# CGL Landau-Fluid Validation

## Validation Tiers

Routine CPU testing covers a quantitative decay case, limiter occupancy with
strict admissibility, explicit-versus-STS reference comparisons, and the live
CGL FOFC regression. The broader scientific suite is a manual tier because it
generates diagnostic CSV files and figures and is intended for interpretation,
not just pass/fail gating.

The manual suite covers:

- collisionless and finite-collision parallel and perpendicular damping;
- isolated field-strength-gradient transport and heat-flux limiting;
- mirror and firehose limiter activity;
- field-aligned and oblique linear waves;
- exact pure-CGL versus CGL-LF eigenmode comparisons.

## Canonical Workflows

Use the workflow driver from the repository root. It builds
`build-cgl-implementation/src/athena` in Release mode when that executable is
not already available.

For a first validated 1D STS damping run with strict diagnostics:

```bash
python3 scripts/cgl_lf_workflow.py quick
```

For the explicit-reference check against STS, including the finite-collision
split case:

```bash
python3 scripts/cgl_lf_workflow.py compare
```

The explicit mode uses the same anisotropy-to-magnetic-moment split lifecycle
as STS. It is deliberately restricted to standalone reference checks.

For a supported two-dimensional AMR run with periodic boundaries, active LF
transport, strict monitoring, and conserved prolongation:

```bash
python3 scripts/cgl_lf_workflow.py amr
```

For the full scientific validation matrix and plots:

```bash
python3 scripts/cgl_lf_workflow.py full
```

The compatibility command `scripts/run_cgl_lf_validation.sh` invokes this
same `full` workflow. Use `--output-dir /path/to/results` with the Python
driver, or `OUTPUT_DIR=/path/to/results` with the compatibility command, to
give a bundle a persistent name.

## Result Bundles

Each executable workflow writes an ignored result bundle below
`build-cgl-implementation/cgl_lf_runs/<timestamp>-<workflow>/` by default:

```text
manifest.json
summary.md
logs/
history/
data/
figures/
```

`manifest.json` records the git revision, executable, source inputs, runtime
overrides, exact commands, result products, and measured checks. `summary.md`
is the readable pass/fail digest. The `amr` summary reports whether refinement
actually occurred, normalized divB, invalid-state count, and normalized
energy residual. The `compare` summary reports maximum differences between
STS and explicit final states.

The workflow also supports regenerating figures from a retained `full` bundle
with `plot`, and rebuilding `summary.md` and diagnostic values with
`summarize`, both using `--output-dir` to select that existing bundle.
Generated bundles are run products and are not committed by default.

## Diagnostics

With LF active, normal MHD history output appends cumulative counters:
`lf_nstage`, `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
`lf_mirror`, `lf_firehs`, and `lf_hardbd`. Strict validation decks require
zero floors, zero nonfinite/nonpositive states, and zero hard-bound
violations; limiter counts may be nonzero when a limiter is intentionally
exercised.

The pre-existing `aam-D` label is retained for compatibility. In ordinary
output and restart state it denotes conserved CGL pressure anisotropy, not
the temporary magnetic-moment representation used internally during LF split
sweeps.

The AMR workflow retains both `.user.hst` and `.mhd.hst` products: user
history provides normalized divB, invalid-state, and anisotropy measures,
while MHD history provides total energy and LF counters. CGL LF rejects
`mesh_refinement/prolong_primitives=true`; conserved prolongation is the
supported AMR path.

## Developer Maintenance

Regenerate exact oblique-background eigenmode input decks only after
intentionally changing the linear reference convention with
`scripts/generate_cgl_lf_eigenmode_inputs.py`, then rerun the `full`
workflow. Routine pass/fail regressions remain under `tst/test_suite/cgl/`.

## Interpretation Limits

The oblique and eigenmode decks compare this ion CGL-LF closure against
linearized reference systems constructed for the same initial-value problem.
They are useful numerical checks, but do not by themselves implement electron
pressure physics or reproduce every historical comparison model.
