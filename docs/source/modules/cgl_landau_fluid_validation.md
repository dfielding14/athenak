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

## Quantitative 1D Example

Build the default binary and run a strict parallel-temperature damping check:

```bash
cmake -S . -B build-cgl-implementation -DCMAKE_BUILD_TYPE=Release
cmake --build build-cgl-implementation -j
build-cgl-implementation/src/athena \
  -i inputs/unit_tests/cgl_lf_quant_parallel.athinput
```

To compare against the verification integrator, use the routine decay input:

```bash
build-cgl-implementation/src/athena \
  -i inputs/tests/cgl_lf_decay.athinput \
  job/basename=cgl_lf_explicit_check \
  mhd/cgl_heat_flux_integrator=explicit \
  time/sts_integrator=none time/nlim=-1
```

The explicit mode uses the same anisotropy-to-magnetic-moment split lifecycle
as STS. It is deliberately restricted to standalone reference checks.

## Extended Workflow

Run the reproducible validation matrix from the repository root:

```bash
scripts/run_cgl_lf_validation.sh
```

By default this writes CSV logs and plots below
`build-cgl-implementation/cgl_lf_validation/`, not into the source tree.
Use `OUTPUT_DIR=/path/to/results` to retain a named validation bundle.

Regenerate exact oblique-background eigenmode input decks after intentionally
changing the linear reference convention:

```bash
python3 scripts/generate_cgl_lf_eigenmode_inputs.py
scripts/run_cgl_lf_validation.sh
```

The plotting script can also be run directly against an archived data
directory:

```bash
python3 scripts/plot_cgl_lf_validation.py \
  --data-dir /path/to/results/data \
  --figure-dir /path/to/results/figures
```

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

## Interpretation Limits

The oblique and eigenmode decks compare this ion CGL-LF closure against
linearized reference systems constructed for the same initial-value problem.
They are useful numerical checks, but do not by themselves implement electron
pressure physics or reproduce every historical comparison model.
