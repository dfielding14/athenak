# CGL Landau-Fluid Validation

## Validation Tiers

Routine CPU testing covers a quantitative decay case, analytic uniform
collisional relaxation, oblique/parallel firehose threshold selection, strict
emergency-bound reporting without automatic correction, limiter occupancy,
explicit-versus-STS reference comparisons, and the live CGL FOFC regression.
It also covers the reduced active paper initializer, deterministic turbulence
seed selection, Alfvenic forcing orientation, nonrelativistic forcing-energy
work, once-per-cycle OU forcing advancement, forcing restart continuation, and
passive-Delta flow independence from diagnostic pressure anisotropy.
The broader scientific suite is a manual tier because it generates diagnostic
CSV files and figures and is intended for interpretation, not just pass/fail
gating.

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

For an archived local accuracy campaign with parallel/perpendicular
resolution and timestep sweeps, finite-collision coefficient limits, and
an extreme-cap STS/explicit exercise:

```bash
python3 scripts/cgl_lf_workflow.py accuracy
```

The `accuracy` bundle is operator evidence, not a reduced or
paper-resolution turbulence result. It records pgen reference errors for
`nx1 = 32, 64, 128`, three STS timestep ratios plus explicit references,
both local and background coefficient modes at three collision frequencies,
pre-cap face ratios above `10^4` in the extreme-cap cases, aligned thermal
decay with mean fields directed along `x`, `y`, `z`, and a periodic oblique
direction, and clean transport shutdown at `|B| <= bfloor`.

Finite-collision LF cases use two collision source updates per cycle, each
over the corresponding LF half-sweep duration. The routine
`cgl_lf_collision_relaxation` regression independently checks the analytic
law `Delta p(t) = Delta p(0) exp(-nu_coll t)`.

The explicit mode uses the same anisotropy-to-magnetic-moment split lifecycle
as STS. It is deliberately restricted to standalone reference checks.

For a supported two-dimensional AMR run with periodic boundaries, active LF
transport, strict monitoring, and conserved prolongation:

```bash
python3 scripts/cgl_lf_workflow.py amr
```

For reduced MKS24-oriented turbulence smoke cases, including active Alfvenic,
active random, and passive-Delta Alfvenic forcing:

```bash
python3 scripts/cgl_lf_workflow.py paper-smoke
```

This workflow archives the explicit parallel-firehose policy, forcing mode,
random seed, correlation time, injection parameter, mode bounds, state
tables, forcing tables, and the physical-shell/common-spectrum settings that
select the documented `abs(k) in (2*pi/L_parallel)[1,3]` and `k^-2` forcing.
Paper decks enable restartable cumulative applied source work. The workflow
checks clean LF safety diagnostics, positive applied work, zero parallel
forcing for the Alfvenic case, and nonzero parallel forcing for the random
case. Focused CPU regressions separately check a one-cycle OU transition,
RK2 source-work identity, and cumulative-work restart continuation. The
workflow also archives the explicit passive-Delta choice; a focused
regression checks that its flow fields are independent of diagnostic initial
anisotropy. It does not qualify paper-resolution or long-time active/passive
forcing statistics, or figure-level analysis.

For short reduced nonlinear product-path qualification, including
`8x8x16`, `16x16x32`, and `32x32x64` resolution variants, two heat-flux
wavenumber variants, and an oblique-firehose-policy variant:

```bash
python3 scripts/cgl_lf_workflow.py paper-convergence
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/paper-convergence/bundle
```

This mode overrides the standard input to `tlim = 0.02` and exists to check
nonlinear startup, snapshot output, and analysis/plot plumbing. It is too
short and too coarse for steady-state convergence or direct MKS24 figure
comparison claims.

Paper-standard and limiter-frequency workflow modes are defined for approved
production campaigns, not local validation:

```bash
python3 scripts/cgl_lf_workflow.py paper-standard \
  --authorize-paper-execution --output-dir /path/to/archive/run
python3 scripts/cgl_lf_workflow.py paper-nulim \
  --authorize-paper-execution --output-dir /path/to/archive/run
```

Their inputs use `192x192x384`, `tlim = 10`, full-field binary snapshots,
checkpoints, and the `t = [8,10]` analysis window. Without
`--authorize-paper-execution` the workflow rejects these expensive modes
before creating a run bundle. Defining these inputs does not constitute a
paper-standard execution or reproduction result.
In CGL primitive snapshots, legacy `eint` is `p_parallel` and the dedicated
`p_perp` output supplies the perpendicular pressure required for paper
anisotropy analysis.

Analyze a retained paper bundle without rerunning simulations:

```bash
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/existing/bundle
```

This action writes `analysis/diagnostics.json`, summarizes reduced histories,
and analyzes retained CGL binary snapshots when present. For standard inputs,
it reads each case's archived `analysis_t_start`/`analysis_t_end` values and
forms common-bin steady-window averages. Current products include
pressure/density PDFs, perpendicular shell spectra, field-projected
`Delta p` and velocity-gradient spectra, the `b b : grad u` strain PDF,
pressure-stress transfer with a real-space closure check, and selected-scale
alignment PDFs. When the run manifest contains the complete LF closure
choices, it also reconstructs a cell-centered heat-flux smoothing proxy,
`integral[-q_parallel b.grad(T_parallel) - q_perp b.grad(T_perp)] dV`,
including regularized/unlimited and cap-active contributions. This product
uses snapshot central differences and is explicitly not the applied
finite-volume face flux or a closed energy budget.

The analyzer runs synthetic binning/gradient/transfer/alignment and
heat-flux-sign checks and renders available history, PDF, spectrum, transfer,
alignment, and heat-flux-proxy figures under `figures/paper/`. Panel-specific
reference-curve comparison and a full anisotropic-pressure/local-work budget
decomposition remain to be completed. For fixed-level runs, the analyzer
also reports the RKL2-applied capped-face heat-flux contractions retained in
`lf_qprwrk` and `lf_qpewrk`; this is not an AMR-corrected or total-energy
budget diagnostic, and the signed value need not equal the positive
cell-centered snapshot proxy. Operator-face heat-flux-cap fractions are separately
summarized from LF history counters over the same selected interval. When the
default staged MKS24 manifest is available, `paper-analyze` also writes
`analysis/reference_provenance.json` with its archive and source-TeX
checksums; provide a different pinned copy with
`--reference-manifest <path>`.
For a standalone LF history file, the same interval calculation is available
through `scripts/analyze_cgl_lf_paper.py --lf-history <path>`.

Regenerate only the readable summary and diagnostic checks for an existing
paper bundle without plotting:

```bash
python3 scripts/cgl_lf_workflow.py paper-summary \
  --output-dir /path/to/existing/bundle
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

`manifest.json` records the git revision, whether the worktree was dirty at
execution, its short status entries, the executable, source inputs, runtime
overrides, exact commands, result products, and measured checks. A dirty
bundle is suitable for implementation testing but is not evidence from an
immutable production source revision. `summary.md` is the readable pass/fail
digest. The `amr` summary reports whether refinement actually occurred,
normalized divB, invalid-state count, and normalized energy residual. The
`compare` summary reports maximum differences between STS and explicit final
states. The `paper-smoke` summary reports forcing orientation, mean beta,
`C_B2`, threshold-volume fractions, and total-energy change while retaining
both state and forcing tables.

The workflow also supports regenerating figures from a retained `full` bundle
with `plot`, and rebuilding `summary.md` and diagnostic values with
`summarize`, both using `--output-dir` to select that existing bundle.
Paper bundles use `paper-analyze` for JSON plus diagnostic plots and
`paper-summary` for summary-only regeneration.
Generated bundles are run products and are not committed by default.

## Paper Reference Staging

Paper-comparison analysis must use a pinned MKS24 reference copy rather than
an unversioned download. Stage the official arXiv `2405.02418v2` source
bundle into an ignored results area:

```bash
python3 scripts/stage_cgl_lf_mks24_reference.py
```

The utility downloads the source archive, safely extracts it, requires
`MKS24.tex`, and writes `manifest.json` containing SHA-256 checksums for the
archive and every extracted source or figure file. For an archive already
transferred to a compute system, use:

```bash
python3 scripts/stage_cgl_lf_mks24_reference.py \
  --archive /path/to/arXiv-2405.02418v2-source.tar.gz \
  --expected-sha256 <checksum-from-an-approved-manifest> \
  --output-dir /path/to/immutable/reference/arXiv-2405.02418v2
```

Do not commit redistributed paper source or figures as code artifacts unless
their redistribution permission is separately established. Archive the
staging manifest or immutable reference path with every quantitative paper
comparison bundle. The staged official `2405.02418v2` source inventory
contains TeX and rendered figure PDFs, not machine-readable curve tables.
Quantitative panel comparison therefore still requires a recorded
digitization procedure or a separately provenance-tracked numerical
reference source.

## Frontier Debug Qualification Preparation

Frontier qualification uses tracked operational tooling but remains a
separate, manually authorized execution tier. Validate the preparation and
accounting rules offline before copying the scripts to the project run root:

```bash
python3 scripts/frontier/cgl_lf_frontier.py self-test
```

For an actual Frontier debug-only qualification campaign, first build an
immutable HIP/MPI executable from a clean committed checkout using:

```bash
scripts/frontier/build_cgl_lf_frontier.sh
```

Then initialize the prescribed project area and prepare a reduced validation
job. The utility writes a retained manifest and concrete `debug`-QOS batch
script; it does not call `sbatch`.

```bash
python3 scripts/frontier/cgl_lf_frontier.py init
python3 scripts/frontier/cgl_lf_frontier.py prepare \
  --campaign qualification --run-name g002_paper_smoke_active \
  --purpose "strict reduced CGL-LF GPU smoke" \
  --acceptance-criterion "clean LF safety diagnostics and finite output" \
  --executable /lustre/orion/ast207/proj-shared/dfielding/CGL/build/<build>/src/athena \
  --input-file /lustre/orion/ast207/proj-shared/dfielding/CGL/repo/athenak-DF/inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput \
  --nodes 1 --walltime 00:20:00 --athena-walltime 00:15:00
python3 scripts/frontier/cgl_lf_frontier.py check-submit \
  --manifest /lustre/orion/ast207/proj-shared/dfielding/CGL/runs/qualification/g002_paper_smoke_active/manifest/prepared_run.json
```

Only after inspecting the printed `sbatch` command and confirming the
preflight should the job be submitted manually. Record the returned job ID
and completed accounting record sequentially:

```bash
python3 scripts/frontier/cgl_lf_frontier.py mark-submitted \
  --manifest <prepared_run.json> --job-id <jobid>
python3 scripts/frontier/cgl_lf_frontier.py record \
  --manifest <prepared_run.json> --job-id <jobid> \
  --result passed --notes "inspected reduced diagnostics"
```

For real runs the utility requires all source, executable, input, and output
locations beneath `/lustre/orion/ast207/proj-shared/dfielding/CGL`, limits
debug walltime to two hours, reserves against the 1000 node-hour testing
budget, checks that no other debug job is queued, and rejects
`paper-standard` and `paper-nulim` inputs. Paper-production simulations must
not be run through this debug-only workflow.

## Diagnostics

With LF active, normal MHD history output appends cumulative counters:
`lf_nstage`, `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
`lf_mirror`, `lf_firehs`, `lf_hardbd`, `lf_qface`, `lf_qprcap`,
`lf_qpr10`, `lf_qpecap`, `lf_qpe10`, `lf_qprwrk`, and `lf_qpewrk`. Strict
validation decks require zero floors, zero nonfinite/nonpositive states, and
zero hard-bound violations; limiter and cap counts may be nonzero when
intentionally exercised.
`lf_mirror` and `lf_firehs` are physical threshold-occupancy counters
evaluated with the selected policy. `lf_hardbd` is an emergency safety
counter and is evaluated whether or not corrective `backup_limiters` are
enabled. The `lf_q*` columns record evaluated LF operator faces and
parallel/perpendicular pre-cap flux ratios above `q_max` and `10*q_max`;
their interval fractions use `lf_qface` as denominator. For fixed-level
meshes, `lf_qprwrk` and `lf_qpewrk` retain restartable cumulative
RKL2-applied owned-face contractions of capped parallel and perpendicular
heat flux with discrete temperature jumps. Refined meshes intentionally
leave these work columns zero because this diagnostic does not include AMR
flux correction. These are signed operator contractions, not a positivity or
total-energy closure condition.

The firehose policy is explicit for threshold-sensitive work:
`cgl_firehose_threshold = parallel` selects the MKS24 production convention
`beta Delta <= -2`, while `oblique` selects `beta Delta <= -1.4` for
comparison studies and backward-compatible operator validation.

The pre-existing `aam-D` label is retained for compatibility. In ordinary
output and restart state it denotes conserved CGL pressure anisotropy, not
the temporary magnetic-moment representation used internally during LF split
sweeps.

The `cgl_lf_paper` user history stores volume integrals named `mass`,
`kinetic`, `magnetic`, `therm_cgl`, `b2`, `b4`, `delta_p`, `abs_dp`,
`beta`, `mirror_vol`, `fire_vol`, `hard_vol`, `nu_eff`, `force_pwr`, and
`force_work`. `force_pwr` is an instantaneous proxy; `force_work` is the
cumulative net energy source actually applied by stage-weighted forcing,
including its zero-net-momentum projection, when the paper input enables
`record_injected_work`. Form means or fractions
using `volume`; for example, `C_B2 = b4*volume/b2^2 - 1`.
Threshold-volume columns use the selected firehose policy, whereas
`hard_vol` is safety-only. For active-Delta bundles, `paper-analyze` compares
the `force_work` interval difference to the conserved MHD `tot-E` interval
difference and reports the global residual for qualification. Passive-Delta bundles
retain applied work but are not labeled as an active-CGL energy budget.

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
