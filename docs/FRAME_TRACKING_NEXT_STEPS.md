# Frame Tracker Production Hardening And Validation Plan

## Purpose

`FrameTracker` provides a non-relativistic moving-frame controller for AthenaK
hydro and MHD calculations. The low-resolution cloud and TRML examples are
wiring validations: they establish that frame-aware boundary conditions,
controller output, restart state, and plotting machinery are connected
correctly. They are not convergence studies or production-science
certifications.

This document is the maintained implementation and validation plan. Results
must be reported as measured quantities with stated tolerances; do not describe
the feature as working perfectly.

## Production Gates

| Gate | Required outcome | Current branch status |
| --- | --- | --- |
| Correctness | Continuous and restart-split runs agree in controller state and final fluid state within precision-aware tolerances. | `data_precision=real` and a strict material TRML restart regression are implemented and pass at test scale. Medium production restart certification is deferred because the serial cloud physical gate fails first. |
| Robustness | CPU, MPI/AMR, MHD, invalid configuration, legacy-state, cloud, and TRML tests run through reproducible commands. | Full local CPU and MPI CPU discovery suites pass, including tracer weighting, scalar AMR parsing, native-precision output, material cloud/TRML smoke coverage, and strict short restart comparison. |
| Observability | Controller state is emitted in machine-readable history data. | Implemented as `<basename>.frame_tracker.hst`; stdout is compatibility-only for archived runs. |
| Efficiency | All active axes share one mesh sampling pass and grouped reductions per controller update. | Measured on clean candidate `78f8e4cd`: serial `x1` overhead is 77.0% lower and serial `all` overhead is 81.4% lower than `cfbcdde7`; see the performance-evidence page and raw CSV. |
| Usability | Canonical parameters, conservative examples, removed-alias migration errors, and supported-physics limits are documented. | Material recipes now use `target=scalar0` and `weight=tracer_mass`; density/temperature recipes are retained as phase diagnostics; legacy restart-state reading remains supported. |
| Publication | GitHub Pages exposes reference pages, slices, CSV summaries, validation limits, and this roadmap. | Material-validation CSV and medium tracer-slice assets are prepared for publication without upgrading the scientific claim; Pages sync remains pending. |

## Interface Contract

### Restart State

New restart output records `state_version=1`. It preserves global controller
timing and reacquisition counters, frame velocity and displacement, and every
mutable per-axis control value required to continue without re-priming filters
or resetting rate/integral state.

Restart files containing only `frame_velocity_x*` and
`frame_displacement_x*` remain readable. This compatibility path initializes
unavailable history conservatively, emits one rank-0 warning, and is explicitly
not exact continuation.

Configuration aliases are no longer accepted. Former spellings fail at startup
with their canonical replacement; this removes hidden configuration paths while
leaving the legacy restart-state reader intact.

### Structured Diagnostics

When a history output is active, frame tracking writes a dedicated
`<basename>.frame_tracker.hst` stream. The maximum three-axis record has
exactly 20 controller columns:

| Scope | Columns |
| --- | --- |
| Each active axis `xN` | `ft_vf_xN`, `ft_dx_xN`, `ft_pos_xN`, `ft_err_xN`, `ft_dv_xN` |
| Shared | `ft_weight`, `ft_misses`, `ft_recov`, `ft_limit`, `ft_skip` |

Values represent global controller state. Rank 0 supplies those values before
the normal history reduction; other ranks contribute zeros.

### Supported Physics

Frame tracking supports non-relativistic hydro or MHD. Set
`tracked_fluid=hydro` or `tracked_fluid=mhd` whenever both eligible fluids
exist. SR, GR, and dynamical-GR coordinates are rejected. General
problem-specific boundaries and sources cannot be transformed automatically;
only shipped frame-aware examples are validated.

## Implementation Record

| Phase | Work item | Release criterion |
| --- | --- | --- |
| 1 | Publish reference docs, example plots, CSVs, and this roadmap through the GitHub Pages docs tree. | Sphinx/MyST build passes and live links render assets. |
| 2 | Version and serialize complete controller state, retain the legacy reader, and preserve full `Real` precision only for stored controller state. | Restart-split test passes at `100 * machine_epsilon` scaled tolerance without changing unrelated runtime parameter formatting. |
| 3 | Emit dedicated history output and migrate plotting to consume it before stdout fallback. | Scripts reproduce tables from `.frame_tracker.hst`. |
| 4 | Add CPU smoke, restart, MHD, MPI/AMR, legacy compatibility, and invalid-input regressions. | Tests run through `tst/run_test_suite.py`. |
| 5 | Fuse multi-axis sampling, group global reductions, and skip no-op boost/timestep kernels. | Clean-revision benchmark passes the declared serial one-axis and three-axis criteria; raw CPU/MPI results are maintained for publication. |
| 6 | Use time-based slew limiting in shipped examples, reject configuration aliases actionably, and publish copy-ready recipes. | New docs and inputs contain only canonical parameter names; legacy restart state is still readable. |
| 7 | Execute scientific validation at useful resolution and decomposition. | Published CSVs contain measured errors, convergence trends, and limits. |

## Automated Test Matrix

| Scenario | Assertion |
| --- | --- |
| Serial controller smoke | Structured history is finite; the controller primes then actuates; no unexpected misses occur. |
| Restart continuity | Uninterrupted and restart-split material runs with `data_precision=real` match final frame state and conservative fluid output to `100 * machine_epsilon` scaled tolerance. |
| Legacy restart keys | Legacy state loads, warns once, and produces finite diagnostics; exact continuation is not claimed. |
| MHD Galilean invariant | Density and magnetic fields are unchanged at the first non-zero boost; momentum and ideal-gas energy match analytic boost updates to `100 * machine_epsilon`. |
| MPI/AMR | One-rank and multi-rank AMR regression results agree within reduction-order tolerance, including an AMR criterion selecting `hydro_w_s00`. |
| Guards | SR, GR, dynamical-GR, missing selected fluids, ambiguous eligible fluids, and removed aliases fail before evolution with actionable messages. |
| Cloud/TRML integrations | Short custom-problem material runs compile, emit native-precision tracer snapshots, reject invalid cloud tracer indices, and execute a strict short TRML material restart comparison. |

## Local Verification Record

The material-tracer, precision-output, scalar-AMR, regression-test, and
follow-up material-recipe/boundary changes were validated on May 24, 2026 in
the feature worktree beginning with feature commit `a893c52f`:

| Command | Result |
| --- | --- |
| `python run_test_suite.py --style` | `2 passed` |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_lwave1d_cpu.py` | `16 passed`; confirms precise restart storage does not perturb existing linear-wave runtime input updates. |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_frame_tracking_cpu.py` | `20 passed`; includes MHD invariant, canonical/removed-input behavior, tracer-mass centroid, invalid scalar AMR index, initialization summary, and three-axis schema. |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_frame_tracking_restart_cpu.py` | `2 passed`; one- and three-axis restart continuation. |
| `python -m pytest test_suite/nr/test_nr_frame_tracking_examples_cpu.py -q` | `2 passed`; centroid material inputs, moving-frame constant cloud inflow, output precision, cloud tracer-index guard, and strict short TRML material restart comparison. |
| `python run_test_suite.py --mpicpu --test test_suite/nr/test_nr_frame_tracking_amr_mpicpu.py` | `3 passed`; serial/MPI AMR comparison including scalar-targeted refinement. |
| `python run_test_suite.py --cpu` | `212 passed, 15 skipped, 61 deselected`; includes cloud and TRML custom-problem integrations. |
| `python run_test_suite.py --mpicpu` | `32 passed, 256 deselected`. |
| `python scripts/diagnose_cloud_frame_discrepancy.py ...` | `34` diagnostic rows plus a localization figure; corrected Cloud discrepancy is localized to the shocked-cloud interaction region. |

No GitHub Actions runs were visible for the exact pushed feature commit when
queried locally; the table above is local evidence. The medium-resolution
serial diagnostics below retain a failed production gate.

## Benchmark Campaign

The benchmark was run after the correctness tests passed. It uses a uniform
zero-velocity material selection in `mode=velocity`, so the enabled and
disabled fluid evolutions are identical and the measured increment isolates
steady sampling/controller overhead. It records Release-build wall-clock
costs for `64 x 64 x 64` cells, `16 x 16 x 16` blocks, 200 cycles, one
discarded warm-up, and seven measured repeats per case.

| Benchmark | `cfbcdde7` | `78f8e4cd` | Acceptance criterion | Result |
| --- | ---: | ---: | --- | --- |
| Single-axis CPU | `3.594 ms/update` | `0.828 ms/update` | No more than 5 percent regression. | Pass; 77.0 percent lower. |
| Three-axis CPU | `5.588 ms/update` | `1.041 ms/update` | At least 25 percent lower. | Pass; 81.4 percent lower. |
| Single-axis MPI CPU, 4 ranks | `0.887 ms/update` | `0.233 ms/update` | Recorded comparison. | Measured; 73.7 percent lower. |
| Three-axis MPI CPU, 4 ranks | `1.443 ms/update` | `0.441 ms/update` | Recorded comparison. | Measured; 69.5 percent lower. |
| GPU | Not run | Not run | Do not infer GPU performance from CPU data. | Pending follow-up evidence. |

The raw data and compact plot are maintained as
`docs/source/_static/frame_tracking_benchmark.csv` and
`docs/source/_static/frame_tracking_benchmark.png`. These measurements do not
replace cloud/TRML transformed-frame physical validation. The benchmark uses
strictly periodic boundaries and therefore does not measure the required
non-periodic ghost-boundary refresh added for restart-correct frame-aware
examples; that path converts ghost slabs only, but still needs a dedicated
boundary-aware timing result before production performance claims.

## Scientific Validation Campaign

The May 24, 2026 campaign replaces threshold-selected material acceptance with
advected tracer mass and tracer-mass centroid. Follow-up diagnosis found that
the initial material recipe's default blended position signal allowed
numerically small positive tracer tails to influence the controller through
its selected-band midpoint. The maintained material recipes now set
`position_signal=centroid`. Density-selected cloud mass and
temperature-selected TRML mass remain mandatory phase-structure diagnostics.
The maintained serial table is
`docs/source/_static/frame_tracking_material_validation_summary.csv`, with
final medium slices in
`docs/source/_static/frame_tracking_material_medium_slices.png`.

| Problem/comparison | Resolution and duration | Measured result | Acceptance result |
| --- | --- | --- | --- |
| Cloud Sedov centroid material tracking on/off | `48 x 16 x 16` and `96 x 32 x 32`, `tlim=0.04` | Tracer-mass relative differences are `7.56e-16` then `0.00`; centroid error decreases to `3.73e-7`; tracer-density L2 decreases from `2.9779e-5` to `2.4588e-5`; no health events. | **Fail:** aggregate lab-frame conserved-state L2 is `9.0074e-5` then `3.1725e-4`, not decreasing. |
| Cloud Sedov localized diagnosis | Same corrected runs | At medium resolution `99.99999%` of aggregate squared error lies at `1 <= x1 < 4`; the `x1 < 0` inflow fraction is `1.22e-7`; integrated energy relative difference improves from `4.51e-7` to `1.05e-9`. | The remaining failure is localized to shocked-cloud field structure, not the inflow boundary. |
| TRML centroid material tracking on/off | `16 x 16 x 32` and `32 x 32 x 64`, `tlim=0.25` | With `max_boost_change_rate=1.00`, tracer mass differences are `6.59e-3` and `7.60e-3`; tracer-density L2 decreases from `6.0174e-3` to `3.1076e-3`; conserved-state L2 decreases from `4.9300e-3` to `2.6596e-3`; no health events. | Pass serial material gate. |
| TRML temperature-window phase diagnostic | Same resolutions and duration | Selected-mass difference is `4.0290e-2` then `6.0872e-2`, with limit events. | Retained phase diagnostic; not a material-retention gate. |

The corrected TRML centroid recipe requires a new low-resolution slew-rate
sweep. Rates `0.02`, `0.05`, `0.10`, and `0.20` incur `9`, `4`, `2`, and `1`
limit events respectively; `1.00` is the first tested no-limit rate and is
confirmed at medium resolution.

The cloud centroid correction reduces the previous low-resolution conserved
discrepancy by nearly two orders of magnitude and resolves the former
density-threshold and band-midpoint ambiguities. It also leaves a precise
blocker: the transformed conserved-state discrepancy in the time-dependent,
shocked Sedov problem does not show the required decreasing resolution trend.
A boundary audit corrected constant lab-inflow velocity transformation; the
Sedov inflow was already transformed correctly. The localized output in
`docs/source/_static/frame_tracking_sedov_cloud_localization.csv` demonstrates
that the remaining discrepancy is in the shocked cloud, not at the boundary.

The restart mismatch diagnosed in the prior campaign was caused by stale
frame-dependent ghost states after post-timestep frame updates. That defect is
fixed. The new `data_precision=real` output path and strict short material
restart regression are implemented, but no medium production restart campaign
is run while the serial cloud physical gate fails.

Required next work is:

1. Determine why small centroid-controlled frame motion produces increasing
   pointwise shocked-cloud field differences even while material and
   integrated energy observables improve.
2. Either make a validated numerical/controller correction or justify a
   discontinuous-flow acceptance metric with additional published evidence.
3. Only after the cloud serial trend passes, execute strict medium restart,
   four-rank MPI, and tracer-refined AMR comparisons for both examples.

Until those results pass, recommend the feature for controlled wiring tests
and method-development work only, not irreversible production science runs.
