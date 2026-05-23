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
| Correctness | Continuous and restart-split runs agree in controller state and final fluid state within precision-aware tolerances. | Small regression fixtures pass locally, but medium TRML split runs exceed the strict continuation tolerance; production gate fails pending investigation. |
| Robustness | CPU, MPI/AMR, MHD, invalid configuration, legacy-state, cloud, and TRML tests run through reproducible commands. | Full local CPU and MPI CPU discovery suites pass, including MHD, three-axis, guards, legacy, cloud, and TRML coverage; medium serial scientific comparisons expose blocking failures. |
| Observability | Controller state is emitted in machine-readable history data. | Implemented as `<basename>.frame_tracker.hst`; stdout is compatibility-only for archived runs. |
| Efficiency | All active axes share one mesh sampling pass and grouped reductions per controller update. | Measured on clean candidate `78f8e4cd`: serial `x1` overhead is 77.0% lower and serial `all` overhead is 81.4% lower than `cfbcdde7`; see the performance-evidence page and raw CSV. |
| Usability | Canonical parameters, conservative examples, removed-alias migration errors, and supported-physics limits are documented. | Implemented in the module reference and recipe/migration page; legacy restart-state reading remains supported. |
| Publication | GitHub Pages exposes reference pages, slices, CSV summaries, validation limits, and this roadmap. | Existing wiring-validation pages are published; recipe, benchmark, and failed-medium-validation evidence are prepared for sync without upgrading the scientific claim. |

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
| Restart continuity | Uninterrupted and restart-split runs match final frame state and conservative fluid output to `100 * machine_epsilon` scaled tolerance. |
| Legacy restart keys | Legacy state loads, warns once, and produces finite diagnostics; exact continuation is not claimed. |
| MHD Galilean invariant | Density and magnetic fields are unchanged at the first non-zero boost; momentum and ideal-gas energy match analytic boost updates to `100 * machine_epsilon`. |
| MPI/AMR | One-rank and multi-rank AMR results agree within reduction-order tolerance and retain valid controller state through refinement. |
| Guards | SR, GR, dynamical-GR, missing selected fluids, ambiguous eligible fluids, and removed aliases fail before evolution with actionable messages. |
| Cloud/TRML integrations | Short custom-problem runs compile and execute with finite structured histories; medium serial comparisons are recorded and fail the declared physical acceptance criteria. |

## Local Verification Record

The interface and regression-test changes were validated on May 23, 2026
before their feature-branch commit:

| Command | Result |
| --- | --- |
| `python run_test_suite.py --style` | `2 passed` |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_lwave1d_cpu.py` | `16 passed`; confirms precise restart storage does not perturb existing linear-wave runtime input updates. |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_frame_tracking_cpu.py` | `17 passed`; includes MHD invariant, canonical/removed-input behavior, initialization summary, and three-axis schema. |
| `python run_test_suite.py --cpu --test test_suite/nr/test_nr_frame_tracking_restart_cpu.py` | `2 passed`; one- and three-axis restart continuation. |
| `python run_test_suite.py --mpicpu --test test_suite/nr/test_nr_frame_tracking_amr_mpicpu.py` | `2 passed`; serial/MPI AMR comparison. |
| `python run_test_suite.py --cpu` | `209 passed, 15 skipped, 60 deselected`; includes cloud and TRML custom-problem integrations. |
| `python run_test_suite.py --mpicpu` | `31 passed, 253 deselected`. |

The style command uses a local warning suppression for a deprecation warning
in the test harness's downloaded `cpplint.py` under Python 3.14; it does not
mask style findings. CI confirmation for the exact pushed commit and
medium-resolution scientific validation are still required.

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
replace cloud/TRML transformed-frame physical validation.

## Scientific Validation Campaign

Medium-resolution serial validation was run on May 23, 2026 with the Release
CPU executable at `505f2df4`; the inputs now expose `enabled = true`
explicitly so the same documented configuration can be switched off at run
time. The common comparison tool writes
`docs/source/_static/frame_tracking_validation_summary.csv`.

| Problem/comparison | Resolution and duration | Measured result | Acceptance result |
| --- | --- | --- | --- |
| Cloud tracking on/off, lab-frame comparison | `96 x 32 x 32`, `tlim=0.04` | Selected-mass relative difference `1.4178e-2`; zero misses and limits. | Fail: mass difference exceeds `1.0e-2`. |
| TRML tracking on/off, lab-frame comparison | `32 x 32 x 64`, `tlim=0.25` | Selected-mass relative difference `6.0872e-2`; three limit events. | Fail: mass difference and no-limit requirements are not met. |
| TRML restart splits at approximately 25, 50, and 75 percent | `32 x 32 x 64`, `tlim=0.25` | Finite output and small controller differences, but strict field/controller rows fail. | Fail: `100 * machine_epsilon` continuation criterion is not met. |
| TRML serial versus 4-rank MPI, uniform grid | `32 x 32 x 64`, `tlim=0.25` | Field/controller comparisons pass `1.0e-10`; three limit events remain. | Numerical comparison passes; health criterion fails. |

The failed serial physical comparisons are already sufficient to block
production-candidate guidance. Cloud restart/MPI/AMR and full parallel
scientific comparisons are deferred until the selected-mass discrepancy and
TRML slew limiting are investigated; running them now could add diagnostics
but cannot satisfy this release gate.

Required next work is:

1. Diagnose transformed tracking-on/off mass divergence and validate the
   selection and boundary/source transformations.
2. Tune or justify the conservative TRML controller limits and rerun the
   medium serial comparison without unexplained limit events.
3. Resolve the strict medium restart-continuity mismatch.
4. Rerun serial, restart, MPI, and AMR comparisons and publish the updated CSV
   and field-norm tables.

Until those results pass, recommend the feature for controlled wiring tests
and method-development work only, not irreversible production science runs.
