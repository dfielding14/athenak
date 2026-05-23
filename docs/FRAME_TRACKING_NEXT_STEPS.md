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
| Correctness | Continuous and restart-split runs agree in controller state and final fluid state within precision-aware tolerances. | One- and three-axis restart equivalence tests pass locally at `100 * machine_epsilon`; GitHub Actions confirmation remains pending. |
| Robustness | CPU, MPI/AMR, MHD, invalid configuration, legacy-state, cloud, and TRML tests run through reproducible commands. | Full local CPU and MPI CPU discovery suites pass, including MHD, three-axis, guards, legacy, cloud, and TRML coverage; medium campaigns remain pending. |
| Observability | Controller state is emitted in machine-readable history data. | Implemented as `<basename>.frame_tracker.hst`; stdout is compatibility-only for archived runs. |
| Efficiency | All active axes share one mesh sampling pass and grouped reductions per controller update. | Implemented; benchmark evidence remains to be collected. |
| Usability | Canonical parameters, conservative examples, removed-alias migration errors, and supported-physics limits are documented. | Implemented in the module reference and recipe/migration page; legacy restart-state reading remains supported. |
| Publication | GitHub Pages exposes reference pages, slices, CSV summaries, validation limits, and this roadmap. | Existing wiring-validation pages are published; recipe/migration updates and new measured results remain to be synced after validation. |

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
| 5 | Fuse multi-axis sampling and group global reductions per update. | Benchmark table records update costs for one and three axes. |
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
| Cloud/TRML integrations | Short custom-problem runs compile, execute, and generate finite structured histories; quantitative medium-resolution CSV summaries remain pending. |

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
mask style findings. CI confirmation for the exact pushed commit, benchmark
evidence, and medium-resolution scientific validation are still required.

## Benchmark Campaign

Run benchmarks only after correctness tests pass. Record wall-clock tracker
update cost, hardware, compiler configuration, mesh/block layout, MPI ranks,
active axes, and target selection.

| Benchmark | Acceptance criterion |
| --- | --- |
| Single-axis CPU | Median tracker-update time does not regress by more than 5 percent relative to the pre-fusion baseline. |
| Three-axis CPU | Median tracker-update time improves by at least 25 percent on a representative multi-block problem. |
| MPI CPU | Grouped reductions preserve the CPU result within reduction-order tolerance and reduce per-update collective overhead. |
| GPU, when available | Record results and kernel timing; do not infer GPU performance from CPU data. |

## Scientific Validation Campaign

1. Run cloud and TRML examples with tracking disabled and enabled, comparing
   fields after transforming moving-frame results back to lab coordinates.
2. Repeat at multiple spatial resolutions and publish convergence tables.
3. Split runs across restarts at more than one controller phase and compare
   final fields and histories.
4. Repeat representative calculations across MPI decompositions and AMR
   layouts.
5. Publish plots, raw CSV data, commands, software revision, tolerances,
   controller limit/miss events, and all observed failures.

Until that campaign is complete, recommend the feature for controlled example
and method-development work, not irreversible production science runs.
