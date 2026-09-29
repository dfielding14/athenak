# WO1 CGL-LF changes and validation

Status: in progress. Base: `8222de3aa`; branch: `c/cgl-lf-wo1`.
The work order controls the design; `CGL_LF_STS_review.md` describes pre-merge
`2aee609` and is used for derivations rather than current bug status.

## Baseline

- Kokkos initialized from the local populated checkout at unchanged `08ceff92`.
- Release CPU Serial build: `cmake -S . -B tst/build -DCMAKE_BUILD_TYPE=Release`,
  then `cmake --build tst/build -j 8` (passed).
- Current built-in CGL CPU suite: 78 passed, 5 pre-existing failures. The failures
  are `test_cgl_lf_paper_production_inputs_explicitly_use_rank_local_io`,
  `test_cgl_lf_stage_i_requires_retained_source_bundle_provenance`,
  `test_cgl_lf_stage_i_authenticates_historical_production_utility`,
  `test_cgl_lf_stage_i_isolates_epoch_and_checks_all_shared_root_jobs`, and
  `test_cgl_lf_stage_i_recovers_ambiguous_atomic_submit`. These concern archived
  campaign inputs/provenance and Linux-specific executable paths, including
  `/usr/bin/python3.11`, unavailable on this macOS host. They are outside WO1.
- Python tests use `/Users/dbf75/.uv/envs/interactive/.venv/bin/python3`.
  Import `test_suite.testutils` from `tst`, then run pytest from `tst/build/src`,
  matching the repository runner's working directory without rebuilding each time.

## Completed tasks

### T-A1: fast magnetosonic speed

Triage: production formula already fixed. `CheckDirectionalLimits` now checks
12 anisotropic closed-form states at densities 0.01, 0.3 and 4, including both
parallel branches. `CheckIndependentEigenvalues` adds three independently
linearized primitive-CGL Jacobian references: 13.64355243397714,
2.6890819006579081 and 0.86379984723740366. Existing expected values unchanged.
Focused standalone regression: 2 passed (double and single precision).
Commit: `3d4716225`.

### T-A2: rates once, walls at split boundaries

Triage: merged code used two half-duration rate applications. The EOS now takes
an explicit `CGLCollisionMode`. The pre-LF and after-RK calls apply only walls;
the post-LF call applies full-cycle rates then walls, after restoring A. Pure CGL
applies the full operator after RK. Wall calls are independent of collision flags.
The post-sweep timestep refresh and the A/magnetic-moment representation ordering
are preserved. No-op projections leave conserved A untouched.

The rate/wall helper split needed by T-C2 is introduced here because A2 requires
it. Configurable thresholds and the extended limiter-map acceptance checks follow
in batch C. Existing hardwall configuration behavior is retained until then.

The new `cgl_collision_once.athinput` and CPU check verify analytic background
relaxation at every one of 20 cycles, with LF enabled and disabled, to 1e-10
relative. A finite-rate mirror test also checks the per-cycle backward-Euler map;
it distinguishes a full kick from two half kicks. Both new checks pass. The other
64 selected physical/workflow regressions pass; 19 campaign/provenance tests are
excluded from this focused run because their baseline is recorded separately.
The full validation workflow also passes all 19 cases, with its output bundle in
`/tmp/cgl-wo1-a2-full`. Changed production C++ and Python pass style checks; four
existing lint warnings remain in the superseded quantitative pgen.

Both decay-reference copies now use the independent continuous 2x2 temperature
amplitude system instead of reproducing the integrator schedule. At nu=10,
k=2*pi, t=0.02 and initial amplitude 1e-4 (before T-D3 coefficient changes), the
parallel-initial expected (Tparallel,Tperp) changes from
(7.65840395e-5,5.36770904e-6) to (7.65863284e-5,5.47261321e-6), and perpendicular-initial
from (1.11553911e-5,8.83615861e-5) to (1.09452264e-5,8.83592076e-5).
The older standalone quantitative copy still encoded doubled full-duration
collisions: its old pairs were (6.79840239e-5,9.76551990e-6) and
(2.02951004e-5,8.36901023e-5), respectively; both now use the same independent
continuous reference above.

### T-B1: pressure-floor consistency

Triage: the perpendicular-energy factor and pressure-floor A writeback were already
fixed. Repaired states still changed by roundoff on a second conversion. The
single-state conversion now returns primitives recovered from the repaired
conserved state and, only when necessary, rounds the corrected energy upward to
make the floor representable. The density-floor writeback is completed with B2.

The existing double/single precision checker now requires bitwise equality of all
conserved and primitive values after a second conversion, for each pressure-floor
branch, zero/threshold magnetic fields, and cancellation-dominated energy. It
checks the internal-energy identity to 1e-14 in double precision. Both precision
tests and the Release build pass. Existing expected physical pressures are
unchanged; the representability correction may raise a floored energy by a few
ULPs. One million additional deterministic scratch states in each precision
passed the admissibility and idempotence checks.

## Remaining task triage

All acceptance outcomes below remain unverified until the corresponding task.

| Tasks | Current code assessment and remaining work |
| --- | --- |
| B1 | Energy factor/A reset exist; literal bitwise floor idempotence needs work. |
| B2 | Density-preserving pressure ratio and density-floor A writeback missing. |
| B3 | CGL face pressure floors and nonfinite FOFC detection missing. |
| B4 | Both-side low-B isotropization and imported extreme A remain. |
| B5 | CGL bfloor positivity/single-precision validation missing. |
| B6 | CGL-aware prolongation exists; verify and omit obsolete fence. |
| B7 | Boolean parsing fixed; backup-only validation and sigma_max missing. |
| C1-C4 | Existing hardwall helpers/policies differ; apply specified parameterized monotone law and align all closure copies. |
| D1-D5 | Face normalization, limited gradients, 2nu perpendicular coefficient, collisional timestep and stiffness factor remain. |
| D6 | Post-RK parabolic reduction exists; add heating/stage-count acceptance. |
| E1 | Representation fence already exists and has regression coverage. |
| E2 | Passive speeds fixed; requested thermal-consistency fence remains. |
| E3 | Spitzer is implemented; verify it and fix missing-units handling. |
| F1 | Type 2 is rejected, so required random-mode support remains. |
| F2 | OU evolves once but force is RK-staged; implement requested single kick and update work diagnostics/primitive refresh. |
| F3 | Energy work exists but uses stale primitives and has secondary-fluid issues. |
| F4 | Generic planar normalization retains the wrong parallel axis. |
| F5 | Dead parameter survives in the legacy paper pgen. |
| G0 | Several files are duplicates; legacy nonlinear wave setups are unique and must be retained. |
| G1-G4 | References/assertions, long wave runs, actual convergence gates and limiter-map tests need strengthening. |
| G5-G6 | Generalize rotated-decay projection to volume-weighted multi-block/MPI and add two-level SMR decay/conservation. |
| G7 | Update current and retained legacy documentation to match demonstrated coverage. |
| P1 | Variable-restricted communication and frozen-B BC skip remain. |
| P2 | Two-variable copies/update exist; flux clear remains. |
| P3 | Reduced refresh region exists; fused temperatures and frozen-B caching remain. |
| P4-P5 | Flat LF kernel and conditional allocation exist; require bitwise/timing verification. |
| P6 | Both side logarithms still evaluated; retain selected-side arithmetic order. |

## Accepted T-D5 correction

The user corrected the inconsistent isothermal/pressure-equilibrium setup. Use
LF-only kinematic/advect evolution, uniform temperature, a normal magnetic field,
one-face density jumps at contrasts 10, 50, 200 and 1000, and paired seeded/unseeded
runs with a 1e-6 temperature seed. Require less than approximately 10-fold growth
of the difference over three cycles at sweep/parabolic-step ratio at least 10;
report stage counts. The old 20-cycle positivity criterion is superseded.

A retained pre-D scratch control already reproduces the instability: peak growth
over three cycles is 1, 2.33e5, 2.17e6 and 7.90e6 for the contrast ladder,
respectively, at seven stages per half-sweep. Positivity alone misses the R=50
failure. This is preliminary control evidence, not D5 completion.
