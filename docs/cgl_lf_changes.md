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
Multilevel runs refresh their coarse representation after each rate/wall call;
the existing 3D AMR coarse-anisotropy check caught a stale-cache regression when
this refresh was initially missing. That repair is folded into the A2 change.

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

### T-B2: preserve anisotropy across density floors

Triage: still present. `SingleC2P_CGLMHD` now recovers the pressure ratio using
the original density, computes pressures from the internal energy after the
density change, and resets A at the new density. `ConsToPrim` writes that A back
with the density. Nonfinite or overflowing anisotropy/ratio logarithms fall back
to isotropy, set a floor flag and repair A without exponentiating an overflow.

The double/float checker covers four sub-floor densities, three pressure ratios,
nonzero momentum, invalid densities, NaN/infinite A and both exponential and
ratio overflow. Ratios are preserved to 1e-12 in double precision, and repaired
states remain bitwise repeatable. Both precision checks and the Release build
pass. The deliberate overflow-policy change replaces the +1000 log-ratio
reference `(ppar,pperp)=(1e-12,3)` with `(2,2)`; the -1000 reference remains
`(6,1e-12)`.

### T-B3: reconstructed pressures and nonfinite FOFC detection

Triage: still present. All six PPMX/WENOZ directional wrappers now floor both CGL
pressures at pfloor while retaining ideal-MHD internal-energy floors. The FOFC
probe records nonfinite incoming internal energy/A before C2P can repair it and
also checks both recovered pressures for finiteness and floor compliance.

The built-in FOFC pgen checks all six reconstruction paths against fixed
pressure-floor values, retains ideal-MHD control cases, and injects NaN and both
infinities into energy and A to exercise the real detector. New PPMX/WENOZ inputs
evolve a three-cell 1000:1 perpendicular-pressure peak for 50 cycles; both remain
finite and positive with zero FOFC events. The same profile also remains finite
with the pre-B3 reconstruction/detector, so the proposed end-to-end NaN negative
control was not reproduced. Direct pre-fix floor and nonfinite-detection controls
do fail, confirming the assertions detect the actual defects.

Release build, C++/Python style checks, all 68 selected physical/workflow CPU
regressions, and the expanded 21-case full workflow pass. The latter bundle is
`/tmp/cgl-wo1-b3-full`. Existing physical expected values are unchanged.

### T-B4: weak-field anisotropy transport

Triage: still present. HLLE now isotropizes only the weak-field side, avoids
division by zero there, and uses the magnetized neighbor's field for weak-upwind
A transport. LLF uses the same reference field. The LLF reference copy is aligned.

The built-in FOFC test exercises real HLLE/LLF fluxes in a downstream finite-volume
cell update, including momentum, total energy and induction, followed by real C2P.
With bfloor=1e-10 and an incoming weak-field mass fraction of 0.1, the old
downstream pressure ratios 729 and 708.157 become 0.729 and 0.708157, inside
[0.5,2], without floors. Reverse flow, magnetized-side anisotropy preservation,
zero/threshold B and both-weak states are checked too. This is a single physical
cell update rather than a many-cycle grid run. Release build, changed C++ style,
and the three built-in FOFC/reconstruction regressions pass.

### T-B5: valid CGL magnetic-field floor

Triage: still present. The constructor rejects nonpositive bfloor and, in float
builds, bfloor cubed below FLT_MIN. The latter message requests an explicit bfloor.
The new focused constructor regression links the real EOS and parameter parser;
it checks zero, negative, default, underflowing and valid floors in both precisions.
Both precision tests, Release build and style checks pass.

A full single-precision build is independently blocked by existing narrowing
errors in `coordinates/cartesian_ks.hpp` and a mixed float/double `std::min` in
`diffusion/hyperviscosity.cpp`. The focused constructor and C2P checks compile the
real affected production sources in float without changing those unrelated files.

### T-B6: primitive prolongation already supported

Triage: already fixed by the merge. `prolong_prims.cpp` uses the CGL projection
helpers and handles both pressures and the A/magnetic-moment representations.
No fence is added. Existing smooth/uniform and strong-anisotropy refinement tests
and the primitive-versus-conserved stage-order test pass on the CPU binary:
`test_cgl_amr_gpu.py -k 'primitive_smooth or primitive_stage_order'` (2 passed).
They exercise actual refinement, require zero AMR repairs, and compare the
no-transfer paths. GPU and MPI execution remain to be checked separately.

### T-B7: limiter parsing and magnetization density floor

Triage: boolean parsing/defaults were already correct. The constructor now reads
backup_limiters independently and rejects backup-only configurations. Enabled
limiters require an explicit limiter_nu_coll; three hardwall AMR inputs now state
1e10 explicitly (their existing hardwall path still ignores this rate until C2).
Both A and magnetic-moment C2P enforce density >= max(dfloor,B^2/sigma_max).

Double/float constructor tests verify false flags, passive parsing, required rates
and backup-only rejection. Double/float C2P tests verify sigma_max with momentum,
energy preservation and the appropriate A/moment invariants. All four focused
tests, the Release build and style checks pass. The combined physical/AMR suite
had 78 passes and one coarse-anisotropy diagnostic failure (crs_d_err), traced to
the A2 post-wall coarse-cache refresh and repaired there. No physical expected
values were changed.

### T-C1: configurable thresholds

Triage: hard-coded thresholds required migration. EOS_Data now holds positive
firehose/mirror magnetic-pressure coefficients (defaults 2/1), backup factors
(defaults 1/2 for firehose/mirror), and backup LF suppression rate (default 1e10
in inverse code time). Shared predicates, limiter calls, AMR projection, LF
diagnostics, paper diagnostics and reference copies use those values. Parameters
are validated; legacy oblique/parallel aliases map to 1.4/2, and conflicts with an
explicit numeric value are rejected. An absent policy now defaults to 2 rather
than the former oblique value 1.4. Workflow provenance records numeric thresholds
and the A2 once-per-cycle schedule.

The finite-rate firehose stress fixture explicitly uses backup factor 2, keeping
its initial -0.85 B^2 state between the -0.7 B^2 soft threshold and the -B^2 wall.
This is necessary because the new factor-1 default would make that initial state
violate its enabled backup wall. Direct helper checks cover default/custom
thresholds and clipped walls. Constructor rejection tests use staged input files,
and numeric-versus-legacy runs produce identical histories. 79 selected physical
CPU regressions and all 21 full-workflow cases pass. Double/float closure and
constructor checks and the AMR projection unit check pass. All numerical decay
reference values remain unchanged in this task.
The repaired A2 coarse refresh passes all 11 AMR integration tests on CPU; the
independent speed and passive-signal checks also pass. C1 retains the existing
enabled-limiter diagnostic scope; C2 completes unconditional fluid-wall handling
through AMR transfers and its corresponding diagnostics.

## Remaining task triage

All acceptance outcomes below remain unverified until the corresponding task.

| Tasks | Current code assessment and remaining work |
| --- | --- |
| C2-C4 | Complete the monotone law, additive LF suppression and explicit input migration. |
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
