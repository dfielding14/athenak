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

Completed and committed through T-G4 (32 of 41 numbered tasks). Candidate checks
below are scratch evidence until integrated into the branch.

| Tasks | Current code assessment and remaining work |
| --- | --- |
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

### T-C2: monotone finite-rate relaxation and hard-wall projection

Triage: the merged implementation mixed an algebraic soft-threshold projection
with limiter relaxation. It now applies exact background decay, followed by the
backward-Euler soft update and the configured hard walls. The fluid firehose wall
is unconditional. LF strictness no longer enables backup limiters implicitly.
Legacy `limiter_hardwall=true` fails with migration instructions; 28 shipped true
settings were removed while preserving their explicit finite rate of 1e10.
Historical campaign contracts remain historical records.

`SingleCollRates_CGLMHD` uses the stable threshold-plus-residual expression.
`SingleCollWalls_CGLMHD` rounds a projected pressure inward when necessary, making
subsequent wall calls exact no-ops. The conserved-A encoder similarly chooses an
admissible representable value toward isotropy; its fallback bisection is used
only if a single ULP is insufficient. Energy is not changed. Primitive recovery
and AMR no longer impose the obsolete soft hardwall; AMR retains positivity,
weak-field, fluid-wall and configured-backup constraints.

The double/single EOS tests cover the analytic map, monotonicity at rate-times-step
0, 1 and 1e10, continuity across the former emergency boundary, positivity,
pressure conservation, exact wall idempotence and encoded-A admissibility on
10,000 deterministic random states. Larger scratch probes used one million wall
states and 100,000 AMR states per precision. The pressure-floor/AMR standalone
suite passes 3 tests; heat-flux policy, constructor and input-matrix checks pass 5.

Changed AMR expectations follow the new physics: the strong-anisotropy fixture
requires only fluid-wall repairs; the uniform firehose state has final
$|A|=\log(22/7)$; the mirror-slope fixture has no AMR soft-threshold projection.
The coarse diagnostic independently clamps the child mean to the coarse hard-wall
interval before comparing at the original tolerance. The LF 3D churn fixture's
pressure amplitude changes from 0.5 to 0.49: the former produces a real
$\Delta p+B^2=-0.005118$ crossing just after derefinement, now caught by the
unconditional diagnostic. At 0.49 it remains admissible without repairs. No
per-stage limiter or diagnostic tolerance was added.

Final C2 integration: all 11 AMR tests, 79 selected CPU regressions and all 21
full-workflow cases pass. Changed C++ passes cpplint. No tolerance was relaxed.

The final shearing regression exposed another fixture outside the new strict
fluid-wall envelope: its collisionless background fails the first post-RK LF
stage with 1016 hard-bound violations on both post-G and final binaries.
Changing only `tst/inputs/cgl_lf_sts_sbox.athinput` background `nu_coll=0` to 30
keeps its nonzero perturbation, field, strict mode, all tolerances, and analytic
magnetic references. The full serial/MPI/explicit/restart test passes on both
binaries, with all 30 output files byte-identical. It completes 38 STS cycles
to t=0.3 and 106 capped/explicit cycles to t=0.04, with no repairs or hard-bound
violations. This is a collisional boundary/restart test; the original
collisionless strict fixture is not claimed to be supported. Evidence:
`/tmp/cgl-wo1-sbox-final/summary.txt`.

### T-C3: additive LF suppression rates

Triage: the shared limiter helper selected a maximum/replacement rate. It now adds
soft and configured backup rates, with background collisions added exactly once
by callers. Threshold equality is admissible and contributes no scattering.
The existing shared calls update both production transport and retained pgen
references. The snapshot transport proxy now uses the same numeric parameters,
configured backup policy and additive rates.

Independent prescribed-rate cases use background 7, soft 11 and backup 101:
totals are 7, 18, 108 and 119, including a fluid-wall-only activation and scaled
magnetic fields. Production safe/fast face fluxes are checked against independently
prescribed totals 18 and 119. The proxy has a separate eight-state rate check.
Standalone heat-flux tests pass in both precisions; all 84 selected CPU regressions
pass. Existing expected helper rates change from backup alone to backup plus soft
(1e10 to 1e10+20 in the fixture, and configured 1234 to 1254).

### T-C4: explicit input thresholds and synchronized references

All 29 current paper inputs explicitly state the five numeric threshold/backup
parameters. Legacy aliases are removed from those inputs. The firehose stress
fixture explicitly retains threshold 1.4 and backup factor 2; the mirror fixture
uses thresholds 2/1. Workflow overrides and parameter examples use numeric keys.
The distinct legacy-policy fixture remains for compatibility checks. References
already use the shared EOS thresholds from C1/C3.

All 85 selected CPU regressions pass, including literal-key checks of every paper
input. A paper-smoke workflow probe passed the active and passive Alfvénic cases;
its random case exposed a pre-existing incompatible driving-type/projection
combination in the workflow override. This forcing configuration is addressed in
batch F. No physical expected value changes in C4.

### T-D1: CT-normal face field and arithmetic cell-magnitude normalization

Triage: still present. `AddHeatFluxes` now receives the frozen CT face field.
All 12 normal, diagnostic and profile face-state paths use its normal component,
averaged transverse components and $\bar B=(|B_L|+|B_R|)/2$. The direction vector
is not renormalized. Inverse-field factors, weak-field cutoff and limiter magnetic
pressure use the same $\bar B$. Both retained grad-B reference copies were updated;
`FaceCParallel` had no magnetic-field averaging to change.

New 1D/2D field-reversal inputs run the LF sweep with the existing kinematic
advection harness, seed 1e-8, guide field 0.03, sweep ratio 10 and 20 cycles.
All eight safe/fast and full/none cases pass, including profile probe paths, with
peak growth at most 0.9982 (1D) and 0.9991 (2D), positive pressures and no repairs.
The 1D negative control reaches 45.27 times the seed after its first half-sweep
and 1.56e6 by the third pre-sweep (peak pressure perturbation 0.0161188).
The specified 2D negative control remains stable at its smaller dimensional dt;
we do not reproduce the review's claimed 1e5 growth per sweep.

All 93 selected CPU checks and 23 full-workflow cases pass. Existing smooth
references/tolerances were not retuned; grad-B RMS error is 1.0126%. Changed C++
and Python pass style checks.

### T-D2: limit transverse temperature slopes

Triage: the 48 transverse temperature expressions still used arithmetic averages.
Every transport, diagnostic and profile path now uses the existing four-slope
van Leer mean. Magnetic-magnitude gradients remain byte-for-byte unchanged.
The existing limiter definitions were moved to a shared diffusion header so
Spitzer and LF use the same implementation.

The new 2D Gaussian hotspot has 100x contrast and runs for 10.02 perpendicular
diffusion times at 30 and 45 degrees. All safe/fast, full/none combinations,
including profile probes, preserve the initial minimum to 5.9e-14 and total energy
to 2.7e-14. Before the fix, parallel/perpendicular minima were 0.57863/0.69414 at
30 degrees and 0.50099/0.60814 at 45 degrees, from initial minima 1. Existing 3D
x/y/z/oblique smooth-decay checks retain their original tolerances (oblique error
0.002405 versus tolerance 0.03). The full 24-case workflow passes.
The eight new Python hotspot regressions pass; the existing 93 selected checks
also pass. The test reads the history writer's truncated `max_energy` label.

A subsequent independent audit found overflow/underflow in the reused van Leer
mean for finite extreme slopes. Its ordinary evaluation order is preserved;
only unusable products/sums use a scaled harmonic mean. The actual header now
passes direct double/single checks for equal maximum-scale, minimum-normal and
subnormal slopes, opposite slopes and extrema. Before this correction several
finite equal-slope means became zero or nonfinite. Both precision regressions
and the integrated Release build pass.

### T-D3: perpendicular BGK collision coefficient

Triage: the denominator still contained one collision frequency. It now contains
$2\nu_{\rm eff}$ in the scaled production ratio and its logarithmic fallback,
both retained quantitative references, the prescribed-rate face-flux reference,
legacy paper diagnostics, snapshot reconstruction and documentation. The code
records the deliberate SHD97/Sharma choice versus Squire et al. eq. 2.7.
The collisionless eigenmode generator requires no change.

The independent strong-collision check verifies $\chi_\perp\nu/c_\parallel^2\to1$.
Double/single closure tests pass, as do 93 selected CPU checks and all 24 workflow
cases. The previously verified eight collisionless hotspot cases are unchanged.
At c=1, nu=10, k=2*pi, t=0.02 and initial amplitude 1e-4, chi_perp changes from
0.07767107945385573 to 0.05594466633444508. The continuous parallel-initial
(Tparallel,Tperp) reference changes from (7.65863283905e-5,5.47261320692e-6) to
(7.65884217595e-5,5.52095332130e-6); the perpendicular-initial pair changes from
(1.09452264138e-5,8.83592076405e-5) to (1.10419066426e-5,8.98859289040e-5).
No tolerance changes were needed.

### T-D4: include background collisions in the stability bound

Triage: the LF timestep still used collisionless chi_parallel. It now includes
background nu, excluding limiter rates as required. The usual arithmetic order
is retained for nu=0; the existing scaled closure handles otherwise-overflowing
intermediates. The fac and cfl multipliers are unchanged.

Five independent Decimal checks pass: nu=0/10/100 and two finite-diffusivity
extremes with overflowing intermediate products. In the ordinary uniform fixture,
allowed timesteps are 0.0615219138504782, 0.08934960839675515 and
0.3397988593132476, respectively; before the fix all three used the first value.
The Release build, 98 selected CPU checks and 25-case full workflow pass. The full
workflow now requests validation CSV output only when the input declares that
parameter, allowing monitoring-only acceptance fixtures without invalid overrides.

### T-D5: conservative density and magnetic stiffness factors

Triage: the LF timestep lacked the face-to-cell stiffness factors. `NewTimeStep`
now receives fresh cell-centered B, scans both faces in each active direction,
and multiplies chi by the maximum arithmetic-face-density ratio times
$\max(1,|B_i|/\bar B_f)$. It does not reuse a stale pre-RK magnitude cache.
The dimensional fac and cfl multiplier are unchanged.

The corrected LF-only density-contact acceptance passes in safe and fast modes
at contrasts 10, 50, 200 and 1000, with 1e-6 seeded/unseeded pairs, three cycles,
sweep ratio 10 and seven stages per half-sweep (2688 cell-stages). Peak difference
growth is at most its initial value; final growth factors are approximately
0.2601, 0.6662, 0.8898 and 0.9745. All repair/nonfinite counters remain zero.
The retained pre-D control grows by 2.33e5 at R=50, 2.17e6 at R=200 and 7.90e6
at R=1000. A separate pre-D5 control with D1 already fixed also fails at R=50.

Independent initial-dt checks match both the density factor $(R+1)/2$ and the
reversal-sheet magnetic factor 1.32333161698 (1D/2D). All five D4 uniform
timesteps are bitwise identical with and without D5. Integrated validation:
106 selected CPU checks and all 26 workflow cases pass; tolerances unchanged.

The smaller D5 step shifts the AMR churn fixture's time-triggered derefinement
from cycle 3 to cycle 4, adding hyperbolic evolution before the LF restart. Its
switch time is adjusted from 2.70e-4 to 2.55e-4 to retain the original three-cycle
transition; pressure amplitude and all strict/conservation checks stay unchanged.
All 11 AMR integration tests pass again with D1-D5 integrated.

### T-D6: refresh the post-sweep stability bound

Triage: already fixed. `Mesh::RefreshSTSParabolicTimeStep` performs the local
minimum and MPI reduction immediately before `Driver::Execute` begins the post
sweep. No production change is needed. The new `timestep_refresh` pgen mode and
input heat a uniform state during RK and independently check its analytic
pressure and timestep, with rank-wide minimum/maximum checks of both stage counts.

Unheated runs retain seven stages in both sweeps (896 cell-stages); heating at
3000 raises the post sweep to eleven (1152 cell-stages), final pressure
8.6902392313 and CFL-scaled LF timestep 6.521748e-5. Serial and one/four-rank MPI
results agree to 2e-12. Removing the existing refresh is a failing negative
control: the heated post sweep incorrectly retains seven stages. Both new CPU
cases, all three MPI regressions, all 11 AMR regressions and the 27-case full
workflow pass. Existing expected values and tolerances are unchanged.

### T-E1: representation-aware boundary guard

Triage: already fixed. The `MHD` constructor rejects inflow/user boundaries for
both explicit and STS LF split integration because IAN temporarily stores
magnetic moment. The message now points to WO2; the existing policy is unchanged.
The regression now includes periodic alongside outflow, reflect and diode.
All ten allowed/rejected boundary cases pass. No reference values change.

### T-E2: disable the inconsistent passive thermal model

Triage: passive HLLE already uses the isothermal signal speed. `CGLMHD` now
rejects `passive=true` before reading its sound speed, explaining review M7 and
the WO2 redesign. Direct unit coverage of the dormant isothermal flux/speed path
and active-CGL CFL remains. Runtime passive tests now check the explicit fence.
Executable workflows omit passive cases and record their reason in the manifest;
historical catalogs and inputs are preserved as disabled references. The legacy
smoke script omits the passive input, and current/retained runbooks mark it disabled.

The Release build, 110 selected CPU regressions and three independent constructor
(double/single precision) and workflow-filter checks pass. No numerical reference
values change; passive runtime evolution is deliberately unavailable until WO2.

### T-E3: validate implemented Spitzer conduction and require units

Triage: the pre-merge unimplemented-flux finding is obsolete. The current
`Conduction::AddHeatFlux` calls the live temperature-dependent flux in all active
directions, and `NewTimeStep` computes its bound. Saturated heat flux already
rejects STS. No obsolete "not implemented" fence is added.

The constructor now requires a units object for Spitzer conductivity before any
conversion can dereference it. Two new regressions check the missing-units error
and a finite, positive, nonzero thermal update with units, for both unsaturated
and saturated explicit Spitzer conduction. All six STS diffusion regressions
pass. Existing reference values and implemented conduction behavior are unchanged.

### T-F1: restore type-2 random forcing

Triage: the merged driver rejected type 2 outright. Its shared `IsDrivenMode`
now includes type 2 in the isotropic shell, and amplitude construction uses the
same isotropic spectrum. Type 2 defaults to the retained power-law, unprojected
random policy; incompatible explicit projection policies fail clearly. Both mode
initializers verify the selected count against allocation. The pre-existing
zero-mode exclusions remain, including for nlow=0.

This affects all driven fluids. Three focused type-2/invalid-policy checks pass,
and the current paper-smoke workflow now passes both enabled cases. Its random
case explicitly selects the compatible random projection; the fenced passive
case is recorded as disabled. Nonzero force is verified here; injected power is
verified with the once-per-cycle schedule in F2. No reference values change.

### T-F2: one full forcing kick before RK

Triage: the merged driver advanced OU coefficients only at their configured
update boundaries, but applied the fluid kick in each RK stage. For every driven
fluid, `IncludeInitializeModesTask` now schedules a single full-dt `AddForcing`
after normalization; the stage registration and RK work recurrence are removed.
`Driver::Execute` reuses its boundary/primitive refresh before `CopyCons` and the
first hyperbolic flux. Initial-cycle forcing is normalized as well.

The shared driver affects hydro, MHD and CGL. Six power tests cover driving types
0/2 and RK1/2/3; every positive-duration step injects dedt=0.1, including the first,
with maximum relative error 5.13e-8. The prior RK2/RK3 schedule gave approximately
0.175/0.16518 after an uninjected first cycle. Three CGL multi-cycle cases verify
that accumulated forcing work matches the measured energy gain. All 126 selected
CGL/turbulence checks pass; changed C++ passes style checks. Modal tcorr/dt_update
behavior is unchanged and receives a direct recurrence check with F3.

### T-F3: exact kinetic work from each fluid's conserved state

Triage: forcing already added energy, but used potentially stale primitive
velocities; the two-fluid branch also reused the primary fluid's velocity and
unconditionally wrote the secondary energy slot. `ApplyForcingWithStep` now
computes each nonrelativistic fluid's exact kinetic-energy increment from its own
pre-kick conserved momentum and density, guarded by its own `is_ideal` flag.
Nonrelativistic net-momentum removal already acts on acceleration before the
kick, so this work includes that correction without a separate energy update.
The existing relativistic source remains algebraically unchanged.

A direct built-in `turb_forcing` pgen deliberately supplies stale primitives,
unequal two-fluid densities/velocities, and isothermal passive-scalar sentinels.
All nine ideal-MHD/CGL/hydro/isothermal and two-fluid combinations pass, with
maximum thermal/scalar/impulse error 1.78e-15. The previous implementation fails
this check. The same fixture verifies nonzero modal amplitudes and the exact OU
hold/update recurrence using tcorr and dt_update. All 23 current turbulence
regressions pass. These shared-driver repairs apply to every driven fluid;
existing physical reference values are unchanged.

### T-F4: consistent parallel axis for planar forcing

Triage: the paper-specific planar policy already used z as its parallel axis,
but the generic policy still used x. `InitializeModes` now uses
$k_\parallel=|k_z|$ and $k_\perp=\sqrt{k_x^2+k_y^2}$ for every type-1 policy,
consistent with its selected modes and force components. This affects all fluids
using the generic planar driver; existing paper-specific behavior is unchanged.
All 25 turbulence regressions pass, including direct nonzero-power checks for
every selected planar and type-2 mode. The planar check fails before the fix and
passes with all 30 selected modes nonzero. No other reference values change.

### T-F5: remove the unused forcing selector

Triage: the built-in paper pgen already has no selector; the retained legacy pgen
still read and validated a parameter that never affected forcing. Its parser,
three input keys and smoke-script override are removed. Both runbooks direct
users to the actual turbulence-driver controls; inferred analysis metadata is
retained. No physics or numerical reference changes.

Building the retained generator exposed two pre-merge integration problems:
its history diagnostic used the removed `Conduction` LF API, and the always-built
turbulence initializer emitted a duplicate `UserProblem` for every custom pgen.
The diagnostic now reads `pmhd->pcgl_lf`; CMake enables the turbulence wrapper
only for `PROBLEM=turb`. These minimal compatibility fixes are necessary to
verify the retained smoke script. Its obsolete analyzer flag/directory argument
is also replaced by the current history arguments and synthetic check.

The default and legacy paper builds pass. All seven enabled legacy smoke cases
run, all seven histories are finite, and the analyzer synthetic check passes.
The output is `analysis/diagnostics.json`; it contains basic finite/time summaries
for the retained legacy history labels, rather than the obsolete summary.json.
Two untouched legacy C++ line-length warnings remain.

### T-G0: remove only superseded local files

All 18 listed files were compared with current equivalents. Five are removed:

| Removed file | Retained replacement |
| --- | --- |
| `src/pgen/unit_tests/cgl_lf_quantitative_test.cpp` | `src/pgen/tests/cgl_landau_fluid.cpp`, covering every old mode plus the new acceptance cases |
| `src/pgen/unit_tests/cgl_fofc_end_to_end_test.cpp` | `src/pgen/tests/cgl_fofc.cpp`, retaining the mutation check and adding reconstruction/NaN/weak-field checks |
| `src/pgen/diffusion_test.cpp` | Built-in CGL LF parallel/perpendicular/grad-B modes, with quantitative assertions |
| `inputs/unit_tests/README-cgl-lf.md` | Current Sphinx CGL LF module and validation guide |
| `docs/cgl_lf_validation.pdf` | `docs/source/_static/cgl_lf_validation.pdf`, the expanded report built from the retained TeX |

The other thirteen remain because they carry unique content: both legacy plan/
runbook documents; all seven legacy paper inputs (four distinct nonlinear wave
setups, active calibration, limiter-off control and disabled passive reference);
the smoke script exercising those inputs; `src/pgen/cgl_lf_paper.cpp` implementing
the four wave modes absent from the built-in turbulence pgen; `diffusion_2d.cpp`
with its two-dimensional perpendicular-temperature setup; and `lpaw_paniso.cpp`
with anisotropic linearly polarized Alfvén initial data. Retained links now name
the built-in quantitative pgen. No numerical state or expected value changes.
The default target does not compile the removed custom pgens; its Release build
and all 27 current full-workflow cases pass after removal.

### T-G1: independent references and synchronized diagnostics

The final copy audit confirms all A-D closure coefficients, configured additive
rates and grad-B face normalization are synchronized. Uniform collisions use an
exact exponential; collisional temperature amplitudes use the independent
continuous 2x2 system; fast speeds use closed limits and offline Jacobian
eigenvalues. Limiter rate sharing is checked by fixed literal totals and fluxes
with independently prescribed frequencies. G2-G6 add time, convergence and
multi-block checks around these references.

The old expanded speed copy is removed from `cgl_fofc_flux_test.cpp`: that wiring
test calls the actual EOS speed, whose independent oracle is tested separately.
Its unequal active/passive flux and equal-state EMF checks pass. The snapshot
heat-flux proxy now reads canonical numeric firehose metadata before falling
back to old archived aliases; previously `model_choices` emitted legacy='none'
and valid new inputs failed analysis. All 11 mechanism-analysis tests pass,
including numeric 2/1.4 and archived parallel/oblique policies at a fixed physical
state. No expected physical values change. The retained legacy cell-centered
heat-flux diagnostics remain proxies, not independent face-flux oracles.

### T-G2: wave tests reject ineffective evolution and require convergence

All nine eigenmode/oblique/field-wave inputs and eigenmode-generator defaults
now use a 1e-3 short-run tolerance, replacing 7.5%, 40% or 25% guards. Each pgen
checker requires the reference evolution to exceed ten tolerance units, so a
frozen or nearly zero-duration state cannot pass. The former initialization-only
multi-block and MPI wave tests now evolve through the normal short duration.
The static-AMR wave already evolved six cycles; its stale 0.75 tolerance is
replaced by 0.001. Its maximum error is 6.05e-6; evolved MPI one/four-rank
errors agree (maximum 2.484e-6). Other zero-step checks measure constructor,
initial-timestep or fixed-flux behavior and remain unchanged.

The new acceptance suite runs each family at 64/128/256 cells to 1.25 periods,
requiring order at least 1.8 for pure/Alfvén waves and at least 1 for LF-dominated
waves. Measured orders are 1.95-2.49, with finest-grid errors at most 5.2e-4.
Three too-short cases and two disabled-LF eigenmodes are rejected. All 23 new
checks and all 27 full-workflow cases pass. Analytic eigenvalues and amplitudes
are unchanged; only test times, resolution coverage and tolerances change.
Optional collisional/entropy Figure 17 eigenbranches remain outside this task.

### T-G3: resolve one e-folding of parallel and perpendicular decay

The four quantitative inputs now use 128 cells, STS ratio 10 and hundreds of
cycles to $\chi k^2t=1$, with amplitude tolerance 0.003 and phase tolerance 1e-6.
An independent Python matrix exponential verifies amplitudes and proves that
swapped diffusivities, half/double collisional perpendicular diffusivity and the
old +nu denominator lie more than ten tolerance units away. Measured errors are
0.000208-0.000271 over 208/208/416/649 cycles; all four new and two existing decay
checks pass.

For initial amplitude 1e-4, collisionless parallel/perpendicular expected final
amplitudes are both 3.678794411714424e-5. The collisional parallel-initial pair is
(1.9221708099735522e-5,1.3784043081777392e-5) at t=0.14484812929181745;
the perpendicular-initial pair is (1.169172828934281e-5,1.7978744606824823e-5)
at t=0.45277409930656082. These replace the t=0.02 values recorded in D3 because
the acceptance duration changes; the reference equations do not.

All 33 accuracy-study cases also pass. Its deliberately coarse 32/64-cell
samples use dx-squared error budgets 0.048/0.012, while 128 cells retains 0.003.
The ratio-1000 few-step probe uses 0.04 (measured maximum error 0.003779).
Those convergence/cost probes are distinct from the strengthened ratio-10
acceptance and do not relax its criteria.

### T-G4: per-cell limiter relaxation oracle

`CheckLimiterStress` now independently reconstructs each initial cell, applies
exact background relaxation, the finite backward-Euler soft map and configured
walls, including the LF pre-wall ordering. It checks final anisotropy and
conserved isotropic pressure to 1e-12. This analytic comparison is used for
uniform LF states or nonuniform, zero-advection pure-CGL states, where spatial
transport cannot invalidate a cellwise collision oracle. The existing nonuniform
LF stress remains an admissibility test of combined transport and limiting.

All 28 new cases pass: mirror/firehose, pure/LF, limiter nu*dt=0/1/1e10,
background nu=0/3, and four initially out-of-wall cases verifying the ordering.
No production formula or existing reference values change.

The one-cycle cellwise oracle is gated to its stated cycle/solver assumptions;
the existing 20-cycle uniform limiter regression remains supported and checks
the independent backward-Euler history at every cycle. The final combined run
caught and repaired an initially overbroad G4 guard that rejected that caller.
