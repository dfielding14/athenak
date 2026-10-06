# WO1 CGL-LF changes and validation

## WO2 follow-ups on Frontier

WO2 starts from `9a4030b09307bc43865d5e597638f8645a6388f8`. The
[implementation and validation report](validation/wo2/README.md) records the
actual CPU/HIP/MPI baseline, retained artifacts, and per-task evidence.

Task 0 adds the verified CCE20/HIP build contract and refreshes physical CGL-LF
corner ghosts after prolongation. This corrects a demonstrated baseline error
at an outflow/coarse-fine intersection. Magnetic boundary updates remain
required. Compact communication is evaluated separately in Task 3.
Conservative A and both accepted strict B4 expected failures are retained.

The following sections preserve the WO1 implementation record.


Status: the 41-task implementation is complete. On 2026-10-06 the user accepted
the T-B4 sharp-contact failure as a known limitation shared with the checked
reference implementations; retain A and the current implementation. This is
an [accepted limitation](validation/wo1/review/weak-field/README.md), not a passed
validation, and no longer blocks WO1 by itself. T-P1 is deferred to WO2 under
the required bitwise rule; its candidate is not retained.
The [closeout](#wo1-closeout-2026-10-06) records the passing CPU/MPI/CUDA compile
checks and the limits of that validation.
Base: `8222de3aa`; branch: `c/cgl-lf-wo1`.
The work order controls the design; `CGL_LF_STS_review.md` describes pre-merge
`2aee609` and is used for derivations rather than current bug status.

## User-visible changes

- Turbulent forcing now applies one full-step kick before RK for every driven
  fluid, including hydro and MHD. Existing driven results can change. Each fluid's
  energy update uses its own conserved momentum; OU coefficient cadence is unchanged.
- Generic planar forcing now uses z as its parallel axis. The paper-specific
  policy already used z.
- The CMake turbulence `UserProblem` wrapper is enabled only for `PROBLEM=turb`,
  allowing other custom problem generators to link without a duplicate definition.
- An unspecified firehose threshold now defaults to $\Lambda_{\rm FH}=2$.
  Set the numeric threshold explicitly to retain another value.
- `limiter_hardwall=true` now fails at startup. Use finite-rate soft relaxation
  and the configured backup walls; the fluid firehose wall is always enforced.
- Passive CGL mode is disabled pending its thermal-model redesign. CGL-LF
  inflow/user boundaries remain rejected while IAN stores magnetic moment.
- LF transverse temperature slopes now use the van Leer limiter shared with
  Spitzer through `diffusion/limiters.hpp`. The shared mean also handles extreme
  finite slopes without overflowing; magnetic-magnitude slopes are unchanged.
- Strict LF checks enforce hard walls at sweep entry and after the scheduled
  end projection. Intermediate crossings remain counted but are not fatal;
  floors and invalid states still fail immediately. Final post-AMR states also
  receive a walls-only projection, without repeating rates.

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

**Review follow-up: multi-cycle failure.** A pressure-balanced transverse
contact on 128 cells, with $B=10^{-12}$ entering $B=1$, `bfloor=1e-10`,
$v_x=\pm10$ and no collisions, passes the first update but exceeds
$p_\perp/p_\parallel=2$ by cycle 5. By cycle 50 the maximum ratio is about
$2.46\times10^{10}$, without EOS floors. The supplied weak-field flux change
therefore does not establish the requested multi-cycle bound. The new
`cgl_weak_field_transport.athinput` and `test_cgl_weak_field_cpu.py` preserve the
reproducer. The 2026-10-06 decision accepts this shared sharp-contact limitation
and retains A; no flux formula or test bound was changed to make the test pass.
See the [reference comparison and decision](validation/wo1/review/weak-field/README.md).

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

## Accepted T-D5 correction

The user corrected the inconsistent isothermal/pressure-equilibrium setup. Use
LF-only kinematic/advect evolution, uniform temperature, a normal magnetic field,
one-face density jumps at contrasts 10, 50, 200 and 1000, and paired seeded/unseeded
runs with a 1e-6 temperature seed. Require less than approximately 10-fold growth
of the difference over three cycles at sweep/parabolic-step ratio at least 10;
report stage counts. The old 20-cycle positivity criterion is superseded.

The retained pre-D control reproduces the instability: peak growth
over three cycles is 1, 2.33e5, 2.17e6 and 7.90e6 for the contrast ladder,
respectively, at seven stages per half-sweep. Positivity alone misses the R=50
failure. The corrected D5 results are recorded below.

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
interval before comparing at the original tolerance. Review follow-up restores
the LF 3D churn fixture's pressure amplitude from 0.49 to its original 0.5.
The completed derefinement/field refresh creates eight violating cells before
LF starts: the measured minimum $\Delta p+B^2$ changes from $+0.06065$ to
$-0.004616$. `AdaptiveMeshRefinement` now applies walls only after that final
field/primitive refresh and restricts the corrected coarse representation.
The new test checks the post-AMR state before LF can repair or evolve it;
total-energy drift is $2.49\times10^{-14}$ at the unchanged conservation tolerance.

Final C2 integration: all 11 AMR tests, 79 selected CPU regressions and all 21
full-workflow cases pass. Changed C++ passes cpplint. No tolerance was relaxed.

Review follow-up also restores the shearing fixture to `nu_coll=0`. At cycle 8,
$t=0.06697$, after-RK walls leave all active cells admissible, but the first
post-LF stage creates 1016 crossings with minimum margin $-9.35\times10^{-7}$.
This is intermediate LF transport, not a missing after-RK wall call. Strict
hard-wall checks now occur at sweep entry and after the prescribed end projection;
every intermediate crossing remains counted. Floors/nonpositive/nonfinite states
still fail at every stage. Safe and fast RKL2 runs reach $t=0.3$ in 39 cycles,
recording 1,172,382 hard-bound stage visits; the explicit reference also completes.
An isolated negative control omitting the end projection fails the exit check,
and invalid initial states still fail at entry. No per-stage projection is added.

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

### T-G5: global oblique-decay projection

Triage: the rotated-temperature projection required a single block. It now sums
cell volumes across all blocks and MPI ranks, then measures mean, sine and cosine
amplitudes. `ProjectRotatedTemperature` and its guards/CSV output in
`src/pgen/tests/cgl_landau_fluid.cpp` changed; the new 64-squared oblique input is
included in the full workflow. CPU and MPI regression checks compare one block,
a 2-by-2 block layout, and one/four ranks for both x and y temperature modes with
B=(sqrt(3)/2, 1/2, 0).

Release CPU/MPI builds pass. The two CPU tests and MPI test (six MPI runs) pass:
one e-folding takes 70/208 cycles, relative amplitude errors are 0.00085305 and
0.00083045 against a 0.005 tolerance, and no repair/admissibility counters fire.
Sampled state is bitwise equal between decompositions; summation-order mean
changes reach 2.78e-14, with Fourier-amplitude changes around 1e-19. The checks
allow bounded roundoff in global reductions. An isolated old-projection negative
control fails by approximately 100% on four ranks. Physical references unchanged.

### T-G6: two-level SMR decay and conservation

Triage: no resolved SMR LF decay acceptance case existed. Added the 40-block,
two-level `cgl_lf_smr_decay_2d.athinput`, CPU/MPI checks in the oblique-decay test
module and MPI suite, and a full-workflow entry. The test verifies that flux
communication/correction, restriction and prolongation execute, runs 70 cycles
for one e-folding, and compares the global amplitude with the continuum decay.

The CPU SMR test and one/four-rank MPI comparison pass: amplitude relative error
0.00205802 (tolerance 0.01), total-energy drift at most 7.11e-15 (tolerance 2e-13),
and all repair/admissibility counters zero. The full workflow passes all 29 cases.
The broader AMR/MPI run exposed two old initialization-only wave callers that
now correctly fail G2's evolution guard; their repair belongs to the G2 commit.
No production physics or reference values changed in G6.

### T-G7: current method documentation and bounded validation claims

Updated `docs/cgl_lf_validation.tex`, the retained reproduction plan/runbook,
and the Sphinx LF module, method primer, code guide, AMR and validation pages.
They now describe the full-cycle rate/wall schedule, numeric thresholds, the
BGK perpendicular coefficient, passive/boundary fences, forcing cadence and
exact kick work, and measured G2-G6 acceptance. The superseded unit-test README
was removed in G0; the current Sphinx validation guide replaces it. No current
claim rests on the archived May 24/25 figures: captions identify their dates,
and retained E02 campaign results are explicitly historical. OU coefficient
updates retain their configured cadence; only the forcing kick occurs each cycle.

The repository report compiles with its real figure assets to 18 pages; every
page was rendered and visually inspected, with no clipped text or figures.
The current PDF is `docs/source/_static/cgl_lf_validation.pdf`. Narrow path/digest
line-break repairs remove all overfull boxes. Sphinx now enables dollar-math
and MathJax, and the new acceptance anchor resolves. A full strict Sphinx build
has only two pre-existing orphan-engineering-document warnings; suppressing
only that warning class yields a successful build. No scientific figure data
or historical campaign measurements were regenerated. QA evidence is retained
in `/tmp/cgl-wo1-g7/evidence.md` and its build logs.

## Performance baseline and acceptance before review follow-up

The final post-G numerical binaries were frozen before P: Release double,
Kokkos Serial, Apple M4 Max ARM64, Clang `-O3`, with a separate MPI build.
G7 only changes documentation. Twelve fixed inputs give 24 serial/one-rank/
four-rank configurations, each repeated three times. Every binary dump,
16-digit history and full-precision restart (including ghost cells) is compared
byte-for-byte. The baseline itself is repeatable with zero differences across
72 runs. Required cases are uniform 1D/2D/3D LF, two-level SMR LF and pure CGL;
additional cases cover shear, low B, density contrast, velocity, outflow, and
nonuniform primitive-prolongation SMR with a passive scalar, including an
interface touching an outflow boundary.

The fixed inputs, hashes, frozen executables, all output files and reproducible
harness are retained under `/tmp/cgl-wo1-perf/`. `harness.py run --label NAME
--serial BINARY --mpi MPI_BINARY --repeats 3 --mpi-ranks 1 4 --compare post-g`
performs the comparison without restaging inputs. Baseline source binaries
include G6; the later G2 caller correction changes no source or staged fixture.
`runs/post-g/provenance.json` records source/compiler/input/executable hashes.

Timings below are medians of three runs with identical output/profiling settings.
Cycle time is the solver's wall-clock timer. Stage time is the sum of existing
exclusive LF/STS compute timers divided by actual RKL stages. Shared transport
timers include RK/initialization, so they are recorded separately in JSON and
are not mislabeled STS-only time. Whole-STS wall time is not available separately;
no new production timing instrumentation was added. CPU timings are indicative,
not GPU benchmarks.

### T-P1: rejected non-bitwise communication optimization

Triage: ordinary parabolic CC exchange still sends all variables. A candidate
added dense offset/count packing to `MeshBoundaryValuesCC::PackAndSendCC` and
`RecvAndUnpackCC`, used the existing active `InitRecv(nvars)` message sizing,
and selected IEN/IAN in MHD LF stage tasks. Shared buffer capacity was retained
for hyperbolic reuse. It also skipped magnetic physical BCs during LF. Shearing
kept full ordinary exchanges because its remap consumes all freshly filled
ghosts; a complete compact shear API would require further work.

Both CPU/MPI builds and style checks passed. Twenty-one of 24 configurations
were byte-identical. The SMR/outflow-boundary case differed reproducibly in
serial and MPI one/four-rank runs: both final binary state files, full restart,
and both histories changed. By cycle 4 this includes active cells, not only
ghosts: maximum absolute energy difference 3.45e-3 and magnetic-component
differences about 2e-3. The total-energy history first differs at cycle 1 and
ends 1.56e-6 relative from the baseline. The entire candidate was removed under the work
order's strict bitwise rule and deferred to WO2. No P1 production change is
retained. The rejected diff and outputs are preserved in
`/tmp/cgl-wo1-p12/P1-tested-rejected.patch` and
`/tmp/cgl-wo1-perf/runs/p1/`. The combined coarse/fine and physical-boundary
provenance needs resolution before this optimization can be retried.

This is a WO2 **correctness investigation**, not only a postponed optimization.
The evidence does not establish whether the candidate's packing was wrong or
the baseline requires refreshed density, momentum, scalar or magnetic ghosts at
the coarse/fine outflow interface. Isolate the two candidate changes (compact
cell exchange and skipped magnetic BCs), then compare those ghost values before
and after boundary filling, prolongation and flux correction through the first
differing cycle. Require an explained dependency and a regression that detects
the active-cell discrepancy before attempting the optimization again. No part
of the rejected candidate is enabled. The durable evidence archive is
[`validation/wo1/`](validation/wo1/ARCHIVE.md).

Baseline and rejected-P1 timings (cycle ms / profiled compute microseconds per
stage; rejected timings are diagnostic, not a claimed speedup):

| Case | Post-G | Rejected P1 |
| --- | --- | --- |
| lf1d-serial-0 | 0.513 / 30.20 | 0.522 / 31.17 |
| lf2d-serial-0 | 17.949 / 1175.40 | 18.227 / 1199.71 |
| lf3d-serial-0 | 222.032 / 14492.02 | 223.315 / 14701.98 |
| smr-serial-0 | 17.653 / 927.78 | 17.171 / 939.89 |
| pure_cgl-serial-0 | 0.414 / n/a | 0.439 / n/a |
| lf2d-mpi-4 | 5.315 / 311.00 | 5.292 / 311.31 |
| lf3d-mpi-4 | 58.974 / 3739.98 | 58.764 / 3763.08 |
| smr-mpi-4 | 5.452 / 245.65 | 5.264 / 246.23 |

### T-P2: remove redundant LF flux clearing

Triage: LF STS register copies already select IEN/IAN and live blocks, and the
update already handles only those two variables. The remaining redundant clear
was removed from LF parabolic task registration in `MHD::AssembleMHDTasks` and
its now-unreachable LF branch removed from `MHD::ClearSTSFlux`. Every live LF
face assigns both updated flux slots before divergence or flux correction;
non-LF clears and A/magnetic-moment conversions are unchanged. The profiling
regression now verifies that the removed clear task has no profile row.

Release CPU/MPI builds, style checks and the profiling regression pass. All
24 configurations, repeated three times, are byte-identical to post-G, including
the SMR/outflow case that rejected P1. Expected values unchanged. Timings
(cycle ms / profiled compute microseconds per stage):

| Case | Post-G | P2 |
| --- | --- | --- |
| lf1d-serial-0 | 0.513 / 30.20 | 0.489 / 27.23 |
| lf2d-serial-0 | 17.949 / 1175.40 | 16.823 / 1093.94 |
| lf3d-serial-0 | 222.032 / 14492.02 | 211.354 / 13710.46 |
| smr-serial-0 | 17.653 / 927.78 | 16.294 / 825.40 |
| pure_cgl-serial-0 | 0.414 / n/a | 0.448 / n/a |
| lf2d-mpi-4 | 5.315 / 311.00 | 4.930 / 283.74 |
| lf3d-mpi-4 | 58.974 / 3739.98 | 55.703 / 3499.44 |
| smr-mpi-4 | 5.452 / 245.65 | 5.044 / 215.16 |

### T-P3: fuse uniform-grid primitive and temperature refresh

Triage: the one-ghost intermediate/full-ghost final refresh already existed.
`CGLLandauFluid::RefreshPrimitives` now recovers pressures and writes temperatures
in the same pass for uniform, non-shearing RKL2 LF. Stage 1 caches the frozen
magnetic magnitude; separate LF-scaled and C2P raw-square-root values preserve
the original rounding. A C2P overload accepts the cached norm while retaining
the original caller interface. The MHD sweep begin/end tasks reset cache state.
Density/velocity stores are retained when density repair is needed, and all
energy/moment repairs and floor counters remain. Multilevel, shear and explicit
reference paths retain their existing refresh. No numerical expression is
reassociated, and no runtime input option is added.

CPU/MPI Release builds, style and both precision C2P regressions pass. All
24 configurations times three repeats are byte-identical to post-G. A separate
uniform 16-block outflow case combines nonuniform density/B, nonzero velocities,
a scalar and ghost output: serial and MPI one/four ranks, three repeats each,
also match all eight output files exactly. Its precompute calls fall from 40
to 8 for 40 stages/four cycles, confirming cache reuse once per half-sweep.
Auxiliary evidence: `/tmp/cgl-wo1-p3-aux/`; primary matrix:
`/tmp/cgl-wo1-perf/runs/p3/`. Expected values unchanged. Timings before/after
this commit (cycle ms / profiled compute microseconds per stage):

| Case | P2 | P3 |
| --- | --- | --- |
| lf1d-serial-0 | 0.489 / 27.23 | 0.474 / 25.85 |
| lf2d-serial-0 | 16.823 / 1093.94 | 16.418 / 1061.74 |
| lf3d-serial-0 | 211.354 / 13710.46 | 204.298 / 13244.26 |
| smr-serial-0 | 16.294 / 825.40 | 16.166 / 819.10 |
| pure_cgl-serial-0 | 0.448 / n/a | 0.436 / n/a |
| lf2d-mpi-4 | 4.930 / 283.74 | 4.821 / 277.25 |
| lf3d-mpi-4 | 55.703 / 3499.44 | 55.044 / 3438.47 |
| smr-mpi-4 | 5.044 / 215.16 | 5.118 / 217.78 |

### T-P4: flat LF STS update already implemented

Triage: already fixed in the merged source. `MHD::STSUpdateU` selects the flat
`mhd_sts_update_cgl_lf_u` kernel for the LF-only STS path on uniform and multilevel
meshes, updating IEN/IAN over live blocks. Its directional divergence order is
unchanged. The remaining team kernel serves other parabolic operators and the
explicit reference path, so it is outside this LF STS task. No production edit
is needed.

The P3 matrix executes this kernel in 1D/2D/3D, SMR and shear and has zero byte
differences in all 24 configurations. After this documentation-only task, both
current executable hashes and all 24 retained output sets were rechecked against
the verified P3 binaries/post-G outputs. Before/after executable identity is
exact; the P3 timing table applies unchanged and no new speedup is attributed
to P4. Repeating identical timing runs would measure host noise only.

### T-P5: conditional STS register allocation already implemented

Triage: already fixed in `MHD::MHD`. Full-size `u_sts0/1/2/rhs` storage is
allocated only for `has_any_parabolic_cell_update`; full-size magnetic registers
only for `has_any_parabolic_field_update`. The cell predicate correctly includes
the explicit LF reference integrator, which also uses the parabolic registers.
Pure CGL retains only the existing one-element placeholders. No production edit
or allocation behavior change is needed.

The pure-CGL and LF baseline/P3 runs cover both allocation branches. As for P4,
after this documentation-only task the CPU/MPI executable hashes and all 24
output sets were rechecked: unchanged executables, byte-identical post-G output.
The before/after measurements are the same P3 table; P5 claims no new speedup
or new memory saving. Expected values unchanged.

### T-P6: evaluate only the selected HLLE anisotropy logarithm

Triage: HLLE still evaluated both side logarithms. `HLLE_CGL` now computes the
selected upwind value after the unchanged mass-flux sign test, preserving the
original `rho*log(...)/rho` order and weak-field reference selection. Its unused
left/right magnetic-moment flux assignments are removed. No reference values
change.

CPU/MPI builds and style pass. A direct old/new real-header comparison covers
32,768 states, three directions, active/passive signal paths, both precisions,
all flux/EMF and pressure-work outputs: zero differences. The full 24-case
serial/MPI matrix times three repeats matches post-G binary/history/full-restart
bytes exactly. P4/P5 changed no executable, so the incremental comparison is
P3 to P6 (cycle ms / profiled compute microseconds per stage):

| Case | P3 | P6 |
| --- | --- | --- |
| lf1d-serial-0 | 0.474 / 25.85 | 0.469 / 25.49 |
| lf2d-serial-0 | 16.418 / 1061.74 | 16.285 / 1052.44 |
| lf3d-serial-0 | 204.298 / 13244.26 | 202.263 / 13095.76 |
| smr-serial-0 | 16.166 / 819.10 | 16.011 / 809.41 |
| pure_cgl-serial-0 | 0.436 / n/a | 0.421 / n/a |
| lf2d-mpi-4 | 4.821 / 277.25 | 4.825 / 276.79 |
| lf3d-mpi-4 | 55.044 / 3438.47 | 54.860 / 3427.13 |
| smr-mpi-4 | 5.118 / 217.78 | 5.109 / 215.66 |

P6's incremental timing changes are small. Across all retained optimizations,
the required LF cases improve solver cycle time by 8.5-9.3% in serial and
6.3-9.2% on four MPI ranks relative to post-G. The tiny pure-CGL control varies
by +1.7%; no LF speedup is inferred from that control. Three repeats and this
one CPU do not establish GPU or production-scale performance.

## Validation and observations before review follow-up

- Release CPU and MPI builds pass. The final 187-check physical CPU selection
  had 186 passes and one G4 caller-guard failure; after repairing that guard,
  all 30 targeted limiter/collision rechecks pass, including the failed test.
  The selection includes both precision EOS/speed/closure/constructor checks,
  AMR, forcing, diffusion, boundary/FOFC and analysis regressions.
- The strengthened acceptance suite passes all 55 checks. Its 28 limiter cases
  also pass after the final guard repair. The full workflow passes all 29 cases;
  the separate G3 accuracy workflow passed all 33 cases.
- The final MPI suite passes 8 checks, including the corrected collisional
  shearing serial/MPI/explicit/restart test. Two GPU-allocation tests are skipped
  on this CPU host; no new GPU validation is claimed.
- All retained custom targets build: transform round trips, FOFC flux checks,
  and the seven-case legacy paper smoke with synthetic analyzer checks pass.
  Logs and diagnostics are under `/tmp/cgl-wo1-final-legacy/`.
- The 18-page report compiles and was visually checked. The final strict
  Sphinx build passes with only the existing orphan-document warning class
  suppressed. The two pre-existing orphan pages are named in G7's QA log.
- Five baseline campaign/provenance tests remain outside this work order, as
  listed at the start; 19 campaign/provenance cases were excluded from the
  focused final CPU selection. Full single-precision application compilation
  remains blocked by the pre-existing coordinate/hyperviscosity errors noted
  in B5; the changed math helpers pass their direct single-precision checks.
- P1 is explicitly deferred to WO2 because it changed active-cell results at
  an SMR/outflow interface. No part of that candidate is retained. P4/P5 were
  already implemented. Primitive prolongation and Spitzer conduction were
  already supported, so their tasks verify support/units instead of adding
  obsolete fences. These deviations follow the work order's triage rules.
- The retained legacy cell-centered heat-flux diagnostics are proxies, not
  independent face-flux references. Historical figures/campaign measurements
  were preserved and labeled; no renewed paper-scale validation is claimed.
- The RKL2 coefficients, stage-count/odd-stage rules, Strang half-sweep
  structure, multidimensional stability factor and mesh CFL multiplier are
  unchanged. No WO2 physics work was started.

Primary final logs are under `/tmp/cgl-wo1-final/`; performance outputs,
input/executable hashes, timings and rejected P1 evidence are under
`/tmp/cgl-wo1-perf/`. Work-order Markdown files and the Kokkos submodule pointer
are unchanged. No branch was pushed and no PR was opened.

Before the review follow-up below, the exact delivered CPU/MPI binaries
were run through the full fixed matrix again: all 24 configurations times three
repeats match post-G, with zero output mismatches or repeat nondeterminism.
`/tmp/cgl-wo1-perf/runs/final/provenance.json` records the final source diff and
these executable SHA-256 hashes:

- CPU: `71698e4a0c36b337e9b00998a02a9f8f7834c10f10131ecaeacc53442f198fb8`
- MPI: `7a92a266a18b740487c0d9470fd6bc4f4ddbc82ef219dd94b643b90964893e3e`

## Review follow-up validation

- Both original fixtures are restored: collisionless shear (`nu_coll=0`) and
  3D AMR pressure amplitude 0.5. Shear now runs 39 cycles to $t=0.3$; capped
  STS and explicit references run 994 cycles to $t=0.04$. Serial/MPI physical,
  magnetic, divergence and timestep assertions keep their original tolerances.
- The shearing restart already permits physical-state differences of $5\times10^{-6}$.
  Its discontinuous hard-bound occupancy count therefore need not be exactly
  equal near a wall. All other integer counters remain exact. Hard-bound visits
  remain positive and at most the stage-cell count; the occupied fraction differs
  by $2.55\times10^{-5}$, below its explicit absolute tolerance $10^{-4}$.
  Full-precision strict exit checks enforce the final wall independently.
- The broader CPU run gives 255 passes and one macOS test-fixture failure:
  `os.sched_getaffinity` is absent on macOS. Both affinity mocks now allow that
  absent attribute. The six-check targeted rerun passes, including this check,
  its sibling and the restored AMR/churn/restart cases. Thus all 256 selected CPU
  checks have passing results. All eight selected MPI checks pass.
- The five earlier campaign failures now pass with portable test-local trusted
  Git/Python paths and the correct campaign-input inventory. Production
  authentication was not changed. They were not solely Linux-path failures.
- Two new B4 grid checks remain strict expected failures, accepted as a known
  limitation on 2026-10-06. Their bounds were not relaxed. Both velocity directions fail;
  donor-cell reconstruction also fails. The isolated pre-B4 headers already give
  pressure ratio $1.11\times10^8$ after one cycle; current headers keep that first
  ratio at 1 but still reach $2.46\times10^{10}$ by cycle 50. No further B4 flux
  redesign is planned; this remains failed validation with an accepted limitation.
- A full float build remains blocked. Investigation additionally found table-reader
  `double*`/`Real*` mismatches and geodesic/unit-system narrowing and range issues.
  Trial portability edits were discarded; no unrelated float repair is retained.
  Changed CGL math helpers retain their passing direct float tests.
- B1 energy rounding runs only on a floor repair whose recovered pressure is
  still below the floor. C2 pressure rounding follows an actual wall clamp;
  encoded-A rounding requires a recovered hard-bound violation and never changes
  total energy. The collision kernel returns without re-encoding when rates and
  walls leave pressures unchanged.
- The final CPU/MPI binaries again match all 24 post-G output configurations,
  each repeated three times, byte for byte. The comparison harness accounts for
  the four added strict wall checkpoints per cycle without counting them as LF
  stages. These concurrent validation runs do not establish a revised speedup;
  the timing tables above describe the pre-review implementation.
- Changed source and tests pass style checks. The 18-page PDF compiles without
  overfull boxes; the two changed pages were rendered and checked. Strict Sphinx
  compilation passes with the same two pre-existing orphan-page warnings suppressed.
- `.github/workflows/cgl-lf.yml` adds GitHub-hosted Linux CPU/MPI checks, explicit
  restored AMR regressions, and a CUDA device-code compile using a pinned toolkit
  container. This Mac has no CUDA/HIP compiler, and the fork has no self-hosted
  runners. CUDA compilation subsequently passed at `7b3345fd2`; no GPU runtime
  validation is claimed.
- The first CUDA 12.6 build found NVCC's extended-lambda access restriction in
  private timestep helpers. Six affected declarations in hydro, MHD and the
  turbulence driver now have public access; their numerical bodies are unchanged.
  The local CPU rebuild is byte-identical to the previously tested binary.
  The initial compiler diagnostic is retained in `validation/wo1/review/cuda/`.
- The next CUDA build exposed a first-capture restriction in the three
  directional CGL flux kernels. `7b3345fd2` captures the existing pressure-work
  Boolean before `if constexpr`, following the surrounding code pattern.
  CUDA compilation and Linux CPU/MPI validation then passed. The local build,
  style, 12 fixed serial cases with three repeats, and six dynamic HLLE cases
  with pressure-work recording on/off preserve the previous results.

Durable provenance, output hashes, rejected P1 patch, inputs and review evidence
are committed under [`validation/wo1/`](validation/wo1/ARCHIVE.md). The P1 issue is
recorded as a WO2 correctness investigation. The two root work-order files remain
untracked and are not included. The accepted B4 limitation no longer requires
draft status by itself; other outstanding validation requirements still apply.

## WO1 closeout, 2026-10-06

WO1 is complete with the documented B4 acceptance decision and the required
P1 deferral. All 41 numbered tasks have individual commits and per-task reports
above. Retain conservative A and the current numerical implementation. The
extreme sharp-contact test remains a strict expected failure with unchanged
bounds; it no longer blocks WO1. The reference comparison is archived beside
the [acceptance note](validation/wo1/review/weak-field/README.md).

The final production code at `7b3345fd262a1cf826e9476fb40ca20a28a15881` passed
[Linux CI run 36881505251](https://github.com/dfielding14/athenak/actions/runs/36881505251):
213 CPU checks passed with the two B4 expected failures, all three additional
AMR checks passed, all seven MPI checks passed, and the complete CUDA 12.6
application compiled. The closeout changes no production source or input.
The only test edit updates B4's expected-failure explanation to the accepted
decision; its assertions and strict marker are unchanged.

The local Release build and changed Python style check pass. Both B4 directions
were rerun: normal pytest reports two expected failures; `--runxfail` exposes
the actual cycle-5 ratio failure, 3.0052200316546487 versus the upper bound 2.
Thus the expected-failure marker is recording the known numerical limitation,
not masking a setup failure. Logs, CI test XML, CUDA build output and provenance
are in [`validation/wo1/review/closeout/`](validation/wo1/review/closeout/).

The validation scope is double-precision CPU/MPI execution and CUDA compilation.
GPU runtime validation was not performed. Full single-precision application
compilation still encounters pre-existing portability errors outside WO1; the
changed CGL math and constructor paths have their direct float checks. These
limitations remain disclosed and do not imply either unsupported target passed.
No new turbulence-validation claim follows from the B4 comparison. P1's rejected
non-bitwise change and broader numerical redesigns remain WO2 work.

PR #21 carries the final branch and live checks. The two root work orders remain
local and untracked; the submodule revision and scientific validation thresholds
are unchanged.
