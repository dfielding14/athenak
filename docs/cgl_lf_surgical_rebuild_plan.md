# CGL LF surgical rebuild plan

Reconstruct the required physics from pre-WO1 source, measuring each correction
and optimization separately. Keep the existing AthenaK structure. A previously
accepted physical requirement does not justify copying its entire implementation.

Branch: `c/cgl-lf-surgical-rebuild`.
Base: `8222de3aae4aaebd886653e7c61846e5f17987b4`.
Status: source-reviewed implementation plan only. No solver/test implementation,
builds, simulations, or numerical tests have been performed on this branch.
Implementation and execution begin on the GPU system, not on this planning machine.

The first deliverable is a correct active solver with attributable cost. The
passive model is a separate scientific deliverable. Optimizations are candidates,
not a list that must all be implemented. Keep an existing implementation when
its independent check passes; add a test only for an uncovered failure or changed
contract. Reuse retained tests without importing the surrounding campaign tools.
The scientific acceptance target is time to a specified physical endpoint at
demonstrated accuracy. Zone-cycles/s and cost per LF evaluation explain that
result; neither alone establishes a faster usable solver.

Read the numbered pieces for source scope, the test map for concrete selectors,
and the execution cadence for when to run them. In particular, the gates below
are not instructions to repeat the full CPU/GPU/MPI matrix after every edit.

## Reference points

| Source | Purpose |
| --- | --- |
| `8222de3aae4aaebd886653e7c61846e5f17987b4` | Exact pre-WO1 implementation baseline. Contains known defects; not a qualified production model. |
| `f2d0a25978e4fad9beae6874f279148c0460533b` | Older historical scaling reference. Do not substitute its throughput for a measurement of the rebuild baseline. |
| `9a4030b09307bc43865d5e597638f8645a6388f8` | Completed WO1, now the tip of `CGL-STS-LF`. Compare against it to isolate WO1. |
| `7a37710f6c224e24e7c7f364e7e0b812b3a9494c` | WO2 release reference. |
| `2afaa0fa1e4985c313c132cd97ba3938e8c2040c` | Preserved newer implementation and evidence, including subsequent numerical fixes. |

Read old files with `git show <revision>:<path>`; do not merge the old branch.
The newer reference contains `WO2_CGL_LF_followups.md`,
`docs/cgl_lf_changes.md`, `docs/validation/wo1/`, `docs/validation/wo2/`,
`docs/source/modules/cgl_lf_timestep.md`, `docs/source/modules/cgl_passive.md`,
and `docs/cgl_lf_performance_audit_20261008.md`.
The original WO1 is local at
`/Users/dbf75/Work/Research/AthenaK/athenak-DF/WO1_CGL_LF_fixes.md`.
Use these as physical requirements and evidence, rechecking applicability to
the current patch. Their historical execution order is superseded by this plan.

WO1's performance harness began after its correctness changes. It does not
clear those changes of performance regressions. Safe LF arithmetic, weighted
fluxes, full diagnostics by default, and per-stage admissibility checks already
exist at the base; their cost cannot simply be attributed to WO1.

## Physics and implementation constraints

- Preserve active conservative A and total energy, with thermal energy
  $U=p_\perp+p_\parallel/2$. LF temporarily evolves magnetic moment and restores A.
- Preserve RKL2 coefficients, odd-stage selection, and LF half-step / RK full-step /
  LF half-step ordering. The specified once-per-cycle collision update is not a
  symmetric collision split; do not claim second-order accuracy for that full model.
- Restore the agreed collision and limiter model: rates once with the full outer
  timestep; hard walls at the specified operator boundaries and after regridding.
  Preserve the mandatory fluid-firehose bound $p_\perp-p_\parallel\ge-B^2$ at
  those wall checkpoints, including configurations without optional limiters.
  Do not add per-RKL-stage pressure clipping. Intermediate wall crossings and
  nonpositive/nonfinite states have different meanings and must remain distinguishable.
- Preserve forcing power, its random-process cadence, and exact kinetic-work
  accounting. Keep checks that prevent invalid states; detailed q diagnostics can
  be off for performance measurements.
- Already present: the corrected range-scaled fast speed, perpendicular-floor
  energy factor, A/magnetic-moment interfaces, CGL primitive prolongation, flat
  two-variable STS updates, and conditional registers. Test these; do not rewrite them.
- Each numbered piece below is a reviewable unit. Separate its independent
  sub-fixes into commits; keep coupled operator/caller changes together when an
  intermediate version would implement the wrong equations. Include the focused
  regression in the same commit. Do not cherry-pick whole historical task files.
- Preserve exceptional arithmetic that handles demonstrated finite states. The
  existing scaled LF arithmetic header and RKL2 coefficient implementation are
  already present at the base and unchanged in the newer reference. Reuse them.
  Shorten duplicated storage and unnecessary passes, not the physical operator.

The collision contract is explicit. With $\Delta=p_\perp-p_\parallel$, background
relaxation gives $\Delta_b=\Delta\exp(-\nu_{\rm coll}\,dt)$. Beyond an enabled soft
threshold $\Delta_*$, use the stable backward-Euler form
$\Delta'=\Delta_*+(\Delta_b-\Delta_*)/(1+\nu_{\rm lim}\,dt)$ at fixed
$(p_\parallel+2p_\perp)/3$, then apply the prescribed walls.
Default soft thresholds are $-B^2$ and $B^2/2$; backup factors are respectively
1 and 2, and backup scattering defaults to $10^{10}$ in code time units.
LF scattering adds background, applicable soft, and enabled backup contributions
once each; equality contributes no extra scattering. Keep numeric threshold
parameters and the intended legacy-policy mapping, not a second limiter model.

Dependency order: establish 0, then local repairs 1-3; 4 owns the coherent
collision model; 5-8 correct face transport; 9-10 qualify its controller; 11
corrects forcing; 12 qualifies refinement. Piece 11 can proceed after 4 without
waiting for LF tuning. Piece 13 follows the relevant corrected reference; 14
does not depend on optional optimizations or multilevel work in 12. The uniform
active milestone is a usable stopping point; add refinement or passive support
when that model is needed. Nonlinear temporal qualification applies to both
models; the chunk implementation in 15 remains conditional on demonstrated need.
Keep separate commits at the performance suspects, even when their functional
checks share a single GPU allocation.

Keep the representation changes visible at existing operator boundaries:

| Model / operator | IEN | IAN | Quantities held fixed during LF |
| --- | --- | --- | --- |
| Active, ordinary evolution | Total energy E | A | — |
| Active, LF sweep | Total energy E | Magnetic moment density | Density, momentum, B; E changes through heat flux |
| Passive, ordinary evolution | J | A | — |
| Passive, LF sweep | Thermal energy U | Magnetic moment density | Density, momentum, B; U changes through heat flux |

Reuse the current EOS conversions, face helpers, task lists and STS controller.
Keep exceptional arithmetic beside the formula it protects. Do not add a generic
state wrapper, representation dispatcher, parallel scheduler or new cache manager
to express this table. Comments should explain the invariant or dependency,
not narrate the loop. Keep physical-oracle tests independent of shared code.

## Implementation sequence

All implementation pieces are pending. Run each piece's focused gate before
stacking dependent changes; batch broader integration at the milestones below.
Record existing failures rather than changing tolerances or weakening fixtures.

**0. Establish the baseline before changing production source.**
Build the exact base and post-WO1 reference with identical qualified compiler,
Kokkos, precision and runtime settings. Retain the documented CCE20/HIP corrective
flags on both builds; verify one minimum-reduction reproducer, not the historical
compiler-search matrix. Record the small baseline checks and paired timings below.
Confirm the existing fast speed against non-unit-density anisotropic limits and
independent eigenvalues using the existing test, strengthened only where needed.
Use current source to locate the real test entry points; old work-order line
numbers and several original findings are stale.
Keep both binaries as controls and build the candidate incrementally in its own
directory. Do not rebuild the immutable controls after each source patch. The
comparison identifies which WO1 costs deserve isolated measurement; it does not
require speculative optimizations before correcting the equations.

**1. Keep unsupported models out of the accepted domain.**
Add the small passive-mode rejection until its redesign is complete; preserve
the existing LF inflow/user-boundary guard. Validate existing CGL floor and rate
parameters at their constructors. New threshold/backup parameters arrive with
their consuming physics in 4, not as accepted-but-ignored options here.
Keep the Spitzer units fix outside this CGL series unless that conduction path
is needed. Scope: EOS initialization and affected input/startup tests.
Gate: invalid and unsupported inputs fail clearly; valid active cases are unchanged.

**2. Repair the remaining EOS floor defects individually.**
Scope: `SingleC2P_CGLMHD` and `SingleC2P_CGLMHDFromMagneticMoment` in
`src/eos/ideal_c2p_mhd.hpp`; repaired-state persistence in
`CGLMHD::ConsToPrim`, `CGLMagneticMomentToPrim`, and the lightweight refresh.
For otherwise valid magnetized states, preserve the pressure ratio under density
floors in the A decoder. In the magnetic-moment decoder, retain the accepted
different density repair: at fixed B and without further pressure-floor repair,
preserve magnetic moment (thus perpendicular pressure), momentum and total
energy; the changed kinetic energy changes parallel pressure. Do not impose
ratio preservation on that representation. Enforce the magnetization ceiling
through the density floor in both decoders, and make repaired states representable
and idempotent. Do not redo the already-correct perpendicular energy factor.
Gate: analytic pressure ratios and energy identity,
repeated C2P with no second change, changed A persisted, focused double/float
checks. A focused float test does not establish whole-application float support.
Extend the existing floor fixture selectively: the later complete fixture also
calls `CheckCollisionMap`, which depends on piece 4. Do not import that dependency
or passive decoder/link machinery into this local repair.

**3. Repair reconstruction and FOFC detection.**
Scope: all directional PPMX and WENO-Z wrappers and the C2P test-floor path;
leave PPM4 unchanged.
Floor both reconstructed CGL pressures; detect nonfinite incoming state before
repair obscures it. Gate: actual reconstruction wrappers, injected NaN/Inf in
E/A, a steep-pressure evolution, and ordinary MHD controls. Change matching
reference consumers only where their old expectation is demonstrably wrong.

**4. Restore thresholds, relaxation, walls and their schedule coherently.**
Scope: `EOS_Data`, `cgl_physics.hpp`, `SingleCollRates_CGLMHD`,
`SingleCollWalls_CGLMHD`, `CGLWallAdmissibleAnisotropy`, `CGLMHD::Collisions`,
and the mode/signature declaration plus the base stub in `src/eos/eos.cpp`.
Update `MHD::AssembleMHDTasks`, `CGLCollisions`, `STSPostSweepCGLCollisions`,
AMR projection, LF face-rate consumers and affected diagnostics/inputs together.
Introduce the agreed soft/backup thresholds and additive scattering, then land
the monotone relaxation map with its full-step schedule. Keep coarse state and
primitives consistent after updates. Incorporate the later canonical active wall
encoding correction here (`71ad25ebc`, regression `4cd46a131`), rather than
reintroducing its known restart failure.
The call sites choose the mode: no LF split means full rates/walls after RK;
split LF means walls after the pre-sweep and RK, then full rates/walls after
the post-sweep, with A restored before those calls. After the final AMR field
refresh, apply walls and refresh the coarse representation. Do not retain an
obsolete soft-threshold projection in C2P or AMR alongside the new finite-rate map.
Honor the specified default thresholds (firehose 2, mirror 1), explicit legacy
equivalence and finite-rate fixtures. Update shipped runnable inputs that use
the removed `limiter_hardwall=true` contract with the operator change; preserve
historical evidence. No input should silently acquire a different physical rate.
Gate: analytic background relaxation with/without LF; an independent finite-rate
map distinguishing one full kick from two half kicks; threshold equality and
additive-rate checks; invariant E/momentum/B; unchanged A on no-op; wall
idempotence; a separate production C2P after collision encoding; real restart.
Give the post-RK timestep one owner. The newer reference registers the ordinary
last-stage `MHD::NewTimeStep` and then recomputes it in `CGLCollisions`, before
the driver consumes the budget. Once the unconditional post-wall CGL hook is
installed, omit only the redundant CGL task registration; keep `NewTimeStep`
callable for initialization and regridding. Keep the post-LF refresh and global
post-RK minimum.
This small dependent deletion can follow 4 immediately; it need not wait for 13.
Use existing heated/wall/restart checks and call counts to verify the ownership,
including zero rates and no LF. Ordinary MHD keeps its last-stage task.

**5. Restore the small weak-field face-transport correction.**
Scope: HLLE/LLF face helpers and their direct/evolution fixtures (`7d1f79557`).
Retain conservative A; use the appropriate neighboring field at the weak-field
face. Gate: first-update contamination, both contact orientations, and smooth
fixed-width interface convergence. Preserve the two documented B4 strict
sharp-contact expected failures with their original predicates. This patch does
not solve that formulation's sharp-contact limitation.

**6. Correct LF face geometry.**
Scope: `BuildCGLLFFaceState`, `AddHeatFluxes`, the header and
`MHD::AddSelectedDiffusionFluxes`, plus reference consumers (`910023601`).
Use CT normal B and the mean of adjacent cell magnitudes
for normalization; do not renormalize this vector to unit length.
Gate: independent face fluxes, reversal-sheet perturbations, frozen-field
conservation, and smooth oblique decay. Use a small explicit reference timestep
while the final stability controller is still pending.

**7. Limit transverse temperature gradients.**
Scope: existing directional face arithmetic and the smallest reusable VL4 helper
(`83f407484`). Leave magnetic-magnitude slopes unlimited.
Gate: hotspot minimum and energy, oblique convergence, and exceptional-slope
checks. Measure this change alone: it replaces centered sums with several
harmonic means at every face and is a leading WO1 cost suspect.
Use the retained `VanLeerLimiter`/`VL4Limiter` arithmetic. Update existing
diagnostic/profile copies of the directional expressions when they consume the
same gradients; do not introduce a new general stencil abstraction.

**8. Correct the perpendicular BGK response.**
Scope: `cgl::PerpendicularHeatFluxRatio` in `src/eos/cgl_physics.hpp`, its normal
and logarithmic paths, `ChiPerp` reference consumers and associated method docs
(`d2ab32f58`). Both safe and fast fluxes use that shared helper. Additive scattering
was completed in 4; do not reimplement it here. Use the intended $2\nu_{\rm eff}$
perpendicular response. Gate: collisionless identity, strong-collision limit
$\chi_\perp\nu/c_\parallel^2\to1$, and resolved coupled decay over an e-folding.
Keep pieces 6-8 at a common conservative reference timestep until piece 9
qualifies the final controller. For their timing screen, enlarge the existing
kinematic `cgl_lf_rotated_decay.athinput` with an oblique field and background
coefficient: zero velocity and fixed density/B keep the pre-9 budget independent
of the evolving pressure. Verify actual dt and RHS counts; selecting explicit
integration alone does not fix dt, and there is no general fixed-dt input at the
base. This is a short patch screen, not another benchmark campaign. Do not use
unrestricted turbulence before the forcing correction to attribute kernel cost.

**9. Implement the LF bound for the final stencil.**
Scope: timestep calculation/callers, existing post-RK refresh, and MPI reduction.
Select the mathematical blocks from `822cad5e8`: `CGLLFVL4DerivativeNorm`, the
local logarithmic helpers, representable-sound-speed fallback, and the row loop
in `CGLLandauFluid::NewTimeStep`. Add face B at the `mhd_newdt.cpp` caller.
Retain CT geometry, density coupling, limiter derivatives and grad-B reverse
coupling. The initial bound reuses existing `tpar_`, `tperp_`, `bmag_` storage
with an explicit fresh one-halo fill: no three new timestep arrays and no new
cache-validity framework. Fused primitive refresh does not exist until 13;
reusing its valid values is a later optimization, not a dependency here.
Include background collisions. Do not first port WO1's provisional D5 bound merely to
replace it, and do not loosen unrelated conduction bounds.
Gate: independent frozen-stencil rows/eigenvalues, staggered checkerboard,
unequal spacing, density-jump seeded-minus-unseeded growth, reversal/grad-B,
heated post-RK refresh, MPI stage agreement, and explicit-reference convergence.
Initially retain the old CFL safety multiplier. The bound covers a frozen local
Jacobian; it does not prove nonlinear positivity, RKL stability, or AMR stability.
Do not substitute a simpler maximum-diffusivity or secant-limiter estimate.
The original coefficient-five proposal underbounded the actual VL4 stencil.
Retained JSON matrices establish old results, not a rerunnable independent test;
recover the actual checker or add one compact independent row/Jacobian check if
the reusable runtime fixtures do not cover the adapted mathematics.
In that same checker, retain real/imaginary eigenvalues and RKL amplification
at the selected step, not just spectral radius versus row norm. Interpret growth
against the explicit reference; a positive or nonnormal frozen operator is not
automatically a timestep bug. Do not turn this into another parameter campaign.

**10. Separate STS safety from advective CFL.**
Only after piece 9 passes, add the parameter in `mesh.hpp`, its fresh/restart
loading in `build_tree.cpp`, and both budget selections in `mesh.cpp`;
explicit diffusion keeps its CFL multiplier. Gate: parameter bounds, analytic
stage selection, existing stability tests, and temporal convergence before
accepting a larger factor such as 0.9. Attribute changes in stage count separately
from kernel speed. Do not alter RKL coefficients or its odd-stage rule.

**11. Restore forcing corrections one at a time.**
Scope: shared turbulence driver, task placement, primitive/ghost refresh and
affected problem generators/fixtures (`475e60ab4`, `951e70277`, `fd9c8a1a1`,
`2769a3c40`). Restore type-2 modes, one full-dt kick, energy from each fluid's own
conserved momentum, and the specified planar axis. Retain the minimum refresh
required by the next flux evaluation; measure its cost before narrowing it.
Gate: fixed-seed force modes, RK1/2/3 kick counts, injected power, unchanged
thermal energy from a pure forcing kick, OU cadence/restart, and affected
ordinary hydro/MHD and two-fluid consumers. Do not silently retune power or seeds.
Keep mode selection/counting changes in `IsDrivenMode`/`InitializeModes` separate
from `ApplyForcingWithStep` energy arithmetic. Move `IncludeInitializeModesTask`
and `AddForcing` with their `MeshBlockPack` registration and `Driver::Execute`
refresh atomically. Replace the old RK-recurrence test expectation with the
physical full-kick expectation. Register the small retained `turb_forcing` pgen
only if needed for the direct work/cadence fixture; do not add a forcing framework.

**12. Repair refinement/boundary correctness before compact communication.**
Two separate fixes: refill physical corners after prolongation (`e1b2ecd00`),
then synchronize shared LF face fluxes near refinement (`a6eef6406`). Reuse the
existing buffers and A/magnetic-moment-aware prolongation; retain magnetic BCs.
Gate: uniform anisotropic oblique-field SMR/outflow state, the original nonuniform
fixture, periodic E and integrated magnetic-moment conservation, both prolongation
modes, 2D/3D, explicit/STS, CPU/HIP and multiple ranks. Preserve AMR churn/restart
and repair-accounting checks. Do not equate positive pressures with conservation.
The original shared-face patch followed compact communication, but need not
depend on it here: add default-off same-face offset/count arguments to flux
pack/unpack and the receive overload; leave coarse/fine messages full-variable
and cell-state exchange unchanged. Match receive/send/unpack activation and
retain existing request completion. Before posting receives, grow both face
buffers if needed to `max(existing capacity, 2 * full face area)`, preserving
their block capacity. Six-variable coarse-face storage provides only 1.5 face
areas in 3D and is insufficient. Cache both original estimates before the
identically ordered symmetric mean. No new completion framework is needed.

**13. Apply only measured, equivalent optimizations.**
Try separately: redundant LF flux-clear removal (`15e7ff619`); fused primitive/
temperature refresh with valid frozen-B reuse (`d02bb714d`); upwind-only HLLE
logarithm (`31b739890`); compact IEN/IAN messages (`6539bfbbc`); then any proven
redundant timestep refresh. Keep full shearing exchange and required magnetic
boundary work. Gate: same-configuration full-precision state/restart equality,
including ghosts/coarse data and scalars, plus paired GPU timing. Keep only a
beneficial or independently justified simplification; no configurable fallback
framework. Handle diagnostic reduction-order differences explicitly.
Only change a path that still contains the waste. For refresh fusion, preserve
the different raw/scaled B norms needed for identical rounding, repair stores,
floor counters, and the final full halo. Cache density/B only while frozen;
invalidate at existing sweep boundaries. For compact exchange, receive sizes,
packing, unpacking and call-site offsets/counts must land together. A simple
deletion can be worthwhile even when timing is within noise; a large new kernel
or cache mechanism needs a repeatable benefit on the target with q diagnostics off.
Start with the existing non-detailed profile of the corrected target. Optimize
the dominant recurring work; do not implement every historical optimization in
the listed order. In particular, reducing bound cost at the expense of extra LF
evaluations can make the whole calculation slower.

**14. Add passive thermodynamics as a separate model series.**
Reuse the physical design of `a2d994a3a` and the later decoding correction
`ab9b543e7`. Scope must include EOS, fluxes, forcing, initialization, histories,
output labels, restart version/representation and support guards. The IEN slot
has a different meaning; active total-energy restarts are not passive J/A data.
Start with the implemented uniform-periodic supported domain; no speculative
AMR/boundary/source extensions. Gate: density/momentum/B identical to native
isothermal flow at matched outer timesteps and forcing; thermal evolution from
independent smooth references; pressure-work balance; wall decode and restart.
Document the absence of irreversible shock heating. Optimize only against the
preceding implementation of this same model, never against the old passive defect.
Keep the native isothermal flow kernel intact: the attempted combined face
implementation previously changed HIP flow bits. Implement thermal additions
behind the passive rejection, then enable the model only when conversions,
fluxes, initialization, outputs and restart contracts are complete. For reference,
$J=\rho\ln(p_\parallel B_{\rm eff}^2/\rho^3)$ is stored in IEN outside LF;
LF uses physical U and magnetic moment instead. A direct source increment to
IEN is therefore not a passive heating operation.
After model qualification, profile before considering narrower thermal
reconstruction or a cheaper wall path. An eager isotropic fallback and bisection
are candidates for simplification, not reasons to alter the thermal model.

**15. Add LF chunking only if nonlinear accuracy or positivity requires it.**
The active/passive milestones below first qualify temporal accuracy without
requiring a new integrator feature. If either target needs smaller LF intervals
at a fixed outer step, use `4af7c9a55` as a small implementation reference;
the retained passive depleted-pressure failure is the known motivating case.
Hold the outer step, forcing and collision cadence fixed; change only LF substeps
within each sweep.
Gate: actual failing-state reproduction and timestep refinement, worst-cell
pressure differences as well as norms, energy and smooth temporal order. Recorded
positive chunked runs still differed locally by about 14%; cap 4 is not an
accuracy-qualified default. Do not merge chunking with a kernel optimization.
Retain the existing smooth refinement, cadence, default-off and restart checks.
Their current test module imports half-sweep-merge fixture helpers; extract only
the needed input/run/comparison pieces, not the merge feature or its campaign.
For the failing state, compare at least three decreasing LF intervals at fixed
outer step on the same grid. Report worst-cell absolute and relative pressure
errors, local thermal-energy-scaled errors, and norms; establish decreasing
temporal error rather than selecting a cap because it survives. The old cap-1
trajectory is not a converged oracle. If refinement remains unresolved, report
that numerical limitation instead of silently lowering forcing or changing walls.

## Deferred work

- LF inflow/user boundaries: implement only when needed, using `738810c88` and
  independent uniform-state/callback tests across all newly supported paths.
- Directional fusion: consider after profiling the corrected workload with q
  diagnostics off. Retain only a measured gain with one readable face-arithmetic
  implementation; do not copy the large duplicate fused kernel by default.
- Half-sweep merging, entropy-variable fluxes, A-to-Q reformulation, new limiter
  models, symmetric collision splitting, and broad float portability are outside
  this reconstruction. No speculative runtime switches or analysis framework.
- Leave unrelated cleanup and the large evidence archive on the reference branch.
  Bring across only the fixture, checker and short derivation needed for a change.

One conditional controller simplification deserves a note, not an implementation
commit yet. The VL4 transverse derivative envelope for finite slopes is bounded
by $\Gamma=8$. If the bound remains a measured bottleneck, compare a scratch
implementation using 8 in the same coupled CT/density/grad-B rows against the
slope-dependent bound. It can remove the timestep's two temperature scratch
reads/fill and derivative arithmetic, but must preserve fresh B and fail-closed
exceptional-state handling. For constant density/B and frozen uniform
$\chi_\parallel$, isotropic pressure, a 45-degree field on a square 2D grid and
equal nonzero slopes, the parallel row increases from $6\chi_\parallel/h^2$ to
$20\chi_\parallel/h^2$, potentially about 1.83 times as many stages asymptotically;
aligned fields have no transverse penalty. This changes the controller, so use
9-10's accuracy/work checks, not 13's bitwise-equivalence claim. Judge actual odd
stage counts and physical-time cost at matched accuracy. Retain only a clear
overall benefit and one implementation, with no new runtime mode. This is a
valid conservative alternative, unlike the rejected secant estimate, but is
not the initial controller choice. An early zero-normal-B shortcut is also
deferred: it must not bypass the current rejection of invalid face sound speeds.

## Testing cadence

Reuse `tst/test_suite/`, `tst/scripts/cgl/`, `test_suite.testutils`, existing
inputs and production pgen fixtures. A gate is evidence for the changed behavior,
not a required new test file. Do not import whole multi-thousand-line test modules
to obtain one regression. Add their needed functions/fixtures to the existing
structure with the owning patch. No whole-repository cleanup or new test framework.

| Change | Immediate check | Broader check or timing trigger |
| --- | --- | --- |
| Parameters, guards, documentation | Existing valid/invalid constructor cases or text review | No throughput timing; no new MPI run |
| Local EOS/reconstruction repair | Relevant arithmetic/wrapper fixture; one affected evolution if not covered directly | Actual HIP reproducer for known backend rounding; broader integration at the active milestone |
| Collision/forcing schedule | Analytic cadence/work check and affected restart | HIP application; one/four ranks when scheduling or reduction changes; a paired timing if full-domain work changes |
| Geometry/VL4/BGK | Focused arithmetic and one defect-sensitive evolution | Common conservative reference timestep until 9; one paired timing for added recurring arithmetic, especially VL4 |
| Timestep controller | Small analytic family, fresh-state and MPI stage-count checks | Compare both work count and physical-time cost; nonlinear case at the active milestone |
| Refinement/communication | One affected one/four-rank comparison | Remaining dimensions/prolongation/restart/shear at the mesh milestone; timing for message optimization |
| Equivalent optimization | Exact primary state, relevant ghosts/coarse data and restart comparison | Paired GPU timing on the affected workload; expand coverage only for paths actually changed |
| Passive model | Identity and independent thermal groups while assembling the model | One model milestone across supported reconstructors, forced restart and MPI; failing-state LF refinement separately |

**Baseline milestone:** on each frozen reference, one focused GPU decay,
one 3D runtime-path case, final-refresh/layout agreement, and actual restart.
Run the existing host fast-speed/floor/closure checks once. These establish
execution and the starting failure list; they are not a new broad qualification
of the known-defective base.

**Active milestone:** after 1-11, run the selected physical suite once: resolved
non-integer-period waves, an e-folding of parallel/perpendicular decay, oblique
transport, finite-rate relaxation, heating/timestep refresh, forcing/restart,
and a resolved active turbulence case with nonzero LF work. Run an actual HIP
wall-encoding fixture; CPU scalar helpers cannot substitute for it. Run the
small safe/fast and diagnostic/profile branch checks once after all face changes,
not after each constant edit. Production timings keep detailed q diagnostics off.
Add one short nonlinear temporal-refinement comparison using the existing
64-by-64 hotspot with `lf_coefficient_mode=local`, kinematic/advect evolution,
zero velocity and no forcing/collisions. Compare full-precision parallel and
perpendicular pressures on the same grid and physical endpoint across at least
three decreasing `sts_max_dt_ratio` values; 1, 0.5, 0.25, with 0.125 as a finer
reference if needed, is a starting ladder, not prevalidated accuracy. Confirm the
actual steps and RHS counts, preserve the existing minima/conservation checks
and its at-least-five nominal diffusion-time duration (`tlim=0.201` in the input).
Report worst-cell absolute/relative and thermal-energy-scaled errors as well as
norms, and require decreasing temporal error. This adds an assertion to the
existing fixture, not a new pgen or chunk feature. It checks nonlinear LF and
representation conversions; it does not certify turbulent threshold switching.

At this milestone, compare once against a qualified executable from `2afaa0fa1`
with the same active physics/settings, merging and chunking off. Reuse a verified
binary or make one additional build, not another full validation matrix. The
pre-WO1 controls locate historical costs; this corrected reference answers
whether the smaller reconstruction improves on the newer solver. Use physical
endpoint checks for intentional arithmetic differences. Exact equivalence of
each optimization remains a comparison with its immediate rebuilt parent.

**Mesh milestone:** after 12, cover both physical-corner cases, all six retained
conservation parameters, actual regrid/restart, repair accounting and shearing.
One/four ranks are the default comparison. Use two-node/16-rank 3D churn only at
final qualification when the allocation supports it, not as a prerequisite to
every local fix. These checks gate multilevel/communication acceptance; they do
not prevent independent uniform-grid work while an unrelated mesh issue is diagnosed.

**Passive milestone:** start with identity and independent thermal evolution;
then run the six existing groups once on HIP, including all supported
reconstructors. Check forced identity/restart on multiple ranks. Do not repeat
the historic hundreds-of-application CPU/HIP/MPI cross-product after each edit.
Require a short nonlinear temporal comparison in the supported passive model;
reuse the smooth refinement and retained depleted-pressure state as applicable.
Keep forcing and rate cadence matched for an LF-only comparison. Reducing the
outer timestep in a forced run changes that composition and is not an isolated
LF refinement. Piece 15 is implemented only when needed to refine LF intervals
independently; a small smooth test or clean termination does not establish
turbulent accuracy. Do not relabel the active kinematic fixture as a passive
test: the current passive guard requires dynamic evolution.

For these scientific milestones, read the existing LF floor/invalid counters,
endpoint pressure extrema, and EOS/FOFC event log. Q diagnostics can stay off;
strict LF admissibility still counts and rejects invalid/floor-repaired stages.
The event log covers repair activity outside LF. Investigate changed or dominant
repair activity before treating the result as accurate physics, without imposing
a universal zero-repair criterion on every turbulent case. Counts may include
ghosts and repeated refreshes; they are not distinct-cell fractions or an energy
budget. Read these in functional runs; add no hot-path instrumentation or new
logging system.

Use inherited scientific tolerances unless the changed equation requires a new
analytic expectation. Final-reference tests sometimes hard-code safety 0.9,
13 heated stages, zero stages for a zero operator, or a timestep cap to hit an AMR
event. Derive the expectation for the controller currently being tested. Preserve
the physical error, conservation, amplitude, topology and strictness assertions;
do not copy a cap blindly or weaken a criterion just to recover a pass.

## Concrete test map

All selectors below were verified in Git source, not executed here. Paths are
relative to `tst/test_suite/`. `L` means `cgl/test_cgl_landau_fluid_cpu.py`,
`LM` means `cgl/test_cgl_landau_fluid_mpicpu.py`, and `A` means
`cgl/test_cgl_lf_acceptance_cpu.py`. Expand a selection as `path::function`.
The newer test source is pinned at `2afaa0fa1`; import only the checks absent
from the base and needed for the current piece.

| Piece | Reusable selectors | Fixture or role |
| --- | --- | --- |
| 0 | `L::test_cgl_lf_quantitative_decay_and_diagnostics`; `L::test_cgl_lf_3d_runtime_modes_exercise_directional_fast_paths`; `L::test_cgl_lf_final_refresh_is_meshblock_layout_independent`; `L::test_cgl_lf_restart_preserves_final_state_and_admissibility` | Already present at the base; select a small representative parameter case first |
| 0 | `cgl/test_cgl_fast_speed_cpu.py::test_cgl_active_hlle_uses_scaled_literature_discriminant` | Existing focused Serial test; strengthen independent anisotropic eigenvalues from WO1 if absent |
| 1, 4 | `cgl/test_cgl_parameters_cpu.py::test_cgl_constructor_parameters` | Import parameter cases only when their consuming implementation lands |
| 2 | `cgl/test_cgl_c2p_pressure_floor_cpu.py::test_cgl_c2p_pressure_floor_energy_consistency` | Existing checker, extended progressively; double and focused float |
| 3, 5 | `L::test_cgl_reconstruction_pressure_peak`; `L::test_cgl_fofc_live_flux_mutation`; `cgl/test_cgl_weak_field_cpu.py::test_weak_field_contact_remains_bounded` | `inputs/unit_tests/cgl_reconstruction_{ppmx,wenoz}.athinput`, built-in FOFC pgen, preserved B4 strict expected failures |
| 4 | `L::test_cgl_collision_rates_once_per_cycle_with_and_without_lf`; `L::test_cgl_lf_face_collision_rates_add`; `A::test_limiter_stress_matches_cellwise_relaxation`; `A::test_limiter_stress_matches_wall_ordering` | `cgl_collision_once.athinput` and existing mirror/firehose unit inputs; independent rate/cadence oracles |
| 4 | `L::test_cgl_collision_refreshes_next_timestep_from_relaxed_state`; `L::test_cgl_lf_restart_with_finite_collision_preserves_corrected_split`; `cgl/test_cgl_amr_walls_cpu.py::test_cgl_amr_walls_after_final_field_refresh` | Inputs under `inputs/tests/`; actual restart and post-regrid state checks |
| 4 | `cgl/test_cgl_c2p_pressure_floor_cpu.py::test_cgl_production_wall_encoding_survives_independent_c2p` | Supplied HIP wall binary; real collision kernel followed by independent C2P; serialized fixture bytes do not replace a driver restart |
| 6, 7 | `L::test_cgl_lf_field_reversal_stability`; `L::test_cgl_lf_hotspot_preserves_minima_and_energy`; `cgl/test_cgl_lf_oblique_decay_cpu.py::test_cgl_lf_oblique_decay_agrees_across_blocks` | Reversal 1D/2D, 30/45-degree hotspot, global oblique decay; start safe/none |
| 7, 8 | `cgl/test_cgl_heat_flux_limiter_cpu.py::test_cgl_heat_flux_limiter_and_perpendicular_closure_endpoints`; `A::test_decay_resolves_closure_coefficients` | Reuse `CheckDiffusionSlopeMeans` and collisional/cancellation checks; coupled two-pressure decay oracle rejects the old coefficient |
| 9 | `L::test_cgl_lf_uniform_collisional_timestep`; `L::test_cgl_lf_density_contact_stability`; `L::test_cgl_lf_staggered_checkerboard_timestep_and_stability`; `L::test_cgl_lf_scaled_face_row_timestep` | Ordinary collision values, seeded/unseeded density jumps, 2D/3D CT counterexample, finite complete-row cold arithmetic |
| 9 | `L::test_cgl_lf_low_field_faces_disable_transport_cleanly`; `L::test_cgl_lf_post_sweep_timestep_refresh`; `LM::test_cgl_lf_post_sweep_timestep_refresh_agrees_across_mpi_ranks` | Zero operator and unheated/heated current-state budget |
| 10 | `L::test_cgl_lf_invalid_sts_safety_is_rejected`; `L::test_cgl_lf_sts_safety_changes_only_sts_budget`; `L::test_cgl_lf_sts_safety_is_loaded_from_restart` | Constructor, explicit-vs-STS semantics, restart parameter loading |
| 11 | In `turb/test_turb_driving_cpu.py`: `test_type_two_has_nonzero_force`, `test_conservative_forcing_kick`, `test_once_per_step_power`; `L::test_cgl_lf_paper_multicycle_forcing_work_matches_full_kicks` | Existing `tst/inputs/turb_driving_edot.athinput`; register the retained `turb_forcing` pgen for direct production tests |
| 12 | `cgl/test_cgl_lf_smr_outflow_mpi_gpu.py::test_cgl_lf_smr_outflow_mpi_gpu[4]`; `cgl/test_cgl_lf_smr_conservation_mpi_gpu.py::test_cgl_lf_smr_conservation_mpi_gpu[2-False-sts]` | Each already compares one/four ranks; other parameters run at the mesh milestone |
| 14 | `cgl/test_cgl_passive_gpu.py::test_cgl_passive_gpu`; `cgl/test_cgl_passive_wall_decode.py::test_passive_wall_survives_canonical_c2p` | Existing groups: identity, heating, linear, advection, restart, fences; `passive_acceptance.py` and `cgl_passive_validation` pgen |
| 15 | `cgl/test_cgl_lf_chunks.py`: `test_default_absent_and_zero_are_identical`, `test_uniform_collision_and_limiter_cadence`, `test_smooth_lf_temporal_refinement_and_conservation`, `test_chunked_restart_matches_uninterrupted_run`, `test_invalid_chunk_configurations_fail_before_evolution` | Prebuilt executable; small fixture dependency extraction required; failing passive checkpoint is an additional accuracy check |

The physical oracle must be independent where it matters. A copied production
Riemann expression only tests wiring; a collision/C2P roundtrip only tests
consistency. Keep analytic threshold/map checks and independent eigenvalue/decay
references alongside them. Reuse the resolved decay checker rather than adding
another permissive smoke test. Do not replay archived JSON and call it a new test.

## Performance comparisons

Freeze two workloads, not a scaling campaign:

- **Uniform control:** a populated-node active case with aligned fields. The
  historical two-node/16-GPU input used 256 cubed cells per GPU and is identified
  in `2afaa0fa1:docs/validation/cgl_lf_performance_20261008/matched-timing.json`.
  Verify the remote file/hash before claiming historical reproduction. If absent
  or the allocation differs, freeze a new one-node case and label it separately.
- **Nonzero 3D LF transport:** make a base-compatible active deck through the
  built-in pgen interface, using the intended production geometry and forcing
  from `2afaa0fa1:inputs/cgl_lf_paper/cgl_lf_physics_benchmark_matched_beta10.athinput`:
  192 by 192 by 384 cells, eight 96 by 96 by 192 blocks on one eight-GPU node,
  PPM4 with three ghost cells, beta 10, seed 271828, `dedt=0.16`, and the full
  3D solenoidal/compressive mode policy. The old revision needs its supported
  parallel-threshold spelling; do not copy later STS/passive options and assume
  they work there. Preserve the specified physical shell and forcing cadence.
  Remove output from the timing interval and verify nonzero transport separately.
  The base's built-in `cgl_lf_paper_smoke_active_beta10.athinput` verifies the pgen
  interface, but its 8 by 8 by 16 grid and different forcing policy are not this
  benchmark. The custom-pgen `cgl_lf_paper_turb_active` deck is another interface
  and planar workload; do not silently substitute it. Short startup timings are
  cost controls, not evidence of developed turbulence or long-run accuracy.

Use a common explicit initial state or a freshly produced compatible active
checkpoint. No shared full-precision active turbulent checkpoint has yet been
verified. The retained late turbulent timing checkpoint is passive J/A and must
not be loaded into pre-WO1 active E/A. New/old physics trajectories can differ;
record those intentional changes instead of treating every endpoint difference
as an implementation error or a same-work kernel comparison.

Default production timing configuration: safe arithmetic, diagnostics none,
weighted STS fluxes, strict checks on. Explicitly put these keys in the timing
input; account for `ATHENAK_CGL_LF_*` environment overrides. Hold reconstruction,
CFL, limits, forcing, block shape and ranks/nodes fixed unless they are the
explicit change under study. Historical fast/physical/non-strict runs remain
separately labeled historical controls. Do not optimize full q diagnostics.
Freeze the other diagnostic settings explicitly too: the nonzero research target
has both `mhd/cgl_lf_record_pressure_work=true` and
`turb_driving/record_injected_work=true`. Keep these for target-workload timings.
Q diagnostics off and output writers removed do not disable their traction
storage, stage reductions, or forcing-energy global reductions. An optional
false/false timing is a separately labeled core-cost control, never a replacement
for the research result. The injected-work flag is validated in restart metadata;
do not toggle it on a resumed checkpoint. Construct any alternate control from
the same initial state with its diagnostic settings fixed from the start.

Use existing cycle/time/elapsed stdout. Time the same completed-cycle interval
after warmup, without output inside it, at fixed diagnostic frequency. For an
ordinary cost-sensitive patch, one alternating baseline/candidate pair screens
the effect. Use a warmup and three alternating pairs to establish the baseline,
resolve an apparent regression, or substantiate a retained speedup. If spread
obscures the result, lengthen the interval; do not repeatedly sample for a favorable
number. Do not run competing jobs concurrently during performance measurement.

Report zone-cycles/s/node, LF RHS evaluations/cycle, elapsed/physical-time advanced,
and whole-step elapsed/RHS. The last metric is amortized total cost, not pure LF
kernel time. At final acceptance report node-seconds to the common physical
endpoint at the demonstrated accuracy. Profile separately because profiling
fences perturb timing: start with one short `cgl_lf_profile=true`,
`cgl_lf_profile_detail=false` run, using bucket call counts and `rank_max_s` to
locate expensive work or imbalance. Do not sum nested buckets or maxima from
different ranks. Detailed profiling replays arithmetic and, on the newer reference,
selects different directional kernels; it is not production kernel timing.
Use the ordinary elapsed interval as the speed authority. Compare bound cost
and LF work together: schematically, step cost is refresh/reduction overhead
plus RHS count times stage cost plus RK/forcing work. Profile only far enough
to select the next small change; do not add instrumentation for every helper.
On a fixed mesh, `lf_nstage` increments are zone-stage evaluations; divide by
active cells to obtain RHS evaluations. Obtain counts/endpoints separately if
their output changes timing. Compare accuracy at a common physical endpoint.
Reserve refinement timings for refinement/communication changes.

For equivalent optimizations, compare full-precision primary fields, relevant
ghost/coarse arrays and actual restart/RNG state on the same backend/toolchain/
rank layout. A float32 snapshot match is insufficient. Existing restart comparators
may mask documented nonphysical padding only; preserve live RNG and state, and
do not normalize unexplained differences away. Report diagnostic-only reduction
rounding separately. A physics correction may legitimately cost time; attribute
that cost and seek a simpler correct implementation without changing the equations.

## GPU handoff

The following are verified source recipes, not commands executed in this planning
session. Use the existing Frontier setup if the same platform is selected; do not
start a toolchain or compiler-flag investigation unless that setup fails.

At `2afaa0fa1`, `scripts/frontier/build_cgl_lf_frontier.sh` records the full build
contract: CPE 25.09, CCE 20.0.0, ROCm 6.4.2, Cray MPICH 9.0.1, gfx90a, and
Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517`. Its module/path cleanup is part
of the recipe. After that setup, the numerical configuration is:

```sh
cmake -S "$SRC_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release -DPROBLEM=built_in_pgens \
  -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_CXX_FLAGS="-fno-cray -mno-daz-ftz" \
  -DCMAKE_EXE_LINKER_FLAGS="-no-pie"
cmake --build "$BUILD_DIR" --parallel 16
```

Apply those flags to both Kokkos and the application. Keep per-reference build
directories; never compare a candidate against a binary with unrecorded flags.
Use `docs/validation/wo2/final/evidence/final-research/run_pytest.sh` at the pinned
ref for the complete accepted MPI/HSA/FI runtime environment. It hardcodes old
binary/output locations: reuse its setup, not its historical destinations.
The GPU MPI tests use:

```sh
export ATHENAK_RUN_MPI_GPU=1
export ATHENAK_MPI_GPU_LAUNCHER="srun --exact --threads-per-core=1 --cpu-bind=threads -c7 --gpus-per-task=1 --gpu-bind=closest"
```

Those wrappers append node/rank counts. Keep functional runs in a fresh workspace
outside the source tree. The retained
`docs/validation/wo2/final/evidence/scripts/pytest_baseline.py` stages the required
layout: import `test_suite.testutils` while in staged `tst/`, then run pytest from
staged `tst/build/src/` with `./athena` and `../../../inputs` available. There is no
repository `conftest.py` that does this automatically. Its small `launch.py` adapter
handles direct and legacy mpirun calls through Slurm without nested launches.
Adapt only source/binary/output selection as needed; do not port a campaign runner.
Avoid the generic build-and-test runners for frozen HIP binaries: they rebuild
the test tree and the generic GPU defaults select CUDA.

CPU-named application tests can use the supplied HIP executable in that layout.
The fast-speed, heat-limiter and ordinary pressure-floor helper tests explicitly
build Serial/HIP-off executables and remain CPU evidence. For the real production
wall check, build `PROBLEM=unit_tests/cgl_c2p_pressure_floor_test` under the same
HIP configuration and provide `ATHENAK_CGL_PRODUCTION_WALL_BINARY` and
`ATHENAK_CGL_PRODUCTION_WALL_LAUNCHER` (the one-rank srun prefix) to its selector.
The passive counterpart uses `unit_tests/cgl_passive_wall_decode_test` and
`ATHENAK_CGL_PASSIVE_WALL_BINARY` / `ATHENAK_CGL_PASSIVE_WALL_LAUNCHER`.
These custom pgens arrive with their owning corrections, not with the baseline.

The original minimum-reduction probe, archived uniform input, and passive failing
checkpoint are remote artifacts; their present existence is unverified here.
The compiler probe is referenced under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/power032/lf-subcycling/compiler-comparison`.
Inspect the retained script and manifest before invoking it. If it is unavailable,
make one small equivalent compiler check during baseline setup, not a new probe
matrix. Missing old timing input need not block the rebuild: use the newly frozen
common workload. Missing failing-state evidence does limit the passive accuracy
claim until that state or a demonstrated equivalent reproducer is recovered.

## Completion record and stopping point

Keep one compact row per implemented piece: commit, behavior changed, reused or
added test, command/backend, result, and timing/counts only where relevant. Keep
raw failing logs and full-precision comparison states outside Git; commit small
fixtures and a concise result summary. A tested correct baseline is enough to
proceed; there is no need for a new provenance database or exhaustive archive.

Investigate an unexplained numerical failure before stacking changes that depend
on it. Do not require a global all-green historical suite before an unrelated
focused correction can proceed. Record known B4 failures and unavailable coverage
explicitly; neither skipped tests nor old retained passes qualify this composition.
Refresh the relevant method/parameter/restart documentation with each user-visible
change; no generated manuscript/PDF rebuild or unrelated formatting pass.

The next GPU session starts with 0: verify environment and needed retained artifacts,
build the two controls, run the small baseline checks, and freeze the two workloads.
Then implement 1-4 with focused checks before moving into the LF face series.
This planning task stops before those builds, tests, fixture ports or solver edits.
