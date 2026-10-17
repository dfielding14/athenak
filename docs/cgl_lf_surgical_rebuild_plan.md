# CGL LF surgical rebuild plan

Reconstruct the required physics from pre-WO1 source, measuring each correction
and optimization separately. Keep the existing AthenaK structure. A previously
accepted physical requirement does not justify copying its entire implementation.

Branch: `c/cgl-lf-surgical-rebuild`.
Base: `8222de3aae4aaebd886653e7c61846e5f17987b4`.
Status: execution has begun on Frontier. Track completed checks, measurements,
commits, and the next piece in the [work order](cgl_lf_surgical_rebuild_work_order.md).
The numbered solver repairs remain pending until recorded there.

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
- Restore Majeski's original two-tier finite relaxation (user correction,
  2026-10-15): background decay, then one selected soft or backup branch once
  with the full outer timestep. No exact anisotropy projection at split
  boundaries, in pressure encoding, or after regridding. Threshold crossings
  remain diagnostic; nonpositive/nonfinite states remain distinct errors.
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

The collision contract follows Majeski's corrected `0e6d5e034` implementation.
With $\Delta=p_\perp-p_\parallel$, background relaxation gives
$\Delta_b=\Delta\exp(-\nu_{\rm coll}dt)$. Select the soft target/rate inside
its backup boundary, or the backup target/rate at or beyond that boundary.
Then $\Delta'=\Delta_*+(\Delta_b-\Delta_*)/(1+\nu_*dt)$ at fixed thermal
energy. Backup replaces soft relaxation; do not clamp the finite residual.
LF scattering adds background to the selected limiter rate, using inclusive
thresholds. The original AthenaK thresholds are soft $[-0.7,+0.5]B^2$ and backup
$[-1,+1]B^2$, with backup rate $10^{10}$. Select these explicitly for the current
campaign; generic parameter defaults and historical paper decks retain their
stated threshold conventions.

Dependency order: establish 0, then local repairs 1-3; 4 owns the coherent
collision model; 5-8 correct face transport; 9-10 qualify its controller; 11
corrects forcing; 12 qualifies refinement. Piece 11 can proceed after 4 without
waiting for LF tuning. Piece 13 follows the relevant corrected reference; 14
does not depend on optional optimizations or multilevel work in 12. The uniform
active milestone is a usable stopping point; add refinement or passive support
when that model is needed. Nonlinear temporal qualification applies to both
models; piece 15 provides optional chunking and a short passive setting screen.
Full passive turbulent acceptance remains with its Tier 3 milestone.
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

Completion is recorded in the work order. Run each piece's focused gate before
stacking dependent changes; batch broader integration at the milestones below.
Record existing failures rather than changing tolerances or weakening fixtures.

**Reassessment after piece 7.** Independent review of the current source and
retained evidence supports the completed production fixes; old reference commits
remain implementation clues, not scientific oracles. The collisional decay
references still compose obsolete collision maps over the entire endpoint and
check only one unsigned amplitude. Replace them in 8 with the continuum coupled
two-pressure solution and check both signed channels. The collisionless oblique
screen cannot detect this coefficient error; Tier 2 mainly samples zero or very
strong scattering. Include an intermediate-rate physical check. Retain the
piece-6 checkpoint for paired cost attribution, without treating it as a
qualified trajectory of the final model. Piece 7 established sampled minima,
conservation and fundamental-mode convergence, not full-field turbulent accuracy.

The snapshot heat-flux proxy still contains older scattering semantics. Do not
silently relabel historical proxy products as the rebuilt operator by changing
one denominator. Either reconcile its model provenance and rate conventions or
exclude it explicitly from rebuilt-model acceptance; applied-face diagnostics
remain the transport qualification. Before accepting a larger STS safety, add
the unforced finite-rate threshold-crossing refinement below. Later imports must
preserve the current initialization, collision and weak-field arithmetic fixes.

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
`CGLMHD::ConsToPrim` and `CGLRefreshPrimFromMagneticMoment`.
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

**4. Restore the original two-tier finite relaxation and its schedule.**
The 2026-10-15 user correction supersedes the earlier scheduled-wall operator.
Use Majeski's corrected `0e6d5e034` rule in `SingleCollRates_CGLMHD` and extend
its branch selection consistently to LF face scattering. Remove exact-wall
helpers, their encoding adjustments, and all walls-only hooks. Preserve
full-step rates once after RK without LF, or after the final post-LF sweep;
retain coarse state, boundary, primitive, and timestep synchronization.
AMR retains pressure positivity and conservative primary fields, without
anisotropy clipping. Finite backup residuals are diagnostic rather than fatal.
Gate: original-formula pressure/energy oracle, disabled-side and threshold
cases, actual active/passive grid encoding and reload, full-cycle cadence,
LF face-rate checks, restart, and regridding without an implicit clamp.
Earlier wall-encoding qualification remains historical evidence for the
superseded operator, not acceptance of this restored model.
Give the post-RK timestep one owner. The newer reference registers the ordinary
last-stage `MHD::NewTimeStep` and then recomputes it in `CGLCollisions`, before
the driver consumes the budget. Once the unconditional post-RK CGL state hook is
installed, omit only the redundant CGL task registration; keep `NewTimeStep`
callable for initialization and regridding. Keep the post-LF refresh and global
post-RK minimum.
This small dependent deletion can follow 4 immediately; it need not wait for 13.
Use existing heated/relaxation/restart checks and call counts to verify the ownership,
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
Use CT normal B and normalize the face vector by
`max(Bbar, |B_face|)`, where Bbar is the mean adjacent cell magnitude.
This preserves reversal attenuation while bounding the direction norm by one;
retain Bbar for focusing, moment conversion and limiter rates. The work order
records the completed early GPU qualification and pending late-state replay.
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
(`d2ab32f58`). Both safe and fast fluxes use that shared helper. Limiter rate selection
belongs to 4; do not reimplement it here. Use the intended $2\nu_{\rm eff}$
perpendicular response. With $\dot\Delta=-\nu\Delta$, this follows
[SHD97, equation 49](https://w3.pppl.gov/~hammett/gyrofluid/papers/snyder_pop97.pdf)
and differs from the printed one-$\nu$ coefficient in
[Squire et al., equation 2.7](https://arxiv.org/html/2303.00468v2#S3.SS1.SSS1).
Document that model distinction. The retained parallel 3+1 response is unchanged;
this is not a claim of exact full Braginskii transport.
Gate: collisionless identity, strong-collision limit
$\chi_\perp\nu/c_\parallel^2\to1$, intermediate collisionality, and resolved
coupled decay over an e-folding. Check both full pressure profiles and signed
Fourier channels against the independent matrix exponential, controlling spatial
and splitting errors with refinement. Preserve both thermal and magnetic-gradient
drivers. Exercise normal and logarithmic paths, including subnormal intermediates
that can yield a wrong finite result and thus evade a nonfinite-result fallback.
Do not form an overflowing literal $2\nu$ when the complete response is representable.
Retain a representable near-canceling flux when doubling the collision term
newly overflows the normalized response denominator.
Keep pieces 6-8 at a common conservative reference timestep until piece 9
qualifies the final controller. For their timing screen, enlarge the existing
kinematic `cgl_lf_rotated_decay.athinput` with an oblique field and background
coefficient: zero velocity and fixed density/B keep the pre-9 budget independent
of the evolving pressure. Verify actual dt and RHS counts; selecting explicit
integration alone does not fix dt, and there is no general fixed-dt input at the
base. This is a short patch screen, not another benchmark campaign. Do not use
unrestricted turbulence before the forcing correction to attribute kernel cost.
After piece 8, repeat the Tier 2 measurement on the same compatible checkpoint
to establish the work-count and cost baseline for pieces 9 and 10.

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
Allocate and fill that scratch on the first initialization/restart call, before
any heat-flux evaluation. Exercise the representable-sound-speed fallback in
both the face flux and bound; they must describe the same admitted state.
Match the flux's grad-B arithmetic: divide coordinate differences by spacing
before projection, then sum x/y/z in the same order. Anisotropy can amplify a
rounding residual lost by an algebraically equivalent grouping of the bound.
The exceptional logarithmic path must cover both its absolute-gradient estimate
and the represented flux gradient. Preserve subnormal background scattering in
the shared parallel response; the bound and flux must retain the same term.
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
Eigenmode amplification checks do not bound nonnormal transient growth.
The rerunnable checker is `test_cgl_lf_bound_cpu.py`; runtime row checks in
`test_cgl_landau_fluid_cpu.py` cover 2D/3D single/split layouts. Compare layout
checks at identical initial data and physical times. This temperature Jacobian
freezes density, field geometry, closure coefficients and cap scales; at zero
VL4 slopes its envelope covers neighboring branches without assuming a unique
derivative. Retain the resolved reversal, closure, hotspot and oblique gates.

**10. Separate STS safety from advective CFL.**
Only after piece 9 passes, add the parameter in `mesh.hpp`, its fresh/restart
loading in `build_tree.cpp`, and both budget selections in `mesh.cpp`;
explicit diffusion keeps its CFL multiplier. Gate: parameter bounds, analytic
stage selection, existing stability tests, and temporal convergence before
accepting a larger factor such as 0.9. Attribute changes in stage count separately
from kernel speed. Do not alter RKL coefficients or its odd-stage rule.
Both budget loops serve all registered STS processes: include an ordinary
conduction or viscosity budget check, not only LF, and retain explicit budgets.
Check zero-operator sentinels before and after safety scaling in both precisions;
the float reduction identity must not become an active finite budget after CFL scaling.
Do not adopt a larger default merely because the parameter exists. Complete the
unforced finite-rate threshold-crossing refinement below before accepting it.
Use the accepted Piece 9 executable and Tier 2 measurements as the immediate
parent. Preserve current-CFL behavior when the new safety is absent; 0.9 is not
an accepted default. Leave an absent key absent in restart parameters, so its
fallback follows the exact current CFL without serializing a rounded default;
validate an explicit override as finite and in `(0, 1]`. The ratio cap uses the
scaled STS budget and can change the outer step, so hold outer steps fixed when
isolating stage-safety effects. Adding a restart-loaded parameter requires regenerating
and relabeling the Tier 2 checkpoint. Keep paired binaries on the same compatible
checkpoint and report each binary's own work count and physical advance.

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

**Development-driver extension (2026-10-15):** selectively adopt development
`c082bea1f`'s final-stage Newtonian kick. Keep OU evolution before RK, normalize
against the live conserved state at the last stage, and reuse that stage's
exchange/recovery instead of the separate pre-RK refresh. Preserve CGL's CT
dependency for primitive restriction, passive J/A, injected-work accounting,
physical mode selection and native restart format. `accel_rms` already exists;
this changes scheduling, not the forcing options or the first-order split
accuracy. Radiation and relativistic task graphs retain their prior placement.
Also port development's remaining ion-neutral initializer fix: compute magnetic
energy from analytic local Bz, avoiding a read of an uninitialized adjacent face.
Qualify the analytic final impulse, power, thermal invariance, restart and both
AMR prolongation modes on CPU/HIP and multiple ranks. Measure the forced Tier 2
state with independent health/q-mode gates and each binary's own work counts;
cross-binary trajectories need not match. Unforced Tier 1 cannot exercise this
change. Keep the frozen milestones intact; the beta-25 campaign remains a plan.
Completed qualification: CPU/HIP jobs `5685579`/`5685580` each pass 37 focused
tests and the four-rank passive forced identity/restart gate. Tier 2 `5685581`
measures +1.14% throughput with +0.23% LF work per cycle over three alternating
pairs, with independent numerical health and bitwise q-mode equivalence. See the
work order and [driver record](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/driver-final-kick/README.md)
for the initializer fix, reporting-only recovery, exact scope and retained evidence.

**12. Repair refinement/boundary correctness before compact communication.**
Refill physical corners after prolongation (`e1b2ecd00`), then synchronize shared
LF face fluxes near refinement (`a6eef6406`) as separate fixes. Reuse the
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
Qualification also requires one boundary operator throughout each sweep: rebuild
magnetic-moment halos before the first LF RHS, then rebuild anisotropy halos
after conversion and relaxation before the next operator or timestep bound. Conservative
prolongation and the nonlinear slot conversion do not commute; mixing their
order produces an STS/explicit error plateau even as the timestep shrinks.
Reuse the existing exchanges at completed communication boundaries and do not
count entry refreshes as RHS evaluations. Add fixed-mesh STS-to-explicit pressure
convergence for both prolongation modes, plus exact shearing restart checks.
These finite checks do not prove a composed nonlinear AMR stability bound.

Time-driven uniform refinement is implemented by the generic `method=time`
criterion (`e750ff130`, repeated-handoff check `9726d30fe`): at `time_start`,
request refinement up to the existing level ceiling without requesting
derefinement. Five Serial tests pass, including original uniform LF checkpoint
loading, refinement at the final time, a later second refinement and reload.
Set `time_start=tlim` to save the refined endpoint before changing rank count;
advance both the threshold and level ceiling for the next handoff. Job `5696626`
passes small active-LF refinement on one GPU followed by an eight-rank restart
and advance; production forcing/conservation handoffs remain pending. PPM4 requires
four ghosts for AMR; primitive prolongation currently excludes pressure-work
recording. Active mesh support does not automatically qualify the later passive
model or LF chunking on AMR.

**13. Apply only measured, equivalent optimizations.**
LF-only flux-clear removal is accepted (`718901fa2`, measurements below).
Select compact IEN/IAN cell messages (`920090906`, adapted from `6539bfbbc`):
prior profiles show cell communication costs more than temperature preparation.
Limit compaction to stage > 1 of strict LF on uniform periodic meshes without
shear/orbital advection. Stage one synchronizes the full state. Owned cells and
halos then apply identical density/scalar repairs using the same frozen B and EOS;
later LF stages change only energy and magnetic moment, preserving that repair.
No extra exchange or reduction is needed. One predicate controls
receive sizes, packing and unpacking. Entry refreshes and all excluded cases
retain full messages, and flux communication is unchanged. No new option or
buffer is needed. The original compact patch passed nine CPU and nine HIP checks,
including exact native conserved/scalar/CT/ghost payloads across full/compact
exchange and one/four ranks.
Accepted Tier 1/2 comparisons improve throughput by 3.09%/8.73% against the
immediate parent, with unchanged work counts and exact developed-state native
fields/ghosts, forcing/RNG and ledgers. No-LF throughput is unchanged (+0.002%).
Tier 1 is `5656691`; Tier 2 uses qualified replays from `5656928` and timings
from `5656954` after scheduler failures interrupted the original allocation.
The work order and workspace `analysis/piece13/README.md` retain all evidence.

Defer primitive/temperature fusion (`d02bb714d`): its small measured preparation
cost does not justify duplicating EOS recovery and adding a magnetic-cache
lifetime. Defer upwind-only HLLE logarithms (`31b739890`) and speculative timestep
refresh removal. Implementing every historical optimization is not the goal.
If revisited, preserve raw/scaled B rounding, all repairs/counters and final halos;
the old logarithm patch must retain piece 5's stable weak-field expression and
neighboring-field rule. Keep only a measured benefit or independently useful
simplification, and report reduction-order differences explicitly.

**14. Add passive thermodynamics as a separate model series.**
Reuse the physical design of `a2d994a3a` and the later decoding correction
`ab9b543e7`. Scope must include EOS, fluxes, forcing, initialization, histories,
output labels, restart version/representation and support guards. The IEN slot
has a different meaning; active total-energy restarts are not passive J/A data.
Start with the implemented uniform-periodic supported domain; no speculative
AMR/boundary/source extensions. Gate: density/momentum/B identical to native
isothermal flow at matched outer timesteps and forcing; thermal evolution from
independent smooth references; pressure-work balance; finite-backup decode and restart.
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
reconstruction; preserve the selected finite-relaxation model.

**Piece 14 model qualification (2026-10-12):** source `276eee6ce` implements
this contract; [the model guide](source/modules/cgl_passive.md) records storage,
forcing, output/restart semantics and support limits. CPU `5657203` and HIP
`5657202`/`5657220`/`5657266` qualify native-flow identity, all five reconstructors,
independent heating/linear/advection references, rotated pressure work, wall
encoding and forced four-rank restart. The nonlinear reference test in `3a2a5409b`
verifies the explicit method's first-order convergence before Richardson
extrapolation; the finest STS orders are 1.999/1.994. The earlier unresolved
reference failure and corrected test-launch failure remain in the evidence.

Active control measurements (`5657269`/`5657270`) give 140.8930 to 139.6849
million zone-cycles/s/node for Tier 1 LF (-0.86%, 26 RHS/cycle), 558.4358 to
604.9033 million without LF (+8.32%), and 57.99294 to 58.07608 million for
Tier 2 (+0.14%, 50.4 RHS/cycle). Tier 2's applications and individual health,
transport and q-mode checks pass, but its inherited exact active-parent gate
fails. Local pressure differences reach 0.38% despite small global budget
residuals, with changed mirror-threshold counts. Short replay `5657582` loads
identical native fields; maximum relative parallel-pressure differences grow
from 1.50e-13 after one cycle to 1.54e-12 after four, then 0.003767 after 41.
Source review finds unchanged active equations; the specific arithmetic origin
remains unproven. Retain the failed comparison and full-precision field
differences. This is a qualified
passive model addition, not a bitwise-equivalent optimization or a demonstrated
passive speedup; turbulent threshold accuracy remains unresolved.

The retained depleted-pressure replay `5657268` still fails at safety 0.3:
cycle 10312, t approximately 4.396717, post-LF stage 25/27, one pressure-floor
event. Its old checkpoint is a stress state, not a current-model physical result.
Piece 15's short continuation and production cost screens select cap 32
without relaxing strict checks. The fresh passive Tier 3 run `5661298` subsequently
reaches t=14 with zero LF floor/invalid or EOS repair events; analysis `5661309`
completes the descriptive [6,14] comparison. The milestone is complete within
this numerical-health scope; the work order records cost and physical findings.
Workspace `analysis/piece14/README.md` retains the qualification and measurements.

**15. LF interval refinement: implemented; passive milestone setting selected.**
Hold the outer step, forcing and collision cadence fixed while refining LF
intervals on the same grid. Select settings against intended bulk observables,
budgets, health and actual work, with a checked reference pair. Report local
pressure differences descriptively. Assess the passive stress state and active
Tier 2 sensitivity separately; survival alone does not qualify an interval,
and the finest available trajectory is not an error oracle.

Source `6ea48ac69` adds about 100 production lines for optional
`time/cgl_lf_max_chunk_ratio`, default zero/off, using the existing RKL recurrence.
Each chunk is at most this ratio times the safety-scaled global LF explicit
bound; the bound is refreshed between chunks.
Logical conversion/wall/rate cadence and the zero-operator collision fallback
are preserved. Initial support is uniform periodic LF-only RKL2; AMR is excluded.
CPU `5658077` and HIP `5658068` pass 16 checks each, including smooth pressure
orders 1.95/2.01, cadence, default-off identity, restart and one/four-rank state.
A native active one-cycle comparison also matches Piece 14 exactly with chunking
off. In active one-cycle study `5658095`, caps 0, 4, 2, 1, 0.5, 0.25, 0.125,
0.0625 and 0.03125 share the outer time, forcing/RNG state and zero repair/invalid
events. Work rises from 48 RHS unchunked to 333 at cap 4 and 25,626 at the finest
cap. Maximum relative parallel-pressure differences against that finest cap are
4.522% unchunked and 0.2404% at cap 4; the two finest caps still differ by
0.2341%. Bulk differences improve, but worst-cell differences are nonmonotone:
no active chunk cap is qualified by this study.
Against that same fine reference, either pressure differs by more than 0.1% in
51,392 cells unchunked (0.363% of the domain), versus 63 cells at cap 4.

Face probe `5658098` preserves the original endpoints, but its full replay-flux
gate fails on small DIR1/2 rounding differences. One exactly reproduced DIR0
face on both runs (block 1, ghost-index k/j/i=132/41/13, post stage 13) nonetheless
confirms different scattering branches: effective rate 0 versus 1e10 and weighted
energy flux 2.31208e-7 versus 3.68616e-16. This is direct local transport evidence,
not qualification of every probe row. A same-binary single-cell +/-1e-13
perturbation at N+27 (`5658089`) produces only about 2e-11 relative pressure
response after 14 cycles. Cross-restart `5658231` isolates the actual accumulated
state difference: either executable amplifies the two N+27 states to the same
8.17e-5 relative parallel-pressure difference after one cycle; swapping only
the executable at fixed state produces at most 1.58e-13. All four zero-cycle loads
and both diagonal restarted native endpoints match the original states exactly.
This identifies state-dependent threshold amplification at the onset, without
establishing a general chaotic-growth law or the earliest roundoff source.

External SSPRK2 experiment `5658198` passes the same 16 focused HIP checks with
smooth orders 2.005-2.018. With identical initial state and outer cadence, its
finest two intervals still differ by 0.1225% in parallel pressure; finest RKL
versus SSP differs by 0.1982%. Bulk differences decrease, but local maxima remain
nonmonotone under both methods. This does not isolate RKL as the cause or qualify
an alternative integrator; SSP remains a scratch executable only.
Only 9-19 of 14,155,776 cells exceed 0.1% in either pressure across the fine
comparisons; none exceeds 1%. These descriptive thresholds are not acceptance
cuts. These findings motivate the bulk-observable criteria below; monotone
worst-cell convergence is not a universal criterion for turbulent physics.

Passive continuation `5658164` completes all eight positive caps without new LF
floor/invalid events; cap zero reproduces the original failure. Density, momentum,
CT, outer time and forcing/RNG match exactly across successful caps. Parallel
pressure differs from cap 1/32 by 4.60% locally at cap 4 and 1.08% at cap 1;
the finest pair still differs by 0.918%. Either pressure differs by more than
0.1% in 5,450 cells (0.0385% of the domain) for that fine pair, versus 15,487
cells at cap 4. Their integrated parallel/perpendicular q-gradient work agrees
to 7.11e-6/1.76e-5 relative. Pressure minima stay positive and stable; this old
deck has no event writer to audit repairs outside LF.
The passive thermal budget is a convergence diagnostic, not an exact discrete
invariant; its residual settles near -1.83e-7 of initial thermal energy.
The subsequent short continuation screen below reuses these one-cycle studies
and keeps the existing physical model and production integrator. The
[Piece 15 workspace record](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/piece15/README.md)
retains the measurements and limitations; frozen active Tier 3 artifacts remain intact.

**Production selection priority (2026-10-12):** prioritize speed and eventual
scaling, accepting local pressure differences when bulk distributions and
integrated transport remain adequate. The short passive continuations use
predeclared tolerances: 1% for pressure/thermal means and normalized temperature
distribution distance, 5% for temperature fluctuation RMS, lower pressure
quantiles and integrated work, and 0.01 for bounded-anisotropy distribution
distance. Reference-pair differences must use at most one third of these budgets.
Strict health and passive mechanical/clock invariants remain required; pointwise
maxima are descriptive. Choose the fastest passing cap using production q-off
timings; report LF stages and chunk counts as communication-work indicators,
without claiming many-node scaling from one node. The external
[continuation criteria](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/piece15/CONTINUATION_CRITERIA.md)
define quantities, matched cadence, reference refinement and the timing screen.

Job `5658763` passes all declared quantities at caps 32, 20, 12, 6, 4 and 2 over
40 matched outer cycles; caps 64 and 128 hit a pressure floor in cycle five.
All successful cases preserve native mechanics/CT, clock and forcing/RNG, with
zero LF floor/invalid and EOS/FOFC events. Cap 32 versus cap 2 differs by at most
0.00633% in pressure means, 0.0984% in temperature fluctuation RMS and 1.403%
in integrated q-gradient work. The 4/2 reference pair uses only 2.273% of a
candidate tolerance, below the one-third allowance. The accepted local tradeoff
is a 10.374% maximum, with 0.01046% of cells differing by more than 1%.
Production q-off timing `5659344` measures 24.717 versus 6.836 million
zone-cycles/s/node at caps 32/4: 3.616 times faster over three alternating pairs.
LF work in the timed interval falls from 411.286 to 119.029 RHS/cycle; chunk
counts over all 41 cycles fall from 88.098 to 12.317/cycle. All eight q-off endpoints match
their full-q physical state, forcing/RNG and non-q ledgers exactly. The job
reuses `5659159`'s full-q companions after correcting only an external runner
assertion about the intentional q-diagnostic reset; the failed job is retained.
Cap 32 at safety 0.3 also completes the fresh passive Tier 3 run to t=14 in
6.807 node-hours including I/O. This extends its clean-health evidence beyond
the short stress window; stationary and spatial turbulent statistics remain
unqualified. Active chunk qualification, AMR
support and many-node scaling remain separate; the production default stays off.

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

`lf_hardbd` counts LF-stage threshold crossings, including permitted finite
backup residuals. `lf_hwproj` remains a zero placeholder for history/restart
compatibility. No anisotropy projection is part of the current model. Check
finite/positive state and conservation separately from threshold occupancy;
zero crossings are only required for fixtures whose physics predicts them.

Reuse `tst/test_suite/`, `tst/scripts/cgl/`, `test_suite.testutils`, existing
inputs and production pgen fixtures. A gate is evidence for the changed behavior,
not a required new test file. Do not import whole multi-thousand-line test modules
to obtain one regression. Add their needed functions/fixtures to the existing
structure with the owning patch. No whole-repository cleanup or new test framework.
Keep validation proportional to the change: reuse accepted evidence unless a
relevant source/configuration change or concrete failure invalidates it. New
diagnostics need a specific unanswered question and a stopping condition; retire
them from the routine workflow once answered. Do not add repeated audit layers,
duplicate readers or new provenance machinery merely to increase confidence.

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
and the Tier 3 matched active run to t=14 with nonzero LF work. Run an actual HIP
finite-backup encoding fixture; CPU scalar helpers cannot substitute for it. Run the
small safe/fast and diagnostic/profile branch checks once after all face changes,
not after each constant edit. Production timings keep detailed q diagnostics off.
Use the existing 64-by-64 hotspot with `lf_coefficient_mode=local`,
kinematic/advect evolution, zero velocity and no forcing/collisions. Its accepted
temporal window is `sts_max_dt_ratio=8,4,2,1` at fixed grid and `tlim=0.201`.
Preserve the original minima/conservation checks and at-least-five nominal
diffusion-time duration. Compare both full-precision pressure fields; report
actual steps/RHS counts, norms and worst-cell absolute/relative and
thermal-energy-scaled differences. Distinguish the five-to-three-stage transition
from the fixed-stage 4,2,1 trend. This checks nonlinear LF and representation
conversions; it does not certify turbulent threshold switching.

The initial 1-to-0.125 ladder is retained as a precision stress study after
finer steps failed the inherited minimum gate in the barely represented Gaussian
tail, and the finest rung also failed energy/moment conservation. Full-pressure
differences still decrease near second order; these are tighter invariant and
minimum failures, not a plateau in the pressure norms. The prospectively fixed
8,4,2,1 window passed every original gate in job `5647113`; fixed-stage indicators
are 1.946–1.992. The 0.25-to-0.125 differences are 1.56–1.66% of the 2-to-1
differences in all retained norms, below the prescribed 10% reference-sensitivity
limit. This is bounded convergence, not a rigorous error bound or a resolution
of the finer failures. The B=1 conditioning control activated the unconditional
fluid-firehose wall and remains invalid; the single B=5 replacement (`5647008`)
passed its diagnostic comparison. Retain both controls and the failed finer
ladder as evidence; no further conditioning sweep is planned.

Before accepting a larger STS safety, also use the existing `limiter_stress`
fixture for one smooth, unforced, finite-rate threshold-crossing case. Keep the
grid, endpoint, physical rates and thresholds fixed across at least three decreasing
outer steps; compare both full pressure fields, means/energy and repair activity.
Verify that face-averaged states actually cross an enabled soft threshold during
the run. Require decreasing errors rather than second order: the once-per-cycle
collision composition is not symmetric. This small full-model test complements
the LF-only hotspot; it does not qualify a forced turbulent trajectory. If an
LF-only refinement at fixed outer step is needed, bring forward conditional 15.

At this milestone, compare once against a qualified executable from `2afaa0fa1`
with the same active physics/settings, merging and chunking off. Reuse a verified
binary or make one additional build, not another full validation matrix. The
pre-WO1 controls locate historical costs; this corrected reference answers
whether the smaller reconstruction improves on the newer solver. Use physical
endpoint checks for intentional arithmetic differences. Exact equivalence of
each optimization remains a comparison with its immediate rebuilt parent.

**Active milestone complete (reviewed 2026-10-12):** the physical suite is
frozen in `build/active-milestone-final` with 33 passes and two known B4 expected
failures. The bounded corrected-reference comparison (`5647118`) is complete;
its near-threshold pressure differences remain descriptive, not an accuracy
oracle. Reuse these results and keep STS safety 0.3.

1. External workflow consolidation is complete. The workspace's
   [workflow README](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/active-milestone/README.md)
   lists the existing Tier 1, Tier 3 and analysis entrypoints and required helpers;
   the work order links the concise milestone report. Other diagnostic/setup/review
   scripts are archived in place. Essential regressions, optional endpoint export,
   frozen inputs/binaries, raw results and failed cases remain intact. Preserve
   provenance paths and reuse working tools; no further consolidation is needed.
2. Tier 1 is complete: job `5647205` passed all 14 invocations in one two-node
   allocation, including three LF/no-LF pairs and the historical control. The
   milestone report records rates, spread and work/cost per RHS. Cancelled job
   `5647148` remains separate partial evidence.
3. Canonical Tier 3 job `5647225` completed fresh to t=14 at safety 0.3:
   47,681 cycles, 3.475833 node-hours including I/O and 47.465 LF RHS/cycle.
   LF floor/invalid and EOS repair counters stayed zero; the [6,14] energy minus
   applied-work residual is 2.23e-11 relative to injected work. Analysis job
   `5647226` completed all four figure groups, now reviewed in the milestone
   report. The scientific interpretation remains inconclusive: one heated
   realization does not establish turbulent convergence or active/passive
   causality. This is evidence for frozen pre-cleanup source `73e2769fc`;
   subsequent cleanup and Piece 12 qualifications remain separate.

Do not repeat completed physical/reference checks or add another review layer
without a relevant change or concrete failure. Investigate a new failure with
the smallest discriminating check. Keep existing acceptance tolerances and
required health checks; consolidation itself needs no new simulation campaign.
The active milestone is closed within that scope; mesh, passive and optional
optimization work retain their separate qualification.

**Authorized code cleanup completed (2026-10-10):** dead-code and environment-override/
replay removal, directional consolidation and LF-only flux-clear removal are accepted.
The work order records gates and separate paired costs: consolidation retains a
3.29% uniform penalty with neutral developed-state cost; clear removal improves
Tier 1/2 throughput by 4.10%/2.95%, with identical native endpoints and work counts.
It records the then-existing shearing restart-counter failure, subsequently
resolved by Piece 12's boundary lifecycle repair.
The frozen Tier 3 executable and analysis dependencies remain intact; both jobs
completed successfully. Defer workspace deletion and additional physics/mode changes.

**Mesh milestone:** after 12, cover both physical-corner cases, all six retained
conservation parameters, actual regrid/restart, repair accounting and shearing.
One/four ranks are the default comparison. Use two-node/16-rank 3D churn only at
final qualification when the allocation supports it, not as a prerequisite to
every local fix. These checks gate multilevel/communication acceptance; they do
not prevent independent uniform-grid work while an unrelated mesh issue is diagnosed.

**Piece 12 completed (2026-10-10):** final source `c8d57ea61` passes 25 checks on
each of CPU/HIP (`5647453`/`5647452`) and both final MPI checks (`5647477`, including
two-node/16-rank churn). Physical-corner refill, shared-face flux reconciliation
and consistent A/moment boundary ordering are qualified together. Both
prolongation modes converge toward refined explicit pressure solutions; native
periodic energy/moment drift is below `7e-16`. The native shear restart matches
all conserved/CT/ghost values and integer counters exactly, including 1522 RHS
evaluations. Its parallel q-work sum differs by one ulp because restart restores
the global accumulator on rank zero; the extra all-floating-bitwise assertion
failed and is retained as such. No solver change is needed for that diagnostic
regrouping. The work order and workspace `analysis/piece12/README.md` retain
evidence and limits. The frozen Tier 3 artifacts remain unchanged.

**Passive milestone:** start with identity and independent thermal evolution;
then run the six existing groups once on HIP, including all supported
reconstructors. Check forced identity/restart on multiple ranks. Do not repeat
the historic hundreds-of-application CPU/HIP/MPI cross-product after each edit.
Require a short nonlinear temporal comparison in the supported passive model;
reuse the smooth refinement and retained depleted-pressure state as applicable.
Keep forcing and rate cadence matched for an LF-only comparison. Reducing the
outer timestep in a forced run changes that composition and is not an isolated
LF refinement. Piece 15 now provides independent LF interval refinement;
a small smooth test or clean termination does not establish turbulent accuracy.
Do not relabel the active kinematic fixture as a passive
test: the current passive guard requires dynamic evolution.
The Tier 3 passive member `5661298` completed t=14 using Piece 15's selected cap
32 at safety 0.3. Analysis `5661309` covers [6,14], including the retained spatial
figures. This closes the milestone after Piece 14 qualified the model but
reproduced the old depleted-pressure failure without chunking. The separate
workspace record `analysis/piece14/TIER3.md` gives health, cost and finite-window
physical findings. The frozen active artifacts remain intact; active chunking
is off, so this is not a comparison at matched LF settings.

For these scientific milestones, read the existing LF floor/invalid counters,
endpoint pressure extrema, and EOS/FOFC event log. Q diagnostics can stay off;
strict LF admissibility rejects invalid states. Finite density and pressure-floor
repairs are counted without aborting (user decision, 2026-10-15); nonfinite raw
density, energy, magnetic moment and decoded LF pressures must be rejected before
repair can hide them. Existing `lf_dfloor`/`lf_pfloor` restart counters are retained;
the per-cycle EOS event log also covers repairs outside LF. Investigate changed
or dominant repair activity before treating the result as accurate physics, without imposing
a universal zero-repair criterion on every turbulent case. Counts may include
ghosts and repeated refreshes; they are not distinct-cell fractions or mass/energy
budget. Read these in functional runs; add no hot-path instrumentation or new
logging system.

Use inherited scientific tolerances unless the changed equation requires a new
analytic expectation. Final-reference tests sometimes hard-code safety 0.9,
13 heated stages, zero stages for a zero operator, or a timestep cap to hit an AMR
event. Derive the expectation for the controller currently being tested. Preserve
the physical error, conservation, amplitude, topology and strictness assertions;
do not copy a cap blindly or weaken a criterion just to recover a pass.

## Concrete test map

The original selectors below were verified in Git source; completion records
identify the imported checks and actual runs. Paths are
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
| 4 | `L::test_cgl_collision_rates_once_per_cycle_with_and_without_lf`; `L::test_cgl_lf_face_collision_rates_are_tiered`; `L::test_limiter_stress_matches_cellwise_relaxation`; `L::test_limiter_stress_matches_finite_backup_relaxation` | `cgl_collision_once.athinput` and existing mirror/firehose unit inputs; independent rate/cadence oracles |
| 4 | `L::test_cgl_sts_viscosity_finishes_with_admissible_pressures` | No-LF viscous heating across the fluid-firehose threshold without clipping; positive pressures, shear decay and conserved totals. The cadence selector also checks that non-LF STS checkpoints add no rate kicks |
| 4 | `L::test_cgl_collision_refreshes_next_timestep_from_relaxed_state`; `L::test_cgl_lf_restart_with_finite_collision_preserves_corrected_split`; `cgl/test_cgl_amr_walls_cpu.py::test_cgl_amr_refresh_preserves_finite_anisotropy` | Inputs under `inputs/tests/`; actual restart and post-regrid state checks |
| 4 | `cgl/test_cgl_c2p_pressure_floor_cpu.py::test_cgl_finite_backup_survives_independent_c2p` | Supplied HIP collision binary; active/passive collision kernels followed by independent C2P; serialized fixture bytes do not replace a driver restart |
| 6, 7 | `L::test_cgl_lf_field_reversal_stability`; `L::test_cgl_lf_hotspot_preserves_minima_and_energy`; `cgl/test_cgl_lf_oblique_decay_cpu.py::test_cgl_lf_oblique_decay_agrees_across_blocks` | Reversal 1D/2D, 30/45-degree hotspot, global oblique decay; start safe/none |
| 7, 8 | `cgl/test_cgl_heat_flux_limiter_cpu.py::test_cgl_heat_flux_limiter_and_perpendicular_closure_endpoints`; `A::test_decay_resolves_closure_coefficients` | Reuse `CheckDiffusionSlopeMeans` and collisional/cancellation checks; coupled two-pressure decay oracle rejects the old coefficient |
| 9 | `L::test_cgl_lf_uniform_collisional_timestep`; `L::test_cgl_lf_density_contact_stability`; `L::test_cgl_lf_staggered_checkerboard_timestep_and_stability`; `L::test_cgl_lf_scaled_face_row_timestep` | Ordinary collision values, seeded/unseeded density jumps, 2D/3D CT counterexample, finite complete-row cold arithmetic |
| 9 | `L::test_cgl_lf_low_field_faces_disable_transport_cleanly`; `L::test_cgl_lf_post_sweep_timestep_refresh`; `LM::test_cgl_lf_post_sweep_timestep_refresh_agrees_across_mpi_ranks` | Zero operator and unheated/heated current-state budget |
| 10 | `diffusion/test_sts_safety_cpu.py`; `diffusion/test_parabolic_budget_cpu.py` | Explicit parameter bounds, exact absent-CFL fallback, initial/refreshed STS versus explicit/advective budgets, real restart, ordinary viscosity and inactive sentinels in both precisions |
| 10 | `cgl/test_cgl_lf_safety_refinement_cpu.py`; `cgl/test_cgl_lf_geometry_cpu.py::test_cgl_lf_sts_reversal_converges_to_explicit` | Full pressure fields at safety 0.3/0.9; finite-rate soft-threshold crossing and reversal against refined references. Scale the ratio cap inversely with safety to match nominal outer-step ladders; verify actual steps and work |
| 11 | In `turb/test_turb_driving_cpu.py`: `test_type_two_has_nonzero_force`, `test_conservative_forcing_kick`, `test_once_per_step_power`; `L::test_cgl_lf_paper_multicycle_forcing_work_matches_full_kicks` | Existing `tst/inputs/turb_driving_edot.athinput`; register the retained `turb_forcing` pgen for direct production tests |
| Active milestone | `cgl/test_cgl_lf_waves_cpu.py`; `cgl/test_cgl_lf_closure_cpu.py`; `cgl/test_cgl_lf_transverse_cpu.py::test_cgl_lf_local_hotspot_temporal_refinement` | Qualified wave/e-fold checks and local-hotspot window 8, 4, 2, 1 (`5647113`); fixed grid/endpoint, full-precision pressures, actual work and decreasing adjacent differences. Retain the finer 1-to-0.125 precision stress failures and reference-sensitivity evidence; no rung is an exact oracle |
| 12 | `cgl/test_cgl_lf_smr_outflow_mpi_gpu.py::test_cgl_lf_smr_outflow_mpi_gpu[4]`; `cgl/test_cgl_lf_smr_conservation_mpi_gpu.py::test_cgl_lf_smr_conservation_mpi_gpu[2-False-sts]` | Each already compares one/four ranks; other parameters run at the mesh milestone |
| 13 | `cgl/test_cgl_lf_compact_mpicpu.py::test_cgl_lf_compact_exchange_matches_full_exchange` | Existing nonuniform 3D scalar fixture; exact full native payload across strict/full exchange and one/four ranks, positive pressures, nonzero LF work and zero repairs; run on CPU/HIP |
| 14 | `cgl/test_cgl_passive_gpu.py::test_cgl_passive_gpu`; `cgl/test_cgl_c2p_pressure_floor_cpu.py::test_cgl_finite_backup_survives_independent_c2p` | Identity, heating, linear, advection, restart, fences and nonlinear temporal refinement; `passive_acceptance.py` and `cgl_passive_validation` pgen; forced four-rank identity/restart through `--suite forced` |
| 15 | `cgl/test_cgl_lf_chunks.py`: `test_default_absent_and_zero_are_identical`, `test_uniform_collision_and_zero_operator_cadence`, `test_smooth_chunk_refinement_and_conservation`, `test_chunked_restart_matches_uninterrupted`, `test_chunked_mpi_state_and_schedule`, `test_invalid_chunk_configurations` | Prebuilt executable; existing passive and conservation fixtures; active threshold and failing passive checkpoint refinements remain separate accuracy checks |

The physical oracle must be independent where it matters. A copied production
Riemann expression only tests wiring; a collision/C2P roundtrip only tests
consistency. Keep analytic threshold/map checks and independent eigenvalue/decay
references alongside them. Reuse the resolved decay checker rather than adding
another permissive smoke test. Do not replay archived JSON and call it a new test.

## Performance comparisons

Use three frozen benchmark tiers. Each answers a different question; none
substitutes for another. Weak/strong scaling sweeps, other MeshBlock sizes and
higher-resolution boxes, including 384 cubed, belong to the separate proposal
`CGL/scaling/PRODUCTION_LIKE_SCALING_EXTENSIONS.md` in the project workspace;
this plan runs none. Only Tier 1 requires two nodes; Tiers 2 and 3 use the same
one-node matched box. Frozen inputs live under the rebuild workspace's `inputs/`.

| Tier | Workload | Allocation | Answers | Run when |
| --- | --- | --- | --- | --- |
| 1. Uniform cost control | Archived scaling geometry and uniform state, 40 cycles | 2 nodes / 16 GPUs, one 256 cubed block/GPU | Cost per LF evaluation; LF overhead relative to no-LF CGL | Any change that adds, removes or alters recurring full-domain work |
| 2. Developed-state restart | Frozen active matched-box checkpoint near t=6; 40-cycle window plus one untimed guard cycle | 1 node / 8 GPUs, eight 96 x 96 x 192 blocks | LF evaluations/cycle and bound, limiter and strict-check cost in a nonuniform state | Controller, STS-safety, LF-interval, wall/limiter or state-dependent kernel changes |
| 3. Matched turbulence box | Matched box from its initial state to t=14 | Same as Tier 2 | Physical behavior and node-hours to the common endpoint | Active and passive milestones and final acceptance |

**Tier 1, uniform cost control.** Use the frozen `uniform-2n.athinput`: the
archived two-node geometry (256 x 1024 x 1024 cells, sixteen 256 cubed blocks,
domain 1 x 4 x 4) and unforced uniform state (rho=1, p_parallel=p_perp=0.6,
Bz=sqrt(0.12), beta 10) with production LF modes. Freeze a no-LF companion
that differs only in basename, `sts_integrator=none` and removal of the heat-flux
keys, matching the archived collisionless variant. Run both in the same
allocation; report the LF/no-LF throughput ratio and absolute rates.
Freeze and hash the no-LF companion in the workspace `inputs/README.md`, and
check it on both controls before counting it as a qualified control.
Keep `uniform-historical-2n.athinput` as a separately labeled historical
control: its bytes, including the original 20-cycle limit and LF modes, match
the archived `cglws_lf_sts_0002n` input (SHA-256
`d1e06b2a9cc1b9d9b1dcea489211d2fa13a8cf386ca37cec34491fc5320ea763`).
Run it when establishing a baseline and at milestones, not for every patch.
Do not shrink this tier below two nodes or change its block size; continuity
with the archive depends on that layout.

This state has zero heat flux and a grid-aligned field. It exercises LF kernels,
exchanges, reductions and the timestep row loop, but does not activate gradient
terms in the bound, heat-flux caps, low-field branches or limiters. It controls
cost per evaluation; its evaluation count does not predict turbulent work.
A change that alters Tier 1 LF evaluations/cycle has changed the controller;
report it and judge the resulting work and cost on Tier 2.

**Tier 2, developed-state restart.** Produce a full-precision active checkpoint
of the Tier 3 box from the prescribed initial state with the current candidate,
q diagnostics off and both work ledgers on from the start. Target t=6, the start
of the analysis window. Produce the checkpoint with the canonical writer blocks
present, so replays can enable or disable them through existing keys. Record the
actual checkpoint time and cycle, producer revision, input/binary/checkpoint
hashes. A periodic restart near t=6 is acceptable
if its actual time is recorded. This is a frozen cost-screening state, not a
physical result or evidence of stationarity. Attempt and qualify it before piece
7. If a pending correction prevents reaching it healthily, retain the failure,
record Tier 2 as unavailable, and qualify the checkpoint after that correction;
do not relax checks or block the needed fix on this prerequisite.

Every invocation loads the same checkpoint. For restart cycle N, use
`time/nlim=N+41`, `time/ndiag=1` and a `time/tlim` beyond the whole interval.
Measure between the stdout timestamps at N+5 and N+40: 35 completed updates
after five warmup cycles. The extra untimed cycle makes N+40 an ordinary loop
diagnostic; final output and problem-finalization work occur before the final
diagnostic and would otherwise contaminate it. Override only keys already in
the checkpoint's parameter record. Disable each existing output writer with
`dcycle=0` if that key is present, otherwise `dt=0`, and verify no files are
written. Delaying the next output alone does not suppress unconditional final
output. Preserve LF diagnostic mode, both ledger flags, forcing parameters and
saved forcing/RNG state. Verify zero new LF floor/invalid events in a separate
short history run from the checkpoint; inspect EOS/FOFC repairs with an event
writer included in the producer deck. With q diagnostics off, the q-work fields
are zero, and stage counts alone do not establish nonzero transport. For that
qualification only, allow a separate untimed full-q replay from the checkpoint;
verify nonzero q work and full-precision physical/forcing endpoint agreement with
a matching q-off replay. Keep the frozen checkpoint and all timing runs q-off.
Obtain LF work counts over the same N+5 to N+40 window in a matching history
replay: delta(lf_nstage)/active_cells gives RHS evaluations, then divide by 35
for evaluations/cycle. In the history replay, set the existing history `dt`
below the smallest cycle timestep (for example, 1e-12) so every cycle is recorded;
verify rows exist at both timing endpoints. Whole-run counts include warmup and
the guard cycle and can bias cost/RHS when stage counts vary.

Parent and candidate must load the same checkpoint with verified compatible
decoding. Regenerate and relabel it after changes to restart contents or newly
loaded parameters, or when accepted physics changes make the old state
unrepresentative. Keep old checkpoints and timings, and remeasure both members
on the replacement; if they cannot read a common checkpoint, establish a new
baseline without claiming a paired speedup across that transition. A passive
J/A checkpoint cannot initialize active E/A. The passive model gets its own
Tier 2 checkpoint from its qualified implementation. Pieces 7, 9 and 10,
bound/refresh optimizations in 13, and any chunking in 15 report work count and
cost here together. Repeat Tier 2 after piece 8 to set the baseline for 9 and 10.

**Tier 3, matched turbulence box.** The canonical Majeski-type workload is
`2afaa0fa1:inputs/cgl_lf_paper/cgl_lf_physics_benchmark_matched_beta10.athinput`.
Derive a separate milestone deck from the geometry and initial state frozen in
`active3d-1n.athinput`; preserve that 40-cycle startup control unchanged:

- Grid: periodic (1,1,2), 192 x 192 x 384, eight 96 x 96 x 192 blocks on one
  eight-GPU node; use the built-in `cgl_lf_paper` turbulence interface.
- Physics and numerics: PPM4 with three ghosts, RK2/RKL2, CFL 0.3, local LF
  coefficients, beta 10, rho=1, B=(0,0,1), velocity=0, p_parallel=p_perp=5;
  retain the canonical forcing, floors, thresholds and finite-rate limiter policy.
- Forcing: seed 271828, `dedt=0.16` per unit volume (nominal total power
  0.32 in volume two, matching the stated MKS24 total), tcorr=2,
  dt_update=0.01, the full 3D solenoidal/compressive amplitude blend
  1/(1+sqrt(2)) and the physical shell pi <= |k| <= 3*pi.

Start fresh with the qualified milestone solver and run to t=14 with no cycle
cap; the pre-piece-7 cost checkpoint cannot supply this trajectory. Restore the
canonical history (0.02), field (0.25, both `mhd_w_bcc` and `turb_force`) and
restart (1.0) cadences, verifying each writer exists in the candidate. The
analysis window is [6,14]. The `cgl_lf_paper_standard_*` MKS24 decks (32 x 32 x 64
blocks, PLM, `dedt=0.32` per volume) are a different benchmark; neither those,
the small smoke deck nor the custom-pgen planar workload substitutes for this box.
Use explicit `sts_safety=0.3` for the active milestone. Piece 10 qualifies the
separate parameter and small 0.3/0.9 refinement cases, but its Tier 2 cost screen
does not establish equal turbulent accuracy. The WO2-era matched deck used
`sts_safety=0.9`, which the frozen startup copy omitted; selecting that value
requires milestone accuracy criteria and evidence first. Freeze the qualified
value for the matched pair and label historical differences.
Run the active member at the active milestone and the passive member at the
passive milestone. Restore canonical `iso_sound_speed=sqrt(5)`, omitted from the
active startup control, and set `mhd/passive` and `problem/passive_delta` together
for each member. Half-sweep merging stays off; chunking requires piece 15. If
passive qualification requires different numerical settings, label them separately.
A numerically matched-pair claim requires independently qualifying and rerunning
active at the common settings; it is not a prerequisite for the passive milestone.

Report:

- node-hours to t=14, including I/O;
- interval rates from periodic `zone-cycles/s=` diagnostics, identifying
  intervals containing I/O and excluding the final diagnostic from evaluation-cost
  comparisons;
- LF evaluations/cycle over time;
- the active energy budget, Delta E minus accumulated injected work;
- LF floor/invalid counters, pressure extrema and EOS/FOFC repair activity;
- the retained figure groups from
  `2afaa0fa1:docs/cgl_lf_matched_physics_benchmark.md` and its retained reports.

The figures are descriptive comparisons with retained references, not a
physical-validation pass criterion. Surviving to t=14 does not establish accuracy.
There is no cooling: thermal energy and beta may evolve, and t=6 does not by
itself establish stationary turbulence. Keep intentional model/trajectory changes
distinct from equivalent-optimization comparisons. Once a milestone run is
accepted, its restart near t=6 may become the next Tier 2 checkpoint.

The 40-cycle `active3d-1n.athinput` startup screen remains an optional quick check
of hyperbolic, forcing and LF evaluation cost on the target geometry. It cannot
replace Tier 2 work counts or Tier 3 physical-time cost and accuracy evidence.

Tier 1 takes a few minutes on two nodes per invocation. Piece 10's checkpoint
producer reached t=6.05 in 1.371 node-hours including I/O (1.400 for the entire
allocation). At safety 0.3, its Tier 2 candidate uses 50.4 LF evaluations/cycle
and about 863 seconds per simulation-time unit near t=6. Timing invocations take
seconds, while the separate full-q qualification replays take about 90 seconds.
The earlier Tier 3 estimate of two node-hours assumed startup work near 38 evaluations
per cycle; budget more for the milestone and measure its actual node-hours to
t=14. A short developed-state interval does not fix the cost of the trajectory.

Default production timing configuration: safe arithmetic, diagnostics none,
weighted STS fluxes, strict checks on. Explicitly put these keys in the timing
input. LF mode environment overrides are removed; clear them when comparing
older binaries. Hold reconstruction,
CFL, limits, forcing, block shape and ranks/nodes fixed unless they are the
explicit change under study. Historical fast/physical/non-strict runs remain
separately labeled historical controls. Do not optimize full q diagnostics.
Freeze the other diagnostic settings explicitly too: the nonzero research target
has both `mhd/cgl_lf_record_pressure_work=true` and
`turb_driving/record_injected_work=true`. Keep these for target-workload timings.
Q diagnostics off and output writers removed do not disable their traction
storage, stage reductions, or forcing-energy global reductions. Tier 1 is the
labeled core-cost control with no forcing and pressure-work recording off;
Tiers 2 and 3 keep both work ledgers on. The injected-work flag is validated in
restart metadata; do not toggle it on a resumed checkpoint. Fix diagnostic and
ledger settings from the initial state.

Use existing cycle/time/elapsed stdout. Time the same completed-cycle interval
after warmup, without output inside it, at fixed diagnostic frequency. For an
ordinary cost-sensitive patch, one alternating baseline/candidate pair screens
the effect. Use a warmup and three alternating pairs to establish the baseline,
resolve an apparent regression, or substantiate a retained speedup. If spread
obscures the result, lengthen the interval; do not repeatedly sample for a favorable
number. Do not run competing jobs concurrently during performance measurement.
The alternating-pair protocol applies to Tiers 1 and 2; Tier 3 is one run per
member at its milestone, repeated only for a changed configuration or a failed gate.

Report zone-cycles/s/node, LF RHS evaluations/cycle, elapsed/physical-time advanced,
and whole-step elapsed/RHS. The last metric is amortized total cost, not pure LF
kernel time. At final acceptance report node-hours to t=14 at the demonstrated
accuracy. Profile separately because profiling fences perturb timing: start with
one short `cgl_lf_profile=true`,
`cgl_lf_profile_detail=false` run, using bucket call counts and `rank_max_s` to
locate expensive work or imbalance. Do not sum nested buckets or maxima from
different ranks. Detailed replay profiling is removed from the rebuilt solver;
`cgl_lf_profile_detail=false` remains accepted for existing decks. Historical
replay timings are not production kernel timings.
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
finite-backup check, build `PROBLEM=unit_tests/cgl_c2p_pressure_floor_test` under
the same HIP configuration and provide `ATHENAK_CGL_COLLISION_BINARY` and
`ATHENAK_CGL_COLLISION_LAUNCHER` (the one-rank srun prefix) to
`test_cgl_finite_backup_survives_independent_c2p`. This fixture covers active
and passive encoding through separate C2P/reload and LF representation changes.

The archived uniform input is verified: the frozen historical copy matches its
SHA-256, recorded in the rebuild workspace's `inputs/README.md` and above.
The compiler probe is referenced under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/power032/lf-subcycling/compiler-comparison`.
Inspect the retained script and manifest before invoking it. If it is unavailable,
make one small equivalent compiler check during baseline setup, not a new probe
matrix. Missing failing-state evidence does limit the passive accuracy
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

Follow the work order for the current next step and retained completion evidence.
Piece 0 describes baseline setup; this benchmark-plan update does not require
repeating completed solver qualification or changing the frozen startup controls.
