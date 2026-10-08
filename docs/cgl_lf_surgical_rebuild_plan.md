# CGL LF surgical rebuild plan

Reconstruct the required physics from pre-WO1 source, measuring each correction
and optimization separately. Keep the existing AthenaK structure. A previously
accepted physical requirement does not justify copying its entire implementation.

Branch: `c/cgl-lf-surgical-rebuild`.
Base: `8222de3aae4aaebd886653e7c61846e5f17987b4`.
Status: branch and plan created; no solver changes or new validation runs.

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

## Implementation sequence

All implementation pieces are pending. Finish the relevant correctness and
measurement gate before advancing; record existing failures rather than making
them disappear by changing tolerances or fixtures.

**0. Establish the baseline before changing production source.**
Build the exact base and post-WO1 reference with identical qualified compiler,
Kokkos, precision and runtime settings. Reproduce the reduction correctness probe
on CCE20/HIP; retain the documented corrective flags on both builds. Record the
baseline tests and paired timings described below. Confirm the existing fast
speed against non-unit-density anisotropic limits and independent eigenvalues.
Use current source to locate the real test entry points; old work-order line
numbers and several original findings are stale.

**1. Keep unsupported models out of the accepted domain.**
Add the small passive-mode rejection until its redesign is complete; preserve
the existing LF inflow/user-boundary guard. Validate CGL magnetic floors and
collision/limiter parameters at their existing constructors. Add the units
requirement if the implemented Spitzer branch is selected; do not add a new
conduction model. Scope: EOS initialization and affected input/startup tests.
Gate: invalid and unsupported inputs fail clearly; valid active cases are unchanged.

**2. Repair the remaining EOS floor defects individually.**
Scope: `src/eos/ideal_c2p_mhd.hpp`, `src/eos/cgl_mhd.cpp`, and existing floor tests.
Preserve pressure anisotropy under density floors, enforce the magnetization
ceiling through the density floor in both A and magnetic-moment decoders, and
make repaired states representable and idempotent. Do not redo the already-correct
perpendicular energy factor. Gate: analytic pressure ratios and energy identity, repeated C2P with no
second change, changed A persisted, finite/extreme inputs, focused double/float
checks. A focused float test does not establish whole-application float support.

**3. Repair reconstruction and FOFC detection.**
Scope: all directional PPMX and WENO-Z wrappers and the C2P test-floor path;
leave PPM4 unchanged.
Floor both reconstructed CGL pressures; detect nonfinite incoming state before
repair obscures it. Gate: actual reconstruction wrappers, injected NaN/Inf in
E/A, a steep-pressure evolution, and ordinary MHD controls. Change matching
reference consumers only where their old expectation is demonstrably wrong.

**4. Restore thresholds, relaxation, walls and their schedule coherently.**
Scope: `EOS_Data`, `cgl_physics.hpp`, EOS collisions, LF scattering, MHD task
callers, AMR projection, affected diagnostics and parameter fixtures.
Introduce the agreed soft/backup thresholds and additive scattering, then land
the monotone relaxation map with its full-step schedule. Keep coarse state and
primitives consistent after updates. Incorporate the later canonical active wall
encoding correction here (`71ad25ebc`, regression `4cd46a131`), rather than
reintroducing its known restart failure.
Gate: analytic background relaxation with/without LF; an independent finite-rate
map distinguishing one full kick from two half kicks; threshold equality and
additive-rate checks; invariant E/momentum/B; unchanged A on no-op; wall
idempotence; a separate production C2P after collision encoding; real restart.
Measure added full-domain passes and timestep refreshes. Remove a superseded
refresh only in a separate equivalence-checked optimization.

**5. Restore the small weak-field face-transport correction.**
Scope: HLLE/LLF face helpers and their direct/evolution fixtures (`7d1f79557`).
Retain conservative A; use the appropriate neighboring field at the weak-field
face. Gate: first-update contamination, both contact orientations, and smooth
fixed-width interface convergence. Preserve the two documented B4 strict
sharp-contact expected failures with their original predicates. This patch does
not solve that formulation's sharp-contact limitation.

**6. Correct LF face geometry.**
Scope: `cgl_landau_fluid.cpp/.hpp`, existing face-field callers, and reference
consumers (`910023601`). Use CT normal B and the mean of adjacent cell magnitudes
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

**8. Correct the perpendicular BGK response and consistent scattering.**
Scope: safe/fast closure arithmetic, logarithmic fallback, and independent
reference consumers (`d2ab32f58`, `a8a25fb68`). Use the intended $2\nu_{\rm eff}$
perpendicular response. Gate: collisionless identity, strong-collision limit
$\chi_\perp\nu/c_\parallel^2\to1$, and resolved coupled decay over an e-folding.
Keep pieces 6-8 at a common conservative explicit-reference timestep until
piece 9 qualifies the final controller. Do not use unrestricted turbulence runs
at these intermediate commits to attribute a kernel performance regression.

**9. Implement the LF bound for the final stencil.**
Scope: timestep calculation/callers, existing post-RK refresh, and MPI reduction.
Use the mathematical work in `822cad5e8`, including CT geometry, neighboring
density, VL4 derivatives and grad-B coupling. Reuse available temperature/B
calculations where valid; avoid the later duplicate cache passes. Include
background collisions. Do not first port WO1's provisional D5 bound merely to
replace it, and do not loosen unrelated conduction bounds.
Gate: independent frozen-stencil rows/eigenvalues, staggered checkerboard,
unequal spacing, density-jump seeded-minus-unseeded growth, reversal/grad-B,
heated post-RK refresh, MPI stage agreement, and explicit-reference convergence.
Initially retain the old CFL safety multiplier. The bound covers a frozen local
Jacobian; it does not prove nonlinear positivity, RKL stability, or AMR stability.

**10. Separate STS safety from advective CFL.**
Only after piece 9 passes, change both initial and refreshed STS budget selection;
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

**12. Repair refinement/boundary correctness before compact communication.**
Two separate fixes: refill physical corners after prolongation (`e1b2ecd00`),
then synchronize shared LF face fluxes near refinement (`a6eef6406`). Reuse the
existing buffers and A/magnetic-moment-aware prolongation; retain magnetic BCs.
Gate: uniform anisotropic oblique-field SMR/outflow state, the original nonuniform
fixture, periodic E and integrated magnetic-moment conservation, both prolongation
modes, 2D/3D, explicit/STS, CPU/HIP and multiple ranks. Preserve AMR churn/restart
and repair-accounting checks. Do not equate positive pressures with conservation.

**13. Apply only measured, equivalent optimizations.**
Try separately: redundant LF flux-clear removal (`15e7ff619`); fused primitive/
temperature refresh with valid frozen-B reuse (`d02bb714d`); upwind-only HLLE
logarithm (`31b739890`); compact IEN/IAN messages (`6539bfbbc`); then any proven
redundant timestep refresh. Keep full shearing exchange and required magnetic
boundary work. Gate: same-configuration full-precision state/restart equality,
including ghosts/coarse data and scalars, plus paired GPU timing. Keep only a
beneficial or independently justified simplification; no configurable fallback
framework. Handle diagnostic reduction-order differences explicitly.

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

**15. Qualify nonlinear LF time integration before selecting a chunk cap.**
If the corrected passive target still exhibits the retained depleted-pressure
failure, use `4af7c9a55` as a small implementation reference. Hold the outer step,
forcing and collision cadence fixed; change only LF substeps within each sweep.
Gate: actual failing-state reproduction and timestep refinement, worst-cell
pressure differences as well as norms, energy and smooth temporal order. Recorded
positive chunked runs still differed locally by about 14%; cap 4 is not an
accuracy-qualified default. Do not merge chunking with a kernel optimization.

## Deferred work

- LF inflow/user boundaries: implement only when needed, using `738810c88` and
  independent uniform-state/callback tests across all newly supported paths.
- Directional fusion: consider after profiling the corrected diagnostic-free
  workload. Retain only a measured gain with one readable face-arithmetic
  implementation; do not copy the large duplicate fused kernel by default.
- Half-sweep merging, entropy-variable fluxes, A-to-Q reformulation, new limiter
  models, symmetric collision splitting, and broad float portability are outside
  this reconstruction. No speculative runtime switches or analysis framework.
- Leave unrelated cleanup and the large evidence archive on the reference branch.
  Bring across only the fixture, checker and short derivation needed for a change.

## Checks and measurements for every piece

Use `tst/test_suite/cgl/`, `tst/scripts/cgl/`, `test_suite.testutils`, and existing
inputs/pgen fixtures. Inspect their working-directory and backend setup before
running; the generic runner's GPU default is CUDA, not Frontier HIP. Recover small
later fixtures from the pinned reference when needed. Keep baseline failures and
known limits explicit; never replace physical accuracy checks with smoke completion.

1. Reproduce the defect or independently establish the expected result first.
   Run the smallest relevant CPU tests and one affected integration case. For
   GPU arithmetic, driver, MPI or boundary changes, also run the affected actual
   HIP/MPI cases before calling the piece qualified. Missing hardware is a
   recorded validation gap, not a pass.
2. Measure baseline and candidate with the same complete build configuration,
   including Kokkos, and record source/input/executable hashes. Match arithmetic,
   diagnostics, strictness, reconstruction, CFL, limits, forcing, block shape and
   rank/node layout unless the patch explicitly changes one of them.
3. Use one large uniform aligned active case and one genuinely 3D active case
   with independently verified nonzero LF transport. Add a refinement timing
   only for relevant pieces. Avoid using planar zero-transport forcing as the
   sole performance case. Reuse checkpoints only while physics/encoding match;
   otherwise initialize explicitly matched primitive states.
4. Use existing cycle/time/elapsed stdout. Time a fixed completed-cycle interval
   after warmup, with no output writers inside it and fixed diagnostic frequency.
   Run a warmup per binary and three alternating sequential baseline/candidate
   pairs on the same allocation. Report individual values, median and spread;
   lengthen noisy intervals instead of selecting favorable repeats.
5. Record zone-cycles/s/node, LF RHS evaluations/cycle, elapsed/physical-time
   advanced, and amortized whole-step elapsed/RHS. The latter is not pure LF
   kernel time. Obtain kernel profiles separately because profiling fences perturb
   timing. Existing `lf_nstage` counts zone-stage evaluations: on a fixed mesh,
   divide its increment by active cells to obtain RHS evaluations. Use a separate
   count/endpoint run when output would affect timing.
6. Compare physics at the same physical time. For optimizations require exact
   primary-state equality on the same backend/toolchain/layout; for deliberate
   operator changes require the independent accuracy/conservation gate. Report
   any diagnostic-only rounding separately. Assess min/max and spatial outliers,
   not just RMS errors. Stability, positivity and accuracy are separate tests.
7. Record one compact result row per commit: requirement, source scope, commands,
   passed/failed/skipped checks, changed expectations, timing/counts, and retained
   limitation. Stop on an unexplained failure or regression; investigate before
   adding the next change. A necessary physics correction may cost time, but its
   cost must be attributed and its implementation reviewed rather than hidden.

After the active correctness pieces pass, run the broader relevant CPU/HIP/MPI
suite once, including resolved wave/decay and refinement/restart cases. Freeze
that reference before optimizations and establish a separate reference when the
passive model is added. No speedup target overrides the physical acceptance gates.
Include a resolved active turbulence check before freezing the active reference
or starting piece 13. After piece 14, qualify the passive target separately using
piece 15; the retained depleted-pressure failure was passive. Small smooth tests
alone do not release either model for the intended turbulent run.

The immediate next action is piece 0: qualify and measure the unmodified
pre-WO1 and post-WO1 sources. None of the retained historical passes certify the
new patch composition in advance.
