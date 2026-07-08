# CGL+LF+STS AMR Code Review and Remediation Plan

- Status: review complete; F1 through F5, F7, and F8 implemented
- Review date: 2026-07-08
- Reviewed source: `40fe047f9ffe74795acea3116caad1cb94e6a796`
- Review range: `f2d0a2597..40fe047f9`
- Companion document: `cgl_amr_implementation_guide_temporary.md`

## Purpose

This document records the code review of the CGL, Landau-fluid (LF),
super-time-stepping (STS), and adaptive mesh refinement (AMR) implementation.
It also defines a detailed, deliberately narrow remediation plan for the four
highest-severity findings.

The review phase was read-only. The F1 changes recorded below were implemented
after the findings and remediation plan were established.

## Review Decision

The four highest-severity numerical corrections are implemented. F1 and F3
intentionally share one guarded task reorder, while F2 and F4 are independent
patches. Focused serial validation and MPI-enabled builds pass; multi-rank and
GPU runtime acceptance still require the target execution environment.

The original four highest-severity findings were:

1. Primitive CGL restriction mutates the live hyperbolic state between the
   conservative update and constrained transport (CT).
2. Active primitive refinement does not conserve momentum, total energy, or
   scalar densities, even when no CGL repair is reported.
3. Coarse CGL thermodynamics are encoded with the pre-CT magnetic field and
   later paired with the post-CT field.
4. The final LF STS primitive refresh leaves outer ghost layers stale before
   the next hyperbolic reconstruction.

These are localized defects. They do not require a rewrite of generic AMR,
the LF flux operator, or the CGL equation of state.

## F1 Implementation Update

F1 was implemented on 2026-07-08 with a narrowly gated task dependency change:

- `src/mhd/mhd_tasks.cpp` defers the existing CGL primitive `RestrictU` and U
  exchange chain until after `CornerE`, CT, and `RecvB_OA`;
- all single-level, non-CGL, and conserved-prolongation configurations retain
  their original task order;
- `src/ion-neutral/ion-neutral_tasks.cpp` applies the same guard to the
  reachable alternate assembler without claiming broader CGL+ion-neutral
  physics support;
- `tst/test_suite/cgl/test_cgl_amr_gpu.py` compares one-cycle primitive and
  conserved AMR modes with refinement suppressed. With no transfer to perform,
  their histories must agree to `1.0e-12`.

The implementation passed serial Debug builds, an MPI-enabled Debug build, the
focused no-regrid comparison, and 2-D, 3-D, and LF+STS refine smoke tests. A
multi-rank runtime check still requires an allocation on the target system.

## Required Invariants

The remediation should be reviewed against these invariants rather than only
against individual test outputs.

### Hyperbolic stage-state invariant

After hyperbolic fluxes have been computed, no AMR helper may mutate the live
fine-grid `u0`, `b0`, `w0`, or `bcc0` state before `CornerE` and CT have
consumed the matching stage state.

### Conservative-transfer invariant

In a no-repair refine or derefine operation, volume-integrated density,
momentum, total energy, magnetic flux, and every scalar density must be
preserved to roundoff. CGL `IAN` is a derived thermodynamic quantity and is
not required to be conserved independently.

### Magnetic-field consistency invariant

Any CGL state encoded into `IAN` must use the same magnetic field that will be
paired with that state by its next consumer. A state must not be encoded using
pre-CT `B` and consumed using post-CT `B`.

### Primitive-ghost invariant

Before any hyperbolic reconstruction, every primitive and cell-centered
magnetic-field ghost cell required by the selected reconstruction method must
have been refreshed from the current conserved and face-centered states.

### Representation invariant

Every use of `IAN` must have an explicit representation:

- normal CGL: conserved anisotropy;
- LF STS: magnetic moment, `p_perp/|B|`;
- private AMR scratch: a named and documented intermediate, never an implicit
  reinterpretation of a live evolution array.

## Findings Map

| ID | Severity | Area | Finding | Disposition |
| --- | --- | --- | --- | --- |
| F1 | Critical | Hyperbolic AMR task order | `RestrictU` mutates live primitives before `CornerE` and CT | Implemented 2026-07-08 |
| F2 | High | Active refinement | Primitive reconstruction is not conservative | Implemented 2026-07-08 |
| F3 | High | Coarse restriction | CGL thermodynamics use pre-CT `B` and are paired with post-CT `B` | Implemented 2026-07-08 by the F1 task reorder |
| F4 | High | LF STS lifecycle | Final refresh leaves reconstruction ghost cells stale | Implemented 2026-07-08 |
| F5 | Medium | Collision lifecycle | Non-LF collision relaxation can invalidate the stored next timestep | Implemented 2026-07-08 |
| F6 | Medium | Physical boundaries | LF fixed-inflow and user boundaries are representation-blind | Implemented 2026-07-08 |
| F7 | Medium | Repair diagnostics | Repair counters omit paths and do not survive restart | Implemented 2026-07-08 |
| F8 | Medium | Failure handling | Nonfinite face fields can remain nonfinite after a reported repair | Implemented 2026-07-08 |
| F9 | Medium | Validation pgen | Multi-block face-field initialization uses incorrect coordinates | Implemented 2026-07-08 |
| F10 | Medium | MPI validation | Quantitative Fourier projections are rank-local | Implemented 2026-07-08 |
| F11 | Low/Medium | Regression coverage | MPI+GPU primitive AMR, restart state, and conservation gates are incomplete | Implemented; worker qualification pending |
| F12 | Low | Style | Changed files do not fully pass the repository C++ style check | Implemented 2026-07-08 |

## Detailed Findings

### F1: Live State Is Mutated Before CornerE and CT

The hyperbolic task order in `src/mhd/mhd_tasks.cpp` is effectively:

```text
Fluxes
  -> RKUpdate
  -> RestrictU
  -> Send/RecvU
  -> CornerE
  -> CT
  -> RestrictB
```

When primitive prolongation and CGL are enabled, `MHD::RestrictU` calls
`ConsToPrim(u0, b0, w0, bcc0, ...)` in place. This has two effects:

- active-zone `w0` and `bcc0` are replaced with values recovered from the
  post-RK conserved state and pre-CT face field;
- CGL floors or hard-wall projection may also modify live `u0`.

`CornerE` then combines those changed primitives with face EMFs and mass fluxes
computed from the earlier stage state. The generic conserved restriction path
does not perform this mutation, so the inconsistency is specific to the CGL
primitive-AMR branch.

Primary locations:

- `src/mhd/mhd_tasks.cpp:82-100` and `:730-763`:
  `MHD::AssembleMHDTasks` and `MHD::RestrictU`;
- `src/eos/cgl_mhd.cpp:169-184`: mutating behavior of `ConsToPrim`;
- `src/mhd/mhd_corner_e.cpp:59-153` and `:206-318`: post-restriction
  consumers of `w0` and `bcc0`.

### F2: Active Primitive Refinement Is Not Conservative

`MeshRefinement::RefineCGLMHDPrimitives` independently prolongs density,
velocity, CGL internal energy `U`, anisotropy `Delta`, and scalar primitives.
It then reconstructs the conserved fields nonlinearly:

```text
momentum       = rho * velocity
total energy   = U + kinetic energy + magnetic energy
scalar density = rho * scalar primitive
```

Zero-mean primitive deviations do not imply zero-mean conserved deviations.
For example,

```text
<rho v> = rho_parent v_parent + <delta_rho delta_v>.
```

Kinetic and magnetic energy introduce additional variance terms. The defect is
therefore present in smooth, admissible states and does not require projection
or a repair event.

A forced-refinement check produced:

| Quantity | Primitive AMR | Conserved AMR comparison |
| --- | ---: | ---: |
| Mesh cells | 576 -> 2304 | 576 -> 2304 |
| Total-energy change | `+9.7538e-05` | `+2.7e-14` |
| Reported CGL repairs | zero | zero |

Primary location:

- `src/mesh/mesh_refinement.cpp:1411-1561`:
  `MeshRefinement::RefineCGLMHDPrimitives`.

### F3: Coarse CGL State Is Encoded Against the Wrong Magnetic Field

The CGL restriction helper currently runs before CT. It restricts the old face
field into `coarse_b0`, averages child `Delta`, and rebuilds coarse `IAN` using
that field. CT and `RestrictB` then update and replace the magnetic field while
the already encoded CGL thermodynamic state remains in `coarse_u0`.

Because both conserved anisotropy and LF magnetic moment depend on magnetic
field magnitude, the final `coarse_u0(IAN)` is generally not the encoding of
the selected coarse thermodynamic state under the final `coarse_b0`.

Primary locations:

- `src/mhd/mhd_tasks.cpp:82-100` and `:730-763`: ordering of `RestrictU`,
  `CT`, and `RestrictB`;
- `src/mesh/mesh_refinement.cpp:1169-1255`:
  `MeshRefinement::RestrictCGLMHDPrimitivesToCons`.

### F4: Final LF Refresh Does Not Cover the Reconstruction Stencil

`MHD::CGLLandauFluidPrimitiveRefresh` refreshes the active domain plus one
ghost cell. That is sufficient for the immediate LF stencil currently used by
intermediate STS stages, but it is not sufficient for the next hyperbolic
reconstruction.

Default PLM reconstructs through `is-1..ie+1` and reads `i-1` and `i+1`, so it
requires valid primitives through `is-2..ie+2`. The current PPM/WENO paths can
reach three ghost cells and FOFC can extend the demand to four. After the final
pre-sweep stage there is no full C2P before the hyperbolic flux task.
Post-sweep collisions also traverse the full ghost allocation.

Primary locations:

- `src/mhd/mhd_sts.cpp:623-650`:
  `MHD::CGLLandauFluidPrimitiveRefresh`;
- `src/mhd/mhd_fluxes.cpp:85-103`: reconstruction call ranges;
- `src/reconstruct/plm.hpp:42-52`: neighbor accesses.

## Remediation Constraints

The first four fixes should remain within the following boundaries:

- do not change generic non-CGL AMR behavior;
- do not change the default value or broad meaning of `prolong_primitives`;
- do not treat `IAN` as a passively prolongated conserved scalar;
- do not change the LF pressure-work support policy;
- do not add a global synchronization when a task dependency or local fence is
  sufficient;
- do not hide conservation errors by loosening regression tolerances;
- do not combine repair-counter cleanup or unrelated pgen corrections with the
  core numerical patches.

## Remediation Plan for F1

### Target behavior

No CGL primitive restriction or EOS recovery may mutate the live stage state
between flux construction and `CornerE`/CT. The existing mutating
`MHD::RestrictU` implementation may run only after those consumers have
finished and the stage magnetic field is complete.

### Implemented approach

F1 and F3 share a task-order seam and should be fixed by one narrowly gated
hyperbolic task-graph branch. When all three conditions are true --
multilevel mesh, `prolong_primitives=true`, and CGL EOS -- defer the existing
full U restriction and exchange chain until CT and orbital-B completion:

```text
RKUpdate
  -> source terms
  -> Send/RecvU_OA
  -> CornerE
  -> electric-field exchange
  -> CT
  -> Send/RecvB_OA
  -> RestrictU
  -> Send/RecvU
  -> Send/RecvU_Shr
  -> RestrictB
  -> Send/RecvB
  -> Send/RecvB_Shr
  -> physical BCs
  -> Prolongate
  -> ConToPrim
```

Implementation steps completed:

1. In `MHD::AssembleMHDTasks`, detect the CGL primitive-AMR hyperbolic path at
   task assembly time.
2. For that path only, make `CornerE` depend on completed orbital-U work rather
   than on the ordinary U boundary-exchange chain.
3. Make `RestrictU` depend on `RecvB_OA`, which follows CT and any orbital
   magnetic remap.
4. Attach the existing U exchange chain after `RestrictU`, then continue with
   the existing B restriction/exchange, boundary, prolongation, and C2P tasks.
5. Leave the current task graph unchanged for non-CGL, non-primitive, and
   single-level configurations.
6. Keep `MHD::RestrictU` and the boundary communication APIs functionally
   unchanged in the first patch. Once deferred, its in-place primitive recovery
   cannot affect `CornerE` and sees the completed stage field.
7. Audit custom task assemblers. Radiation cannot be combined with CGL under
   the current coordinate guards and its MHD ordering is already post-CT. The
   reachable ion-neutral task list received the corresponding guarded reorder.

This task-only change is preferred to adding a second CGL scalar communication
path. It also avoids a new intermediate representation in `coarse_u0`.

### Files changed

- `src/mhd/mhd_tasks.cpp`;
- `src/ion-neutral/ion-neutral_tasks.cpp`;
- `tst/test_suite/cgl/test_cgl_amr_gpu.py`.

### Verification

The checked-in regression runs the same multilevel CGL problem for one cycle
with refinement suppressed, once with primitive prolongation selected and once
with conserved prolongation selected. Because no coarse/fine transfer occurs,
the selection must not affect the solution. The test compares time, timestep,
all global conserved histories, anisotropy, kinetic energy, and magnetic energy
to `1.0e-12`.

The pre-fix executable fails this test with a timestep difference of about
`1.9e-09` and directional magnetic-energy differences up to `1.0e-07`. The
patched executable produces bitwise-identical history files. Additional serial
smoke tests exercised 2-D refinement, 3-D refinement, and LF+STS refinement.
Both serial and MPI-enabled Debug builds pass; a multi-rank runtime check is
still pending an allocation.

### Acceptance criteria

- No invocation of CGL primitive `RestrictU` occurs before `CornerE`, CT, and
  `RecvB_OA` on the guarded path.
- Smooth two- and three-dimensional CT results do not depend on whether the
  CGL AMR restriction helper ran.
- Same-rank and cross-rank U exchanges both carry the post-CT finalized state.
- Non-CGL and non-primitive AMR task paths remain unchanged.

## Remediation Plan for F2

### Implementation update

F2 was implemented on 2026-07-08 by making generic `RefineCC` authoritative
for density, momentum, total energy, and scalar densities. The CGL completion
now prolongs and limits only `Delta`, then rebuilds `IAN` at fixed conserved
state. A separately counted fallback may change primary fields only when the
generic child state is itself inadmissible.

The focused single-refinement regression conserves mass, all momentum
components, and total energy to roundoff with zero CGL repairs. A real passive
scalar with correlated density and concentration slopes also conserves scalar
mass exactly. One-level, mixed-topology, and deep-AMR audits found no change to
the selected thermodynamics across `RepairAMRFC`.

### Target behavior

For an admissible no-repair parent state, child averages must exactly recover
the parent conserved density, momentum, total energy, and every scalar density.
Only the derived CGL `IAN` encoding may differ from a generic component-wise
prolongation.

### Proposed implementation

Refactor `RefineCGLMHDPrimitives` into a conservative baseline plus a CGL
thermodynamic completion:

1. Use the existing generic conserved prolongation for the complete conserved
   array. Its temporary component-wise `IAN` result will be overwritten by the
   CGL completion pass; `IDN`, `IM1`, `IM2`, `IM3`, `IEN`, and all scalar
   densities remain authoritative.
2. Keep the existing divergence-preserving face-field prolongation and CT
   repair. Use the provisional child field for admissibility planning, but
   defer final `IAN` encoding until the repaired child faces are available.
3. Decode the parent CGL state into `U_parent` and `Delta_parent`.
4. Prolong only the thermodynamic anisotropy candidate `Delta` with the common
   sibling slope factor. Do not independently prolong velocity and then
   rebuild momentum; do not convert scalar primitives back to densities.
5. For each child, derive

   ```text
   rho = child conserved density
   v   = child momentum / rho
   U   = child total energy - child kinetic energy - child magnetic energy
   ```

6. Determine the admissible `Delta` interval at fixed `rho`, momentum, total
   energy, and magnetic field. Verify first that the zero-slope state is
   admissible for every sibling, then scale the sibling `Delta` deviations by
   one common factor until every child is admissible. If the zero-slope state
   is not admissible because child `B` changes the interval, project each
   child's `Delta` at fixed `U`; do not change the generic conserved fields.
7. Encode `IAN` from the accepted child `U`, `Delta`, and magnetic field. In
   the normal no-repair case, write only `IAN`; leave the conservatively
   prolonged fields byte-for-byte unchanged.
8. Treat an invalid conservative child, such as subfloor density or
   insufficient internal energy, as a separate fallback. First attempt a
   sibling-coupled conservative slope reduction or bounded redistribution.
   Only if no admissible conservative solution exists may the existing
   projection alter density or energy, and that change must be explicitly
   counted as a conservation-changing repair, distinct from ordinary slope
   scaling, and exposed in conservation diagnostics. A strict/debug mode
   should fail with the offending state rather than silently take this path.
9. Preserve the accepted child `Delta` until AMR face exchange and
   `RepairAMRFC` have completed. Recompute `U` using the unchanged conservative
   state and the repaired final child field, then perform the final `IAN`
   encoding. This prevents a second pre-/post-field mismatch in newly created
   blocks.

The key separation is:

```text
generic AMR owns conserved interpolation
CGL AMR owns Delta admissibility and IAN reconstruction
```

### Avoided shortcuts

The following are not acceptable fixes:

- subtracting the measured energy error from one arbitrary child;
- accepting a looser energy tolerance;
- conserving only density while continuing to rebuild momentum and energy
  from independently prolonged primitives;
- skipping projection and allowing invalid child pressures;
- conserving scalar primitives instead of scalar densities.

### Files expected to change

- `src/mesh/mesh_refinement.cpp`;
- `MeshRefinement::AdaptiveMeshRefinement`, where final face repair and the
  last primitive recovery are ordered;
- possibly `src/eos/cgl_amr_projection.hpp` for a fixed-conserved-state
  `Delta` projection helper;
- the MHD allocation only if a small `Delta` scratch view is needed across
  face exchange and `RepairAMRFC`;
- `tst/test_suite/cgl/test_cgl_amr_gpu.py`;
- the current-churn and passive-scalar CGL AMR inputs if an explicit single
  regrid checkpoint is needed.

### Verification matrix

Add forced single-refinement tests in one, two, and three dimensions with:

- simultaneous nonzero density and velocity slopes;
- nonuniform magnetic fields;
- nonzero parallel and perpendicular pressure slopes;
- at least one passive scalar with a nonuniform concentration;
- projection disabled by construction in the smooth case;
- a second case that requires common `Delta` slope scaling;
- a stress case that requires the documented repair fallback.

At the output immediately before and after regridding, assert:

- density conservation to roundoff;
- all momentum components to roundoff;
- total energy to roundoff in the no-repair and slope-scaling cases;
- every scalar density to roundoff;
- magnetic flux and divergence within the existing CT tolerances;
- zero repair counters for the smooth case;
- the precise expected repair category for the fallback case.

Repeat the current-triggered churn test with at least one refine/derefine cycle.
The previously observed `9.7538e-05` refinement energy jump must disappear; it
must not merely move to derefinement or restart.

### Acceptance criteria

- No-repair regrids conserve all generic conserved variables to roundoff.
- Common `Delta` slope scaling does not change those conserved totals.
- Passive scalar mass is conserved independently of density variation.
- Child pressures remain finite, positive, and inside enabled hard-wall bounds.
- Final `IAN` remains consistent with the field after `RepairAMRFC`.
- An infeasible primary state cannot be reported as an ordinary slope repair.
- The implementation continues to use the existing generic face-field AMR
  machinery.

## Remediation Plan for F3

### Implementation update

F3 is implemented by the same post-CT task dependency committed for F1. A
separate numerical rewrite is neither necessary nor desirable: serial audits
covering one-level, mixed-topology, and deep refinement found no additional
thermodynamic change across `RepairAMRFC`, while a late active-state rewrite
would occur after boundary communication.

The F3-specific regression decodes every restricted `coarse_u0(IAN)` using the
current `coarse_b0` and compares its `Delta` with the final child average. The
pre-fix task order produces a maximum error of `1.0016e-03`; the corrected
post-CT order reduces it to `4.4e-15` in the focused two-dimensional run.

### Target behavior

Coarse `IAN` must be finalized using the completed stage magnetic field. No
consumer may observe a coarse CGL state whose thermodynamic encoding and
magnetic field come from different stage times.

### Proposed implementation

Use the same guarded task-graph change specified for F1. Deferring the existing
restriction and full U exchange is sufficient and avoids any new
representation or communication protocol:

1. Complete `CornerE`, electric-field communication, CT, and `RecvB_OA`.
2. Run the existing `MHD::RestrictU`. Its internal
   `RestrictFC(b0, coarse_b0)` now sees completed stage faces.
3. In `RestrictCGLMHDPrimitivesToCons`, decode final child states, average
   `Delta`, and encode coarse `IAN` using that final restricted field.
4. Run the existing full U exchange. Fine senders pack finalized `coarse_u0`
   for coarser neighbors, and coarse receivers obtain the completed CGL state.
5. Run the ordinary `RestrictB` and B exchange. This initially recomputes the
   same restricted field already formed inside `RestrictU`; keep the redundant
   operation in the first surgical patch rather than combining numerical and
   cleanup changes.
6. Continue through physical boundaries and primitive prolongation only after
   both completed U and B exchanges.

The ownership audit shows why finalizing only a sender's `coarse_u0(IAN)` after
the current early `SendU` cannot work:

- `PackAndSendCC` copies buffer values immediately; same-rank paths also copy
  through receive buffers and are not aliases;
- the coarse receiver has aggregate U and B but does not own the fine children
  needed to recover the final average `Delta`;
- `Prolongate` cannot reconstruct that missing thermodynamic information.

Therefore, if early U communication were retained, a second fine-to-coarse CGL
scalar exchange would be mandatory. Deferring the existing full exchange is
smaller and reuses the established buffer/request lifecycle.

The LF parabolic graph already orders `STSUpdateB` before `RestrictU`. Confirm
that property with a task-trace test, but do not change that graph as part of
this fix unless the test exposes a separate ordering problem.

### Files expected to change

- `src/mhd/mhd_tasks.cpp` for the guarded task dependency;
- optionally an assertion or explanatory comment in `MHD::RestrictU`;
- no functional boundary-communication or mesh-refinement change should be
  required for the primary fix;
- a focused AMR restriction test and the MPI CGL AMR regression suite.

### Verification

Add a manufactured restriction test in which CT changes both the magnitude and
direction of child magnetic fields:

1. choose child states with known final `U` and `Delta`;
2. apply a nonzero CT magnetic update;
3. restrict and finalize the coarse CGL state;
4. decode the final coarse `IAN` using final coarse `B`;
5. require the decoded `Delta` to equal the selected restricted/projected
   target, not the value implied by the old magnetic field.

Run this test:

- in two and three dimensions;
- with ordinary anisotropy representation;
- inside LF STS with magnetic-moment representation where that task path is
  supported;
- on one rank and across a fine/coarse MPI boundary;
- with projection inactive and with one controlled projected state.

### Acceptance criteria

- Final coarse `IAN` decodes correctly under final coarse `B`.
- The result is independent of MPI ownership and block decomposition.
- U buffer packing occurs after post-CT CGL restriction on both same-rank and
  cross-rank fine/coarse boundaries.
- No pending or scratch representation is introduced.
- The F1 live-state purity acceptance test remains satisfied.

## Remediation Plan for F4

### Implementation update

F4 was implemented on 2026-07-08 by retaining the active-plus-one refresh for
intermediate LF stages and refreshing the full allocated primitive and
cell-centered magnetic-field extent at the final stage of both split sweeps.
No task dependency, communication path, or intermediate-stage stencil changed.

The regression compares the same periodic 24-by-24 problem as one MeshBlock
and as nine MeshBlocks. Before the fix, the one-cycle total energies differed
by `6.3185e-11`; after the fix, the difference is `5.3291e-15`, and all tested
conserved quantities agree within `5.0e-13`. Conserved- and
primitive-prolongation LF AMR smoke runs remain admissible and repair-free.

### Target behavior

The last LF STS stage of each split sweep must leave `w0` and `bcc0` valid over
the full ghost region needed by the next consumer. Intermediate stages should
retain the smallest safe refresh extent to avoid an unnecessary performance
regression.

### Proposed implementation

1. Keep the current active-plus-one refresh extent for intermediate LF STS
   stages, provided the LF stencil audit confirms that extent is sufficient.
   The conserved communication immediately before refresh already fills the
   full ghost allocation; the purpose here is to bring the corresponding
   primitive ghosts up to date.
2. On `stage == pdrive->sts.nstages`, refresh cell-centered magnetic fields and
   CGL primitives across the full allocated ghost region:

   ```text
   x1: 0 .. nx1 + 2*ng - 1
   x2: 0 .. nx2 + 2*ng - 1 when multidimensional, otherwise 0 .. 0
   x3: 0 .. nx3 + 2*ng - 1 when three-dimensional, otherwise 0 .. 0
   ```

3. Run `RefreshCellCenteredBFromFace` over exactly the same range before
   `CGLRefreshPrimFromMagneticMoment`.
4. Perform the full refresh before
   `EndCGLLandauFluidSTSSweep` changes the `IAN` representation. The primitive
   pressures remain valid after the slot conversion.
5. Ensure that post-sweep collisions see the same fully refreshed `w0/bcc0`.
6. Apply the rule to the final stage of both pre and post split sweeps. The pre
   sweep feeds the hyperbolic integrator; the post sweep feeds output,
   regridding, restart, and the next cycle.
7. Add a short comment at the bounds selection explaining why the final stage
   differs from intermediate stages. Avoid reconstruction-method-specific
   hard-coded halo widths when a full final refresh is simpler and safer.

### Files expected to change

- `src/mhd/mhd_sts.cpp`;
- `tst/test_suite/cgl/test_cgl_amr_gpu.py`;
- `tst/test_suite/cgl/test_cgl_landau_fluid_mpicpu.py` for decomposition
  coverage.

### Verification matrix

Use the same global periodic LF problem with different MeshBlock sizes so that
the number of internal boundaries changes while the physical resolution does
not. Compare after one complete pre-sweep plus hyperbolic stage and after a
full split cycle.

Required cases:

1. **Stage-extent sentinel:** use `nghost=4` and RKL2 with at least three
   stages. Poison outer `w0` halos. Intermediate stages may leave layers beyond
   the LF one-cell requirement untouched; the final stage of both sweeps must
   refresh through the outermost ghost.
2. **Default production path:** run PLM with `nghost=2` using the same periodic
   global problem on one large MeshBlock and several smaller MeshBlocks.
   Compare interface-adjacent active fields after one cycle, not only history
   totals.
3. **Maximum stencil:** run WENO-Z or PPM with FOFC and `nghost=4`. Compare all
   active fields, emphasizing cells within four zones of block boundaries.
4. **Collision path:** enable nontrivial collision relaxation and compare a
   fixed multi-block layout on one rank and multiple MPI ranks after both
   final sweeps.
5. **AMR boundary path:** use a fixed two-level topology with primitive AMR,
   distributed first on one rank and then across ranks. Exercise both
   same-level and coarse/fine face-B-to-`bcc0` refreshes.
6. **GPU path:** repeat a short multi-block case on device to exercise the
   full-extent refresh kernels and their ordering.

Use a smooth state with nonuniform CGL pressure and a demonstrably nonzero LF
update. Require zero floors, hard-wall repairs, and nonfinite events. In
addition to field norms, retain the sentinel or a checksum proving that the
outer primitive layer changed after final conserved ghost communication. A
smooth aggregate history comparison alone may not expose a stale-halo defect.

### Acceptance criteria

- Hyperbolic results are invariant, within the scheme's normal roundoff, under
  MeshBlock decomposition.
- Every reconstruction-required ghost primitive is refreshed after the final
  LF stage.
- Collision-enabled runs do not consume stale outer ghosts.
- Intermediate-stage LF performance is not materially changed.
- No additional global fence is introduced unless required by a measured data
  dependency.

## Patch and Review Sequence

Keep the work reviewable by separating numerical behavior from diagnostics and
test infrastructure.

1. **Reproduction tests:** add focused failing tests for F1 through F4 without
   changing production behavior.
2. **Task-state and field-consistency patch:** fix F1 and F3 together or as two
   consecutive commits with a shared final design. The combined result must be
   tested before merge.
3. **Conservative active-refinement patch:** implement F2 without changing
   generic AMR operators.
4. **Final LF refresh patch:** implement F4 locally in the STS refresh task.
5. **Validation expansion:** run serial CPU, MPI CPU, GPU, and MPI+GPU tests,
   including refine/derefine and restart-through-regrid.
6. **Secondary findings:** address F5 through F12 in separately scoped patches.

Each production patch should include its regression and should pass the
repository style checks before the next numerical change is added.

## Required Final Validation for F1-F4

The combined change should not be accepted until all of the following pass:

- clean debug and release builds;
- existing CGL and LF CPU tests;
- existing CGL MPI tests;
- checked-in GPU CGL AMR tests;
- a real MPI+GPU primitive-AMR run;
- static coarse/fine-interface LF comparison;
- forced single refinement with conservation assertions;
- forced refine/derefine churn with conservation assertions at each topology
  change;
- restart immediately before and immediately after a regrid;
- two- and three-dimensional CT task-order tests;
- repair counters checked for both zero-repair and controlled-repair cases;
- `git diff --check`, C++ style checks, and Python lint.

The validation report should record topology transitions and conservation
jumps at the transition itself. Comparing only initial and final totals can
allow a refinement error to be partially hidden by a later derefinement error.

## Secondary Findings Backlog

The following items remain real review findings but should not be mixed into
the four core numerical patches:

### F5: timestep refresh after ordinary collisions

F5 was implemented on 2026-07-08 without moving the existing collision task.
After ordinary non-LF relaxation, `MHD::CGLCollisions` now calls the shared
state-based MHD timestep estimator. The later mesh-level timestep selection
continues to own the MPI minimum, growth limit, and final-time clipping.

The focused regression uses a uniform state in which isotropization increases
the CGL fast speed. Before the fix, the stored next timestep remained
`1.4940357616679920e-02`; after the fix it matches the analytic post-collision
value `1.3574280209859372e-02` to roundoff. The existing collision finalizer
also verifies that relaxation advanced exactly one physical timestep.

### F6: LF-aware physical and user boundaries

F6 was implemented on 2026-07-08 as a configuration-time guard. Both STS and
explicit LF split integration now reject fixed-inflow and user boundaries on
any active mesh face. During LF stages `IAN` temporarily stores magnetic
moment, while fixed-inflow data has the ordinary conserved-state contract and
user callbacks receive no active-representation argument. There is therefore
no truthful opt-in until one of those boundary contracts is extended.

The guard deliberately permits periodic/shear exchange, reflection, outflow,
diode, and vacuum boundaries because those paths copy, sign-adjust unrelated
components, or zero `IAN` without interpreting it as anisotropy. Focused tests
cover STS inflow rejection, explicit user-boundary rejection, and one-cycle
STS/explicit smokes for reflection, outflow, and diode boundaries.

### F7: trustworthy repair accounting

F7 was implemented on 2026-07-08. The fixed-width cumulative counters now
cover repair-bearing CGL projections that participate in coarse-boundary
packing, derefinement, active-refinement coarse stencils and fine children,
and fine/coarse boundary prolongation. Rectangular scratch cells outside the
actual axial interpolation stencil are not counted, and restriction scratch
is counted only for a real coarse-boundary or derefinement transfer.

The counters are cumulative projection-event counts, not unique physical-cell
counts. A cell projected on multiple stages or boundary exchanges therefore
contributes once per participating invocation. Boundary kernels accumulate
into a persistent device buffer and synchronize it only for history or restart
output. The shared projection path also retains nonfinite and density repair
bits from the original conserved density and uses the repaired density for
passive-scalar conversion.

Restart files carry a versioned array of 12 exact `uint64_t` values. Shared
MPI checkpoints store the global sum and restore that baseline on rank zero
only, so later history reductions do not multiply it by the rank count;
rank-local checkpoints preserve each rank's local values. Files without the
version marker retain the legacy layout and initialize the counters to zero.

Focused tests cover retained conserved-density repair bits, zero accounting
when no transfer occurs, active-refinement accounting, static coarse/fine
boundary accounting, and exact restoration of a nonzero 12-counter vector.

### F8: nonfinite face-field handling

F8 was implemented on 2026-07-08 with an explicit fail-fast contract. The
shared CGL AMR projection helper now aborts before projection or state writes
when a face-derived magnetic component, squared magnitude, or magnitude is
nonfinite. It no longer replaces a local copy of the field with zero or
reports that substitution as a thermodynamic repair while the evolved face
field remains invalid.

A focused CPU death test covers NaN and infinite components as well as finite
components whose squared magnitude overflows. Its finite controls verify that
ordinary projection remains unchanged and that nonfinite thermodynamic inputs
continue through the separately reported repair path.

### F9: quantitative pgen field initialization

F9 was implemented on 2026-07-08. The transverse face-field initializers in
the quantitative CGL LF pgen now compute cell-center coordinates from each
MeshBlock's physical bounds, matching the already-correct primitive
initializer. A zero-cycle fast-eigenmode regression uses four MeshBlocks and
checks the supplied eigenvector; before the fix its measured `By` Fourier
amplitude was effectively zero and failed with relative error one.

### F10: quantitative pgen MPI projections

F10 was implemented on 2026-07-08. MPI-capable primitive and cell-field
projections now all-reduce the raw weighted mean numerator before dividing by
the global domain length. They then calculate local sine/cosine numerators
using that global mean, all-reduce both harmonics together, and normalize only
after the reduction. Every rank therefore receives the same complete
projection before running the quantitative checks.

The MPI regression runs a zero-cycle pure-CGL oblique wave on four 32-cell
MeshBlocks with both one and four ranks. Its constant backgrounds detect a
rank-local mean while its transverse velocity perturbation detects rank-local
Fourier amplitudes.

### F11: regression completeness

F11 now includes full-field restart comparisons, conservation checks at every
history sample through refine/derefine, repair assertions for pure CGL, LF,
and passive-scalar cases, and a clean 24-cubed three-dimensional LF churn. The
restart regression sorts MeshBlocks by logical location and requires every
primitive and cell-centered magnetic field to replay bit-for-bit from both a
checkpoint one cycle before derefinement and the terminal checkpoint directly
after derefinement.

The opt-in MPI+GPU module runs the restart/conservation case on four ranks and
the three-dimensional LF churn on 16 ranks, with one GPU per rank. It is
skipped unless `ATHENAK_RUN_MPI_GPU=1` and requires an allocation-specific
launcher in `ATHENAK_MPI_GPU_LAUNCHER`, so ordinary CPU and single-GPU CI do
not pretend to provide multi-GPU evidence. Target-worker qualification remains
required before F11 is considered accepted.

The first full refine/derefine conservation gate exposed a separate stale-halo
defect after `RepairAMRFC`: the repair ran after boundary exchange, so the next
step could reconstruct from neighboring ghost fields and primitives created
from the pre-repair face field. AMR now reruns the established boundary and
primitive initialization after face repair. The focused churn regression
checks mass, all momentum components, and total energy at every history row;
the previous first-coarse-step errors were as large as `8.7e-7`.

### F12: style

F12 was completed on 2026-07-08 without a broad formatting diff. All C++ files
touched by F5--F11 pass the repository `cpplint` filter and custom whitespace,
brace, pragma, and permission checks. The Python files touched by those
patches pass the configured 90-column `flake8` check; the final style-only
changes were confined to pre-existing violations in a file already modified
by F5, F6, and F9.

## Positive Review Results

No concrete defect was found in the following reviewed areas:

- explicit CGL slot-representation guards;
- LF begin/end representation conversions;
- weighted LF `IEN` and `IAN` flux/reflux handling;
- cycle-end AMR regrid placement;
- post-LF collision followed by timestep refresh;
- GPU value captures in the timestep kernels;
- load-balance host-buffer ordering;
- `U`/`Delta` inversion and pressure-floor interval construction;
- mirror and firehose bound formulas;
- one-, two-, and three-dimensional parent/child index mapping.

These paths should still be retained in regression coverage because the F1-F4
changes touch adjacent lifecycle and AMR code.

## Review Evidence

The review included:

- a clean serial debug build of the reviewed source;
- focused current-triggered refine/derefine runs;
- a forced single-refinement conservation comparison;
- smooth LF oblique-wave comparisons across block layouts;
- inspection of the existing Frontier MPI+GPU AMR campaign artifacts;
- `git diff --check` and changed-file C++ linting;
- source-level audits of task ordering, CGL representation transitions,
  restriction/refinement kernels, reconstruction stencils, boundary paths,
  restart state, and test assertions.

The accepted external campaign remains useful integration evidence, but its
aggregate final-state gates do not detect the F1 task-state mutation or the F2
energy jump at the instant of refinement.
