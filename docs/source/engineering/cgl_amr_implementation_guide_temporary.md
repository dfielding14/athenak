# Temporary CGL AMR Implementation Guide

This is a temporary engineering guide for hardening CGL-aware AMR support.
Do not add this page to the Sphinx toctree or public documentation index. It is
intended to be removed once the implementation, tests, and permanent docs are
complete.

## Purpose

The initial CGL AMR primitive-transfer path now exists. CGL variables are
carried through the MHD AMR arrays and buffers, primitive fine/coarse boundary
prolongation is CGL-aware, and active MeshBlock creation/deletion has CGL
wrappers that rebuild the thermodynamic state instead of treating `IAN` as a
passive scalar.

The stronger active-block transfer path now exists: active block refinement
interpolates thermodynamic variables such as `U` and `Delta`, scale-limits
children before creating inadmissible states, and centralizes final
projection/repair accounting in a reusable helper. This temporary guide remains
until the validation matrix below is complete and permanent documentation is
updated.

This guide describes a staged implementation plan for:

1. CGL-aware primitive boundary prolongation.
2. Admissibility-preserving CGL interpolation/projection for AMR transfers.
3. Tests and rollout controls.

The goal is to preserve the existing generic AMR machinery where possible while
adding CGL-specific thermodynamic decoding, projection, and reconstruction at
the physics boundaries.

## Current State

CGL uses a sixth MHD slot:

- `IEN` / `IPR`: total energy in conserved arrays, parallel pressure in primitive arrays.
- `IAN` / `IPP`: conserved CGL anisotropy in conserved arrays, perpendicular pressure in primitive arrays.

Important source paths:

- `src/athena.hpp`: defines `IAN`, `IPP`, and the CGL one-dimensional state fields.
- `src/eos/ideal_c2p_mhd.hpp`: defines `CGLConservedAnisotropy`,
  `SingleC2P_CGLMHD`, `SingleP2C_CGLMHD`, and magnetic-moment conversion helpers.
- `src/eos/cgl_mhd.cpp`: full-array CGL primitive recovery, pressure floors,
  hard-wall projection, collisions, and LF slot conversion.
- `src/mhd/mhd_tasks.cpp`: AMR task flow, including conserved or primitive prolongation.
- `src/bvals/prolong_prims.cpp`: current primitive boundary conversion path.
- `src/mesh/mesh_refinement.cpp`: actual block refinement and restriction.

Current implemented support:

- `prolong_primitives=true` is supported for pure CGL and CGL LF/STS AMR.
- The primitive boundary converters in `src/bvals/prolong_prims.cpp` handle
  CGL `IPP` and `IAN` and accept an explicit magnetic-moment flag for LF STS
  stages.
- `src/mesh/mesh_refinement.cpp` has CGL-specific block refine/restrict
  wrappers that transfer `U`/`Delta`, apply floors/hard-wall projection, and
  rebuild conserved `IAN`.
- `src/eos/cgl_amr_projection.hpp` centralizes single-cell AMR projection,
  anisotropy vs LF magnetic-moment encoding, and per-cell repair masks.
- `divb_amr` user histories expose CGL AMR repair counters for smooth/stress
  validation.
- During a CGL Landau-fluid STS sweep, `IAN` temporarily stores magnetic moment
  `p_perp / |B|`; CGL AMR code paths now branch on that representation instead
  of assuming conserved anisotropy.
- CGL LF pressure-work recording remains disabled with AMR primitive
  prolongation because that diagnostic path has not been audited across
  fine/coarse flux communication.
- AMR MPI load-balance buffers are host-staged for GPU portability.

Current limitations:

- Fine/coarse ghost-boundary primitive prolongation is CGL-aware but still uses
  component-wise primitive interpolation followed by projection; the stronger
  `U`/`Delta` slope-scaling path is currently implemented for active block
  creation/regrid refinement.
- The repair counters are cumulative history diagnostics, not a full event-log
  system with per-step reset semantics.
- Runtime controls are still the generic `prolong_primitives` flag plus the
  existing LF pressure-work guard; the CGL-specific controls suggested below
  are not implemented yet.

## Design Principle

Keep generic AMR numerics generic:

- use existing face-centered magnetic-field restriction/prolongation and CT repair;
- use existing conservative transfer for density, momentum, total energy, and scalars;
- add CGL-specific wrappers that decode and reconstruct the thermodynamic state.

Do not interpolate CGL `IAN` as a passive scalar. Decode the CGL thermodynamics,
interpolate/project the pressure state, and recompute `IAN`.

## Phase 0: Guard Unsafe Modes

Status: mostly complete, but the guard has changed from the original broad
rejection.

Implemented behavior:

- Ordinary CGL consumers call representation guards before using `IAN` as
  conserved anisotropy.
- LF STS callbacks require the temporary magnetic-moment representation when
  they operate inside the split LF sweep.
- CGL LF AMR with `prolong_primitives=true` is allowed only through the
  CGL-aware primitive and magnetic-moment paths.
- CGL LF pressure-work recording is rejected with AMR primitive prolongation.

Current pressure-work guard message:

```text
CGL Landau-fluid pressure-work recording is not supported with AMR primitive
prolongation. Set <mhd>/cgl_lf_record_pressure_work = false for LF/STS AMR runs.
```

## Phase 1: CGL-Aware Primitive Boundary Prolongation

Status: implemented for both ordinary CGL anisotropy and the LF STS
magnetic-moment representation.

### ConsToPrimCoarseBndry

The MHD overload in `src/bvals/prolong_prims.cpp` now does this.

Required changes:

1. Capture `eos.is_cgl`.
2. For CGL, load:

   ```cpp
   u.mu = cons(m, IAN, k, j, i);
   ```

3. For CGL, call:

   ```cpp
   SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used,
                    tfloor_used, bfloor_used);
   ```

4. If `eos.hardwall_lim` is enabled and `|B| > bfloor`, apply the same
   hard-wall projection used by `CGLMHD::ConsToPrim`.
5. Store both CGL primitive pressures:

   ```cpp
   prim(m, IPR, k, j, i) = w.e;   // p_parallel
   prim(m, IPP, k, j, i) = w.pp;  // p_perp
   ```

6. Keep scalar conversion unchanged.

### PrimToConsFineBndry

The MHD overload in `src/bvals/prolong_prims.cpp` now does this.

Required changes:

1. Capture full `eos`, not only `gamma`.
2. For CGL, load:

   ```cpp
   w.e  = prim(m, IPR, k, j, i);
   w.pp = prim(m, IPP, k, j, i);
   ```

3. Enforce minimum admissibility before calling CGL P2C:

   ```cpp
   w.d  = fmax(w.d, eos.dfloor);
   w.e  = fmax(w.e, eos.pfloor);
   w.pp = fmax(w.pp, eos.pfloor);
   ```

4. Apply hard-wall projection if configured.
5. Call:

   ```cpp
   SingleP2C_CGLMHD(w, eos.bfloor, u);
   ```

6. Store:

   ```cpp
   cons(m, IAN, k, j, i) = u.mu;
   ```

### Invariants

- Primitive boundary mode must never prolong conserved `IAN` directly.
- Fine-side `IAN` must be reconstructed using the already-prolonged fine
  face-centered magnetic field.
- CGL primitive pressures must be positive before any logarithm is evaluated.
- LF primitive prolongation must pass the explicit magnetic-moment flag during
  LF STS stages. Do not let a generic CGL consumer interpret magnetic moment as
  conserved anisotropy.
- LF pressure-work recording remains disabled for AMR primitive prolongation
  until that diagnostic path is audited.

## Phase 2: CGL AMR Projection Helper

Status: implemented for active block restriction/refinement in
`src/eos/cgl_amr_projection.hpp`.

The reusable CGL AMR projection helper operates on a single cell state and is
callable from boundary prolongation, restriction, refinement, and post-regrid
repair paths.

Recommended thermodynamic variables:

```text
U     = p_perp + 0.5 * p_parallel
Delta = p_perp - p_parallel
```

Here `U` is CGL internal energy density and `Delta` is pressure anisotropy.

Projection algorithm:

```text
rho = max(rho, dfloor)
B_eff = max(|B|, bfloor)
U = E - 0.5 * |m|^2 / rho - 0.5 * |B|^2
U = max(U, 1.5 * pfloor)

Delta_min = 3 * pfloor - 2 * U
Delta_max = U - 1.5 * pfloor

if hardwall firehose:
    Delta_min = max(Delta_min, firehose_threshold * B^2)

if hardwall mirror:
    Delta_max = min(Delta_max, 0.5 * B^2)

if |B| <= bfloor:
    Delta = 0
else:
    Delta = clamp(Delta_candidate, Delta_min, Delta_max)

p_parallel = (2/3) * (U - Delta)
p_perp     = (2/3) * U + (1/3) * Delta
```

Then store:

```text
normal CGL: IAN = CGLConservedAnisotropy(rho, p_parallel, p_perp, B_eff)
LF STS:     IAN = p_perp / B_eff
```

The helper should report whether it changed density, internal energy,
anisotropy, or hard-wall bounds so AMR-specific repair counters can be added.

## Phase 3: Admissibility-Preserving Prolongation

Status: implemented for active MeshBlock creation in
`MeshRefinement::RefineCGLMHDPrimitives`. Fine/coarse ghost-boundary primitive
prolongation is still the older CGL-aware component-wise primitive path.

For true CGL active-block prolongation, use a slope-scaling approach before
final clipping.

Procedure:

1. Build parent primitive thermodynamics from the coarse conserved state.
2. Compute limited slopes for `rho`, velocity, `U`, and `Delta`.
3. Generate child candidate states.
4. Check all children for:
   - finite density and pressures;
   - positive `p_parallel` and `p_perp`;
   - hard-wall firehose and mirror bounds when enabled;
   - valid low-`B` fallback behavior.
5. If any child violates admissibility, scale all child deviations from the
   parent toward zero by a common factor.
6. Apply final projection only if needed.
7. Recompute conserved variables, including `IAN`.

The common slope-scaling step preserves sibling averages when possible. Final
projection may change the thermodynamic state; count those changes as CGL AMR
repairs.

## Phase 4: Restriction And Regrid Transfers

Status: active block restriction/refinement uses the `U`/`Delta` projection
helper. Restriction preserves the generic conservative averages first, averages
child `Delta` as the thermodynamic anisotropy candidate, projects the coarse
state, and rebuilds `IAN` in either anisotropy or magnetic-moment
representation.

Restriction should preserve density, momentum, total energy, and magnetic flux
first. Treat `IAN` as a derived thermodynamic invariant.

Recommended restriction path:

1. Use generic restriction for `IDN`, `IM1`, `IM2`, `IM3`, and `IEN`.
2. Restrict or reconstruct the thermodynamic anisotropy from child states using
   `U` and `Delta`.
3. Project the coarse state into the admissible CGL domain.
4. Recompute `IAN`.

Exact conservation of `IAN` is lower priority than an admissible CGL state.
If desired later, implement a bounded residual redistribution for `IAN`, but
drop residuals that cannot fit without violating floors or hard-wall bounds.

For actual block creation in `MeshRefinement::RefineCC`, do not rely on generic
component-wise interpolation of `IAN`. Add a CGL-specific path that either:

- refines thermodynamic variables and reconstructs conserved state; or
- applies the CGL projection immediately after generic refinement.

The first option is more accurate. The second is smaller but should be marked as
a transitional implementation.

## Phase 5: Landau-Fluid STS Representation

During LF STS, `IAN` stores magnetic moment:

```text
IAN = p_perp / |B|
```

This means ordinary CGL anisotropy projection cannot be blindly applied inside
the LF STS sweep.

Rules:

- Normal hyperbolic AMR transfers should require anisotropy representation.
- LF parabolic-stage AMR transfers must know when `IAN` is magnetic moment.
- Boundary/user BCs used during LF stages must not interpret `IAN` as conserved
  anisotropy.
- LF primitive prolongation is allowed through the explicit magnetic-moment
  path, but pressure-work recording remains disabled for AMR primitive
  prolongation.

## Suggested Runtime Controls

The current implementation still uses the existing generic
`prolong_primitives` flag. The following CGL-specific controls remain a future
cleanup if we want finer rollout modes than the generic flag provides.

Suggested controls:

```ini
<mhd>
cgl_amr_prolongation = conserved        # conserved | primitive_boundary | primitive_all
cgl_amr_admissibility = strict          # strict | floor_project | limiter_project
```

Suggested behavior:

- `conserved`: generic conserved AMR transfer.
- `primitive_boundary`: CGL-aware coarse/fine ghost-boundary primitive prolongation only.
- `primitive_all`: boundary plus block refine/restrict/regrid CGL thermodynamic transfer.
- `strict`: fail on any AMR-created inadmissible state.
- `floor_project`: project floors only.
- `limiter_project`: project floors and selected mirror/firehose hard walls.

Add AMR-specific CGL repair counters separate from LF counters.

## Validation Plan

### Stage 0: Preserve Existing Behavior

Keep current tests passing:

- CGL LF conserved-prolongation AMR smoke.
- CGL LF AMR MPI reproducibility.
- CGL LF turbulence-driving AMR restart test.
- Negative test rejecting LF pressure-work recording with AMR primitive
  prolongation.

### Stage 1: Unit-Level Projection Tests

Focused projection coverage now exists through AMR stress decks and the GPU
test module. Keep coverage for:

- constant-state round trip;
- isotropic CGL round trip;
- low-`B` fallback near `bfloor`;
- pressure floors;
- mirror and firehose hard-wall bounds;
- slope scaling preserving sibling averages when no final projection is needed.

The retained stress decks intentionally exercise low-`B`, firehose, mirror,
and slope-scaling repair counters while smooth decks require zero CGL AMR
repairs.

### Stage 2: Pure CGL Primitive Boundary AMR

Pure-CGL primitive AMR inputs now exist for uniform and current-triggered
refinement. Keep them as smoke tests and extend them as the stronger
projection/slope-scaling helper lands.

Add a pure-CGL `divb_amr` input with:

```ini
<mhd>
eos = cgl
cgl_heat_flux = none

<mesh_refinement>
prolong_primitives = true
```

Acceptance criteria:

- `bad_state == 0`
- finite `abs_anis`
- `max_ndiv < 1.0e-12` for CPU 2D
- no unexpected density or pressure repairs
- MPI one-rank and four-rank histories agree within roundoff

### Stage 3: True Refine/Derefine Tests

Exercise active block creation, deletion, and restart:

- moving refinement pattern;
- refine/derefine through anisotropic pressure structure;
- restart immediately after regrid;
- compare against a uniform-grid or fixed-refinement reference.

Acceptance criteria:

- no nonfinite conserved or primitive state;
- no new `div B` regression;
- mass, momentum, and total energy budgets match existing AMR tolerances;
- CGL AMR repair counters are zero for smooth states.

### Stage 4: LF AMR Extension

Initial LF AMR primitive support exists. The remaining work is to harden it:

- keep LF AMR tests in both conserved and primitive modes;
- exercise magnetic-moment representation through refine/derefine churn;
- require strict LF counters to remain clean.

### Stage 5: 3D/GPU And Physics Tests

Added focused input decks and a GPU test module for:

- small 3D CGL AMR smoke;
- GPU version if CI resources permit;
- passive CGL AMR case;
- limiter stress case near mirror/firehose thresholds.
- CGL wave crossing a coarse/fine interface with a smooth reference comparison.

### Stage 6: Final Production-Readiness Gate

Before calling CGL LF/STS AMR production-ready, run one final validation set
that separates AMR-transfer correctness from ordinary resolution-history
differences in turbulence.

Required final validation cases:

- one meaningful multi-node GPU LF+STS AMR churn run;
- one current-triggered refine/derefine case, not just refine-up;
- restart-through-regrid validation;
- comparison against uniform-grid or fixed-refinement references for at least
  one smooth problem;
- no unexpected CGL AMR repair counters in smooth cases;
- documented tolerances for mass, momentum, energy, `divB`, and LF counters.

Address these items explicitly before declaring the implementation ready:

- which runs exercised true derefinement, not only block creation;
- whether restart files were written immediately after regrid and then used to
  continue the run;
- which smooth reference problem was used, what reference resolution or fixed
  refinement pattern it used, and what convergence/error metric was applied;
- whether any CGL AMR repair counter fired in a smooth case, and if so whether
  the event is physically justified or a blocker;
- the accepted tolerances for mass, each momentum component, total energy,
  `divB`, LF cap/NaN/floor counters, and CGL AMR repair counters;
- any known differences between the AMR and uniform-grid turbulence histories
  that are resolution-history effects rather than AMR-transfer failures.

The multi-node GPU LF+STS churn run should exercise the same code paths we
would use in production: primitive AMR transfer enabled, Landau-fluid STS
active, real block creation and deletion, and enough evolution after regrid to
exercise repeated fine/coarse synchronization. It does not need to be a long
turbulence campaign, but it must be long enough to make one-step regrid-only
success irrelevant.

The current-triggered case should demonstrate both sides of the trigger. The
history needs to show refinement when `|curl B|` is high and derefinement after
the current falls, the threshold changes, or the tagged structure moves away.

The restart-through-regrid test should write a restart immediately after a
regrid event, resume from that restart, and compare the continued history
against an uninterrupted run over the same interval.

Smooth reference comparisons should use a problem where pointwise or normed
errors are meaningful. Turbulence comparisons should not be used as the only
reference because phase divergence can hide a transfer bug or falsely suggest
one.

A useful turbulent-box version is:

1. Use the existing uniform high-resolution run as the reference.
2. Run a matched uniform low-resolution case to the same final time.
3. Run a low-to-high case: start at half resolution, globally refine at an
   early time such as `0.1T` to `0.25T`, then continue to the final time.
4. Optionally repeat with a later refinement time, such as `0.5T`, to expose
   expected resolution-history sensitivity.
5. Run a current-triggered AMR version to the same final time.

Do not require pointwise field agreement in turbulent runs after substantial
evolution; late-time phase differences can be physical. Compare budgets,
admissibility counters, LF counters, `divB`, spectra, PDFs, energy partitions,
anisotropy distributions, and qualitative slice morphology. For smooth
reference problems, use stricter convergence and conservation tolerances.

## Rollout Order

1. Keep the existing unsafe-mode guards current.
2. Keep CGL-aware primitive boundary converters tested.
3. Keep pure-CGL primitive-boundary tests passing.
4. Implement reusable CGL AMR projection helper.
5. Replace the transitional block refine/restrict/regrid transfer with
   `U`/`Delta` slope scaling.
6. Add repair counters and, if useful, docs for runtime controls.
7. Extend tests to LF representation, 3D, GPU, restart, and physics convergence.
8. Pass the final production-readiness gate above.
9. Remove this temporary guide after the implementation and permanent docs land.

## Non-Goals For The First Patch

- Do not change default AMR behavior.
- Do not rewrite generic AMR interpolation kernels.
- Do not claim exact conservation of CGL `IAN` across AMR transfers.
- Do not add this temporary page to the public docs index.
