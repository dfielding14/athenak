# Temporary CGL AMR Implementation Guide

This is a temporary engineering guide for implementing CGL-aware AMR support.
Do not add this page to the Sphinx toctree or public documentation index. It is
intended to be removed once the implementation, tests, and permanent docs are
complete.

## Purpose

The current CGL AMR path is only partially supported. CGL variables are carried
through the generic MHD AMR arrays and buffers, but the primitive prolongation
path is not CGL-aware, and generic conserved interpolation does not preserve CGL
thermodynamic admissibility.

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

Current limitations:

- `prolong_primitives=true` is rejected for CGL Landau-fluid AMR, but pure CGL
  can still enter the primitive path.
- The primitive boundary converters use ideal-MHD C2P/P2C and do not handle
  `IPP` or `IAN` correctly.
- Generic conserved refinement/restriction treats `IAN` like an ordinary scalar.
  That is not admissibility preserving because CGL `IAN` is a nonlinear
  thermodynamic variable.
- During a CGL Landau-fluid STS sweep, `IAN` temporarily stores magnetic moment
  `p_perp / |B|`, not conserved anisotropy.

## Design Principle

Keep generic AMR numerics generic:

- use existing face-centered magnetic-field restriction/prolongation and CT repair;
- use existing conservative transfer for density, momentum, total energy, and scalars;
- add CGL-specific wrappers that decode and reconstruct the thermodynamic state.

Do not interpolate CGL `IAN` as a passive scalar. Decode the CGL thermodynamics,
interpolate/project the pressure state, and recompute `IAN`.

## Phase 0: Guard Unsafe Modes

Before adding new functionality, make failure modes explicit.

Required behavior:

- Keep rejecting CGL LF AMR with `prolong_primitives=true` until a
  magnetic-moment-aware primitive path exists.
- Add a fatal guard or clearly opt-in control for pure CGL
  `prolong_primitives=true` until the CGL-aware converter is implemented.

Suggested message:

```text
CGL AMR primitive prolongation requires the CGL-aware primitive boundary path.
Use conserved prolongation or enable the experimental CGL AMR primitive mode.
```

## Phase 1: CGL-Aware Primitive Boundary Prolongation

This phase fixes coarse/fine ghost-boundary primitive prolongation. It does not
yet make active block creation fully CGL-aware.

### ConsToPrimCoarseBndry

Modify the MHD overload in `src/bvals/prolong_prims.cpp`.

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

Modify the MHD overload in `src/bvals/prolong_prims.cpp`.

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
- LF primitive prolongation remains disabled until the `IAN` magnetic-moment
  representation is handled explicitly.

## Phase 2: CGL AMR Projection Helper

Add a reusable CGL AMR projection helper. The helper should operate on a single
cell state and should be callable from boundary prolongation, restriction,
refinement, and post-regrid repair paths.

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

For true CGL prolongation, use a slope-scaling approach before final clipping.

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
- Keep LF primitive prolongation disabled until this representation split has
  explicit tests.

## Suggested Runtime Controls

Do not silently change the meaning of the existing generic
`prolong_primitives` flag. Add CGL-specific controls while this is experimental.

Suggested controls:

```ini
<mhd>
cgl_amr_prolongation = conserved        # conserved | primitive_boundary | primitive_all
cgl_amr_admissibility = strict          # strict | floor_project | limiter_project
```

Suggested behavior:

- `conserved`: current validated mode.
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
- Negative test rejecting LF primitive prolongation.

### Stage 1: Unit-Level Projection Tests

Add focused tests for:

- constant-state round trip;
- isotropic CGL round trip;
- low-`B` fallback near `bfloor`;
- pressure floors;
- mirror and firehose hard-wall bounds;
- slope scaling preserving sibling averages when no final projection is needed.

### Stage 2: Pure CGL Primitive Boundary AMR

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

Only after pure CGL is stable:

- add LF AMR tests in conserved mode with the new projection helpers available;
- add LF primitive-boundary mode only if magnetic-moment representation is
  explicitly handled;
- require strict LF counters to remain clean.

### Stage 5: 3D/GPU And Physics Tests

Add:

- small 3D CGL AMR smoke;
- GPU version if CI resources permit;
- CGL wave crossing a coarse/fine interface;
- passive CGL AMR case;
- limiter stress case near mirror/firehose thresholds.

## Rollout Order

1. Guard unsafe primitive-prolongation modes.
2. Implement CGL-aware primitive boundary converters.
3. Add pure-CGL primitive-boundary tests.
4. Implement reusable CGL AMR projection helper.
5. Add CGL-specific block refine/restrict/regrid transfer.
6. Add repair counters and docs for runtime controls.
7. Extend tests to LF representation, 3D, GPU, restart, and physics convergence.
8. Remove this temporary guide after the implementation and permanent docs land.

## Non-Goals For The First Patch

- Do not make LF primitive prolongation work in the first patch.
- Do not change default AMR behavior.
- Do not rewrite generic AMR interpolation kernels.
- Do not claim exact conservation of CGL `IAN` across AMR transfers.
- Do not add this temporary page to the public docs index.
