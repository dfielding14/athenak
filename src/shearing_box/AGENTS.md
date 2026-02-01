# AGENTS.md

## Purpose
This directory implements shearing-box physics: orbital advection in x2, shearing
periodic boundary conditions on x1 faces, source terms for Hydro/MHD, and the
rotating-frame E-field correction used by MHD.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core interfaces
- `shearing_box.hpp`: abstract base classes for orbital advection and shearing box
  boundaries (`OrbitalAdvection`, `ShearingBoxBoundary`) plus CC/FC derived classes,
  shared buffer type (`ShearingBoxBoundaryBuffer`), and task-ID struct.
- `remap_fluxes.hpp`: conservative remap helpers (`DCRemapFlx`, `PLMRemapFlx`).

### Orbital advection
- `orbital_advection.cpp`: base-class constructor, `maxjshift`, MPI communicator.
- `orbital_advection_cc.cpp`: CC pack/send and recv/unpack with conservative remap.
- `orbital_advection_fc.cpp`: FC pack/send and recv/unpack for B1/B3 and CT update.

### Shearing-box boundaries
- `shearing_box.cpp`: base-class constructor, boundary MB discovery, MPI setup,
  and `FindTargetMB` mapping.
- `shearing_box_cc.cpp`: CC shearing-box pack/send/recv/unpack.
- `shearing_box_fc.cpp`: FC (magnetic field) shearing-box pack/send/recv/unpack.
- `shearing_box_tasks.cpp`: MPI post/clear routines for orbital advection and
  shearing boundaries.

### Source terms
- `shearing_box_srcterms.cpp`: shearing-box momentum/energy sources and the
  rotating-frame E-field correction for 2D MHD.

---

## Activation and Inputs

### Input block
- Presence of `<shearing_box>` activates the module in `SourceTerms` and in
  Hydro/MHD constructors.
- Parameters read in `src/srcterms/srcterms.cpp`:
  - `qshear`
  - `omega0`

### Mesh constraints (enforced in `src/mesh/mesh.cpp`)
- `ix1_bc` and `ox1_bc` must both be `shear_periodic` if either is.
- `shear_periodic` is forbidden on x2 and x3 boundaries.
- 1D is not allowed (requires 2D or 3D).
- Shearing box is not compatible with SMR/AMR (`multilevel`).
- If `shear_periodic` is set but `<shearing_box>` is missing, the code aborts.

### Object construction
- `Hydro` and `MHD` constructors allocate orbital advection and shearing-box
  boundary objects when `<shearing_box>` exists:
  - Hydro: `OrbitalAdvectionCC`, `ShearingBoxBoundaryCC`.
  - MHD: `OrbitalAdvectionCC/FC`, `ShearingBoxBoundaryCC/FC`.

### Dimensional gating
- Most shearing-box tasks are guarded by:
  `pmesh->three_d || psrc->shearing_box_r_phi`.
- `shearing_box_r_phi` defaults to `false` in `SourceTerms` and is not set by any
  input parameter in the current codebase.

---

## Data Layout and Buffers

### Shared buffer container
- `ShearingBoxBoundaryBuffer` owns:
  - `vars`: 5D device buffer for data.
  - `flux`: 5D device buffer (declared, but unused in this directory).
  - MPI request arrays for `vars` (and `flux`, unused).

### Orbital advection buffers
- CC (`OrbitalAdvectionCC`):
  - `sendbuf[n].vars`, `recvbuf[n].vars` shape:
    `(nmb, nvar, nx3, ng + maxjshift, nx1)`.
- FC (`OrbitalAdvectionFC`):
  - `sendbuf[n].vars`, `recvbuf[n].vars` shape:
    `(nmb, 2, nx3+1, ng + maxjshift, nx1+1)`.
  - Stores only B3 and B1 (B2 is not exchanged).
- `maxjshift` computed in `OrbitalAdvection` constructor as:
  `int(cfl_no * max(|x1min|, |x1max|)) + 1`.

### Shearing-box boundary buffers
- CC (`ShearingBoxBoundaryCC`):
  - `sendbuf[n].vars`, `recvbuf[n].vars` shape:
    `(nmb_x1bndry, nx2+2*ng, nvar, nx3+2*ng, ng)`.
- FC (`ShearingBoxBoundaryFC`):
  - `sendbuf[n].vars`, `recvbuf[n].vars` shape:
    `(nmb_x1bndry, nx2+2*ng, 3, nx3+2*ng, ng)`.

### Boundary metadata
- `nmb_x1bndry(0/1)`: number of MeshBlocks touching inner/outer x1 boundaries.
- `x1bndry_mbgid(0/1, m)`: GIDs of boundary MeshBlocks (host + device views).
- `yshear`: current shear distance in x2 computed in `InitRecv`.
- MPI uses dedicated communicators duplicated from `MPI_COMM_WORLD`:
  `comm_orb_advect` for orbital advection and `comm_sbox` for shearing boundaries.

---

## Execution Flow (Hydro/MHD Tasks)

### Orbital advection (OA)
- `InitRecv` (OA): posts MPI receives for x2-face neighbors.
- `SendU_OA`/`SendB_OA`: packs and sends buffers.
- `RecvU_OA`/`RecvB_OA`: tests receives, then applies remap and shift.
- OA tasks run only at the last explicit stage
  (`stage == pdrive->nexp_stages`).

### Shearing-box boundary exchange
- `InitRecv` (shearing): computes `yshear` and posts MPI receives for x1 faces.
- `SendU_Shr`/`SendB_Shr`: applies fractional shift and sends integer-shifted
  buffers to target MBs.
- `RecvU_Shr`/`RecvB_Shr`: waits/tests for receives and copies into x1 ghost zones.

### Source terms and E-field
- `SourceTerms::ShearingBox` is called from Hydro/MHD source-term tasks.
- `SourceTerms::SBoxEField` is called from `MHD::EFieldSrc` in 2D only.

### Stage ordering (stagen task list)
- Hydro order (subset): `CopyCons` -> `Fluxes` -> `SendFlux` -> `RecvFlux` ->
  `RKUpdate` -> `HydroSrcTerms` -> `SendU_OA` -> `RecvU_OA` -> `RestrictU` ->
  `SendU` -> `RecvU` -> `SendU_Shr` -> `RecvU_Shr` -> `ApplyPhysicalBCs` ->
  `Prolongate` -> `ConToPrim` -> `NewTimeStep`.
- MHD order (subset): `CopyCons` -> `Fluxes` -> `SendFlux` -> `RecvFlux` ->
  `RKUpdate` -> `MHDSrcTerms` -> `SendU_OA` -> `RecvU_OA` -> `RestrictU` ->
  `SendU` -> `RecvU` -> `SendU_Shr` -> `RecvU_Shr` -> `CornerE` ->
  `EFieldSrc` -> `SendE` -> `RecvE` -> `CT` -> `SendB_OA` -> `RecvB_OA` ->
  `RestrictB` -> `SendB` -> `RecvB` -> `SendB_Shr` -> `RecvB_Shr` ->
  `ApplyPhysicalBCs` -> `Prolongate` -> `ConToPrim` -> `NewTimeStep`.

---

## Algorithm Notes

### Orbital advection (CC)
- Communicates only x2-face buffers (neighbor indices 8 and 12 in `nghbr`).
- Recv/unpack stage:
  - Constructs a 1D scratch array with left boundary, interior, right boundary.
  - Computes `yshear = -qom * x1 * dt` and integer `joffset`.
  - Uses `DCRemapFlx` or `PLMRemapFlx` for the fractional shift.
  - Updates `a(m,...)` with integer shift plus conservative flux difference.

### Orbital advection (FC)
- Packs only B3 (x3f) and B1 (x1f) into buffers.
- Recv/unpack computes effective EMFs from remap fluxes and integer shifts:
  - `emfx` from B3, `emfz` from B1.
- Updates face-centered fields using CT:
  - `B1` and `B3` updated in multi-D.
  - `B2` always updated (includes 3D term if `three_d`).

### Shearing-box boundaries (CC/FC)
- `yshear = qshear * omega0 * Lx * time` in `ShearingBoxBoundary::InitRecv`.
- For each boundary MB, `joffset = int(yshear / dx2)`; split into:
  - Case 1: `jr < ng` (3 target MBs).
  - Case 2: `jr < nx2 - ng` (2 target MBs).
  - Case 3: otherwise (3 target MBs, wraparound).
- Fractional shift is applied in `PackAndSend*` using remap fluxes; integer shift
  is handled by slicing/redistributing buffers to target MBs.
- `FindTargetMB` wraps target locations in x2 using the logical-location tree.

### Source terms
- `SourceTerms::ShearingBox` adds qshear/omega0-dependent terms to momentum and
  energy.
- Computation uses primitives `w0` (not conserved `u0`).
- In 2D, the component pairing depends on `shearing_box_r_phi`.
- `SBoxEField` adds rotating-frame E-field terms in 2D only; 3D is TODO.

---

## Constraints and Cautions
- Shearing box is incompatible with mesh refinement (`multilevel`).
- Shear-periodic BCs are only supported on x1 faces and require 2D/3D meshes.
- OA and shearing-boundary tasks only run when
  `three_d || shearing_box_r_phi`; `shearing_box_r_phi` is not set by input.
- Remap routines support only `dc` and `plm`. If the global reconstruction
  method is `ppm4`, `ppmx`, or `wenoz`, the remap switch falls through and
  no remap fluxes are computed.
- Many shearing-box routines assume all MeshBlocks share the same `nx2`.
- `ShearingBoxBoundaryBuffer::flux` is allocated but unused in this module.

---

## Related Areas
- `src/mesh/mesh.cpp`: shearing-box boundary checks and AMR restriction.
- `src/srcterms/srcterms.cpp`: parses `<shearing_box>` and sets `qshear`, `omega0`.
- `src/hydro/hydro.cpp` and `src/mhd/mhd.cpp`: allocate OA and shearing BC objects.
- `src/hydro/hydro_tasks.cpp` and `src/mhd/mhd_tasks.cpp`: schedule OA/BC tasks.
- `src/pgen/mri2d.cpp`: 2D shearing-box problem setup and checks.

---

## Extension Points
- Add PPM/WENO remap support in `remap_fluxes.hpp` and update the switches in
  orbital advection and shearing-boundary pack/unpack.
- Implement 3D rotating-frame E-field updates in `SBoxEField`.
- Relax the uniform-`nx2` and no-AMR assumptions if shearing boxes are extended
  to refined meshes.
