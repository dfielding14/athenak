# AGENTS.md

## Purpose
This directory implements boundary exchange and physical boundary conditions for
mesh variables and particles, including SMR/AMR prolongation and flux correction.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Boundary value framework
- `bvals.hpp`: enums, buffer structs, base class, derived classes for CC/FC, particles.
- `bvals.cpp`: constructors/destructors, buffer init (`InitializeBuffers`), MPI
  communicator setup.

### Task list helpers
- `bvals_tasks.cpp`: `InitRecv`, `ClearRecv`, `ClearSend`, `ClearFluxRecv`,
  `ClearFluxSend` (generic MPI posting/clearing).

### Cell-centered data path
- `bvals_cc.cpp`: `PackAndSendCC` / `RecvAndUnpackCC` for CC data; handles
  same/coarse/fine neighbors plus Z4c same-level coarse append.
- `buffs_cc.cpp`: `InitSendIndices` / `InitRecvIndices` for CC buffers, including
  flux-correction index ranges.
- `flux_correct_cc.cpp`: `PackAndSendFluxCC` / `RecvAndUnpackFluxCC` and
  `InitFluxRecv` for CC flux correction (faces to coarser levels).

### Face-centered data path
- `bvals_fc.cpp`: `PackAndSendFC` / `RecvAndUnpackFC` for face-centered fields,
  including three-component packing.
- `buffs_fc.cpp`: `InitSendIndices` / `InitRecvIndices` for FC buffers; per-component
  index ranges with overlapping face handling for multilevel.
- `flux_correct_fc.cpp`: `PackAndSendFluxFC` / `RecvAndUnpackFluxFC` for edge EMFs;
  `InitFluxRecv`, plus `SumBoundaryFluxes`, `ZeroFluxesAtBoundaryWithFiner`,
  `AverageBoundaryFluxes` helpers.

### Prolongation and primitive conversions
- `prolongation.cpp`: `FillCoarseInBndryCC/FC` and `ProlongateCC/FC` using
  `mesh/prolongation.hpp` and `mesh/restriction.hpp` (Z4c high-order options).
- `prolong_prims.cpp`: `ConsToPrimCoarseBndry` and `PrimToConsFineBndry` (Hydro/MHD
  variants) for primitive-variable prolongation.

### Particle boundary exchange
- `bvals_part.cpp`: `ParticlesBoundaryValues` tasks, send/destroy list handling,
  MPI pack/unpack, and GID updates.

### Physics-specific BCs
- `physics/hydro_bcs.cpp`: `HydroBCs` for reflect/inflow/outflow/diode/vacuum.
- `physics/bfield_bcs.cpp`: `BFieldBCs` for face-centered B with the same BC flags.
- `physics/radiation_bcs.cpp`: `RadiationBCs` (inflow/outflow).
- `physics/z4c_bcs.cpp`: `Z4cBCs` with extrapolation order 2/3/4 and reflection rules.

---

## Core Data Structures
- `BoundaryFace` and `BoundaryFlag` (in `bvals.hpp`) define face identifiers and BC
  types.
- `MeshBoundaryBuffer` stores index ranges (same/coarse/fine/prolong/flux), buffer
  sizes, Kokkos storage (`vars`, `flux`), and MPI request arrays.
- `MeshBoundaryValues` owns 56 send/recv buffers, inflow state arrays
  (`u_in`, `b_in`, `i_in`), MPI communicators, and generic task helpers.
- `MeshBoundaryValuesCC` / `MeshBoundaryValuesFC` implement CC or FC packing,
  prolongation, and flux correction.
- `particles::ParticlesBoundaryValues` manages per-particle send/recv lists and MPI
  buffers.

---

## Boundary Exchange Flow (high level)
1. `InitializeBuffers` uses `NeighborIndex` and mesh indices to set buffer ranges and
   allocate storage. Inflow arrays are allocated when the domain is not strictly
   periodic.
2. Task lists call `InitRecv` (and `InitFluxRecv` during flux correction) to post MPI
   receives.
3. `PackAndSendCC` or `PackAndSendFC` fill send buffers (or copy directly into recv
   buffers for same-rank neighbors) and post MPI sends.
4. `RecvAndUnpackCC/FC` (and `RecvAndUnpackFlux*`) wait/test for completion and unpack
   into ghost zones or coarse buffers.
5. `ClearSend/ClearRecv` (and `ClearFlux*`) finalize MPI requests after the stage.

Notes:
- CC flux corrections are exchanged only on faces with finer neighbors. FC flux
  corrections are exchanged on faces and edges for same-level and finer neighbors,
  then summed/averaged with helper routines.
- For Z4c with multilevel, same-level exchanges append coarse data for higher-order
  prolongation.

---

## Physical Boundary Conditions
- Applied via static helpers in `MeshBoundaryValues`: `HydroBCs`, `BFieldBCs`,
  `RadiationBCs`, `Z4cBCs`.
- `BoundaryFlag` values are stored per MeshBlock (`mb_bcs`) and per mesh (`mesh_bcs`),
  with periodic and shear_periodic handled outside or skipped.
- Z4c BCs use one-sided extrapolation (order 2/3/4) and also update coarse arrays when
  multilevel is enabled.

---

## Extension Points and Cautions
- To add a new BC type, extend `BoundaryFlag` and update each BC helper to keep
  behavior consistent across physics modules.
- When adding new pack/unpack paths, preserve the `NeighborIndex` ordering and buffer
  index conventions in `buffs_cc.cpp` and `buffs_fc.cpp`.
- Prolongation paths depend on coarse buffers being up to date; if you add new
  variables, ensure `FillCoarseInBndry` and `ConsToPrim`/`PrimToCons` coverage.
- Particle boundary handling is separate from mesh variables; update both when
  changing mesh BC semantics.
