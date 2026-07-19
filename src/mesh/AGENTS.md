<!-- BEGIN build-memory-table -->
# Mesh and AMR Navigation

## Scope

This subtree owns the global mesh topology, lightweight MeshBlock metadata,
MeshBlockPack construction, neighbor discovery, load balancing, and static/adaptive mesh
refinement. Physics state remains in the modules attached to `MeshBlockPack`; ordinary
boundary exchange and physical boundary conditions live in `src/bvals/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change root-grid geometry, boundary flags, or mesh-wide state | `src/mesh/mesh.hpp`, `src/mesh/mesh.cpp` | `src/main.cpp`, representative `inputs/` files |
| Change initial or restart topology construction | `src/mesh/build_tree.cpp` | `src/outputs/restart.cpp`, `src/pgen/pgen.cpp` |
| Change logical tree or neighbor discovery | `src/mesh/meshblock_tree.*`, `src/mesh/meshblock.cpp` | `src/mesh/nghbr_index.hpp`, `src/bvals/AGENTS.md` |
| Change pack ownership or physics registration | `src/mesh/meshblock_pack.*` | Module constructors and `src/AGENTS.md` |
| Add or change AMR criteria | `src/mesh/refinement_criteria.*` | `<amr_criterionN>` examples in `inputs/` and AMR regression tests |
| Change refinement, derefinement, or state migration | `src/mesh/mesh_refinement.*` | `src/mesh/prolongation.hpp`, `src/mesh/restriction.hpp`, every affected physics array |
| Change redistribution or MPI AMR transfers | `src/mesh/load_balance.cpp` | `Mesh::LoadBalance`, restart/sharding tests, MPI AMR tests |

## Important flow

1. `main.cpp` constructs `Mesh`, then `BuildTreeFromScratch` or
   `BuildTreeFromRestart` creates the logical tree and rank assignment.
2. `MeshBlockPack` groups the rank-local blocks for Kokkos execution;
   `MeshBlock` stores IDs, levels, extents, boundary flags, and neighbors.
3. `Mesh::AddCoordinatesAndPhysics` attaches coordinates and input-selected physics, then
   creates refinement criteria after the physics arrays exist.
4. Adaptive refinement checks criteria, updates the tree, computes a new load balance,
   migrates/restricts/prolongs each evolved array, rebuilds pack metadata and physics,
   then reinitializes boundaries, primitives, and timestep limits.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/rad/test_rad_lwave1d_amr_cpu.py`:
  focused CPU AMR evolution and cell-centered boundary path.
- `cd tst && python3 run_test_suite.py --test test_suite/sr/test_sr_lwave2d_amr_mpicpu.py`:
  MPI AMR coverage across relativistic hydro and MHD.
- Use the nearest physics-specific AMR test under `tst/test_suite/` when changing how a
  particular module migrates state.

## Local constraints

- `MeshBlock` is deliberately lightweight; do not move bulk cell data out of the physics
  arrays owned by `MeshBlockPack` without changing the execution model.
- Adding an evolved or restart-persistent array requires auditing
  `RedistAndRefineMeshBlocks`, AMR MPI packing, restart write/read layout, and pack rebuilds.
- Cell-centered and face-centered data use different restriction, prolongation, and
  repair paths. MHD face fields require the post-AMR divergence-preserving repair path.
- Keep host/device `DualView` synchronization explicit when host topology decisions feed
  device kernels.
<!-- END build-memory-table -->
