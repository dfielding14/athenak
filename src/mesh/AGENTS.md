# AGENTS.md

## Purpose
This directory defines AthenaK's mesh hierarchy: the global `Mesh`, per-block metadata,
AMR tree/refinement logic, neighbor discovery, and MPI load balancing. It is responsible
for building the MeshBlock layout (from scratch or restart), tracking logical locations,
and maintaining neighbor relationships used by boundary exchange and AMR.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core mesh types
- `mesh.hpp` / `mesh.cpp`: `Mesh` class and global mesh state (domain size, boundary
  flags, meshblock counts, time/cycle, restart metadata).
- `meshblock.hpp` / `meshblock.cpp`: `MeshBlock` metadata (gid/level, physical extents,
  boundary flags, neighbor lists) stored in DualArrays for host/device use.
- `meshblock_pack.hpp` / `meshblock_pack.cpp`: `MeshBlockPack` container that groups
  MeshBlocks, owns `Coordinates`, and constructs physics modules + task lists.

### Logical topology and refinement
- `meshblock_tree.hpp` / `meshblock_tree.cpp`: `MeshBlockTree` nodes (logical AMR tree)
  used for neighbor lookup and layout bookkeeping.
- `mesh_refinement.hpp` / `mesh_refinement.cpp`: `MeshRefinement` (SMR/AMR orchestration,
  refine/derefine decisions, data prolongation/restriction, AMR communications).

### Build and redistribution
- `build_tree.cpp`: `Mesh::BuildTreeFromScratch` and `Mesh::BuildTreeFromRestart`.
- `load_balance.cpp`: MPI load balancing for uniform meshes and AMR redistribution.

### AMR helpers
- `prolongation.hpp`: inline prolongation operators for CC/FC data.
- `restriction.hpp`: inline restriction operators (templated by NGHOST).
- `nghbr_index.hpp`: neighbor indexing scheme for faces/edges/corners.

---

## How the Mesh Gets Built

### New runs
`Mesh::BuildTreeFromScratch` (in `build_tree.cpp`) does the following:
1. Computes root-grid layout and logical levels.
2. Builds the `MeshBlockTree` and expands it for SMR/AMR refinement regions.
3. Allocates `cost_eachmb`, `rank_eachmb`, `lloc_eachmb`, and rank offsets.
4. Calls `Mesh::LoadBalance` to assign MeshBlocks to MPI ranks.
5. Constructs `MeshBlockPack`, `MeshBlock`s, and neighbor lists.
6. Instantiates `MeshRefinement` if multilevel refinement is enabled.

### Restarts
`Mesh::BuildTreeFromRestart` rehydrates mesh metadata from the restart file, rebuilds
the logical tree, assigns MeshBlocks to ranks, and then follows the same pack/block
construction path.
- For `single_file_per_rank` restarts, pass the mode flag through all `IOWrapper`
  calls (including `GetPosition(single_file_per_rank)`), otherwise MPI-enabled
  builds can hit MPI-IO position errors while reading rank-local files.
- Restart metadata (`rank_eachmb`, `gids_eachrank`, `nmb_eachrank`) is consumed
  downstream by `src/pgen/pgen.cpp` to remap source-rank data when restoring
  per-rank restart files onto a possibly different runtime rank layout.

---

## Key Invariants and Rules
- `MeshBlockPack::AddCoordinates` must be called before `AddPhysics` (physics constructors
  use coordinate data).
- `MeshBlock::SetNeighbors` uses `nghbr_index.hpp` to map face/edge/corner neighbors.
- AMR requires MeshBlock dimensions divisible by 2 and a configured
  `<mesh_refinement>/max_nmb_per_rank`.
- The Mesh constructor enforces boundary consistency (e.g., periodic pairs, shearing-box
  constraints) and rejects incompatible options (e.g., shearing box with refinement).

---

## Extension Points
- **New refinement criteria:** add checks in `MeshRefinement::CheckForRefinement` and
  wire parameters in `mesh_refinement.cpp`.
- **Custom load-balancing cost model:** update the cost initialization in
  `BuildTreeFromScratch` or extend `Mesh::LoadBalance`.
- **Neighbor logic:** adjust `NeighborIndex` or `MeshBlock::SetNeighbors` if the
  connectivity scheme changes.

---

## Related Areas
- Physics modules are instantiated in `MeshBlockPack::AddPhysics` using `<hydro>`,
  `<mhd>`, `<radiation>`, `<particles>`, `<z4c>`, etc. input blocks.
- Boundary exchanges and AMR ghost-zone handling also depend on the neighbor mapping
  and refinement operators defined here.
