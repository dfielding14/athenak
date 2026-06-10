# Module: Mesh

## Overview

The mesh package turns user input into a block-structured AMR hierarchy.  It
constructs the space-filling curve, distributes MeshBlocks across MPI ranks,
creates GPU-friendly MeshBlockPacks, and manages refinement/coarsening through
`MeshRefinement`.

All heavy data (hydro, MHD, GR state) lives in the physics modules; the mesh layer
primarily tracks topology and geometry.

## Source Layout

`src/mesh/`

| File | Role |
| --- | --- |
| `mesh.hpp/.cpp` | Public `Mesh` class and high-level driver hooks |
| `meshblock.hpp/.cpp` | Lightweight per-block metadata |
| `meshblock_pack.hpp/.cpp` | Pack container and physics module wiring |
| `meshblock_tree.hpp/.cpp` | Octree used for AMR bookkeeping |
| `mesh_refinement.hpp/.cpp` | Refinement criteria, prolongation/restriction |
| `build_tree.cpp` | Root grid construction and MPI partitioning |
| `load_balance.cpp` | Space-filling-curve load balancing |
| `prolongation.hpp`, `restriction.hpp` | Inter-grid transfer operators |

## Mesh Class (selected fields)

```cpp
class Mesh {
  RegionSize  mesh_size;    // physical size of the root grid
  RegionIndcs mesh_indcs;   // cell indices on the root grid
  RegionIndcs mb_indcs;     // indices including ghost zones for each MeshBlock
  BoundaryFlag mesh_bcs[6]; // physical boundary conditions

  int nmb_total;            // total MeshBlocks across all ranks
  int nmb_thisrank;         // MeshBlocks owned by this rank
  int nmb_maxperrank;       // AMR safety limit (input-controlled)

  float *cost_eachmb;       // per-block weighting used by LoadBalance
  int   *rank_eachmb;       // owning rank for each block
  LogicalLocation *lloc_eachmb; // logical coordinates in the AMR tree

  MeshBlockPack *pmb_pack;  // primary pack for this rank
  MeshRefinement *pmr;      // optional AMR controller
};
```

`BuildTreeFromScratch()` (or the restart variant) allocates these buffers, creates
the initial pack, and fills neighbor metadata via `MeshBlock::SetNeighbors`.

## MeshBlock & MeshBlockPack

```cpp
class MeshBlock {
 public:
  DualArray1D<int>          mb_gid;    // global ID
  DualArray1D<int>          mb_lev;    // logical AMR level
  DualArray1D<RegionSize>   mb_size;   // physical extents
  DualArray2D<BoundaryFlag> mb_bcs;    // boundary flags (6 faces)
  DualArray2D<NeighborBlock> nghbr;    // neighbor descriptors (up to 56)
  DualArray1D<bool>         newly_created; // true for blocks created this cycle
  int nnghbr;                          // max neighbors per block (cached)
};

class MeshBlockPack {
 public:
  Mesh *pmesh;
  int nmb_thispack;
  MeshBlock     *pmb;
  Coordinates   *pcoord;
  hydro::Hydro  *phydro;
  mhd::MHD      *pmhd;
  z4c::Z4c      *pz4c;
  adm::ADM      *padm;
  numrel::NumericalRelativity *pnr;
  // ... other optional physics modules (particles, radiation, source terms, etc.)
};
```

Physics modules are allocated on demand inside `MeshBlockPack::AddPhysics()` based
on the presence of blocks in the input deck.

## Input Blocks

### `<mesh>`

Required for every run.

| Parameter | Description |
| --- | --- |
| `nx1`, `nx2`, `nx3` | active cell counts on the root grid |
| `nghost` | ghost-zone depth (default 2) |
| `x#min`, `x#max` | domain extents |
| `i#_bc`, `o#_bc` | boundary conditions (`periodic`, `outflow`, `reflect`, `inflow`, `vacuum`, `shear_periodic`, `user`) |

### `<meshblock>`

Optional override of the block resolution. When omitted, each MeshBlock inherits
the root-grid resolution (`nx#`).  Values must evenly divide the corresponding
root-grid counts.

### `<mesh_refinement>`

Controls AMR.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `refinement` | `"none"` | `"none"`, `"static"`, or `"adaptive"` |
| `num_levels` | 1 | Maximum AMR level (root level + `num_levels` - 1) |
| `max_nmb_per_rank` | required for AMR | Upper bound on blocks per rank (avoids GPU OOM) |
| `ncycle_check` | 1 | Cycles between refinement checks |
| `refinement_interval` | 5 | Minimum cycles between successive refinements |
| `prolong_primitives` | false | When true, prolong primitive variables instead of conserved |
| `dens_max`, `ddens_max`, `dpres_max`, `dvel_max` | optional thresholds | Trigger flags for density/gradient-based refinement |

Problem generators set `MeshRefinement::refine_flag` and `derefine_flag`
directly—there is no Athena++-style `UserRefinementCondition`.

## AMR Workflow

1. **Flagging:** problem generators populate `pmr->refine_flag.h_view(m)` (1 =
   refine, -1 = derefine). Convenience helpers such as
   `MeshRefinement::CheckRefinementCondition` exist in individual generators.
2. **Sync:** flags are mirrored to the device and evaluated every `ncycle_check`
   cycles (respecting `refinement_interval`).
3. **Tree Update:** `MeshRefinement::Apply` creates/destroys blocks, updates
   `newly_created`, and repacks MeshBlockPack.
4. **Load Balance:** `LoadBalance` adjusts the space-filling curve when the cost
   distribution changes significantly.
5. **Transfer Operators:** `prolongation.hpp` and `restriction.hpp` perform
   conservative interpolation/averaging; AMR-aware constrained transport is handled
   inside the MHD module.

## Load Balancing

`float *cost_eachmb` stores a user-adjustable cost for each block (default 1).  The
`LoadBalance` helper computes a prefix sum along the space-filling curve and
assigns contiguous segments to ranks, ensuring `max_nmb_per_rank` is respected.

## Memory / Performance Notes

- All field arrays live in Kokkos `View`s with layout `(nmb, nvar, nk, nj, ni)` so
  the `i` index is fastest varying.
- Pack size is determined automatically from the space-filling curve segment; GPU
  performance is best with several blocks per pack (8–32 typical).
- Ghost-zone exchange uses nonblocking MPI and overlaps with compute via the task
  system.

## Boundary Conditions

`Mesh::GetBoundaryFlag()` recognizes: `periodic`, `outflow`, `reflecting`,
`inflow`, `shear_periodic`, `vacuum`, and `user`. Custom boundary logic is
implemented inside the physics modules or problem generators via the boundary
value package in `src/bvals/`.

## Practical Tips

- Always specify `max_nmb_per_rank` when AMR is enabled—AthenaK enforces this at
  runtime.
- If you shrink `<meshblock>` sizes, ensure they still divide the global counts and
  that each pack maintains enough work for the GPU.
- Use the event log (`outputs/eventlog.cpp`) to monitor `nfofc` and other counters
  that originate in the mesh/refinement layer.

### Ghost Cells
- Default is 2 ghost cells (`nghost=2`)
- Higher-order methods may need 3-4 ghost cells
- Set in input file, not hardcoded

### MeshBlock Indexing
- Indices are local to each MPI rank
- Use global ID (`mb_gid`) for unique identification

### AMR Notes
- MeshBlockPacks rebuilt after refinement
- Load balancing automatic after AMR
- Refinement controlled via problem generators

## Example Usage

### Basic Mesh Setup
```ini
<mesh>
nx1 = 256
x1min = -1.0
x1max = 1.0
ix1_bc = periodic
ox1_bc = periodic

<meshblock>
nx1 = 64   # 256/64 = 4 MeshBlocks in x1
```

### AMR Configuration
```ini
<mesh_refinement>
refinement = adaptive
num_levels = 3
ncycle_check = 10
refinement_interval = 5
```

## See Also
- [Coordinates Module](coordinates.md)
- [Boundary Values Module](boundaries.md)
- Source: `src/mesh/`
