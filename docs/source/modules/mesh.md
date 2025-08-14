# Module: Mesh

## Overview
The Mesh module provides the fundamental spatial discretization framework for AthenaK, implementing block-structured adaptive mesh refinement (AMR) with MPI domain decomposition.

## Source Location
`src/mesh/`

## Key Components

| File | Purpose | Key Classes/Functions |
|------|---------|----------------------|
| `mesh.hpp/cpp` | Global mesh management | `Mesh` class |
| `meshblock.hpp/cpp` | Individual mesh block data | `MeshBlock` class |
| `meshblock_pack.hpp/cpp` | GPU-optimized block container | `MeshBlockPack` class |
| `meshblock_tree.hpp/cpp` | Octree structure for AMR | `MeshBlockTree` class |
| `mesh_refinement.hpp/cpp` | Refinement/derefinement logic | `MeshRefinement` class |
| `build_tree.cpp` | Tree construction | `BuildTreeFromScratch()`, `BuildTreeFromRestart()` |
| `load_balance.cpp` | MPI load balancing | `LoadBalance()` |
| `prolongation.hpp` | Coarse-to-fine interpolation | Prolongation operators |
| `restriction.hpp` | Fine-to-coarse averaging | Restriction operators |

## Data Structures

### Mesh Class
```cpp
class Mesh {
  // Core data members (from mesh.hpp)
  RegionSize mesh_size;      // Physical size of mesh
  RegionIndcs mesh_indcs;     // Cell indices
  RegionIndcs mb_indcs;       // MeshBlock cell indices
  BoundaryFlag mesh_bcs[6];   // Boundary conditions
  
  int nmb_total;             // Total MeshBlocks across all ranks
  int nmb_thisrank;          // MeshBlocks on this rank
  int nmb_maxperrank;        // Max MBs per device (memory limit)
  
  float *cost_eachmb;        // Cost of each MeshBlock
  int *rank_eachmb;          // Rank of each MeshBlock
  LogicalLocation *lloc_eachmb; // Logical locations
  
  MeshBlockPack* pmb_pack;   // Container for MeshBlocks
  MeshRefinement *pmr;        // AMR data/functions (if enabled)
};
```

### MeshBlock Class
```cpp
class MeshBlock {
  // Lightweight data structure (from meshblock.hpp)
  DualArray1D<int> mb_gid;           // Global ID
  DualArray1D<int> mb_lev;           // Logical level
  DualArray1D<RegionSize> mb_size;   // Physical size
  DualArray2D<BoundaryFlag> mb_bcs;  // Boundary conditions
  DualArray2D<NeighborBlock> nghbr;  // Neighbor data
  
  // Note: Field data (u, w) stored in physics modules
  // e.g., phydro->u0, pmhd->u0, etc.
};
```

### MeshBlockPack Class
```cpp
class MeshBlockPack {
  Mesh *pmesh;               // Pointer to parent Mesh
  int nmb_thispack;          // Number of MBs in this pack
  
  MeshBlock* pmb;            // MeshBlocks in this pack
  Coordinates* pcoord;       // Coordinate system
  
  // Physics modules (allocated as needed)
  hydro::Hydro *phydro;
  mhd::MHD *pmhd;
  z4c::Z4c *pz4c;
  particles::Particles *ppart;
  // ... etc
};
```

## Configuration Parameters

### Mesh Block (`<mesh>`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nx1, nx2, nx3` | int | required | Number of cells in domain |
| `x1min, x1max` | Real | required | x1 domain boundaries |
| `x2min, x2max` | Real | required | x2 domain boundaries |
| `x3min, x3max` | Real | required | x3 domain boundaries |
| `ix1_bc, ox1_bc` | string | required | Inner/outer x1 boundary conditions |
| `ix2_bc, ox2_bc` | string | required | Inner/outer x2 boundary conditions |
| `ix3_bc, ox3_bc` | string | required | Inner/outer x3 boundary conditions |

### MeshBlock Block (`<meshblock>`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nx1, nx2, nx3` | int | required | Cells per MeshBlock |

### Mesh Refinement Block (`<mesh_refinement>`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `refinement` | string | "none" | "none", "static", or "adaptive" |
| `num_levels` | int | 1 | Maximum number of refinement levels |
| `ncycle_check` | int | 1 | Cycles between AMR checks |
| `refinement_interval` | int | 5 | Cycles between successive refinements |

## Key Functions

### Tree Construction
```cpp
// From build_tree.cpp
void Mesh::BuildTreeFromScratch(ParameterInput *pin);
void Mesh::BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile, 
                                bool single_file_per_rank);
```
**Note**: No generic `BuildTree()` function exists.

### Load Balancing
```cpp
// From mesh.hpp (private member)
void LoadBalance(float *clist, int *rlist, int *slist, int *nlist, int nb);
```
Redistributes MeshBlocks across MPI ranks using a space-filling curve.

### MeshBlockPack Creation
```cpp
// Created in mesh.cpp constructor
MeshBlockPack(Mesh *pm, int igids, int igide);
```
Groups MeshBlocks for GPU efficiency. Pack size typically 1-32 blocks.

## AMR Operations

### Refinement Control
AMR is controlled through problem generators, not a user-defined function:

```cpp
// In problem generator (pgen/*.cpp)
// Set refinement flags during initialization or evolution
if (need_refinement) {
  pmb->pmr->refine_flag.h_view(m) = 1;  // Flag for refinement
}
```

**Note**: There is NO `UserRefinementCondition()` function in AthenaK.

### Prolongation and Restriction
- **Prolongation**: Conservative linear interpolation (in `prolongation.hpp`)
- **Restriction**: Volume-weighted averaging (in `restriction.hpp`)
- Both operators preserve conservation

### Recent AMR + div(B) Improvements (2024)

#### MHD AMR div(B) Preservation
1. **Constrained Transport**: Extended to AMR interfaces
2. **Face-centered field prolongation**: Preserves divergence
3. **Test files verified to exist**:
   - `inputs/mhd/test_divb_minimal.athinput`
   - `inputs/mhd/test_divb_ranks.athinput`

## Memory Layout

### Array Layout
- **Kokkos Views**: All arrays use LayoutRight for GPU efficiency
- **Field arrays**: `(nmb, nvar, nk, nj, ni)` where:
  - nmb: MeshBlock index within pack
  - nvar: Variable index
  - nk, nj, ni: z, y, x indices

### MeshBlockPack Optimization
- Groups multiple MeshBlocks for coalesced GPU memory access
- Typical size: 1-32 blocks per pack
- All blocks in pack execute same kernels

## Boundary Conditions

Supported types (from `mesh.cpp` GetBoundaryFlag()):
- `periodic`: Periodic boundary
- `outflow`: Zero gradient
- `reflecting`: Reflecting (symmetry)
- `inflow`: Fixed inflow values
- `vacuum`: Vacuum boundary
- `user`: User-defined in problem generator

## Performance Considerations

### GPU Optimization
- MeshBlockPacks provide coalesced memory access
- Pack size tuned for GPU occupancy
- Shared coordinate arrays within packs

### MPI Communication
- Non-blocking MPI for neighbor exchange
- Load balancing via space-filling curve
- Ghost cell exchange overlapped with computation

## Common Issues

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