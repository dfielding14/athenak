# AMR-Compatible Turbulence Driver Implementation

## Overview

This document describes the modifications made to the turbulence driver in AthenaK to support Adaptive Mesh Refinement (AMR). The turbulence driver implements randomly forced turbulence using an Ornstein-Uhlenbeck stochastic process with Cartesian basis functions.

## Key Challenges with AMR

When AMR is enabled, the mesh structure changes dynamically:
1. **MeshBlock count changes**: Refinement creates new MeshBlocks (8 children from 1 parent), derefinement destroys MeshBlocks (8 children merge to 1 parent)
2. **MeshBlock positions change**: After refinement, MeshBlocks may have different spatial locations
3. **Array sizes must adapt**: All arrays indexed by MeshBlock must be resized
4. **Basis functions need recalculation**: Position-dependent basis functions must be recomputed for new MeshBlock locations

## Implementation Details

### 1. Dynamic Array Resizing (`ResizeArrays` function)

The `ResizeArrays` function (lines 602-806 in `turb_driver.cpp`) handles all array resizing when the mesh changes:

```cpp
void TurbulenceDriver::ResizeArrays(int new_nmb)
```

#### Key Features:
- **Refinement Detection**: Tracks `nmb_created` and `nmb_deleted` counters from the mesh refinement module
- **Conditional Resizing**: Only resizes when MeshBlock count actually changes
- **Comprehensive Array Updates**: Resizes all turbulence-related arrays:
  - Force arrays: `force`, `force_tmp1`, `force_tmp2` 
  - Mode arrays: `xcos`, `xsin`, `ycos`, `ysin`, `zcos`, `zsin`

### 2. MeshBlock Tracking

The turbulence driver maintains several tracking variables:

```cpp
int current_nmb = 0;      // Current number of MeshBlocks
int last_nmb_created = 0; // Last known nmb_created value
int last_nmb_deleted = 0; // Last known nmb_deleted value
```

These track refinement events and prevent unnecessary array operations.

### 3. Global Coordinate Recalculation

After resizing, all position-dependent arrays are recalculated using **global coordinates** (lines 708-744):

```cpp
// Recalculate for ALL MeshBlocks with global coordinates
par_for("xsin/xcos", DevExeSpace(),0,new_nmb-1,0,mode_count-1,is,ie,
KOKKOS_LAMBDA(int m, int n, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v_global = CellCenterX(i-is, nx1, x1min, x1max);
    Real k1v = kx_mode_.d_view(n);
    xsin_(m,n,i) = sin(k1v*x1v_global);
    xcos_(m,n,i) = cos(k1v*x1v_global);
});
```

This ensures basis functions are correct for each MeshBlock's actual position in the domain.

### 4. Force Preservation During Refinement

The implementation preserves forcing continuity across refinement by:
1. **Preserving mode amplitudes**: The `aka` and `akb` arrays (Fourier amplitudes) are NOT resized or modified
2. **Recomputing with same amplitudes**: After resizing, forces are recalculated using the same amplitudes but new basis functions (lines 768-785)
3. **Maintaining temporal coherence**: The Ornstein-Uhlenbeck evolution continues uninterrupted

### 5. Array Dimensions

All major arrays now include MeshBlock dimension:

- **Before AMR support**: 
  - `force[3][nk][nj][ni]` - force components per cell
  - `xcos[mode_count][ni]` - cosine basis per mode and x-position

- **After AMR support**:
  - `force[nmb][3][nk][nj][ni]` - includes MeshBlock index
  - `xcos[nmb][mode_count][ni]` - basis functions per MeshBlock

### 6. Task Integration

The turbulence driver integrates with the task system through:

```cpp
TaskStatus TurbulenceDriver::CheckResize(Driver *pdrive, int stage) {
    int nmb = pmy_pack->nmb_thispack;
    if (nmb != current_nmb || force.extent(0) != nmb) {
        ResizeArrays(nmb);
    }
    return TaskStatus::complete;
}
```

This task runs before force updates to ensure arrays are correctly sized.

## Initialization and Configuration

### Constructor Modifications

The constructor (lines 33-72) initializes AMR tracking:

```cpp
// Set initial MeshBlock count
current_nmb = nmb;

// Initialize refinement tracking variables  
if (pm->adaptive && pm->pmr != nullptr) {
    last_nmb_created = pm->pmr->nmb_created;
    last_nmb_deleted = pm->pmr->nmb_deleted;
} else {
    last_nmb_created = 0;
    last_nmb_deleted = 0;
}
```

### Mode Selection

The driver uses Cartesian Fourier modes with wavenumbers `(kx, ky, kz)`. Mode counting accounts for the selected wavenumber range (lines 142-168).

## Critical Implementation Decisions

### 1. Global vs Local Coordinates
**Decision**: Use global coordinates for all basis functions
**Rationale**: Ensures forcing pattern continuity across MeshBlock boundaries

### 2. Array Resizing Strategy
**Decision**: Resize all arrays when MeshBlock count changes
**Rationale**: Simpler than complex reindexing schemes, minimal performance impact

### 3. Force Amplitude Preservation
**Decision**: Keep `aka` and `akb` arrays unchanged during refinement
**Rationale**: Maintains forcing statistics and temporal evolution

### 4. Refinement Detection Method
**Decision**: Track both MeshBlock count and refinement counters
**Rationale**: Catches all refinement events reliably

## Performance Considerations

1. **Memory Overhead**: Arrays scale with number of MeshBlocks × modes × cells
2. **Recomputation Cost**: Basis functions recalculated for all MeshBlocks after refinement
3. **Parallel Efficiency**: All recalculations use Kokkos parallel_for for GPU acceleration

## Testing and Validation

The implementation was tested with:
- **Problem generator**: `turb_timed_amr` 
- **Input file**: `turb_timed_amr_stage1.athinput`
- **Refinement trigger**: Time-based refinement at t=2.0
- **MeshBlock scaling**: 8 blocks → 64 blocks on refinement

## Known Issues and Limitations

1. **Memory Scaling**: Large mode counts with many MeshBlocks can consume significant memory
2. **Load Balancing**: No special handling for load imbalance after refinement

## Future Improvements

1. **Lazy Recalculation**: Only recalculate basis functions for affected MeshBlocks
2. **Memory Optimization**: Consider memory pooling for temporary arrays
3. **Interpolation**: Add force interpolation for refined cells instead of recalculation

## Summary

The AMR-compatible turbulence driver successfully maintains forced turbulence across dynamic mesh refinement by:
- Dynamically resizing all arrays when MeshBlock count changes
- Recalculating position-dependent basis functions using global coordinates
- Preserving forcing amplitudes to maintain statistical properties
- Integrating cleanly with the AthenaK task system

This implementation allows turbulence simulations to benefit from AMR's computational efficiency while maintaining physically consistent forcing throughout the domain.