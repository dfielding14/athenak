# Turbulence Driver AMR Derefinement Issue Analysis

## Problem Summary

The turbulence driver in AthenaK fails during AMR derefinement because:

1. **Array Size Mismatch**: Turbulence driver arrays (`force`, `force_tmp1`, `force_tmp2`, etc.) are allocated based on the initial number of MeshBlocks (`nmb`) in the constructor
2. **No Dynamic Resizing**: These arrays are never resized when the mesh undergoes refinement or derefinement
3. **Invalid Indexing**: When MeshBlocks are destroyed during derefinement, the turbulence driver still tries to access array indices based on the new MeshBlock count, causing out-of-bounds access

## Root Cause

In `turb_driver.cpp` constructor (lines 54-63):
```cpp
// allocate memory for force registers
int nmb = pmy_pack->nmb_thispack;  // Initial MeshBlock count
...
Kokkos::realloc(force, nmb, 3, ncells3, ncells2, ncells1);
Kokkos::realloc(force_tmp1, nmb, 3, ncells3, ncells2, ncells1);
Kokkos::realloc(force_tmp2, nmb, 3, ncells3, ncells2, ncells1);
```

These arrays are sized for the initial mesh configuration and never updated.

## Impact

During the `Generate()` task (line 767+):
```cpp
int nmb = pmy_pack->nmb_thispack;  // Current MeshBlock count (may have changed!)
...
par_for("force_update", DevExeSpace(),0,nmb-1,ks,ke,js,je,is,ie,
KOKKOS_LAMBDA(int m, int k, int j, int i) {
    force_(m,0,k,j,i) = ...  // May access out-of-bounds if nmb > original size
```

## Diagnostic Test

Created `turb_amr_wave_test.cpp` that:
- Implements a traveling wave refinement criterion
- Forces periodic refinement and derefinement as the wave moves
- Exposes array sizing issues when MeshBlock count changes

## Proposed Fixes

### Option 1: Dynamic Array Resizing
Add a method to resize turbulence arrays when mesh changes:
```cpp
void TurbulenceDriver::ResizeArrays(int new_nmb) {
  if (new_nmb != force.extent(0)) {
    // Resize all arrays
    Kokkos::realloc(force, new_nmb, 3, ncells3, ncells2, ncells1);
    Kokkos::realloc(force_tmp1, new_nmb, 3, ncells3, ncells2, ncells1);
    // ... etc
    
    // Reinitialize if needed
    Initialize();
  }
}
```

### Option 2: Over-allocate Arrays
Allocate arrays with a safety margin:
```cpp
int max_nmb = nmb * 8;  // Allow up to 3 levels of refinement
Kokkos::realloc(force, max_nmb, 3, ncells3, ncells2, ncells1);
```

### Option 3: Use Dynamic Views
Switch to dynamically-sized Kokkos Views that can handle resizing.

## Testing Instructions

1. Build with the diagnostic test:
```bash
cmake -B build -DPROBLEM=turb_amr_wave_test -DAthena_ENABLE_MPI=ON
cd build && make -j8
```

2. Run the test:
```bash
cd build/src
./athena -i ../../inputs/hydro/turb_amr_wave_demo.athinput
```

3. Analyze results:
```bash
python ../../vis/python/diagnose_turb_amr_derefinement.py .
```

## Expected Behavior

Without fix:
- Segmentation fault or incorrect force values when derefinement occurs
- Force field discontinuities at refinement boundaries
- Array index out-of-bounds errors

With fix:
- Smooth force field evolution across refinement/derefinement
- Consistent array sizes with MeshBlock count
- No crashes or discontinuities