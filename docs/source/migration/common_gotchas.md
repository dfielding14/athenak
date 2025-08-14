# Common Gotchas

## Migration Issues from Athena++

### 1. Extra MeshBlock Index
**Issue**: All arrays have an extra dimension for MeshBlock index
```cpp
// Athena++
prim(IDN,k,j,i)

// AthenaK  
u0(m,IDN,k,j,i)  // 'm' is MeshBlock index, variable index is second
```

### 2. Device Memory
**Issue**: Accessing device memory from host code
```cpp
// WRONG - will crash
Real value = u0(0,IDN,0,0,0);  // Can't access from CPU

// CORRECT
Kokkos::deep_copy(host_view, device_view);
Real value = host_view(0,IDN,0,0,0);
```

### 3. Lambda Capture
**Issue**: Capturing by reference in Kokkos kernels
```cpp
// WRONG
Real factor = 2.0;
par_for(..., KOKKOS_LAMBDA(int i) {
  data(i) *= factor;  // Undefined behavior
});

// CORRECT
const Real factor = 2.0;  // Make const
par_for(..., KOKKOS_LAMBDA(int i) {
  data(i) *= factor;  // OK
});
```

### 4. Boundary Conditions
**Issue**: Ghost zones indexed differently
- Athena++: Negative indices for ghost zones
- AthenaK: Explicit ghost zone ranges

### 5. MPI + GPU
**Issue**: MPI not GPU-aware
```bash
# May need to set
export MPICH_GPU_SUPPORT_ENABLED=1
```

## Performance Gotchas

### 1. Small MeshBlocks on GPU
**Issue**: MeshBlocks too small for GPU efficiency
```ini
<meshblock>
nx1 = 32  # Better for GPU (default is mesh nx1)
```

### 2. Memory Transfers
**Issue**: Frequent CPU-GPU transfers
- Minimize deep_copy operations
- Keep data on device

### 3. Task Granularity
**Issue**: Too many small tasks
- Batch operations when possible
- Use MeshBlockPacks

## Common Errors

### Build Errors
```bash
# Kokkos not found
git submodule update --init

# CUDA not detected
export CUDA_ROOT=/usr/local/cuda
```

### Runtime Errors
```bash
# Segfault on GPU
# Check for race conditions
# Verify array bounds

# Wrong results
# Check reduction operations
# Verify atomic operations
```

## See Also
- [Migration Guide](from_athena_plus_plus.md)
- [Troubleshooting Guide](../troubleshooting.md)