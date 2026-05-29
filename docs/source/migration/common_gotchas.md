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

### 3. Lambda Capture and Device-Accessible State
**Issue**: Capturing host-only state or host stack references in device kernels

`KOKKOS_LAMBDA` captures ordinary local variables by value. Scalars used as
read-only kernel parameters do not need to be declared `const`.
```cpp
// CORRECT - factor is copied into the kernel closure
Real factor = 2.0;
par_for(..., KOKKOS_LAMBDA(int i) {
  data(i) *= factor;
});

// WRONG for portable device execution - factor_ptr points to host stack storage
Real *factor_ptr = &factor;
par_for(..., KOKKOS_LAMBDA(int i) {
  data(i) *= *factor_ptr;
});
```

For outer Kokkos kernels, use values captured by copy for scalar parameters and
device-accessible Kokkos Views for array data; avoid reference capture such as
`[&]`. Nested team/vector lambdas are a separate case: reference capture may be
used only where permitted by Kokkos and where the computation is also correct
under capture-by-copy semantics.

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
- [Migration Guide](index.md)
- [Troubleshooting Guide](../troubleshooting.md)
