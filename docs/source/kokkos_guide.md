# Kokkos Programming Guide for AthenaK

## Overview

This guide explains Kokkos concepts essential for understanding and developing AthenaK code. Kokkos enables performance portability across CPUs and GPUs with a single code base. This guide focuses on practical concepts needed for AthenaK development.

## Core Concepts

### 1. Execution and Memory Spaces

Kokkos separates **where code runs** (execution space) from **where data lives** (memory space):

```cpp
// Execution Spaces
Kokkos::Serial         // Single CPU thread
Kokkos::OpenMP         // CPU with OpenMP threads  
Kokkos::Cuda           // NVIDIA GPUs
Kokkos::HIP            // AMD GPUs

// Memory Spaces
Kokkos::HostSpace      // CPU RAM
Kokkos::CudaSpace      // GPU memory
Kokkos::HIPSpace       // AMD GPU memory
```

In AthenaK, we use aliases:
```cpp
using DevExeSpace = Kokkos::DefaultExecutionSpace;  // GPU if available, else CPU
using HostExeSpace = Kokkos::Serial;               // Always CPU
using DevMemSpace = DevExeSpace::memory_space;     // Matches execution space
```

### 2. Kokkos Views (Arrays)

Views are Kokkos' multi-dimensional arrays that know their memory space:

```cpp
// AthenaK type aliases for Views
using DvceArray1D = Kokkos::View<Real*, LayoutWrapper, DevMemSpace>;
using DvceArray5D = Kokkos::View<Real*****, LayoutWrapper, DevMemSpace>;

// Creating Views
DvceArray5D u("conserved", nmb, nvar, nk, nj, ni);  // Device memory
auto u_host = Kokkos::create_mirror_view(u);        // Host mirror

// Accessing data
u(m, n, k, j, i) = value;  // Device access (inside kernels)
u_host(m, n, k, j, i) = value;  // Host access (outside kernels)
```

**Key Points:**
- Views manage memory lifetime automatically (reference counting)
- Access uses parentheses `()`, not brackets `[]`
- Views can be passed by value to kernels (they're lightweight handles)

### 3. Device-Host Communication

Data must be explicitly copied between host and device:

```cpp
// Create device and host views
DvceArray3D device_data("data", nz, ny, nx);
auto host_data = Kokkos::create_mirror_view(device_data);

// Initialize on host
for (int k=0; k<nz; k++) {
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      host_data(k,j,i) = initial_value;
    }
  }
}

// Copy host → device
Kokkos::deep_copy(device_data, host_data);

// Process on device
par_for("process", DevExeSpace(), 0, nx,
  KOKKOS_LAMBDA(int i) {
    device_data(0,0,i) *= 2.0;
  });

// Copy device → host to read results
Kokkos::deep_copy(host_data, device_data);
```

**Important**: Never access host memory from device kernels or vice versa!

### 4. Parallel Patterns

AthenaK uses Kokkos parallel patterns extensively:

#### Parallel For
Most common pattern for independent iterations:
```cpp
// 1D loop
par_for("label", DevExeSpace(), 0, n,
  KOKKOS_LAMBDA(int i) {
    array(i) = i * 2.0;
  });

// 3D nested loop
par_for("label", DevExeSpace(), 
  {0, nz}, {0, ny}, {0, nx},
  KOKKOS_LAMBDA(int k, int j, int i) {
    u(k,j,i) = initial_value;
  });

// AthenaK's typical 6D loop (MeshBlock + 3D space)
par_for("hydro", DevExeSpace(),
  0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u(m,IDN,k,j,i) = density;
  });
```

#### Parallel Reduce
For reductions (sum, min, max):
```cpp
Real total_mass = 0.0;
par_reduce("mass", DevExeSpace(),
  0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i, Real &sum) {
    sum += u(m,IDN,k,j,i) * vol(m,k,j,i);
  }, total_mass);
```

### 5. Lambda Captures and KOKKOS_LAMBDA

The `KOKKOS_LAMBDA` macro ensures code works on both CPU and GPU:

```cpp
// KOKKOS_LAMBDA expands to:
// - [=] on CPU (capture by value)
// - [=] __device__ on GPU (device function)
```

#### What Can Be Captured?

✅ **Safe to Capture:**
- Kokkos Views (by value)
- Simple POD types (int, Real, bool)
- Enums and constants

```cpp
DvceArray3D data("data", nz, ny, nx);
Real coefficient = 2.0;
int offset = 5;

par_for("safe", DevExeSpace(), 0, nx,
  KOKKOS_LAMBDA(int i) {
    data(0,0,i) = coefficient * (i + offset);  // All safe
  });
```

❌ **Cannot Capture:**
- Raw pointers to host memory
- STL containers (vector, map, etc.)
- Class objects with virtual functions
- Host-only Views

```cpp
// WRONG - Will crash on GPU!
std::vector<Real> host_vec(100);
Real* host_ptr = new Real[100];

par_for("bad", DevExeSpace(), 0, nx,
  KOKKOS_LAMBDA(int i) {
    value = host_vec[i];    // ERROR: STL container
    value = host_ptr[i];    // ERROR: Host pointer
  });
```

### 6. Memory Management Best Practices

#### View Allocation
```cpp
// Allocate once, reuse many times
DvceArray3D scratch;  // Declare outside loops

// In initialization
Kokkos::realloc(scratch, nz, ny, nx);  // Allocate/resize

// Reuse in time loop without reallocation
for (int cycle=0; cycle<ncycles; cycle++) {
  par_for(..., KOKKOS_LAMBDA(int i) {
    scratch(k,j,i) = compute_value();
  });
}
```

#### Minimize Host-Device Transfers
```cpp
// BAD: Transfer every timestep
for (int n=0; n<1000; n++) {
  Kokkos::deep_copy(device_data, host_data);  // Slow!
  process_on_device(device_data);
  Kokkos::deep_copy(host_data, device_data);  // Slow!
}

// GOOD: Keep data on device
Kokkos::deep_copy(device_data, host_data);  // Once
for (int n=0; n<1000; n++) {
  process_on_device(device_data);  // All on device
}
Kokkos::deep_copy(host_data, device_data);  // Once
```

### 7. Scratch Memory

For temporary workspace inside kernels:
```cpp
// Team-level scratch memory
using ScratchPad = Kokkos::View<Real*, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
size_t scratch_size = ScratchPad::shmem_size(scratch_level_size);

par_for_outer("team", TeamPolicy<>(DevExeSpace(), nmb, Kokkos::AUTO),
  scratch_size, KOKKOS_LAMBDA(TeamMember member) {
    ScratchPad scratch(member.team_scratch(0), scratch_level_size);
    
    // Use scratch array within team
    par_for_inner(member, 0, ni, [&](int i) {
      scratch(i) = compute_temp_value(i);
    });
  });
```

### 8. Common AthenaK Patterns

#### MeshBlock Pack Pattern
All physics modules work on packs of MeshBlocks:
```cpp
void HydroUpdate(MeshBlockPack *pmbp) {
  auto &u = pmbp->phydro->u;  // 5D array: (nmb,nvar,nk,nj,ni)
  int nmb = pmbp->nmb_thispack;
  
  // Process all blocks in parallel
  par_for("hydro", DevExeSpace(),
    0, nmb-1, pmbp->ks, pmbp->ke, 
    pmbp->js, pmbp->je, pmbp->is, pmbp->ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      // m = MeshBlock index in pack
      u(m,IDN,k,j,i) = update_density();
    });
}
```

#### Task Pattern
Tasks operate on entire MeshBlock packs:
```cpp
TaskStatus MyTask(Driver *pdrive, int stage) {
  auto &pmbp = pdrive->pmesh->pmb_pack;
  auto &u = pmbp->phydro->u;
  
  // Launch kernel for all blocks
  par_for("task", DevExeSpace(), 0, pmbp->nmb_thispack-1,
    KOKKOS_LAMBDA(int m) {
      // Work on block m
    });
  
  // Synchronize if needed
  Kokkos::fence();
  
  return TaskStatus::complete;
}
```

## Performance Optimization

### GPU-Specific Optimizations

#### 1. Coalesced Memory Access
Keep fastest-varying index innermost:
```cpp
// GOOD: Coalesced access on GPU
par_for(..., KOKKOS_LAMBDA(int k, int j, int i) {
  u(m,n,k,j,i) = value;  // i is innermost → contiguous
});

// BAD: Strided access on GPU  
par_for(..., KOKKOS_LAMBDA(int i, int j, int k) {
  u(m,n,k,j,i) = value;  // k varies fastest but not innermost
});
```

#### 2. Minimize Divergence
Avoid complex branching in kernels:
```cpp
// BAD: Divergent branches
KOKKOS_LAMBDA(int i) {
  if (complex_condition(i)) {
    // Long computation A
  } else {
    // Long computation B  
  }
}

// BETTER: Separate kernels
if (global_condition) {
  par_for("A", ..., KOKKOS_LAMBDA(int i) { /* A */ });
} else {
  par_for("B", ..., KOKKOS_LAMBDA(int i) { /* B */ });
}
```

#### 3. Avoid Atomics When Possible
Design algorithms to avoid atomic operations:
```cpp
// BAD: Atomic updates (serializes)
KOKKOS_LAMBDA(int i) {
  Kokkos::atomic_add(&global_sum, local_value(i));
}

// GOOD: Use reduction
par_reduce(..., KOKKOS_LAMBDA(int i, Real &sum) {
  sum += local_value(i);
}, global_sum);
```

### CPU-Specific Optimizations

#### Vectorization
Structure loops for SIMD:
```cpp
// Enable vectorization hints
#pragma ivdep
for (int i=is; i<=ie; i++) {
  u(i) = v(i) + w(i);
}
```

## Debugging Tips

### 1. Test on CPU First
```bash
# Build for CPU only
cmake -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=OFF ..
```

### 2. Use Bounds Checking
```bash
# Enable bounds checking
cmake -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON ..
```

### 3. Add Fences for Debugging
```cpp
par_for("kernel1", ...);
Kokkos::fence();  // Force synchronization
// Check results here
par_for("kernel2", ...);
```

### 4. Print from Limited Threads
```cpp
KOKKOS_LAMBDA(int i) {
  if (i == 0) {  // Only thread 0 prints
    printf("Debug: value = %f\n", value);
  }
}
```

## Common Pitfalls and Solutions

### Pitfall 1: Accessing Host Memory from Device
```cpp
// WRONG
Real host_array[100];
par_for(..., KOKKOS_LAMBDA(int i) {
  value = host_array[i];  // CRASH on GPU!
});

// CORRECT
DvceArray1D device_array("arr", 100);
// ... initialize ...
par_for(..., KOKKOS_LAMBDA(int i) {
  value = device_array(i);  // Safe
});
```

### Pitfall 2: Race Conditions
```cpp
// WRONG: Race condition
par_for(..., KOKKOS_LAMBDA(int i) {
  shared_sum += array(i);  // Multiple threads write
});

// CORRECT: Use reduction
par_reduce(..., KOKKOS_LAMBDA(int i, Real &sum) {
  sum += array(i);  // Thread-safe
}, total);
```

### Pitfall 3: Forgetting MeshBlock Index
```cpp
// WRONG: Missing 'm' dimension
u(IDN,k,j,i) = density;

// CORRECT: Include MeshBlock index
u(m,IDN,k,j,i) = density;
```

### Pitfall 4: Wrong Execution Space
```cpp
// WRONG: Host execution for device data
par_for("loop", HostExeSpace(), ...,
  KOKKOS_LAMBDA(int i) {
    device_view(i) = 0;  // May fail
  });

// CORRECT: Match execution and memory
par_for("loop", DevExeSpace(), ...,
  KOKKOS_LAMBDA(int i) {
    device_view(i) = 0;  // Correct
  });
```

## AthenaK-Specific Conventions

### Naming Conventions
- `Dvce*` prefix: Device memory arrays
- `Host*` prefix: Host memory arrays  
- `par_for`: Kokkos parallel for wrapper
- `par_reduce`: Kokkos parallel reduce wrapper
- `KOKKOS_LAMBDA`: Device-compatible lambda

### Standard Indices
```cpp
m   // MeshBlock index in pack (0 to nmb-1)
n   // Variable index (IDN, IM1, IM2, IM3, IEN)
k   // Z-direction index
j   // Y-direction index  
i   // X-direction index
```

### Memory Layout
AthenaK uses column-major (Fortran) layout:
- Rightmost index varies fastest in memory
- Optimal access pattern: `(m,n,k,j,i)` with `i` in innermost loop

## Quick Reference Card

```cpp
// Essential includes
#include <Kokkos_Core.hpp>
#include "athena.hpp"

// Type aliases
using Real = double;
using DvceArray5D = Kokkos::View<Real*****, LayoutWrapper, DevMemSpace>;

// Common operations
DvceArray5D u("array", nmb, nvar, nz, ny, nx);  // Allocate
Kokkos::realloc(u, ...);                        // Resize
Kokkos::deep_copy(dest, src);                   // Copy
auto u_host = Kokkos::create_mirror_view(u);    // Host mirror

// Parallel patterns
par_for("label", DevExeSpace(), 0, n,           // 1D loop
  KOKKOS_LAMBDA(int i) { /*...*/ });

par_for("label", DevExeSpace(),                 // 6D nested
  0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u(m,n,k,j,i) = value;
  });

par_reduce("label", DevExeSpace(), 0, n,        // Reduction
  KOKKOS_LAMBDA(int i, Real &sum) {
    sum += array(i);
  }, total);

Kokkos::fence();  // Synchronize all kernels
```

## Further Reading

- [Kokkos Documentation](https://kokkos.github.io/kokkos-core-wiki/)
- [Kokkos Tutorials](https://github.com/kokkos/kokkos-tutorials)
- [AthenaK Examples](examples/shock_tube.md)
- [Migration Guide](migration/index.md)

## Summary

Kokkos enables AthenaK to run efficiently on both CPUs and GPUs with a single code base. The key concepts are:

1. **Execution/Memory Spaces**: Separate where code runs from where data lives
2. **Views**: Smart multi-dimensional arrays 
3. **Device-Host Communication**: Explicit data movement
4. **Lambda Captures**: Understand what can be safely captured
5. **Parallel Patterns**: Use par_for and par_reduce
6. **MeshBlock Packs**: Process multiple blocks together

Master these concepts to write efficient, portable AthenaK code!