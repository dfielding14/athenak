# Migration Pattern: [Pattern Name]

## Overview
[Brief description of what aspect of code is being migrated]

## Key Differences

| Aspect | Athena++ | AthenaK |
|--------|----------|---------|
| [Aspect 1] | [Old approach] | [New approach] |
| [Aspect 2] | [Old approach] | [New approach] |
| Memory Model | CPU-centric | Device-aware (Kokkos) |
| Parallelism | OpenMP/MPI | Kokkos/MPI |

## Side-by-Side Comparison

<div class="migration-comparison">

<div>

### Athena++ Implementation
```cpp
// Athena++ style code
for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      // Direct array access
      cons(IDN,k,j,i) = value;
      cons(IM1,k,j,i) = momentum;
    }
  }
}
```

**Key Points:**
- Explicit nested loops
- Direct array indexing
- CPU-only execution
- Sequential memory access

</div>

<div>

### AthenaK Implementation
```cpp
// AthenaK style code
auto &indcs = pmbp->pmesh->mb_indcs;
int &is = indcs.is, &js = indcs.js, &ks = indcs.ks;
int &ie = indcs.ie, &je = indcs.je, &ke = indcs.ke;

par_for("UpdateCons", DevExeSpace(), 
    0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, 
                  const int j, const int i) {
      // Kokkos View access
      cons(m,IDN,k,j,i) = value;
      cons(m,IM1,k,j,i) = momentum;
    });
```

**Key Points:**
- Single `par_for` call
- Lambda function for kernel
- Device-compatible execution
- MeshBlock pack dimension

</div>

</div>

## Step-by-Step Migration

### 1. Identify Loop Structure
```cpp
// Original Athena++ pattern
for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      // Loop body
    }
  }
}
```

### 2. Convert to par_for
```cpp
// Step 1: Replace loops with par_for
par_for("KernelName", DevExeSpace(), 
    ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int k, int j, int i) {
      // Loop body (modified for device)
    });
```

### 3. Handle MeshBlock Packs
```cpp
// Step 2: Add MeshBlock dimension if needed
par_for("KernelName", DevExeSpace(),
    0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      // Access with MeshBlock index m
    });
```

### 4. Update Array Access
```cpp
// Athena++ array access
AthenaArray<Real> &u = phydro->u;
u(IDN,k,j,i) = value;

// AthenaK array access
auto u = pmbp->phydro->u;  // Kokkos View
u(m,IDN,k,j,i) = value;    // Note extra dimension
```

## Data Structure Changes

### Arrays and Views

| Athena++ | AthenaK | Notes |
|----------|---------|-------|
| `AthenaArray<Real>` | `DvceArray5D<Real>` | Device array for conserved vars |
| `Real***` | `DvceArray3D<Real>` | 3D device array |
| `new Real[n]` | `Kokkos::View<Real*>` | Dynamic allocation |

### Memory Spaces

```cpp
// Athena++: CPU memory only
Real *data = new Real[size];

// AthenaK: Explicit memory space
DualArray1D<Real> data("label", size);
auto d_data = data.d_view;  // Device view
auto h_data = data.h_view;  // Host view
```

## Common Pitfalls

### ❌ Incorrect: Capturing Host Data
```cpp
// WRONG: Captures host pointer in device lambda
Real *host_ptr = some_host_data;
par_for("Bad", DevExeSpace(), 0, n,
    KOKKOS_LAMBDA(int i) {
      value = host_ptr[i];  // Will crash on GPU!
    });
```

### ✅ Correct: Use Device Views
```cpp
// RIGHT: Use device-compatible view
auto device_view = some_device_data;
par_for("Good", DevExeSpace(), 0, n,
    KOKKOS_LAMBDA(int i) {
      value = device_view(i);  // Safe on GPU
    });
```

### ❌ Incorrect: Wrong Index Order
```cpp
// WRONG: Forgetting MeshBlock index
u(IDN,k,j,i) = value;  // Missing 'm' dimension
```

### ✅ Correct: Include All Dimensions
```cpp
// RIGHT: Include MeshBlock pack index
u(m,IDN,k,j,i) = value;  // Correct for MeshBlockPack
```

## Performance Tips

1. **Coalesced Access**: Ensure innermost loop accesses contiguous memory
2. **Avoid Atomics**: Design algorithms to avoid atomic operations when possible
3. **Minimize Divergence**: Keep conditional logic simple in kernels
4. **Use Scratch Memory**: Leverage Kokkos scratch spaces for temporary data

## Verification Checklist

- [ ] All loops converted to `par_for` or equivalent
- [ ] Array accesses updated for new dimensions
- [ ] No host pointers captured in device lambdas
- [ ] Memory spaces explicitly managed
- [ ] Synchronization points added where needed
- [ ] Reductions use Kokkos reduction patterns
- [ ] Boundary conditions updated for new structure

## Related Patterns

- [Array Management](arrays.md)
- [Loop Parallelization](loops.md)
- [Memory Spaces](memory.md)
- [Reductions](reductions.md)

## Examples in Codebase

- Simple loop: `src/hydro/hydro_update.cpp`
- Complex kernel: `src/mhd/mhd_fluxes.cpp`
- Reduction: `src/hydro/hydro_newdt.cpp`