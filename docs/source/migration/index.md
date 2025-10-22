# Athena++ to AthenaK Migration Guide

## Overview

AthenaK is a complete rewrite of Athena++ using Kokkos for performance portability. While the physics and algorithms remain similar, the implementation details have changed significantly to support GPU execution. This guide helps you port existing Athena++ code to AthenaK.

## Quick Reference

| Concept | Athena++ | AthenaK |
|---------|----------|---------|
| **Arrays** | `AthenaArray<Real>` | `DvceArray5D<Real>` (Kokkos::View) |
| **Loops** | Nested for loops | `par_for` with lambdas |
| **Memory** | CPU only | Host/Device spaces |
| **MeshBlock Access** | Direct member access | Through MeshBlockPack |
| **Parallelism** | OpenMP threads | Kokkos execution spaces |
| **Boundary Conditions** | Function pointers | Task-based |
| **Problem Generators** | Similar structure | Add MeshBlock dimension |

## Major Architectural Changes

### 1. Kokkos-Based Parallelism
- All computational loops use Kokkos parallel patterns
- Memory managed through Kokkos Views
- Explicit host/device memory spaces

### 2. MeshBlockPack Design
- Multiple MeshBlocks grouped for GPU efficiency
- Extra dimension in all field arrays
- Enables coalesced memory access

### 3. Task-Based Execution
- All operations expressed as tasks with dependencies
- Automatic overlap of computation and communication
- Fine-grained parallelism

## Migration Workflow

The typical migration workflow from Athena++ to AthenaK:

1. **Identify Loops** - Find all computational loops
2. **Convert Arrays** - Change to Kokkos Views
3. **Update Access Patterns** - Add MeshBlock dimension
4. **Add par_for** - Replace loops with Kokkos patterns
5. **Handle Device/Host** - Manage memory spaces
6. **Test & Validate** - Verify correctness

## Common Patterns

### Loop Conversion

<div class="migration-comparison">

<div>

**Athena++**
```cpp
for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      u(IDN,k,j,i) = rho;
      u(IM1,k,j,i) = rho * v1;
    }
  }
}
```

</div>

<div>

**AthenaK**
```cpp
par_for("InitDensity", DevExeSpace(),
  0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u(m,IDN,k,j,i) = rho;
    u(m,IM1,k,j,i) = rho * v1;
  });
```

</div>

</div>

### Array Declaration

<div class="migration-comparison">

<div>

**Athena++**
```cpp
AthenaArray<Real> u;
u.NewAthenaArray(NVAR,nk,nj,ni);
u(IDN,k,j,i) = value;
```

</div>

<div>

**AthenaK**
```cpp
DvceArray5D<Real> u;
Kokkos::realloc(u, nmb,NVAR,nk,nj,ni);
u(m,IDN,k,j,i) = value;
```

</div>

</div>

### Problem Generator Structure

<div class="migration-comparison">

<div>

**Athena++**
```cpp
void MeshBlock::ProblemGenerator(
  ParameterInput *pin) {
  // Direct MeshBlock access
  auto &u = pmb->phydro->u;
  
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        // Initialize single MeshBlock
        prim(IDN,k,j,i) = 1.0;
      }
    }
  }
}
```

</div>

<div>

**AthenaK**
```cpp
void ProblemGenerator::UserProblem(
  ParameterInput *pin, const bool restart) {
  // Access through MeshBlockPack
  auto &mbp = pmy_mesh_->pmb_pack;
  
  par_for("pgen", DevExeSpace(),
    0, mbp->nmb_thispack-1,
    mbp->is, mbp->ie,
    mbp->js, mbp->je,
    mbp->ks, mbp->ke,
    KOKKOS_LAMBDA(int m, int i, int j, int k) {
      // Initialize all MeshBlocks in pack
      prim(IDN,m,k,j,i) = 1.0;
    });
}
```

</div>

</div>

### Task System Structure

<div class="migration-comparison">

<div>

**Athena++**
```cpp
TaskStatus MyTask(MeshBlock *pmb, 
                  int stage) {
  // Task operates on single MeshBlock
  auto &u = pmb->phydro->u;
  
  // Do work on this MeshBlock
  for (int k=ks; k<=ke; k++) {
    // Process cells
  }
  
  return TaskStatus::complete;
}
```

</div>

<div>

**AthenaK**
```cpp
TaskStatus MyTask(Driver *pdrive, 
                  int stage) {
  // Task operates on all MeshBlocks
  auto &mbp = pdrive->pmesh->pmb_pack;
  
  par_for("task", DevExeSpace(),
    0, mbp->nmb_thispack-1,
    mbp->ks, mbp->ke,
    KOKKOS_LAMBDA(int m, int k) {
      // Process all blocks in parallel
    });
  
  return TaskStatus::complete;
}
```

</div>

</div>

## Key Concepts to Master

### 1. **Kokkos Views** 
Replace all arrays with Kokkos Views for device compatibility

### 2. **Lambda Captures**
Understand what can be safely captured in device lambdas

### 3. **Memory Spaces**
Manage data movement between host and device

### 4. **Execution Policies**
Choose appropriate execution spaces (Serial, OpenMP, CUDA, etc.)

### 5. **MeshBlock Packing**
Work with multiple blocks simultaneously

## Input File Changes

Most Athena++ input files work with minimal modifications:

### Parameters That Remain The Same
```ini
<time>
tlim = 1.0
integrator = rk2    # Same integrator options

<mesh>
nx1 = 256           # Same resolution parameters
x1min = 0.0         # Same domain bounds
x1max = 1.0         

<hydro>             # Same physics modules
gamma = 1.4
reconstruct = plm
rsolver = hllc
```

### New Capabilities
```ini
<particles>         # New: Lagrangian particles
nspecies = 1
pusher = boris

<z4c>               # New: Numerical relativity
use_z4c = true

<dyn_grmhd>         # New: Dynamical GRMHD
use_dyngr = true
```

### Changed Parameters
```ini
# Athena++
<mesh>
refinement = adaptive

# AthenaK
<mesh>
refinement = adaptive
max_level = 3       # Must specify max level
```

## Additional Resources

- [**Common Gotchas**](common_gotchas.md) - Pitfalls and solutions
- [**Kokkos Concepts**](../kokkos_guide.md) - Understanding Kokkos for AthenaK

## Performance Tips

### GPU Optimization
- **Coalesced Access**: Keep innermost loop over i for contiguous memory
- **Minimize Divergence**: Avoid complex conditionals in kernels
- **Batch Operations**: Process multiple MeshBlocks together
- **Avoid Atomics**: Design algorithms to avoid atomic operations

### Memory Management
- **Minimize Transfers**: Keep data on device as long as possible
- **Use Scratch Space**: Leverage Kokkos scratch memory for temporaries
- **Explicit Sync**: Add fence() operations where needed
- **View Reuse**: Reuse allocated Views rather than reallocating

## Common Pitfalls and Solutions

### ❌ Problem: Capturing Host Pointers
```cpp
// WRONG - Will crash on GPU
Real *host_data = new Real[n];
par_for(..., KOKKOS_LAMBDA(int i) {
  value = host_data[i];  // Error!
});
```

### ✅ Solution: Use Device Views
```cpp
// CORRECT - Device-compatible
DvceArray1D<Real> device_data("data", n);
par_for(..., KOKKOS_LAMBDA(int i) {
  value = device_data(i);  // Safe
});
```

### ❌ Problem: Missing MeshBlock Index
```cpp
// WRONG - Missing 'm' dimension
u(IDN,k,j,i) = density;
```

### ✅ Solution: Include All Dimensions
```cpp
// CORRECT - Include MeshBlock index
u(m,IDN,k,j,i) = density;
```

### ❌ Problem: Wrong Execution Space
```cpp
// WRONG - Using host execution for device data
par_for("Loop", HostSpace(), ...,
  KOKKOS_LAMBDA(int i) {
    device_array(i) = value;  // May fail
  });
```

### ✅ Solution: Match Execution and Memory Spaces
```cpp
// CORRECT - Device execution for device data
par_for("Loop", DevExeSpace(), ...,
  KOKKOS_LAMBDA(int i) {
    device_array(i) = value;  // Correct
  });
```

## Validation Checklist

Before considering your port complete:

- [ ] All loops converted to Kokkos patterns
- [ ] No raw pointers in device kernels
- [ ] Memory spaces explicitly managed
- [ ] MeshBlock pack dimension added
- [ ] Boundary conditions updated
- [ ] No race conditions in parallel regions
- [ ] Reductions use Kokkos patterns
- [ ] Code tested on both CPU and GPU
- [ ] Performance comparable or better
- [ ] Results match Athena++ baseline

## Getting Help

### Resources
- [Kokkos Documentation](https://kokkos.github.io/kokkos-core-wiki/)
- [AthenaK Examples](../examples/shock_tube.md)
- [API Reference](../reference/api_reference.md)

### Support
- GitHub Issues: Report bugs or ask questions
- Wiki: Additional tutorials and guides
- Example Problems: Study existing implementations

## Next Steps

1. Review [Common Gotchas](common_gotchas.md) to avoid common pitfalls
2. Study the [Kokkos Guide](../kokkos_guide.md) to understand device programming
3. Try migrating a simple problem first (e.g., shock tube)
4. Test on CPU before enabling GPU
5. Profile and optimize for your target architecture