# Kokkos Guide

AthenaK uses Kokkos for portable device execution and storage. This page
documents the public helper types and loop conventions implemented in
`src/athena.hpp`; inspect nearby production kernels before adding a new one.

## Public Types

The project defines:

```cpp
using DevExeSpace = Kokkos::DefaultExecutionSpace;
using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;

template <typename T>
using DvceArray5D = Kokkos::View<T *****, LayoutWrapper, DevMemSpace>;
template <typename T>
using HostArray5D = Kokkos::View<T *****, LayoutWrapper, HostMemSpace>;
template <typename T>
using DualArray1D = Kokkos::DualView<T *, LayoutWrapper, DevMemSpace>;
```

`LayoutRight` makes the last array index contiguous. Fluid field arrays
normally place the MeshBlock index first and cell index `i` last, for example
`u(m, IDN, k, j, i)`.

## Device And Host Data

An array in `DvceArray*` storage must be used from device-capable kernels.
Prepare or inspect its data on the host through a host view and an explicit
copy:

```cpp
DvceArray1D<Real> values("values", n);
auto host_values = Kokkos::create_mirror_view(values);
for (int i = 0; i < n; ++i) {
  host_values(i) = 0.0;
}
Kokkos::deep_copy(values, host_values);
```

For state intentionally modified on both sides of the device boundary,
AthenaK also defines `DualArray*` templates. Follow the `modify_*` and
`sync_*` discipline of the owning module before reading a dual view from the
other memory space.

## AthenaK Loop Helpers

`src/athena.hpp` implements `par_for` overloads for one through five logical
indices. Each helper flattens its inclusive ranges into a Kokkos
`RangePolicy`, then supplies reconstructed indices to a device lambda. A
typical pack-wide kernel therefore has this shape:

```cpp
auto &u0 = pmbp->phydro->u0;
auto &indcs = pmbp->pmesh->mb_indcs;
int nmb = pmbp->nmb_thispack;

par_for("initialize", DevExeSpace(), 0, nmb - 1,
        indcs.ks, indcs.ke, indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m, IDN, k, j, i) = 1.0;
  });
```

Use values and Kokkos views that are device-accessible inside
`KOKKOS_LAMBDA`; do not dereference ordinary host pointers or access standard
library containers there.

AthenaK also provides `par_for_outer` and `par_for_inner` for team/vector
kernels with scratch storage. Reuse an established production pattern when a
new kernel needs that layout.

## Porting And Validation Rules

1. Get block and index metadata from the `MeshBlockPack` / owning `Mesh`, as
   existing problem generators do.
2. Add the MeshBlock index to field accesses and keep `i` as the
   fastest-varying cell index.
3. Avoid host-device copies in the timestep path unless the algorithm requires
   them and the cost is measured.
4. Compile the relevant problem selection and execute a short input deck on
   the enabled backend before treating the port as complete.

See [Migration](migration/index.md), [Problem Generators](modules/pgen.md),
and [Building](building.md) for the corresponding public workflow.
