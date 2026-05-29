# Migration Gotchas

## Pack-Wide Arrays

AthenaK field storage normally includes a MeshBlockPack index:

```cpp
u0(m, IDN, k, j, i)
```

Use `pmy_mesh_->pmb_pack` and `pmy_mesh_->mb_indcs` in a problem generator,
following an implemented generator in `src/pgen/`. Do not mechanically retain
single-block Athena++ indexing.

## Device Accessibility

`DvceArray*` data is device storage. Fill or inspect it on the host by using a
host mirror plus `Kokkos::deep_copy`, or write it in a device lambda:

```cpp
DvceArray1D<Real> data("data", n);
auto host = Kokkos::create_mirror_view(data);
host(0) = 1.0;
Kokkos::deep_copy(data, host);
```

Ordinary host pointers and standard-library containers must not be
dereferenced in a `KOKKOS_LAMBDA` intended for GPU execution.

## Input Compatibility

- A custom `src/pgen/name.cpp` generally requires a build configured with
  `-DPROBLEM=name`; only generators registered in the built-in dispatcher can
  be selected with `<problem>/pgen_name`.
- SR and GR fluid configurations reject `eos = isothermal`.
- Refined cases that require at least three ghost zones must choose an even
  value, normally `nghost >= 4`.
- The public particle module accepts `particle_type = cosmic_ray` and
  `pusher = drift`; branch records describing Boris or star particles are not
  public runtime configuration.

## Backend Builds

CUDA builds use the bundled wrapper:

```bash
cmake -S . -B build-cuda -DKokkos_ENABLE_CUDA=ON \
  -DCMAKE_CXX_COMPILER=${PWD}/kokkos/bin/nvcc_wrapper
```

HIP builds must select a ROCm compiler:

```bash
export ROCM_PATH=/opt/rocm
cmake -S . -B build-hip -DKokkos_ENABLE_HIP=ON \
  -DCMAKE_CXX_COMPILER=${ROCM_PATH}/bin/hipcc
```

Add an appropriate `Kokkos_ARCH_*` setting for the actual device and validate
with a short shipped input before assessing performance.

See [Migration](index.md), [Building](../building.md), and
[Troubleshooting](../troubleshooting.md).
