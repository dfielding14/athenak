# Troubleshooting Guide

This page collects quick diagnostics for common build and runtime issues. Start with the section that matches your symptoms.

## Build Problems

### "Cannot find Kokkos" / missing headers
```bash
git submodule update --init --recursive
```
Make sure you configure from the repository root (`cmake -S . -B build`). Delete the build directory and reconfigure if the cache still references an old path.

### CUDA / HIP toolchain not detected
```bash
export CUDAToolkit_ROOT=/usr/local/cuda   # adjust to installation path
export PATH="$CUDAToolkit_ROOT/bin:$PATH"
export ROCM_PATH=/opt/rocm
```
For CUDA, re-run CMake with `-DKokkos_ENABLE_CUDA=ON`,
`-DCMAKE_CXX_COMPILER=${PWD}/kokkos/bin/nvcc_wrapper`, and the matching
`Kokkos_ARCH_*` option. For HIP, use `-DKokkos_ENABLE_HIP=ON`,
`-DCMAKE_CXX_COMPILER=${ROCM_PATH}/bin/hipcc`, and the matching AMD
`Kokkos_ARCH_*` option.

### MPI compiler wrappers unavailable
```bash
cmake -S . -B build-mpi -DAthena_ENABLE_MPI=ON -DMPI_CXX_COMPILER=mpicxx
```
Ensure your MPI module is loaded and CMake sees the same compiler that Kokkos was built with.

## Runtime Diagnostics

### Inspect parsed input before running
```bash
./build/src/athena -i input.athinput -n
```
This parses the file, applies overrides, and dumps the resolved parameters without
starting the simulation. It exits before mesh and physics objects are constructed,
so a short actual run is still required to test compatibility checks.

### Crash or segmentation fault
- Rebuild with debug information (`cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug`).
- Run under `gdb --args ./build-debug/src/athena -i input.athinput` to obtain a backtrace.
- Check that all problem-generator arrays respect mesh/ghost extents.

### Simulation stops because `dt` becomes tiny
- Inspect the solution immediately before the timestep collapse for nonphysical
  states or rapidly growing velocities/wave speeds.
- Reduce to a shipped problem or a smaller reproduction before changing
  resolution or numerical methods.
- Inspect history output (`*.hst`) for unstable trends in global conserved
  quantities; inspect field output to diagnose spatial gradients.

### NaNs or negative density/pressure
- Verify `<hydro>` / `<mhd>` parameters (EOS, floors) match the problem.
- Try a more diffusive reconstruction (`reconstruct = plm`) or solver (`rsolver = hlle`).
- To test first-order flux correction (FOFC), add `fofc = true` to the active
  `<hydro>` or `<mhd>` block and set `<mesh>/nghost >= 3` when using `plm`;
  higher-order reconstruction with FOFC requires `nghost >= 4`.
- SMR/AMR requires an even `nghost`, so a refined case that otherwise needs
  at least three ghost zones must use `nghost >= 4`.

### MHD `div(B)` growth
- Confirm `<mhd>` is enabled (constrained transport is always active for MHD runs).
- To test FOFC, add `fofc = true` under `<mhd>` and satisfy the ghost-zone
  requirement above; FOFC is a run-wide MHD setting, not an AMR-level switch.
- Reduce the case to a shipped MHD input and record whether mesh refinement is
  involved before reporting a divergence-control issue.

### Parallel job hangs or underperforms
- Verify each executable was built with matching MPI/OpenMP settings.
- Use `mpirun --bind-to core --map-by socket` (OpenMPI) or equivalent binding flags.
- Enable Kokkos profiling to spot load imbalance (adjust the library path to your tool installation):
  ```bash
  export KOKKOS_PROFILE_LIBRARY=/path/to/libkokkos-tools.so
  ./build/src/athena -i input.athinput
  ```

## Output & Restart Tips

- History tables (`file_type = hst`) are helpful for tracking conserved quantities.
- Restart dumps (`file_type = rst`) are written under `rst/`; resume with
  `./build/src/athena -r rst/basename.00001.rst` or the latest later
  checkpoint. The `00000` checkpoint is generally the initialization state.
- The `-t hh:mm:ss` command-line flag sets AthenaK's wall-clock limit. Set it
  below the scheduler limit and combine it with a restart output stream so
  AthenaK has time to finalize and write the last checkpoint.

## When to Ask for Help

- Search the documentation (`make -C docs html` then open `docs/build/html/search.html`).
- Reproduce the problem with a bundled example if possible.
- Open an issue on GitHub with your input file, run command, compiler/toolchain versions, and relevant logs.

Related resources: [Configuration](configuration.md), [Common Gotchas](migration/common_gotchas.md), [Migration Guide](migration/index.md).
