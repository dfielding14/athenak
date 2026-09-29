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
export CUDA_HOME=/usr/local/cuda          # adjust to installation path
export PATH="$CUDA_HOME/bin:$PATH"
export HIP_PATH=/opt/rocm
```
Re-run CMake with the appropriate flags (`-DKokkos_ENABLE_CUDA=ON` or `-DKokkos_ENABLE_HIP=ON`) and set the matching `Kokkos_ARCH_*` option.

### MPI compiler wrappers unavailable
```bash
cmake -S . -B build-mpi   -DAthena_ENABLE_MPI=ON   -DKokkos_ENABLE_MPI=ON   -DMPI_CXX_COMPILER=mpicxx
```
Ensure your MPI module is loaded and CMake sees the same compiler that Kokkos was built with.

## Runtime Diagnostics

### Validate the input before running
```bash
./build/src/athena -i input.athinput -n
```
This parses the file, applies overrides, and prints any warnings without starting the simulation.

### Crash or segmentation fault
- Rebuild with debug information (`cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug`).
- Run under `gdb --args ./build-debug/src/athena -i input.athinput` to obtain a backtrace.
- Check that all problem-generator arrays respect mesh/ghost extents.

### Simulation stops because `dt` becomes tiny
- Lower `time/cfl_number` to improve stability.
- Increase resolution around discontinuities (refine or reduce MeshBlock size).
- Inspect the history output (`*.hst`) for rapidly growing gradients.

### NaNs or negative density/pressure
- Verify `<hydro>` / `<mhd>` parameters (EOS, floors) match the problem.
- Try a more diffusive reconstruction (`reconstruct = plm`) or solver (`rsolver = hlle`).
- Enable first-order flux correction (FOFC) in problematic regions.

### MHD `div(B)` growth
- Confirm `<mhd>` is enabled (constrained transport is always active for MHD runs).
- Enable FOFC (`hydro/fofc = true` or `mhd/fofc = true`) on troublesome AMR levels.
- Check prolongation settings in `<mesh_refinement>` to ensure divergence-free interpolation.

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
- Restart dumps (`file_type = rst`) allow `./build/src/athena -r basename.NNNNN.rst` to resume runs; ensure the restart build was configured with the same modules.
- The `-t hh:mm:ss` command-line flag sets a wall-clock limit—combine it with a restart output stream to survive queue timeouts.

## When to Ask for Help

- Search the documentation (`make html` → `docs/build/html/search.html`).
- Reproduce the problem with a bundled example if possible.
- Open an issue on GitHub with your input file, run command, compiler/toolchain versions, and relevant logs.

Related resources: [Configuration](configuration.md), [Common Gotchas](migration/common_gotchas.md), [Migration Guide](migration/index.md).
