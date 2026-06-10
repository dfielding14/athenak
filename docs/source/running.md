# Running Simulations

All commands assume you are in the repository root and have already built AthenaK in `build/` (adjust the build directory name if you created separate MPI/CUDA builds).

## Basic Execution

```bash
./build/src/athena -i inputs/hydro/sod.athinput
```

You can override parameters at runtime using `block/name=value` assignments:

```bash
./build/src/athena -i inputs/hydro/sod.athinput mesh/nx1=512 time/tlim=1.0
```

## Parallel Execution

### MPI

```bash
mpirun -np 16 ./build-mpi/src/athena -i inputs/hydro/sod.athinput
```

`build-mpi` must be configured with `Athena_ENABLE_MPI=ON` and `Kokkos_ENABLE_MPI=ON`.

### OpenMP

```bash
export OMP_NUM_THREADS=8
./build/src/athena -i inputs/hydro/sod.athinput
```

### Hybrid MPI + OpenMP

```bash
export OMP_NUM_THREADS=4
mpirun -np 4 ./build-mpi/src/athena -i inputs/hydro/sod.athinput
```

## Command-Line Options

| Flag | Meaning |
|------|---------|
| `-i <file>` | Input file (required for new runs) |
| `-r <file>` | Restart from checkpoint |
| `-d <dir>` | Use directory as runtime/output root |
| `-n` | Parse the input and exit (no execution) |
| `-m` | Write mesh structure information and exit |
| `-c` | Show compiled configuration and exit |
| `-t hh:mm:ss` | Wall-clock time limit; triggers a final dump when reached |
| `block/name=value` | Override parameter after parsing the input |

Run `./build/src/athena -h` to see the full usage message.

## Monitoring Progress

AthenaK prints cycle/time information and writes an `eventlog` file (if requested in the input):

```
cycle=0 time=0.0000e+00 dt=1.0000e-03
cycle=100 time=1.0000e-01 dt=1.0234e-03
```

Historical data can also be captured via `<output*>` blocks with `file_type = hst`.

## Output Files

| File type | Description | Typical files |
|-----------|-------------|----------------|
| `vtk` | Structured VTK for ParaView/VisIt | `basename.block*.out*.vtk` |
| `tab` | ASCII tables (one file per output) | `basename.out*.tab` |
| `hst` | History diagnostics | `basename.hst` |
| `bin` | Binary mesh dumps | `basename.block*.out*.bin` |
| `athdf` | HDF5 snapshots (if enabled) | `basename.out*.athdf` |
| `rst` | Restart checkpoints | `basename.NNNNN.rst` |

For parallel I/O, set `single_file_per_rank = true` inside the output block.

## Performance Tips

- **GPU runs**: keep several MeshBlocks per device (8–32 is typical). Control device selection with `CUDA_VISIBLE_DEVICES`/`HIP_VISIBLE_DEVICES`.
- **CPU runs**: enable OpenMP support at build time (`Athena_ENABLE_OPENMP=ON`, `Kokkos_ENABLE_OPENMP=ON`) and set `OMP_PROC_BIND=spread` / `OMP_PLACES=cores` for better affinity.
- **Restart cadence**: configure `<output_restart>` with a modest `dt` or `ncycle` so long jobs can recover after hitting the `-t` wall-clock limit.

## Troubleshooting

| Symptom | Suggested actions |
|---------|-------------------|
| CFL timestep becomes tiny | Lower `time/cfl_number` or refine the mesh to resolve steep gradients |
| Out-of-memory errors | Reduce MeshBlock size (`mesh/nx*`), use more MPI ranks, or decrease output frequency |
| NaNs in solution | Check initial conditions, equation-of-state parameters, and try a more diffusive reconstruction/solver |
| Cannot open restart | Ensure the restart build matches the original configuration (MPI count, physics modules). Use `./build/src/athena -n -r file.rst` to validate |

For additional diagnostics, rerun with `-v 2` (verbose logging) or enable Kokkos profiling (`export KOKKOS_PROFILE_LIBRARY=...`).
