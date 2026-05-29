# Running Simulations

All commands assume you are in the repository root and have already built AthenaK in `build/` (adjust the build directory name if you created separate MPI/CUDA builds).

## Basic Execution

```bash
./build/src/athena -i inputs/hydro/sod.athinput
```

Use `-d <dir>` to place generated output beneath a run directory. The input or
restart path is opened before AthenaK changes into that directory, so specify
those paths relative to the directory from which you launch the executable:

```bash
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod
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

`build-mpi` must be configured with `Athena_ENABLE_MPI=ON`; top-level CMake
enables Kokkos MPI support as part of that configuration.

### OpenMP

```bash
cmake -S . -B build-omp -DAthena_ENABLE_OPENMP=ON
cmake --build build-omp
export OMP_NUM_THREADS=8
./build-omp/src/athena -i inputs/hydro/sod.athinput
```

### Hybrid MPI + OpenMP

```bash
cmake -S . -B build-mpi-omp -DAthena_ENABLE_MPI=ON -DAthena_ENABLE_OPENMP=ON
cmake --build build-mpi-omp
export OMP_NUM_THREADS=4
mpirun -np 4 ./build-mpi-omp/src/athena -i inputs/hydro/sod.athinput
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

AthenaK prints cycle/time information during a run. An output block with
`file_type = log` writes event-counter diagnostics to `basename.log` when
events such as floor applications or FOFC activity occur:

```
elapsed=0.0000e+00 cycle=0 time=0.0000e+00 dt=1.0000e-03
elapsed=1.0000e+00 cycle=100 time=1.0000e-01 dt=1.0234e-03
```

Historical data can also be captured via `<output*>` blocks with `file_type = hst`.

## Output Files

| File type | Description | Typical files |
|-----------|-------------|----------------|
| `vtk` | Legacy VTK mesh data for ParaView/VisIt | `vtk/basename.id.00000.vtk` |
| `tab` | ASCII sliced tables | `tab/basename.id.00000.tab` |
| `hst` | History diagnostics, one file per contributing physics module | `basename.hydro.hst`, `basename.mhd.hst` |
| `bin` | AthenaK binary mesh dumps | `bin/basename.id.00000.bin` |
| `rst` | Restart checkpoints | `rst/basename.00000.rst` |

Additional registered output types are `log`, `pvtk`, `trk`, `cbin`, `pdf`,
`cart`, and `sph`; see [Outputs](modules/outputs.md). For `bin`, `cbin`, or
`rst`, `single_file_per_rank = true` selects per-rank files.

## Performance Tips

- **GPU runs**: control device selection with `CUDA_VISIBLE_DEVICES` or
  `HIP_VISIBLE_DEVICES`; tune MeshBlock size using measurements for the target problem.
- **CPU runs**: enable OpenMP support at build time with
  `Athena_ENABLE_OPENMP=ON`; top-level CMake enables Kokkos OpenMP support.
- **Restart cadence**: configure one numbered output block with
  `file_type = rst` and a suitable `dt` or `dcycle` so a wall-clock-limited job
  writes resumable state.

## Troubleshooting

| Symptom | Suggested actions |
|---------|-------------------|
| CFL timestep becomes tiny | Inspect the solution and diagnostics for nonphysical states or rapidly growing signal speeds; reduce the case before changing numerical parameters |
| Out-of-memory errors | Reduce global cell counts (`mesh/nx*`) or use more MPI ranks; separately measure the effect of MeshBlock decomposition (`meshblock/nx*`), since smaller blocks add ghost-zone overhead; check whether failure coincides with field-output allocation |
| NaNs in solution | Check initial conditions, equation-of-state parameters, and try a more diffusive reconstruction/solver |
| Cannot open restart | Verify the checkpoint path under `rst/` and use `./build/src/athena -n -r rst/file.rst` to parse its saved parameters |

For additional diagnostics, inspect `./build/src/athena -c` and enable an
installed Kokkos Tools library with `export KOKKOS_PROFILE_LIBRARY=/path/to/library`.
