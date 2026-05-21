# Lagrangian Monte-Carlo Thermodynamic Tracers

`particle_type = lagrangian_mc` adds Monte-Carlo tracer particles that jump between
cell centers with probabilities set by the RK-accumulated density flux through the
faces of the current cell. The tracer state is sampled from the cell-centered fluid
state after the fluid update, so the path-history output follows each particle's
thermodynamic history at the requested output cadence.

This implementation is intended for ideal-gas Hydro and MHD runs. It preserves the
existing drift-particle paths and rejects non-ideal or isothermal thermo-history
requests with a fatal input error.

## Runtime Blocks

Enable the particle type with:

```ini
<particles>
particle_type = lagrangian_mc
pusher        = lagrangian_mc
random_seed   = 12345
```

Tracer seeding is controlled by one or more `<tracer_seedN>` blocks. Each block defines
a schedule, a spatial region, an optional thermodynamic mask, and a deterministic
sampling seed:

```ini
<tracer_seed1>
id              = 1
start_time      = 0.0
end_time        = 0.1
cadence         = 0.01      # <= 0 means one-shot
count_per_event = 1000
weight          = mass      # mass or volume
region          = box       # all, box, sphere, slab
x1min           = -0.5
x1max           =  0.5
x2min           = -0.5
x2max           =  0.5
target          = temperature
target_min      = 1.0e-2
target_max      = 1.0e2
seed            = 24680
```

Schedule fields:

| Field | Meaning |
| --- | --- |
| `id` | Integer stored in particle history output as `seed_id`. |
| `start_time` | First seed event time. Events due at the initial time fire after pgen setup. |
| `end_time` | Last allowed seed event time. |
| `cadence` | Time between seed events. `cadence <= 0` means one-shot. |
| `count_per_event` | Global number of particles requested per event. |
| `weight` | `mass` samples cells in proportion to `rho*dV`; `volume` samples by `dV`. |
| `seed` | Deterministic schedule seed. Same input and decomposition give reproducible tags. |

Region fields:

| Region | Required or useful fields |
| --- | --- |
| `all` | Uses the full mesh. |
| `box` | `x1min/x1max`, `x2min/x2max`, and `x3min/x3max`. |
| `sphere` | `center1`, `center2`, `center3`, and `radius`. |
| `slab` | `slab_axis`, `slab_min`, and `slab_max`. |

The optional `target` selector may be `density`, `temperature`, `pressure`, `entropy`,
or `scalarN`. Add `target_min` and/or `target_max` to keep only cells within a
thermodynamic range. New particles are inserted at selected cell centers, which matches
the cell-jump Monte-Carlo pusher.

Initial seeding fires after the problem generator and primitive-variable initialization,
before initial outputs. Timed seeding fires at the end of `after_timeintegrator`, after
existing particles move and communicate.

## Thermodynamic History Output

Use `file_type = prtcl_thermo_history` for append-only path histories:

```ini
<output1>
file_type = prtcl_thermo_history
dt        = 0.01
```

Each record contains:

```text
time, cycle, tag, seed_id, x1, x2, x3, gid,
rho, pressure, temperature, specific_entropy, internal_energy,
v1, v2, v3, scalar0, scalar1, ...
```

The file is binary and append-only. Read it directly with:

```bash
python scripts/read_prtcl_thermo_history.py \
  prtcl_thermo_history/<basename>.prtcl_thermo_history.thp \
  --npz tracers.npz
```

## Restarts, MPI, and AMR

Restart files persist particle real arrays, integer arrays, `next_tracer_tag`, and each
schedule's `next_time`, `event_index`, and completion flag. Restarted runs therefore
continue tags monotonically and do not duplicate seed events that fired before the
restart.

For MPI restarts, write restart files with:

```ini
<output2>
file_type            = rst
dt                   = 0.02
single_file_per_rank = true
```

The restart command may point at rank 0's file; the restart reader maps the per-rank
files internally. After AMR or load balancing changes the MeshBlock layout, tracers are
remapped to the owning rank and snapped back to valid cell centers on the new mesh.

## Example Inputs

| Input | Purpose |
| --- | --- |
| `inputs/particles/lagrangian_mc_thermo.athinput` | Serial Hydro tracer seeding, timed box schedule, history output, and restart append. |
| `inputs/particles/lagrangian_mc_thermo_amr.athinput` | Hydro AMR tracer run with an initial all-mesh seed, a timed spherical temperature-masked seed, restart output, and MPI-safe restart files. |
| `inputs/particles/lagrangian_mc_thermo_mhd.athinput` | MHD linear-wave smoke test with `lagrangian_mc` tracers and thermo-history output. |

## Validation Run on 2026-05-20

The branch was validated on macOS with Open MPI 5.0.9. The release MPI executable
reported:

```text
Problem generator:          built_in_pgens
Floating-point precision:   double
MPI parallelism:            ON
OpenMP parallelism:         OFF
```

Build and syntax checks:

| Check | Command | Result |
| --- | --- | --- |
| Serial release build | `cmake --build build -j 8` | Passed. |
| MPI release build | `cmake --build build_mpi -j 8` | Passed. |
| ASAN build | `cmake --build build_asan -j 8` | Passed. |
| Python reader syntax | `python3 -m py_compile scripts/read_prtcl_thermo_history.py` | Passed. |
| Whitespace check | `git diff --check` | Passed. |

Runtime checks:

| Check | Command | Key result |
| --- | --- | --- |
| Serial Hydro plus restart append | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo.athinput -d run_serial_hydro` then restart from `lagrangian_mc_thermo.00001.rst` to `time/tlim=0.03` | 104 history records, 32 unique tags, tag range `0..31`, finite data, no duplicate restarted seed event, and seed 2 remained in the requested box. |
| Serial MHD | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo_mhd.athinput -d run_serial_mhd` | 48 history records, 16 unique tags, tag range `0..15`, finite data, stopped at runtime `time=0.02`, cycle 4. |
| Serial Hydro AMR | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo_amr.athinput -d run_serial_amr` | 10 final MeshBlocks, 6 AMR-created MeshBlocks, 136 records, 40 unique tags, tag range `0..39`, finite data, seed 2 satisfied `radius <= 0.3` and `temperature >= 0.1`. |
| Two-rank MPI no-particle AMR smoke | `mpirun -np 2 ./build_mpi/src/athena -i inputs/tests/linear_wave_hydro_amr.athinput -d run_mpi_noparticles mesh/nx1=32 mesh/nx2=32 meshblock/nx1=16 meshblock/nx2=16 time/nlim=1 time/tlim=0.001 output1/dt=1.0 output2/dt=1.0 output3/dt=1.0` | Passed with 4 MeshBlocks, no particle updates, and no AMR load-balancing errors. |
| Two-rank MPI Hydro AMR plus restart | `mpirun -np 2 ./build_mpi/src/athena -i inputs/particles/lagrangian_mc_thermo_amr.athinput -d run_mpi_amr` then restart from rank 0's `lagrangian_mc_thermo_amr.00001.rst` to `time/tlim=0.04` | Forward run ended with 10 MeshBlocks and 40 unique tags. Restarted history had 216 records, still 40 unique tags with tag range `0..39`, and no duplicated timed seed events. |
| Four-rank MPI Hydro AMR | `mpirun -np 4 ./build_mpi/src/athena -i inputs/particles/lagrangian_mc_thermo_amr.athinput -d run_mpi_amr4` | Passed with 10 final MeshBlocks, 6 AMR-created MeshBlocks, 8 MeshBlocks communicated for load balancing, 136 records, 40 unique tags, finite data, and valid seed 2 mask checks. |
| ASAN Hydro AMR | `ASAN_OPTIONS=abort_on_error=0:detect_leaks=0 ./build_asan/src/athena -i inputs/particles/lagrangian_mc_thermo_amr.athinput -d run_asan_amr` | Passed with 10 final MeshBlocks, 136 records, 40 unique tags, and no AddressSanitizer report. |

History record counts by output time:

| Case | Counts |
| --- | --- |
| Serial Hydro restart | `0: 16`, `0.0124834403634: 24`, `0.02: 32`, `0.03: 32`. |
| Serial MHD | `0: 16`, `0.0124918535307: 16`, `0.02: 16`. |
| AMR forward runs | `0: 24`, `0.0124834403634: 32`, `0.02496688271: 40`, `0.03: 40`. |
| Two-rank AMR restart | `0: 24`, `0.0124834403634: 32`, `0.02496688271: 40`, `0.03: 40`, `0.0312086052974: 40`, `0.04: 40`. |

During four-rank testing, a real MPI dynamic-seeding bug was found and fixed: ranks that
received zero newly seeded particles were skipping the collective particle-count update.
The count update now runs collectively after every seeding pass, and the four-rank AMR
case above exercises that fix.
