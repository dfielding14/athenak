# Lagrangian Monte-Carlo Thermodynamic Tracers

This source-adjacent note stays next to the runnable tracer inputs. The
GitHub Pages-ready module page is
`docs/source/modules/lagrangian_mc_tracers.md`; use that file when publishing this
feature into the Sphinx/MyST documentation tree on `origin/gh-pages`.

`particle_type = lagrangian_mc` adds Monte-Carlo tracer particles that jump between
cell centers with probabilities set by the RK-accumulated density flux through the
faces of the current cell. Tracers carry only transport state in the particle arrays:
position, tag, owning MeshBlock, seed id, creation time, and AMR bookkeeping. Requested
fluid variables are sampled from the current cell-centered fluid state when history
output is written.

This implementation is intended for ideal-gas Hydro and MHD runs. It preserves the
existing drift-particle paths and rejects non-ideal or isothermal thermo-history
requests with a fatal input error.

## Particle Setup

Enable the particle type with:

```ini
<particles>
particle_type    = lagrangian_mc
pusher           = lagrangian_mc
random_seed      = 12345
track_variables  = default, mach, sound_speed, vmag, scalar0
```

`track_variables` is the default history-output field list. Individual
`prtcl_thermo_history` outputs may override it with their own `variables` list.

`default` expands to:

```text
density, pressure, temperature, entropy, internal_energy, v1, v2, v3
```

## Tracer Fields

The same tracer-field registry is used for `track_variables`, per-output `variables`,
and seed-block `target` masks. Names are case-insensitive, comma or whitespace
separated, and stored in output files using the canonical names below.

| Canonical name | Accepted aliases | Meaning |
| --- | --- | --- |
| `density` | `rho` | Cell density. |
| `pressure` | `p` | Ideal-gas pressure. |
| `temperature` | `T` | `pressure/density`. |
| `entropy` | `s`, `specific_entropy` | `log(pressure/density^gamma)`. |
| `internal_energy` | `eint` | Internal-energy primitive slot. |
| `v1`, `v2`, `v3` | `vx`, `vy`, `vz` | Fluid velocity components. |
| `vmag` | `velocity_magnitude` | Magnitude of the fluid velocity. |
| `sound_speed` | `cs` | Ideal-gas sound speed. |
| `mach` | `mach_number` | `vmag/sound_speed`. |
| `scalarN` | none | Passive scalar index `N`, for example `scalar0`. |

MHD-only fields:

| Canonical name | Accepted aliases | Meaning |
| --- | --- | --- |
| `b1`, `b2`, `b3` | `bx`, `by`, `bz` | Cell-centered magnetic-field components. |
| `bmag` | `magnetic_field_magnitude` | Magnetic-field magnitude. |
| `magnetic_pressure` | `pmag` | `0.5*bmag^2`. |
| `beta` | `plasma_beta` | Gas pressure divided by magnetic pressure. |
| `alfven_speed` | `va` | `bmag/sqrt(density)`. |

Unknown field names, unavailable MHD fields in Hydro runs, and out-of-range passive
scalar indices are fatal input errors.

Problem-specific derived quantities are not parsed from runtime expressions in this
version. For custom diagnostics, store the quantity in a passive scalar and track
`scalarN`. A later extension can add pgen-enrolled tracer-field callbacks for quantities
such as `cooling_time` or `tcool_over_tff`.

## Seeding Schedules

Tracer seeding is controlled by one or more `<tracer_seedN>` blocks. Each block defines
a schedule, a spatial region, an optional field mask, and a deterministic sampling seed:

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
target          = mach
target_min      = 1.0
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

`target`, `target_min`, and `target_max` use the tracer-field registry above. New
particles are inserted at selected cell centers, which matches the cell-jump
Monte-Carlo pusher.

Initial seeding fires after the problem generator and primitive-variable initialization,
before initial outputs. Timed seeding fires at the end of `after_timeintegrator`, after
existing particles move and communicate.

## Thermodynamic History Output

Use `file_type = prtcl_thermo_history` for append-only path histories:

```ini
<output1>
file_type = prtcl_thermo_history
dt        = 0.01
variables = density, temperature, mach, scalar0
```

If `variables` is omitted, the output uses `<particles>/track_variables`. If that is
also omitted, it uses `default`.

Every record contains the fixed metadata fields followed by the requested variables in
the exact order stored in the file header:

```text
time, cycle, tag, seed_id, x1, x2, x3, gid, <requested variables...>
```

The binary `.thp` file starts with a versioned schema header containing the canonical
variable names. Each appended block repeats the record sizes and output time. Appending
is intentionally strict: if a restart or rerun tries to append to an existing `.thp`
file with a different variable list, real precision, or schema version, AthenaK exits
with a clear fatal error rather than silently corrupting the file.

Read history files with:

```bash
python scripts/read_prtcl_thermo_history.py \
  prtcl_thermo_history/<basename>.prtcl_thermo_history.thp \
  --npz tracers.npz
```

The reader returns arrays by column name. For the example above, the keys include
`density`, `temperature`, `mach`, and `scalar0`. The reader also supports the previous
fixed-layout version-1 files for compatibility.

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
| `inputs/particles/lagrangian_mc_thermo.athinput` | Serial Hydro tracer seeding, scalar tracking, derived-Mach seed mask, history output with `variables = density, temperature, mach, scalar0`, and restart append. |
| `inputs/particles/lagrangian_mc_thermo_amr.athinput` | Hydro AMR tracer run using `<particles>/track_variables = density, temperature, mach`, an initial all-mesh seed, a timed spherical Mach-masked seed, restart output, and MPI-safe restart files. |
| `inputs/particles/lagrangian_mc_thermo_mhd.athinput` | MHD linear-wave smoke test with `variables = bmag, beta, alfven_speed, mach`. |

## Validation Run on 2026-05-21

Refresh this section when changing the tracer registry, binary history schema, restart
logic, or particle communication code.

Build and static checks:

| Check | Command | Result |
| --- | --- | --- |
| Serial release build | `cmake --build build -j 8` | Passed. The build reported only pre-existing VLA warnings in `src/mesh/mesh.cpp`. |
| MPI release build | `cmake --build build_mpi -j 8` | Passed with the same pre-existing VLA warnings. |
| Reader unit test | `python3 scripts/test_prtcl_thermo_history_reader.py` | Passed. The synthetic v2 file returned named columns `density`, `temperature`, `mach`, and `scalar0` with 3 records. |
| Reader syntax | `python3 -m py_compile scripts/read_prtcl_thermo_history.py scripts/test_prtcl_thermo_history_reader.py` | Passed. |
| Whitespace check | `git diff --check` | Passed. |

Runtime checks:

| Check | Command | Result |
| --- | --- | --- |
| Serial Hydro custom variables | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo.athinput -d run_validation_serial_hydro` | Passed with columns `density`, `temperature`, `mach`, `scalar0`; 72 records; 32 unique tags; tag range `0..31`; finite data; `scalar0 = 0.25` for all records. |
| Serial Hydro same-schema restart append | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo.athinput -d run_validation_serial_hydro_restart output3/dt=0.01`, then `./build/src/athena -r run_validation_serial_hydro_restart/rst/rank_00000000/lagrangian_mc_thermo.00001.rst -d run_validation_serial_hydro_restart time/tlim=0.03` | Passed. The appended file has 136 records, 32 unique tags, tag range `0..31`, and time counts `0:16`, `0.0124834403634:24`, `0.02:32`, `0.0249668822191:32`, `0.03:32`. |
| Schema mismatch append | Restart into `run_validation_serial_hydro_restart` with `output1/variables=density,temperature` | Failed as intended with `existing particle thermo file schema does not match this output` and exit status 1. |
| MHD derived fields | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo_mhd.athinput -d run_validation_serial_mhd` | Passed with columns `bmag`, `beta`, `alfven_speed`, `mach`; 48 records; 16 unique tags; all requested variables finite. |
| MPI Hydro AMR custom variable list | `mpirun -np 2 ./build_mpi/src/athena -i inputs/particles/lagrangian_mc_thermo_amr.athinput -d run_validation_mpi_amr` | Passed. The run ended with 10 MeshBlocks after creating 6 by AMR; history columns are `density`, `temperature`, `mach`; 136 records; 40 unique tags; tag range `0..39`; all requested variables finite; seed-2 particles satisfy the spherical mask with maximum recorded radius `0.2974911369`. |
| MPI Hydro AMR restart append | `mpirun -np 2 ./build_mpi/src/athena -r run_validation_mpi_amr/rst/rank_00000000/lagrangian_mc_thermo_amr.00001.rst -d run_validation_mpi_amr time/tlim=0.04` | Passed. The appended file has 216 records, still 40 unique tags with tag range `0..39`, and time counts `0:24`, `0.0124834403634:32`, `0.0249668827100:40`, `0.03:40`, `0.0312086052974:40`, `0.04:40`. |
| Hydro request for MHD-only field | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo.athinput -d run_validation_bad_mhd_field output1/variables=bmag` | Failed as intended with `output1/variables: field 'bmag' requires MHD` and exit status 1. |
| Unknown field name | `./build/src/athena -i inputs/particles/lagrangian_mc_thermo.athinput -d run_validation_bad_unknown_field output1/variables=not_a_field` | Failed as intended with `output1/variables: unknown field 'not_a_field'` and exit status 1. |
