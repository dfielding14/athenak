<!--
GitHub Pages integration checklist:
1. Copy this file to docs/source/modules/lagrangian_mc_tracers.md on origin/gh-pages.
2. Add modules/lagrangian_mc_tracers to the Modules toctree in docs/source/index.md.
3. Add a row for this page under Physics Modules in docs/source/modules/index.md.
4. Add a short link from docs/source/modules/particles.md and, optionally,
   docs/source/modules/outputs.md.
5. Rebuild with: cd docs && make html.
-->

# Module: Lagrangian Monte-Carlo Thermodynamic Tracers

## Overview

`particle_type = lagrangian_mc` adds Monte-Carlo tracer particles that follow the
mass flux of the fluid. Each tracer sits at a cell center and jumps to neighboring
cell centers with probabilities derived from the RK-accumulated density flux through
the current cell faces.

The tracer particle arrays store only transport state: position, tag, owning
MeshBlock, seed id, creation time, and AMR bookkeeping. Requested thermodynamic and
derived variables are sampled from the cell-centered Hydro or MHD state when
`prtcl_thermo_history` output is written. This keeps restarts independent of the
chosen diagnostic variable list.

This feature is intended for ideal-gas Hydro and MHD. Non-ideal and isothermal
thermodynamic-history requests should fail during input validation.

## Source Location

| Path | Role |
| --- | --- |
| `src/particles/particles_lagrangian_mc.cpp` | Seeding schedules, stochastic pusher, AMR remapping, and restart persistence. |
| `src/particles/tracer_fields.hpp` | Tracer-field registry interface and field metadata. |
| `src/particles/tracer_fields.cpp` | Field-name parser and host-side field evaluator. |
| `src/outputs/prtcl_thermo_history.cpp` | Versioned append-only binary history output. |
| `scripts/read_prtcl_thermo_history.py` | Python reader for v1 and v2 history files. |
| `scripts/test_prtcl_thermo_history_reader.py` | Synthetic v2 reader regression test. |
| `inputs/particles/lagrangian_mc_thermo*.athinput` | Hydro, AMR, and MHD smoke-test inputs. |

## Quick Start

Minimal Hydro setup:

```ini
<hydro>
eos         = ideal
gamma       = 1.6666666666666667
nscalars    = 1

<particles>
particle_type    = lagrangian_mc
pusher           = lagrangian_mc
random_seed      = 12345
track_variables  = default, mach, sound_speed, vmag, scalar0

<tracer_seed1>
id              = 1
start_time      = 0.0
end_time        = 0.0
cadence         = -1.0
count_per_event = 1000
weight          = mass
region          = all
seed            = 24680

<output1>
file_type = prtcl_thermo_history
dt        = 0.01
variables = density, temperature, mach, scalar0
```

Run the serial smoke test:

```bash
./build/src/athena \
  -i inputs/particles/lagrangian_mc_thermo.athinput \
  -d run_lmc_hydro
```

Read the history file:

```bash
python scripts/read_prtcl_thermo_history.py \
  run_lmc_hydro/prtcl_thermo_history/lagrangian_mc_thermo.prtcl_thermo_history.thp
```

## Runtime Blocks

### `<particles>`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `particle_type` | string | required | Set to `lagrangian_mc`. |
| `pusher` | string | required | Set to `lagrangian_mc`. |
| `random_seed` | integer | `12345` | Base seed used by tracer operations. |
| `track_variables` | string list | `default` | Default tracer-field list for `prtcl_thermo_history` outputs that do not specify their own `variables`. |
| `next_tracer_tag` | integer | `0` | Restart-managed next global tracer tag. Users normally do not set this manually. |

### `<tracer_seedN>`

Add one or more numbered seed blocks, for example `<tracer_seed1>`,
`<tracer_seed2>`, and so on.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `id` | integer | block order | Stored as `seed_id` in history output. |
| `start_time` | Real | `0.0` | First event time. Events due at the initial time fire after problem-generator setup. |
| `end_time` | Real | `start_time` | Last allowed event time. |
| `cadence` | Real | `-1.0` | Time between seed events. `cadence <= 0` means one-shot. |
| `count_per_event` | integer | required | Global number of particles requested per event. |
| `weight` | string | `mass` | `mass` samples cells by `rho*dV`; `volume` samples by `dV`. |
| `region` | string | `all` | One of `all`, `box`, `sphere`, or `slab`. |
| `seed` | integer | `0` | Deterministic seed for this schedule. |
| `target` | string | none | Optional tracer-field selector used as a seed mask. |
| `target_min` | Real | none | Minimum accepted target value. |
| `target_max` | Real | none | Maximum accepted target value. |

Region parameters:

| Region | Parameters |
| --- | --- |
| `all` | No extra parameters. |
| `box` | `x1min`, `x1max`, `x2min`, `x2max`, `x3min`, `x3max`. |
| `sphere` | `center1`, `center2`, `center3`, `radius`. |
| `slab` | `slab_axis`, `slab_min`, `slab_max`. |

New particles are inserted at selected cell centers. This matches the cell-jump
Monte-Carlo transport model and avoids subcell interpolation ambiguity.

### `<outputN>` with `file_type = prtcl_thermo_history`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `file_type` | string | required | Set to `prtcl_thermo_history`. |
| `dt` or `dcycle` | Real or integer | required | Output cadence. |
| `variables` | string list | `<particles>/track_variables` | Per-output tracer-field list. |
| `id` | string | `prtcl_thermo_history` | Output id used in the `.thp` filename. |

## Tracer Field Registry

The same field names are used by:

- `<particles>/track_variables`
- `<outputN>/variables` for `prtcl_thermo_history`
- `<tracer_seedN>/target`

Names are case-insensitive and may be separated by commas or whitespace. Aliases
are normalized to the canonical names stored in the v2 `.thp` file header.

### Hydro and MHD Fields

| Canonical name | Accepted aliases | Definition |
| --- | --- | --- |
| `density` | `rho` | Cell density. |
| `pressure` | `p` | Ideal-gas pressure, `(gamma - 1)*eint`. |
| `temperature` | `T` | `pressure/density`. |
| `entropy` | `s`, `specific_entropy` | `log(pressure/density^gamma)`. |
| `internal_energy` | `eint` | Internal-energy primitive slot. |
| `v1`, `v2`, `v3` | `vx`, `vy`, `vz` | Fluid velocity components. |
| `vmag` | `velocity_magnitude` | `sqrt(v1^2 + v2^2 + v3^2)`. |
| `sound_speed` | `cs` | Ideal-gas sound speed. |
| `mach` | `mach_number` | `vmag/sound_speed`. |
| `scalarN` | none | Passive scalar with zero-based index `N`, for example `scalar0`. |

`default` expands to:

```text
density, pressure, temperature, entropy, internal_energy, v1, v2, v3
```

### MHD-Only Fields

| Canonical name | Accepted aliases | Definition |
| --- | --- | --- |
| `b1`, `b2`, `b3` | `bx`, `by`, `bz` | Cell-centered magnetic-field components. |
| `bmag` | `magnetic_field_magnitude` | `sqrt(b1^2 + b2^2 + b3^2)`. |
| `magnetic_pressure` | `pmag` | `0.5*bmag^2`. |
| `beta` | `plasma_beta` | Gas pressure divided by magnetic pressure. |
| `alfven_speed` | `va` | `bmag/sqrt(density)`. |

Input validation is deliberately strict. Unknown fields, MHD-only fields in Hydro
runs, and out-of-range `scalarN` requests are fatal errors.

## History File Format

`prtcl_thermo_history` writes append-only binary `.thp` files under
`prtcl_thermo_history/`.

Each record contains fixed metadata followed by the requested variables in schema
order:

```text
time, cycle, tag, seed_id, x1, x2, x3, gid, <requested variables...>
```

The v2 file header stores:

- magic string
- schema version
- `Real` size
- number of requested fields
- newline-separated canonical field names

Each appended block stores:

- block magic string
- schema version
- number of records
- integer and real record widths
- cycle and time
- integer record payload
- real record payload

Append behavior is strict. If a restart or rerun attempts to append to an existing
`.thp` with a different variable list, schema version, or real precision, AthenaK
exits with a fatal schema-mismatch error instead of appending incompatible records.

## Restarts, MPI, and AMR

Restart files persist:

- particle real and integer arrays
- global next tracer tag
- each seed schedule's next fire time
- each seed schedule's event index and completion flag

Restarted runs continue tag assignment monotonically and do not duplicate seed
events that fired before the restart. For MPI restarts, write restart files with
`single_file_per_rank = true` and restart from rank 0's file; the restart reader
maps the companion rank files internally.

After AMR or load balancing changes the MeshBlock layout, tracers are remapped to
the owning rank and snapped back to valid cell centers on the new mesh.

## Example Inputs

| Input | Purpose |
| --- | --- |
| `inputs/particles/lagrangian_mc_thermo.athinput` | Serial Hydro seeding, scalar tracking, derived-Mach seed mask, and custom history variables. |
| `inputs/particles/lagrangian_mc_thermo_amr.athinput` | Hydro AMR plus MPI-safe restart output using `<particles>/track_variables = density, temperature, mach`. |
| `inputs/particles/lagrangian_mc_thermo_mhd.athinput` | MHD smoke test with `bmag`, `beta`, `alfven_speed`, and `mach`. |

## Reader Usage

Use `scripts/read_prtcl_thermo_history.py` to inspect or convert the binary file:

```bash
python scripts/read_prtcl_thermo_history.py path/to/file.thp --npz tracers.npz
```

The reader returns arrays keyed by column name. For v2 files, the data keys come
from the stored schema, so a file written with:

```ini
variables = density, temperature, mach, scalar0
```

returns `density`, `temperature`, `mach`, and `scalar0` arrays in addition to the
fixed metadata arrays. The reader also supports previous fixed-layout v1 files.

## Validation

The feature branch validation covered:

| Check | Result |
| --- | --- |
| Serial build | `cmake --build build -j 8` passed. |
| MPI build | `cmake --build build_mpi -j 8` passed. |
| Reader regression | `python3 scripts/test_prtcl_thermo_history_reader.py` passed. |
| Whitespace | `git diff --check` passed. |
| Serial Hydro custom variables | 72 records, 32 unique tags, finite `density`, `temperature`, `mach`, and `scalar0`. |
| Same-schema restart append | Passed; appended history retained 32 unique tags. |
| Changed-schema append | Failed as intended with a schema-mismatch fatal error. |
| Serial MHD derived fields | 48 records with finite `bmag`, `beta`, `alfven_speed`, and `mach`. |
| Two-rank MPI AMR | 136 records, 40 unique tags, finite custom fields, and valid spherical seed-mask locations. |
| Two-rank MPI AMR restart append | 216 appended records, still 40 unique tags, no duplicated timed seed events. |

## Current Limitations

- The built-in fields are intentionally finite and explicit; there is no runtime
  expression parser.
- Problem-specific derived fields should be stored in passive scalars for now and
  tracked as `scalarN`.
- A future extension can add pgen-enrolled tracer-field callbacks for quantities
  such as `cooling_time` or `tcool_over_tff`.
- The history output samples complete trajectories at the requested output cadence,
  not every RK substage.

## See Also

- [Particles](particles.md)
- [Outputs](outputs.md)
- [Hydro](hydro.md)
- [MHD](mhd.md)
- Source-adjacent notes: `inputs/particles/lagrangian_mc_thermo.md`
