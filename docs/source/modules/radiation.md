# Module: Radiation

The public Radiation module evolves angle-resolved radiation intensities on a
geodesic angular mesh. It does not implement the M1-closure moment interface
described by earlier draft documentation.

## Public Representation

`src/radiation/radiation.hpp` owns:

| Member | Role |
| --- | --- |
| `GeodesicGrid *prgeo` | Angular mesh selected by `nlevel` |
| `i0`, `i1` | Intensity arrays indexed by radiation angle and cell |
| `iflx` | Spatial face fluxes of intensity |
| `divfa` | Angular-flux divergence when `angular_fluxes = true` |
| tetrad arrays | Coordinate/tetrad transformations used by the transport kernels |

Radiation construction currently requires `<coord>/general_rel = true`; flat
spacetime tests use `<coord>/minkowski = true`. Radiation may run alone or
alongside either Hydro or MHD, but not both simultaneously.

## Implementation Files

| File | Responsibility |
| --- | --- |
| `src/radiation/radiation.cpp`, `radiation.hpp` | State, parameter parsing, angular-grid allocation |
| `src/radiation/radiation_fluxes.cpp` | Spatial and angular transport fluxes |
| `src/radiation/radiation_update.cpp` | Intensity update |
| `src/radiation/radiation_source.cpp` | Radiation-fluid coupling |
| `src/radiation/radiation_tetrad.cpp` | Tetrad geometry |
| `src/radiation/radiation_tasks.cpp` | Coupled and uncoupled task sequencing |

## `<radiation>` Parameters

| Parameter | Requirement/default | Behavior |
| --- | --- | --- |
| `nlevel` | Required | Geodesic angular mesh level |
| `rotate_geo` | Default `true` | Rotation option for the angular mesh |
| `angular_fluxes` | Default `true` | Enable angular flux-divergence update |
| `n_0_floor` | Default `0.1` | Angular transport floor used by the implementation |
| `reconstruct` | Default `plm` | `dc`, `plm`, `ppm4`, `ppmx`, or `wenoz`; high-order choices require `nghost >= 3` |
| `fixed_fluid` | Default `false` | Holds an accompanying fluid fixed when selected |

When `<hydro>` or `<mhd>` is present, radiation-fluid coupling is enabled by
default through `rad_source = true`. In that case:

| Parameter | Requirement/default |
| --- | --- |
| `rad_source` | Default `true` with a fluid module |
| `kappa_s` | Required when `rad_source = true` |
| `power_opacity` | Default `false` |
| `kappa_a`, `kappa_p` | Required when `rad_source = true` and `power_opacity = false` |
| `compton` | Default `false`; enabling it requires a `<units>` block |
| `arad` | Required without `<units>`; computed from units when present |
| `affect_fluid` | Default `true` |

## Beam Source Term

The current beam interface belongs in `<rad_srcterms>`, as documented in
[Source Terms](srcterms.md):

```ini
<rad_srcterms>
rad_beam = true
dii_dt   = 1.0
pos_1    = 3.91
pos_2    = 0.0
pos_3    = 0.0
dir_1    = 0.0
dir_2    = 1.0
dir_3    = 0.0
width    = 0.7
spread   = 10.0
```

`inputs/radiation/bh_beam.athinput` contains that current source-term block,
but its `<problem>` block omits the built-in selector required by
`src/pgen/tests/rad_beam.cpp`. To validate the deck with the default build,
work from a copy and add:

```ini
<problem>
pgen_name = rad_beam
```

For example, after adding that selector to a copy:

```bash
./build/src/athena -i my_bh_beam.athinput -d run-beam time/nlim=1
```

With `basename = beam_audit` in the audit copy, a verified one-cycle run
writes binary radiation-coordinate output:

```text
run-beam/bin/beam_audit.rad_coord.00000.bin
run-beam/bin/beam_audit.rad_coord.00001.bin
```

The older `inputs/radiation/beam.athinput` and `inputs/radiation/snake.athinput`
place beam parameters in blocks that are not consumed by the current public
beam-source constructor; they should not be used as runnable beam guidance
without correction.

## Task Coupling

`Radiation::AssembleRadTasks()` inserts transport, source, boundary,
prolongation, and update tasks. When an evolving Hydro or MHD module is
present, the radiation tasks are sequenced with that module and apply
`RadFluidCoupling`; when no fluid evolves, the radiation transport tasks run
without a fluid update.

## Outputs

Radiation-derived output selections include values beginning with `rad_coord`,
`rad_fluid`, `rad_hydro`, and `rad_mhd` where the corresponding module state
exists; output resolution is implemented in
`src/outputs/basetype_output.cpp` and `src/outputs/derived_variables.cpp`.

## See Also

- [Source Terms](srcterms.md)
- [Outputs](outputs.md)
- [Public Input Parameters](../reference/input_parameters.md)
