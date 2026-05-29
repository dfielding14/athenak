# Example: Driven Turbulence

This public example runs continuously forced, compressible hydrodynamic
turbulence. It uses the shipped input deck `inputs/hydro/turb.athinput`, the
custom problem generator `src/pgen/turb.cpp`, and the forcing implementation
in `src/srcterms/turb_driver.cpp`.

```{warning}
The public source tree does not ship the former AMR demonstration decks,
staged MHD turbulence workflow, or tiled-driving parameters described by older
draft documentation. Do not add those settings to a public input deck.
```

## Build And Validate

The turbulence problem generator is not in the default built-in problem list.
Select it when configuring the executable:

```bash
cmake -S . -B build-turb -DPROBLEM=turb
cmake --build build-turb

# Fast initialization/forcing check
./build-turb/src/athena -i inputs/hydro/turb.athinput -d run-turb time/nlim=1
```

The short check writes history output and VTK field dumps:

```text
run-turb/Turb.hydro.hst
run-turb/vtk/Turb.hydro_w.00000.vtk
run-turb/vtk/Turb.hydro_w.00001.vtk
```

For the supplied full calculation, omit the cycle-limit override:

```bash
./build-turb/src/athena -i inputs/hydro/turb.athinput -d run-turb-full
```

## Supplied Configuration

The supplied deck defines:

| Item | Supplied value | Meaning |
| --- | --- | --- |
| Mesh | `64 x 64 x 64`, periodic | Single periodic three-dimensional domain |
| Hydro | Ideal EOS, `plm`, `hllc`, `gamma = 1.4` | Dynamical hydrodynamic solver |
| Runtime | `tlim = 10.0`, `cfl_number = 0.3` | Full supplied evolution interval |
| Output | `hst` and `vtk` (`hydro_w`) | History and primitive field output |

The public `src/pgen/turb.cpp` initializes a uniform fluid and installs its
turbulent history function. For hydro it reads optional `<problem>/d_n`
(default `1.0`); the public input relies on the default.

## Driving Parameters

Presence of `<turb_driving>` creates the public `TurbulenceDriver`. Its
parameters are read by `src/srcterms/turb_driver.cpp`:

| Parameter | Default | Purpose |
| --- | --- | --- |
| `nlow` | `1` | Lowest forced mode index |
| `nhigh` | `2` | Highest forced mode index |
| `driving_type` | `0` | Selects the driver's supported forcing geometry |
| `expo` | `5.0/3.0` | Spectral exponent used by isotropic driving |
| `exp_prp` | `5.0/3.0` | Perpendicular spectral exponent |
| `exp_prl` | `0.0` | Parallel spectral exponent |
| `dedt` | `0.0` | Energy-injection normalization |
| `tcorr` | `0.0` | Ornstein-Uhlenbeck correlation time |

The shipped input sets `tcorr = 0.5`, `dedt = 0.1`, `nlow = 1`, and
`nhigh = 2`.

## MHD Extension Boundary

`src/pgen/turb.cpp` can initialize MHD when a custom input contains an
`<mhd>` block and the executable is built with `-DPROBLEM=turb`. The generator
then accepts optional `<problem>/beta`, `<problem>/d_i`, `<problem>/d_n`, and
`<problem>/ifield` (`1` for a zero-net-flux vertical field or `2` for a
uniform vertical field). No public MHD turbulence input deck is shipped, so
this is an implementation entry point rather than a validated worked example.

## Next Steps

- [Source Terms](../modules/srcterms.md) describes how forcing is attached to
  the task list.
- [Visualization Utilities](../tools/visualization.md) describes public
  readers and plotting scripts for generated output.
- [Configuration](../configuration.md) documents the shared mesh, runtime, and
  output controls.
