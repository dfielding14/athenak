# Module: Diffusion

The public diffusion implementation contributes explicit viscous,
conductive, or resistive terms to Hydro/MHD updates. These helpers are
constructed from parameters in the active `<hydro>` or `<mhd>` block; there
is no separate public diffusion input block.

## Implemented Paths

| Source | Activated by | Role |
| --- | --- | --- |
| `src/diffusion/viscosity.*` | `viscosity` in `<hydro>` or `<mhd>` | Isotropic kinematic shear viscosity |
| `src/diffusion/conduction.*` | `conductivity` or `tdep_conductivity` in `<hydro>` or `<mhd>` | Heat flux and conductive timestep bound |
| `src/diffusion/resistivity.*` | `ohmic_resistivity` in `<mhd>` | Ohmic electric field and resistive timestep bound |
| `src/diffusion/current_density.hpp` | Resistive path | Current-density helpers |

## Inputs

```ini
<hydro>
eos = ideal
rsolver = hlle
viscosity = 1.0e-3
conductivity = 1.0e-2

<mhd>
eos = ideal
rsolver = hlle
ohmic_resistivity = 1.0e-3
```

| Parameter | Module(s) | Verified parser behavior |
| --- | --- | --- |
| `viscosity` | Hydro, MHD | Required value when viscosity helper is constructed |
| `conductivity` | Hydro, MHD | Constant thermal conductivity; default `0.0` within constructed helper |
| `tdep_conductivity` | Hydro, MHD | Boolean enabling temperature-dependent conductivity |
| `cond_ceiling` | Hydro, MHD | Ceiling for temperature-dependent conductivity |
| `sat_hflux` | Hydro, MHD | Boolean saturated-heat-flux option |
| `ohmic_resistivity` | MHD only | Required Ohmic coefficient when resistivity helper is constructed |

Thermal conduction currently requires an ideal-gas fluid; the constructor
rejects conduction with an isothermal EOS.

`sat_hflux` is meaningful only with `tdep_conductivity = true`. With constant
conductivity, enabling it disables the usual conduction timestep bound without
selecting the saturated-flux calculation path.

DynGRMHD bypasses the standard viscosity and conduction flux-task paths.
Resistivity has only partial routing through its corner electric-field
calculation. Do not treat the standard diffusion settings on this page as
supported DynGRMHD behavior.

## Timestep Consequences

Viscosity and resistivity compute explicit diffusion timestep limits using
the smallest cell spacing and dimensional prefactors in their constructors.
Conduction supplies its own timestep limit from the configured conductivity
path. Large coefficients can therefore reduce the allowed simulation
timestep substantially; the public module does not advertise an implicit
diffusion alternative.

## Shipped Evidence

- `inputs/hydro/viscosity.athinput` supplies a Hydro viscosity example.
- `inputs/mhd/resistivity.athinput` supplies an MHD Ohmic-resistivity example.

Inspect those decks and run a short cycle-limited case before adapting a
diffusion configuration to a new problem.

## See Also

- [Hydrodynamics](hydro.md)
- [Magnetohydrodynamics](mhd.md)
- [Configuration](../configuration.md)
