# Module: Equations Of State

Equation-of-state selection is performed by fluid constructors. The public
standard Hydro and MHD paths expose ideal and isothermal choices; the
dynamical-GR MHD coupling uses a separate primitive-solver policy selected
through `<mhd>/dyn_eos`.

## Standard Hydro/MHD EOS

| Fluid block | `eos` value | Implementation |
| --- | --- | --- |
| `<hydro>` | `ideal` | `ideal_hyd.cpp`, or SR/GR variants selected by coordinates |
| `<hydro>` | `isothermal` | `isothermal_hyd.cpp`; non-relativistic only |
| `<mhd>` | `ideal` | `ideal_mhd.cpp`, or SR/GR variants selected by coordinates |
| `<mhd>` | `isothermal` | `isothermal_mhd.cpp`; non-relativistic only |

```ini
<hydro>
eos = ideal
gamma = 1.4
```

```ini
<mhd>
eos = isothermal
iso_sound_speed = 1.0
```

The Hydro/MHD constructors reject `eos = isothermal` in special- or
general-relativistic coordinate modes. Relevant controls also include
`gamma_max` in SR/GR ideal-fluid constructors and fluid floors parsed by the
selected EOS implementation.

## Dynamical GRMHD Primitive-Solver EOS

When `<mhd>` is paired with `<adm>` or `<z4c>`,
`src/dyn_grmhd/dyn_grmhd.cpp` builds a primitive-solver EOS policy. Z4c uses
the staged matter-feedback route; an ADM-only path is not equivalent to a
Z4c matter-coupled spacetime evolution.

| Required `<mhd>/dyn_eos` | Public policy path |
| --- | --- |
| `ideal` | Ideal-gas primitive solver |
| `piecewise_poly` | Piecewise-polytrope primitive solver |
| `compose` | CompOSE table primitive solver |
| `hybrid` | Hybrid EOS primitive solver |

The required `<mhd>/dyn_error` currently accepts only `reset_floor`. Further
policy-specific controls are parsed in `src/eos/primitive_solver_hyd.hpp` and
the corresponding policy implementation; use a shipped `inputs/dyngr/` deck
as the starting point rather than inferring a generic table interface.

Dynamic GRMHD currently also requires the ordinary `<mhd>/eos = ideal`
setting so that the MHD state includes the energy storage read by the
dynamical-GRMHD update path.

## Source Map

| Source area | Role |
| --- | --- |
| `src/eos/eos.hpp`, `eos.cpp` | Shared standard-EOS data/interface |
| `src/eos/ideal_*`, `src/eos/isothermal_*` | Standard Hydro/MHD conversions |
| `src/eos/noop_dyngrmhd.cpp` | Base MHD placeholder when DynGRMHD owns primitive conversion |
| `src/eos/primitive_solver_hyd.hpp` | DynGRMHD primitive-solver wiring and common controls |
| `src/eos/primitive-solver/` | Policy implementations for ideal, piecewise, CompOSE, and hybrid routes |

## See Also

- [Hydrodynamics](hydro.md)
- [Magnetohydrodynamics](mhd.md)
- [Dynamical GRMHD](dyn_grmhd.md)
