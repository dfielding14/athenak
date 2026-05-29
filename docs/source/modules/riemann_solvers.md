# Module: Riemann Solvers

Hydro and MHD constructors parse `rsolver` from their fluid input block and
map the same public keyword to a non-relativistic, SR, or GR implementation
according to the configured coordinates. Users do not select internal enum
suffixes such as `hlle_sr` in an input deck.

## Hydro Keywords

| Regime | Accepted `rsolver` keywords | Source directory |
| --- | --- | --- |
| Non-relativistic dynamic | `llf`, `hlle`, `hllc`, `roe` | `src/hydro/rsolvers/` |
| Non-relativistic kinematic | `advect` | `src/hydro/rsolvers/advect_hyd.hpp` |
| Special relativistic dynamic | `llf`, `hlle`, `hllc` | SR headers under `src/hydro/rsolvers/` |
| General relativistic dynamic | `llf`, `hlle` | GR headers under `src/hydro/rsolvers/` |

For non-relativistic Hydro, `hllc` is rejected with isothermal EOS. SR/GR
kinematic Hydro is rejected.

## MHD Keywords

| Regime | Accepted `rsolver` keywords | Source directory |
| --- | --- | --- |
| Non-relativistic dynamic | `llf`, `hlle`, `hlld` | `src/mhd/rsolvers/` |
| Non-relativistic kinematic | `advect` | `src/mhd/rsolvers/advect_mhd.hpp` |
| Special relativistic dynamic | `llf`, `hlle` | SR headers under `src/mhd/rsolvers/` |
| General relativistic dynamic | `llf`, `hlle` | GR headers under `src/mhd/rsolvers/` |
| DynGRMHD with `<adm>` or `<z4c>` | `llf`, `hlle` | `src/dyn_grmhd/rsolvers/` |

The public non-relativistic HLLD source includes ideal-gas and isothermal
paths. SR/GR kinematic MHD is rejected.

## Input Examples

```ini
<hydro>
eos = ideal
rsolver = hllc
```

```ini
<mhd>
eos = ideal
rsolver = hlld
```

For a relativistic configuration the input still uses, for example,
`rsolver = hlle`; static/SR/GR MHD coordinates choose the ordinary MHD
variant, while DynGRMHD independently dispatches its `llf_dyngr` or
`hlle_dyngr` implementation for `<mhd>` paired with `<adm>` or `<z4c>`.

## Validation

The shipped Sod deck uses a built-in Hydro problem generator; the shipped
Orszag-Tang deck exercises MHD. See [Quickstart](../quickstart.md) for
executed commands and observed output paths.

## See Also

- [Hydrodynamics](hydro.md)
- [Magnetohydrodynamics](mhd.md)
- [Reconstruction](reconstruction.md)
