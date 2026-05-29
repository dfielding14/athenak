# Module: Reconstruction

Hydro selects reconstruction through `<hydro>/reconstruct`; MHD and DynGRMHD
select it through `<mhd>/reconstruct`; radiation selects it through its
radiation block. Public reconstruction kernels are inline implementations in
`src/reconstruct/` for Cartesian-like uniform spacing.

## Available Methods

| Input keyword | Source | Implemented method | Ghost zones without FOFC |
| --- | --- | --- | --- |
| `dc` | `dc.hpp` | Piecewise constant donor cell | `>= 2` |
| `plm` | `plm.hpp` | Limited piecewise linear | `>= 2` |
| `ppm4` | `ppm.hpp` | Colella-Woodward PPM limiter path | `>= 3` |
| `ppmx` | `ppm.hpp` | Colella-Sekora extrema-preserving limiter path | `>= 3` |
| `wenoz` | `wenoz.hpp` | Fifth-order WENO-Z stencil | `>= 3` |

Hydro/MHD first-order flux correction changes the requirements to
`nghost >= 3` for `plm` and `nghost >= 4` for `ppm4`, `ppmx`, or `wenoz`.
Refined meshes require even `nghost`, so their effective minimum is four when
three would otherwise suffice.

## Configuration

```ini
<hydro>
reconstruct = plm
```

```ini
<mhd>
reconstruct = ppm4
```

```ini
<radiation>
reconstruct = plm
```

The constructor of the active module validates the keyword and ghost-zone
requirements. Choice of method is problem-dependent; this documentation does
not promise a universal accuracy or robustness ranking.

## See Also

- [Hydrodynamics](hydro.md)
- [Magnetohydrodynamics](mhd.md)
- [Configuration](../configuration.md)
