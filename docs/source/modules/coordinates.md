# Module: Coordinates

## Overview

The Coordinates module owns the geometry metadata needed by every MeshBlockPack.
AthenaK currently supports

* **Flat Cartesian** coordinates (default for Newtonian and SR runs).
* **Cartesian Kerr–Schild** coordinates for general-relativistic problems (either
  static background metrics or fully dynamical evolutions with ADM/Z4c).

Curvilinear grids (spherical, cylindrical, etc.) are not implemented in this branch.

## Source Layout

`src/coordinates/`

| File | Role | Notes |
| --- | --- | --- |
| `coordinates.hpp/.cpp` | Core `Coordinates` class, relativistic flags, source terms | Exposes the `CoordSrcTerms` overloads used by hydro/MHD |
| `cartesian_ks.hpp` | Metric helpers for Minkowski / Kerr–Schild space-times | Provides `ComputeMetricAndInverse` and derivatives |
| `adm.hpp/.cpp` | ADM helper class for dynamic GR runs | Populates metric slices for the Z4c/ADM modules |
| `excision.cpp` | Builds excision masks around BH horizons | Used when `coord_data.bh_excise` is true |
| `cell_locations.hpp` | Static utilities for cell/face locations | Used widely by problem generators/tasks |

## Relativistic Modes

The constructor inspects the input deck and sets:

- `is_special_relativistic` – `<coord>/special_rel = true`.
- `is_general_relativistic` – `<coord>/general_rel = true` (static metric).
- `is_dynamical_relativistic` – automatically enabled when either `<adm>` or `<z4c>`
  blocks are present alongside `<hydro>`/`<mhd>`.

If both SR and GR flags are requested simultaneously the code aborts.

When GR or dynamical GR is active, the module loads additional parameters into
`coord_data`:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `minkowski` | `false` | Treat metric as flat even in GR mode |
| `a` | — | Dimensionless BH spin (|a|<1); ignored for Minkowski |
| `excise` | `true` (if `minkowski = false`) | Enable excision masks near the horizon |
| `dexcise`, `pexcise` | — | Density/pressure floors applied inside the excision region |
| `flux_excise_r` | `1.0` without radiation; not parsed when radiation is enabled | Non-radiation input radius inside which FOFC is forced; radiation forces the horizon radius |
| `excision_scheme` | `"fixed"` (dynamic GR only) | `"fixed"` builds masks once; `"lapse"` recomputes masks from the lapse |
| `excise_lapse` | `0.25` | Lapse threshold for the `"lapse"` scheme (dynamic GR only) |

When the radiation module is active, it forces both `flux_excise_r` and the
internal `rexcise` radius to the horizon radius. Without radiation, the
public constructor defaults both to `1.0`, and only `flux_excise_r` is a
parsed input parameter.

## Coordinate Source Terms

`Coordinates::CoordSrcTerms` is overloaded:

```cpp
void CoordSrcTerms(const DvceArray5D<Real> &prim,
                   const EOS_Data &eos,
                   Real dt,
                   DvceArray5D<Real> &cons);

void CoordSrcTerms(const DvceArray5D<Real> &prim,
                   const DvceArray5D<Real> &bcc,
                   const EOS_Data &eos,
                   Real dt,
                   DvceArray5D<Real> &cons);
```

These routines rebuild the ADM 4-metric (via `ComputeMetricAndInverse`), compute
metric derivatives, and add geometric source terms to conserved variables in
the GR paths. Special-relativistic configurations do not require curvature
source terms.

## Excision Support

When `coord_data.bh_excise` is true the constructor allocates two boolean masks:

```cpp
DvceArray4D<bool> excision_floor;  // mark cells to floor primitives
DvceArray4D<bool> excision_flux;   // mark cells to force FOFC
```

For fixed excision the masks are filled once via `SetExcisionMasks`.  For lapse-based
excision the masks are recomputed each cycle through `UpdateExcisionMasks()` using
the current Z4c lapse.

Cells flagged in `excision_floor` are reset to `(dexcise, pexcise, v=0)` during
primitive recovery.  The flux mask forces first-order fluxes (`FOFC`) inside
`flux_excise_r` to stabilize the solution near the horizon.

## Cell Location Utilities

`cell_locations.hpp` provides grid-spacing helpers; the coordinate module uses these
to evaluate metrics, and problem generators frequently call them directly:

```cpp
Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
Real x1f = FaceX(i-is, nx1, x1min, x1max);
```

All grids are currently uniform Cartesian; non-uniform spacing is not supported.

## Usage Examples

### Cartesian / Newtonian (default)

```ini
<mesh>
nx1 = 256
x1min = -1.0
x1max =  1.0
# No <coord> block required.
```

### Special Relativistic Hydro/MHD

```ini
<coord>
special_rel = true

<hydro>
gamma_max = 10.0
```

### Static Kerr–Schild Background

```ini
<coord>
general_rel = true
minkowski  = false
a          = 0.7
excise     = true
dexcise    = 1.0e-8
pexcise    = 1.0e-10
flux_excise_r = 1.0   # optional; defaults to 1.0 without radiation
```

### Dynamical GRMHD With Z4c

This activation fragment is not a complete runnable deck; begin from a
shipped `inputs/dyngr/` input and add the relevant Z4c and problem-generator
controls. Dynamical relativistic fluid coordinates require MHD, not Hydro.

```ini
<z4c>
# add Z4c evolution parameters from the selected starting deck

<mhd>
eos       = ideal
dyn_eos   = ideal
dyn_error = reset_floor
rsolver   = hlle

<coord>
a = 0.5
excise = true
dexcise = 1.0e-8
pexcise = 1.0e-10
excision_scheme = lapse
excise_lapse    = 0.2
```

## Limitations & Recommendations

- Only Cartesian grids are available.  Spherical or cylindrical coordinates are not
  implemented; mimic spherical problems by embedding them in a Cartesian domain.
- Make sure `dexcise`/`pexcise` floors are compatible with the EOS used around the
  horizon.
- When enabling GR without excision (`excise = false`), ensure the mesh resolves the
  horizon to avoid C2P failures.

## Related Modules

- [Z4c Module](z4c.md) – time-evolves the dynamical metric (hands the lapse/shift to the coordinates module).
- [Dyn GRMHD Module](dyn_grmhd.md) – couples fluid evolution to the evolving space-time.
- [Mesh Module](mesh.md) – describes MeshBlockPack layout and AMR integration.
- Ghost zones must extend beyond excision region
- Use appropriate boundary conditions

### Relativistic Speeds
- Requires smaller CFL number (typically 0.2-0.4)
- Check Lorentz factor limits

## See Also
- [Mesh Module](mesh.md) - Grid structure
- [Z4c Module](z4c.md) - Numerical relativity
- [DynGRMHD Module](dyn_grmhd.md) - Relativistic MHD
- Source: `src/coordinates/`
