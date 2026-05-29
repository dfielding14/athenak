# Module: Shearing Box

The public shearing-box code provides shearing-periodic x1 boundaries,
orbital-advection helpers, and shearing-box source terms for Hydro/MHD
problems. The operational shifted-boundary/orbital-advection path documented
here is currently three-dimensional.

## Configuration

```ini
<mesh>
ix1_bc = shear_periodic
ox1_bc = shear_periodic

<shearing_box>
qshear = 1.5
omega0 = 1.0
stratified = false
```

| Parameter | Requirement/default | Parsed by |
| --- | --- | --- |
| `qshear` | Required | `ShearingBox` and `OrbitalAdvection` constructors |
| `omega0` | Required | `ShearingBox` and `OrbitalAdvection` constructors |
| `stratified` | Default `false` | `ShearingBox` constructor |

## Constraints

The mesh constructor enforces that:

- `shear_periodic` appears on both x1 faces together.
- A shearing-periodic mesh is two- or three-dimensional.
- A `<shearing_box>` block is present when those boundary flags are used.
- Shearing-box configurations are not currently compatible with SMR/AMR.
- `shear_periodic` is rejected on x2 and x3 faces.

The parser permits a two-dimensional shearing-periodic mesh, and shearing-box
source terms include 2D branches. However, the shifted-boundary and
orbital-advection communication/remap paths execute only for
three-dimensional runs unless the unimplemented radial-azimuthal path is
enabled. Treat 2D shearing-periodic run workflows as unsupported.

## Implementation Map

| Source | Responsibility |
| --- | --- |
| `src/shearing_box/shearing_box.*` | Parameters and boundary buffer coordination |
| `shearing_box_cc.cpp`, `shearing_box_fc.cpp` | Cell- and face-centered shifted boundaries |
| `shearing_box_srcterms.cpp` | Fluid/MHD shearing-box source updates |
| `orbital_advection*.cpp` | Orbital-advection operations |
| `remap_fluxes.hpp` | Conservative remap helper |
| `shearing_box_tasks.cpp`, `orbital_advection_tasks.cpp` | Task-list wrappers |

Hydro constructs its cell-centered shearing/orbital objects whenever a
`<shearing_box>` block exists. MHD constructs both cell- and face-centered
objects so magnetic fields participate in the shearing-boundary workflow.

## Verified Example

`inputs/shearing_box/mri3d_unstratified.athinput` uses the built-in
`pgen_name = mri3d` route and was executed for one cycle during this
documentation audit:

```bash
./build/src/athena -i inputs/shearing_box/mri3d_unstratified.athinput \
  -d run-mri3d time/nlim=1
```

See [MRI In A Shearing Box](../examples/mri_turbulence.md) for outputs and
the current non-runnable status of the custom `mri2d` deck.

## See Also

- [Magnetohydrodynamics](mhd.md)
- [Boundaries](boundaries.md)
- [MRI In A Shearing Box](../examples/mri_turbulence.md)
