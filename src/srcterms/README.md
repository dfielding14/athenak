# Source Terms

## External Gravity

Add an `<external_gravity>` block to enable a fixed external gravitational potential for
hydro or MHD. No `<hydro_srcterms>` or `<mhd_srcterms>` block is required.

Example:

```text
<external_gravity>
model     = point_mass
G         = 1.0
mass      = 1.0
softening = 0.01
x1_origin = 0.0
x2_origin = 0.0
x3_origin = 0.0
```

The source term applies `a = -grad(phi)`. Built-in models use analytic acceleration;
`model=user` finite-differences the user potential hook. For ideal-gas fluids, total
energy is updated by `rho v . a`. When external gravity is enabled, default hydro/MHD
history output gains `grav-PE`, the volume integral of `rho phi`.

Built-in `model` values:

- `uniform`: Cartesian linear potential from constant acceleration components `g1`,
  `g2`, and `g3`.
- `cartesian_harmonic`: `0.5*(omega1^2*x^2 + omega2^2*y^2 + omega3^2*z^2)`.
- `spherical_harmonic`: `0.5*omega^2*r^2`.
- `point_mass` or `plummer`: `-G*mass/sqrt(r^2 + softening^2)`.
- `nfw`: spherical NFW potential with `rho_scale` and `r_scale`.
- `logarithmic_halo`: cylindrical logarithmic halo with `v0`, `core_radius`, and `q`.
- `miyamoto_nagai`: cylindrical disk with `disk_mass`, `disk_a`, and `disk_b`.
- `outer_cgm`: spherical outer CGM background with `rho_mean` and `r_outer`.
- `gotham`: NFW plus Miyamoto-Nagai disk plus outer CGM background, matching the
  useful form previously used in gotham problem generators.
- `user`: calls `external_gravity::UserExternalGravityPotential(x1, x2, x3)` in
  `src/pgen/pgen.cpp`; edit that function for a problem-specific compiled-in potential.

All coordinate-dependent models are centered on `x1_origin`, `x2_origin`, `x3_origin`.
The spherical models use `r = sqrt(x^2 + y^2 + z^2)`, and the cylindrical models use
`R = sqrt(x^2 + y^2)`.

Optional boundary tapering can smoothly turn the acceleration off near domain faces:

```text
<external_gravity>
model         = uniform
g1            = -1.0
taper_gravity = true
taper_width_x1_inner = 0.1
taper_width_x1_outer = 0.1
```

For hydrostatic hydro atmospheres, set a mesh face to `hydrostatic_gravity` and add
`<external_gravity_boundary>`. The boundary extrapolates conserved variables into ghost
zones using the same potential as the source term:

```text
<mesh>
ix1_bc = hydrostatic_gravity
ox1_bc = hydrostatic_gravity

<external_gravity_boundary>
velocity    = no_inflow
sound_speed = 1.0
```

The hydrostatic boundary currently supports hydro. MHD runs that request
`hydrostatic_gravity` fail at setup with a clear error; for MHD use gravity tapering or a
problem-specific user boundary condition.
