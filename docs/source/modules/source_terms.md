# Source Terms

## External Gravity

The external-gravity source term lets a hydro or MHD problem evolve in a fixed,
user-prescribed gravitational potential. Add an `<external_gravity>` block to the input
file to turn it on. No `<hydro_srcterms>` or `<mhd_srcterms>` block is required.

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

When the block is enabled, hydro and MHD source-term objects are created automatically.
The source updates momentum with `rho a`, where `a = -grad(phi)`. For ideal-gas fluids it
also updates total energy by `rho v dot a`. The default hydro/MHD history output gains
`grav-PE`, the volume integral of `rho phi`.

### Built-In Models

All coordinate-dependent models are centered on `x1_origin`, `x2_origin`, and
`x3_origin`. Spherical models use `r = sqrt(x^2 + y^2 + z^2)`. Cylindrical disk/halo
models use `R = sqrt(x^2 + y^2)`. Built-in models use analytic accelerations.

| `model` | Main parameters | Form |
| --- | --- | --- |
| `uniform` | `g1`, `g2`, `g3` | Constant Cartesian acceleration. |
| `cartesian_harmonic` | `omega1`, `omega2`, `omega3` | `0.5*(omega1^2*x^2 + omega2^2*y^2 + omega3^2*z^2)`. |
| `spherical_harmonic` | `omega` | `0.5*omega^2*r^2`. |
| `point_mass`, `plummer` | `G`, `mass`, `softening` | `-G*mass/sqrt(r^2 + softening^2)`. |
| `nfw` | `G`, `rho_scale`, `r_scale` | Spherical NFW potential. |
| `logarithmic_halo` | `v0`, `core_radius`, `q` | Flattened logarithmic halo. |
| `miyamoto_nagai` | `G`, `disk_mass`, `disk_a`, `disk_b` | Miyamoto-Nagai disk. |
| `outer_cgm` | `G`, `rho_mean`, `r_outer` | Spherical outer-CGM background. |
| `gotham` | NFW, disk, and outer-CGM parameters | Sum of NFW, Miyamoto-Nagai disk, and outer-CGM terms. |
| `user` | `user_fd_step` | Calls the compiled-in user potential hook. |

Parameter parsing is model-specific: irrelevant parameters are not required and are not
added as defaults. Common aliases such as `M`, `rs`, `rho_s`, `r_soft`, and `eps` are
accepted where they match the active model.

### User Potential Hook

For a custom potential, set:

```text
<external_gravity>
model = user
user_fd_step = 1.0e-5
```

Then edit the hook in `src/pgen/pgen.cpp`:

```cpp
namespace external_gravity {

KOKKOS_FUNCTION
Real UserExternalGravityPotential(Real x1, Real x2, Real x3) {
  return /* phi(x1, x2, x3) */;
}

} // namespace external_gravity
```

The source term finite-differences this hook with `user_fd_step`. Built-in models should
be preferred when possible because their accelerations are analytic and cheaper.

### Gravity Tapering

Open boundaries can be made more robust by tapering the acceleration smoothly to zero
near selected domain faces. The physical potential used in `grav-PE` is not tapered.

```text
<external_gravity>
model         = uniform
g1            = -1.0
taper_gravity = true
taper_width_x1_inner = 0.1
taper_width_x1_outer = 0.1
```

Use `taper_width` to set all faces at once, or use face-specific widths:
`taper_width_x1_inner`, `taper_width_x1_outer`, `taper_width_x2_inner`,
`taper_width_x2_outer`, `taper_width_x3_inner`, and `taper_width_x3_outer`.

### Hydrostatic Boundaries

For hydrostatic hydro atmospheres, use the `hydrostatic_gravity` mesh boundary and an
`<external_gravity_boundary>` block:

```text
<mesh>
ix1_bc = hydrostatic_gravity
ox1_bc = hydrostatic_gravity

<external_gravity_boundary>
velocity    = no_inflow
sound_speed = 1.0
```

The hydrostatic boundary extrapolates conserved variables into ghost zones with

```text
rho_g = rho_a exp[-(phi_g - phi_a)/cs^2]
P_g   = cs^2 rho_g
```

where `a` is the adjacent active cell and `g` is the ghost cell. It uses the same
`Potential()` implementation as the source term and history output. This is the preferred
boundary for an isothermal atmosphere in a fixed potential.

Boundary options:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `velocity` | `no_inflow` | One of `copy`, `zero_normal`, `no_inflow`, or `reflect_normal`. |
| `sound_speed` | EOS/local state | Isothermal sound speed used in the exponential extrapolation. |
| `density_floor` | EOS floor | Optional ghost-zone density floor. |
| `pressure_floor` | EOS floor | Optional ghost-zone pressure floor for ideal gas. |
| `max_exponent` | `60.0` | Clamp on the exponential argument to prevent overflow. |

The hydrostatic boundary currently supports hydro. MHD setup fails clearly if a face uses
`hydrostatic_gravity`; use tapering, `diode`/`outflow`, or a problem-specific user
boundary for MHD equilibria that must account for magnetic energy.

### Example Problem Generators

Two built-in problem generators are available for external-gravity workflows:

| `pgen_name` | Use |
| --- | --- |
| `external_gravity_hydrostatic` | Initializes an isothermal hydrostatic atmosphere in the configured potential. |
| `external_gravity_disk` | Initializes a rotating hydro disk using the configured acceleration for circular speed. |

Example inputs live in:

```text
inputs/hydro/external_gravity/
inputs/mhd/external_gravity/
```

The regression input `tst/inputs/external_gravity_hydrostatic.athinput` is the minimal
tested example for `hydrostatic_gravity` boundaries.
