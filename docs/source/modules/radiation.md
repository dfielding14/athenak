# Module: Radiation

## Overview
The Radiation module solves frequency-integrated radiative transport using a
discrete ordinates ($S_n$) scheme on a geodesic angular grid.  Specific
intensities are evolved along a set of directions and coupled to the fluid
through implicit source terms.  Both special and general relativistic
spacetimes are supported.

**Limitations**

- Requires a general relativistic coordinate system.
- Adaptive mesh refinement is currently unsupported.
- Only one of hydrodynamics or MHD may be enabled at a time.

## Source Location
`src/radiation/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `radiation.hpp/cpp` | Core radiation class | Initialization, intensity evolution |
| `radiation_fluxes.cpp` | Flux computation along angular directions | `CalculateFluxes()` |
| `radiation_source.cpp` | Matter coupling | Emission, absorption, scattering |
| `radiation_update.cpp` | Implicit update of intensities | `ImplicitUpdate()` |
| `radiation_tetrad.hpp/cpp` | Tetrad formalism | Frame transformations |
| `radiation_tasks.cpp` | Task management | Task registration |
| `radiation_opacities.hpp` | Opacity models | κ calculations |

## Equations Solved
For each geodesic direction $n^\mu$ the module evolves the specific intensity
$I$:

\[
\frac{1}{c}\frac{\partial I}{\partial t} + n^i \nabla_i I =
-(\kappa_a+\kappa_s)I + \kappa_a B + \kappa_s J ,
\]

where $B$ is the equilibrium intensity and $J$ is the angle-averaged
intensity.  Moments such as radiation energy density and flux are computed from
the intensities when needed for coupling to the fluid or for diagnostic output.

## Configuration Parameters

From `<radiation>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `arad` | Real | required | Radiation constant |
| `kappa_a` | Real | required | Absorption opacity |
| `kappa_p` | Real | required | Planck mean opacity |
| `kappa_s` | Real | required | Scattering opacity |
| `power_opacity` | bool | false | Use power-law opacity |
| `compton` | bool | false | Enable Compton scattering (requires `<units>`) |
| `beam_source` | bool | false | Collimated beam source |
| `rad_source` | bool | true | Couple radiation to fluid |
| `affect_fluid` | bool | true | Radiation acts back on fluid |
| `fixed_fluid` | bool | false | Disable fluid evolution |
| `nlevel` | int | required | Geodesic grid refinement level |
| `rotate_geo` | bool | true | Rotate geodesic grid each step |
| `angular_fluxes` | bool | true | Evolve angular flux terms |
| `reconstruct` | string | plm | Spatial reconstruction method |
| `n_0_floor` | Real | 0.1 | Flux limiter floor |

## Execution Flow

```{mermaid}
flowchart TD
    Start[Radiation Tasks] --> Source[Calculate Sources]
    Source --> Opacity[Update Opacities]
    Opacity --> Flux[Calculate Fluxes]
    Flux --> Implicit[Implicit Update]
    Implicit --> Couple[Fluid Coupling]
    Couple --> BC[Apply BCs]
    BC --> Done[Complete]
```

## Tetrad Formalism (GR)

The module constructs orthonormal tetrads to transform between coordinate and
fluid frames.  These tetrads are used to compute moments and perform Lorentz
transforms when coupling to the fluid.

## Implicit Solver

Radiation–matter coupling is solved implicitly using a Newton–Raphson
iteration.  The scheme forms a Jacobian of residuals with respect to the
intensities and updates them until convergence.

## Variable Arrays

### Radiation Variables
- `i0`: specific intensity per angle
- `i1`: reconstruction scratch array
- `iflx`: face fluxes per angle
- `divfa`: angular divergence terms

### Geometry and Coupling
- `tet_c`, `tetcov_c`: tetrad and metric data
- `na`: angular face areas
- `norm_to_tet`: transformation from fluid velocity to tetrad frame

Derived moments are available via the `rad_coord` (coordinate frame) and
`rad_fluid` (fluid frame) derived variables.

## Common Test Problems

### Radiation Diffusion
```ini
<radiation>
arad = 1.0
kappa_a = 100.0  # Optically thick
kappa_p = 100.0
kappa_s = 0.0
nlevel = 1
```

### Radiation Beam
```ini
<radiation>
beam_source = true
kappa_a = 0.1
kappa_p = 0.1
kappa_s = 0.0
nlevel = 1
```

## Performance Considerations

- Implicit solves dominate cost in optically thick regimes.
- Memory usage scales with the number of angular directions.
- Reduce the Courant number when coupling is stiff.

## Common Issues

- Use sufficiently small time steps for convergence.
- Ensure required opacities and `nlevel` are specified.

## See Also
- [Hydro Module](hydro.md)
- Source: `src/radiation/radiation.cpp`

