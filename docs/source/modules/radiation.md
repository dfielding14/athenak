# Module: Radiation

## Overview
The Radiation module implements radiation transport using the M1 closure scheme, supporting both optically thin and thick regimes with full special and general relativistic treatments.

## Source Location
`src/radiation/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `radiation.hpp/cpp` | Core radiation class | Initialization, closures |
| `radiation_fluxes.cpp` | Radiation flux calculation | `CalculateFluxes()` |
| `radiation_source.cpp` | Source terms | Emission, absorption |
| `radiation_update.cpp` | Implicit update | `ImplicitUpdate()` |
| `radiation_tetrad.hpp/cpp` | Tetrad formalism | Frame transformations |
| `radiation_tasks.cpp` | Task management | Task registration |
| `radiation_opacities.hpp` | Opacity models | κ calculations |

## Equations Solved

### M1 Radiation Transport
$$\frac{\partial E_r}{\partial t} + \nabla \cdot \mathbf{F}_r = -\kappa_a(E_r - aT^4) - \kappa_s \cdot \mathbf{F}_r/c$$
$$\frac{\partial \mathbf{F}_r}{\partial t} + \nabla \cdot \mathbf{P}_r = -\kappa_a \mathbf{F}_r - (\kappa_a + \kappa_s)\mathbf{F}_r/c$$

Where:
- $E_r$: Radiation energy density
- $\mathbf{F}_r$: Radiation flux
- $\mathbf{P}_r$: Radiation pressure tensor
- $\kappa_a$: Absorption opacity
- $\kappa_s$: Scattering opacity

## Configuration Parameters

From `<radiation>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `arad` | Real | required | Radiation constant |
| `kappa_a` | Real | required | Absorption opacity |
| `kappa_s` | Real | required | Scattering opacity |
| `power_opacity` | bool | false | Power-law opacity |
| `compton` | bool | false | Compton scattering |
| `beam_source` | bool | false | Beam source |
| `affect_fluid` | bool | true | Radiation affects fluid |
| `fixed_fluid` | bool | false | Fixed background |
| `angular_fluxes` | bool | true | Angular flux transport |
| `reconstruct` | string | plm | Reconstruction method |
| `n_0_floor` | Real | 0.1 | Flux limiter floor |

## M1 Closure Relations

### Eddington Tensor
Closure function:
$$\chi = \frac{3 + 4f^2}{5 + 2\sqrt{4-3f^2}}$$

Where $f = |\mathbf{F}_r|/(cE_r)$ is the flux factor

### Pressure Tensor
$$P_{ij} = D_{ij} E_r$$
$$D_{ij} = \frac{1-\chi}{2} \delta_{ij} + \frac{3\chi-1}{2} n_i n_j$$

## Opacity Models

### Constant Opacity
$$\kappa = \kappa_0 = \text{const}$$

### Power-Law Opacity
$$\kappa = \kappa_0 \left(\frac{\rho}{\rho_0}\right)^a \left(\frac{T}{T_0}\right)^b$$

### Rosseland Mean
Pre-tabulated opacity tables supported

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

### Tetrad Basis
Orthonormal tetrad:
- $e^\mu_0$: Timelike (fluid frame)
- $e^\mu_1$: Spacelike (radial)
- $e^\mu_2$: Spacelike (theta)
- $e^\mu_3$: Spacelike (phi)

### Frame Transformations
Lab to fluid frame:
$$E_{\text{fluid}} = \Lambda^0_0 E_{\text{lab}} - \Lambda^0_i F^i_{\text{lab}}$$
$$F^i_{\text{fluid}} = \Lambda^i_0 E_{\text{lab}} - \Lambda^i_j F^j_{\text{lab}}$$

## Implicit Solver

### Newton-Raphson Iteration
```cpp
while (error > tolerance) {
  // Jacobian matrix
  // J = ∂R/∂U
  
  // Solve linear system
  // ΔU = -J^(-1) * R
  
  // Update
  // U_new = U_old + ΔU
}
```

Update equations:
$$\mathbf{J} = \frac{\partial \mathbf{R}}{\partial \mathbf{U}}$$
$$\Delta \mathbf{U} = -\mathbf{J}^{-1} \cdot \mathbf{R}$$
$$\mathbf{U}_{\text{new}} = \mathbf{U}_{\text{old}} + \Delta \mathbf{U}$$

## Source Terms

### Emission
$$S_{\text{emit}} = \kappa_a \cdot a \cdot T^4$$

### Absorption
$$S_{\text{abs}} = -\kappa_a \cdot E_r$$

### Scattering
$$S_{\text{scat}} = -\kappa_s \cdot \mathbf{F}_r/c$$

## Beam Sources

For testing and laser applications:
```cpp
// Collimated beam
F_beam = I_0 * exp(-r²/w²) * direction
```

## Variable Arrays

### Radiation Variables
- `rad.er`: Radiation energy density
- `rad.fr1`: Radiation flux x1
- `rad.fr2`: Radiation flux x2
- `rad.fr3`: Radiation flux x3

### Radiation Pressure (computed)
- `pr11, pr12, pr13`
- `pr22, pr23`
- `pr33`

## Common Test Problems

### Radiation Diffusion
```ini
<radiation>
arad = 1.0
kappa_a = 100.0  # Optically thick
kappa_s = 0.0
```

### Radiation Beam
```ini
<radiation>
beam_source = true
kappa_a = 0.1  # Optically thin
```

### Radiation Shock
```ini
<radiation>
arad = 1.0
kappa_a = 1.0
affect_fluid = true
```

## Performance Considerations

### Implicit Solver Cost
- Dominates in optically thick regime
- Use subcycling for efficiency
- Consider reduced speed of light

### Memory Usage
- 4 variables (er, fr1, fr2, fr3)
- Plus coupling terms
- Jacobian storage for implicit

## Common Issues

### Convergence Problems
- Reduce CFL for radiation
- Adjust Newton-Raphson tolerance
- Check opacity values

### Negative Energy
- Use flux limiters
- Check reconstruction
- Ensure proper floor values

## See Also
- [Hydro Module](hydro.md)
- Source: `src/radiation/radiation.cpp`