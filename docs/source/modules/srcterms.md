# Module: Source Terms

## Overview
The Source Terms module implements external sources and forcing including turbulence driving, cooling, gravity, rotation, and user-defined source terms.

## Source Location
`src/srcterms/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `srcterms.hpp/cpp` | Core source terms class | Management, registration |
| `turb_driver.hpp/cpp` | Turbulence forcing | Fourier/SFB driving |
| `srcterms_newdt.cpp` | Timestep constraints | Source term CFL |
| `cooling_tables.hpp` | Cooling functions | Tabulated cooling |
| `ismcooling.hpp` | ISM cooling curves | Temperature-dependent |

## Turbulence Driving

### Configuration (`<turb_driving>` block)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `turb_flag` | int | 2 | 0=off, 1=initial, 2=continuous |
| `dedt` | Real | 0.0 | Energy injection rate |
| `tcorr` | Real | 0.0 | Correlation time |
| `rseed` | int | -1 | Random seed |
| `basis_type` | int | 0 | 0=Fourier, 1=SFB |
| `driving_type` | int | 0 | 0=solenoidal, 1=compressive |
| `nlow` | int | 1 | Min wavenumber |
| `nhigh` | int | 3 | Max wavenumber |

### Fourier Driving
Random phases in k-space:
$$f(\mathbf{k}) = A(\mathbf{k}) \cdot \exp(i\phi(\mathbf{k}))$$

Force in real space:
$$\mathbf{F}(\mathbf{x}) = \text{FFT}^{-1}[f(\mathbf{k})]$$

### Spherical Fourier-Bessel (SFB)
For spherical geometries:
$$f(r,\theta,\phi) = \sum_{nlm} A_{nlm} \cdot j_l\left(\frac{k_{ln} r}{r_0}\right) \cdot Y_{lm}(\theta,\phi)$$

where $k_{ln} = x_{ln}/r_0^{\text{turb}}$ ($x_{ln}$ = spherical Bessel roots)

#### SFB-specific parameters:
| Parameter | Type | Description |
|-----------|------|-------------|
| `lmax` | int | Maximum spherical harmonic degree |
| `nmax` | int | Maximum radial mode index |
| `r0_turb` | Real | Outer radius of driven region |

### Implementation
```cpp
// In turb_driver.cpp
void TurbulenceDriver::Generate() {
  // Generate random amplitudes
  // Apply projection (solenoidal/compressive)
  // Normalize to target power
}
```

## Cooling Functions

### ISM Cooling
Temperature-dependent cooling:
$$\Lambda(T) = \begin{cases}
\Lambda_0 \left(\frac{T}{T_0}\right)^\alpha & \text{for } T < T_{\text{peak}} \\
\Lambda_1 \left(\frac{T}{T_1}\right)^\beta & \text{for } T > T_{\text{peak}}
\end{cases}$$

### Tabulated Cooling
```cpp
// From cooling_tables.hpp
struct CoolingTable {
  Real GetCoolingRate(T, n);
  Real GetHeatingRate(T, n);
};
```

## Gravity

### Point Source
$$\mathbf{F}_{\text{grav}} = -\frac{GM}{r^2} \hat{\mathbf{r}}$$

### Uniform Field
$$\mathbf{F}_{\text{grav}} = -g \hat{\mathbf{z}}$$

### NFW Potential
$$\Phi = -\frac{GM}{r} \ln\left(1 + \frac{r}{r_s}\right)$$

## Rotation (Coriolis/Centrifugal)

### Rotating Frame
$$\mathbf{F}_{\text{rot}} = -2\boldsymbol{\Omega} \times \mathbf{v} - \boldsymbol{\Omega} \times (\boldsymbol{\Omega} \times \mathbf{r})$$

### Shearing Box
From `<shearing_box>` block:
$$\mathbf{F}_{\text{shear}} = 2\Omega v_y \hat{\mathbf{x}} - 2\Omega v_x \hat{\mathbf{y}}$$

## User-Defined Sources

**Note**: User-defined source terms are implemented directly in problem generators, not through a registration mechanism.

### Implementation
```cpp
// In problem generator, add source terms in the appropriate task
void MyProblem::AddSources(MeshBlockPack *pmbp, Real dt) {
  par_for("user_source", DevExeSpace(),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Add custom source terms directly
    cons(IEN,m,k,j,i) += dt * heating_rate;
    cons(IM1,m,k,j,i) += dt * force_x;
  });
}
```

## Timestep Constraints

### Cooling Time
$$\Delta t_{\text{cool}} = \frac{e}{|de/dt|_{\text{cooling}}}$$

### Orbital Time
$$\Delta t_{\text{orbit}} = \frac{2\pi}{\Omega}$$

## Execution Flow

```{mermaid}
flowchart TD
    Start[Source Terms] --> Check{Which Sources?}
    Check -->|Turbulence| Turb[Generate Forcing]
    Check -->|Cooling| Cool[Apply Cooling]
    Check -->|Gravity| Grav[Add Gravity]
    Check -->|User| User[User Sources]
    
    Turb --> Apply[Apply to Momentum]
    Cool --> Apply2[Apply to Energy]
    Grav --> Apply
    User --> Apply3[Apply to All]
    
    Apply --> Done[Update Complete]
    Apply2 --> Done
    Apply3 --> Done
```

## Common Applications

### Driven Turbulence
```ini
<turb_driving>
turb_flag = 2       # Continuous driving
dedt = 1.0          # Power input
tcorr = 0.5         # Correlation time
driving_type = 0    # Solenoidal
```

### Cooling Flow
```ini
<problem>
cooling_type = ism  # ISM cooling curve
T_floor = 1e4       # Temperature floor
```

### Gravitational Collapse
```ini
<problem>
grav_type = point_mass
GM = 1.0
softening = 0.1
```

## Performance Considerations

### Turbulence Driving
- FFT operations expensive
- Cache forcing patterns
- Update every N timesteps

### Cooling
- Table lookups can be slow
- Use interpolation
- Subcycle if needed

## Common Issues

### Energy Conservation
- Source terms can violate conservation
- Monitor total energy
- Use conservative formulation

### Stiff Cooling
- Cooling time << dynamical time
- Use implicit integration
- Apply temperature floor

### Force Imbalance
- Check force normalization
- Verify momentum conservation
- Monitor angular momentum

## Testing

Source term tests:
- Energy injection rate
- Cooling equilibrium
- Orbital dynamics

## See Also
- [Hydro Module](hydro.md)
- [MHD Module](mhd.md)
- [Shearing Box Module](shearing_box.md)
- Source: `src/srcterms/srcterms.cpp`