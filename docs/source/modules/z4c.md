# Module: Z4c (Numerical Relativity)

## Overview
The Z4c module evolves Einstein's field equations using the Z4 formulation with constraint damping, enabling simulations of dynamical spacetimes including black holes, neutron stars, and gravitational waves.

## Source Location
`src/z4c/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `z4c.hpp/cpp` | Core Z4c evolution | System initialization |
| `z4c_calcrhs.cpp` | Right-hand side calculation | `CalculateRHS()` |
| `z4c_adm.cpp` | ADM decomposition | 3+1 split |
| `z4c_gauge.cpp` | Gauge conditions | Lapse, shift evolution |
| `z4c_amr.cpp` | AMR criteria | Refinement triggers |
| `z4c_wave_extr.cpp` | Wave extraction | GW signals |
| `z4c_Sbc.cpp` | Sommerfeld BCs | Outgoing waves |
| `tmunu.hpp/cpp` | Stress-energy tensor | Matter coupling |
| `compact_object_tracker.cpp` | BH/NS tracking | Apparent horizons |

## Evolution Variables

### Z4c Variables (19 total)
```cpp
// Conformal metric (6 components)
gxx, gxy, gxz, gyy, gyz, gzz

// Extrinsic curvature (6 components)
kxx, kxy, kxz, kyy, kyz, kzz

// Constraint variable
theta

// Conformal connection (3 components)
gamx, gamy, gamz

// Lapse and shift (4 components)
alpha, betax, betay, betaz
```

## Configuration Parameters

From `<z4c>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_z4c` | bool | true | Use Z4c (vs BSSN) |
| `chi_psi_power` | Real | -4.0 | Conformal factor power |
| `damp_kappa1` | Real | 0.0 | Constraint damping $\kappa_1$ |
| `damp_kappa2` | Real | 0.0 | Constraint damping $\kappa_2$ |
| `lapse_oplog` | Real | 2.0 | 1+log slicing power |
| `shift_eta` | Real | 2.0 | Shift damping $\eta$ |
| `shift_Gamma` | Real | 1.0 | Gamma-driver strength |
| `diss` | Real | 0.0 | Kreiss-Oliger dissipation |
| `chi_div_floor` | Real | -1000.0 | $\chi$ divergence floor |
| `eps_floor` | Real | 1e-12 | Small number floor |

## Evolution Equations

### Z4c System
$$\partial_t \chi = \frac{2}{3} \chi(\alpha K - \partial_i \beta^i)$$
$$\partial_t \tilde{g}_{ij} = -2\alpha \tilde{K}_{ij} + \mathcal{L}_\beta \tilde{g}_{ij}$$
$$\partial_t \tilde{K}_{ij} = \chi(\tilde{R}_{ij} + \nabla_i\nabla_j\alpha) - 2\tilde{K}_{ik}\tilde{K}^k_j + \alpha K \tilde{K}_{ij} + \mathcal{L}_\beta \tilde{K}_{ij}$$
$$\partial_t \Theta = \alpha(R - K_{ij}K^{ij} + K^2) - \kappa_1(2+\kappa_2)\alpha\Theta$$
$$\partial_t \tilde{\Gamma}^i = 2\alpha(\tilde{\Gamma}^i_{jk}\tilde{K}^{jk} - \tilde{g}^{ij}\partial_j K) + \mathcal{L}_\beta \tilde{\Gamma}^i - \kappa_1\alpha\tilde{\Gamma}^i$$

## Gauge Conditions

### Lapse Evolution (1+log)
$$\partial_t \alpha = -2\alpha K \left(\frac{\alpha}{\alpha_0}\right)^{\text{oplog}}$$

### Shift Evolution (Gamma-driver)
$$\partial_t \beta^i = \frac{3}{4} B^i$$
$$\partial_t B^i = \partial_t \tilde{\Gamma}^i - \eta B^i$$

## Constraint Monitoring

### Hamiltonian Constraint
$$\mathcal{H} = R - K_{ij}K^{ij} + K^2 - 16\pi \rho$$

### Momentum Constraint
$$\mathcal{M}^i = \nabla_j K^{ij} - \partial^i K - 8\pi S^i$$

## Implementation Flow

```{mermaid}
flowchart TD
    Init[Initialize Metric] --> ADM[ADM Decomposition]
    ADM --> Conformal[Conformal Transform]
    Conformal --> RHS[Calculate RHS]
    RHS --> Deriv[Spatial Derivatives]
    Deriv --> Ricci[Ricci Tensor]
    Ricci --> Source[Add Sources]
    Source --> Damp[Constraint Damping]
    Damp --> Update[RK Update]
    Update --> Gauge[Update Gauge]
    Gauge --> BC[Boundary Conditions]
    BC --> Extract[Wave Extraction]
```

## Wave Extraction

### Newman-Penrose Scalars
Weyl scalar $\Psi_4$:
$$\Psi_4 = \text{Weyl\_contraction}(\text{tetrad})$$

Gravitational wave strain:
$$h_+ - ih_\times = \int\int \Psi_4 \, dt \, dt$$

### Extraction Spheres
```ini
<z4c>
nrad_wave_extraction = 5    # Number of radii
extraction_nlev = 10        # Angular resolution
waveform_dt = 1.0          # Output frequency
```

## AMR Criteria

### Truncation Error
From z4c_amr.cpp:
$$\text{error} = |\nabla^2 \chi| \cdot \Delta x^2$$
$$\text{if } (\text{error} > \text{threshold}) \text{ then refine}$$

### Configuration
```ini
<z4c_amr>
method = chi_gradient
chi_min = 0.2        # Refine if $\chi$ < threshold
dchi_max = 0.1       # Refine if $|\nabla\chi|$ > threshold
```

## Black Hole Tracking

### Apparent Horizon Finding
```cpp
// compact_object_tracker.cpp
// Tracks BH properties:
- Mass (irreducible, ADM)
- Spin (magnitude, direction)
- Position
- Horizon area
```

## Common Initial Data

### Single Puncture
```cpp
// z4c_one_puncture.cpp
punc_ADM_mass = 1.0
punc_center_x1 = 0.0
```

### Binary Black Holes
```cpp
// z4c_two_puncture.cpp
// Uses TwoPunctures library
```

### Neutron Stars
```cpp
// lorene_bns.cpp, sgrid_bns.cpp
// Import from LORENE/SGRID
```

## Boundary Conditions

### Sommerfeld (Outgoing)
$$\partial_t f + v^i \partial_i f + \frac{f-f_0}{r} = 0$$

### Constraint-Preserving
Custom BCs that preserve constraints

## Performance Optimization

### Vectorization
- Derivative stencils vectorized
- Ricci tensor computation optimized

### GPU Considerations
```cpp
// All kernels use Kokkos
par_for("z4c_rhs", DevExeSpace(), ...
KOKKOS_LAMBDA(int m, int k, int j, int i) {
  // RHS kernel
});
```

## Common Issues

### Constraint Violations
- Increase constraint damping (κ₁, κ₂)
- Check boundary conditions
- Verify initial data

### Gauge Instabilities
- Adjust η for shift damping
- Modify lapse condition
- Add gauge driver

### AMR Issues
- Prolongation can violate constraints
- Use conservative variables
- Add buffer zones

## Testing

Standard tests in `src/pgen/tests/`:
- `z4c_linear_wave.cpp`: GW convergence
- `gr_monopole.cpp`: Trumpet solution

## Usage Examples

### Binary Black Hole
```ini
<z4c>
damp_kappa1 = 0.1
damp_kappa2 = 0.0
lapse_oplog = 2.0
shift_eta = 1.0

<z4c_amr>
method = chi_gradient
chi_min = 0.3
```

### Neutron Star
```ini
<z4c>
use_z4c = true
chi_psi_power = -4.0
damp_kappa1 = 0.02
shift_eta = 2.0
```

## See Also
- [DynGRMHD Module](dyn_grmhd.md) - GRMHD coupling
- [Coordinates Module](coordinates.md) - Metric components
- Source: `src/z4c/z4c.cpp`