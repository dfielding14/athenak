# Module: DynGRMHD (Dynamical General Relativistic MHD)

## Overview
The DynGRMHD module couples magnetohydrodynamics with dynamical spacetimes from the Z4c module, enabling simulations of neutron star mergers, black hole accretion, and other strong-gravity phenomena.

## Source Location
`src/dyn_grmhd/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `dyn_grmhd.hpp/cpp` | Core GRMHD class | System initialization |
| `dyn_grmhd_fluxes.cpp` | GRMHD flux calculation | `CalculateFluxes()` |
| `dyn_grmhd_fofc.cpp` | First-order flux correction | Robustness near horizons |
| `dyn_grmhd_util.hpp` | Utility functions | Metric operations |

## Evolution Variables

### Conservative Variables
$$D = \sqrt{\gamma} \cdot \rho \cdot W$$ (Conserved density)
$$S_i = \sqrt{\gamma} \cdot \rho h^* \cdot W^2 v_i$$ (Conserved momentum)
$$\tau = \sqrt{\gamma} \cdot (\rho h^* W^2 - P - D)$$ (Conserved energy)
$$B^i = \sqrt{\gamma} \cdot B^i$$ (Magnetic field)

Where:
- $W$: Lorentz factor
- $h^*$: Specific enthalpy including magnetic
- $\gamma$: Determinant of spatial metric

## Configuration Parameters

From `<mhd>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dyn_eos` | string | required | EOS type (ideal, compose) |
| `dyn_error` | string | required | Error handling policy |
| `rsolver` | string | required | Riemann solver (llf, hlle) |
| `fofc` | bool | false | First-order flux correction |
| `fofc_method` | string | llf | FOFC Riemann solver |
| `dmp_M` | Real | 1.2 | Density floor multiplier |
| `enforce_maximum` | bool | true | Enforce ceilings |
| `fixed` | bool | false | Fixed metric (Cowling) |

## Conservative to Primitive Conversion

### Valencia Formulation
Iterative solver for primitives:
- Given: $D$, $S_i$, $\tau$, $B^i$
- Find: $\rho$, $v^i$, $P$

2D Newton-Raphson in $(W, P)$:
$$f_1 = S^2 - (D+\tau+P)^2W^2 + P^2$$
$$f_2 = \tau + D - \rho h W^2 + P$$

Iterate until convergence

### Recovery Policies
- `reset_floor`: Reset to floor values
- `reset_prim`: Keep previous primitives
- `kill_zone`: Set to vacuum

## Metric Coupling

### From Z4c Module
- Lapse function: $\alpha = \text{z4c.alpha}$
- Shift vector: $\beta^i = \text{z4c.beta}^i$
- Spatial metric: $\gamma_{ij} = \text{z4c.gamma}_{ij}$
- Extrinsic curvature: $K_{ij} = \text{z4c.K}_{ij}$

### Stress-Energy Tensor
Feedback to spacetime:
$$T^{\mu\nu} = (\rho h + b^2)u^\mu u^\nu + Pg^{\mu\nu} - b^\mu b^\nu$$

Source for Z4c:
$$S_{\text{Z4c}} = -8\pi\sqrt{\gamma} T^{\mu\nu}$$

## Execution Flow

```{mermaid}
flowchart TD
    Start[DynGRMHD Tasks] --> Metric[Get Metric from Z4c]
    Metric --> C2P[Conservative to Primitive]
    C2P --> Check{C2P Success?}
    Check -->|Yes| Flux[Calculate Fluxes]
    Check -->|No| Recovery[Apply Recovery]
    Recovery --> Flux
    Flux --> FOFC{FOFC Needed?}
    FOFC -->|Yes| FirstOrder[First-Order Flux]
    FOFC -->|No| Update[Conservative Update]
    FirstOrder --> Update
    Update --> Source[Add Source Terms]
    Source --> Tmunu[Calculate T^μν]
    Tmunu --> Feedback[Feedback to Z4c]
```

## Riemann Solvers

### Local Lax-Friedrichs (LLF)
Simple, robust:
$$\mathbf{F} = \frac{1}{2}(\mathbf{F}_L + \mathbf{F}_R - \lambda_{\text{max}}(\mathbf{U}_R - \mathbf{U}_L))$$

### HLLE
Better accuracy:
$$\mathbf{F} = \begin{cases}
\mathbf{F}_L & \text{if } 0 < S_L \\
\mathbf{F}_R & \text{if } S_R < 0 \\
\frac{S_R\mathbf{F}_L - S_L\mathbf{F}_R + S_LS_R(\mathbf{U}_R-\mathbf{U}_L)}{S_R-S_L} & \text{otherwise}
\end{cases}$$

## First-Order Flux Correction

### Shock Detection
Detect troubled cells:
$$\text{use\_fofc} = \begin{cases}
\text{true} & \text{if } \frac{|\nabla P|}{P} > \text{threshold} \text{ or } W > W_{\text{max}} \\
\text{false} & \text{otherwise}
\end{cases}$$

### FOFC Application
```cpp
if (use_fofc) {
  // Use first-order flux
  flux = LLF_flux(U_L, U_R);
}
```

## Floor and Ceiling Values

### Atmosphere Treatment
```cpp
ρ_atm = ρ_floor * r^(-power)
P_atm = P_floor * r^(-power)
```

### Velocity Ceiling
```cpp
if (W > W_max) {
  // Reduce velocity
  v = v * sqrt(1 - 1/W_max²)/|v|;
}
```

## Common Applications

### Binary Neutron Star Merger
```ini
<mhd>
dyn_eos = compose  # Nuclear EOS
rsolver = hlle
fofc = true
dmp_M = 1.5
```

### Black Hole Accretion
```ini
<mhd>
dyn_eos = ideal
gamma = 4.0/3.0
rsolver = llf
enforce_maximum = true
```

### Magnetar Formation
```ini
<mhd>
dyn_eos = compose
fixed = false  # Dynamical spacetime
fofc = true
```

## Performance Considerations

### C2P Iterations
- Most expensive operation
- Cache converged values
- Use good initial guesses

### GPU Optimization
```cpp
// All kernels GPU-ready
par_for("dyngrmhd_c2p", DevExeSpace(),
KOKKOS_LAMBDA(int m, int k, int j, int i) {
  // C2P kernel
});
```

## Common Issues

### C2P Failures
- Near black hole horizons
- In low-density regions
- At refinement boundaries

**Solutions:**
- Apply appropriate floors
- Use recovery policies
- Enable FOFC

### Constraint Violations
- Monitor Hamiltonian constraint
- Check momentum constraint
- Verify div(B) = 0

### Instabilities
- Reduce CFL factor
- Increase resolution
- Check boundary conditions

## Testing

Test problems:
- TOV star (equilibrium)
- BNS merger (dynamics)
- Black hole accretion

## See Also
- [Z4c Module](z4c.md) - Spacetime evolution
- [MHD Module](mhd.md) - Magnetohydrodynamics
- [EOS Module](eos.md) - Equations of state
- Source: `src/dyn_grmhd/dyn_grmhd.cpp`