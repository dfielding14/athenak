# Module: Equations of State (EOS)

## Overview
The EOS module provides thermodynamic closures for hydrodynamics and MHD, supporting ideal gas, isothermal, and tabulated equations of state for both Newtonian and relativistic regimes.

## Source Location
`src/eos/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `eos.hpp/cpp` | Base EOS class | Interface definitions |
| `ideal_hyd.cpp` | Ideal gas hydro | `ConsToPrim()`, `PrimToCons()` |
| `ideal_mhd.cpp` | Ideal gas MHD | Magnetic pressure included |
| `isothermal_hyd.cpp` | Isothermal hydro | Constant temperature |
| `isothermal_mhd.cpp` | Isothermal MHD | No energy equation |
| `ideal_srhyd.cpp` | Special relativistic | SR transformations |
| `ideal_grhyd.cpp` | General relativistic | GR with ideal gas |
| `primitive-solver/` | Advanced solvers | Tabulated EOS support |

## Supported EOS Types

### Ideal Gas
$$P = (\gamma - 1) \cdot e_{\text{int}}$$
$$e_{\text{int}} = E_{\text{total}} - \frac{1}{2}\rho v^2$$
$$T = \frac{P}{\rho R_{\text{gas}}}$$
$$c_s^2 = \frac{\gamma P}{\rho}$$

### Isothermal
$$P = c_s^2 \cdot \rho$$
$$T = \text{constant}$$
No energy equation needed

### Polytropic
$$P = K \cdot \rho^\gamma$$

### Tabulated (Nuclear)
From table lookup:
$$P = P(\rho, T, Y_e)$$
$$e = e(\rho, T, Y_e)$$

## Configuration

### Ideal Gas
```ini
<hydro>  # or <mhd>
eos = ideal
gamma = 1.4  # Adiabatic index
```

### Isothermal
```ini
<hydro>  # or <mhd>
eos = isothermal
iso_sound_speed = 1.0
```

## Conservative to Primitive Conversion

### Newtonian (ideal_hyd.cpp)
```cpp
void ConsToPrim(cons, prim) {
  // ρ = cons[IDN];
  // v = cons[IMX]/ρ;
  // E = cons[IEN];
  // P = (γ-1)*(E - 0.5*ρ*v²);
  
  // prim[IDN] = ρ;
  // prim[IVX] = v;
  // prim[IPR] = P;
}
```

Conversion equations:
$$\rho = \text{cons}[\text{IDN}]$$
$$\mathbf{v} = \frac{\text{cons}[\text{IMX}]}{\rho}$$
$$P = (\gamma - 1)\left(E - \frac{1}{2}\rho v^2\right)$$

### Relativistic (ideal_srhyd.cpp)
Newton-Raphson solver:
$$D = \gamma\rho$$ (Conserved density)
$$\mathbf{S} = \gamma^2\rho h\mathbf{v}$$ (Conserved momentum)
$$\tau = \gamma^2\rho h - P - D$$ (Conserved energy)

Solve for primitives iteratively

## Primitive Solver

### Noble Solver (2D)
2D Newton-Raphson in $(W, v^2)$:
$$W = \gamma^2\rho h$$ (Lorentz factor times enthalpy)

### Palenzuela Solver (3D)
```cpp
// Full 3D solver for GRMHD
// Robust for extreme conditions
```

## Sound Speed Calculations

### Ideal Gas
$$c_s^2 = \frac{\gamma P}{\rho}$$

### Isothermal
$$c_s^2 = \text{constant}$$ (input parameter)

### Relativistic
$$c_s^2 = \left(\frac{\partial P}{\partial \rho}\right)_s + \frac{P}{\rho^2}\left(\frac{\partial P}{\partial e}\right)_\rho$$

## Execution Flow

```{mermaid}
flowchart TD
    Cons[Conservative Vars] --> Check{EOS Type}
    Check -->|Ideal| IdealC2P[Algebraic C2P]
    Check -->|Isothermal| IsoC2P[Direct C2P]
    Check -->|Tabulated| TabC2P[Table Lookup]
    Check -->|Relativistic| RelC2P[Iterative C2P]
    
    IdealC2P --> Prim[Primitive Vars]
    IsoC2P --> Prim
    TabC2P --> Prim
    RelC2P --> Prim
    
    Prim --> Thermo[Thermodynamics]
    Thermo --> Sound[Sound Speed]
```

## Floor Values

### Density Floor
```cpp
// if (ρ < ρ_floor) ρ = ρ_floor;
```

Floor condition:
$$\rho = \max(\rho, \rho_{\text{floor}})$$

### Pressure Floor
```cpp
if (P < P_floor) P = P_floor;
```

### Temperature Floor
```cpp
if (T < T_floor) T = T_floor;
```

## Relativistic Considerations

### Maximum Lorentz Factor
```ini
<hydro>
gamma_max = 10.0  # Limit Lorentz factor
```

### Velocity Limit
```cpp
// if (v² > v_max²) {
//   v = v * v_max/|v|;
// }
```

Velocity limiting:
$$\mathbf{v} = \begin{cases}
\mathbf{v} & \text{if } v^2 \leq v_{\text{max}}^2 \\
\mathbf{v} \cdot \frac{v_{\text{max}}}{|\mathbf{v}|} & \text{if } v^2 > v_{\text{max}}^2
\end{cases}$$

## Tabulated EOS Support

### COMPOSE Format
```cpp
// eos_compose.cpp
// Nuclear EOS tables
struct ComposeEOS {
  Real GetPressure(ρ, T, Ye);
  Real GetEnergy(ρ, T, Ye);
  Real GetEntropy(ρ, T, Ye);
};
```

### Piecewise Polytrope
```cpp
// piecewise_polytrope.cpp
// Multiple polytropic segments
P = K_i * ρ^(Γ_i) for ρ_i < ρ < ρ_{i+1}
```

## Common Applications

### Shock Tubes
```ini
<hydro>
eos = ideal
gamma = 1.4  # Air
```

### Astrophysical Flows
```ini
<hydro>
eos = ideal
gamma = 5.0/3.0  # Monatomic gas
```

### Isothermal Turbulence
```ini
<mhd>
eos = isothermal
iso_sound_speed = 1.0
```

### Neutron Stars
```ini
<dyn_grmhd>
dyn_eos = compose  # Tabulated nuclear
table_path = /path/to/eos.h5
```

## Performance Optimization

### Caching
```cpp
// Cache frequently used values
struct EOSCache {
  Real cs2;  // Sound speed squared
  Real h;    // Specific enthalpy
  Real dpde; // Pressure derivative
};
```

### Vectorization
```cpp
#pragma omp simd
for (int i = is; i <= ie; ++i) {
  prim(IPR,k,j,i) = (gm1)*e_int(i);
}
```

## Common Issues

### Negative Pressure
- Check energy conservation
- Apply pressure floor
- Reduce CFL number

### C2P Convergence Failure
- Check initial guess
- Increase iteration limit
- Apply velocity ceiling

### Table Interpolation
- Ensure table bounds
- Check extrapolation
- Verify table resolution

## Testing

Standard tests:
- Sod shock tube (ideal gas)
- Isothermal shock (isothermal)
- Relativistic shock (SR/GR)

## See Also
- [Hydro Module](hydro.md)
- [MHD Module](mhd.md)
- Source: `src/eos/eos.cpp`