# Module: Riemann Solvers

## Overview
Riemann solvers compute numerical fluxes at cell interfaces by solving the Riemann problem between left and right states, essential for capturing shocks and discontinuities.

## Source Locations
- Hydro: `src/hydro/rsolvers/`
- MHD: `src/mhd/rsolvers/`
- Radiation: Built into module
- GRMHD: `src/dyn_grmhd/rsolvers/`

## Available Solvers

### Hydrodynamics

| Solver | File | Description | Properties |
|--------|------|-------------|------------|
| `llf` | `llf_hyd.hpp` | Local Lax-Friedrichs | Simple, diffusive |
| `hlle` | `hlle_hyd.hpp` | Harten-Lax-van Leer-Einfeldt | Robust, moderate diffusion |
| `hllc` | `hllc_hyd.hpp` | HLLE + Contact | Accurate, less diffusive |
| `roe` | `roe_hyd.hpp` | Roe's linearization | Accurate, can fail |
| `advect` | `advect_hyd.hpp` | Pure advection | Testing only |

### MHD

| Solver | File | Description | Properties |
|--------|------|-------------|------------|
| `llf` | `llf_mhd.hpp` | Local Lax-Friedrichs | Most robust |
| `hlle` | `hlle_mhd.hpp` | HLLE for MHD | Good balance |
| `hlld` | `hlld_mhd.hpp` | HLLE + Discontinuities | Most accurate |
| `advect` | `advect_mhd.hpp` | Pure advection | Testing only |

### Relativistic

| Solver | File | Description | Properties |
|--------|------|-------------|------------|
| `llf_sr` | `llf_srhyd.hpp` | Special relativistic LLF | Simple, robust |
| `hlle_sr` | `hlle_srhyd.hpp` | Special relativistic HLLE | Better accuracy |
| `hllc_sr` | `hllc_srhyd.hpp` | Special relativistic HLLC | Best accuracy |
| `llf_gr` | `llf_grhyd.hpp` | General relativistic LLF | For curved spacetime |
| `hlle_gr` | `hlle_grhyd.hpp` | General relativistic HLLE | GR with accuracy |

## Solver Details

### Local Lax-Friedrichs (LLF)
```cpp
// Maximum wave speed
λ_max = max(|v_L| + cs_L, |v_R| + cs_R)

// Flux
F = 0.5*(F_L + F_R - λ_max*(U_R - U_L))
```

**Pros:** Simple, always stable
**Cons:** Most diffusive

### HLLE
```cpp
// Wave speeds
S_L = min(v_L - cs_L, v_R - cs_R)
S_R = max(v_L + cs_L, v_R + cs_R)

// Flux
if (0 <= S_L) F = F_L
else if (S_R <= 0) F = F_R
else F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L))/(S_R - S_L)
```

**Pros:** Good balance of accuracy and robustness
**Cons:** Misses contact discontinuities

### HLLC (Hydro)
```cpp
// Additional middle wave
S_M = (pressure weighted average velocity)

// Four-state solution
if (0 <= S_L) F = F_L
else if (0 <= S_M) F = F*_L
else if (0 <= S_R) F = F*_R
else F = F_R
```

**Pros:** Captures all waves accurately
**Cons:** More complex, can fail in extreme conditions

### HLLD (MHD)
```cpp
// Five waves for MHD
// Fast, Alfven, Contact, Alfven, Fast
// Seven-state solution
```

**Pros:** Most accurate for MHD
**Cons:** Complex, expensive

### Roe Solver
```cpp
// Linearization around Roe average
Ũ = Roe_Average(U_L, U_R)

// Eigendecomposition
A(Ũ) = R Λ L

// Flux
F = 0.5*(F_L + F_R) - 0.5*R|Λ|L(U_R - U_L)
```

**Pros:** Very accurate for smooth flows
**Cons:** Can produce negative pressures, carbuncle instability

## Configuration

### Hydro
```ini
<hydro>
rsolver = hllc  # Recommended for most problems
```

### MHD
```ini
<mhd>
rsolver = hlld  # Best for MHD
# rsolver = llf  # If stability issues
```

### Relativistic
```ini
<hydro>
rsolver = hlle_sr  # Special relativity
# rsolver = hlle_gr  # General relativity
```

## Performance Comparison

| Solver | Speed | Accuracy | Robustness |
|--------|-------|----------|------------|
| LLF | Fastest | Low | Excellent |
| HLLE | Fast | Good | Very Good |
| HLLC | Medium | High | Good |
| HLLD | Slow | Highest | Good |
| Roe | Medium | High | Fair |

## Implementation Pattern

```cpp
template<typename T>
KOKKOS_INLINE_FUNCTION
void RiemannSolver(
    const T &prim_l,  // Left primitives
    const T &prim_r,  // Right primitives
    T &flux,          // Output flux
    const Real bx     // Normal B-field (MHD)
) {
    // 1. Compute wave speeds
    // 2. Evaluate flux formula
    // 3. Add source terms if needed
}
```

## Special Considerations

### Strong Shocks
- LLF or HLLE most robust
- HLLC/HLLD may need FOFC

### Contact Discontinuities
- HLLC/HLLD required for accuracy
- HLLE smears contacts

### Magnetic Reconnection
- HLLD recommended
- Maintains div(B)=0 better

### Relativistic Flows
- Use specialized SR/GR solvers
- Standard solvers fail for v→c

## Common Issues

### Negative Pressure
- Switch to more diffusive solver
- Enable first-order flux correction
- Check reconstruction method

### Carbuncle Instability
- Grid-aligned shock instability
- Use HLLE or add dissipation
- Avoid Roe solver

### Slow Convergence
- HLLD expensive for MHD
- Consider HLLE for testing
- Profile to identify bottleneck

## Testing

Standard tests for each solver:
- Sod shock tube
- Strong rarefaction
- Contact discontinuity
- Relativistic shocks

## See Also
- [Hydro Module](hydro.md)
- [MHD Module](mhd.md)
- [Reconstruction Module](reconstruction.md)
- Source: `src/hydro/rsolvers/`, `src/mhd/rsolvers/`