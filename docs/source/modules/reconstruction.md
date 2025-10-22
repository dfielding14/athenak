# Module: Reconstruction

## Overview
The Reconstruction module provides spatial reconstruction methods for computing interface values from cell-centered data, crucial for high-order accuracy in finite volume schemes.

## Source Location
`src/reconstruct/`

## Reconstruction Methods

| File | Method | Order | Description |
|------|--------|-------|-------------|
| `dc.hpp` | DC | 1st | Piecewise constant (donor cell) |
| `plm.hpp` | PLM | 2nd | Piecewise linear with limiters |
| `ppm.hpp` | PPM4/PPMX | 4th/3rd | PPM4: 4th order accurate, PPMX: extrema preserving |
| `wenoz.hpp` | WENO-Z | 5th | Weighted ENO with Z-indicator |

## Method Details

### Donor Cell (DC)
```cpp
// First-order - most diffusive but robust
u_L = u_i
u_R = u_i
```

### Piecewise Linear (PLM)
```cpp
// Second-order with minmod limiter
δu = minmod(u_i - u_{i-1}, u_{i+1} - u_i)
u_L = u_i - 0.5*δu
u_R = u_i + 0.5*δu
```

### Piecewise Parabolic (PPM4/PPMX)
```cpp
// PPM4: 4th-order accurate PPM from Colella & Woodward
// PPMX: PPM with Colella-Sekora extrema preserving limiters
// Both construct parabola through cell
// Apply different limiters to prevent oscillations
```

### WENO-Z
```cpp
// Fifth-order weighted ENO
// Uses multiple stencils with smooth weights
// Z-indicator improves accuracy at discontinuities
```

## Slope Limiters

### MinMod
```cpp
minmod(a,b) = 0.5*(sign(a)+sign(b))*min(|a|,|b|)
```

### Van Leer
```cpp
vanleer(a,b) = 2ab/(a+b) if ab>0, else 0
```

### MC Limiter
```cpp
mc(a,b) = minmod(2a, 2b, 0.5*(a+b))
```

## Implementation

### Interface Template
```cpp
template<typename T>
KOKKOS_INLINE_FUNCTION
void ReconstructionX1(
    const T &q,      // Input array
    T &ql,          // Left states
    T &qr,          // Right states
    const int il,   // Start index
    const int iu    // End index
);
```

### Characteristic Variables
For systems of equations, reconstruction in characteristic variables improves accuracy:

```cpp
// Transform to characteristic variables
W = L * Q  // L = left eigenvectors

// Reconstruct W
Reconstruct(W, W_L, W_R)

// Transform back
Q_L = R * W_L  // R = right eigenvectors
Q_R = R * W_R
```

## Usage in Physics Modules

### Hydro
```ini
<hydro>
reconstruct = ppm4  # Options: dc, plm, ppm4, ppmx, wenoz
```

### MHD
```ini
<mhd>
reconstruct = plm  # Often more stable than ppm
```

### Radiation
```ini
<radiation>
reconstruct = plm  # Balance accuracy/robustness
```

## Performance Considerations

### Vectorization
- All methods vectorized over i-direction
- SIMD operations for slope calculations

### GPU Optimization
```cpp
// Reconstruction kernel
par_for("reconstruct", DevExeSpace(),
KOKKOS_LAMBDA(int k, int j, int i) {
  // Reconstruction at (i,j,k)
});
```

### Memory Access
- Stencil width affects cache usage
- DC: 1 cell, PLM: 3 cells, PPM4/PPMX: 5 cells, WENOZ: 6 cells

## Accuracy vs Robustness

| Method | Accuracy | Robustness | Cost |
|--------|----------|------------|------|
| DC | Low | Excellent | Minimal |
| PLM | Good | Good | Low |
| PPM4 | High | Moderate | Medium |
| PPMX | High | Good | Medium |
| WENOZ | Highest | Good | High |

## Common Issues

### Oscillations
- Use more diffusive limiter
- Reduce to lower order near shocks
- Check CFL number

### Carbuncle Instability
- Common with PPM4 in strong shocks
- Switch to PLM or PPMX (more robust)
- Enable first-order flux correction (FOFC)

## Testing

Convergence tests demonstrate order of accuracy:
- Smooth problems: Full order
- Discontinuous: Reduced to first order

## See Also
- [Hydro Module](hydro.md)
- [Riemann Solvers](riemann_solvers.md)
- Source: `src/reconstruct/`