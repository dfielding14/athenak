# Module: Coordinates

## Overview
The Coordinates module handles coordinate systems and metrics for Cartesian and general relativistic spacetimes. **Note: AthenaK currently only supports Cartesian coordinates and Kerr-Schild coordinates for black hole spacetimes. Spherical polar and cylindrical coordinates are NOT implemented.**

## Source Location
`src/coordinates/`

## Key Components

| File | Purpose | Key Classes/Functions |
|------|---------|----------------------|
| `coordinates.hpp` | Base coordinate class | `Coordinates` |
| `coordinates.cpp` | Coordinate implementations | `CoordSrcTerms` |
| `cartesian_ks.hpp` | Kerr-Schild coordinates | `ComputeMetricAndInverse` |
| `adm.hpp/cpp` | ADM formalism for GR | `ADM` class |
| `excision.cpp` | Black hole excision | `SetExcisionMasks` |
| `cell_locations.hpp` | Cell position utilities | `CellCenterX`, `FaceX` |

## Supported Coordinate Systems

### 1. Cartesian (Default)
- Standard flat space coordinates (x, y, z)
- Uniform grid spacing
- Used for all non-relativistic simulations
- No coordinate singularities

### 2. Kerr-Schild (General Relativity Only)
- Used for black hole spacetimes
- Horizon-penetrating coordinates
- Requires `general_rel = true` or dynamic GR (Z4c/ADM)
- Supports spinning black holes via spin parameter `a`

## Configuration Parameters

From `<coord>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `special_rel` | bool | false | Enable special relativity |
| `general_rel` | bool | false | Enable general relativity |
| `minkowski` | bool | false | Use flat Minkowski spacetime (GR mode) |
| `a` | Real | - | Black hole spin parameter (-1 < a < 1) |
| `excise` | bool | true | Enable excision inside black hole |
| `dexcise` | Real | - | Density floor in excision region |
| `pexcise` | Real | - | Pressure floor in excision region |
| `excision_scheme` | string | "fixed" | Excision method: "fixed" or "lapse" |
| `excise_lapse` | Real | 0.25 | Lapse threshold for excision (lapse method) |

## Cell Location Utilities

### Getting Cell Positions
```cpp
// From cell_locations.hpp
Real x1v = CellCenterX(i-is, nx1, x1min, x1max);  // Cell center in x1
Real x1f = FaceX(i-is, nx1, x1min, x1max);        // Cell face in x1
```

### Index Calculations
```cpp
// Linear spacing assumed
Real dx = (x1max - x1min) / nx1;
Real x_center = x1min + (i + 0.5) * dx;
Real x_face = x1min + i * dx;
```

## General Relativistic Coordinates

### Kerr-Schild Metric
For black holes, AthenaK uses Kerr-Schild coordinates which are regular at the horizon:

```cpp
// From cartesian_ks.hpp
ComputeMetricAndInverse(x1, x2, x3, is_minkowski, spin, glower, gupper);
```

The metric components include:
- Lapse function α
- Shift vector βⁱ
- Spatial metric γᵢⱼ
- Full 4-metric gμν

### Coordinate Source Terms
For relativistic hydrodynamics/MHD, geometric source terms are computed:

```cpp
// From coordinates.cpp
void CoordSrcTerms(prim, eos, dt, cons);
```

These account for:
- Christoffel symbols
- Metric derivatives
- Curved spacetime effects

## Excision for Black Holes

### Excision Region
Inside the black hole horizon, cells are excised (masked) from computation:

```cpp
// Check if inside excision region
if (excision_floor(m,k,j,i)) {
  // Reset to floor values
  prim(m,IDN,k,j,i) = dexcise;
  prim(m,IEN,k,j,i) = pexcise;
  prim(m,IVX,k,j,i) = 0.0;
  // ...
}
```

### Excision Methods
1. **Fixed**: Excise at fixed coordinate radius
2. **Lapse**: Excise where lapse α < threshold

## Usage Examples

### Basic Cartesian Setup
```ini
# Default - no coordinate block needed
<mesh>
nx1 = 256
x1min = -1.0
x1max = 1.0
# Automatically uses Cartesian coordinates
```

### Special Relativity
```ini
<coord>
special_rel = true  # Enables SR

<hydro>
gamma_max = 10.0   # Lorentz factor limit
```

### Kerr Black Hole
```ini
<coord>
general_rel = true
a = 0.9           # Dimensionless spin
excise = true
dexcise = 1.0e-8  # Floor density
pexcise = 1.0e-10 # Floor pressure
excision_scheme = lapse
excise_lapse = 0.2
```

### Dynamic GR (with Z4c)
```ini
<z4c>
# Enables dynamic spacetime evolution

<coord>
# general_rel automatically set by z4c
a = 0.7           # Initial black hole spin
excise = true
```

## Integration with Physics Modules

### Hydro/MHD
- Provides metric for covariant equations
- Computes geometric source terms
- Handles Lorentz factor calculations

### Z4c/ADM
- Uses ADM decomposition of spacetime
- Evolves metric components dynamically
- Coordinates module provides initial metric

## Important Limitations

### No Curvilinear Coordinates
- **Spherical polar coordinates NOT available**
- **Cylindrical coordinates NOT available**  
- All simulations use Cartesian grid
- For spherical problems, use Cartesian with appropriate domain

### Workarounds for Spherical Problems
```ini
# Simulate sphere in Cartesian coordinates
<mesh>
nx1 = 128
x1min = -1.0
x1max = 1.0
nx2 = 128  
x2min = -1.0
x2max = 1.0
nx3 = 128
x3min = -1.0
x3max = 1.0

# Apply spherical initial conditions in problem generator
```

## Common Issues

### Black Hole Excision
- Must have sufficient resolution near horizon
- Ghost zones must extend beyond excision region
- Use appropriate boundary conditions

### Relativistic Speeds
- Requires smaller CFL number (typically 0.2-0.4)
- Check Lorentz factor limits

## See Also
- [Mesh Module](mesh.md) - Grid structure
- [Z4c Module](z4c.md) - Numerical relativity
- [DynGRMHD Module](dyn_grmhd.md) - Relativistic MHD
- Source: `src/coordinates/`