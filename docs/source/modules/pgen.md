# Module: Problem Generators

## Overview
Problem Generators (pgen) define initial conditions and problem-specific behavior for simulations. Each problem generator is a self-contained setup for a specific physical scenario.

## Source Location
`src/pgen/`

## Structure

### Categories
- **Test Problems**: `tests/` subdirectory
- **Astrophysical**: Main directory
- **Benchmarks**: Standard problems

## Key Components

Every problem generator must define:
```cpp
void ProblemGenerator(ParameterInput *pin, const bool restart) {
  // Initialize simulation
  // Set initial conditions
  // Register user functions
}
```

## Standard Test Problems

### Hydrodynamics Tests (`tests/`)

| File | Test | Purpose |
|------|------|---------|
| `linear_wave.cpp` | Linear waves | Convergence testing |
| `shock_tube.cpp` | Sod shock tube | Riemann solver validation |
| `lw_implode.cpp` | 2D implosion | Multi-dimensional |
| `advection.cpp` | Passive advection | Transport accuracy |

### MHD Tests

| File | Test | Purpose |
|------|------|---------|
| `cpaw.cpp` | Circularly polarized Alfven | Wave propagation |
| `orszag_tang.cpp` | Orszag-Tang vortex | MHD turbulence |
| `field_loop.cpp` | Field loop advection | div(B) preservation |

### Relativity Tests

| File | Test | Purpose |
|------|------|---------|
| `gr_bondi.cpp` | Bondi accretion | GR hydro |
| `gr_monopole.cpp` | Magnetic monopole | GRMHD |
| `z4c_linear_wave.cpp` | GW propagation | NR evolution |

## Astrophysical Problems

### Accretion Disks

| File | Description | Physics |
|------|-------------|---------|
| `gr_torus.cpp` | Relativistic torus | GRMHD equilibrium |
| `mri2d.cpp` | 2D MRI | Magnetorotational instability |
| `mri3d.cpp` | 3D MRI | Full MRI turbulence |

### Stellar Problems

| File | Description | Physics |
|------|-------------|---------|
| `dyngr_tov.cpp` | TOV star | Neutron star equilibrium |
| `lorene_bns.cpp` | Binary NS (LORENE) | Initial data import |
| `sgrid_bns.cpp` | Binary NS (SGRID) | Initial data import |

### ISM/CGM Problems

| File | Description | Physics |
|------|-------------|---------|
| `turb.cpp` | Driven turbulence | Supersonic turbulence |
| `sfb_turb.cpp` | SFB turbulence | Spherical forcing |
| `blast.cpp` | Blast waves | Explosions |
| `rt.cpp` | Rayleigh-Taylor | Instabilities |
| `kh.cpp` | Kelvin-Helmholtz | Shear instabilities |

### Shocks

| File | Description | Physics |
|------|-------------|---------|
| `shock_cloud.cpp` | Shock-cloud | Cloud crushing |
| `shu_osher.cpp` | Shu-Osher | High Mach shock |
| `cshock.cpp` | C-type shock | Ion-neutral |

## Problem Generator Structure

### Basic Template
```cpp
namespace pgen {

void ProblemGenerator(ParameterInput *pin, const bool restart) {
  // 1. Read parameters
  Real param = pin->GetReal("problem", "param");
  
  // 2. Get mesh pointer
  auto pm = pmy_mesh_;
  
  // 3. Initialize over all MeshBlocks
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  
  // 4. Set initial conditions
  par_for("pgen", DevExeSpace(), 0, nmb-1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Set primitive variables
    prim(IDN,m,k,j,i) = density;
    prim(IVX,m,k,j,i) = velocity_x;
    prim(IPR,m,k,j,i) = pressure;
  });
  
  // 5. Register user functions (optional)
  RegisterUserHistoryFunctions();
  RegisterUserSourceTerms();
  RegisterUserBoundaryConditions();
}

} // namespace pgen
```

### Reading Parameters
```cpp
// From <problem> block
Real x0 = pin->GetReal("problem", "x0");
int n = pin->GetInteger("problem", "n");
std::string type = pin->GetString("problem", "type");
bool flag = pin->GetBoolean("problem", "flag");
```

### Coordinate Access
```cpp
// Get coordinates
Real x1 = pco->x1v(i);  // x or r
Real x2 = pco->x2v(j);  // y or theta
Real x3 = pco->x3v(k);  // z or phi
```

## User-Defined Functions

### History Functions
```cpp
Real UserHistory(MeshBlockPack *pmbp) {
  // Calculate diagnostic quantity
  return total_energy;
}

// Register in ProblemGenerator
pmbp->pmesh->RegisterUserHistoryFunction(UserHistory);
```

### Source Terms
```cpp
void UserSource(MeshBlockPack *pmbp, Real dt) {
  // Add source terms
}

// Register
pmbp->RegisterUserSourceTerm(UserSource);
```

### Boundary Conditions
```cpp
void UserBoundary(MeshBlockPack *pmbp) {
  // Apply custom BCs
}

// Register
pmbp->RegisterUserBoundary(INNER_X1, UserBoundary);
```

## Building with Custom Problem

### CMake Configuration
```bash
cmake -B build -DPROBLEM=my_problem
```

### File Naming
- File: `src/pgen/my_problem.cpp`
- Build flag: `-DPROBLEM=my_problem`

## Common Patterns

### Random Perturbations
```cpp
// Initialize random number generator
Random rng(rseed);

// Add perturbations
prim(IDN,m,k,j,i) = rho0 * (1.0 + amp*rng.ran());
```

### Reading External Data
```cpp
// Read profile from file
std::ifstream file(filename);
std::vector<Real> data;
// ... read data ...

// Interpolate to grid
Real value = interpolate(data, x);
```

### Wave Initialization
```cpp
// Sinusoidal wave
Real phase = kx*x + ky*y + kz*z;
prim(IDN,m,k,j,i) = rho0 + amp*sin(phase);
```

## Testing Problem Generators

### Convergence Test
```bash
# Run at different resolutions
./athena -i test_N32.athinput
./athena -i test_N64.athinput
./athena -i test_N128.athinput

# Check convergence rate
```

### Restart Test
```bash
# Run with restart
./athena -i test.athinput
./athena -r test.00100.rst
```

## Common Issues

### Incorrect Initial Conditions
- Check coordinate system
- Verify unit conversions
- Test with simple cases

### Boundary Artifacts
- Check BC compatibility
- Verify ghost zone filling
- Test uniform flow

### Performance
- Use Kokkos parallel_for
- Minimize serial sections
- Profile initialization

## See Also
- [Configuration Guide](../configuration.md)
- [Input Parameters](../reference/input_parameters.md)
- Source: `src/pgen/pgen.cpp`