# Module: MHD

## Overview
The MHD module solves the ideal magnetohydrodynamics equations using constrained transport (CT) to maintain the divergence-free constraint $\nabla \cdot \mathbf{B} = 0$.

## Source Location
`src/mhd/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `mhd.hpp/cpp` | Core MHD class | Constructor, initialization |
| `mhd_fluxes.cpp` | Flux calculations | `CalculateFluxes()` |
| `mhd_ct.cpp` | Constrained transport | `CT()` |
| `mhd_corner_e.cpp` | Corner electric fields | `CornerE()` |
| `mhd_update.cpp` | Conservative update | `Update()` |
| `mhd_newdt.cpp` | Timestep calculation | `NewTimeStep()` |
| `mhd_tasks.cpp` | Task registration | Task management |
| `mhd_fofc.cpp` | First-order flux correction | FOFC implementation |

## Equations Solved

### Ideal MHD Equations (Conservative Form)
$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0$$
$$\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot \left(\rho \mathbf{v}\mathbf{v} - \mathbf{B}\mathbf{B} + \left(p + \frac{B^2}{8\pi}\right)\mathbb{I}\right) = 0$$
$$\frac{\partial E}{\partial t} + \nabla \cdot \left(\left(E + p + \frac{B^2}{8\pi}\right)\mathbf{v} - \frac{\mathbf{B}(\mathbf{B} \cdot \mathbf{v})}{4\pi}\right) = 0$$
$$\frac{\partial \mathbf{B}}{\partial t} - \nabla \times (\mathbf{v} \times \mathbf{B}) = 0$$

Where:
- $\rho$: density
- $\mathbf{v}$: velocity
- $\mathbf{B}$: magnetic field
- $p$: gas pressure
- $E$: total energy density = $\rho\varepsilon + \frac{\rho v^2}{2} + \frac{B^2}{8\pi}$

## Configuration Parameters

From `<mhd>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `eos` | string | required | Equation of state (ideal only) |
| `gamma` | Real | - | Adiabatic index |
| `reconstruct` | string | plm | Reconstruction (dc, plm, ppm4, ppmx, wenoz) |
| `rsolver` | string | required | Riemann solver |
| `fofc` | bool | false | First-order flux correction |
| `nscalars` | int | 0 | Number of passive scalars |
| `resistivity` | Real | - | Ohmic resistivity coefficient |
| `viscosity` | Real | - | Kinematic viscosity coefficient |
| `conductivity` | Real | - | Thermal conductivity coefficient |

## Riemann Solvers

### Non-Relativistic MHD
| Solver | Description | Waves | Speed | Accuracy |
|--------|-------------|-------|-------|----------|
| `llf` | Local Lax-Friedrichs | All | Fast | Diffusive |
| `hlle` | Harten-Lax-van Leer-Einfeldt | 2 | Medium | Good |
| `hlld` | HLLE with all discontinuities | 7 | Slow | Excellent |

### Relativistic MHD
- `llf_sr`: Special relativistic LLF
- `hlle_sr`: Special relativistic HLLE  
- `llf_gr`: General relativistic LLF
- `hlle_gr`: General relativistic HLLE

**Note**: Roe solver is not yet implemented for MHD (commented out in code).

## Reconstruction Methods

| Method | Order | Description | Ghost Zones |
|--------|-------|-------------|-------------|
| `dc` | 1st | Donor cell (piecewise constant) | 2 |
| `plm` | 2nd | Piecewise linear (minmod limiter) | 2 (3 with FOFC) |
| `ppm4` | 4th | Piecewise parabolic (4th order) | 3 (4 with FOFC) |
| `ppmx` | 4th | Extremum-preserving PPM | 3 (4 with FOFC) |
| `wenoz` | 5th | Weighted ENO-Z | 3 (4 with FOFC) |

## Constrained Transport (CT)

The CT algorithm maintains $\nabla \cdot \mathbf{B} = 0$ to machine precision:

### Key Features
1. **Face-centered B-fields**: Magnetic field components stored at cell faces
2. **Edge-centered E-fields**: Electric fields computed at cell edges
3. **Stokes theorem update**: Magnetic flux through faces updated via line integral

### CT Algorithm Flow
```cpp
// 1. Compute edge-centered electric fields
CornerE(bcc, e3x1f, e1x2f, e2x3f);

// 2. Update face-centered magnetic fields using Stokes theorem
// b1f_new = b1f_old - dt/dy*(e3_upper - e3_lower);
// b2f_new = b2f_old - dt/dz*(e1_upper - e1_lower);  
// b3f_new = b3f_old - dt/dx*(e2_upper - e2_lower);

// 3. Average to cell centers for primitive variables
// bcc = 0.5*(b_face_left + b_face_right);
```

Update equations:
$$B^{n+1}_{1,f} = B^n_{1,f} - \frac{\Delta t}{\Delta y}(E_{3,\text{upper}} - E_{3,\text{lower}})$$
$$B^{n+1}_{2,f} = B^n_{2,f} - \frac{\Delta t}{\Delta z}(E_{1,\text{upper}} - E_{1,\text{lower}})$$
$$B^{n+1}_{3,f} = B^n_{3,f} - \frac{\Delta t}{\Delta x}(E_{2,\text{upper}} - E_{2,\text{lower}})$$

## Variable Arrays

### Conservative Variables (`cons`)
- `cons(IDN)`: Density $\rho$
- `cons(IM1,IM2,IM3)`: Momentum $\rho\mathbf{v}$
- `cons(IEN)`: Total energy $E$
- Face-centered B-field arrays

### Primitive Variables (`prim`)
- `prim(IDN)`: Density $\rho$
- `prim(IVX,IVY,IVZ)`: Velocity $\mathbf{v}$
- `prim(IPR)`: Gas pressure $p$
- `prim(IBX,IBY,IBZ)`: Magnetic field $\mathbf{B}$ (cell-centered)

### Face-Centered Fields
- `b1f`: $B_1$ at x1-faces
- `b2f`: $B_2$ at x2-faces
- `b3f`: $B_3$ at x3-faces

## Execution Flow

```{mermaid}
flowchart TD
    Start[MHD Tasks] --> Recv[Start Boundary Recv]
    Recv --> FaceB[Face B to Cell B]
    FaceB --> P2C[Primitive to Conservative]
    P2C --> Recon[Reconstruction]
    Recon --> Riemann[Riemann Solver]
    Riemann --> EMF[Corner EMF]
    EMF --> CT[Constrained Transport]
    CT --> Wait[Wait for Boundaries]
    Wait --> Update[Conservative Update]
    Update --> C2P[Conservative to Primitive]
    C2P --> BC[Apply BCs]
    BC --> Done[Task Complete]
```

## AMR and div(B) Preservation

### Recent Improvements (2024)
The module includes fixes for maintaining $\nabla \cdot \mathbf{B} = 0$ with AMR:

1. **Prolongation**: Uses divergence-preserving prolongation operators
2. **Restriction**: Maintains magnetic flux conservation
3. **Flux correction**: Ensures consistency at refinement boundaries

From recent commits:
- Fixed div(B) issues during AMR derefinement
- Added comprehensive test suite (`tst/regression/mhd_amr_divb.py`)
- Verified with field loop advection tests

## First-Order Flux Correction (FOFC)

For strong shocks and to maintain positivity:
```cpp
if (fofc && shock_detected) {
  // Use first-order fluxes
  flux = flux_first_order;
}
```

## Common Test Problems

### 1D Tests
- **Brio-Wu shock tube**: Standard MHD Riemann problem
- **Ryu-Jones tests**: Suite of 1D MHD shocks

### 2D Tests  
- **Orszag-Tang vortex**: MHD turbulence and reconnection
- **Field loop advection**: Tests $\nabla \cdot \mathbf{B}$ preservation
- **MHD rotor**: Strong torsional Alfvén waves

### 3D Tests
- **MHD blast wave**: Spherical explosion in magnetized medium
- **Turbulence**: Driven MHD turbulence

## Performance Considerations

### GPU Optimization
- Face-centered fields require careful indexing
- CT electric field calculation vectorized
- MeshBlockPacks for coalesced memory access

### Typical Performance
- ~3-5×10⁶ zone-cycles/sec on V100 GPU
- ~80% efficiency compared to hydro

## Common Issues and Solutions

### Negative Pressure
- Reduce CFL number (typically 0.4 for MHD)
- Enable FOFC for strong shocks
- Check initial magnetic field strength

### $\nabla \cdot \mathbf{B}$ Growth
- Ensure using CT (always enabled in AthenaK)
- Check boundary conditions
- Verify AMR prolongation/restriction if using refinement

### Carbuncle Instability
- Use HLLD solver (most accurate)
- Add small resistivity
- Increase resolution

## See Also
- [Hydro Module](hydro.md) - Hydrodynamics solver
- [Reconstruction Module](reconstruction.md) - Reconstruction methods
- [Riemann Solvers Module](riemann_solvers.md) - Solver details
- Source: `src/mhd/`