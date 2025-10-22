# Module: Hydrodynamics

## Overview
The Hydro module solves the Euler equations for inviscid fluid dynamics, supporting both Newtonian and relativistic hydrodynamics with various equations of state.

## Source Location
`src/hydro/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `hydro.hpp/cpp` | Core hydro class | Constructor, initialization |
| `hydro_fluxes.cpp` | Flux calculations | `CalculateFluxes()` |
| `hydro_update.cpp` | Conservative update | `Update()` |
| `hydro_newdt.cpp` | Timestep calculation | `NewTimeStep()` |
| `hydro_tasks.cpp` | Task registration | Task management |
| `hydro_fofc.cpp` | First-order flux correction | FOFC implementation |

## Equations Solved

### Euler Equations (Conservative Form)

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0$$
$$\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot (\rho \mathbf{v} \mathbf{v} + P\mathbb{I}) = 0$$
$$\frac{\partial E}{\partial t} + \nabla \cdot ((E+P)\mathbf{v}) = 0$$

Where:
- $\rho$: density
- $\mathbf{v}$: velocity
- $P$: pressure
- $E$: total energy density

## Configuration Parameters

From `<hydro>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `eos` | string | required | Equation of state (ideal, isothermal) |
| `reconstruct` | string | plm | Reconstruction (dc, plm, ppm4, ppmx, wenoz) |
| `rsolver` | string | required | Riemann solver |
| `gamma` | Real | - | Adiabatic index (ideal EOS) |
| `iso_sound_speed` | Real | - | Sound speed (isothermal) |
| `fofc` | bool | false | First-order flux correction |
| `nscalars` | int | 0 | Number of passive scalars |
| `viscosity` | Real | - | Kinematic viscosity coefficient |
| `conductivity` | Real | - | Thermal conductivity coefficient |

## Riemann Solvers

| Solver | Description | Speed | Accuracy |
|--------|-------------|-------|----------|
| `llf` | Local Lax-Friedrichs | Fast | Diffusive |
| `hlle` | Harten-Lax-van Leer-Einfeldt | Medium | Good |
| `hllc` | HLLE with Contact | Slow | Excellent |
| `roe` | Roe's approximate solver | Medium | Very good |
| `advect` | Pure advection (testing) | Fast | Limited |

### Relativistic Solvers
- `llf_sr`: Special relativistic LLF
- `hlle_sr`: Special relativistic HLLE
- `hllc_sr`: Special relativistic HLLC
- `llf_gr`: General relativistic LLF
- `hlle_gr`: General relativistic HLLE

## Reconstruction Methods

| Method | Order | Description | Ghost Zones |
|--------|-------|-------------|-------------|
| `dc` | 1st | Donor cell (piecewise constant) | 2 |
| `plm` | 2nd | Piecewise linear (minmod limiter) | 2 (3 with FOFC) |
| `ppm4` | 4th | Piecewise parabolic (4th order) | 3 (4 with FOFC) |
| `ppmx` | 4th | Extremum-preserving PPM | 3 (4 with FOFC) |
| `wenoz` | 5th | Weighted ENO-Z | 3 (4 with FOFC) |

## Execution Flow

```{mermaid}
flowchart TD
    Start[Hydro Tasks] --> Recv[Start Boundary Recv]
    Recv --> P2C[Primitive to Conservative]
    P2C --> Recon[Reconstruction]
    Recon --> Riemann[Riemann Solver]
    Riemann --> Flux[Calculate Fluxes]
    Flux --> Wait[Wait for Boundaries]
    Wait --> Update[Conservative Update]
    Update --> C2P[Conservative to Primitive]
    C2P --> BC[Apply BCs]
    BC --> Done[Task Complete]
```

## Variable Arrays

### Conservative Variables (`cons`)
- `cons(IDN)`: Density $\rho$
- `cons(IM1)`: Momentum $\rho v_1$
- `cons(IM2)`: Momentum $\rho v_2$
- `cons(IM3)`: Momentum $\rho v_3$
- `cons(IEN)`: Total energy $E$

### Primitive Variables (`prim`)
- `prim(IDN)`: Density $\rho$
- `prim(IV1)`: Velocity $v_1$
- `prim(IV2)`: Velocity $v_2$
- `prim(IV3)`: Velocity $v_3$
- `prim(IPR)`: Pressure $P$

## Key Methods

### Flux Calculation
```cpp
// hydro_fluxes.cpp L100-200
void Hydro::CalculateFluxes(MeshBlockPack *pmbp) {
  // Reconstruct primitives at faces
  Reconstruction(prim_left, prim_right);
  
  // Solve Riemann problem
  RiemannSolver(prim_left, prim_right, flux);
  
  // Store fluxes for update
  StoreFluxes(flux);
}
```

### Conservative Update
```cpp
// hydro_update.cpp L50-100
void Hydro::Update(MeshBlockPack *pmbp, Real dt) {
  // Update conserved variables
  cons_new = cons_old - dt/dx * (flux_right - flux_left);
}
```

### Timestep Calculation
```cpp
// hydro_newdt.cpp L30-80
Real Hydro::NewTimeStep() {
  // CFL condition
  // dt = cfl_number * min(dx/(|v|+cs));
}
```

CFL condition:

$$\Delta t = \text{CFL} \cdot \min\left(\frac{\Delta x}{|v| + c_s}\right)$$

## Equations of State

### Ideal Gas

$$P = (\gamma - 1)\left(E - \frac{1}{2}\rho v^2\right)$$
$$c_s = \sqrt{\frac{\gamma P}{\rho}}$$
```

### Isothermal
```cpp
P = cs² * ρ
// No energy equation needed
```

## Boundary Conditions

Standard boundary types:
- `reflecting`: Zero normal velocity
- `outflow`: Zero gradient
- `periodic`: Periodic wrap
- `user`: User-defined in problem generator

## First-Order Flux Correction (FOFC)

Ensures robustness near shocks:
```cpp
if (fofc && shock_detected) {
  // Use first-order flux
  flux = flux_first_order;
}
```

## Performance Optimization

### Vectorization
- Reconstruction vectorized over i-direction
- Riemann solver operates on pencils

### GPU Optimization
```cpp
par_for("hydro_fluxes", DevExeSpace(),
        0, nmb-1, ks, ke, js, je, is, ie+1,
KOKKOS_LAMBDA(int m, int k, int j, int i) {
  // Flux kernel
});
```

## Common Issues and Solutions

### Negative Pressure
- Check CFL number (reduce if needed)
- Enable FOFC for strong shocks
- Check initial conditions

### Carbuncle Instability
- Use HLLC or Roe solver
- Add small viscosity
- Refine grid

## Usage Examples

### Shock Tube
```ini
<hydro>
eos = ideal
gamma = 1.4
reconstruct = ppm
rsolver = hllc
```

### Isothermal Turbulence
```ini
<hydro>
eos = isothermal
iso_sound_speed = 1.0
reconstruct = plm
rsolver = hlle
```

### Relativistic Flow
```ini
<hydro>
eos = ideal
gamma = 4.0/3.0
rsolver = hlle_sr  # Special relativistic
```

## Testing

Standard test problems in `src/pgen/tests/`:
- `linear_wave.cpp`: Linear wave convergence
- `shock_tube.cpp`: Sod shock tube
- `lw_implode.cpp`: 2D implosion test

## See Also
- [MHD Module](mhd.md) - Magnetohydrodynamics
- [EOS Module](eos.md) - Equations of state
- [Riemann Solvers](riemann_solvers.md)
- Source: `src/hydro/hydro.cpp`