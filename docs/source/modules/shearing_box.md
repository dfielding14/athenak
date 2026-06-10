# Module: Shearing Box

## Overview
The Shearing Box module implements the local shearing box approximation for studying disk dynamics, including orbital advection, Coriolis forces, and shearing-periodic boundary conditions.

## Source Location
`src/shearing_box/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `shearing_box.hpp/cpp` | Core shearing box | Setup, parameters |
| `orbital_advection.cpp` | FARGO algorithm | Orbital advection |
| `orbital_advection_cc.cpp` | Cell-centered advection | Scalar transport |
| `orbital_advection_fc.cpp` | Face-centered advection | B-field transport |
| `shearing_box_cc.cpp` | Cell-centered boundaries | Shearing-periodic |
| `shearing_box_fc.cpp` | Face-centered boundaries | B-field boundaries |
| `shearing_box_srcterms.cpp` | Source terms | Coriolis, tidal |
| `shearing_box_tasks.cpp` | Task management | Integration tasks |
| `remap_fluxes.hpp` | Flux remapping | Conservative remap |

## Configuration Parameters

From `<shearing_box>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `omega0` | Real | required | Orbital frequency $\Omega_0$ |
| `qshear` | Real | required | Shear parameter $q = -d\ln\Omega/d\ln r$ |

## Shearing Box Equations

### Momentum Equation
$$\frac{\partial \mathbf{v}}{\partial t} + \mathbf{v} \cdot \nabla \mathbf{v} = -\frac{\nabla P}{\rho} + 2q\Omega_0^2 x\mathbf{e}_x - 2\boldsymbol{\Omega}_0 \times \mathbf{v}$$

Where:
- First source: Tidal force
- Second source: Coriolis force

### Shear Flow
$$\mathbf{v}_0 = -q\Omega_0 x\mathbf{e}_y$$

## Orbital Advection (FARGO)

### Algorithm
Split velocity into background + perturbation:
$$\mathbf{v} = \mathbf{v}_0 + \delta\mathbf{v}$$

- Advect by background analytically
- Advect by perturbation numerically

### Implementation
Shift in y-direction:
$$y_{\text{shift}} = -q\Omega x \Delta t$$

Remap using high-order interpolation:
$$u_{\text{new}}(y) = u_{\text{old}}(y - y_{\text{shift}})$$

## Shearing-Periodic Boundaries

### Boundary Matching
At time $t$, the shift is:
$$\Delta y(t) = q\Omega L_x t$$

Periodic with shift:
$$u(x=0, y) = u(x=L_x, y-\Delta y)$$
$$u(x=L_x, y) = u(x=0, y+\Delta y)$$

### Radial Boundary Communication
1. Pack data at x-boundaries
2. Apply y-shift
3. Exchange via MPI
4. Unpack with interpolation

## Super Time Stepping

RKL2 super time stepping is compatible with shearing-periodic boundaries on
uniform grids. The pre-sweep uses the boundary displacement at time `t`,
while the post-sweep uses the displacement at `t + dt`.
After every STS stage, Hydro and MHD state is communicated and remapped
across the radial shearing boundaries.

Orbital advection remains in the ordinary hyperbolic step; it is not applied
inside each STS stage. Cell-centered diffusive fluxes on the two radial
faces are sheared and averaged before the STS divergence update. For MHD,
edge EMFs on the two radial faces are likewise sheared and averaged before
constrained transport updates the face-centered magnetic field, both in the
ordinary update and for Ohmic-resistive STS stages. The combined STS and
shearing-box path does not currently support SMR or AMR.

## Source Terms

### Coriolis Force
$$\mathbf{F}_{\text{Coriolis}} = 2\Omega_0(v_y\mathbf{e}_x - v_x\mathbf{e}_y)$$

### Tidal Force
$$\mathbf{F}_{\text{tidal}} = 2q\Omega_0^2 x\mathbf{e}_x$$

## Execution Flow

```{mermaid}
flowchart TD
    Start[Shearing Box] --> Orbital[Orbital Advection]
    Orbital --> Interior[Compute Interior]
    Interior --> Boundary[Shearing Boundaries]
    Boundary --> Remap[Remap with Shift]
    Remap --> Source[Add Source Terms]
    Source --> Coriolis[Coriolis Force]
    Source --> Tidal[Tidal Force]
    Coriolis --> Done[Complete]
    Tidal --> Done
```

## Common Applications

### MRI (Magnetorotational Instability)
```ini
<shearing_box>
omega0 = 1.0
qshear = 1.5  # Keplerian

<mhd>
# Vertical field for MRI
```

### Planetary Migration
```ini
<shearing_box>
omega0 = 1.0
qshear = 1.5

<problem>
# Add planet potential
```

### Disk Turbulence
```ini
<shearing_box>
omega0 = 1.0
qshear = 1.5

<turb_driving>
turb_flag = 2  # Drive turbulence
```

## Diagnostics

### Angular Momentum
$$\mathbf{L} = \int \rho(\mathbf{x} \times \mathbf{v}) dV$$

### Reynolds Stress
$$\alpha_{\text{Re}} = -\frac{\langle \rho v_x \delta v_y \rangle}{\langle P \rangle}$$

### Maxwell Stress
$$\alpha_{\text{Max}} = \frac{\langle B_x B_y \rangle}{4\pi \langle P \rangle}$$

## Implementation Details

### FARGO Speed-up
Only advect perturbations:
$$\text{CFL}_{\text{FARGO}} = \Delta t \cdot \frac{|\delta v|}{\Delta x} \ll \text{CFL}_{\text{standard}}$$

### Conservative Remapping
Ensure mass conservation during remap:
$$\int \rho_{\text{new}} \, dV = \int \rho_{\text{old}} \, dV$$

## Common Issues

### Radial Boundary Artifacts
- Check shearing-periodic implementation
- Verify interpolation order
- Test with linear shear wave

### Epicyclic Oscillations
- Natural at frequency 2Ω for q=1.5
- Not numerical artifact
- Physical disk response

### MRI Resolution
- Need ~10 cells per wavelength
- λ_MRI = 2π|vₐ|/Ω
- Under-resolved → decay

## Testing

Run the focused STS and shearing-box regressions from `tst/`:

```bash
python run_test_suite.py --mpicpu --test test_suite/sbox/test_sbox_sts_mpicpu.py
```

This compares STS against capped explicit integration, compares serial and
MPI decompositions with remote shearing partners, verifies passive-scalar
conservation and net magnetic fluxes, checks resistive `div B`, and exercises
a two-rank restart.

## See Also
- [MHD Module](mhd.md)
- [Source Terms Module](srcterms.md)
- Source: `src/shearing_box/shearing_box.cpp`
