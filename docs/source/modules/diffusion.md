# Module: Diffusion

## Overview
The Diffusion module implements physical diffusion processes including viscosity, thermal conduction, and resistivity using explicit or implicit time integration schemes.

## Source Location
`src/diffusion/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `viscosity.hpp/cpp` | Kinematic/dynamic viscosity | Momentum diffusion |
| `conduction.hpp/cpp` | Thermal conduction | Heat diffusion |
| `resistivity.hpp/cpp` | Ohmic resistivity | Magnetic diffusion |
| `current_density.hpp` | Current calculation | $\mathbf{J} = \nabla \times \mathbf{B}$ |

## Configuration Parameters

### Viscosity
```ini
<mhd>  # or <hydro>
viscosity = 0.01  # Kinematic viscosity coefficient
```

### Thermal Conduction
```ini
<mhd>  # or <hydro>
conductivity = 0.1  # Thermal conductivity coefficient
cond_ceiling = 1.0e10  # Optional ceiling value
```

### Resistivity
```ini
<mhd>
ohmic_resistivity = 0.001   # Resistivity coefficient η
```

## Diffusion Equations

### Viscosity (Navier-Stokes)

The viscous force per unit mass (for constant kinematic viscosity $\nu$):

$$\frac{\partial \mathbf{v}}{\partial t} = ... + \nu \nabla^2 \mathbf{v} + \frac{\nu}{3} \nabla(\nabla \cdot \mathbf{v})$$

**Note on bulk viscosity**: AthenaK implements only shear (dynamic) viscosity with the standard assumption that bulk viscosity $\zeta = 0$. This is why the coefficient of the divergence term is $\frac{2\mu}{3}$ (or $\frac{2\nu}{3}$ in kinematic form) rather than $\zeta + \frac{2\mu}{3}$. This assumption is valid for monatomic gases and is commonly used in astrophysical simulations.

### Thermal Conduction

The energy equation with thermal conduction:

$$\frac{\partial E}{\partial t} + \nabla \cdot \mathbf{F}_E = ... - \nabla \cdot \mathbf{q}$$

where the heat flux is:

$$\mathbf{q} = -\kappa \nabla T$$

Note: In the code, $\kappa$ represents thermal **conductivity** (not diffusivity), with units of energy/(length·time·temperature).

### Resistivity (Magnetic Diffusion)

The induction equation with Ohmic resistivity:

$$\frac{\partial \mathbf{B}}{\partial t} = \nabla \times (\mathbf{v} \times \mathbf{B}) - \nabla \times (\eta \mathbf{J})$$

where the current density is:

$$\mathbf{J} = \nabla \times \mathbf{B}$$

For constant resistivity $\eta$, this simplifies to:

$$\frac{\partial \mathbf{B}}{\partial t} = \nabla \times (\mathbf{v} \times \mathbf{B}) + \eta \nabla^2 \mathbf{B}$$

## Numerical Implementation

### Heat Flux Implementation

At x1-faces, the code computes:

$$q_{x,i+1/2} = -\kappa \frac{T_i - T_{i-1}}{\Delta x}$$

where temperature is computed from:
- **Ideal gas**: $T = (\gamma - 1) e / \rho$ where $e$ is internal energy per unit mass
- **With ITM variable**: $T$ is stored directly

The energy flux is then updated:

$$F^E_{i+1/2} \mathrel{-}= q_{x,i+1/2}$$

### Resistive Electric Field Implementation

The code uses Ohm's law:

$$\mathbf{E} = -\mathbf{v} \times \mathbf{B} + \eta \mathbf{J}$$

where:
- The inductive term $-\mathbf{v} \times \mathbf{B}$ is computed in the MHD Riemann solver
- The resistive term $\eta \mathbf{J}$ is added separately

The current density components at cell edges are computed using:

$$J_x = \frac{\partial B_z}{\partial y} - \frac{\partial B_y}{\partial z}$$

$$J_y = \frac{\partial B_x}{\partial z} - \frac{\partial B_z}{\partial x}$$

$$J_z = \frac{\partial B_y}{\partial x} - \frac{\partial B_x}{\partial y}$$

Note: The code omits the $1/\mu_0$ factor by using normalized units where $\mu_0 = 1$.

### Stability Constraint
Diffusive CFL condition:

$$\Delta t < \frac{0.5 \Delta x^2}{D}$$

Where $D = \max(\nu, \kappa, \eta)$

## Viscous Stress Tensor

### Definition
The viscous stress tensor for a Newtonian fluid (with zero bulk viscosity):

$$\tau_{ij} = \mu\left(\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i} - \frac{2}{3}\delta_{ij} \nabla \cdot \mathbf{v}\right)$$

where:
- $\mu = \rho \nu$ is the dynamic (shear) viscosity
- $\nu$ is the kinematic viscosity (input parameter `viscosity`)
- The $-\frac{2}{3}\delta_{ij} \nabla \cdot \mathbf{v}$ term ensures the stress tensor is traceless (zero bulk viscosity)

### Viscous Flux Implementation
The code computes the viscous momentum flux directly. For example, at x1-faces:

$$F^{visc}_{x,1} = -\tau_{11} = -\mu\left(\frac{4}{3}\frac{\partial v_x}{\partial x} - \frac{2}{3}\frac{\partial v_y}{\partial y} - \frac{2}{3}\frac{\partial v_z}{\partial z}\right)$$

$$F^{visc}_{y,1} = -\tau_{12} = -\mu\left(\frac{\partial v_y}{\partial x} + \frac{\partial v_x}{\partial y}\right)$$

$$F^{visc}_{z,1} = -\tau_{13} = -\mu\left(\frac{\partial v_z}{\partial x} + \frac{\partial v_x}{\partial z}\right)$$


## Reynolds Numbers

### Kinematic Reynolds

$$\text{Re} = \frac{LV}{\nu}$$

### Magnetic Reynolds

$$\text{Rm} = \frac{LV}{\eta}$$

### Prandtl Number

$$\text{Pr} = \frac{\nu}{\kappa}$$

## Common Applications

### Viscous Shear Flow
```ini
<hydro>
viscosity = 0.1
```

### MHD Reconnection
```ini
<mhd>
ohmic_resistivity = 0.001
```

### Thermal Diffusion
```ini
<hydro>
conductivity = 0.01
```

## Performance Considerations

- Explicit diffusion limits timestep
- Consider implicit methods for large diffusion
- Diffusion terms add ~20% computational cost

## See Also
- [MHD Module](mhd.md)
- [Hydro Module](hydro.md)
- Source: `src/diffusion/`