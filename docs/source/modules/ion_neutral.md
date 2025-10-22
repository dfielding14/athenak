# Module: Ion-Neutral

## Overview
The Ion-Neutral module implements two-fluid MHD for partially ionized plasmas, including drag forces, ionization/recombination, and separate evolution of ion and neutral fluids.

## Source Location
`src/ion-neutral/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `ion-neutral.hpp/cpp` | Core two-fluid class | System setup |
| `ion-neutral_tasks.cpp` | Task registration | Evolution tasks |

## Configuration Parameters

From `<ion-neutral>` block:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `drag_coeff` | Real | required | Ion-neutral drag coefficient |
| `ionization_coeff` | Real | 0.0 | Ionization rate |
| `recombination_coeff` | Real | 0.0 | Recombination rate |

## Two-Fluid Equations

### Ion Fluid
$$\frac{\partial \rho_i}{\partial t} + \nabla \cdot (\rho_i \mathbf{v}_i) = S_i$$
$$\frac{\partial (\rho_i \mathbf{v}_i)}{\partial t} + \nabla \cdot (\rho_i \mathbf{v}_i \mathbf{v}_i + P_i) = \mathbf{J} \times \mathbf{B} + \mathbf{F}^{\text{drag}}$$
$$\frac{\partial E_i}{\partial t} + \nabla \cdot ((E_i + P_i)\mathbf{v}_i) = \mathbf{J} \cdot \mathbf{E} + Q^{\text{drag}}$$

### Neutral Fluid
$$\frac{\partial \rho_n}{\partial t} + \nabla \cdot (\rho_n \mathbf{v}_n) = S_n$$
$$\frac{\partial (\rho_n \mathbf{v}_n)}{\partial t} + \nabla \cdot (\rho_n \mathbf{v}_n \mathbf{v}_n + P_n) = -\mathbf{F}^{\text{drag}}$$
$$\frac{\partial E_n}{\partial t} + \nabla \cdot ((E_n + P_n)\mathbf{v}_n) = -Q^{\text{drag}}$$

## Drag Force

$$\mathbf{F}^{\text{drag}} = \alpha^{\text{drag}} \cdot \frac{\rho_i \rho_n}{\rho_i + \rho_n} \cdot (\mathbf{v}_i - \mathbf{v}_n)$$

## Source Terms

### Ionization
$$S_i = +\alpha_{\text{ion}} \cdot \rho_n$$
$$S_n = -\alpha_{\text{ion}} \cdot \rho_n$$

### Recombination
$$S_i = -\alpha_{\text{rec}} \cdot \rho_i^2$$
$$S_n = +\alpha_{\text{rec}} \cdot \rho_i^2$$

## See Also
- [MHD Module](mhd.md)
- Source: `src/ion-neutral/ion-neutral.cpp`