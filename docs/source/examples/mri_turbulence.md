# Example: MRI Turbulence

Magnetorotational instability (MRI) in a shearing box - the fundamental mechanism for angular momentum transport in accretion disks.

## Problem Generators
- **2D MRI**: `src/pgen/mri2d.cpp`
- **3D MRI**: `src/pgen/mri3d.cpp`

## Available Input Files
- **2D Setup**: `inputs/mhd/mri2d.athinput`
- **3D Setup**: `inputs/mhd/mri3d.athinput`
- **High Resolution**: `inputs/mhd/mri3d_hires.athinput`

## Physics

The MRI occurs when a weak magnetic field destabilizes a differentially rotating disk:
- Linear growth rate: $\gamma = 0.75 \Omega$ for optimal wavelength
- Critical wavelength: $\lambda_{\text{MRI}} = 2\pi v_A/\Omega$
- Requires: $\lambda_{\text{MRI}} > 6-8 \Delta x$ for resolution

## Running the Simulation

```bash
# Build with MRI problem
cmake -B build -DPROBLEM=mri3d
make -C build -j8

# Run 3D MRI
./build/src/athena -i inputs/mhd/mri3d.athinput
```

## Complete Input File

```ini
<shearing_box>
omega0 = 1.0          # Orbital frequency
qshear = 1.5          # Keplerian shear

<mhd>
eos = ideal
gamma = 1.66667       # 5/3
rsolver = hlld

<problem>
beta = 400.0          # Plasma beta
amp = 0.025           # Perturbation amplitude
ifield = 1            # 1=vertical, 2=toroidal

<meshblock>
nx1 = 32              # Radial
nx2 = 64              # Azimuthal  
nx3 = 32              # Vertical
```

## Analysis

Monitor volume-averaged stresses:

Reynolds stress:
$$\alpha_{\text{Re}} = -\frac{\langle \rho \delta v_x \delta v_y \rangle}{\langle P \rangle}$$

Maxwell stress:
$$\alpha_{\text{Max}} = \frac{\langle B_x B_y \rangle}{4\pi \langle P \rangle}$$

Total stress ($\alpha$-parameter):
$$\alpha = \alpha_{\text{Re}} + \alpha_{\text{Max}}$$

## See Also
- [Shearing Box Module](../modules/shearing_box.md)
- [MHD Module](../modules/mhd.md)