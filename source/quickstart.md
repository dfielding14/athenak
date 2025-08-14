# Quick Start Guide

Get AthenaK running in 5 minutes!

## Prerequisites

- C++17 compiler (GCC 8+, Clang 10+, or NVCC 11+)
- CMake 3.13+
- MPI (optional, for parallel runs)
- Git

## 1. Clone and Build

```bash
# Clone the repository
git clone https://github.com/IAS-Astrophysics/athenak.git
cd athenak

# Create build directory
mkdir build && cd build

# Configure (basic CPU build)
cmake ..

# Build (use -j8 for parallel compilation)
make -j8
```

## 2. Run Your First Simulation

```bash
# From the build directory
cd src

# Run a shock tube test
./athena -i ../../inputs/tests/shock_tube.athinput

# Output will be in the current directory
ls *.vtk
```

## 3. Visualize Results

Open the VTK files in ParaView or VisIt:
```bash
# If ParaView is installed
paraview ShockTube.*.vtk
```

## Build Options

### MPI Parallel Build
```bash
cmake -DAthena_ENABLE_MPI=ON ..
make -j8

# Run with MPI
mpirun -np 4 ./athena -i input.athinput
```

### GPU Build (NVIDIA)
```bash
cmake -DKokkos_ENABLE_CUDA=ON \
      -DKokkos_ARCH_VOLTA70=ON \
      -DCMAKE_CXX_COMPILER=$KOKKOS_PATH/bin/nvcc_wrapper ..
make -j8
```

### Debug Build
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j8
```

## Input File Basics

AthenaK uses `.athinput` files to configure simulations:

```ini
# Example: 1D Shock Tube
<job>
basename = shock_tube

<mesh>
nx1 = 256        # Resolution
x1min = 0.0      # Domain
x1max = 1.0
ix1_bc = outflow # Boundary conditions
ox1_bc = outflow

<time>
cfl_number = 0.4
tlim = 0.25

<hydro>
gamma = 1.4
reconstruct = plm
rsolver = hllc

<problem>
pgen_name = shock_tube
```

## Common Problem Generators

| Problem | Description | Input File |
|---------|-------------|------------|
| `shock_tube` | 1D Riemann problems | `inputs/tests/shock_tube.athinput` |
| `blast` | Spherical blast wave | `inputs/tests/blast.athinput` |
| `linear_wave` | Linear wave convergence | `inputs/tests/linear_wave_hydro.athinput` |
| `kh` | Kelvin-Helmholtz instability | `inputs/hydro/kh.athinput` |
| `turb` | Driven turbulence | `inputs/hydro/turb.athinput` |

## Next Steps

### Learn More
- [Building Guide](building.md) - Detailed build options
- [Configuration](configuration.md) - Input file parameters
- [Module Guides](modules/index.md) - Physics modules

### Try Examples
```bash
# MHD linear waves
./athena -i ../../inputs/tests/linear_wave_mhd.athinput

# 3D blast wave
./athena -i ../../inputs/tests/blast_3d.athinput

# Turbulence
./athena -i ../../inputs/hydro/turb.athinput
```

### Customize
1. Copy an existing problem generator from `src/pgen/`
2. Modify for your problem
3. Build with `-DPROBLEM=your_problem`

## Troubleshooting

### Build Errors

**Issue**: CMake can't find MPI
```bash
# Specify MPI compiler
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
```

**Issue**: CUDA not detected
```bash
# Set CUDA path
export CUDA_HOME=/usr/local/cuda
cmake -DKokkos_ENABLE_CUDA=ON ..
```

### Runtime Errors

**Issue**: Timestep too small
- Reduce CFL number in input file
- Check for unresolved features

**Issue**: NaN in solution
- Check initial conditions
- Ensure proper EOS parameters
- Try first-order reconstruction

## Getting Help

- [Documentation](index.rst)
- [GitHub Issues](https://github.com/IAS-Astrophysics/athenak/issues)
- [Wiki](https://github.com/IAS-Astrophysics/athenak/wiki)

## Quick Reference Card

### Essential Commands
```bash
# Build
cmake .. && make -j

# Run
./athena -i input.athinput

# Check configuration
./athena -c

# Parse input without running
./athena -n -i input.athinput

# Show mesh structure
./athena -m -i input.athinput
```

### Key Input Parameters
```ini
<mesh>
nx1, nx2, nx3     # Resolution
x1min, x1max      # Domain size
ix1_bc, ox1_bc    # Boundaries

<time>
cfl_number        # Stability (0.1-0.4)
tlim              # End time
dt                # Fixed timestep

<hydro>/<mhd>
gamma             # Adiabatic index
reconstruct       # plm, ppm, wenoz
rsolver           # hllc, hlld, roe
```