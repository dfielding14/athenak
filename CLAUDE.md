# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AthenaK is a high-performance astrophysical simulation framework for solving fluid dynamics, magnetohydrodynamics (MHD), and numerical relativity problems. It's a complete rewrite of Athena++ using Kokkos for performance portability across CPUs and GPUs.

## Build Commands

### Basic CPU Build
```bash
mkdir build
cd build
cmake ..
make -j8
```

### MPI-Enabled Build
```bash
cmake -B build -DAthena_ENABLE_MPI=ON
cd build && make -j8
```

### GPU Build (CUDA)
```bash
cmake -B build -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON \
      -DCMAKE_CXX_COMPILER=/path/to/kokkos/bin/nvcc_wrapper
cd build && make -j8
```

### Debug Build
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cd build && make -j8
```

### Custom Problem Generator
```bash
cmake -B build -DPROBLEM=your_problem_name
cd build && make -j8
```

## Test Commands

### Run All Tests
```bash
cd tst
python run_tests.py
```

### Run Specific Test Suite
```bash
cd tst
python run_tests.py hydro mhd  # Run only hydro and mhd tests
```

### Run Single Test
After building, from the build directory:
```bash
cd build/src
./athena -i ../../inputs/tests/linear_wave_hydro.athinput
```

## Code Style and Linting

### C++ Style Check
```bash
cd tst/scripts/style
bash check_athena_cpp_style.sh
```

Key style requirements:
- No tabs (spaces only)
- Line length limit: 90 characters
- Single closing brace per line
- No trailing whitespace

### Python Style Check
```bash
flake8 tst/ vis/
```

## Architecture and Key Components

### Core Source Structure (`src/`)
- **mesh/**: AMR mesh infrastructure with MeshBlock as the fundamental unit
- **hydro/**: Hydrodynamics solver for Euler equations
- **mhd/**: Magnetohydrodynamics solver
- **driver/**: Main simulation driver and time-stepping control
- **pgen/**: Problem generators (initial conditions) - custom problems go here
- **coordinates/**: Coordinate systems (Cartesian, spherical, cylindrical)
- **outputs/**: Output formats (VTK, HDF5, binary) 
- **bvals/**: Boundary values and MPI communication
- **reconstruct/**: Spatial reconstruction methods (PLM, PPM, WENOZ)

### Advanced Physics Modules
- **dyn_grmhd/**: General relativistic MHD in dynamical spacetimes
- **radiation/**: Radiation transport solver
- **z4c/**: Numerical relativity using Z4c formalism
- **particles/**: Lagrangian particle tracking
- **srcterms/**: Source terms including turbulence driving (currently working on SFB turbulence driver)

### Key Design Patterns
1. **Task-Based Execution**: Uses TaskList for managing computational tasks
2. **Kokkos Views**: All array data uses Kokkos::View for performance portability
3. **MeshBlock Structure**: Domain decomposed into MeshBlocks for AMR and parallelization
4. **Physics Modules**: Each physics solver (hydro, MHD, etc.) is a separate module with standardized interfaces

### Problem Generator Development
When creating a new problem generator:
1. Create file in `src/pgen/` following existing patterns
2. Must define `void ProblemGenerator(ParameterInput *pin, const bool restart)`
3. Access mesh data through `pm->pmb_pack` 
4. Use Kokkos parallel constructs for initialization
5. Build with `-DPROBLEM=your_problem_name`

### Current Development Focus
