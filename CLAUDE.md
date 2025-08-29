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

### Verify Installation
```bash
# Quick verification test (recommended after build)
cd tst
python run_tests.py hydro/hydro_linwave
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
- Left-justified pragmas
- No executable permissions on source files

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

Example template in: `src/pgen/tests/linear_wave_hydro.cpp`

## Input File Format

AthenaK uses `.athinput` files with hierarchical parameter blocks:

```
<mesh>
nghost = 2        # Number of ghost cells
nx1 = 64          # Zones in X1-direction

<time>
cfl_number = 0.3  # CFL stability parameter
tlim = 5.0        # Time limit

<hydro>
reconstruct = plm # Spatial reconstruction (plm, ppm, wenoz)
rsolver = llf     # Riemann solver (llf, hllc, hlld, etc.)
```

Reference input files in `inputs/tests/` directory.

## Useful Development Commands

### Analyze Code Without Running
```bash
./athena -c   # Display configuration
./athena -m   # Output mesh structure  
./athena -n   # Parse input file only
```

### Visualization Tools
```bash
# Python analysis scripts in vis/python/
python vis/python/plot_slice.py output.vtk
```

## CMake Build Options

Key configuration flags:
- `Athena_SINGLE_PRECISION`: Use single precision (default: OFF)
- `Athena_ENABLE_MPI`: Enable MPI parallelization (default: OFF)
- `Athena_ENABLE_OPENMP`: Enable OpenMP (default: OFF)
- `PROBLEM`: Custom problem generator name (default: built_in_pgens)
- `CMAKE_BUILD_TYPE`: Debug/Release/RelWithDebInfo (default: Release)

## Dependencies

- **Required**: C++17 compiler, CMake 3.13+
- **Kokkos**: Automatically fetched and built (performance portability layer)
- **Optional**: MPI, OpenMP, GSL (for specific problems)

## Athena++ Reference Directory

The `athena++/` directory contains the original Athena++ codebase for reference. This is the pre-Kokkos version that AthenaK was rewritten from. Use it for:

### Migration and Porting
- **Problem generators**: Compare `athena++/src/pgen/` with `src/pgen/` to understand porting patterns
- **Algorithm reference**: Original implementations without Kokkos abstractions
- **Input file compatibility**: Many `.athinput` files work with minor modifications

### Key Migration Patterns
When porting from Athena++ to AthenaK:
1. **Arrays**: `AthenaArray<Real>` → `Kokkos::View` with appropriate memory space
2. **Loops**: Nested for loops → `par_for` with Kokkos lambdas
3. **MeshBlock access**: Direct member access → `pmb_pack` interface
4. **Physics modules**: Similar structure but with Kokkos execution policies
5. **Boundary conditions**: Function pointers → Task-based BC application

### Common Gotchas
- **Memory spaces**: Host vs Device memory - all computation arrays must be device-accessible
- **Lambda captures**: Use `KOKKOS_LAMBDA` for device-compatible lambdas
- **Reductions**: Replace manual reductions with Kokkos reduction patterns
- **Random access**: Array indexing patterns must be GPU-friendly (coalesced access)

### Useful Comparisons
```bash
# Compare problem generator implementations
diff athena++/src/pgen/blast.cpp src/pgen/blast.cpp

# Check for equivalent physics modules
ls athena++/src/hydro/ src/hydro/

# Compare input file parameters
diff athena++/inputs/hydro/athinput.blast inputs/tests/blast_hydro.athinput
```

## Current Development Focus

The current branch (`dev/particles_turb_divb`) implements:
- Spherical Fourier-Bessel turbulence driving in `src/srcterms/turb_driver.hpp/cpp`
- Particle tracking enhancements
- Divergence-free magnetic field verification tests

## References

- [AthenaK Wiki](https://github.com/IAS-Astrophysics/athenak/wiki)
- [Athena++ Documentation](https://github.com/PrincetonUniversity/athena/wiki) (for similar concepts)
- Code papers: Stone et al. (2024), Zhu et al. (2024), Fields et al. (2024)
- Athena++ paper: Stone et al. (2020, ApJS 249, 4) - Original framework design

- memorize Do not build or test the code unless I explicitly tell you to