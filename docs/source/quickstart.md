# 5-Minute Quickstart

Get AthenaK running in 5 minutes with this step-by-step guide.

## Step 1: Clone and Build (2 minutes)

```bash
# Clone the repository
git clone https://github.com/IAS-Astrophysics/athenak.git
cd athenak

# Pull bundled dependencies (e.g. Kokkos)
git submodule update --init --recursive

# Configure and build out-of-source
cmake -S . -B build
cmake --build build -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

## Step 2: Run Your First Simulation (1 minute)

```bash
# Run from the build tree
cd build/src

# Classic Sod shock tube test
./athena -i ../../inputs/hydro/sod.athinput
```

## Step 3: Inspect the Output (1 minute)

```bash
# List output files
ls Sod.*

# Typical files:
# Sod.block0.out1.00000.tab   (slice of primitive variables)
# Sod.block0.out2.00000.hst   (history diagnostics)

# Peek at the tabular output
head Sod.block0.out1.00000.tab
```

## Step 4: Visualize (1 minute)

The default Sod problem writes tabular and history data. To view spatial fields:

```bash
# Option A: plot the tab output with your favourite tool (Python/gnuplot)
python - <<'PY'
import numpy as np
data = np.loadtxt("Sod.block0.out1.00000.tab")
print("Columns:", data.shape[1])
print("Density head:", data[:5,0])
PY

# Option B: run a deck that emits VTK files
./athena -i ../../inputs/mhd/field_loop.athinput
ls Loop*.vtk
```

Then open the generated VTK files in ParaView or VisIt:
1. Open ParaView
2. File → Open → Select all `*.vtk` files (either from the optional run above or your own deck)
3. Click Apply
4. Choose variable to visualize (density, pressure, etc.)

## What Just Happened?

You simulated a **shock tube** - a fundamental test of fluid dynamics:

```{mermaid}
flowchart LR
    subgraph Initial["Initial Condition (t=0)"]
        LEFT[High Pressure<br/>ρ=1.0, P=1.0]
        RIGHT[Low Pressure<br/>ρ=0.125, P=0.1]
        LEFT ---|Discontinuity| RIGHT
    end
    
    subgraph Evolution["Evolution (t>0)"]
        SHOCK[Shock Wave →]
        CONTACT[Contact<br/>Discontinuity →]
        RAREFACTION[← Rarefaction<br/>Wave]
    end
    
    Initial --> Evolution
```

## Next: Customize Your Simulation

### Change Resolution
Edit the input file:
```ini
<mesh>
nx1 = 512  # Increase from 256 to 512
```

### Change Output Frequency
```ini
<output1>
dt = 0.01  # Output every 0.01 time units
```

### Try Different Physics

#### MHD Shock Tube (Brio-Wu)
```bash
./athena -i ../../inputs/mhd/bw.athinput
```

#### 2D MHD Blast Wave
```bash
./athena -i ../../inputs/mhd/blast_mhd.athinput
```

## Common Quick Tasks

### Run in Parallel (MPI)
```bash
# Configure a build with MPI support
cmake -S . -B build-mpi -DAthena_ENABLE_MPI=ON
cmake --build build-mpi -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

mpirun -np 4 ./build-mpi/src/athena -i inputs/hydro/sod.athinput
```

### Run on GPU
```bash
# Configure a CUDA build (set architecture to match your device)
cmake -S . -B build-cuda \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_AMPERE80=ON
cmake --build build-cuda -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

./build-cuda/src/athena -i inputs/hydro/sod.athinput
```

### Restart from Checkpoint
```bash
# Create restart file
./athena -i ../../inputs/hydro/sod.athinput

# Restart from last checkpoint (replace with the filename written in your run dir)
./athena -r Sod.00010.rst
```

### Change Problem Type
```bash
# Run a different built-in problem without rebuilding
./athena -i ../../inputs/hydro/blast.athinput

# Need a custom problem generator? Configure a fresh build directory:
cmake -S . -B build-blast -DPROBLEM=blast
cmake --build build-blast -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

## Quick Input File Reference

### Minimal Input File
```ini
<time>
tlim = 1.0        # End time
integrator = rk2  # Time integration

<mesh>
nx1 = 256         # Resolution
x1min = 0.0       # Domain
x1max = 1.0

<hydro>
eos = ideal       # Equation of state
gamma = 1.4       # Adiabatic index
rsolver = hllc    # Riemann solver

<output1>
file_type = vtk   # Output format
dt = 0.1          # Output interval
```

## Troubleshooting

### Build Errors
```bash
# Missing MPI
cmake -DAthena_ENABLE_MPI=OFF ..

# Missing Kokkos
git submodule update --init --recursive
```

### Runtime Errors
```bash
# Check input file syntax
./athena -i input.athinput -n

# Increase verbosity
./athena -i input.athinput -v 2
```

### Performance Issues
```bash
# Check task timing
./athena -i input.athinput -t

# Profile with Kokkos
export KOKKOS_PROFILE=1
./athena -i input.athinput
```

## What's Next?

### Learn More About:
- **[System Overview](overview.md)** - Understand the architecture
- **[Configuration Guide](configuration.md)** - Detailed input options
- **[Physics Modules](modules/index.md)** - Available physics
- **[Problem Generators](modules/pgen.md)** - Set up new problems

### Try These Examples:
1. **Kelvin-Helmholtz Instability** - Shear flow instability
2. **Orszag-Tang Vortex** - MHD turbulence
3. **Rayleigh-Taylor Instability** - Buoyancy-driven mixing
4. **MRI** - Magnetorotational instability

## Quick Command Reference

| Command | Purpose |
|---------|---------|
| `./athena -i file.athinput` | Run simulation |
| `./athena -r file.rst` | Restart from checkpoint |
| `./athena -n` | Parse input only (no run) |
| `./athena -t` | Show timing information |
| `./athena -v 2` | Verbose output |
| `./athena -h` | Show help |

---

**Ready for more?** → [Full Configuration Guide](configuration.md)
