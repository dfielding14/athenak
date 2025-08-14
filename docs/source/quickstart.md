# 5-Minute Quickstart

Get AthenaK running in 5 minutes with this step-by-step guide.

## Step 1: Clone and Build (2 minutes)

```bash
# Clone the repository
git clone https://github.com/IAS-Astrophysics/athenak.git
cd athenak

# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make -j8
```

## Step 2: Run Your First Simulation (1 minute)

```bash
# Navigate to executable
cd src

# Run classic Sod shock tube test
./athena -i ../../inputs/hydro/sod.athinput
```

## Step 3: View the Output (1 minute)

```bash
# List output files
ls *.vtk

# Output will be:
# shock.block0.out1.00000.vtk
# shock.block0.out1.00001.vtk
# ...
```

## Step 4: Visualize (1 minute)

Open the VTK files in ParaView or VisIt:
1. Open ParaView
2. File → Open → Select all `*.vtk` files
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
mpirun -np 4 ./athena -i input.athinput
```

### Run on GPU
```bash
# Build with CUDA support
cmake -DKokkos_ENABLE_CUDA=ON ..
make -j8

# Run on GPU (automatic detection)
./athena -i input.athinput
```

### Restart from Checkpoint
```bash
# Create restart file
./athena -i input.athinput

# Restart from last checkpoint
./athena -r shock.00010.rst
```

### Change Problem Type
```bash
# Build with different problem generator
cmake -DPROBLEM=blast ..
make -j8
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