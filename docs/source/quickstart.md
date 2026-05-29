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
# Classic Sod shock tube test, with outputs kept in one run directory
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod
```

## Step 3: Inspect the Output (1 minute)

```bash
# List output files
ls run-sod/tab/Sod.hydro_w.*.tab run-sod/Sod.hydro.hst

# Typical files:
# run-sod/tab/Sod.hydro_w.00000.tab  (primitive-variable slice)
# run-sod/Sod.hydro.hst              (history diagnostics)

# Peek at the tabular output
head run-sod/tab/Sod.hydro_w.00000.tab
```

## Step 4: Visualize (1 minute)

The default Sod problem writes tabular and history data. To view spatial fields:

```bash
# Option A: plot the tab output with your favourite tool (Python/gnuplot)
python - <<'PY'
import numpy as np
data = np.loadtxt("run-sod/tab/Sod.hydro_w.00000.tab")
print("Columns:", data.shape[1])
print("First rows:", data[:5])
PY

# Option B: run one short cycle of a built-in MHD deck that emits VTK files
./build/src/athena -i inputs/mhd/orszag_tang.athinput -d run-ot time/nlim=1
ls run-ot/vtk/OrszagTang.mhd_w.*.vtk run-ot/vtk/OrszagTang.mhd_bcc.*.vtk
```

Then open the generated VTK files in ParaView or VisIt:
1. Open ParaView
2. File → Open → Select all `*.vtk` files (either from the optional run above or your own deck)
3. Click Apply
4. Choose a variable to visualize

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
Override the shipped deck without modifying it:
```bash
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod-512 mesh/nx1=512
```

### Change Output Frequency
```bash
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod-frequent output1/dt=0.005
```

### Try Different Physics

#### MHD Shock Tube (Brio-Wu)
```bash
./build/src/athena -i inputs/mhd/bw.athinput -d run-bw
```

#### Orszag-Tang MHD Vortex
```bash
./build/src/athena -i inputs/mhd/orszag_tang.athinput -d run-ot-full
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
The Sod deck writes tables and history data, not restart checkpoints. Add a
restart stream to a run-specific input deck:

```bash
cp inputs/hydro/sod.athinput my_sod_restart.athinput
```

```ini
<output3>
file_type = rst
dt        = 0.1
```

Then run and resume from a produced checkpoint:

```bash
./build/src/athena -i my_sod_restart.athinput -d run-sod-rst
./build/src/athena -r run-sod-rst/rst/Sod.00001.rst -d run-sod-resumed
```

The `00000` dump is written during initialization; select a later dump to
demonstrate continuation from evolved state.

### Change Problem Type
```bash
# Run a different built-in problem without rebuilding
./build/src/athena -i inputs/hydro/lw_implode.athinput -d run-implode

# A new source file src/pgen/my_problem.cpp is selected at build time:
cmake -S . -B build-my-problem -DPROBLEM=my_problem
cmake --build build-my-problem -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Some shipped decks correspond to a source file selected at build time rather
than a `pgen_name` handled by the default executable. For example, run the
field-loop or blast deck only after configuring `-DPROBLEM=field_loop` or
`-DPROBLEM=blast`, respectively.

## Quick Input File Reference

Start from a shipped deck such as `inputs/hydro/sod.athinput`: it includes the
required mesh extents, boundary conditions, `<meshblock>`, `<time>`, physics,
problem-generator, and output blocks. See the
[Configuration Guide](configuration.md) before authoring a deck from scratch.

## Troubleshooting

### Build Errors
```bash
# Missing Kokkos
git submodule update --init --recursive

# Inspect compiled options or parse an input without running
./build/src/athena -c
./build/src/athena -i inputs/hydro/sod.athinput -n
```

### Runtime Errors
```bash
# For a longer-running deck, set an AthenaK wall-clock limit below the queue allocation
./build/src/athena -i long_running.athinput -t 00:05:00

# Load an installed Kokkos Tools profiling library
export KOKKOS_PROFILE_LIBRARY=/path/to/libkokkos-tools.so
./build/src/athena -i inputs/hydro/sod.athinput
```

## What's Next?

### Learn More About:
- **[System Overview](overview.md)** - Understand the architecture
- **[Configuration Guide](configuration.md)** - Detailed input options
- **[Physics Modules](modules/index.md)** - Available physics
- **[Problem Generators](modules/pgen.md)** - Set up new problems

### Try These Examples:
1. `inputs/hydro/sod.athinput` - hydrodynamic shock tube
2. `inputs/mhd/bw.athinput` - MHD shock tube
3. `inputs/mhd/orszag_tang.athinput` - MHD vortex with VTK output
4. `inputs/hydro/lw_implode.athinput` - 2D hydrodynamic implosion with VTK output

## Quick Command Reference

| Command | Purpose |
|---------|---------|
| `./build/src/athena -i file.athinput` | Run a new simulation |
| `./build/src/athena -r rst/file.rst` | Restart from a checkpoint |
| `./build/src/athena -i file.athinput -n` | Parse input and exit |
| `./build/src/athena -c` | Show compiled configuration and exit |
| `./build/src/athena -m -i file.athinput` | Write mesh structure and exit |
| `./build/src/athena -i file.athinput -t hh:mm:ss` | Set wall-clock limit for final output |
| `./build/src/athena -h` | Show help |

---

**Ready for more?** → [Full Configuration Guide](configuration.md)
