# AthenaK System Overview

## What is AthenaK?

AthenaK is a complete rewrite of Athena++ using Kokkos for performance portability. It solves the equations of **hydrodynamics**, **magnetohydrodynamics (MHD)**, and **general relativity** on CPUs and GPUs with a single codebase.

## How It Works - The Complete Flow

```{mermaid}
flowchart TD
    subgraph Input["1️⃣ INPUT"]
        CONFIG[Configuration File<br/>.athinput]
        PROB[Problem Setup<br/>Initial Conditions]
    end

    subgraph Core["2️⃣ CORE SYSTEM"]
        MESH[Mesh Creation<br/>Domain Decomposition]
        BLOCKS[MeshBlocks<br/>Local Compute Units]
        PACK[MeshBlockPacks<br/>GPU Optimization]
    end

    subgraph Physics["3️⃣ PHYSICS"]
        HYDRO[Hydrodynamics<br/>Fluid Evolution]
        MHD[MHD<br/>Magnetic Fields]
        GR[General Relativity<br/>Curved Spacetime]
    end

    subgraph Evolution["4️⃣ TIME EVOLUTION"]
        TASKS[Task System<br/>Dependency Management]
        RK[Runge-Kutta<br/>Time Integration]
        AMR[AMR<br/>Adaptive Refinement]
    end

    subgraph Output["5️⃣ OUTPUT"]
        VTK[Visualization<br/>ParaView/VisIt]
        DATA[Data Files<br/>Analysis]
        RESTART[Checkpoints<br/>Restart Files]
    end

    CONFIG --> MESH
    PROB --> MESH
    MESH --> BLOCKS
    BLOCKS --> PACK
    
    PACK --> HYDRO
    PACK --> MHD
    PACK --> GR
    
    HYDRO --> TASKS
    MHD --> TASKS
    GR --> TASKS
    
    TASKS --> RK
    RK --> AMR
    AMR --> VTK
    AMR --> DATA
    AMR --> RESTART

    style Input fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
    style Core fill:#fff3e0,stroke:#f57c00,stroke-width:2px
    style Physics fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
    style Evolution fill:#e8f5e9,stroke:#388e3c,stroke-width:2px
    style Output fill:#fce4ec,stroke:#c2185b,stroke-width:2px
    
    click CONFIG "configuration.html" "Go to Configuration Guide"
    click PROB "modules/pgen.html" "Go to Problem Generators"
    click MESH "modules/mesh.html" "Go to Mesh Module"
    click BLOCKS "modules/mesh.html" "Go to Mesh Module"
    click HYDRO "modules/hydro.html" "Go to Hydro Module"
    click MHD "modules/mhd.html" "Go to MHD Module"
    click GR "modules/z4c.html" "Go to Z4c Module"
    click TASKS "modules/tasklist.html" "Go to TaskList Module"
    click AMR "modules/mesh.html#adaptive-mesh-refinement" "Go to AMR Documentation"
    click VTK "modules/outputs.html" "Go to Outputs Module"
    click DATA "modules/outputs.html" "Go to Outputs Module"
    click RESTART "modules/outputs.html#restart-files" "Go to Restart Documentation"
```

## Quick Start Path

### 1. Choose Your Physics
```ini
# In your .athinput file
<hydro>           # Fluid dynamics
eos = ideal       # Equation of state
gamma = 1.4       # Adiabatic index

<mhd>             # Magnetohydrodynamics
eos = ideal       # Equation of state

<radiation>       # Radiation transport
arad = 1.0        # Radiation constant
```

### 2. Set Up Your Domain
```ini
<mesh>
nx1 = 256         # Resolution in x
x1min = -1.0      # Domain boundaries
x1max = 1.0
```

### 3. Run Your Simulation
```bash
./build/src/athena -i problem.athinput
```

## Core Concepts

### MeshBlocks - The Fundamental Unit
- Each MeshBlock is a small 3D grid (e.g., 16³ cells)
- Distributed across MPI ranks for parallel execution
- Grouped into MeshBlockPacks for GPU efficiency

### Task-Based Execution
- Physics operations are tasks with dependencies
- Automatically schedules and overlaps computation
- Enables efficient CPU/GPU execution

### Driver Orchestration
- `src/main.cpp` parses the CLI, loads the `.athinput`, and constructs the `Driver`
- The `Driver` builds global task lists (before-stage, stage, after-stage) and executes them every timestep
- Physics modules register their own tasks (e.g., hydro fluxes, radiation updates) so the driver can honour data dependencies

### Adaptive Mesh Refinement (AMR)
- Dynamically refines/coarsens grid
- Follows features of interest
- Maintains conservation

## Execution Pipeline (Code-Level View)

1. **Startup** (`src/main.cpp`)
   - Parses command-line arguments (`-i`, overrides, restart flags)
   - Loads the `.athinput` via `ParameterInput`
   - Constructs the global `Mesh` and its `MeshBlockPack`
2. **Driver Initialisation** (`Driver::Initialize`)
   - Each enabled physics module (`hydro`, `mhd`, `radiation`, etc.) attaches to the MeshBlockPack and registers tasks
   - Output streams are configured from `<output*>` blocks
3. **Time Loop** (`Driver::Execute`)
   - For every stage: execute task lists (communication, fluxes, source terms, AMR restriction/prolongation)
   - Update diagnostics, wall-clock checks, and AMR triggers
4. **Outputs / Restart** (`outputs/` modules)
   - Streams write according to their cadence (`dt`, `ncycle`, `ndump`)
   - Restart dumps (`file_type = rst`) capture the full MeshBlock state for `-r` restarts

## Where to Go Next

### 🎯 By Goal

| I want to... | Start here |
|-------------|------------|
| Run a simulation quickly | [Quickstart Guide](quickstart.md) |
| Understand the physics | [Physics Modules](modules/index.md) |
| Optimize performance | [Task System](modules/tasklist.md) |
| Set up a new problem | [Problem Generators](modules/pgen.md) |
| Analyze output | [Output Formats](modules/outputs.md) |

### 📚 By Physics

| Physics Type | Module | Key Features |
|-------------|---------|--------------|
| Fluid Dynamics | [Hydro](modules/hydro.md) | Shock capturing, multiple Riemann solvers |
| MHD | [MHD](modules/mhd.md) | Constrained transport, div(B)=0 |
| General Relativity | [Z4c](modules/z4c.md) | Black holes, neutron stars, gravitational waves |
| Radiation | [Radiation](modules/radiation.md) | M1 closure, implicit solver |
| Particles | [Particles](modules/particles.md) | Lagrangian tracking, cosmic rays |

### 🔧 By Task

| Task | Documentation |
|------|---------------|
| Configure simulation | [Input Parameters](reference/input_parameters.md) (340 parameters) |
| Build for GPU | [Building Guide](building.md) |
| Debug performance | [Task Execution](flowcharts/runtime.md) |
| Migrate from Athena++ | [Migration Guide](migration/index.md) |

## System Architecture Deep Dive

For complete understanding of how components interact:
→ **[Full System Architecture](flowcharts/system_architecture.md)**

## Key Design Principles

### 1. Performance Portability
- **Single source code** for CPUs and GPUs
- **Kokkos abstractions** hide hardware details
- **Automatic optimization** for target architecture

### 2. Modularity
- **Physics modules** are independent
- **Easy to extend** with new physics
- **Clean interfaces** between components

### 3. Scalability
- **MPI parallelization** across nodes
- **GPU acceleration** within nodes
- **AMR** for efficient resolution

### 4. Robustness
- **Conservative schemes** preserve physics
- **Error checking** throughout
- **Restart capability** for long runs

## Example: Shock Tube Simulation

Here's how the system processes a simple shock tube:

1. **Input**: Read initial discontinuity from `.athinput`
2. **Mesh**: Create uniform grid of MeshBlocks
3. **Physics**: Initialize left/right states in Hydro module
4. **Evolution**: 
   - Reconstruct variables at faces (PLM/PPM)
   - Solve Riemann problem (HLLC solver)
   - Update conserved variables
   - Apply boundary conditions
5. **Output**: Write VTK files for visualization
