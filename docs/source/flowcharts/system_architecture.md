# AthenaK System Architecture

## Complete System Overview

```{mermaid}
flowchart TD
    subgraph Input["📥 Input Layer"]
        direction LR
        ATHINPUT[athinput file]
        RESTART[restart file]
        EXTDATA[external data]
    end

    subgraph Core["⚙️ Core Infrastructure"]
        MAIN[main.cpp<br/>Entry Point]
        MESH[Mesh Module<br/>Domain Decomposition]
        DRIVER[Driver Module<br/>Time Integration]
        TASKS[TaskList<br/>Execution Manager]
        COORDS[Coordinates<br/>Geometry]
    end

    subgraph PhysicsRow1["🔬 Physics Modules (Primary)"]
        direction LR
        HYDRO[Hydrodynamics<br/>Euler Equations]
        MHD[MHD<br/>Maxwell+Fluid]
        RAD[Radiation<br/>Transport]
    end

    subgraph PhysicsRow2["🔬 Physics Modules (Advanced)"]
        direction LR
        Z4C[Z4c<br/>Relativity]
        GRMHD[DynGRMHD<br/>Relativistic MHD]
        PARTICLES[Particles<br/>Lagrangian]
        IONNEUTRAL[Ion-Neutral<br/>Two-Fluid]
    end

    subgraph Numerical["🔢 Numerical Methods"]
        direction LR
        RECON[Reconstruction<br/>PLM/PPM/WENOZ]
        RIEMANN[Riemann<br/>Solvers]
        EOS[EOS<br/>Thermodynamics]
        DIFF[Diffusion<br/>Viscosity/Resistivity]
    end

    subgraph Support["🛠️ Support Systems"]
        direction LR
        BVALS[Boundaries<br/>MPI Comm]
        OUTPUTS[Outputs<br/>I/O Manager]
        SRCTERMS[Source Terms<br/>Forces/Cooling]
        SHEARBOX[Shearing Box<br/>Disks]
        PGEN[Problem<br/>Generators]
    end

    subgraph Output["💾 Output Layer"]
        direction LR
        VTK[VTK<br/>Visualization]
        BIN[Binary<br/>Analysis]
        RST[Restart<br/>Checkpoints]
        HST[History<br/>Time Series]
        PART[Particle<br/>Tracking]
    end

    %% Vertical flow - main pipeline
    Input --> MAIN
    MAIN --> MESH
    MESH --> DRIVER
    DRIVER --> TASKS
    TASKS --> PhysicsRow1
    TASKS --> PhysicsRow2
    PhysicsRow1 --> Numerical
    PhysicsRow2 --> Numerical
    Numerical --> Support
    Support --> OUTPUTS
    OUTPUTS --> Output

    %% Core connections
    MESH -.-> COORDS
    COORDS -.-> PhysicsRow1
    COORDS -.-> PhysicsRow2
    
    %% Physics details
    Z4C --> GRMHD
    
    %% Support connections
    TASKS -.-> BVALS
    PhysicsRow1 -.-> SRCTERMS
    MHD -.-> SHEARBOX
    MESH -.-> PGEN
    EXTDATA -.-> PGEN
    PARTICLES -.-> PART

    %% Click handlers for navigation
    click MESH "../modules/mesh.html" "Go to Mesh Module"
    click DRIVER "../modules/driver.html" "Go to Driver Module"
    click TASKS "../modules/tasklist.html" "Go to TaskList Module"
    click COORDS "../modules/coordinates.html" "Go to Coordinates Module"
    click HYDRO "../modules/hydro.html" "Go to Hydro Module"
    click MHD "../modules/mhd.html" "Go to MHD Module"
    click RAD "../modules/radiation.html" "Go to Radiation Module"
    click Z4C "../modules/z4c.html" "Go to Z4c Module"
    click GRMHD "../modules/dyn_grmhd.html" "Go to DynGRMHD Module"
    click PARTICLES "../modules/particles.html" "Go to Particles Module"
    click IONNEUTRAL "../modules/ion_neutral.html" "Go to Ion-Neutral Module"
    click RECON "../modules/reconstruction.html" "Go to Reconstruction Module"
    click RIEMANN "../modules/riemann_solvers.html" "Go to Riemann Solvers"
    click EOS "../modules/eos.html" "Go to EOS Module"
    click DIFF "../modules/diffusion.html" "Go to Diffusion Module"
    click BVALS "../modules/boundaries.html" "Go to Boundaries Module"
    click OUTPUTS "../modules/outputs.html" "Go to Outputs Module"
    click SRCTERMS "../modules/srcterms.html" "Go to Source Terms"
    click SHEARBOX "../modules/shearing_box.html" "Go to Shearing Box"
    click PGEN "../modules/pgen.html" "Go to Problem Generators"

    %% Styling
    classDef input fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef core fill:#fff3e0,stroke:#e65100,stroke-width:3px
    classDef physics fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef numerical fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px
    classDef support fill:#fce4ec,stroke:#880e4f,stroke-width:2px
    classDef output fill:#f3e5f5,stroke:#4a148c,stroke-width:2px

    class ATHINPUT,RESTART,EXTDATA input
    class MAIN,MESH,DRIVER,TASKS,COORDS core
    class HYDRO,MHD,RAD,Z4C,GRMHD,PARTICLES,IONNEUTRAL physics
    class RECON,RIEMANN,EOS,DIFF numerical
    class BVALS,OUTPUTS,SRCTERMS,SHEARBOX,PGEN support
    class VTK,BIN,RST,HST,PART output
```

## Execution Flow

```{mermaid}
flowchart TD
    subgraph Init["🚀 Initialization Phase"]
        START([Start Program])
        PARSE[Parse Input<br/>Read .athinput]
        BUILD[Build Mesh<br/>Create Grid]
        DECOMP[Domain Decomposition<br/>MPI Distribution]
        PGEN_INIT[Problem Generator<br/>Set Initial Conditions]
        
        START --> PARSE
        PARSE --> BUILD
        BUILD --> DECOMP
        DECOMP --> PGEN_INIT
    end

    subgraph Evolution["🔄 Main Evolution Loop"]
        LOOP{{"Time < tlim?"}}
        
        subgraph TimeStep["Time Step"]
            STAGE[RK Stage<br/>Integration]
            TASKS_EXEC[Execute Tasks<br/>Physics Modules]
            UPDATE[Update Variables<br/>Conservative → Primitive]
            NEWDT[Calculate dt<br/>CFL Condition]
        end
        
        subgraph Checks["Periodic Checks"]
            OUTPUT_CHECK{{"Output Time?"}}
            WRITE[Write Files<br/>VTK/Binary/History]
            AMR_CHECK{{"Refine Mesh?"}}
            REFINE[Mesh Refinement<br/>AMR Operations]
        end
        
        PGEN_INIT --> LOOP
        LOOP -->|Yes| STAGE
        STAGE --> TASKS_EXEC
        TASKS_EXEC --> UPDATE
        UPDATE --> NEWDT
        NEWDT --> OUTPUT_CHECK
        OUTPUT_CHECK -->|Yes| WRITE
        OUTPUT_CHECK -->|No| AMR_CHECK
        WRITE --> AMR_CHECK
        AMR_CHECK -->|Yes| REFINE
        AMR_CHECK -->|No| LOOP
        REFINE --> LOOP
    end

    subgraph Finalize["🏁 Finalization"]
        FINAL_OUT[Final Output<br/>Save Results]
        CLEANUP[Cleanup<br/>Free Memory]
        END([End Program])
        
        LOOP -->|No| FINAL_OUT
        FINAL_OUT --> CLEANUP
        CLEANUP --> END
    end

    style Init fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
    style Evolution fill:#fff9c4,stroke:#f57c00,stroke-width:2px
    style TimeStep fill:#ffffff,stroke:#666,stroke-width:1px
    style Checks fill:#ffffff,stroke:#666,stroke-width:1px
    style Finalize fill:#efebe9,stroke:#5d4037,stroke-width:2px
    
    click BUILD "../modules/mesh.html" "Go to Mesh Module"
    click PGEN_INIT "../modules/pgen.html" "Go to Problem Generators"
    click TASKS_EXEC "../modules/tasklist.html" "Go to TaskList Module"
    click WRITE "../modules/outputs.html" "Go to Outputs Module"
    click REFINE "../modules/mesh.html" "Go to Mesh Module"
```

## Task Execution Detail

```{mermaid}
flowchart TD
    subgraph TaskStage["Single RK Stage"]
        START_STAGE[Stage Start] --> RECV[Start MPI Recv]
        
        subgraph PhysicsTasks["Physics Tasks (Parallel)"]
            RECV --> HYDRO_TASK[Hydro Tasks]
            RECV --> MHD_TASK[MHD Tasks]
            RECV --> RAD_TASK[Radiation Tasks]
            
            HYDRO_TASK --> FLUX_H[Calculate Fluxes]
            MHD_TASK --> FLUX_M[Calculate Fluxes]
            MHD_TASK --> CT[Constrained Transport]
            RAD_TASK --> FLUX_R[Calculate Fluxes]
        end
        
        FLUX_H --> WAIT[Wait for MPI]
        FLUX_M --> WAIT
        FLUX_R --> WAIT
        CT --> WAIT
        
        WAIT --> UPDATE_ALL[Update All Variables]
        UPDATE_ALL --> BC[Apply BCs]
        BC --> STAGE_END[Stage Complete]
    end
    
    style PhysicsTasks fill:#c8e6c9
    
    click HYDRO_TASK "../modules/hydro.html" "Go to Hydro Module"
    click MHD_TASK "../modules/mhd.html" "Go to MHD Module"
    click RAD_TASK "../modules/radiation.html" "Go to Radiation Module"
    click CT "../modules/mhd.html" "Go to MHD Module"
    click BC "../modules/boundaries.html" "Go to Boundaries Module"
```

## Data Flow

```{mermaid}
flowchart LR
    subgraph MeshBlock["MeshBlock Data"]
        CONS[Conservative Variables]
        PRIM[Primitive Variables]
        BFIELD[Magnetic Field]
        METRIC[Metric Data]
    end

    subgraph Computation["Computation"]
        C2P[Cons to Prim]
        RECONSTRUCT[Reconstruction]
        RIEMANN_SOLVE[Riemann Solver]
        FLUXES[Fluxes]
    end

    subgraph Communication["MPI Communication"]
        PACK[Pack Buffers]
        SEND[MPI Send/Recv]
        UNPACK[Unpack Buffers]
    end

    CONS --> C2P
    C2P --> PRIM
    PRIM --> RECONSTRUCT
    RECONSTRUCT --> RIEMANN_SOLVE
    RIEMANN_SOLVE --> FLUXES
    FLUXES --> CONS

    PRIM --> PACK
    PACK --> SEND
    SEND --> UNPACK
    UNPACK --> PRIM

    style MeshBlock fill:#ffecb3
    style Computation fill:#b2dfdb
    style Communication fill:#d1c4e9
    
    click RECONSTRUCT "../modules/reconstruction.html" "Go to Reconstruction Module"
    click RIEMANN_SOLVE "../modules/riemann_solvers.html" "Go to Riemann Solvers"
    click PACK "../modules/boundaries.html" "Go to Boundaries Module"
    click SEND "../modules/boundaries.html" "Go to Boundaries Module"
```

## Memory Layout

```{mermaid}
flowchart TD
    subgraph Global["Global Memory"]
        MESH_TREE[Mesh Tree]
        TASK_LIST[Task List]
        PARAMS[Parameters]
    end

    subgraph PerRank["Per MPI Rank"]
        MB_PACK[MeshBlockPack]
        
        subgraph PerBlock["Per MeshBlock"]
            COORDS_DATA[Coordinates]
            CONS_VARS[Conservative Vars]
            PRIM_VARS[Primitive Vars]
            FLUX_ARRAYS[Flux Arrays]
            SCRATCH[Scratch Space]
        end
        
        MPI_BUFFERS[MPI Buffers]
        OUTPUT_BUFFERS[Output Buffers]
    end

    subgraph Device["Device Memory (GPU)"]
        KOKKOS_VIEWS[Kokkos Views]
        KERNELS[Compute Kernels]
    end

    MESH_TREE --> MB_PACK
    MB_PACK --> PerBlock
    PerBlock --> KOKKOS_VIEWS
    KOKKOS_VIEWS --> KERNELS

    style Global fill:#e1bee7
    style PerRank fill:#c5cae9
    style Device fill:#b3e5fc
    
    click MESH_TREE "../modules/mesh.html" "Go to Mesh Module"
    click TASK_LIST "../modules/tasklist.html" "Go to TaskList Module"
    click COORDS_DATA "../modules/coordinates.html" "Go to Coordinates Module"
```

## Module Dependencies

```{mermaid}
graph TD
    MAIN[main.cpp] --> MESH[Mesh]
    MESH --> COORDS[Coordinates]
    MESH --> BVALS[Boundaries]
    MESH --> DRIVER[Driver]
    
    DRIVER --> TASKS[TaskList]
    DRIVER --> OUTPUTS[Outputs]
    
    TASKS --> HYDRO[Hydro]
    TASKS --> MHD[MHD]
    TASKS --> RAD[Radiation]
    TASKS --> Z4C[Z4c]
    TASKS --> PARTICLES[Particles]
    
    HYDRO --> EOS[EOS]
    HYDRO --> RECON[Reconstruction]
    HYDRO --> RIEMANN[Riemann]
    
    MHD --> EOS
    MHD --> RECON
    MHD --> RIEMANN
    MHD --> DIFF[Diffusion]
    
    Z4C --> GRMHD[DynGRMHD]
    GRMHD --> EOS
    
    HYDRO --> SRCTERMS[Source Terms]
    MHD --> SRCTERMS
    MHD --> SHEARBOX[Shearing Box]
    
    MESH --> PGEN[Problem Generator]
    PGEN --> HYDRO
    PGEN --> MHD

    style MAIN fill:#ffcdd2,stroke:#b71c1c,stroke-width:3px
    style DRIVER fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style MESH fill:#fff3e0,stroke:#e65100,stroke-width:2px
    
    click MESH "../modules/mesh.html" "Go to Mesh Module"
    click COORDS "../modules/coordinates.html" "Go to Coordinates Module"
    click BVALS "../modules/boundaries.html" "Go to Boundaries Module"
    click DRIVER "../modules/driver.html" "Go to Driver Module"
    click TASKS "../modules/tasklist.html" "Go to TaskList Module"
    click OUTPUTS "../modules/outputs.html" "Go to Outputs Module"
    click HYDRO "../modules/hydro.html" "Go to Hydro Module"
    click MHD "../modules/mhd.html" "Go to MHD Module"
    click RAD "../modules/radiation.html" "Go to Radiation Module"
    click Z4C "../modules/z4c.html" "Go to Z4c Module"
    click PARTICLES "../modules/particles.html" "Go to Particles Module"
    click EOS "../modules/eos.html" "Go to EOS Module"
    click RECON "../modules/reconstruction.html" "Go to Reconstruction Module"
    click RIEMANN "../modules/riemann_solvers.html" "Go to Riemann Solvers"
    click DIFF "../modules/diffusion.html" "Go to Diffusion Module"
    click GRMHD "../modules/dyn_grmhd.html" "Go to DynGRMHD Module"
    click SRCTERMS "../modules/srcterms.html" "Go to Source Terms"
    click SHEARBOX "../modules/shearing_box.html" "Go to Shearing Box"
    click PGEN "../modules/pgen.html" "Go to Problem Generators"
```

## Performance Scaling

```{mermaid}
flowchart LR
    subgraph Hardware["Hardware Detection"]
        CPU[CPU Cores]
        GPU[GPU Available?]
        MPI[MPI Ranks]
    end

    subgraph Execution["Execution Strategy"]
        SERIAL[Serial]
        OPENMP[OpenMP]
        CUDA[CUDA]
        MPI_EXEC[MPI Parallel]
    end

    subgraph Optimization["Optimizations"]
        VECTOR[Vectorization]
        CACHE[Cache Blocking]
        OVERLAP[Comm/Comp Overlap]
        LOAD_BAL[Load Balance]
    end

    CPU --> OPENMP
    GPU -->|Yes| CUDA
    GPU -->|No| OPENMP
    MPI --> MPI_EXEC

    OPENMP --> VECTOR
    CUDA --> CACHE
    MPI_EXEC --> OVERLAP
    MPI_EXEC --> LOAD_BAL

    style Hardware fill:#e8eaf6
    style Execution fill:#c5e1a5
    style Optimization fill:#ffe0b2
    
    click MPI_EXEC "../modules/boundaries.html" "Go to Boundaries/MPI Module"
    click LOAD_BAL "../modules/mesh.html" "Go to Mesh Module"
```

## See Also
- [Module Index](../modules/index.md)
- [Configuration Guide](../configuration.md)
- [Running Simulations](../running.md)