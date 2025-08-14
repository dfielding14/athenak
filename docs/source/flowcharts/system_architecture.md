# AthenaK System Architecture

## Complete System Overview

```{mermaid}
flowchart TB
    subgraph Input["Input Layer"]
        ATHINPUT[athinput file]
        RESTART[restart file]
        EXTDATA[external data]
    end

    subgraph Core["Core Infrastructure"]
        MAIN[main.cpp]
        MESH[Mesh Module]
        DRIVER[Driver Module]
        TASKS[TaskList]
        COORDS[Coordinates]
    end

    subgraph Physics["Physics Modules"]
        HYDRO[Hydrodynamics]
        MHD[Magnetohydrodynamics]
        RAD[Radiation]
        Z4C[Z4c Relativity]
        GRMHD[DynGRMHD]
        PARTICLES[Particles]
        IONNEUTRAL[Ion-Neutral]
    end

    subgraph Numerical["Numerical Methods"]
        RECON[Reconstruction]
        RIEMANN[Riemann Solvers]
        EOS[Equations of State]
        DIFF[Diffusion]
    end

    subgraph Support["Support Systems"]
        BVALS[Boundaries/MPI]
        OUTPUTS[Outputs]
        SRCTERMS[Source Terms]
        SHEARBOX[Shearing Box]
        PGEN[Problem Generators]
    end

    subgraph Output["Output Layer"]
        VTK[VTK Files]
        BIN[Binary Files]
        RST[Restart Files]
        HST[History Files]
        PART[Particle Files]
    end

    %% Input connections
    ATHINPUT --> MAIN
    RESTART --> MAIN
    EXTDATA --> PGEN

    %% Core connections
    MAIN --> MESH
    MESH --> DRIVER
    DRIVER --> TASKS
    MESH --> COORDS
    
    %% Physics initialization
    DRIVER --> HYDRO
    DRIVER --> MHD
    DRIVER --> RAD
    DRIVER --> Z4C
    Z4C --> GRMHD
    DRIVER --> PARTICLES
    DRIVER --> IONNEUTRAL

    %% Numerical methods used by physics
    HYDRO --> RECON
    HYDRO --> RIEMANN
    HYDRO --> EOS
    MHD --> RECON
    MHD --> RIEMANN
    MHD --> EOS
    MHD --> DIFF
    RAD --> RECON
    GRMHD --> RIEMANN
    GRMHD --> EOS

    %% Support systems
    TASKS --> BVALS
    DRIVER --> OUTPUTS
    HYDRO --> SRCTERMS
    MHD --> SRCTERMS
    MHD --> SHEARBOX
    MESH --> PGEN

    %% Output generation
    OUTPUTS --> VTK
    OUTPUTS --> BIN
    OUTPUTS --> RST
    OUTPUTS --> HST
    PARTICLES --> PART

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
flowchart LR
    subgraph Init["Initialization Phase"]
        START[Start] --> PARSE[Parse Input]
        PARSE --> BUILD[Build Mesh]
        BUILD --> DECOMP[Domain Decomposition]
        DECOMP --> PGEN_INIT[Problem Generator]
    end

    subgraph Evolution["Evolution Loop"]
        PGEN_INIT --> LOOP{t < tlim?}
        LOOP -->|Yes| STAGE[RK Stage]
        STAGE --> TASKS_EXEC[Execute Tasks]
        TASKS_EXEC --> UPDATE[Update Variables]
        UPDATE --> NEWDT[Calculate dt]
        NEWDT --> OUTPUT_CHECK{Output?}
        OUTPUT_CHECK -->|Yes| WRITE[Write Files]
        OUTPUT_CHECK -->|No| AMR_CHECK{Refine?}
        WRITE --> AMR_CHECK
        AMR_CHECK -->|Yes| REFINE[Mesh Refinement]
        AMR_CHECK -->|No| LOOP
        REFINE --> LOOP
    end

    subgraph Finalize["Finalization"]
        LOOP -->|No| FINAL_OUT[Final Output]
        FINAL_OUT --> CLEANUP[Cleanup]
        CLEANUP --> END[End]
    end

    style Init fill:#e3f2fd
    style Evolution fill:#fff9c4
    style Finalize fill:#efebe9
    
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