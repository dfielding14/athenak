# AthenaK System Architecture

## Complete System Overview

```{mermaid}
flowchart TD
    subgraph INPUT[Input Layer]
        direction LR
        ATHINPUT[athinput file]
        RESTART[restart file]
        EXTDATA[external data]
    end

    subgraph CORE[Core Infrastructure]
        MAIN[main.cpp - Entry Point]
        MESH[Mesh Module - Domain Decomposition]
        DRIVER[Driver Module - Time Integration]
        TASKS[TaskList - Execution Manager]
        COORDS[Coordinates - Geometry]
    end

    subgraph PHYSICS1[Physics Modules - Primary]
        direction LR
        HYDRO[Hydrodynamics - Euler Equations]
        MHD[MHD - Maxwell+Fluid]
        RAD[Radiation - Transport]
    end

    subgraph PHYSICS2[Physics Modules - Advanced]
        direction LR
        Z4C[Z4c - Relativity]
        GRMHD[DynGRMHD - Relativistic MHD]
        PARTICLES[Particles - Lagrangian]
        IONNEUTRAL[Ion-Neutral - Two-Fluid]
    end

    subgraph NUMERICAL[Numerical Methods]
        direction LR
        RECON[Reconstruction - PLM/PPM/WENOZ]
        RIEMANN[Riemann - Solvers]
        EOS[EOS - Thermodynamics]
        DIFF[Diffusion - Viscosity/Resistivity]
    end

    subgraph SUPPORT[Support Systems]
        direction LR
        BVALS[Boundaries - MPI Comm]
        OUTPUTS[Outputs - I/O Manager]
        SRCTERMS[Source Terms - Forces/Cooling]
        SHEARBOX[Shearing Box - Disks]
        PGEN[Problem - Generators]
    end

    subgraph OUTPUT[Output Layer]
        direction LR
        VTK[VTK - Visualization]
        BIN[Binary - Analysis]
        RST[Restart - Checkpoints]
        HST[History - Time Series]
        PART[Particle - Tracking]
    end

    %% Vertical flow - main pipeline
    ATHINPUT --> MAIN
    RESTART --> MAIN
    MAIN --> MESH
    MESH --> DRIVER
    DRIVER --> TASKS
    TASKS --> HYDRO
    TASKS --> MHD
    TASKS --> Z4C
    HYDRO --> RECON
    MHD --> RIEMANN
    Z4C --> GRMHD
    RECON --> BVALS
    RIEMANN --> OUTPUTS
    OUTPUTS --> VTK
    OUTPUTS --> BIN

    %% Core connections
    MESH -.-> COORDS
    COORDS -.-> HYDRO
    COORDS -.-> MHD
    
    %% Physics details
    Z4C --> GRMHD
    
    %% Support connections
    TASKS -.-> BVALS
    MHD -.-> SRCTERMS
    MHD -.-> SHEARBOX
    MESH -.-> PGEN
    EXTDATA -.-> PGEN
    PARTICLES -.-> PART

    %% Click handlers for navigation
    click MESH "../modules/mesh.html"
    click DRIVER "../modules/driver.html"
    click TASKS "../modules/tasklist.html"
    click COORDS "../modules/coordinates.html"
    click HYDRO "../modules/hydro.html"
    click MHD "../modules/mhd.html"
    click RAD "../modules/radiation.html"
    click Z4C "../modules/z4c.html"
    click GRMHD "../modules/dyn_grmhd.html"
    click PARTICLES "../modules/particles.html"
    click RECON "../modules/reconstruction.html"
    click RIEMANN "../modules/riemann_solvers.html"
    click EOS "../modules/eos.html"
    click DIFF "../modules/diffusion.html"
    click BVALS "../modules/boundaries.html"
    click OUTPUTS "../modules/outputs.html"
    click SRCTERMS "../modules/srcterms.html"
    click SHEARBOX "../modules/shearing_box.html"
    click PGEN "../modules/pgen.html"

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
    subgraph INIT[Initialization Phase]
        START([Start Program])
        PARSE[Parse Input - Read .athinput]
        BUILD[Build Mesh - Create Grid]
        DECOMP[Domain Decomposition - MPI Distribution]
        PGEN_INIT[Problem Generator - Set Initial Conditions]
        
        START --> PARSE
        PARSE --> BUILD
        BUILD --> DECOMP
        DECOMP --> PGEN_INIT
    end

    subgraph EVOLUTION[Main Evolution Loop]
        LOOP{Time less than tlim?}
        
        subgraph TIMESTEP[Time Step]
            STAGE[RK Stage - Integration]
            TASKS_EXEC[Execute Tasks - Physics Modules]
            UPDATE[Update Variables - Conservative to Primitive]
            NEWDT[Calculate dt - CFL Condition]
        end
        
        subgraph CHECKS[Periodic Checks]
            OUTPUT_CHECK{Output Time?}
            WRITE[Write Files - VTK/Binary/History]
            AMR_CHECK{Refine Mesh?}
            REFINE[Mesh Refinement - AMR Operations]
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

    subgraph FINALIZE[Finalization]
        FINAL_OUT[Final Output - Save Results]
        CLEANUP[Cleanup - Free Memory]
        END([End Program])
        
        LOOP -->|No| FINAL_OUT
        FINAL_OUT --> CLEANUP
        CLEANUP --> END
    end

    style INIT fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
    style EVOLUTION fill:#fff9c4,stroke:#f57c00,stroke-width:2px
    style TIMESTEP fill:#ffffff,stroke:#666,stroke-width:1px
    style CHECKS fill:#ffffff,stroke:#666,stroke-width:1px
    style FINALIZE fill:#efebe9,stroke:#5d4037,stroke-width:2px
    
    click BUILD "../modules/mesh.html"
    click PGEN_INIT "../modules/pgen.html"
    click TASKS_EXEC "../modules/tasklist.html"
    click WRITE "../modules/outputs.html"
    click REFINE "../modules/mesh.html"
```

## Task Execution Detail

```{mermaid}
flowchart TD
    subgraph TaskStage[Single RK Stage]
        START_STAGE[Stage Start] --> RECV[Start MPI Recv]
        
        subgraph PhysicsTasks[Physics Tasks - Parallel]
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
    
    click HYDRO_TASK "../modules/hydro.html"
    click MHD_TASK "../modules/mhd.html"
    click RAD_TASK "../modules/radiation.html"
    click CT "../modules/mhd.html"
    click BC "../modules/boundaries.html"
```

## Data Flow

```{mermaid}
flowchart LR
    subgraph MeshBlock[MeshBlock Data]
        CONS[Conservative Variables]
        PRIM[Primitive Variables]
        BFIELD[Magnetic Field]
        METRIC[Metric Data]
    end

    subgraph Computation[Computation]
        C2P[Cons to Prim]
        RECONSTRUCT[Reconstruction]
        RIEMANN_SOLVE[Riemann Solver]
        FLUXES[Fluxes]
    end

    subgraph Communication[MPI Communication]
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
    
    click RECONSTRUCT "../modules/reconstruction.html"
    click RIEMANN_SOLVE "../modules/riemann_solvers.html"
    click PACK "../modules/boundaries.html"
    click SEND "../modules/boundaries.html"
```

## Memory Layout

```{mermaid}
flowchart TD
    subgraph Global[Global Memory]
        MESH_TREE[Mesh Tree]
        TASK_LIST[Task List]
        PARAMS[Parameters]
    end

    subgraph PerRank[Per MPI Rank]
        MB_PACK[MeshBlockPack]
        
        subgraph PerBlock[Per MeshBlock]
            COORDS_DATA[Coordinates]
            CONS_VARS[Conservative Vars]
            PRIM_VARS[Primitive Vars]
            FLUX_ARRAYS[Flux Arrays]
            SCRATCH[Scratch Space]
        end
        
        MPI_BUFFERS[MPI Buffers]
        OUTPUT_BUFFERS[Output Buffers]
    end

    subgraph Device[Device Memory - GPU]
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
    
    click MESH_TREE "../modules/mesh.html"
    click TASK_LIST "../modules/tasklist.html"
    click COORDS_DATA "../modules/coordinates.html"
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
    
    click MESH "../modules/mesh.html"
    click COORDS "../modules/coordinates.html"
    click BVALS "../modules/boundaries.html"
    click DRIVER "../modules/driver.html"
    click TASKS "../modules/tasklist.html"
    click OUTPUTS "../modules/outputs.html"
    click HYDRO "../modules/hydro.html"
    click MHD "../modules/mhd.html"
    click RAD "../modules/radiation.html"
    click Z4C "../modules/z4c.html"
    click PARTICLES "../modules/particles.html"
    click EOS "../modules/eos.html"
    click RECON "../modules/reconstruction.html"
    click RIEMANN "../modules/riemann_solvers.html"
    click DIFF "../modules/diffusion.html"
    click GRMHD "../modules/dyn_grmhd.html"
    click SRCTERMS "../modules/srcterms.html"
    click SHEARBOX "../modules/shearing_box.html"
    click PGEN "../modules/pgen.html"
```

## Performance Scaling

```{mermaid}
flowchart LR
    subgraph Hardware[Hardware Detection]
        CPU[CPU Cores]
        GPU[GPU Available?]
        MPI[MPI Ranks]
    end

    subgraph Execution[Execution Strategy]
        SERIAL[Serial]
        OPENMP[OpenMP]
        CUDA[CUDA]
        MPI_EXEC[MPI Parallel]
    end

    subgraph Optimization[Optimizations]
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
    
    click MPI_EXEC "../modules/boundaries.html"
    click LOAD_BAL "../modules/mesh.html"
```

## See Also
- [Module Index](../modules/index.md)
- [Configuration Guide](../configuration.md)
- [Running Simulations](../running.md)