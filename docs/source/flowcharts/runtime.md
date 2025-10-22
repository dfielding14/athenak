# AthenaK Runtime Flow

This interactive flowchart shows the complete execution flow of an AthenaK simulation from startup to completion. Click on any node to navigate to detailed documentation.

## Main Execution Flow

```{mermaid}
flowchart TD
    Start([Start]) --> InitMPI[Initialize MPI (if enabled)]
    InitMPI --> InitKokkos[Initialize Kokkos]
    InitKokkos --> ParseArgs[Parse Command Line]
    ParseArgs --> ReadInput[Read .athinput (ParameterInput)]

    ReadInput --> BuildMesh[Construct Mesh & MeshBlocks]
    BuildMesh --> BuildPacks[Assemble MeshBlockPacks]
    BuildPacks --> AddPhysics{Register Physics Modules}

    AddPhysics --> |Hydro| AddHydro[Attach Hydro]
    AddPhysics --> |MHD| AddMHD[Attach MHD]
    AddPhysics --> |Radiation| AddRad[Attach Radiation]
    AddPhysics --> |Particles| AddPart[Attach Particles]
    AddPhysics --> |GR| AddGR[Attach Z4c / Dyn GR]

    AddHydro --> CreateDriver
    AddMHD --> CreateDriver
    AddRad --> CreateDriver
    AddPart --> CreateDriver
    AddGR --> CreateDriver

    CreateDriver[Create Driver & Outputs] --> DriverInit[Driver::Initialize]
    DriverInit --> ProblemInit[Problem Generator Sets ICs]
    DriverInit --> BuildTaskList[Build Task Lists]

    BuildTaskList --> TimeLoop{Time Evolution Loop}

    TimeLoop --> CheckTime{t < tlim?}
    CheckTime --> |Yes| ExecuteTasks[Execute Task Lists]
    CheckTime --> |No| Finalize[Finalize]

    ExecuteTasks --> UpdateDt[New Timestep]
    UpdateDt --> AMRCheck{AMR Triggered?}

    AMRCheck --> |Yes| RefineGrid[Refine / Derefine Grid]
    AMRCheck --> |No| OutputCheck
    RefineGrid --> Redistribute[Load Balance]
    Redistribute --> OutputCheck

    OutputCheck{Output due?}
    OutputCheck --> |Yes| WriteOutput[Write Outputs / Restart]
    OutputCheck --> |No| TimeLoop
    WriteOutput --> TimeLoop

    Finalize --> CleanupKokkos[Finalize Kokkos]
    CleanupKokkos --> CleanupMPI[Finalize MPI]
    CleanupMPI --> End([End])

    click BuildMesh "../modules/mesh.html" "Go to Mesh Module"
    click BuildPacks "../modules/mesh.html" "Go to Mesh Module"
    click AddHydro "../modules/hydro.html" "Go to Hydro Module"
    click AddMHD "../modules/mhd.html" "Go to MHD Module"
    click AddRad "../modules/radiation.html" "Go to Radiation Module"
    click AddPart "../modules/particles.html" "Go to Particles Module"
    click AddGR "../modules/z4c.html" "Go to Z4c Module"
    click CreateDriver "../modules/driver.html" "Go to Driver Module"
    click DriverInit "../modules/driver.html" "Go to Driver Module"
    click BuildTaskList "../modules/tasklist.html" "Go to TaskList Module"
    click ExecuteTasks "../modules/tasklist.html" "Go to TaskList Module"
    click RefineGrid "../modules/mesh.html" "Go to Mesh Module"
    click WriteOutput "../modules/outputs.html" "Go to Outputs Module"
```


## Task Execution Detail

```{mermaid}
flowchart LR
    subgraph TaskList[Task List Execution]
        Start2([Start Tasks]) --> CheckDep{Check Dependencies}
        CheckDep --> |Ready| Launch[Launch Kernel]
        CheckDep --> |Waiting| Wait[Wait]
        
        Launch --> Kernel{Kernel Type}
        Kernel --> |Hydro| HydroKernel[Hydro Update]
        Kernel --> |MHD| MHDKernel[MHD Update]
        Kernel --> |BC| BCKernel[Apply BCs]
        Kernel --> |Flux| FluxKernel[Compute Fluxes]
        
        HydroKernel --> Complete
        MHDKernel --> Complete
        BCKernel --> Complete
        FluxKernel --> Complete
        
        Complete[Mark Complete] --> UpdateDep[Update Dependencies]
        UpdateDep --> NextTask{More Tasks?}
        NextTask --> |Yes| CheckDep
        NextTask --> |No| EndTasks([End Tasks])
        
        Wait --> CheckDep
    end
    
    click HydroKernel "../modules/hydro.html" "Go to Hydro Module"
    click MHDKernel "../modules/mhd.html" "Go to MHD Module"
    click BCKernel "../modules/boundaries.html" "Go to Boundaries Module"
```

## AMR Workflow

```{mermaid}
flowchart TD
    subgraph AMR[Adaptive Mesh Refinement]
        StartAMR([Check Refinement]) --> EvalCrit[Evaluate Criteria]
        EvalCrit --> MarkBlocks[Mark Blocks]
        
        MarkBlocks --> Refine{Action?}
        Refine --> |Refine| SplitBlock[Split MeshBlock]
        Refine --> |Derefine| MergeBlock[Merge MeshBlocks]
        Refine --> |None| Skip[Skip]
        
        SplitBlock --> Prolongate[Prolongate Data]
        MergeBlock --> Restrict[Restrict Data]
        
        Prolongate --> UpdateTree
        Restrict --> UpdateTree
        Skip --> UpdateTree
        
        UpdateTree[Update Tree] --> UpdateNeighbors[Update Neighbors]
        UpdateNeighbors --> LoadBalance[Redistribute Blocks]
        LoadBalance --> EndAMR([AMR Complete])
    end
    
    click SplitBlock "../modules/mesh.html" "Go to Mesh Module"
    click MergeBlock "../modules/mesh.html" "Go to Mesh Module"
    click LoadBalance "../modules/mesh.html" "Go to Mesh Module"
```

## Key Components

### Initialization Phase
- **MPI**: Sets up parallel communication infrastructure
- **Kokkos**: Initializes performance-portable execution space
- **Parameter Input**: Reads configuration from `.athinput` files
- **Mesh Construction**: Builds domain decomposition and AMR tree

### Physics Modules
Each physics module adds its own:
- Conserved and primitive variables
- Flux computation routines
- Source terms
- Boundary conditions
- Task dependencies

### Time Evolution
The main time loop uses a task-based execution model:
1. Build task list with dependencies
2. Execute ready tasks in parallel
3. Update dependencies as tasks complete
4. Synchronize at stage boundaries

### Output System
Multiple output formats supported:
- VTK for visualization (ParaView/VisIt)
- Binary for efficient native I/O
- Restart files for checkpointing
- History for time series diagnostics
- Particle outputs (pvtk, trk, ppd)
- Coarsened binary for reduced data

### Adaptive Mesh Refinement
Dynamic mesh adaptation:
- User-defined refinement criteria
- Prolongation/restriction operators
- Load balancing across MPI ranks
- Flux correction at interfaces

## Performance Considerations

The runtime is optimized for:
- **GPU Efficiency**: MeshBlockPacks group blocks for coalesced access
- **Task Parallelism**: Overlapping computation and communication
- **Memory Locality**: Data structures optimized for cache/GPU memory
- **Load Balance**: Dynamic redistribution of MeshBlocks
- **State Preservation**: AMR rebuilds use module-specific hooks (e.g., `EnsureBasisSize` in the turbulence driver) to resize without reinitialising physics state

## See Also

- [Task System Module](../modules/tasklist.md)
- [Mesh Module](../modules/mesh.md)
- [Output Module - I/O](../modules/outputs.md)
- [Back to Documentation Home](../index.md)
