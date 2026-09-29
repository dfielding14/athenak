# AthenaK Runtime Flow

The runtime sequence below mirrors the control flow in `src/main.cpp` and the driver implementation. Nodes link to the relevant documentation for deeper inspection.

## Main Program Sequence

```{mermaid}
flowchart TB
    classDef setup fill:#e3f2fd,stroke:#1e88e5,stroke-width:2px,color:#0d47a1;
    classDef io fill:#ede7f6,stroke:#512da8,stroke-width:2px,color:#311b92;
    classDef mesh fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#1b5e20;
    classDef sim fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#e65100;
    classDef output fill:#fce4ec,stroke:#ad1457,stroke-width:2px,color:#880e4f;

    START([Start]):::setup --> MPI_INIT["Initialize MPI (if enabled)"]:::setup
    MPI_INIT --> KOKKOS_INIT["Initialize Kokkos"]:::setup
    KOKKOS_INIT --> PARSE["Parse CLI arguments"]:::setup
    PARSE --> VALIDATE{"Input or restart specified?"}:::io
    VALIDATE -->|No| ERROR["Abort: missing -i or -r"]:::io
    VALIDATE -->|Yes| LOAD_RST{"Restart file supplied?"}:::io
    LOAD_RST -->|Yes| LOAD_RESTART["Load restart header\n(ParameterInput::LoadFromFile)"]:::io
    LOAD_RST -->|No| SKIP_RST["Skip restart load"]:::io
    LOAD_RESTART --> LOAD_INPUT
    SKIP_RST --> LOAD_INPUT["Load athinput\n(ParameterInput::LoadFromFile)"]:::io
    LOAD_INPUT --> OVERRIDES["Apply CLI overrides\n(ModifyFromCmdline)"]:::io
    OVERRIDES --> DUMP_CHECK{"-n requested?"}:::io
    DUMP_CHECK -->|Yes| PAR_DUMP["Dump parameters and exit"]:::io
    DUMP_CHECK -->|No| MESH_STEP["Construct Mesh"]:::mesh

    MESH_STEP --> TREE_BUILD{"Restart run?"}:::mesh
    TREE_BUILD -->|No| BUILD_SCRATCH["Build tree from scratch"]:::mesh
    TREE_BUILD -->|Yes| BUILD_RESTART["Build tree from restart"]:::mesh
    BUILD_SCRATCH --> MESH_POST
    BUILD_RESTART --> MESH_POST["Set restart metadata"]:::mesh

    MESH_POST --> MESH_OPTION{"-m requested?"}:::mesh
    MESH_OPTION -->|Yes| WRITE_MESH["Write mesh structure and exit"]:::mesh
    MESH_OPTION -->|No| ADD_PHYS["Mesh::AddCoordinatesAndPhysics"]:::mesh

    ADD_PHYS --> PROBLEM_INIT["Construct ProblemGenerator\n(ICs or restart)"]:::sim
    PROBLEM_INIT --> DRIVER_BUILD["Create Driver"]:::sim
    PROBLEM_INIT --> OUTPUT_BUILD["Create Outputs manager"]:::output
    DRIVER_BUILD --> INITIALIZE["Driver::Initialize"]:::sim
    INITIALIZE --> EXECUTE["Driver::Execute"]:::sim
    EXECUTE --> FINALIZE["Driver::Finalize"]:::sim
    FINALIZE --> CLEANUP["Destroy Driver/Outputs/Mesh/ParameterInput"]:::sim
    CLEANUP --> KOKKOS_FINI["Finalize Kokkos"]:::setup
    KOKKOS_FINI --> MPI_FINI["Finalize MPI"]:::setup
    MPI_FINI --> END([End]):::setup

    click LOAD_RESTART "../modules/mesh.html" "Mesh restart handling"
    click LOAD_INPUT "../configuration.html" "Configuration guide"
    click MESH_STEP "../modules/mesh.html" "Mesh module"
    click ADD_PHYS "../modules/tasklist.html#assembly" "Coordinates and physics assembly"
    click PROBLEM_INIT "../modules/pgen.html" "Problem generator"
    click DRIVER_BUILD "../modules/driver.html" "Driver overview"
    click OUTPUT_BUILD "../modules/outputs.html" "Outputs manager"
```

## Driver Time Integration

```{mermaid}
flowchart LR
    classDef loop fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#e65100;
    classDef stage fill:#f3e5f5,stroke:#6a1b9a,stroke-width:2px,color:#4a148c;
    classDef check fill:#ede7f6,stroke:#512da8,stroke-width:2px,color:#311b92;

    INIT["Driver::Initialize"]:::loop --> CYCLE{"Cycle loop\n(t < tlim, n < nlim)"}:::loop
    CYCLE -->|No| DONE(["Driver::Finalize"]):::loop
    CYCLE -->|Yes| STAGE_LOOP{"For each RK/IMEX stage"}:::loop

    STAGE_LOOP --> POST_RECV["Post boundary receives"]:::stage
    POST_RECV --> DO_TASKS["Execute TaskList stage\n(module kernels)"]:::stage
    DO_TASKS --> POST_UPDATE["Update state & apply BCs"]:::stage
    POST_UPDATE --> NEWDT["Hydro/MHD/Radiation dt estimates"]:::stage
    NEWDT --> StageNext{"More stages?"}:::stage
    StageNext -->|Yes| STAGE_LOOP
    StageNext -->|No| OUTPUT_CHECK{"Output due?"}:::check

    OUTPUT_CHECK -->|Yes| WRITE["Outputs::MakeOutputs"]:::check
    OUTPUT_CHECK -->|No| AMR_CHECK{"Adaptive refinement triggered?"}:::check
    WRITE --> AMR_CHECK
    AMR_CHECK -->|Yes| ADAPT["Mesh::AdaptiveRefinement\n+ load balance"]:::check
    AMR_CHECK -->|No| ADVANCE
    ADAPT --> ADVANCE["Advance time & counters"]:::loop
    ADVANCE --> CYCLE

    click DO_TASKS "../modules/tasklist.html" "Task scheduling"
    click NEWDT "../modules/hydro.html#timestep-control" "Timestep control"
    click WRITE "../modules/outputs.html" "Outputs"
    click ADAPT "../modules/mesh.html#adaptive-mesh-refinement" "AMR workflow"
```

## AMR Workflow Snapshot

```{mermaid}
flowchart TB
    classDef amr fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#1b5e20;

    START_AMR[["AMR trigger"]]:::amr --> EVAL["Evaluate criteria\n(MeshRefinement::Evaluate)"]:::amr
    EVAL --> MARK["Mark MeshBlocks\n(refine / derefine)"]:::amr
    MARK --> ACTION{"Action"}:::amr
    ACTION -->|Refine| SPLIT["Split block\n(MeshBlockTree::Refine)"]:::amr
    ACTION -->|Derefine| MERGE["Merge blocks\n(MeshBlockTree::Derefine)"]:::amr
    ACTION -->|None| SKIP[["No change"]]:::amr
    SPLIT --> PROLONG["Prolongate data\n(module-specific operators)"]:::amr
    MERGE --> RESTRICT["Restrict data"]:::amr
    PROLONG --> SYNC_NGHBR[Update neighbour metadata]:::amr
    RESTRICT --> SYNC_NGHBR
    SKIP --> SYNC_NGHBR
    SYNC_NGHBR --> LB[Load balance / redistribute packs]:::amr
    LB --> AMR_DONE[[Return to Driver loop]]:::amr

    click PROLONG "../modules/mesh.html#prolongation" "Prolongation operators"
    click RESTRICT "../modules/mesh.html#restriction" "Restriction operators"
```

The main loop diagram follows the logic in `Driver::Execute`, while the AMR snapshot corresponds to the refinement helpers in `Mesh` and associated module hooks. Together with the system architecture overview, these charts provide a complete, code-accurate picture of how AthenaK evolves a simulation.
