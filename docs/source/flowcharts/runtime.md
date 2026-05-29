# Runtime Flow

These diagrams summarize control flow implemented in `src/main.cpp`,
`src/driver/driver.cpp`, and `src/mesh/mesh_refinement.cpp`.

## Startup And Shutdown

```{mermaid}
flowchart TD
    A["MPI init when enabled"] --> B["Kokkos::initialize"]
    B --> C["Parse CLI options and paths"]
    C --> D["Load input or restart parameters"]
    D --> E["Apply block/name=value overrides"]
    E --> F{"-n requested?"}
    F -->|Yes| G["Dump parsed parameters and exit"]
    F -->|No| H["Construct Mesh and mesh tree"]
    H --> I{"-m requested?"}
    I -->|Yes| J["Report mesh structure and exit"]
    I -->|No| K["Add coordinates and enabled physics"]
    K --> L["Construct ProblemGenerator"]
    L --> M["Construct Driver and Outputs"]
    M --> N["Driver::Initialize"]
    N --> O["Driver::Execute"]
    O --> P["Driver::Finalize"]
    P --> Q["Delete runtime objects"]
    Q --> R["Kokkos::finalize; MPI_Finalize when enabled"]
```

`Driver::Initialize()` sets boundary values and primitives, computes the
initial timestep for non-static problems, and writes initial outputs for a new
run. Restart runs do not emit the new-run initialization outputs.

## Time Loop

```{mermaid}
flowchart TD
    LOOP{"time/cycle/wall limit not reached"} -->|"Continue"| DIAG["Print cycle diagnostics"]
    DIAG --> PRE["Execute before_timeintegrator tasks"]
    PRE --> STAGE["For each explicit stage: before_stagen, stagen, after_stagen"]
    STAGE --> POST["Execute after_timeintegrator tasks"]
    POST --> ADV["Advance time and cycle counters"]
    ADV --> OUT["Write due output streams"]
    OUT --> AMR{"Adaptive mesh active?"}
    AMR -->|Yes| REF["MeshRefinement::AdaptiveMeshRefinement"]
    AMR -->|No| DT["Mesh::NewTimeStep"]
    REF --> DT
    DT --> LOOP
    LOOP -->|"Stop"| FINAL["Driver::Finalize writes final outputs"]
```

Outputs are directed by the driver after a completed cycle; they are not
TaskList nodes. AMR is applied after due outputs, and the next timestep is
computed after refinement or derefinement has completed.

## Adaptive Refinement Entry Point

When `<mesh_refinement>/refinement = adaptive`, the driver calls
`MeshRefinement::AdaptiveMeshRefinement()`. That routine begins with
`CheckForRefinement()`, changes the tree and pack allocation if required,
reestablishes module timestep estimates as needed, and returns before
`Mesh::NewTimeStep()` selects the next global step.

See [Mesh](../modules/mesh.md), [Driver](../modules/driver.md),
[Outputs](../modules/outputs.md), and [TaskList](../modules/tasklist.md) for
the relevant interfaces and configuration.
