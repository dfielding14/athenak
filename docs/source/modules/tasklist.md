# Module: TaskList

## Overview
The TaskList module implements AthenaK's task-based execution system, managing computational dependencies and enabling efficient parallel execution on CPUs and GPUs.

## Source Location
`src/tasklist/`

## Key Components

| File | Purpose | Key Classes/Functions |
|------|---------|----------------------|
| `task_list.hpp` | Core task system definitions | `TaskList`, `TaskFunction` |
| `numerical_relativity.hpp` | NR-specific task management | `NumericalRelativityTasks` |
| `numerical_relativity.cpp` | NR task implementations | Task registration |

## Task System Architecture

### Task States
- **Incomplete**: Task not yet ready to run
- **Ready**: All dependencies satisfied
- **Running**: Currently executing
- **Complete**: Successfully finished

### Task Dependencies
Tasks can depend on:
- Completion of previous tasks
- Boundary communication
- Data availability
- Iteration count

## Task Registration

```cpp
// Example from hydro_tasks.cpp
TaskID none;
id.irecv = tl["before_stagen"]->AddTask(&Hydro::InitRecv, this, none);
id.copyu = tl["stagen"]->AddTask(&Hydro::CopyCons, this, none);
id.flux  = tl["stagen"]->AddTask(&Hydro::Fluxes, this, id.copyu);
id.recvf = tl["stagen"]->AddTask(&Hydro::RecvFlux, this, id.flux);
id.rkupdt = tl["stagen"]->AddTask(&Hydro::RKUpdate, this, id.recvf);
```

## Execution Flow

```{mermaid}
flowchart TD
    Start[Driver::Execute] --> Check{Check Tasks}
    Check --> Ready[Ready Tasks]
    Ready --> Launch[Launch Kernel]
    Launch --> Device{Device Type}
    Device -->|GPU| GPU[Kokkos::parallel_for<Cuda>]
    Device -->|CPU| CPU[Kokkos::parallel_for<OpenMP>]
    GPU --> Complete[Mark Complete]
    CPU --> Complete
    Complete --> Update[Update Dependencies]
    Update --> Check
```

## Task Categories

### Integration Stages
- `"before_stagen"`: Pre-stage setup
- `"stagen"`: Main computation stage
- `"after_stagen"`: Post-stage operations

**Note**: Task lists are accessed via string keys, not enums. The actual stage number is passed as an argument to task functions.

### Physics Tasks
Each physics module registers its tasks:

| Module | Key Tasks | Dependencies |
|--------|-----------|--------------|
| Hydro | StartRecv, Fluxes, Update | Boundary data |
| MHD | CT_Update, CornerE | Electric fields |
| Radiation | Source, Fluxes, Update | Opacity data |
| Particles | Push, Deposit, Interpolate | Field data |

## Performance Optimization

### Task Granularity
- Fine-grain: One task per MeshBlock operation
- Coarse-grain: Batch multiple blocks per task
- Adaptive: Runtime adjustment based on hardware

### Load Balancing
```cpp
// Task execution with load balancing
if (ntasks_per_rank > threshold) {
  // Use finer granularity
  ExecuteTasksParallel();
} else {
  // Sequential execution
  ExecuteTasksSerial();
}
```

## Configuration

No direct input parameters. Task behavior controlled by:
- Hardware detection (CPU vs GPU)
- MeshBlock count
- MPI rank distribution

## Integration with Other Modules

### Driver Module
- Driver calls `ExecuteTaskList()` for each RK stage
- Manages task list for entire simulation

### Physics Modules
- Each module adds its tasks during initialization
- Tasks operate on MeshBlockPack data

### MPI Communication
- Boundary tasks handle MPI message passing
- Overlaps computation with communication

## Common Task Patterns

### Boundary Exchange
```cpp
// Pattern for boundary communication
id.irecv = tl["before_stagen"]->AddTask(&Module::InitRecv, this, none);
id.compute = tl["stagen"]->AddTask(&Module::Compute, this, none);
id.recvu = tl["stagen"]->AddTask(&Module::RecvU, this, id.compute);
id.bcs = tl["stagen"]->AddTask(&Module::ApplyPhysicalBCs, this, id.recvu);
```

### Multi-Stage Update
```cpp
// Tasks execute in dependency order
id.copyu = tl["stagen"]->AddTask(&Module::CopyCons, this, none);
id.flux = tl["stagen"]->AddTask(&Module::Fluxes, this, id.copyu);
id.update = tl["stagen"]->AddTask(&Module::RKUpdate, this, id.flux);
// Driver calls this task list for each RK stage
```

## Debugging Tasks

### Task Timing
Enable timing output to identify bottlenecks:
```cpp
// In driver output
Task: Hydro::Fluxes - 0.045s (23%)
Task: MHD::CT_Update - 0.032s (16%)
```

### Dependency Visualization
Tasks can output dependency graphs for analysis.

## Example: Adding Custom Task

```cpp
// In your physics module constructor
void MyModule::MyModule(MeshBlockPack *ppack) {
  // Get task list map
  auto &tl = ppack->pmesh->pmb_pack->tl_map;
  
  // Register tasks with dependencies
  TaskID none;
  id.compute = tl["stagen"]->AddTask(&MyModule::MyComputation, this, none);
  id.update = tl["stagen"]->AddTask(&MyModule::MyUpdate, this, id.compute);
}

// Task implementation
TaskStatus MyModule::MyComputation(MeshBlockPack *ppack) {
  // Kokkos parallel execution
  auto &indcs = ppack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  
  par_for("MyComp", DevExeSpace(), 0, nmb-1,
          is, ie, js, je, ks, ke,
  KOKKOS_LAMBDA(int m, int i, int j, int k) {
    // Computation kernel
  });
  
  return TaskStatus::complete;
}
```

## See Also
- [Driver Module](driver.md)
- [Architecture Overview](../flowcharts/runtime.md)
- Source: `src/tasklist/task_list.hpp`