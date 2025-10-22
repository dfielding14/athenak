# Module: Driver

## Overview
The Driver module controls the main time evolution loop, managing time integration schemes and coordinating task execution.

## Source Location
`src/driver/`

## Key Components

| File | Purpose | Key Classes/Functions |
|------|---------|----------------------|
| `driver.hpp` | Driver class definition | `Driver` class |
| `driver.cpp` | Main evolution loop | `Execute()`, `ExecuteTaskList()` |

## Time Integration Schemes

Configured via `integrator` parameter in `<time>` block:

| Integrator | Order | Stages | Description | Lines in driver.cpp |
|------------|-------|--------|-------------|---------------------|
| `rk1` | 1st | 1 | Forward Euler | L97-105 |
| `rk2` | 2nd | 2 | Heun's method (SSPRK2, default) | L105-118 |
| `rk3` | 3rd | 3 | Strong stability preserving (SSPRK3) | L118-135 |
| `rk4` | 4th | 4 | RK4()4[2S] from Ketcheson (2010) | L135-167 |
| `imex2` | 2nd | 3 | IMEX-SSP2(3,2,2) Pareschi & Russo | L167-194 |
| `imex3` | 3rd | 4 | IMEX-SSP3(4,3,3) Pareschi & Russo | L194-236 |
| `imex+` | 2nd | 3 | IMEX(2,3,2) Krapp et al. (2024) | L236-263 |

**Note**: Error message at L267 incorrectly lists only `[rk1,rk2,rk3,imex2,imex3]` but `rk4` and `imex+` are implemented.

## Configuration Parameters

### Time Block (`<time>`)

| Parameter | Type | Default | Description | Set in |
|-----------|------|---------|-------------|--------|
| `evolution` | string | required | Evolution type (static/kinematic/dynamic) | driver.cpp |
| `integrator` | string | rk2 | Time integration scheme | driver.cpp |
| `cfl_number` | Real | required | CFL stability factor | build_tree.cpp L306, L502 |
| `tlim` | Real | required | Simulation end time | driver.cpp |
| `nlim` | int | -1 | Maximum cycles (-1 = unlimited) | driver.cpp |
| `ndiag` | int | 1 | Cycles between diagnostics | driver.cpp |

**Note**: `cfl_number` is read in `mesh/build_tree.cpp`, not in `driver.cpp`. There is no fixed `dt` parameter.

## Evolution Types

Set via `evolution` parameter:

| Type | Description | Enum Value |
|------|-------------|------------|
| `static` | Single evaluation, no time integration | TimeEvolution::tstatic |
| `kinematic` | Fixed velocity field, passive evolution | TimeEvolution::tkinematic |
| `dynamic` | Full time-dependent evolution | TimeEvolution::tevolve |

## Key Methods

### Main Evolution Loop
```cpp
// driver.cpp L375
void Driver::Execute(Mesh *pmesh, ParameterInput *pin, Outputs *pout) {
  // Main time evolution loop
  while (pmesh->time < pmesh->tlim && 
         pmesh->ncycle < pmesh->nlim) {
    // Execute stages
    for (int stage = 1; stage <= nstages; ++stage) {
      ExecuteTaskList(pmesh, tl_name, stage);
    }
    // Update time
    pmesh->time += pmesh->dt;
    pmesh->ncycle++;
  }
}
```

### Task List Execution
```cpp
// driver.cpp L279
void Driver::ExecuteTaskList(Mesh *pm, std::string tl, int stage) {
  // Execute task list for given stage
  pm->pmb_pack->tl_map[tl]->Execute(pm, stage);
}
```

## Task Execution Stages

For each integration stage, the following task lists are executed:

1. `before_timeintegrator` - Prepare for time integration
2. `before_stageN` - Stage-specific preparation (N = stage number)
3. `stageN` - Main stage computation
4. `after_stageN` - Stage cleanup
5. `after_timeintegrator` - Finalize timestep

## Time Step Calculation

New timestep is calculated in `Mesh::NewTimeStep()` (mesh.cpp L571-621):
```cpp
dt = cfl_no * min(dt_hydro, dt_mhd, dt_z4c, dt_rad, ...);
```

## Performance Monitoring

The driver tracks (driver.cpp L440-500):
- Wall clock time per cycle
- Zone-cycles per second  
- Time spent in MPI communication
- Diagnostic output interval

## Integration with Other Modules

- **Mesh**: Provides spatial discretization and stores `cfl_no`, `dt`, `time`
- **TaskList**: Manages and executes computational tasks
- **Physics Modules**: Register tasks with task lists
- **Outputs**: Called by driver at specified intervals
- **AMR**: Mesh refinement checked between cycles

## Example Usage

### Standard Evolution
```ini
<time>
evolution = dynamic   # Full time evolution
integrator = rk2      # 2nd order Heun's method
cfl_number = 0.4      # CFL factor
tlim = 10.0          # Run to t=10
nlim = -1            # No cycle limit
ndiag = 10           # Diagnostics every 10 cycles
```

### Higher-Order Integration
```ini
<time>
evolution = dynamic
integrator = rk3      # 3rd order SSP-RK
cfl_number = 0.3     # More conservative CFL for higher order
tlim = 1.0
```

### Static Problem (No Evolution)
```ini
<time>
evolution = static    # Single evaluation only
tlim = 0.0           # Not used but required
cfl_number = 1.0     # Not used but required
```

## Common Issues

### Timestep Too Small
- CFL number may be too large for the problem
- Check for unresolved features or shocks
- Consider enabling FOFC in physics modules

### Integration Instability
- Reduce `cfl_number` (typical: 0.4 for RK2, 0.3 for RK3)
- Use SSP integrators (rk2, rk3) for problems with shocks
- Check boundary conditions

### Performance
- Adjust MeshBlock size for better GPU utilization
- Monitor load balance with `ndiag` diagnostics
- Consider lower-order integrator if accuracy permits

## See Also
- [TaskList Module](tasklist.md) - Task management system
- [Mesh Module](mesh.md) - Timestep calculation
- Source: `src/driver/driver.cpp`