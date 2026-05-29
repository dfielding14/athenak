# Module: TaskList

## Role in AthenaK
The TaskList module provides the lightweight scheduler that the driver uses to orchestrate
per-stage work across every `MeshBlockPack`. Each physics package contributes tasks to the
shared lists, while the `TaskList` class enforces ordering constraints, retries
communication-heavy work until it completes. The scheduler
itself is intentionally minimal: it executes tasks sequentially on the host while each task
invokes Kokkos kernels or MPI routines as needed.

## Source Layout
| File | Responsibility | Key APIs |
|------|----------------|----------|
| `task_list.hpp` | Core scheduler, task container, and dependency tracking | `TaskID`, `Task`, `TaskList` |
| `numerical_relativity.hpp/cpp` | Queues NR-specific tasks and resolves cross-physics dependencies | `numrel::NumericalRelativity` |

## Task Graph Lifecycle
- `MeshBlockPack::MeshBlockPack` creates the canonical task lists (`before_timeintegrator`,
  `after_timeintegrator`, `before_stagen`, `stagen`, `after_stagen`) and stores them in
  `tl_map`, keyed by string so modules can opt-in by name (`src/mesh/meshblock_pack.cpp:43`,
  `src/mesh/meshblock_pack.hpp:79`).
- During `MeshBlockPack::AddPhysics`, standalone Hydro or MHD registers its
  own tasks in the shared `tl_map`. Evolving radiation-fluid configurations
  (`fixed_fluid = false`), ion-neutral, and numerical relativity
  configurations instead assemble combined graphs for constructed fluid
  modules; radiation with `fixed_fluid = true` uses a transport-only branch.
  Returned `TaskID`s still allow dependent extensions such as turbulence
  forcing (`src/mesh/meshblock_pack.cpp:113`,
  `src/mesh/meshblock_pack.cpp:140`, `src/mesh/meshblock_pack.cpp:168`,
  `src/radiation/radiation_tasks.cpp:33`, `src/mesh/meshblock_pack.cpp:227`).
- At runtime the driver calls `Driver::ExecuteTaskList` for each stage, which repeatedly
  invokes `TaskList::DoAvailable` on every pack until the list reports completion
  (`src/driver/driver.cpp:273`, `src/tasklist/task_list.hpp:145`).

### Task Map Keys
| Key | Typical Responsibility |
|-----|------------------------|
| `before_timeintegrator` | Operator-split work once per cycle, including current turbulence mode initialization and one forcing application. |
| `before_stagen` | Stage-wide preparation such as posting boundary receives (`src/hydro/hydro_tasks.cpp:50`). |
| `stagen` | Main body of the per-stage update (fluxes, Runge–Kutta updates, source terms). |
| `after_stagen` | Cleanup once all packs finish the stage (clearing outstanding MPI requests). |
| `after_timeintegrator` | Post-timestep maintenance (rarely populated today). |

The driver forwards the current stage index to every task, allowing functions to gate work
on the first or final stage (`src/hydro/hydro_tasks.cpp:126`).

## Task and List Semantics
### Status Codes
Tasks return one of three `TaskStatus` values (`src/tasklist/task_list.hpp:29`):
- `complete`: the task finished and any dependent work may proceed.
- `incomplete`: dependencies remain external (e.g., non-blocking MPI has not completed);
  `TaskList::DoAvailable` will retry the task on the next sweep. MPI unpack routines use
  this to wait without busy-waiting (`src/bvals/bvals_cc.cpp:263`).
- `fail`: declared as a task return value, but `TaskList::DoAvailable()` currently
  marks only `complete` and otherwise continues reporting `running`. A new task
  should not rely on scheduler-level propagation of `fail` without implementing
  and testing that behavior.

`TaskListStatus` records list-wide progress. `DoAvailable` currently returns only
`running` or `complete`, while the `stuck` and `nothing_to_do` enumerators remain reserved
for future use (`src/tasklist/task_list.hpp:30`, `src/tasklist/task_list.hpp:156`).

### Dependency Encoding
`TaskID` wraps a 64-bit bitset that records which tasks have completed
(`src/tasklist/task_list.hpp:36`). When `AddTask` appends a new task, it assigns the next
bit, so dependencies are expressed as bitwise combinations of upstream IDs. `TaskList`
tracks a second bitset, `tasks_completed_`, that it updates whenever a task reports
`complete` (`src/tasklist/task_list.hpp:121`).

### Scheduler Loop
Calling `TaskList::Reset` clears completion flags before a new timestep or stage
(`src/tasklist/task_list.hpp:139`). `DoAvailable` then walks the list in order, executing
any task whose dependency bits are satisfied. Tasks returning `complete` are marked
finished immediately; tasks returning `incomplete` are retried on the next pass
(`src/tasklist/task_list.hpp:145`). Since there is no internal threading or GPU queueing,
the degree of parallelism is entirely determined by each task’s implementation (usually a
Kokkos kernel).

## Registering Tasks
Physics packages build their task graphs during construction and retain the resulting
`TaskID`s to make later insertions deterministic. Hydrodynamics provides the canonical
example:

```cpp
// src/hydro/hydro_tasks.cpp:47
TaskID none(0);
id.irecv = tl["before_stagen"]->AddTask(&Hydro::InitRecv, this, none);
id.copyu = tl["stagen"]->AddTask(&Hydro::CopyCons, this, none);
id.flux  = tl["stagen"]->AddTask(&Hydro::Fluxes, this, id.copyu);
id.sendf = tl["stagen"]->AddTask(&Hydro::SendFlux, this, id.flux);
id.recvf = tl["stagen"]->AddTask(&Hydro::RecvFlux, this, id.sendf);
id.rkupdt = tl["stagen"]->AddTask(&Hydro::RKUpdate, this, id.recvf);
```

The returned IDs are stored in `HydroTaskIDs` so that other systems (e.g., turbulence
forcing) can add work at named dependency points (`src/hydro/hydro.hpp:39`).
Every task implements `TaskStatus f(Driver*, int stage)` so it can access shared driver
state such as the current stage count or CFL buffers.

### Inserting Work Between Existing Tasks
`TaskList::InsertTask` splices a task before a target `TaskID` and rewrites downstream
dependencies to preserve ordering (`src/tasklist/task_list.hpp:197`). The turbulence
driver uses this to add forcing immediately before the Runge-Kutta update, dependent on
the corresponding flux task:

```cpp
// src/srcterms/turb_driver.cpp:446
tl->InsertTask(&TurbulenceDriver::AddForcing, this,
               pmy_pack->phydro->id.flux, pmy_pack->phydro->id.rkupdt);
```

### Numerical Relativity Choreography
`numrel::NumericalRelativity` maintains three queues (start/run/end). Each queued task
records optional dependencies on other NR tasks and on physics features such as Z4c or
DynGRMHD. `AssembleNumericalRelativityTasks` walks the queues, adds tasks when all
dependencies have been satisfied, and aborts with a helpful report if it detects a cycle
or a missing prerequisite (`src/tasklist/numerical_relativity.cpp:112`).

## Guidelines for New Tasks
- Always thread safety through Kokkos or MPI inside the task body; the scheduler itself is
  single-threaded.
- Return `TaskStatus::incomplete` while waiting on non-blocking communication or host
  events so subsequent tasks do not race.
- Capture the `TaskID` you receive from `AddTask`; downstream modules may rely on it to
  hook additional behaviour.
- Prefer adding tasks during module construction so the graph is immutable during the run.
  If dynamic behaviour is required, reuse `InsertTask` to preserve dependency ordering.

## Debugging and Inspection
- `TaskList::PrintIDs` and `PrintDependencies` dump the current graph (`src/tasklist/task_list.hpp:135`).
- `TaskID::PrintID` and `TaskList::GetIDLastTask` aid in debugging dynamic insertions
  (`src/tasklist/task_list.hpp:55`, `src/tasklist/task_list.hpp:133`).
- There is no built-in timing or dependency visualisation; use the driver diagnostics or
  instrument individual tasks if detailed profiling is required.

## See Also
- [Driver Module](driver.md) – invokes task lists each stage.
- [Mesh Module](mesh.md) – owns `MeshBlockPack` construction and timestep control.
- [Source Terms Module](srcterms.md) – additional examples of inserting specialised tasks.
