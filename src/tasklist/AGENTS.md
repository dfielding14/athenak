# AGENTS.md

## Purpose
This directory defines AthenaK's task-list infrastructure and the numerical
relativity (NR) task-queue system. It provides the generic `TaskList` scheduler
(`TaskID`, `Task`, `TaskList`) and the `NumericalRelativity` builder that queues
DynGRMHD and Z4c tasks into the `before_stagen`, `stagen`, and `after_stagen`
lists.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core task-list framework
- `task_list.hpp`: header-only definitions for `TaskID`, `Task`, and `TaskList`,
  plus `TaskStatus`/`TaskListStatus` enums and the dependency-checking logic.

### Numerical relativity task queuing
- `numerical_relativity.hpp` / `numerical_relativity.cpp`: `NumericalRelativity`
  class, `TaskName` enums for DynGRMHD/Z4c tasks, dependency filtering, and
  queue assembly into concrete `TaskList` instances.

---

## Task System Core (`task_list.hpp`)

### TaskStatus and TaskListStatus
- `TaskStatus`: `fail`, `complete`, `incomplete`.
- `TaskListStatus`: `running`, `stuck`, `complete`, `nothing_to_do`.
- `TaskList::DoAvailable` only returns `running` or `complete`.

### TaskID
- Bitset wrapper with `NUMBER_TASKID_BITS = 64` (64 dependency bits).
- Constructor `TaskID(id)`:
  - `id == 0` clears all bits (no dependencies).
  - otherwise sets bit `id-1`.
- `CheckDependencies(dep)` returns true when `bitfld_` contains all bits in `dep`.
- `SetComplete(rhs)` ORs `rhs` into the current bitfield.
- Supports `|`, `&`, `^`, and equality operators for dependency composition.

### Task
- Holds a `TaskID` (`myid_`), dependency `TaskID` (`dep_`), completion flag, and a
  `std::function<TaskStatus(Driver*, int)>`.
- `operator()(Driver*, int)` calls the bound function.
- `ChangeDependency(id, newdep)` replaces dependency bits `id` with `newdep` when
  `dep_` includes `id`.

### TaskList
- Stores tasks in a `std::list<Task>` and tracks completion with a `TaskID` bitset
  (`tasks_completed_`).
- `AddTask` overloads:
  - free/static function pointer
  - member function pointer + object
  - `std::function`
  - All assign `TaskID(size+1)` and append to the list.
- `InsertTask`:
  - Inserts a task *before* the task with ID `loc`.
  - New task ID is `TaskID(size+1)` (end-of-list index, not position).
  - Rewrites dependencies for **all other tasks**: any task whose dependency set
    includes the *dependency bitset of `loc`* (not the ID of `loc`) has that
    dependency replaced with the new task ID.
  - Returns `TaskID(0)` if `loc` is not found.
- `Reset` clears `tasks_completed_` and marks each task incomplete.
- `DoAvailable` executes any task whose dependencies are satisfied and that is not
  marked complete, and marks it complete only if it returns `TaskStatus::complete`.
  Tasks returning `fail` or `incomplete` remain pending.
- `IsComplete` checks `tasks_completed_` against all task IDs.

---

## NumericalRelativity Task Queues (`numerical_relativity.*`)

### Core types
- `TaskName` enum defines the NR task names used for dependency resolution:

  **DynGRMHD task names** (require `pmy_pack->pdyngr != nullptr`):
  - `MHD_Recv`, `MHD_CopyU`, `MHD_Flux`, `MHD_SetTmunu`, `MHD_SendFlux`,
    `MHD_RecvFlux`, `MHD_ExplRK`, `MHD_AddSrc`, `MHD_RestU`, `MHD_SendU`,
    `MHD_RecvU`, `MHD_EField`, `MHD_SendE`, `MHD_RecvE`, `MHD_CT`, `MHD_RestB`,
    `MHD_SendB`, `MHD_RecvB`, `MHD_BCS`, `MHD_Prolong`, `MHD_C2P`, `MHD_Newdt`,
    `MHD_ClearS`, `MHD_ClearR`, `MHD_NTASKS`

  **Z4c task names** (require `pmy_pack->pz4c != nullptr`):
  - `Z4c_Recv`, `Z4c_IRecvW`, `Z4c_CopyU`, `Z4c_CalcRHS`, `Z4c_SomBC`,
    `Z4c_ExplRK`, `Z4c_SendU`, `Z4c_RestU`, `Z4c_RecvU`, `Z4c_Newdt`, `Z4c_BCS`,
    `Z4c_Prolong`, `Z4c_AlgC`, `Z4c_Z4c2ADM`, `Z4c_Excise`, `Z4c_ADMC`,
    `Z4c_ClearS`, `Z4c_ClearR`, `Z4c_Weyl`, `Z4c_RestW`, `Z4c_SendW`,
    `Z4c_RecvW`, `Z4c_ProlW`, `Z4c_ClearSW`, `Z4c_ClearRW`, `Z4c_Wave`, `Z4c_PT`,
    `Z4c_NTASKS`

- `PhysicsDependency` enum:
  - `Phys_None`, `Phys_MHD`, `Phys_Z4c`.
  - `NeedsPhysics(task)` maps `task < MHD_NTASKS` to `Phys_MHD`,
    `task < Z4c_NTASKS` to `Phys_Z4c`, otherwise `Phys_None`.
- `TaskLocation` enum selects queue: `Task_Start`, `Task_Run`, `Task_End`.
- `QueuedTask` stores the task name, string label, queued/added state, assigned
  `TaskID`, dependency list, and bound function.

### Queueing and dependency filtering
- `QueueTask` overloads accept free/static functions or member functions with the
  signature `TaskStatus(Driver*, int)`.
- A task is queued only if all **required** dependencies are compatible with
  active physics (`DependenciesMet` + `DependencyAvailable`).
- `AddExtraDependencies` appends **optional** dependencies when their physics
  requirements are met.
- `QueueTask` does not validate ordering or cross-queue availability; missing or
  cyclic dependencies are detected later during list assembly.

### Assembly into TaskLists
- `AssembleNumericalRelativityTasks` (public) performs:
  - `pmy_pack->pdyngr->QueueDynGRMHDTasks()` when `pdyngr != nullptr`.
  - `pmy_pack->pz4c->QueueZ4cTasks()` when `pz4c != nullptr`.
  - Assembles three task lists: `before_stagen` (start queue), `stagen` (run queue),
    and `after_stagen` (end queue).
- `AssembleNumericalRelativityTasks(list, queue)` (private):
  - Iteratively adds tasks whose dependencies are already added **within the same
    queue**, building the dependency bitset from queued task IDs.
  - Returns `false` if a full pass adds no tasks (cyclic or missing dependencies).
- On assembly failure, it prints added vs missing tasks and aborts.

---

## Cautions and Constraints
- `NUMBER_TASKID_BITS` is 64; there is no guard against creating more than 64 tasks
  in a `TaskList`.
- `TaskList::DoAvailable` only returns `running` or `complete`; `stuck` and
  `nothing_to_do` are defined but not used.
- `NumericalRelativity` constructor accepts `ParameterInput *pin` but does not use it
  in the current implementation.
- `HasDependency` is declared in `numerical_relativity.hpp` but not defined in this
  directory.
