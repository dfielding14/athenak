# Module: Driver

## Role in AthenaK
The Driver owns AthenaK’s top-level time-integration loop. It decides whether the run is
static or time-evolving, configures the Runge–Kutta/ImEx coefficients, schedules the
shared task lists each stage, coordinates adaptive mesh refinement, and emits diagnostics
and outputs. Its responsibilities are split across `Initialize`, `Execute`, and `Finalize`,
each invoked from `main.cpp` after the mesh and physics modules have been constructed
([src/main.cpp:320], [src/driver/driver.hpp:22]).

## Source Layout
| File | Responsibility | Key APIs |
|------|----------------|----------|
| `driver.hpp` | Driver class declaration; integrator state, timers, counters | `Driver`, `ExecuteTaskList` |
| `driver.cpp` | Construction, lifecycle methods, scheduler, wall-clock handling | `Driver::Initialize`, `Driver::Execute`, `Driver::Finalize` |

## Execution Lifecycle
### Initialization
1. **Boundary priming** – `InitBoundaryValuesAndPrimitives` walks every enabled physics
   module, forcing the boundary communication tasks to run once with sentinel stage values
   (e.g., `stage = -1` suppresses flux receives, `stage = -4` clears shear buffers)
   ([src/driver/driver.cpp:549]).
2. **Module-specific timestep prep** – Any active hydro, MHD, radiation, or Z4c module is
   asked to compute its stage-local limit via `NewTimeStep(this, nexp_stages)`
   ([src/driver/driver.cpp:305]).
3. **Global timestep** – The mesh reduces those stage limits alongside diffusion and
   particle constraints in `Mesh::NewTimeStep`, clamping the first cycle to `tlim`
   ([src/mesh/mesh.cpp:580]).
4. **Initial outputs** – For fresh runs, every configured output writes once so the initial
   state is captured ([src/driver/driver.cpp:326]).
5. **ImEx scratch allocation** – If the ion-neutral module is present, the driver allocates
   `impl_src` and enforces that an ImEx integrator was selected ([src/driver/driver.cpp:340]).

### Main Loop
`Driver::Execute` exits early for static runs (currently a TODO) and otherwise advances
while time, cycle count, and optional wall-clock limits permit
([src/driver/driver.cpp:371]). Each cycle:
- Prints diagnostic timing every `ndiag` steps ([src/driver/driver.cpp:378]).
- Runs the shared task lists in order:
  `before_timeintegrator` @ stage 0, then for each explicit stage
  `before_stagen`, `stagen`, `after_stagen`, and finally
  `after_timeintegrator` @ stage 1 ([src/driver/driver.cpp:384]).
- Updates conserved time/cycle counters and the running meshblock / particle statistics
  ([src/driver/driver.cpp:396]).
- Triggers outputs when either their cadence `dt` or `dcycle` threshold is met
  ([src/driver/driver.cpp:413]).
- Invokes AMR (if enabled) and recomputes the next global timestep after refinement
  ([src/driver/driver.cpp:428]).
- Synchronises the wall-clock timer across ranks when a limit was provided
  ([src/driver/driver.cpp:374], [src/driver/driver.cpp:533]).

### Finalization
After the loop the driver flushes all outputs, calls any problem-specific final hook, and
prints aggregate performance metrics (zone-cycles, particle updates, AMR statistics). The
diagnostics are only emitted on rank 0 and include the reason for termination
([src/driver/driver.cpp:446]).

## Time Integration Options
The `<time>/integrator` string determines the explicit/implicit stage counts and CFL bound
([src/driver/driver.cpp:86]).

| Keyword | Order | Explicit stages (`nexp_stages`) | Implicit stages (`nimp_stages`) | CFL limit | Notes |
|---------|-------|---------------------------------|----------------------------------|-----------|-------|
| `rk1` | 1 | 1 | 0 | 1.0 | Forward Euler baseline ([src/driver/driver.cpp:91]). |
| `rk2` | 2 | 2 | 0 | 1.0 | SSPRK(2,2) Heun scheme ([src/driver/driver.cpp:99]). |
| `rk3` | 3 | 3 | 0 | 1.0 | SSPRK(3,3) ([src/driver/driver.cpp:112]). |
| `rk4` | 4 | 4 | 0 | 1.3925 | Low-storage RK4()4[2S]; uses `delta` weights for register updates ([src/driver/driver.cpp:129]). |
| `imex2` | 2 | 2 | 3 | 1.0 | Pareschi & Russo IMEX-SSP2; explicit leg matches RK2 ([src/driver/driver.cpp:162]). |
| `imex3` | 3 | 3 | 4 | 1.0 | IMEX-SSP3 with RK3 explicit leg ([src/driver/driver.cpp:188]). |
| `imex+` | 2 | 3 | 2 | 1.0 | Krapp et al. (2024) IMEX(2,3,2); explicit stages have non-RK2 weights ([src/driver/driver.cpp:230]). |

> The fatal error message printed for unknown integrators omits `rk4` and `imex+`; be aware
> the warning is stale ([src/driver/driver.cpp:258]).

## Driver Configuration
| Setting | Source | Default | Effect |
|---------|--------|---------|--------|
| `time/evolution` | input file | – (required) | Chooses `static`, `kinematic`, or `dynamic` evolution; anything else aborts ([src/driver/driver.cpp:68]). |
| `time/integrator` | input file | `rk2` | Only read for non-static runs; selects table entry above ([src/driver/driver.cpp:86]). |
| `time/tlim` | input file | – (required) | Simulation end time tested by the main loop ([src/driver/driver.cpp:274], [src/driver/driver.cpp:378]). |
| `time/nlim` | input file | `-1` | Maximum cycles; negative keeps running indefinitely ([src/driver/driver.cpp:88]). |
| `time/ndiag` | input file | `1` | Frequency for `OutputCycleDiagnostics` prints ([src/driver/driver.cpp:89], [src/driver/driver.cpp:513]). |
| `time/cfl_number` | input file | – (required) | Stored on the mesh and applied when reducing module timesteps ([src/mesh/build_tree.cpp:304], [src/mesh/mesh.cpp:586]). |
| `-t hh:mm:ss` | CLI flag | unset | Optional wall-clock limit broadcast each loop ([src/main.cpp:172], [src/driver/driver.cpp:374]). |

Static mode currently short-circuits without executing any task lists (marked TODO)
([src/driver/driver.cpp:371]).

## Task Scheduling Mechanics
`ExecuteTaskList` wraps each shared `TaskList`, clearing completion flags, running tasks
until dependencies are satisfied, and marking the pack complete before moving on
([src/driver/driver.cpp:273]). Every task sees the stage index that triggered it, which is
why physics modules use the integer to distinguish first/final stages or special sentinels.
The driver does not attempt to interleave task lists across packs; concurrency is provided
inside each task via Kokkos kernels.

## Boundary Priming Helper
`InitBoundaryValuesAndPrimitives` mirrors the actions performed during a live stage,
calling the same task wrappers with sentinel stages so that ghost zones, AMR restriction,
and shearing-box buffers start in a consistent state ([src/driver/driver.cpp:549]). The
routine is also reused by mesh refinement to seed newly created blocks.

## Diagnostics and Performance Counters
- `OutputCycleDiagnostics` prints the elapsed wall time, cycle, simulation time, and `dt`
  every `ndiag` cycles ([src/driver/driver.cpp:513]).
- `Execute` tallies meshblock and particle updates each cycle, then reports aggregate
  zone-cycles/second and particle updates/second during `Finalize`
  ([src/driver/driver.cpp:400], [src/driver/driver.cpp:493]).
- When AMR is enabled, `Finalize` also reports created/deleted meshblocks and the average
  load-balancing efficiency ([src/driver/driver.cpp:483]).

## Limitations & Notes
- Static runs (`time/evolution = static`) currently perform initialization and finalization
  but skip any analysis in `Execute` because the branch is still a placeholder
  ([src/driver/driver.cpp:371]).
- Ion-neutral two-fluid runs *must* use an ImEx integrator; the driver enforces this during
  initialization before allocating implicit scratch space ([src/driver/driver.cpp:340]).
- The global CFL number is read during mesh construction, so ensure `<time>/cfl_number` is
  present even if only the driver section is being edited ([src/mesh/build_tree.cpp:304]).

## See Also
- [TaskList Module](tasklist.md) – scheduler internals referenced by the driver.
- [Mesh Module](mesh.md) – timestep reduction, AMR, and mesh-wide counters.
- [Hydrodynamics Module](hydro.md) – example of per-stage task registration used by the driver.
