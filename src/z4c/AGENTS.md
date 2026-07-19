<!-- BEGIN build-memory-table -->
# Z4c Numerical Relativity Navigation

## Scope

This subtree owns Z4c/BSSN spacetime state, RHS and gauge evolution, algebraic and ADM
constraints, Z4c-specific AMR behavior, compact-object tracking, waveform extraction,
CCE, and horizon dumps. Cross-module task dependency assembly lives in
`src/tasklist/numerical_relativity.*`; dynamical matter evolution lives in
`src/dyn_grmhd/` and supplies stress-energy through `Tmunu`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change evolved fields, options, or allocation | `src/z4c/z4c.hpp`, `src/z4c/z4c.cpp` | Output/restart variable layouts and pgens |
| Change evolution equations or finite-difference order | `src/z4c/z4c_calcrhs.cpp` | `src/utils/finite_diff.hpp`, ghost-zone configuration, convergence tests |
| Change gauge, RK update, or timestep | `z4c_gauge.cpp`, `z4c_update.cpp`, `z4c_newdt.cpp` | Driver integrators and task ordering |
| Change task ordering or matter coupling | `src/z4c/z4c_tasks.cpp` | `src/tasklist/numerical_relativity.*`, `src/dyn_grmhd/` |
| Change boundaries or AMR interpolation | `src/bvals/physics/z4c_bcs.cpp`, `src/z4c/z4c_amr.*` | `src/bvals/AGENTS.md`, `src/mesh/AGENTS.md` |
| Change ADM conversion or constraints | `src/z4c/z4c_adm.cpp`, `src/coordinates/adm.*` | `src/z4c/tmunu.*`, derived outputs |
| Change wave or compact-object diagnostics | `z4c_wave_extr.cpp`, `compact_object_tracker.*`, `horizon_dump.*`, `cce/` | Input blocks, restart state, output readers |
| Change initial data | Z4c pgens in `src/pgen/` | `inputs/z4c/`, `tst/inputs/`, external-library CMake branches if used |

## Important flow

- `MeshBlockPack::AddPhysics` creates Z4c and ADM state, optional dynamical GRMHD and
  `Tmunu`, then asks `NumericalRelativity` to assemble one dependency graph from the
  tasks queued by each active module.
- A stage copies Z4c state, calculates the RHS, applies boundary RHS handling and RK
  update, exchanges/restricts/prolongates state, enforces algebraic constraints, converts
  to ADM variables, and updates the timestep.
- End-of-integrator tasks calculate ADM constraints and Weyl fields, complete their AMR
  exchange, extract waves, update compact-object trackers, and run CCE/horizon outputs.
- Matter-coupled runs make stress-energy production and spacetime RHS dependencies
  explicit in `numerical_relativity.*`; do not infer ordering from source-file order.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/nr/test_nr_lwave1d_cpu.py`:
  CPU numerical-relativity and dynamical-matter evolution.
- `cd tst && python3 run_test_suite.py --test test_suite/z4c/test_z4c_lwave2d_amr_mpicpu.py`:
  Z4c finite-difference convergence with MPI AMR.
- GPU-specific changes should additionally use the closest `_gpu.py` case under
  `tst/test_suite/z4c/` or `tst/test_suite/nr/` with the required Kokkos architecture
  flags.

## Local constraints

- The configured finite-difference order, required ghost width, templated RHS/constraint
  kernels, and Z4c boundary interpolation order must remain consistent.
- When adding or moving a task, declare its required and optional dependencies in the NR
  queue; missing or cyclic dependencies abort task-list assembly.
- Changes to evolved field count/order propagate to ADM aliases, boundary buffers, AMR,
  outputs, and restart serialization.
- Keep diagnostics that depend on fully exchanged Weyl/ADM data after their communication
  and prolongation tasks.
<!-- END build-memory-table -->
