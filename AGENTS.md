# AGENTS.md

## Purpose
This file orients automated agents working in AthenaK. It summarizes the codebase
layout and the expected approach to implementing changes. Use it as a map, but
verify details in the source before editing or documenting behavior.

---

## Codebase Overview

### What AthenaK is
AthenaK is a high-performance astrophysical simulation framework that targets
multi-physics magnetohydrodynamics, radiation transport, and numerical relativity
on block-structured AMR meshes. It is a Kokkos-based rewrite of Athena++ to run
efficiently on CPUs and GPUs.

### Execution flow (high level)
1. `src/main.cpp` initializes MPI/Kokkos, parses the input deck via
   `ParameterInput`, and constructs the `Mesh`.
2. `Mesh` owns the domain decomposition/AMR metadata and builds `MeshBlockPack`
   containers for GPU-friendly execution.
3. `MeshBlockPack::AddPhysics` instantiates only the physics modules requested
   in the input file.
4. `Driver` selects task lists and advances the simulation through integration
   stages, coordinating module updates and diagnostics.
5. `Outputs` schedules and writes history/restart/mesh/particle diagnostics
   before shutdown in `main.cpp`.

### Source tree map (core areas)
- `src/main.cpp` + `src/parameter_input.*`: startup, CLI, and input parsing.
- `src/AGENTS.md`: module map and top-level source orientation.
- `src/mesh/`: `Mesh`, `MeshBlock`, AMR tree/refinement, and `MeshBlockPack`.
- `src/driver/`: time integration, task scheduling, diagnostics.
- `src/tasklist/`: dependency-driven task lists for physics modules.
- `src/outputs/`: output formats and scheduling (`Outputs`, `OutputParameters`).
- `src/hydro/`, `src/mhd/`, `src/radiation/`, `src/dyn_grmhd/`, `src/z4c/`,
  `src/particles/`, `src/srcterms/`: physics packages.
- `src/reconstruct/` + `src/*/rsolvers/`: reconstruction and Riemann solvers.
- `src/bvals/`: boundary conditions and MPI exchanges.
- `src/eos/`: equations of state.
- `src/coordinates/`: geometry and metric helpers.
- `src/pgen/`: problem generators and regression scenarios.
- `inputs/`: runtime input decks and format notes (see `inputs/AGENTS.md`).
- `vis/python/`: Python post-processing utilities; see `vis/python/AGENTS.md`.
- `tst/`: regression tests and harness (`tst/run_tests.py`).

### Key abstractions
- `athena.hpp`: central typedefs, variable indices, Kokkos layouts.
- `ParameterInput`: block/parameter parsing for `.athinput` files.
- `MeshBlockPack`: batch container for per-block physics and task lists.
- `TaskList`: dependency graph for update stages.
- Kokkos `View` + `par_for`: required for performance-portable kernels.

---

## Implementation & Change Guidelines

### Approach to changes
- Be direct and pragmatic: small, incremental changes that build and test cleanly.
- Prefer clarity over cleverness; if it needs a long explanation, simplify.
- Study nearby code patterns first; mirror existing interfaces and data layouts.
- Avoid premature abstractions; add only what the current change needs.

### Adding or extending physics
- New modules typically require:
  - Constructors and storage in `MeshBlockPack`.
  - Registration in `MeshBlockPack::AddPhysics`.
  - TaskList entries for staging updates.
  - Input parsing via `ParameterInput` blocks.
  - Output variable wiring if new fields should be written.
- Keep kernels in Kokkos `par_for` with established memory layouts.

### Problem generators
- Implement `ProblemGenerator::UserProblem(ParameterInput*, bool restart)` in
  `src/pgen/`.
- Use `pmy_mesh_` and `pmy_mesh_->pmb_pack` for initialization.
- Always use Kokkos parallel patterns for loops.

### Documentation
- Docs in `docs/` are useful but may be stale. Verify claims in code before
  changing behavior or documentation.
- If editing docs, follow the audit workflow in `docs/AGENT_PRIMER.md`.

### PIC PR2 Coupling Status (Current)
- PR2 particle-to-MHD coupling is runtime opt-in under the `<particles>` block.
- Defaults preserve legacy behavior:
  - `couple_moments_to_mhd = false`
  - `couple_j_to_efield_representation = cell_centered`
  - fluid feedback toggles remain off by default.
- Coupled ordering now follows a staged split:
  - particle push/deposition wrappers are inserted in `stagen` (stage 1 only)
  - particle migration communication runs in `after_timeintegrator` when coupled.
- PR2 guards remain strict:
  - coupling is limited to single-fluid MHD task paths
  - radiation+MHD, hydro/ion-neutral, and NR (`adm`/`z4c`) compositions are
    fatal in coupled mode
  - `edge_staggered` current representation and fluid momentum/energy feedback
    are limited to non-relativistic MHD paths.
- Staged PIC runtime controls are now parsed under `<particles>` for the
  PR5 test-suite expansion path (`pic_background_mode`, `pic_feedback_mode`,
  `pic_deltaf_mode`, `pic_expanding_box_mode`, and related guard knobs). In
  the current stage:
  - `pic_background_mode=passive_mhd` is active and freezes fluid evolution in
    MHD update/source/CT tasks while preserving MHD state containers for PIC
    test-particle workflows.
  - `pic_background_mode=no_mhd` is active and uses a particle-owned magnetic
    field carrier (`pic_no_mhd_bx/by/bz`) for Boris pushers when no `<mhd>`
    block is present.
  - Boris CR pushers now use a midpoint E+B sequence (`cE = -u x B`) and
    store per-step particle delta channels used by coupled fluid feedback in
    `MHDSrcTerms`/`EFieldSrc` when `pic_feedback_mode=coupled`.
  - Step-5 diagnostics outputs are available for midpoint/feedback checks:
    `prtcl_dpxdt`, `prtcl_dpydt`, `prtcl_dpzdt`, `prtcl_dedt`, `prtcl_ebdot`.
  - Entity-mirroring deposit parity tests are now present in
    `particles/pic_entity_deposit_mink` and
    `particles/pic_entity_deposit_reflect` with serial/MPI decomposition checks.
  - An EM-vacuum-style analytic convergence anchor is present as
    `particles/pic_em_vacuum_wave` (linear-wave adaptation in test-particle mode).
  Other staged controls remain parse/validation hooks until their
  corresponding implementation steps.

---

## Build, Test, and Style

### Build (examples)
```bash
cmake -B build
cmake --build build -j8
```

### Tests
```bash
cd tst
python run_tests.py
```

### Style checks
```bash
cd tst/scripts/style
bash check_athena_cpp_style.sh
flake8 tst/ vis/
```

### Style reminders
- No tabs; 90-character line limit.
- One closing brace per line.
- Keep code boring and obvious.

---

## When in doubt
Prefer the source code over documentation. If behavior is unclear, stop and ask
for clarification rather than guessing.
