# AthenaK PIC PR1 Agent Implementation Guide

This guide is a companion to:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENTS.md`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`

## 1. Instruction Priority

1. Always follow `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENTS.md`.
2. Use `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
   as the PR1 scope/spec.
3. Use this guide for implementation style and source-selection policy.

If this guide conflicts with root `AGENTS.md`, root `AGENTS.md` wins.

## 2. Core Rule: Copy/Adapt First, Do Not Invent First

For PR1 PIC work, prefer adapting existing code from:
1. AthenaK (same module first, then similar module).
2. Entity (algorithm reference for PIC-specific behavior).

Do not start from blank implementations when equivalent patterns exist.

## 3. Source Selection Order

For each task, use this order:

1. Same-file AthenaK precedent
- Example: if editing CC boundary comm, first mirror structure in
  `src/bvals/bvals_cc.cpp` and `src/bvals/buffs_cc.cpp`.

2. Same-subsystem AthenaK precedent
- Example: for task wiring and communication clear ordering, mirror
  `src/particles/particles_tasks.cpp`.

3. Cross-subsystem AthenaK precedent
- Example: for additive accumulation style, mirror existing additive patterns in
  AthenaK boundary/flux code.

4. Entity reference
- Use Entity for PIC algorithm semantics, ordering, and synchronize behavior.
- Map Entity concepts onto AthenaK containers/tasks; do not transplant APIs
  blindly.

## 4. Required Implementation Workflow Per Change

For every non-trivial function:

1. Identify source anchor(s)
- Record exact source files used for adaptation (AthenaK and/or Entity).

2. Port structure first
- Copy control flow and data-access pattern before micro-optimizing.

3. Adapt types/indexing second
- Replace with AthenaK types (`DvceArray*`, task wrappers, boundary objects).
- Keep index and ghost-zone conventions consistent with AthenaK.

4. Preserve invariants
- Maintain existing AthenaK task ordering and communication completion rules.
- Keep compatibility for non-PIC code paths by default.

5. Validate immediately
- Compile after each major file group.
- Run the smallest relevant test before moving on.

## 5. PR1-Specific Guardrails

1. Do not add unnecessary abstractions.
- No new framework layer for PR1.
- Add only targeted enums/methods needed by the handoff plan.

2. Preserve existing behavior for non-PIC paths.
- New CC communication modes must default to current behavior.
- Existing hydro/mhd/radiation/z4c call sites should remain unchanged unless
  required by interface updates.

3. Keep task chain conservative.
- For moment communication, keep clear ordering as:
  `ClearRecvMoments -> ClearSendMoments`.

4. Keep runtime scope guard strict.
- `deposit_moments=true` requires PR1-supported configuration only.

## 6. Workstream-Specific Copy Targets

### WS-A (CC communication core)
- Primary AthenaK templates:
  - `src/bvals/buffs_cc.cpp`
  - `src/bvals/bvals_cc.cpp`
  - `src/bvals/bvals_tasks.cpp`
- Entity semantic reference:
  - synchronize/additive current communication behavior.

### WS-B (particles scaffolding and guards)
- Primary AthenaK templates:
  - `src/particles/particles.cpp`
  - `src/hydro/hydro.cpp`
  - `src/mhd/mhd.cpp`
- Follow existing AthenaK constructor/destructor allocation style.

### WS-C (deposition kernels and task DAG)
- Primary AthenaK templates:
  - `src/particles/particles_tasks.cpp`
  - `src/outputs/derived_variables.cpp` (`prtcl_d` atomic binning style)
- Entity reference:
  - deposition kernel ordering and old/new state usage.

### WS-D (output observability)
- Primary AthenaK templates:
  - `src/outputs/basetype_output.cpp` (derived var registration pattern)
  - `src/outputs/derived_variables.cpp` (derived var computation pattern)
  - `src/outputs/outputs.hpp` (valid variable names list)
- Required WS-D checklist:
  - update `NOUTPUT_CHOICES` and append
    `prtcl_rho/prtcl_jx/prtcl_jy/prtcl_jz` in `src/outputs/outputs.hpp`
  - add one derived registration block per new variable in
    `src/outputs/basetype_output.cpp`
  - enforce guard in `BaseTypeOutput`: these outputs require
    `particles/deposit_moments=true`
  - implement one `ComputeDerivedVariable()` branch per new variable in
    `src/outputs/derived_variables.cpp`
  - map each branch directly to `Particles::moments` component indices
    (`IMOM_RHO/JX/JY/JZ`)
  - keep `prtcl_d` unchanged and avoid refactoring unrelated output logic

### WS-E (tests and inputs)
- Primary AthenaK templates:
  - `tst/scripts/mhd/divb_amr.py` (run/analyze structure)
  - `tst/run_tests.py` discovery expectations
- Required WS-E checklist:
  - add `inputs/tests/pic_deposit_conservation.athinput` with a built-in
    problem generator (`pgen_name=linear_wave`) and PR1-supported particle
    deposition settings (`deposit_moments=true`, `deposit_order=1`,
    periodic boundaries, no AMR/multilevel/shearing-box)
  - include output variables `prtcl_rho`, `prtcl_jx`, `prtcl_jy`, `prtcl_jz`,
    and `prtcl_d` in the test deck
  - add `tst/scripts/particles/__init__.py` for test discovery
  - add `tst/scripts/particles/pic_deposit_conservation.py`:
    `run()` executes serial and MPI cases, `analyze()` validates global
    deposited `Q/J` integrals against expected charge/current invariants
  - add `tst/scripts/particles/pic_decomp_invariance.py`:
    compare global deposited `Q/J` integrals across at least two decomposition
    choices (and at least one MPI configuration)
  - parse mesh outputs in analysis using
    `vis/python/bin_convert_new.py` helpers for binary files
  - include explicit negative guard checks for unsupported configuration paths
    required by PR1 guardrails
  - ensure deterministic thresholds and clear logger diagnostics for measured
    values and errors
  - keep WS-E scope limited to test deck and regression scripts
  - keep tests deterministic and decomposition-aware.

### WS-F (PR1 closeout and merge readiness)
- Primary AthenaK templates:
  - `tst/run_tests.py` (targeted/package regression invocation)
  - `tst/scripts/style/check_athena_cpp_style_changed.sh` (C++ changed-file style gate)
  - `tst/scripts/style/check_python_style_changed.sh` (Python changed-file style gate)
  - existing `tst/scripts/particles/*.py` (PIC regression evidence format)
- Required WS-F checklist:
  - keep WS-F scope to validation hardening and closeout evidence only
  - do not add PR2 coupling or new PIC runtime features in WS-F
  - re-run MPI-enabled targeted PIC tests:
    `particles/pic_deposit_conservation` and
    `particles/pic_decomp_invariance`
  - re-run MPI-enabled `particles` package
  - run at least one non-MPI sanity execution path so serial behavior remains
    validated
  - run changed-file C++/Python style gates; fix only PR1-introduced issues
  - preserve deterministic pass/fail thresholds for global deposited `Q/J`
    integrals and guard-failure checks
  - produce a concise closeout summary with commands, outcomes, and residual
    risks carried to PR2/PR3

## 7. What Not To Do

1. Do not rewrite large existing subsystems to "clean up" style.
2. Do not change unrelated physics modules in PR1.
3. Do not silently alter default communication semantics.
4. Do not invent new naming/style patterns if AthenaK already has one.
5. Do not skip output wiring and tests for deposited moments.

## 8. Minimal Evidence Each Agent Must Provide

In each PR/patch description include:

1. Files changed.
2. AthenaK source patterns copied/adapted from.
3. Entity source semantics mirrored (if applicable).
4. Behavior preserved vs behavior intentionally changed.
5. Compile/test commands run and outcome.
6. For WS-E specifically:
   - serial and MPI run evidence
   - measured global deposited `Q/J` values and tolerance checks
   - decomposition-invariance comparison output.
7. For WS-F specifically:
   - repeated MPI evidence for targeted tests and package run
   - non-MPI sanity evidence
   - style gate outputs
   - explicit statement that no PR2-coupling files were changed

## 9. Definition of Done (per workstream)

1. Compiles on its own branch.
2. Does not break existing defaults.
3. Matches handoff scope for that stream only.
4. Includes at least one concrete validation relevant to that stream.
5. Is ready to merge in the order defined in
   `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`.
