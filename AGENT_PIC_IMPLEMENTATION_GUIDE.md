# AthenaK PIC PR1-PR5 Agent Implementation Guide

This guide is a companion to:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENTS.md`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_PR4_REVIEW_BRIEF.md`

## 1. Instruction Priority

1. Always follow `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENTS.md`.
2. Use `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
   as the active PIC scope/spec.
3. Use this guide for implementation style and source-selection policy.

If this guide conflicts with root `AGENTS.md`, root `AGENTS.md` wins.

## 1.1 Current Audit Snapshot (As-Built, 2026-02-07 Working Tree)

This section records the audited implementation state from source, not intent.

1. Runtime knobs and guards are implemented in particles constructor:
- coupling/representation/feedback knobs:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:168`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:216`
- PR1/PR2 guardrails:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:218`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:332`
- edge-current storage allocation in coupled edge mode:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:355`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:364`.

2. Coupled task ordering and insertion are implemented:
- migration moved to `after_timeintegrator` in coupled mode:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:52`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:63`
- stage insertion anchor follows feedback order:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:76`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:85`
- edge conversion insertion depends on both wrappers and `CornerE`:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:158`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:180`.

3. Moment wrappers and representation conversion are implemented:
- stage gating:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:16`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:25`
- edge-current zeroing:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:77`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:83`
- deterministic CC->edge conversion kernel:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:239`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:329`.

4. MHD coupling paths are implemented:
- fluid feedback in `MHDSrcTerms`:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:275`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:326`
- E-field source coupling branch by representation:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:427`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:513`
- alternate fluid-feedback placement in `EFieldSrc`:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:515`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:566`.

5. Regression coverage is present for both representation/order branches:
- coupled/uncoupled + guard + feedback checks:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:235`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:760`
- decomposition invariance across `cell_centered`/`edge_staggered` and
  `mhd_src_terms`/`efield_src`:
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:245`
  through
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:364`.

6. Entity alignment references used for audit:
- ordering invariants:
  `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:99`
  through
  `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:124`
- current-to-E coupling coefficient path:
  `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:563`
  through
  `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:612`
  and
  `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/ampere_mink.hpp:142`
  through
  `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/ampere_mink.hpp:205`.

## 2. Core Rule: Copy/Adapt First, Do Not Invent First

For PIC PR1/PR2 work, prefer adapting existing code from:
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

### WS-G (PR2 current-to-MHD coupling)
- Primary AthenaK templates:
  - `src/mhd/mhd_tasks.cpp` (`MHD::EFieldSrc` integration point)
  - `src/diffusion/resistivity.cpp` (edge-centered additive E-field loops)
  - `src/srcterms/turb_driver.cpp` (`TaskList::InsertTask` usage pattern)
  - `src/particles/particles_tasks.cpp` (task-chain wiring conventions)
  - `src/particles/particles_moments.cpp` (moment wrappers + guards)
- Required WS-G checklist:
  - keep coupling runtime-opt-in only with
    `particles/couple_moments_to_mhd=false` default behavior preserved
  - require `deposit_moments=true` for coupling mode and fail early otherwise
  - keep PR1 defaults unchanged for non-coupled runs and existing WS-E/WS-F tests
  - move moment scheduling toward stage-level parity by inserting moment tasks
    in `stagen` immediately before `MHD::EFieldSrc`
  - preserve deterministic PR2 behavior by executing stage-inserted moment
    deposition on stage 1 only
  - add particle-current contribution in `MHD::EFieldSrc` by mapping deposited
    `IMOM_JX/IMOM_JY/IMOM_JZ` into additive updates of `efld`
  - preserve existing shearing-box `EFieldSrc` behavior; particle term is
    additive, not a replacement
  - keep PR2 split explicit: no fluid momentum/energy feedback terms in WS-G
  - add deterministic test controls for non-zero current (e.g. CR velocity init
    knobs with zero defaults) to support coupled-field regression checks
  - add new regression deck/scripts for coupled-field behavior plus
    decomposition invariance with MPI coverage
  - include explicit negative checks for unsupported coupling configurations
  - keep WS-G scope limited to PR2 current-to-E coupling and tests only
    (no PR3 AMR/restart/non-periodic expansion).

### WS-G Post-Review Lock-Ins (2026-02-06)
- Add an explicit current-to-E coupling coefficient runtime knob in PR2 follow-up
  work, with default `1.0` so current behavior is preserved.
- Keep current WS-G ordering as a transitional PR2 state; defer full Entity-style
  particle-communication ordering parity to Step H.
- Expand WS-G negative tests to cover all guard branches
  (`radiation+MHD`, `ion-neutral/hydro`, and `adm/z4c`).
- Harden stage insertion by checking every `TaskList::InsertTask` return value,
  not only first/final insertions.
- Add determinism/linearity regression coverage for coupled mode (`cr_v*` and
  optional `deposit_qscale` sweeps).
- Canonical follow-up details live in
  `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md` sections
  `9.12` and `9.13`.
- Post-patch status note:
  explicit coupling coefficient and expanded WS-G tests are implemented.

### WS-H (PR2 fluid-feedback split + ordering movement)
- Primary AthenaK templates:
  - `src/mhd/mhd_tasks.cpp` (`MHDSrcTerms` and `EFieldSrc` additive source style)
  - `src/particles/particles_tasks.cpp` (task placement and insertion anchors)
  - `src/particles/particles.cpp` (runtime guards and option parsing)
- Required WS-H checklist:
  - keep all new behavior opt-in; preserve uncoupled/default paths
  - add independent momentum and energy feedback controls with coefficients
  - keep stage-1-only feedback application in coupled mode
  - enforce compatibility guards:
    - feedback requires `couple_moments_to_mhd=true`
    - energy feedback requires ideal MHD EOS
    - relativistic paths are rejected for fluid feedback in PR2
  - move coupled particle migration communication to
    `after_timeintegrator` while leaving uncoupled ordering unchanged
  - add regression checks proving each feedback branch changes the intended
    conserved quantity and leaves `Q/J` invariants intact.

### WS-I (PR2 representation parity slice + ordering anchor fix)
- Primary AthenaK templates:
  - `src/particles/particles_moments.cpp` (deterministic stencil conversion)
  - `src/particles/particles_tasks.cpp` (`InsertTask` dependency wiring)
  - `src/mhd/mhd_tasks.cpp` (representation branch in `EFieldSrc`)
- Required WS-I checklist:
  - add explicit representation selector
    `particles/couple_j_to_efield_representation` with default
    `cell_centered`
  - implement particle-owned edge-current arrays for converted `J`
  - zero edge-current arrays in the same lifecycle as moments
  - insert conversion immediately before `MHD::EFieldSrc` and ensure its
    dependency is after both:
    - inserted moment wrappers
    - `MHD::CornerE`
  - preserve existing `couple_j_to_efield_coeff` scaling path
  - add guard tests for unsupported `edge_staggered` relativistic use
  - extend decomposition tests over representation + feedback-order combinations.

### Current Deferred Items After WS-I
1. Full Entity-style charge-conserving staggered current deposition
   (trajectory-based deposit directly to staggered locations) is not yet
   implemented; AthenaK currently uses CC deposit + deterministic conversion.
2. PR3 scope is split into three small PRs:
   - PR3a: restart fidelity for coupled moment/edge-current state
   - PR3b: AMR/multilevel-safe deposited-moment handling
   - PR3c: non-periodic boundary policy for deposited moments/currents.

### PR3 Split Execution Rules (PR3a/PR3b/PR3c)
1. Keep PR3a/PR3b/PR3c scope-isolated.
   - PR3a must not include AMR or boundary-policy implementation work.
   - PR3b must not include restart-format redesign beyond what AMR needs.
   - PR3c must not include AMR algorithm work or restart-format expansion.
2. Preserve PR2 defaults and opt-in behavior.
   - default uncoupled paths remain unchanged
   - coupled paths change only within the active PR3x scope.
3. Require independent validation per PR3x.
   - PR3a: restart A/B equivalence (serial+MPI) for both representations
   - PR3b: multilevel/AMR conservation and decomposition invariance
   - PR3c: non-periodic boundary policy tests and guard coverage.
4. Maintain traceability in both docs.
   - each PR3x merge updates this guide and
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
     with as-built file:line references.

### PR3b As-Built Snapshot (Current Branch State)
1. Multilevel moment storage and wrapper API are now present in particle state/task
   wiring.
   - `src/particles/particles.hpp`:
     `ParticlesTaskIDs::{rest_mom, prol_mom}` and
     `Particles::{RestrictMoments, ProlongateMoments}` declarations.
2. Deposited moments allocate coarse mirrors when `mesh->multilevel` is enabled.
   - `src/particles/particles.cpp`: multilevel allocation of `coarse_moments`.
3. Moment wrapper order now follows AthenaK multilevel CC pattern:
   - `DepositMoments -> RestrictMoments -> SendMoments -> RecvMoments ->
     ClearRecvMoments -> ClearSendMoments -> ProlongateMoments`.
   - applied in both:
     - baseline `before_timeintegrator` path
     - coupled `stagen` insertion path.
4. Multilevel deterministic regression added:
   - input deck:
     `inputs/tests/pic_mhd_coupling_multilevel.athinput`
   - test:
     `tst/scripts/particles/pic_mhd_coupling_multilevel.py`
   - coverage includes:
     - `cell_centered` and `edge_staggered`
     - `mhd_src_terms` and `efield_src`
     - serial vs MPI decomposition invariance.

### PR3c As-Built Snapshot (Current Branch State)
1. Scope completion
   - PR3c is implemented as a non-periodic boundary-policy extension for deposited
     moments/currents only.
   - No restart-format expansion (PR3a scope) and no AMR algorithm expansion
     (PR3b scope) were added in this slice.
2. Implemented boundary policy and insertion points
   - Constructor policy guard implementation:
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:28`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:96`
     - invocation in constructor guard path:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:308`
   - Supported BC classes for deposited moments:
     - `periodic`
     - `outflow`
     - `reflect`
   - Unsupported BC classes fail fast with explicit face/flag diagnostics.
3. Implemented ordering contract
   - Added physical-BC wrapper:
     - declaration/task IDs:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:64`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:161`
     - implementation:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:253`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:265`
     - baseline/coupled task wiring:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:56`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:59`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:175`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:183`
   - Effective coupled ordering is:
     deposit/restrict/send/recv/clear -> apply moment BCs -> prolongate ->
     convert (edge mode) -> MHD source consumption.
4. Representation parity status
   - `cell_centered` and `edge_staggered` coupled modes run under the same
     boundary-policy contract.
   - Unsupported boundary classes are rejected before advance.
5. Regression evidence
   - New tests/decks:
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_mhd_coupling_nonperiodic.athinput`
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_mhd_coupling_guard_nonperiodic_unsupported.athinput`
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_nonperiodic.py`
   - Updated guards:
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py`
     - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_deposit_conservation.py`
   - Focused review runs passed:
     - `python run_tests.py particles/pic_mhd_current_coupling --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
     - `python run_tests.py particles/pic_mhd_coupling_nonperiodic --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
6. Entity alignment status
   - Sequencing parity remains aligned with Entity invariants:
     - `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:210`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:219`
     - `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:114`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:127`
   - Boundary taxonomy remains intentionally AthenaK-native.

### PR4 As-Built Status and PR4e Closeout Plan
1. As-built implementation status (`PR4a` through `PR4d`)
   - Runtime selector and defaults:
     - deposition-mode enum/default:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:40`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:111`
     - parser wiring:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:265`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:277`
   - Guardrails:
     - `direct_staggered` requires `edge_staggered`:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:384`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:392`
     - PR5-B removes the strict-periodic-only guard and applies explicit
       direct edge-current physical BCs before `MHD::EFieldSrc`:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:267`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:1181`
   - Dataflow and comm objects:
     - edge-comm pointer and task IDs:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:67`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:71`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp:145`
     - allocation/destruction:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:460`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:463`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:483`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:484`
   - Stage ordering and direct trajectory use:
     - direct mode old-position capture stage gate:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:82`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:104`
     - direct edge-current sync insertion before `MHD::EFieldSrc`:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:198`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:265`
   - Direct deposition kernel:
     - kernel entry:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:263`
     - first-order zig-zag style weighting and edge-current atomics:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:300`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:588`
   - Comm semantics and coupling integration:
     - additive edge-current exchange via `SumBoundaryFluxes`:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:748`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:751`
     - conversion path retained and gated to `cc_convert`:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:267`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:291`
       and
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:799`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:800`
     - MHD consumer path unchanged:
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:437`
       through
       `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp:476`.
2. As-built PR4 validation coverage
   - coupled direct-mode positive and guard coverage:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:392`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:585`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:761`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:774`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:868`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:1141`
   - decomposition parity now includes `edge_direct`:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:23`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:26`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:362`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:512`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:526`
   - multilevel parity now includes `edge_direct`:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:21`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:148`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:231`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:266`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:306`
   - restart fidelity now includes `edge_direct`:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:19`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:176`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:207`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:272`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:288`
   - deck defaults expose deposition-mode knob for CLI overrides:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_mhd_restart_fidelity.athinput:60`
     and
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_mhd_coupling_multilevel.athinput:72`.
   - continuity-residual oracle coverage is now active for both coupled and
     decomposition suites:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:28`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:165`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:605`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:872`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:1172`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:18`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:164`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:388`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:543`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:590`.
   - continuity residual measurement is now robust to duplicate terminal outputs
     at identical time/cycle:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:176`
     and
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:175`.
3. Entity-reference contract for PR4 closeout
   - ordering and pipeline references:
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:210`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/AGENTS.md:219`
     and
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:114`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:132`
   - direct trajectory deposition references:
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:521`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:560`
     and
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:170`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:404`
   - deferred higher-order reference surface:
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:406`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:560`
     and
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/particle_shapes.hpp:906`
     through
     `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/particle_shapes.hpp:997`.
4. PR4e remaining work
   - Completed in this step: direct-mode continuity residual oracle
     (`d(rho)/dt + div(J)`) in:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:605`
     through
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:1183`
     and
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:388`
     through
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_decomp.py:603`.
   - Completed in this step: explicit direct-vs-conversion closeout checks on
     multilevel + restart suites:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:266`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_multilevel.py:306`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:272`,
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_restart_fidelity.py:288`
   - Final compatibility decision:
     - retain `cc_convert` as the default compatibility/debug mode for PR4.
     - keep `direct_staggered` runtime opt-in until higher-order direct
       deposition and broader boundary coverage are implemented and rebaselined.
   - Explicitly document current direct-mode boundary restriction at:
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:394`
     through
     `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:400`.
5. PR4e exit criteria
   - All coupled suites pass with MPI:
     - `particles/pic_mhd_current_coupling`
     - `particles/pic_mhd_coupling_decomp`
     - `particles/pic_mhd_restart_fidelity`
     - `particles/pic_mhd_coupling_multilevel`
     - `particles/pic_mhd_coupling_nonperiodic`
   - style gates pass:
     - `bash tst/scripts/style/check_athena_cpp_style_changed.sh`
     - `bash tst/scripts/style/check_python_style_changed.sh`
   - handoff/docs contain as-built anchors and explicit compatibility-policy outcome.

### PR5 Entity-Parity Plan (Tracked Checklist)

1. Objective and contract
   - PR5 extends direct edge-current deposition beyond first-order and follows
     Entity behavior as closely as AthenaK architecture allows.
   - Preserve ordering, trajectory-deposition, and coupling semantics unless
     divergence is explicitly documented and justified.
   - Keep temporary divergences behind runtime guards and tracked checklist
     items.
   - Entity references:
     - `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:114`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:132`
     - `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:521`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/engines/srpic.hpp:560`
     - `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:406`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/currents_deposit.hpp:560`
       and
       `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/particle_shapes.hpp:906`
       through
       `/Users/dbf75/Work/Research/AthenaK/entity/src/kernels/particle_shapes.hpp:997`

2. Tracked checklist
   - PR5-A: higher-order direct deposition kernel parity
     - [x] Port higher-order Entity-style trajectory deposition into AthenaK
       direct path without changing PR4 first-order behavior.
       - Status (2026-02-08, commit `56244bc9`): periodic direct mode supports
         `deposit_order=2`.
     - [x] Add runtime support for direct-mode deposition orders beyond order-1
       with clear fatal guards for unsupported combinations.
       - Status (2026-02-08, commit `56244bc9`): runtime guards permit
         `deposit_order={1,2}` only for coupled `direct_staggered`; other
         unsupported combinations remain fatal.
     - [x] Document exact AthenaK-to-Entity mapping for kernel primitives
       (weights, trajectory split, current accumulation).
       - Status (2026-02-08, baseline `4d2236be`): completed in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
         Section 9.20.3.3.
   - PR5-B: boundary-policy expansion for direct mode
     - [x] Implement conservation-safe direct trajectory handling for supported
       non-periodic policies (`reflect`, `outflow`).
     - [x] Remove strict-periodic direct-mode guard after PR5-B validation is
       green.
       - Status (2026-02-08, commit `f8d2642d`): constructor periodic guard removed;
         direct edge-current BC task added in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_moments.cpp:1181`
         and inserted before `MHD::EFieldSrc` in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_tasks.cpp:267`.
   - PR5-C: test and parity hardening
     - [x] Extend direct-mode coverage (higher-order and non-periodic) in:
       `pic_mhd_current_coupling.py`, `pic_mhd_coupling_decomp.py`,
       `pic_mhd_coupling_multilevel.py`, `pic_mhd_restart_fidelity.py`,
       `pic_mhd_coupling_nonperiodic.py`.
       - Status (2026-02-08, commit `f8d2642d`): non-periodic direct
         `direct_staggered` + `deposit_order=2` cases added in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_nonperiodic.py:227`
         and
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_coupling_nonperiodic.py:290`;
         obsolete periodic-only guard removed in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_mhd_current_coupling.py:827`.
     - [x] Keep continuity-residual oracle (`d(rho)/dt + div(J)`) as required.
     - [x] Add explicit Entity parity classification per axis: aligned,
       intentional divergence, unplanned divergence.
       - Status (2026-02-08, commit `f8d2642d`): matrix recorded in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
         Section 9.20.3.1; no unplanned divergence found in audited PR5 axes.
     - [x] Record baseline numeric envelopes for decomposition, restart,
       multilevel, and direct-vs-conversion deltas.
       - Status (2026-02-08, commit `f8d2642d`): envelope values recorded in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
         Section 9.20.3.2.
   - PR5-D: default-mode decision gate
     - [x] Evaluate switching default from `cc_convert` to
       `direct_staggered` only after PR5-A/B/C pass.
       - Status (2026-02-08, baseline `4d2236be`): adopt a
         context-aware default in
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:265`
         through
         `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp:274`
         (`edge_staggered -> direct_staggered`,
         `cell_centered -> cc_convert`).
     - [x] If any parity/stability/performance gate fails, keep `cc_convert`
       default and retain `direct_staggered` as opt-in with explicit rationale.
       - Status (2026-02-08): full MPI matrix is green and no blocking parity
         divergence was found; compatibility-driven context-aware default
         selected instead of a global direct default.

3. PR5 merge gates
   - Required MPI tests (existing five plus new PR5 additions) pass:
     - `particles/pic_mhd_current_coupling`
     - `particles/pic_mhd_coupling_decomp`
     - `particles/pic_mhd_restart_fidelity`
     - `particles/pic_mhd_coupling_multilevel`
     - `particles/pic_mhd_coupling_nonperiodic`
   - Changed-file style checks pass:
     - `bash tst/scripts/style/check_athena_cpp_style_changed.sh`
     - `bash tst/scripts/style/check_python_style_changed.sh`
   - PR5 closeout docs contain:
     - completed checklist state
     - explicit Entity-parity matrix
     - explicit default-mode decision outcome with evidence.
     - latest matrix log reference:
       `/tmp/pr5_closeout_matrix_20260207_212254.log`.

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
8. For WS-G specifically:
   - explicit evidence that coupling is opt-in and defaults are unchanged
   - serial and MPI coupled-field test evidence
   - decomposition-invariance evidence for coupled runs
   - measured `Q/J` diagnostics and coupled-field deltas used for pass/fail
   - explicit statement that fluid momentum/energy feedback was not added
     in WS-G scope.
9. For WS-H specifically:
   - evidence for independent momentum and energy feedback toggles
   - evidence that coupled communication ordering moved while default ordering
     remained unchanged
   - guard-coverage evidence for invalid feedback compositions/EOS.
10. For WS-I specifically:
   - evidence for both current-representation modes
     (`cell_centered`, `edge_staggered`)
   - evidence that conversion is inserted before `MHD::EFieldSrc` with
     deterministic dependencies
   - decomposition invariance evidence across representation and feedback-order
     combinations.
11. For PR4 specifically:
   - explicit mapping of AthenaK kernel/control-flow choices to Entity source
     references
   - continuity-residual and decomposition-invariance evidence for direct
     deposition mode
   - A/B comparison evidence versus converted-edge mode on shared test decks
   - explicit statement of default behavior and compatibility-mode status.

## 9. Definition of Done (per workstream)

1. Compiles on its own branch.
2. Does not break existing defaults.
3. Matches handoff scope for that stream only.
4. Includes at least one concrete validation relevant to that stream.
5. Is ready to merge in the order defined in
   `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`.
