# MHD-PIC Active Development Tracker

> Temporary working document for `PIC_development`, simplified 2026-07-11.
> Detailed findings from the original review remain available in Git history at
> `9364e0ff3:MHD_PIC_REVIEW_TRACKER.md`.

## Purpose and baseline

This tracker records work that affects the next planned MHD-PIC science runs. It follows
the project principles in [ethos.md](ethos.md): physical consistency is non-negotiable,
while validation, hardening, and optimization should be proportional to demonstrated
risk and intended use.

- Branch: `PIC_development`
- Current evaluated HEAD: `9364e0ff37`
- Near-term workflows: uniform-grid Bell, non-relativistic shocks, and turbulent boxes
- Detailed roadmap: `MHD_PIC_NEXT_STEPS_GUIDE.md` in the shared PIC workspace

Status meanings:

- **ACTIVE**: work on the immediate path.
- **NEXT**: small, concrete cleanup to do after or alongside active work.
- **LIMITATION**: known and documented; not on the critical path for current runs.
- **BACKLOG**: act only when profiling, a failure, or a planned configuration triggers it.
- **CLOSED**: fixed and covered by focused evidence.

Closed work is summarized here rather than carried as hundreds of lines of completed
checklists:

- `TURB-1`–`TURB-4`: **CLOSED**. Correct normalization root, RK-consistent work,
  active-dimension projection, and signed isotropic modes. Implementation `a458bf08e`;
  focused validation `75da80082`.
- `MIG-1`: **CLOSED**. Unified migration/destruction compaction with survivor, restart,
  and ledger tests. Implementation `a458bf08e`; focused validation `75da80082`.
- `STAB-1`: **CLOSED**. Post-feedback admissibility check before EOS repair in
  `9364e0ff3`.
- `DT-1`, `DT-2`: **CLOSED**. Global configured-species gyro bound and pre-injection
  timestep bound in `9364e0ff3`.

These fixes should remain covered by their focused regressions. They do not need to be
requalified through a new platform/precision matrix for every subsequent change.

## 1. Active now

### `ROBUST-4` — Use the correct MPI datatype for `Real`

- Status: **ACTIVE — quick fix**
- Scope: `src/srcterms/turb_driver.cpp`

Three turbulence reductions still pass `Real` buffers to `MPI_DOUBLE`. This is correct
only in double-precision builds and can overrun buffers in single precision.

Do:

1. Replace those datatypes with `MPI_ATHENA_REAL`.
2. Run the focused turbulence regression in the normal build.
3. Run one MPI single-precision turbulence smoke to exercise the repaired path.

Done when the focused tests pass and no hard-coded `MPI_DOUBLE` remains for `Real` buffers
in the turbulence driver.

### `HALL-1` — Implement the complete large-scale CR-Hall closure

- Status: **ACTIVE — highest physics priority**
- Scope: particle moments/pusher/tasks, MHD CT and energy update, model documentation

The current production coupling uses ideal-MHD induction. The existing
`current_to_ct_experimental` mode is a source-isolation experiment with a free
coefficient; it is not the physical CR-Hall model.

Do:

1. Write the compact paper-to-code map for signed CR charge, current, background-ion
   charge density, electric field, force, time level, and storage location.
2. Derive the Hall EMF from the deposited `J_cr - q_cr u_g` and the physical electron
   charge denominator. Do not use a freely tunable Hall-strength coefficient.
3. Use the same midpoint full electric field in the particle push and constrained
   transport.
4. Keep gas momentum exchange, particle energy exchange, and the Hall-related gas energy
   flux mutually consistent and counted exactly once.
5. Provide one atomic physical choice, `pic_cr_hall_mode=full|off`. New coupled science
   decks use `full`; `off` is retained for legacy reproduction and controlled comparison.
6. Add only `max(|R|)` and `max(Lambda)` as new history diagnostics.
7. Update the model contract and nearby comments where signs, units, or centering are not
   obvious from the code.

Focused validation:

- one algebra/zero-limit test;
- one manufactured uniform Hall-EMF test;
- one existing periodic exchange test with Hall enabled; and
- the Bell comparison in `BELL-1` below.

Done when these focused tests pass on CPU/MPI and one Frontier GPU smoke exercises the
full path. Do not build a multidimensional Hall parameter matrix before using it.

### `BELL-1` — Compact physical Bell qualification

- Status: **ACTIVE after the first `HALL-1` implementation**
- Scope: corrected volume-aware current setup and linear Bell problem

Do:

1. Verify deposited `J_cr/c` directly from raw particle/grid output for one corrected
   volume-aware case.
2. Run one Hall-off linear eigenmode and recover growth rate, dominant wavenumber, and
   polarization/helicity.
3. Run one Hall-on case near `Lambda ~ 1` and recover the Hall-shifted growth,
   wavenumber, real frequency, and polarization.
4. Repeat the Hall-on case at one higher useful resolution.
5. Vary particles per cell only if the first mode fit is visibly noise-limited.

Done when the two physical branches agree with their dispersion relations and the single
resolution repeat supports the measured quantities. Proceed directly to a modest nonlinear
pilot; add another sensitivity only when that pilot identifies one.

### `SHOCK-1` — Choose one reproducible shock and injection workflow

- Status: **ACTIVE after the compact Bell test**
- Scope: `src/pgen/tests/pic_parallel_shock.cpp` and the active shock deck

The existing immediate-injection/early-removal history and delayed-injection history are
not physically equivalent because removing particles does not restore what was subtracted
from the gas. Stop treating both as interchangeable science setups.

Do:

1. Select and document one startup history; delayed injection after initial shock
   formation is the current preferred simple choice.
2. Track the area-averaged shock position with one smoothed density- or pressure-gradient
   criterion.
3. Maintain one swept-mass and injection ledger. Add local surface reconstruction only if
   later shock corrugation makes the area-averaged tracker inadequate.
4. Run a no-CR planar shock, a low-efficiency full-Hall shock, and one matched Hall-off
   comparison.
5. Choose one resolution or particle-count repeat from what the pilot shows is limiting.

Done when shock speed/compression are sensible, the injection and gas-particle exchange
ledgers close to the expected numerical accuracy, and upstream current and magnetic growth
support a quantitative Bell interpretation.

## 2. Small correctness cleanups

These are narrow fixes, not new validation programs.

### `CLEAN-1` — Stop silently rewriting invalid particle weights

- Status: **NEXT**

Several push/deposit paths replace a nonpositive particle weight with one. That silently
changes mass, charge, current, and feedback.

Validate weights once at the narrowest common construction/load boundary and fail clearly
if an invalid weight reaches a supported full-f workflow. Remove the fallback assignments.
Add one invalid-weight regression; do not add a per-cycle full-particle audit.

### `CLEAN-2` — Prevent physical use of the legacy Bell normalization

- Status: **NEXT**

The old Q023 relation omits root-cell volume and includes the artificial light speed in
the current target. Active science decks must use the corrected volume-aware relation.

Keep the old generator only if it is explicitly labeled and guarded as
`legacy_nonphysical_mechanics_only`; otherwise remove it and its active decks. One startup
guard test is enough.

### `CLEAN-3` — Repair default regression discovery

- Status: **NEXT**

`pic_paper_smooth_tsc_oracle.py` is a helper but is discovered as a regression despite
lacking `run()` and `analyze()`. Rename/reroute it as a helper or add the normal
interface, then confirm the focused particle suite reaches its real tests.

## 3. Known limitation

### `AMR-1` — Coupled refinement-interface deposition is not conservative

- Status: **LIMITATION — uniform grids are the supported science path**

Receiver-resolution TSC is evaluated independently on each refinement level without a
cross-interface partition of unity. Integrated gas feedback can therefore differ from the
opposite particle change. Current exact-conservation shock paths correctly reject AMR/SMR.

Policy now:

- Do not claim exact gas-particle conservation for coupled AMR/SMR runs.
- Do not make AMR qualification a gate for uniform-grid Bell, shock, or turbulence work.
- Keep the limitation visible in the model/deck documentation.

If a planned production run needs AMR, resume with the smallest useful sequence:

1. one planar static-refinement crossing;
2. particle count plus momentum/energy and first-moment closure;
3. one rank-split repeat; and
4. only the additional corner, dynamic-refinement, restart, or GPU case used by that run.

## 4. Profile- or failure-triggered backlog

These observations are real, but they are not prerequisites for the current uniform-grid
science program.

### Performance

- AMR deposition performs host mirroring, global communication, and repeated allocation.
  Address this only if `AMR-1` becomes active.
- Ordinary migration uses global metadata and allocation churn. Profile a
  production-shaped run before redesigning it.
- Shock injection is globally replicated and repeatedly reallocates particle arrays.
  Optimize it if the shock pilot shows that injection materially limits runtime.
- `pic_sort_interval` is parsed but not implemented. Do not implement a sorting
  framework until profiling justifies it; remove or reject nonzero settings in active
  decks so inputs do not claim nonexistent behavior.

### Robustness

- Some rank-local error exits could strand MPI peers. Replace a specific path when it is
  encountered or touched; do not begin a whole-code fatal-path audit.
- `TaskStatus::fail` is not propagated by the generic task-list driver and could
  busy-loop. Track this as a general AthenaK issue and fix it if a reachable task
  failure is observed.
- Do not add a comprehensive per-cycle particle validator unless an actual corruption
  demonstrates the need.

## Deletion rule

Delete this tracker when the active items and small cleanups are closed or transferred to
the ordinary project backlog, and the AMR limitation is documented in the permanent model
contract. Deletion does not require completing speculative optimization, every platform
combination, or the entire future AMR program.
