# MHD-PIC Active Development Tracker

> Temporary tracker for `PIC_development`, refreshed 2026-07-11. The committed base is
> `e37799eb85`; the CR-Hall implementation was developed from that base on
> `PIC_development`. Detailed findings from the original review remain in Git history at
> `9364e0ff3:MHD_PIC_REVIEW_TRACKER.md`.

This file tracks only work that can affect the next uniform-grid MHD-PIC science runs. It
follows [ethos.md](ethos.md): physical consistency is mandatory, while extra validation,
hardening, and optimization require a demonstrated reason. Validation artifacts are under
`/lustre/orion/ast207/proj-shared/dfielding/PIC/validation`.

Status meanings: **CLOSED** has decisive focused evidence; **ACTIVE** is on the immediate
execution path; **NEXT** is a small concrete cleanup; **LIMITATION** is documented and
outside the supported current workflow; **BACKLOG** requires a profile, observed failure,
or planned use case.

## 1. Active now

### Completed foundation

- **BASELINE / ROBUST-4 — CLOSED.** The focused normal-precision
  `particles/pic_turbulent_dynamo_smoke` passed in job `4972237`. The repaired
  single-precision build and two-rank turbulence smoke passed in job `4972471`.
  `turb_driver.cpp` now reduces `Real` buffers with `MPI_ATHENA_REAL`; no hard-coded
  `MPI_DOUBLE` remains in that driver.
- **HALL-1 — CLOSED for the supported uniform-grid model.**
  `pic_cr_hall_mode=full|off` is an atomic
  choice across the midpoint particle field, CT edge EMF, gas-particle momentum and energy
  exchange, and Hall gas-energy flux. The full path uses deposited
  `K_cr - Q_cr u_g` and the physical denominator
  `Q_e = alpha_i rho + Q_cr`; the artificial CR light speed and the legacy experimental
  coefficient do not set the Hall strength. The only new history quantities are maximum
  `|R|` and maximum `Lambda`. Signs, normalization, staggering, and time centering are
  recorded in [MHD_PIC_CR_HALL_CODE_MAP.md](MHD_PIC_CR_HALL_CODE_MAP.md) and the permanent
  model contract.
- **CONS-1 — CLOSED.** The original centered-route run `4972459` passed, and the
  limited-face requalification job `4972934` produced a terminal `PASS`.
  Hall-off results are invariant when the legacy experimental coefficient changes from
  zero to seven. The full mode produces a distinct particle response while conserving
  total momentum to `9.96e-7` and total energy to `9.67e-6` in the focused test.
- **BELL-LIN-1 — CLOSED.** Immutable-source job `4973137` passed all 12 raw-current,
  Hall-parameter, growth, frequency, wavelength, polarization, and resolution checks for
  the final limited-face route. The Hall-off growth and frequency errors are `1.16%` and
  `0.069%`. At `R=0.01` and `Lambda=1`, the full-Hall growth and real-frequency errors are
  `0.70%` and `0.41%`; doubling the parallel resolution changes them by only `0.54%` and
  `0.13%`. The measured deposited values are `J_cr/c = (4*pi,0,0)` and
  `Q_cr/c = 0.04*pi`, independently reduced from raw particle moment outputs. The first
  face-route job `4972780` had exposed a late high-wavenumber branch because it selected
  raw donor cells before interface reconstruction. Applying the Appendix-B order—limited
  PLM interface reconstruction followed by density-flux upwind selection—removed that
  resolution failure.

Earlier turbulence, migration, admissibility, and timestep blockers (`TURB-1`--`TURB-4`,
`MIG-1`, `STAB-1`, `DT-1`, and `DT-2`) remain closed by their focused regressions. They do
not need a new validation matrix for each Hall run.

### `BELL-NL-1` — Modest nonlinear Bell pilot

- Status: **ACTIVE — prepared, qualification pending**

Three superseded centered-route pilots lost gas admissibility at nearly the same physical
time: the base run, a four-times-higher particle-count run, and a four-times-smaller-CFL
run. This ruled out particle noise and timestep size and triggered the matched face
induction/energy route. Run one clean corrected-route pilot from `t=0`; confirm terminal
nonlinear growth, finite states, and sensible `R` and `Lambda`. Add no further variation
unless that report exposes a new sensitivity.

### `SHOCK-1` — Controlled delayed-injection shocks

- Status: **ACTIVE after `BELL-NL-1` — deck and analysis prepared, runs pending**

Use delayed CR injection after the shock forms, one area-averaged shock tracker, and one
swept-mass/injection ledger. Run the no-CR, full-Hall, and matched Hall-off cases. Check
shock speed and compression, ledger closure, upstream current and predicted Bell scale,
magnetic growth, and CR acceleration. Choose one resolution or particle-count repeat from
the pilot evidence; do not pre-build a parameter matrix.

### `BOX-128-1` — Return to the turbulent box

- Status: **ACTIVE after the shock comparison — prepared, run pending**

Run one canonical `128^3` turbulent box with the full closure and analyze it as a
self-contained worker job. Add a full-size Hall-off comparison only if the measured
`Lambda` makes the comparison scientifically useful.

## 2. Small correctness cleanups

These are narrow fixes, not new validation programs.

### `CLEAN-1` — Reject invalid particle weights

- Status: **NEXT**

Validate weights once at the common construction/load boundary for supported full-f
workflows and remove fallbacks that silently replace a nonpositive weight with one. Add
one invalid-weight regression, not a per-cycle particle audit.

### `CLEAN-2` — Quarantine the legacy Bell normalization

- Status: **NEXT**

Active science decks must use the volume-aware current relation. Keep the old Q023 setup
only if it is explicitly labeled and guarded as nonphysical legacy mechanics; otherwise
remove it and its active decks. One startup guard test is sufficient.

### `CLEAN-3` — Repair focused test discovery

- Status: **NEXT**

Stop discovering `pic_paper_smooth_tsc_oracle.py` as a standalone regression without
`run()` and `analyze()`. Treat it as a helper or give it the normal interface, then confirm
the focused particle suite reaches its actual tests.

## 3. Known limitations

### `AMR-1` — Coupled refinement-interface deposition is not conservative

- Status: **LIMITATION — uniform grids are the supported science path**

Receiver-resolution TSC is evaluated independently on each refinement level without a
cross-interface partition of unity, so coupled AMR/SMR runs cannot claim exact integrated
gas-particle exchange. The current exact-conservation shock path correctly rejects this
configuration.

Do not make AMR qualification a gate for uniform-grid Bell, shock, or turbulence work.
If a planned science run needs AMR, begin with one planar static-refinement crossing and
measure particle count, momentum/energy exchange, and first-moment closure. Add only the
rank split, restart, corner, dynamic-refinement, or GPU case that the intended run uses.

The target full-Hall claim is deliberately limited to uniform Cartesian, full-f,
non-relativistic ideal MHD with the VL2/TSC coupling path. Other model combinations remain
unsupported rather than silently approximated.

## 4. Profile- or failure-triggered backlog

- **Performance — BACKLOG.** Profile a production-shaped nonlinear Bell or shock run
  before rewriting deposition, migration, or injection. Address host mirroring, global
  metadata/communication, allocation churn, or replicated shock injection only when the
  profile shows material cost. Do not implement particle sorting merely because
  `pic_sort_interval` exists; reject nonzero settings until a measured need justifies it.
- **Robustness — BACKLOG.** Fix a rank-local fatal path when it is encountered or touched,
  rather than auditing the entire code. Revisit generic `TaskStatus::fail` propagation if
  a reachable task failure demonstrates a hang. Add broader particle validation only in
  response to observed corruption.
- **Sensitivity — BACKLOG.** Particle-count, timestep, decomposition, precision, and broad
  platform sweeps are triggered by a failed fit, noisy pilot, or production requirement;
  they are not standing gates.

**Deletion rule.** Delete this tracker when the nonlinear Bell, controlled shock, and
`128^3` full-Hall steps are complete and the three small cleanups are closed or transferred
to an ordinary backlog. AMR and speculative optimization do not block deletion once their
limitations are permanent documentation.
