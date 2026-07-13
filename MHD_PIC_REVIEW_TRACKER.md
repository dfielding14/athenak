# MHD-PIC Active Development Tracker

> Temporary tracker for `PIC_development`, refreshed 2026-07-12. The validated candidate
> is being consolidated above `e15c820f4`: `978344bcd` contains the conservative
> CR-Hall update, `b808d3400` the narrow correctness cleanups, and `683e120ca` the
> analysis and shock workflow. Detailed review history is in Git; this file tracks only
> work relevant to the next uniform-grid science runs.

Follow [ethos.md](ethos.md): physical consistency is mandatory, but validation,
hardening, and optimization must answer a concrete risk. Validation artifacts live under
`/lustre/orion/ast207/proj-shared/dfielding/PIC/validation`.

Status meanings: **CLOSED** has focused evidence; **ACTIVE** is on the immediate path;
**NEXT** is a narrow cleanup; **LIMITATION** is outside the supported workflow; and
**BACKLOG** requires a profile, observed failure, or planned use case.

## 1. Active now

### Candidate CR-Hall closure requalification

- Status: **CLOSED for the supported uniform-grid model**
- The candidate consistently routes the exact deposited stage-2 particle impulse through
  the gas corrector and places the CR-Hall induction and energy fluxes in the conservative
  face-flux/FOFC update. Source-aware FOFC tests the same PIC feedback applied to the live
  state.
- Jobs `4972934` (focused conservation) and `4973137` (12-check linear Bell suite) remain
  useful historical evidence for the previous route; they do **not** qualify the current
  discrete candidate.
- Candidate linear job `4975948` passed all 12 unchanged Hall-off, full-Hall,
  raw-moment, polarization, wavelength, and resolution checks. Full-Hall growth and
  frequency were `5.5254` and `2.5701`; the doubled-resolution values were `5.5502`
  and `2.5787`.
- Exact-final-source core job `4978017` passed the normal-precision turbulence and
  focused Hall-conservation regressions, all seven safety guards, the static task-order
  oracle, and the Q019 source-closure validator. The full-Hall total momentum and energy
  residuals remained `8.33e-7` and `5.18e-6`. Single-precision two-rank GPU smoke
  `4975966` also passed.

### `BELL-NL-1` — Nonlinear Bell onset pilot

- Status: **CLOSED for onset only**
- The Hall-off source probe `4975816` lost gas admissibility at `t=1.41828`, directly
  after the exact PIC source update. The source-aware full-Hall run `4975854` failed near
  `t=1.5255`; halving the CFL in `4975875` reproduced the same trajectory and failed near
  `t=1.5251`. A `160x64x64` repeat (`4975930`) failed earlier, near `t=1.364`, in the same
  20 coherent tiled filament cores; its peak raw-output
  `B_perp,rms/B0` was `0.742` at `t=1.35`.
- The bounded onset job `4975940` passed all six unchanged checks at `t=1.5`, reached peak
  `B_perp,rms/B0=1.043`, and recovered the exact initial `R`, `Lambda`, current, and charge.
  This qualifies nonlinear onset, **not** saturation or convergence.
- The evidence does not support adding a timestep limiter or more general robustness
  machinery now. Treat deep filament-core thermal collapse as the present coarse-pilot
  boundary and move on after focused candidate requalification.

### `SHOCK-1` — Controlled delayed-injection shocks

- Status: **CLOSED as an engineering pilot**
- Job `4976035` completed the no-CR, Hall-off, and full-Hall trio through `t=100`; the
  compact completeness gate passed. Measured shock speeds were `10.06`, `9.98`, and
  `9.98`, with terminal compression ratios `3.99`, `4.01`, and `4.01`.
- The two CR ledgers agree exactly: injected mass `422.3286` plus the `4.89e-4`
  reservoir gives swept mass `422329.09` at `eta=1e-3`. The terminal full-Hall upstream
  drive current is `0.06520`, predicting `lambda_Bell=192.93` code units, or about
  `64` cells at `dx=3`.
- This first shock is deliberately an engineering proxy, not an acceleration or
  nonlinear-Bell result. Its terminal upstream `Lambda_parallel=0.0651` and peak
  domain-wide `Lambda=0.150` make Hall shifts small; `p_max/p_inj=1.0091` and no particle
  weight exceeds `2 p_inj`. The current is smooth and the Bell scale is already well
  resolved, so neither a resolution nor particle-count repeat is triggered. After the
  turbulent-box comparison, run one longer/larger full-Hall successor at unchanged
  resolution and injection physics; do not launch a speculative matrix.

### `SHOCK-2` — Long full-Hall shock

- Status: **CLOSED — coherent 2D Bell-in-a-shock result**
- This one `5000 by 128`, `t=400` full-Hall run kept the pilot's `dx=3`, injection
  physics, particle weight, seed, solver, and coupling parameters. It added upstream
  extent and duration rather than a control matrix.
- The first submission, `4977517`, stopped during CMake after 14 seconds because the
  archive retained a stale `kokkos/.git` worktree pointer; no simulation ran. The
  portable archive compiled cleanly in debug build proof `4977596`.
- Job `4977605` then reached `t=165.3` before stochastic carrier-local gas subtraction
  would have pushed one cell below the pressure floor. It stopped before clipping; the
  healthy partial shock isolated source localization, rather than resolution, particle
  noise, or Hall strength, as the demonstrated limiter.
- Job `4977996` reran the same science case with the existing exact surface-averaged gas
  sink and a 31-cell downstream stencil. It reached `t=400` with all simulation,
  completeness, injection-ledger, and shock gates passing. The shock remained stable at
  `v_sh=9.941` and compression `4.009`.
- Hash-bound post-analysis job `4978079` measured `N_Bell=6.64`, terminal upstream
  `B_perp,rms/B0=0.185`, and a 488.8-code-unit (`162.9`-cell) dominant transverse mode
  with helicity `-0.998`. The near-shock transverse field reaches `0.673 B0`, and the
  morphology shows the expected magnetized filaments and underdense cavities. Together
  the sustained current, amplification, scale selection, circular polarization, and
  morphology establish a coherent **2D** Bell instability in the shock precursor.
- The terminal mean/local-maximum `Lambda` values are only `0.059/0.203`, with a mean
  Hall wavelength shift below `0.1%`; a matched Hall-off shock is not warranted. The
  shortest local Bell scale still has `20.6` cells and the observed mode has `162.9`, so
  no resolution repeat is triggered. The smooth current and coherent morphology with
  4.78 million particles do not trigger a particle-count repeat. The run also shows
  emerging acceleration (`pmax/pinj=2.45`, `p99/pinj=1.78`), without making a mature
  spectrum claim.

### `SHOCK-3` — Modest 3D Bell shock

- Status: **ACTIVE — throughput gate before one 3D science run**
- First run one production-shaped 3D probe only through `t=47`, just past delayed
  injection. The exact surface-averaged source preparation currently repeats global host
  work on every rank, so this probe measures the real scaling risk before committing the
  larger allocation.
- If throughput is acceptable, run one full-Hall 3D case through `t=350` to test whether
  the polarization, magnetic amplification, cavities, filaments, and CR transport survive
  without the 2D geometry. If throughput is not acceptable, profile and fix only the
  measured source-preparation bottleneck; do not launch a parameter matrix.

### `BOX-128-1` — Return to the turbulent box

- Status: **CLOSED — matched full-Hall/Hall-off comparison passed**
- Full-Hall job `4976062` reached `t=10` and passed all five gates. Maximum driven-energy
  residual was `1.92e-4`, total-momentum norm `2.09e-6`, and `divB=4.98e-13`.
  It amplified `B_rms` by about `100`, ending at `0.09998` with `v_rms=0.4710`.
- Every sample in the final 30% had maximum `Lambda >= 0.1`; the tail median, p90, and
  peak were `41.93`, `82.55`, and `277.41`. This decisively triggered one matched
  Hall-off control under the predeclared rule. Debug job `4976612` reached `t=10` and
  passed simulation, analysis, and comparison gates with the identical source archive,
  input, seed, build, and eight-GPU layout, changing only basename and Hall mode.
- Full Hall ends with `B_rms` `9.21%` above Hall-off and magnetic energy `19.22%`
  higher; the final-30%-median magnetic-energy difference is `6.57%`. Both spectra peak
  at mode `9`; velocity RMS differs by only `-2.00%` at the final snapshot and CR energy
  by `+0.042%`. This is a measurable finite-duration Hall response, not a convergence
  result. Add no broader box matrix unless a later science claim requires it.

The earlier baseline jobs (`4972237`, `4972471`) remain useful history; current-candidate
coverage is supplied by `4978017` and `4975966`. These checks need not be repeated for
every subsequent science deck.

## 2. Small correctness cleanups

These are narrow fixes, not new validation programs.

- **`CLEAN-1` — CLOSED: invalid CR macro weights are rejected.** One device reduction
  plus coordinated MPI maximum validates the assembled population after fresh problem
  setup and after restart problem setup. Nonfinite and nonpositive weights fail there;
  push, deposition, and turbulent-history consumers now use the stored weight directly
  instead of silently replacing it with one. The compact safety test covers NaN, zero,
  negative, boundary-call, and no-repair semantics. Post-cleanup worker job `4976135`
  compiled the full tree and passed the normal-precision turbulence smoke plus all seven
  focused safety guards.
- **`CLEAN-2` — CLOSED: legacy Bell normalization quarantined.** Fresh and restart
  dispatch now default-deny the historical Q023/Q029 generators whose normalization omits
  root-cell volume and includes an artificial-`C` factor. Reproducing those mechanics
  requires the explicit `allow_legacy_nonphysical_bell_mechanics=true` opt-in; physical
  runs are directed to Q043 or the corrected Q023 J-over-c path. The focused safety suite
  passes without modifying the hash-bound historical decks.
- **`CLEAN-3` — CLOSED: focused test discovery repaired.** `run_tests.py` now treats
  `*_oracle.py` modules as helpers, so broad discovery skips the standalone
  `pic_paper_smooth_tsc_oracle.py` while retaining the compiled interface regression that
  imports it. Discovery, the four existing oracle tests, and the standalone command pass.

## 3. Known limitations

- **`AMR-1` — LIMITATION: coupled refinement-interface deposition is not conservative.**
  Receiver-resolution TSC is evaluated independently on each refinement level without a
  cross-interface partition of unity. Uniform Cartesian, full-f, non-relativistic ideal
  MHD with VL2/TSC is therefore the supported science path; the exact-conservation shock
  path correctly rejects coupled AMR/SMR.
- AMR does not gate uniform-grid Bell, shock, or turbulence work. If a planned science run
  needs it, start with one planar static-refinement crossing and add only cases required by
  that run.

## 4. Profile- or failure-triggered backlog

- **Performance — MEASURED BACKLOG.** The full-Hall box spent `21.4%` of wall time in
  particle migration; output was below `1%` and particle/mesh balance was already
  excellent. The remaining `77.7%` stage block cannot be split because synchronized
  particle kernel timers were disabled. If runtime becomes limiting, first run one short
  restart with migration substep/byte counters, then optimize only the demonstrated
  communication or polling cost. No output, load-balance, or Hall rewrite is justified.
- **Robustness — BACKLOG.** Fix reachable rank-local fatal paths, task-failure propagation,
  or particle corruption when encountered; do not conduct a broad speculative audit.
- **Sensitivity — BACKLOG.** Particle-count, timestep, decomposition, precision, and broad
  platform sweeps require a failed fit, noisy pilot, or production need. The nonlinear
  evidence above does not justify further timestep machinery now.

**Deletion rule.** Delete this tracker after the modest 3D item is completed or transferred
to the permanent science guide. The 2D Bell, controlled-shock, matched-box, and cleanup
requirements are now closed; AMR and speculative optimization do not block deletion once
their limitations are documented permanently.
