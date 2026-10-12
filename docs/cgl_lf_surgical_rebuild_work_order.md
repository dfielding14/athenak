# CGL LF surgical rebuild work order

Deliver a correct uniform active solver with attributable cost, then qualify mesh
and passive support separately. Follow the
[rebuild plan](cgl_lf_surgical_rebuild_plan.md) for equations, source scope,
test selectors, and GPU recipes; this checklist tracks execution.

Branch: `c/cgl-lf-surgical-rebuild`.
Controls: pre-WO1 `8222de3aa`, post-WO1 `9a4030b09`.
Corrected reference: `2afaa0fa1`.
Here pre-WO1 is this rebuild's starting implementation; post-WO1 is the completed
historical first work order, not the candidate on this branch. The older scaling
archive uses a third revision, `f2d0a259`.

**Status:** pieces 0–14 are implemented and qualified within their recorded
scopes. The frozen active milestone is complete; the passive turbulent milestone
remains open. Piece 15's default-off implementation and focused checks pass;
turbulent interval qualification remains open. Piece 14's active-parent bitwise
comparison failed (details below).
**Workspace:** `/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild`.
**Next:** Set physical-observable and pressure-depletion tolerances, then use
short matched-cadence passive continuations to select an LF interval before
passive Tier 3. Reuse the completed one-cycle studies. Larger refresh/cache and
logarithm changes remain deferred. The
[workflow](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/active-milestone/README.md)
retains reproduction commands and links to the completed milestone evidence.
Reuse accepted physics and reference evidence unless a relevant change or concrete
failure requires repeating it.
Keep the plan's piece numbers; execute forcing (11) immediately after collisions
(4). The uniform active milestone precedes optional optimization and later models.
Passive work (14) does not depend on mesh work (12) or optimizations (13).
Piece 15 addresses the reproduced passive pressure loss and active threshold sensitivity.

**Working rules**

- Before each edit, trace all affected callers, tests, inputs, diagnostics, and
  documentation. Reuse existing helpers and import only required reference code.
- Keep independent fixes in separate commits with their focused regressions.
  Land coupled equations, scheduling, and callers together. Pass the focused
  gate before stacking dependent changes; broaden checks at milestones.
- Preserve active E/A, thermal energy `U = p_perp + p_parallel/2`, temporary LF
  magnetic moment, RKL2 coefficients, odd stages, and the two LF half-sweeps.
  Reuse the existing safe arithmetic and already-correct numerical machinery.
- Use independent physical oracles. Preserve inherited scientific tolerances;
  record known B4 failures and investigate new failures before dependent work.
- Keep one compact completion record below. Store raw logs and full-precision
  comparison states outside Git. Update this status and the next piece as work lands.

**Execution checklist**

- [x] **0 — Controls and baseline.** Verify the pinned Frontier environment and
  retained artifacts. Build both controls once, separately from the candidate,
  with matching precision/Kokkos and `-fno-cray -mno-daz-ftz` for application and
  Kokkos, plus `-no-pie`. Run one minimum-reduction probe, host arithmetic checks,
  and the prescribed GPU decay/3D/layout/restart checks. Record starting failures.
  Freeze the aligned uniform and eight-GPU `192 x 192 x 384` active workloads;
  verify nonzero LF transport and compatible active initial data. A passive J/A
  checkpoint cannot initialize active E/A. Stage frozen HIP tests outside the
  source tree using the plan's layout; avoid generic CUDA/rebuild runners.

- [x] **1 — Guards and parameters.** Reject passive mode pending redesign;
  preserve LF boundary guards and validate existing floors/rates. Gate: invalid
  inputs fail clearly and valid active inputs retain their behavior.

- [x] **2 — EOS repairs.** Fix density/magnetization floors and persist repaired
  state. Preserve pressure ratio in the A decoder; preserve magnetic moment,
  momentum, and energy in the moment decoder. Gate: independent pressure/energy
  identities, idempotence, stored A, and focused double/float checks.

- [x] **3 — Reconstruction and FOFC.** Repair both CGL pressure floors in every
  directional PPMX/WENO-Z wrapper and detect incoming nonfinite E/A before repair.
  Reject nonfinite recovered velocities and repaired A as well as pressures.
  Gate: affine/quadratic wrapper tests, NaN/Inf injection through the live FOFC
  predictor/correction path, steep-pressure evolution, and MHD controls.

- [x] **4 — Collisions and walls.** Implement thresholds, relaxation, additive LF
  scattering, canonical wall encoding, and all scheduling/AMR/input consumers
  coherently. Apply full-dt rates once per cycle and walls at prescribed checkpoints,
  including non-LF parabolic sweeps and the mandatory fluid-firehose wall without
  optional limiters. Preserve small pressures under tiny inward rate kicks. Update
  admissibility: count intermediate wall crossings; reject invalid/floor-repaired
  stages and enforce walls at checkpoints without inflating RHS counts. Gate:
  independent rate/cadence/equality checks, invariants/no-op/idempotence, production
  HIP wall encoding followed by C2P, and real restart. Then separately remove the
  redundant CGL timestep task after verifying unconditional post-wall ownership,
  including zero rates and no LF.

- [x] **11 — Forcing.** Separate mode selection, conserved-momentum work, and
  scheduling fixes; move task registration and required refresh together. Gate:
  RK1/2/3 full-step kick counts, power, unchanged thermal energy, OU cadence/restart,
  and affected hydro/MHD/two-fluid consumers. Preserve seeds and prescribed power.

- [x] **5 — Weak-field faces.** Repair neighboring-field use in HLLE/LLF while
  retaining A. Gate: first-update contamination, both contact orientations, smooth
  interface convergence, and unchanged predicates for the two B4 expected failures.

- [x] **6 — LF geometry.** Use CT normal B and mean adjacent cell magnitudes;
  do not renormalize to unit length. Gate: independent face fluxes, reversals,
  frozen-field conservation, and oblique decay.

- [x] **7 — Transverse slopes.** Apply VL4 to temperature gradients, retaining
  unlimited magnetic-magnitude slopes and consistent diagnostic/profile copies.
  Gate: hotspot minima/energy, oblique convergence, exceptional slopes, and a
  separate paired timing of the added arithmetic.

- [x] **8 — Perpendicular closure.** Correct the `2*nu_eff` response in normal
  and logarithmic paths and reference consumers, explicitly identifying the
  SHD/BGK convention. Replace stale whole-endpoint collision-map references with
  the independent coupled solution and both signed amplitudes. Gate: collisionless
  identity, intermediate/strong-collision limits, subnormal intermediates,
  newly overflowing response denominators near cancellation, and full-field
  resolved coupled decay. Across 6–8 use the common
  conservative reference step and oblique background-coefficient kinematic screen;
  verify actual dt/RHS counts instead of assuming explicit integration fixes dt.
  Repeat Tier 2 afterward to set the work-count and cost baseline for 9 and 10.

- [x] **9 — Final-stencil timestep bound.** Adapt the coupled CT/density/VL4/grad-B
  rows directly, with fresh halo values in existing scratch storage and safe
  exceptional arithmetic. Retain the old safety multiplier. Gate: rerunnable
  independent rows/eigenvalues/RKL amplification, checkerboards, density jumps,
  reversals, heated refresh, MPI agreement, and explicit-reference convergence.
  Retained closure/hotspot/oblique gates and Tier 1/2 measurements pass; the
  developed-state slowdown is recorded below, separately from uniform savings.

- [x] **10 — STS safety.** Separate STS safety from advective/explicit CFL,
  including restart loading and both budget selections. Gate: parameter bounds,
  analytic odd-stage counts, stability, and temporal convergence before accepting
  a larger factor. Record work-count changes separately from evaluation cost.
  Include an ordinary STS consumer and the unforced finite-rate soft-threshold
  crossing case: require full-pressure-field convergence under outer-step
  refinement before accepting a larger safety. Preserve first-order split semantics.

- [x] **Active milestone.** The physical suite is frozen in
  `build/active-milestone-final`: source `73e2769fc`, selected test source
  `08189c1e2`, 33 passes and two known B4 expected failures. The fixed 8, 4, 2, 1
  local-hotspot window (`5647113`) passes the original physical gates and
  convergence checks. Retain the finer ladder's minimum/conservation failures
  (`5646938`) and invalid B=1 control (`5646989`); B=5 conditioning control
  `5647008` passes its diagnostic criteria. Corrected-reference comparison
  `5647118` is accepted as a bounded descriptive comparison, with localized
  near-threshold differences retained and no equal-accuracy speedup claim.
  Tier 1 job `5647205` completed all 14 invocations on two nodes: production
  LF/no-LF rates 135.512/557.546M zone-cycles/s/node (ratio 0.243051), 26 LF
  RHS/cycle, and rate spans below 0.04%. The historical control is separately
  labeled in the report. Cancelled job `5647148` remains partial evidence.
  Tier 3 job `5647225` completed fresh to t=14 at safety 0.3: 47,681 cycles,
  3.475833 node-hours including I/O and 47.465 LF RHS/cycle. LF floor/invalid
  and EOS repair counters stayed zero; the [6,14] energy minus applied-work
  residual is 2.23e-11 relative to injected work. Analysis job `5647226`
  completed all four figure groups, reviewed on 2026-10-12. The scientific
  interpretation remains inconclusive: one heated realization does not establish
  turbulent convergence or active/passive causality. The run used frozen
  pre-cleanup source `73e2769fc`; later cleanup and Piece 12 evidence are separate.

  External workflow consolidation is complete: three existing entrypoints and
  their required helpers are listed in the workflow above. One-off and duplicate
  audit/setup scripts are archived in place, preserving evidence and provenance
  paths. The [milestone report](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/active-milestone/active-milestone-report.md)
  records the completed evidence and limits. Reuse completed physics and reference
  work unless a relevant change or concrete failure requires a rerun.

- [x] **12 — Mesh support.** Separately repair physical corners after prolongation
  and synchronize shared LF face fluxes. Grow both buffers before receives, with
  sufficient 3D full-face capacity and existing request completion preserved.
  Gate: energy/moment conservation, both prolongation modes, 2D/3D explicit/STS,
  CPU/HIP, one/four ranks, actual regrid/restart, repair accounting, and shearing.
  Rebuild moment halos before the first LF RHS and anisotropy halos after split
  checkpoints; conversion and prolongation do not commute. Require fixed-mesh
  STS/explicit pressure convergence without extra diagnostic RHS counts.
  Two-node churn and native shear restart qualification are complete. CPU/HIP
  each pass 25 checks; the final MPI gate passes both checks. All native shear
  fields/ghosts and integer counters match exactly; one q-work sum differs by
  one ulp after restart regrouping. See the completion record and workspace
  [Piece 12 report](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/piece12/README.md).

  Requested follow-up: a monotonic physical-time target level for uniform
  refinement in the turbulence problem, using the existing AMR hook. Qualify
  forcing/restart and conservation through refinement; passive and chunked LF
  support require their own checks. This callback is not implemented yet.

- [x] **13 — Optional optimizations.** Profile the corrected target and select
  the smallest useful removal/fusion/message change. Compact communication follows
  mesh correctness. Gate each retained change on full-precision primary, ghost,
  coarse, scalar, and restart/RNG equivalence plus paired GPU timing. A simple
  deletion may stand within timing noise; larger machinery needs repeatable benefit.

  **Compact exchange (`920090906`):** later strict LF stages on uniform periodic
  meshes send only IEN/IAN; the first stage and all excluded cases retain full
  messages. No new option or buffer. CPU `5656649` and HIP `5656646` each pass
  nine checks, including exact native scalar/conserved/CT/ghost payloads on
  one/four ranks. Tier 1 `5656691` passes all 18 invocations: 136.122→140.333M
  zone-cycles/s/node (+3.09%), 26 RHS/cycle unchanged, no-LF +0.002%.
  Tier 2 qualified replays `5656928` plus timing continuation `5656954` give
  53.350→58.008M zone-cycles/s/node (+8.73%; all three pairs +8.59–8.88%).
  Native fields/ghosts, forcing/RNG, ledgers, time/dt and work are exact;
  both builds retain 50.4 RHS/cycle and pass the separate full-q checks.
  Retain the change: 22 net production lines and one 51-line regression.
  Failed Slurm allocations `5656757`/`5656928` remain recorded separately;
  successful replays were reused without repeating physics checks.
  Defer the larger refresh/cache fusion and
  upwind-log changes; details and artifacts are in the
  [Piece 13 report](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/piece13/README.md).

  **Authorized cleanup (2026-10-10):** dead code `e763f1d10` passes three host and
  seven HIP checks (`5647261`); environment-override/replay removal `1e8c02c37`
  passes 20 HIP checks (`5647272`). Directional consolidation
  `7f5c7a83f` + `a60e9586b` removes 279 net lines and passes 20 HIP checks (`5647286`).
  Tier 1 `5647287`: parent/candidate 138.080→133.536M zone-cycles/s/node
  (−3.29%), both 26 RHS/cycle, no-LF unchanged. Tier 2 `5647288`:
  52.0043→52.0795M (+0.145%, within variation), both 50.4 RHS/cycle;
  native conserved/CT/ghost state, RNG and ledgers are bitwise identical, and
  both full-q qualification checks pass. Retain consolidation for maintenance
  and neutral developed-state cost, explicitly accepting the uniform penalty.
  First-draft Tier 1 `5647277` (−2.50%) is superseded evidence; `5647278` was cancelled.
  LF-only clear removal `718901fa2` removes 42 net lines. Gate `5647316` exposed
  a stale zero-intermediate-crossing assertion, reproduced on its parent (`5647337`);
  `6c36d9fc4` replaces it with the actual endpoint wall check. Gate `5647341`
  passes 16 HIP checks and the shearing field/magnetic/endpoint comparisons, but
  the shearing restart still reports 1520 versus 1522 LF RHS evaluations. Parent
  `5647358` and the pre-consolidation binary (`5647376`) reproduce that exact
  failure; the assertion was retained and Piece 12 subsequently repaired the
  boundary lifecycle. Both restart paths now perform 1522 LF evaluations.
  Accept this documented existing failure for the clear
  comparison; raw failed records remain intact. Clear Tier 1 `5647364` passes all
  18 invocations: 132.362→137.793M zone-cycles/s/node (+4.10%), both 26 RHS/cycle;
  no-LF changes by −0.043%. Tier 2 `5647365` has identical native endpoints,
  RNG/ledgers and 50.4 RHS/cycle; both full-q qualifications pass with nonzero
  transport and zero floor/invalid events. Three timing pairs give
  51.5675→53.0907M zone-cycles/s/node (+2.95%, repeat spans below 0.13%).
  Retain the clear removal. Earlier dependent timings were cancelled without
  running. Artifacts are in workspace
  `analysis/cleanup/`. Frozen Tier 3 binary `5fa3c0df…`, jobs `5647225`/`5647226`,
  shared helpers and analysis dependencies remain unchanged. No workspace cleanup,
  physics changes, exceptional-arithmetic fallback removal or fast-mode removal.

- [x] **14 — Passive model implementation and qualification.** Source `276eee6ce`
  evolves J/A outside LF and U/moment inside LF on native isothermal flow,
  including forcing, initialization, outputs, restart encoding and support guards.
  Model scope and absent shock heating are documented in
  [passive thermodynamics](source/modules/cgl_passive.md). CPU `5657203` and HIP
  `5657202`/`5657220`/`5657266` pass identity, thermal/work references, independent
  wall decoding, all five reconstructors and forced four-rank restart. The
  resolved nonlinear oracle (`3a2a5409b`) measures STS orders 1.999/1.994;
  failed reference/launch evidence remains retained. Initial scope is uniform
  periodic. Workspace `analysis/piece14/README.md` records the full evidence.
  Active control timings (`5657269`/`5657270`) show Tier 1 LF throughput
  -0.86%, no-LF +8.32%, and Tier 2 +0.14%, with unchanged LF work counts.
  **Tier 2's exact parent-state gate failed:** all applications and health/q-mode
  checks pass, but local pressure differences reach 0.38% and mirror counts
  differ. These are descriptive costs, not an equivalent-optimization result;
  small integrated budgets do not resolve the local accuracy question. Replay
  `5657582` confirms identical loading, then pressure differences of order 1e-13
  after one cycle and 1e-12 after four. Piece 15 identifies the onset's local
  amplification; the earliest arithmetic source remains unresolved.
  **The passive milestone is still open:** stress replay `5657268` hits one
  pressure floor at post-LF stage 25/27 near cycle 10312, t=4.396717, safety 0.3.
  This old state is a stress probe, not a current-model physical result. Tier 3
  is prepared but unlaunched; qualify smaller LF intervals in Piece 15 first.

- [ ] **15 — LF chunking implemented; turbulent qualification open.**
  Source `6ea48ac69` adds about 100 production lines, default off, for uniform
  periodic LF-only RKL2; AMR is excluded. CPU `5658077` and HIP `5658068` pass
  16 checks each, including smooth orders 1.95/2.01, cadence and restart; the
  default-off active one-cycle native state matches Piece 14 exactly. Active
  fixed-outer study `5658095` preserves forcing/RNG and records no repairs, but
  the finest two caps still differ by 0.2341% in parallel pressure locally, with
  nonmonotone worst-cell refinement. Work grows from 48 to 25,626 RHS per cycle.
  No turbulent cap is qualified. The face probe confirms a local abrupt
  scattering/transport branch while retaining its failed full replay gate;
  cross-restart `5658231` shows the onset follows the accumulated starting-state
  difference, with only roundoff-sized differences between executables at fixed
  state. External SSPRK2 `5658198` also passes 16 focused HIP checks but retains
  nonmonotone local refinement differences; only 9-19 cells exceed 0.1% in either
  pressure in the fine comparisons. No production integrator is added. Select
  physical-observable and depletion tolerances before a production cap.
  Passive `5658164` completes all eight positive caps without new LF floor/invalid
  events; cap zero reproduces the failure. Cap 4 differs from the finest endpoint
  by 4.60% in parallel pressure locally; the finest pair differs by 0.918%, with
  5,450 cells (0.0385%) above 0.1% in either pressure. Integrated q-gradient work
  agrees to about 1e-5 relative. The old deck has no event writer for repairs
  outside LF. Keep passive Tier 3 unlaunched. Measurements and limitations are in the
  [Piece 15 workspace record](/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild/analysis/piece15/README.md).

**Measurement and acceptance**

Use the plan's three frozen benchmark tiers; scaling sweeps and higher-resolution
boxes are separate work.

- **Tier 1, uniform cost control:** two nodes, sixteen 256 cubed blocks,
  frozen `uniform-2n` plus its no-LF companion, 40 cycles. Run the byte-identical
  historical input at baselines and milestones only. Use Tier 1 for cost per LF
  evaluation on any patch that changes recurring full-domain work; report the
  LF/no-LF throughput ratio and absolute rates.
  Freeze/hash the no-LF companion in `inputs/README.md` and check both controls.
- **Tier 2, developed-state restart:** matched box, one node, frozen active
  checkpoint near t=6. Measure between five and 40 completed cycles after restart
  (35 updates), with one extra untimed cycle to exclude finalization. Disable
  writers using existing checkpoint keys and verify no output. Count LF work over
  the same 35 updates in a separate history replay. Produce with canonical writer
  blocks present, and enable per-cycle history by setting its existing `dt` below
  the smallest cycle timestep. Keep both ledgers, LF diagnostic mode and forcing
  state unchanged for timings. Qualify actual transport with a separate untimed
  full-q replay and check its physical/forcing endpoint against a q-off replay.
  Include an event writer to inspect EOS/FOFC repairs. Tier 2 measures work count and cost
  in a nonuniform state; require it for 7, 9, 10, bound/refresh work in 13, and 15.
  Attempt and qualify the checkpoint before piece 7; retain and diagnose any
  failure without blocking the required correction. Repeat Tier 2 after 8.
  Parent/candidate share the same compatible checkpoint; regenerate, relabel and
  rebaseline when restart contents or model changes require it. Passive uses its
  own qualified checkpoint; passive J/A cannot initialize active E/A.
- **Tier 3, matched box to t=14:** fresh initial state on one node, active at the
  active milestone, passive at the passive milestone, canonical history/field/restart
  output.
  It provides the physical comparison and node-hours including I/O to the common
  endpoint. Use explicit STS safety 0.3 for the active milestone; Piece 10's small
  refinement examples and cost screen do not select historical 0.9 for turbulent
  accuracy. Chunking requires piece 15.

Use safe/none/weighted LF with strict checks and explicit environment-override
control. Keep pressure-work and injected-work recording enabled on Tiers 2 and 3;
Tier 1 is the core-cost control without them or forcing. For Tiers 1 and 2, use
warmup and matching output-free intervals: one alternating pair screens a patch;
three establish a baseline or support a speedup. Tier 3 is one run per member
at its milestone, repeated for changed settings or failed gates. Run timing jobs
without competition; profile separately with detail off initially.

Report the tier, zone-cycles/s/node, LF evaluations/cycle, elapsed/physical time,
and amortized whole-step time/RHS. A change in Tier 1 evaluations/cycle is a
controller change; judge it on Tier 2. Final acceptance reports node-hours to
t=14 at demonstrated accuracy. Figures are descriptive; reaching t=6 or t=14
does not prove stationarity or accuracy. Set additional endpoint-error criteria
before selecting performance settings. The frozen-row bound does not prove
nonlinear or AMR stability; hotspot convergence does not certify turbulent
threshold switching; once-per-cycle collisions do not make the full split second
order.

Directional fusion, the constant derivative-envelope controller, and boundary
extensions remain conditional under the plan. Half-sweep merging, A-to-Q,
new limiter models, symmetric collision splitting, and broad float support remain
outside scope. Optional work may be recorded as deferred with its reason.

**Completion record**

Append one row per completed piece or independent sub-fix, including any deferred
candidate. Keep commands/backend and evidence paths sufficient to resume the work.

| Piece | Commit and behavior | Check and command/backend | Result or limitation | Timing/work and evidence |
| --- | --- | --- | --- | --- |
| 0 / fast speed | `16f3762b9`; reuse fixture-only `3d4716225` analytic limits and independent eigenvalues | Existing `test_cgl_fast_speed_cpu.py` selector, CCE20 Serial double/float | 2 passed; production formula unchanged | Workspace `runs/fast-speed-5640622/results.xml` and `fixture.patch` |
| 0 / controls | Exact `8222de3aa` and `9a4030b09`, Kokkos `08ceff92b`, qualified CCE20/HIP flags | `analysis/baseline.sbatch`, job `5640622`; reduction probe and four GPU selectors per control; six Serial helper cases | Both probes match analytic dt; 8 GPU and 6 host cases passed | Workspace `runs/baseline-5640622/`; binary hashes, compile commands, XML results, and logs retained |
| 0 / timings | Frozen uniform 2-node and active 3D 1-node inputs; both controls unchanged | `analysis/timing.sbatch`, job `5640674`; warmup, three alternating pairs, separate detail-off profiles and transport checks | All 22 simulations passed; both active checks have nonzero LF work and zero repair/invalid events. Initial report parser failed on one-based history labels; corrected analysis passed without rerunning simulations | Mean zone-cycles/s/node pre→post: uniform 105.154M→116.708M (42 RHS/cycle); active 70.743M→73.712M (38 RHS/cycle). Startup cost only; full intervals, spreads, elapsed/physical-time and inferred whole-step s/RHS in workspace `README.md` and `runs/timing-5640674/summary.json` |
| 0 / archive context | Historical scaling `f2d0a259`, job `4925508`, compared with fresh controls | Recomputed all three from timestamps 5→20 using the archived mean-per-cycle estimator | Archive 103.496M; fresh pre-WO1 105.151M; fresh post-WO1 116.712M zone-cycles/s/node. Runtime/compiler settings differ; no attribution of that gap to a single change | Workspace `analysis/scaling-comparison.json`; original `scaling/current_consistent_primary_f2d0a259` CSV, manifest, input, and log |
| Logging | `82665ed37`; append interval `zone-cycles/s=` to periodic diagnostics using completed global block updates | `analysis/candidate-logging.sbatch`, job `5640730`; qualified HIP build/probe; two ranks, `ndiag=1/3`, cycle-5 restart | Probe and all 3 checks passed; rate covers all ranks and includes interval I/O/regridding. AMR accounting reviewed, evolving-AMR run not repeated | Workspace `build/candidate` and `runs/logging-5640730/`; source patch, commands, logs, and results retained |
| 1 / startup domain | `ad8a420d1`; reject passive evolution; validate finite positive density/pressure/B floors, magnetization ceiling and LF scales, and finite nonnegative T/entropy floors and rates | Qualified candidate HIP build; `analysis/step1.sbatch` and `step1-startup-retry.sbatch`, jobs `5640759`/`5640764` | 54 startup + 12 active/boundary GPU checks and 1 Serial arithmetic check passed. Three passive evolution tests explicitly deferred; active/scalar coverage retained. Initial new fixtures failed on absent command-line keys, then passed after writing keys into the input; same executable throughout | Workspace `runs/step1-5640759/` and `runs/step1-5640764/`. Both immutable control hashes verified. No timing repeat for constructor-only guards; float representability and passive AMR remain unqualified |
| 2 / custom-pgen build | `f7e1ce884`; enable the turbulence custom entry point only for `PROBLEM=turb` | Full custom-pgen HIP link and zero-step execution, job `5640869` | Existing duplicate-symbol defect exposed by the new production storage fixture; only the entry-point guard imported, no passive machinery | Workspace `runs/piece2-5640866/` preserves failed link; corrected build in `runs/piece2-5640869/` |
| 2 / density and A storage | `5aed2e006`; decode the pressure ratio at original positive density, re-encode/store A at repaired density; preserve moment density-repair contract and nonpositive-density fallback | Existing `test_cgl_c2p_pressure_floor_cpu.py` in Serial double/float; custom `PROBLEM=unit_tests/cgl_c2p_storage_test` HIP binary with `inputs/unit_tests/cgl_transform_check.athinput` | 2 host cases passed, including 15 new density states; all 3 production wrappers preserve analytic invariants and exact second-call state. Initial fixture include order corrected before successful HIP check | `analysis/piece2.sbatch` and `piece2-density-retry.sbatch`; jobs `5640840`/`5640869`; original control hashes unchanged |
| 2 / magnetization ceiling | `798336118`; impose `max(dfloor, B²/sigma_max)` in both scalar decoders, retaining their distinct density-repair invariants | `analysis/piece2.sbatch`, job `5640878`; Serial double/float and all 3 HIP storage wrappers | 2 host cases passed; density- and sigma-triggered device states match analytic pressures, stored fields, counters and exact repeat conversion | Workspace `runs/piece2-5640878/`; original control hashes unchanged; combined piece-2 timing follows the rounding repair |
| 2 / representable repairs | `ddb912c55`; return pressures decoded from repaired storage; round E upward only when needed, and round moment upward if its decoded perpendicular pressure misses the floor | `analysis/piece2-final.sbatch`, job `5640915`; 2 Serial precision configurations and 21 HIP cases through all 3 production wrappers; specialized `PROBLEM=turb` object also checked | Passed density/sigma/pressure repairs, weak fields, lost thermal energy, moment-division rounding and exact repeated states. Retained the existing ±1000-log decoder behavior and all inherited numerical tolerances; removed the float repeated-floor allowance | Workspace `runs/piece2-5640915/`. New unresolved-energy assertion in `5640890` corrected to compare against incoming E: rounded K+M can already leave positive thermal energy. Focused float checks do not qualify whole-application float evolution |
| 2 / integration | Same `ddb912c55` candidate; no rebuild after focused gates | `analysis/piece2-regression-retry.sbatch`, job `5640969`; qualified compiler probe plus frozen decay, 3D runtime, layout and restart HIP selectors | Probe ratio 1; 4 GPU tests passed. Job `5640915` had stopped before these checks because a copied probe helper lacked its companion source; corrected staging only | Workspace `runs/piece2-5640969/`; candidate SHA-256 `acc48e3cf58c78e64bee9075ba858e0870d27928571cae4f32658d25e63993d1`, identical to the completed build in `5640915`; original control hashes unchanged |
| 2 / timing screen | `ddb912c55` candidate versus saved piece-1 binary `8b9e48ec…`; same frozen physics and qualified flags | `analysis/piece2-timing.sbatch`, job `5640972`; one pair per workload after warmup, separate detail-off work-count profiles and active transport checks | All 14 invocations passed. Matched printed timesteps and RHS/cycle (uniform 42, active 38); identical recorded active transport diagnostics, with zero floor/invalid/hard-bound events | Million zone-cycles/s/node before→after: uniform 104.944→104.967 (+0.02%); active 70.892→70.520 (−0.53%). Single-pair cost screen only; no speedup/accuracy claim. Full elapsed/physical time and whole-step s/RHS in workspace `runs/piece2-timing-5640972/summary.json` and `README.md` |
| 2 / ratio underflow | `ef42d98e9`; split an underflowing pressure-ratio exponential into two factors so representable small pressures survive | `analysis/piece2-sweep.sbatch`, jobs `5640984`/`5640989`; bounded grid plus seeded scalar sweep in both precisions. Permanent boundary fixture and 4 HIP extreme-ratio cases in job `5640993` | Before: two analytic underflow failures per precision; after: all 106,252 states per precision pass pressure/energy, primary-field and exact-repeat checks. Added 16 boundary states per precision without relaxing inherited tolerances | Workspace `runs/piece2-sweep-5640984/` preserves failures; `runs/piece2-sweep-5640989/` preserves successful sweep, source, seed and commands. Ordinary arithmetic unchanged |
| 2 / scalar persistence and final gate | `c8f31d67e`; refresh scalar fractions in the lightweight moment decoder when density is repaired | `analysis/piece2-followup.sbatch`, job `5640993`; 2 Serial precision configurations, 24 HIP wrapper cases, 8 read-only FOFC floor checks, 4 device ratio cases, and 4 retained LF selectors | All passed, including scalar mass/fractions, low-B writeback and exact repeated arrays. Independent whole-piece review found no further persistence defect. Original controls and saved initial piece-2 binary unchanged | Workspace `runs/piece2-followup-5640993/`; candidate SHA-256 `4a0881e6522c94941847687859daa35f0c9576ea9b85ac05b75726475a3405d5`. Whole-application float, nonfinite FOFC detection and later physics remain unqualified |
| 2 / follow-up cost screen | `c8f31d67e` candidate versus saved initial piece-2 binary `acc48e3c…` | `analysis/piece2-followup-timing.sbatch`, job `5640998`; same frozen inputs, one measured pair after warmup, separate profiles and transport checks | All 14 invocations passed; matching printed timesteps, RHS/cycle and recorded active transport, with zero repair/invalid/hard-bound events | Million zone-cycles/s/node before→after: uniform 103.932→103.873 (−0.06%, 42 RHS/cycle); active 70.870→69.743 (−1.59%, 38 RHS/cycle). Single pair, no repeatable performance claim. Full metrics in workspace `runs/piece2-followup-timing-5640998/summary.json` and `README.md` |
| 3 / incoming FOFC state | `51b24a3b1`; detect nonfinite rho/momentum/E/A/B before local C2P can hide it, and reject nonfinite or below-floor recovered pressures in the test-only path | `analysis/piece3.sbatch`, job `5641014`; extended production storage HIP fixture | All 59 states passed: NaN/±Inf in all nine incoming fields at finite/zero B, valid controls including extreme ratios, and finite-input pressure overflow. Exact input/primitive/B preservation and counters checked; earlier piece-2 device cases also passed | Workspace `runs/piece3-5641014/device.log`; ordinary C2P repair policy unchanged. Saved piece-2 and original control hashes verified |
| 3 / reconstruction and evolution | `51d70fc86`; all directional PPMX/WENO-Z wrappers use pfloor for both CGL pressures; retain PPM4 and ordinary-MHD behavior | Same job; `test_cgl_reconstruction_floor_cpu.py` double/float plus custom `cgl_reconstruction_floor_test` HIP pgen; pressure-peak, live FOFC and RJ2a selectors | 72 wrapper calls / 1,008 analytic face-variable checks passed per backend/precision. Both 50-cycle peaks stayed finite/positive with zero EOS/FOFC events; live flux/EMF/traction mutation passed; ordinary-MHD PPMX/WENO-Z shock-tube convergence passed at 128/256 cells | 2 host and 5 application pytest cases passed, plus both custom HIP fixtures. Workspace `runs/piece3-5641014/`; qualified flags checked. No throughput rerun under the plan's local-repair cadence; broad physical/MPI/mesh acceptance and whole-application float remain later gates |
| 3 / recovered-state overflow | `2c9be7d9d`; reject nonfinite repaired A and recovered velocities in the read-only FOFC path | `analysis/piece3-repro.sbatch`, job `5641076`, then `piece3-followup.sbatch`, job `5641081`; production storage HIP fixture | Four missing flags reproduced before the fix: weak-B repaired-A overflow and velocity overflow in all three directions, with finite incoming fields and pressures. Afterward all 63 states pass with exactly 59 flags, unchanged input arrays and unchanged EOS floor counters; prior piece-2 device checks also pass | Workspace `runs/piece3-repro-5641076/` preserves the failures; `runs/piece3-followup-5641081/` preserves the passing run. Normal C2P repair policy unchanged |
| 3 / adversarial final gate | `eae8ebbaa`; add admissible quadratic pressure undershoots and real predictor-to-FOFC E/A injections | Serial double/float in `5641076`; HIP reconstruction and five application pytest cases in `5641081` | All pass: 144 wrapper calls / 2,016 analytic comparisons per backend/precision; 18 directional NaN/±Inf flux injections detect both adjacent cells, correct their faces/EMFs/tractions, preserve primary state and reset flags. Both 50-cycle peaks and ordinary-MHD convergence remain clean; inherited tolerances unchanged | Tested sources verified against commits; original controls, piece-2 control and saved initial piece-3 binary unchanged. Candidate SHA-256 `124ad5270cc1a16c269009f862604f504481c6dbe536e66b3ba89adf6d51aaea`. No performance claim or extension of later physics/MPI/mesh/full-application float qualification |
| 4 / strict LF recovery | `d1c22bb44`; reject nonfinite E/moment before weak-field recovery can hide it | Actual lightweight EOS kernel on HIP: healthy control and NaN/±Inf in E/moment, strict checking enabled/disabled | All 14 executions passed expected acceptance/rejection; retained production-storage/FOFC cases passed | Workspace `runs/piece4-qualification-5641169/strict-refresh.xml`, per-case logs and `storage-device.log` |
| 4 / collisions and walls | `037cbb0d4`; full-cycle relaxation, additive LF scattering, canonical wall encoding/cache refresh, scheduled walls and 28 input migrations | Serial pressure/collision and additive-rate helpers in double/float; AMR projection in double; qualified HIP wall, cadence, finite-rate map, restart, heated MPI and AMR gates | Passed 49 application + 1 static case in `5641164`; 80 startup + 9 application + 1 static case in `5641169`. HIP fixture passes six collision/restart and six LF-conversion states. No inherited scientific tolerance relaxed | Workspace `runs/piece4-host-5641111/`, `piece4-host-retry-5641126/`, `piece4-gpu-retry-5641164/`, `piece4-qualification-5641169/`. Stale pressure-cache defect reproduced in `5641155` before correction. Four controls unchanged; host checks do not qualify whole-application float |
| 4 / checkpoint acceptance | Same physics commit; retain original grad-B oracle and add independent intermediate fluid-wall crossing regression | Qualified HIP, `analysis/piece4-final-gates.sbatch`, job `5641184` | Both passed: original grad-B relative error 0.01013132 ≤ 0.1; wall case records 36 crossings with admissible final pressures, conserved energy and unchanged density/momentum/B. 896 cell-stage samples = 64 cells × 14 LF evaluations | Workspace `runs/piece4-final-gates-5641184/`. Exploratory wall run used the wrong linear-response oracle; dedicated wall test keeps that oracle unchanged. Unrelated paper inventory test already expected 20 decks against 29 at `808ffe1bf`; failure retained in `5641155`, proof in `analysis/piece4-preexisting-paper-inventory.json` |
| 4 / timestep ownership and cost | `5d60ba528`; remove only the redundant ordinary CGL stage task; retain post-wall, post-LF, initialization/regrid and ordinary-MHD owners | `analysis/piece4-timing.sbatch`, job `5641186` (11m54s, exit 0); six paired HIP state comparisons, ten heated/MPI/AMR/wall/restart/MHD regressions, then 21 sequential cost/profile/transport invocations | All passed. Full-precision primary/ghost state, histories and restart U/A/staggered-B payloads match exactly, including zero rates, no LF and disabled sweeps. Over 40 cycles, reductions drop 121→81 with unchanged LF RHS (42 uniform, 38 active per cycle). Active diagnostics have nonzero transport and zero repair/invalid/hard events | Million zone-cycles/s/node piece 3→piece-4 physics→final: uniform 104.114→100.117→100.236; active 71.041→69.953→70.274. One measured comparison after warmup, no repeatable speedup/accuracy claim. Workspace `runs/piece4-timing-5641186/summary.json`, `equivalence/results.json`, `pytest/results.xml` and `committed-source.json`; full metrics in workspace `README.md` |
| 4 / backup-wall coverage | `866056091`; exercise configured asymmetric backup walls and a finite kick through the existing production fixture | Job `5641224`; both HIP wall policies | Both passed independent C2P, canonical cache, invariant primary fields, serialized recovery and repeated-wall checks. This job separately exposed tiny-rate scalar failures; it was not an all-passing qualification | Workspace `runs/piece4-wall-probe-5641224/`; existing inherited tolerances retained |
| 4 / non-LF parabolic checkpoints | `d0f08b844`; apply walls after every CGL STS sweep, with full rates still once post-RK without LF or post-sweep with LF | Job `5641282`; 24 cadence selectors (RK1/2/3 × background/soft × STS/explicit/disabled-STS/viscosity), viscous evolution and five other selected regressions: 30 cases passed | Frozen control left all 64 viscously heated cells outside the fluid wall (minimum margin −0.002884); fixed minimum is zero, with saved rho/momentum/E/B bitwise identical. Task trace and analytic rates exclude extra kicks | Workspace `runs/piece4-device-repro-5641282/viscous-fixed/` and `pytest/results.xml`; control failure in `runs/piece4-wall-probe-5641224/viscous-control/`. The same job intentionally reproduced the tiny-rate HIP failure before its arithmetic fix |
| 4 / tiny-rate pressure preservation | `5fd6342c3`; retain accurate exchange with expm1/BE fractions, update the smaller pressure directly and recover the larger from target anisotropy | Final job `5641333`: frozen refined double/float sweep, both host precisions, both HIP wall policies and 67 application cases | All pass, including unchanged strict monotonicity, 16 new scalar tiny-kick states and two actual HIP states per policy. Double/float encoded-state and second-C2P failures fall 131/209→0; 380,563/264,268 changed maps are clean. No floor/energy repair added and no inherited tolerance changed | Failures retained in `5641224` (scalar) and `5641282` (actual HIP). First fix `5641316` failed only retained float monotonicity; final target reconstruction fixes it. Workspace `runs/piece4-followup-5641333/committed-source.json` verifies all 658 tested files; final binary frozen in `build/piece4-final/`, SHA-256 `83b25c4f…`. Whole-application float remains unqualified |
| 4 / follow-up cost screen | `5fd6342c3` candidate versus initial piece-4 control `c95e61b58`; final candidate frozen separately in `build/piece4-final/` | `analysis/piece4-followup-timing.sbatch`, job `5641349` (8m47s, exit 0); 14 sequential warmup/measured/profile/transport invocations | All passed. Matching time/dt sequences, 81 timestep reductions over 40 cycles and unchanged LF RHS/cycle (42 uniform, 38 active). Nonzero active transport and zero floor/invalid/hard events; all binary/input/source/provenance hashes verified | Million zone-cycles/s/node before→after: uniform 100.679→100.509 (−0.17%); active 70.298→70.021 (−0.39%). Single-pair screen on frozen zero-background-rate decks; no repeatable performance/accuracy claim. Full metrics in workspace `runs/piece4-followup-timing-5641349/summary.json` and `README.md` |
| 11 / isotropic modes | `4926cd2ff`; restore type-2 defaults/unprojected policy and share isotropic selection/normalization; validate generated mode counts | Qualified HIP, staged job `5642148`, 4 cases | All pass: nonzero type-2 force, incompatible-policy rejection and exact fixed-seed type-0/type-2 unprojected equivalence. Type-0/1 enumeration and RNG order retained | Workspace `runs/piece11-qualification-5642148/modes/`; old constructor rejection reproduced in `runs/piece11-repro-5642143/` |
| 11 / conserved work | `16624c806`; compute each fluid's kinetic work from its own conserved density/momentum; independent ideal-EOS guards | Same staged job, 9 direct HIP fixture cases: hydro/MHD/CGL and four two-fluid EOS combinations | All pass with deliberately stale primitives: correct impulse and unchanged thermal energy, density, CGL A, staggered B and isothermal scalar sentinel. Fixture also checks OU initial/held/update coefficients | Workspace `runs/piece11-qualification-5642148/energy/`; pre-fix failures retained in `5642143`. Relativistic source preserved algebraically; SR evolution not requalified |
| 11 / planar axis | `f8686f7bf`; use z as the parallel axis for every planar policy | Same staged job, 2 selected-mode HIP cases with power-law index shells | All selected modes have finite nonzero amplitude; no mode enumeration or RNG retuning | Workspace `runs/piece11-qualification-5642148/planar/`; old planar zero-mode failure retained in `5642143` |
| 11 / full-cycle forcing | `d26d9c156`; atomically move the full-dt kick before CopyCons, retain required boundary/primitive refresh, and accumulate actual kick work | Final stage of `analysis/piece11-qualification.sbatch`, job `5642148`; exact failed-launch retry `analysis/piece11-launch-retry.sbatch`, job `5642155` | 51 selected application cases covered: 50 passed initially; one Slurm socket timeout occurred before launch, then passed unchanged. RK1/2/3 first/every-cycle power, IMEX two-fluid power, OU/restart, AMR, one/four ranks and retained LF/MHD gates pass without tolerance changes | `committed-source.json` verifies all 661 tested files against the final source; frozen `build/piece11-final` SHA-256 `ee8544b3…`. Pre-fix 18 failures in `5642143` include zero first-cycle power and RK2/RK3 overinjection; no broad float or resolved-turbulence qualification |
| 11 / cost screen | Same qualified `d26d9c156` executable; full refresh retained, identical helper and placement to historical post-WO1 | `analysis/piece11-timing.sbatch`, job `5642157` (11m59s, exit 0): 21 sequential warmup/measured/profile/transport invocations, same allocation and frozen flags/inputs | All pass; LF RHS/cycle remains 42 uniform / 38 active; candidate and piece 4 retain 81 timestep reductions versus historical 121. Active transport is nonzero with zero floor/invalid/hard events. Injected/prescribed work is 1.000000000006 for candidate, versus legacy piece-4 1.142694176408 | Million zone-cycles/s/node historical post-WO1→piece 4→piece 11: uniform 116.673→100.798→100.863; active 73.826→70.059→69.523. Piece 11 changes +0.06% / −0.76% versus piece 4; historical gaps remain −13.55% / −5.83%. Single measured samples after warmup, not repeatable small differences or equal-accuracy claims. All hashes verified; full physical-time/work/profile metrics in workspace `runs/piece11-timing-5642157/summary.json` and `README.md` |
| 11 / workflow and archive identities | `e1634b8e7`; record full-cycle collisions and keep historical Figure-13 fixture inputs/source labels pinned; method docs in `44c6da641` | `analysis/piece11-metadata-local.sbatch`, job `5642250`; old fixture checksum mismatch reproduced, then all 78 analysis cases pass | Current forcing/collision labels are explicitly rejected as the old model; archived input/model checksums remain unchanged. No solver change or rebuild. Initial `5642233` was stopped during slow synthetic-fixture creation on Lustre; unchanged tests pass with compute-node temporary files | Workspace `runs/piece11-metadata-local-5642250/` retains XML/logs and verified committed source; interrupted-run note retained in `runs/piece11-metadata-5642233/`. Archive checks do not qualify the rebuilt model's long-time scientific results |
| 5 / weak-field faces | `6ddfdf125`; HLLE/LLF use the neighboring full field for weak-donor A, retain magnetized pressure tensors, and isotropize weak states at fixed thermal energy | `analysis/piece5-check.sbatch qualification`, job `5643328`; qualified HIP double: direct flux/update/C2P, both first-cycle orientations, six smooth evolutions and four CGL/MHD regressions; Serial reference | 8 candidate pytest cases, 1 baseline smooth case and the Serial reference pass; 2 original B4 strict xfails retained. First-cycle ratio is 0.695673–1; analytic direct updates need no EOS repair. Fully magnetized smooth errors converge at orders 1.98/2.00, with all six endpoints exactly equal to baseline. Sharp contacts still violate ratio bounds at cycle 5 | Workspace `runs/piece5-qualification-5643328/summary.json` and `committed-source.json` verify all 663 tested files; frozen `build/piece5-final/src/athena`, SHA-256 `fc423707…`. Baseline failures and interrupted I/O-heavy run retained in `5643286`/`5643292`; completed baseline rerun in final job. No throughput, broad float, turbulence or MPI/restart requalification claim |
| 6 / LF face geometry | `fa840a367`; CT normal field and mean adjacent cell magnitudes in all 12 transport/diagnostic/profile paths, with matching inverse-field and scattering factors and both grad-B references | `analysis/piece6-check.sbatch qualification`, job `5643545`; qualified HIP double: 23 pytest cases, retained grad-B run and Serial reference syntax check pass. Baseline `5643521` retains 7 scientific failures / 8 passing controls | Independent reversal/CT oracles pass; 1D STS sampled seed growth drops from 1.61e6 to approximately 1, with frozen primary fields and conserved energy/moment. Evolved-pressure spatial orders 1.96/2.00, explicit temporal orders 1.002/1.0015; one-e-fold aligned/oblique decay passes. All 13 evolution cases retain exact baseline time/dt/RHS histories | `runs/piece6-qualification-5643545/{summary,matched-work,committed-source}.json`; all 666 source files match Git; frozen `build/piece6-final` SHA-256 `244a4160…`. Paired one-GPU cost job `5643604`: 43.166→43.438M zone-cycles/s over cycles 5–35, excluding final host validation; 6 RHS/cycle, 81 timestep reductions unchanged. One pair only; no speedup or general stability claim. Later slope/closure/controller and mesh qualification remain pending |
| Tier 1 / no-LF companion | Frozen `uniform-no-lf-2n.athinput`, SHA-256 `76c26d47…`; only basename, disabled STS and removal of heat-flux keys differ from uniform control | `analysis/tier1-no-lf-check.sbatch`, job `5643975`; exact pre-WO1 and post-WO1 controls, sequential two-node/sixteen-GPU checks | Both pass 40 cycles with the archived geometry, finite constant dt, no solver output and verified hashes. No numeric LF counter exists in no-LF stdout; absent LF constructor key and source establish no LF process | Workspace `inputs/README.md` and `runs/tier1-no-lf-5643975/`; qualified control only, no paired LF/no-LF performance claim |
| Tier 2 / piece-6 checkpoint | Qualified `fa840a367` binary; q-off, both ledgers, canonical writers plus EOS event log; checkpoint t=6.000234496394158, cycle 19334, SHA-256 `016f9822…` | `analysis/tier2-generate.sbatch` / `tier2-replay.sbatch`, jobs `5643996` / `5644607`; one node/eight GPUs; `analysis/qualify_tier2.py` self-check, native restart/history/field audit, separate full-q probe and output-free timings | Generation and replays pass: zero LF/EOS repair or invalid events, finite positive pressure snapshots, completed hard-bound volume zero; generation energy residual 1.21e-11. Nonzero transport probe matches q-off primary/forcing/live-RNG endpoints bitwise. Healthy cost state, not stationarity or turbulent-accuracy qualification | `runs/tier2-generate-5643996/{checkpoint.json,generation/qualification.json}` and `runs/tier2-replay-5644607/{summary,q-mode-comparison}.json`. Three reference repetitions: 74.951M zone-cycles/s/node mean (0.288% full span), 34 RHS/cycle over cycles 19339–19374; 605.039 s/physical-time and 0.00555492 whole-step s/RHS. Generation allocation 1.056 node-hours including preflight/I/O; no paired speedup claim |
| 7 / transverse slopes | `38b9b5cea`; shared robust VL4 for transverse temperatures in all directions and diagnostic/profile copies; magnetic-magnitude slopes unchanged | HIP double: baseline `5644704` reproduces 12 failures; qualification `5644802` plus byte-identical continuation `5644898` passes 38 pytest cases, grad-B and active conduction controls; host double/float exceptional slopes | Eight hotspots retain initial minima to roundoff, energy drift <=6.11e-14; oblique fundamental converges at orders 2.22/2.17, but under-damps at coarse resolution. Inherited oblique grid 24→32 retains 3% tolerance. Thermal drift uses the conserved energy scale; raw drift and interrupted/failed gates retained | Frozen `build/piece7-final` SHA `86d9cded…`; 669 source files verified. `analysis/piece7-{tier1,tier2,oblique}.sbatch`, jobs `5644939`/`5644984`/`5644998`: Tier 1 101.308→102.837M zone-cycles/s/node, 42 RHS/cycle; Tier 2 74.724→73.780M, 34 RHS/cycle, cost/RHS +1.28% across three pairs; oblique cost unchanged within one-pair screen. Workspace README/JSON/plot retain full rates, work, q-mode identity and provenance. Final bound, turbulent accuracy, AMR and broad float application remain unqualified |
| 8 / perpendicular BGK closure | `476134da8`; shared normal/logarithmic `2*nu_eff` response with subnormal-safe ordering and updated reference consumers. Reassessment `8c7c5ef3a`; legacy-only snapshot proxy `977fa2aa9` (91 analysis checks) | Reproduction `5646185` isolates host/HIP/coupled-decay failures. Qualification `5646195` passes 21 pytest cases, standalone reference compile and both inherited collisional STS decks; identical-binary continuation `5646222` repairs only a workspace table-reader assumption and passes collisionless identity | Independent signed two-pressure continuum: resolved errors 8–20% → <0.48%, approximately fourfold refinement improvement, energy drift <=8.89e-16, zero safety events. Collisionless primary fields and all raw history columns match Piece 7 exactly; native A is not directly exported by that table control | Frozen `build/piece8-final` SHA `d6e9ddef…`, all 833 source files verified. Jobs `5646223`/`5646257`/`5646260`: Tier 1 102.783→102.637M zone-cycles/s/node (42 RHS/cycle); Tier 2 74.248→74.166M (34 RHS/cycle, cost/RHS +0.11%, within observed variation); oblique 43.407→43.369M on one GPU (6 RHS/cycle, one pair). q-full/none restart state and live RNG match bitwise. Workspace README/JSON/plot preserve physics, timings and reader failure. No threshold-switching, turbulent-trajectory, AMR or broad float qualification |
| 8 follow-up / response overflow | `70e1002bc`; retry half of a newly overflowing denominator before logarithmic cancellation; retain ordinary driver arithmetic. Coupled-decay gates now check full pressure profiles | Reproduction `5646292` fails host double/float and HIP endpoints; qualification `5646298` plus same-binary launch retry `5646308` passes all 21 solver pytest cases, standalone compile and both inherited STS decks. Fresh 91 analysis checks pass with captured stdout/JUnit | Resolved full-field error <0.476%; checkerboard negative control fails despite unchanged fundamental amplitudes. Tier 2 original/fixed conserved/CT fields including ghosts, forcing and live RNG match bitwise; q-full/none identity and zero new floor/invalid repairs pass | Frozen `build/piece8-final-v2` SHA `25cd995e…`, all 833 source files verified. Tier 1 `5646316`: 102.559→101.791M, 42 RHS/cycle, cost/RHS +0.755%. Tier 2 `5646352`: 74.436→74.665M, 34 RHS/cycle, cost/RHS −0.306% within variation. Oblique `5646361`: one pair +0.059% cost, 6 RHS/cycle. Original scheduler failures retained; no speedup claim. Arbitrary-exponent arithmetic beyond the introduced overflow boundary and previous physical scope limits remain unqualified |
| 9 / final-stencil bound | `f8000cc01`; coupled CT/density/VL4/grad-B temperature rows, background collisions, fresh halo in existing scratch and shared representable-sound-speed fallback; current-CFL safety retained | Parent reproduction `5646412`: 8 failures / 1 passing control. Main qualification `5646455`: 65 passes with independent frozen Jacobian/finite-difference/spectrum/RKL checks and HIP evolution | Rejected coefficient-five envelope 54.8585 < spectral radius 72.9081; new row norm 99.7224 <= envelope 99.8419. Positive/complex eigenvalues retained and RKL eigenmode refinement judged against exp(dt*lambda); nonnormal transients are not bounded. Frozen coefficients only, not nonlinear/AMR stability | Early `5646418` fixture failures traced to degenerate slopes, an unintended magnetization floor, unequal layout endpoints and excessive reference timestep; corrected fixtures preserve scientific tolerances. Full diagnoses under workspace `analysis/piece9-preflight/` |
| 9 / grad-B arithmetic | `486bbb570`; match flux division/projection and x/y/z summation order so anisotropy cannot amplify a residual omitted by the bound | Fresh final qualification `5646500`: all 66 cases pass across nine gates, including new source-derived cancellation regression, actual 2D/3D single/split runtime bounds, heated refresh, MPI/restart and retained closure/hotspot checks | Reversal spatial errors decrease about fourfold; STS full-profile errors converge with explicit-reference uncertainty <4.41% of the finest STS error. Closure full-field error <0.483%; hotspot energy error <=1.03e-13, positive pressures and zero repair/invalid events. Cancellation regression is host arithmetic evidence, not a reproduced GPU instability | Final binary SHA `413bcff4…`; all 923 application and 1,455 Kokkos source files verified against Git. Workspace `runs/piece9-qualification-5646500/` and `analysis/piece9-{final-math-summary,qualification-5646500-physics-summary}.json` retain spectra, rerunnable matrix checks, full-precision errors and provenance |
| 9 / convergence and cost | Same qualified `486bbb570` binary frozen in `build/piece9-final`; current-CFL safety and retained checkpoint unchanged | Oblique continuation `5646516`: finest Fourier-amplitude error 0.599%, spatial orders 2.22/2.17, half-step contamination 0.234% of intermediate spatial error. Total 67 physics cases pass. Sequential Tier 1/2 jobs `5646534`/`5646545`: all 31 invocations pass, three measured pairs per workload plus separate profiles | Tier 1 LF 100.827→130.202M zone-cycles/s/node, 42→26 RHS/cycle; no-LF 558.887→558.731M. Tier 2 74.011→50.824M, 34→50.457 RHS/cycle: throughput −31.33%, amortized whole-step s/RHS 0.00562546→0.00552006. Uniform savings do not predict developed-state cost | Tier 2 s/physical-time 612.815→892.387 using each history; full-q/q-off conserved/CT/forcing/live-RNG endpoints match bitwise, nonzero transport and no new floor/invalid events. Heat-flux profile s/RHS rises 20.0% / 11.6% in Tiers 1/2; new bound cost retained explicitly. Workspace README and per-run summaries/reviews retain spreads, physical advances, nested profiles and provenance. No equal-accuracy turbulence or general speedup claim |
| 9 follow-up / exceptional arithmetic | `4b5ea7f67` retains the represented grad-B residual in the cold bound; `7011864b9` preserves subnormal parallel scattering; `afb2be2d6` makes the local-speed fixture satisfy E/A recovery and active pressure walls. Final archived source `171e81eb7` | Old-source reproduction `5646561`: intended cold-row failure and parallel double/float failures, with a passing legacy device control. Final HIP qualification `5646611`: all 69 tests pass across nine gates, including 14 signed device closure cases and safe/fast admitted zero-cycle local-speed checks | Old cold row 132874.47 was below a capped reverse-coupling row 177087.31; corrected row 177165.96 covers it. The double parallel response was 14.04% above the cold coefficient because scattering vanished. These are local arithmetic counterexamples, not evolved-instability reproductions | Binary SHA `3ffebd3e…`; all 923 application and 1,455 Kokkos files match Git. All 89 retained numerical artifacts match `5646500` byte-for-byte. Preliminary `5646569` is preserved separately. Workspace `analysis/piece9-reaudit/` and final qualification review retain evidence; frozen-Jacobian, nonlinear/AMR and full-float scope limits remain |
| 9 follow-up / retained results and cost | Exact qualified binary frozen in `build/piece9-final-v2`; original Piece 9 freeze and all old timings retained | Oblique `5646631` passes with all eight numerical files byte-identical to `5646516`; final total 70 tests. Tier 1/2 `5646643`/`5646649`: all 31 invocations pass, three measured pairs per workload, separate profiles and transport replay | Tier 1 LF 130.649→136.862M zone-cycles/s/node (+4.76%), no-LF 559.014→559.716M; 26 RHS/cycle unchanged. Tier 2 51.189→51.422M (+0.46%, within observed variation), 50.457 RHS/cycle unchanged; whole-step s/RHS 0.00548095→0.00545585 | Tier 2 native conserved/CT/ghost arrays, forcing/live RNG, ledgers and time/dt/work histories match exactly; q-full/q-off identity, nonzero transport and zero floor/invalid events pass. Intermediate limiter/hard-bound detections remain nonzero. Timings/profile scopes, spreads and independent reconciliation are in workspace summaries/reviews. No isolated source attribution or general speedup claim |
| 10 / independent STS safety | `28d748f49`; preserve absent-key current-CFL behavior, load/validate explicit safety on fresh starts and restarts, scale only STS budgets and retain inactive sentinels. Final source `ad99fad72` | Qualification `5646694`: all 50 tests pass across eight gates, including host double/float budgets, actual HIP startup/restart/refresh, ordinary diffusion, MPI, closure and hotspots. Exact binary frozen as `build/piece10-final` | Safety 0.3/0.9 reversal and finite-rate mirror-crossing full pressure fields converge at matched nominal outer budgets. Larger safety reduces coarse reversal RHS 20→12 but raises coarse pressure error about 42–45%; both pass unchanged tolerances. Threshold reference differences are 33.35–34.72% of finest reported errors: coarse first-order evidence, not turbulent accuracy | Workspace `analysis/piece10/physics-5646694-independent-review.json`, runtime review and seven-matrix RKL probe retain raw reconciliation and limits. Exploratory `5646678` fixture/provenance issues are retained; final source archive matches all 926 application and 1,455 Kokkos files. Tier 2 generation `5646710`, same-safety comparison `5646841` and separate 0.3/0.9 cost screen `5646854` pass. No larger default adopted |
| 10 / checkpoint and paired cost | New canonical producer with explicit safety 0.3; frozen t=6.000305031020128, cycle 19315, SHA `bb475c16…`; old checkpoint and timings retained | Generation `5646710` reaches t=6.05 in 1.371 node-hours including I/O. Comparisons `5646841`/`5646854` each pass 14 invocations with three alternating measured pairs and separate profiles | Parent/candidate at 0.3 match native conserved/CT/ghost state, forcing/live RNG, ledgers and all q-off history columns exactly. Rates 52.053→51.934M zone-cycles/s/node (−0.229%, within observed variation), both 50.4 RHS/cycle. Safety 0.3→0.9 reduces work 50.4→30 RHS/cycle, rates 51.723→75.130M (+45.255%), and s/physical-time 866.138→596.288 using each history | Both safeties pass q-mode identity, nonzero transport, positivity and zero new LF floor/invalid or EOS repair checks. Intermediate hard-bound crossings remain nonzero; strict checkpoints and completed-state hard-bound volume are clean. The legacy `lf_hwproj` value stays zero and is not a projection measurement. Safety screen cost/RHS rises 5.43027→6.28058 ms as fewer RHS share outer-step work; separate heat-flux ms/RHS 1.70897→1.70781. Workspace completion/review records retain raw reconciliation, spreads and profile scopes. No equal-accuracy turbulent speedup claim; explicit 0.3 remains the active-milestone setting. Tier 1 unnecessary without recurring full-domain changes; active milestone remains next |
| 10 follow-up / independent review | Production source and frozen executable unchanged (`ad99fad72`) | Fresh connected-code review; all 50 retained JUnit results and 12 output-free timing invocations independently reconciled; no new simulations | No solver defect found. Raw threshold snapshots give successive-difference temporal orders 0.963–1.039, independent of the shared fine reference; reversal errors and stage/work changes reproduce. Restart, generic diffusion and zero-LF collision paths reviewed | Corrected only the completion record's warmup wording: five completed updates, then 35 measured updates, then one guard (41 total). Timing arithmetic was already correct. Original report retained as `analysis/piece10-completion-original.json`; `analysis/piece10-reaudit.json` binds this review and correction. Previous qualification limits and explicit 0.3 milestone setting remain |
| 12 / mesh and boundary consistency | `de75846b9` corner refill; `76ca6cb4f` shared LF fluxes; `af2aa9480` boundary lifecycle; final source `c8d57ea61` | CPU `5647453` and HIP `5647452`: 25/25 each. Final MPI `5647477`: 2/2, including two-node/16-rank AMR churn. External instrumented CPU `5647461`: 230 refreshes preserve active E/B/CT bitwise | Both prolongation modes show decreasing STS/explicit pressure errors; native energy/moment drift <7e-16. Shear restart primary/CT/ghost state and all integer counters are exact, with 1522 RHS on both paths. One qpar-work sum differs by 1 ulp from MPI accumulator regrouping; the extra bitwise-all-diagnostics assertion failed, and that job status is retained | Workspace `analysis/piece12/README.md` and `qualification-summary.json` retain parent negatives, intermediate failures and final evidence. Local row bounds do not establish composed nonlinear AMR stability; face q-work diagnostics precede reconciliation. No performance claim or new uniform-periodic device work; frozen Tier 3 executable/dependencies remain intact |
