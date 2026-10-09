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

**Status:** pieces 0–4 complete.
**Workspace:** `/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild`.
**Next:** piece 11 — forcing mode, conserved-momentum work and scheduling.
Keep the plan's piece numbers; execute forcing (11) immediately after collisions
(4). The uniform active milestone precedes optional optimization and later models.
Passive work (14) does not depend on mesh work (12) or optimizations (13).
Bring conditional chunking (15) forward if active/passive qualification requires it.

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

- [ ] **11 — Forcing.** Separate mode selection, conserved-momentum work, and
  scheduling fixes; move task registration and required refresh together. Gate:
  RK1/2/3 full-step kick counts, power, unchanged thermal energy, OU cadence/restart,
  and affected hydro/MHD/two-fluid consumers. Preserve seeds and prescribed power.

- [ ] **5 — Weak-field faces.** Repair neighboring-field use in HLLE/LLF while
  retaining A. Gate: first-update contamination, both contact orientations, smooth
  interface convergence, and unchanged predicates for the two B4 expected failures.

- [ ] **6 — LF geometry.** Use CT normal B and mean adjacent cell magnitudes;
  do not renormalize to unit length. Gate: independent face fluxes, reversals,
  frozen-field conservation, and oblique decay.

- [ ] **7 — Transverse slopes.** Apply VL4 to temperature gradients, retaining
  unlimited magnetic-magnitude slopes and consistent diagnostic/profile copies.
  Gate: hotspot minima/energy, oblique convergence, exceptional slopes, and a
  separate paired timing of the added arithmetic.

- [ ] **8 — Perpendicular closure.** Correct the `2*nu_eff` response in normal
  and logarithmic paths and reference consumers. Gate: collisionless identity,
  strong-collision limit, and resolved coupled decay. Across 6–8 use the common
  conservative reference step and oblique background-coefficient kinematic screen;
  verify actual dt/RHS counts instead of assuming explicit integration fixes dt.

- [ ] **9 — Final-stencil timestep bound.** Adapt the coupled CT/density/VL4/grad-B
  rows directly, with fresh halo values in existing scratch storage and safe
  exceptional arithmetic. Retain the old safety multiplier. Gate: rerunnable
  independent rows/eigenvalues/RKL amplification, checkerboards, density jumps,
  reversals, heated refresh, MPI agreement, and explicit-reference convergence.

- [ ] **10 — STS safety.** Separate STS safety from advective/explicit CFL,
  including restart loading and both budget selections. Gate: parameter bounds,
  analytic odd-stage counts, stability, and temporal convergence before accepting
  a larger factor. Record work-count changes separately from evaluation cost.

- [ ] **Active milestone.** Run the plan's selected physical suite once, including
  resolved waves/decay, finite-rate relaxation, forcing/restart, HIP wall encoding,
  and resolved active turbulence. Extend the local-coefficient `64 x 64` hotspot
  to at least three temporal refinements at the same grid/endpoint, preserving
  duration and conservation checks. Report worst-cell absolute/relative and
  thermal-energy-scaled pressure errors, norms, extrema, and repair activity.
  Require decreasing temporal errors.
  Compare once with qualified `2afaa0fa1` at matched active physics, merging and
  chunking off. This qualifies the first usable deliverable within tested limits.

- [ ] **12 — Mesh support.** Separately repair physical corners after prolongation
  and synchronize shared LF face fluxes. Grow both buffers before receives, with
  sufficient 3D full-face capacity and existing request completion preserved.
  Gate: energy/moment conservation, both prolongation modes, 2D/3D explicit/STS,
  CPU/HIP, one/four ranks, actual regrid/restart, repair accounting, and shearing.
  Reserve two-node churn for final mesh qualification.

- [ ] **13 — Optional optimizations.** Profile the corrected target and select
  the smallest useful removal/fusion/message change. Compact communication follows
  mesh correctness. Gate each retained change on full-precision primary, ghost,
  coarse, scalar, and restart/RNG equivalence plus paired GPU timing. A simple
  deletion may stand within timing noise; larger machinery needs repeatable benefit.

- [ ] **14 — Passive model.** Complete J/A outside LF and U/moment inside LF,
  including EOS, fluxes, forcing, initialization, outputs, restart encoding, and
  guards before enabling support. Preserve the native isothermal flow kernel.
  Gate: matched-timestep flow identity, independent thermal evolution/work, HIP
  wall decoding, forced MPI restart, all supported reconstructors, and nonlinear
  temporal refinement. Initial scope is uniform periodic; document absent shock heating.

- [ ] **15 — Conditional LF chunking.** Implement only if active/passive accuracy
  or positivity needs smaller LF intervals at fixed outer dt. Reproduce the failing
  state, preserve forcing/collision cadence, and demonstrate decreasing pressure
  errors across at least three LF intervals, with conservation and restart checks.
  Surviving with cap 4 is insufficient; an unconverged cap-1 run is not an oracle.

**Measurement and acceptance**

Use safe/none/weighted LF with strict checks and explicit environment-override
control. Keep pressure-work and injected-work recording enabled on research
timings. Use warmup and matching output-free intervals: one alternating pair
screens a patch; three establish a baseline or support a speedup. Run timing
jobs without competition; profile separately with detail off initially.

Report zone-cycles/s/node, LF evaluations/cycle, elapsed/physical time, and
amortized whole-step time/RHS. Final acceptance reports node-seconds to the common
endpoint at demonstrated accuracy. Set additional endpoint-error criteria before
selecting performance settings. The frozen-row bound does not prove nonlinear or
AMR stability; hotspot convergence does not certify turbulent threshold switching;
once-per-cycle collisions do not make the full split second order.

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
