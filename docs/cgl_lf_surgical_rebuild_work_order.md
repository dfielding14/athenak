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

**Status:** piece 0 and requested throughput logging complete; solver repairs pending.
**Workspace:** `/lustre/orion/ast207/proj-shared/dfielding/CGL/surgical-rebuild`.
**Next:** piece 1 — guards and parameters.
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

- [ ] **1 — Guards and parameters.** Reject passive mode pending redesign;
  preserve LF boundary guards and validate existing floors/rates. Gate: invalid
  inputs fail clearly and valid active inputs retain their behavior.

- [ ] **2 — EOS repairs.** Fix density/magnetization floors and persist repaired
  state. Preserve pressure ratio in the A decoder; preserve magnetic moment,
  momentum, and energy in the moment decoder. Gate: independent pressure/energy
  identities, idempotence, stored A, and focused double/float checks.

- [ ] **3 — Reconstruction and FOFC.** Repair both CGL pressure floors in every
  directional PPMX/WENO-Z wrapper and detect incoming nonfinite E/A before repair.
  Gate: wrapper tests, NaN/Inf injection, steep-pressure evolution, and MHD controls.

- [ ] **4 — Collisions and walls.** Implement thresholds, relaxation, additive LF
  scattering, canonical wall encoding, and all scheduling/AMR/input consumers
  coherently. Apply full-dt rates once per cycle and walls at prescribed checkpoints,
  including the mandatory fluid-firehose wall without optional limiters. Update
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
