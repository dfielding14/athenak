# Work Order 2: CGL + Landau-fluid + RKL2 STS follow-ups (result-changing performance and deferred redesigns)

**Repository:** AthenaK, branch `CGL-STS-LF`. WO1 was merged through [PR #21](https://github.com/dfielding14/athenak/pull/21) on 2026-10-06 at `ac33b6b04b8093c2d65e87a30db93148fa65714d`; its final implementation/evidence commit is `543a1d7ec967a592f155fdd0d2213a5dfbe20b5c`.

**Execution handoff:** carry out WO2 on the actual GPU machine. This document updates the work order only; no WO2 implementation or GPU runtime validation has been performed. Start from the merged WO1 state and establish the device baseline below before changing code.

## 0. Starting state, rules, and GPU baseline

### Decisions and current triage

- **Retain conservative A.** The extreme unresolved sharp-contact failure is an accepted limitation shared with the checked Majeski and archived Squire implementations. Keep both B4 tests as `xfail(strict=True)` with their existing numerical bounds. Do not replace A with Q, redesign B4, or treat the old Q prototype as active work. See the [accepted decision](docs/validation/wo1/review/weak-field/README.md) and [reference comparison](docs/validation/wo1/review/weak-field-reference/STATUS_UPDATE.md). The Squire archive is not a verified publication revision; smooth-contact convergence is not general turbulence validation.
- **P1 is absent from production.** Its compact exchange and skipped magnetic-boundary updates were rejected after unexplained active-cell changes at a coarse/fine outflow interface. Diagnose the separate Task 0-P1 below before retrying that optimization or restructuring the multilevel path.
- **Task 6 is already implemented.** CGL primitive prolongation handles both A and magnetic-moment representations through the CGL AMR projection helpers. WO1 verified CPU behavior; WO2 must verify it on the GPU and preserve it through Task 3, not add a duplicate implementation or remove a nonexistent fence.
- **Tasks 4 and 5 remain fenced.** Passive CGL and LF inflow/user boundaries require the redesigns below. The current boundary guard covers explicit and STS LF, including multilevel configurations.
- **Task 9 needs triage.** Post-AMR wall projection, log-space A encoding, and overflow-safe heat-flux cap arithmetic already exist. Initial-state wall handling and the remaining optional items still need inspection.
- **Paths have changed.** The LF implementation is `CGLLandauFluid` in `src/diffusion/cgl_landau_fluid.cpp`; locate functions by name rather than old line numbers. Existing uniform-grid fast paths and full-path optimizations must be preserved where applicable.

For every task, record **already done**, **adapted**, or **implemented**, with source and validation evidence. The historical review described `2aee609`; its bug-status claims are not the current implementation. The local root WO1/review Markdown files are not shipped with this branch. This order is self-contained; use the committed [WO1 per-task report](docs/cgl_lf_changes.md), [evidence index](docs/validation/wo1/ARCHIVE.md), and current source as the authoritative starting point.

### Ground rules

1. Create a new `c/` branch from the merged baseline. Make one commit per executed numbered task, naming the task ID. Keep builds and applicable tests passing; the two accepted B4 expected failures remain explicit exceptions.
2. Triage before editing. Reuse existing helpers and tests, keep duplicated formulas/diagnostics in lockstep, prefer independent analytic or numerical references, and list changed expected values. No unrelated cleanup or new runtime switches beyond those specified here.
3. Retain the verified RKL2 coefficients and odd-stage rule. Task 1 alone changes the STS stability budget/CFL factor; Task 2 alone introduces optional merged sweeps, with the default off. Preserve the WO1 collision schedule and wall semantics unless a specified task explicitly changes them and verifies the new split.
4. Most WO2 tasks may change results by design. Report the size and cause of each change, with accuracy/convergence evidence against the frozen WO1 baseline. That permission does not excuse unexplained discrepancies. Where bitwise equality is required, compare on the same backend, compiler/options, Kokkos revision, hardware, and rank layout; do not demand CPU/GPU byte identity.
5. Use the actual GPU for runtime correctness and performance measurements. Record source/binary/input hashes, device model, driver/runtime/compiler, Kokkos revision, precision, MPI implementation, rank/device binding, block geometry, warm-up, and repeats. Use at least three timed repeats with fixed inputs and outputs; synchronize device work appropriately. Report median wall time per cycle, STS stage counts, and available kernel/communication timings without double-counting nested regions. Keep baseline and candidate measurements on the same machine/allocation.
6. Preserve scientific tolerances and the collisionless shear/AMR fixtures. Do not turn on collisions, reduce amplitudes, disable strict diagnostics, or weaken checks to hide a failure. Diagnose baseline failures before attributing them to WO2. Keep full-float portability work within optional Task 9; CUDA compilation alone is not GPU execution.

### GPU baseline before edits

WO1's [final CI run](https://github.com/dfielding14/athenak/actions/runs/37514510289) passed 213 CPU checks plus two B4 expected failures, three additional AMR checks, seven MPI checks, and CUDA 12.6 compilation. **GPU runtime validation remains to be done here.** Full application single-precision compilation has pre-existing portability failures; focused float checks of the changed WO1 math/EOS paths passed.

1. Confirm the checkout contains the WO1 merge, record its exact SHA, and initialize Kokkos at the committed revision (`git submodule update --init --recursive`). Inspect the actual GPU/allocation and select the correct Kokkos backend/architecture; do not blindly reuse CI's `AMPERE80` setting. Use separate CPU and GPU build directories or checkouts, with Release double precision initially.
2. Run the existing WO1 acceptance workflows using an explicitly supplied GPU executable. `scripts/run_cgl_lf_validation.sh` now delegates to `scripts/cgl_lf_workflow.py full`; use `--no-build` so it cannot silently create a CPU executable. For an already configured GPU build in `tst/build`, for example:

   ```sh
   python3 scripts/cgl_lf_workflow.py full \
     --build-dir "$PWD/tst/build" \
     --athena-bin "$PWD/tst/build/src/athena" --no-build \
     --output-dir "$CGL_WO2_RUNS/wo1-gpu-full"
   python3 scripts/cgl_lf_workflow.py paper-smoke \
     --build-dir "$PWD/tst/build" \
     --athena-bin "$PWD/tst/build/src/athena" --no-build \
     --output-dir "$CGL_WO2_RUNS/wo1-gpu-smoke"
   ```

   Set `CGL_WO2_RUNS` to a persistent results directory first. These commands are handoff instructions, not claimed GPU passes. Include field reversal, corrected density-jump, hot-spot, limiter-map, collision-schedule, fast-speed, resolved wave/decay, and conservation checks. The corrected density test uses uniform temperature, normal B, contrasts 10/50/200/1000, paired $10^{-6}$ seeded/unseeded runs, and less than about tenfold difference growth over three cycles; the old isobaric/20-cycle wording is superseded.
3. Run the existing CGL GPU AMR tests (`tst/test_suite/cgl/test_cgl_amr_gpu.py`) on that GPU executable, including primitive prolongation, 3D churn, and restart through regridding. Use the repository's `tst/test_suite/testutils.py` setup and test working-directory conventions. Retain CPU reference checks and run the collisionless shear, oblique-decay, timestep-refresh, and forcing MPI regressions with appropriate GPU/MPI builds and rank/device binding. Record skipped coverage explicitly.
4. When the allocation supports it, run the opt-in `test_cgl_amr_mpi_gpu.py` with `ATHENAK_RUN_MPI_GPU=1` and a valid `ATHENAK_MPI_GPU_LAUNCHER`. That harness appends Slurm-style `-N`, `-n`, and `--ntasks-per-node` flags; use a compatible launcher and allocation. If only one GPU is available, complete single-GPU validation and report the multi-GPU gap; do not count skipped tests as passes.
5. Generate and freeze fresh WO1 binary, history, and full-precision restart outputs on this machine before any WO2 edit. Use the fixed inputs in `docs/validation/wo1/inputs/`, including uniform 1D/2D/3D LF, SMR, pure CGL, shear, and `smr_outflow.athinput`, and retain a resolved turbulence benchmark. The committed archive stores hashes and summaries, not all original raw outputs; old Mac `/tmp` paths and CPU hashes are not portable GPU references. Adapt a working copy of the retained harness, never overwrite the historical archive. Check repeatability before using byte comparisons.

**Task order:** GPU baseline → diagnose 0-P1 on the frozen baseline → 1 → 2 → 3; then remaining work in 4–7 as appropriate; 8–9 remain optional and last. Resolve 0-P1 before accepting any communication/full-path optimization. Task 6 is validation of existing support. Tasks 1 and 2 are performance priorities, but their expected savings must be measured on this GPU.

## 0-P1. Resolve the rejected communication optimization

**Question:** during LF stages, can we exchange only IEN/IAN and skip magnetic physical-boundary updates without losing data needed at refinement/physical boundaries? The two changes were tested together in WO1, so their individual effects remain unknown.

**Evidence:** the candidate passed 21 of 24 configurations but failed SMR/outflow in serial and MPI on one/four ranks. By cycle 4, active-cell energy differs by up to $3.45\times10^{-3}$ and magnetic components by about $2\times10^{-3}$; total-energy history first differs at cycle 1. This is not merely a ghost-output mismatch. The original candidate is fully removed. See the [rejected patch](docs/validation/wo1/P1-tested-rejected.patch), [fixed input](docs/validation/wo1/inputs/smr_outflow.athinput), [difference analysis](docs/validation/wo1/runs/p1/difference-analysis.json), and [T-P1 report](docs/cgl_lf_changes.md#t-p1-rejected-non-bitwise-communication-optimization).

**Work:**
- Reproduce the discrepancy with the frozen input and baseline. Historical patches may need adaptation to the merged source; do not reapply the combined patch blindly.
- Isolate two candidates: compact IEN/IAN cell exchange with magnetic BCs retained, and skipped magnetic BCs with full cell exchange retained. Keep diagnostic controls in scratch builds; do not add production runtime switches.
- Trace the first differing ghost and active values through send/receive packing, physical boundary filling, restriction/prolongation, and flux correction. Inspect density, momentum, scalars, cell/face magnetic fields, and A/μ representation. Determine the actual dependency before proposing a correction.
- If the baseline has a real boundary defect, fix and validate that defect separately from the performance optimization. If the candidate omitted necessary work, retain that work. Keep shearing-boundary dependencies explicit.

**Acceptance:** an explained cause and a cheap regression that detects it, with serial and MPI coverage and actual GPU execution. A pure communication optimization must preserve results bit-for-bit on the same execution configuration. If a separately justified correctness fix changes the reference solution, report the changed quantities and conservation/convergence evidence before freezing a new baseline. Measure transfer volume and GPU timing only for a validated candidate. Do not waive the discrepancy because WO2 allows other result-changing designs; if a safe optimization is not established, retain full communication and report the finding.

---

## 1. Per-cell Gershgorin timestep bound for the LF operator, and `sts_safety`

**Why:** the historical turbulence setup used about 30 STS stages against 2 RK stages per cycle. Treat this as motivation and measure the actual GPU baseline before quoting a stage-count or speedup claim.
- The STS dt is multiplied by the advective `cfl_number` in `src/mesh/mesh.cpp` (`Real process_dt = (cfl_no)*(process.ExplicitDt());`). That adds a factor of √(1/cfl) to the stage count, about 1.8–2× at cfl = 0.25–0.3.
- The `cfl_no` factor cannot simply be removed while the estimate ignores real stiffness. WO1 T-D5 added a conservative factor for this.
- **Do not loosen `fac = 1/(2·ndim)` on its own.** An arbiter showed it is the sharp worst case for grid-scale divergence-free field patterns, where the spectral radius equals 4d·χ/h².

**Change:**
1. Replace the per-cell estimate in `CGLLandauFluid::NewTimeStep` (`src/diffusion/cgl_landau_fluid.cpp`) with a rigorous **row-sum (Gershgorin, ∞-norm) bound of the discrete operator as implemented after WO1**. This replaces both `fac` and the WO1 T-D5 factor R_i for the LF branch.
   - Starting point from the review arbiter:
     λ_i ≤ h⁻² Σ_faces (K_f/ρ_i)·(2 b_n² + |b_n| Σ_t |b_t|)
     - K_f = ρ_f χ_f is the face conductivity for each of the T∥ and T⊥ parts.
     - For anisotropic spacing, use the actual Δx per direction, not h_min.
   - For the μ/p∥ coupling, include the factor B_i/B̄_f and the (1 − B_i/B̄_f) cross term that feeds q⊥ into p∥.
   - Include the grad-B term's contribution.
   - With cross-face cancellations, the bound should be about 4–6·χ/h² for smooth fields: 4.0 field-aligned, 5.0 at 45° in 2D.
   - **You must derive the exact expression for the final WO1 stencil**, which includes the T-D1 face field and the T-D2 limited transverse gradients.
2. Verify numerically, e.g. with a Python replica of the frozen linearized operator, that the bound is ≥ the true spectral radius. Test at least these cases in 2D and 3D, with and without the T-D2 limiter active:
   - uniform b̂ at 0°, 30° and 45°;
   - a random grid-scale field;
   - a tanh current sheet 1–2 cells wide with a 1–3% guide field;
   - density jumps up to 10³;
   - a turbulence-like field (e.g. from a WO1 turbulence smoke run).

   Report the ratio bound/λ_max for each case.
3. Add `<time>/sts_safety` (default **0.9**, must lie in (0, 1]) and use it instead of `cfl_no` for STS-managed processes in `Mesh::NewTimeStep`. Keep `cfl_no` for explicit parabolic processes. Also use it in WO1's `Mesh::RefreshParabolicTimeStep`.
   - Keep a real margin: RKL2 grows catastrophically just past the stability edge (a 5% overshoot gives |R| = 72 at s = 15).
4. Remove WO1's conservative T-D5 factor from the LF branch, since the new bound subsumes it. Do not touch `fac` for the other (isotropic) conduction branches.

**Acceptance:**
- All WO1 stability tests still pass: field reversal, density jump, hot spot, turbulence smoke.
- Eigenmode, decay and convergence tests stay within tolerance.
- Report the stage-count reduction per sweep and the wall-clock speedup on these inputs:
  - `inputs/cgl_lf_paper/cgl_lf_paper_turb_active.athinput`
  - one wave input
  - one SMR input

## 2. Merge back-to-back half-sweeps (optional mode)

**Why:** The post sweep of cycle n is immediately followed by the pre sweep of cycle n+1. RKL2 stages scale as √(duration), so one sweep of (dt_n + dt_{n+1})/2 costs about 1/√2 of two separate half-sweeps, roughly 30% fewer stages. RKL2 is still the right integrator family: RKL1 is first order, and RKG2 has a smaller stability interval.

**Change:** add `<time>/sts_merge_half_sweeps` (default **false**). When it is true:
- Skip the post half-sweep unless an output, restart dump, AMR step or tlim/nlim stop falls between cycles.
- At the start of the next cycle, run a merged pre-sweep of (dt_prev + dt)/2, with its stage count computed for that length.
- Before any output, restart, AMR or final step, run the real post half-sweep so every output is at a consistent Strang state.

**Constraints:** the WO1 collision schedule ("rates once, walls everywhere") must still hold. Specifically:
- Collision rates run exactly once per cycle with that cycle's dt.
- The hard-wall projection runs after every sub-operator: each sweep (merged or not) and each hyperbolic step.

Document where the end-of-cycle rates call moves when the post sweep is skipped, and justify the splitting order. Also check the interplay with the fresh-dt refresh (WO1 T-D6) and with the dt change between cycles.

**Acceptance:**
- With the option off, the results are **bitwise identical** to before.
- With it on:
  - smooth tests converge at the same order;
  - outputs agree with the unmerged run to within the splitting error, which should shrink as O(dt²) when dt is halved;
  - the stage-count and wall-clock savings are reported;
  - restart and AMR produce correct, consistent states.

## 3. Restructure the multilevel "full" STS path

**Why:** SMR and AMR runs use the full path (`src/mhd/mhd_sts.cpp`, about 173-273; `src/mhd/mhd_tasks.cpp`, about 92-131). WO1 T-P2 removed its bitwise-safe waste. What remains is roundoff-level: it converts A↔μ in active cells at **every stage** and runs the full `ConToPrim` every stage, which recomputes the frozen bcc.

**Change:** make the full path mirror the two-var path:
- convert once per sweep at Begin and End, on all cells including ghosts and coarse buffers;
- use the lightweight primitive refresh;
- restrict flux correction to IEN and IAN;
- handle restriction and prolongation of the μ representation correctly during the sweep. μ = p⊥/|B| is a linear conserved density with B frozen, so conservative restriction and prolongation of μ are valid. Check prolongation specifically.

Reuse the existing representation-aware AMR projection helpers and Task 6 coverage. Merge the two code paths only if it simplifies the validated implementation; do not create a second projection system.

**Acceptance:**
- The SMR LF tests from WO1 T-G6 agree with the old path to roundoff-level differences, with the magnitude reported. Include the resolved 0-P1 coarse/fine outflow case so stale-ghost changes cannot be mistaken for roundoff.
- Total E and ∫μ dV are conserved.
- The per-stage data volume reduction is reported.

## 4. Passive-mode redesign (removes the WO1 T-E2 fence)

**Why (review M7):** passive mode takes its momentum flux from isothermal MHD, but the energy flux keeps the full CGL stress work. The recovered thermal pressures therefore pick up a spurious v·(∇p_iso − ∇·P_CGL) term. Integrated over a periodic box, that term exactly cancels the Δp viscous heating a passive-Δp run is meant to measure.

**Change:**
- Redesign passive mode so that ρ, v and B follow isothermal MHD exactly, while p⊥ and p∥ follow the CGL equations (Squire et al. 2023 eqs. 2.3–2.4), including LF heat fluxes, collisions and limiters.
- Recommended approach: evolve the two double-adiabatic invariants as conserved scalars, i.e. A (already present) plus a second one, e.g. ρ·ln(p∥B²/ρ³), in the IEN slot or a new slot.
  - p⊥ and p∥ are recovered from these invariants in the c2p.
  - The LF sweep and collisions convert through the primitives: invariants → (p∥, p⊥) → (E_th, μ) → sweep → back.
  - The energy slot then no longer carries total energy in passive mode. Document this, and make the history outputs reflect it.
- Consider shock-heating semantics: invariant advection gives no heating at shocks. That is acceptable for passive-Δp turbulence, but document it.
- Remove the T-E2 fatal error, and re-enable `inputs/cgl_lf_paper/cgl_lf_paper_turb_passive.athinput`.

**Acceptance:**
- A passive run's ρ, v and B are **bitwise identical** to an isothermal-MHD run with the same forcing seed.
- In a smooth linear test, passive p⊥ and p∥ match the linear CGL(-LF) evolution for the prescribed flow.
- The box-integrated Δp viscous heating is no longer cancelled.

## 5. Representation-aware LF inflow and user boundary conditions (removes the WO1 T-E1 fence)

**Current guard:** `src/mhd/mhd.cpp` rejects inflow/user boundaries for explicit and STS LF, including multilevel runs. Validate every path released by removing this guard; do not assume the current full path is an available unfenced reference.

**Why (review M4):** during the uniform-grid sweep, IAN holds μ. Inflow BCs copy a stored conserved vector whose IAN is A (`src/bvals/physics/hydro_bcs.cpp:59-62`), and `PrimToCons`-based user BCs also write A.

**Change:**
- After `ApplyPhysicalBCs` in the two-var parabolic task list, convert A→μ **only on ghost slabs of faces whose BC is inflow or user**. Copy-type BCs (outflow, reflect, diode) already carry μ. Periodic faces are handled by the exchange.
- Alternatively, give BC functions a representation flag. Pick whichever is cleaner, and document it.
- Remove the T-E1 fatal error.

**Acceptance:**
- A 1D CGL-LF run with `ix1_bc = inflow` (isotropic inflow state) matches a representation-correct full-path reference to roundoff. Establish that reference explicitly because the baseline currently rejects it; also check an analytic uniform-state case. Cover explicit LF and multilevel paths before lifting their restrictions.
- Ghost p⊥/p∥ equals the inflow state.
- Add an equivalent test with a `PrimToCons`-based user BC.

## 6. Verify existing CGL-aware primitive prolongation on the GPU

**Triage: already implemented and tested on CPU in WO1 T-B6.** `src/bvals/prolong_prims.cpp` uses CGL projection helpers for both pressures and explicitly selects the anisotropy (A) or magnetic-moment (μ) slot. No WO1 primitive-prolongation fence was added. The old assumption that prolongation always happens outside STS is obsolete.

**Work:** reuse the existing smooth/uniform, strong-anisotropy, primitive-versus-conserved stage-order, churn, and restart tests on the GPU. Preserve representation-aware coarse/fine buffers during Task 3; do not replace the current implementation with the old proposed A-only branches. Add code only for a demonstrated gap.

**Acceptance:** existing tests pass on the actual GPU and, where available, MPI+GPU. Verify fine-ghost A/μ consistency in the representation active at the tested stage, pressure admissibility, repair accounting, and conservation. For smooth SMR data, retain the original O(Δx²) agreement check against `prolong_primitives = false` and report the measured error/convergence. Record this task as already done if no new implementation is needed.

## 7. Fuse the three direction flux kernels (performance, roundoff-level)

**Why:** after WO1 T-D1 and T-D2 settle the kernel contents, the x1, x2 and x3 LF flux kernels (`src/diffusion/cgl_landau_fluid.cpp`) each re-stream the same cell arrays. A fused kernel computing all face fluxes per cell (or per pencil) should roughly halve flux-kernel traffic. Hoist only frozen geometry such as b̂ and B̄. Recompute temperature- and limiter-dependent χ when its inputs change; profile the current GPU kernels before claiming a traffic or timing reduction.

**Acceptance:** differences are roundoff-level (report the magnitude), and speedups are reported for 2D and 3D.

## 8. (Optional) Entropy-consistent LF face flux (review M10)

**Why:** the continuous LF subsystem has an exact H-theorem:
σ = χ⊥ρG²/T⊥² + χ∥ρ(∇∥T∥)²/(2T∥²) ≥ 0, with G = ∇∥T⊥ − T⊥(1 − T⊥/T∥)∇∥B/B.

The grad-B term is exactly what makes the q⊥ part a perfect square. The discrete face flux can violate this at grid scale.

**Change:** write q⊥ at the face in entropy variables:
q⊥,f = χ⊥ρ_f T⊥,L T⊥,R b̂ [ (1/T⊥,R − 1/T⊥,L)/Δx + ū (|B|_R − |B|_L)/(Δx B̄) ]
with ū = mean(1/T⊥ − 1/T∥). Use B̄ in the μ flux, and do the same for q∥. The linearization stays unchanged at second order, and the cap preserves the sign.

**Acceptance:**
- In 1D, discrete entropy production is ≥ 0 face by face on randomized states.
- Convergence order and eigenmode accuracy are unchanged.
- The WO1 hot-spot and reversal tests still pass.

## 9. (Optional) Backlog of review low/info items
Pick these up as capacity allows. These are historical review leads: verify the current behavior before implementing each item. Each selected item needs its own small test.
- **Exact kinetic mirror threshold.** Offer the p⊥/p∥-weighted form (p⊥/p∥ − 1 > 1/β⊥) as an option.
- **Smooth ν_eff switching.** Avoid discontinuous ν_eff jumps at thresholds, e.g. with a smooth ramp, and evaluate it once per sweep instead of per stage.
- **Symmetric collision split** (optional mode) for second-order collisional runs. Limiters are currently Lie-split from the Strang composition. Post-AMR wall projection already exists and must be preserved. Check initial-state wall handling separately; do not assume it was implemented with the post-AMR fix.
- **`sts_max_dt_ratio` default.** Consider a finite default. It is currently unbounded, and RKL2 damps stiff modes only weakly at very large s.
- **Single-precision robustness.** Log-space A encoding and overflow-safe heat-flux cap arithmetic already exist; verify their coverage rather than duplicating them. Full application float compilation has separate pre-existing EOS/table/coordinate/units/diffusion portability failures recorded in the WO1 report. Triage those only if this optional item is selected, and distinguish focused float checks from a working full float application.
- **The `hfpow` history diagnostic** in `cgl_lf_paper` is q·∇T, not a heat-flux power, and uses a different stencil from the evolved fluxes. Rename it or recompute it from the actual face fluxes.
- **`kh.cpp`** sets CGL primitives with ideal-gas semantics. Fix it, or fence it for CGL.
- **Explicit viscosity or resistivity with CGL** is currently allowed silently. Fence it, or document the support.

---

## Final deliverables
- A branch from the merged WO1 baseline with one commit per executed task, all applicable tests passing, and the two accepted B4 strict expected failures retained. Record already-completed and unselected optional tasks explicitly.
- A final report with each task's triage, measured result change, acceptance/convergence evidence, and performance change on the actual GPU. Include device/toolchain/allocation metadata, ranks, stage counts, timing definitions/repeats, skipped coverage, and retained input/output hashes.
- An updated `docs/cgl_lf_changes.md`, or equivalent, listing the new parameters (`sts_safety`, `sts_merge_half_sweeps`) and behaviour changes.
- Docs and README updates for removed fences and the new passive-mode semantics.
- A durable GPU baseline and comparison archive, including the 0-P1 diagnosis/regression and CPU/GPU/MPI distinctions. Do not present compilation, skipped tests, smoke tests, or accepted B4 failures as broader physics validation.
