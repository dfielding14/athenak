# AthenaK MHD-PIC Program Guide

> Tracked working guide for the PIC_development branch, revised 2026-07-20.
> Part I is the active uniform-grid implementation and science program.
> Part II is the deferred static/adaptive mesh-refinement roadmap.
> Update this guide when evidence changes a decision; do not turn it into a
> speculative requirements catalog.

## 1. Objective, priorities, and project standard

The primary scientific objective is to study kinetic cosmic rays in
three-dimensional turbulent MHD boxes. The shock calculations have a narrower
purpose: they validate the MHD-PIC method in a real, self-consistent Bell
precursor and provide a clear proof of concept for the methods section of the
eventual paper.

The program is divided into two parts.

### Part I — Uniform-grid implementation and science

Part I is active:

1. preserve the physically complete full-Hall MHD-PIC closure;
2. make its supported uniform-grid path fast and operationally reliable;
3. run one efficiently designed 3D Bell shock, plus at most one
   evidence-driven numerical repeat;
4. evolve uniform-grid MHD turbulence to a declared saturated state;
5. introduce CR particles into that common saturated state; and
6. run the production uniform-grid CR-turbulence boxes.

The goal is a rock-solid supported workflow, not a promise that every AthenaK
configuration works with MHD-PIC. Rock solid means that the selected model is
physically coherent, conserves the intended quantities, restarts cleanly,
retains particles, completes MPI communication, writes manageable output, and
runs efficiently in production-shaped jobs.

### Part II — Static and adaptive mesh refinement

Part II is a detailed future plan and does not gate Part I. It covers three
different capabilities:

1. **time-triggered whole-domain refinement**, in which every block is refined
   after a chosen time;
2. **static mesh refinement**, especially fixed fine regions covering a shock
   and its precursor; and
3. **adaptive mesh refinement**, in which those regions move or change size.

The first item is exactly the global resolution staircase discussed in project
planning. It uses AMR machinery for a transition such as
\(128^3\rightarrow256^3\), but once every block has refined, the leaf mesh is
spatially uniform again. Each settled epoch is therefore scientifically a
uniform-grid calculation. It is distinct from persistent spatial AMR.

### Working rules

- Physical signs, units, normalization, time centering, staggering, momentum
  exchange, energy exchange, and CT must describe one consistent equation set.
- Make the smallest coherent change that solves a demonstrated problem.
- Match AthenaK style and reuse its task, data, communication, restart, and
  output paths.
- Profile production-shaped work before optimizing it.
- Repeat the workload that exposed a bottleneck; do not create a general
  performance matrix.
- Use one decisive physics calculation and one uncertainty-triggered repeat.
- Run long simulations and substantial analysis as self-contained worker jobs.
- Preserve immutable source, input, executable, and analysis provenance.
- Fix observed failures and direct physical inconsistencies. Do not pre-build
  defenses against every imaginable state.

The fuller rationale is in the
[MHD-PIC project ethos](docs/source/engineering/pic_project_ethos.md).

## 2. Current validated baseline

The source repository is
[/ccs/home/dfielding/athenak-pic](/ccs/home/dfielding/athenak-pic) on branch
PIC_development. The validated reusable uniform-grid reference is tagged
**pic-science-candidate-20260713-v2**. Its principal component commits are:

- **978344bcd**: conservative CR-Hall coupling;
- **b808d3400**: narrow correctness cleanups;
- **683e120ca**: analysis and shock workflow;
- **de65bc447**: reusable qualified Frontier build template;
- **b5422ca57**: scalable uniform shock-source preparation and profiling;
- **6600c0806**: particle-migration completion; and
- **225038f76**: conditional same-level redistribution for localized uniform
  shocks.

That tag is immutable. The previous **pic-science-candidate-20260713** tag at
**7b784e3b8** remains the pre-performance reference.

The completed Mignone-R2-style 2D3V shock is frozen separately at
**pic-nonrelshock-t3000-20260720** (**3e03007e6**). That tag contains the exact
safe-carrier injection source, production and smoke decks, focused tests, and
setup specification used by the archived \(t=3000\) campaign. The exact
Frontier executable, input, source patch, output, final restart, figures, and
hash manifests are retained under
[`NonRelShock`](/lustre/orion/ast207/proj-shared/dfielding/PIC/NonRelShock).
This problem-specific tag is immutable evidence for the shock result; it does
not replace the reusable v2 baseline for unrelated workflows.

| Work item | Status | Decisive evidence or next action |
| --- | --- | --- |
| Paper-to-code full-Hall map | Complete | [CR-Hall paper-to-code map](docs/source/engineering/pic_cr_hall_code_map.md) |
| Full-Hall linear physics | Complete | Job 4975948 passed all 12 Hall-off, full-Hall, moment, wavelength, frequency, polarization, and resolution checks |
| Exact source/core qualification | Complete | Job 4978017 passed conservation, source closure, normal precision, task ordering, and safety checks |
| MPI single-precision smoke | Complete | Two-rank GPU job 4975966 passed |
| Nonlinear Bell onset | Bounded pass | Job 4975940 reached \(B_{\perp,\mathrm{rms}}/B_0=1.043\); no saturation claim |
| Compact controlled 2D shock | Complete | Jobs 4976035 and 4977996 plus analysis 4978079 established a stable shock and coherent Bell precursor |
| Mignone-R2-style 2D3V shock | Complete and frozen | The full-Hall uniform run completed continuously from \(t=0\) to 3000; see the [setup specification](docs/source/engineering/pic_mignone_r2_shock_setup.md) and frozen campaign archive |
| Matched 128-cubed boxes | Engineering comparison complete | Full-Hall 4976062 and Hall-off 4976612 completed through \(t=10\) and exposed a measurable Hall response |
| 3D shock scaling probe | Complete | Job 4980740 passed and exposed shock-source duplication and early particle imbalance |
| Uniform-grid performance | PERF-U1 through PERF-U4 complete | Source preparation and migration improved; unsafe Frontier atomics cut representative-box cost by 13%; conditional shock redistribution delivered a further 1.40x speedup; stop general optimization |
| Focused correctness cleanups | Complete | Invalid macro weights are rejected, legacy nonphysical Bell mechanics require explicit opt-in, and focused test discovery is repaired |
| Production 3D shock | Next; candidate frozen | Finalize the best-science-per-node-hour uniform case; the old narrow deck remains a benchmark, not the automatic science target |
| Production turbulent boxes | Primary science after 3D proof | Saturate MHD first, save one common checkpoint, then introduce CRs |
| SMR/AMR | Deferred Part II | Full Hall remains fail-closed on multilevel meshes |

This guide is now the single durable program record. Temporary execution
trackers should not be recreated unless a new active campaign genuinely needs
one.

### 2.1 Frozen Mignone-R2-style shock through \(t=3000\)

The current method-section shock reference is the uniform 2D3V, full-Hall
Mignone-R2-style campaign documented in the
[shock setup specification](docs/source/engineering/pic_mignone_r2_shock_setup.md).
It used a \(120000\times3000\,c/\omega_{pi}\) domain on a
\(11520\times288\) mesh, a parallel upstream field, Alfvénic Mach number 30,
and shock-created CRs. The injection ledger uses the swept-mass trigger, a
cheap \(p>30p_0\) carrier prefilter, and an exact carrier-local affordability
transaction. An unaffordable carrier remains pending; the run does not clip
the subtraction, redistribute it to another cell, or invoke the legacy global
floor throttle.

The archived history is strictly monotonic from \(t=0\) through 3000 with 3001
rows. It contains 101 fluid/moment snapshots, 51 full particle snapshots, 31
publication triptychs, and one verified final restart. The Release,
double-precision MPI+HIP executable used `-munsafe-fp-atomics`; its SHA-256 and
the exact deck, environment, and source-patch hashes are recorded in the
archive validation manifest. This establishes a long, restartable,
production-scale uniform-grid shock workflow and preserves the magnetic
morphology and CR-distribution evidence needed for the method comparison. It
does not by itself establish a 3D result or a converged maximum-energy spectrum.

### 2.2 Earlier compact Bell qualification

The qualified long shock reached \(t=400\) on a \(5000\) by \(128\) mesh. It
measured:

- shock speed \(9.941\) and terminal compression \(4.009\);
- a 7836-code-unit current-bearing precursor;
- integrated Bell growth opportunity \(N_{\mathrm{Bell}}=6.64\);
- terminal upstream \(B_{\perp,\mathrm{rms}}/B_0=0.185\);
- near-shock transverse field reaching \(0.673B_0\);
- dominant wavelength \(488.8\), or 162.9 cells at \(dx=3\);
- helicity \(-0.998\) at the spectral peak;
- magnetized filaments and underdense cavities; and
- emerging acceleration with \(p_{99}/p_{\mathrm{inj}}=1.78\) and
  \(p_{\max}/p_{\mathrm{inj}}=2.45\).

The terminal mean and local-maximum parallel Hall parameters were only 0.059
and 0.203, and the mean predicted Hall wavelength shift was below 0.1 percent.
A matched Hall-off 3D shock is therefore not part of the validation plan.

These compact measurements remain the quantitative Bell scale, handedness,
and Hall-strength benchmark used to design the 3D proof run. The frozen
\(t=3000\) Mignone-R2-style campaign supersedes it as the long-duration 2D
method reference.

### 2.3 Completed 3D throughput probe

Production-shaped job 4980740 ran the current
[3D benchmark deck](/ccs/home/dfielding/athenak-pic/inputs/publication/pic_parallel_shock_cr_hall_3d_full.athinput)
through \(t=47\) on 32 Frontier nodes:

- mesh \(5000\times128\times128=81.92\) million cells;
- 1808 cycles and 115.82 simulation seconds;
- pre-injection cost 0.05766 seconds per cycle;
- early post-injection cost 0.14010 seconds per cycle;
- 4.349 million particles at completion;
- 58.42 GB of tracked AthenaK allocations globally;
- exact mesh balance and 3.35 percent particle-load efficiency; and
- no memory, particle-validity, or physical failure.

Every rank currently gathers and sorts the complete \(128\times128\) shock
surface, constructs 507904 stencil keys, and repeats every injected tag's
kinematics. That duplicated source preparation accounts for about 0.042
seconds per cycle in the early post-injection fit.

Projection from the completed 2D history gives approximately 563 million
particles at \(t=350\) for the old narrow deck. Its present planning cost is
roughly 112 to 128 node-hours including I/O, with about 245 GB of retained
output. This proves feasibility and provides a benchmark. It does not make the
384-wide transverse domain, only about 0.8 of the measured nonlinear Bell
wavelength, the correct science design.

## 3. Physical contract that Part I must preserve

### 3.1 Model and references

[Bai et al. (2015)](https://arxiv.org/abs/1412.1087) is the primary equation
source. [Mignone et al. (2018)](https://arxiv.org/abs/1804.01946) provides an
independent conservative implementation cross-check.
[Sun and Bai (2023)](https://arxiv.org/abs/2304.10568) supplies the
Athena++-style VL2/Boris/TSC pattern and refinement/performance context, while
intentionally omitting CR Hall from induction.

The supported Part I model is non-relativistic ideal MHD coupled to relativistic,
full-orbit CR superparticles on a uniform Cartesian mesh. It assumes:

- massless thermal electrons enforcing charge neutrality;
- gas mass dominated by thermal ions;
- dynamically negligible CR mass density;
- scales much larger than the thermal-ion inertial length; and
- full-f particles using the paper MHD-PIC VL2/TSC integrator.

Conventional Hall MHD, electron pressure and inertia, resistivity, delta-f,
expanding boxes, relativistic MHD, and multilevel meshes are outside this
supported full-Hall path.

### 3.2 AthenaK code-unit equations

Use the signed deposited moments

\[
Q_{\rm cr}=q_{\rm cr}/c,\qquad
\boldsymbol K_{\rm cr}=\boldsymbol J_{\rm cr}/c,
\]

and define

\[
Q_e=\alpha_i\rho+Q_{\rm cr},\qquad
\boldsymbol D=\boldsymbol K_{\rm cr}-Q_{\rm cr}\boldsymbol u_g,\qquad
\boldsymbol v_H=\boldsymbol D/Q_e.
\]

The stored electric field is \(c\boldsymbol E\):

\[
\boldsymbol{cE}_0=-\boldsymbol u_g\times\boldsymbol B,\qquad
\boldsymbol{cE}_H=-\boldsymbol v_H\times\boldsymbol B,\qquad
\boldsymbol{cE}=\boldsymbol{cE}_0+\boldsymbol{cE}_H.
\]

The particle force and power are

\[
\boldsymbol F_{\rm cr}
=Q_{\rm cr}\boldsymbol{cE}+\boldsymbol K_{\rm cr}\times\boldsymbol B,
\qquad
P_{\rm cr}=\boldsymbol K_{\rm cr}\cdot\boldsymbol{cE}.
\]

The gas receives the opposite particle momentum and kinetic-energy changes.
Its total energy also receives the conservative Hall flux

\[
\boldsymbol F_{E,H}=\boldsymbol{cE}_H\times\boldsymbol B.
\]

AthenaK has absorbed the usual factor of \(4\pi\). No extra physical light
speed, artificial particle light speed, or adjustable Hall-strength
coefficient belongs in these grid equations.

The model selection is atomic:

- **pic_cr_hall_mode=full** enables the full pusher, CT, feedback, and Hall
  energy-flux closure;
- **pic_cr_hall_mode=off** removes the complete Hall correction while retaining
  gas-particle momentum and energy exchange.

The physical mode uses the deposited difference
\(\boldsymbol K_{\rm cr}-Q_{\rm cr}\boldsymbol u_g\) directly. It does not
divide by noisy \(Q_{\rm cr}\). The background-ion charge-to-mass parameter is
a physical normalization, not a tunable Hall coefficient.

### 3.3 Time centering and exchange

The full-f particle chain runs after the MHD conserved-state copy and before
MHD fluxes in both VL2 stages:

1. Stage 1 deposits base-time signed charge and current. A scratch half kick
   predicts midpoint particle momentum without changing the true state or
   returning a gas impulse.
2. The stage-1 grid predictor uses the predicted Hall drift for induction and
   Hall energy flux. Analytic particle feedback advances the gas with the
   opposite sign.
3. Stage 2 deposits midpoint charge and current from midpoint position and
   scratch-predicted momentum. That field time-centers the true Boris push.
4. The true push deposits its realized momentum and kinetic-energy rates,
   DPDT and DEDT. The gas receives their opposites exactly once.
5. The stage-2 grid corrector derives the Hall electric field from the exact
   discrete particle impulse rather than adding another analytic force.
6. Hall induction and \((\boldsymbol{cE}_H\times\boldsymbol B)_n\) enter the
   face fluxes before FOFC and RK. CornerE constructs the edge field through
   the normal CT path.
7. FOFC tests the composite flux-plus-particle-source state. If it replaces a
   face flux, Hall induction and Hall energy are restored together from the
   same donor; the deposited particle impulse is not clipped or redistributed.

Predicted Hall drift belongs to particle time centering and the stage-1
predictor. Deposited DPDT closes the stage-2 grid update. Do not combine or
double count them.

### 3.4 Artificial light speed and diagnostics

The artificial particle light speed controls only relativistic particle
kinematics. It does not divide charge-to-mass ratio, deposited current,
background-ion normalization, Hall drift, or Lorentz force.

The two global regime diagnostics are

\[
\max |R|=\max\left|Q_{\rm cr}/Q_e\right|,\qquad
\max\Lambda=\max |\boldsymbol v_H|/v_A.
\]

They characterize the simulated regime. They are not automatic warnings or
abort criteria.

### 3.5 Injection and CR initialization must be physically explicit

Every production problem must state whether a newly introduced CR population:

- represents external CRs added to the modeled system; or
- receives mass, momentum, and energy removed from the gas.

The shock uses a conservative swept-mass/injection ledger and returns the
opposite gas changes. A turbulent-box CR population may represent an externally
specified population, but that choice must be explicit in the deck, analysis,
and paper. Performance work must never alter these semantics.

## 4. Part I — Uniform-grid implementation and science program

### 4.1 Part I completion criterion

Part I is complete when:

- the measured uniform-grid bottlenecks are fixed or demonstrated adequate;
- one new immutable descendant candidate passes compact correctness and
  production-reliability gates;
- one interpretable 3D shock validates the dimensional method claim;
- the saturated-MHD-to-CR box workflow is physically explicit and reproducible;
- at least one retained uniform-grid CR-turbulence science box is complete; and
- costs, figures, limitations, and provenance are documented.

SMR or AMR does not block this completion rule or promotion of the uniform model.

### 4.2 PERF-U1 — Shock-source preparation

The first performance change targets the duplicated source preparation exposed
by job 4980740.

Implementation objective:

- avoid allgathering and sorting the complete surface on every rank every cycle;
- cache geometry that is invariant between topology changes;
- construct and apply only locally owned source work, using compact reductions
  or routing for the exact global ledger;
- retain equal-area surface averaging and the declared downstream stencil; and
- preserve deterministic injected-particle and gas-source semantics.

Validation is deliberately compact:

1. run the closest focused injection/source regression;
2. run the unchanged normal-precision core check if the coupling path changed;
3. repeat the identical 32-node \(t=47\) benchmark; and
4. compare injected mass, reservoir, momentum, energy, shock trajectory,
   particle validity, and throughput.

Keep the change only if physics is unchanged within the appropriate reduction
tolerance and throughput improves materially. Do not create a permanent
performance acceptance bureaucracy.

**Completed result.** Focused job 4981305 passed exact serial/MPI injection
parity, restart controls, and the normal-precision Hall core. The identical
32-node repeat 4981353 reduced post-injection cost from 0.14010 to 0.12917
seconds per cycle (7.8 percent) and the incremental injection-window overhead
by 13.0 percent. Counts were identical, invalid records remained zero, and the
maximum scaled history difference was `5.6e-14`.

### 4.3 PERF-U2 — Particle-migration completion

The full-Hall \(128^3\) box spent 21.4 percent of total runtime in particle
migration and made roughly \(8.95\times10^8\) migration-wrapper calls. Mesh
and particle populations were already balanced, and output was below one
percent of runtime.

Instrument only enough to confirm the completion path, then replace unnecessary
MPI polling with the smallest existing-style completion mechanism. Benchmark one
fixed-cycle production-shaped box and check:

- particle count, identity, and ownership;
- MPI completion without hangs;
- unchanged gas-particle conservation and Hall observables; and
- end-to-end rather than kernel-only improvement.

Do not rewrite output, Hall kernels, deposition, or the box decomposition from
this profile.

**Completed result.** The receive-completion task now uses one blocking
completion at the terminal point of its strictly ordered migration chain rather
than repeatedly returning incomplete to the task scheduler. Focused migration,
restart-ledger, and telemetry checks passed in job 4981735. Two sequential
treatment/control benchmarks ran from identical full-Hall `128^3` checkpoints
on the same node: job 4981821 at `t=5` and job 4981838 at `t=10`, with 300
cycles per executable. Every retained history row was bitwise identical, all
12,582,912 particles remained present, and invalid records remained zero.

The late-state comparison reduced receive-wrapper calls by 96.6 percent,
fitted cycle cost by 3.1 percent, median interval cost by 3.4 percent, and total
driver time by 5.4 percent. The migration section itself improved by only 1.6
percent. The original uninterrupted run's extreme polling count was not
reproduced after restart, so this is recorded as a modest end-to-end
improvement, not recovery of the full 21.4 percent migration fraction. Keep the
small change and move on; reopen migration only if a target-experiment profile
again shows a material cost.

### 4.4 PERF-U3 — Conditional same-level particle redistribution

The post-injection shock profile demonstrated particle localization as a real
limiter. AthenaK's existing nonadaptive same-level redistribution machinery was
therefore extended to uniform meshes without enabling multilevel physics:

- use actual per-block particle counts and a measured particle cost;
- redistribute only at a sparse cadence;
- include migration cost in the end-to-end comparison; and
- remove or disable it if the net gain is negligible.

This is smaller and more directly useful than beginning AMR solely to obtain
load balancing.

**Completed result.** The 32-node, 100-cycle profile `4982464` measured `3.93%`
particle-load efficiency. Critical-rank push plus deposition occupied `26.1%`
of driver time, while particle-array resizing occupied only `0.55%`. Matched
job `4982553` applied one redistribution at the beginning of the identical
window. Driver time fell from 16.248 to 11.583 seconds (`1.40x` speedup), the
maximum particle count per rank fell by `37.9%`, and critical-rank push and
deposition each fell by about `51%`. All `8,098,700` particles remained,
invalid records were zero, and histories agreed to `5.6e-14` maximum scaled
difference.

The tracked particle-storage high-water mark rose by `34.3%`, to `0.575` GiB
per rank. Retain this as an opt-in, sparsely applied tool for localized uniform
shocks and remeasure memory for the final geometry. Leave it off for balanced
turbulent boxes. The full-Hall multilevel prohibition remains unchanged; this
work does not begin Part II.

### 4.5 PERF-U4 — Target-experiment profiles and follow-on optimization

The two target profiles cover the workloads we actually intend to run:

1. **Localized shock workload:** use the optimized final-geometry probe just
   after injection and one particle-rich continuation milestone. Separate MHD,
   shock-source preparation, particle gather/push/deposition, migration,
   particle imbalance, global reductions, and I/O.
2. **Volume-filling turbulent-box workload:** use a short fixed-cycle restart
   from saturated MHD after CR introduction at the intended particle density.
   Separate MHD/Hall work, gather/push/deposition, particle sorting/locality,
   migration, host-device movement, output, and load balance.

Those profiles decide, rather than automatically require, the next improvements:

- tune particle sorting cadence or data locality for a volume-filling box;
- use contiguous/intermediate TSC storage only if measured particle density and
  memory cost make it profitable;
- batch or fuse particle kernels only if launch overhead is material;
- remove repeated host staging or global reductions where they dominate;
- use same-level particle redistribution for the localized shock; and
- choose MeshBlock, rank, and node layout from measured end-to-end scaling.

Change and benchmark one demonstrated limiter at a time. Retain an optimization
only when the production-shaped workload improves without changing physics or
creating unreasonable memory growth.

**First target-box result.** Synchronized 100-cycle profile 4981885 used the
full-Hall `128^3`, 12,582,912-particle `t=10` checkpoint. Particle balance was
99.8 percent and invalid records were zero. Deposition consumed 38.9 percent of
driver time, push/gather 16.5 percent, and migration 13.7 percent, making the
particle path about 69 percent of this workload. This rules out box load
balancing as the next optimization target and identifies deposition/data
locality as the main remaining opportunity.

One deliberately small deposition experiment skipped the redundant shape
renormalization pass when the uniform-grid TSC stencil was complete. Focused
normal-precision conservation and two-GPU single-precision turbulence checks
passed in jobs 4981910 and 4981942. The same-checkpoint production A/B 4981964
retained bitwise-identical histories and every particle, but fitted and median
cycle costs worsened by 2.3 and 1.8 percent. The experiment was reverted. Do not
revive it without different evidence; any next deposition change must address
atomic contention or particle/cell locality and earn an end-to-end improvement.

The next atomic-contention test required no source rewrite. Job 4982226 rebuilt
the same snapshot for Frontier/gfx90a with `-munsafe-fp-atomics` and compared it
with the safe build for 300 cycles from the identical `t=10` checkpoint on the
same node. Fitted cycle cost, median cycle cost, and total driver time each fell
by about 13.3 percent; stage time fell by 15.2 percent. Both runs reached cycle
17317, retained all 12,582,912 particles with zero invalid records, and had
identical histories at written precision. Job 4982256 then passed the existing
HIP double-precision full-Hall conservation test and two-GPU single-precision
turbulence smoke. Use the flag in the qualified Frontier ROCm 6.2.4/gfx90a
production build profile. Because its floating-point and memory assumptions are
platform-specific, repeat this compact qualification after a ROCm or GPU
architecture change rather than treating it as a portable source default.

**Target-shock result.** Job `4982464` profiled 100 fixed post-injection cycles
on 32 nodes with `8,098,700` particles and zero invalid records. Source
preparation consumed `18.6%` of driver time, push plus deposition `26.1%`, and
migration `8.9%`. The inclusive append operation was `1.11%`, including only
`0.55%` in particle-array resizing. Although resizing copied an estimated 128
GB summed across ranks over the window, it was not a wall-time limiter, so
capacity-managed particle growth was rejected. Localization was material, and
the single redistribution experiment documented in PERF-U3 earned a `1.40x`
end-to-end speedup.

The two intended workload classes have now selected the useful changes. Stop
general performance optimization. Reopen sorting, data layout, kernels,
staging, migration, or capacity growth only if a final science-shaped run
demonstrates a material limiter.

### 4.6 Production reliability and focused requalification

Reliability work follows the supported workflow, not every possible
configuration.

After each coherent implementation batch:

1. run the focused normal-precision core/conservation checks;
2. run the MPI single-precision turbulence smoke when that path changed;
3. repeat the production-shaped workload that motivated the change; and
4. inspect its concise result, particle count, CT divergence, conservation
   ledger, memory high-water mark, and terminal state.

Before the long science runs, exercise one real checkpoint/restart continuation
using the intended decomposition and output path. It must preserve:

- particle identity, weight, species, momentum, birth information, and ownership;
- forcing and problem-generator state needed by the run;
- injection reservoirs and cumulative source ledgers;
- Hall mode and normalization fingerprints;
- gas-plus-particle momentum and energy accounting; and
- the next-step result within the intended numerical reproducibility.

Fix crashes, hangs, particle loss, invalid records, source overdraw, silent
physics changes, and material conservation or CT regressions. Do not add broad
platform, precision, timestep, and particle-number matrices without evidence.

### 4.7 Frozen uniform-science candidate

Never move or overwrite `pic-science-candidate-20260713` or its qualified
successor, `pic-science-candidate-20260713-v2`. The performance batch was split
into independent shock-source, migration, redistribution, build, and
documentation commits. The compact reliability gates passed, and the v2 tag:

1. records the before/after production timings and physics comparisons;
2. archives a reproducible Frontier build configuration including
   `-munsafe-fp-atomics`;
3. contains the focused correctness, restart, and MPI smoke coverage; and
4. is the exact candidate to use for the redesigned shock and boxes.

The reusable Frontier build entry point is
[`build_frontier.sh`](/ccs/home/dfielding/athenak-pic/build_frontier.sh), added
in commit `de65bc447`. Profiling build `4982429` qualified its Release,
double-precision, MPI/HIP, gfx90a configuration with
`-munsafe-fp-atomics`. Exact-final source job `4982681` used the same profile
and passed the normal-precision full-Hall conservation test plus the two-GPU MPI
single-precision turbulence smoke. Its source snapshot SHA-256 is
`472096c3b26ec574212bbb97dc8bf1c158045332fec78ae34bc26173127c223d`.

The later tag `pic-nonrelshock-t3000-20260720` freezes the exact demonstrated
safe-carrier shock source at commit `3e03007e6`. Its archived executable and
campaign products completed the uniform Mignone-R2-style run through
\(t=3000\). Preserve that tag for shock reproduction; keep the reusable v2 tag
as the pre-shock baseline and create a new candidate only after a coherent
cross-cutting implementation batch passes the focused requalification above.

### 4.8 SHOCK-U1 — Smart 3D method validation

#### Scientific role

The 3D shock is not intended to become a major parameter campaign. It must show
that the qualified 2D Bell scale and handedness, magnetic amplification,
cavities, filaments, shock corrugation, CR transport, and emerging acceleration
remain interpretable without the geometric restrictions of 2D.

[van Marle et al. (2019)](https://arxiv.org/abs/1909.06931) remains the closest
large 3D MHD-PIC comparison. Its principal lesson for this program is that a
moderately wide 3D domain can establish tubular morphology and CR escape without
requiring an unlimited grid.

#### Provisional best-science-per-node-hour design

PERF-U1 through PERF-U4 are measured and the qualified uniform candidate is
frozen. Choose the final deck from this working design:

- full Hall, delayed injection, and the same qualified physical shock model;
- \(dx\) near 6, resolving the injection gyroradius with about 16 cells;
- about 81 cells across the measured 2D dominant wavelength;
- longitudinal extent set by the measured precursor, roughly 160 to 190
  injection gyroradii rather than copied blindly from the old deck;
- transverse width roughly 12 to 16 injection gyroradii, sufficient for several
  Bell structures and nonlinear tubes/cavities;
- 8 to 12 independently sampled downstream particles per cell unless a short
  current-noise measurement requires more;
- duration roughly \(t=700\) to 1000, or the first time at which the declared
  morphology, acceleration, and escape measurements are decisive; and
- sparse field, moment, tracked-particle, checkpoint, and final-particle output
  designed around the analysis rather than convenience.

This corresponds approximately to 90 to 200 million cells depending on the
selected length and width. The post-performance planning range is 500 to 1500
Frontier node-hours for the fiducial. The entire fiducial-plus-repeat program
should remain inside a 3000-node-hour envelope.

#### Execution gate

1. Run one short production-shaped cost and current-noise probe with the final
   geometry, sampling, executable, and output definitions.
2. Project runtime, memory, particle count, restart size, and total retained
   output from that probe.
3. Run the fiducial if it fits the approved budget.
4. Add one resolution or particle-statistics repeat only if the result makes
   that uncertainty limiting.

Do not add a Hall-off 3D shock unless the measured upstream \(\Lambda\) becomes
large enough to make it a physical comparison. Do not begin a Mach-number,
obliquity, injection-efficiency, or shock-parameter grid as an automatic sequel.

#### Required measurements and claim boundary

Measure:

- shock trajectory, compression, stability, and corrugation;
- current-bearing precursor width and integrated Bell growth opportunity;
- measured versus locally predicted Bell scale;
- magnetic amplification and spectrum;
- handedness and polarization;
- 3D tubes, cavities, and filaments;
- CR residence, escape, spectrum, acceleration efficiency, and maximum momentum;
- gas-CR momentum and energy exchange; and
- memory, throughput, particle balance, and output cost.

The intended claim is method validation and physically interpretable 3D proof
of concept. It is not a converged maximum-energy spectrum or a definitive
state-of-the-art shock survey.

### 4.9 BOX-U1 — Saturated-MHD-to-CR turbulent boxes

The turbulent boxes are the principal science program.

#### Phase A: MHD dynamo

For each retained resolution:

1. run forced MHD without kinetic CRs on the final uniform science mesh;
2. retain the forcing realization, OU amplitudes/state, RNG state, physical
   forcing band, and normalization;
3. declare saturation or statistical stationarity using a compact sustained
   interval of magnetic and kinetic energies, RMS quantities, and spectra;
4. continue long enough that the selected statistics are not a transient
   crossing; and
5. save one provenance-bound saturated MHD checkpoint as the common initial
   condition.

Part I does not use mesh refinement during this run. If Part II later uses a
time-triggered global resolution staircase, saturation and the retained
checkpoint must still be established after the final-resolution cascade has
settled.

#### Phase B: CR introduction

Introduce CRs only from the common saturated MHD state. Before running, state:

- distribution, species, charge sign, artificial light speed, and physical
  normalization;
- particle number or target particles per final-resolution cell;
- whether CR energy and momentum are externally supplied or removed from gas;
- full-Hall background-ion normalization;
- any tracked-particle selection; and
- the post-injection adjustment interval excluded from stationary statistics.

Run one compact particle-loaded restart smoke that checks initialization
accounting, CT, total momentum and energy behavior, particle retention, Hall
diagnostics, output, and restart continuity. Then profile the particle-rich box
before selecting the final node layout.

#### Phase C: science runs

Run the full-Hall fiducial from the shared saturated checkpoint. Add a Hall-off
control only when the measured \(\Lambda\) distribution or a specific paper
claim makes it informative. Add physical comparisons one at a time.

For each retained run measure:

- kinetic and magnetic energies, RMS quantities, spectra, and intermittency;
- CR energy, spectrum, anisotropy, transport, acceleration, and tracked
  trajectories where useful;
- gas-CR momentum and energy transfer;
- \(R\) and \(\Lambda\) distributions rather than maxima alone when science
  interpretation requires them;
- correlations of CR current/charge with magnetic and velocity structure; and
- conservation, CT divergence, runtime, memory, and particle balance.

Increase uniform resolution deliberately. Each resolution is a separate
production run from an appropriate saturated MHD checkpoint; do not treat a
short resolution change as stationary science.

### 4.10 Analysis, publication products, and promotion

Long simulations perform compact analysis on worker nodes and retain only the
raw outputs required to reproduce the selected PNG figures and tables. For
each science result preserve:

- immutable source and executable identity;
- exact input and runtime overrides;
- checkpoint ancestry;
- concise JSON measurements;
- publication-quality PNG figures;
- analyzer identity; and
- the one numerical repeat required by the claim.

After final code review and uniform-grid science validation, prepare
PIC_development for merge into PIC with a clean history, concise physical
description, reproducible evidence, measured performance, and explicit
multilevel limitation. Do not merge without user approval. Part II refinement
does not have to be complete before promoting a correct uniform implementation.

## 5. Worker-job and evidence policy

Every long simulation or substantial analysis is a self-contained worker job.
A useful worker job:

1. creates a run directory from the scheduler job ID;
2. verifies an immutable source snapshot and records environment, source HEAD,
   input, launcher, executable, and analyzer hashes;
3. builds on node-local storage when appropriate;
4. runs without an interactive watcher;
5. writes bounded outputs and checkpoints;
6. performs the declared analysis on a worker node;
7. writes a concise JSON result and sharp PNG figures; and
8. creates PASS only after simulation and required analysis succeed.

Do not use timestamps to associate outputs. Use the explicit run directory and
known basename. Occasional targeted scheduler checks are fine; do not sit in a
polling loop. Partial outputs are diagnostic evidence, not a pass.

For development decisions retain source, build, input, terminal state, concise
measurements, and the figures supporting the decision. Science results also
retain the raw data necessary to reproduce their figures and the targeted
repeat needed by the claim.

## 6. Part II — Static and adaptive mesh refinement

> **Status: deferred.** Part II begins after the uniform program, or earlier
> only if a production cost model shows that refinement is required for
> memory/output capacity or should improve end-to-end runtime materially.

### 6.1 Capability order

| Capability | Purpose | PIC active during refinement? | Current status |
| --- | --- | ---: | --- |
| Time-triggered whole-domain refinement | Cheap MHD dynamo spin-up before CR introduction | No | Existing MHD mechanism; uniform PIC handoff unqualified |
| Static mesh refinement | Keep a fixed shock and its Bell precursor highly resolved | Yes | Full-Hall coupling unsupported |
| Adaptive mesh refinement | Move or resize the shock/precursor region | Yes | Full-Hall coupling unsupported |

Implement in this order. Whole-domain MHD refinement is the simplest special
case. Static refinement establishes the multilevel PIC interface physics.
Dynamic AMR must reuse that operator and add only repeated topology change,
tagging, ownership, and balancing.

### 6.2 Time-triggered whole-domain refinement

The user's normal time criterion is exactly the global resolution staircase.
AthenaK already contains this pattern in
[turb_timed_amr.cpp](/ccs/home/dfielding/athenak-pic/src/pgen/turb_timed_amr.cpp):
after **t_refine**, every block is flagged.

The checked-in two-level example performs one global transition. A deliberate
\(128^3\rightarrow256^3\rightarrow512^3\) sequence needs a separate trigger
for each level or staged jobs; one threshold with additional allowed levels
would continue flagging blocks at later eligible checks until the maximum.

The preferred turbulence use is:

1. evolve the dynamo without CR particles;
2. trigger one global factor-two refinement at each selected time;
3. preserve the physical forcing band, OU amplitudes/state, phases, RNG state,
   and correlation time;
4. allow newly available high-wavenumber modes to fill after every transition;
5. establish saturation and retained statistics at the final resolution;
6. produce a qualified uniform-grid state handoff; and
7. introduce CRs only in the genuinely uniform Part I PIC run.

The handoff must preserve conserved MHD state, face-centered magnetic flux,
\(\nabla\cdot\boldsymbol B\), forcing state, and physical time. Do not merely
remove the multilevel flag from an AMR restart. Current full Hall rejects any
multilevel Mesh, even when all leaves happen to occupy the same final level.
Qualify either a supported uniform-restart conversion or an explicit
full-resolution state import.

A single small \(64^3\rightarrow128^3\) MHD comparison is sufficient initially.
Check volume-integrated conserved quantities, magnetic flux and divergence,
forcing continuity, low-wavenumber spectral continuity, and the settled
final-resolution spectrum. This test does not establish MHD-PIC AMR.

The method saves burn-in, not the final science interval. A half-resolution 3D
turnover has roughly one eighth the cells and twice the timestep, about one
sixteenth of the mesh update cost before overhead. The final high-wavenumber
cascade is absent immediately after prolongation, so no post-refinement output
is stationary until an empirically measured settling interval has passed.

### 6.3 Current full-Hall multilevel limitation

Full-Hall MHD-PIC deliberately fails closed on multilevel meshes. AthenaK has
useful infrastructure for MHD prolongation/restriction, flux and CT correction,
particle owner remapping, load-cost accounting, topology restart, and bounded
AMR proxy tests. That infrastructure is not a qualified multilevel full-Hall
coupling.

The current receiver-resolution TSC path deposits using each owning block's
resolution and then applies generic exchange/restriction/prolongation. It does
not build one particle shape over the actual composite leaf mesh. Near a
coarse/fine boundary, the cross-interface particle weight sum need not equal
one. That is incompatible with the exact particle/gas momentum and energy
closure required for production.

Part II needs one clear composite-grid deposition/interpolation policy. Do not
retain several experimental policies unless a science case requires them.

### 6.4 Conservative composite-grid PIC operator

The particle stencil must be evaluated over the actual leaf cells on both sides
of a refinement interface. It must provide:

- partition of unity for signed deposited charge;
- correct physical cell-volume normalization;
- one cross-level geometry for \(Q_{\rm cr}\), \(\boldsymbol K_{\rm cr}\),
  DPDT, and DEDT;
- smooth second-order-consistent interpolation/deposition across the planar
  interface;
- equality between deposited DPDT/DEDT and recorded particle momentum and
  kinetic-energy changes; and
- the exact opposite gas momentum and energy changes using leaf-cell volumes.

First-spatial-moment consistency belongs in the derivation and oracle where
practical. Do not repair a bad interface stencil afterward with a global
rescaling.

The full-Hall VL2 chronology remains unchanged across levels: stage 1 deposits
the predictor, stage 2 deposits midpoint moments and realized DPDT/DEDT, gas
receives the opposite realized exchange, and artificial particle light speed
does not enter the grid Hall closure.

### 6.5 Hall CT and total-energy reflux

Hall terms must use AthenaK's multilevel correction paths rather than a second
interface algorithm:

- Hall induction enters face/edge CT correction so coarse magnetic flux agrees
  with the restricted fine update;
- the matched Hall energy flux
  \(\boldsymbol{cE}_H\times\boldsymbol B\) is refluxed through total energy;
- induction and Hall energy are accepted or replaced together under FOFC;
- the stage-2 corrector remains based on deposited DPDT; and
- the face-centered field remains divergence-free through prolongation,
  refluxing, restart, and topology change.

A manufactured planar case must establish flux compatibility, bounded
divergence, and total gas-plus-particle energy closure before a shock is run.

### 6.6 Particle ownership and mesh lifetime

Every physical macro-particle has exactly one leaf-block owner. Refinement
changes owner, not position, momentum, species, weight, identity, birth time,
tracked identity, old position, or exchange channels.

The first implementation does not split or merge particles. Refined cells will
have fewer particles per cell; add a statistically sound splitting policy only
if a production calculation demonstrates unacceptable noise.

Qualify:

- coarse-to-fine and fine-to-coarse crossings;
- ownership changes caused by refinement and derefinement;
- MPI migration when a new owner lies on another rank;
- preservation of all staged particle state;
- restart after a completed topology transition;
- injection weights based on physical leaf volume or surface area; and
- particle-aware load balance based on measured cost.

Topology changes occur only between committed full timesteps after deposition
and migration complete, never inside a VL2 stage.

### 6.7 Static refinement for shocks

Begin coupled refinement with one factor-two planar static layout. Use a
shock-frame geometry where possible; otherwise make the fine region wide enough
for the full shock trajectory.

The fine region covers:

- the shock transition;
- the complete current-bearing Bell precursor;
- the nonlinear cavities and filaments being measured; and
- enough downstream volume for validation diagnostics.

Refining only the discontinuity is physically insufficient because the CR
current drives Bell growth upstream.

The shock tracker and injection ledger must traverse the composite leaf surface
exactly once. Covered coarse faces must not be counted, and physical particle
normalization must not depend on refinement level.

After the planar operator passes, compare one short static-refinement shock with
a uniform calculation at the same finest resolution. Compare shock trajectory,
compression, precursor current, Bell scale and growth, magnetic amplification,
global exchange ledgers, particle retention, and magnetic divergence.

### 6.8 Dynamic AMR

Dynamic AMR begins only if a fixed fine slab wastes enough cells to matter or
cannot follow the required precursor. Reuse the qualified static operator and
add:

- a two-level shock/precursor tag;
- a physical buffer around the tracked shock and current-bearing precursor;
- simple hysteresis to avoid mesh thrashing;
- repeated owner remapping, MPI redistribution, and restart;
- the 3D edge/corner cases required by the production geometry; and
- particle-plus-mesh load balancing.

Start with the shock position and a buffered physical precursor extent. Do not
begin with a collection of noisy current, vorticity, density, and magnetic
thresholds. Add a field-based trigger only when the geometric rule misses
dynamically relevant structure.

### 6.9 Minimal qualification ladder

| Step | Decisive check |
| --- | --- |
| 1. MHD global refinement | One \(64^3\rightarrow128^3\) transition: conservation, forcing continuity, bounded divergence, and post-settling spectra |
| 2. Planar particle interface | A few particle positions on both sides: identity, signed deposited sums, and first moment |
| 3. Rank-split interface | The same crossing with interface and owner transition across MPI ranks |
| 4. Coupled full-Hall interface | Exact gas-particle exchange, Hall CT/energy reflux, and bounded divergence |
| 5. Restart | Restart after a completed crossing/topology transition and recover the next-step result |
| 6. Static shock | One short SMR shock versus a uniform finest-resolution reference |
| 7. Dynamic shock | One buffered moving region versus the qualified static case |
| 8. Frontier performance gate | One production-shaped GPU/MPI profile showing useful end-to-end savings |

Do not create an orientation, precision, particle-count, and refinement-level
matrix. Start with the planar orientation used by the shock. Add a 3D
edge/corner case when the 3D target requires it and a statistical repeat only
when noise limits the comparison.

### 6.10 Refinement performance gate

AthenaK advances multilevel MHD with a global timestep set by the finest cells.
Refinement reduces evolved cells but does not give coarse regions larger
timesteps.

For one factor-two level in 3D, if fraction \(f\) of the volume is refined, the
mesh-cell work relative to a uniform finest mesh is approximately

\[
f+\frac{1-f}{8},
\]

before boundary, communication, reconstruction, and load-balance overhead.
Particle work generally does not decrease merely because the mesh is refined.

If fraction \(p\) of runtime is particle/global work that refinement does not
remove and \(h\) is refinement overhead, a useful first estimate is

\[
\frac{T_{\rm refined}}{T_{\rm uniform}}
\simeq p+(1-p)\left(f+\frac{1-f}{8}\right)+h.
\]

Before promoting refinement, record fine-volume fraction, finest-grid cycle
count, active leaf cells, particle count and particle wall fraction, memory,
output size, and projected end-to-end speedup. Begin Part II for a real capacity
need or a worthwhile total benefit, approximately a factor of two or larger,
not for a mesh-only headline.

### 6.11 Part II completion criterion

Part II is complete only when one immutable candidate:

- retains the Part I full-Hall equation and staging contract;
- conserves signed deposition and exact gas-particle momentum/energy across
  refinement interfaces;
- refluxes Hall induction and energy consistently;
- preserves particle identity and ownership through crossing, topology change,
  MPI redistribution, and restart;
- keeps magnetic divergence bounded;
- reproduces one uniform finest-resolution shock within the intended accuracy;
  and
- demonstrates a measured production benefit.

Until then, the parser continues to reject full-Hall multilevel science runs.

## 7. Decision rules

- If a focused correctness or production benchmark passes, proceed.
- If a result exposes one clear numerical sensitivity, run one targeted repeat.
- If a physical signal is weak but the run is healthy, change the demonstrated
  physical limitation rather than the implementation.
- If signs, conservation, normalization, CT, or the selected equation set are
  wrong, stop and fix them before scaling.
- If a job fails after writing partial output, use it diagnostically but require
  a clean rerun for evidence.
- If optimized uniform performance makes the science affordable, keep Part II
  deferred.
- If refinement becomes necessary, implement whole-domain MHD refinement,
  static coupled refinement, and dynamic refinement in that order.

The immediate next action is to finalize the efficient 3D proof-shock deck and
run its short cost/current-noise pilot from the frozen v2 candidate. If the
projection fits the declared budget, run the proof shock and then move directly
into saturated-MHD-to-CR turbulent-box science. Reopen performance work only in
response to a measured limiter in one of those target runs.
