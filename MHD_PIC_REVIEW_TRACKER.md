# MHD-PIC Uniform Science Tracker

> Temporary execution tracker for `PIC_development`, refreshed 2026-07-13.
> The physical contract, completed evidence, execution details, and deferred
> refinement roadmap live in
> [MHD_PIC_NEXT_STEPS_GUIDE.md](MHD_PIC_NEXT_STEPS_GUIDE.md).
> This file is the concise source of truth for active work and the next decision.

Follow [ethos.md](ethos.md). Physical consistency is mandatory. Reliability and
performance work must protect or accelerate a workflow we intend to use; they
must not grow into speculative hardening or a general benchmark matrix.

## Project decision

Part I uses the supported uniform Cartesian MHD-PIC model and follows this order;
the first item is complete:

1. make the uniform-grid particle coupling fast and production-reliable;
2. run one efficiently designed 3D Bell shock as validation and method proof;
3. drive uniform-grid MHD turbulence to saturation at the science resolution;
4. inject CR particles into the common saturated state and begin the main
   turbulent-box science.

The shock is a validation result for the method section, not the primary science
campaign. Static mesh refinement (SMR) and adaptive mesh refinement (AMR) are
deferred to Part II of the guide and do not gate Part I.

## 1. Active Part I work

### `PERF-U1` — Remove replicated shock-source preparation

- Status: **COMPLETE**
- Production-shaped scaling probe `4980740` completed through `t=47` on 32
  Frontier nodes. Injection increased cost from `0.05766` to `0.14010` seconds
  per cycle. The run remained physically healthy.
- The qualified uniform-grid path now caches planar geometry, constructs only
  locally owned stencil work, and routes one canonical tag-derived particle
  stream instead of repeating it on every rank. The legacy/AMR path is unchanged.
- Focused parity, restart-control, and normal-precision Hall tests passed in job
  `4981305`. The identical 32-node repeat `4981353` reduced post-injection cost
  from `0.14010` to `0.12917` seconds per cycle (7.8 percent) and incremental
  injection-window overhead by 13.0 percent. Particle counts were identical,
  invalid records remained zero, and history differences were at roundoff
  (`5.6e-14` maximum scaled difference).

### `PERF-U2` — Improve particle-migration completion

- Status: **COMPLETE — MODEST RETAINED IMPROVEMENT**
- The representative full-Hall `128^3` box spent `21.4%` of total wall time in
  particle migration and made roughly `8.95e8` migration-wrapper calls. The
  evidence points to repeated MPI completion polling, not particle imbalance;
  the box particle balance was `99.85%`.
- The terminal receive task now waits once for its posted real and integer
  receives instead of returning into a task-list busy poll. Focused migration,
  restart-ledger, and telemetry tests passed in job `4981735`.
- Same-node, same-checkpoint treatment/control jobs `4981821` (`t=5`) and
  `4981838` (`t=10`) completed 300 cycles per case. Histories were bitwise
  identical, all `12,582,912` particles remained present, and invalid records
  remained zero. At `t=10`, wrapper calls fell by `96.6%`, fitted and median
  cycle cost fell by `3.1%` and `3.4%`, and total driver time fell by `5.4%`;
  migration-section time itself fell by only `1.6%`.
- The extreme call count in the original uninterrupted run was not reproduced
  after restart, so it is not used to claim a `21.4%` recovery. The small,
  terminal-task change is retained as a modest end-to-end win; further migration
  work is deferred unless a target profile again exposes it as material.

### `PERF-U3` — Particle-aware uniform-grid redistribution

- Status: **COMPLETE — RETAINED FOR LOCALIZED UNIFORM SHOCKS**
- The 100-cycle, 32-node post-injection profile `4982464` measured only `3.93%`
  particle-load efficiency. Push plus deposition occupied `26.1%` of the
  critical-rank driver time, while array resizing occupied only `0.55%`.
  Capacity-managed particle growth is therefore not warranted.
- The existing nonadaptive same-level redistribution path was extended to the
  uniform mesh without enabling multilevel physics. Full Hall remains on its
  supported uniform-grid model, and the option remains off by default.
- Matched job `4982553` redistributed once at the start of the same 100-cycle
  window. Driver time fell from `16.248` to `11.583` seconds (a `1.40x`
  speedup), maximum particles per rank fell by `37.9%`, and push and deposition
  maxima fell by about `51%`. Particle counts were identical, invalid records
  remained zero, and histories agreed to a maximum scaled difference of
  `5.6e-14`.
- The tracked particle-storage high-water mark rose by `34.3%`, to only `0.575`
  GiB per rank in this case. Use a sparse, measured cadence for localized shocks
  and remeasure memory for the final geometry. Do not enable it for an already
  balanced turbulent box.

### `PERF-U4` — Profile the target shock and box workloads

- Status: **COMPLETE — TARGET PROFILES QUALIFIED; STOP GENERAL OPTIMIZATION**
- The localized shock and volume-filling CR-loaded box were profiled with
  synchronized timers separating source preparation, append/resize,
  gather/push/deposition, migration, imbalance, memory, and driver time.
- Synchronized 100-cycle box profile `4981885` retained all `12,582,912`
  particles with `99.8%` balance and zero invalid records. Deposition consumed
  `38.9%` of driver time, push/gather `16.5%`, and migration `13.7%`; particle
  work therefore accounted for about `69%`. Box load balancing is not indicated.
- A complete-stencil TSC renormalization shortcut passed focused normal- and
  single-precision checks (`4981910`, `4981942`) but was rejected after the
  production A/B `4981964`: histories were bitwise identical, while fitted and
  median cycle costs worsened by `2.3%` and `1.8%`. The source change was
  reverted.
- The flag-only Frontier/gfx90a A/B `4982226` rebuilt the same source snapshot
  with `-munsafe-fp-atomics` and ran the same 300-cycle `t=10` checkpoint on one
  node. Fitted, median, and driver costs fell by `13.3%`, `13.2%`, and `13.3%`;
  stage cost fell by `15.2%`. Both cases reached cycle `17317`, retained all
  `12,582,912` particles with zero invalid records, and produced identical
  written histories. HIP double-precision full-Hall conservation and the
  two-GPU single-precision turbulence smoke then passed in `4982256`.
- Retain `-munsafe-fp-atomics` in the Frontier ROCm 6.2.4/gfx90a production
  build profile. It is a platform-specific build choice, not an MHD-PIC
  algorithm change; requalify it when the GPU architecture or ROCm toolchain
  changes.
- The optimized 3D shock profile `4982464` retained `8,098,700` particles with
  zero invalid records. Source preparation consumed `18.6%` of driver time,
  push plus deposition `26.1%`, migration `8.9%`, and particle-array resizing
  only `0.55%`. Source preparation's inclusive append component was just
  `1.11%`, so neither capacity growth nor a broader source rewrite is justified.
- The one evidence-driven follow-up was the successful same-level
  redistribution in `4982553`, recorded under `PERF-U3`. Both intended workload
  classes have now been profiled. Stop general optimization and reopen sorting,
  layout, kernels, staging, or migration only if the final shock or box exposes
  a material limiter.

### `RELIABILITY-U1` — Qualify the production uniform-grid path

- Status: **COMPLETE FOR THE PERFORMANCE BATCH**
- Exact-final job `4982681` passed the normal-precision full-Hall conservation
  test and the two-GPU MPI single-precision turbulence smoke with the optimized
  Frontier build and `-munsafe-fp-atomics`.
- Preserve the normal-precision focused core regression and MPI
  single-precision GPU turbulence smoke after each coherent implementation
  batch. Performance changes must not alter the supported physical model,
  discrete gas-particle exchange, Hall closure, or CT behavior.
- Exercise checkpoint/restart, bounded output, and the intended decomposition in
  one production-shaped worker workflow before the long shock and boxes. Measure
  memory and I/O at the same time; do not create separate campaigns for them.
- Fix observed crashes, hangs, particle loss, corruption, conservation failures,
  or silent physics changes. Do not build a Cartesian product of platforms,
  precisions, particle counts, and decompositions.

### `SHOCK-U1` — Efficient 3D proof-of-concept shock

- Status: **NEXT — FINALIZE THE 3D PILOT**
- Replace the current narrow `5000 by 128 by 128` geometry-check target with one
  science-per-node-hour design. The provisional design space is:
  - `dx` near `6`, giving about 16 cells per injection gyroradius and about 81
    cells across the measured 2D dominant Bell mode;
  - a transverse domain near `12` to `16` injection gyroradii, large enough for
    several Bell structures and genuinely 3D tubes and cavities;
  - approximately `8` to `12` downstream particles per cell unless a short
    current-noise measurement requires more; and
  - evolution to roughly `t=700` to `1000`, or the first time that the declared
    morphology, acceleration, and escape measurements become decisive.
- The run supports a method-section claim that the qualified 2D amplification,
  polarization, cavities, filaments, and CR transport survive in 3D. It is not
  intended to establish a converged high-energy CR spectrum or a broad shock
  parameter survey.
- The working estimate after focused performance work is `500` to `1500`
  Frontier node-hours for the fiducial. Keep the fiducial plus one
  evidence-driven repeat within a `3000` node-hour planning envelope.
- Repeat resolution or particle statistics only if that uncertainty limits the
  method claim. The measured shock Hall parameter remains too small to justify a
  matched Hall-off 3D shock.

### `BOX-U1` — Saturated-turbulence CR campaign

- Status: **PRIMARY SCIENCE, AFTER `SHOCK-U1`**
- Drive the MHD dynamo without CR particles on the final uniform science mesh.
  Before launching, define saturation/statistical stationarity using a compact
  interval of magnetic and kinetic energies, RMS quantities, and spectra.
- Save a reproducible saturated MHD checkpoint and use it as the common initial
  condition for every matched MHD-PIC box. Introduce CRs only after saturation,
  then allow a declared adjustment interval before collecting stationary CR
  statistics.
- Full Hall is the default physical model. Add a Hall-off control only when the
  measured `Lambda` distribution or a specific scientific claim makes the
  comparison informative.
- Measure CR transport and acceleration, magnetic and velocity statistics,
  gas-CR energy and momentum exchange, Hall-regime diagnostics, and spectra.
  Increase resolution, particle number, duration, or physical parameter coverage
  in response to the science question, not as a general code-validation matrix.

## 2. Decisive completed evidence

- **Frozen uniform candidate:** tag `pic-science-candidate-20260713-v2` contains
  the qualified full-Hall implementation, performance corrections, Frontier
  build template, focused tests, decks, and documentation. The previous
  `pic-science-candidate-20260713` tag remains immutable historical evidence.
- **Full CR-Hall closure:** linear suite `4975948`, exact-source/core job
  `4978017`, and two-rank single-precision GPU smoke `4975966` passed.
- **Nonlinear Bell onset:** job `4975940` qualified onset through the declared
  bounded regime. The deep coarse-pilot filament collapse remains outside that
  claim and does not motivate speculative source repair.
- **Controlled 2D shock:** jobs `4976035` and `4977996`, with analysis `4978079`,
  established a stable shock and coherent Bell amplification, scale selection,
  circular polarization, filaments, and cavities. The long run reached
  `N_Bell=6.64`, helicity `-0.998`, and `pmax/pinj=2.45`.
- **Turbulent-box Hall response:** full-Hall `4976062` and matched Hall-off
  `4976612` completed through `t=10`. The comparison found a measurable magnetic
  response and supplied the representative uniform-box profile.
- **Correctness cleanups:** invalid CR macro weights are rejected, legacy
  nonphysical Bell mechanics require explicit opt-in, and focused test discovery
  is repaired.

These results are not rerun after every deck or documentation edit. Requalify
after a change to the coupling physics, discrete conservation, particle ownership,
or a relevant production hot path.

## 3. Part I decision rules

- Uniform Cartesian, full-f, non-relativistic ideal MHD with VL2/TSC and full
  CR-Hall coupling is the supported science model.
- Optimize a measured production-shaped bottleneck, then repeat that workload.
- Preserve physical results while improving throughput; never trade correctness
  for a favorable timing number.
- Use one fiducial calculation and one evidence-triggered repeat, not a general
  parameter matrix.
- Long runs and substantial analysis execute as self-contained worker jobs.
  Occasional targeted scheduler checks are fine; do not poll long jobs in a loop.
- Physics runs use an immutable tagged commit and retain the exact input,
  executable identity, concise report, and analysis needed for the claim.

## 4. Part II handoff: mesh refinement

### `REFINE-1` — Known limitation, not a Part I gate

Coupled refinement-interface deposition is not yet conservative. The current
receiver-resolution TSC path evaluates levels independently without a cross-level
partition of unity. The exact full-Hall science path therefore remains uniform-grid
only.

Part II of the guide separately specifies:

1. time-triggered whole-domain refinement for a uniform resolution staircase;
2. conservative planar two-level particle crossing and deposition;
3. exact gas-particle momentum and energy closure across levels;
4. Hall CT/EMF and Hall-energy-flux synchronization and refluxing;
5. particle ownership, migration, restart, refinement, and derefinement;
6. fixed nested refinement covering a shock and its Bell precursor; and
7. dynamic criteria and load balancing only after the fixed path is qualified.

A time criterion that refines the entire turbulent box is exactly the “global
resolution staircase” discussed in planning. After each event all leaf blocks are
again at one level, so each epoch is scientifically a uniform-grid calculation.
If used later, the MHD turbulence must settle at the final resolution before CR
injection and before collecting final-resolution statistics.

## Deletion rule

Delete this tracker after the 3D proof shock and first saturated-state MHD-PIC
box are complete and their durable evidence has moved into the guide or permanent
documentation. If Part II becomes active later, create a new refinement tracker
rather than keeping this Part I dashboard indefinitely.
