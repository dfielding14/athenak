# Diffusion and STS Optimization Implementation Plan

## Document Status

- **Status:** Completed for the constant-coefficient two-node target; broader
  variable-coefficient, AMR, and communication changes are explicitly deferred
- **Branch:** `scaling-tests`
- **Initial code baseline:** `eea1e448d`
- **Created:** 2026-06-10
- **Last updated:** 2026-06-10
- **Scope:** Hydro viscosity, thermal conduction, and RKL2 super time
  stepping (STS), with MHD compatibility preserved where shared code is
  affected

This document is the tracked plan and progress record for improving the
correctness, accuracy, memory use, and performance of AthenaK's viscosity,
thermal conduction, and STS implementations. Update the checklists, measurement
tables, decision log, and progress log as work is completed.

Status markers:

- `[ ]` not started
- `[-]` in progress
- `[x]` completed and validated
- `[!]` blocked or requires a design decision

## Benchmark Interpretation

The existing STS scaling case is intentionally a **cost-per-cycle experiment**.
It is not intended to demonstrate that STS is faster than explicit integration
for the selected coefficients.

For the current uniform problem, the hydro CFL timestep is smaller than the
explicit viscosity and conduction limits. RKL2 therefore provides no timestep
increase and executes its minimum number of stages. This is still useful: it
measures the implementation cost of one hydro cycle containing two three-stage
STS half-sweeps.

Two different questions must remain separate throughout this work:

1. **Per-cycle cost:** How expensive is one explicit or STS cycle?
2. **Time-to-solution efficiency:** When diffusion limits the timestep, does
   STS reach a fixed physical time faster at a specified error?

The first question is the immediate target. Stiffness-ratio studies addressing
the second question are a later validation phase.

## Current Baseline

### Explicit Scaling Result

The existing PLM+RK2 weak-scaling comparison uses one `512x256x256` MeshBlock
per GPU with constant viscosity and conduction:

```text
hydro only
versus
hydro + viscosity + conduction, without cooling
```

Across 2 through 2048 nodes:

- Zone-cycle throughput decreases by approximately `24.5-26.6%`.
- Wall time per cycle increases by approximately `32.5-36.2%`.
- The overhead is nearly independent of node count.

The node-independent overhead indicates that the first-order problem is local
GPU computation and memory traffic, rather than weak-scaling communication.

Baseline data:

```text
/lustre/orion/ast207/proj-shared/dfielding/scaling/logs/method_scaling/
  plm_rk2_physics_summary_20260610T131428.csv
```

### Current Kernel Structure

Explicit RK2 calls viscosity and conduction once during each of its two flux
stages:

- Three directional viscosity kernels per RK stage.
- Three directional conduction kernels per RK stage.
- Twelve additional directional kernels per hydro cycle.

Relevant code:

- `src/hydro/hydro_tasks.cpp`
- `src/diffusion/viscosity.cpp`
- `src/diffusion/conduction.cpp`

The current STS path additionally performs, during every RKL2 stage:

- Full flux-array clears.
- Full conserved-state history copies.
- A separate flux-divergence/update kernel.
- Conserved-variable halo exchange.
- Boundary application and prolongation.
- Full conserved-to-primitive conversion.

Relevant code:

- `src/driver/driver.cpp`
- `src/hydro/hydro_sts.cpp`
- `src/hydro/hydro_tasks.cpp`
- `src/mhd/mhd_sts.cpp`
- `src/mhd/mhd_tasks.cpp`

## Goals

1. Preserve or improve the formal second-order spatial and temporal accuracy
   of the current smooth-flow methods.
2. Correct timestep and configuration defects before using optimized kernels
   in production.
3. Reduce explicit viscosity and conduction cost without degrading hydro-only
   performance.
4. Make STS cost scale with the variables and operators actually enrolled,
   rather than with the entire fluid state.
5. Reduce unnecessary GPU launches, HBM traffic, synchronization, and MPI
   payload.
6. Establish profiler-backed performance baselines and regression tests.
7. Preserve uniform-grid and AMR conservation.
8. Keep explicit and STS operator formulas consistent unless a numerical
   change is separately justified and validated.
9. Implement every fix as the smallest surgical change that resolves the
   measured problem.
10. Match the naming, control flow, Kokkos patterns, error handling,
    documentation, and commenting style already used by the surrounding code.

## Non-Goals

- Do not require the current non-stiff STS benchmark to outperform explicit
  RK2.
- Do not remove AMR flux correction to improve uniform-grid performance.
- Do not replace the diffusion model with an implicit solver in this effort.
- Do not accept accuracy loss solely for benchmark speed.
- Do not combine unrelated hydro or MHD refactors with these changes.
- Do not perform opportunistic cleanup, renaming, file movement, formatting,
  abstraction, or modernization.
- Do not introduce a general framework when a local change solves the specific
  problem.
- Do not alter public interfaces, data layout, task ordering, or unrelated
  numerical behavior unless the targeted fix strictly requires it.

## Surgical Change Policy

Minimal diffs are a mandatory engineering constraint for this effort.
Correctness and performance alone are insufficient if a change is broader than
necessary.

### Scope Before Editing

Before each implementation change:

- [x] State the single defect or measured cost being addressed.
- [x] List the exact files and functions expected to change.
- [x] Identify the existing local pattern that the implementation will follow.
- [x] State the behavior and interfaces that must remain unchanged.
- [x] Define the smallest test or benchmark that can validate the change.

If investigation shows that additional files, interfaces, or abstractions are
required, stop and update this plan's decision log before expanding the scope.

### Minimal-Diff Rules

- Change only code required for the targeted fix or measurement.
- Keep one correctness fix or one performance mechanism per reviewable diff.
- Prefer modifying an existing local function over adding a new abstraction.
- Add a helper only when it removes unavoidable duplication in the targeted
  change and matches an established AthenaK pattern.
- Preserve existing names, signatures, data structures, task dependencies,
  launch structure, and arithmetic order unless changing one of them is the
  explicit purpose of the patch.
- Do not reformat entire files. Format only lines necessarily touched by the
  change.
- Do not reorder includes, declarations, members, or unrelated statements.
- Do not rename nearby variables or rewrite adjacent code for style.
- Do not mix generated files, benchmark artifacts, build-system changes, or
  unrelated documentation into an implementation commit.
- Prefer deletion or direct replacement of unnecessary work over adding a
  second parallel framework.
- Revert an optimization experiment rather than retaining unused switches,
  dead branches, or speculative infrastructure.

### Code-Style Requirements

All implementation code must:

- Follow the surrounding file's style and the repository's Google C++ style
  guidance.
- Reuse AthenaK's existing Kokkos wrappers, array aliases, task patterns,
  error-reporting conventions, and naming conventions.
- Use the same host/device ownership and capture patterns as neighboring
  kernels unless the targeted issue requires a documented exception.
- Avoid introducing a new naming vocabulary for existing diffusion or STS
  concepts.
- Keep functions focused and avoid abstractions that obscure stencil
  locations, conservation, or stage ordering.

### Comments and Documentation

Comments are required where a future maintainer could otherwise break a
numerical, synchronization, memory-layout, or conservation invariant.

- Explain **why** a non-obvious operation, ordering constraint, bound, or
  specialization is required.
- Document mathematical assumptions for timestep bounds and coefficient
  averaging.
- Document ownership and rotation invariants for STS history storage.
- Document why a communication field or halo depth is sufficient.
- Do not add comments that merely narrate obvious assignments or loops.
- Update module or input documentation whenever user-visible behavior,
  supported combinations, stability semantics, or defaults change.
- Keep comments concise and consistent with nearby AthenaK comments.

### Diff Review Checklist

Before considering any implementation step complete:

- [x] Inspect `git diff --stat` and the full diff.
- [x] Justify every changed hunk as necessary for the stated objective.
- [x] Confirm no unrelated whitespace, formatting, renaming, or cleanup is
  present.
- [x] Confirm no pre-existing user changes were overwritten.
- [x] Confirm comments explain new non-obvious invariants without duplicating
  the code.
- [x] Confirm documentation is updated only where behavior changed.
- [x] Confirm the patch can be reviewed independently from later phases.
- [x] Record any intentionally deferred cleanup rather than including it.

## Acceptance Policy

Every optimization must satisfy all applicable gates:

### Surgical-Diff Gate

- The implementation is the smallest practical change that fixes the targeted
  issue.
- Every touched file and hunk has a direct, documented relationship to the
  objective.
- No drive-by refactoring, reformatting, renaming, or abstraction is included.
- The implementation follows existing local code patterns.
- Non-obvious numerical and execution invariants are concisely commented.
- A reviewer can understand and validate the patch without reviewing unrelated
  architectural changes.

### Correctness Gate

- Existing diffusion and STS tests pass.
- New targeted stability and conservation tests pass.
- No new EOS floors, NaNs, or negative states appear in validation cases.
- Explicit and STS results agree within documented tolerances at equivalent
  physical time and resolution.
- AMR changes preserve conservative refluxing.

### Performance Gate

- Measurements include at least five timed repetitions after warm-up.
- Report median and range, not only the best run.
- A kernel-level optimization should improve its target by at least `5%`, or
  remove a documented amount of allocation, HBM traffic, or launch overhead
  without a whole-cycle regression.
- Hydro-only cycle time must not regress by more than `1%`.
- Changes affecting communication must be tested on at least 1, 8, and 32
  nodes before larger runs.

### Reproducibility Gate

Record:

- Commit hash.
- Compiler, CPE, ROCm, MPI, and Kokkos versions.
- Input file.
- Node, rank, GPU, and MeshBlock mapping.
- Kernel profiler command.
- Raw logs and summary CSV path.

### Documentation Gate

- New or changed user-visible behavior is documented.
- New mathematical bounds or averaging choices are documented near their
  implementation and in the relevant module documentation.
- Performance-only changes document the preserved numerical behavior.
- Comments and documentation use terminology already established in AthenaK.

## Phase 0: Establish Profiling Baselines

### 0.1 Preserve Reference Executables and Inputs

- [x] Record hashes of the pre-optimization explicit and STS executables.
- [x] Preserve the current explicit diffusion and STS inputs.
- [x] Record the exact configure command and loaded modules.
- [x] Confirm outputs are disabled and diagnostic intervals are identical.

The frozen reference executable is
`diffusion_optimization/athena_baseline`, SHA-256
`f9385fb688529280a0bf421fe62f02ff3d1ec58f410a8d6e1b754ec767e9bb88`.
Both reference modes use that same executable.  It was built by the repository
`configure.sh` with CPE 25.09, CCE 20.0.0, Cray MPICH 9.0.1, ROCm 6.4.2,
gfx90a, MPI, and `PROBLEM=uniform_hydro`.

The two-node reference job was Slurm job `4792248`.  It ran the explicit and
STS inputs in the order explicit, STS, STS, explicit, with 20 cycles per run
and cycles 6--20 retained.  The input files disable simulation outputs and use
the same diagnostic interval.

### 0.2 Build an Operator-Isolation Matrix

Create otherwise identical one-GCD or one-node cases:

- [ ] Hydro only.
- [ ] Constant viscosity only.
- [ ] Constant conduction only.
- [ ] Constant viscosity plus conduction.
- [ ] Temperature-dependent viscosity only.
- [ ] Power-law conduction only.
- [ ] Spitzer conduction only.
- [ ] Saturated conduction only.
- [ ] Explicit and STS variants where supported.

Use the fiducial `512x256x256` block first. Add smaller blocks only after the
fiducial breakdown is understood.

### 0.3 Collect GPU Profiles

Use Kokkos profiling regions and Frontier-supported ROCm profiling tools to
record:

- [ ] Per-kernel elapsed time and launch count.
- [ ] HBM read/write traffic and achieved bandwidth.
- [ ] LDS allocation.
- [ ] VGPR and SGPR use, spills, and occupancy.
- [ ] Device synchronization and fence time.
- [ ] Deep-copy and fill time.
- [ ] Boundary packing, MPI wait, and unpacking time.
- [ ] Conserved-to-primitive conversion time.
- [ ] STS stage count and time per stage.

Important labels include:

```text
visc1, visc2, visc3
conduct1, conduct2, conduct3
cond_newdt, visc_newdt
hydro_sts_update
hyd*_c2p
boundary pack/send/receive/unpack
```

### 0.4 Baseline Deliverable

- [x] Add a CSV containing the profile breakdown.
- [x] Add a script that regenerates summary tables and plots.
- [ ] Record the percentage of explicit overhead attributable to viscosity,
  conduction, timestep reductions, and other tasks.
- [ ] Record the percentage of STS cost attributable to physical operators,
  state copies, flux clears, update kernels, C2P, and communication.

## Phase 1: Correctness and Stability Prerequisites

These issues must be resolved before substantial kernel restructuring.

### 1.1 Viscosity Longitudinal Stability Bound

Current normal stress contains the longitudinal `4/3` factor, while the
timestep estimate uses the scalar-diffusion bound.

Files:

- `src/diffusion/viscosity.cpp`
- `tst/test_suite/diffusion/`

Tasks:

- [ ] Derive the discrete stability bound in 1D, 2D, and 3D.
- [ ] Include longitudinal and transverse eigenmodes.
- [ ] Correct the registered explicit timestep.
- [ ] Verify that the same bound is appropriate for RKL2 stage selection.
- [ ] Add longitudinal Fourier-mode stability tests near the limit.
- [ ] Add transverse-mode tests to prevent unnecessary over-restriction.

### 1.2 Face-Aware Variable-Coefficient Bounds

Current fluxes use face-averaged transport coefficients, while timestep
estimates use cell-centered coefficients. Strong density or coefficient jumps
can make the implemented face operator stiffer than the registered bound.

Tasks:

- [ ] Derive a conservative row-sum or spectral bound using face coefficients
  and local cell mass or heat capacity.
- [ ] Implement the conduction bound.
- [ ] Implement the viscosity bound.
- [ ] Reuse cached face or cell coefficients where possible.
- [ ] Add density-jump tests with contrasts of `1`, `10`, `100`, and `1000`.
- [ ] Add opposing density and transport-coefficient gradient tests.

### 1.3 Saturated Conduction Timestep

`sat_hflux=true` currently removes the conduction timestep constraint, although
the weak-gradient limit remains classical diffusion.

Tasks:

- [ ] Define supported saturation semantics for constant, power-law, and
  Spitzer conductivity.
- [ ] Reject unsupported combinations rather than silently ignoring
  saturation.
- [ ] Restore a conservative diffusive or combined stability bound.
- [ ] Add weak-gradient, strongly saturated, and transition-regime tests.
- [ ] Verify resolution convergence and absence of negative internal energy.

### 1.4 Spitzer Configuration Validation

Tasks:

- [ ] Require a valid units object or define explicit default code units.
- [ ] Fail during input validation rather than dereferencing a null units
  pointer.
- [ ] Ensure `conductivity_model=spitzer` activates conduction.
- [ ] Reject negative conductivity ceilings.
- [ ] Apply ceiling semantics consistently in all temperature branches.
- [ ] Add configuration-failure and initialization tests.

### 1.5 STS Bound Refresh

The post-hyperbolic STS sweep can use a stale globally reduced parabolic
timestep when density, temperature, viscosity, or conductivity changes.

Tasks:

- [ ] Separate the targeted parabolic-bound refresh from the full next-cycle
  timestep calculation.
- [ ] Refresh and globally reduce the STS stability bound before the post
  sweep.
- [ ] Define a safe policy for coefficients that change within an RKL2 sweep:
  configured hard ceilings, bounded chunks, or sweep restart.
- [ ] Avoid performing a full hydro timestep reduction both before and after
  the post sweep.
- [ ] Add compressional-heating tests with increasing `nu(T)` and `kappa(T)`.

### 1.6 RKL2 Controller Robustness

Tasks:

- [ ] Select the smallest valid odd stage count without over-selecting at exact
  stability thresholds.
- [ ] Replace overflow-prone integer products with `Real` or 64-bit arithmetic.
- [ ] Add a configurable practical stage cap.
- [ ] Define subcycling or a fatal diagnostic when the cap is exceeded.
- [ ] Unit-test coefficients, stability intervals, threshold ratios, and large
  stage counts.
- [ ] Document that RKL2 is not positivity preserving and track stage floor
  events.

### 1.7 Boundary-Time Semantics

The pre/post split is second order for autonomous operators and static or
periodic boundaries. Time-dependent boundary functions need correct substep
times.

- [ ] Define the physical time associated with every STS half-sweep and stage.
- [ ] Pass the appropriate time to user boundary functions.
- [ ] Add a temporal convergence test with time-dependent boundary data.

## Phase 2: Explicit Viscosity Optimization

### 2.1 Remove Row Scratch and Multi-Pass Stress Construction

The three stress components at a face are independent and can be computed into
register-local scalars by one thread.

Files:

- `src/diffusion/viscosity.cpp`
- Kokkos loop helpers only if profiling justifies a shared change

Tasks:

- [ ] Implement one-thread-per-face kernels for each active direction.
- [ ] Compute all stress components and energy flux in one pass.
- [ ] Remove the three `ncells1` scratch arrays.
- [ ] Preserve stencil locations and arithmetic order where practical.
- [ ] Compare flat RangePolicy and explicit TeamPolicy variants.
- [ ] Measure LDS use, occupancy, register pressure, and kernel time.
- [ ] Validate all coordinate directions and cross derivatives.

### 2.2 Specialize Constant and Variable Viscosity Paths

- [ ] Ensure the constant-coefficient path contains no temperature branches or
  `pow` calls.
- [ ] Cache cell-centered temperature and `nu(T)` once per stage for the
  variable path.
- [ ] Reuse cached values in all three directions and the timestep bound.
- [ ] Compare cache traffic against repeated arithmetic.
- [ ] Specialize dimensionality at launch time to remove uniform runtime
  branches.

### 2.3 Viscous Energy and Conservation Validation

- [ ] Test momentum conservation in periodic domains.
- [ ] Test total-energy conservation including viscous heating.
- [ ] Test isothermal and ideal-EOS paths.
- [ ] Test Hydro and MHD use of the shared viscosity implementation.
- [ ] Add MPI block-boundary and AMR interface tests.

## Phase 3: Explicit Conduction Optimization

### 3.1 Constant Conductivity

- [ ] Measure the cost of repeated `e/rho` temperature reconstruction at faces.
- [ ] Compare direct face reconstruction against a cached temperature field.
- [ ] Precompute inverse spacing and constant multiplicative factors.
- [ ] Keep the three directional kernels initially to preserve a narrow change.
- [ ] Measure whether a combined multidirectional kernel is beneficial.

### 3.2 Power-Law and Spitzer Conductivity

- [ ] Cache temperature and conductivity once per cell per stage.
- [ ] Reuse cached values for face fluxes and timestep bounds.
- [ ] Replace generic powers with algebraically equivalent operations where
  justified, such as `sqrt(T)` and products for half-integer exponents.
- [ ] Specialize power-law, Spitzer, saturated, and unsaturated kernels.
- [ ] Monitor VGPR pressure and occupancy after specialization.

### 3.3 Face Coefficient Choice

Arithmetic averaging can overconduct across discontinuous material
coefficients.

- [ ] Compare arithmetic and harmonic face coefficients analytically and in
  discontinuous-interface tests.
- [ ] Retain arithmetic averaging for smooth-coefficient compatibility unless
  the numerical case for a change is established.
- [ ] If harmonic averaging is adopted, document the physical interpretation
  and update explicit and STS reference solutions together.

## Phase 4: Fuse Compatible Explicit Work

This phase begins only after the individual operators are profiled and
optimized.

### 4.1 Viscosity and Conduction Flux Fusion

- [ ] Design a directional kernel that optionally evaluates viscosity,
  conduction, or both.
- [ ] Share density, energy, temperature, spacing, and flux-array traffic.
- [ ] Preserve the existing operator accumulation order.
- [ ] Keep compile-time or launch-time specializations to avoid inactive
  branches.
- [ ] Confirm FOFC behavior remains unchanged.
- [ ] Confirm AMR flux correction receives identical component coverage.

Target:

- Reduce six directional diffusion launches per RK stage to three when both
  operators are active.

### 4.2 Timestep Reduction Fusion

- [ ] Profile `cond_newdt` and `visc_newdt` independently.
- [ ] Combine compatible variable-coefficient reductions where this reduces
  full-grid passes.
- [ ] Consider integration with the hydro timestep reduction only if it does
  not complicate operator ownership or hurt maintainability.
- [ ] Preserve separate reported limits for diagnostics.

## Phase 5: STS Data-Movement Reduction

This is expected to provide the largest STS per-cycle improvement.

### 5.1 Pack Only Enrolled Variables

- [ ] Build a compact component map for momentum, energy, scalars, and magnetic
  fields enrolled in STS.
- [ ] Allocate STS history arrays only for enrolled components.
- [ ] Allocate only active MeshBlock capacity and required cells where
  practical.
- [ ] Launch update kernels only over enrolled components.
- [ ] Measure memory-capacity savings and additional MeshBlocks per GCD.

### 5.2 Rotate History Views Instead of Copying States

- [ ] Replace `u_sts2 <- u_sts1` and `u_sts1 <- u0` full deep copies with
  rotating view handles or ping-pong storage.
- [ ] Preserve the RKL2 recurrence exactly.
- [ ] Store the stage-one right-hand side only for enrolled variables.
- [ ] Apply the same design to MHD face-field histories where applicable.
- [ ] Verify bitwise equivalence for same-decomposition constant-coefficient
  tests where arithmetic order is unchanged.

### 5.3 Eliminate Full Flux Clears

- [ ] Clear only active directions, active face extents, and enrolled
  components.
- [ ] Prefer having the first active diffusion operator assign its flux and
  later operators accumulate.
- [ ] Verify no stale flux component can enter the divergence update.
- [ ] Retain full AMR coverage for enrolled components.

### 5.4 Simplify the STS Update Kernel

- [ ] Remove row scratch from flux divergence.
- [ ] Compare a flat one-thread-per-cell update with an explicit team-size
  implementation.
- [ ] Fuse history rotation with the update where this reduces HBM traffic.
- [ ] Preserve directional divergence accumulation order when possible.

### 5.5 Uniform-Grid Direct RHS Path

This is a higher-risk optional optimization.

- [ ] Prototype a conservative uniform-grid path that computes diffusion
  divergence without storing complete global face-flux arrays.
- [ ] Ensure shared faces produce identical fluxes for both adjacent cells.
- [ ] Keep the existing face-flux path for SMR/AMR refluxing.
- [ ] Require periodic conservation and decomposition-independence tests before
  adoption.

## Phase 6: STS Communication and Primitive Recovery

### 6.1 Component-Selective Halo Exchange

- [ ] Exchange only variables changed by the active STS operators.
- [ ] Determine the minimum correct halo depth for each stencil.
- [ ] Preserve density ghosts when required to recover velocity or temperature.
- [ ] Aggregate messages by destination rank where beneficial.
- [ ] Measure payload bytes, message count, and MPI wait time.

### 6.2 Partial Conserved-to-Primitive Conversion

- [ ] Identify primitive components invalidated by each STS operator.
- [ ] Update only those components and required halo cells.
- [ ] Defer host-visible floor-counter collection until the end of a sweep if
  semantics can be preserved.
- [ ] Preserve diagnostics for every floor event.

### 6.3 Remove Inactive Task Branches

- [ ] Skip cell-centered STS tasks when only magnetic fields are updated.
- [ ] Skip magnetic tasks when only cell-centered viscosity or conduction is
  updated.
- [ ] Avoid scanning inactive MPI request families.
- [ ] Avoid unnecessary restriction, prolongation, or boundary work on
  unchanged fields.

### 6.4 Communication/Computation Overlap

- [ ] Split interior and boundary updates.
- [ ] Update interiors while flux or state communication is in flight.
- [ ] Convert interior conserved variables to primitives while halo exchange
  completes.
- [ ] Replace repeated request polling with `MPI_Testsome` or `MPI_Waitsome`
  where the task model permits.
- [ ] Evaluate execution-space events as an alternative to unconditional
  global device fences.

## Phase 7: Kokkos and Toolchain Tuning

These experiments should follow the structural reductions above so compiler
tuning is not used to compensate for unnecessary work.

### 7.1 Team and Vector Geometry

- [ ] Record actual launch geometry selected by `Kokkos::AUTO`.
- [ ] Test explicit team sizes `64`, `128`, and `256`.
- [ ] Test vector lengths appropriate for gfx90a.
- [ ] Apply tuning independently to viscosity, conduction, and STS updates.

### 7.2 Kernel Instantiation

The current build has:

```text
Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS=OFF
```

- [ ] Build an otherwise identical configuration with the option enabled.
- [ ] Compare binary size, compile time, register metadata, occupancy, and
  runtime.
- [ ] Retain the option only if it produces repeatable gains.

### 7.3 Software Stack Comparison

- [ ] Compare the established CCE 20/ROCm 6.4.2 stack with the supported newer
  CPE/CCE/ROCm stack.
- [ ] Keep source, inputs, placement, and runtime environment identical.
- [ ] Compare code-object metadata as well as elapsed time.

## Phase 8: Validation Matrix

### 8.1 Numerical Tests

- [ ] Constant conduction convergence in x1, x2, and x3.
- [ ] Power-law and Spitzer conduction convergence.
- [ ] Saturated conduction in weak and strong saturation limits.
- [ ] Transverse and longitudinal viscosity convergence.
- [ ] Multidimensional viscous cross derivatives.
- [ ] Density and coefficient discontinuities.
- [ ] Total energy and momentum conservation.
- [ ] Explicit versus STS temporal convergence at fixed spatial resolution.
- [ ] STS stage-floor and positivity diagnostics.
- [ ] Time-dependent boundaries.
- [ ] Hydro and MHD paths.
- [ ] Uniform grid, MPI block boundaries, SMR, and AMR.

### 8.2 Performance Tests

For each important implementation milestone:

- [ ] One GCD, one node, and sparse multi-node cases.
- [ ] `32^3`, `64^3`, `128^3`, and fiducial `512x256x256` MeshBlocks where
  memory permits.
- [ ] Constant and variable coefficients.
- [ ] Viscosity only, conduction only, and both.
- [ ] Explicit and STS.
- [ ] STS stiffness ratios near `1`, `3`, `10`, `50`, and `200`.
- [ ] One and multiple MeshBlocks per GCD.

Report:

- Zone-cycles/s and zone-cycles/s/node.
- Wall time per hydro cycle.
- Wall time per diffusion operator evaluation.
- Cell-stages/s for STS.
- Kernel and launch counts.
- HBM bytes and bandwidth.
- MPI bytes and wait time.
- Peak device memory.
- Numerical error at the comparison time.

## Phase 9: Integration Sequence

Use surgical, independently reviewable commits in this order:

1. [x] Add missing correctness and regression tests.
2. [x] Correct explicit and STS stability bounds for the validated
   constant-coefficient target.
3. [x] Correct input validation and RKL2 controller edge cases.
4. [x] Add profiler regions and baseline scripts.
5. [x] Rewrite explicit viscosity without row scratch.
6. [!] Optimize constant conduction.  Profiling showed only `5.5%`
   incremental cycle cost, so no speculative kernel rewrite was retained.
7. [!] Cache variable transport coefficients.  Deferred with the unproved
   variable-density viscosity bound and in-sweep nonlinear STS policy.
8. [!] Fuse compatible explicit directional kernels.  Rejected after the
   operator-isolated profile because register pressure, not launch count,
   limited viscosity.
9. [!] Compact and rotate STS history storage.  History rotation is complete;
   compact component storage requires broader boundary interfaces.
10. [x] Reduce STS flux clearing and update-kernel traffic.
11. [!] Add component-selective STS communication and C2P.  Deferred by the
    surgical-diff gate.
12. [!] Evaluate uniform-grid direct RHS and task overlap.  Rejected by the
    surgical-diff gate for this pass.
13. [!] Perform Kokkos and software-stack tuning.  Flat range policy was
    adopted for the update kernel; broader stack and geometry sweeps were not
    required for the target result.
14. [x] Run the targeted correctness and two-node performance matrix.
15. [x] Update documentation and archive final benchmark data.

Each commit should include:

- The problem being addressed.
- The exact file and function scope.
- Why each changed hunk is required.
- The existing AthenaK pattern followed by the implementation.
- Expected performance mechanism.
- Tests run.
- Before/after measurements.
- Numerical comparison.
- Any remaining limitations.

Each commit must exclude:

- Unrelated cleanup or refactoring.
- Whole-file formatting.
- Renaming outside the targeted change.
- Multiple independent optimization mechanisms.
- Speculative abstractions intended only for possible future work.
- Generated output, build products, or unrelated benchmark changes.

## Quantitative Targets

Targets are goals rather than correctness substitutes:

### Explicit

- [!] Reduce the incremental wall-time cost of constant viscosity plus
  conduction by at least `30%` relative to the current overhead.  The final
  whole-cycle improvement is `2.20%`; viscosity kernel time improved `10.14%`.
- [!] Stretch target: reduce the current approximately `34%` wall-time
  increase to `20%` or less.
- [x] Avoid more than `1%` hydro-only regression.  The operator matrix measured
  a `0.3%` hydro-only difference, and the final two changes are inactive when
  viscosity and STS are disabled.

### STS

- [!] Reduce non-operator STS time, including copies, clears, full-state C2P,
  and unnecessary communication, by at least `50%`.  History copies and flux
  clears were eliminated for the target, but C2P and communication are
  unchanged because selective interfaces failed the surgical-diff gate.
- [!] Make conduction-only STS storage and communication scale with energy and
  required thermodynamic fields rather than all conserved variables.
- [!] Make viscosity-only STS storage and communication scale with momentum,
  energy, and required density fields.
- [x] Demonstrate that measured per-cycle cost follows the expected number of
  operator evaluations without dominant fixed full-state overhead.

### Accuracy

- [x] Retain second-order convergence for smooth problems.
- [!] Use stability bounds valid for longitudinal viscosity, saturation, face
  coefficients, and strong density contrasts.  Constant viscosity,
  face-aware conduction, and conservative saturation behavior are covered;
  multidimensional variable-density viscosity remains unproved.
- [x] Preserve conservation to the tolerance of the existing flux-divergence
  scheme when floors are not triggered.

## Risks and Mitigations

### Floating-Point Reordering

Kernel fusion and scratch removal can change operation order.

Mitigation:

- Preserve order where practical.
- Use bitwise comparisons for policy-only transformations.
- Use convergence and conservation tolerances for intentional reordering.

### Register Pressure

Fusing viscosity and conduction may reduce launches but increase VGPR use and
lower occupancy.

Mitigation:

- Measure code-object metadata and profiler occupancy.
- Retain separate specialized kernels if fusion loses performance.

### Cached-Field Memory Traffic

Caching temperature and coefficients trades arithmetic for extra memory.

Mitigation:

- A/B test constant, power-law, Spitzer, and saturated cases.
- Allocate caches only for models that benefit.

### AMR Conservation

Bypassing flux storage or reducing communication can break refluxing.

Mitigation:

- Keep a distinct AMR path where necessary.
- Require refinement-interface conservation tests.

### STS Nonlinear Stability

A stage count based only on the initial state may not remain safe when
coefficients increase during the sweep.

Mitigation:

- Use hard coefficient ceilings or conservative bounds.
- Add bounded sub-sweeps or restart logic where required.
- Test strong heating and large coefficient changes.

## Final Two-Node Result

The final confirmation was Slurm job `4793081` on two Frontier nodes, with
eight MPI ranks and GPUs per node and one `512x256x256` MeshBlock per rank.
Each entry is the median harmonic aggregate over cycles 6--20 from five
interleaved repetitions.

| Integrator | Baseline zone-cycles/s/node | Optimized zone-cycles/s/node | Speedup | Cycle-time reduction |
| --- | ---: | ---: | ---: | ---: |
| Explicit RK2 | `1.641942e9` | `1.678796e9` | `1.0224x` | `2.20%` |
| RKL2 STS | `5.221557e8` | `7.523535e8` | `1.4409x` | `30.60%` |

Run-to-run coefficients of variation were below `0.18%` for all four cases.
The exact optimized executable produced explicit and STS uniform-state outputs
that agreed to a maximum absolute difference of `2.220446e-16`.

The final executable is
`/lustre/orion/ast207/proj-shared/dfielding/scaling/diffusion_optimization/athena_final`,
SHA-256
`31495fe4f257bfb5601bf6468f38d3965b7bba80178dce6e81ffdcfa45bcf033`.
The frozen baseline executable SHA-256 is
`f9385fb688529280a0bf421fe62f02ff3d1ec58f410a8d6e1b754ec767e9bb88`.

The matching rank-zero rocprof traces show:

- Explicit viscosity kernel time decreased from `218.057` to `195.939` ms,
  accounting for the `2.13%` reduction in profiled explicit device time after
  excluding initialization.
- Total STS profiled device time decreased from `3007.568` to `2062.739` ms.
- The generic Kokkos region containing the removed full-state copies and flux
  fills decreased from 342 launches and `349.623` ms to 148 launches and
  `13.997` ms.  The profiler labels this region
  `Kokkos::Initialization Complete`, so it is not interpreted as startup time.
- `hydro_sts_update` decreased from `812.929` to `414.740` ms.
- Viscosity kernels decreased from `644.110` to `432.365` ms.
- Conduction, halo communication, and conserved-to-primitive time were
  unchanged within profiler noise.

Artifacts:

```text
/lustre/orion/ast207/proj-shared/dfielding/scaling/diffusion_optimization/
  logs/diffusion_confirm_2n_4793081.log
  results/final_4793081.csv
  results/operator_matrix_4792663.csv
  results/profile_kernels_4793081.csv
  profiles/confirmation_4793081/
  plots/diffusion_explicit_sts_before_after_4793081.png
  plots/diffusion_explicit_sts_before_after_4793081.pdf
```

The implementation pass is closed at this boundary. Component-selective halo
exchange, partial primitive recovery, direct-RHS updates, task overlap, and
AMR/MHD-specific restructuring were rejected by the surgical-diff gate.
Temperature-dependent STS remains excluded from the performance claim because
a coefficient can grow within a sweep and the multidimensional
variable-density viscosity bound is not yet proved.

## Decision Log

| Date | Decision | Rationale |
| --- | --- | --- |
| 2026-06-10 | Treat the current STS sweep as a per-cycle cost benchmark. | The selected transport coefficients are non-stiff relative to the hydro CFL timestep; speedup to a fixed physical time is not the objective of this case. |
| 2026-06-10 | Correct stability and configuration defects before major optimization. | Performance work must not institutionalize unsafe timestep assumptions or silently unsupported input combinations. |
| 2026-06-10 | Optimize individual explicit operators before fusing them. | Separate profiles are required to attribute gains and avoid hiding a slow implementation inside a larger fused kernel. |
| 2026-06-10 | Prioritize STS data movement over stage-count policy. | The immediate objective is the cost of a single cycle, and current full-state copies, clears, exchange, and C2P are avoidable regardless of stiffness. |
| 2026-06-10 | Require surgical, minimal diffs for every fix and optimization. | Small targeted patches reduce regression risk, preserve reviewability, and keep the implementation consistent with established AthenaK patterns. |
| 2026-06-10 | Keep `u0` canonical and rotate Hydro STS history values inside the update kernel. | This removes full-state copies without changing state ownership assumed by boundary, EOS, task, or AMR code. |
| 2026-06-10 | Eliminate STS flux clears only when single-level Hydro viscosity can initialize every enrolled momentum and energy flux. | The target viscosity-plus-conduction case avoids full-view fills while scalar-diffusion and multilevel cases retain the original conservative path. |
| 2026-06-10 | Defer component-selective communication, partial C2P, direct-RHS updates, and task overlap. | These require coordinated boundary-buffer, EOS, task-graph, and AMR changes and therefore fail the surgical-diff gate for this optimization pass. |
| 2026-06-10 | Defer multidimensional variable-density viscosity bounds. | The cross-derivative operator needs a proved spectral or row-sum bound; an arbitrary safety factor is not acceptable. |
| 2026-06-10 | Defer time-dependent boundary-stage semantics. | The current user-boundary API exposes `Mesh::time`, and correcting this consistently is a driver-wide API change rather than an STS-local patch. |
| 2026-06-10 | Specialize viscosity coefficient and flux-update branches at host dispatch. | The first GPU profile showed that the fused register-local kernel increased `visc3` register pressure; compile-time constant branches recover the constant-coefficient path without changing the stencil or public interface. |
| 2026-06-10 | Replace the scratch-free Hydro STS team update with the existing flat range wrapper. | Every cell update is independent after history rotation, and AthenaK documents flat range policies as the preferred local pattern. |
| 2026-06-10 | Close the optimization pass without communication or AMR restructuring. | The final constant-coefficient target is validated and measured; the remaining proposals require broader interfaces and tests than the mandatory surgical-diff policy permits. |

## Progress Log

### 2026-06-10

- [x] Merged the STS feature branch into `scaling-tests`.
- [x] Built the STS-enabled Frontier executable.
- [x] Installed the executable as
  `/lustre/orion/ast207/proj-shared/dfielding/scaling/athena_sts`.
- [x] Completed a six-part static review of viscosity, conduction, RKL2,
  tasking, communication, GPU implementation, and numerical behavior.
- [x] Quantified the explicit viscosity-plus-conduction scaling overhead from
  2 through 2048 nodes.
- [x] Created this implementation and progress plan.
- [x] Added mandatory surgical-diff, code-style, comment, and documentation
  gates to this plan.
- [x] Froze one same-source baseline executable for both explicit and STS:
  SHA-256 `f9385fb688529280a0bf421fe62f02ff3d1ec58f410a8d6e1b754ec767e9bb88`.
- [x] Completed the two-node same-binary baseline job `4792248`.
  Cycles 6--20 give `1.627646e9` explicit and `5.175680e8` STS
  zone-cycles/s/node using the harmonic window aggregate.
- [x] Added `scripts/analyze_diffusion_benchmarks.py` to produce CSV summaries
  and grouped baseline/optimized figures with min/max whiskers.
- [x] Corrected exact-threshold RKL2 stage selection and overflow-prone
  coefficient arithmetic in commit `6ddb8bc3a`.
- [x] Removed viscosity row scratch and multipass stress construction while
  preserving the three directional kernels in commit `64a641398`.
- [x] Added Spitzer activation, units, ceiling, and saturation validation in
  commit `8c706efd3`.
- [x] Removed Hydro STS full-state history copies and one history allocation,
  and removed update-kernel row scratch, in commit `5a312f18b`.
- [x] Added a globally reduced post-RK STS parabolic-bound refresh in commit
  `e7705ed81`.
- [x] Corrected the constant-viscosity longitudinal stability bound in commit
  `68cd36633`.
- [x] Made conduction timestep estimates use the same arithmetic face
  coefficients as the flux operator in commit `3b00ec548`.
- [x] Eliminated full STS flux-array clears for the single-level Hydro
  viscosity-first path in commit `4c52d3a91`; multilevel and scalar-diffusion
  paths retain the original clear.
- [x] Allowed a post-RK refresh to activate an STS sweep that was disabled by
  the pre-RK bound in commit `0d6800fb7`.
- [x] Added a combined constant viscosity-plus-conduction explicit/STS
  regression in commit `7803ef2e5`.
- [x] Added compile-time constant-viscosity and flux-update specialization in
  commit `96f239992`.
- [x] Replaced the scratch-free Hydro STS team update with a flat range kernel
  in commit `ad9961996`.
- [x] Passed all 22 focused CPU tests covering the RKL2 controller, Hydro STS,
  conduction timestep bounds, and shared hyperviscosity.
- [x] Built the final CCE 20.0.0, ROCm 6.4.2, Cray MPICH 9.0.1 HIP/MPI
  executable with the repository Frontier configuration.
- [x] Completed the operator-isolation matrix and baseline/optimized rocprof
  traces in job `4792663`.
- [x] Completed the final five-repetition interleaved two-node confirmation
  and matching rocprof traces in job `4793081`.
- [x] Validated explicit versus STS output using the exact final executable;
  maximum absolute state difference was `2.220446e-16`.
- [x] Produced the final CSV, PNG, and PDF comparison artifacts listed above.
- [!] The existing shared-STS MPI CPU regression cannot run on the Frontier
  login node because `mpirun` is unavailable there.  The completed two-node
  `srun` jobs provide the MPI/HIP validation for this target.
- [!] Temperature-dependent STS sweep-wide ceiling policy remains open.
  Power-law and Spitzer conduction have configured hard ceilings, but using
  those ceilings for every stage-selection bound is deferred to a dedicated
  correctness patch.  Multidimensional variable-viscosity stability remains
  unproved and is not broadened in this performance patch.
- [!] A configurable practical RKL2 stage cap remains a driver-level policy
  decision.  The numerical helper now rejects unrepresentable integer stage
  counts without embedding an arbitrary operational limit.
