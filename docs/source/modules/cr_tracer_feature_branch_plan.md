---
orphan: true
---

# CR Tracer Feature-Branch Hardening Plan

## Purpose

This plan covers improvements that fit naturally into the current
`feature/CR_tracers` branch without changing the particle architecture, AMR
ownership model, or MHD div(B) repair path.  The goal is to make the current
branch faster, easier to validate, and more useful to users while preserving the
behavior already covered by the regression suite.

Changes in this document should remain compatible with the existing CPU and
local MPI validation workflow.  They should not require GPU execution,
production-scale communication redesign, or a new particle data model.

## Scope

In scope:

- small performance fixes in existing particle kernels and host paths,
- stronger runtime guardrails that detect unsafe particle motion,
- additional low-cost diagnostics derived from fields already stored on
  particles,
- clearer validation commands and readback output,
- tests that extend the existing serial and two-rank MPI coverage,
- documentation updates for users and GitHub Pages integration.

Out of scope:

- device-side AMR remap,
- redesign of particle MPI exchange,
- new field-gather architecture,
- high-order CT-aware interpolation,
- particle-aware global load balancing,
- relativistic or electric-field pushers.

Those deeper changes belong in the follow-up architecture branch.

## Design Constraints

- Preserve the current `particle_type = cosmic_ray` input interface.
- Preserve existing `drift`, `boris + tsc`, and `boris + lin` behavior unless a
  new optional mode is explicitly selected.
- Keep `<particles>/check_consistency = false` as the production default.
- Keep all expensive checks opt-in or test-only.
- Do not modify `RefineFC()`, `RestrictFC()`, or `RepairAMRFC()`.
- Keep div(B) AMR regression inputs in the acceptance run.
- Keep new documentation in Markdown so it can be copied directly into the
  GitHub Pages branch.

## Phase 0: Baseline and Guardrails

### Objectives

Establish a measured baseline before making additional changes and make sure
future work does not weaken the already validated AMR behavior.

### Tasks

1. Record baseline timing for:
   - serial drift smoke,
   - serial Boris AMR with `tsc`,
   - serial Boris AMR with `lin`,
   - two-rank Boris AMR stress,
   - reduced `df` and `dxh` output.
2. Add a lightweight particle performance log line behind a test/debug option:
   - particles pushed,
   - particles sent,
   - particles received,
   - particles remapped after AMR,
   - histogram particles processed.
3. Keep the current inspector output stable so tests and users see the same
   count validation.
4. Add a short performance-baseline subsection to the particle docs with the
   exact commands used to collect timing.

### Acceptance Criteria

- Existing particle CPU and MPI suites still pass.
- Existing div(B) AMR compatibility inputs still pass.
- New timing/logging path is disabled by default.
- Documentation states that timings are local reference numbers, not portable
  performance guarantees.

## Phase 1: Low-Risk Performance Cleanup

### Objectives

Remove avoidable work that is already visible in the current implementation,
without changing algorithms.

### Tasks

1. Remove unused local cell-index calculations in the drift pusher path.
2. In host AMR remap, copy only data needed for destination discovery:
   - `IPX`,
   - `IPY`,
   - `IPZ`,
   - `PGID`.
3. Avoid full particle-array host mirrors in remap when only positions and
   ownership are needed.
4. Reuse temporary host buffers in remap when the existing object lifetime makes
   that practical and clean.
5. In output diagnostics, avoid repeated allocation when the existing output
   object can safely own a reusable histogram buffer.
6. Keep any reusable buffers rank-local and output-local to avoid cross-output
   state coupling.

### Acceptance Criteria

- Numerical outputs are unchanged except for non-semantic ordering or timing
  logs.
- `prst`, `ppd`, `df`, and `dxh` inspector checks still pass.
- AMR stress still reports both created and deleted MeshBlocks.
- No new production path depends on `<particles>/check_consistency = true`.

## Phase 2: Runtime Safety Checks for Motion and Ownership

### Objectives

Catch conditions that the current exchange path assumes are rare or prevented
by the timestep, especially large particle jumps across more than one
MeshBlock.

### Tasks

1. Add an optional `<particles>/check_motion_bounds` setting, defaulting off.
2. When enabled, check after each push that a particle did not move farther than
   the supported ownership-update stencil.
3. Report a clear error if a particle displacement can cross more than one
   MeshBlock in any direction before exchange.
4. Add a documentation note explaining that this is a guardrail, not a
   production subcycling method.
5. Extend consistency checks with a cheaper mode split if useful:
   - `none`,
   - `counts`,
   - `local`,
   - `full`.
6. Keep `check_consistency = true` as an alias for the full check to preserve
   current inputs.

### Acceptance Criteria

- Default production behavior is unchanged.
- Existing tests using `check_consistency = true` still run the full invariant
  set.
- A targeted test input can trigger the motion-bound failure in a controlled
  way.
- Error messages identify the particle species, tag, displacement, and allowed
  bound.

## Phase 3: Branch-Local Accuracy Improvements

### Objectives

Improve diagnostics and validation around numerical accuracy without replacing
the pusher or field-gather architecture.

### Tasks

1. Add an optional Boris uniform-field accuracy test:
   - known gyro period,
   - known gyro radius,
   - energy conservation check,
   - phase-error estimate.
2. Add a periodic-wrap displacement test:
   - force particles across periodic boundaries,
   - verify accumulated `IPDX`, `IPDY`, and `IPDZ` remain physically
     consistent,
   - verify `PGID` remains valid after exchange.
3. Add a local note to the docs clarifying the current meaning of
   `interpolation = lin` and `interpolation = tsc`.
4. Add a warning in docs that higher-order or CT-aware field gather is planned
   for the follow-up branch.

### Acceptance Criteria

- Uniform-field Boris test passes with a documented tolerance.
- Periodic-wrap displacement test passes in serial and MPI if practical.
- No change is made to the default pusher equations in this branch.

## Phase 4: Additional Low-Cost Diagnostics

### Objectives

Add useful diagnostics that can be computed from existing particle variables and
the already sampled magnetic field, without requiring a new field-gather layer.

### Tasks

1. Add reduced pitch-angle moments per species:
   - `<mu>`,
   - `<mu^2>`,
   - anisotropy proxy `3 <mu^2> - 1`.
2. Add displacement moments per species:
   - `<Delta x>`,
   - `<Delta y>`,
   - `<Delta z>`,
   - `<Delta x^2>`,
   - `<Delta y^2>`,
   - `<Delta z^2>`.
3. Add optional scalar displacement histogram `drh` for `|Delta x|`.
4. Add optional parallel-displacement histogram using existing `IPDB`.
5. Teach `cr_tracer_inspect.py` to summarize the new moments and validate
   their species totals when present.
6. Keep new outputs opt-in to avoid changing current output volume.

### Acceptance Criteria

- Existing output formats remain readable.
- New output files include enough header metadata for the inspector to validate
  bin counts and species counts.
- MPI-reduced diagnostics agree with particle counts to integer precision for
  count-like quantities.
- Documentation includes one minimal example of each new diagnostic.

## Phase 5: Documentation and GitHub Pages Integration

### Objectives

Keep the current feature branch easy to review and easy to publish to the docs
site.

### Tasks

1. Update `docs/source/modules/particles.md` with:
   - branch-local new options,
   - updated validation commands,
   - new diagnostic examples,
   - known limitations.
2. Keep `src/particles/README.md` as the implementation-oriented reference.
3. Add a short future-architecture link pointing to the deeper plan.
4. Avoid Sphinx-only syntax in the new planning pages unless already used by
   the GitHub Pages branch.
5. Rebuild the docs in a temporary GitHub Pages worktree before final commit if
   the Pages branch is available locally.

### Acceptance Criteria

- Markdown renders correctly in the local docs build or the Pages worktree.
- New docs can be copied into `docs/source/modules/` on `gh-pages`.
- Existing `particles.md` remains a user-facing how-to, not just a design log.

## Phase 6: Final Validation

### Required Commands

```bash
cd tst
python run_test_suite.py --test test_suite/particles/test_particles_cr_cpu.py --cpu
python run_test_suite.py --test test_suite/particles/test_particles_cr_mpicpu.py --mpicpu
python run_test_suite.py --style
git diff --check
```

Also rerun:

```bash
./athena -i inputs/tests/divb_amr_1d.athinput
./athena -i inputs/tests/divb_amr_2d.athinput
./athena -i inputs/tests/divb_amr_3d.athinput
```

### Acceptance Criteria

- All current particle tests pass.
- Any new particle tests pass in CPU and MPI modes.
- The moving AMR stress test still creates and deletes MeshBlocks.
- Total and per-species particle counts remain conserved.
- No invalid `PGID`, out-of-block position, duplicate tag, or malformed restart
  file is accepted when full consistency checks are enabled.
- div(B) AMR regressions remain within existing thresholds.
- Documentation includes exact validation results from the final run.

## Recommended Branch Exit Criteria

The current feature branch should be considered ready after this plan when:

- the branch has no known correctness regressions,
- branch-local performance cleanup is complete,
- guardrails detect large unsupported particle jumps,
- low-cost CR diagnostics are documented and tested,
- the docs clearly state which deeper improvements are intentionally deferred,
- the final commit and push include the updated validation record.
