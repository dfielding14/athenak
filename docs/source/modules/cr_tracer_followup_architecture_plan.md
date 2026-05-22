---
orphan: true
---

# CR Tracer Follow-Up Architecture Plan

## Purpose

This plan covers deeper computational and numerical improvements that should be
implemented in a follow-up branch after `feature/CR_tracers` is accepted.  These
changes are larger because they touch core particle architecture, communication
scaling, field interpolation policy, and possibly the load balancer.

The goal is to turn the correctness-first CR tracer feature into a production
path suitable for large AMR, MPI, and eventually GPU runs.

## Scope

In scope:

- device-side AMR particle remap,
- scalable particle MPI exchange,
- device-side particle compaction,
- modular field-gather policies,
- optional higher-order interpolation,
- optional particle subcycling,
- richer CR physics diagnostics,
- particle-aware load balancing,
- performance benchmarks at increasing particle and rank counts.

Out of scope:

- changing the accepted behavior of the original `feature/CR_tracers` tests
  without explicit migration notes,
- modifying the div(B) AMR repair algorithm unless a particle change exposes a
  genuine integration bug,
- requiring all users to adopt higher-order methods by default,
- replacing the current particle restart format without a compatibility reader.

## Architectural Principles

- Keep low-cost default methods available.
- Add higher-order and higher-cost methods as explicit runtime options.
- Separate pusher logic from field-gather logic.
- Keep AMR remap independent from MHD face-centered repair.
- Prefer device-resident particle operations for production paths.
- Keep debug consistency checks separate from production kernels.
- Preserve restart and output compatibility where possible.

## Phase 0: Benchmarks and Scaling Targets

### Objectives

Quantify the bottlenecks before replacing core paths.

### Tasks

1. Add particle microbenchmarks for:
   - push throughput,
   - AMR remap throughput,
   - particle MPI exchange throughput,
   - histogram output throughput,
   - restart write/read throughput.
2. Sweep:
   - particles per rank,
   - MeshBlocks per rank,
   - AMR level count,
   - migration fraction,
   - species count,
   - histogram bin count.
3. Record metrics:
   - particles pushed per second,
   - particles remapped per second,
   - bytes exchanged per particle,
   - MPI messages per exchange,
   - output time per particle,
   - host-device copy volume.
4. Add a benchmark documentation page with machine, compiler, MPI, and Kokkos
   configuration.

### Acceptance Criteria

- Baseline measurements identify the dominant cost at small, medium, and large
  particle counts.
- Benchmarks can be run without changing production inputs.
- Results include enough metadata to compare across machines.

## Phase 1: Device-Side AMR Remap

### Objectives

Remove the host loop and full particle mirror copies from AMR remapping.

### Design

Use a compact device-readable leaf lookup table.  Each live leaf MeshBlock is
represented by:

- finest-level logical lower index,
- finest-level logical upper index,
- `gid`,
- owning rank,
- optional local MeshBlock index.

A particle position is wrapped, converted to finest-level logical coordinates,
and mapped to a leaf by either:

- binary search over sorted logical keys,
- a hash table keyed by finest-level cell or block prefix,
- a flattened mesh-tree representation suitable for device traversal.

### Tasks

1. Define a device lookup structure owned by `Mesh` or a small AMR lookup
   helper.
2. Rebuild the lookup only after AMR changes the leaf list.
3. Implement a device kernel that updates particle `PGID` from position.
4. Preserve the existing MPI send path after `PGID` reassignment.
5. Add a debug-only validation path that compares device lookup against the
   current host tree lookup for selected particles.
6. Add tests with refinement, derefinement, and periodic boundaries.

### Acceptance Criteria

- AMR remap no longer copies all particle arrays to host.
- Device and host lookup agree in debug validation.
- Current AMR stress tests pass unchanged.
- Remap time scales approximately with particle count, not particle count times
  MeshBlock count.

## Phase 2: Scalable Particle MPI Exchange

### Objectives

Replace global metadata collection and host-heavy exchange bookkeeping with a
scalable sparse communication path.

### Design Options

Option A: sorted sparse point-to-point exchange.

- Sort or bin outgoing particles by destination rank.
- Exchange send counts only with ranks that have traffic.
- Use nonblocking sends and receives for packed particle buffers.

Option B: neighborhood collectives.

- Build a communication graph from neighboring ranks or recent particle traffic.
- Use `MPI_Neighbor_alltoallv` where the MPI implementation supports it well.

Option C: rank-count all-to-all.

- Use an all-to-all count exchange but avoid all-gathering detailed send tuples.
- This is simpler and may be sufficient for moderate rank counts.

### Tasks

1. Choose one exchange design based on benchmark results.
2. Pack particle real and integer fields into fewer messages per destination.
3. Keep the old exchange path under a compile-time or runtime fallback until
   validation is complete.
4. Add stress tests with sparse and dense particle migration.
5. Track migration counts and bytes exchanged in diagnostics.

### Acceptance Criteria

- Exchange metadata cost does not grow as a dominant all-rank bottleneck for
  sparse migration.
- Particle counts and tags remain conserved under high migration.
- The fallback path can be used to bisect communication regressions.

## Phase 3: Device-Side Compaction and Buffer Management

### Objectives

Remove serial host-side particle removal and reduce allocation churn during
exchange and diagnostics.

### Tasks

1. Implement device-side compaction after sends:
   - survivor mask,
   - prefix sum,
   - scatter survivors to compacted slots.
2. Add buffer capacity management for receive growth.
3. Reuse send, receive, and diagnostic buffers across exchanges where safe.
4. Add tests with:
   - no sends,
   - all particles sent,
   - mixed sends and receives,
   - receive count larger than send count,
   - send count larger than receive count.

### Acceptance Criteria

- Exchange no longer uses a serial host loop to fill particle holes.
- Particle order is documented as unstable unless a test explicitly requires
  otherwise.
- Reallocation occurs only when capacity is insufficient.

## Phase 4: Modular Field Gather and Optional Higher Order Methods

### Objectives

Separate the Boris pusher from magnetic-field interpolation so new methods can
be added without duplicating pusher logic.

### Runtime Interface

Keep:

```text
<particles>
pusher        = boris
interpolation = tsc
```

Add optional values such as:

```text
interpolation = trilinear
interpolation = pcs
interpolation = ct_aware
```

### Methods

1. `lin_legacy`:
   - preserves the current face-component interpolation behavior,
   - useful for backward comparisons.
2. `trilinear`:
   - gathers a consistent point-centered magnetic field with full 3D
     interpolation.
3. `tsc`:
   - keeps current triangular-shaped-cloud gather.
4. `pcs`:
   - optional piecewise-cubic spline gather for smoother fields.
5. `ct_aware`:
   - uses the face-centered CT representation to gather a field that is more
     geometrically consistent with the divergence-controlled MHD state.

### Tasks

1. Define a `FieldGather` policy interface.
2. Move current `lin` and `tsc` gather code behind that interface.
3. Add `trilinear` as the first new non-legacy method.
4. Add convergence tests in analytic magnetic fields.
5. Add a uniform-field orbit test that compares interpolation methods.
6. Add documentation describing cost, accuracy, and recommended use.

### Acceptance Criteria

- Existing `lin` and `tsc` results can be reproduced.
- New `trilinear` mode has demonstrably better smooth-field convergence than
  the legacy linear gather.
- The pusher update is not duplicated for each gather method.

## Phase 5: Particle Subcycling and Time-Centered Fields

### Objectives

Improve orbit accuracy and guarantee that particle motion remains compatible
with ownership updates.

### Tasks

1. Add particle subcycling based on:
   - gyro period,
   - cell crossing fraction,
   - MeshBlock crossing fraction,
   - user maximum substeps.
2. Add a timestep diagnostic reporting the active particle constraint.
3. Add optional midpoint or predictor-corrector field sampling.
4. Explore time interpolation of MHD fields if old/new states are available in
   the task flow.
5. Keep the default mode simple until convergence tests justify a new default.

### Acceptance Criteria

- Large particle velocities no longer violate ownership-update assumptions.
- Uniform-field orbit phase error decreases with substep count.
- The method remains stable across AMR and MPI exchange.
- Runtime failure modes are clear when the requested subcycle limit is too low.

## Phase 6: Richer CR Diagnostics and Moment Outputs

### Objectives

Make the tracer particles useful for CR transport analysis, not only position
tracking.

### Diagnostics

Add optional reduced diagnostics for:

- `dN/dp`,
- `dN/dE`,
- `dN/dlogE`,
- `f(mu, p)`,
- `f(mu, E)`,
- `f(v_parallel, v_perp)`,
- `<mu>`, `<mu^2>`, and anisotropy,
- parallel and perpendicular displacement moments,
- diffusion tensor estimates,
- magnetic moment proxy `v_perp^2 / B`,
- CR number density deposited to mesh,
- CR current or flux,
- CR pressure tensor,
- particles per MeshBlock and particles per rank.

### Tasks

1. Define output metadata for multidimensional histograms.
2. Add inspector support for each new diagnostic family.
3. Add deterministic sample-output mode selected by species and tag hash.
4. Add tests that check normalization, finite values, and MPI reduction.
5. Document which diagnostics use stored pusher fields and which resample the
   current MHD fields.

### Acceptance Criteria

- Histogram and moment sums are validated by the inspector.
- Diagnostics are globally reduced by default where that is the documented
  behavior.
- Output volume can be controlled by user options.
- Tests cover both one-species and multi-species cases.

## Phase 7: Particle-Aware Load Balancing

### Objectives

Prevent CR tracer cost from becoming concentrated on a small number of ranks in
AMR runs.

### Tasks

1. Add particle count to MeshBlock cost estimates.
2. Add runtime parameters for particle load weight.
3. Track particles per MeshBlock and per rank before and after rebalancing.
4. Add stress tests where AMR refinement and particle concentration move across
   the domain.
5. Compare mesh-only and particle-aware balancing in benchmark outputs.

### Acceptance Criteria

- Particle-rich regions are distributed more evenly when particle weighting is
  enabled.
- Mesh-only load balancing remains available.
- Particle migration caused by rebalancing preserves total and per-species
  counts.

## Phase 8: Restart and Output Format Robustness

### Objectives

Make particle outputs easier to parse, validate, and extend.

### Tasks

1. Add explicit binary record metadata where missing:
   - magic value,
   - version,
   - real size,
   - integer size,
   - species count,
   - field count,
   - particle count.
2. Keep a compatibility reader for current `prst`, `ppd`, `df`, and `dxh`
   outputs.
3. Add optional selected-field particle dumps.
4. Add checksums or record-length validation for restart files.
5. Update `cr_tracer_inspect.py` to report format versions.

### Acceptance Criteria

- Current branch outputs remain readable.
- New outputs are self-describing enough for independent analysis scripts.
- Restart read errors report the exact malformed record or length mismatch.

## Phase 9: Final Validation and Production Readiness

### Required Validation

1. Existing CR tracer CPU suite.
2. Existing CR tracer MPI CPU suite.
3. New architecture benchmarks.
4. New field-gather convergence tests.
5. New subcycling orbit tests.
6. MPI AMR stress at more than two ranks.
7. GPU build and smoke test if the target production path includes GPUs.
8. div(B) AMR 1D, 2D, and 3D compatibility inputs.
9. Documentation build for the GitHub Pages branch.

### Acceptance Criteria

- Production path avoids host-side AMR remap.
- MPI exchange no longer relies on global detailed send metadata.
- Particle push, remap, exchange, and diagnostics have benchmark coverage.
- Higher-order methods are optional and documented with cost/accuracy guidance.
- Existing low-cost methods remain available for comparison and regression.
- Documentation includes migration notes from `feature/CR_tracers`.

## Recommended Branch Exit Criteria

The follow-up architecture branch should be considered complete when:

- the device-side remap is the default production AMR path,
- the old remap remains available only as a debug or fallback path,
- particle exchange scales beyond two-rank testing,
- at least one improved interpolation method is validated,
- subcycling or an equivalent motion bound protects ownership updates,
- CR diagnostics support transport analysis beyond raw position dumps,
- docs and inspector utilities make the new outputs straightforward to use.
