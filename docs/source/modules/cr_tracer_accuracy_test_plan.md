---
orphan: true
---

# CR Tracer Accuracy Test Implementation Plan

## Purpose

This plan turns the proposed CR tracer accuracy ladder into an implementation
roadmap.  The goal is to demonstrate what the particle pusher, field-gather
methods, AMR remap, MPI exchange, and CR diagnostics are actually accurate for,
not only that they preserve particle counts.

The tests should produce both automated pass/fail regression checks and
publication-quality documentation artifacts.  Each test problem should generate
machine-readable outputs for CI and human-readable figures for the GitHub Pages
documentation.

## Scope

In scope:

- analytic and manufactured-field tests for the CR tracer pusher,
- resolution and timestep convergence studies,
- comparisons among `lin_legacy`, `trilinear`, and `tsc`,
- placeholders for future `pcs` and `ct_aware` methods,
- AMR and MPI decomposition invariance tests,
- quantitative diagnostics for pitch angle, energy, magnetic moment, drift, and
  distribution functions,
- qualitative trajectory, phase-space, and field-error figures for the docs.

Out of scope for the first implementation pass:

- changing the default interpolation method,
- adding `pcs` or `ct_aware` unless those methods have already landed,
- GPU performance or GPU-specific convergence claims,
- stochastic scattering physics beyond deterministic tracer motion,
- production-scale turbulence campaigns.

## Shared Test Infrastructure

### New Inputs

Add focused inputs under `inputs/particles/accuracy/`.  Keep each input small
enough for CI, then allow command-line overrides for larger documentation runs.

Recommended inputs:

- `cr_uniform_gyro.athinput`
- `cr_uniform_gyro_amr.athinput`
- `cr_linear_gather.athinput`
- `cr_manufactured_divb_free.athinput`
- `cr_smooth_orbit_reference.athinput`
- `cr_magnetic_mirror.athinput`
- `cr_gradb_curvature_drift.athinput`
- `cr_amr_boundary_convergence.athinput`
- `cr_mpi_decomposition_invariance.athinput`
- `cr_isotropic_ensemble.athinput`
- `cr_frozen_turbulent_field.athinput`

### New Test Helpers

Add a Python helper module, for example
`tst/test_suite/particles/cr_accuracy_utils.py`, with routines for:

- reading `prst`, `ppd`, `psamp`, `pmom`, `pspec`, and `pspec2`,
- matching particles by species and tag,
- comparing trajectories across runs,
- fitting convergence slopes from log-log error data,
- computing reference gyro-orbit and magnetic-mirror quantities,
- producing reusable JSON summaries for documentation plots.

The existing `scripts/particles/cr_tracer_inspect.py` should remain the
readback utility for users.  Test-only fitting code can live in the test suite,
but any generally useful plot or summary function should be placed under
`scripts/particles/`.

### Documentation Artifacts

Create a docs artifact directory:

```text
docs/source/modules/figures/cr_tracer_accuracy/
```

Each accuracy test should write:

- one JSON summary with exact numerical errors and run metadata,
- one or more PNG or SVG figures for the documentation,
- a short Markdown results block that can be pasted into the docs page.

The documentation page should report:

- run date,
- branch and commit,
- compiler, MPI, and Kokkos backend,
- input file and command-line overrides,
- tested interpolation methods,
- tested resolutions or timestep values,
- fitted convergence slopes,
- quantitative pass/fail tolerances,
- qualitative interpretation of what the figure shows.

## Common Metrics

Use these metrics consistently across the ladder.

| Metric | Definition | Used by |
|--------|------------|---------|
| `err_x` | Euclidean position error against analytic or reference solution. | Orbit, mirror, drift, AMR tests. |
| `err_v` | Euclidean velocity error against analytic or reference solution. | Orbit and reference tests. |
| `err_B_L1`, `err_B_L2`, `err_B_Linf` | Sampled magnetic-field error norms. | Gather-only and manufactured-field tests. |
| `speed_drift` | `(max(|v|) - min(|v|)) / mean(|v|)`. | Boris pusher tests. |
| `phase_error` | Gyro phase difference after a fixed number of periods. | Uniform and smooth orbit tests. |
| `mu_error` | Error in pitch-angle cosine `mu = v dot B / (|v||B|)`. | Mirror and ensemble tests. |
| `magnetic_moment_drift` | Drift in `v_perp^2 / B` relative to initial value. | Mirror and smooth-field tests. |
| `drift_error` | Error in mean perpendicular drift against analytic or reference value. | Grad-B and curvature tests. |
| `distribution_error` | Error in reduced histograms or moments. | Ensemble and turbulent-field tests. |
| `decomposition_delta` | Per-tag difference between one-rank and multi-rank results. | MPI invariance tests. |

Convergence fits should use at least three grid spacings or timesteps.  Four
points are preferred for documentation figures because they reveal pre-asymptotic
or saturation behavior.

## Figure Standards

Each test should produce at least one quantitative figure and, when useful, one
qualitative figure.

Quantitative figure types:

- log-log error versus `dt` with fitted slope,
- log-log error versus `dx` with fitted slope,
- invariant drift versus time,
- per-method bar chart of final error,
- table of fitted slopes and final errors,
- rank-decomposition difference histogram.

Qualitative figure types:

- particle trajectory overlaid on analytic trajectory,
- phase-space projection,
- pitch-angle distribution,
- field-error slice,
- AMR-level map with particle trajectory,
- histogram comparison for `df`, `pspec`, and `pspec2`.

Implemented documentation context figures:

- initial smooth and multi-mode prescribed magnetic-field slices, showing
  `delta Bz` and normalized in-plane perturbation streamlines,
- projected representative trajectories for one smooth-field particle and six
  multi-mode-field particles, with start/end markers.

All figure captions should state which error source is isolated.  For example,
do not caption a smooth-orbit result as an interpolation convergence result
unless the timestep error has already been made subdominant.

## Implementation Phases

### Phase 1: Clean Analytic Tests

Implement tests 1 through 4 first.  These isolate pusher and gather accuracy
without AMR, MPI, or long-time transport complications.

Acceptance gate:

- analytic gyro orbit passes,
- gather-only linear field identifies exact and non-exact interpolation paths,
- manufactured-field convergence slopes are measured and documented.

### Phase 2: Coupled Orbit and Invariant Tests

Implement tests 5 through 7 next.  These combine pusher and gather errors and
measure physically meaningful CR trajectory properties.

Acceptance gate:

- timestep and resolution errors can be separated,
- magnetic moment and drift diagnostics are finite and interpretable,
- method recommendations are supported by plots, not only by text.

### Phase 3: AMR, MPI, and Ensemble Tests

Implement tests 8 through 11 last.  These are more production-like and should
reuse the lower-level reference data from phases 1 and 2.

Acceptance gate:

- AMR boundaries do not destroy convergence,
- MPI decomposition changes do not change per-tag final states beyond tolerance,
- distribution diagnostics remain normalized and useful.

## Test 1: Uniform-B Single-Particle Gyro Orbit

### Purpose

Isolate the Boris pusher time accuracy in the simplest possible field.
Interpolation should not matter in this test.

### Setup

- Uniform magnetic field, for example `B = (0, 0, B0)`.
- One particle with nonzero perpendicular and parallel velocity.
- Periodic domain.
- No AMR for the first run.
- Run for an integer number of gyroperiods and for a non-integer number of
  gyroperiods.

### Sweeps

- `dt / tgyro = 1/8, 1/16, 1/32, 1/64`.
- `subcycle = false` and `subcycle = true`.
- `lin_legacy`, `trilinear`, and `tsc`, mainly to confirm they agree.

### Quantitative Outputs

- final `err_x` and `err_v`,
- `phase_error` after fixed time,
- `speed_drift`,
- fitted convergence slope of `phase_error` versus `dt`.

Expected behavior:

- `speed_drift` should be near roundoff,
- phase and position errors should converge approximately second order with
  timestep or subcycle size,
- all interpolation methods should agree to roundoff in uniform `B`.

### Figures

- Quantitative: log-log `phase_error` versus `dt / tgyro` with fitted slope.
- Qualitative: orbit projection in the perpendicular plane, showing analytic
  and numerical trajectories.

### Regression Criteria

- speed conservation below a strict tolerance,
- fitted phase-error slope above a minimum threshold, for example `1.8`,
- no method-dependent discrepancy in uniform `B`.

## Test 2: Uniform-B Mesh, AMR, and MPI Crossing Orbit

### Purpose

Show that MeshBlock ownership updates, AMR remap, MPI exchange, and periodic
wrapping do not introduce trajectory jumps.

### Setup

- Same analytic uniform-field orbit as test 1.
- Choose velocity, domain, and runtime so the particle crosses:
  - MeshBlock boundaries,
  - periodic boundaries,
  - MPI rank boundaries,
  - static or moving AMR boundaries.
- Compare to a no-AMR one-rank reference run with the same physical parameters.

### Sweeps

- one rank, two ranks, four ranks,
- no AMR, static AMR, moving AMR,
- low and high subcycle limits.

### Quantitative Outputs

- per-tag `err_x` and `err_v` versus the reference,
- maximum trajectory jump at ownership changes,
- count of sends, receives, and AMR remaps,
- final `phase_error`.

Expected behavior:

- convergence should match test 1 once timestep is controlled,
- errors should not spike at ownership changes,
- particle counts, species counts, and tags must be conserved.

### Figures

- Quantitative: final trajectory error versus timestep for no-AMR and AMR/MPI.
- Qualitative: trajectory over an AMR-level map, with boundary-crossing points
  marked.

### Regression Criteria

- no invalid `PGID`,
- no out-of-block positions after exchange or remap,
- AMR/MPI trajectory error remains within a small factor of the one-rank
  reference at the same timestep.

## Test 3: Gather-Only Linear Field Test

### Purpose

Isolate spatial field-gather accuracy without particle motion.

### Setup

- Place fixed particles at off-center positions in a linear field, for example:

```text
Bx = Bx0 + a*y
By = By0 + a*x
Bz = Bz0
```

- Use zero velocity or one-step runs that sample fields without meaningful
  orbit evolution.
- Compare sampled particle `B` fields against the analytic values.

### Sweeps

- resolutions `8^3`, `16^3`, `32^3`, `64^3`,
- particle positions near cell centers, faces, edges, and corners,
- all available interpolation modes.

### Quantitative Outputs

- `err_B_L1`, `err_B_L2`, and `err_B_Linf`,
- exactness flags for linear fields,
- per-component error tables.

Expected behavior:

- `trilinear` should be exact or near roundoff for a truly linear point-field
  setup,
- `lin_legacy` may show component-placement artifacts,
- `tsc` should remain consistent and smooth.

### Figures

- Quantitative: field-sample error versus resolution for each method.
- Qualitative: scatter plot of particle locations colored by field-error
  magnitude.

### Regression Criteria

- `trilinear` error below a near-roundoff tolerance for the manufactured linear
  field,
- no method returns non-finite sampled fields,
- documented legacy-method behavior remains stable.

## Test 4: Manufactured Smooth Divergence-Free Field

### Purpose

Measure pure spatial convergence in a smooth, divergence-free magnetic field.

### Setup

Define the magnetic field from a vector potential, for example:

```text
A = (A0 sin(k y) sin(k z),
     A0 sin(k z) sin(k x),
     A0 sin(k x) sin(k y))
B = B0 + curl(A)
```

This gives an analytic field with known point values.  Particles can be fixed
or moved negligibly so the test measures gather accuracy, not orbit error.

### Sweeps

- resolutions `16^3`, `32^3`, `64^3`, `128^3` for documentation runs,
- smaller `8^3`, `16^3`, `32^3` set for CI,
- all available interpolation methods.

### Quantitative Outputs

- `err_B_L1`, `err_B_L2`, `err_B_Linf`,
- fitted convergence slope for each method,
- per-component error norms.

Expected behavior:

- `trilinear` and `tsc` should show approximately second-order convergence in
  smooth fields,
- future `pcs` should exceed second order in the smooth-field regime,
- future `ct_aware` should be judged on both pointwise error and geometric
  consistency with face-centered fields.

### Figures

- Quantitative: log-log field-error norms versus `dx`, with fitted slopes.
- Qualitative: 2D slice of field-error magnitude for each method.

### Regression Criteria

- fitted slopes above method-specific minimums,
- no loss of convergence at ordinary cell boundaries,
- field-error summaries are written to JSON for docs.

## Test 5: Smooth-Field Orbit Against High-Resolution Reference

### Purpose

Combine pusher and gather errors in a nonuniform smooth field, while retaining a
controlled reference.

### Setup

- Use the manufactured divergence-free field from test 4.
- Launch particles with several pitch angles and gyro radii.
- Build a reference run using high spatial resolution and small subcycled
  timestep.

### Sweeps

- resolution sweep at fixed small timestep,
- timestep sweep at fixed high resolution,
- interpolation-method sweep.

### Quantitative Outputs

- final `err_x` and `err_v` against reference,
- `phase_error`,
- energy or speed conservation,
- pitch-angle error,
- magnetic moment proxy drift.

Expected behavior:

- timestep sweep should recover Boris/subcycling time convergence,
- resolution sweep should recover gather-method convergence until reference or
  timestep error dominates,
- errors should behave like `C_t dt^2 + C_x dx^p`.

### Figures

- Quantitative: two-panel convergence plot, one panel for timestep and one for
  resolution.
- Qualitative: representative 3D or projected trajectory against the
  high-resolution reference.

### Regression Criteria

- reference comparison is deterministic,
- fitted slopes are stable under excluding the coarsest point,
- no method produces secular speed drift above tolerance.

## Test 6: Divergence-Free Magnetic Mirror

### Purpose

Measure adiabatic invariant behavior and pitch-angle dynamics in a physically
recognizable CR problem.

### Setup

Use a smooth mirror-like field that is analytically divergence-free to leading
order, for example:

```text
Bz = B0 + a*z^2
Bx = -a*x*z
By = -a*y*z
```

Launch particles with several pitch angles, including trapped and passing
orbits.

### Sweeps

- pitch angles spanning mirror and non-mirror regimes,
- resolution sweep,
- subcycle/timestep sweep,
- interpolation-method sweep.

### Quantitative Outputs

- mirror point location error,
- bounce period error,
- `mu` evolution,
- magnetic moment proxy drift,
- trapped versus passing classification.

Expected behavior:

- better gather methods should reduce magnetic-moment drift,
- subcycling should reduce bounce-phase error,
- mirror classification should converge with resolution and timestep.

### Figures

- Quantitative: magnetic moment drift versus time for each method.
- Quantitative: mirror-point error versus resolution.
- Qualitative: `z(t)` and `mu(t)` curves for trapped and passing particles.

### Regression Criteria

- trapped particles remain trapped in converged runs,
- magnetic moment drift decreases with timestep and resolution,
- no interpolation method creates nonphysical speed growth.

## Test 7: Curvature and Grad-B Drift

### Purpose

Detect long-time transport bias that may not appear in one-orbit tests.

### Setup

Use a smooth field with controlled curvature or gradient strength.  The
preferred reference is a high-resolution full-orbit run, but where analytic
guiding-center drift estimates are reliable, compare to those as well.

### Sweeps

- weak, moderate, and strong field gradients,
- particle gyro radius relative to gradient scale,
- interpolation-method sweep,
- long runtime covering many gyroperiods.

### Quantitative Outputs

- mean perpendicular drift velocity,
- drift direction error,
- diffusion tensor proxy from displacement moments,
- long-time speed and magnetic-moment drift.

Expected behavior:

- high-order or smoother gather methods should reduce secular drift bias in
  smooth fields,
- noisy or inconsistent gathers can give acceptable short-time errors but poor
  long-time transport.

### Figures

- Quantitative: drift velocity error versus resolution.
- Quantitative: cumulative perpendicular displacement versus time.
- Qualitative: trajectory projected onto field-gradient coordinates.

### Regression Criteria

- drift sign and direction agree with reference,
- drift error decreases with resolution,
- long-time invariant drift stays within documented method-specific bounds.

## Test 8: AMR Refinement-Boundary Convergence

### Purpose

Demonstrate that AMR boundaries do not introduce first-order trajectory or
field-sampling artifacts.

### Setup

- Use a smooth manufactured field and orbit from tests 4 and 5.
- Run with:
  - uniform grid,
  - static AMR patch,
  - moving AMR patch,
  - refine plus derefine pattern.
- Force trajectories to cross coarse-fine boundaries repeatedly.

### Sweeps

- base resolution and maximum AMR level,
- static versus moving AMR,
- interpolation-method sweep,
- one-rank and multi-rank configurations.

### Quantitative Outputs

- trajectory error binned by distance to AMR boundary,
- field-sample error near and far from AMR boundaries,
- global convergence against uniform high-resolution reference,
- particle count and tag conservation.

Expected behavior:

- AMR should not change the asymptotic convergence order,
- error may be larger near refinement boundaries but should remain bounded and
  decrease with resolution,
- no jumps should occur when particles cross refinement boundaries.

### Figures

- Quantitative: near-boundary and far-from-boundary error versus resolution.
- Qualitative: AMR level map with particle trajectory and error markers.

### Regression Criteria

- refinement and derefinement both occur in the stress version,
- total and per-species particle counts are conserved,
- near-boundary error does not show first-order convergence unless documented as
  an expected limitation of a method.

## Test 9: MPI Decomposition Invariance

### Purpose

Prove that parallel decomposition does not change particle trajectories or
diagnostics beyond documented floating-point tolerances.

### Setup

- Use deterministic particle initialization and fixed AMR schedule.
- Run identical physical setup with one, two, and four MPI ranks.
- Sort particles by species and tag before comparison.

### Sweeps

- no AMR and AMR,
- low and high migration fraction,
- one species and two species,
- all supported interpolation methods.

### Quantitative Outputs

- per-tag `err_x`, `err_v`, and sampled-field differences,
- rank-count difference in reduced diagnostics,
- send/receive/remap counts,
- histogram of decomposition deltas.

Expected behavior:

- one-, two-, and four-rank runs should agree within tolerance,
- reduced diagnostics should be identical for integer histograms and
  close-to-identical for floating moments,
- particle rank ownership can differ, but physical state cannot.

### Figures

- Quantitative: histogram of per-particle differences by rank count.
- Quantitative: table of max and RMS decomposition deltas.
- Qualitative: rank ownership map or particles-per-rank plot.

### Regression Criteria

- no lost or duplicated tags,
- reduced histogram sums match expected species counts,
- max decomposition delta remains below a stated tolerance.

## Test 10: Distribution-Level Isotropic Ensemble

### Purpose

Validate CR ensemble diagnostics and show that the method preserves a simple
stationary distribution.

### Setup

- Many particles in uniform or weakly varying magnetic field.
- Fixed speed and isotropic pitch-angle distribution.
- Multiple species with known counts.

### Sweeps

- particle count,
- species count,
- histogram bin count,
- interpolation method,
- one-rank and multi-rank reductions.

### Quantitative Outputs

- `f(mu,p)`, `f(mu,E)`, and `f(v_parallel,v_perp)` normalization,
- `<mu>`, `<mu^2>`, anisotropy,
- pressure tensor proxy,
- flux/current proxy,
- diffusion tensor proxy for controlled displacement cases.

Expected behavior:

- for an isotropic ensemble, `<mu>` should be near zero and `<mu^2>` near
  `1/3`,
- reduced integer histograms should sum exactly to species counts,
- MPI reductions should match serial totals.

### Figures

- Quantitative: histogram residuals relative to expected isotropic distribution.
- Quantitative: table of moments by species and rank count.
- Qualitative: 2D `pspec2` image for `f(mu,p)` or
  `f(v_parallel,v_perp)`.

### Regression Criteria

- histogram counts are exactly conserved,
- moment errors fall with particle count as sampling noise permits,
- selected-field `psamp` rows are deterministic by species and tag.

## Test 11: Frozen Turbulent-Field Reference Test

### Purpose

Provide a production-like accuracy demonstration where method recommendations
can be made honestly.

### Setup

- Use a prescribed smooth turbulent magnetic field or a frozen MHD snapshot.
- Freeze the fluid state and evolve tracer particles through the field.
- Compare lower-resolution runs against a high-resolution reference or a
  filtered-reference hierarchy.

### Sweeps

- grid resolution,
- particle count,
- interpolation method,
- subcycle setting,
- optional AMR versus uniform grid.

### Quantitative Outputs

- pitch-angle evolution,
- displacement moments,
- diffusion tensor proxy,
- `df`, `pspec`, and `pspec2` differences from reference,
- method cost versus error.

Expected behavior:

- `tsc` may suppress grid-scale noise in ensemble diagnostics,
- `trilinear` may preserve sharper local structure,
- higher-order methods should help only where the field is smooth enough and
  the wider stencil is well behaved near AMR boundaries.

### Figures

- Quantitative: cost-error plot for each method.
- Quantitative: diffusion tensor components versus time.
- Qualitative: particle trajectories over a magnetic-field slice.
- Qualitative: side-by-side `pspec2` images for each method.

### Regression Criteria

- documentation run establishes method tradeoffs,
- CI version checks only reduced, deterministic quantities,
- no method recommendation is made without both error and cost evidence.

## Method Recommendation Logic

The docs should present recommendations as a decision table, not as a single
winner.

| Method | Best use | Accuracy expectation | Cost and risk |
|--------|----------|----------------------|---------------|
| `lin_legacy` | Backward comparison and regression. | Stable legacy behavior, not the preferred accuracy path. | Lowest conceptual migration cost. |
| `trilinear` | Smooth pointwise field sampling. | Second-order smooth-field convergence when time error is controlled. | Cheap and easy to interpret. |
| `tsc` | Ensemble transport where smoothness matters. | Good noise suppression and stable histogram behavior. | Wider stencil and more smoothing. |
| future `pcs` | Very smooth fields needing higher spatial order. | Expected higher-than-second-order convergence in smooth regions. | Wider stencil, AMR-boundary complexity. |
| future `ct_aware` | Geometric consistency with face-centered MHD fields. | Should improve behavior tied to constrained-transport geometry. | More complex implementation and validation. |
| subcycling | High gyro frequency or large cell-crossing distance. | Reduces time and ownership errors; should recover pusher convergence. | More particle pushes and exchanges. |

## Proposed Test Files

Split tests into CI-sized tests and documentation-sized tests.

CI tests:

- `tst/test_suite/particles/test_particles_cr_accuracy_cpu.py`
- `tst/test_suite/particles/test_particles_cr_accuracy_mpicpu.py`

Documentation driver:

- `scripts/particles/cr_tracer_accuracy_study.py`

The documentation driver should run larger sweeps, write JSON summaries, and
make figures.  The CI tests should run small versions that assert slopes and
invariants quickly.

## Documentation Page Structure

Add a final user-facing page after the tests are implemented:

```text
docs/source/modules/cr_tracer_accuracy.md
```

Suggested sections:

1. What Accuracy Is Being Tested
2. Summary Recommendation Table
3. Pusher Time Accuracy
4. Field-Gather Spatial Accuracy
5. AMR and MPI Accuracy
6. Distribution and Transport Diagnostics
7. Method Selection Guidance
8. Reproducibility Commands
9. Known Limits and GPU Follow-Up

The page should embed figures from
`docs/source/modules/figures/cr_tracer_accuracy/` and link back to this plan.

## Acceptance Criteria For The Full Accuracy Suite

The full suite is complete when:

- all 11 test problems have small CI versions,
- at least tests 1 through 5 have documented convergence figures,
- AMR and MPI tests compare per-tag states across decompositions,
- distribution tests validate `df`, `pspec`, `pspec2`, and `pmom`,
- every figure is regenerated by a documented command,
- the docs include both qualitative interpretation and quantitative slopes,
- recommendations for `lin_legacy`, `trilinear`, `tsc`, and subcycling are
  based on the measured results,
- future `pcs` and `ct_aware` rows are present but clearly marked TODO until
  those methods exist,
- the docs can be copied into `docs/source/modules/` on `origin/gh-pages`
  without changing the main `particles.md` toctree entry.
