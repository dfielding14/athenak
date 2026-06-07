# Ito Tracer Review Remediation, Verification, and Test Guide

## Purpose

This document turns the review of `feature/ito_tracers` into an executable
remediation plan. It covers:

- confirmed correctness defects;
- numerical and physics-model limitations;
- missing regression and validation coverage;
- performance and scaling improvements;
- documentation corrections; and
- unrelated branch issues found during the review.

The reviewed baseline is commit
`b25a406c0fdc60a051e53873900e24580c3d688f`. The implementation is based on
Moseley, Teyssier, and Abel,
[arXiv:2604.23041](https://arxiv.org/abs/2604.23041).

The most important distinction in this guide is:

> The current implementation matches the Monte Carlo jump mean and variance
> separately in each coordinate. It does not match the full multidimensional
> finite-step covariance tensor.

No documentation, test name, or scientific claim should blur that distinction.

## Completion Definition

The remediation is complete only when all of the following are true:

1. The intended mathematical contract is written explicitly and tested.
2. The full finite-step displacement covariance is either implemented or
   explicitly excluded from the supported claim.
3. Particle tags remain unique and deterministic past the signed 32-bit range.
4. AMR and Runge-Kutta behavior are validated with moment-level tests rather
   than smoke tests.
5. Restarted trajectories match uninterrupted trajectories.
6. Boundary, high-CFL, invalid-input, single-precision, MPI, and GPU cases have
   explicit outcomes.
7. Performance changes are measured against a recorded baseline.
8. The repository style gate passes.
9. User documentation reflects the tested behavior and is validated against
   the live `origin/gh-pages` documentation tree.

## Priority Summary

| Priority | Finding | Disposition |
| --- | --- | --- |
| P0 | Multidimensional finite-step covariance is omitted | Fix or narrow all claims before scientific use |
| P0 | 64-bit next-tag counter is truncated into 32-bit particle tags | Convert tag storage and I/O to 64-bit |
| P1 | AMR transfer acts on nonlinear derived coefficients | Establish and test a level-consistent moment transfer |
| P1 | RK-weighted net flux is converted to probabilities after stage summation | Define intended stage semantics and add reversing-flux tests |
| P1 | Existing tests do not establish paper-level fidelity | Add covariance, distribution, restart, AMR, and turbulence tests |
| P1 | Exact upper periodic boundary is not wrapped | Use half-open periodic bounds and add an exact-boundary test |
| P1 | Mass-changing source terms are not represented by flux tracers | Reject, document, or implement source coupling |
| P2 | Coefficient construction and communication have avoidable cost | Optimize only after correctness changes are frozen |
| P2 | Seeding and history output are host-heavy and root-gathered | Add scalable paths and benchmarks |
| P2 | Thermodynamic history uses containing-cell sampling | Add optional CIC sampling |
| P2 | Ito-2 cannot reproduce MC distribution shape in general | Add Ito-3 or retain an explicit limitation |
| P0 branch gate | C++ and Python style checks fail in cooling code | Repair before merge |
| P1 branch gate | Single-precision build fails in existing coordinate code | Fix or explicitly exclude single precision |
| P2 docs | Ito page is not present on live `origin/gh-pages` | Integrate and strictly build the live docs tree |
| P2 integration | Branch combines Ito with cooling, perturbation, divB, and MC work | Split or validate dependencies deliberately |

## Observed Review Baseline

The following results were reproduced on the reviewed branch tip.

Passed:

- Release serial build;
- Release MPI build;
- two serial Ito CPU tests;
- the two-rank MPI/AMR Ito test;
- thermodynamic-history reader test;
- `git diff --check`; and
- repository connectivity check.

Failed or not established:

- the repository style gate fails in branch-added cooling files;
- the repository does not build in single precision with AppleClang because of
  pre-existing narrowing conversions;
- no GPU run was performed;
- no paper-equivalent square-pulse or 3D turbulence validation exists; and
- the Ito page is not published on the live `origin/gh-pages` branch.

Independent reproductions:

- a 2D diagonal-flow run matched both marginal variances but measured
  correlation `-0.000278`, compared with the finite-step MC expectation
  `-0.111111`;
- a two-particle seed beginning at `INT_MAX` produced tags
  `[2147483647, -2147483648]`.

## Phase 0: Freeze the Baseline

Before changing behavior, record the current branch state and preserve
reproducible baselines.

### Required records

Record:

- branch and commit SHA;
- compiler, CMake flags, Kokkos backend, MPI implementation, and precision;
- serial and MPI test results;
- current one-step 1D, 2D, and 3D moment measurements;
- current AMR smoke-test output;
- particle updates per second for representative particle counts; and
- memory use for saved fluxes and Ito coefficient arrays.

### Baseline commands

```bash
git status --short --branch
git rev-parse HEAD
git diff --check

cmake -S . -B /tmp/ito-serial -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/ito-serial -j 8

cmake -S . -B /tmp/ito-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=ON
cmake --build /tmp/ito-mpi -j 8
```

The repository test runner assumes it is launched from `tst/`:

```bash
cd tst
python run_test_suite.py \
  --cpu \
  --test test_suite/particles/test_particles_ito_cpu.py

python run_test_suite.py \
  --mpicpu \
  --test test_suite/particles/test_particles_ito_mpicpu.py
```

### Baseline acceptance

- Builds and existing Ito tests must pass before behavioral changes.
- New tests that expose known defects should initially fail for the expected
  reason.
- Baseline performance results must be stored in a small machine-readable
  table so later optimizations can be evaluated.

## Phase 1: Correct the Multidimensional Covariance

### Exact paper locations

The relevant passage is Section 2.2, printed page 7 of
[arXiv:2604.23041v2](https://arxiv.org/abs/2604.23041v2):

| Paper location | Assessment |
| --- | --- |
| Equations (28)-(29) | Correct 1D drift and central variance. The `-C_-^2` term is the diagonal mean-product subtraction. |
| Equation (45) | Uses `E[dX_i dX_j]/(2 dt)` while calling the result the covariance matrix. This is the raw second-moment rate, not the finite-step central covariance. |
| Equation (46) | Correctly observes that the raw MC cross moment is zero, but does not subtract `E[dX_i] E[dX_j]`. |
| Text after equation (46) | The conclusion that the dimensions are independent is not true for the finite-step MC displacement when more than one mean component is nonzero. |
| Equations (47)-(48) | The vector SDE and drift can remain unchanged. |
| Equation (49) | For a non-diagonal tensor, replace the componentwise square root with a factor satisfying `Sigma Sigma^T = 2 K`. |
| Equation (59) and footnote 4 | The finite-step update can remain unchanged once `Sigma` is corrected. The choice `delta t = Delta t` is why the omitted finite-step cross covariance does not vanish in the implemented step. |
| Equations (60)-(61) | The stay-local estimate should be re-derived if it is retained for a non-diagonal diffusion tensor. |

Important qualification: equations (45)-(46) are consistent with the
diagonal **raw** second Kramers-Moyal coefficient in the strict
`Delta t -> 0` limit. The defect arises when that quantity is described as a
finite-step covariance and used to support full finite-step multidimensional
moment-matching claims.

### Mathematical contract

For one finite Monte Carlo tracer step, let coordinate `i` have left and right
jump probabilities `p_i-` and `p_i+`, cell spacing `h_i`, and

```text
C_i+ = p_i+ + p_i-
C_i- = p_i+ - p_i-
```

Only one of all possible face jumps can occur during the MC step. Therefore,

```text
E[dX_i]       = h_i C_i-
E[dX_i^2]     = h_i^2 C_i+
E[dX_i dX_j]  = 0                         for i != j
```

The finite-step central covariance is consequently

```text
Cov[dX_i,dX_i] = h_i^2 (C_i+ - C_i-^2)
Cov[dX_i,dX_j] = -h_i h_j C_i- C_j-       for i != j
```

The current diagonal diffusion coefficients implement only the first line.
Independent coordinate kicks force the second line to zero.

### Recommended implementation

Represent the full symmetric covariance or diffusion tensor:

```text
Q11 Q22 Q33 Q12 Q13 Q23
```

where `Q = Cov[dX]` for one fluid step, or equivalently
`K = Q/(2 dt)` if the code retains diffusion-rate units.

The coefficient field then contains:

```text
U1 U2 U3 Q11 Q22 Q33 Q12 Q13 Q23
```

Recommended code targets:

- `src/particles/particles.hpp`
- `src/particles/particles.cpp`
- `src/particles/particles_lagrangian_ito.cpp`
- the Ito coefficient boundary-communication calls
- CPU and MPI particle tests

At the particle position:

1. CIC-interpolate the drift and six independent covariance components.
2. Symmetrize the interpolated matrix explicitly.
3. Verify finiteness and positive semidefiniteness within precision-scaled
   tolerance.
4. Compute a robust 2x2 or 3x3 factor `L` with `Q = L L^T`.
5. Draw a vector of independent zero-mean, unit-variance bounded variables.
6. Advance with:

```text
dX = U dt + L xi
```

A Cholesky-like factorization must handle positive-semidefinite matrices with
zero pivots. Negative eigenvalues or pivots larger than roundoff tolerance
must fail closed. Tiny negative values caused by roundoff may be clipped only
after recording the tolerance policy in code and tests.

### Important implementation choice

There are two defensible contracts:

1. **Full finite-step moment matching.** Implement the off-diagonal central
   covariance above.
2. **Per-coordinate marginal matching only.** Keep independent kicks, but
   rename the method and narrow every claim and test accordingly.

The first contract is recommended if the feature is described as matching the
first two multidimensional displacement moments of the MC kernel.

### Required covariance tests

Add deterministic statistical tests with at least 65,536 particles for
multidimensional cases:

| Case | Required checks |
| --- | --- |
| 2D diagonal positive flow | both means, both variances, negative cross covariance |
| 2D mixed-sign flow | covariance sign follows `-mean_x mean_y` |
| 2D one stationary axis | cross covariance is zero |
| 3D flow | all 3 means, 3 variances, and 3 cross covariances |
| zero flow | exactly zero drift and covariance |
| outgoing probability near one | factorization remains stable |
| rank-deficient covariance | no NaNs and expected zero stochastic direction |
| single precision | tolerance-scaled covariance agreement |

For each covariance component, use a sampling-error-based tolerance, not a
fixed percentage alone. The test should report the expected value, measured
value, standard error, and normalized residual.

### Covariance acceptance

- All means and covariance components agree with the analytical finite-step MC
  kernel within five estimated standard errors.
- No materially negative covariance eigenvalue survives validation.
- Results are deterministic for fixed tags, cycle, and random seed.
- Serial and MPI decompositions give identical per-tag trajectories.

## Phase 2: Make Tracer Tags Truly 64-bit

### Confirmed defect

`next_tracer_tag` is 64-bit, but `PTAG` is stored in
`DvceArray2D<int> prtcl_idata`. The seeding path explicitly casts the next tag
to `int`. Tags become negative immediately after `INT_MAX` and eventually
repeat after the 32-bit space wraps.

### Recommended representation

Use a dedicated unsigned 64-bit tag array:

```cpp
DvceArray1D<std::uint64_t> prtcl_tag;
```

Keeping GID, seed ID, and small state fields in the existing integer array
limits the blast radius and avoids doubling every integer particle property.

Update all paths that create, move, resize, serialize, or output particles:

- initial tag creation;
- scheduled seeding;
- append/reallocation;
- MPI packing and unpacking;
- AMR remapping;
- restart write/read;
- VTK and tracked-particle output;
- thermodynamic-history output;
- RNG key construction; and
- Python history readers.

### File-format policy

- Increment the particle restart format version.
- Increment the thermodynamic-history format version.
- Store tags as explicit `uint64`.
- Preserve readers for existing restart/history versions.
- Reject unsupported narrowing conversions rather than silently truncating.
- Add a 64-bit parameter-input path if `next_tracer_tag` remains user-settable.

### Required tag tests

Test tags around:

```text
INT_MAX - 1
INT_MAX
INT_MAX + 1
UINT32_MAX
UINT32_MAX + 1
2^63
UINT64_MAX - 1
```

Verify:

- uniqueness;
- exact restart round trip;
- exact history-reader values;
- MPI migration;
- AMR remapping;
- deterministic RNG streams;
- no collision between old and newly seeded particles; and
- a precise fail-closed error instead of unsigned wrap at `UINT64_MAX`.

## Phase 3: Define Runge-Kutta Flux Semantics

### Open numerical issue

The code sums final-RK-weighted signed stage fluxes and then applies
`max(flux, 0)` when constructing jump probabilities.

Because:

```text
max(sum(stage_flux), 0) != sum(max(stage_flux, 0))
```

opposing stage fluxes can cancel before the probability moments are computed.
This exactly represents the signed mass change of the final finite-volume
update, but it is not generally the same stochastic process as composing
stage-wise MC transition kernels.

### Required decision

Document which object the tracer is intended to match:

- the net final finite-volume mass transfer over the full timestep; or
- the composition of the stage-level mass-transfer kernels.

Do not change the algorithm until this contract is settled. A stage-wise
implementation can be wrong if low-storage RK stage states and weights are
treated as ordinary positive substeps.

### Required RK tests

Add manufactured or instrumented tests with:

- constant-sign stage fluxes;
- sign reversal between stages;
- near cancellation;
- compression with outward flux through multiple faces;
- RK1, RK2, and RK3;
- FOFC activation; and
- AMR flux correction.

Record both:

- moments derived from the current net-flux kernel; and
- moments from an explicitly composed stage-level reference process.

The final design must state which reference is authoritative.

## Phase 4: Establish AMR Moment Fidelity

### Current risk

The code computes nonlinear `u` and `kappa` fields and then applies generic
cell-centered restriction and prolongation. In general:

```text
restrict(f(probabilities)) != f(restrict(probabilities))
```

where `f` contains the nonlinear `C_-^2` term.

### Recommended design study

Compare these approaches:

1. Transfer derived drift and covariance components.
2. Transfer raw first and second displacement moments, then construct central
   covariance on the destination level.
3. Transfer corrected mass fluxes and density, then reconstruct probabilities
   and moments on the destination level.

Approach 2 is the most direct way to preserve the intended statistical
quantities. Approach 3 is preferable if level-local consistency with the
finite-volume update can be demonstrated.

### Required AMR tests

Use static and adaptive refinement:

- uniform flow crossing a stationary coarse/fine interface;
- diagonal 2D flow crossing a refinement corner;
- nonuniform compressible flow;
- particles seeded on both sides of the interface;
- repeated refine/derefine cycles;
- MPI ownership changes during refinement;
- restart immediately before and after refinement.

Measure:

- per-level mean and covariance;
- discontinuity in coefficients across the interface;
- particle density relative to gas density;
- dependence on which side of the interface particles were seeded;
- convergence when both levels are globally refined.

The existing AMR smoke test should remain, but it is not sufficient evidence
of moment fidelity.

## Phase 5: Fix Boundary and Source-Term Semantics

### Exact periodic upper boundary

Particle domains are otherwise treated as half-open intervals
`[xmin, xmax)`, but periodic wrapping checks `x > xmax`.

Change periodic wrapping to include `x == xmax`. Add tests that place a
particle exactly at every active upper boundary, including corners, then run:

- GID reassignment;
- CIC interpolation;
- MPI migration;
- restart output; and
- another Ito push.

### Mass-changing source terms

Flux tracers follow mass exchanged through cell faces. A source that changes
cell mass without a corresponding tracer operation breaks tracer-to-gas mass
proportionality.

Choose one policy:

1. reject Ito and MC mass-flux tracers when a mass source is enabled;
2. document that mass-changing sources are unsupported and provide a runtime
   hook for modules to declare them; or
3. implement tracer birth, death, or weight changes consistent with the source.

Add a manufactured density-source test. It should either fail at startup with
a precise message or demonstrate the intended tracer source coupling.

## Phase 6: Strengthen Restart and Decomposition Reproducibility

Current restart tests prove only that the restart command succeeds.

Add paired runs:

```text
A: cycles 0 -> N without restart
B: cycles 0 -> K, restart, then K -> N
```

Compare, sorted by 64-bit tag:

- position bit patterns;
- GID;
- seed ID and creation time;
- schedule state and next tag;
- tracked thermodynamic variables;
- final history blocks.

Run the comparison for:

- serial;
- 2 and 4 MPI ranks;
- rank-count change across restart, if supported;
- uniform and AMR meshes;
- RK2 and RK3;
- particles crossing periodic boundaries.

If bitwise identity is not a supported contract, specify and test the exact
numerical tolerance and explain the source of nondeterminism.

## Phase 7: Add Physics Validation

### Square-pulse reference

Port the paper's square-pulse test in the closest AthenaK-compatible form.
Compare MC and Ito-2 using:

- mean;
- variance or sheet width;
- skewness and kurtosis;
- binned chi-square or another predeclared distribution test;
- convergence with particle number;
- convergence with spatial resolution.

The expected result should not claim that Ito-2 reproduces the full MC
distribution. The paper reports agreement in width but disagreement in
distribution shape. A correct regression should preserve that known
distinction.

### Three-dimensional turbulence

Add an expensive validation tier modeled on the paper's turbulent test:

- gas and tracer density correlation;
- PDF of tracer-to-gas density ratio;
- column-density comparison;
- tracer and gas power spectra;
- dependence on particles per cell;
- comparison among MC, Ito-2, and any future Ito-3 implementation.

This test belongs in a nightly, release, or publication-validation tier rather
than every pull request.

### Resolution and solver dependence

Run at multiple:

- resolutions;
- CFL numbers;
- reconstruction methods;
- Riemann solvers; and
- RK integrators.

The interpretation must state that the current method follows the numerical
mass flux and its scheme-dependent mixing. It is not a unique reconstruction
of physical fluid parcel trajectories.

### Ito-3 decision

If distributional agreement with MC is a requirement, implement and validate
Ito-3. If not, keep Ito-3 rejected and state explicitly that Ito-2 matches only
the selected lower-order statistics.

## Phase 8: Improve Thermodynamic Sampling

The transport coefficients use CIC interpolation, while history output samples
the single containing cell. This causes discontinuous thermodynamic histories
when a continuous trajectory crosses a cell boundary.

Add a selectable history sampling mode:

```ini
<outputN>
particle_field_sampling = cell
```

or:

```ini
particle_field_sampling = cic
```

Retain containing-cell sampling for discrete phase membership and exact cell
diagnostics. Use CIC for smooth path histories.

Tests must cover:

- uniform fields;
- linear gradients, where CIC has an exact expected value;
- discontinuities;
- periodic boundaries;
- coarse/fine interfaces;
- Hydro and MHD derived variables.

## Phase 9: Optimize After Correctness

Do not optimize the diagonal implementation and then redesign storage for the
full covariance tensor. Freeze the corrected coefficient contract first.

### Coefficient work

- Store and communicate only components required by the active dimension.
- Consider storing per-step displacement moments to avoid repeated `dt`
  scaling.
- Replace full-array zeroing when every active interior element is overwritten.
- Fuse coefficient validation into a reduction with one host synchronization.
- Avoid drawing inactive-coordinate random values.
- Measure the cost of covariance factorization separately.

### Particle migration

- Reuse send-list capacity instead of reallocating to `npart` each step.
- Avoid global count collectives when no consumer needs updated global counts.
- Audit all count and displacement products for 32-bit MPI overflow.

### Seeding

- Avoid copying the complete fluid state to the host for each seed event.
- Compute eligible-cell weights and prefix sums on device.
- Append into reserved particle capacity instead of copying all old particles.
- Batch seed events due on the same step.

### History output

- Evaluate fields on device.
- Avoid gathering all particle records to rank zero.
- Use sharded per-rank output or collective MPI-IO.
- Use 64-bit counts and offsets for large outputs.
- Document ordering guarantees.

### Benchmark matrix

Benchmark:

```text
particles: 10^4, 10^5, 10^6, 10^7 where feasible
dimensions: 2D, 3D
mesh: uniform, AMR
ranks: 1, 2, 8, 32+
backends: serial/OpenMP, CUDA or HIP
methods: MC, old Ito-2, corrected Ito-2
```

Record:

- wall time per particle update;
- coefficient-build time;
- communication time;
- migration time;
- output time;
- peak device and host memory;
- strong-scaling efficiency.

## Phase 10: Repair Unrelated Branch Gates

### Style failures

The current repository style gate fails in branch-added cooling files:

- C++ line-length failures in `src/srcterms/cooling.cpp`,
  `src/srcterms/cooling.hpp`, and `src/srcterms/srcterms.cpp`;
- Python indentation and formatting failures in
  `tools/validate_cooling_table.py`,
  `tst/test_suite/cooling/plot_cooling_tests.py`, and
  `tst/test_suite/cooling/__init__.py`.

Run:

```bash
cd tst
python run_test_suite.py --style
```

The Ito documentation must not claim that the style gate passes until this
command passes on the exact branch tip.

### Single-precision build

The branch currently fails an AppleClang single-precision build because of
narrowing conversions in `src/coordinates/cartesian_ks.hpp`.

Run:

```bash
cmake -S . -B /tmp/ito-single \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_SINGLE_PRECISION=ON
cmake --build /tmp/ito-single -j 8
```

This failure is not introduced by Ito, but single precision cannot be listed
as tested until the repository builds and the Ito probability/covariance tests
pass in that configuration.

### Latent multiple-pack particle count issue

`Mesh::UpdateParticleCounts()` loops over `nmb_packs_thisrank` but repeatedly
uses `pmb_pack` without advancing to another pack. The current tree appears to
use one pack per rank, so this is latent. Either:

- assert that only one pack is supported; or
- implement correct iteration before enabling multiple packs.

### Documentation publication

The Ito page is not currently present on `origin/gh-pages`.

Validate publication in a temporary detached worktree based on the live remote:

```bash
git worktree add --detach /tmp/athenak-ito-gh-pages origin/gh-pages
```

Integrate:

- `docs/source/modules/ito_tracers.md`;
- figures;
- `docs/source/index.md`;
- `docs/source/modules/index.md`;
- relevant links from particle and output pages.

Then run:

```bash
cd /tmp/athenak-ito-gh-pages/docs
make clean html SPHINXOPTS="-W --keep-going"
```

Verify the rendered page and links in a browser before publication.

### Branch integration scope

The reviewed branch is 28 commits and 151 files ahead of the reviewed
`origin/main` baseline. It includes cooling, initial perturbations, divB AMR,
MC tracer, and Ito tracer work.

Before merge, choose one strategy:

1. split Ito and its minimum MC-tracer dependency into a focused branch;
2. merge the prerequisite feature branches first and rebase Ito onto that
   result; or
3. retain the combined branch and run the complete test surface for every
   included subsystem.

Do not treat passing Ito tests as evidence that the other bundled features are
merge-ready.

## Documentation Corrections

Update the module page and technical note to state:

- current or corrected multidimensional covariance behavior;
- whether the claim is per-coordinate or full-tensor;
- the chosen RK stage semantics;
- AMR validation scope;
- source-term limitations;
- restart reproducibility guarantee;
- expected Ito-2 distribution-shape limitation;
- tested precision and backends;
- exact validation date and commit.

Remove historical claims such as "the style gate passed" when they are not
true for the current branch tip.

## Final Test Matrix

### Pull-request tier

- Release serial build.
- Release MPI build.
- Style gate.
- `git diff --check`.
- 1D marginal moments.
- 2D and 3D full covariance.
- tag overflow and format compatibility.
- exact periodic upper boundary.
- serial restart equivalence.
- two-rank MPI migration and restart.
- static-AMR interface moments.
- invalid probability and unsupported-input rejection.
- history reader compatibility.

### Nightly tier

- adaptive-AMR refine/derefine cycles.
- four-rank decomposition reproducibility.
- single precision.
- GPU build and Ito tests.
- high-particle-count migration.
- output scaling.
- square-pulse distribution comparison.
- multiple CFL, solver, and reconstruction combinations.

### Release/publication tier

- 3D turbulent validation.
- resolution and particle-number convergence.
- MC versus Ito-2 statistical comparison.
- Ito-3 comparison if implemented.
- performance and memory report.
- strict `gh-pages` build, link check, and browser review.

## Final Acceptance Checklist

- [ ] The mathematical contract is approved.
- [ ] Full covariance is implemented, or all claims say per-coordinate only.
- [ ] Covariance tests pass in 2D and 3D.
- [ ] 64-bit tags pass overflow, migration, restart, and output tests.
- [ ] RK stage semantics are documented and tested.
- [ ] AMR moment tests pass.
- [ ] Exact upper-boundary tests pass.
- [ ] Mass-changing source behavior is explicit.
- [ ] Restart equivalence passes.
- [ ] Square-pulse validation records the expected Ito-2 limitation.
- [ ] Single precision has a documented pass or explicit exclusion.
- [ ] CPU, MPI, and GPU support claims match actual tests.
- [ ] Style passes.
- [ ] Performance regressions are within an approved budget.
- [ ] Documentation is correct on the feature branch.
- [ ] The page builds and renders in the live `gh-pages` tree.
