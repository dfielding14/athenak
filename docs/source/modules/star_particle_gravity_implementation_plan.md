# Star Particle Gravity Implementation Plan

This document lays out the design and implementation plan for a part 2 star-particle
feature branch that builds on `feature/star-particles`.  The goal is to add
gravitational interactions between star particles without losing the formation,
accretion, particle accounting, boundary migration, tests, and documentation already
present in the standalone star-particle branch.

This is intentionally an implementation plan, not user-facing feature documentation.
The user documentation should be written after the behavior and input contract are
implemented and tested.

Implementation status for `feature/star-particle-gravity`: the branch implements this
milestone with direct and replicated-tree star-star gravity, constant and point-mass
external accelerations, a pgen callback for custom external accelerations, KDK and RK4
integrators, gravity diagnostics, timestep limiting, momentum-aware accretion, and
CPU/MPI regression coverage.

## Starting Point

The existing star-particle branch provides:

- `particle_type = star` in the `<particles>` block.
- Star-particle file initialization.
- Dynamic star formation from dense gas cells.
- Optional removal of gas mass from the host cell when a star forms.
- Gas accretion from the host cell and neighboring cells.
- Particle mass stored in `IPMASS`.
- Particle creation time stored in `IPT_CREATE`.
- Particle migration and outflow removal through the particle boundary path.
- Particle count and mass history output.
- CPU and MPI regression tests for initialization, formation, accretion, and
  outflow removal.

The current particle pusher is not a gravitational dynamics integrator.  It only
drifts positions from existing velocities.  Adding star-star gravity therefore requires
a real force calculation, a real particle integrator, timestep control, diagnostics,
and new tests.

## Recommended Branch

Create a new branch from the current star-particle feature branch:

```bash
git switch feature/star-particles
git switch -c feature/star-particle-gravity
```

This should not start from `origin/main`, because the part 2 work depends directly on
the star-particle type, mass field, formation/accretion logic, and boundary fixes.

## Goals

The first complete gravity branch should provide:

1. Softened Newtonian gravity between star particles.
2. Explicit external accelerations for star particles.
3. A required user-specified gravitational constant.
4. At least two particle integrator choices.
5. A direct pairwise force backend for correctness and small systems.
6. A scalable large-N backend behind the same force interface.
7. MPI-correct particle force calculation and migration.
8. Gravity-aware particle timestep control.
9. Momentum-aware accretion when gas is added to a star particle.
10. Diagnostics for conservation and timestep behavior.
11. CPU and MPI tests that validate force accuracy, orbit behavior, conservation,
    formation plus gravity, and accretion plus gravity.
12. Developer documentation and user-facing documentation once the feature works.

## Non-Goals For The First Gravity Branch

The following are important, but should not be allowed to block the first correct
star-star gravity implementation:

- Full gas self-gravity.
- A general mesh Poisson solver.
- Star-gas gravitational backreaction.
- Automatically inferring particle gravity from hydro/MHD pgen source terms.
- A fully distributed tree or FMM implementation.
- Collisions, mergers, stellar feedback, sink-particle feedback, or stellar evolution.
- Relativistic gravity.

The implementation should be structured so these can be added later without rewriting
the star-star gravity path.

## Major Design Decisions

### 1. Physics Scope

Recommended first scope:

- Implement star-star gravity and explicitly configured external acceleration.
- Do not apply gravitational forces to gas in the first branch.
- Do not include gas mass in the gravity source.
- Document that formation/accretion couple gas mass to particle mass, but gas does not
  source gravity and does not feel star-particle gravity.

Later extensions can add:

- star-to-gas gravitational source terms;
- additional gas-to-star gravity models beyond the explicit external acceleration
  interface;
- full gas/star self-gravity with a Poisson solver;
- feedback or sink-particle models.

Rationale: the current code has constant-acceleration source terms, but no general
Poisson solver.  Star-star gravity can be implemented, tested, and used without
building a full gravity module for the gas.

### 2. Input Contract

Add a dedicated `<star_gravity>` block.  Gravity is a distinct physical module and
should not be hidden among formation and accretion parameters in `<particles>`.

Proposed input:

```ini
<particles>
particle_type = star
pusher = gravity

<star_gravity>
enabled = true
grav_constant = 1.0
force_method = direct
integrator = kdk
softening_type = plummer
softening_length = 0.01
timestep_mode = acceleration
timestep_eta = 0.02
max_timestep = -1.0
periodic_mode = domain
tree_theta = 0.7
```

Parameter meanings:

| Parameter | Recommended choices | Meaning |
| --- | --- | --- |
| `enabled` | `true`, `false` | Enables star-star gravity. |
| `grav_constant` | positive real | Code-unit gravitational constant. Required when enabled. |
| `force_method` | `direct`, `tree` | Force backend. Start with `direct`; add `tree` next. |
| `integrator` | `kdk`, `rk4` initially | Particle time integration method. |
| `softening_type` | `plummer` initially | Gravitational softening kernel. |
| `softening_length` | non-negative real | Softening length in code units. |
| `timestep_mode` | `fixed`, `acceleration`, `pair_orbit` | Gravity timestep limiter. |
| `timestep_eta` | positive real | Accuracy parameter for non-fixed timestep modes. |
| `max_timestep` | real, negative means disabled | Optional absolute cap on particle timestep. |
| `periodic_mode` | `domain`, `minimum_image`, `error` | How to treat domain periodicity for star-star gravity. |
| `tree_theta` | positive real | Barnes-Hut opening angle for `force_method = tree`. |

Important validation rules:

- `grav_constant` must be present and positive when gravity is enabled.
- `softening_length` must be non-negative.
- `tree_theta` must be positive when the tree backend is selected.
- Gravity should require `particle_type = star`.
- Gravity should require a 3D mesh, matching the current star-particle contract.
- Unknown `force_method`, `integrator`, `softening_type`, `timestep_mode`, or
  `periodic_mode` values should fail early with clear errors.

### 3. Force Law And Softening

Recommended first force law:

```text
a_i = sum_j G m_j r_ij / (|r_ij|^2 + eps^2)^(3/2)
```

where `r_ij = x_j - x_i` and `eps = softening_length`.

For diagnostics, the matching Plummer potential energy is:

```text
U = - sum_{i<j} G m_i m_j / sqrt(|r_ij|^2 + eps^2)
```

Design choices:

- `eps = 0` should be allowed only if explicitly requested and tested.
- The code should protect against exact zero separation producing non-finite
  accelerations.
- If two particles occupy the same position and `eps = 0`, fail or skip the singular
  pair with a diagnostic depending on the configured policy.  The safer first choice is
  to require positive softening for gravity tests and documented production examples.

### 4. Periodic Boundaries

There are several physically distinct choices:

| Mode | Meaning | Recommendation |
| --- | --- | --- |
| `domain` | Use positions as stored; boundary conditions only affect particle migration. | Default first implementation. |
| `minimum_image` | For periodic dimensions, use the nearest periodic image. | Add after direct backend is stable. |
| `error` | Abort if any mesh dimension is periodic. | Useful for conservative validation runs. |

Do not claim to implement true periodic self-gravity unless an Ewald-like or mesh-based
periodic gravity method is actually added.

### 5. Force Backend Architecture

Add a force-interface layer instead of embedding one force loop directly in the pusher.

Recommended files:

- `src/particles/star_gravity.hpp`
- `src/particles/star_gravity.cpp`
- `src/particles/star_gravity_direct.cpp`
- `src/particles/star_gravity_tree.cpp` once the tree backend is implemented
- `src/pgen/tests/star_gravity_test.cpp`

Recommended data structures:

```cpp
enum class StarGravityForceMethod {direct, tree};
enum class StarGravityIntegrator {kdk, rk4};
enum class StarGravitySoftening {plummer};
enum class StarGravityTimestepMode {fixed, acceleration, pair_orbit};
enum class StarGravityPeriodicMode {domain, minimum_image, error};

struct StarGravityConfig {
  bool enabled;
  Real grav_constant;
  StarGravityForceMethod force_method;
  StarGravityIntegrator integrator;
  StarGravitySoftening softening;
  StarGravityTimestepMode timestep_mode;
  StarGravityPeriodicMode periodic_mode;
  Real softening_length;
  Real timestep_eta;
  Real max_timestep;
  Real tree_theta;
};

struct StarGravityDiagnostics {
  Real kinetic_energy;
  Real potential_energy;
  Real total_mass;
  Real max_acceleration;
  Real min_pair_separation;
  Real momentum_x;
  Real momentum_y;
  Real momentum_z;
  Real angular_momentum_x;
  Real angular_momentum_y;
  Real angular_momentum_z;
  Real center_of_mass_x;
  Real center_of_mass_y;
  Real center_of_mass_z;
  Real gravity_dt;
};
```

The `Particles` class should own gravity configuration and scratch arrays, either
directly or through a small `StarGravity` helper.  A helper class is cleaner because
gravity will otherwise make `Particles` too large.

### 6. Direct Force Backend

The direct backend should be implemented first and treated as the correctness oracle.

Recommended MPI strategy:

1. Every rank packs its local star particles into compact host arrays:
   `x`, `y`, `z`, `vx`, `vy`, `vz`, `mass`, `tag`.
2. Every rank participates in collective count exchange, even if it owns zero particles.
3. Use `MPI_Allgatherv` to build a global particle snapshot on each rank.
4. Copy the global snapshot to device memory if the force loop runs on device.
5. Compute accelerations only for locally owned particles.
6. Reduce diagnostics over all ranks.

This gives:

- exact pairwise forces for the chosen softened law;
- simple MPI semantics;
- deterministic tests against serial runs;
- a baseline for tree-force error tests.

Cost:

- memory is `O(N_total)` per rank;
- work is `O(N_local * N_total)`.

That is acceptable for the correctness backend, not for very large production runs.

### 7. Large-N Backend

The first scalable backend should be a Barnes-Hut tree with a replicated global tree:

1. Gather the global particle snapshot exactly as in the direct backend.
2. Build a tree on each rank from the global snapshot.
3. Compute approximate accelerations for local particles using opening angle
   `tree_theta`.
4. Compare tree accelerations against direct accelerations in dedicated tests.

This is not fully distributed, but it changes force evaluation from `O(N^2)` to roughly
`O(N log N)` while preserving a simple MPI model.  It is a realistic intermediate step
before a distributed tree or particle-mesh method.

Tree implementation choices:

- Start with monopole moments only.
- Use an axis-aligned octree over the simulation domain or over the particle bounding
  box.
- Store node total mass, center of mass, bounding size, child indices, and particle
  index range.
- Avoid recursion on device unless the target backend supports it cleanly.
- Consider building the tree on host first, then optimizing later.

Validation requirements:

- tree force agrees with direct force to a documented tolerance for a random particle
  distribution;
- error decreases as `tree_theta` decreases;
- total momentum drift remains controlled for symmetric systems;
- tree backend works with zero-particle ranks.

### 8. Particle Integrators

Implement at least two integrators.

#### KDK Leapfrog

This should be the default production integrator.

Recommended task sequence:

1. Form stars.
2. Accrete gas onto stars.
3. Compute acceleration at time `n`.
4. Kick velocities by `0.5*dt`.
5. Drift positions by `dt`.
6. Run particle boundary migration and outflow removal.
7. Recompute acceleration at time `n+1` after migration.
8. Kick velocities by `0.5*dt`.
9. Update diagnostics and gravity timestep.

This requires splitting the current single `Push` task into gravity-aware subtasks, or
adding gravity subtasks around the existing migration tasks.

Reasons to prefer KDK:

- time-centered velocities;
- good long-term energy behavior for fixed masses;
- simple diagnostics;
- standard for collisionless particle dynamics.

Important caveat:

- Star formation and accretion change particle number and particle mass, so the system
  is not exactly Hamiltonian.  Conservation tests should use runs with no formation and
  no accretion when evaluating integrator quality.

#### RK4

Add RK4 as a high-accuracy comparison integrator.

Use cases:

- short-term convergence tests;
- comparing against analytic two-body motion;
- catching leapfrog staging bugs.

Limitations:

- not symplectic;
- more force evaluations per step;
- less attractive for long integrations.

Implementation note:

- RK4 can update positions and velocities in one task and then use the existing
  migration path after the final update.
- Intermediate states should use scratch arrays, not overwrite live particle arrays.

#### Optional Later Integrators

Consider later:

- DKD leapfrog;
- velocity Verlet;
- adaptive subcycling;
- Hermite integrator with jerk estimates for close N-body systems.

Do not add these before the KDK and RK4 paths are correct and tested.

### 9. Timestep Control

Gravity must update `Particles::dtnew`, because `Mesh::NewTimeStep()` already includes
the particle timestep in the global timestep minimum.

Recommended modes:

| Mode | Formula | Use |
| --- | --- | --- |
| `fixed` | use `<particles>/dt` or `<star_gravity>/max_timestep` | Simple controlled tests. |
| `acceleration` | `dt = eta * sqrt(eps / max(|a|))` | Softened force default. |
| `pair_orbit` | `dt = eta * min(sqrt(r_ij^3 / (G*(m_i+m_j))))` | Better for close pairs. |

Rules:

- If `softening_length = 0`, acceleration mode needs a different length scale.  Use
  the minimum pair separation with a floor, or require `pair_orbit`.
- `max_timestep > 0` should cap all gravity timestep estimates.
- All ranks must participate in timestep reductions.
- The history output should record the final gravity timestep estimate.

### 10. Formation And Accretion With Gravity

Formation behavior should stay physically conservative:

- when gas density exceeds the threshold, create a star particle;
- if `star_remove_gas_on_formation = true`, remove gas mass from the host cell;
- the actual particle mass is the mass removed from gas;
- do not create a star if no mass can be removed above the configured floor.

Accretion should be upgraded for gravity:

Current behavior increases particle mass by gas accretion.  With gravity, accretion
should also update particle momentum:

```text
v_new = (m_old*v_old + dm*v_gas_removed) / (m_old + dm)
```

where `dm` and `v_gas_removed` are accumulated over the host and neighboring cells.

Recommended behavior:

- Add `star_accretion_conserve_momentum = true` with default `true`.
- Preserve the existing mass-only behavior only through an explicit opt-out.
- Reduce gas conserved mass, momentum, energy, and scalars consistently.
- Add tests that verify total gas-plus-star momentum is conserved during accretion when
  no other forces are active.

### 11. Task Ordering

The existing star-particle task order is:

1. Form stars.
2. Accrete stars.
3. Push particles.
4. Migrate particles.

Gravity-aware ordering should become:

1. Form stars.
2. Accrete stars.
3. Prepare/gather gravity state.
4. First gravity kick or RK stage.
5. Drift/update positions.
6. Migrate particles.
7. Rebuild/gather gravity state if the integrator needs post-migration forces.
8. Final gravity kick.
9. Update diagnostics and timestep.

The task graph must handle:

- no particles on one or more MPI ranks;
- particles formed on only one rank;
- particles removed by outflow on only one rank;
- particles crossing rank boundaries during a gravity step;
- formation and accretion before the gravity force is computed;
- all MPI collectives executed by all ranks.

### 12. Diagnostics And History Output

Extend particle history output with gravity-specific columns when star gravity is
enabled.

Recommended columns:

| Column | Meaning |
| --- | --- |
| `grav-dt` | Gravity timestep estimate. |
| `star-ekin` | Total star-particle kinetic energy. |
| `star-epot` | Total softened star-star potential energy. |
| `star-etot` | `star-ekin + star-epot`. |
| `star-px`, `star-py`, `star-pz` | Total star-particle momentum. |
| `star-lx`, `star-ly`, `star-lz` | Total star-particle angular momentum. |
| `star-com-x`, `star-com-y`, `star-com-z` | Center of mass. |
| `star-amax` | Maximum particle acceleration magnitude. |
| `star-rmin` | Minimum pair separation. |

For tree gravity, diagnostics should either:

- compute potential energy directly for test-sized runs; or
- clearly label the diagnostic as approximate if it uses tree approximations.

For regression tests, the direct backend should provide the reference diagnostics.

### 13. Restart And Output

Before implementation, inspect restart support for particles.  The gravity branch should
not silently lose star-particle state on restart.

Implementation status: mesh restart support was inspected in `src/outputs/restart.cpp`
and the restart constructor in `src/pgen/pgen.cpp`.  Mesh restart files still do not
serialize particle arrays, so this branch adds a first-class particle sidecar restart
output:

- `file_type = rst_prtcl` writes `rst_prtcl/<basename>.<file_number>.rst_prtcl`;
- `star_init = restart` loads the sidecar named by `particle_restart_file`;
- the binary header records magic/version, `Real` and `int` sizes, mesh time/cycle,
  total particles, particle array dimensions, maximum tag, formed/accreted totals, and a
  basename token;
- load validates sidecar time/cycle against the mesh restart before restoring particles;
- particle real/int arrays are restored exactly;
- `PGID` is recomputed from current mesh geometry, so MPI rank count can change;
- restored tags are preserved and `RefreshNextStarTag()` is used instead of
  reassigning tags;
- formed/accreted totals are restored once globally so history continuity survives MPI
  reductions.

The restart regression tests now cover full-versus-restarted agreement, mismatch
rejection, dynamic formation counter/tag continuation, and MPI repartitioning.

Particle outputs should include enough information for tests and analysis:

- positions;
- velocities;
- mass;
- tag;
- creation time;
- optionally acceleration if useful for debugging.

### 14. Tests

The gravity branch needs tests beyond simple run completion.

#### Direct Force Tests

Add a two-particle force test:

- two particles at known separation;
- known masses;
- known `G`;
- known softening length;
- expected acceleration from the Plummer formula;
- verify equal and opposite momentum tendency.

#### Two-Body Orbit Tests

Add an equal-mass circular binary:

- particles start at `x = +/- d/2`;
- velocities produce a circular orbit;
- no formation;
- no accretion;
- direct force backend;
- KDK integrator.

Validate:

- center of mass remains fixed;
- total momentum remains near zero;
- angular momentum is conserved to the expected tolerance;
- total energy drift scales down with timestep;
- final separation remains near the initial separation for a short run.

Add an eccentric orbit test after circular orbit is stable.

#### Integrator Comparison Tests

Run the same two-body setup with:

- `integrator = kdk`;
- `integrator = rk4`.

Validate:

- both remain bounded;
- RK4 and KDK agree over short integrations;
- KDK has acceptable long-term energy behavior for fixed-mass runs.

#### MPI Tests

Add MPI tests for:

- one particle on each rank;
- both particles on one rank and zero particles on another rank;
- particles crossing rank boundaries during an orbit;
- serial versus two-rank agreement for direct gravity;
- outflow removal during a gravity run.

All MPI tests should be written to catch collective mismatches.  A rank with zero local
particles must still participate in all gravity collectives.

#### Formation Plus Gravity Tests

Add a setup with:

- one file-loaded star particle;
- one dense gas cell that forms a second star;
- gravity enabled immediately after formation.

Validate:

- the dense cell loses mass if `star_remove_gas_on_formation = true`;
- global particle count updates correctly;
- the new particle participates in gravity on the same step or documented next step;
- total star mass equals initial star mass plus removed gas mass.

#### Accretion Plus Gravity Tests

Add a setup with:

- one or two star particles;
- gas accretion enabled;
- nonzero gas velocity;
- gravity optionally disabled for the conservation subtest.

Validate:

- particle mass increases by the removed gas mass;
- particle velocity changes according to momentum-conserving accretion;
- total gas-plus-star momentum is conserved when gravity is disabled;
- gravity acceleration changes when accretion changes particle mass.

#### Tree Backend Tests

After the tree backend exists:

- random particle distribution direct-versus-tree force comparison;
- decreasing error as `tree_theta` decreases;
- two-body tree path agrees with direct path when the opening criterion forces exact
  or near-exact behavior;
- MPI tree path handles zero-particle ranks.

### 15. Documentation

Update or add:

- `docs/source/modules/star_particles.md` with a short pointer to the gravity feature
  once it exists.
- `docs/source/modules/star_particle_gravity.md` as the user-facing gravity page.
- `src/particles/README.md` with developer notes on gravity architecture.
- Input examples under `inputs/particles/`.
- Test inputs under `tst/inputs/particles/`.

The user documentation should include:

- the physical model;
- all input parameters;
- recommended defaults;
- examples for direct and tree gravity;
- timestep guidance;
- softening guidance;
- MPI and scaling notes;
- known limitations;
- how formation and accretion interact with gravity.

### 16. Implementation Phases

#### Phase 0: Branch And Baseline

- Create `feature/star-particle-gravity` from `feature/star-particles`.
- Confirm the existing star-particle CPU and MPI tests still pass.
- Record baseline behavior and current limitations.

#### Phase 1: Input And Architecture

- Add `<star_gravity>` to allowed input blocks.
- Add gravity enums and config parsing.
- Add a `StarGravity` helper or equivalent.
- Add scratch arrays for acceleration and global snapshots.
- Add clear fatal errors for invalid parameter combinations.
- Add no-op gravity path when disabled.

#### Phase 2: Direct Backend

- Implement global particle snapshot gathering.
- Implement direct softened force computation.
- Implement acceleration and diagnostic reductions.
- Add force-law unit/regression test.

#### Phase 3: KDK Integrator

- Split gravity integration into task-list-safe subtasks.
- Add KDK kick/drift/kick behavior.
- Ensure particle migration happens between drift and final kick when needed.
- Update `dtnew` from gravity timestep estimates.
- Add circular binary tests.

#### Phase 4: Diagnostics

- Extend particle history output.
- Add energy, momentum, angular momentum, center-of-mass, `rmin`, `amax`, and
  gravity timestep columns.
- Add tests that parse these diagnostics.

#### Phase 5: Momentum-Conserving Accretion

- Modify accretion to accumulate removed gas momentum.
- Update particle velocity using mass-weighted momentum conservation.
- Add tests for gas-plus-star momentum conservation.

#### Phase 6: RK4 Integrator

- Add RK4 scratch state.
- Add short-run comparison tests against KDK and analytic expectations.
- Document that RK4 is not the recommended long-term collisionless integrator.

#### Phase 7: MPI Hardening

- Add zero-particle-rank tests.
- Add cross-rank orbit or drift tests.
- Add outflow-removal-with-gravity tests.
- Run serial and MPI tests repeatedly enough to catch collective-ordering issues.

#### Phase 8: Tree Backend

- Implement a replicated Barnes-Hut tree.
- Add direct-versus-tree force tests.
- Add tree orbit smoke tests.
- Document performance envelope and limitations.

#### Phase 9: Documentation And Final Validation

- Write user documentation.
- Update developer README.
- Add example inputs.
- Copy `docs/source/modules/star_particles.md`,
  `docs/source/modules/star_particle_gravity.md`, and
  `docs/source/_static/star_particles/` into a temporary `origin/gh-pages` worktree.
- Register both star-particle pages in the live `docs/source/modules/index.md` module
  table and hidden toctree during the `gh-pages` integration.
- Keep the live `docs/source/modules/particles.md` cosmic-ray drift page as the base
  particle-module entry point.  The star-particle pages document the standalone feature
  without hiding the existing public workflow.
- Validate the integrated Pages tree with
  `cd docs && make clean html SPHINXOPTS="-W --keep-going"`.
- Run style checks.
- Run CPU and MPI tests.
- Build docs if the docs infrastructure supports it.
- Commit and push the branch only after validation passes.

#### Phase 10: Review Hardening

- Preserve every co-located particle in Barnes-Hut leaf buckets.
- Pass RK4 stage positions and times to user external-acceleration callbacks.
- Recompute adaptive gravity timestep estimates from the configured particle timestep
  ceiling so the timestep can recover after forces weaken.
- Evaluate accretion stencils from a replicated star snapshot and reduce removed gas
  mass and momentum back to owners so the stencil crosses MeshBlock and MPI boundaries.
- Write `rst_prtcl` sidecars before mesh restart files regardless of input-block order.
- Reject malformed sidecars whose basename metadata, bookkeeping values, or complete
  byte count do not match the binary header.
- Allow large tree runs to set `exact_diagnostics = false` when exact pairwise history
  diagnostics would otherwise dominate runtime.
- Treat upper root-mesh bounds as half-open during particle migration so stars landing
  exactly on `xmax` are removed or wrapped correctly.
- Reject particle updates that cross more than one MeshBlock in any direction instead
  of indexing beyond the immediate-neighbor table.

### 17. Acceptance Criteria

The branch is not complete until:

- Existing star-particle tests still pass.
- New gravity tests pass in CPU mode.
- New gravity tests pass in MPI CPU mode.
- Direct two-particle forces match the analytic softened force.
- KDK binary orbit conserves center of mass, momentum, angular momentum, and energy to
  documented tolerances.
- RK4 agrees with KDK over a short integration.
- Accretion conserves gas-plus-star momentum when gravity is disabled.
- Formation-created stars participate in gravity correctly.
- Zero-particle MPI ranks do not deadlock.
- Outflow particle removal works during a gravity run.
- Exact upper-boundary outflow removal and oversized migration rejection are covered.
- Gravity timestep control limits `Mesh::dt`.
- Docs explain inputs, examples, limitations, and test coverage.
- `git diff --check` passes.
- The full relevant test suite has been run and recorded in the final branch summary.

## Open Questions To Decide Before Coding

1. Should `particle_type = star` use `pusher = gravity`, or should gravity be enabled
   solely by `<star_gravity>/enabled = true`?
2. Should `softening_length = 0` be allowed, or should gravity require positive
   softening?
3. Should the first implementation abort on periodic mesh boundaries unless
   `periodic_mode` is explicitly specified?
4. Should the replicated-tree backend be considered sufficient for the first large-N
   path, or should a distributed tree/particle-mesh method be required before merge?
5. Should particle restart support be a hard requirement for this branch?
6. Should accretion momentum conservation be introduced in this branch even though it
   changes behavior relative to the existing mass-only accretion path?
7. Should gravitational potential energy diagnostics use direct summation even when the
   tree backend is selected, at least for small test runs?

## Recommended First Milestone

The first milestone should be a scientifically honest star-gravity implementation:

- explicit `G`;
- Plummer softening;
- direct allgather force backend;
- replicated Barnes-Hut tree backend;
- explicit external acceleration interface;
- KDK integrator;
- RK4 cross-check integrator;
- gravity timestep limiter;
- diagnostics;
- formation plus gravity;
- momentum-conserving accretion;
- CPU and MPI tests.

The direct backend should remain permanently as the correctness oracle for tests and
debugging.
