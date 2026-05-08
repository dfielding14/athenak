# TRML Lagrangian MC Particle Tracking Test Plan

This plan defines the validation suite for TRML `lagrangian_mc` tracer particles,
with emphasis on thermodynamic histories, AMR bookkeeping, physical-boundary
particle exits, inflow injection, restart continuity, and MPI/rank crossing.

The goal is to test the particle machinery with small deterministic problems
before trusting it in a production TRML run. A full turbulent mixing layer is too
complicated for first-line verification because transport, AMR, cooling, and
stochastic sampling all vary at once.

## Scope

The tests should verify these behaviors:

1. Tracked particles write correct thermodynamic histories in `.trk` output.
2. The Monte Carlo pusher moves particles according to the saved density fluxes.
3. Particles that cross non-periodic outflow boundaries are removed cleanly.
4. New tracers can be injected through physical inflow boundaries without tag
   reuse or bad mesh ownership.
5. Scheduled tracking can reserve stable `.trk` rows for a deterministic initial
   sample plus later physical top-boundary injections.
6. AMR refinement and derefinement keep particle position, `gid`, and level
   bookkeeping consistent.
7. Restarted runs produce the same tracked-particle state as uninterrupted runs.
8. Serial and MPI runs agree for the same deterministic setup.

The current `.trk` record contains fixed columns:

```text
tag x y z vx vy vz rho press temp eint scalar0 gid level active
edot_cool edot_heat edot_net dedt_rad_mass dTdt_rad tcool
entropy ln_entropy divv dTdt_ad T_mix_scalar T_minus_T_mix_scalar
T_label_mix T_minus_T_label_mix gradT_mag grad_scalar_mag strain_mag
```

For legacy tracking, row number and `tag` both identify the low-tag particle.
For `track_selection = scheduled_injection`, row number is the stable tracking
slot and the `tag` column is the current physical particle identity. Unfilled
scheduled slots must have `active = 0` and `tag = -1`.

For TRML runs with `nscalars >= 2`, `scalar1` is initialized as the passive
mixing-only thermal label sampled in the track file as `T_label_mix`.  It is
advected but not radiatively cooled, so `T - T_label_mix` is the cleanest
particle-local diagnostic of cooling beyond passive mixing, up to compression
and shock terms.

## Testing Principles

Use a purpose-built deterministic test problem, not a full TRML simulation, for
the primary validation. The test should use analytic cell fields that make
indexing, level, and stale-cell mistakes obvious. A useful pattern is:

```text
rho     = 10 + i + 10*j + 100*k + 1000*level
temp    = 1 + 0.01*i + 0.1*j + k
scalar0 = x + 2*y + 3*z
```

The Python oracle should reconstruct the expected cell from each particle
position and AMR level, then compare the expected values against the `.trk`
columns. Exact tag-by-tag checks should be preferred wherever possible.
Statistical Monte Carlo checks should be secondary and used only where exact
realizations are intentionally not fixed.

The test suite should keep runtime short. Small meshes, short stop times, and
fixed seeds are preferred because failures must be fast to reproduce.

## Proposed Test Cases

### 1. Static Thermodynamic Sampling

Setup:

- Zero velocity and zero saved mass flux.
- Fixed particle positions.
- Analytic `rho`, `temp`, `press`, `eint`, velocity, and `scalar0` fields.
- Track all initialized particles.

Assertions:

- `active = 1` for every initialized tracked particle.
- Positions do not change between outputs.
- `rho`, `press`, `temp`, `eint`, `scalar0`, and velocity columns match the
  analytic oracle for each particle cell.
- `gid` and `level` match the owning MeshBlock.

This isolates track-output sampling from transport.

### 2. Deterministic Monte Carlo Push

Setup:

- Uniform density.
- Controlled one-direction density flux across selected faces.
- Fixed `random_seed`.
- A small number of particles placed in known cells.

Assertions:

- The Python oracle reproduces the particle random draw from `tag`, cycle, and
  seed, then predicts the exact tag-by-tag move.
- Particle positions move to the expected neighboring cell centers.
- Thermodynamic history columns after the move match the destination cell.
- A zero-flux control group remains fixed.

This tests the stochastic pusher without AMR or physical boundaries.

### 3. AMR Refine and Derefine Bookkeeping

Setup:

- Force a small region to refine at a known time.
- Later force the same region to derefine.
- Place tracked particles in cells that will become children and later parents.
- Keep analytic fields level-dependent so wrong level sampling is obvious.

Assertions:

- Particles remain active unless they physically leave the domain.
- Particle `gid` points to a valid MeshBlock after refinement and derefinement.
- The `level` column changes as expected when particles are remapped.
- Thermodynamic columns match the analytic oracle at the current level.
- Logs contain no orphan-particle or invalid-particle messages.

This is the critical AMR-specific validation.

### 4. Non-Periodic Boundary Exit

Setup:

- Use an outflow physical boundary in one direction.
- Place particles in the boundary-adjacent cells.
- Drive outward mass flux through the boundary.

Assertions:

- The oracle predicts which tracked tags leave the domain.
- Exited tracked tags appear with `active = 0` in `.trk` output after removal.
- Active particle count decreases by the expected amount.
- No exited particle retains a valid active `gid` or level.
- No out-of-bounds cell index is reported.

This verifies lifecycle handling for particles that leave the computational box.

### 5. Inflow Injection

Setup:

- Use an inward physical boundary flux.
- Set `inject_at_inflow = true`.
- Set `mass_per_particle` so the expected injected count is an integer for at
  least one output interval.
- Set track `nparticles` larger than the initial particle count so newly injected
  tags are observable.

Assertions:

- New tags are unique and contiguous from the previous `next_tag`.
- No initial tag is reused.
- Injected particles start in the correct boundary cells on the correct AMR
  level.
- Injected tracked particles have `active = 1` after injection.
- Their thermodynamic history columns match the injection cell.
- `max_inject_per_step`, when set, caps the number of injected particles without
  corrupting tag order.

This verifies both particle creation and track-output visibility for particles
that did not exist at initialization.

### 6. Scheduled Top-Injection Tracking

Setup:

- Use `track_selection = scheduled_injection`.
- Set `track_initial_fraction` to a small deterministic value or zero.
- Schedule injected slots before a controlled top inflow begins.
- Keep another physical inflow face active before the top inflow starts.

Assertions:

- Pre-scheduled slots remain `active = 0`, `tag = -1` before
  `track_inject_t_start`.
- Slots that become due while no eligible top particles exist remain inactive.
- Non-selected inflow faces do not fill scheduled injected slots.
- The first later physical particles injected through `track_inject_face` fill
  the due slots, with stable row indices and visible non-negative tags.
- Thermodynamic columns match the cells occupied by those tracked particles.

This verifies that scheduled tracking follows physical top-boundary tracers
without forcing diagnostic-only particles.

### 7. Restart Continuity

Setup:

- Run the same deterministic case in two ways:
  - one uninterrupted run from `t = 0` to `t = t_end`;
  - one split run from `t = 0` to restart time, then restart to `t_end`.
- Include at least one injected particle before restart and one after restart.

Assertions:

- The final tag set is identical between uninterrupted and restarted runs.
- Positions, active flags, `gid`, `level`, and thermodynamic columns agree for
  each tracked tag.
- Injected tags after restart continue from the correct `next_tag`.

This test should fail if `next_tag` is not persisted or reconstructed correctly.

### 8. MPI and Rank-Crossing Parity

Setup:

- Run the deterministic push case in serial and with two MPI ranks.
- Arrange at least one particle to cross a rank-owned MeshBlock boundary.
- Keep the same seeds and output cadence.

Assertions:

- The tracked-tag histories are identical between serial and two-rank runs.
- Particle communication does not change tag order, active flags, or sampled
  thermodynamic fields.
- No particle is duplicated or lost at rank boundaries.

This validates the communication path separately from AMR remapping.

## Highest-Value Combined Test

After the individual tests pass, run one combined stress test:

- AMR enabled.
- A refinement event before restart and a derefinement event after restart.
- Inflow injection enabled.
- At least one particle exits through an outflow boundary.
- At least one particle crosses an MPI rank boundary.
- `nparticles` in the track output larger than the initial particle count.

Acceptance for this combined test:

- No duplicate tags.
- No invalid active `gid` or level.
- Expected exit tags become inactive.
- Expected injected tags become active and visible in `.trk`.
- Continuous and restarted runs agree.
- Serial and MPI runs agree where the domain decomposition is not intended to
  change the deterministic random sequence.

This is the strongest single validation of the fragile bookkeeping paths.

## Acceptance Criteria

The particle tracking implementation should be considered ready for TRML use
when all of the following are true:

1. `.trk` files parse cleanly and advertise the expected fixed columns.
2. Active flags are exact for initialized, exited, and injected particles.
3. No duplicate tags appear at any output time.
4. No active particle has an invalid `gid` or negative level.
5. Thermodynamic history columns match the analytic oracle to single-precision
   output tolerance.
6. AMR refinement and derefinement update particle ownership and level correctly.
7. Restarted and uninterrupted runs agree tag-by-tag.
8. Serial and two-rank MPI runs agree for deterministic cases.
9. Logs contain no orphan-particle, invalid-particle, or out-of-bounds warnings.

## Implementation Sketch

Recommended file layout:

```text
TRML_testing/
  docs/
    particle_tracking_test_plan.md
  particle_tests/
    athinput.particle_static
    athinput.particle_mc_push
    athinput.particle_exit_x1
    athinput.particle_inflow_injection
    athinput.particle_amr_refine_derefine
    run_particle_tracking_tests.py
    read_trk.py
    oracle.py
src/pgen/
  particle_tracking_test.cpp
```

The initial executable harness is implemented under `TRML_testing/particle_tests/`.
It uses the custom pgen `src/pgen/particle_tracking_test.cpp`, which should be
built with:

```text
cmake -S . -B TRML_testing/particle_tests/build_particle_tracking -D PROBLEM=particle_tracking_test
cmake --build TRML_testing/particle_tests/build_particle_tracking -j 4
```

The harness can build this automatically and run the implemented checks with:

```text
python3 TRML_testing/particle_tests/run_particle_tracking_tests.py --jobs 4
```

Currently implemented cases:

1. `static`: analytic thermodynamic history sampling in 3D.
2. `mc_push`: exact tag-by-tag one-cycle Monte Carlo x1 transport.
3. `exit_x1`: removal of particles through a non-periodic outflow boundary.
4. `inflow_injection`: creation and tracking of new particles at an inflow face.
5. `amr_refine_derefine`: serial AMR refinement and derefinement bookkeeping.

Remaining planned extensions:

1. Restart continuity, including injected-particle `next_tag` continuity.
2. Serial versus two-rank MPI parity for rank-crossing trajectories.

The test runner should:

1. Build or locate the AthenaK executable.
2. Run each tiny input in a fresh scratch directory.
3. Parse `.trk` output records.
4. Compare each tracked tag against the Python oracle.
5. Print a compact pass/fail table.
6. Save failed-run logs and parsed comparison tables for debugging.

Generated test outputs should stay out of version control.

## Phased Rollout

Phase 1:

- Implement the `.trk` reader.
- Add the static thermodynamic sampling test.

Phase 2:

- Add deterministic Monte Carlo transport.
- Add non-periodic boundary exit.

Phase 3:

- Add inflow injection.
- Add AMR refinement and derefinement.

Phase 4:

- Add restart continuity.
- Add serial versus two-rank MPI parity.

The suite should become part of the normal local validation path before running
large TRML particle simulations.
