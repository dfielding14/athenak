# Temporary MHD-PIC Implementation Review and Remediation Tracker

> **Temporary working document.** Update this file as fixes, tests, and
> qualification runs land. Delete it only after every required item in the
> final deletion checklist is complete.

## Baseline and scope

- Review branch: `PIC_development`
- Reviewed commit: `f8a56172983a1f009662af4217dcb7b52032d7dc`
- Review date: 2026-07-08
- Review type: source/static review, implementation remediation, and focused runtime tests
- Fresh AthenaK build or simulation performed during review: serial and MPI-enabled
  Debug builds plus bounded serial/MPI-singleton smoke runs
- Primary implementation reviewed:
  - `src/particles/`
  - `src/mhd/mhd_tasks.cpp` and `src/mhd/mhd_update.cpp`
  - `src/bvals/bvals_part.cpp` and `src/bvals/bvals_mom.cpp`
  - `src/pgen/tests/pic_parallel_shock.cpp`
  - Bell problem generators and campaign contracts
  - `src/srcterms/turb_driver.cpp`
  - Current turbulent-box, Bell, and parallel-shock inputs and tests

This review concentrates on coding practice, numerical stability, efficiency,
particle-fluid conservation, physical correctness, and the robustness of the
current campaign choices.

## Status convention

- **OPEN**: confirmed issue; no accepted fix yet.
- **IN PROGRESS**: implementation or validation is actively underway.
- **FIXED, NEEDS VALIDATION**: source fix exists but the required evidence is incomplete.
- **CLOSED**: fix, regression coverage, and relevant physical validation are complete.
- **ACCEPTED LIMITATION**: deliberately retained, clearly scoped, and protected from
  unsupported scientific claims.

Severity convention:

- **P0**: stop physical interpretation or production use of the affected workflow.
- **P1**: fix or explicitly qualify before science use.
- **P2**: engineering, performance, maintainability, or secondary robustness issue.

## Executive assessment

The uniform-grid, full-f `paper_mhd_pic_vl2_tsc` core is conceptually sound.
The explicit-midpoint chronology, relativistic Boris kick, TSC gather/deposit,
and opposite gas/particle momentum and kinetic-energy exchange are mutually
consistent on uniform grids.

The implementation is not yet generally science-ready. The immediate blockers
are the turbulence driver, mixed particle migration/destruction compaction,
AMR conservation, and feedback/timestep admissibility. The shock and Bell
campaigns also have important physical-qualification gates that must remain
separate from code-mechanics tests.

Until the corresponding items below are closed:

- Treat current turbulent-box outputs as exploratory.
- Do not trust multi-rank shock survivor populations when migration and physical
  escape can occur on the same rank and stage.
- Do not claim exact particle-fluid conservation on AMR/SMR meshes.
- Treat current nonlinear Bell runs as engineering candidates rather than
  qualifying physical evidence.

## Recommended remediation order

### 0. Protect interpretation of existing results

- [ ] Mark current turbulent-box results as exploratory pending `TURB-1` through
  `TURB-4`.
- [ ] Record which shock outputs could have encountered simultaneous migration and
  physical escape pending `MIG-1`.
- [ ] Keep AMR exact-conservation and nonlinear Bell claims fail-closed.

### 1. Fix active correctness blockers

- [x] Fix `TURB-1`: constant-energy normalization root.
- [x] Fix `TURB-2`: midpoint/VL2 forcing energy update.
- [x] Fix `TURB-3`: inactive-dimension forcing modes and solenoidal projection.
- [x] Fix `TURB-4`: signed isotropic Fourier-mode enumeration.
- [x] Fix `MIG-1`: unified migration/destruction compaction.
- [x] Add focused regressions for all five fixes before rerunning campaigns.

### 2. Add stability and admissibility protection

- [ ] Implement `STAB-1`: feedback source positivity/source-timescale control.
- [ ] Implement `DT-1`: globally safe gyrofrequency bound.
- [ ] Implement `DT-2`: post-injection timestep validation.
- [ ] Add fail-fast accounting for source-induced gas floors.

### 3. Resolve the AMR conservation policy

- [ ] Choose between a conservative cross-level deposition implementation and an
  explicitly accepted non-conservative `paper_smooth` limitation.
- [ ] Implement and validate the choice in `AMR-1`.
- [ ] Run matched uniform/AMR interface tests before using AMR shock outputs for
  conservation-sensitive conclusions.

### 4. Improve shock physical fidelity

- [ ] Freeze one startup injection history under `SHOCK-1`.
- [ ] Implement or validate measured shock-surface and mass-flux tracking under
  `SHOCK-2`.
- [ ] Add and monitor a CR-Hall applicability diagnostic under `SHOCK-3`.
- [ ] Re-run controlled planar reproduction tests before nonlinear shock/Bell claims.

### 5. Qualify the Bell workflow

- [ ] Prevent accidental physical use of the legacy normalization under `BELL-1`.
- [ ] Close the Q043 raw deposited-current matrix under `BELL-2`.
- [ ] Close corrected linear Bell growth, wavelength, phase, and polarization tests.
- [ ] Close PPC, resolution, timestep, decomposition, rigidity, and artificial-light-
  speed convergence gates before nonlinear use.

### 6. Scale and harden the implementation

- [ ] Address the highest-impact performance items in `PERF-1` through `PERF-4`.
- [ ] Address validation and MPI failure behavior in `ROBUST-1` through `ROBUST-4`.
- [ ] Make all source-local and default regression entry points green.

## Detailed findings

### TURB-1 — Incorrect constant-energy normalization root

- Status: **FIXED, NEEDS VALIDATION**
- Severity: **P0**
- Affected workflows: all current turbulent boxes with `constant_edot=true`
- Evidence: `src/srcterms/turb_driver.cpp:908-966`

The scale factor solves

\[
m_0 s^2 + m_1 s = \dot E,
\]

whose positive root always begins with `-m1/(2*m0)`. The `m1 < 0` branch instead
uses `+m1/(2*m0)`, so it does not solve the requested energy-injection equation
when fluid velocity and force are anti-correlated.

Required work:

- [x] Replace the branch with a numerically stable form of the correct positive root.
- [x] Add manufactured tests with positive, zero, and negative `m1`.
- [x] Verify measured `Delta E / Delta t` against configured `dedt` in serial and MPI.
- [x] Verify restart continuation preserves the same forcing normalization.

Current evidence: serial and true two-rank MHD/MHD-PIC runtime budgets pass. A
cycle-one restart continuation reproduces the uninterrupted final force exactly and
the final history to `1e-14`. An existing 2-D campaign checkpoint carrying the old
nonnegative default bounds also loads successfully through the compatibility upgrade;
its trajectory intentionally changes after the next corrected signed-mode refresh.
Runtime forcing with a deliberately negative `m1` remains the closure gap.

Closure criterion: the measured injection agrees with the requested value within
a documented discretization/roundoff tolerance for both signs of `m1`.

### TURB-2 — Midpoint/VL2 forcing double-counts quadratic kinetic energy

- Status: **FIXED, NEEDS VALIDATION**
- Severity: **P0**
- Affected workflows: all turbulent boxes using `paper_mhd_pic_vl2_tsc`
- Evidence:
  - midpoint weights: `src/driver/driver.cpp:159-167`
  - force update: `src/srcterms/turb_driver.cpp:1067-1092`
  - stage weight application: `src/srcterms/turb_driver.cpp:1348-1365`

Stage 2 evaluates work using the midpoint velocity and then adds another finite-kick
term, `0.5*a^2*dt^2`. For constant acceleration this gives

\[
\Delta E = \rho\left(v_0\cdot a\,\Delta t + a^2\Delta t^2\right)
\]

instead of

\[
\Delta E = \rho\left(v_0\cdot a\,\Delta t + \tfrac12 a^2\Delta t^2\right).
\]

When `v dot a` is small, the injected energy can approach twice the requested value.

Required work:

- [x] Make staged forcing obey the selected RK tableau.
- [x] Use midpoint work without the extra stage-2 quadratic term, or implement a
  consistent once-per-cycle operator-split exact kick.
- [ ] Add a production-path uniform constant-acceleration one-cycle invariant test.
- [x] Test MHD-only and MHD-PIC modes to ensure selecting PIC does not change the
  intended forcing power.

The implementation now applies only `rho*v.a` in every multi-stage RK source update;
the finite-kick quadratic term is retained only for RK1 and explicit standalone
impulses. A host tableau oracle covers uniform acceleration, while serial and true
two-rank stochastic Heun MHD and explicit-midpoint MHD-PIC runtime budgets pass. A
production-path manufactured constant-acceleration test remains required for closure.

Closure criterion: momentum and energy match the analytic constant-force solution
through the expected order, and measured forcing power is independent of the PIC
integrator selection.

### TURB-3 — Current 2-D forcing is not solenoidal

- Status: **FIXED, NEEDS VALIDATION**
- Severity: **P0** for 2-D turbulent boxes
- Evidence:
  - nonzero default `kz`: `src/srcterms/turb_driver.cpp:92-97`
  - mode generation: `src/srcterms/turb_driver.cpp:570-595`
  - flattened z dependence: `src/srcterms/turb_driver.cpp:429-448`
  - 3-D projection: `src/srcterms/turb_driver.cpp:647-667`

The 2-D basis removes spatial z dependence but the amplitude projection still uses
the nonzero three-dimensional wavevector. Therefore a nominally solenoidal force
generally has nonzero two-dimensional divergence.

Required work:

- [x] Force inactive-dimensional mode numbers to zero.
- [x] Project using only active spatial wavevector components.
- [x] Decide and document whether an out-of-plane acceleration component is retained
  in 2D3V.
- [x] Add an emitted-force spectral-divergence regression.

The supported isotropic driver retains the out-of-plane acceleration as a valid
2D3V solenoidal component. The inconsistent legacy `driving_type=1` path is now
rejected rather than exposed as a physical option. Emitted 2D3V and 3-D fields pass
spectral-divergence checks; a dedicated 1-D emitted-field case remains for closure.

Closure criterion: `sol_fraction=1` produces divergence consistent with roundoff and
the selected discrete derivative in every supported dimensionality.

### TURB-4 — Fourier mode enumeration is angularly biased

- Status: **FIXED, NEEDS VALIDATION**
- Severity: **P0** for isotropic-turbulence claims
- Evidence: `src/srcterms/turb_driver.cpp:210-231,570-588`
- Reference implementation: `src/srcterms/initial_perturbations.cpp:272-279,414-444`

The driver enumerates only nonnegative Cartesian mode components. The real-field
conjugate supplies `-k`, but mixed-sign directions such as `(1,-1,0)` remain absent.
This is not an isotropic shell.

Required work:

- [x] Adopt signed canonical-half-space enumeration.
- [x] Preserve deterministic seeding and restart behavior.
- [x] Test mode counts and absence of duplicate conjugate pairs.
- [x] Test the unweighted mode angular tensor for isotropy.
- [ ] Test stochastic forcing covariance over enough realizations.

Default isotropic bounds are the complete signed cube and asymmetric explicit bounds
are rejected. Runtime FFT checks confirm mixed-sign modes in 2-D and 3-D, and a
cycle-one restart is bitwise identical in the emitted force. Physical-k isotropy is
currently claimed only for equal active-axis box/tile lengths; stochastic covariance
and broader decomposition evidence remain required for closure.

Closure criterion: shell angular moments are isotropic within a documented finite-mode
tolerance and reproducible across decompositions and restarts.

### MIG-1 — Mixed migration and destruction corrupt particle identity

- Status: **FIXED, NEEDS VALIDATION**
- Severity: **P0** for multi-rank nonperiodic particle runs
- Evidence: `src/bvals/bvals_part.cpp:976-1035`

Send holes and destruction holes are compacted in separate passes. The send-hole pass
can copy a particle that is itself in `destroylist`. For example, with ten particles,
`send={1}`, `destroy={9}`, and no receives, particle 9 is copied into slot 1 and survives,
while shrinking the arrays drops valid particle 8. The escape ledger records the original
destruction before this corruption and can therefore appear correct.

Required work:

- [x] Merge send and destruction indices into one unique sorted hole set.
- [x] Fill receives and compact survivors against that unified set exactly once.
- [x] Move compaction to a device kernel or otherwise avoid per-particle deep copies.
- [x] Add a two-rank test with simultaneous inter-rank migration and physical escape.
- [x] Compare complete tag/source/species/state inventories, not only particle counts.

The exhaustive planner oracle and serial/two-rank end-to-end identity tests pass. In
the MPI fixture both ranks send and physically destroy particles in the same migration
call; all 180 survivor tags, metadata fields, positions, and velocities match serial
exactly. Restart, repeated crossing, and ledger parity remain closure requirements.

Closure criterion: exact survivor identity and boundary ledgers agree across serial,
MPI, restart, and repeated boundary-crossing cases.

### AMR-1 — `paper_smooth` is not conservative across refinement interfaces

- Status: **OPEN POLICY DECISION**
- Severity: **P1**
- Evidence:
  - implementation: `src/particles/particles_moments.cpp:629-933`
  - frozen raw totals: `tst/scripts/particles/pic_paper_smooth_tsc_oracle.py:117-126`
  - explicit no-renormalization policy:
    `tst/scripts/particles/pic_paper_smooth_tsc_oracle.py:251-263`

Receiver-resolution TSC is evaluated independently on each level without cross-interface
normalization. Representative 1-D totals include `37/32` and `707/800`, so integrated
gas momentum and energy feedback do not necessarily equal the opposite particle change.

The exact shock conservation ledger correctly rejects AMR/SMR in
`src/pgen/tests/pic_parallel_shock.cpp:4178-4243`; this protection must remain until the
policy is resolved.

Required decision:

- [ ] Implement a conservative cross-level partition/reflux/adjoint scheme; **or**
- [ ] Retain `paper_smooth` as an accepted non-conservative model with explicit error
  bounds and prohibit exact-conservation claims.

Required validation:

- [ ] Check partition of unity for every interface position and supported dimension.
- [ ] Check total gas plus particle momentum and energy closure.
- [ ] Test dynamic AMR, SMR, MPI decomposition, restart, and interface corners.
- [ ] Compare Bell growth and shock precursor structure with a matched uniform grid.

### STAB-1 — No feedback source-timescale or positivity constraint

- Status: **OPEN**
- Severity: **P1**
- Evidence:
  - feedback update: `src/mhd/mhd_tasks.cpp:432-443`
  - particle timestep: `src/particles/particles.cpp:1909-2034`
  - energy floor: `src/eos/ideal_c2p_mhd.hpp:43-56`

The final gas update subtracts the full deposited particle impulse without checking that
the resulting gas state has positive internal energy. FOFC only repairs flux updates and
cannot preempt source-driven failures. Subsequent EOS floors inject energy and invalidate
exact conservation.

Required work:

- [ ] Form a trial post-feedback gas state.
- [ ] Add a source CFL, controlled subcycling, or a fail-fast admissibility threshold.
- [ ] Record source-induced density, pressure, energy, and temperature floor events.
- [ ] Make conservation-qualified runs abort on any such event.
- [ ] Stress-test high macro-particle loading, low beta, strong currents, and shock cells.

### DT-1 — Gyro bound combines only rank-local extrema

- Status: **OPEN**
- Severity: **P1**
- Evidence: `src/particles/particles.cpp:1976-2030` and `src/mesh/mesh.cpp:641-650`

A rank-local maximum particle `|q/mc|` is multiplied by a rank-local maximum magnetic
field. Only the resulting timestep is globally minimized. A particle and a strong field
on different ranks therefore never form the required worst-case product.

Required work:

- [ ] Use the maximum configured species `|q/mc|`, or reduce the extrema separately.
- [ ] Include any ghost-field region that can be sampled during the next half drift.
- [ ] Add a two-rank high-q/m/strong-B boundary-crossing regression.

### DT-2 — First injected shock cohort can bypass particle timestep limits

- Status: **OPEN**
- Severity: **P1**
- Evidence:
  - empty-particle early return: `src/particles/particles.cpp:1929-1931`
  - callback ordering: `src/driver/driver.cpp:507-520,617-618`
  - injection append: `src/pgen/tests/pic_parallel_shock.cpp:3292-3350,4000-4012`

The timestep is chosen before the next cycle's injection. If no particles existed during
selection, the first cohort is pushed using an MHD-only timestep.

Required work:

- [ ] Include configured injection velocity and q/m bounds in timestep selection; or
- [ ] Recompute/validate the timestep after injection and before the stage-1 push.
- [ ] Add an intentionally restrictive first-cohort test.

### SHOCK-1 — Incompatible startup injection histories

- Status: **OPEN POLICY DECISION**
- Severity: **P1**
- Evidence:
  - particle removal: `src/pgen/tests/pic_parallel_shock.cpp:2142-2213`
  - production-style t=0 injection and t=45 removal:
    `inputs/publication/pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput:134-140`
  - delayed static injection: `inputs/q011_section54_static_dx3_final_v1_vl2_tsc.athinput:157-164`

Deleting early particles does not restore gas mass, momentum, or energy already removed
at injection. A run that injects at t=0 and deletes the cohort at t=45 is therefore not
equivalent to a run that begins injection at t=45.

Required work:

- [ ] Freeze one science policy; delayed injection after shock formation is preferred.
- [ ] Run a matched startup-history comparison.
- [ ] Prevent cross-comparison of histories unless the resulting shock-state difference
  is explicitly quantified.

### SHOCK-2 — Injection follows an analytic rather than measured shock surface

- Status: **OPEN / CONTROLLED-MODE LIMITATION**
- Severity: **P1**
- Evidence: `src/pgen/tests/pic_parallel_shock.cpp:736-769,2758-2783,3003-3018`

Particle placement and injected mass use an analytic planar surface, fixed model speed,
fixed upstream density, and planar area. This is suitable for a short Section 5.4-style
controlled reproduction but not for a corrugated or CR-modified nonlinear shock.

Required work:

- [ ] Track the local and area-averaged shock surface.
- [ ] Measure upstream surface-normal mass flux.
- [ ] Gate injection validity on alignment with the detected front.
- [ ] Report the effective injected fraction using measured swept mass.
- [ ] Compare analytic-surface and tracked-surface pilots.

### SHOCK-3 — Paper mode omits CR-Hall induction physics

- Status: **ACCEPTED LIMITATION, NEEDS DIAGNOSTIC**
- Severity: **P1** when the Hall parameter is not small
- Evidence: `src/particles/particles.hpp:554-559` and
  `src/particles/particles.cpp:1302-1312`

Paper mode correctly reproduces the ideal-MHD induction choice and rejects the
experimental direct-current CT path. The resulting shock model is physically applicable
only while the CR-Hall correction remains demonstrably small.

Required work:

- [ ] Define and output local/volume/surface Hall applicability measures.
- [ ] Freeze an acceptance threshold.
- [ ] Stop or mark runs nonqualifying if the threshold is exceeded.
- [ ] Keep the experimental Hall extension outside paper-mode claims until independently
  validated.

### BELL-1 — Legacy nonphysical current normalization remains runnable

- Status: **OPEN**
- Severity: **P1**
- Evidence:
  - legacy relation: `src/pgen/tests/q023_paper_bell_linear.cpp:121-136`
  - corrected relation: `src/pgen/tests/q043_bell_current_volume_aware.cpp:88-116`
  - supersession record:
    `tst/publication/readiness/q043_bell_current_normalization_supersession_2026-06-06.json`

The legacy generator omits root-cell volume and incorrectly includes the artificial CR
light speed in the deposited-current target. The repository documents the invalidation,
but the generator and decks remain directly runnable.

Required work:

- [ ] Require an explicit `legacy_nonphysical_mechanics_only` opt-in, or reject legacy
  generators in science launch tooling.
- [ ] Ensure every active Bell deck uses the volume-aware `J_CR/c` relation.
- [ ] Print derived deposited current and root-cell volume at startup.

### BELL-2 — Corrected Bell mechanics lack qualifying runtime evidence

- Status: **OPEN**
- Severity: **P1**
- Evidence:
  - Q043 volume-aware source and oracle
  - Q019 deck fields `q043_independent_raw_cycle_one_oracle_bound=false`,
    `q023_independent_linear_predecessor_bound=false`, and `launch_authorized=false`
  - Q019 checked-in manifest currently drifts from generated content

The corrected Q043 and current Q019 source arithmetic are well designed, but host harnesses
that include the production source are self-consistency tests rather than independent
runtime evidence.

Required work:

- [ ] Complete the Q043 deposited-current matrix across dimension, resolution, PPC,
  decomposition, and artificial light speed.
- [ ] Independently recompute `J_CR/c` from raw particle and grid outputs.
- [ ] Complete corrected linear Bell growth, wavelength, phase, and polarization tests.
- [ ] Close convergence gates before nonlinear saturation work.
- [ ] Regenerate and review manifests only after the source and deck set are frozen.

### PERF-1 — AMR paper deposition abandons GPU residency

- Status: **OPEN**
- Severity: **P2**, potentially run-limiting
- Evidence: `src/particles/particles_moments.cpp:662-868` and
  `src/bvals/bvals_mom.cpp:142-196`

Each deposition stage mirrors all particle data to the host, builds receiver records in a
serial host loop, performs global collectives, reallocates a device record array, and copies
records back. Paper mode does this twice per timestep.

Required work:

- [ ] Replace global exchange with sparse neighbor/device-aware transport.
- [ ] Build and compact records on device.
- [ ] Retain capacity rather than reallocating exact sizes every stage.
- [ ] Benchmark strong/weak scaling on the intended Frontier topology.

### PERF-2 — Ordinary migration uses global metadata and allocation churn

- Status: **OPEN**
- Severity: **P2**
- Evidence: `src/bvals/bvals_part.cpp:618-645,680-797,909-913,976-1038`

Migration globally gathers sender descriptors, repeats the global particle-count gather,
reallocates full-size send/destroy lists and exact message buffers, and uses many small
deep copies during compaction.

Required work:

- [ ] Use sparse handshakes or neighbor collectives.
- [ ] Remove the redundant particle-count allgather.
- [ ] Introduce retained geometric-capacity buffers.
- [ ] Use a single device compaction kernel shared with `MIG-1`.

### PERF-3 — Shock injection is globally replicated and repeatedly reallocates

- Status: **OPEN**
- Severity: **P2**
- Evidence: `src/pgen/tests/pic_parallel_shock.cpp:2835-2851,3060-3111,3292-3348`

Every rank receives the full shock-cell list, loops over the global injected population,
and grows particle arrays to the exact new count each cycle.

Required work:

- [ ] Use distributed prefix selection and direct owner construction.
- [ ] Maintain particle-array capacity separately from live particle count.
- [ ] Benchmark injection cost versus total particle count and rank count.

### PERF-4 — Advertised sorting and load-balancing controls are ineffective

- Status: **OPEN**
- Severity: **P2**
- Evidence:
  - `pic_sort_interval` parse only: `src/particles/particles.cpp:749-755`
  - uniform-mesh balancing exclusion: `src/driver/driver.cpp:642-647`

`pic_sort_interval` is parsed and restart-fingerprinted but never used. Bell and shock decks
set it nonzero. Particle-weighted balancing is also inactive for uniform
`refinement=none` Bell decks even when their cost is nonzero.

Required work:

- [ ] Implement stable device sorting/binning by MeshBlock/cell/tag on the requested cadence,
  or reject nonzero values.
- [ ] Document the mesh classes on which particle-aware balancing is active.
- [ ] Decide whether periodic redistribution is needed for uniform particle runs.

### ROBUST-1 — Invalid particle state can be silently accepted or rewritten

- Status: **OPEN**
- Severity: **P2**
- Evidence: `src/particles/particles.cpp:219-225,412-475,572-645,1589-1613`,
  `src/particles/particles_pushers.cpp:507-515`, and
  `src/particles/particles_moments.cpp:702-705,1012-1016,1066-1070`

Several physical inputs lack complete finite/positivity checks, and nonpositive particle
weights are silently replaced by one in paper push/deposition paths.

Required work:

- [ ] Centralize finite and admissibility validation for configuration, initialization,
  injection, migration, and restart.
- [ ] Fail closed on invalid weights, species, GID, state, or macro mass.
- [ ] Add an optional per-cycle debug state audit.

### ROBUST-2 — Rank-local exits can strand MPI peers

- Status: **OPEN**
- Severity: **P2**
- Evidence: error paths in `src/particles/particles_moments.cpp:690-818` and
  `src/bvals/bvals_part.cpp:680-1089`

Several locally detected errors call `std::exit` before other ranks reach collective
operations. Use collective validation where appropriate and the repository's MPI-aware
fatal path for unrecoverable failures.

Required work:

- [ ] Audit PIC and turbulence fatal paths.
- [ ] Replace unsafe local exits with collective checks or `AbortOnFatalError()`.
- [ ] Add an injected-error MPI test that terminates cleanly rather than hanging.

### ROBUST-3 — `TaskStatus::fail` can produce an infinite busy loop

- Status: **OPEN**
- Severity: **P2**
- Evidence: `src/tasklist/task_list.hpp:144-157` and `src/driver/driver.cpp:369-385`

Only `complete` is handled specially. A task returning `fail` remains runnable forever,
and the driver has no progress or timeout detection.

Required work:

- [ ] Propagate `fail` to the driver and terminate through the MPI-aware fatal path.
- [ ] Add no-progress detection for task lists.
- [ ] Test a manufactured boundary-communication failure.

### ROBUST-4 — Turbulence MPI reductions assume double precision

- Status: **OPEN**
- Severity: **P2**
- Evidence: `src/srcterms/turb_driver.cpp:890-895,943-946,1228-1233`

`Real` arrays are reduced using `MPI_DOUBLE`. Single-precision builds can read and write
beyond the buffers.

Required work:

- [ ] Replace `MPI_DOUBLE` with `MPI_ATHENA_REAL`.
- [ ] Build and run the turbulence unit tests in single and double precision.

## Validation baseline from the review

These results describe the reviewed baseline and should be updated as fixes land.

| Validation | Baseline result | Interpretation |
|---|---:|---|
| Paper-smooth topology/raw-TSC tests | 11 passed | Confirms current routing and intentionally non-unit interface totals; not a conservation qualification |
| Q043 corrected-current host tests | 12 passed | Strong source-local arithmetic check; not independent runtime evidence |
| Q019 design tests | 23 passed, 1 failed, 5 skipped | Checked-in manifest drift |
| Q023 host tests | 3 passed, 1 failed | Stale artifact hash |
| Q011 shock preparation tests | 1 passed, 4 failed | Stale source-fragment contracts and hashes |
| Q019 checked-in deck CLI validation | failed | Checked-in manifest drifted |
| Fresh serial Debug build (`PROBLEM=turb`) | passed | Full executable build and `-c` startup identity |
| Fresh MPI-enabled Debug build (`PROBLEM=turb`) | passed | Full executable build plus true two-rank MHD and MHD-PIC runs |
| Blocker host tests | 4 passed | Stable root, RK energy maps, signed modes/projection, and exhaustive compaction identity |
| Focused AthenaK regression harness | 2 of 2 passed | Fresh serial build ran turbulence and migration/destruction smokes |
| Turbulence runtime smoke | serial/MPI2 passed | 3-D MHD-PIC/MHD-only budgets, 2D3V, force FFT, restart parity, and legacy-bound upgrade |
| Migration/destruction runtime smoke | serial/MPI2 passed | Exact 180-of-192 survivor inventory and full state/metadata parity |

Additional test-infrastructure issue:

- `tst/run_tests.py` discovers every non-helper module and calls `run()` and `analyze()`.
  `tst/scripts/particles/pic_paper_smooth_tsc_oracle.py` is discovered but implements
  neither function. The default all-suite GPU CI entry point is therefore structurally
  broken until the oracle is renamed/routed as a helper or given the regression interface.

## Strengths and invariants to preserve

The following should be protected with focused regression tests while addressing the
open findings:

- Paper-mode fail-closed composition guards for coupled feedback, TSC order, ghost width,
  RK2, source placement, ideal induction, and unit conservative coefficients.
- Explicit-midpoint/VL2 particle chronology on uniform grids.
- Relativistic Boris rotation using the correct midpoint gamma.
- Stable relativistic kinetic-energy evaluation.
- Macro-particle mass, charge, and q/m semantics with root-cell-volume normalization.
- Matching uniform-grid TSC gather and deposition support.
- Exact opposite full-f particle/gas momentum and kinetic-energy deltas on uniform grids.
- Explicit 2D3V handling and rejection of unsupported one-dimensional particle runs.
- Deterministic shock tags and velocity sampling.
- Once-per-cycle shock particle creation with RK-weighted gas transaction replay.
- Restart fingerprints, ownership validation, boundary ledgers, and overflow checks.
- Corrected Q043/current-Q019 volume-aware Bell current relation.
- Existing fail-closed labels that distinguish preparation, engineering, and qualifying
  science artifacts.

## Decisions and update log

Record substantive decisions and validation evidence here. Link commits, tests, run
directories, or retained evidence as appropriate.

| Date | Item | Status change | Decision or evidence |
|---|---|---|---|
| 2026-07-08 | Initial review | Created | Static/source review at `f8a56172983a1f009662af4217dcb7b52032d7dc`; no implementation edits or fresh binary run |
| 2026-07-08 | `TURB-1`–`TURB-4` | OPEN -> FIXED, NEEDS VALIDATION | Stable normalization, RK-consistent work, dimension-aware signed modes, strict isotropic bounds, legacy-restart migration, host tests, and serial/MPI2 runtime smokes pass; the remaining item-specific closure tests stay open |
| 2026-07-08 | `MIG-1` | OPEN -> FIXED, NEEDS VALIDATION | Unified validated compaction and one batched survivor kernel; exhaustive oracle plus exact serial/MPI2 survivor-state parity pass; restart/repeated-crossing/ledger evidence remains |

## Final deletion checklist

Delete this document only when all of the following are true:

- [ ] Every P0 item is **CLOSED**.
- [ ] Every P1 item is **CLOSED** or has a reviewed **ACCEPTED LIMITATION** with runtime
  guards and claim boundaries.
- [ ] Required P2 reliability issues are closed; deferred performance work is captured
  in permanent project tracking rather than silently abandoned.
- [ ] The complete particle regression suite has a green CPU and GPU run.
- [ ] MPI, restart, AMR/SMR, and single/double-precision targeted tests are green.
- [ ] Turbulent-box energy injection, divergence, isotropy, and paired MHD/PIC tests pass.
- [ ] Mixed migration/escape survivor-identity tests pass.
- [ ] Uniform and selected multilevel gas-plus-particle conservation tests pass under the
  adopted AMR policy.
- [ ] Corrected Bell raw-current and linear-physics matrices pass independently.
- [ ] Shock injection startup, tracked-surface, and applicability tests pass.
- [ ] Existing affected campaign outputs have been rerun or clearly retired.
- [ ] Scientific claims and figure provenance reference only qualified replacement runs.
- [ ] The update log contains the final qualifying commit and evidence locations.
