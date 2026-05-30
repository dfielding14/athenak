---
orphan: true
---

# CR Tracer Relativistic Acceleration Handoff

## Verdict Boundary

This branch implements an explicit opt-in **passive** relativistic cosmic-ray
tracer model.  It advances full-orbit test particles with authoritative
`w = gamma v` state and the Higuera-Cary electromagnetic map.  It does not feed
particle momentum or energy back into the fluid.

The workstation qualification verdict is:

- `CPU/MPI QUALIFIED`: yes, within the bounded scope below;
- `GPU QUALIFIED`: no, unavailable locally and explicitly unclaimed;
- `MERGE READY`: pending the final Phase-9 public-overlay privacy seal,
  `RG-010`, and `CP-7 Final Cold Review`.  The strict public documentation
  overlay has passed as a development-preview render, not as a stable-site
  publication.

## Branch Facts

| Item | Value |
| --- | --- |
| Branch | `feature/CR_tracers_relativistic_acceleration` |
| Required prerequisite | `feature/CR_tracers_followup_architecture` |
| Frozen prerequisite SHA | `64a4d1be8da1c22d1328cc47280195b3747fa0ab` |
| Intended integration target | `origin/development` |
| Refreshed integration-target SHA | `c6a73b08e60807f8b925164c5e7edd5cb820c8ae` |
| Accepted Phase-8 source SHA | `addd12d4e26f7d8b275165b6be7b364d39f22a43` |
| Accepted Phase-8 evidence seal | `aa8663e8a5d49e26c206363d028b52d0e350a91f` |
| Historical removed Phase-9 tips | `720fb8fc`, `89084c66`, `fa39416a`, and `5e031387` were replaced after cold-review holds. |
| Sanitized Phase-9 pre-seal candidate SHA | `5c2d39b5c36a48e3b552eb14406a039edf1d29f8` |
| Pushed Phase-9 portability-rebound seal SHA | `aa66f6c27531116e12554631281c8f2ed07d93c6` |
| Documentation-durability correction candidate SHA | `9a270755f2739398d61023fa7a950add2dd550e0` |
| Pushed documentation-durability seal SHA reviewed on rebound | `88d631e4943648fe83f0624cb30291fa52ab4296` |
| Pushed capture-completeness seal SHA held on cold review | `1a7086add5fffd55356109b99e6a66fcd0b43486` |
| Final public-overlay privacy review target | Resolve and record the exact pushed branch-tip SHA after this documentation-only reseal.  A commit cannot embed its own SHA; `CP-7` must name the reviewed immutable tip explicitly. |
| Next reviewed packet range | `64a4d1be8da1c22d1328cc47280195b3747fa0ab..refs/remotes/origin/feature/CR_tracers_relativistic_acceleration` resolved after the final documentation-only reseal |
| Decision ledger | `cr_tracer_relativistic_acceleration_ledger.md` |

## Physical Model

The reduced-model relations are:

```text
gamma = sqrt(1 + |w|^2 / c_model^2)
v = w / gamma
kinetic_energy_model = (gamma - 1) c_model^2
```

For the solver-coupled selector, the code trilinearly gathers Newtonian
cell-centered fluid velocity and magnetic field and forms:

```text
cE = -u_fluid x B
```

Each Higuera-Cary substep updates `w` and accumulates direct sampled-field work:

```text
alpha_s dt cE dot (w_old + w_new) / (gamma_old + gamma_new)
```

`c_model` is an explicit finite positive code-unit model speed.  `alpha_s` is an
explicit finite nonzero signed normalized force coefficient.  The sign is not
the legacy `IPM < 0` alignment sentinel.

## Qualified Scope

The qualified solver-coupled path is deliberately narrow:

- `pusher = relativistic_hc`;
- `relativistic_field_source = mhd_ideal`;
- `relativistic_temporal_sampling = frozen_tn`;
- Newtonian coordinates;
- ideal MHD;
- explicit trilinear interpolation;
- 3-D strictly periodic boundaries;
- serial uniform-level exactly-one-MeshBlock execution;
- `c_model = 1.0`;
- passive particles only.

The qualified prescribed-test path is mechanically separate:

- `relativistic_field_source = prescribed_test`;
- uniform prescribed `B` and `cE`;
- analytical acceleration, deceleration, drift, orbit, convergence, and work
  closure;
- periodic MPI rank ladder `1`, `2`, `4`, and `8`;
- same-level multiblock migration;
- static SMR;
- adaptive refine/derefine;
- both `device_table` and `host_tree` multilevel lookup implementations;
- same-topology restart continuation before and after migration;
- deterministic rank-count-change rejection.

The prescribed-test migration opening proves infrastructure behavior.  It does
not widen solver-coupled sampled-MHD execution.

## State And Restart

Legacy `drift` and `boris` particles retain their 14-real-field and
3-integer-field layout.  `relativistic_hc` appends:

| Field | Meaning |
| --- | --- |
| `IPWX`, `IPWY`, `IPWZ` | Authoritative `w = gamma v`. |
| `IPCEX`, `IPCEY`, `IPCEZ` | Sampled or prescribed `cE`. |
| `IPWORK` | Accumulated direct sampled-field work. |
| `IPALPHA` | Signed normalized force coefficient. |

Relativistic restart uses typed-v2 particle shards and a manifest paired with
the native mesh restart and `.rst.rmeta` witness.  The bundle binds cycle, time,
timestep, rank count, topology, byte counts, checksums, and checkpoint
identity.  Changed-rank continuation rejects with a deterministic diagnostic
until redistribution support is designed and reviewed separately.

## Diagnostics

The retained relativistic diagnostic schema is deliberately explicit and
narrow.  Headers declare:

- `schema=akcr_particle_output_v1`;
- `mode=relativistic_hc`;
- `units=code_model`;
- `c_model`;
- `alpha_s`;
- field-specific units and layouts.

Relativistic `pspec`, `pspec2`, and `psamp` require canonical field names.
Ambiguous legacy aliases reject rather than silently changing meaning.
Field-relative diagnostics fail closed when sampled `|B| = 0`.

The selected-sample output can expose `wx`, `wy`, `wz`, `wmag`, `gamma`,
`kinetic_energy_model`, `cex`, `cey`, `cez`, `work`, `alpha_s`, and
`r_larmor_over_dx_min`.

## Accepted Qualification

Phase-8 migration evidence is sealed under:

```text
docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase8_*
```

The accepted source-bound migration package records:

| Gate | Result |
| --- | --- |
| Serial and MPI builds | Passed |
| Frozen migration registration | `62/62` |
| Deliberate oracle mutations | `17/17` rejected |
| Isolated migration runtime cases | `21` |
| MPI rank ladder | `1`, `2`, `4`, `8` |
| Empty-particle-rank witness | Passed |
| Static-SMR device-table cross-level witness | Passed |
| Static-SMR host-tree cross-level witness | Passed |
| Adaptive AMR remap modes | `device_table`, `host_tree` |
| Restart-after-AMR continuations | Exact parity for both remap modes |
| Same-topology continuation before and after migration | Exact parity |
| Changed-rank continuation | Exact deterministic rejection |
| Parser contract | `172/172` |
| Legacy CPU suites | `20/20`, `12/12` |
| Legacy MPI suites | `4/4`, `5/5` |
| Durable raw replay tarball | Extracted and decoded offline |
| Retained typed-v2 checkpoints decoded offline | `29` |
| Path-sorted SHA-256 evidence manifest | Verified |

Earlier sealed phases retain the helper, parser, gather, pusher, coupled-field,
subcycle, restart, and diagnostic evidence.  Phase 9 replays those analyzers on
the final branch candidate and archives the refreshed outputs.

## Phase-9 Evidence Index

The final workstation evidence lives under:

```text
docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/
```

| Surface | Retained artifact |
| --- | --- |
| Complete replay command inventory | `phase9_qualification_command_log.md` |
| Candidate provenance snapshot | `phase9_branch_snapshot.txt` |
| Source-bound provenance envelope | `phase9_provenance_envelope.json` |
| Path-sorted evidence manifest | `phase9_evidence_sha256.txt` |
| Manifest verification | `phase9_evidence_verify.log` |
| Relocated-root manifest verification | `phase9_evidence_relocated_verify.log` |
| Serial release, MPI release, debug-bounds, MPI debug-bounds builds | `phase9_release_*`, `phase9_mpi_release_*`, `phase9_debug_*`, `phase9_mpi_debug_*` |
| Dedicated pusher and coupled builds | `phase9_pusher_release_*`, `phase9_coupled_release_*`, `phase9_coupled_debug_*` |
| Analytical replays and metrics | `phase9_pusher_runtime*`, `phase9_coupled_runtime*`, `phase9_subcycle_runtime*` |
| Restart and output replays | `phase9_restart_runtime*`, `phase9_diagnostics_runtime*`, `phase9_all_formats*` |
| Reader and writer adversarial probes | `phase9_*negative_controls.log`, `phase9_diagnostics_writer_adversarial.log` |
| Release and debug-bounds migration qualification | `phase9_migration_runtime.log`, `phase9_migration_debug_runtime.log`, and paired metrics |
| Repository parser and legacy harness | `phase9_parser_contract.log`, `phase9_legacy_*.log`, `phase9_style.log` |
| Deep AMR `divB` regression | `phase9_deep_amr_divb.log` |
| Strict public documentation overlay | `phase9_gh_pages_overlay_build.log` |
| Frozen public baseline and overlay inventory | `phase9_gh_pages_overlay_inventory.txt` |
| Executed public solver-coupled example | `phase9_public_solver_coupled_example_runtime.log` |
| Complete public solver-coupled input | `inputs/particles/cr_tracer_relativistic_mhd_ideal_example.athinput` |
| Tracked typed-v2 restart resume template | `inputs/particles/cr_tracer_relativistic_prescribed_restart_resume.athinput` |
| Executed typed-v2 restart resume template | `phase9_public_restart_resume_template_runtime.log` |
| Commit-range whitespace audit | `phase9_commit_range_diff_check.log` |
| Retained-artifact privacy audit | `phase9_artifact_privacy_scan.log` |
| Provenance-envelope verification | `phase9_provenance_envelope_verify.log` |
| Workstation accelerator disposition | `phase9_gpu_disposition.log` |
| Unsupported single-precision build attempt | `phase9_single_precision_build.log` and cache |
| Refreshed development merge-tree audit | `phase9_merge_tree_origin_development.log` |

The deterministic archived-artifact replay is retained as
`phase9_all_formats_portable_replay_bundle.tar.gz`.  It includes the qualified
release binary, runtime artifacts, input, criteria, required inspector scripts,
a deterministic portable manifest, and a replay wrapper.  Fresh extraction
and strict reinspection are recorded in
`phase9_all_formats_portable_replay.log`.

## Affected Source Inventory

The implementation is intentionally distributed across a bounded set of
surfaces:

- `src/athena.hpp`: appended relativistic particle-state indices;
- `src/particles/`: parser fences, state helpers, Higuera-Cary push, gather,
  subcycle refresh, migration lookup diagnostics, and restart codec;
- `src/pgen/part_random.cpp`: explicit relativistic initialization and typed-v2
  continuation application;
- `src/outputs/`: paired typed-v2 restart publication and versioned
  relativistic diagnostics;
- `scripts/particles/`: inspectors, deterministic analyzers, and adversarial
  controls;
- `inputs/unit_tests/` and `tst/test_suite/particles/`: bounded runtime decks
  and parser or legacy regression coverage;
- `docs/source/modules/`: public parameter contract, implementation guide,
  living ledger, and handoff documents.

## Integration Status

The refreshed merge-tree audit against
`origin/development@c6a73b08e60807f8b925164c5e7edd5cb820c8ae` reports manual
integration conflicts in:

```text
src/main.cpp
src/mesh/mesh.cpp
src/mesh/mesh.hpp
src/outputs/outputs.cpp
src/outputs/outputs.hpp
src/outputs/restart.cpp
src/pgen/pgen.cpp
```

Do not merge this stacked branch directly into `development` without a
dedicated reviewed integration step.  Preserve the prerequisite branch first,
resolve each overlap explicitly, rerun the complete qualification matrix, and
repeat the cold review on the integrated candidate.

## Known Limitations

- No GPU qualification is claimed.  Use
  [CR Tracer Relativistic GPU Testing Handoff](cr_tracer_relativistic_gpu_testing_handoff.md).
- Solver-coupled `mhd_ideal` MPI, SMR, and AMR execution remains fail-closed.
- Nonperiodic boundaries remain fail-closed.
- Changed-rank typed-v2 restart redistribution remains unsupported.
- Coupled temporal sampling is `frozen_tn`; stage-coupled or time-centered
  widening is deferred.
- The model is passive and full-orbit.  It is not MHD-PIC, SRMHD, guiding
  center, or scientific production-readiness certification.
- Single-precision qualification remains blocked by inherited repository
  narrowing errors.
- Multi-node scaling and performance optimization remain follow-up work.

## Deferred Follow-Ups

Keep follow-ups separate and bounded:

| Proposed branch | Entry gate | Scope |
| --- | --- | --- |
| `feature/CR_tracers_relativistic_acceleration_gpu_qualification` | Accepted handoff SHA on a configured CUDA or HIP runner | GPU backend qualification and profiling only. |
| `feature/CR_tracers_relativistic_coupled_mpi_amr` | Accepted GPU-independent handoff plus new preregistration | Solver-coupled MPI, SMR, and AMR widening only. |
| `feature/CR_tracers_relativistic_stage_coupling` | Accepted coupled MPI/AMR branch plus temporal oracle | Stage-coupled temporal sampling only. |
| `feature/CR_tracers_relativistic_restart_redistribution` | Accepted typed-v2 same-topology handoff plus redistribution design | Changed-rank typed-v2 restart redistribution only. |
| `feature/CR_tracers_relativistic_nonperiodic_boundaries` | Accepted boundary semantics decision | Nonperiodic boundary design only. |
| `feature/CR_tracers_guiding_center_study` | Full-orbit applicability trigger reached | Guiding-center applicability study only. |
| `feature/CR_tracers_feedback_model_design` | Separate scientific requirement and model review | Passive-to-feedback model design only. |
| `feature/CR_tracers_relativistic_acceleration_integration` | Current `development` tip frozen and conflict plan reviewed | Reconcile the stacked prerequisite and this branch onto `development`, then rerun the full matrix. |

## Review Rule

Every follow-up must preserve the decision-ledger discipline used here:

1. Record the choice and rejected alternatives.
2. Freeze acceptance criteria before widening source.
3. Stop on failed or ambiguous evidence.
4. Ask fresh independent reviewers to try to falsify the claimed milestone.
5. Recheck fixes before recording `PROCEED`.
