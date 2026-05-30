---
orphan: true
---

# CR Tracer Relativistic Acceleration Phase 0 Ledger

## Ledger Status

This is the Phase 0 implementation ledger required by
`cr_tracer_relativistic_acceleration_implementation_guide.md`.  It freezes the
audited baseline, opens the mandatory decision records, and records missing
evidence without authorizing production-code edits.

The statuses in this initial draft are intentionally conservative:

- `selected`: a Phase 0 working selection defensible from the current code and
  guide.  It is not an accepted checkpoint and does not open a production edit
  set by itself.
- `proposed`: deliberately unresolved pending a spike, derivation, test,
  independent review, or migration audit.

## Control Header

| Control | Current value | Update rule |
| --- | --- | --- |
| Ledger date | `2026-05-30` | Record the date of each accepted update. |
| Working branch | `feature/CR_tracers_relativistic_acceleration` | Stop if the working branch changes unexpectedly. |
| Working branch HEAD at latest accepted commit | `be24fa12ca9b330f1da20c489c826a231bf01d85` | Refresh after every accepted commit. |
| Frozen feature base | `64a4d1be8da1c22d1328cc47280195b3747fa0ab` from `feature/CR_tracers_followup_architecture` | Change only through an accepted `DR-000` update. |
| Intended eventual integration target | `origin/development` at `c6a73b08e60807f8b925164c5e7edd5cb820c8ae` | Refresh the target SHA and merge-tree audit after target updates and before handoff. |
| Current implementation phase | `Phase 3: shared trilinear gather helper and prescribed-field diagnostics` | Advance to acceleration only after helper and diagnostic-only gather tests, independent review, and `RG-004` record `PROCEED`. |
| Allowed write manifest | `src/particles/trilinear_gather.hpp`, `src/particles/particles_pushers.cpp`, narrowly scoped Phase-3 diagnostic pgen and input files, gather-only tests, Phase-3 evidence, and this ledger only | Keep particle widths, restart payloads, outputs, MHD tasks, AMR algorithms, subcycling, displacement, Boris algebra, and positive-cycle `relativistic_hc` abort unchanged. |
| Last accepted checkpoint | Phase-2 fail-closed parser contract at `be24fa12ca9b330f1da20c489c826a231bf01d85` | Record each milestone commit before opening the next phase. |
| Open blocking findings | Phase-3 helper review corrections, true MPI launch parity, diagnostic-only runtime integration, and actual ghost-fill, periodic, block-boundary, SMR, forced-AMR evidence; downstream `BF-009`, `BF-014`, backend qualification, layout, restart, output, solver-coupling, and migration obligations remain open for their owning phases | Keep the live list current.  Do not proceed while a blocker for the next edit set remains open. |
| CPU/MPI qualification status | Phase 2 accepted: focused parser `49 passed` plus `nlim=0` smoke and positive-cycle abort; full CPU `269 passed, 15 skipped, 67 deselected`; MPI CPU `38 passed, 313 deselected`; style `2 passed`; single precision remains blocked by inherited repository narrowing errors | Rerun and archive qualification at each accepted milestone. |
| GPU qualification status | Not started; unavailable on this workstation | GPU qualification is optional for `MERGE READY` on this workstation, but it remains a separate unqualified residual risk.  Do not claim `GPU QUALIFIED` without accelerator evidence on an accepted SHA. |
| Merge-ready status | Not eligible | Requires the guide's checkpoint sequence, including `CP-4 Solver Coupling` and `CP-7 Final Cold Review`. |

## Historical Branch And Base Facts At Ledger Creation

The labels in this original Phase-0 table are historical snapshots at ledger
creation.  Use the refreshed control header and the latest bound
`evidence/phase*_branch_snapshot.txt` archive for live checkpoint provenance.

| Fact | Recorded value | Evidence |
| --- | --- | --- |
| Current branch | `feature/CR_tracers_relativistic_acceleration` | `git branch --show-current` |
| Current local HEAD | `3984b57ec8832b3e3e8140e08188a8a05e773136` | `git rev-parse HEAD` |
| Current remote tracking tip | `origin/feature/CR_tracers_relativistic_acceleration` at the same SHA | `git status --short --branch` |
| Current HEAD purpose | Guide-only commit: `docs: add relativistic CR acceleration implementation guide` | `git show --no-patch --format=... HEAD` |
| Required prerequisite branch | `feature/CR_tracers_followup_architecture` | Guide control header |
| Frozen prerequisite tip | `64a4d1be8da1c22d1328cc47280195b3747fa0ab` | Parent of current HEAD and `git merge-base HEAD 64a4d1be` |
| Integration target | `origin/development` | `DR-000` working selection |
| Integration-target SHA | `c6a73b08e60807f8b925164c5e7edd5cb820c8ae` | `git rev-parse origin/development` |
| Merge base with integration target | `c3cfc63fba22ae3794a0b87fdccca2c1c5a3b883` | `git merge-base HEAD origin/development` |
| Dirty state at first audit command | Clean | Initial `git status --short --branch` |
| Dirty state later in the same audit | Unexpected untracked files appeared; see `BF-001` | `git status --porcelain=v1 --untracked-files=all` |

The unexpected untracked files were not removed or modified during this ledger
draft:

```text
docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
scripts/particles/__pycache__/cr_relativistic_pusher_spike.cpython-313.pyc
scripts/particles/__pycache__/cr_tracer_inspect.cpython-313.pyc
scripts/particles/cr_relativistic_pusher_spike.py
tst/test_log.txt
tst/test_suite/__pycache__/__init__.cpython-313.pyc
tst/test_suite/__pycache__/testutils.cpython-313.pyc
tst/test_suite/particles/__pycache__/__init__.cpython-313.pyc
tst/test_suite/particles/__pycache__/cr_accuracy_utils.cpython-313.pyc
tst/test_suite/particles/__pycache__/test_particles_cr_accuracy_mpicpu.cpython-313-pytest-8.4.1.pyc
tst/test_suite/particles/__pycache__/test_particles_cr_mpicpu.cpython-313-pytest-8.4.1.pyc
vis/python/__pycache__/athena_read.cpython-313.pyc
```

The standalone pusher-spike script and JSON result appeared after the first
ledger draft during read-only verification.  They are concurrent, unreviewed
candidate artifacts.  This ledger does not claim ownership of them and does
not treat them as accepted Phase 0 evidence.

## Merge-Tree Audit Against `origin/development`

The required read-only audit was run against the selected target SHA:

```bash
git merge-tree --write-tree \
  3984b57ec8832b3e3e8140e08188a8a05e773136 \
  c6a73b08e60807f8b925164c5e7edd5cb820c8ae
```

The command returned exit status `1` and synthetic tree
`a29705847e632b98012f835ad5cd3014fcad4730`.

### Content Conflicts

| Path | Required reconciliation focus |
| --- | --- |
| `src/mesh/mesh.cpp` | Particle integration hooks versus development integration repairs and merged features. |
| `src/mesh/mesh.hpp` | Particle state owned by `Mesh` versus development mesh additions. |
| `src/outputs/outputs.cpp` | Particle output registration versus development output integration. |
| `src/outputs/outputs.hpp` | Particle output declarations versus development output declarations. |
| `src/pgen/pgen.cpp` | Particle random problem-generator registration versus development pgen integration. |

### Additional Same-Ancestry Overlaps

The following paths changed on both sides of the merge base and currently
merge without content conflicts.  They remain mandatory review and retest
surfaces:

| Path | Required reconciliation focus |
| --- | --- |
| `src/CMakeLists.txt` | Pgen and source registration. |
| `src/mesh/mesh_refinement.cpp` | Particle-aware AMR load balancing versus development refinement repairs. |
| `src/outputs/basetype_output.cpp` | Particle output exemptions versus development output handling. |
| `src/pgen/pgen.hpp` | Problem-generator declarations. |
| `tst/scripts/mhd/mhd_divb_amr.py` | Deep AMR regression behavior and acceptance thresholds. |

### Selected Conflict Strategy

Keep the feature stacked on
`feature/CR_tracers_followup_architecture@64a4d1be` while Phase 0 evidence and
small branch-local commits are produced.  Do not opportunistically rebase
before the baseline is archived.  Before an integration or merge-ready claim:

1. Refresh `origin/development` and rerun this audit.
2. Resolve conflicts in a dedicated integration commit or a clean successor
   branch with one reviewed claim.
3. Review all ten overlapping paths, not only the five content-conflict paths.
4. Rerun the fixed-energy particle suites, new relativistic ladder, deep AMR
   `divB` regression, restart checks, output checks, and documentation overlay.
5. Rerun affected checkpoints after conflict resolution.

## Current Code Baseline

### Particle State And Runtime

| Surface | Audited behavior | Relativistic obligation |
| --- | --- | --- |
| `src/athena.hpp:65` | Three integer indices (`PGID`, `PTAG`, `PSP`) and fourteen real indices (`IPX` through `IPDB`). | Version any new state meaning.  Do not scatter new index assumptions. |
| `src/particles/particles.cpp:200` | Runtime pusher selection supports `drift` and legacy `boris`; interpolation supports `tsc`, `lin_legacy`, and `trilinear`. | Add a new explicit opt-in mode.  Preserve legacy `boris`. |
| `src/particles/particles_pushers.cpp:212` | Legacy Boris gathers `B`, rotates stored velocity, and updates displacement. | Keep this fixed-energy control unchanged.  Add a separate reviewed relativistic path. |
| `src/particles/particles_pushers.cpp:221` | Positive `IPM` enters the magnetic rotation denominator; negative `IPM` aligns velocity with sampled `B`. | Resolve reciprocal, sign, sentinel, and migration semantics in `DR-005`. |
| `src/particles/particles_pushers.cpp:325` | Existing subcycling bounds cell crossing, block crossing, and magnetic gyro angle once before the substep loop. | Redesign for momentum-changing electric kicks and refresh the outer limit. |
| `src/particles/particles_tasks.cpp:31` | Particle push and exchange execute in `before_timeintegrator`. | Do not claim dynamic-MHD temporal order until scheduling is qualified separately. |
| `src/pgen/part_random.cpp:535` | `dtnew` is initialized once under fixed-speed and uniform-mesh assumptions. | Add a reviewed accelerated-particle timestep refresh policy. |

### Migration And Restart

| Surface | Audited behavior | Relativistic obligation |
| --- | --- | --- |
| `src/bvals/bvals_part.cpp:608` | MPI pack/unpack loops iterate over runtime `nidata` and `nrdata`, with integer and real buffers separated. | Appended state can migrate generically, but new layouts still need MPI and AMR tests. |
| `src/bvals/bvals_part.cpp:763` | Compaction copies runtime `nidata` and `nrdata`. | Verify no stale derived cache survives compaction. |
| `src/particles/particles.cpp:759` | AMR remap wraps particle coordinates periodically in device-table and host-tree paths. | Restrict the first relativistic mode to deliberate periodic semantics. |
| `src/particles/particles.cpp:147` | Restart sizing assumes `3 + 17*nparticle` `Real` values and derives the shard from the current MPI rank. | Replace fixed assumptions with an explicit new schema and selected rank-count policy. |
| `src/pgen/part_random.cpp:221` | Restart load decodes each particle as seventeen `Real` values and casts the first three back to integers. | Reject unsafe conversion and document each supported version. |
| `src/outputs/rst_prtcl.cpp:80` | Writer emits a three-`Real` header plus seventeen `Real` values per particle, casting integer identity fields through `Real`. | New relativistic records must serialize integers as integers. |
| `src/outputs/rst_prtcl.cpp:117` | Writer adds a `.pmeta` sidecar with magic, version, field counts, length, and checksum. | Preserve useful validation ideas, but do not mistake this sidecar for a complete extensible restart schema. |
| `scripts/particles/cr_tracer_inspect.py:19` | Inspector still uses `RESTART_FIELDS = 17`; metadata is optional and the binary payload remains `legacy_prst`. | Decode every retained version explicitly and reject malformed or unsafe records. |

### Existing Output Registration

`src/outputs/outputs.cpp:259` through `src/outputs/outputs.cpp:291` registers the
particle-facing outputs that must be classified in `DR-010`: `pvtk`, `trk`,
`df`, `dxh`, `drh`, `dparh`, `pmom`, `pspec`, `pspec2`, `psamp`, `ppd`, and
`prst`.

### Known Output Hazards

| Finding | Code evidence | Consequence |
| --- | --- | --- |
| Legacy spectrum `p` is speed, not momentum. | `src/outputs/df_prtcl.cpp:456` | Do not reuse the name for relativistic momentum silently. |
| Legacy spectrum `E` and `logE` derive from `0.5*v^2`. | `src/outputs/df_prtcl.cpp:458` | Add explicit relativistic kinetic-energy naming or versioned semantics. |
| Legacy joint spectra derive from speed, `0.5*v^2`, `v_parallel`, and `v_perp`. | `src/outputs/df_prtcl.cpp:674` | Preserve legacy interpretation or version the new mode's outputs. |
| `psamp p` aliases speed. | `src/outputs/df_prtcl.cpp:872` | Do not label it as momentum in the new mode. |
| `psamp energy` is `0.5*v^2`. | `src/outputs/df_prtcl.cpp:912` | Add an unambiguous relativistic field. |
| `psamp mass` emits overloaded `IPM`. | `src/outputs/df_prtcl.cpp:856` | Do not expose a misleading physical mass label after charge-semantics work. |
| Tracked output repeatedly writes `data[0]` and omits `sizeof(float)` in tag offsets. | `src/outputs/track_prtcl.cpp:140` | Repair separately or exclude `trk` from relativistic acceptance evidence. |
| Particle VTK loops over all integer fields but names only `PGID` and `PTAG`. | `src/outputs/vtk_prtcl.cpp:188` | Repair separately or exclude `pvtk`; `PSP` is structurally suspect. |

## Historical-Plan Reconciliation

The prior plans remain evidence.  They are not silently redefined as completed
relativistic behavior.

| Historical item | Source plan | Classification | Code evidence | Test evidence | Action |
| --- | --- | --- | --- | --- | --- |
| Preserve legacy `drift`, `boris + tsc`, and `boris + lin` behavior. | `cr_tracer_feature_branch_plan.md`, Design Constraints | implemented | `src/particles/particles.cpp:200`, `src/particles/particles_pushers.cpp:486` | Fresh Phase 0 rerun pending. | Retain as mandatory compatibility controls. |
| Device-table AMR remap with host-tree fallback. | `cr_tracer_followup_architecture_plan.md`, Phase 1 | implemented | `src/particles/particles.cpp:759` | Fresh AMR baseline pending. | Retain and qualify appended state. |
| Rank-count all-to-all metadata exchange with fallback. | `cr_tracer_followup_architecture_plan.md`, Phase 2 | implemented | `src/particles/particles.cpp:135`, `src/bvals/bvals_part.cpp:603` | Fresh MPI baseline pending. | Retain and qualify new state. |
| Device-side compaction and runtime-sized migration buffers. | `cr_tracer_followup_architecture_plan.md`, Phase 3 | implemented | `src/bvals/bvals_part.cpp:690`, `src/bvals/bvals_part.cpp:763` | Fresh empty-rank and migration baseline pending. | Retain; add stale-cache checks if derived state is stored. |
| Modular magnetic gather policies. | `cr_tracer_followup_architecture_plan.md`, Phase 4 | implemented | `src/particles/particles_pushers.cpp:48`, `src/particles/particles_pushers.cpp:88`, `src/particles/particles_pushers.cpp:142` | Existing fixed-energy accuracy report; fresh rerun pending. | Reuse policy structure.  Do not use `lin_legacy` for the first relativistic gather without a matching fluid gather decision. |
| Future `pcs` and `ct_aware` gathers. | `cr_tracer_accuracy_test_plan.md`, acceptance criteria | deferred | No implementation in current enum. | None. | Keep out of first relativistic feature. |
| Fixed-energy particle subcycling for cell, block, and gyro limits. | `cr_tracer_followup_architecture_plan.md`, Phase 5 | implemented | `src/particles/particles_pushers.cpp:325` | Fresh subcycle controls pending. | Retain for legacy mode; supersede with acceleration-aware bounds in the new mode. |
| Reduced CR diagnostics and selected samples. | `cr_tracer_followup_architecture_plan.md`, Phase 6 | implemented | `src/outputs/df_prtcl.cpp`, `src/outputs/dxhist_prtcl.cpp` | Inspector baseline pending. | Preserve legacy semantics and add explicit relativistic names. |
| Particle-aware AMR load balancing. | `cr_tracer_followup_architecture_plan.md`, Phase 7 | implemented | `src/mesh/mesh_refinement.cpp:67`, `src/mesh/mesh_refinement.cpp:467` | Fresh weighted-AMR baseline pending. | Retain and include in later AMR qualification. |
| Restart metadata, checksum, and inspector format reporting. | `cr_tracer_followup_architecture_plan.md`, Phase 8 | retained | `src/outputs/rst_prtcl.cpp:117`, `scripts/particles/cr_tracer_inspect.py:131` | Fresh round-trip baseline pending. | Preserve useful validation behavior, but replace the fixed legacy payload for the new mode. |
| Twelve-test fixed-energy accuracy ladder. | `cr_tracer_accuracy_test_plan.md`; `cr_tracer_accuracy.md` | implemented | `tst/test_suite/particles/test_particles_cr_accuracy_cpu.py`, `tst/test_suite/particles/test_particles_cr_accuracy_mpicpu.py` | Historical report in `cr_tracer_accuracy.md`; fresh Phase 0 rerun pending. | Retain unchanged and add a separate relativistic ladder. |
| GPU qualification on the same accepted SHA. | `cr_tracer_gpu_testing_handoff.md` | deferred | No local accelerator evidence. | Not available on this workstation. | Keep as separate unqualified residual risk.  It is optional for workstation merge-ready, not silently passed. |
| Particle feedback, current deposition, Hall terms, and coupled PIC behavior. | Relativistic implementation guide scope exclusions | not applicable | Current feature baseline is passive. | None required for first feature. | Keep out of scope. |

## Normalized Equations

### Physical-Units Baseline

The initial production feature is passive Newtonian ideal-MHD sampling.  It
must not take a curl inside the particle pusher.  The physical baseline is:

```text
cE_particle = -u_fluid,particle x B_particle

w = p / m = gamma v
gamma = sqrt(1 + |w|^2 / c^2)
v = w / gamma

dw/dt = (q / (m c)) [cE_particle + v x B_particle]

K = (gamma - 1) m c^2
dK/dt = (q / c) cE_particle dot v
```

The curl belongs to the MHD induction update:

```text
dB/dt = -curl(E)
```

### Code-Unit Contract Still To Derive

The implementation must not copy physical symbols mechanically into code.
Before a production edit, `DR-001` must derive:

```text
cE_code = -u_fluid,code x B_code
dw_code/dt = alpha_code [cE_code + v_code(w_code) x B_code]
gamma = sqrt(1 + |w_code|^2 / C_state^2)
v_code = w_code / gamma
dK_model/dt = beta_code cE_code dot v_code
```

`DR-001` must state how `alpha_code`, `beta_code`, `K_model`, and `C_state`
map to physical units and existing AthenaK MHD units.  If configurable reduced
light speed `C` replaces physical `c`, record it as a model approximation and
derive where it enters the state relation, force prefactor, field mapping, and
work identity.

### Existing Legacy Magnetic-Only Normalization

For positive `IPM`, the current magnetic-only Boris rotation is consistent
with the normalized legacy relation:

```text
dv/dt = (1 / IPM) v x B
```

That observation does not settle the new mode's physical charge convention.
Current negative `IPM` values are sentinel-like: they align stored velocity
with sampled `B`.  `IPM == 0` would divide by zero in the legacy rotation.
`DR-005` must resolve migration, charge sign, reciprocal safety, and sentinel
behavior explicitly.

## Compatibility Matrix

This matrix is the selected initial fail-fast boundary for the new opt-in
relativistic mode.  It does not restrict existing legacy modes.

| Combination | Initial relativistic-mode disposition | Evidence required to widen |
| --- | --- | --- |
| Passive one-way MHD-to-particle coupling | Selected scope | Keep deposition and feedback absent in source review. |
| Newtonian ideal MHD | Selected production model after qualification | Gather, pointwise ideal-field, prescribed-versus-coupled, and solver-coupling tests. |
| Prescribed analytical fields | Selected test-harness-only model | Mechanically isolate from production ideal-MHD construction. |
| `time/evolution = static` or explicitly frozen background | Selected first qualification target | Kernel and coupled frozen-background convergence. |
| Kinematic evolution | Initially reject | Separate temporal-centering qualification. |
| Dynamic MHD evolution | Initially reject for production claim | `CP-4 Solver Coupling`, including time-dependent manufactured convergence. |
| Three spatial dimensions with 3-vector momentum and field | Selected first supported dimensional contract | Dedicated 3-D analytical and infrastructure ladder. |
| Two spatial dimensions with 3-vector momentum and field (`2D3V`) | Initially reject | Separate explicit contract and validation. |
| One spatial dimension | Reject | Current particle module already rejects 1-D. |
| Strictly periodic particle coordinates | Selected first boundary contract | Periodic wrap and migration tests. |
| Reflecting, inflow, outflow, user, or shear-periodic particle boundaries | Reject | Separate boundary design and review. |
| Active turbulence forcing | Initially reject | Separate evolving-background temporal qualification. |
| Frozen turbulent snapshot | Eligible after static-background qualification | Frozen-field comparison and convergence. |
| Ion-neutral coupling | Reject | Separate model derivation and scheduling review. |
| Resistivity or non-ideal electric fields | Reject | Separate physical model and field-source decision. |
| Orbital advection or shearing box | Reject | Separate frame and boundary derivation. |
| SRMHD | Reject | Separate relativistic fluid-frame field transformation design. |
| GRMHD or dynamical GR | Reject | Separate covariant model design. |
| Particle feedback, current deposition, Hall corrections | Reject | Separate coupled-feature branch. |
| GPU backend | Optional for merge-ready on this workstation; unqualified residual risk | Run accelerator build, analytical, MPI, AMR, restart, and diagnostic gates on the accepted SHA before claiming `GPU QUALIFIED`. |

## Output Audit Classification

| Output | Legacy meaning | Initial relativistic-mode classification | Required action |
| --- | --- | --- | --- |
| `ppd` | Species and position dump. | Preserve. | Validate unchanged positions and species semantics. |
| `prst` | Legacy three-`Real` header plus seventeen-`Real` particle records; optional `.pmeta` sidecar. | Version. | Define a new schema with typed integers, metadata, checksum, shard manifest, and explicit legacy policy. |
| `df` | Pitch-angle cosine histogram from stored velocity and sampled `B`. | Preserve legacy; extend deliberately. | Compute from documented derived physical velocity in the new mode and validate semantics. |
| `dxh` | Component displacement histograms. | Preserve. | Validate that physical displacement remains unchanged in meaning. |
| `drh` | Scalar displacement histogram. | Preserve. | Validate that physical displacement remains unchanged in meaning. |
| `dparh` | Parallel displacement histogram using accumulated `IPDB`. | Extend or version. | Define how changing sampled `B` and acceleration affect interpretation. |
| `pmom` | Reduced moments including velocity and speed-squared moments. | Version or add explicit relativistic sibling fields. | Do not silently reinterpret legacy columns as momentum or relativistic energy. |
| `pspec` | `p` means speed; `E` means `0.5*v^2`; `logE` follows that proxy. | Version or reject ambiguous quantities in the new mode. | Add explicit momentum and relativistic kinetic-energy quantity names. |
| `pspec2` | Joint histograms based on `mu`, speed, `0.5*v^2`, `v_parallel`, and `v_perp`. | Version or add explicit relativistic quantities. | Preserve legacy quantity names only with legacy meanings. |
| `psamp` | Selected raw and derived fields; `p` aliases speed, `energy` is `0.5*v^2`, and `mass` emits `IPM`. | Extend and repair naming. | Add explicit momentum, `gamma`, relativistic energy, sampled `cE`, work, and charge-control names.  Do not retain misleading aliases for new semantics. |
| `trk` | Position and velocity record for selected tags. | Intentionally unsupported as acceptance evidence until repaired. | Repair repeated-buffer and byte-offset defects in a separately reviewed prerequisite or isolate the output. |
| `pvtk` | Legacy particle positions and integer scalars. | Intentionally unsupported as acceptance evidence until repaired. | Repair integer-field naming and audit new-state coverage separately. |
| Inspector | Decodes legacy restart, `ppd`, histograms, moments, and samples. | Extend. | Decode every retained restart and output version with explicit names and failure checks. |

## Live Blocking Findings

| ID | Severity | Finding | Blocking scope | Required closure evidence |
| --- | --- | --- | --- | --- |
| `BF-001` | blocking | The tree was clean at the first status check, then untracked Python caches, `tst/test_log.txt`, and concurrent standalone pusher-spike artifacts appeared during the read-only audit and verification. | Any accepted baseline claim. | Attribute or remove the artifacts deliberately outside this ledger-only edit, then rerun and archive `git status --porcelain=v1 --untracked-files=all`. |
| `BF-002` | blocking | Fresh fixed-energy CPU, MPI, accuracy, style, deep AMR, restart, and inspector baseline evidence has not been run and archived on this branch. | `CP-0 Baseline`. | Complete the baseline evidence placeholders below on a clean agreed tree. |
| `BF-003` | blocking | Independent physics, integration, and test-adversary reviewers have not reviewed Phase 0. | `CP-0`, `CP-1`, and `RG-001`. | Archive prompts, findings, dispositions, and rechecks. |
| `BF-004` | blocking | The exact code-unit normalization, physical-versus-reduced light-speed model, and `IPM` migration semantics remain unresolved. | State-helper and pusher edits. | Accept `DR-001` and `DR-005` after derivation and review. |
| `BF-005` | blocking | A concurrent standalone Boris-versus-Higuera-Cary candidate spike appeared after the first ledger draft, but no preregistered, audited, independently reviewed pusher-selection artifact has been accepted. | Pusher selection and any pusher edit. | Audit the candidate artifact, preregister acceptance criteria before using results for selection, justify included and excluded methods, archive metrics and independent references, and accept `DR-003`. |
| `BF-006` | blocking | The restart reader and writer retain fixed seventeen-`Real` records, cast integer identity through `Real`, and derive restart shard selection from the current MPI rank. | Restart edit set and merge-ready claim. | Accept and implement `DR-009`; test typed integer safety, sidecars, shard policy, corruption, continuation, and rank-count behavior. |
| `BF-007` | blocking | Solver coupling currently pushes particles before time-integrator stages.  A static-field midpoint gather does not establish dynamic-MHD temporal order. | Dynamic-MHD support claim and `CP-4 Solver Coupling`. | Accept `DR-007`; archive time-dependent manufactured convergence and scheduling review. |
| `BF-008` | blocking | `trk` and `pvtk` have confirmed pre-existing structural defects. | Use of these formats as acceptance evidence. | Repair in isolated reviewed work or explicitly keep unsupported for the first relativistic mode. |

## Residual Risks That Are Not Phase 0 Blockers

| ID | Risk | Current disposition |
| --- | --- | --- |
| `RR-001` | GPU behavior is unqualified because this workstation supplies CPU and MPI execution only. | GPU qualification is optional for workstation merge-ready but remains a separate unqualified residual risk. |
| `RR-002` | Pointwise `cE = -u x B` gathered from cell-centered quantities is not the CT edge electromotive force. | Use the pointwise construction first; document the tradeoff and require convergence evidence. |
| `RR-003` | A full-orbit method may be unresolved or too costly when `r_L / dx` or gyro-angle resolution is poor. | Add applicability diagnostics and a revisit trigger; do not hide a guiding-center model inside this branch. |
| `RR-004` | Five merge-tree content conflicts and five additional overlaps create integration risk. | Follow the selected `DR-000` strategy and rerun dependent gates after reconciliation. |

## Decision Records

## DR-000: Base, Target, Stacking, And Backend Qualification Strategy

- Status: selected
- Question: What base, eventual integration target, stacking strategy, and required backend verdicts govern the branch?
- Selected option: freeze `feature/CR_tracers_followup_architecture@64a4d1be`; keep the current branch stacked while Phase 0 is archived; target `origin/development@c6a73b08`; resolve integration overlaps in a dedicated reviewed step; refresh merge-tree after target changes and before handoff.
- Backend qualification: `CPU/MPI QUALIFIED` is required.  `GPU QUALIFIED` is optional for `MERGE READY` on this workstation, but unavailable GPU evidence remains an explicit unqualified residual risk.
- Evidence and reasoning: current branch ancestry, exact merge-tree audit, existing GPU handoff, and workstation limits.
- Compatibility and migration impact: all dependent checkpoints must rerun after rebase or conflict resolution.
- Failure signals: changed base or target without a ledger update; unreviewed conflict resolution; accelerator claim without evidence.
- Revisit trigger: target SHA changes, intended deployment backend changes, or a reviewer requires earlier integration.

## DR-001: Physical Normalization And Reduced-Light Model

- Status: proposed
- Question: What exact code-unit equations map physical `cE`, `q/m` or `m/q`, and physical `c` or reduced `C` into the pusher and work diagnostic?
- Proposed direction: adopt the physical-units baseline above, then derive the AthenaK code-unit map explicitly.  Do not replace `c` with configurable `C` mechanically.
- Evidence and reasoning: the guide requires explicit factors; current legacy code only establishes the magnetic-only `1/IPM` coefficient.
- Compatibility and migration impact: controls state conversion, force prefactor, initialization, diagnostics, and restart metadata.
- Failure signals: implicit units, mixed `E` and `cE`, or a work identity not derived from the accepted model.
- Revisit trigger: completed dimensional derivation and physics review.

## DR-002: Authoritative Particle State

- Status: selected
- Question: Which evolving particle quantity is authoritative?
- Selected option: use mass-normalized momentum `w = p/m = gamma*v` as the authoritative relativistic state; derive `gamma` and physical velocity.  Treat cached derived fields as proposed until synchronization cost and invariants are reviewed.
- Evidence and reasoning: directional momentum is required by the force update; storing only `gamma` is insufficient; derived velocity avoids making two mutable representations authoritative.
- Compatibility and migration impact: exact appended-versus-replaced layout and legacy conversion remain proposed pending `DR-009` and `DR-010`.
- Failure signals: stale velocity cache, ambiguous restart authority, or a state helper that reconstructs `|v| >= C_state`.
- Revisit trigger: helper spike shows a compelling compatibility or performance reason to cache a derived quantity.

## DR-003: Relativistic Electromagnetic Pusher

- Status: proposed
- Question: Which reviewed second-order pusher should be integrated?
- Proposed options: standard relativistic Boris, Vay, and Higuera-Cary, with any exclusion justified from primary references and spike evidence.
- Selected option: none.  Pusher selection remains proposed pending the preregistered comparison spike.
- Evidence and reasoning: familiarity is not sufficient evidence; crossed-field drift and boosted-frame cancellation are mandatory discriminators.
- Compatibility and migration impact: blocks pusher algebra, position split, work quadrature, and detailed subcycling decisions.
- Failure signals: selection before spike review, copied formulas without primary-reference audit, or missing independent oracle.
- Revisit trigger: spike metrics and independent physics review are archived.

## DR-004: Explicit Runtime Compatibility Mode

- Status: selected
- Question: Should relativistic acceleration replace legacy `boris` silently?
- Selected option: add a new explicit opt-in runtime mode and preserve legacy `boris` as the fixed-energy compatibility control.
- Evidence and reasoning: historical plans require existing `drift`, `boris + tsc`, and `boris + lin` behavior to remain available.
- Compatibility and migration impact: existing decks must run unchanged.
- Failure signals: old input activates new semantics or old output meanings change silently.
- Revisit trigger: none for the first feature.

## DR-005: Charge, Species, And `IPM` Semantics

- Status: proposed
- Question: What exactly does `IPM` mean, how are charge signs represented, and what happens to the legacy negative sentinel?
- Proposed direction: introduce explicit new-mode charge semantics without silently reinterpreting legacy `IPM`; retain legacy sentinel behavior only inside legacy mode unless an intentional migration is accepted.
- Evidence and reasoning: current positive `IPM` is used as a denominator; current negative `IPM` aligns velocity with `B`.
- Compatibility and migration impact: affects initialization, state helpers, restart conversion, spectra, samples, and failure injection.
- Failure signals: unnoticed reciprocal or sign change, `IPM == 0`, or overloaded sentinel values in the new mode.
- Revisit trigger: units derivation and migration matrix are reviewed.

## DR-006: Spatial Gather And Ideal-Field Construction

- Status: selected
- Question: How should the first production mode sample fields?
- Selected option: gather fluid velocity and magnetic field separately with the same documented non-legacy policy, then construct pointwise `cE_particle = -u_fluid,particle x B_particle`.  Reject `lin_legacy` initially because it lacks a matching clean fluid-velocity gather.
- Proposed detail: choose `trilinear`, `tsc`, or a reviewed selectable subset after gather-only comparisons.
- Evidence and reasoning: current gather-policy structure is reusable; CT edge EMF staggering and time availability differ.
- Failure signals: mixed gather policies, wrong sign, unexpected `E dot B`, or silent use of the CT edge field.
- Revisit trigger: gather-only convergence and AMR-boundary evidence.

## DR-007: Temporal Centering And Solver Coupling

- Status: selected
- Question: What temporal claim is allowed before stage-coupled evidence exists?
- Selected option: qualify a frozen-background reference path first, with spatial midpoint gather from a clearly stated grid state.  Do not claim second-order dynamic-MHD coupling from static-field evidence.
- Evidence and reasoning: current particle tasks run in `before_timeintegrator`.
- Compatibility and migration impact: static or frozen backgrounds can be qualified before dynamic MHD; evolving backgrounds remain rejected initially.
- Failure signals: dynamic-MHD claim without time-dependent manufactured convergence.
- Revisit trigger: `CP-4 Solver Coupling` evidence.

## DR-008: Dimensional Semantics

- Status: selected
- Question: What spatial dimensionality is supported first?
- Selected option: restrict the first relativistic mode to 3-D.  Keep `2D3V` proposed as a later explicitly validated widening.
- Evidence and reasoning: current particles support 2-D and 3-D, but the relativistic field, boundary, and diagnostic contract should not inherit lower-dimensional assumptions silently.
- Compatibility and migration impact: legacy 2-D behavior remains unchanged.
- Failure signals: new mode runs in 2-D without an accepted widening decision.
- Revisit trigger: dedicated `2D3V` analytical and migration tests.

## DR-009: Restart Schema And Rank-Count Policy

- Status: proposed
- Question: What typed, versioned schema is written; which legacy records load; and are MPI rank-count changes supported or rejected?
- Proposed direction: add an extensible schema with typed integers, explicit field dictionary, magic, version, checksum, checkpoint time and cycle, topology policy, and shard manifest.  Choose either manifest-based old-shard discovery plus redistribution or deterministic rank-count-change rejection.
- Evidence and reasoning: current payload is fixed seventeen-`Real`; current readers derive shard choice from the current MPI rank; the `.pmeta` sidecar is useful but incomplete.
- Compatibility and migration impact: legacy conversion is allowed only where precision and semantics are safe.
- Failure signals: magic record sizes, lossy identity conversion, stale sidecar acceptance, or undocumented rank-count behavior.
- Revisit trigger: restart design review and corruption-test plan.

## DR-010: Diagnostic Migration And Naming

- Status: selected
- Question: How are legacy outputs preserved without ambiguous relativistic reuse?
- Selected option: retain legacy names with legacy meanings; add explicit relativistic momentum, `gamma`, kinetic-energy, sampled-`cE`, work, and applicability names; version or reject ambiguous new-mode requests.
- Evidence and reasoning: current `pspec p`, `psamp p`, legacy energy, and `psamp mass` names are not valid relativistic definitions.
- Compatibility and migration impact: use the output audit classification above.  Isolate `trk` and `pvtk` until repaired.
- Failure signals: silent column reinterpretation or analysis scripts unable to identify units.
- Revisit trigger: output dictionary and restart schema review.

## DR-011: Acceleration-Aware Subcycling

- Status: proposed
- Question: Which limits protect changing-momentum trajectories?
- Proposed direction: bound cell crossing, block crossing, relativistic gyro angle, electric kick, fractional momentum change, optional fractional kinetic-energy change, and total subcycle count.  State whether constraints are recomputed after kicks or conservatively bounded.
- Evidence and reasoning: current loop computes fixed-energy bounds once before subcycling.
- Compatibility and migration impact: retain legacy policy for legacy mode.
- Failure signals: accelerated particle violates ownership assumptions or cap handling hangs.
- Revisit trigger: pusher spike and subcycle-limit unit tests.

## DR-012: Failure Policy And Validity Domain

- Status: selected
- Question: Which inputs and states fail fast?
- Selected option: reject nonfinite state, invalid model speed, unsupported dimensions and mode combinations, `IPM == 0`, reciprocal overflow or underflow, superluminal reconstructed velocity, unsupported charge or sentinel cases, and sampled states outside the accepted model domain.
- Proposed detail: preregister exact warning-versus-error thresholds for electric-dominated and near-limit states.
- Evidence and reasoning: silent invalid state can produce plausible trajectories.
- Failure signals: NaN propagation, silent clipping, or unsupported combination execution.
- Revisit trigger: accepted `DR-001`, `DR-005`, and failure-injection tests.

## DR-013: Electric-Field Scope

- Status: selected
- Question: Which electric-field sources belong in the first production feature?
- Selected option: production support is passive Newtonian ideal MHD only, using sampled pointwise `cE = -u x B`.  Prescribed fields are test-harness-only.  Non-ideal, resistive, feedback, Hall, SRMHD, and GRMHD fields remain out of scope.
- Evidence and reasoning: bounded passive scope is the guide's initial scientific contract.
- Failure signals: production path acquires feedback or non-ideal behavior without a reopened scope decision.
- Revisit trigger: separate follow-up branch.

## DR-014: Position Integration Split

- Status: proposed
- Question: Which position update is consistent with the selected pusher and gather timing?
- Proposed options: drift-kick-drift, kick-drift-kick, or another reviewed second-order split.
- Selected option: none pending `DR-003` and `DR-007`.
- Evidence and reasoning: the existing magnetic path uses a midpoint spatial gather, but that does not settle the relativistic split.
- Failure signals: first-order trajectory convergence or inconsistent gather timing.
- Revisit trigger: pusher spike and analytical position-oracle results.

## DR-015: Explicit Work Accumulation

- Status: proposed
- Question: Which pusher-consistent explicit quadrature closes against kinetic-energy change?
- Proposed direction: accumulate sampled-field work directly from the accepted update, never from `Delta K`; preregister closure tolerances before inspecting results; decide whether accumulated work is serialized.
- Evidence and reasoning: work closure is required for acceleration and deceleration validation.
- Compatibility and migration impact: affects diagnostics, restart, and analytical tests.
- Failure signals: tautological diagnostic or post-hoc tolerance.
- Revisit trigger: accepted pusher algebra and units derivation.

## DR-016: Outer Timestep Refresh

- Status: proposed
- Question: Where is the mesh-level particle timestep constraint recomputed after momentum changes?
- Proposed direction: add a reviewed refresh point so the next outer step consumes current particle limits.
- Evidence and reasoning: `part_random.cpp` initializes `dtnew` once using fixed-speed assumptions.
- Compatibility and migration impact: retain legacy timestep behavior unless explicitly changed.
- Failure signals: next outer step consumes stale pre-acceleration limits.
- Revisit trigger: task-order audit and accelerated timestep tests.

## DR-017: Particle Boundary Semantics

- Status: selected
- Question: Which coordinate boundaries does the first relativistic mode support?
- Selected option: `strictly_periodic` only, with fail-fast rejection otherwise.
- Evidence and reasoning: current AMR remap wraps coordinates periodically; nonperiodic relativistic behavior is not implemented deliberately.
- Compatibility and migration impact: legacy modes remain unchanged.
- Failure signals: silent periodic wrapping under a requested nonperiodic boundary.
- Revisit trigger: separately reviewed boundary implementation.

## DR-018: Compatibility Matrix

- Status: selected
- Question: Which physics-module and evolution combinations are initially allowed?
- Selected option: enforce the compatibility matrix above for the new opt-in mode.  Start with 3-D, strictly periodic, passive Newtonian ideal-MHD static or explicitly frozen backgrounds; reject every unqualified widening.
- Evidence and reasoning: the current code contains multiple solver modules and scheduling paths whose semantics cannot be inherited silently.
- Compatibility and migration impact: legacy modes remain available.
- Failure signals: unlisted combination executes or unsupported mode is only documented rather than rejected.
- Revisit trigger: separate qualification evidence for a named widening.

## Baseline Evidence Placeholders

Historical fixed-energy evidence exists in `cr_tracer_accuracy.md`: the
2026-05-23 report records `12` CPU accuracy passes, `5` MPI CPU accuracy passes,
style passes, and a deep AMR `divB` pass.  Those results are inherited context,
not a substitute for fresh Phase 0 artifacts on an agreed clean tree.

## Evidence: Branch And Dirty-State Freeze

- Commit: `3984b57ec8832b3e3e8140e08188a8a05e773136`
- Runtime command: `git status --porcelain=v1 --untracked-files=all`
- Expected result: agreed clean baseline, or an explicitly attributed and accepted artifact list.
- Measured result: `HOLD`; unexpected untracked files appeared during audit.
- Acceptance threshold: no unexplained dirty-tree entry.
- Pass or fail: fail
- Artifact paths: pending archive
- Notes: do not remove another worker's artifacts without attribution.

## Evidence: Exact Merge-Tree Audit

- Commit: `3984b57ec8832b3e3e8140e08188a8a05e773136`
- Runtime command: `git merge-tree --write-tree HEAD c6a73b08e60807f8b925164c5e7edd5cb820c8ae`
- Expected result: enumerate all conflicts and same-ancestry overlaps.
- Measured result: five content conflicts and five additional overlaps, recorded above.
- Acceptance threshold: accepted `DR-000` conflict strategy and complete affected-path inventory.
- Pass or fail: pass for inventory; integration remains pending
- Artifact paths: pending command-log archive
- Notes: rerun after target updates and before handoff.

## Evidence: Legacy Particle CPU Suite

- Commit: pending clean baseline SHA
- Build configuration: pending compiler, Kokkos backend, and build command
- Runtime command: `cd tst && python run_test_suite.py --test test_suite/particles/test_particles_cr_cpu.py --cpu`
- Expected result: all existing CPU particle regressions pass unchanged.
- Measured result: concurrent unreviewed candidate files appeared during ledger verification; no accepted spike result yet
- Acceptance threshold: no failure and no relaxed assertion.
- Pass or fail: pending
- Artifact paths: pending
- Notes: archive stdout, stderr, executable SHA, and configuration.

## Evidence: Legacy Particle MPI CPU Suite

- Commit: pending clean baseline SHA
- Build configuration: pending compiler, MPI version, Kokkos backend, rank counts, and build command
- Runtime command: `cd tst && python run_test_suite.py --test test_suite/particles/test_particles_cr_mpicpu.py --mpicpu`
- Expected result: all existing MPI particle regressions pass unchanged.
- Measured result: pending
- Acceptance threshold: no failure; include empty-rank and migration coverage.
- Pass or fail: pending
- Artifact paths: pending
- Notes: record supported rank counts.

## Evidence: Fixed-Energy Accuracy CPU Suite

- Commit: pending clean baseline SHA
- Build configuration: pending
- Runtime command: `cd tst && python run_test_suite.py --test test_suite/particles/test_particles_cr_accuracy_cpu.py --cpu`
- Expected result: all twelve fixed-energy CPU accuracy controls pass unchanged.
- Measured result: pending
- Acceptance threshold: existing preregistered suite thresholds.
- Pass or fail: pending
- Artifact paths: pending
- Notes: preserve the historical report as comparison context.

## Evidence: Fixed-Energy Accuracy MPI CPU Suite

- Commit: pending clean baseline SHA
- Build configuration: pending
- Runtime command: `cd tst && python run_test_suite.py --test test_suite/particles/test_particles_cr_accuracy_mpicpu.py --mpicpu`
- Expected result: all fixed-energy MPI accuracy controls pass unchanged.
- Measured result: pending
- Acceptance threshold: existing preregistered suite thresholds.
- Pass or fail: pending
- Artifact paths: pending
- Notes: archive decomposition and empty-rank evidence.

## Evidence: Style And Diff Hygiene

- Commit: pending clean baseline SHA
- Build configuration: not applicable
- Runtime command: `cd tst && python run_test_suite.py --style`; then `git diff --check`
- Expected result: style and whitespace checks pass.
- Measured result: pending for style; pre-ledger `git diff --check` passed.
- Acceptance threshold: no style or whitespace failure.
- Pass or fail: pending
- Artifact paths: pending
- Notes: rerun after the ledger draft.

## Evidence: Deep AMR `divB` Compatibility

- Commit: pending clean baseline SHA
- Build configuration: pending
- Runtime command: `cd tst && python run_tests.py mhd/mhd_divb_amr`
- Expected result: existing 1-D, 2-D, and 3-D normalized-divergence regression passes.
- Measured result: pending
- Acceptance threshold: existing thresholds; do not relax them.
- Pass or fail: pending
- Artifact paths: pending
- Notes: this is required again after development conflict resolution.

## Evidence: Restart And Inspector Baseline

- Commit: pending clean baseline SHA
- Build configuration: pending CPU and MPI configurations
- Runtime command: pending exact existing restart-producing particle cases plus `python scripts/particles/cr_tracer_inspect.py ...`
- Expected result: legacy round trip and inspector checks pass; sidecar behavior is recorded.
- Measured result: pending
- Acceptance threshold: exact counts, tags, species, record length, and checksum behavior.
- Pass or fail: pending
- Artifact paths: pending
- Notes: distinguish C++ load behavior from inspector-only metadata validation.

## Evidence: Relativistic Pusher Comparison Spike

- Commit: pending isolated spike artifact
- Build configuration: pending
- Runtime command: pending
- Expected result: compare at least two credible pushers against independent references for uniform `B`, crossed fields, force cancellation, high `gamma`, convergence, reversibility where expected, long-time stress, boosted-frame cancellation, finite-state limits, and device-suitable algebra.
- Measured result: pending
- Acceptance threshold: preregister before results inspection.
- Pass or fail: pending
- Artifact paths: pending
- Notes: the concurrent candidate compares relativistic Boris and Higuera-Cary in a standalone script.  Audit it before reuse; pusher selection remains proposed.

## Evidence: Normalized-Equation Review

- Commit: ledger SHA pending
- Build configuration: not applicable
- Runtime command: independent physics-review prompt from the guide
- Expected result: reviewer derives factors independently and checks `E` versus `cE`, `q/m` versus `m/q`, physical `c` versus reduced `C`, work identity, and legacy `IPM`.
- Measured result: pending
- Acceptance threshold: no unresolved algebra or units finding.
- Pass or fail: pending
- Artifact paths: pending
- Notes: required before accepting `DR-001`, `DR-003`, `DR-005`, `DR-014`, or `DR-015`.

## RG-001: Is The Feature Still Passive And Bounded?

- Date: `2026-05-29`
- Commit range reviewed: `64a4d1be..3984b57e`
- Evidence reviewed: guide, source inventory, historical plans, merge-tree audit, output audit, and dirty-tree status.
- Assumptions that still hold: the first feature can remain passive; an explicit opt-in runtime mode can preserve legacy fixed-energy behavior; pointwise ideal-MHD field construction is a bounded first production model.
- Assumptions weakened or falsified: the tree did not remain clean during audit; the refreshed overlap inventory adds `tst/scripts/mhd/mhd_divb_amr.py`; existing restart metadata is not a complete extensible schema.
- New risks: dynamic-MHD temporal order, restart identity safety, output ambiguity, integration conflict resolution, and unavailable GPU qualification.
- Blocking findings: `BF-001` through `BF-008`.
- Decision: stop
- Plan changes: keep Phase 0 open.  Do not open Phase 1 state-helper edits.
- Tests added or changed: none; ledger-only draft.
- Independent reviewer findings and responses: pending.
- Next smallest reviewable increment: attribute the dirty-tree artifacts, archive fresh baseline evidence, run the pusher spike, complete normalized-equation review, and update the proposed decision records.

## CP-0 Baseline

- Commit: `3984b57ec8832b3e3e8140e08188a8a05e773136` plus this uncommitted ledger draft
- Allowed edit set reviewed: guide, ledger, and read-only source inventory only
- Evidence reviewed: branch facts, exact merge-tree audit, historical-plan reconciliation, affected-file inventory, output audit, and inherited fixed-energy report.
- Reviewer verdicts: pending physics, integration, and test-adversary review.
- Assumptions falsified: clean tree remained stable during audit.
- Scope drift found: none in this ledger draft.
- Unresolved risks: code-unit normalization, `IPM`, pusher choice, restart schema, timestep refresh, work quadrature, output schema, integration conflicts, and GPU residual risk.
- Blocking findings: `BF-001` through `BF-008`.
- Evidence added: merge-tree inventory and baseline placeholders.
- Choices made since previous checkpoint, including retained defaults: selected records `DR-000`, `DR-002`, `DR-004`, `DR-006`, `DR-007`, `DR-008`, `DR-010`, `DR-012`, `DR-013`, `DR-017`, and `DR-018`.
- Decision records still proposed and why they do not authorize current edits: `DR-001`, `DR-003`, `DR-005`, `DR-009`, `DR-011`, `DR-014`, `DR-015`, and `DR-016` require derivation, spike evidence, or migration review.
- Verdict: HOLD
- Next allowed edit set: ledger updates and archived Phase 0 evidence only.

## CP-1 Physics

- Commit: pending
- Allowed edit set reviewed: no state-layout or pusher edit is open.
- Evidence reviewed: physical-units baseline only.
- Reviewer verdicts: pending independent physics and numerics review.
- Assumptions falsified: none assessed yet.
- Scope drift found: none.
- Unresolved risks: code-unit map, charge semantics, reduced-light model, pusher selection, position split, and work quadrature.
- Blocking findings: `BF-003`, `BF-004`, and `BF-005`.
- Evidence added: normalized-equation scaffold and pusher-spike placeholder.
- Choices made since previous checkpoint, including retained defaults: momentum-like authoritative state selected provisionally in `DR-002`.
- Decision records still proposed and why they do not authorize current edits: `DR-001`, `DR-003`, `DR-005`, `DR-014`, and `DR-015`.
- Verdict: HOLD
- Next allowed edit set: equation derivation, isolated spike artifacts, reviewer records, and ledger updates only.

## Unresolved Choices Summary

The following choices remain intentionally unresolved:

1. AthenaK code-unit mapping, physical `c` versus reduced `C`, and exact work identity.
2. Relativistic pusher selection after the preregistered spike.
3. New-mode charge sign, species, reciprocal, and legacy negative-`IPM` sentinel policy.
4. Exact appended-versus-replaced particle storage layout and whether any derived value is cached.
5. Typed restart schema, safe legacy conversion boundary, paired-checkpoint policy, shard manifest, and MPI rank-count-change policy.
6. Exact accelerated subcycle bounds and whether they are recomputed after kicks.
7. Position integration split.
8. Explicit sampled-field work quadrature and serialization policy.
9. Outer timestep refresh ownership and task position.
10. Final relativistic output dictionary and whether `trk` and `pvtk` are repaired or remain intentionally unsupported.
11. Exact first production gather subset: `trilinear`, `tsc`, or both after gather-only evidence.

## Phase 0 Review Update 1

This section is append-only evidence.  It supersedes the initial draft where
stated, but preserves the earlier `HOLD` record so later reviewers can see
which assumptions changed.

### Attributed Working-Tree Snapshot

After the baseline particle suites, style rerun, and deliberate cache cleanup,
the only expected untracked files are the three owned Phase 0 artifacts:

```text
docs/source/modules/cr_tracer_relativistic_acceleration_ledger.md
docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
scripts/particles/cr_relativistic_pusher_spike.py
```

The script and JSON are attributed to the isolated Phase 0 pusher experiment.
The ledger is attributed to the Phase 0 architecture audit.  Python cache
directories and `tst/test_log.txt` are generated test noise and must be removed
before each accepted dirty-state snapshot.

Snapshot checksums before the spike-v2 expansion:

```text
8473999dbfba444b214c1cc610c6a0f98790f0101da68ff84ad60a34186b3eef  scripts/particles/cr_relativistic_pusher_spike.py
b5c863b7b813271f6cb46b292463c03fc2f08e45996997f158a708493e3f9ed5  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
```

`BF-001` is closed for attribution.  Refresh the snapshot and checksums after
the expanded spike is regenerated.

### Independent Phase 0 Reviews

Two fresh reviewers independently inspected the guide, ledger, current source
surfaces, and isolated pusher experiment:

| Review | Disposition | Required response |
| --- | --- | --- |
| Physics and numerics | `HOLD` before production edits | Accept the normalized contract below; keep Higuera-Cary provisional until spike-v2 preregisters new criteria and closes missing adversarial cases; use a robust scaled or `hypot`-style root evaluation in production. |
| Architecture, integration, and validation | `HOLD` before production edits | Archive fresh baseline evidence; expand the positional-consumer inventory; repair the legacy particle-restart timestep semantic defect; refresh the relativistic particle timestep after kicks and after AMR. |

The review prompts and full findings remain in the implementation thread.  The
required changes are normalized below instead of being silently folded into
the initial draft.

### Additional Blocking Findings

| ID | Severity | Finding | Blocking scope | Required closure evidence |
| --- | --- | --- | --- | --- |
| `BF-009` | blocking | The legacy particle-restart writer stores `pm->dt` in header slot `1`, while the C++ loader restores slot `1` into `ppart->dtnew`; `Mesh::NewTimeStep()` later consumes that value as a particle timestep bound. | Restart qualification and merge-ready claim. | Repair the semantic mismatch deliberately and add a regression in which mesh `dt` differs from particle `dtnew`. |
| `BF-010` | blocking | The initial Higuera-Cary comparison was exploratory: results were inspected before expanded acceptance criteria were preregistered, Vay was not implemented or explicitly excluded, and boosted-cancellation, phase-volume, resonance, pure-electric, and explicit-work cases were incomplete. | `DR-003`, pusher edits, and `CP-1 Physics`. | Regenerate a spike-v2 artifact whose new-case criteria are emitted before execution; either implement Vay from a primary source or archive a bounded literature-backed exclusion; rerun independent review. |
| `BF-011` | blocking | The positional particle-state consumer inventory in the initial draft omitted raw tuple slicing in Python studies and tests, mesh-derived `prtcl_d` / `prtcl_all`, and the exact `dparh` interpretation. | State-layout edit set and `CP-2 Data Contract`. | Keep the legacy prefix stable, append new-mode fields, update every positional consumer deliberately, and rerun consumer search plus tests. |

### Expanded Particle-State Consumer Inventory

| Surface | Positional or semantic dependency | Required treatment |
| --- | --- | --- |
| `src/athena.hpp` | Defines the shared `IPX` through `IPDB` prefix. | Preserve indices `0` through `13`; append new-mode fields only. |
| `src/particles/particles.cpp` | Allocates fixed legacy width and inspects restart size. | Select width by explicit mode and versioned restart schema. |
| `src/pgen/part_random.cpp` | Initializes and reloads every legacy real slot positionally. | Initialize appended relativistic state and route restart decoding through the versioned schema. |
| `src/outputs/rst_prtcl.cpp` | Writes fixed seventeen-`Real` records and casts identity integers through `Real`. | Preserve legacy-v1 writes only for legacy mode; add typed v2 records for the new mode. |
| `scripts/particles/cr_tracer_inspect.py` | Decodes fixed records and rounds identity values. | Decode v1 and typed v2 explicitly; reject malformed metadata, lengths, checksum, and unsafe identity records. |
| `scripts/particles/cr_tracer_accuracy_study.py` | Slices position, velocity, field, and displacement tuples by fixed offsets. | Preserve prefix behavior and add explicit helpers for appended relativistic fields. |
| `tst/test_suite/particles/cr_accuracy_utils.py` | Slices the same prefix in regression helpers. | Preserve prefix behavior and add explicit relativistic helpers. |
| `src/outputs/df_prtcl.cpp` | Interprets legacy velocity, energy proxies, moments, and samples. | Preserve legacy meanings; add explicit relativistic names or reject ambiguous requests. |
| `src/outputs/dxhist_prtcl.cpp` | Uses displacement prefix including `IPDB`. | Preserve physical displacement; document `IPDB` as accumulated projection onto each sampled midpoint field. |
| `src/outputs/derived_variables.cpp` and `src/outputs/basetype_output.cpp` | Provide mesh-derived `prtcl_d` / `prtcl_all`. | Audit as legacy density-style derived products; do not claim relativistic phase-space coverage. |
| `src/bvals/bvals_part.cpp` | Copies runtime-sized real and integer arrays. | No layout-specific rewrite expected; prove appended-state migration with MPI and AMR tests. |

### Accepted Normalized Code-Unit Contract

`DR-001` is accepted for the first bounded model:

```text
cE = -u_fluid x B
gamma = sqrt(1 + |w|^2 / C_model^2)
v = w / gamma
dw/dt = alpha_s [cE + v x B]
k_model = (gamma - 1) C_model^2
dk_model/dt = alpha_s cE dot v
Delta W = alpha_s h cE_mid dot (w_n + w_n+1) / (gamma_n + gamma_n+1)
```

Here `w = (p/m) / V0`, `C_model = c_model / V0`, and `alpha_s` is a separate
finite signed species coefficient.  For a physical-light mapping,
`C_model = c_phys / V0` and
`alpha_s = q_s B0 T0 / (m_s c_phys)`.  A reduced `C_model` is an explicit model
approximation.  Changing `C_model`, rescaling `alpha_s`, or doing both are
distinct models and must never happen implicitly.

The discrete work increment is accumulated directly from the sampled field,
not reconstructed from kinetic-energy change.

### Accepted Initial Design Choices

These choices supersede narrower unresolved statements above:

| Record | Accepted first-feature choice |
| --- | --- |
| `DR-002` | Preserve the legacy real prefix. Append authoritative `w`, sampled `cE`, explicit accumulated work, signed `alpha_s`, and applicability diagnostics. Derive `gamma`; do not serialize a mutable `gamma` cache. A compatibility velocity shadow may be synchronized only through one invariant helper because existing legacy-prefix consumers require physical velocity. |
| `DR-004` | Add explicit pusher mode `relativistic_hc`; preserve legacy `boris` unchanged. |
| `DR-005` | Add separate finite signed new-mode `alpha_s`; allow both signs with tests; reject zero initially. Keep legacy `IPM` and its negative alignment sentinel legacy-only. |
| `DR-006` | First qualify trilinear interpolation only. Gather Newtonian primitive fluid velocity and cell-centered `B` with the same weights, then construct pointwise `cE = -u_fluid x B`. Reject `lin_legacy`; defer TSC widening. |
| `DR-007` | Qualify fields frozen at `t^n`, spatially regathered at each particle-substep midpoint. Use `time/evolution = kinematic` with an explicit frozen-background acknowledgement because the driver does not execute particle tasks for `static`. Reject dynamic-MHD claims until a separately reviewed stage-coupled widening. |
| `DR-008` | Support 3-D only initially. |
| `DR-012` | Fail on nonfinite state, `C_model <= 0`, invalid or zero `alpha_s`, unsupported dimensions or modules, `|v| >= C_model`, `|u_fluid| >= C_model`, and production `|cE| >= C_model |B|` when `B != 0`. Warn at a preregistered near-limit threshold. |
| `DR-013` | Production support is passive Newtonian ideal-MHD pointwise reconstruction only. Prescribed uniform fields are `part_random` test-harness-only. Never reuse CT edge EMFs. |
| `DR-014` | Use drift-kick-drift: half drift with `v(w_n)`, midpoint gather, full Higuera-Cary candidate momentum push, half drift with `v(w_n+1)`, then accumulate physical displacement. |
| `DR-015` | Accumulate and restart the explicit sampled-field `Delta W` quadrature from the accepted contract. |
| `DR-016` | Refresh the next outer particle timestep bound after relativistic kicks and again after AMR topology changes before the driver calls `Mesh::NewTimeStep()`. |
| `DR-017` | Support strictly periodic coordinates only initially. |
| `DR-018` | Initial runtime boundary is 3-D, strictly periodic, passive, frozen-background Newtonian ideal-gas MHD. Reject missing MHD, isothermal EOS, unacknowledged kinematic or any dynamic evolution, active forcing, ion-neutral, resistivity, orbital advection, shearing box, SRMHD, GRMHD, feedback, Hall paths, and unlisted widenings. |

### Still-Proposed Decisions

| Record | Required closure before dependent edit |
| --- | --- |
| `DR-003` | Accept Higuera-Cary only after spike-v2 evidence and independent re-review, including robust-root handling. |
| `DR-009` | Accept typed-v2 restart layout, inline magic/version/field counts/checksum/time/cycle/timestep/rank topology, sidecar manifest, and deterministic MPI rank-count-change rejection after restart-specialist review. |
| `DR-010` | Accept the exact explicit relativistic output dictionary and constructor-level rejection of every ambiguous unsupported request after output-specialist review. |
| `DR-011` | Accept a conservative first subcycle implementation after adversarial review: `C_model` crossing bound, midpoint regather, strict cap abort, relativistic gyro-angle bound, electric-impulse bound, post-kick refresh, and diagnostic reporting. |

### Fresh Baseline Evidence Recorded So Far

| Gate | Command summary | Result |
| --- | --- | --- |
| Legacy particle CPU plus accuracy | One clean CPU build followed by `pytest -q test_particles_cr_cpu.py test_particles_cr_accuracy_cpu.py` in the repository test harness. | `32 passed in 6.73s` |
| Legacy particle MPI CPU plus accuracy | One clean MPI CPU build followed by `pytest -q test_particles_cr_mpicpu.py test_particles_cr_accuracy_mpicpu.py` in the repository test harness. | `9 passed in 2.19s` |
| Style and whitespace | `cd tst && python run_test_suite.py --style`; `git diff --check` scoped to Phase 0 artifacts. | First style run exposed a spike-only whitespace defect. The expression was clarified. Rerun: `2 passed in 14.48s`; scoped whitespace check passed. |
| Deep AMR `divB` | `cd tst && python run_tests.py mhd/mhd_divb_amr --log_file /tmp/cr-rel-baseline-divb-amr-rerun.log` | Passed: `1 out of 1 test passed`; elapsed `138 s`. An earlier parallel attempt was invalidated when the style harness cleaned the shared build directory during execution. |

### Phase 0 Update-1 Verdict

`RG-001`, `CP-0 Baseline`, and `CP-1 Physics` remain `HOLD`.

The next allowed edits remain ledger updates and isolated evidence artifacts.
Do not open production source edits until the deep-AMR rerun, spike-v2 artifact,
dirty-state refresh, and reviewer recheck are archived.

## Phase 0 Review Update 2

This append-only update records the specialist audits requested by Update 1.
It does not authorize production edits by itself.  The candidate opening-gate
verdict remains subject to a fresh independent recheck after the final spike-v2
checksums are archived.

### Spike-V2 Scope And Selection Evidence

The expanded standalone experiment compares three methods:

| Method | Primary-reference basis | Phase 0 role |
| --- | --- | --- |
| Relativistic Boris | Standard comparison map | Fixed comparison control |
| Vay | Vay (2008), equations 9 through 12 | Crossed-field cancellation comparator |
| Higuera-Cary | Higuera and Cary (2017), equations 20, 21, and 24 | Candidate production pusher |

The script emits its preregistered acceptance criteria before executing the
measurements.  The acceptance set now includes pure-electric acceleration and
deceleration, low-speed agreement, explicit sampled-field work closure,
crossed-field force cancellation, a finite-difference momentum-map Jacobian
probe, a resonance-neighborhood scan, convergence, reversibility, high-gamma
behavior, and long-time finite-state controls.  The root evaluation uses an
overflow-resistant `hypot` form and the conjugate positive-root expression when
needed.  The artifact remains deliberately bounded: it is normalized
`c = q/m = 1` evidence, not AthenaK integration evidence.

Final checksums after the overflow-resistant rerun:

```text
e78198f0771283b548a0b105e8e5ffe8e34c8f6a6a62849e288ac43b1fdc0e8c  scripts/particles/cr_relativistic_pusher_spike.py
29bdd9140d6cd665c31298cca38b68c3c307468c94ad2c0f0037e14bb0d415e7  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
```

`BF-010` is candidate-closed pending independent recheck.  The first
production implementation selects Higuera-Cary because it
retains the crossed-field cancellation behavior required here while also
retaining phase-volume behavior in the bounded Jacobian probe.  Vay remains an
important comparison method but is not selected because the bounded
resonance-neighborhood probe exposes larger energy error and the Jacobian probe
does not support the selected volume-preserving contract.

### Restart-Schema Specialist Audit

The restart specialist independently confirmed:

- legacy v1 casts integer identifiers through native `Real`;
- slot `1` writes mesh `dt` but is read into particle `dtnew`;
- the current C++ loader does not validate the `.pmeta` checksum;
- rank topology is inferred from filename replacement rather than declared;
- particle restart and mesh restart are not paired checkpoints;
- constructor sizing and `part_random` loading duplicate the decoder contract.

`DR-009` is accepted with this first-feature policy:

1. Preserve byte-compatible v1 writes for legacy `drift` and `boris`.
   Document v1 slot `1` as mesh `dt`; never restore it into particle `dtnew`.
   Recompute the legacy particle bound after load.  Validate finite integral
   counts and identifiers, checked allocation sizes, exact representability,
   and `int32` range before conversion.
2. Reject v1 input for `relativistic_hc`; it cannot recover authoritative
   momentum-like state, work, signed coefficient, or checkpoint pairing.
3. Add typed-v2 relativistic-only shards with canonical little-endian encoding,
   no serialized native structs, a fixed `144`-byte header, fixed typed records,
   explicit schema version, flags, endian marker, integer and real field
   counts, saved rank topology, particle counts, mesh cycle, mesh time, mesh
   `dt`, particle `dtnew`, checkpoint ID, byte counts, payload checksum, and
   header checksum.
4. Store checked `i32le` identifiers and `f64le` real values.  Preserve the
   legacy runtime prefix, append authoritative `w`, sampled midpoint `cE`,
   accumulated work, and signed `alpha_s`, and validate then resynchronize the
   compatibility velocity shadow from `w` after load.  Derive `gamma`.
5. Treat instantaneous topology-dependent applicability values as recomputed
   diagnostics, not authoritative restart fields.
6. Publish shards through temporary-file rename and publish one mandatory
   manifest last.  Require exact shard coverage, relative paths without
   traversal, a paired mesh checkpoint ID, same mesh time/cycle/`dt`, and a
   deterministic `reject_rank_count_change` topology policy.
7. Centralize decoding so constructor sizing, loading, production checks, and
   inspector parity cannot drift.

The exact typed-v2 header and record offsets from the specialist report are the
implementation acceptance oracle.  Any layout change requires a new reviewed
schema-minor or schema-major decision rather than an unrecorded rewrite.

### Diagnostic-Schema Specialist Audit

The output specialist independently confirmed:

- `trk` can overrun or corrupt tracked-subset output and uses an
  accelerator-unsafe host counter;
- `pvtk` emits a species integer payload without a matching VTK scalar header
  and casts identifiers through float;
- legacy `pspec p`, `E`, `logE`, and their aliases remain nonrelativistic
  proxies;
- legacy `psamp p`, `energy`, and `mass` are ambiguous for the new mode;
- the inspector does not yet reject duplicate or malformed sample columns.

`DR-010` is accepted with these rules:

1. Keep `IPVX`, `IPVY`, and `IPVZ` as synchronized physical-velocity shadows.
   Add explicit relativistic sample names for `w`, `gamma`, reduced-model
   kinetic energy, sampled midpoint `cE`, direct accumulated work, signed
   `alpha_s`, and applicability metrics.  Do not serialize mutable `gamma`.
2. Require explicit canonical quantities for relativistic `psamp`, `pspec`,
   and `pspec2`.  Reject missing quantities, legacy aliases, and duplicate
   canonical fields after alias normalization rather than silently changing
   their meanings.
3. Retain `ppd`, `df`, `dxh`, and `drh` as legacy-compatible physical-space
   products.  Retain `dparh` only with the documented midpoint
   `dx dot Bhat` accumulation meaning.  Treat `pmom` as a physical-velocity
   transport diagnostic with explicit schema and mode metadata.
4. Isolate `trk` and `pvtk` repairs as small reviewed compatibility edits.
   After repair they remain visualization outputs, not relativistic acceptance
   evidence.
5. Extend the inspector to retain schema and quantity metadata, reject
   duplicate or malformed columns, decode typed v2, and compare direct work
   increments against reduced-model kinetic-energy increments.

### Subcycling And Timestep Specialist Audit

The adversarial numerics specialist independently confirmed:

- cached sampled particle `B` is stale by construction and cannot drive a new
  force or timestep decision;
- a pre-step velocity bound is unsafe when acceleration can increase the second
  half drift;
- intermediate exchange requires an MPI-global schedule, including neutral
  participation from empty ranks;
- cap clipping is unacceptable for the new mode;
- a field-dependent refresh at the end of `RemapAfterAMR()` is too early,
  because post-AMR boundary repair and primitive reconstruction have not yet
  completed.

`DR-011` is accepted conservatively:

```text
Bmax = max ||B||
Umax = max ||u_fluid||
Emax = Umax Bmax

w_floor     = max(0, ||w_n|| - |alpha_s| Emax |dt|)
gamma_floor = sqrt(1 + w_floor^2 / C_model^2)

Ncell  = ceil(C_model |dt| / (f_cell dx_min))
Nblock = ceil(C_model |dt| / (f_block Lblock_min))
Ngyro  = ceil(|alpha_s| Bmax |dt| / (theta_max gamma_floor))
NE     = ceil(|alpha_s| Emax |dt| / delta_w_max)

nsub = MPI_Allreduce(max(Ncell, Nblock, Ngyro, NE), MPI_MAX)
```

Compute the grid envelopes over the exact trilinear-accessible frozen-field
storage region after boundary fill.  Use global leaf-cell and leaf-block
minimum lengths after AMR.  Require `0 < f_cell, f_block <= 0.5`.
Re-gather fields at every substep midpoint.  Treat cached sampled fields as
diagnostics only.  Defer fractional-momentum and fractional-energy subcycle
bounds because they become singular near rest without an independently
reviewed scale.  If the global request exceeds the configured cap, abort
globally before the first drift; do not clip.

`DR-016` is accepted with one shared dirty-aware
`RefreshRelativisticTimestepBound()`:

```text
dtnew = Ncap min(
  f_cell dx_min / C_model,
  f_block Lblock_min / C_model,
  delta_w_max / (alpha_max Emax),
  theta_max / (alpha_max Bmax))
```

Zero-field terms contribute infinity.  Mark the bound dirty after remap and
refresh after initial boundary fill and primitive reconstruction, after the
final post-kick particle exchange, after AMR boundary repair and primitive
reconstruction, and defensively before `Mesh::NewTimeStep()` consumes it.

### Selected Phase-0 Decisions

The following formerly proposed decisions are now candidate-accepted for the
first bounded implementation:

| Record | Selected first-feature choice | Deliberately deferred alternative |
| --- | --- | --- |
| `DR-003` | Higuera-Cary with robust positive-root evaluation | Vay remains an oracle comparator; standard relativistic Boris remains a compatibility comparison |
| `DR-009` | Legacy v1 preserved for legacy pushers; paired typed-v2 manifest-backed relativistic restart; deterministic rank-count-change rejection | Rank-count redistribution and wider identifiers require a separate reviewed widening |
| `DR-010` | Explicit canonical relativistic output dictionary and constructor-level rejection of ambiguous aliases | Silent reinterpretation is forbidden |
| `DR-011` | Pessimistic MPI-global subcycle count from grid envelopes and speed cap; midpoint re-gather every substep; strict abort on cap overflow | Per-rank or per-particle adaptivity is deferred because collective exchange requires a common schedule |
| `DR-016` | Dirty-aware shared outer-bound refresh at explicit lifecycle points | Cached sampled fields and initialization-only bounds are rejected |

### Phase-0 Candidate Gate Record

The opening-gate candidate is now:

- `CP-0 Baseline`: candidate `PROCEED`; fresh CPU, MPI CPU, style, whitespace,
  deep-AMR, dirty-state attribution, consumer inventory, target strategy, and
  specialist audits are archived.
- `CP-1 Physics`: candidate `PROCEED`; the normalized equations, charge
  semantics, direct-work identity, spike-v2 pusher choice, robust root policy,
  and bounded oracle plan are archived.
- `RG-001`: candidate `PROCEED`; the first branch remains passive, 3-D,
  strictly periodic, Newtonian ideal-MHD, explicitly frozen-background, and
  opt-in.  Feedback, dynamic-MHD order claims, non-ideal fields, SRMHD, GRMHD,
  nonperiodic boundaries, TSC widening, and rank-count redistribution remain
  outside scope.

The next permitted edit set remains Phase 0 artifacts only until a fresh
physics reviewer and a fresh architecture/test-adversary reviewer inspect this
update and record `PROCEED`.

## Phase 0 Review Update 3

This append-only correction supersedes stale Update-2 artifact hashes and
narrows several candidate statements to the evidence actually established.
It also archives the final baseline reruns and reviewer-disposition closures.
No production source edit is authorized until fresh independent physics,
architecture, and test-adversary rechecks inspect this update and record
`PROCEED`.

### Spike-V2 Final Artifact Set And Scope Correction

The preregistered criteria are now loaded from a separately frozen manifest.
The spike validates the exact criterion identifiers, metrics, operators, and
thresholds before running measurements, rejects artifact-path aliasing, and
rechecks the manifest checksum before emitting results.  A deliberate
mismatched-manifest probe was rejected before JSON or transcript emission.
All preregistered criteria `P0-01` through `P0-17` passed.

The current final artifact checksums are:

```text
9cc8811a692e8dfdff87b93941fee400ec7d2cb4076a389995e53322860ad9dc  scripts/particles/cr_relativistic_pusher_spike.py
a8659125ecd37f21eb9eaff4dab869b806e31352b72bf3bcaea087c9b3d0af0c  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_preregistered_criteria.json
8458b28006555afe12a3a2f44d593512c70dc2f289d372c4eff0465d1cdec523  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
c4879698fd40a7d2afef2d36e75ed4f9962ddd34521cd32315aee6aebf77da46  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike_transcript.txt
```

Update 2's earlier script and result checksums are superseded.  Its
`overflow-resistant` wording is also too broad.  The bounded Python spike
establishes a cancellation-resistant positive-root evaluation within the
preregistered parameter domain.  Production code must still use scaled or
`hypot`-style norm evaluation and explicitly tested finite-value guards before
claiming robust large-magnitude behavior.

Selected acceptance-oracle metrics include:

| Probe | Measured | Required |
| --- | ---: | ---: |
| Uniform-`B` gamma drift | `1.811201292677658e-15` | `<= 1e-12` |
| Uniform-`B` phase error | `0.002284001491163727` | `<= 0.003` |
| Ordinary crossed-field drift versus ODE | `0.00038356340601482334` | `<= 0.0005` |
| High-gamma phase error | `0.04054875030393133` | `<= 0.05` |
| High-gamma relative endpoint error | `0.04054597243181514` | `<= 0.05` |
| Near-limit departure | `7.850375460079695e-20` | `<= 1e-12` |
| Higuera-Cary direct-work single-step residual | `5.941427905220564e-16` | `<= 2e-15` |
| Higuera-Cary direct-work accumulated residual | `5.756506382681437e-14` | `<= 2e-12` |
| Trapezoid-work finest residual | `6.801550711532656e-09` | `<= 1e-4` |
| Trapezoid-work convergence slope | `2.0005454567407606` | `>= 1.8` |

For the normalized reduced model, the direct discrete-work oracle is:

```text
Delta W = alpha_s h cE_mid dot(w_n + w_np1) / (gamma_n + gamma_np1)
```

The separately accumulated trapezoid-work estimate remains a diagnostic
convergence probe, not a replacement for the exact Higuera-Cary discrete-work
identity.

### Bounded Runtime Contract Correction

Update 2's static-driver language was incomplete because the current static
driver does not execute particle task lists.  The first implementation
contract is instead:

```text
time/evolution = kinematic
particles/relativistic_background = frozen
hydro/mhd rsolver = advect
```

The constructor must reject dynamic evolution and must require explicit
frozen-background acknowledgement.  The accepted first implementation remains
3-D, strictly periodic, passive, Newtonian ideal-gas MHD only.  All previously
listed widenings remain rejected.

For `DR-016`, define:

```text
alpha_max = max_s |alpha_s|
```

The dirty-aware outer timestep helper uses this absolute species maximum.

### Baseline Evidence Archive

The final archived Phase-0 baseline evidence is:

| Gate | Archived evidence | Final result |
| --- | --- | --- |
| Legacy particle CPU plus accuracy | `evidence/phase0_cpu_particle_baseline.log` | `32 passed in 6.79s` |
| Named serial drift restart and inspector path | `evidence/phase0_cpu_restart_inspector_named.log` | `1 passed in 0.28s` |
| Legacy particle MPI CPU plus accuracy | `evidence/phase0_mpicpu_particle_baseline.log` | `9 passed in 2.17s` |
| Named MPI drift restart round trip and inspector path | `evidence/phase0_mpicpu_restart_inspector_named.log` | `1 passed in 0.37s` |
| Style | `evidence/phase0_style.log` | `2 passed in 12.52s` |
| Deep AMR `divB` | `evidence/phase0_deep_amr.log` | `1 out of 1 test passed`; elapsed `139 s` |
| Effective artifact whitespace | `evidence/phase0_untracked_whitespace_check.log` | Passed |
| Dirty-state, branch tips, submodule pin | `evidence/phase0_branch_snapshot.txt` | Archived |
| Evidence checksums | `evidence/phase0_evidence_sha256.txt` | Archived |

Generated logs and CMake-cache evidence were mechanically normalized to remove
trailing whitespace and blank end-of-file noise before the effective
whitespace check.

The refreshed merge-tree probe against target development is archived at
`evidence/phase0_merge_tree.log`.  Its synthetic-tree hash is
`a04bf754c9923049b89e1299f67a721b7ed6652f`.  Update 1's older synthetic-tree
hash is superseded.  The durable direct conflict paths remain:

```text
src/mesh/mesh.cpp
src/mesh/mesh.hpp
src/outputs/outputs.cpp
src/outputs/outputs.hpp
src/pgen/pgen.cpp
```

Additional overlapping paths requiring explicit final integration review are:

```text
src/CMakeLists.txt
src/mesh/mesh_refinement.cpp
src/outputs/basetype_output.cpp
src/pgen/pgen.hpp
tst/scripts/mhd/mhd_divb_amr.py
```

### Consumer-Inventory Correction

The direct-index inventory must include the existing CPU assertions in
`tst/test_suite/particles/test_particles_cr_cpu.py` and the positional slices
in `scripts/particles/cr_tracer_accuracy_study.py`, in addition to the C++
consumers already recorded.  Any Phase-2 state-layout edit must either preserve
the legacy prefix or update every inventoried consumer deliberately.

### Candidate Decision Status Clarification

`DR-003`, `DR-009`, `DR-010`, `DR-011`, and `DR-016` are directionally
accepted first-feature contracts for their dependent phase gates.  This does
not authorize an unreviewed downstream implementation.  Each dependent phase
must still implement the recorded choice, add the prescribed tests, pause at
its reflection point, and obtain its own fresh review.

### Prior Reviewer Disposition Closures

| Prior reviewer concern | Disposition in this update |
| --- | --- |
| Spike criteria were embedded in the executable artifact | Separate frozen criteria manifest, exact contract validation, checksum recheck, and deliberate mismatch rejection added |
| Root claim exceeded bounded evidence | Claim narrowed; scaled production norms and finite guards remain mandatory |
| Static mode did not execute particle task lists | First bounded contract changed to kinematic plus explicit frozen-background acknowledgement |
| `alpha_max` sign semantics were implicit | Recorded as `max_s |alpha_s|` |
| Direct-index test and script consumers were omitted | Inventory corrected explicitly |
| Specialist decisions sounded globally authorized | Status narrowed to directionally accepted contracts with dependent phase gates |
| Archived baseline evidence was not fully enumerated | Final CPU, MPI CPU, named restart, inspector, style, AMR, merge-tree, branch, whitespace, and checksum artifacts listed |

### Phase-0 Recheck Request

The opening-gate candidate remains:

- `CP-0 Baseline`: candidate `PROCEED`
- `CP-1 Physics`: candidate `PROCEED`
- `RG-001`: candidate `PROCEED`

The only permitted next edits remain Phase-0 artifacts until fresh independent
physics, architecture, and test-adversary reviewers inspect Update 3 and each
record `PROCEED`.

## Phase 0 Review Update 4

This append-only correction records and closes the fresh Update-3 adversarial
findings.  It supersedes the invalid general frozen-background runtime claim
from Update 3, hardens the spike manifest contract, and expands the
direct-consumer inventory.  The production-code gate remains closed until the
revised artifacts receive independent rechecks.

### Fresh Reviewer Findings And Dispositions

| ID | Severity | Finding | Disposition |
| --- | --- | --- | --- |
| `BF-012` | blocking | Update 3 incorrectly described `time/evolution = kinematic`, `rsolver = advect`, and an acknowledgement flag as a generally frozen Newtonian ideal-gas MHD background.  The driver still executes MHD RK and CT stages, while the `advect` solver is explicitly isothermal-only. | Closed by the qualification-boundary correction below.  Prescribed frozen fields remain mechanically distinct through Phase 4a.  Solver-coupled behavior is deferred to `CP-4 Solver Coupling`; no general frozen-background production support is claimed. |
| `BF-013` | blocking | The spike verified criterion identifiers, metrics, and operators, but did not pin criterion requirements and thresholds.  Its path-alias guard also missed existing hard links. | Closed by digest-pinning the complete manifest bytes and normalized full criterion contract, checking filesystem identity as well as canonical paths, and archiving threshold-tampering plus direct, symlink, and hard-link negative controls. |
| `BF-014` | blocking for Phase-2 layout edits | The direct-consumer closure table omitted several direct or semantic prefix consumers, including `pos_prtcl.cpp`. | Closed for the opening helper phase by the explicit inventory delta below.  The Phase-2 layout gate must re-audit and deliberately treat every listed consumer. |

### Qualification-Boundary Correction

Update 3's following combination is superseded as a general qualification
contract:

```text
time/evolution = kinematic
particles/relativistic_background = frozen
hydro/mhd rsolver = advect
```

An acknowledgement flag cannot freeze grid fields.  For `kinematic`
evolution, the driver executes the MHD time-integrator task lists, including
flux calculation, Runge-Kutta updates, constrained transport, boundary repair,
primitive reconstruction, and timestep refresh.  The `advect` MHD solver is a
pure-advection solver documented as isothermal-only, so it is not an
ideal-gas qualification contract.

The corrected staged contract is:

1. Phase 1 adds state helpers only and changes no runtime behavior.
2. Phase 2 may add an explicit opt-in parser contract, but it must not silently
   authorize solver-coupled production use.
3. Phase 3 and Phase 4a qualify a mechanically distinct prescribed-field
   harness with fields frozen over the outer step and spatially regathered at
   particle-substep midpoints.
4. Phase 4b may open experimental solver coupling only after `CP-3 Minimal
   Serial` and `RG-005`.
5. `CP-4 Solver Coupling` must separately accept the sampled MHD state, task
   schedule, supported EOS and evolution mode, temporal-order statement,
   prescribed-versus-coupled comparison, and time-dependent manufactured
   convergence evidence before any solver-coupled production claim.

No initial production solver-coupling mode is selected yet.  Credible later
options include a reviewed particle-only execution mode that suppresses MHD
stage evolution, narrowly qualified invariant analytical kinematic
backgrounds, or explicitly evolving Newtonian ideal-MHD backgrounds with
separately measured temporal accuracy.  Choosing among them is a later gated
decision, not an implementation shortcut.

This correction keeps `RG-001` passive and bounded: helper and prescribed-field
work remain one-way particle work with no fluid feedback, current deposition,
Hall correction, resistivity widening, SRMHD, or GRMHD transformation.

### Spike Manifest Integrity Hardening

The standalone spike now:

- pins the complete frozen manifest byte digest;
- pins a normalized full-contract digest over criterion identifier,
  requirement, metric, operator, and threshold;
- continues to validate the explicit metric/operator shape before execution;
- rejects canonical-path aliasing and existing filesystem-identity aliasing;
- repeats the manifest and alias checks immediately before artifact emission.

The archived negative-control transcript is:

```text
evidence/phase0_spike_manifest_negative_controls.log
```

It records passing rejection probes for:

```text
threshold tampering
direct manifest/results alias
symlink manifest/results alias
hard-link manifest/results alias
direct results/transcript alias
symlink results/transcript alias
hard-link results/transcript alias
```

The full spike was rerun after hardening.  All `P0-01` through `P0-17`
criteria passed unchanged.  Current artifact hashes are:

```text
0d09caed911a958e7bd3f25b20c5f3e6f45fbc35015c4f9ffb7dc139eb6fa0d5  scripts/particles/cr_relativistic_pusher_spike.py
a8659125ecd37f21eb9eaff4dab869b806e31352b72bf3bcaea087c9b3d0af0c  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_preregistered_criteria.json
26b47c2e4b8a246687180e24f3275de7ddb31612f945070ff1ca7e9c6e5b89a4  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike.json
feaa2760283386d725d8c593ea21a09fb14c3b9bc44c6bcb7d89ecc213f2a228  docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_pusher_spike_transcript.txt
4709f227026f070c60164f3a59d2722354b4a26bb1acba2debdfd987ff866cf5  docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase0_spike_manifest_negative_controls.log
```

Update 3's script, result, and transcript hashes are superseded.

### Direct-Consumer Inventory Delta

The Phase-2 state-layout review must preserve the existing prefix and
deliberately treat the following additional direct or semantic consumers:

| Surface | Dependency | Required treatment |
| --- | --- | --- |
| `src/particles/particles_pushers.cpp` | Legacy Boris, drift, subcycling, sampled-`B`, displacement, and `IPM` paths consume prefix fields directly. | Keep legacy paths unchanged. Add a separate reviewed relativistic path and explicit helper-mediated shadow synchronization. |
| `src/outputs/pos_prtcl.cpp` | Copies runtime-sized real and integer particle arrays into position output storage. | Audit raw appended-state exposure and metadata deliberately before treating the format as a relativistic diagnostic. |
| `src/outputs/track_prtcl.cpp` | Reads legacy position and velocity shadows directly and has pre-existing structural defects. | Repair in an isolated reviewed Phase-7 change or keep excluded from relativistic acceptance evidence. |
| `src/outputs/vtk_prtcl.cpp` | Copies runtime-sized arrays and has an incomplete integer-field naming contract. | Repair in an isolated reviewed Phase-7 change or keep excluded from relativistic acceptance evidence. |
| `tst/test_suite/particles/test_particles_cr_cpu.py` | Contains direct legacy-prefix assertions through restart summaries. | Preserve the legacy prefix and add separate explicit relativistic assertions. |

The previously recorded consumers remain binding: `src/athena.hpp`,
`src/particles/particles.cpp`, `src/pgen/part_random.cpp`,
`src/outputs/rst_prtcl.cpp`, `src/outputs/df_prtcl.cpp`,
`src/outputs/dxhist_prtcl.cpp`, `src/outputs/derived_variables.cpp`,
`src/outputs/basetype_output.cpp`, `src/bvals/bvals_part.cpp`,
`scripts/particles/cr_tracer_inspect.py`,
`scripts/particles/cr_tracer_accuracy_study.py`, and
`tst/test_suite/particles/cr_accuracy_utils.py`.

### Phase-0 Revised Recheck Request

The candidate opening gate is now:

- `CP-0 Baseline`: candidate `PROCEED`
- `CP-1 Physics`: candidate `PROCEED`
- `RG-001`: candidate `PROCEED` for helper and prescribed-field work only

If the independent rechecks concur, the next allowed edit set is Phase 1 state
helpers only.  Solver-coupled runtime support remains closed until its
downstream checkpoint.

## Phase 0 Review Update 5: Accepted Opening Checkpoint

This append-only checkpoint records the fresh Update-4 rechecks and opens the
smallest Phase-1 helper-only edit set.  It does not authorize runtime
activation, particle-layout edits, prescribed-field integration, or
solver-coupled behavior.

### Final Manifest Negative-Control Expansion

The archived negative-control log now exercises all three artifact-pair
relationships for direct paths, symlinks, and hard links:

```text
manifest / results JSON
manifest / transcript
results JSON / transcript
```

It also retains the threshold-tampering rejection probe.  The superseding log
hash is:

```text
e2258557b66f1f8b126e60b35b68be637a90fa55937312a4f018dccb70cecc3e  docs/source/modules/figures/cr_tracer_relativistic_acceleration/evidence/phase0_spike_manifest_negative_controls.log
```

Update 4's earlier negative-control-log hash is superseded.

### Fresh Independent Recheck Record

| Reviewer lane | Verdict | Findings and disposition |
| --- | --- | --- |
| Physics and numerics | `PROCEED` for Phase 1 helpers only | Recomputed and matched complete-manifest and normalized full-contract digests.  Reconfirmed accepted reduced-model equations, Higuera-Cary map, direct-work identity, bounded-root wording, and pusher-selection rationale.  Production helpers must use scaled norms, finite guards, and later signed-`alpha_s` tests. |
| Architecture and runtime integration | `PROCEED` for Phase 1 helpers only | Reconfirmed that Update 4 retracts the invalid `kinematic + advect` general frozen ideal-gas claim.  Solver coupling remains closed until `CP-4`.  Layout, parser, pusher, gather, restart, output, timestep, and MHD changes remain outside the next write manifest. |
| Adversarial testing | `PROCEED` for Phase 1 helpers only | Reproduced threshold-tampering rejection and alias rejection.  Reconfirmed current checksums, whitespace, branch snapshot, and consumer-inventory delta.  Required a fresh Phase-2 direct-consumer audit before layout edits and retained `BF-009` as a downstream restart blocker. |

### CP-0 Baseline

- Commit reviewed: `3984b57ec8832b3e3e8140e08188a8a05e773136`
- Allowed edit set reviewed: Phase-0 ledger, spike, criteria manifest, results,
  transcript, and archived evidence only
- Evidence reviewed: branch snapshot, prerequisite and target tips, Kokkos
  pin, merge-tree archive, CPU and MPI CPU fixed-energy suites, named serial
  and MPI restart-inspector paths, style, deep-AMR `divB`, spike artifacts,
  manifest-integrity negative controls, effective whitespace audit, and
  evidence checksum index
- Assumptions falsified: `kinematic + advect` is not a general frozen
  ideal-gas MHD production contract; path-resolution equality alone does not
  reject hard-link aliases; the first manifest loader did not pin thresholds
- Scope drift found: none after Update 4
- Unresolved risks: target-development merge conflicts; unavailable GPU
  qualification; downstream legacy restart timestep defect `BF-009`;
  downstream solver-coupled temporal model; mandatory Phase-2 direct-consumer
  re-audit
- Choices made since previous checkpoint: keep initial work prescribed-field
  and helper-only; defer solver coupling; pin complete manifest bytes and full
  normalized criterion contracts; reject direct, symlink, and hard-link
  artifact aliases
- Verdict: `PROCEED`
- Next allowed edit set: Phase-1 state helpers only

### CP-1 Physics

- Commit reviewed: `3984b57ec8832b3e3e8140e08188a8a05e773136`
- Evidence reviewed: accepted normalized equations, primary-reference
  transcriptions, spike criteria manifest, hardened spike script, result JSON,
  transcript, and manifest-integrity controls
- Assumptions falsified: none after Update 4
- Scope drift found: none
- Unresolved risks: production large-magnitude behavior still requires scaled
  norms and finite guards; AthenaK gather, staggering, subcycling, MPI, AMR,
  restart, and output paths remain unqualified
- Choices made since previous checkpoint: retain Higuera-Cary as the bounded
  production candidate; retain relativistic Boris and Vay as comparison
  oracles; keep direct discrete work distinct from trapezoid-work diagnostic
- Verdict: `PROCEED`
- Next allowed edit set: Phase-1 state helpers only

### RG-001: Passive And Bounded Opening Gate

- Date: 2026-05-30
- Commit range reviewed: Phase-0 artifacts on top of
  `3984b57ec8832b3e3e8140e08188a8a05e773136`
- Evidence reviewed: `CP-0`, `CP-1`, Update-4 correction, and all fresh
  independent rechecks above
- Assumptions that still hold: the first branch can remain passive, one-way,
  3-D, strictly periodic, and bounded; helper and prescribed-field
  qualification can proceed without solver coupling
- Assumptions weakened or falsified: general frozen-background support cannot
  be represented by `kinematic + advect` plus an acknowledgement flag
- New risks: a later solver-coupled mode requires an explicit new selection and
  temporal-accuracy evidence
- Blocking findings: none for Phase 1 helpers
- Decision: `proceed`
- Plan changes: Phase 1 remains helper-only.  Prescribed-field work remains
  mechanically distinct through Phase 4a.  Solver coupling remains closed
  until `CP-4 Solver Coupling`.
- Tests added or changed: spike manifest-integrity negative controls expanded
  to all artifact pairs and alias classes
- Independent reviewer findings and responses: `BF-012`, `BF-013`, and
  `BF-014` closed for this gate; downstream obligations retained explicitly
- Next smallest reviewable increment: add a header-only device-capable
  relativistic state helper and isolated custom-pgen host/device-parity harness

### Phase-1 Write Manifest

The next edit set is intentionally narrower than the guide's default Phase-1
surface:

```text
src/particles/relativistic_state.hpp
src/pgen/unit_tests/cr_relativistic_state_test.cpp
inputs/unit_tests/cr_relativistic_state_test.athinput
docs/source/modules/cr_tracer_relativistic_acceleration_ledger.md
Phase-1 evidence artifacts
```

Do not edit `src/athena.hpp`, `src/particles/particles.hpp`,
`src/particles/particles.cpp`, `src/particles/particles_pushers.cpp`, restart
code, outputs, `src/CMakeLists.txt`, or runtime inputs during Phase 1.  The
existing custom-pgen hook accepts nested `PROBLEM=unit_tests/...` paths without
a CMake change.

## Phase 1 Review Update 1: Velocity-Shadow Representability Domain

The first Phase-1 physics review found that reconstructing a compatibility
physical-velocity shadow from arbitrarily large finite `w` is not a
deterministic operation near floating-point saturation.  Mathematically valid
large-`gamma` states can round to `|v| = C_model`, and nearby values can
alternate between acceptance and rejection if the only policy is a post-hoc
strict-speed check.  This append-only decision record is accepted before
revising the helper.

## DR-019: Deterministic Velocity-Shadow Representability Limit

- Status: accepted
- Date: 2026-05-30
- Author: implementation agent after fresh Phase-1 physics review
- Commit: working tree on top of `04090c1a`
- Question: Which finite momentum-like states can be accepted while the first
  implementation retains a synchronized physical-velocity compatibility
  shadow?
- Constraints: the shadow must remain strictly subluminal after rounding;
  acceptance must be deterministic and monotonic rather than dependent on
  rounded reconstruction accidents; the Phase-0 `gamma approximately 1000`
  regime must remain representable where backend precision permits; clipping
  is forbidden.
- Options considered:
  1. Accept every finite `w` and rely only on a strict post-reconstruction
     `|v| < C_model` check.
  2. Clip rounded velocities below `C_model`.
  3. Define a conservative precision-dependent representability domain for the
     compatibility shadow and reject states beyond it explicitly.
  4. Remove the compatibility shadow immediately.
- Selected option: option 3.  Define:

  ```text
  gamma_shadow_max = 1 / sqrt(8 epsilon_Real)
  ```

  `VelocityFromW()`, `ValidateWState()`, and `WFromVelocity()` reject states
  above this cap with an explicit
  `velocity_shadow_unrepresentable` status.  `GammaFromW()` and
  `KineticEnergyFromW()` may still evaluate larger finite momentum states when
  their requested derived quantity remains finite.  Retain the strict
  reconstructed-speed check as a defensive invariant inside the accepted
  domain.
- Evidence and reasoning: near the ultrarelativistic limit,
  `1 - |v|/C_model` is approximately `1/(2 gamma^2)`.  The selected cap reserves
  an approximate four-`epsilon_Real` gap below `C_model`, avoiding the
  non-monotonic saturation region while retaining the Phase-0 high-gamma
  oracle in default double precision and, narrowly, around `gamma = 1000` in
  single precision.
- Why the other options were not selected: option 1 is non-deterministic near
  saturation; option 2 silently changes the physical state; option 4 would
  expand Phase 1 into an unreviewed layout and output migration.
- Compatibility and migration impact: later runtime parsing must reject
  momentum or velocity initializations that cannot produce a trustworthy
  compatibility shadow.  A later momentum-only layout may revisit this cap.
- Validation required: test accepted high momentum, asymptotic gamma and
  kinetic-energy behavior, deterministic acceptance below the cap, rejection
  above the cap, direct invalid-`C_model` conversion cases, nonfinite
  validation inputs, and host/device parity for all status classes.
- Failure signals: non-monotonic status below the selected cap, accepted
  rounded `|v| >= C_model`, clipped values, backend disagreement, or inability
  to retain the Phase-0 high-gamma regime.
- Revisit trigger: remove or redesign the physical-velocity compatibility
  shadow; widen supported precision backends; or qualify a science case above
  the selected cap.
- Independent reviewers: fresh Phase-1 physics reviewer requested recheck
  after implementation
- Follow-up actions: revise the helper and parity harness; record the inherited
  single-precision repository build blocker separately; rerun helper and
  compatibility evidence.

## Phase 1 Review Update 2: Candidate Helper Checkpoint

This append-only update records the post-`DR-019` helper implementation,
evidence, and initial reviewer-finding closures.  Phase 2 remains closed until
fresh physics and integration/test rechecks record `PROCEED`.

### Phase-1 Implementation Scope

The helper increment adds only:

```text
src/particles/relativistic_state.hpp
src/pgen/unit_tests/cr_relativistic_state_test.cpp
inputs/unit_tests/cr_relativistic_state_test.athinput
```

No shared runtime pusher, particle state width, state index, parser, restart,
output, gather, task-ordering, AMR, timestep, MHD, or CMake surface changed.
The nested custom-pgen build uses the existing
`PROBLEM=unit_tests/cr_relativistic_state_test` hook.

The header-only `particles::relativistic` helper now provides device-capable:

```text
ScaledNorm3
MaxVelocityShadowGamma
GammaFromW
VelocityFromW
KineticEnergyFromW
WFromVelocity
ValidateWState
```

The implementation uses scaled norms for large-magnitude state, the
cancellation-resistant kinetic-energy form
`|w| * (|w| / (gamma + 1))`, explicit finite-state status returns, strict
sub-luminal reconstruction, and the accepted `DR-019` compatibility-shadow
representability cap.  It clips no value.

### Phase-1 Test Coverage

The isolated host/device parity harness covers:

- zero momentum and `gamma = 1`;
- nonrelativistic kinetic-energy behavior;
- moderate vector reconstruction;
- large-magnitude scaled-norm and `GammaFromW()` behavior where naive squaring
  overflows;
- deterministic compatibility-shadow rejection for the overflow-prone state;
- accepted `gamma approximately 1000` velocity-shadow reconstruction;
- high-momentum gamma and kinetic-energy asymptotes;
- accepted state at `0.9 * gamma_shadow_max`;
- deterministic rejection at `2 * gamma_shadow_max`;
- `v -> w -> v` and `w -> v -> w` round trips;
- strict `|v| = C_model` and `|v| > C_model` rejection;
- zero, negative, NaN, and infinite `C_model` rejection for direct gamma and
  conversion paths;
- NaN and infinite momentum, velocity, and validation input rejection;
- `DevExeSpace` parity for success and rejection status classes.

### Phase-1 Evidence Archive

| Gate | Archived evidence | Result |
| --- | --- | --- |
| Double-precision CPU helper harness | `evidence/phase1_helper_cpu.log` | Passed: `CR relativistic state helper tests passed` |
| Single-precision CPU qualification attempt | `evidence/phase1_helper_cpu_single_precision_blocked.log` | Blocked before the new helper by inherited narrowing errors in existing repository source including `src/coordinates/cartesian_ks.hpp` and `src/dyn_grmhd/dyn_grmhd.cpp`; not a pass |
| Legacy particle CPU plus accuracy | `evidence/phase1_legacy_cpu.log` | `32 passed in 6.68s` |
| Legacy particle MPI CPU plus accuracy | `evidence/phase1_legacy_mpicpu.log` | `9 passed in 2.20s` |
| Style | `evidence/phase1_style.log` | `2 passed in 13.15s` |
| Branch snapshot | `evidence/phase1_branch_snapshot.txt` | Archived |
| Evidence hashes | `evidence/phase1_evidence_sha256.txt` | Archived |

GPU execution remains unavailable locally and explicitly unqualified.
Single-precision qualification remains blocked by inherited repository source
outside the accepted Phase-1 write manifest.  Neither backend is silently
reported as passing evidence.

### Initial Phase-1 Reviewer Findings And Closures

| Finding | Closure |
| --- | --- |
| Strict reconstructed-speed validation became non-monotonic near floating-point saturation. | Accepted and implemented `DR-019` deterministic `gamma_shadow_max` domain with explicit `velocity_shadow_unrepresentable` status. |
| Required high-momentum coverage was incomplete. | Added gamma and kinetic-energy asymptotes, accepted high-momentum velocity reconstruction, and accepted/rejected shadow-cap boundary cases. |
| Conversion rejection coverage was incomplete. | Added direct and parity cases for zero, negative, NaN, and infinite conversion `C_model`, plus nonfinite validation inputs. |
| `w -> v -> w` round trip was absent. | Added explicit momentum round-trip assertion. |
| Style defects remained. | Added explicit include and wrapped long lines; style rerun passes. |
| Generated caches and logs remained. | Removed generated `test_log.txt`, Python caches, and the invalid non-MPI diagnostic artifact; mechanically normalized retained evidence logs. |

### RG-002: Is The Chosen State Still The Right One?

- Date: 2026-05-30
- Commit range reviewed: working tree on top of `04090c1a`
- Evidence reviewed: helper diff, `DR-019`, double-precision helper harness,
  inherited single-precision blocker, legacy CPU and MPI CPU regressions,
  style, whitespace, branch snapshot, and evidence hashes
- Assumptions that still hold: momentum-like `w` remains the authoritative
  relativistic state candidate; `gamma` should be derived; helper code can
  remain header-only and device-capable; no mutable `gamma` cache is justified
- Assumptions weakened or falsified: a synchronized physical-velocity shadow
  cannot represent arbitrarily high finite `w` deterministically in finite
  precision
- New risks: runtime parsing and pusher integration must propagate the explicit
  `velocity_shadow_unrepresentable` failure; a future momentum-only layout may
  revisit the cap; single-precision and GPU qualification remain unavailable
- Blocking findings: none candidate-closed for Phase 1; independent rechecks
  pending
- Decision: candidate `proceed`
- Plan changes: retain a compatibility velocity shadow only inside the
  reviewed deterministic representability domain; carry `DR-019` into Phase-2
  parser validation
- Tests added or changed: high-momentum asymptotes, deterministic shadow-cap
  boundaries, both conversion round trips, invalid conversion inputs,
  nonfinite validation parity, style rerun
- Independent reviewer findings and responses: initial physics and
  integration/test reviews returned `HOLD`; every local Phase-1 blocker is
  dispositioned above; request fresh recheck
- Next smallest reviewable increment: after recheck only, Phase-2 opt-in
  runtime contract and layout audit

### Phase-1 Candidate Done Claim

- Allowed edit set reviewed: narrow Phase-1 helper manifest only
- Evidence added: listed Phase-1 archive above
- Scope drift found: none
- Unresolved risks: inherited single-precision blocker; unavailable GPU
  execution; downstream direct-consumer layout re-audit; downstream restart
  `BF-009`; downstream solver-coupling gate
- Verdict: candidate `PROCEED`, pending fresh independent physics and
  integration/test rechecks
- Next allowed edit set if accepted: Phase-2 runtime contract only

## Phase 1 Review Update 3: Separate Inverse Velocity-Parsing Domain

The fresh physics recheck found that the accepted forward compatibility-shadow
domain was not closed under inverse conversion near
`MaxVelocityShadowGamma()`.  A rounded velocity shadow can remain strictly
sub-luminal and valid for output while the inverse `v -> w` map is too
ill-conditioned for reliable initialization parsing.

### DR-020: Keep Broad Forward Shadows But Narrow Inverse Parsing

- Status: accepted for Phase-1 final recheck
- Date: 2026-05-30
- Problem: the broad forward `w -> v` output-shadow domain from `DR-019`
  remains deterministic, but using the same cap for `v -> w` parsing permits
  large conditioning error near `|v| / C_model = 1` and does not guarantee a
  stable high-gamma round trip
- Options considered:
  1. Lower the forward shadow cap until the whole compatibility-shadow domain
     is round-trip safe.
  2. Retain the reviewed forward output-shadow cap and introduce a separate
     conservative inverse velocity-parsing cap.
  3. Remove inverse velocity parsing and require momentum-only initialization
     immediately.
  4. Accept the ill-conditioned inverse domain and rely on downstream
     tolerances.
- Selected option: option 2
- Exact rule:
  - `w -> v` output reconstruction retains
    `gamma <= 1 / sqrt(8 * epsilon_Real)`;
  - `v -> w` initialization parsing requires
    `beta < beta_parse_max`, where
    `gamma_parse_max = 0.5 * epsilon_Real^(-1/4)` and
    `beta_parse_max = sqrt((gamma_parse_max - 1) *
    (gamma_parse_max + 1)) / gamma_parse_max`;
  - the inverse guard is strict at `beta_parse_max`;
  - `WFromVelocity()` retains a defensive decoded-gamma guard;
  - values are rejected with `velocity_shadow_unrepresentable`; no value is
    clipped.
- Evidence and reasoning: the inverse map has worsening conditioning as
  `beta -> 1`.  The selected precision-aware inverse cap bounds the leading
  round-trip sensitivity scale `epsilon_Real * gamma^2` at approximately
  `sqrt(epsilon_Real) / 4`, while preserving a materially relativistic
  velocity-initialization range.  Forward shadows beyond the narrower inverse
  cap remain output-only compatibility values.
- Why the other options were not selected: option 1 unnecessarily narrows
  deterministic output reconstruction; option 3 removes a useful bounded
  compatibility path before the parser contract is implemented; option 4
  accepts a known unstable initialization map.
- Compatibility and migration impact: Phase-2 parsing must document that
  velocity initialization is supported only below the inverse parsing cap.
  Higher-gamma initialization must use the authoritative momentum-like state
  after the later data-contract gate.  No legacy runtime behavior changes in
  Phase 1.
- Validation required: direct `WFromVelocity()` cases immediately below, at,
  and above `beta_parse_max`; a high-gamma `w -> v -> w` round trip inside the
  parsing domain with an explicit tolerance; reverse-cap rejection under
  `DevExeSpace`; all six `StateStatus` classes under host/device parity.
- Failure signals: accepted `beta >= beta_parse_max`, clipped values,
  backend disagreement, unstable interior round trip, or accidental narrowing
  of the broader forward output-shadow domain.
- Revisit trigger: remove velocity initialization, remove the compatibility
  shadow, widen precision/backend support, or qualify a science case that
  requires velocity initialization above the selected inverse cap.
- Independent reviewers: physics and integration/test rechecks requested
  after implementation.

### Fresh-Recheck Finding Closure Delta

| Finding | Closure |
| --- | --- |
| Forward output shadows were not closed under high-gamma inverse parsing. | Added accepted `DR-020`, `MaxVelocityParsingGamma()`, `MaxVelocityParsingBeta()`, and strict inverse parsing rejection. |
| Reverse cap-edge tests were absent. | Added direct below, at, and above parsing-boundary tests, a high-gamma interior round trip, and DevExeSpace parity for reverse-cap statuses. |
| Host/device parity omitted `nonfinite_result`. | Added an explicit finite-input `GammaFromW()` overflow case to the parity table. |
| DevExeSpace parity omitted nonfinite `WFromVelocity()` inputs and high-momentum energy. | Added NaN/infinite velocity inputs and high-momentum kinetic energy to the parity table. |
| Phase-1 evidence hashes omitted the uncommitted implementation. | Regenerate the final evidence index with helper, harness, and input hashes before final recheck. |
| Phase-1 branch snapshot lacked full checkpoint provenance. | Regenerate the final snapshot using the Phase-0 provenance structure before final recheck. |
| Ledger control header was stale. | Refreshed the control header to the Phase-1 done-claim boundary and latest accepted Phase-0 commit. |
| Mixed-precision `std::nextafter()` arguments would collapse immediate float neighbors onto the inverse parsing boundary. | Passed typed `Real` zero and one arguments so the strict below, at, and above parsing-boundary tests remain precision-generic after the inherited single-precision repository blocker is repaired. |

### Phase-1 Final-Recheck Evidence Refresh

The final-recheck archive is regenerated after `DR-020` and now binds:

```text
src/particles/relativistic_state.hpp
src/pgen/unit_tests/cr_relativistic_state_test.cpp
inputs/unit_tests/cr_relativistic_state_test.athinput
evidence/phase1_branch_snapshot.txt
evidence/phase1_helper_cpu.log
evidence/phase1_helper_cpu_single_precision_blocked.log
evidence/phase1_legacy_cpu.log
evidence/phase1_legacy_mpicpu.log
evidence/phase1_style.log
evidence/phase1_untracked_whitespace_check.log
```

`evidence/phase1_branch_snapshot.txt` records timestamp, branch, accepted local
HEAD, tracking tip, frozen prerequisite, integration target, Kokkos pin, and
the full short branch status.  The SHA-256 index verifies every bound file.

## Phase 1 Review Update 4: Accepted State-Helper Milestone

The final post-`DR-020` physics/numerics and integration/test rechecks both
record `PROCEED`.  Phase 1 may commit.  The next open edit set is Phase-2
parser-only runtime-contract work; layout widening, restart schema changes,
output changes, prescribed-field gathering, particle acceleration, production
solver coupling, subcycling redesign, and MPI/AMR migration qualification
remain closed behind their owning gates.

### Final Independent Rechecks

| Reviewer | Verdict | Verified closure |
| --- | --- | --- |
| Physics/numerics | `PROCEED` | Typed-`Real` inverse-boundary neighbors, split forward-shadow and inverse-parsing domains, strict rejection without clipping, live isolated helper pass, SHA bindings, snapshot, style, legacy regressions, whitespace hygiene, inherited single-precision blocker, and explicit GPU residual risk. |
| Integration/test | `PROCEED` | Phase-1 scope discipline, typed-neighbor repair, isolated helper rerun, `2 passed in 13.15s` style rerun, legacy `32` CPU and `9` MPI CPU passes, full provenance snapshot, bound SHA index, whitespace hygiene, inherited single-precision disposition, and historical-label clarification. |

### RG-002 Final Disposition

- Decision: `PROCEED`
- Authoritative relativistic state: momentum-like `w`
- Derived values: `gamma`, physical velocity, and kinetic energy
- Cached mutable `gamma`: rejected as unnecessary synchronization risk
- Compatibility-shadow policy: broad deterministic `w -> v` output domain
  from `DR-019`; narrower conservative `v -> w` initialization-parsing domain
  from `DR-020`
- Backend status: double-precision CPU helper qualified; single-precision
  full-build qualification blocked before helper compilation by inherited
  repository narrowing errors; GPU execution unavailable and unqualified
- Next allowed edit set: Phase-2 explicit opt-in parser contract only, with
  legacy particle layout and runtime semantics preserved

## Phase 2 Opening Update 1: Fail-Closed Parser Contract

The accepted Phase-1 milestone is committed at
`0687973f3df831381f8a6ad71ffd3729febdc0fa`.  Phase 2 opens only the explicit
runtime parser contract and its tests.  The particle layout remains the legacy
`nrdata = 14`, `nidata = 3` prefix.  No relativistic particle execution is
authorized in this phase.

### DR-021: Prescribed-Test-Only Parser Boundary

- Status: selected for Phase-2 implementation and independent review
- Date: 2026-05-30
- Problem: expose an explicit relativistic opt-in without silently changing
  legacy `boris`, widening storage, reading incompatible restart payloads, or
  allowing a selectable-but-unimplemented pusher to fall through as a no-op
- Options considered:
  1. Alias the new mode onto legacy `boris`.
  2. Add `relativistic_hc` only to the enum and defer validation.
  3. Add a parser-visible `relativistic_hc` contract that is construction-only,
     prescribed-test-only, and explicitly aborts attempted execution.
  4. Widen layout, restart, outputs, gathers, and pusher behavior together.
- Selected option: option 3
- Exact parser contract:
  - preserve legacy `drift` and `boris` parsing and runtime semantics;
  - require explicit `pusher = relativistic_hc`;
  - require finite `c_model > 0`;
  - require finite signed `alpha_s != 0`;
  - require `nspecies = 1`;
  - require a 3-D strictly periodic mesh and Newtonian coordinates;
  - require ideal MHD;
  - require explicit `interpolation = trilinear`;
  - require `relativistic_field_source = prescribed_test`;
  - require `time/evolution = kinematic` as a particle-task scheduling vehicle
    for the future mechanically isolated prescribed-field harness;
  - use the constructor-required `<mhd>/rsolver = advect` in parser fixtures;
  - reject legacy particle restart before restart-file access;
  - reject unqualified coupled physics blocks and MHD diffusion;
  - abort any attempted `relativistic_hc` particle push with a clear
    unimplemented-or-unqualified message.
- Physics clarification: `kinematic` does not freeze MHD fields.  This parser
  contract authorizes only the planned prescribed-field harness.  Sampled-MHD
  production coupling remains closed until `CP-4 Solver Coupling`.  The
  `advect` solver selection is constructor plumbing for the Phase-2 fixture,
  not qualified positive-cycle ideal-MHD evolution.
- Why the other options were not selected: option 1 silently changes legacy
  semantics; option 2 creates a silent no-op risk in the existing default push
  switch; option 4 violates the reviewed phase boundary and prevents isolated
  parser qualification.
- Compatibility and migration impact: no legacy layout, restart payload,
  output schema, gather, AMR, MPI, or particle timestep behavior changes.
- Validation required: legacy deck regressions, successful construction-only
  parsing for both `alpha_s` signs, clear rejection of each unsupported input
  class, restart rejection before file access, and explicit positive-cycle
  execution abort.
- Failure signals: old deck behavior changes, missing validation path,
  accepted unsupported combination, restart file access before new-mode
  rejection, or positive-cycle silent no-op.
- Revisit trigger: Phase-3 prescribed-field harness implementation or a
  separately reviewed compatibility widening.

### Phase-2 Layout Audit Before Layout Edits

The fresh direct-consumer sidecar audit confirmed that Phase 2 can and should
remain parser-only:

- Keep `ParticlesIndex`, `nrdata = 14`, and `nidata = 3` unchanged.
- Keep legacy restart bytes and positional tuple consumers unchanged.
- Keep outputs unchanged.
- Keep MPI/AMR transport unchanged; it is width-generic but not yet qualified
  for appended relativistic state.
- Treat hard-coded 17-`Real` restart payloads, legacy velocity/energy diagnostic
  aliases, `IPM` sentinel semantics, timestep refresh, `IPDB` displacement
  semantics, and known `ptrack` / `pvtk` findings as downstream gates.

## Phase 2 Review Update 2: Superseding Fail-Closed Boundary

Independent adversarial and documentation reviews held the first Phase-2
candidate.  This append-only update supersedes the incomplete `DR-021` parser
list above without rewriting the historical record.

### DR-022: Reject Ambiguous Staged Configuration

- Status: selected for Phase-2 closure and independent re-review
- Date: 2026-05-30
- Problem: parser-only qualification is unsafe if an input deck can retain a
  staged relativistic key, a stale architecture spelling, an output request, or
  an unqualified MHD widening that is silently ignored.
- Options considered:
  1. Ignore unknown or future-facing keys until the implementation phase that
     consumes them.
  2. Warn and continue.
  3. Fail closed until each spelling, output, and widening is deliberately
     implemented and qualified.
- Selected option: option 3.
- Superseding parser boundary:
  - reject every `output*` block for `relativistic_hc` before initial output
    emission; Phase 2 does not serialize legacy particle artifacts under the
    new mode;
  - reject `relativistic_field_source`, `c_model`, `alpha_s`, and the stale
    `relativistic_background` key when `pusher != relativistic_hc`;
  - reject the superseded `particles/relativistic_background` spelling under
    `relativistic_hc`; select `relativistic_field_source` explicitly instead;
  - reject configured MHD diffusion coefficients and nonzero MHD passive-scalar
    counts;
  - retain the earlier rejection of coupled physics blocks, unsupported
    dimensions, nonperiodic boundaries, relativistic coordinates, static or
    dynamic scheduling, legacy particle restart, non-ideal MHD, non-trilinear
    interpolation, non-prescribed field sources, invalid `c_model`, invalid
    `alpha_s`, and every attempted positive-cycle execution.
- Historical-spelling clarification:
  `particles/relativistic_background = frozen` in earlier architecture notes is
  not an accepted input.  It is superseded by the mechanically distinct
  `particles/relativistic_field_source = prescribed_test` Phase-2 parser
  fixture.  Neither spelling authorizes solver-coupled field sampling.
- Why the other options were not selected: warnings and silent ignores allow a
  deck to appear configured for a relativistic contract while executing a
  different semantic surface.
- Compatibility impact: legacy `drift` and `boris` decks remain valid when they
  do not contain staged relativistic-only keys.  Existing legacy layout,
  restart bytes, output schemas, gathers, MPI/AMR transport, and timestep logic
  remain unchanged.
- Validation expansion:
  - direct new-mode and legacy-mode rejection of the stale background spelling;
  - direct legacy-mode rejection of staged relativistic-only keys;
  - rejection of every output block before file creation;
  - absent-MHD and static-scheduling rejection;
  - direct passive-scalar, Ohmic resistivity, viscosity, conductivity, and
    time-dependent-conductivity rejection;
  - active turbulence forcing and shearing-box/orbital-advection rejection;
  - direct SR and GR coordinate rejection after selecting compatible dynamic
    MHD constructor plumbing;
  - inherited fail-closed rejection for ion-neutral and dynamical-GR inputs
    whose module constructors reject the incompatible deck before particle
    construction.
- Evidence status: focused parser archive refreshed after this update; CPU,
  MPI CPU, style, whitespace, snapshot, and checksum archives must be rebound
  before the Phase-2 done claim.
- Revisit trigger: Phase-3 prescribed-field harness implementation, any output
  schema opening, or any solver-coupled widening.

## Phase 2 Done Claim: Accepted Parser-Only Milestone

Phase 2 is accepted as a bounded parser-only milestone.  It does not authorize
relativistic particle execution, layout widening, restart changes, output
changes, prescribed-field gathering, sampled-MHD coupling, acceleration,
subcycling changes, or MPI/AMR migration claims for appended state.

### Accepted Review Dispositions

| Review | Disposition | Closure |
| --- | --- | --- |
| Adversarial parser matrix | `PROCEED` | Verified missing-MHD, static scheduling, active forcing, shearing-box/orbital-advection, SR, GR, inherited ion-neutral and dynamical-GR fail-closed paths, nonzero passive scalars, all four diffusion arms, stale spelling, output rejection, and legacy controls. |
| `RG-003` documentation gate | `PROCEED` | README and append-only `DR-022` clearly separate constructor plumbing, prescribed-test scope, legacy data and output meanings, stale-spelling rejection, and the positive-cycle abort. |
| Runtime integration done claim | `PROCEED` | Verified parser ordering before output construction and restart-file access, reachability of direct guards, unchanged legacy widths, appended enum ordering, bounded out-of-tree compatibility impact, evidence hashes, and live whitespace hygiene. |

### Accepted Evidence Archive

| Artifact | Result |
| --- | --- |
| `evidence/phase2_parser_contract_cpu.log` | `49 passed`; direct `nlim = 0` construction succeeded; direct positive-cycle execution rejected explicitly. |
| `evidence/phase2_particle_cpu.log` | Full CPU collection: `269 passed, 15 skipped, 67 deselected`. |
| `evidence/phase2_particle_mpicpu.log` | Full MPI CPU collection: `38 passed, 313 deselected`. |
| `evidence/phase2_style.log` | `2 passed`. |
| `evidence/phase2_diff_check.log` | Empty; `git diff --check` passed. |
| `evidence/phase2_untracked_whitespace_check.log` | Every retained untracked fixture, test, and runtime log passed the trailing-whitespace scan. |
| `evidence/phase2_branch_snapshot.txt` | Captures the intended tracked and untracked Phase-2 surface before commit. |
| `evidence/phase2_evidence_sha256.txt` | Binds the production sources, README, fixture, test, ledger, runtime logs, diff check, whitespace report, and branch snapshot. |

### Reflection Gate RG-003: Proceed

- Status: `PROCEED`
- Accepted scope: explicit constructor-only `relativistic_hc` opt-in with
  fail-closed unsupported combinations and unchanged legacy drift/Boris
  semantics.
- Assumptions weakened: stale architecture spellings and staged future keys
  cannot remain silently ignored; every such ambiguity must fail closed until
  deliberately implemented.
- Residual risks: single-precision repository build remains blocked by inherited
  narrowing errors; GPU execution remains unavailable and unqualified; the
  prescribed-field gather helper, runtime diagnostic-only integration, state
  widening, pusher, timestep, restart, outputs, and MPI/AMR migration remain
  downstream gates.
- Next permitted edit set: Phase-3 mechanically isolated trilinear gather helper
  and custom test harness only.  Preserve the positive-cycle
  `relativistic_hc` abort and keep production sampled-MHD coupling closed.

## Phase 3 Opening Update 1: Helper-First Gather Boundary

Phase 2 is committed at `be24fa12`.  Independent architecture and adversarial
test sidecars both held any broad Phase-3 done claim and selected a narrower
first edit set.

### DR-023: Shared Pointwise Trilinear Gather Helper

- Status: selected for helper implementation and independent review
- Date: 2026-05-30
- Problem: interpolate Newtonian primitive fluid velocity and cell-centered
  magnetic field without permitting silent stencil drift, CT-edge reuse,
  face-centered-field reuse, precomputed-field gathering, or premature
  acceleration.
- Options considered:
  1. Interpolate `u_fluid` and `B` through separate copied kernels.
  2. Gather a precomputed grid electric field or CT edge electromotive force.
  3. Extract one device-capable cell-centered trilinear stencil object, reuse
     its exact indices and weights for both vectors, then construct pointwise
     normalized `cE = -u_fluid x B`.
- Selected option: option 3.
- Why the other options were not selected: duplicated interpolation can drift;
  precomputed grid fields and CT edge EMFs answer a different discretization
  question and require a separate reviewed widening.
- Initial helper-only write manifest:
  `src/particles/trilinear_gather.hpp`,
  `src/particles/particles_pushers.cpp`, a custom unit-test pgen, test inputs,
  test utilities, gather-only tests, this ledger, and Phase-3 evidence files.
- Protected surfaces: do not change particle widths, `ParticlesIndex`, restart,
  outputs, MHD tasks, AMR algorithms, subcycling, displacement, Boris algebra,
  or the positive-cycle `relativistic_hc` abort.
- Helper qualification contract:
  - gather arbitrary cell-centered vectors through one reusable stencil;
  - gather primitive `w0(IVX:IVZ)` and cell-centered `bcc0(IBX:IBZ)` with that
    same stencil when the diagnostic-only integration layer opens;
  - name the normalized electric quantity `cE`, not unqualified `E`;
  - keep `trilinear` as the only relativistic gather policy;
  - require post-boundary-repair lifecycle points for ghost-zone access.
- Helper-only validation ladder: host/device parity; uniform fields; affine
  exactness; smooth nonlinear convergence; parallel-field cancellation;
  basis-vector sign table; `cE dot B` debug invariant; noncommutation trap;
  periodic ghost probes; poisoned-ghost negative control; synthetic two-block
  interface; synthetic coarse/fine views; MPI launch parity.
- Diagnostic-only runtime ladder before `RG-004`: actual ghost fill,
  block-boundary continuity, periodic wrap, MPI decomposition including empty
  ranks, SMR behavior, and forced refine/derefine migration with identity
  preservation and no stale sampled state.
- Stop conditions: any source edit outside the manifest; separate `u` and `B`
  stencils; conserved-momentum velocity input; CT-edge or face-centered-field
  reuse; out-of-range stencil access; analytical-sign failure; unexplained MPI,
  block-boundary, periodic, SMR, or AMR disagreement; or removal of the
  positive-cycle runtime abort.
- Revisit trigger: helper review findings or the diagnostic-only runtime
  integration checkpoint.

## Phase 3 Update 2: Independent Helper Review Hold And Correction Scope

The first helper implementation compiled and passed an isolated serial harness,
but an independent reviewer returned `HOLD`.  This is intentionally recorded
before the fixes are accepted.

### Helper Review Finding Record

- Review status: `HOLD`
- Date: 2026-05-30
- Findings:
  1. Host-side parity checks dereferenced device views and were valid only on
     host-accessible serial memory.
  2. The nonlinear convergence ladder varied `u_fluid` and `B` together rather
     than isolating the two required manufactured cases.
  3. Synthetic ghost and interface fixtures did not qualify actual AthenaK
     periodic repair, MPI exchange, SMR prolongation, or AMR remapping.
  4. Boundary probes did not assert stencil index bounds explicitly.
  5. MPI launch parity had not yet been run through a true MPI-enabled build.
  6. The control header still described the accepted Phase-1 checkpoint.
  7. Lower-dimensional helper branches remain outside the first 3-D-only
     relativistic contract.
- Selected correction scope:
  - mirror device fields for host-oracle comparisons;
  - replace captured profile callables with a device-capable profile enum;
  - split nonlinear convergence into varying-`u_fluid`/uniform-`B` and
    uniform-`u_fluid`/varying-`B` sequences;
  - assert every derived stencil index stays inside the allocated synthetic
    ghost extent;
  - retain synthetic interface checks as helper-contract evidence only;
  - run a true `Athena_ENABLE_MPI=ON` launch ladder;
  - keep actual periodic, SMR, and AMR claims open for a separate test-only
    runtime diagnostic harness;
  - refresh the live control header to the accepted Phase-2 SHA.
- Why the scope remains narrow: none of the review findings requires a particle
  layout, restart, output, MHD-task, AMR-algorithm, subcycling, displacement,
  Boris-algebra, or positive-cycle relativistic-execution change.
- Initial corrected serial result:
  - varying-`u_fluid` RMS order: `1.95697`, then `1.98919`;
  - varying-`B` RMS order: `1.95708`, then `1.98922`.
- Open before helper acceptance: true MPI launch results, independent rereview,
  archived evidence, and the separate runtime diagnostic ladder.

## Phase 3 Update 3: Corrected Helper Validation In Progress

- Date: 2026-05-30
- Corrected serial release harness: `PASS`
- Corrected smooth manufactured sequences:
  - varying-`u_fluid`, uniform-`B`: RMS errors `3.89056e-03`,
    `1.00209e-03`, `2.52406e-04`; observed orders `1.95697`, `1.98919`;
  - uniform-`u_fluid`, varying-`B`: RMS errors `2.93677e-03`,
    `7.56365e-04`, `1.90510e-04`; observed orders `1.95708`, `1.98922`.
- True MPI configuration:
  `-D PROBLEM=unit_tests/cr_relativistic_gather_test
  -D Athena_ENABLE_MPI=ON`.
- True MPI launch ladder:
  - rank `1`: `PASS`;
  - ranks `2`: initial one-MeshBlock fixture rejected as intended by AthenaK,
    then `PASS` with `mesh/nx1=8`;
  - ranks `4`: initial one-MeshBlock fixture rejected as intended by AthenaK,
    then `PASS` with `mesh/nx1=16`.
- Claim boundary: these launches qualify helper execution parity only.  They do
  not qualify production MPI particle collection, actual periodic ghost repair,
  SMR prolongation, or AMR remapping.
- Open before helper-only acceptance: debug-bounds build, independent rereview,
  and archived evidence.

## Phase 3 Helper-Only Checkpoint Closure

- Date: 2026-05-30
- Status: `PROCEED` for the mechanically isolated helper slice only.
- Reviewer loop:
  1. Initial review returned `HOLD` for host dereferences of device views,
     combined convergence norms, missing explicit bounds assertions, missing
     true MPI launches, and a stale control header.
  2. Rereview returned `HOLD` for one remaining host dereference of a device
     `cell_cE` view and then for an evidence manifest that did not bind untracked
     helper files.
  3. Final rereview returned `PROCEED` after moving the last gather into a device
     kernel and directly hashing the helper source files, harness input, and
     retained logs.
- Accepted evidence:
  - `evidence/phase3_helper_serial_release.log`;
  - `evidence/phase3_helper_serial_debug_bounds.log`;
  - `evidence/phase3_helper_mpi_launch_parity.log`;
  - `evidence/phase3_helper_diff_check.log`;
  - `evidence/phase3_helper_sha256.txt`.
- Accepted narrow implementation:
  - reuse one cell-centered trilinear stencil for arbitrary gathered vectors;
  - preserve the legacy trilinear Boris magnetic gather through that helper;
  - gather primitive fluid velocity and cell-centered magnetic field with the
    same stencil in helper-only tests;
  - construct explicitly named pointwise `cE = -u_fluid x B`;
  - keep primitive-velocity gather and `cE` construction unused by production
    particle execution.
- Still open before `RG-004`: actual repaired periodic ghosts, live
  block-boundary behavior, production MPI collection including particle-empty
  ranks, static refinement, forced AMR migration shadow evidence, and an
  independent runtime-diagnostic review.
