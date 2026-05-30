---
orphan: true
---

# CR Tracer Relativistic Acceleration Living Ledger

## Ledger Status

This is the living implementation ledger required by
`cr_tracer_relativistic_acceleration_implementation_guide.md`.  It began as the
Phase-0 audited baseline and now records each accepted implementation milestone,
decision, reflection gate, independent review, and retained evidence seal.

The historical Phase-0 entries use these intentionally conservative statuses:

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
| Working branch HEAD at latest accepted commit | `aa8663e8a5d49e26c206363d028b52d0e350a91f` | Refresh after every accepted commit. |
| Frozen feature base | `64a4d1be8da1c22d1328cc47280195b3747fa0ab` from `feature/CR_tracers_followup_architecture` | Change only through an accepted `DR-000` update. |
| Intended eventual integration target | `origin/development` at `c6a73b08e60807f8b925164c5e7edd5cb820c8ae` | Refresh the target SHA and merge-tree audit after target updates and before handoff. |
| Current implementation phase | `Phase 9: GPU disposition, documentation, public overlay, and final handoff` | Keep the temporal claim explicitly first-order for evolving fields until a separately reviewed stage-coupled widening exists. |
| Allowed write manifest | Phase-9 documentation and retained qualification evidence only: `docs/source/modules/particles.md`, bounded handoff documents, append-only ledger updates, and evidence artifacts | Do not widen solver-coupled MPI/SMR/AMR execution, nonperiodic boundaries, changed-rank redistribution, or physics semantics in the handoff phase. |
| Last accepted checkpoint | `CP-6 AMR And MPI PROCEED`: bounded periodic `prescribed_test` MPI, same-level multiblock, static-SMR, adaptive-AMR, and same-topology continuation qualification accepted on `addd12d4` after source and evidence rebounds | Record each milestone commit before opening the next phase. |
| Open blocking findings | Corrected Phase-9 evidence seal and fresh cold-review rebound remain open; public docs preview overlay and workstation GPU disposition are complete | Keep the live list current.  Do not proceed while a blocker for the next edit set remains open. |
| CPU/MPI qualification status | Phase-8 accepted source `addd12d4` plus Phase-9 candidate replay: frozen registration `62/62`; release and debug-bounds migration replay `62/62` across `21` cases each; mutation controls `17/17`, `88/88`, `8/8`, and coupled-oracle controls rejected; release and debug restart, diagnostic, subcycle, coupled, and all-format replays passed; parser `172/172`; legacy CPU plus accuracy `20/20` and `12/12`; legacy MPI CPU plus accuracy `4/4` and `5/5`; style `2/2`; deep AMR `divB` `1/1`; bounded periodic `prescribed_test` MPI/SMR/AMR opening accepted; solver-coupled `mhd_ideal` MPI/SMR/AMR remains fail-closed; single precision remains blocked before CR source by inherited repository narrowing errors | Rerun and archive qualification at each accepted milestone. |
| GPU qualification status | Workstation disposition archived; unavailable locally and explicitly unclaimed | GPU qualification is optional for `MERGE READY` on this workstation, but it remains a separate unqualified residual risk.  Do not claim `GPU QUALIFIED` without accelerator evidence on an accepted SHA. |
| Merge-ready status | Not eligible | Requires Phase-9 qualification, public docs overlay, `RG-010`, and `CP-7 Final Cold Review`. |

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

## Phase 3 Update 4: Test-Only Live-Grid Diagnostic Harness

### DR-024: Finalizer-Only Runtime Sampling Harness

- Status: selected for Phase-3 diagnostic qualification
- Date: 2026-05-30
- Problem: exercise actual repaired `w0` and `bcc0` storage without opening
  positive-cycle `relativistic_hc`, changing particle records, routing through
  production outputs, or editing MHD and AMR tasks.
- Options considered:
  1. Add a production particle diagnostic task beside `Push()`.
  2. Sample during AMR remapping or before boundary repair.
  3. Add one compile-selected custom pgen whose finalizer samples live MHD
     fields after initialization or the final repaired AMR hierarchy.
- Selected option: option 3.
- Reasoning: the existing finalizer runs after initialization ghost repair and,
  for positive-cycle legacy shadow tests, after final MHD stages and AMR repair.
  It observes the real integration surface while remaining mechanically
  distinct from future production sampled-MHD coupling.
- Implementation boundary:
  - call existing `PartRandom()` initialization;
  - prescribe active-cell primitive velocity and convert it back to conserved
    state through the existing EOS;
  - retain face-centered magnetic initialization from `PartRandom()`, because
    later repair reconstructs cell-centered `bcc0`;
  - gather live `w0(IVX:IVZ)` and `bcc0(IBX:IBZ)` through the reviewed helper;
  - construct diagnostic-only `cE = -u_fluid x B`;
  - globally gather rows from every MPI rank, including particle-empty ranks,
    sort by `(species, tag)`, and emit explicit TSV metadata;
  - use a separate analyzer to reconstruct analytic values and compare
    decomposition-invariant physical rows.
- Deliberately excluded: particle-record widening, restart, production output
  registration, MHD task edits, solver-coupled acceleration, dynamic-MHD order
  claims, and positive-cycle `relativistic_hc`.

### DR-025: Separate Hierarchy Oracles By Question

- Status: selected
- Date: 2026-05-30
- Observation: an affine velocity profile is not periodic.  Using it as a
  strict live-grid oracle on a periodic domain correctly exposed wrapped-ghost
  discontinuities near the domain boundary.
- Selected split:
  - retain affine exactness across synthetic same-level and coarse/fine helper
    views, where the interpolation coordinate contract is isolated;
  - use a smooth periodic manufactured `u_fluid` and the existing
    `sinusoidal_divb_free` magnetic profile for live periodic and SMR
    qualification, with a preregistered `4e-2` pointwise profile-error bound;
  - use deterministic tags immediately on both sides of periodic faces,
    same-level interfaces, and each coarse/fine interface plane;
  - keep dynamic AMR as a separate moving-particle lifecycle shadow: require
    exact identity preservation, nonzero motion for every tag, final repaired
    sampled-field construction, serial/MPI parity, and recorded mesh churn,
    but do not claim a static analytic MHD-field oracle after the background
    evolves.
- Rejected alternative: loosening tolerances around an incompatible affine
  periodic fixture.  That would hide a test-design error.

## Phase 3 Update 5: Runtime-Diagnostic Review Hold And Corrections

- Date: 2026-05-30
- Independent runtime reviews: `HOLD`, honored before `RG-004`.
- Findings:
  1. Initial live-grid decks used uniform fields, which could mask wrong-neighbor
     and stale-ghost reads.
  2. The analyzer checked values but did not bind deterministic face,
     same-level-interface, coarse/fine-interface, or exact AMR identity
     coverage.
  3. Particle-empty-rank evidence inferred absence from rows instead of binding
     independent per-rank counts.
  4. Decomposition comparison did not structurally validate the reference TSV
     and omitted lifecycle coordinates.
  5. The dynamic AMR shadow proved hierarchy churn but did not require
     particle motion.
  6. Test-only failures used rank-local `exit`, which could strand MPI peers.
  7. The retained archive omitted periodic MPI ranks `1` and `2`, CMake caches,
     exact command transcripts, and a branch snapshot.
- Selected corrections:
  - add smooth periodic manufactured `u_fluid` and use the existing smooth
    divergence-free periodic magnetic profile;
  - add deterministic all-axis periodic and same-level probes plus deterministic
    all-axis coarse/fine probes;
  - emit and validate explicit rank-count metadata;
  - require exact contiguous `(species, tag)` identity sets and a final cycle;
  - validate both decomposition-comparison operands structurally and compare
    `cycle` and `time`;
  - move every AMR-shadow particle and require motion reconstructed
    independently from deterministic tag-random initialization;
  - use `MPI_Abort` for test-harness failures after MPI initialization;
  - retain serial and MPI CMake caches, exact set-`-x` transcripts, all compared
    TSVs, negative controls, legacy regressions, style, branch snapshot, diff
    check, and direct SHA-256 manifest.
- Deliberate boundary: the dynamic AMR shadow is lifecycle and
  decomposition-parity evidence.  It is not a dynamic-MHD temporal-order claim.
  The static periodic and SMR manufactured cases own analytic repaired-grid
  field qualification.

### Corrected Runtime Validation In Progress

- Serial periodic construction-only diagnostic:
  all-axis periodic-face, corner, and same-level-interface deterministic probes,
  nonuniform periodic `u_fluid` and `B`, `64` rows, bounded analytic analyzer
  `PASS`; observed max absolute errors are `2.2877e-3` for `u_fluid` and
  `2.3893e-2` for `B`, below the preregistered `4e-2` bound.
- Serial SMR construction-only diagnostic:
  deterministic probes on both sides of all six coarse/fine planes, level-1
  static refined region, nonuniform periodic `u_fluid` and `B`, `512` rows,
  bounded analytic analyzer `PASS`; observed max absolute errors are
  `5.0446e-4` for `u_fluid` and `1.9321e-3` for `B`.
- Serial dynamic AMR shadow:
  moving legacy drift particles, exact contiguous identity set, independently
  reconstructed nonzero motion for every tag, `56` MeshBlocks created and `21`
  deleted, `64` rows, final sampled-field construction analyzer `PASS`.
- MPI periodic construction-only diagnostic:
  ranks `1`, `2`, and `4`, explicit rank-count metadata, `64` rows each,
  bounded analytic analyzer and serial physical-row comparison `PASS`.
- MPI particle-empty-rank diagnostic:
  ranks `6`, uneven `8`-MeshBlock decomposition, explicit rank counts
  `0,0,0,0,2,0`, `2` global rows, particle-empty ranks participate in
  collection, bounded analytic analyzer `PASS`.
- MPI SMR construction-only diagnostic:
  ranks `4`, explicit rank-count metadata, deterministic coarse/fine probes,
  `512` rows, bounded analytic analyzer and serial physical-row comparison
  `PASS`.
- MPI dynamic AMR shadow:
  moving legacy drift particles, ranks `4`, `56` MeshBlocks created, `21`
  deleted, `11` communicated for load balancing, exact contiguous identity set,
  independently reconstructed nonzero motion for every tag, `64` rows, final
  sampled-field construction analyzer and serial physical-row comparison
  `PASS`.
- Analyzer adversarial controls:
  collapsed periodic probes, collapsed SMR probes, shifted AMR tags, false
  serial-as-six-rank evidence, duplicated comparison reference rows, and
  falsified rank-count metadata all reject as intended.
- Fresh shared-path regressions:
  - legacy CPU, fixed-energy accuracy, and parser contract: `81 passed`;
  - legacy MPI CPU and fixed-energy accuracy: `9 passed`.
- Open before `RG-004`: style refresh, retained evidence manifest, branch
  snapshot, and fresh independent rereview of the corrected runtime package.

## Phase 3 Update 6: Adversarial Oracle Hold And Rebound

- Date: 2026-05-30
- Independent runtime rereviews: `HOLD`, honored before `RG-004`.
- Findings:
  1. A forged serial TSV could pad `rank_counts` with zero entries and satisfy
     the empty-rank oracle without a separately emitted participation witness.
  2. AMR serial/MPI comparison bound equal lifecycle times but did not bind the
     final time independently.  A coherent time rewrite of both operands could
     pass.
  3. The negative-control transcript invoked mutated `/tmp` TSVs without
     retaining the fixture-construction commands.
  4. The control header and tail still described completed style, manifest, and
     branch-snapshot packaging as pending.
- Selected corrections:
  - gather the communicator rank independently of particle-row counts and emit
    an explicit `rank_participants` witness with `world_size`;
  - require contiguous participant identities and exact agreement between the
    participant witness and rank-count width;
  - add `--expected-time` with an explicit tolerance and require it for the
    AMR serial and MPI lifecycle-shadow validations;
  - regenerate the negative-control transcript with traced fixture
    construction commands, including padded serial rank counts, malformed
    participant witness, and coherent serial/MPI `time=999` mutations;
  - rerun every retained serial and MPI runtime case through the corrected
    producer and analyzer;
  - refresh full style, CMake caches, build transcript, diff check, direct
    untracked-artifact whitespace scan, branch snapshot, and SHA-256 manifest.
- Rebound validation:
  - serial periodic, serial SMR, and serial dynamic-AMR-shadow analyzers pass;
  - MPI periodic ranks `1`, `2`, and `4`, MPI-`6` particle-empty-rank,
    MPI-`4` SMR, and MPI-`4` dynamic-AMR-shadow analyzers pass;
  - the six-rank artifact records rank counts `0,0,0,0,2,0` and an independent
    participant witness `0,1,2,3,4,5`;
  - AMR validations bind cycle `3` and final time
    `0.025970518172472944`;
  - `8` traced adversarial rejection cases are retained, including both newly
    identified attacks;
  - full style passes: `2 passed`.
- Deliberate boundary: the participant witness and retained `mpirun -np 6`
  transcript make the MPI-empty-rank claim auditable.  The TSV remains an
  evidence artifact, not a cryptographically unforgeable attestation.
- Open before `RG-004`: refreshed direct SHA-256 manifest and fresh independent
  rereview of this rebound package.

### DR-026: Phase-4a Prescribed-Test Layout And Initialization Candidate

- Status: candidate-selected; opens only after `RG-004` records `PROCEED`
- Date: 2026-05-30
- Question: What is the smallest explicit state-layout and initialization
  contract needed for prescribed-field-only relativistic acceleration?
- Options considered:
  1. Replace legacy prefix velocity slots with authoritative momentum.
  2. Append authoritative momentum, prescribed sampled fields, direct work,
     signed coefficient, and stored applicability metrics immediately.
  3. Append only authoritative momentum, prescribed sampled fields, direct
     work, and signed coefficient; defer recomputed applicability metrics to
     the diagnostics phase.
- Selected option: option 3.
- Reasoning:
  - preserve real-data prefix indices `0..13` unchanged;
  - append `IPWX=14`, `IPWY=15`, `IPWZ=16`, `IPCEX=17`, `IPCEY=18`,
    `IPCEZ=19`, `IPWORK=20`, and `IPALPHA=21`;
  - use `nrdata=22` only for `relativistic_hc`; retain legacy `nrdata=14` and
    common `nidata=3`;
  - derive `gamma`; do not cache or serialize it;
  - synchronize `IPVX`, `IPVY`, and `IPVZ` compatibility shadows from
    authoritative `w` through one invariant helper after initialization and
    every kick;
  - retain deterministic legacy-prefix `IPM` initialization but never read it
    in relativistic algebra;
  - defer topology-dependent or instantaneously recomputable applicability
    metrics until the diagnostics phase rather than expanding authoritative
    restart state prematurely.
- Initialization contract:
  - add explicit test-harness selector
    `problem/relativistic_initial_state = velocity | w`;
  - support `problem/w0x`, `problem/w0y`, and `problem/w0z` for authoritative
    momentum initialization;
  - support `problem/cE0x`, `problem/cE0y`, and `problem/cE0z` only for the
    mechanically isolated `prescribed_test` harness;
  - reuse `problem/B0x`, `problem/B0y`, and `problem/B0z` as the prescribed
    magnetic-field source, stored in legacy-prefix `IPBX`, `IPBY`, and `IPBZ`;
  - prohibit routing Phase 4a through `w0`, `bcc0`, the Phase-3 live-grid
    diagnostic gather, or CT edge EMFs;
  - reject `relativistic_initial_state`, `w0x`, `w0y`, `w0z`, `cE0x`, `cE0y`,
    and `cE0z` when `pusher != relativistic_hc`;
  - reject mixed initialization modes, nonfinite state, inverse parsing beyond
    the velocity-to-`w` cap, and forward shadows beyond the `w`-to-velocity
    cap.
- Phase-4a write manifest reopened exactly for:
  - `src/athena.hpp`;
  - `src/particles/particles.cpp`;
  - `src/pgen/part_random.cpp`;
  - `src/particles/particles_pushers.cpp`;
  - one narrowly scoped invariant helper;
  - prescribed-field kernel tests and decks;
  - this ledger and Phase-4a evidence artifacts.
- Still excluded during Phase 4a: restart writers and decoders, particle
  outputs, MHD tasks, AMR infrastructure, live-grid coupling, CT EMF reuse,
  subcycling redesign, outer-timestep redesign, and solver-coupled execution.
- Compatibility impact: legacy drift and Boris records remain byte-for-byte
  runtime-width compatible.  Restart and outputs remain rejected for
  `relativistic_hc` until their owning phases.
- Validation required: isolated prescribed-field kernel harness, positive and
  negative `alpha_s`, state-cap rejections, shadow synchronization, direct-work
  closure, reversibility, convergence, and retained legacy regressions.
- Failure signals: any live-MHD gather in Phase 4a, any read of `IPM` in
  relativistic algebra, any mutable `gamma` cache, any ambiguous output or
  restart widening, or any legacy record-width change.
- Revisit trigger: `RG-005`, restart-schema implementation, diagnostics-schema
  implementation, or solver-coupling design.

## RG-004: Pointwise Ideal Construction And Prescribed-Test Opening

- Date: 2026-05-30
- Commit range reviewed: Phase-3 helper commit
  `56dcbf3f3c0120171604af4db9300bcde86b2397` plus the bound runtime-diagnostic
  package and Update-6 rebound corrections.
- Evidence reviewed:
  - retained direct SHA-256 manifest;
  - serial periodic, SMR, and dynamic-AMR-shadow runtime transcripts;
  - MPI periodic ranks `1`, `2`, `4`, particle-empty rank-`6`, SMR rank-`4`,
    and dynamic-AMR-shadow rank-`4` runtime transcripts;
  - full style, legacy CPU, legacy MPI CPU, CMake-cache, branch-snapshot,
    whitespace, and diff-check evidence;
  - reconstructible negative-control transcript with `8` traced rejection
    cases.
- Assumptions that still hold:
  - gather Newtonian primitive fluid velocity and cell-centered `B` with the
    same trilinear stencil before constructing pointwise
    `cE = -u_fluid x B`;
  - do not substitute CT edge EMFs or add a particle-side curl;
  - retain a mechanically distinct prescribed-field harness before opening
    solver coupling;
  - dynamic-AMR evidence is a lifecycle shadow, not dynamic-MHD temporal-order
    qualification.
- Assumptions weakened or falsified:
  - row-count metadata alone was not a sufficient empty-rank oracle;
  - pairwise lifecycle-time equality alone was not a sufficient AMR-time
    oracle;
  - unexplained temporary negative-control fixtures were not reviewable
    evidence.
- New risks:
  - retained TSVs are auditable evidence artifacts, not cryptographic
    attestations.  Their narrow claims depend on the bound producer,
    transcripts, and direct manifest.
- Blocking findings: none after the Update-6 rebound and focused DR-026
  recheck.
- Independent reviewer findings and responses:
  - fresh physics/oracle reviewer: `PROCEED`; replayed all nine positive
    analyzers and checked the participant witness, expected-time binding, and
    reconstructible controls;
  - fresh integration-boundary reviewer: initial `HOLD` until DR-026 named the
    reopened files, prescribed `B` source, live-grid prohibition, and
    fail-closed legacy key policy;
  - focused integration-boundary recheck: `PROCEED` after all three DR-026
    amendments.
- Decision: `PROCEED`.
- Plan changes:
  - close Phase 3;
  - open Phase 4a only for the mechanically isolated prescribed-test
    relativistic Higuera-Cary push and the exact DR-026 write manifest;
  - keep solver-coupled execution closed behind `CP-3 Minimal Serial`,
    `RG-005`, and later `CP-4 Solver Coupling`.
- Next smallest reviewable increment: append the explicit Phase-4a state
  fields, add one invariant helper for velocity shadows, initialize
  prescribed-test state fail-closed, integrate the prescribed HC map, and add
  isolated analytical kernel tests.  Do not edit restart, outputs, subcycling,
  outer timestep, MHD tasks, or AMR infrastructure.

## DR-027: Phase-4a Prescribed-Test HC Implementation Choices

- Status: candidate-selected; pending `RG-005`
- Date: 2026-05-30
- Question: Which bounded implementation choices make the first executable
  relativistic push reviewable without leaking into solver coupling or later
  schema work?
- Options considered:
  1. Use direct vector arithmetic and rely on post-operation finite checks.
  2. Add checked scalar, scale, add, dot, and cross helpers and reject
     arithmetic range violations before partial particle-state mutation.
  3. Introduce a general vector-arithmetic abstraction across legacy and
     relativistic pushers immediately.
- Selected option: option 2.
- Reasoning:
  - direct post-operation checks can observe overflow only after it occurs;
  - narrowly scoped checked helpers keep the legacy path unchanged;
  - checked dot and cross operations avoid naive intermediate overflow;
  - the negative-`sigma` Higuera-Cary root uses the conjugate form;
  - the runtime path computes candidate work, positions, displacements,
    `|B|`, and `IPDB` first, then commits the accepted state.
- State-invariant choice:
  - `IPWX`, `IPWY`, and `IPWZ` remain authoritative;
  - one helper stores accepted `w` and synchronizes `IPVX`, `IPVY`, and
    `IPVZ` compatibility shadows after initialization and every accepted kick;
  - no relativistic algebra reads legacy `IPM`;
  - a paired production-path poison test changes `IPM` while holding the
    final-limited outer step fixed and requires identical physics fields.
- Initialization choice:
  - require explicit complete `B0x`, `B0y`, and `B0z` prescribed tuples;
  - require explicit `cE0x`, `cE0y`, and `cE0z`;
  - for `relativistic_initial_state = velocity`, require
    `particle_velocity = uniform`, explicit complete `v0x`, `v0y`, and `v0z`,
    and no scalar `v0`;
  - for `relativistic_initial_state = w`, require exactly complete `w0x`,
    `w0y`, and `w0z`, with no velocity initialization tuple;
  - retain legacy defaults unchanged outside `relativistic_hc`.
- Harness choice:
  - use a helper-only custom pgen for algebra, root branches, host/device
    parity, and deterministic failure statuses;
  - use a distinct production-path custom pgen that calls `PartRandom`,
    installs a finalizer, emits witness TSV data, and never invokes the HC
    helper directly;
  - validate the production path against independent closed-form and RK4
    oracles in Python;
  - poison finite live `w0` and `bcc0` arrays in a paired production run and
    require identical physics fields, proving Phase-4a isolation from sampled
    MHD arrays.
- Convergence choice and inflection:
  - the first attempted ladder varied mesh CFL values over a final-limited
    `0.04` window; retained traces showed every case took the same single
    `0.04` step, so that ladder was rejected as non-discriminating;
  - probe runs then established an observed refinement axis using paired
    mesh and particle CFL values `0.08`, `0.04`, `0.02`, and `0.01` over a
    `0.2` window;
  - the accepted ladder measures a combined momentum-plus-displacement slope
    of `1.9996235673293787`.
- Documentation reopening:
  - reopen `src/particles/README.md` only to replace stale constructor-only
    language with the bounded prescribed-test execution contract and the
    appended in-memory state layout;
  - keep restart and output semantics explicitly gated.
- Compatibility impact: legacy prefix indices, legacy widths, restart bytes,
  output implementations, MHD tasks, task scheduling, and subcycling remain
  unchanged.
- Alternatives retained:
  - revisit a broader arithmetic abstraction only if later solver-coupled or
    subcycling code demonstrates real duplication;
  - revisit the authoritative state layout only at restart and diagnostic
    schema gates;
  - retain prescribed and solver-coupled evidence as separate classifications.
- Revisit trigger: `RG-005`, solver-coupled Phase 4b design, acceleration-aware
  subcycling, restart schema, or any reviewer finding that the prescribed
  harness is not mechanically isolated.

## RG-005 Candidate: Does The Prescribed-Test Pusher Merit Integration?

- Date: 2026-05-30
- Status: candidate `PROCEED`; independent physics, test-adversary, and
  integration-boundary done-claim reviews still required.
- Scope reviewed:
  - append-only relativistic in-memory state with legacy prefix preservation;
  - arbitrary-`c_model` checked Higuera-Cary helper;
  - prescribed-test-only production path;
  - explicit direct sampled-field work quadrature;
  - fail-closed initialization, output, restart, solver-coupling, and
    subcycling boundaries;
  - README correction and preregistered Phase-4a criteria.
- Prescribed-test evidence observed before review:
  - helper-only serial release and debug runs pass;
  - helper-only host/device parity covers accepted and rejected probes,
    zero field, pure electric field, magnetic reversibility, signed
    `alpha_s`, low-speed limit, high-`gamma`, and conjugate-root execution;
  - parser contract reports `73 passed`;
  - production metrics satisfy all preregistered `P4A-01` through `P4A-16`
    thresholds;
  - direct work closure is `2.168404344971009e-19` in pure `cE`,
    `2.168404344971009e-19` for non-unit `c_model`, and
    `6.776263578034403e-21` in crossed fields;
  - uniform-`B` phase error is `1.1196392133608346e-05`;
  - high-`gamma` relative energy drift is `0.0`;
  - legacy-`IPM` independence and poisoned-live-grid isolation errors are
    both `0.0`;
  - production convergence slope is `1.9996235673293787`;
  - refreshed Phase-3 serial periodic, SMR, and dynamic-AMR-shadow gather
    witnesses pass after their input-deck initialization amendments.
- Retained regression status at record time:
  - legacy core CR CPU suite: `20 passed`;
  - additional accuracy, MPI, style, evidence binding, and independent
    done-claim reviews still running or pending.
- Scope still closed:
  - sampled-MHD solver coupling;
  - production-path fluid-velocity cancellation;
  - dynamic-MHD temporal-order claims;
  - relativistic subcycling and outer-timestep redesign;
  - relativistic restart and output schemas;
  - MPI, AMR, and accelerator qualification for the new momentum-changing
    runtime path.
- Decision: hold the gate open until retained regressions, direct evidence
  binding, and all three fresh reviews complete.

## RG-005 Review Hold And Correction Scope

- Date: 2026-05-30
- Status: `HOLD`; correction and rebound validation in progress.
- Independent review outcome:
  - physics reviewer: `HOLD`;
  - test-adversary reviewer: `HOLD`;
  - integration-boundary reviewer: `HOLD`.
- Blocking findings accepted:
  1. the documented serial-only Phase-4a scope was reachable under MPI and on
     multilevel meshes because the production constructor lacked a fail-closed
     rank and mesh-level guard;
  2. the frozen prescribed-field tuple was reachable with nonuniform
     `B_profile` selections even though the pusher intentionally reads only the
     stored particle tuple;
  3. the helper conjugate-root check asserted branch entry and magnetic-energy
     conservation but did not bind an independently evaluated momentum result;
  4. retained production-path evidence lacked expected-failure subprocess
     cases, a non-unit-`c_model` magnetic rotation case, a high-`gamma` phase
     oracle, and independent Python reconstruction of velocity shadows;
  5. the original frozen JSON bound `P4A-01` through `P4A-16` but did not bind
     all supplemental metrics emitted by the broadened analyzer;
  6. the control-header write manifest omitted the DR-027 README reopening and
     the accepted Phase-3 runtime-diagnostic commit;
  7. direct Phase-4a SHA-256 evidence binding was still pending.
- Correction choices:
  - reject `global_variable::nranks != 1` and `Mesh::multilevel` in the
    production `relativistic_hc` constructor until Phase 8 qualification;
  - reject `B_profile != uniform` in the Phase-4a prescribed-test path;
  - preserve the accepted Phase-3 gather witnesses through one explicit
    `problem/relativistic_gather_diagnostic_only = true` exception restricted
    to `pgen_name = cr_relativistic_runtime_gather_test` with `time/nlim = 0`;
    this exception may exercise MPI and refinement construction but can never
    authorize a particle push;
  - keep the original `P4A-01` through `P4A-16` criteria intact and add
    explicitly labeled `P4A-R*` reviewer-driven gate-close criteria rather
    than rewriting history;
  - retain production expected-failure transcripts for initialization,
    electric-range, magnetic-range, and output-shadow-cap rejection paths;
  - bind the corrected sources, analyzer, criteria, logs, runtime case
    transcripts, and control snapshot through a direct SHA-256 manifest after
    rebound reviews complete.
- Alternatives rejected:
  - broadening Phase 4a to MPI, SMR, and AMR qualification would collapse the
    staged migration boundary and skip the guide's dedicated Phase 8 review;
  - silently treating nonuniform `B_profile` as a live-grid poison would leave
    a misleading user-facing input contract;
  - relabeling reviewer-driven criteria as preregistered would make the
    evidence history dishonest.
- Revisit trigger: corrected build, parser, runtime-oracle, MPI rejection,
  retained-regression, direct-binding, and independent rebound-review results.

## RG-005 Rebound Candidate And CP-3 Minimal Serial Evidence

- Date: 2026-05-30
- Status: rebound candidate `PROCEED`; final test-adversary seal recheck pending.
- Corrected parser evidence:
  - `92 passed`;
  - production `relativistic_hc` rejects MPI, SMR, and AMR;
  - production `relativistic_hc` rejects nonuniform `B_profile`;
  - the Phase-3 diagnostic-only exception rejects a normal pgen, rejects
    `time/nlim = 1` under the registered gather pgen, and remains rejected by
    legacy pushers.
- Corrected prescribed-test production evidence:
  - helper-only serial release and debug runs pass;
  - all emitted runtime metrics are self-enforced against the complete
    `38`-criterion JSON surface;
  - reviewer-driven `P4A-R*` additions bind non-unit-`c_model` magnetic
    rotation, high-`gamma` phase, independent negative-`sigma` momentum and
    trajectory, four production failure paths, separate convergence-component
    slopes, finest-step component errors, independently reconstructed
    `gamma`, and velocity-shadow consistency;
  - direct-work closure remains `2.168404344971009e-19` in pure `cE`,
    `2.168404344971009e-19` for non-unit `c_model`, and
    `6.776263578034403e-21` in crossed fields;
  - production convergence slope remains `1.9996235673293787`;
  - legacy-`IPM` independence and poisoned-live-grid isolation errors remain
    `0.0`.
- Corrected fail-closed evidence:
  - two-rank production launch rejects with the explicit serial-uniform-level
    diagnostic on both ranks;
  - retained production subprocesses reject unrepresentable initialization,
    electric arithmetic range, magnetic arithmetic range, and output
    velocity-shadow range.
- Replayed Phase-3 diagnostic evidence:
  - serial periodic `64` rows, SMR `512` rows, and dynamic-AMR shadow `64`
    rows pass;
  - MPI periodic ranks `1`, `2`, `4`, and `6` with an empty-particle rank,
    SMR rank `4`, and dynamic-AMR-shadow rank `4` pass;
  - traced replay transcripts retain every analyzer command and successful
    row-count result.
- Retained compatibility evidence:
  - legacy CR CPU `20 passed`;
  - legacy CR accuracy CPU `12 passed`;
  - legacy CR MPI CPU `4 passed`;
  - legacy CR accuracy MPI CPU `5 passed`;
  - style `2 passed`;
  - tracked and authored-untracked whitespace checks pass.
- Independent rebound review status:
  - physics reviewer: `PROCEED`;
  - integration-boundary reviewer: `PROCEED`;
  - test-adversary reviewer: initial rebound `HOLD` narrowed to one missing
    parser negative case, stale generic artifacts, absent traced analyzer
    transcripts, and the deferred direct SHA-256 seal; all four corrections
    are now implemented and awaiting final seal recheck.
- Reflection:
  - the Phase-4a prescribed-field pusher merits integration as a minimal
    serial checkpoint if and only if the direct evidence manifest verifies and
    the final test-adversary seal recheck returns `PROCEED`;
  - solver-coupled ideal-MHD execution remains closed until experimental
    Phase 4b records its own temporal-model decision and separate evidence.

## RG-005 Final: Prescribed-Test Pusher Merits Minimal-Serial Integration

- Date: 2026-05-30
- Status: `PROCEED`.
- Accepted checkpoint: `CP-3 Minimal Serial`.
- Direct seal:
  - `evidence/phase4a_evidence_sha256.txt`;
  - path-sorted manifest covers corrected sources, criteria, harnesses, decks,
    analyzer, generic regressions, parser and style logs, MPI rejection,
    nested prescribed-runtime cases, nested Phase-3 diagnostic replay
    artifacts, whitespace report, and control snapshot;
  - independent test-adversary replay verified every recorded digest.
- Independent rebound verdicts:
  - physics reviewer: `PROCEED`;
  - integration-boundary reviewer: `PROCEED`;
  - test-adversary seal reviewer: `PROCEED`.
- Accepted evidence summary:
  - helper-only release and debug HC harnesses pass;
  - prescribed production analyzer passes all `38/38` bound criteria;
  - parser contract passes `92` cases;
  - explicit two-rank production probe rejects before execution;
  - Phase-3 diagnostic-only serial and MPI replay ladder passes with traced
    analyzer commands;
  - retained legacy CR CPU `20`, accuracy CPU `12`, legacy MPI CPU `4`, MPI
    accuracy CPU `5`, and style `2` suites pass;
  - tracked and authored-untracked whitespace checks pass.
- Scope classification:
  - accepted: mechanically isolated serial uniform-level prescribed-field HC
    execution and the named `nlim = 0` Phase-3 gather-diagnostic exception;
  - still closed: solver-coupled ideal MHD, positive-cycle MPI, SMR, AMR,
    acceleration-aware subcycling, outer-timestep refresh, restart, outputs,
    GPU qualification, and scientific production-readiness claims.
- Reflection decision:
  - keep Higuera-Cary as the selected pusher;
  - do not reopen `DR-003`;
  - open experimental Phase 4b only as a separate grid-frozen-at-`t^n`
    solver-coupling increment with its own preregistered evidence and reviewers.

## DR-028: Experimental Solver-Coupled Frozen-`t^n` Ideal-MHD Path

- Date: 2026-05-30
- Status: selected for bounded Phase-4b implementation after `RG-005 PROCEED`
  and `CP-3 Minimal Serial`.
- Question: what is the smallest honest solver-coupled acceleration increment
  that can follow the mechanically isolated prescribed-field checkpoint?
- Options considered:
  1. reuse CT edge EMFs or add a stage-centered MHD coupling immediately;
  2. gather the repaired Newtonian primitive velocity and cell-centered
     magnetic field from the live grid at the particle spatial midpoint,
     construct checked pointwise `cE = -u_fluid x B`, and freeze that sample
     for one Higuera-Cary update;
  3. leave all live-grid execution closed until the later subcycling and
     timestep redesigns are complete.
- Selected option: option 2.
- Reasoning:
  - the accepted Phase-3 trilinear helper already applies one stencil to
    `u_fluid` and cell-centered `B`;
  - the accepted Higuera-Cary helper already consumes a single sampled field
    tuple and commits only after candidate validation;
  - this adds the missing live ideal-MHD force without changing CT ownership,
    MHD RK stages, restart bytes, outputs, or relativistic subcycling;
  - explicit `SaveMHDState -> Push` task dependency documents the chosen
    repaired-grid `t^n` scheduling boundary instead of relying on insertion
    order;
  - dynamic evolving-field accuracy is expected to be first order globally
    under this temporal model, even though the spatial drift-kick-drift update
    remains midpoint based.
- Runtime contract:
  - select
    `<particles>/relativistic_field_source = mhd_ideal` and
    `<particles>/relativistic_temporal_sampling = frozen_tn`;
  - require dynamic, Newtonian, three-dimensional, strictly periodic,
    serial, uniform-level, passive ideal-gas MHD with trilinear gather;
  - reject diffusion, source blocks, CT-EMF reuse aliases, restart, outputs,
    and relativistic subcycling until their owning phases qualify them;
  - preserve `prescribed_test` as a mechanically separate kinematic harness.
- Atomicity requirement:
  - gather and validate live `u_fluid` and `B`;
  - construct `cE` through checked arithmetic;
  - compute the full Higuera-Cary and drift-kick-drift candidate;
  - commit sampled `B`, sampled `cE`, authoritative `w`, velocity shadows,
    explicit work, positions, displacement, and `IPDB` only after every check
    succeeds.
- Alternatives retained:
  - stage-coupled or time-centered reconstruction remains a separate widening
    if later physics requirements demand second-order evolving-field accuracy;
  - CT EMF reuse remains excluded because edge-centered solver EMFs are not
    the selected particle-facing pointwise ideal field;
  - production MPI, SMR, AMR, acceleration-aware subcycling, and outer
    timestep refresh remain deferred to their dedicated phases.
- Inflection point:
  - stop and reassess before `CP-4 Solver Coupling` if the live-grid poison
    test does not prove overwrite, the cancellation invariant fails, or a
    manufactured evolving-field ladder does not show the preregistered
    frozen-`t^n` first-order trend.

### DR-028 Independent Scheduling-Audit Corrections

- Date: 2026-05-30
- Status: accepted corrections before Phase-4b qualification.
- Findings:
  1. the first parser patch introduced
     `<particles>/relativistic_mhd_temporal_mode`, which drifted from the
     preregistered `<particles>/relativistic_temporal_sampling` spelling;
  2. rejecting `<mhd_srcterms>` alone did not reject
     `<problem>/user_srcs = true`, which executes through the MHD source task;
  3. the particle constructor rejects legacy particle-only restart input but
     cannot observe command-line full-mesh `-r` restart selection;
  4. the live gather indexed `PGID - gids` without a local-range guard.
- Corrections:
  - use only `<particles>/relativistic_temporal_sampling = frozen_tn`;
  - reject the unregistered temporal spelling rather than treating it as an
    alias;
  - reject user source enrollment and orphan source blocks in the bounded
    coupled contract;
  - reopen `src/main.cpp` narrowly to reject full-mesh restart selection for
    `relativistic_hc` before mesh construction;
  - fail closed before a live gather if a particle does not map to a local
    MeshBlock index.
- Alternative retained:
  - restart support remains Phase-6 work; this correction rejects unsupported
    state restoration and does not introduce serialization behavior.

## RG-005 Phase-4b Candidate: Does The Experimental Frozen-`t^n` Coupling Merit CP-4?

- Date: 2026-05-30
- Status: candidate `PROCEED`; fresh physics, test-adversary, and
  integration-boundary done-claim reviews in progress.
- Scope reviewed:
  - explicit `mhd_ideal` field source with one accepted
    `relativistic_temporal_sampling = frozen_tn` selector;
  - live primitive-velocity and cell-centered-`B` midpoint gather;
  - checked pointwise `cE = -u_fluid x B`;
  - explicit `SaveMHDState -> Particles::Push` scheduling edge;
  - atomic sampled-field and accepted-state commit;
  - fail-closed serial uniform-level, three-dimensional, strictly periodic,
    Newtonian ideal-MHD execution envelope.
- Preregistered coupled-runtime evidence:
  - all `20/20` Phase-4b criteria pass;
  - stale sampled-field poison overwrite error is `0.0`;
  - prescribed-versus-coupled parity error is `0.0`;
  - `cE dot B` error is `0.0`;
  - fluid-velocity cancellation momentum and trajectory errors are
    `2.7755575615628914e-17` and `3.602423855278284e-18`;
  - acceleration and deceleration direct-work closure errors are
    `1.6940658945086007e-21` and `6.071532165918825e-17`;
  - spatial-midpoint field error is `3.469446951953614e-17`, with a
    `0.003922322702763681` separation from stale old-position sampling;
  - the dynamic-MHD frozen-`t^n` ladder measures slope
    `0.9872005890039331` and finest error `4.726247195889718e-06`.
- Temporal-qualification honesty:
  - the accepted ladder is a periodic dynamic-MHD self-convergence test
    against a finer reference run;
  - it is not a closed-form manufactured evolving-field solution;
  - therefore it supports only the experimental first-order frozen-`t^n`
    temporal classification and does not justify stage-centered, nonlinear
    production-science, or broad dynamic-MHD accuracy claims.
- Parser and retained prescribed evidence:
  - corrected parser matrix passes `110` cases;
  - prescribed Phase-4a production analyzer remains `38/38`.
- Review-driven corrections already applied:
  - normalized the temporal-selector spelling;
  - rejected obsolete aliases, user source enrollment, orphan source blocks,
    and full-mesh restart;
  - added a local MeshBlock-index guard before live gather;
  - added a defensive standalone-MHD task dependency guard.
- Scope still closed:
  - relativistic subcycling and outer-timestep refresh;
  - restart serialization and outputs;
  - positive-cycle MPI, SMR, AMR, GPU, and scientific production-readiness
    claims.

## DR-029: Phase-4b Rebound After Independent Done-Claim Review

- Date: 2026-05-30
- Status: corrections implemented; rebound evidence replay and fresh independent
  review pending.
- Trigger: the first Phase-4b done-claim review returned `HOLD`.
- Findings:
  1. the coupled midpoint gather could index beyond the allocated ghost stencil
     when the particle midpoint left the supported local range;
  2. same-level multi-MeshBlock ownership and periodic-wrap behavior had not
     been qualified for solver-coupled live gathers;
  3. the coupled timestep still inherited a legacy cap and had not qualified
     arbitrary `c_model`;
  4. the original evolving-field ladder used fine-reference self convergence
     rather than an independently integrated manufactured field;
  5. the stale-field poison criterion did not prove at runtime that poisoned
     per-particle fields existed before the live gather overwrote them;
  6. explicit negative probes were missing for local-GID, superluminal-fluid,
     and stencil-range rejection;
  7. two README passages still described `prescribed_test` as the only
     relativistic allocation and output-rejection path;
  8. the first rebound stencil guard still constructed integer stencil
     indices before rejecting an extreme finite midpoint coordinate;
  9. the first rebound one-MeshBlock fence unintentionally narrowed the
     retained `prescribed_test` Phase-4a surface;
  10. the manufactured source exception trusted an input-controlled
      `pgen_name` instead of the compiled test-harness identity.
  11. first-order manufactured convergence did not distinguish the intended
      frozen-`t^n` left-end sample from an incorrect right-end sample.
- Options considered:
  1. qualify migration, wrap, arbitrary-`c_model`, and acceleration-aware
     timesteps immediately inside Phase 4b;
  2. narrow Phase 4b to one MeshBlock and coupled `c_model = 1.0`, add a
     pre-gather stencil guard, and defer those widenings to their owning
     phases;
  3. close solver coupling again and return to the Phase-4a prescribed-only
     checkpoint.
- Selected option: option 2.
- Reasoning:
  - one MeshBlock removes the unqualified same-level ownership transition
    without discarding the live-grid force path;
  - coupled `c_model = 1.0` keeps the opening normalized until Phase 5 replaces
    the inherited timestep model with acceleration-aware bounds;
  - an explicit stencil-range guard fails closed before dereferencing live MHD
    arrays;
  - this preserves Phase 4b as a bounded solver-coupling checkpoint rather than
    importing Phase-5 and Phase-8 obligations prematurely.
- Implemented corrections:
  - require exactly one MeshBlock for positive-cycle solver-coupled
    `relativistic_field_source = mhd_ideal` execution until Phase-8 migration
    qualification;
  - require coupled `<particles>/c_model = 1.0` until the Phase-5 timestep
    redesign;
  - validate local `PGID - gids` and the complete cell-centered trilinear
    coordinate range before constructing integer stencil indices or gathering
    live-grid fields;
  - scope the exactly-one-MeshBlock execution fence to `mhd_ideal`, preserving
    the retained same-level multi-MeshBlock `prescribed_test` harness;
  - bind the manufactured source exception to the compiled
    `unit_tests/cr_relativistic_coupled_runtime_test` harness as well as its
    runtime pgen selector;
  - add a test-generator-only manufactured MHD source that evolves uniform
    fluid velocity in time, while the analyzer compares production results
    against an independent RK4 integration;
  - emit and consume a runtime pre-push witness showing poisoned stored `B` and
    `cE` values before live-grid overwrite;
  - emit final-live gathered fields under explicit `final_live_*` names and
    bind `final_live_cE_dot_b`;
  - bind a one-step manufactured orientation witness: the committed particle
    `cE` and momentum must match the pre-source `t^n` field, the final live
    field must match the separately evolved right-end state, and a deliberate
    right-end oracle must be rejected;
  - add runtime probes for status `207`, `201`, ordinary status `208`, and
    positive and negative extreme finite-coordinate status `208` preflights;
  - correct the README allocation and output wording.
- Rebound criteria:
  - the reviewer-rebound surface now contains `34` criteria binding `32`
    unique metrics;
  - the criteria file is relabeled explicitly as a post-HOLD reviewer-rebound
    freeze before corrected replay rather than as the original pre-harness
    preregistration;
  - the original periodic self-convergence ladder remains as a compatibility
    witness;
  - the independent manufactured ladder is the controlling evolving-field
    temporal witness.
- Alternatives retained:
  - acceleration-aware subcycling and timestep refresh remain Phase-5 work;
  - restart schema remains Phase-6 work;
  - explicit outputs remain Phase-7 work;
  - same-level migration, periodic wrap, MPI, SMR, and AMR remain Phase-8
    work.
- Inflection point:
  - do not accept `CP-4 Solver Coupling` unless a fresh independent reviewer
  replay verifies the narrowed parser contract, the manufactured ladder, all
  five runtime negative probes, the executable restart fence, and the
  direct evidence seal.

## CP-4 Solver Coupling: Experimental Frozen-`t^n` Ideal-MHD Opening

- Date: 2026-05-30
- Verdict: `PROCEED`
- Accepted scope:
  - one-way passive Newtonian ideal-MHD sampling only;
  - explicit `relativistic_field_source = mhd_ideal`;
  - explicit `relativistic_temporal_sampling = frozen_tn`;
  - live primitive velocity and cell-centered magnetic field gathered at the
    particle spatial midpoint from the repaired `t^n` state;
  - pointwise normalized `cE = -u_fluid x B`;
  - serial, uniform-level, exactly-one-MeshBlock, 3-D strictly periodic
    execution with normalized `c_model = 1.0`.
- Accepted evidence:
  - `34/34` rebound criteria bind `32` unique metrics and pass;
  - the one-step manufactured orientation witness binds committed sampled
    `cE` to `t^n` with error `2.4747374745180785e-18`, committed momentum with
    error `0.0`, final-live right-end `cE` with error
    `1.0842021724855044e-19`, and a nontrivial left-versus-right separation of
    `0.0035902646142032483`;
  - a deliberate right-end oracle is rejected with sampled-field mismatch
    `3.590265e-03`;
  - the independent evolving-field RK4 ladder measures slope
    `0.9902663432088565` and finest error `1.0848032296554225e-05`;
  - multi-cycle accumulated-work closure error is
    `1.4405658740543337e-16`;
  - release and Kokkos bounds-check coupled analyzers pass;
  - parser contract replay passes `114` cases;
  - retained prescribed Phase-4a analyzer passes;
  - fixed-energy compatibility regressions pass CPU `20 + 12`, MPI CPU
    `4 + 5`, style `2`, and whitespace checks;
  - the executable two-rank coupled fence and full-mesh-restart fence reject
    unsupported paths;
  - status probes reject invalid GID `207`, superluminal fluid `201`,
    ordinary stencil range `208`, and positive and negative extreme finite
    coordinates `208`;
  - the final path-sorted Phase-4b direct evidence seal authenticates `160`
    files and verifies cleanly.
- Independent review disposition:
  - initial three-reviewer `HOLD` exposed the unsafe integer conversion,
    overbroad one-MeshBlock fence, spoofable manufactured-source exception,
    evidence-provenance gaps, and missing cumulative controls;
  - post-rebound review exposed replay-packaging gaps and the temporal-
    orientation oracle gap;
  - the final narrow reviewer returned `PROCEED` after verifying the sealed
    package and orientation-sensitive mutation rejection.
- Deferred residual risks:
  - the shared trilinear helper still depends on callers to preflight
    representable coordinates;
  - relativistic subcycling and outer timestep refresh remain Phase-5 work;
  - restart remains Phase-6 work;
  - diagnostics remain Phase-7 work;
  - same-level migration, periodic wrap, MPI, SMR, and AMR remain Phase-8
    work;
  - GPU and public documentation overlay remain Phase-9 work.
- Next permitted edit set:
  - Phase 5 only: `src/particles/particles_pushers.cpp`,
    `src/particles/particles_tasks.cpp`, `src/pgen/part_random.cpp`,
    timestep-focused tests and evidence, and conditionally
    `src/mesh/mesh.cpp` after an explicit ownership review.

## DR-030: Phase-5 Shared Bound Ownership And Reopened Write Manifest

- Date: 2026-05-30
- Status: accepted before Phase-5 production edits
- Question: where should the accepted `DR-011` global subcycle schedule and
  `DR-016` dirty-aware outer timestep refresh live, and which previously
  closed files must reopen to implement the lifecycle contract without
  silently widening physics scope?
- Selected implementation:
  - keep the schedule and shared refresh implementation in
    `src/particles/particles_pushers.cpp`;
  - retain one pessimistic MPI-global subcycle count for the full outer step
    from the accepted grid-envelope formula; do not introduce per-rank,
    per-particle, or post-kick adaptive collective schedules;
  - spatially re-gather the frozen-`t^n` fields at every relativistic particle
    substep midpoint through the already accepted coupled pusher;
  - use `src/particles/particles_tasks.cpp` to refresh after the completed MHD
    time integrator and use `src/mesh/mesh.cpp` as the defensive final consumer
    check immediately before the particle `dtnew` value is consumed;
  - mark the bound dirty at the end of particle push and after particle AMR
    remap, but refresh AMR-dependent fields only in
    `src/mesh/mesh_refinement.cpp` after face-field repair and primitive
    reconstruction;
  - reopen `src/particles/particles.cpp` only for new-mode parser validation,
    dirty-state initialization, and the too-early AMR-remap dirty marker;
  - expose the last selected subcycle count as a diagnostic harness witness,
    without changing legacy output meanings.
  - compute coupled envelopes over unique active leaf cells.  The qualified
    mode is strictly periodic, so boundary-filled ghost cells are duplicates;
    including pre-fill harness ghost extrapolations would make initialization
    order affect the physical envelope.
- Alternatives rejected for this phase:
  - refreshing inside `RemapAfterAMR()` is too early because repaired magnetic
    fields and reconstructed primitives are not available yet;
  - relying only on an after-integrator task is insufficient because
    initialization and future AMR callers need a defensive consumer check;
  - adding fractional momentum or kinetic-energy constraints is deferred as
    accepted in `DR-011`: both need an independently reviewed nonsingular
    scale near rest.
- Phase-5 reopened write manifest:
  - `src/particles/particles.hpp`;
  - `src/particles/particles.cpp`;
  - `src/particles/particles_pushers.cpp`;
  - `src/particles/particles_tasks.cpp`;
  - `src/mesh/mesh.cpp`;
  - `src/mesh/mesh_refinement.cpp`;
  - `src/pgen/part_random.cpp`;
  - `src/pgen/unit_tests/cr_relativistic_coupled_runtime_test.cpp`;
  - timestep-focused scripts, criteria, evidence, and this append-only ledger.
- Inflection point:
  - stop before `RG-006` if current-field consumption cannot be demonstrated
    at initialization, after MHD evolution, and after the AMR repair hook, or
    if the global schedule requires cap clipping rather than a clear abort.

## DR-031: Phase-5 Reviewer Rebound And Test-Harness Isolation

- Date: 2026-05-30
- Status: accepted after independent HOLD findings and corrected replay
- Independent HOLD findings:
  - exact-integer ratios in the inherited `StepsForRatio()` helper used `<=`
    and requested one excess substep;
  - an initial draft refreshed `dtnew` inside `Push()` before the final
    task-list exchange;
  - the coupled test pgen refreshed before the driver's initial boundary fill,
    masking the required defensive-consumer recomputation;
  - broad activation tests did not bind every isolated component count;
  - the first overflow probe tested integer saturation, not checked-product
    overflow;
  - the standalone prescribed-field kernel oracle inherited the new outer
    refresh and lost its historical controlled timestep ladder.
- Corrected choices:
  - use a mathematically exact ceiling helper and saturate only unreasonable
    nonfinite or integer-range requests;
  - mark dirty at the end of `Push()` without recomputing there;
  - let the coupled pgen mark dirty only, so `Driver::Initialize()` fills
    boundaries and reconstructs primitives before `Mesh::NewTimeStep()`
    defensively refreshes and consumes the initial bound;
  - bind exact crossing, gyro, electric-kick, gamma-floor, and total schedule
    counts in the Phase-5 analyzer;
  - require a true finite-input electric-kick overflow abort before integer
    conversion;
  - preserve the standalone prescribed-field kernel oracle through a narrowly
    compiled `unit_tests/cr_relativistic_pusher_runtime_test`-only timestep
    override.  This bypass is test-only and does not apply to coupled or
    production execution.
- Mutation controls:
  - reject eight weakened-source variants: exact-ceil regression, removed cap
    preflight, lost checked multiplication, lost gamma floor, lost midpoint
    regather placement, lost post-integrator refresh, premature initialization
    refresh, and refresh inside `Push()`.
- Alternatives retained:
  - fractional-momentum and fractional-kinetic-energy controls remain rejected
    pending an independently reviewed nonsingular near-rest scale;
  - per-rank and per-particle adaptive schedules remain deferred because the
    accepted collective exchange requires one shared schedule;
  - actual MPI, same-level migration, SMR, and AMR execution qualification
    remain Phase-8 work.

## RG-006: Is The Temporal Model Still Defensible?

- Date: 2026-05-30
- Verdict: `PROCEED` for the bounded passive first-order model
- Reflection:
  - Phase-5 substeps spatially re-gather the accepted frozen-`t^n` primitive
    velocity and cell-centered magnetic field at every particle midpoint;
  - the schedule remains one pessimistic full-outer-step MPI-global envelope,
    so electric acceleration cannot silently weaken the gyro bound during the
    outer step;
  - the next-step `dtnew` cache is refreshed only after current fields exist,
    with lifecycle hooks after the MHD integrator, after AMR repair and
    primitive reconstruction, and defensively before mesh consumption;
  - this does not widen the Phase-4b temporal claim beyond first-order
    evolving-field accuracy.
- Revisit trigger:
  - if target science requires stage-centered evolving fields or materially
    tighter schedules, stop and design a separately reviewed stage-coupled
    follow-up rather than mixing stage fields into this passive path.

## CP-5 Acceleration-Aware Subcycling And Outer Refresh

- Date: 2026-05-30
- Verdict: `PROCEED`
- Accepted scope:
  - retain the Phase-4b passive, one-MeshBlock, serial, uniform-level,
    normalized ideal-MHD execution fence;
  - add one pessimistic full-outer-step acceleration-aware subcycle schedule
    from model-speed crossing, MeshBlock-bound, relativistic gyro-angle, and
    electric-kick constraints;
  - abort before first drift when the global request exceeds the configured
    cap; never clip relativistic requests;
  - spatially re-gather frozen-`t^n` coupled fields at every substep midpoint;
  - refresh dirty outer bounds after current fields exist and defensively at
    mesh consumption.
- Accepted evidence:
  - release and Kokkos bounds-check Phase-5 analyzers pass all `36/36`
    criteria, including `12` lifecycle source-contract assertions;
  - all `8/8` deliberate weakening mutations are rejected;
  - retained prescribed Phase-4a and coupled Phase-4b analyzers pass;
  - parser contract replay passes `117` cases;
  - fixed-energy compatibility regressions pass CPU `20 + 12`, MPI CPU
    `4 + 5`, style `2`, and whitespace checks;
  - the explicit two-rank relativistic coupled execution fence still rejects
    MPI before Phase-8 qualification.
- Independent review disposition:
  - the architecture review corrected ownership, rejected early remap refresh,
    preserved the `c_model == 1` coupled fence, and prevented legacy parser
    tightening;
  - the test-adversary HOLD exposed the inherited exact-ceil bug, masked
    initialization consumption, incomplete isolated-bound checks, weak
    overflow probe, and missing mutation controls;
  - the lifecycle HOLD removed the premature refresh inside `Push()`;
  - both fresh post-rebound reviewers returned `PROCEED`.
- Deferred residual risks:
  - same-level migration, periodic wrap, MPI, SMR, and AMR execution remain
    Phase-8 work;
  - restart schema remains Phase-6 work;
  - explicit production diagnostics remain Phase-7 work;
  - GPU and public documentation overlay remain Phase-9 work.
- Next permitted edit set:
  - Phase 6 only: versioned relativistic restart writer, reader, inspector,
    restart-focused tests and evidence, plus append-only ledger updates.

## DR-032: Phase-6 Paired Typed-V2 Restart Boundary

- Date: 2026-05-30
- Status: accepted candidate pending final restart-specialist replay
- Question: how should relativistic particle state cross a restart boundary
  without weakening legacy pushers, silently truncating identity, or allowing
  an unpaired particle shard to masquerade as a complete checkpoint?
- Selected implementation:
  - preserve legacy v1 seventeen-`Real` particle restart payloads for legacy
    pushers through one centralized decoder;
  - require `relativistic_hc` to use a paired native mesh `.rst`, mesh
    `.rst.rmeta` witness, rank shard, and typed-v2 manifest;
  - encode every relativistic particle identity as three little-endian
    `int32` values: `PGID`, `PTAG`, and `PSP`;
  - encode all authoritative and compatibility real state as twenty-two
    little-endian `float64` values in the exact `ParticlesIndex` order:
    `IPX`, `IPVX`, `IPY`, `IPVY`, `IPZ`, `IPVZ`, `IPM`, `IPBX`, `IPBY`,
    `IPBZ`, `IPDX`, `IPDY`, `IPDZ`, `IPDB`, `IPWX`, `IPWY`, `IPWZ`,
    `IPCEX`, `IPCEY`, `IPCEZ`, `IPWORK`, and `IPALPHA`;
  - keep `w = gamma v` authoritative, keep physical velocity as an invariant-
    checked compatibility shadow, and reconstruct `gamma` diagnostically
    rather than serializing another authoritative field;
  - publish a fixed-width `144`-byte shard header with magic, schema,
    endian marker, record widths, saved rank topology, local and global
    counts, mesh cycle, mesh time, mesh timestep, particle timestep,
    checkpoint ID, payload bytes, payload checksum, header checksum, manifest
    requirement, and semantic configuration fingerprint;
  - bind the mesh witness and particle manifest to the native mesh checkpoint
    byte count, checksum, rank-to-MeshBlock topology hash, exact five-digit
    checkpoint number, cycle, time, timestep, and deterministic checkpoint ID;
  - deterministically reject MPI rank-count change for typed-v2 restart.
- Alternatives rejected:
  - reusing legacy seventeen-`Real` records would omit relativistic state and
    continue routing integer identity through floating-point storage;
  - silently reconstructing relativistic momentum from legacy velocity is an
    ambiguous model migration and is not labeled backward compatible;
  - implementing rank-count redistribution inside this phase would mix a
    restart-schema change with unqualified Phase-8 migration work.

## DR-033: Phase-6 Restart-Specialist Rebound

- Date: 2026-05-30
- Status: corrected candidate replayed after two independent HOLD rounds
- Initial HOLD findings:
  - the shard could be renamed before paired mesh-witness preflight;
  - paired output blocks lacked one coordinator for cadence and counter
    equality;
  - restart state semantic validation depended too heavily on optional debug
    checks;
  - the semantic fingerprint omitted accepted subcycle controls;
  - five-digit writer bounds and strict parser tests were incomplete.
- Second HOLD findings:
  - mesh-witness validation occurred after native mesh restart consumption;
  - the C++ reader did not bind particle and mesh checkpoint numbers exactly;
  - the Python inspector accepted checkpoint identities and filename widths
    that C++ rejected;
  - semantic particle validation could occur after native mesh output had
    already started overwriting a same-number checkpoint;
  - the runtime qualifier did not mutate native mesh bytes, fabricate a
    checkpoint ID, exercise overwidth shard filenames, or compare regenerated
    native mesh state payloads.
- Third HOLD findings:
  - the runtime qualifier recorded the criteria SHA-256 digest but did not
    enforce an immutable expected digest or the complete registered criterion
    ID sequence;
  - Python `int()` accepted underscore-separated numeric spellings that the
    strict C++ parser rejected, allowing the offline inspector to certify a
    bundle that Athena refused to restore.
- Corrected choices:
  - preflight any present mesh witness before native restart parameter reads
    and require the witness again for paired relativistic restoration before
    mesh-tree construction;
  - validate relativistic particle semantics in mesh restart
    `LoadOutputData()` before native mesh output staging and write;
  - parse particle and mesh file numbers strictly and require exact equality;
  - make the Python inspector recompute the deterministic checkpoint ID and
    enforce the same typed-v2 five-digit filename contract;
  - compare native mesh state bytes after the `<par_end>` parameter terminator,
    because restart-only input flags intentionally change the native
    parameter dump while the serialized mesh state remains byte-identical;
  - bind the corrected runtime replay to a frozen criteria file with schema,
    status, exact ordered criterion IDs, archived SHA-256 digest, and a
    hard-coded expected SHA-256 assertion;
  - apply a strict decimal lexical grammar before Python typed-v2 integer and
    floating-point conversion so the inspector cannot bless runtime-rejected
    numeric spellings.

## RG-007: Is The Migration Contract Honest?

- Date: 2026-05-30
- Verdict: `PROCEED` candidate pending final restart-specialist replay
- Reflection:
  - legacy v1 remains supported for legacy pushers with centralized validation
    of header length, record count, finite state, optional metadata, and
    checksum;
  - `relativistic_hc` does not claim legacy-v1 conversion because reconstructing
    the authoritative momentum, sampled `cE`, accumulated work, signed
    coupling, or semantic fingerprint would be lossy or ambiguous;
  - typed-v2 restart requires a complete paired checkpoint and rejects stale,
    missing, mismatched, corrupted, future-version, truncated, semantically
    incompatible, and rank-count-changed inputs;
  - rank-count redistribution remains explicit Phase-8 follow-up work rather
    than an undocumented partial feature.
- Revisit trigger:
  - if Phase 8 opens rank-count redistribution or legacy-to-relativistic
    conversion, stop and write a new reviewed migration design before coding.

## CP-6 Versioned Relativistic Restart And Inspection Contract

- Date: 2026-05-30
- Verdict: `PROCEED`
- Accepted scope:
  - preserve validated legacy v1 particle restart behavior for legacy pushers;
  - require `relativistic_hc` to use paired native mesh restart, mesh witness,
    typed-v2 rank shard, and typed-v2 manifest artifacts;
  - serialize typed identity and the full twenty-two-real relativistic state
    without treating derived `gamma` as authoritative;
  - reject typed-v2 MPI rank-count changes deterministically pending Phase 8;
  - require strict C++ and Python parsing, checksums, topology witness,
    checkpoint-number equality, deterministic checkpoint identity, and
    semantic configuration equality.
- Accepted evidence:
  - release serial and MPI CPU builds pass;
  - paired continuation regenerates byte-identical typed-v2 particle shards
    and byte-identical native mesh state payloads after the native parameter
    dump terminator;
  - the runtime qualifier passes `21/21` frozen criteria, `13/13` copied-
    artifact runtime corruption controls, and `1/1` inspector numeric-grammar
    control;
  - direct typed-v2 inspector replay passes;
  - parser contract replay passes `120` cases;
  - fixed-energy compatibility regressions pass CPU `20 + 12` and MPI CPU
    `4 + 5`;
  - the two-rank relativistic execution fence still rejects MPI before
    Phase-8 qualification;
  - style and whitespace gates pass;
  - a tampered criteria copy is rejected by immutable SHA-256 assertion;
  - the path-sorted Phase-6 evidence seal verifies every archived artifact.
- Independent review disposition:
  - three restart-specialist HOLD rounds exposed shard-publication ordering,
    paired-output coordination, semantic-fingerprint coverage, mesh-witness
    preflight ordering, mesh-particle number binding, Python/C++ parser drift,
    missing runtime controls, mutable criteria acceptance, and stale evidence
    sealing;
  - each finding was corrected and mutation-replayed;
  - the final narrow restart-specialist audit returned `PROCEED` with no
    remaining Phase-6 blocker.
- Deferred residual risks:
  - relativistic diagnostic semantics remain Phase-7 work;
  - MPI execution, same-level migration, SMR, AMR migration, restart after
    migration, and any rank-count redistribution remain Phase-8 work;
  - GPU, performance, public documentation overlay, and final handoff remain
    Phase-9 work.
- Next permitted edit set:
  - Phase 7 only: particle-facing output implementations, output metadata,
    inspector output validation, diagnostics-focused tests and evidence,
    `src/particles/README.md`, and append-only ledger updates;
  - keep migration infrastructure, MPI widening, AMR algorithms, GPU
    optimization, and public documentation overlay closed.

## DR-034: Phase-7 Explicit Diagnostic Allowlist And Semantic Names

- Date: 2026-05-30
- Status: accepted before Phase-7 production edits
- Question: which particle-facing outputs can be interpreted honestly for
  `relativistic_hc`, and which legacy names must remain isolated rather than
  silently changing meaning?
- Selected implementation:
  - allow `ppd`, `df`, `dxh`, `drh`, `dparh`, `pmom`, `pspec`, `pspec2`,
    `psamp`, and paired `prst/rst` only;
  - preserve `ppd` species and position payload meaning;
  - preserve `df` pitch-angle cosine from physical velocity shadow and sampled
    magnetic field;
  - preserve `dxh` physical component displacement, `drh` physical
    displacement norm, and define `dparh` explicitly as the accumulated
    midpoint sum of `dx dot Bhat`, not final displacement projected onto final
    sampled `B`;
  - preserve `pmom` physical-velocity transport moments with explicit schema,
    mode, and units metadata;
  - add canonical relativistic spectrum names: `speed`, `wmag`,
    `kinetic_energy_model`, `log10_kinetic_energy_model`, `mu`, and an
    explicitly named velocity-based magnetic-moment proxy;
  - add canonical relativistic joint-spectrum names: `mu_speed`, `mu_wmag`,
    `mu_kinetic_energy_model`, and physical-velocity `vpar_vperp`;
  - add canonical samples for momentum, `gamma`, model kinetic energy,
    sampled `cE`, accumulated work, signed `alpha_s`, and applicability
    diagnostics while retaining clearly named physical fields;
  - reject ambiguous relativistic aliases and defaults: legacy `pspec p`
    means speed, legacy energy means `0.5 v^2`, legacy `psamp p` means speed,
    and legacy `psamp mass` emits overloaded `IPM`, which is not physical
    mass for `relativistic_hc`;
  - reject `trk`, `pvtk`, generic `prtcl_d`, and advertised-but-unimplemented
    generic `prtcl_all` for the new mode.
- Alternatives rejected:
  - silently redefining old `p`, `E`, `energy`, or `mass` labels would make
    mixed legacy and relativistic analyses scientifically ambiguous;
  - repairing tracked-particle, VTK, and generic deposition infrastructure
    inside the diagnostics phase would mix independent correctness defects
    into the relativistic semantic opening.
- Inflection point:
  - stop before `RG-008` if a metadata-only reader cannot identify every
    retained field's model meaning and units convention without reading the
    implementation patch.

## DR-035: Phase-7 Metadata And Histogram-Safety Rebound

- Date: 2026-05-30
- Status: corrected candidate replayed; independent rebound reviews pending
- Initial metadata-only HOLD findings:
  - retained `df`, `dxh`, `drh`, and `dparh` headers did not identify bin
    geometry, reduction policy, or species count;
  - `psamp` did not state units for sampled `cE`, `gamma`, signed `alpha_s`,
    or `r_larmor_over_dx_min`;
  - the tracked qualifier exercised only `psamp`, one `pspec`, and paired
    restart, so it could not bind the full retained diagnostic dictionary;
  - the tracked runtime artifact bundle did not expose `ppd`, `df`, `dxh`,
    `drh`, `dparh`, `pmom`, `pspec2`, or paired mesh restart artifacts to the
    independent scientist.
- Fresh code-adversary HOLD findings:
  - `df`, `dxh`, `drh`, `dparh`, and `pspec` accepted zero bins; the first four
    also accepted empty or non-finite ranges before device atomics;
  - the inspector skipped any binary prefix byte equal to `#`, allowing a
    valid first little-endian histogram bin count of `35` to be misread as a
    comment;
  - scalar histogram and `pmom` inspector paths discarded emitted metadata;
  - the shared legacy histogram helper no longer checked an optional `pspec`
    spectrum normalization;
  - spectrum-derived values could reach floating-point-to-integer bin
    conversion without one generic finite-state guard.
- Corrected choices:
  - emit explicit histogram geometry and reduction metadata for every retained
    scalar histogram; emit complete sampled-field and moment-family units;
  - make the tracked diagnostics fixture produce the full retained dictionary
    and paired typed-v2 restart bundle in one run;
  - extend the frozen Phase-7 criteria from `27` to `33` checks, adding both
    applicability-ratio checks, metadata coverage, inspector round trips, and
    paired restart artifact inspection;
  - reject nonpositive bin counts, non-finite bounds, and empty histogram
    ranges before allocation for every retained histogram family;
  - reject any non-finite derived histogram value before integer bin
    conversion and device mutation;
  - parse the final fixed-width binary histogram payload from the record end,
    then accept only recognized textual metadata and layout prefixes;
  - add metadata-preserving inspector record APIs while retaining the
    list-returning compatibility wrappers;
  - restore optional legacy `pspec` normalization through
    `validate_histograms()`;
  - add a leading-count-`35` payload regression, metadata mismatch mutation,
    and invalid histogram-geometry constructor tests.
- Alternatives rejected:
  - keeping metadata only in C++ comments would not satisfy the scientist-facing
    contract;
  - widening the parser's comment skipper with more byte heuristics would keep
    binary/text ambiguity;
  - limiting geometry validation to relativistic mode would leave the shared
    output constructors unsafe for legacy callers.
- Inflection point:
  - keep Phase 8 closed until a fresh metadata-only scientist and fresh
    adversarial source reviewer both return `PROCEED` on the replayed bundle.

## DR-036: Phase-7 Second Rebound For Offline Semantics And Safe Binning

- Date: 2026-05-30
- Status: corrected candidate replayed; second independent rebound pending
- Metadata-only scientist HOLD findings:
  - the offline inspector normalized forbidden relativistic sample aliases and
    accepted forbidden relativistic spectrum and joint-spectrum aliases;
  - scientist-facing documentation mixed legacy particle-only restart wording
    with the paired relativistic typed-v2 restart contract;
  - the qualifier checked artifact existence for intermediate typed-v2
    checkpoints but decoded only the latest checkpoint;
  - general non-unit `c_model` formulas and work quadrature remained implicit;
  - the archived full-format report lacked replay provenance.
- Source-adversary HOLD findings:
  - finite out-of-range histogram values were converted to integer bins before
    clamping, so an extreme normalized coordinate could reach undefined
    floating-point-to-integer conversion;
  - `df` and `pspec2` zero-field fallbacks could hide non-finite dependencies;
  - the binary-prefix regression covered `#` but not whitespace-leading
    payloads;
  - Phase-7 source hashes and evidence seal had not yet been archived.
- Corrected choices:
  - reject relativistic legacy aliases in offline `psamp`, `pspec`, and
    `pspec2` parsing, and reject sample metadata that appears after columns;
  - describe legacy and relativistic restart paths separately in the README,
    including the typed-v2 shard, manifest, mesh restart, and witness bundle;
  - document `gamma = sqrt(1 + |w|^2 / c_model^2)`, `v = w / gamma`, model
    kinetic energy, and the direct substep work quadrature;
  - add a non-unit `c_model = 2` work-closure runtime replay and extend the
    frozen criteria from `33` to `43` checks;
  - decode every typed-v2 checkpoint during retained-dictionary qualification;
  - make offline all-checkpoint duplicate detection checkpoint-aware so
    temporal history remains inspectable without weakening same-checkpoint
    duplicate rejection;
  - branch on histogram range before integer conversion and reject
    non-finite dependencies before legitimate zero-field fallbacks;
  - reject non-finite histogram spans before allocation;
  - add forbidden-alias, metadata-order, whitespace-leading binary payload,
    and overflow-range controls;
  - include replay root, binary, input, criteria path, criteria digest,
    metrics path, generator label, decoded shard count, and decoded format
    counts in the full-format archived report.
- Alternatives rejected:
  - trusting writer-side rejection alone would allow stale or hand-edited
    artifacts to be silently reinterpreted offline;
  - post-cast clamping is not a valid overflow defense;
  - documenting only the normalized `c_model = 1` formulas would leave the
    accepted prescribed-test parameter contract underspecified.
- Inflection point:
  - keep Phase 8 source widening closed until second-rebound reviewers return
    `PROCEED` and the path-sorted Phase-7 source snapshot and evidence seal
    verify.

## DR-037: Phase-7 Spectrum Cardinality Rebound

- Date: 2026-05-30
- Status: corrected candidate replayed; metadata-only scientist returned `PROCEED`
- Metadata-only HOLD finding:
  - retained `pspec` and `pspec2` headers described quantity, axes, bins, range,
    reduction policy, schema, mode, and units but omitted species cardinality;
  - the README summary of rejected relativistic `psamp` aliases omitted
    `magnetic_moment`, although both writer and inspector rejected it.
- Corrected choices:
  - emit `nspecies` in retained spectrum and joint-spectrum headers;
  - validate optional `nspecies` metadata in the offline readers, preserving
    legacy artifact compatibility while binding the new schema;
  - require `nspecies = 1` in the tracked relativistic runtime qualifier;
  - add spectrum and joint-spectrum species-cardinality mutation controls;
  - align the README alias summary with the executable policy.
- Alternative rejected:
  - relying on out-of-band caller knowledge of `nspecies` would violate the
    metadata-only scientist gate.
- Inflection point:
  - rerun the Phase-7 source and metadata-only rebounds after evidence replay;
    do not open Phase 8 widening on a partially self-describing spectrum
    schema.

## DR-038: Phase-7 Third Rebound For Zero-Field Diagnostic Fallbacks

- Date: 2026-05-30
- Status: corrected candidate replayed; provenance rebound pending
- Fresh source-adversary HOLD findings:
  - retained `pspec` `mu` and `velocity_magnetic_moment_proxy` paths could
    replace non-finite dependencies with a finite zero when `B = 0`;
  - retained `psamp` `mu`, `vpar`, `vperp`,
    `velocity_magnetic_moment_proxy`, and `r_larmor_over_dx_min` paths could
    apply the same zero-field fallback before validating all required
    dependencies;
  - constructor coverage bound representative invalid histogram geometries
    but did not require finite non-overflowing spans for every retained
    histogram family;
  - the archived all-format metadata report named an ad hoc generator label
    rather than a tracked replay tool.
- Corrected choices:
  - validate `pspec` dependencies before applying zero-field fallbacks;
  - return a non-finite sentinel from `psamp` derived-field evaluation whenever
    a required velocity, field, momentum, alpha, or cell-width dependency is
    invalid, so the retained writer fails closed before emitting a row;
  - require finite positive MeshBlock cell widths before evaluating the
    applicability ratio;
  - extend constructor regressions to cover finite overflowing ranges for
    every retained histogram family and bind the dependency-order contract;
  - replace the ad hoc all-format metadata export with a tracked deterministic
    replay tool.
- Alternatives rejected:
  - relying on the final derived value to be non-finite is insufficient because
    a legitimate zero-field fallback can mask an invalid dependency;
  - treating the all-format report as a sealed but manually generated artifact
    would weaken replayability at the scientist-facing gate.
- Inflection point:
  - rerun both CPU builds, the complete Phase-7 evidence replay, legacy
    regressions, and fresh source and metadata-only rebounds before recording
    `RG-008` or opening Phase 8.

## DR-039: Phase-7 Provenance Reseal After Diagnostic Rebound

- Date: 2026-05-30
- Status: corrected candidate resealed; superseded by `DR-040` source rebound
- Fresh metadata-only scientist findings:
  - the tracked all-format exporter reproduces the archived JSON and the
    scientist-facing semantics are adequate;
  - the Phase-7 checksum manifest and source snapshot still described the
    pre-`DR-038` candidate and therefore failed direct verification;
  - the new deterministic exporter must be added to the source snapshot.
- Corrected choices:
  - regenerate the path-sorted source snapshot after the final `DR-038` replay,
    including `scripts/particles/cr_relativistic_all_formats_inspect.py`;
  - regenerate the evidence checksum manifest and archive direct verification;
  - retain the preregistration manifest's pre-replay status string and frozen
    SHA-256 intentionally: rewriting a preregistration artifact after replay
    would erase the evidence that criteria were fixed before the candidate was
    accepted.
- Alternative rejected:
  - changing the frozen criteria status string to describe the replayed state
    would make the accepted replay validate against a post hoc manifest.
- Inflection point:
  - require a narrow metadata-only reseal recheck before recording `RG-008`.

## DR-040: Phase-7 Retained-Writer Adversarial Rebound

- Date: 2026-05-30
- Status: corrected candidate replayed; superseded by `DR-041` source rebound
- Fresh source-adversary HOLD findings:
  - `pmom` used naive magnetic-field squaring and had no relativistic failure
    channel, so finite extreme fields could emit a scientifically wrong
    pitch-angle moment;
  - `ppd` could downcast a finite `Real` position to non-finite `float32`;
  - `psamp` validated derived values while streaming fields, so a late failure
    could leave a malformed partial data row;
  - source-substring checks alone did not bind these writer failure paths.
- Corrected choices:
  - use scaled norms and normalized magnetic-field components when constructing
    retained pitch-angle quantities;
  - validate every relativistic `pmom` contribution before atomics and fail the
    output when any retained moment is non-finite;
  - validate `ppd` float32 representability collectively before opening the
    artifact;
  - evaluate and validate a complete `psamp` row before writing its prefix;
  - add executable adversarial output probes rather than relying only on source
    substring assertions.
- Alternatives rejected:
  - retaining naive dot products divided by extreme magnitudes would preserve
    avoidable intermediate overflow;
  - allowing a partially written diagnostic followed by a fatal exit would
    create an artifact that offline tools cannot interpret safely;
  - applying `ppd` validation after opening the file would leave an ambiguous
    partial artifact.
- Inflection point:
  - rerun both CPU builds, executable adversarial output probes, complete
    Phase-7 replay, legacy regressions, source rebound, metadata-only rebound,
    source snapshot, and evidence reseal before recording `RG-008`.

## DR-041: Phase-7 Retained-Writer Numerical Second Rebound

- Date: 2026-05-30
- Status: corrected candidate replayed; superseded by `DR-042` source rebound
- Fresh source-adversary HOLD findings:
  - retained `df` still computed pitch angle through a raw velocity-field dot
    product divided by `|v||B|`, so finite extreme fields could overflow an
    intermediate denominator and emit a scientifically wrong finite bin;
  - sampled `r_larmor_over_dx_min` squared the momentum magnitude and parallel
    component, then multiplied `alpha_s*dx_min` before dividing, so avoidable
    intermediate overflow could silently turn a representable applicability
    ratio into zero;
  - the first executable writer probe did not bind `df`, `pmom`
    reject-before-artifact behavior, or the representable applicability-ratio
    edge case.
- Corrected choices:
  - route retained `df` through the normalized-component pitch-angle helper
    already shared by the other relativistic outputs;
  - obtain perpendicular momentum from a normalized-field cross product,
    avoiding square-subtract cancellation and overflow;
  - evaluate the three-factor positive applicability ratio with mantissa and
    exponent decomposition, preserving representable subnormal results without
    admitting non-finite final values;
  - extend the tracked writer probe to six cases covering extreme-field
    `pmom`, `pmom` rejection before artifact creation, extreme-field `df`,
    late-field `psamp` row atomicity, subnormal applicability-ratio retention,
    and pre-open `ppd` rejection.
- Alternatives rejected:
  - sequential division alone would still depend on factor ordering and could
    overflow or underflow an intermediate even when the final ratio is
    representable;
  - retaining squared perpendicular-momentum reconstruction would keep both
    avoidable overflow and cancellation near field-aligned states.
- Inflection point:
  - replay the corrected Phase-7 evidence and require a fresh source rebound
    before generating the final source snapshot or recording `RG-008`.

## DR-042: Phase-7 Field-Relative Diagnostic Consolidation

- Date: 2026-05-30
- Status: corrected candidate replayed; superseded by `DR-043` source rebound
- Fresh source-adversary HOLD findings:
  - physical-velocity perpendicular diagnostics still reconstructed
    `v_perp^2` by subtracting nearly equal squares, creating avoidable
    cancellation for almost field-aligned particles;
  - `pspec mu` and `pspec2 mu_wmag` rejected representable diagnostics when an
    unused raw `speed2` intermediate overflowed;
  - zero-field pitch-angle and Larmor applicability outputs were encoded as an
    ordinary zero even though the field-relative diagnostic is undefined.
- Corrected choices:
  - compute physical `v_perp` through the normalized-field cross-product
    helper and reuse it for sampled fields, `pspec`, and `pspec2`;
  - validate only intermediates required by the selected spectrum quantity;
  - fail closed for relativistic field-relative diagnostics when the sampled
    magnetic-field magnitude is zero, including pitch angle, field-relative
    velocity decomposition, physical velocity magnetic-moment proxy, and
    `r_larmor_over_dx_min`;
  - change the Phase-7 prescribed diagnostic fixture from `B = 0` to a nonzero
    field aligned with the initial momentum and electric field, retaining the
    analytical acceleration and deceleration trajectories while making every
    retained field-relative diagnostic scientifically defined;
  - extend the tracked writer adversarial probe to bind cancellation,
    scoped-intermediate validation, and zero-field rejection.
- Alternatives rejected:
  - documenting a zero-field sentinel would still conflate an undefined field
    direction with a resolved zero in histograms and moments;
  - retaining square-subtract reconstruction for physical velocity would leave
    an avoidable cancellation defect even after the momentum-side repair.
- Inflection point:
  - rebuild both CPU binaries, replay all Phase-7 evidence and legacy suites,
    and require another fresh source rebound before generating the final
    snapshot or recording `RG-008`.

## DR-043: Phase-7 Proxy Scaling And Offline Histogram Rebound

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - physical-velocity magnetic-moment proxy evaluation used
    `v_perp*(v_perp/B)`, so a tiny positive `B` could overflow the intermediate
    quotient even when the documented final `v_perp^2/B` value remained
    representable;
  - offline histogram readers accepted impossible negative int32 bins when the
    remaining bins preserved the expected aggregate count.
- Corrected choices:
  - evaluate `v_perp^2/B` through a device-callable exponent-scaled helper
    built from Kokkos `logb` and `exp2`, normalizing the final mantissa before
    applying its exponent;
  - reuse the helper in both device-side `pspec` and host-side `psamp`;
  - reject negative bins in the shared histogram record decoder before
    retained-family reshaping or aggregate validation;
  - add tiny-`B` executable proxy probes for `pspec` and `psamp`, and mutation
    controls for negative bins in `df`, `dxh`, `drh`, `dparh`, `pspec`, and
    `pspec2`.
- Alternatives rejected:
  - quotient-first and square-first evaluation each fail in a different
    representable regime;
  - aggregate-count validation alone cannot distinguish a physical histogram
    from cancelling negative and positive bin corruption.
- Inflection point:
  - rebuild, replay all Phase-7 evidence, and require fresh source and
    metadata-only rebounds before final snapshot generation or `RG-008`.

## DR-044: Phase-7 Offline Schema Identity Rebound

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - the retained `ppd` reader ignored its header particle count and rounded
    fractional species identities;
  - the retained `pmom` reader accepted headerless payloads, rounded
    fractional species and count identities, and silently selected a suffix
    when the latest block contained extra rows;
  - retained histogram readers decoded arbitrary layout declarations without
    checking the family-specific binary layout.
- Corrected choices:
  - bind the `ppd` header particle count to payload cardinality and require
    finite nonnegative integral species identities;
  - require a recognizable latest `pmom` header, exact latest-block row count,
    and finite nonnegative integral species and particle-count identities;
  - reject contradictory histogram layouts for every retained histogram
    family, require layouts for `relativistic_hc` histograms, and retain
    compatibility with legacy histograms that predate layout declarations;
  - extend mutation controls for count mismatches, fractional identities,
    headerless moments, extra moment rows, contradictory family layouts, and
    missing relativistic layout metadata.
- Alternatives rejected:
  - rounding stored identities would convert malformed artifacts into
    plausible scientific records;
  - requiring layouts for every historical histogram would unnecessarily
    break legacy inspection even when the older binary layout is unambiguous.
- Inflection point:
  - replay all Phase-7 evidence and require fresh source and metadata-only
    rebounds before final snapshot generation or `RG-008`.

## DR-045: Phase-7 Declared-Modern Metadata Rebound

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - offline readers still accepted artifacts that declared `relativistic_hc`
    while omitting required schema, units, geometry, or semantic metadata;
  - `pmom` still admitted a legacy-width row even when its metadata declared
    the modern relativistic schema;
  - the shared metadata parser silently discarded malformed standalone tokens;
  - the first layout validator incorrectly flattened the intentionally
    distinct `df` and `dxh` writer declarations to a generic scalar layout.
- Corrected choices:
  - validate required common and family-specific metadata whenever an artifact
    declares either `mode=relativistic_hc` or `schema=akcr_particle_output_v1`;
  - require modern `pmom` rows to carry the complete 26-column retained
    transport-moment payload while preserving legacy 14-column inspection;
  - reject malformed standalone metadata tokens and empty values centrally;
  - bind `df`, `dxh`, scalar histograms, `pspec`, and `pspec2` to their actual
    emitted family layouts, while retaining positive compatibility controls
    for layout-free legacy artifacts;
  - add mutation controls for fractional `pmom` species, malformed metadata,
    underspecified modern artifacts, modern `pmom` legacy-width rows, and
    explicit legacy no-layout acceptance.
- Alternatives rejected:
  - relying only on the runtime qualifier's complete-metadata checks would
    leave the reusable offline reader fail-open for independently supplied
    artifacts;
  - erasing semantic layout suffixes from the writer would discard useful
    retained-format documentation to simplify the inspector.
- Inflection point:
  - replay the full Phase-7 evidence set and require fresh source and
    metadata-only rebounds before final snapshot generation or `RG-008`.

## DR-046: Phase-7 Offline Identity Closure

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - direct modern `pmom` inspection accepted duplicate, reordered, or
    out-of-range species rows even though the higher-level validation wrapper
    rejected them later;
  - direct modern `psamp` inspection accepted negative row identities and a
    stored row rank inconsistent with the shard path;
  - the first `DR-045` provenance refresh was necessarily superseded by these
    source-adversary findings.
- Corrected choices:
  - require modern `pmom` row `n` to describe species `n` inside the shared
    record reader, not only in a higher-level optional validator;
  - require modern `psamp` rank, packet-global ID, tag, and species identities
    to be nonnegative and require row rank to match the shard path;
  - add direct mutation controls for duplicate, reordered, and out-of-range
    modern moment species, every sampled-row negative identity, and sampled
    shard-rank mismatch.
- Alternative rejected:
  - leaving identity checks only in scientist-specific wrappers would make the
    reusable shared reader unsafe for independent offline consumers.
- Inflection point:
  - replay the complete affected evidence set, regenerate provenance from the
    corrected candidate, and require another fresh source rebound before
    recording `RG-008`.

## DR-047: Phase-7 Sample-Shard And Spectrum-Semantic Closure

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - direct modern `psamp` inspection inferred rank zero for unsharded files and
    accepted arbitrary path substrings rather than the writer's canonical
    `psamp/rank_########/` shard structure;
  - direct modern `psamp` inspection retained but did not enforce header rank
    bounds, selector geometry, SplitMix64 selection, selected species, or the
    writer's signed-integer identity range;
  - direct modern spectrum readers required unit fields but did not bind those
    units to the declared relativistic quantity or require the logarithmic
    spectrum floor.
- Corrected choices:
  - require exact canonical modern sample shard paths, a valid header rank
    count, and shard rank strictly below the declared rank count;
  - parse and validate modern selector metadata, replay `sample_selected()` for
    every retained row, and require sampled identities to fit the writer's
    signed-integer contract;
  - bind `pspec` and `pspec2` units to their declared relativistic quantities
    and require a finite positive `log10_floor` for logarithmic spectra;
  - add direct mutations for malformed shard paths, rank bounds, selector
    geometry and predicate drift, integer upper bounds, contradictory spectrum
    units, and missing logarithmic floor metadata.
- Alternatives rejected:
  - preserving permissive sample path inference for modern artifacts would
    admit files the retained writer cannot produce;
  - treating units as decorative metadata would leave scientifically
    contradictory spectra readable as valid records.
- Inflection point:
  - replay the complete affected Phase-7 evidence set and require another
    hostile source rebound before regenerating provenance or recording
    `RG-008`.

## DR-048: Phase-7 Retained-Format Token And Floor-Unit Closure

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source- and metadata-adversary HOLD findings:
  - direct `ppd` and modern `psamp` readers accepted a leading integer prefix
    from duplicated or fractional header-count tokens instead of requiring the
    writer's single complete integer token;
  - direct modern `ppd` inspection accepted integral float32 species values
    outside the signed writer identity domain, and the writer did not reject a
    species identity if its float32 serialization lost exactness;
  - logarithmic relativistic spectra declared the floor value but not the
    floor's units, and the user-facing retained-format documentation did not
    state the exact floored logarithm.
- Corrected choices:
  - require exactly one complete nonnegative integer `particles` token in
    `ppd` headers and one complete positive integer `nranks` token in modern
    `psamp` headers;
  - reject modern `ppd` species identities outside the signed writer domain
    and reject writer-side float32 identity serialization unless the stored
    value exactly round-trips to the integer species identity;
  - emit and bind `log10_floor_units=code_velocity_squared`, document
    `log10(max(kinetic_energy_model, log10_floor))`, and add a positive
    writer-emitted logarithmic-spectrum probe plus missing and contradictory
    floor-unit mutation controls.
- Alternatives rejected:
  - accepting a syntactically valid integer prefix would admit malformed
    artifacts that the writer cannot emit;
  - treating logarithmic-floor units as inferable from the quantity name would
    leave a scientific interpretation contract implicit.
- Inflection point:
  - replay the complete Phase-7 evidence set, require fresh hostile source and
    metadata rebounds, and regenerate provenance only from the accepted
    corrected candidate.

## DR-049: Phase-7 Reserved-Token Downgrade Closure

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh source-adversary HOLD findings:
  - direct `ppd` and modern `psamp` header parsing admitted punctuation-prefixed
    lookalike labels such as `x-particles` and `x-nranks`;
  - malformed reserved `schema` and `mode` values could silently bypass the
    modern retained-format validator and downgrade an artifact to permissive
    legacy inspection.
- Corrected choices:
  - recognize reserved header-count labels only at line start or after
    whitespace;
  - treat the presence of either reserved `schema` or `mode` metadata key as a
    request for strict modern validation, then reject any contradictory value;
  - add direct prefixed-label and malformed-marker downgrade mutations for
    `ppd`, `psamp`, and `pspec`.
- Alternative rejected:
  - retaining exact-value-only modern detection would allow a one-character
    marker corruption to disable the fail-closed contract.
- Inflection point:
  - replay the affected offline suite, regenerate final Phase-7 provenance,
    and require a fresh hostile rebound before recording `RG-008`.

## DR-050: Phase-7 Digest And Boundary-Probe Closure

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh metadata-adversary HOLD findings:
  - the deterministic all-format exporter bound runtime metrics to binary and
    input paths but not their bytes;
  - the writer-emitted logarithmic-spectrum probe checked floor units and
    cardinality without checking the emitted floor value or expected zero-state
    bin;
  - the writer's float32 `ppd` identity branch and parser boundaries needed
    more direct controls and the README needed the retained serialization rule.
- Corrected choices:
  - record and verify SHA-256 digests for the binary, input, and criteria in
    Phase-7 runtime metrics and all-format export provenance;
  - bind the zero-state emitted logarithmic spectrum to
    `log10_floor=1e-30`, its declared units, and the exact occupied bin;
  - factor the production float32 species predicate into a `constexpr` helper,
    compile-bind exact and lossy boundaries with `static_assert`, execute a
    positive large exact float32 reader boundary, and add missing, zero,
    nonfinite, and out-of-range parser controls;
  - document that relativistic `ppd` creation fails before artifact creation
    when float32 serialization loses finite positions or exact signed species
    identity.
- Alternative rejected:
  - adding a hidden runtime-only species override to manufacture an impossible
    high species ID would widen production input semantics solely for a test;
    compile-time assertions bind the exact production predicate without such a
    seam.
- Inflection point:
  - replay the full Phase-7 gate, reseal source and evidence manifests, and
    require final bounded hostile source and provenance rebounds before
    recording `RG-008`.

## DR-051: Phase-7 Writer-Evidence Enumeration

- Date: 2026-05-30
- Status: corrected candidate replay pending
- Fresh metadata-adversary HOLD finding:
  - the archived writer-adversarial evidence summarized `19/19` without naming
    the probes, so a later auditor could not verify writer-specific closure
    from the evidence bundle alone.
- Corrected choice:
  - enumerate every writer-adversarial PASS in the archived log, record the
    exact zero-state logarithmic floor, units, and occupied bin, identify
    pre-artifact rejection probes explicitly, and archive the compile-bound
    float32 species predicate boundaries as a distinct twentieth check.
- Alternative rejected:
  - treating a single aggregate PASS count as sufficient would make the sealed
    evidence depend on rereading the test implementation.
- Inflection point:
  - rerun the enumerated writer probe, regenerate provenance, and require the
    final bounded metadata rebound before recording `RG-008`.

## RG-008: Scientists Can Interpret The Retained Outputs

- Date: 2026-05-30
- Status: `PROCEED`
- Scope accepted:
  - relativistic retained writers and direct offline readers now bind explicit
    schemas, units, layouts, counts, identities, logarithmic-floor semantics,
    and digest-bound replay provenance;
  - legacy readers retain their documented compatibility behavior without
    silently accepting malformed modern declarations.
- Evidence:
  - frozen `43`-criterion Phase-7 runtime qualification passes;
  - `88` direct parser, metadata, and provenance mutations reject;
  - enumerated writer-adversarial evidence passes `20/20`;
  - parser contract passes `172` cases;
  - legacy serial and MPI ladders pass `20 + 12 + 4 + 5` cases;
  - style passes `2` cases and `git diff --check` is empty;
  - deterministic all-format export binds binary, input, and criteria
    SHA-256 digests.
- Independent review disposition:
  - source-adversary `HOLD` findings from `DR-045` through `DR-050` were
    corrected and the final bounded source rebound returned `PROCEED`;
  - metadata-adversary evidence-enumeration `HOLD` is corrected by `DR-051`
    and requires one final metadata-only reseal verification before Phase-8
    source edits.
- Residual boundary:
  - MPI, block migration, SMR, AMR, and changed-rank restart policy remain
    closed until the separately preregistered Phase-8 gate passes.

## DR-052: Phase-8 Narrow Migration Qualification Opening

- Date: 2026-05-30
- Status: qualification replay pending
- Selected choices:
  - preserve the serial uniform-level fence for
    `relativistic_field_source=mhd_ideal`;
  - remove that parent fence only from the prescribed-test path so the frozen
    Phase-8 oracle can exercise generic particle migration, MPI exchange, SMR,
    adaptive AMR remap, periodic wrapping, and momentum evolution;
  - preserve the typed-v2 changed-rank restart rejection policy and report its
    exact diagnostic before evaluating the derived topology hash.
- Evidence before edit:
  - frozen `57`-criterion Phase-8 registration validates with SHA-256
    `c11ee798522629e4cc4aa14af3c45bf3df730d84a1eb4f114609651c1d75c645`;
  - `14` Phase-8 oracle mutation controls reject;
  - pre-fence smoke confirms one-block and multiblock serial parity, owner
    changes, periodic wraps, expected MPI/SMR/AMR fences, and the preempted
    changed-rank diagnostic.
- Alternatives rejected:
  - opening solver-coupled `mhd_ideal` migration in the same edit would combine
    generic transport qualification with live-grid coupling scope;
  - adding special migration code before the existing generic loops fail would
    widen the patch without evidence.
- Required review surfaces:
  - inspect `src/bvals/bvals_part.cpp`, `src/mesh/mesh_refinement.cpp`, and
    `Particles::RemapAfterAMR()` even if the generic record loops already cover
    the appended relativistic fields.
- Inflection point:
  - stop if the frozen runtime oracle exposes state drift, stale appended
    fields, missing off-rank witnesses, or changed-rank diagnostic ambiguity.

## DR-053: Phase-8 Static-SMR Oracle Population Correction

- Date: 2026-05-30
- Status: corrected qualifier replay pending
- First post-fence qualifier HOLD:
  - the static-SMR fixture requested `8 / (15 * 512)` particles per cell;
  - serial initialization truncated once and created `8` particles, while the
    four MPI packs truncated independently across `3 + 4 + 4 + 4` leaf blocks
    and created only `7`;
  - the qualifier therefore rejected the MPI static-SMR case before comparing
    transport state.
- Corrected choice:
  - keep the preregistered base, adaptive, continuation, and empty-rank ladders
    at `8` particles;
  - give the static-SMR and matching fine-uniform reference leg a separate
    deterministic `15`-particle identity set using one particle per static-SMR
    leaf block, which is exactly additive across serial and four MPI packs.
- Alternative rejected:
  - weakening identity checks or accepting a seven-particle MPI artifact would
    hide fixture decomposition drift instead of testing migration.
- Frozen-criteria disposition:
  - retain the registered `57` acceptance criteria unchanged; the correction
    repairs the oracle population used to evaluate the existing static-SMR
    parity criteria.
- Inflection point:
  - rerun the complete Phase-8 qualifier and stop again if any physical-state
    comparison, exchange witness, adaptive remap witness, or restart policy
    assertion fails.

## DR-054: Phase-8 Adaptive Oracle Timestep Correction

- Date: 2026-05-30
- Status: corrected qualifier replay pending
- Second post-fence qualifier HOLD:
  - the adaptive fixture capped `tlim` below the default relativistic outer
    timestep bound, so it completed one outer step;
  - the moving-box oracle refined after that step but could not execute a
    second refinement pass to witness derefinement.
- Corrected choice:
  - retain the registered three-cycle adaptive oracle;
  - give the adaptive legs and their matching fine-uniform reference the same
    fixture-local `subcycle_max_steps=4` and
    `subcycle_gyro_fraction=0.0005`;
  - this makes the gyro-limited outer timestep smaller than the finest-cell
    crossing bound while exercising real relativistic subcycling.
- Alternatives rejected:
  - a production-only timestep override would widen runtime semantics solely
    for the qualification fixture;
  - staging an AMR restart continuation would add checkpoint orchestration
    when the existing outer-bound contract can express the required schedule;
  - allowing adaptive and fine-uniform legs to take different outer steps
    would weaken their physical-state comparison.
- Frozen-criteria disposition:
  - retain the registered `57` acceptance criteria unchanged; this repairs the
    execution schedule used to evaluate the existing adaptive AMR criteria.
- Inflection point:
  - rerun the complete Phase-8 qualifier and stop if AMR churn, off-rank remap,
    physical-state parity, or restart policy evidence fails.

## DR-055: Phase-8 Continuation Reference Schedule Correction

- Date: 2026-05-30
- Status: corrected qualifier replay pending
- Third post-fence qualifier HOLD:
  - the first continuation oracle compared an uninterrupted `tlim=0.8` run
    taking one outer step with a restarted path taking `0.6 + 0.2`;
  - the resulting small state difference was expected integrator partition
    sensitivity, not evidence that typed-v2 restart restoration lost state.
- Corrected choice:
  - add a dedicated MPI-4 continuation source instead of overloading the rank
    ladder;
  - apply fixture-local `subcycle_max_steps=4` and
    `subcycle_cell_fraction=0.4` to source, uninterrupted reference, and both
    restart continuations so they all advance through the same `0.1` outer
    steps;
  - retain the post-migration restart source, exact typed-v2 inspection, owner
    migration assertion, same-rank policy, and topology-hash equality checks.
- Alternatives rejected:
  - loosening the state tolerance would conceal an invalid reference schedule;
  - requiring the relativistic integrator to be invariant under arbitrary
    outer-step repartitioning would test a stronger and unrelated property;
  - changing production restart state solely for this mismatch would patch
    code before demonstrating a code defect.
- Frozen-criteria disposition:
  - retain the registered `57` acceptance criteria unchanged; the dedicated
    split reference repairs the schedule used to evaluate restart parity.
- Inflection point:
  - rerun the complete Phase-8 qualifier and stop if the schedule-matched
    post-migration restart path differs from its uninterrupted reference.

## DR-056: Phase-8 Narrow MHD-Ideal Fence Diagnostic Update

- Date: 2026-05-30
- Status: corrected qualifier replay pending
- Fourth post-fence qualifier HOLD:
  - the retained `mhd_ideal` MPI negative control rejected correctly, but its
    expected-diagnostic list still named only the removed parent fence and the
    older one-MeshBlock restriction;
  - the narrowed production fence now reports that `mhd_ideal` requires a
    serial uniform-level mesh.
- Corrected choice:
  - add the exact narrowed `mhd_ideal` runtime-fence diagnostic to the retained
    negative-control allow-list;
  - preserve the historical parent-fence and one-MeshBlock diagnostics so the
    pre-fence and constructor-level rejection paths remain inspectable.
- Alternative rejected:
  - matching a broad substring such as `currently requires` would stop binding
    the intentional scope boundary expressed by the diagnostic.
- Frozen-criteria disposition:
  - retain the registered `57` acceptance criteria unchanged; the production
    rejection was correct and only its expected text contract lagged.
- Inflection point:
  - rerun the complete Phase-8 qualifier and stop if any retained negative
    control unexpectedly executes or reports an unrecognized scope boundary.

## DR-057: Phase-8 Parser Contract Multilevel Opening

- Date: 2026-05-30
- Status: corrected regression replay pending
- Post-qualification regression HOLD:
  - the parser suite still required every prescribed-test relativistic
    multilevel constructor to reject with the removed parent-fence diagnostic;
  - after the intentional Phase-8 opening, its incomplete adaptive fixture
    proceeded far enough to fail for a missing AMR criterion instead.
- Corrected choice:
  - replace that obsolete rejection with a valid constructor-level static-SMR
    acceptance case for `relativistic_field_source=prescribed_test`;
  - retain the separate solver-coupled `mhd_ideal` multilevel rejection test
    unchanged.
- Alternative rejected:
  - preserving the old rejection would contradict the qualified MPI, SMR, and
    AMR migration scope;
  - accepting an unrelated missing-criterion failure would make the parser
    suite pass for the wrong reason.
- Inflection point:
  - replay the parser contract and broader legacy regressions before recording
    `RG-009`.

## DR-058: Phase-8 Evidence-Adversary Strengthening

- Date: 2026-05-30
- Status: strengthened qualifier replay pending
- Independent evidence-adversary `HOLD` findings:
  - the first passing qualifier restarted after same-level MPI migration but
    did not restart after adaptive refinement and derefinement;
  - continuation schedules were not parsed and allowed floating-point
    cleanup cycles near zero timestep;
  - AMR churn parsing scanned setup output as well as runtime output;
  - static-SMR multilevel lookup was inferred from hierarchy plus traffic
    instead of an explicit source-to-destination refinement-level transition;
  - the passing package was not yet sealed with durable replay artifacts and
    SHA-256 provenance.
- Corrected choices:
  - add schedule-matched adaptive restart continuations after AMR churn for
    both `device_table` and `host_tree`;
  - stop the same-level continuation fixtures by exact cycle count, parse the
    intended positive-timestep schedule, and reject near-zero cleanup cycles;
  - accept create/delete churn only after the setup-complete marker;
  - add a diagnostic-only `particle_multilevel_lookup` performance counter
    that counts actual old-level versus destination-level transitions in both
    multilevel lookup implementations and bind `P8-20` to it;
  - add binary, MPI binary, input, criteria, script, and source-HEAD digests to
    the runtime report, then create a durable Phase-8 seal after replay.
- Alternatives rejected:
  - inferring cross-level behavior from MPI sends could pass for same-level
    movement inside an SMR hierarchy;
  - weakening restart tolerances would conceal reference-schedule defects;
  - preserving temporary paths as the only replay evidence would make the
    qualification package non-durable.
- Frozen-criteria disposition:
  - retain the registered `57` acceptance criteria unchanged; the new checks
    strengthen witnesses and add defense-in-depth metrics without weakening a
    threshold.
- Inflection point:
  - rebuild, replay the strengthened qualifier and mutation controls, archive
    durable evidence, and require a fresh evidence-adversary recheck before
    recording `CP-6` and `RG-009`.

## DR-059: Phase-8 Cycle-Parser Token Grammar Correction

- Date: 2026-05-30
- Status: corrected qualifier replay pending
- Strengthened-qualifier HOLD:
  - the new cycle-schedule regular expression used a literal-space exclusion,
    so its timestep capture consumed a newline and the next `Particle` token;
  - static-SMR cross-level evidence and adaptive restart-after-AMR legs had
    already executed successfully before this parser stop.
- Corrected choice:
  - parse elapsed, time, and timestep values as non-whitespace tokens.
- Alternative rejected:
  - trimming the captured value after a permissive match would keep accepting
    malformed schedule lines.
- Inflection point:
  - rerun script syntax, mutation controls, and the complete strengthened
    Phase-8 qualifier.

## DR-060: Phase-8 Diagnostic Isolation and Host-Tree Witness Amendment

- Date: 2026-05-30
- Status: corrected qualifier replay and fresh source rebound passed
- Fresh source-rebound `HOLD` findings:
  - the first cross-level counter allocated and copied the diagnostic level
    table and dereferenced source levels even when performance logging was
    disabled;
  - the device-table static-SMR fixture asserted a cross-level transition, but
    the host-tree implementation was only exercised indirectly by adaptive
    continuations.
- Corrected choices:
  - allocate the device level table and read source levels only when
    `particles/log_performance=true`;
  - validate source GID ranges before diagnostic source-level lookups in both
    the device-table and host-tree implementations;
  - add a dedicated MPI-4 static-SMR `host_tree` fixture with exact identity,
    field-parity, off-rank-send, and cross-level-transition assertions;
  - amend the registered criteria from `57` to `62` checks and add a mutation
    control that must reject loss of the host-tree transition count.
- Alternatives rejected:
  - leaving the counter active outside diagnostics would impose qualification
    overhead and new reads on ordinary multilevel migration;
  - treating host-tree observation in an archived log as sufficient would
    allow the executable oracle to pass after a future host-tree regression;
  - inferring host-tree behavior from adaptive AMR traffic alone would not
    prove an actual source-to-destination refinement-level transition.
- Registration disposition:
  - the criteria amendment strengthens the oracle after independent rebound
    review; it does not loosen an existing threshold or widen production
    scope;
  - bind the amended manifest with SHA-256
    `add2ab12fbaf7f669e8e4eecc9e8f9cb0d581f06947c7e39876f934239708b74`.
- Inflection point:
  - rerun registration, mutation controls, both builds, and the complete
    Phase-8 qualifier; stop again if either lookup implementation lacks an
    explicit cross-level witness or if the source rebound remains on hold.

## Phase-8 DR-060 Source Rebound Recheck

- Date: 2026-05-30
- Reviewer verdict: `PROCEED`
- Reviewer scope:
  - diagnostic-only level-tracking isolation;
  - source-index bounds checks;
  - multilevel empty-rank MPI collective cadence;
  - explicit static-SMR device-table and host-tree cross-level assertions;
  - bounded `prescribed_test` opening and retained `mhd_ideal` closure;
  - typed-v2 rank-count rejection ordering;
  - amended criteria/schema coherence.
- Recheck evidence:
  - frozen registration passed with `62` ordered criteria;
  - mutation controls rejected `17` deliberate weakenings;
  - serial and MPI builds passed;
  - complete corrected runtime replay passed `62/62` criteria across `21`
    isolated runtime cases;
  - device-table and host-tree static-SMR fixtures each recorded one explicit
    cross-level transition and exact physical-state parity.
- Disposition:
  - accept the corrected source rebound for a source-bound milestone commit;
  - rebuild and rerun the Phase-8 evidence package on that accepted source SHA
    before recording `CP-6` and `RG-009`.

## CP-6 AMR And MPI Seal

- Date: 2026-05-30
- Accepted source SHA:
  `addd12d4e26f7d8b275165b6be7b364d39f22a43`
- Verdict: `PROCEED`.
- Bounded production opening:
  - `relativistic_hc` plus `relativistic_field_source=prescribed_test` is open
    for periodic MPI, same-level multiblock, static SMR, and adaptive AMR
    qualification;
  - solver-coupled `relativistic_field_source=mhd_ideal` remains fail-closed for
    MPI, SMR, and AMR;
  - nonperiodic execution, changed-rank redistribution, GPU execution, and
    scientific production-readiness claims remain deferred.
- Accepted-SHA evidence:
  - serial and MPI builds passed;
  - frozen amended registration validated `62` ordered criteria with SHA-256
    `add2ab12fbaf7f669e8e4eecc9e8f9cb0d581f06947c7e39876f934239708b74`;
  - mutation controls rejected `17` deliberate weakenings;
  - complete replay passed `62/62` criteria across `21` isolated runtime cases;
  - MPI rank ladder `1`, `2`, `4`, and `8` passed with particle-empty-rank
    evidence;
  - periodic wrap, nonzero prescribed field, momentum-change, and work
    witnesses passed;
  - static-SMR device-table and host-tree fixtures each recorded exact physical
    parity, off-rank movement, and one explicit cross-level transition;
  - adaptive moving-box refinement and derefinement passed for both remap
    implementations with off-rank sends and fine-uniform reference agreement;
  - restart-after-AMR continuations passed for both remap implementations with
    exact physical-state agreement;
  - schedule-matched same-topology continuations before and after migration
    passed exactly;
  - changed-rank typed-v2 continuation rejected with the exact rank-count
    diagnostic before the derivative topology diagnostic;
  - parser contract passed `172/172`;
  - retained legacy particle suites passed `20/20`, `12/12`, `4/4`, and `5/5`.
- Durable replay package:
  - `phase8_replay_bundle.tar.gz` retains selected static-SMR, adaptive-AMR,
    restart-after-AMR, same-topology continuation, and changed-rank rejection
    roots;
  - its sorted archive inventory is retained;
  - fresh extraction and offline typed-v2 decode passed for `29` retained
    checkpoints.

## DR-061: Phase-8 Durable Replay Decoder Completeness

- Date: 2026-05-30
- Status: corrected verifier replay and evidence rebound passed
- Fresh evidence-adversary `HOLD` finding:
  - the durable replay tarball correctly retained the changed-rank rejection
    root, but the retained offline decoder enumerated only the positive static,
    adaptive, and same-topology continuation roots;
  - the archive therefore contained `29` valid typed-v2 checkpoints while the
    retained decoder proved only `28`.
- Corrected choice:
  - add `mpi_rank_count_change_rejection` to the retained offline decode list;
  - update the CP-6 candidate count to `29`;
  - regenerate the direct path-sorted SHA-256 manifest and verification log;
  - require a fresh evidence rebound after the corrected verifier replay.
- Alternative rejected:
  - accepting the adversary's independent decode without correcting the
    retained executable evidence would leave the seal weaker than the archive
    it describes.
- Inflection point:
  - stop again if the changed-rank root fails typed-v2 decoding or if the fresh
    evidence rebound identifies another durability or provenance gap.

## Phase-8 Evidence-Adversary Rebound Sequence

- Date: 2026-05-30
- First evidence-adversary verdict: `HOLD`
  - the durable tarball retained the changed-rank rejection root, but the
    retained offline decoder omitted that root and proved only `28` of its `29`
    typed-v2 checkpoints.
- Correction:
  - add the omitted changed-rank rejection root to the retained decoder;
  - regenerate the replay verification log and evidence manifest.
- Second evidence-adversary verdict: `HOLD`
  - corrected replay durability, provenance, executable assertions, and
    retained logs passed;
  - the living-ledger header still described a stale Phase-5 state, and
    `DR-061` still said pending after its verifier replay passed.
- Correction:
  - refresh the living-ledger control header for the accepted Phase-8 scope;
  - mark `DR-061` complete;
  - require one final ledger-focused read-only rebound before committing the
    Phase-8 seal.

## RG-009: Is Infrastructure Noise Hiding Physics Error?

- Date: 2026-05-30
- Verdict: `PROCEED`.
- Reassessment:
  - kernel-only, gather, coupled, subcycle, restart, diagnostic, parser, and
    legacy milestones remained sealed before the bounded migration opening;
  - new migration-only failures were isolated in fixture scheduling,
    diagnostic expectations, parser grammar, and witness strength rather than
    hidden by relaxed physical tolerances;
  - every discovered fixture defect stopped the replay and produced a decision
    record before qualification resumed;
  - accepted adaptive and restart continuations now agree exactly, so the
    migration infrastructure does not introduce a detectable physics drift in
    the qualified prescribed-test scope.
- Deferred risks:
  - solver-coupled `mhd_ideal` MPI, SMR, and AMR widening remains closed;
  - GPU qualification remains unavailable locally and unclaimed;
  - changed-rank redistribution remains rejected pending separate design;
  - public documentation overlay and final cold review remain Phase-9 work.
- Inflection point:
  - require a fresh evidence-adversary `PROCEED` before sealing CP-6;
  - if the evidence review finds a source-bound provenance, replay durability,
    or executable-oracle gap, reopen Phase 8 rather than advancing to Phase 9.

## DR-062: Phase-9 Workstation Backend Disposition

- Date: 2026-05-30
- Status: selected and archived
- Choice:
  - keep `GPU QUALIFIED` explicitly unclaimed on this workstation;
  - allow workstation `MERGE READY` assessment to proceed through the optional
    backend policy selected in `DR-000`;
  - archive a separate GPU testing handoff that must run on an accepted SHA
    before any accelerator-readiness claim.
- Evidence:
  - the Apple M4 Max workstation exposes Metal graphics but no configured CUDA
    or HIP Kokkos backend;
  - `nvcc`, `hipcc`, `nvidia-smi`, and `rocminfo` are unavailable;
  - accepted release, MPI release, debug-bounds, and MPI debug-bounds caches
    enable Kokkos Serial and disable CUDA, HIP, and SYCL.
- Alternative rejected:
  - treating Metal presence as GPU qualification would claim an unconfigured
    backend that AthenaK did not build or execute;
  - blocking the workstation handoff on unavailable optional hardware would
    contradict the Phase-0 backend policy without adding executable evidence.
- Residual risk:
  - accelerator compilation, execution, CPU/GPU differential replay,
    device-aware MPI transport, and profiling remain unqualified follow-up
    work.

## DR-063: Phase-9 Public Documentation Overlay Scope

- Date: 2026-05-30
- Status: selected and validated
- Choice:
  - validate the public page update in a temporary detached
    `origin/gh-pages` worktree;
  - overlay `particles.md`, the complete stacked CR tracer markdown set, and
    the inherited CR tracer accuracy figures required by linked pages;
  - run `make clean html SPHINXOPTS="-W --keep-going"` under `docs/`.
- Correction during validation:
  - the first invocation ran from the repository root, but the public branch
    keeps its Sphinx `Makefile` under `docs/`;
  - the second invocation exposed missing inherited CR accuracy figures;
  - the corrected stacked overlay included those prerequisite images and
    passed with warnings treated as errors.
- Alternatives rejected:
  - validating only the new relativistic markdown files would not represent
    the linked public integration surface;
  - accepting a warning-tolerant build would hide broken public assets.

## DR-064: Phase-9 Unsupported Single-Precision Build Disposition

- Date: 2026-05-30
- Status: archived inherited blocker
- Choice:
  - rerun and retain a final `Athena_SINGLE_PRECISION=ON` build attempt;
  - classify single precision as unsupported for this handoff because the
    build fails before CR source in inherited coordinate headers;
  - keep single-precision repair outside this bounded feature branch.
- Evidence:
  - configuration succeeded;
  - compilation failed with pre-existing C++ narrowing errors in
    `src/coordinates/cartesian_ks.hpp`, including initializer-list conversions
    from `double` to `Real` when `Real` is `float`.
- Alternative rejected:
  - widening Phase 9 into repository-wide coordinate-header cleanup would
    obscure the CR feature boundary and invalidate the final cold-review
    surface.

## Phase-9 Candidate Qualification Before Final Cold Review

- Date: 2026-05-30
- Candidate parent seal:
  `aa8663e8a5d49e26c206363d028b52d0e350a91f`
- Runtime and build results:
  - serial release, MPI release, debug-bounds, and MPI debug-bounds general
    builds passed;
  - dedicated release pusher, release coupled, and debug-bounds coupled builds
    passed;
  - prescribed pusher, solver-coupled, subcycle, restart, and diagnostics
    analyzers passed, including debug-bounds replays for coupled, subcycle,
    restart, and diagnostics paths;
  - release and debug-bounds prescribed-test migration qualification each
    passed the frozen `62/62` contract across `21` isolated runtime cases;
  - migration mutation controls rejected `17/17`, diagnostic mutation controls
    rejected `88/88`, subcycle controls rejected `8/8`, coupled-oracle
    mutations rejected as expected, and retained writer adversarial probes
    passed `20/20`;
  - release and debug-bounds all-format inspection passed;
  - deterministic release all-format analysis rerun from the retained tarball
    matched the original report byte-for-byte;
  - parser contract passed `172/172`;
  - legacy CPU suites passed `20/20` and `12/12`;
  - legacy MPI suites passed `4/4` and `5/5`;
  - style passed `2/2`;
  - inherited deep AMR `divB` regression passed `1/1`;
  - strict public `origin/gh-pages` overlay build passed with warnings as
    errors.
- Integration result:
  - refreshed merge-tree audit against
    `origin/development@c6a73b08e60807f8b925164c5e7edd5cb820c8ae`
    remains intentionally unmerged and reports seven manual content-conflict
    paths;
  - preserve the stacked prerequisite and use a dedicated reviewed
    integration step before merging to `development`.
- Inflection point:
  - commit and push this Phase-9 candidate documentation and evidence;
  - ask fresh cold reviewers with no implementation ownership to try to
    falsify physics, migration, evidence durability, public documentation, and
    handoff claims;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` until blocking reviewer
    findings are closed.

## DR-065: Phase-9 First Cold-Review HOLD And Replacement Seal

- Date: 2026-05-30
- Status: corrected replacement packet prepared for rebound
- Candidate reviewed:
  `2a18501df4ab64a60283f18a63e6b090c8516b1a`
- Independent-review result:
  - physics and migration source audits found no implementation blocker within
    the bounded passive scope;
  - physics, migration, evidence-adversary, and handoff reviewers each returned
    `HOLD` on the release packet;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` on `2a18501d`.
- Public-documentation corrections:
  - add the required `relativistic_initial_state = velocity` selector to the
    solver-coupled example;
  - document `subcycle_electric_kick_max`, `relativistic_initial_state`,
    direct `w0x`, `w0y`, `w0z` initialization, prescribed `cE0x`, `cE0y`,
    `cE0z`, and the relativistic `B_profile = uniform` restriction;
  - state that non-strict cap clipping is legacy-only and that
    `relativistic_hc` requires `subcycle_strict = true`;
  - separate legacy independent particle restart from typed-v2 paired mesh and
    particle restart, including changed-rank rejection;
  - gate public stable-site publication on the later reviewed integration
    branch;
  - relabel the implementation-guide control header as the historical Phase-0
    creation snapshot and point readers to this living ledger;
  - add concrete bounded follow-up branch names and entry gates;
  - add the missing GPU all-format replay command and require debug-bounds
    accelerator replays.
- Evidence-packet corrections:
  - archive exact analytical replay commands instead of ellipses;
  - add one source-bound Phase-9 provenance envelope covering source SHA,
    source tree, binaries, analyzers, inputs, criteria, metrics, frozen
    `development`, and frozen `gh-pages`;
  - replace the nonportable all-format tarball with a deterministic portable
    package containing a prefix-mapped qualified release binary, runtime
    artifacts, input, criteria, required inspectors, manifest, and replay
    wrapper;
  - prove the portable replay from a fresh extraction path;
  - redact local checkout and workstation identity from retained text
    artifacts and record a package privacy scan;
  - reduce the GPU disposition artifact to backend facts only;
  - rerun the Pages overlay from immutable
    `origin/gh-pages@4833aa9341e19861297e330ff02aabfd8001935c` and archive an
    overlay inventory;
  - regenerate the merge-tree artifact in the intended
    `development <- candidate` direction;
  - normalize retained text-log trailing whitespace and archive a real
    commit-range `git diff --check`.
- Additional executable check:
  - run the public solver-coupled documentation example with the general
    release binary and retain its successful runtime log.
- Alternatives rejected:
  - treating rendering success as proof of a runnable public example;
  - accepting opaque one-line replay logs without a cross-artifact provenance
    envelope;
  - retaining a package that requires an ephemeral external binary path;
  - publishing host-specific identity or workstation inventory in the handoff;
  - force-recording `PROCEED` before a fresh rebound reviews the corrected
    immutable seal.
- Inflection point:
  - seal and push a corrected immutable packet;
  - rerun a fresh physics, migration, evidence-adversary, and handoff cold
    review against that exact SHA;
  - stop again on any blocking finding.

## DR-066: Phase-9 First Rebound Seal-Metadata Correction

- Date: 2026-05-30
- Status: corrected packet prepared for second rebound
- First sanitized evidence seal reviewed:
  `89084c66b149e9b310a06e13eb4465c58240a4aa`
- Fresh handoff-rebound verdict: `HOLD`.
- Findings:
  - the handoff named the corrected documentation candidate
    `720fb8fc193a8d598dccbf780cee3777ebba8bc9` but did not distinguish the
    immutable evidence-seal tip `89084c66b149e9b310a06e13eb4465c58240a4aa`;
  - the merge-tree artifact remained bound to the documentation candidate
    rather than the first sanitized evidence-seal tip;
  - the Pages overlay inventory preceded final handoff and ledger edits and
    therefore did not match the exact sealed overlay set;
  - the GPU follow-up handoff required debug-bounds coupled replay without
    giving a dedicated coupled-debug build;
  - the empty whitespace-audit artifact did not state its intended Phase-9
    packet range or disclose inherited archive noise in the older stacked
    range.
- Corrected choices:
  - preserve separate labels for the documentation candidate and each
    immutable evidence-seal tip;
  - rerun the merge-tree artifact in the intended `development <- seal-tip`
    direction;
  - rerun the frozen-`gh-pages` preview overlay after the complete seal-doc
    update and archive the matching inventory;
  - add an explicit dedicated coupled GPU debug-bounds build and exact binary
    selection for required debug replays;
  - record the scoped Phase-9 `git diff --check` command and disclose that
    older stacked retained evidence still contains historical whitespace.
- Alternative rejected:
  - treating a two-stage evidence seal as self-explanatory would force future
    reviewers to infer candidate roles from Git history and leave the handoff
    ambiguous.
- Inflection point:
  - reseal and push the corrected packet;
  - require fresh read-only rebound review against the new immutable tip;
  - stop again on any blocking finding.

## DR-067: Phase-9 Second Rebound Portability And Public-Replay Correction

- Date: 2026-05-30
- Status: corrected packet prepared for sanitized-history replacement
- Second rebound evidence seal reviewed:
  `5e031387e66224b0e9dc4462fbf4d9a7ee01c9df`
- Independent-review result:
  - one integration reviewer returned `PROCEED`;
  - physics, evidence-adversary, and documentation reviewers returned `HOLD`;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` on `5e031387`.
- Findings:
  - the outer Phase-9 manifest and its verification log used absolute local
    checkout paths, so they were not offline-portable and leaked workstation
    identity despite the retained privacy claim;
  - the provenance envelope retained one stale Pages-build digest;
  - the public documentation described the solver-coupled selectors but did
    not ship or link the complete runnable one-cycle deck;
  - the typed-v2 restart example exposed the internally injected paired-mesh
    key without showing the required full-mesh `-r` invocation;
  - the retained command inventory omitted package construction, manifest
    generation, envelope verification, privacy scanning, and final outer
    verification commands;
  - the archived whitespace check stopped at the earlier candidate;
  - portable replay reports remained semantically equivalent but encoded
    extraction-local paths.
- Corrected choices:
  - publish a complete one-cycle solver-coupled input under
    `inputs/particles/` and execute that public path in the retained smoke;
  - document typed-v2 resume as paired `-r` plus particle-shard input and mark
    `relativistic_paired_mesh_restart_file` as internal plumbing;
  - normalize portable report prefixes to `$PACKAGE_ROOT`;
  - regenerate the outer manifest with repository-relative paths, verify it
    from a relocated checkout root, and run the privacy scan only after all
    retained text artifacts exist;
  - archive package-construction, deterministic metadata normalization,
    internal-manifest, outer-manifest, envelope-verification, privacy-scan,
    and final verification commands;
  - bind the final whitespace audit to the latest pre-seal candidate and ask
    the rebound reviewers to rerun it through the later seal tip;
  - rewrite the Phase-9 documentation history so the leaked local identity is
    removed from the branch history rather than merely corrected at the tip.
- Alternatives rejected:
  - retaining an absolute-path manifest because it verifies in the original
    checkout would not prove offline replay;
  - documenting an internal restart key as a user parameter would preserve an
    operationally incomplete resume procedure;
  - accepting path-varying portable reports would weaken deterministic replay
    evidence without providing useful information.
- Inflection point:
  - construct a sanitized replacement candidate from the accepted Phase-8
    evidence seal;
  - regenerate and verify the complete Phase-9 envelope and portable packet;
  - force-update the remote branch with an exact lease against the leaked tip;
  - require a fresh exact-SHA cold review before recording closure.

## DR-068: Phase-9 Portability-Rebound Documentation Durability Correction

- Date: 2026-05-30
- Status: documentation-only correction prepared for final rebound
- Portability-rebound evidence seal reviewed:
  `aa66f6c27531116e12554631281c8f2ed07d93c6`
- Independent-review result:
  - evidence durability reviewer returned `PROCEED`;
  - physics, integration, and handoff reviewers returned `HOLD`;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` on `aa66f6c2`.
- Findings:
  - the handoff Branch Facts table still presented removed pre-rewrite Phase-9
    commits as the live lineage and left the top-level verdict stale;
  - the typed-v2 resume command named an untracked input deck;
  - the replay inventory initialized `$REPO_ROOT` after first use, quoted the
    prefix-map flags so shell expansion would not occur, and omitted producer
    commands or retained-output capture for several generated artifacts;
  - the retained recipe ran its privacy scan before creating all final text
    verification artifacts;
  - the envelope listed ephemeral workstation binary hashes that its retained
    verifier could not reconstruct independently.
- Corrected choices:
  - distinguish historical removed Phase-9 tips from the live sanitized
    candidate and portability-rebound seal;
  - ship a tracked serial prescribed-test typed-v2 resume template and
    document paired native mesh `-r` plus particle-shard cycle selection;
  - define `$REPO_ROOT` and `$EVID` once before use, expand prefix-map flags,
    and retain producer commands for snapshots, overlay artifacts, public
    smoke, merge-tree output, and whitespace output;
  - run the privacy scan after the outer manifest and both verification logs
    exist, then rerun the two manifest verification commands to prove that the
    scan output remained stable;
  - narrow the reconstructable envelope to the archived portable binary,
    retained analyzers, criteria, inputs, result artifacts, and bundle.
- Alternative rejected:
  - preserving unverifiable ephemeral binary hashes would overstate the
    durable packet;
  - documenting an untracked illustrative resume deck would keep the public
    procedure non-runnable.
- Inflection point:
  - reseal the documentation-only correction;
  - rebuild the strict frozen Pages overlay including both public inputs;
  - require another exact-SHA cold review before recording closure.

## DR-069: Phase-9 Documentation-Durability Rebound Capture Correction

- Date: 2026-05-30
- Status: documentation-only capture-completeness correction prepared for
  final cold review
- Documentation-durability correction candidate:
  `9a270755f2739398d61023fa7a950add2dd550e0`
- Documentation-durability seal reviewed:
  `88d631e4943648fe83f0624cb30291fa52ab4296`
- Independent-review result:
  - physics and evidence-durability reviewers returned `PROCEED`;
  - integration and documentation reviewers returned `HOLD`;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` on `88d631e4`.
- Findings:
  - the handoff still named `aa66f6c2` as the live portability-rebound seal and
    did not identify the later correction candidate `9a270755` or the reviewed
    documentation-durability seal `88d631e4`;
  - the replay inventory did not bind retained stdout, stderr, metrics, and
    CMake-cache artifacts for the general builds, dedicated builds, release
    and debug analytical replays, migration replays, adversarial probes,
    sequential repository harness, portable build, or unsupported
    single-precision attempt;
  - the retained evidence-local public solver deck copy and standalone empty
    portable-binary privacy scan also lacked explicit producer commands.
- Corrected choices:
  - preserve immutable historical candidate and seal labels, then require
    `CP-7` to resolve and record the exact pushed tip after this
    capture-completeness reseal because a commit cannot embed its own SHA;
  - replace implied capture with explicit output redirection, metrics paths,
    cache copies, and producer commands for every retained Phase-9 artifact;
  - retain the standalone empty binary privacy-scan artifact as an explicit
    successful no-match output rather than deleting an already sealed
    evidence surface;
  - rerun the frozen `gh-pages` overlay, bounded whitespace audit, provenance
    envelope, repository-relative outer manifest, relocated-root manifest
    verification, and final privacy audit after the correction reaches its
    final bytes.
- Alternatives rejected:
  - relying on filename conventions and prose descriptions would force the
    next agent to reconstruct evidence-capture behavior by inference;
  - deleting retained evidence surfaces solely to shorten the recipe would
    weaken continuity with the already reviewed durable packet.
- Inflection point:
  - seal and push the capture-completeness correction;
  - require fresh exact-SHA cold review of the pushed immutable tip;
  - stop again on any blocking finding before recording `RG-010`, `CP-7`, or
    `MERGE READY`.

## DR-070: Phase-9 Capture Rebound Public-Overlay Privacy Correction

- Date: 2026-05-30
- Status: documentation-only public-overlay privacy correction prepared for
  final cold review
- Capture-completeness seal reviewed:
  `1a7086add5fffd55356109b99e6a66fcd0b43486`
- Independent-review result:
  - physics and integration reviewers returned `PROCEED`;
  - evidence-durability and documentation reviewers returned `HOLD`;
  - do not record `RG-010`, `CP-7`, or `MERGE READY` on `1a7086ad`.
- Findings:
  - the portable replay wrapper and README were sealed and copied into the
    package but lacked explicit producer commands in the retained replay
    inventory;
  - the frozen public overlay copied
    `figures/cr_tracer_accuracy/cr_tracer_accuracy_summary.json`, which
    retained `54` obsolete absolute workstation paths;
  - the privacy audit checked retained `phase9_*` artifacts and the portable
    package but did not scan the complete source set copied into the public
    overlay.
- Corrected choices:
  - archive exact here-document producers for `phase9_portable_replay.py` and
    `phase9_portable_replay_README.md`;
  - normalize the obsolete public accuracy-summary worktree prefix to
    `$REPO_ROOT` without changing numerical results;
  - extend the privacy scan to `particles.md`, all copied `cr_tracer*.md`
    files, the complete `figures/cr_tracer_accuracy` directory, and both
    tracked relativistic public inputs before sealing;
  - rerun the strict frozen Pages overlay, bounded whitespace audit,
    provenance envelope, repository-relative outer manifest, relocated-root
    verification, and final privacy scan after the correction reaches its
    final bytes.
- Alternatives rejected:
  - accepting a privacy scan scoped only to retained evidence would ignore
    the actual public preview surface;
  - deleting the inherited accuracy summary from the overlay would hide a
    public artifact rather than make it portable;
  - treating package wrapper files as hand-authored exceptions would leave the
    replay recipe incomplete.
- Inflection point:
  - seal and push the public-overlay privacy correction;
  - require another fresh exact-SHA cold review;
  - stop again on any blocking finding before recording `RG-010`, `CP-7`, or
    `MERGE READY`.
