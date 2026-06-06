# AthenaK MHD-PIC Manuscript Claim-Evidence Ledger

## Purpose and Use

This ledger is the actionable manuscript-facing bridge between claim IDs,
immutable evidence, open work, and allowed prose. It does not replace the
machine-readable claims registry, qualification manifests, production ledger,
or named reviewer dispositions.

For every result-bearing revision:

1. Confirm the claim ID and exact intended statement.
2. Bind immutable evidence, analysis, uncertainty, scope, and exclusions.
3. Record the named reviewer disposition.
4. Change manuscript wording only to the strength supported by that disposition.
5. Preserve failed, limited, rejected, and contradictory evidence.

Allowed status values:

- `draftable_method_context`: may support bounded methods explanation, not a
  scientific result.
- `blocked_result`: result wording and result figures are prohibited.
- `limited`: a named reviewer accepted only a narrower statement.
- `rejected`: evidence does not support the intended claim.
- `qualified_pending_manuscript_review`: scientific gate closed; wording,
  consistency, and render review remain.
- `release_ready`: evidence, wording, provenance, and final review all close.
- `out_of_target_scope`: not selected for this target manuscript; the
  underlying claim remains open and no result wording is admissible.

Unless an immutable signed disposition says otherwise, every claim below remains
open or blocked.

## Current Production and Evidence Snapshot

Snapshot date: 2026-06-06. Recheck before every substantive manuscript revision.

| Item | Current evidence | Manuscript consequence |
|---|---|---|
| Active registered-science allowlist | Production policy has `registered_science_slices: []` | No current production policy authorizes a qualifying Q011 or nonlinear Bell campaign |
| Frontier accounting | Latest reconciled cumulative consumed node-hours: `1.3019444444444446`; project cap: `10000` | Budget headroom exists, but no campaign-specific node-hour sizing or authorization follows from it |
| Active science freeze | An older clean candidate is authorized by active policy | It does not bind the current timing/retention source candidate |
| Target manuscript scope | Core implementation, corrected Bell linear/nonlinear, Section 5.4 shock, reproducibility, and measured performance | Standalone oscillation, CRSI, CRPAI, driven-box, physical-damping, and calibrated-transport results are outside this target paper; their claims remain open |
| Bell current normalization | Generic deposition includes statistical/delta-f weights, physical-density scaling, shape, configured species charge, and recipient-cell volume. The uniform full-f Q043 specialization is `deposited J_CR/c = PPC * deposit_qscale * species_charge * v_CR / V_root_cell = 2 B_g k_0`; the corrected source-local matrix uses nonunit species mass and separately binds `species_charge/species_mass`. Historical preparation targeted `2 B_g C k_0` and held `PPC * deposit_qscale` fixed without the root-cell-volume factor | Historical volume-blind/artificial-`C` Bell evidence changes current with dimension/resolution and is invalid; corrected source/decks/raw-output oracle exist source-locally, but the complete registered runtime matrix and reruns remain required |
| Q043 registered-execution admission | Independent review found self-attested control-plane provenance, incomplete candidate source closure, no reconciliation-time digest inventory of raw output, and no scientific validation of cycle-zero deposited moments | The provisional admission cannot support registered or scientific evidence; a hardened fail-closed successor and negative-test matrix are required |
| Q043 registered launch preparation | A deterministic fail-closed successor generates 132 launch candidates, 132 policy candidates, exact raw-output topology interfaces, and one aggregate review fragment under a 22 node-hour and 182804545536-byte ceiling | Planning evidence only; it grants no launch, policy, science, or publication authority and remains blocked on the hardened installed reconciliation/admission path, final clean candidate, live preflight, and an empty user queue |
| Section 5.4 pressure selection | Human-only receipt selected `ps_p0_1p00`, `problem_ps_p0=1`; receipt SHA-256 `e5c1492cfc67d5cfad0110e7d772bf75f9d6d2e1fd90b0d5d8dd3338965905cc` | May state only that `p0=1.0` is the provenance-first baseline |
| Q011 future campaign policy | Three grid variants, eight paired seeds, 24 expected baseline attempts; nominal-slot/observed-time rules and all-attempt retention are frozen in a policy-only successor | Planning evidence only; no execution authorization or result |
| Q011 resource-scaling plan | Excluded preproduction pilots and a three-seed paired core with resource-gated complete-triad expansion are designed under a 500-node-hour pilot ceiling | Does not change the historical eight-seed preregistration; measured sizing and authorization remain open |
| Q011 storage projection | Approximately `10.553535598859627 TB` planning envelope for 24 baseline attempts | Planning value only; actual sizes, retry capacity, capacity confirmation, and reservation remain open |
| Q009 dynamic AMR/load-balance pilot | A reviewed fail-closed five-case contract specifies paired load-balance off/on migration, restart through `t>45`, and a held-out full-geometry engineering pilot under 100 node-hours and 200 GiB | Contract evidence only; deck materialization, trusted telemetry/admission, live policy slices, execution, and scientific-AMR qualification remain open |
| Q019 nonlinear Bell | The historical source-local foundation and analyzer inherit the invalid volume-blind/artificial-`C` normalization. A corrected volume-aware source candidate exists outside the integrated primary, but independent review found that the HR large-box pair is an exact tiled replication, the FR 3D gyroradius is under-resolved, particle noise overwhelms the intended seed, finite-rigidity linear behavior is unvalidated, and trusted raw-evidence bindings are incomplete | Historical Q019 artifacts and the reviewed candidate remain ineligible; production requires a non-tiled box-sufficiency pair, 3D resolution/PPC controls, controlled seed-versus-noise treatment, a finite-rigidity predecessor, and trusted reconciliation before any nonlinear result or figure is admissible |
| Q022 nonlinear Bell comparison | Bounded equation/normalization/parameter maps bind the retained references and establish that the target-paper Bell benchmark is linear. No standalone nonlinear-saturation dataset or numeric acceptance tolerance is supplied by that benchmark; extraction, execution bindings, discrepancy closure, and external review remain open | The bounded map may define comparison scope, but no nonlinear agreement, independent-comparison, or scoped-superiority claim is admissible |
| Primary-source literature map | Bell (2004), Bai et al. (2015), Sun & Bai (2023), Riquelme & Spitkovsky (2009), Gargaté et al. (2010), Zacharegkas et al. (2024), Stone & Gardiner (2009), Birdsall & Langdon, Mignone et al. (2018), and van Marle et al. (2018) are bound in the manuscript bibliography | Supports method lineage and bounded comparison context only; it supplies no AthenaK result, cross-code agreement, universal nonlinear-saturation tolerance, or Hall-dominated validation |
| Manuscript scaffold | Contains explicit open gates and TODOs | It is architecture, not evidence |

## Central Result Gate A: Nonlinear Bell Saturation

### Intended Claim

`CLAIM-PROD-BELL-NONLINEAR-NOHALL-001`: qualify corrected-normalization
nonlinear Bell evolution for the exact tested non-Hall AthenaK regime.

### Current Disposition

`blocked_result`

### Current Admissible Evidence

- Source and deposited-current audit confirms the exact AthenaK closure:
  `deposited J_CR/c = PPC * deposit_qscale * species_charge * v_CR /
  V_root_cell`, where `V_root_cell` is determined by global root-mesh extents
  and counts. The shorthand `species_charge=q/(mc)` applies only when the
  configured species mass is one.
- The historical Bell preparation targeted `J_CR/c = 2 B_g C k_0` rather than
  `2 B_g k_0` and held `PPC * deposit_qscale` fixed without the
  root-cell-volume factor. Artificial-`C` multiplication changes the physical
  Bell mode; the volume omission makes the current change with dimension or
  root resolution. A pure meshblock/rank decomposition change must leave the
  deposited current invariant.
- The historical Bell preparation/proxy evidence and historical Q019
  foundation/analyzer that inherit that normalization are invalid for
  qualification. A separately named corrected volume-aware Q019 campaign
  design exists, but it grants no execution or claim authority and is not
  qualifying evidence.
- Corrected Q043 source, decks, and a cycle-one raw-output oracle exist
  source-locally and use nonunit species mass, but the complete registered
  runtime matrix is unobserved and the first registered-execution admission
  failed independent provenance review.
- The bounded Q022 reference map establishes that the retained target-paper
  Bell benchmark is linear and supplies no standalone nonlinear-saturation
  dataset or numeric nonlinear acceptance tolerance.
- The mapped nonlinear Bell literature spans fixed-current, fully kinetic, and
  hybrid regimes. Comparisons are admissible only for matched dimensionless
  parameters and shared observables; no mapped source supplies a universal
  saturation amplitude or numerical acceptance tolerance.
- The production-readiness plan defines the required nonlinear saturation
  qualification contract.
- Corrected Bell linear and nonlinear campaigns do not yet exist as qualifying
  evidence.

### Prohibited Wording

- AthenaK reproduces or predicts nonlinear Bell saturation.
- Historical volume-blind or artificial-`C`-multiplied Bell artifacts validate
  the Bell mode or may seed a corrected nonlinear qualification campaign.
- AthenaK establishes a saturation amplitude, time, mechanism, morphology, or
  energy partition.
- A non-Hall result applies to Hall-dominated shock-front conditions.
- A mapped nonlinear study supplies a universal saturation tolerance across
  fixed-current, MHD-PIC, hybrid-PIC, and fully kinetic regimes.
- AthenaK agrees with, improves upon, or is competitive with another method or
  code for nonlinear Bell behavior.

### Required Evidence Package

| Required item | Acceptance condition | Current status | Next action |
|---|---|---|---|
| Corrected Bell current normalization | Source, deck, analyzer, and direct cycle-one deposited-current oracle bind `PPC * deposit_qscale * species_charge * v_CR / V_root_cell = 2 B_g k_0`; reject artificial-`C` multiplication; and demonstrate invariant physical current across the registered dimension, root-resolution, meshblock/rank-decomposition, and artificial-`C` matrix | Corrected source-local source/decks/oracle exist and focused tests plus one serial cycle-one smoke pass; complete registered matrix and hardened provenance admission remain absent | Harden admission, execute the complete registered raw-output matrix, retain immutable inventories, and obtain independent recomputation/review |
| Q003/Q004/Q005 prerequisites | Required mechanics, coupling, and corrected paper-faithful linear Bell gates close | Open; historical Bell linear evidence invalidated | Complete the corrected qualifying linear Bell convergence, MPI/GPU, independent-recompute, and review matrix |
| Exact physical regime | Corrected fiducial normalization, parameters, dimensional scope, exclusions, and artifact root frozen before qualifying output | Corrected volume-aware successor design exists but is not execution-frozen; corrected implementation and evidence are absent | Implement, review, pilot, and freeze the separately named corrected Q019 successor before registration |
| Qualifying decks | Exact corrected deck matrix and checksums frozen | Historical source-local foundation invalid for qualification; corrected successor design exists but corrected deck matrix is absent | Implement and review the separately named corrected Q019 source and exact pilot matrix |
| Nonlinear analyzer | Production analyzer, source digest, outputs, failure behavior, and independent recompute contract frozen | Bounded source-local analyzer exists; production and independent-recompute contracts absent | Extend only after required raw diagnostics and pilot-frozen criteria close |
| Saturation windows and criteria | Theory-derived or excluded-pilot-derived windows, estimators, uncertainty, tolerances, outlier rule, and stop conditions frozen before qualifying inspection | Open | Run only explicitly excluded pilots, then freeze successor |
| Sensitivity design | Particle count, timestep, resolution, box size, dimensionality, and seed requirements bound; any sequential design has maximum seeds and node-hour ceiling | Proposed Q019 staged design; not frozen for execution | Complete independent physics/resource review and excluded pilots |
| Resource model | Nodes, walltime, memory, output size, node-hour ceiling, retry ceiling, and storage envelope measured and reviewed | Planning ceiling exists; measured model absent | Run excluded sizing pilots and produce estimator |
| Registered execution | Clean candidate, executable, decks, analyzer, policy slices, launch contracts, and retention bound | Unauthorized | Install/promote reviewed successor and execute serially |
| Raw evidence | Every qualifying, failed, outlier, and replacement attempt archived with inventories | Absent | Enforce registered retention during execution |
| Primary observables | Amplification, dominant wavelength, spectra, morphology, energy partition, saturation amplitude, and saturation time measured with uncertainty | Absent | Analyze only after complete registered dataset |
| Q022/Q028 comparison | Reference extraction, equation/normalization/parameter maps, defensible comparison criteria, discrepancy ledger, and comparison runs close | Bounded central/reference-specific maps complete; target-paper Bell benchmark is linear and supplies no standalone nonlinear numeric tolerance; execution comparison and external review remain open | Bind corrected execution evidence, retain the explicit no-standalone-nonlinear-tolerance boundary, and close the discrepancy ledger and review |
| Independent recomputation | Reviewer-owned implementation regenerates the primary metric table from archived raw output | Absent | Produce independent artifact after execution |
| Named review | Reviewer accepts exact claim scope and limitations | Pending external review | Assign reviewer after evidence package is complete |

### Required Manuscript Insertions After Closure

- Exact regime and non-Hall applicability statement.
- Exact volume-aware deposited-`J_CR/c` normalization, direct-current
  dimension/resolution/decomposition/artificial-`C` invariance results, and
  explicit inventory of superseded historical artifacts.
- Registered simulation matrix and resource/provenance table.
- Linear-to-nonlinear transition and saturation metrics with uncertainty.
- Sensitivity trends without post-inspection selection.
- Figure showing only reviewer-accepted amplification, spectra, morphology,
  energy partition, and saturation observables.
- Independent-comparison table and discrepancy ledger.
- Explicit statement that the result cannot be transferred to Hall-dominated
  shock-front conditions.

### Fail-Closed Outcomes

If the expected asymptotic regime is not reached, sensitivities do not close,
the result leaves its preregistered validity envelope, or comparisons expose an
unresolved material discrepancy, record `limited`, `rejected`, or
`blocked_result`. Do not tune windows, drop seeds, or broaden tolerances after
inspection. Any use of a historical volume-blind or artificial-`C`-multiplied
Bell artifact keeps the corrected linear and nonlinear claims blocked.

## Central Result Gate B: Section 5.4 Parallel Shock

### Intended Claim

`CLAIM-PAPER-SHOCK-001`: reproduce the specified Section 5.4 non-relativistic
parallel-shock morphology, magnetic amplification, downstream spectra, and
coarse/AMR/fine comparison within the exact paper-mode and injection
prescription.

### Current Disposition

`blocked_result`

### Current Admissible Evidence

- A human-only immutable receipt selected `problem/ps_p0=1.0` as the
  provenance-first baseline. Bai et al. (2015) explicitly specify `P0=T0=1`
  and motivate its irrelevance only while `T0 << v0^2`; Sun & Bai (2023)
  Section 5.4 does not state `P0`.
- Pressure pilots and their analysis are engineering calibration only.
- Source-local preparation, bounded injection/provenance diagnostics, restart
  guards, and future campaign policies may support methods and planned-analysis
  explanation.
- A preproduction scaling design defines excluded strong-scaling,
  reduced-transverse full-time, and held-out validation pilots under a hard
  500-node-hour ceiling.
- The future policy freezes three variants, eight paired seeds, 24 expected
  baseline attempts, nominal output slots, observed committed times, retained
  raw products, and all-attempt failure retention.

### Prohibited Wording

- `p0=1.0` is pressure-independent, physically preferred, or demonstrated to be
  optimal.
- AthenaK reproduces the Section 5.4 shock result.
- AthenaK establishes shock morphology, magnetic amplification, spectra,
  acceleration slope, AMR agreement, load-balance benefit, or production
  performance.
- The injection prescription resolves thermal-pool injection.
- The result generalizes to oblique shocks.

### Required Evidence Package

| Required item | Acceptance condition | Current status | Next action |
|---|---|---|---|
| Stable source baseline | Curated clean commit equals reviewed remote tip; focused and control-plane regressions pass against exact source | Open | Freeze concurrent work, commit, push, and validate |
| Q003/Q004/Q009 prerequisites | Pusher, coupling, multi-rank AMR migration/restart, HIP, and scientific-AMR gates close for required scope | Open | Complete prerequisite registered simulations |
| Pressure decision binding | Final campaign binds the immutable human receipt and selected case exactly | Receipt complete; campaign binding absent | Bind receipt during immutable campaign materialization |
| Campaign registration | Exact clean commit, executable, per-variant decks, analyzers, Orion root, policy slices, and independent recompute plan bound | Unauthorized | Materialize, review, install, and promote successor |
| Resource sizing | Per-variant nodes, walltime, throughput, memory, node-hours, output/checkpoint size, retry ceiling, and capacity confirmation frozen | Staged preproduction design exists; measurements absent | Run registered excluded coarse/fine/AMR timing and I/O pilots |
| Short qualification matrix | Short coarse, fine, AMR, and restart-AMR pilots pass and remain distinct from production evidence | Open | Register and execute eligible short pilots |
| Full qualifying matrix | Coarse, three-level AMR, and fine variants complete for eight paired seeds under registered `normal`/`batch` policy | Absent | Execute serially only after prerequisites close |
| Time semantics | Each required nominal slot binds one canonical observed committed time and common cycle; terminal `t=1200` requirements pass | Policy frozen; runtime evidence absent | Validate every attempt manifest and raw inventory |
| Retention | Every raw mesh, particle, restart, completion, failure, outlier, and replacement artifact retained and inventoried | Policy frozen; runtime evidence absent | Enforce retention and size retry contingency |
| Morphology and amplification | Preregistered `t=500` products and ideal-surface/detected-front distinction pass review | Absent | Generate from complete registered raw data |
| Spectra and slope | Downstream spectra at required slots, fixed bins, exclusions, overflow accounting, uncertainty, and late-tail fit pass | Absent | Generate from provenance-bound particles |
| AMR comparison | Paired-seed AMR versus fine residuals and honest speed/memory comparison pass | Absent | Complete paired analysis and uncertainty |
| Restart and provenance | Registered restart carrier, particle provenance, cohort spectra, MPI/GPU, load balance, runtime, and memory evidence close | Open | Execute and review bound evidence |
| Independent recomputation | Reviewer-owned implementation regenerates primary metrics from archived raw output | Absent | Produce after full dataset is immutable |
| Named review | Reviewer accepts exact reproduction claim and limitations | Pending external review | Assign reviewer after evidence package is complete |

### Required Manuscript Insertions After Closure

- Exact campaign registration, source/executable/deck/analyzer bindings, and
  selected pressure receipt.
- Table of coarse/AMR/fine paired-seed runs, failures, replacements, resources,
  and retained evidence.
- Morphology and profile figure at preregistered slots, distinguishing ideal
  injection surface from detected front.
- Upstream magnetic-amplification and downstream-spectrum figures with
  uncertainty, exclusions, overflow accounting, and fit details.
- Paired AMR-versus-fine residual table and a carefully scoped interpretation.
- Restart, provenance, MPI/GPU, load-balance, runtime, memory, and storage
  evidence for the same qualifying campaign.
- Explicit limitations: simplified injection, no thermal-pool injection claim,
  no oblique-shock generality, and no pressure-independence claim.

### Fail-Closed Outcomes

Any missing required slot, ambiguous committed time, failed terminal completion,
missing raw artifact, post-inspection threshold change, unregistered
replacement, unresolved prerequisite, or failed independent recomputation keeps
the claim blocked or causes a limited/rejected disposition.

## Supporting Claim Ledger

These claims support the narrowed target methods-paper argument but do not
substitute for either central-result gate.

| Claim ID or topic | Intended manuscript role | Current admissible evidence | Required closure before result wording | Status |
|---|---|---|---|---|
| `CLAIM-PAPER-GYRO-001` | Support core particle-mechanics verification; not a standalone physical centerpiece | Bounded analytical convergence and selected registered GPU oracle | Q-003, Q-025, Q-026; final figure/table, uncertainty, provenance, named review | `blocked_result` |
| Paper coupling methods evidence | Explain opposite-impulse exchange on its qualified uniform-mesh scope and ideal-MHD induction isolation | Bounded host and selected registered GPU evidence described in production plan | Exact source audit, full intended decomposition/GPU scope, explicit boundary/AMR/injection/removal exclusions, and claim-specific review | `draftable_method_context` |
| `CLAIM-PAPER-BELL-LINEAR-001` | Establish the mandatory corrected linear predecessor to nonlinear Bell | Historical preparation/proxy evidence is invalid because it omitted `V_root_cell`, fixed `PPC * deposit_qscale` across dimension/resolution changes, and multiplied deposited `J_CR/c` by artificial `C` | Corrected volume-aware `J_CR/c` source/decks/analyzer; dimension/resolution/decomposition oracle; Q-003, Q-004, Q-005, Q-023, Q-025, Q-026; complete registered thin-2D/3V one-dimensional-variation plus 2D/3D rerun, independent recompute, and named review | `blocked_result` |
| Performance and scalability | State measured cost, memory, communication, and load balance | Bounded observability scaffolding only for broad claims | Registered measurements, exact machine/runtime mapping, repeated measures, uncertainty, and review | `blocked_result` |
| `CLAIM-RELEASE-PAPER-MHD-PIC-001` | Final production-ready bounded paper-mode release claim | Individual bounded gates and selected registered slices only | Final selected profile gates, Q-014 terminal sign-off, and release review | `blocked_result` |

## Claims Outside Target Paper Scope

These claims remain open in the claim registry but are not planned as
standalone results in this target paper. They must not re-enter the manuscript
without an explicit scope revision and their original evidence gates.

| Claim ID or topic | Target-manuscript disposition | Claim consequence |
|---|---|---|
| `CLAIM-PAPER-OSCILLATION-001` | `out_of_target_scope` | Remains open; no standalone oscillation result or figure |
| `CLAIM-PAPER-CRSI-LINEAR-001` | `out_of_target_scope` | Remains open; no CRSI result or figure |
| `CLAIM-PAPER-CRPAI-LINEAR-001` | `out_of_target_scope` | Remains open; no CRPAI result or figure |
| `CLAIM-PAPER-CRPAI-DRIVEN-001` | `out_of_target_scope` | Remains open; no driven expanding/compressing-box result or figure |
| Physical-damping and calibrated-transport extensions | `out_of_target_scope` | Remain separately gated; no transport-calibration result or figure |

## Figure and Table Evidence Actions

| Planned manuscript object | Claim binding | Current status | Evidence required before insertion |
|---|---|---|---|
| Method and task-order schematic | Methods context | Open | Exact release-source trace and methods reviewer approval |
| Runtime-mode/applicability table | Methods context | Open | Verified runtime identities, equations, exclusions, and source bindings |
| Verification hierarchy/simulation matrix | Multiple claims | Open | Reconcile every row with final claim disposition and immutable registrations |
| Corrected Bell linear dispersion figure | `CLAIM-PAPER-BELL-LINEAR-001` | Blocked | Corrected volume-aware `J_CR/c` contract, dimension/resolution/decomposition oracle, complete registered qualifying rerun matrix, analytical comparison, uncertainty, provenance, independent recompute, and review |
| Bell nonlinear saturation bundle | `CLAIM-PROD-BELL-NONLINEAR-NOHALL-001` | Blocked central result | Full Gate A evidence package |
| Section 5.4 shock bundle | `CLAIM-PAPER-SHOCK-001` | Blocked central result | Full Gate B evidence package |
| Frontier performance/scaling bundle | Performance claims | Blocked | Registered measurement campaign and reviewed scope |
| Final claim-evidence table | All result-bearing claims | Open | Final dispositions, exact artifact links, limitations, and named reviewers |
| Reproducibility table | All result-bearing claims | Open | Exact commits, binaries, decks, analyzers, commands, inventories, and archive disposition |

For each inserted object, record:

- manuscript label and caption;
- claim ID;
- exact immutable raw-artifact root and inventory digest;
- generating script and digest;
- command or workflow needed to reconstruct it;
- selections, exclusions, units, normalization, and uncertainty;
- reviewer identity and disposition;
- supported conclusion and explicit unsupported inference.

## Ordered Manuscript-Evidence Critical Path

1. Freeze the exact release-source candidate and rerun the full required
   verification/control-plane suites.
2. Close Q003/Q004 and the campaign-specific prerequisites shared by Bell and
   shock work.
3. Supersede every historical volume-blind/artificial-`C` Bell binding, close
   the dimension/resolution/decomposition deposited-current oracle, then close
   corrected paper-faithful linear Bell Q005 before any nonlinear Bell
   qualifying campaign.
4. Close multi-rank/HIP/scientific-AMR Q009 requirements before the Q011 full
   shock campaign.
5. Freeze campaign-specific resource models, retry ceilings, storage envelopes,
   decks, analyzers, statistical designs, and independent-recompute plans.
6. Install and promote reviewed control-plane/policy successors with exact
   registered-science slices.
7. Execute, reconcile, archive, and independently recompute the qualifying
   campaigns without selecting away failures or outliers.
8. Complete Q022/Q028 Bell comparison work and record discrepancies.
9. Obtain named scientific and provenance reviewer dispositions.
10. Insert only accepted results, figures, tables, limitations, and provenance
    into the TeX manuscript.
11. Run scientific-fidelity, argument, pedagogy, citation, equation,
    figure/table, adversarial, consistency, and rendered-object reviews.

## Open Questions

| Question | Classification | Required disposition |
|---|---|---|
| What exact corrected-normalization non-Hall Bell regime and dimensional scope should the Q019 successor target? | Blocks Bell claim | Freeze before qualifying execution |
| Does the corrected Q019 box-sufficiency comparison use a genuinely non-tiled matched realization with controlled access to modes below the small-box fundamental? | Blocks Bell claim | Demonstrate translation-symmetry rejection and preregister the paired spectral comparison |
| Does the finite-rigidity Q019 branch resolve the evolving gyroradius and distinguish controlled deterministic seeding from particle-noise-seeded evolution? | Blocks Bell claim | Close a finite-rigidity early-time predecessor and freeze 3D resolution, PPC, quiet-start/noise, and early-growth gates |
| Has every historical Bell source, deck, analyzer, foundation, design, and artifact inheriting the volume omission, fixed `PPC * deposit_qscale`, or artificial-`C` multiplication been inventoried and excluded? | Blocks Bell claims | Complete supersession audit before corrected execution |
| Does the direct deposited-current oracle hold the physical current fixed across dimension, root resolution, and meshblock/rank decomposition? | Blocks Bell claims | Close before corrected Bell execution |
| Does the hardened Q043 admission bind an independently trusted control-plane event, the actual built executable and deposition/output source closure, reconciliation-time raw-output digests, and valid cycle-zero moments? | Blocks Bell claims | Close before any Q043 registered evidence is admissible |
| What excluded pilots are sufficient to freeze Bell saturation windows and resource sizing? | Blocks Bell claim | Review pilot design before launch |
| What node, walltime, retry, node-hour, and storage ceilings should Q011 register? | Blocks shock execution | Measure and approve before policy promotion |
| Does every full Q011 attempt fit one compliant `normal` allocation? | Blocks shock execution | Demonstrate by reviewed timing/sizing evidence |
| Can the production AMR particle-moment policy close a per-particle partition-of-unity conservation oracle without losing the intended receiver-resolution smoothness, or must exact conservation remain limited to uniform meshes? | Blocks the scope of shock conservation claims | Evaluate a separately named conservative policy and retain `paper_smooth` as a distinct, explicitly non-exact reproduction option |
| What is the final exact release runtime identity and portability scope? | Blocks broad main claim | Resolve from final closed gates |
| What final archive/data-availability statement is supportable under Orion-only retention? | Blocks release wording | Resolve through Q-026/external review |
| Which additional prior-method citations are required for AMR, reproducibility, GPU portability, and the relativistic Boris method? | Requires verification | Extend the verified bibliography without inventing citations |

## Evidence Acceptance Record Template

Use one copy of this block for each proposed result before manuscript insertion.

```text
Claim ID:
Proposed manuscript statement:
Exact scope and exclusions:
Evidence class:
Immutable run/manifest paths:
Raw-artifact inventory digests:
Analysis and figure/table script digests:
Uncertainty and tolerance record:
Independent recomputation artifact:
Contradictory, failed, and outlier evidence:
Named reviewer:
Reviewer disposition:
Allowed wording:
Prohibited broader inference:
Manuscript locations:
Final consistency/render review:
```
