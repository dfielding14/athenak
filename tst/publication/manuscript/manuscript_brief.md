# AthenaK MHD-PIC ApJ Methods-Paper Manuscript Brief

## Control Status

- **Task type:** Manuscript planning and argument-development contract.
- **Target venue:** The Astrophysical Journal (ApJ), methods-paper framing.
- **Article subtype and formal limits:** `[AUTHOR DECISION / VERIFY AGAINST CURRENT APJ REQUIREMENTS]`.
- **Current manuscript source:** `tst/publication/manuscript/athenak_mhd_pic_apj_methods_scaffold.tex`.
- **Companion operational tracker:** `tst/publication/manuscript/claim_evidence_ledger.md`.
- **Status:** Gate 1 planning artifact with a narrowed target-paper scope. It
  does not authorize result wording, certify evidence, or replace a claim
  manifest.
- **Default reviewer state:** `pending external review`.

This brief follows
`/lustre/orion/ast207/proj-shared/dfielding/PIC/writing_guide.md`. Truth,
evidence integrity, explicit scope, and reader orientation take precedence over
fluency or apparent completeness.

## Authority Hierarchy

Use the following sources in this order when drafting or adjudicating a claim:

1. Immutable production evidence, registered run manifests, raw-artifact
   inventories, deterministic analyses, and signed reviewer dispositions.
2. `tst/publication/readiness/claims_registry.json` for stable claim IDs and
   dispositions.
3. `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md` for controlling gates,
   evidence boundaries, and authorized interpretations.
4. Source-controlled physical-model, algorithm, AMR, provenance, and toolchain
   contracts under `docs/source/engineering/`, after verification against the
   exact release source.
5. The current TeX scaffold for manuscript architecture and visible gates.
6. Historical status records and engineering proxies, which may orient future
   work but may not independently support a manuscript result.

The human pressure-selection receipt at
`/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/q011_section54_pressure_selection_receipt.json`
is authoritative only for the decision to use `problem/ps_p0=1.0` as the
provenance-first Section 5.4 baseline. It does not establish pressure
independence, physical optimality, or a qualifying shock result.

## Central Question

Under what documented physical and numerical conditions can AthenaK's core
MHD-PIC implementation support corrected linear and nonlinear Bell
calculations and a Section 5.4 parallel shock with quantitatively verified
accuracy, reproducibility, measured performance, and explicitly bounded
scientific claims?

## Proposed Main Claim

The intended methods-paper claim is:

> Within a precisely stated non-Hall MHD-PIC applicability envelope, AthenaK
> provides an auditable core implementation whose corrected Bell linear and
> nonlinear behavior, Section 5.4 parallel-shock reproduction, reproducibility,
> and measured performance are evaluated through analytical, numerical,
> registered-campaign, and independent-review gates.

This is a proposed claim, not a verified release statement. Before submission,
it must be narrowed to the claim IDs that have closed, the exact portability
matrix demonstrated, and the limitations supported by final evidence.

## Target Paper Scope

The target paper is narrowed to:

- the core AthenaK MHD-PIC physical model, discrete implementation, mesh/AMR,
  restart, provenance, and implementation-level verification needed to support
  its result claims;
- corrected Bell linear qualification followed by corrected nonlinear Bell
  qualification;
- the Section 5.4 parallel-shock campaign;
- reproducibility and evidence provenance; and
- measured performance on the exact demonstrated hardware/software matrix.

Standalone multispecies-oscillation, CRSI, CRPAI, driven expanding/compressing
box, physical-damping, and calibrated-transport results are outside this target
paper. Their open claim IDs remain open; narrowing manuscript scope does not
close, reject, or weaken their independent evidence gates.

## Intended Contribution

The paper should make the implementation and its validation logic auditable to
an astrophysical-computation reader. Its contribution is not merely that a set
of test problems runs. It should explain:

- the physical model and explicit applicability envelope;
- the discrete particle and MHD coupling design;
- how AthenaK mesh, AMR, restart, provenance, and portable execution constraints
  shape the method;
- why the validation hierarchy separates core implementation checks, corrected
  Bell qualification, Section 5.4 reproduction, nonlinear production
  qualification, and cross-code comparison;
- which physical results close, remain limited, or remain unsupported.

Contribution relative to prior methods and codes remains `[CITATION NEEDED /
VERIFY]`. Do not claim novelty, superiority, or cross-code agreement before the
relevant sources, equation maps, comparison data, and reviewer dispositions are
verified.

## Central-Result Gates

### Gate A: Nonlinear Bell Saturation

`CLAIM-PROD-BELL-NONLINEAR-NOHALL-001` is an explicit central-result gate.

The exact AthenaK deposited-current closure is
`deposited J_CR/c = PPC * deposit_qscale * species_charge * v_CR /
V_root_cell = 2 B_g k_0`, where `V_root_cell` is determined from the global
root-mesh extents and counts. Replacing `species_charge` by `q/(mc)` is valid
only for the present unit-species-mass decks. The historical Bell preparation instead targeted
`2 B_g C k_0` and treated fixed `PPC * deposit_qscale` as sufficient without
the root-cell-volume factor. Multiplication by artificial light speed `C`
changes the physical Bell mode, while fixed `PPC * deposit_qscale` makes the
physical current change with dimension or root resolution. Correct deposition
must also be invariant under a pure meshblock/rank decomposition change. This
confirmed volume-aware current-normalization defect invalidates the historical
Bell linear preparation/proxy evidence for qualification and invalidates any
nonlinear foundation or campaign design that inherits it. It does not itself
establish a corrected Bell result.

Until Q-019, Q-023, Q-025, Q-026, and Q-028 close with a named reviewer:

- no historical volume-blind or artificial-`C`-multiplied Bell artifact may
  support a linear or nonlinear result;
- no corrected nonlinear Bell campaign may qualify before the cycle-one
  deposited-current dimension/resolution/decomposition/artificial-`C` oracle
  and corrected linear predecessor close;
- no nonlinear Bell saturation result may appear in the abstract, conclusions,
  results narrative, or reader-facing figure bundle;
- no claim may transfer a non-Hall result to Hall-dominated shock-front
  conditions;
- preparation records, excluded pilots, or plausible-looking saturation may not
  be described as production qualification.

Required result scope, if the gate closes: the exact tested non-Hall physical
regime, corrected volume-aware deposited-`J_CR/c` normalization, direct-current
dimension/resolution/decomposition oracle, induction model, dimensions,
parameter matrix, uncertainty, numerical sensitivity, comparison scope, and
unresolved discrepancies.

### Gate B: Section 5.4 Parallel Shock

`CLAIM-PAPER-SHOCK-001` is an explicit central-result gate.

Until Q-003, Q-004, Q-009, Q-011, Q-016, Q-023, Q-025, Q-026, and Q-027 close
with a named reviewer:

- no Section 5.4 shock morphology, amplification, spectrum, acceleration-slope,
  AMR-agreement, load-balance, or performance conclusion may appear as a result;
- the pressure pilots must remain labeled engineering calibration;
- `problem/ps_p0=1.0` must be described only as the human-selected baseline;
- no claim may imply thermal-pool injection physics or oblique-shock generality.

Required result scope, if the gate closes: the exact paper-mode implementation,
coarse/AMR/fine paired-seed matrix, nominal-slot and observed-committed-time
semantics, retained failures, restart behavior, provenance, MPI/GPU scope,
analysis criteria, and limitations.

## Dominant Uncertainty

The dominant uncertainty is the transition from bounded local and selected
registered GPU evidence to complete long-horizon physical qualification,
broader MPI/GPU portability, nonlinear saturation, Section 5.4 shock results,
and production-scale performance.

For the two central-result gates specifically:

- **Bell nonlinear:** a source-local fixed-current-like foundation deck,
  bounded analyzer, and Q019 campaign design exist, but they inherit the
  invalid historical volume-blind, artificial-`C` current normalization and are
  ineligible for qualification. Corrected linear qualification, a
  dimension/resolution/decomposition deposited-current oracle, a superseding
  nonlinear design and foundation, qualifying generator and deck matrix,
  excluded-pilot-frozen saturation windows, comparison maps, numeric tolerances,
  measured resource model, and registered campaign remain open.
- **Section 5.4 shock:** the pressure baseline is selected and the future
  qualifying policy is detailed, and a staged resource-scaling plan now exists,
  but production registration, stable clean source, prerequisite simulations,
  measured runtime sizing, qualifying execution, and result review remain open.

## Applicability and Claim Boundaries

The manuscript must distinguish the following:

- `paper_mhd_pic` or its exact final registered runtime identity from historical
  engineering modes and separately named extensions;
- paper reproduction from AthenaK production qualification;
- corrected volume-aware Bell `J_CR/c` normalization from the invalid
  historical volume-blind, artificial-`C`-multiplied normalization and all
  artifacts that inherit it;
- non-Hall Bell behavior from any Hall extension or Hall-dominated inference;
- energetic-particle injection prescriptions from resolved thermal-pool
  injection physics;
- bounded verification from broad MPI/HIP, resilience, or performance claims;
- preparation, preregistration, and launch authorization from qualifying
  scientific results;
- AthenaK-selected release criteria from tolerances stated by a reference paper;
- Orion-only evidence retention from an institutional durable archive.

## Intended Audience and Reader Model

Primary readers are computational plasma astrophysicists, cosmic-ray transport
researchers, and developers or users of performance-portable astrophysical
simulation codes. Assume familiarity with MHD, particle methods, and numerical
validation, but do not assume familiarity with AthenaK task ordering, runtime
identities, the claim registry, or the difference between the paper-mode and
extension-mode contracts.

The manuscript should help the reader build this mental model:

1. The physical approximation defines what questions the method may answer.
2. The discrete coupling and task ordering define what must be conserved and
   verified.
3. Mesh, AMR, restart, provenance, and accelerator execution can invalidate an
   otherwise correct local operator.
4. Physical claims therefore require a hierarchy of increasingly broad evidence.
5. Bell nonlinear saturation and the Section 5.4 shock are decisive
   campaign-level demonstrations, not decorative examples.

## Conceptual Bottlenecks

- Why an MHD-PIC model is useful while remaining physically narrower than full
  kinetic plasma modeling.
- Why ideal-MHD induction and a separately named Hall extension must not be
  conflated.
- Why AthenaK's deposited current includes `1/V_root_cell`, why fixed
  `PPC * deposit_qscale` changes physical current with dimension or resolution,
  why artificial-`C` multiplication changes the physical Bell mode, and why all
  affected Bell evidence must be rerun.
- How particle interpolation/deposition, feedback, and MHD evolution compose
  into a conservative update.
- Why AMR receiver resolution, particle ownership, and mesh transitions are
  scientific-validity issues rather than implementation details.
- Why a passing engineering proxy, serial/MPI overlay, or short GPU run does not
  establish physical qualification.
- Why preregistered windows, seeds, exclusions, retained failures, and
  independent recomputation are necessary for nonlinear and shock claims.
- For Section 5.4, why the nominal output slot and actual committed simulation
  time are separate but jointly bound quantities.

## Section-Level Argument and Take-Home Messages

| Manuscript section | Argument job | Required take-home message |
|---|---|---|
| Abstract | State the verified method result and its dominant limitation | Include only closed claims; omit Bell nonlinear and shock results while their gates are open |
| Introduction | Define the scale-separation problem and why an auditable MHD-PIC implementation is needed | The paper is about a bounded method and evidence hierarchy, not implementation chronology |
| Scope, runtime identity, and claim discipline | Establish what physical model is being tested | Runtime identity and exclusions are part of every claim |
| Governing model and notation | Define equations, normalization, and approximation boundaries | The reader can identify what is evolved, coupled, omitted, and separately extended |
| Numerical method and AthenaK integration | Explain discrete updates, ordering, mesh behavior, restart, and provenance | Correctness depends on the composed algorithm, not only isolated formulas |
| Verification and qualification design | Explain evidence classes, preregistration, uncertainty, and review | Different claims require different evidence; proxies cannot close physics gates |
| Core implementation verification | Establish only the mechanics, coupling, mesh/AMR, restart, and portability evidence needed by the target claims | Each result has an oracle, scope, uncertainty, artifact chain, and limitation |
| Corrected Bell linear result | Establish the mandatory corrected predecessor to nonlinear Bell | Historical volume-blind or artificial-`C`-multiplied Bell evidence is inadmissible; corrected normalization, invariance oracles, and reruns are required |
| Nonlinear Bell result | Present a central production result only if Gate A closes | Any accepted result is restricted to the demonstrated non-Hall regime |
| Section 5.4 shock result | Present a central paper-reproduction result only if Gate B closes | Morphology, spectra, AMR comparison, and provenance must come from one registered evidence chain |
| Performance and scalability | Quantify cost and portability only from registered measurements | Performance claims require exact hardware, mapping, problem size, repeated measures, and uncertainty |
| Discussion and limitations | State what the evidence supports and where it stops | Limitations are part of the method result, not afterthoughts |
| Reproducibility and provenance | Connect every result to immutable evidence | A manuscript result must be reconstructable and reviewable |
| Conclusions | State the strongest verified contribution and next decisive test | Do not summarize intended or open claims as achieved |

## Planned Reader-Facing Figures and Tables

The current scaffold's planned bundle is retained as a proposal. Each item needs
an immutable evidence binding before release.

Priority methods items:

1. Method and task-order schematic, verified against exact release source.
2. Runtime-mode and applicability-envelope table.
3. Verification-hierarchy and simulation-matrix table.
4. Core particle mechanics, coupling, mesh/AMR, and restart/provenance
   verification figures needed by the target claims and supported by closed
   evidence.
5. Final reproducibility and claim-evidence tables.

Explicitly gated central-result items:

- Bell nonlinear saturation figure and quantitative table: prohibited until
  Gate A closes.
- Section 5.4 shock morphology, profiles, spectra, and AMR-comparison bundle:
  prohibited until Gate B closes.
- Frontier production performance figures: prohibited until registered
  performance evidence and its scope close.

Every figure and table must record raw inputs, generating command or script,
analysis revision, selections and exclusions, units, normalization, uncertainty,
claim ID, and reviewer disposition.

## Required Terminology and Notation Discipline

- Use `AthenaK` for the code and `MHD-PIC` for the bounded method.
- Name the exact runtime identity for every physical claim.
- State the full deposited-current closure, including `V_root_cell`; never
  describe a volume-blind or artificial-`C`-multiplied historical Bell artifact
  as qualifying evidence.
- Use `paper reproduction`, `production qualification`, `engineering proxy`,
  `preparation`, and `cross-code comparison` as distinct evidence classes.
- Use `non-Hall` wherever a Bell result excludes Hall induction.
- For Section 5.4, distinguish the `ideal injection surface` from the
  independently detected shock front.
- Distinguish `nominal_slot_time` from `observed_committed_time`.
- Define every symbol, normalization, sign convention, unit, and source binding
  before quantitative use.
- Do not use `state of the art` without exact compared methods, regimes,
  observables, references, and limitations.

## Allowed and Prohibited Drafting Actions

Allowed:

- Draft explanatory methods prose from verified model and implementation
  contracts, while retaining `[VERIFY]` markers where source-to-prose audit is
  incomplete.
- Draft section transitions, reader orientation, limitation statements, and
  evidence requirements.
- Insert a result only after the companion ledger points to immutable evidence,
  uncertainty, scope, and named reviewer acceptance.

Prohibited:

- Inventing citations, results, values, tolerances, novelty, comparisons, or
  performance claims.
- Treating source-controlled plans, preregistrations, or manuscript assertions
  as result evidence.
- Removing a visible gate because a result appears plausible.
- Broadening bounded evidence across dimensions, regimes, hardware, runtime
  modes, injection physics, or Hall/non-Hall boundaries.
- Editing thresholds or windows after inspecting qualifying output without a
  versioned successor and new qualifying dataset.

## Current Admissible Statements

The following statements may orient methods drafting, subject to source audit:

- AthenaK contains an MHD-PIC implementation with documented paper-mode and
  extension boundaries.
- The repository contains bounded analytical, numerical, restart, provenance,
  and selected registered GPU evidence.
- A human-only production receipt selected `problem/ps_p0=1.0` as the
  provenance-first Section 5.4 baseline.
- The current Section 5.4 pressure pilots are engineering calibration only.
- The historical Bell linear preparation omitted the root-cell-volume factor
  from its current contract, held `PPC * deposit_qscale` fixed across
  dimension/resolution changes, and multiplied the deposited `J_CR/c` closure
  by artificial `C`; it and all Bell plans or foundations inheriting that
  normalization are invalid for qualification.
- Corrected Bell linear and nonlinear campaigns require superseding contracts,
  a dimension/resolution/decomposition deposited-current oracle, complete
  reruns, immutable evidence, independent recomputation, and named review
  before they support a result.
- A Q011 preproduction resource-scaling plan recommends a three-seed paired
  core followed by resource-gated expansion; it does not modify the historical
  eight-seed preregistration or authorize production.
- Bell nonlinear saturation and the Section 5.4 shock result remain open.
- The active production policy currently authorizes no registered science
  slices.

These statements do not close a scientific claim.

## Open Author Decisions

- Confirm the ApJ article subtype, target length, supplemental-material strategy,
  and practical figure/table budget.
- Approve or narrow the proposed main claim after reviewing the final set of
  closed claim IDs.
- Preserve the narrowed target scope: core implementation, corrected Bell
  linear/nonlinear, Section 5.4 shock, reproducibility, and measured
  performance. If either central-result gate remains open, revise the target
  paper scope explicitly rather than weakening a gate or presenting an open
  result.
- Approve the final terminology for the exact runtime identity.
- Decide the final archive and data-availability statement after Q-026 review.
- Assign named scientific, methods, provenance, figure/table, and final release
  reviewers.

## Known Risks

- The source and campaign-control work remains active; prose may become stale
  unless it binds an exact release candidate.
- A large inventory of bounded evidence can obscure the central argument unless
  the manuscript preserves evidence hierarchy and scope.
- The central Bell and shock campaigns may remain blocked by prerequisites,
  resource sizing, authorization, or inconclusive physics.
- The historical Bell current-normalization defect invalidates
  inherited Bell plans and artifacts; failing to supersede every affected
  binding could silently reintroduce a current that changes with dimension or
  resolution, or the wrong physical Bell mode.
- Cross-code comparison maps and numeric tolerances are incomplete.
- Orion-only retention remains a terminal durability-review limitation.
- The manuscript can overstate portability or performance if selected GPU
  oracles are generalized beyond their registered matrix.

## Required Deliverables Before Release

- Approved manuscript brief and maintained claim-evidence ledger.
- Exact authoritative source, executable, deck, analyzer, run, and artifact
  provenance for every quantitative statement.
- Verified equations and notation sheet.
- Figure/table reconstruction records and captions with explicit limitations.
- Named reviewer dispositions for every result-bearing claim.
- Independent scientific, grounding, provenance, consistency, citation,
  figure/table, adversarial, and rendered-manuscript reviews.
- A final release-gate report listing closed, limited, rejected, and open claims.
