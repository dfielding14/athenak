# AthenaK MHD-PIC Production Readiness and Sun & Bai (2023) Reproduction Plan

## Document Status

| Item | Value |
| --- | --- |
| Purpose | Canonical implementation, verification, publication-reproduction, and Frontier qualification plan for the AthenaK MHD-PIC module |
| Primary reference | Sun & Bai, *The Magnetohydrodynamic-Particle-In-Cell Module in Athena++: Implementation and Code Tests*, source manuscript in `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` |
| Initial review baseline | Historical local `c/pic-review` checkout at base commit `ca7e43f0d0d6`, reviewed with an extensive dirty working tree containing in-progress PIC changes |
| Plan-revision baseline | Clean `PIC` branch at commit `1f7aeff19e55a909d2c8dd9c26d2eafe09673ef6`, tracking `origin/PIC` |
| Production verdict at review time | **Not production-ready. Paper reproduction is blocked by model-level implementation gaps.** |
| Supersedes | All scope, terminology, pass criteria, operational instructions and sign-off language in `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md`; that retired note is historical evidence only |
| Required execution root on Frontier | `/lustre/orion/ast207/proj-shared/dfielding/PIC` |
| Frontier scheduling constraint | Use only the `debug` QOS, one job in any state at a time, no job longer than 2 hours |
| Total Frontier testing budget | 1000 node-hours maximum, tracked before and after every submitted job |

This plan is written as a production gate, not as a list of desirable enhancements.
Agents must not declare the MHD-PIC implementation production-ready, nor claim
reproduction of Sun & Bai (2023), until every blocking gate below is satisfied
with archived evidence.

## Executive Verdict

The current implementation contains valuable scaffolding, tests, diagnostics,
publication tooling, and a shock problem generator. It is not yet an
implementation of the physical and numerical model required by the reference
paper. Four defects are immediately disqualifying:

1. The coupled path injects deposited particle current into the CT electric
   field in `src/mhd/mhd_tasks.cpp`, while the paper explicitly uses ideal-MHD
   induction, `cE = -u x B`, and neglects the CR-induced Hall term. Several
   current legacy paper-like input decks enable this field modification while leaving
   conservative momentum and energy feedback disabled.
2. The particle pusher advances non-relativistic velocities with
   `0.5*m*v^2`, while the paper advances mass-normalized momentum with an
   artificial light speed `C`, relativistic `gamma`, and relativistic energy.
   The current parser rejects non-default `pic_cr_light_speed`, although the
   paper's validation cases require large, case-specific values of `C`.
3. `pic_deltaf_mode=on` is currently a deterministic quiet-start option, not
   the paper's evolving delta-f weight method. It cannot support physical
   reproduction of the CR streaming or pressure-anisotropy experiments.
4. The expanding/compressing-box hook applies particle velocity scaling only.
   The paper requires a coupled MHD and particle coordinate transformation with
   gas and magnetic-field source terms and analytic validation.

There are additional high-risk reliability issues: AMR reconstructs child
MeshBlock and coordinate objects while retaining physics modules and selected
boundary helpers, but the lifetime and buffer-refresh contract is not proved;
the load-balancing path does not yet demonstrate particle-cost weighting; and
the refinement-boundary deposition policy must be made explicit because the
paper knowingly trades individual conservation for smooth feedback at
resolution interfaces. Restart output also lacks a qualified atomic-publication,
completion-marker and last-known-good recovery contract, which blocks
production operation even after the physical-model defects are repaired.

The correct next step is not a larger physics run. The correct next step is to
repair the model contract and establish stringent small analytical tests. Only
then should Frontier be used for GPU, AMR, scaling, and publication simulations.

## Non-Negotiable Rules For Future Agents

1. Treat this document as a living release gate. Do not silently reduce a
   tolerance, omit a failed diagnostic, or reinterpret a test as passing.
2. Prefer the paper-faithful model for reproduction. Any alternate MHD-PIC
   formulation must have a separate runtime mode, derivation, documentation,
   and independent validation; it must never be conflated with the paper mode.
3. Use only the canonical machine-readable claim classes defined in the claims
   registry below. In particular, keep `unit/regression`, `engineering_proxy`,
   `physics_validation`, `sun_bai_2023_reproduction`,
   `athenak_production_mode`, `cross_code_comparison`,
   `scoped_state_of_the_art`, and `unsupported` separate.
4. Existing proxy tests are useful guardrails only. A proxy that grows, remains
   stable, or has serial/MPI parity is not a substitute for the paper's
   dispersion relation, conservation budget, spectrum, or shock diagnostic.
5. Any current dirty-tree changes must be reviewed and curated into small,
   justified commits before being used as a production baseline.
6. On Frontier, all simulation, build, log, metric, ledger, and artifact paths
   must remain below `/lustre/orion/ast207/proj-shared/dfielding/PIC`.
7. On Frontier, submit only `debug` jobs. Never submit a second job while any
   previous `debug` job is queued or running, never request more than 2 hours,
   and never exceed the cumulative 1000 node-hour budget.
8. Treat `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md` as retired historical
   evidence only. It must never authorize a run, close a gate, or define a
   scientific claim.
9. Do not use the phrase `state of the art` without naming the compared methods,
   codes, physical regimes, metrics, reference versions and limitations. A
   scoped comparative conclusion is required instead.
10. Do not start a production-scale Frontier campaign under the `debug` QOS.
    The currently authorized Frontier phase is short correctness and
    qualification work only. Any future production campaign requires an
    explicit plan revision and execution authorization.

## Claims Registry And Evidence Boundaries

This registry is the controlling interpretation layer for all implementation,
testing, publication and release work. A passing test supports only the claim
class explicitly assigned to it. Future agents must update this registry when a
new runtime mode, algorithm, experiment or external comparison is added.

| Claim class | Allowed claim | Mandatory evidence before use | Explicitly insufficient evidence |
| --- | --- | --- | --- |
| `unit/regression` | A narrow implementation invariant remains intact | Deterministic test, stated invariant, failure artifact and supported configurations | Successful execution, trend sign or serial/MPI parity alone |
| `engineering_proxy` | A development workflow or broad behavior trend is useful for detecting regressions | Proxy label in test name, manifest and figure; documented limitation; no promotion into a physics result table | Positive growth, stable output, finite values, clean restart or plot generation |
| `physics_validation` | A documented equation or algorithm is accurate in a stated regime | Analytical, manufactured or frozen external oracle; convergence; tolerance rationale; CPU/MPI/GPU scope; raw artifacts | Qualitative similarity or a single-resolution overlay |
| `sun_bai_2023_reproduction` | A named Sun and Bai (2023) benchmark or figure is quantitatively reproduced | Paper-faithful mode; frozen paper parameters; deterministic analysis; uncertainty-aware metrics; archived raw output and figure reconstruction | Paper-like input names, proxies, or results from an extension mode |
| `athenak_production_mode` | A documented AthenaK mode is reliable for its stated applicability envelope | Closed physical oracles; conservation; AMR/restart/boundary/decomposition tests; portability; resilience; usability; performance evidence; limitations | Paper reproduction alone |
| `cross_code_comparison` | AthenaK agrees with an independent code or published result for named observables and regime | Frozen reference version/data; mapping of equations and normalization; uncertainty-aware comparison; discrepancy ledger | Shared ancestry, copied kernels, visual resemblance or unmatched physics |
| `scoped_state_of_the_art` | AthenaK matches or extends named contemporary capabilities for a precisely bounded use case | Closed production-mode gate plus independent literature/code comparison matrix, nonlinear saturation evidence where claimed, performance measurements and explicit exclusions | An unqualified global claim that the implementation is `state of the art` |
| `unsupported` | A capability is not qualified | Explicit parser rejection where unsafe, documentation, and known-limitations entry | Silent fallback, inert flag or ambiguous mode name |

### Claim Registration Procedure

Before any new result is described as validation, reproduction, production or
state-of-the-art evidence:

1. Add a stable claim ID such as `CLAIM-PAPER-BELL-LINEAR-001`.
2. State the runtime mode, equations, dimensionality, boundary conditions,
   parameter range and intentionally excluded physics.
3. Link the exact test or campaign IDs, metrics, predeclared tolerances,
   uncertainty method, raw artifact root and analysis revision.
4. Record whether the evidence is analytical, manufactured, paper reproduction,
   external published data, or a frozen independent-code comparison.
5. Name the claim reviewer and disposition: `open`, `qualified`, `limited`,
   `superseded`, or `rejected`.
6. Reject any wording broader than the archived evidence.

### Initial Claim Instance Register

This table is intentionally populated before implementation work resumes. Add
rows rather than silently broadening an existing claim. Store immutable evidence
under `$PIC_ROOT/manifests/claims/<claim_id>/` and link compact repository-side
reports from the final sign-off bundle.

Each row inherits a required manifest pointer of
`$PIC_ROOT/manifests/claims/<claim_id>/claim.json`. That manifest must expose the
registered test and campaign IDs, tolerances, artifact roots, analysis
revision, evidence links and reviewer. The initial reviewer is `TBD`; no row
may move to `qualified` until a named reviewer signs the manifest.

| Claim ID | Intended claim and regime | Runtime mode and exclusions | Required gates | Initial disposition |
| --- | --- | --- | --- | --- |
| `CLAIM-PAPER-GYRO-001` | Reproduce Sun and Bai Section 5.1 relativistic gyro-motion | `paper_mhd_pic` particle mechanics; no extension substitutions | Q-003, Q-025, Q-026 | `open` |
| `CLAIM-PAPER-BELL-LINEAR-001` | Reproduce Sun and Bai Section 5.2 Bell real and imaginary dispersion in 1D/2D/3D | `paper_mhd_pic`; excludes Hall-extension claims | Q-003, Q-004, Q-005, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-PAPER-OSCILLATION-001` | Reproduce Sun and Bai Section 5.3 oscillation on uniform, SMR and AMR grids | `paper_mhd_pic`; requires audited AMR policy | Q-004, Q-006, Q-009, Q-025, Q-026 | `open` |
| `CLAIM-PAPER-SHOCK-001` | Reproduce Sun and Bai Section 5.4 parallel-shock results | `paper_mhd_pic`; excludes thermal-pool injection and oblique-shock generality | Q-003, Q-004, Q-009, Q-011, Q-016, Q-023, Q-025, Q-026, Q-027 | `open`; full run may remain authorization-blocked |
| `CLAIM-PAPER-CRSI-LINEAR-001` | Reproduce Sun and Bai Section 5.5 CRSI dispersion and spectra | `paper_mhd_pic` with true delta-f | Q-003, Q-004, Q-007, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-PAPER-CRPAI-LINEAR-001` | Reproduce Sun and Bai Section 5.6 prolate/oblate CRPAI dispersion and polarization | `paper_mhd_pic` with true delta-f | Q-003, Q-004, Q-007, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-PAPER-CRPAI-DRIVEN-001` | Reproduce Sun and Bai Section 5.7 driven-box CRPAI evolution | `paper_mhd_pic` with true delta-f and full box equations; numerical-damping balance is not a calibrated transport coefficient | Q-007, Q-008, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-PROD-BELL-NONLINEAR-NOHALL-001` | Qualify nonlinear Bell evolution for the tested non-Hall AthenaK regime | `paper_mhd_pic`; explicitly excludes Hall-dominated shock-front inference | Q-019, Q-023, Q-025, Q-026, Q-028 | `open` |
| `CLAIM-PROD-CRSI-SATURATION-001` | Qualify bounded nonlinear CRSI growth, quasi-linear evolution and saturation | `paper_mhd_pic` with true delta-f; physical damping is separate | Q-020, Q-023, Q-025, Q-026, Q-031 | `open` |
| `CLAIM-PROD-CRPAI-SATURATION-001` | Qualify bounded driven CRPAI saturation behavior | `paper_mhd_pic` with true delta-f and full box equations; no calibrated transport claim | Q-021, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-PROD-CRPAI-UNDRIVEN-SATURATION-001` | Qualify bounded undriven CRPAI nonlinear evolution and saturation behavior | `paper_mhd_pic` with true delta-f; separate from Section 5.6 linear reproduction and from physical-damping transport calibration | Q-041, Q-023, Q-025, Q-026 | `open` |
| `CLAIM-EXT-HALL-BELL-001` | Qualify a CR-induced Hall Bell extension if implemented | Separately named `extended_mhd_pic` mode only | Q-029 | `unsupported` until implemented |
| `CLAIM-EXT-CRSI-IN-DAMPING-001` | Qualify ion-neutral-damped CRSI if implemented | Separately named extension; excludes ordinary paper-mode claim | Q-032 | `unsupported` until implemented |
| `CLAIM-EXT-CRPAI-TRANSPORT-001` | Calibrate CRPAI saturated-state effective scattering with physical damping if implemented | Separately named adaptive-delta-f and physical-damping extension | Q-033 | `unsupported` until implemented |
| `CLAIM-XCODE-ATHENA-BELL-SHOCK-001` | Compare matched Bell and shock observables against the frozen Athena/Bai reference family | Matched observable subset; Hall-sensitive regimes must be separated | Q-022, Q-028, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-BELL-SATURATION-001` | Compare nonlinear Bell saturation against independent kinetic and hybrid-PIC literature | Bounded parameter overlap; report model differences explicitly | Q-022, Q-028, Q-040 | `open` |
| `CLAIM-XCODE-BAI2019-CRSI-001` | Compare nonlinear CRSI against Bai et al. (2019) | `paper_mhd_pic` true delta-f overlap only | Q-022, Q-031, Q-040 | `open` |
| `CLAIM-XCODE-PLUTO-COUPLING-001` | Compare matched conservative-coupling behavior against PLUTO | Matched-equation subset only | Q-022, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-AMRVAC-SHOCK-001` | Compare overlapping AMR-shock observables against MPI-AMRVAC publications | Matched observable subset only | Q-022, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-GIZMO-RSOL-DECISION-001` | Decide whether a bounded Ji-Hopkins/GIZMO RSOL comparison is informative, or exclude it with rationale | Do not conflate GIZMO RSOL equations with the Sun and Bai artificial-`C` contract | Q-039 | `open` |
| `CLAIM-XCODE-ENTITY-MICRO-001` | Compare bounded shared particle-numerics micro-oracles against frozen Entity Toolkit source | Shared-kernel pusher/deposition/restart subset only; excludes EM evolution and gas feedback | Q-015, Q-022 | `open` |
| `CLAIM-STATEART-CRPAI-SCATTERING-001` | Make a scoped contemporary CRPAI scattering comparison if AthenaK implements the required extension | Requires physical damping, adaptive delta-f comparison and explicit parameter envelope | Q-033, Q-035 | `unsupported` until implemented |
| `CLAIM-BUNDLE-SUN-BAI-2023-001` | Qualify the complete, paper-faithful Sun and Bai (2023) reproduction bundle | Exactly the `sun_bai_2023_reproduction` profile; no extension substitutions and no production-readiness implication | All `sun_bai_2023_reproduction` profile gates, then terminal Q-014 disposition recorded in the signed final manifest | `open` |
| `CLAIM-RELEASE-PAPER-MHD-PIC-001` | Qualify the bounded production-ready AthenaK `paper_mhd_pic` release | Exactly the documented paper-mode applicability envelope; unsupported extensions remain excluded | All Production-ready `paper_mhd_pic` profile gates, then terminal Q-014 disposition recorded in the signed final manifest | `open` |

### Explicit Unsupported-Capability Register

Unsupported claims must remain visible until a separately named implementation
and its qualification gates are complete.

| Capability or claim | Handling | Closure artifact |
| --- | --- | --- |
| Full electromagnetic PIC physics | `rename_proxy`, `docs_limit` | Reclassified manifest, figure watermark and known-limitations entry; current EM-vacuum, two-stream and Weibel names must not imply this capability |
| Self-consistent Langmuir-wave physics | `rename_proxy`, `docs_limit` | Rename the uniform-`B` orbit-frequency anchor and document its narrow oracle |
| Physical two-stream or Weibel fidelity | `rename_proxy`, `docs_limit` | Reclassified engineering-proxy manifest and figures |
| Hall-dominated Bell or shock-front behavior | `extension_gate` | Explicit unsupported disposition unless `CLAIM-EXT-HALL-BELL-001` is implemented and qualified |
| Injection from the thermal pool | `docs_limit` | Section 5.4 reproduction manifest documents its simplified prescription |
| Relativistic MHD background fluid | `parser_reject`, `docs_limit` | Parser test and supported-mode table |
| Physical damping and calibrated CR transport coefficients | `extension_gate` | Explicit unsupported disposition unless separately implemented and qualified |
| Oblique-shock generality | `parser_reject`, `docs_limit` | Parser test or unsupported-mode documentation until independently implemented |
| Production-scale Frontier qualification under the current `debug`-only authorization | `docs_limit`, `extension_gate` | Blocked disposition pending explicit revised production-campaign authorization |

## Review Basis And Evidence Status

### Material Reviewed

| Material | Scope Used In This Plan |
| --- | --- |
| `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` | Equations, algorithms, all benchmark and reproduction requirements, performance and AMR claims |
| `src/particles/particles.cpp`, `particles.hpp`, `particles_pushers.cpp`, `particles_moments.cpp`, `particles_tasks.cpp` | Runtime controls, stored particle variables, pushing, deposition, task ordering, partial staged features |
| `src/mhd/mhd_tasks.cpp` | Field coupling and fluid feedback paths |
| `src/mesh/mesh_refinement.cpp`, `src/mesh/load_balance.cpp` | AMR rebuild and load-balancing concerns |
| `src/pgen/tests/pic_parallel_shock.cpp` | Shock initialization, injection, feedback subtraction, refinement behavior |
| `inputs/tests/pic_*.athinput` and `tst/scripts/particles/pic_*.py` | Current tests and proposed publication workflows |
| `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md` | Retired historical record of local validation and large-machine workflow ideas; absorbed and superseded here |
| Entity Toolkit wiki, <https://entity-toolkit.github.io/wiki/> | Comparative reference for documented particle methods, output controls, checkpoints, diagnostics and Frontier practices; not the governing MHD-PIC model |
| `/Users/dbf75/Work/Research/AthenaK/entity` at local commit `a59065fc` | Source-level comparison of particle storage/pushers, current deposition tests, filtering, outputs, timers, parameter/restart contracts and example problems |
| OLCF Frontier User Guide | Frontier node layout, GPU-aware MPI setup, task/GPU binding, and `debug` QOS restrictions |

The checked-in reference-paper directory in the reviewed working tree contains
the manuscript source and figures rather than a rendered PDF. The TeX source
is the operative reference used here. Before publication sign-off, preserve the
exact paper source archive, generated PDF, checksum, and citation metadata in
the qualification artifact set.

### Existing Evidence That May Be Reused Only As Historical Context

The earlier validation note records the following results from the in-progress
local review workspace:

| Historical action | Reported outcome | Proper interpretation |
| --- | --- | --- |
| `python3 run_tests.py particles --cmake=-DCMAKE_BUILD_TYPE=Debug --cmake=-DAthena_ENABLE_MPI=ON` | 25/25 default particle tests passed | Basic local regression evidence only; rerun on a curated baseline |
| Publication proxy archive run with MPI | `overall_status: ok` | Tooling and proxy evidence only; not paper physics validation |
| Metrics and plotting workflow | Completed using a noninteractive Matplotlib cache/backend configuration | Artifact-generation workflow evidence only |

No production sign-off may rely on those runs until they are repeated on clean,
identified commits after the blocking physical-model repairs.

### Evidence Still Required

Every one of the following is currently required:

- A clean production candidate commit series, with no unexplained local changes.
- Analytical unit and convergence verification for the paper-faithful pusher and
  coupling algorithm.
- Exact delta-f and expanding-box implementations with analytical tests.
- AMR ownership, child-state and boundary-buffer refresh stress tests, plus a
  documented AMR deposition contract.
- Quantitative reproduction of every benchmark in the paper.
- Frontier HIP/MPI short-run validation and, after separately approved
  production-campaign authorization, controlled scaling results.
- Restart, decomposition, boundary, robustness, usability, and documentation
  qualification on the final model.
- Frozen, provenance-recorded Entity comparison material for any shared-kernel
  differential test used as production evidence.
- A populated claims registry that limits every release, reproduction and
  comparative statement to its archived evidence.
- Nonlinear saturation qualification for each instability regime presented as
  a production capability, including seed, particle-count, timestep,
  resolution, domain-size and dimensionality sensitivity where applicable.
- Independent cross-code and literature comparisons for any scoped
  state-of-the-art claim.
- A statistical qualification report for every stochastic campaign.
- Resilience, portability, archive-integrity and licensing reports for the final
  release candidate.

## Reference-Paper Reproduction Contract

### Governing Numerical Model

The paper's MHD-PIC model is a required external contract. A paper-reproduction
runtime mode must satisfy all of these properties:

| Paper requirement | Consequence for AthenaK |
| --- | --- |
| Particle state is mass-normalized momentum `p/m`, with `v = (p/m)/gamma`, `gamma = sqrt(1 + (p/m)^2/C^2)`, and energy `(gamma - 1) C^2` | Store and update momentum consistently; implement configurable artificial light speed `C`; diagnose relativistic energy |
| Particle push uses the Boris method with TSC interpolation/deposition and second-order integration | Validate orbit phase, energy error, temporal convergence, interpolation order, and deposition consistency |
| Ideal-MHD induction remains `cE = -u x B`; the CR Hall term is neglected in the tested formulation | Paper mode must not directly add CR current to the CT electric field; feedback must enter through the documented force and energy exchange |
| Stage-two particle momentum and energy changes are subtracted from gas to conserve total momentum and energy | Require cell-resolved conservative exchange and global conservation diagnostics |
| Maximum particle cell crossing is bounded (`Nmax = 2` in the paper) and gyro-angle step obeys `theta_max <= 0.3` | Enforce runtime timestep checks and fail fast with diagnostic output |
| Delta-f evolves particle weights from the background distribution and deposits weighted perturbation moments | A quiet start alone is not delta-f; weight storage, evolution, deposition, restart, and diagnostics are mandatory |
| Expanding/compressing box transforms particle and MHD equations | Implement gas, magnetic-field, and particle terms, with analytical box tests |
| Static/refined meshes and dynamic load balancing are supported | Correctness must persist across SMR/AMR boundaries, migration, repartitioning, and restart |

### Required Paper Experiments

The following matrix is the minimum reproduction set. A row is complete only
when the implemented physics matches the paper, input parameters are recorded,
raw outputs are archived, the analysis is deterministic, and quantitative
acceptance thresholds are documented before the qualifying run.

| Paper experiment | Required result | Current state at review | Blocking work |
| --- | --- | --- | --- |
| Section 5.1 gyro-motion | Relativistic orbit, energy and phase accuracy for the paper's `C`, `v0`, `Omega`, and timestep constraints | Non-relativistic Boris proxy exists | Implement momentum/`C` model and analytical oracle |
| Section 5.2 Bell instability | Measured phase and growth rate versus analytical dispersion in 1D, 2D, and 3D | Bell-like positive-growth proxy exists; coupling differs from paper | Repair coupling and add analytical dispersion comparison |
| Section 5.3 gas plus electron/positron oscillation | Correct oscillation frequency and equivalent behavior on uniform, SMR, and AMR grids | Multi-species parity proxy exists | Repair coupling/AMR, use paper frequency oracle and convergence |
| Section 5.4 non-relativistic parallel shock acceleration | Paper domain, injection, AMR, spectra, morphology, timing and load-balancing behavior | Shock generator and input scaffolding exist but pusher/coupling/`C` are not paper-faithful | Complete model, injection audit, GPU production run and analysis |
| Section 5.5 CR gyro-resonant streaming instability | Polarization-resolved spectra and growth rates against analytical prediction with true delta-f | Growth proxy toggles quiet start | Implement true delta-f and dispersion/spectral oracle |
| Section 5.6 CR pressure anisotropy instability | Prolate/oblate polarization branch selection and quantitative growth against theory | Branch proxy exists | Implement true delta-f, relativistic model and oracle |
| Section 5.7 driven expanding/compressing-box CRPAI | Correct anisotropy evolution, spectra and distribution evolution under box driving | Particle-only slope proxy exists | Implement full expanding box and paper diagnostics |
| Appendix A circularly polarized Alfven wave | Analytic expansion/compression amplitude and phase response | No demonstrated paper-faithful MHD box test | Implement MHD box terms and test |
| Appendix B expanding-box gyro-motion | Analytic gamma/phase history with prescribed expansion rate | No demonstrated paper-faithful result | Implement relativistic expanding-box pusher and test |
| Optimization, scaling and load balance claims | Measured sorting/intermediate-array effects, weak scaling and shock AMR load distribution | Controls appear staged or unproven | Implement/validate production performance mechanisms |

### Nonlinear Saturation Qualification Contract

Linear agreement is mandatory but insufficient for production claims involving
instability evolution. The campaigns below are required after the corresponding
linear oracle passes. They must use predeclared analysis windows, multiple
seeds where sampling noise is present, particle-count and resolution scans, and
an independently regenerated metric table from archived raw output.

| Campaign | Required observables | Minimum sensitivity studies | Required conclusion boundary |
| --- | --- | --- | --- |
| Bell nonlinear saturation | Linear-to-nonlinear transition, magnetic amplification, mode spectrum, filament/cavity morphology, energy partition, saturation amplitude and time history | Particle count, grid resolution, timestep, box size, dimensionality, seed and paper-mode versus separately named Hall-extension mode if the latter is implemented | State exactly which Bell regime and induction model are qualified; do not transfer a non-Hall result to Hall-dominated shock conditions |
| CRSI nonlinear saturation | Polarization-resolved spectra, wave-energy evolution, pitch-angle distribution, particle-distribution evolution, effective scattering or diffusion behavior, isotropization trend and saturation level | Matched reduced nonlinear full-f versus true-delta-f runs, archived delta-f weight distributions and validity envelope; particle count, timestep, resolution, momentum-bin resolution, box size and seed | Claim only the validated drift, distribution and wave regimes inside the demonstrated delta-f envelope; separate linear reproduction from nonlinear transport evidence |
| CRPAI nonlinear saturation | Prolate and oblate branch evolution, wave spectra, pressure anisotropy, effective scattering rate, distribution evolution, quasi-steady or saturated state and driving-rate response | Undriven and driven campaigns; matched reduced nonlinear full-f versus true-delta-f runs; archived delta-f weight distributions and validity envelope; expansion/compression sign; driving rate; particle count; timestep; resolution; box size; seed; and damping model if varied | Register undriven nonlinear saturation, driven saturation and physical-damping transport calibration as separate claims; characterize resolution and effective-damping trends and register an algorithm-specific bounded result if the state remains numerically controlled |

For each campaign, archive failed and outlier seeds rather than selecting only
clean examples. If the expected asymptotic regime cannot be reached under the
currently authorized Frontier policy, record the blocked evidence explicitly
and do not close the associated claim.

For every nonlinear delta-f campaign, archive the full particle-weight
distribution through saturation, predeclare a quantitative validity envelope
for weight excursions and effective represented distribution changes, and run
matched reduced nonlinear full-f cases. If the qualifying trajectory leaves the
envelope or disagrees with the reduced full-f controls outside its predeclared
tolerance, the disposition is `limited` or `rejected`; linear full-f agreement
alone cannot close Q-020, Q-021 or Q-041.

Keep three saturation targets separate:

1. `sun_bai_2023_reproduction`: linear Bell, linear CRSI, linear CRPAI and the
   preliminary driven-box CRPAI result from the bundled manuscript.
2. `athenak_production_mode`: bounded nonlinear Bell and CRSI evolution plus
   separately registered driven and undriven CRPAI saturation behavior in
   explicitly tested regimes.
3. `transport_calibration_extension`: physical-damping and adaptive-delta-f
   extensions capable of measuring effective scattering rates. This is a
   separate future capability, not an implication of reproducing the manuscript.

### Nonlinear Campaign Preregistration Record

Before generating qualifying nonlinear output, freeze a campaign record with:

| Field | Minimum requirement |
| --- | --- |
| Claim and reference | Claim ID, exact reference setup, governing mode, normalization, exclusions and artifact root |
| Seeds | Fixed qualifying seed list and count, or a formally preregistered sequential rule with confidence level, maximum seeds, pilot-reuse policy, inspection schedule, stopping boundaries, multiplicity control and node-hour ceiling |
| Particle count | At least 3 levels spanning the intended production choice |
| Timestep | At least 3 levels or gyro-angle/cell-crossing settings spanning the intended production choice |
| Grid resolution | At least 3 levels for convergence or a justified asymptotic-trend assessment |
| Domain size | At least 2 sizes and more if the dominant mode or saturation statistic remains box-sensitive |
| Dimensionality | 1D/2D/3D wherever the claim crosses dimensionality; otherwise state the restricted dimension |
| Distribution resolution | Momentum-bin and pitch-angle resolution scans for CRSI/CRPAI claims |
| Driving and damping | Driving-rate scan for driven CRPAI; damping-model and rate scan only for a separately named physical-damping extension |
| Analysis freeze | Saturation-window rule, outlier rule, estimators, intervals or bootstrap method, tolerance per observable, script checksum and failure-artifact path |

Plausible-looking saturation does not pass a gate. Acceptance requires the
predeclared quantitative criteria and uncertainty intervals to close for every
observable associated with the claim.

### Independent Cross-Code And Literature Comparison Matrix

Sun and Bai (2023) reproduction is the primary paper contract. It is not by
itself sufficient for a scoped state-of-the-art conclusion. The matrix below is
mandatory comparative work. Future agents must freeze exact reference
citations, code versions when available, input translations, normalization
mappings, extracted data and extraction uncertainties before running AthenaK
comparisons.

| Reference family | AthenaK comparison target | Required observables | Required handling |
| --- | --- | --- | --- |
| Bai et al. (2015), arXiv:1412.1087 | Bell growth and nonlinear evolution; shock precursor behavior; CR-induced Hall-effect applicability | Dispersion, amplification, spectra, morphology, shock diagnostics and parameter mapping | Implement Hall coupling only as a separately named, derived and validated extension; otherwise state the excluded regime explicitly |
| Riquelme and Spitkovsky (2009), arXiv:0810.4565 | Independent kinetic nonlinear Bell/CRCD saturation mechanism | Linear growth and wavelength, turbulence, dominant-wavelength evolution, plasma acceleration, CR deflection and saturation level | Compare only bounded overlaps; use as an independent nonlinear reference rather than as a paper-mode equation oracle |
| Gargaté et al. (2010), arXiv:1002.1701 | Independent hybrid-PIC nonlinear Bell saturation mechanism | Growth, wavelength, magnetic amplification, background-plasma response, CR perpendicular-energy transfer and parameter dependence | Record hybrid-model differences and compare only overlapping observables |
| Zacharegkas et al. (2022), arXiv:2210.08072 | Contemporary hybrid nonlinear Bell saturation survey | Saturation magnetic pressure, CR anisotropic pressure, amplification scaling and parameter envelope | Use for a scoped saturation comparison with exact equation and parameter mapping |
| Mignone et al. (2018), arXiv:1804.01946 | Independent PLUTO conservative-coupling comparison | Linear modes, feedback conservation, representative instability metrics and any matched shock diagnostics | Use matched equations and normalization only; differences become ledger findings, not tolerance adjustments |
| van Marle, Casse and Marcowith (2018), DOI:10.1093/mnras/stx2509 | Independent MPI-AMRVAC AMR-shock and particle-acceleration comparison | Shock morphology, precursor structure, spectra, AMR behavior and documented setup differences | Compare only overlapping regimes; record code-model differences and digitization uncertainty |
| Bai et al. (2019), arXiv:1902.10219 | True delta-f nonlinear CRSI comparison | Growth, polarization spectra, distribution evolution, quasi-linear diffusion, pitch-angle evolution, 90-degree crossing behavior, isotropization, saturation and resolution dependence | Reproduce a bounded published case before claiming nonlinear CR transport capability |
| Plotnikov, Ostriker and Bai (2021), arXiv:2102.11878 | Ion-neutral-damped CRSI extension comparison | Damping-rate dependence, unstable bandwidth, saturation and isotropization | Use only for a separately named physical-damping extension |
| Sun, Bai and Zhao (2024), arXiv:2409.08592 | Adaptive-delta-f, ion-neutral-friction CRPAI saturated-state comparison | Anisotropy, spectra, effective scattering rate, quasi-steady state, driving-rate and friction scaling | Use only for a separately named transport-calibration extension; do not claim calibrated transport from numerical damping |
| Ji and Hopkins (2022), arXiv:2111.14704, plus frozen GIZMO documentation | Reduced-speed-of-light and meshless-code comparison or documented exclusion | Steady-state RSOL invariance, normalization and applicability differences | Record whether a bounded comparison is useful; do not conflate GIZMO RSOL equations with the Sun and Bai artificial-`C` reproduction contract |
| Entity Toolkit frozen snapshot | Shared particle-numerics micro-oracles only | Relativistic pusher cases, compatible trajectory deposition, continuity and restart payload behavior | Never use Entity electromagnetic-field evolution as an oracle for `paper_mhd_pic` gas feedback or induction |

### Pinned External Reference Baseline

The references above are the initial mandatory baseline, not a closed list.
Before a qualifying comparison, record the exact version used, retrieval date,
source URL or DOI, checksum for any archived document or extracted data, and
license or redistribution basis. Re-search the literature when implementation
work reaches Phase 7A and add relevant newer work through the plan-revision
procedure.

The primary reproduction contract is the archived local source for
`arXiv:2304.10568v1`. At this revision,
`docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` has SHA-256
`8f99cec40b8c9fad7f6011c0c32f14724cc0dddabb86db95a7f1687ce04cbdae`.
No rendered PDF is currently checked in. Before qualification, generate the
rendered PDF from the archived source in a recorded environment, store its
SHA-256 in the sign-off manifest, and preserve the exact arXiv-version citation.

Every comparison report must include:

1. exact citation, source URL or DOI, version or commit, and retrieval date;
2. governing-equation overlap and intentional mismatch table;
3. unit and normalization mapping;
4. raw reference data provenance, including digitization method and uncertainty;
5. AthenaK commit, manifest, input checksum and analysis checksum;
6. quantitative residuals and a discrepancy ledger;
7. a scoped conclusion stating what the comparison does and does not qualify.

### Extracted Paper Parameters To Freeze In Reproduction Decks

Agents must recheck these values against the archived manuscript source when
freezing final input decks, then preserve the checked source and input
checksum in the run manifest. The table prevents paper-named tests from
quietly substituting convenient engineering parameters.

| Paper case | Parameters and required setup extracted from the manuscript |
| --- | --- |
| Gyro-motion, Section 5.1 | Use transrelativistic particles with artificial light speed and initial particle speed `C = v0 = 10 u0`; computational normalization includes `u0 = B_g = q/(mc) = Omega = 1`; cubic box length `500 u0/Omega` with 32 cells; impose `theta_max = 0.3` |
| Bell instability, Section 5.2 | Measure real and imaginary parts of the mode frequency against the analytical dispersion relation in 1D, 2D and 3D; the setup uses `Omega = 1e-6 k0 U_A` and `C = 1e3 v_CR` to suppress competing resonant behavior |
| Gas/electron/positron oscillation, Section 5.3 | Compare measured oscillation frequency against `Omega * sqrt(1 + 2 m_e n0/rho)` and require uniform-grid, SMR and AMR results to agree |
| Parallel shock, Section 5.4 | Reflecting left boundary; magnetic field parallel to `x`; inject monoenergetic isotropic supra-thermal particles at the ideal shock surface with `p/m = sqrt(10) u0`, subtracting mass, momentum and energy from gas; two-dimensional domain `(48 x 3.12) * 10^3 c/omega_pi`; `M_A = 30`, gas `gamma = 5/3`, injection fraction `eta = 1e-3`, `C = 1e4 U_A0`; exclude particles injected before `45 Omega_0^-1` where required by the paper analysis; AMR cell sizes from 12 to 3 with refinement/derefinement curvature thresholds 1.0/0.1; compare outputs near `t = 500` and `t = 1200` |
| CR streaming instability, Section 5.5 | Use true delta-f; drift speed `v_d = 2 U_A`; CR mass-density ratio `m n_CR/rho0 = 1e-4`; use the paper's momentum binning and particle count (eight logarithmic bins with 256 particles per cell per bin); measure polarization-resolved spectra and fitted linear growth |
| CR pressure anisotropy instability, Section 5.6 | Use true delta-f; prolate/oblate cases with `xi = 0.99` and `xi = 1.01`; `C = 3e4 U_A = 100 p0/m`, `kappa = 1.75`, and spatial resolution `dx = 20 U_A/Omega`; measure branch polarization and analytical growth |
| Driven CRPAI, Section 5.7 | Use full expanding/compressing-box plus true delta-f; initialize isotropically with `a1 = a^2`, `a2 = a3 = a`; `C = 300 U_A`, `kappa = 1.25`, `dot(a) = 1e-5 Omega_0`; evolve the scaled background distribution and compare spectra plus distribution evolution, not merely a scalar slope sign |
| CPAW box validation, Appendix A | Reproduce the analytic amplitude and phase of the circularly polarized Alfven wave under expansion/compression before accepting any MHD box claim |
| Gyro-motion box validation, Appendix B | Validate analytic gamma and phase histories under the specified expansion case, including the manuscript's `dot(a) = 0.01 Omega_0` test |

### Reproduction Output Contract

For each experiment, archive all of:

- Exact git commit, submodule status, CMake cache, executable checksum, loaded
  modules, reviewed redacted allowlisted environment capture, Slurm script,
  input deck and analysis script. Never archive unrestricted environment dumps.
- Raw AthenaK output, stdout/stderr, scheduler accounting, restart checkpoints,
  and machine-readable diagnostics.
- An immutable metrics JSON or CSV containing analytical reference values,
  measured quantities, errors, resolution, timestep and pass/fail criteria.
- Reconstructed figures with scripts and deterministic plotting environment.
- A short Markdown report stating whether the run reproduces the paper and,
  if not, what failed without altering acceptance criteria after looking at the
  answer.

## Findings From The Current MHD-PIC Implementation

### Blocking Finding Register

| ID | Severity | Finding | Source evidence reviewed | Consequence | Release gate |
| --- | --- | --- | --- | --- | --- |
| PIC-P0-001 | P0 | Coupled particle current is added directly to CT electric fields, but paper mode neglects CR Hall induction; legacy paper-like decks also leave momentum/energy feedback off | `src/mhd/mhd_tasks.cpp` field-source and feedback paths; `inputs/tests/pic_parallel_shock_section54_hpc.athinput` and related PIC decks; paper Section 2 | Current coupled runs solve a different and generally non-conservative model | Implement an explicit paper-faithful coupling mode, forbid incompatible toggles, and pass conservation plus Bell/oscillation oracles |
| PIC-P0-002 | P0 | Pusher advances non-relativistic velocity and energy; configurable artificial light speed is rejected | `src/particles/particles_pushers.cpp`; `src/particles/particles.cpp`; paper Sections 2 and 5 | Cannot reproduce the required orbit, shock, CRSI or CRPAI physics | Implement relativistic momentum Boris with `C`; pass analytic orbit and energy/convergence suites |
| PIC-P0-003 | P0 | `pic_deltaf_mode=on` is quiet-start sampling rather than evolving delta-f weights | `src/particles/particles.cpp`, `particles_moments.cpp`, `src/particles/AGENTS.md`; paper Section 4.4 and Sections 5.5-5.7 | Published delta-f physics and low-noise instability tests are not represented | Implement weight evolution/deposition/restart; reproduce CRSI/CRPAI theory |
| PIC-P0-004 | P0 | Expanding-box code updates particle velocities only, with no corresponding MHD coordinate-source evolution | `src/particles/particles_pushers.cpp`; absence of matching MHD implementation; paper Section 3 and appendices | Driven anisotropy results do not validate the paper model | Implement full box equations; pass CPAW, gyro-motion and driven CRPAI tests |
| PIC-P0-005 | P0 | AMR retains the `MeshBlockPack` object but reconstructs child MeshBlock/coordinate objects while preserving physics modules and particle-owned boundary helpers; the lifetime, neighbor-state and buffer-refresh contract is not proved | `src/mesh/mesh_refinement.cpp`; `src/particles/particles.hpp`; particle and mesh boundary-helper construction paths | Static inspection does not prove a dangling `Particles::pmy_pack` pointer, but stale child state, stale communication metadata or invalid AMR scientific results remain release-blocking risks | Audit ownership and refresh requirements explicitly; repair any demonstrated gap; pass repeated refine/derefine/migration/load-balance/restart tests under memory checking and GPU execution |
| PIC-P0-006 | P0 | Restart output is written directly to its final pathname without a qualified atomic-publication, completion-marker and last-known-good recovery contract; short writes and close failures are not demonstrated to abort safely | `src/outputs/restart.cpp`; `src/outputs/io_wrapper.cpp`; `src/driver/driver.cpp` | A walltime event, I/O fault or interrupted checkpoint may destroy or expose an invalid recovery point and invalidate an expensive campaign | Implement atomic restart publication and fail-closed I/O handling; pass the crash-consistency and I/O failure matrix |
| PIC-P1-001 | P1 | Refinement-boundary deposition policy has not been reconciled with the paper's smoothness-versus-conservation choice | `src/particles/particles_moments.cpp`; paper Section 4.3 | AMR results may be smooth but non-reproducing, or conservative but physically different from the paper | Specify policies, validate both errors, select paper policy for reproduction |
| PIC-P1-002 | P1 | Task-stage ordering and conservative delta exchange are not yet proved against the paper's second-order method | `src/particles/particles_tasks.cpp`; MHD source paths; paper Section 2.3 | Good-looking tests may mask order loss or incorrect exchange timing | Create stage-contract tests and convergence/conservation gates |
| PIC-P1-003 | P1 | Shock setup contains a minimal directed injection construction and cannot yet satisfy the full paper setup | `src/pgen/tests/pic_parallel_shock.cpp`; shock input decks; paper Section 5.4 | Shock acceleration plots cannot yet be publication evidence | Establish exact injection distribution, units and diagnostic contract after P0 work |
| PIC-P1-004 | P1 | Sorting/intermediate-array/load-cost controls are not demonstrated as effective production mechanisms | particle runtime controls; load-balance code; paper Section 4 | Performance and load-balance claims are unsupported | Implement and benchmark, or clearly remove claims and unsupported inputs |
| PIC-P1-005 | P1 | Existing publication-named tests frequently assert sign/parity rather than analytical paper oracles | `tst/scripts/particles/pic_*publication.py` and proxy suites | Passing suite can overstate physical confidence | Reclassify proxies and add quantitative reference tests |
| PIC-P1-006 | P1 | The `pic_entity_deposit_*` tests use Entity-oriented names but currently demonstrate AthenaK integrated-moment and decomposition parity checks rather than a provenance-frozen differential comparison with Entity kernels | `tst/scripts/particles/pic_entity_deposit_mink.py`, `pic_entity_deposit_reflect.py`; Entity `src/kernels/currents_deposit.hpp` and `src/kernels/tests/deposit.cpp` | An apparently external cross-check can be overinterpreted as independent validation | Either rename the tests as internal deposit regressions or add a bounded, frozen Entity-reference comparison for shared deposition behavior |
| PIC-P1-007 | P1 | Orion is a purge-eligible working filesystem rather than durable evidence storage, no OLCF resource is guaranteed as permanent institutional retention, and unrestricted environment capture can preserve sensitive values | OLCF Frontier and storage guidance; historical artifact layout and script templates | Qualification evidence may disappear or archives may expose inappropriate environment data | Mirror operational artifacts to approved nearline storage, export final bundles to an institutional or approved off-site long-term archive with checksums and restore drills, and archive a redacted environment allowlist only |

### What May Be Retained

The current effort should not be discarded wholesale. Subject to review and
repair, the following components are productive foundations:

- The staged test directory and publication artifact workflow.
- Particle diagnostics for current, feedback deltas, work and output fields.
- Restart and decomposition-test scaffolding.
- Shock pgen framework, AMR deck scaffolding and analysis pipeline.
- Serial/MPI parity checks and boundary-focused regression tests.

Retaining these foundations is conditional: their names, documentation and
pass criteria must accurately reflect whether they validate infrastructure,
proxy behavior or actual paper physics.

## Entity Toolkit Comparative Audit And Adoption Decisions

### Review Scope And Boundary

The Entity Toolkit wiki and local Entity source tree were reviewed as an
independent implementation reference. The local comparison used commit
`a59065fc`; that checkout was locally ahead of its configured upstream at the
time of review, so future differential tests must freeze the exact Entity
source snapshot, commit identifier and checksum used to construct an oracle.

Entity is a relativistic Vlasov-Maxwell PIC code with electromagnetic field
evolution. Sun and Bai require a hybrid MHD-PIC model whose paper-tested mode
retains ideal-MHD induction and couples particles to gas through conservative
momentum and energy exchange. Entity is therefore useful for shared particle
numerics, verification patterns and operational discipline. It is not a
physics oracle for the AthenaK `paper_mhd_pic` field-coupling algorithm.

### Transferable Requirements To Adopt

| Entity practice observed in documentation/source | Why it matters here | Required AthenaK action and gate |
| --- | --- | --- |
| Relativistic particle state and analytical pusher checks: Entity stores spatial four-velocity-like momentum state, advances a gamma-dependent Boris/Vay update and tests uniform-field gyro behavior | AthenaK currently lacks the paper-required mass-normalized momentum and configurable artificial `C` implementation | Use Entity only as a test-design reference while implementing the Sun and Bai `p/m`, `gamma(C)` and energy definitions; add kernel-level orbit, phase and energy tests before end-to-end paper tests |
| Previous/current particle positions and trajectory-deposition tests, including discrete current-divergence checks and multiple shape orders | AthenaK already has Entity-named deposit regressions and direct-staggered extension paths, but their independent reference status is incomplete | Create a frozen differential/reference harness for the deposition behavior genuinely shared with Entity, limited to AthenaK-supported orders and compatible current-diagnostic or extension modes; require local stencil, continuity, boundary and MPI checks |
| Per-species particle tracking and payload/provenance fields | Shock injection, delta-f weights, migration and restart require traceable particle histories | Define persistent species ID, injection cohort/source, birth time, tracking ID and delta-f state metadata; verify migration, AMR and restart preserve these fields |
| On-run spectra plus aggregate field/current/particle statistics | Paper shock, CRSI and CRPAI claims depend on distributions and spectra, not only scalar proxies | Provide species- and cohort-resolved weighted spectra with bin edges and definitions in output; cross-check in-run reductions against independent postprocessing |
| Explicit checkpoint metadata and parameter lifecycle distinction between immutable physical setup and mutable run controls | Restarts must not silently change the simulated equations | Version the PIC restart schema; freeze physical mode, species, pusher, `C`, deposition policy, delta-f, box and AMR-coupling settings across restart; reject mismatches with deterministic diagnostics while recording allowed cadence/output changes |
| Named kernel/substep timers and reported particle memory footprint | Frontier performance, particle-aware balancing and capacity estimates require evidence rather than flags | Record pusher, deposition, coupling, communication, sorting, AMR/load-balance, output and checkpoint timers plus memory/count telemetry by species and mesh level; archive them in Frontier manifests |
| Tested bin/tag sorting utilities and periodic particle maintenance controls | Spatial ordering may improve GPU deposition/push locality, but can alter ordering and expose determinism issues | Treat Entity-style sorting as a benchmark candidate only; enable in production only after physics-invariance, restart, CPU/GPU and cost-benefit tests pass |

### Deliberate Non-Adoptions

| Entity capability or practice | Decision for this plan | Reason and required handling |
| --- | --- | --- |
| Electromagnetic Maxwell/current-to-field evolution | Do not import into `paper_mhd_pic` | It changes the governing model; paper mode must use the derived MHD feedback contract. Any future extension needs separate equations, names and tests |
| Digital filtering/smoothing of deposited current | Do not enable for paper reproduction by default | Filtering can change conservation and instability spectra. A separately named extension would require preservation, boundary/corner, growth-rate and convergence studies |
| Entity shock or Bell example decks | Do not treat as reproduction inputs | They exercise a full kinetic electromagnetic model rather than the required hybrid MHD-PIC experiment; retain only transferable diagnostics or bookkeeping ideas |
| Vay/GCA/radiative-cooling/GR machinery | Out of scope for production readiness of this paper mode | They do not close any current Sun and Bai reproduction blocker and would expand the validation surface |
| Entity-specific Frontier recommendation concerning GPU-aware MPI | Do not copy as an AthenaK default | Communication behavior is application dependent. Keep the official OLCF-supported baseline, then perform a short AthenaK A/B correctness/performance comparison and record the chosen environment |

### Consequences For Validation

1. Any Entity-derived oracle must state which quantity is shared between the
   two codes and which physics is intentionally excluded.
2. Existing `pic_entity_deposit_*` cases remain infrastructure regressions
   unless a frozen source/formula provenance record and direct differential
   comparison are added.
3. Entity comparisons cannot close a paper-physics gate involving gas feedback,
   induction, delta-f evolution, expanding-box MHD evolution or paper
   benchmark parameters.
4. Differences discovered by a valid shared-kernel comparison must be entered
   in the findings ledger before a threshold or algorithm is changed.

## Production-Quality Architecture Contract

Before adding further simulations, write a concise design specification under
`docs/` and link it from the particle module documentation. It must define:

### Runtime Modes

| Mode | Intended use | Allowed equations | Forbidden combinations |
| --- | --- | --- | --- |
| `test_particle` | Particle orbit validation in prescribed fluid/field background | No MHD feedback; may support prescribed/ideal background E and B | Any claim of coupled instability or shock backreaction |
| `paper_mhd_pic` | Reproduction of Sun & Bai model | Ideal-MHD induction, conservative CR-to-gas momentum/energy exchange, configurable `C`, optional true delta-f and full box equations | Direct CR-current CT injection; disabled feedback for coupled runs |
| `extended_mhd_pic` | Future physically derived extensions such as Hall-current induction, if desired | Only equations explicitly derived and tested for that extension | Being used under paper-reproduction test names |
| `passive_mhd` / `no_mhd` | Engineering tests and controlled particle experiments | Clearly documented restricted behavior | Publication-level coupled-physics claims |

### Configuration Validation

The input parser must reject configurations whose physics claims cannot be
true. At minimum:

- `paper_mhd_pic` requires momentum and energy exchange to be enabled as a
  coherent algorithm, not independent unchecked toggles.
- `paper_mhd_pic` rejects any direct-current-to-CT induction option.
- Experiments requiring relativistic/artificial-`C` operation reject the
  non-relativistic pusher.
- Delta-f publication cases reject quiet-start-only mode.
- Expanding-box publication cases reject particle-only scaling.
- AMR publication cases reject builds without the selected tested deposition
  policy and an audited post-rebuild ownership, child-state and buffer-refresh
  contract.
- Inputs print the active physical model, units, `C`, timestep limits,
  deposition kernel, AMR policy, feedback algorithm and restart compatibility
  at startup and store that metadata with outputs.

### Conservation And Diagnostics

Production output must permit an independent audit of:

- Particle and gas momentum in each direction and globally.
- Particle kinetic energy, gas total energy and total exchange residual.
- Deposited charge/current or momentum-change moments by level and across
  refinement boundaries.
- Magnetic divergence, field energy and electric work terms.
- Particle timestep constraints, cell-crossing count and gyro-angle maximum.
- Population counts, injection/removal bookkeeping, migration and load cost.
- Delta-f weight statistics and background/perturbation decomposition.
- Expanding-box scale factors and source-term energy/momentum accounting.
- Restart metadata sufficient to assert bitwise or tolerance-level continuity.
- Tracked particle provenance and species/cohort-resolved spectra sufficient to
  audit injection, acceleration, migration and weighted distributions.
- Kernel timers, particle memory and mesh-level population/load information
  sufficient to distinguish physical failure from capacity or imbalance.

## Implementation Workstreams And Acceptance Gates

The phases below are ordered. Agents may develop independent tests in parallel,
but must not start expensive paper or Frontier campaigns before the preceding
blocking gates pass.

## Phase 0: Baseline Curation And Reproducibility Control

### Objective

Turn the current exploratory workspace into a reviewable series of commits and
establish the evidence system used by every later phase.

### Required Work

1. Inventory every modified and untracked PIC-related file in the current
   workspace. Classify changes as:
   - physical-model implementation;
   - correctness fix;
   - diagnostics/restart/output;
   - regression/proxy test;
   - publication tooling;
   - temporary log/artifact that must not be committed.
2. Preserve raw exploratory outputs outside the source tree or under an ignored
   artifact root with a manifest; never mix them into implementation commits.
3. Split retained changes into narrow commits with a stated invariant and a
   matching test. Drop or rewrite misleading publication labels where tests
   are merely proxies.
4. Commit or otherwise provenance-lock the exact paper source archive used for
   reproduction, subject to repository licensing and policy.
5. Add a machine-readable validation manifest schema containing commit,
   executable checksum, test ID, physical mode, parameters, expected oracle,
   measured metrics, pass criteria, resources and artifact directory.
6. Create an external-artifact inventory covering bundled manuscript source and
   figures, rendered PDF, external code snapshots, digitized curves, reference
   datasets, plotting dependencies, citation metadata, checksums,
   redistribution status and allowed archive location.
7. Record whether each external artifact may be committed, archived privately,
   or referenced only through checksum and retrieval instructions.

### Acceptance Gate P0

- Clean candidate branch exists and its retained changes have code review.
- No generated run artifacts are accidentally versioned.
- Every future validation run can write a complete provenance manifest.
- External-artifact licensing and provenance inventory exists.
- This plan is updated with the curated baseline commit hash.

## Phase 1: Paper-Faithful Equation And Interface Specification

### Objective

Make it impossible to confuse staged infrastructure with the paper's physical
model.

### Required Work

1. Write down the discrete state variables and units for particles, gas,
   magnetic field, electric field, current and coupling source terms.
2. Define the exact stage sequence corresponding to the paper's second-order
   scheme: interpolation time, push time, deposited moments, gas feedback,
   boundary exchange, migration and CT updates.
3. Define how feedback is deposited in uniform, SMR and AMR cases and how any
   deviation from paper behavior is exposed as a separate policy.
4. Implement the runtime-mode and validation contract above.
5. Refactor ambiguous flags. A coupled-paper mode should be selected as one
   coherent model, not assembled from booleans that permit unphysical mixes.
6. Ensure help text, input examples, startup logs and module documentation use
   precise names: `proxy`, `engineering`, `paper-reproduction`, or `extension`.

### Tests

- Parser failure tests for every forbidden mode combination.
- Metadata output/restart round-trip tests.
- Stage trace test on a tiny deterministic problem recording which source
  terms are applied at each integrator substep.

### Acceptance Gate P1

- A reviewer can identify the solved equations from a run's input and metadata
  alone.
- Paper test inputs cannot launch under a partial or extension model.

## Phase 2: Relativistic Momentum Pusher And Artificial Light Speed

### Objective

Implement the particle mechanics required by the paper.

### Required Work

1. Replace or supplement velocity-only CR state with mass-normalized momentum
   state suitable for `gamma(C)` and relativistic energy.
2. Implement configurable `pic_cr_light_speed` with unit checks and valid
   range checking. Preserve any explicitly non-relativistic engineering mode
   only under a distinct name.
3. Implement the relativistic Boris update and midpoint interpolation required
   by the paper; define exact field-time centering.
4. Rework all feedback delta channels to use momentum and relativistic kinetic
   energy consistently.
5. Update restart format, output variables, diagnostics and conversion
   utilities with backward-compatibility handling or a deliberate versioned
   restart break.
6. Enforce and record particle timestep constraints for maximum crossing count
   and gyro-angle.
7. Use the frozen Entity pusher implementation and analytical test pattern only
   as a comparative implementation aid: document the mapping from Entity's
   normalized relativistic state to AthenaK's paper-defined `p/m` and
   artificial `C`, and do not substitute Entity's unit conventions.

### Tests

| Test | Purpose | Required pass condition |
| --- | --- | --- |
| Zero-field ballistic motion | State and position update sanity | Exact/tight-roundoff trajectory and restart parity |
| Uniform magnetic-field gyro-motion | Paper Section 5.1 oracle | Second-order phase convergence, bounded energy error, correct `C` dependence |
| Entity-reference pusher microtest | Independent shared-kernel comparison | A frozen reference formulation agrees for matched normalized cases, while AthenaK-specific `C` scaling is validated separately against the paper equations |
| Uniform electric-field energy gain | Momentum-energy consistency | Work equals kinetic-energy change within convergence tolerance |
| Midpoint crossed-field orbit | Field interpolation and centering | Error decreases at designed order |
| Timestep guard tests | Safety | Runs exceeding crossing/gyro limits abort with deterministic reason |
| CPU/MPI/GPU parity | Portability | Metrics agree within predefined roundoff/reduction tolerances |

### Acceptance Gate P2

- Section 5.1 and Appendix B particle-only analytical baselines pass for the
  paper parameters on CPU and Frontier GPU execution.
- No paper-required input is silently forced back to `C=1`.

## Phase 3: Conservative Paper-Faithful MHD Coupling

### Objective

Implement the feedback model actually described in the paper.

### Required Work

1. Disable direct deposited-current modification of CT electric fields in
   `paper_mhd_pic`. If retained as an extension, isolate it under a new mode
   and document its equations and applicability.
2. Deposit the particle momentum and energy changes from the completed push
   into gas source terms with the paper's sign convention and time centering.
3. Specify whether mass exchange occurs in injected/removed particle
   scenarios, and ensure injection subtraction is conservative.
4. Validate boundary handling for deposited exchange before/after boundary
   synchronization, including periodic, reflecting and outflow boundaries.
5. Add global diagnostic residuals that are independent of the update kernel
   implementation.
6. Keep any Entity-derived trajectory/current-deposition comparison explicitly
   separated from paper-mode gas feedback and CT induction qualification.

### Tests

| Test | Scope | Required evidence |
| --- | --- | --- |
| Single particle plus uniform gas exchange | Local sign/time centering | Exact expected momentum and energy transfer |
| Many-particle periodic exchange | Global conservation | Residual convergence and roundoff floor |
| Cell-boundary and physical-boundary crossing | Shape/boundary deposition | Conservation and symmetry properties |
| Frozen Entity-reference trajectory deposition | Shared deposition mechanics in diagnostic/extension paths only | AthenaK-supported shape orders match the recorded local-stencil/continuity reference without being used to validate paper-mode induction |
| Serial/MPI/decomposition invariance | Communication correctness | Same metrics for multiple decompositions |
| Bell linear mode | Coupled physical oracle | Growth rate and phase versus analytical dispersion |
| Gas/electron/positron oscillation | Backreaction and grid support | Frequency versus analytic value on uniform/SMR/AMR |

### Acceptance Gate P3

- No paper input modifies induction through an unapproved current source.
- Total momentum and energy exchange meets declared conservation tolerances.
- Paper Sections 5.2 and 5.3 pass their analytical metrics before proceeding.

## Phase 4: True Delta-F Implementation

### Objective

Implement the delta-f method described in the reference paper rather than
using quiet-start sampling as a surrogate.

### Required Work

1. Define supported background distribution functions `f0`, their parameters
   and their transformation under ordinary and expanding-box dynamics.
2. Store the initial invariant/background information and evolving particle
   weight required for `w = 1 - f0(current)/f(initial)` or the precisely
   documented equivalent discretization.
3. Deposit perturbation moments using delta-f weights while retaining
   diagnostics for full-f versus perturbation quantities.
4. Define conservation implications and error diagnostics, because the paper
   notes that delta-f no longer provides the same machine-level conservation
   property as full-f exchange.
5. Restart all delta-f state exactly and make incompatible restart files fail
   clearly.
6. Keep quiet-start sampling as a separately named variance-reduction option,
   available independently of delta-f.

### Tests

- Weight remains zero or analytically constant for an equilibrium orbit.
- Controlled perturbation deposition matches a direct numerical reference.
- Full-f and delta-f agree on signal while measured noise reduction is
  reproducible for an appropriate linear test.
- Restart and MPI decomposition preserve weight statistics and measured growth.
- CRSI and CRPAI tests compare polarization-resolved growth rates and spectra
  directly with the paper's analytical prediction and input distributions.

### Acceptance Gate P4

- Sections 5.5 and 5.6 are quantitatively reproduced with true delta-f.
- Test names and documentation no longer describe quiet-start-only runs as
  delta-f physics.

## Phase 5: Full Expanding And Compressing Box

### Objective

Implement the paper's coordinate transformation consistently for particles,
gas and magnetic fields.

### Required Work

1. Represent time-dependent box scale factors and rates unambiguously and
   store them in run metadata and restart output.
2. Add the MHD density, momentum, energy and magnetic-field evolution/source
   terms required by the paper's equations.
3. Implement particle momentum transformation and its centered coupling to the
   Boris advance.
4. Implement the delta-f background transformation used by driven CRPAI.
5. Audit CFL/timestep behavior and conservation accounting in changing-volume
   coordinates.

### Tests

| Test | Paper target | Required result |
| --- | --- | --- |
| Expanding/compressing circularly polarized Alfven wave | Appendix A | Amplitude and phase follow analytic solution with convergence |
| Expanding-box gyro-motion | Appendix B | Gamma and phase follow analytic history |
| Static-rate-zero recovery | Baseline consistency | Identical result to ordinary paper mode |
| Driven CRPAI | Section 5.7 | Reproduce anisotropy development, polarization and spectra/distribution evolution |
| Restart across changing scale factor | Reliability | Continuous metrics and exact restored coordinate state |

### Acceptance Gate P5

- Appendices A and B pass and the Section 5.7 experiment is reproducible.
- Particle-only box scaling is no longer presented as full expanding-box MHD-PIC.

## Phase 6: AMR, SMR, Boundary, Migration, Restart And Load Balance

### Objective

Make mesh evolution and distributed execution scientifically trustworthy.

### Required Work

1. Audit the AMR lifetime and refresh contract precisely. `MeshBlockPack`
   itself is retained, while child MeshBlock and coordinate objects are
   reconstructed and physics modules plus particle-owned boundary helpers
   persist. Inventory every retained pointer, view, neighbor table, MPI request,
   buffer and cached dimension. Add explicit post-AMR refresh or reconstruction
   methods wherever required, and document why any retained object is safe.
2. Test refine, derefine, load-balance, particle migration and restart in
   combination, including repeated transitions.
3. Decide and document deposition behavior at refinement interfaces:
   - `paper_smooth`: reproduce the paper's locally smooth but not individually
     conservative feedback behavior;
   - optional `conservative`: an AthenaK extension, separately documented and
     validated if retained.
4. Implement meaningful particle-aware load costs and document how they are
   computed, reduced and used by the balancer.
5. Verify boundary-condition support for every production mode and reject
   unsupported combinations before execution.

### Tests

| Test family | Required configurations | Pass criteria |
| --- | --- | --- |
| Forced AMR lifetime stress | Repeated refine/derefine with particles and coupling | No invalid memory, no missing particles, invariant diagnostics |
| Refinement deposition characterization | Interface traversal under each policy | Documented smoothness and conservation errors; paper policy matches intended result |
| Uniform/SMR/AMR oscillation | Section 5.3 setup | Frequency and conservation within oracle tolerance |
| Boundary exchange | periodic, reflecting, outflow as supported | Conservation/sign/symmetry and MPI parity |
| Restart | no-MHD engineering, paper coupled uniform, AMR coupled, delta-f, box modes | Continuation agrees with uninterrupted run |
| Load balance | Particle-dominated imbalanced case and shock | Cost correlates with measured time; no physical regression |
| GPU memory/lifetime | Frontier HIP build and AMR stress | No runtime faults or unexplained numerical corruption |

### Acceptance Gate P6

- AMR coupled results are permitted only after ownership/lifetime,
  neighbor-state, communication-buffer refresh and deposition-policy tests pass.
- Section 5.3 AMR behavior and load-balancing requirements are supported by
  archived metrics, not merely by successful completion.

## Reliability, Crash-Consistency And I/O Failure Contract

Restart equivalence is insufficient if a checkpoint interrupted during writing
can replace the last usable recovery point. Treat crash consistency and
fail-closed artifact handling as production requirements.

### Required Implementation Work

1. Publish restart checkpoints atomically: write a temporary or `.partial`
   artifact, complete all rank writes, close successfully, record completion
   metadata and checksum, then promote the artifact to its restartable name.
2. Preserve at least one prior completed checkpoint until the new checkpoint is
   verified and promoted.
3. For per-rank restart mode, publish a completion manifest only after every
   member exists with the expected checksum and size.
4. Make short writes, close failures, unwritable paths, incomplete per-rank
   sets, malformed headers, payload corruption and checksum mismatches fail
   deterministically with actionable diagnostics.
5. Make analysis pipelines fail closed when required raw outputs, metrics,
   checksums or figure inputs are absent or incomplete.
6. Record failure state, last-known-good restart, exit code and artifact
   completeness in the run manifest.

### Required Failure Matrix

| Test | Required result |
| --- | --- |
| Soft Athena `-t` stop | Final checkpoint is complete and restartable; continued metrics match an uninterrupted control within tolerance |
| Slurm pre-timeout pilot | Tested wrapper leaves a valid restart before scheduler termination |
| Forced termination during compute | Previous completed restart remains usable |
| Forced termination during checkpoint write | Incomplete checkpoint is not selectable; previous completed checkpoint remains usable |
| Truncated shared MPI restart | Deterministic rejection with precise diagnostic |
| Missing or truncated per-rank restart member | Deterministic rejection before evolution begins |
| Corrupt header, particle payload or checksum | Deterministic rejection |
| Unwritable output directory and bounded simulated short write | Nonzero exit, preserved prior checkpoint and manifest failure record |
| Partial analysis artifact set | Qualification pipeline fails closed instead of producing a passing summary |

### Timeout Safety Policy

Athena's internal `-t` limit stops after a completed cycle and final output still
requires time. Before adopting a production walltime margin, measure worst-case
cycle duration and worst-case checkpoint duration for the selected profile and
problem size. Require:

```text
scheduler_walltime - athena_internal_walltime
  > measured_worst_cycle_time + measured_worst_checkpoint_time + safety_margin
```

OLCF supports `#SBATCH --signal=B:USR1@300`, but a signal is useful only when a
tested wrapper traps it and requests a safe stop or checkpoint. Qualify any such
wrapper with short debug-QOS pilots before relying on it.

## Phase 7: Parallel Shock Reproduction

### Objective

Reproduce the paper's non-relativistic parallel shock acceleration benchmark
with a correct underlying model.

### Required Work

1. Audit the shock problem generator against Section 5.4 line by line:
   reflecting-wall geometry, upstream flow, parallel magnetic field, Mach
   number, particle injection timing, spatial distribution, isotropic
   monoenergetic velocity distribution, mass/momentum/energy removal from gas,
   AMR thresholds, output cadence and load-balancing scheme.
2. Replace any minimal engineering injection distribution with the exact
   paper-reproduction distribution; retain the engineering version only under
   a distinct input/problem identifier.
3. Use the paper-required artificial light speed and all dimensional
   conversions; prohibit launch if the pusher is not in `paper_mhd_pic`.
4. Establish coarse, fine and AMR inputs with frozen provenance and an analysis
   script measuring the quantities shown in the paper, including spectra and
   spatial morphology at prescribed times.
5. Validate checkpoint/restart equivalence in short, compliant debug pilots.
   Do not chain a full production shock campaign through the `debug` QOS; if a
   full run exceeds permitted debug use, mark it blocked pending an explicitly
   revised execution authorization.
6. Store particle provenance and species/injection-cohort spectra needed to
   audit injected populations and acceleration histories, with independent
   offline reconstruction of the same binned distributions.

### Required Shock Qualification Matrix

| Run | Purpose | Required comparison |
| --- | --- | --- |
| Short uniform coarse | Early shock setup sanity | Shock structure, injection and conservation diagnostic |
| Short uniform fine | Resolution response | Compare with coarse and convergence/trend expectation |
| Short AMR | Refine/derefine and load-balance sanity | Agreement with fine representation in resolved region |
| Full coarse/fine/AMR paper runs | Section 5.4 reproduction | Morphology, spectra and reported-time figures/metrics |
| Short restart AMR pilot | Restart correctness under a compliant debug job | Same diagnostics as uninterrupted short control within tolerance; not authorization to chain a full run |

### Acceptance Gate P7

- Section 5.4 plots and quantitative metrics are reproducible from archived raw
  outputs and analysis scripts.
- AMR speed/memory benefit and physical equivalence are both stated honestly.

## Phase 7A: Nonlinear Saturation And Independent Comparison

### Objective

Establish the evidence needed for any nonlinear-instability, transport,
comparative or scoped state-of-the-art statement.

### Required Work

1. Complete the nonlinear Bell, CRSI and separately registered driven and
   undriven CRPAI campaigns defined above using preregistered statistical
   designs and the final qualified runtime modes.
2. Close the independent cross-code and literature comparison matrix for every
   capability named in a comparative statement.
3. Record discrepancies as findings. Do not hide disagreement by changing
   windows, dropping seeds or broadening tolerances after inspection.
4. Separate paper reproduction from AthenaK extension results in manifests,
   figures and prose.
5. Identify regimes that remain blocked by execution authorization, model
   exclusions or insufficient evidence and mark the corresponding claims
   `limited` or `unsupported`.
6. Close the Ji-Hopkins/GIZMO RSOL decision before interpreting comparative
   evidence: either preregister a bounded comparison with exact metrics and
   equation mapping, or archive a reviewer-approved exclusion explaining why
   the distinct RSOL formulation is not informative for the claimed profile.

### Acceptance Gate P7A

- Every claimed nonlinear capability has archived saturation, sensitivity and
  uncertainty evidence.
- Every scoped state-of-the-art statement names the compared references,
  equations, regime, metrics and limitations.
- No unresolved cross-code discrepancy is hidden behind a passing proxy or
  paper-reproduction result.

## Phase 8: Performance, GPU Portability, Usability And Documentation

### Objective

Convert a physically correct model into a maintainable production module.

### Required Work

1. Determine whether the paper's intermediate-array and sorting strategies are
   implemented, useful and appropriate for Kokkos/HIP. Remove inert runtime
   flags or implement them with benchmark evidence.
2. Implement particle-aware load-balancing metrics and performance logging.
3. Benchmark kernel time, communication, sorting, deposition, AMR overhead,
   memory use and output overhead on representative Frontier sizes.
4. Add required telemetry for per-stage timers, per-species particle counts
   and memory, mesh-level particle distribution and load imbalance, including
   these fields in run manifests and qualification metrics.
5. Confirm deterministic/reproducible reduction expectations; document where
   GPU/MPI order permits tolerance-level rather than bitwise agreement.
6. Retain the OLCF-supported GPU-aware MPI environment as the initial Frontier
   baseline, then run a short correctness/performance A/B comparison with any
   Entity-motivated communication alternative before altering production
   settings.
7. Provide user documentation:
   - governing equations and supported modes;
   - parameter reference with valid/invalid combinations;
   - paper-reproduction quick start;
   - Frontier runbook and accounting requirements;
   - output/analysis definitions;
   - troubleshooting and failure interpretation.
8. Make failure messages specific: unsupported mode, timestep constraint,
   missing GPU-aware MPI support, AMR policy mismatch and restart incompatibility.
9. Execute the reliability, resilience and portability matrix below on the
   release candidate. Resolve findings before sign-off rather than treating
   successful physics runs as a substitute.
10. Complete the archive-integrity and licensing gate, including an artifact
    restore drill in a fresh directory.
11. Populate the claims registry with the exact supported production envelope,
    comparative conclusions and known exclusions.

### Acceptance Gate P8

- GPU runs are physically equivalent to CPU references within documented
  tolerances.
- Performance claims include raw job accounting and reproducible scripts.
- A new user can launch a supported small validation run without reading source.
- Restart, interruption, corrupt-checkpoint and output-failure behavior is
  deterministic and documented.
- Archive restore and licensing audits are complete.

## Phase 9: Final Qualification And Production Sign-Off

### Objective

Provide a defensible yes/no production decision.

### Required Sign-Off Bundle

- Clean release candidate tag or commit and complete change review.
- Host-side unit/regression results and build logs.
- Frontier HIP/MPI build and runtime environment record.
- Complete analytical validation metrics.
- All paper benchmark artifacts and reconstructed figures.
- AMR/restart/decomposition/boundary/load-balance robustness reports.
- Performance and resource-usage report.
- Known limitations and unsupported-mode documentation.
- Populated claims registry with reviewer dispositions and scoped wording.
- Statistical qualification report for every stochastic production claim.
- Nonlinear Bell, CRSI and CRPAI saturation reports for every corresponding
  production capability claimed by the release.
- Independent cross-code and literature comparison report for every scoped
  state-of-the-art statement.
- Reliability, resilience and portability matrix report.
- Archive-integrity, restore-drill and licensing report.
- Cumulative Frontier node-hour ledger, demonstrating budget compliance.
- Final checklist signed in this document with links to immutable artifact paths.

### Release-Profile Closure Matrix

Do not use one profile's closure to imply another. Record the selected profile
in the release manifest and close every mandatory gate listed for it.

| Release or evidence profile | Mandatory gates |
| --- | --- |
| `sun_bai_2023_reproduction` evidence bundle | Q-001 through Q-009, Q-011 through Q-013, Q-016, Q-018, Q-023, Q-025 through Q-027, Q-034, Q-036 through Q-038 |
| Production-ready `paper_mhd_pic` | All `sun_bai_2023_reproduction` gates plus Q-010, Q-017, Q-019 through Q-021, Q-024, Q-028, Q-030, Q-031, Q-040 and Q-041 |
| Scoped state-of-the-art statement | Production-ready `paper_mhd_pic` plus Q-022, Q-028 through Q-035, Q-039 and Q-040; unsupported extensions must close through explicit exclusion rather than omission |
| Optional `extended_mhd_pic` capability | Production-ready `paper_mhd_pic` plus an `implemented_and_qualified` extension-specific outcome: Q-029 for Hall Bell, Q-032 for ion-neutral-damped CRSI, or Q-033 for CRPAI transport calibration. An `excluded_as_unsupported` outcome closes scope bookkeeping but cannot qualify the extension |

`Q-014` is the terminal aggregate sign-off report, not a child gate in any
profile. Its lifecycle is deliberately two-step to avoid circular sign-off:

1. After every selected-profile child gate closes, freeze a review-ready
   aggregate claim manifest with disposition `pending_terminal_review`.
2. Q-014 reviews that immutable candidate, records the terminal disposition,
   signs and archives the final manifest, and updates the aggregate claim row.

A paper-reproduction bundle may be archived while a long run remains
execution-authorization blocked, but it must be labeled `blocked`, not
`qualified`.

### Definition Of Production-Ready

The implementation is production-ready only when:

1. The code solves a documented, physically coherent MHD-PIC model for each
   enabled runtime mode.
2. `paper_mhd_pic` quantitatively reproduces the required results of Sun & Bai
   (2023), or any discrepancy is demonstrated to be a documented AthenaK
   correction with independently validated consequences.
3. Accuracy, convergence, conservation, restart, AMR, MPI, GPU and
   performance evidence is complete and repeatable.
4. Unsupported options fail explicitly and supported options are usable from
   documentation alone.
5. All critical and high-severity findings are closed with evidence; no waiver
   may hide a failed physical invariant.
6. Every comparative or state-of-the-art statement is scoped through the claims
   registry and supported by the required independent evidence.
7. The archived qualification bundle can regenerate critical metrics and paper
   figures in a fresh environment.

## Test And Verification Strategy

### Test Taxonomy

| Test class | Runtime scale | Purpose | May support production sign-off? |
| --- | --- | --- | --- |
| Unit/oracle | seconds to minutes | Equations, state conversions, deposition, parser, restart serialization | Yes, as component evidence |
| Regression | minutes | Prevent known failures and preserve supported behavior | Yes, as reliability evidence |
| Proxy | minutes to hours | Detect broad behavior trends or pipeline breakage | No, unless promoted with an analytical contract |
| Physics validation | minutes to hours | Compare measured physical values with analytical reference | Yes |
| Paper reproduction | hours and multi-run campaigns | Reproduce documented paper benchmark results | Yes, mandatory |
| Performance/scaling | Frontier jobs | Demonstrate production feasibility and balance | Yes, mandatory for production operation |

The human-readable taxonomy above maps to the claims registry as follows:
`Unit/oracle` and `Regression` map to `unit/regression`; `Proxy` maps to
`engineering_proxy`; `Physics validation` maps to `physics_validation`; `Paper
reproduction` maps only to `sun_bai_2023_reproduction`; and
`Performance/scaling` is supporting evidence for `athenak_production_mode`.
Never create a second spelling for a machine-readable claim class.

### Required Reclassification Of Existing Publication-Named Cases

Several current cases have useful regression value but names or notes that can
be overread as physical qualification. Phase 0 must rename them or change their
manifest classification before any release artifact bundle is produced.

| Current case family | Required classification now | Reason | Promotion requirement |
| --- | --- | --- | --- |
| `pic_em_vacuum_wave*` | `engineering_proxy` or renamed MHD linear-wave regression | Current deck is an MHD linear-wave adaptation with inactive neutral drift particles, not a full electromagnetic vacuum-wave PIC benchmark | Add an explicit governing-equation oracle and rename the test to match what is actually solved |
| `pic_langmuir_frequency_proxy` and manifest publication alias | Narrow `physics_validation` after renaming | No-MHD uniform-`B` Boris frequency anchor checks a component gyrofrequency oracle, but it is not a self-consistent electrostatic Langmuir-wave validation and not the paper's relativistic gyro reproduction | Rename as nonrelativistic uniform-`B` gyrofrequency validation; keep its claim narrowly scoped |
| `pic_two_stream_growth*` | `engineering_proxy` | Current AthenaK-adapted control is not a full kinetic two-stream reproduction | Implement and validate an explicitly supported physical model before promotion |
| `pic_weibel_growth*` | `engineering_proxy` | Current transverse-current trend control is not a full electromagnetic Weibel reproduction | Implement and validate an explicitly supported physical model before promotion |
| `pic_bell_growth*` | `engineering_proxy` until rewritten | Positive growth and rank parity do not reproduce the analytical Bell dispersion relation, and the current coupling differs from paper mode | Require paper-faithful coupling, real and imaginary dispersion, convergence, dimensionality and nonlinear saturation gates |
| `pic_multispecies_backreaction*` | `engineering_proxy` until frequency oracle closes | Oscillation and parity scaffolding are useful, but paper-grid physics must be measured quantitatively | Require analytical frequency, conservation and uniform/SMR/AMR agreement |
| `pic_crsi_deltaf*` | `engineering_proxy` until true delta-f exists | Current `pic_deltaf_mode=on` is quiet-start sampling only | Implement evolving weights, theory comparison, spectra and nonlinear saturation |
| `pic_crpai_polarization*` | `engineering_proxy` until true delta-f exists | Branch-sign separation is not quantitative CRPAI reproduction | Implement evolving weights, theory comparison, spectra and nonlinear saturation |
| `pic_expanding_box_anisotropy*` | `engineering_proxy` until full box equations exist | Particle-only velocity scaling is not expanding-box MHD-PIC | Implement MHD plus particle box equations and appendix oracles |
| `pic_entity_deposit_*` | `unit/regression` unless frozen differential oracle is added | Entity-oriented names do not alone demonstrate an independent Entity comparison | Freeze source/formula provenance and run bounded shared-kernel differential checks |
| `pic_amr_shock_lb*` and `run_pic_shock_scan.py` Orszag-Tang scans | `engineering_proxy` and stress testing | Clean execution, output generation and `problem/ot_mach` scans do not reproduce the paper parallel shock | Use distinct paper Section 5.4 decks and quantitative shock analysis |

Every manifest, figure bundle, README and report must preserve this
classification. Proxy figures must use a visible `engineering_proxy` label and must never
appear in a paper-reproduction table. The checked-in exploratory manifest is not
a paper-reproduction manifest: no currently checked-in case may be cited as
`sun_bai_2023_reproduction` evidence until it is promoted through the claims
registry with the required physical oracle. Checked-in tooling must use
evidence-class names such as `ArtifactCase`, `entity_core_engineering`,
`extended_benchmark_engineering`, and `extended_engineering`. Retain historical
identifiers containing `publication` only where archive or CLI compatibility
requires them; annotate them explicitly as legacy and unqualified. Plotting and
metric tools must not silently fall back from a reproduction artifact to a
proxy artifact while retaining a reproduction label.

Checked-in enforcement lives in
`tst/publication/test_pic_artifact_taxonomy.py`. Keep it passing as a Phase-0
gate whenever artifact tooling, legacy aliases, deck labels, checksums,
watermarks, scan defaults or HPC emit-only rules change.

### Host And Small-Scale Validation Ladder

Run each stage on a clean release candidate and archive logs:

1. Compile with warnings enabled and run fast serial unit/oracle tests.
2. Run Debug plus MPI regression:

   ```bash
   cd tst
   python3 run_tests.py particles \
     --cmake=-DCMAKE_BUILD_TYPE=Debug \
     --cmake=-DAthena_ENABLE_MPI=ON
   ```

3. Run Release plus MPI regression and quantitative physics tests.
4. Run sanitizer or debug-memory checking builds for host-supported paths,
   especially forced AMR lifetime and restart tests.
5. Run the publication artifact pipeline only after each physical model gate
   relevant to its cases has passed.

The historical 25/25 result from the exploratory workspace is a useful starting
reference, not an excuse to skip the clean-baseline rerun.

### Mandatory Analytical Oracles

| Component | Oracle |
| --- | --- |
| Relativistic particle pusher | Analytic uniform-B orbit, phase and energy |
| Shared trajectory/deposition extensions | Frozen Entity-reference stencil and discrete-continuity comparison, restricted to compatible AthenaK modes |
| Feedback coupling | Global momentum/energy exchange and paper oscillation frequency |
| Bell mode | Analytical dispersion relation, growth and phase |
| Delta-f | Controlled equilibrium/perturbation weight result plus CRSI/CRPAI dispersion |
| Expanding box | CPAW and gyro analytical solutions |
| AMR | Uniform/SMR/AMR physics agreement plus documented interface error |
| Shock | Paper diagnostic and figure reconstruction, with conservation histories |

### Tolerance Policy

Every new test must state:

- measured quantity and units;
- analytical/reference value;
- expected discretization-order trend if applicable;
- absolute and relative tolerance and why it is appropriate;
- allowed CPU/MPI/GPU reduction variation;
- failure artifact location.

Never relax a threshold merely because a run failed. A threshold change requires
a documented physical or numerical justification, repeated resolution evidence,
and review in the finding/change ledger.

### Statistical Qualification Protocol

Stochastic simulations require a declared statistical design before qualifying
data are generated. A single clean seed is exploratory evidence only.

For every stochastic physics-validation, paper-reproduction, nonlinear or
cross-code campaign:

1. Register the target observables, estimator, units, fit window, exclusion
   rules, seed list, parameter grid, reference values and tolerance rationale
   before inspecting qualifying output.
2. Use a fixed-sample design by default. Freeze the qualifying seed list,
   sample count, endpoint hierarchy and multiplicity-control policy across
   observables and sensitivity scans before inspecting qualifying outputs.
   Exploratory pilot seeds may inform the design only if they are excluded from
   the qualifying estimate.
   If compute cost requires a sequential design, preregister the confidence
   level, maximum seed count, pilot-reuse policy, inspection schedule,
   alpha-spending or equivalent error-control rule, stopping boundaries,
   multiplicity treatment across observables, and campaign node-hour ceiling.
   Stop only at a predeclared boundary; never add optional seeds informally
   because an interval is inconvenient.
3. Perform particle-count, timestep and resolution sensitivity studies. Add
   box-size, dimensionality, momentum-bin, driving-rate and damping-model scans
   where they can affect the claimed regime.
4. Report central estimates and uncertainty intervals. Use bootstrap or another
   justified method for noisy fitted quantities and preserve the resampling
   configuration.
5. Archive every attempted seed, including failures and outliers. Exclusions
   require a predeclared rule or a ledger finding with physical justification.
6. Freeze reference-data extraction, digitization scripts and uncertainty when
   comparing published figures without machine-readable source data.
7. Recompute critical metrics with an independent analysis implementation or
   reviewer-owned script before closing a claim.
8. Distinguish bitwise reproducibility, deterministic metric reproducibility,
   and statistical agreement. State which level each CPU/MPI/GPU comparison
   requires.
9. Do not tune thresholds, select windows or discard seeds after looking at the
   qualifying answer. Any necessary revision invalidates the candidate dataset
   and requires a registered rerun.

### Statistical Qualification Report Schema

Store one machine-readable report and one short Markdown interpretation for
each stochastic claim. Include:

| Field group | Required content |
| --- | --- |
| Identity | Claim ID, campaign ID, commit, executable checksum, physical mode, input and analysis checksums |
| Preregistration | Frozen observables, estimators, windows, exclusions, seed list, scan grid, reference values and tolerances |
| Sampling rule | Fixed qualifying sample count plus endpoint hierarchy and multiplicity control, or preregistered sequential confidence level, maximum seeds, pilot-reuse policy, inspection schedule, stopping boundaries, multiplicity control and node-hour ceiling |
| Samples | Every attempted seed and parameter point, completion state, artifact path and exclusion disposition |
| Estimates | Per-run metrics, aggregate central estimates, uncertainty intervals, bootstrap or interval method and configuration |
| Sensitivity | Particle-count, timestep, resolution, box-size, dimensionality and problem-specific scan results |
| Independent check | Reviewer-owned or independently implemented critical-metric recomputation |
| Disposition | Pass, fail, limited or rerun-required with linked findings and claim-register update |

### Reliability, Resilience And Portability Matrix

Production readiness requires controlled failure behavior, not only successful
physics runs.

| Gate family | Required configurations | Required evidence |
| --- | --- | --- |
| Build portability | Debug and Release; warnings-enabled host build; supported CPU serial/OpenMP if enabled; CPU MPI; Frontier HIP/MPI | Build logs, compiler and dependency versions, warning disposition and identical active-mode metadata |
| Host memory safety | Sanitizer or supported debug-memory tooling on small particle, coupling, AMR and restart cases | No invalid access, use-after-free, leak affecting campaign operation or undefined-behavior finding |
| GPU execution safety | Frontier HIP smoke, pusher/coupling oracles and repeated AMR stress | No device faults, hangs, unexplained corruption or CPU/GPU metric discrepancy outside declared tolerance |
| Decomposition robustness | Multiple MeshBlock layouts and MPI rank counts, including rank-boundary migration | Invariant physical metrics within declared deterministic or reduction-order policy |
| Forced walltime recovery | Deliberately interrupted short runs with checkpoint cadence exercised | Resume instructions work; continued metrics agree with uninterrupted control |
| Restart integrity | Uniform, AMR, delta-f, box and shock-pilot restarts; permitted runtime-control change; forbidden physical-setting change | Schema version recorded; immutable settings rejected on mismatch; continuation equivalence demonstrated |
| Output interruption | Interrupted or incomplete output write where safely testable | Failure is diagnosed; previous valid restart remains usable; incomplete artifact is not silently accepted |
| Corrupt or incomplete checkpoint | Truncated or checksum-mismatched checkpoint fixture | Deterministic rejection with actionable diagnostic; no undefined continuation |
| Capacity and storage failure | Bounded low-space or write-failure fixture where safe; manifest preflight | Clean abort, preserved previous checkpoint and explicit artifact status |
| Analysis reproducibility | Fresh environment and independent critical-metric recomputation | Metrics and figures regenerate from archived raw data with recorded tool versions |

### Mandatory Portability Execution Matrix

Archive compiler identity, module resolution, Kokkos configuration, CMake cache,
binary checksum, selected environment profile, run manifest and metrics for
every row.

Before running this matrix, publish a supported-toolchain declaration for the
release candidate. It must list the supported host compiler, MPI, Kokkos and
backend combinations; state whether OpenMP is supported; and list intentionally
excluded backends as `unsupported` rather than silently omitting them. If
OpenMP is enabled for the release, the two OpenMP rows below are mandatory. If
it is excluded, replace them with a reviewed unsupported disposition.

| Build/runtime | Required scope |
| --- | --- |
| Host Debug serial | Unit and analytical-oracle suite |
| Host Release serial | Unit/oracle suite plus representative quantitative physics validations |
| Host Debug OpenMP, if supported | Unit/oracle suite, thread-count variation and deterministic/reduction-order policy |
| Host Release OpenMP, if supported | Representative quantitative physics validations and thread-count parity |
| Host Debug MPI `np=2,4` | Decomposition, boundary, migration, restart and failure-path checks |
| Host Release MPI `np=2,4` | Quantitative physics parity and representative stochastic metrics |
| Host ASan/UBSan or supported memory tooling | Particle push/deposition, AMR, migration, restart and output-failure fixtures |
| Frontier HIP/MPI one rank | Basic device execution, analytical oracles and memory behavior |
| Frontier HIP/MPI eight ranks on one node | Rank/GPU mapping, decomposition and physical parity |
| Frontier HIP/MPI multi-node debug pilot | Communication, restart, AMR and load-balance qualification |
| Shared-file and per-rank restart modes | Complete restart and crash-consistency matrix |

Gate ownership is hierarchical: `Q-036` closes the restart-publication and I/O
failure implementation contract; `Q-024` aggregates broader resilience
evidence, including `Q-036`; and `Q-012` is the release-level reliability
report that closes only after `Q-024` plus MPI, decomposition, boundary and
restart-continuation evidence close. Likewise, `Q-038` closes the Frontier
environment-profile decision and is a child of the broader portability gate
`Q-025`.

### Archive Integrity And Licensing Gate

Before the first third-party ingestion, archive transfer or Frontier submission,
classify licensing, redistribution, access-control and sensitive-data handling.
Repeat the review before release sign-off:

1. Inventory every committed or archived third-party paper source, figure,
   digitized dataset, script, kernel excerpt and reference output.
2. Record citation, source URL or DOI, retrieval date, checksum, license or
   redistribution basis, and whether the artifact may be committed, archived
   privately, or referenced only by retrieval instructions.
3. Preserve the exact Sun and Bai manuscript source archive and rendered PDF
   used for qualification with checksums, subject to redistribution policy.
4. Preserve frozen Entity comparison provenance without importing incompatible
   code or licensing assumptions into AthenaK.
5. Version the validation-manifest schema and include checksums for executable,
   input deck, analysis code, raw reference data and compact result table.
6. Run an archive restore drill: reconstruct at least one representative
   analytical result, one stochastic metric table and one paper figure from the
   recorded artifact bundle in a fresh directory.
7. Add a release checklist entry signed by a reviewer confirming that archived
   evidence is retrievable and legally retainable.
8. Scrub paths, job names, scripts, logs and allowlisted environment values for
   sensitive or controlled strings before submission and before export.

### Publication Artifact Workflow

The previous large-machine note established a useful pattern: archive raw
simulation output, compute metrics, and generate figures noninteractively. Keep
that structure but require physical qualification first. On systems without a
display, set a deterministic Matplotlib environment, for example:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/athenak_pic_mplconfig"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/athenak_pic_xdg"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
```

Artifacts from proxy campaigns must be labeled `proxy` and must not appear in
a paper-reproduction result table as though they validate the paper.

## Living Findings And Plan-Revision Procedure

This document must change as evidence changes. Agents must not defer plan
updates until the end of a large campaign.

### Finding Ledger Format

For every newly discovered issue, add an entry to the current finding register
or its successor with:

| Field | Required content |
| --- | --- |
| ID and severity | Stable `PIC-P0/P1/P2-###` identifier |
| Discovery date and commit | Exact code state examined |
| Invariant at risk | Physical, numerical, reliability, performance or usability claim |
| Evidence | Source locations, reproducer, metrics and artifact paths |
| Consequence | Which tests, paper claims or production uses are blocked |
| Proposed remedy | Smallest defensible fix or design decision |
| Verification gate | Exact test/metric required for closure |
| Status | Open, implementing, verifying, closed or rejected with reason |

### Revision Rules

1. Add a finding immediately when a claimed invariant is not demonstrated.
2. Before implementing a design-changing fix, update the affected phase and
   acceptance test in this plan.
3. When a failure exposes an inadequate test, improve the test before closing
   the implementation finding.
4. When a test is promoted from proxy to physics validation, document its
   analytical oracle and tolerance before using its result.
5. When a finding is closed, link the fixing commit and archived validation
   artifacts; never delete the historical issue.
6. After each Frontier job, update node-hour accounting and the phase evidence
   table before submitting another job.
7. At each release-candidate point, reread the paper reproduction matrix and
   confirm that no requirement was lost during implementation.

### Recommended Evidence Layout In The Repository

Keep source-controlled instructions and compact metrics under `tst/publication/`.
Keep large machine output outside git in the Frontier execution root:

```text
/lustre/orion/ast207/proj-shared/dfielding/PIC/
  source/                       # clean checkout(s), identified by commit
  build/<commit>/<config>/      # build trees and CMake caches
  bin/<commit>/<config>/        # immutable executable copies/checksums
  inputs/<campaign>/            # exact submitted input decks
  jobs/<campaign>/              # exact sbatch scripts and submission metadata
  logs/build/                   # compiler logs
  logs/slurm/                   # scheduler stdout/stderr
  runs/<campaign>/<run_id>/     # simulation raw output and restarts
  metrics/<campaign>/           # compact machine-readable analyses
  figures/<campaign>/           # generated plots
  manifests/<campaign>/         # provenance manifests
  ledger/node_hours.jsonl       # authoritative append-only accounting ledger
  ledger/node_hours.csv         # derived RFC-4180 human-readable index
  ledger/mirror_receipts.jsonl  # immutable non-recursive mirror receipts
```

### Archive Retention And Restore Procedure

Orion is the working filesystem, not durable evidence storage. OLCF documents
that Orion work areas are not backed up and are purge-eligible. Before closing a
qualification gate:

1. Build a sign-off bundle containing raw outputs required for reproduction,
   compact metrics, figures, inputs, analysis scripts, logs, scheduler
   accounting, executable checksum, CMake cache, selected environment profile,
   ledger snapshot and claim-register entry.
2. Generate a SHA-256 manifest for the bundle and verify it before transfer.
3. Export required evidence to Kronos as a nearline operational mirror where
   approved, and to an institutional or otherwise approved off-site long-term
   archive for final sign-off. No OLCF filesystem is the sole permanent record.
   Record archive destination, transfer timestamp, checksum verification,
   retention owner, migration procedure and any access restriction.
4. Keep compact source-controlled metadata pointing to the immutable archive
   location. Do not rely on an Orion path as the sole long-term evidence record.
5. Restore a representative bundle into a fresh directory and regenerate at
   least one analytical metric, one stochastic summary and one paper figure
   before final sign-off.
6. Use OLCF data-transfer guidance and data-transfer nodes for substantial
   transfers. Recheck current retention policy before every archive operation.

For the node-hour ledger, durability is continuous rather than a final-gate
action. The authoritative ledger is append-only JSONL; `node_hours.csv` is a
derived human-readable index. Serialize each primary event as canonical JSON
with UTF-8 encoding, lexicographically sorted keys, no insignificant whitespace
and a trailing newline. Compute `event_sha256` over that serialization while
omitting only the `event_sha256` field itself. After every reservation, job-ID
attachment, cancellation and accounting reconciliation, append a hash-chained
primary event and mirror the event plus sequence head to a preflighted approved
durable target. Use a project-approved DTN transfer or approved remote `rsync`,
`scp` or Globus workflow; record the exact transport and target before the first
reservation. Each primary event contains `sequence_number`,
`previous_event_sha256`, `event_sha256`, `event_type`, identity, accounting and
artifact fields. After transport succeeds, append a separate immutable
`mirror_ack` receipt to `mirror_receipts.jsonl`. A receipt refers to the
mirrored primary event's SHA-256 and contains `mirror_destination`,
`mirror_transport`, `mirror_acknowledged_utc` and `mirror_ack_sha256`.
Serialize receipts with the same canonical-JSON rule and compute
`mirror_ack_sha256` while omitting only that field. Receipts do not participate
in the primary hash chain and are not themselves mirrored, avoiding recursive
acknowledgement. The durable target's primary-chain head and the Orion
primary-chain head are the heads compared during preflight; the local receipt
stream proves completed transport. Never mutate a primary event to add
acknowledgement state. CSV output uses RFC 4180 quoting and the fixed column
order defined below, projecting the latest matching receipt fields when present.
Before creating a reservation, require the validator to confirm that the Orion
ledger and durable sequence head agree and that a test acknowledgement can
still reach the mirror target. If recovery is ambiguous or the mirror is
unavailable, block all new submissions until a manual `sacct` reconciliation
and reviewed ledger repair are complete.

## Frontier Operating Procedure

### Official Frontier Constraints And Configuration

Consult the current OLCF Frontier User Guide before every substantial campaign:

<https://docs.olcf.ornl.gov/systems/frontier_user_guide.html>

At the time this plan was prepared, the guide states that:

- A Frontier node provides eight GPU-visible AMD MI250X GCDs.
- The common one-GPU-per-rank placement uses eight MPI ranks per node, seven
  CPU cores per rank and `--gpus-per-task=1 --gpu-bind=closest`.
- GPU-aware Cray MPICH requires `craype-accel-amd-gfx90a`, `rocm` and
  `MPICH_GPU_SUPPORT_ENABLED=1`; Cray compiler wrappers with HIP sources may
  require explicit ROCm include/link flags.
- The `debug` QOS permits only one user job in any state and rejects walltimes
  exceeding 2 hours. It is intended for short debugging, not production job
  chaining.
- Orion project-work storage is not backed up and is purge-eligible; required
  evidence must be mirrored to approved nearline storage and exported to an
  institutional or approved off-site long-term archive.
- The default Frontier GPU mode is `HSA_XNACK=0`; setting `HSA_XNACK=1` changes
  page-migration behavior and must be treated as an experimental profile until
  benchmark evidence supports promotion.

This project deliberately restricts all currently authorized Frontier jobs to
`debug`. Because OLCF identifies that QOS as short, non-production use and
prohibits production job chaining, do not attempt to evade the policy by
manually chaining production restart segments. Use `debug` only for compliant
short qualification/debug runs. If complete paper reproduction or a production
scaling experiment cannot be completed within the permitted use of `debug`,
stop, record the blocked work and obtain explicit revised authorization before
using any other execution policy.

### Authorized Debug-QOS Qualification Scope

The only currently authorized Frontier work is bounded short-run qualification
under `#SBATCH -q debug`. These jobs may:

- verify HIP/MPI compilation, GPU placement and startup metadata;
- run analytical pusher, coupling, delta-f and expanding-box oracles;
- run short AMR, migration, load-balance, boundary and restart stress tests;
- run short Bell, oscillation, CRSI and CRPAI qualification windows;
- run bounded shock setup, output, checkpoint and single-resolution pilots;
- measure small controlled performance A/B comparisons needed to choose a
  recorded environment setting.

These jobs may not:

- chain restart segments to approximate a production allocation;
- run a full shock reproduction that requires production-scale time or repeated
  restart chaining;
- run a scaling study whose purpose is production performance characterization;
- exceed the one-job, two-hour or cumulative 1000 node-hour constraints;
- use another QOS, partition or account policy without explicit revised
  authorization.

### Future Approved Production Campaign Boundary

Full shock reproduction, nonlinear saturation runs that exceed compliant debug
use, and controlled production scaling are a separate future campaign. They are
blocked, not waived. Before any such submission, future agents must:

1. Obtain explicit user and site-policy authorization for the execution policy.
2. Amend this plan with the approved QOS, partition, walltime, node-count
   ceiling, concurrency rule and revised budget envelope.
3. Preserve the same `/lustre/orion/ast207/proj-shared/dfielding/PIC` execution
   root, provenance manifests and node-hour ledger unless explicitly revised.
4. Freeze the qualifying commit, binary checksum, deck matrix, statistical
   design, predeclared metrics and stop conditions before submission.
5. Run the authorized production jobs only after all prerequisite small
   correctness gates are closed.
6. Record the resulting accounting and artifacts separately from debug-QOS
   evidence so no short pilot is mistaken for a production result.

### Mandatory Resource And Accounting Policy

1. Set `PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC` for every
   Frontier action. Do not run simulations in another project directory.
2. Use only `#SBATCH -q debug` and a time request no larger than `02:00:00`.
3. Submit at most one job at a time. Direct `sbatch` use is prohibited. Before
   each validated-wrapper submission, run:

   ```bash
   squeue -u "$USER" -h -o "%i %q %T %j"
   ```

   If any `debug` job is listed in any state, do not submit another.
4. Calculate the maximum additional budget before submission:

   ```text
   maximum_node_hours = requested_nodes * requested_walltime_hours
   ```

   Refuse submission if:

   ```text
   cumulative_consumed_node_hours
     + currently_reserved_node_hours
     + maximum_node_hours > 1000
   ```

5. Reserve maximum requested node-hours in the locked ledger before submission.
   After the job ends, obtain actual accounting with `sacct`, reconcile the
   reservation and append the terminal state. The consumed value is:

   ```text
   billed_nodes = max(requested_nodes, scheduler_reported_allocated_nodes)
   consumed_node_hours = billed_nodes * elapsed_seconds / 3600
   ```

6. Keep the authoritative ledger under `$PIC_ROOT/ledger/node_hours.jsonl`;
   regenerate `$PIC_ROOT/ledger/node_hours.csv` as a derived RFC-4180 index
   after every acknowledged append. Also copy its current state into each
   campaign manifest and durably mirror every
   sequence-numbered mutation as required by the archive procedure. The ledger
   columns must include:

   ```text
   sequence_number,previous_event_sha256,event_sha256,event_type,timestamp,reservation_id,
   submission_id,job_id,git_commit,campaign,test_id,requested_nodes,
   scheduler_reported_allocated_nodes,billed_nodes,
   requested_walltime,reserved_node_hours,elapsed_seconds,consumed_node_hours,
   cumulative_consumed_node_hours,state,reconciled,artifact_dir,
   mirror_destination,mirror_transport,mirror_acknowledged_utc,mirror_ack_sha256,
   notes
   ```

7. Use the budget deliberately. Begin with one-node correctness tests and only
   grow when a previous result satisfies its gate. A failed setup should cost
   no more than one small debug job before it is corrected.

### Frontier Environment Profiles

Store the reusable environment setup as
`$PIC_ROOT/jobs/frontier_pic_environment.sh` and source it from build/run
scripts. Verify module versions when the OLCF software stack changes. Start from
the minimum supported profile. Treat page migration and communication tuning as
experiments until matched AthenaK A/B evidence supports promotion.

```bash
#!/bin/bash
set -euo pipefail

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1

PIC_FRONTIER_PROFILE="${PIC_FRONTIER_PROFILE:-frontier_minimum_supported}"
export PIC_FRONTIER_PROFILE

unset HSA_XNACK
unset MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED
unset MPICH_OFI_NIC_POLICY
unset MPICH_GPU_IPC_CACHE_MAX_SIZE
unset MPICH_MPIIO_HINTS
unset MPICH_OFI_NUM_CQ_ENTRIES
unset FI_MR_CACHE_MONITOR
unset FI_CXI_RX_MATCH_MODE

case "$PIC_FRONTIER_PROFILE" in
  frontier_minimum_supported)
    ;;
  frontier_xnack1_experimental)
    export HSA_XNACK=1
    export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
    ;;
  frontier_ofi_tuned_experimental)
    export MPICH_OFI_NIC_POLICY=GPU
    export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
    export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
    export MPICH_OFI_NUM_CQ_ENTRIES=131072
    export FI_MR_CACHE_MONITOR=kdreg2
    export FI_CXI_RX_MATCH_MODE=software
    ;;
  *)
    printf "Unsupported PIC_FRONTIER_PROFILE=%s\n" "$PIC_FRONTIER_PROFILE" >&2
    exit 1
    ;;
esac

record_pic_environment() {
  local name value
  for name in PIC_FRONTIER_PROFILE MPICH_ENV_DISPLAY MPICH_VERSION_DISPLAY \
      MPICH_GPU_SUPPORT_ENABLED MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED \
      MPICH_OFI_NIC_POLICY MPICH_GPU_IPC_CACHE_MAX_SIZE MPICH_MPIIO_HINTS \
      MPICH_OFI_NUM_CQ_ENTRIES FI_MR_CACHE_MONITOR FI_CXI_RX_MATCH_MODE \
      HSA_XNACK ROCM_PATH; do
    if value="$(printenv "$name")"; then
      printf "%s=%s\n" "$name" "$value"
    else
      printf "%s=<unset>\n" "$name"
    fi
  done
}
```

| Profile | Purpose | Rule |
| --- | --- | --- |
| `frontier_minimum_supported` | Initial qualifying baseline | Site-supported modules, GPU-aware MPICH and explicit rank/GPU binding only |
| `frontier_xnack1_experimental` | Managed-memory experiment | Compare correctness, runtime and memory against the default `HSA_XNACK=0` behavior before promotion |
| `frontier_ofi_tuned_experimental` | Communication experiment | Add OFI/CXI/MPI-IO settings as one controlled profile; benchmark against minimum supported baseline |
| `frontier_selected_production` | Future approved production profile | Define only after matched evidence and explicit future-campaign authorization |

Entity documentation reports application-specific Frontier communication
experience. Treat that information as a prompt for a bounded AthenaK benchmark,
not as permission to change the minimum environment above. After correctness
gates pass, run matched short A/B cases with the OLCF-supported GPU-aware MPI
baseline and any proposed alternative, checking physical metrics, runtime and
memory before selecting a recorded production setting. Archive the redacted
allowlist emitted by `record_pic_environment`; do not archive unrestricted
`env | sort` output because environment variables can contain sensitive values.

### Improved Frontier Build Script Template

The example supplied with the original validation request builds an unrelated
problem generator and leaves configuration provenance incomplete. MHD-PIC
qualification must build the built-in PIC problem generators, store the cache
and log, and make the resulting executable immutable by commit/config identity.

Create `$PIC_ROOT/jobs/build_frontier_pic.sh`:

```bash
#!/bin/bash
set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
SRC_DIR="${PIC_ROOT}/source/athenak-DF"
ENV_FILE="${PIC_ROOT}/jobs/frontier_pic_environment.sh"

source "$ENV_FILE"

cd "$SRC_DIR"
GIT_COMMIT="$(git rev-parse --short=12 HEAD)"
CONFIG=hip-mpi-release-paper-pic
BUILD_DIR="${PIC_ROOT}/build/${GIT_COMMIT}/${CONFIG}"
BIN_DIR="${PIC_ROOT}/bin/${GIT_COMMIT}/${CONFIG}"
LOG_DIR="${PIC_ROOT}/logs/build"

mkdir -p "$BUILD_DIR" "$BIN_DIR" "$LOG_DIR"

cmake -S "$SRC_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
  -DPROBLEM=built_in_pgens \
  2>&1 | tee "${LOG_DIR}/${GIT_COMMIT}.${CONFIG}.configure.log"

cmake --build "$BUILD_DIR" --parallel 32 \
  2>&1 | tee "${LOG_DIR}/${GIT_COMMIT}.${CONFIG}.build.log"

cp -p "${BUILD_DIR}/src/athena" "${BIN_DIR}/athena" 2>/dev/null || \
  cp -p "${BUILD_DIR}/athena" "${BIN_DIR}/athena"

sha256sum "${BIN_DIR}/athena" \
  > "${BIN_DIR}/athena.sha256"
cp -p "${BUILD_DIR}/CMakeCache.txt" "${BIN_DIR}/CMakeCache.txt"
git -C "$SRC_DIR" status --short --branch \
  > "${BIN_DIR}/git_status.txt"
module -t list 2> "${BIN_DIR}/modules.txt"
record_pic_environment > "${BIN_DIR}/environment.allowlist.txt"

printf "Built %s for %s at %s\n" "$GIT_COMMIT" "$CONFIG" "${BIN_DIR}/athena"
```

Notes:

- Use `Kokkos_ARCH_AMD_GFX90A`, the explicit gfx90a architecture spelling in
  current Kokkos, rather than relying on legacy naming.
- The ROCm include and HIP link flags are retained because the Frontier guide
  describes them for Cray compiler-wrapper GPU-aware builds. Record any future
  toolchain simplification in this plan and in the manifest.
- Run compiles on an appropriate Frontier login/build workflow according to
  OLCF policy; run simulations only through allocated compute resources.
- Do not qualify a build if `git_status.txt` contains unexplained modified or
  untracked source/config files.

### Improved Debug-QOS Job Template

Create one immutable job script per test/campaign under `$PIC_ROOT/jobs/`.
This template uses one GPU per MPI rank as recommended by the Frontier guide:

Before submission, implement and run
`$PIC_ROOT/jobs/create_pre_submit_manifest.py`. It must create an immutable
`pre_submit_manifest.json` under
`$PIC_ROOT/manifests/<campaign>/<submission_id>/` before the validator reserves
node-hours. Generate a new UUID `submission_id` for every submission, including
reruns. Snapshot the executable, input deck, job script, environment profile,
analysis scripts, timeout artifact and configuration into that immutable
submission directory. Freeze the commit, checksums, resource request, control-
plane schema version and every snapshot path, including the compute-node
verifier, wrapper, validator and reconciler. Never point an allocation at a
mutable shared job script, input, verifier or executable path.
The validator and snapshotted compute-node script must recompute and verify checksum
values, not merely test that fields exist. Runtime state belongs in a separate
append-only status stream; never mutate the pre-submit manifest from inside an
allocation. Direct `sbatch` use is prohibited: submit only through the validated
wrapper below.

Freeze the submission wrapper, validator, reconciler, ledger initializer and
schema under an immutable checksummed control-plane version directory. Invoke
that frozen wrapper path, record every control-plane checksum and schema version
in each pre-submit manifest and ledger event, and reject checksum drift before
reservation. Do not execute mutable scripts directly from `$PIC_ROOT/jobs/`.

The timeout-margin evidence artifact is mandatory and expires after any code,
toolchain, environment-profile, mesh, PPC, rank-layout, restart-mode,
checkpoint-cadence or output-mode change. It records profile checksum, module
stack, mesh size, PPC, node/rank count, restart mode, checkpoint cadence, output
mode, sample count, maximum measured cycle and checkpoint times, selected
safety margin, measurement timestamp and expiry rule. Reject stale artifacts
and unresolved `REPLACE_*` placeholders.

```bash
#!/bin/bash
#SBATCH -J PIC_GYRO_P2_1N
#SBATCH -A AST207
#SBATCH -o /lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log
#SBATCH -t 00:20:00
#SBATCH -p batch
#SBATCH -q debug
#SBATCH -N 1

set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
ENV_FILE="${PIC_ROOT}/jobs/frontier_pic_environment.sh"
GIT_COMMIT=REPLACE_WITH_COMMIT
CONFIG=hip-mpi-release-paper-pic
CAMPAIGN=paper_pusher_validation
TEST_ID=gyro_section51_p2_1n
SUBMISSION_ID="${PIC_SUBMISSION_ID:?validated wrapper must export PIC_SUBMISSION_ID}"
RESERVATION_ID="${PIC_RESERVATION_ID:?validated wrapper must export PIC_RESERVATION_ID}"
SLURM_JOB_KEY="${SLURM_JOB_ID:-${SLURM_JOBID:?Slurm job ID is required}}"
RUN_ID="${TEST_ID}.${SUBMISSION_ID}.${SLURM_JOB_KEY}"
NNODES="${SLURM_JOB_NUM_NODES}"
NRANKS=$((8 * NNODES))
NTHREADS=7
ATHENA_WALLTIME=REPLACE_FROM_TIMEOUT_MARGIN_MANIFEST

SNAPSHOT_DIR="${PIC_ROOT}/manifests/${CAMPAIGN}/${SUBMISSION_ID}/snapshot"
PRE_SUBMIT_MANIFEST="${PIC_ROOT}/manifests/${CAMPAIGN}/${SUBMISSION_ID}/pre_submit_manifest.json"
ENV_FILE="${SNAPSHOT_DIR}/frontier_pic_environment.sh"
ATHENA="${SNAPSHOT_DIR}/athena"
INPUT="${SNAPSHOT_DIR}/${TEST_ID}.athinput"
SNAPSHOT_VERIFIER="${SNAPSHOT_DIR}/verify_compute_node_snapshot.py"
OUTDIR="${PIC_ROOT}/runs/${CAMPAIGN}/${RUN_ID}"
RUNTIME_DIR="${PIC_ROOT}/runtime_status/${CAMPAIGN}/${RUN_ID}"
RUNTIME_STATUS="${RUNTIME_DIR}/events.jsonl"

record_exit() {
  status="$?"
  set +e
  mkdir -p "$RUNTIME_DIR"
  printf '{"event":"job_exit","exit_status":%s,"finished_utc":"%s"}\n' \
    "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$RUNTIME_STATUS"
}
trap record_exit EXIT

preflight_failed() {
  reason="$1"
  set +e
  mkdir -p "$RUNTIME_DIR"
  printf '{"event":"preflight_failed","reason":"%s","finished_utc":"%s"}\n' \
    "$reason" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$RUNTIME_STATUS"
  exit 1
}

mkdir -p "$OUTDIR" "$RUNTIME_DIR" "${PIC_ROOT}/logs/slurm" || \
  preflight_failed runtime_directory_creation_failed
printf '{"event":"job_start","submission_id":"%s","reservation_id":"%s","started_utc":"%s"}\n' \
  "$SUBMISSION_ID" "$RESERVATION_ID" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  >> "$RUNTIME_STATUS"
test -r "$PRE_SUBMIT_MANIFEST" || preflight_failed manifest_unreadable
test -r "$SNAPSHOT_VERIFIER" || preflight_failed snapshot_verifier_unreadable
python3 "$SNAPSHOT_VERIFIER" \
  --manifest "$PRE_SUBMIT_MANIFEST" \
  --submission-id "$SUBMISSION_ID" \
  --reservation-id "$RESERVATION_ID" || preflight_failed snapshot_checksum_mismatch
source "$ENV_FILE" || preflight_failed environment_source_failed
export SLURM_EXPORT_ENV=ALL
export OMP_NUM_THREADS="$NTHREADS"
module -t list 2> "$RUNTIME_DIR/modules.txt" || preflight_failed module_capture_failed
record_pic_environment > "$RUNTIME_DIR/environment.allowlist.txt" || \
  preflight_failed environment_capture_failed

run_status=0
srun -N "$NNODES" -n "$NRANKS" -c "$NTHREADS" \
  --gpus-per-task=1 --gpu-bind=closest \
  "$ATHENA" \
  -i "$INPUT" \
  -d "$OUTDIR" \
  -t "$ATHENA_WALLTIME" \
  job/basename="$RUN_ID" || run_status="$?"
exit "$run_status"
```

Adjust job name, walltime, node count, campaign and input only after checking
that the worst-case budget remains within the approved cap.
`ATHENA_WALLTIME` is a required substitution: the generated job script must
insert an internal limit justified by the selected
profile/problem-size timeout-margin artifact in `pre_submit_manifest.json`.
Use shorter walltimes for quick correctness tests. Reconcile `sacct` after the
batch job reaches a terminal state; an in-job `sacct` command is not sufficient.

### Submission And Ledger Procedure

Because only one debug job can exist at once, submission is an explicit serial
workflow. Implement and test
`$PIC_ROOT/jobs/validate_and_reserve_debug_job.py` and
`$PIC_ROOT/jobs/reconcile_debug_job.py` before the first Frontier submission.
The validator must parse the immutable job script and pre-submit manifest, then reject
submission unless all policy checks pass:

- job script requests exactly `#SBATCH -q debug`;
- requested walltime is `<=02:00:00`;
- requested nodes and worst-case node-hours are parsed successfully;
- consumed plus currently reserved plus requested node-hours is `<=1000`;
- `squeue` shows no `debug` job for the user in any state;
- no previous reservation remains unreconciled;
- executable, input, analysis and job-script checksums exist and recompute to
  the frozen values;
- the timeout-margin artifact matches the selected profile and problem size;
- a restricted-export smoke fixture proves that rank 0 observes the approved
  allowlisted profile after the job sources its snapshot and sets
  `SLURM_EXPORT_ENV=ALL`;
- the Orion ledger sequence head matches the durably mirrored sequence head;
- the durable mirror target passes an acknowledgement preflight;
- the submission UUID is new and every runtime dependency was snapshotted under
  its immutable manifest directory;
- paths, job names, scripts, logs and allowlisted environment values pass the
  sensitive-string scrub;
- executable, input, job, output, log and manifest paths remain under
  `/lustre/orion/ast207/proj-shared/dfielding/PIC`;
- the ledger can be locked and a reservation row can be written before
  submission.

The submission wrapper must fail closed:

```bash
#!/bin/bash
set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
JOB_SCRIPT="$1"
RUN_MANIFEST="$2"
LEDGER_JSONL="${PIC_ROOT}/ledger/node_hours.jsonl"
LEDGER_CSV="${PIC_ROOT}/ledger/node_hours.csv"
CONTROL_PLANE_VERSION=REPLACE_WITH_CHECKSUMMED_CONTROL_PLANE_VERSION
CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}"
VALIDATOR="${CONTROL_PLANE_DIR}/validate_and_reserve_debug_job.py"
RECONCILER="${CONTROL_PLANE_DIR}/reconcile_debug_job.py"
PENDING_SUBMISSION_FILE="${PIC_ROOT}/ledger/pending_submission.json"

reservation_id=""
job_id=""
attached=0

cleanup_unattached_job() {
  status="$?"
  set +e
  if [[ -n "$job_id" && "$attached" -eq 0 ]]; then
    scancel "$job_id"
    python3 "$RECONCILER" \
      --ledger-jsonl "$LEDGER_JSONL" \
      --ledger-csv "$LEDGER_CSV" \
      --job-id "$job_id" \
      --run-manifest "$RUN_MANIFEST" \
      --state submission_attach_failed
  elif [[ -n "$reservation_id" && -z "$job_id" ]]; then
    python3 "$VALIDATOR" cancel-reservation \
      --ledger-jsonl "$LEDGER_JSONL" \
      --ledger-csv "$LEDGER_CSV" \
      --reservation-id "$reservation_id"
  fi
  exit "$status"
}
trap cleanup_unattached_job ERR INT TERM

mkdir -p "${PIC_ROOT}/ledger"
test -r "$LEDGER_JSONL" || {
  printf "Missing initialized mirrored ledger: %s\n" "$LEDGER_JSONL" >&2
  printf "Run the reviewed initialize_debug_ledger.py genesis procedure first.\n" >&2
  exit 1
}

reservation_id="$(
  python3 "$VALIDATOR" reserve \
    --job-script "$JOB_SCRIPT" \
    --run-manifest "$RUN_MANIFEST" \
    --ledger-jsonl "$LEDGER_JSONL" \
    --ledger-csv "$LEDGER_CSV" \
    --required-qos debug \
    --max-walltime 02:00:00 \
    --max-node-hours 1000
)"

submission_id="$(python3 "$VALIDATOR" submission-id \
  --run-manifest "$RUN_MANIFEST" \
  --reservation-id "$reservation_id")"
snapshot_job_script="$(python3 "$VALIDATOR" snapshot-path \
  --run-manifest "$RUN_MANIFEST" \
  --reservation-id "$reservation_id" \
  --kind job-script)"
pending_tmp="${PENDING_SUBMISSION_FILE}.tmp.$$"
printf '{"reservation_id":"%s","submission_id":"%s","state":"reserved_not_attached"}\n' \
  "$reservation_id" "$submission_id" > "$pending_tmp"
mv "$pending_tmp" "$PENDING_SUBMISSION_FILE"

if ! job_id="$(sbatch --parsable \
    --comment "pic-reservation=${reservation_id}" \
    --export "PIC_RESERVATION_ID=${reservation_id},PIC_SUBMISSION_ID=${submission_id}" \
    "$snapshot_job_script")"; then
  python3 "$VALIDATOR" cancel-reservation \
    --ledger-jsonl "$LEDGER_JSONL" \
    --ledger-csv "$LEDGER_CSV" \
    --reservation-id "$reservation_id"
  exit 1
fi
printf '{"reservation_id":"%s","submission_id":"%s","job_id":"%s","state":"submitted_not_attached"}\n' \
  "$reservation_id" "$submission_id" "$job_id" > "$pending_tmp"
mv "$pending_tmp" "$PENDING_SUBMISSION_FILE"
python3 "$VALIDATOR" attach-job-id \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --reservation-id "$reservation_id" \
  --job-id "$job_id"
attached=1
rm -f "$PENDING_SUBMISSION_FILE"
trap - ERR INT TERM
printf "Submitted %s with reservation %s\n" "$job_id" "$reservation_id"

printf "After completion run:\n"
printf "python3 %q --ledger-jsonl %q --ledger-csv %q --job-id %q --run-manifest %q\n" \
  "$RECONCILER" "$LEDGER_JSONL" "$LEDGER_CSV" "$job_id" "$RUN_MANIFEST"
```

Before the first submission, initialize an explicit genesis primary event,
mirror it, append its receipt and derive the CSV index through the reviewed
`initialize_debug_ledger.py` workflow. Never create an empty ledger implicitly.

After completion, run the reconciler to archive `sacct`, replace reserved usage
with actual usage, and record job state. A human or agent must still inspect the
run identity, artifact completeness, failure state and cumulative total before
the next submission. The validator must refuse a new reservation while any prior
row remains unreconciled, while the durable ledger head is stale, or while
`pending_submission.json` exists without a reviewed attachment, cancellation and
reconciliation. Because shell traps cannot close the post-`sbatch` process-loss
window, recovery must query `squeue` and `sacct` for the reservation token stored
in the Slurm comment, attach or cancel the matching job, reconcile usage and
clear the pending marker only after the durable mirror acknowledges the
transition.
Every reservation, attachment, cancellation and reconciliation mutation must be
durably mirrored before the wrapper or reconciler reports success. If
`scancel`, `sacct` reconciliation or durable mirroring cannot be confirmed,
block later submissions and require manual recovery.

### Authorized Debug-QOS Qualification Budget Envelope

The budget below is deliberately conservative. It uses maximum requested
node-hours, not expected elapsed cost. Do not spend a later tier while its
prerequisite phase gate is open.

| Tier | Prerequisite | Suggested job pattern | Maximum requested node-hours | Cumulative ceiling after tier |
| --- | --- | --- | ---: | ---: |
| F0 Build/run smoke and GPU mapping | Clean baseline plus parser tests | Up to 6 x 1 node x 0.25 h | 1.5 | 1.5 |
| F1 Relativistic pusher and conservative coupling oracles | P2/P3 implementation complete | Up to 16 x 1 node x 0.5 h | 8 | 9.5 |
| F2 Delta-f and box analytical tests | P4/P5 implementation complete | Up to 16 x 1 node x 1 h | 16 | 25.5 |
| F3 AMR/restart/decomposition/load balance validation | P6 implementation complete | Up to 20 jobs varying 1-4 nodes x <=1 h, budget-capped | 50 | 75.5 |
| F4 Bell, oscillation, CRSI, CRPAI short qualification cases | P3-P6 analytical gates closed | Up to 24 compliant debug jobs varying 1-8 nodes x <=2 h, budget-capped | 200 | 275.5 |
| F5 Shock debug pilots and single-job resolution checks | P7 setup review complete | Up to 12 compliant debug jobs varying 1-8 nodes x <=2 h, budget-capped | 120 | 395.5 |
| F6 Full shock reproduction, long nonlinear saturation and controlled scaling | All correctness gates closed | **Blocked:** not authorized under the current debug-only/non-production policy; revise plan and obtain approval before submission | 0 | 395.5 |
| Reserve | Unexpected defects or compliant debug-QOS reruns only until a future plan revision explicitly reallocates it | Must be justified in ledger/plan before use | 604.5 | 1000 |

This envelope is not permission to spend the allotment. Stop as soon as enough
evidence exists to decide a gate. If a one-node test fails, fix it before
allocating more nodes.

The reserve is not implicit authorization for production work. Any future
approved production envelope must be added as a new revision with its own
authorization record and cumulative accounting that preserves the 1000
node-hour cap unless the user explicitly changes that cap.

## Detailed Qualification Matrix

Future agents should expand the table with commit hashes, artifact locations
and measured values as work progresses.

| Gate | Capability | Required test/evidence | Status at initial review |
| --- | --- | --- | --- |
| Q-001 | Clean provenance | Curated baseline, manifest, executable checksum | Open |
| Q-002 | Physical-mode configuration | Parser rejection and metadata tests | Open |
| Q-003 | Relativistic pusher and `C` | Section 5.1 analytic convergence CPU/GPU | Blocked by PIC-P0-002 |
| Q-004 | Paper coupling | Conservation tests and no direct CT-current induction | Blocked by PIC-P0-001 |
| Q-005 | Bell | Section 5.2 dispersion and phase | Blocked by Q-003/Q-004 |
| Q-006 | Oscillation plus grid refinement | Section 5.3 uniform/SMR/AMR frequency | Blocked by Q-004 and AMR audit |
| Q-007 | Delta-f | Weight/deposition tests, CRSI/CRPAI reproduction | Blocked by PIC-P0-003 |
| Q-008 | Expanding box | CPAW, gyro and driven CRPAI tests | Blocked by PIC-P0-004 |
| Q-009 | AMR lifetime and boundary policy | Forced rebuild/restart/interface tests | Blocked by PIC-P0-005/PIC-P1-001 |
| Q-010 | Load balance/performance | Particle-cost implementation and Frontier measurements | Open |
| Q-011 | Shock | Section 5.4 coarse/fine/AMR reproduction | Blocked by Q-003/Q-004/Q-009 |
| Q-012 | Reliability aggregate report | Closes after Q-024/Q-036 plus MPI/decomposition/restart-continuation/boundary evidence | Partial scaffolding only |
| Q-013 | Usability/docs | Supported-mode docs and runbook tested by clean launch | Open |
| Q-014 | Terminal sign-off aggregate | After every selected-profile gate closes, freeze a review-ready aggregate manifest with `pending_terminal_review`; then review, record the terminal disposition, sign and archive the final manifest. Never list Q-014 as its own child | Open |
| Q-015 | Entity-derived shared-kernel comparison | Frozen source/formula provenance, relativistic pusher microtest and bounded trajectory-deposition differential checks | Open; existing Entity-named deposit regressions are partial scaffolding only |
| Q-016 | Particle provenance and spectra | Persistent tracking/cohort metadata plus in-run/offline spectral agreement through restart/migration | Open |
| Q-017 | Performance observability and Frontier communication choice | Stage timers, memory/load telemetry and recorded GPU-aware MPI A/B decision | Open |
| Q-018 | Claims registry | Stable claim IDs, evidence links, limitations and reviewer dispositions for every release statement | Open |
| Q-019 | Bell nonlinear saturation | Amplification, spectra, morphology, energy partition, saturation and sensitivity matrix | Blocked by Q-003/Q-004/Q-005 |
| Q-020 | CRSI nonlinear saturation | True delta-f spectra, distribution evolution, scattering/diffusion, saturation, sensitivity matrix, matched reduced nonlinear full-f controls and archived weight-validity envelope | Blocked by Q-007 |
| Q-021 | Driven CRPAI nonlinear saturation | Driven-box branch evolution, anisotropy, scattering, saturation, sensitivity matrix, effective-damping trends, matched reduced nonlinear full-f controls and archived weight-validity envelope | Blocked by Q-007/Q-008 |
| Q-022 | Independent-comparison preregistration | Frozen references, equation mappings, normalization, observables, tolerances, uncertainty methods and discrepancy-ledger schema | Open; prerequisite for comparison runs |
| Q-023 | Statistical qualification | Preregistered estimators, seeds, scans, intervals, exclusions and independent metric recomputation | Open |
| Q-024 | Resilience aggregate report | Forced walltime, interrupted output, corrupt checkpoint, restart schema and bounded storage-failure tests; closes after Q-036 | Blocked by PIC-P0-006 |
| Q-025 | Portability aggregate report | Debug/Release, warnings, host memory checking, declared CPU/OpenMP/MPI scope and Frontier HIP/MPI matrix; closes after Q-038 | Open |
| Q-026 | Archive integrity, licensing and sensitive-data handling | Pre-ingestion and pre-export classification, scrubbed submission artifacts, checksums, redistribution basis, manifest schema and fresh-directory restore drill | Open |
| Q-027 | Frontier debug-QOS authorization boundary | Short-run qualification ledger complete; prohibited production work remains blocked pending revised authorization | Open |
| Q-028 | Independent Bell nonlinear comparisons | Non-Hall paper-mode campaign compared against Bai et al., Riquelme-Spitkovsky, Gargaté et al. and Zacharegkas et al. where regimes overlap: amplification, wavelength evolution, spectra, cavities, filaments, energy transfer and saturation time/mechanism | Blocked by Q-019/Q-022/Q-023 |
| Q-029 | Hall-extension qualification or exclusion | Separately named derived mode with linear/nonlinear Bell and shock-front tests. Record exactly one outcome: `implemented_and_qualified` or `excluded_as_unsupported`; only the first qualifies an extension release | Unsupported unless extension is implemented |
| Q-030 | Independent matched-code comparisons | Close explicit sub-gates Q-030-A Athena/Bai matched Bell-shock observables, Q-030-P PLUTO matched conservative coupling and Q-030-M MPI-AMRVAC overlapping AMR-shock observables, with frozen mappings, quantitative residuals and discrepancy reports | Blocked by Q-004/Q-009/Q-011/Q-019/Q-022/Q-023 |
| Q-030-A | Athena/Bai matched Bell-shock comparison | Frozen Athena/Bai mapping, matched Bell and shock observables, quantitative residuals and discrepancy report | Blocked by Q-011/Q-019/Q-022/Q-023 |
| Q-030-P | PLUTO matched conservative-coupling comparison | Frozen PLUTO mapping, matched conservative-coupling observables, quantitative residuals and discrepancy report | Blocked by Q-004/Q-022/Q-023 |
| Q-030-M | MPI-AMRVAC overlapping AMR-shock comparison | Frozen MPI-AMRVAC mapping, overlapping AMR-shock observables, quantitative residuals and discrepancy report | Blocked by Q-009/Q-011/Q-022/Q-023 |
| Q-031 | Bai et al. 2019 nonlinear CRSI comparison | Growth, saturation, spectra, diffusion, pitch-angle evolution, 90-degree crossing, isotropization and resolution dependence | Blocked by Q-020/Q-022/Q-023 |
| Q-032 | Ion-neutral-damped CRSI extension or exclusion | Plotnikov et al. matched comparison with damping-rate dependence. Record exactly one outcome: `implemented_and_qualified` or `excluded_as_unsupported`; only the first qualifies an extension release | Unsupported unless extension is implemented |
| Q-033 | CRPAI transport-calibration extension or exclusion | Adaptive-delta-f, physical-damping, `nu_eff`, anisotropy, spectra, quasi-steady-state and scaling comparison against Sun, Bai and Zhao. Record exactly one outcome: `implemented_and_qualified` or `excluded_as_unsupported`; only the first qualifies an extension release | Unsupported unless extension is implemented |
| Q-034 | Unsupported-capability review | Execute each register row's `parser_reject`, `rename_proxy`, `docs_limit` or `extension_gate` handling and archive its closure artifact | Open |
| Q-035 | Scoped state-of-the-art sign-off | Reviewer-approved wording tied to qualified claim IDs, exact references, regimes, metrics, performance evidence and exclusions | Blocked by Q-018/Q-023/Q-034/Q-040 |
| Q-036 | Crash consistency and I/O failure handling | Atomic restart publication, completion markers, last-known-good recovery, timeout, interrupted-write, truncated-restart and write-failure tests | Blocked by PIC-P0-006 |
| Q-037 | Durable evidence export | Orion-to-nearline mirror plus institutional or approved off-site long-term archive transfer, checksum verification, retention owner, migration procedure and restore drill | Open |
| Q-038 | Frontier environment profile selection | Minimum-supported baseline plus controlled XNACK/OFI A/B evidence, redacted allowlist capture and restricted-submit-export rank-0 propagation smoke test | Open |
| Q-039 | GIZMO/RSOL comparison decision | Before P7A interpretation, archive a bounded Ji-Hopkins/GIZMO comparison with preregistered metrics or a reviewer-approved documented exclusion from qualification scope | Open |
| Q-040 | Independent-comparison aggregate closure | For non-Hall production Bell close Q-028; for nonlinear CRSI close Q-031; for matched code comparisons close Q-030-A, Q-030-P and Q-030-M. Extension exclusions are allowed only for Q-029/Q-032/Q-033 and the bounded Q-039 GIZMO decision. Archive the discrepancy ledger and scoped conclusions | Blocked by Q-022 and the named child gates |
| Q-041 | Undriven CRPAI nonlinear saturation | Undriven branch evolution, anisotropy, spectra, scattering, saturation, effective-damping trends, sensitivity matrix, matched reduced nonlinear full-f controls and archived weight-validity envelope, distinct from Section 5.6 linear reproduction and physical-damping calibration | Blocked by Q-007 |

## Immediate Agent Handoff: First Actions

Future implementation agents should execute the following sequence:

1. Read this plan, the paper source and all `PIC-P0-*` evidence locations.
2. Confirm the current `PIC` branch, commit and clean-worktree state. Treat the
   retired large-machine note as historical evidence only.
3. Create or update the findings ledger, claims registry and provenance
   manifest framework.
4. Reclassify publication-named proxies in manifests, READMEs, figures and
   reports before generating new qualification artifacts.
5. Write the paper-mode equations/interface specification and obtain review
   before changing more coupling code.
6. Implement `C`-aware relativistic momentum pushing and its analytical tests.
7. Replace or isolate the current current-to-CT field coupling and implement
   conservative paper feedback; run analytical coupling tests.
8. Implement true delta-f and full expanding-box physics with oracle tests.
9. Resolve AMR ownership, buffer-refresh, deposition-policy and load-cost
   correctness.
10. Implement atomic restart publication and fail-closed I/O handling, then
    pass the crash-consistency failure matrix.
11. Only after all small analytical gates pass, build and run Frontier validation
   under the debug-only, budget-tracked procedure above.
12. Reproduce every paper result, execute registered nonlinear and cross-code
    campaigns within authorization, and complete the production sign-off bundle.
13. If long saturation, full shock or scaling work exceeds compliant debug-QOS
    use, stop and obtain a revised approved production-campaign policy rather
    than chaining jobs or treating pilots as final evidence.

## Initial Change Log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-05-25 | Created comprehensive production-readiness and paper-reproduction plan | Expanded narrow large-machine validation scope after full code/paper review exposed model-level blockers and Frontier execution requirements |
| 2026-05-25 | Added Entity Toolkit comparative audit and bounded adoption gates | Adopt shared particle verification, provenance, spectra, restart and observability practices without importing incompatible electromagnetic PIC physics into paper mode |
| 2026-05-30 | Retired the narrow large-machine note as historical-only evidence; corrected AMR wording; added claims registry, nonlinear saturation, cross-code, statistical, resilience, portability, archive/licensing and Frontier authorization-boundary gates | A fresh source and plan audit showed that paper reproduction plus exploratory proxies alone could not support production or scoped state-of-the-art claims |
| 2026-05-30 | Enforced engineering-proxy classification in checked-in artifact tooling and hardened the plan after independent physics, taxonomy and Frontier-operability reviews | Legacy publication identifiers now carry explicit unqualified metadata and visible figure labels; comparison gates, sampling rules, release profiles and fail-closed submission accounting are executable specifications rather than implied policy |
| 2026-05-30 | Added immutable Frontier control-plane snapshots, JSONL accounting with non-recursive mirror receipts, dual Slurm-ID guards, inert scan shells, root-document retirement banners and checksummed standalone-helper lineage companions | Adversarial execution review found that documentation-only policy was insufficient while mutable paths, ambient exports, copied shell scripts or partially checksummed quicklooks could bypass provenance and submission boundaries |

## Source Pointers For The Initial Review

These pointers identify the initial audit basis. Line numbers may move as fixes
are implemented; future findings must record the commit examined.

| Topic | Initial source locations |
| --- | --- |
| Coupled E-field and feedback source paths | `src/mhd/mhd_tasks.cpp` |
| Particle controls and artificial-light-speed guard | `src/particles/particles.cpp` |
| Particle state and ownership | `src/particles/particles.hpp` |
| Non-relativistic Boris and box scaling | `src/particles/particles_pushers.cpp` |
| Moment deposition and particle weighting | `src/particles/particles_moments.cpp` |
| Particle/MHD task ordering | `src/particles/particles_tasks.cpp` |
| AMR child MeshBlock/coordinate reconstruction and retained-module refresh | `src/mesh/mesh_refinement.cpp` |
| Particle shock generator | `src/pgen/tests/pic_parallel_shock.cpp` |
| Restart publication and I/O error handling | `src/outputs/restart.cpp`, `src/outputs/io_wrapper.cpp` |
| Internal walltime stop and final output | `src/driver/driver.cpp` |
| Existing particle and publication tests | `inputs/tests/pic_*.athinput`, `tst/scripts/particles/pic_*.py`, `tst/publication/` |
| Paper manuscript source | `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` |
| Entity comparative documentation | <https://entity-toolkit.github.io/wiki/> |
| Entity local source snapshot used for comparative audit | `/Users/dbf75/Work/Research/AthenaK/entity` at `a59065fc`; especially `src/framework/containers/particles.h`, `src/kernels/particle_pusher_sr.hpp`, `src/kernels/currents_deposit.hpp`, `src/kernels/tests/pusher.cpp`, `src/kernels/tests/deposit.cpp`, `src/kernels/digital_filter.hpp`, `src/framework/domain/output.cpp`, `src/framework/domain/stats.cpp`, `src/output/checkpoint.cpp`, `src/engines/srpic.hpp` and `src/engines/engine_run.cpp` |

This document is the entry point for all future MHD-PIC production work. A
simulation result is not a production result until its governing mode, oracle,
provenance, resource accounting and archived evidence satisfy this plan.
