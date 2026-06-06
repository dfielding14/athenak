# AthenaK MHD-PIC Production Readiness and Sun & Bai (2023) Reproduction Plan

## Document Status

| Item | Value |
| --- | --- |
| Purpose | Canonical implementation, verification, publication-reproduction, and Frontier qualification plan for the AthenaK MHD-PIC module |
| Primary reference | Sun & Bai, *The Magnetohydrodynamic-Particle-In-Cell Module in Athena++: Implementation and Code Tests*, source manuscript in `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` |
| Initial review baseline | Historical local `c/pic-review` checkout at base commit `ca7e43f0d0d6`, reviewed with an extensive dirty working tree containing in-progress PIC changes |
| Plan-revision baseline | Clean `PIC` branch at commit `3bcd3f21e4f5439e655d977854f03b287fa814b9`, tracking `origin/PIC` |
| Production verdict at review time | **Not production-ready. Paper reproduction is blocked by model-level implementation gaps.** |
| Supersedes | All scope, terminology, pass criteria, operational instructions and sign-off language in `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md`; that retired note is historical evidence only |
| Required execution root on Frontier | Run every Frontier simulation under `/lustre/orion/ast207/proj-shared/dfielding/PIC`. No alternate simulation directory is authorized. |
| Frontier scheduling constraint | Prefer the `debug` QOS only for eligible short non-production work when the user has no `debug` job in any state; otherwise use the `normal` QOS on the `batch` partition. Submit at most one AthenaK PIC job at a time across both QOS classes. |
| Total Frontier testing budget | 10000 node-hours maximum, tracked before and after every submitted job. Stop and ask the user for permission before increasing the cap. |
| Authorized target profiles | Complete the production-ready `paper_mhd_pic` profile and implement plus qualify the Hall Bell, ion-neutral-damped CRSI and adaptive-delta-f physical-damping CRPAI `extended_mhd_pic` capabilities. Explicit exclusion is not sufficient for those three selected extension gates. |
| Review disposition placeholder | Record gate and terminal dispositions as `pending external review` until a named external reviewer is assigned. No pending disposition may be promoted to `qualified`. |
| OLCF-side evidence strategy | Use `/lustre/orion/ast207/proj-shared/dfielding/PIC` as the user-selected sole bulk-evidence root for simulation output, immutable sign-off bundles, private reference-artifact staging and restore-drill evidence. Use backed-up AST207 Project Home only for the small append-only ledger/control-plane mirror with `filesystem_copy`. Kronos, DTN and Globus are removed from the execution boundary. Orion-only retention carries an explicit durability risk and must not be represented as an institutional archive. |

This plan is written as a production gate, not as a list of desirable enhancements.
Agents must not declare the MHD-PIC implementation production-ready, nor claim
reproduction of Sun & Bai (2023), until every blocking gate below is satisfied
with archived evidence.

## Executive Verdict

The initial review found four immediately disqualifying implementation defects:
paper and Hall induction were conflated, particles used non-relativistic
velocity state, delta-f was quiet-start only, and expanding-box support scaled
particles without the corresponding MHD map. The current worktree lands bounded
repairs for those defects and adds isolated host oracles. It is still not
production-ready: paper reproduction, extension qualification, portability,
resilience, archival, and external-review gates remain open.

The current implementation now provides:

1. Explicit `paper_mhd_pic` and `extended_mhd_pic` identities. Paper mode keeps
   ideal-MHD induction, while the opt-in Hall-current CT source remains an
   experimental extension with a host odd-in-`alpha_H` source-isolation oracle.
2. Mass-normalized momentum state, configurable artificial light speed `C`,
   relativistic Boris mechanics, relativistic kinetic-energy feedback, and a
   host gyro/restart smoke.
3. Separately named quiet-start and evolving physical delta-f paths, weighted
   deposition, VTK diagnostics, restart fingerprints, plus an extension-only
   adaptive global bi-kappa fit with exact host fit and restart-state checks.
4. Particle and MHD expanding-box transforms with explicit scale-factor laws,
   uniform-MHD invariants, and an extension-only reduced static-neutral
   ion-neutral friction map with an exact manufactured-source oracle.

High-risk work remains. AMR ownership refresh and the explicit `paper_smooth`
refinement-interface policy have bounded serial coverage but still need the
repeated MPI/GPU lifetime matrix. The selected physical-damping adaptive-delta-f
CRPAI composition now passes a bounded serial-host endpoint-source and
damping-order oracle; physical campaign qualification remains open.
Particle-aware load cost exists but needs Frontier measurement. Atomic restart
publication, checksummed completion markers, local corruption guards, bounded
short-write/seek/storage-failure injection, interrupted-writer preservation, and
soft wallclock continuation parity pass locally; per-rank MPI, node-loss,
scheduler-pretimeout, and Frontier-filesystem drills remain open.

The bounded local analytical, conservation, restart and portability tranches
are archived. The reviewed Orion-only storage policy, paired immutable control
plane and mirrored ledger genesis were initialized by predecessor snapshots.
The historically installed `6f3458ca` control-plane generation is chronology
only and must not launch a new Frontier job. New Frontier submissions remain
prohibited until the reviewed successor is installed on both paired roots,
promoted, bound to a canonical clean-candidate freeze and authorized by a
campaign-specific immutable registration. After that transition, run only
registered, budget-tracked Frontier campaigns.

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
6. Run every Frontier simulation under
   `/lustre/orion/ast207/proj-shared/dfielding/PIC`. No alternate simulation
   directory is authorized. All build, log, metric, ledger, runtime-status and
   artifact paths must remain below that same root.
7. On Frontier, prefer `debug` only for eligible short non-production work when
   the user has no `debug` job in any state. Otherwise submit to the `normal`
   QOS on the `batch` partition. Keep AthenaK PIC submissions serial across both
   QOS classes and never exceed the cumulative 10000 node-hour budget without
   stopping and obtaining explicit user permission.
8. Treat `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md` as retired historical
   evidence only. It must never authorize a run, close a gate, or define a
   scientific claim.
9. Do not use the phrase `state of the art` without naming the compared methods,
   codes, physical regimes, metrics, reference versions and limitations. A
   scoped comparative conclusion is required instead.
10. Do not start a production-scale Frontier campaign under the `debug` QOS.
    Use `normal` on `batch` for registered work that is not eligible for
    `debug`, including longer paper-reproduction, nonlinear-saturation and
    controlled-scaling jobs after their prerequisite gates close. Any request
    for another partition, QOS, account policy or a budget above 10000 node-hours
    requires an explicit plan revision and user authorization.
11. Make PIC implementation work fit AthenaK's existing style framework and
    design choices wherever possible. Depart only when the AthenaK-native design
    would materially reduce physical accuracy or computational efficiency, and
    document the evidence, scope and review disposition for every exception.

## AthenaK Style And Design Compatibility Contract

PIC code must read as an AthenaK module, not as a separately designed codebase
embedded inside AthenaK. This requirement covers architecture, interfaces,
execution patterns, naming, diagnostics, documentation and formatting.

### Default Implementation Rules

1. Study the nearest AthenaK implementation before editing. Reuse established
   `MeshBlockPack`, module-constructor, `ParameterInput`, `TaskList`, boundary
   exchange, output, restart and problem-generator patterns where they satisfy
   the physical contract.
2. Keep hot kernels in the established Kokkos execution model. Match nearby
   `View` ownership, layout, `par_for`, host/device synchronization and MPI
   staging patterns unless measured evidence requires a scoped alternative.
3. Prefer existing AthenaK ownership boundaries and helper APIs over new PIC-only
   abstractions. Add an abstraction only when it removes demonstrated complexity
   or is required for accuracy, reliability or performance.
4. Keep runtime parameters explicit, narrowly scoped and validated at startup.
   Follow nearby AthenaK naming and failure-message conventions. Reject unsafe
   combinations rather than silently selecting a fallback.
5. Integrate outputs, diagnostics, timers and restart state through AthenaK's
   existing extension surfaces. Do not introduce parallel ad hoc formats when an
   established output or manifest mechanism can represent the required data.
6. Keep changes incremental and reviewable. Separate physical-model changes,
   performance changes and broad refactors unless they cannot be validated
   independently.
7. Treat Entity Toolkit and external-code patterns as comparison inputs, not as
   automatic design replacements. Translate a useful idea into AthenaK's
   conventions unless a documented exception is justified.
8. Follow repository formatting rules for every touched source file: no tabs,
   90-character C++ line limit, one closing brace per line, concise comments and
   ASCII text unless an existing file requires otherwise.

### Exception Policy

An AthenaK-style deviation is allowed only when the default design would
materially reduce physical accuracy or computational efficiency. Before merging
an exception:

1. Record the rejected AthenaK-native design, the proposed deviation and the
   smallest affected surface in the findings ledger.
2. For an accuracy exception, archive an analytical, manufactured or
   independently reproduced comparison showing the material accuracy loss.
3. For an efficiency exception, archive matched CPU and relevant GPU benchmarks
   showing the material runtime, scaling or memory regression. Do not accept a
   speculative performance argument.
4. Confirm that the exception preserves restart, AMR, MPI, GPU, diagnostic and
   usability contracts.
5. Obtain code review and record the reviewer disposition in the release
   manifest. Revisit the exception if AthenaK's common infrastructure changes.

### Style And Architecture Verification

For each implementation batch:

1. Run targeted formatting and lint checks on every touched C++ and Python file.
2. Run `bash tst/scripts/style/check_athena_cpp_style.sh`, record the existing
   repository-wide baseline debt separately, and introduce no new violations.
3. Run the relevant regression suite and the smallest representative GPU test
   for any hot-kernel or data-layout change.
4. Review the diff specifically for avoidable PIC-only infrastructure,
   inconsistent ownership, duplicated helpers, ambiguous parameters and
   deviations from nearby AthenaK patterns.
5. Add any approved exception and its archived evidence to the Q-042
   architecture/style-conformance report.

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
revision, evidence links and reviewer. The initial reviewer is
`pending external review`; no row
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
| `CLAIM-EXT-HALL-BELL-001` | Qualify a CR-induced Hall Bell extension | Separately named `extended_mhd_pic` mode only | Q-029 | `open`; implementation and qualification required by authorized scope |
| `CLAIM-EXT-CRSI-IN-DAMPING-001` | Qualify ion-neutral-damped CRSI | Separately named extension; excludes ordinary paper-mode claim | Q-032 | `open`; implementation and qualification required by authorized scope |
| `CLAIM-EXT-CRPAI-TRANSPORT-001` | Calibrate CRPAI saturated-state effective scattering with physical damping | Separately named adaptive-delta-f and physical-damping extension | Q-033 | `open`; implementation and qualification required by authorized scope |
| `CLAIM-XCODE-ATHENA-BELL-SHOCK-001` | Compare matched Bell and shock observables against the frozen Athena/Bai reference family | Matched observable subset; Hall-sensitive regimes must be separated | Q-022, Q-028, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-BELL-SATURATION-001` | Compare nonlinear Bell saturation against independent kinetic and hybrid-PIC literature | Bounded parameter overlap; report model differences explicitly | Q-022, Q-028, Q-040 | `open` |
| `CLAIM-XCODE-BAI2019-CRSI-001` | Compare nonlinear CRSI against Bai et al. (2019) | `paper_mhd_pic` true delta-f overlap only | Q-022, Q-031, Q-040 | `open` |
| `CLAIM-XCODE-PLUTO-COUPLING-001` | Compare matched conservative-coupling behavior against PLUTO | Matched-equation subset only | Q-022, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-AMRVAC-SHOCK-001` | Compare overlapping AMR-shock observables against MPI-AMRVAC publications | Matched observable subset only | Q-022, Q-030, Q-040 | `open` |
| `CLAIM-XCODE-GIZMO-RSOL-DECISION-001` | Decide whether a bounded Ji-Hopkins/GIZMO RSOL comparison is informative, or exclude it with rationale | Do not conflate GIZMO RSOL equations with the Sun and Bai artificial-`C` contract | Q-039 | `open` |
| `CLAIM-XCODE-ENTITY-MICRO-001` | Compare bounded shared particle-numerics micro-oracles against frozen Entity Toolkit source | Shared-kernel pusher/deposition/restart subset only; excludes EM evolution and gas feedback | Q-015, Q-022 | `open` |
| `CLAIM-STATEART-CRPAI-SCATTERING-001` | Make a scoped contemporary CRPAI scattering comparison for the selected extension | Requires physical damping, adaptive delta-f comparison and explicit parameter envelope | Q-033, Q-035 | `open`; implementation and qualification required by authorized scope |
| `CLAIM-BUNDLE-SUN-BAI-2023-001` | Qualify the complete, paper-faithful Sun and Bai (2023) reproduction bundle | Exactly the `sun_bai_2023_reproduction` profile; no extension substitutions and no production-readiness implication | All `sun_bai_2023_reproduction` profile gates, then terminal Q-014 disposition recorded in the signed final manifest | `open` |
| `CLAIM-RELEASE-PAPER-MHD-PIC-001` | Qualify the bounded production-ready AthenaK `paper_mhd_pic` release | Exactly the documented paper-mode applicability envelope; unsupported extensions remain excluded | All Production-ready `paper_mhd_pic` profile gates, then terminal Q-014 disposition recorded in the signed final manifest | `open` |
| `CLAIM-RELEASE-EXTENDED-MHD-PIC-001` | Qualify the authorized AthenaK `extended_mhd_pic` release bundle | Production-ready `paper_mhd_pic` plus implemented-and-qualified Hall Bell, ion-neutral-damped CRSI and adaptive-delta-f physical-damping CRPAI extensions; exclusion is not sufficient | Production-ready `paper_mhd_pic` profile gates, Q-029, Q-032, Q-033, then terminal Q-014 disposition recorded in the signed final manifest | `open`; implementation and qualification required by authorized scope |

### Explicit Unsupported-Capability Register

Unsupported claims must remain visible until a separately named implementation
and its qualification gates are complete.

| Capability or claim | Handling | Closure artifact |
| --- | --- | --- |
| Full electromagnetic PIC physics | `rename_proxy`, `docs_limit` | Reclassified manifest, figure watermark and known-limitations entry; current EM-vacuum, two-stream and Weibel names must not imply this capability |
| Self-consistent Langmuir-wave physics | `rename_proxy`, `docs_limit` | Rename the uniform-`B` orbit-frequency anchor and document its narrow oracle |
| Physical two-stream or Weibel fidelity | `rename_proxy`, `docs_limit` | Reclassified engineering-proxy manifest and figures |
| Hall-dominated Bell or shock-front behavior | `extension_gate` | Implement and qualify `CLAIM-EXT-HALL-BELL-001`; reject Hall claims until that selected extension closes |
| Injection from the thermal pool | `docs_limit` | Section 5.4 reproduction manifest documents its simplified prescription |
| Relativistic MHD background fluid | `parser_reject`, `docs_limit` | Parser test and supported-mode table |
| Physical damping and calibrated CR transport coefficients | `extension_gate` | Implement and qualify the selected physical-damping extensions; reject calibrated-transport claims until their gates close |
| Oblique-shock generality | `parser_reject`, `docs_limit` | Parser test or unsupported-mode documentation until independently implemented |
| Frontier work outside `debug`-preferred, `normal`-fallback scheduling on `batch`, or beyond the 10000 node-hour cap | `docs_limit`, `extension_gate` | Blocked disposition pending explicit revised execution authorization and user permission |

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
| Historical Entity checkout at unavailable local commit `a59065fc`; frozen accessible replacement `/ccs/home/dfielding/entity` at commit `512998c471bf3fdec292cb4a64150c4f0aeea539`, tree `38d511d3d317d0866fab33201a7355815c26bedf`, archive SHA-256 `c92f5fac17cdbf90fc0c6adcc0b9445fa5b9f4f78c1b6a53c97aa6520423ce1d` and manifest `tst/publication/readiness/entity_snapshot.json` | Source-level shared particle storage/pusher, deposition, filtering, output, timer, parameter and restart comparisons may use only the frozen bounded subset; Entity electromagnetic evolution remains out of scope as an AthenaK MHD-PIC oracle |
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
- Frontier HIP/MPI short-run validation and registered `normal`-QOS controlled
  scaling results after prerequisite gates close.
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
| Section 5.1 gyro-motion | Relativistic orbit, energy and phase accuracy for the paper's `C`, `v0`, `Omega`, and timestep constraints | Momentum/`C` mechanics and bounded host analytical convergence pass | Freeze paper deck; run clean-candidate GPU and portability matrix |
| Section 5.2 Bell instability | Measured phase and growth rate versus analytical dispersion in 1D, 2D, and 3D | Paper coupling isolation passes; existing Bell deck remains an engineering proxy | Freeze paper-faithful decks and add analytical dispersion, convergence, MPI and GPU comparison |
| Section 5.3 gas plus electron/positron oscillation | Correct oscillation frequency and equivalent behavior on uniform, SMR, and AMR grids | Coupling conservation and serial AMR refresh smoke pass; existing oscillation fixtures remain proxies | Freeze long-horizon paper decks; close frequency, AMR-policy, MPI and GPU matrix |
| Section 5.4 non-relativistic parallel shock acceleration | Paper domain, injection, AMR, spectra, morphology, timing and load-balancing behavior | Isotropic shock-surface injection scaffold and restart fingerprints landed; a frozen Section 5.4 preparation-only deck and fail-closed analyzer contract now reject incomplete calibration, while the executed fixtures remain engineering proxies | Complete calibration and injection audit, bind provenance/spectra outputs, then run the clean GPU pilot and registered paper campaign |
| Section 5.5 CR gyro-resonant streaming instability | Polarization-resolved spectra and growth rates against analytical prediction with true delta-f | True delta-f mechanics, deposition and restart guards pass | Freeze theory/spectral oracle and run clean MPI/GPU campaign |
| Section 5.6 CR pressure anisotropy instability | Prolate/oblate polarization branch selection and quantitative growth against theory | True delta-f and relativistic mechanics pass bounded host tests | Freeze branch theory oracle and run clean MPI/GPU campaign |
| Section 5.7 driven expanding/compressing-box CRPAI | Correct anisotropy evolution, spectra and distribution evolution under box driving | Bounded comoving-flux box repair, the narrow active-MHD endpoint-normalized adaptive-source plus reduced-damping launch oracle and serial uninterrupted-versus-restarted endpoint parity pass; physical transport calibration remains open | Close clean-candidate MPI/Frontier restart slices and qualify the driven physical campaign |
| Appendix A circularly polarized Alfven wave | Analytic expansion/compression amplitude and phase response | Bounded one-step expanding/compressing CPAW closure and oblique-divB oracles pass | Run full history and resolution-convergence campaign with changing-volume accounting |
| Appendix B expanding-box gyro-motion | Analytic gamma/phase history with prescribed expansion rate | Bounded relativistic expanding-box history convergence passes | Freeze manuscript deck and run clean portability matrix |
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
No rendered PDF is checked into Git. The bundled archive omits `mnras.cls` and
`mnras.bst`, although its `readme.txt` lists both files. On 2026-05-30 a private
reconstruction retrieved and checksummed the CTAN MNRAS class, bibliography
style and TeX dependencies, rendered a 20-page PDF with SHA-256
`86e1f53ee335685a833b65b000e105d72c83c9bcea3a7b4a0ccaaa4deabaf3d6`,
archived the private Orion bundle and passed a fresh-directory checksum
restore drill. The bundle manifest SHA-256 is
`6103d0da9c5782e22b702fd7ecf50cc04ae75126b5d11e1aa12bea41684611ed`;
see `tst/publication/readiness/q026_private_render_restore_2026-05-30.json`.
External redistribution review and the external disposition on user-selected
Orion-only retention risk remain open.

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

- Exact git commit and tree, clean-candidate manifest checksum, parent source
  archive checksum, canonical recursive-submodule source-bundle checksum,
  structured clean pinned submodule archive attestations, CMake cache,
  executable checksum, loaded modules, reviewed redacted allowlisted
  environment capture, Slurm script, input deck and analysis script. Never
  archive unrestricted environment dumps.
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
| PIC-P0-001 | P0 | Initial defect: engineering coupling added deposited current directly to CT while paper mode requires ideal-MHD induction and no CR Hall source. Explicit paper and extension identities are now separated; paper mode suppresses direct CT current and extension Hall is opt-in | `src/mhd/mhd_tasks.cpp`; `src/particles/particles.hpp`; paper Section 2 | Implementation landed; conservation, Bell, oscillation, GPU and decomposition qualification remain open | Pass Q-004 |
| PIC-P0-002 | P0 | Initial defect: pusher advanced non-relativistic velocity and energy and rejected configurable `C`. Explicit paper/extension modes now store `p/m`, derive relativistic velocity and kinetic energy, and apply a `C`-aware Boris rotation | `src/particles/particles_pushers.cpp`; `src/particles/particles.cpp`; `tst/scripts/particles/pic_relativistic_gyro_paper.py`; paper Sections 2 and 5 | Host analytical smoke and restart parity pass; convergence, GPU and portability qualification remain open | Pass Q-003 |
| PIC-P0-003 | P0 | Initial defect: `pic_deltaf_mode=on` was quiet-start only. Physical paper/extension mode now stores initial `f0`, evolves perturbation weights, deposits weighted perturbation moments, fingerprints restart state and emits VTK diagnostics | `src/particles/particles.cpp`, `particles_moments.cpp`, `particles_pushers.cpp`; `tst/scripts/particles/pic_relativistic_gyro_paper.py`; paper Section 4.4 and Sections 5.5-5.7 | Mechanics smoke passes; CRSI/CRPAI theory, spectra and nonlinear qualification remain open | Pass Q-007 |
| PIC-P0-004 | P0 | Initial defect: expanding-box code updated particle velocity only. The repaired bounded path applies exact scale-factor particle half steps, comoving drift, one exact post-RK MHD conserved-variable map per cycle, stores raw face arrays as divergence-preserving comoving magnetic fluxes, derives physical MHD views, maps edge EMFs before CT, deposits physical-volume moments, applies non-delta-f conservative feedback in the final physical source frame, provides a separately identified endpoint-normalized adaptive delta-f analytic source path and emits built-in physical-volume MHD history | `src/particles/particles_pushers.cpp`; `src/particles/particles_moments.cpp`; `src/mhd/mhd_tasks.cpp`; `src/mhd/mhd_fluxes.cpp`; `src/outputs/history.cpp`; `tst/scripts/particles/pic_mhd_expanding_box_uniform.py`; `tst/scripts/particles/pic_mhd_expanding_box_coupled_conservation.py`; `tst/scripts/particles/pic_mhd_expanding_box_adaptive_damping_smoke.py`; `tst/scripts/particles/pic_mhd_expanding_box_oblique_divb.py`; paper Section 3 and appendices | Uniform, non-delta-f final-frame conservation, adaptive-delta-f source-normalization and damping-order, oblique-divB and serial uninterrupted-versus-restarted endpoint host oracles pass; unsupported compositions remain fail-closed; full Appendix-A convergence, clean-candidate MPI/Frontier restart slices, driven CRPAI and GPU qualification remain open | Pass Q-008 |
| PIC-P0-005 | P0 | AMR retains the `MeshBlockPack` object while reconstructing child MeshBlock/coordinate objects. Cross-rank particle migration already resolves new geometric ownership; retained particles now receive an explicit post-AMR refresh hook that validates MeshBlock-sized capacities and refreshes the stable pack pointer | `src/mesh/mesh_refinement.cpp`; `src/mesh/load_balance.cpp`; `src/particles/particles.cpp`; particle and mesh boundary-helper construction paths | Serial shock-rich refine smoke, a six-transition serial refine/derefine stress with five restart continuations, coupled PIC/MHD reflecting/outflow Debug-plus-UBSan lifetime stress and repaired x2-inflow Debug-plus-host-ASan/UBSan lifetime stress pass; multi-rank MPI migration/load-balance/restart and HIP evidence remain open | Complete Q-009 matrix |
| PIC-P0-006 | P0 | Initial defect: restart output wrote directly to its final pathname and ignored close failures. Restart output now writes `.partial`, checks every sequential header append and offset seek/write, syncs, closes, atomically promotes, writes checksummed completion markers and manifests, and validates completed artifacts before consumption. Schema-2 rank-sharded checkpoints add a shared nonzero best-effort nonce to detect accidental mixed-shard assembly; schema-1 and schema-absent legacy reads intentionally remain supported with the weaker pre-nonce guarantee. | `src/outputs/restart.cpp`; `src/outputs/io_wrapper.cpp`; `src/outputs/restart_utils.cpp`; `src/main.cpp`; `tst/scripts/particles/pic_restart_safety_guards.py`; `tst/scripts/particles/restart_fault_injector.c` | Serial-host parity, checksum/payload/schema rejection, short-write, seek-failure, `/dev/full`, interrupted-writer and soft-wallclock continuation drills pass; MPI per-rank, node-loss, scheduler-pretimeout and Frontier-filesystem pilots remain open | Pass Q-036 |
| PIC-P1-001 | P1 | Refinement-boundary deposition policy has not been reconciled with the paper's smoothness-versus-conservation choice | `src/particles/particles_moments.cpp`; paper Section 4.3 | AMR results may be smooth but non-reproducing, or conservative but physically different from the paper | Specify policies, validate both errors, select paper policy for reproduction |
| PIC-P1-002 | P1 | Task-stage ordering and conservative delta exchange are not yet proved against the paper's second-order method | `src/particles/particles_tasks.cpp`; MHD source paths; paper Section 2.3 | Good-looking tests may mask order loss or incorrect exchange timing | Create stage-contract tests and convergence/conservation gates |
| PIC-P1-003 | P1 | Shock setup now contains a source-local manuscript ideal-surface model with full-sphere isotropic monoenergetic injection, a separately retained finite-Mach engineering surface option, frame-tracking controls, schema-6 continuation fingerprints and a schema-7 bounded provenance/spectrum successor. A recursively read-only reduced serial-host runtime audit checks actual PVTK shock-injected provenance, clamped surface placement, monoenergetic surface-relative speed and bounded full-sphere statistics; physical calibration remains open | `src/pgen/tests/pic_parallel_shock.cpp`; shock input decks; `tst/scripts/particles/pic_parallel_shock_restart_controls.py`; `tst/scripts/particles/pic_q016_particle_provenance.py`; `tst/publication/analyze_q011_injection_distribution_runtime_local.py`; paper Section 5.4 | Shock acceleration plots cannot yet be publication evidence | Close physical units and macro-mass calibration, then close the preregistered AMR/MPI/GPU campaign |
| PIC-P1-004 | P1 | Sorting/intermediate-array controls remain to be benchmarked. AMR load balancing now supports an opt-in `pic_load_balance_cost_per_particle` term computed from geometric post-AMR particle ownership | particle runtime controls; `src/mesh/mesh_refinement.cpp`; load-balance code; paper Section 4 | Performance and load-balance claims remain unsupported until Frontier measurements close | Benchmark sorting, intermediate-array and particle-cost choices; record selected production settings |
| PIC-P1-005 | P1 | Existing publication-named tests frequently assert sign/parity rather than analytical paper oracles | `tst/scripts/particles/pic_*publication.py` and proxy suites | Passing suite can overstate physical confidence | Reclassify proxies and add quantitative reference tests |
| PIC-P1-006 | P1 | The `pic_entity_deposit_*` tests use Entity-oriented names but currently demonstrate AthenaK integrated-moment and decomposition parity checks rather than a provenance-frozen differential comparison with Entity kernels | `tst/scripts/particles/pic_entity_deposit_mink.py`, `pic_entity_deposit_reflect.py`; Entity `src/kernels/currents_deposit.hpp` and `src/kernels/tests/deposit.cpp` | An apparently external cross-check can be overinterpreted as independent validation | Either rename the tests as internal deposit regressions or add a bounded, frozen Entity-reference comparison for shared deposition behavior |
| PIC-P1-007 | P1 | Orion is a purge-eligible working filesystem rather than durable evidence storage, no OLCF resource is guaranteed as permanent institutional retention, and unrestricted environment capture can preserve sensitive values | OLCF Frontier and storage guidance; historical artifact layout and script templates | Qualification evidence may disappear or archives may expose inappropriate environment data | The user selected Orion-only bulk-evidence retention and removed Kronos from scope. Keep the durability risk explicit, run checksummed Orion restore drills, archive a redacted environment allowlist only, and require external disposition before terminal sign-off |
| PIC-P1-008 | P1 | The opt-in `extended_mhd_pic` Hall-current path now has an explicit experimental normalization and a host manufactured-source oracle: Hall-off and `+/- alpha_H` one-cycle runs produce nonzero magnetic increments with odd residual `1.25e-5` and cosine `-0.9999999999`. This proves source isolation only; the derived CR-Hall normalization and scientific envelope remain open | `docs/source/engineering/pic_mhd_model_contract.md`; `src/mhd/mhd_tasks.cpp`; `tst/scripts/particles/pic_extended_hall_ct_smoke.py`; Bai et al. (2015) | Treating the source smoke as Bell or shock-front qualification would overstate the implementation | Derive and preregister the extension model, then close linear/nonlinear Bell, shock-front, GPU and decomposition matrices under Q-029 |
| PIC-P1-009 | P1 | The selected ion-neutral damping extension now has a bounded `extended_mhd_pic` implementation: final-physical-frame transverse ion momentum is multiplied by `exp(-nu_in dt)` and ideal-MHD energy loses the removed ion kinetic energy. The host oracle measures factor `0.9686827210` with density/longitudinal residual zero and transverse/energy residuals below `1.2e-7` | `docs/source/engineering/pic_mhd_model_contract.md`; `src/mhd/mhd_tasks.cpp`; `tst/scripts/particles/pic_ion_neutral_friction_smoke.py`; Plotnikov et al. (2021) | The static-neutral manufactured source is not a matched damped-CRSI dispersion, saturation, GPU, or decomposition result | Preregister the reduced-model applicability envelope and close the Plotnikov comparison matrix under Q-032 |
| PIC-P1-010 | P1 | The selected adaptive delta-f extension now has an extension-only global bi-kappa moment fit for `xi` and `p0`, fixed cadence, normalized fitted background, restart schema version 7 persistence, and a bounded active-MHD expanding-box endpoint-normalized analytic-source plus reduced-damping launch composition. Deterministic host oracles match both fitted values exactly, preserve uninterrupted-versus-restarted particle state without an unintended refit and pass coupled-carrier endpoint parity | `docs/source/engineering/pic_mhd_model_contract.md`; `src/particles/particles_tasks.cpp`; `src/particles/particles_pushers.cpp`; `tst/scripts/particles/pic_adaptive_deltaf_smoke.py`; `tst/scripts/particles/pic_mhd_expanding_box_adaptive_damping_smoke.py`; `tst/scripts/particles/pic_mhd_expanding_box_adaptive_damping_restart.py`; Sun, Bai and Zhao (2024) | The bounded x1-parallel fit, launch mechanics and serial restart parity do not establish physical-damping CRPAI transport calibration, effective scattering, saturation, GPU, or MPI qualification | Close clean-candidate MPI/Frontier restart slices and the physical-damping, `nu_eff`, anisotropy, spectra, quasi-steady-state and scaling matrix under Q-033 |

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
independent implementation reference. The historical local comparison used
commit `a59065fc`; that object is not present in the accessible replacement
checkout. The replacement source is frozen in
`tst/publication/readiness/entity_snapshot.json` from
`/ccs/home/dfielding/entity` at commit
`512998c471bf3fdec292cb4a64150c4f0aeea539`, tree
`38d511d3d317d0866fab33201a7355815c26bedf`, and clean-tree archive SHA-256
`c92f5fac17cdbf90fc0c6adcc0b9445fa5b9f4f78c1b6a53c97aa6520423ce1d`.
Any future differential test must cite that manifest and remain within its
bounded shared-kernel scope.

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
| `paper_test_particle` | Particle orbit validation in prescribed fluid/field background | No MHD feedback; may support prescribed/ideal background E and B | Any claim of coupled instability or shock backreaction |
| `paper_mhd_pic` | Reproduction of Sun & Bai model | Ideal-MHD induction, conservative CR-to-gas momentum/energy exchange, configurable `C`, optional true delta-f and full box equations | Direct CR-current CT injection; disabled feedback for coupled runs |
| `extended_mhd_pic` | Selected derived extensions: Hall-current induction, ion-neutral-damped CRSI and adaptive-delta-f physical-damping CRPAI transport calibration | Only equations explicitly derived and tested for the selected extension | Being used under paper-reproduction test names or silently substituting for `paper_mhd_pic` |
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
5. Validate checkpoint/restart equivalence in short, compliant `debug` pilots.
   Do not chain a full production shock campaign through the `debug` QOS. Run a
   registered full shock campaign under `normal` on `batch` after prerequisites
   close and only while the cumulative reservation remains within the approved
   node-hour cap.
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
- AthenaK architecture/style-conformance report, including any evidence-backed
  exceptions.
- Archive-integrity, restore-drill and licensing report.
- Cumulative Frontier node-hour ledger, demonstrating budget compliance.
- Final checklist signed in this document with links to immutable artifact paths.

### Release-Profile Closure Matrix

Do not use one profile's closure to imply another. Record the selected profile
in the release manifest and close every mandatory gate listed for it.

| Release or evidence profile | Mandatory gates |
| --- | --- |
| `sun_bai_2023_reproduction` evidence bundle | Q-001 through Q-009, Q-011 through Q-013, Q-016, Q-018, Q-023, Q-025 through Q-027, Q-034, Q-036 through Q-038 and Q-042 |
| Production-ready `paper_mhd_pic` | All `sun_bai_2023_reproduction` gates plus Q-010, Q-017, Q-019 through Q-021, Q-024, Q-028, Q-030, Q-031, Q-040 through Q-042 |
| Scoped state-of-the-art statement | Production-ready `paper_mhd_pic` plus Q-022, Q-028 through Q-035, Q-039 and Q-040; unsupported extensions must close through explicit exclusion rather than omission |
| Authorized `extended_mhd_pic` bundle | Production-ready `paper_mhd_pic` plus `implemented_and_qualified` outcomes for Q-029 Hall Bell, Q-032 ion-neutral-damped CRSI and Q-033 CRPAI transport calibration. `excluded_as_unsupported` does not satisfy the selected release scope |

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
| `pic_crsi_deltaf*` | `engineering_proxy` until quantitative delta-f qualification closes | Evolving weights now exist, but existing cases still lack theory comparison, spectra and nonlinear saturation evidence | Add quantitative theory, spectra and nonlinear saturation |
| `pic_crpai_polarization*` | `engineering_proxy` until quantitative delta-f qualification closes | Branch-sign separation is not quantitative CRPAI reproduction | Add theory comparison, spectra and nonlinear saturation |
| `pic_expanding_box_anisotropy*` | `engineering_proxy` until full box qualification closes | Particle and MHD transforms now exist and a separate uniform-MHD appendix invariant oracle passes, but the trend proxy is not a driven-CRPAI reproduction | Add CPAW and driven-CRPAI oracles |
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
  runs/<campaign>/<submission-id>/ # exact immutable simulation output and restarts
  metrics/<campaign>/           # compact machine-readable analyses
  figures/<campaign>/           # generated plots
  manifests/<campaign>/         # provenance manifests
  ledger/node_hours.jsonl       # authoritative append-only accounting ledger
  ledger/node_hours.csv         # derived RFC-4180 human-readable index
  ledger/mirror_receipts.jsonl  # immutable non-recursive mirror receipts
```

### Archive Retention And Restore Procedure

The user selected Orion as the sole bulk-evidence root and explicitly removed
Kronos from the execution design. Orion-only retention is a documented
durability risk, not an institutional or approved off-site archive. Before
closing a qualification gate:

1. Build a sign-off bundle containing raw outputs required for reproduction,
   compact metrics, figures, inputs, analysis scripts, logs, scheduler
   accounting, executable checksum, CMake cache, selected environment profile,
   ledger snapshot and claim-register entry.
2. Generate a SHA-256 manifest for the bundle and verify it before transfer.
3. Store required evidence below the authorized Orion root, record bundle
   destination, timestamp, checksum verification and any access restriction,
   and keep the Orion-only durability-risk disposition visible for review.
4. Keep compact source-controlled metadata pointing to the immutable Orion
   bundle location. Do not describe Orion-only retention as durable archival.
5. Restore a representative bundle into a fresh directory and regenerate at
   least one analytical metric, one stochastic summary and one paper figure
   before final sign-off.
6. Recheck the selected Orion root and current retention policy before every
   archive operation.

For the node-hour ledger, durability is continuous rather than a final-gate
action. The authoritative ledger is append-only JSONL; `node_hours.csv` is a
derived human-readable index. Serialize each primary event as canonical JSON
with UTF-8 encoding, lexicographically sorted keys, no insignificant whitespace
and a trailing newline. Compute `event_sha256` over that serialization while
omitting only the `event_sha256` field itself. After every reservation, job-ID
attachment, cancellation and accounting reconciliation, append a hash-chained
primary event and mirror the event plus sequence head to the preflighted Project
Home ledger mirror using mounted `filesystem_copy`. Record the exact target
before the first reservation. Each primary event contains `sequence_number`,
`previous_event_sha256`, `event_sha256`, `event_type`, identity, accounting and
artifact fields. After transport succeeds, append a separate immutable
`mirror_ack` receipt to `mirror_receipts.jsonl`. A receipt refers to the
mirrored primary event's SHA-256 and contains `mirror_destination`,
`mirror_transport`, `mirror_acknowledged_utc` and `mirror_ack_sha256`.
Serialize receipts with the same canonical-JSON rule and compute
`mirror_ack_sha256` while omitting only that field. Receipts do not participate
in the primary hash chain and are not themselves mirrored, avoiding recursive
acknowledgement. The Project Home mirror's primary-chain head and the Orion
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
- The `batch` partition is the default partition for production work and the
  default QOS is `normal`.
- The `debug` QOS permits only one user job in any state and rejects walltimes
  exceeding 2 hours. It is intended for short non-production debugging, not
  production work or job chaining.
- Orion project-work storage is not backed up and is purge-eligible; required
  evidence must be mirrored to approved nearline storage and exported to an
  institutional or approved off-site long-term archive.
- The default Frontier GPU mode is `HSA_XNACK=0`; setting `HSA_XNACK=1` changes
  page-migration behavior and must be treated as an experimental profile until
  benchmark evidence supports promotion.

The user selected Orion-only bulk-evidence retention despite the durability
warning above. This removes the storage block on registered Frontier execution,
but it leaves terminal durable-retention disposition open for external review.

This project uses a `debug`-preferred, `normal`-fallback scheduling policy on the
`batch` partition. Because OLCF identifies `debug` as short non-production use
and prohibits production work and job chaining under that QOS, do not attempt to
evade the policy by manually chaining restart segments. Prefer `debug` only for
eligible short qualification and debugging runs when the user has no `debug`
job in any state. Use `normal` on `batch` when the `debug` slot is occupied or
when the job is not eligible for `debug`, including registered longer
paper-reproduction, nonlinear-saturation and controlled-scaling work after its
prerequisite gates close. Keep AthenaK PIC submissions serial under both QOS
classes.

### Authorized Frontier Execution Scope And QOS Selection

Eligible short non-production jobs should request `#SBATCH -q debug` when the
user has no `debug` job in any state. These jobs may:

- verify HIP/MPI compilation, GPU placement and startup metadata;
- run analytical pusher, coupling, delta-f and expanding-box oracles;
- run short AMR, migration, load-balance, boundary and restart stress tests;
- run short Bell, oscillation, CRSI and CRPAI qualification windows;
- run bounded shock setup, output, checkpoint and single-resolution pilots;
- measure small controlled performance A/B comparisons needed to choose a
  recorded environment setting.

No `debug` job may:

- chain restart segments to approximate a production allocation;
- run a full shock reproduction that requires production-scale time or repeated
  restart chaining;
- run a scaling study whose purpose is production performance characterization;
- exceed the one-job or two-hour `debug` constraints;
- be used for production work, production scaling or restart chaining;
- cause the cumulative reserved-plus-consumed usage to exceed the 10000
  node-hour project testing cap.

Use `#SBATCH -q normal` on the `batch` partition when a short eligible job cannot
use `debug` because the user's `debug` slot is occupied, or when a registered job
is not eligible for `debug`. Record exactly one QOS-selection reason in the
pre-submit manifest: `debug_available`, `debug_slot_occupied`,
`debug_ineligible_production`, `debug_ineligible_walltime`, or
`normal_required_by_registered_campaign`. A `normal` submission is not a waiver:
all provenance, prerequisite, accounting, artifact and serial-submission checks
remain mandatory. Recheck the current OLCF scheduling policy before every
substantial campaign and record `site_policy_checked_utc` in the pre-submit
manifest; never assume that the site's `normal`-QOS limits are static.

### Authorization Boundary For Any Expanded Campaign

Registered paper reproduction, nonlinear saturation and controlled scaling may
run under `normal` on `batch` after their prerequisite gates close. Before any
submission outside the scheduling, account or budget envelope in this plan,
future agents must:

1. Stop and obtain explicit user and site-policy authorization for the expanded
   execution policy. Always stop and ask the user before raising the cumulative
   cap above 10000 node-hours.
2. Amend this plan with the approved QOS, partition, walltime, node-count
   ceiling, concurrency rule and revised budget envelope.
3. Preserve the same `/lustre/orion/ast207/proj-shared/dfielding/PIC` execution
   root, provenance manifests and node-hour ledger unless explicitly revised.
4. Freeze the qualifying commit, binary checksum, deck matrix, statistical
   design, predeclared metrics and stop conditions before submission.
5. Run registered large jobs only after all prerequisite small correctness gates
   are closed.
6. Record `debug` and `normal` artifacts distinctly so no short pilot is
   mistaken for a production result.

### Mandatory Resource And Accounting Policy

1. Set `PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC` for every
   Frontier action. Run every simulation under that exact root. Do not launch,
   execute, continue or write simulation output in another project directory.
2. Use the `batch` partition. Prefer `#SBATCH -q debug` only for eligible short
   non-production work when the user has no `debug` job in any state; otherwise
   use `#SBATCH -q normal`. A `debug` request must not exceed `02:00:00`.
3. Submit at most one AthenaK PIC job at a time across both QOS classes. Direct
   `sbatch` use is prohibited. Before each validated-wrapper submission, run:

   ```bash
   squeue -u "$USER" -h -o "%i %P %q %T %j %k"
   ```

   If any user `debug` job is listed in any state, select `normal` rather than
   submitting another `debug` job. If any AthenaK PIC job is listed in any state,
   do not submit another PIC job under either QOS.
4. Calculate the maximum additional budget before submission:

   ```text
   maximum_node_hours = requested_nodes * requested_walltime_hours
   ```

   Refuse submission if:

   ```text
   cumulative_consumed_node_hours
     + currently_reserved_node_hours
     + maximum_node_hours > 10000
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
   submission_id,job_id,git_commit,campaign,test_id,partition,qos,qos_selection_reason,
   queue_snapshot_sha256,site_policy_checked_utc,requested_nodes,
   scheduler_reported_allocated_nodes,billed_nodes,
   requested_walltime,reserved_node_hours,elapsed_seconds,consumed_node_hours,
   cumulative_consumed_node_hours,state,reconciled,artifact_dir,
   mirror_destination,mirror_transport,mirror_acknowledged_utc,mirror_ack_sha256,
   notes
   ```

7. Use the budget deliberately. Begin with one-node correctness tests and only
   grow when a previous result satisfies its gate. A failed setup should cost
   no more than one small `debug` job before it is corrected. Stop and ask the
   user for permission before increasing the cap above 10000 node-hours.

### Frontier Environment Profiles

Store the reusable environment setup only inside the installed checksummed
control plane as
`$PIC_ROOT/control_plane/<reviewed-digest>/frontier_pic_environment.sh` and
source it with a checked return status from build/run scripts. Verify module
versions when the OLCF software stack changes. Start from the minimum supported
profile. Treat page migration and communication tuning as experiments until
matched AthenaK A/B evidence supports promotion.

```bash
#!/bin/bash
PIC_FRONTIER_PROFILE="${PIC_FRONTIER_PROFILE:-frontier_minimum_supported}"
export PIC_FRONTIER_PROFILE

case "$PIC_FRONTIER_PROFILE" in
  frontier_minimum_supported|frontier_xnack1_experimental|frontier_ofi_tuned_experimental)
    ;;
  *)
    printf 'Unsupported PIC_FRONTIER_PROFILE=%s\n' "$PIC_FRONTIER_PROFILE" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

if ! module reset \
    || ! module load PrgEnv-amd/8.6.0 \
    || ! module load amd/6.2.4 \
    || ! module load rocm/6.2.4 \
    || ! module load craype/2.7.33 \
    || ! module load cray-mpich/8.1.31 \
    || ! module load cray-pmi/6.1.15 \
    || ! module load cray-libsci/24.11.0 \
    || ! module load craype-accel-amd-gfx90a; then
  printf 'Failed to load Frontier PIC module profile\n' >&2
  return 1 2>/dev/null || exit 1
fi

if module is-loaded darshan-runtime && ! module unload darshan-runtime; then
  printf 'Failed to unload inactive darshan-runtime module\n' >&2
  return 1 2>/dev/null || exit 1
fi

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export SLURM_EXPORT_ENV=ALL
export ROCM_PATH=/opt/rocm-6.2.4

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
esac

record_pic_environment() {
  local name value
  printf 'PIC_FRONTIER_PROFILE=%s\n' "$PIC_FRONTIER_PROFILE"
  printf 'HSA_XNACK=%s\n' "${HSA_XNACK:-0}"
  for name in MPICH_ENV_DISPLAY MPICH_VERSION_DISPLAY \
      MPICH_GPU_SUPPORT_ENABLED MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED \
      MPICH_OFI_NIC_POLICY MPICH_GPU_IPC_CACHE_MAX_SIZE MPICH_MPIIO_HINTS \
      MPICH_OFI_NUM_CQ_ENTRIES FI_MR_CACHE_MONITOR FI_CXI_RX_MATCH_MODE \
      OMP_NUM_THREADS SLURM_EXPORT_ENV ROCM_PATH; do
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
| `frontier_selected_production` | Selected large-run profile | Define only after matched environment-profile evidence and before any registered large `normal`-QOS campaign |

Entity documentation reports application-specific Frontier communication
experience. Treat that information as a prompt for a bounded AthenaK benchmark,
not as permission to change the minimum environment above. After correctness
gates pass, run matched short A/B cases with the OLCF-supported GPU-aware MPI
baseline and any proposed alternative, checking physical metrics, runtime and
memory before selecting a recorded production setting. Archive the redacted
allowlist emitted by `record_pic_environment`; do not archive unrestricted
`env | sort` output because environment variables can contain sensitive values.

### Installed Frontier Build Profile Workflow

The example supplied with the original validation request builds an unrelated
problem generator and leaves configuration provenance incomplete. MHD-PIC
qualification must build the built-in PIC problem generators, store the cache
and log, and make the resulting executable immutable by commit/config identity.
The installed checksummed control plane is the build authority. Do not create a
campaign build script that reassembles CMake flags, output paths, provenance
sidecars, or profile JSON. Source the installed environment profile with a
checked status and invoke only its build-profile writer:

```bash
set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
SRC_DIR=/ccs/home/dfielding/athenak-pic
CONTROL_PLANE_VERSION=<reviewed-control-plane-digest>
CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}"
ENV_FILE="${CONTROL_PLANE_DIR}/frontier_pic_environment.sh"
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
CONTROL_PLANE=("$PYTHON" -I "${CONTROL_PLANE_DIR}/run_control_plane.py")
source "$ENV_FILE" || exit $?

GIT_COMMIT_FULL="$(git -C "$SRC_DIR" rev-parse HEAD)"

"${CONTROL_PLANE[@]}" write_orion_build_profile.py \
  --source-root "$SRC_DIR" \
  --expected-git-commit "$GIT_COMMIT_FULL" \
  --profile-id hip-mpi-release-paper-pic
```

The writer accepts only `--source-root`, `--expected-git-commit`, and
`--profile-id hip-mpi-release-paper-pic`. It verifies the installed control
plane and authorized clean source closure, derives every Orion artifact and log
path, requires fresh build, bin, and log paths, and materializes a fresh local
detached checkout plus recursive local submodule checkouts under the derived
build directory. It then invokes a closed direct `/usr/bin/cmake` configure and
build argv: Release, explicit double precision through
`Athena_SINGLE_PRECISION=OFF`, MPI, HIP, `Kokkos_ARCH_ZEN3`,
`Kokkos_ARCH_AMD_GFX90A`, the Cray
`/opt/cray/pe/craype/2.7.33/bin/CC` wrapper, `/opt/rocm-6.2.4` include and HIP
link flags, built-in problem generators, and parallel build width 32. Future
toolchain changes belong in the installed writer and its review, not in an
operator-authored campaign recipe.

The writer captures the exact configure and build argv in
`build-invocations.json`; empty `git_status.preconfigure.txt` and post-build
`git_status.txt` captures from the fresh checkout; configure/build logs;
`CMakeCache.txt`; `modules.txt`; the fixed `toolchain.txt`; recursive
`submodule_status.txt`; the redacted `environment.allowlist.txt`; and the
retained minimal `build-environment.json`. It binds those eleven provenance
inputs, the clean source and recursive-submodule closure, and the executable
digest in a schema-v3 `build_profile.json`, then publishes the adjacent
`profile_receipt.json` exclusively under Orion.

Run compiles on an appropriate Frontier login/build workflow according to OLCF
policy; run simulations only through allocated compute resources. The closed
argv, hashes, and fresh detached checkout reduce stale or mixed-build mistakes;
they do not cryptographically prove that the compiler honored the recorded
command or establish compiler semantics. Reproducibility and qualification
remain separate gates. Descriptor-relative publication below pinned authorized
parents closes ancestor-swap pathname redirection during publication; it does
not protect writable artifacts against a malicious same-UID process.

Create one candidate freeze from the helper-published profile:

```bash
GIT_COMMIT="${GIT_COMMIT_FULL:0:12}"
CONFIG=hip-mpi-release-paper-pic
BIN_DIR="${PIC_ROOT}/bin/${GIT_COMMIT}/${CONFIG}"

"${CONTROL_PLANE[@]}" create_clean_candidate_freeze.py \
  --source-root "$SRC_DIR" \
  --executable "${BIN_DIR}/athena" \
  --build-profile "${BIN_DIR}/build_profile.json" \
  --build-profile-id "$CONFIG" \
  --prepared-artifact-inventory \
    tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json
```

Regenerate and review the committed prepared-artifact inventory before freezing
the clean candidate. Freeze creation reads that source-relative JSON from the
generated `source.tar` and revalidates every listed PIC deck and publication
analyzer SHA-256 from archived bytes.

### Improved Frontier Job Template

Create one immutable job script per test/campaign under `$PIC_ROOT/jobs/`.
This template uses one GPU per MPI rank as recommended by the Frontier guide:

Before submission, run `create_pre_submit_manifest.py` through
`${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}/run_control_plane.py`
from the installed checksummed control-plane snapshot. Keep each campaign's
immutable job script, queue snapshot, timeout artifact and pre-submit
configuration under `$PIC_ROOT/jobs/<campaign>/`. The manifest creator must
create an immutable `pre_submit_manifest.json` under
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
append-only status stream; record the manifest digest in the reservation chain,
export that exact digest into the scheduled job and reject a compute-node
manifest mismatch. Never mutate the pre-submit manifest from inside an
allocation. Direct `sbatch` use is prohibited: submit only through the validated
wrapper below.

The manifest creator must resolve `REPLACE_WITH_VALIDATED_QOS` before snapshot:
choose `debug` only when the registered job is eligible for short
non-production use and `squeue` shows no user `debug` job in any state;
otherwise choose `normal`. Snapshot the queue observation, selected QOS and
QOS-selection reason. The validator must recompute that decision immediately
before reservation and fail closed if a race or policy change makes the
snapshotted selection stale.

Freeze the submission wrapper, validator, reconciler, ledger initializer and
schema under an immutable checksummed control-plane version directory. Invoke
that frozen wrapper path, record every control-plane checksum and schema version
in each pre-submit manifest and ledger event, and reject checksum drift before
reservation. Production installation in Orion and Project Home is permitted
only from reviewed clean tracked control-plane source files. Installation,
clean-candidate freeze and manifest-snapshot creation must retain a pinned
authorized-parent descriptor while staging, publish with descriptor-relative
rename, sync the pinned parent and fail closed if the lexical parent no longer
names that directory. This closes ancestor-swap pathname redirection during
publication, but it is not an integrity boundary against a malicious process
running as the same Unix UID. Do not execute mutable scripts directly from
`$PIC_ROOT/jobs/`.

The shell block below is a legacy directive-and-argv reference only. The
successor control plane must translate its Athena invocation into a closed
structured `launch_contract`; it must not execute this mutable shell body.
The installed trampoline captures each
`<action_id>.environment.allowlist.txt` through an anchored inherited file
descriptor immediately before `srun`. Before workload execution, the wrapper
syncs that descriptor, changes it to mode `0400`, syncs it again, syncs its
pinned directory descriptor, closes both descriptors and unsets their
environment bindings.

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
#SBATCH -q REPLACE_WITH_VALIDATED_QOS
#SBATCH -N 1

set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PATH=/usr/bin:/bin
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
ENV_FILE="${PIC_ROOT}/jobs/frontier_pic_environment.sh"
GIT_COMMIT=REPLACE_WITH_COMMIT
CONFIG=hip-mpi-release-paper-pic
CAMPAIGN=paper_pusher_validation
TEST_ID=gyro_section51_p2_1n
SUBMISSION_ID="${PIC_SUBMISSION_ID:?validated wrapper must export PIC_SUBMISSION_ID}"
RESERVATION_ID="${PIC_RESERVATION_ID:?validated wrapper must export PIC_RESERVATION_ID}"
MANIFEST_SHA256="${PIC_MANIFEST_SHA256:?validated wrapper must export PIC_MANIFEST_SHA256}"
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
OUTDIR="${PIC_ROOT}/runs/${CAMPAIGN}/${SUBMISSION_ID}"
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
"$PYTHON" -E -s "$SNAPSHOT_VERIFIER" \
  --manifest "$PRE_SUBMIT_MANIFEST" \
  --submission-id "$SUBMISSION_ID" \
  --reservation-id "$RESERVATION_ID" \
  --manifest-sha256 "$MANIFEST_SHA256" || preflight_failed snapshot_checksum_mismatch
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

Submission is an explicit serial AthenaK PIC workflow across both allowed QOS
classes. Install and test
`${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}/run_control_plane.py validate_and_reserve_frontier_job.py`
and
`${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}/run_control_plane.py reconcile_frontier_job.py`
before the first Frontier submission. Execute reservation, submission and
reconciliation only through that installed checksummed snapshot; keep immutable
campaign job scripts and submission metadata under `$PIC_ROOT/jobs/<campaign>/`.
The validator must parse the immutable job script and pre-submit manifest, then reject
submission unless all policy checks pass:

- job script requests exactly `#SBATCH -p batch` and exactly one allowed QOS,
  `#SBATCH -q debug` or `#SBATCH -q normal`;
- the pre-submit manifest records the selected QOS, queue snapshot, one allowed
  QOS-selection reason and the timestamp of the current OLCF scheduling-policy
  check;
- a `debug` request is a short non-production task, requests walltime
  `<=02:00:00`, and `squeue` shows no existing user `debug` job in any state;
- a `normal` request is selected when `debug` is occupied or ineligible and its
  registered walltime, node count, campaign bounds and current site-policy
  limits are satisfied;
- requested nodes and worst-case node-hours are parsed successfully;
- consumed plus currently reserved plus requested node-hours is `<=10000`;
- `squeue` shows no existing AthenaK PIC job for the user in any state across
  either allowed QOS, identified by the `pic-reservation=` Slurm comment and
  snapshotted job metadata;
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

The installed submission wrapper is the only operational implementation; do
not fork it into campaign shell. Its fail-closed transaction is:

```text
reserve and durably mirror worst-case node-hours
persist reserved_not_submitted marker
durably persist scheduler_dispatch_started marker
sbatch --parsable --hold ... installed launch_trampoline.py
durably persist scheduler_job_id_received marker
verify the held job with scontrol
durably persist submitted_not_attached marker
append and durably mirror scheduler-ID attachment
scontrol release <job-id>
```

Every registered run artifact directory is exactly
`$PIC_ROOT/runs/<campaign>/<submission-id>`. Reservation rejects an existing
artifact directory and a launch creates that exact directory exclusively.
Reruns use a new submission ID.

Before the first submission, initialize an explicit genesis primary event,
mirror it, append its receipt and derive the CSV index through the reviewed
`initialize_frontier_ledger.py` workflow. Never create an empty ledger implicitly.

After completion, run the reconciler to archive `sacct`, replace reserved usage
with actual usage, and record job state. A human or agent must still inspect the
run identity, artifact completeness, failure state and cumulative total before
the next submission. The validator must refuse a new reservation while any prior
row remains unreconciled, while the durable ledger head is stale, or while
`pending_submission.json` exists without reviewed attachment, cancellation and
reconciliation. The durable marker states are `reserved_not_submitted`,
`scheduler_dispatch_started`, `scheduler_job_id_received` and
`submitted_not_attached`. Only failures known to precede
`scheduler_dispatch_started` may automatically cancel an unused reservation.
Recovery must query `squeue` and `sacct` for the reservation token stored in the
Slurm comment, attach or cancel the matching job, reconcile terminal usage even
if the job reached a terminal state before attachment, and clear the pending
marker only after the durable mirror acknowledges the transition.
Every reservation, attachment, cancellation and reconciliation mutation must be
durably mirrored before the wrapper or reconciler reports success. If
`scancel`, `sacct` reconciliation or durable mirroring cannot be confirmed,
block later submissions and require manual recovery.

### Authorized Frontier Testing Budget Envelope

The budget below is deliberately conservative. It uses maximum requested
node-hours, not expected elapsed cost. Do not spend a later tier while its
prerequisite phase gate is open.

| Tier | Prerequisite | Suggested job pattern | Maximum requested node-hours | Cumulative ceiling after tier |
| --- | --- | --- | ---: | ---: |
| F0 Build/run smoke and GPU mapping | Clean baseline plus parser tests | Up to 6 x 1 node x 0.25 h | 1.5 | 1.5 |
| F1 Relativistic pusher and conservative coupling oracles | P2/P3 implementation complete | Up to 16 x 1 node x 0.5 h | 8 | 9.5 |
| F2 Delta-f and box analytical tests | P4/P5 implementation complete | Up to 16 x 1 node x 1 h | 16 | 25.5 |
| F3 AMR/restart/decomposition/load balance validation | P6 implementation complete | Up to 20 jobs varying 1-4 nodes x <=1 h, budget-capped | 50 | 75.5 |
| F4 Bell, oscillation, CRSI, CRPAI short qualification cases | P3-P6 analytical gates closed | Up to 24 compliant `debug`-preferred, `normal`-fallback jobs varying 1-8 nodes x <=2 h, budget-capped | 200 | 275.5 |
| F5 Shock pilots and single-job resolution checks | P7 setup review complete | Up to 12 compliant `debug`-preferred, `normal`-fallback jobs varying 1-8 nodes x <=2 h, budget-capped | 120 | 395.5 |
| F6 Full shock reproduction, long nonlinear saturation and controlled scaling | All correctness gates closed | Registered `normal`-QOS jobs on `batch`; freeze each campaign matrix and statistical design before use | 2400 | 2795.5 |
| F7 Selected extension qualification | Production-ready `paper_mhd_pic` child gates and extension-specific host oracles closed | Registered Hall Bell, ion-neutral-damped CRSI and adaptive-delta-f physical-damping CRPAI comparison campaigns under `normal` on `batch` where required | 3000 | 5795.5 |
| Reserve | Unexpected defects, justified reruns and registered follow-up tests only | Must be justified in ledger/plan before use | 4204.5 | 10000 |

This envelope is not permission to spend the allotment. Stop as soon as enough
evidence exists to decide a gate. If a one-node test fails, fix it before
allocating more nodes.

The reserve is not permission to skip prerequisite gates or broaden campaign
scope silently. Any future expanded envelope must be added as a new revision
with its own authorization record and cumulative accounting. Stop and ask the
user for permission before increasing the 10000 node-hour cap.

## Detailed Qualification Matrix

Future agents should expand the table with commit hashes, artifact locations
and measured values as work progresses.

| Gate | Capability | Required test/evidence | Status at initial review |
| --- | --- | --- | --- |
| Q-001 | Clean provenance | Curated baseline, manifest, executable checksum | Verifying: historical `origin/PIC` baseline remains pinned at `3bcd3f21`; exact clean source commit `4cceb5d4`, tree `933f2e3a`, pinned `kokkos`, source bundle `e11dc8fb`, immutable clean-candidate manifest `3174515d` and executable `bea2a418` are archived under freeze `31be2cd6`; terminal release aggregate remains open |
| Q-002 | Physical-mode configuration | Parser rejection and metadata tests | Verifying: host parser contract accepts five bounded named-mode compositions, rejects 28 invalid combinations and checks runtime identity; registered F0 Frontier metadata capture and parser/startup smoke pass. Bounded clean-candidate F2 eight-rank one-node Athena runtime-metadata v1 job `4746310` completed and passed its inspected-artifact analyzer, but remains rejected before evidence acceptance because its launch contract omitted `-d` and Athena error-history output escaped the registered run tree. Repaired v2 from `d59b8ebf` routes that output below the artifact tree and requires it during descriptor-pinned offline recomputation. Attested v2 job `4746316` completed on one Frontier node in 9 seconds, reconciled `0.0025` node-hours and passed descriptor-pinned no-write replay with inventory SHA `ca673515` and result SHA `95dc6e40`; preserve that bounded engineering evidence while external review and broader F2 coverage remain open |
| Q-003 | Relativistic pusher and `C` | Section 5.1 analytic convergence CPU/GPU | Verifying: host continuum gyro scan passes three `C` values and four timesteps with approximately second-order convergence; registered clean-candidate Frontier F1 gyro-v3 one-rank GPU oracle passes at cycle 2 with 64 particles and maximum Boris velocity error `5.808176634092277e-09`, including visible-GPU and HIP-link evidence; the full GPU matrix remains open |
| Q-004 | Paper coupling | Conservation tests and no direct CT-current induction | Verifying: host paper-mode task-stage trace, particle-plus-fluid conservation and coefficient-invariant ideal-MHD induction isolation pass; registered clean-candidate Frontier F1 paper-coupling-v2 one-rank GPU conservation and coefficient-invariance oracle passes; the multi-rank MPI/GPU and paper-campaign matrix remain open |
| Q-005 | Bell | Section 5.2 dispersion and phase | Verifying: bounded serial-host extractor validation passes and the reused Bell fixture is now quantitatively classified as an engineering proxy because its growth and phase miss the deck-bound no-Hall analytical branch; source-local `Q023-PAPER-BELL-LINEAR` preparations now freeze Section 5.2 paper-value 1D/2D/3D geometry and normalization, a dedicated right-polarized eigenmode generator and a fail-closed combined magnetic-plus-fluid-velocity raw-trace analytical-dispersion candidate requiring signed phase, approved materialized geometry and exact raw-artifact provenance. The deterministic source-local materializer has staged and recursively frozen the fixed 405-deck epsilon, dimensionality, resolution, local-CFL and PPC matrix only as an ignored descriptor-relative direct child below `tst/.codex`, with exact inventory recorded in `q023_paper_bell_linear_materialized_variants_local_2026-05-30.json`, while retaining centered loading as an explicit review boundary. The physical 1D preparation uses an exact transverse-invariant thin 2D3V carrier because the particle module intentionally supports 2D/3D meshes. Bounded Debug `mhd_w_bcc` cycle-zero initialization passes for all three dimensions and a 2D schema-7 restart continuation advances one cycle. The retained cycle-zero projected fluid-velocity-to-magnetic ratio is a nonqualifying source-local initialization diagnostic only; the source-local analyzer separately freezes the manuscript-literal spatial sine-fit `delta u_y` phase and volume-averaged `|delta u|` growth contract. Clean freeze, qualifying variants, convergence, MPI, GPU and external review remain open |
| Q-006 | Oscillation plus grid refinement | Section 5.3 uniform/SMR/AMR frequency | Verifying: bounded serial-host extractor validation passes, reused uniform/SMR proxies have broad frequency parity and the nested level-2 AMR fixture remains finite short-horizon interface smoke. `q006_paper_multispecies_oscillation_source_local_preparation_2026-05-30.json` adds a separately named registered guarded source-local generator and exact source-local decks. `q006_paper_multispecies_oscillation_runtime_local_2026-05-30.json` adds the narrowly guarded exact-isothermal full-f momentum-only runtime successor and recursively read-only Orion uniform startup, uniform evolution, SMR evolution, deterministic audited-AMR evolution and audited-AMR restart snapshot recomputations. Accepted long-horizon uniform/SMR/true-AMR residual tolerances, true-AMR policy qualification, MPI, GPU, Frontier and external review remain open |
| Q-007 | Delta-f | Weight/deposition tests, CRSI/CRPAI reproduction | Verifying: evolving state/deposition/restart/VTK mechanics smoke passes; `q007_paper_deltaf_linear_source_local_preparation_2026-05-30.json` retains the narrowly guarded exact-isothermal `paper_mhd_pic` true-delta-f parser allowance and separately named CRSI/CRPAI thin-2D3V carriers. The source-local successor implements eight logarithmic momentum bins with IPWT macro-multiplicity weighting, a deterministic antipodal angular sampler, a deterministic four-branch wave carrier, raw startup MHD-cell and particle-loading validation for CRSI and both CRPAI carriers, static full-`Q1 + Q2` dispersion-oracle preparation and a recursively read-only Orion two-cycle CRSI mechanics replay from a pinned executable. The finite shell quadrature, angular sampler, wave seed and discrete mode set remain explicit source-local conventions rather than recovered paper-run provenance. Long-horizon growth fitting, paper-run provenance review, CRPAI handedness-label review and runtime replay, MPI, GPU, Frontier and external review remain open |
| Q-008 | Expanding box | CPAW, gyro and driven CRPAI tests | Verifying: particle mechanics, uniform MHD appendix invariants, bounded gyro-history convergence, one-step expanding/compressing CPAW closure, comoving-flux physical-view repair, oblique native-divB, built-in physical-volume MHD history, non-delta-f final-frame conservative source accounting and selected adaptive-delta-f endpoint-source plus reduced-damping ordering pass locally; `q008_expanding_box_restart_resilience_bounded_local_2026-05-30.json` adds a passing uninterrupted-versus-restarted coupled-carrier endpoint regression; `q008_expanding_box_cpaw_history_convergence_preparation_2026-05-30.json` records the retained Orion-backed ten-case serial exponential-profile Appendix-A CPAW four-crossing-time preparation matrix plus a recursively read-only 22-case extension with literal linear and reciprocal-linear convergence preparations, LLF/HLLE/HLLD sensitivity, x1/x2/x3 axis-aligned carriers and a fail-closed unsupported Roe probe. Publication thresholds, expanded dimensional convergence, solver-sensitivity qualification, driven CRPAI, MPI, GPU, Frontier and external review remain open |
| Q-009 | AMR lifetime and boundary policy | Forced rebuild/restart/interface tests | Verifying: explicit particle refresh, retained-state inventory, `paper_smooth` interface-policy decision and serial refine smoke pass are archived in `q009_amr_lifetime_policy_successor_2026-05-30.json`; `q009_repeated_amr_lifetime_bounded_local_2026-05-30.json` adds six serial refine/derefine transitions through five restart continuations with stable particle identity and ownership refresh; `q009_coupled_boundary_lifetime_bounded_local_2026-05-30.json` adds serial Debug and UBSan coupled PIC/MHD reflecting/outflow lifetime stress; preserved blocker `q009_coupled_inflow_lifetime_bounded_local_2026-05-30.json` exposed an uninitialized MHD inflow reservoir, restart `ProblemGenerator` leak and `ParticleMeshBlockOffset` overflow path; repaired successor `q009_coupled_inflow_lifetime_repaired_bounded_local_2026-05-30.json` records passing Debug and host ASan/UBSan six-transition x2-inflow replay after reservoir initialization, stack scratch storage and guarded offset conversion; multi-rank MPI migration/restart, HIP and scientific-AMR qualification remain open |
| Q-010 | Load balance/performance | Particle-cost implementation and Frontier measurements | Verifying: opt-in particle-aware AMR cost implementation landed; Frontier measurements open |
| Q-011 | Shock | Section 5.4 coarse/fine/AMR reproduction | Blocked by Q-003/Q-004/Q-009: source-local manuscript ideal-surface, unique half-open surface carrier, upstream-relative swept-mass budget, full-sphere isotropic monoenergetic injection, restart-persisted runtime removal of the startup cohort, frame controls and schema-6 continuation fingerprints pass bounded guards while the finite-Mach surface remains an explicitly separate engineering option; schema-7 initial-versus-shock-injected provenance, restart preservation, MeshBlock migration and independently reconstructed weighted spectra pass the bounded Q-016 successor; `q011_parallel_shock_section54_paper_preparation_2026-05-30.json` freezes the Section 5.4 preparation-only 2D3V deck, three-level AMR sizing, `M_A=30`, `C/U_A0=10000` and continuous injection with runtime removal of particles born before `45 Omega0^-1`, and records a recursively read-only Orion three-stage restart-before, crossing and restart-after cutoff diagnostic. `q011_injection_distribution_runtime_local_2026-05-30.json` adds an executed recursively read-only bounded serial-host audit of actual PVTK shock-injected provenance, clamped shock-surface placement, monoenergetic surface-relative speed and full-sphere statistics. Physical calibration, AMR, MPI, GPU and campaign evidence remain open |
| Q-012 | Reliability aggregate report | Closes after Q-024/Q-036 plus MPI/decomposition/restart-continuation/boundary evidence | Partial scaffolding only |
| Q-013 | Usability/docs | Supported-mode docs and runbook tested by clean launch | Verifying: clean-launch runbook, supported-toolchain declaration, runtime model contract, AMR lifetime/interface policy and bounded Q-016 provenance/spectrum page landed; isolated `docs/requirements.txt` environment renders all 58 Sphinx pages with warnings promoted to errors; controlled clean-candidate Frontier F0 launch and registered F1 gyro-v3 plus paper-coupling-v2 slices pass, while the broader production-launch matrix remains open |
| Q-014 | Terminal sign-off aggregate | After every selected-profile gate closes, freeze a review-ready aggregate manifest with `pending_terminal_review`; then review, record the terminal disposition, sign and archive the final manifest. Never list Q-014 as its own child | Open |
| Q-015 | Entity-derived shared-kernel comparison | Frozen source/formula provenance, relativistic pusher microtest and bounded trajectory-deposition differential checks | Verifying: frozen Entity Boris and particle-shape exact-overlap differential passes with explicit non-overlap boundaries; full current-kernel equivalence, MPI and GPU remain open |
| Q-016 | Particle provenance and spectra | Persistent tracking/cohort metadata plus in-run/offline spectral agreement through restart/migration | Verifying: schema-7 persistent `gid`, tag, species, source, macro-weight, birth-time and delta-f fields; strict typed PVTK decode; tracked-row metadata; and independently reconstructed weighted species/source/birth-cohort spectra pass a bounded serial restart and MeshBlock-migration regression archived in `q016_particle_provenance_spectra_local_2026-05-30.json`. The MPI-ready harness merges schema-compatible same-cycle optional `gid` slices before comparison and rejects missing or mismatched `POINT_DATA`, duplicate sections, non-float vectors, trailing unknown content and duplicate tags; the synthetic `gid`-slice fixture is not MPI evidence. Controlled two-rank, one-node Orion-shared MPI decomposition replay passes in `q016_particle_provenance_spectra_mpi2_orion_local_2026-05-30.json`; registered clean-candidate Frontier HIP parity, repeated AMR/HIP lifetime and preregistered shock-campaign evidence remain open |
| Q-017 | Performance observability and Frontier communication choice | Stage timers, memory/load telemetry and recorded GPU-aware MPI A/B decision | Verifying: bounded driver-level task-list, output-publication, AMR/load-balance, throughput and safely derivable rank-load telemetry are archived in `q017_driver_performance_observability_2026-05-30.json`; `q017_particle_observability_local_2026-05-30.json` adds passing serial and one-rank MPI host regressions for wrapper-boundary particle timers, fixed-record resident bytes by species and absolute logical level, and direct particle-view allocated-byte snapshots; successor `q017_owned_helper_allocation_high_water_local_2026-05-30.json` adds final snapshots and observed strict rank-local high-water margins for the bounded AthenaK-owned particle, boundary and MeshRefinement Kokkos-view spans visible at owned allocation boundaries, including migration lists and counters sampled while live and explicitly excluding allocator overhead, runtime caching, untracked subsystems and any temporally concurrent multi-rank peak claim; Frontier synchronized-timer overhead, tracked GPU measurement, total GPU-memory accounting, GPU-aware-MPI A/B, multi-rank scaling and retained manifest evidence remain open |
| Q-018 | Claims registry | Stable claim IDs, evidence links, limitations and reviewer dispositions for every release statement | Verifying: stable IDs, a per-claim applicability, limitation and current-link scaffold, and the separate fail-closed `pic_qualification_manifest.py` validator/freezer are archived; registered-F1 gyro-v3 and paper-coupling-v2 terminal qualification manifests are refreshed through historical bounded-F2 projection policy SHA `446db26a` and linked as immutable pending-external-review evidence without promoting any claim, while strict operational policy SHA `36228012` remains live; the paper-Bell link now exposes its source-local paper-value preregistration profile and the paper-shock link exposes the fail-closed Q-009 inflow blocker record; remaining qualifying bundle links and named external reviewer dispositions remain open |
| Q-019 | Bell nonlinear saturation | Amplification, spectra, morphology, energy partition, saturation and sensitivity matrix | Blocked by Q-003/Q-004/Q-005 |
| Q-020 | CRSI nonlinear saturation | True delta-f spectra, distribution evolution, scattering/diffusion, saturation, sensitivity matrix, matched reduced nonlinear full-f controls and archived weight-validity envelope | Blocked by Q-007 |
| Q-021 | Driven CRPAI nonlinear saturation | Driven-box branch evolution, anisotropy, scattering, saturation, sensitivity matrix, effective-damping trends, matched reduced nonlinear full-f controls and archived weight-validity envelope | Blocked by Q-007/Q-008 |
| Q-022 | Independent-comparison preregistration | Frozen references, equation mappings, normalization, observables, tolerances, uncertainty methods and discrepancy-ledger schema | Verifying: fail-closed inventory, paper-source anchor, AthenaK paper-mode equation profile, bounded Entity overlap map, map schema, tolerance schema, observable families, tolerance policy and discrepancy schema are frozen; non-Entity public source references and GIZMO interpretation artifacts are privately archived in Orion with post-copy verification under manifest SHA `7028552e` in `q022_external_reference_private_ingest_2026-05-30.json`; successor dataset-provenance, per-comparison equation-map and tolerance-table sidecars fail closed with `pending external review` placeholders and no fabricated measurements, synthetic complete test-only records exercise each schema ready path, and the GIZMO route has a provisional documented-exclusion candidate; authorized Orion extraction, numeric thresholds, external review and comparison runs remain open |
| Q-023 | Statistical qualification | Preregistered estimators, seeds, scans, intervals, exclusions and independent metric recomputation | Verifying: `q023_local_preregistration_profiles_2026-05-30.json` freezes eight qualifying seeds, two excluded pilot seeds, fixed-sample Holm/BCa policy, locally available grids, defensible windows and thresholds, fail-closed exclusions, checksummed local analyzer/input candidates and the independent raw-artifact recomputation plan; source-local `Q023-PAPER-BELL-LINEAR` preparation adds checksummed Section 5.2 paper-value 1D/2D/3D candidates, the dedicated eigenmode generator, a strict 15-row combined magnetic-plus-fluid-velocity analytical-grid candidate requiring signed phase, approved materialized geometry and exact raw-artifact provenance, bounded three-dimension Debug `mhd_w_bcc` initialization smoke and a 2D schema-7 restart continuation. The deterministic source-local materializer has staged and recursively frozen the fixed 405-deck epsilon, dimensionality, resolution, local-CFL and PPC preparation matrix below `tst/.codex`, with exact inventory recorded in `q023_paper_bell_linear_materialized_variants_local_2026-05-30.json`; it refuses qualification and authorization claims, and centered loading remains an explicit review boundary. The physical 1D case is carried on an exact transverse-invariant thin 2D3V mesh. The projected fluid-velocity cross-check is explicitly nonqualifying, and the source-local analyzer separately freezes the manuscript-literal spatial sine-fit `delta u_y` phase and volume-averaged `|delta u|` growth contract. All 12 campaign drafts remain `locally_preregistered_run_blocked`; clean-candidate and Frontier bindings, qualifying seed and selected time-limit variants, paper-shock fit interval, theory- or excluded-pilot-derived nonlinear regimes, external-reference tolerances, exact plotting lock, reviewer-owned recomputation artifact and external review remain open |
| Q-024 | Resilience aggregate report | Forced walltime, interrupted output, corrupt checkpoint, restart schema and bounded storage-failure tests; closes after Q-036 | Verifying: schema-6 hardening parity plus the schema-7 Q-016 provenance successor, canonical manifest binding, wrong-path, checksum, payload, selector, adaptive-state, star-potential and shock continuation-control guards pass; bounded serial-host schema mutation, short-write, seek-failure, `/dev/full`, interrupted-writer residue and soft-wallclock continuation-parity drills pass; per-rank MPI, node-loss, scheduler-pretimeout and Frontier-filesystem matrix remain open |
| Q-025 | Portability aggregate report | Debug/Release, warnings, host memory checking, declared CPU/OpenMP/MPI scope and Frontier HIP/MPI matrix; closes after Q-038 | Verifying: host Debug/Release serial, Debug/Release MPI compile, Debug/Release OpenMP compile, OpenMP 1/2-thread byte-identity, warnings-enabled build, GNU ASan bounded run and UBSan checks pass; registered Frontier HIP/MPI build, clean-candidate one-rank parser/startup launch, gyro-v3 Release GPU runtime and paper-coupling-v2 Release GPU runtime pass after the device-capture portability fix; the rank-launched Athena MPI and broader HIP matrix remain open |
| Q-026 | Archive integrity, licensing and sensitive-data handling | Pre-ingestion and pre-export classification, scrubbed submission artifacts, checksums, redistribution basis, manifest schema and fresh-directory restore drill | Verifying: paper source, bounded Entity snapshot and private 20-page rendered PDF checksums are frozen; Orion private-bundle restore drill passes; export review and the external disposition on Orion-only durability risk remain open |
| Q-027 | Frontier QOS-selection and budget boundary | `debug`-preferred, `normal`-fallback selection is recorded and policy-compliant; serialized PIC submission ledger is complete; cumulative reserved-plus-consumed usage remains `<=10000` node-hours; expanded work remains blocked pending revised authorization and user permission | Verifying: mirrored ledger genesis event `849bf340` and predecessor digests through compute-snapshot successor `e8e47ead` remain immutable chronology. Registered F0 job `4745842` passed and reconciled. Retry predecessor `6002c80e` was paired-installed and promoted with policy SHA `ee3923fb` and promotion SHA `433d0afc`; bounded gyro attempt `4746123` completed but was scientifically rejected before result publication because the v1 analyzer required empty stderr while Frontier emitted reviewed Cray MPICH diagnostics. Historical successor `cbc6fb50` paired-installed and promoted with policy SHA `58f85503` and promotion SHA `05dfbf32`, but a late independent review reproduced generated-inventory close/reopen replacement and offline initial-loader identity substitutions before a v2 reservation or scheduler submission; one fail-closed pre-submit manifest remains as chronology. Superseded staged repairs `d55d064d`, `d8cc7d61` and `d538af91` closed those gaps and then exposed narrower post-check replacement windows. Historical repaired successor `4ccde8df` retains workload-payload descriptors through a closing namespace and byte sweep, rebinds directory identities afterward and rechecks the retained generated-inventory descriptor. Fresh exploit retest, curated commit `6fa926a9`, paired immutable install inventory SHA `ea45cdc5`, existing-genesis migration, reviewed policy SHA `3895221e` and promotion SHA `016cbb66` pass. One operator queue-format mismatch manifest failed closed before ledger intent. Attested gyro-v2 job `4746245` completed and reconciled at `0.0030555555555555557` node-hours, but analysis rejected before result publication after exposing deterministic cycle-three output closure and emitted particle-VTK-schema omissions in the staged parser fixture. Paired gyro-v3 and paper-v2 parser binding from commit `c878dd6c` promoted with policy SHA `0538da6d` and promotion SHA `f84fc7b7`. Attested gyro-v3 job `4746290` completed, reconciled at `0.0033333333333333335` node-hours, passed offline analysis and initially froze pending-review qualification manifest SHA `9f7a5d64`. The first paper-coupling manifest creation failed closed before publication or ledger intent when the submission scrubber matched ordinary analyzer identifiers. Scrub-safe commit `ae9fbb97` retained the unused coupling-v2 authorization ID, promoted mirrored policy SHA `559a4b1e` and promotion SHA `84065c8a`, and attested paper-coupling-v2 job `4746297` completed, reconciled at `0.006111111111111111` node-hours, passed offline analysis and froze pending-review qualification manifest SHA `6262332a`. After that active-policy replacement, gyro-v3 replayed through current-policy successor qualification manifest SHA `17265e02`; preserve initial manifest `9f7a5d64` as immutable chronology and use `17265e02` for terminal review. Rejected F2 v1 job `4746310` preserved its escaped-output defect as chronology. Repaired F2 v2 policy SHA `446db26a` and promotion SHA `b772da00` passed paired mirroring; fail-closed queue-snapshot retries remained pre-reservation, and attested v2 job `4746316` completed, reconciled `0.0025` node-hours and passed descriptor-pinned no-write replay with contained output `output/f2_multirank_runtime_metadata-errs.dat`. Strict registered-lifecycle successor `6f3458ca` then passed its retained 354-test matrix and two fresh independent adversarial reviews, paired-installed inventory SHA `10bba8e7`, promoted policy SHA `36228012` with promotion SHA `a2e04161`, and imported the reviewed seven-row Q-016 direct-`srun` accounting-only tranche while the queue was empty. All three ledger streams are coherent at 56 records with cumulative consumption `0.11055555555555553` node-hours, no active reservation and no PIC pending marker. Registered F1 closure and bounded F2 engineering evidence pass; imported Q-016 rows are budget accounting only; broader authorization and external review remain open |
| Q-028 | Independent Bell nonlinear comparisons | Non-Hall paper-mode campaign compared against Bai et al., Riquelme-Spitkovsky, Gargaté et al. and Zacharegkas et al. where regimes overlap: amplification, wavelength evolution, spectra, cavities, filaments, energy transfer and saturation time/mechanism | Blocked by Q-019/Q-022/Q-023 |
| Q-029 | Hall-extension qualification | Separately named derived mode with linear/nonlinear Bell and shock-front tests. The authorized target requires an `implemented_and_qualified` outcome | Verifying: opt-in experimental CT source and host odd-in-`alpha_H` manufactured-source smoke pass; predecessor `q029_hall_extension_source_local_candidate_2026-05-30.json` freezes the exact implemented additive CT-current normalization and dimensionless `chi_H` form as synthetic-contract chronology. Successor `q029_hall_bell_linear_source_local_preparation_2026-05-30.json` adds a separately named runnable source-local 1D/2D/3D launch preparation, a guarded shared Q-023 seed-carrier header, exact runtime geometry guards, a warnings-as-errors compiled host contract and 45 epsilon-specific positive-`chi_H` materializations. Fresh combined `mhd_w_bcc` 1D/2D/3D startups plus exact serial 2D uninterrupted-versus-restart continuation array parity pass locally. `q029_hall_bell_q022_prerequisite_successor_2026-05-30.json` adds a dedicated fail-closed `XCMP-EXT-HALL-BELL` Bai provenance route with empty equation-map and tolerance placeholders. `q029_hall_bell_linear_raw_extractor_source_local_2026-05-30.json` adds pinned-decoder Orion replay of digest-verified byte copies, rejecting symlinks and special inventory entries while retaining combined magnetic and fluid-velocity geometric seed-carrier projections and explicitly omitting a physical Bai oracle, reviewed coefficient grid and numeric tolerance. This remains nonqualifying launch, restart and raw-projection preparation only; reviewed physical Bai mapping, extracted comparison data, numeric tolerances, qualifying Hall-dispersion analysis, clean-candidate timestep freeze, linear/nonlinear Bell, shock-front, GPU and decomposition qualification remain open |
| Q-030 | Independent matched-code comparisons | Close explicit sub-gates Q-030-A Athena/Bai matched Bell-shock observables, Q-030-P PLUTO matched conservative coupling and Q-030-M MPI-AMRVAC overlapping AMR-shock observables, with frozen mappings, quantitative residuals and discrepancy reports | Blocked by Q-004/Q-009/Q-011/Q-019/Q-022/Q-023 |
| Q-030-A | Athena/Bai matched Bell-shock comparison | Frozen Athena/Bai mapping, matched Bell and shock observables, quantitative residuals and discrepancy report | Blocked by Q-011/Q-019/Q-022/Q-023 |
| Q-030-P | PLUTO matched conservative-coupling comparison | Frozen PLUTO mapping, matched conservative-coupling observables, quantitative residuals and discrepancy report | Blocked by Q-004/Q-022/Q-023 |
| Q-030-M | MPI-AMRVAC overlapping AMR-shock comparison | Frozen MPI-AMRVAC mapping, overlapping AMR-shock observables, quantitative residuals and discrepancy report | Blocked by Q-009/Q-011/Q-022/Q-023 |
| Q-031 | Bai et al. 2019 nonlinear CRSI comparison | Growth, saturation, spectra, diffusion, pitch-angle evolution, 90-degree crossing, isotropization and resolution dependence | Blocked by Q-020/Q-022/Q-023 |
| Q-032 | Ion-neutral-damped CRSI extension qualification | Plotnikov et al. matched comparison with damping-rate dependence. The authorized target requires an `implemented_and_qualified` outcome | Verifying: reduced static-neutral exact friction map, manufactured-source oracle and bounded Alfvén-envelope damping scan pass; `q032_reduced_static_neutral_runtime_local_2026-05-30.json` adds a separately named exact thin-2D3V source-local runtime carrier and retained Orion-local control-versus-damped mechanics probe with density and longitudinal residual zero and transverse/energy residuals below `1.0e-7`. `q032_plotnikov_reduced_map_applicability_derivation_2026-05-30.md` separates the implemented transverse momentum sink `exp(-nu_in dt)` from the conditional high-frequency wave-amplitude envelope `exp(-nu_in t/2)` and wave-energy envelope `exp(-nu_in t)`, while leaving the matched equation map, phase-scrambling disposition, extracted data and thresholds fail-closed. Shared parser and wrapper guards now admit this map only on the Newtonian single-fluid MHD source-task path and reject SR, GR, dynamical-GR, radiation, hydro and alternate ion-neutral task-list compositions. The MPI-ready harness admits an exact two-way x1 decomposition, but registered compute-node replay remains open. Preserve `q032_q033_two_rank_mpi_host_successor_2026-05-30.json` only as rejected chronology because OLCF policy forbids parallel login-node launches. Plotnikov-matched damped CRSI, nonlinear, registered MPI, Frontier HIP/GPU and review qualification remain open |
| Q-033 | CRPAI transport-calibration extension qualification | Adaptive-delta-f, physical-damping, `nu_eff`, anisotropy, spectra, quasi-steady-state and scaling comparison against Sun, Bai and Zhao. The authorized target requires an `implemented_and_qualified` outcome | Verifying: extension-only global bi-kappa fit, schema-7 restart state, parser guards, exact host restart oracle, adaptive-cadence sensitivity scan and the bounded active-MHD expanding-box endpoint-normalized analytic-source plus reduced-damping composition smoke pass locally; `q033_expanding_box_adaptive_damping_restart_resilience_bounded_local_2026-05-30.json` adds a passing coupled-carrier restart-continuity regression without unintended adaptive refit, and `q033_crpai_transport_calibration_q022_prerequisite_successor_2026-05-30.json` carries the exact blocked Q-022 Sun-Bai-Zhao provenance, equation-map and empty tolerance-table placeholder hashes through that synthetic bundle. `q033_crpai_transport_runtime_local_2026-05-30.json` adds a separately named deterministic thin-2D3V source-local runtime carrier with seeded transverse modes, antipodal bounded prolate particle pairs and immutable initialized-plus-post-step extraction diagnostics. The MPI-ready harness supports an opt-in widened two-block carrier, but registered compute-node replay remains open. Preserve `q032_q033_two_rank_mpi_host_successor_2026-05-30.json` only as rejected chronology because OLCF policy forbids parallel login-node launches. The bounded fixtures do not show variance reduction; reference extraction, physical calibration, effective scattering, saturation, registered MPI, Frontier HIP/GPU and external review remain open |
| Q-034 | Unsupported-capability review | Execute each register row's `parser_reject`, `rename_proxy`, `docs_limit` or `extension_gate` handling and archive its closure artifact | Verifying: local parser-reject, rename-proxy, docs-limit and extension-gate classifications are archived in `q034_unsupported_capability_local_closure.md`; successor active-MHD expanding-box fail-closed compositions are archived in `q034_unsupported_capability_successor_2026-05-30.json`; scientific extension qualification and external review remain open |
| Q-035 | Scoped state-of-the-art sign-off | Reviewer-approved wording tied to qualified claim IDs, exact references, regimes, metrics, performance evidence and exclusions | Blocked by Q-018/Q-023/Q-034/Q-040: an explicitly unapproved fail-closed wording draft is archived in `q035_scoped_state_of_the_art_draft_2026-05-30.json`; it is not manuscript text or a result |
| Q-036 | Crash consistency and I/O failure handling | Atomic restart publication, completion markers, last-known-good recovery, timeout, interrupted-write, truncated-restart and write-failure tests | Verifying: atomic writer, completion markers, manifests, startup checksum validation, truncated-artifact rejection and unwritable-target guards pass locally; deterministic serial-host schema mutation, bounded short-write, seek-failure, `/dev/full`, killed-writer residue and soft-wallclock continuation-parity drills pass; MPI per-rank publication, node-loss, scheduler-pretimeout and Frontier-filesystem matrix remain open |
| Q-037 | Durable evidence export | Record the user-selected Orion-only bulk-evidence deviation, checksum verification, retention risk, access policy and restore drills; obtain external disposition before terminal sign-off | Verifying: Orion is selected as the sole bulk-evidence root, the private PDF restore drill passes and Kronos is removed from scope; Orion-only retention is not durable archival, and external reviewer disposition remains open |
| Q-038 | Frontier environment profile selection | Minimum-supported baseline plus controlled XNACK/OFI A/B evidence, redacted allowlist capture and restricted-submit-export rank-0 propagation smoke test | Verifying: registered F0 baseline captures `HSA_XNACK=0`, `MPICH_GPU_SUPPORT_ENABLED=1` and `SLURM_EXPORT_ENV=ALL`, and maps eight ranks to eight distinct `ROCR_VISIBLE_DEVICES`; hardened F1 dirty-candidate gyro chronology plus registered clean-candidate gyro-v3 and paper-coupling-v2 slices capture the same allowlist, a visible rank-0 GPU, `libamdhip64`, `libmpi_amd` and `libmpi_gtl_hsa`; source-alias successor `6cbbbbd6` preserves the closed profile selector, strict Frontier value validation, Bash-startup sanitization, `sbatch --export=NIL`, scheduler `env -i` boundary, descriptor-based redacted allowlist capture and canonical 19-component `MODULEPATH`, then closes the `/ccs/home/...` source-alias clean-build blocker. Recovery successor `f2ad817a` preserves those controls and narrowly accepts scheduler account spellings `AST207` and `ast207` while rejecting unrelated account drift. Historical compute-snapshot successor `e8e47ead` preserves the same controls and passes a fresh structured clean-candidate F0 launch with `frontier_minimum_supported`, `HSA_XNACK=0`, `MPICH_GPU_SUPPORT_ENABLED=1`, `SLURM_EXPORT_ENV=ALL` and the canonical 19-component `MODULEPATH`. Exact clean freeze `31be2cd6` is authorized; matched XNACK/OFI A/B and Athena communication evidence remain open |
| Q-039 | GIZMO/RSOL comparison decision | Before P7A interpretation, archive a bounded Ji-Hopkins/GIZMO comparison with preregistered metrics or a reviewer-approved documented exclusion from qualification scope | Verifying: fail-closed local decision scaffold prohibits manuscript interpretation; Ji-Hopkins, official GIZMO documentation and cosmic-ray notes are privately archived in Orion with post-copy checksum verification, and the public repository HEAD is frozen; mapped bounded comparison or reviewed exclusion and external disposition remain open |
| Q-040 | Independent-comparison aggregate closure | For non-Hall production Bell close Q-028; for nonlinear CRSI close Q-031; for matched code comparisons close Q-030-A, Q-030-P and Q-030-M. The selected extension scope requires implemented-and-qualified outcomes for Q-029/Q-032/Q-033. Archive the discrepancy ledger and scoped conclusions; Q-039 still permits a reviewed bounded exclusion | Blocked by Q-022 and the named child gates |
| Q-041 | Undriven CRPAI nonlinear saturation | Undriven branch evolution, anisotropy, spectra, scattering, saturation, effective-damping trends, sensitivity matrix, matched reduced nonlinear full-f controls and archived weight-validity envelope, distinct from Section 5.6 linear reproduction and physical-damping calibration | Blocked by Q-007 |
| Q-042 | AthenaK architecture and style conformance | Review every implementation batch against nearby AthenaK patterns; run targeted formatting/lint checks and the repository style baseline; archive evidence for any accuracy- or efficiency-driven exception; introduce no new style violations or avoidable PIC-only infrastructure | Verifying: predecessor host architecture sidecars plus `q042_schema7_architecture_docs_successor_2026-05-30.json` cover the schema-6, shock-control, comoving-flux, schema-7 provenance/spectrum, repeated-AMR, warning-free Sphinx and bounded serial-host standard-MHD compatibility tranches, including one-cycle blast-AMR, shearing and orbital-advection launches; focused host regressions, targeted touched-file C++ lint and `git diff --check` pass; predecessor control-plane digests are preserved as chronology. Source-alias digest `6cbbbbd6` retains the outer-anchor serialization and canonical-MODULEPATH controls, closes the exact clean-build mount-alias blocker and publishes authorized immutable freeze `31be2cd6`. Recovery successor `f2ad817a` adds narrow scheduler-account normalization, immutable cross-generation terminal-recovery handoffs and historical predecessor-inventory verification. Historical compute-snapshot successor `e8e47ead` adds a dedicated numeric scheduler-ID compute launch boundary with descriptor-pinned exact-byte ledger semantics, authorized-root traversal for new and inherited parent descriptors, transient pathname-ABA rejection and retained login-side writer-lock discipline; it passes 261 control-plane tests, 76 PIC publication tests, 5 Frontier publication tests, 8 readiness-registry tests, static checks, a live 28-record snapshot probe and three independent reviews. Fresh structured clean-candidate F0 execution passes through installed job `4745842`. Non-PIC MPI/multilevel runtime slices and external review remain open |

Historical Q-027/Q-038/Q-042 successor `e8e47ead` preserved the
source-alias predecessor `6cbbbbd6`, its exact clean build-profile publication
and immutable clean-candidate freeze `31be2cd6` with source bundle `e11dc8fb`
and executable `bea2a418`. The first structured F0 parser-contract submission
through `6cbbbbd6` created held job `4745523`, then rejected scheduler account
spelling `ast207` before releasing the hold. The trap cancelled the job before
execution, retained the durable marker and blocked later PIC submissions.
After scheduler age-out purged live `scontrol` and `squeue` state, reviewed
successor `f2ad817a` published immutable mirrored recovery handoff `06392534`,
appended synthetic attachment event `f518c3d8` and zero-consumption
reconciliation `a13311b3`, and removed the pending marker. Its paired immutable
install, historical predecessor-inventory verification, 250 control-plane
tests, 76 PIC publication tests, 5 Frontier publication tests, static checks
and three independent reviews pass. Its fresh structured F0 retry attached and
released job `4745755`, then failed before Athena when the compute-node Project
Home mount returned `OSError 524` for the login-side mirrored writer lock.
Reviewed successor `e8e47ead` retains writer locks for login-side reads and
mutations while its dedicated numeric scheduler-ID trampoline uses an
authorized-root-traversed descriptor-pinned exact-byte read-only snapshot.
Its paired immutable install, 261 control-plane tests, static checks and three
independent adversarial reviews pass. Mirrored policy SHA `cfe6610a` and
promotion SHA `6d6bf1cd` were live at that transition. Fresh structured F0 job `4745842` completed,
passed parser analysis and reconciled `0.0025` node-hours. The coherent
31-record Orion, receipt and Project Home streams report cumulative
consumption `0.07833333333333331` node-hours and no active reservation.
Applicable registered-science prerequisites remain open.

Staged registered-F1 successor `6002c80e` preserves the active `e8e47ead`
chronology and narrows the next execution boundary to two exact one-attempt
clean-candidate slices: GPU relativistic gyro and GPU paper-mode coupling. Its
structured artifact freeze, immutable inventory, descriptor-pinned offline
analysis and qualification-time no-write recomputation are locally implemented.
Fresh adversarial review found two residual gaps: qualification closed the
inventory, result and receipt descriptors before recomputation, and the coupling
analyzer ignored extra files outside each expected coefficient prefix. Both are
repaired: qualification now retains and rechecks the three evidence descriptors
plus their directory ancestry, and coupling closure rejects every extra file
beneath `output/`. A subsequent launch-publication retest found that the
trampoline released nested-directory identity and its newly created `analysis/`
descriptor before final verification. The repaired successor retains both
through final publication, rejects byte-identical nested replacement and
requires the original empty owner-only analysis directory. The exact staged
bindings are independently recomputed by the readiness registry. A final
operability retest then found ambient restrictive-`umask` dependence during
artifact creation. The repaired launch lifecycle scopes a deterministic
`umask 022` across artifact-root, workload and publication creation and restores
the inherited value on exit. Paired production install, active-policy promotion
and both registered F1 executions remain pending clean commit curation. Fresh
independent launch-publication and qualification-domain exploit retests pass for
superseded pre-activation install `140e9a29`. An adjacent-path review then
rejected promotion before any live-policy mutation and staged `6002c80e`: it
publishes `analysis/` through a descriptor-pinned staging rename, rejects empty
structured subtrees and anchors qualification traversal at the trusted PIC root
so replaceable `runs/` and manifest `snapshot/` components cannot become trust
roots. Fresh launch-publication and qualification-domain exploit retests pass
for `6002c80e`. Curated commit `c1be6cab`, identical immutable paired install
with inventory SHA `5144e836`, reviewed policy SHA `ee3923fb` and active
promotion SHA `433d0afc` pass. The existing paired genesis anchors are
preserved without append; this successor does not require another anchor
migration. Both registered F1 executions remain pending.

Registered gyro attempt `4746123` then completed and reconciled
`0.0030555555555555557` node-hours, but scientific analysis rejected the
immutable run before result publication: the required Frontier profile emits
reviewed Cray MPICH version and environment diagnostics to stderr while the
first analyzer required empty stderr. Preserve that run as unqualified
chronology. Retry registrations bind the exact reviewed Cray MPICH
informational transcript SHA-256, locally archive its immutable-inventory
provenance in
`q027_frontier_f1_failed_gyro_mpich_stderr_provenance_2026-05-30.json`, reject
any stderr byte drift and rotate both one-attempt authorization IDs before
retry promotion.

Before retry promotion, an adjacent-path review staged successor `cbc6fb50`.
It preserves the active `6002c80e` chronology, anchors launch and qualification
below the site-owned serialization root, retains observed artifact-directory
ancestry and workload-created directory identities through final publication,
rejects regular-file namespace substitution during freeze, and explicitly
records the irreducible same-account process-isolation
prerequisite from `mkdirat` through no-follow descriptor and retained-ancestry
binding, plus workload-created descendant isolation through the post-action
recursive-capture handoff. Fresh exploit retest, clean commit `559d5e54`,
paired install inventory SHA `f5268721`, reviewed retry-policy SHA `58f85503`
and promotion SHA `05dfbf32` pass. Before either v2 submission, a late
independent review reproduced generated-inventory replacement between creation
close and read-only reopen plus offline initial-loader nested-directory and
regular-file identity substitutions. Staged repair successor `d55d064d`
retains the generated inventory creation descriptor through final publication
verification and binds offline pre-open metadata to opened descriptors while
retaining the observed inventory bytes, payload-file identities and published
analysis-result descriptors through receipt publication. Fresh launch and
offline-analysis exploit retests pass. A later launch-publication review then
reproduced byte-identical top-level and nested frozen workload-payload namespace
replacement after inventory freeze. Intermediate successor `d8cc7d61`
retained those identities but an exact review reproduced payload and inventory
replacement after their individual checks while later verification continued.
Intermediate successor `d538af91` retained descriptors through a closing
namespace and byte sweep but exact review reproduced nested-directory
transplant between directory rebinding and the payload sweep. Repaired
successor `4ccde8df` rebinds directory identities afterward and rechecks the
generated inventory. Fresh exact-successor retest, curated commit `6fa926a9`,
paired immutable install inventory SHA `ea45cdc5`, existing-genesis migration,
reviewed policy SHA `3895221e` and promotion SHA `016cbb66` pass.
One corrected-controller gyro manifest then failed closed before reservation
because the operator froze the attestation template's seven-field queue format
instead of the validator's six-field snapshot format. Fresh attested gyro-v2
job `4746245` completed and reconciled `0.0030555555555555557` node-hours, but
offline analysis rejected before result publication: Athena deterministically
published cycle-zero through cycle-three particle VTK artifacts while the
staged analyzer registered only cycle two, and read-only replay then exposed an
abbreviated synthetic parser fixture missing emitted `POINT_DATA`, `cr_source`,
`macro_weight` and `birth_time` records. Preserve that run as unqualified
chronology. The parser repair registers both exact deterministic output trees,
requires the emitted particle-VTK schema, measures cycle two and rotates only
the consumed gyro authorization to `f1-clean-gyro-mpich-stderr-v3`. Unused
paper-coupling v2 keeps its authorization ID while receiving the repaired
analyzer binding before its first serialized execution boundary. Commit
`c878dd6c` was promoted with mirrored policy SHA `0538da6d` and promotion SHA
`f84fc7b7`. Attested gyro-v3 job `4746290` then completed, reconciled and
initially froze pending-review qualification manifest SHA `9f7a5d64`. The first paper-coupling
manifest creation failed closed before manifest publication or ledger intent
because the submission scrubber matched ordinary analyzer identifiers.
Scrub-safe commit `ae9fbb97` retained the unused coupling-v2 authorization ID
and promoted mirrored policy SHA `559a4b1e` with promotion SHA `84065c8a`.
Attested paper-coupling-v2 job `4746297` completed, reconciled and froze
pending-review qualification manifest SHA `6262332a`. After that active-policy
replacement, gyro-v3 replayed through current-policy successor qualification
manifest SHA `17265e02`; preserve initial gyro manifest `9f7a5d64` as immutable
chronology and use `17265e02` for terminal review. The ledger is coherent
at 43 records with cumulative consumption `0.09388888888888887` node-hours,
zero active reservations and no pending PIC marker. Source-local accepted
closure index
`q027_frontier_f1_accepted_closure_source_local_evidence_2026-05-30.json`
binds both live execution projections for clone-local review. The opt-in live
registry additionally validates exact mirrored ledger chains and replays both
terminal manifests through descriptor-pinned bytes.
Bounded F2 v1 job `4746310` remains rejected chronology because its launch
contract allowed Athena error history to escape the registered artifact tree.
Repaired F2 v2 policy SHA `446db26a` and promotion SHA `b772da00` route that
output below the run tree. Attested job `4746316` completed, reconciled
`0.0025` node-hours and passed descriptor-pinned no-write replay with inventory
SHA `ca673515`, analysis SHA `95dc6e40`, receipt SHA `1685be75` and terminal
qualification SHA `09e71dbd`. The final ledger is coherent at 49 records with
cumulative consumption `0.09944444444444443` node-hours, zero active
reservations and no pending PIC marker. After the F2 policy replacement, both
registered-F1 executions replayed through current-policy successor
qualification manifests: gyro-v3 SHA `d77b88f1` and paper-coupling-v2 SHA
`df0876df`. Preserve the older F1 qualification freezes as immutable chronology.
Source-local accepted closure index
`q027_frontier_f2_accepted_closure_source_local_evidence_2026-05-30.json`
binds accepted F2 evidence plus fail-closed queue-format and transient-drift
chronology for clone-local review.

Strict registered-lifecycle successor `6f3458ca` from operational commit
`fe2d7f7c` then closed adjacent replay holes found during independent review:
registered transition arithmetic and job-ID ownership are enforced, Slurm
walltime parsing rounds conservatively, terminal-recovery handoff bytes are
bound to every referencing record, purged zero-execution snapshots match the
issuance verifier and compute-side read-only snapshots pin and recheck mirrored
handoff files. Its retained 354-test matrix, static checks and two fresh
independent adversarial reviews pass. The identical read-only paired install
has inventory SHA `10bba8e7`; reviewed policy SHA `36228012` and promotion SHA
`a2e04161` are active. With the trusted Frontier queue empty, the reviewed
Q-016 direct-`srun` accounting-only authorization imported seven rows for jobs
`4746332`, `4746335`, `4746336`, `4746337`, `4746341`, `4746342` and
`4746343`. All three ledger streams are coherent at 56 records with cumulative
consumption `0.11055555555555553` node-hours, zero active reservations and no
pending PIC marker. Preserve
`q027_manual_frontier_accounting_activation_2026-05-30.json` as the current
operational activation record; these imported rows are ineligible for
scientific evidence.

## Immediate Agent Handoff: First Actions

Future implementation agents should execute the following successor sequence:

1. Read this plan, the paper source and all `PIC-P0-*` evidence locations.
2. Confirm the current `PIC` branch, commit and clean-worktree state. Treat the
   retired large-machine note as historical evidence only.
3. Preserve predecessor readiness records as chronology and add successor
   sidecars for schema-6 restart, schema-7 provenance, shock controls, comoving-flux expanding-box
   repair, unsupported-capability handling and architecture review.
4. Preserve the local independent adversarial PASS for historical successor
   control-plane digest `12750ea6` as chronology. Outer-anchor digest
   `80c0797b` passed its combined-domain clone and hardlink retest, then was
   activated in a proved-quiescent cutover. Its first clean build-profile
   preflight rejected caller-dependent Frontier `MODULEPATH` provenance before
   artifact creation. Reviewed successor `c8002a1d` from operational
   commit `37ccee56` closed that drift, but the exact clean-build retry then
   rejected the authenticated `/ccs/home/...` source alias before artifact
   creation. Reviewed successor `6cbbbbd6` from operational commit
   `5f1458e3` is paired-installed and promoted with its existing-genesis anchor
   migration complete. Its exact clean build and immutable freeze `31be2cd6`
   pass. Its first structured F0 submission produced held job `4745523`, then
   failed closed on scheduler account spelling before execution. Historical
   recovery successor `f2ad817a` from operational commit `3df7deee` publishes
   immutable cross-generation handoff `06392534`, reconciles that purged job
   at zero node-hours and preserves the exact freeze. Its fresh structured F0
   retry attached job `4745755`, then failed before Athena on compute-node
   `OSError 524`. Historical successor `e8e47ead` from operational commit
   `a110e38c` uses the reviewed descriptor-pinned read-only compute snapshot,
   preserves login-side writer locks and passes fresh structured F0 job
   `4745842`; preserve it while closing applicable registered-science
   prerequisites.
5. Preserve the curated clean source commit series, exact clean candidate
   `31be2cd6`, executable digest `bea2a418`, historical clean-freeze
   authorization promotion `3e7e8c7f` and historical predecessor promotion
   `6d6bf1cd`.
6. Preserve freshly retested registered-F1 successor `6002c80e`, curated commit
   `c1be6cab`, identical paired-install inventory SHA `5144e836`, reviewed
   policy SHA `ee3923fb` and historical promotion SHA `433d0afc`. The paired
   genesis anchors were preserved without append.
7. Preserve historical successor `cbc6fb50`, clean commit
   `559d5e54`, paired install inventory SHA `f5268721`, reviewed policy SHA
   `58f85503` and promotion SHA `05dfbf32` as chronology only; preserve one
   fail-closed v2 pre-submit manifest, with no reservation or scheduler
   submission created under it. Preserve staged repair successor `d55d064d`
   as superseded chronology after late frozen workload-payload identity review;
   preserve `d8cc7d61` and `d538af91` as superseded post-check-window
   chronology; preserve freshly retested historical repaired successor `4ccde8df`,
   curated commit `6fa926a9`, paired-install inventory SHA `ea45cdc5`,
   reviewed policy SHA `3895221e` and promotion SHA `016cbb66`; preserve one
   fail-closed operator queue-format manifest and reconciled unqualified gyro-v2
   job `4746245`. Preserve promoted parser policy SHA `0538da6d` and promotion
   SHA `f84fc7b7`, scrub-safe commit `ae9fbb97`, replacement policy SHA
   `559a4b1e` and promotion SHA `84065c8a`. Preserve completed gyro-v3 job
   `4746290` and paper-coupling-v2 job `4746297`, their serialized ledger
   reconciliations, snapshotted offline-analysis publications, no-write
   qualification replays and immutable pending-review qualification manifests
   `17265e02` and `6262332a`; preserve initial gyro freeze `9f7a5d64` as
   immutable chronology after the coupling-policy replacement. Preserve
   clone-local accepted-closure projection
   `q027_frontier_f1_accepted_closure_source_local_evidence_2026-05-30.json`
   alongside the immutable Orion evidence. Preserve rejected F2 v1 job
   `4746310`, repaired F2 policy SHA `446db26a`, promotion SHA `b772da00`,
   fail-closed pre-reservation queue-snapshot chronology and accepted F2 v2 job
   `4746316`. Use current-policy F1 qualification manifests `d77b88f1` and
   `df0876df`, bounded F2 qualification manifest `09e71dbd` and clone-local F2
   closure projection
   `q027_frontier_f2_accepted_closure_source_local_evidence_2026-05-30.json`
   for terminal review while retaining older F1 freezes as chronology.
   Preserve strict registered-lifecycle successor `6f3458ca`, source commit
   `fe2d7f7c`, policy SHA `36228012`, promotion SHA `a2e04161`, paired inventory
   SHA `10bba8e7` and Q-016 accounting-only activation sidecar
   `q027_manual_frontier_accounting_activation_2026-05-30.json`. The current
   coherent mirrored ledger has 56 records, cumulative consumption
   `0.11055555555555553` node-hours, zero active reservations and no pending PIC
   marker.
8. Preserve the completed Orion-local Q-016 two-rank MPI decomposition replay
   and close remaining local preparation where feasible: Q-017 kernel
   telemetry, Q-018 immutable claim links, Q-022/Q-023
   extraction/reference-threshold scaffolding and the Q-039 mapped or
   reviewed-exclusion route. Preserve the rejected direct-PALS Q-032/Q-033
   login-host attempt only as policy-invalid chronology, then replay it solely
   through a registered compute-node allocation. Execute the remaining Q-009
   MPI/HIP stress and Q-016 HIP
   parity slices only as separately registered Frontier work.
9. Preserve predecessor control plane `e8e47ead`, the original
   genesis anchor, immutable recovery handoff `06392534` and mirrored historical
   predecessor promotion `6d6bf1cd` as chronology. After the registered-F1 cutover, use only
   the promoted installed successor runner for subsequent
   registered Frontier work and keep every reservation serialized through the
   audited ledger path.
10. Only after all local and control-plane gates pass, build and run Frontier validation
   under the `debug`-preferred, `normal`-fallback, budget-tracked procedure above.
11. Reproduce every paper result, execute registered nonlinear and cross-code
   campaigns within authorization, and complete the production sign-off bundle.
12. Use `normal` on `batch` for registered long saturation, full shock or
    controlled-scaling work after prerequisites close. Stop and ask the user for
    permission before exceeding 10000 cumulative node-hours or expanding beyond
    the authorized QOS, partition, account or campaign envelope.

### Current Q011 Stage-4 no-science publication boundary

The current Q011 action is to publish the reviewed `problem/ps_p0=1.0`
pressure selection, not to launch a qualifying campaign. The live paired
controller `930a04d1` and candidate-only policy retain an empty
registered-science allowlist, closed admission smoke, no pending submission,
and no incomplete manual-accounting marker. Preserve the exact earlier human
selection and rationale verbatim. As a separate scientific scope
interpretation, the selected case is a provenance-first baseline because it
matches Bai et al. (2015)'s explicit P0=T0=1 normalization; the short pressure
pilots do not establish pressure independence or a physics-preferred Bell/DSA
baseline.

Before production publication, require the exact Stage-4 source commit to be
pushed, independently rereviewed, and passed by the clean-snapshot
pressure-gate worker. Run machine reanalysis, human sealing, publication, and
any reconciliation only from one authenticated read-only Git archive of that
exact commit. Machine reanalysis must stop without creating a reviewer
attestation or candidate receipt. Publication requires a separately supplied
post-reanalysis human decision record that binds the exact sealed reanalysis
and has a strictly later review timestamp. The software authenticates the
decision record, not human authorship, so the operational human stop remains
mandatory. Require the sealed Stage-4 preparation-source attestation, human
decision, candidate publication authorization, and controller-state
attestation to prove that preparation, sealing, and publication used the same
replacement-ref-disabled authenticated publisher archive.
Accept the result only after independently verifying the canonical schema-v3
receipt, sealed reanalysis and reviewer attestations, launch-prohibited
controller-state attestation, inode-bound success seal, absent publication
guard, unchanged empty allowlist, and unchanged active policy/promotion.
Stage-4 publication grants no science authority.

## Initial Change Log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-05-25 | Created comprehensive production-readiness and paper-reproduction plan | Expanded narrow large-machine validation scope after full code/paper review exposed model-level blockers and Frontier execution requirements |
| 2026-05-25 | Added Entity Toolkit comparative audit and bounded adoption gates | Adopt shared particle verification, provenance, spectra, restart and observability practices without importing incompatible electromagnetic PIC physics into paper mode |
| 2026-05-30 | Retired the narrow large-machine note as historical-only evidence; corrected AMR wording; added claims registry, nonlinear saturation, cross-code, statistical, resilience, portability, archive/licensing and Frontier authorization-boundary gates | A fresh source and plan audit showed that paper reproduction plus exploratory proxies alone could not support production or scoped state-of-the-art claims |
| 2026-05-30 | Enforced engineering-proxy classification in checked-in artifact tooling and hardened the plan after independent physics, taxonomy and Frontier-operability reviews | Legacy publication identifiers now carry explicit unqualified metadata and visible figure labels; comparison gates, sampling rules, release profiles and fail-closed submission accounting are executable specifications rather than implied policy |
| 2026-05-30 | Added immutable Frontier control-plane snapshots, JSONL accounting with non-recursive mirror receipts, dual Slurm-ID guards, inert scan shells, root-document retirement banners and checksummed standalone-helper lineage companions | Adversarial execution review found that documentation-only policy was insufficient while mutable paths, ambient exports, copied shell scripts or partially checksummed quicklooks could bypass provenance and submission boundaries |
| 2026-05-30 | Replaced the initial `debug`-only Frontier rule with a `debug`-preferred, `normal`-fallback policy on `batch`; raised the cumulative testing cap to 4000 node-hours with a mandatory stop-and-ask boundary | The user authorized normal-QOS fallback when the single-user `debug` slot is unavailable and expanded the tracked cluster-testing budget while preserving serial PIC submissions and fail-closed accounting |
| 2026-05-30 | Added the AthenaK architecture/style-conformance contract and Q-042 | PIC implementation work should reuse AthenaK's existing framework and design choices unless archived evidence shows that a scoped deviation is necessary for material accuracy or efficiency |
| 2026-05-30 | Expanded the authorized target to include qualified Hall Bell, ion-neutral-damped CRSI and adaptive-delta-f physical-damping CRPAI extensions; raised the cumulative testing cap to 10000 node-hours; recorded `pending external review`; selected the then-current Project Home plus Kronos OLCF-side archive design pending exact path and off-site destination freeze | Preserve the superseded intermediate decision accurately: the user requested completion of the optional extensions, authorized the larger tracked Frontier budget and deferred named external review and final institutional-retention selection |
| 2026-05-30 | Implemented explicit paper momentum state, `C`-aware Boris mechanics, paper-versus-Hall induction separation, physical delta-f state/deposition/restart diagnostics, particle-plus-MHD expanding-box transforms, restart schema fingerprints, and atomic checksummed restart publication | Close initial implementation defects without prematurely closing analytical, AMR, GPU, Frontier, resilience or review gates |
| 2026-05-30 | Added a uniform-MHD expanding-box appendix invariant oracle, corrected multistage over-expansion by applying the exact MHD box map once after RK, made retained-particle AMR refresh explicit, and added opt-in geometric particle-aware AMR load costs | Convert box and AMR implementation assumptions into executable checks while keeping portability, MPI and Frontier measurements open |
| 2026-05-30 | Froze the accessible Entity replacement tree and bounded shared-kernel hashes, recorded the manuscript PDF render dependency blocker, and established the repository-wide C++ style baseline | Make Q-015, Q-026 and Q-042 prerequisites explicit and machine-readable without promoting incomplete comparison, archive or style gates |
| 2026-05-30 | Added an explicit `extended_mhd_pic` Hall-current CT manufactured-source smoke and documented its experimental `cE_CT = cE_ideal + alpha_H P_edge[J_CR]` normalization | Verify host source isolation with nonzero odd-in-`alpha_H` magnetic increments while keeping derived Hall Bell, nonlinear, shock-front, GPU and decomposition qualification open under Q-029 |
| 2026-05-30 | Added the bounded reduced static-neutral ion-neutral friction map, extension-only adaptive global bi-kappa delta-f fit, restart schema version 6 fitted-state persistence, and direct host manufactured-source/restart oracles | Land selected-extension mechanics without promoting them to damped-CRSI or CRPAI transport qualification before the Q-032/Q-033 comparison matrices close |
| 2026-05-30 | Recorded bounded host Debug serial, Debug MPI and Release MPI compilation plus serial runtime evidence; recalibrated only the short AMR multispecies interface-smoke floor after measuring the corrected nearest-center TSC response against `origin/PIC` | Preserve a nonzero AMR liveness regression without misclassifying a one-turn engineering proxy as an oscillation-frequency or MPI qualification oracle |
| 2026-05-30 | Added parser-contract, paper task-stage, paper-coupling conservation, continuum-gyro, Entity exact-overlap, preregistration, GIZMO fail-closed decision, extension-envelope, style-conformance, OpenMP and sanitizer evidence; staged the repaired Frontier control plane without installing or initializing it | Close locally executable evidence slices while preserving the then-pending authenticated Kronos, Frontier runtime, scientific-campaign and external-review boundaries |
| 2026-05-30 | Reconstructed the omitted MNRAS render dependencies in a private environment, archived a checksummed 20-page PDF bundle in Project Home and passed a fresh-directory restore drill | Close the local Q-026 PDF retrieval blocker without authorizing external redistribution, then-pending Kronos mirroring or terminal archive sign-off |
| 2026-05-30 | Archived bounded Q-005/Q-006 proxy characterization, Q-008 expanding-box gyro/CPAW oracles and Q-032/Q-033 extension scans with explicit non-closure limitations | Preserve quantitative local evidence while preventing engineering proxies and bounded fixtures from being promoted to paper, Frontier or external-review qualification |
| 2026-05-30 | Bound Project Home ledger mirroring to reviewed `filesystem_copy`, separated the then-pending authenticated Kronos `globus` or `dtn_rsync` bulk transfer selection, and broadened Project Home's non-permanent role to include private reference-artifact staging and restore drills | Resolve transport and archive-role ambiguity before any ledger genesis, reservation or Frontier submission |
| 2026-05-30 | Removed Kronos, DTN and Globus from the execution design at user direction; selected Orion as the sole bulk-evidence root with explicit durability risk; initialized mirrored ledger genesis `849bf340` with paired immutable snapshot `cde67c16`; installed hardened active snapshot `36761ee2`; copied and restore-verified the private PDF bundle from Orion | Unlock registered Frontier execution without misrepresenting Orion-only retention as institutional archival |
| 2026-05-30 | Passed registered Frontier F0 HIP/MPI admission job `4744225` after preserving and reconciling bounded bootstrap failures; fixed HIP device captures in the parallel-shock initializer; promoted paired control-plane snapshot `1635c9f6` with reservation-ledger and scheduled-job manifest-digest anchoring | Admit subsequent registered science work only after a real Frontier build, rank/GPU mapping and parser/startup smoke while preventing self-consistent post-reservation manifest replacement |
| 2026-05-30 | Passed registered Frontier F1 dirty-candidate one-rank GPU relativistic-gyro job `4744232` after a preserved preliminary pass and independent evidence-hardening review; captured visible-GPU binding, HIP/GPU-aware-MPI linked libraries, scheduled manifest-digest equality and the cycle-2 Boris oracle | Preserve bounded candidate evidence while keeping the clean-candidate freeze, full GPU matrix, MPI communication, performance and scientific-claim promotion gates open |
| 2026-05-30 | Advanced PIC restart schema to version 6, bound continuation-sensitive selectors and coefficients plus adaptive and star-potential state, made restart sync/close cleanup unconditional, repaired zero-particle-rank AMR load-cost collectives, and fingerprinted parallel-shock continuation controls | Close clean-freeze blockers exposed by independent implementation review while retaining MPI, restart-failure and Frontier qualification gates |
| 2026-05-30 | Repaired active-MHD expanding-box face storage into divergence-preserving comoving magnetic fluxes with derived physical views, mapped edge EMFs, directional divergence and CFL metrics; added oblique-divB and fail-closed guard regressions; rebuilt the then-current 56 Sphinx pages warning-free in an isolated declared dependency environment | Preserve bounded host correctness while keeping unqualified active-MHD compositions, full Appendix-A histories, coupled driven CRPAI, MPI, GPU and external review open |
| 2026-05-30 | Added a per-claim Q-018 link scaffold, compact Sun-Bai source anchor, AthenaK paper-mode equation profile, bounded Entity overlap map, Q-022/Q-023 schemas, per-campaign Q-023 draft bundle, plotting-lock candidate and an explicitly unapproved Q-035 wording draft | Freeze every locally available comparison and claim-control artifact while keeping missing external datasets, numerical thresholds, exact plotting lock, immutable qualifying bundles and reviewer dispositions fail-closed |
| 2026-05-30 | Archived the Q-009 retained-state inventory and selected the locally smooth `paper_smooth` cell-centered AMR-interface deposition path; kept the optional conservative policy and direct-staggered trajectory-current candidate outside qualified production scope | Make the AMR lifetime and refinement-interface policy explicit without overstating the existing one-shot host proxies as the required repeated MPI/HIP stress matrix |
| 2026-05-30 | Paused all new Frontier submissions after independent adversarial review reproduced mutable post-reservation snapshot replacement, forged scheduler attachment/reconciliation, launcher substitution/TOCTOU and interrupted-ledger recovery gaps in the staged successor control plane | Keep Orion unchanged and fail closed until the successor control plane is repaired, independently re-audited and promoted from a clean source candidate |
| 2026-05-30 | Added a separate fail-closed qualification-manifest validator/freezer that rejects proxy evidence, dirty source candidates, unknown claim IDs, escaped or missing files and checksum drift; landed bounded Q-017 driver-level task-list, output, AMR/load, throughput and safely derivable rank-load telemetry with explicit residuals | Prevent exploratory artifacts from being mislabeled as review-ready evidence and expose low-overhead observability without overstating driver-only instrumentation as complete performance qualification |
| 2026-05-30 | Repaired three successor control-plane defects found by independent adversarial review: preserved execute mode on snapshotted Athena, retained ledger-path roles when checking caller arguments, and bound receipt provenance to the configured Project Home mirror plus reviewed filesystem transport; local full-suite re-audit passes and fresh independent re-audit is in progress | Keep Orion unchanged and submissions paused until the repaired successor passes an independent exploit retest and is frozen from a clean source candidate |
| 2026-05-30 | Advanced PIC restart schema to version 7 with persistent source and birth-time provenance, explicit typed particle outputs, tracked-row metadata and bounded independently reconstructed weighted spectra; added a six-transition serial Q-009 refine/derefine lifetime stress through restart continuations; staged checksummed public comparison sources and GIZMO interpretation references for private Orion ingest | Close locally implementable Q-009/Q-016/Q-022 slices without promoting bounded host regressions or retrieved source references into publication evidence |
| 2026-05-30 | Added physical-volume particle deposition and built-in MHD history semantics for active expanding MHD; mapped non-delta-f conservative particle EM feedback into the final physical source frame; admitted the separate endpoint-normalized adaptive delta-f analytic source path under a narrow cell-centered envelope; added final-frame conservation and adaptive-delta-f reduced-damping source-ordering smokes | Remove the local Q-008/Q-033 composition-launch blocker without promoting bounded host mechanics into Appendix-A, CRPAI transport, MPI/HIP or external-review qualification |
| 2026-05-30 | Checked every sequential restart-header append and stdio offset seek; added deterministic `LD_PRELOAD` short-write, seek-failure and interrupted-writer injection plus `/dev/full`, schema-mutation and soft-wallclock continuation drills | Close bounded serial-host Q-024/Q-036 fault-injection slices while preserving MPI per-rank, node-loss, scheduler-pretimeout and Frontier-filesystem gates |
| 2026-05-30 | Passed the repaired successor Frontier control plane through its 98-test local suite, Python and shell static checks and independent adversarial retest at staged digest `12750ea6`; froze local Q-023 profiles for eight qualifying seeds plus two excluded pilots; repaired standard-MHD custom generators and isothermal timestep compatibility; corrected local-staging wording for Q-022/Q-039 | Preserve an honest local successor boundary before clean-candidate commit, paired Orion and Project Home promotion, registered Frontier execution and external review |
| 2026-05-30 | Hardened the next Frontier control-plane candidate after fresh independent reviews: closed qualification platform relabeling and run-artifact lineage gaps, bound reconciled node-hours, anchored explicit-genesis mirrored-ledger reads, rejected permissive candidate JSON and portable-layout aliases, made immutable trees and their parents durable, made qualification and replacement publication roll back after parent-directory sync failures, synced ledger first-create and CSV replacement, preserved trusted lexical Project Home mount aliases below an anchored root, anchored output-parent creation and atomic writes below trusted roots, bound staged-policy digest equality, serialized mirrored-policy promotion, added an installed Orion build-profile writer and exposed the three closed pre-selection Frontier environment profiles with redacted capture; its 118-test local suite, 16-test ledger suite, 19-test qualification suite, 199-test publication discovery, warning-free docs and static checks pass at staged digest `d8820c33` | Preserve the live `3e933edd` install as chronology and pause all submissions until this exact candidate passes fresh independent adversarial retest, paired install, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Repaired four additional activation blockers reproduced by a fresh independent audit of `d8820c33`: require an explicitly open active policy before ledger genesis creation, validate immutable genesis anchors before mirror-repair writes, require coherent mirrored state before CSV publication, and close plus unset the runtime-allowlist descriptor before executing the workload; the repaired digest `121f9f8f` passes its 125-test control-plane suite, 22-test ledger suite, 20-test qualification suite, 8-test readiness-registry suite, 213-test publication discovery and static checks | Preserve the live `3e933edd` install as chronology and pause all submissions until `121f9f8f` passes fresh independent adversarial retest, paired install, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Preserved `121f9f8f` as superseded chronology after an independent audit found residual writeability-preflight ordering, interrupted fresh-bootstrap and completed-cancellation marker recovery gaps; hardened those paths in recovery successor digest `7b722ab2` and added focused regressions | Keep Orion unchanged and submissions paused until `7b722ab2` passes the full local gate, fresh independent adversarial retest, paired install, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Preserved `7b722ab2` as superseded chronology after an independent audit found stranded completed-attachment markers, mutable-source mutating-repair entry points and inherited Bash-startup descriptor duplication; hardened those paths in recovery-provenance successor digest `7ea4aa8b` and expanded exact interrupted-bootstrap prefix coverage | Keep Orion unchanged and submissions paused until `7ea4aa8b` passes the full local gate, fresh independent adversarial retest, paired install, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Preserved `7ea4aa8b` as superseded chronology after an adjacent-path sweep found that reconciliation did not verify its paired Project Home install before scheduler accounting queries or mirrored ledger append; hardened that path in accounting-provenance successor digest `1843ce72` and added a fail-before-query regression | Keep Orion unchanged and submissions paused until `1843ce72` passes the full local gate, fresh independent adversarial retest, paired install, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Preserved `1843ce72` as superseded chronology after an independent audit found bootstrap-import, mutable-source-install, ambient-Git, scheduler-environment, trampoline-path and module-provenance activation blockers; hardened those paths in bootstrap-scheduler successor digest `7397038b`, added the isolated exact-inventory runner, reviewed HEAD-blob installation, hermetic Git reads, `env -i` Frontier scheduler commands, `sbatch --export=NIL`, descriptor-anchored artifact finalization, exact canonical module provenance and an explicit paired-install lifecycle assertion, then passed its 198-test control-plane-plus-ledger suite, static checks, warning-free 58-page Sphinx build, targeted touched-file C++ lint and 274-test publication discovery | Keep Orion unchanged and submissions paused until `7397038b` passes fresh independent adversarial retest, clean commit curation, paired install, anchor migration, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Integrated bounded local Q-008/Q-033 expanding-box coupled-carrier restart parity, Q-009 coupled reflecting/outflow Debug-plus-UBSan lifetime stress, Q-011 frozen Section 5.4 preparation deck and fail-closed analyzer, Q-017 particle-wrapper telemetry and Q-022 fail-closed comparison map/tolerance/provenance scaffolds | Close every locally executable preparation slice without promoting bounded host mechanics, empty comparison thresholds or unexecuted paper campaigns into publication evidence |
| 2026-05-30 | Preserved `0d07f132` as superseded chronology after an independent combined-domain clone probe showed that replaceable Orion and Project Home PIC generations could still bypass inner locks; staged outer-anchor successor digest `80c0797b` serializes ledger mutation and policy promotion through a descriptor lock on site-owned `/lustre/orion/ast207`, adds direct contention regressions, and passes its 208-test control-plane suite, 36-test focused ledger suite, paired temporary install, 286-test publication discovery and independent clone/hardlink retest | Activate only during a proved-quiescent cutover because predecessor snapshots do not take the outer lock; keep submissions paused until paired install, anchor migration, active-policy promotion and clean-candidate science freeze |
| 2026-05-30 | Activated outer-anchor snapshot `80c0797b` in a proved-quiescent cutover: paired immutable Orion and Project Home install, one-time genesis-anchor migration, coherent active-policy promotion and idempotent post-promotion anchor validation pass while all three ledger streams remain at 22 records | Preserve `q027_outer_anchor_paired_activation_2026-05-30.json`; keep science submissions paused until the exact clean-candidate build and freeze plus applicable prerequisites close |
| 2026-05-30 | Preserved live `80c0797b` after its first clean build-profile attempt rejected caller-dependent Frontier `MODULEPATH` provenance before artifact creation; staged runtime-hardening digest `c8002a1d` resets a reviewed seed before loading, publishes the exact reviewed 19-component value after module work and passes its 212-test control-plane suite, 288-test publication discovery, static checks and independent stripped-wrapper audit | Keep science submissions paused until paired successor install, active-policy promotion and exact clean-candidate freeze |
| 2026-05-30 | Activated runtime-hardening snapshot `c8002a1d`: paired immutable Orion and Project Home install, idempotent existing-genesis anchor migration, coherent policy SHA `46bac279` promotion and post-promotion validation pass while all three ledger streams remain at 22 records | Keep science submissions paused until the exact clean-candidate build and freeze plus applicable prerequisites close |
| 2026-05-30 | Preserved live `c8002a1d` after its exact clean build-profile retry rejected clean initialized `kokkos` because the authorized `/ccs/home/...` source mount alias was compared directly with its resolved `/autofs/...` path before artifact creation; staged source-alias hardening digest `6cbbbbd6` normalizes only after exact lexical authentication, preserves the second writer authorization handoff, rejects arbitrary checkout aliases including zero-submodule roots and passes its 217-test control-plane suite, 293-test publication discovery, static checks, actual-worktree probe and two independent reviews | Keep science submissions paused until paired successor install, active-policy promotion and exact clean-candidate freeze |
| 2026-05-30 | Activated source-alias hardening snapshot `6cbbbbd6`: paired immutable Orion and Project Home install, idempotent existing-genesis anchor migration, coherent policy SHA `6533172b` promotion and post-promotion validation pass while all three ledger streams remain at 22 records | Keep science submissions paused until the exact clean-candidate build and freeze plus applicable prerequisites close |
| 2026-05-30 | Passed exact clean HIP/MPI build-profile publication and immutable clean-candidate freeze `31be2cd6` through live `6cbbbbd6`: archived clean source commit `4cceb5d4`, tree `933f2e3a`, pinned `kokkos`, source bundle `e11dc8fb`, build profile `353412be`, receipt `f6a9ac20` and executable `bea2a418` | Authorize only this exact clean freeze in the reviewed policy while keeping registered science submissions paused until their applicable prerequisites close |
| 2026-05-30 | Promoted the exact-freeze authorization through installed `6cbbbbd6`: mirrored policy SHA `e8909bf5`, promotion SHA `3e7e8c7f`, idempotent anchor validation and unchanged 22-record Orion, Project Home and receipt streams pass | Preserve freeze `31be2cd6`; keep registered science submissions paused until their applicable prerequisites close |
| 2026-05-30 | Activated narrow structured clean-candidate F0 parser-contract admission policy through installed `6cbbbbd6`: mirrored policy SHA `1d525be1`, promotion SHA `2460cc0c`, exact template, deck, profile, analyzer, executable and trusted launch-contract bindings pass with unchanged 22-record streams | Submit and reconcile only the authorized F0 admission smoke; keep registered science submissions paused |
| 2026-05-30 | Preserved the first structured clean-candidate F0 retry through `6cbbbbd6` as fail-closed chronology: held Slurm job `4745523` was cancelled before execution after scheduler account spelling canonicalized to `ast207`; the durable pending marker remained until reviewed recovery | Do not guess terminal scheduler state or clear an unresolved marker after scheduler age-out; require an explicit immutable recovery boundary |
| 2026-05-30 | Activated paired immutable recovery successor `f2ad817a` after its 250-test control-plane suite, 76-test PIC publication suite, 5-test Frontier publication suite, static checks and three independent reviews pass; published mirrored handoff `06392534`, appended synthetic attachment `f518c3d8` and zero-consumption reconciliation `a13311b3`, removed the pending marker, and promoted mirrored policy SHA `2f01e58f` with promotion SHA `51c9d6a4` while all three ledger streams remain coherent at 25 records | Resume only a fresh authorized structured F0 parser-contract retry through the installed successor; keep registered science submissions paused |
| 2026-05-30 | Preserved fresh F0 compute-startup failure job `4745755` after Frontier returned `OSError 524` for the login-side Project Home writer lock before Athena; activated paired immutable successor `e8e47ead` from commit `a110e38c` with dedicated numeric scheduler-ID descriptor-pinned exact-byte read-only compute snapshots, authorized-root traversal for new and inherited parent descriptors, transient pathname-ABA rejection and retained login-side writer locks; passed its 261-test control-plane suite, 76-test PIC publication suite, 5-test Frontier publication suite, 8-test readiness registry, static checks and three independent reviews; promoted mirrored policy SHA `cfe6610a` with promotion SHA `6d6bf1cd`; completed fresh structured F0 job `4745842`, parser analysis and reconciliation while all three ledger streams remain coherent at 31 records | Close the real Frontier compute-mount startup blocker without weakening writer serialization; preserve narrow F0 authorization and keep registered-science submissions paused until applicable prerequisites close |
| 2026-05-30 | Staged registered-F1 successor `140e9a29` with exact one-attempt clean-candidate gyro and paper-coupling slices, structured immutable artifact inventories, descriptor-pinned analyzer execution and qualification-time no-write recomputation; repaired fresh independent-review findings by retaining inventory/result/receipt descriptors plus ancestry across recomputation, rejecting every extra coupling file beneath `output/`, retaining launch-time nested-directory identities, keeping the original empty owner-only `analysis/` descriptor through final publication and scoping deterministic `umask 022` creation with restoration; added independent readiness-registry recomputation of the control-plane inventory and every staged slice binding; passed fresh independent launch-publication and qualification-domain exploit retests | Keep registered-science submissions paused until clean commit curation, paired immutable install, reviewed lifecycle update and active-policy promotion pass |
| 2026-05-30 | Preserved paired immutable `140e9a29` as a never-promoted pre-activation install after an adjacent-path review found an `analysis/` mkdir/open substitution window, unbound empty structured subtrees and replaceable qualification traversal trust roots; staged repaired successor `6002c80e` with descriptor-pinned staging rename, empty-subtree rejection, trusted-PIC-root ancestry and focused regressions; passed fresh launch-publication and qualification-domain exploit retests | Keep live predecessor `e8e47ead` active and registered-science submissions paused until `6002c80e` passes clean commit curation, paired immutable install, reviewed lifecycle update and active-policy promotion |
| 2026-05-30 | Curated registered-F1 successor `6002c80e` at commit `c1be6cab`, paired-installed identical immutable inventory SHA `5144e836`, preserved existing genesis anchors without append and promoted reviewed mirrored policy SHA `ee3923fb` with active-promotion SHA `433d0afc` after empty-queue and zero-reservation validation | Execute only registered slices `f1-clean-gyro-v1` and `f1-clean-paper-coupling-v1` serially; reconcile, analyze and qualify each before broader work |
| 2026-05-30 | Reconciled registered gyro attempt `4746123` at `0.0030555555555555557` node-hours but rejected scientific result publication because required Cray MPICH display diagnostics populate stderr; staged exact reviewed-transcript SHA-256 validation in both F1 analyzers and rotated both one-attempt authorization IDs to `v2` | Preserve `4746123` as unqualified immutable chronology; promote the reviewed retry policy before submitting either `v2` slice |
| 2026-05-30 | Staged adjacent-path retry successor `cbc6fb50`: anchored launch and qualification below the site-owned serialization root, retained observed directory ancestry and workload-created directory identities through final publication, retained the pre-submit manifest through qualification recomputation, rejected regular-file namespace substitution during freeze, recorded the irreducible same-account process-isolation prerequisite from `mkdirat` through no-follow descriptor and retained-ancestry binding plus post-action descendant capture, and passed fresh launch, qualification and exact-stderr exploit retests | Preserve active `6002c80e`; keep submissions paused until clean commit curation, paired immutable install and reviewed retry-policy promotion |
| 2026-05-30 | Activated reviewed retry successor `cbc6fb50` from clean commit `559d5e54`: paired immutable Orion and Project Home inventory SHA `f5268721`, existing genesis-anchor validation without append, mirrored policy SHA `58f85503`, promotion SHA `05dfbf32` and idempotent post-promotion validation pass while all three ledger streams remain at 34 records | Archive and review the required same-account process-isolation attestation, then execute only the two `v2` registered F1 slices serially |
| 2026-05-30 | Staged late-review repair successor `d55d064d`: retained the generated immutable-inventory exclusive-creation descriptor through final publication verification; bound offline nested-directory, inventory-file and payload-file stat/open identities; retained exact inventory bytes and each published analysis-result descriptor through receipt publication after independent review reproduced accepted substitutions in active `cbc6fb50` before any v2 reservation or scheduler submission; fresh launch and offline-analysis exploit retests pass | Preserve active `cbc6fb50` and the fail-closed pre-submit manifest as immutable chronology only; keep submissions paused until clean commit curation, paired immutable install and reviewed retry-policy promotion |
| 2026-05-30 | Preserved one fail-closed `cbc6fb50` v2 pre-submit manifest as chronology after its immutable helper binding differed from the active reviewed retry-policy helper binding and reservation stopped before ledger intent or scheduler submission; separately added a changed-queue regression proving the same no-intent boundary; preserved `d55d064d` as superseded after a later review reproduced byte-identical top-level and nested frozen workload-payload namespace replacement; staged repaired successor `d8cc7d61` with retained payload identities through final inventory publication verification and direct regressions | Keep submissions paused until the exact `d8cc7d61` candidate passes fresh independent exploit retest, clean commit curation, paired immutable install and reviewed retry-policy promotion |
| 2026-05-30 | Preserved `d8cc7d61` as superseded chronology after exact review reproduced top-level payload, nested payload and generated-inventory replacement after their individual final checks while later verification continued; staged repaired successor `d538af91` with retained descriptor ownership through a closing payload namespace and byte sweep plus a final generated-inventory recheck and direct regressions | Keep submissions paused until the exact `d538af91` candidate passes fresh independent exploit retest, clean commit curation, paired immutable install and reviewed retry-policy promotion |
| 2026-05-30 | Preserved `d538af91` as superseded chronology after exact review reproduced nested-directory transplant between retained directory rebinding and the closing payload-byte sweep; staged repaired successor `4ccde8df` with a second retained directory namespace sweep after payload verification and a direct transplant regression; fresh exact exploit retest passes | Keep submissions paused until clean commit curation, paired immutable install and reviewed retry-policy promotion |
| 2026-05-30 | Curated repaired registered-F1 successor `4ccde8df` at commit `6fa926a9`, paired-installed byte-identical immutable Orion and Project Home inventory SHA `ea45cdc5`, reran live zero-reservation and no-marker validation, and preserved existing genesis anchors without ledger append | Keep submissions paused until reviewed policy promotion and the required same-account process-isolation attestations pass |
| 2026-05-30 | Promoted repaired registered-F1 successor `4ccde8df` after a fresh live zero-reservation, no-marker and rejected-UUID-absence boundary: mirrored reviewed policy SHA `3895221e`, promotion SHA `016cbb66`, idempotent genesis-anchor validation and unchanged 34-record streams pass | Archive and review the required same-account process-isolation attestation, then execute only the two `v2` registered F1 slices serially |
| 2026-05-30 | Preserved corrected-controller manifest `13f053df` after the installed wrapper rejected its operator-attestation seven-field queue snapshot before ledger intent; regenerated with the validator's exact six-field queue format, attested and reconciled gyro-v2 job `4746245` at `0.0030555555555555557` node-hours, then rejected result publication after immutable output exposed deterministic cycle-three closure and emitted particle-VTK-schema omissions in the staged gyro analyzer | Preserve both fail-closed events as unqualified chronology; bind both exact output-tree schema repairs, rotate only the consumed gyro authorization to `v3`, preserve the unused paper-v2 authorization ID with its repaired analyzer binding, promote the reviewed policy and keep submissions paused |
| 2026-05-30 | Curated parser-bound retry commits `81dba519` and `c878dd6c`, added deterministic paper-coupling linear-wave error-table closure, archived reviewed pre-promotion attestation `49d94859`, and promoted paired mirrored policy SHA `0538da6d` with promotion SHA `f84fc7b7` through installed successor `4ccde8df` while all three ledger streams remained unchanged at 37 records | Execute only attested gyro-v3 and unused paper-coupling-v2 qualification slices serially |
| 2026-05-30 | Attested gyro-v3 job `4746290` under policy SHA `0538da6d`, reconciled `0.0033333333333333335` node-hours, published passing cycle-two analytical result and receipt, replayed the snapshotted analyzer through the qualification gate and froze pending-review qualification manifest SHA `9f7a5d64` | Preserve immutable gyro-v3 evidence and execute only the serialized unused paper-coupling-v2 slice |
| 2026-05-30 | Preserved failed paper-coupling UUID `1e12dd79` after manifest creation rejected ordinary parser identifier substrings through the submission scrubber before manifest publication or ledger intent; curated scrub-safe commit `ae9fbb97`, archived pre-promotion attestation `a45a004a`, and promoted paired mirrored policy SHA `559a4b1e` with promotion SHA `84065c8a` while all three ledger streams remained unchanged at 40 records | Retain the unused paper-coupling-v2 authorization ID and retry only its attested serialized execution |
| 2026-05-30 | Attested paper-coupling-v2 job `4746297` under policy SHA `559a4b1e`, reconciled `0.006111111111111111` node-hours, published passing conservation and coefficient-invariance result plus receipt, replayed the snapshotted analyzer through the qualification gate and froze pending-review qualification manifest SHA `6262332a`; all three ledger streams remain coherent at 43 records with cumulative consumption `0.09388888888888887` node-hours, no active reservation and no pending PIC marker | Preserve registered F1 closure as pending external review evidence and continue only separately authorized remaining local and campaign gates |
| 2026-05-30 | Replayed gyro-v3 evidence through the final active policy after the paper-coupling scrub-safe promotion and froze current-policy successor qualification manifest SHA `17265e02`; preserved initial gyro freeze `9f7a5d64` as immutable chronology | Use gyro successor `17265e02` and paper-coupling manifest `6262332a` for terminal external review |
| 2026-05-30 | Linked registered-F1 gyro-v3 qualification manifest SHA `17265e02` and paper-coupling-v2 qualification manifest SHA `6262332a` into the Q-018 claim registry as immutable pending-external-review evidence | Preserve claim dispositions as open until all required child gates and named external reviews close |
| 2026-05-30 | Archived the accepted gyro-v3 and paper-coupling-v2 manifests, inventories, analytical results, receipts and terminal mirrored reconciliation records into source-local fixture index `q027_frontier_f1_accepted_closure_source_local_evidence_2026-05-30.json`; hardened the opt-in registry to validate exact mirrored ledger chains and descriptor-pinned terminal-manifest bytes against the accepted-execution projection | Preserve the source-local review bundle as a clone-local projection of immutable Orion evidence while keeping Orion as the sole bulk root |
| 2026-05-30 | Preserved rejected bounded F2 v1 job `4746310` after its launch contract omitted `-d` and allowed Athena error history to escape the registered artifact tree; archived escaped-output SHA `df1a5f2a` in Orion, repaired routing at `d59b8ebf`, rotated authorization at `2db6b060`, archived pre-promotion and pre-submit attestations, promoted mirrored policy SHA `446db26a` with promotion SHA `b772da00`, retained transient queue-snapshot retries as fail-closed pre-reservation chronology, and reconciled attested v2 job `4746316` at `0.0025` node-hours | Preserve v1 as rejected chronology and v2 inventory SHA `ca673515`, analysis SHA `95dc6e40`, receipt SHA `1685be75` and contained output SHA `df1a5f2a` as bounded F2 engineering evidence pending external review |
| 2026-05-30 | Replayed the F2 v2 snapshotted analyzer with its descriptor-pinned no-write mode against the reconciled run tree; exact inventory SHA `ca673515` and result SHA `95dc6e40` passed while the source worktree remained free of escaped `f2_multirank_runtime_metadata-errs.dat` output. All three ledger streams remain coherent at 49 records with cumulative consumption `0.09944444444444443` node-hours, no active reservation and no pending PIC marker | Continue only separately scoped local prerequisites and preregistered campaign work; do not broaden the bounded Q-002 engineering slice into a paper claim |
| 2026-05-30 | Replayed accepted gyro-v3 and paper-coupling-v2 evidence through final bounded-F2 projection policy SHA `446db26a`, froze qualification manifests SHA `d77b88f1` and `df0876df`, validated those plus bounded-F2 manifest SHA `09e71dbd`, and archived F2 terminal artifacts, attestations, ledger projections and fail-closed queue-snapshot retries in source-local index `q027_frontier_f2_accepted_closure_source_local_evidence_2026-05-30.json` | Use the historical bounded-F2 F1 projections and source-local archive for clone-local review while preserving prior F1 freezes and rejected F2 retries as immutable chronology; strict operational policy `36228012` remains live |
| 2026-05-30 | Preserved the fail-closed Q-009 x2-inflow chronology, initialized the linear-wave MHD inflow reservoir across restart, replaced leaked restart scratch buffers with stack storage, guarded particle MeshBlock offset conversion before integer casts, replayed reflecting/outflow and x2-inflow lifetime siblings under Debug and host ASan/UBSan, corrected Q-016 optional `gid`-slice wording and sampled Q-017 migration lists plus counters while live | Close bounded local broader-boundary and host-sanitizer prerequisites without promoting them to MPI, HIP, scientific-AMR, total GPU-memory or external-review qualification |
| 2026-05-30 | Hardened registered Frontier ledger lifecycle replay in strict successor `6f3458ca`, passed a retained 354-test matrix and two fresh independent adversarial reviews, paired-installed identical read-only inventories, promoted policy SHA `36228012` with promotion SHA `a2e04161`, and serialized the reviewed seven-row Q-016 direct-`srun` accounting-only import while the queue was empty. All mirrored ledger streams now contain 56 rows with cumulative consumption `0.11055555555555553` node-hours, zero active reservations and no pending PIC marker | Preserve `q027_manual_frontier_accounting_activation_2026-05-30.json` as current operational chronology; imported rows consume budget but remain ineligible for scientific evidence |
| 2026-05-30 | Replaced the Q-023 Bell source-local smoke decks' magnetic-only `mhd_bcc` output with combined `mhd_w_bcc`, retained fresh immutable 1D/2D/3D and restart-smoke Orion artifacts below `/lustre/orion/ast207/proj-shared/dfielding/PIC`, added fail-closed velocity-field extraction plus projected velocity-to-magnetic ratio diagnostics, and hardened the opt-in live registry to reject paired manual-accounting markers | Treat the new velocity ratio as a nonqualifying source-local diagnostic only; preregister the paper-literal Bell estimator, campaign tolerances and external review before any scientific qualification submission |
| 2026-05-30 | Bound the exact blocked Q-022 Sun-Bai-Zhao provenance, equation-map and empty tolerance-table placeholders into the Q-033 synthetic CRPAI transport bundle and added fail-closed checksum and semantic regressions in `q033_crpai_transport_calibration_q022_prerequisite_successor_2026-05-30.json` | Preserve the synthetic analyzer as launch-blocked nonqualifying preparation only; reference extraction, reviewed physical estimator choices, dedicated runtime generator, MPI/HIP matrix and external review remain open |
| 2026-05-30 | Replaced the Q-029 Hall-Bell source-local smoke decks' magnetic-only output with combined `mhd_w_bcc`, retained fresh immutable 1D/2D/3D Orion-local startups and exact serial 2D uninterrupted-versus-restart continuation array parity from carrier commit `731b550e`, and refreshed the Hall preregistration wording to distinguish the bound positive launch-preparation grid from the open reviewed physical grid | Preserve the Q-029 tranche as nonqualifying local mechanics evidence only; reviewed Bai mapping, physical Hall coefficients, raw Hall-dispersion analyzer, clean timestep freeze, MPI/HIP and scientific campaigns remain open |
| 2026-05-30 | Added a warnings-as-errors Q-029 compiled host contract for `CurrentDensity`, `ChiH`, the positive prepared-`chi_H` whitelist boundary and the shared Q-023 carrier reuse path | Keep the compiled contract as nonqualifying source-local implementation evidence; it does not replace the reviewed Bai mapping or Hall-Bell dispersion campaign |
| 2026-05-30 | Added full sorted-inventory digests for the immutable Q-023 and Q-029 Orion-local smoke trees, narrowed the retained Q-023 fluid-velocity statement to its cycle-zero ratio diagnostic, and archived a direct-PALS two-rank Q-032/Q-033 login-host attempt | Preserve Q-023/Q-029 inventory-bound source-local evidence; reject the Q-032/Q-033 direct-PALS outputs as policy-invalid chronology because OLCF forbids parallel login-node launches, and replay only through a registered compute-node allocation |
| 2026-05-30 | Froze the source-local Q-023 manuscript-literal Bell velocity estimator: spatial sine-fit `delta u_y` phase over fixed `pi/(k0*U_A)` intervals plus exponential growth of volume-averaged `|delta u|`, retaining fitted traces for independent phase replay | Close the locally decidable Section 5.2 estimator gap without promoting preparation fixtures into scientific qualification before clean binding, exact variants, registered execution, independent recomputation and external review |
| 2026-05-30 | Added a dedicated fail-closed Q-022 `XCMP-EXT-HALL-BELL` route, empty Bai equation-map and tolerance placeholders, a Q-029 prerequisite-binding successor and Q023 inventory entries for the existing Hall-Bell 1D/2D/3D preparation decks plus positive-`chi_H` grid materializer | Preserve the Q-029 source-local mechanics tranche without inventing a physical Bai mapping, extracted measurements, numerical thresholds, raw Hall-dispersion result or execution authorization |
| 2026-05-30 | Added a bounded Q-032 Plotnikov reduced-map applicability derivation separating the implemented `exp(-nu_in dt)` transverse source sink from the conditional high-frequency `exp(-nu_in t/2)` wave-amplitude and `exp(-nu_in t)` wave-energy envelopes, and recorded the particle phase-scrambling decision as unresolved for a matched campaign | Preserve the source-map mechanics evidence without silently promoting a conditional asymptote into a matched Plotnikov equation map, extracted dataset, tolerance table or qualification claim |
| 2026-05-30 | Repaired the Q-011 Section 5.4 preparation path to select the manuscript ideal surface `u_sh'=(Gamma-1)u0/2` and full-sphere isotropic monoenergetic injection while retaining the previous finite-Mach shock-speed estimate as an explicitly separate engineering option; added a source-local derivation regression and claim link | Close the locally decidable paper-surface and sampler mapping defect without promoting the unexecuted preparation deck into shock-campaign evidence; executed distribution audit, thermodynamic and macro-mass calibration, AMR/MPI/GPU runs, independent recompute and external review remain open |
| 2026-05-30 | Added a deterministic source-local Q-023 Bell variant materializer that stages the fixed 405-deck epsilon, dimensionality, resolution, local-CFL and PPC preparation matrix only below `tst/.codex`, refusing qualification and authorization claims while retaining centered loading as a review boundary | Advance exact Section 5.2 preparation mechanics without treating staged decks as a clean freeze, authorized run or qualifying output |
| 2026-05-30 | Added a Q-029-only raw extractor that replays the exact recursively read-only Orion Hall-Bell preparation tree and retains combined magnetic and fluid-velocity geometric seed-carrier projections with exact provenance | Preserve locally decidable raw schema and replay mechanics without inventing a Bai oracle, reviewed coefficient grid, numeric tolerance or Hall-dispersion qualification |
| 2026-05-30 | Hardened the Q-023 materializer to reserve only a direct `tst/.codex` child and create every payload through retained directory descriptors with the manifest published last; hardened Q-029 replay to reject decoder injection and special inventory entries, pin the decoder and numerical runtime, decode verified byte copies and narrow its inventory wording | Retain both paths as source-local nonqualifying preparation tools; clean-candidate binding, registered execution, physical mapping and external review remain open |
| 2026-05-30 | Added and registered the guarded Q-006 Section 5.3 source-local preparation generator with exact uniform, SMR and deterministic audited-AMR decks plus an integrated host-build pass | Preserve the ideal-EOS compatibility boundary and run retained long-horizon, true-AMR, MPI and GPU qualification separately before any paper-frequency claim |
| 2026-05-30 | Added and registered the bounded Q-033 thin-2D3V source-local runtime successor with deterministic seeded transverse modes, antipodal bounded prolate particle pairs and immutable initialized-plus-post-step extraction diagnostics | Preserve the runtime carrier as nonqualifying mechanics evidence only; physical transport calibration, Q-022 closure, MPI, Frontier HIP and external review remain open |
| 2026-05-30 | Repaired the Q-011 Section 5.4 source-local runtime contract to use upstream-relative swept mass, a unique half-open shock-surface carrier cell, modeled-surface particle placement and one-time restart-persisted removal of the startup injected cohort | Keep the frozen deck preparation-only until executed injection-distribution, thermodynamic, macro-mass, AMR, MPI, GPU, recompute and review gates close |
| 2026-05-30 | At this chronology point, added a retained Orion-backed bounded serial Q-008 exponential-profile Appendix-A CPAW four-crossing-time history and resolution-convergence preparation matrix with static mode-on/off parity and exact source bindings | Preserve the passing ten-case preparation matrix without promoting it into Appendix-A qualification; the 22-case successor below supersedes the locally executable profile, supported-solver and dimensional-carrier preparation scope |
| 2026-05-30 | At this chronology point, added a narrowly guarded exact-isothermal Q-007 `paper_mhd_pic` true-delta-f parser allowance, separately named CRSI and CRPAI cycle-zero carriers, exact source-local deck freezes and paper-literal static kappa, resonant-scale, `Q2` and low-density signed-branch mappings | Preserve source-local preparation only; the weighted-loader, four-branch-carrier and static-oracle successor below supersedes the locally executable preparation scope |
| 2026-05-30 | At this chronology point, archived a recursively read-only Orion Q-011 three-stage reduced restart-boundary diagnostic that preserves the startup cohort before cutoff, removes it while crossing cutoff and retains newly injected particles after cutoff; retained one stopped oversized-extent override only as rejected chronology | Preserve bounded nonqualifying runtime-state evidence without promoting the preparation deck into shock-campaign reproduction; the executed distribution-audit successor below supersedes the locally executable distribution item |
| 2026-05-30 | Executed the deterministic Q-023 Bell source-local materializer and recursively froze its ignored direct-child `tst/.codex` staging tree: 405 decks, 407 files, zero writable entries and exact request, manifest and inventory hashes recorded in `q023_paper_bell_linear_materialized_variants_local_2026-05-30.json` | Preserve auditable Section 5.2 source-local preparation mechanics without treating the small deck-staging tree as bulk output, a clean-candidate binding, authorized execution or qualifying evidence |
| 2026-05-30 | Executed a bounded direct serial-host Q-011 shock-injection audit and froze its Orion payload recursively read-only: actual PVTK provenance, clamped ideal-surface placement, monoenergetic surface-relative speed and bounded full-sphere statistics pass for 309 injected particles | Close the locally executable distribution-audit item without promoting the reduced carrier into physical calibration, AMR, MPI, GPU, Frontier, recomputation, review or paper-shock campaign evidence |
| 2026-05-30 | Executed and recursively froze a 22-case Orion-backed Q-008 Appendix-A CPAW extension: literal linear and reciprocal-linear convergence preparations, LLF/HLLE/HLLD sensitivity, x1/x2/x3 axis-aligned carriers and a zero-output fail-closed unsupported Roe probe pass coarse serial-host guardrails | Close locally executable profile, supported-solver and dimensional-carrier preparation gaps without promoting coarse bounds into publication thresholds, expanded dimensional convergence, driven CRPAI, MPI, GPU, Frontier or external-review qualification |
| 2026-05-30 | Added the narrowly guarded Q-006 exact-isothermal full-f momentum-only runtime successor and recursively froze bounded Orion uniform startup, uniform evolution, SMR evolution, deterministic audited-AMR evolution and audited-AMR restart snapshot recomputations | Close the locally executable Section 5.3 isothermal mechanics blocker without promoting short-horizon carriers into accepted long-horizon residual tolerances, true-AMR policy qualification, MPI, GPU, Frontier or external-review evidence |
| 2026-05-30 | Extended the Q-007 source-local true-delta-f preparation with eight-bin IPWT-weighted loading, deterministic antipodal angular sampling, a deterministic four-branch wave carrier, static full-`Q1 + Q2` oracle preparation and a pinned-executable recursively read-only two-cycle CRSI mechanics replay | Keep finite quadrature, angular sampler, seed and discrete-mode choices explicit as source-local conventions; do not promote them into recovered paper-run provenance, long-horizon growth-fit qualification, CRPAI handedness closure, MPI, GPU, Frontier or external-review evidence |
| 2026-05-30 | Rebuilt and recursively froze the final integrated serial-host provenance pin below Orion, replayed final Q-006, Q-007, Q-008 and Q-011 bounded local artifacts from that exact pin, retained raw invocations and strict root-relative-path inventories, and added direct shared-dependency bindings where the local evidence depends on shared PIC registration or extraction code | Preserve the final local mechanics evidence as nonqualifying preparation only; no Frontier submission was authorized or performed, and every scientific campaign, portability and external-review gate remains open |
| 2026-05-30 | Hardened the shared Orion verifier to reject noncanonical archive names, non-ELF pinned executables and archives unrelated to their retained dependency manifests; closed Q-006 parser rejection return codes and wrapper schemas; closed Q-007 raw startup grid metadata, wrapper schemas, empty-directory topology and direct-CLI imports; and froze the Q-008 v5 self-contained retained script, deck, reader and analysis-environment successor | Preserve the local serial-host evidence as bounded nonqualifying preparation while preventing provenance receipts, reshaped startup payloads, hidden claim keys or writable-worktree analysis dependencies from silently passing |

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
| Historical Entity source snapshot used for comparative audit | Unavailable local checkout at `a59065fc`; retain as historical provenance only |
| Frozen accessible Entity replacement snapshot | `/ccs/home/dfielding/entity` at `512998c471bf3fdec292cb4a64150c4f0aeea539`, tree `38d511d3d317d0866fab33201a7355815c26bedf`, archive SHA-256 `c92f5fac17cdbf90fc0c6adcc0b9445fa5b9f4f78c1b6a53c97aa6520423ce1d`; bounded file manifest `tst/publication/readiness/entity_snapshot.json`, especially `src/framework/containers/particles.h`, `src/kernels/pushers/sr.hpp`, `src/kernels/currents_deposit.hpp`, `tests/kernels/pusher.cpp`, `tests/kernels/deposit.cpp`, `src/kernels/digital_filter.hpp`, `src/output/stats.cpp`, `src/output/checkpoint.cpp`, `src/engines/srpic/srpic.hpp`, `src/engines/engine.hpp` and `src/framework/domain/metadomain_stats.cpp` |

This document is the entry point for all future MHD-PIC production work. A
simulation result is not a production result until its governing mode, oracle,
provenance, resource accounting and archived evidence satisfy this plan.
