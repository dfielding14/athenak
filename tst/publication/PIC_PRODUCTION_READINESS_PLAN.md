# AthenaK MHD-PIC Production Readiness and Sun & Bai (2023) Reproduction Plan

## Document Status

| Item | Value |
| --- | --- |
| Purpose | Canonical implementation, verification, publication-reproduction, and Frontier qualification plan for the AthenaK MHD-PIC module |
| Primary reference | Sun & Bai, *The Magnetohydrodynamic-Particle-In-Cell Module in Athena++: Implementation and Code Tests*, source manuscript in `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` |
| Review baseline | Local `c/pic-review` checkout at base commit `ca7e43f0d0d6`, reviewed with an extensive dirty working tree containing in-progress PIC changes |
| Production verdict at review time | **Not production-ready. Paper reproduction is blocked by model-level implementation gaps.** |
| Supersedes | The scope and operational content of `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md`; that narrower note remains useful historical evidence but is not sufficient for qualification |
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
   current paper-style input decks enable this field modification while leaving
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

There are additional high-risk reliability issues: particles retain a
`MeshBlockPack` pointer through AMR pack reconstruction without a visible
post-AMR pointer refresh; the load-balancing path does not yet demonstrate
particle-cost weighting; and the refinement-boundary deposition policy must be
made explicit because the paper knowingly trades individual conservation for
smooth feedback at resolution interfaces.

The correct next step is not a larger physics run. The correct next step is to
repair the model contract and establish stringent small analytical tests. Only
then should Frontier be used for GPU, AMR, scaling, and publication simulations.

## Non-Negotiable Rules For Future Agents

1. Treat this document as a living release gate. Do not silently reduce a
   tolerance, omit a failed diagnostic, or reinterpret a test as passing.
2. Prefer the paper-faithful model for reproduction. Any alternate MHD-PIC
   formulation must have a separate runtime mode, derivation, documentation,
   and independent validation; it must never be conflated with the paper mode.
3. Keep three evidence classes separate:
   - `unit/regression`: detects implementation failures;
   - `physics validation`: proves equations and algorithms against analytical
     or externally reproducible references;
   - `publication reproduction`: reruns the paper experiments and reconstructs
     quantitative claims and figures.
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
| `tst/publication/PIC_LARGE_MACHINE_VALIDATION.md` | Previously documented local validation and large-machine workflow ideas, absorbed and expanded here |
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
- AMR pointer/lifetime safety tests and a documented AMR deposition contract.
- Quantitative reproduction of every benchmark in the paper.
- Frontier HIP/MPI validation and controlled scaling results.
- Restart, decomposition, boundary, robustness, usability, and documentation
  qualification on the final model.

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
  modules, environment variables, Slurm script, input deck and analysis script.
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
| PIC-P0-001 | P0 | Coupled particle current is added directly to CT electric fields, but paper mode neglects CR Hall induction; paper-style decks also leave momentum/energy feedback off | `src/mhd/mhd_tasks.cpp` field-source and feedback paths; `inputs/tests/pic_parallel_shock_section54_hpc.athinput` and related PIC decks; paper Section 2 | Current coupled runs solve a different and generally non-conservative model | Implement an explicit paper-faithful coupling mode, forbid incompatible toggles, and pass conservation plus Bell/oscillation oracles |
| PIC-P0-002 | P0 | Pusher advances non-relativistic velocity and energy; configurable artificial light speed is rejected | `src/particles/particles_pushers.cpp`; `src/particles/particles.cpp`; paper Sections 2 and 5 | Cannot reproduce the required orbit, shock, CRSI or CRPAI physics | Implement relativistic momentum Boris with `C`; pass analytic orbit and energy/convergence suites |
| PIC-P0-003 | P0 | `pic_deltaf_mode=on` is quiet-start sampling rather than evolving delta-f weights | `src/particles/particles.cpp`, `particles_moments.cpp`, `src/particles/AGENTS.md`; paper Section 4.4 and Sections 5.5-5.7 | Published delta-f physics and low-noise instability tests are not represented | Implement weight evolution/deposition/restart; reproduce CRSI/CRPAI theory |
| PIC-P0-004 | P0 | Expanding-box code updates particle velocities only, with no corresponding MHD coordinate-source evolution | `src/particles/particles_pushers.cpp`; absence of matching MHD implementation; paper Section 3 and appendices | Driven anisotropy results do not validate the paper model | Implement full box equations; pass CPAW, gyro-motion and driven CRPAI tests |
| PIC-P0-005 | P0 | AMR rebuild reconstructs the pack while the particle module stores a pack pointer without a verified refresh | `src/mesh/mesh_refinement.cpp`; `src/particles/particles.hpp` | Potential stale-pointer use and invalid AMR scientific results | Repair lifetime/update path and pass forced refine/derefine/migration/restart tests under memory checking and GPU execution |
| PIC-P1-001 | P1 | Refinement-boundary deposition policy has not been reconciled with the paper's smoothness-versus-conservation choice | `src/particles/particles_moments.cpp`; paper Section 4.3 | AMR results may be smooth but non-reproducing, or conservative but physically different from the paper | Specify policies, validate both errors, select paper policy for reproduction |
| PIC-P1-002 | P1 | Task-stage ordering and conservative delta exchange are not yet proved against the paper's second-order method | `src/particles/particles_tasks.cpp`; MHD source paths; paper Section 2.3 | Good-looking tests may mask order loss or incorrect exchange timing | Create stage-contract tests and convergence/conservation gates |
| PIC-P1-003 | P1 | Shock setup contains a minimal directed injection construction and cannot yet satisfy the full paper setup | `src/pgen/tests/pic_parallel_shock.cpp`; shock input decks; paper Section 5.4 | Shock acceleration plots cannot yet be publication evidence | Establish exact injection distribution, units and diagnostic contract after P0 work |
| PIC-P1-004 | P1 | Sorting/intermediate-array/load-cost controls are not demonstrated as effective production mechanisms | particle runtime controls; load-balance code; paper Section 4 | Performance and load-balance claims are unsupported | Implement and benchmark, or clearly remove claims and unsupported inputs |
| PIC-P1-005 | P1 | Existing publication-named tests frequently assert sign/parity rather than analytical paper oracles | `tst/scripts/particles/pic_*publication.py` and proxy suites | Passing suite can overstate physical confidence | Reclassify proxies and add quantitative reference tests |

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
  policy and the AMR pack-lifetime fix.
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

### Acceptance Gate P0

- Clean candidate branch exists and its retained changes have code review.
- No generated run artifacts are accidentally versioned.
- Every future validation run can write a complete provenance manifest.
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

### Tests

| Test | Purpose | Required pass condition |
| --- | --- | --- |
| Zero-field ballistic motion | State and position update sanity | Exact/tight-roundoff trajectory and restart parity |
| Uniform magnetic-field gyro-motion | Paper Section 5.1 oracle | Second-order phase convergence, bounded energy error, correct `C` dependence |
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

### Tests

| Test | Scope | Required evidence |
| --- | --- | --- |
| Single particle plus uniform gas exchange | Local sign/time centering | Exact expected momentum and energy transfer |
| Many-particle periodic exchange | Global conservation | Residual convergence and roundoff floor |
| Cell-boundary and physical-boundary crossing | Shape/boundary deposition | Conservation and symmetry properties |
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

1. Repair particle-module references across `MeshBlockPack` reconstruction.
   Add a clearly owned post-AMR update method or redesign ownership so stale
   references cannot occur.
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

- AMR coupled results are permitted only after pointer/lifetime and deposition
  policy tests pass.
- Section 5.3 AMR behavior and load-balancing requirements are supported by
  archived metrics, not merely by successful completion.

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

## Phase 8: Performance, GPU Portability, Usability And Documentation

### Objective

Convert a physically correct model into an maintainable production module.

### Required Work

1. Determine whether the paper's intermediate-array and sorting strategies are
   implemented, useful and appropriate for Kokkos/HIP. Remove inert runtime
   flags or implement them with benchmark evidence.
2. Implement particle-aware load-balancing metrics and performance logging.
3. Benchmark kernel time, communication, sorting, deposition, AMR overhead,
   memory use and output overhead on representative Frontier sizes.
4. Confirm deterministic/reproducible reduction expectations; document where
   GPU/MPI order permits tolerance-level rather than bitwise agreement.
5. Provide user documentation:
   - governing equations and supported modes;
   - parameter reference with valid/invalid combinations;
   - paper-reproduction quick start;
   - Frontier runbook and accounting requirements;
   - output/analysis definitions;
   - troubleshooting and failure interpretation.
6. Make failure messages specific: unsupported mode, timestep constraint,
   missing GPU-aware MPI support, AMR policy mismatch and restart incompatibility.

### Acceptance Gate P8

- GPU runs are physically equivalent to CPU references within documented
  tolerances.
- Performance claims include raw job accounting and reproducible scripts.
- A new user can launch a supported small validation run without reading source.

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
- Cumulative Frontier node-hour ledger, demonstrating budget compliance.
- Final checklist signed in this document with links to immutable artifact paths.

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
  ledger/node_hours.csv         # mandatory accounting ledger
```

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

This project deliberately restricts all currently authorized Frontier jobs to
`debug`. Because OLCF identifies that QOS as short, non-production use and
prohibits production job chaining, do not attempt to evade the policy by
manually chaining production restart segments. Use `debug` only for compliant
short qualification/debug runs. If complete paper reproduction or a production
scaling experiment cannot be completed within the permitted use of `debug`,
stop, record the blocked work and obtain explicit revised authorization before
using any other execution policy.

### Mandatory Resource And Accounting Policy

1. Set `PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC` for every
   Frontier action. Do not run simulations in another project directory.
2. Use only `#SBATCH -q debug` and a time request no larger than `02:00:00`.
3. Submit at most one job at a time. Before each `sbatch`, run:

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
   cumulative_consumed_node_hours + maximum_node_hours > 1000
   ```

5. After the job ends, obtain actual accounting with `sacct` and append one
   ledger row. The consumed value is:

   ```text
   consumed_node_hours = allocated_nodes * elapsed_seconds / 3600
   ```

6. Keep the ledger under `$PIC_ROOT/ledger/node_hours.csv`; also copy its
   current state into each campaign manifest. The ledger columns must include:

   ```text
   timestamp,job_id,git_commit,campaign,test_id,nodes,requested_walltime,
   elapsed_seconds,node_hours,cumulative_node_hours,state,artifact_dir,notes
   ```

7. Use the budget deliberately. Begin with one-node correctness tests and only
   grow when a previous result satisfies its gate. A failed setup should cost
   no more than one small debug job before it is corrected.

### Frontier Environment File

Store the reusable environment setup as
`$PIC_ROOT/jobs/frontier_pic_environment.sh` and source it from build/run
scripts. Verify module versions when the OLCF software stack changes.

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
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
```

Keep experimental communication variables in the manifest. If a result depends
on a nondefault environment setting, run a confirming comparison before treating
it as production evidence.

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
env | sort > "${BIN_DIR}/environment.txt"

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
RUN_ID="${TEST_ID}.${SLURM_JOB_ID}"
NNODES="${SLURM_JOB_NUM_NODES}"
NRANKS=$((8 * NNODES))
NTHREADS=7
ATHENA_WALLTIME=00:18:00

source "$ENV_FILE"
export OMP_NUM_THREADS="$NTHREADS"

ATHENA="${PIC_ROOT}/bin/${GIT_COMMIT}/${CONFIG}/athena"
INPUT="${PIC_ROOT}/inputs/${CAMPAIGN}/${TEST_ID}.athinput"
OUTDIR="${PIC_ROOT}/runs/${CAMPAIGN}/${RUN_ID}"
MANIFEST_DIR="${PIC_ROOT}/manifests/${CAMPAIGN}/${RUN_ID}"

mkdir -p "$OUTDIR" "$MANIFEST_DIR" "${PIC_ROOT}/logs/slurm"
cp -p "$INPUT" "$MANIFEST_DIR/input.athinput"
cp -p "${PIC_ROOT}/jobs/${CAMPAIGN}/${TEST_ID}.sbatch" \
  "$MANIFEST_DIR/job.sbatch"
cp -p "${PIC_ROOT}/bin/${GIT_COMMIT}/${CONFIG}/athena.sha256" \
  "$MANIFEST_DIR/athena.sha256"
module -t list 2> "$MANIFEST_DIR/modules.txt"
env | sort > "$MANIFEST_DIR/environment.txt"

srun -N "$NNODES" -n "$NRANKS" -c "$NTHREADS" \
  --gpus-per-task=1 --gpu-bind=closest \
  "$ATHENA" \
  -i "$INPUT" \
  -d "$OUTDIR" \
  -t "$ATHENA_WALLTIME" \
  job/basename="$RUN_ID"

sacct -j "$SLURM_JOB_ID" \
  --format=JobIDRaw,JobName,State,ElapsedRaw,AllocNodes,AllocTRES \
  > "$MANIFEST_DIR/sacct.txt"
```

Adjust job name, walltime, node count, campaign and input only after checking
that the worst-case budget remains within the approved cap. Use shorter
walltimes for quick correctness tests.

### Submission And Ledger Procedure

Because only one debug job can exist at once, submission is an explicit serial
workflow:

```bash
#!/bin/bash
set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
JOB_SCRIPT="$1"
LEDGER="${PIC_ROOT}/ledger/node_hours.csv"

mkdir -p "${PIC_ROOT}/ledger"
touch "$LEDGER"

if squeue -u "$USER" -h -o "%q" | grep -qx "debug"; then
  printf "A debug-QOS job is already queued or running. Do not submit.\n" >&2
  exit 1
fi

printf "Review node-hour budget in %s before submitting %s\n" \
  "$LEDGER" "$JOB_SCRIPT"
sbatch "$JOB_SCRIPT"
```

After completion, manually inspect the output and append a ledger row based on
the `sacct` record. A helper script may automate parsing, but a human/agent must
still verify job state, run identity, artifact path and cumulative total before
the next submission.

### Escalating Frontier Campaign And Budget Envelope

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
| F6 Full shock reproduction and controlled scaling | All correctness gates closed | Not authorized under the current debug-only/non-production policy; revise plan and obtain approval before submission | 0 | 395.5 |
| Reserve | Unexpected defects, compliant reruns, or future explicitly authorized campaign | Must be justified in ledger/plan before use | 604.5 | 1000 |

This envelope is not permission to spend the allotment. Stop as soon as enough
evidence exists to decide a gate. If a one-node test fails, fix it before
allocating more nodes.

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
| Q-012 | Reliability | MPI/decomposition/restart/boundary and failure tests | Partial scaffolding only |
| Q-013 | Usability/docs | Supported-mode docs and runbook tested by clean launch | Open |
| Q-014 | Production sign-off | Complete artifact bundle and finding closure | Open |

## Immediate Agent Handoff: First Actions

Future implementation agents should execute the following sequence:

1. Read this plan, the paper source and all `PIC-P0-*` evidence locations.
2. Inspect the existing dirty PIC review workspace without assuming its changes
   should be merged. Preserve user work and curate changes intentionally.
3. Create or update a findings ledger and provenance manifest framework.
4. Write the paper-mode equations/interface specification and obtain review
   before changing more coupling code.
5. Implement `C`-aware relativistic momentum pushing and its analytical tests.
6. Replace or isolate the current current-to-CT field coupling and implement
   conservative paper feedback; run analytical coupling tests.
7. Implement true delta-f and full expanding-box physics with oracle tests.
8. Resolve AMR pack ownership, deposition policy and load-cost correctness.
9. Only after all small analytical gates pass, build and run Frontier validation
   under the debug-only, budget-tracked procedure above.
10. Reproduce every paper result and complete the production sign-off bundle.

## Initial Change Log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-05-25 | Created comprehensive production-readiness and paper-reproduction plan | Expanded narrow large-machine validation scope after full code/paper review exposed model-level blockers and Frontier execution requirements |

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
| AMR pack reconstruction | `src/mesh/mesh_refinement.cpp` |
| Particle shock generator | `src/pgen/tests/pic_parallel_shock.cpp` |
| Existing particle and publication tests | `inputs/tests/pic_*.athinput`, `tst/scripts/particles/pic_*.py`, `tst/publication/` |
| Paper manuscript source | `docs/reference_paper/arXiv-2304.10568v1/mnras_template.tex` |

This document is the entry point for all future MHD-PIC production work. A
simulation result is not a production result until its governing mode, oracle,
provenance, resource accounting and archived evidence satisfy this plan.
