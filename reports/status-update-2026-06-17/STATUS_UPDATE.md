# Project status update: Frontier validation handoff for stochastic scalar forcing

- Date: 2026-06-17
- Exact timestamp: 2026-06-17T09:44:43-04:00
- Project or repository: `/Users/dbf75/.codex/worktrees/scalar-ou-forcing/athenak-DF`
- Report profile: standard, because this is a multi-file scientific/HPC handoff with an intentionally unverified final tree
- Status: partially complete
- Branch: `c/scalar-ou-forcing`
- Implementation base: `3c2829697616dee3f64eae6e205df047ce030e25`; the transferred state is the tip of `c/scalar-ou-forcing`
- Transfer state: scalar-generator consolidation, Frontier integration, forcing normalization optimization, documentation, inputs, and tests are packaged together on the branch
- Agent identifier: Codex GPT-5
- Data or simulations analyzed: no simulation data analyzed for the current consolidated tree
- Compute environment: implementation prepared on macOS; final build and all post-consolidation validation are intentionally deferred to Frontier

# Tier 0: What happened and why it matters

This branch adds stochastic Ornstein--Uhlenbeck forcing for a passive scalar evolving in a frozen velocity field. The source is confined to a configurable Fourier shell, has zero mass-weighted mean, supports either prescribed source RMS or prescribed scalar-variance injection, is restartable, and exposes history diagnostics for scalar mean, variance, forcing, injection, and gradient dissipation.

The scalar problem generators were also consolidated. The canonical implementation is now `src/pgen/scalar_mixing.cpp`, derived from the newer exact-shell implementation. The redundant `scalar_mixing_old.cpp` and `scalar_mixing_perfect_powerlaw.cpp` identities were removed. Builds, inputs, tests, and documentation now target `PROBLEM=scalar_mixing`.

The canonical generator retains the features relevant to the autonomous-scalar campaign: projection, 2D stream-function, and 3D Clebsch velocity construction; selectable `exact_shell` or `statistical` spectrum contracts; configurable solenoidal fraction where supported; divergence-free scalar face velocities; flexible scalar initial conditions; OU scalar forcing hooks; and scalar forcing history diagnostics.

The scalar forcing normalization was changed for Frontier scalability. One forcing refresh now computes all five required raw moments in one Kokkos reduction and performs one `MPI_Allreduce`, followed by one device kernel that projects out the mean and applies the normalization. This replaces two device reductions, two global collectives, and separate projection/scaling kernels.

Frontier-specific build and launch paths are present. `configure_make_frontier.sh` builds a Release HIP/MPI executable for `gfx90a`. `scripts/run_frontier_scalar_ou.sbatch` requests eight MPI ranks and eight GPU GCDs per node, binds one GCD per rank, and writes binary and restart data using one file per node.

The final consolidated tree has not been built or tested. This is deliberate: the owner explicitly requested that no further testing be performed on the local machine. Earlier tests and local builds applied to a pre-consolidation state and must not be treated as validation of the current tree.

The next action is a clean Frontier build, followed by the focused validation matrix in this report. Do not begin a production campaign until the HIP build, scalar forcing normalization, restart continuity, exact-shell spectrum checks, and one-node/multi-node I/O smoke tests pass.

# Tier 1: How the work was done

The scalar equation is evolved in conservative form for $q=\rho\theta$:

$$
\partial_t q+\nabla\cdot(q\boldsymbol{u})
=\nabla\cdot(\rho\kappa\nabla\theta)+\rho f.
$$

For the intended uniform-density frozen-flow problem, this becomes

$$
\partial_t\theta+\boldsymbol{u}\cdot\nabla\theta
=\chi\nabla^2\theta+f,
$$

with $\chi=\kappa_{\rm code}$. The stochastic source is rendered from OU-evolved Fourier coefficients and projected as

$$
f=s\left(f_0-\langle f_0\rangle_\rho\right).
$$

The zero-mode projection preserves the mass-weighted scalar mean. Diffusion removes scalar variance, so statistical stationarity is diagnosed from the balance between forcing injection and scalar dissipation rather than from mean preservation alone.

Two normalizations are implemented:

- `source_rms`: impose $\langle f^2\rangle_\rho^{1/2}$.
- `variance_rate`: choose the positive scale factor that injects a configured scalar-variance increment over one timestep.

The updated normalization reduces the raw moments

$$
M,\quad \int \rho\theta\,dV,\quad \int \rho f_0\,dV,\quad
\int \rho\theta f_0\,dV,\quad \int \rho f_0^2\,dV
$$

in one device pass. The projected covariance and variance are then reconstructed algebraically. This reduces synchronization cost but introduces a numerical point that must be checked on Frontier: the projected source variance is evaluated as a difference of raw moments and may be sensitive to cancellation in extreme nonuniform-density cases.

The current production input is `inputs/hydro/scalar_mixing_stochastic_forcing.athinput`. It uses a $128^2$ periodic scalar-only problem, statistical frozen-velocity spectrum, low-wave-number OU scalar forcing, source-RMS normalization, history output, per-node binary output, and per-node restart output.

The canonical Frontier workflow is:

1. build with `configure_make_frontier.sh`;
2. validate the executable and focused tests in a GPU allocation;
3. run a short one-node smoke test through `scripts/run_frontier_scalar_ou.sbatch`;
4. restart from its per-node checkpoint;
5. repeat the smoke test on at least two nodes;
6. only then stage longer statistical-steady-state runs.

# Tier 2: Detailed methods, implementation, and validation

## 2.1 Problem definition

The objective is a statistically stationary passive-scalar cascade in a prescribed autonomous velocity field. The source must not drive the scalar mean, and the forcing/diffusion budget must be measurable well enough to diagnose steady state, anomalous dissipation, and Yaglom-law scaling.

## 2.2 Data model and assumptions

- AthenaK stores passive scalar density $q=\rho\theta$.
- The intended physical interpretation assumes uniform frozen density.
- The velocity field is generated once and then held fixed by `scalar_only=true`.
- Scalar forcing modes are global Fourier modes; their complex OU coefficients and RNG state are authoritative restart data.
- History data are additive volume integrals and are reduced by the existing history-output machinery.
- The primary Frontier mapping is eight MPI ranks per node and one MI250X GCD per rank.

## 2.3 Mathematical definitions

The mass-weighted mean and scalar variance are

$$
\bar{\theta}=\frac{\int\rho\theta\,dV}{\int\rho\,dV},
\qquad
Z=\frac{1}{2}\left\langle(\theta-\bar{\theta})^2\right\rangle_\rho.
$$

The stationary budget is

$$
\frac{dZ}{dt}=\mathcal I-\varepsilon_\theta,
\qquad
\mathcal I=\left\langle(\theta-\bar{\theta})f\right\rangle_\rho,
\qquad
\varepsilon_\theta=\chi\langle|\nabla\theta|^2\rangle_\rho.
$$

The history columns needed to reconstruct these terms are documented in `docs/stochastic_scalar_forcing.tex`.

## 2.4 Algorithm and implementation

Relevant implementation surfaces:

- `src/srcterms/scalar_driver.cpp` and `.hpp`
  - OU coefficient update;
  - modal rendering;
  - zero-mean projection;
  - `source_rms` and `variance_rate` normalization;
  - Runge--Kutta stage application;
  - AMR-aware basis rebuilding;
  - restart metadata.
- `src/pgen/scalar_mixing.cpp`
  - canonical frozen-velocity/scalar problem generator;
  - exact-shell or statistical velocity spectrum;
  - divergence-free scalar face velocity;
  - scalar initialization and legacy scalar source modes;
  - stochastic forcing history diagnostics.
- `src/outputs/restart.cpp` and `src/pgen/pgen.cpp`
  - scalar OU restart write/read path.
- `src/outputs/basetype_output.cpp`
  - `variable=scalar_force` output.
- `configure_make_frontier.sh`
  - Release, MPI, HIP, ZEN3, and VEGA90A build.
- `scripts/run_frontier_scalar_ou.sbatch`
  - Frontier task/GPU binding and run-directory setup.

The generator consolidation intentionally removed the legacy JSON generator-diagnostics implementation and its three dedicated post-processing scripts. The retained exact-shell regression reads the actual solver output field, which is the stronger validation target for the current campaign.

## 2.5 Validation

All rows below refer to the consolidated tree transferred on `c/scalar-ou-forcing`.

| Purpose | Method | Expected result | Actual result | Tolerance | Status | Caveat |
| --- | --- | --- | --- | --- | --- | --- |
| HIP/MPI compilation | Run `./configure_make_frontier.sh` on Frontier | Release executable at `build_scalar_mixing/src/athena` | Not run | Clean build | not run | Highest-priority gate |
| Canonical generator resolution | Build with `-DPROBLEM=scalar_mixing` | CMake compiles `src/pgen/scalar_mixing.cpp` | Not run | Exact target match | not run | Rename occurred after prior local builds |
| Source mean projection | Run scalar source-RMS regression | `rf_s0/rho` remains near zero | Not run | Existing test uses absolute $2\times10^{-15}$ | not run | HIP/MPI reduction order may expose a larger residual |
| Source-RMS normalization | Run scalar source-RMS regression | RMS equals 0.2 | Not run | Existing test uses relative $2\times10^{-14}$ | not run | Do not relax before recording observed error |
| Variance-rate normalization | Run variance-rate regression | Measured rate equals 0.04 | Not run | Existing test uses relative $2\times10^{-13}$ | not run | Tests the raw-moment algebra |
| Mean preservation | Run both normalization tests | Scalar mean is unchanged | Not run | Existing test uses absolute $2\times10^{-15}$ | not run | Required for zero-mode correctness |
| Restart continuity | Split and resume a run | Scalar and rendered force are bitwise identical | Not run | Bitwise equality | not run | Must test shared and per-node restart paths |
| MPI global normalization | Run at least two ranks | Global source mean is zero and RMS target is met | Not run | Same as serial test initially | not run | Existing pytest hardcodes `mpirun`; Frontier should use `srun` or a launcher override |
| Exact-shell velocity spectrum | Run `hydro/scalar_mixing_spectrum` regression | Shell error, divergence, leakage, RMS, and mean pass | Not run | Stream shell $5\times10^{-7}$; projection shell $2\times10^{-6}$; divergence $2\times10^{-7}$; leakage $5\times10^{-8}$ | not run | Validates solver-effective velocity, not just coefficients |
| One-node production smoke | Submit short batch run | Eight ranks/GCDs complete and outputs are readable | Not run | No hangs/errors; expected files present | not run | Also validates module tuple and binding |
| Per-node restart smoke | Restart the one-node checkpoint | Continued run completes with preserved OU state | Not run | Exact continuation where compared | not run | Critical for long Frontier campaigns |
| Multi-node collective/I/O smoke | Repeat on two nodes | No collective mismatch or I/O failure | Not run | Completion with sane history | not run | Exercises one `MPI_Allreduce` per forcing refresh and per-node shards |

## 2.6 Results

There are no Frontier results yet. The robust implementation conclusions are limited to source inspection:

- the current source has one canonical scalar-mixing generator;
- the forcing normalization contains one five-moment device reduction and one MPI global reduction per refresh;
- the Frontier build and batch scripts target `scalar_mixing`;
- the production input requests per-node binary and restart output;
- focused regression tests exist for forcing normalization, restart continuity, MPI projection, and exact-shell velocity spectra.

## 2.7 Failures and discarded approaches

- The repository previously exposed `scalar_mixing_old` and `scalar_mixing_perfect_powerlaw`, while many build commands still requested a missing `scalar_mixing.cpp`. This ambiguity was removed.
- The legacy JSON generator-diagnostics path was not transplanted into the canonical generator. Its dedicated analysis scripts were removed because they depended on outputs the retained generator did not produce.
- Local post-consolidation testing was explicitly not performed.
- Earlier local builds were made before the generator rename and are superseded. Their build directories were removed and must not be used as evidence.

## 2.8 Remaining risks

1. The canonical generator rename has not been compiled.
2. Five simultaneous Kokkos reducers have not been compiled with the Frontier HIP toolchain in this tree.
3. The raw-moment variance formula may expose cancellation at extreme density contrast or pathological AMR weighting.
4. The very tight CPU-oriented regression tolerances may be too strict for a different parallel reduction order on HIP. Any adjustment must be based on measured residuals and budget closure, not convenience.
5. `test_scalar_driving_mpicpu.py` invokes `mpirun` directly. Frontier validation should use `srun`; either run an equivalent command manually or make the launcher configurable before using that test.
6. The default module tuple must be checked against the currently available Frontier compatibility matrix before building.
7. `RUN_DIR` should be set explicitly if the local `$MEMBERWORK/$SLURM_JOB_ACCOUNT` layout differs from the batch script assumption.
8. The Frontier agent must fetch `origin/c/scalar-ou-forcing` and record its exact tip before building so all validation evidence is tied to one transferred tree.

## 2.9 Recommended next steps

1. Fetch and check out `origin/c/scalar-ou-forcing` without changing scientific parameters.
2. Build on Frontier using the supplied build script.
3. Run the forcing normalization and restart tests in a one-GPU allocation.
4. Adapt the MPI regression launcher to `srun` and run it on two GCDs.
5. Run the exact-shell velocity regression.
6. Run one-node and two-node production smoke tests with aggressive output/restart cadence.
7. Record all logs, module versions, job IDs, history files, and restart manifests.
8. Only after all gates pass, begin a steady-state forcing run and verify $\langle\mathcal I\rangle_t\simeq\langle\varepsilon_\theta\rangle_t$.

# Tier 3: Reproducibility, audit trail, and handoff

## 3.1 Repository state

At report time:

```text
branch: c/scalar-ou-forcing
implementation base: 3c2829697616dee3f64eae6e205df047ce030e25
tracking after transfer: origin/c/scalar-ou-forcing
commits since origin/autonomous_scalar: 5
state: all changes described by this report are included in the branch tip
```

Major current changes:

- canonicalized `src/pgen/scalar_mixing.cpp`;
- deleted `src/pgen/scalar_mixing_old.cpp`;
- deleted `src/pgen/scalar_mixing_perfect_powerlaw.cpp`;
- optimized `src/srcterms/scalar_driver.cpp`;
- canonicalized `configure_make_frontier.sh`;
- added `scripts/run_frontier_scalar_ou.sbatch`;
- updated scalar forcing and exact-shell inputs/tests;
- updated the scalar forcing and scalar mixing documentation;
- removed legacy JSON-diagnostics post-processing scripts.

## 3.2 Commands and scripts

### Frontier module compatibility check

Before building, compare the default tuple in `configure_make_frontier.sh` with the current OLCF documentation:

- [Frontier user guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)
- [Frontier job scripts and binding examples](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#running-jobs)

The current defaults are:

```text
cpe/25.09
cray-mpich/9.0.1
rocm/6.4.2
cce/20.0.0
craype-accel-amd-gfx90a
```

Override all four versioned modules together if OLCF requires another compatible tuple.

### Build

```bash
cd /path/to/athenak-DF
./configure_make_frontier.sh
```

Expected executable:

```text
build_scalar_mixing/src/athena
```

### Configuration check on a compute node

```bash
salloc -A PROJECT_ID -p batch -N 1 -t 00:30:00
srun -N 1 -n 1 -c 7 \
  --gpus-per-task=1 --gpu-bind=closest \
  ./build_scalar_mixing/src/athena -c
```

Confirm:

```text
Problem generator: scalar_mixing
MPI parallelism: ON
```

### Scalar forcing tests

Use a Frontier Python environment containing `pytest` and `numpy`. Run the single-process test inside a one-GPU allocation:

```bash
srun -N 1 -n 1 -c 7 \
  --gpus-per-task=1 --gpu-bind=closest \
  bash -lc '
    cd /path/to/athenak-DF/build_scalar_mixing/src
    python3 -m pytest -q ../../tst/test_suite/scalar/test_scalar_driving_cpu.py
  '
```

For the MPI test, do not rely on the hardcoded `mpirun` path without checking Frontier behavior. Prefer a small patch that makes the launcher configurable, then use an `srun` launcher with two ranks and one GCD per rank.

### Exact-shell velocity regression

The retained regression entry point is:

```text
tst/scripts/hydro/scalar_mixing_spectrum.py
```

Run it through `tst/run_tests.py` with the same HIP/MPI CMake options as the production build. Preserve the complete build and test log.

### One-node smoke test

```bash
cd /path/to/athenak-DF
RUN_DIR="$MEMBERWORK/PROJECT_ID/scalar_ou/smoke-$(date +%Y%m%d-%H%M%S)" \
sbatch -A PROJECT_ID -N 1 -t 00:30:00 \
  scripts/run_frontier_scalar_ou.sbatch \
  time/nlim=6 time/tlim=1.0 \
  output1/dt=1.0e-6 output2/dt=1.0e-6 \
  output3/dt=1.0e-6 output4/dt=1.0e-6
```

Inspect the Slurm log, history file, per-node binary shards, restart manifest, and node payload shard.

### Two-node smoke test

Repeat the command with `-N 2`. The script will launch 16 ranks and bind one GCD per rank.

## 3.3 Compute accounting

- Frontier jobs used so far: none
- Frontier node-hours used so far: 0
- Recommended initial allocation: one node for 30--60 minutes
- Recommended multi-node smoke: two nodes for 15--30 minutes
- Estimates above are planning values, not measured runtimes

## 3.4 Output inventory

| Output path | Description | Status | Size | Needed for future work | Regenerable |
| --- | --- | --- | --- | --- | --- |
| `reports/status-update-2026-06-17/STATUS_UPDATE.md` | This Frontier handoff | available | small | yes | yes |
| `docs/stochastic_scalar_forcing.tex` | Mathematical setup, diagnostics, anomalous dissipation, Yaglom law, Frontier notes | included in transfer | small | yes | yes |
| `inputs/hydro/scalar_mixing_stochastic_forcing.athinput` | Canonical OU forcing production input | included in transfer | small | yes | yes |
| `configure_make_frontier.sh` | Canonical Frontier build script | included in transfer | small | yes | yes |
| `scripts/run_frontier_scalar_ou.sbatch` | Frontier launch script | included in transfer | small | yes | yes |
| Frontier build log | HIP/MPI compiler evidence | not available | not available | yes | yes |
| Frontier test logs | Normalization/restart/spectrum evidence | not available | not available | yes | yes |
| Frontier smoke outputs | History, binary shards, restart manifest/payload | not available | not available | yes | yes |

## 3.5 Known issues

- No final-tree build exists.
- No Frontier job has run.
- The MPI pytest launcher is not Frontier-native.
- Legacy JSON initialization diagnostics were removed; downstream users of those deleted scripts must migrate to solver-output-based diagnostics.
- The retained figure directory `docs/figures/scalar_mixing_perfect_powerlaw/` keeps its historical name because it contains an existing binary artifact; this does not indicate a second active generator.

## 3.6 Continuation instructions

The Frontier agent should:

1. fetch `origin/c/scalar-ou-forcing` and record `git rev-parse HEAD`;
2. record `git status --short`, branch, commit, module list, compiler version, ROCm version, and Cray MPICH version;
3. perform no cleanup or refactor before the first build;
4. build with `configure_make_frontier.sh`;
5. stop immediately on a compile failure and report the first complete compiler error;
6. run the validation matrix in the listed order;
7. preserve failed outputs and logs;
8. do not loosen numerical tolerances until the measured residual and expected HIP reduction-order effect are documented;
9. report job IDs and output paths in a follow-up status update;
10. declare production readiness only after restart and two-node smoke tests pass.

## 3.7 Document quality assurance

The standard-profile status-report audit completed with zero warnings. The report was then reviewed directly against the current repository status. Independent subagent review was not performed because delegation was not requested for this turn.
