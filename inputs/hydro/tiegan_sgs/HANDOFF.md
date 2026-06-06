# Project status update: Tiegan SGS turbulence and high-resolution handoff

- Date: 2026-06-06
- Exact timestamp: 2026-06-06T09:52:06-04:00
- Project or repository: AthenaK, `dfielding14/athenak`
- Report profile: standard, because this is a simulation-pilot and supercomputer handoff
- Status: validated locally but not yet production-ready at `16384 x 16384`
- Branch: `Tiegan_SGS`
- Commit: `e70ae0bd206718a5616d5558d42598b7bbe555ef` before this handoff document
- Worktree state: clean before this handoff document
- Agent identifier: Codex
- Data or simulations analyzed: completed `512 x 512` drag-0.25 and drag-0.025 pilots
- Compute environment: local macOS CPU, Release build, MPI enabled, Kokkos Serial

# Tier 0: What happened and why it matters

The project goal is to generate fully two-dimensional, compressible isothermal
turbulence data for subgrid-scale modeling. AthenaK evolves a high-resolution
density and velocity field while writing non-overlapping square-filtered density,
Favre velocity, and the three independent components of the two-dimensional SGS
stress tensor.

The implementation, inputs, tests, and pilot configurations are committed on the
`Tiegan_SGS` branch. Immediately before this handoff document was added, the clean
local branch and the live remote branch both pointed to commit
`e70ae0bd206718a5616d5558d42598b7bbe555ef`.

The simulations are genuinely two-dimensional: `nx3 = 1`, the forcing contains no
$k_z$ modes or out-of-plane force, and completed pilots retained exactly zero third
momentum and third-component kinetic energy. Periodic boundaries are correct in the
degenerate direction; reflecting boundaries are neither required nor desirable.

The forcing is a global Ornstein-Uhlenbeck process sampled from a narrow Fourier
annulus. The sparse-annulus implementation avoids enumerating the complete Cartesian
mode volume, and therefore avoids the setup cost that made large forcing wave numbers
impractical. It also avoids the repeated spatial pattern produced by tiled forcing.

Uniform linear Rayleigh friction is used as the standard large-scale sink. The
drag-0.25 `512 x 512` pilot reached a stable RMS velocity of `0.24452` and retained
only `1.31%` of its final velocity-spectrum energy in the box mode. Reducing the drag
to `0.025` at fixed injection produced a strong condensate and a final RMS velocity
of `0.45581`; `23.6%` of its final spectral energy occupied the box mode.

The current pilots do not show a resolved $k^{-3}$ forward-enstrophy range. Their
late-time velocity spectra are approximately $k^{-4.2}$ for drag `0.25` and
$k^{-3.8}$ to $k^{-4.1}$ for drag `0.025`. The likely primary limitation is that
`512 x 512` with forcing at $k_f=16$ has too little resolved forward scale separation
before WENOZ/Roe implicit dissipation becomes important. Lower drag changes the
inverse cascade dramatically but does not recover a clean $k^{-3}$ regime.

A `16384 x 16384` run is scientifically attractive because a balanced forcing choice
$k_f=\sqrt{16384}=128$ leaves substantial scale separation on both sides of the
forcing and resolves the forcing wavelength with 128 cells. It must not be launched
as an unbenchmarked one-shot production job. Sparse forcing render cost, output
volume, restart cost, MeshBlock size, and time to statistical stationarity all need
machine-specific measurements first.

The safest next action is a staged supercomputer campaign: reproduce a small MPI
case, benchmark short `2048 x 2048` and `4096 x 4096` runs, measure forcing and I/O
costs, verify restart identity, then submit a short `16384 x 16384` scaling and output
pilot. Only after those gates pass should the long production allocation begin.

# Tier 1: How the work was done

## Objective and physical model

The intended fiducial science case is low-Mach, two-dimensional, isothermal
turbulence with $c_s=1$ and target RMS Mach number $\mathcal{M}_{\rm rms}=0.1$.
Mach `0.25` has been used for faster local commissioning and drag calibration.

| Quantity | Definition or current choice |
| --- | --- |
| Domain | Periodic square with $L_{\rm box}=1$ |
| Active dimensions | `nx1 = nx2`, `nx3 = 1` |
| Equation of state | Isothermal, `iso_sound_speed = 1.0` |
| Integrator | `rk3`, CFL `0.4` |
| Reconstruction and solver | `wenoz`, `roe` |
| Driving | Solenoidal global sparse Fourier annulus |
| OU correlation time | $t_{\rm corr}=t_{\rm eddy}$ |
| Large-scale sink | Uniform Rayleigh drag, $d\boldsymbol{v}/dt=-\alpha\boldsymbol{v}$ |
| Forcing normalization | Constant energy injection rate, `normalization = edot` |
| Baseline dissipation | WENOZ/Roe implicit numerical dissipation; no explicit viscosity |

For forcing centered on box mode $k_f$,

$$
L_{\rm drive} = \frac{L_{\rm box}}{k_f},
\qquad
t_{\rm eddy} = \frac{L_{\rm drive}}{v_{\rm rms}}.
$$

The initial drag and injection calibration is

$$
\alpha = \frac{v_{\rm rms}}{L_{\rm box}},
\qquad
\dot{E} = \alpha v_{\rm rms}^2
         = \frac{v_{\rm rms}^3}{L_{\rm box}}.
$$

This balance is not a Mach-number thermostat. The achieved RMS velocity must be
measured from the history output and the drag or injection retuned if necessary.

## SGS data products

For a square, non-overlapping top-hat filter, an overbar denotes the cell average
inside one filter box and a tilde denotes Favre filtering:

$$
\widetilde{v_i} = \frac{\overline{\rho v_i}}{\overline{\rho}},
\qquad
\tau_{ij} = \overline{\rho v_i v_j}
            - \overline{\rho}\,\widetilde{v_i}\widetilde{v_j}.
$$

Each `hydro_sgs_2d` coarsened-binary output contains
`dens`, `velx`, `vely`, `tau_xx`, `tau_xy`, and `tau_yy`. Every coarsening factor
must divide both active MeshBlock dimensions. Filters are local to MeshBlocks, so
the largest usable factor is limited by the chosen MeshBlock extent.

## Completed pilots

| Run | Drag | Duration | Final RMS velocity | Final box-mode fraction | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `mach025_512_k16_mpi8_steady` | `0.25` | 40 forcing eddies | `0.24452` | `1.31%` | Stable target-Mach calibration without a box condensate |
| `mach025_512_k16_mpi8_drag0025` | `0.025` | 40 initial forcing eddies | `0.45581` | `23.6%` | Nonstationary condensate growth; not a target-Mach steady state |

Both pilots completed successfully on eight MPI ranks, produced 801 full-resolution
snapshots and all requested SGS outputs, and retained zero out-of-plane velocity.

## Proposed `16384 x 16384` production candidate

This is a candidate configuration to benchmark and refine, not a submission-ready
input.

| Parameter | Fiducial candidate | Reason |
| --- | --- | --- |
| Resolution | `16384 x 16384 x 1` | Large dual-cascade scale separation |
| MeshBlock baseline | `512 x 512 x 1` | 1024 MeshBlocks; permits filters through factor 512 |
| Alternative MeshBlock | `256 x 256 x 1` | 4096 MeshBlocks; may expose more parallelism but caps filters at 256 |
| Target Mach | `0.1` | Fiducial low-Mach SGS case |
| Sound speed | `1.0` | Then target $v_{\rm rms}=0.1$ |
| Forcing peak | `npeak = 128` | $\sqrt{N_{\rm res}}$ balanced-scale choice |
| Forcing annulus | `127 <= |k| <= 129` | Narrow annulus analogous to the validated pilots |
| Sparse complex modes | benchmark `64`, `128`, and `256` | Cost and angular isotropy are both unresolved at $k_f=128$ |
| `tcorr` | `0.078125` | One forcing-scale eddy time |
| `dt_update` baseline | `0.00078125` | `tcorr / 100`, matching pilot practice |
| Drag rate | `0.1` | Initial Mach-0.1 calibration |
| `dedt` | `0.001` | Initial Mach-0.1 calibration |
| Filter factors | `4, 8, 16, 32, 64, 128, 256, 512` with 512-square blocks | Broad SGS hierarchy; all factors divide each active block dimension |

At this forcing scale, $\alpha t_{\rm eddy}=1/128$, so uniform drag acts weakly
during one forcing-scale turnover even though it controls the large-scale energy
budget. Nevertheless, the drag remains a possible influence on direct-cascade
spectra and should be included in the scientific uncertainty budget.

## Staged high-resolution campaign

1. Reproduce the focused CPU and MPI tests on the target machine.
2. Reproduce a short committed `512 x 512` or `1024 x 1024` case and verify output
   parsing, zero $v_3$, restart continuity, and SGS stresses.
3. Benchmark `2048 x 2048` and `4096 x 4096` short runs with candidate MeshBlock
   sizes and sparse-mode counts. Separate evolution, forcing-update, restart, and
   output timings.
4. Run a short `16384 x 16384` pilot with low output cadence. Confirm memory headroom,
   load balance, forcing-update cost, I/O bandwidth, and restart viability.
5. Freeze the production input and scheduler script only after the short pilot.
6. Begin production statistics only after stationarity diagnostics pass.

# Tier 2: Detailed methods, implementation, and validation

## 2.1 Problem definition

The project needs high-resolution training targets for SGS modeling across a wide
range of filter widths. The high-resolution state supplies density and in-plane
velocity. AthenaK computes the Favre-filtered fields and exact filter-scale stresses
on the fly, avoiding the need to retain every full-resolution state solely to
construct SGS labels later.

The immediate scientific questions are:

- Can a low-Mach compressible two-dimensional flow produce a useful dual-cascade
  steady state without a box-scale condensate?
- How do the exact SGS stresses vary with filter width and flow scale?
- Does increasing resolution reveal a resolved forward-enstrophy interval that is
  absent in the `512 x 512` pilots?
- How sensitive are the results to drag, implicit dissipation, compressibility, and
  sparse forcing-mode count?

## 2.2 Data model and assumptions

The current coarsening operator is a non-overlapping local box average. It is not a
sliding convolution and it does not average across MeshBlock boundaries. Aligning
filter widths with MeshBlock dimensions makes the local boxes form a consistent
global tiling, but the filter origin remains tied to the mesh decomposition.

The `cbin` outputs use the final Favre velocities and stresses, not raw moments.
The regression test independently reconstructs the expected values from a
full-resolution conserved-state output.

The high-resolution binary output currently saves the full primitive hydro state.
The separate analysis workflow assembles all uniform MeshBlocks into global arrays
before plotting or applying NumPy FFTs. That workflow is appropriate for the local
pilots but may require a large-memory analysis node or a distributed FFT replacement
for `16384 x 16384` data.

## 2.3 Mathematical definitions

The history-file RMS velocity is

$$
v_{\rm rms}
= \sqrt{\frac{2(K_{x}+K_{y}+K_{z})}{M}},
$$

where $K_i$ is the volume-integrated kinetic energy in component $i$ and $M$ is the
total mass. A fully two-dimensional run requires $K_z=0$ for the complete run.

The local analysis scripts use shell-integrated spectra,

$$
E_f(k) =
\sum_{k-1/2 \leq |\boldsymbol{n}| < k+1/2}
|\widehat{f}(\boldsymbol{n})|^2.
$$

With this convention, a Kolmogorov inverse-energy range is proportional to
$k^{-5/3}$ and a classical forward-enstrophy range is proportional to $k^{-3}$.

## 2.4 Algorithm and implementation

The main project-specific implementation surfaces are:

| Area | Files and behavior |
| --- | --- |
| Fully 2D SGS outputs | `src/outputs/derived_variables.cpp`, `src/outputs/basetype_output.cpp`, `src/outputs/coarsened_binary.cpp`, `src/outputs/outputs.hpp` |
| Sparse annulus driving | `src/srcterms/turb_driver.cpp`, `src/srcterms/turb_driver.hpp` |
| Rayleigh drag | `src/srcterms/srcterms.cpp`, `src/srcterms/srcterms.hpp`, `src/srcterms/srcterms_newdt.cpp` |
| Restart state | `src/outputs/restart.cpp` and turbulence-driver restart metadata |
| Project inputs | `inputs/hydro/tiegan_sgs/` |
| Regression tests | `tst/test_suite/turb/` and `tst/inputs/` |

Sparse-annulus setup chooses approximately equal-angle lattice modes from one
Fourier half-plane and obtains the conjugate half-plane implicitly. This avoids
scanning a mode cube or disk up to `nhigh`.

The force rendering step still loops over all selected modes for every active cell.
Its leading cost is therefore approximately proportional to
$N_{\rm cells}N_{\rm modes}$ each time the force is refreshed. This cost has not
been benchmarked at `16384 x 16384` and is a required production gate.

## 2.5 Validation

| Purpose | Method | Expected result | Actual result | Tolerance | Status | Caveat |
| --- | --- | --- | --- | --- | --- | --- |
| Confirm branch publication | Compared `HEAD`, `origin/Tiegan_SGS`, fetched tip, and `git ls-remote` | All SHAs equal | All were `e70ae0bd206718a5616d5558d42598b7bbe555ef` before this document | Exact | passed | Must recheck after handoff commit |
| Verify SGS math and 2D contract | Ran focused CPU turbulence regressions | Direct Favre reconstruction matches; $v_3=0$ | 10 focused CPU tests passed | SGS `rtol=5e-6`, `atol=5e-8`; exact zero checks | passed | Current local Release/MPI executable |
| Verify MPI turbulence path | Ran focused MPI regression | MPI run completes and normalization matches | 1 MPI test passed | RMS forcing `rel=2e-6` | passed | Two-rank regression, not production scale |
| Verify candidate $k_f=128$ mode construction | Ran zero-cycle inputs with annulus `127` to `129` | Construct each requested sparse set | Constructed `64`, `128`, and `256` modes | Exact requested count | passed | Does not benchmark render cost or isotropy |
| Verify high-drag pilot completion | Inspected exit status, run log, and history | Reaches `t=10`, zero exit, $K_z=0$ | Reached `t=10`, exit `0`, $v_{\rm rms}=0.24452$, $K_z=0$ | Exact completion and zero $K_z$ | passed | Local eight-rank run |
| Verify low-drag comparison completion | Inspected exit status, run log, and history | Reaches `t=10`, zero exit, $K_z=0$ | Reached `t=10`, exit `0`, $v_{\rm rms}=0.45581$, $K_z=0$ | Exact completion and zero $K_z$ | passed | Not statistically stationary |
| Establish forward spectral behavior | Fit late-time shell-integrated velocity spectra | Identify resolved slope if present | High drag about `-4.1` to `-4.3`; low drag about `-3.7` to `-4.1` over tested bands | Multiple fit bands | inconclusive | A slope alone does not establish flux |
| Validate GPU production build | Not run | Target-machine GPU build and focused tests pass | No evidence available | Not available | not run | Required before GPU production |
| Validate `16384 x 16384` scaling and I/O | Not run | Stable memory, acceptable forcing and output cost | No evidence available | Machine-specific | not run | Required before production |

The CPU tests were run from the existing Release/MPI build with a temporary Python
environment containing `pytest`, NumPy, and h5py. The default local Python lacked
the complete test dependency set; that was an environment issue rather than a code
failure.

## 2.6 Results

The drag-0.25 pilot provides the current best calibrated state. Its final RMS
velocity is within about `2.2%` of the target `0.25`; its final box-mode fraction is
small; and its third velocity remains zero.

The low-drag comparison demonstrates that reducing uniform friction by a factor of
ten at fixed injection does not preserve the target Mach number. It permits rapid
large-scale accumulation and produces a strong condensate during the measured
interval.

The forward spectra do not yet support a robust classical enstrophy-cascade claim.
For the high-drag run, late-time fitted slopes are `-4.194` over `24 <= k <= 48`
and `-4.095` over `32 <= k <= 64`. For the low-drag run, the corresponding slopes
are `-3.807` and `-3.750`. The fact that both runs settle into steep tails soon
above the forcing supports insufficient resolved scale separation and broad implicit
dissipation as the leading explanation.

## 2.7 Failures and discarded approaches

- Reflecting boundaries in the degenerate direction were considered but are not
  needed. The correct fully two-dimensional contract uses `nx3 = 1`, periodic
  boundaries, no $k_z$ modes, and zero initial third momentum.
- Tiled forcing was tested as a route to higher forcing wave numbers but was rejected
  for science production because it repeats a smaller spatial realization.
- Tenfold lower drag was tested and produced a strong condensate rather than a
  comparable lower-friction steady state.
- Saving high-cadence full-resolution data at `16384 x 16384` is rejected as the
  default production plan because it creates several tebibytes of output.

## 2.8 Remaining risks

- The sparse-annulus force render is $O(N_{\rm cells}N_{\rm modes})$ per refresh.
  It may dominate at `16384 x 16384`.
- `64` modes were effective for the $k_f=16$ pilot, but angular sampling and cost
  must be reevaluated for $k_f=128$.
- Current pilots use implicit WENOZ/Roe dissipation. A resolved constant-flux
  forward-enstrophy range has not been demonstrated.
- Uniform drag is physically standard but can influence spectral slopes.
- A start-from-rest high-resolution run may require many forcing-scale eddy times
  to reach large-scale statistical stationarity because the drag time is much longer
  than $t_{\rm eddy}$ at $k_f=128$.
- Coarsening factors cannot exceed the active MeshBlock dimensions and do not cross
  MeshBlock boundaries.
- Current `cbin`, full-resolution binary, and restart I/O have not been scaled on a
  parallel filesystem.
- The local NumPy spectral workflow is not distributed and may have high peak memory
  use at `16384 x 16384`.
- No machine-specific scheduler script, GPU build recipe, or allocation estimate is
  committed because the target supercomputer has not been specified.

## 2.9 Recommended next steps

Use two distinct resolution studies rather than asking one run to answer every
question:

1. Diagnose the missing forward $k^{-3}$ range with a fixed-$k_f=16$ resolution
   ladder. Increasing resolution while holding forcing physics fixed makes movement
   of the dissipative break directly interpretable.
2. Use the balanced $k_f=\sqrt{N_{\rm res}}$ design for the eventual broad-range SGS
   production dataset, after forcing-mode and I/O benchmarks.
3. Add spectral enstrophy or potential-enstrophy flux diagnostics before making a
   cascade claim. A flux plateau is more diagnostic than a fitted slope.
4. Consider an ordinary explicit-viscosity comparison at moderate resolution to
   separate physical and implicit dissipation. Do not silently change the baseline
   ILES production model.

# Tier 3: Reproducibility, audit trail, and handoff

## 3.1 Repository state

Before this handoff document:

```text
repository: git@github.com:dfielding14/athenak.git
branch:     Tiegan_SGS
HEAD:       e70ae0bd206718a5616d5558d42598b7bbe555ef
remote:     e70ae0bd206718a5616d5558d42598b7bbe555ef
worktree:   clean
```

After cloning or pulling, use `git log -1 --oneline` to record the later commit that
adds this handoff document.

The scientific analysis scripts and generated simulation data intentionally live
outside the AthenaK branch under the local research directory
`~/Work/Research/Tiegan_SGS`. They must be transferred or versioned separately if
they are needed on the supercomputer.

## 3.2 Commands and scripts

Clone and build a CPU/MPI baseline:

```bash
git clone --branch Tiegan_SGS git@github.com:dfielding14/athenak.git
cd athenak
cmake -S . -B build-tiegan-sgs \
  -DPROBLEM=turb \
  -DAthena_ENABLE_MPI=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-tiegan-sgs -j
```

Add the target machine's Kokkos backend and architecture flags according to that
machine's supported AthenaK build recipe. Do not assume the local Kokkos Serial
configuration is appropriate for production.

Run a committed MPI commissioning input:

```bash
mpirun -np 8 build-tiegan-sgs/src/athena \
  -d RUN_DIR \
  -i inputs/hydro/tiegan_sgs/mach025_512_k16_mpi8_steady.athinput
```

On a scheduler-managed machine, replace `mpirun` with the site-supported launcher.

The focused tests used for this handoff were:

```bash
cd build-tiegan-sgs/src
PYTHONPATH=../../tst python -m pytest -q \
  ../../tst/test_suite/turb/test_turb_sgs_2d_cpu.py \
  ../../tst/test_suite/turb/test_linear_drag_cpu.py \
  ../../tst/test_suite/turb/test_turb_driving_cpu.py

PYTHONPATH=../../tst python -m pytest -q \
  ../../tst/test_suite/turb/test_turb_driving_mpicpu.py
```

The target environment needs `pytest`, NumPy, and h5py for these focused tests.

Create the first `16384 x 16384` candidate by copying the Mach-0.1 input and changing
the parameters in the proposed-candidate table. Do not submit the long run until a
short candidate input has passed the scaling gates. Keep the finalized input and
machine scheduler script together in a run-specific directory and record the exact
Git commit in that directory.

The local analysis entry points are:

```bash
python analysis/plot_sgs_slices.py --run-dir RUN_DIR
python analysis/calculate_power_spectra.py --run-dir RUN_DIR
python analysis/plot_power_spectra.py --run-dir RUN_DIR
python analysis/make_high_resolution_movie.py --run-dir RUN_DIR --overwrite
```

These scripts are in the separate local research directory, not in this branch.

## 3.3 Compute accounting

The completed local `512 x 512` runs used eight MPI ranks:

| Run | Cycles | Recorded CPU time | Zone-cycles per CPU second |
| --- | ---: | ---: | ---: |
| drag `0.25` | `21042` | `396.555 s` | `1.391e7` |
| drag `0.025` | `25246` | `492.261 s` | `1.344e7` |

These timings are not a reliable `16384 x 16384` allocation estimate. The explicit
hydrodynamic time step shrinks with resolution, sparse-force rendering has its own
scaling, and production I/O can dominate. Measure a machine-specific scaling curve
before requesting or consuming a large allocation.

## 3.4 Output inventory

The local research products are not committed:

| Local path | Description | Status | Approximate size |
| --- | --- | --- | ---: |
| `~/Work/Research/Tiegan_SGS/analysis` | Slice, movie, and spectral analysis scripts | available locally | `272 KiB` |
| `~/Work/Research/Tiegan_SGS/tiegan-sgs-m025-512-k16-mpi8-steady` | Drag-0.25 pilot, plots, spectra, movie | complete | `6.2 GiB` |
| `~/Work/Research/Tiegan_SGS/tiegan-sgs-m025-512-k16-mpi8-drag0025` | Drag-0.025 comparison, plots, spectra, movie | complete | `5.2 GiB` |

The following `16384 x 16384` estimates scale the observed `512 x 512` output layout
and assume 32-bit field output. They exclude filesystem overhead and are planning
estimates, not measured production I/O.

| Product | Estimated size |
| --- | ---: |
| One full-resolution primitive snapshot | `4.0 GiB` |
| One complete SGS set for factors `4` through `512` | about `0.50 GiB` |
| 801 full-resolution snapshots | about `3.13 TiB` |
| 801 complete SGS sets | about `0.39 TiB` |
| 41 restarts, scaled from the pilot | about `0.34 TiB` |
| High-cadence full-resolution total | about `3.86 TiB` |
| Reduced plan with 41 full-resolution snapshots and 801 SGS sets | about `0.90 TiB` |

The recommended production default is to decouple full-resolution and SGS cadence.
Retain high-cadence SGS output only if it is scientifically necessary, and save
full-resolution states much less frequently. Confirm that restart retention policy,
filesystem quotas, and postprocessing capacity agree with the final plan.

## 3.5 Known issues

- No committed `16384 x 16384` input or scheduler script exists yet.
- No GPU or target-supercomputer validation has been run.
- No high-resolution forcing-render benchmark exists.
- No parallel-filesystem I/O benchmark exists.
- No spectral flux diagnostic exists in the current analysis workflow.
- No resolved $k^{-3}$ forward range has been demonstrated.
- No initializer exists to promote a statistically steady lower-resolution state to
  a higher-resolution production mesh.

## 3.6 Continuation instructions

Before the first expensive submission:

1. Pull `Tiegan_SGS` and record the exact commit.
2. Build with the target machine's supported MPI and Kokkos configuration.
3. Run the focused CPU or GPU tests and the MPI regression.
4. Run a short committed commissioning input and verify all output readers.
5. Create a machine-specific candidate input and scheduler script.
6. Benchmark candidate MeshBlock sizes and sparse mode counts.
7. Measure output and restart time with the intended production cadence.
8. Confirm that all requested coarsening factors divide the active MeshBlock sizes.
9. Confirm exact zero third momentum and third kinetic energy.
10. Confirm restart continuity before the long run.
11. Freeze the input, scheduler script, Git commit, random seed, and output policy.
12. Start production only after memory, runtime, I/O, and stationarity gates pass.

No independent subagent review was run because delegation was not requested. The
report was directly audited against the repository state, completed pilot outputs,
focused test results, and the standard project-status report contract.
