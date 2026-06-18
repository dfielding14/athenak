# Project status update: Tiegan SGS turbulence and high-resolution handoff

- Date: 2026-06-18
- Exact timestamp: 2026-06-18T10:16:48-04:00
- Project or repository: AthenaK, `dfielding14/athenak`
- Report profile: standard, because this is a simulation-pilot and supercomputer handoff
- Status: explicit-viscosity production inputs configured; supercomputer convergence not run
- Branch: `Tiegan_SGS`
- Commit: `a03d55141b36ddc25e822f0cb33465d89e208a40` before the scale-separation retuning
- Worktree state: modified by the scale-separation input and test update
- Agent identifier: Codex
- Data or simulations analyzed: completed `512 x 512` drag-0.25 and drag-0.025 pilots
- Compute environment: local macOS CPU, Release build, MPI enabled, Kokkos Serial

# Tier 0: What happened and why it matters

The project goal is to generate fully two-dimensional, compressible isothermal
turbulence data for subgrid-scale modeling. AthenaK evolves a high-resolution
density and velocity field while writing non-overlapping square-filtered density,
Favre velocity, and the three independent components of the two-dimensional SGS
stress tensor.

The implementation, inputs, tests, and pilot configurations live on the `Tiegan_SGS`
branch. Immediately before the scale-separation retuning, the clean local branch and
the live remote branch both pointed to commit
`a03d55141b36ddc25e822f0cb33465d89e208a40`.

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

The science configuration now requires ordinary explicit viscosity. The same physical
coefficient, $\nu=8.4\times10^{-7}$, is used in committed `12288 x 12288` and
`16384 x 16384` inputs. With Mach `0.1`, $k_f=16$, and energy injection `0.001`, the
estimated viscous length is sampled by `7.66` and `10.21` cells, respectively, while
the estimated dissipation mode remains a factor `15.96` above the forcing peak.

Resolved viscosity changes the forcing-scale argument. Keeping $k_f=128$ would put
a conservatively resolved viscous cutoff too close to the forcing to leave a useful
direct-cascade interval. The production pair therefore uses $k_f=16$, close to the
logarithmic midpoint between the box mode and the estimated viscous cutoff.

Roe remains the baseline Riemann solver. AthenaK's HLLC implementation is ideal-gas
only and explicitly rejects an isothermal EOS. HLLE supports isothermal hydro but is
more diffusive, so it is reserved as a robustness fallback rather than mixed into a
calculation intended to have explicit viscosity control the small-scale dissipation.

The safest next action is a staged supercomputer campaign: reproduce the focused MPI
case, run a short downscaled parse/output check, benchmark the committed `12288 x 12288`
input, and then submit a short `16384 x 16384` scaling and restart pilot. The long
production run begins only after the same-physics spectra and dissipation budgets show
that the `16384 x 16384` cutoff is converged and explicit-viscosity dominated.

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
| Baseline dissipation | Uniform isotropic kinematic shear viscosity, $\nu=8.4\times10^{-7}$ |

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

For a low-Mach solenoidal forcing band centered on physical wavenumber
$q_f=2\pi k_f/L_{\rm box}$, the design estimate is

$$
\eta_f \simeq q_f^2\dot{E},
\qquad
\ell_\nu = \left(\frac{\nu^3}{\eta_f}\right)^{1/6},
\qquad
k_\nu = \frac{L_{\rm box}}{2\pi\ell_\nu}.
$$

For the production parameters, $\eta_f\simeq10.1065$,
$\ell_\nu\simeq6.2331\times10^{-4}$, and $k_\nu\simeq255.3$. These are design
estimates. The actual run must measure the enstrophy injection and dissipation budget.

As a reference point rather than a proof, the `16384^2` explicit-viscosity calculation
reported by Bernard et al. used $\nu=10^{-6}$ and quoted a dissipation length about
nine grid cells wide. The present candidate has comparable sampling at `10.21` cells.
See [Bernard et al. (2006)](https://arxiv.org/abs/nlin/0602017).

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
They had no explicit viscosity and remain historical ILES commissioning evidence.

## Proposed `16384 x 16384` production candidate

This committed candidate must pass target-machine scaling and the `12288^2` convergence
comparison before it is treated as production-ready.

| Parameter | Fiducial candidate | Reason |
| --- | --- | --- |
| Resolution | `16384 x 16384 x 1` | Large dual-cascade scale separation |
| MeshBlock baseline | `512 x 512 x 1` | 1024 MeshBlocks; permits filters through factor 512 |
| Alternative MeshBlock | `256 x 256 x 1` | 4096 MeshBlocks; may expose more parallelism but caps filters at 256 |
| Target Mach | `0.1` | Fiducial low-Mach SGS case |
| Sound speed | `1.0` | Then target $v_{\rm rms}=0.1$ |
| Explicit viscosity | `8.4e-7` | Same physical coefficient at `12288^2` and `16384^2` |
| Forcing peak | `npeak = 16` | Near the log midpoint of box and resolved viscous cutoff |
| Forcing annulus | `15 <= |k| <= 17` | Validated narrow global annulus |
| Sparse complex modes | `64` | Validated angular sampling at $k_f=16$ |
| `tcorr` | `0.625` | One forcing-scale eddy time |
| `dt_update` baseline | `0.00625` | `tcorr / 100`, matching pilot practice |
| Drag rate | `0.1` | Initial Mach-0.1 calibration |
| `dedt` | `0.001` | Initial Mach-0.1 calibration |
| Estimated $\ell_\nu/\Delta x$ | `10.21` | Resolved viscous-scale sampling |
| Estimated $k_\nu/k_f$ | `15.96` | Target direct-cascade scale separation of about 16 |
| Filter factors | `4, 8, 16, 32, 64, 128, 256, 512` with 512-square blocks | Broad SGS hierarchy; all factors divide each active block dimension |

At this forcing scale, $\alpha t_{\rm eddy}=1/16$, so uniform drag acts weakly
during one forcing-scale turnover even though it controls the large-scale energy
budget. Nevertheless, the drag remains a possible influence on direct-cascade
spectra and should be included in the scientific uncertainty budget.

## Staged high-resolution campaign

1. Reproduce the focused CPU and MPI tests on the target machine.
2. Run the committed `12288 x 12288` input and verify output parsing, zero $v_3$,
   restart continuity, SGS stresses, and the viscous dissipation budget.
3. Benchmark `512 x 512` and `256 x 256` MeshBlocks. Separate evolution, forcing,
   restart, and output timings.
4. Run a short `16384 x 16384` pilot. Confirm memory headroom, load balance, forcing
   cost, I/O bandwidth, and restart viability.
5. Compare `12288 x 12288` and `16384 x 16384` spectra, fluxes, and integrated viscous
   dissipation over their common resolved range.
6. Begin production statistics only after convergence and stationarity diagnostics pass.

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
- Does the same-physics `12288^2`/`16384^2` pair demonstrate an explicitly viscous,
  converged forward-enstrophy interval absent from the `512 x 512` ILES pilots?
- How sensitive are the results to drag, compressibility, and residual numerical
  dissipation above the physical viscous cutoff?

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

AthenaK implements the requested dissipation as a uniform, isotropic Newtonian shear
viscosity added directly to the hydrodynamic fluxes. In an isothermal calculation,
the removed kinetic energy is not retained as thermal energy, consistent with the
assumed instantaneous isothermal cooling.

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
| Explicit shear viscosity | `src/diffusion/viscosity.cpp`, enabled by `<hydro>/viscosity` |
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

HLLC is not an alternative for these inputs: AthenaK rejects HLLC with an isothermal
EOS. Roe supports isothermal hydro and is less diffusive than HLLE. Since the explicit
viscous cutoff is predicted near mode `255`, far below the numerical grid cutoff,
Roe is retained so numerical diffusion is minimized rather than deliberately increased.

## 2.5 Validation

| Purpose | Method | Expected result | Actual result | Tolerance | Status | Caveat |
| --- | --- | --- | --- | --- | --- | --- |
| Confirm branch baseline | Compared `HEAD`, `origin/Tiegan_SGS`, and the live remote before editing | All SHAs equal | All were `a03d55141b36ddc25e822f0cb33465d89e208a40` | Exact | passed | Must recheck after retuning commit |
| Verify SGS math, viscosity path, and 2D contract | Ran focused CPU turbulence regressions with isothermal viscosity enabled | Direct Favre reconstruction matches; $v_3=0$ | 13 focused CPU tests passed | SGS `rtol=5e-6`, `atol=5e-8`; exact zero checks | passed | Current local Release/MPI executable |
| Verify MPI turbulence path | Ran forcing and four-rank viscous-SGS regressions | MPI normalization matches; viscous SGS output retains $v_3=0$ | 2 MPI tests passed | RMS forcing `rel=2e-6`; exact zero checks | passed | Local regression scale, not production scale |
| Verify production input execution | Ran both inputs at downscaled `512^2` with zero cycles and the `16384^2` input for two cycles | Inputs initialize and the explicit-viscosity task list advances | Both initialized; the smoke run reached cycle 2 | Exact successful execution | passed | Not a performance or physics run |
| Verify viscosity resolution contract | Parsed both production inputs and evaluated the design estimate | Same $\nu$; at least 7.5 and 10 cells per $\ell_\nu$; $k_\nu/k_f\simeq16$ | `7.66` and `10.21` cells; $k_\nu/k_f=15.96$ | Configuration assertions | passed | Uses estimated, not measured, enstrophy injection |
| Verify solver compatibility | Attempted isothermal HLLC initialization | Unsupported combination is rejected | Fatal rejection matched the expected message | Exact message | passed | Roe remains the baseline |
| Verify high-drag pilot completion | Inspected exit status, run log, and history | Reaches `t=10`, zero exit, $K_z=0$ | Reached `t=10`, exit `0`, $v_{\rm rms}=0.24452$, $K_z=0$ | Exact completion and zero $K_z$ | passed | Local eight-rank run |
| Verify low-drag comparison completion | Inspected exit status, run log, and history | Reaches `t=10`, zero exit, $K_z=0$ | Reached `t=10`, exit `0`, $v_{\rm rms}=0.45581$, $K_z=0$ | Exact completion and zero $K_z$ | passed | Not statistically stationary |
| Establish forward spectral behavior | Fit late-time shell-integrated velocity spectra | Identify resolved slope if present | High drag about `-4.1` to `-4.3`; low drag about `-3.7` to `-4.1` over tested bands | Multiple fit bands | inconclusive | A slope alone does not establish flux |
| Validate GPU production build | Not run | Target-machine GPU build and focused tests pass | No evidence available | Not available | not run | Required before GPU production |
| Validate `16384 x 16384` scaling and I/O | Not run | Stable memory, acceptable forcing and output cost | No evidence available | Machine-specific | not run | Required before production |

The CPU tests were run from the existing Release/MPI build with the persistent local
analysis Python environment containing `pytest`, NumPy, and h5py.

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

The explicit-viscosity production design is a calculated configuration, not a new
simulation result. Its predicted viscous length is `6.2331e-4`, corresponding to
`7.66` cells at `12288^2` and `10.21` cells at `16384^2`. The forcing-scale Reynolds
number is approximately `7440`, so the selected viscosity remains small at the
driving scale while being deliberately resolvable at the dissipation scale.

## 2.7 Failures and discarded approaches

- Reflecting boundaries in the degenerate direction were considered but are not
  needed. The correct fully two-dimensional contract uses `nx3 = 1`, periodic
  boundaries, no $k_z$ modes, and zero initial third momentum.
- Tiled forcing was tested as a route to higher forcing wave numbers but was rejected
  for science production because it repeats a smaller spatial realization.
- Tenfold lower drag was tested and produced a strong condensate rather than a
  comparable lower-friction steady state.
- The previous $k_f=128$ `16384^2` proposal was discarded after imposing the
  resolved-viscosity requirement. A conservative viscous cutoff would leave too
  little scale separation above that forcing band.
- The `8192^2` convergence input was discarded after retuning to $k_\nu/k_f\simeq16$;
  it would sample the estimated viscous length with only `5.11` cells.
- HLLC was considered but is unavailable for AthenaK's isothermal hydro equations.
- Saving high-cadence full-resolution data at `16384 x 16384` is rejected as the
  default production plan because it creates several tebibytes of output.

## 2.8 Remaining risks

- The sparse-annulus force render is $O(N_{\rm cells}N_{\rm modes})$ per refresh.
  It may dominate at `16384 x 16384`.
- The completed pilots use implicit WENOZ/Roe dissipation. The newly committed
  explicit-viscosity pair has not yet been run at full resolution.
- The estimated $\ell_\nu$ uses the low-Mach relation
  $\eta_f\simeq(2\pi k_f/L)^2\dot E$; the measured enstrophy budget may shift it.
- AthenaK's Roe solver has no dedicated all-speed low-Mach correction. WENOZ reduces
  interface jumps in smooth flow, but convergence must still show that Roe's residual
  numerical dissipation is subdominant to the explicit viscous budget.
- A resolved constant-flux forward-enstrophy range has not been demonstrated.
- Uniform drag is physically standard but can influence spectral slopes.
- A start-from-rest high-resolution run may require many forcing-scale eddy times
  to reach large-scale statistical stationarity because the drag time is much longer
  than $t_{\rm eddy}$ at $k_f=16$.
- Coarsening factors cannot exceed the active MeshBlock dimensions and do not cross
  MeshBlock boundaries.
- Current `cbin`, full-resolution binary, and restart I/O have not been scaled on a
  parallel filesystem.
- The local NumPy spectral workflow is not distributed and may have high peak memory
  use at `16384 x 16384`.
- No machine-specific scheduler script, GPU build recipe, or allocation estimate is
  committed because the target supercomputer has not been specified.

## 2.9 Recommended next steps

1. Run the committed `12288^2` case with $\nu=8.4\times10^{-7}$ and measure the actual
   energy and enstrophy injection and viscous dissipation rates.
2. Run the same physical model at `16384^2`. Require agreement over the shared
   resolved range and confirm that the viscous rolloff does not move materially.
3. Add spectral enstrophy or potential-enstrophy flux diagnostics before making a
   cascade claim. A flux plateau is more diagnostic than a fitted slope.
4. Treat the older inviscid inputs as pipeline and historical comparison cases only.

# Tier 3: Reproducibility, audit trail, and handoff

## 3.1 Repository state

Before the scale-separation retuning:

```text
repository: git@github.com:dfielding14/athenak.git
branch:     Tiegan_SGS
HEAD:       a03d55141b36ddc25e822f0cb33465d89e208a40
remote:     a03d55141b36ddc25e822f0cb33465d89e208a40
worktree:   clean
```

After cloning or pulling, use `git log -1 --oneline` to record the later commit that
retunes the explicit-viscosity production pair.

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

Run the lower-resolution member of the explicit-viscosity pair first:

```bash
mpirun -np RANKS build-tiegan-sgs/src/athena \
  -d RUN_DIR \
  -i inputs/hydro/tiegan_sgs/mach010_12288_k16_viscous.athinput
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
  ../../tst/test_suite/turb/test_turb_sgs_2d_mpicpu.py \
  ../../tst/test_suite/turb/test_turb_driving_mpicpu.py
```

The target environment needs `pytest`, NumPy, and h5py for these focused tests.

The matching high-resolution input is
`inputs/hydro/tiegan_sgs/mach010_16384_k16_viscous.athinput`. Do not submit its full
40-turnover duration until a short candidate job has passed the scaling, restart,
and output gates. Keep the machine scheduler script with the run and record the exact
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
| 801 complete SGS sets | about `0.39 TiB` |
| 41 full-resolution snapshots | about `0.16 TiB` |
| 11 restarts, scaled from the pilot | about `0.09 TiB` |
| Committed `16384^2` output plan | about `0.64 TiB` |
| Matching `12288^2` output plan | about `0.36 TiB` |

The committed inputs decouple full-resolution and SGS cadence: SGS products are
written 20 times per forcing turnover, primitive states once per turnover, and
restarts every four turnovers. Confirm filesystem quotas and measured I/O bandwidth
before retaining the full schedule.

## 3.5 Known issues

- No machine-specific scheduler script exists yet.
- No GPU or target-supercomputer validation has been run.
- No high-resolution forcing-render benchmark exists.
- No parallel-filesystem I/O benchmark exists.
- No spectral flux diagnostic exists in the current analysis workflow.
- No resolved $k^{-3}$ forward range has been demonstrated.
- The explicit-viscosity resolution estimate has not yet been checked against a
  measured enstrophy injection and dissipation budget.
- No initializer exists to promote a statistically steady lower-resolution state to
  a higher-resolution production mesh.

## 3.6 Continuation instructions

Before the first expensive submission:

1. Pull `Tiegan_SGS` and record the exact commit.
2. Build with the target machine's supported MPI and Kokkos configuration.
3. Run the focused CPU or GPU tests and the MPI regression.
4. Run a short `mach010_12288_k16_viscous.athinput` job and verify all output readers.
5. Create a machine-specific scheduler script without changing the committed physics.
6. Benchmark candidate MeshBlock sizes with the fixed 64-mode forcing.
7. Measure output and restart time with the intended production cadence.
8. Confirm that all requested coarsening factors divide the active MeshBlock sizes.
9. Confirm exact zero third momentum and third kinetic energy.
10. Confirm restart continuity before the long run.
11. Measure the enstrophy injection and viscous dissipation scales in the `12288^2` run.
12. Run the short `16384^2` convergence pilot and compare common resolved scales.
13. Freeze the input, scheduler script, Git commit, random seed, and output policy.
14. Start production only after memory, runtime, I/O, convergence, and stationarity gates pass.

No independent subagent review was run because delegation was not requested. The
report was directly audited against the repository state, completed pilot outputs,
focused test results, and the standard project-status report contract.
