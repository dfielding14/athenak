# Example: Localized And AMR Turbulence Driving

These short `turb` problem-generator calculations demonstrate the driver
controls and write the applied acceleration with `variable=turb_force`. They
are validation examples, not statistically converged turbulence experiments.

The default mode is `normalization=accel_rms`, with an explicit `accel_rms`
strength. The power-controlled examples select `normalization=edot`
explicitly. All cases use exact half kicks before and after the fluid step.

## Build And Run

From the repository root:

```bash
cmake -S . -B build-turb -DPROBLEM=turb -DCMAKE_BUILD_TYPE=Release
cmake --build build-turb -j

mkdir -p run/turb_edot_fixed run/turb_tiled_include run/turb_amr_exclude
build-turb/src/athena -d run/turb_edot_fixed \
  -i inputs/hydro/turb_driving/constant_edot_fixed_grid.athinput
build-turb/src/athena -d run/turb_tiled_include \
  -i inputs/hydro/turb_driving/tiled_localized_include.athinput
build-turb/src/athena -d run/turb_amr_exclude \
  -i inputs/hydro/turb_driving/amr_localized_exclude.athinput
```

The `-d` option keeps each case's `bin/` and history files in its run
directory.

## Cases

### Constant Energy Injection

`constant_edot_fixed_grid.athinput` uses a $32^3$ periodic uniform grid and

```text
normalization = edot
dedt          = 0.1
```

It supplies the homogeneous baseline for checking requested energy injection.
In this periodic box, with no other energy sources, the total energy should
increase by `dedt * volume * dt` each driven cycle, including cycle zero when
`tdriv_start=0`. The kick itself preserves thermal energy; later fluid
evolution can convert kinetic energy to heat.

### Tiled Included Region

`tiled_localized_include.athinput` uses

```text
normalization = accel_rms
accel_rms     = 0.3
tile_nx       = 2
tile_ny       = 2
localization  = include
sigma_x1      = 0.18
sigma_x2      = 0.18
sigma_x3      = 0.28
```

The same realization repeats over $2\times2\times1$ tiles and is
concentrated near the origin.

### Adaptive Excluded Region

`amr_localized_exclude.athinput` uses location-based adaptive refinement
around the origin, the same $2\times2\times1$ tile pattern, and
`localization=exclude`. Its refined central region is suppressed while
driving remains active outside. It selects `normalization=edot` with `dedt=0.1`.

## Projected Acceleration

The committed figure predates the full-timestep kick update and illustrates
the spatial controls; it is not validation of the current integration. It
was generated from cycle-12 `turb_force` dumps. Each panel
shows $L_3^{-1}\int|\boldsymbol{a}|\,dx_3$; the cyan contour in the AMR
panel marks the level-0/level-1 interface.

![Projections of turbulence-driving examples](../_static/turbulence_driving_projections.png)

Checks to make when rerunning the committed inputs are:

| Case | Check |
| --- | --- |
| Constant energy injection | Per-cycle total-energy increase equals `dedt * volume * dt` within numerical error. |
| Tiled included region | Volume-weighted RMS of `turb_force` equals `accel_rms=0.3`; acceleration is concentrated near the selected region. |
| AMR excluded region | Refinement does not add statistical weight to fine cells; acceleration is suppressed near the selected region. |

These energy and RMS checks concern the applied force. The raw modal OU
coefficients are subsequently localized, mean-subtracted, and normalized.
Exact injected work is separate from the temporal accuracy of the coupled
fluid evolution. The current driver uses symmetric splitting; the
single-kick measurements below document its first-order baseline.

## Symmetric-Kick Validation And Cost

The 2026-10-12 symmetric implementation passed 44 serial regressions and two
two-rank MPI regressions. The tests cover exact kick work, hydro/MHD temporal
convergence, event timing, strong initial forcing, STS viscosity, two-fluid
integration, AMR, and restart continuation. They also reject a long tail of
tiny steps near the final forcing event.

All three committed examples also ran for 12 cycles with the production
`turb` problem generator, increasing only history-output precision for the
check. Maximum per-cycle work discrepancies were $4.3\times10^{-10}$ for the
uniform power-controlled box and $2.9\times10^{-10}$ for the AMR case, against
an initial box energy of approximately 1000. The included-region RMS was
0.300000000002 when measured from its binary force output.

The fixed-grid, frozen-forcing experiment described below isolates temporal
error by comparing both methods with the same fine-step symmetric reference
at CFL 0.003125. At CFL 0.4, both methods take 83 steps to time 0.2:

| Quantity | Single-kick relative error | Symmetric relative error | Reduction |
| --- | ---: | ---: | ---: |
| Density perturbation | $9.306\times10^{-3}$ | $3.241\times10^{-5}$ | $287\times$ |
| Momentum | $1.0303\times10^{-2}$ | $3.202\times10^{-5}$ | $322\times$ |

These are errors relative to the small physical perturbation, not the mean
density, and exclude spatial discretization error. RK3 density/momentum
self-convergence orders change from approximately 1.00 to 2.00. Independent
RK2 and MHD cases also converge at approximately second order; the MHD tests
use a nonzero transverse magnetic field.

To reach comparable density and momentum errors in this smooth hydro test,
the single-kick method needed CFL 0.00123 and 26,894 steps, versus 83 symmetric
steps at CFL 0.4. This illustrates the potential gain when temporal error
controls the calculation; it is not a speedup prediction for turbulent statistics.

Runtime comparisons used commit `de43aaca9` as the single-kick baseline and
the symmetric turbulence-driver source with SHA-256 prefix `7e7ae5246fac`.
Both used AppleClang 21, Release builds, Kokkos Serial, an Apple M4 Max,
and the uniform ideal-MHD regression fixture: $32^3$ cells, $16^3$
MeshBlocks, RK3, PLM/HLLE, strength 0.1, and no outputs. Reported times are
median executable integration times from three alternating repetitions;
two-rank MPI runs used the same machine.

With one frozen realization (`dt_update=1000`), 100 cycles isolate the cost
per cycle:

| Backend / normalization | Single kick (s) | Symmetric (s) | Extra cost |
| --- | ---: | ---: | ---: |
| Serial / `accel_rms` | 2.711 | 3.240 | 19.5% |
| Serial / `edot` | 2.703 | 3.621 | 33.9% |
| MPI, 2 ranks / `accel_rms` | 1.415 | 1.675 | 18.4% |
| MPI, 2 ranks / `edot` | 1.412 | 1.830 | 29.6% |

With `dt_update=0.01` and `tcorr=0.25`, both methods were also run to the same
physical time 0.5. This comparison includes the additional steps needed to
land on OU events and respect the first kick's CFL margin:

| Backend / normalization | Single kick (s) | Symmetric (s) | Extra cost | Steps, single → symmetric |
| --- | ---: | ---: | ---: | ---: |
| Serial / `accel_rms` | 2.426 | 3.262 | 34.5% | 89 → 100 |
| Serial / `edot` | 2.827 | 4.963 | 75.5% | 104 → 139 |
| MPI, 2 ranks / `accel_rms` | 1.272 | 1.689 | 32.8% | 89 → 100 |
| MPI, 2 ranks / `edot` | 1.429 | 2.530 | 77.0% | 104 → 139 |

The extra work is a first-half kick, a boundary/primitive refresh, a second
normalization, and a global CFL reduction. Fixed-power forcing can require
another normalization after the timestep changes. No additional cell-sized
arrays are allocated. Balancing event remainders reduces avoidable timestep
recovery; the CFL calculation reuses the fluid's existing stability bound.

For RMS driving, the measured cost buys a clear second-order coupling
improvement while preserving exact kick work. Whether this changes a
particular turbulent statistic still requires timestep and sampling studies.
The OU experiment has an additional event-alignment benefit, so its error
reduction cannot be attributed solely to symmetric splitting. Continuous-OU
accuracy, GPU performance, and large-rank scaling were not measured.
Fixed-power normalization retains the singular-state qualification explained
in {doc}`../modules/turbulence_driving`.

## Single-Kick Baseline

The 2026-10-12 single-kick validation at `de43aaca9` passed 35 serial regressions and
two additional two-rank MPI regressions. All three production examples above
also ran successfully: both power-controlled cases matched their per-cycle
energy budgets, and the included-region case retained RMS acceleration 0.3.

Matched AppleClang/Serial timings used the uniform MHD regression fixture,
a $32^3$ grid with $16^3$ MeshBlocks, RK3, PLM/HLLE, 100 cycles, and no outputs.
Both modes used strength 0.1 and one frozen modal realization
(`dt_update=1000`). Three repetitions alternated old/new executable order.
The table reports median executable integration time against `b6035d254`:

| Normalization | Before (s) | Single kick (s) | Ratio |
| --- | ---: | ---: | ---: |
| `accel_rms` | 2.611135 | 2.531720 | 0.970 |
| `edot` | 2.625087 | 2.559827 | 0.975 |

A separate time-refinement experiment used a periodic one-dimensional unit
box with 128 cells, uniform ideal hydro ($\rho=p=1$, $\gamma=5/3$), RK3/PLM/HLLC,
and one compressive mode (`nlow=nhigh=1`, `spectrum=power_law`,
`min_kx=max_kx=1`, `max_ky=max_kz=0`, `sol_fraction=0`, `rseed=2718`). With
`accel_rms=0.01`, `dt_update=1`, `tcorr=1e6`, and final time 0.2, the modal field
is frozen. Halving CFL through 0.4, 0.2, 0.1, and 0.05 gave density/momentum
self-convergence orders 0.997--1.001 on the same spatial grid; CFL 0.003125
provided an additional reference. The applied acceleration varied by less
than $6\times10^{-12}$. This confirms the expected first-order deterministic
coupling and provides the baseline for the symmetric-splitting comparison.
It does not establish stochastic convergence or GPU performance.

## Earlier Modal-State Benchmark

Three Release runs of each input were compared against the initial feature
prototype at commit `a84c7ffd`; the table reports median executable CPU time.
These measurements predate the full-timestep kick update and have not been
repeated for the current integration.
These short total-runtime measurements include integration and output work,
not only driver kernels.

| Case | Prototype CPU time (s) | Modal-state CPU time (s) | Ratio |
| --- | ---: | ---: | ---: |
| Constant energy injection | `0.12035` | `0.11882` | `0.987` |
| Tiled included region | `0.22149` | `0.24983` | `1.128` |
| AMR excluded region | `0.27642` | `0.28859` | `1.044` |

The correctness change removes the prototype's two cell-sized OU scratch
arrays (`force_tmp1` and `force_tmp2`) and replaces them with two modal
innovation arrays. For these inputs, at double precision, the removed cell
storage is `2.93 MiB`, `5.86 MiB`, and `40.50 MiB`, respectively; the new
modal innovations use `1.31 KiB` in all three cases. The measured tiled-case
overhead is retained as evidence for any future kernel-level optimization
rather than hiding it behind unprofiled changes.

Regenerate the figure after running the examples:

```bash
python vis/python/plot_turbulence_driving_projection.py \
  --panel "Constant Edot=run/turb_edot_fixed/bin/turb_edot_fixed.force.00004.bin" \
  --panel "Tiled + Gaussian include=run/turb_tiled_include/bin/turb_tiled_include.force.00004.bin" \
  --panel "AMR + Gaussian exclude=run/turb_amr_exclude/bin/turb_amr_exclude.force.00004.bin" \
  --output docs/source/_static/turbulence_driving_projections.png
```

## Regression Coverage

The regression fixture is distinct from the `turb` examples. From the
repository root, build it and run the CPU tests from the binary directory:

```bash
cmake -S . -B build-turb-tests -DPROBLEM=tests/turb_forcing -DCMAKE_BUILD_TYPE=Release
cmake --build build-turb-tests -j
cd build-turb-tests/src
python3 -m pytest ../../tst/test_suite/turb/test_turb_driving_cpu.py
```

The turbulence regression tests exercise:

- conserved-state kick work and thermal-energy preservation with stale primitives;
- per-cycle fixed-power injection for RK1--RK4, including the first cycle;
- second-order frozen-forcing convergence for RK2/RK3 in hydro and MHD;
- forcing start/stop and OU-event alignment, with one modal realization per interval;
- a forcing-aware timestep for a strong first half kick;
- exact final-time arrival without a long tail of tiny event-clipped steps;
- per-cycle fixed-power injection with RKL2 viscosity sweeps in hydro and MHD;
- combined two-fluid work and momentum for ImEx2 and ImEx2+;
- default RMS acceleration and zero-target forcing;
- zero power with opposing flow and cancellation-resistant power normalization;
- fixed-RMS normalization with tiling and inclusion localization;
- AMR execution with exclusion localization;
- exact rendered-force and final-history continuation across fixed-grid and AMR restarts;
- MPI execution of the AMR/localized case and periodic MHD energy budget;
- rejection of removed parameter names and singular parabolic spectra.

The driver uses precision-correct MPI reductions. The current upstream source
tree does not build a complete single-precision executable with AppleClang
because unchanged coordinate and tabulated-EOS translation units contain
single-precision type/narrowing errors; those blockers are outside this
feature.

The equations, restart contract, and full canonical parameter reference are
in {doc}`../modules/turbulence_driving`.
