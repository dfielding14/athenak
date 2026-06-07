# Ito-2 Old-versus-Corrected Turbulence Validation

## Question

Does restoring the finite-step multidimensional covariance change the turbulence
results enough to matter?

The comparison is deliberately paired:

- identical gas initial conditions;
- identical particle positions and random seeds;
- identical resolution, solver, integrator, CFL number, outputs, and decomposition;
- legacy diagonal Ito-2 in one executable;
- full-covariance Ito-2 in the other.

Particles do not back-react on the gas. The old and corrected gas fields must therefore
agree before any tracer comparison is interpreted.

The current executable can perform the same comparison without maintaining two
binaries by overriding:

```text
particles/ito_covariance_model=published_diagonal
particles/ito_covariance_model=full_finite_step
```

Use this same-executable form for performance measurements. It isolates the
covariance model from unrelated implementation and compiler changes.

## Experiment tiers

| Input | Grid | Particles | End time | Purpose |
| --- | ---: | ---: | ---: | --- |
| `ito_tracers_turbulence_preflight.athinput` | $16^3$ | 8,192 | 8 cycles | CFL qualification and failure detection |
| `ito_tracers_turbulence_pilot.athinput` | $32^3$ | 131,072 | $0.05L/c_s$ | tractable local end-to-end pilot |
| `ito_tracers_turbulence_production.athinput` | $64^3$ | 4,194,304 | $0.2L/c_s$ | paper-scale comparison with 16 particles per cell |

All three start from the same deterministic, solenoidal Fourier velocity field with
$\sigma_{\rm rms}=15c_s$, then decay isothermally. The gas seed is `8675309`; the
particle transport and seeding seeds are `271828` and `161803`.

The production grid, particle count, Mach number, and end time follow the turbulence
test in Moseley, Teyssier, and Abel, arXiv:2604.23041v2, Section 3.3. This is not a
bit-for-bit reproduction of their RAMSES initial snapshot. AthenaK uses PLM plus HLLE
here because its HLLC solver is unavailable for an isothermal EOS. It is an
AthenaK-native, deterministic decaying-turbulence test designed to isolate the covariance
change.

## CFL qualification

AthenaK's RK2 driver has a formal single-operator stability limit of

$$
C_{\rm CFL}\leq 1.
$$
The categorical MC kernel requires the sum of all outward face probabilities
to remain at or below one. Write

$$
Q=D-mm^{\mathsf T},
\qquad
D_{ii}=h_i^2C_{+,i},
\qquad
m_i=h_iC_{-,i}.
$$
For positive diagonal $D$,

$$
Q\succeq 0
\quad\Longleftrightarrow\quad
m^{\mathsf T}D^{-1}m
=
\sum_i\frac{C_{-,i}^2}{C_{+,i}}
\leq 1.
$$
For one-way transport in each active direction,

$$
\frac{C_{-,i}^2}{C_{+,i}}=C_i,
$$
the covariance condition reduces to

$$
C_x+C_y+C_z\leq1.
$$
This motivates an initial three-dimensional guard

$$
C_{\rm CFL}\leq\frac{p_{\rm target}}{3}.
$$
The corrected code enforces this guard from the first step using
`particles/ito_probability_target`. With the default $p_{\rm target}=0.99$,
it is

$$
C_{\rm CFL}=0.99/3=0.33.
$$
This is not a universal stability proof. A divergent cell can have outward
flux through both faces of one coordinate, so the sum over faces need not be
bounded by three times the configured CFL number. The measured outgoing
probability tightens later steps predictively, and every current step still
fails closed if the realized finite-volume fluxes violate the categorical
condition.

The preflight sweep uses

```text
0.25, 0.30, 0.32, 0.329, 0.33, 0.331, 0.34
```

The eight-cycle preflight identifies whether the initial guard changes the gas
timestep. It must be followed by a full-duration qualification using the
production gas setup and a reduced tracer count. On the June 7, 2026 run,
`0.324` completed to $t=0.2$, while `0.3245` failed closed. The production
comparison therefore uses `0.324`, within 0.15 percent of the observed
stability edge. This value is empirical and flow-specific.

## Build the two executables

Use the documented baseline commit for the old implementation:

```bash
git worktree add --detach /tmp/athenak-ito-old e26bf0ec4
git -C /tmp/athenak-ito-old submodule update --init --recursive

cmake -S /tmp/athenak-ito-old -B /tmp/athenak-ito-old-build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/athenak-ito-old-build -j 8
```

Build the corrected working tree separately:

```bash
cmake -S . -B /tmp/athenak-ito-corrected-build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/athenak-ito-corrected-build -j 8
```

For production MPI runs, add `-DAthena_ENABLE_MPI=ON` to both configurations. Do not
compare binaries built with different precision, reconstruction options, or Kokkos
backends.

The runner records:

- executable SHA-256 and `athena -c` configuration;
- input SHA-256;
- exact command and parameter overrides;
- repository SHA, dirty state, and a binary Git patch;
- hostname and selected environment variables;
- wall, user, and system time;
- exit status and complete AthenaK log.

## Run the CFL preflight

```bash
python3 scripts/run_ito_turbulence_validation.py preflight \
  --old-exe /tmp/athenak-ito-old-build/src/athena \
  --corrected-exe /tmp/athenak-ito-corrected-build/src/athena \
  --input inputs/particles/ito_tracers_turbulence_preflight.athinput \
  --output-root /tmp/ito-turbulence-preflight
```

`preflight_summary.json` reports every lane and selects the highest passing CFL at or
below $0.33$. A lane passes only when:

1. both executables finish without fatal text or non-finite output;
2. gas and tracer mesh outputs exist at the same time;
3. gas density is finite and positive;
4. particle density is finite and non-negative; and
5. old and corrected gas fields agree to the requested tolerance.

Use the short preflight value for the pilot. Before production, repeat the
production gas evolution with fewer tracers through the full end time and
bracket the failure threshold. Do not call `0.33` a mathematical maximum.

## Run the local pilot

```bash
python3 scripts/run_ito_turbulence_validation.py run-pair \
  --old-exe /tmp/athenak-ito-old-build/src/athena \
  --corrected-exe /tmp/athenak-ito-corrected-build/src/athena \
  --input inputs/particles/ito_tracers_turbulence_pilot.athinput \
  --output-root /tmp/ito-turbulence-pilot \
  --cfl 0.33
```

The pilot is successful when the full runner and analysis pipeline completes. It is not
large enough for a publication-level statement.

## Run the production comparison

Example with eight MPI ranks:

```bash
python3 scripts/run_ito_turbulence_validation.py run-pair \
  --old-exe /tmp/athenak-ito-old-mpi-build/src/athena \
  --corrected-exe /tmp/athenak-ito-corrected-mpi-build/src/athena \
  --input inputs/particles/ito_tracers_turbulence_production.athinput \
  --output-root runs/ito-turbulence-production-seed-161803 \
  --launcher "mpirun -np 8" \
  --cfl 0.324
```

Repeat at least three paired particle seeds while leaving the gas seed unchanged:

```bash
python3 scripts/run_ito_turbulence_validation.py run-pair \
  --old-exe /tmp/athenak-ito-old-mpi-build/src/athena \
  --corrected-exe /tmp/athenak-ito-corrected-mpi-build/src/athena \
  --input inputs/particles/ito_tracers_turbulence_production.athinput \
  --output-root runs/ito-turbulence-production-seed-271829 \
  --launcher "mpirun -np 8" \
  --cfl 0.324 \
  --override particles/random_seed=314159 \
  --override tracer_seed1/seed=271829
```

Use the same two overrides for both lanes of each pair. Compare paired metric
differences across seeds; do not estimate uncertainty by treating grid cells as
independent samples.

For a same-executable mode comparison, run the production input twice with the
two covariance overrides above, then pass the published run as `--old-run` and
the full run as `--corrected-run` to
`scripts/analyze_ito_turbulence_validation.py`.

## Diagnostics

The final paper-facing tracer field is deposited from particle positions using periodic
cell-centered CIC. The intermediate native `prtcl_d` output uses NGP and is retained for
health checks and time evolution.

The analysis writes:

| Output | Contents |
| --- | --- |
| `summary.json` | snapshot metadata, all scalar metrics, corrected-minus-old changes |
| `metrics.csv` | long-form scalar metric table |
| `ratio_pdf.csv` | PDF of $\log_{10}(\rho_t/\rho_g)$ |
| `power_spectra.csv` | gas and tracer density spectra |
| `analysis_arrays.npz` | normalized 3D fields, columns, PDFs, and spectra |
| `column_density_comparison.png` | gas, old, corrected, ratios, and difference |
| `tracer_gas_ratio_pdf.png` | old and corrected ratio PDFs |
| `density_power_spectra.png` | spectra and relative spectral error |

The scalar metrics include:

$$
r(\rho_g,\rho_t),\quad
R^2_{\rm fit},\quad
R^2_{1:1},
$$
$$
\mu\!\left[\log_{10}\frac{\rho_t}{\rho_g}\right],\quad
\sigma\!\left[\log_{10}\frac{\rho_t}{\rho_g}\right],
$$
column-density correlation and residual norms, and mean absolute tracer/gas spectral
error in large-, intermediate-, and small-scale bands.

## Interpretation

The experiment is invalid if the gas fields differ appreciably. With the same
decomposition and unchanged gas code, the target is

$$
\frac{\|\rho_{g,\rm corrected}-\rho_{g,\rm old}\|_2}
{\|\rho_{g,\rm old}\|_2}\leq10^{-12}.
$$
For tracer metrics, report the signed paired change:

$$
\Delta M=M_{\rm corrected}-M_{\rm old}.
$$
The correction is scientifically material when the paired change is reproducible across
particle seeds and is comparable to, or larger than, the differences among tracer
methods reported in the paper. A visible single-seed plot is not enough.

At minimum, report:

- $\Delta R^2_{1:1}$;
- change in the mean and width of the log density-ratio PDF;
- change in column-density residual norms;
- change in spectral error by scale band;
- Jensen-Shannon divergence between old and corrected ratio PDFs;
- runtime and memory change; and
- any CFL-dependent failure or covariance-factorization rejection.

## Script checks

```bash
python3 -m py_compile \
  scripts/analyze_ito_turbulence_validation.py \
  scripts/run_ito_turbulence_validation.py \
  scripts/test_ito_turbulence_validation.py

python3 scripts/test_ito_turbulence_validation.py -v
```

The unit tests cover CIC mass conservation, exact identical-field metrics, Fourier-mode
recovery, PDF divergence, particle VTK parsing, and command construction.
