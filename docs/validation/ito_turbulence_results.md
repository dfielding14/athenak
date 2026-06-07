# Ito-2 Published-Diagonal versus Full-Covariance Turbulence Results

## Bottom line

The full finite-step model matters mathematically and changes individual tracer
realizations substantially. In this production-like turbulence test it does
**not** appreciably change the headline tracer-gas correlation or density-ratio
PDF. It produces a small, reproducible, scale-dependent change in the spectra.

The full implementation is also more expensive. The historical two-executable
comparison measured 1.94x. After unrolling the hot-path covariance
factorization, a repeated same-executable comparison measured 1.47x. The published
diagonal model is therefore the default; full finite-step covariance is an
explicit runtime choice.

## Run identity

| Item | Value |
| --- | --- |
| Date | June 7, 2026 |
| Old source | `e26bf0ec4724a9204bb74e8efcf4d0b8f2dfcd16` |
| Old MPI executable SHA-256 | `f953fc80fc1c70275dd634aeb2e672f1c4bd881caef272426fd5472763cf3416` |
| Corrected MPI executable SHA-256 | `7554d0a6b4f36a4724122186c97766a9d92204c02da14ed037ae3d9e7c13b2f5` |
| Configuration | double precision, MPI, Kokkos Serial |
| MPI ranks | 8 |
| Production input | `inputs/particles/ito_tracers_turbulence_production.athinput` |
| Resolution | $64^3$ |
| Particle count | 4,194,304, or 16 per cell |
| Gas seed | `8675309` |
| Production particle seeds | transport `271828`, seeding `161803` |
| Integrator / solver | RK2, PLM, HLLE, isothermal Hydro |
| End time | $0.2L/c_s$ |
| Selected CFL | `0.324` |
| Runtime parameter | `particles/ito_covariance_model` |

This is an AthenaK-native deterministic decaying-turbulence experiment. It
matches the paper's production resolution, particle count, Mach number, and
end time, but it is not a bit-for-bit replay of the RAMSES snapshot in
Moseley, Teyssier, and Abel.

## CFL qualification

The short eight-cycle sweep passed with identical gas fields through `0.33`.
At `0.331` the corrected first-step guard changed the gas timestep, so the
paired gas fields diverged as expected.

That short sweep was not sufficient. On the full $t=0.2$ production gas
evolution:

| CFL | Outcome |
| ---: | --- |
| 0.320 | completed |
| 0.322 | completed |
| 0.323 | completed |
| 0.324 | completed |
| 0.3245 | failed closed |
| 0.325 | failed closed |
| 0.330 | failed closed |

The production value `0.324` is within 0.15 percent of the observed failure
edge. This is an empirical, flow-specific stability result. The
`ito_probability_target / active_dimensions` guard is not a general proof
because a divergent cell can have outward mass flux through both faces of a
coordinate.

## Three-seed pilot

Three independent $32^3$, 131,072-particle pairs were run to
$0.05L/c_s$ at CFL `0.33`. Values are corrected minus old; uncertainties are
paired standard errors across the three particle seeds.

| Metric | Mean change | Standard error |
| --- | ---: | ---: |
| Density Pearson $r$ | $+1.636\times10^{-3}$ | $8.29\times10^{-5}$ |
| Density $R^2_{1:1}$ | $-1.882\times10^{-3}$ | $1.25\times10^{-3}$ |
| Std. $\log_{10}(\rho_t/\rho_g)$ | $+1.508\times10^{-3}$ | $5.78\times10^{-3}$ |
| Column Pearson $r$ | $-6.753\times10^{-4}$ | $5.63\times10^{-4}$ |
| Large-scale spectral error | $-6.249\times10^{-3}$ | $1.66\times10^{-3}$ |
| Mid-scale spectral error | $+1.861\times10^{-2}$ | $5.62\times10^{-3}$ |
| Small-scale spectral error | $+1.209\times10^{-2}$ | $3.17\times10^{-3}$ |
| PDF Jensen-Shannon divergence | $1.015\times10^{-3}$ bits | $5.74\times10^{-5}$ |

The pilot therefore predicts small changes in integrated statistics and a
more consistent change in scale-dependent power.

## Production result

The old and corrected gas densities are bit-identical:

| Gas comparison | Value |
| --- | ---: |
| Relative $L_2$ difference | 0 |
| Maximum absolute difference | 0 |

Tracer-gas statistics:

| Metric | Old Ito-2 | Corrected Ito-2 | Corrected - old |
| --- | ---: | ---: | ---: |
| Density Pearson $r$ | 0.976168 | 0.975919 | -0.000249 |
| Density $R^2_{1:1}$ | 0.949262 | 0.948116 | -0.001146 |
| Density normalized $L_1$ | 0.149927 | 0.148817 | -0.001110 |
| Density normalized $L_2$ | 0.273741 | 0.275489 | +0.001748 |
| Mean $\log_{10}(\rho_t/\rho_g)$ | -0.002532 | +0.000038 | +0.002570 |
| Std. $\log_{10}(\rho_t/\rho_g)$ | 0.122601 | 0.119919 | -0.002683 |
| Column Pearson $r$ | 0.994382 | 0.994533 | +0.000151 |
| Column normalized $L_1$ | 0.035382 | 0.034907 | -0.000475 |
| Large-scale spectral error | 0.036589 | 0.028718 | -0.007872 |
| Mid-scale spectral error | 0.109345 | 0.118127 | +0.008782 |
| Small-scale spectral error | 0.464526 | 0.473717 | +0.009191 |

The old-versus-corrected density-ratio PDF Jensen-Shannon divergence is
`1.699e-4` bits.

The normalized old-versus-corrected tracer-field difference is much larger
than the summary-statistic changes:

| Field difference | Value |
| --- | ---: |
| Normalized $L_1$ | 0.138946 |
| Normalized $L_2$ | 0.198309 |

That is expected. Correlating the coordinate kicks changes each deterministic
tag trajectory, even when ensemble diagnostics remain close.

## Historical runtime

| Quantity | Old Ito-2 | Corrected Ito-2 | Ratio |
| --- | ---: | ---: | ---: |
| Wall time | 24.87 s | 48.33 s | 1.943 |
| Child user time | 192.48 s | 378.58 s | 1.967 |
| Child system time | 3.47 s | 3.78 s | 1.089 |

Peak RSS is not reported because `RUSAGE_CHILDREN` exposes a cumulative peak
that cannot be attributed safely to one lane of a sequential pair.

## Final same-executable runtime

Both modes were rerun twice from the same optimized executable
(`SHA-256 74289f3b4d7bc77617fab7573d9cc2e5238e792a766ace90a28e04a9c817bcb6`)
with the same production input, eight MPI ranks, and CFL `0.324`. The sequence
was published, full, full, published.

The final verification executable has SHA-256
`0f3a0cc342d634597aed5e29272ed2c15aa62cfd9ef758135dea823c0fb051c4`.
The changes between those binaries are confined to legacy-restart inference,
reserved-value validation, and tests; the coefficient construction,
factorization, and particle-push hot paths are identical.

| Quantity | `published_diagonal` | `full_finite_step` | Ratio |
| --- | ---: | ---: | ---: |
| Wall-time samples | 25.18, 25.30 s | 37.35, 37.07 s | 1.465-1.483 |
| Mean wall time | 25.24 s | 37.21 s | 1.474 |

The optimized factorization reduced the full-mode cost substantially from the
pre-optimization same-executable run. Absolute times varied between sessions,
but the final repeated mode ratio remained close to the earlier 1.48x result.
The science metrics from the final mode comparison agree with the historical
old-versus-corrected result to roundoff.

A separate 65,536-particle one-step check at CFL `0.495` compared
`published_diagonal` against the archived published executable. The timesteps
were identical; particle coordinates were bit-identical in `x1` and `x3`, with
a maximum `x2` difference of `1.11e-16`.

## Interpretation

This is not a pure covariance ablation. The kick vector is built from
independent bounded-uniform variables, not a rotationally invariant Gaussian.
Replacing the diagonal factor with a full factor therefore changes fourth and
higher moments of the continuous kick law along with its cross covariance. The
numbers below measure the complete `published_diagonal` versus
`full_finite_step` sampler change.

For the reported turbulence observables, that change is not headline-changing
at this resolution and particle count:

- correlation changes by $2.5\times10^{-4}$;
- $R^2_{1:1}$ changes by $1.1\times10^{-3}$;
- the density-ratio PDF width changes by about 2.2 percent;
- the PDF divergence is very small.

The spectral response is more visible. Relative to the old error, the
large-scale error improves by about 22 percent, the mid-scale error worsens by
about 8 percent, and the small-scale error worsens by about 2 percent. The
three pilot seeds show the same sign for these spectral changes.

The practical conclusion is:

1. use `published_diagonal` for compatibility with the paper and for production
   runs where the measured observable is insensitive to the joint kick tensor;
2. use `full_finite_step` when claiming finite-step vector covariance fidelity
   or interpreting scale-by-scale tracer power;
3. do not expect the main turbulence correlation/PDF conclusions to move much
   in this setup; and
4. budget about 47 percent more local wall time for the full mode on this
   benchmark. GPU performance remains unmeasured.

## Artifacts

- `ito_turbulence_production_summary.json`
- `ito_turbulence_covariance_modes.json`
- `ito_turbulence_pilot_ensemble.json`
- `ito_turbulence_pilot_ensemble.csv`
- `ito_turbulence_column_density_comparison.png`
- `ito_turbulence_tracer_gas_ratio_pdf.png`
- `ito_turbulence_density_power_spectra.png`
