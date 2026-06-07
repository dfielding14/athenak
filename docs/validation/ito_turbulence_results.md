# Ito-2 Old-versus-Corrected Turbulence Results

## Bottom line

The covariance correction matters mathematically and changes individual tracer
realizations substantially. In this production-like turbulence test it does
**not** appreciably change the headline tracer-gas correlation or density-ratio
PDF. It produces a small, reproducible, scale-dependent change in the spectra.

The corrected implementation is also more expensive. On the eight-rank local
production run it took 1.94 times as long as the old diagonal implementation.

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

## Runtime

| Quantity | Old Ito-2 | Corrected Ito-2 | Ratio |
| --- | ---: | ---: | ---: |
| Wall time | 24.87 s | 48.33 s | 1.943 |
| Child user time | 192.48 s | 378.58 s | 1.967 |
| Child system time | 3.47 s | 3.78 s | 1.089 |

Peak RSS is not reported because `RUSAGE_CHILDREN` exposes a cumulative peak
that cannot be attributed safely to one lane of a sequential pair.

## Interpretation

For the reported turbulence observables, the covariance correction is not a
headline-changing effect at this resolution and particle count:

- correlation changes by $2.5\times10^{-4}$;
- $R^2_{1:1}$ changes by $1.1\times10^{-3}$;
- the density-ratio PDF width changes by about 2.2 percent;
- the PDF divergence is very small.

The spectral response is more visible. Relative to the old error, the
large-scale error improves by about 22 percent, the mid-scale error worsens by
about 8 percent, and the small-scale error worsens by about 2 percent. The
three pilot seeds show the same sign for these spectral changes.

The practical conclusion is:

1. fix the covariance because the old multidimensional process is
   mathematically wrong;
2. do not expect the main turbulence correlation/PDF conclusions to move much
   in this setup;
3. rerun any claim that depends on scale-by-scale tracer power; and
4. budget roughly a factor of two in local CPU runtime for the corrected
   implementation until the factorization and coefficient path are optimized.

## Artifacts

- `ito_turbulence_production_summary.json`
- `ito_turbulence_pilot_ensemble.json`
- `ito_turbulence_pilot_ensemble.csv`
- `ito_turbulence_column_density_comparison.png`
- `ito_turbulence_tracer_gas_ratio_pdf.png`
- `ito_turbulence_density_power_spectra.png`
