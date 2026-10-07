# Reviewed reference: beta ten, seed 271828

Figures have been redrawn with `dbfplot`: compact panels, Dark2 colors,
inward ticks, and PDF/240-dpi PNG exports. Logarithmic PDFs, mathematical
labels and illustrative spectral slopes are retained. The pressure panels
identify fluctuation spectra, thermal–magnetic correlation and normalized
residual variance. [Strict figure audits](figure-audit.json) and
[rendering provenance](figure-rendering.json) are separate from the original
numerical analysis. Measurements and scientific findings below are unchanged.

“Range of 2-unit time averages” means the pointwise minimum and maximum of
six averaged curves over [6,8], [8,10], [10,12], [12,14], [14,16] and [16,18].
It measures temporal variability, not a confidence interval; the blocks can
remain correlated. The solid PDF curve averages the full [6,18] interval.

**Numerical integrity is consistent. Physical evidence is useful but remains
inconclusive for quantitative paper agreement or a causal magneto-immutability
claim.** The corrected box completed t=18, including two actual restarts. The
fixed [6,18] interval provides six duration-two blocks; their variability is
descriptive, with correlation and heating limitations retained. The
[scientific review](scientific-review.md) assesses all four figure groups and
also compares three duration-four blocks. No further extension is warranted
merely to improve agreement.

The achieved Mach proxy is **0.21366 ± 0.01626**, vector delta-B RMS/B0
**0.80664 ± 0.07136**, and volume-mean beta **9.2424 ± 1.4971** (block SD).
Resolved parallel-strain/transverse-gradient power is **0.004826**, while
parallel fluctuations supply **44.35%** of the integrated, mean-subtracted
local velocity variance. This is a persistent directional hierarchy in one
active realization; it does not establish an active/passive causal contrast.

Global perpendicular pressure compensation is modest and variable:
**C=−0.15199 ± 0.13213**, **R=0.87297 ± 0.10487**. The
[three-snapshot signed-band audit](pressure-episode/README.md) finds that the
t=14 global positive correlation is dominated by forcing scales, while the
resolved smaller-scale band remains anticorrelated at the inspected times.
These selected snapshots are not a second time-averaged ensemble.

The [mirror-tail audit](mirror-tail/README.md) retains **three cells** beyond
the float32 rounding envelope in one snapshot, with maximum X=1.00022452.
These small tails are compatible with finite soft-mirror relaxation; the
saved state cannot determine the exact pre-relaxation excess. No configured
mirror hard wall was violated, and no firehose crossing beyond the envelope
is observed. Threshold occupancy, strict exceedance and near-band residence
remain distinct.

The original restart exposed a three-cell collision-wall encoding defect.
The [numerical correction and checks](numerical-checks.md) document its
surgical fix, unchanged conservative-A formulation/strict thresholds, passing
regressions and preserved B4 sharp-contact limitations. For the corrected
[6,18] interval, actual applied work is **7.68** and energy change minus work
is **−5.72e−12**; floor/nonfinite/nonpositive counters stay zero.

Use the [reproduction guide](../../../cgl_lf_physics_benchmark.md) and
[reference README](README.md) for exact commands, costs and provenance.
The primary comparisons are [Squire 2023](https://arxiv.org/html/2303.00468v2)
and [MKS24](https://arxiv.org/html/2405.02418v2), with their forcing, threshold,
normalization and resolution differences stated in the guide. One box does
not establish spatial convergence, a precise cascade slope, or LF coefficient
accuracy; retain independent physical damping checks.

The following automatically generated report is preserved with trailing
whitespace removed. Its directional labels describe the measured diagnostics;
the joint scientific assessment above and linked review qualify their scope.

---

# CGL-LF single-run physics benchmark

**Inconclusive** for a finite-window descriptive comparison.

Requested interval [6.0, 18.0]; retained interval [6.0, 18.0]; 50 snapshots and 6 complete blocks of duration 2.0.

Block bands show temporal variability, not independent-snapshot confidence intervals. No cooling is present: evolving thermal energy/beta are reported rather than assumed stationary.

Bracketing snapshot times span [5.75034954004, 18.0]; maximum gap 0.2503838736000006. Diagnostic values (not cell fields) are linearly interpolated to covered requested endpoints. Snapshot blocks are anchored at t=6; this equals the requested start only when its coverage is available. No extrapolation is performed.


| Figure group | PDF |
| --- | --- |
| [marginality](marginality.png) | [vector PDF](marginality.pdf) |
| [pressure_balance](pressure_balance.png) | [vector PDF](pressure_balance.pdf) |
| [gradients](gradients.png) | [vector PDF](gradients.pdf) |
| [spectra_energy](spectra_energy.png) | [vector PDF](spectra_energy.pdf) |

| Directional finding | Classification | Reason |
| --- | --- | --- |
| marginality | inconclusive | report near-band residence and strict/rounding-robust crossings separately; snapshot precision and projection cadence do not establish full-time admissibility |
| pressure_balance | consistent | negative covariance and R<1 indicate compensation relative to uncorrelated fields; extent and variability require review, with no paper amplitude tolerance |
| gradients | inconclusive | resolved/cutoff tensor ratios and parallel velocity are measured; one active box cannot establish suppression caused by anisotropy feedback |
| spectra_energy | inconclusive | Parseval and measured forcing ledger available as stated; thermal drift and cascade extent require finite-window review, without a target slope |

Analysis normalization checks: `{"classification": "consistent", "max_parseval_relative_error": 4.188727624076505e-15, "max_PDF_integral_error": 1.2245759961615477e-12, "scope": "analysis normalization/finite-field checks, not simulation physics acceptance"}`. These do not classify simulation health or replace scientific review.

**Simulation integrity: consistent.**

Numerical floor/nonfinite/nonpositive counters: `{"lf_dfloor": {"available": true, "full_retained_min": 0.0, "full_retained_max": 0.0, "last": 0.0, "window_increment": 0.0}, "lf_pfloor": {"available": true, "full_retained_min": 0.0, "full_retained_max": 0.0, "last": 0.0, "window_increment": 0.0}, "lf_nonfin": {"available": true, "full_retained_min": 0.0, "full_retained_max": 0.0, "last": 0.0, "window_increment": 0.0}, "lf_nonpos": {"available": true, "full_retained_min": 0.0, "full_retained_max": 0.0, "last": 0.0, "window_increment": 0.0}}`.

Unprojected LF stage crossings and reserved restart field (not automatically failures): `{"lf_hardbd": {"available": true, "full_retained_min": 0.0, "full_retained_max": 3550875530.0, "last": 3550875530.0, "window_increment": 2922212058.8112035, "meaning": "cumulative active-cell hard-bound crossings at unprojected LF stages; repeated cell-stage events, not unique cells or a time-integrated volume fraction"}, "lf_hwproj": {"available": true, "full_retained_min": 0.0, "full_retained_max": 0.0, "last": 0.0, "window_increment": 0.0, "instrumentation": "reserved_uninstrumented", "audited_simulation_revision": "7a37710f6c224e24e7c7f364e7e0b812b3a9494c", "audited_simulation_revisions": ["7a37710f6c224e24e7c7f364e7e0b812b3a9494c", "71ad25ebce73d33db048defd8f585a7dce0528c2"], "meaning": "post-WO2 audited source initializes and serializes this field but never increments it; retained values are preserved for provenance and do not establish wall-projection activity or its absence"}}`.

**Projection-count limitation:** in the audited post-WO2 simulation revision, `lf_hwproj` is reserved/uninstrumented. Its raw value and increment do not measure wall-projection activity, and zero does not demonstrate absence of projections. `lf_hardbd` and nonfinite/nonpositive counters count active-cell stage checks; density/pressure-floor counters count EOS refresh events including refreshed halo cells. All are instrumented cumulative counts, preserved across restart and history output.

Post-operator hard-bound history: `{"available": true, "min_fraction": 0.0, "max_fraction": 0.0, "meaning": "strict hard-bound predicate at retained history times; distinct from unprojected LF stage crossings"}`.


| Measurement | Time mean | Block SD |
| --- | ---: | ---: |
| Mach_isotropic_proxy | 0.213656 | 0.0162595 |
| deltaB_rms_over_B0 | 0.806641 | 0.0713566 |
| beta_volume_mean | 9.2424 | 1.49711 |
| u_parallel_fraction | 0.449312 | 0.0268027 |
| mirror_strict | 0.00792503 | 0.00469342 |
| firehose_strict | 0.00220907 | 0.0018063 |
| mirror_near | 0.032882 | 0.0194274 |
| firehose_near | 0.0113955 | 0.00955712 |
| pressure_correlation | -0.151989 | 0.132133 |
| pressure_normalized_residual_variance | 0.872968 | 0.10487 |
| S_parallel_rms | 0.547608 | 0.0829221 |
| induction_rms | 0.4911 | 0.0817845 |

Pressure spectra alone cannot establish compensation; use the signed correlation and normalized residual together. Local S_parallel is not b·grad(u_parallel). The latter includes field-direction curvature; S_parallel-div(u) is the ideal compressible induction proxy. Parallel velocity is retained explicitly. Compare the resolved and cutoff bands separately. The conservative resolved guide is Nyquist/4 (eight cells per wavelength); scalar resolved ratios also restrict full |k|. Plotted k_perp shells sum all k_parallel, so small k_perp alone does not imply fully resolved gradients. Green shading marks the projection of the physical forcing shell onto k_perp, including zero; no universal spectral slope is prescribed.

Gradient band measurements: `{"resolved_strain_to_perpendicular_gradient_power": 0.00482573674039363, "cutoff_strain_to_perpendicular_gradient_power": 0.0825060677175919}`.

Actual applied forcing energy budget: `{"available": true, "interval": [6.0, 18.0], "requested_interval": [6.0, 18.0], "requested_window_covered": true, "actual_applied_work": 7.6800000000001205, "actual_mean_total_power": 0.64000000000001, "conserved_total_energy_change": 7.679999999994401, "residual_E_minus_work": -5.7198690228688065e-12, "relative_residual": 7.447746123526975e-13, "sampling": "linear interpolation only if explicit endpoints are between retained history rows"}`.

Retained dedt=0.32 is nominal power per domain volume; nominal total power=0.64. The measured accumulated work remains authoritative.
Realized Helmholtz solenoidal acceleration-power fraction: `{"mean": 0.5620842180121426, "block_starts": [6.0, 8.0, 10.0, 12.0, 14.0, 16.0], "block_duration": 2.0, "block_means": [0.6677452686912362, 0.6832273964609474, 0.7625585867826209, 0.421251661688755, 0.1964747301041403, 0.6412476643451565], "block_count": 6, "block_sd": 0.21250922693072058, "block_min": 0.1964747301041403, "block_max": 0.7625585867826209, "effective_window": [6.0, 18.0], "linear_slope": -0.025086573478371458, "fitted_change_over_window": -0.3010388817404575}`.
Nominal ratio of expected isotropic innovation powers: 0.5000000000000001. This is not the expectation of an instantaneous fraction or a target for this finite OU realization. Ratio of time-mean solenoidal to total acceleration powers (a different statistic): 0.5022943689732536. Neither statistic is an energy-injection partition.

| Temporal diagnostic | Integral correlation time | Caution |
| --- | ---: | --- |
| deltaB_rms_over_B0 | 1.1258318122988702 | Finite-window estimate; inspect drift |
| beta_volume_mean | 1.6539664330698884 | Finite-window estimate; inspect drift |
| S_parallel_rms | 0.7068672032012219 | Finite-window estimate; inspect drift |
| u_parallel_rms | 0.638406326184646 | Finite-window estimate; inspect drift |

- deltaB_rms_over_B0: block duration is less than twice the empirical positive-sequence integral correlation time.

- beta_volume_mean: block duration is less than twice the empirical positive-sequence integral correlation time.

| History drift | Fitted change in window | Block means |
| --- | ---: | --- |
| kinetic | -0.028779558688034217 | [0.5573199082871366, 0.5811362925139936, 0.48198229743123405, 0.5693089003835458, 0.5115918258136806, 0.5397739033557087] |
| magnetic | -0.07809632149006067 | [1.6436198846949839, 1.5858249590222486, 1.7775073597608282, 1.8211415321629998, 1.5322826039353925, 1.5770787187162112] |
| therm_cgl | 7.786875880168389 | [18.279060207012417, 19.593038748450795, 20.78051034278854, 21.929549567432062, 23.556125570231433, 24.76314737791554] |
| beta_volume_mean | 4.133604803696176 | [7.860418866545128, 8.570160721456665, 8.167555714525049, 8.570684373141622, 11.24582723252871, 11.042836931554705] |

Empirical autocorrelation and linear trends are retained in metrics.json. Positive-sequence correlation times are finite-window estimates, not evidence that OU-sized blocks are independent. Fewer than four complete blocks warrants extending this same realization before quoting variability; a longer heated run is still a finite-time ensemble.

Thresholds: mirror X=1.0, firehose X=-2.0; symmetric near halfwidth 0.05. Strict exceedance, inclusive solver history switches, and near-threshold bands are separate. Snapshot payload precision and post-operator/source/projection cadence limit admissibility claims; the dynamic half-ULP float32 sensitivity envelope is a representation bound, not altered physical thresholds or a bound on dynamical error.

One box does not establish active/passive causality, spatial convergence, a universal spectral slope, or LF coefficient correctness. The directional classification above uses explicitly stated coverage prerequisites and signed compensation evidence, not a fitted paper-image tolerance.

## Definitions

- **X**: 2*(p_perp-p_parallel)/B^2; AthenaK magnetic pressure is B^2/2.
- **pressure_input**: primitive binary 'eint' is p_parallel, not thermal energy.
- **PDF**: sum(cell volume in bin)/(domain volume * bin width); shared edges include every retained cell.
- **time_average**: trapezoidal integration of snapshot diagnostics over actual retained times, divided by covered duration.
- **blocks**: contiguous physical-time blocks; linearly interpolated diagnostic endpoints; sample SD/range of block means, not an iid confidence interval.
- **derivatives**: second-order centered periodic differences on the uniform full 3D mesh.
- **S_parallel**: b_i b_j partial_j u_i, b=B/|B| (local instantaneous direction).
- **projected_gradients**: G_ij=partial_j u_i, P_ij=delta_ij-b_i*b_j; parallel/perpendicular components are b.G.b, P.G.b, b.G.P, P.G.P, projected after differentiating u; no derivatives of b enter these four components.
- **parallel_velocity**: u_parallel=u dot b; this includes the retained bulk velocity.
- **parallel_velocity_derivative**: b dot grad(u dot b) = S_parallel + u dot [(b dot grad)b] in the continuum; discretization need not obey the product rule exactly.
- **induction**: S_parallel-div(u) reconstructs ideal material D ln|B|/Dt; not a measured time derivative or resistive/numerical induction budget.
- **pressure_balance**: corr(delta p_perp, delta(B^2/2)); R=mean[(delta p_perp+delta(B^2/2))^2]/(var(p_perp)+var(B^2/2)); exact compensation gives corr=-1,R=0.
- **FFT**: F=fftn(field)/N; shell power=sum_shell |F|^2 / dk; integral over all shells (including k_perp=0) equals the stated real-space mean square.
- **shells**: k_perp=sqrt(kx^2+ky^2) relative to the initial z guide field; physical radians/length; bins [n*dk,(n+1)*dk), dk=min(2*pi/Lx,2*pi/Ly).
- **kinetic_spectrum**: FFT of sqrt(rho/2)*(u-u_bulk), u_bulk=<rho*u>/<rho>; no further mean removal; integral=mean[rho*|u-u_bulk|^2/2].
- **magnetic_spectrum**: FFT of (B-<B>)/sqrt(2); integral=mean[|B-<B>|^2/2].
- **pressure_spectra**: FFT of each scalar minus its volume mean; unnormalized physical pressure units; integral=variance.
- **gradient_spectra**: FFT of each reconstructed scalar/vector minus its component means; signed fields are squared only by the power spectrum; spectra do not determine signs.
- **Mach**: u_rms about volume-mean velocity / sqrt(gamma*<p_iso>/<rho>), p_iso=(p_parallel+2*p_perp)/3; isotropic-pressure proxy, not a CGL characteristic-wave Mach number.
- **beta**: volume mean of 2*p_iso/B^2; also report 2*<p_iso>/<B^2>, which differs in general.
- **deltaB**: sqrt(<|B-<B>|^2>)/B0; mean field removed separately at each snapshot.
- **forcing**: actual accumulated applied forcing work is user-history force_work; force_pwr is instantaneous rho*u dot f and is not its exact quadrature.
- **force_decomposition**: nonzero Fourier acceleration modes: P_compressive=sum |k dot fhat|^2/k^2; P_solenoidal=P_total-P_compressive; this is acceleration power, not injected-energy partition.
- **nominal_force_mixture**: expected_solenoidal_power_fraction and expected_solenoidal_fraction retain the nominal ratio of expected isotropic innovation powers, 2*s^2/[2*s^2+(1-s)^2]; this is not the expectation of an instantaneous fraction or a target for one finite OU realization.
- **precision**: primitive/force snapshots may be float32; strict threshold crossings are descriptive and not tight full-precision admissibility tests.
- **X_rounding_envelope**: half the larger adjacent float32 gap for each stored pressure/B component; bound |delta X| <= [2*(delta p_perp+delta p_parallel)+|X|*delta(B^2)]/[B^2-delta(B^2)], delta(B^2)=sum(2*|Bi|*delta Bi+delta Bi^2); no dynamical-error claim.
- **limiter_history**: mirror_vol/fire_vol are inclusive threshold predicates; with nu_coll=0, both soft limiters enabled and backups disabled, nu_eff/(limiter_nu_coll*volume) measures the strict soft-rate fraction at history sampling; hard_vol counts strict physical firehose violation even with backups off.
- **LF_counters**: lf_nstage counts cumulative active-cell LF-stage checks; lf_hardbd counts cumulative hard-bound crossings in those unprojected checks, with repeated cells counted repeatedly. nonfin/nonpos also inspect active cells; dfloor/pfloor accumulate EOS refresh events including refreshed halo cells. All persist restart and do not reset at history output. In audited post-WO2 simulation revisions 7a37710f6 and 71ad25ebc, lf_hwproj is a reserved, uninstrumented field: its retained value does not measure wall projections or their absence. Physical stage crossings are not automatically numerical failure.

## Provenance

Simulation revision from retained launch metadata: `71ad25ebce73d33db048defd8f585a7dce0528c2`. Analysis checkout revision: `64d652317bec173d1d33e23906d0e93f21463409` (separate from simulation provenance).

[metrics.json](metrics.json) includes embedded effective input, file hashes, analysis/reader hashes, runtime metadata, sampling, block means, Parseval errors, and restart-time deduplication audits.
