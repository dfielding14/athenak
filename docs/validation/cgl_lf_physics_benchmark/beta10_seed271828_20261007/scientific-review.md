# Final independent scientific review: corrected canonical [6,18]

**Recommendation: retain this as the completed finite-window reference and do not extend further merely to improve agreement.** The extra interval captures the strongly compressive forcing episode and its subsequent recovery. Six blocks support a descriptive record of this realization; they do not support stationary thermal statistics, independent-block confidence intervals or quantitative agreement with paper figures.

**Overall scientific classification: inconclusive for a quantitative or causal paper comparison, with directionally consistent pressure/gradient signatures and consistent numerical/energy evidence.** No new numerical defect is established. The three-cell mirror tail described below is a real nonzero diagnostic result, compatible with the configured finite-rate mirror relaxation; it must remain visible.

Reviewed the final [metrics](metrics.json), [report](report.md), all four PNGs, the [signed pressure-band audit](pressure-episode/README.md), and the [mirror-tail audit](mirror-tail/README.md). Derived scalar/block values and provenance are retained alongside this review. This review uses existing metrics, images and the separately authorized one-snapshot audit; it makes no source change or simulation run.

## Adequacy and block dependence

The requested [6,18] interval is fully covered by **50 bracketing snapshots**, maximum gap 0.250384, with six duration-two blocks beginning at 6,8,10,12,14,16. Histories and force snapshots cover the same requested interval. No undefined metric, conflicting duplicate or abandoned-future sample is reported.

| Block | Solenoidal acceleration fraction | Pressure C | Pressure R | S RMS |
| --- | ---: | ---: | ---: | ---: |
| [6,8] | 0.6677 | -0.2018 | 0.8334 | 0.5402 |
| [8,10] | 0.6832 | -0.1595 | 0.8778 | 0.4971 |
| [10,12] | 0.7626 | -0.3733 | 0.6936 | 0.5333 |
| [12,14] | 0.4213 | +0.0010 | 0.9909 | 0.6242 |
| [14,16] | 0.1965 | -0.0419 | 0.9579 | 0.6588 |
| [16,18] | 0.6412 | -0.1364 | 0.8841 | 0.4321 |

The interval now includes both the compressive episode and a later return to more solenoidal forcing with reduced strain. These changes are part of the retained sample, not grounds for selecting a more favorable window. The realized forcing mixture is strongly variable even though the prescribed parameters did not change.

Empirical positive-sequence correlation times are approximately **1.126 for delta-B**, **1.654 for beta**, **0.707 for S RMS**, and **0.638 for parallel-velocity RMS**. Duration-two blocks are shorter than twice the first two estimates. The ACF is mean-subtracted, not detrended, and uses bracketing data resampled through approximately 17.75; its nominal effective-sample counts are not measured independent sample sizes.

Pairing adjacent duration-two means gives only **three duration-four blocks**, with no field reanalysis:

| Interval | C | R | S RMS | Solenoidal fraction |
| --- | ---: | ---: | ---: | ---: |
| [6,10] | -0.1807 | 0.8556 | 0.5186 | 0.6755 |
| [10,14] | -0.1861 | 0.8422 | 0.5788 | 0.5919 |
| [14,18] | -0.0892 | 0.9210 | 0.5454 | 0.4189 |

Their sample SD is 0.0545 for C and 0.0421 for R. All three longer intervals retain modest compensation, while its magnitude varies. Reblocking smooths the episode; it neither proves independence nor supplies enough blocks for precise uncertainty estimates. Known heating and finite-sample variability should be reported, not removed by endless extension. The predeclared progression [6,10]→[6,14]→[6,18] has now served its purpose.

## Interpretation of the four groups

**Magnetic strength and marginality: qualitatively compatible, precise residence comparison inconclusive.** Mean |B|/B0 is **1.2784**. Mean near-band fractions are **3.288% mirror** and **1.140% firehose**, with large block variability. Float32 strict fractions are **0.7925%** and **0.2209%**, respectively; most are ambiguous under the storage-rounding envelope. Double-history mirror-inclusive/strict-soft-rate means are about **1.5865%**, and firehose-inclusive **0.03757%**. These are different predicates and checkpoints; the firehose equality-counter sensitivity documented in the earlier review still applies.

Mirror crossings outside the rounding envelope are **not zero**: the time-mean volume fraction is **3.53247e-8**, entirely from **three of 1,769,472 cells** at **t=14.250230767343458**, instantaneous fraction **1.695421e-6**. The largest stored X is **1.0002245168**. The three X−1 values are 2.24517e-4, 2.13621e-4 and 4.02660e-6; excesses beyond their individual rounding envelopes are 2.33900e-5, 1.16724e-6 and 8.64381e-9. Two cells have weak B and local beta near 2965/3151, magnifying X's pressure-rounding sensitivity. The third has beta about 62.9.

The source audit confirms finite monotone mirror relaxation of positive excess by `1/(1+nu*dt)`, with no mirror upper hard wall when backups are off. A small residual above X=1 is therefore allowed by this model; the snapshot is after the post-STS rates/walls. The measured tail is compatible with finite soft relaxation, but the pre-relaxation state and exact source balance cannot be reconstructed from this float32 snapshot. No configured mirror hard wall was violated. This event is distinct from the earlier strict firehose restart defect. It is not legitimate to label these three cells rounding-only, to impose a new hard mirror criterion, or to infer continuous-time admissibility from snapshots. **No firehose crossing outside the rounding envelope is observed**, and sampled post-operator hard-volume fraction remains zero.

**Pressure balance: consistent with modest, variable compensation; not strong global balance.** Mean instantaneous C is **-0.15199**, R **0.87297**. These are time means of normalized diagnostics, not ratios formed from pooled variances. The perpendicular/magnetic pressure-variance ratio is **8.95**, and **93.0%** of perpendicular-pressure variance lies in the first two perpendicular shells, which also contain all parallel wavenumbers. Global statistics are consequently sensitive to large-scale forcing.

The signed-band audit resolves the apparent t=14 contradiction: forcing-band C=+0.517/R=1.120, while resolved smaller scales `3*pi<|k|<=24*pi` have C=-0.670/R=0.337. The resolved band remains anticorrelated at the inspected t=10,12,14 snapshots. This is signed evidence beyond matched auto-spectra, but only for those three instants, not a new time-averaged cross-spectrum. Positive global correlation during the episode does not show that pressure compensation disappeared throughout the resolved spectrum.

**Gradients: a persistent directional hierarchy, causality inconclusive.** Resolved strain/transverse-gradient power is **0.004826** overall; all six block ratios lie between **0.00354 and 0.00623**. The perpendicular-cutoff ratio is **0.08251**, so scale bands must remain separate. Integrated mean-subtracted local parallel velocity variance is **0.234464**, versus summed perpendicular variance **0.294152**. Their square roots are **0.48422** and **0.54236**, with parallel fluctuations contributing **44.35%**; these are square roots of time-mean variances, not means of RMS. The hierarchy coexists with substantial parallel motion. S, b.grad(u.b), and S−div(u) remain distinct diagnostics; the large difference between the first two includes curvature and finite-difference product-rule residual. There is no matched passive control establishing causation.

**Spectra and achieved regime: finite-time evidence, broad-cascade claims inconclusive.** Mean Mach proxy is **0.21366**, vector delta-B RMS/B0 **0.80664**, and volume-mean beta **9.2424**. The kinetic/magnetic spectra have comparable power away from the forcing peak but steeply roll off. Forcing maximum 3pi to the eight-cell guide 24pi spans less than a decade. Perpendicular shell sums, one resolution and centered derivatives cannot establish convergence or a broad inertial range.

Thermal-energy block means increase from **18.279 to 24.763**; beta-history means increase nonmonotonically from **7.860 to 11.043**. The duration-four Mach means decline 0.2327→0.2104→0.1978; this proxy also changes as thermal pressure heats, so its decline cannot alone establish weaker velocity fluctuations. Stationary thermal pressure is not required in this uncooled setup.

## Energy, provenance and final retention

Applied work is **7.68000000000012** over [6,18], total power **0.64000000000001**. Conserved-energy change minus work is **-5.720e-12**, relative **7.45e-13**, across both restarts. Mean instantaneous solenoidal acceleration fraction is **0.56208**, while the ratio of mean powers is **0.50229**. The latter's proximity to nominal 0.5 does not erase the episode or validate a forcing implementation by agreement: nominal 0.5 is a ratio of expected innovation powers, and neither measured statistic is an injected-energy partition. The deliberate dedt=.32 per-volume convention gives total .64 for this volume-two box.

All three corrected segments returned zero, ending at t=18/cycle32604, with numerical revision `71ad25ebce73d33db048defd8f585a7dce0528c2` and binary SHA `bd699cf711d0c2a21037690d28e8639b140f0ed96ad19c2d0406636c331527a3`. Retained lineage boundaries are t=10 and14, without duplicate or discarded-future data; prior aggregate manifests remain unchanged. Final analysis revision is `64d652317bec173d1d33e23906d0e93f21463409`, script SHA `a8241542d0d3cc5e4ecc83dd17b4cd7819467f9c2079eb707ebf231a5df812aa`.

Density/pressure-floor, nonfinite and nonpositive counters remain zero. About **2.922e9** unprojected LF hard-bound events correspond to **0.3759% of active cell-stage checks** during [6,18], not unstable volume or unique cells. Reserved `lf_hwproj` cannot count projections. Analysis Parseval/PDF normalization errors are about 4.19e-15/1.22e-12.

Retain the four final PNG/PDF groups, metrics/report, exact three-segment metadata and commands, input/binary/source provenance, earlier [6,10]/[6,14] reports, pressure-band and mirror-tail audits, and the original failed checkpoint/continuation as a negative control. The reusable guide and primary-reference definitions remain essential: [Squire Fig. 3](https://arxiv.org/html/2303.00468v2) uses X/2 and an active/passive comparison; [MKS24](https://arxiv.org/html/2405.02418v2)'s forcing/occupancy conventions differ and Fig.13b is a beta-100 limiter scan. No image-derived tolerance, quantitative paper reproduction, LF-coefficient validation, active/passive causality or resolution convergence is claimed.
