# Signed pressure covariance across the late forcing episode

**The sampled global pressure-correlation reversal is dominated by forcing-scale modes. Resolved smaller-scale compensation persists in all three inspected snapshots.** This narrows the interpretation of the [6,14] review: the positive late global correlation is not evidence that pressure compensation disappeared throughout the resolved spectrum, and does not by itself indicate a solver defect. It remains useful to complete the predeclared extension to characterize the evolving forcing-scale episode.

Only the existing snapshots nearest t=10,12,14 were inspected, at actual times **10**, **12.000347469555237**, **14**. These are instantaneous measurements, not a new time-averaging window. No source, input or reusable figure group was changed.

## Exact definitions

Use `a=p_perp-<p_perp>` and `b=B²/2-<B²/2>`, with uniform-cell volume means. The existing benchmark reader loads float32 fields into float64. With `A=fftn(a)/N` and `Bhat=fftn(b)/N`, each band has:

- `Va=sum |A|²`, `Vb=sum |Bhat|²`, signed covariance `Q=sum Re[A*conj(Bhat)]`;
- `C=Q/sqrt(Va*Vb)`;
- `R=sum |A+Bhat|²/(Va+Vb)=1+2Q/(Va+Vb)`.

Both conjugate modes are included. These are integrated powers/covariances in pressure squared; there is no division by shell width or mode count. The masks use the full physical three-dimensional wavevector, including pure parallel modes, rather than perpendicular shells:

| Band | Definition |
| --- | --- |
| Forcing | `0<|k|<=3*pi`; the box's smallest nonzero magnitude is pi, so this contains the imposed pi..3pi shell (38 full FFT modes). |
| Resolved smaller scales | `3*pi<|k|<=min_i(pi*N_i/L_i)/4=24*pi`; this is called `resolved_subforcing` in the JSON and denotes smaller wavelengths than forcing. |
| Remaining high k | `|k|>24*pi`, through every represented FFT mode, including high parallel modes at small perpendicular k. |

The three bands partition all nonzero modes. The result also retains their full-resolved union, all nonzero modes, full grid, and the negligible DC residual from mean subtraction. Bounds follow the reusable analyzer's physical grids and Nyquist/4 convention.

## Results

| Time | Band | Signed covariance Q | C | R |
| --- | --- | ---: | ---: | ---: |
| 10 | Forcing | -0.0043684 | -0.17150 | 0.88721 |
| 10 | Resolved smaller scales | -0.0113964 | -0.66109 | 0.37612 |
| 10 | Remaining high k | -0.0002190 | -0.96236 | 0.03774 |
| 10 | All nonzero | -0.0159838 | -0.35440 | 0.72069 |
| 12.000347 | Forcing | +0.0012770 | +0.02904 | 1.01545 |
| 12.000347 | Resolved smaller scales | -0.0108884 | -0.70583 | 0.30821 |
| 12.000347 | Remaining high k | -0.0001290 | -0.95381 | 0.04641 |
| 12.000347 | All nonzero | -0.0097403 | -0.14747 | 0.90113 |
| 14 | Forcing | +0.0581913 | +0.51704 | 1.12017 |
| 14 | Resolved smaller scales | -0.0215815 | -0.66962 | 0.33736 |
| 14 | Remaining high k | -0.0006166 | -0.95484 | 0.04516 |
| 14 | All nonzero | +0.0359932 | +0.17669 | 1.06956 |

Forcing-band perpendicular-pressure variance grows **0.067908→0.152602→0.955189**, whereas its magnetic-pressure variance is **0.009554→0.012673→0.013261**. At t=14, forcing-band positive covariance exceeds the global covariance; the remaining modes partially cancel it. The resolved smaller-scale covariance is more negative at t=14 than at t=10, with C staying near -0.67 and R around one third. This signed result supplies evidence that matched auto-spectra alone could not provide.

The strong high-k anticorrelation is recorded but is not physical convergence evidence: this band lies beyond the conservative resolved guide, has small variance, and can be affected by numerical dissipation. The resolved band itself spans only a factor eight in wavenumber and is a band-integrated measurement, not proof of per-mode cancellation or an inertial range. Three snapshots do not establish a stationary ensemble or a causal attribution to the forcing; they show where the measured covariance resides.

## Verification and provenance

[results.json](results.json) retains all band variances, residual variances, covariances, correlations, mode counts, exact definitions and source/file hashes. Full Fourier sums reproduce the stored snapshot Va, Vb, Q, C and R, and the direct real-space computation independently matches the same values. The disjoint-band variances/covariances sum to the all-nonzero values. The verification bound is solely a float64 arithmetic check (`256*epsilon*max(1,|reference|)`), not a new physical acceptance tolerance; actual errors are in the JSON and [derived-values.json](derived-values.json).

The single CPU-only step used eight allocated cores, zero GPUs and `--exact --overlap` inside allocation **5629677**, concurrent with the unchanged simulation. It completed successfully. [run.sh](run.sh), [run.log](run.log), and [analyze_pressure_episode.py](analyze_pressure_episode.py) preserve the command and script. The current analyzer SHA and each snapshot hash were checked against the retained [6,14] metrics before computing. There was no numerical application run, parameter change, paper-image tolerance or additional figure group.
