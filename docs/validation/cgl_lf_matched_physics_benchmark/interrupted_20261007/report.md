# Interrupted matched CGL-LF experiment, 7 October 2026

**This is a reviewed diagnostic of a failed/interrupted experiment, not a completed physical benchmark or a reference for statistical acceptance.** The passive member failed during a pre-LF RKL stage after t≈2.5077. The user then stopped the active member; its last history is t=6.16005 and last field snapshot t=6.00019. Neither member supplies the planned matched [6,14] average.

Both members used 192×192×384, PPM4/HLLE, RK2/RKL2, conservative A, seed 271828, and dedt=0.32 per volume (total 0.64). The box has L=(1,1,2), initial rho=B0=1 and p_parallel=p_perp=5 (beta0=10); time is in L_perp/v_A0. Closure k_L=2π, nu_coll=0, finite limiter rate 10¹⁰ and thresholds X=−2,+1 were unchanged. The passive flow is isothermal MHD with c_iso=√5 and diagnostic CGL pressures. Physics, closure, forcing prescription and sampling were matched; realized adaptive timesteps and acceleration amplitudes differ. Numerical source revision: `ab9b543e7e6a972ebfc026d0c98ea3dfee9cb55b`. Simulation and analysis revisions are recorded separately in [provenance.json](provenance.json).

## Figures and sampling

| Diagnostic | Raster | Vector |
| --- | --- | --- |
| Shared startup PDFs | [PNG](paired_pdfs.png) | [PDF](paired_pdfs.pdf) |
| Shared startup threshold occupancy | [PNG](paired_occupancy.png) | [PDF](paired_occupancy.pdf) |
| Shared startup energy, pressure and strain spectra | [PNG](paired_spectra.png) | [PDF](paired_spectra.pdf) |
| Shared startup signed pressure balance versus scale | [PNG](paired_pressure_balance_scale.png) | [PDF](paired_pressure_balance_scale.pdf) |
| Later active energy, spectra and achieved regime | [PNG](active_2to6_spectra_energy.png) | [PDF](active_2to6_spectra_energy.pdf) |
| Later active flow geometry | [PNG](active_2to6_gradients.png) | [PDF](active_2to6_gradients.pdf) |
| Later active pressure spectra and whole-box statistics | [PNG](active_2to6_pressure_balance.png) | [PDF](active_2to6_pressure_balance.pdf) |
| Passive pressure minima and spatial variability | [PNG](passive_pressure_tails.png) | [PDF](passive_pressure_tails.pdf) |

The paired interval is **[0.5,2.5]**, using ten bracketing snapshots per member and four 0.5-unit blocks. The later active interval is **[2,6]**, using 18 bracketing snapshots and four 1-unit blocks. Snapshot cadence is 0.25; histories are 0.02. Diagnostics are interpolated to interval endpoints between bracketing outputs; fields are not interpolated. Bands show pointwise minimum/maximum of time-block means. These short blocks are correlated relative to forcing tcorr=2; they are not independent samples or confidence intervals. For example, active startup Mach proxy block means rise 0.098→0.132→0.207→0.257.

## Why the shared startup anisotropy PDF looks this way

**The threshold spikes are compatible with the prescribed stiff limiter, but the strongly positive, mirror-dominated active distribution is not yet agreement with the published steady-state result.** The [right panel of Shared startup PDFs](paired_pdfs.png) combines a large startup transient with accumulation near the limiter boundaries. Its shape should not be presented as a validated, universal turbulence distribution.

### What is actually plotted

The abscissa is the local, instantaneous

$$X=\frac{2(p_\perp-p_\parallel)}{|\mathbf B|^2},$$

using the cell's magnetic strength, not $B_0$ or an average beta. The ordinate is probability **density per unit $X$**, with volume weighting and subsequent trapezoidal time averaging over $[0.5,2.5]$:

$$P_i=\left\langle\frac{\sum_{c:\,X_c\in[e_i,e_{i+1})}V_c}{V\,(e_{i+1}-e_i)}\right\rangle_t,
\qquad \sum_i P_i\,\Delta X_i=1.$$

Both runs use the same 128 linear bins, with $\Delta X=0.023438$, spanning all retained values. There is no clipping to $[-2,1]$ by the analysis. The PDF integrals differ from one by at most $1.24\times10^{-12}$. A density above one is normal: the last bin has height 9.41 for active and 12.25 for passive, but contains respectively **22.06% and 28.71%** of the time-volume sample. That bin spans approximately $[0.97660,1.00004]$; plotting at its center, $X=0.98832$, does not mean the mirror threshold moved. The log ordinate emphasizes small populations and cannot be read as probability per plotted vertical interval.

| Measured population over [0.5,2.5] | Active | Passive |
| --- | ---: | ---: |
| Last histogram bin, near the mirror boundary | 22.06% | 28.71% |
| First histogram bin, near the firehose boundary | 0.0100% | 3.430% |
| Mirror band, $\lvert X-1\rvert\le0.05$ | 25.73% | 32.10% |
| Firehose band, $\lvert X+2\rvert\le0.05$ | 0.0146% | 4.971% |
| Outside both near-threshold bands | 74.26% | 62.92% |

Endpoint-bin populations and near-band populations are different measurements; neither is the fraction strictly beyond a threshold. The large endpoint densities coexist with substantial interior volume. Bin width changes peak heights, so use integrated fractions when comparing plots. These values and their input hashes are retained in [startup-x-pdf-audit.json](startup-x-pdf-audit.json).

### Why there are boundary peaks and a positive bias

The actual finite-rate limiter relaxes excess anisotropy toward a threshold, not all the way to isotropy. With $d=p_\perp-p_\parallel$ and crossed threshold $d_*$, its update is

$$d_{\rm new}=d_*+\frac{d_{\rm old}-d_*}{1+\nu_{\rm lim}\Delta t}.$$

Here $\nu_{\rm lim}=10^{10}$ and typical outer timesteps are order $10^{-4}$, so an overshoot is compressed into a very narrow region beside the threshold at output time. Continued driving can consequently accumulate a large population there. The unconditional fluid firehose wall also enforces $X\ge-2$ at scheduled wall checkpoints; disabling backups does not remove that wall. Intermediate LF stages may cross it, and float32 output can straddle it by rounding. This makes the two numerical boundaries asymmetric. With background collisions zero, stable interior states have no background relaxation toward $X=0$. These are effects of the prescribed fluid closure, not resolved kinetic instability structures. In passive runs the same limiter acts on the evolved tracer pressures despite their absence from the momentum force.

The **positive-side bias is demonstrably associated with an evolving startup sample**. Successive half-unit blocks have the following mirror-band fractions:

| Time block | Active | Passive |
| --- | ---: | ---: |
| [0.5,1.0] | 58.01% | 62.62% |
| [1.0,1.5] | 24.17% | 28.54% |
| [1.5,2.0] | 9.49% | 20.18% |
| [2.0,2.5] | 11.25% | 17.07% |

At the individual $t\simeq0.5$ snapshot, **93.89% of the active box** is mirror-near. This early population contributes strongly to the averaged spike. Initialization is isotropic, at rest, with uniform $B$; magnetic fluctuations subsequently grow. The active block-mean $\langle|B|\rangle/B_0$ rises from 1.116 to 1.200 across these blocks. Startup field amplification and the closer positive boundary ($+1$ versus $-2$) provide a plausible physical explanation for the skew, but this is an **inference**, not a measured decomposition of its cause. In compressible CGL, near isotropy and before heat fluxes/scattering,

$$\frac{D(p_\perp-p_\parallel)}{Dt}\simeq p\left(3S_\parallel-\nabla\!\cdot\mathbf u\right),\qquad
S_\parallel=b_i b_j\partial_j u_i.$$

Thus field amplification alone does not determine the pressure PDF: density changes, heat transport, the existing mean anisotropy and limiter history also matter. The marginal PDF cannot establish their individual contributions. The entire averaging interval is only one forcing correlation time; the shaded block ranges quantify drift as well as fluctuations. For context, the later active-only [2,6] sample has mirror/firehose near fractions **4.53%/0.215%**, substantially below startup. That does not supply a matched passive comparison or establish stationarity.

### Is this consistent with the papers?

**Partly in mechanism and shape; not established as a matched physical result.** [MKS24 §3.8 and Figure 13b](https://arxiv.org/html/2405.02418v2#S3.SS8) show a rising active interior and threshold spikes for the $10^{10}$ limiter, qualitatively resembling this plot. A central maximum is therefore not mandatory. That figure is an **Alfvénically forced $\beta=100$ limiter scan**, not this $\beta_0=10$ mixed-forcing experiment. Its $\beta\Delta$ corresponds exactly to our $X$ when $\beta$ uses local $p_\parallel$ and $\Delta=p_\perp/p_\parallel-1$. [MKS24 §3.3 and Figure 2b](https://arxiv.org/html/2405.02418v2#S3.SS3) identify very large first-crossing-time occupancy as undeveloped turbulence; their statistical windows last at least two time units and start after $t=6$. Our startup window cannot test that steady-state result. Their threshold-exceedance diagnostic also must not be equated with our finite-width near bands or post-operator, float32 snapshot crossings.

[Squire 2023 §4.1 and Figures 3–4](https://arxiv.org/html/2303.00468v2#S4.SS1) describe active interiors concentrated near zero, passive interiors much flatter, and markedly reduced active threshold populations. Our active interior's strong positive bias is a real difference from that description, not something to dismiss because endpoint peaks are allowed. Their anisotropy coordinate is $Y=4\pi\Delta p/B_{\rm physical}^2=X/2$ in our magnetic units, so comparison requires $P_X(X)=P_Y(X/2)/2$ and thresholds $Y=-1,+1/2$. Figure 4 also includes a different firehose limiter; it should not be overlaid without accounting for it.

Forcing projection, total injected power (ours is 0.64), changing beta and limiter/output timing further limit direct comparison. The figure's sharply bounded shape is **consistent with the configured limiter**, the strong temporal evolution makes the steady-state PDF comparison **inconclusive**, and the separate reproduced LF failure remains **concerning**. After the numerical issue is resolved, assess any persistent late-time mirror bias against the relevant reference regime and temporal convergence; neither a rising interior alone nor resemblance to Figure 13b establishes success or failure. Check definitions and sampling first, achieved forcing/regime second, resolution third, and implementation as indicated by the numerical evidence. Do not modify thresholds or forcing merely to obtain a desired PDF.

## Physical evidence

**Consistent directional signatures, with insufficient evidence for stationary validation.** The active startup distribution of magnetic strength is narrower, field-aligned strain is lower, and parallel velocities remain substantial. These are the contrasts motivated by [Squire 2023, Figure 3 and §3.1.3](https://arxiv.org/html/2303.00468v2). The numerical failure prevents promoting them to an accepted turbulence reference.

| Time-averaged startup quantity | Active | Passive |
| --- | ---: | ---: |
| Magnetic-strength distribution width | 0.11785 | 0.25682 |
| delta B_rms/B0 | 0.60444 | 0.74354 |
| Thermal Mach proxy | 0.17326 | 0.17953 |
| Parallel velocity power fraction | 0.5472 | 0.52507 |
| S_parallel rms | 0.39895 | 1.0226 |
| Near mirror fraction | 0.25729 | 0.32105 |
| Near firehose fraction | 0.00014565 | 0.04971 |

The magnetic-strength width is the standard deviation of the time-volume pooled |B|/B0 distribution, computed from exact field moments; it differs from the vector fluctuation delta B_rms/B0. Active/passive width ratio is 0.459. Resolved S_parallel/perpendicular-gradient power ratios are 0.0063644 and 0.041956; at the cutoff these become 0.09862 and 0.1053. Small-scale numerical behavior must not be equated with the resolved contrast. S_parallel=b_i b_j ∂_j u_i differs from b·∇(u·b); the latter includes field curvature. The supplement also retains S_parallel−∇·u, the ideal compressible induction proxy.

Near-threshold fractions mean |X−1|≤0.05 and |X+2|≤0.05, where X=2(p_perp−p_parallel)/B². They are separate from strict X>1/X<−2 counts. Almost all strict snapshot crossings lie within float32 rounding envelopes; robust firehose exceedance is zero in both startup averages, while robust mirror fractions are only order 10⁻⁸. Inclusive history switches differ again. None of these snapshot distinctions excuses the later pressure-floor event.

**Pressure compensation is clearest as a function of scale.** At k_perp/(2π)=8.5, startup C/R are −0.964/0.0365 active and −0.836/0.174 passive. In the later active interval, whole-box C=-0.253 and R=0.779, whereas that same spectral shell has C=−0.941 and R=0.0590. The full-|k| resolved filter gives essentially the same shell values. Here C=Re(P_ab)/√(P_aa P_bb), R=(P_aa+P_bb+2Re(P_ab))/(P_aa+P_bb), a=δp_perp, b=δ(B²/2), with powers averaged before ratios. This supplies phase-sensitive evidence beyond matching amplitudes, relevant to [MKS24 Figure 6](https://arxiv.org/html/2405.02418v2). Passive CGL p_perp is diagnostic; the actual passive momentum pressure is 5rho and its spectrum is labeled separately.

**A broad inertial range remains inconclusive.** Later active E_K k_perp^(5/3), normalized at k_perp/(2π)=2.5, is 1.000,0.737,0.717,0.398,0.216 at shell centers 2.5,4.5,8.5,16.5,23.5. There is a short approximately level interval followed by substantial roll-off, not a wide established cascade. The startup pair is much steeper and still developing. Curves are shell sums divided by Δk, summing all k_parallel; the zero transverse shell is omitted in paired log plots. The eight-cell guide and shaded high-k band are heuristic scale markers, not measured convergence limits. No slope tolerance is extracted from paper images. The old 96 PLM run also has a different averaging epoch, so this is not a controlled resolution/reconstruction comparison.

Later active averages are Mach proxy 0.268 (block SD 0.034), delta B_rms/B0=0.759 (0.094), and volume-mean beta=7.357 (0.237). The proxy uses c_th²=(5/3)〈p_iso〉/〈rho〉; it is not a CGL wave speed. Passive startup isothermal Mach is 0.240, using its fixed dynamical sound speed. No forcing retuning was performed.

## Numerical integrity and energy

**Concerning: the passive evolution failed.** The fatal condition was one pressure-floor refresh event at pre-sweep stage 13/15, with finite raw pressure becoming negative. The 2254 reported hard-wall crossings were intermediate stage activity, not the reason for the strict stop. Floor counters include refreshed halos, while nonpositive/hard-bound checks inspect active cells; post-refresh nonpositive=0 is not proof that the raw state was admissible. [Retained forensic evidence](/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/floor-diagnosis/evidence.md) records the original failure and checkpoint.

All saved primitive states remained positive. At passive t=2.5003003, global minimum p_parallel=2.5140, minimum p_perp=3.7129 and minimum U=5.2292, with mean U=8.5792 and spatial SD 1.3292. Across all 11 outputs the minimum U is 4.589. A [diagnostic replay](/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched/floor-diagnosis/reconstruction-same-settings.log) reproduced the event at cycle 5800, t=2.507714514270586, with raw p_perp=−0.00875815. Its LF entry cell had p_parallel≈0.1573 and p_perp≈0.04455: this much colder state was not a cold tail already visible in the saved active cells at t=2.5003. Its formation between outputs needs implementation/stability diagnosis; the snapshot audit alone does not identify the cause.

All retained floor/nonfinite/nonpositive history counters are zero before the fatal step, and retained strict hard-wall volume is zero. Active histories also show zero such counters through user interruption. Analysis Parseval errors are at most 7.3×10⁻¹⁵ and PDF normalization errors at most 1.3×10⁻¹²; these validate normalization, not the physical implementation.

Measured forcing work supplies total power 0.64 in both startup intervals. Active ΔE−ΔW residuals are 3.11e-12 over [0.5,2.5] and -1.25e-11 over [2,6]. These use the actual conserved-energy/forcing ledgers. Passive thermal energy rises 15→17.1583 by the last output; KE+ME+passive U does not obey the active conserved-energy budget. Passive pressure-work diagnostics are diagnostic thermal-stress contractions, not forces applied to passive momentum.

Equal expected forcing innovation power is prescribed; measured instantaneous solenoidal acceleration fractions agree between the matched saved frames to 1.4e-11. Their time means are 0.64093/0.64094, while ratios of integrated powers are 0.62964/0.61651 because forcing amplitudes depend on each flow. Neither statistic is an injected-energy partition.

The active launch was canceled before its metadata finalizer, leaving empty history inventories and no return code. An [analysis-only metadata copy](active-analysis-metadata.json) lists the actual retained histories without inventing successful completion or editing raw metadata. [The amendment record](history-inventory-amendment.json) hashes the original/corrected cached metrics and confirms unchanged snapshots, FFTs, PDFs, scalar measurements and forcing products. The complete cached metrics and original backups remain outside Git.

## Cost and reproduction

Allocation 5631202 reserved 16 nodes for 3926 s: **17.449 node-hours**, or 139.591 logical GPU-device-hours at 8 devices/node, including idle reservation after passive failure. Canonical steps used 8 nodes each for 1221 s passive and 3922 s active; a separate diagnostic replay used 8 nodes for 332 s. Allocation and step costs overlap and must not be added. The queued follow-up was canceled before starting. See [costs.json](costs.json).

The paired analysis took 1332.64 s and the later active analysis 1145.75 s, concurrently on one CPU node, with maximum process RSS about 9.36 and 8.41 GiB. Pressure-tail inspection took 21.48 s and history amendment/replot 36.10 s. No GPU computation was used for analysis.

Run on an allocated compute node, using the recorded environment and source revision:

```bash
S=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2
M=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/matched
O="$M/interrupted-analysis"
PY=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark/dbfplot-refresh/venv/bin/python
mkdir -p "$O/reproduction-runtime"
cd "$O/reproduction-runtime"
source "$M/../scripts/runtime_cpu.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$PY" "$S/scripts/compare_cgl_lf_physics_benchmark.py" "$M/active-to14" "$M/passive-to14" --time-start 0.5 --time-end 2.5 --block-duration 0.5 --output-dir "$O/paired-startup"
"$PY" "$S/scripts/analyze_cgl_lf_physics_benchmark.py" "$M/active-to14" --time-start 2 --time-end 6 --block-duration 1 --output-dir "$O/active-interrupted"
"$PY" "$O/restore_interrupted_history_inventory.py"
"$PY" "$O/inspect_passive_pressure_tails.py"
```

The retained amendment is specific to this interrupted launch and runs only after both analysis commands complete; original execution records and raw metadata remain in O. Full command arrays, environment, source hashes, snapshot provenance and machine-readable measurements are in [metrics.json](metrics.json), [provenance.json](provenance.json), and the external analysis plan. The [benchmark guide](../../../cgl_lf_matched_physics_benchmark.md) controls future canonical launches; this failed/interrupted directory must not be promoted to a statistical reference. Existing LF damping checks and the accepted B4 sharp-contact limitation remain separate.
