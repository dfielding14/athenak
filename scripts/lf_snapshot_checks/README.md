# Offline active CGL closure checks

Run this bundle on existing developed active $\beta_{i0}=16$ snapshots. It compares the nominal fixed-$k_L$ closure with a finite-reach Snyder–Hammett–Dorland (SHD) field-line integral, an eight-pole approximation, and path-dependent scattering for positive background collision rates. It operates on frozen states. Output quantities describe the instantaneous heat-flux contribution to pressure evolution.

## Cluster quick start

Use a compute node with a **CPU** allocation and a Python environment with NumPy, SciPy, Numba, and Matplotlib. There is no MPI launch or GPU kernel in the analyzer. Use one process per snapshot and set Numba threads to the CPUs assigned to it. Local verification used NumPy 2.2.6, SciPy 1.15.2, Numba 0.61.2, and Matplotlib 3.10.8. The bundle is self-contained. The standard install command on a host with package access is `python -m pip install -r requirements.txt`.

1. Select full-volume native AthenaK 1.1 `mhd_w_bcc` snapshots. A shared file or a **quoted glob containing every rank at one output time** works. Use the same resolved mesh and comparable developed turnover times in both runs, preferably at least two times per run. Use clean output directories. The reader verifies active CGL, periodic boundaries, full coverage and a uniform leaf level, including uniformly refined meshes. For a restart-only archive, see below.
2. Set `LA` and `ELL` to **physical code lengths**, with the convention documented. Use the matched MHD control's measured Alfvén scale, e.g. the separation where the previously adopted velocity increment equals the reference Alfvén speed. Use the measured parallel viscous scale for each analysis when available. If no cutoff can be identified in the zero-background run, use the paired collisional cutoff as an explicitly labeled comparison scale. If only a design target is available, label it as such. Convert box fractions to code lengths. Distinguish $1/k$ from $2\pi/k$. The code never infers a viscous scale from the input collision rate.
3. Stage one snapshot at a time. The assembled six float32 fields need **24 bytes per cell**, about 10.1 GiB at $768^3$ or 24 GiB at $1024^3$, plus working memory. `--cache-dir` places that array in a private temporary directory on scratch and removes it after the run. Prefer node-local scratch with room for the array, sufficient RAM for its working set, and enough CPUs. Random rays into a disk-backed array larger than available page cache can be very slow. This is a single-node postprocessor; choose a manageable developed snapshot. The code accepts unsliced 3D periodic output and rejects mixed refinement levels.
4. From the AthenaK repository root, follow [the GPU-agent handoff](GPU_AGENT.txt) for the pull, per-snapshot runs, and final comparison commands. Select fresh per-snapshot output names for each pair of times, and keep the combined results on persistent storage.

The default sample is 64 common anchors across separations and methods. The optical calculation uses the first eight of those anchors and all required neighboring faces. Increase `--anchors` and `--optical-anchors` and repeat with a second `--seed` if scatter is comparable to the differences. `--optical-anchors 0` skips optical work for an initial runtime estimate. Include it in the final collisional comparison by using a positive count.

## What to compare

The two **measurement scales** $L_A$ and $\ell_\nu$ set the separation between pairs of cell centers. At $0.5$, $1$, and $2$ times each scale, the script samples straight periodic chords parallel to the field at each anchor and at random azimuth perpendicular to it. The intervals describe local-field structure functions. They are distinct from a Fourier band or a field-line-traced separation. Each row includes its realized median separation and angle. Treat rows below four coarse-direction cell spacings, longer than half the shortest box side, or with zero-length pairs as flagged. Near the resolution threshold, inspect the recorded separation distribution and rounding angles.

The **integration reach** $R$ sets how far each heat-flux calculation traces in each direction from a face. By default the runner includes both measurement lengths and $0.25L$, $0.5L$, $L$, $2L$, with $L$ the shortest box side. These controls are independent. Override with `--reaches` in units of $L$, e.g. `--reaches .125 .25 .5 1 2 4 8`. A finite trace is always a comparison reference. If the last doubling still changes a result appreciably, report it as reach-sensitive. Repeat with longer traces if feasible. Both step sizes, by default half and one-quarter of the smallest cell spacing, use identical faces and reaches. `--step-cells` changes the coarse step.

For each endpoint the runner reconstructs a single pair of signed fluxes on each of six shared faces. At fixed density and field, set $\Delta p=p_\perp-p_\parallel$. The diagnostic is
$$
\dot E_q=-\nabla\cdot[\boldsymbol b(q_\perp+q_\parallel/2)],\qquad
\dot\mu_q=-\nabla\cdot[\boldsymbol b q_\perp/B],\qquad
\dot{\Delta p}_q=3 B\dot\mu_q-2\dot E_q .
$$
With $\delta_\ell$ denoting an endpoint difference, the **signed transfer** is $\langle\delta_\ell\Delta p\,\delta_\ell\dot{\Delta p}_q\rangle$. At frozen pair geometry this is half the LF contribution to the time derivative of the pressure-anisotropy second-order structure function. Negative values mean that heat transport reduces those increments. Also compare the **relative RMS difference** in $\delta_\ell\dot{\Delta p}_q$, sign changes, and signed heat-flux increment differences. Pointwise disagreement can cancel in the signed transfer. Batch and paired standard errors are descriptive sampling-scatter estimates. A second snapshot or independent seed provides a stronger repeatability check.

| CSV family | Comparison |
|---|---|
| receiver | Nominal fixed-$k_L$ and short SHD reaches against the longest receiver-frozen SHD reach, with the source run's face scattering |
| uniform kernel | Eight positive fitted poles compared against SHD **at the identical reach and uniform inverse scattering length** for each channel. This isolates kernel fit error. Both omit intermittent scattering in this control. |
| ray step | Coarse versus fine step at the longest reach, on the same faces |
| path scattering | Path-integrated optical depth versus receiver-frozen scattering at equal reach on matched faces, only for positive background scattering |
| optical reach and step | Half reach and coarser optical step against full reach at finer step |

The rational filter evaluates exponential kernels on the same finite rays. It is not a benchmark of a global elliptic solve. The default angular fit band $kL$ runs from $0.1$ to the axial grid Nyquist; the per-channel sampled fit error is saved. For a band-sensitivity check use `--fit-band KMIN_L KMAX_L`; keep physical reach and faces fixed. Inspect the softest pole's decay length $1/\kappa_{\min}$ relative to the longest trace. For optical comparisons the default full reach is $0.5L$; use `--optical-reach` to extend it separately and inspect the half/full comparison. The $\nu_{\rm coll}=0$ run can have local limiter scattering, but a zero background invalidates the strictly positive optical-depth coordinate.

Both channels include the free-streaming cap after combining the kernel contributions. The perpendicular channel includes the receiving-face magnetic drive proportional to $(1-p_\perp/p_\parallel)\nabla_\parallel\ln B$. The trace calculation holds thermal prefactors at the receiving face. The optical variant varies the rate along the path, including finite soft/backup branches. These are specified nonlinear extensions of the uniform-coefficient closure. Native binary output supplies cell-centered magnetic fields; adjacent cell averages approximate the unavailable staggered normal fields. Consequently the fixed-$k_L$ result is a reconstructed baseline.

## Outputs to return

Return the combined `scales.csv`, `scale_comparison.png` and `scale_transfer.png`, with `metadata.json` from each snapshot. The two plots show local-parallel results, including a transfer curve around each target scale. The CSV also has perpendicular chords, half/double separations, actual lengths, sample counts, reach/step checks, paired uncertainties and resolution flags. Per-snapshot `samples.npz` retains sampled pairs, faces, fluxes and tendencies for follow-up. Report each scale's source, snapshot times, case settings, and any under-resolved or unconverged comparisons. CPU postprocessing timings cannot establish production GPU solver speed.

## Restart-only inputs

Use a compatible AthenaK executable, the source rank/file layout as required, and a **fresh output directory**. With a pre-existing full `mhd_w_bcc` bin block, e.g. `output2`, request a zero-cycle continuation:

```bash
./athena -r "$RST" -d "$NEW_EXPORT_DIR" time/nlim=0 \
  output2/file_number=0 output2/single_file_per_rank=false
```

If no such block exists, use `-i restart_to_bin.athinput` with an otherwise unused output number. Run AthenaK with the site's usual launcher for the required ranks. Finalization may write other configured outputs into that fresh directory. Verify zero meshblock-cycles and that the binary time and cycle equal the checkpoint. The existing-block route reproduced all nine primitive and magnetic arrays bitwise in a local active cooling restart export.

## Smoke test and analytic checks

```bash
python smoke.py /path/to/new/smoke-output
```

The included $24^3$ native pilots are solely for portability and regression testing. The two example lengths in the smoke script are test separations rather than measured turbulence scales. The smoke test exercises both collision cases, the plots, the rank-sharded reader, and signed fluxes and tendencies independently saved from full-grid pilot calculations. `benchmarks/` contains the original standalone helical, scattering-barrier, and shifted-solve checks, runnable individually. The benchmarks write their JSON results beside their scripts.

Physics references: [Snyder, Hammett & Dorland, §§VI–VII](https://w3.pppl.gov/~hammett/gyrofluid/papers/snyder_pop97.pdf) and [Dimits, Joseph & Umansky](https://doi.org/10.1063/1.4876617).
