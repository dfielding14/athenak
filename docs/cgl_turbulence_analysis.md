# Snapshot turbulence comparisons

`scripts/analyze_cgl_turbulence.py` calculates PDFs, energy/pressure spectra and
local-field strain diagnostics once per snapshot. `scripts/plot_cgl_turbulence.py`
replots six figure groups from small NPZ caches. They use the existing Athena
binary reader, NumPy/SciPy and `~/dbfplot`. The Frontier environment expects the
checkout at `~/athenak-cgl`; there are no new dependencies.

```bash
python scripts/analyze_cgl_turbulence.py --self-check
python scripts/plot_cgl_turbulence.py --self-check
python scripts/analyze_cgl_turbulence.py manifest.json
python scripts/plot_cgl_turbulence.py manifest.json --time 14 --output figures_t14
python scripts/plot_cgl_turbulence.py manifest.json --window 6.25 14 --output figures_window
```

Use one process per compute node for large fields. Under Slurm, `SLURM_PROCID`
and `SLURM_NTASKS` partition independent snapshots; `SLURM_CPUS_PER_TASK` sets
FFT threads. Cached results are skipped when source size/mtime and analyzer hash
match. Batch scripts should retain the exact analyzer and manifest they run.
Do not run two writers for the same cache simultaneously. Plotting does not read
large simulation fields.

## Manifest

The calculation reads `tasks`, each containing absolute `source` and `cache`
paths and a nominal output `time`. Plotting reads `cases`, each containing
`name`, `label`, snapshot NPZ `files`, and one `history` NPZ with `ht`, `mirror`,
`fire`, `backup` in percent. Relative plot paths are resolved against the
manifest. Optional `reference`, `k_label`, `kmax`, and `dissipation_start` control
comparison figures. Keep case files on a common time grid. A changing resolution
requires separate windows, since unlike Fourier grids cannot be averaged.
Different resolutions can be plotted together when their physical low-k shells
match; ratios stop at the common cutoff. Average each resolution separately.
Keep comparisons below the coarser axis Nyquist limit: the 384³/768³ figures
stop at kLx/(2 pi)=128, and quantitative band comparisons use 4–32.

The short discovery scripts and actual campaign manifests are retained in
[CGL/analysis/turbulence_20261017](/lustre/orion/ast207/proj-shared/dfielding/CGL/analysis/turbulence_20261017/README.md).
They list the accepted continuation chains, actual output times and missing
cadence slots. Regenerating a manifest picks up new outputs; rerunning analysis
calculates only uncached snapshots.

## Definitions and checks

- Assemble complete same-level blocks, checking logical tiling and physical
  bounds. Uniform AMR leaves use their true resolution, not root-grid Nx.
  Mixed-level grids are rejected rather than silently interpolated.
- Fourier shells use physical wavenumbers and are labeled `k Lx/(2 pi)`.
  Shell width is `Lx/max(Lx,Ly,Lz)`; spectra are densities in that coordinate.
  Real-FFT multiplicities and even-grid Nyquist derivatives are explicit.
  Parseval checks run on every spectrum.
- Kinetic spectra are half the variance spectrum of `sqrt(rho) u`; magnetic
  spectra are half the variance spectrum of cell-centered B. Their means are
  removed, including the guide field. Pressure spectra are variance densities
  without the half. These integrals differ from total history energies.
- Differentiate u first, then project its gradient using local b. This gives
  `S=b_i b_j partial_j u_i`, `Q=S-div(u)/3` and four local parallel/perpendicular
  gradient components without differentiating b. The four component powers,
  including their zero modes, add after integration, not shell by shell.
  S/Q and local `u_parallel=u.b` spectra remove their means.
- B PDFs are densities per `d(|B|/B0)` despite logarithmic bins/axes.
  Beta-anisotropy is `2(p_perp-p_parallel)/B^2`; probability outside the PDF
  range [-3,3] is recorded rather than renormalized away.
- Trapezoidal time averages require matched times and both requested endpoints.
  Temporal scatter is descriptive, not an independent-sampling confidence interval.

Self-checks cover rectangular Fourier modes, Parseval, odd/even transforms,
four isolated gradient projections, varying-field integrated closure, refined
block tiling, mixed-level rejection and uneven-time averaging. Actual field
analysis checks positive pressures/density, finite fields, nonzero local B,
spectral normalization and projected-gradient closure.

The full-grid CPU FFT route is intended for the present uniform grids through
768³, one process per node. The memory estimate rejects 1536³ before loading its
fields; that size needs distributed FFTs or an explicitly chosen coarsening
policy. No downsampling is used in the present analysis. Job `5697979` measured
275/282 seconds and 66.73 GiB process peak RSS on two actual 768³ snapshots,
using one process per node with 16 CPU threads; both passed all spectral checks.
Job `5697978` processed the 110 smaller Majeski snapshots in 4m07s on four nodes.

## Current evidence

The [completed Majeski comparison](/lustre/orion/ast207/proj-shared/dfielding/CGL/cgl_majeski_beta10_matched_20261016/analysis_20261017/RESULTS.md)
retains 110 canonical snapshots, full histories, endpoint figures and a
29-time [6.25,14] field average. Missing restart-boundary snapshots are explicit;
no field is synthesized. This pair used frozen pre-fix source `795c2e180c`,
with active cap 0 and passive cap 32. New passive runs default to cap 0;
the completed figures retain their original settings and provenance. The
[corrected cap comparison](/lustre/orion/ast207/proj-shared/dfielding/CGL/cgl_majeski_beta10_matched_20261016/passive_cap_comparison_20261017/RESULTS.md)
supports faster bulk comparisons but shows setting-dependent populations near
the limiter thresholds. Exact inequalities in history count those populations;
they do not measure the size of a threshold excursion.
The active beta-16 AMR continuations use bounded-face source `197bafc0f` and have
separate manifests/results. Preserve that distinction when interpreting the
comparison or selecting the next qualification run.
Their [matched resolution comparison](/lustre/orion/ast207/proj-shared/dfielding/CGL/cgl_beta16_uniform_amr_20261016/handoff10_to14/analysis_20261017/README.md)
waits for all four t=14 endpoints, compares [11,14], and retains [10,11] separately
as refinement adjustment. Further refinement remains disabled pending that review.

## Pressure-balance supplement

`scripts/analyze_cgl_pressure_balance.py` computes the signed cross spectrum
between perpendicular pressure and magnetic pressure, using two FFTs per field.
Give it a manifest with the same `tasks`/`cases` structure, new cache paths, and
only the matched times to average; it leaves the full turbulence caches intact.
Its plotting needs no history file. Use absolute cache paths in this supplement.

```bash
python scripts/analyze_cgl_pressure_balance.py --self-check
python scripts/analyze_cgl_pressure_balance.py pressure_manifest.json
python scripts/analyze_cgl_pressure_balance.py pressure_manifest.json --plot
```

The signed co-spectrum is `Re[p_perp_hat conjugate((B²/2)_hat)]`, with the same
physical shells and real-FFT weights as the auto-spectra. The band correlation
integrates cross power and normalizes by the square root of the two integrated
auto-powers. It is calculated from time-averaged spectra, not by averaging the
correlation coefficient. The total-pressure residual is
`E(p_perp+B²/2)/(E(p_perp)+E(B²/2))`: zero means exact cancellation, one means
uncorrelated fields, and two means equal in-phase fluctuations. Both phase and
amplitude therefore matter. The default reported band is `4 <= k Lx/(2 pi) <= 32`.

The summary separately records snapshot correlations, temporal scatter and
early/late windows; these are not independent-sample confidence intervals.
Self-checks cover exact anti-phase/in-phase/quadrature fields and a rectangular
box Fourier mode. Every field checks Parseval and total-pressure spectral
closure. This is a pressure-balance diagnostic, not a complete perpendicular
force balance: magnetic tension and inertia are not evaluated here.

Job `5698013` qualified all 58 Majeski fields in the common [6.25,14] window
and two full-resolution 768³ fields. The latter took 43/50 seconds and 30.85 GiB
peak RSS per process, using two 16-thread processes on one node. Their pressure
auto-spectra exactly match the existing full-analysis caches. The complete job,
including plots, took 1m38s on one node; no gradient FFTs were repeated.
