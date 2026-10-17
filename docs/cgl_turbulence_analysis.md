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
no field is synthesized. This pair used frozen pre-fix source `795c2e180c`.
The active beta-16 AMR continuations use bounded-face source `197bafc0f` and have
separate manifests/results. Preserve that distinction when interpreting the
comparison or selecting the next qualification run.
