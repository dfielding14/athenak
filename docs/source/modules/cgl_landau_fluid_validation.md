# CGL Landau-Fluid Validation

## Validation Tiers

Routine CPU testing covers a quantitative decay case, analytic uniform
collisional relaxation, oblique/parallel firehose threshold selection, strict
emergency-bound reporting without automatic correction, limiter occupancy,
explicit-versus-STS reference comparisons, and the live CGL FOFC regression.
It also covers the reduced active paper initializer, deterministic turbulence
seed selection, Alfvenic forcing orientation, nonrelativistic forcing-energy
work, once-per-cycle OU forcing advancement, forcing restart continuation, and
passive-Delta flow independence from diagnostic pressure anisotropy.
AMR, restart, and live FOFC checks also exercise the optional retained
RK-integrated CGL pressure-traction work path.
The broader scientific suite is a manual tier because it generates diagnostic
CSV files and figures and is intended for interpretation, not just pass/fail
gating.

## Illustrated Validation Note

The illustrated method note originally developed on `CGL-STS-LF` is now
maintained in this branch as `docs/cgl_lf_validation.tex`. It includes the
five method-validation figures plus representative reduced MKS24-oriented
analysis products, with production and reference-admission boundaries called
out explicitly.

[Download the illustrated CGL-LF validation and MKS24 analysis note
(PDF)](../_static/cgl_lf_validation.pdf)

The execution plan for expanding that note into a weak-guide-field turbulence
manuscript is maintained in
`docs/cgl_lf_weak_guide_manuscript_plan.md`. Manuscript prose and claim
discipline are governed by the recovered `docs/writing_style_guide.md`.

The manual suite covers:

- collisionless and finite-collision parallel and perpendicular damping;
- isolated field-strength-gradient transport and heat-flux limiting;
- mirror and firehose limiter activity;
- field-aligned and oblique linear waves;
- exact pure-CGL versus CGL-LF eigenmode comparisons.

## Canonical Workflows

Use the workflow driver from the repository root. It builds
`build-cgl-implementation/src/athena` in Release mode when that executable is
not already available.

For a first validated 1D STS damping run with strict diagnostics:

```bash
python3 scripts/cgl_lf_workflow.py quick
```

For the explicit-reference check against STS, including the finite-collision
split case:

```bash
python3 scripts/cgl_lf_workflow.py compare
```

For an archived local accuracy campaign with parallel/perpendicular
resolution and timestep sweeps, finite-collision coefficient limits, and
an extreme-cap STS/explicit exercise:

```bash
python3 scripts/cgl_lf_workflow.py accuracy
```

The `accuracy` bundle is operator evidence, not a reduced or
paper-resolution turbulence result. It records pgen reference errors for
`nx1 = 32, 64, 128`, three STS timestep ratios plus explicit references,
both local and background coefficient modes at three collision frequencies,
pre-cap face ratios above `10^4` in the extreme-cap cases, aligned thermal
decay with mean fields directed along `x`, `y`, `z`, and a periodic oblique
direction, and clean transport shutdown at `|B| <= bfloor`.

Finite-collision LF cases use two collision source updates per cycle, each
over the corresponding LF half-sweep duration. The routine
`cgl_lf_collision_relaxation` regression independently checks the analytic
law `Delta p(t) = Delta p(0) exp(-nu_coll t)`.

The explicit mode uses the same anisotropy-to-magnetic-moment split lifecycle
as STS. It is deliberately restricted to standalone reference checks.

For a supported two-dimensional AMR run with periodic boundaries, active LF
transport, strict monitoring, and conserved prolongation:

```bash
python3 scripts/cgl_lf_workflow.py amr
```

For reduced MKS24-oriented turbulence smoke cases, including active Alfvenic,
active random, and passive-Delta Alfvenic forcing:

```bash
python3 scripts/cgl_lf_workflow.py paper-smoke
```

This workflow archives the explicit parallel-firehose policy, forcing mode,
random seed, correlation time, injection parameter, mode bounds, state
tables, forcing tables, and the physical-shell/common-spectrum settings that
select the documented `abs(k) in (2*pi/L_parallel)[1,3]` and `k^-2` forcing.
Paper decks enable restartable cumulative RK-integrated applied source work.
They also enable `cgl_lf_record_pressure_work`, which retains total and
anisotropic applied hyperbolic CGL pressure-traction work in `lf_cpwrk` and
`lf_cawrk`.
The workflow checks clean LF safety diagnostics, positive applied work, active
global energy/work residual closure, zero parallel
forcing for the Alfvenic case, and nonzero parallel forcing for the random
case. Focused CPU regressions separately check a one-cycle OU transition,
multi-cycle RK2 companion source-work identity, and cumulative-work restart continuation. The
workflow also archives the explicit passive-Delta choice; a focused
regression checks that its flow fields are independent of diagnostic initial
anisotropy. It does not qualify paper-resolution or long-time active/passive
forcing statistics, or figure-level analysis.

For short reduced nonlinear product-path qualification, including
`8x8x16`, `16x16x32`, and `32x32x64` resolution variants, two heat-flux
wavenumber variants, and an oblique-firehose-policy variant:

```bash
python3 scripts/cgl_lf_workflow.py paper-convergence
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/paper-convergence/bundle
```

This mode overrides the standard input to `tlim = 0.02` and exists to check
nonlinear startup, snapshot output, and analysis/plot plumbing. It is too
short and too coarse for steady-state convergence or direct MKS24 figure
comparison claims.

Paper-scale workflow modes are defined for approved production campaigns,
not local validation:

```bash
python3 scripts/cgl_lf_workflow.py paper-standard \
  --authorize-paper-execution --output-dir /path/to/archive/run
python3 scripts/cgl_lf_workflow.py paper-nulim \
  --authorize-paper-execution --output-dir /path/to/archive/run
python3 scripts/cgl_lf_workflow.py paper-heat-flux \
  --authorize-paper-execution --output-dir /path/to/archive/run
python3 scripts/cgl_lf_workflow.py paper-compressive \
  --authorize-paper-execution --output-dir /path/to/archive/run
python3 scripts/cgl_lf_workflow.py paper-scale-separation \
  --authorize-paper-execution --output-dir /path/to/archive/run
```

The standard, limiter, heat-flux, and compressive inputs use `192x192x384`; the
scale-separation workflow adds `96x96x192` and `384x384x768` around its
reused standard-resolution case. All use `tlim = 10`, full-field binary
snapshots, checkpoints, and the `t = [8,10]` analysis window. Without
`--authorize-paper-execution` the workflow rejects these expensive modes
before creating a run bundle. The nine `paper-standard` definitions include
the eight active/passive, Alfvenic/random beta-10/beta-100 histories plotted
in MKS24 Figure 2(b), plus the active Alfvenic beta-1 case. The two
`paper-heat-flux` definitions add the nonnominal active beta-10 Figure 12
heat-flux cases; the nominal active and passive comparisons reuse standard
definitions. The two `paper-compressive` definitions add the active random
beta-1 and beta-100 sonic-correlation cases required by Figure 3; its other
four plotted cases reuse standard definitions. The two
`paper-scale-separation` definitions add the nonstandard
active Alfvenic beta-10 Figure 11 resolutions; the `n_perp = 192` comparison
reuses its standard definition. Defining these inputs does not constitute
paper-standard execution or a reproduction result.
In CGL primitive snapshots, legacy `eint` is `p_parallel` and the dedicated
`p_perp` output supplies the perpendicular pressure required for paper
anisotropy analysis.

Analyze a retained paper bundle without rerunning simulations:

```bash
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/existing/bundle
```

This action writes `analysis/diagnostics.json`, summarizes reduced histories,
and analyzes retained CGL binary snapshots when present. For standard inputs,
it reads each case's archived `analysis_t_start`/`analysis_t_end` values and
forms common-bin steady-window averages. Current products include
pressure/density PDFs, Figure 2(a)-coordinate joint pressure-density PDFs,
perpendicular shell spectra of total and compressive velocity, normalized
density fluctuations, separate parallel/perpendicular/magnetic pressures,
field-projected
`Delta p` and velocity-gradient spectra, the `b b : grad u` strain PDF,
pressure-stress transfer with a real-space closure check, selected-scale
alignment PDFs and their peak-versus-`k_perp` curve, and a cell-centered CGL
pressure-work split,
`integral[p_perp div(u) - Delta p (b b : grad u)] dV`. The anisotropic
stress term is compared with the direct pressure-transfer integral and is
labelled as active feedback or passive-only interpretation when that model
choice is archived. The pressure-transfer output also records the MKS24
dimensionless normalization
`T_Delta_p/T_total`, where the paper approximates
`T_total ~= E_K (2 pi u_rms/L_perp)`. For windows with at least two retained snapshots, the
ensemble record also reports cadence-limited trapezoidal time-integral
estimates of these pressure powers. When the run manifest contains the complete LF closure
choices, the analyzer also reconstructs a cell-centered heat-flux smoothing proxy,
`integral[-q_parallel b.grad(T_parallel) - q_perp b.grad(T_perp)] dV`,
including regularized/unlimited and cap-active contributions. This product
uses snapshot central differences; its optional snapshot-time quadrature and
the pressure-work quadrature are explicitly not the applied finite-volume
face flux or a closed energy budget.
When requested with `--eddy-samples <count>`, the analyzer also estimates
three-point second-order structure functions conditioned within 15 degrees
of the three-point local mean magnetic field, solves
`S2(ell_perp) = S2(ell_parallel)`, and reports
`eddy_anisotropy.velocity_perp` and `eddy_anisotropy.magnetic_perp`.
It bins `|ell|/L_perp` inside each angular cone, retaining that estimator
definition in the diagnostic JSON.
The deterministic sample count, logarithmic-bin count, and seed are retained
through `--eddy-samples`, `--eddy-bins`, and `--eddy-seed`.

The joint pressure-density output uses the plotted Figure 2(a) coordinates
`x = <p> delta rho/<rho>` and `y = delta p_parallel` or `delta p_perp`,
where `<p> = <(2 p_perp + p_parallel)/3>`, and renders one two-panel PDF per
analyzed case. The optional reference manifest can compare sampled surfaces
against these products by bilinear interpolation. The pinned raster is
decoded through its labeled colorbar into sampled surfaces. For this beta-10
panel, the pinned source fixes the pressure conversion: it states that `p0`
scales with `beta0` at fixed `B0` and reports `p0 = 100` for its beta-100
limiter runs, so the plotted run has donor `p0 = 10`; the matching AthenaK
deck has `p0 = 5`. The extractor admits the surfaces with pressure scale
`s = 0.5` and applies the two-dimensional PDF Jacobian to density and
uncertainty.

The compressive spectrum implements the Figure 3 quantity
`E_{khat dot u}(k_perp)` by Fourier-projecting velocity onto the full
three-dimensional wavevector before binning by `k_perp`. The thermal-pressure
products and AthenaK magnetic-pressure product `B^2/2` supply the fields
needed for Figure 6(a), and `spectra.density_fluctuation` supplies a
mean-normalized density diagnostic for Figure 4(b). These products do not
admit the plotted MKS24 spectral curves: their ordinate transforms remain
unqualified. MKS24 and its cited donor source state shell-binned definitions
and plotted spectral units, but do not state the discrete Fourier
normalization required to convert digitized absolute ordinates into these
AthenaK products. Figure 3's previously missing sonic-correlation and beta-1
random cases are now definition-only entries in the guarded
`paper-compressive` workflow; they have not been executed.

The analyzer runs synthetic binning/gradient/transfer/alignment and
heat-flux-sign/pressure-work checks and renders available history, PDF,
compressive/pressure-spectrum, joint-pressure-density, transfer, alignment,
heat-flux-proxy, and
pressure-work figures
under `figures/paper/`, including eddy anisotropy when requested. When an
externally prepared reference-data manifest
is supplied, it validates data/source checksums and reported pointwise
uncertainties, computes uncertainty-normalized residuals for supported
snapshot products, sampled joint-PDF surfaces, and threshold-volume history
series, and renders curve and surface comparison figures. Figure 2(b)
vector-curve extraction is
available below, together with Figure 5(b) normalized eddy anisotropy,
dimensionless Figure 7 lower-panel transfer,
Figure 8 selected-shell alignment PDFs, Figure 9, Figure 11 lower-panel, and
Figure 12 alignment, and Figure 13(b),(d) limiter curves. Dimensional Figure
7 upper spectra, Figure 11 upper spectra, Figure 12 lower spectra, and Figure
13(a),(c) curves remain excluded because their complete spectral or strain
ordinate transforms are not qualified by the pressure-coordinate conversion
used only for Figure 2(a). The remaining admission boundary requires
machine-readable author/archive data, donor diagnostic code, or an explicit
discrete Fourier-normalization statement; digitizing those absolute curves
alone is insufficient.
Figure 8 is admitted only through its reviewed raster contract: its embedded
RGB heatmaps are decoded against their labeled linear colorbar, selected
shell slices are checked against the paper's per-`k_perp` unit-normalization
statement and renormalized, and an absolute PDF-density uncertainty is
retained in the manifest. Figure 4(a)'s normalized-density PDF is also
admitted through its plotted `rho/<rho>` coordinate. Figure 5(b) is admitted
through its normalized length coordinates and the opt-in conditioned
structure functions. Figure 2(a)'s matching joint product and sampled-surface
comparison route are implemented; its labeled raster/color mapping and
source-derived beta-10 pressure conversion emit admitted surfaces while
retaining donor-coordinate samples for audit. Figure 3 and Figure 6(a)
now have matching emitted spectral fields, and Figure 4(b) has an explicit
normalized-density spectral product; their plotted spectral ordinates remain
unadmitted pending qualified transformations and authorized case execution.
Remaining dimensional or spectral-normalization MKS24
panels, exact
time-integrated/production local budget closure, and production comparisons
remain to be completed. For retained LF histories, the analyzer
also reports the RKL2-applied capped-face heat-flux contractions retained in
`lf_qprwrk` and `lf_qpewrk`. On refined meshes these contractions count
coarse/fine interfaces from their fine-side closure faces; they are not a
total-energy budget diagnostic, and their signed value need not equal the
positive cell-centered snapshot proxy. Operator-face heat-flux-cap fractions are separately
summarized from LF history counters over the same selected interval. When the
selected history contains `lf_cpwrk` and `lf_cawrk`, it additionally reports
the explicit-RK-applied CGL pressure-traction ledger. That ledger is evaluated
from the retained traction after the same AMR flux correction used by
momentum; passive-Delta decks record zero applied pressure work.
When the
default staged MKS24 manifest is available, `paper-analyze` also writes
`analysis/reference_provenance.json` with its archive and source-TeX
checksums; provide a different pinned copy with
`--reference-manifest <path>`.
Provide numerical or digitized panel curves separately with
`--reference-curves <path/to/curves.json>`.
For a standalone LF history file, the same interval calculation is available
through `scripts/analyze_cgl_lf_paper.py --lf-history <path>`.

Regenerate only the readable summary and diagnostic checks for an existing
paper bundle without plotting:

```bash
python3 scripts/cgl_lf_workflow.py paper-summary \
  --output-dir /path/to/existing/bundle
```

For the full scientific validation matrix and plots:

```bash
python3 scripts/cgl_lf_workflow.py full
```

The compatibility command `scripts/run_cgl_lf_validation.sh` invokes this
same `full` workflow. Use `--output-dir /path/to/results` with the Python
driver, or `OUTPUT_DIR=/path/to/results` with the compatibility command, to
give a bundle a persistent name.

## Result Bundles

Each executable workflow writes an ignored result bundle below
`build-cgl-implementation/cgl_lf_runs/<timestamp>-<workflow>/` by default:

```text
manifest.json
summary.md
logs/
history/
data/
figures/
```

`manifest.json` records the git revision, whether the worktree was dirty at
execution, its short status entries, the executable, source inputs, runtime
overrides, exact commands, result products, and measured checks. A dirty
bundle is suitable for implementation testing but is not evidence from an
immutable production source revision. `summary.md` is the readable pass/fail
digest. The `amr` summary reports whether refinement actually occurred,
normalized divB, invalid-state count, and normalized energy residual. The
`compare` summary reports maximum differences between STS and explicit final
states. The `paper-smoke` summary reports forcing orientation, mean beta,
`C_B2`, threshold-volume fractions, and total-energy change while retaining
both state and forcing tables.

The workflow also supports regenerating figures from a retained `full` bundle
with `plot`, and rebuilding `summary.md` and diagnostic values with
`summarize`, both using `--output-dir` to select that existing bundle.
Paper bundles use `paper-analyze` for JSON plus diagnostic plots and
`paper-summary` for summary-only regeneration.
Generated bundles are run products and are not committed by default.

## Paper Reference Staging

Paper-comparison analysis must use a pinned MKS24 reference copy rather than
an unversioned download. Stage the official arXiv `2405.02418v2` source
bundle into an ignored results area:

```bash
python3 scripts/stage_cgl_lf_mks24_reference.py
```

The utility downloads the source archive, safely extracts it, requires
`MKS24.tex`, and writes `manifest.json` containing SHA-256 checksums for the
archive and every extracted source or figure file. For an archive already
transferred to a compute system, use:

```bash
python3 scripts/stage_cgl_lf_mks24_reference.py \
  --archive /path/to/arXiv-2405.02418v2-source.tar.gz \
  --expected-sha256 <checksum-from-an-approved-manifest> \
  --output-dir /path/to/immutable/reference/arXiv-2405.02418v2
```

Do not commit redistributed paper source or figures as code artifacts unless
their redistribution permission is separately established. Archive the
staging manifest or immutable reference path with every quantitative paper
comparison bundle. The staged official `2405.02418v2` source inventory
contains TeX and rendered figure PDFs, not machine-readable curve tables.
For Figure 2(a), decode its labeled raster and apply its source-qualified
beta-10 pressure-coordinate conversion:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig2a.py \
  /path/to/arXiv-2405.02418v2/source/fig2a.pdf \
  /path/to/arXiv-2405.02418v2/source/MKS24.tex \
  /path/to/arXiv-2405.02418v2/digitized_fig2a_v2
```

The utility refuses either input if its SHA-256 differs from the pinned
`2405.02418v2` sources, reconstructs the tiled source-resolution chart, and
decodes its labeled linear colorbar. It retains two donor-coordinate CSV
surfaces with absolute joint-PDF density uncertainty `0.5`, masking only
black plotted overlays, and emits `surfaces.json` plus two
AthenaK-coordinate surfaces. The source-derived pressure scale is `s = 0.5`,
so `x_A = s x_paper`, `y_A = s y_paper`,
`z_A = z_paper/s^2`, and `z_uncertainty_A = 2.0`. Supply `surfaces.json` to
`paper-analyze --reference-curves` for the matching standard active beta-10
case.

For Figure 2(b), extract the plotted vector paths from the pinned source PDF
into an ignored or archival results area:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig2b.py \
  /path/to/arXiv-2405.02418v2/source/fig2b.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig2b_v1
```

The utility refuses a PDF whose SHA-256 is not the pinned Figure 2(b)
checksum, copies the source figure alongside its generated manifest, and
writes the eight active/passive unstable-volume histories as CSV curves. It
records an absolute `0.0025` plotted-line digitization uncertainty and omits
vertices hidden above the panel's `0.8` vertical limit.

For Figure 4(a), extract the eight normalized-density PDF paths from the
pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig4a.py \
  /path/to/arXiv-2405.02418v2/source/fig4a.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig4a_v1
```

This extractor maps the plotted `rho/<rho>` abscissa to
`pdf.density_fluctuation` by subtracting one, emits the active/passive
Alfvenic/random beta-10/beta-100 curves, and records five-percent relative
PDF ordinate uncertainty. It omits plotted-boundary-clipped vertices. The
companion Figure 4(b) `E_rho` spectrum remains excluded because the plotted
ordinate does not declare a `rho/<rho>-1` spectral normalization.

For Figure 5(b), extract the four normalized eddy-length curves from the
pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig5b.py \
  /path/to/arXiv-2405.02418v2/source/fig5b.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig5b_v1
```

This extractor maps blue/orange curves to
`eddy_anisotropy.velocity_perp`/`eddy_anisotropy.magnetic_perp`, and
solid/dash-dotted curves to the random/Alfvenic active beta-10 standard
cases. Both plotted axes are normalized lengths; the analysis normalization
uses `L_perp = sqrt(Lx Ly)`, which is unity in the shared `[1,1,2]` domain.
It records five-percent relative parallel-length uncertainty and omits
plot-boundary-clipped vertices. A comparison must explicitly enable the
costlier structure-function sampling, for example:

```bash
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/existing/standard/bundle \
  --eddy-samples 2000000 --eddy-bins 24 --eddy-seed 731 \
  --reference-curves /path/to/arXiv-2405.02418v2/digitized_fig5b_v1/curves.json
```

For Figure 7, extract the directly comparable lower-panel transfer-ratio
curves from the pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig7.py \
  /path/to/arXiv-2405.02418v2/source/fig7.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig7_v1
```

This extractor emits the active random, active Alfvenic, and passive
Alfvenic beta-10 curves as
`pressure_transfer.transfer_normalized_by_total`, with absolute `0.025`
digitization uncertainty. The plotted upper-panel gradient spectra remain
excluded pending the paper-to-AthenaK dimensional conversion.

For Figure 8, extract the active/passive beta-10 Alfvenic selected-shell
alignment PDFs from the pinned raster heatmaps:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig8.py \
  /path/to/arXiv-2405.02418v2/source/fig8.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig8_v1
```

The extractor decodes the embedded RGB panels through the figure's labeled
linear colorbar, samples shells `8,16,32,64,128`, and emits ten
`alignment.<shell>` curves. It requires each decoded slice to integrate
within `0.03` of unity before renormalizing it consistently with the paper's
stated per-`k_perp` normalization, and records absolute probability-density
uncertainty `0.05`. A comparison bundle must analyze the emitted shell set
with `paper-analyze --alignment-shells 8,16,32,64,128`.

For Figure 9, extract the eight standard-case peak-alignment paths from the
pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig9.py \
  /path/to/arXiv-2405.02418v2/source/fig9.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig9_v1
```

This extractor maps active/passive Alfvenic/random beta-10/beta-100 cases to
`alignment_peak.cos_theta`, uses absolute `0.005` ordinate uncertainty,
omits vertices clipped at the horizontal plot boundaries, and records its
recommended alignment shells for `paper-analyze --alignment-shells`.

For Figure 11, extract the three scale-separation lower-panel peak-alignment
paths from the pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig11.py \
  /path/to/arXiv-2405.02418v2/source/fig11.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig11_v1
```

This extractor maps active Alfvenic beta-10 cases at `n_perp = 96`, `192`,
and `384` to `alignment_peak.cos_theta`, using absolute `0.005` ordinate
uncertainty. It omits horizontal-boundary-clipped vertices and records the
alignment shells required by `paper-analyze`. The upper-panel kinetic spectra
are not emitted because their dimensionful scale is excluded pending the
paper-to-AthenaK conversion.

For Figure 12, extract the alignment-peak paths from the pinned source PDF:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig12.py \
  /path/to/arXiv-2405.02418v2/source/fig12.pdf \
  /path/to/arXiv-2405.02418v2/digitized_fig12_v2
```

This extractor maps four active/passive heat-flux cases to
`alignment_peak.cos_theta`, using absolute `0.005` top-panel uncertainty. It
omits vertices clipped at the horizontal plot boundaries and records its
recommended alignment shells. The lower-panel
`spectra.grad_parallel_delta_p` curves are not emitted because their
dimensionful scale cannot be compared until the paper-to-AthenaK code-unit
conversion is qualified. At the standard domain size, a Figure 12 comparison
can request the admitted shell centers directly; the workflow archives that
selection in its analysis-product metadata:

```bash
python3 scripts/cgl_lf_workflow.py paper-analyze \
  --output-dir /path/to/existing/bundle \
  --alignment-shells 8,12,16,17,20,24,27,28,32,36,37,40,44,47,48,52,\
56,57,60,64,67,68,72,76,77,80,84,87,88,92 \
  --reference-curves /path/to/arXiv-2405.02418v2/digitized_fig12_v2/curves.json
```

For Figure 13, extract the directly comparable dimensionless `beta Delta` PDF
and normalized pressure-transfer paths from the four checksum-pinned panel
PDFs:

```bash
python3 scripts/digitize_cgl_lf_mks24_fig13.py \
  /path/to/arXiv-2405.02418v2/source \
  /path/to/arXiv-2405.02418v2/digitized_fig13_v3
```

The extractor maps panel (b) to `pdf.beta_delta` for the three
`paper-nulim` cases plus the passive beta-100 comparison and panel (d) to
`pressure_transfer.transfer_normalized_by_total` for those same cases. It
records five-percent relative uncertainty for the panel-(b) PDFs, absolute
`0.025` uncertainty for panel-(d) transfer ratios, and omits
plotted-boundary-clipped vertices. Panels (a) and (c) are copied and
checksum-recorded for audit context but are not emitted as comparison curves:
panel (a) is dimensionful and panel (c) is plotted as
`bhat bhat : grad(u) / <B^2>`.

Other quantitative panels still require a recorded digitization procedure,
a qualified coordinate transform, or a separately provenance-tracked
numerical reference source. Figure 2(a) now has a matching joint
pressure-density product and admitted sampled surfaces, with raw
donor-coordinate raster samples retained for audit. Figure 4(b)'s
`E_rho` spectrum must not be
substituted for `spectra.density_fluctuation` without a matching ordinate
transformation. For the remaining absolute spectral and strain panels, the
manuscripts' shell-spectrum definitions and stated plotted units do not
provide the missing discrete Fourier normalization; admission therefore
requires machine-readable original data, donor diagnostic code, or an
explicit normalization statement in addition to any digitization.

Optional reference-data manifests passed through repeated
`paper-analyze --reference-curves` options
have `schema_version = 1`, a `provenance` object, and one or more `curves`
and/or `surfaces` entries. For `method = "digitized"`, provenance must name
either one source
figure and its SHA-256 digest or a nonempty `source_figures` list of such
pairs, together with the digitization tool and uncertainty procedure.
Each CSV curve must contain ordered `x`, `y`, and positive `y_uncertainty`
columns; its entry records `data_sha256`, target analysis `case`, supported
`product` (for example `spectra.grad_parallel_delta_p` or
`pressure_transfer.transfer_normalized_by_total` or
`alignment.<shell>` or `alignment_peak.cos_theta` or
`eddy_anisotropy.velocity_perp` or `history.unstable_fraction`), and optional
`interpolation` (`linear` or `loglog`). The analyzer fails closed on missing
uncertainty, checksum mismatches, or out-of-domain coordinates.
Each sampled surface CSV must contain `x`, `y`, `z`, and positive
`z_uncertainty` columns; its entry selects
`pressure_density_joint.parallel` or `.perpendicular` and supports
`interpolation = "bilinear"`. The renderer writes reference, analyzed, and
uncertainty-normalized residual views for each supplied surface.
Comparison product identifiers must be unique across all supplied manifests.
Dimensionful paper curves must additionally be withheld from such a manifest
until their paper-to-AthenaK normalization is explicitly established.

## Frontier Debug Qualification Preparation

Frontier qualification uses tracked operational tooling but remains a
separate, manually authorized execution tier. Validate the preparation and
accounting rules offline before copying the scripts to the project run root:

```bash
python3 scripts/frontier/cgl_lf_frontier.py self-test
```

For an actual Frontier debug-only qualification campaign, first build an
immutable HIP/MPI executable from a clean committed checkout using:

```bash
scripts/frontier/build_cgl_lf_frontier.sh
```

Then initialize the prescribed project area and prepare a reduced validation
job. The utility writes a retained manifest and concrete `debug`-QOS batch
script; it does not call `sbatch`.

```bash
python3 scripts/frontier/cgl_lf_frontier.py init
python3 scripts/frontier/cgl_lf_frontier.py prepare \
  --campaign qualification --run-name g002_paper_smoke_active \
  --purpose "strict reduced CGL-LF GPU smoke" \
  --acceptance-criterion "clean LF safety diagnostics and finite output" \
  --executable /lustre/orion/ast207/proj-shared/dfielding/CGL/build/<build>/src/athena \
  --input-file /lustre/orion/ast207/proj-shared/dfielding/CGL/repo/athenak-DF/inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput \
  --nodes 1 --walltime 00:20:00 --athena-walltime 00:15:00
python3 scripts/frontier/cgl_lf_frontier.py check-submit \
  --manifest /lustre/orion/ast207/proj-shared/dfielding/CGL/runs/qualification/g002_paper_smoke_active/manifest/prepared_run.json
```

Only after inspecting the printed `sbatch` command and confirming the
preflight should the job be submitted manually. Record the returned job ID
and completed accounting record sequentially:

```bash
python3 scripts/frontier/cgl_lf_frontier.py mark-submitted \
  --manifest <prepared_run.json> --job-id <jobid>
python3 scripts/frontier/cgl_lf_frontier.py record \
  --manifest <prepared_run.json> --job-id <jobid> \
  --result passed --notes "inspected reduced diagnostics"
```

For real runs the utility requires all source, executable, input, and output
locations beneath `/lustre/orion/ast207/proj-shared/dfielding/CGL`, limits
debug walltime to two hours, reserves against the 1000 node-hour testing
budget, checks that no other debug job is queued, and rejects
`paper-standard`, `paper-nulim`, `paper-heat-flux`, `paper-compressive`, and
`paper-scale-separation` inputs.
Paper-production simulations must not be run through this debug-only workflow.

## Diagnostics

With LF active, normal MHD history output appends cumulative counters:
`lf_nstage`, `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
`lf_mirror`, `lf_firehs`, `lf_hardbd`, `lf_qface`, `lf_qprcap`,
`lf_qpr10`, `lf_qpecap`, `lf_qpe10`, `lf_qprwrk`, `lf_qpewrk`, and
`lf_hwproj`. Strict
validation decks require zero floors, zero nonfinite/nonpositive states, and
zero hard-bound violations; limiter and cap counts may be nonzero when
intentionally exercised.
`lf_mirror` and `lf_firehs` are physical threshold-occupancy counters
evaluated with the selected policy. `lf_hardbd` is an emergency safety
counter and is evaluated whether or not corrective `backup_limiters` are
enabled. `lf_hwproj` is the cumulative number of primitive-refresh
applications of `limiter_hardwall = true`, including refreshed support/ghost
states; it records constraint activity rather than an admissibility failure
and is not a normalized active-cell occupancy. The `lf_q*` columns record evaluated LF operator faces and
parallel/perpendicular pre-cap flux ratios above `q_max` and `10*q_max`;
their interval fractions use `lf_qface` as denominator. Restart files retain
the baseline for every cumulative LF diagnostic column, including the face
and cap counters, so segmented-run intervals can be analyzed continuously.
`lf_qprwrk` and `lf_qpewrk` retain cumulative
RKL2-applied owned-face contractions of capped parallel and perpendicular
heat flux with discrete temperature jumps. At a coarse/fine interface, the
fine-side closure faces own the contraction because their fluxes are
restricted into the coarse update. These are signed operator contractions,
not a positivity or total-energy closure condition.

The firehose policy is explicit for threshold-sensitive work:
`cgl_firehose_threshold = parallel` selects the MKS24 production convention
`beta Delta <= -2`, while `oblique` selects `beta Delta <= -1.4` for
comparison studies and backward-compatible operator validation. MKS24
hard-wall decks additionally set `limiter_hardwall = true`, which projects
pressures to the selected threshold while preserving CGL internal energy and
replaces finite-rate limiter pressure relaxation.

The pre-existing `aam-D` label is retained for compatibility. In ordinary
output and restart state it denotes conserved CGL pressure anisotropy, not
the temporary magnetic-moment representation used internally during LF split
sweeps.

The `cgl_lf_paper` user history stores volume integrals named `mass`,
`kinetic`, `magnetic`, `therm_cgl`, `b2`, `b4`, `delta_p`, `abs_dp`,
`beta`, `mirror_vol`, `fire_vol`, `hard_vol`, `nu_eff`, `force_pwr`, and
`force_work`. `force_pwr` is an instantaneous proxy; `force_work` is the
cumulative net energy source actually applied by stage-weighted forcing,
including its zero-net-momentum projection, when the paper input enables
`record_injected_work`. It is updated as a companion variable with the same
explicit-RK recurrence used by the source-modified conserved state. Form means or fractions
using `volume`; for example, `C_B2 = b4*volume/b2^2 - 1`.
Threshold-volume columns use the selected firehose policy, whereas
`hard_vol` is safety-only. For active-Delta bundles, `paper-analyze` compares
the `force_work` interval difference to the conserved MHD `tot-E` interval
difference and reports the global residual for qualification; the reduced
workflow rejects a relative residual of `1.0e-8` or larger. Passive-Delta
bundles retain applied work but are not labeled as an active-CGL energy budget.

The AMR workflow retains both `.user.hst` and `.mhd.hst` products: user
history provides normalized divB, invalid-state, and anisotropy measures,
while MHD history provides total energy and LF counters. CGL LF rejects
`mesh_refinement/prolong_primitives=true`; conserved prolongation is the
supported AMR path.

## Developer Maintenance

Regenerate exact oblique-background eigenmode input decks only after
intentionally changing the linear reference convention with
`scripts/generate_cgl_lf_eigenmode_inputs.py`, then rerun the `full`
workflow. Routine pass/fail regressions remain under `tst/test_suite/cgl/`.

## Interpretation Limits

The oblique and eigenmode decks compare this ion CGL-LF closure against
linearized reference systems constructed for the same initial-value problem.
They are useful numerical checks, but do not by themselves implement electron
pressure physics or reproduce every historical comparison model.
