# PIC Engineering Readiness Fast Path

This path answers whether the current implementation is ready for exploratory
Bell and parallel-shock simulations. It is intentionally separate from the
publication qualification and claim-signoff program.

## Current Engineering Status (2026-06-26)

| Gate | Frontier job | Result |
|---|---:|---|
| ER1 corrected current | 4905833 | Pass: 26/26 cases in 44 s |
| ER2 compact linear Bell | 4905859 | Pass: growth and phase agree with theory |
| ER3 compact nonlinear Bell | 4905865 | Pass: peak `Bperp_rms/B0 = 2.0746` |
| ER4 shock transport/AMR/coupling | 4905864 | Pass: triad and coupled RK2 smoke |
| ER4 coupled shock transport | 4905888 | Pass: strong shock, injection, and transport |

The ER2 measurement used 89 snapshots. The measured normalized growth was
`0.90587175` versus `0.91651514` expected and the measured phase was
`-0.39613559` versus `-0.4` expected; both fits had `R^2 > 0.99999998`.
The ER3 coarse run retained 49 matched ten-product snapshots, reached a peak
`Bperp_rms/B0 = 2.07458316`, and ended at `1.14421976`.

ER4 proves shock evolution, particle injection/transport, dynamic AMR, and the
coupled deposition/feedback path. The compact coupled run retained 21 outputs
through `t=20` and completed with 129,069 particles. The measured shock front
was `x=199.5` versus the prescribed ideal position `x=200`, and the downstream
density approached the strong-shock value of four. Exact tag matching for the
initial 8,600-particle cohort found a maximum momentum increase of only
`0.077%`; this run does not resolve diffusive shock acceleration.

## Milestones

1. **ER0: build and mechanics**
   - Freeze one clean HIP/MPI Release executable.
   - Pass focused pusher, coupling, restart, and GPU/MPI tests.
2. **ER1: corrected Bell current**
   - Run the 26-case Q043 engineering matrix in one Frontier allocation.
   - Require every exact-byte deposited-current oracle and the cross-case
     current-spread check to pass.
3. **ER2: linear Bell**
   - Run 1D, 2D, and 3D baseline, fine-resolution, and MPI cases.
   - Require signed phase, analytical growth-rate agreement, finite
     conservation residuals, and serial/MPI consistency.
4. **ER3: nonlinear Bell**
   - Run deterministic-seed, particle-noise-seed, and enlarged/high-resolution
     pilots.
   - Require a resolved linear interval, departure from exponential growth,
     magnetic amplification, and bounded energy transfer.
5. **ER4: parallel shock**
   - Pass compact uniform, dynamic-AMR, and restart-continuation mechanics.
   - Run one paired coarse/AMR/fine shock triad and retain raw particle spectra.

ER2 and ER4 preparation may proceed in parallel after ER1 and the common ER0
mechanics gates pass. The full Q043 matrix, expanded seed sets, archive drills,
external comparisons, and named review remain later publication work.

## ER1 Commands

Validate the deterministic case selection locally:

```bash
python3 -B -m tst.publication.q043_engineering_readiness_v1 \
  --validate-selection
```

Preview the exact manifest:

```bash
python3 -B -m tst.publication.q043_engineering_readiness_v1 \
  --print-manifest
```

Submit the one-allocation Frontier run:

```bash
sbatch tst/publication/frontier_q043_engineering_readiness_v1.sh
```

The terminal evidence is
`engineering_readiness_summary.json` below the selected output root. A passing
record has `engineering_ready: true`; it does not authorize a scientific or
publication claim and does not replace the 132-case Q043 campaign.

## Bell And Shock Commands

```bash
sbatch tst/publication/frontier_q023_linear_bell_engineering_v1.sh
sbatch tst/publication/frontier_q019_nonlinear_bell_engineering_v1.sh
sbatch tst/publication/frontier_pic_shock_engineering_readiness_v1.sh
sbatch tst/publication/frontier_q011_coupled_shock_transport_engineering_v1.sh
```

The Bell jobs write `quicklook.json` under their job-specific
`engineering_readiness/` roots. The shock mechanics gate writes an artifact
manifest, CSV/JSON summaries, and comparison plot. The coupled transport job
runs a one-node uniform reduction of the Section 5.4 production deck and writes
time-resolved particle and exact-cohort statistics to `quicklook.json`.

## Publication-Style Figures

The checked-in figure generator produces two-column-width vector PDFs,
400-DPI PNGs, manuscript-ready captions, and a SHA-256 provenance manifest:

```bash
python3 -B -m tst.publication.make_bell_shock_publication_figures_v1 \
  --linear-root /lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/linear-bell-4905859 \
  --nonlinear-root /lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/nonlinear-bell-4905865 \
  --shock-root /lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/shock-acceleration-4905888 \
  --output-root /lustre/orion/ast207/proj-shared/dfielding/PIC/publication_figures/bell_shock_20260626
```

The figures are publication-quality presentations of compact engineering runs.
They support corrected linear Bell growth, coherent finite-amplitude Bell
response, shock formation, coupled injection, and particle transport. They do
not support turbulent Bell saturation or diffusive shock acceleration claims.
