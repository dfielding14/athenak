# Final Traceability Record

Final pass: 2026-06-05

Path/status update: 2026-07-10. This is the immutable 76-product campaign
record. The later 69-product campaign and current reviewed source are described
in `tools/gotham_pdf_rebuild/FRONTIER_HANDOFF.md`.

This file records the final state of the explanatory write-ups after the
second writing pass, all required reviews, and the completed Frontier release
sequence. It distinguishes remaining scientific claims from the now-completed
production release-engineering gates.

## Final Documents

| Document | Pages | SHA-256 |
|---|---:|---|
| `corruption_root_cause/corruption_root_cause.pdf` | 6 | `2ba5a6f0024dc6b408ccad058c7a858ade42917594b13a850d53f3a6abc09ab7` |
| `streaming_amr_reducer/streaming_amr_reducer.pdf` | 6 | `29e6bbe1e61a7fbb998e1253cdbae7283b1af4b6b3950674cc28252f07baa5ed` |
| `frontier_release_strategy/frontier_release_strategy.pdf` | 6 | `637ac7c47c68c79be8af08c0f1f6f5bf80df263d4d565bec3385195b4a813b05` |

## Final Generation Checks

The quantitative figures and metrics were regenerated with:

```bash
module load python/3.7-anaconda3
PYTHONPATH=/lustre/orion/ast207/proj-shared/gotham/analysis \
python docs/explanatory_writeups/scripts/generate_figures.py
```

Result: success. The generator read existing archived sparse PDFs, rebuilt
test PDFs, logs, manifests, fixed-record preheaders, and file sizes. It did
not read full cube payloads or submit a scheduled job.

Each LaTeX document was then compiled twice with:

```bash
pdflatex -interaction=nonstopmode -halt-on-error DOCUMENT.tex
```

Result: three six-page PDFs. Final logs contain no overfull-box,
undefined-reference, fatal-error, or emergency-stop messages.

Additional final checks:

- `git diff --check`: pass.
- Required explanatory sections: 13 of 13 present in each document.
- Required role notes: 7 of 7 present in each document, 21 total.
- Quantitative figures: 6 PDF and 6 PNG versions.
- Diagram prompts: 6.
- No original simulation data was written or overwritten.

## Completed Frontier Release Sequence

The Frontier release sequence and complete production campaign completed on
2026-06-05.

Frozen artifacts:

- Release identity:
  `frontier-release-identity.7a517368fe384bc3e793dca9bea96ff0042a5a6995815be1ceb25ab5abd9cf4f.json`
- Source-bundle SHA-256:
  `0cbb0da68c2505280b846f22b6424f753e67fdae81b5d91684076df17425b885`
- HIP executable SHA-256:
  `ea53c3cf0f2ed8a949a6068f94da3d62f6bc230a23069e61fe4e5585ebefb692`
- Reducer source SHA-256:
  `ffce7ab54296a293382bf722e13e7b36aa530ee3110c02319cf09a6b88b87607`
- Validator SHA-256:
  `fda65cbffd3dfbbff4b378e5b385632bea4b11e03fc5676105d8dae9bf2e41d5`
- Performance-budget SHA-256:
  `1f24c3104164cb216e08532b9d5d23cbe869a797b4c4164f2a3fd179ce45cf4a`

All immutable artifacts and SHA-256 sidecars are under:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release
```

| Gate | Job | Nodes | Shards | Reducer elapsed | Cells/s | Validator result |
|---|---:|---:|---:|---:|---:|---|
| golden original | `4767325` | 128 | 1024 | 5.97 s | 14.48B | 31 pass, 0 fail |
| golden science | `4767326` | 128 | 1024 | 6.21 s | 13.93B | 22 pass, 0 fail |
| golden all | `4767327` | 128 | 1024 | 7.11 s | 12.17B | 42 pass, 0 fail |
| 256-node canary | `4767328` | 256 | 2048 | 4.39 s | 14.14B | 41 pass, 0 fail |
| no-control canary | `4767329` | 64 | 512 | 4.20 s | 3.07B | 26 pass, 0 fail |

Every run validated geometry-to-logical-key consistency and passed all
predeclared performance-budget checks. The all-product golden processed
86,467,149,824 cells in 329,846 MeshBlocks. The production campaign then used
the same frozen identity and five immutable pass reports to complete all 161
snapshots and their plots.

The immutable production campaign report is:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release/production-campaign.8fd47b0c5cc8b405779e37fc865685aa255088bb6aeaa00e389de953ee49388a.json
```

- Report SHA-256: `8fd47b0c5cc8b405779e37fc865685aa255088bb6aeaa00e389de953ee49388a`
- Calculations: 161 of 161, 76 products each
- Plots: 161 of 161, 222 groups each, 35,742 PNGs total
- Cells processed: 11,099,956,051,968
- Source payload bytes read: 355,202,319,844,512
- Reducer elapsed range: 2.79 to 71.11 seconds
- Validator result: zero errors, zero missing calculations, zero missing plots

Those historical production outputs were moved to
`/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/delete_after_rebin_review_20260606/production`.
Their plots were moved to
`/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production`.
The live PDF production tree now contains the later 69-product campaign, and
the live plot production path is absent.

The plot count above is the immutable release-gated legacy campaign. On
2026-06-06, the visible plot tree was replaced by 26,565 descriptively named
PNGs and 1,155 QuickTime-compatible MP4 movies. The current plot and movie
audit manifests are:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production/plot_campaign_manifest.json
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production/movies/movie_manifest.json
```

The final report also records the audited correction for
`res_8pc_highmdot/phase1/00000`: nodes `00000200` through `00000249` are stale
May 7 shards, while the contiguous May 12 nodes `00000000` through `00000199`
form the validated full snapshot. Corrected job `4769652` validated 25,192
unique AMR leaves and full domain-volume closure before publication.

## Quantitative Evidence Inventory

| Figure | Source class | Main use |
|---|---|---|
| `corruption_root_cause/figures/output29_corruption_signature.pdf` | Eight archived sparse output29 shards plus matching rebuilt test PDF | Shows the mixed-survival fingerprint: the full PDF is broken while the radial marginal survives |
| `corruption_root_cause/figures/archived_control_errors.pdf` | Archived output29/output31 controls plus rebuilt test PDFs | Separates the broken-order differential from the later-order control |
| `streaming_amr_reducer/figures/andes_runtime_throughput.pdf` | Completed Andes real-shard logs | Shows measured partial-input runtime and throughput; does not predict MI250X performance |
| `streaming_amr_reducer/figures/same_weight_total_closure.pdf` | Completed rebuilt-test manifest | Shows internal same-weight reduction closure; it is not an external correctness oracle |
| `frontier_release_strategy/figures/validation_ladder.pdf` | Completed pre-Frontier validation artifacts plus live fixed-record inventory | Shows the scale boundary that existed when the figure was generated |
| `frontier_release_strategy/figures/production_inventory.pdf` | Five production manifests | Shows 161 snapshots, including 21 without mapped archived controls |

Exact paths and metrics used to generate the figures are in
`figure_metrics.json`. The figures intentionally preserve the pre-Frontier
state at generation time. The previously pending golden inventory of 1024
shards, 329,846 MeshBlocks, and 86,467,149,824 cells is now also a completed
reducer result, recorded in the immutable golden pass reports.

## Review and Verification Record

Every document has:

- `solution_advocate.md`
- `logic_technical_audit.md`
- `evidence_figure_audit.md`
- `visual_diagram_plan.md`
- `red_team_review.md`
- `human_editor_pass.md`
- `reproducibility_audit.md`

The reproducibility review also recorded:

- Complete synthetic/reference, plotter, validator, and release-gate suite at
  the 2026-06-05 pass: 17 passed, 3 skipped.
- Fresh one-block CPU all-product run: 76 products.
- CPU/GPU comparison: all 76 payloads matched, with global maximum relative
  difference `1.0749e-14`.
- Completed Andes eight-shard job `3318742`: exit `0:0`.

These checks support the local mathematics, reducer implementation, and
partial real-input path. The complete Frontier sequence above independently
exercised whole-input invariants, 1024- and 2048-rank behavior, archived
controls, no-control validation, and MI250X performance.

The 2026-07-10 review rebuilt the final source with serial GNU/Cray MPI and
Frontier HIP CCE 20/ROCm 6.4.2 backends. The full current executable suite
reported 58 passed and one login-node launcher skip.

## Claims Still Open

1. A same-state bad-line/fixed-line AthenaK A/B run is still needed for the
   cleanest causal proof of the archive corruption.
2. The reported all-five-simulation, twelve-of-fourteen archive scope lacks a
   preserved completed raw audit product.
3. The current separate oracle covers every implemented product, but most
   science products still lack an external archived truth product.

## Resource Statement

No Slurm job, Frontier job, external API call, or image-generation call was
used while generating the explanatory documents and figures. The later
release-engineering phase used Frontier jobs `4767103` and `4767325` through
`4767329`, followed by the complete 161-snapshot production calculation and
plot campaigns. A later 69-product release used jobs `4780097` through
`4780099`, `4780107`, and `4780108` before its complete 161-snapshot science
campaign. Their current handoff is recorded separately as noted above.
