# Resource Ledger

## Policy

No expensive new computation was needed for these explanatory documents.
Figures were generated from existing processed products, archived sparse PDF
shards, manifests, and completed-job logs.

## Resources Used

| Date | Resource | Purpose | Cost / scope | Output |
|---|---|---|---|---|
| 2026-06-05 | Andes login-node Python with `python/3.7-anaconda3` | Generate quantitative explanatory figures | Lightweight reads of eight archived sparse PDF shards per control, existing rebuilt dense PDFs, logs, and static manifests; no Slurm job | `figure_metrics.json` and figure PDFs/PNGs |
| 2026-06-05 | Local TeX Live `pdflatex` | Compile stand-alone explanatory PDFs | Local document compilation only | PDF documents beside each `.tex` source |
| 2026-06-05 | Subagents | Advocate, visual-design, audit, red-team, and editing passes | Text review only | Per-document `notes/*.md` |
| 2026-06-05 | Andes login-node CPU | Reproducibility audit | Lightweight synthetic/reference suite (`6 passed`), one-block 76-product CPU run in `/tmp`, and CPU/GPU payload comparison; no scheduler job | Results recorded in per-document reproducibility audits |
| 2026-06-05 | Andes login-node Python | Final figure regeneration after audit revisions | Same existing eight-shard controls, logs, and manifests; no full-cube reads | Final figure PDFs/PNGs and refreshed `figure_metrics.json` |
| 2026-06-05 | Local TeX Live and shell checks | Final pass-two build and traceability check | Each document compiled twice; section, review-note, figure, prompt, warning, hash, and `git diff --check` audits | Three final six-page PDFs and `FINAL_TRACEABILITY.md` |

## Explicitly Not Run

- No full-cube reads.
- No new reducer jobs.
- No Frontier jobs.
- No external APIs.
- No image-generation API calls. Detailed diagram prompts and TikZ schematics were used instead.
