# Reproducibility / Traceability Audit

Audit date: 2026-06-05
Draft: `frontier_release_strategy.tex`

## Verdict

**FAIL on exact figure-artifact traceability and an underdefined performance
gate; PASS on inventory, completed-step evidence, pending-state honesty, and
compilation.**

The inventory, completed validation steps, exact golden-snapshot identity,
release-script guards, and still-pending Frontier state are traceable. The
document does not pretend that the golden run has completed.

Two release-traceability issues remain:

1. "Acceptable throughput" is not defined quantitatively, so the performance
   part of the release gate is not yet reproducible as a pass/fail decision.
2. The checked-in `frontier_release_strategy.pdf` is older than the current
   `.tex` source. The current source compiles, but the PDF beside it must be
   regenerated before release.
3. The current `validation_ladder` figure is not exactly reproducible from the
   current generator in the audited Python environment. Its underlying metrics
   reproduce exactly, but the current PNG pixels and PDF bytes differ.

## Figure And Input Traceability

| Item | Status | Trace |
|---|---|---|
| `figures/validation_ladder.pdf` and `.png` | FAIL exact artifact reproduction | Files exist and the underlying metrics reproduce, but the current PNG is not pixel-identical to output from the current generator. |
| `figures/production_inventory.pdf` and `.png` | PASS | Files exist. Built from the five static TSV manifests. |
| Inventory totals | PASS | Direct TSV count gives `161` rows, `140` mapped controls, and `21` rows without controls. |
| Golden identity | PASS | `frontier_golden_8pc.sbatch`, `validate_real_8pc.py`, and `jobs/manifests/res_8pc.tsv` agree on source `00028`, archived PDF `00062`, exact time `10.520002963096713`, cycle `202829`, 1024 shards/ranks, 128 nodes, and domain volume `400^3`. |
| Lightweight regeneration | PASS metrics / FAIL one figure artifact | Two staged regenerations reproduced inventory and validation-ladder data exactly. `production_inventory.png` is pixel-identical; `validation_ladder.png` is not. |
| Projected golden cell count | PASS as labeled estimate | `85,899,345,920` is computed as the measured eight-shard cell count times 128. The figure and metrics explicitly label it as an extrapolation, not a completed measurement. |

## Claim Traceability

Supported:

- The local synthetic/reference suite was rerun and reported
  `6 passed in 153.55s`.
- A fresh local CPU all-product one-block run matched the preserved K80
  all-product output for all 76 payloads:
  `global_max_abs=1.27898e-13`, `global_max_rel=1.0749e-14`.
- One real shard and eight real shards completed through the reducer.
- `sacct -j 3318742` reports the four-GPU/eight-shard Slurm job as
  `COMPLETED` with exit code `0:0`.
- The eight-shard log explicitly says global AMR closure was skipped for
  intentionally limited input.
- The static manifests contain 161 snapshots, 140 controls, 21 snapshots
  without controls, and 51 4 pc snapshots requiring 256 nodes.
- Job scripts default to plan mode, enforce confirmation tokens, refuse
  existing output by default, and require Frontier for submission.
- No Frontier HIP executable exists in the worktree, and no golden Frontier
  result or complete `validate_real_8pc.py` JSON report was found. The draft
  correctly labels these as pending.

Underdefined:

- **"Acceptable `original`, `science`, and `all` throughput on MI250X."** No
  numerical time, throughput, scaling, or node-hour threshold defines
  acceptable. Correctness gates are reproducible; the performance gate is not.
  Before production, record explicit thresholds and the exact benchmark
  command/output used to judge them.

No unsupported completed-Frontier claim was found. The missing Frontier
results are presented as missing, which is correct. Exact provenance for the
current `validation_ladder` artifact remains unresolved.

## Compilation

Status: **PASS with a layout warning; checked-in PDF stale.**

The current source compiled twice with `pdflatex -halt-on-error` to temporary
storage and produced a five-page PDF. The final log had no fatal errors and one
overfull box. At audit time:

- `frontier_release_strategy.tex` modification time:
  `2026-06-05 07:49:02 -0400`
- checked-in `frontier_release_strategy.pdf` modification time:
  `2026-06-05 07:39:19 -0400`

## No-Overwrite And Resource Check

- PASS: figure regeneration, synthetic tests, CPU comparison, and document
  compilation wrote only beneath `/tmp/gotham-repro-audit.T5cDY9`.
- PASS: both staged regenerations wrote only to temporary storage and reproduced
  the scientific metrics exactly.
- FAIL exact generated-artifact traceability: generated target figures changed
  concurrently during the audit, and the current `validation_ladder` artifact
  does not match output from the current generator.
- PASS: manifests, preserved outputs, logs, source files, and archived inputs
  were read only.
- PASS: no new Slurm job, Frontier job, external API, or image-generation
  service was used by this audit.
- PARTIAL: `RESOURCE_LEDGER.md` records the original lightweight work but not
  this audit's local test, one-block CPU comparison, or the concurrent
  target-figure changes. They are recorded here because this auditor was
  instructed not to edit other files.

## Commands Attempted

**PASS: staged figure regeneration**

```bash
mkdir -p /tmp/gotham-repro-audit.T5cDY9/docs/explanatory_writeups/scripts
cp docs/explanatory_writeups/scripts/generate_figures.py \
  /tmp/gotham-repro-audit.T5cDY9/docs/explanatory_writeups/scripts/
ln -s /ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools \
  /tmp/gotham-repro-audit.T5cDY9/tools
module load python/3.7-anaconda3
PYTHONPATH=/lustre/orion/ast207/proj-shared/gotham/analysis \
python /tmp/gotham-repro-audit.T5cDY9/docs/explanatory_writeups/scripts/generate_figures.py
```

Result: inventory and validation-ladder data exactly matched the checked-in
metrics. The current production-inventory PNG matched pixel-for-pixel; the
current validation-ladder figure did not match the current generator output.

**FAIL: exact current-figure comparison**

```bash
python - <<'PY'
# Load each target and staged PNG with matplotlib.image.imread and compare
# array shape and pixels exactly.
PY
```

Result: `production_inventory.png` matched exactly.
`validation_ladder.png` had the same shape but non-identical pixels.

**PASS: direct manifest inventory**

```bash
awk -F'\t' 'FNR>1 {n++; if ($9=="none") none++; else mapped++}
END {printf "rows=%d mapped=%d none=%d\n", n,mapped,none}' \
  tools/gotham_pdf_rebuild/jobs/manifests/*.tsv
```

Result: `rows=161 mapped=140 none=21`.

**PASS: Slurm accounting**

```bash
sacct -j 3318742 -X --noheader \
  -o JobIDRaw,JobName,State,ExitCode,Elapsed,NNodes,NTasks,AllocTRES
```

Result: job `3318742` is `COMPLETED`, exit code `0:0`, elapsed `00:06:27`.

**PASS: locate pending Frontier state**

```bash
test -x build-gotham-pdf-frontier/gotham_pdf_rebuild
find /lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild \
     /lustre/orion/ast207/proj-shared/brent/gotham/logs \
     -type f \( -iname '*golden*' -o -iname '*frontier*' \)
```

Result: Frontier executable absent; no golden Frontier result found. This is
consistent with the draft's pending status.

**PASS: current-source compilation**

```bash
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory=/tmp/gotham-repro-audit.T5cDY9/compile-latest/frontier_release_strategy \
  frontier_release_strategy.tex
```

Run twice from the document directory. Result: five-page PDF, no fatal error.

## Pass-Two Follow-up

After this audit, the quantitative figures were regenerated from the final
script and the checked-in six-page document was compiled twice from the revised
source.
The exact pending golden inventory is now derived from fixed-record preheaders
and file sizes. The gate remains intentionally incomplete: numeric performance
budgets, frozen artifact identity, immutable validator report, and enforced
production dependency are still required before execution.
