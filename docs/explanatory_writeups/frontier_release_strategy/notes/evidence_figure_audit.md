# Evidence / Figure Audit: Frontier Release Strategy

Audit date: 2026-06-05

## Verdict

The static production inventory is well supported. The five manifests contain
161 rows, 140 mapped archived controls, and 21 no-control uniform-phase rows.
A read-only audit found all 140 representative control payloads and verified
that their embedded float64 times exactly match the manifest times.

The validation ladder correctly communicates that the complete Frontier run is
pending and logically different from the partial tests. The pass-one final
cell count was a low projection because the 1024 source shards do not all have
equal block counts. The final pass-two figure uses the exact fixed-record live
input inventory instead.

The golden-first strategy is supported. The golden run itself, Frontier HIP
build, 1024-rank reduction, complete AMR report, and MI250X profile remain
unverified.

## Audit Method

- Read both figures, the generator, metrics manifest, static production
  manifests, job scripts, validator, reducer source, logs, and Slurm accounting.
- Verified all 140 mapped representative control files and exact embedded
  times.
- Counted live golden source shards and fixed-size block records.
- Reran the current synthetic/reference suite:
  `6 passed in 77.59s`.
- Did not regenerate figures and did not edit the `.tex` source.

## Figure Records

### `validation_ladder`

- Figure paths:
  - `docs/explanatory_writeups/frontier_release_strategy/figures/validation_ladder.pdf`
  - `docs/explanatory_writeups/frontier_release_strategy/figures/validation_ladder.png`
- Status: newly generated explanatory figure from pre-existing logs/manifests
  plus the live fixed-record inventory of the pending input.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
- Inputs and source artifacts:
  - One-block GPU all-product artifact:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_one_block_all2/rebuild_manifest.json`
  - One-real-shard log:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/logs/gpu_8pc_shard.3318736.out`
  - Eight-real-shard log:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/logs/gpu_8pc_8shard.3318742.out`
  - Golden job configuration:
    `tools/gotham_pdf_rebuild/jobs/frontier_golden_8pc.sbatch`
  - Golden inventory row:
    `tools/gotham_pdf_rebuild/jobs/manifests/res_8pc.tsv`
  - Metrics record:
    `docs/explanatory_writeups/figure_metrics.json`
- Plotted values:
  - One block, all products: `262,144` cells, passed
  - One real 8 pc shard: `83,886,080` cells, passed
  - Eight real 8 pc shards: `671,088,640` cells, passed
  - Complete 1024-shard golden: exact live input inventory
    `86,467,149,824` cells, pending
- Independent source-record count for the pending golden input:
  - Live source files: `1024`
  - Fixed records: `329,846` MeshBlocks
  - MeshBlock size: `64^3` cells
  - Cells implied by live file sizes:
    `86,467,149,824`
  - The pass-one projection was low by `567,803,904` cells, or
    `0.6566700824%`, because 753 shards contain 320 blocks, one contains 326,
    and 270 contain 328.
- What it supports:
  - Real data, GPU execution, and a four-rank partial MPI reduction have
    completed.
  - The complete golden run is still pending.
  - The complete run is the first planned run that activates whole-input
    logical-key and summed-volume checks.
- What it does not support:
  - It does not show Frontier execution or performance.
  - It does not show a completed 1024-rank reduction.
  - It does not prove exact physical spatial coverage.
  - The exact plotted cell count is an input inventory, not a completed
    reducer result.
- Caveats and interpretation flags:
  - The final caption correctly labels the final value as a live fixed-record
    inventory rather than a measured completed run.
  - `figure_metrics.json` hardcodes the one-block passed status and `64^3`
    count without recording its source artifact path. The artifact exists, but
    the metrics file is incomplete provenance for that stage.
  - The partial eight-shard log explicitly says:
    `skipping global AMR closure for intentionally limited input`.

### `production_inventory`

- Figure paths:
  - `docs/explanatory_writeups/frontier_release_strategy/figures/production_inventory.pdf`
  - `docs/explanatory_writeups/frontier_release_strategy/figures/production_inventory.png`
- Status: newly generated explanatory figure from pre-existing static TSV
  manifests.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
- Input manifests:
  - `tools/gotham_pdf_rebuild/jobs/manifests/res_4pc.tsv`
  - `tools/gotham_pdf_rebuild/jobs/manifests/res_4pc_highmdot.tsv`
  - `tools/gotham_pdf_rebuild/jobs/manifests/res_8pc.tsv`
  - `tools/gotham_pdf_rebuild/jobs/manifests/res_8pc_highmdot.tsv`
  - `tools/gotham_pdf_rebuild/jobs/manifests/res_8pc_lowmdot.tsv`
- Mapping source:
  - `tools/gotham_pdf_rebuild/jobs/add_pdf_identity_to_manifests.py`
  - `/lustre/orion/ast207/proj-shared/gotham/analysis/catalog/generated/output_time_lookup.csv`
- Independently verified counts:
  - `res_4pc`: 29 snapshots, 28 mapped, 1 without control
  - `res_4pc_highmdot`: 28 snapshots, 28 mapped, 0 without control
  - `res_8pc`: 49 snapshots, 29 mapped, 20 without control
  - `res_8pc_highmdot`: 30 snapshots, 30 mapped, 0 without control
  - `res_8pc_lowmdot`: 25 snapshots, 25 mapped, 0 without control
  - Total: 161 snapshots, 140 mapped, 21 without control
  - All 21 no-control rows are `phase2_uniform`.
  - All 140 mapped representative control payloads currently exist and their
    embedded float64 times exactly equal the manifest `output_time`.
  - 51 manifest rows request 256 Frontier nodes.
- What it supports:
  - The static campaign inventory and per-simulation mapped/no-control counts.
  - The existence and exact embedded time of a representative archived
    output31 control for all 140 mapped rows.
  - The need for a no-control canary and internal invariants.
- What it does not support:
  - It does not prove every listed source cube currently has every expected
    shard.
  - It does not prove every archived control is scientifically intact.
  - It does not prove production jobs will run successfully at the listed node
    counts.
- Caveats and interpretation flags:
  - This is a static manifest inventory, not a live full-archive completeness
    scan.
  - "Matching archived PDF identity" means nearest same-phase control within
    the mapping script's time tolerance, with the exact selected payload time
    stored afterward. It is not a causal proof that the cube and control are
    otherwise equivalent.

## Important Numerical and Release Claims

| Claim | Audit result | Evidence / caveat |
|---|---|---|
| Golden identity `res_8pc/phase2/00028 -> PDF 00062 @ 10.520002963096713` | Verified | Static manifest row, archived payload time, and first-eight rebuild artifacts agree. |
| Expected source cycle `202829` | Verified from real-run logs/manifests | Both real-shard logs report cycle 202829. |
| Golden source shards `1024` and cells `86,467,149,824` | Verified live | Live file count is 1024; fixed-record preheaders and file sizes imply 329,846 MeshBlocks and the stated exact input cell count. |
| Frontier ranks `1024`, nodes `128` | Verified as planned configuration only | Hard-coded in `frontier_golden_8pc.sbatch`; not executed. |
| Domain volume `400^3` | Verified as source/configured expectation | Embedded source mesh spans `[-200,200]` in x/y/z; validator expects `400.0**3`. No complete closure report exists. |
| Global volume tolerance `2e-9` | Verified | Implemented in reducer and validator. |
| Six passing tests | Verified by audit-time rerun | `6 passed in 77.59s`; this is not an immutable release artifact. |
| One-block CPU/GPU near-roundoff agreement | Partially verified | Preserved CPU/GPU artifacts share 56 products and have maximum normalized L1 `4.345444850314735e-15`. The GPU artifact has 76 products, so 20 current products lack a paired preserved CPU/GPU comparison. |
| One-shard and eight-shard runs completed | Verified | Logs, manifests, and Slurm accounting agree. |
| Four-GPU/eight-shard exit `0:0` | Verified | Slurm job `3318742` completed with exit `0:0`. |
| Real-data controls pass by orders of magnitude | Verified for the eight-shard subset | The measured errors are thousands to millions of times below their corresponding validator limits. |
| Some 4 pc snapshots require 256 nodes | Verified | 51 static manifest rows request 256 nodes. |
| Frontier HIP has not been built/run | Supported by known artifact tree | No `build-gotham-pdf-frontier/gotham_pdf_rebuild` executable or Frontier result artifact was found. This is absence from the known workspace, not proof about every external location. |
| No complete AMR report / MI250X profile | Supported by known artifacts | No such report or profile was found. |
| Jobs default to plan mode and require confirmation | Verified in scripts | `common.sh`, golden job, production worker, and submitter implement these guards. |
| Current golden gate is advisory | Verified | The golden job performs a lightweight manifest check and tells the operator to run `validate_real_8pc.py`; production submission does not require a passing immutable validator report. |

## Required Corrections or Explicit Caveats

1. Resolved in pass two: the golden cell count was replaced with the live
   fixed-record count `86,467,149,824` and labeled as an input inventory.
2. Limit "one-block CPU/GPU agreement" to the 56 products with paired preserved
   artifacts unless a current 76-product paired run is recorded.
3. Treat the 161-row inventory as static planning evidence, not proof of live
   full-source completeness.
4. Continue to label the 1024 ranks and 128 nodes as planned, not measured.
5. Do not describe the golden strategy as enforced until production submission
   requires a frozen executable/source identity and passing validator report.

## Pass-Two Resolution

The final validation ladder uses the exact live input inventory: 1024 shards,
329,846 MeshBlocks, and 86,467,149,824 cells. The final document also calls the
gate advisory, requires a frozen artifact and immutable validator report,
narrows the AMR guarantee to implemented operational checks, and names the
256-node and no-control canaries. Frontier execution, objective performance
budgets, generic production-control validation, and enforcement remain open.
