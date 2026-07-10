# GOTHAM PDF Reducer Jobs

These scripts operate the standalone `gotham_pdf_rebuild` reducer on Andes and
Frontier. They do not build the reducer.

All jobs default to `RUN_MODE=plan`, which invokes the reducer with `--dry-run`
and does not create output. Execution is always explicit:

```bash
RUN_MODE=execute sbatch andes_cpu_smoke.sbatch
RUN_MODE=execute sbatch andes_k80_correctness.sbatch
RUN_MODE=execute sbatch andes_8pc_limited.sbatch
```

The full Andes run has an additional snapshot-specific confirmation:

```bash
RUN_MODE=execute \
CONFIRM_FULL_8PC=res_8pc/phase2/00028 \
sbatch andes_8pc_full.sbatch
```

It defaults to the six currently schedulable K80 nodes. The task count follows
the actual allocation, so a different node count can be requested with
`sbatch --nodes=N` when Andes GPU node health changes.

The full Andes job defaults to `PRODUCTS=original` because Tesla K80
double-precision atomics are a correctness platform, not a production
performance platform. Override with `PRODUCTS=all` only intentionally.

## Frontier Production

The exact golden release gate is hard-coded in `frontier_golden_8pc.sbatch`.
It automatically runs the full validator, enforces the predeclared performance
budget, and publishes a content-addressed immutable pass report.

The latest pre-review release sequence passed on 2026-06-09:

```text
golden-original   4780097
golden-science    4780098
golden-all        4780099
canary-256        4780107
canary-no-control 4780108
```

The frozen release identity is:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release/frontier-release-identity.04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b.json
```

That identity and its five pass reports cover the exact 69-science/83-complete
source used for the June 8--9 campaign. Subsequent review changes mean a future
run must create a new identity and repeat the release sequence. Exact artifact
paths and campaign status are recorded in `../FRONTIER_HANDOFF.md`.

For a future identity, build the HIP executable, freeze it with
`release_gate.py freeze`, then submit the original repair set first:

```bash
RUN_MODE=execute \
CONFIRM_FRONTIER_GOLDEN=res_8pc/phase2/00028 \
RELEASE_IDENTITY=/path/to/frontier-release-identity.<sha256>.json \
sbatch frontier_golden_8pc.sbatch
```

Submit `PRODUCTS=science`, `PRODUCTS=all`, and both explicit canaries only
after the preceding immutable pass report exists. Production execution
requires the frozen identity plus all five exact pass reports; both the
submitter and every worker verify them before proceeding.

The five manifests are an explicit inventory of all 161 available full-volume
`hydro_w` snapshots. Shard counts were measured from the archive on June 4,
2026. Each row also carries the matching archived PDF sequence and exact
float64 PDF time when one exists; the 21 uniform-phase rows without archived
PDFs are explicitly marked `none`. Node counts use one MPI rank per Frontier
GPU and round up to eight ranks per node.

The submitter groups each simulation's snapshot indices by required node count.
It prints commands by default and refuses to submit unless run on Frontier:

```bash
SIMULATION=res_8pc ./submit_frontier_production.sh
```

Submit dry-run array jobs:

```bash
SIMULATION=res_8pc SUBMIT_MODE=submit RUN_MODE=plan \
./submit_frontier_production.sh
```

Submit production array jobs:

```bash
# First export RELEASE_IDENTITY and all five *_PASS variables exactly as
# recorded in ../FRONTIER_HANDOFF.md.
SIMULATION=res_8pc \
SUBMIT_MODE=submit \
RUN_MODE=execute \
PRODUCTS=all \
ARRAY_CONCURRENCY=1 \
CONFIRM_FRONTIER_PRODUCTION=GOTHAM_PDF_REBUILD_PRODUCTION \
./submit_frontier_production.sh
```

Repeat the production command for:

```text
res_4pc
res_4pc_highmdot
res_8pc
res_8pc_lowmdot
res_8pc_highmdot
```

`ARRAY_CONCURRENCY=1` is the conservative default. Increase it only after
confirming filesystem and queue behavior with the first production snapshots.

## Output Guards

The default output root is:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild
```

Every executing job atomically claims a previously nonexistent output
directory and refuses any existing target. Setting
`ALLOW_EXISTING_OUTPUT=YES` bypasses that guard and should only be used for an
intentional rerun.

After execution, jobs validate `rebuild_manifest.json`, exact processed shard
counts, positive block/cell/byte counts, finite product sums, and the existence
of every declared output file. The K80 correctness job also compares one-rank
and four-rank output payloads with a streaming numerical comparator.

The full Andes job validates complete reconstruction of one 8 pc snapshot. The
implemented partial real-data controls compare reconstructed `output31` and the
preserved radial marginal of `output29` against archived online PDFs. The
complete 1024-shard Frontier run, full validator, and both canaries have passed;
their immutable reports remain the required production release gate.

## Assumed Executables

```text
/ccs/home/dfielding/athenak-gotham-pdf-rebuild/build-gotham-pdf-cpu/gotham_pdf_rebuild
/ccs/home/dfielding/athenak-gotham-pdf-rebuild/build-gotham-pdf-andes-cuda/gotham_pdf_rebuild
/ccs/home/dfielding/athenak-gotham-pdf-rebuild/build-gotham-pdf-frontier/gotham_pdf_rebuild
```

Override `EXECUTABLE`, `ARCHIVE_ROOT`, `OUTPUT_ROOT`, or other explicit
variables at submission time when needed.
