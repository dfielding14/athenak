# Frontier handoff: GOTHAM PDF reconstruction

## Current status

Source worktree:

```text
/ccs/home/dfielding/athenak-gotham-pdf-rebuild
branch: gotham-pdf-rebuild
base: origin/gotham-1.0 at db953863b
```

The 69-science/83-complete source was frozen on 2026-06-09 and passed all five
Frontier release gates under identity
`04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b`.
The subsequent science campaign completed all 161 inventory snapshots. A
read-only audit on 2026-07-10 checked every inventory mapping, manifest,
69-product ID list, shard count, output number and time, geometry/volume guard,
finite product sum, and nonempty header and payload; it found zero errors.

The reviewed source now differs from that frozen identity: it fixes the
physical-unit form of AthenaK's tiny cooling-time regularizer, makes the
release gate require an exact current source inventory, records the bundled
cooling table in that inventory, and corrects batch-log destinations. The
regularizer affects only the zero/tiny net-rate limit, but the source identity
has changed. No existing release artifact authorizes a future run of this
tree; freeze a new executable and repeat all goldens and canaries first.

The current 69-product science PDFs are under:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/production
```

The audit found 161 manifests: 29 `res_4pc`, 28 `res_4pc_highmdot`, 49
`res_8pc`, 25 `res_8pc_lowmdot`, and 30 `res_8pc_highmdot`. Together they
record 11,099,956,051,968 cells and 355,202,319,844,512 source payload bytes.
The manifest elapsed times and inventory node counts imply 125.55 reducer
node-hours. There is no live production plot tree for this campaign.

The earlier 76-product PDFs and their 26,565-PNG/1,155-movie descriptive plot
suite were moved to the `delete_after_rebin_review_20260606` trees. Their
immutable final campaign report remains a historical baseline and passed with
zero errors:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release/production-campaign.8fd47b0c5cc8b405779e37fc865685aa255088bb6aeaa00e389de953ee49388a.json
```

The historical PDFs and plots are respectively under
`pdf_rebuild/delete_after_rebin_review_20260606/production` and
`pdf_rebuild_plots/delete_after_rebin_review_20260606/production` within the
shared analysis-products tree.

## Latest pre-review Frontier release record

Frozen release identity:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release/frontier-release-identity.04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b.json
```

- Identity SHA-256: `04bbaa1e76bd25becff45d2b6e2975a64e477fdd71905def1eb70fd8ad260b9b`
- Source-bundle SHA-256: `82d118e7a31c0b29c7f5a4241015e0355244d1353f3a812e84fbde8bd536fa59`
- HIP executable SHA-256: `fc498bdf445c8f1481f807125fbdee9817ec1bab675959e49a02eb3e76e2c96b`
- Reducer source SHA-256: `c515c6babecddae984c71661c34f2b8d338c5d48c426f6adf1bc8706c4f1ce80`
- Validator SHA-256: `c96ee445404deb31dee8aa9744d3f71cf624deb108c375a0c1c220b7bb896e09`
- Performance-budget SHA-256: `9020a5b80adc9b5a48dec983e5bd877135bf9afc17d4a1ec40d716c6c925d1a1`

| Gate | Job | Shards | Reducer elapsed | Cells/s | Validator | Immutable pass report |
|---|---:|---:|---:|---:|---:|---|
| golden original | `4780097` | 1024 | 8.07 s | 10.71B | 31 pass, 0 fail | `golden-original.3f8db00b746fc6950dd4de4322280a2a7df919f67f4fb68379e657b1d5f234bc.json` |
| golden science | `4780098` | 1024 | 14.33 s | 6.03B | 27 pass, 0 fail | `golden-science.4c01864ef1247630495f912a1df552eba00c9f7fa7e5f0b0814dabac3afd3366.json` |
| golden all | `4780099` | 1024 | 15.34 s | 5.64B | 44 pass, 0 fail | `golden-all.576717ae481dab43d74a452fe0544494cb30507766f91eb2cb7e0fc950c4d413.json` |
| 256-node canary | `4780107` | 2048 | 14.87 s | 4.17B | 43 pass, 0 fail | `canary-256.898b5601037d3dceffb5b965161fe4fa0f4579a71d355945e86ff1db81401163.json` |
| no-control canary | `4780108` | 512 | 13.01 s | 0.99B | 28 pass, 0 fail | `canary-no-control.2de3350ba8d520bae44c86d1c33995816bb574fefc064af34bf05354c31d2693.json` |

All five runs validated geometry-to-logical-key consistency and passed every
predeclared performance-budget check. Pass reports and SHA-256 sidecars are in
`/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release`.
They are retained as evidence for the exact pre-review identity, not as a gate
for the changed source tree.

## Verified on Andes

- CPU standalone build succeeds.
- Reviewed serial source build plus the complete synthetic/reference, plotter,
  validator, and release-gate suite: `58 passed, 1 skipped` on 2026-07-10. The
  skip is the multi-rank launcher test because no `mpiexec`/`mpirun` is
  available on the login node.
- Reviewed Frontier HIP source build succeeds with CCE 20.0.0, ROCm 6.4.2,
  Cray MPICH 9.0.1, and bundled Kokkos 4.7.2 in the sanitized module shell
  documented in `README.md`.
- Real production layout: `64^3` cells per MeshBlock.
- K80 CUDA build succeeds after the documented one-line bundled-Kokkos Kepler
  patch.
- CPU/GPU one-block all-product agreement is near roundoff:
  representative normalized L1 errors `0` to `3.7e-15`.
- One complete real `res_8pc/phase2/00028` shard:
  - 320 MeshBlocks
  - 83,886,080 cells
  - 2.500 GiB payload
  - original repair set completed on one K80 in 99.69 seconds
- Eight complete real `res_8pc/phase2/00028` shards on four K80 GPUs:
  - 2,560 MeshBlocks
  - 671,088,640 cells
  - 20.000 GiB payload
  - original repair set completed in 381.63 seconds
  - Slurm job `3318742`, exit code `0:0`

Archived one-shard control comparisons at exact PDF sequence `00062`:

| Comparison | normalized L1 | normalized total error |
|---|---:|---:|
| output29 preserved radial marginal | `1.86e-11` | `1.35e-12` |
| output31 full 4D | `2.35e-7` | `1.36e-12` |
| output31 `(costheta,r)` marginal | `3.59e-12` | `1.36e-12` |
| output31 temperature marginal | `2.02e-7` | `1.36e-12` |
| output31 radial-velocity marginal | `1.63e-12` | `1.36e-12` |

The rebuilt output29 angle coordinate populates all 16 interior bins. Its full
comparison to archived output29 fails strongly, as expected, because the
archived angle axis is the corrupted quantity.

The four-GPU/eight-shard result independently confirmed the MPI reduction path:

| Comparison | normalized L1 | normalized total error |
|---|---:|---:|
| output29 preserved radial marginal | `1.60e-11` | `8.87e-13` |
| output31 full 4D | `3.87e-7` | `9.03e-13` |
| output31 `(costheta,r)` marginal | `7.23e-12` | `3.33e-12` |
| output31 temperature marginal | `2.07e-7` | `9.03e-13` |
| output31 radial-velocity marginal | `7.11e-8` | `9.03e-13` |

The eight-shard rebuilt output29 also populates all 16 interior angular bins.
Its intentionally failing full comparison to the corrupted archive has
normalized L1 `1.92`, while its total still agrees to `8.87e-13`.

## Correctness protections already implemented

- Exact GOTHAM schema and type validation.
- Exact complete embedded-header digest agreement across shards.
- Fixed-record divisibility and per-block index/geometry validation.
- Non-finite or non-positive hydro-state rejection.
- Exact expected contiguous `node_XXXXXXXX` shard manifest option.
- Global unique AMR leaf-key and ancestor-overlap validation.
- Full-snapshot domain-volume closure at relative tolerance `2e-9`.
- Exact AthenaK PDF boundary semantics and C-order strides.
- Explicit output-number and exact output-time mapping.
- Separate output tree, `.partial` payload publication, manifest written last.
- Independent Python oracle and full-snapshot real-control validator.

## Frontier golden sequence

This sequence completed successfully using the frozen identity recorded above.
The commands below are retained as the rerun procedure for a future identity.

Build the HIP executable using `README.md`, then run only the original repair
set for the golden complete snapshot using `jobs/frontier_golden_8pc.sbatch`:

```text
simulation: res_8pc
phase: phase2
cube sequence: 00028
PDF sequence: 00062
exact PDF time: 10.520002963096713
source shards/ranks: 1024
Frontier nodes: 128
```

Use one rank per GCD and `--expected-shards 1024`. Run
`validate_real_8pc.py` immediately afterward. Production is blocked unless it
passes all checks, including full AMR closure and archived global controls.

```bash
cd /ccs/home/dfielding/athenak-gotham-pdf-rebuild
FRONTIER_EXE="$PWD/build-gotham-pdf-frontier/gotham_pdf_rebuild"
FRONTIER_PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
RELEASE_IDENTITY=$(
  "$FRONTIER_PYTHON" tools/gotham_pdf_rebuild/jobs/release_gate.py freeze \
    --repo-root "$PWD" \
    --executable "$FRONTIER_EXE" \
    --artifact-dir /lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release
)
export RELEASE_IDENTITY

cd tools/gotham_pdf_rebuild/jobs
RUN_MODE=execute \
PRODUCTS=original \
CONFIRM_FRONTIER_GOLDEN=res_8pc/phase2/00028 \
RELEASE_IDENTITY="$RELEASE_IDENTITY" \
sbatch frontier_golden_8pc.sbatch
```

Only after the original golden pass report exists, submit the science and all
goldens separately. The production release gate requires pass reports for all
three product sets. Record each submitted job ID; the complete visual gate
must use the `PRODUCTS=all` output.

```bash
RUN_MODE=execute \
PRODUCTS=science \
CONFIRM_FRONTIER_GOLDEN=res_8pc/phase2/00028 \
RELEASE_IDENTITY="$RELEASE_IDENTITY" \
sbatch frontier_golden_8pc.sbatch

RUN_MODE=execute \
PRODUCTS=all \
CONFIRM_FRONTIER_GOLDEN=res_8pc/phase2/00028 \
RELEASE_IDENTITY="$RELEASE_IDENTITY" \
sbatch frontier_golden_8pc.sbatch
```

After the all-product golden validation passes, render it before releasing
production. The renderer reads only rebuilt dense PDF products; it does not
reread full cubes or modify reducer output. Run this step in a separate Andes
shell against the shared Orion output. The Python module names and venv below
are Andes-specific; no Frontier plotting environment is required.

```bash
cd /ccs/home/dfielding/athenak-gotham-pdf-rebuild
module reset
module load gcc/9.3.0 python/.3.11-anaconda3
source /lustre/orion/ast207/proj-shared/gotham/venv_gotham/bin/activate

ALL_JOB_ID=4780099
ALL_GOLDEN=/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/validation/frontier_golden/res_8pc/phase2/00028/all/"$ALL_JOB_ID"
PLOT_DIR=/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/validation/frontier_golden/res_8pc/phase2/00028/all/"$ALL_JOB_ID"

python "$PWD/tools/gotham_pdf_rebuild/plot_rebuilt_pdfs.py" \
  "$ALL_GOLDEN" \
  --output-dir "$PLOT_DIR" \
  --render-profile production \
  --workers 8 \
  --skip-existing
```

Inspect the generated `index.html`, especially the repaired `output29`
angular distribution, signed inflow/outflow products, and science products
with large reported underflow/overflow fractions. Read the per-product
metadata before interpreting amplitudes: panels are independently autoscaled
bin-integrated histograms, all axes are restricted to interior bins, and a
recorded visual-support mask omits at most `1e-4` of absolute panel weight.
Radial transport panels retaining radius are shell-normalized; other
transport panels remain raw integrated moments. No plot used for the release
gate may carry a `PARTIAL / UNVALIDATED REDUCTION` watermark.

## Frontier performance result

Predeclared numeric budgets live in `jobs/performance_budgets.json` and are
included in the frozen source bundle. All three golden product sets and both
canaries passed elapsed-time, node-hour, throughput, planned-buffer, node-count,
product-set, and timeout-margin checks. The largest planned peak buffer was
1.150 GiB/rank, below the 1.5 GiB/rank budget.

The current direct-global-atomic kernel is fast enough for production on the
measured golden and canary snapshots. Future optimization is optional and must
preserve the frozen release behavior:

```text
--products original
--products science
--products all
```

If future profiling identifies an atomic-bound workload:

1. Group products sharing axes and reuse transformed bin indices.
2. Use team/LDS-private accumulation for hot 2D products.
3. Keep direct global atomics for sparse high-dimensional products.
4. Double-buffer pinned host reads and device chunks to overlap `pread`, copy,
   and kernel execution.
5. Benchmark product-by-product GPU-aware Cray MPICH reductions against the
   current host-staged fallback.

Do not apply these optimizations to the repair set without rerunning the
synthetic oracle and full archived-control validator.

## Production inventory

The five `jobs/manifests/*.tsv` files contain 161 full cubes:

- 140 rows have matching archived PDF number and exact time.
- 21 uniform-phase rows have no archived PDF control and are marked `none`.
- Node count is `ceil(shards/8)`.
- Array concurrency defaults to one.

The production submitter and worker default to plan mode and require an
explicit confirmation token. Keep concurrency at one for the initial
production submission even though the full-snapshot and 256-node canary
throughput checks passed.

On 2026-06-05, the historical immutable-gated 76-product campaign completed 29
`res_4pc`, 28 `res_4pc_highmdot`, 49 `res_8pc`, 25 `res_8pc_lowmdot`, and 30
`res_8pc_highmdot` snapshots. It processed 11,099,956,051,968 cells and
355,202,319,844,512 source payload bytes. Reducer elapsed times ranged from
2.79 to 71.11 seconds.

The immutable release-gated plot campaign completed 35,742 production-profile
PNG plot groups and 161 snapshot `index.html` pages. The final campaign
validator checked every inventory mapping, rebuild manifest, product file,
plot metadata file, PNG count, source-view entry, and provenance checksum
before publishing the content-addressed pass report above.

The later user-approved scientific plot regeneration replaced those legacy
PNGs with 26,565 descriptively named PNGs and 1,155 H.264/YUV420p MP4 movies.
Its audit manifests are:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production/plot_campaign_manifest.json
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production/movies/movie_manifest.json
```

The June 8--9 radius-revised campaign repopulated the live PDF `production`
tree with the 69-product science catalog. The 2026-07-10 audit described at the
top of this handoff found all 161 expected snapshots and zero metadata or file
errors; summed reducer elapsed time was 2,762.91 seconds. No content-addressed
campaign report was found for this newer campaign, and the rebuild manifests
do not embed a release identity. Frontier worker logs show the
`04bbaa1e...` identity and five pass reports were verified before execution,
but a future campaign should also attest that identity directly in each output
manifest or a final campaign report.

One source-inventory correction was required for
`res_8pc_highmdot/phase1/00000`. The original directory contains 200 current
May 12 shards plus 50 stale May 7 shards. The audited source view includes
only contiguous nodes `00000000` through `00000199`; job `4769652` then
validated 25,192 unique AMR leaves and full domain-volume closure before
publishing the result. The source-view entry-metadata SHA-256 is
`35817ad7c167caf1a6a43b94a0aa62f2d7a82ff8530612c1df74192f9f92c400`.

Fail-closed orchestration attempts that did not publish a rebuild manifest are
retained in the Slurm logs. They include empty `MANIFEST_PATH` submissions,
the intentional all-250-shard duplicate-leaf rejection, and a symlink
precheck rejection. The immutable final report validates only the 161
successfully published snapshots.

The exact fail-closed operator environment for the historical 76-product
campaign was:

```bash
cd /ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs
RELEASE_DIR=/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/release
export RELEASE_IDENTITY="$RELEASE_DIR/frontier-release-identity.7a517368fe384bc3e793dca9bea96ff0042a5a6995815be1ceb25ab5abd9cf4f.json"
export GOLDEN_ORIGINAL_PASS="$RELEASE_DIR/golden-original.fdea5f6d147b25894dea642f0257a50710a7662d9c6f15317482f83978c95b30.json"
export GOLDEN_SCIENCE_PASS="$RELEASE_DIR/golden-science.f6cbc719c4ca7c98126363c15785c652819c0fb282de2cdcccb9b03be7a45d12.json"
export GOLDEN_ALL_PASS="$RELEASE_DIR/golden-all.7fe62b04a7d6d3110fa84818f82c17cd64639b0b153609ed5cf579e895b0de9b.json"
export CANARY_256_PASS="$RELEASE_DIR/canary-256.f8bb514afc911aca0b2c8a5c8ad066e6cd0c32539e280a92a7abe1deefb67eab.json"
export CANARY_NO_CONTROL_PASS="$RELEASE_DIR/canary-no-control.5bbb1603b2a5dba32100d2adf7381db96859c6beeb1613f7ced64a231caa67cf.json"

for simulation in res_4pc res_4pc_highmdot res_8pc res_8pc_lowmdot res_8pc_highmdot; do
  SIMULATION="$simulation" \
  SUBMIT_MODE=print \
  RUN_MODE=execute \
  PRODUCTS=all \
  ARRAY_CONCURRENCY=1 \
  CONFIRM_FRONTIER_PRODUCTION=GOTHAM_PDF_REBUILD_PRODUCTION \
  ./submit_frontier_production.sh
done
```

This exact block passed on 2026-06-05 with `SUBMIT_MODE=print`, and the
production campaign subsequently completed under the same frozen identity and
five pass reports. A future rerun must preserve the audited source correction
above or deliberately refreeze and repeat the release sequence.

## Known reconstruction limits

- Saved fluid fields are float32, so repaired PDFs cannot be bitwise-identical
  to online double-precision PDFs.
- Magnetic, particle, event-history, ionization, and exact shock-history
  products cannot be reconstructed from `hydro_w`.
- The added cooling products use the exact bundled CGM cooling tables and the
  source-term shielding/heating convention from `gotham-1.0`. They are signed
  net cooling-minus-heating diagnostics, matching AthenaK's saved
  `cooling_time` field rather than a positive-only radiative timescale. The
  cooling-weighted temperature products use cell-integrated net cooling
  luminosity in `erg s^-1`; the 3D cooling-rate product remains a local
  source-rate density in `erg s^-1 cm^-3`.
- The Kokkos Kepler source patch is Andes-only. Frontier HIP does not need it.
