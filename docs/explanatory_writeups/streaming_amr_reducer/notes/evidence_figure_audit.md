# Evidence / Figure Audit: Streaming Native-AMR Reducer

Audit date: 2026-06-05

## Verdict

The figures support the narrow claims they should support: the real-data
MPI+Kokkos path completed on Andes, and same-weight original products close
tightly within one eight-shard run. The source code, current test suite, and
existing artifacts strongly support the streaming architecture and original
repair catalog.

They do not establish Frontier performance, complete-snapshot correctness,
exact physical coverage, or external scientific validity of the 62-product
science catalog.

## Audit Method

- Read both figures, the generator, metrics manifest, reducer source, validator,
  tests, manifests, runtime logs, and Slurm accounting.
- Independently recomputed throughput, same-weight spreads, control metrics,
  product counts, and histogram memory.
- Reran the current test suite read-only:
  `GOTHAM_PDF_REBUILD_EXE=$PWD/build-gotham-pdf-cpu/gotham_pdf_rebuild python -m pytest -q tools/gotham_pdf_rebuild/tests`
  Result: `6 passed in 77.59s`.
- Did not regenerate figures and did not edit the `.tex` source.

## Figure Records

### `andes_runtime_throughput`

- Figure paths:
  - `docs/explanatory_writeups/streaming_amr_reducer/figures/andes_runtime_throughput.pdf`
  - `docs/explanatory_writeups/streaming_amr_reducer/figures/andes_runtime_throughput.png`
- Status: newly generated explanatory figure from pre-existing job logs.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
- Input logs:
  - `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/logs/gpu_8pc_shard.3318736.out`
  - `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/logs/gpu_8pc_8shard.3318742.out`
- Slurm accounting checked with:
  - `sacct -j 3318736,3318742`
- Independently verified values:
  - One-rank / one-K80 / one-shard run:
    - 320 blocks
    - 83,886,080 cells
    - 2.50002622604 GiB
    - 99.685631315 s reducer elapsed
    - 0.0250791031 GiB/s
    - 0.841506232 million cells/s
  - Four-rank / four-K80 / eight-shard run:
    - 2,560 blocks
    - 671,088,640 cells
    - 20.0002098083 GiB
    - 381.626483859 s reducer elapsed
    - 0.0524078141 GiB/s
    - 1.758495986 million cells/s
  - Aggregate rate ratio: `2.08970049`
  - Slurm job `3318742`: `COMPLETED`, exit code `0:0`
- What it supports:
  - The fixed-record reader, Kokkos derive/bin kernel, FP64 accumulation, MPI
    reduction, dense publication, and manifest path completed on real data.
  - Four K80 ranks processed and reduced eight real shards.
  - Aggregate throughput improved by about `2.09x` from the one-rank case.
- What it does not support:
  - It does not measure Frontier or MI250X performance.
  - It does not validate a complete snapshot.
  - It does not demonstrate strong scaling or efficient weak scaling.
  - It does not measure peak memory.
  - It does not measure the 62-product science or 76-product all-catalog path;
    both plotted real-shard runs used the 14-product `original` catalog.
- Caveats and interpretation flags:
  - "End-to-end" is the reducer's measured main path from just before shard
    processing through reduction and publication. It excludes process launch,
    Kokkos/MPI initialization, initial shard discovery, and product setup.
  - The one-shard Slurm allocation exposed four K80s, but the reducer step used
    one rank and one GPU; Slurm accounting confirms that distinction.
  - The figure is evidence of execution, not production efficiency.

### `same_weight_total_closure`

- Figure paths:
  - `docs/explanatory_writeups/streaming_amr_reducer/figures/same_weight_total_closure.pdf`
  - `docs/explanatory_writeups/streaming_amr_reducer/figures/same_weight_total_closure.png`
- Status: newly generated explanatory figure from a pre-existing rebuild
  manifest.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
- Input manifest:
  - `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/rebuild_manifest.json`
- Threshold source:
  - `tools/gotham_pdf_rebuild/validate_real_8pc.py`
- Independently verified plotted values:
  - `edot_sph`: `1.0532546461132228e-13`
    from output29/output31
  - `mass`: `6.415880453357775e-13`
    from output23/25/26/27/32/33/34/35/36
  - `mdot_sph`: exactly `0.0`
    from output28/output30
- What it supports:
  - Within this original-product eight-shard run, changing histogram axes does
    not materially change the total deposited value for products sharing the
    same computed weight.
  - It is a useful deposit/reduction bookkeeping invariant.
- What it does not support:
  - It does not externally validate any weight formula.
  - It does not validate the science catalog.
  - It does not validate AMR completeness.
  - These are separate products from the same kernel and same run, not
    statistically or implementation-independent measurements.
- Caveats and interpretation flags:
  - The figure generator normalizes spread by `max(abs(product total))`.
    `validate_real_8pc.py` normalizes its same-weight check by the maximum
    **absolute bin sum**. The dashed `2e-9` validator line is therefore not
    exactly the same metric as the bars, although the plotted metric is more
    conservative for signed products.
  - The zero `mdot_sph` spread is invisible on a logarithmic axis. The final
    pass-two figure labels it explicitly.
  - The pass-one figure title's word "Independent" was too strong. The final
    title describes this as a shared-weight closure check.

## Important Numerical and Implementation Claims

| Claim | Audit result | Evidence / caveat |
|---|---|---|
| Five simulations and 161 snapshots | Verified as static manifest inventory | Five TSV manifests contain 161 rows. This does not itself prove every live source snapshot is complete. |
| Original catalog uses 0.292 GiB/rank | Verified | Current executable reports 39,176,280 FP64 bins and `0.291886` GiB. |
| All 76 products use 0.325 GiB/rank | Verified | Current executable reports 43,593,352 FP64 bins and `0.324796` GiB. |
| Science catalog contains 62 products | Verified | Current executable `--products science --list-products` reports 62 science products. |
| Dense 8 pc cube estimate is `1.25e14` cells and 4 PB for eight float32 fields | Verified analytical estimate | Assumes a uniform 400 kpc cube with 0.008 kpc cells and excludes coordinates, masks, metadata, and working storage. |
| Current suite records six passing tests | Verified by audit-time rerun | `6 passed in 77.59s`. No prior immutable test report was found. |
| Independent oracle coverage | Verified from the reference catalog and tests | The suite checks all 14 original products and 20 of 62 science products against the separate Python reference. The remaining 42 science products do not have independent formula-oracle coverage. |
| Oracle covers "real 64^3 block layouts" | Misleading wording | The oracle's 64-cubed test is synthetic and validates one science product. Separate real-shard runs exercise actual 64-cubed blocks, but not against the Python oracle cell by cell. |
| Synthetic one-rank/two-rank invariance | Verified for science products | `test_science_products_are_mpi_invariant` checks the science catalog, not the original catalog. |
| Production reader accepts reconstructed files | Partially verified | The current test suite reads representative output29; the figure workflow also successfully reads rebuilt output29 and output31. This is not an all-product reader audit. |
| One block with all 76 products runs on GPU | Verified | `gpu_one_block_all2/rebuild_manifest.json` records 76 products and 262,144 processed cells; Slurm job `3318735` completed. |
| CPU/GPU one-block agreement is near roundoff | Verified only for the 56 products common to preserved CPU/GPU artifacts | Comparing `one_block_all` to `gpu_one_block_all2` gives maximum normalized L1 `4.345444850314735e-15`. The GPU artifact has 76 products; the preserved CPU artifact has 56, so 20 current products lack a paired CPU/GPU artifact. |
| Publication is fail-closed | Correct only under a manifest-required consumer contract | Product headers/payloads are renamed before the final manifest. A consumer that ignores the manifest can see a partial product set. |
| Global AMR checks prove exact physical coverage | Correctly not claimed in the revised draft | Code checks contiguous shard names, header consistency, unique logical keys, ancestor/descendant pairs, and summed volume. It lacks geometry-to-key validation. |

## Required Corrections or Explicit Caveats

1. State clearly that both real-shard throughput bars are for the 14-product
   `original` catalog, not science/all.
2. Do not call the same-weight products independent evidence.
3. Do not imply the dashed same-weight threshold uses exactly the same
   normalization as the plotted bars.
4. Split the oracle/64-cubed claim: synthetic oracle coverage and real
   64-cubed execution are different evidence.
5. Keep CPU/GPU "all-product" agreement limited to the 56 products with paired
   preserved artifacts unless a current 76-product CPU/GPU comparison is
   generated.

## Pass-Two Resolution

The final document makes the 14/14 original and 20/62 science oracle boundary
explicit, labels same-weight closure as an internal deposit/reduction
invariant, and separates synthetic oracle coverage from real-shard execution.
The final closure figure labels the exact-zero `mdot_sph` spread and no longer
uses "Independent" in its title. Complete-snapshot validation, exact physical
coverage, and external validation of the remaining science definitions remain
open.
