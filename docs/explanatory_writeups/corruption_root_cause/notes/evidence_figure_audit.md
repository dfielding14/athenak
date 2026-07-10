# Evidence / Figure Audit: Corruption Root Cause

Audit date: 2026-06-05

## Verdict

The eight-shard quantitative figures are real, reproducible from existing
artifacts, and interpreted correctly at the level of the main conclusion:
`output29` preserves its signed radial marginal while its angular distribution
is corrupted, and `output31` is a strong surviving control for this specific
bug.

The source-level mechanism is also strongly supported. The available AthenaK
history and bundled Kokkos implementation show exactly the destructive
reallocation sequence described in the draft.

The weakest important claim is the archive-wide scope. The available full
radius-audit job failed before completion and produced no CSV or Markdown
summary. Therefore, "audited latest snapshots across all five simulations;
twelve of fourteen affected" is not independently traceable to a completed
machine-readable audit artifact, even though the stream definitions and root
cause make it plausible.

## Audit Method

- Read the current draft, `figure_metrics.json`, and the figure generator.
- Independently recomputed the output29/output31 metrics from the eight
  archived sparse shards and rebuilt dense payloads.
- Inspected the source commits and bundled Kokkos `realloc` implementation.
- Inspected the archive-radius audit script, output directory, Slurm logs, and
  Slurm accounting.
- Did not regenerate figures and did not edit the `.tex` source.

## Figure Records

### `output29_corruption_signature`

- Figure paths:
  - `docs/explanatory_writeups/corruption_root_cause/figures/output29_corruption_signature.pdf`
  - `docs/explanatory_writeups/corruption_root_cause/figures/output29_corruption_signature.png`
- Status: newly generated explanatory figure from pre-existing data products.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
  - Reproducible command:
    `PYTHONPATH=/lustre/orion/ast207/proj-shared/gotham/analysis python docs/explanatory_writeups/scripts/generate_figures.py`
  - The exact original invocation is not preserved in a command log.
- Reader:
  - `/lustre/orion/ast207/proj-shared/gotham/analysis/gotham_analysis/readers/pdf.py`
- Archived inputs:
  - Header:
    `/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/pdf/pdf_coord_abscostheta_coord_r/node_00000000/gotham.header.pdf`
  - Payloads:
    `/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/pdf/pdf_coord_abscostheta_coord_r/node_00000000..00000007/gotham.00062.pdf`
- Rebuilt inputs:
  - Header:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/output29/gotham.header.pdf`
  - Payload:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/output29/gotham.00062.pdf`
  - Manifest:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/rebuild_manifest.json`
- Independently verified values:
  - Shape: `(18, 130)`
  - Archived and rebuilt embedded time:
    `10.520002963096713`
  - Full normalized L1 error:
    `1.9175953424932073`
  - Full normalized total error:
    `8.869061191580815e-13`
  - Signed radial-marginal normalized L1 error:
    `1.6037961007158223e-11`
  - Signed radial-marginal normalized total error:
    `8.868511479098382e-13`
  - Archived absolute angular support occupies only bin index `1`.
  - Rebuilt absolute angular support occupies all 16 interior bins,
    indices `1` through `16`.
- What it supports:
  - The rebuilt and archived eight-shard subsets carry nearly identical signed
    radial energy-flux profiles.
  - Their full joint distributions and angular support disagree strongly.
  - The total signed energy-flux weight is nearly unchanged despite the joint
    distribution failure.
- What it does not support:
  - It does not prove complete-snapshot AMR closure.
  - It does not prove all five simulations have the same corruption.
  - It does not independently isolate the hydro state, radial coordinate,
    weight formula, and reduction as separately correct; it shows their
    combined radial result is consistent.
- Caveats and interpretation flags:
  - The lower-right panel is an **absolute-weight angular-support marginal**,
    not the signed angular marginal. Calling it simply "the angular marginal"
    is imprecise.
  - The rebuilt manifest records eight processed shards but does not record the
    exact shard path list or original command line. Matching it to archived
    nodes `00000000..00000007` relies on the reducer's sorted-then-truncated
    shard selection, the output directory name, and the job context.
  - Zero `coord_abscostheta` maps to the first interior bin because its minimum
    is zero. It does **not** map to the underflow bin. Radius zero does map to
    radius underflow because the radius minimum is `0.1`.

### `archived_control_errors`

- Figure paths:
  - `docs/explanatory_writeups/corruption_root_cause/figures/archived_control_errors.pdf`
  - `docs/explanatory_writeups/corruption_root_cause/figures/archived_control_errors.png`
- Status: newly generated explanatory figure from pre-existing control data.
- Generator:
  - `docs/explanatory_writeups/scripts/generate_figures.py`
- Metrics manifest:
  - `docs/explanatory_writeups/figure_metrics.json`
- Output29 inputs: same archived and rebuilt inputs listed above.
- Output31 archived inputs:
  - Header:
    `/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/pdf/pdf_coord_costheta_coord_r_temperature_vel_sph_r/node_00000000/gotham.header.pdf`
  - Payloads:
    `/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/pdf/pdf_coord_costheta_coord_r_temperature_vel_sph_r/node_00000000..00000007/gotham.00062.pdf`
- Output31 rebuilt inputs:
  - Header:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/output31/gotham.header.pdf`
  - Payload:
    `/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/gpu_8pc_first8_original/output31/gotham.00062.pdf`
- Threshold source:
  - `tools/gotham_pdf_rebuild/validate_real_8pc.py`
- Independently verified normalized L1 values:
  - Output29 full, expected broken: `1.9175953424932073`
  - Output29 signed radial marginal: `1.6037961007158223e-11`
  - Output31 full 4D: `3.869265168184719e-07`
  - Output31 `(costheta, radius)` marginal: `7.2314517494673725e-12`
  - Output31 temperature marginal: `2.0686361120776285e-07`
  - Output31 radial-velocity marginal: `7.111634336486642e-08`
- What it supports:
  - The reducer strongly disagrees with the known-broken full output29.
  - The same reducer agrees closely with the preserved output29 radial
    marginal and multiple output31 controls.
  - This pattern is difficult to explain as a generic reader failure or a
    globally bad source state.
- What it does not support:
  - It does not validate a complete 1024-shard output31.
  - It does not validate products other than output29/output31.
  - It does not by itself prove the source-line change caused the archive bug.
- Caveats and interpretation flags:
  - The dashed `2e-3` line is only the full-validator limit for output31 full
    4D. The other green bars have different, tighter limits:
    output29 radial `2e-5`, output31 geometry `5e-5`, and output31 1D
    marginals `5e-4`.
  - "Complete rebuilt and archived 4D controls" is potentially misleading.
    The comparison is full-dimensional but only for the matched eight-shard
    subset.

## Important Numerical and Causal Claims

| Claim | Audit result | Evidence / caveat |
|---|---|---|
| Fourteen stream definitions, twelve radius-first | Verified for the latest visible phase directory of each of the five simulations | Shallow archive inspection found 14 streams and 12 `pdf_coord_r_coord_abscostheta*` streams in each simulation. This verifies layout, not corruption state. |
| Twelve streams affected across all five latest snapshots | Not fully traceable | `PDF_OUTPUT_CORRUPTION_SUMMARY.md` states this, but no completed CSV exists. The full audit job `3318069` failed after 81 of 238 streams because one partition manifest was incomplete. Jobs `3318065` and `3318066` were cancelled. |
| Full historical archive audited | Explicitly not supported | The draft now correctly says the stronger all-history claim remains untested. |
| Root-cause source history | Verified | Commit `ec009a27b87029db6f35a95045ae624ff5c0d1e4`, dated 2026-05-02, made the inner allocation unconditional. Commit `db953863b6c0e06d893b4a9738d4c111279e42fe`, dated 2026-06-02, removed it. |
| `Kokkos::realloc` destroys prior view contents | Verified | Bundled `kokkos/core/src/Kokkos_CopyViews.hpp` describes `realloc` as resizing while discarding old data; same-size reallocation initializes the view to zero. |
| Output29 fails differently because angle is first | Verified from code and data fingerprint | The inner reallocation shrinks extent 0 from `nmb_max` to `nmb`; the following radius calculation triggers the common capacity reallocation and rewrites radius after angle is lost. |
| Output29 radial marginal `E_L1=1.60e-11` | Verified | Independent recomputation from the raw eight-shard controls. |
| Output31 full 4D `E_L1=3.87e-7` | Verified | Independent recomputation from the raw eight-shard controls. |
| Control totals agree at roughly `1e-12` | Verified for the displayed eight-shard controls | Output29/output31 normalized total errors are approximately `8.87e-13` to `3.33e-12`, depending on marginal. |

## Required Corrections or Explicit Caveats

1. Do not present the all-five latest-snapshot corruption classification as
   backed by a completed archive-audit artifact. The available full audit
   failed and wrote no result files.
2. Call the lower-right output29 panel an absolute-weight angular-support
   marginal.
3. Distinguish radius zero, which lands in radius underflow, from
   `abs(costheta)=0`, which lands in the first interior angular bin.
4. Do not imply the `2e-3` dashed line is the threshold for every control bar.
5. Treat the eight-shard source identity as strongly inferred rather than
   perfectly recorded because the rebuild manifest omits exact input shard
   paths.

## Pass-Two Resolution

The final document presents the root cause as strongly supported rather than
demonstrated, includes the two-reallocation `output29` mechanism, calls
`output31` a differential negative control, and labels the all-five,
twelve-of-fourteen scope as reported rather than independently reproduced.
The same-state AthenaK A/B run and completed archive-wide audit artifact remain
open.
