# Logic / Mathematical / Technical Audit

## Actionable findings

1. **High: show the second reallocation step; without it, the draft does not
   actually explain why `output29` loses angle.**

   The draft currently says that the inner `coord_abscostheta` reallocation
   erases earlier rows and then writes a valid angle row
   (`corruption_root_cause.tex:26-31`, `:61-63`). That explains the radius
   collapse in radius-first streams. By itself, it predicts that `output29`,
   whose order is `(coord_abscostheta, coord_r)`, should retain both angle and
   radius. The missing step is:

   - the bad inner allocation uses extent `nmb`, while the normal capacity
     allocation uses `nmb_max`;
   - PDF output calls `ComputeDerivedVariable()` for every later PDF dimension
     and variable weight, including names that are not derived
     (`src/outputs/pdf.cpp:92-103`);
   - on a later call, the entry capacity guard can see
     `derived_var.extent(0) < nmb_max` and reallocate again
     (`src/outputs/derived_variables.cpp:40-65`);
   - that second reallocation erases the freshly written angle row, after which
     later rows such as radius are recomputed.

   The historical bad line is visible at
   `ec009a27b87029db6f35a95045ae624ff5c0d1e4:src/outputs/derived_variables.cpp:1619-1621`.
   The archived input embedded in
   `/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/bin/node_00000000/gotham.hydro_w.00028.bin`
   sets `max_nmb_per_rank = 75`, so the `nmb` versus `nmb_max` distinction is
   not hypothetical.

   **Required revision:** replace the one-step cartoon with a two-step sequence:
   inner reallocation erases rows before angle; a subsequent entry-guard
   reallocation erases angle; later requested rows are then recomputed. State
   that the second erasure requires `nmb < nmb_max`.

2. **High: make the survival rule explicit rather than saying only that "a
   marginal over an uncorrupted axis can survive."**

   The precise rule is order-dependent:

   - rows computed before the bad `coord_abscostheta` call are erased by the
     inner reallocation;
   - `coord_abscostheta` is written, but can be erased by the next
     `ComputeDerivedVariable()` entry reallocation;
   - rows computed after that second reallocation survive to binning;
   - stored mass/volume weights are outside `derived_var`, while original
     variable flux weights are computed late enough to survive.

   This rule explains `output29`: its final radius row and `edot_sph` weight
   survive, all angle values collapse to zero, and summing over angle removes
   the wrong-bin assignment. The current formal section
   (`corruption_root_cause.tex:72-85`) accommodates the observation but does not
   derive this specific result.

   **Required revision:** state the row-lifetime rule and walk through
   `output29` call-by-call. Use the archived ordering recorded in
   `tools/gotham_pdf_rebuild/tests/reference.py:83-97`.

3. **Medium: qualify the statement that the bug does not change weights.**

   `corruption_root_cause.tex:77-84` correctly says that the bin coordinate can
   change without changing the weight, but this is not a general property of
   the bug. A variable weight may itself live in `derived_var` and may be
   erased or recomputed depending on call order. It is true for the fourteen
   GOTHAM products because their mass/volume weights are external to
   `derived_var`, and their `mdot_sph`/`edot_sph` weights are computed after the
   destructive calls. PDF output explicitly computes a variable weight after
   all axes (`src/outputs/pdf.cpp:97-103`).

   **Required revision:** say "for this fourteen-product catalog, the final
   weight rows survive" rather than implying a general guarantee.

4. **Medium: the proposed regression test is underspecified and can miss the
   `output29` failure mode.**

   A test that merely contains radius and `coord_abscostheta`
   (`corruption_root_cause.tex:198-202`) can pass if angle is last and no later
   call triggers the entry guard, or if `nmb == nmb_max`. It also does not
   distinguish the two destructive reallocations.

   **Required revision:** require two ordered products, `(coord_r,
   coord_abscostheta, trailing_dimension)` and `(coord_abscostheta, coord_r)`
   with a trailing variable weight, while forcing `nmb_thispack <
   nmb_maxperrank`. Compare both full histograms to an oracle, not merely
   "nontrivial population."

5. **Medium: distinguish the archive-wide claim from the eight-shard figure.**

   The "all five simulations" and "twelve of fourteen streams" claims are
   supported by the raw sparse-shard audit summarized in
   `/lustre/orion/ast207/proj-shared/gotham/analysis/PDF_OUTPUT_CORRUPTION_SUMMARY.md:42-65`.
   The quantitative figures in this document use only the first eight
   `res_8pc/phase2` shards, as recorded in
   `docs/explanatory_writeups/figure_metrics.json:48-79`. The draft currently
   moves between those scopes without naming the archive-wide audit artifact.

   **Required revision:** cite the archive-wide audit separately wherever the
   five-simulation/twelve-stream claim is made. Keep the eight-shard comparison
   labeled as the numerical repair-path control.

6. **Low: call `output31` a differential control, not a fully independent
   control.**

   `output31` avoids the implicated `coord_abscostheta` routine, so it is a
   strong control for this root-cause hypothesis. It still shares the PDF
   writer, radius, temperature, radial-velocity, and `edot_sph` machinery with
   affected products. "Independent control" at
   `corruption_root_cause.tex:169` is therefore stronger than the implementation
   supports.

   **Required revision:** use "independent of the implicated
   `coord_abscostheta` path" or "differential control."

## Claims that are technically sound

- The bad line was introduced by `ec009a27...` and removed by `db953863...`.
- `Kokkos::realloc` is destructive for the prior contents relied upon here.
- The histogram identity in `corruption_root_cause.tex:72-75` is appropriate
  when `w_p` is understood as a per-volume weight field.
- The eight-shard metrics are internally consistent:
  `output29` full `E_L1 = 1.917595...`, radial marginal
  `E_L1 = 1.6038e-11`, and total error `8.87e-13`
  (`docs/explanatory_writeups/figure_metrics.json:4-17`).
- Renormalization or relabeling cannot recover a lost cell-to-angle
  association.

## Recommended claim boundary

The root cause is demonstrated for the observed GOTHAM archive pattern, but the
document should demonstrate the complete two-reallocation mechanism rather
than relying on the reader to infer it. The eight-shard rebuild validates the
repair semantics and fingerprint; the archive-wide sparse audit establishes
scope.
