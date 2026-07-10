# Logic / Mathematical / Technical Audit

## Actionable findings

1. **High: separate architecture validation from product-definition
   validation; the independent oracle does not cover all 62 science products.**

   The reducer defines 62 science products
   (`tools/gotham_pdf_rebuild/gotham_pdf_rebuild.cpp:645-817`), but the
   independent Python reference defines and checks only 20
   (`tools/gotham_pdf_rebuild/tests/reference.py:112-145`;
   `tools/gotham_pdf_rebuild/tests/test_reducer.py:168-173`). All 14 original
   products are independently checked
   (`tools/gotham_pdf_rebuild/tests/test_reducer.py:240-248`). The one-block
   all-product CPU/GPU runs establish execution and backend consistency for 76
   products, but both backends execute the same C++ definitions.

   `streaming_amr_reducer.tex:152-158` is individually accurate, but the status
   and bottom line can be read as a correctness claim for every science
   product.

   **Required revision:** state explicitly: "independent oracle: 14/14 original
   and 20/62 science products." The remaining 42 science products have
   execution/backend consistency and some internal-identity coverage, but no
   independent formula oracle. Either add independent cases for the remaining
   weight/axis families or label them provisional.

2. **High: the global AMR checks do not mathematically prove exact domain
   coverage.**

   The desired invariant in `streaming_amr_reducer.tex:71-76` is correct. The
   implementation checks nonnegative logical keys and positive extents
   (`gotham_pdf_rebuild.cpp:1038-1072`), duplicate keys,
   ancestor/descendant pairs, and scalar volume closure
   (`gotham_pdf_rebuild.cpp:1394-1419`). It does not check that physical
   extents match logical keys, that all blocks lie within the domain, or that
   distinct same-level geometries do not overlap. A compensating gap and
   out-of-domain or geometrically inconsistent block can satisfy the current
   scalar volume equation.

   **Required revision:** call this an operational closure check under the
   trusted AthenaK shard-format assumption, not proof that the accepted leaf
   set covers the domain exactly once. For a stronger guarantee, validate each
   block's geometry against its logical location and domain bounds.

3. **High: correct the memory model and narrow "bounded memory" to bounded
   field-data memory.**

   The formula at `streaming_amr_reducer.tex:78-85` counts one FP64 histogram
   but omits:

   - the full host histogram mirror created before MPI reduction
     (`gotham_pdf_rebuild.cpp:1557-1564`);
   - the local `BlockKey` vector retained for every processed block
     (`gotham_pdf_rebuild.cpp:1291-1295`, `:1742-1746`);
   - root-only global key storage plus the `unordered_set` used during AMR
     validation (`gotham_pdf_rebuild.cpp:1369-1396`);
   - root-only per-product receive storage (`gotham_pdf_rebuild.cpp:1575-1581`).

   The 0.292/0.325 GiB figures are the size of one device histogram, not peak
   per-rank memory. Snapshot size does not enter the chunked field-data term,
   but block metadata does scale with snapshot size and root validation gathers
   global metadata.

   **Required revision:** use a phase-wise peak estimate, approximately the
   maximum of processing memory, global-validation metadata, and
   reduction/publication memory. Say "bounded field-data memory" rather than
   unqualified "bounded memory."

4. **Medium: same-weight total closure is an internal checksum, not an
   independent correctness test.**

   Every cell is assigned to exactly one bin because underflow and overflow are
   retained, and products sharing a weight call the same
   `CellValues::weighted_value()` implementation
   (`gotham_pdf_rebuild.cpp:970-1017`, `:1270-1278`). Therefore same-weight
   totals can agree even if every axis is wrong or the shared physical weight
   formula is wrong. The check is useful for detecting dropped contributions,
   reduction errors, or publication errors. It does not validate bin semantics
   or the weight definition.

   The figure generator title "Independent products"
   (`docs/explanatory_writeups/scripts/generate_figures.py:316`) overstates this
   relationship. Its normalization also differs from the full validator:
   max absolute signed total in the figure generator
   (`generate_figures.py:299-305`) versus max sum of absolute bin weights in
   `validate_real_8pc.py:665-681`.

   **Required revision:** call the plot a conservative internal accounting
   checksum and do not present it as independent evidence.

5. **Medium: define `w_p` and the units of histogram values explicitly.**

   In `streaming_amr_reducer.tex:62-66`, the notation is correct only if `w_p`
   means a per-volume scalar such as density or flux density. The implementation
   function named `weighted_value()` already multiplies by cell volume
   (`gotham_pdf_rebuild.cpp:970-1015`), then atomically adds that returned value
   (`:1277-1278`). A literal mapping from the equation to that function would
   double-count volume.

   Also, these products are weighted histograms, not normalized probability
   densities. A radial marginal of `volume * edot_sph` is an integral of energy
   flux density over a shell volume; it is not directly a surface luminosity
   without the intended radial normalization.

   **Required revision:** define `w_p` as the pre-volume field and state the
   units/normalization convention.

6. **Medium: qualify "complete run" and the missing-shard guarantee.**

   Global mesh validation runs whenever both limiting options are zero
   (`gotham_pdf_rebuild.cpp:1745-1749`). Exact contiguous shard validation runs
   only when `--expected-shards` is nonzero
   (`gotham_pdf_rebuild.cpp:1099-1115`, `:1711-1716`). Thus the executable does
   not infer completeness merely because a run is unbounded.

   **Required revision:** state that a production-complete run must supply
   `--expected-shards`; otherwise "reject missing shards" is not guaranteed.

7. **Medium: publication is fail-closed only for consumers that require the
   manifest as a commit marker.**

   Each product header and payload is renamed separately, and products are
   published one at a time before the manifest is written
   (`gotham_pdf_rebuild.cpp:1451-1554`). A crash can leave a subset of
   final-named product files without a manifest. The individual renames are
   atomic; the output tree is not atomically published as a unit.

   **Required revision:** change `streaming_amr_reducer.tex:159-160` to say that
   the manifest is the required commit marker and incomplete trees without it
   must be rejected.

## Claims that are technically sound

- Direct accumulation over AMR leaf cells is the correct additive architecture
  for these weighted histograms.
- The 8 pc uniform-grid estimate is correct:
  `(400 / 0.008)^3 = 1.25e14` cells and eight float32 fields are about 4 PB
  decimal (`streaming_amr_reducer.tex:87-92`).
- C-order strides and underflow/overflow semantics match AthenaK
  (`gotham_pdf_rebuild.cpp:820-859`, `:887-897`).
- The measured Andes rates and 2.09x aggregate throughput ratio agree with
  `docs/explanatory_writeups/figure_metrics.json:159-179`.
- The draft correctly leaves complete-snapshot correctness and MI250X
  performance open.

## Recommended claim boundary

The native-AMR streaming architecture is strongly supported. The fourteen
repair products have a strong correctness case. The complete 62-product
science catalog, strict O(1) memory claim, and proof-level AMR coverage claim
are not yet supported at the same level.
