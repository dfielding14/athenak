# Logic / Mathematical / Technical Audit

## Actionable findings

1. **Critical: the golden gate is a sound policy but is not currently
   machine-enforced before production.**

   The golden job validates only the basic rebuild manifest and then prints a
   reminder to run `validate_real_8pc.py`
   (`tools/gotham_pdf_rebuild/jobs/frontier_golden_8pc.sbatch:93-95`). It does
   not run the validator or write a durable golden-pass marker. The production
   submitter and worker require only a generic confirmation token
   (`submit_frontier_production.sh:46-50`;
   `frontier_production_array.sbatch:50-54`) and do not check a golden report.

   Therefore the prediction at `frontier_release_strategy.tex:102-103` that a
   failed golden run blocks every production array is procedural, not enforced.

   **Required revision:** say this explicitly. For an actual release gate,
   automatically run the validator, save its JSON, record executable/source
   identity, and require a matching pass artifact in the production submitter.

2. **High: "exact output number and exact float64 archived PDF time" is not
   fully checked by the validator.**

   The golden job hard-codes output `00062` and time
   `10.520002963096713` (`frontier_golden_8pc.sbatch:25-30`), so that script is
   correctly configured. However, `validate_real_8pc.py` does not assert that
   `manifest["output_number"] == "00062"` or that the manifest output time
   equals the archived exact time. It checks that rebuilt payload times equal
   the manifest exactly and that source time is within a six-significant-digit
   tolerance of archived time
   (`validate_real_8pc.py:534-585`). Its manifest checks cover source sequence,
   cycle, shard counts, and source path, but not output identity
   (`validate_real_8pc.py:337-373`).

   **Required revision:** do not claim the validator enforces item 4 at
   `frontier_release_strategy.tex:81`. Add explicit expected output-number and
   exact archived-output-time checks before calling the gate complete.

3. **High: "acceptable throughput" is undefined, so the performance gate cannot
   produce an objective release decision.**

   The release list requires acceptable `original`, `science`, and `all`
   throughput (`frontier_release_strategy.tex:87-88`), but no minimum
   cells/second, maximum wall time, maximum node-hours, peak-memory limit, or
   filesystem criterion is stated. The golden job's two-hour Slurm limit
   (`frontier_golden_8pc.sbatch:6`) is only a timeout, not a justified acceptance
   threshold.

   **Required revision:** define pass/fail budgets before the run. At minimum:
   maximum wall time and node-hours for each product set, required completion
   margin below timeout, and a rule for when optimization is mandatory.

4. **High: science/all release evidence is weaker than the draft implies.**

   `validate_real_8pc.py` requires original products `output29` and `output31`,
   so it cannot validate a `PRODUCTS=science` run by itself
   (`validate_real_8pc.py:872-880`). The independent oracle covers only 20 of
   the 62 science products (`tools/gotham_pdf_rebuild/tests/reference.py:112-145`),
   while all 62 are defined in
   `tools/gotham_pdf_rebuild/gotham_pdf_rebuild.cpp:645-817`. A science-only
   golden run is therefore primarily a performance benchmark.

   **Required revision:** require the `all` run to pass the full validator and
   compare the science-product payloads from `science` and `all` runs. Label the
   remaining 42 science definitions as not independently oracle-validated.

5. **Medium: "prove global AMR uniqueness and domain-volume closure" overstates
   what the implementation checks.**

   The complete run does activate checks that partial runs skip, which is a
   genuine logical distinction. But the implementation checks duplicate
   logical keys, ancestor/descendant pairs, and scalar volume closure
   (`gotham_pdf_rebuild.cpp:1394-1419`); it does not validate physical geometry
   against logical keys or domain bounds. The result is a strong operational
   gate under a trusted AthenaK-format assumption, not a mathematical proof of
   exact spatial coverage.

   **Required revision:** replace "prove" at
   `frontier_release_strategy.tex:118-119` with "exercise and pass the global
   operational checks," or strengthen the implementation.

6. **Medium: specify the 256-node canary and its acceptance criteria.**

   The draft correctly notes that one 128-node/1024-rank golden run does not
   exercise the 256-node/2048-rank production scale
   (`frontier_release_strategy.tex:160-166`). The current "first 256-node
   snapshot" is not an explicit reproducible gate. The manifests contain
   concrete candidates, including
   `res_4pc/phase2a/00003` at 2048 shards and 256 nodes
   (`tools/gotham_pdf_rebuild/jobs/manifests/res_4pc.tsv`).

   The existing `validate_real_8pc.py` is hard-coded to the 8 pc golden
   identity, so this canary also needs a generic or 4 pc-specific validator.

   **Required revision:** name the exact second canary, required product set,
   validator/control path, throughput threshold, and release condition for the
   remaining 256-node rows.

7. **Medium: label the domain-volume units.**

   The table gives expected domain volume as `400^3`
   (`frontier_release_strategy.tex:63`) and the validator uses the same bare
   value (`validate_real_8pc.py:18-24`). A volume requires units. If the saved
   coordinate unit is kpc, write `400^3 kpc^3`; otherwise label it explicitly as
   code-volume units.

8. **Low: distinguish validation thresholds from predictions.**

   The statement that full-golden errors should be "similar" to partial
   real-data errors (`frontier_release_strategy.tex:93-97`) is vague. The
   validator already defines concrete pass limits, including
   `2e-5` for the `output29` radial marginal and `2e-3` for `output31` full 4D
   (`validate_real_8pc.py:47-63`).

   **Required revision:** list those limits as release criteria. Treat similarity
   to the partial results as a useful diagnostic, not the pass/fail definition.

## Claims that are technically sound

- The remaining risks are predominantly common-mode, so multiplying the run
  count before testing the assembled Frontier path is unsound.
- The golden identity, source cycle, shard/rank count, and node count match the
  job and validator:
  `res_8pc/phase2/00028`, cycle `202829`, 1024 shards/ranks, 128 nodes.
- The absolute-weight denominators in `E_L1` and `E_total`
  (`frontier_release_strategy.tex:68-74`) are appropriate for signed flux
  histograms.
- The figure clearly labels the projected 1024-shard cell count as an
  extrapolation rather than a measured run
  (`docs/explanatory_writeups/figure_metrics.json:155-157`, `:181-201`).
- The inventory counts of 161 snapshots, 140 mapped controls, and 21 without
  controls agree with `figure_metrics.json:118-154`.
- One successful golden run should not be treated as universal proof; the draft
  correctly retains per-snapshot validation.

## Recommended claim boundary

Golden-first is the correct release strategy. At present it is a documented
human procedure, not a fail-closed software gate. Production should remain
blocked until the pass artifact, exact output-identity checks, quantitative
performance budgets, and explicit 256-node canary are defined and enforced.
