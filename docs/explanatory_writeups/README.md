# GOTHAM PDF Reconstruction: Explanatory Write-ups

These documents turn the PDF-corruption investigation and repair prototype
into three stand-alone arguments. They are meant to be read in order. The
first explains what failed, the second explains why the repair has its current
shape, and the third explains what still has to happen before the repair can
run across all five simulations.

Historical-scope note (2026-07-10): the PDFs and their review notes preserve
the evidence state when they were written, before the Frontier gates ran.
`FINAL_TRACEABILITY.md` records the subsequent 76-product release. The current
69-product campaign and reviewed source status are authoritative in
`tools/gotham_pdf_rebuild/FRONTIER_HANDOFF.md`; paths printed inside the older
documents may now refer to the `delete_after_rebin_review_20260606` archive.

## Recommended Reading Order

1. [Why the GOTHAM PDF Files Broke](corruption_root_cause/corruption_root_cause.pdf)

   This is the root-cause argument. It connects one destructive
   `Kokkos::realloc` to the unusual archive fingerprint: the reported twelve collapsed
   radius-first streams, a partially salvageable `output29`, and an intact
   `output31` control.

2. [Why the Repair Streams Native AMR Cells](streaming_amr_reducer/streaming_amr_reducer.pdf)

   This explains the preferred repair architecture. The desired PDFs are
   additive reductions over AMR leaves, so the reducer streams bounded chunks
   directly into GPU histograms instead of constructing a global dense cube.

3. [Why Frontier Production Starts with One Golden Snapshot](frontier_release_strategy/frontier_release_strategy.pdf)

   This explains the release gate that was subsequently implemented and
   completed. The complete Frontier sequence exercised whole-input AMR
   invariants, 1024- and 2048-rank behavior, archived and no-control
   validation, and MI250X performance before production.

The editable LaTeX source, quantitative figures, diagram prompts, and review
notes live beside each PDF. [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md) gives a
compact comparison. [RESOURCE_LEDGER.md](RESOURCE_LEDGER.md) records the
compute and external resources used for this explanatory phase.
[FINAL_TRACEABILITY.md](FINAL_TRACEABILITY.md) records the final build checks,
completed Frontier jobs, immutable artifact hashes, and unresolved claims.

## Current Judgment

The leading scientific conclusion is strongly supported: the archived
corruption comes from shared derived-coordinate storage being reallocated
during `coord_abscostheta`, including a second capacity realloc in the
`output29` order. A same-state bad-line/fixed-line AthenaK reproduction remains
the clean causal proof. The leading repair design is strongly supported:
stream native AMR leaves into MPI+Kokkos histograms with explicit invariants.

The remaining unresolved issues are scientific proof and validation breadth.
A same-state AthenaK A/B run is still missing, and the science catalog has
less independent external validation than the original products.

The original release-engineering requirements were completed. A frozen HIP
executable and source bundle passed the full 1024-shard golden sequence, a
256-node canary, and a no-control canary. The 76-product campaign then
calculated all 161 snapshots and rendered all 35,742 production-profile PNG
plot groups. The immutable all-pass campaign report and audited stale-shard
correction are recorded in [FINAL_TRACEABILITY.md](FINAL_TRACEABILITY.md).

That 35,742-PNG count describes the immutable release-gated legacy plot
campaign. On 2026-06-06, the visible plot tree was replaced by a
user-approved scientific suite containing 26,565 descriptively named PNGs and
1,155 QuickTime-compatible movies. The reducer outputs and immutable release
report did not change. A later 69-product radius-revised science campaign
repopulated the live PDF production tree for all 161 snapshots; its live plot
tree has not yet been regenerated.
