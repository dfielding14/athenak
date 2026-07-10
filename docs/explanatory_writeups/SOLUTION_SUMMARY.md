# Solution Summary

This table summarizes the pre-Frontier argument captured by the explanatory
PDFs. The status column has been updated for the completed release sequence;
current artifact paths and source status are in
`tools/gotham_pdf_rebuild/FRONTIER_HANDOFF.md`.

| Solution / conclusion | What it explains | Best evidence | Biggest weakness | Best figure | Status |
|---|---|---|---|---|---|
| Shared-buffer reallocation is the archive root cause | Why the reported twelve streams collapse, `output29` is partly salvageable, and `output31` survives | The two-reallocation source mechanism predicts the observed mixed archive and the rebuilt controls reproduce that fingerprint | No preserved same-state bad-line/fixed-line AthenaK A/B run; full historical audit incomplete | `corruption_root_cause/figures/output29_corruption_signature.pdf` | Strongly supported causal diagnosis |
| Stream native AMR leaves directly into histograms | How to rebuild useful additive PDFs without a dense global cube | Separate oracle tests, bounded-memory real-shard GPU runs, and surviving original controls | Current checks do not prove arbitrary physical AMR coverage; science catalog lacks external validation | `streaming_amr_reducer/figures/andes_runtime_throughput.pdf` | Strongly supported architecture and original repair path |
| Gate production on one complete Frontier golden snapshot | How to avoid multiplying a common-mode error across 161 snapshots | The completed 1024-shard goldens, 2048-shard and no-control canaries, and 161-snapshot campaigns exercised the intended ladder | The newer 69-product campaign lacks a content-addressed final campaign report, and any reviewed-source rerun needs a new identity | `frontier_release_strategy/figures/validation_ladder.pdf` | Enforced and executed; repeat for each new source identity |

## Leading Answer

The root cause is selected strongly enough to act on. For the repair itself,
the native-AMR streaming reducer plus an enforced Frontier golden gate is the
leading path. These are complementary conclusions, not competing designs.

## Still Open

The same-state AthenaK A/B reproduction and a preserved full raw audit of the
historical archive scope remain open. The science catalog has exhaustive
internal oracle coverage but less independent external physical validation
than the original products. The reviewed source also needs a new frozen
identity and release sequence before any future production run.
