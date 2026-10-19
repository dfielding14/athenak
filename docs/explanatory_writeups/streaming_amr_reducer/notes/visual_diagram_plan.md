# Visual Diagram Plan

## Visual Language

Blue means trusted data, purple means GPU work, green means validated output,
red means rejected design, and gray means infrastructure. Schematics explain
the mechanism; they are not performance evidence.

## Opening Cartoon

Show AMR shards flowing through bounded `pread` chunks, GPU derivation and
binning, rank-local histograms, MPI reduction, and atomic publication. Cross
out the global dense cube below the main path.

Prompt: `diagram_prompts/streaming_dataflow.md`.

## Formal Schematic

Show coarse and fine leaf cells contributing with their own physical volumes.
Make the forbidden double-counting case obvious: an ancestor and its
descendants may not both enter the sum.

Prompt: `diagram_prompts/amr_volume_weighting.md`.

## Evidence Pair

Use `andes_runtime_throughput.pdf` to show the real path executes and
`same_weight_total_closure.pdf` to show internal accounting consistency.
