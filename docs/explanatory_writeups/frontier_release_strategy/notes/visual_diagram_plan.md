# Visual Diagram Plan

## Visual Language

Green means passed, amber means pending gate, red means reject and revise, and
gray means not yet attempted. Schematics are decision aids, not evidence.

## Opening Cartoon

Show the validation ladder from oracle through Andes tests to the complete
Frontier golden snapshot, named canaries, and production. Put the frozen
artifact and immutable validator report beside the golden gate. The gate must
have a visible failure path back to revision rather than forward to production.

Prompt: `diagram_prompts/golden_release_gate.md`.

## Risk Diagram

Contrast a partial shard test with a complete snapshot. The partial side can
test parsing and local reduction. Only the complete side can test contiguous
identity, global leaf uniqueness, ancestor exclusion, and volume closure.

Prompt: `diagram_prompts/partial_vs_complete.md`.

## Evidence Pair

Use `validation_ladder.pdf` to mark the exact untested boundary and
`production_inventory.pdf` to show the 161-snapshot exposure, including the 21
snapshots without archived controls.
