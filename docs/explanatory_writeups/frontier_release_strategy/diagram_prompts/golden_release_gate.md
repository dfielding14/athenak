# Golden Release Gate

**Purpose:** Explain why the complete Frontier snapshot sits between partial tests and production.

**Image-generation prompt:**

> Horizontal staged-release diagram for an HPC reconstruction campaign. Six gates from left to right: separate synthetic tests, Andes K80 partial real-data tests, frozen executable plus exact 1024-shard 8 pc golden snapshot, MI250X performance profile with predeclared budgets, 256-node and no-control canaries, full 161-snapshot production array. Green checks on completed early gates, amber locks on pending gates, gray production array waiting at the end. Show a red rollback arrow from any failed gate back to code revision. Clearly show that production submission requires an immutable passing validator report tied to the frozen executable. Include labels: “agreed formulas?”, “real file format?”, “whole-input invariants and archived controls?”, “fast enough at scale?”, “variant canaries?”, “release”.

**Caption:** The golden run is the first whole-input test; the canaries cover scale and no-control variants before production.

**How to read it:** Each stage removes a distinct risk. Passing an earlier stage does not answer the next question.

**Evidence status:** Schematic, not evidence.
