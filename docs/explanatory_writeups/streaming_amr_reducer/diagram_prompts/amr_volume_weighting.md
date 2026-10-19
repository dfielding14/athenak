# AMR Volume Weighting

**Purpose:** Explain why native leaves can be combined without uniform-grid resampling.

**Image-generation prompt:**

> Scientific cartoon of a two-dimensional AMR mesh standing in for a three-dimensional cube. One coarse region is refined into smaller cells. Arrows from cells enter histogram bins, with arrow thickness proportional to physical cell volume. Show that refined children together occupy the same volume as their replaced parent, and that parent and children must never both be counted. Include equations H_b += Delta V_i w_i and sum leaf Delta V_i = domain volume. Add a red overlap warning for parent plus children.

**Caption:** Every accepted AMR leaf contributes according to physical volume; completeness and non-overlap make the sum correct.

**How to read it:** Small refined cells contribute less individually. Their total replaces, rather than supplements, the parent.

**Evidence status:** Schematic, not evidence.
