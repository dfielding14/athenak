# Shared Buffer Erasure

**Purpose:** Show how a valid set of derived coordinates was silently erased while calculating `coord_abscostheta`.

**Image-generation prompt:**

> Clean technical cartoon on a white background, designed like a blackboard explanation. Show two horizontal timelines. Timeline A, “radius first”: radius is written into a shared buffer sized nmb_max; a red inner Kokkos::realloc for coord_abscostheta shrinks it to nmb and invalidates radius; the angle row is written. Timeline B, “angle first / output29”: the inner realloc writes angle into an nmb-sized buffer; the next coordinate call notices nmb is smaller than nmb_max, performs a second red capacity realloc, invalidates angle, and then writes radius. Below, show correct weights entering wrong coordinate bins. Use blue for valid values, red for destructive reallocation, gray for invalidated rows. Clearly label that non-preservation is guaranteed while the observed archived result is zero/underflow. Minimal labels, scientific diagram, no decoration.

**Caption:** The corruption was caused by resizing shared derived-coordinate storage, not by bad simulation state.

**How to read it:** Follow the rows through the red reallocation. The weights still exist, but earlier coordinates do not.

**Evidence status:** Schematic, not evidence. The evidence is the source diff and archived/rebuilt comparison.
