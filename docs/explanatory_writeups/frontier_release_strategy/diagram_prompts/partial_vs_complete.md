# Partial Versus Complete Validation

**Purpose:** Make clear why eight successful shards cannot certify a 1024-shard snapshot.

**Image-generation prompt:**

> Side-by-side comparison titled “Partial test” and “Golden whole-input test”. Left panel: eight highlighted shards out of a field of 1024 boxes feed a reducer. Checks for parsing, GPU kernel, MPI combination, and matched archived comparison are green; contiguous shard identity, global logical-key checks, and summed-volume checks are gray. Right panel: all 1024 shards highlighted, those operational checks active, outputs compared to archived controls, plus an amber note that geometry-to-key validation is still required to prove arbitrary physical coverage. At bottom: “Partial success says the method works on sampled input. Whole-input success passes the implemented global invariants; it is not universal proof.”

**Caption:** Partial validation exercises the computational path; whole-input validation activates the implemented global AMR invariants.

**How to read it:** Look at the gray checks on the partial side. Those are the remaining common-mode risks.

**Evidence status:** Schematic, not evidence.
