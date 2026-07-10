# Mixed Survival Timeline

**Purpose:** Explain why output29 retains a radial marginal while output31 remains intact.

**Image-generation prompt:**

> Three-lane technical timeline. First lane: twelve radius-first products compute radius, then coord_abscostheta performs an inner realloc from nmb_max to nmb and invalidates radius, producing observed radius underflow. Second lane: output29 computes coord_abscostheta first in the shrunken nmb buffer; the next radius call's common capacity check reallocates back to nmb_max, invalidates angle, and then writes radius, so the radial marginal survives. Third lane: output31 never calls coord_abscostheta and ends in an intact four-dimensional negative-control histogram. Green checks mark validated quantities, red marks invalidated coordinates, amber marks a marginal that remains usable. White background, compact scientific style.

**Caption:** Axis order and coordinate dependencies produce a mixed archive rather than one uniform failure.

**How to read it:** A correct marginal does not imply a correct joint distribution.

**Evidence status:** Schematic, not evidence.
