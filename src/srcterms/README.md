# Source Terms

This directory contains AthenaK source-term implementations.

`initial_perturbations.hpp/cpp` provides one-time initial-condition
perturbations from an `<initial_perturbations>` input block. It can perturb
density, velocity, or face-centered magnetic fields when the corresponding RMS
amplitude is positive. Magnetic perturbations are generated from an
edge-centered vector potential and added as a discrete curl so
constrained-transport `div(B)` is preserved to roundoff.

The regression examples in `inputs/tests/initial_perturbations.athinput` and
`inputs/tests/initial_perturbations_2d.athinput` write full-grid VTK snapshots.
The snapshots can be turned into density, velocity-decomposition, magnetic, and
`divB` documentation figures with `scripts/plot_initial_perturbations_example.py`.
