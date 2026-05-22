# Source Terms

This directory contains AthenaK source-term implementations.

`initial_perturbations.hpp/cpp` provides one-time initial-condition
perturbations from an `<initial_perturbations>` input block. It can perturb
density, velocity, or face-centered magnetic fields when the corresponding RMS
amplitude is positive. Magnetic perturbations are generated from an
edge-centered vector potential and added as a discrete curl so
constrained-transport `div(B)` is preserved to roundoff.

The regression example in `inputs/tests/initial_perturbations.athinput` writes a
full-grid VTK snapshot that can be turned into documentation figures with
`scripts/plot_initial_perturbations_example.py`.
