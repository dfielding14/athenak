# Full single-precision application build: unresolved portability work

A fresh AppleClang CPU Serial build with `Athena_SINGLE_PRECISION=ON` and no
narrowing-warning suppression reproduced the pre-existing build failure.
Correcting the first type errors uncovered further blockers in unrelated
relativistic EOS, unit conversion, and geodesic-grid code. All trial repairs
were removed; production files retain their pre-review contents.

The cascade includes:

- Narrowing aggregate initializers in `coordinates/cartesian_ks.hpp`,
  `dyn_grmhd/dyn_grmhd.cpp`, `geodesic-grid/geodesic_grid.cpp`, and
  `eos/primitive-solver/unit_system.cpp`.
- Mixed float/double arguments to `std::min` in `diffusion/hyperviscosity.cpp`.
- Table-reader `double*` values assigned to `Real*` in `eos_compose.cpp` and
  `eos_hybrid.cpp`.
- `double[]` inputs passed to a `Real*` interface in `piecewise_polytrope.cpp`.

Unit conversions also require a float range audit: merely silencing narrowing
would not validate the large CGS products or tiny conversion factors. This is
broader than the CGL fixes and remains a validation limitation, not a passing
full-float build.

The configure and three compiler logs retain the successive failures. The
`float-partial-repair-not-retained.patch` is diagnostic evidence only. The
already passing standalone float CGL math tests do not substitute for this
application-level build.

Reproduction from the repository root:

```sh
cmake -S . -B tst/build-single -DAthena_SINGLE_PRECISION=ON -DCMAKE_CXX_FLAGS=
cmake --build tst/build-single -j 6
```
