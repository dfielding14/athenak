# Example: Localized and AMR Turbulence Driving

These examples use the `turb` problem generator with hydro and write the
applied acceleration field through `variable = turb_force`.  They are short
demonstration runs rather than statistically converged turbulence experiments.

## Build

From the repository root:

```bash
cmake -S . -B build-turb -DPROBLEM=turb -DCMAKE_BUILD_TYPE=Release
cmake --build build-turb -j
```

Run each example from a separate output directory because AthenaK creates
`bin/` and history outputs in the current working directory:

```bash
mkdir -p run/turb_edot_fixed run/turb_tiled_include run/turb_amr_exclude
(cd run/turb_edot_fixed && ../../build-turb/src/athena \
  -i ../../inputs/hydro/turb_driving/constant_edot_fixed_grid.athinput)
(cd run/turb_tiled_include && ../../build-turb/src/athena \
  -i ../../inputs/hydro/turb_driving/tiled_localized_include.athinput)
(cd run/turb_amr_exclude && ../../build-turb/src/athena \
  -i ../../inputs/hydro/turb_driving/amr_localized_exclude.athinput)
```

## Cases

### Constant-Edot Fixed Grid

`constant_edot_fixed_grid.athinput` uses a \(32^3\) periodic grid and
`constant_edot = true`.  This is the baseline for a statistically homogeneous
force field normalized to a requested energy injection rate.

### Tiled Included Region

`tiled_localized_include.athinput` repeats one realization in
\(2\times2\times1\) tiles.  It selects `constant_edot = false`,
`accel_rms = 0.3`, and multiplies the field by a three-dimensional Gaussian
centred at the origin.  The resulting acceleration is concentrated in the
central region while the substructure records the tile scale.

### Adaptive Excluded Region

`amr_localized_exclude.athinput` enables adaptive refinement using a location
criterion about the origin.  It also uses a \(2\times2\times1\) tile pattern,
but applies `localization = exclude`; the central refined region is quiet and
the forcing is retained outside it.  In a validation run to cycle 12, the
input created 56 MeshBlocks and completed with 88 active MeshBlocks.

## Projected Acceleration

The figure below was generated from the final `turb_force` binary dump of each
included input.  Every panel displays
\(L_3^{-1}\int |\boldsymbol{a}|\,dx_3\).  The cyan outline in the AMR panel
marks the level-0/level-1 transition in the projected mesh.

![Projections of turbulence-driving examples](../_static/turbulence_driving_projections.png)

The short validation runs used for this figure produced the following checks:

| Case | Measured check |
| --- | --- |
| Constant-Edot fixed grid | Least-squares energy growth rate after forcing onset: `0.0992` for requested `dedt = 0.1`. |
| Tiled Gaussian include | Volume-weighted `a_rms = 0.29999999`; requested `accel_rms = 0.3`; centre/outer projected-amplitude ratio `13.67`. |
| AMR Gaussian exclude | `88` final MeshBlocks at maximum physical level `1`; centre/outer projected-amplitude ratio `0.56`. |

For the localized runs, the centre/outer ratio is evaluated on the cycle-12
projection shown above: the mean in `r < 0.12` divided by the mean in
`r > 0.32`, where `r = sqrt(x1^2 + x2^2)`.

Regenerate the figure after running the examples:

```bash
python vis/python/plot_turbulence_driving_projection.py \
  --panel "Constant Edot=run/turb_edot_fixed/bin/turb_edot_fixed.force.00004.bin" \
  --panel "Tiled + Gaussian include=run/turb_tiled_include/bin/turb_tiled_include.force.00004.bin" \
  --panel "AMR + Gaussian exclude=run/turb_amr_exclude/bin/turb_amr_exclude.force.00004.bin" \
  --output docs/source/_static/turbulence_driving_projections.png
```

## Adapting the Inputs

- Use `constant_edot = true` when the scientific control variable is the
  injected power density; use `constant_edot = false` with an explicit
  `accel_rms` when the acceleration field itself is the controlled amplitude.
- Change `localization = include` to concentrate driving in a cloud or disc,
  and `localization = exclude` to leave a protected central region undriven.
- Keep tile counts divisible into the root-grid cell counts.  Tiling is
  evaluated in physical coordinates and is compatible with refined blocks.
- For adaptive runs, choose `max_nmb_per_rank` using the intended refinement
  ceiling rather than the initial number of blocks.

The detailed equations and complete parameter reference are in
{doc}`../modules/turbulence_driving`.
