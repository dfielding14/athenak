# Frame-Tracking Example Problem Generators

`src/pgen/cloud_crushing.cpp` is a custom problem generator for a cold cloud in
a warm ISM. It is intended to be built with:

```bash
cmake -S . -B build_cloud_crushing -DPROBLEM=cloud_crushing
```

The setup requires a `<hydro_srcterms>` block with `ism_cooling = true` and a
`<units>` block so the ISM equilibrium and Sedov boundary can be evaluated in
cgs-backed code units.

The problem generator supports two inner-`x1` boundary modes:

| Mode | Description |
|------|-------------|
| `sedov` | Samples a Taylor-von Neumann-Sedov blast profile at the left boundary. |
| `constant` | Fills the left boundary with a constant user-specified state. |

When `<frame_tracking>` is enabled, the boundary uses `FrameDisplacement(axis)`
to convert grid-frame sample points back to lab-frame positions and
`FrameVelocity(axis)` to transform lab-frame velocities into the moving grid
frame.

`src/pgen/TRML_frame_tracking.cpp` is a compact turbulent radiative mixing-layer
example. It is intended to be built with:

```bash
cmake -S . -B build_trml_frame_tracking -DPROBLEM=TRML_frame_tracking
```

The setup initializes a pressure-balanced hot/cold interface, optional local
top-hat cooling, a passive scalar marking the cold phase, and optional user
`x3` boundaries. Those user boundaries use `FrameDisplacement()` to sample the
same lab-frame layer as the grid frame moves and subtract `FrameVelocity()` from
the lab-frame boundary velocity.

# Self-Gravity Test Problem Generators

The self-gravity feature adds built-in problem generators under
`src/pgen/tests/`:

| `pgen_name` | Source | Purpose |
|-------------|--------|---------|
| `gravity` | `jeans_wave.cpp` | Periodic Jeans-wave setup for hydro or MHD self-gravity. |
| `binary_gravity` | `binary_gravity.cpp` | Two dense spheres for isolated-potential tests. |
| `be_collapse` | `be_collapse.cpp` | Bonnor-Ebert-like collapse with Jeans AMR refinement. |

The lightweight regression inputs are `inputs/tests/selfgravity.athinput` and
`inputs/tests/selfgravity_mhd.athinput`. The full solver and parameter notes are
in [Self-Gravity and Multigrid](self_gravity.md).
