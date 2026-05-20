# Cloud-Crushing Problem Generator

`src/pgen/cloud_crushing.cpp` is a custom problem generator for a cold cloud in
a warm ISM. It is intended to be built with:

```bash
cmake -S . -B build_cloud_crushing -DPROBLEM=cloud_crushing
```

The setup requires a `<hydro>` block with `ism_cooling = true` and a `<units>`
block so the ISM equilibrium and Sedov boundary can be evaluated in cgs-backed
code units.

The problem generator supports two inner-`x1` boundary modes:

| Mode | Description |
|------|-------------|
| `sedov` | Samples a Taylor-von Neumann-Sedov blast profile at the left boundary. |
| `constant` | Fills the left boundary with a constant user-specified state. |

When `<frame_tracking>` is enabled, the boundary uses `FrameDisplacement(axis)`
to convert grid-frame sample points back to lab-frame positions and
`FrameVelocity(axis)` to transform lab-frame velocities into the moving grid
frame.
