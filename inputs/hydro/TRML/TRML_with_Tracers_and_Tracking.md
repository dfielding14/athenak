# Simple TRML with tracers and frame tracking

Build this problem with:

```bash
cmake -S . -B build_trml -DPROBLEM=simple_TRML
cmake --build build_trml -j
```

The companion input `TRML_with_Tracers_and_Tracking.athinput` combines four pieces:

- `simple_TRML.cpp` supplies the pressure-balanced shear layer, exact cooling update,
  x3 reservoirs, and user history diagnostics.
- `<initial_perturbations>` supplies the one-time, reproducibly seeded velocity field.
  The pgen-local perturbation amplitude is zero to prevent applying two perturbations.
- `<frame_tracking>` follows the hot-side edge of the cold-material scalar
  (`0.05 <= scalar0 <= 0.2`) along x3. Its PD derivative uses the measured
  interface-position rate rather than the velocity of gas flowing through the interface.
- `particle_type=lagrangian_mc` follows the mass flux and samples thermodynamic fields.

## Reproducible local-run inputs

These standalone inputs record the complete configurations used for the local
parameter runs, without requiring command-line parameter overrides:

- `TRML_chi56p234_mach0p5_xi1e3_t150_16x16x32.athinput`: `xi=1000`,
  `16 x 16 x 32`, through `t=150`.
- `TRML_chi56p234_mach0p5_xi1e2_t75_48x48x96.athinput`: `xi=100`,
  `48 x 48 x 96`, through `t=75`.
- `TRML_chi56p234_mach0p5_xi1e2_t75_48x48x96_zero_gradient_vx.athinput`:
  the same `xi=100` run with zero-gradient `vx` at both x3 reservoirs.

All three use eight root MeshBlocks for one-block-per-rank launches on eight MPI ranks,
and retain the same 3,072-particle split and 51-event injection schedule.

## X3 reservoir velocity condition

Both x3 boundaries always hold density, pressure, and cold-material fraction at their
reservoir values. By default, `problem/zero_gradient_vx=true` copies `vx` from the
boundary-adjacent active cell into every ghost layer. Set it to `false` to hold `vx`
at the imposed shear values instead. The `vy` and `vz` treatment is unchanged. The
standalone inputs for the completed fixed-`vx` runs keep `false` explicitly so those
results remain reproducible.

## Frame and particle coordinates

Lagrangian Monte-Carlo tracers store grid-frame cell-center positions. They do not carry
particle velocities. A frame-controller update is an instantaneous velocity-coordinate
change, so it must not translate these particle positions or write the unused drift-
particle velocity slots. The after-integrator task dependency is instead:

1. advance the frame displacement and apply the fluid boost;
2. move Monte-Carlo tracers with the saved mass fluxes from the completed timestep;
3. migrate particles and seed any events due at the new time.

The boosted fluid fluxes naturally control tracer motion on the following timestep. This
avoids applying the frame change twice.

`prtcl_thermo_history` records grid-frame coordinates. The aligned frame history records
`ft_dx_x3` and `ft_vf_x3`. Reconstruct lab-frame quantities at a common time with

```text
x3_lab = x3_grid + ft_dx_x3
v3_lab = v3_grid + ft_vf_x3
```

The canonical input gives the history and particle-history outputs the same cadence so
they can be joined exactly by time or cycle. Velocity-carrying particle species have not
been validated with this frame tracker; the statement above is specific to
`lagrangian_mc` tracers.

## Tracer population and injection

The canonical run contains 3,072 tracers in total:

- seed ID 1 places 777 tracers uniformly by volume throughout the initial domain at
  `t=0` (25.3% of the total population);
- seed ID 2 places 45 tracers uniformly by area in the top boundary slab at
  each of 51 times, `t=0, 0.4, ..., 20` (2,295 tracers, or 74.7%).

The top schedule uses the slab `0.984375 <= x3 <= 1.0`, which is one root-grid cell
thick on the committed 128-cell x3 mesh. In the default uniform run this is the
uppermost active-cell layer; if that slab is refined it contains the corresponding fine
cells. Particles are seeded at active-cell centers, not in ghost cells. If `time/tlim`,
the x3 domain, or the root x3 resolution is changed, update `tracer_seed2/end_time`,
`slab_min`, and (if needed) `cadence` together. Keeping
`(end_time - start_time) / cadence` integral includes both endpoints.

## Passive scalar convention

The pgen stores cold-material fraction in `scalar0` as a conserved scalar:
`rho * cold_fraction`. The inner x3 reservoir supplies fraction one and the outer
reservoir supplies fraction zero. When frame tracking is active, both reservoirs are
transformed from their lab velocities using the current frame velocity.

The canonical frame tracker volume-weights cells with cold fraction from 0.05 to 0.2.
This band marks the visible hot-side density/temperature front and is present in the
smoothed initial condition. `velocity_signal=position_rate` differentiates the filtered
band centroid, avoiding the persistent offset produced by using the mean velocity of
gas that crosses a cooling interface.

## Restart convention

MPI runs with Monte-Carlo tracers must use `single_file_per_rank=true` for restart
output. Restart files preserve particle tags and seed schedules as well as the complete
frame-controller state. Initial perturbations are applied only to new runs and are not
replayed after restart.

The staged HIP/MPI validation procedure is in
[`docs/TRML_WITH_TRACERS_AND_TRACKING_FRONTIER_HANDOFF.md`](../../../docs/TRML_WITH_TRACERS_AND_TRACKING_FRONTIER_HANDOFF.md).
