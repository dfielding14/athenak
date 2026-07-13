# Historical extended TRML with tracers and frame tracking

Current minimal-pgen documentation lives in [`README.md`](README.md). This file
records the older tracer/frame-tracking workflow preserved in
`simple_TRML_extended.cpp`; it is not the reference for `simple_TRML` or the
current `TRML_with_Tracers_and_Tracking.athinput`.

Build this problem with:

```bash
cmake -S . -B build_trml -DPROBLEM=simple_TRML_extended
cmake --build build_trml -j
```

The historical configuration combined four pieces:

- `simple_TRML_extended.cpp` supplies the pressure-balanced shear layer, exact cooling update,
  upper x3 hot reservoir, and user history diagnostics.
- `<initial_perturbations>` supplies the one-time, reproducibly seeded velocity field.
  The pgen-local perturbation amplitude is zero to prevent applying two perturbations.
- `<frame_tracking>` follows the hot-side edge of the cold-material scalar
  (`0.05 <= scalar0 <= 0.2`) along x3. Its PD derivative uses the measured
  interface-position rate rather than the velocity of gas flowing through the interface.
- `particle_type=lagrangian_mc` follows the mass flux and samples thermodynamic fields.

## Historical input

The historical integrated run used the `xi=100`, `48 x 48 x 96`,
50-shear-time setup,
including eight root MeshBlocks for one-block-per-rank launches on eight MPI ranks,
the low-bandwidth controller, single-precision full-volume and slice outputs, and
the 3,072-particle/51-event tracer schedule.

## Low-bandwidth frame control

The frame should remove secular interface drift, not force every resolved turbulent
excursion back to exactly zero. The stable input uses `tau_avg=1`, `tau_relax=5`,
`tau_vel=2.5`, and `max_boost_change_rate=0.02` from `t=0`. These ratios give a
nominally critically damped PD response while placing its bandwidth below the resolved
turbulent fluctuations. Proportional feedback naturally strengthens as the interface
moves farther from center, so the same controller can acquire the initial interface
without a timed gain switch, mode transition, or controller-memory reset.

There is no parameter-free way to distinguish secular drift from turbulence: that is a
choice of scale. For another problem, choose one tracking timescale that is longer than
the turbulent correlation time but shorter than the time for the interface to approach a
boundary, then keep it fixed for the entire run. A useful critically damped starting
ratio is `tau_avg = 0.2 tau_relax` and `tau_vel = 0.5 tau_relax`. This is a bandwidth
choice, not a simulation-phase or absolute-time trigger.

In the `xi=100`, `48 x 48 x 96` validation used to choose this controller, the
tracked-band offset peaked at 0.136 and returned within one cell of center by
`t=19.5`. Over `t=20.25--74.75`, its detrended cooling standard deviation was
0.0163 versus 0.0798 for the aggressive controller, and the fraction of cooling
power in the controller-artifact band at frequencies 1.5--2.0 fell from 80.8% to
2.35%.

## X3 boundary conditions

Both x3 boundaries use frame-aware user reservoirs. The lower reservoir supplies
fixed cold gas with fixed pressure and cold-material fraction one; the upper reservoir
supplies fixed hot gas with fixed pressure and cold-material fraction zero. The
canonical `problem/zero_gradient_vx=false` holds `vx` at the imposed cold- and
hot-side shear values. Set it to `true` to copy `vx` from each boundary-adjacent
active cell instead.

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

See [TRML_TRACER_PARTICLE_OUTPUT_GUIDE.md](TRML_TRACER_PARTICLE_OUTPUT_GUIDE.md)
for the complete column reference, MPI/restart behavior, integrity checks, and
temperature-evolution plotting examples.

The historical input gave the history and particle-history outputs the same cadence so
they can be joined exactly by time or cycle. Velocity-carrying particle species have not
been validated with this frame tracker; the statement above is specific to
`lagrangian_mc` tracers.

## Tracer population and injection

The canonical run contains 3,072 tracers in total:

- seed ID 1 places 777 tracers uniformly by volume throughout the initial domain at
  `t=0` (25.3% of the total population);
- seed ID 2 places 45 tracers uniformly by area in the top boundary slab at
  each of 51 times, `t=0, t_shear, ..., 50 t_shear`, where
  `t_shear=1.5491933384829668` (2,295 tracers, or 74.7%).

The top schedule uses the slab `0.9791666666666666 <= x3 <= 1.0`, which is one
root-grid cell thick on the committed 96-cell x3 mesh. In the default uniform run this
is the uppermost active-cell layer; if that slab is refined it contains the corresponding
fine cells. Particles are seeded at active-cell centers, not in ghost cells. If `time/tlim`,
the x3 domain, or the root x3 resolution is changed, update `tracer_seed2/end_time`,
`slab_min`, and (if needed) `cadence` together. Keeping
`(end_time - start_time) / cadence` integral includes both endpoints.

## Passive scalar convention

The pgen stores cold-material fraction in `scalar0` as a conserved scalar:
`rho * cold_fraction`. The lower reservoir supplies fraction one, while the upper
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
