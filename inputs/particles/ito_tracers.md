# Second-Moment Itô Mass-Flux Tracers

The runnable inputs in this directory exercise the Itô-2 implementation:

| Input | Purpose |
| --- | --- |
| `ito_tracers.athinput` | Uniform isothermal advection used to verify the analytical drift and variance, RK weighting, restart behavior, and continuous trajectories. |
| `ito_tracers_amr.athinput` | Ideal-gas AMR and MPI smoke test for coefficient exchange, migration, restart output, and preservation of subcell positions. |

Enable the method with:

```ini
<particles>
particle_type = lagrangian_ito
pusher        = ito2
```

Only Itô-2 is implemented. `pusher = ito3`, `ito_order` other than `2`, and
non-uniform kick distributions are rejected explicitly. See
`docs/source/modules/ito_tracers.md` for the equations, supported configuration,
and implementation details.
