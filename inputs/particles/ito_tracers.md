# Second-Moment Itô Mass-Flux Tracers

The runnable inputs in this directory exercise the Itô-2 implementation:

| Input | Purpose |
| --- | --- |
| `ito_tracers.athinput` | Uniform isothermal advection used to verify the analytical drift and variance, RK weighting, restart behavior, and continuous trajectories. |
| `ito_tracers_amr.athinput` | Ideal-gas AMR and MPI smoke test for coefficient exchange, migration, restart output, and preservation of subcell positions. |
| `ito_tracers_1d_sheet.athinput` | Thin-domain sheet advection used for the 1D-style trajectory and displacement-distribution figure. |
| `ito_tracers_2d_cloud.athinput` | Compact cloud under diagonal advection used for the 2D continuous-path figure. |

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

## Behavior Figures

The 1D-style sheet test shows continuous trajectories and an ensemble
displacement distribution with the MC jump process's drift and variance:

![Ito-2 1D-style sheet behavior](../../docs/source/modules/figures/ito_tracers_1d_sheet.png)

The 2D cloud test shows diagonal drift, stochastic spreading, and continuous
subcell paths across the mesh:

![Ito-2 2D cloud behavior](../../docs/source/modules/figures/ito_tracers_2d_cloud.png)
