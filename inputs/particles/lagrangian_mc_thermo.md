# Lagrangian Monte-Carlo Thermodynamic Tracers

`particle_type = lagrangian_mc` adds Monte-Carlo tracer particles that jump between cell
centers with probabilities set by the RK-accumulated density flux through the faces of the
current cell. The pusher is intended for ideal-gas Hydro and MHD runs.

Enable the particle type with:

```ini
<particles>
particle_type = lagrangian_mc
pusher        = lagrangian_mc
random_seed   = 12345
```

Tracer seeding is controlled by one or more `<tracer_seedN>` blocks. Each block defines a
schedule, a spatial region, an optional thermodynamic mask, and deterministic sampling
seed:

```ini
<tracer_seed1>
id              = 1
start_time      = 0.0
end_time        = 0.1
cadence         = 0.01      # <= 0 means one-shot
count_per_event = 1000
weight          = mass      # mass or volume
region          = box       # all, box, sphere, slab
x1min           = -0.5
x1max           =  0.5
target          = temperature
target_min      = 1.0e-2
target_max      = 1.0e2
seed            = 24680
```

`sphere` regions use `center1`, `center2`, `center3`, and `radius`; `slab` regions use
`slab_axis`, `slab_min`, and `slab_max`. `target` may be `density`, `temperature`,
`pressure`, `entropy`, or `scalarN`. New particles are inserted at selected cell centers.
Initial seeding fires after the problem generator and primitive-variable initialization,
before initial outputs. Timed seeding fires after the fluid update, existing tracer
motion, and particle communication.

Use `file_type = prtcl_thermo_history` for append-only path histories:

```ini
<output1>
file_type = prtcl_thermo_history
dt        = 0.01
```

The binary records contain `time, cycle, tag, seed_id, x1, x2, x3, gid, rho, p, T, s,
eint, v1, v2, v3, scalars`. Read them with:

```bash
python scripts/read_prtcl_thermo_history.py prtcl_thermo_history/<basename>.prtcl_thermo_history.thp --npz tracers.npz
```

Restart files persist particle arrays, `next_tracer_tag`, and each schedule's next fire
time. In MPI runs, use `single_file_per_rank = true` for restart outputs with
`lagrangian_mc` particles.
