# Module: Particles

The public particles module is currently a narrow Lagrangian drift path, not
a general catalogue of particle species and pushers.

## Public Interface

`src/particles/particles.cpp` constructs particles only when `<particles>` is
present and accepts a two- or three-dimensional mesh at initialization. The
publicly supportable path documented here is a three-dimensional, strictly
periodic drift run: the shipped pusher and migration path access z state
unconditionally and wrap periodic boundaries without general physical-boundary
handling.

| Parameter | Requirement/default | Supported public value |
| --- | --- | --- |
| `particle_type` | Required | `cosmic_ray` |
| `pusher` | Required | `drift` |
| `ppc` | Default `1.0` | Particles per cell used to allocate the pack |
| `assign_tag` | Default `index_order` | `index_order` or `rank_order` |

Values such as `star` particle types, Boris pushers, gravity pushers, and
species/mass-spectrum controls are not accepted by the public constructor.

## Shipped Example

`inputs/particles/random_particle_drift.athinput` uses:

```ini
<particles>
particle_type = cosmic_ray
ppc = 0.01
pusher = drift

<output1>
file_type = pvtk
variable = prtcl_all
dt = 0.01
```

The deck is three-dimensional and also defines mesh particle-density output
with `variable = prtcl_d`.

Do not treat two-dimensional meshes or non-periodic physical boundaries as
validated particle workflows.

## Implementation Map

| Source | Role |
| --- | --- |
| `src/particles/particles.hpp`, `particles.cpp` | State allocation, accepted inputs, tag creation |
| `src/particles/particles_pushers.cpp` | Implemented drift update |
| `src/particles/particles_tasks.cpp` | Push and communication tasks |
| `src/bvals/bvals_part.cpp` | Particle boundary communication |

The documented three-dimensional path allocates two integer arrays and six
real components per particle for the accepted `cosmic_ray` drift workflow.

## See Also

- [Outputs](outputs.md)
- [Mesh](mesh.md)
