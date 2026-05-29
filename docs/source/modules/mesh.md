# Module: Mesh

The mesh layer constructs the root domain, partitions it into MeshBlocks,
tracks static/adaptive refinement, and owns the `MeshBlockPack` containers
through which physics modules operate. The public implementation is in
`src/mesh/`.

## Public Objects

| Source | Public responsibility |
| --- | --- |
| `mesh.hpp`, `mesh.cpp` | `Mesh`, root-domain parsing, boundary validation, physics attachment |
| `meshblock.hpp`, `meshblock.cpp` | Per-pack MeshBlock IDs, levels, extents, and neighbors |
| `meshblock_pack.hpp`, `meshblock_pack.cpp` | Coordinates, physics pointers, and task lists for a pack |
| `meshblock_tree.*`, `build_tree.cpp` | Static/refined tree construction |
| `mesh_refinement.*` | Adaptive refinement driver and transfer coordination |
| `refinement_criteria.*` | Built-in adaptive criteria parsed from `<amr_criterion*>` |
| `load_balance.cpp` | MeshBlock redistribution |

`MeshBlockPack::AddPhysics()` creates a module only when its input block is
present, for example `<hydro>`, `<mhd>`, `<radiation>`, or `<particles>`.

## Root Mesh And MeshBlocks

Every run supplies `<mesh>`. Required domain fields are `nx1`, `nx2`, `nx3`,
`x1min`/`x1max`, `x2min`/`x2max`, `x3min`/`x3max`, and boundary flags for
active directions. `nghost` defaults to `2`.

```ini
<mesh>
nx1 = 256
nx2 = 1
nx3 = 1
nghost = 2
x1min = -1.0
x1max = 1.0
x2min = -0.5
x2max = 0.5
x3min = -0.5
x3max = 0.5
ix1_bc = periodic
ox1_bc = periodic

<meshblock>
nx1 = 64
```

`<meshblock>/nx*` sets block cell dimensions; omitted values inherit the
corresponding root values. Each block dimension must divide the global mesh
dimension. Reducing block dimensions at fixed global resolution changes
decomposition and ghost-zone overhead, not physical resolution.

## Mesh Checks

The constructor enforces these constraints:

| Rule | Consequence |
| --- | --- |
| Each active direction has at least 4 cells | `nx1 >= 4`, plus `nx2`/`nx3 >= 4` when active |
| An `x1-x3` two-dimensional mesh is unsupported | `nx2 = 1`, `nx3 > 1` is rejected |
| `nghost >= 2` | Increase it for high-order reconstruction or FOFC |
| SMR/AMR requires even `nghost` | Refined cases that need three ghost zones must use at least four |
| Periodic boundaries are paired on a direction | Both faces must be periodic together |
| `shear_periodic` is x1-only and paired | It also requires a `<shearing_box>` block and is incompatible with mesh refinement |

For hydro and MHD reconstruction-specific ghost requirements, see
[Configuration](../configuration.md).

## Static And Adaptive Refinement

`<mesh_refinement>/refinement` accepts `none`, `static`, or `adaptive`.

| Parameter | Use |
| --- | --- |
| `num_levels` | Adaptive only; maximum adaptive refinement-level count, default `1` |
| `max_nmb_per_rank` | Required for adaptive refinement; caps allocated MeshBlocks per rank |
| `ncycle_check` | Adaptive checking cadence; default `1` |
| `refinement_interval` | Minimum cycles between adaptive changes; default `5` |
| `prolong_primitives` | Multilevel transfer option; use primitive rather than conserved variables in SMR or AMR prolongation when `true` |

Static refinement regions are specified with `<refined_region*>` blocks; see a
shipped example such as `inputs/grhydro/gr_fm_torus_smr.athinput`. Static
refinement derives its constructed levels from those regions and does not
require `num_levels` or `max_nmb_per_rank`.

Adaptive refinement requires one or more `<amr_criterion*>` blocks. The
initial tree may also contain `<refined_region*>` blocks before dynamic
criterion checks begin. The current criteria parser accepts:

| `method` | Additional fields |
| --- | --- |
| `min_max`, `slope`, `second_deriv` | `variable`; optional `value_min`, `value_max` |
| `location` | Optional `location_x1`, `location_x2`, `location_x3`, `location_rad` |
| `user` | Selected problem generator must enroll `user_ref_func` |

For example, the public adaptive blast input uses:

```ini
<mesh_refinement>
refinement = adaptive
num_levels = 2
refinement_interval = 3
max_nmb_per_rank = 1024

<amr_criterion0>
method = slope
variable = hydro_w_d
value_max = 0.1
```

The older threshold keys `dens_max`, `ddens_max`, `dpres_max`, and
`dvel_max` are not the interface parsed by the public refinement criteria
implementation.

## Problem-Generator Hook

Problem generators may enroll a custom adaptive criterion through
`ProblemGenerator::user_ref_func`. Public source uses that hook in several
relativity/problem-generator files. An `<amr_criterion*>` with
`method = user` invokes that enrolled function; it is not an automatic
replacement for supplying an implementation.

## See Also

- [Public Input Parameters](../reference/input_parameters.md)
- [Boundaries](boundaries.md)
- [Problem Generators](pgen.md)
