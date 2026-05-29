# Example: Blast Wave

The shipped blast inputs exercise the custom problem generator
`src/pgen/blast.cpp` with hydro or MHD, with uniform-grid and adaptive-mesh
variants.

## Inputs

| Input deck | Physics | Mesh behavior | Output |
| --- | --- | --- | --- |
| `inputs/hydro/blast_hydro.athinput` | Hydro | Uniform two-dimensional grid | History and `hydro_w` VTK |
| `inputs/hydro/blast_hydro_amr.athinput` | Hydro | Adaptive mesh using an `<amr_criterion0>` density-slope criterion | History and binary `hydro_w` |
| `inputs/mhd/blast_mhd.athinput` | MHD | Uniform-grid magnetic blast | History and VTK |
| `inputs/mhd/blast_mhd_amr.athinput` | MHD | Adaptive mesh magnetic blast | See the shipped deck |

## Build And Run

`blast.cpp` is a build-selected custom generator; it is not in the default
`pgen_name` dispatcher.

```bash
cmake -S . -B build-blast -DPROBLEM=blast
cmake --build build-blast

# Short public hydro validation run
./build-blast/src/athena -i inputs/hydro/blast_hydro.athinput \
  -d run-blast time/nlim=1
```

This one-cycle command was validated against the public source tree and
writes:

```text
run-blast/Blast.hydro.hst
run-blast/vtk/Blast.hydro_w.00000.vtk
run-blast/vtk/Blast.hydro_w.00001.vtk
```

For the complete supplied run, omit `time/nlim=1`. Switch to one of the MHD
or AMR input decks without rebuilding; they use the same selected generator.

## Verified Controls

The generator reads required `<problem>/inner_radius`,
`<problem>/outer_radius`, and `<problem>/prat`. It also accepts public optional
controls including `pn_amb`, `dn_amb`, `pi_amb`, `di_amb`, `drat`, `b_amb`,
and `coordinates`, as defined in `src/pgen/blast.cpp`.

The uniform hydro input sets `nghost = 3` because it uses
`reconstruct = ppm4`. The AMR hydro input instead uses `plm` and defines:

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

## Validation Use

These decks are appropriate for checking shock propagation, symmetry and AMR
behavior. Quantitative Sedov-Taylor comparison requires defining the intended
energy normalization and measuring the output; it is not asserted as an
automatic result of the supplied run.

## See Also

- [Configuration](../configuration.md)
- [Problem Generators](../modules/pgen.md)
- [Visualization Utilities](../tools/visualization.md)
