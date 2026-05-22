# Module: Particles

## Role in AthenaK

The Particles module implements Lagrangian cosmic-ray tracer particles on top of
the existing AthenaK task list, boundary exchange, and AMR infrastructure.  The
current supported particle type is `cosmic_ray`, with drift and magnetic Boris
pushers, per-particle species tags, restart output, compact position output,
pitch-angle histograms, displacement histograms, and AMR-aware particle remap.

Particle evolution is owned by `particles::Particles` and is scheduled through
the standard mesh task flow.  Particle exchange reuses the boundary-value
communication path, so particles that cross MPI-rank ownership after motion or
AMR are packed, sent, unpacked, and counted through the same rank metadata used
by MeshBlocks.

## Source Location

| Path | Purpose |
|------|---------|
| `src/particles/particles.hpp` | Particle class state, runtime options, and debug-check API. |
| `src/particles/particles.cpp` | Construction, restart sizing, AMR remap, and consistency checks. |
| `src/particles/particles_pushers.cpp` | Drift and Boris particle pushers. |
| `src/particles/particles_tasks.cpp` | Particle task registration. |
| `src/bvals/bvals_part.cpp` | Particle MPI exchange. |
| `src/outputs/*_prtcl.cpp` | Particle `ppd`, `prst`, `df`, `dxh`, VTK, and tracking outputs. |
| `src/pgen/part_random.cpp` | Standard random-particle test/problem generator. |

## Minimal Input

Enable cosmic-ray tracer particles with a `<particles>` block and select the
standard built-in particle generator with `<problem>/pgen_name = part_random`.

```text
<particles>
particle_type     = cosmic_ray
ppc               = 0.125
nspecies          = 2
pusher            = drift
check_consistency = true

<problem>
pgen_name      = part_random
prtcl_rst_flag = 0
```

For magnetic Boris pushing, add an `<mhd>` block and select an interpolation
scheme:

```text
<mhd>
eos         = ideal
reconstruct = plm
rsolver     = hlld
gamma       = 1.666666667

<particles>
particle_type     = cosmic_ray
ppc               = 0.125
nspecies          = 2
pusher            = boris
interpolation     = tsc
min_mass          = 1.0
mass_log_spacing  = 2.0
cfl_part          = 0.1
check_consistency = true

<problem>
pgen_name = part_random
B0z       = 1.0
```

The shipped smoke and stress inputs are:

| Input | Purpose |
|-------|---------|
| `inputs/particles/random_particle_drift.athinput` | Drift pusher smoke and restart tests. |
| `inputs/particles/cr_tracer_boris_amr.athinput` | 3-D periodic MHD Boris AMR smoke test. |
| `inputs/particles/cr_tracer_boris_amr_stress.athinput` | Moving refine/derefine MPI AMR stress test. |

## Runtime Parameters

### `<particles>`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `particle_type` | required | Must be `cosmic_ray`. |
| `ppc` | `1.0` | Particles per cell per species.  Fractional values are allowed. |
| `nspecies` | `1` | Number of particle species. |
| `pusher` | required | `drift` or `boris`. |
| `interpolation` | `tsc` | For `boris`: `tsc` or `lin`. |
| `min_mass` | `1.0` | Species-zero gyro parameter used by the Boris update. |
| `mass_log_spacing` | `1.0` | Multiplicative spacing between species gyro parameters. |
| `cfl_part` | `0.05` | Particle timestep cap used by `part_random`. |
| `assign_tag` | `index_order` | Initial tag assignment; `rank_order` is also supported. |
| `check_consistency` | `false` | Enables expensive host-side invariant checks for tests and debugging. |

### `<problem>` for `part_random`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pgen_name` | required for built-in builds | Use `part_random` for the standard pgen. |
| `prtcl_rst_flag` | `0` | If nonzero, load particle data from `prtcl_res_file`. |
| `prtcl_res_file` | required for particle restart | Per-rank `prst` path; MPI ranks replace `rank_00000000` with their rank directory. |
| `particle_refinement` | `none` | `moving_boxes` enrolls the AMR stress-test refinement pattern. |
| `B0x`, `B0y`, `B0z` | `0, 0, 1` | Uniform magnetic field initialized by `part_random`. |

## Particle Data Layout

Cosmic-ray tracer particles allocate 14 real fields and 3 integer fields.

Integer fields:

| Field | Meaning |
|-------|---------|
| `PGID` | Global MeshBlock ID currently owning the particle. |
| `PTAG` | Particle tag, unique within each species. |
| `PSP` | Particle species index. |

Real fields:

| Field | Meaning |
|-------|---------|
| `IPX`, `IPY`, `IPZ` | Particle position. |
| `IPVX`, `IPVY`, `IPVZ` | Particle velocity. |
| `IPM` | Particle gyro parameter. |
| `IPBX`, `IPBY`, `IPBZ` | Magnetic field sampled by the pusher. |
| `IPDX`, `IPDY`, `IPDZ` | Accumulated displacement. |
| `IPDB` | Accumulated displacement parallel to the sampled magnetic field. |

## Pushers

| Pusher | Requirements | Notes |
|--------|--------------|-------|
| `drift` | `<particles>` only | Advances particles with stored velocities. |
| `boris` + `tsc` | `<mhd>` block | Uses cell-centered magnetic fields and triangular-shaped-cloud weights. |
| `boris` + `lin` | `<mhd>` block | Uses linearly interpolated face-centered magnetic fields. |

The Boris pusher requires MHD because it samples magnetic fields.  The
interpolation choice is made through `<particles>/interpolation`.

## AMR Remapping

AMR changes the live leaf MeshBlock list and can also move MeshBlocks between
MPI ranks.  After AMR, `Particles::RemapAfterAMR()` remaps particles by
position onto the new leaf list before handing off-rank particles to the normal
particle MPI exchange.

The remap path is:

1. Copy local particle positions and integer IDs to host memory.
2. Wrap periodic particle coordinates into the mesh domain.
3. Convert each wrapped position to a finest-level logical location using
   `Mesh::max_level`.
4. Descend the mesh tree with `MeshBlockTree::FindLeafContaining()`.
5. Assign the returned leaf MeshBlock `gid` to `PGID`.
6. Build the particle send list for particles whose new `PGID` belongs to
   another MPI rank.
7. Run the existing particle MPI exchange and update particle counts.

This avoids the previous host-side scan over every MeshBlock for every particle.
When `check_consistency = true`, the tree result is checked against the returned
MeshBlock's physical bounds.

The particle remap is separate from face-centered MHD AMR repair.  It does not
modify `RefineFC()`, `RestrictFC()`, or `RepairAMRFC()`, and it is called only
after mesh metadata and physics arrays have been rebuilt.

## Consistency Checks

`<particles>/check_consistency = true` enables fail-fast checks after problem
generation, particle restart load, AMR remap, particle MPI exchange, and
particle restart output.

The checks validate:

- every particle `PGID` maps to a live MeshBlock on the owning rank,
- every particle position lies inside its assigned MeshBlock,
- total and per-species particle counts remain conserved,
- tags are unique within each species,
- particle real fields are finite,
- particle restart files contain exactly the declared header and data length.

These checks copy particle data to host and use MPI reductions where needed.
They are intended for regression tests and debugging, not production runs.

## Outputs

| File type | Output directory | Description |
|-----------|------------------|-------------|
| `ppd` | `ppd/` | Shared compact particle species and position dump. |
| `prst` | `prst/rank_XXXXXXXX/` | Per-rank particle restart dump. |
| `df` | `df/` | Pitch-angle cosine histogram. |
| `dxh` | `dxh/` | Displacement histograms for `IPDX`, `IPDY`, and `IPDZ`. |
| `pvtk` | run directory | Particle VTK output. |
| `trk` | run directory | Tracked-particle output. |

`prst` files store a three-`Real` header `(time, dt, particles_this_rank)`,
followed by 17 `Real` values per particle: `PGID`, `PTAG`, `PSP`, and the 14
real particle fields.

For `df` and `dxh`, `reduce = true` is the default.  MPI runs write one
globally reduced histogram from rank zero.  The sum of each `df` species
histogram equals that species particle count.  Each `dxh` species/component
histogram also sums to that species count.

Set `df_single_file_per_rank = 1` or `dxhist_single_file_per_rank = 1` to
preserve per-rank histogram output.  Those options force `reduce = false`.

## Particle Restart

Particle restart is independent from full mesh restart.  To reload particle
data:

```text
<problem>
prtcl_rst_flag = 1
prtcl_res_file = prst/rank_00000000/cr_tracer_amr.00002.prst
```

In MPI, each rank replaces `rank_00000000` with its own rank directory.  The
mesh configuration must remain compatible with the saved particle `PGID`
values.

## Readback Utility

Use `scripts/particles/cr_tracer_inspect.py` to inspect and validate particle
outputs:

```bash
python scripts/particles/cr_tracer_inspect.py run-dir \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8
```

The utility reads `prst` and `ppd`, reports totals by rank and species, verifies
restart header/data lengths, checks finite particle data, checks unique tags
within each species, and validates reduced `df`/`dxh` histogram sums when
`--histogram-nbin` is supplied.

The particle regression tests import this utility directly so command-line and
automated checks use the same parser.

## Regression Suite

Particle tests live in `tst/test_suite/particles/`.

```bash
cd tst
python run_test_suite.py --test test_suite/particles/test_particles_cr_cpu.py --cpu
python run_test_suite.py --test test_suite/particles/test_particles_cr_mpicpu.py --mpicpu
python run_test_suite.py --style
git diff --check
```

The CPU suite covers:

- serial drift smoke,
- serial particle restart write/read,
- Boris AMR with `tsc`,
- Boris AMR with `lin`,
- reduced `df`/`dxh` histogram sums,
- multi-species `ppd` and `prst` count/tag checks.

The MPI suite covers:

- two-rank Boris AMR particle migration,
- moving-box AMR refine/derefine stress,
- total and per-species conservation across rank exchange,
- MPI particle restart write/read,
- consistency-check validation during AMR remap and MPI exchange.

Local validation for this branch on 2026-05-22:

| Command | Result |
|---------|--------|
| `test_particles_cr_cpu.py --cpu` | 4 passed. |
| `test_particles_cr_mpicpu.py --mpicpu` | 3 passed. |
| `python run_test_suite.py --style` | 2 passed. |
| `git diff --check` | Passed. |
| `-D PROBLEM=part_random` build | Passed. |

The explicit two-rank stress run reported `105 MeshBlocks created, 77 deleted
by AMR`.  The inspector reported 1024 total particles, rank counts `0=520
1=504`, species counts `0=512 1=512`, and valid reduced `df`/`dxh` histograms.

Div(B) AMR compatibility inputs were rerun against the same built-in pgen
build:

| Input | AMR result |
|-------|------------|
| `inputs/tests/divb_amr_1d.athinput` | 69 MeshBlocks created, 50 deleted. |
| `inputs/tests/divb_amr_2d.athinput` | 1311 MeshBlocks created, 1095 deleted. |
| `inputs/tests/divb_amr_3d.athinput` | 3633 MeshBlocks created, 2842 deleted. |

## Integration Notes

The GitHub Pages branch already includes this page in the module toctree as
`docs/source/modules/particles.md`.  To publish these docs there, copy or merge
this file over the stale page on `origin/gh-pages`; no `index.md` toctree edit
is needed.
