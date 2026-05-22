# Cosmic-Ray Tracer Particles

This directory contains the AthenaK particle module and the cosmic-ray tracer
extension.  The implementation is built on the existing particle task list and
particle boundary exchange path, with additional fields, pushers, outputs, AMR
remapping support, and test-oriented validation checks for charged tracer
particles.

## Quick Start

Use `pgen_name = part_random` for the standard built-in particle generator.  A
custom `-D PROBLEM=part_random` build is still supported, but normal tests and
examples use the built-in pgen path.

Minimal drifting CR tracer input:

```text
<particles>
particle_type     = cosmic_ray
ppc               = 0.125
nspecies          = 2
pusher            = drift
check_consistency = true

<problem>
pgen_name = part_random

<output1>
file_type = prst
dt        = 0.01

<output2>
file_type = ppd
dt        = 0.01
```

Minimal Boris/MHD AMR setup:

```text
<mesh_refinement>
refinement = adaptive
num_levels = 2

<amr_criterion0>
method = location

<mhd>
eos     = ideal
rsolver = hlld

<particles>
particle_type     = cosmic_ray
ppc               = 0.125
nspecies          = 2
pusher            = boris
interpolation     = tsc
check_consistency = true

<problem>
pgen_name = part_random
B0z       = 1.0
```

The moving refine/derefine AMR stress example is:

```bash
./athena -i inputs/particles/cr_tracer_boris_amr_stress.athinput
```

For MPI:

```bash
mpirun -np 2 ./athena -i inputs/particles/cr_tracer_boris_amr_stress.athinput
```

Validate the latest particle outputs from a run directory with:

```bash
python scripts/particles/cr_tracer_inspect.py . \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8
```

Particle restart reads are enabled with:

```text
<problem>
prtcl_rst_flag = 1
prtcl_res_file = prst/rank_00000000/basename.00042.prst
```

For MPI particle restarts, the input may name the rank-zero file.  The reader
replaces `rank_00000000` with each rank's local restart directory name.  The
mesh configuration must be compatible with the saved particle `PGID` values; a
full mesh restart is separate from this particle-only `prst` path.

## Input Blocks

Enable particles with a `<particles>` block.  The cosmic-ray tracer mode uses

```text
<particles>
particle_type     = cosmic_ray
ppc               = 0.125
nspecies          = 2
pusher            = boris
interpolation     = tsc
min_mass          = 1.0
mass_log_spacing  = 2.0
cfl_part          = 0.1
check_consistency = false
```

The supported particle settings are:

- `particle_type = cosmic_ray`: selects the cosmic-ray tracer data layout.
- `ppc`: number of particles per cell per species.  Values below one are
  allowed.
- `nspecies`: number of particle species.
- `pusher = drift`: moves particles with stored velocities.
- `pusher = boris`: uses the magnetic Boris update and requires an `<mhd>`
  block.
- `interpolation = tsc`: for `pusher = boris`, interpolates the cell-centered
  magnetic field with triangular-shaped-cloud weights.
- `interpolation = lin`: for `pusher = boris`, linearly interpolates the
  face-centered magnetic field components.
- `min_mass`: species-zero gyro parameter used by the Boris update.
- `mass_log_spacing`: multiplicative spacing between species masses.
- `cfl_part`: particle timestep cap used by `part_random`.
- `assign_tag = index_order` or `rank_order`: controls initial particle tags.
- `check_consistency = true`: enables host-side invariant checks intended for
  tests and debugging.  Leave it off for production runs.

`inputs/particles/cr_tracer_boris_amr.athinput` is the branch smoke-test input
for a 3D MHD, Boris-pushed, AMR cosmic-ray tracer run.

## Particle Data

Cosmic-ray tracer particles allocate 14 real fields and 3 integer fields.

Integer fields:

- `PGID`: global MeshBlock ID containing the particle.
- `PTAG`: particle tag, unique within each species.
- `PSP`: species index.

Real fields:

- `IPX`, `IPY`, `IPZ`: particle position.
- `IPVX`, `IPVY`, `IPVZ`: particle velocity.
- `IPM`: particle gyro parameter.
- `IPBX`, `IPBY`, `IPBZ`: magnetic field sampled by the pusher.
- `IPDX`, `IPDY`, `IPDZ`: accumulated displacement.
- `IPDB`: accumulated displacement parallel to the sampled magnetic field.

## AMR Behavior

AMR changes the global leaf MeshBlock list and can also move MeshBlocks between
MPI ranks.  After each AMR update, `Particles::RemapAfterAMR()` remaps every
particle by position onto the new leaf MeshBlock list.

The remap procedure is:

1. Copy particle positions and IDs to host memory.
2. Wrap periodic coordinates back into the mesh domain.
3. Convert the wrapped position to a finest-level logical location using
   `Mesh::max_level`.
4. Use `MeshBlockTree::FindLeafContaining()` to descend directly to the owning
   leaf MeshBlock.
5. Update `PGID` to the new MeshBlock ID.
6. Use the existing particle MPI boundary path for particles whose new block is
   owned by another rank.
7. Recompute per-rank and total particle counts with
   `Mesh::UpdateParticleCounts()`.

This replaces the previous host-side scan over every leaf MeshBlock for every
particle.  When `check_consistency = true`, the tree result is validated against
the returned MeshBlock's physical bounds.

The AMR particle remap is intentionally separate from the face-centered MHD AMR
repair used by the div(B) fix.  It does not modify `RefineFC()`,
`RestrictFC()`, or `RepairAMRFC()`, and the particle remap is called only after
the mesh AMR metadata and physics arrays have been rebuilt.

## Outputs

The inherited particle outputs are still available:

- `pvtk`: particle VTK output.
- `trk`: tracked-particle output.

The cosmic-ray tracer branch adds four particle-oriented output types:

- `ppd`: compact species and position dump.  Each particle contributes
  `(species, x, y, z)` as floats.
- `prst`: per-rank particle restart dump.  The header stores
  `(time, dt, particles_this_rank)` as `Real`, followed by 17 `Real` values per
  particle containing `PGID`, `PTAG`, `PSP`, and the 14 real particle fields.
- `df`: pitch-angle cosine histogram.  Use `nbin`, `vmin`, and `vmax` to set the
  binning range.
- `dxh`: displacement histograms for `IPDX`, `IPDY`, and `IPDZ`.  Use `nbin`,
  `vmin`, and `vmax` to set the binning range.

For `df` and `dxh`, `reduce = true` is the default.  MPI runs write one
globally reduced histogram from rank zero, so each species sum in `df` equals
the total particle count for that species, and each displacement component in
`dxh` sums to the total particle count for that species.  Set
`df_single_file_per_rank = 1` or `dxhist_single_file_per_rank = 1` to preserve
the older per-rank output mode; those settings force `reduce = false`.

## Consistency Checks

`<particles>/check_consistency = true` fails fast if any of these invariants are
violated:

- every particle `PGID` maps to a live MeshBlock on the owning rank,
- every particle position lies inside its assigned MeshBlock,
- total and per-species particle counts remain conserved after initialization,
  AMR remap, MPI exchange, restart load, and restart output,
- tags are unique within each species,
- particle restart files have exactly the header and data length implied by the
  stored particle count.

The checks copy particle data to host and use MPI reductions where needed.  They
are intended for regression tests and debugging, not for production runs.

## Readback Utility

`scripts/particles/cr_tracer_inspect.py` reads `prst`, `ppd`, `df`, and `dxh`
outputs.  It reports totals by rank and species, verifies particle restart file
lengths, checks finite positions and velocities, checks unique tags within each
species, and can assert expected total/species counts.

Typical commands:

```bash
python scripts/particles/cr_tracer_inspect.py run-dir/prst
python scripts/particles/cr_tracer_inspect.py run-dir \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8
```

The regression tests import the same utility directly, so command-line and test
validation use the same parser and count checks.

## Local Validation

The regression suite exercises serial drift, serial Boris AMR with `tsc` and
`lin` interpolation, restart write/read, reduced `df`/`dxh` histograms,
multi-species tags/counts, two-rank MPI AMR migration, and a moving-box AMR
refine/derefine stress run.

Required commands:

```bash
cd tst
python run_test_suite.py --test test_suite/particles/test_particles_cr_cpu.py --cpu
python run_test_suite.py --test test_suite/particles/test_particles_cr_mpicpu.py --mpicpu
python run_test_suite.py --style
git diff --check
```

Local results from 2026-05-22:

- `test_particles_cr_cpu.py --cpu`: 4 passed.
- `test_particles_cr_mpicpu.py --mpicpu`: 3 passed.
- `--style`: 2 passed.
- `git diff --check`: passed.
- Custom `-D PROBLEM=part_random` build: passed.

The explicit two-rank stress run:

```bash
mpirun -np 2 build-cr-robust-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_stress.athinput \
  -d run-cr-stress-mpi-final
python scripts/particles/cr_tracer_inspect.py run-cr-stress-mpi-final \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8
```

Stress result: `105 MeshBlocks created, 77 deleted by AMR`; inspector result:
1024 total particles, rank counts `0=520 1=504`, species counts `0=512 1=512`,
and valid reduced `df`/`dxh` histogram totals.

The div(B) AMR compatibility inputs also passed:

```bash
./athena -i inputs/tests/divb_amr_1d.athinput
./athena -i inputs/tests/divb_amr_2d.athinput
./athena -i inputs/tests/divb_amr_3d.athinput
```

Div(B) results:

- 1D: `cycle=24`, `Current number of MeshBlocks = 27`,
  `69 MeshBlocks created, 50 deleted by AMR`.
- 2D: `cycle=20`, `Current number of MeshBlocks = 252`,
  `1311 MeshBlocks created, 1095 deleted by AMR`.
- 3D: `cycle=12`, `Current number of MeshBlocks = 855`,
  `3633 MeshBlocks created, 2842 deleted by AMR`.

## Known Validation Limits

The validation above is CPU and local MPI focused.  It does not replace a large
multi-node run, and no GPU build is required for this hardening pass.
