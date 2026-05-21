# Cosmic-Ray Tracer Particles

This directory contains the AthenaK particle module and the cosmic-ray tracer
extension.  The implementation is built on the existing particle task list and
particle boundary exchange path, with additional fields, pushers, outputs, and
AMR remapping support for charged tracer particles.

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
```

The supported particle settings are:

- `particle_type = cosmic_ray`: selects the cosmic-ray tracer data layout.
- `ppc`: number of particles per cell per species.  Values below one are
  allowed.
- `nspecies`: number of particle species.  Tags are unique across species.
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

`inputs/particles/cr_tracer_boris_amr.athinput` is the branch smoke-test input
for a 3D MHD, Boris-pushed, AMR cosmic-ray tracer run.

## Particle Data

Cosmic-ray tracer particles allocate 14 real fields and 3 integer fields.

Integer fields:

- `PGID`: global MeshBlock ID containing the particle.
- `PTAG`: unique particle tag.
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
3. Find the new leaf MeshBlock whose logical location contains the particle.
4. Update `PGID` to the new MeshBlock ID.
5. Use the existing particle MPI boundary path for particles whose new block is
   owned by another rank.
6. Recompute per-rank and total particle counts with
   `Mesh::UpdateParticleCounts()`.

This AMR particle remap is intentionally separate from the face-centered MHD AMR
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
  binning range.  `df_single_file_per_rank = 1` writes separate rank files.
- `dxh`: displacement histograms for `IPDX`, `IPDY`, and `IPDZ`.  Use `nbin`,
  `vmin`, and `vmax` to set the binning range.

For `df` and `dxh`, the default shared file stores each rank's local histogram
at a rank-specific offset.  It is not reduced across ranks before writing.

Particle restart reads are enabled with:

```text
<problem>
prtcl_rst_flag = 1
prtcl_res_file = prst/rank_00000000/basename.00042.prst
```

For MPI restarts, the input may name the rank-zero file.  The reader replaces
`rank_00000000` with each rank's local restart directory name.

## Local Validation

These results were produced on 2026-05-20 from branch `feature/CR_tracers`
based on `main` plus merge `c3cfc63f`, using CMake 4.3.2 and Open MPI 5.0.9.
The local validation used CPU builds.

Build commands:

```bash
cmake -S . -B build-cr-tracers -D PROBLEM=part_random
cmake --build build-cr-tracers --target athena -j 6

cmake -S . -B build-cr-tracers-mpi -D PROBLEM=part_random -DAthena_ENABLE_MPI=ON
cmake --build build-cr-tracers-mpi --target athena -j 6

cmake -S . -B build-cr-tracers-builtin -D PROBLEM=built_in_pgens
cmake --build build-cr-tracers-builtin --target athena -j 6
```

Serial drift smoke test:

```bash
build-cr-tracers/src/athena \
  -i inputs/particles/random_particle_drift.athinput \
  -d run-cr-tracers-smoke \
  mesh/nx1=16 mesh/nx2=16 mesh/nx3=16 \
  meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=8 \
  particles/ppc=0.125 time/nlim=2 time/tlim=0.02 \
  output1/dt=0.02 output2/dt=0.02
```

Result: completed at `time=2.000000e-02`, `cycle=1`,
`MeshBlock-cycles = 8`, and
`particle-updates/cpu_second = 2.299832e+06`.

Serial Boris plus AMR smoke test:

```bash
build-cr-tracers/src/athena \
  -i inputs/particles/cr_tracer_boris_amr.athinput \
  -d run-cr-tracers-boris-amr
```

Result: completed at `time=3.000000e-02`, `cycle=3`; AMR ended with
`Current number of MeshBlocks = 64` and
`56 MeshBlocks created, 0 deleted by AMR`.

MPI Boris plus AMR test:

```bash
mpirun -np 2 build-cr-tracers-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr.athinput \
  -d run-cr-tracers-mpi-amr
```

Result: completed at `time=3.000000e-02`, `cycle=3`; AMR ended with
`Current number of MeshBlocks = 64`,
`56 MeshBlocks created, 0 deleted by AMR`, and
`load balancing efficiency = 1.000000e+00`.

The MPI restart headers show that particles migrated across ranks while the
total count stayed fixed at 1024:

```text
rank 0: 512 -> 515 -> 516 particles
rank 1: 512 -> 509 -> 508 particles
```

Div(B) AMR regression inputs were also rerun with the built-in pgen build:

```bash
build-cr-tracers-builtin/src/athena \
  -i inputs/tests/divb_amr_1d.athinput \
  -d run-divb-amr-1d
build-cr-tracers-builtin/src/athena \
  -i inputs/tests/divb_amr_2d.athinput \
  -d run-divb-amr-2d
build-cr-tracers-builtin/src/athena \
  -i inputs/tests/divb_amr_3d.athinput \
  -d run-divb-amr-3d
```

Results:

- 1D: `cycle=24`, `Current number of MeshBlocks = 27`,
  `69 MeshBlocks created, 50 deleted by AMR`.
- 2D: `cycle=20`, `Current number of MeshBlocks = 252`,
  `1311 MeshBlocks created, 1095 deleted by AMR`.
- 3D: `cycle=12`, `Current number of MeshBlocks = 855`,
  `3633 MeshBlocks created, 2842 deleted by AMR`.

Style and consistency checks:

```bash
git diff --check
awk 'length($0)>100 {print FILENAME ":" FNR ":" length($0)}' <changed files>
```

`git diff --check` passed.  The changed C++ and input files had no lines longer
than 100 characters at the time this documentation was written.  `cpplint` was
not available in the local Python environment, so no `cpplint` run was included.

## Known Validation Limits

The local validation above confirms serial operation, local two-rank MPI AMR
particle migration, and compatibility with the current div(B) AMR regression
inputs.  It does not replace a large multi-node run, and no GPU build was tested
in this branch validation pass.
