# Cosmic-Ray Tracer Particles

This directory contains the AthenaK particle module and the cosmic-ray tracer
extension.  The implementation is built on the existing particle task list and
particle boundary exchange path, with additional fields, pushers, outputs, AMR
remapping support, reduced diagnostics, and test-oriented validation checks for
charged tracer particles.

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
check_motion_bounds = false

<problem>
pgen_name = part_random
particle_position = random
particle_velocity = random

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
particle_load_weight = 0.0

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
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
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
check_consistency_mode = none
check_motion_bounds = false
log_performance = false
amr_remap = device_table
validate_amr_lookup = false
amr_lookup_max_cells = 20000000
exchange_mode = alltoall_counts
subcycle = false
subcycle_max_steps = 8
subcycle_cell_fraction = 0.5
subcycle_meshblock_fraction = 0.5
subcycle_gyro_fraction = 0.25
subcycle_strict = true
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
- `interpolation = trilinear`: for `pusher = boris`, uses full
  multidimensional trilinear interpolation of the cell-centered magnetic field.
- `interpolation = lin` or `lin_legacy`: for `pusher = boris`, uses the legacy
  component-wise face-centered magnetic-field gather.
- `min_mass`: species-zero gyro parameter used by the Boris update.
- `mass_log_spacing`: multiplicative spacing between species masses.
- `cfl_part`: particle timestep cap used by `part_random`.
- `assign_tag = index_order` or `rank_order`: controls initial particle tags.
- `check_consistency = true`: legacy alias for full host-side invariant checks.
- `check_consistency_mode = none`, `counts`, `local`, or `full`: selects how
  much consistency checking to do.  Leave it at `none` for production runs.
- `check_motion_bounds = true`: fails after a push if a particle has moved
  outside the old one-neighbor MeshBlock stencil.
- `log_performance = true`: prints lightweight particle push, exchange, remap,
  and diagnostic counters.
- `amr_remap = device_table`: remaps particles after AMR with a flattened
  finest-level logical lookup table.  Use `host_tree` as a debug fallback.
- `validate_amr_lookup = true`: compares lookup-table remap results against the
  host mesh-tree lookup.
- `amr_lookup_max_cells`: maximum flattened AMR lookup entries before the code
  falls back to `host_tree`.
- `exchange_mode = alltoall_counts`: exchanges MPI rank counts without
  all-gathering detailed send metadata.  Use `allgather` as a compatibility
  fallback.
- `subcycle = true`: enables opt-in particle subcycling within each mesh
  timestep.
- `subcycle_max_steps`: maximum allowed substeps.  With
  `subcycle_strict = true`, the run fails if the constraints need more.
- `subcycle_cell_fraction`: per-substep cell-crossing fraction limit.
- `subcycle_meshblock_fraction`: per-substep MeshBlock-crossing fraction limit.
- `subcycle_gyro_fraction`: per-substep Boris gyro-angle limit.
- `subcycle_strict`: if false, caps subcycling at `subcycle_max_steps` instead
  of failing immediately.

`inputs/particles/cr_tracer_boris_amr.athinput` is the branch smoke-test input
for a 3D MHD, Boris-pushed, AMR cosmic-ray tracer run.
`inputs/particles/cr_tracer_boris_amr_stress.athinput` is the moving-box AMR
refine/derefine stress test used by MPI validation and output checks.
`inputs/particles/cr_tracer_boris_uniform.athinput` is a deterministic
one-particle uniform-field Boris accuracy test.
`inputs/particles/cr_tracer_boris_amr_perf.athinput` and
`inputs/particles/cr_tracer_drift_amr_perf.athinput` are CPU/MPI performance
comparison inputs for end-to-end Boris AMR and remap-focused AMR timing.

Useful `part_random` problem settings:

- `particle_position = random`, `center`, `meshblock_center`, `fixed`, or
  `tag_random`.  The `fixed` mode uses `particle_x1`, `particle_x2`, and
  `particle_x3`.  `tag_random` generates deterministic global positions from
  species and tag.
- `particle_velocity = random`, `uniform`, `isotropic`, or
  `isotropic_tag_random`.  The `uniform` mode uses `v0x`, `v0y`, and `v0z`.
  `isotropic_tag_random` generates deterministic directions from species and
  tag while preserving `v0`.
- `particle_seed` selects the deterministic `tag_random` sequence.  These
  modes are intended for MPI-decomposition and resolution comparisons; initial
  ownership is assigned through the normal remap/exchange path.
- `B_profile = uniform`, `linear_cross`, `sinusoidal_divb_free`, `mirror`,
  `gradb`, or `turbulent`.  These prescribed accuracy fields use `Bgrad`, or
  `Bamp` and `Bwave_number`, as applicable.

## Particle Subcycling

Subcycling is disabled by default.  Enable it when particle velocities or Boris
gyro motion would otherwise make the field gather stale or require a smaller
gyro-angle step:

```text
<particles>
subcycle = true
subcycle_max_steps = 64
subcycle_cell_fraction = 0.5
subcycle_meshblock_fraction = 0.5
subcycle_gyro_fraction = 0.25
subcycle_strict = true
```

The particle module computes the required substep count from cell crossing,
MeshBlock crossing, and gyro-angle constraints, then takes the maximum across
MPI ranks.  Between intermediate substeps it updates `PGID` and runs the same
particle MPI exchange path used by the normal task list, so the next field
gather uses the current owner.  The final substep is followed by the normal
particle `NewGID` and exchange tasks.

`subcycle_strict = true` is recommended for tests and production runs because
it fails with the required cell, MeshBlock, and gyro substep counts when the
configured maximum is too low.  `subcycle_strict = false` is useful for
experiments, but should be paired with `check_motion_bounds = true`.

## Field Gather Checks

The CPU regression suite includes a smooth-field check using
`B_profile = linear_cross`.  In that test, `trilinear` recovers the analytic
point-centered field to roundoff for a linear cross-field, while `lin_legacy`
shows the expected cross-direction interpolation error that decreases with
resolution.  This is a pusher/gather test, not an MHD evolution test.

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

The default remap procedure is:

1. Build a flattened finest-level logical lookup table from the live AMR leaf
   list.
2. In a Kokkos kernel, wrap periodic coordinates back into the mesh domain.
3. Convert the wrapped position to a finest-level logical table index using
   `Mesh::max_level`.
4. Update `PGID` to the table's leaf MeshBlock ID.
5. Build the particle send list for particles whose new block is owned by
   another rank.
6. Exchange particle MPI metadata with `MPI_Alltoall` rank counts.
7. Pack, send, receive, unpack, and compact migrated particles.
8. Recompute per-rank and total particle counts with
   `Mesh::UpdateParticleCounts()`.

This replaces the previous host-side scan over every leaf MeshBlock for every
particle.  Set `amr_remap = host_tree` to use the host mesh-tree fallback.  Set
`exchange_mode = allgather` to use the older detailed metadata exchange path.
When `validate_amr_lookup = true`, or when local/full consistency checks are
active, the lookup-table result is validated against
`MeshBlockTree::FindLeafContaining()` and the returned MeshBlock's physical
bounds.

The post-push particle `NewGID` task also uses position-based lookup on
multilevel meshes, so AMR ownership after particle motion is checked against
the live leaf list.  Uniform meshes keep the cheaper neighbor-stencil update.

The AMR particle remap is intentionally separate from the face-centered MHD AMR
repair used by the div(B) fix.  It does not modify `RefineFC()`,
`RestrictFC()`, or `RepairAMRFC()`, and the particle remap is called only after
the mesh AMR metadata and physics arrays have been rebuilt.

## Particle-Aware Load Balancing

AMR uses mesh-only load balancing by default.  Set
`<mesh_refinement>/particle_load_weight` to a positive value to add CR tracer
particles to the post-AMR MeshBlock cost estimate:

```text
<mesh_refinement>
particle_load_weight = 0.001
```

The load balancer uses the pre-AMR particle owners to estimate the new costs.
If an old MeshBlock is refined, its particles are divided evenly among the new
children for costing only; the later particle AMR remap still assigns every
particle from its actual position.  If MeshBlocks are derefined, the old child
particle counts are assigned to the new parent.  `particle_load_weight = 0.0`
preserves the original mesh-only behavior and is the right baseline for
performance comparisons.

## Outputs

The inherited particle outputs are still available:

- `pvtk`: particle VTK output.
- `trk`: tracked-particle output.

The cosmic-ray tracer branch adds particle-oriented output types:

- `ppd`: compact species and position dump.  Each particle contributes
  `(species, x, y, z)` as floats.
- `prst`: per-rank particle restart dump.  The header stores
  `(time, dt, particles_this_rank)` as `Real`, followed by 17 `Real` values per
  particle containing `PGID`, `PTAG`, `PSP`, and the 14 real particle fields.
- `df`: pitch-angle cosine histogram.  Use `nbin`, `vmin`, and `vmax` to set the
  binning range.
- `pspec`: per-species spectrum.  Set `quantity` to `p`, `E`, `logE`, `mu`, or
  `magnetic_moment`; use `nbin`, `vmin`, and `vmax` to set the binning range.
- `pspec2`: per-species two-dimensional spectrum.  Set `quantity` to `mu_p`,
  `mu_E`, or `vpar_vperp`; use `nbin1`, `nbin2`, `vmin1`, `vmax1`, `vmin2`,
  and `vmax2` to set the binning range.
- `dxh`: displacement histograms for `IPDX`, `IPDY`, and `IPDZ`.  Use `nbin`,
  `vmin`, and `vmax` to set the binning range.
- `drh`: scalar displacement histogram for
  `sqrt(IPDX^2 + IPDY^2 + IPDZ^2)`.
- `dparh`: parallel-displacement histogram using `IPDB`.
- `pmom`: reduced per-species moments including count, `<mu>`, `<mu^2>`,
  `3<mu^2>-1`, displacement moments, parallel-displacement moments, `<v^2>`,
  flux/current proxies `<v_i>`, pressure-tensor proxies `<v_i v_j>`, and
  diffusion-tensor proxies `<dx_i dx_j>`.
- `psamp`: per-rank text sample of selected particle fields.  The sample is
  deterministic in species and tag:
  `splitmix64(species, tag) % sample_stride == sample_offset`.

For `df`, `pspec`, `dxh`, `drh`, `dparh`, and `pmom`, `reduce = true` is the
default.  MPI runs write one globally reduced diagnostic from rank zero, so
each histogram species sum equals the total particle count for that species.
Set `df_single_file_per_rank = 1`, `pspec_single_file_per_rank = 1`,
`dxhist_single_file_per_rank = 1`, `drh_single_file_per_rank = 1`,
`dparh_single_file_per_rank = 1`, or `pmom_single_file_per_rank = 1` to
preserve per-rank output.  Those settings force `reduce = false`.

The `pspec` quantities use stored particle velocities and the sampled pusher
magnetic field.  `p` is the speed proxy, `E` is `0.5 p^2`, `logE` is
`log10(E)`, `mu` is the pitch-angle cosine, and `magnetic_moment` is the
diagnostic proxy `v_perp^2 / B`.  The output does not resample the current MHD
field during the diagnostic write.

`pspec2` provides the reduced joint distributions `f(mu,p)`, `f(mu,E)`, and
`f(v_parallel,v_perp)`.  It uses the same stored particle velocities and
sampled pusher fields as `pspec`.

`psamp` is for selected-field analysis and debugging, not restart.  It always
writes `rank`, `pgid`, `tag`, and `species`, then the user-selected `fields`.
Stored fields are `x`, `y`, `z`, `vx`, `vy`, `vz`, `mass`, `bx`, `by`, `bz`,
`dx`, `dy`, `dz`, and `dpar`.  Derived fields are `speed`, `energy`, `bmag`,
`mu`, `vpar`, `vperp`, and `magnetic_moment`.  Use `sample_stride = 1` for all
particles, or a larger stride for a stable subset.  Use `species = -1` for all
species.

## Consistency Checks

`<particles>/check_consistency = true` is a legacy alias for
`<particles>/check_consistency_mode = full`.  The available modes are:

- `none`: no extra host-side checks,
- `counts`: conserved total and per-species counts,
- `local`: count checks plus local ownership, bounds, finite-value, and
  duplicate-tag checks,
- `full`: local checks plus global duplicate-tag checks across MPI ranks.

The full mode fails fast if any of these invariants are violated:

- every particle `PGID` maps to a live MeshBlock on the owning rank,
- every particle position lies inside its assigned MeshBlock,
- total and per-species particle counts remain conserved after initialization,
  AMR remap, MPI exchange, restart load, and restart output,
- tags are unique within each species,
- particle restart files have exactly the header and data length implied by the
  stored particle count.

The checks copy particle data to host and use MPI reductions where needed.  They
are intended for regression tests and debugging, not for production runs.

`<particles>/check_motion_bounds = true` is a cheaper exchange guardrail.  It
detects a pushed particle that has moved outside the old one-neighbor MeshBlock
stencil.  It is a diagnostic check; use `subcycle = true` when large particle
velocities need intermediate ownership updates or smaller Boris gyro-angle
steps within one mesh timestep.

## Readback Utility

`scripts/particles/cr_tracer_inspect.py` reads `prst`, `ppd`, `psamp`, `df`,
`pspec`, `pspec2`, `dxh`, `drh`, `dparh`, and `pmom` outputs.  It reports totals by rank
and species, reports particles per MeshBlock `PGID`, reports restart format
metadata, validates optional `*.prst.pmeta` sidecar payload sizes and checksums,
verifies particle restart file lengths, checks finite positions and velocities,
checks unique tags within each species, validates reduced histogram sums,
reports sampled-field totals and field names, and can assert expected
total/species counts.

Typical commands:

```bash
python scripts/particles/cr_tracer_inspect.py run-dir/prst
python scripts/particles/cr_tracer_inspect.py run-dir \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
```

The regression tests import the same utility directly, so command-line and test
validation use the same parser and count checks.

## Performance Logging

Set `<particles>/log_performance = true` to print particle counters from the
push, MPI exchange, AMR remap, and diagnostic-output paths.  For example:

```bash
./athena -i inputs/particles/random_particle_drift.athinput \
  particles/log_performance=true \
  particles/ppc=0.001953125 \
  problem/particle_position=center \
  problem/particle_velocity=uniform \
  time/nlim=1
```

The local 2026-05-22 smoke run reported:

```text
Particle performance push: pushed=1 sent=0 received=0 remapped=0 diagnostics=0
```

Use these counters to compare runs on the same machine; do not treat them as
portable performance guarantees.

## CPU/MPI Performance Comparison

The full CPU/MPI benchmark write-up is in
`docs/source/modules/cr_tracer_cpu_mpi_performance.md`.  The 2026-05-22 local
comparison used a Release MPI build with the Kokkos Serial backend on an Apple
M4 Max.  It includes the original old-scan versus tree-remap comparison and a
four-mode follow-up matrix for `device_table`/`host_tree` and
`alltoall_counts`/`allgather`.  GPU performance was not tested on this machine
and remains a TODO.

Median results from three runs:

| Benchmark | MPI ranks | Old scan/full-copy CPU s | Current CPU s | Speedup |
|-----------|-----------|--------------------------|---------------|---------|
| Remap-focused drift AMR | 1 | 0.5686 | 0.00386 | 147.4x |
| Remap-focused drift AMR | 2 | 0.4210 | 0.00540 | 77.9x |
| Remap-focused drift AMR | 4 | 0.2786 | 0.01233 | 22.6x |
| Boris+MHD AMR end-to-end | 1 | 1.1082 | 0.7974 | 1.39x |
| Boris+MHD AMR end-to-end | 2 | 0.7437 | 0.5433 | 1.37x |
| Boris+MHD AMR end-to-end | 4 | 0.4789 | 0.3451 | 1.39x |

The mesh-tree lookup is the dominant performance fix.  The later position-only
host-copy optimization is neutral within noise for the small CPU Boris case
after the tree lookup is present, but it removes unnecessary AMR remap data
movement and is the right direction for GPU follow-up.

For the follow-up architecture slice on CPU/MPI, `device_table` is
1.13x-1.63x faster than `host_tree` in the remap-focused drift matrix and
neutral to 1.02x faster in the Boris+MHD matrix.  `alltoall_counts` is neutral
within local noise in the drift, exchange, and Boris+MHD matrices.

Use the benchmark harness for repeatable local sweeps:

```bash
python scripts/particles/cr_tracer_benchmark.py \
  --athena build-cr-robust-mpi/src/athena \
  --cases push,amr_remap,exchange,histogram,histogram_per_rank,restart,boris_amr \
  --modes default,host_tree,allgather,fallback \
  --ranks 1,2,4 \
  --repeats 3 \
  --mpi-n1 \
  --output-json cr_tracer_bench_cpu_mpi.json \
  --output-md cr_tracer_bench_cpu_mpi.md
```

Reduced `df`/`pspec`/`pspec2`/`dxh`/`drh`/`dparh`/`pmom` diagnostic output had no
significant runtime penalty in the local four-rank MPI test.  It reduced each
diagnostic from four per-rank files to one global file, which is why it remains
the default.

## Local Validation

The regression suite exercises serial drift, serial Boris AMR with `tsc` and
`lin` interpolation, restart write/read, uniform-field Boris speed
conservation, periodic wrap displacement accounting, the large-motion guard,
reduced `df`/`pspec`/`dxh`/`drh`/`dparh` histograms, `pmom` moments,
reduced `pspec2` joint spectra, deterministic `psamp` readback,
smooth-field `trilinear` convergence,
multi-species tags/counts, two-rank MPI AMR migration, and a four-rank
moving-box AMR refine/derefine stress run.

Required commands:

```bash
cd tst
python run_test_suite.py --test test_suite/particles/test_particles_cr_cpu.py --cpu
python run_test_suite.py --test test_suite/particles/test_particles_cr_mpicpu.py --mpicpu
python run_test_suite.py --style
git diff --check
```

Local results from 2026-05-22:

- `test_particles_cr_cpu.py --cpu`: 20 passed.
- `test_particles_cr_mpicpu.py --mpicpu`: 4 passed.
- `--style`: 2 passed.
- `git diff --check`: passed.
- Custom `-D PROBLEM=part_random` build: passed.

The explicit four-rank stress run:

```bash
mpirun -np 4 build-cr-robust-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_stress.athinput \
  -d run-cr-stress-mpi-final \
  particles/check_consistency_mode=full \
  particles/validate_amr_lookup=true \
  mesh_refinement/particle_load_weight=0.001
python scripts/particles/cr_tracer_inspect.py run-cr-stress-mpi-final \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
```

Stress result: `105 MeshBlocks created, 77 deleted by AMR`; inspector result:
1024 total particles, rank counts `0=238 1=259 2=267 3=260`, species counts
`0=512 1=512`, valid `AKPRST-v1` restart metadata, and valid reduced
`df`/`pspec`/`pspec2`/`dxh`/`drh`/`dparh` histogram totals.

The div(B) AMR compatibility inputs also passed:

```bash
./athena -i inputs/tests/divb_amr_1d.athinput
./athena -i inputs/tests/divb_amr_2d.athinput
./athena -i inputs/tests/divb_amr_3d.athinput
```

Run the deeper normalized-divergence regression from the repository test
directory:

```bash
cd tst
python run_tests.py mhd/mhd_divb_amr
```

Div(B) results:

- 1D: `cycle=24`, `Current number of MeshBlocks = 27`,
  `69 MeshBlocks created, 50 deleted by AMR`.
- 2D: `cycle=20`, `Current number of MeshBlocks = 252`,
  `1311 MeshBlocks created, 1095 deleted by AMR`.
- 3D: `cycle=12`, `Current number of MeshBlocks = 855`,
  `3633 MeshBlocks created, 2842 deleted by AMR`.

The normalized-divergence regression passed its 1D, 2D, and 3D
physical-level-1-through-5 cases.  The deepest 3D case completed with `30,829`
live MeshBlocks after `50,673` creations and `19,908` deletions while staying
within the scripted normalized `div(B)` tolerances.

## Known Validation Limits

The validation above is CPU and local MPI focused.  It does not replace a large
multi-node run, and no GPU build is required for this hardening pass.
