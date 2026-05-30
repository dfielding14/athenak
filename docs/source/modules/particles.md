# Module: Particles

## Role in AthenaK

The Particles module implements Lagrangian cosmic-ray tracer particles on top of
the existing AthenaK task list, boundary exchange, and AMR infrastructure.  The
current supported particle type is `cosmic_ray`, with drift and magnetic Boris
pushers, per-particle species tags, restart output, compact position output,
pitch-angle histograms, spectra, displacement histograms, reduced particle moments, and
AMR-aware particle remap.

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
| `src/outputs/*_prtcl.cpp` | Particle `ppd`, `prst`, `df`, `pspec`, `dxh`, `drh`, `dparh`, `pmom`, VTK, and tracking outputs. |
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
| `inputs/particles/cr_tracer_boris_uniform.athinput` | One-particle uniform-field Boris accuracy test. |
| `inputs/particles/cr_tracer_boris_amr_perf.athinput` | End-to-end CPU/MPI Boris AMR performance comparison. |
| `inputs/particles/cr_tracer_drift_amr_perf.athinput` | Remap-focused CPU/MPI AMR performance comparison. |

## Runtime Parameters

### `<particles>`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `particle_type` | required | Must be `cosmic_ray`. |
| `ppc` | `1.0` | Particles per cell per species.  Fractional values are allowed. |
| `nspecies` | `1` | Number of particle species. |
| `pusher` | required | `drift` or `boris`. |
| `interpolation` | `tsc` | For `boris`: `tsc`, `trilinear`, `lin`, or `lin_legacy`. |
| `min_mass` | `1.0` | Species-zero gyro parameter used by the Boris update. |
| `mass_log_spacing` | `1.0` | Multiplicative spacing between species gyro parameters. |
| `cfl_part` | `0.05` | Particle timestep cap used by `part_random`. |
| `assign_tag` | `index_order` | Initial tag assignment; `rank_order` is also supported. |
| `check_consistency` | `false` | Legacy boolean alias for `check_consistency_mode = full`. |
| `check_consistency_mode` | `none` | `none`, `counts`, `local`, or `full`. |
| `check_motion_bounds` | `false` | Fails if a pushed particle leaves the old one-neighbor update stencil. |
| `log_performance` | `false` | Prints lightweight particle push, remap, exchange, and diagnostic counters. |
| `amr_remap` | `device_table` | `device_table` uses a flattened finest-level logical lookup table; `host_tree` is the debug fallback. |
| `validate_amr_lookup` | `false` | Compares device-table remap results to the host mesh-tree lookup. |
| `amr_lookup_max_cells` | `20000000` | Maximum flattened AMR lookup entries before falling back to `host_tree`. |
| `exchange_mode` | `alltoall_counts` | `alltoall_counts` exchanges rank counts only; `allgather` preserves the older detailed metadata path. |
| `subcycle` | `false` | Enables opt-in particle subcycling during each particle push. |
| `subcycle_max_steps` | `8` | Maximum allowed substeps before strict mode fails or non-strict mode caps. |
| `subcycle_cell_fraction` | `0.5` | Limits per-substep motion to this fraction of the local cell width. |
| `subcycle_meshblock_fraction` | `0.5` | Limits per-substep motion to this fraction of the local MeshBlock size. |
| `subcycle_gyro_fraction` | `0.25` | Limits Boris gyro angle per substep in radians. |
| `subcycle_strict` | `true` | If true, fail when the requested constraints need more than `subcycle_max_steps`. |

### `<problem>` for `part_random`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pgen_name` | required for built-in builds | Use `part_random` for the standard pgen. |
| `prtcl_rst_flag` | `0` | If nonzero, load particle data from `prtcl_res_file`. |
| `prtcl_res_file` | required for particle restart | Per-rank `prst` path; MPI ranks replace `rank_00000000` with their rank directory. |
| `particle_refinement` | `none` | `moving_boxes` enrolls the AMR stress-test refinement pattern. |
| `particle_position` | `random` | `random`, `center`, `meshblock_center`, `fixed`, or `tag_random`; `fixed` uses `particle_x1`, `particle_x2`, and `particle_x3`, while `tag_random` is deterministic across MPI decompositions. |
| `particle_x1`, `particle_x2`, `particle_x3` | `0, 0, 0` | Fixed particle position used when `particle_position = fixed`. |
| `particle_velocity` | `random` | `random`, `uniform`, `isotropic`, or `isotropic_tag_random`; `uniform` uses `v0x`, `v0y`, and `v0z`, while `isotropic_tag_random` is deterministic across MPI decompositions. |
| `particle_seed` | `0` | Seed combined with species and tag for `tag_random` and `isotropic_tag_random`. |
| `v0x`, `v0y`, `v0z` | `1, 0, 0` | Uniform particle velocity used when `particle_velocity = uniform`. |
| `B0x`, `B0y`, `B0z` | `0, 0, 1` | Uniform magnetic field initialized by `part_random`. |
| `B_profile` | `uniform` | `uniform`, `linear_cross`, `sinusoidal_divb_free`, `mirror`, `gradb`, or `turbulent`; the nonuniform profiles are prescribed accuracy fields. |
| `Bgrad` | `1.0` | Gradient strength for `linear_cross`, `mirror`, and `gradb` manufactured fields. |
| `Bamp`, `Bwave_number` | `0.05`, `1.0` | Amplitude and wave number for `sinusoidal_divb_free` and `turbulent` profiles. |

### `<mesh_refinement>` for Particle AMR Runs

| Parameter | Default | Description |
|-----------|---------|-------------|
| `particle_load_weight` | `0.0` | Additional AMR load-balancing cost per particle.  Zero preserves mesh-only balancing. |

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
| `boris` + `trilinear` | `<mhd>` block | Uses full multidimensional trilinear interpolation of cell-centered magnetic fields. |
| `boris` + `lin` or `lin_legacy` | `<mhd>` block | Uses legacy component-wise interpolation of face-centered magnetic fields. |

The Boris pusher requires MHD because it samples magnetic fields.  The
interpolation choice is made through `<particles>/interpolation`.  All Boris
modes share the same pusher update; only the field-gather policy changes.  Use
`lin` for backward comparisons, `trilinear` for a low-cost point-centered
gather, and `tsc` for the smoother default.

The regression suite includes a smooth-field convergence check using
`B_profile = linear_cross`.  In that test, `trilinear` recovers the analytic
point-centered field to roundoff for a linear cross-field, while `lin` shows
the expected cross-direction interpolation error that decreases with
resolution.  The test is intentionally a pusher/gather check, not an MHD
evolution test.

## Particle Subcycling

Subcycling is off by default so existing CR tracer inputs preserve their
previous timestep behavior.  Set `<particles>/subcycle = true` when particle
motion in one mesh timestep would make the field gather stale or when a Boris
orbit needs a smaller gyro-angle step.

When enabled, AthenaK computes a global substep count from three constraints:

- `subcycle_cell_fraction`: maximum particle displacement as a fraction of the
  local cell width,
- `subcycle_meshblock_fraction`: maximum particle displacement as a fraction of
  the owning MeshBlock size,
- `subcycle_gyro_fraction`: maximum Boris gyro angle per substep.

The active substep count is the maximum over local particles and MPI ranks.
Between intermediate substeps, the particle module updates `PGID` and runs the
same particle MPI exchange used by the normal task list.  This keeps the next
field gather tied to the current owning MeshBlock instead of waiting until the
end of the full mesh timestep.  The final substep still hands off to the normal
`NewGID`/exchange tasks.

With `subcycle_strict = true`, the run fails fast if the constraints require
more than `subcycle_max_steps`.  With `subcycle_strict = false`, the run caps at
`subcycle_max_steps`; pair that mode with `check_motion_bounds = true` while
debugging.  When `log_performance = true`, rank zero prints the selected
substep count and active constraint.

## AMR Remapping

AMR changes the live leaf MeshBlock list and can also move MeshBlocks between
MPI ranks.  After AMR, `Particles::RemapAfterAMR()` remaps particles by
position onto the new leaf list before handing off-rank particles to the normal
particle MPI exchange.

The default remap path is:

1. Build a flattened finest-level logical lookup table from the live AMR leaf
   list.
2. In a Kokkos kernel, wrap periodic particle coordinates into the mesh domain.
3. Convert each wrapped position to a finest-level logical table index using
   `Mesh::max_level`.
4. Assign the table's leaf MeshBlock `gid` to `PGID`.
5. Build the particle send list for particles whose new `PGID` belongs to
   another MPI rank.
6. Exchange rank-count metadata with `MPI_Alltoall`.
7. Pack, send, receive, unpack, and compact particles.
8. Update particle counts.

This avoids the previous host-side scan over every MeshBlock for every particle.
Set `<particles>/amr_remap = host_tree` to use the host mesh-tree fallback.
Set `<particles>/exchange_mode = allgather` to use the older detailed metadata
exchange path.  When `validate_amr_lookup = true`, or when local/full
consistency checks are active, device-table results are checked against
`MeshBlockTree::FindLeafContaining()` and the returned MeshBlock bounds.

The particle remap is separate from face-centered MHD AMR repair.  It does not
modify `RefineFC()`, `RestrictFC()`, or `RepairAMRFC()`, and it is called only
after mesh metadata and physics arrays have been rebuilt.

The same position-based lookup is also used by the post-push particle
`NewGID` task on multilevel meshes, so AMR ownership after particle motion is
checked against the live leaf list rather than an immediate-neighbor stencil.
Uniform meshes keep the cheaper neighbor-stencil update.

## Particle-Aware Load Balancing

AMR load balancing remains mesh-only by default.  Set
`<mesh_refinement>/particle_load_weight` to a positive value to add particle
counts to the MeshBlock cost estimate used after refinement and derefinement.
The cost model is:

```text
cost = 1 + particle_load_weight * estimated_particles_in_block
```

The estimate uses the pre-AMR particle owners.  When a MeshBlock is refined,
its particle count is split evenly among the new children for load-balancing
purposes; after the mesh has been rebuilt, `RemapAfterAMR()` still remaps each
particle by its actual position.  When MeshBlocks are derefined, child particle
counts are assigned to the new parent.  This makes particle-rich regions more
expensive without changing particle ownership, particle counts, or the div(B)
AMR repair path.

Recommended use:

- leave `particle_load_weight = 0.0` for baseline mesh-only comparisons,
- start with a small value such as `0.001` or `0.01` for particle-heavy AMR
  tests,
- compare particles per rank before and after rebalancing using restart or
  inspector output,
- keep consistency checks enabled while tuning the weight.

## Consistency Checks

`<particles>/check_consistency = true` is a legacy alias for
`<particles>/check_consistency_mode = full`.  Consistency checks run after
problem generation, particle restart load, AMR remap, particle MPI exchange, and
particle restart output.

The modes are:

- `none`: no extra host-side checks,
- `counts`: conserved total and per-species counts,
- `local`: count checks plus local `PGID`, rank, position, finite-value, and
  local duplicate-tag checks,
- `full`: local checks plus global duplicate-tag checks across MPI ranks.

The full checks validate:

- every particle `PGID` maps to a live MeshBlock on the owning rank,
- every particle position lies inside its assigned MeshBlock,
- total and per-species particle counts remain conserved,
- tags are unique within each species,
- particle real fields are finite,
- particle restart files contain exactly the declared header and data length.

These checks copy particle data to host and use MPI reductions where needed.
They are intended for regression tests and debugging, not production runs.

`<particles>/check_motion_bounds = true` is a cheaper guardrail for particle
exchange assumptions.  It fails immediately after a push if a particle has
moved outside the old one-neighbor MeshBlock stencil.  This option detects
large timesteps; use `subcycle = true` when large particle velocities require
intermediate ownership updates or smaller Boris gyro-angle steps.

## Outputs

| File type | Output directory | Description |
|-----------|------------------|-------------|
| `ppd` | `ppd/` | Shared compact particle species and position dump. |
| `prst` | `prst/rank_XXXXXXXX/` | Per-rank particle restart dump. |
| `df` | `df/` | Pitch-angle cosine histogram. |
| `pspec` | `pspec/` | Per-species one-dimensional spectrum. |
| `pspec2` | `pspec2/` | Per-species two-dimensional spectrum. |
| `dxh` | `dxh/` | Displacement histograms for `IPDX`, `IPDY`, and `IPDZ`. |
| `drh` | `drh/` | Scalar displacement histogram for `sqrt(IPDX^2 + IPDY^2 + IPDZ^2)`. |
| `dparh` | `dparh/` | Parallel-displacement histogram using `IPDB`. |
| `pmom` | `pmom/` | Reduced per-species pitch-angle, displacement, flux, and tensor moments. |
| `psamp` | `psamp/rank_XXXXXXXX/` | Deterministic selected-field particle sample dump. |
| `pvtk` | run directory | Particle VTK output. |
| `trk` | run directory | Tracked-particle output. |

`prst` files store a three-`Real` header `(time, dt, particles_this_rank)`,
followed by 17 `Real` values per particle: `PGID`, `PTAG`, `PSP`, and the 14
real particle fields.

For `df`, `pspec`, `dxh`, `drh`, `dparh`, and `pmom`, `reduce = true` is the default.
MPI runs write one globally reduced diagnostic from rank zero.  The sum of each
histogram species row equals that species particle count.  Each `dxh`
species/component histogram also sums to that species count.

Set `df_single_file_per_rank = 1` or `dxhist_single_file_per_rank = 1` to
preserve per-rank histogram output.  Those options force `reduce = false`.
The spectrum histogram uses `pspec_single_file_per_rank = 1`, scalar histograms
use `drh_single_file_per_rank = 1` and `dparh_single_file_per_rank = 1`, and
particle moments use `pmom_single_file_per_rank = 1`.

`pspec` uses stored particle velocity and sampled pusher fields.  The `p`
quantity is the tracer speed proxy `sqrt(vx^2 + vy^2 + vz^2)`, `E` is the
nonrelativistic proxy `0.5 p^2`, `logE` is `log10(E)`, `mu` is the pitch-angle
cosine against the stored sampled magnetic field, and `magnetic_moment` is the
proxy `v_perp^2 / B`.  These are diagnostics of the tracer state; they do not
resample the current MHD fields during output.

`pspec2` writes reduced two-dimensional histograms.  Supported `quantity`
values are `mu_p` for `f(mu,p)`, `mu_E` for `f(mu,E)`, and `vpar_vperp` for
`f(v_parallel,v_perp)`.  Use `nbin1`, `nbin2`, `vmin1`, `vmax1`, `vmin2`, and
`vmax2` to control binning.  The histogram layout is species-major, then first
coordinate bin, then second coordinate bin.  As with `pspec`, reduced MPI
histograms are the default.

`pmom` writes globally reduced per-species means.  In addition to count,
`<mu>`, `<mu^2>`, anisotropy, displacement means, parallel-displacement means,
and `<v^2>`, it includes `<v_i>` as a CR flux/current proxy,
`<v_i v_j>` as a pressure-tensor proxy, and cross-displacement moments
`<dx_i dx_j>` as a diffusion-tensor proxy.  These are tracer moments, not
coupled CR-fluid source terms.

`psamp` writes a text sample of selected particle fields to one file per rank.
It is intended for debugging and lightweight analysis, not for restart.  The
selection is deterministic:

```text
splitmix64(species, tag) % sample_stride == sample_offset
```

Use `species = -1` to sample all species, or set `species` to one species
index.  `sample_stride = 1` writes every local particle; larger values write a
stable subset.  `fields` is a comma-separated list.  Supported stored fields
are `x`, `y`, `z`, `vx`, `vy`, `vz`, `mass`, `bx`, `by`, `bz`, `dx`, `dy`,
`dz`, and `dpar`.  Supported derived fields are `speed`, `energy`, `bmag`,
`mu`, `vpar`, `vperp`, and `magnetic_moment`.  Every row always includes
`rank`, `pgid`, `tag`, and `species`, so those do not need to be listed in
`fields`.

Example:

```text
<output9>
file_type     = psamp
dt            = 0.02
species       = -1
sample_stride = 4
sample_offset = 0
fields        = x,y,z,vx,vy,vz,bx,by,bz,dx,dy,dz,dpar,mass,mu,bmag
```

The file is self-describing: each block records the simulation time, cycle,
rank count, selection parameters, sample count, field count, and column names.

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
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
```

The utility reads `prst`, `ppd`, and `psamp`, reports totals by rank, species,
and MeshBlock `PGID`, reports restart format metadata, validates optional
`*.prst.pmeta` sidecar payload sizes and checksums, verifies restart
header/data lengths, checks finite particle data, checks unique tags within
each species, validates reduced `df`/`pspec`/`pspec2`/`dxh`/`drh`/`dparh`
histogram sums when `--histogram-nbin`, `--spectrum-nbin`, or
`--joint-spectrum-nbin` is supplied, summarizes `pmom` files when present, and
reports sampled-field totals and field names when `psamp` output is present.

The particle regression tests import this utility directly so command-line and
automated checks use the same parser.

## Performance Logging

Set `<particles>/log_performance = true` to print low-overhead counters from
the particle push, MPI exchange, AMR remap, and diagnostic-output paths.  A
minimal smoke command is:

```bash
./athena -i inputs/particles/random_particle_drift.athinput \
  particles/log_performance=true \
  particles/ppc=0.001953125 \
  problem/particle_position=center \
  problem/particle_velocity=uniform \
  time/nlim=1
```

On the 2026-05-22 local validation build, this reported:

```text
Particle performance push: pushed=1 sent=0 received=0 remapped=0 diagnostics=0
```

These counters are intended for local regression and scaling comparisons; they
are not portable performance guarantees.

## CPU/MPI Performance Comparison

The local CPU/MPI comparison is documented in
[CR Tracer CPU/MPI Performance Comparison](cr_tracer_cpu_mpi_performance.md).
The 2026-05-22 benchmark used a Release MPI build with the Kokkos Serial
backend on an Apple M4 Max and three repeated runs per case.  The detailed page
includes the compared implementation states, fix-by-fix interpretation,
rank-scaling notes, a four-mode follow-up matrix for `device_table` versus
`host_tree` and `alltoall_counts` versus `allgather`, reproducibility commands,
bottleneck summary, and a GPU follow-up checklist.

Key results:

| Benchmark | MPI ranks | Old scan/full-copy CPU s | Current CPU s | Speedup |
|-----------|-----------|--------------------------|---------------|---------|
| Remap-focused drift AMR | 1 | 0.5686 | 0.00386 | 147.4x |
| Remap-focused drift AMR | 2 | 0.4210 | 0.00540 | 77.9x |
| Remap-focused drift AMR | 4 | 0.2786 | 0.01233 | 22.6x |
| Boris+MHD AMR end-to-end | 1 | 1.1082 | 0.7974 | 1.39x |
| Boris+MHD AMR end-to-end | 2 | 0.7437 | 0.5433 | 1.37x |
| Boris+MHD AMR end-to-end | 4 | 0.4789 | 0.3451 | 1.39x |

The later position-only AMR host-copy optimization was neutral within noise in
the small CPU Boris benchmark after tree lookup was already present, but it
removes unnecessary remap data movement and is expected to matter more on GPU
builds.  GPU timing remains a TODO.

The same comparison found no significant local runtime penalty for reduced
`df`/`pspec`/`pspec2`/`dxh`/`drh`/`dparh`/`pmom` output under four-rank MPI.
The default reduced mode wrote one file per diagnostic; per-rank mode wrote
one file per rank per diagnostic.

The follow-up architecture matrix is intentionally CPU/MPI-only.  On this
machine the `device_table` AMR remap path is 1.13x-1.63x faster than
`host_tree` in the remap-focused drift benchmark and neutral to 1.02x faster in
the Boris+MHD benchmark.  The `alltoall_counts` metadata path is neutral within
local noise in the drift, exchange, and Boris+MHD matrices, but remains the
default because it removes detailed all-gather metadata from the exchange
design.  GPU timing remains a separate TODO.

The benchmark harness is:

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

It runs existing particle inputs with command-line overrides and records
AthenaK CPU time, wall time, particle update rate, AMR activity, restart
readback totals, and particle exchange message/byte estimates.

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
- smooth-field `trilinear` gather convergence against legacy `lin`,
- uniform-field Boris speed conservation,
- periodic drift displacement with wrapping,
- large-motion guard failure,
- reduced `df`/`pspec`/`dxh`/`drh`/`dparh` histogram sums,
- reduced `pspec2` joint-spectrum sums,
- reduced `pmom` count, flux, and tensor checks,
- deterministic selected-field `psamp` readback,
- multi-species `ppd` and `prst` count/tag checks.

The MPI suite covers:

- two-rank Boris AMR particle migration,
- four-rank moving-box AMR refine/derefine stress,
- total and per-species conservation across rank exchange,
- MPI particle restart write/read,
- reduced `pmom`, `pspec`, `pspec2`, `drh`, and `dparh` validation under MPI,
- deterministic selected-field `psamp` readback under MPI,
- consistency-check validation during AMR remap and MPI exchange.

Local validation for this branch on 2026-05-22:

| Command | Result |
|---------|--------|
| `test_particles_cr_cpu.py --cpu` | 20 passed. |
| `test_particles_cr_mpicpu.py --mpicpu` | 4 passed. |
| `python run_test_suite.py --style` | 2 passed. |
| `python run_tests.py mhd/mhd_divb_amr` | Passed through 3D physical AMR level 5. |
| `git diff --check` | Passed. |
| `-D PROBLEM=part_random` build | Passed. |

The explicit four-rank stress run reported `105 MeshBlocks created, 77 deleted
by AMR`.  The inspector reported 1024 total particles, rank counts
`0=238 1=259 2=267 3=260`, species counts `0=512 1=512`, valid
`AKPRST-v1` restart metadata, and valid reduced
`df`/`pspec`/`pspec2`/`dxh`/`drh`/`dparh` histograms.

Div(B) AMR compatibility inputs were rerun against the same built-in pgen
build:

| Input | AMR result |
|-------|------------|
| `inputs/tests/divb_amr_1d.athinput` | 69 MeshBlocks created, 50 deleted. |
| `inputs/tests/divb_amr_2d.athinput` | 1311 MeshBlocks created, 1095 deleted. |
| `inputs/tests/divb_amr_3d.athinput` | 3633 MeshBlocks created, 2842 deleted. |

The dedicated normalized-divergence regression,
`python run_tests.py mhd/mhd_divb_amr`, also passed its 1D, 2D, and 3D
physical-level-1-through-5 cases.  The deepest 3D case reached `30,829` live
MeshBlocks after `50,673` creations and `19,908` deletions while remaining
within the test's normalized `div(B)` tolerances.

## Integration Notes

The GitHub Pages branch already includes this page in the module toctree as
`docs/source/modules/particles.md`.  To publish these docs there, copy or merge
this file and the linked CR tracer Markdown pages into
`docs/source/modules/` on `origin/gh-pages`; no `index.md` toctree edit is
needed because the supporting pages are reached from this page.

These linked pages summarize validation evidence and recommended next CR tracer
work:

- [CR Tracer CPU/MPI Performance Comparison](cr_tracer_cpu_mpi_performance.md)
  records the local CPU and MPI performance evidence for the AMR remap and
  diagnostic-output changes.
- [CR Tracer GPU Testing Handoff](cr_tracer_gpu_testing_handoff.md)
  gives the exact accelerator validation checklist that still needs to run on a
  GPU machine.
- [CR Tracer Accuracy Validation](cr_tracer_accuracy.md)
  records the 12-test CPU/MPI accuracy ladder, current quantitative results,
  and documentation figures.
- [CR Tracer Accuracy Test Implementation Plan](cr_tracer_accuracy_test_plan.md)
  lays out the proposed analytic, AMR, MPI, and ensemble accuracy tests and the
  documentation figures they should produce.
- [CR Tracer Feature-Branch Hardening Plan](cr_tracer_feature_branch_plan.md)
  covers changes that fit naturally in `feature/CR_tracers`.
- [CR Tracer Follow-Up Architecture Plan](cr_tracer_followup_architecture_plan.md)
  covers deeper remap, communication, interpolation, and diagnostics work for a
  subsequent feature branch.
- [CR Tracer Relativistic Acceleration Implementation Guide](cr_tracer_relativistic_acceleration_implementation_guide.md)
  defines the gated implementation and validation workflow for an opt-in
  passive relativistic acceleration follow-up.
