# Particles in a frozen CGL checkpoint

Use the executable built from this CGL branch to load its own checkpoint and
advance test particles in the saved magnetic field. This preserves the CGL
restart representation, Landau-fluid diagnostics, turbulence restart payload,
and loaded mesh hierarchy.

This path currently accepts 3D periodic Newtonian single-fluid MHD checkpoints.
Live CGL evolution with particles and frozen runs started without a checkpoint
are rejected.

The Boris pusher integrates $d\mathbf v/dt=\mathbf v\times\mathbf B/m_g$,
where `min_mass` specifies the gyro parameter $m_g$ in code units. The local
gyroradius is $r_L=m_g v_\perp/|\mathbf B|$; choose particle speed and gyro
parameters to represent the desired CR rigidities in the simulation units.
It has no electric force, feedback on the fluid, or stochastic scattering.
The saved fluid velocity does not enter this equation. The magnetic field is
sampled from the loaded mesh, including its initial ghost fill; `problem/B0*`
does not replace it.

## Seed a particle run

Build with the usual CGL configuration and built-in problem generators; no
separate `PROBLEM=part_random` build is required. For example:

```sh
./build/src/athena -r /absolute/path/cgl.00100.rst \
  -i inputs/particles/frozen_cgl.athinput -d cr-frozen \
  job/basename=cr-frozen time/tlim=1.2
```

`time/tlim` is the absolute final time and must exceed the time in the checkpoint.
The overlay keeps time evolution enabled so particles advance. Its
`particles/frozen_mhd=true` stops fluid updates, LF/STS, collisions, forcing,
and adaptive refinement while preserving the loaded fields and hierarchy.
Do not use `time/evolution=static` for a particle advance.

The supplied overlay seeds deterministic positions and isotropic directions.
`particles/ppc` sets particles per active cell per species; the global count
must be nonzero. Empty ranks are supported. `particles/nspecies`, `min_mass`, and
`mass_log_spacing` set $m_s=m_{\min}(\mathrm{mass\_log\_spacing})^s$.
The default `assign_tag=index_order` supports uneven rank populations;
`rank_order` requires equal initial particle counts on every rank.
Masses must be positive. The speed is
$\sqrt{v_{0x}^2+v_{0y}^2+v_{0z}^2}$ from `problem/v0x`, `v0y`, and `v0z`.
`particle_velocity=uniform` uses those vector components directly.
`particle_position=fixed` uses `problem/particle_x1`, `particle_x2`, and
`particle_x3`; `tag_random` samples the full domain reproducibly.

The example uses TSC gathering and no subcycling. `particles/cfl_part` limits
the gyro step using the actual loaded field. A separate motion limit uses
particle speed and cell spacing. Check
trajectory convergence before relaxing the timestep or enabling subcycling for
a production run. Particle speed conservation alone does not check gyro phase.
Transport measurements also need convergence with magnetic-field resolution;
interpolation does not supply unresolved magnetic fluctuations.
Use `particles/check_consistency_mode=full` for short validation runs; the
example leaves these host-side checks disabled during production pushes.

The overlay adds particle checkpoints as `<output99>`. Choose another unused
output block if the source checkpoint already uses that number. Existing fluid
outputs remain configured; adjust their cadence if they are unnecessary for the
particle study.

For position and velocity tracks, add an unused output block with `file_type=trk`,
`dt`, and `nparticles`. It tracks the first `nparticles` global tags. Particle
checkpoints also retain species, gyro parameter, and accumulated displacement.
Tags are grouped by species, so a short track list can contain only the lightest
bin. For balanced samples in frozen Boris runs, use `particles_per_species`
instead of `nparticles`:

```ini
<output98>
file_type = trk
particles_per_species = 1024
dt = 0.001
```

This writes buffered `rich_v1` files under `trk/rank_*/`, with named fields for
compact output tag, original global tag, species, mass, position, velocity,
endpoint magnetic field, unwrapped displacement, and signed displacement along
the local field. The endpoint field uses the pusher's interpolation. It requires
`assign_tag=index_order`, enough particles in every bin, and at most 2²⁴ total
particles so the float32 tag fields remain exact. Particle checkpoints retain
the full ensemble in double precision. Track cadence controls saved samples,
not the Boris timestep; resolve the gyroperiod when studying individual events.

Restart initialization emits no initial outputs. For an exact initial track and
particle sample, run a separate zero-cycle replay with the same seed and combine
its output with the advancing run. Read actual saved times rather than assuming
outputs land exactly on the requested cadence.

For a wide range of gyro parameters, per-particle gyro subcycling avoids making
every bin take the lightest bin's number of substeps. The 768³ pilot used these
overrides to the supplied overlay:

```ini
<particles>
ppc = 1e-4
nspecies = 19
min_mass = 0.00025
mass_log_spacing = 1.41421356
cfl_part = 0.3
subcycle = true
subcycle_max_steps = 64
subcycle_gyro_fraction = 0.05
subcycle_per_particle_gyro = true
exchange_gyro_only_substeps = false
```

The outer timestep still limits travel to half the finest cell, with ownership
exchange after each outer step. The gyro bound uses the maximum saved field;
each mass bin takes enough inner steps to satisfy that bound. Validate with
`check_consistency_mode=full`, `check_motion_bounds=true`, and
`log_performance=true` on short runs, then disable those checks for timing.

## Resume particles

Keep the matching fluid checkpoint and every rank's `.prst` file from the same
time. Resume with the same rank count and mesh layout. The rank-zero particle
path is rewritten to each rank's corresponding directory:

```sh
mpirun -np 2 ./build/src/athena -r /absolute/path/cr-frozen.00120.rst \
  -i inputs/particles/frozen_cgl.athinput -d cr-resumed \
  particles/seed_on_restart=false problem/prtcl_rst_flag=1 \
  problem/prtcl_res_file=/absolute/path/prst/rank_00000000/cr-frozen.00120.prst \
  time/tlim=1.4
```

Match checkpoint **times**, not filename numbers: fluid and particle outputs
may have different cadences. Keep the original particle configuration when
resuming. This path does not redistribute a particle checkpoint to a new rank
count.

## Focused checks

Run the CPU integration checks against a double-precision executable with the
built-in `divb_amr` generator. They create their own small CGL checkpoint:

```sh
ATHENAK_TEST_EXE=/absolute/path/build/src/athena PYTHONPATH=tst \
  python3 -m pytest -q tst/test_suite/cgl/test_cgl_cr_cpu.py \
  tst/test_suite/cgl/test_cgl_cr_smr_cpu.py \
  tst/test_suite/cgl/test_cgl_cr_tracks_cpu.py
```

For an MPI-enabled executable and a working `mpirun`:

```sh
ATHENAK_TEST_EXE=/absolute/path/build-mpi/src/athena PYTHONPATH=tst \
  python3 -m pytest -q tst/test_suite/cgl/test_cgl_cr_mpicpu.py \
  tst/test_suite/cgl/test_cgl_cr_tracks_mpicpu.py
```

The checks cover frozen active and passive CGL fields, a particle-free evolving
control, actual-field gathering and gyro-phase convergence with both interpolation
methods, both subcycling paths, multi-species checkpoint continuation, binary
tracks and VTK output, mismatched checkpoint rejection, and ownership across MPI and
static coarse/fine boundaries.

Initial port validation on 2026-10-17 used Release, double-precision Kokkos Serial builds,
with MPI enabled and separately with MPI disabled and bounds checking enabled.
All nine checks passed with the MPI build; all eight serial checks passed with
bounds checking.
A separate build of the original `3293b6188` produced a byte-identical final
fluid checkpoint and identical histories to the port with particles disabled,
including CGL/LF evolution, forcing, and refinement from two to sixteen blocks.
The subsequent qualification passed all nine serial checks and all eleven
Frontier HIP/MPI checks, including uneven-rank tag validation and analytic
two-mass gyro motion with per-particle subcycling. These checks do not establish
turbulent transport convergence.

The native pilot used the completed active ν=0 and ν=20 snapshots at t=14:
768³ cells in 128 uniformly refined 192×96×192 blocks, on 16 Frontier nodes
with eight ranks/GPUs per node. Each snapshot seeded 858,496 particles, with
45,184 in each of 19 bins spanning gyro parameters 0.00025–0.128 at unit speed.
Original, zero-cycle loaded, and pushed fluid/face/ghost payloads were bitwise
identical, as were topology, forcing/RNG state, LF ledgers and AMR counters.
Particle tags, bin membership, mass, ownership and accumulated displacement
passed their checks. Maximum speed error remained below 9×10⁻¹⁴ through the
2048-cycle timing runs.

For ν=0 at particle time 0.02, halving the gyro fraction from 0.05 to 0.025 changed the
worst bin's mean-square displacement by 0.0017% and its velocity second-moment
components by at most 9.2×10⁻⁶. The predeclared limits were 1% and 0.01,
respectively; this qualifies a short push, not long-time diffusion coefficients.

Two alternating ν=0 timing pairs measured 17.04–17.08 million outer particle
updates/s/node with per-particle subcycling, versus 3.16 million with global
subcycling: a 5.39–5.41× gain. The ν=20 fast run measured 17.11 million.
These are output-free intervals after warmup, excluding checkpoint loading and
final output; an outer update can contain multiple Boris steps. No node-count
scaling sweep was performed. Inputs, binaries, logs and per-bin results are
retained in
`/lustre/orion/ast207/proj-shared/dfielding/CGL/cgl_cr_frozen_768_20261017/`.
