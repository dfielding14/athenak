# AGENTS.md

## Purpose
This directory implements the particle module: data storage for particle arrays,
particle initialization, pushers (drift, RK4 gravity, Boris), and integration with
particle boundary exchange handled in `src/bvals`.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core class and setup
- `particles.hpp`: `particles::Particles` class definition, pusher/type enums, task IDs.
- `particles.cpp`: constructor, input parsing, initialization for cosmic rays and stars,
  tag assignment, species setup, and allocation of particle arrays.

### Pushers and field interpolation
- `particles_pushers.cpp`: pusher implementations (drift, RK4 gravity, Boris) plus
  `InterpolateLinear` and `InterpolateTSC` for B-field sampling.

### Task list wiring
- `particles_tasks.cpp`: inserts particle tasks into the `before_timeintegrator`
  and `after_timeintegrator` task lists, and wires coupled-stage insertion points.
- `particles_moments.cpp`: moment deposition/communication wrappers and the
  deterministic cell-centered to edge-current conversion task used by PR2
  `edge_staggered` coupling mode.

### Data structs
- `particles_data_structs.hpp`: `ParticleLocationData` and `ParticleMessageData` used
  by particle MPI exchange in `src/bvals/bvals_part.cpp`.

---

## Data Layout and Indices
Particles store two Kokkos 2D arrays with index order `(var, particle)`:
- `prtcl_rdata` (Real): continuous particle properties
- `prtcl_idata` (int): integer properties

Index constants are defined in `athena.hpp`:
- **Integer indices** (for `prtcl_idata`):
  - `PGID`: owning MeshBlock global ID
  - `PTAG`: unique tag
  - `PSP` (cosmic rays) or `NSN` (stars): species or SN count (same index value)
- **Real indices** (for `prtcl_rdata`):
  - `IPX, IPVX, IPY, IPVY, IPZ, IPVZ`: position and velocity
  - Cosmic rays (`CRParticlesIndex`): `IPM` (q/m), `IPBX/IPBY/IPBZ` (sampled B),
    `IPDX/IPDY/IPDZ` (displacement), `IPDB` (parallel displacement),
    `IPEX/IPEY/IPEZ` (sampled midpoint `cE`), `IPDPX/IPDPY/IPDPZ` (per-step
    momentum-rate feedback channels), `IPDE` (per-step energy-rate channel),
    `IPEBDOT` (midpoint frozen-in orthogonality diagnostic `cE dot B`),
    `IPWT` (relative macro-particle weight)
  - Stars (`StarParticlesIndex`): `IPT_CREATE`, `IPMASS`, `IPT_NEXT_SN`

The number of real slots (`nrdata`) is chosen at construction time:
- Cosmic rays:
  - Always provisioned as `nrdata = IPWT + 1` so midpoint E+B, coupled feedback
    diagnostics, and per-particle macro weights are available in all staged PIC
    modes.
- Stars: `nrdata = 9` (positions, velocities, and three star fields).

`nidata` is always 3 (`PGID`, `PTAG`, and `PSP` or `NSN`).

---

## Particle Types

### Cosmic rays (`particle_type = cosmic_ray`)
- Particle count: `ppc` (particles per cell). Counts are computed per global
  MeshBlock from cumulative `floor(gid * ppc * nx1 * nx2 * nx3)` differences
  so total particle count is independent of MPI rank decomposition, including
  fractional `ppc`.
- Species:
  - `nspecies` in `<particles>` block.
  - Species properties in blocks `species0`, `species1`, ... with `mass` and
    `charge`.
  - Optional per-species drift overrides in each `speciesN` block:
    `vx0`, `vy0`, `vz0` (fallback to global `cr_vx0/cr_vy0/cr_vz0`).
  - `IPM` stores charge-to-mass ratio for each particle.
- Initialization:
  - `cr_distribution = center` (default) maps particles over MeshBlock cell
    centers; it does not collapse every particle in a block onto the block
    center.
  - `cr_distribution = random` uses deterministic per-global-block hashing
    controlled by `pic_random_seed` to place uniformly in each block. Species
    are assigned round-robin per block, and particles with the same
    `position_index = pinmb/nspecies` share a position across species.
  - `IPWT` is initialized as local cell volume divided by root-cell volume, so
    ppc-created particles carry refinement-consistent macro weights while
    preserving root-level `deposit_qscale` semantics.
- Optional displacement tracking via `track_displacement` updates `IPD*` and `IPDB`.
- Boris pushers support both MHD-carried fields (`coupled` / `passive_mhd`) and
  particle-owned no-MHD carriers (`pic_background_mode=no_mhd`).

### Stars (`particle_type = star`)
- **3D only**. Constructor aborts if `three_d` is false.
- Particles are loaded from `star_particle_file` with 8 columns per line:
  `x y z vx vy vz t_create mass` (comments starting with `#` are skipped).
- Each rank reads the file and keeps only particles inside its local MeshBlocks.
- Additional fields:
  - `NSN` (integer) counts supernovae for that particle.
  - `IPT_NEXT_SN` is initialized via `GetNthSNTime` from `utils/sn_scheduler.hpp`.

---

## Pushers and Interpolation

### Pusher selection (`<particles>` block)
Supported strings in code:
- `drift` -> `PushDrift`
- `rk4_gravity` -> `PushStars`; valid only with `particle_type=star`
- `boris_lin` -> `PushCosmicRays` with linear interpolation
- `boris_tsc` -> `PushCosmicRays` with TSC interpolation
- Boris pushers are valid only with `particle_type=cosmic_ray`.

Note: the enum also lists `leap_frog`, `lagrangian_tracer`, and `lagrangian_mc`,
but these are not wired in the constructor.

### Pushers
- **Drift**: updates positions with velocities using `mesh->dt`.
- **RK4 gravity** (stars): integrates in an analytic potential using finite-difference
  gradients of `GravPot`. Parameters are read from the `<potential>` block:
  `r_scale`, `rho_scale`, `mass_gal`, `scale_gal`, `z_gal`, `r_200`, `rho_mean`.
  Step size uses `particles:grav_dx` (default `1e-6`).
- **Boris** (cosmic rays): midpoint E+B Boris sequence:
  first drift by `dt/2`, interpolate midpoint `u` and `B`, compute
  `cE = -u x B`, apply full-step Boris momentum update, then second drift by
  `dt/2`. Per-step `dp/dt`, `dE/dt`, and `cE dot B` diagnostics are stored in
  the CR payload for deposition and regression checks.

### Field interpolation
- `InterpolateLinear`: trilinear (or bilinear in 2D) interpolation of midpoint
  carrier fields (`B` and fluid velocity `u`).
- `InterpolateTSC`: triangular-shaped cloud weighting over a 3x3x3 stencil for
  midpoint carrier fields (`B` and fluid velocity `u`).

---

## Moment Deposition and PR2 Coupling Controls

### Deposition controls
- `deposit_moments` (default `false`): enables particle moment deposition.
  In Boris CR paths this includes `rho/J` and midpoint diagnostics channels
  (`E dot B`, and in coupled mode `dp/dt`, `dE/dt`).
- `deposit_order`:
  - default and all non-direct paths: `1` only
  - coupled `direct_staggered` path: `1` and `2` supported
  - coupled `direct_staggered` deposition is trajectory-based and aborts if the
    old-to-new particle shape support shifts by more than one cell in any active
    dimension; lower `time/cfl_number` or `pic_max_cell_cross` if that guard trips.
- `deposit_qscale`: macro-charge scaling.
- Restart files persist the particle record for every active particle mode.
  When `deposit_moments=true`, restart files also persist deposited moment
  arrays; coupled edge-staggered runs additionally persist `j_edge_x*e`.

### PR2 E-field coupling controls
- `couple_moments_to_mhd` (default `false`): opt-in particle-to-MHD coupling.
- `couple_j_to_efield_coeff` (default `1.0`): explicit current-to-E coupling
  coefficient applied in MHD `EFieldSrc`.
- `couple_j_to_efield_representation`:
  - `cell_centered` (default): uses deposited CC current directly.
  - `edge_staggered`: stores current on the edge-centered layout consumed by
    `MHD::EFieldSrc`.

### PR2 fluid feedback controls
- `couple_moments_momentum_to_mhd` (default `false`)
- `couple_moments_energy_to_mhd` (default `false`)
- `couple_moments_momentum_coeff` / `couple_moments_energy_coeff`
- `couple_fluid_feedback_order`:
  - `mhd_src_terms` (default)
  - `efield_src` (parity experiment mode).

### Deterministic CR initialization knobs
- `cr_distribution=center` maps initialized positions deterministically onto cell
  centers across each MeshBlock. It no longer collapses all particles in a block
  onto the MeshBlock center; extra particles wrap over cell centers.
- `cr_vx0`, `cr_vy0`, `cr_vz0` default to `0.0` and are used by PR2 regression
  tests to make integrated current expectations deterministic.
- `speciesN/vx0`, `speciesN/vy0`, `speciesN/vz0` optionally override those
  global drifts for species-specific beam initialization without changing
  legacy single-drift defaults.

### Staged PIC runtime controls (Step 1-4 status)
- `pic_background_mode`: `coupled` (default), `passive_mhd`, `no_mhd`.
- `pic_feedback_mode`:
  - default `coupled` for `pic_background_mode=coupled`
  - default `test_particle` for `pic_background_mode=passive_mhd`
  - accepted values: `coupled`, `test_particle`
- `pic_interp_scheme`: `tsc` (default and currently only valid value).
- `pic_enable_2d3v`: required for Boris pushers on 2D meshes; the reduced
  2D/2V Lorentz-force path is not implemented.
- `pic_cr_light_speed` is reserved for a future reduced-speed-of-light pusher;
  only the default value `1.0` is accepted.
- `pic_max_cell_cross` (default `2`) and `pic_theta_max` (default `0.3`)
  with positivity guards. `pic_max_cell_cross` must not exceed the smallest
  active MeshBlock dimension because particle exchange is nearest-neighbor.
- `pic_deltaf_mode`: `off` (default) or `on`; `on` requires
  `pic_deltaf_f0` to be explicitly set to `kappa_iso` or `uniform_quiet`.
  - staged behavior: with `cr_distribution=random`, `pic_deltaf_mode=on`
    applies deterministic low-discrepancy quiet-start particle placement for
    reduced sampling noise in proxy instability tests. It does not yet apply a
    full delta-f particle-weight evolution or otherwise use `pic_deltaf_f0` in
    moment deposition.
- `pic_sort_interval` (default `0`, must be `>= 0`).
- `pic_random_seed` (default `0`, must be `>= 0`) controls deterministic
  `cr_distribution=random` particle placement.
- `pic_intermediate_arrays`: `auto` (default) or `off`.
- `pic_expanding_box_mode`: `off` (default) or `on`.
- `pic_expansion_rate_x1/x2/x3` default to `0.0`; non-zero expansion rates
  require `pic_expanding_box_mode=on`.
  - staged behavior: when expanding-box mode is on, Boris pushers apply
    per-component half-step source scaling around the Boris update using
    `exp(-0.5*dt*pic_expansion_rate_xi)`.
- `pic_no_mhd_bx`, `pic_no_mhd_by`, `pic_no_mhd_bz` define uniform no-MHD
  Boris background fields for the particle-owned carrier.
- `passive_mhd` currently requires:
  - active `<mhd>` block
  - `pic_feedback_mode=test_particle`
  - coupling toggles disabled (`couple_moments_to_mhd`,
    `couple_moments_momentum_to_mhd`, `couple_moments_energy_to_mhd`).
- `no_mhd` currently requires:
  - `pic_feedback_mode=test_particle`
  - coupling toggles disabled (`couple_moments_to_mhd`,
    `couple_moments_momentum_to_mhd`, `couple_moments_energy_to_mhd`)
  - Boris pushers use `pic_no_mhd_bcc0` instead of `pmhd->bcc0`; the
    particle-owned carrier is allocated to AMR `max_nmb_per_rank` capacity.

### Runtime guards (constructor)
- Coupling requires `deposit_moments=true` and an active `<mhd>` block.
- Coupled mode is rejected for `radiation+MHD`, hydro/ion-neutral, and
  numerical-relativity (`adm`/`z4c`) compositions.
- Energy feedback requires ideal MHD EOS.
- `edge_staggered` and fluid feedback branches are restricted to
  non-relativistic MHD in PR2.
- `pic_feedback_mode=test_particle` explicitly rejects particle-to-MHD coupling
  toggles in the current staged implementation.
- Boris pushers now require either active MHD fields or
  `pic_background_mode=no_mhd` with a valid no-MHD carrier path.

---

## Boundary Exchange (MPI)
Particle communication is implemented in `src/bvals/bvals_part.cpp` via
`ParticlesBoundaryValues`, but is driven by particle tasks in this directory.

High-level flow:
1. **Push**: advance particle positions/velocities.
2. **NewGID**: detect boundary crossings, update `PGID`, and build send/destroy lists.
3. **Count**: share send counts across ranks.
4. **InitRecv**: post non-blocking receives for particle buffers.
5. **SendP**: pack and send particle data (real + int buffers).
6. **RecvP**: unpack received particles, fill holes, destroy out-of-domain particles.
7. **ClearRecv/ClearSend**: finalize MPI requests.

Notes from `bvals_part.cpp`:
- `BoundaryFlag::reflect` is handled immediately after each pusher by mirroring
  the particle position and flipping the normal velocity before deposition or
  exchange.
- Any particle that still exits through a non-periodic physical face during
  exchange (`outflow`, `inflow`, `user`, or an unhandled physical crossing) is
  marked for destruction.
- Periodic boundaries wrap positions back into the global domain.
- Arrays are resized when receives exceed sends, and compacted after destruction.
- `mesh->CountParticles()` is called after exchanges to update global counts.

---

## Task List Integration
`Particles::AssembleTasks` always wires:
- `before_timeintegrator` push/deposition chain:
  `SaveOldPositions -> Push -> ZeroMoments -> InitRecvMoments -> DepositMoments ->
  SendMoments -> RecvMoments -> ClearRecvMoments -> ClearSendMoments`

Particle migration communication depends on coupling mode:
- uncoupled/default: `NewGID -> SendCnt -> InitRecv -> SendP -> RecvP ->
  ClearRecv -> ClearSend` in `before_timeintegrator`
- coupled: same migration chain moved to `after_timeintegrator`.

In coupled mode, moment wrappers are also inserted into `stagen` on stage 1:
- insertion anchor is selected by `couple_fluid_feedback_order`
  (`MHD::MHDSrcTerms` vs `MHD::EFieldSrc`)
- if `couple_j_to_efield_representation=edge_staggered` and
  `couple_j_deposition_mode=cc_convert`, `ConvertCoupledCurrentRepresentation`
  is inserted immediately before `MHD::EFieldSrc` and depends on both `CornerE`
  and wrapper completion.
- if `couple_j_deposition_mode=direct_staggered`, direct edge-current
  synchronization and physical-BC tasks are inserted before `MHD::EFieldSrc`
  instead of the CC conversion task. Physical edge-current BCs cover
  `periodic`, `reflect`, and `outflow`; direct mode rejects `inflow`, so use
  `couple_j_deposition_mode=cc_convert` for inflow-boundary coupled runs.

---

## Outputs and Diagnostics (references)
- `outputs/vtk_prtcl.cpp` and `outputs/track_prtcl.cpp` consume `prtcl_rdata` and
  `prtcl_idata` for particle dumps and tracked particles.
- `trk` output selects particles by nonnegative `PTAG < <output>/nparticles` and
  writes big-endian float32 position/velocity rows at tag-derived offsets.
- Derived variable `prtcl_d` (in `outputs/derived_variables.cpp`) bins particle
  counts onto the mesh.

---

## Extension Points and Cautions
- Any new particle type or pusher must update:
  - enum(s) and constructor parsing in `particles.cpp`.
  - data layout decisions for `nrdata` and `nidata`.
  - push logic in `particles_pushers.cpp`.
- Boris pushers support no-MHD mode only through the explicit
  `pic_background_mode=no_mhd` carrier path; keep constructor guards strict.
- Star particles assume 3D and use an external file for initialization; keep the
  file format consistent.
- Particle boundary exchange relies on `PGID` and the neighbor index scheme; keep
  GID updates and periodic wrapping consistent with mesh BCs.
