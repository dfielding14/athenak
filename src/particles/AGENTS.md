# AGENTS.md

## Purpose
This directory implements the particle module: data storage for particle arrays,
particle initialization, pushers (drift, RK4 gravity, Boris), and integration with
particle boundary exchange handled in `src/bvals`.
See `../../AGENTS.md` for repository-wide conventions and workflow.
For active MHD-PIC work use `MHD_PIC_NEXT_STEPS_GUIDE.md`; durable model and
implementation documents live under `docs/source/engineering/`.

---

## Key Files and Responsibilities

### Core class and setup
- `particles.hpp`: `particles::Particles` class definition, pusher/type enums, task IDs.
- `particles.cpp`: constructor, input parsing, initialization for cosmic rays and stars,
  tag assignment, species setup, and allocation of particle arrays.

### Pushers and field interpolation
- `particles_pushers.cpp`: pusher implementations (drift, RK4 gravity, Boris)
  plus `InterpolateLinearFields` for carrier-field sampling.
- `field_interpolation.hpp`: shared `InterpolateTSCFields` implementation used
  by Boris pushers and shock injection.

### Task list wiring
- `particles_tasks.cpp`: inserts particle tasks into the `before_timeintegrator`
  and `after_timeintegrator` task lists, and wires coupled-stage insertion points.
- `particles_moments.cpp`: moment deposition/communication wrappers and the
  deterministic cell-centered to edge-current conversion task retained for
  staggered-current diagnostics and engineering paths.

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
  - `PCRSOURCE` (cosmic rays): persistent source cohort (`initial` or
    `shock_injected`)
- **Real indices** (for `prtcl_rdata`):
  - `IPX, IPVX, IPY, IPVY, IPZ, IPVZ`: position plus velocity for drift/star
    particles; Boris cosmic rays always use the second slot of each pair for
    mass-normalized momentum `p/m`
  - Cosmic rays (`CRParticlesIndex`): `IPM` (configured charge-to-mass coupling
    coefficient), `IPBX/IPBY/IPBZ` (sampled B),
    `IPDX/IPDY/IPDZ` (displacement), `IPDB` (parallel displacement),
    `IPEX/IPEY/IPEZ` (sampled midpoint `cE`), `IPDPX/IPDPY/IPDPZ` (per-step
    momentum-rate feedback channels), `IPDE` (per-step energy-rate channel),
    `IPEBDOT` (midpoint frozen-in orthogonality diagnostic `cE dot B`),
    `IPWT` (relative macro-particle weight), `IPF0` (initial analytic delta-f
    background value), `IPDFWT` (evolving delta-f perturbation weight),
    `IPT_BIRTH` (persistent creation time)
  - Stars (`StarParticlesIndex`): `IPT_CREATE`, `IPMASS`, `IPT_NEXT_SN`

The number of real slots (`nrdata`) is chosen at construction time:
- Cosmic rays:
  - Always provisioned as `nrdata = IPT_BIRTH + 1` so midpoint E+B, coupled
    feedback diagnostics, per-particle macro weights, delta-f state, and birth
    time are available for every cosmic-ray configuration.
- Stars: `nrdata = 9` (positions, velocities, and three star fields).

`nidata` is 4 for cosmic rays (`PGID`, `PTAG`, `PSP`, and `PCRSOURCE`) and 3
for stars (`PGID`, `PTAG`, and `NSN`).

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
  - `IPM` stores `species_charge/species_mass`; the Boris update consumes that
    code-normalized coupling coefficient directly.
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
  - CR populations validate all active `IPWT` values once after
    fresh or restart problem setup. Values must be finite and positive; push,
    deposition, and history consumers never repair invalid weights locally.
- Optional displacement tracking via `track_displacement` updates `IPD*` and `IPDB`.
- Boris pushers support both MHD-carried fields (`coupled` / `passive_mhd`) and
  particle-owned no-MHD carriers (`pic_background_mode=no_mhd`).

### Stars (`particle_type = star`)
- **3D only**. Constructor aborts if `three_d` is false.
- A `<units>` block is required because the analytic potential converts code
  units through `Units`.
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
  `cE`, apply the Boris momentum update, then second drift by `dt/2`. The
  standard ideal carrier uses `cE = -u x B`; full CR-Hall uses `-(u+v_H) x B`
  for its particle predictor. With `time/integrator=vl2`, stage 1 makes only a
  scratch half-kick predictor and stage 2 performs the real midpoint full kick.
  The realized stage-2 `dp/dt`, `dE/dt`, and `cE dot B` diagnostics are stored
  in the CR payload for deposition and conservative feedback.

### Field interpolation
- `InterpolateLinearFields`: trilinear (or bilinear in 2D) interpolation of midpoint
  carrier fields (`B` and fluid velocity `u`).
- `InterpolateTSCFields`: triangular-shaped cloud weighting over a 3x3x3 stencil for
  midpoint carrier fields (`B` and fluid velocity `u`).

---

## Moment Deposition and Coupling Controls

### Deposition controls
- `deposit_moments`: enables particle moment deposition and defaults to `false`.
  In Boris CR paths this includes `rho/J` and midpoint diagnostics channels
  (`E dot B`, and in coupled mode `dp/dt`, `dE/dt`).
- `deposit_order`:
  - defaults to `1`
  - `time/integrator=vl2` requires `2` for TSC deposition
  - coupled `direct_staggered` path: `1` and `2` supported
  - coupled `direct_staggered` deposition is trajectory-based and aborts if the
    old-to-new particle shape support shifts by more than one cell in any active
    dimension; lower `time/cfl_number` or `pic_max_cell_cross` if that guard trips.
- `deposit_qscale`: macro-charge scaling.
- Restart files persist the particle record for every active particle mode.
  When `deposit_moments=true`, restart files also persist deposited moment
  arrays; coupled edge-staggered runs additionally persist `j_edge_x*e`.

### Coupling and retained current-representation controls
- `couple_moments_to_mhd` (default `false`): opt-in particle-to-MHD coupling.
- `couple_j_to_efield_coeff` (default `1.0`): retained compatibility coefficient;
  current MHD induction does not consume it in `EFieldSrc`.
- `couple_j_to_efield_representation`:
  - `cell_centered` (default): keeps deposited current in cell-centered moments.
  - `edge_staggered`: allocates edge-current arrays populated by CC conversion or
    direct staggered deposition for retained diagnostics/engineering paths.
- VL2 rejects `edge_staggered` and `direct_staggered` as direct-current CT
  options. Use `pic_cr_hall_mode=off|full`; only `full` adds CR-current-dependent
  induction through MHD face fluxes, `CornerE`, and CT.

### Fluid feedback controls
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
- `cr_vx0`, `cr_vy0`, `cr_vz0` default to `0.0` and are used by focused regression
  tests to make integrated current expectations deterministic.
- `speciesN/vx0`, `speciesN/vy0`, `speciesN/vz0` optionally override those
  global drifts for species-specific beam initialization without changing
  existing single-drift defaults.

### Explicit PIC runtime controls
- `<time>/integrator=vl2` selects the two-stage VL2/TSC coupling algorithm.
  It requires cosmic rays with `pusher=boris_tsc`, a coupled MHD background,
  `pic_feedback_mode=coupled`, `deposit_moments=true`, `deposit_order=2`,
  `couple_moments_to_mhd=true`, conservative momentum feedback, and
  `couple_fluid_feedback_order=mhd_src_terms`. Ideal MHD also requires energy
  feedback; the guarded exact-isothermal paper paths use momentum-only feedback.
  These
  settings are explicit:
  deposition and every particle-to-gas coupling toggle default to `false`, and
  `deposit_order` defaults to `1`.
- Boris cosmic rays always store `p/m`. `pic_cr_initial_state=velocity|momentum`
  controls only whether initializer components are converted from velocity or
  accepted as momentum. `pic_cr_light_speed` is the positive artificial light
  speed used by the relativistic Boris kinematics.
- Fresh inputs containing obsolete `pic_physical_mode` are fatal. A legacy
  restart header may contain that field, but it is warned about and ignored in
  favor of the explicit controls documented here.
- `pic_background_mode=coupled|passive_mhd|no_mhd` selects the field carrier;
  `pic_feedback_mode=coupled|test_particle` independently selects whether
  particle feedback is admitted. `pic_interp_scheme=tsc` is currently the only
  accepted explicit interpolation scheme.
- `pic_enable_2d3v=true` is required for Boris pushers on 2D meshes; there is no
  reduced 2D/2V Lorentz-force path.
- `pic_cr_hall_mode=off|full` selects the CR-Hall closure. `full` requires the
  explicit VL2/TSC full-f contract, uniform-grid ideal MHD, positive
  `pic_background_ion_q_over_mc`, and 3V on non-3D meshes. Its predictor uses
  deposited charge/current; its corrector uses the realized stage-2 particle
  impulse for CT, the matched Hall energy flux, and opposite gas feedback. See
  `docs/source/engineering/pic_cr_hall_code_map.md` for equations and signs.
- `pic_wave_damping_mode=off|ion_neutral_friction` uses a non-negative
  `pic_ion_neutral_collision_rate`; the active option applies the reduced
  static-neutral transverse damping and matched ideal-MHD energy sink.
- `pic_max_cell_cross` (default `2`) cannot exceed the smallest active
  MeshBlock dimension because migration is nearest-neighbor.
  `pic_theta_max` (default `0.3`) limits the relativistic gyro step. The Boris
  timestep scans active and ghost magnetic fields and uses each particle's
  current Lorentz factor; only a globally empty population uses the configured
  species fallback.
- `pic_deltaf_mode=off|quiet_start|physical`. Enabled modes require an explicit
  `pic_deltaf_f0`; `quiet_start` changes sampling only, while `physical` evolves
  and deposits the delta-f weight. The optional
  `pic_deltaf_adapt_mode=global_bikappa_moments_experimental` is limited to its
  guarded expanding-box bi-kappa configuration, and restart schema 8 preserves
  its fitted state and cadence bucket.
- `pic_sort_interval` (default `0`, must be `>= 0`).
- `pic_random_seed` (default `0`, must be `>= 0`) controls deterministic
  `cr_distribution=random` particle placement.
- `pic_load_balance_cost_per_particle` (default `0.0`, must be `>= 0`) adds an
  opt-in particle-count contribution to each post-AMR MeshBlock load cost.
- `pic_intermediate_arrays`: `auto` (default) or `off`.
- `pic_expanding_box_mode`: `off` (default) or `on`.
- `pic_expansion_law`: `linear` (default), `reciprocal_linear`, or
  `exponential`.
- `pic_expansion_rate_x1/x2/x3` default to `0.0`; non-zero expansion rates
  require `pic_expanding_box_mode=on`, and enabled scale factors must remain
  finite and positive through `<time>/tlim`.
  - When expanding-box mode is on, Boris pushers apply exact scale-factor
    half-step ratios around the Boris update and drift in comoving coordinates.
    MHD tasks rescale gas conserved variables, retain raw face-centered arrays
    as divergence-preserving comoving magnetic fluxes, derive physical face
    fields for MHD consumers, and map edge EMFs before constrained transport.
  - Active-MHD expanding-box mode is deliberately fail-closed for unqualified
    compositions, including AMR, non-periodic boundaries, relativistic
    coordinates, nonideal MHD, fluid/radiation/relativity/turbulence blocks,
    staggered/direct current deposition, Hall-current induction, nonzero
    analytic delta-f background current, and user-defined history callbacks.
    The qualified non-delta-f coupled split uses cell-centered `cc_convert` moments,
    conservative momentum and energy feedback, `mhd_src_terms` ordering,
    physical-volume deposition, final-frame particle EM impulses, built-in
    physical-volume MHD `hst`, and optional reduced ion-neutral wave damping.
    The admitted adaptive physical delta-f path retains a separate
    endpoint-normalized analytic `rho E + J x B` source path; nonadaptive
    expanding-MHD delta-f feedback remains rejected.
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
- After AMR child reconstruction, `Particles::UpdateAfterAMR` refreshes the
  retained pack pointer and validates particle-owned MeshBlock array capacity.

### Runtime guards (constructor)
- Coupling requires `deposit_moments=true` and an active `<mhd>` block.
- Coupled mode is rejected for `radiation+MHD`, hydro/ion-neutral, and
  numerical-relativity (`adm`/`z4c`) compositions.
- Energy feedback requires ideal MHD EOS; the exact-isothermal VL2 paper paths
  are explicitly guarded momentum-only exceptions with energy feedback off.
- `edge_staggered` and fluid feedback branches are restricted to
  non-relativistic MHD.
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

Boundary notes:
- Ordinary pushers handle `BoundaryFlag::reflect` before deposition or exchange
  by mirroring position and flipping the normal state component. VL2 deposits
  its predictor or realized impulse first, then reflects in the staged
  half-drift before exchange.
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
  `AdaptDeltaF -> SaveOldPositions -> Push -> ZeroMoments -> InitRecvMoments ->
  DepositMoments -> RestrictMoments -> SendMoments -> RecvMoments ->
  ClearRecvMoments -> ClearSendMoments -> ApplyMomentPhysicalBCs ->
  ProlongateMoments`. VL2 wrappers are stage-gated, so this stage-0 chain
  does not perform its coupled update.

Particle migration communication depends on coupling mode:
- uncoupled/default: `NewGID -> SendCnt -> InitRecv -> SendP -> RecvP ->
  ClearRecv -> ClearSend` in `before_timeintegrator`
- coupled non-VL2 runs: the same migration chain moves to
  `after_timeintegrator`
- VL2: the chain runs in `after_stagen`, so ownership is corrected after
  both the midpoint (stage 1) and endpoint (stage 2) half drifts.

In coupled mode, moment wrappers are also inserted into `stagen`:
- non-VL2 runs retain the stage-1 insertion anchor selected by
  `couple_fluid_feedback_order` (`MHD::MHDSrcTerms` vs `MHD::EFieldSrc`)
- full-f VL2 inserts its staged chain after `MHD::CopyCons` and before
  `MHD::Fluxes`: optional full-Hall predictor deposition/`v_H` construction ->
  `Push` -> `SaveOldPositions` -> moment deposition/synchronization ->
  staged half drift. Hall-off runs the generic moment wrappers in both stages.
  Full Hall instead uses separately synchronized `cr_hall_moments` for its
  stage-1 analytic source and runs the generic `moments` wrappers only in stage
  2, where they deposit the realized particle `dp/dt` and `dE/dt` used by the
  exact gas and grid correctors.
- if `couple_j_to_efield_representation=edge_staggered` and
  `couple_j_deposition_mode=cc_convert`, `ConvertCoupledCurrentRepresentation`
  is inserted before `MHD::EFieldSrc` to populate the retained edge-current
  representation; current `EFieldSrc` does not turn that array into an EMF.
- if `couple_j_deposition_mode=direct_staggered`, direct edge-current
  synchronization and physical-BC tasks are inserted before `MHD::EFieldSrc`
  instead of the CC conversion task. These arrays are retained engineering and
  diagnostic state, not the full-Hall CT route. Physical edge-current BCs cover
  `periodic`, `reflect`, and `outflow`; direct mode rejects `inflow`, so use
  `couple_j_deposition_mode=cc_convert` for inflow-boundary coupled runs.

---

## Outputs and Diagnostics (references)
- `outputs/vtk_prtcl.cpp` and `outputs/track_prtcl.cpp` consume `prtcl_rdata` and
  `prtcl_idata` for particle dumps and tracked particles.
- Particle VTK output emits explicit `gid`, `ptag`, `species`, `cr_source`,
  `macro_weight`, `birth_time`, `deltaf_f0`, and `deltaf_weight` scalar
  diagnostics.
- `trk` output selects particles by nonnegative `PTAG < <output>/nparticles` and
  writes big-endian float32 position/velocity rows at tag-derived offsets.
- Derived variable `prtcl_d` (in `outputs/derived_variables.cpp`) bins particle
  counts onto the mesh.
- Final Q-017 stdout telemetry reports wrapper-boundary `adaptive_deltaf`,
  `push`, `deposition`, and `migration` elapsed timers, fixed-record resident
  bytes by species and absolute logical mesh-refinement level, a direct
  particle-view allocated-byte snapshot, and tracked AthenaK-owned Kokkos-view
  final snapshots plus rank-local high-water marks. Tracked high-water samples
  include particle migration lists and counters while live, moment-boundary
  helpers, and visible MeshRefinement helpers. They include view padding but
  exclude allocator overhead, runtime caching, untracked subsystems and any
  temporally concurrent aggregate multi-rank peak; they are not total
  GPU-memory measurements.
- Set `<particles>/pic_q017_sync_kernel_timers=true` only for dedicated timing
  runs. It fences measured particle boundaries so device elapsed times are
  interpretable, but the synchronization intentionally perturbs execution.

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
