# Module: Particles

## Overview
The Particles module implements Lagrangian particle tracking with pusher algorithms,
interpolation schemes, and particle-mesh coupling. The module supports two
**mutually exclusive** particle types: star particles and cosmic ray particles.

The coupled runtime interface and its qualification boundary are specified in
[MHD-PIC Runtime Model Contract](../engineering/pic_mhd_model_contract.md).
Particle state, integration chronology, coupling, and optional physics are
selected by explicit controls rather than a bundled physical-mode selector.

## Source Location
`src/particles/`

## Key Components

| File | Purpose | Key Functions |
|------|---------|---------------|
| `particles.hpp/cpp` | Core particle class | Initialization, type selection, memory management |
| `particles_pushers.cpp` | Integration algorithms | `rk4_gravity()`, `boris_lin()`, `boris_tsc()` |
| `particles_tasks.cpp` | Task registration | Particle evolution tasks |
| `bvals/bvals_part.cpp` | Boundary communication | MPI particle exchange |

## Particle Types (Mutually Exclusive)

### Star Particles
- **Purpose**: Represent stellar objects in gravitational potentials
- **Memory**: 9 real data slots per particle
- **Compatible Pushers**: `rk4_gravity`, `drift`
- **Incompatible**: Boris pushers (will cause fatal error)

### Cosmic Ray Particles
- **Purpose**: Charged particles in magnetic fields
- **Memory**: 26 real data slots per particle, including sampled fields,
  displacement/feedback diagnostics, per-particle macro weight, delta-f state,
  and birth time
- **Compatible Pushers**: `drift`, `boris_lin`, `boris_tsc`

## Configuration Parameters

### Common Parameters (`<particles>` block)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `particle_type` | string | Yes | Either `star` or `cosmic_ray` |
| `pusher` | string | Yes | Integration algorithm (see below) |
| `pic_interp_scheme` | string | No | PIC interpolation scheme; currently only `tsc` is supported |

For cosmic rays, `boris_lin` and `boris_tsc` always store mass-normalized
momentum `p/m`; `drift` stores velocity. `<particles>/pic_cr_initial_state`
controls only whether the initializer supplies velocity (converted to `p/m`
for a Boris pusher) or already supplies momentum. It does not select the stored
state representation.

`<time>/integrator=vl2` selects the Sun--Bai staged predictor-corrector
chronology. It requires cosmic-ray `boris_tsc`, coupled MHD background and
feedback, moment deposition with `deposit_order=2`, conservative momentum
feedback, the qualified energy-feedback composition, and
`couple_fluid_feedback_order=mhd_src_terms`. `<time>/integrator=rk2` remains
Heun and does not imply staged PIC coupling.

### Star Particle Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `star_particle_file` | string | Yes | ASCII file with initial positions |
| `grav_dx` | Real | No | Finite difference spacing for gravity (default: 1e-6) |

Star-particle runs require a `<units>` block because the analytic potential
converts code units through the configured unit system.

**Required `<potential>` block parameters:**
```ini
r_scale     # Radial scale
rho_scale   # Density scale
mass_gal    # Galaxy mass
scale_gal   # Galaxy scale length
z_gal       # Vertical scale
r_200       # Virial radius
rho_mean    # Mean density
```

### Cosmic Ray Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nspecies` | int | 1 | Number of species |
| `ppc` | Real | 1.0 | Particles per cell; must be non-negative; fractional counts are distributed by global MeshBlock ID |
| `cr_distribution` | string | center | `center` places particles on deterministic cell centers; `random` uses deterministic per-block pseudorandom positions |
| `assign_tag` | string | index_order | Tagging method |
| `deposit_moments` | bool | false | Deposit particle charge/current moments |
| `deposit_qscale` | Real | 1.0 | Root-grid macro-charge and macro-mass scale |
| `pic_background_mode` | string | coupled | `coupled`, `passive_mhd`, or `no_mhd` |
| `pic_feedback_mode` | string | background-dependent | `coupled` or `test_particle`; defaults to `coupled` for a coupled background and `test_particle` otherwise |
| `pic_enable_2d3v` | bool | false | Required for Boris pushers on 2D meshes; keeps `vz` and `Bz` active when `nx3=1` |
| `pic_cr_light_speed` | Real | 1.0 | Positive artificial CR light speed used by momentum-state cosmic-ray Boris pushers |
| `pic_cr_initial_state` | string | velocity | Input interpretation only: `velocity` converts initialized values to stored `p/m` for Boris pushers; `momentum` treats them as already `p/m` |
| `pic_cr_hall_mode` | string | off | `off` or `full`; `full` selects the atomic CR-Hall closure and requires `time/integrator=vl2`, ideal MHD, uniform-grid full-f coupling, and positive `pic_background_ion_q_over_mc` |
| `pic_background_ion_q_over_mc` | Real | 0.0 | Background-ion charge-to-mass normalization; it must be positive when `pic_cr_hall_mode=full`, fixes the physical electron-charge denominator, and is not a Hall-strength knob |
| `pic_wave_damping_mode` | string | off | `off` or `ion_neutral_friction`; the reduced static-neutral transverse-friction map requires coupled Newtonian single-fluid MHD |
| `pic_ion_neutral_collision_rate` | Real | 0.0 | Non-negative `nu_in`; must be positive when reduced ion-neutral friction is enabled |
| `pic_max_cell_cross` | int | 2 | Particle cell-crossing timestep limit before global CFL scaling; must not exceed the smallest active MeshBlock dimension |
| `pic_theta_max` | Real | 0.3 | Boris gyro-angle timestep limit before global CFL scaling |
| `pic_random_seed` | int | 0 | Seed for deterministic `cr_distribution=random` placement |
| `pic_load_balance_cost_per_particle` | Real | 0.0 | Non-negative opt-in particle contribution to each post-AMR MeshBlock load cost |
| `pic_deltaf_mode` | string | off | `off`, `quiet_start`, or `physical`; `quiet_start` changes sampling only, while `physical` evolves and deposits delta-f weights |
| `pic_deltaf_f0` | string | empty | Required when delta-f is enabled; accepted backgrounds are `uniform`, `uniform_quiet`, `kappa_iso`, `kappa_drift`, and `kappa_aniso` |
| `pic_deltaf_p0`, `pic_deltaf_kappa` | Real | 1.0, 1.25 | Positive momentum scale and kappa index for analytic kappa backgrounds |
| `pic_deltaf_drift_x1/x2/x3` | Real | 0.0 | Analytic background drift components |
| `pic_deltaf_aniso_x1/x2/x3` | Real | 1.0 | Positive analytic background anisotropy scales |
| `pic_deltaf_background_rho` | Real | 0.0 | Non-negative analytic background charge density used by physical delta-f gas feedback |
| `pic_deltaf_background_jx/jy/jz` | Real | 0.0 | Analytic background current used by physical delta-f gas feedback |
| `pic_deltaf_adapt_mode` | string | off | `off` or `global_bikappa_moments_experimental`; the adaptive x1-parallel global fit is restricted to its bounded physical-delta-f expanding-box configuration |
| `pic_deltaf_adapt_interval` | Real | 0.0 | Non-negative adaptive-fit cadence; must be positive when adaptive delta-f is enabled |
| `pic_expanding_box_mode` | string | off | `off` or `on`; enables particle and MHD expanding-box transformations |
| `pic_expansion_law` | string | linear | `linear`, `reciprocal_linear`, or `exponential` scale-factor law |
| `pic_expansion_rate_x1/x2/x3` | Real | 0.0 | Directional rates; enabled laws must remain finite and positive through `<time>/tlim` |

Each `<speciesN>` block defines at least `mass` and `charge`. Species masses
must be positive.

When `couple_j_deposition_mode=direct_staggered`, the deposited edge current is
constructed from each particle trajectory over the step. The implementation
requires the old-to-new shape support to shift by no more than one cell in each
active dimension and aborts if a step is too large. Reduce `time/cfl_number` or
`pic_max_cell_cross` if that guard is triggered. Direct edge-current physical
BCs support periodic, reflecting, and outflow boundaries; use
`couple_j_deposition_mode=cc_convert` for inflow-boundary coupled runs.

## Particle Data Arrays

### Index Definitions
```cpp
enum ParticlesIndex {
  // Integer indices
  PGID=0,  // Global MeshBlock ID
  PTAG=1,  // Particle tag
  PSP=2,   // Species index
  PCRSOURCE=3, // Persistent CR source cohort

  // Real data indices
  IPX=0, IPVX=1,  // X position, velocity or p/m component
  IPY=2, IPVY=3,  // Y position, velocity or p/m component
  IPZ=4, IPVZ=5,  // Z position, velocity or p/m component

  // Cosmic ray only (indices 6-25)
  IPM=6,          // Charge/mass ratio
  IPBX=7, IPBY=8, IPBZ=9,  // B-field components
  IPDX=10, IPDY=11, IPDZ=12, IPDB=13, // Displacements
  IPEX=14, IPEY=15, IPEZ=16,          // Electric-field samples
  IPDPX=17, IPDPY=18, IPDPZ=19,       // Feedback momentum rate
  IPDE=20, IPEBDOT=21,                // Feedback energy diagnostics
  IPWT=22,                            // Relative macro-particle weight
  IPF0=23,                            // Initial analytic delta-f background value
  IPDFWT=24,                          // Evolving delta-f perturbation weight
  IPT_BIRTH=25                        // Persistent creation time
};
```

The shared `IPVX/IPVY/IPVZ` names are retained for restart-layout
compatibility. They store mass-normalized momentum `p/m` for cosmic-ray Boris
pushers and velocity otherwise. Current deposition, timestep checks, VTK
output, and tracked-particle output derive physical velocity from that
pusher-selected representation.

Cosmic-ray particles created from `ppc` carry `IPWT = cell_volume/root_cell_volume`
so refined-cell particles keep the same physical density represented by a
root-level `deposit_qscale`. Manually injected particles may set `IPWT=1.0`
when `deposit_qscale` is already the intended macro-particle mass/charge.
Every cosmic-ray particle also carries a persistent integer `PCRSOURCE` cohort
(`initial` or `shock_injected`) and real `IPT_BIRTH` creation time.

`pic_deltaf_mode=quiet_start` selects deterministic low-discrepancy placement
for `cr_distribution=random` tests without changing deposited moments. With
`pic_deltaf_mode=physical`, `IPF0` stores the initial analytic background,
`IPDFWT` evolves as `1 - f0(t,x,p)/f(0,x0,p0)`, and moment deposition uses the
perturbation weight. Particle VTK output exposes named `gid`, `ptag`, `species`,
`cr_source`, `macro_weight`, `birth_time`, `deltaf_f0`, and `deltaf_weight`
scalars. The bounded Q-016 helper resolves weighted physical-speed spectra by
species, source, and birth-time cohort. Its `full_f` mode accumulates
`macro_weight`; its explicit `delta_f_perturbation` mode accumulates the signed
sampled perturbation `macro_weight * deltaf_weight` without adding the analytic
background represented by `deltaf_f0`. See
{doc}`../engineering/pic_q016_particle_provenance_spectra` for the bounded
workflow and qualification limits.

`pic_deltaf_adapt_mode=global_bikappa_moments_experimental` inserts a global
weighted bi-kappa fit before particle pushing when a new cadence bucket begins.
The fitted `xi`, fitted `p0`, and cadence bucket persist in restart schema
version 8. This host-verified extension is a bounded mechanics path, not a
qualified CRPAI transport calibration.

During adaptive refinement, particle migration resolves the post-AMR owning
MeshBlock geometrically before rank reassignment completes. The retained
particle module then validates its MeshBlock-sized arrays against the rebuilt
pack. `pic_load_balance_cost_per_particle > 0` adds a particle-count term to
the base fluid cost of each post-AMR MeshBlock before load balancing.
The retained-state inventory and the target `paper_smooth`
receiver-resolution TSC refinement-interface deposition policy are documented
in {doc}`../engineering/pic_amr_lifetime_and_interface_policy`. The generic
restriction/exchange/prolongation path remains supporting infrastructure until
the fine/coarse oracle passes. The existing direct-staggered trajectory-current
path is an experimental candidate, not a qualified conservative AMR
gas-feedback policy.

In expanding-box mode with active MHD, raw face-centered arrays are
divergence-preserving comoving magnetic fluxes. MHD consumers derive physical
face fields and CT maps edge EMFs before updating those fluxes. This bounded
path rejects unqualified compositions, including AMR, non-periodic
boundaries, relativistic coordinates, nonideal MHD, fluid/radiation/relativity
or turbulence blocks, staggered/direct current deposition, nonzero analytic
delta-f background current, the full CR-Hall closure, and user-defined history
callbacks. The qualified coupled source split uses cell-centered `cc_convert`
moments, conservative momentum and energy feedback, `mhd_src_terms` ordering,
final-frame particle EM impulses, and physical-volume moment normalization.
The admitted adaptive physical delta-f extension uses a separate
endpoint-normalized analytic `rho E + J x B` source path and is not an
opposite-impulse conservation claim. Built-in MHD `hst` output uses physical
cell volumes and physical magnetic fields. The reduced static-neutral damping
map is permitted only with one of those admitted cell-centered source splits.

## Pusher Algorithms

### Compatibility Matrix

| Pusher | Star Particles | Cosmic Rays | Description |
|--------|---------------|-------------|-------------|
| `drift` | yes | yes | Simple drift motion |
| `rk4_gravity` | yes | no | RK4 in gravitational potential |
| `boris_lin` | no | yes | Boris pusher with linear field interpolation |
| `boris_tsc` | no | yes | Boris pusher with TSC field interpolation |

### RK4 Gravity Pusher (Stars)
4th-order Runge-Kutta integration in gravitational potential:
```cpp
// Evaluates gravitational acceleration via finite differences
// Uses potential parameters from <potential> block
```

### Boris Pushers (Cosmic Rays)
Relativistic particle motion in electromagnetic fields:
- `boris_lin`: Linear interpolation of fields
- `boris_tsc`: Triangular-shaped cloud (TSC) interpolation

## Usage Examples

### Star Particle Simulation
```ini
<job>
problem_id = stellar_dynamics

<particles>
particle_type = star
pusher = rk4_gravity
star_particle_file = initial_stars.txt
grav_dx = 1.0e-6

<units>
length_cgs = 1.0
mass_cgs   = 1.0
time_cgs   = 1.0
mu         = 1.0

<potential>
r_scale = 8.0      # kpc
rho_scale = 1.0e7  # Msun/kpc^3
mass_gal = 1.0e12  # Msun
scale_gal = 3.5    # kpc
z_gal = 0.5        # kpc
r_200 = 200.0      # kpc
rho_mean = 1.0e6   # Msun/kpc^3
```

### Cosmic Ray Simulation
```ini
<job>
problem_id = cosmic_ray_propagation

<particles>
particle_type = cosmic_ray
pusher = boris_tsc
pic_interp_scheme = tsc
nspecies = 1
ppc = 8.0
pic_cr_light_speed = 1.0
pic_max_cell_cross = 2
pic_theta_max = 0.3

<species0>
mass = 1.0
charge = 1.0
```

## Star Particle File Format
ASCII file with 8 columns (3D only):
```
# x y z vx vy vz t_create mass
0.1 0.2 0.3 1.0 0.5 0.0 0.0 1.0e6
0.2 0.3 0.4 0.8 0.6 0.1 1.0 2.0e6
```

## Safety Features

### Type-Specific Memory Allocation
- Star particles: Exactly 9 slots (no waste)
- Cosmic rays: 26 slots, including field samples, feedback diagnostics,
  relative macro-particle weight, delta-f state, and birth time

### Runtime Protection
1. **Initialization check**: Fatal error if incompatible pusher selected
2. **Pusher verification**: Boris pushers require cosmic-ray particles; RK4
   gravity requires star particles
3. **Clear error messages**: Explains exactly what went wrong

Example error:
```
### FATAL ERROR in particles.cpp at line 198
Boris pushers are incompatible with star particles; use
<particles>/particle_type=cosmic_ray for boris_lin/boris_tsc, or use
drift/rk4_gravity for star particles.
```

## Implementation Notes

### Mutual Exclusivity
- Only ONE particle type per simulation
- Types cannot be mixed in same run
- Particle type cannot change across a restart
- This is by design for safety and efficiency

### Memory Layout
- Real data: `prtcl_rdata[nrdata][nparticles]`
- Integer data: `prtcl_idata[nidata][nparticles]`
- Kokkos Views for GPU portability
- PIC restart files store particle records, `nrdata`, and `nidata` for every
  active particle mode. Runs with deposited moments also persist the moment
  arrays, and edge-current coupling persists its staggered current arrays.
  Changing the particle record layout intentionally rejects older incompatible
  PIC restart files. Schema version 8 fingerprints continuation-sensitive PIC
  selectors and coefficients, adaptive delta-f state, and star-potential
  parameters so incompatible continuation attempts fail before evolution
  resumes.

### Boundary Communication
- The MPI exchange implementation migrates particles between MeshBlocks and
  ranks.
- The AMR ownership-refresh path geometrically reconstructs post-refinement
  owners before rank reassignment completes.
- Qualification remains open for repeated MPI migration/restart, forced
  empty-rank load balance, refine/derefine lifetime, memory checking, and GPU
  execution.

## Limitations

1. Star particles require 3D simulations
2. Particle types cannot be mixed
3. No dynamic type switching
4. Boris pushers only work with cosmic rays

## Performance Considerations

- Star particles: 72 bytes per particle (9x8)
- Cosmic rays: 208 bytes per particle for real data (26x8), plus integer data
- All pushers use Kokkos parallel_for for GPU acceleration

## References

- Stone et al. (2020) for Athena++ particle framework
