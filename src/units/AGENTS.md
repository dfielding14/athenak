# AGENTS.md

## Purpose
This directory defines the `units::Units` class that manages code-unit <-> cgs
conversions and exposes physical constants. It is constructed early in
`MeshBlockPack::AddPhysics` when a `<units>` block exists and is used by
multiple physics modules for opacity, cooling, gravity, and scale conversions.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files

### `units.hpp`
- Declares `units::Units`.
- Defines static `constexpr` cgs constants and conversion helpers.

### `units.cpp`
- Implements `Units` constructor/destructor and all conversion functions.
- Implements the GR-specific override for length/mass/time scales.

---

## Construction and Lifetime

### Creation
`MeshBlockPack::AddPhysics` creates `punit` **only** when a `<units>` block
exists; otherwise `punit` is `nullptr`. (`src/mesh/meshblock_pack.cpp`)

### Constructor inputs (`<units>` block)
The constructor reads:
- `length_cgs` (default `1.0`)
- `mass_cgs` (default `1.0`)
- `time_cgs` (default `1.0`)
- `mu` (mean molecular weight, default `1.0`)

### GR override
If `coord/general_rel=true`, the constructor **overrides** the MLT scales using:
- `density_cgs` (required)
- `bhmass_msun` (required, in solar masses)

Computed overrides:
- `length_cgs_ = G * M_bh / c^2`
- `mass_cgs_   = density_cgs * length_cgs_^3`
- `time_cgs_   = length_cgs_ / c`

These values replace any `length_cgs`/`mass_cgs`/`time_cgs` values passed in the
`<units>` block.

---

## Stored Scales and Conventions

### Base scales (cgs per code unit)
`length_cgs_`, `mass_cgs_`, `time_cgs_` are stored as **cgs per code unit**.
Thus:
- **code -> cgs**: multiply by `*_cgs()`
- **cgs -> code**: multiply by the reciprocal provided by `cm()`, `g()`, etc.

### Derived scales (cgs per code unit)
All derived scale functions return **cgs per code unit**:
- `velocity_cgs = length_cgs / time_cgs`
- `density_cgs  = mass_cgs / length_cgs^3`
- `energy_cgs   = mass_cgs * velocity_cgs^2`
- `pressure_cgs = energy_cgs / length_cgs^3`
- `temperature_cgs = velocity_cgs^2 * mu * atomic_mass_unit_cgs / k_boltzmann_cgs`

---

## Conversion Helpers (cgs -> code)

These return **code units per cgs unit** (divide by the scale above):
- Length: `cm()`, `pc()`, `kpc()`
- Mass: `g()`, `msun()`, `atomic_mass_unit()`
- Time: `s()`, `yr()`, `myr()`
- Velocity: `cm_s()`, `km_s()`
- Density: `g_cm3()`
- Energy: `erg()`
- Pressure: `dyne_cm2()`
- Temperature: `kelvin()`

---

## Physical Constants (cgs and code units)

### Static cgs constants (from `units.hpp`)
Available as `Units::<name>_cgs`:
- Base units: `cm_cgs`, `pc_cgs`, `kpc_cgs`, `g_cgs`, `msun_cgs`,
  `atomic_mass_unit_cgs`, `s_cgs`, `yr_cgs`, `myr_cgs`,
  `cm_s_cgs`, `km_s_cgs`, `g_cm3_cgs`, `erg_cgs`, `dyne_cm2_cgs`, `kelvin_cgs`
- Physical constants: `k_boltzmann_cgs`, `grav_constant_cgs`,
  `speed_of_light_cgs`, `rad_constant_cgs`,
  `electron_rest_mass_energy_cgs`
- Opacity coefficients: `rosseland_coef_cgs`,
  `planck_minus_rosseland_coef_cgs`

### Code-unit constants (computed)
Provided as member functions:
- `k_boltzmann()`:
  `k_boltzmann_cgs / (energy_cgs / temperature_cgs)`
- `grav_constant()`:
  `grav_constant_cgs * density_cgs * time_cgs^2`
- `speed_of_light()`:
  `speed_of_light_cgs / velocity_cgs`

---

## Direct Usage in the Codebase

Files with direct references to `punit` or `Units`:
- `src/mesh/meshblock_pack.cpp` (construction)
- `src/mesh/meshblock_pack.hpp` (pointer declaration)
- `src/mesh/mesh_refinement.cpp` (preserve `punit` across AMR rebuild)
- `src/radiation/radiation.cpp` (compute `arad` when units enabled)
- `src/radiation/radiation_source.cpp`
  (density/temperature scales, opacity coefficients, electron energy)
- `src/diffusion/conduction.cpp` (temperature and kappa scaling)
- `src/srcterms/srcterms.cpp`, `src/srcterms/srcterms_newdt.cpp`
  (cooling/heating scaling)
- `src/particles/particles.cpp`, `src/particles/particles_pushers.cpp`
  (gravity constant and time scaling)
- `src/pgen/cgm_*.cpp`, `src/pgen/cgm_static.cpp`
  (gravity constant and unit conversions)

---

## Cautions

- `punit` is `nullptr` if the `<units>` block is omitted. Any module that uses
  `punit` must be guarded by configuration checks (some features, e.g. Compton
  in radiation, explicitly require units).
- In GR mode (`coord/general_rel=true`), `density_cgs` and `bhmass_msun` are
  mandatory because the constructor calls `GetReal` (no default).
