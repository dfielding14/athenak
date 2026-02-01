# AGENTS.md

## Purpose
This directory provides source-term infrastructure and turbulence driving used by
Hydro, MHD, and Radiation modules. It includes cooling/heating models, shearing-box
terms, beam sources, and the Ornstein-Uhlenbeck turbulence driver.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### SourceTerms
- `srcterms.hpp` / `srcterms.cpp`: `SourceTerms` class, parameter parsing, and source
  term implementations for hydro/MHD/radiation.
- `srcterms_newdt.cpp`: computes timestep constraints from active cooling terms.
- `ismcooling.hpp`: `ISMCoolFn` cooling curve (tabulated + fits).
- `cooling_tables.hpp`: CGM cooling tables (large static arrays, loaded to device).

### Turbulence driving
- `turb_driver.hpp` / `turb_driver.cpp`: `TurbulenceDriver` class, mode selection,
  Ornstein-Uhlenbeck evolution, and forcing application.
- `README.md`: turbulence driver overview (verify against code before relying on it).

---

## How SourceTerms Are Used
`SourceTerms` instances are created by physics modules using their own input blocks:
- Hydro: `new SourceTerms("hydro", ...)`
- MHD: `new SourceTerms("mhd", ...)`
- Radiation: `new SourceTerms("radiation", ...)` (beam source support)

The same class handles multiple term types; only enabled terms are applied.

---

## SourceTerms Features and Parameters

Parameters are read from the block passed at construction (`<hydro>`, `<mhd>`,
or `<radiation>`), except where noted.

### Constant acceleration
- `const_accel` (bool)
- `const_accel_val`, `const_accel_dir` (1..3)

### ISM cooling (optically thin)
- `ism_cooling` (bool)
- `hrate` (heating rate)
- `cooling_dt_factor` (multiplies dt constraint)
- `t_start_ism_cooling` (delay start time)

### CGM cooling
- `cgm_cooling` (bool) -> incompatible with `ism_cooling`
- Heating profile parameters: `hrate`, `hscale_norm`, `hscale_height`,
  `hscale_radius`, `hscale_alpha`
- `T_max` temperature ceiling
- Uses tables from `cooling_tables.hpp` loaded in `Initialize()`

### Relativistic cooling
- `rel_cooling` (bool)
- `crate_rel`, `cpower_rel`

### Beam source (radiation)
- `beam_source` (bool)
- `dii_dt` (intensity injection rate)

### Shearing box (block `<shearing_box>`)
- `qshear`, `omega0`
- Source term functions have hydro and MHD overloads, plus `SBoxEField` for 2D.

### Units dependency
Cooling/heating terms use `pmy_pack->punit` conversions. Enable `<units>` when using
ISM/CGM cooling or radiation terms that rely on physical units.

---

## TurbulenceDriver Overview
The turbulence driver is not a `SourceTerms` subclass. It is instantiated when the
`<turb_driving>` block is present (via `MeshBlockPack::AddPhysics`) and injects a
stochastic forcing field.

Key behaviors (see `turb_driver.cpp`):
- Mode selection between `nlow` and `nhigh`, with optional `npeak`/`kpeak` shaping.
- Ornstein-Uhlenbeck evolution with `tcorr`, `dt_turb_update`.
- Solenoidal/compressive mix via `sol_fraction`.
- Optional spatial windowing (`x_turb_scale_height`, centers, etc.).
- Optional tiled driving (`tile_driving`, `tile_nx/ny/nz` or `tile_factor`).
- `turb_flag`: 1 = decaying, 2 = continuously driven.

Tasks:
- `IncludeInitializeModesTask` inserts mode setup and forcing update into
  `before_timeintegrator`.
- `IncludeAddForcingTask` inserts forcing between RK update and source terms.

### Initial turbulence kick
If an `<initial_turb>` block exists, `main.cpp` creates a temporary
`TurbulenceDriver` to apply a one-time impulse at startup.

---

## Timestep Constraints
`SourceTerms::NewTimeStep` computes dt limits from active cooling terms and applies
`cooling_dt_factor` when any cooling is enabled. Hydro/MHD call this each cycle.

---

## Extension Points
- Add a new source term in `srcterms.cpp` and wire it into `SourceTerms::NewTimeStep`.
- Add or modify cooling tables in `cooling_tables.hpp` (be mindful of size).
- Extend `TurbulenceDriver` by adding new forcing spectra or spatial windows and
  updating task insertion logic as needed.
