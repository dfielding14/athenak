# AGENTS.md

## Purpose
This directory implements the Z4c numerical relativity module: evolution of the
Z4c variables, ADM conversions, constraint monitoring, gauge handling, optional
Sommerfeld boundary terms, AMR drivers, Weyl scalar extraction, and compact
object tracking. It also defines the `Tmunu` stress-energy container used to
couple matter to the Z4c RHS.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core class and data layout
- `z4c.hpp`: `z4c::Z4c` definition, variable indices/names, data arrays, task
  entry points, and options.
- `z4c.cpp`: constructor/destructor, option parsing, allocations, spherical
  grids for wave extraction, and algebraic constraint enforcement.

### Evolution and gauge
- `z4c_calcrhs.cpp`: RHS evaluation for the Z4c system, with matter coupling
  (via `Tmunu`), gauge terms, constraint damping, and KO dissipation.
- `z4c_update.cpp`: explicit RK update (`u0` from `u1` and `u_rhs`).
- `z4c_newdt.cpp`: time step estimate (min cell size).
- `z4c_gauge.cpp`: pre-collapsed lapse initialization from ADM `psi4`.

### ADM conversions and constraints
- `z4c_adm.cpp`: `ADMToZ4c`, `Z4cToADM`, and ADM constraint diagnostics.

### Boundary conditions and task scheduling
- `z4c_tasks.cpp`: Z4c task queue for the NumericalRelativity task list.
- `z4c_Sbc.cpp`: Sommerfeld boundary RHS terms on outflow/diode (and optional
  user) boundaries.

### AMR and diagnostics
- `z4c_amr.hpp` / `z4c_amr.cpp`: Z4c-specific refinement drivers.
- `z4c_calculate_weyl_scalars.cpp`: computation of Psi4 (real/imag) from ADM
  data into `u_weyl`.
- `z4c_wave_extr.cpp`: spherical interpolation, spin-weighted harmonic
  decomposition, and waveform output.
- `compact_object_tracker.hpp` / `compact_object_tracker.cpp`: tracker for
  puncture motion using shift (and optional fluid velocity for NS).

### Matter coupling
- `tmunu.hpp` / `tmunu.cpp`: `Tmunu` container for `E`, `S_i`, `S_ij` and storage
  in `u_tmunu`.

---

## Data Layout and Indices

### Z4c evolved variables (`z4c::Z4c::nz4c`)
Indices are defined in `z4c.hpp` and names in `Z4c_names`:
- `chi`
- conformal metric `g_dd` (`gxx`, `gxy`, `gxz`, `gyy`, `gyz`, `gzz`)
- `Khat`
- traceless extrinsic curvature `A_dd` (`Axx`, `Axy`, `Axz`, `Ayy`, `Ayz`, `Azz`)
- conformal connection `Gam_u` (`Gamx`, `Gamy`, `Gamz`)
- `Theta` (zeroed when `use_z4c=false` for BSSN)
- lapse `alpha`
- shift `beta_u` (`betax`, `betay`, `betaz`)

### Constraint variables (`z4c::Z4c::ncon`)
Stored in `u_con` and exposed via `con.*`:
- `C` (constraint monitor)
- `H` (Hamiltonian constraint)
- `M` (momentum norm)
- `Z` (Z-constraint norm)
- `Mx`, `My`, `Mz` (momentum components)

### Weyl scalars
- `u_weyl` holds two components: real and imaginary Psi4 (`weyl.rpsi4`,
  `weyl.ipsi4`).

### Matter stress-energy (`Tmunu`)
`u_tmunu` stores:
- `S_dd` (stress tensor, symmetric)
- `E` (energy density)
- `S_d` (momentum density)

---

## Task Flow (NumericalRelativity list)

`Z4c::QueueZ4cTasks` enqueues the following (see `z4c_tasks.cpp`):

**Start**
- `InitRecv` (post boundary receives for `u0`)
- `InitRecvWeyl` (schedule Weyl exchange when waveform output is due)

**Run**
- `CopyU` (stage bookkeeping for RK)
- `CalcRHS<NGHOST>` (NGHOST=2/3/4)
  - Optional dependency on `MHD_SetTmunu` when matter is present.
- `Z4cBoundaryRHS` (Sommerfeld RHS on selected boundaries)
- `ExpRKUpdate` (update `u0` from `u1` and `u_rhs`)
- `RestrictU` (AMR/SMR restriction)
- `SendU` / `RecvU` (boundary exchange)
- `ApplyPhysicalBCs` (calls `MeshBoundaryValues::Z4cBCs`; user BCs after)
- `Prolongate` (AMR/SMR prolongation)
- `EnforceAlgConstr` (det(g)=1, trace(A)=0) when `pdyngr` or last RK stage
- `ConvertZ4cToADM` (update ADM variables) when `pdyngr` or last RK stage
- `UpdateExcisionMasks` (if BH excision enabled in coordinates)
- `NewTimeStep` (last RK stage only; min cell size)

**End**
- `ClearSend` / `ClearRecv`
- `ADMConstraints_` (computes `u_con` on last stage)
- Weyl pipeline (only if wave extraction enabled and scheduled):
  `CalcWeylScalar` -> `RestrictWeyl` -> `SendWeyl` -> `RecvWeyl` ->
  `ProlongateWeyl` -> `ClearSendWeyl` -> `ClearRecvWeyl` ->
  `CalcWaveForm`
- `TrackCompactObjects` (if trackers configured)

---

## RHS and Gauge Details (from `z4c_calcrhs.cpp`)

- Uses finite-difference operators from `utils/finite_diff.hpp` (`Dx`, `Dxx`,
  `Dxy`, `Diss`) with template `NGHOST` order.
- Matter coupling is included if `ptmunu` is present:
  - Adds `E`, `S_i`, and `S_ij` contributions to RHS terms.
- Gauge parameters are read from `<z4c>` and used as:
  - Lapse: `rhs.alpha = lapse_advect * Lalpha - f * alpha * Khat`, with
    `f = lapse_oplog*lapse_harmonicf + lapse_harmonic*alpha`.
  - Shift: `rhs.beta_u = shift_ggamma * Gam + shift_advect * Lbeta - shift_eta * beta`
    plus optional harmonic terms (`shift_alpha2Gamma`, `shift_H`).
- `chi_div_floor` bounds `chi` in divisions (`chi_guarded`).
- `use_z4c=false` multiplies `rhs.Theta` by zero (BSSN mode).
- KO dissipation is added via `Diss<NGHOST>` and the precomputed `diss` factor.

---

## ADM Conversion and Constraints

### `ADMToZ4c<NGHOST>` (`z4c_adm.cpp`)
- Computes conformal metric, `chi`, `Khat`, and traceless `A_dd` from ADM
  `g_dd` and `K_dd`.
- Computes `Gam_u` via derivatives of the inverse conformal metric.
- Calls `AlgConstr` to enforce det(g)=1 and trace(A)=0.

### `Z4cToADM`
- Builds `adm.psi4`, `adm.g_dd`, and `adm.vK_dd` from Z4c variables.
- Uses `(Khat + 2*Theta)` in the ADM extrinsic curvature.

### `ADMConstraints<NGHOST>`
- Computes Hamiltonian, momentum, Z-constraint, and monitor `C^2` on interior
  zones only.
- Includes matter terms via `Tmunu` when present.

---

## Boundary Conditions

### Physical BCs (in `bvals/physics/z4c_bcs.cpp`)
`MeshBoundaryValues::Z4cBCs` applies:
- Reflecting, outflow/diode/vacuum (one-sided extrapolation), inflow, and
  periodic/shear-periodic handling.
- Extrapolation order set by `opt.extrap_order` (clamped to 2..4 and <= nghost).

### Sommerfeld RHS (`z4c_Sbc.cpp`)
- `Z4cBoundaryRHS` applies Sommerfeld-type RHS corrections to `Theta`, `Khat`,
  `Gam_u`, and `A_dd` on outflow/diode boundaries.
- On `BoundaryFlag::user`, Sommerfeld is applied only if `z4c/user_Sbc=true`.

### User BCs
If any mesh face uses `BoundaryFlag::user`, `ApplyPhysicalBCs` calls the
problem-generator `user_bcs_func` after `Z4cBCs`.

---

## AMR: `Z4c_AMR`

Configured via `<z4c_amr>`:
- `method`: `trivial`, `tracker`, `chi`, or `dchi`
- `chi_min` (for `chi`), `dchi_max` (for `dchi`)
- `radius_N_rad`, `radius_N_reflevel`: enforce minimum level inside spheres

`Refine()` dispatches to the selected criterion and always applies the
radius-based minimum refinement. The AMR driver is typically called from
problem generators (e.g., puncture setups) via a user refinement callback.

---

## Weyl Scalars and Wave Extraction

### Weyl scalars (`Z4cWeyl<NGHOST>`)
- Computes Psi4 (real/imag) from ADM variables into `u_weyl`.
- Executed at end-of-step when wave extraction is scheduled.

### Wave extraction (`WaveExtr`)
- Uses `SphericalGrid` objects created in the Z4c constructor.
- Interpolates `u_weyl` to spheres and decomposes into spin-weighted spherical
  harmonics (s=-2, l=2..8 => 77 modes).
- Writes to `waveforms/rpsi4_real_####.txt` and `waveforms/rpsi4_imag_####.txt`.
- Scheduling:
  - `nrad_wave_extraction=0` disables.
  - `waveform_dt` and `last_output_time` gate output; `InitRecvWeyl` updates
    `last_output_time` when output is due.

---

## Compact Object Tracker

Configured via `<z4c>`:
- `co_N_type`: `BH`/`BlackHole` or `NS`/`NeutronStar`
- `co_N_x`, `co_N_y`, `co_N_z` (initial position)
- `co_N_reflevel`, `co_N_radius` (AMR hints)
- `co_N_out_every` (output cadence)
- `filename` (prefix) and `<job>/basename` for output naming

The tracker interpolates shift (and for NS, fluid velocity and metric) using
`LagrangeInterpolator`, advances the position, and writes a per-object file.

---

## Input Parameters Summary

### `<z4c>`
- `chi_psi_power`, `chi_div_floor`, `diss`, `eps_floor`
- `damp_kappa1`, `damp_kappa2`
- `lapse_harmonicf`, `lapse_harmonic`, `lapse_oplog`, `lapse_advect`
- `shift_Gamma`, `shift_advect`, `shift_alpha2Gamma`, `shift_H`, `shift_eta`
- `use_z4c` (false disables Theta, BSSN-like behavior)
- `user_Sbc` (Sommerfeld on user boundaries)
- `extrap_order` (2..4, clamped to <= nghost)
- Wave extraction:
  - `nrad_wave_extraction`, `extraction_nlev`,
    `extraction_radius_1`, `extraction_radius_2`, ...
  - `waveform_dt`
- Compact object tracker:
  - `filename`, `co_N_type`, `co_N_x/y/z`, `co_N_reflevel`, `co_N_radius`,
    `co_N_out_every`

### `<z4c_amr>`
- `method` (`trivial`, `tracker`, `chi`, `dchi`)
- `chi_min`, `dchi_max`
- `radius_N_rad`, `radius_N_reflevel`

---

## Cautions

- `NewTimeStep` uses only the minimum grid spacing; there is no wave-speed
  estimate for Z4c.
- `eps_floor` is read but not used in the active `AlgConstr` implementation
  (the older guarded path is commented out).
- `psi_out` is allocated with a hard-coded mode count (77 for l=2..8). Changing
  `lmax` in `WaveExtr` requires updating this allocation.
