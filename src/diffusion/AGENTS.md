# AGENTS.md

## Purpose
This directory implements diffusion-style physics for Hydro and MHD: thermal
conduction, isotropic shear viscosity, and Ohmic resistivity. It also provides
current-density helpers used in resistive electric-field calculations.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Conduction
- `conduction.hpp` / `conduction.cpp`: `Conduction` class, isotropic and
  temperature-dependent conductivity, optional saturated heat flux, and dt
  constraints.

### Viscosity
- `viscosity.hpp` / `viscosity.cpp`: `Viscosity` class implementing isotropic
  Navier-Stokes shear viscosity and viscous timestep estimates.

### Resistivity
- `resistivity.hpp` / `resistivity.cpp`: `Resistivity` class for Ohmic diffusion,
  resistive electric fields (`OhmicEField`) and resistive energy flux
  (`OhmicEnergyFlux`).

### Current density helper
- `current_density.hpp`: inline `CurrentDensity` helper that computes edge-centered
  current density from face-centered magnetic fields, used by resistivity kernels.

---

## How Diffusion Terms Are Used
- Hydro constructs `Viscosity` and `Conduction` from `<hydro>` options.
- MHD constructs `Viscosity`, `Conduction`, and `Resistivity` from `<mhd>` options.
- Flux additions are applied during flux construction in the Hydro/MHD task flow.
- Each module maintains its own `dtnew` for diffusion stability checks, used by
  `NewTimeStep` in Hydro/MHD.

---

## Conduction Details
- **Inputs** (from block passed to constructor, `hydro` or `mhd`):
  - `conductivity` (Real) -> constant kappa.
  - `tdep_conductivity` (bool) -> enable temperature-dependent kappa.
  - `cond_ceiling` (Real) -> cap on kappa when tdep is enabled.
  - `sat_hflux` (bool) -> apply harmonic saturation of heat flux.
- **Modes**:
  - `IsotropicHeatFlux`: constant kappa, isotropic heat conduction.
  - `TempDependentHeatFlux`: Parker/Spitzer law via `KappaTemp`, with optional
    saturation using local gradients and pressure.
- **Units**: temperature-dependent kappa uses `Units` to convert to cgs and back.
- **Timestep**: `NewTimeStep` computes conduction-limited dt (skipped if saturated).
- **Constraints**: conduction only supports ideal-gas EOS (checked in constructor).

---

## Viscosity Details
- **Inputs**: `viscosity` (nu) from `<hydro>` or `<mhd>`.
- **Implementation**: isotropic shear viscosity adds momentum and energy fluxes
  using velocity gradients and shared-memory scratch arrays for 1D pencils.
- **Timestep**: dt is computed once per pack at construction based on grid spacing.

---

## Resistivity Details
- **Inputs**: `ohmic_resistivity` from `<mhd>`.
- **OhmicEField**: adds `eta * J` to edge-centered electric fields, using
  `CurrentDensity` to compute J on edges. Handles 1D/2D/3D cases explicitly.
- **OhmicEnergyFlux**: adds Poynting flux from resistive E-fields to energy fluxes.
- **Timestep**: dt computed at construction based on grid spacing and eta.

---

## Extension Points and Cautions
- New diffusion operators should mirror existing flux-addition patterns and update
  per-module dt constraints.
- Resistivity currently assumes Ohmic diffusion only; Hall/ambipolar would require
  new current-density and E-field formulations.
- Temperature-dependent conduction relies on ideal EOS and `Units`; keep those
  assumptions consistent if adding new EOS models.
- Scratch-array usage in viscosity and resistivity kernels assumes 1D pencil
  patterns; keep access ranges aligned with existing mesh indices.
