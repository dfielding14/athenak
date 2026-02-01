# AGENTS.md

## Purpose
This directory implements equation-of-state (EOS) support and conserved <->
primitive conversions for Hydro and MHD in Newtonian, SR, and GR regimes. It
also contains the policy-based primitive solver framework used by dynamical GR
modules and advanced EOS tables.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core EOS interface
- `eos.hpp` / `eos.cpp`: `EquationOfState` base class, `EOS_Data` container, and
  derived-class declarations for isothermal/ideal, hydro/MHD, SR/GR variants.

### Ideal and isothermal EOS implementations
- `ideal_hyd.cpp`, `ideal_mhd.cpp`: nonrelativistic ideal-gas C2P/P2C.
- `isothermal_hyd.cpp`, `isothermal_mhd.cpp`: nonrelativistic isothermal C2P/P2C.
- `ideal_srhyd.cpp`, `ideal_srmhd.cpp`: SR ideal-gas C2P/P2C with root finding.
- `ideal_grhyd.cpp`, `ideal_grmhd.cpp`: GR ideal-gas C2P/P2C using metric transforms.

### Inline conversion helpers
- `ideal_c2p_hyd.hpp`: inline single-state C2P/P2C for hydro (Newtonian, SR, GR).
- `ideal_c2p_mhd.hpp`: inline single-state C2P/P2C for MHD (Newtonian, SR, GR).

### Primitive solver framework (policy-based)
- `primitive_solver_hyd.hpp`: adapter that plugs the policy-based primitive solver
  into AthenaK data structures (GRMHD-focused). Includes FOFC and excision handling.
- `primitive-solver/`: standalone EOS + primitive solver framework, including:
  - `primitive_solver.hpp`: root finding + primitive solve logic.
  - `eos.hpp`, `eos_policy_interface.hpp`, `error_policy_interface.hpp`.
  - `idealgas.hpp`, `piecewise_polytrope.*`, `eos_compose.*`.
  - `unit_system.*`, `reset_floor.hpp`, `ps_types.hpp`, `ps_error.hpp`.

---

## EOS_Data and Base Class
- `EOS_Data` stores gamma, isothermal sound speed, floors, and a Lorentz-factor
  ceiling. It also provides inline wave-speed helpers for Newtonian/SR/GR flows.
- `EquationOfState` supplies virtual `ConsToPrim` and `PrimToCons` interfaces for
  hydro and MHD signatures; derived classes implement only the relevant overloads.

---

## Conversion Flow and Floor Handling
- Each C2P path uses single-state helpers (`SingleC2P_*`) inside Kokkos loops.
- Floors are applied to density, pressure/internal energy, and sometimes entropy or
  temperature. Usage is tracked in `Mesh::ecounter` (dfloor, efloor, tfloor, vceil,
  c2p failures).
- `only_testfloors` is used by FOFC: if floors or failures are detected, the FOFC
  flag is set and the update path switches to first-order corrections.

---

## Relativistic Variants
- SR C2P for ideal gas uses the Galeazzi et al. root solve (false position), with
  velocity ceilings enforced via `gamma_max`.
- GR C2P converts to local SR variables using metric transforms, performs SR solve,
  and transforms back (see `TransformToSRHyd/MHD` in `ideal_c2p_*`).
- GR variants integrate excision masks (`Coordinates::excision_floor/flux`) to
  bypass or reset states inside the excised region.

---

## Policy-Based Primitive Solver (GRMHD)
The policy-based solver is separate from `EquationOfState` and is used where
advanced EOS tables or GRMHD diagnostics are needed.

- `PrimitiveSolverHydro` configures an EOS policy (IdealGas, PiecewisePolytrope,
  or CompOSE) and error policy, then exposes `ConsToPrim`/`PrimToCons` on AthenaK
  arrays.
- Uses ADM metric fields for densitization and converts between AthenaK arrays and
  the solver's primitive/conserved layouts.
- Supports floors, error caps, and detailed failure logging. Integrates excision
  handling in GR cases.

---

## Configuration Inputs
Common parameters parsed from `<hydro>` or `<mhd>` blocks:
- `eos`: `ideal` or `isothermal`
- `gamma`, `iso_sound_speed`
- Floors: `dfloor`, `pfloor`, `tfloor`, `sfloor`
- Relativistic ceiling: `gamma_max`

Primitive-solver-specific (for policy-based EOS):
- `c2p_tol`, `c2p_iter`, `c2perrs`
- EOS-specific parameters (e.g., CompOSE tables, unit system)

---

## Extension Points and Cautions
- When adding a new EOS, implement both C2P and P2C paths and ensure floors are
  accounted for with proper counters.
- Keep SR/GR transformations consistent with `Coordinates` and metric conventions;
  many routines assume `sqrt(-g)=1` in Kerr-Schild.
- If you change primitive variable definitions (e.g., use temperature vs energy),
  update `eos_data.use_e/use_t` and any reconstruction or solver logic that assumes
  internal energy.
- The policy-based primitive solver is independent of the legacy EOS classes; keep
  behaviors aligned if both are supported for the same physics module.
