# Passive CGL thermodynamics

The supported scope is uniform periodic Newtonian MHD with dynamic
evolution and the HLLE solver, including first-order flux correction, local
Landau-fluid heat flux, collisions/limiters, and turbulence forcing.

With `<mhd>/passive=true`, density, momentum and magnetic fields follow the
isothermal MHD system with the required positive `iso_sound_speed`. Thermal
pressures do not contribute to the force or advective wave speed. Their
primitive slots retain `p_parallel` and `p_perp`, while the conserved slots are

```
IEN: J = rho [ln(p_parallel) + 2 ln(B_eff) - 3 ln(rho)]
IAN: A = rho [ln(p_perp) - ln(p_parallel) + 2 ln(rho) - 3 ln(B_eff)]
B_eff = max(|B|, bfloor)
```

Both quantities are transported by the isothermal mass flux. The flow uses
the exact native isothermal flux kernel; a second reconstruction supplies
thermal scalar fluxes without rewriting flow or electric-field components. In smooth ideal
flow this gives the double-adiabatic material laws
`p_parallel B² / rho³ = constant` and `p_perp / (rho B) = constant`.
The thermal energy density `U = p_perp + p_parallel/2` consequently obeys the
CGL pressure work, including the nonzero `Delta_p bb:grad(v)` contribution
in incompressible flow. It is not obtained by subtracting isothermal kinetic
and magnetic energy from a stored total energy.

LF substeps temporarily store **U in IEN** and `p_perp/B_eff` in IAN. Both
slots return to J/A before hyperbolic evolution, histories, and restarts.
Collisions and limiter relaxation change both invariants at fixed physical U;
no-op maps preserve the conserved bits. Below the magnetic floor, the thermal
pair is isotropized at fixed U. Pressure floors may add thermal energy and are
reported through the existing repair counters. Thermodynamic repairs never
change momentum or trigger a flow flux correction. Density floors and FOFC
flow decisions use the same arithmetic as the isothermal EOS.

Identical flow evolution requires identical initial flow, reconstruction,
forcing realization, boundaries, and **outer timestep sequence**. LF explicit
stability limits or `sts_max_dt_ratio` can shorten that sequence. For a pure
flow-identity comparison with STS, set `sts_max_dt_ratio=-1` so its optional
cycle cap is inactive; explicit comparisons must keep the LF stability bound
above the shared advective step, for example with a weak heat-flux fixture.
LF transport does not modify density, momentum, or magnetic fields; finite
density-floor repairs remain recorded exceptions to the frozen-density state.

Use `<time>/cgl_lf_max_chunk_ratio = 0` for new passive runs. This is the
runtime default: LF remains enabled, with one RKL2 sequence per half-cycle
sweep. A positive value such as `32` subdivides each sweep and refreshes the
LF bound between chunks; it does not cap the heat-flux amplitude or stage count.
The corrected-geometry [continuation comparison](https://github.com/dfielding14/athenak/blob/c/cgl-lf-surgical-rebuild/docs/cgl_turbulence_analysis.md)
supports cap 0 for quick passive comparisons, while instability-threshold
populations remain sensitive to the choice. Existing inputs and restarts
with an explicit positive cap retain it; override that existing key to adopt 0.

Passive histories label the conserved integrals `cgl-J` and `cgl-A`, and add
`thermal-U` for physical thermal energy. Kinetic and magnetic columns keep
their usual meaning. Passive conserved binary fields use `cgl_J` and `cgl_A`;
primitive `eint` retains the legacy meaning `p_parallel`, with `p_perp` exposed
separately. Recorded forcing work is the kinetic-energy change. Passive
pressure-work recording uses the centered physical thermal stress as a
diagnostic; that stress is not applied as a force to the isothermal fluid.
For periodic boundaries, the thermal budget is `Delta integral(U) + W = 0`,
where W is the recorded mechanical pressure work. This is a convergence
relation, not an exact finite-step conservation identity for advected J/A.

Full-precision restarts carry `passive_restart_encoding=1`. Active/passive
mode changes and incompatible encoding versions on restart are rejected.
The supported built-in initializers are `turb`, `cgl_lf_paper`, and the regression
initializer `cgl_passive_validation`. Other
built-in problem generators are rejected. Custom problem generators must
initialize physical primitives through the
EOS `PrimToCons` interface, which honors the current J/A or U/mu
representation for every destination buffer. Custom source terms must likewise
respect that contract: a direct IEN energy increment would corrupt J.
The unredefined `mhd_sgs` energy-moment output is rejected for passive runs.

AMR/SMR, nonperiodic physical boundaries, viscosity, hyperviscosity,
resistivity, cooling, shearing boxes, kinematic evolution, and coupled fluids
or radiation remain explicitly unsupported. Active CGL's mesh and boundary
qualification does not establish passive support, including uniformly
increasing refinement over time.

The passive model transports material invariants through shocks. It does
not add irreversible shock heating; interpreting it as a control for a
shock-heated active calculation requires this distinction. Smooth pressure
work, heat flux, and collisional relaxation remain part of the model.

Acceptance tests are `test_cgl_passive_{cpu,gpu}.py` and
`test_cgl_c2p_pressure_floor_cpu.py::test_cgl_finite_backup_survives_independent_c2p`.
The latter checks active/passive finite residuals through independent recovery.
The shared `passive_acceptance.py` exposes
`--suite forced` for focused multiple-rank flow/restart checks. They compare
full-precision flow and restart arrays, timestep sequences, reconstructors,
forcing, floors, independent linear modes, thermal advection and pressure work.
The nonlinear temporal case uses stationary dynamic isothermal flow and
refines local LF against a three-timestep, Richardson-extrapolated explicit
reference, checking the explicit method's first-order convergence and reference
uncertainty. These smooth checks do not certify accuracy in depleted turbulent
cells. The retained depleted-pressure stress failure, subsequent LF interval
refinement, and passive turbulent milestone used the earlier limiter operator;
see the [work order](https://github.com/dfielding14/athenak/blob/c/cgl-lf-surgical-rebuild/docs/cgl_lf_surgical_rebuild_work_order.md) for their evidence
and the separate qualification of the restored finite two-tier model.
