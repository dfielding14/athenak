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
above the shared advective step, for example with a weak heat-flux fixture. LF substeps themselves do not modify density, momentum,
or magnetic fields.

Passive histories label the conserved integrals `cgl-J` and `cgl-A`, and add
`thermal-U` for physical thermal energy. Kinetic and magnetic columns keep
their usual meaning. Passive conserved binary fields use `cgl_J` and `cgl_A`;
primitive `eint` retains the legacy meaning `p_parallel`, with `p_perp` exposed
separately. Recorded forcing work is the kinetic-energy change. Passive
pressure-work recording uses the centered physical thermal stress as a
diagnostic; that stress is not applied as a force to the isothermal fluid.

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
or radiation remain explicitly unsupported. Task5's
physical-boundary support applies to active CGL and does not remove the
independent passive boundary restriction.

The passive model transports material invariants through shocks. It does
not add irreversible shock heating; interpreting it as a control for a
shock-heated active calculation requires this distinction. Smooth pressure
work, heat flux, and collisional relaxation remain part of the model.

Permanent acceptance tests are `test_cgl_passive_{cpu,gpu,mpicpu,mpi_gpu}.py`.
They compare full-precision flow and restart arrays, actual timestep sequences,
all five reconstructors, forcing, floors, independent linear modes, thermal
advection, and the nonzero periodic pressure-work derivative. The MPI GPU
variant uses the existing `ATHENAK_RUN_MPI_GPU=1` opt-in and launcher contract.
Each test retains inputs, logs, complete restart payloads, and JSON results
under its pytest temporary directory.
