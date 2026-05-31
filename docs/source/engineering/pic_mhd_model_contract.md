# MHD-PIC Runtime Model Contract

## Scope

AthenaK exposes several particle and MHD-PIC execution paths. They are not
interchangeable. A run that claims to reproduce Sun & Bai (2023),
arXiv:2304.10568v1, must select an explicit paper mode and must not assemble a
model from independent legacy booleans.

The runtime identity is selected with `<particles>/pic_physical_mode`.

| Mode | Purpose | CR state | Induction policy | Gas feedback |
| --- | --- | --- | --- | --- |
| `engineering` | Backward-compatible development and proxy tests | Historical velocity slots | Legacy opt-in current-to-CT source remains available | Legacy opt-in policies |
| `paper_test_particle` | Particle-only analytical tests of paper mechanics | Mass-normalized momentum `p/m` | No CR-current CT source | Disabled |
| `paper_mhd_pic` | Sun & Bai paper reproduction | Mass-normalized momentum `p/m` | Frozen-in `cE = -u x B`; no CR Hall term | Conservative momentum and kinetic-energy deltas for ideal MHD; exact-isothermal delta-f and the Q-006 runtime-local full-f carrier use momentum-only feedback |
| `extended_mhd_pic` | Separately named extensions requiring separate qualification | Mass-normalized momentum `p/m` | Extension-specific, never implied by paper mode | Extension-specific and recorded |

`engineering` is not a paper-reproduction mode. It exists to preserve the
historical interface while migration tests are written. New publication decks
must choose one of the explicit paper or extension modes.

## Particle State

The shared particle payload slots `IPVX`, `IPVY`, and `IPVZ` retain their
historical names for restart-layout compatibility. Their meaning depends on the
runtime mode:

| Mode family | Slot meaning |
| --- | --- |
| `engineering` | Coordinate velocity components |
| `paper_test_particle`, `paper_mhd_pic`, `extended_mhd_pic` | Mass-normalized momentum components `p/m` |

In a momentum-state mode, the derived quantities are

```text
gamma = sqrt(1 + |p/m|^2 / C^2)
v     = (p/m) / gamma
ekin  = (gamma - 1) C^2
```

where `C = <particles>/pic_cr_light_speed` is the artificial CR light speed.
Particle motion, current deposition, VTK output, tracked-particle output, and
cell-crossing timestep checks use derived velocity. Restart files preserve the
stored state and must record enough metadata to reject incompatible reads.

The `drift` and `rk4_gravity` pushers continue to use velocity slots.
Momentum-state paper modes require a CR Boris pusher.

## Paper Equations

For each CR super-particle, paper mode integrates

```text
d x / dt       = v
d (p/m) / dt   = (q/mc) (cE + v x B)
cE             = -u x B
```

The configured charge-to-mass slot stores the normalized `q/(mc)` factor used
by the AthenaK units. The relativistic Boris rotation evaluates its magnetic
rotation with the Lorentz factor after the first electric half-kick.

For `paper_mhd_pic`, the MHD induction update remains the ideal-MHD constrained
transport update. Deposited CR current must not be added directly to the CT
electric field. After a completed particle push, an ideal-MHD gas receives the
negative of the deposited CR momentum and relativistic kinetic-energy changes.
Exact-isothermal paper delta-f uses momentum-only feedback. The separately
named `q006_paper_multispecies_oscillation_runtime_local` generator admits the
same momentum-only contract for its bounded full-f Section 5.3 mechanics
carrier; other exact-isothermal full-f paper-mode compositions fail closed.

## Stage Ordering

The target second-order paper sequence is:

1. Deposit initial-position CR charge and current moments needed by the gas
   predictor.
2. Advance CR positions from the initial state to the midpoint with initial
   derived velocities.
3. Advance the MHD predictor and apply CR source terms.
4. Interpolate midpoint gas velocity and magnetic field to each particle.
5. Form midpoint `cE = -u x B` and advance `p/m` with the relativistic Boris
   update.
6. Record per-particle momentum and relativistic kinetic-energy deltas.
7. Advance CR positions from midpoint to the final state with final derived
   velocities.
8. Deposit the recorded deltas and subtract them from the gas update.
9. Apply physical boundaries, particle migration, moment synchronization, and
   constrained transport in the documented task order.

Implementation tests must trace the actual task sequence. This document records
the required sequence, not an assertion that every step is already complete.

## Mesh And Boundary Contract

Particle state, old positions, delta channels, deposited moments, and
particle-owned communication helpers must remain valid across:

- periodic, reflecting, and outflow boundaries;
- MeshBlock migration between ranks;
- SMR fine/coarse interfaces;
- AMR refine/derefine reconstruction;
- restart and decomposition changes.

The retained-state inventory and the refinement-interface policy decision are
recorded in
{doc}`pic_amr_lifetime_and_interface_policy`. Paper-mode AMR uses the
`paper_smooth` cell-centered restriction/exchange/prolongation path. The
optional `conservative` interface policy is not retained as a qualified
production mode.

Paper-mode deposition uses TSC interpolation/deposition unless a separately
named extension explicitly documents a different policy. A particle timestep
must satisfy both the configured maximum cell crossing bound and the configured
gyro-angle bound.

## Delta-F And Expanding Box

Quiet-start placement remains available as the separately named
`pic_deltaf_mode=quiet_start` variance-reduction technique. In explicit
paper/extension modes, `pic_deltaf_mode=physical` stores the initial analytic
background value, evolves `1 - f0(t,x,p)/f(0,x0,p0)`, deposits perturbation
moments, fingerprints the model in restart schema version 7, and emits VTK
weight diagnostics. Quantitative CRSI and CRPAI qualification remains a
release gate.

`pic_expanding_box_mode=on` evolves directional scale factors, applies
half-step particle momentum transforms around the Boris push, drifts in
comoving coordinates, rescales gas conserved variables, stores raw
face-centered arrays as divergence-preserving comoving magnetic fluxes,
derives physical face fields for MHD consumers, maps edge EMFs before CT, and
fingerprints its law and rates in restart schema version 7. Active-MHD
expansion is intentionally restricted to parser-qualified compositions until
each additional source package is mapped and tested. The bounded non-delta-f
conservative source split deposits physical-volume moments, records the
particle EM impulse in the final physical frame, and applies its opposite gas
impulse after the gas expansion map. The adaptive physical delta-f extension
uses a separate endpoint-normalized analytic `rho E + J x B` feedback path;
its host oracle checks source normalization and ordering, not opposite-impulse
conservation or transport calibration. Both admitted paths require
cell-centered `cc_convert` deposition, momentum and energy feedback,
`mhd_src_terms` ordering, and Hall mode off. Built-in MHD `hst` output
integrates physical cell volumes and physical magnetic fields; user-defined
history callbacks remain fail-closed because their volume contract is not
defined.
Appendix-level MHD and particle analytical qualification remains a release
gate.

## Extension Boundary

The paper formulation neglects the CR-induced Hall term. Any Hall-capable
implementation belongs to `extended_mhd_pic`, with its equation, normalization,
applicability envelope, tests, and claims registered separately. Entity Toolkit
may be used as a frozen comparative reference for relativistic particle
mechanics and compatible deposition micro-oracles; it is not an oracle for
AthenaK paper-mode induction or gas feedback.

The currently implemented `pic_cr_hall_mode=current_to_ct_experimental`
extension is a deliberately narrow source experiment:

```text
cE_CT = cE_ideal + alpha_H P_edge[J_CR]
alpha_H = <particles>/couple_j_to_efield_coeff
```

where `P_edge` is either the selected cell-centered-to-edge conversion or the
explicit staggered-current representation. This mode requires
`pic_physical_mode=extended_mhd_pic` and coupled particle moments. The host
manufactured-source smoke compares Hall-off, `+alpha_H`, and `-alpha_H` runs
and requires the magnetic-field increments to be nonzero and odd in
`alpha_H`.

This source-isolation oracle does not establish a derived CR-induced Hall
normalization, a Hall Bell dispersion relation, a nonlinear applicability
envelope, or shock-front validity. Those remain separate extension release
gates.

### Reduced Ion-Neutral Friction

The `pic_wave_damping_mode=ion_neutral_friction` extension implements the
bounded high-frequency, static-neutral reduction used for ion-neutral-damped
wave experiments:

```text
p_ion,perp(t + dt) = p_ion,perp(t) exp(-nu_in dt)
```

where `nu_in = <particles>/pic_ion_neutral_collision_rate`. For non-expanding
runs, the exact map is applied once after the explicit RK source update. For
expanding runs, it is applied once after the endpoint physical-frame feedback
map. Ideal-MHD total energy loses the removed transverse ion kinetic energy.
This mode requires
`pic_physical_mode=extended_mhd_pic`, a positive collision rate, and an active
coupled MHD background.

This is not the repository's general two-fluid `ion-neutral` module. It assumes
static neutrals and does not establish a CRSI dispersion comparison, a damping
envelope, GPU parity, or nonlinear saturation. The host manufactured-source
oracle verifies the exact exponential momentum factor and energy sink only.

### Adaptive Delta-F

The `pic_deltaf_adapt_mode=global_bikappa_moments_experimental` extension fits
the global bi-kappa reference state at a fixed
`pic_deltaf_adapt_interval`. In an x1-parallel frame it computes

```text
xi = (2/pi) sum(weight p_perp) / sum(weight |p_parallel|)
p0 = A(kappa) sum(weight sqrt(xi^4 p_parallel^2 + xi^2 p_perp^2))
                  / sum(weight)
A(kappa) = sqrt(pi kappa) (kappa - 1) Gamma(kappa - 1/2)
           / (2 Gamma(kappa + 1))
```

and evaluates the fitted background with the corresponding `xi^4`, box-volume,
and fitted-`p0` normalization. Restart schema version 7 preserves the fitted
`xi`, fitted `p0`, and cadence bucket so a restart does not introduce an
unrequested refit.

The bounded implementation requires `extended_mhd_pic`, physical
`kappa_aniso` delta-f, expanding-box mode, `kappa > 1`, zero configured drift,
unit configured anisotropy scales, and a positive fit interval. Its host oracle
checks the closed-form two-species fit, parser guards, restart fingerprint
rejection, and exact uninterrupted-versus-restarted state. It does not establish
CRPAI transport calibration, a saturated-state scattering rate, or GPU/MPI
qualification.
