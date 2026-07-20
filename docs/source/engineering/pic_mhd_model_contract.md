# MHD-PIC Runtime Model Contract

## Scope

AthenaK exposes several particle and MHD-PIC execution paths. They are not
interchangeable. Runtime identity is described by the controls that actually
change the algorithm; there is no umbrella `pic_physical_mode` selector.

| Concern | Explicit control |
| --- | --- |
| VL2 predictor/corrector chronology | `<time>/integrator=vl2` |
| Particle representation | Cosmic rays with a Boris pusher store mass-normalized momentum `p/m`; drift and star particles store velocity |
| CR input interpretation | `<particles>/pic_cr_initial_state=velocity|momentum`; velocity input is converted once during initialization |
| Background evolution | `<particles>/pic_background_mode=coupled|passive_mhd|no_mhd` |
| Gas feedback | `<particles>/pic_feedback_mode=coupled|test_particle` plus the explicit deposition and feedback switches |
| CR-Hall closure | `<particles>/pic_cr_hall_mode=off|full` |

Fresh inputs containing the obsolete umbrella selector fail with a migration
error so that old implicit defaults cannot silently change the physics.
Schema-8 restart headers may still contain the old field; restart reads ignore
that field and validate the stored state against the explicit model metadata.
Legacy VL2 checkpoints must be restarted with `<time>/integrator=vl2`.

New coupled science decks use cosmic-ray Boris-TSC particles, second-order TSC
deposition, explicit momentum and energy feedback, and
`<time>/integrator=vl2`. Full Hall is one atomic physical choice; there is no
free-strength or experimental Hall mode.

## Particle State

The shared particle payload slots `IPVX`, `IPVY`, and `IPVZ` retain their
historical names for restart-layout compatibility. Cosmic rays advanced by a
Boris pusher store mass-normalized momentum `p/m`; drift and star particles
retain coordinate velocity. In a momentum-state path, the derived quantities
are

```text
gamma = sqrt(1 + |p/m|^2 / C^2)
v     = (p/m) / gamma
ekin  = (gamma - 1) C^2
```

where `C = <particles>/pic_cr_light_speed` is the artificial CR light speed.
Particle motion, current deposition, VTK output, tracked-particle output, and
cell-crossing timestep checks use derived velocity. Restart files preserve the
stored state and must record enough metadata to reject incompatible reads.

The `drift` and `rk4_gravity` pushers continue to use velocity slots. The
`pic_cr_initial_state` selector changes only how a CR's supplied initial values
are interpreted before conversion; it does not change the stored runtime state.

## Paper Equations

For each CR super-particle, the VL2/TSC path integrates

```text
d x / dt       = v
d (p/m) / dt   = (q/mc) (cE + v x B)
cE             = -(u + v_H) x B
```

Hall-off sets `v_H=0`. Full Hall uses the signed deposited moments

```text
Q_CR = q_CR/c
K_CR = J_CR/c
Q_e  = alpha_i rho + Q_CR
v_H  = (K_CR - Q_CR u)/Q_e
alpha_i = <particles>/pic_background_ion_q_over_mc
```

The particle pusher TSC-interpolates `u_g`, the predicted `v_H`, and `B` to
each particle. Full Hall uses an initial-current predictor at stage 1 and a
midpoint-current predictor at stage 2. The stage-1 grid update uses the same
predicted Hall electric field

```text
cE_H = -v_H x B.
```

After the stage-2 Boris kick, the deposited particle momentum-rate density is
the exact discrete Lorentz force `DPDT = dP_CR/dt`. The grid corrector therefore
uses the algebraically equivalent but exactly coupled form

```text
cE_H = -DPDT/(alpha_i rho)
```

for both induction and Hall energy transport. It does not reconstruct a second
stage-2 drift from charge density. The stage-2 gas source likewise uses the
deposited `DPDT` and relativistic kinetic-energy rate `DEDT` directly. Gas
momentum and energy receive the opposite particle exchange. Artificial
`pic_cr_light_speed` affects particle kinematics but does not enter the Hall
closure. See the
[CR-Hall paper-to-code map](pic_cr_hall_code_map.md) for the signed
normalization and exact staging.

The configured charge-to-mass slot stores the normalized `q/(mc)` factor used
by the AthenaK units. The relativistic Boris rotation evaluates its magnetic
rotation with the Lorentz factor after the first electric half-kick.

For the active `<time>/integrator=vl2` model, deposited CR current is never
added directly to the final edge EMF. Hall-off retains ideal-MHD CT; full Hall
adds the derived Hall correction to the face induction fluxes and adds the
matched `(cE_H x B)` term to the face total-energy flux before FOFC and the RK
update. Limited PLM states are selected by the ordinary density-flux sign and
use the exact staggered face-normal magnetic field. `CornerE` then constructs
the edge field from the resulting face terms. Exact-isothermal paper delta-f
uses momentum-only feedback. The separately named
`q006_paper_multispecies_oscillation_runtime_local` generator admits the same
momentum-only contract for its bounded full-f Section 5.3 mechanics carrier;
other exact-isothermal full-f compositions fail closed.

## Stage Ordering

The full-f particle chain runs after `MHD::CopyCons` and before `MHD::Fluxes`
on both VL2 stages. This makes the predictor or exact corrector data available
to the conservative face update and to FOFC. Hall-off deposits the generic
`moments` array in both stages. Full Hall deposits its dedicated
`cr_hall_moments` predictor in both stages, while its generic `moments` wrappers
run only in stage 2 to supply the realized `DPDT` and `DEDT` corrector.

1. At stage 1, deposit initial-position `Q_CR` and `K_CR`, synchronize their
   ghost contributions, and construct the initial `v_H`.
2. Use that drift in a scratch Boris half-kick that predicts momentum without
   changing the true particle state. Drift the true particle from the initial
   position to the midpoint with its initial derived velocity.
3. Add the predicted `-v_H x B` induction term and its matched Hall energy flux
   to MHD face fluxes before FOFC and the predictor RK update. Apply analytic
   gas feedback from `Q_CR cE + K_CR x B` and `K_CR dot cE` with the opposite
   sign. In 2D/3D, `CornerE` uses the same predicted cell Hall field; in 1D it
   copies the corrected face EMFs directly. CT advances the staggered magnetic
   field.
4. After midpoint particle boundary handling, deposit `Q_CR` and `K_CR` from
   the scratch-predicted momentum and construct the midpoint `v_H`.
5. TSC-interpolate midpoint `u_g`, predicted `v_H`, and `B`; perform the true
   full-step Boris kick at the midpoint; record momentum and relativistic
   kinetic-energy rates; deposit and synchronize `DPDT` and `DEDT`; then drift
   the particle from the midpoint to the endpoint with its final velocity.
6. Reconstruct the exact corrector `cE_H=-DPDT/(alpha_i rho)` to faces. Add its
   induction components and matched Hall energy flux before FOFC and the final
   RK update. Apply `-DPDT` and `-DEDT` to the gas in `MHDSrcTerms` with the
   ordinary stage weight.
7. In 2D/3D, `CornerE` uses the exact cell corrector for its Hall contribution;
   in 1D it copies the already corrected face EMFs directly. CT updates the
   staggered field, and full-Hall gas admissibility is checked against that
   updated field before conserved-to-primitive conversion. Physical boundaries,
   communication, and endpoint particle migration then complete the stage.

FOFC tests the composite flux-plus-particle-source update: its trial conserved
state receives the same stage-1 analytic or stage-2 exact feedback source as
the live state. If FOFC replaces an ordinary face flux with its first-order
fallback, it restores the Hall induction and energy terms together from the
donor cell selected by the replaced mass-flux sign before clearing the flags.
The deposited particle impulse itself is not clipped or redistributed.

Hall-off follows the same full-f predictor/corrector chronology but omits the
Hall-current predictor, Hall drift, and Hall face terms.

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
{doc}`pic_amr_lifetime_and_interface_policy`. Paper-mode AMR targets the
`paper_smooth` receiver-resolution TSC policy. The existing cell-centered
restriction/exchange/prolongation path is supporting infrastructure, not a
qualified substitute at fine/coarse interfaces. The optional `conservative`
interface policy is not retained as a qualified production mode.

Paper-mode deposition uses TSC interpolation/deposition unless a separately
named extension explicitly documents a different policy. A particle timestep
must satisfy both the configured maximum cell crossing bound and the configured
gyro-angle bound.

The first full-Hall implementation is deliberately uniform-grid and full-f.
AMR/SMR, delta-f, expanding boxes, `<mhd>/eos` values other than `ideal`, and
relativistic MHD fail closed in full mode. Viscosity, resistivity, and conduction
are outside the intended full-Hall qualification scope but are not currently
rejected by the parser. This is a known scope boundary, not a claim that the Hall
equations are optional on supported uniform-grid science runs.

Full Hall with `<mhd>/fofc=true` requires `<mesh>/nghost >= 3`; the expanded
stencil supplies the neighboring states needed by expanded Hall face
reconstruction and donor-cell replacement.

## Delta-F And Expanding Box

Quiet-start placement remains available as the separately named
`pic_deltaf_mode=quiet_start` variance-reduction technique. In explicit
paper/extension modes, `pic_deltaf_mode=physical` stores the initial analytic
background value, evolves `1 - f0(t,x,p)/f(0,x0,p0)`, deposits perturbation
moments, fingerprints the model in restart schema version 8, and emits VTK
weight diagnostics. Quantitative CRSI and CRPAI qualification remains a
release gate.

`pic_expanding_box_mode=on` evolves directional scale factors, applies
half-step particle momentum transforms around the Boris push, drifts in
comoving coordinates, rescales gas conserved variables, stores raw
face-centered arrays as divergence-preserving comoving magnetic fluxes,
derives physical face fields for MHD consumers, maps edge EMFs before CT, and
fingerprints its law and rates in restart schema version 8. Active-MHD
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

Sun & Bai (2023) used Hall-off induction, but the complete large-scale MHD-PIC
system in Bai et al. (2015) includes the CR-induced Hall term. AthenaK keeps the
VL2/TSC algorithm identity and selects the physical closure explicitly with
`pic_cr_hall_mode=off|full`. `full` is not a free-coefficient experiment. It
has a fixed signed equation, uniform-grid
applicability envelope, and two regime diagnostics, `max|R|` and `max Lambda`.
The frozen uniform candidates passed the focused conservation, Bell, restart,
and MPI qualification summarized in `MHD_PIC_NEXT_STEPS_GUIDE.md`.

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
This mode requires a positive collision rate and an active coupled MHD
background.

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
and fitted-`p0` normalization. Restart schema version 8 preserves the fitted
`xi`, fitted `p0`, and cadence bucket so a restart does not introduce an
unrequested refit.

The bounded implementation requires physical `kappa_aniso` delta-f,
expanding-box mode, `kappa > 1`, zero configured drift,
unit configured anisotropy scales, and a positive fit interval. Its host oracle
checks the closed-form two-species fit, parser guards, restart fingerprint
rejection, and exact uninterrupted-versus-restarted state. It does not establish
CRPAI transport calibration, a saturated-state scattering rate, or GPU/MPI
qualification.
