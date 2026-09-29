# CGL Method Physics Primer

This page explains the physics model behind AthenaK's non-relativistic CGL MHD
and optional Landau-fluid heat flux. It is a method guide, not a production
claim: validation status and campaign workflows live in
[CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md), and the code
map lives in [CGL Landau-Fluid Code Guide](cgl_landau_fluid_code_guide.md).

## Model Lineage

The implementation follows the reduced-fluid model class used in the modern
magneto-immutability literature:

- Chew-Goldberger-Low double-adiabatic MHD supplies separate parallel and
  perpendicular pressures.
- Landau-fluid closures approximate collisionless, field-aligned heat fluxes
  using a prescribed parallel closure scale.
- Microinstability regulation represents unresolved mirror and firehose
  scattering through finite-rate soft limiters, with energy-preserving fluid
  and optional backup walls.

The Squire et al. and Majeski et al. turbulence studies use this model class to
study how pressure anisotropy changes high-beta MHD turbulence. AthenaK is an
independent implementation in that model lineage. Agreement with those papers
must be established through the validation and analysis workflows; the method
definition alone is not a reproduction claim.

## Pressure Tensor

CGL MHD evolves two thermal pressures:

$$
p_\parallel,\qquad p_\perp ,
$$

defined relative to the local magnetic-field direction

$$
\boldsymbol{b} = \boldsymbol{B}/|\boldsymbol{B}|.
$$

AthenaK uses

$$
\Delta p = p_\perp - p_\parallel
$$

and the pressure tensor

$$
\mathbf{P} = p_\perp \mathbf{I} - \Delta p\,\boldsymbol{b}\boldsymbol{b}.
$$

The internal energy density associated with the two pressures is

$$
u_{\rm th} = p_\perp + \frac{1}{2}p_\parallel .
$$

At high beta, a fractional pressure anisotropy of order $\beta^{-1}$ can still
matter dynamically because magnetic tension is modified by the anisotropic
stress. This is why mirror and firehose bounds are central to the model even
when $p_\parallel$ and $p_\perp$ differ only weakly in relative terms.

## Double-Adiabatic Limit

Without heat fluxes, collisions, limiters, explicit sources, or numerical
floors, the CGL thermodynamics conserve the double-adiabatic invariants

$$
\frac{p_\perp}{\rho B},\qquad
\frac{p_\parallel B^2}{\rho^3}.
$$

AthenaK stores a conserved anisotropy variable equivalent to

$$
A = \rho \ln\left(
  \frac{p_\perp}{p_\parallel}\frac{\rho^2}{B^3}
\right),
$$

where $B=|\boldsymbol{B}|$. This is the sixth MHD conserved component outside
Landau-fluid sweeps. The code guide describes the temporary magnetic-moment
reinterpretation used during LF transport.

## Magneto-Immutability

The model is useful for high-beta turbulence because pressure anisotropy
responds directly to flow motions that change magnetic-field strength. In ideal
MHD,

$$
\frac{D\ln B}{Dt} =
\boldsymbol{b}\boldsymbol{b}:\nabla\boldsymbol{u}
- \nabla\cdot\boldsymbol{u}.
$$

For incompressible motions the first term is the field-strength-changing strain.
For compressible motions the divergence term must be retained; do not interpret
$\boldsymbol{b}\boldsymbol{b}:\nabla\boldsymbol{u}$ alone as $D\ln B/Dt$.

The pressure work decomposes as

$$
\mathbf{P}:\nabla\boldsymbol{u}
= p_\perp \nabla\cdot\boldsymbol{u}
- \Delta p\,\boldsymbol{b}\boldsymbol{b}:\nabla\boldsymbol{u}.
$$

The anisotropic term can resist the motions that would change $|B|$. The
Majeski/Squire comparisons therefore emphasize diagnostics such as
field-strength PDFs, pressure-anisotropy PDFs, pressure and magnetic spectra,
pressure-work transfer, and suppression of field-aligned strain.

## Landau-Fluid Heat Flux

The optional LF closure adds field-aligned heat transport for the parallel and
perpendicular temperatures. AthenaK evaluates local gradients of

$$
T_\parallel = p_\parallel/\rho,\qquad
T_\perp = p_\perp/\rho,
$$

projects them along $\boldsymbol{b}$, and uses a prescribed positive
`lf_k_parallel` as the closure wavenumber magnitude. This parameter is not a
measured turbulent wavenumber and not a nonlocal solve; smaller
`lf_k_parallel` gives larger heat-flux coefficients and stronger smoothing.

In the collisionless normal range, the implemented ratios are equivalent to

$$
\frac{q_\parallel^L}
{\sqrt{8/\pi}\,c_\parallel p_\parallel}
= -\frac{\rho\,\nabla_\parallel T_\parallel}
{p_\parallel\,|k_\parallel|},
$$

and

$$
\frac{q_\perp^L}
{\sqrt{2/\pi}\,c_\parallel p_\perp}
= -\frac{1}{|k_\parallel|}
\left[
  \frac{\rho\,\nabla_\parallel T_\perp}{p_\perp}
  - \left(1-\frac{p_\perp}{p_\parallel}\right)
    \nabla_\parallel\ln B
\right].
$$

Background collisions and active soft/backup limiter rates add at each face,
reducing the heat flux. The finite-collision coefficients are

$$
\chi_\parallel=\frac{8c_\parallel^2}{\sqrt{8\pi}c_\parallel|k_\parallel|+(3\pi-8)\nu_{\rm eff}},
\qquad
\chi_\perp=\frac{2c_\parallel^2}{\sqrt{2\pi}c_\parallel|k_\parallel|+2\nu_{\rm eff}}.
$$

The perpendicular form follows SHD97/Sharma and the BGK moment result; the
$+\nu$ printed in Squire et al. (2023) is deliberately not used. The code
supports either local $c_\parallel$ or a configured background `lf_c_parallel0`.

## Free-Streaming Cap

The unlimited heat flux is smoothly capped by a free-streaming-scale limiter.
Conceptually,

$$
q = \frac{q^L q_{\max}}{q_{\max} + |q^L|},
$$

with

$$
q_{\parallel,\max} = \sqrt{8/\pi}\,c_\parallel p_\parallel,\qquad
q_{\perp,\max} = \sqrt{2/\pi}\,c_\parallel p_\perp.
$$

This cap protects the fluid closure from unphysical large local fluxes. It
should not be described as a kinetic calculation of the heat flux.

## Collisions And Instability Regulation

AthenaK supports a background anisotropy-relaxation frequency `nu_coll` plus
optional mirror and firehose limiters. The code uses the dimensionless
anisotropy coordinate

$$
\frac{2\Delta p}{B^2}.
$$

The mirror threshold is

$$
\frac{2\Delta p}{B^2}=1,
$$

corresponding to $\Delta p = B^2/2$. The parallel-fluid firehose threshold is

$$
\frac{2\Delta p}{B^2}=-2,
$$

corresponding to $\Delta p=-B^2$. AthenaK also supports an oblique-firehose
activation policy at $2\Delta p/B^2=-1.4$ for legacy or diagnostic studies.

Finite-rate limiters relax anisotropy monotonically toward the selected soft
threshold after exact background decay. Rates advance once per full timestep,
after the LF post sweep or after RK without LF. The defaults are numeric
`mirror_threshold=1` and `firehose_threshold=2`; the legacy oblique alias maps
to 1.4. Enabling either limiter requires an explicit `limiter_nu_coll`.

At a scheduled wall update, the fluid bound $\Delta p\geq-B^2$ always applies.
Enabled backup walls default to $[-B^2,B^2]$ and preserve
$p_\perp+p_\parallel/2$; they are independent of which soft limiter is enabled.
LF additionally applies walls after its pre sweep and after RK. The old
`limiter_hardwall=true` soft projection during primitive recovery is rejected;
a stiff finite rate approaches the soft threshold asymptotically.

These are fluid-scale regulation models, not resolved mirror, firehose,
ion-Larmor-radius, or kinetic phase-space physics.

## Active And Passive Runs

Active CGL applies the anisotropic pressure tensor in the momentum and energy
fluxes. The former passive-Delta model is disabled because its thermal-energy
equation is inconsistent; `passive=true` fails at construction until the WO2
redesign. Its isothermal signal-speed paths remain unit tested. Archived
active/passive comparisons are retained as historical evidence, not as a
currently executable controlled experiment.

## Turbulence Comparisons

Majeski/Squire-style turbulence studies typically use a periodic elongated box
with a guide field, modal Ornstein-Uhlenbeck forcing, active/passive contrasts,
heat-flux-scale scans, and limiter/collisionality scans. The most relevant
observables are not just energies; they include pressure-anisotropy occupancy,
pressure and magnetic spectra, pressure-work transfer, field/strain alignment,
and the field-strength distribution.

For AthenaK, the forcing implementation and production input details are
documented separately in [Turbulence Driving](turbulence_driving.md) and
[CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md). A private
Majeski note reviewed during development is useful qualitative sanity evidence,
but it is not a public archival citation or a substitute for reproducible
validation products.

## What This Method Does Not Claim

Do not describe the CGL-LF implementation as:

- a kinetic or Vlasov calculation;
- a resolved mirror/firehose simulation;
- a finite-Larmor-radius model;
- a proof of continuum convergence;
- a completed reproduction of a specific paper without the corresponding
  retained run and analysis evidence;
- an irreversible heating measurement based only on signed pressure-work
  diagnostics.

Use this primer for the physics model, the
[CGL Landau-Fluid Code Guide](cgl_landau_fluid_code_guide.md) for the
implementation map, and
[CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md) for evidence
requirements.
