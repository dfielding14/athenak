# Scalar Mixing Velocity-Method Notes

`scalar_mixing` now supports three divergence-free velocity initialization
methods:

- `turb_velocity_method = projection`: projection of random vector Fourier modes
- `turb_velocity_method = stream_2d`: 2D stream-function construction
- `turb_velocity_method = clebsch`: 3D Clebsch construction

The older boolean alias `turb_use_stream_function` is still accepted for
backward compatibility, but it is deprecated. On `nx3 = 1` meshes it maps to
`stream_2d`; on 3D meshes it maps to `clebsch`.

`scalar_mixing` now always installs the divergence-free staggered face velocity
used by the scalar-only transport operator. Projection mode is always
solenoidal; `turb_sol_frac` and `divfree_scalar_flux` are no longer accepted by
this problem generator.

## 2D Convention

The target `turb_expo` is interpreted as the shell-integrated velocity spectrum
`E_u(k) ~ k^{-expo}`.

- Projection mode in true 2D uses per-mode velocity amplitudes
  `|u_k| ~ k^{-(expo+1)/2}`.
- 2D stream-function mode uses
  `u = (dpsi/dy, -dpsi/dx, 0)`,
  so `|psi_k| ~ |u_k| / k ~ k^{-(expo+3)/2}`.

This keeps the velocity spectrum target consistent between projection and
stream-function initialization on `nx3 = 1` meshes.

## 3D Clebsch Model

In 3D, the Clebsch path uses

`u = grad(phi1) x grad(phi2)`.

For Fourier modes `p` and `q`, the exact velocity coefficient is

`u_hat(p + q) = (p x q) phi1_hat(p) phi2_hat(q)`.

The implementation now mirrors the retained-mode notebook construction rather
than solving a separate shell-fit problem. The user specifies `turb_alpha`, and
AthenaK converts that to

- `s_phi = 3 + 2*alpha` for `alpha >= 0`
- `s_phi = 3 + alpha` for `alpha < 0`
- target velocity shell slope `s_u = 2*alpha + 1`

It then:

- builds two independent retained scalar catalogs over `[turb_nlow, turb_nhigh]`
- draws Gaussian scalar coefficients shell by shell
- rescales each retained shell so `E_phi(k) ~ k^{-s_phi}` is exact up to roundoff
- forms an edge-centered vector potential
  `A = 0.5 (phi1 grad(phi2) - phi2 grad(phi1))`
  and takes its discrete curl onto scalar-transport faces
- normalizes the realized velocity field to zero mean and the requested `turb_v_rms`

Because the velocity is a quadratic function of the two scalar fields and the
retained catalog is thinned by `turb_k_crit`, the realized velocity slope is
only approximate in practice.

## Leakage Definition

For the 3D Clebsch generator, out-of-band leakage is reported as

`sum_{K < nlow or K > nhigh} E_u(K) / sum_{K >= 1} E_u(K)`.

The code does not run a separate quality gate anymore. Instead the diagnostics
sidecar exports the retained scalar shell spectra, the latent scalar slices, and
the normalization metadata needed to audit the realization directly.

## Input Contract

- On `nx3 = 1` meshes, `stream_2d` always uses the 2D construction and forces
  `v3 = 0`, regardless of x3 boundary flags.
- The 3D Clebsch sidecar keeps `phi1`, `phi2`, their shell spectra, the face
  divergence diagnostics, and the normalization metadata. It no longer records
  fit starts, sector responses, quality-gate state, or a stale retained-velocity
  shell spectrum.
