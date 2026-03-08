# Scalar Mixing Stream-Function Notes

`scalar_mixing` now supports two divergence-free velocity initialization methods:

- `turb_use_stream_function = false`: projection of random vector Fourier modes
- `turb_use_stream_function = true`: stream-function construction

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

In 3D, the stream-function path uses a Clebsch representation

`u = grad(phi1) x grad(phi2)`.

For Fourier modes `p` and `q`, the exact velocity coefficient is

`u_hat(p + q) = (p x q) phi1_hat(p) phi2_hat(q)`.

The implementation fits nonnegative shell weights `A_s` for the scalar
potentials, with shells `s = 1..turb_nhigh`, by building a deterministic shell
response tensor `T[K,s,t]` from the actually accepted mode catalogs. This means
the fit includes the retained `turb_k_crit` importance sampling and spectral
boosting used during mode selection.

The initial guess uses the local-triad asymptotic scaling:

- the number of contributing pairs at fixed output shell grows like `k^2`
- `|p x q|^2` contributes another `k^4`
- with `|phi_k| ~ k^{-beta}`, this gives `E_u(k) ~ k^{8 - 4 beta}`

Matching `E_u(k) ~ k^{-expo}` then suggests
`beta ~= (expo + 8) / 4`, which seeds the nonlinear shell fit.

## Leakage Definition

For the 3D Clebsch fit, out-of-band leakage is reported as

`sum_{K < nlow or K > nhigh} E_u(K) / sum_{K >= 1} E_u(K)`.

The solver penalizes leakage outside `[turb_nlow, turb_nhigh]` while matching
the requested in-band slope.

For regression on small boxes, low integer shell sums are noisy and can distort
the apparent slope even when the retained Fourier coefficients follow the
requested law. The added regression therefore fits ensemble-averaged mode
energy versus the exact discrete `|k|` values inside the requested band.

## Input Contract

- `turb_sol_frac` is ignored in stream-function mode because the construction is
  intrinsically solenoidal.
- On `nx3 = 1` meshes, stream-function mode always uses the 2D construction and
  forces `v3 = 0`, regardless of x3 boundary flags.
