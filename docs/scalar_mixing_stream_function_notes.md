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

The code no longer treats any single continuum asymptotic law as authoritative
for the 3D Clebsch fit. Instead it:

- builds the exact retained-mode shell-response tensor;
- records sector-resolved contributions (`local_local`, `low_high`,
  `high_low`, `high_high_cancel`);
- tries three provisional seed families internally:
  - the older local-triad guide,
  - a symmetric Camillo-style guide,
  - and a flat shell-power start;
- chooses the best fit by final loss, leakage, and in-band slope error.

## Leakage Definition

For the 3D Clebsch fit, out-of-band leakage is reported as

`sum_{K < nlow or K > nhigh} E_u(K) / sum_{K >= 1} E_u(K)`.

The solver penalizes leakage outside `[turb_nlow, turb_nhigh]` while matching
the requested in-band slope.

The 3D stream path now also applies a quality gate:

- deterministic fit must converge, have `abs(fitted_slope + expo) <= 0.5`, and
  leakage `<= 0.05`;
- the chosen realized spectrum must satisfy the same slope and leakage limits;
- if either check fails, initialization aborts unless
  `turb_stream_allow_poor_fit = true`.

The regression now relies on this diagnostics gate rather than on an exact-`|k|`
mode-slope fit for 3D Clebsch, because small discrete boxes can have strong
radius-by-radius nulls even when the shell-integrated retained response is
acceptable.

## Input Contract

- `turb_sol_frac` is ignored in stream-function mode because the construction is
  intrinsically solenoidal.
- `turb_stream_allow_poor_fit` is an expert escape hatch for audits. The default
  remains fail-closed.
- On `nx3 = 1` meshes, stream-function mode always uses the 2D construction and
  forces `v3 = 0`, regardless of x3 boundary flags.
