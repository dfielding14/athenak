# CGL Landau-Fluid Unit Inputs

`lf_k_parallel` is the Snyder-Hammett-Dorland closure wavenumber `|k_parallel|`.
It is not a conductivity.

`lf_coefficient_mode = local` is the default physics mode. It evaluates the
closure coefficient from the live face state with `c_parallel = sqrt(p_parallel/rho)`.

`lf_coefficient_mode = background` is for controlled linear/reference tests. It
freezes only the coefficient value `c_parallel = lf_c_parallel0`; the pressure,
density, magnetic-field gradients, and anisotropy state still come from the live
solution.

CGL Landau-fluid heat flux requires `<mhd>/eos = cgl` and
`<mhd>/conductivity_integrator = sts`.

The LF face fluxes are capped with the Squire et al. (2023) equation 3.2 form,
`q = q_L*q_max/(q_max + abs(q_L))`. The current implementation does not yet
include collision-frequency-dependent suppression of the closure coefficients.

The `cgl_lf_paper_oblique_wave` and `cgl_pure_paper_oblique_wave` inputs use
the Figure 17 background state from Squire et al. (2023) and compare against a
linearized offline reference for the same initialized perturbation.
