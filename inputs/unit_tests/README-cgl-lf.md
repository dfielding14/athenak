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
`<mhd>/cgl_heat_flux = landau_fluid`. Select `cgl_heat_flux_integrator = sts`
with `<time>/sts_integrator = rkl2`, or `explicit` with `sts_integrator = none`.

The perpendicular LF coefficient follows the SHD97/Sharma03 BGK convention:
`chi_perp = 2*c_parallel^2/(sqrt(2*pi)*c_parallel*lf_k_parallel + 2*nu_eff)`.
Here background collisions decay `p_perp-p_parallel` at rate `nu_coll`, and
`nu_eff` also includes the applicable configured limiter scattering. This
deliberately differs from the single collision frequency printed in Squire
et al. (2023), equation 2.7. The collisionless coefficient is unchanged.
The resulting unlimited flux is then capped with the equation 3.2 form,
`q = q_L*q_max/(q_max + abs(q_L))`.

The `cgl_lf_paper_oblique_wave` and `cgl_pure_paper_oblique_wave` inputs use
the Figure 17 background state from Squire et al. (2023), initialize a
transverse-velocity perturbation rather than a CGL eigenmode, and compare
against a linearized offline initial-value reference for that same perturbation.

The `cgl_*_paper_eigen_{alfven,slow,fast}` inputs initialize exact eigenvectors
of the same linearized Figure 17 CGL system. Regenerate those inputs with
`scripts/generate_cgl_lf_eigenmode_inputs.py`; the script computes the pure-CGL
or CGL-LF linear matrix, identifies the Alfvén, slow-like, and fast-like
positive-frequency branches, and writes the complex eigenvectors and
eigenvalues into AthenaK input parameters.
