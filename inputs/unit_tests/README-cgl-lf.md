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

The LF face fluxes use the finite-collision Squire et al. (2023) coefficients
when `<mhd>/nu_coll` or active mirror/firehose limiter scattering is present.
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

CGL instability thresholds are configured under `<mhd>`. In AthenaK units the
magnetic pressure is $B^2/2$ and $\Delta p = p_\perp-p_\parallel$. The positive
coefficients `firehose_threshold` (default 2) and `mirror_threshold` (default 1)
give soft thresholds $-\Lambda_{\rm FH}B^2/2$ and $+\Lambda_{\rm M}B^2/2$.
`firehose_backup_factor` (default 1) and `mirror_backup_factor` (default 2)
multiply these soft thresholds; the firehose backup wall is clipped at $-B^2$.
Thresholds must be finite and positive; backup factors must be finite and at
least 1. `limiter_backup_nu` (default $10^{10}$) is a finite, nonnegative LF
heat-flux suppression frequency in inverse code time, separate from pressure
projection. The fluid firehose wall at $-B^2$ applies independently of limiter
flags. Backup-wall diagnostics apply only when the effective backup policy is on.

The legacy `cgl_firehose_threshold = oblique|parallel` maps to numeric
`firehose_threshold = 1.4|2.0` when the numeric key is absent; conflicting explicit
values are rejected. `limiter_hardwall` retains its selected-soft-threshold
projection behavior. The threshold and backup factors do not enable a limiter.
