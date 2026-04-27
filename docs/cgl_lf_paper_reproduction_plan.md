# CGL-LF Paper Reproduction and Stress-Test Plan

This plan targets the CGL-Landau-fluid method on the `CGL-STS-LF` branch after
adding the Squire heat-flux cap and finite-collision heat-flux coefficients.
The goals are to reproduce the paper results as closely as this fluid model can
support, then to push the implementation with deliberately difficult tests.

AthenaK uses magnetic fields normalized so that `vA = B/sqrt(rho)`. When
translating paper parameters that use `B/sqrt(4*pi*rho)`, use `beta = 2*p/B^2`
for an isotropic pressure in code units. The heat-flux coefficients and caps are
unchanged because they depend on `c_parallel`, pressure, `lf_k_parallel`, and
collision frequency, not directly on the magnetic pressure normalization.

## Local Runs Already Completed

The current quantitative pgen already contains a Squire Figure 17-style oblique
linear initial-value setup. It seeds a transverse velocity perturbation on the
Figure 17 background and compares the result with an offline integration of the
same linearized initial-value problem; it is not an eigenmode initialization.
I reran the CGL-LF and pure-CGL variants with validation CSV output under
`/tmp/athenak_cgl_lf_repro_data`.

| Input | Result |
| --- | --- |
| `inputs/unit_tests/cgl_lf_paper_oblique_wave.athinput` | Passed against the linear oblique IVP reference; `vy_rel_err=1.382402e-06`, `By_rel_err=2.483205e-06`, `p_parallel_rel_err=2.119553e-06`, `p_perp_rel_err=1.097882e-06`. |
| `inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput` | Passed against the pure-CGL linear IVP reference; `vy_rel_err=1.355613e-06`, `By_rel_err=2.483684e-06`, `p_parallel_rel_err=2.120868e-06`, `p_perp_rel_err=1.099899e-06`. |
| `inputs/unit_tests/cgl_pure_paper_eigen_{alfven,slow,fast}.athinput` | Exact pure-CGL eigenvectors for the Figure 17 background; all three passed with maximum tracked-component relative errors below `2.2e-6`. |
| `inputs/unit_tests/cgl_lf_paper_eigen_{alfven,slow,fast}.athinput` | Exact CGL-LF eigenvectors for the Figure 17 background; all three passed with maximum tracked-component relative errors below `2.9e-5`. |
| `inputs/unit_tests/cgl_lf_quant_parallel_collisional.athinput` | Passed finite-collision `q_parallel` decay; relative error `1.564782e-04`. |
| `inputs/unit_tests/cgl_lf_quant_perp_collisional.athinput` | Passed finite-collision `q_perp` decay; relative error `5.596518e-05`. |
| `inputs/unit_tests/cgl_lf_flux_limiter.athinput` | Passed large-gradient heat-flux cap check; unlimited flux reached `3.622256e+02 q_max`. |
| `inputs/unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput` | Passed limiter-collisional heat-flux suppression; capped collisionless unlimited flux was `6.312215e+01 q_max`, while limiter suppression reduced `q_parallel` to `3.506286e-03` and `q_perp` to `2.500361e-03` of the collisionless value. |

## Squire et al. 2023 Targets

### Linear CGL-LF Waves

Reproduce Appendix A/Figure 17 at higher fidelity than the current smoke input.

- Background: `rho0=1`, `p_parallel0=p_perp0=5`, `B0=(1,sqrt(2),0.5)`,
  `k=2*pi*xhat`, and `lf_k_parallel=|k|`.
- Runs: pure CGL, collisionless CGL-LF, `nu_coll=10`, and `nu_coll=1e10`.
- Modes: Alfvén, slow-like, and fast-like eigenvector inputs are implemented for
  pure CGL and CGL--LF.  Entropy-like and anisotropy-like eigenvectors remain
  useful follow-up branches for a full Figure 17 scan.
- Resolutions: `Nx=128,256,512,1024`.
- Acceptance: recover the complex frequency from the time series and show the
  expected first- to second-order convergence trends by mode. Strongly damped
  LF modes may converge closer to first order because the heat-flux operator is
  split and limited.

Implementation status: `src/pgen/unit_tests/cgl_lf_quantitative_test.cpp` now
initializes supplied complex eigenvectors and compares against their
\(\exp(\lambda t)\) evolution.  Regenerate the inputs with
`scripts/generate_cgl_lf_eigenmode_inputs.py`.

### Driven Alfvénic Turbulence

Reproduce the method-level Squire turbulence trends rather than a single bitwise
figure match.

- Domain: fully periodic `[Lx,Ly,Lz]=[1,1,2]` with `B0` along `z`.
- Resolution ladder: local smoke `48x48x96`, standard `192x192x384`, high
  `384x384x768` if resources allow.
- Parameters: `beta0=1,10,100`, `rho0=1`, isotropic initial pressure
  `p0=0.5*beta0*B0^2` in AthenaK units.
- Heat flux: `lf_k_parallel=4*pi/L_parallel`, local coefficient mode.
- Limiters: hard-wall reference with `limiter_nu_coll=1e10*vA/L_perp`, plus
  finite-limiter scans.
- Forcing: Ornstein-Uhlenbeck large-scale incompressible velocity forcing,
  perpendicular to the mean field, with `t_corr=L_parallel/vA`, power over the
  lowest few modes, and an energy-injection rate tuned to
  `delta B_perp/B0 ~= 0.5`.
- Duration: at least `10 L_perp/vA`; use restarts for longer saturated
  statistics.
- Controls: active CGL-LF, pure CGL without LF, ideal/isothermal MHD, and
  passive-delta if a passive anisotropy force switch is added.

Required diagnostics:

- Brazil plots of local `beta` versus pressure anisotropy.
- Active limiter volume fraction and time spent beyond mirror/firehose
  thresholds.
- `delta B^2` compressibility statistic `C_B2`.
- Kinetic/magnetic/pressure spectra; `p_parallel`, `p_perp`, `Delta p`, and
  `B^2` spectra.
- Field-parallel versus field-perpendicular gradient spectra for `u_parallel`
  and `Delta p`.
- Anisotropic-pressure work and heat-flux power, split into capped and uncapped
  regions.
- Cascade efficiency from injected power minus large-scale anisotropic/heat-flux
  losses.

## Majeski, Kunz, and Squire 2023 Targets

This paper is partly hybrid-kinetic. AthenaK CGL-LF cannot reproduce Larmor-scale
mirror/firehose fluctuations directly, so the reproducible target is the
large-scale CGL-LF envelope with a limiter-scattering surrogate for
microinstability feedback.

### Linear Magnetosonic Dispersion

- Target: Landau-fluid CGL-MHD dispersion as a function of
  `nu/(k_parallel*v_th,i)`.
- Parameters from the paper figure: `k_perp=4*|k_parallel|`, `beta_i0=16`,
  and `T_e=T_i0`.
- Runs: scan `nu_coll/(k_parallel*v_th,i)` from collisionless to Braginskii-like
  values; compare measured real frequency and damping rate for NP/entropy,
  anisotropy, oblique ion-acoustic, and fast branches.
- Caveat: if the current EOS has no explicit isothermal electron pressure, the
  exact `T_e=T_i0` curve needs an EOS extension or the comparison should be
  labelled as the cold-electron analogue.

### Non-Propagating Mode

- Geometry: 2D periodic domain elongated along `B0`, with one wavelength
  `lambda_parallel x lambda_perp`.
- Obliquity: `k_perp/k_parallel=4` and `8`.
- Amplitudes: `alpha=|delta B_parallel|/B0` below, near, and above the mirror
  threshold: `0.2,0.3,0.4,0.7,0.9`.
- Baseline beta: `beta_i0=16`; add `beta_i0=64` for asymptotic behavior.
- Fluid surrogate: use CGL-LF with mirror limiter enabled and scan
  `limiter_nu_coll`; compare with limiter disabled.
- Metrics: NP-mode Fourier amplitude versus time, pressure-anisotropy extrema,
  mirror-limiter occupancy localized to `delta B_parallel<0`, and decay
  transition across `alpha ~= 0.3`.

### Perpendicular Fast Wave

- Geometry: 1D perpendicular propagation with compressive `B` and density
  perturbations.
- Betas: `beta=16,100`.
- Amplitude scan: straddle `delta B/B0 ~= 2/beta`.
- Metrics: firehose and mirror limiter occupancy over phase, shock steepening
  time, effective adiabatic response, and convergence toward MHD-like dynamics
  when limiter scattering is strong.

### Oblique Ion-Acoustic Wave

- Geometry: 2D oblique mode with `k_perp >> k_parallel`.
- Amplitude scan: around `|delta B_parallel|/B0 ~= 1/beta`.
- Collisional scan: verify the non-propagating band near
  `nu/(k_parallel*v_th,i) in [2, 0.75*sqrt(beta)]` where applicable.
- Metrics: phase speed, damping rate, pressure anisotropy phase, active limiter
  fraction, and sign symmetry between mirror and firehose portions of the wave.

## Majeski, Kunz, and Squire 2024 Targets

This is the full high-beta self-organization turbulence problem. It should share
the Squire turbulence infrastructure but add the diagnostics emphasized in the
Majeski 2024 paper.

- Domain: `[1,1,2]`, `B0` along the long direction.
- Standard resolution: `n_perp=192`, `n_parallel=384`; low-resolution smoke
  `48x48x96`; convergence `96x96x192` and `384x384x768`.
- Betas: `beta0=1,10,100`.
- Forcing modes: Alfvénic incompressible forcing and random forcing; add sonic
  correlation-time runs to deliberately excite compressive fluctuations.
- Heat flux: `lf_k_parallel=4*pi/L_parallel`.
- Limiter scan: hard-wall `limiter_nu_coll`, finite `limiter_nu_coll`, and
  limiter disabled for controlled failures.
- Runtime: `t_final >= 10 L_perp/vA`.
- Acceptance: active-delta runs should suppress mirror/firehose volume fraction,
  suppress `grad_parallel u_parallel` and `grad_parallel Delta p`, preserve an
  MHD-like Alfvénic inertial range, and show heat-flux power becoming subleading
  when field lines self-organize toward isotherms.

## Demanding Physics Tests

These should become separate CI/smoke, nightly, and HPC-regression tiers.

1. **Limiter-threshold face scan**: sweep pressure anisotropy through ordinary
   and backup mirror/firehose thresholds on a single face, verifying `nu_eff`,
   heat-flux suppression, cap sign, and no discontinuous energy-path failure.
2. **Collision asymptotics**: for `nu_eff << c*k` and `nu_eff >> c*k`, verify
   `chi_perp` and `chi_parallel` approach the collisionless and Braginskii-like
   limits with both local and background coefficient modes.
3. **Rotated-field invariance**: run the same field-aligned temperature wave in
   `x`, `y`, `z`, and oblique directions; assert identical decay after
   projecting onto `bhat`.
4. **Extreme heat-flux cap**: use a discontinuous temperature profile with
   `|q_L|/q_max > 10^4`; assert bounded flux, sign preservation, positive
   pressures, and no STS overshoot.
5. **Low-B transition**: advect a temperature gradient through regions where
   `B` approaches `bfloor`; assert no NaNs and controlled fallback behavior.
6. **High-beta large-amplitude Alfvén wave**: use `beta=100`, amplitudes below
   and above the interruption threshold, and verify the effective Alfvén speed
   and pressure-anisotropy feedback.
7. **Perpendicular fast-wave shock**: measure shock formation delay across the
   `delta B/B0 ~= 2/beta` threshold with limiter scattering on/off.
8. **Oblique NP/IAW branch tracking**: initialize exact eigenvectors and verify
   the measured complex frequency across a two-dimensional `(nu,beta)` scan.
9. **STS robustness sweep**: scan `sts_max_dt_ratio`, resolution, and limiter
   activity for the strongest heat-flux gradients; compare against a small-step
   non-STS reference.
10. **Conservation budget**: for every nonlinear problem, track injected power,
    kinetic/magnetic/internal energy, anisotropic work, collisional heating, and
    heat-flux transport so the residual remains at truncation-error level.
11. **Active/passive-delta divergence**: add a passive-delta mode and show that
    active CGL-LF suppresses threshold occupancy relative to identical MHD-like
    fields.
12. **Limiter collisionality locality**: create adjacent stable/unstable regions
    with identical thermal gradients and verify only unstable faces receive
    limiter-suppressed heat flux.

## Implementation Slices

1. Extend the current exact-eigenvector inputs to entropy-like and
   anisotropy-like branches, then scan resolution.
2. Add a CGL-specific turbulence pgen or extend `turb.cpp` so it initializes
   `p_parallel`, `p_perp`, total energy, and CGL diagnostics correctly.
3. Add history diagnostics for anisotropy thresholds, heat-flux cap activity,
   limiter activity, `Delta p` work, and heat-flux power.
4. Add offline Python analysis for Fourier mode fitting, spectra, Brazil plots,
   transfer functions, and run summaries.
5. Add passive-delta and optional isothermal-electron support if exact Majeski
   comparisons require them.
6. Split runs into quick local smoke, nightly quantitative, and HPC paper
   reproduction tiers.
