# CGL-LF Paper Validation Runbook

This runbook is for the dedicated `cgl_lf_paper` problem generator. Build it with:

```bash
cmake -S . -B build-cgl-lf-paper -DPROBLEM=cgl_lf_paper -DCMAKE_BUILD_TYPE=Release
cmake --build build-cgl-lf-paper -j
```

AthenaK uses `vA = B/sqrt(rho)`, so isotropic beta is `beta = 2*p/B^2`.
The Squire heat-flux caps remain
`q_parallel,max = sqrt(8/pi)*c_parallel*p_parallel` and
`q_perp,max = sqrt(2/pi)*c_parallel*p_perp`.

## Implemented Pgen Modes

- `mode = turbulence`: periodic `[1,1,2]` boxes with `B0` along `x3`.
- `mode = np_mode`: cold-electron analogue of the Majeski non-propagating mode.
- `mode = fast_wave`: perpendicular compressive fast-wave analogue.
- `mode = oblique_iaw`: oblique ion-acoustic-wave analogue.
- `mode = linear_wave_scan`: small-amplitude transverse wave seed for dispersion scans.

The pgen initializes `rho`, `u`, face-centered `B`, cell-centered `B`, `p_parallel`,
`p_perp`, and the CGL conserved state. The pressure default is
`p_parallel0 = p_perp0 = 0.5*beta0*B0^2`.

## Passive-Delta Controls

Use the existing CGL EOS switch:

```text
<mhd>
passive = true
iso_sound_speed = sqrt(0.5*beta0*B0^2/rho0)
```

In passive mode, `p_parallel` and `p_perp` still evolve as CGL/LF thermodynamic
variables, while the CGL Riemann solver removes anisotropic-pressure feedback from
the momentum flux.

## Diagnostics

Set `<problem>/user_hist = true`. The pgen writes a `.user.hst` file with:

- volume, mass, kinetic/magnetic/CGL thermal energies.
- `B^2` and `B^4` for `C_B2 = <B^4>/<B^2>^2 - 1`.
- volume-weighted beta, `Delta p`, and `abs(Delta p)`.
- ordinary and backup mirror/firehose limiter volume.
- mean effective collision rate from background plus active limiter scattering.
- heat-flux cap activity fractions for `|q_L|/qmax > 1` and `> 10`.
- a cell-centered heat-flux power proxy using the same capped LF coefficients.

Full derived field dumps are produced through existing `mhd_w_bcc` VTK output; the
offline script derives pressure anisotropy and beta from those primitive and
cell-centered magnetic fields when reduced dumps are converted to analysis arrays.

## Local Smoke

Use:

```bash
scripts/run_cgl_lf_paper_smoke.sh
```

The smoke script builds the custom pgen if needed, runs reduced-size active,
passive, limiter-disabled, NP-mode, and fast-wave cases, then writes
`summary.json` with `scripts/analyze_cgl_lf_paper.py`.

## Tiered Runs

Tier 1 local smoke:

- Use the shipped `inputs/cgl_lf_paper/*.athinput` files.
- Override to small grids for quick local checks, as the smoke script does.
- Acceptance: finite histories, no unexpected floors/NaNs, analysis completes.

Tier 2 nightly:

- `96x96x192`, `t_final = 3 L_perp/vA`.
- Run `beta0 = 1, 10, 100`, active/passive, Alfvénic and `driving_type=2`
  random forcing, hard-wall and finite-limiter scans.
- Acceptance: active CGL-LF reduces limiter occupancy relative to passive-delta,
  heat-flux cap fractions stay bounded, and energy residuals remain small.

Tier 3 paper grade:

- `192x192x384`, `t_final >= 10 L_perp/vA`.
- Squire turbulence: active CGL-LF, passive-delta, pure CGL, and MHD controls.
- Majeski turbulence: Alfvénic, random, and sonic-correlation-time forcing.
- Majeski waves: NP `beta=16,64`, `k_perp/k_parallel=4,8`, amplitudes
  `0.2-0.9`; fast-wave scans around `2/beta`; oblique IAW scans around `1/beta`.

Tier 4 convergence/HPC:

- `384x384x768`, selected Tier 3 cases only.
- Use sparse full dumps and frequent history/reduced diagnostics.

## Current Scope

The validation suite is an ion CGL-LF suite. Majeski comparisons that rely on
explicit isothermal electrons should be labelled cold-electron analogues until an
electron-pressure EOS extension is added.
