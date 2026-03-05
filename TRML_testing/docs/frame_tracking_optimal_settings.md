# Optimal Frame-Tracking Configuration for TRML

## Executive Recommendation
Use the same conservative PD controller shape across resolutions, but scale only the cadence with resolution.

Core controller (all resolutions):
- `ft_mode = pd`
- `ft_tau_avg = 0.04`
- `ft_tau_relax = 2.6`
- `ft_tau_vel = 0.09`
- `ft_temp_lo = 0.050`
- `ft_temp_hi = 0.066`
- `ft_max_abs_boost = 0.012`
- `ft_max_boost_change = 0.002`
- `ft_tau_int = 12.0`
- `ft_int_max_abs = 0.0015`
- `ft_int_leak_tau = 8.0`
- `ft_int_unsat_only = true`

Resolution-dependent cadence:
- `dx = 1/128`  -> `ft_apply_every = 2`
- `dx = 1/256`  -> `ft_apply_every = 4`
- `dx = 1/512`  -> `ft_apply_every = 8`

## Domain Recommendation for Low-`xi`
For `xi \lesssim 10`, use a symmetric vertical domain
- `x3min = -1.0`
- `x3max = 1.0`

to increase headroom and reduce early top-boundary interaction.

This recommendation should be treated as default for low-`xi` production runs.

## Evidence Basis
This recommendation combines outcomes from:
1. `2026-03-02_ftpd_xi100_xi1000`
2. `2026-03-03_ftpd_followup_refine`
3. resolution-scaling follow-up runs (`dx=1/64` to `1/512`) with cadence scaling checks

Main findings:
- The conservative PD shape above remained the most robust option across `xi=10,100,1000` when combined with cadence scaling.
- Keeping `ft_apply_every` fixed while increasing resolution degraded robustness; scaling cadence with resolution stabilized behavior.
- Low-`xi` cases are disproportionately sensitive to vertical headroom; extending to `z \in [-1,1]` improves safety margin.

## Input Templates Added
The following templates encode the recommendations above with `T_peak = T_cutoff = 10^0.5` and `chi=10^1.75`:
- `inputs/hydro/TRML/TRML_frametracking_tophatcooling.athinput` (default 1/128)
- `inputs/hydro/TRML/TRML_frametracking_tophatcooling_dx128.athinput`
- `inputs/hydro/TRML/TRML_frametracking_tophatcooling_dx256.athinput`
- `inputs/hydro/TRML/TRML_frametracking_tophatcooling_dx512.athinput`
