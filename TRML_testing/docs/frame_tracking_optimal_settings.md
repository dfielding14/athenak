# Optimal Frame-Tracking Configuration for TRML

## Recommended Production Configuration
The following configuration is the best compromise between interface centering, robustness, and compute efficiency across the completed campaign set:

- `ft_apply_every = 2`
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

## Why This Point Is Preferred
Two large campaign suites were completed and analyzed in detail:

1. `2026-03-02_ftpd_xi100_xi1000` (broad initial sweep)
2. `2026-03-03_ftpd_followup_refine` (targeted refinement near the robust region)

The follow-up suite outperformed the initial suite on reliability and controller behavior:

- Initial suite status: `43 ok`, `5 no_window_samples`
- Follow-up suite status: `48 ok`, `0 no_window_samples`
- Mean skip fraction reduced from `0.1168` to `0.0187`
- Mean cap fraction reduced from `0.0941` to `0.0121`

Within the follow-up suite, the strict quality winner (`c22`) was more aggressive and reached edge-contact earlier in interface-range diagnostics. The selected point (`c13`) had slightly weaker raw quality score but materially better conservative behavior with similar runtime cost.

Relative to follow-up baseline (`c00_baseline`), selected point (`c13`) gives:

- `quality_max`: `0.025455 -> 0.021376` (about `-16.0%`)
- `quality_mean`: `0.017243 -> 0.016342` (about `-5.2%`)
- `cost_mean`: `0.747725 -> 0.707248` (about `-5.4%`)
- conservative edge-contact metric `tmin_worst`: `15.7127 -> 23.9389` (about `+52.4%`, more conservative)

## Compute-Cost Guidance
`ft_apply_every=2` is retained deliberately.

The direct campaign evidence indicates that making updates less frequent (`3` or `4`) reduced quality strongly while yielding only modest runtime changes:

- `apply_every=3`: quality degradation roughly `+34%` to `+40%` on paired metrics
- `apply_every=4`: quality degradation roughly `+56%` to `+76%` on paired metrics
- runtime savings remained comparatively small

Given this, larger cadence values (e.g. `20`) are expected to be substantially less accurate for interface centering, especially at high `xi`.

---

## Appendix A: Test Matrix and Campaign Structure

### A1. Campaign 1: Broad Parameter Sweep
- Campaign: `2026-03-02_ftpd_xi100_xi1000`
- 48 total runs (`xi=100` and `1000`, 24 parameter variants each)
- Fixed physics:
  - `chi = 10^1.75` (`density_contrast=56.234132519`)
  - `tlim=150`
- Grid/decomposition:
  - `128 x 128 x 192`
  - meshblock `64 x 64 x 96` (8 meshblocks)
- Purpose: identify robust region around baseline PD control.

### A2. Campaign 2: Follow-Up Refinement
- Campaign: `2026-03-03_ftpd_followup_refine`
- 48 total runs (`xi=100` and `1000`, 24 targeted variants)
- Same physics/numerics as campaign 1
- Purpose: refine around high-performing region and improve cross-`xi` robustness.

### A3. Ranking Metrics
Ranking used a weighted quality functional over `t >= 20`:

- primary: `abs(mean(mean_z_ft))`
- penalties: skip fraction, cap fraction, and `std(mean_z_ft)`
- pairwise ranking across `xi=100` and `xi=1000`
- conservative cross-check from history diagnostics (`interface_z_ranges`) including edge-contact timing

### A4. Conservative Diagnostic Review
Interface-range diagnostics (`interface_z_ranges.png`) were reviewed for edge-reaching behavior. This review was used to break ties among near-best quantitative candidates and avoid overly aggressive settings that reached boundaries too quickly.

### A5. Key Output Artifacts
Primary evidence is in campaign analysis products, including:

- `ranking_combined.csv`
- `pairwise_robust_summary.csv`
- `pair_quality_conservative.csv`
- per-case diagnostics from `plot_run_diagnostics.py`

These are located under the campaign analysis directories in:
`/lustre/orion/ast207/proj-shared/dfielding/TRML/frame_tracking/campaigns/`.
