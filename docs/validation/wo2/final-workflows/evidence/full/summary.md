# CGL-LF full workflow summary

- Status: **passed**
- Created UTC: `2026-10-06T23:14:06.780769+00:00`
- Git revision: `738810c8807ac25f7e871f486d36c0d93c65aca2`
- Dirty worktree at execution: `True`
- Executable: `/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/launch_final.py`

## Cases

| Case | Status | LF safety | Firehose threshold |
| --- | --- | --- | --- |
| `cgl_lf_smr_decay_2d` | passed | clean | `2.0` |
| `cgl_lf_oblique_decay_2d` | passed | clean | `2.0` |
| `cgl_lf_timestep_refresh` | passed | clean | `2.0` |
| `cgl_lf_density_contact` | passed | clean | `2.0` |
| `cgl_lf_uniform_timestep` | passed | clean | `2.0` |
| `cgl_lf_hotspot` | passed | clean | `2.0` |
| `cgl_lf_field_reversal_1d` | passed | clean | `2.0` |
| `cgl_lf_field_reversal_2d` | passed | clean | `2.0` |
| `cgl_reconstruction_ppmx` | passed | not applicable | `n/a` |
| `cgl_reconstruction_wenoz` | passed | not applicable | `n/a` |
| `cgl_collision_once` | passed | not applicable | `n/a` |
| `cgl_lf_quant_parallel` | passed | clean | `2.0` |
| `cgl_lf_quant_parallel_collisional` | passed | clean | `2.0` |
| `cgl_lf_quant_perp` | passed | clean | `2.0` |
| `cgl_lf_quant_perp_collisional` | passed | clean | `2.0` |
| `cgl_lf_quant_grad_b` | passed | clean | `2.0` |
| `cgl_lf_flux_limiter` | passed | clean | `2.0` |
| `cgl_lf_limiter_heat_flux_suppression` | passed | clean | `2.0` |
| `cgl_lf_limiter_mirror` | passed | clean | `2.0` |
| `cgl_lf_limiter_firehose` | passed | clean | `1.4` |
| `cgl_lf_field_wave` | passed | clean | `2.0` |
| `cgl_lf_paper_oblique_wave` | passed | clean | `2.0` |
| `cgl_pure_paper_oblique_wave` | passed | not applicable | `n/a` |
| `cgl_pure_paper_eigen_alfven` | passed | not applicable | `n/a` |
| `cgl_pure_paper_eigen_slow` | passed | not applicable | `n/a` |
| `cgl_pure_paper_eigen_fast` | passed | not applicable | `n/a` |
| `cgl_lf_paper_eigen_alfven` | passed | clean | `2.0` |
| `cgl_lf_paper_eigen_slow` | passed | clean | `2.0` |
| `cgl_lf_paper_eigen_fast` | passed | clean | `2.0` |

## LF Diagnostics

- `cgl_lf_smr_decay_2d`: stages=`1863680`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_oblique_decay_2d`: stages=`3940352`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_timestep_refresh`: stages=`1280`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_density_contact`: stages=`2688`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_uniform_timestep`: stages=`0`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
- `cgl_lf_hotspot`: stages=`25116672`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`3.029876e-02`/`0.000000e+00`, perpendicular=`4.084426e-02`/`0.000000e+00`.
- `cgl_lf_field_reversal_1d`: stages=`17920`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_field_reversal_2d`: stages=`286720`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_quant_parallel`: stages=`118528`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_quant_parallel_collisional`: stages=`118528`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_quant_perp`: stages=`236288`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_quant_perp_collisional`: stages=`369408`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_quant_grad_b`: stages=`384`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_flux_limiter`: stages=`896`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`9.687500e-01`/`9.687500e-01`, perpendicular=`9.687500e-01`/`9.687500e-01`.
- `cgl_lf_limiter_heat_flux_suppression`: stages=`896`, mirror occupancy=`1.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_limiter_mirror`: stages=`640`, mirror occupancy=`1.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_limiter_firehose`: stages=`640`, mirror occupancy=`0.000000e+00`, firehose occupancy=`1.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_field_wave`: stages=`1152`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_paper_oblique_wave`: stages=`4608`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_paper_eigen_alfven`: stages=`18432`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_paper_eigen_slow`: stages=`18432`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `cgl_lf_paper_eigen_fast`: stages=`18432`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.

Generated output belongs to this result bundle and is not a source artifact.
