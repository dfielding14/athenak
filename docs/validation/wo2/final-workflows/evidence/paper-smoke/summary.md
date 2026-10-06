# CGL-LF paper-smoke workflow summary

- Status: **passed**
- Created UTC: `2026-10-06T23:14:21.568315+00:00`
- Git revision: `738810c8807ac25f7e871f486d36c0d93c65aca2`
- Dirty worktree at execution: `True`
- Executable: `/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/launch_final.py`

## Cases

| Case | Status | LF safety | Firehose threshold |
| --- | --- | --- | --- |
| `paper_smoke_active_alfvenic` | passed | clean | `2.0` |
| `paper_smoke_active_random` | passed | clean | `2.0` |
| `paper_smoke_passive_alfvenic` | passed | clean | `2.0` |

## Paper Smoke Diagnostics

- `paper_smoke_active_alfvenic`: forcing=`alfvenic_z_perpendicular`, perpendicular force squared=`3.285833e+01`, parallel force squared=`0.000000e+00`, mean beta=`9.999967e+00`, C_B2=`6.007781e-10`, mirror/firehose fractions=`0.000000e+00/0.000000e+00`, energy delta/work/residual=`6.400000e-03`/`6.400000e-03`/`-1.776357e-15`, passed.
- `paper_smoke_active_random`: forcing=`isotropic_random`, perpendicular force squared=`2.145116e+01`, parallel force squared=`1.178701e+01`, mean beta=`1.000004e+01`, C_B2=`4.435392e-06`, mirror/firehose fractions=`0.000000e+00/0.000000e+00`, energy delta/work/residual=`6.400000e-03`/`6.400000e-03`/`-5.329071e-15`, passed.
- `paper_smoke_passive_alfvenic`: forcing=`alfvenic_z_perpendicular`, perpendicular force squared=`1.280000e+02`, parallel force squared=`0.000000e+00`, mean beta=`9.999924e+00`, C_B2=`5.677829e-10`, mirror/firehose fractions=`0.000000e+00/0.000000e+00`, energy delta/work/residual=`6.353132e-03`/`6.400000e-03`/`-4.686815e-05`, passed.

## LF Diagnostics

- `paper_smoke_active_alfvenic`: stages=`12288`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `paper_smoke_active_random`: stages=`12288`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.
- `paper_smoke_passive_alfvenic`: stages=`6144`, mirror occupancy=`0.000000e+00`, firehose occupancy=`0.000000e+00`, unsafe counters=`0/0/0/0/0`.
  hard-wall projections: `0`.
  heat-flux cap face fractions (`>1`, `>10`): parallel=`0.000000e+00`/`0.000000e+00`, perpendicular=`0.000000e+00`/`0.000000e+00`.

Generated output belongs to this result bundle and is not a source artifact.
