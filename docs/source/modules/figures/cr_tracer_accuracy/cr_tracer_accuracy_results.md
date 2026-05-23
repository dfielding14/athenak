---
orphan: true
---

# CR Tracer Accuracy Results

- profile: `docs`
- git revision: `8c041244a5c6910fe2f22125ffea1fd7cdd96071`
- git worktree dirty during generation: `True`
- MPI: `mpiexec (Open MPI) 5.0.9`
- distributed documentation ranks: `8`

## uniform_gyro
- fitted slope: `2.000`
- figure: `docs/source/modules/figures/cr_tracer_accuracy/uniform_gyro_phase_convergence.png`

## uniform_amr_mpi
- figure: `docs/source/modules/figures/cr_tracer_accuracy/uniform_amr_mpi_activity.png`

## linear_gather
- figure: `docs/source/modules/figures/cr_tracer_accuracy/linear_gather_error.png`

## manufactured_gather
- fitted slope: `1.756`
- figure: `docs/source/modules/figures/cr_tracer_accuracy/manufactured_gather_convergence.png`

## smooth_orbit_reference
- figure: `docs/source/modules/figures/cr_tracer_accuracy/smooth_orbit_reference_error.png`

## magnetic_mirror
- figure: `docs/source/modules/figures/cr_tracer_accuracy/magnetic_mirror_moment.png`

## gradb_drift
- figure: `docs/source/modules/figures/cr_tracer_accuracy/gradb_drift_response.png`

## amr_boundary
- figure: `docs/source/modules/figures/cr_tracer_accuracy/amr_boundary_activity.png`

## mpi_decomposition
- figure: `docs/source/modules/figures/cr_tracer_accuracy/mpi_decomposition_delta.png`

## isotropic_ensemble
- figure: `docs/source/modules/figures/cr_tracer_accuracy/isotropic_ensemble_moments.png`

## frozen_turbulent
- figure: `docs/source/modules/figures/cr_tracer_accuracy/frozen_turbulent_moments.png`

## pitch_angle_decorrelation
- figure: `docs/source/modules/figures/cr_tracer_accuracy/pitch_angle_decorrelation.png`

## qualitative_figures
- uniform_gyro: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_uniform_gyro.png` (`t = 0` to `6.40`, `1` tracks, `1` ranks)
- uniform_amr_mpi: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_uniform_amr_mpi.png` (`t = 0` to `0.64`, `4` tracks, `8` ranks)
- linear_gather: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_linear_gather.png` (`t = 0` to `0.01`, `1` tracks, `1` ranks)
- manufactured_gather: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_manufactured_gather.png` (`t = 0` to `0.01`, `1` tracks, `1` ranks)
- smooth_orbit_reference: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_smooth_orbit_reference.png` (`t = 0` to `1.60`, `1` tracks, `1` ranks)
- magnetic_mirror: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_magnetic_mirror.png` (`t = 0` to `18.00`, `1` tracks, `1` ranks)
- gradb_drift: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_gradb_drift.png` (`t = 0` to `1.60`, `1` tracks, `1` ranks)
- amr_boundary: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_amr_boundary.png` (`t = 0` to `0.64`, `4` tracks, `8` ranks)
- mpi_decomposition: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_mpi_decomposition.png` (`t = 0` to `1.00`, `4` tracks, `8` ranks)
- isotropic_ensemble: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_isotropic_ensemble.png` (`t = 0` to `0.80`, `6` tracks, `1` ranks)
- frozen_turbulent: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_frozen_turbulent.png` (`t = 0` to `3.20`, `5` tracks, `8` ranks)
- pitch_angle_decorrelation: `docs/source/modules/figures/cr_tracer_accuracy/qualitative_pitch_angle_decorrelation.png` (`t = 0` to `1.28`, `3` tracks, `8` ranks)

## executions

| case | ranks | wall seconds |
|------|------:|-------------:|
| uniform_gyro_cfl_0.04 | 1 | 0.044 |
| uniform_gyro_cfl_0.02 | 1 | 0.049 |
| uniform_gyro_cfl_0.01 | 1 | 0.055 |
| uniform_amr_mpi | 8 | 0.229 |
| linear_gather_lin | 1 | 0.044 |
| linear_gather_trilinear | 1 | 0.047 |
| linear_gather_tsc | 1 | 0.043 |
| manufactured_gather_16 | 1 | 0.045 |
| manufactured_gather_32 | 1 | 0.056 |
| manufactured_gather_64 | 1 | 0.123 |
| manufactured_gather_128 | 1 | 0.628 |
| smooth_orbit_16 | 1 | 0.058 |
| smooth_orbit_32 | 1 | 0.166 |
| smooth_orbit_64 | 1 | 1.713 |
| smooth_orbit_128 | 1 | 26.223 |
| smooth_orbit_timestep_0.0125 | 1 | 0.518 |
| smooth_orbit_timestep_0.00625 | 1 | 0.924 |
| smooth_orbit_timestep_0.0015625 | 1 | 3.390 |
| magnetic_mirror | 1 | 9.023 |
| magnetic_mirror_passing_control | 1 | 3.033 |
| magnetic_mirror_convergence_16 | 1 | 0.105 |
| magnetic_mirror_convergence_32 | 1 | 0.375 |
| magnetic_mirror_convergence_64 | 1 | 2.585 |
| magnetic_mirror_convergence_128 | 1 | 19.412 |
| gradb_16 | 1 | 0.133 |
| gradb_32 | 1 | 0.507 |
| gradb_64 | 1 | 3.395 |
| gradb_128 | 1 | 32.713 |
| amr_boundary_16 | 8 | 0.399 |
| amr_boundary_32 | 8 | 0.756 |
| amr_boundary_64 | 8 | 4.959 |
| amr_boundary_uniform_128 | 8 | 9.565 |
| mpi_decomposition_n1 | 1 | 0.065 |
| mpi_decomposition_n2 | 2 | 0.087 |
| mpi_decomposition_n4 | 4 | 0.086 |
| mpi_decomposition_n8 | 8 | 0.111 |
| isotropic_ensemble | 1 | 0.049 |
| frozen_turbulent | 8 | 0.279 |
| pitch_angle_structured_32 | 8 | 0.294 |
| pitch_angle_structured_64 | 8 | 1.431 |
| pitch_angle_structured_128 | 8 | 19.639 |
| pitch_angle_uniform_control | 8 | 0.308 |
| qualitative_uniform_gyro | 1 | 0.145 |
| qualitative_uniform_amr_mpi | 8 | 0.223 |
| qualitative_linear_gather | 1 | 0.050 |
| qualitative_manufactured_gather | 1 | 0.052 |
| qualitative_smooth_orbit_reference | 1 | 0.136 |
| qualitative_magnetic_mirror | 1 | 9.049 |
| qualitative_gradb_drift | 1 | 0.253 |
| qualitative_amr_boundary | 8 | 0.215 |
| qualitative_mpi_decomposition | 8 | 0.121 |
| qualitative_isotropic_ensemble | 1 | 0.079 |
| qualitative_frozen_turbulent | 8 | 0.369 |
| qualitative_pitch_angle_decorrelation | 8 | 0.173 |
