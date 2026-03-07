# `scalar_mixing_test` Calibration Notes

The verifier in `vis/python/scalar_mixing_test_verify.py` was calibrated against fresh runs from the user pgen `src/pgen/scalar_mixing_test.cpp` with `scalar_only = true`. The analytic reference now evolves the exact AthenaK cell-centered initial condition, so the numerical and analytic fields match to roundoff at `t = 0`.

## Diffusive Final-Time Trend

2D reference family, `t = 1.0`, `kappa = 1/128`, `cfl = 0.4`

| Resolution | L1 | L2 | Linf | Centroid | Width Rel |
| --- | ---: | ---: | ---: | ---: | ---: |
| `64^2` | `9.266e-4` | `1.848e-3` | `1.000e-2` | `1.056e-3` | `8.061e-3` |
| `128^2` | `1.962e-4` | `3.804e-4` | `1.851e-3` | `2.763e-4` | `1.242e-3` |
| `256^2` | `4.648e-5` | `8.814e-5` | `4.209e-4` | `7.037e-5` | `1.836e-4` |

3D reference family, `t = 1.0`, `kappa = 1/128`, `cfl = 0.4`

| Resolution | L1 | L2 | Linf | Centroid | Width Rel |
| --- | ---: | ---: | ---: | ---: | ---: |
| `48^3` | `6.168e-4` | `1.907e-3` | `2.451e-2` | `6.214e-3` | `3.371e-2` |
| `64^3` | `1.929e-4` | `5.778e-4` | `7.605e-3` | `1.204e-3` | `8.085e-3` |
| `80^3` | `1.184e-4` | `3.445e-4` | `4.244e-3` | `8.057e-4` | `4.435e-3` |

The underresolved `32^3` and `40^3` 3D diffusive runs were numerically unstable and were not used for thresholds. The shipped regression thresholds are anchored to the stable `48^3` case.

## Pure-Advection Smoke

Pure advection was checked by overriding `hydro/scalar_diffusivity = 0.0` on the reference inputs.

- 2D at `cfl = 0.2`: `mass_rel_error = 1.20e-10`, `centroid_error = 1.65e-4`
- 3D at `cfl = 0.2`: `mass_rel_error = 7.65e-10`, `centroid_error = 8.12e-4`

At the more aggressive `cfl = 0.4`, the 3D pure-advection top-hat becomes numerically unstable, so the regression smoke keeps the pure-advection overrides at `cfl = 0.2`.

## Plotting

The verifier can also generate scalar-slice and convergence figures directly from binary dumps.

Time-evolution slices from a single run:

```bash
python3 vis/python/scalar_mixing_test_verify.py \
  --input inputs/tests/scalar_mixing_blob_diffusion_2d.athinput \
  --glob 'tmp/scalar_mixing_test_smoke/2d/bin/scalar_blob_diff_2d.hydro_w.*.bin' \
  --plot-slices \
  --plot-dir tmp/scalar_mixing_test_plots
```

Resolution/error convergence plots from multiple final dumps:

```bash
python3 vis/python/scalar_mixing_test_verify.py \
  --convergence-case inputs/tests/scalar_mixing_blob_diffusion_2d.athinput \
    'tmp/scalar_mixing_test_calibration/2d_64/bin/scalar_blob_cal_2d_64.hydro_w.*.bin' \
  --convergence-case inputs/tests/scalar_mixing_blob_diffusion_2d.athinput \
    'tmp/scalar_mixing_test_calibration/2d_128/bin/scalar_blob_cal_2d_128.hydro_w.*.bin' \
  --convergence-case inputs/tests/scalar_mixing_blob_diffusion_2d.athinput \
    'tmp/scalar_mixing_test_calibration/2d_256/bin/scalar_blob_cal_2d_256.hydro_w.*.bin' \
  --plot-dir tmp/scalar_mixing_test_plots
```
