# Module: Source Terms

This page documents source-term behavior present in the public AthenaK
implementation. Source-term objects are constructed from
`<hydro_srcterms>`, `<mhd_srcterms>`, or `<rad_srcterms>` blocks. Random
turbulence forcing is a separate `TurbulenceDriver` constructed when a
`<turb_driving>` block exists.

```{note}
CGM cooling tables, an `<initial_turb>` impulse, tiled turbulence controls, and
AMR-specific turbulence-resize tasks are not present in the public
implementation baseline. They are not supported interfaces documented here.
```

## Implementation Files

| File | Public responsibility |
| --- | --- |
| `src/srcterms/srcterms.hpp`, `srcterms.cpp` | Flags and kernels for fluid and radiation source terms. |
| `src/srcterms/srcterms_newdt.cpp` | Timestep constraints for implemented cooling terms. |
| `src/srcterms/ismcooling.hpp` | Temperature-dependent ISM cooling coefficient. |
| `src/srcterms/turb_driver.hpp`, `turb_driver.cpp` | Ornstein-Uhlenbeck random forcing driver. |

Shearing-box Coriolis/tidal updates are implemented in
`src/shearing_box/shearing_box_srcterms.cpp`, not by `SourceTerms`.

## Fluid And Radiation Source Blocks

The same fluid flags can be placed in `<hydro_srcterms>` or
`<mhd_srcterms>`, according to the evolved fluid module. Radiation beam
injection belongs in `<rad_srcterms>`.

| Flag | Block | Required parameters when enabled | Behavior |
| --- | --- | --- | --- |
| `const_accel = true` | fluid | `const_accel_val`, `const_accel_dir` (`1`, `2`, or `3`) | Adds uniform momentum acceleration; ideal-gas states also receive the corresponding energy work. |
| `ism_cooling = true` | fluid | `hrate` and a `<units>` block | Applies the public optically thin ISM heating/cooling energy term using unit conversion. |
| `rel_cooling = true` | fluid | `crate_rel`; optional `cpower_rel` (default `1.0`) | Applies the implemented relativistic momentum and energy loss term. |
| `rad_beam = true` | radiation | `dii_dt`, `pos_1`, `pos_2`, `pos_3`, `dir_1`, `dir_2`, `dir_3`, `width`, `spread` | Injects radiation intensity in the beam mask. |

Constant acceleration is exercised by
`inputs/hydro/rt2d.athinput` and `inputs/mhd/rt2d-mhd.athinput`.
Radiation beam source-term configuration is present in
`inputs/radiation/bh_beam.athinput`.

### Cooling And Timestep Control

`ISMCoolFn()` combines the Koyama and Inutsuka analytic low-temperature
expression, interpolation through the bundled SPEX cooling table, and a
high-temperature power-law tail. `SourceTerms::NewTimeStep()` computes a
cellwise energy-over-net-cooling limit when `ism_cooling` is enabled. It also
computes an energy-over-cooling limit using `crate_rel` and `cpower_rel` when
`rel_cooling` is enabled. No public `cooling_dt_factor`, start-time flag, CGM
table flag, or metallicity parameter is parsed by this implementation.

### Example Constant Acceleration Block

```ini
<hydro_srcterms>
const_accel     = true
const_accel_val = -0.1
const_accel_dir = 2
```

## Turbulence Driver

When `<turb_driving>` is present, `MeshBlockPack::AddPhysics()` creates one
`TurbulenceDriver`. `InitializeModes` and one `AddForcing` call are registered
in `before_timeintegrator`, which the driver executes once per cycle before
the stage loop. The driver also inserts `AddForcing` into the per-stage task
sequence for hydro, MHD, or ion-neutral evolution.

The implemented driver evolves a random force through an
Ornstein-Uhlenbeck update. For `tcorr <= 1e-6` it uses the white-noise limit;
otherwise the correlation factors are
`exp(-dt/tcorr)` and `sqrt(1 - exp(-2*dt/tcorr))`. The normalization uses
`dedt` as the requested energy-injection rate.

### Public `<turb_driving>` Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `nlow` | `1` | Lower integer mode magnitude included in the forcing catalogue. |
| `nhigh` | `2` | Upper integer mode magnitude included in the forcing catalogue. |
| `driving_type` | `0` | Mode-selection branch used by the implemented driver (`0` isotropic; `1` uses perpendicular/parallel exponents). |
| `expo` | `5.0/3.0` | Spectral exponent for `driving_type = 0`. |
| `exp_prp` | `5.0/3.0` | Perpendicular spectral exponent for `driving_type = 1`. |
| `exp_prl` | `0.0` | Parallel spectral exponent for `driving_type = 1`. |
| `dedt` | `0.0` | Energy injection rate used when normalizing the force. |
| `tcorr` | `0.0` | Correlation time; values at or below `1e-6` select white noise. |

The shipped turbulence input `inputs/hydro/turb.athinput` supplies the runtime
configuration. Select its problem-specific source file
`src/pgen/turb.cpp` at build time with:

```bash
cmake -S . -B build-turb -DPROBLEM=turb
cmake --build build-turb
./build-turb/src/athena -i inputs/hydro/turb.athinput -d run-turb time/nlim=1
```

Remove the `time/nlim=1` override for the full supplied calculation.

## Related Pages

- [Hydrodynamics](hydro.md) and [MHD](mhd.md) describe the fluid modules that
  call fluid source terms.
- [Shearing Box](shearing_box.md) describes orbital and tidal updates kept
  outside this module.
- [Development Record: CGM Cooling Flow](../cgm_cooling_flow_metals.md)
  records excluded feature documentation without presenting it as supported.
