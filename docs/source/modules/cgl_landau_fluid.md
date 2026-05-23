# CGL Landau-Fluid Heat Flux

## Scope

AthenaK supports non-relativistic anisotropic MHD with the CGL equation of
state and an optional Landau-fluid (LF) heat-flux closure. The feature is
implemented as a distinct MHD parabolic process, `mhd/cgl_heat_flux`, rather
than as ordinary isotropic thermal conduction.

Use CGL with the HLLE solver:

```ini
<time>
sts_integrator = rkl2

<mhd>
eos = cgl
passive = false
rsolver = hlle
gamma = 1.666666666666667
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
lf_k_parallel = 32.0
lf_coefficient_mode = local
```

`lf_k_parallel` is the closure wavenumber magnitude, not a conductivity.
Set `lf_coefficient_mode = background` with a positive `lf_c_parallel0` to
hold the closure speed fixed while gradients and evolved variables remain
local.

## State And Update

CGL evolves density, momentum, total energy, and a sixth MHD state storing
conserved pressure anisotropy outside the LF sweep. Primitive arrays use
`IPR` for parallel pressure and `IPP` for perpendicular pressure.

The LF flux advances total energy and magnetic moment. During each RKL2 STS
sweep the MHD task graph performs:

1. Convert conserved anisotropy to magnetic moment.
2. Evaluate field-aligned LF energy and magnetic-moment fluxes.
3. Advance those two conserved quantities through all STS stages.
4. Refresh CGL primitives between stages in the magnetic-moment
   representation.
5. Convert magnetic moment back to conserved anisotropy at sweep completion.

This lifecycle prevents ordinary hyperbolic fluxes, output, and restart state
from interpreting magnetic moment as pressure anisotropy.

## Closure Controls

| Parameter | Default | Meaning |
| --- | --- | --- |
| `eos` | required | Set to `cgl` for this feature. |
| `passive` | `false` | Use passive isotropic-wave speeds when `true`. |
| `cgl_heat_flux` | absent | Set to `landau_fluid` to enable LF transport. |
| `cgl_heat_flux_integrator` | `sts` | Currently only `sts` is supported. |
| `lf_k_parallel` | required | Positive closure wavenumber magnitude. |
| `lf_coefficient_mode` | `local` | `local` or `background`. |
| `lf_c_parallel0` | required in background mode | Positive fixed parallel thermal speed. |
| `nu_coll` | `0.0` | Background anisotropy-relaxation frequency. |
| `mirror_limiter` | `false` | Enable mirror-limiter relaxation. |
| `firehose_limiter` | `false` | Enable firehose-limiter relaxation. |
| `limiter_nu_coll` | `0.0` | Limiter relaxation frequency. |
| `backup_limiters` | `false` | Enable hard backup bounds for active limiters. |

## Current Restrictions

- CGL is not available for SR, GR, or dynamical-GR MHD.
- CGL dynamic runs use `rsolver = hlle`; LLF and HLLD are rejected.
- Ordinary `<mhd>/conductivity` is rejected with `eos = cgl`.
- CGL LF is currently STS-only and cannot be combined with another MHD
  process selecting STS in the same run.

## Verification

Focused unit problems in `inputs/unit_tests/` exercise CGL transforms, CGL
FOFC, LF parallel and perpendicular decay, magnetic-field-gradient coupling,
flux limiting, limiter suppression, and a field-aligned wave. The LF
quantitative pgen is `src/pgen/unit_tests/cgl_lf_quantitative_test.cpp`.

Run the LF regression suite from `tst/`:

```bash
python run_tests.py cgl/cgl_landau_fluid \
    --cmake=-DPROBLEM=unit_tests/cgl_lf_quantitative_test
```

See also [Super Time Stepping](super_time_stepping.md) and
[Magnetohydrodynamics](mhd.md).
