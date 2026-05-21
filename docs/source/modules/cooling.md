# Cooling Source Terms

Cooling and heating source terms are configured through a standalone `<cooling>`
block. Legacy `<hydro>/ism_cooling`, `<mhd>/ism_cooling`,
`<hydro_srcterms>/ism_cooling`, `<mhd_srcterms>/ism_cooling`, and the
corresponding `cgm_cooling` parameters are not supported.

## Basic block

```text
<cooling>
enabled = true
cooling_model = table
heating_model = none
units = code
cooling_density = mass_density
cooling_table = path/to/table.tbl
timestep = true
timestep_factor = 0.1
history = true
```

`cooling_model` may be `none`, `ism`, `cgm`, `table`, `powerlaw`,
`piecewise_powerlaw`, or `user`. `heating_model` may be `none`, `constant`,
`table`, `powerlaw`, `piecewise_powerlaw`, or `user`.

Cooling is applied as

```text
edot_cool = q_cool^2 Lambda
```

and heating is applied as

```text
edot_heat = q_heat Gamma
```

where `q_cool` and `q_heat` are selected by `cooling_density` and
`heating_density`.

## Units And Density

`units` sets the default value units for built-in analytic cooling and heating:

```text
units = code | cgs
```

Tables are self-describing and carry their own `value_units` and axis units.
The final volumetric source term is converted to code units before updating the
energy equation.

Density selectors are explicit:

```text
cooling_density = mass_density | number_density | hydrogen_number_density
heating_density = mass_density | number_density | hydrogen_number_density
```

For cgs number densities the `<cooling>` block must provide the composition
needed by that conversion:

```text
mu = 0.62                         # required for number_density
hydrogen_mass_fraction = 0.75     # required for hydrogen_number_density
# or
mu_H = 1.3333333333333333         # equivalent to 1 / hydrogen_mass_fraction
```

For cgs temperature axes or models, provide either `temperature_mu` or `mu`.

## Table Format

Tables are fixed-format ASCII files. Irregular source data should be
preprocessed into this format before runtime.

```text
ATHENAK_COOLING_TABLE 1
ndim 2
value_kind Lambda
value_units cgs
value_scale log10
bounds zero
data_order c_row_major_last_axis_fastest
axis0 temperature cgs log10 2.0 8.0 121
axis1 density cgs log10 -8.0 2.0 101 density_kind=hydrogen_number_density
data
...
```

Table interpolation is logarithmic along each axis by default. An `axisN` line
therefore has the form

```text
axisN kind units [linear|log10] xmin xmax n [axis_options]
```

If the optional scale token is omitted, the axis is treated as `log10`. For a
logarithmic axis, `xmin` and `xmax` are the lower and upper `log10` coordinates,
not the raw physical or code-unit values. Add `linear` on any axis that should
use linear coordinates; axes may be mixed within one table.

Required header keys:

- `ndim`: 1, 2, or 3.
- `value_kind`: `Lambda` for cooling tables, `Gamma` for heating tables,
  `modifier` for dimensionless modifier tables.
- `value_units`: `code` or `cgs`.
- `value_scale`: `linear` or `log10`.
- `bounds`: `zero`, `clamp`, or `fatal`; default behavior is zeroing
  out-of-range values.
- `data_order`: must be `c_row_major_last_axis_fastest`.
- `axisN`: one line per axis.

Allowed axes are `temperature`, `density`, and `scalar`. Density axes may add
`density_kind=mass_density`, `density_kind=number_density`, or
`density_kind=hydrogen_number_density`. Scalar axes may add `scalar_index=N`.

## Power Laws

Power-law cooling/heating uses

```text
cooling_model = powerlaw
cooling_powerlaw_axis = temperature
cooling_powerlaw_axis_units = code
cooling_powerlaw_value_units = code
cooling_reference_axis = 1.0
cooling_reference_value = 1.0e-2
cooling_slope = -1.0
```

Piecewise power laws use ordered break points and one more slope than break:

```text
cooling_model = piecewise_powerlaw
cooling_breaks = 1.0e4, 1.0e6
cooling_slopes = 2.0, -0.7, 0.5
```

The curve is continuous across breaks and normalized by
`cooling_reference_axis` and `cooling_reference_value`.

The same key pattern is available with the `heating_` prefix.

## CGM Mode And Modifiers

`cooling_model = cgm` blends separate PIE and CIE Lambda tables:

```text
cooling_model = cgm
cgm_pie_table = path/to/pie.tbl
cgm_cie_table = path/to/cie.tbl
shielding_model = cgm_approx
```

With `shielding_model = cgm_approx`, the PIE/CIE blend fraction is estimated
from hydrogen number density, temperature, cell length, and a cross section.
The tunable parameters are:

```text
shielding_density = hydrogen_number_density
shielding_cross_section_cgs = 1.0e-17
shielding_transition_temperature_cgs = 8.0e3
shielding_transition_width_cgs = 1.5e3
shielding_length_code = 0.0
shielding_apply_to_heating = true
```

`shielding_length_code = 0.0` means use the local x1 cell width.

Cooling and heating can also be multiplied by reusable modifiers:

```text
cooling_modifier_model = constant
cooling_modifier = 0.5

heating_modifier_model = powerlaw
heating_modifier_powerlaw_axis = temperature
heating_modifier_reference_axis = 1.0
heating_modifier_reference_value = 1.0
heating_modifier_slope = -0.5
```

Modifier models may be `none`, `constant`, `table`, `powerlaw`, or
`piecewise_powerlaw`.

## Timestep Constraint

Cooling timestep control is opt-in:

```text
timestep = true
timestep_factor = 0.1
```

Temperature and density filters require explicit units and are enabled by
separate booleans:

```text
timestep_use_temperature_bounds = true
timestep_temperature_units = cgs
timestep_temperature_min = 2.0e4
timestep_temperature_max = 1.0e9

timestep_use_density_bounds = true
timestep_density_units = cgs
timestep_density_kind = hydrogen_number_density
timestep_density_min = 1.0e-6
timestep_density_max = 1.0e2
```

There are no sentinel values such as `-1` for disabled bounds.

## History Output

When `history = true`, the standard Hydro or MHD history output gets interval
averaged cooling rates:

```text
history = true
history_gross = true
history_net = true
```

`cool_gross` is the gross radiative cooling rate integrated over space and over
the interval since the previous history output, then divided by that elapsed
simulation time. `cool_net` is the corresponding cooling-minus-heating rate.

## User Cooling Hook

`cooling_model = user` or `heating_model = user` delegates the whole cooling
source update to problem-generator launcher hooks. In `UserProblem()`, set

```cpp
user_cooling_func = MyCoolingSource;
user_cooling_timestep_func = MyCoolingTimeStep;  // required only if timestep=true
```

The source hook signature is declared in `src/pgen/pgen.hpp`:

```cpp
void MyCoolingSource(MeshBlockPack *pmbp, const DvceArray5D<Real> &w0,
                     const EOS_Data &eos_data,
                     const cooling::RuntimeData &runtime, const Real bdt,
                     const Real history_bdt,
                     DvceArray5D<Real> &u0,
                     Real &gross_energy, Real &net_energy);
```

The hook should launch its own device work and update `u0` with `bdt`.
It should return domain-local gross and net cooling energies in code units
using `history_bdt`.  `bdt` and `history_bdt` differ for multistage RK
integrators because earlier source increments are blended before they reach the
final step solution.

Most pgens should include `src/srcterms/cooling_hooks.hpp` and use its shared
templated launchers
instead of duplicating the unit conversion and history bookkeeping:

```cpp
cooling::ApplyCoolingWithEvaluator(pmbp, w0, eos_data, runtime,
                                   MyCoolingEvaluator{}, bdt, history_bdt, u0,
                                   gross_energy, net_energy);

cooling::CoolingNewDtWithEvaluator(pmbp, w0, eos_data, runtime, timestep,
                                   MyCoolingEvaluator{}, dtnew);
```

The timestep hook receives a `cooling::TimestepData` object so the shared helper
can apply the same temperature and density selection bounds as the built-in
models.

## Cooling Regression Test Suite

The dedicated cooling regression suite lives in `tst/test_suite/cooling/` and
uses the built-in `cooling_test` problem generator in
`src/pgen/tests/cooling_test.cpp`.  The pgen initializes a uniform 1D ideal-gas
state with zero velocity, optional passive scalar zero, and optional uniform
magnetic field.  There are no fluid gradients, so any energy change and any
cooling history signal must come from the cooling source term itself.

This pgen is intentionally not a physical flow test.  It is a source-term test
surface with analytic expectations:

- `density = 1`, `pressure = 1`, and `gamma = 5/3` give code temperature
  `T = p/rho = 1`.
- The domain volume is one, so a code-unit volumetric rate is also the
  domain-integrated rate.
- With `cooling_density = mass_density` and `rho = 1`, a table or power law
  returning `Lambda = 1.0e-2` should produce `cool_gross = 1.0e-2`.
- If heating returns `Gamma = 2.0e-3`, `cool_net` should be `8.0e-3`.
- The difference in total energy over the first evolved history interval must
  match the reported `cool_net` rate.

The suite covers:

- Hydro and MHD source-term paths.
- Constant power-law cooling under `rk1`, `rk2`, and `rk3`.
- 1D, 2D, and 3D cooling tables.
- Default logarithmic table axes and mixed logarithmic/linear table axes.
- Piecewise power-law cooling.
- Constant and table heating.
- Constant and table multiplicative cooling modifiers.
- User-defined pgen cooling hooks through `cooling_hooks.hpp`.
- CGM PIE/CIE table blending with shielding disabled.
- CGM shielding with explicit `mu` and `hydrogen_mass_fraction`.
- Cooling history disabled mode.
- Table bounds behavior for `zero` and `clamp`.
- Cooling timestep control and temperature-bound exclusion.
- Clean failures for removed legacy cooling keys, fatal table bounds, missing
  cgs `mu`, missing cgs hydrogen fraction, and cgs timestep-temperature bounds
  without a temperature composition.

During the review that produced the suite, the multistage RK history
accumulation was found to be overcounting.  The source update must still use the
stage coefficient `beta[stage-1]*dt`, but the history accumulator must use the
weight with which that stage contributes to the final RK solution.  For example,
`rk2` contributes one half of each stage to the final solution, not raw source
substeps with weights `1` and `1/2`.  The driver now passes a separate
`history_bdt` to cooling source launchers, and the test
`test_constant_powerlaw_history_weights` verifies that `rk1`, `rk2`, and `rk3`
all report the same interval-integrated cooling rate for a constant source.

To run the suite:

```bash
cd tst
python run_test_suite.py --test test_suite/cooling/test_cooling_cpu.py --cpu
```

Local result for the implementation documented here:

```text
collected 26 items
../../test_suite/cooling/test_cooling_cpu.py ..........................  [100%]
26 passed in 0.46s
```

AthenaK may write a final duplicate-time history row after the accumulator has
been reset.  The regression tests read raw history files and check the first
strictly evolved interval, which is the interval with the source contribution.
