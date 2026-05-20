# General Cooling Source Term Implementation Plan

## Purpose

Build a new, mergeable feature branch that replaces the current hard-wired ISM
and CGM cooling paths with a general `<cooling>` source-term block. The new
system should support table-driven cooling, analytic power-law cooling,
piecewise power-law cooling, reusable modifiers, optional cooling timestep
limits, and interval-integrated cooling history diagnostics.

This is a design and implementation plan. It is not the final user
documentation.

## Current State

Cooling currently lives inside `SourceTerms` and is selected with flags in the
fluid block:

- `<hydro>/ism_cooling`
- `<hydro>/cgm_cooling`
- `<hydro>` or `<mhd>/rel_cooling`

The ISM path uses `ISMCoolFn(T)` from `src/srcterms/ismcooling.hpp` plus a
uniform heating rate. The CGM path uses hard-coded tables in
`src/srcterms/cooling_tables.hpp`, metallicity from passive scalar zero, a
temperature ceiling, a spatially varying heating profile, and an approximate
self-shielding blend between PIE and CIE cooling. Cooling timestep limits are
duplicated in `src/srcterms/srcterms_newdt.cpp`.

TRML carries a separate pgen-level cooling/heating implementation with its own
timestep and history accumulation logic. That code provides useful patterns,
especially for reducing integrated cooling rates, but should not remain the
general mechanism.

## Design Goals

1. Make `<cooling>` the only supported interface for radiative cooling.
2. Remove legacy `<hydro>/ism_cooling` and `<hydro>/cgm_cooling` input modes.
3. Preserve ISM-like and CGM-like physics as named models inside `<cooling>`,
   not as legacy fluid-block flags.
4. Support cooling tables with one, two, or three dimensions.
5. In the first pass, allow table axes drawn only from temperature, density,
   and passive scalar or metallicity.
6. Exclude pressure axes and coordinate axes from the first pass.
7. Treat table values as cooling coefficients `Lambda`, not direct volumetric
   cooling rates.
8. Let users explicitly choose the density measure used in cooling and heating
   rates.
9. Make unit handling explicit for model axes, model values, and timestep
   bounds.
10. Support reusable cooling, heating, and shielding modifiers.
11. Provide a pgen hook for custom user-defined cooling/heating evaluators.
12. Accumulate gross and net cooling over every source update, then output
    interval-averaged rates in history output.

## Non-Goals For V1

- No pressure-dependent cooling tables.
- No coordinate-dependent cooling tables.
- No arbitrary-dimensional tables beyond three dimensions.
- No direct volumetric `Edot` cooling tables.
- No implicit composition assumptions when cgs number densities are requested.
- No hidden fallback to old `<hydro>/ism_cooling` or `<hydro>/cgm_cooling`
  parameters.
- No temperature-bin history in the first mergeable version unless the core
  implementation proves straightforward. The accumulator should be extensible
  enough to add up to five bins later.

## Input Interface

The canonical interface is a standalone block:

```ini
<cooling>
enabled = true

cooling_model = table                 # ism | cgm | table | powerlaw | piecewise_powerlaw | user
heating_model = constant              # none | constant | table | powerlaw | piecewise_powerlaw | user

units = cgs                           # default model unit system: cgs | code

cooling_density = hydrogen_number_density
heating_density = hydrogen_number_density

temperature_mu = 0.62
hydrogen_mass_fraction = 0.75

cooling_bounds = zero                 # zero | clamp | fatal
heating_bounds = zero                 # zero | clamp | fatal

timestep = true
timestep_factor = 0.1

timestep_use_temperature_bounds = true
timestep_temperature_units = cgs
timestep_temperature_min = 2.0e4
timestep_temperature_max = 1.0e9

timestep_use_density_bounds = false
timestep_density_units = cgs
timestep_density_kind = hydrogen_number_density

history = true
history_gross = true
history_net = true
```

The exact key names can change during implementation, but the semantic
requirements should not: each ambiguous quantity needs an explicit unit and
density convention.

## Legacy Input Policy

The old fluid-block cooling keys should not remain valid aliases.

If the parser sees any of these keys, it should fail early with a migration
message:

```text
<hydro>/ism_cooling
<hydro>/cgm_cooling
<mhd>/ism_cooling
<mhd>/cgm_cooling
```

The message should direct users to `<cooling>` and the relevant named model.
This keeps old inputs from silently running with different semantics.

## Rate Definitions

Table values are `Lambda`, never direct volumetric `Edot`.

The general source term is

```text
Edot_cool = q_cool^2 * Lambda
Edot_heat = q_heat   * Gamma
Edot_net  = Edot_cool - Edot_heat
```

where `q_cool` and `q_heat` are explicitly selected density measures.

Allowed density selectors for V1:

```text
mass_density
number_density
hydrogen_number_density
```

For `mass_density`, `Lambda` has units consistent with
`Edot = rho^2 Lambda`. In cgs, that means
`erg cm^3 g^-2 s^-1`.

For `number_density` or `hydrogen_number_density`, `Lambda` has the usual cgs
cooling-function units `erg cm^3 s^-1`.

For heating, `Gamma` is per selected heating density. In cgs:

- with mass density: `Gamma` has units `erg g^-1 s^-1`
- with number density: `Gamma` has units `erg s^-1` per particle
- with hydrogen number density: `Gamma` has units `erg s^-1` per hydrogen

## Units And Composition

The cooling system must not infer number density conventions from mass density
without explicit composition metadata.

For cgs evaluation:

```text
rho_cgs = rho_code * density_unit

number_density:
n = rho_cgs / (mu * m_u)

hydrogen_number_density:
n_H = X * rho_cgs / m_u
    = rho_cgs / (mu_H * m_u), where mu_H = 1 / X
```

Required composition parameters:

- `number_density` requires `mu`
- `hydrogen_number_density` requires either `hydrogen_mass_fraction` (`X`) or
  `mu_H`
- if both `X` and `mu_H` are provided, they must satisfy `mu_H = 1/X` within
  tolerance
- cgs temperature derived from internal energy requires `temperature_mu`

The implementation may allow a single `mu` to default `temperature_mu` and
`number_density` conversion, but this must be explicit in the docs and visible
in the parsed runtime state. Hydrogen density should always require `X` or
`mu_H`.

The unit contract should be centralized in one small helper layer. Source
kernels should receive already prepared conversion factors instead of repeating
unit algebra in multiple places.

## Table Format

V1 should use a narrow AthenaK-native table format instead of a general CSV or
HDF5 reader. The format should be easy to parse once at startup and cheap to
evaluate on device.

Recommendation: require rectilinear, uniformly spaced axes in the stored
coordinate. If a source dataset is irregular, users should preprocess it into
this format before running AthenaK.

Example:

```text
ATHENAK_COOLING_TABLE 1
ndim 3
value_kind Lambda
value_units cgs
value_scale log10
bounds zero
data_order c_row_major_last_axis_fastest

axis0 temperature cgs log10 2.0 9.0 352
axis1 density cgs log10 -8.0 0.0 81 hydrogen_number_density
axis2 scalar code linear 0.0 2.0 16 scalar_index=0 scalar_scale=raw

data
...
```

Rules:

- `ndim` must be `1`, `2`, or `3`.
- `value_kind` must be `Lambda`.
- each axis declares kind, units, scale, lower edge, upper edge, and number of
  points
- allowed axis kinds are `temperature`, `density`, and `scalar`
- allowed axis scales are `linear` and `log10`
- table values can be stored as `linear` or `log10`
- flattening is C row-major with the last axis fastest
- comments are allowed only in clearly documented positions or with a single
  comment character such as `#`
- the parser validates the exact number of values
- the parser validates finite values and positive dimensional sizes

Uniform axes allow device interpolation to compute indices directly instead of
looping over axis arrays.

## Interpolation And Bounds

Supported table dimensions:

- 1D: linear interpolation
- 2D: bilinear interpolation
- 3D: trilinear interpolation

Interpolation is in the table's stored coordinate. For example, a `log10`
temperature axis and `log10` value scale means the evaluator interpolates in
`log10(T)` and returns `10^interpolated_value`.

Out-of-range behavior is user selectable:

```text
zero
clamp
fatal
```

Default: `zero`.

`fatal` should not call `exit` inside a device kernel. The kernel should reduce
an out-of-range flag and enough diagnostic context to the host, then the host
should abort with a clear message after the kernel completes.

## Analytic Power-Law Models

V1 should support:

- single power law
- piecewise power law
- multi-piece piecewise power law

The first implementation should make temperature the independent variable for
analytic power laws unless there is a concrete need to support density or scalar
power laws immediately. If generalized analytic axes are added, they should use
the same axis-kind and unit machinery as tables.

Example:

```ini
cooling_model = piecewise_powerlaw
cooling_powerlaw_axis = temperature
cooling_powerlaw_axis_units = cgs
cooling_powerlaw_reference_temperature = 1.0e5
cooling_powerlaw_reference_lambda = 1.0e-22
cooling_powerlaw_breaks = 1.0e4,1.0e5,1.0e6
cooling_powerlaw_slopes = 2.0,-0.7,0.5,1.0
```

Heating should use the same infrastructure, with `Gamma` replacing `Lambda`.

## Named ISM And CGM Models

The old physical models can remain available only as `<cooling>` models:

```ini
<cooling>
enabled = true
cooling_model = ism
heating_model = constant
```

and

```ini
<cooling>
enabled = true
cooling_model = cgm
```

The implementation should avoid preserving old parameter names unless they are
still the clearest names in the new block.

The current CGM behavior includes several separable pieces:

- PIE table cooling
- CIE table cooling
- low-temperature analytic cooling
- metallicity scaling
- self-shielding blend between PIE and CIE
- heating profile
- temperature ceiling

The new design should express as many of these as reusable model components or
modifiers as practical.

## Modifiers

Cooling, heating, and shielding modifiers should be reusable, not trapped
inside the CGM model.

Conceptually:

```text
cooling = cooling_modifier * base_cooling_model
heating = heating_modifier * base_heating_model
shielding may modify cooling, heating, or both
```

The current CGM shielding is not a pure multiplicative modifier because it
blends PIE and CIE cooling:

```text
Lambda = (1 - f_shield) * Lambda_CIE + f_shield * Lambda_PIE
Gamma  = (1 - f_shield) * Gamma
```

So the modifier framework should support both:

- multiplicative modifiers
- blend modifiers between two cooling evaluators

V1 can start with the modifiers needed to reproduce the current CGM behavior,
then add more general variants after the base infrastructure is stable.

## Timestep Constraint

The timestep constraint lives in `<cooling>`, not in pgens.

Base formula:

```text
dt_cool = timestep_factor * e_int / abs(Edot_net)
```

Only cells passing the timestep-selection bounds participate.

Bounds must not use sentinel values such as `-1` to mean disabled. Use explicit
booleans:

```ini
timestep_use_temperature_bounds = true
timestep_temperature_units = cgs
timestep_temperature_min = 2.0e4
timestep_temperature_max = 1.0e9

timestep_use_density_bounds = false
timestep_density_units = cgs
timestep_density_kind = hydrogen_number_density
```

Temperature and density bound units are independent because a user may want, for
example, cgs temperature bounds and code-density bounds.

If cgs density bounds use `number_density` or `hydrogen_number_density`, the
same composition requirements apply as for cooling and heating rates.

## History Accumulation

When enabled, history output should include both gross and net cooling rates:

```text
cooling_gross_rate = integral_interval integral_volume Edot_cool dV dt / Delta t
cooling_net_rate   = integral_interval integral_volume (Edot_cool - Edot_heat) dV dt / Delta t
```

These are interval-integrated rates, not instantaneous samples at the history
output time.

Implementation requirements:

- accumulate during every cooling source update
- use the actual source substep weight `bdt`
- accumulate gross cooling energy and net cooling energy separately
- perform device reductions for pack-local contributions
- perform MPI reductions where needed
- divide by elapsed simulation time since the previous history output
- reset accumulators after history writes

Kokkos reductions are preferred over atomics for the global totals. Reductions
are cleaner and map naturally to the existing MPI reduction pattern.

Temperature-bin diagnostics can be added later by extending the accumulator to
hold a small fixed number of bins, likely no more than five.

## User-Defined Pgen Cooling Hook

The system should allow a problem generator to provide custom cooling/heating
rates instead of using the built-in table or power-law evaluators.

Because the evaluator runs inside Kokkos device kernels, this should not be a
raw host callback or virtual function. A portable design is to provide common
templated cooling kernels and let a pgen register a small host launcher that
instantiates those kernels with a pgen-defined device evaluator.

Sketch:

```cpp
struct CoolingCellState {
  Real rho_code;
  Real eint_code;
  Real temperature;
  Real cooling_density;
  Real heating_density;
  Real scalar0;
};

struct CoolingRates {
  Real lambda;
  Real gamma;
};

struct MyPgenCoolingEvaluator {
  KOKKOS_INLINE_FUNCTION
  CoolingRates operator()(const CoolingCellState &state) const {
    CoolingRates rates;
    // User-defined device-safe calculation.
    return rates;
  }
};
```

The core cooling module should provide shared templated entry points:

```cpp
ApplyCoolingWithEvaluator(evaluator, runtime_data, w0, u0, bdt);
CoolingNewDtWithEvaluator(evaluator, runtime_data, w0);
```

The pgen hook registers launch functions that call those templates with its
evaluator. This keeps the core implementation responsible for unit conversion,
bounds, timestep selection, and history accumulation, while avoiding device
function-pointer portability issues.

The first implementation should include a small spike or prototype to confirm
this pattern works with the enabled Kokkos backends.

## Proposed Code Organization

Add a focused cooling submodule under `src/srcterms`:

```text
src/srcterms/cooling.hpp
src/srcterms/cooling.cpp
src/srcterms/cooling_table.hpp
src/srcterms/cooling_table.cpp
src/srcterms/cooling_units.hpp
src/srcterms/cooling_evaluators.hpp
src/srcterms/cooling_history.hpp
```

Expected responsibilities:

- `cooling.hpp/cpp`: owns parsed runtime state and source-term entry points
- `cooling_table.*`: table file parsing and host/device table storage
- `cooling_units.hpp`: density, temperature, Lambda, Gamma, and Edot conversion
- `cooling_evaluators.hpp`: device-callable built-in evaluators and interpolation
- `cooling_history.hpp`: accumulators and history-output bridge

`SourceTerms` should own an optional cooling object and call it from the same
task slots where ISM/CGM cooling currently run. It should no longer own CGM
tables directly.

## Implementation Phases

### Phase 1: Branch And Scaffolding

1. Create a feature branch from current `master`.
2. Add the new cooling files.
3. Add `<cooling>` parsing with `enabled=false` default.
4. Add fatal checks for old fluid-block `ism_cooling` and `cgm_cooling`.
5. Wire an empty cooling object through `SourceTerms` without changing runtime
   behavior when `<cooling>` is absent.

### Phase 2: Unit And State Layer

1. Implement the explicit unit/conversion parser.
2. Implement composition validation.
3. Implement `CoolingCellState` construction from Hydro and MHD primitives.
4. Add validation for ambiguous cgs number-density requests.
5. Add clear fatal messages for missing `mu`, `temperature_mu`, `X`, or `mu_H`.

### Phase 3: Built-In Analytic Models

1. Implement constant heating.
2. Implement single power-law cooling/heating.
3. Implement piecewise and multi-piece power-law cooling/heating.
4. Add bounds handling.
5. Add small deterministic tests or problem setups that can verify expected
   energy changes in one zone or a uniform box.

### Phase 4: Table Loader And Interpolation

1. Implement the fixed table parser.
2. Store table metadata and values in device-accessible arrays.
3. Implement 1D, 2D, and 3D interpolation.
4. Implement `zero`, `clamp`, and `fatal` bounds behavior.
5. Add a minimal example table for testing.

### Phase 5: Source Update

1. Replace the old ISM/CGM source update paths with calls to the new cooling
   object.
2. Keep old ISM/CGM physics only as named `<cooling>` models.
3. Ensure Hydro and MHD both work.
4. Verify energy updates use primitive inputs and update conserved energy.
5. Confirm behavior with both internal-energy and temperature EOS paths.

### Phase 6: Timestep Constraint

1. Move cooling timestep logic into the new cooling object.
2. Implement explicit timestep bound units and density kinds.
3. Reduce minimum `dt_cool` over valid cells.
4. Make out-of-range fatal handling host-side.
5. Verify the result feeds into the existing global timestep reduction.

### Phase 7: History Accumulators

1. Add gross and net cooling accumulators.
2. Accumulate `Edot * bdt * dV` during source updates.
3. Add history columns when `<cooling>/history=true`.
4. Divide by elapsed time since the previous history write.
5. Reset accumulators after output.
6. Validate against a uniform-box analytic cooling/heating rate.

### Phase 8: Reusable Modifiers

1. Factor the current CGM shielding/heating pieces into reusable components.
2. Implement multiplicative modifiers first.
3. Implement blend modifiers only where needed to preserve CGM-like PIE/CIE
   behavior.
4. Keep modifier parsing conservative and explicit.

### Phase 9: Pgen User Hook

1. Add pgen registration fields for user cooling launcher functions.
2. Implement templated common kernels that accept a device evaluator.
3. Add a simple pgen example using `cooling_model=user`.
4. Confirm the pattern compiles on the target Kokkos backends.
5. Document the required evaluator signature and restrictions.

### Phase 10: Documentation And Cleanup

1. Update `docs/source/modules/srcterms.md`.
2. Add full `<cooling>` user documentation.
3. Add table-format documentation.
4. Add migration notes for old ISM and CGM inputs.
5. Remove obsolete `cooling_tables.hpp` ownership from `SourceTerms` once the
   new table path is active.

## Validation Plan

Minimum validation before merge:

- build succeeds with Hydro and MHD enabled
- run with no `<cooling>` block preserves baseline behavior
- old `<hydro>/ism_cooling` and `<hydro>/cgm_cooling` inputs fail with clear
  migration errors
- uniform-box constant-heating test matches analytic energy change
- uniform-box power-law cooling test matches analytic energy change for one
  timestep
- table interpolation test matches hand-computed 1D, 2D, and 3D values
- timestep bound selection ignores excluded cells
- history gross and net rates match interval-integrated source updates
- cgs `number_density` fails without `mu`
- cgs `hydrogen_number_density` fails without `X` or `mu_H`
- table out-of-range behavior is verified for `zero`, `clamp`, and `fatal`

## Open Decisions

1. Whether analytic piecewise power laws should support only temperature in V1,
   or all V1 axis kinds.
2. Whether the table file should be fully self-describing or whether some unit
   metadata should be supplied only by the input file. The recommendation above
   is to make the table self-describing and validate it against the input.
3. How exactly to expose ISM and CGM named models in the new input block while
   avoiding old parameter names that imply legacy compatibility.
4. Whether temperature-bin history should be included in the first branch or
   only reserved in the accumulator design.
