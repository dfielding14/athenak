# External Gravity: Implementation Plan and User Documentation

This document is written as a MyST Markdown page that can be integrated directly into
the GitHub Pages documentation source tree.  The intended final location is:

```text
docs/source/modules/external_gravity.md
```

To expose it in the current Pages documentation, add a row for **External Gravity** to
`docs/source/modules/index.md` near the existing **Source Terms** entry, and optionally
link to it from `docs/source/modules/srcterms.md`.

## Purpose

The external-gravity feature should let users turn on a fixed gravitational potential
with a single input block:

```ini
<external_gravity>
model = point_mass
mass  = 1.0
G     = 1.0
```

The feature belongs in `src/srcterms` because it is a source term on the fluid
conserved variables, not a problem-generator-specific force.  Problem generators should
only be needed when a user wants custom initial conditions or a custom compiled-in
potential.

The final feature should provide:

- simple built-in potentials for common Cartesian, cylindrical, and spherical cases;
- a clear `pgen.cpp` hook for user-defined potentials;
- robust hydro/MHD source-term updates;
- `grav-PE` in the default history output when external gravity is enabled;
- gravity-aware boundary conditions for hydrostatic or force-balanced states;
- documented examples and regression tests covering the supported usage.

## Current Branch Status

The current `feature/external-gravity` branch has the first implementation pass:

| Area | Current state | Needed before calling this complete |
|------|---------------|--------------------------------------|
| Input block | `<external_gravity>` is accepted and can enable hydro/MHD source terms without `<hydro_srcterms>` or `<mhd_srcterms>`. | Keep, but make parsing more model-specific and less noisy. |
| Built-in potentials | Uniform, Cartesian harmonic, spherical harmonic, point mass/Plummer, NFW, logarithmic halo, Miyamoto-Nagai, outer CGM, gotham composite, and user hook. | Add analytic accelerations and model-specific tests. |
| Force application | Source term computes acceleration by face-centered finite differences of the potential. | Add analytic acceleration path for built-ins and keep finite differencing as fallback. |
| History | Hydro/MHD history gains `grav-PE` when external gravity is enabled. | Test passive scalars, MHD, MPI, AMR, and restart behavior. |
| Documentation | Short `src/srcterms/README.md` exists. | Replace/augment with this full Pages-ready page and example pages. |
| Examples | One regression input exists in `tst/inputs/external_gravity.athinput`. | Add user-facing example inputs and useful pgens. |
| Boundary conditions | Not implemented. | Must add gravity-aware boundary conditions before hydrostatic examples are reliable. |

## Runtime Architecture

### Source Files

| File | Role |
|------|------|
| `src/srcterms/external_gravity.hpp` | Defines the potential model enum, parameter struct, device-callable potential function, and user-potential hook declaration. |
| `src/srcterms/external_gravity.cpp` | Parses `<external_gravity>` and applies the source term to fluid conserved variables. |
| `src/srcterms/srcterms.hpp/cpp` | Owns the external-gravity flag and stores parsed potential data in `SourceTerms`. |
| `src/hydro/hydro.cpp` | Constructs `SourceTerms` when either `<hydro_srcterms>` or `<external_gravity>` exists. |
| `src/mhd/mhd.cpp` | Same construction logic for MHD. |
| `src/outputs/history.cpp` | Adds `grav-PE = ∫ rho phi dV` to default hydro/MHD history output. |
| `src/pgen/pgen.cpp` | Defines `external_gravity::UserExternalGravityPotential(x1, x2, x3)` for `model=user`. |
| `src/parameter_input.cpp` | Registers `<external_gravity>` as a valid input block. |

### Execution Flow

1. `ParameterInput::CheckBlockNames()` accepts `<external_gravity>`.
2. `Hydro` or `MHD` constructors create a `SourceTerms` object if the block exists.
3. `SourceTerms` parses the block once and stores an `ExternalGravityData` object.
4. During each hydro/MHD source-term task, `SourceTerms::ExternalGravity()` evaluates
   the potential at cell faces, computes

   ```math
   \mathbf{a} = -\nabla \Phi ,
   ```

   and applies

   ```math
   \Delta \mathbf{m} = \Delta t\,\rho\,\mathbf{a}.
   ```

5. For ideal-gas fluids, the source also applies the kinetic-work term

   ```math
   \Delta E = \Delta t\,\rho\,\mathbf{a}\cdot\mathbf{v}.
   ```

6. History output adds

   ```math
   E_{\rm grav} = \int \rho \Phi\,dV
   ```

   as `grav-PE`.

## User Interface

### Minimal Input

```ini
<external_gravity>
model = uniform
g1    = 0.0
g2    = 0.0
g3    = -1.0
```

No `<hydro_srcterms>` or `<mhd_srcterms>` block is needed.  The feature should work for
hydro and MHD in Newtonian/non-relativistic mode.

### Coordinate Convention

All models are expressed in the mesh coordinate variables `x1`, `x2`, and `x3`.
For models using cylindrical or spherical coordinates, define

```math
x = x_1 - x_{1,0}, \qquad y = x_2 - x_{2,0}, \qquad z = x_3 - x_{3,0},
```

```math
R = \sqrt{x^2 + y^2}, \qquad r = \sqrt{x^2 + y^2 + z^2}.
```

The origin is controlled by:

```ini
x1_origin = 0.0
x2_origin = 0.0
x3_origin = 0.0
```

### Built-In Models

#### Uniform Acceleration

```ini
<external_gravity>
model = uniform
g1    = 0.0
g2    = 0.0
g3    = -1.0
```

Potential:

```math
\Phi = -(g_1 x + g_2 y + g_3 z).
```

Acceleration:

```math
\mathbf{a} = (g_1, g_2, g_3).
```

#### Cartesian Harmonic Potential

```ini
<external_gravity>
model  = cartesian_harmonic
omega1 = 1.0
omega2 = 0.0
omega3 = 0.0
```

Potential:

```math
\Phi = \frac{1}{2}\left(\omega_1^2 x^2 + \omega_2^2 y^2 + \omega_3^2 z^2\right).
```

Acceleration:

```math
\mathbf{a} = -(\omega_1^2 x,\omega_2^2 y,\omega_3^2 z).
```

This is useful for local box tests where the gravitational field is linear in one or more
directions.

#### Spherical Harmonic Potential

```ini
<external_gravity>
model = spherical_harmonic
omega = 1.0
```

Potential:

```math
\Phi = \frac{1}{2}\omega^2 r^2.
```

Acceleration:

```math
\mathbf{a} = -\omega^2 \mathbf{r}.
```

#### Point Mass / Plummer Potential

```ini
<external_gravity>
model     = point_mass
G         = 1.0
mass      = 1.0
softening = 0.01
```

Potential:

```math
\Phi = -\frac{GM}{\sqrt{r^2 + \epsilon^2}}.
```

Use `model=plummer` for the same functional form when the softening is physically a
Plummer scale rather than a numerical regularization.

#### NFW Halo

```ini
<external_gravity>
model     = nfw
G         = 1.0
rho_scale = 1.0
r_scale   = 1.0
```

Potential:

```math
\Phi = -4\pi G\rho_s r_s^2\,\frac{\ln(1+r/r_s)}{r/r_s}.
```

The implementation uses a regularized small-`r` series for
`\ln(1+x)/x`.

#### Logarithmic Halo

```ini
<external_gravity>
model       = logarithmic_halo
v0          = 1.0
core_radius = 1.0
q           = 1.0
```

Potential:

```math
\Phi = \frac{1}{2}v_0^2 \ln\left(r_c^2 + R^2 + \frac{z^2}{q^2}\right).
```

This is a cylindrical potential.  `q < 1` flattens the halo in the vertical direction.

#### Miyamoto-Nagai Disk

```ini
<external_gravity>
model     = miyamoto_nagai
G         = 1.0
disk_mass = 1.0
disk_a    = 1.0
disk_b    = 0.1
```

Potential:

```math
\Phi = -\frac{GM_{\rm d}}
{\sqrt{R^2 + \left(a + \sqrt{z^2 + b^2}\right)^2}}.
```

This is useful for disk-galaxy examples and for the gotham-style CGM setup.

#### Outer CGM Background

```ini
<external_gravity>
model    = outer_cgm
G        = 1.0
rho_mean = 1.0e-4
r_outer  = 5.0
```

Potential:

```math
\Phi = 4\pi G\rho_{\rm mean}
\left[\frac{4}{3}r_{\rm outer}^{3/2}r^{1/2} + \frac{1}{6}r^2\right].
```

This follows the outer-background form used in the gotham problem generators.

#### Gotham Composite

```ini
<external_gravity>
model     = gotham
G         = 1.0
rho_scale = 1.0
r_scale   = 1.0
disk_mass = 1.0
disk_a    = 1.0
disk_b    = 0.1
rho_mean  = 1.0e-4
r_outer   = 5.0
```

Potential:

```math
\Phi = \Phi_{\rm NFW} + \Phi_{\rm MN} + \Phi_{\rm outer}.
```

This should replace the pgen-local duplicated galaxy potentials used in the gotham
branches.

#### User-Defined Potential

```ini
<external_gravity>
model = user
```

Then edit the hook in `src/pgen/pgen.cpp`:

```cpp
namespace external_gravity {

KOKKOS_FUNCTION
Real UserExternalGravityPotential(Real x1, Real x2, Real x3) {
  return 0.5*(x1*x1 + x2*x2 + x3*x3);
}

} // namespace external_gravity
```

The function must be device-callable.  It should avoid host-only state, STL containers,
and non-Kokkos-compatible calls.

## Boundary Conditions: Required Next Feature

External gravity is not complete without boundary conditions that understand the
gravitational field.

The problem is simple: if the gravitational acceleration remains nonzero at a physical
boundary, a standard `outflow` or zero-gradient boundary copies interior pressure and
density into ghost zones without balancing gravity.  A nominally hydrostatic atmosphere
then develops spurious boundary forces, inflows, shocks, or mass loss.

A robust implementation needs one of two strategies:

1. Make the gravitational acceleration go smoothly to zero near the boundary.
2. Extrapolate thermodynamic state variables into ghost zones so hydrostatic equilibrium
   can extend through the boundary.

Both should be supported because they solve different problems.

### Strategy A: Taper the Gravitational Field Near Boundaries

This is the simplest option for simulations where the boundary should behave like an
ordinary outflow region.

Add optional taper parameters to `<external_gravity>`:

```ini
<external_gravity>
model = nfw

taper_gravity       = true
taper_width_x1_inner = 0.1
taper_width_x1_outer = 0.1
taper_width_x2_inner = 0.0
taper_width_x2_outer = 0.0
taper_width_x3_inner = 0.0
taper_width_x3_outer = 0.0
taper_function      = smoothstep
```

The acceleration should be multiplied by a smooth mask `w(x)` that approaches zero at
the selected boundary and one in the active domain interior.  For example, with
`s = d / width`,

```math
w(s) = 3s^2 - 2s^3,\qquad 0 \le s \le 1.
```

The potential-energy history should still use the physical potential unless the user
explicitly requests a tapered effective potential.  The docs must state which convention
is used.

Pros:

- easy to implement;
- works with existing boundary types;
- avoids ghost-zone hydrostatic extrapolation complexity.

Cons:

- changes the physical force near the boundary;
- not appropriate when the boundary is meant to represent a continuation of a
  stratified atmosphere.

### Strategy B: Hydrostatic Gravity Boundary

This is necessary for the proposed hydrostatic-atmosphere examples and for any run where
the boundary should continue a stratified equilibrium.

Add a new boundary mode, either as a new mesh boundary flag or as a user boundary helper:

```ini
<mesh>
ix1_bc = hydrostatic_gravity
ox1_bc = hydrostatic_gravity
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = hydrostatic_gravity
ox3_bc = hydrostatic_gravity

<external_gravity_boundary>
mode              = hydrostatic
thermo_closure    = isothermal
temperature_floor = 1.0e-12
density_floor     = 1.0e-12
pressure_floor    = 1.0e-12
velocity_mode     = no_inflow
```

If adding a new top-level block is too much, these keys can live in
`<external_gravity>` with a `boundary_` prefix.  A separate
`<external_gravity_boundary>` block is cleaner and easier to document.

#### Isothermal Closure

For an isothermal atmosphere with sound speed `c_s`, hydrostatic equilibrium gives

```math
\frac{dP}{dn} = -\rho \frac{d\Phi}{dn}, \qquad P = \rho c_s^2.
```

Between an active reference cell `a` and a ghost cell `g`,

```math
\rho_g = \rho_a \exp\left[-\frac{\Phi_g-\Phi_a}{c_s^2}\right],
```

```math
P_g = c_s^2 \rho_g.
```

This is the most important first implementation because it is stable, easy to verify,
and covers many atmosphere tests.

#### Polytropic / Adiabatic Closure

For a barotropic polytrope,

```math
P = K\rho^\gamma,
```

hydrostatic equilibrium can be written as

```math
h + \Phi = {\rm constant},
```

where

```math
h = \frac{\gamma}{\gamma-1}\frac{P}{\rho}.
```

Using the nearest active cell as the reference,

```math
h_g = h_a + \Phi_a - \Phi_g.
```

Then

```math
\rho_g =
\left[\frac{\gamma-1}{\gamma K}h_g\right]^{1/(\gamma-1)}.
```

If `h_g <= 0`, the boundary must apply the configured density/pressure floors and should
optionally warn once.

#### General Discrete Closure

For maximum robustness, also provide a discrete integration mode:

```ini
thermo_closure = integrate_pressure
```

For each ghost layer, integrate normal to the boundary:

```math
P_g = P_a - \frac{1}{2}(\rho_a + \rho_g)(\Phi_g-\Phi_a),
```

with an EOS-specific update for `rho_g`.  This can be solved by a few Newton or fixed
point iterations.  This is more flexible but should be implemented after the isothermal
and polytropic closures.

### Boundary Velocity Treatment

Hydrostatic thermodynamic extrapolation is not enough.  The boundary also has to choose
velocities consistently.

Recommended modes:

| `velocity_mode` | Behavior |
|-----------------|----------|
| `copy` | Copy all velocity components from the nearest active cell. |
| `zero_normal` | Set the normal component to zero and copy tangential components. |
| `no_inflow` | Copy velocities, but zero the normal component if it points into the domain. |
| `reflect_normal` | Reflect the normal velocity and copy tangential components. |

For hydrostatic-atmosphere examples, use `zero_normal` or `no_inflow`.

### MHD Boundary Treatment

For MHD, the hydrostatic boundary should update density and gas pressure, but magnetic
fields must still satisfy the existing face-centered boundary machinery.

Recommended first implementation:

- use existing magnetic-field boundary behavior (`outflow`, `reflect`, or `user`) for
  face-centered fields;
- apply hydrostatic extrapolation only to cell-centered fluid primitives;
- recompute conserved variables after primitive ghost zones are filled;
- document that magnetostatic equilibria require a dedicated MHD boundary if magnetic
  pressure/tension contributes to the force balance.

### Implementation Location

The boundary implementation should live near the existing physical boundary routines:

```text
src/bvals/physics/hydro_bcs.cpp
src/bvals/physics/bfield_bcs.cpp
src/bvals/bvals.hpp
```

External-gravity-specific helper functions can live in:

```text
src/srcterms/external_gravity.hpp
src/srcterms/external_gravity.cpp
```

The boundary code should call a shared device-callable potential/acceleration API so it
uses the exact same gravitational field as the source term and history output.

### Boundary Tests

Boundary support should not be accepted without tests:

1. Isothermal atmosphere in uniform gravity, 1-D, with `ix1_bc=hydrostatic_gravity` and
   `ox1_bc=hydrostatic_gravity`; verify density and pressure stay stationary.
2. Same atmosphere with ordinary `outflow`; verify the test detects larger boundary
   drift, proving the hydrostatic BC is doing real work.
3. Spherical atmosphere in point-mass or Plummer potential; verify radial profiles remain
   stable.
4. AMR version of the isothermal atmosphere; verify prolongation plus physical boundary
   filling does not create edge artifacts.
5. MPI version with domain decomposition; verify physical and internal boundaries do not
   diverge.
6. MHD smoke test with hydrostatic thermodynamics and simple copied magnetic fields.

## Efficiency Roadmap

The current implementation finite-differences `Phi` at six face points per cell per
source-term call.  That is simple and general, but it is not the efficient path for
production runs.

### Target Design

Use two APIs:

```cpp
KOKKOS_INLINE_FUNCTION
Real Potential(const ExternalGravityData &grav, Real x1, Real x2, Real x3);

KOKKOS_INLINE_FUNCTION
void Acceleration(const ExternalGravityData &grav, Real x1, Real x2, Real x3,
                  Real &a1, Real &a2, Real &a3);
```

For built-in models, `Acceleration()` should use analytic derivatives.  For `model=user`,
the default can finite-difference `UserExternalGravityPotential()`, with an optional
future user acceleration hook.

### Optional Caching

For static potentials, add cached arrays:

```text
phi(m,k,j,i)
grav_accel(m,3,k,j,i)
```

These can be allocated on the pack or source-term object and rebuilt when:

- the mesh pack is created;
- AMR changes the pack layout;
- load balancing changes local MeshBlocks;
- user explicitly requests recomputation.

Benefits:

- source updates become one array read instead of repeated potential evaluations;
- history output uses the same cached `phi`;
- hydrostatic boundaries can use the same `phi` values in active and ghost zones.

Costs:

- more memory;
- must handle AMR and restarts carefully;
- requires ghost-zone-compatible coordinate evaluation.

The first optimized pass should add analytic accelerations.  Caching is a second pass
unless profiling shows potential evaluation dominates.

## Example Pgens and Inputs

The final feature should include both simple input-only examples and useful pgens.

### Input-Only Examples

Place user-facing examples under:

```text
inputs/hydro/external_gravity/
inputs/mhd/external_gravity/
```

Recommended examples:

| File | Purpose |
|------|---------|
| `uniform_gravity.athinput` | Constant acceleration in a box; demonstrates source-term behavior and `grav-PE`. |
| `harmonic_box.athinput` | Cartesian or spherical harmonic potential. |
| `point_mass.athinput` | Plummer-softened point mass. |
| `nfw_halo.athinput` | Static spherical halo. |
| `miyamoto_nagai_disk.athinput` | Cylindrical disk potential. |
| `gotham_composite.athinput` | NFW plus Miyamoto-Nagai plus outer CGM background. |
| `hydrostatic_isothermal.athinput` | Requires the hydrostatic gravity boundary condition. |

### Useful Pgens

Recommended pgens:

#### `external_gravity_hydrostatic.cpp`

Initializes an isothermal atmosphere satisfying

```math
\rho = \rho_0 \exp\left[-\frac{\Phi-\Phi_0}{c_s^2}\right].
```

This should be the primary validation pgen for gravity-aware boundaries.

Required inputs:

```ini
<problem>
pgen_name = external_gravity_hydrostatic
rho0      = 1.0
cs        = 1.0
x1_ref    = 0.0
x2_ref    = 0.0
x3_ref    = 0.0
```

#### `external_gravity_orbit.cpp`

Initializes a cold ring or cloud in a spherical potential with the circular velocity
computed from the same acceleration API used by the source term.

This validates signs, normalization, and cylindrical/spherical coordinate conventions.

#### `external_gravity_disk.cpp`

Initializes a rotating disk in a Miyamoto-Nagai or logarithmic-halo potential.  This is a
useful example for users who want realistic galaxy simulations.

#### `external_gravity_gotham_cgm.cpp`

A stripped-down CGM atmosphere using `model=gotham`.  The goal is to replace duplicated
potential code in gotham pgens with the generic source-term infrastructure.

#### `external_gravity_custom.cpp`

Demonstrates `model=user` by defining a compact custom potential in `pgen.cpp` or a
small dedicated pgen file.  This should be deliberately minimal so users can copy it.

## Testing Plan

### Current Test

The current branch includes:

```text
tst/test_suite/nr/test_nr_external_gravity_cpu.py
tst/inputs/external_gravity.athinput
```

It checks:

- `<external_gravity>` alone activates source terms;
- hydro history contains `grav-PE`;
- `grav-PE` matches a discrete Cartesian harmonic potential sum;
- uniform acceleration changes total momentum by `g * t`.

### Required Robustness Tests

Add tests in this order:

1. **Model parser tests**: invalid model names, negative masses, zero scale radii, bad
   logarithmic-halo core radius, and disabled block behavior.
2. **Analytic acceleration tests**: each built-in model compared against analytic
   acceleration at representative points.
3. **History tests**: hydro and MHD, ideal and isothermal EOS, with and without passive
   scalars.
4. **Boundary tests**: isothermal hydrostatic atmosphere with gravity-aware boundaries.
5. **AMR tests**: hydrostatic and force-update tests with refinement enabled.
6. **MPI tests**: at least one uniform-gravity and one hydrostatic-boundary test under
   domain decomposition.
7. **Restart tests**: run half, restart, and compare against a continuous run.
8. **User hook test**: compile a simple custom potential and verify source/history
   behavior.

### Validation Commands

Baseline local checks:

```bash
cd tst
python run_test_suite.py --test test_suite/nr/test_nr_external_gravity_cpu.py --cpu
python run_test_suite.py --style
```

After adding AMR/MPI coverage:

```bash
cd tst
python run_test_suite.py --test test_suite/nr/test_nr_external_gravity_hydrostatic_cpu.py --cpu
python run_test_suite.py --test test_suite/nr/test_nr_external_gravity_mpicpu.py --mpicpu
```

## Robustness Checklist

Before considering the feature review-ready:

- `<external_gravity>` works without legacy source-term blocks.
- Every built-in model has documented equations and input parameters.
- Built-in models use analytic acceleration where possible.
- `model=user` is documented with a minimal working example.
- The code rejects unsupported relativistic or coordinate-system combinations clearly.
- `grav-PE` works for hydro and MHD.
- `grav-PE` indexing is safe with passive scalars and does not exceed
  `NREDUCTION_VARIABLES`.
- Hydrostatic gravity boundaries exist and are tested.
- Gravity tapering exists or is explicitly deferred with documented limitations.
- Restart behavior is tested.
- AMR behavior is tested.
- MPI behavior is tested.
- Example pgens are present and documented.
- The Pages docs include this page and link it from the source-terms module page.

## Recommended Implementation Order

1. Refactor `external_gravity` into model-specific parsing and add analytic acceleration.
2. Add the hydrostatic gravity boundary condition for isothermal atmospheres.
3. Add `external_gravity_hydrostatic.cpp` and the corresponding input/test.
4. Add AMR and MPI tests for the hydrostatic setup.
5. Add remaining example inputs for the built-in potentials.
6. Add orbit/disk/gotham/custom example pgens.
7. Integrate this page into the GitHub Pages documentation branch.
8. Run the full targeted validation suite and then open the PR.

## Notes for Future Extensions

- Time-dependent external potentials should be a separate feature.  The present design
  assumes a static potential.
- Self-gravity should remain separate.  This module is for externally prescribed
  potentials only.
- Particle gravity can reuse the same potential/acceleration API, but particle orbit
  integration should be handled in `src/particles`.
- Magnetostatic boundaries require additional terms beyond hydrostatic pressure balance
  and should not be hidden behind the first hydrostatic boundary implementation.
