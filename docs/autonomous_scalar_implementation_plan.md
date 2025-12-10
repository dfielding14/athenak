# Autonomous Scalar Evolution Implementation Plan

## Overview

This document outlines the implementation of a simulation mode where:
- Velocity and density fields are **frozen** (time-independent)
- Passive scalars are **evolved** via advection (and eventually diffusion)
- Initial velocity field has a specified **power-law spectrum**

This enables studying turbulent mixing/diffusion in a statistically stationary flow.

---

## Design Principles

**These principles govern all implementation decisions:**

1. **Simple, minimally invasive code**
   - Prefer small, targeted changes over sweeping modifications
   - Reuse existing infrastructure wherever possible
   - Avoid introducing new abstractions unless absolutely necessary
   - If a feature can be achieved by adding 5 lines instead of 50, use 5 lines

2. **Match the existing codebase style**
   - Follow the formatting, naming conventions, and patterns already in use
   - Study similar existing code before implementing new code
   - When in doubt, copy the approach used elsewhere in AthenaK

3. **Clearly commented and documented**
   - Every new function gets a descriptive header comment (following existing style)
   - Non-obvious logic gets inline comments explaining the "why"
   - Input parameters documented in the implementation plan and code

4. **NEVER GUESS**
   - If uncertain about the correct approach, **STOP and ask for clarification**
   - If unsure how existing code works, **read it thoroughly first**
   - If a design decision has multiple valid options, **ask the user to choose**
   - Incorrect assumptions waste more time than asking questions

---

## Physics Background

### Power Spectrum Convention

We use the standard definition where the 3D energy spectrum E(k) relates to the structure function:

For E(k) ~ k^{-α}:
- Second-order structure function: S₂(ℓ) = ⟨|δv(ℓ)|²⟩ ~ ℓ^{α-1}  (for 1 < α < 3)
- Velocity increment: δv(ℓ) ~ ℓ^{(α-1)/2}

**Kolmogorov turbulence:**
- E(k) ~ k^{-5/3}  →  α = 5/3
- δv(ℓ) ~ ℓ^{1/3}  ✓

### Fourier Amplitude Scaling

For 3D isotropic turbulence, E(k) = 4πk² |û(k)|²

If Fourier amplitude |û(k)| ~ k^{-β}, then:
- E(k) ~ k^{2-2β}

For E(k) ~ k^{-5/3}:
- 2 - 2β = -5/3  →  β = 11/6

In the existing turb_driver code (`spect_form=2`):
```cpp
norm = 1.0/pow(kiso, (expo+2.0)/2.0);
```
So β = (expo+2)/2, meaning **expo = 5/3** gives Kolmogorov spectrum.

---

## Implementation Phases

### Phase 1: Add `scalar_only` Flag to Hydro ✓ COMPLETE

**Goal:** Allow freezing hydro variables while evolving scalars.

- [x] **1.1** Add member variable to `Hydro` class
  - File: `src/hydro/hydro.hpp`
  - Add: `bool scalar_only;`

- [x] **1.2** Read parameter from input file
  - File: `src/hydro/hydro.cpp` (constructor)
  - Add: `scalar_only = pin->GetOrAddBoolean("hydro", "scalar_only", false);`

- [x] **1.3** Modify RK update loop to skip hydro variables
  - File: `src/hydro/hydro_update.cpp`
  - Location: Line ~50, the `par_for_outer` loop
  - Change:
    ```cpp
    // OLD: par_for_outer(..., 0, nvar-1, ...)
    // NEW:
    int nstart = scalar_only ? nhydro : 0;
    par_for_outer("h_update", DevExeSpace(), scr_size, scr_level,
                  0, nmb1, nstart, nvar-1, ks, ke, js, je, ...)
    ```

- [x] **1.4** Verify `scalar_only` is accessible in RKUpdate
  - Member variable accessible directly, no lambda capture needed

**Test 1.1:** ✓ Compile and run with `scalar_only = false` (default behavior unchanged)

**Test 1.2:** ✓ Run with `scalar_only = true`:
  - Initialize uniform density, uniform velocity, scalar blob
  - Verify: density and velocity remain EXACTLY constant (max_change=0.00e+00)
  - Verify: scalar advects with the flow

---

### Phase 2: Timestep Calculation for Scalar-Only Mode ✓ COMPLETE

**Goal:** Use velocity-only CFL when hydro is frozen (no sound speed contribution).

- [x] **2.1** Modify timestep calculation
  - File: `src/hydro/hydro_newdt.cpp`
  - Location: Line ~54 where `time_evolution == TimeEvolution::kinematic` is checked
  - Change logic to also use kinematic branch when `scalar_only == true`:
    ```cpp
    if (pdrive->time_evolution == TimeEvolution::kinematic || scalar_only) {
      // velocity-only CFL
    }
    ```

- [x] **2.2** Ensure `scalar_only` flag is accessible in `NewTimeStep`
  - It's a member of `Hydro`, so `this->scalar_only` works directly

**Test 2.1:** ✓ Compare timestep with `scalar_only=true` vs `evolution=kinematic`
  - Both produce dt=1.250e-2 (identical)

**Test 2.2:** ✓ Verify timestep is NOT affected by sound speed when `scalar_only=true`
  - scalar_only=false: dt=5.456e-3 (includes cs)
  - scalar_only=true:  dt=1.250e-2 (velocity only, ~2.3x larger)

---

### Phase 3: Velocity Field Initialization with Specified Spectrum ✓ COMPLETE

**Goal:** Create function to generate a turbulent velocity field with:
- Power-law spectrum E(k) ~ k^{-α}
- Specified kmin, kmax
- Specified mean velocity (default 0)
- Specified velocity dispersion (v_rms)
- Solenoidal (divergence-free) option

**Implementation:** Added `InitTurbulentVelocity()` function directly to `turbulent_box.cpp` (simpler than creating separate module).

- [x] **3.1** Implement `InitTurbulentVelocity()` function
  - File: `src/pgen/turbulent_box.cpp`
  - Fourier-based velocity field generation with power-law spectrum

- [x] **3.2** Parameters (read from `<problem>` block):
  ```
  turb_vel_init  = true/false  # Enable turbulent velocity
  turb_v_rms     = 1.0         # Target velocity dispersion
  turb_nlow      = 1           # Min wavenumber index
  turb_nhigh     = 4           # Max wavenumber index
  turb_expo      = 1.6667      # Power-law exponent (5/3 for Kolmogorov)
  turb_sol_frac  = 1.0         # Solenoidal fraction (1.0 = div-free)
  turb_rseed     = 12345       # Random seed
  ```

- [x] **3.3** Core algorithm implemented:
  1. Count and generate Fourier modes in range [nlow, nhigh]
  2. Set amplitude |û(k)| ~ k^{-(expo+2)/2} so E(k) ~ k^{-expo}
  3. Apply solenoidal projection using sol_frac parameter
  4. Sum modes to real-space velocity using Kokkos parallel_for
  5. Compute mean momentum with MPI reduction
  6. Subtract mean and scale to target v_rms

- [x] **3.4** MPI-safe momentum subtraction
  - Kokkos::parallel_reduce for local sums
  - MPI_Allreduce for global sums

- [x] **3.5** v_rms scaling implemented
  - Compute current_vrms = sqrt(⟨v²⟩ - ⟨v⟩²)
  - Scale velocity by (target_v_rms / current_vrms)

**Test 3.1:** ✓ Generate velocity field, verify:
  - Mean velocity ≈ 0 (10^-16 machine precision)
  - v_rms = 1.0000 (matches target exactly)

**Test 3.2:** (Future) Compute power spectrum of generated field

**Test 3.3:** (Future) Verify divergence-free when sol_frac=1.0

---

### Phase 4: Integrate with Problem Generator ✓ COMPLETE

**Goal:** Modify turbulent_box.cpp to use new velocity initialization.

**Note:** Phases 3 and 4 were combined - initialization is already in turbulent_box.cpp.

- [x] **4.1** Add input parameters to turbulent_box.cpp
  - Read turb_vel_init boolean to switch modes
  - If false: use uniform velocity (vx0, vy0, vz0)
  - If true: call InitTurbulentVelocity()
  ```
  <problem>
  turb_vel_init  = true
  turb_v_rms     = 1.0
  turb_nlow      = 1
  turb_nhigh     = 4
  turb_expo      = 1.6667    # 5/3
  sol_fraction = 1.0
  turb_rseed  = 12345
  ```

- [x] **4.2** Call velocity initialization in ProblemGenerator
  - After density/energy initialized
  - Before scalar initialization
  - Energy updated with kinetic energy after velocity set

- [x] **4.3** Initialize scalars as in current turbulent_box.cpp
  - Kept existing scalar initialization (random dye concentration)
  - Refactored into separate parallel_for for clarity

- [x] **4.4** Turbulence driving not needed for scalar_only mode
  - No turb_driving block needed in input file
  - Velocity initialized once and frozen

**Test 4.1:** ✓ Full simulation test (32^3, 3D)
  - scalar_only=true, turb_vel_init=true
  - Verified: max(|Δρ|) = 0, max(|Δv|) = 0
  - Velocity field frozen, scalars evolve

**Test 4.2:** (Future) Conservation test

**Test 4.3:** (Future) Qualitative visualization test

---

### Phase 5: Add Passive Scalar Diffusion

**Goal:** Enable scalar diffusion in addition to advection.

- [ ] **5.1** Check existing diffusion infrastructure
  - File: `src/diffusion/` directory
  - Understand how viscosity/resistivity are implemented

- [ ] **5.2** Add scalar diffusion coefficient
  - File: `src/hydro/hydro.hpp`
  - Add: `Real scalar_diffusivity;` or array for multiple scalars

- [ ] **5.3** Read from input
  - `<hydro>` block: `scalar_diffusivity = 0.01`
  - Or per-scalar: `scalar_diff_0 = 0.01`, etc.

- [ ] **5.4** Implement diffusion operator
  - Options:
    a. Explicit: add ∇²s term in source terms
    b. Use existing diffusion framework
  - Compute: ∂s/∂t += D∇²s

- [ ] **5.5** Add diffusion timestep constraint
  - dt_diff ~ dx² / (2*D*ndim)
  - Take minimum with advection dt

- [ ] **5.6** Handle boundary conditions for diffusion
  - Periodic: automatic
  - Other BCs: may need ghost zone treatment

**Test 5.1:** Pure diffusion test (v=0)
  - Initialize Gaussian scalar blob
  - Verify it spreads as expected: σ² ~ σ₀² + 2Dt

**Test 5.2:** Advection-diffusion test
  - Uniform velocity + diffusion
  - Verify blob advects AND spreads

**Test 5.3:** Diffusion coefficient test
  - Run with different D values
  - Verify spreading rate scales with D

---

## File Summary

| File | Changes |
|------|---------|
| `src/hydro/hydro.hpp` | Add `scalar_only`, `scalar_diffusivity` |
| `src/hydro/hydro.cpp` | Read new parameters |
| `src/hydro/hydro_update.cpp` | Conditional loop start |
| `src/hydro/hydro_newdt.cpp` | Kinematic CFL for scalar_only |
| `src/pgen/turb_init.hpp` | New file: velocity init declaration |
| `src/pgen/turb_init.cpp` | New file: velocity init implementation |
| `src/pgen/turbulent_box.cpp` | Use new velocity init |
| `src/pgen/CMakeLists.txt` | Add new source files |

---

## Input File Template

```
<comment>
Autonomous scalar evolution in frozen turbulent velocity field

<job>
problem_id = turb_scalar

<time>
evolution    = kinematic
cfl_number   = 0.4
tlim         = 10.0
nlim         = -1

<mesh>
nghost = 2
nx1    = 128
x1min  = 0.0
x1max  = 1.0
ix1_bc = periodic
ox1_bc = periodic

nx2    = 128
x2min  = 0.0
x2max  = 1.0
ix2_bc = periodic
ox2_bc = periodic

nx3    = 128
x3min  = 0.0
x3max  = 1.0
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = 32
nx2 = 32
nx3 = 32

<hydro>
eos         = isothermal
iso_sound_speed = 1.0
nscalars    = 1
scalar_only = true
# scalar_diffusivity = 0.01  # Uncomment for Phase 5

<problem>
# Velocity field parameters
v_rms        = 1.0
nlow         = 1
nhigh        = 10
expo         = 1.6666667   # 5/3 for Kolmogorov
sol_fraction = 1.0
turb_rseed   = 12345

# Scalar initialization
scalar_init  = blob        # Options: blob, gradient, random
blob_radius  = 0.1
blob_center_x = 0.5
blob_center_y = 0.5
blob_center_z = 0.5

<output1>
file_type   = hdf5
variable    = hydro_w
dt          = 0.1
```

---

## Testing Checklist Summary

| Phase | Test | Description | Pass Criterion |
|-------|------|-------------|----------------|
| 1 | 1.1 | Default behavior unchanged | Existing tests pass |
| 1 | 1.2 | scalar_only freezes hydro | max(Δρ), max(Δv) < ε_mach |
| 2 | 2.1 | Kinematic CFL | dt matches kinematic mode |
| 3 | 3.1 | v_rms scaling | Measured v_rms = target |
| 3 | 3.2 | Power spectrum | E(k) ~ k^{-5/3} verified |
| 3 | 3.3 | Divergence-free | max(∇·v) < ε_mach |
| 4 | 4.1 | Full frozen evolution | Hydro frozen, scalars evolve |
| 4 | 4.2 | Scalar conservation | ∫ρs dV conserved |
| 5 | 5.1 | Pure diffusion | Gaussian spreading correct |
| 5 | 5.2 | Advection-diffusion | Combined behavior correct |

---

## Execution Order

1. Phase 1 (scalar_only flag) - foundation
2. Phase 2 (timestep) - required for correct CFL
3. Phase 3 (velocity init) - core new functionality
4. Phase 4 (pgen integration) - put it all together
5. Phase 5 (diffusion) - final enhancement

Estimated effort: Phases 1-2 are quick (~1-2 hours each). Phase 3 is the bulk of the work (~1 day). Phase 4 is integration (~half day). Phase 5 depends on existing infrastructure (~half to full day).
