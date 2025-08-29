# AMR + Turbulence Volume Scaling Fix Plan

## Problem Summary
The simulation experiences catastrophic failure around t≈1.035 when running with both AMR and turbulence driving:
- dt plummets from ~1e-3 to ~1e-19
- Total energy explodes from ~1e5 to ~1e22  
- Momentum components blow up (e.g., x-momentum reaches -8.35e6)

### Timing Analysis
- AMR triggers at t=0.5, refining from 8 to 64 blocks
- Simulation runs normally from t=0.5 to t=0.92
- Gradual dt decrease starts at t≈0.92
- Catastrophic failure at t≈1.035
- Notable: tcorr=1.0 in input file, suggesting issue may relate to turbulence mode update

## Root Cause Hypothesis
The turbulence forcing normalization in `UpdateForcing()` accounts for varying cell volumes when calculating energy injection rates, but the force application in `AddForcing()` doesn't properly scale for different refinement levels. This creates an imbalance where refined cells receive disproportionate forcing.

## Diagnostic Plan (Do First)

### 1. Track Force Magnitudes by Refinement Level
Add diagnostics in `turb_driver.cpp` to monitor:
```cpp
// In UpdateForcing() - after force normalization
if (global_variable::my_rank == 0) {
  // Track force magnitude on base level blocks vs refined blocks
  // Log: block level, force magnitude, cell volume
}
```

### 2. Monitor Energy Injection Rates
Track actual vs target energy injection:
- Per refinement level
- Before and after normalization
- At mode update times (when InitializeModes is called)

### 3. Log Turbulence Mode Updates
Add logging to identify when new random modes are generated:
```cpp
// In InitializeModes()
if (global_variable::my_rank == 0) {
  std::cout << "### InitializeModes called at t=" << pm->time 
            << ", updating random amplitudes" << std::endl;
}
```

## Fix Implementation Options

### Option A: Force Array Contains Acceleration (force per unit mass)
If this is the case:
- Force should NOT be scaled in `AddForcing()`
- But normalization in `UpdateForcing()` needs to account for cell volumes properly

### Option B: Force Array Contains Force Density (force per unit volume)  
If this is the case:
- Scale force by `cell_volume/base_volume` in `AddForcing()`
- Keep current normalization in `UpdateForcing()`

### Determining Which Option:
Examine the units and usage:
- In `AddForcing()`: `u0(m,IM1,k,j,i) += den*a1*bdt`
- This suggests `a1` is acceleration (force/mass), supporting Option A
- But normalization uses volumes, suggesting Option B intent

## Proposed Fix

### Step 1: Add Volume Scaling (if Option B is correct)
```cpp
// In AddForcing() - before applying force
Real Lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
Real Ly = pm->mesh_size.x2max - pm->mesh_size.x2min;  
Real Lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
// Calculate base cell volume
Real base_vol = (Lx*Ly*Lz) / (mesh_nx1 * mesh_nx2 * mesh_nx3);

// Get current cell volume
Real cell_vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

// Scale force to maintain constant force per unit volume
Real vol_scale = cell_vol / base_vol;
Real a1 = force_(m,0,k,j,i) * vol_scale;
Real a2 = force_(m,1,k,j,i) * vol_scale;
Real a3 = force_(m,2,k,j,i) * vol_scale;
```

### Step 2: Add Safety Limiter
```cpp
// Prevent catastrophic acceleration
Real a_mag = sqrt(a1*a1 + a2*a2 + a3*a3);
Real cs = 1.0;  // Approximate sound speed
Real a_max = cs / bdt;  // Maximum safe acceleration
if (a_mag > a_max) {
  Real scale = a_max / a_mag;
  a1 *= scale;
  a2 *= scale;
  a3 *= scale;
  if (global_variable::my_rank == 0) {
    std::cout << "WARNING: Limiting acceleration from " << a_mag 
              << " to " << a_max << std::endl;
  }
}
```

### Step 3: Monitor dt Evolution
Add check for rapid dt decrease as early warning:
```cpp
static Real last_dt = 0.0;
if (last_dt > 0.0 && dt < 0.1 * last_dt) {
  std::cout << "WARNING: dt dropped by >10x from " << last_dt 
            << " to " << dt << std::endl;
}
last_dt = dt;
```

## Testing Strategy

### 1. Short Diagnostic Run
- Run until t=1.1 with full diagnostics
- Compare force values on level 0 vs level 1 blocks
- Check if energy injection rate is consistent

### 2. Verify Fix
- Implement volume scaling based on diagnostic results
- Run past t=1.035 to verify no crash
- Check that energy remains bounded

### 3. Production Test
- Remove diagnostics
- Run full simulation with fix
- Verify physical results are reasonable

## Files to Modify

1. **src/srcterms/turb_driver.cpp**
   - `UpdateForcing()`: Add diagnostics, fix normalization if needed
   - `AddForcing()`: Add volume scaling and safety limiter
   - `InitializeModes()`: Add logging

2. **src/srcterms/turb_driver.hpp**
   - Add any new member variables for tracking

## Notes for Implementation

- The key is determining whether the force array represents acceleration or force density
- The timing near tcorr=1.0 suggests the issue may be triggered by mode updates
- Safety limiters are essential to prevent future instabilities
- Consider making volume scaling a runtime parameter for testing

## Current Status
- AMR + turbulence integration is working (arrays resize correctly)
- Volume scaling code was temporarily added but needs proper testing
- This plan should be followed when ready to fix the physics issue