# AMR + Turbulence + MHD Fixes Summary

## Critical Fixes Applied

### 1. div(B) Conservation Fix (mesh_refinement.cpp/hpp)
**Problem**: When AMR refines blocks adjacent to already-fine blocks, the prolongated face-centered B-field values don't match existing fine values, violating div(B)=0 and causing catastrophic instabilities.

**Solution**: Integrated `FixRefinedFCBoundaries()` function from `fix-amr-divb-issue` branch that:
- Copies correct B-field values from already-fine neighbors to newly-refined blocks
- Maintains div(B)=0 constraint after AMR refinement
- Includes diagnostic output showing boundary fixes

**Files Modified**:
- `src/mesh/mesh_refinement.hpp`: Added function declaration (line 122-123)
- `src/mesh/mesh_refinement.cpp`: Added function implementation (lines 1196-1391) and function call (line 645)

### 2. Derived Variable Reallocation Fix (derived_variables.cpp)
**Problem**: After AMR increases mesh blocks from 8 to 64, derived variable arrays weren't properly reallocated, causing memory access faults when writing outputs like `mhd_divb`.

**Solution**: Fixed reallocation condition to check if nmb (number of mesh blocks) has changed:
- Changed from: `if (derived_var.extent(4) <= 1)`
- Changed to: `if (derived_var.extent(0) != nmb || derived_var.extent(4) <= 1)`
- Applied to all derived variable allocations in the file

**File Modified**:
- `src/outputs/derived_variables.cpp`: Updated all reallocation conditions (10+ locations)

### 3. div(B) Loop Limits Fix (derived_variables.cpp) 
**Problem**: Ghost zone loop limits for div(B) calculation used `else if` incorrectly, preventing proper k-limits in 3D.

**Solution**: Changed from `else if (three_d)` to `if (three_d)` to properly set k-limits in 3D problems.

**File Modified**:
- `src/outputs/derived_variables.cpp`: Lines 969-976

## Testing Recommendations

1. **Compile the code** with the standard build process
2. **Test with your input file** including both:
   ```
   <output5>
   file_type   = bin
   variable    = mhd_divb
   dt          = 0.1
   slice_x3    = 0.0
   
   <output6>
   file_type   = bin  
   variable    = turb_force
   dt          = 0.1
   slice_x3    = 0.0
   ```

3. **Monitor for**:
   - No memory access faults after AMR triggers at t=0.5
   - Stable dt evolution (no catastrophic drops)
   - div(B) diagnostic output showing boundary fixes
   - Energy conservation (no exponential growth)

## Expected Improvements

1. **Stability**: Should eliminate catastrophic instability at t≈1.035
2. **div(B) Conservation**: Magnetic field divergence maintained at machine precision
3. **Output Reliability**: Both `mhd_divb` and `turb_force` outputs should work correctly after AMR
4. **Diagnostic Information**: Will see messages like "[FixDivB] Fixed N boundaries" when AMR occurs

## Verification Steps

To verify the fixes are working:
1. Check that the simulation runs past t=1.035 without crashes
2. Confirm output files are generated for both mhd_divb and turb_force
3. Look for "[FixDivB]" messages in output when refinement occurs at t=0.5
4. Monitor that dt remains around 1e-3 (not dropping to 1e-14 or lower)

## Future Considerations

If instabilities persist, revisit `AMR_TURBULENCE_FIX_PLAN.md` for volume scaling implementation (currently not applied, waiting to see if div(B) fix resolves the issue).