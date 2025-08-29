# div(B) Conservation Fix Comparison

## Summary
The `fix-amr-divb-issue` branch adds a critical fix for div(B) conservation during AMR that is **NOT present** in the current branch. This explains why div(B) is not conserved when using AMR in the current branch.

## The Problem
When AMR refines a coarse meshblock adjacent to an already-fine meshblock, the prolongation operator creates face-centered magnetic field values at the shared boundary that may not match the existing fine values. This mismatch violates the divergence-free constraint (div(B) = 0) for magnetic fields.

## The Fix (in `fix-amr-divb-issue` branch)

### 1. New Function Added
The fix-amr-divb-issue branch adds a new function `FixRefinedFCBoundaries()` in `mesh_refinement.cpp`:

**Location**: src/mesh/mesh_refinement.cpp:1195-1385
**Declaration**: src/mesh/mesh_refinement.hpp:116

### 2. Function Purpose
`FixRefinedFCBoundaries()` corrects face-centered B-field values after refinement by:
- Identifying blocks that were just refined (refine_flag > 0)
- Finding their neighbors that were already fine (refine_flag == 0) 
- Copying face-centered B-field values from fine neighbors to newly refined blocks
- This overwrites the incorrectly prolongated values with correct values from fine neighbors

### 3. Where It's Called
In `RedistAndRefineMeshBlocks()` at line 639, right after `RefineFC()`:
```cpp
if (pmhd != nullptr) {
  RefineCC(new_to_old, pmhd->u0, pmhd->coarse_u0);
  RefineFC(new_to_old, pmhd->b0, pmhd->coarse_b0);
  // IMPORTANT: Fix div(B) errors at fine/coarse boundaries after refinement
  FixRefinedFCBoundaries(new_to_old, pmhd->b0);  // <-- THIS IS THE FIX
}
```

## Key Implementation Details

### Algorithm
1. Loop over all MeshBlocks on this rank that were just refined
2. For each refined block, check all other blocks on same rank at same level
3. If a neighbor was already fine (not refined), copy its boundary face values
4. Handle all 6 faces (±x1, ±x2, ±x3) as appropriate for dimensionality

### Limitations
- Only handles same-rank neighbors (cross-rank boundaries handled by normal boundary communication)
- Includes diagnostic output showing number of boundaries fixed

## Comparison with Current Branch

| Aspect | Current Branch | fix-amr-divb-issue Branch |
|--------|---------------|---------------------------|
| Has FixRefinedFCBoundaries() | ❌ No | ✅ Yes |
| Maintains div(B)=0 after AMR | ❌ No | ✅ Yes |
| Issue #595 addressed | ❌ No | ✅ Yes |
| Diagnostic output | ❌ No | ✅ Yes (boundary fix count) |

## Related Commits
- `93b4071e`: "Fix div(B) issue in AMR with magnetic fields" - Main implementation
- `0c907093`: "Complete div(B) test suite and documentation with style fixes"
- `98ff81d7`: "Add comprehensive MHD+AMR test suite for div(B) verification"

## Integration Path
To fix the div(B) conservation issue in the current branch:

1. **Cherry-pick the fix commit**:
   ```bash
   git cherry-pick 93b4071e
   ```

2. **Or manually integrate**:
   - Add `FixRefinedFCBoundaries()` declaration to mesh_refinement.hpp
   - Add the function implementation to mesh_refinement.cpp
   - Call it after `RefineFC()` in `RedistAndRefineMeshBlocks()`

## Why This Matters for Turbulence + AMR
The catastrophic instability at t≈1.035 could be related to div(B) violations accumulating over time:
- div(B) errors act as spurious magnetic monopoles
- These create unphysical forces and energy injection
- Combined with turbulence forcing, this could lead to exponential growth
- The timing (after multiple turbulence correlation times) suggests cumulative effect

## Recommendation
**Integrate the div(B) fix from the `fix-amr-divb-issue` branch before addressing volume scaling issues**. The div(B) conservation violation is likely the root cause of the catastrophic instability, not the volume scaling.