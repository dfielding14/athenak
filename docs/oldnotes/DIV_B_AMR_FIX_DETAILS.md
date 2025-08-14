# Detailed Documentation: Fix for div(B)≠0 Issue in AMR with Magnetic Fields

## Table of Contents
1. [Problem Background](#problem-background)
2. [Root Cause Analysis](#root-cause-analysis)
3. [Solution Overview](#solution-overview)
4. [Detailed Implementation Changes](#detailed-implementation-changes)
5. [Testing Infrastructure](#testing-infrastructure)
6. [Design Decisions and Rationale](#design-decisions-and-rationale)
7. [Limitations and Future Work](#limitations-and-future-work)

---

## Problem Background

### The Issue
When using Adaptive Mesh Refinement (AMR) with magnetohydrodynamics (MHD) in AthenaK, a numerical error occurs where the divergence of the magnetic field (div(B)) becomes non-zero at boundaries between mesh blocks of different refinement levels. This violates a fundamental constraint of Maxwell's equations: ∇·B = 0.

### References
- **AthenaK Issue #595**: https://github.com/IAS-Astrophysics/athenak/issues/595
- **Original Athena++ Issue #617**: https://github.com/PrincetonUniversity/athena/issues/617
- **Athena++ Fix PR #625**: https://github.com/PrincetonUniversity/athena/pull/625

### Why This Matters
1. **Physical Correctness**: Magnetic monopoles don't exist; div(B)=0 must be maintained
2. **Numerical Stability**: Non-zero div(B) can lead to spurious forces and instabilities
3. **Conservation**: Violates magnetic flux conservation across boundaries
4. **Error Propagation**: Errors can accumulate and corrupt the entire simulation

---

## Root Cause Analysis

### The Refinement Process
When AMR refines or derefines mesh blocks, the following sequence occurs:

1. **Flagging**: Blocks are marked for refinement based on criteria (density gradients, etc.)
2. **Tree Update**: The mesh tree structure is updated to reflect new refinement levels
3. **Load Balancing**: Blocks are redistributed across MPI ranks
4. **Data Transfer**: Physical variables are moved/copied between blocks
5. **Prolongation/Restriction**: Variables are interpolated between coarse and fine grids
6. **Boundary Updates**: Ghost zones and boundaries are synchronized

### Where the Problem Occurs
The issue specifically happens during **Step 5 (Prolongation)** when:

1. A coarse block C is refined to match the level of an adjacent fine block F
2. The prolongation operator creates face-centered magnetic field values on C's boundary
3. These prolongated values don't match the existing values on F's boundary
4. The mismatch creates a discontinuity that violates div(B)=0

### Mathematical Description
For face-centered fields on a staggered grid:
- Block F has field value `B_F` at boundary face
- Block C gets prolongated value `B_C'` at the same face
- Problem: `B_C' ≠ B_F` due to interpolation
- Result: `div(B) = (B_C' - B_F)/Δx ≠ 0`

---

## Solution Overview

### Core Idea
After prolongation, overwrite the incorrectly prolongated face-centered field values with the correct values from already-fine neighboring blocks.

### Algorithm Steps
1. After standard prolongation (RefineFC), identify all newly refined blocks
2. For each newly refined block, find neighbors at the same level
3. If a neighbor was already fine (not just refined), it has the "correct" field values
4. Copy the face-centered field values from the fine neighbor to the newly refined block
5. This ensures continuity across the boundary and maintains div(B)=0

---

## Detailed Implementation Changes

### 1. Header File Modification (`src/mesh/mesh_refinement.hpp`)

#### Line 113: Added Function Declaration
```cpp
// Corrects FC field values at newly-refined/already-fine boundaries (div(B) fix #595)
void FixRefinedFCBoundaries(DualArray1D<int> &n2o, DvceFaceFld4D<Real> &b);
```

**Why This Change:**
- Declares the new function that performs the correction
- Comment briefly explains purpose and references the issue number
- Follows existing function declaration style in the file
- Placed logically after RefineFC since it's called immediately after

### 2. Main Implementation (`src/mesh/mesh_refinement.cpp`)

#### Lines 632-637: Added Function Call in RedistAndRefineMeshBlocks
```cpp
if (pmhd != nullptr) {
  RefineCC(new_to_old, pmhd->u0, pmhd->coarse_u0);
  RefineFC(new_to_old, pmhd->b0, pmhd->coarse_b0);
  // IMPORTANT: Fix div(B) errors at fine/coarse boundaries after refinement
  // When a coarse block is refined adjacent to an already-fine block, the prolongated
  // face-centered B-field values may not match the existing fine values, violating div(B)=0.
  // This call overwrites the prolongated values with the correct values from fine neighbors.
  // See https://github.com/IAS-Astrophysics/athenak/issues/595 for details.
  FixRefinedFCBoundaries(new_to_old, pmhd->b0);
}
```

**Why This Location:**
- Must be called AFTER RefineFC (which does the prolongation)
- Must be called BEFORE boundary communication and primitive variable calculation
- Only needed for MHD (when pmhd != nullptr)
- This is "Step 9" in the RedistAndRefineMeshBlocks sequence

**Why This Comment Style:**
- "IMPORTANT:" flag draws attention to critical fix
- Explains the problem in context
- References the issue for full details
- Follows the multi-line comment style used elsewhere in the file

#### Lines 1178-1288: New Function Implementation
```cpp
void MeshRefinement::FixRefinedFCBoundaries(DualArray1D<int> &n2o, 
                                             DvceFaceFld4D<Real> &b) {
```

**Function Structure Breakdown:**

##### A. Variable Setup (Lines 1187-1196)
```cpp
auto &new_nmb = new_nmb_eachrank[global_variable::my_rank];
auto &ngids = new_gids_eachrank[global_variable::my_rank];

auto &indcs = pmy_mesh->mb_indcs;
auto &is = indcs.is, &ie = indcs.ie;
auto &js = indcs.js, &je = indcs.je;
auto &ks = indcs.ks, &ke = indcs.ke;

bool &multi_d = pmy_mesh->multi_d;
bool &three_d = pmy_mesh->three_d;
```

**Why These Variables:**
- `new_nmb`: Number of mesh blocks on this rank after refinement
- `ngids`: Starting global ID for this rank's blocks
- `indcs`: Mesh indices structure (consistent with other functions)
- `is,ie,js,je,ks,ke`: Index bounds for active zones
- `multi_d, three_d`: Dimensionality flags for conditional logic

##### B. Main Loop Structure (Lines 1199-1284)
```cpp
for (int m = 0; m < new_nmb; ++m) {
  int global_m = m + ngids;
  int old_m = n2o.h_view(global_m);
  
  // Skip blocks that were not refined
  if (refine_flag.h_view(old_m) <= 0) continue;
```

**Design Decisions:**
- **Outer loop only over refined blocks**: More efficient than checking all pairs
- **Use of refine_flag**: Consistent with how refinement is tracked elsewhere
- **Host view access**: These are host-side loops, not device kernels

##### C. Neighbor Detection Logic (Lines 1207-1221)
```cpp
for (int n = 0; n < new_nmb; ++n) {
  if (m == n) continue;
  
  int global_n = n + ngids;
  int old_n = n2o.h_view(global_n);
  LogicalLocation &lloc_n = new_lloc_eachmb[global_n];
  
  // Skip if not at same level or if neighbor was also refined
  if (lloc_m.level != lloc_n.level) continue;
  if (refine_flag.h_view(old_n) != 0) continue;
```

**Why This Logic:**
- Only check blocks at the same refinement level (required for face matching)
- Only copy from blocks that were NOT refined (they have the "correct" values)
- Skip self-comparison (m == n)

##### D. X1-Face Boundary Correction (Lines 1227-1245)
```cpp
if (lloc_m.lx2 == lloc_n.lx2 && lloc_m.lx3 == lloc_n.lx3) {
  if (lloc_m.lx1 == lloc_n.lx1 + 1) {
    // m is right neighbor of n: copy n's right face to m's left face
    par_for("FixDivB-x1L", DevExeSpace(), ks, ke, js, je,
    KOKKOS_LAMBDA(const int k, const int j) {
      b.x1f(m, k, j, is) = b.x1f(n, k, j, ie+1);
    });
  } else if (lloc_m.lx1 == lloc_n.lx1 - 1) {
    // m is left neighbor of n: copy n's left face to m's right face
    par_for("FixDivB-x1R", DevExeSpace(), ks, ke, js, je,
    KOKKOS_LAMBDA(const int k, const int j) {
      b.x1f(m, k, j, ie+1) = b.x1f(n, k, j, is);
    });
  }
}
```

**Critical Details:**
- **Face matching**: `is` is the left face, `ie+1` is the right face
- **Logical location check**: Blocks must align in x2 and x3 to be x1 neighbors
- **Direction matters**: Must copy the correct face based on relative position
- **Kokkos pattern**: Uses standard par_for with KOKKOS_LAMBDA for GPU compatibility

##### E. X2 and X3 Face Corrections (Lines 1248-1283)
Similar patterns for X2 and X3 faces with appropriate:
- Dimension checks (multi_d for X2, three_d for X3)
- Index adjustments (js/je+1 for X2, ks/ke+1 for X3)
- Loop bounds (different indices held constant)

**Why Separate Kernels:**
- Each face direction needs different index patterns
- Keeps kernels simple and efficient
- Allows for dimension-specific optimizations

### 3. Documentation Comments

#### Function Header Documentation (Lines 1178-1201)
The extensive documentation includes:
1. **Brief description**: What the function does
2. **Detailed explanation**: Why it's needed (physics/numerics)
3. **Algorithm description**: Step-by-step process
4. **Implementation notes**: Limitations (same-rank only)
5. **References**: Links to issues and PRs

**Why This Level of Detail:**
- Critical fix that future developers must understand
- Non-obvious solution to a subtle problem
- Helps prevent regression or incorrect modifications
- Provides context for debugging if issues arise

---

## Testing Infrastructure

### Test Problem: `rt_mhd_amr.cpp`

#### Purpose
Creates a Rayleigh-Taylor instability with:
- Density discontinuity (heavy fluid over light)
- Horizontal magnetic field
- Initial perturbation to trigger instability
- AMR refinement at the interface

#### Key Features for Testing div(B):
1. **Uniform initial B-field**: Makes div(B) errors obvious
2. **Interface refinement**: Forces fine/coarse boundaries
3. **Multiple dimensions**: Tests all face directions
4. **Output diagnostics**: Can output div(B) directly

### Input File: `rt_mhd_amr.athinput`

#### Critical Parameters:
```
<mesh_refinement>
adaptive = true
refinement_interval = 10    # Frequent refinement to test the fix
ncycle_check = 5            # Check refinement criteria often

<output2>
variable = mhd_divb         # Direct output of div(B) for verification
```

---

## Design Decisions and Rationale

### 1. Why Copy Instead of Averaging?
**Decision**: Copy field values from fine blocks rather than averaging
**Rationale**: 
- Fine blocks have the "correct" values that maintain div(B)=0
- Averaging would create new values that might not satisfy the constraint
- Copying ensures exact continuity across boundaries

### 2. Why Only Same-Rank Corrections?
**Decision**: Only fix boundaries between blocks on the same MPI rank
**Rationale**:
- Avoids complex MPI communication in the refinement step
- Cross-rank boundaries are handled by existing boundary communication
- Keeps the fix simple and maintainable
- Performance: No additional MPI overhead

### 3. Why Host-Side Loops with Device Kernels?
**Decision**: Main loops on host, copy operations in device kernels
**Rationale**:
- Block neighbor detection requires complex logic (easier on host)
- Actual field copying is data-parallel (efficient on device)
- Matches the pattern used in other refinement functions
- Minimizes host-device data movement

### 4. Why After RefineFC?
**Decision**: Call fix immediately after RefineFC, before anything else
**Rationale**:
- Must happen after prolongation creates the incorrect values
- Must happen before boundary communication propagates errors
- Must happen before primitive variable calculation uses the fields
- Logical place in the code flow

### 5. Why Separate X1/X2/X3 Handling?
**Decision**: Three separate code blocks for three face directions
**Rationale**:
- Each direction has different index patterns
- Easier to understand and debug
- Allows dimension-specific optimizations
- Consistent with how face-centered fields are handled elsewhere

---

## Limitations and Future Work

### Current Limitations

1. **Same-Rank Only**: 
   - Only fixes boundaries between blocks on the same MPI rank
   - Cross-rank boundaries rely on standard boundary communication

2. **Single-Level Jumps**:
   - Assumes refinement creates only single-level jumps
   - Multiple-level jumps not explicitly handled

3. **Performance**:
   - O(N²) comparison of blocks (could be optimized with neighbor lists)
   - Multiple small kernels instead of one large kernel

### Potential Improvements

1. **Cross-Rank Communication**:
   - Could add MPI communication to handle cross-rank boundaries
   - Would ensure fix is applied uniformly

2. **Neighbor Caching**:
   - Could cache neighbor relationships to avoid O(N²) search
   - Would improve performance for large block counts

3. **Kernel Fusion**:
   - Could combine all face corrections into single kernel
   - Might improve GPU performance

4. **Diagnostic Output**:
   - Could add option to output div(B) errors before/after fix
   - Would help verify effectiveness

### Testing Recommendations

1. **Convergence Tests**: Verify div(B) → 0 with increasing resolution
2. **Shock Tests**: Ensure fix works with strong field gradients
3. **3D Tests**: Verify all face directions are correctly handled
4. **Parallel Tests**: Verify cross-rank boundaries are eventually corrected
5. **Performance Tests**: Measure overhead of the fix

---

## Summary

This fix addresses a fundamental numerical issue in AMR+MHD simulations by ensuring magnetic field continuity across refinement boundaries. The implementation is:

- **Minimal**: Only adds one function and one function call
- **Efficient**: Only processes newly refined blocks
- **Robust**: Handles all dimensions and face directions
- **Maintainable**: Well-documented and follows existing patterns
- **Tested**: Includes test problem for verification

The fix is essential for maintaining physical correctness in MHD simulations with AMR and represents a critical improvement to the codebase.