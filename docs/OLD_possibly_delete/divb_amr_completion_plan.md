# Face-Field Correction Implementation Checklist

## Phase 1: Fix Foundation Issues

### [x] 1.1 Fix BuildLocationMaps() Method
**File:** `src/mesh/mesh.cpp`
- [x] ~~Fix pointer assignment to use array indexing~~ (Not needed - MeshBlock is single object with arrays)
- [x] Understood MeshBlock structure - single object with DualArrays
- [x] BuildLocationMaps stores base pointer correctly

### [x] 1.2 Fix Neighbor Access Pattern
**File:** `src/mesh/mesh.cpp` in `ApplyFaceFieldCorrection()`
- [x] Fix neighbor array access pattern
- [x] Use proper 2D array indexing: `nghbr.h_view(m, n)`
- [x] Ensure gid access uses correct array index
- [x] Fixed newly_created to be DualArray1D<bool> instead of single bool

## Phase 2: Implement Face Orientation Detection

### [x] 2.1 Add Face Type Detection Logic
**File:** `src/mesh/mesh.cpp` in `ApplyFaceFieldCorrection()`
- [x] Implement face type detection using neighbor indices
- [x] Added get_face_info lambda function
- [x] Determine inner vs outer face correctly

### [x] 2.2 Calculate Buffer Sizes Correctly
- [x] For X1 faces: `buf_size = (je-js+1) * (ke-ks+1)`
- [x] For X2 faces: `buf_size = (ie-is+1) * (ke-ks+1)`
- [x] For X3 faces: `buf_size = (ie-is+1) * (je-js+1)`
- [x] Handle 2D case properly

## Phase 3: Add Kokkos Memory Management

### [x] 3.1 Create Host Mirrors for Magnetic Fields
**File:** `src/mesh/mesh.cpp`
- [x] Add host mirror declarations for all three face fields
- [x] Add deep_copy before packing to get device data to host

### [x] 3.2 Add Deep Copy After Unpacking
- [x] After modifying host data, copy back to device
- [x] Deep copies for all three face field components

## Phase 4: Implement Data Packing (Coarse Side)

### [x] 4.1 Pack X1 Face Data
**Location:** Phase 2 of `ApplyFaceFieldCorrection()`
- [x] Determine face index correctly
- [x] Pack data from host mirror with proper indexing
- [x] Send via MPI_Isend with correct tags

### [x] 4.2 Pack X2 Face Data
- [x] Determine face index for X2 faces
- [x] Pack b0.x2f data with correct ordering

### [x] 4.3 Pack X3 Face Data
- [x] Determine face index for X3 faces  
- [x] Pack b0.x3f data with correct ordering

## Phase 5: Implement Data Unpacking (Fine Side)

### [x] 5.1 Unpack X1 Face with 2×2 Replication
**Location:** Phase 3 of `ApplyFaceFieldCorrection()`
- [x] Calculate fine face index
- [x] Apply 2×2 replication pattern
- [x] Track receive buffer correspondence with RecvInfo struct

### [x] 5.2 Unpack X2 Face with 2×2 Replication
- [x] Implemented pattern for X2 faces
- [x] Correct index mapping from coarse to fine

### [x] 5.3 Unpack X3 Face with 2×2 Replication
- [x] Implemented pattern for X3 faces
- [x] Correct index mapping from coarse to fine

### [x] 5.4 Deep Copy Back to Device
- [x] Added Kokkos::deep_copy for all three face fields
- [x] Updates device arrays after unpacking

## Phase 6: Update Prolongation Operators

### [ ] 6.1 Add Coarse Neighbor Detection
**File:** `src/mesh/prolongation.hpp`
- [ ] Add function to check if face has coarse neighbor
- [ ] Pass neighbor information to prolongation functions

### [ ] 6.2 Modify ProlongFCSharedX1Face
- [ ] Add parameter for skip flag
- [ ] Check if face should be skipped (has coarse neighbor)
- [ ] Return early if face was corrected by FFC

### [ ] 6.3 Modify ProlongFCSharedX2Face
- [ ] Similar modifications as X1

### [ ] 6.4 Modify ProlongFCSharedX3Face
- [ ] Similar modifications as X1

## Phase 7: Handle Edge Cases

### [ ] 7.1 Handle Non-MPI Case
- [ ] Implement local (same-rank) transfers without MPI
- [ ] Direct memory copy for local neighbors

### [ ] 7.2 Handle 2D Problems
- [ ] Check if problem is 2D (ke == ks)
- [ ] Adjust replication pattern accordingly

### [ ] 7.3 Handle Domain Boundaries
- [ ] Check for physical boundary conditions
- [ ] Skip FFC for domain boundary faces

## Verification Points

### After Phase 1:
- [x] Code compiles without errors
- [x] BuildLocationMaps creates correct mappings

### After Phase 3:
- [x] Host mirrors created successfully
- [x] Deep copies execute without errors

### After Phase 4:
- [x] Face data packed correctly
- [x] MPI sends/receives posted properly

### After Phase 5:
- [x] Face values replicated correctly
- [x] Device arrays updated
- [x] **COMPILATION SUCCESSFUL**

### After Phase 6:
- [ ] Prolongation skips corrected faces
- [ ] No double-correction of faces

## Notes and Issues

### Current Issues Found:
1. BuildLocationMaps stores same pointer for all blocks
2. Neighbor array access pattern incorrect
3. Missing Kokkos host mirrors
4. Face orientation detection incomplete

### Key Physics Requirements:
- Magnetic flux conservation: Φ_coarse = Σ Φ_fine
- Normal B continuous across interfaces
- div(B) = 0 preserved to machine precision

### Implementation Strategy:
- Fix foundation issues first (Phase 1)
- Add memory management infrastructure (Phase 3)
- Implement core logic (Phases 4-5)
- Integrate with existing code (Phase 6)
- Handle special cases (Phase 7)