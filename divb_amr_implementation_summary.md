# Face-Field Correction for AMR div(B) Implementation Summary

## Implementation Status: PARTIAL SUCCESS ✓

### Completed Steps

#### 1. Added LogicalLocationHash and Equality Operator ✓
**File:** `src/mesh/mesh.hpp`
- Added `#include <unordered_map>`
- Implemented `LogicalLocationHash` struct with hash function
- Added equality operator for `LogicalLocation`
- Enables O(1) neighbor lookups using hash maps

#### 2. Added Tracking Structures ✓
**Files:** `src/mesh/mesh.hpp`, `src/mesh/meshblock.hpp`
- Added `bool newly_created` flag to MeshBlock class
- Added `location_maps_` vector of hash maps to Mesh class (per refinement level)
- Added public methods `BuildLocationMaps()` and `ApplyFaceFieldCorrection()`

#### 3. Implemented BuildLocationMaps() Method ✓
**File:** `src/mesh/mesh.cpp`
- Builds per-level hash maps for fast neighbor lookup
- Maps LogicalLocation → MeshBlock pointer
- Called after mesh structure changes

#### 4. Implemented ApplyFaceFieldCorrection() Method (Skeleton) ✓
**File:** `src/mesh/mesh.cpp`
- Three-phase implementation structure:
  - Phase 1: Post receives on fine blocks
  - Phase 2: Pack and send from coarse blocks
  - Phase 3: Wait and apply corrections
- MPI and non-MPI versions included
- **Note:** Current implementation is a skeleton - face data packing/unpacking needs refinement

#### 5. Added MPI Infrastructure ✓
**File:** `src/mesh/mesh_refinement.hpp`
- Added `static constexpr int ffc_tag = 0x50` for FFC MPI messages
- Added face-field correction buffers (send/receive)
- Added MPI request vectors for asynchronous communication

#### 6. Modified AMR Flow ✓
**File:** `src/mesh/mesh_refinement.cpp`
- Added `BuildLocationMaps()` call after refinement
- Added `ApplyFaceFieldCorrection()` call for MHD problems
- Added logic to mark newly created blocks using `newtoold` array

#### 7. Compilation Test ✓
- Code compiles successfully with blast problem generator
- No compilation errors introduced

### Incomplete/TODO Items

#### 1. Complete Face-Field Correction Implementation
The current `ApplyFaceFieldCorrection()` is a skeleton that needs:
- Proper identification of face orientations (X1, X2, X3)
- Actual magnetic field data packing from device to host
- Correct unpacking and application to fine faces (2×2 replication)
- Kokkos deep_copy for device/host transfers

#### 2. Update Prolongation to Skip Corrected Faces
**File:** `src/mesh/prolongation.hpp`
- Need to modify `ProlongFCSharedX1Face`, `ProlongFCSharedX2Face`, `ProlongFCSharedX3Face`
- Add checks to skip faces that have been corrected by FFC
- Requires passing neighbor information to prolongation functions

#### 3. Testing and Validation
- Test with actual AMR problems (e.g., field loop with refinement)
- Verify div(B) = 0 across refinement boundaries
- Test MPI parallel execution
- Benchmark performance impact

## Key Technical Achievements

1. **Infrastructure Ready:** All necessary data structures and method signatures are in place
2. **AMR Integration:** Face-field correction is properly integrated into the AMR workflow
3. **Compilation Clean:** No errors introduced, code builds successfully
4. **MPI Framework:** Asynchronous communication structure is established

## Next Steps for Full Implementation

1. **Complete FFC Logic:** Fill in the actual magnetic field data handling in `ApplyFaceFieldCorrection()`
2. **Neighbor Detection:** Properly identify coarse-fine interfaces using neighbor information
3. **Face Data Transfer:** Implement correct packing/unpacking based on face orientation
4. **Prolongation Updates:** Modify prolongation operators to respect corrected faces
5. **Extensive Testing:** Validate with standard MHD AMR test problems

## Physics Impact

When fully implemented, this will:
- Maintain div(B) = 0 to machine precision across AMR boundaries
- Preserve magnetic flux conservation in refined regions
- Enable robust MHD simulations with adaptive mesh refinement
- Match the approach used in Athena++ PR #625

## Code Quality

The implementation follows AthenaK conventions:
- Uses Kokkos for performance portability
- Respects existing MPI communication patterns
- Maintains separation between host and device data
- Integrates cleanly with existing AMR infrastructure