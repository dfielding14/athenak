# Face-Field Correction for AMR div(B) Implementation Plan

## Overview
This plan implements face-field correction (FFC) to maintain div(B) = 0 across AMR refinement boundaries in AthenaK. The approach is based on the Athena++ PR #625 and adapted for AthenaK's Kokkos-based architecture.

## Physics Background
When AMR creates fine/coarse boundaries, the magnetic field faces at these boundaries must satisfy flux conservation to maintain div(B) = 0. The face-field correction ensures that the area-averaged magnetic flux through coarse faces equals the sum of fluxes through the corresponding fine faces.

## Implementation Steps

### Step 1: Add LogicalLocationHash and Equality Operator
**File:** `src/mesh/mesh.hpp`
- Add `#include <unordered_map>` 
- Add LogicalLocationHash struct after LogicalLocation definition
- Add equality operator for LogicalLocation
- This enables fast O(1) neighbor lookups using hash maps

### Step 2: Add Tracking Structures to Mesh Classes
**File:** `src/mesh/mesh.hpp`
- Add `bool newly_created` flag to MeshBlock class (line ~50)
- Add `std::vector<std::unordered_map<LogicalLocation, MeshBlock*, LogicalLocationHash>> location_maps_` to Mesh class private members
- Add public methods `BuildLocationMaps()` and `ApplyFaceFieldCorrection()` to Mesh class

### Step 3: Implement BuildLocationMaps Method
**File:** `src/mesh/mesh.cpp`
- Create new method that builds per-level hash maps of LogicalLocation -> MeshBlock*
- Iterate through all MeshBlocks and insert into appropriate level map
- Called after any mesh structure change (refinement/load balancing)

### Step 4: Implement ApplyFaceFieldCorrection Method
**File:** `src/mesh/mesh_refinement.cpp`
- Three-phase implementation:
  1. **Receive posting (fine side):** Post MPI_Irecv for face data from coarse neighbors
  2. **Send packing (coarse side):** Pack and send face-centered B-field to fine neighbors
  3. **Unpacking (fine side):** Wait for receives and apply coarse values to fine faces
- Handle all three face orientations (X1, X2, X3)
- Use host buffers for MPI, then copy to device with Kokkos::deep_copy

### Step 5: Add MPI Tags and Buffers
**File:** `src/mesh/mesh_refinement.hpp`
- Add `static constexpr int ffc_tag = 0x50;` for FFC MPI messages
- Add temporary send/receive buffers for face data exchange

### Step 6: Modify AMR Flow
**File:** `src/mesh/mesh_refinement.cpp` in `AdaptiveMeshRefinement()` method
- After UpdateMeshBlockTree, call BuildLocationMaps()
- Mark newly created blocks with `newly_created = true`
- After redistribution, call BuildLocationMaps() again
- Insert ApplyFaceFieldCorrection() before prolongation
- Sequence: Restrict → FFC → Prolongation

### Step 7: Update Prolongation to Skip Corrected Faces
**File:** `src/mesh/prolongation.hpp`
- Modify ProlongFCSharedX1Face, ProlongFCSharedX2Face, ProlongFCSharedX3Face
- Add checks to skip faces that have coarse neighbors (already corrected by FFC)
- Only prolongate faces that don't border coarse blocks

### Step 8: Integration with MHD Module
**File:** `src/mhd/mhd.cpp` or relevant task file
- Ensure FFC is called in the correct task sequence
- May need to add new task ID if using task-based execution

## Key Technical Considerations

### Memory Management
- Use Kokkos::View for device arrays
- Host mirrors for MPI communication
- Proper deep_copy between host/device

### MPI Communication Pattern
- Non-blocking receives posted first
- Pack/send from coarse blocks
- Wait and unpack on fine blocks
- Handle both local (same rank) and remote (different rank) neighbors

### Face Value Replication
- Coarse face value → 2×2 fine faces (2D)
- Coarse face value → 2×2×2 fine faces (3D)
- Ensures flux conservation

### Neighbor Identification
- Use location_maps_ for O(1) local neighbor lookup
- Check neighbor levels to identify coarse/fine boundaries
- Handle edge cases (domain boundaries, same-level neighbors)

## Testing Strategy
1. Compile with basic changes first
2. Test with simple AMR problem (e.g., field loop)
3. Monitor div(B) across refinement boundaries
4. Verify flux conservation
5. Test with MPI parallelization

## Expected Outcomes
- div(B) remains machine precision across AMR boundaries
- Magnetic flux conservation maintained
- Compatible with existing prolongation/restriction operators
- Minimal performance impact

## Risk Mitigation
- Make changes incrementally
- Test compilation after each major step
- Use existing AthenaK patterns for Kokkos/MPI
- Preserve existing functionality when FFC not needed