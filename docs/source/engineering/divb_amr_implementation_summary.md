# Face-Field Correction for AMR div(B) Implementation Summary

## Implementation Status: COMPLETE ✔️

### Key Features
- **Per-level location maps** map each `LogicalLocation` to its exact `MeshBlock` for fast neighbor lookup.
- **Face-field correction** executes in three phases:
  1. Post non-blocking receives and record metadata for each coarse neighbor.
  2. Pack only the needed coarse-face slices directly from device subviews and send them.
  3. Replicate coarse faces onto the fine mesh on the device and clear `newly_created` flags.
- **Serial fallback** applies the same replication logic without MPI for single-rank runs.
- **AMR divergence test** (`tst/scripts/mhd/divb_amr.py`) with input deck `inputs/mhd/ffc_divb.athinput` ensures ∇⋅B remains near machine precision across refinement boundaries.

### Remaining Work
- Additional divergence monitoring and error handling may further harden the implementation.
