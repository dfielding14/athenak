# Passive CGL (Task4)

The standalone canonical candidate passed **240 application cases**: 80 each
on CPU, HIP1 and HIP4, with 65 assertions per backend. Each backend compared
full-precision flow/face-field bits, timestep sequences and next timestep;
all five reconstructors, floors, FOFC, weak fields, forcing and LF were covered.
Independent heating, linear modes and material advection passed their
refinement gates. Forced collisional and zero-rate restarts were byte exact.

Final release gates use the complete Task0–Task7 composition. Their packaged
and original/3D paper outcomes are listed in `summary.json`; pending entries
must not be interpreted as completed acceptance. The full source/binary
manifests under `WO2/final-research` identify the integrated executable.

## Failures retained as evidence

The first HIP flow attempt failed exact equality. The same-input actual HIP
face helper produced 2,258 mismatches despite CPU exactness. The corrected
candidate runs the exact native isothermal kernel first, then thermal-only
scalar fluxes. It has zero mismatches among 18,432 HIP flow values; the legacy
CGL negative control has 11,650 mismatches. `hip-face-failure.log` preserves
the original failure and the helper result JSON preserves the passing control.

The first forced HIP restart preserved flow, B, time and dt but changed 63,488
full-state thermal bits. Canonical C2P after both collision wrappers fixed it;
all final CPU1/HIP1/HIP4 restarts then passed exactly. The original failure
log is retained here; full states and intermediate probes remain in the
scratch delivery archive. No tolerance was relaxed.

Paper repeatability initially exposed uninitialized restart metadata. The raw
failed checker JSON is preserved beside an adjudicated JSON. Exact physical
comparison uses the pre-existing validated restart contract: unused root mesh
coarse indices, dormant startup RNG fields only at time/cycle zero with the
verified seed/cache state, and four non-member RNG alignment bytes. All live
RNG members, diagnostics and complete state remain compared; exact excluded
byte ranges and original hashes are recorded per file.

The first final HIP4 packaged run used an incomplete Frontier runtime
contract and faulted on GPU memory access during forcing. Exclusive probes
reproduced this with both fused and unfused executables, so fusion was not
necessary. Restoring the exact accepted HSA/NIC/IPC/FI environment gave six
consecutive unchanged-binary probe passes and a fresh full HIP4 qualification.
The failed runs, raw environment dumps, node/binary identity and control
states are retained. The specific responsible flag or runtime mechanism is
not claimed to be isolated. No solver or tolerance workaround was applied.

Root's final compatibility followups add explicit Real casts (identity for
the accepted double build), extract the literal unchanged pressure-traction
helper to avoid an unrelated GR float include dependency, and link the actual
Kokkos runtime in focused tests. Root's followup evidence records those
source/build-only changes and double-arithmetic identity proofs.

The original paper deck has planar legacy forcing and cannot validate LF
activity. It is retained unchanged as the requested input. A separately
labeled 3D forcing policy tests nonzero LF work. Additional output blocks and
explicit unchanged defaults support provenance. Concurrent run times are
functional evidence only and are not used as performance measurements.

## Scope

See [the passive model](../../../source/modules/cgl_passive.md) for J/A,
physical thermal U, LF U/mu, initialization and restart contracts. Scope is
uniform periodic Newtonian HLLE flow with LF, collisions, limiters and
forcing. Unsupported source/mesh/boundary combinations and unredefined SGS
energy moments have explicit fences. Historical paper-scale analysis catalogs
remain intact; unvalidated passive analysis consumers retain their own fence.

The full-precision inputs, logs and restart states are retained under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/task4-research` and bundled
in its delivery archive. `archive-manifest.json` identifies that artifact.
