# Task 6: representation validation and conservative shared LF fluxes

Triage: **adapted**. Existing A/mu-aware primitive prolongation is verified on
the GPU; a separately demonstrated shared-face conservation gap is corrected
and covered by permanent CPU/GPU/MPI regressions.

The existing primitive-prolongation implementation correctly handles A and mu on
an actual GPU. The validation also exposed an independent baseline conservation
bug: separately prolonged corner magnetic ghosts produce different transverse
LF flux estimates on two copies of the same physical face. A separate correction
now synchronizes same-level IEN/IAN fluxes during multilevel LF stages. It retains
the existing coarse/fine correction, projections, magnetic exchanges, and physical
boundary fills. Task3's earlier bitwise proof remains a separate result.

## Baseline defect and cause

[baseline-defect.md](baseline-defect.md) records the original diagnosis. In a periodic,
zero-velocity, frozen-field 32-cell-root smooth SMR case at time 0.002, the old
baseline changes integrated energy by -4.445033718880609e-8 and mu by
-3.0149784668864754e-8 with conserved prolongation. Primitive prolongation changes
them by -1.0930872873515796e-8 and +2.8896616210971615e-9. The immutable baseline and
physical-refill-only reference give identical normalized restart states. Every
LF and AMR repair counter stays zero.

The stage ledger first locates the drift at STSUpdateU. Active face B does not
change. Pairing shared faces then isolates the residual at same-level faces near
a refinement corner; coarse/fine corrected faces cancel to rounding. Direct
stencil inspection finds different reconstructed Bcc and pressures in ghost
copies of the same transverse physical cell. These cells enter the LF stencil.

## Correction

The generic cell-centered flux API accepts an optional same-level variable range;
ordinary callers keep it disabled. Multilevel LF stages select IEN/IAN. Both
original face estimates are packed before unpacking. Identical estimates remain
unchanged; differing estimates use a symmetric, overflow-safe half-plus-half
mean with the operands ordered identically on both sides, including under FMA
contraction. No new runtime switch is added.

The existing communicator and send/receive completion slots cover these faces.
Only the requested full-face buffer columns grow when necessary, before posting
receives; 3D regression coverage exercises this larger capacity. The existing
mesh startup guard rejects multilevel shearing boxes, so remapped shear faces
cannot enter this path. Supporting refined shear remains outside this change.

## CPU, GPU, MPI and stage checks

All 24 standalone regression cases pass: CPU and HIP, one and four MPI ranks,
2D conserved/primitive prolongation with both STS and explicit LF, and 3D
conserved/primitive prolongation with STS. One/four-rank final states agree
exactly. Across these cases the largest absolute integrated energy or mu drift
is 4.440892098500626e-16. Every LF and AMR repair counter is zero; active density,
momentum, scalar, and Bcc remain unchanged; the thermal state evolves.

The combined Task0/1/3/6/5 scratch binaries also pass the same 24 CPU/HIP and
one/four-rank checks, with maximum absolute E/mu drift 4.440892098500626e-16 and
exact one/four-rank state agreement. The tested CPU SHA256 is
`ec446b3c50935cfe1d11b2b563d820395cba7ee1af71dcc2186bc80535ba1eee`; HIP is
`f47cea69f5d454dfb9bcc62913140a15bdc993a48c61c7ec82caba9e34167775`.

The permanent checker retains a 5e-12 absolute conservation bound for both
energy and mu, verifies finite positive states and divergence control, and
compares full-precision restart states. The old baseline fails this checker.
The tolerance has not been relaxed.

The corrected HIP trace covers conserved and primitive prolongation across ten
stages each. Same-level shared faces cancel exactly, and the first coarse/fine
face-sum residual is at most 4.4e-23. Across all twenty post-receive stages the
largest global weighted E/mu flux divergence is 2.8587361969832636e-21. The
stage-integrated energy stays identical; mu changes at most 2.220446049250313e-16.
Active face B never changes.

Independent thermodynamic checks include the active cells and one-cell halo at
Begin, primitive refresh, and End, using the representation recorded at that
stage. All 48 snapshots pass A/mu, energy, momentum, scalar, pressure-admissibility,
and face-to-cell-B consistency checks; the maximum normalized residual is
6.293385664684257e-16 against the unchanged 2e-12 bound.

## Smooth primitive/conserved agreement

Actual HIP runs use root resolutions 32, 64, and 128 with the same physical
refinement geometry and matched final time 0.002. Both prolongation choices start
with identical active states. The comparison uses full-precision restart data.

| Quantity | Initial fine-ghost mean order 32→64 / 64→128 | Final active volume-L1 order 32→64 / 64→128 |
| --- | ---: | ---: |
| Energy | 1.948 / 1.975 | 1.933 / 1.895 |
| mu | 1.889 / 1.952 | 1.812 / 1.992 |
| Parallel pressure | 1.900 / 1.955 | 1.893 / 1.951 |
| Perpendicular pressure | 1.886 / 1.950 | 1.815 / 1.992 |

Final energy volume-L1 differences are 5.267512582403898e-5,
1.37962865030654e-5, and 3.710170344601316e-6. Final mu differences are
4.6684066647311384e-5, 1.3290974846510119e-5, and 3.341866941542884e-6.
These are approximately second-order integral agreement between the two
prolongation choices, not an exact analytic solution error. Max-norm agreement
is not uniformly second order: final energy max-norm orders are 0.553 and 0.664,
and final fine-ghost mean energy orders are 1.837 and 1.526. All raw norms are
retained; no stronger convergence claim is made. Every resolution passes the
independent direct-conservation, pressure, frozen-field and complete repair-
counter audit. The largest direct integral drift across all six smooth runs is
2.886579864025407e-15.

## Provenance and archive

The [artifact index](artifact-index.json) preserves absolute raw paths and
SHA256 hashes. Raw outputs remain under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/p1-research`.
The archive retains all 48 permanent-check records, [causal trace summaries](causal-trace-summary.json),
[negative baseline controls](negative-control.json), [all convergence norms](smooth-convergence.json),
and the [complete smooth health audit](smooth-health-conservation.json).
The negative control uses the permanent test's exact restart reader and its
unchanged 5e-12 bound. It demonstrates rejection of the original physical defect.

Isolated links include scratch forwarding overloads for old baseline callers;
the combined full builds use the production signatures. Each original build
manifest, referenced by the artifact index, retains complete source/object hashes
and compile/link commands. No scratch wrapper or instrumentation enters production.
The following paths are relative to the retained raw evidence root, unless
provided as files in this archive.

## Review artifacts

- Production patch after Task0 refill: `task6-sync.after-refill.patch`.
- Production patch after Task3: `task6-sync.after-task3.patch`.
- Permanent regression patch: `task6-sync.regression.patch`.
- Standalone CPU/HIP evidence: `runs/task6-sync-{cpu,hip}/results.json`.
- Combined candidate evidence: `runs/task6-integration-{cpu,hip}/results.json`.
- Corrected trace evidence: `runs/task6-sync-trace-hip/*/{stage-ledger,face-balance,representation}.json`.
- Smooth data, all error norms and health audit:
  `runs/task6-sync-convergence-hip/{convergence,health-conservation}.json`.
- Each build manifest and run command records source/binary hashes and commands.

These validation runs were allowed to share the allocation and are not timing
measurements. This correctness correction makes no throughput claim.
