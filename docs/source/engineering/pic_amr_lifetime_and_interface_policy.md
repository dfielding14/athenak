# PIC AMR Lifetime And Interface Policy

This document records the bounded local Q-009 audit for particle-bearing AMR.
It is an implementation inventory and policy decision, not a claim that the
full AMR qualification matrix has passed.

## Lifetime Inventory

`MeshBlockPack` is retained across adaptive reconstruction. Its child
`MeshBlock` and coordinate objects are deleted and rebuilt, then the rebuilt
`MeshBlock` neighbor table is regenerated from the current tree and rank map.
Physics-module pointers are preserved while those children are reconstructed.

| Retained surface | Local handling after AMR reconstruction | Residual evidence requirement |
| --- | --- | --- |
| `MeshBlockPack` identity | Retained; rebuilt child pointers are installed in place | Repeated refine/derefine under host memory checking and Frontier HIP |
| Child `MeshBlock` and coordinates | Reconstructed; `SetNeighbors` rebuilds neighbor state | Boundary and MPI migration stress |
| MHD module | `MHD::UpdateAfterAMR` refreshes the retained pack view | Coupled AMR continuation matrix |
| Particle module | `Particles::UpdateAfterAMR` refreshes the retained pack pointer and validates MeshBlock-sized moment, coarse-moment, edge-current and no-MHD field capacities | Repeated coupled AMR, restart and GPU lifetime stress |
| Particle position history | Resized when needed after AMR so direct trajectory-current deposition can retain old-position storage | Direct-staggered extension stress |
| Particle boundary helper | Retains the stable `Particles` object; particle send/receive buffers are resized from the current record layout and current message counts | MPI migration and decomposition changes |
| Moment and edge-current boundary helpers | Retain the stable pack pointer and allocate MeshBlock request/buffer capacity against `max_nmb_per_rank` | Repeated multilevel communication stress |
| AMR and load-balance request arrays | AMR transfers are cleared before child reconstruction; load-balance particle request vectors are transient to each migration operation | MPI interruption and sanitizer runs |
| Particle-aware load costs | Post-AMR geometric ownership assigns particle-count contributions before balancing | Correlate selected costs with measured Frontier time |

This inventory does not make persistent raw child pointers acceptable. New PIC
or MHD state added behind a MeshBlock-sized view must either be provisioned
against `max_nmb_per_rank`, reconstructed, or refreshed explicitly after AMR.

## Refinement-Interface Deposition Policy

The production paper-mode policy is named `paper_smooth`. In the current
implementation it maps to the existing cell-centered moment path:

1. Deposit cell-centered particle moments.
2. Restrict deposited moments into coarse storage on multilevel meshes.
3. Exchange neighboring values with additive receive accumulation.
4. Fill coarse boundary state and prolongate fine moment ghosts.
5. Use the synchronized cell-centered moments for the paper-mode gas feedback
   path.

This is the locally smooth paper policy. It is not advertised as individually
conservative particle feedback at every refinement-interface crossing.

The optional `conservative` AMR-interface policy is **not retained as a
qualified production mode**. The existing
`couple_j_deposition_mode=direct_staggered` machinery is an experimental
trajectory-current candidate for separately scoped edge-current work. It has
cell-crossing guards and boundary handling, but it does not by itself establish
a conservative AMR gas-feedback policy. Any future `conservative` policy must
be separately named, derived and validated before use in a production claim.

## Current Evidence Boundary

The checked-in AMR shock/load-balance smoke exercises serial refinement and
conditionally exercises MPI migration and particle tracking when an
MPI-enabled runtime is available. The refinement-boundary characterization
records bounded SMR and AMR-proxy behavior. These are useful local guards.

The bounded serial Q-009 successor now alternates six refine/derefine
transitions through five restart continuations while verifying stable particle
identity and refreshed ownership. Q-009 remains open until coupled-boundary,
host memory-checking, controlled MPI and Frontier HIP lifetime evidence is
archived. Section 5.3 AMR scientific qualification also remains a separate
gate.
