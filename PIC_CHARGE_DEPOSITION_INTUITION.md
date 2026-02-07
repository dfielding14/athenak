# PIC Charge Deposition: Physics Intuition First

## Why charge deposition exists
In particle-in-cell (PIC), particles carry charge and move continuously, while
fields live on a mesh. Charge deposition is the bridge between those two
pictures.

At each time step, particles must tell the grid:
1. How much charge density `rho` exists in each region.
2. How much current density `J` is flowing through the grid.

Once `rho` and `J` are on the mesh, field updates can respond self-consistently.

## The non-negotiable physics idea: local charge conservation
The key physical constraint is the continuity equation:

`d(rho)/dt + div(J) = 0`

Interpretation: charge cannot appear or disappear; it can only move.

If deposition violates this relation (even slightly and repeatedly), the field
solver must absorb inconsistency, which can show up as noise, wrong wave content,
or long-time drift.

## What "shape clouds" mean physically
A particle is not represented as an infinitesimal point on the grid. Instead,
its charge/current is distributed with a finite shape (a "cloud").

Physical intuition:
1. The grid only resolves finite scales, so sub-cell point structure is not
   physically meaningful at mesh resolution.
2. A finite shape reduces grid aliasing and unphysical high-k noise.

Common shapes (increasing smoothness/support):
1. NGP (nearest-grid-point): all weight to one cell/node.
2. CIC (cloud-in-cell): linear split over adjacent cells/nodes.
3. TSC (triangular-shaped cloud): wider, smoother quadratic weighting.

## Is TSC the "industry standard"?
Short answer: it is common, but not universal.

In practice, production PIC codes use different combinations depending on goals:
1. Low cost / robustness: CIC-like schemes.
2. Better noise properties: TSC or higher-order shapes.
3. Strict conservation priority: charge-conserving current deposition algorithms
   (for example, Esirkepov-style trajectory-based deposition), often paired with
   a chosen shape order.

Important distinction:
- "TSC" describes the interpolation/deposition shape order.
- "Charge-conserving trajectory deposition" describes how current is deposited so
  continuity is satisfied discretely.

These are related but not the same decision.

## Why direct staggered, trajectory-based current deposition matters
For electromagnetic PIC, current naturally couples to staggered field locations
(in Yee-like layouts). Depositing current directly to those staggered locations
using particle trajectory information (old -> new position) improves physical
consistency:
1. Better discrete charge conservation.
2. Cleaner coupling to field updates.
3. Less reliance on post-hoc conversions.

## Current AthenaK approach (today)
AthenaK currently uses a transitional approach:
1. Deposit `rho/J` in a cell-centered particle-moment array.
2. In `edge_staggered` coupling mode, convert cell-centered `J` to edge values
   using deterministic local averaging.
3. Use those values in coupled E-field source updates.

This is intentionally incremental and testable, but it is not yet full
Entity-style direct staggered charge-conserving deposition.

## Why not jump straight to Entity-style everywhere immediately?
Mostly engineering risk control, not a rejection of the physics.

A full direct-deposition change touches multiple correctness-sensitive systems at
once:
1. Restart fidelity for coupled moment/edge-current state.
2. AMR/multilevel handling for deposited quantities.
3. Non-periodic boundary policy for deposited moments/currents.
4. Task ordering and MPI decomposition invariance.

AthenaK has been sequencing these as separate gates so each source of error can
be isolated, tested, and validated.

## Where this is heading
The planned PR4 direction is to add Entity-style direct staggered,
charge-conserving trajectory deposition as an explicit path, with side-by-side
validation against the current converted-edge path before any default change.

That keeps physical fidelity moving forward without losing determinism and
regression confidence.

## Mental model summary
1. Shape order (CIC/TSC/...) controls smoothness/noise and effective particle
   footprint.
2. Charge-conserving current deposition controls whether continuity is respected
   discretely.
3. High-quality PIC needs both choices to be coherent with the field layout and
   time update.
