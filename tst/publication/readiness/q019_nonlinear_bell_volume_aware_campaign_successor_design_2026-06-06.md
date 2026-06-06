# Q-019 Volume-Aware Nonlinear Bell Campaign Successor

Date: 2026-06-06

Status: **source-local design only; all physical pilots, production execution,
qualification, and publication claims remain prohibited**

This design supersedes the historical Q-019 nonlinear Bell campaign for every
physical use. The historical JSON and Markdown remain byte-preserved
chronology. Their volume-blind `PPC*deposit_qscale=2e9` rule, artificial-`C`
current normalization, time windows, and resource estimates are not inherited.

The authoritative adversarial audit changes the predecessor order. Q-019
depends first on the **passed raw cycle-one**
`Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE` oracle. A Q-023 source-local or
cycle-zero oracle is not a substitute. The exact Q043 pass receipt, raw
cycle-one inventory, independent recompute, and review bindings must enter the
future Q-019 execution preregistration. Only after that oracle passes may
Q-019 depend on the separately corrected
`Q023-PAPER-BELL-LINEAR-JOVERC` linear Bell predecessor.

## Required Current Closure

Every corrected deck must satisfy and directly verify from raw `prtcl_j`:

```text
PPC * deposit_qscale * species_charge * v_CR / V_root_cell
  = J_CR/c
  = 2 B0 k0
```

Here `V_root_cell` includes all three root-grid axes, including invariant
carrier axes. `deposit_qscale` is therefore recomputed for every root-cell
volume and PPC. The artificial CR light speed is absent from the target and
from the `deposit_qscale` derivation. The general closure uses the configured
`species_charge` coefficient. Writing that coefficient as bare `q/(mc)` is
allowed only when `species_mass=1` is explicitly bound.

No nonlinear physical pilot may run until the passed Q043 raw cycle-one oracle
is exactly bound and the separately corrected
`Q023-PAPER-BELL-LINEAR-JOVERC` predecessor then closes growth, signed phase,
polarization, convergence, retained provenance, independent recompute, and
review.

## Separate Scientific Branches

The first branch is
`Q019-HR-JOVERC-FIXED-CURRENT-LIKE-NOHALL`: a self-consistent,
high-rigidity, periodic-box experiment that may be called
*fixed-current-like* only if the lab-frame CR current and CR momentum remain
invariant while the measured gas-frame current
`J_CR,gas = J_CR,lab - rho_CR u_gas` declines consistently with gas
acceleration. It is not an exact locked-current experiment.

The high-rigidity cold centered beam uses **PPC=1** in production. Higher PPC
only duplicates colocated, identical phase-space samples while
`deposit_qscale` decreases inversely with PPC. It increases cost without
sampling a distribution, so it is useful for the deposition oracle and
debugging but is not a nonlinear particle-convergence test.

The second branch is
`Q019-FR-JOVERC-SELF-CONSISTENT-NOHALL`: a separately named finite-rigidity
campaign with a reviewed distribution function, nonduplicated particle
sampling, a meaningful PPC ladder, and a reviewed initial CR momentum-flux
tensor. It may test CR deflection and anisotropic momentum-flux saturation
only after its own source, mapping, convergence, and evidence gates close.
Its periodic box is an undriven analogue, not a reproduction of driven CR
boundaries or a kinetic-background-ion calculation.

## Geometry And Numerical Controls

All corrected nonlinear decks align `B0`, `J_CR`, and `x1`. The primary
ensemble and numerical sensitivities use 2D3V. Any 3D saturation claim
requires a paired, elongated, axis-aligned fiducial with
`L1 > L2 = L3` and a control that doubles every extent at the same cell size,
seed, solver, physics, and terminal normalized time. A missing member,
fundamental-box-mode saturation, or box-sensitive mechanism classification
fails the 3D claim.

The fiducial numerical method is RK2, PLM, LLF, and Boris--TSC. Required
paired sensitivities change only the MHD Riemann solver (`LLF -> HLLD`), only
the reconstruction (`PLM -> WENOZ`, with its required ghost-zone change),
the cell size, or the timestep. Each alternative must first pass source-local
runtime, restart, and output compatibility.

## Diagnostics And Residuals

Retained evidence must include raw MHD fields, raw deposited particle moments,
sparse particle state, and schema-7 restarts. Analysis must report lab-frame
and gas-frame current, CR momentum and kinetic energy, full CR
momentum-flux/pressure tensors, sampled-field gyroradii, gas acceleration,
relative drift, magnetic spectra, dominant scale, and filament/cavity
statistics.

Energy partitions include parallel and transverse gas bulk kinetic energy,
gas thermal energy, guide and transverse magnetic energy, total MHD energy,
CR kinetic energy, and cumulative CR-to-MHD transfer. Conservation residuals
must be reported against several denominators separately: initial total
energy, cumulative transferred energy, MHD energy change, and magnetic-energy
gain, with analogous momentum normalizations. A single favorable denominator
is insufficient.

## Pilots, Seeds, And Resources

Pilot and qualifying seed namespaces are disjoint. Seeds are paired across
every required sensitivity and box control. The exact qualifying inventory,
seed count, geometry, terminal time, output cadence, saturation algorithm,
numeric tolerances, and resource allocation freeze only after excluded pilots
and before any qualifying output is inspected.

The project-wide Frontier cap remains 10,000 node-hours. This design allocates
none of it. A future registered pilot preregistration must impose a hard cap
no larger than the lesser of 500 node-hours or ten percent of then-unreserved
project budget. Production authorization must include complete paired groups,
the 3D box pair when claimed, Q-011 obligations, and a failure reserve. Budget
pressure may trigger a reviewed redesign; it may not silently remove seeds,
sensitivities, or the 3D box control.

## Literature Boundary

Bell (2004) supplies the nonresonant current-driven context. The
high-rigidity branch tests the plasma-acceleration and current-frame
discriminants associated with the constant-current limit of Riquelme &
Spitkovsky (2009), without claiming locked particles. The finite-rigidity
branch can test self-consistent current evolution and isotropization in the
spirit of Gargate et al. (2010), but AthenaK's ideal-MHD background cannot
test kinetic-ion demagnetization. It can test an undriven analogue of the
anisotropic momentum-flux relation of Zacharegkas et al. (2024), but not
their driven-boundary setup. Sun & Bai (2023) remains the linear AthenaK
context only through the corrected, volume-aware predecessor.
