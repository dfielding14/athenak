# Q022 central reference-mapping successor v1

This additive successor maps only the two central paper comparisons:

1. corrected, volume-aware, no-Hall nonlinear Bell; and
2. the selected `p0=1` Q011 parallel shock against Bai et al. (2015) and
   Sun and Bai Section 5.4.

Historical Q022 files remain byte-preserved. The maps use only the locally
retained Bai et al. and Sun-and-Bai PDFs and exact AthenaK readiness/deck
contracts. They grant no launch, execution, qualification, claim, publication,
or policy authority.

## Corrected nonlinear Bell

Sun and Bai Section 5.2, PDF page 9, Equation (27), defines
`k0=j_CR/(2*B_g*c)`. AthenaK raw `prtcl_j` is already deposited `J_CR/c`, so
the corrected target is

```text
PPC * deposit_qscale * species_charge * v_CR / V_root_cell
  = J_CR/c
  = 2 * B0 * k0
```

Artificial particle light speed `C` is not a current-conversion factor. The
target-paper Bell benchmark is linear. It supplies no standalone nonlinear
saturation dataset or numeric nonlinear acceptance tolerance. The corrected
Q019 branches therefore may use the target-paper normalization and linear
predecessor, but may not claim a direct target-paper nonlinear saturation
reproduction.

## Section 5.4 shock

The Q011 production-science deck exactly carries the Sun-and-Bai textual setup
values for domain, AMR ladder, Mach number, adiabatic index, injection
efficiency, injection momentum, artificial particle light speed, ideal shock
surface, gas subtraction, and removal of particles injected before
`45 Omega0^-1`.

Sun and Bai Section 5.4 does not explicitly state upstream `P0`. The selected
`p0=1` baseline is instead anchored to Bai et al. (2015), PDF page 8, which
explicitly sets `P0=T0=1` and states that this choice is irrelevant while
thermal pressure is much smaller than ram pressure. This provenance distinction
must remain explicit.

Bai et al. dynamically identifies the shock and uses a more involved injection
and compensation treatment; Sun and Bai and Q011 use an ideal injection
surface. Q011 matches Sun and Bai's `eta=1e-3` and smaller domain, not Bai's
fiducial R2 setup. Q011 is also no-Hall.

## Open comparison inputs

No numeric comparison tolerance is frozen by this successor. Approximate text,
undigitized figures, and existing Q011 engineering gates are not promoted into
Q022 reference-derived tolerances. Figure extraction with uncertainty,
branch-specific nonlinear Bell execution contracts, tolerance freeze, and
external scientific review remain open and fail closed.
