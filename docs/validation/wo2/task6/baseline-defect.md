# Periodic LF conservation defect isolated during Task6

The corrected-compiler, immutable baseline loses total energy and magnetic moment
in a periodic, frozen-flow LF sweep on a smooth refined mesh. This defect is
independent of the physical-boundary refill fix and the Task3 communication
optimization. No production source was changed for this diagnosis.

The probe reuses `divb_amr` with zero velocity, periodic boundaries, smooth density
and divergence-free magnetic field, and two touching refined regions. Both
conserved and primitive prolongation were tested to the fixed time 0.002.
Full-precision restart integrals confirm the history result: on the 32-cell root
mesh, conserved prolongation changes energy by -4.445033718880609e-8 and magnetic
moment by -3.0149784668864754e-8. Primitive prolongation changes them by
-1.0930872873515796e-8 and +2.8896616210971615e-9. The unmodified old baseline and
refill-only candidate have identical normalized restart states. All LF floor,
nonfinite, positivity, instability/hard-wall counters and all twelve AMR repair
counters remain zero.

A stage ledger isolates the first change to `STSUpdateU`, not A/mu conversion,
restriction, prolongation, or magnetic evolution. No active face field changes
in the traced first cycle. For conserved prolongation, the first post-receive
global weighted flux divergence is 1.023527494183719e-9 in energy and
1.1131340939532e-10 in magnetic moment; the corresponding update changes the
integrals by their negatives to rounding.

Pairing every physical block face locates the residual entirely at same-level
faces adjoining the refinement corner. Coarse/fine corrected fluxes cancel to
about 2e-23. The two significant unmatched shared faces are:

- x=.25, y in [31/64,32/64], between block IDs 5 and 8;
- y=.5, x in [16/64,17/64], between block IDs 8 and 15.

The initial hypothesis was independently prolonged transverse corner ghosts.
Direct pre-flux state inspection confirms it and narrows the mechanism to the
magnetic reconstruction. Physical fine cell (15,32) has identical ghost density
in blocks 5/8/15 but different Bcc and parallel pressure (1.002882896,
1.005042352, and 1.005023450). Block 15's ghost copy of active block 5 cell (15,31)
also has different B_y (.221240235 versus .222438189), producing p_parallel
1.001193254 rather than 1.0. These cells enter the first transverse LF stencil.
The ordinary cell-centered flux correction synchronizes coarse/fine interfaces
only, leaving these same-level flux estimates independent.

The proposed separate correctness fix synchronizes IEN/IAN face fluxes between
same-level neighbors during multilevel LF sweeps, retaining the existing
coarse/fine correction and all boundary/projection work. It caches both original
estimates and forms an identical overflow-safe symmetric mean. Equal inputs are
preserved exactly. Operand order is fixed on both sides even under FMA
contraction. The existing mesh startup guard rejects shearing boxes with
refinement, so no remapped shear face enters this path. Generic MPI send/receive
completion already covers these request slots.

At the diagnosis checkpoint, CPU/HIP scratch binaries compiled and linked, and
runtime acceptance was still pending. Subsequent validation is recorded in
`TASK6_REPORT.md`: standalone and combined CPU/HIP one/four-rank checks, causal
traces, and smooth-resolution audits passed. The permanent regression
deliberately keeps a 5e-12 absolute conservation tolerance for
both energy and magnetic moment, requires every repair counter to stay zero,
and checks frozen active fields and one/four-rank identity. It correctly rejects
the old baseline; the tolerance has not been weakened.

Evidence under this directory:

- `runs/task6-ledger-hip/trace-{conserved,primitive}/stage-ledger.json`
- `runs/task6-ledger-hip/trace-{conserved,primitive}/face-balance.json`
- `runs/task6-ledger-hip/trace-conserved/shared-face-stencil.json`
- `runs/task6-smooth-kinematic-hip/convergence.json`
- `task6-sync.after-refill.patch` and `task6-sync.after-task3.patch`
- `task6-sync.regression.patch` and `task6_regression_candidate/`

Task3's earlier bitwise proof against the refill-only reference remains intact;
the new synchronization is a separate change that intentionally corrects the
baseline's flux mismatch.
