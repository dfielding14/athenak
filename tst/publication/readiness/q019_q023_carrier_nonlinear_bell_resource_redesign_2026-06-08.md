# Q-019 Q023-Carrier Nonlinear Bell Resource Redesign

Date: 2026-06-08

Status: **versioned source-local redesign; execution and scientific claims
remain prohibited**

The registered finite-rigidity physical-pilot design is computationally
infeasible under its 500 node-hour cap. Its exact global particle timestep
requires 5.7 million to 782 million cycles per case, with 68.8 million cycles
in Stage 1 alone. This is a preregistered resource-overrun condition, so the
matrix is superseded for execution rather than silently reduced.

The replacement uses the corrected Q023 carrier:

- `rho0=B0=U_A=lambda0=P0=1`, `k0=2 pi`;
- `epsilon=0.4`, `v_CR=2.5`;
- `Omega0=(q/mc)B0=10^-6 k0`;
- `J_CR/c=2 B0 k0`, imposed with the exact root-cell-volume-aware macro weight;
- `C=2500`;
- one centered cold-beam macro-particle per root cell.

The current closure forces `rho_CR/rho0=8e5` and `k0 r_g0=2.5e6`. This is an
explicitly artificial large-inertia current carrier, not a physical
finite-density CR population. It is also not an exactly locked external
current. Any eventual publication statement must therefore be limited to an
ideal-MHD Bell calculation with a self-consistently advanced current carrier
whose lab-frame current, momentum, and energy are measured to remain
invariant.

The pilot matrix contains:

1. Two short 2D linear/onset seeds through `tau=12`.
2. Two 2D nonlinear-window seeds through `tau=30`.
3. Paired resolution, timestep, Riemann-solver, reconstruction, and
   `q/m={10^-4,10^-5,10^-6} k0` controls at fixed current.
4. Short 3D fiducial and large-box cost/mode pilots.

Every pilot is bounded by a design safety envelope `B/B0=10`, which is used
only for timestep and resource planning. It is not a saturation threshold.
The unmocked estimator gives fewer than 20,000 cycles for every row and keeps
the largest pilot at 4,194,304 root cells. Actual node counts, seconds per
cycle, memory, terminal time, saturation window, and acceptance tolerances
must be frozen from excluded Frontier measurements before qualifying output.

Required prerequisites remain the passed registered Q043 deposited-current
matrix followed by the passed registered Q023 linear matrix. No launch,
qualification, fixed-current-like classification, saturation claim, or
publication authority is created by this record.

## Native host validation

The source-local redesign was compiled on 2026-06-08 with the pinned Kokkos
revision `08ceff92bcf3a828844480bc1e6137eb74028517`, Release mode, MPI off, and
OpenMP off. The resulting executable SHA-256 was
`db4d58b979789fc8ee7280eacf2978ee95c9d76bc9a33851f57c58ea543198ba`.

The first exact Stage-1 deck initially failed closed because Athena inserts
the runtime default `particles/pic_boundary_conservation_ledger=0` before the
problem generator computes its immutable semantic identity. The Python
canonicalizer did not yet model that field. The modeled defaults, all 91
Q019 deck identities, the compiled case registry, historical manifests, and
runtime-controller overlays were regenerated together. No checksum bypass or
runtime exception was introduced.

The corrected exact deck `q019-q023-carrier-s1-onset-s0` then completed on one
serial host rank through `t=1.909859317102744` in 445 cycles, compared with the
3,082-cycle conservative planning estimate. It advanced 8,192 active cells
and 8,192 particles in 7.31 s wall time with 54,120 kB peak resident memory.
This is a source-local functional and resource sanity check only. It is not a
Frontier benchmark, registered evidence, or execution authorization.
