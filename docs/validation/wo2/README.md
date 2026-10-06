# WO2 CGL-LF implementation and Frontier validation

The implementation branch is `c/cgl-lf-wo2`, based on
`9a4030b09307bc43865d5e597638f8645a6388f8` (including the WO1 merge
`ac33b6b04b8093c2d65e87a30db93148fa65714d`). The original working tree and its
uncommitted engineering-guide change are preserved. All builds, tests, traces,
inputs, raw outputs, and scheduler logs are retained under
`/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2`.

Conservative A and both strict B4 expected failures remain unchanged. The
unresolved sharp-contact limitation is not treated as a passing test.

The [Task 1 report](task1/README.md) records the new timestep bound, independent
proof checks, changed expectations and paired GPU measurements.

The [Task 3 report](task3/README.md) records compact LF cell/flux messages and
CPU/HIP one/four-rank bitwise state preservation. With one scalar the selected
message payload shrinks by 71.43%; paired GPU timings show no consistent
throughput gain. Existing once-per-sweep conversion and primitive refresh remain.

## Task 0: establish the actual GPU baseline

Triage: **implemented**, including diagnosis of two compiler-contract failures
before attributing any numerical difference to WO2.

The reference uses double precision, Release, MPI, AMD MI250X (`gfx90a`),
ROCm 6.4.2, CCE 20.0.0, Cray MPICH 9.0.1, driver 6.16.13, and the committed
Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517`. Each HIP rank receives one GPU
and seven CPU cores with `--gpu-bind=closest --cpu-bind=threads`. The CPU
reference uses the same source revision and MPI toolchain.

Two independent build issues initially produced genuine baseline failures:

1. CCE's startup FTZ/DAZ mode flushed a representable subnormal intermediate
   in the existing extreme-coefficient timestep regression. Relinking the
   unchanged CPU objects with `-mno-daz-ftz` restored the reference timestep
   `0.0097915166977773482`, previously replaced by the `1e30` sentinel.
2. CCE 20's HIP optimizer miscompiled a Kokkos minimum reduction. On the
   ordinary density-contrast-200 input, device per-cell calculations gave the
   correct stiffness, but the reduction discarded its minimum. Neither
   `-ffp-model=precise` nor `strict` repaired this. `-fno-cray` did, matching
   the [upstream Kokkos issue](https://github.com/kokkos/kokkos/issues/9476)
   and [documented Cray/HIP workaround](https://kokkos.org/kokkos-core-wiki/known-issues.html#cray-clang-with-hip).
   The initial cycle timestep became `3.8259896673182957e-5`; the seeded
   temperature deviation after three cycles fell from `1.092832` to
   `8.898016e-7`, with the original `1e-6` seed and scientific bounds intact.

The corrected HIP reference compiles all application and Kokkos objects with
`-fno-cray -mno-daz-ftz` and links with `-no-pie`. Inherited CCE 18 include-path
overrides are cleared when selecting CCE 20. No reduction or physics formula
was changed to compensate for the compiler.

Frozen corrected binaries:

| Backend | Retained executable | SHA-256 |
| --- | --- | --- |
| CPU | `baseline/bin/athena-cpu-noftz` | `f2a05e83b0404bc4f49788b114dc27408f23e47dcc3dcfe8afd87b96b0189f09` |
| HIP | `baseline/bin/athena-hip-nocray` | `adcdadd805c63fd774ef69bd7673c9b0ab52e77cf37e2e3913cff2b890c18150` |

Completed baseline coverage:

- CPU: 213 passed and two accepted strict B4 expected failures. This is the
  union of the full invocation and documented harness reruns, recorded by
  `baseline/cpu-noftz-aggregate.json`. Two reruns repair missing analysis-path
  symlinks. Nineteen Stage I utility fixtures use an exact-source import hook
  that redirects only their canonical-root constant to an unused test path;
  otherwise their placement beneath the real production root triggers the
  utility's lock guard. Production source, safeguards, and assertions are
  unchanged. The first failures and rerun provenance are retained.
- HIP: all 29 full-workflow cases, both active paper smoke cases, and all
  11 AMR GPU checks pass. The broad CPU-named suite also passes with the HIP
  application launcher: 194 passed and two accepted B4 expected failures in
  2265.14 seconds. Tests which build their own Serial helper/application remain
  CPU coverage; the suite name alone is not GPU evidence. Passive mode remains
  explicitly fenced.
- CPU MPI: seven passed. HIP MPI: seven passed, including collisionless shear,
  oblique decay, timestep refresh, and forcing/restart coverage.
- HIP MPI AMR: four-rank restart/conservation and two-node, 16-rank 3D churn
  pass. The two-node run is job `5628668` on `frontier[06610,06651]`.
- Fixed WO1 outputs: all 12 inputs, 18 selected rank/case combinations,
  one warm-up plus three repeats each on job `5628672`, `frontier05888`.
  `matrix/baseline-nocray/results.json` records zero physical repeatability
  failures, raw output hashes, stage counts, solver timing, and launch-inclusive
  timing. These retained inputs enable profiling; they are not uninstrumented
  performance results. Nested profile regions must not be added together.

The application writes nine uninitialized, unused coarse-index integers in the
root `RegionIndcs` restart header. Raw restart hashes therefore vary between
identical baseline runs. `p1-research/restart_compare.py` validates the header
layout and excludes exactly those 36 bytes for state comparisons, retaining
each raw hash and byte range. It never rewrites an output or excludes state,
time, cycle, diagnostics, or forcing data. All binary fields and history files
repeat byte-for-byte. This pre-existing serialization issue is not silently
described as raw restart byte identity.

The original-grid custom-pgen active turbulence input (48 x 48 x 96 cells,
24 x 24 x 48 blocks, unchanged forcing seed/physics/tlim) also passes. One
64-cycle warm-up and three unprofiled repeats give median solver wall time
0.18398578125 seconds/cycle and nine stages per half-sweep. A separate profile
run retains synchronized kernel timing. The full original tlim=1 run completes
620 cycles in 107.411 solver seconds with zero LF repair/wall counters and
maximum relative mass drift 1.1102230246251565e-16. This is original-duration
stability/throughput evidence, not a turbulence resolution-convergence study.
A later full-precision audit established that this original input is exactly
planar and its LF q-work is zero. It therefore does not validate nonzero heat
transport. Task1 includes a separately labeled 3D forcing pair with verified
nonzero LF transport; see the linked report.

The turbulence harness initially failed its profiling command because the input
had not declared the two profiling parameters. The corrected input declares their
unchanged false defaults; only the separate profiling run turns profiling on.
Both original failures and the corrected run are retained. All histories, binary
fields, and final restart state repeat exactly. Initial turbulence restarts also
serialize dormant, uninitialized RNG fields. An explicit opt-in comparison validates
native forcing metadata version 3, 24 metadata/config integers, time/cycle/update
count zero, idum=-1, and iset=0, then excludes only 272 bytes of unused shuffle
state plus the unused 8-byte Gaussian cache. Seed/cache flags, all live forcing state, and physical fields remain compared.
Task1 adds a separate opt-in for the four bytes of native RNG alignment padding;
its ABI proof and exact ranges are recorded in that report. Source inspection proves
those dormant fields are initialized before their first use; mutation probes verify
that changes outside the two ranges still fail comparison. No output is rewritten.
The default comparison still excludes only the nine root coarse-index integers.

Evidence is retained in `turbulence-research/runs/baseline-1rank-recheck`, the
separate `startup-rng-optin-audit.json` reports and their summary, plus the
custom-binary build manifest. The custom baseline executable SHA-256 is
`2904ce1fe800d7b82c32f9163eade045f4e75773582eb2bf567ab57c4453f0ae`.
The complete baseline gate precedes the production changes below; speedups are
reported only after paired candidate comparisons.

## Task 0-P1: magnetic boundary dependency and a separate baseline defect

Triage: **implemented diagnosis and separate correctness fix**. The compact
communication optimization is evaluated separately in Task 3.

With the original outflow/SMR input, compacting only cell-centered IEN/IAN
exchange preserves every physical output byte on CPU and HIP, on one and four
ranks. All 110 one-rank and 440 four-rank fine/coarse state snapshots agree.
Skipping magnetic physical BCs alone reproduces the rejected WO1 drift:
maximum active energy difference `3.45158577e-3`, Bx `1.76364183e-3`, and By
`1.99840963e-3` at cycle 4. Total-energy history first changes at cycle 1.

The first difference is in face-B corner ghosts of blocks 2 and 4 where
physical and coarse/fine boundaries meet. Removing the magnetic fill leaves
zeros there. Primitive refresh then changes cell B and pressures; the next
LF stage reads them and changes active energy and magnetic moment. Later
hyperbolic evolution propagates the error into active B. Magnetic BCs must
remain, including their shearing-boundary dependencies.

The trace also exposed a separate defect in the existing baseline:
`ApplyPhysicalBCs -> Prolongate -> ConToPrim` can leave physical-corner ghosts
stale after prolongation changes their transverse donors. The very first LF
flux reads these ghosts. The fix retains the original fill for coarse-stencil
construction and fills the physical boundaries again after prolongation,
before conversion or the next LF stencil. Supported copy/reflection BCs
preserve whichever A or magnetic-moment representation is currently active.

Independent analytic validation uses a uniform state with nonzero oblique B,
velocity, anisotropy, and a scalar on the same 25-block refinement geometry.
The baseline develops a perpendicular-pressure error `5.7459e-5` and velocity
error `7.7263e-6` after four cycles. The fix preserves the state on CPU/HIP,
one/four ranks, through all active cells and the complete one-cell LF halo.
Double-precision trace residuals stay below `5e-13`; the permanent output-based
regression allows `2e-7` solely because binary field output is float32. Its
eight CPU/GPU/serial/MPI cases at initialization and cycle 4 all pass.
The baseline fails the same checker, so this is not a comparison that merely
reproduces the implementation.

The original nonuniform input changes by design after this correctness fix:
maximum active differences at cycle 4 are energy `2.49624e-4`, parallel
pressure `1.59323e-4`, perpendicular pressure `1.74940e-4`, Bx `2.80142e-6`,
and density `1.72853e-6`. CPU/HIP and one/four ranks agree to field-output
precision. Time and cycle are unchanged. HIP total-energy history changes by
`3.8769e-7`; the maximum normalized divergence remains approximately
`1.7e-15`. This outflow problem does not require constant box energy.

Evidence is in `p1-research/runs/{cpu,hip}-corrected`,
`p1-research/runs/refill-{cpu,hip}`, `original-refill-quantification.json`,
the analytic trace/checker reports, and the permanent regression JUnit files
under `baseline/task0-refill-{cpu,hip}-cpu`. No tolerance or collision setting
was changed to obtain these results.

## Additional baseline conservation gap for Task 6


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

CPU/HIP scratch binaries compile and link. Runtime conservation, MPI, GPU, and
smooth-convergence acceptance of this fix remain pending. The new permanent
regression draft deliberately keeps a 5e-12 absolute conservation tolerance for
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

## Task 2: optional merged half-sweeps

Triage: **adapted and implemented**, with a narrow eligibility gate that preserves
noncommuting collision schedules and strict failure semantics. The mode is off
by default and falls back outside uniform periodic collisionless active-CGL LF.
The [Task2 report](task2/README.md) records CPU/HIP/MPI default/fallback identity,
smooth second-order agreement, synchronization/rollback checks and paired GPU
timings including snapshot overhead. It does not claim a speedup for the original
paper turbulence deck, whose limiter physics makes it ineligible.
