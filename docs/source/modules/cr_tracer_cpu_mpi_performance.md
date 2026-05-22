---
orphan: true
---

# CR Tracer CPU/MPI Performance Comparison

This page records the local CPU and MPI performance comparison for the CR tracer
particle hardening work on `feature/CR_tracers` and the first follow-up
architecture slice.  It is intentionally limited to the hardware available on
2026-05-22: CPU execution with the Kokkos Serial backend and Open MPI.  GPU
timing is still required before making claims about accelerator performance.

## Scope

This page is a CPU/MPI performance record, not a general performance claim.  All
numeric results below were measured on this local machine with the Kokkos Serial
backend and Open MPI.  GPU timing, GPU-aware MPI behavior, and accelerator
copy-volume profiling are deliberately left as follow-up work.

The reported comparison runs all use an MPI-enabled executable.  The
follow-up matrix uses `mpiexec -n 1` even for the one-rank rows so the command
path is consistent across rank counts.  No final timing claim in this page
comes from a separate non-MPI CPU build.

The comparison is split into two layers:

1. The first CR tracer feature branch, `feature/CR_tracers`.
2. The follow-up architecture slice that sits on top of that branch.

The first layer tests the branch-local performance fixes from
`feature/CR_tracers`:

- replacing the host-side `O(N_particles * N_meshblocks)` AMR destination scan
  with mesh-tree lookup,
- reducing `Particles::RemapAfterAMR()` host mirrors from all particle real and
  integer fields to only `IPX`, `IPY`, `IPZ`, and `PGID`,
- preserving MPI particle exchange after AMR remap,
- making reduced `df`, `pspec`, `pspec2`, `dxh`, `drh`, `dparh`, and `pmom`
  diagnostics the MPI default while keeping per-rank diagnostic files as an
  option,
- splitting consistency checks into `none`, `counts`, `local`, and `full`.

The most important result from the first layer is that the mesh-tree AMR remap
is a large CPU/MPI win when remapping is exposed directly, and it remains a
measurable end-to-end win in a real Boris+MHD AMR run.  The later position-only
host copy is neutral within noise for the small CPU Boris benchmark; it is
still the right implementation because it removes unnecessary data movement
from the AMR remap path.

The second layer tests the first follow-up architecture slice:

- `amr_remap = device_table`, a flattened finest-level logical AMR lookup table
  used from a Kokkos remap kernel,
- `exchange_mode = alltoall_counts`, which exchanges rank particle counts
  without all-gathering detailed send tuples,
- reusable particle MPI send/receive buffer capacity,
- prefix-sum particle compaction after exchange.

For the second layer, the CPU/MPI matrix compares four runtime modes:

| Mode | AMR lookup | MPI exchange metadata |
|------|------------|-----------------------|
| Default | `device_table` | `alltoall_counts` |
| Host-tree only | `host_tree` | `alltoall_counts` |
| Allgather only | `device_table` | `allgather` |
| Both fallback | `host_tree` | `allgather` |

That matrix isolates the table lookup and count-exchange switches from one
another.  The reusable buffers and prefix-sum compaction are not runtime
switches, so this page treats them as part of the current implementation rather
than claiming an isolated speedup for each one.

The comparison does not measure:

- GPU execution,
- multi-node MPI,
- OpenMP, pthreads, CUDA, HIP, or SYCL Kokkos backends,
- GPU-aware MPI,
- device memory traffic,
- a pre-buffer-reuse or pre-prefix-compaction implementation with all other
  follow-up changes held fixed.

The first benchmark block should therefore be read as the CPU/MPI evidence for
the first CR tracer feature branch, and as the baseline that the follow-up
architecture branch must beat.  The final benchmark block should be read as the
CPU/MPI evidence that the first follow-up architecture slice is at least neutral
and sometimes modestly faster on this machine while preserving the same AMR
activity and particle correctness.

## How To Read This Comparison

The performance fixes in this branch do not all have the same kind of
evidence.  Some were measured against a direct old-code comparator, some are
runtime-selectable and can be compared in one executable, and some are
structural changes that are always enabled in the current implementation.  The
tables below keep those categories separate so the documentation does not
overstate what was measured.

| Evidence type | Fixes covered | How the comparison was made | Strong claim allowed |
|---------------|---------------|------------------------------|----------------------|
| Direct old-code comparator | Old global AMR scan and full particle host mirror. | A temporary worktree restored the old scan/full-copy path while keeping the same pgen, inputs, compiler, MPI build, and AMR pattern. | The branch is faster than the old scan/full-copy path for these CPU/MPI inputs. |
| Intermediate branch comparator | Position-only AMR host mirror after the tree lookup was already present. | Current branch numbers were compared with commit `b12e02c1`, which used tree lookup but still copied full particle arrays. | The copy reduction is neutral in the small CPU Boris benchmark and removes unnecessary remap data movement. |
| Runtime switch matrix | `device_table` versus `host_tree`; `alltoall_counts` versus `allgather`. | One executable was run with four override combinations for the same inputs and rank counts. | The default follow-up choices are mostly neutral to faster locally; the four-rank remap-focused case shows the smallest table-lookup gain. |
| Semantics comparison | Reduced diagnostics versus per-rank diagnostic files. | The same four-rank stress input was run with reduced output and rank-local output. | Reduced output improves usability and file count without a measured local runtime penalty. |
| Always-on structural changes | Reusable exchange buffers and prefix-sum compaction. | These are included in every current run and were not isolated behind runtime switches. | They preserve CPU/MPI behavior while replacing allocation and host-serial structures with more scalable ones. |

For CPU/MPI timing, the safest metric for code-path changes is AthenaK
`cpu time used`.  Wall time is also reported because short MPI tests are what
users actually run, but wall time includes process launch, filesystem setup,
and local scheduling noise.  The documentation therefore treats a speedup as
strong only when it is visible in AthenaK CPU time and the AMR activity is the
same.

For GPU performance, this page makes no claim.  Every GPU-relevant conclusion
is marked as expected follow-up work and must be repeated on an accelerator
machine with profiler evidence.

## Compared Implementations

Three implementation states are useful for interpreting the first-layer
numbers.

| Name | AMR destination lookup | AMR host mirror | MPI particle exchange | Purpose |
|------|------------------------|-----------------|-----------------------|---------|
| Old scan/full-copy | Host loop over every live MeshBlock for every particle | Full real and integer particle arrays | Existing path | Direct baseline for the expensive original remap behavior. |
| Tree/full-copy | Host mesh-tree lookup by finest logical location | Full real and integer particle arrays | Existing path | Intermediate state used to isolate the lookup algorithm from copy volume. |
| Feature branch current | Host mesh-tree lookup by finest logical location | Positions plus `PGID` only | Existing path | Branch implementation on `feature/CR_tracers`. |

The first feature branch still performed the AMR remap lookup on the host, so
the term "tree" in the first-layer comparison does not mean the follow-up
device-table lookup.  It means that each particle is mapped to the owning leaf through
`MeshBlockTree::FindLeafContaining()` instead of through a global leaf scan.

The follow-up branch adds the `device_table` runtime mode.  With the Serial
backend used here, that table lookup is still a CPU Kokkos kernel.  The name
describes the data structure and execution path needed by device backends; it
does not mean a GPU was used in these measurements.

## Performance-Fix Map

| Fix | Expected CPU/MPI effect | Measured CPU/MPI result | Current conclusion |
|-----|-------------------------|-------------------------|--------------------|
| Mesh-tree AMR lookup | Removes the dominant `N_particles * N_meshblocks` search during AMR remap. | 22.6x-147.4x faster AthenaK CPU time in the remap-focused drift AMR benchmark, depending on rank count. | Keep as the default; this is the main branch-local performance fix. |
| Position-only AMR host mirror | Reduces remap data movement and avoids copying velocity, mass, field sample, displacement, tag, and species data when only position and destination `PGID` are needed. | Neutral within noise in the small CPU Boris benchmark once tree lookup is already present. | Keep it anyway.  It simplifies the remap data path and is expected to matter more on GPU/device backends. |
| Position-based multilevel ownership update | Replaces the immediate-neighbor particle `NewGID` stencil under multilevel meshes with the same position-based lookup strategy used by AMR remap. | Fixed a four-rank AMR out-of-block particle failure; included in the current timing tables. | Keep for robustness.  On the Serial backend, fixed AMR/table overhead becomes visible in the four-rank remap-focused microbenchmark. |
| Existing MPI exchange after remap | Avoids destabilizing correctness by keeping the validated particle migration path. | End-to-end Boris+MHD AMR is 1.37x-1.39x faster than the old scan/full-copy path. | Correct for the feature branch, but this is still a known scaling target for the follow-up branch. |
| Flattened `device_table` AMR lookup | Replaces host tree traversal inside the particle loop with an indexed finest-level lookup. | 1.13x-1.63x relative to `host_tree` in the drift AMR matrix and neutral to 1.02x faster in the Boris+MHD matrix. | Keep as the follow-up default.  The CPU gain is workload-dependent, but this is the right structure for GPU/device execution. |
| `alltoall_counts` exchange metadata | Avoids all-gathering detailed send tuples when only per-rank receive counts are needed. | Neutral within local CPU/MPI noise in the drift, exchange, and Boris+MHD matrices. | Keep as the default because it removes a detailed all-gather scaling hazard without a measured local penalty. |
| Reusable particle MPI buffers | Avoids reallocating send/receive buffers on every exchange when capacity is already sufficient. | Included in the current follow-up runs, but not isolated by a runtime switch. | Keep.  Add a dedicated allocation-count/profile benchmark if this becomes a suspected bottleneck. |
| Prefix-sum survivor compaction | Replaces serial hole-filling after particle exchange with a Kokkos survivor mask and scan. | Included in the current follow-up runs, but not isolated by a runtime switch. | Keep.  The Serial backend may not show much benefit, but the implementation removes a host-serial structure needed for scalable backends. |
| Reduced MPI diagnostics by default | Replaces rank-local histogram files with one globally reduced diagnostic file by default. | No significant local four-rank runtime penalty; file count drops from four files to one file per diagnostic. | Keep as the user-facing default; retain per-rank mode for debugging. |
| Split consistency-check modes | Makes cheap count checks available without forcing full host-side invariant checks. | The small four-rank stress run shows only noise-level overhead for all modes. | Use `counts`, `local`, and `full` in tests/debugging; keep production default at `none`. |

The CPU/MPI evidence is strongest for the mesh-tree lookup, the
`device_table` lookup, and diagnostic-output semantics, because those rows have
direct timing comparisons that isolate the change under study.  The evidence is
intentionally more conservative for reusable buffers and compaction: they are
performance-oriented architecture changes, but the present local benchmark does
not assign them a standalone speedup number.

## Fix-By-Fix CPU/MPI Evidence

This section is the detailed audit of each performance fix.  It separates what
was directly measured from what is supported by structure, correctness tests, or
scaling logic.

### 1. Mesh-Tree AMR Remap Versus Global Scan

The original AMR remap path found a particle's owning MeshBlock by scanning the
global leaf list on the host.  That made the search cost scale like
`N_particles * N_leaf_meshblocks` for every remap.  The feature branch replaces
that search with a logical-location lookup through the AMR tree.

Direct CPU/MPI comparator:

- temporary old-code worktree at the same production commit,
- same Release MPI build, pgen, inputs, compiler, and moving AMR pattern,
- old path restored only for the destination search and full-copy mirror,
- one, two, and four local MPI ranks.

Measured result:

- `147.4x`, `77.9x`, and `22.6x` lower AthenaK CPU time in the remap-focused
  drift AMR case at one, two, and four ranks,
- `1.39x`, `1.37x`, and `1.39x` lower AthenaK CPU time in the end-to-end
  Boris+MHD AMR case at one, two, and four ranks,
- identical MeshBlock creation and deletion counts between old and new paths.

Conclusion:

The mesh-tree lookup is the highest-confidence CPU/MPI performance fix in the
branch.  It removes the pathological global scan without changing the AMR
workload or particle migration semantics.

Remaining limits:

- the local machine only reaches four MPI ranks in this comparison,
- the direct old-code comparator was a temporary timing worktree, not a
  supported runtime mode,
- GPU/device behavior is not measured here.

### 2. Position-Only AMR Host Mirror

After the tree lookup exists, the AMR remap only needs particle positions to
compute destination ownership and `PGID` to store the answer.  The feature
branch therefore stops mirroring velocity, mass, sampled field, displacement,
tag, and species arrays for this remap operation.

Direct CPU/MPI comparator:

- commit `b12e02c1`, which already had tree lookup but still copied full
  particle arrays,
- same small Boris+MHD AMR CPU/MPI benchmark,
- one, two, and four local MPI ranks.

Measured result:

- `0.99x`, `1.01x`, and `1.00x` relative CPU-time ratios at one, two, and four
  ranks,
- no measurable local CPU speedup once the global scan was already removed.

Conclusion:

This is a data-movement cleanup rather than a demonstrated local CPU speedup.
It should stay because it makes the remap contract explicit and removes
unnecessary transfer volume from the path that must eventually run efficiently
on accelerator memory spaces.

Remaining limits:

- this benchmark is too small to expose CPU memory-bandwidth effects,
- no host/device copy-volume measurement is possible on this machine,
- larger particle arrays should be timed before treating the CPU result as
  universal.

### 3. Position-Based Multilevel Ownership After Pushes

The first implementation still used an immediate-neighbor `NewGID` stencil for
post-push particle ownership updates.  That assumption is too weak for
multilevel AMR and MPI migration.  The follow-up branch now uses
position-based ownership lookup on multilevel meshes, matching the AMR remap
logic.

CPU/MPI evidence:

- the four-rank AMR stress run failed before this change with an out-of-block
  particle under `check_consistency_mode = full`,
- the same four-rank stress run now passes with `amr_remap = device_table`,
  `validate_amr_lookup = true`, and 105 MeshBlocks created and 77 deleted,
- the explicit `host_tree`/`allgather` fallback stress also passes.

Measured performance status:

- included in the current benchmark tables,
- not isolated behind a runtime switch,
- fixed AMR/table overhead is visible in the four-rank remap-focused Serial
  backend case, where local particle work is very small.

Conclusion:

This is primarily a robustness fix.  It is performance-relevant because an
incorrect cheap ownership shortcut is not an acceptable optimization for AMR.
The next CPU/MPI optimization target is reducing or caching the table rebuild
cost when many small rank-local remaps occur.

### 4. Flattened `device_table` AMR Lookup

The follow-up branch adds `amr_remap = device_table`, which maps finest-level
logical cells to owning leaf `gid`s and uses an indexed Kokkos remap kernel.
On this machine that still runs on the Kokkos Serial backend, but it exercises
the same data-structure direction needed by GPU backends.

Runtime CPU/MPI comparator:

- same executable,
- `device_table` versus `host_tree`,
- `alltoall_counts` held fixed when isolating lookup,
- one, two, and four local MPI ranks.

Measured result:

- drift AMR matrix: `1.63x`, `1.25x`, and `1.13x` faster than `host_tree` at
  one, two, and four ranks,
- Boris+MHD AMR matrix: neutral to `1.02x` faster than `host_tree`.

Conclusion:

The default should remain `device_table`, but the CPU/MPI result is
workload-dependent.  The table is useful when lookup work is exposed, neutral
in the end-to-end Boris case, and can be slower on this Serial-backend
microbenchmark when rebuild overhead dominates.

Remaining limits:

- no threaded CPU backend was measured,
- no GPU memory-space behavior was measured,
- the present table rebuild strategy should be profiled at larger rank counts
  and particle counts.

### 5. `alltoall_counts` Exchange Metadata

The follow-up exchange path avoids making detailed send tuples globally
visible.  Ranks exchange receive counts directly, then keep the particle payload
exchange in the existing point-to-point path.

Runtime CPU/MPI comparator:

- same executable,
- `alltoall_counts` versus `allgather`,
- `device_table` held fixed when isolating exchange metadata,
- drift AMR, exchange-focused subcycled drift, and Boris+MHD AMR cases.

Measured result:

- neutral within local CPU/MPI noise in the drift matrix,
- neutral within local CPU/MPI noise in the exchange-focused case,
- neutral within local CPU/MPI noise in the Boris+MHD matrix.

Conclusion:

This is a scaling-structure fix, not a locally measured speedup.  It should
remain the default because it removes a detailed all-gather pattern without
adding a measurable four-rank CPU penalty.

Remaining limits:

- the local test only reaches four ranks,
- a cluster-scale CPU/MPI run is needed to quantify the expected metadata
  benefit,
- MPI profiling would be needed to split metadata exchange from payload
  exchange and waits.

### 6. Reusable MPI Buffers

Particle exchange now keeps send and receive buffer capacity and grows buffers
only when the required exchange exceeds current capacity.

CPU/MPI evidence:

- all current correctness and performance runs include this path,
- no behavior changes were seen in CPU and MPI regression tests,
- no standalone timing is claimed.

Conclusion:

This is an allocation-control and scaling hygiene fix.  It should stay, but the
current documentation should not assign it a speedup.  A stronger future CPU
comparison would restore per-exchange allocation in a temporary build and
record allocation counts, peak capacity, and exchange timing.

### 7. Prefix-Sum Particle Compaction

The follow-up branch replaces host-serial survivor packing after exchange with
a Kokkos mask and prefix-sum compaction pattern.

CPU/MPI evidence:

- all current correctness and performance runs include this path,
- CPU and MPI stress tests pass with AMR migration, subcycling, and restarts,
- no runtime switch isolates this change.

Conclusion:

This removes a host-serial structure and aligns the CPU path with the execution
model needed by threaded and GPU backends.  On the Kokkos Serial backend, the
honest local claim is behavioral preservation and better architecture, not a
measured standalone speedup.

### 8. Reduced Diagnostics By Default

The branch changes `df`, `pspec`, `pspec2`, `dxh`, `drh`, `dparh`, and `pmom`
diagnostics so MPI runs write globally reduced files by default, while retaining
per-rank files for debugging.

Runtime CPU/MPI comparator:

- four-rank stress input,
- reduced output versus per-rank output,
- three repeated runs.

Measured result:

- reduced default: `0.03779 s` CPU, `0.12066 s` wall,
- per-rank files: `0.03850 s` CPU, `0.12424 s` wall,
- file count changes from four files per diagnostic to one file per diagnostic
  at four ranks.

Conclusion:

The reduced default is a usability improvement with no measured local
performance penalty.  It also makes histogram totals directly comparable to
global species counts, which simplifies automated validation.

Remaining limits:

- the local filesystem and small output size are not a parallel-I/O stress
  test,
- large bin counts, many species, and high rank counts should be remeasured.

### 9. Split Consistency-Check Modes

The branch separates particle consistency checks into `none`, `counts`,
`local`, and `full`.  This allows tests to choose the invariant strength they
need without forcing full host-side validation in production runs.

Runtime CPU/MPI comparator:

- four-rank stress input,
- diagnostics disabled,
- mode sweep over `none`, `counts`, `local`, and `full`.

Measured result:

- `counts` is `1.01x` the CPU time of `none`,
- `local` is `1.03x` the CPU time of `none`,
- `full` is `1.02x` the CPU time of `none`,
- all differences are close to local noise for this small stress test.

Conclusion:

`counts` is cheap enough for frequent regression use.  `local` and `full`
should remain debug/test modes because their host copies and global tag
gathering should scale with particle count and rank count.

## Measurement Method

### Test Environment

| Item | Value |
|------|-------|
| Date | 2026-05-22 |
| Machine | Apple M4 Max, 16 physical/logical CPU cores |
| OS | macOS Darwin 25.5.0, arm64 |
| Compiler | Apple clang 21.0.0 |
| MPI | Open MPI 5.0.9 |
| Kokkos backend | Serial |
| Build type | `Release` |
| Production commit | `ac909458` |
| Tree-remap comparison commit | `b12e02c1` |
| Follow-up matrix branch | `feature/CR_tracers_followup_architecture` |
| Current full-harness sweep | 240 runs, three repeats, one/two/four local MPI ranks |

Build command:

```bash
cmake -S . -B build-cr-perf-mpi \
  -D Athena_ENABLE_MPI=ON \
  -D PROBLEM=part_random \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build-cr-perf-mpi -j 8
```

The direct old-scan comparator was a temporary worktree at `ac909458` with
`Particles::RemapAfterAMR()` changed back to the old global MeshBlock scan and
full particle-array host mirror.  That keeps the pgen, AMR stress pattern, MPI
exchange code, compiler, and input files identical while isolating the remap
lookup/copy change.

That temporary comparator is not a supported runtime mode.  It was used only to
answer the historical question "how much did the old scan cost?"  The supported
fallback for users is `particles/amr_remap=host_tree`, which keeps the tree
lookup and only bypasses the newer flattened lookup table.

The follow-up architecture matrix was measured on the current follow-up branch
with runtime overrides.  The four modes differ only in `<particles>/amr_remap`
and `<particles>/exchange_mode`; all four use the current reusable buffer and
prefix-sum compaction implementation.  This is why the matrix can isolate the
lookup and exchange switches, but cannot isolate buffer reuse or compaction.

Each table reports the median of three runs.  `cpu time used` and
`particle-updates/cpu_second` are AthenaK's own runtime counters; `wall` is the
outer Python subprocess wall time.  The absolute numbers are local benchmark
data, not portable performance guarantees.

The benchmark design keeps the global problem fixed while changing the number
of MPI ranks.  These are therefore useful local strong-scaling comparisons, but
they are not clean microarchitectural profiles.  At these small problem sizes,
process launch, MPI synchronization, filesystem metadata, and AMR bookkeeping
are visible in wall time.  AthenaK CPU time is usually the better metric for
the branch-local code change, while wall time is the better reminder of what a
short test actually costs to run.

The tables below intentionally include both metrics:

- `CPU s`: AthenaK's own measured CPU time.
- `Wall s`: Python subprocess elapsed time including launch/setup overhead.
- `CPU speedup`: old or comparison CPU time divided by current CPU time.
- `Wall speedup`: old or comparison wall time divided by current wall time.

The medians were used because the runs are short enough that a single slow
launch or filesystem event can skew a mean.

For the follow-up matrix, even the one-rank entries were launched with
`mpiexec -n 1` so that every entry exercised the MPI build and the same command
path.  Wall times for the shortest drift cases should therefore be interpreted
with caution: MPI process launch and filesystem setup can dominate a run whose
AthenaK CPU time is only a few milliseconds.

### Benchmark Inputs

| Input | Purpose |
|-------|---------|
| `inputs/particles/cr_tracer_drift_amr_perf.athinput` | Remap-focused AMR benchmark.  It uses drift particles with zero velocity, 64^3 cells, 8^3-cell MeshBlocks, `ppc = 0.125`, two species, moving AMR boxes, and no diagnostic outputs.  This exposes AMR remap and load-balance particle migration without MHD/Boris work. |
| `inputs/particles/cr_tracer_boris_amr_perf.athinput` | End-to-end Boris+MHD AMR benchmark.  It uses 32^3 cells, 8^3-cell MeshBlocks, `ppc = 1`, two species, moving AMR boxes, and no diagnostic outputs. |
| `inputs/particles/cr_tracer_boris_amr_stress.athinput` | Validated debug/output benchmark input.  It is smaller, includes the CR diagnostics, and is the target used for consistency-check and output-semantics timing. |

Representative commands:

```bash
./build-cr-perf-mpi/src/athena \
  -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
  -d run-cr-drift-perf

mpiexec -n 4 ./build-cr-perf-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_perf.athinput \
  -d run-cr-boris-perf-n4
```

## AMR Remap Microbenchmark

This benchmark isolates the AMR remap path by using zero-velocity drift
particles.  All runs produced the same AMR activity: 1078 MeshBlocks created,
742 deleted, and 848 live MeshBlocks at termination.

| MPI ranks | Old scan/full-copy CPU s | Tree/current CPU s | CPU speedup | Old wall s | Tree wall s | Wall speedup |
|-----------|--------------------------|--------------------|-------------|------------|-------------|--------------|
| 1 | 0.5686 | 0.00386 | 147.4x | 0.6177 | 0.05953 | 10.4x |
| 2 | 0.4210 | 0.00540 | 77.9x | 0.4841 | 0.06874 | 7.0x |
| 4 | 0.2786 | 0.01233 | 22.6x | 0.3503 | 0.08285 | 4.2x |

Interpretation:

- The old scan is dominated by repeatedly testing particle positions against
  the global MeshBlock list.
- The tree path removes that scaling and makes MPI overhead, load balancing,
  and process launch a larger fraction of the small benchmark.
- The CPU-time speedup decreases with rank count because the old scan work is
  divided across ranks while fixed MPI overhead becomes more visible.
- The wall-clock speedup is smaller than the AthenaK CPU-time speedup because
  these short runs include process launch and setup overhead.

Rank-by-rank details:

- At one rank, the lookup change is almost a pure algorithmic comparison:
  `0.5686 s` becomes `0.00386 s`, so the old global scan is more than two
  orders of magnitude slower than the tree lookup.
- At two ranks, the old scan cost is partially divided across ranks, but the
  current tree path exposes MPI and AMR metadata overhead.  The CPU speedup is
  still `77.9x`.
- At four ranks, the current remap CPU time rises to `0.01233 s` because the
  useful local remap work is extremely small and fixed overheads dominate.
  Even in that regime the tree path is still `22.6x` faster by AthenaK CPU
  time and `4.2x` faster by subprocess wall time.

The important scaling point is not that the current path is perfectly flat with
rank count.  It is that the pathological multiplication by global leaf count is
gone.  The remaining cost is now small enough that MPI launch, AMR bookkeeping,
and particle exchange become the next visible bottlenecks.

## Boris+MHD End-To-End AMR Benchmark

This benchmark includes MHD evolution, Boris particle pushing, AMR, particle
migration, and no diagnostic output.  All runs produced the same AMR activity:
588 MeshBlocks created, 476 deleted, and 176 live MeshBlocks at termination.

| MPI ranks | Old scan/full-copy CPU s | Current CPU s | CPU speedup | Old wall s | Current wall s | Wall speedup |
|-----------|--------------------------|---------------|-------------|------------|----------------|--------------|
| 1 | 1.1082 | 0.7974 | 1.39x | 1.1742 | 0.8778 | 1.34x |
| 2 | 0.7437 | 0.5433 | 1.37x | 0.8469 | 0.6505 | 1.30x |
| 4 | 0.4789 | 0.3451 | 1.39x | 0.5828 | 0.4512 | 1.29x |

The end-to-end speedup is smaller than in the remap microbenchmark because MHD,
Boris interpolation, and particle push work remain unchanged.  This is the more
realistic expectation for production CPU runs: the fast remap removes a clear
AMR bottleneck, but it does not change the cost of the physical update.

The current implementation removes these fractions of the old end-to-end CPU
time:

| MPI ranks | CPU time removed | Fraction of old CPU time removed |
|-----------|------------------|----------------------------------|
| 1 | 0.3108 s | 28.0% |
| 2 | 0.2004 s | 26.9% |
| 4 | 0.1338 s | 27.9% |

That is a useful production-scale interpretation.  The remap optimization does
not turn the Boris pusher or MHD solver into a different algorithm, but it does
remove roughly one quarter to one third of the local CPU time in this AMR
workload.  The slight decline in removed fraction with rank count is expected: smaller
per-rank particle and MeshBlock counts reduce the old scan burden while
communication and synchronization take a larger share of each short run.

The current path also preserves the same AMR behavior as the old-scan
comparator.  The identical created/deleted/live MeshBlock counts are important:
the benchmark is not faster because it refined less, derefined less, or skipped
particle migration.

## Position-Only Host-Copy Check

The `b12e02c1` implementation used the mesh-tree lookup but still mirrored the
full particle arrays during AMR remap.  The first current implementation
mirrored only positions and `PGID`.  On the small CPU Boris benchmark this copy
optimization was neutral within run-to-run noise when isolated before the later
four-rank AMR ownership fix:

| MPI ranks | Tree/full-copy CPU s | Current CPU s | Ratio |
|-----------|----------------------|---------------|-------|
| 1 | 0.8220 | 0.8325 | 0.99x |
| 2 | 0.5750 | 0.5692 | 1.01x |
| 4 | 0.3741 | 0.3743 | 1.00x |

This does not mean the copy change is unnecessary.  It means the CPU benchmark
is dominated by other work once the algorithmic scan has been removed.  The
position-only copy still reduces remap data movement, keeps the host path
honest, and is expected to matter more for GPU builds where host/device copies
are much more expensive.  GPU timing is a required follow-up.

The current full branch also includes the position-based multilevel ownership
update after particle pushes.  That later robustness fix is included in the
current end-to-end Boris timings above, so the table in this section should be
read only as the historical isolation of copy volume.

For CPU/MPI, this result should be interpreted narrowly:

- It does show that the copy reduction did not harm the small Boris+MHD AMR
  case.
- It does not prove that copy volume is irrelevant for larger CPU particle
  counts.
- It does not measure host-device transfer cost on accelerator backends.
- It does not replace the follow-up goal of moving remap lookup and compaction
  onto the device.

The practical reason to keep the fix is that `RemapAfterAMR()` only needs
position data to find a destination block and `PGID` to record it.  Copying
unrelated particle fields obscures that contract and makes future device-side
work harder to reason about.

## Consistency-Check Overhead

The consistency-check timing used the validated AMR stress input with four MPI
ranks and diagnostic outputs disabled.  The run is deliberately small, so the
numbers measure smoke-test overhead rather than large-production overhead.

| Mode | CPU s | Wall s | CPU ratio to `none` |
|------|-------|--------|---------------------|
| `none` | 0.03758 | 0.12170 | 1.00x |
| `counts` | 0.03806 | 0.12071 | 1.01x |
| `local` | 0.03863 | 0.12058 | 1.03x |
| `full` | 0.03848 | 0.11991 | 1.02x |

Interpretation:

- On the small stress input, all consistency modes are close to noise-level
  overhead.
- `counts` is cheap enough for frequent regression use.
- `local` and `full` copy particle data to host; `full` also gathers duplicate
  tag keys across MPI ranks.  Their cost should scale with particle count and
  rank count, so they should remain debug/test options rather than production
  defaults.
- High-particle synthetic benchmarks should usually run with
  `check_consistency_mode = none` unless the benchmark is specifically testing
  the checks.

The modes are performance tools as much as correctness tools:

| Mode | Intended use | Scaling concern |
|------|--------------|-----------------|
| `none` | Production and high-particle benchmarks. | No extra particle validation cost. |
| `counts` | Cheap regression guard for total and per-species conservation. | MPI reductions scale with species count and rank count, not directly with particle count. |
| `local` | Debugging invalid ownership, bounds, finite values, and local duplicate tags. | Copies local particle data to host and scans it. |
| `full` | Strongest test mode, including duplicate-tag checks across ranks. | Adds global tag gathering; use in small and medium tests, not routine production. |

For future performance tests, `counts` is the best default when the goal is to
guard conservation without measuring the checker itself.  `local` and `full`
should be timed separately when checker overhead is under study.

## Diagnostic Output Reduction

The output comparison used the current stress input with four MPI ranks and
three repeated runs.  The reduced default and per-rank mode had similar runtime
on this local filesystem, but very different file counts:

| Mode | CPU s | Wall s | `df` files | `pspec` files | `pspec2` files | `dxh` files | `drh` files | `dparh` files | `pmom` files |
|------|-------|--------|------------|---------------|----------------|-------------|-------------|---------------|--------------|
| Reduced default | 0.03779 | 0.12066 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Per-rank files | 0.03850 | 0.12424 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |

The performance difference is not significant at this scale.  The reduced
default is still the better user-facing behavior because each diagnostic has one
global file, histogram sums directly match global species counts, and analysis
scripts do not need to merge rank-local histograms.  Per-rank diagnostic files
remain useful for debugging rank-local particle migration and can be enabled
with:

```text
df_single_file_per_rank = 1
pspec_single_file_per_rank = 1
pspec2_single_file_per_rank = 1
dxhist_single_file_per_rank = 1
drh_single_file_per_rank = 1
dparh_single_file_per_rank = 1
pmom_single_file_per_rank = 1
```

This comparison is deliberately about output semantics and local overhead, not
about peak parallel I/O performance.  The measured filesystem is local and the
run produces small files.  The reduced path is still the better default because
it gives a single physically meaningful global distribution.  Per-rank output
is an inspection mode: it answers where particles live and migrate, not what
the global distribution function is.

Expected scaling behavior:

- Reduced histograms add MPI reductions over histogram bins and species.
- Per-rank histograms avoid those reductions but increase file count linearly
  with rank count.
- On a shared filesystem, many rank-local files can become more expensive than
  the small reduction, especially for frequent diagnostic dumps.
- For very large histograms, many species, or very high rank counts, the
  reduction cost should be remeasured.

The current CPU/MPI result says that reduced diagnostics are not a measurable
local penalty at four ranks.  It does not say that reduction cost is always
zero.

## Bottleneck Summary

The CPU/MPI evidence points to this current cost hierarchy:

| Path | Before fixes | Current state | Remaining CPU/MPI question |
|------|--------------|---------------|----------------------------|
| AMR destination lookup | Dominant in remap-heavy AMR runs because of the global leaf scan. | No longer dominant in the tested CPU/MPI runs.  The follow-up `device_table` path is faster than `host_tree` in the drift matrix and neutral to modestly faster in Boris+MHD. | Recheck at larger particle counts and more ranks. |
| AMR host data movement | Mixed with full particle mirror traffic. | Reduced to positions and `PGID`; not measurable in the small CPU Boris test. | Recheck with larger particle arrays and GPU host/device profiling. |
| Particle MPI exchange metadata | Detailed send information was globally visible. | `alltoall_counts` avoids detailed all-gather metadata and is neutral locally. | Test at higher rank count where all-gather metadata should matter more. |
| Particle MPI buffers | Buffers could be rebuilt around each exchange. | Buffers grow only when capacity is insufficient. | Add allocation-count or profiler evidence if exchange becomes noisy. |
| Particle compaction | Hole handling was host/serial in structure. | Survivor-mask and prefix-sum compaction now use the Kokkos execution path. | Isolate against the old compaction path only if profiling shows compaction cost. |
| Boris gather/push | Unchanged by the first performance branch. | Modular gather keeps the same default `tsc` method while adding `trilinear` coverage. | Profile interpolation choice separately; do not mix with AMR remap timing. |
| Diagnostics | Rank-local files were easier to generate but harder to analyze. | Reduced output default improves usability with no local four-rank penalty measured. | Rebenchmark with larger histograms, species counts, and production filesystems. |

This hierarchy is why the first follow-up work concentrated on AMR lookup,
particle exchange metadata, buffer reuse, and compaction before making deeper
numerical changes.  The first feature branch fixed the largest
correctness-adjacent bottleneck without destabilizing particle migration.  The
follow-up slice then changes the architecture in the places that matter for
scaling, even where the local CPU/MPI speedup is modest because the old global
scan bottleneck was already gone.

## Follow-Up Architecture Matrix

This section measures the current follow-up defaults against independent
fallbacks for the two runtime-switchable performance fixes.  The same Release
MPI build and inputs were used.  Each entry is the median of three runs.

The four tested modes are:

| Mode | Command-line overrides |
|------|------------------------|
| Default | `particles/amr_remap=device_table particles/exchange_mode=alltoall_counts` |
| Host-tree only | `particles/amr_remap=host_tree particles/exchange_mode=alltoall_counts` |
| Allgather only | `particles/amr_remap=device_table particles/exchange_mode=allgather` |
| Both fallback | `particles/amr_remap=host_tree particles/exchange_mode=allgather` |

All four modes include reusable MPI buffers and prefix-sum compaction.  These
tables therefore compare lookup and exchange metadata choices, not buffer or
compaction implementations.

### Remap-Focused Drift AMR

All runs produced the same AMR activity: 1078 MeshBlocks created, 742 deleted,
and 848 live MeshBlocks at termination.

CPU-time medians:

| MPI ranks | Default | Host-tree only | Allgather only | Both fallback |
|-----------|---------|----------------|----------------|---------------|
| 1 | 0.004072 | 0.006629 | 0.004080 | 0.006650 |
| 2 | 0.005806 | 0.007252 | 0.005700 | 0.007295 |
| 4 | 0.009824 | 0.011068 | 0.010006 | 0.010738 |

Wall-time medians:

| MPI ranks | Default | Host-tree only | Allgather only | Both fallback |
|-----------|---------|----------------|----------------|---------------|
| 1 | 0.06514 | 0.06748 | 0.06387 | 0.06624 |
| 2 | 0.07091 | 0.07078 | 0.07219 | 0.07242 |
| 4 | 0.08172 | 0.08399 | 0.08287 | 0.08133 |

CPU-time speedups:

| MPI ranks | Table vs host-tree | `alltoall_counts` vs `allgather` | Default vs both fallback |
|-----------|--------------------|----------------------------------|--------------------------|
| 1 | 1.63x | 1.00x | 1.63x |
| 2 | 1.25x | 0.98x | 1.26x |
| 4 | 1.13x | 1.02x | 1.09x |

Interpretation:

- The table lookup is the visible improvement in the remap-focused drift
  benchmark: it is 1.13x-1.63x faster than the host tree by AthenaK CPU time
  across one, two, and four local MPI ranks.
- At four ranks on the Serial backend, the table benefit is smaller than at one
  and two ranks.  The useful local particle work is very small there, so fixed
  AMR, table rebuild, MPI, and launch overheads become a larger fraction of the
  short run.
- The exchange metadata switch is essentially neutral in this drift case.  The
  rank-to-rank sign changes are expected because the benchmark is dominated by
  AMR remapping and load-balance migration, not by repeated high-volume
  particle boundary exchange.
- The one-rank wall-time table is noisier than the CPU-time table because the
  AthenaK work is only a few milliseconds.  It should not be used to argue that
  a slower CPU path is faster in practice.
- The default mode and both-fallback mode generate the same refinement and
  derefinement pattern, so the speedup is not caused by changing the AMR
  workload.

### Boris+MHD AMR

All runs produced the same AMR activity: 588 MeshBlocks created, 476 deleted,
and 176 live MeshBlocks at termination.

CPU-time medians:

| MPI ranks | Default | Host-tree only | Allgather only | Both fallback |
|-----------|---------|----------------|----------------|---------------|
| 1 | 0.812301 | 0.823076 | 0.812418 | 0.822398 |
| 2 | 0.545783 | 0.557509 | 0.556394 | 0.553186 |
| 4 | 0.345134 | 0.349544 | 0.346115 | 0.350253 |

Wall-time medians:

| MPI ranks | Default | Host-tree only | Allgather only | Both fallback |
|-----------|---------|----------------|----------------|---------------|
| 1 | 0.89522 | 0.90836 | 0.89408 | 0.90756 |
| 2 | 0.65570 | 0.66544 | 0.66384 | 0.66275 |
| 4 | 0.45300 | 0.45701 | 0.45572 | 0.46105 |

CPU-time speedups:

| MPI ranks | Table vs host-tree | `alltoall_counts` vs `allgather` | Default vs both fallback |
|-----------|--------------------|----------------------------------|--------------------------|
| 1 | 1.01x | 1.00x | 1.01x |
| 2 | 1.02x | 1.02x | 1.01x |
| 4 | 1.01x | 1.00x | 1.01x |

Interpretation:

- In the end-to-end Boris+MHD workload, the follow-up default is only modestly
  faster than the fallback modes.  The largest original bottleneck, the global
  MeshBlock scan, has already been removed before this matrix starts.
- The table lookup is neutral to 1.02x faster than `host_tree`; most of the run
  time is MHD evolution, Boris gather/push work, AMR bookkeeping, and particle
  exchange rather than AMR destination lookup alone.
- `alltoall_counts` is neutral relative to `allgather` by AthenaK CPU time in
  this local matrix.  The reason to keep it as the default is not a large
  four-rank speedup here; it is the removal of a detailed all-gather metadata
  path that scales worse with MPI size.
- The default-vs-both-fallback comparison is neutral to 1.02x faster by CPU time.
  That is the appropriate expectation for this first architecture slice on CPU:
  it improves structure and removes scaling hazards without changing the cost
  of the physical update.

### Exchange-Focused Subcycled Drift

The exchange case uses high-velocity drift particles with subcycling enabled so
particles cross MPI ownership boundaries during a single mesh timestep.  It is
not an AMR benchmark, so the `amr_remap` switch is intentionally irrelevant
here.  The useful comparison is `alltoall_counts` against `allgather`.

| MPI ranks | Default CPU s | Host-tree CPU s | Allgather CPU s | Both fallback CPU s | Sent particles | Payload bytes | Messages |
|-----------|---------------|-----------------|-----------------|---------------------|----------------|---------------|----------|
| 2 | 0.003491 | 0.003388 | 0.003416 | 0.003464 | 1272 | 157728 | 334 |
| 4 | 0.003096 | 0.003075 | 0.003217 | 0.003188 | 2558 | 317192 | 672 |

Interpretation:

- The all-to-all-counts path is neutral within local noise at two and four
  ranks.  The measured CPU times are only a few milliseconds, so process
  launch and synchronization noise are larger than the metadata difference.
- The reported bytes and messages are particle payload estimates from
  `LogPerformance()`.  They do not include the old detailed all-gather metadata
  cost, so this table should not be read as a full MPI profiling trace.
- The implementation change is still important for scaling: the new path asks
  each rank only how many particles it will receive from each other rank,
  instead of making detailed send tuples globally visible.
- A higher-rank CPU run on a cluster is the right follow-up for quantifying the
  metadata benefit.  On this local machine, the honest result is that the new
  path preserves behavior and does not add measurable overhead.

### Full Default-Mode Harness Sweep

The current benchmark harness also ran all six benchmark cases in the default
mode.  This table is not a replacement for the more controlled comparisons
above; it is a coverage map showing what the default path costs across the main
particle workflows.

| Case | Ranks | CPU s | Wall s | Particle updates/s | Sent | Bytes | Messages | Extra check |
|------|-------|-------|--------|--------------------|------|-------|----------|-------------|
| `push` | 1 | 0.004257 | 0.06497 | 2.69e7 | 0 | 0 | 0 | no AMR |
| `push` | 2 | 0.005274 | 0.07261 | 2.17e7 | 1604 | 198896 | 28 | no AMR |
| `push` | 4 | 0.002728 | 0.08074 | 4.20e7 | 3227 | 400148 | 134 | no AMR |
| `amr_remap` | 1 | 0.004072 | 0.06514 | 4.83e7 | 0 | 0 | 0 | AMR 1078/742/848 |
| `amr_remap` | 2 | 0.005806 | 0.07091 | 3.39e7 | 15817 | 1961308 | 6 | AMR 1078/742/848 |
| `amr_remap` | 4 | 0.009824 | 0.08172 | 2.00e7 | 44452 | 5512048 | 18 | AMR 1078/742/848 |
| `exchange` | 2 | 0.003491 | 0.06879 | 1.17e6 | 1272 | 157728 | 334 | subcycled |
| `exchange` | 4 | 0.003096 | 0.07769 | 1.32e6 | 2558 | 317192 | 672 | subcycled |
| `histogram` | 1 | 0.086640 | 0.14922 | 4.73e4 | 0 | 0 | 0 | AMR 63/35/36 |
| `histogram` | 2 | 0.053556 | 0.12797 | 7.65e4 | 146 | 18104 | 20 | AMR 63/35/36 |
| `histogram` | 4 | 0.037789 | 0.12066 | 1.08e5 | 602 | 74648 | 72 | AMR 63/35/36 |
| `histogram_per_rank` | 1 | 0.086513 | 0.14962 | 4.73e4 | 0 | 0 | 0 | AMR 63/35/36 |
| `histogram_per_rank` | 2 | 0.053678 | 0.12862 | 7.63e4 | 146 | 18104 | 20 | AMR 63/35/36 |
| `histogram_per_rank` | 4 | 0.038500 | 0.12424 | 1.06e5 | 602 | 74648 | 72 | AMR 63/35/36 |
| `restart` | 1 | 0.008488 | 0.06954 | 3.86e6 | 0 | 0 | 0 | restart total 16384 |
| `restart` | 2 | 0.006924 | 0.08412 | 4.73e6 | 410 | 50840 | 8 | restart total 16384 |
| `restart` | 4 | 0.003014 | 0.08475 | 1.09e7 | 772 | 95728 | 34 | restart total 16384 |
| `boris_amr` | 1 | 0.812301 | 0.89522 | 6.45e5 | 0 | 0 | 0 | AMR 588/476/176 |
| `boris_amr` | 2 | 0.545783 | 0.65570 | 9.61e5 | 52095 | 6459780 | 46 | AMR 588/476/176 |
| `boris_amr` | 4 | 0.345134 | 0.45300 | 1.52e6 | 179330 | 22236920 | 230 | AMR 588/476/176 |

The generated 80-row Markdown table and JSON data are intentionally treated as
benchmark artifacts rather than source documentation.  They should be
regenerated on the target machine with the harness command below, because the
absolute timings are machine, MPI, filesystem, and launch-policy dependent.

### Non-Isolated Follow-Up Fixes

Two follow-up performance fixes are intentionally not assigned direct speedup
numbers in the matrix:

- Reusable MPI buffers are always enabled in the current branch.  The right
  direct comparison would be a temporary build that restores per-exchange
  allocation while keeping all other code unchanged, plus allocation counters
  or profiler evidence.
- Prefix-sum compaction is also always enabled in the current branch.  On the
  Serial backend, the Kokkos scan may not beat a small hand-written host loop in
  every short run.  Its purpose is to remove a host-serial structure and make
  the CPU path match the execution model needed by threaded and GPU backends.

Those fixes are still valuable, but the honest CPU/MPI claim is architectural
neutrality plus better scaling structure, not a separately measured local
speedup.

Additional follow-up features are covered by correctness and usability tests
rather than the timing matrix:

- Particle subcycling is opt-in.  The CPU suite includes a high-velocity drift
  case that fails the motion guard without subcycling, passes with subcycling,
  and fails clearly when `subcycle_max_steps` is too small.  It also checks that
  subcycling reduces uniform-field Boris phase error.
- The MPI CPU suite includes a high-velocity drift case that performs
  intermediate ownership exchange across two ranks during one particle push.
- Particle-aware AMR load balancing is opt-in through
  `<mesh_refinement>/particle_load_weight`; the MPI AMR stress test runs with a
  nonzero weight while conserving particle totals.
- `cr_tracer_inspect.py` now reports particles per MeshBlock `PGID` and restart
  format metadata in addition to rank/species totals, spectrum checks, and
  histogram checks.  New `*.prst.pmeta` sidecars record the particle restart
  payload size and checksum while preserving the legacy binary `prst` layout.

The CPU/MPI acceptance evidence for this slice is:

- four-rank AMR stress with default modes passed and reported 105 MeshBlocks
  created and 77 deleted,
- `cr_tracer_inspect.py` reported 1024 total particles and 512 particles per
  species,
- the explicit `host_tree`/`allgather` fallback passed the same stress and
  inspector checks,
- the CPU particle suite passed with 20 tests,
- the MPI CPU particle suite passed with 4 tests,
- the style suite passed,
- `git diff --check` passed.

## Reproducibility Checklist

To reproduce the CPU/MPI comparison on this machine or a similar CPU-only
machine:

1. Build a Release MPI executable with the same pgen:

   ```bash
   cmake -S . -B build-cr-perf-mpi \
     -D Athena_ENABLE_MPI=ON \
     -D PROBLEM=part_random \
     -D CMAKE_BUILD_TYPE=Release
   cmake --build build-cr-perf-mpi -j 8
   ```

2. Run each input at one, two, and four ranks, three times each.  Use
   `mpiexec -n 1` for the one-rank entry if you want to match the follow-up
   matrix exactly:

   ```bash
   mpiexec -n 1 ./build-cr-perf-mpi/src/athena \
     -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
     -d run-cr-drift-n1 \
     particles/amr_remap=device_table \
     particles/exchange_mode=alltoall_counts

   mpiexec -n 2 ./build-cr-perf-mpi/src/athena \
     -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
     -d run-cr-drift-n2 \
     particles/amr_remap=device_table \
     particles/exchange_mode=alltoall_counts

   mpiexec -n 4 ./build-cr-perf-mpi/src/athena \
     -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
     -d run-cr-drift-n4 \
     particles/amr_remap=device_table \
     particles/exchange_mode=alltoall_counts
   ```

   Repeat the same pattern for
   `inputs/particles/cr_tracer_boris_amr_perf.athinput`.

3. Repeat each rank/input pair with these three fallback override sets:

   ```text
   particles/amr_remap=host_tree particles/exchange_mode=alltoall_counts
   particles/amr_remap=device_table particles/exchange_mode=allgather
   particles/amr_remap=host_tree particles/exchange_mode=allgather
   ```

4. Record:

   - AthenaK `cpu time used`,
   - subprocess wall time,
   - `MeshBlocks created`,
   - `MeshBlocks deleted`,
   - final live MeshBlock count,
   - particle counts from `scripts/particles/cr_tracer_inspect.py` when outputs
     are enabled.

5. Compare medians, not individual fastest runs.

6. Keep correctness validation separate from performance timing unless the
   goal is to measure checker overhead.

## Benchmark Harness

The reusable benchmark driver is:

```bash
python scripts/particles/cr_tracer_benchmark.py \
  --athena build-cr-robust-mpi/src/athena \
  --cases push,amr_remap,exchange,histogram,histogram_per_rank,restart,boris_amr \
  --modes default,host_tree,allgather,fallback \
  --ranks 1,2,4 \
  --repeats 3 \
  --mpi-n1 \
  --output-json cr_tracer_bench_cpu_mpi.json \
  --output-md cr_tracer_bench_cpu_mpi.md
```

The harness runs existing inputs with command-line overrides, so benchmark
sweeps do not require editing production inputs.  It records machine/compiler
metadata available from the local environment, AthenaK CPU time, subprocess
wall time, particle update rate, MeshBlock activity, restart readback totals,
and lightweight particle counters.

The supported benchmark cases are:

| Case | Purpose |
|------|---------|
| `push` | Drift push throughput without AMR. |
| `amr_remap` | Zero-velocity AMR remap throughput. |
| `exchange` | High-velocity subcycled drift exchange throughput. |
| `histogram` | Reduced histogram, spectrum, and moment output throughput. |
| `histogram_per_rank` | Per-rank histogram, spectrum, and moment output throughput. |
| `restart` | Particle restart write and inspector read throughput. |
| `boris_amr` | End-to-end Boris, MHD, AMR, and particle migration. |

The supported runtime modes are:

| Mode | Overrides |
|------|-----------|
| `default` | `amr_remap=device_table`, `exchange_mode=alltoall_counts` |
| `host_tree` | `amr_remap=host_tree`, `exchange_mode=alltoall_counts` |
| `allgather` | `amr_remap=device_table`, `exchange_mode=allgather` |
| `fallback` | `amr_remap=host_tree`, `exchange_mode=allgather` |

For particle exchange, `LogPerformance()` now reports estimated MPI message and
byte counts.  The byte estimate is based on packed particle real and integer
fields sent by each rank.  Host-device copy volume is intentionally not
reported on this CPU-only machine; it is part of the GPU follow-up checklist.

## GitHub Pages Integration

This file is written as a standalone Markdown page so it can be copied directly
into the GitHub Pages documentation tree.  In this development branch it is
marked `orphan: true` because the branch only carries module-page fragments,
not the full Pages `index` and toctree.  The user-facing particle page links to
it with a relative Markdown link, so the link remains valid when all CR tracer
module pages live side by side in `docs/source/modules/`.

For the Pages branch, copy or merge these files together:

- `docs/source/modules/particles.md`,
- `docs/source/modules/cr_tracer_cpu_mpi_performance.md`,
- `docs/source/modules/cr_tracer_gpu_testing_handoff.md`.

The implementation-plan pages are useful review artifacts but are not required
for the user-facing CR tracer documentation.  If the Pages branch already has a
module toctree entry for `particles.md`, no extra toctree entry is required for
this performance page unless the maintainers want it to appear as a separate
sidebar item.

Before publishing, run the normal Pages documentation build from the Pages
worktree and check that:

- the relative link from `particles.md` to this page resolves,
- the relative link from this page to the GPU handoff resolves,
- wide benchmark tables render acceptably on desktop,
- the GPU handoff remains clearly labeled as TODO work rather than completed
  validation.

## Conclusions

The CPU/MPI evidence supports keeping the current branch design:

- The mesh-tree AMR lookup is the main performance fix.  It removes the old
  scan bottleneck and gives 22.6x-147.4x CPU-time speedup in a remap-focused
  test.
- In a realistic Boris+MHD AMR test, the full branch is 1.37x-1.39x faster by
  AthenaK CPU time than the old scan/full-copy path.
- The position-only host copy is not measurable in the small CPU Boris case,
  but it removes unnecessary remap data movement and is the correct direction
  for future GPU/device-side work.
- The reduced diagnostic default improves usability and file count without a
  measurable local runtime penalty.
- The consistency-check modes have low overhead in the validated small stress
  test, but `local` and `full` should remain debug/test modes because they copy
  and gather particle data.
- The position-based multilevel ownership update fixes a four-rank AMR
  out-of-block particle failure that the old immediate-neighbor stencil missed.
- The follow-up `device_table` path is faster than `host_tree` in the
  remap-focused drift matrix and neutral to 1.02x faster in the Boris+MHD
  matrix.
- The follow-up `alltoall_counts` path is neutral within local noise in the
  drift, exchange, and Boris+MHD matrices.  It should remain the default
  because it removes detailed all-gather metadata from the exchange design.
- Reusable buffers and prefix-sum compaction are included in the current runs
  but are not independently isolated by runtime switches.  The right claim is
  that they preserve CPU/MPI behavior while removing scaling hazards.
- The next CPU/MPI performance targets are higher-rank exchange scaling, larger
  particle-count scaling, allocation/profiler evidence for buffer reuse,
  compaction profiling, and interpolation/pusher profiling separated from AMR
  remap timing.

## GPU Follow-Up

These results do not cover GPU execution.  This machine only provides CPU and
local MPI testing, so accelerator claims must wait for a GPU machine.

The GPU follow-up should repeat the same comparison with CUDA, HIP, or SYCL
builds on a machine with accelerators and should add device-copy timing for:

- AMR remap before and after the position-only host copy,
- `device_table` against `host_tree` on real accelerator memory spaces,
- `alltoall_counts` against `allgather` with GPU-aware MPI enabled and
  disabled, if the machine supports both modes,
- reusable buffers and prefix-sum compaction against profiler counters for
  allocation, host/device synchronization, and kernel time,
- `df`, `pspec`, `dxh`, `drh`, `dparh`, and `pmom` diagnostics,
- `check_consistency_mode = counts`, `local`, and `full`,
- any future device-side or tree-cached remap implementation.

The minimum GPU TODO list is:

1. Build the same branch with the production GPU backend and MPI enabled.
2. Run the remap-focused drift AMR benchmark at one GPU and at multiple MPI
   ranks if the machine supports GPU-aware MPI.
3. Run the Boris+MHD AMR benchmark with the same inputs and particle counts.
4. Repeat the reduced-diagnostic and per-rank-diagnostic comparison.
5. Time `check_consistency_mode = none`, `counts`, `local`, and `full`.
6. Record host-device copy volume or profiler counters for AMR remap and
   diagnostics.
7. Repeat the four-mode matrix from this page:
   default, host-tree only, allgather only, and both fallback.
8. Confirm that any device-side remap implementation agrees with the host tree
   validation path.
9. Document the GPU model, driver/runtime, Kokkos backend, MPI implementation,
   compiler, rank/GPU mapping, and whether GPU-aware MPI was enabled.

The detailed execution checklist is in
[CR Tracer GPU Testing Handoff](cr_tracer_gpu_testing_handoff.md).  Until that
work is done, this page should be cited only for CPU/MPI behavior.
