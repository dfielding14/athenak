---
orphan: true
---

# CR Tracer CPU/MPI Performance Comparison

This page records the local CPU and MPI performance comparison for the CR tracer
particle hardening work on `feature/CR_tracers`.  It is intentionally limited to
the hardware available on 2026-05-22: CPU execution with the Kokkos Serial
backend and Open MPI.  GPU timing is still required before making claims about
accelerator performance.

## Scope

The comparison targets the branch-local performance fixes:

- replacing the host-side `O(N_particles * N_meshblocks)` AMR destination scan
  with mesh-tree lookup,
- reducing `Particles::RemapAfterAMR()` host mirrors from all particle real and
  integer fields to only `IPX`, `IPY`, `IPZ`, and `PGID`,
- preserving MPI particle exchange after AMR remap,
- making reduced `df`, `dxh`, `drh`, `dparh`, and `pmom` diagnostics the MPI
  default while keeping per-rank diagnostic files as an option,
- splitting consistency checks into `none`, `counts`, `local`, and `full`.

The most important result is that the mesh-tree AMR remap is a large win when
remapping is exposed directly, and it remains a measurable end-to-end win in a
real Boris+MHD AMR run.  The later position-only host copy is neutral within
noise for the small CPU Boris benchmark; it is still the right implementation
because it removes unnecessary data movement from the AMR remap path.

## Test Environment

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

Each table reports the median of three runs.  `cpu time used` and
`particle-updates/cpu_second` are AthenaK's own runtime counters; `wall` is the
outer Python subprocess wall time.  The absolute numbers are local benchmark
data, not portable performance guarantees.

## Benchmark Inputs

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
| 1 | 0.5686 | 0.00511 | 111.3x | 0.6177 | 0.04997 | 12.4x |
| 2 | 0.4210 | 0.01261 | 33.4x | 0.4841 | 0.07574 | 6.4x |
| 4 | 0.2786 | 0.01696 | 16.4x | 0.3503 | 0.08825 | 4.0x |

Interpretation:

- The old scan is dominated by repeatedly testing particle positions against
  the global MeshBlock list.
- The tree path removes that scaling and makes MPI overhead, load balancing,
  and process launch a larger fraction of the small benchmark.
- The CPU-time speedup decreases with rank count because the old scan work is
  divided across ranks while fixed MPI overhead becomes more visible.
- The wall-clock speedup is smaller than the AthenaK CPU-time speedup because
  these short runs include process launch and setup overhead.

## Boris+MHD End-To-End AMR Benchmark

This benchmark includes MHD evolution, Boris particle pushing, AMR, particle
migration, and no diagnostic output.  All runs produced the same AMR activity:
588 MeshBlocks created, 476 deleted, and 176 live MeshBlocks at termination.

| MPI ranks | Old scan/full-copy CPU s | Current CPU s | CPU speedup | Old wall s | Current wall s | Wall speedup |
|-----------|--------------------------|---------------|-------------|------------|----------------|--------------|
| 1 | 1.1082 | 0.8325 | 1.33x | 1.1742 | 0.9007 | 1.30x |
| 2 | 0.7437 | 0.5692 | 1.31x | 0.8469 | 0.6740 | 1.26x |
| 4 | 0.4789 | 0.3743 | 1.28x | 0.5828 | 0.4801 | 1.21x |

The end-to-end speedup is smaller than in the remap microbenchmark because MHD,
Boris interpolation, and particle push work remain unchanged.  This is the more
realistic expectation for production CPU runs: the fast remap removes a clear
AMR bottleneck, but it does not change the cost of the physical update.

## Position-Only Host-Copy Check

The `b12e02c1` implementation used the mesh-tree lookup but still mirrored the
full particle arrays during AMR remap.  The current implementation mirrors only
positions and `PGID`.  On the small CPU Boris benchmark this later copy
optimization is neutral within run-to-run noise:

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

## Diagnostic Output Reduction

The output comparison used the stress input promoted to 32^3 cells, `ppc = 1`,
four MPI ranks, and two diagnostic dumps.  The reduced default and per-rank
mode had indistinguishable runtime on this local filesystem, but very different
file counts:

| Mode | CPU s | Wall s | `df` files | `dxh` files | `drh` files | `dparh` files | `pmom` files |
|------|-------|--------|------------|-------------|-------------|---------------|--------------|
| Reduced default | 0.2091 | 0.3356 | 1 | 1 | 1 | 1 | 1 |
| Per-rank files | 0.2068 | 0.3315 | 4 | 4 | 4 | 4 | 4 |

The performance difference is not significant at this scale.  The reduced
default is still the better user-facing behavior because each diagnostic has one
global file, histogram sums directly match global species counts, and analysis
scripts do not need to merge rank-local histograms.  Per-rank diagnostic files
remain useful for debugging rank-local particle migration and can be enabled
with:

```text
df_single_file_per_rank = 1
dxhist_single_file_per_rank = 1
drh_single_file_per_rank = 1
dparh_single_file_per_rank = 1
pmom_single_file_per_rank = 1
```

## Conclusions

The CPU/MPI evidence supports keeping the current branch design:

- The mesh-tree AMR lookup is the main performance fix.  It removes the old
  scan bottleneck and gives 16x-111x CPU-time speedup in a remap-focused test.
- In a realistic Boris+MHD AMR test, the full branch is 1.28x-1.33x faster by
  AthenaK CPU time than the old scan/full-copy path.
- The position-only host copy is not measurable in the small CPU Boris case,
  but it removes unnecessary remap data movement and is the correct direction
  for future GPU/device-side work.
- The reduced diagnostic default improves usability and file count without a
  measurable local runtime penalty.
- The consistency-check modes have low overhead in the validated small stress
  test, but `local` and `full` should remain debug/test modes because they copy
  and gather particle data.

## GPU Follow-Up

These results do not cover GPU execution.  The GPU follow-up should repeat the
same comparison with CUDA/HIP/SYCL builds on a machine with accelerators and
should add device-copy timing for:

- AMR remap before and after the position-only host copy,
- `df`, `dxh`, `drh`, `dparh`, and `pmom` diagnostics,
- `check_consistency_mode = counts`, `local`, and `full`,
- a future device-side or tree-cached remap implementation.
