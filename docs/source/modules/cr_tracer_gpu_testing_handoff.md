---
orphan: true
---

# CR Tracer GPU Testing Handoff

## Purpose

This document is the exact follow-up checklist for validating the CR tracer
architecture branch on a GPU machine.  The current machine only supports CPU
and local MPI testing, so no accelerator performance or correctness claim is
complete until the steps below have been run elsewhere.

The GPU validation should be run after the CPU/MPI branch passes locally.  Use
the same branch and commit that passed CPU/MPI tests so any GPU failure can be
attributed to backend, memory-space, MPI, or scaling behavior rather than to
unverified CPU changes.

## Required Machine Metadata

Record this information at the top of the GPU result note or pull-request
comment:

- hostname and machine type,
- GPU model and number of GPUs,
- CPU model,
- compiler and version,
- MPI implementation and version,
- Kokkos backend: CUDA, HIP, SYCL, or other,
- GPU driver and runtime version,
- rank-to-GPU mapping,
- whether GPU-aware MPI is enabled,
- branch name and commit hash,
- exact CMake command.

## Build Commands

CUDA example:

```bash
cmake -S . -B build-cr-gpu-mpi \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  -D PROBLEM=part_random \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build-cr-gpu-mpi -j 8
```

HIP example:

```bash
cmake -S . -B build-cr-hip-mpi \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_HIP=ON \
  -D PROBLEM=part_random \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build-cr-hip-mpi -j 8
```

Add site-specific architecture flags required by the machine, for example a
`Kokkos_ARCH_*` flag or compiler wrapper.  Record every added flag.

## Correctness Smoke Tests

Run the moving AMR stress test first.  It exercises the device-table AMR remap,
particle MPI exchange, device-style compaction, reduced diagnostics, restart
output, particle-aware load-balancing cost hooks, and refine/derefine particle
conservation.

```bash
rm -rf run-cr-gpu-smoke
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_stress.athinput \
  -d run-cr-gpu-smoke \
  particles/check_consistency_mode=full \
  particles/validate_amr_lookup=true \
  mesh_refinement/particle_load_weight=0.001

python scripts/particles/cr_tracer_inspect.py run-cr-gpu-smoke \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
```

Acceptance criteria:

- run exits with status 0,
- AMR reports `MeshBlocks created > 0`,
- AMR reports `MeshBlocks deleted > 0`,
- inspector reports 1024 total particles,
- inspector reports species counts `0=512 1=512`,
- no invalid `PGID`, out-of-block particle, duplicate tag, or malformed restart
  is accepted,
- a separate high-velocity drift smoke also passes with `particles/subcycle=true`,
- reduced `df`, `pspec`, `dxh`, `drh`, `dparh`, and `pmom` files are finite
  and have the expected totals.

Then run the subcycling smoke.  This mirrors the local CPU/MPI test that forces
intermediate ownership updates during one particle push.

```bash
rm -rf run-cr-gpu-subcycle
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/random_particle_drift.athinput \
  -d run-cr-gpu-subcycle \
  job/basename=cr_gpu_subcycle \
  mesh/nx1=32 mesh/nx2=8 mesh/nx3=8 \
  meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=8 \
  particles/ppc=0.0009765625 \
  particles/nspecies=1 \
  particles/check_consistency_mode=local \
  particles/check_motion_bounds=true \
  particles/subcycle=true \
  particles/subcycle_max_steps=256 \
  problem/particle_position=center \
  problem/particle_velocity=uniform \
  problem/v0x=100.0 problem/v0y=0.0 problem/v0z=0.0 \
  time/nlim=1 time/tlim=0.125 \
  output1/file_type=prst output1/dt=0.01

python scripts/particles/cr_tracer_inspect.py run-cr-gpu-subcycle \
  --expected-total 2 \
  --expected-species-counts 2
```

## Backend Comparison Runs

Run these paired tests to separate backend issues from new architecture issues.

New default path:

```bash
rm -rf run-cr-gpu-default
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_stress.athinput \
  -d run-cr-gpu-default \
  particles/amr_remap=device_table \
  particles/exchange_mode=alltoall_counts \
  particles/check_consistency_mode=full \
  particles/validate_amr_lookup=true
```

Fallback path:

```bash
rm -rf run-cr-gpu-fallback
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_stress.athinput \
  -d run-cr-gpu-fallback \
  particles/amr_remap=host_tree \
  particles/exchange_mode=allgather \
  particles/check_consistency_mode=full
```

Validate both outputs:

```bash
python scripts/particles/cr_tracer_inspect.py run-cr-gpu-default \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8

python scripts/particles/cr_tracer_inspect.py run-cr-gpu-fallback \
  --expected-total 1024 \
  --expected-species-counts 512,512 \
  --histogram-nbin 8 \
  --spectrum-nbin 8 \
  --joint-spectrum-nbin 8,8
```

The default and fallback runs do not need bitwise-identical diagnostic moments,
but they must conserve total and per-species particle counts and pass restart
and histogram validation.

## Performance Benchmarks

Run the same benchmark harness used in the local CPU/MPI documentation.  Use
one rank/GPU first, then a rank count that maps cleanly to the available GPUs.
The CPU/MPI documentation compares four runtime modes; repeat the same matrix
on GPU before making accelerator claims:

| Mode | Overrides |
|------|-----------|
| Default | `particles/amr_remap=device_table particles/exchange_mode=alltoall_counts` |
| Host-tree only | `particles/amr_remap=host_tree particles/exchange_mode=alltoall_counts` |
| Allgather only | `particles/amr_remap=device_table particles/exchange_mode=allgather` |
| Both fallback | `particles/amr_remap=host_tree particles/exchange_mode=allgather` |

The CPU/MPI result says the new defaults are neutral to modestly faster on this
machine.  The GPU result still needs to answer the separate question: whether
the device-table remap, reduced host mirrors, reusable buffers, and prefix-sum
compaction reduce host/device copies, synchronization, and kernel time on real
accelerator memory spaces.  It also needs to verify that subcycled ownership
exchange and particle-aware AMR load balancing do not introduce backend-specific
failures.

Harness command template:

```bash
python scripts/particles/cr_tracer_benchmark.py \
  --athena build-cr-gpu-mpi/src/athena \
  --cases push,amr_remap,exchange,histogram,histogram_per_rank,restart,boris_amr \
  --modes default,host_tree,allgather,fallback \
  --ranks 1,2 \
  --repeats 3 \
  --mpi-n1 \
  --output-json cr_tracer_bench_gpu_mpi.json \
  --output-md cr_tracer_bench_gpu_mpi.md
```

Increase `--ranks` to match the GPU allocation and rank/GPU mapping.  Keep the
JSON output as the machine-readable artifact and use the Markdown output in the
PR or benchmark note.

Manual remap-focused AMR benchmark:

```bash
rm -rf run-cr-gpu-drift-n1
./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
  -d run-cr-gpu-drift-n1 \
  particles/check_consistency_mode=none

rm -rf run-cr-gpu-drift-n2
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_drift_amr_perf.athinput \
  -d run-cr-gpu-drift-n2 \
  particles/check_consistency_mode=none
```

Boris+MHD AMR benchmark:

```bash
rm -rf run-cr-gpu-boris-n1
./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_perf.athinput \
  -d run-cr-gpu-boris-n1 \
  particles/check_consistency_mode=none

rm -rf run-cr-gpu-boris-n2
mpiexec -n 2 ./build-cr-gpu-mpi/src/athena \
  -i inputs/particles/cr_tracer_boris_amr_perf.athinput \
  -d run-cr-gpu-boris-n2 \
  particles/check_consistency_mode=none
```

For each benchmark, repeat at least three times and report medians for:

- AthenaK `cpu time used`,
- outer wall time,
- `particle-updates/cpu_second`,
- MeshBlocks created,
- MeshBlocks deleted,
- final live MeshBlock count,
- total particles and species counts if outputs are enabled.

## Profiling Requirements

Use the site profiler appropriate to the backend.  Examples include Nsight
Systems for CUDA, rocprof for HIP, and the Kokkos profiling tools where
available.

Profile at least:

- `Particles::RemapAfterAMR()` device-table kernel,
- host-tree fallback remap,
- particle MPI exchange pack/unpack,
- compaction after sends exceed receives,
- intermediate exchange during particle subcycling,
- AMR load balancing with `particle_load_weight > 0`,
- reduced diagnostics, including `df`, `pspec`, `dxh`, `drh`, `dparh`, and
  `pmom`,
- consistency modes `counts`, `local`, and `full`.

Record:

- host-to-device copy volume,
- device-to-host copy volume,
- remap kernel time,
- multilevel post-push `NewGID` lookup-table time,
- particle pack/unpack kernel time,
- MPI wait time if available,
- allocation/reallocation count if visible,
- diagnostic reduction time.

## Failure Triage

If the default GPU path fails but the fallback path passes:

1. Rerun the default path with `particles/validate_amr_lookup=true`.
2. Check whether the first failure is remap lookup, particle exchange, or
   compaction.
3. Rebuild with `-D CMAKE_BUILD_TYPE=Debug` and Kokkos bounds checking if the
   site can tolerate the overhead.
4. Rerun the smallest failing input with one rank and then two ranks.
5. Save stdout, stderr, the input file, and the exact CMake command.

If both default and fallback fail:

1. Confirm the branch still passes the CPU/MPI stress test.
2. Check whether the pgen or output path is using unsupported host-only logic.
3. Confirm that GPU-aware MPI configuration matches the site recommendation.
4. Try a one-rank run to separate GPU execution from MPI exchange.

## Completion Criteria

GPU validation is complete only when:

- the GPU+MPI build succeeds,
- the default device-table path passes the AMR stress run,
- the host-tree/allgather fallback passes the same stress run,
- restart and diagnostic outputs pass `cr_tracer_inspect.py`,
- remap-focused and Boris+MHD performance medians are recorded,
- profiler evidence records host-device copy volume and remap kernel time,
- a high-velocity subcycle smoke passes on GPU,
- an AMR stress run with `particle_load_weight > 0` passes on GPU,
- results are written into the CPU/MPI performance page or a linked GPU
  performance page.
