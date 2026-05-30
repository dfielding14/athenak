---
orphan: true
---

# CR Tracer Relativistic GPU Testing Handoff

## Status

`GPU QUALIFIED` is **not claimed** for the relativistic CR tracer acceleration
branch.  The local qualification workstation is an Apple M4 Max system with the
Kokkos Serial backend.  It exposes Metal graphics, but it does not provide a
configured CUDA or HIP Kokkos backend, `nvcc`, `hipcc`, `nvidia-smi`, or
`rocminfo`.

GPU qualification is optional for the workstation `MERGE READY` decision and
remains an explicit residual risk.  Run this handoff on the accepted branch tip
before making any `GPU QUALIFIED` claim.

## Scope

The accelerator follow-up must test the already reviewed passive full-orbit
implementation without widening physics:

- `pusher = relativistic_hc`;
- authoritative `w = gamma v` state;
- Higuera-Cary push;
- `prescribed_test` analytical execution;
- serial solver-coupled `mhd_ideal` execution with
  `relativistic_temporal_sampling = frozen_tn`;
- acceleration-aware subcycling;
- typed-v2 restart;
- retained diagnostics and inspector round trips;
- bounded periodic MPI, static-SMR, and adaptive-AMR migration for
  `prescribed_test`.

Do not use GPU qualification to open solver-coupled `mhd_ideal` MPI, SMR, or AMR
execution, nonperiodic boundaries, changed-rank restart redistribution,
stage-coupled temporal sampling, or fluid feedback.

## Backend Record

Before building, archive:

```bash
git rev-parse HEAD
git status --short
cmake --version
command -v nvcc || true
command -v hipcc || true
nvidia-smi || true
rocminfo || true
```

Record the exact compiler, GPU model, driver, MPI implementation, Kokkos
backend, architecture flags, node count, rank placement, and build cache.

## Build Matrix

Use the backend appropriate to the target machine.  CUDA examples:

```bash
cmake -S . -B build-cr-rel-gpu-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON
cmake --build build-cr-rel-gpu-release -j 4

cmake -S . -B build-cr-rel-gpu-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build build-cr-rel-gpu-debug -j 4

cmake -S . -B build-cr-rel-gpu-mpi-release \
  -D CMAKE_BUILD_TYPE=Release \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON
cmake --build build-cr-rel-gpu-mpi-release -j 4
```

For HIP, replace `Kokkos_ENABLE_CUDA` with `Kokkos_ENABLE_HIP` and record the
architecture selection explicitly.

Build the dedicated analytical binaries as well:

```bash
cmake -S . -B build-cr-rel-gpu-pusher \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_pusher_runtime_test \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON
cmake --build build-cr-rel-gpu-pusher -j 4

cmake -S . -B build-cr-rel-gpu-coupled \
  -D CMAKE_BUILD_TYPE=Release \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON
cmake --build build-cr-rel-gpu-coupled -j 4

cmake -S . -B build-cr-rel-gpu-coupled-debug \
  -D CMAKE_BUILD_TYPE=Debug \
  -D PROBLEM=unit_tests/cr_relativistic_coupled_runtime_test \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_DEBUG=ON \
  -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build build-cr-rel-gpu-coupled-debug -j 4
```

## Required Replays

Run the tracked analyzers with fresh work directories:

```bash
python3 scripts/particles/cr_relativistic_pusher_runtime_inspect.py \
  --binary build-cr-rel-gpu-pusher/src/athena \
  --input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4a_preregistered_criteria.json \
  --work-dir run-cr-rel-gpu-pusher \
  --metrics cr-rel-gpu-pusher-metrics.json

python3 scripts/particles/cr_relativistic_coupled_runtime_inspect.py \
  --binary build-cr-rel-gpu-coupled/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --prescribed-input inputs/unit_tests/cr_relativistic_pusher_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase4b_preregistered_criteria.json \
  --work-dir run-cr-rel-gpu-coupled \
  --metrics cr-rel-gpu-coupled-metrics.json

python3 scripts/particles/cr_relativistic_subcycle_runtime_inspect.py \
  --binary build-cr-rel-gpu-coupled/src/athena \
  --input inputs/unit_tests/cr_relativistic_coupled_runtime_test.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase5_preregistered_criteria.json \
  --work-dir run-cr-rel-gpu-subcycle \
  --metrics cr-rel-gpu-subcycle-metrics.json
```

Run the restart and diagnostic analyzers with the general GPU binary:

```bash
python3 scripts/particles/cr_relativistic_restart_runtime_inspect.py \
  --binary build-cr-rel-gpu-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_restart_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase6_preregistered_criteria.json \
  --work-dir run-cr-rel-gpu-restart \
  --json cr-rel-gpu-restart-metrics.json

python3 scripts/particles/cr_relativistic_diagnostics_runtime_inspect.py \
  --binary build-cr-rel-gpu-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --work-dir run-cr-rel-gpu-diagnostics \
  --metrics cr-rel-gpu-diagnostics-metrics.json

python3 scripts/particles/cr_relativistic_all_formats_inspect.py \
  --runtime-root run-cr-rel-gpu-diagnostics/acceleration_uninterrupted \
  --binary build-cr-rel-gpu-release/src/athena \
  --input inputs/unit_tests/cr_relativistic_diagnostics_runtime.athinput \
  --criteria docs/source/modules/figures/cr_tracer_relativistic_acceleration/cr_relativistic_phase7_preregistered_criteria.json \
  --metrics cr-rel-gpu-diagnostics-metrics.json \
  > cr-rel-gpu-all-formats.json
```

Repeat the coupled, subcycle, restart, diagnostic, and all-format replays with
the debug-bounds accelerator binaries.  Keep release and debug-bounds reports
separate so a release-only pass cannot satisfy the accelerator exit contract.
Use `build-cr-rel-gpu-coupled-debug/src/athena` for the coupled and subcycle
debug-bounds replays, and `build-cr-rel-gpu-debug/src/athena` for restart,
diagnostic, and all-format debug-bounds replays.  A dedicated pusher-debug
replay and MPI migration debug-bounds replay are recommended strengthening
checks; they become mandatory if release/debug differential behavior appears
or if the target backend's launch configuration differs between build types.

## MPI And AMR Gate

Run the bounded prescribed-test migration qualifier with the MPI GPU binary:

```bash
python3 scripts/particles/cr_relativistic_migration_runtime_inspect.py \
  --binary build-cr-rel-gpu-release/src/athena \
  --mpi-binary build-cr-rel-gpu-mpi-release/src/athena \
  --work-dir run-cr-rel-gpu-migration \
  --json cr-rel-gpu-migration-metrics.json
```

This must preserve the existing `62/62` acceptance contract, including:

- MPI rank ladder `1`, `2`, `4`, and `8`;
- particle-empty-rank witness;
- periodic wrap;
- static-SMR cross-level transitions through both `device_table` and
  `host_tree`;
- adaptive refine/derefine with off-rank remap sends for both remap modes;
- restart-after-AMR continuation for both remap modes;
- deterministic changed-rank restart rejection.

Run a GPU MPI smoke where the machine supports device-aware MPI.  If MPI stages
through host memory, record that fact and keep the result distinct from a
device-aware transport claim.

## Differential Review

Compare CPU and GPU metrics on the same accepted SHA:

- analytical errors and registered thresholds;
- acceleration and deceleration work closure;
- sampled-field signs and manufactured-field convergence;
- subcycle activation and outer-timestep bounds;
- restart byte counts, checksums, and decoded state;
- output metadata and all-format inspector JSON;
- migration identity, cross-level lookup counts, AMR churn, and restart
  continuation parity.

Do not explain a backend mismatch as floating-point noise until the compared
quantity, expected tolerance, execution backend, and reproducibility have been
written down.

## Performance Pass

Profile before optimizing.  Record:

- total runtime;
- particle updates per second;
- remap sends and receives;
- diagnostic counters;
- kernel time for gather, pusher, subcycle-envelope, AMR lookup, pack, unpack,
  and restart output;
- host-device transfers;
- device memory footprint;
- MPI transport mode.

Keep any optimization in a separate commit with CPU/GPU differential replay.

## Exit Contract

Claim `GPU QUALIFIED` only after:

1. Release and debug accelerator builds pass.
2. Analytical, coupled, subcycle, restart, diagnostic, and all-format replays
   pass on one accepted SHA.
3. GPU MPI prescribed-test migration passes where supported.
4. CPU/GPU differentials are within written tolerances.
5. Profiling evidence is archived before optimization.
6. A fresh independent reviewer audits the accelerator evidence and records
   `PROCEED`.
