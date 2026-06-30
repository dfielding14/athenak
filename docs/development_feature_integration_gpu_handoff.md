# Development Feature Integration GPU Handoff

Date: 2026-06-30

## Branch State

Integration branch: `c/development-feature-integration`

Worktree used for integration:
`/Users/dbf75/.codex/worktrees/development-feature-integration/athenak-DF`

Code integration head before this handoff document:
`71420545ac5f329b9402f41f15640118856ab487`

Base: `origin/development` at `462a0dbcd`

Merged feature refs:

| Feature | Integrated ref |
| --- | --- |
| `feature/initial-perturbations` | `origin/feature/initial-perturbations` at `f96baef2e` |
| `feature/cooling` | `origin/feature/cooling` at `3eaaaf224` |
| `feature/turbulence-driving` | `origin/feature/turbulence-driving` at `6fa601cb9` |
| `feature/monte-carlo-thermo-tracers` | `origin/feature/monte-carlo-thermo-tracers` at `75fc522d6` |
| `feature/frame-tracker` | `origin/feature/frame-tracker` at `52a3427fc` |

First-parent integration commits:

```text
71420545a Merge frame tracker hardening into development integration
088cb34f0 Merge Monte-Carlo thermo tracers into development integration
1e29c3d35 Merge turbulence driving into development integration
896bc519d Merge cooling hardening into development integration
b29bd365f Merge initial perturbation hardening into development integration
462a0dbcd origin/development
```

Local checks already run:

```bash
git diff --cached --check
rg -n "^(<<<<<<<|=======|>>>>>>>)" .github docs inputs scripts src tst
```

Both checks passed before the frame-tracker merge commit.

## Pre-Push Local Compile And Smoke Audit

After the handoff was first written, a local compile audit found a missing
closing brace in `ProblemGenerator::ProblemGenerator(ParameterInput*, Mesh*,
IOWrapper, FileShardMode)` after the shared restart data block. The pre-push fix
restores that restart scope before MC particle restart data is read.

These local CPU/macOS builds passed after that fix:

```bash
cmake -S . -B /tmp/athenak_integration_build_builtin -DPROBLEM=built_in_pgens
cmake --build /tmp/athenak_integration_build_builtin -j 8

cmake -S . -B /tmp/athenak_integration_build_cloud -DPROBLEM=cloud_crushing
cmake --build /tmp/athenak_integration_build_cloud -j 8

cmake -S . -B /tmp/athenak_integration_build_trml -DPROBLEM=TRML_frame_tracking
cmake --build /tmp/athenak_integration_build_trml -j 8

cmake -S . -B /tmp/athenak_integration_build_turb -DPROBLEM=turb
cmake --build /tmp/athenak_integration_build_turb -j 8
```

These local serial smoke runs also passed:

```text
inputs/particles/lagrangian_mc_thermo.athinput
inputs/hydro/frame_tracking_smoke.athinput
inputs/tests/initial_perturbations.athinput
inputs/hydro/cloud_crushing_snr.athinput
inputs/hydro/TRML/TRML_frame_tracking_material.athinput
inputs/hydro/turb_driving/constant_edot_fixed_grid.athinput
```

GPU validation has not been run locally; the remaining validation target is the
GPU-machine matrix below.

## GPU Particle Counter Fix

The TRML-specific commit requested by the user was:

```text
d37a6551e0e20349c5e41b2487144e593624f91c Fix GPU particle counter allocation
```

This exact commit is not an ancestor of the integration branch because it also
contains a TRML-specific `configure.sh` change. The source fix is present via:

```text
75fc522d6da74a2fa0087ea544e7f07e66916cb9 Backport GPU particle counter allocation fix
```

Audit evidence:

```text
d37 source patch-id: b96164dfbeb615cb13d881c1b5ab60eb55db2dde
75f source patch-id: b96164dfbeb615cb13d881c1b5ab60eb55db2dde
75fc-to-HEAD source diff in affected files: 0
```

Affected source files:

```text
src/bvals/bvals.cpp
src/bvals/bvals.hpp
src/bvals/bvals_part.cpp
src/outputs/track_prtcl.cpp
```

The fix replaces host stack counters used inside GPU kernels with device-resident
counters for particle send lists and tracked particle output.

## Checkout Audit For GPU Agent

Run these first on the GPU machine:

```bash
git status --short
git rev-parse HEAD
git merge-base --is-ancestor origin/feature/cooling HEAD
git merge-base --is-ancestor origin/feature/turbulence-driving HEAD
git merge-base --is-ancestor origin/feature/monte-carlo-thermo-tracers HEAD
git merge-base --is-ancestor origin/feature/initial-perturbations HEAD
git merge-base --is-ancestor origin/feature/frame-tracker HEAD
git merge-base --is-ancestor 75fc522d6da74a2fa0087ea544e7f07e66916cb9 HEAD
git diff --name-only 75fc522d6 HEAD -- \
  src/bvals/bvals.cpp src/bvals/bvals.hpp \
  src/bvals/bvals_part.cpp src/outputs/track_prtcl.cpp
```

Expected: clean status, all `merge-base` commands exit 0, and the final diff
prints nothing.

## Build Setup

Use the site's normal Kokkos GPU toolchain. For NVIDIA/CUDA, the repo test
wrapper is usable because `tst/run_test_suite.py --gpu` adds
`Kokkos_ENABLE_CUDA=On`.

For AMD/HIP systems such as Frontier, do not use `run_test_suite.py --gpu`
directly because it hardcodes CUDA. Configure directly instead, adapting the
compiler wrapper and architecture to the site environment:

```bash
cmake -S . -B build_gpu_hip \
  -D CMAKE_BUILD_TYPE=Release \
  -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_HIP=ON \
  -D Kokkos_ARCH_AMD_GFX90A=ON \
  -D CMAKE_CXX_COMPILER=hipcc

cmake --build build_gpu_hip -j 8
export ATHENA=$PWD/build_gpu_hip/src/athena
```

If the site requires Cray wrappers, replace the compiler line with the local
wrapper, for example `-D CMAKE_CXX_COMPILER=CC`.

## Validation Matrix

Run from the repository root unless noted. Use `srun`/launcher options required
by the GPU system.

### 1. Built-In GPU Smoke Tests

CUDA systems:

```bash
cd tst
python run_test_suite.py --test test_suite/cooling/test_cooling_gpu.py --gpu
python run_test_suite.py --test test_suite/nr/test_nr_lwave3d_amr_gpu.py --gpu
cd ..
```

HIP systems should either port the wrapper flags locally or run equivalent
direct executable commands with the HIP build.

Acceptance: tests pass without Kokkos allocation, bounds, or illegal-address
failures.

### 2. Cooling GPU Path

```bash
rm -rf run_gpu_cooling
srun -n 1 "$ATHENA" \
  -i tst/inputs/cooling_test_hydro.athinput \
  -d run_gpu_cooling \
  job/basename=cooling_gpu_direct \
  cooling/cooling_model=table \
  cooling/cooling_table=tst/inputs/tables/cooling_lambda_3d.tbl \
  problem/scalar0=1.0
```

Acceptance: run exits 0 and writes `run_gpu_cooling/cooling_gpu_direct.hydro.hst`
with finite cooling history columns.

### 3. Turbulence Driving GPU Path

```bash
rm -rf run_gpu_turbulence
srun -n 1 "$ATHENA" \
  -i inputs/hydro/turb_driving/constant_edot_fixed_grid.athinput \
  -d run_gpu_turbulence \
  time/nlim=12
```

Acceptance: run exits 0, writes force and primitive binary outputs, and does not
report nonfinite acceleration or restart-state errors.

### 4. Initial Perturbations GPU Path

```bash
rm -rf run_gpu_initial_perturbations
srun -n 1 "$ATHENA" \
  -i inputs/tests/initial_perturbations.athinput \
  -d run_gpu_initial_perturbations
```

Acceptance: run exits 0 at cycle 0, writes `InitialPerturbations.mhd_w_bcc`
VTK output, and reports no Kokkos bounds or illegal-address failures.

### 5. Monte-Carlo Thermodynamic Tracers And Particle Counter Fix

Serial GPU:

```bash
rm -rf run_gpu_mc_serial
srun -n 1 "$ATHENA" \
  -i inputs/particles/lagrangian_mc_thermo.athinput \
  -d run_gpu_mc_serial

python scripts/read_prtcl_thermo_history.py \
  run_gpu_mc_serial/prtcl_thermo_history/lagrangian_mc_thermo.prtcl_thermo_history.thp \
  --npz run_gpu_mc_serial/tracers.npz
```

MPI+AMR GPU:

```bash
rm -rf run_gpu_mc_amr
srun -n 2 "$ATHENA" \
  -i inputs/particles/lagrangian_mc_thermo_amr.athinput \
  -d run_gpu_mc_amr

srun -n 2 "$ATHENA" \
  -r run_gpu_mc_amr/rst/rank_00000000/lagrangian_mc_thermo_amr.00001.rst \
  -d run_gpu_mc_amr \
  time/tlim=0.04

python scripts/read_prtcl_thermo_history.py \
  run_gpu_mc_amr/prtcl_thermo_history/lagrangian_mc_thermo_amr.prtcl_thermo_history.thp \
  --npz run_gpu_mc_amr/tracers.npz
```

Acceptance: both runs exit 0, restart append succeeds, `.thp` files are
nonempty and readable, tracer columns are finite, and there are no GPU failures
from particle send-list or tracked-particle counters.

### 6. Frame Tracking GPU Path

Built-in pgen smoke:

```bash
rm -rf run_gpu_frame_smoke
srun -n 1 "$ATHENA" \
  -i inputs/hydro/frame_tracking_smoke.athinput \
  -d run_gpu_frame_smoke \
  time/nlim=20
```

Frame-aware cloud material input:

```bash
rm -rf run_gpu_cloud_material
srun -n 1 "$ATHENA" \
  -i inputs/hydro/cloud_crushing_material_tracking.athinput \
  -d run_gpu_cloud_material \
  mesh/nx1=32 mesh/nx2=16 mesh/nx3=16 \
  meshblock/nx1=16 meshblock/nx2=8 meshblock/nx3=8 \
  time/nlim=80 time/tlim=0.004 time/cfl_number=0.01 \
  output1/dt=2.0e-5 output2/dt=5.0e-5 output3/dt=0.002 \
  frame_tracking/start_time=0.0 \
  frame_tracking/diagnostic_every=5
```

Frame-aware TRML material input:

```bash
rm -rf run_gpu_trml_material
srun -n 1 "$ATHENA" \
  -i inputs/hydro/TRML/TRML_frame_tracking_material.athinput \
  -d run_gpu_trml_material \
  mesh/nx1=32 mesh/nx2=32 mesh/nx3=64 \
  meshblock/nx1=16 meshblock/nx2=16 meshblock/nx3=32 \
  time/nlim=80 time/tlim=0.05 \
  output1/dt=0.01 output2/dt=0.025 output3/dt=0.025
```

Acceptance: runs exit 0, each writes `<basename>.frame_tracker.hst`, `data_precision=real`
binary output works for the material inputs, restart files are produced, and
there are no frame-state or boundary-refresh errors.

### 7. Combined Stress Run

This is the highest-value gate because it exercises AMR, MPI, particles,
restart append, and the GPU counter fix together.

```bash
rm -rf run_gpu_combined_stress
srun -n 2 "$ATHENA" \
  -i inputs/particles/lagrangian_mc_thermo_amr.athinput \
  -d run_gpu_combined_stress \
  time/nlim=12 time/tlim=0.05 \
  output1/dt=0.01 output2/dt=0.025
```

Acceptance: run exits 0 under the GPU backend, `.thp` output remains readable,
AMR changes complete, and the log contains no Kokkos illegal memory access,
allocation, or host-device counter errors.

## Report Back

For each gate, report:

- exact commit SHA tested,
- GPU backend and architecture,
- build command and launcher command,
- pass/fail status,
- first failing command and the relevant log tail if any,
- whether particle `.thp`, frame tracker `.hst`, restart, and binary outputs were produced.
