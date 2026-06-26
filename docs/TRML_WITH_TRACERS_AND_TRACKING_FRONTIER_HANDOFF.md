# TRML with tracers and tracking: Frontier handoff

This is the stopping point for local integration of
`TRML_with_Tracers_and_Tracking`. The combined problem has passed the serial,
MPI, AMR, and restart checks listed below on macOS. The HIP build and larger GPU
runs in this document have **not** been run on Frontier.

The canonical input is
`inputs/hydro/TRML/TRML_with_Tracers_and_Tracking.athinput`; build it with
`PROBLEM=simple_TRML`.

## Integrated behavior

- `simple_TRML.cpp` supplies the pressure-balanced shear layer, exact cooling,
  frame-aware x3 reservoirs, the cold-material scalar, and 24 user-history
  diagnostics.
- `<initial_perturbations>` supplies the single seeded velocity perturbation.
- `<frame_tracking>` follows the conserved cold-material scalar along x3.
- `lagrangian_mc` particles sample the fluid and move using saved mass fluxes.
- Of 3,072 canonical tracers, 777 start throughout the volume and 2,295 enter
  through a one-root-cell-thick top boundary slab in 51 equal batches from `t=0` through
  `t=20`.
- A frame update precedes the Monte-Carlo particle push. These particles have
  positions but no physical velocity state, so the frame boost is not applied a
  second time to their positions. Lab coordinates are reconstructed with
  `x3_lab = x3_grid + ft_dx_x3`; sampled velocities use
  `v3_lab = v3_grid + ft_vf_x3`.

Do not generalize that last statement to velocity-carrying particle species;
their interaction with frame tracking is not certified by this integration.

## Local evidence already collected

The following checks passed on the integration branch:

| Gate | Result |
| --- | --- |
| C++ and Python style | 2 passed |
| Combined population, RK2/RK3/RK4 cooling, serial restart, and AMR regression | 6 passed |
| Frame-tracker CPU regression | 20 passed |
| Initial-perturbation Hydro/MHD regression | 2 passed |
| Frame restart regression | 2 passed |
| Shipped frame-aware problem examples | 2 passed |
| Frame MPI/AMR regression | 3 passed |
| Thermodynamic-history reader regression | passed |

The combined manual matrix additionally covered:

- one-rank versus two-rank uniform runs;
- two-rank uninterrupted versus restart-split runs;
- serial versus two-rank adaptive runs, including refinement from 8 to 64
  MeshBlocks;
- exact particle tags and grid-frame positions across the comparisons;
- fluid, particle thermodynamic, and controller differences at approximately
  `1e-11` or smaller;
- 96 unique test tags with no repeated seed events after restart.

The principal reproducible local commands are:

```bash
cd tst
python run_test_suite.py --cpu \
  --test test_suite/nr/test_nr_trml_tracers_tracking_cpu.py
python run_test_suite.py --mpicpu \
  --test test_suite/nr/test_nr_frame_tracking_amr_mpicpu.py
```

## Frontier preflight

Start from an unmodified, published integration commit and record its identity
outside the checkout:

```bash
git fetch origin
git switch TRML_with_Tracers_and_Tracking
git pull --ff-only
git submodule update --init
test -z "$(git status --porcelain)"

export TRML_HANDOFF="$(cd .. && pwd -P)/trml-frontier-handoff-$(date +%Y%m%d_%H%M%S)"
mkdir -p "${TRML_HANDOFF}"
git rev-parse HEAD | tee "${TRML_HANDOFF}/tested-commit.txt"
git submodule status | tee "${TRML_HANDOFF}/tested-submodules.txt"
```

The module tuple below matches the repository's June 2026 Frontier build
pattern. Treat the four versioned modules as one compatible tuple. If OLCF has
retired it, replace all four with a currently supported tuple rather than
mixing versions.

```bash
module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0
module load cmake
module unload darshan-runtime

module -t list 2>&1 | tee "${TRML_HANDOFF}/build-modules.txt"
command -v CC
command -v hipcc

: "${CRAY_LD_LIBRARY_PATH:?missing Frontier library path}"
: "${ROCM_PATH:?missing ROCm path}"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export FI_MR_CACHE_MONITOR=kdreg2
```

Configure and build the exact combined pgen:

```bash
cmake -S . -B build_trml_frontier \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
  -DPROBLEM=simple_TRML \
  2>&1 | tee "${TRML_HANDOFF}/cmake.log"

cmake --build build_trml_frontier --parallel 16 \
  2>&1 | tee "${TRML_HANDOFF}/build.log"
```

The configure log must report MPI and HIP, the executable must be
`build_trml_frontier/src/athena`, and the build must finish without warnings
that the requested Kokkos architecture was ignored.

## Allocation setup

Use eight ranks per Frontier node and one rank per MI250X GCD. Inside an
allocation, define:

```bash
export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_VERSION_DISPLAY=1
export FI_MR_CACHE_MONITOR=kdreg2

export REPO="$(pwd -P)"
export EXE="${REPO}/build_trml_frontier/src/athena"
export INPUT="${REPO}/inputs/hydro/TRML/TRML_with_Tracers_and_Tracking.athinput"
export TRML_RUN_ROOT="${MEMBERWORK}/${SLURM_JOB_ACCOUNT}/trml-tracers-tracking/${SLURM_JOB_ID}"
mkdir -p "${TRML_RUN_ROOT}"
```

All launches use `--cpus-per-task=7 --threads-per-core=1 --gpus-per-task=1
--gpu-bind=closest`. Keep restart runs at the same rank count. With
`single_file_per_rank=true`, pass rank 0's restart path; AthenaK maps the
companion rank files internally.

## Frontier gate F0: one-node GPU smoke

Run this first in a one-node allocation. It retains all four integrated
features while reducing particle count and limiting the run to eight steps.

```bash
export F0="${TRML_RUN_ROOT}/F0-smoke"
mkdir -p "${F0}"
srun -N 1 -n 8 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -i "${INPUT}" -d "${F0}" \
  job/basename=F0 \
  time/nlim=8 time/tlim=1.0 \
  frame_tracking/diagnostic_every=1000 \
  tracer_seed1/count_per_event=48 \
  tracer_seed2/end_time=0.01 tracer_seed2/cadence=0.005 \
  tracer_seed2/count_per_event=16 \
  output1/dt=1.0e-20 output2/dt=1.0e-20 \
  output3/variable=hydro_u output3/dt=1.0e-20 \
  output4/dt=10 output5/dt=1.0e-20 \
  2>&1 | tee "${F0}/launch.log"
```

F0 passes only if:

- the run reaches cycle 8 without a HIP, MPI, Kokkos, or finite-value failure;
- `F0.frame_tracker.hst` has positive `ft_weight`, zero `ft_misses`, and
  nonzero `ft_dv_x3` and `ft_dx_x3` after controller actuation;
- the particle history contains 96 unique tags, seed IDs 1 and 2, finite
  coordinates and thermodynamic fields, and nonnegative MeshBlock IDs;
- the final conserved binary contains finite `dens`, `mom1`, `mom2`, `mom3`,
  `ener`, and `r_00` fields;
- all eight `rst/rank_*` directories contain the same final restart index.

Inspect the particle stream with:

```bash
python "${REPO}/scripts/read_prtcl_thermo_history.py" \
  "${F0}/prtcl_thermo_history/F0.thermo.thp" \
  --npz "${F0}/F0-thermo.npz"
```

Do not proceed if F0 fails.

## Frontier gate F1: uninterrupted versus restart-split

Use one node and the canonical 64 x 64 x 128 mesh. Run an uninterrupted
32-cycle reference, then a 16 + 16 cycle restart using the same eight ranks.

```bash
export CONTINUOUS="${TRML_RUN_ROOT}/F1-continuous"
export SPLIT="${TRML_RUN_ROOT}/F1-split"
mkdir -p "${CONTINUOUS}" "${SPLIT}"

common_f1_overrides=(
  time/tlim=1.0
  frame_tracking/diagnostic_every=1000
  tracer_seed2/end_time=0.06 tracer_seed2/cadence=0.0012
  output1/dt=1.0e-20 output2/dt=1.0e-20
  output3/variable=hydro_u output3/dt=1.0e-20
  output4/dt=10 output5/dt=1.0e-20
)

srun -N 1 -n 8 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -i "${INPUT}" -d "${CONTINUOUS}" \
  job/basename=F1Continuous time/nlim=32 "${common_f1_overrides[@]}" \
  2>&1 | tee "${CONTINUOUS}/launch.log"

srun -N 1 -n 8 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -i "${INPUT}" -d "${SPLIT}" \
  job/basename=F1Split time/nlim=16 "${common_f1_overrides[@]}" \
  2>&1 | tee "${SPLIT}/initial-launch.log"

restart_file="$(find "${SPLIT}/rst/rank_00000000" -name 'F1Split.*.rst' \
  -print | sort | tail -n 1)"
test -n "${restart_file}"

srun -N 1 -n 8 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -r "${restart_file}" -i "${INPUT}" -d "${SPLIT}" \
  job/basename=F1Split time/nlim=32 "${common_f1_overrides[@]}" \
  2>&1 | tee "${SPLIT}/restart-launch.log"
```

Do not change the particle variable list, precision, or rank count between the
two halves. Compare final conserved state and controller state with:

```bash
python "${REPO}/scripts/compare_frame_tracking_validation.py" \
  --reference-dir "${CONTINUOUS}" \
  --candidate-dir "${SPLIT}" \
  --binary-id state \
  --output "${TRML_RUN_ROOT}/F1-restart.csv" \
  --problem simple_TRML \
  --resolution 64x64x128 \
  --tracking-mode scalar0-tracer-mass-centroid \
  --restart-split 16+16 \
  --ranks 8 \
  --amr-mode uniform \
  --comparison-reference uninterrupted \
  --comparison-kind restart \
  --axis x3 \
  --selection scalar \
  --scalar-field r_00 \
  --state conserved \
  --require-real-output \
  --target-min 0 --target-max 1 \
  --tolerance 1.0e-10
```

F1 passes if every CSV row is `pass`, particle tags and seed IDs match exactly,
the final tag set has 3,072 members, all 51 top-injection events are represented,
and no seed event is duplicated across the restart boundary.

## Frontier gate F2: MPI plus adaptive refinement

After F1 passes, use one node and enable the committed scalar refinement
criterion:

```bash
export F2="${TRML_RUN_ROOT}/F2-amr"
mkdir -p "${F2}"
srun -N 1 -n 8 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -i "${INPUT}" -d "${F2}" \
  job/basename=F2 \
  mesh_refinement/refinement=adaptive \
  time/nlim=32 time/tlim=1.0 \
  tracer_seed2/end_time=0.03 tracer_seed2/cadence=0.0006 \
  output3/variable=hydro_u output3/dt=1.0e-20 \
  output4/dt=10 output5/dt=1.0e-20 \
  2>&1 | tee "${F2}/launch.log"
```

F2 passes if refinement actually increases the MeshBlock count, all 3,072 tags
remain unique and assigned to valid MeshBlocks, frame history remains finite
with zero missed samples, and a same-rank restart can advance the refined state
for at least eight additional cycles without reseeding particles.

## Frontier gate F3: larger multi-node run

Only after F0-F2 pass, request four nodes and use a 256 x 256 x 512 mesh with
32 x 32 x 64 MeshBlocks. This produces 512 root MeshBlocks, or 16 per GPU
across 32 ranks:

```bash
export F3="${TRML_RUN_ROOT}/F3-4node"
mkdir -p "${F3}"
srun -N 4 -n 32 --ntasks-per-node=8 --cpus-per-task=7 \
  --threads-per-core=1 --gpus-per-task=1 --gpu-bind=closest --unbuffered \
  "${EXE}" -i "${INPUT}" -d "${F3}" \
  job/basename=F3 \
  mesh/nx1=256 mesh/nx2=256 mesh/nx3=512 \
  meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=64 \
  time/nlim=100 time/tlim=1.0 \
  output4/dt=10 \
  2>&1 | tee "${F3}/launch.log"
```

First run F3 with uniform refinement. An adaptive four-node campaign is a
separate performance/science decision after the uniform run establishes memory
headroom and stable step time.

Record the commit, module list, allocation, complete launch line, wall time,
maximum MeshBlock count, restart inventory, and all validation artifacts with
the run. A passing F3 establishes platform readiness; it does not by itself
establish physical convergence or a production resolution.

## Stop conditions

Stop at the first failed gate. In particular, do not spend a multi-node
allocation after any of these observations:

- GPU compilation failure or a host-only Kokkos build;
- any non-finite fluid, frame, or particle value;
- missing or duplicate particle tags;
- a nonzero `ft_misses` count;
- particle history schema mismatch on restart;
- restart disagreement above the F1 tolerance;
- AMR particles with invalid MeshBlock IDs;
- incomplete per-rank restart sets.

Those failures are integration evidence and should be diagnosed before tuning
the physical setup or increasing resolution.
