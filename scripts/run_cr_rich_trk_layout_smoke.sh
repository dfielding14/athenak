#!/usr/bin/env bash
# Run a tiny GPU smoke for rich tracked-particle output layouts.
#
#SBATCH -A AST207
#SBATCH -J CR_RICH_TRK
#SBATCH -o /lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing/%x.%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH --qos=debug
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1

set -euo pipefail

REPO=${AMR_CR_RICH_TRK_REPO:-/ccs/home/dfielding/athenak-cr-tracers-followup-architecture}
EXE=${AMR_CR_RICH_TRK_EXE:-${REPO}/build-frontier/src/athena}
INPUT=${AMR_CR_RICH_TRK_INPUT:-${REPO}/inputs/particles/cr_pusher_uniform_b.athinput}
TESTING_ROOT=${AMR_CR_RICH_TRK_TESTING_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing}
STAMP=${AMR_CR_RICH_TRK_STAMP:-$(date +%Y%m%d_%H%M%S)}
RUN_ROOT=${AMR_CR_RICH_TRK_RUN_ROOT:-${TESTING_ROOT}/cr_rich_trk_layout_${STAMP}}
RANKS=${AMR_CR_RICH_TRK_RANKS:-8}

mkdir -p "${RUN_ROOT}"

if command -v module >/dev/null 2>&1; then
  module restore
  module load cpe/24.07 PrgEnv-amd cray-mpich/8.1.30 craype-accel-amd-gfx90a amd/6.2.0 rocm/6.2.0
  module unload darshan-runtime || true
fi

export MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-1}
export MPICH_SMP_SINGLE_COPY_MODE=${MPICH_SMP_SINGLE_COPY_MODE:-NONE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export HSA_ENABLE_DEBUG=0
ulimit -c 0 || true

test -x "${EXE}"
test -f "${INPUT}"

git -C "${REPO}" rev-parse HEAD >"${RUN_ROOT}/source.commit"
git -C "${REPO}" status --short >"${RUN_ROOT}/source.status"
sha256sum "${EXE}" "${INPUT}" "${REPO}/scripts/analyze_cr_pusher_accuracy.py" \
  >"${RUN_ROOT}/bindings.sha256"
env | sort >"${RUN_ROOT}/environment.txt"

run_mode() {
  local mode=$1
  local per_rank=$2
  local per_node=$3
  local run_dir=${RUN_ROOT}/${mode}
  rm -rf "${run_dir}"
  mkdir -p "${run_dir}"

  local command=(
    "${EXE}" -i "${INPUT}" -d "${run_dir}"
    "job/basename=trk_${mode}"
    "meshblock/nx1=8"
    "meshblock/nx2=8"
    "meshblock/nx3=8"
    "time/tlim=0.125"
    "time/nlim=1000"
    "time/ndiag=20"
    "particles/log_performance=false"
    "particles/check_consistency_mode=none"
    "output1/dt=0.125"
    "output1/ncycle=1"
    "output1/buffer_size=0"
    "output1/single_file_per_rank=${per_rank}"
    "output1/single_file_per_node=${per_node}"
    "output1/validate_global_tags=true"
    "output2/dt=-1.0"
  )
  printf '%q ' "${command[@]}" >"${run_dir}/command.txt"
  printf '\n' >>"${run_dir}/command.txt"

  /usr/bin/time -p -o "${run_dir}/run.time" \
    srun -N1 -n"${RANKS}" --ntasks-per-node="${RANKS}" \
      --cpus-per-task=1 --gpus-per-task=1 --gpu-bind=closest \
      --kill-on-bad-exit=1 \
      "${command[@]}" >"${run_dir}/stdout.txt" 2>"${run_dir}/stderr.txt"
}

run_mode shared false false
run_mode rank true false
run_mode node false true

python3 - "${RUN_ROOT}" "${REPO}" "${RANKS}" <<'PY'
import sys
from pathlib import Path

run_root = Path(sys.argv[1])
repo = Path(sys.argv[2])
ranks = int(sys.argv[3])
sys.path.insert(0, str(repo / "scripts"))
from analyze_cr_pusher_accuracy import read_trk_file  # noqa: E402

expected_counts = {
    "shared": 1,
    "rank": ranks,
    "node": 1,
}

summary = []
for mode, expected_file_count in expected_counts.items():
    run_dir = run_root / mode
    if mode == "shared":
        files = sorted((run_dir / "trk").glob("*.trk"))
    elif mode == "rank":
        files = sorted((run_dir / "trk").glob("rank_*/*.trk"))
    else:
        files = sorted((run_dir / "trk").glob("node_*/*.trk"))
    if len(files) != expected_file_count:
        raise SystemExit(
            f"{mode}: expected {expected_file_count} trk files, found {len(files)}")

    frames = []
    for path in files:
        data = path.read_bytes()
        if b"trk_format=rich_v1" not in data or b"nfields=18" not in data:
            raise SystemExit(f"{mode}: {path} is not a rich_v1 nfields=18 file")
        if b"fields=tag,time,x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag" not in data:
            raise SystemExit(f"{mode}: {path} has the wrong fields header")
        frames.extend(read_trk_file(path))

    latest_time = max(frame.time for frame in frames)
    latest = [frame for frame in frames if abs(frame.time - latest_time) < 1.0e-7]
    expected_tracks = max(frame.ntracked for frame in latest)
    keys = []
    for frame in latest:
        keys.extend(frame.particles)
    if len(keys) != expected_tracks:
        raise SystemExit(
            f"{mode}: latest frame has {len(keys)} records, expected {expected_tracks}")
    if len(set(keys)) != len(keys):
        raise SystemExit(f"{mode}: latest frame has duplicate track keys")
    summary.append(f"{mode}: files={len(files)} latest_time={latest_time} records={len(keys)}")

(run_root / "verification.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))
PY

touch "${RUN_ROOT}/.complete"
printf 'Rich TRK layout smoke complete: %s\n' "${RUN_ROOT}"
