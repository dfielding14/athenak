#!/bin/bash
# Build and run a compact 3D finite-rigidity Bell saturation experiment.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-bell-saturation
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
pic_root="${PIC_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC}"
deck="${repo_root}/inputs/publication/pic_bell_nonlinear_saturation_engineering_v1.athinput"
basename="bell_nonlinear_saturation_engineering_v1"
output_root="${PIC_OUTPUT_ROOT:-${pic_root}/engineering_readiness/bell-saturation-${SLURM_JOB_ID}}"
build_root="${output_root}/build"

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"
export OMP_NUM_THREADS=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0

mkdir -p "${output_root}" "${build_root}"
git -C "${repo_root}" status --short >"${output_root}/source.status"
git -C "${repo_root}" rev-parse HEAD >"${output_root}/source.commit"

/usr/bin/time -p -o "${output_root}/configure.time" \
  /usr/bin/cmake -S "${repo_root}" -B "${build_root}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DAthena_SINGLE_PRECISION=OFF \
    -DAthena_ENABLE_MPI=ON \
    -DKokkos_ENABLE_HIP=ON \
    -DKokkos_ARCH_ZEN3=ON \
    -DKokkos_ARCH_AMD_GFX90A=ON \
    -DCMAKE_CXX_COMPILER=/opt/cray/pe/craype/2.7.33/bin/CC \
    -DCMAKE_CXX_FLAGS=-I/opt/rocm-6.2.4/include \
    '-DCMAKE_EXE_LINKER_FLAGS=-L/opt/rocm-6.2.4/lib -lamdhip64' \
    -DPROBLEM=built_in_pgens \
    >"${output_root}/configure.stdout" 2>"${output_root}/configure.stderr"

/usr/bin/time -p -o "${output_root}/build.time" \
  /usr/bin/cmake --build "${build_root}" --parallel 32 \
    >"${output_root}/build.stdout" 2>"${output_root}/build.stderr"

if [[ -x "${build_root}/src/athena" ]]; then
  executable="${build_root}/src/athena"
elif [[ -x "${build_root}/athena" ]]; then
  executable="${build_root}/athena"
else
  printf 'Frontier build produced no Athena executable\n' >&2
  exit 1
fi

sha256sum "${executable}" "${deck}" \
  "${repo_root}/tst/publication/bell_saturation_engineering_analysis_v1.py" \
  >"${output_root}/bindings.sha256"

/usr/bin/time -p -o "${output_root}/run.time" \
  srun -N1 -n8 --ntasks-per-node=8 --cpus-per-task=1 \
    --gpus-per-task=1 --gpu-bind=closest \
    "${executable}" -i "${deck}" -d "${output_root}" \
    >"${output_root}/athena.stdout" 2>"${output_root}/athena.stderr"

cd "${repo_root}"
/opt/cray/pe/python/3.11.7/bin/python3 -B -m \
  tst.publication.bell_saturation_engineering_analysis_v1 \
  "${output_root}" "${basename}" --output "${output_root}/saturation_analysis.json" \
  --require-saturation

find "${output_root}/bin" -type f -print0 | sort -z | \
  xargs -0 sha256sum >"${output_root}/raw_outputs.sha256"
printf '%s\n' 'execution_complete=true' 'sustained_saturation_gate_passed=true' \
  >"${output_root}/engineering_readiness.status"
