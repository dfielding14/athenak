#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module unload darshan-runtime rocm craype-accel-amd-gfx90a
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export TMPDIR="$PWD/final-research/release-build-cpu/tmp"
export MPICH_GPU_SUPPORT_ENABLED=0 PYTHONDONTWRITEBYTECODE=1
module -t list > final-research/release-cpu-link-modules.txt 2>&1
python3 - <<'PY'
from pathlib import Path
import subprocess,shlex,hashlib,json,os,datetime
root=Path.cwd();build=root/'final-research/release-build-cpu';out=root/'final-research/release-bin/athena-cpu';out.parent.mkdir(exist_ok=True);assert not out.exists()
def hashes():
 return {str(p.relative_to(build)):hashlib.sha256(p.read_bytes()).hexdigest() for pat in ('*.o','*.a') for p in sorted(build.rglob(pat))}
before=hashes();original=build/'src/athena';oldsha=hashlib.sha256(original.read_bytes()).hexdigest()
cmd=shlex.split((build/'src/CMakeFiles/athena.dir/link.txt').read_text());assert cmd[cmd.index('-o')+1]=='athena';cmd[cmd.index('-o')+1]=str(out)
subprocess.run(cmd,cwd=build/'src',check=True)
after=hashes();assert before==after;assert hashlib.sha256(original.read_bytes()).hexdigest()==oldsha
needed=subprocess.check_output(['readelf','-d',str(out)],text=True)
assert 'libamdhip' not in needed and 'libmpi_gtl_hsa' not in needed
record={'created':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_binary':str(original),'original_sha256':oldsha,'binary':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'command':cmd,'cwd':str(build/'src'),'module_log':'final-research/release-cpu-link-modules.txt','object_and_archive_hashes':before,'objects_and_original_binary_unchanged':True,'dynamic_section':needed,'reason':'Original release CPU link inherited ROCm library injection; relink identical CPU objects under CPU-only modules.'}
(root/'final-research/release-cpu-link-manifest.json').write_text(json.dumps(record,indent=2)+'\n');print(record['binary'],record['sha256'])
PY
