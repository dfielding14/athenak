#!/usr/bin/env python3
"""Retain and launch one member/segment of the matched CGL-LF experiment.

Run inside an allocation, or supply --job-id for an existing allocation.
The executable's revision comes from its retained build manifest, not HEAD.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time


SOURCE = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4*1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def override_input(text, overrides):
    for item in overrides:
        lhs, value = item.split("=", 1)
        section, key = lhs.split("/", 1)
        match = re.search(rf"(?ms)(^<{re.escape(section)}>\s*\n)(.*?)(?=^<|\Z)", text)
        if match is None:
            raise ValueError("Unknown override section: "+section)
        body = match.group(2)
        expression = rf"(?m)^{re.escape(key)}\s*=.*$"
        body = (re.sub(expression, lambda _: key+" = "+value, body)
                if re.search(expression, body) else body+key+" = "+value+"\n")
        text = text[:match.start()]+match.group(1)+body+text[match.end():]
    return text


def parameters(text):
    result, section = {}, None
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            section = line[1:-1]
        elif section and "=" in line:
            key, value = line.split("=", 1)
            result[section+"/"+key.strip()] = value.strip()
    return result


def same_value(a, b):
    if a is None or b is None:
        return a == b
    booleans = {"true": "true", "1": "true", "false": "false", "0": "false"}
    if a.lower() in ("true", "false") or b.lower() in ("true", "false"):
        return booleans.get(a.lower(), a) == booleans.get(b.lower(), b)
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-14, abs_tol=0)
    except ValueError:
        return a == b


def restart_input(path, intended, overrides):
    """Check physical/numerical identity before allowing a checkpoint resume."""
    with Path(path).open("rb") as stream:
        raw = stream.read(4*1024*1024)
    marker = b"<par_end>"
    if marker not in raw:
        raise ValueError("Restart has no parameter header within the first 4 MiB")
    header = raw.split(marker, 1)[0].decode("utf-8")
    old, wanted = parameters(header), parameters(intended)
    administrative = {"time/tlim", "time/nlim", "time/ndiag"}
    mismatches = []
    for key, value in wanted.items():
        if key.startswith(("output", "job/")) or key in administrative:
            continue
        if not same_value(old.get(key), value):
            mismatches.append(f"{key}: checkpoint={old.get(key)!r}, requested={value!r}")
    if mismatches:
        raise ValueError("Checkpoint configuration does not match this member:\n"+"\n".join(mismatches))
    # The code reads the restart header plus CLI overrides, not the proposed deck.
    return override_input(header, overrides)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--mode", choices=("active", "passive"), required=True)
    parser.add_argument("--input", type=Path, default=SOURCE/"inputs/cgl_lf_paper/cgl_lf_physics_benchmark_matched_beta10.athinput")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--build-manifest", type=Path, required=True)
    parser.add_argument("--build-cache", type=Path, required=True)
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--ranks-per-node", type=int, default=8)
    parser.add_argument("--job-id", default=os.environ.get("SLURM_JOB_ID"))
    parser.add_argument("--nodelist", help="Optional Slurm node list within the allocation")
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--wall-time", help="AthenaK clean-checkpoint wall limit HH:MM:SS")
    parser.add_argument("--classification", default="matched scientific experiment")
    parser.add_argument("--set", action="append", default=[], dest="overrides")
    args = parser.parse_args()
    if args.nodes < 1 or args.ranks_per_node < 1 or not args.job_id:
        parser.error("positive node/rank counts and an allocation job ID are required")
    run, binary = args.run_dir.resolve(), args.executable.resolve()
    manifest = json.loads(args.build_manifest.read_text())
    revision = manifest.get("revision")
    if not revision or manifest.get("binary_sha256") != sha(binary):
        raise ValueError("Build manifest must identify the executable revision and matching SHA256")
    cache = dict(re.findall(r"(?m)^([A-Za-z0-9_]+):[^=\n]+=(.*)$", args.build_cache.read_text()))
    if cache.get("PROBLEM") != "built_in_pgens":
        raise ValueError("This benchmark requires a retained PROBLEM=built_in_pgens build")
    for key, value in (("Kokkos_ENABLE_HIP", "ON"), ("Athena_ENABLE_MPI", "ON"),
                       ("Athena_SINGLE_PRECISION", "OFF")):
        if cache.get(key) != value:
            raise ValueError(f"This GPU launcher requires {key}={value} in the retained build cache")
    if run.exists() and any(run.iterdir()):
        raise ValueError("Use a new empty segment directory; existing data are never overwritten")
    # The numerical source must still match this executable's recorded revision.
    # Custom unit-test pgens are not compiled into PROBLEM=built_in_pgens.
    for diff in (["diff", revision, "--"], ["diff", "--cached", "--"]):
        changed = subprocess.check_output(["git", "-C", str(SOURCE), *diff,
            "src", "CMakeLists.txt", "cmake", ":(exclude)src/pgen/unit_tests"], text=True)
        if changed:
            raise ValueError("Numerical source differs from the retained executable; build and record it first")
    passive = "true" if args.mode == "passive" else "false"
    for item in args.overrides:
        key, value = item.split("=", 1)
        if key in ("mhd/passive", "problem/passive_delta") and not same_value(value, passive):
            raise ValueError("--set conflicts with the selected --mode: "+item)
    overrides = ["mhd/passive="+passive, "problem/passive_delta="+passive, *args.overrides]
    effective = override_input(args.input.read_text(), overrides)
    if args.restart:
        effective = restart_input(args.restart, effective, overrides)
    run.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.input, run/"input.athinput")
    (run/"effective.athinput").write_text(effective)
    shutil.copy2(args.build_manifest, run/"binary_manifest.json")
    shutil.copy2(args.build_cache, run/"CMakeCache.txt")
    shutil.copy2(__file__, run/"launch-script.py")
    launch = ["srun", "--jobid", str(args.job_id), "--exact", "-N", str(args.nodes),
        "-n", str(args.nodes*args.ranks_per_node), "--ntasks-per-node", str(args.ranks_per_node),
        "--threads-per-core=1", "--cpu-bind=threads", "-c", "7", "--gpus-per-task=1", "--gpu-bind=closest"]
    if args.nodelist:
        launch.extend(["--nodelist", args.nodelist])
    initial = ["-r", str(args.restart.resolve())] if args.restart else ["-i", "input.athinput"]
    wall = ["-t", args.wall_time] if args.wall_time else []
    command = [*launch, str(binary), *initial, *wall, *overrides]
    record = {"schema_version": 1, "classification": args.classification, "member": args.mode,
        "simulation": {"revision": revision, "dirty": False,
            "revision_basis": "Exact executable SHA from retained build manifest; numerical source comparison passed before launch.",
            "executable": str(binary), "executable_sha256": sha(binary),
            "input_path": "effective.athinput", "input_sha256": sha(run/"effective.athinput"),
            "canonical_input_path": "input.athinput", "canonical_input_sha256": sha(run/"input.athinput"),
            "build_backend": "HIP MPI double; PROBLEM=built_in_pgens; compiler="
                +cache.get("CMAKE_CXX_COMPILER", "unknown")+"; version="
                +cache.get("KOKKOS_CXX_COMPILER_VERSION", "unknown"),
            "build_manifest_path": "binary_manifest.json", "build_manifest_sha256": sha(run/"binary_manifest.json"),
            "build_cache_path": "CMakeCache.txt", "build_cache_sha256": sha(run/"CMakeCache.txt")},
        # Discover histories if cancellation prevents the final metadata save.
        # Finalization below adds exact inventories, which take precedence.
        "outputs": {"snapshot_glob": "**/*.mhd_w_bcc.*.bin", "forcing_glob": "**/*.turb_force.*.bin",
            "user_history_glob": "**/*.user.hst", "mhd_history_glob": "**/*.mhd.hst"},
        "launch": {"command": command, "shell_command": shlex.join(command), "overrides": overrides,
            "cwd": str(run), "nodes": args.nodes, "ranks": args.nodes*args.ranks_per_node,
            "gpus": args.nodes*args.ranks_per_node, "slurm_job_id": str(args.job_id),
            "requested_nodelist": args.nodelist,
            "restart": str(args.restart.resolve()) if args.restart else None,
            "restart_sha256": sha(args.restart) if args.restart else None,
            "source_checkout_revision": subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip(),
            "launcher_sha256": sha(run/"launch-script.py"),
            "environment": {k: v for k, v in os.environ.items() if k.startswith(("MPICH_", "FI_", "HSA_", "OMP_")) or k in ("LOADEDMODULES", "LD_LIBRARY_PATH")},
            "started_utc": datetime.now(timezone.utc).isoformat()},
        "planned_analysis": {"initial_window": [6, 14], "block_duration": 2,
            "procedure": "Use the same interval for both members; inspect four two-unit blocks and extend both to18 if needed. No tuning to a paper slope or desired contrast."}}
    def save():
        (run/"benchmark_metadata.json").write_text(json.dumps(record, indent=2)+"\n")
    save()
    (run/"run-command.txt").write_text(shlex.join(command)+"\n")
    begin = time.monotonic()
    with (run/"run.log").open("w") as stream:
        result = subprocess.run(command, cwd=run, stdout=stream, stderr=subprocess.STDOUT)
    record["launch"].update({"ended_utc": datetime.now(timezone.utc).isoformat(),
        "returncode": result.returncode, "launch_wall_seconds": time.monotonic()-begin,
        "binary_hash_unchanged": sha(binary) == record["simulation"]["executable_sha256"]})
    for key, pattern in (("user_history", "*.user.hst"), ("mhd_history", "*.mhd.hst")):
        record["outputs"][key] = [str(p.relative_to(run)) for p in sorted(run.rglob(pattern))]
    save()
    print(json.dumps(record["launch"], indent=2))
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
