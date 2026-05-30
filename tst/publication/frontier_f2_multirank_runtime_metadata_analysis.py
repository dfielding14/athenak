#!/usr/bin/env python3
"""Validate the bounded Frontier F2 eight-rank runtime-metadata slice."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from frontier_f0_structured_smoke_analysis import _parse_allowlist
from frontier_f0_structured_smoke_analysis import _read_regular
from frontier_f0_structured_smoke_analysis import _sha256
from frontier_f0_structured_smoke_analysis import _write_result


EXPECTED_RANKS = set(range(8))
EXPECTED_DEVICES = set(range(8))
GPU_PREFLIGHT = re.compile(
    r"^PIC trusted GPU launch: rank=(?P<rank>[0-9]+) "
    r"host=(?P<host>\S+) ROCR_VISIBLE_DEVICES=(?P<device>[0-9]+) "
    r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$"
)
RUNTIME_MODEL_PREFIX = "PIC runtime model: physical_mode=extended_mhd_pic "
EXPECTED_RUNTIME_TOKENS = {
    "background=coupled",
    "feedback=coupled",
    "induction=ideal_mhd_only",
    "deposition=tsc",
    "deltaf=off",
    "expanding_box=off",
    "wave_damping=off",
    "restart_schema=7",
}


def _decode(path: Path, *, require_nonempty: bool) -> tuple[bytes, str]:
    data = _read_regular(path, require_nonempty=require_nonempty)
    try:
        return data, data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"F2 runtime artifact is not UTF-8: {path}") from error


def _rank_gpu_bindings(stdout: str) -> list[dict[str, object]]:
    bindings = []
    for line in stdout.splitlines():
        match = GPU_PREFLIGHT.fullmatch(line)
        if match is not None:
            bindings.append(
                {
                    "rank": int(match.group("rank")),
                    "host": match.group("host"),
                    "rocr_visible_device": int(match.group("device")),
                }
            )
    ranks = {binding["rank"] for binding in bindings}
    devices = {binding["rocr_visible_device"] for binding in bindings}
    hosts = {binding["host"] for binding in bindings}
    if len(bindings) != len(EXPECTED_RANKS) or ranks != EXPECTED_RANKS:
        raise ValueError("F2 trusted GPU preflight does not contain exactly ranks 0 through 7")
    if devices != EXPECTED_DEVICES:
        raise ValueError("F2 trusted GPU preflight does not bind exactly devices 0 through 7")
    if len(hosts) != 1:
        raise ValueError("F2 trusted GPU preflight does not remain on one Frontier node")
    return sorted(bindings, key=lambda binding: int(binding["rank"]))


def _require_runtime_model(stdout: str) -> str:
    models = [
        line for line in stdout.splitlines() if line.startswith(RUNTIME_MODEL_PREFIX)
    ]
    if len(models) != 1:
        raise ValueError("F2 Athena stdout does not contain exactly one runtime-model line")
    tokens = set(models[0].split())
    missing = EXPECTED_RUNTIME_TOKENS - tokens
    if missing:
        raise ValueError(f"F2 Athena runtime model is missing expected tokens: {sorted(missing)}")
    return models[0]


def analyze(artifact_dir: Path) -> dict[str, object]:
    allowlist_data, allowlist = _parse_allowlist(
        artifact_dir / "f2-runtime-metadata.environment.allowlist.txt"
    )
    stdout_data, stdout = _decode(
        artifact_dir / "athena_stdout.txt", require_nonempty=True
    )
    stderr_data, _ = _decode(
        artifact_dir / "athena_stderr.txt", require_nonempty=False
    )
    if "Number of parallel ranks = 8" not in stdout.splitlines():
        raise ValueError("F2 Athena stdout does not report eight parallel ranks")
    bindings = _rank_gpu_bindings(stdout)
    runtime_model = _require_runtime_model(stdout)
    return {
        "schema_version": 1,
        "status": "pass",
        "evidence_class": "frontier_f2_clean_candidate_multirank_runtime_metadata",
        "launch_contract": "trusted_trampoline_athena_argv_v1",
        "parser_runtime_smoke": "pass",
        "runtime_profile": allowlist["PIC_FRONTIER_PROFILE"],
        "parallel_ranks": 8,
        "hosts": sorted({str(binding["host"]) for binding in bindings}),
        "rank_gpu_bindings": bindings,
        "runtime_model": runtime_model,
        "runtime_artifacts": {
            "f2-runtime-metadata.environment.allowlist.txt": _sha256(allowlist_data),
            "athena_stdout.txt": _sha256(stdout_data),
            "athena_stderr.txt": _sha256(stderr_data),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    _write_result(args.artifact_dir / "analysis.json", analyze(args.artifact_dir))


if __name__ == "__main__":
    main()
