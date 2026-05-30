#!/usr/bin/env python3
"""Validate artifacts from the registered Frontier F0 HIP/MPI smoke."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


MAPPING_PATTERN = re.compile(
    r"^rank=(?P<rank>\d+) host=(?P<host>\S+) "
    r"ROCR_VISIBLE_DEVICES=(?P<rocr>\S+) "
    r"GPU_DEVICE_ORDINAL=(?P<ordinal>\S+) "
    r"HIP_VISIBLE_DEVICES=(?P<hip>\S+)$"
)


def require_nonempty(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"Missing F0 artifact: {path}")
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"Empty F0 artifact: {path}")
    return text


def analyze(artifact_dir: Path) -> dict[str, object]:
    require_nonempty(artifact_dir / "snapshot_verification.txt")
    environment_text = require_nonempty(artifact_dir / "environment.allowlist.txt")
    require_nonempty(artifact_dir / "modules.txt")
    require_nonempty(artifact_dir / "rocm_smi.txt")
    require_nonempty(artifact_dir / "athena_parser_stdout.txt")
    mapping_text = require_nonempty(artifact_dir / "gpu_mapping.txt")

    mappings = []
    for line in mapping_text.splitlines():
        match = MAPPING_PATTERN.match(line)
        if match is None:
            raise ValueError(f"Malformed GPU mapping line: {line}")
        mappings.append(match.groupdict())
    ranks = sorted(int(mapping["rank"]) for mapping in mappings)
    if ranks != list(range(8)):
        raise ValueError(f"Expected GPU mapping ranks 0..7, found {ranks}")
    if any(mapping["rocr"] == "unset" for mapping in mappings):
        raise ValueError("ROCR_VISIBLE_DEVICES was not assigned for every MPI rank")
    visible_devices = {mapping["rocr"] for mapping in mappings}
    if len(visible_devices) != 8:
        raise ValueError(
            f"Expected eight distinct ROCR_VISIBLE_DEVICES bindings, found {visible_devices}"
        )

    environment = dict(
        line.split("=", 1)
        for line in environment_text.splitlines()
        if "=" in line
    )
    if environment.get("MPICH_GPU_SUPPORT_ENABLED") != "1":
        raise ValueError("MPICH_GPU_SUPPORT_ENABLED must be 1")
    if environment.get("SLURM_EXPORT_ENV") != "ALL":
        raise ValueError("SLURM_EXPORT_ENV must be ALL for compute steps")

    return {
        "schema_version": 1,
        "status": "pass",
        "gpu_mapping_rank_count": len(mappings),
        "gpu_mapping_ranks": ranks,
        "parser_smoke": "pass",
        "snapshot_verification": "pass",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    result = analyze(args.artifact_dir)
    (args.artifact_dir / "analysis.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
