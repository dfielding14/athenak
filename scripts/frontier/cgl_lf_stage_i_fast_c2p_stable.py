#!/usr/bin/env python3
"""Launch fresh R14/R15 C2P and heat-flux stability successors on Frontier."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


SCRIPT_PATH = Path(__file__).resolve()
CORRECTED_LAUNCHER = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_corrected.py")
SUPPORTED_CASES = ("R14", "R15")
IDENTITY_SCHEMA = "athenak-cgl-c2p-stable-campaign-identity"
HISTORICAL_IDENTITY = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1/campaign-identity.json"
)
HISTORICAL_IDENTITY_SHA256 = (
    "fad99d9651e6da155bb5dfae1b81a8675cce232fb17e31a36e9be36d62810643"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def requested_root(arguments: list[str]) -> Path:
    for index, argument in enumerate(arguments):
        if argument == "--root":
            if index + 1 >= len(arguments):
                raise RuntimeError("--root requires a value")
            return Path(arguments[index + 1]).expanduser().resolve()
        if argument.startswith("--root="):
            return Path(argument.split("=", 1)[1]).expanduser().resolve()
    raise RuntimeError("--root is required for the C2P-stability successor launcher")


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"expected a JSON object: {path}")
    return value


def load_key_values(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key] = value
    return result


def require_binding(
    identity: dict[str, object], label: str
) -> tuple[Path, str]:
    artifacts = identity.get("artifacts")
    if not isinstance(artifacts, dict):
        raise RuntimeError("campaign identity lacks artifacts")
    value = artifacts.get(label)
    if not isinstance(value, dict):
        raise RuntimeError(f"campaign identity lacks {label} binding")
    raw_path = value.get("path")
    expected_sha = value.get("sha256")
    if not isinstance(raw_path, str) or not isinstance(expected_sha, str):
        raise RuntimeError(f"campaign identity has invalid {label} binding")
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise RuntimeError(f"campaign identity {label} does not exist: {path}")
    actual_sha = sha256_file(path)
    if actual_sha != expected_sha:
        raise RuntimeError(
            f"campaign identity {label} SHA differs: {actual_sha} != {expected_sha}"
        )
    return path, expected_sha


def load_identity(root: Path) -> tuple[Path, dict[str, object]]:
    identity_path = root / "campaign-identity.json"
    identity = load_json(identity_path)
    if (
        identity.get("schema") != IDENTITY_SCHEMA
        or identity.get("schema_version") != 1
    ):
        raise RuntimeError("campaign identity schema differs")
    if Path(str(identity.get("campaign_root"))).resolve() != root:
        raise RuntimeError("campaign identity root differs")
    if identity.get("selected_cases") != list(SUPPORTED_CASES):
        raise RuntimeError("campaign identity must select exactly R14 and R15")
    return identity_path, identity


def load_corrected_launcher():
    name = "_cgl_lf_stage_i_fast_c2p_stable_base"
    spec = importlib.util.spec_from_file_location(name, CORRECTED_LAUNCHER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import corrected launcher: {CORRECTED_LAUNCHER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def restrict_validate_provenance(corrected) -> None:
    base_validate_provenance = corrected.validate_provenance

    def successor_validate_provenance(
        root: Path, cases=SUPPORTED_CASES
    ) -> dict[str, object]:
        return base_validate_provenance(root, cases)

    corrected.validate_provenance = successor_validate_provenance


def configure_successor(root: Path):
    identity_path, identity = load_identity(root)
    source_path, _ = require_binding(identity, "source_identity")
    historical_path, historical_sha = require_binding(
        identity, "historical_campaign_identity"
    )
    build_environment_path, _ = require_binding(identity, "build_environment")
    executable_path, executable_sha = require_binding(identity, "executable")
    matrix_path, matrix_sha = require_binding(identity, "matrix")
    r14_input, r14_input_sha = require_binding(identity, "R14_input")
    r15_input, r15_input_sha = require_binding(identity, "R15_input")

    source = Path(str(identity.get("source"))).expanduser().resolve()
    if source_path != source / "source-identity.json":
        raise RuntimeError("source-tree binding is not rooted in the declared source")
    source_identity = load_json(source_path)
    if (
        source_identity.get("schema") != "athenak-cgl-source-purpose"
        or source_identity.get("schema_version") != 1
        or source_identity.get("repair") != "cgl-c2p-and-heat-flux-stability"
    ):
        raise RuntimeError("source-purpose identity differs")
    if (
        historical_path != HISTORICAL_IDENTITY
        or historical_sha != HISTORICAL_IDENTITY_SHA256
    ):
        raise RuntimeError("historical campaign identity binding differs")
    historical_identity = load_json(historical_path)
    if historical_identity.get("campaign_id") != (
        "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
    ):
        raise RuntimeError("historical campaign ID differs")
    if SCRIPT_PATH.parents[2] != source:
        raise RuntimeError("launcher is not executing from the declared frozen source")
    if matrix_path.parent.parent.parent != source:
        raise RuntimeError("matrix is not rooted in the declared frozen source")
    if r14_input.parent.parent.parent != source or r15_input.parent.parent.parent != source:
        raise RuntimeError("case inputs are not rooted in the declared frozen source")
    source_revision = identity.get("source_revision")
    campaign_id = identity.get("campaign_id")
    if not isinstance(source_revision, str) or len(source_revision) != 40:
        raise RuntimeError("campaign identity has invalid source revision")
    if not isinstance(campaign_id, str) or not campaign_id:
        raise RuntimeError("campaign identity has invalid campaign ID")
    if root.name != campaign_id:
        raise RuntimeError("campaign identity ID differs from the root directory")
    build_environment = load_key_values(build_environment_path)
    if build_environment.get("git_revision") != source_revision:
        raise RuntimeError("build environment source revision differs")
    if Path(build_environment.get("source_dir", "")).resolve() != source:
        raise RuntimeError("build environment source directory differs")
    if (
        Path(build_environment.get("build_dir", "")).resolve()
        != executable_path.parent.parent
    ):
        raise RuntimeError("build environment executable directory differs")

    corrected = load_corrected_launcher()
    corrected.SCRIPT_PATH = SCRIPT_PATH
    corrected.CAMPAIGN_ID = campaign_id
    corrected.CAMPAIGN_ROOT = root
    corrected.FROZEN_SOURCE = source
    corrected.SOURCE_REVISION = source_revision
    corrected.MATRIX_RELATIVE = matrix_path.relative_to(source)
    corrected.MATRIX_SHA256 = matrix_sha
    corrected.ATHENA = executable_path
    corrected.ATHENA_SHA256 = executable_sha
    corrected.IDENTITY = identity_path
    corrected.IDENTITY_SHA256 = sha256_file(identity_path)
    corrected.ACTIVE_CASES = SUPPORTED_CASES
    corrected.FINITE_LIMITER_DIAGNOSTIC_CASES = frozenset(SUPPORTED_CASES)
    corrected.FINITE_LIMITER_VARIANT = (
        "finite_limiter_c2p_heat_flux_stable_successor_nonfatal"
    )
    corrected.DEFAULT_CASE_NODES = {"R14": 24, "R15": 24}
    corrected.MAX_USEFUL_NODES = {"R14": 27, "R15": 27}
    corrected.CASE_INPUT_SHA256 = {
        "R14": r14_input_sha,
        "R15": r15_input_sha,
    }
    restrict_validate_provenance(corrected)

    base_batch_script_text = corrected.batch_script_text

    def successor_batch_script_text(manifest: dict[str, object]) -> str:
        text = base_batch_script_text(manifest)
        expected = "#SBATCH -J cglc_"
        if text.count(expected) != 1:
            raise RuntimeError("corrected launcher job-name prefix differs")
        return text.replace(expected, "#SBATCH -J cgls_", 1)

    base_prepare_segment = corrected.prepare_segment

    def successor_prepare_segment(*args, **kwargs) -> Path:
        segment = base_prepare_segment(*args, **kwargs)
        manifest_path = corrected.fast.segment_manifest(segment)
        manifest = corrected.fast.load_json(manifest_path)
        manifest.update(
            {
                "lineage_relationship": "independent_fresh_successor",
                "successor_launch_origin": "fresh_t0_c2p_stable_successor"
                if int(manifest["sequence"]) == 0
                else "c2p_stable_successor_continuation",
                "historical_campaign_identity": str(historical_path),
                "historical_campaign_identity_sha256": historical_sha,
                "repair": "cgl-c2p-and-heat-flux-stability",
            }
        )
        corrected.fast.write_json(manifest_path, manifest)
        return segment

    base_validate_segment = corrected.validate_segment

    def successor_validate_segment(segment: Path) -> dict[str, object]:
        manifest = base_validate_segment(segment)
        sequence = int(manifest["sequence"])
        expected = {
            "lineage_relationship": "independent_fresh_successor",
            "successor_launch_origin": "fresh_t0_c2p_stable_successor"
            if sequence == 0
            else "c2p_stable_successor_continuation",
            "historical_campaign_identity": str(historical_path),
            "historical_campaign_identity_sha256": historical_sha,
            "repair": "cgl-c2p-and-heat-flux-stability",
        }
        for key, value in expected.items():
            if manifest.get(key) != value:
                raise RuntimeError(f"successor segment {key} differs: {segment}")
        return manifest

    corrected.batch_script_text = successor_batch_script_text
    corrected.prepare_segment = successor_prepare_segment
    corrected.validate_segment = successor_validate_segment
    corrected.configure(root)
    corrected.validate_common_provenance(root)
    return corrected


def main() -> int:
    try:
        root = requested_root(sys.argv[1:])
        corrected = configure_successor(root)
        return corrected.main()
    except (OSError, RuntimeError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
