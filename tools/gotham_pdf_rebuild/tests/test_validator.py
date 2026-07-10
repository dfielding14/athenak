"""Focused release-gate validator tests."""

from __future__ import annotations

import importlib.util
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "validate_real_8pc", TOOLS_DIR / "validate_real_8pc.py"
)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def test_control_payload_time_uses_exact_manifest_output_time() -> None:
    source_time = 10.52
    output_time = 10.520002963096713
    report = {"checks": [], "failures": []}
    manifest = {"source_time": source_time, "output_time": output_time}
    rebuilt = {
        product_id: {"time": output_time} for product_id in VALIDATOR.EXPECTED_CONTROLS
    }
    archived = {
        product_id: {"time": output_time} for product_id in VALIDATOR.EXPECTED_CONTROLS
    }

    VALIDATOR.add_control_time_checks(report, manifest, rebuilt, archived)

    checks = {check["name"]: check for check in report["checks"]}
    assert checks["rebuilt_source_time_matches_archived_sequence_00062"]["status"] == "pass"
    exact_check = checks["required_rebuilt_control_times_match_manifest"]
    assert exact_check["status"] == "pass"
    assert exact_check["metrics"]["manifest_output_time"] == output_time
    assert exact_check["metrics"]["maximum_absolute_time_error"] == 0.0


def test_manifest_checks_require_exact_output_identity(tmp_path: Path) -> None:
    archived_root = tmp_path / "phase2"
    manifest = {
        "source_sequence": VALIDATOR.EXPECTED_SOURCE_SEQUENCE,
        "source_cycle": VALIDATOR.EXPECTED_SOURCE_CYCLE,
        "shards_available": VALIDATOR.EXPECTED_SHARDS,
        "shards_processed": VALIDATOR.EXPECTED_SHARDS,
        "source_input_dir": str((archived_root / "bin").resolve()),
        "output_number": "99999",
        "output_time": VALIDATOR.EXPECTED_OUTPUT_TIME + 1.0,
        "geometry_to_logical_key_validated": True,
        "expected_domain_volume": VALIDATOR.EXPECTED_DOMAIN_VOLUME,
        "summed_leaf_volume": VALIDATOR.EXPECTED_DOMAIN_VOLUME,
    }
    report = {"checks": [], "failures": []}

    VALIDATOR.add_manifest_checks(report, manifest, archived_root)

    checks = {check["name"]: check for check in report["checks"]}
    assert checks["manifest_output_number_matches_required_identity"]["status"] == "fail"
    assert checks["manifest_output_time_matches_required_identity"]["status"] == "fail"
