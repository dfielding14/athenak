"""Focused tests for the corrected-production direct-fast Stage I launcher."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_corrected.py"
EXPECTED_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
)
EXPECTED_SOURCE = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-corrected-0c406312f"
)
EXPECTED_REVISION = "0c406312fa35d5e1c7041d80333b0ea24b0127ae"
EXPECTED_MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
EXPECTED_ATHENA = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/build/"
    "frontier-hip-0c406312fa35-cpe25.09-cce20-rocm6.4.2/src/athena"
)
EXPECTED_ATHENA_SHA256 = "0f032379c4d829fc86353cc734a716f6ef6d6eda8951a5b79a3a97ad671e9b4b"
ACTIVE_CASES = (
    "R02",
    "R03",
    "R04",
    "R05",
    "R10",
    "R11",
    "R12",
    "R13",
    "R14",
    "R15",
    "R16",
    "R17",
)


def load_launcher():
    name = "cgl_lf_stage_i_fast_corrected_for_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def corrected():
    return load_launcher()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    names = list(columns)
    lines = [
        "# "
        + " ".join(f"[{index}]={name}" for index, name in enumerate(names, start=1))
    ]
    lines.extend(
        " ".join(format(value, ".17g") for value in row)
        for row in zip(*(columns[name] for name in names))
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_restart_group(directory: Path, rank_count: int, time: float) -> Path:
    payload = f"<time>\ntime = {time:.17g}\n<par_end>\n".encode() + b"\0payload"
    rank_zero = directory / "rank_00000000/fixture.00001.rst"
    for rank in range(rank_count):
        path = directory / f"rank_{rank:08d}/fixture.00001.rst"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    return rank_zero


def test_corrected_defaults_are_exact_active_only_and_aggressive(corrected):
    assert corrected.CAMPAIGN_ROOT == EXPECTED_ROOT
    assert corrected.FROZEN_SOURCE == EXPECTED_SOURCE
    assert corrected.SOURCE_REVISION == EXPECTED_REVISION
    assert corrected.MATRIX_SHA256 == EXPECTED_MATRIX_SHA256
    assert corrected.ATHENA == EXPECTED_ATHENA
    assert corrected.ATHENA_SHA256 == EXPECTED_ATHENA_SHA256
    assert corrected.ACTIVE_CASES == ACTIVE_CASES
    assert set(corrected.DEFAULT_CASE_NODES) == set(ACTIVE_CASES)
    assert corrected.DEFAULT_CASE_NODES["R02"] == 24
    assert corrected.DEFAULT_CASE_NODES["R12"] == 24
    assert corrected.DEFAULT_CASE_NODES["R16"] == 3
    assert corrected.DEFAULT_CASE_NODES["R17"] == 128
    assert corrected.FINITE_LIMITER_DIAGNOSTIC_CASES == {"R14", "R15"}
    assert not set(ACTIVE_CASES) & {"R06", "R07", "R08", "R09"}


def test_campaign_root_is_exact_and_cannot_escape(corrected, tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)

    assert corrected.require_campaign_root(root) == root.resolve()
    with pytest.raises(corrected.CorrectedFastError, match="must be exactly"):
        corrected.require_campaign_root(tmp_path / "other")
    with pytest.raises(corrected.CorrectedFastError, match="escapes corrected campaign"):
        corrected.require_beneath(tmp_path / "legacy/restart.rst", root, "restart")


def test_node_overrides_are_selected_active_and_meshblock_limited(corrected):
    assert corrected.parse_node_overrides(["R02=27", "R17=216"], ["R02", "R17"]) == {
        "R02": 27,
        "R17": 216,
    }
    for value, selected in (
        ("R06=4", ["R06"]),
        ("R16=4", ["R16"]),
        ("R17=217", ["R17"]),
        ("R02=0", ["R02"]),
        ("R02:4", ["R02"]),
        ("R02=4", ["R03"]),
    ):
        with pytest.raises(corrected.CorrectedFastError):
            corrected.parse_node_overrides([value], selected)


def test_configure_erases_legacy_seeds_and_installs_only_corrected_cases(
    corrected, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)

    assert corrected.configure(root, {"R17": 192}) == root.resolve()
    assert corrected.fast.DEFAULT_ROOT == root.resolve()
    assert corrected.fast.FROZEN_SOURCE == EXPECTED_SOURCE
    assert corrected.fast.ATHENA == EXPECTED_ATHENA
    assert corrected.fast.CASE_SEEDS == {}
    assert set(corrected.fast.CASE_NODES) == set(ACTIVE_CASES)
    assert corrected.fast.CASE_NODES["R17"] == 192
    assert corrected.fast.seed_for_case("R03") == (0.0, None, None)


def test_source_revision_requires_exact_clean_commit(corrected, monkeypatch):
    monkeypatch.setattr(corrected.FROZEN_SOURCE.__class__, "is_dir", lambda _self: True)
    replies = iter([EXPECTED_REVISION, ""])
    monkeypatch.setattr(corrected, "git_read", lambda _args: next(replies))
    corrected.require_source_revision()

    replies = iter(["f" * 40])
    monkeypatch.setattr(corrected, "git_read", lambda _args: next(replies))
    with pytest.raises(corrected.CorrectedFastError, match="revision mismatch"):
        corrected.require_source_revision()

    replies = iter([EXPECTED_REVISION, " M src/eos/eos.hpp"])
    monkeypatch.setattr(corrected, "git_read", lambda _args: next(replies))
    with pytest.raises(corrected.CorrectedFastError, match="is dirty"):
        corrected.require_source_revision()


def test_batch_script_reauthenticates_corrected_provenance_and_root(corrected, tmp_path):
    root = tmp_path / "campaign"
    run_dir = root / "runs/E03-forcing-policy/R02/fast_s000_t0_to_t10"
    manifest = {
        "case_id": "R02",
        "sequence": 0,
        "nodes": 24,
        "input": str(EXPECTED_SOURCE / "inputs/R02.athinput"),
        "input_sha256": corrected.CASE_INPUT_SHA256["R02"],
        "restart": None,
        "restart_sha256": None,
        "output_dir": str(run_dir / "output"),
        "run_dir": str(run_dir),
        "slurm_log": str(root / "logs/%x.%j.log"),
        "run_basename": "E03_corrected_R02",
        "target_time": 10.0,
        "root": str(root),
        "fast_script": str(LAUNCHER),
        "fast_script_sha256": "a" * 64,
    }

    text = corrected.batch_script_text(manifest)

    assert "#SBATCH -J cglc_R02_s000" in text
    assert f'require_sha {EXPECTED_ATHENA_SHA256} "$ATHENA" executable' in text
    assert f'require_sha {EXPECTED_MATRIX_SHA256} "$MATRIX" stage_i_matrix' in text
    assert f'test "$source_revision" = {EXPECTED_REVISION}' in text
    assert 'status --porcelain --untracked-files=all' in text
    assert '"$CAMPAIGN_ROOT"/runs/*' in text
    assert '"$ATHENA" -i ' in text
    assert corrected.FINITE_LIMITER_OVERRIDE not in text
    assert "athenak-cgl-e03-9e075422" not in text


def test_prepare_segment_forces_fresh_t0_and_corrected_in_root_continuations(
    corrected, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)
    corrected.configure(root)
    monkeypatch.setattr(corrected, "validate_common_provenance", lambda _root: None)
    observed: list[tuple[object, ...]] = []

    def fake_prepare(root_arg, case, sequence, start_time, restart, restart_sha):
        observed.append((root_arg, case, sequence, start_time, restart, restart_sha))
        segment = root_arg / corrected.RUNS_RELATIVE / case["id"] / f"fast_s{sequence:03d}"
        write_json(
            corrected.fast.segment_manifest(segment),
            {
                "case_id": case["id"],
                "sequence": sequence,
                "start_time": start_time,
                "restart": str(restart) if restart else None,
                "restart_sha256": restart_sha,
            },
        )
        return segment

    monkeypatch.setattr(corrected, "_BASE_PREPARE_SEGMENT", fake_prepare)
    case = {"id": "R02", "name": "two", "input": "inputs/two.athinput"}

    initial = corrected.prepare_segment(root, case, 0, 0.0, None, None)
    initial_manifest = corrected.fast.load_json(corrected.fast.segment_manifest(initial))
    assert initial_manifest["launch_origin"] == "fresh_t0_corrected_production"
    assert initial_manifest["fresh_lineage_root"] is True
    assert initial_manifest["legacy_restart_permitted"] is False

    with pytest.raises(corrected.CorrectedFastError, match="start fresh from t=0"):
        corrected.prepare_segment(
            root,
            case,
            0,
            0.5,
            tmp_path / "legacy/restart.rst",
            "a" * 64,
        )
    with pytest.raises(corrected.CorrectedFastError, match="escapes corrected campaign"):
        corrected.prepare_segment(
            root,
            case,
            1,
            0.5,
            tmp_path / "legacy/restart.rst",
            "a" * 64,
        )

    restart = root / corrected.RUNS_RELATIVE / "R02/fast_s000/output/rst/rank_00000000/r.rst"
    continuation = corrected.prepare_segment(root, case, 1, 0.5, restart, "a" * 64)
    continuation_manifest = corrected.fast.load_json(
        corrected.fast.segment_manifest(continuation)
    )
    assert continuation_manifest["launch_origin"] == "corrected_campaign_continuation"
    assert continuation_manifest["fresh_lineage_root"] is False
    assert len(observed) == 2


def test_r14_r15_are_explicit_nonfatal_finite_limiter_diagnostics(
    corrected, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)
    corrected.configure(root)
    monkeypatch.setattr(corrected, "validate_common_provenance", lambda _root: None)

    def fake_prepare(root_arg, case, sequence, start_time, restart, restart_sha):
        segment = root_arg / corrected.RUNS_RELATIVE / case["id"] / f"fast_s{sequence:03d}"
        write_json(
            corrected.fast.segment_manifest(segment),
            {
                "case_id": case["id"],
                "sequence": sequence,
                "start_time": start_time,
                "restart": restart,
                "restart_sha256": restart_sha,
            },
        )
        script = segment / "manifest/run.sbatch"
        script.write_text(
            '#!/bin/bash\n"$ATHENA" -i input -d output time/tlim=10.0\n',
            encoding="utf-8",
        )
        return segment

    monkeypatch.setattr(corrected, "_BASE_PREPARE_SEGMENT", fake_prepare)
    for case_id in ("R14", "R15"):
        case = {"id": case_id, "name": case_id, "input": f"inputs/{case_id}.athinput"}
        segment = corrected.prepare_segment(root, case, 0, 0.0, None, None)
        manifest = corrected.fast.load_json(corrected.fast.segment_manifest(segment))
        text = (segment / "manifest/run.sbatch").read_text(encoding="utf-8")
        assert manifest["variant"] == corrected.FINITE_LIMITER_VARIANT
        assert manifest["strict_admissibility"] is False
        assert manifest["command_line_overrides"] == [corrected.FINITE_LIMITER_OVERRIDE]
        assert manifest["continuation_policy"] == "finite_progress_complete_terminal_products"
        assert text.count(corrected.FINITE_LIMITER_OVERRIDE) == 1


def test_validate_segment_rejects_legacy_or_nonfresh_initial_lineage(
    corrected, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)
    corrected.configure(root)
    segment = root / corrected.RUNS_RELATIVE / "R02/fast_s000_t0_to_t10"
    manifest = {
        "root": str(root.resolve()),
        "campaign_root": str(root.resolve()),
        "campaign_id": corrected.CAMPAIGN_ID,
        "source": str(corrected.FROZEN_SOURCE),
        "source_revision": corrected.SOURCE_REVISION,
        "campaign_identity": str(corrected.IDENTITY),
        "campaign_identity_sha256": corrected.IDENTITY_SHA256,
        "matrix_sha256": corrected.MATRIX_SHA256,
        "executable": str(corrected.ATHENA),
        "executable_sha256": corrected.ATHENA_SHA256,
        "input_sha256": corrected.CASE_INPUT_SHA256["R02"],
        "legacy_restart_permitted": False,
        "case_id": "R02",
        "nodes": 24,
        "ranks": 192,
        "sequence": 0,
        "start_time": 0.0,
        "restart": None,
        "restart_sha256": None,
        "launch_origin": "fresh_t0_corrected_production",
        "fresh_lineage_root": True,
        "variant": "corrected_production_strict",
        "strict_admissibility": True,
        "command_line_overrides": [],
        "continuation_policy": "strict_base_scientific_sanity",
    }
    write_json(corrected.fast.segment_manifest(segment), manifest)
    script = segment / "manifest/run.sbatch"
    script.write_text("#!/bin/bash\nstrict\n", encoding="utf-8")

    assert corrected.validate_segment(segment)["source_revision"] == EXPECTED_REVISION

    manifest["source_revision"] = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
    write_json(corrected.fast.segment_manifest(segment), manifest)
    with pytest.raises(corrected.CorrectedFastError, match="source_revision differs"):
        corrected.validate_segment(segment)

    manifest["source_revision"] = corrected.SOURCE_REVISION
    manifest["start_time"] = 0.5
    manifest["restart"] = "/legacy/restart.rst"
    manifest["restart_sha256"] = "a" * 64
    write_json(corrected.fast.segment_manifest(segment), manifest)
    with pytest.raises(corrected.CorrectedFastError, match="not fresh from t=0"):
        corrected.validate_segment(segment)


def test_diagnostic_analysis_continues_on_finite_complete_products_not_hard_bound_zero(
    corrected, tmp_path, monkeypatch
):
    segment = tmp_path / "campaign/runs/E03-forcing-policy/R14/fast_s000"
    output = segment / "output"
    manifest = {
        "case_id": "R14",
        "ranks": 8,
        "start_time": 0.0,
        "target_time": 10.0,
    }
    write_history(
        output / "fixture.mhd.hst",
        {
            "time": [0.0, 0.5],
            "mass": [1.0, 1.0],
            "lf_dfloor": [0.0, 0.0],
            "lf_pfloor": [0.0, 0.0],
            "lf_nonfin": [0.0, 0.0],
            "lf_nonpos": [0.0, 0.0],
            "lf_hardbd": [0.0, 42.0],
        },
    )
    write_history(
        output / "fixture.user.hst",
        {"time": [0.0, 0.5], "mass": [1.0, 1.0]},
    )
    write_restart_group(output / "rst", 8, 0.5)
    monkeypatch.setattr(corrected, "validate_segment", lambda _segment: manifest)

    result = corrected.analyze_segment(segment)

    assert result["base_strict_passed"] is False
    assert result["strict_lf_failure_maxima"]["lf_hardbd"] == 42.0
    assert result["continuation_gate"]["hard_bound_zero_required"] is False
    assert result["continuation_gate"]["terminal_snapshot_required"] is False
    assert result["terminal_snapshot"] is None
    assert result["terminal_snapshot_error"] is not None
    assert result["passed"] is True
    saved = corrected.fast.load_json(segment / "manifest/fast_analysis.json")
    assert saved["passed"] is True


def test_preflight_authenticates_matrix_executable_and_every_selected_input(
    corrected, tmp_path, monkeypatch
):
    root = tmp_path / "campaign"
    source = tmp_path / "source"
    executable = tmp_path / "athena"
    monkeypatch.setattr(corrected, "CAMPAIGN_ROOT", root)
    monkeypatch.setattr(corrected, "FROZEN_SOURCE", source)
    monkeypatch.setattr(corrected, "ATHENA", executable)
    monkeypatch.setattr(corrected, "require_source_revision", lambda: None)
    observed: list[tuple[Path, str, str]] = []
    monkeypatch.setattr(
        corrected.fast,
        "require_sha",
        lambda path, expected, label: observed.append((path, expected, label)),
    )
    monkeypatch.setattr(
        corrected.fast,
        "campaign_cases",
        lambda: {
            "R02": {"id": "R02", "input": "inputs/R02.athinput"},
            "R17": {"id": "R17", "input": "inputs/R17.athinput"},
        },
    )
    corrected.configure(root)

    result = corrected.validate_provenance(root, ["R02", "R17"])

    assert result["active_cases"] == ["R02", "R17"]
    assert result["case_nodes"] == {"R02": 24, "R17": 128}
    assert observed == [
        (corrected.IDENTITY, corrected.IDENTITY_SHA256, "corrected campaign identity"),
        (source / corrected.MATRIX_RELATIVE, corrected.MATRIX_SHA256, "Stage I matrix"),
        (executable, corrected.ATHENA_SHA256, "corrected executable"),
        (
            source / "inputs/R02.athinput",
            corrected.CASE_INPUT_SHA256["R02"],
            "R02 corrected frozen input",
        ),
        (
            source / "inputs/R17.athinput",
            corrected.CASE_INPUT_SHA256["R17"],
            "R17 corrected frozen input",
        ),
    ]


def test_submit_is_direct_and_only_runs_when_explicitly_requested(
    corrected, tmp_path, monkeypatch
):
    segment = tmp_path / "segment"
    write_json(
        corrected.fast.segment_manifest(segment),
        {"case_id": "R02", "job_id": None},
    )
    script = segment / "manifest/run.sbatch"
    script.write_text("#!/bin/bash\n", encoding="utf-8")
    calls: list[list[str]] = []

    def fake_run(argv, **_kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout="12345;frontier\n")

    monkeypatch.setattr(corrected.fast.subprocess, "run", fake_run)
    assert corrected.fast.submit_segment(segment) == "12345"
    assert calls == [["/usr/bin/sbatch", "--parsable", str(script)]]
