"""Focused tests for corrected-production campaign identity isolation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
TOOL = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_corrected_identity.py"
REVISION = "0c406312fa35d5e1c7041d80333b0ea24b0127ae"


def load_tool():
    name = "cgl_lf_stage_i_corrected_identity_tests"
    spec = importlib.util.spec_from_file_location(name, TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def identity_tool():
    return load_tool()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def campaign(tmp_path: Path) -> dict[str, object]:
    artifacts = {}
    for name in ("audit", "eos", "executable", "matrix"):
        path = tmp_path / "authority" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{name} authority\n", encoding="utf-8")
        artifacts[name] = path
    return {
        "id": "mks24-stage-i-fast-corrected-E03",
        "root": tmp_path / "runs/mks24-stage-i-fast-corrected/E03-forcing-policy",
        "artifacts": artifacts,
        "run_paths": ["selected/R17", "selected/R02", "selected/R02"],
    }


def write_campaign(identity_tool, campaign: dict[str, object]) -> Path:
    return identity_tool.write_identity(
        campaign["id"],
        campaign["root"],
        REVISION,
        campaign["artifacts"],
        campaign["run_paths"],
    )


def test_write_and_validate_binds_artifacts_and_selected_paths(
    identity_tool, campaign: dict[str, object]
) -> None:
    identity_path = write_campaign(identity_tool, campaign)
    root = Path(campaign["root"]).resolve()
    retained = json.loads(identity_path.read_text(encoding="utf-8"))

    assert identity_path == root / "campaign-identity.json"
    assert retained["campaign_kind"] == "corrected-production"
    assert retained["campaign_root"] == str(root)
    assert retained["source_revision"] == REVISION
    assert retained["selected_run_paths"] == ["selected/R02", "selected/R17"]
    for name, artifact in campaign["artifacts"].items():
        assert retained["artifacts"][name] == {
            "path": str(Path(artifact).resolve()),
            "sha256": sha256(Path(artifact)),
        }

    result = identity_tool.validate_identity(
        identity_path, [root / "selected/R02", "selected/R17"]
    )
    assert result["result"] == "pass"
    assert result["requested_run_paths"] == ["selected/R02", "selected/R17"]


@pytest.mark.parametrize("run_path", ["../legacy/R02", "/tmp/legacy/R02"])
def test_write_rejects_run_path_outside_corrected_root(
    identity_tool, campaign: dict[str, object], run_path: str
) -> None:
    campaign["run_paths"] = [run_path]
    with pytest.raises(identity_tool.IdentityError, match="escapes corrected campaign"):
        write_campaign(identity_tool, campaign)


def test_validate_rejects_unselected_and_outside_requested_paths(
    identity_tool, campaign: dict[str, object]
) -> None:
    identity_path = write_campaign(identity_tool, campaign)
    with pytest.raises(identity_tool.IdentityError, match="not selected by identity"):
        identity_tool.validate_identity(identity_path, ["selected/R03"])
    with pytest.raises(identity_tool.IdentityError, match="escapes corrected campaign"):
        identity_tool.validate_identity(identity_path, ["../legacy/R02"])


def test_validate_rejects_artifact_drift(
    identity_tool, campaign: dict[str, object]
) -> None:
    identity_path = write_campaign(identity_tool, campaign)
    eos = Path(campaign["artifacts"]["eos"])
    eos.write_text("changed EOS\n", encoding="utf-8")
    with pytest.raises(identity_tool.IdentityError, match="EOS|eos artifact SHA-256"):
        identity_tool.validate_identity(identity_path)


def test_validate_rejects_malformed_schema_marker(
    identity_tool, campaign: dict[str, object]
) -> None:
    identity_path = write_campaign(identity_tool, campaign)
    retained = json.loads(identity_path.read_text(encoding="utf-8"))
    retained["schema_version"] = True
    identity_path.write_text(
        json.dumps(retained, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(identity_tool.IdentityError, match="schema or schema_version"):
        identity_tool.validate_identity(identity_path)


def test_validate_rejects_selected_symlink_escape(
    identity_tool, campaign: dict[str, object], tmp_path: Path
) -> None:
    campaign["run_paths"] = ["selected/R02"]
    identity_path = write_campaign(identity_tool, campaign)
    root = Path(campaign["root"])
    outside = tmp_path / "legacy/R02"
    outside.mkdir(parents=True)
    selected = root / "selected"
    selected.mkdir()
    (selected / "R02").symlink_to(outside, target_is_directory=True)

    with pytest.raises(identity_tool.IdentityError, match="escapes corrected campaign"):
        identity_tool.validate_identity(identity_path)


def test_existing_identity_is_idempotent_but_cannot_be_replaced(
    identity_tool, campaign: dict[str, object]
) -> None:
    identity_path = write_campaign(identity_tool, campaign)
    original = identity_path.read_bytes()
    assert write_campaign(identity_tool, campaign) == identity_path
    assert identity_path.read_bytes() == original

    campaign["id"] = "different-corrected-campaign"
    with pytest.raises(identity_tool.IdentityError, match="refusing to replace"):
        write_campaign(identity_tool, campaign)


def test_cli_write_then_validate(
    identity_tool, campaign: dict[str, object], capsys
) -> None:
    artifacts = campaign["artifacts"]
    write_args = [
        "write",
        "--campaign-id",
        campaign["id"],
        "--campaign-root",
        str(campaign["root"]),
        "--source-revision",
        REVISION,
        "--executable",
        str(artifacts["executable"]),
        "--eos",
        str(artifacts["eos"]),
        "--matrix",
        str(artifacts["matrix"]),
        "--audit",
        str(artifacts["audit"]),
        "--run-path",
        "selected/R02",
        "--run-path",
        "selected/R17",
    ]
    assert identity_tool.main(write_args) == 0
    written = json.loads(capsys.readouterr().out)
    identity_path = Path(written["identity"]["path"])

    assert identity_tool.main(
        ["validate", "--identity", str(identity_path), "--run-path", "selected/R17"]
    ) == 0
    validated = json.loads(capsys.readouterr().out)
    assert validated["requested_run_paths"] == ["selected/R17"]
