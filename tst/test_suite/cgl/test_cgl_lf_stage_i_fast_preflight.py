"""Regression tests for the final direct-fast Stage I analysis preflight."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
PREFLIGHT_PATH = (
    REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_preflight.py"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def preflight():
    return load_module("cgl_lf_stage_i_fast_preflight_test", PREFLIGHT_PATH)


@pytest.fixture(scope="module")
def acceptance(preflight):
    return preflight.load_acceptance_module()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def create_archive(root: Path) -> Path:
    archive = root / "source.tar.gz"
    source = root / "source/MKS24.tex"
    archive.parent.mkdir(parents=True, exist_ok=True)
    source.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(b"bound source archive\n")
    source.write_text("bound source\n", encoding="utf-8")
    manifest = root / "manifest.json"
    write_json(
        manifest,
        {
            "schema_version": 1,
            "archive": {"path": archive.name, "sha256": sha256(archive)},
            "extracted_file_sha256": {
                "source/MKS24.tex": sha256(source),
            },
            "source_tex": "source/MKS24.tex",
        },
    )
    return manifest


class PolicyAcceptance:
    """Delegate verification primitives while supplying one fixture policy."""

    def __init__(self, acceptance, policy):
        self.acceptance = acceptance
        self.policy = policy
        self.calls: list[tuple[Path, Path]] = []

    def load_validated_policy(self, criteria: Path, review: Path):
        self.calls.append((criteria, review))
        return self.policy

    def __getattr__(self, name: str):
        return getattr(self.acceptance, name)


def build_fixture(tmp_path: Path, acceptance):
    main_manifest = create_archive(tmp_path / "reference")
    verified_manifest = create_archive(tmp_path / "reference-verified")

    product_dir = main_manifest.parent / "digitized_fixture"
    data = product_dir / "fixture.csv"
    data.parent.mkdir(parents=True, exist_ok=True)
    data.write_text(
        "x,y,y_uncertainty\n0,1,0.1\n1,2,0.2\n",
        encoding="utf-8",
    )
    product_manifest = product_dir / "curves.json"
    write_json(
        product_manifest,
        {
            "schema_version": 1,
            "curves": [
                {
                    "id": "fixture_product",
                    "case": "fixture_case",
                    "product": "fixture.metric",
                    "data_file": data.name,
                    "data_sha256": sha256(data),
                    "interpolation": "linear",
                }
            ],
        },
    )

    panel_status = {
        "reference_manifests": {
            "fixture_manifest": {
                "path": "digitized_fixture/curves.json",
                "sha256": sha256(product_manifest),
            }
        },
        "reference_product_bindings": {
            "fixture_product": {
                "kind": "curve",
                "case": "fixture_case",
                "product": "fixture.metric",
                "data_file": data.name,
                "data_sha256": sha256(data),
                "reference_manifest": "fixture_manifest",
            }
        },
        "panels": [
            {
                "id": "fixture_panel",
                "disposition": "comparison",
                "reference_products": ["fixture_product"],
            }
        ],
    }
    stage_manifest = tmp_path / "stage_manifest.json"
    write_json(stage_manifest, {"panel_status": panel_status})
    criteria = tmp_path / "criteria.json"
    review = tmp_path / "review.json"
    criteria.write_text("{}\n", encoding="utf-8")
    review.write_text("{}\n", encoding="utf-8")
    policy = {
        "criteria": {"comparison_panels": [{"id": "fixture_panel"}]},
        "criteria_binding": acceptance.regular_file_binding(criteria, "criteria"),
        "review_binding": acceptance.regular_file_binding(review, "review"),
        "manifest": {"panel_status": panel_status},
        "verified_sources": {
            "acceptance_utility": acceptance.regular_file_binding(
                Path(acceptance.__file__), "acceptance utility"
            ),
            "stage_i_manifest": acceptance.regular_file_binding(
                stage_manifest, "Stage I manifest"
            ),
            "reference_archive_manifest": acceptance.regular_file_binding(
                main_manifest, "reference archive manifest"
            ),
            "verified_reference_archive_manifest": acceptance.regular_file_binding(
                verified_manifest, "verified reference archive manifest"
            ),
        },
    }
    return {
        "policy": policy,
        "criteria": criteria,
        "review": review,
        "data": data,
        "verified_source": verified_manifest.parent / "source/MKS24.tex",
    }


def test_preflight_reuses_policy_verification_and_is_deterministic(
    tmp_path: Path, preflight, acceptance
) -> None:
    fixture = build_fixture(tmp_path, acceptance)
    delegated = PolicyAcceptance(acceptance, fixture["policy"])
    before = {
        str(path.relative_to(tmp_path)): path.read_bytes()
        for path in sorted(tmp_path.rglob("*"))
        if path.is_file()
    }

    first = preflight.build_preflight(
        fixture["criteria"], fixture["review"], acceptance=delegated
    )
    second = preflight.build_preflight(
        fixture["criteria"], fixture["review"], acceptance=delegated
    )
    after = {
        str(path.relative_to(tmp_path)): path.read_bytes()
        for path in sorted(tmp_path.rglob("*"))
        if path.is_file()
    }

    assert first == second
    assert before == after
    assert first["result"] == "pass"
    assert first["read_only"] is True
    assert first["fail_closed"] is True
    assert first["summary"] == {
        "verified_archive_count": 2,
        "verified_archive_file_count": 4,
        "verified_reference_manifest_count": 1,
        "verified_reference_product_count": 1,
    }
    assert delegated.calls == [
        (fixture["criteria"], fixture["review"]),
        (fixture["criteria"], fixture["review"]),
    ]
    assert "Result: **pass**." in preflight.render_markdown(first)


def test_preflight_fails_closed_when_archive_member_is_missing(
    tmp_path: Path, preflight, acceptance
) -> None:
    fixture = build_fixture(tmp_path, acceptance)
    fixture["verified_source"].unlink()
    delegated = PolicyAcceptance(acceptance, fixture["policy"])

    with pytest.raises(preflight.PreflightError, match="does not resolve"):
        preflight.build_preflight(
            fixture["criteria"], fixture["review"], acceptance=delegated
        )


def test_preflight_fails_closed_when_reference_product_digest_drifts(
    tmp_path: Path, preflight, acceptance
) -> None:
    fixture = build_fixture(tmp_path, acceptance)
    fixture["data"].write_text(
        "x,y,y_uncertainty\n0,1,0.1\n1,3,0.2\n",
        encoding="utf-8",
    )
    delegated = PolicyAcceptance(acceptance, fixture["policy"])

    with pytest.raises(acceptance.AcceptanceError, match="data SHA-256 differs"):
        preflight.build_preflight(
            fixture["criteria"], fixture["review"], acceptance=delegated
        )


def test_preflight_rejects_reference_manifest_escape(
    tmp_path: Path, preflight, acceptance
) -> None:
    fixture = build_fixture(tmp_path, acceptance)
    outside = tmp_path / "outside.json"
    write_json(outside, {"schema_version": 1})
    declaration = fixture["policy"]["manifest"]["panel_status"][
        "reference_manifests"
    ]["fixture_manifest"]
    declaration["path"] = "../outside.json"
    declaration["sha256"] = sha256(outside)
    delegated = PolicyAcceptance(acceptance, fixture["policy"])

    with pytest.raises(preflight.PreflightError, match="escapes its bound archive"):
        preflight.build_preflight(
            fixture["criteria"], fixture["review"], acceptance=delegated
        )


def test_command_reports_deterministic_failure_for_missing_policy(
    tmp_path: Path, preflight, capsys
) -> None:
    criteria = tmp_path / "missing-criteria.json"
    review = tmp_path / "missing-review.json"

    status = preflight.main(
        [
            "--criteria",
            str(criteria),
            "--criteria-review",
            str(review),
            "--format",
            "json",
        ]
    )
    record = json.loads(capsys.readouterr().out)

    assert status == 2
    assert record["result"] == "fail"
    assert record["read_only"] is True
    assert record["fail_closed"] is True
    assert str(criteria) in record["error"]
    assert list(tmp_path.iterdir()) == []
