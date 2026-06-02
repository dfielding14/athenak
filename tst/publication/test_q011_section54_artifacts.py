#!/usr/bin/env python3
"""Focused tests for immutable Q-011 Section 5.4 derived artifacts."""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest

from . import q011_section54_artifacts as artifacts


_RAW_DIGEST = "a" * 64
_ANALYZER_DIGEST = "b" * 64


def _pending() -> dict[str, object]:
    return {
        "status": "pending_external_review",
        "reviewer": None,
        "reviewed_at_utc": None,
        "notes": "Named external review remains explicitly deferred.",
    }


def _manifest() -> dict[str, object]:
    return artifacts.build_derived_manifest(
        campaign_id="q011-section54-fixture",
        raw_artifact_inventory_sha256=_RAW_DIGEST,
        analyzer_bindings={"tst/publication/analyze_fixture.py": _ANALYZER_DIGEST},
        reviewer_disposition=_pending(),
    )


def _make_writable(root: Path) -> None:
    for directory, names, filenames in os.walk(root):
        os.chmod(directory, 0o755)
        for name in names:
            os.chmod(Path(directory) / name, 0o755)
        for name in filenames:
            os.chmod(Path(directory) / name, 0o644)


class Q011Section54DerivedArtifactTests(unittest.TestCase):
    def test_publish_verify_and_refuse_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            receipt = artifacts.publish_derived_bundle(
                root,
                manifest=_manifest(),
                artifacts={
                    "metrics/summary.json": artifacts.canonical_json_bytes({"pass": True}),
                    "figures/morphology.png": b"png-fixture",
                },
            )
            self.assertEqual(
                artifacts.verify_published_derived_bundle(
                    root, expected_inventory_sha256=receipt["inventory_sha256"]
                ),
                _manifest(),
            )
            self.assertFalse(root.stat().st_mode & 0o222)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "already exists"):
                artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)

    def test_verify_rejects_tamper_and_tree_addition(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(
                root, manifest=_manifest(), artifacts={"metrics/result.json": b"{}\n"}
            )
            _make_writable(root)
            (root / "metrics/result.json").write_bytes(b'{"drift": true}\n')
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "checksum drifted"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            (root / "undeclared.txt").write_text("unexpected", encoding="utf-8")
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "tree closure"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_publication_rejects_unsafe_reserved_and_nonfinite_members(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for name in ("../escape", "/absolute", artifacts.INVENTORY_NAME):
                with self.subTest(name=name), self.assertRaises(artifacts.DerivedArtifactError):
                    artifacts.publish_derived_bundle(
                        Path(directory) / name.replace("/", "_"),
                        manifest=_manifest(),
                        artifacts={name: b"payload"},
                    )
        with self.assertRaises(artifacts.DerivedArtifactError):
            artifacts.canonical_json_bytes({"not_finite": float("nan")})

    def test_reviewer_disposition_requires_named_terminal_reviewer(self) -> None:
        accepted = {
            "status": "accepted",
            "reviewer": "Named Reviewer",
            "reviewed_at_utc": "2026-06-02T00:00:00Z",
            "notes": "Reviewed.",
        }
        self.assertEqual(artifacts.validate_reviewer_disposition(accepted), accepted)
        for disposition in (
            {**accepted, "reviewer": ""},
            {**accepted, "reviewed_at_utc": None},
            {**_pending(), "reviewer": "Premature Name"},
            {**_pending(), "status": "unknown"},
        ):
            with self.subTest(disposition=disposition), self.assertRaises(
                artifacts.DerivedArtifactError
            ):
                artifacts.validate_reviewer_disposition(disposition)

    def test_verify_rejects_duplicate_inventory_keys(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            inventory = root / artifacts.INVENTORY_NAME
            inventory.write_text(
                '{"record_type":"q011_section54_derived_artifact_inventory",'
                '"schema_version":1,"schema_version":1,"members":[]}',
                encoding="utf-8",
            )
            os.chmod(inventory, 0o444)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "duplicate JSON key"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_numeric_schema_aliases_fail_closed(self) -> None:
        with self.assertRaisesRegex(artifacts.DerivedArtifactError, "identity drifted"):
            artifacts.validate_derived_manifest({**_manifest(), "schema_version": True})

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            inventory = root / artifacts.INVENTORY_NAME
            value = json.loads(inventory.read_text(encoding="utf-8"))
            value["schema_version"] = True
            inventory.write_bytes(artifacts.canonical_json_bytes(value))
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "schema drifted"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)


if __name__ == "__main__":
    unittest.main()
