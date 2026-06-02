#!/usr/bin/env python3
"""Focused adversarial tests for Q-011 pressure-pilot raw-bundle publication."""

from __future__ import annotations

import base64
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Iterator
import unittest
from unittest.mock import patch

from tst.publication import analyze_q011_section54_pressure_pilot as pilot
from tst.publication import analyze_q011_section54_pressure_pilot_case as case_verifier
from tst.publication import publish_q011_section54_pressure_pilot_bundle as publisher
from tst.publication import test_analyze_q011_section54_pressure_pilot as pilot_fixture
from tst.publication.frontier_control_plane.control_plane_common import (
    PRODUCTION_RUNTIME_LOADED_MODULES,
)
from tst.publication.frontier_control_plane.control_plane_common import (
    PRODUCTION_RUNTIME_MODULEFILES,
)
from tst.publication.frontier_control_plane.control_plane_common import (
    PRODUCTION_RUNTIME_MODULEPATH,
)


PUBLICATION_DIR = Path(__file__).resolve().parent


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _put(root: Path, relative: str, payload: bytes) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _allowlist() -> bytes:
    values = {
        **{key: "<unset>" for key in case_verifier._RUNTIME_ALLOWLIST_KEYS},
        **case_verifier._EXPECTED_BASELINE,
        "LOADEDMODULES": ":".join(PRODUCTION_RUNTIME_LOADED_MODULES),
        "_LMFILES_": ":".join(PRODUCTION_RUNTIME_MODULEFILES),
        "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH,
    }
    return "".join(
        f"{key}={values[key]}\n" for key in case_verifier._RUNTIME_ALLOWLIST_KEYS
    ).encode("utf-8")


def _stderr() -> bytes:
    return base64.b64decode(
        (
            PUBLICATION_DIR / "readiness/frontier_mpich_diagnostic_stderr.base64"
        ).read_bytes().replace(b"\n", b""),
        validate=True,
    )


def _stdout() -> bytes:
    return (
        "PIC trusted GPU launch: rank=0 host=frontier00001 "
        "ROCR_VISIBLE_DEVICES=0 linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
    ).encode("ascii") + pilot_fixture._stdout()


def _restart(root: Path, case_id: str, index: int) -> None:
    suffix = f"{index:05d}"
    artifact_relative = f"output/rst/{case_id}.{suffix}.rst"
    manifest_relative = artifact_relative + ".manifest"
    artifact = f"terminal restart {case_id} {suffix}\n".encode("ascii")
    manifest = {
        "schema": "ATHENAK_RESTART_MANIFEST_V1",
        "members": [
            {
                "path": f"rst/{case_id}.{suffix}.rst",
                "size": len(artifact),
                "fnv1a64": pilot_fixture._fnv1a64(artifact),
            }
        ],
    }
    manifest_payload = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    _put(root, artifact_relative, artifact)
    _put(root, artifact_relative + ".complete", pilot_fixture._marker(artifact))
    _put(root, manifest_relative, manifest_payload)
    _put(
        root,
        manifest_relative + ".complete",
        pilot_fixture._marker(manifest_payload),
    )


def _publish_inventory(root: Path) -> None:
    (root / "analysis").mkdir(mode=0o700)
    (root / "analysis").chmod(0o700)
    records = []
    for path in sorted(root.rglob("*")):
        if "analysis" in path.relative_to(root).parts or not path.is_file():
            continue
        payload = path.read_bytes()
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "sha256": _sha256(payload),
                "size": len(payload),
            }
        )
        path.chmod(0o444)
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_dir() and path != root / "analysis":
            path.chmod(0o555)
    inventory = root / "artifact_inventory.json"
    inventory.write_text(
        json.dumps({"schema_version": 1, "files": records}, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    inventory.chmod(0o444)
    root.chmod(0o555)


def _raw_tree(
    root: Path,
    case_id: str,
    *,
    extra: bool = False,
    false_stdout_checksum: bool = False,
    legacy_allowlist_name: bool = False,
    omit_stdout_checksum: bool = False,
) -> None:
    root.mkdir()
    action_id = case_id if legacy_allowlist_name else case_verifier._action_id(case_id)
    _put(root, f"{action_id}.environment.allowlist.txt", _allowlist())
    stdout = _stdout()
    _put(root, "athena_stdout.txt", stdout)
    stdout_digest = "0" * 64 if false_stdout_checksum else _sha256(stdout)
    if not omit_stdout_checksum:
        _put(root, "athena_stdout.sha256", (stdout_digest + "\n").encode("ascii"))
    _put(root, "athena_stderr.txt", _stderr())
    _ps_p0, argv_value = case_verifier._CASE_BY_ID[case_id]
    for index, time in enumerate(case_verifier._TIMES):
        sources = case_verifier._snapshot_source_paths(case_id, index)
        _put(root, sources["mhd_w_bcc"], pilot_fixture._binary(time, argv_value, pilot._MHD_FIELDS))
        _put(root, sources["bmag"], pilot_fixture._binary(time, argv_value, ("bmag",)))
        _put(root, sources["prtcl_jx"], pilot_fixture._binary(time, argv_value, ("prtcl_jx",)))
        _put(root, sources["j2"], pilot_fixture._binary(time, argv_value, ("j2",)))
        _put(root, sources["prtcl_all"], pilot_fixture._particle_vtk(time))
        _restart(root, case_id, index)
    if extra:
        _put(root, "output/undeclared.txt", b"must fail closure\n")
    _publish_inventory(root)


def _make_writable(root: Path) -> None:
    if not os.path.lexists(root):
        return
    for directory, names, filenames in os.walk(root):
        base = Path(directory)
        base.chmod(0o755)
        for name in names:
            (base / name).chmod(0o755)
        for name in filenames:
            (base / name).chmod(0o644)


def _freeze_existing(root: Path) -> None:
    for path in sorted(root.rglob("*")):
        path.chmod(0o444 if path.is_file() else 0o555)
    root.chmod(0o555)


def _descriptor_sha256(descriptor: dict[str, object]) -> str:
    return _sha256(case_verifier.canonical_json_bytes(descriptor))


@contextmanager
def _verified_raw_cases() -> Iterator[tuple[Path, dict[str, Path], dict[str, str]]]:
    with tempfile.TemporaryDirectory() as directory:
        pic_root = Path(directory) / "pic"
        base = pic_root / "publication"
        runs = pic_root / "runs"
        base.mkdir(parents=True)
        runs.mkdir()
        roots = {}
        digests = {}
        for case_id in case_verifier.CASE_IDS:
            root = runs / f"raw-{case_id}"
            _raw_tree(root, case_id)
            descriptor = case_verifier.publish_case_descriptor(root, case_id)
            roots[case_id] = root
            digests[case_id] = _descriptor_sha256(descriptor)
        try:
            yield base, roots, digests
        finally:
            for root in roots.values():
                _make_writable(root)


class Q011Section54PressurePilotBundlePublicationTests(unittest.TestCase):
    def test_case_descriptor_requires_execution_tranche_action_id_and_stdout_checksum(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "raw"
            case_id = case_verifier.CASE_IDS[0]
            _raw_tree(root, case_id)
            descriptor = case_verifier.analyze_case_tree(root, case_id)
            allowlist = "q011-pressure-ps-p0-1p00.environment.allowlist.txt"
            self.assertIn(allowlist, descriptor["runtime_artifacts"])
            self.assertIn("athena_stdout.sha256", descriptor["runtime_artifacts"])
            terminal = descriptor["manifest_case"]["terminal_restart"]
            self.assertEqual(
                set(terminal),
                {"time", "manifest", "manifest_complete", "members"},
            )
            _make_writable(root)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "raw"
            case_id = case_verifier.CASE_IDS[0]
            _raw_tree(root, case_id, false_stdout_checksum=True)
            with self.assertRaisesRegex(
                case_verifier.PressurePilotCaseError, "does not bind"
            ):
                case_verifier.analyze_case_tree(root, case_id)
            _make_writable(root)

        for variant in ("legacy allowlist name", "omitted stdout checksum"):
            with self.subTest(variant=variant), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "raw"
                case_id = case_verifier.CASE_IDS[0]
                _raw_tree(
                    root,
                    case_id,
                    legacy_allowlist_name=variant == "legacy allowlist name",
                    omit_stdout_checksum=variant == "omitted stdout checksum",
                )
                with self.assertRaisesRegex(ValueError, "omits required path"):
                    case_verifier.analyze_case_tree(root, case_id)
                _make_writable(root)

    def test_case_descriptor_and_aggregate_publication_preserve_exact_closure(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt_path = base / "publication-receipt.json"
            receipt = publisher.publish_pressure_pilot_bundle(
                bundle,
                receipt_path=receipt_path,
                analysis_result_path=base / "aggregate-analysis.json",
                case_artifact_dirs=roots,
                case_descriptor_sha256=digests,
                authorized_pic_root=base.parent,
            )
            result = publisher.verify_published_pressure_pilot_bundle(
                bundle,
                receipt["manifest_sha256"],
                authorized_publication_root=base,
            )
            self.assertEqual(result["status"], "pass_engineering_calibration_only")
            manifest = pilot._manifest_schema(
                (bundle / pilot.MANIFEST_NAME).read_bytes()
            )
            actual = {
                path.relative_to(bundle).as_posix()
                for path in bundle.rglob("*")
                if path.is_file()
            }
            self.assertEqual(actual, pilot._declared_paths(manifest))
            self.assertNotIn("artifact_inventory.json", actual)
            self.assertFalse(any("environment.allowlist" in path for path in actual))
            self.assertFalse(any(path.endswith("-errs.dat") for path in actual))
            self.assertFalse(any(path.endswith("athena_stderr.txt") for path in actual))
            publication = json.loads(receipt_path.read_text(encoding="utf-8"))
            self.assertEqual(
                receipt["receipt_sha256"], _sha256(receipt_path.read_bytes())
            )
            analysis_result = Path(receipt["analysis_result_path"])
            self.assertEqual(
                receipt["analysis_result_sha256"],
                _sha256(analysis_result.read_bytes()),
            )
            self.assertEqual(
                publication["aggregate_analysis"],
                {
                    "path": str(analysis_result),
                    "sha256": receipt["analysis_result_sha256"],
                },
            )
            self.assertEqual(
                set(publication["source_bindings"]),
                {
                    "registered_execution_preregistration",
                    "publisher",
                    "aggregate_analyzer",
                },
            )
            self.assertEqual(
                [record["case_id"] for record in publication["raw_cases"]],
                list(case_verifier.CASE_IDS),
            )
            self.assertEqual(
                {
                    record["case_id"]: record["descriptor_sha256"]
                    for record in publication["raw_cases"]
                },
                digests,
            )
            for path in [bundle, *bundle.rglob("*")]:
                self.assertFalse(path.stat().st_mode & 0o222)
            self.assertFalse(receipt_path.stat().st_mode & 0o222)
            self.assertFalse(analysis_result.stat().st_mode & 0o222)
            retained = publisher.verify_published_pressure_pilot_receipt(
                receipt_path, authorized_pic_root=base.parent
            )
            self.assertEqual(retained["receipt_sha256"], receipt["receipt_sha256"])
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "already exists"
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=base / "second-publication-receipt.json",
                    analysis_result_path=base / "second-aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            _make_writable(bundle)

    def test_case_verifier_rejects_raw_payload_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "raw"
            case_id = case_verifier.CASE_IDS[0]
            _raw_tree(root, case_id)
            path = root / case_verifier._snapshot_source_paths(case_id, 0)["bmag"]
            path.chmod(0o644)
            path.write_bytes(path.read_bytes() + b"forged")
            path.chmod(0o444)
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                case_verifier.analyze_case_tree(root, case_id)
            _make_writable(root)

    def test_case_verifier_rejects_inventory_closed_but_unregistered_raw_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "raw"
            case_id = case_verifier.CASE_IDS[0]
            _raw_tree(root, case_id, extra=True)
            with self.assertRaisesRegex(
                case_verifier.PressurePilotCaseError, "tree closure drifted"
            ):
                case_verifier.analyze_case_tree(root, case_id)
            _make_writable(root)

    def test_publisher_recomputes_and_rejects_descriptor_tamper(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            case_id = case_verifier.CASE_IDS[0]
            descriptor_path = roots[case_id] / case_verifier.CASE_DESCRIPTOR_PATH
            descriptor_path.chmod(0o644)
            descriptor = json.loads(descriptor_path.read_text(encoding="utf-8"))
            descriptor["runtime_profile"] = "forged"
            descriptor_path.write_bytes(case_verifier.canonical_json_bytes(descriptor))
            descriptor_path.chmod(0o444)
            digests[case_id] = _sha256(descriptor_path.read_bytes())
            with self.assertRaisesRegex(
                case_verifier.PressurePilotCaseError, "value drifted"
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )

    def test_final_verifier_rejects_payload_tamper_and_tree_addition(self) -> None:
        for variant in ("payload", "addition"):
            with self.subTest(variant=variant), _verified_raw_cases() as (base, roots, digests):
                bundle = base / "bundle"
                receipt = publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
                _make_writable(bundle)
                if variant == "payload":
                    path = bundle / "cases/ps_p0_1p00/stdout.txt"
                    path.write_bytes(path.read_bytes() + b"forged\n")
                else:
                    _put(bundle, "cases/ps_p0_1p00/undeclared.txt", b"forged\n")
                _freeze_existing(bundle)
                with self.assertRaises(publisher.PressurePilotPublicationError):
                    publisher.verify_published_pressure_pilot_bundle(
                        bundle,
                        receipt["manifest_sha256"],
                        authorized_publication_root=base,
                    )
                _make_writable(bundle)

    def test_publisher_rejects_incomplete_case_set(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            del roots[case_verifier.CASE_IDS[-1]]
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "source tree set drifted"
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )

    def test_publisher_rejects_off_root_output_and_raw_case(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            outside = base.parent.parent / "outside"
            outside.mkdir()
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "outside the authorized PIC publication root",
            ):
                publisher.publish_pressure_pilot_bundle(
                    outside / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            roots[case_verifier.CASE_IDS[0]] = base
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "outside the authorized PIC runs root",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )

    def test_retained_receipt_rejects_source_binding_tamper(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt_path = base / "publication-receipt.json"
            publisher.publish_pressure_pilot_bundle(
                base / "bundle",
                receipt_path=receipt_path,
                analysis_result_path=base / "aggregate-analysis.json",
                case_artifact_dirs=roots,
                case_descriptor_sha256=digests,
                authorized_pic_root=base.parent,
            )
            receipt_path.chmod(0o600)
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["source_bindings"]["publisher"]["sha256"] = "0" * 64
            receipt_path.write_bytes(publisher._canonical_json_bytes(receipt))
            receipt_path.chmod(0o444)
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "source bindings drifted"
            ):
                publisher.verify_published_pressure_pilot_receipt(
                    receipt_path, authorized_pic_root=base.parent
                )

    def test_analyzer_and_publisher_require_registered_source_tranche(self) -> None:
        with _verified_raw_cases() as (base, roots, digests), patch.object(
            pilot.execution,
            "validate_source_tranche",
            side_effect=pilot.execution.ContractError("source drift"),
        ):
            with self.assertRaisesRegex(pilot.execution.ContractError, "source drift"):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            with self.assertRaisesRegex(pilot.execution.ContractError, "source drift"):
                pilot.analyze_pressure_pilot_bundle(base, "0" * 64)

    def test_descriptor_relative_rename_fails_closed_without_atomic_no_replace(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), patch.object(publisher.os, "rename") as rename:
                    with self.assertRaisesRegex(
                        publisher.PressurePilotPublicationError,
                        "requires atomic no-replace rename support",
                    ):
                        publisher._rename_no_replace_at(
                            descriptor, "source", "destination"
                        )
                rename.assert_not_called()
                self.assertTrue((parent / "source").is_file())
                self.assertFalse((parent / "destination").exists())
            finally:
                os.close(descriptor)


if __name__ == "__main__":
    unittest.main()
