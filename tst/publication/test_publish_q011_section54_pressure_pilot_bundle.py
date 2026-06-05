#!/usr/bin/env python3
"""Focused adversarial tests for Q-011 pressure-pilot raw-bundle publication."""

from __future__ import annotations

import base64
from contextlib import contextmanager
import errno
import fcntl
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import shutil
import stat
import tarfile
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
        metadata = pilot_fixture._approved_snapshot_metadata(case_id, index)
        observed_time = float(metadata["observed_time_omega0_inverse"])
        cycle = int(metadata["cycle"])
        _put(
            root,
            sources["mhd_w_bcc"],
            pilot_fixture._binary(
                observed_time,
                argv_value,
                pilot_fixture._ATHENAK_MHD_W_BCC_FIELDS,
                cycle=cycle,
            ),
        )
        _put(root, sources["bmag"], pilot_fixture._binary(observed_time, argv_value, ("bmag",), cycle=cycle))
        _put(root, sources["prtcl_jx"], pilot_fixture._binary(observed_time, argv_value, ("prtcl_jx",), cycle=cycle))
        _put(root, sources["j2"], pilot_fixture._binary(observed_time, argv_value, ("j2",), cycle=cycle))
        _put(root, sources["prtcl_all"], pilot_fixture._particle_vtk(observed_time, cycle=cycle))
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


def _rewrite_publication_seal(base: Path, receipt: Path) -> None:
    publication_descriptor = os.open(base, os.O_RDONLY | os.O_DIRECTORY)
    try:
        seal = (
            base.parent
            / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
            / publisher._publication_seal_name(receipt.name)
        )
        seal.chmod(0o600)
        seal.write_bytes(
            publisher._publication_seal_payload(
                publication_descriptor,
                receipt.name,
                _sha256(receipt.read_bytes()),
                publisher._file_identity_at(
                    publication_descriptor, receipt.name, "test receipt"
                ),
            )
        )
        seal.chmod(0o444)
    finally:
        os.close(publication_descriptor)


def _write_forged_publication_namespace_seal(base: Path, receipt: Path) -> None:
    publication_descriptor = os.open(base, os.O_RDONLY | os.O_DIRECTORY)
    try:
        seal = base / publisher._publication_seal_name(receipt.name)
        seal.write_bytes(
            publisher._publication_seal_payload(
                publication_descriptor,
                receipt.name,
                _sha256(receipt.read_bytes()),
                publisher._file_identity_at(
                    publication_descriptor, receipt.name, "test receipt"
                ),
            )
        )
        seal.chmod(0o444)
    finally:
        os.close(publication_descriptor)


@contextmanager
def _reviewed_source_archive() -> Iterator[tuple[Path, str]]:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "source.tar"
        static = publisher._static_source_bindings()
        required = {
            record["path"]: record["sha256"]
            for record in static["reviewed_source_closure"]
        }
        authorization = static["postrun_aggregate_source_authorization"]
        required[authorization["path"]] = authorization["sha256"]
        with tarfile.open(path, mode="w", pax_headers={"comment": "a" * 40}) as archive:
            for relative in sorted(required):
                payload = (publisher.pilot.REPO_ROOT / relative).read_bytes()
                member = tarfile.TarInfo(relative)
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
        path.chmod(0o444)
        yield path, _sha256(path.read_bytes())


@contextmanager
def _worker_source_snapshot(archive: Path) -> Iterator[Path]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory) / "snapshot"
        root.mkdir()
        with tarfile.open(archive, mode="r:") as source:
            source.extractall(root)
        _freeze_existing(root)
        with patch.object(publisher, "EXECUTING_SOURCE_ROOT", root), patch.dict(
            os.environ,
            {
                "PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH": str(archive),
                publisher.WORKER_SOURCE_SNAPSHOT_ROOT_ENV: str(root),
            },
            clear=True,
        ):
            yield root
        _make_writable(root)


@contextmanager
def _verified_raw_cases() -> Iterator[tuple[Path, dict[str, Path], dict[str, str]]]:
    with tempfile.TemporaryDirectory() as directory:
        pic_root = Path(directory) / "pic"
        base = pic_root / "publication"
        runs = pic_root / "runs"
        base.mkdir(parents=True)
        (pic_root / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY).mkdir()
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
    def test_public_receipt_verifier_exposes_no_marker_bypass_flags(self) -> None:
        self.assertEqual(
            list(inspect.signature(publisher.verify_published_pressure_pilot_receipt).parameters),
            ["receipt_path", "authorized_pic_root"],
        )

    def test_guard_disarm_rejects_replacement_inode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            receipt_name = "receipt.json"
            guard_name = publisher._publication_guard_name(receipt_name)
            moved = root / "moved-guard"
            try:
                identity = publisher._arm_publication_guard_at(descriptor, receipt_name)
                (root / guard_name).rename(moved)
                (root / guard_name).write_bytes(publisher._PUBLICATION_GUARD_PAYLOAD)
                (root / guard_name).chmod(0o444)
                with self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError,
                    "changed during publication",
                ):
                    publisher._disarm_publication_guard_at(
                        descriptor, receipt_name, identity
                    )
                self.assertTrue((root / guard_name).is_file())
                self.assertTrue(moved.is_file())
            finally:
                os.close(descriptor)

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
                    "postrun_aggregate_source_authorization",
                    "registered_execution_preregistration",
                    "historical_v2_execution_preregistration",
                    "reviewed_source_closure",
                    "runtime_source_archive",
                },
            )
            self.assertEqual(
                publication["source_bindings"]["runtime_source_archive"],
                {
                    "execution_mode": "direct_api_nonproduction_only",
                    "git_commit": None,
                    "archive_sha256": None,
                    "verified_source_closure_sha256": None,
                },
            )
            self.assertEqual(publication["consumption_rule"], publisher.CONSUMPTION_RULE)
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
            success_seal = (
                base.parent
                / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                / publisher._publication_seal_name(receipt_path.name)
            )
            self.assertTrue(success_seal.is_file())
            self.assertFalse(success_seal.stat().st_mode & 0o222)
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

    def test_direct_bundle_verifier_rejects_production_namespace(self) -> None:
        production_root = Path(os.path.abspath(publisher.AUTHORIZED_PUBLICATION_ROOT))
        drifted_analyzer_root = production_root.parent / "drifted-analyzer-authority"

        def canonical(path: Path, _label: str) -> Path:
            return Path(os.path.abspath(path))

        with patch.object(
            pilot, "AUTHORIZED_PUBLICATION_ROOT", drifted_analyzer_root
        ), patch.object(
            publisher, "_canonical_existing_directory", side_effect=canonical
        ), patch.object(
            publisher, "_open_absolute_directory"
        ) as open_directory, self.assertRaisesRegex(
            publisher.PressurePilotPublicationError,
            "requires receipt consumption",
        ):
            publisher.verify_published_pressure_pilot_bundle(
                production_root / "must-not-open",
                "0" * 64,
                authorized_publication_root=production_root,
            )
        open_directory.assert_not_called()

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
            receipt["source_bindings"]["reviewed_source_closure"][0][
                "sha256"
            ] = "0" * 64
            receipt_path.write_bytes(publisher._canonical_json_bytes(receipt))
            receipt_path.chmod(0o444)
            _rewrite_publication_seal(base, receipt_path)
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "source binding drifted: reviewed_source_closure",
            ):
                publisher.verify_published_pressure_pilot_receipt(
                    receipt_path, authorized_pic_root=base.parent
                )

    def test_mocked_worker_runtime_archive_binding_is_retained_and_strict(self) -> None:
        with _reviewed_source_archive() as (archive, digest), _worker_source_snapshot(
            archive
        ), patch.object(
            publisher,
            "_trusted_source_archive",
            return_value=("a" * 40, archive.read_bytes()),
        ):
            bindings = publisher._source_bindings(require_verified_archive=True)
            binding = bindings["runtime_source_archive"]
        self.assertEqual(
            binding,
            {
                "execution_mode": "worker_extracted_git_archive_head_verified",
                "git_commit": "a" * 40,
                "archive_sha256": digest,
                "verified_source_closure_sha256": publisher._source_closure_sha256(
                    {key: value for key, value in bindings.items() if key != "runtime_source_archive"}
                ),
            },
        )
        with patch.dict(
            os.environ,
            {
                "PIC_PRESSURE_PUBLICATION_SOURCE_COMMIT": "a" * 40,
                "PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_SHA256": "b" * 64,
            },
            clear=True,
        ), self.assertRaisesRegex(
            publisher.PressurePilotPublicationError,
            "production publication requires a verified worker source archive",
        ):
            publisher._source_bindings(require_verified_archive=True)
        with patch.dict(os.environ, {}, clear=True):
            direct = publisher._source_bindings()
        with self.assertRaisesRegex(
            publisher.PressurePilotPublicationError,
            "direct API archive binding drifted",
        ):
            publisher._validate_retained_source_bindings(
                direct, require_verified_archive=True
            )
        with patch.dict(
            os.environ, {"PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH": "/absent/archive"}
        ), self.assertRaisesRegex(
            OSError, "No such file or directory"
        ):
            publisher._source_bindings(require_verified_archive=True)

    def test_synthetic_pax_archive_forgery_is_rejected(self) -> None:
        with _reviewed_source_archive() as (archive, _digest), _worker_source_snapshot(
            archive
        ), self.assertRaisesRegex(
            publisher.PressurePilotPublicationError,
            "differs from trusted repository HEAD archive",
        ):
            publisher._source_bindings(require_verified_archive=True)

    def test_production_api_rejects_direct_trusted_checkout_execution(self) -> None:
        with _reviewed_source_archive() as (archive, _digest), patch.dict(
            os.environ,
            {
                "PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH": str(archive),
                publisher.WORKER_SOURCE_SNAPSHOT_ROOT_ENV: str(
                    publisher.EXECUTING_SOURCE_ROOT
                ),
            },
            clear=True,
        ), patch.object(
            publisher,
            "_trusted_source_archive",
            return_value=("a" * 40, archive.read_bytes()),
        ), self.assertRaisesRegex(
            publisher.PressurePilotPublicationError,
            "must not execute from the trusted checkout",
        ):
            publisher._source_bindings(require_verified_archive=True)

    def test_receipt_is_the_only_accepted_consumption_marker(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            original = publisher._require_same_directory_at
            observed_pre_receipt_window = False

            def inspect_pre_receipt_window(
                parent_descriptor: int, name: str, descriptor: int, label: str
            ) -> None:
                nonlocal observed_pre_receipt_window
                original(parent_descriptor, name, descriptor, label)
                if label != "pressure-pilot public build tree" or observed_pre_receipt_window:
                    return
                observed_pre_receipt_window = True
                self.assertTrue(bundle.is_dir())
                self.assertFalse(receipt.exists())
                self.assertTrue(
                    (base / publisher._publication_guard_name(receipt.name)).is_file()
                )
                with self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError, "fail-closed guard"
                ):
                    publisher.consume_published_pressure_pilot_bundle(
                        receipt, authorized_pic_root=base.parent
                    )

            with patch.object(
                publisher,
                "_require_same_directory_at",
                side_effect=inspect_pre_receipt_window,
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_pre_receipt_window)
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )

    def test_production_publication_without_measured_archive_fails_before_exposure(
        self,
    ) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            with patch.object(
                publisher, "AUTHORIZED_PIC_ROOT", base.parent
            ), patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "production publication requires a verified worker source archive",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertFalse((base / "bundle").exists())
            self.assertFalse((base / "publication-receipt.json").exists())
            self.assertFalse((base / "aggregate-analysis.json").exists())

    def test_analyzer_and_publisher_require_registered_execution_digest(self) -> None:
        with _verified_raw_cases() as (base, roots, digests), patch.object(
            pilot, "REGISTERED_EXECUTION_PREREGISTRATION_SHA256", "0" * 64
        ):
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "registered-execution preregistration SHA-256 drifted",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
        with pilot_fixture._fixture() as fixture, patch.object(
            pilot, "REGISTERED_EXECUTION_PREREGISTRATION_SHA256", "0" * 64
        ):
            with self.assertRaisesRegex(
                pilot.PilotAnalysisError,
                "registered-execution preregistration SHA-256 drifted",
            ):
                pilot.analyze_pressure_pilot_bundle(
                    fixture[0],
                    fixture[2],
                    authorized_publication_root=fixture[0].parent,
                )

    def test_late_aggregate_failure_retains_guarded_public_artifacts(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            original = publisher._verify_pressure_pilot_bundle_at
            calls = 0

            def fail_after_all_public_renames(*args: object, **kwargs: object) -> object:
                nonlocal calls
                calls += 1
                if calls == 3:
                    raise publisher.PressurePilotPublicationError(
                        "injected late aggregate failure"
                    )
                return original(*args, **kwargs)

            with patch.object(
                publisher,
                "_verify_pressure_pilot_bundle_at",
                side_effect=fail_after_all_public_renames,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(bundle.is_dir())
            self.assertTrue(receipt.is_file())
            self.assertTrue(result.is_file())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            self.assertFalse(
                any(".staging-" in path.name for path in base.iterdir())
            )

    def test_late_publication_root_clone_is_rejected_by_retained_identity(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            anchor = base.parent / (base.name + ".descriptor-anchor")
            original = publisher._read_stable_readonly_regular_at
            cloned = False

            def clone_root(
                parent_descriptor: int, name: str, label: str
            ) -> bytes:
                nonlocal cloned
                payload = original(parent_descriptor, name, label)
                if label == "pressure-pilot aggregate analysis result" and not cloned:
                    cloned = True
                    os.rename(base, anchor)
                    shutil.copytree(anchor, base)
                return payload

            with patch.object(
                publisher,
                "_read_stable_readonly_regular_at",
                side_effect=clone_root,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(cloned)
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "fail-closed guard",
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_failure_does_not_delete_substituted_public_bundle(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            moved = base / "bundle.moved-original"
            original = publisher._verify_pressure_pilot_bundle_at
            calls = 0

            def substitute_before_failure(
                *args: object, **kwargs: object
            ) -> object:
                nonlocal calls
                calls += 1
                if calls == 3:
                    os.rename(bundle, moved)
                    shutil.copytree(moved, bundle)
                    raise publisher.PressurePilotPublicationError(
                        "injected post-substitution aggregate failure"
                    )
                return original(*args, **kwargs)

            with patch.object(
                publisher,
                "_verify_pressure_pilot_bundle_at",
                side_effect=substitute_before_failure,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(bundle.is_dir())
            self.assertTrue(moved.is_dir())
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_rejects_receipt_substitution_before_publication(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            moved = base / "publication-receipt.json.moved-original"
            original = publisher._publish_publication_seal_at

            def substitute_before_seal(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> None:
                os.rename(receipt, moved)
                shutil.copy2(moved, receipt)
                return original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=substitute_before_seal,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            self.assertFalse(
                (
                    base.parent
                    / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                    / publisher._publication_seal_name(receipt.name)
                ).exists()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_coordinated_replacement_failure_leaves_receipt_invalidated(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            original = publisher._verify_pressure_pilot_bundle_at
            calls = 0

            def replace_every_public_artifact(
                *args: object, **kwargs: object
            ) -> object:
                nonlocal calls
                calls += 1
                if calls == 3:
                    for path in (bundle, receipt, result):
                        moved = base / f"{path.name}.moved-original"
                        os.rename(path, moved)
                        if moved.is_dir():
                            shutil.copytree(moved, path)
                        else:
                            shutil.copy2(moved, path)
                    raise publisher.PressurePilotPublicationError(
                        "injected coordinated aggregate replacement failure"
                    )
                return original(*args, **kwargs)

            with patch.object(
                publisher,
                "_verify_pressure_pilot_bundle_at",
                side_effect=replace_every_public_artifact,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_rename_commit_is_reconciled_after_wrapper_raise(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            original = publisher._rename_no_replace_at
            observed_committed_seal = False

            def raise_after_seal_commit(
                parent_descriptor: int, source_name: str, destination_name: str
            ) -> None:
                nonlocal observed_committed_seal
                original(parent_descriptor, source_name, destination_name)
                if (
                    destination_name
                    != publisher._publication_seal_name(receipt.name)
                    or observed_committed_seal
                ):
                    return
                observed_committed_seal = True
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )
                raise publisher.PressurePilotPublicationError(
                    "injected post-commit seal-wrapper failure"
                )

            with patch.object(
                publisher, "_rename_no_replace_at", side_effect=raise_after_seal_commit
            ):
                published = publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_committed_seal)
            self.assertEqual(published["receipt_path"], str(receipt))
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )

    def test_success_seal_helper_wrapper_failure_rearms_guard(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            original = publisher._publish_publication_seal_at
            observed_committed_seal = False

            def raise_after_seal_helper_commit(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> None:
                nonlocal observed_committed_seal
                original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )
                observed_committed_seal = True
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )
                raise publisher.PressurePilotPublicationError(
                    "injected post-commit seal-helper failure"
                )

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=raise_after_seal_helper_commit,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_committed_seal)
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_is_committed_while_guard_remains_armed(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            original = publisher._publish_publication_seal_at
            observed_guard = False

            def inspect_guard_before_seal(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> tuple[int, int]:
                nonlocal observed_guard
                publisher._file_identity_at(
                    publication_descriptor,
                    publisher._publication_guard_name(receipt_name),
                    "receipt publication guard",
                )
                observed_guard = True
                return original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=inspect_guard_before_seal,
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_guard)
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )

    def test_guard_disarm_commit_is_reconciled_after_wrapper_raise(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            original = publisher._disarm_publication_guard_at

            def raise_after_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                original(parent_descriptor, receipt_name, guard_identity)
                raise publisher.PressurePilotPublicationError(
                    "injected post-commit guard-disarm wrapper failure"
                )

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=raise_after_disarm,
            ):
                published = publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertEqual(published["receipt_path"], str(receipt))
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )

    def test_guard_unlink_without_parent_sync_is_reconciled_durably(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            guard_unlinked = False
            observed_retry_sync = False
            publication_descriptor: int | None = None
            original_fsync = publisher._fsync_descriptor
            original_require_same_file = publisher._require_same_file_at
            original_require_guard_absent = publisher._require_publication_guard_absent_at

            def unlink_then_raise(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                nonlocal guard_unlinked, publication_descriptor
                publication_descriptor = parent_descriptor
                guard_name = publisher._publication_guard_name(receipt_name)
                publisher._require_same_file_at(
                    parent_descriptor,
                    guard_name,
                    guard_identity,
                    "receipt publication guard",
                )
                os.unlink(guard_name, dir_fd=parent_descriptor)
                guard_unlinked = True
                raise OSError("injected failure before parent sync")

            def record_fsync(descriptor: int) -> None:
                nonlocal observed_retry_sync
                if guard_unlinked and descriptor == publication_descriptor:
                    observed_retry_sync = True
                original_fsync(descriptor)

            def require_same_file_after_sync(
                parent_descriptor: int,
                name: str,
                identity: tuple[int, int],
                label: str,
            ) -> None:
                if guard_unlinked:
                    self.assertTrue(observed_retry_sync)
                original_require_same_file(parent_descriptor, name, identity, label)

            def require_guard_absent_after_sync(
                parent_descriptor: int, receipt_name: str, label: str
            ) -> None:
                if guard_unlinked:
                    self.assertTrue(observed_retry_sync)
                original_require_guard_absent(parent_descriptor, receipt_name, label)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=unlink_then_raise,
            ), patch.object(
                publisher,
                "_fsync_descriptor",
                side_effect=record_fsync,
            ), patch.object(
                publisher,
                "_require_same_file_at",
                side_effect=require_same_file_after_sync,
            ), patch.object(
                publisher,
                "_require_publication_guard_absent_at",
                side_effect=require_guard_absent_after_sync,
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_retry_sync)
            self.assertIsNotNone(publication_descriptor)
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )

    def test_publication_holds_transaction_lock_through_guard_disarm(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            acceptance = base.parent / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
            anchor = publisher._publication_transaction_anchor(base.parent)
            observed_lock = False
            original = publisher._disarm_publication_guard_at

            def inspect_lock_before_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                nonlocal observed_lock
                for path in (anchor, acceptance):
                    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
                    try:
                        with self.assertRaises(BlockingIOError):
                            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    finally:
                        os.close(descriptor)
                observed_lock = True
                original(parent_descriptor, receipt_name, guard_identity)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=inspect_lock_before_disarm,
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_lock)

    def test_descriptor_close_failure_does_not_retain_transaction_lock(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            publication = root / "publication"
            transaction_anchor = root / "transaction-anchor"
            publication.mkdir()
            transaction_anchor.mkdir()
            publication_descriptor = os.open(
                publication, os.O_RDONLY | os.O_DIRECTORY
            )
            transaction_descriptor = os.open(
                transaction_anchor, os.O_RDONLY | os.O_DIRECTORY
            )
            fcntl.flock(publication_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(transaction_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            original_close = publisher.os.close
            injected_error = OSError(errno.EIO, os.strerror(errno.EIO))

            def close_then_raise(descriptor: int) -> None:
                original_close(descriptor)
                if descriptor == publication_descriptor:
                    raise injected_error

            with patch.object(publisher.os, "close", side_effect=close_then_raise):
                observed_error = publisher._close_descriptors(
                    (publication_descriptor, transaction_descriptor)
                )

            self.assertIs(observed_error, injected_error)
            for path in (publication, transaction_anchor):
                descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                finally:
                    os.close(descriptor)

    def test_success_seal_rejects_byte_identical_different_inode_collision(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            acceptance = base.parent / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
            seal_name = publisher._publication_seal_name(receipt.name)
            original = publisher._rename_no_replace_at

            def collide_with_byte_identical_seal(
                parent_descriptor: int, source_name: str, destination_name: str
            ) -> None:
                if destination_name != seal_name:
                    original(parent_descriptor, source_name, destination_name)
                    return
                payload = publisher._read_stable_readonly_regular_at(
                    parent_descriptor,
                    source_name,
                    "receipt staged durable success seal",
                )
                descriptor = os.open(
                    destination_name,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                    0o600,
                    dir_fd=parent_descriptor,
                )
                try:
                    os.write(descriptor, payload)
                    os.fsync(descriptor)
                    os.fchmod(descriptor, 0o444)
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
                os.unlink(source_name, dir_fd=parent_descriptor)
                raise publisher.PressurePilotPublicationError(
                    "injected byte-identical different-inode seal collision"
                )

            with patch.object(
                publisher,
                "_rename_no_replace_at",
                side_effect=collide_with_byte_identical_seal,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue((acceptance / seal_name).is_file())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_post_disarm_receipt_substitution_fails_closed(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            moved = base / "publication-receipt.json.moved-original"
            original = publisher._disarm_publication_guard_at

            def substitute_after_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                original(parent_descriptor, receipt_name, guard_identity)
                os.rename(receipt, moved)
                shutil.copy2(moved, receipt)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=substitute_after_disarm,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=receipt,
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                publisher.consume_published_pressure_pilot_bundle(
                    receipt, authorized_pic_root=base.parent
                )

    def test_retained_receipt_schema_version_rejects_json_boolean_alias(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            publisher.publish_pressure_pilot_bundle(
                base / "bundle",
                receipt_path=receipt,
                analysis_result_path=base / "aggregate-analysis.json",
                case_artifact_dirs=roots,
                case_descriptor_sha256=digests,
                authorized_pic_root=base.parent,
            )
            receipt.chmod(0o600)
            record = json.loads(receipt.read_text(encoding="utf-8"))
            record["schema_version"] = True
            receipt.write_bytes(publisher._canonical_json_bytes(record))
            receipt.chmod(0o444)
            _rewrite_publication_seal(base, receipt)
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "receipt identity drifted",
            ):
                publisher.verify_published_pressure_pilot_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_retained_receipt_requires_durable_success_seal(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            receipt = base / "publication-receipt.json"
            publisher.publish_pressure_pilot_bundle(
                base / "bundle",
                receipt_path=receipt,
                analysis_result_path=base / "aggregate-analysis.json",
                case_artifact_dirs=roots,
                case_descriptor_sha256=digests,
                authorized_pic_root=base.parent,
            )
            (
                base.parent
                / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                / publisher._publication_seal_name(receipt.name)
            ).unlink()
            _write_forged_publication_namespace_seal(base, receipt)
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "durable success seal"
            ):
                publisher.verify_published_pressure_pilot_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_public_build_substitution_is_not_followed_or_deleted(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            bundle = base / "bundle"
            receipt = base / "publication-receipt.json"
            result = base / "aggregate-analysis.json"
            original = publisher._require_same_directory_at
            observed_build_names = []

            def substitute_public_build(
                parent_descriptor: int, name: str, descriptor: int, label: str
            ) -> None:
                if label == "pressure-pilot public build tree" and not observed_build_names:
                    observed_build_names.append(name)
                    os.rename(
                        name,
                        name + ".descriptor-anchor",
                        src_dir_fd=parent_descriptor,
                        dst_dir_fd=parent_descriptor,
                    )
                    os.mkdir(name, dir_fd=parent_descriptor)
                original(parent_descriptor, name, descriptor, label)

            with patch.object(
                publisher,
                "_require_same_directory_at",
                side_effect=substitute_public_build,
            ), self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "reviewed reconciliation required",
            ):
                publisher.publish_pressure_pilot_bundle(
                    bundle,
                    receipt_path=receipt,
                    analysis_result_path=result,
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(bundle.is_dir())
            self.assertFalse(receipt.exists())
            self.assertFalse(result.exists())
            self.assertEqual(len(observed_build_names), 1)
            substituted = base / observed_build_names[0]
            anchored = base / (observed_build_names[0] + ".descriptor-anchor")
            self.assertTrue(substituted.is_dir())
            self.assertTrue(anchored.is_dir())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            _make_writable(anchored)

    def test_lustre_file_rename_fallback_uses_atomic_link_and_preserves_identity(
        self,
    ) -> None:
        class UnsupportedRenameat2:
            argtypes: object = None
            restype: object = None

            def __call__(self, *_args: object) -> int:
                return -1

        class UnsupportedLibc:
            renameat2 = UnsupportedRenameat2()

        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            identity = os.stat(parent / "source").st_ino
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=UnsupportedLibc()
                ), patch.object(
                    publisher.ctypes, "get_errno", return_value=errno.EINVAL
                ), patch.object(publisher.os, "rename") as rename:
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                rename.assert_not_called()
                self.assertFalse((parent / "source").exists())
                self.assertEqual(os.stat(parent / "destination").st_ino, identity)
                self.assertEqual(
                    (parent / "destination").read_text(encoding="utf-8"),
                    "retain source\n",
                )
            finally:
                os.close(descriptor)

    def test_lustre_file_rename_fallback_rejects_collision(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            (parent / "destination").write_text("retain collision\n", encoding="utf-8")
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError, "collided"
                ):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertEqual(
                    (parent / "source").read_text(encoding="utf-8"),
                    "retain source\n",
                )
                self.assertEqual(
                    (parent / "destination").read_text(encoding="utf-8"),
                    "retain collision\n",
                )
            finally:
                os.close(descriptor)

    def test_lustre_file_rename_fallback_reconciles_post_link_wrapper_failure(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            identity = os.stat(parent / "source").st_ino
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            original_link = publisher.os.link

            def link_then_raise(*args: object, **kwargs: object) -> None:
                original_link(*args, **kwargs)
                raise OSError(errno.EIO, os.strerror(errno.EIO))

            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), patch.object(publisher.os, "link", side_effect=link_then_raise):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertFalse((parent / "source").exists())
                self.assertEqual(os.stat(parent / "destination").st_ino, identity)
            finally:
                os.close(descriptor)

    def test_lustre_file_rename_fallback_reconciles_post_unlink_wrapper_failure(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            identity = os.stat(parent / "source").st_ino
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            original_unlink = publisher.os.unlink

            def unlink_then_raise(path: object, *args: object, **kwargs: object) -> None:
                original_unlink(path, *args, **kwargs)
                if path == "source":
                    raise OSError(errno.EIO, os.strerror(errno.EIO))

            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), patch.object(publisher.os, "unlink", side_effect=unlink_then_raise):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertFalse((parent / "source").exists())
                destination = os.stat(parent / "destination")
                self.assertEqual(destination.st_ino, identity)
                self.assertEqual(destination.st_nlink, 1)
            finally:
                os.close(descriptor)

    def test_lustre_file_rename_fallback_reconciles_terminal_check_wrapper_failure(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            identity = os.stat(parent / "source").st_ino
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            original_require_absent = publisher._require_absent_at
            injected = False

            def require_absent_then_raise(
                parent_descriptor: int, name: str, label: str
            ) -> None:
                nonlocal injected
                original_require_absent(parent_descriptor, name, label)
                if name == "source" and not injected:
                    injected = True
                    raise RuntimeError("injected terminal-check wrapper failure")

            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), patch.object(
                    publisher,
                    "_require_absent_at",
                    side_effect=require_absent_then_raise,
                ):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertTrue(injected)
                self.assertFalse((parent / "source").exists())
                destination = os.stat(parent / "destination")
                self.assertEqual(destination.st_ino, identity)
                self.assertEqual(destination.st_nlink, 1)
            finally:
                os.close(descriptor)

    def test_lustre_file_fallback_failure_retains_canonical_and_staging_links(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            source = parent / "source"
            destination = parent / "destination"
            source.write_text("retain for reconciliation\n", encoding="utf-8")
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            original_unlink = publisher.os.unlink

            def reject_source_unlink(
                path: object, *args: object, **kwargs: object
            ) -> None:
                if path == "source":
                    raise OSError(errno.EIO, os.strerror(errno.EIO))
                original_unlink(path, *args, **kwargs)

            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), patch.object(
                    publisher.os, "unlink", side_effect=reject_source_unlink
                ), self.assertRaises(OSError):
                    publisher._rename_no_replace_at(
                        descriptor, source.name, destination.name
                    )
                self.assertTrue(source.is_file())
                self.assertTrue(destination.is_file())
                self.assertEqual(source.stat().st_ino, destination.stat().st_ino)
                self.assertEqual(source.stat().st_nlink, 2)
            finally:
                os.close(descriptor)

    def test_retained_reader_rejects_interrupted_two_link_commit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            source = parent / "source"
            destination = parent / "destination"
            source.write_text("sealed\n", encoding="utf-8")
            source.chmod(0o444)
            os.link(source, destination)
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError, "singly linked"
                ):
                    publisher._read_stable_readonly_regular_at(
                        descriptor, destination.name, "interrupted commit"
                    )
            finally:
                os.close(descriptor)

    def test_lustre_directory_rename_fallback_is_rejected_without_overwrite(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            source = parent / "source"
            source.mkdir()
            (source / "member").write_text("sealed\n", encoding="utf-8")
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError,
                    "exclusive final-name creation",
                ):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertTrue(source.is_dir())
                self.assertFalse((parent / "destination").exists())
            finally:
                os.close(descriptor)

    def test_lustre_file_fallback_publishes_complete_bundle_and_seal(self) -> None:
        with _verified_raw_cases() as (base, roots, digests), patch.object(
            publisher.ctypes, "CDLL", return_value=object()
        ):
            receipt = base / "publication-receipt.json"
            publisher.publish_pressure_pilot_bundle(
                base / "bundle",
                receipt_path=receipt,
                analysis_result_path=base / "aggregate-analysis.json",
                case_artifact_dirs=roots,
                case_descriptor_sha256=digests,
                authorized_pic_root=base.parent,
            )
            publisher.consume_published_pressure_pilot_bundle(
                receipt, authorized_pic_root=base.parent
            )
            self.assertFalse(any(".staging-" in path.name for path in base.iterdir()))

    def test_lustre_file_fallback_runs_under_transaction_lock(self) -> None:
        with _verified_raw_cases() as (base, roots, digests):
            acceptance = base.parent / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
            anchor = publisher._publication_transaction_anchor(base.parent)
            observed_locked_commits = 0
            original = publisher._link_no_replace_file_at

            def inspect_lock_during_link_commit(
                parent_descriptor: int,
                source_name: str,
                destination_name: str,
                identity: tuple[int, int],
            ) -> None:
                nonlocal observed_locked_commits
                for path in (anchor, acceptance):
                    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
                    try:
                        with self.assertRaises(BlockingIOError):
                            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    finally:
                        os.close(descriptor)
                observed_locked_commits += 1
                original(
                    parent_descriptor,
                    source_name,
                    destination_name,
                    identity,
                )

            with patch.object(
                publisher.ctypes, "CDLL", return_value=object()
            ), patch.object(
                publisher,
                "_link_no_replace_file_at",
                side_effect=inspect_lock_during_link_commit,
            ):
                publisher.publish_pressure_pilot_bundle(
                    base / "bundle",
                    receipt_path=base / "publication-receipt.json",
                    analysis_result_path=base / "aggregate-analysis.json",
                    case_artifact_dirs=roots,
                    case_descriptor_sha256=digests,
                    authorized_pic_root=base.parent,
                )
            self.assertGreaterEqual(observed_locked_commits, 3)

    def test_lustre_file_rename_fallback_requires_isolated_parent(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").write_text("retain source\n", encoding="utf-8")
            os.chmod(parent, 0o770)
            descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with patch.object(
                    publisher.ctypes, "CDLL", return_value=object()
                ), self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError,
                    "same-account isolated parent",
                ):
                    publisher._rename_no_replace_at(
                        descriptor, "source", "destination"
                    )
                self.assertTrue((parent / "source").is_file())
                self.assertFalse((parent / "destination").exists())
            finally:
                os.close(descriptor)


if __name__ == "__main__":
    unittest.main()
