#!/usr/bin/env python3
"""Focused adversarial tests for the read-only Q011 review-packet verifier."""

from __future__ import annotations

import ast
import builtins
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Iterator
import types
import unittest
from unittest.mock import patch

from tst.publication.frontier_control_plane import q011_pressure_review_packet_verifier as verifier


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _source_closure_sha256(records: list[dict[str, str]]) -> str:
    return _sha256(
        json.dumps(
            records,
            separators=(",", ":"),
            sort_keys=True,
            allow_nan=False,
        ).encode("utf-8")
    )


def _json_clone(value: object) -> object:
    return json.loads(json.dumps(value))


_AUTHORIZED_ACTIVE_DECK_BINDING = {
    "path": "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
    "sha256": "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1",
}
_AUTHORIZED_COMMON_OVERRIDES = (
    "mesh/nx1=100",
    "mesh/x1max=1200",
    "mesh/nx2=20",
    "mesh/x2max=240",
    "mesh_refinement/refinement=none",
    "mesh_refinement/num_levels=1",
    "time/tlim=60",
    "time/nlim=4096",
    "time/ndiag=50",
    "problem/ps_enable_curvature_amr=false",
    "problem/ps_feedback_diag_dcycle=50",
    "output1/variable=mhd_w_bcc",
    "output1/id=mhd_w_bcc",
    "output1/dt=15",
    "output2/dt=15",
    "output3/dt=15",
    "output4/dt=15",
    "output5/dt=15",
    "output6/dt=15",
)
_AUTHORIZED_RUNTIME_PROFILE = "frontier_minimum_supported"
_AUTHORIZED_PARALLEL_RANKS = 1
_AUTHORIZED_ROCR_VISIBLE_DEVICE = 0


def _make_tree_writable(root: Path) -> None:
    for directory, names, filenames in os.walk(root):
        base = Path(directory)
        base.chmod(0o755)
        for name in names:
            (base / name).chmod(0o755)
        for name in filenames:
            (base / name).chmod(0o644)


class _PublishedPacket:
    def __init__(self) -> None:
        self._temporary = tempfile.TemporaryDirectory(
            prefix="q011-pressure-review-packet-verifier-"
        )
        self.root = Path(self._temporary.name) / "PIC"
        self.publication = self.root / "publication"
        self.acceptance = self.root / "publication_acceptance"
        self.runs = self.root / "runs"
        self.packet_root = self.publication / "review-packet"
        self.aggregate_bundle = self.publication / "aggregate-bundle"
        self.aggregate_analysis = self.publication / "aggregate-analysis.json"
        self.aggregate_receipt = self.publication / "aggregate-receipt.json"
        self.packet_receipt = self.publication / "review-packet-receipt.json"
        self.publication.mkdir(parents=True)
        self.acceptance.mkdir()
        self.runs.mkdir()
        self._build()

    def close(self) -> None:
        if self.root.exists():
            for path in sorted(
                self.root.rglob("*"),
                key=lambda candidate: len(candidate.parts),
                reverse=True,
            ):
                if path.is_symlink():
                    path.unlink()
                elif path.is_dir():
                    path.chmod(0o755)
                else:
                    path.chmod(0o644)
            self.root.chmod(0o755)
        self._temporary.cleanup()

    def _source_bindings(self) -> dict[str, object]:
        return {
            "postrun_aggregate_source_authorization": {
                "path": "tst/publication/readiness/source-authorization.json",
                "sha256": "1" * 64,
            },
            "registered_execution_preregistration": {
                "path": "tst/publication/readiness/registered-execution.json",
                "sha256": "2" * 64,
            },
            "historical_v2_execution_preregistration": {
                "path": "tst/publication/readiness/historical-v2.json",
                "sha256": "3" * 64,
            },
            "reviewed_source_closure": [
                {
                    "role": "aggregate_publisher",
                    "path": "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
                    "sha256": "4" * 64,
                },
                {
                    "role": "review_packet_renderer",
                    "path": "tst/publication/render_q011_section54_pressure_pilot_review_packet.py",
                    "sha256": "5" * 64,
                },
            ],
            "runtime_source_archive": {
                "execution_mode": "direct_api_nonproduction_only",
                "git_commit": None,
                "archive_sha256": None,
                "verified_source_closure_sha256": None,
            },
        }

    def _write_readonly(self, path: Path, payload: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        path.chmod(0o444)

    def _write_canonical_readonly(self, path: Path, value: object) -> bytes:
        payload = _canonical_json_bytes(value)
        self._write_readonly(path, payload)
        return payload

    def _freeze_tree(self, root: Path) -> None:
        for path in sorted(root.rglob("*"), key=lambda candidate: len(candidate.parts), reverse=True):
            path.chmod(0o555 if path.is_dir() else 0o444)
        root.chmod(0o555)

    def _raw_cases(
        self,
        aggregate_bundle: Path,
    ) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        cases: list[dict[str, object]] = []
        manifest_cases: list[dict[str, object]] = []
        for index, (case_id, pressure, argv_value) in enumerate(
            zip(
                verifier.RAW_CASE_IDS,
                verifier.RAW_CASE_PRESSURES,
                verifier.RAW_CASE_ARGV_VALUES,
            )
        ):
            artifact_dir = self.runs / f"raw-{case_id}" / "attempt"
            artifact_dir.mkdir(parents=True)
            raw_payloads: dict[str, bytes] = {}
            aggregate_payloads: dict[str, bytes] = {}
            raw_inventory: dict[str, dict[str, object]] = {}
            bundle_members: list[dict[str, object]] = []

            def bind_raw(source: str, payload: bytes) -> None:
                if source in raw_payloads:
                    self.assert_payload_equal(raw_payloads[source], payload, source)
                    return
                raw_payloads[source] = payload
                raw_inventory[source] = {
                    "path": source,
                    "sha256": _sha256(payload),
                    "size": len(payload),
                }

            def bind_bundle(source: str, target: str, payload: bytes) -> dict[str, str]:
                bind_raw(source, payload)
                aggregate_payloads[target] = payload
                bundle_members.append(
                    {
                        "path": target,
                        "sha256": _sha256(payload),
                        "size": len(payload),
                        "source_path": source,
                    }
                )
                return {"path": target, "sha256": _sha256(payload)}

            snapshots: list[dict[str, object]] = []
            for snapshot_index, time in enumerate(verifier.RAW_CASE_TIMES):
                suffix = f"{snapshot_index:05d}"
                snapshot: dict[str, object] = {"time": time}
                for kind, directory, extension in (
                    ("mhd_w_bcc", "bin", "bin"),
                    ("bmag", "bin", "bin"),
                    ("prtcl_jx", "bin", "bin"),
                    ("j2", "bin", "bin"),
                    ("prtcl_all", "pvtk", "part.vtk"),
                ):
                    target = f"cases/{case_id}/{directory}/{case_id}.{kind}.{suffix}.{extension}"
                    source = f"output/{directory}/{case_id}.{kind}.{suffix}.{extension}"
                    snapshot[kind] = bind_bundle(
                        source,
                        target,
                        f"{case_id}:{kind}:{suffix}\n".encode("ascii"),
                    )
                snapshots.append(snapshot)

            stdout_payload = f"{case_id}: stdout\n".encode("ascii")
            stdout = bind_bundle(
                "athena_stdout.txt",
                f"cases/{case_id}/stdout.txt",
                stdout_payload,
            )
            restart_prefix = f"{case_id}.00004.rst"
            restart_bindings = {}
            for name, suffix in (
                ("manifest", ".manifest"),
                ("manifest_complete", ".manifest.complete"),
                ("artifact", ""),
                ("complete", ".complete"),
            ):
                filename = restart_prefix + suffix
                restart_bindings[name] = bind_bundle(
                    f"output/rst/{filename}",
                    f"cases/{case_id}/rst/{filename}",
                    f"{case_id}:{name}\n".encode("ascii"),
                )
            manifest_case = {
                "case_id": case_id,
                "ps_p0": pressure,
                "overrides": [
                    *_AUTHORIZED_COMMON_OVERRIDES,
                    f"problem/ps_p0={argv_value}",
                ],
                "snapshots": snapshots,
                "stdout": stdout,
                "terminal_restart": {
                    "time": verifier.RAW_CASE_TIMES[-1],
                    "manifest": restart_bindings["manifest"],
                    "manifest_complete": restart_bindings["manifest_complete"],
                    "members": [
                        {
                            "artifact": restart_bindings["artifact"],
                            "complete": restart_bindings["complete"],
                        }
                    ],
                },
            }
            manifest_cases.append(manifest_case)

            allowlist = f"q011-pressure-{case_id.replace('_', '-')}.environment.allowlist.txt"
            runtime_payloads = {
                allowlist: b"PIC_FRONTIER_PROFILE=frontier_minimum_supported\n",
                "athena_stdout.txt": stdout_payload,
                "athena_stdout.sha256": (_sha256(stdout_payload) + "\n").encode("ascii"),
                "athena_stderr.txt": b"frontier diagnostic stderr\n",
            }
            for path, payload in runtime_payloads.items():
                bind_raw(path, payload)
            runtime_artifacts = {
                path: _sha256(payload) for path, payload in runtime_payloads.items()
            }
            inventory_payload = self._write_canonical_readonly(
                artifact_dir / verifier.RAW_CASE_INVENTORY_NAME,
                {
                    "schema_version": 1,
                    "files": [raw_inventory[path] for path in sorted(raw_inventory)],
                },
            )
            descriptor = {
                "schema_version": 1,
                "record_type": verifier.RAW_CASE_RECORD_TYPE,
                "evidence_class": verifier.AGGREGATE_EVIDENCE_CLASS,
                "qualification_effect": verifier.AGGREGATE_QUALIFICATION_EFFECT,
                "launch_contract": "trusted_trampoline_athena_argv_v1",
                "case_id": case_id,
                "ps_p0": pressure,
                "argv_value": argv_value,
                "artifact_inventory_sha256": _sha256(inventory_payload),
                "runtime_artifacts": runtime_artifacts,
                "runtime_profile": _AUTHORIZED_RUNTIME_PROFILE,
                "parallel_ranks": _AUTHORIZED_PARALLEL_RANKS,
                "rank_gpu_bindings": [
                    {
                        "host": "frontier00001",
                        "rank": 0,
                        "rocr_visible_device": _AUTHORIZED_ROCR_VISIBLE_DEVICE,
                    }
                ],
                "manifest_case": manifest_case,
                "bundle_members": sorted(bundle_members, key=lambda member: member["path"]),
            }
            descriptor_payload = self._write_canonical_readonly(
                artifact_dir / verifier.RAW_CASE_DESCRIPTOR_PATH,
                descriptor,
            )
            for path, payload in raw_payloads.items():
                self._write_readonly(artifact_dir / path, payload)
            for path, payload in aggregate_payloads.items():
                self._write_readonly(aggregate_bundle / path, payload)
            self._freeze_tree(artifact_dir)
            (artifact_dir / "analysis").chmod(0o700)
            cases.append(
                {
                    "case_id": case_id,
                    "artifact_dir": str(artifact_dir),
                    "descriptor_path": verifier.RAW_CASE_DESCRIPTOR_PATH,
                    "descriptor_sha256": _sha256(descriptor_payload),
                    "artifact_inventory_sha256": _sha256(inventory_payload),
                    "runtime_artifacts": runtime_artifacts,
                }
            )
        return cases, manifest_cases

    def assert_payload_equal(self, left: bytes, right: bytes, label: str) -> None:
        if left != right:
            raise AssertionError(f"fixture payload collision: {label}")

    def _publication_identity(self) -> dict[str, int]:
        metadata = self.publication.stat()
        return {"device": metadata.st_dev, "inode": metadata.st_ino}

    def _seal_record(self, receipt: Path) -> dict[str, object]:
        payload = receipt.read_bytes()
        metadata = receipt.stat()
        return {
            "schema_version": 1,
            "record_type": verifier.SUCCESS_SEAL_RECORD_TYPE,
            "publication_root_identity": self._publication_identity(),
            "receipt_name": receipt.name,
            "receipt_sha256": _sha256(payload),
            "receipt_identity": {
                "device": metadata.st_dev,
                "inode": metadata.st_ino,
            },
        }

    def seal_path(self, receipt: Path) -> Path:
        return self.acceptance / f".{receipt.name}.publication-success"

    def guard_path(self, receipt: Path) -> Path:
        return self.publication / f".{receipt.name}.publication-invalid"

    def _replace_seal(self, receipt: Path) -> None:
        seal = self.seal_path(receipt)
        if seal.exists():
            seal.chmod(0o644)
            seal.unlink()
        self._write_canonical_readonly(seal, self._seal_record(receipt))

    def _build(self) -> None:
        self.aggregate_bundle.mkdir()
        self._write_canonical_readonly(self.aggregate_analysis, {"status": "pass"})
        source_bindings = self._source_bindings()
        raw_cases, manifest_cases = self._raw_cases(self.aggregate_bundle)
        aggregate_manifest = {
            "schema_version": 1,
            "record_type": verifier.AGGREGATE_MANIFEST_RECORD_TYPE,
            "evidence_class": verifier.AGGREGATE_EVIDENCE_CLASS,
            "qualification_effect": verifier.AGGREGATE_QUALIFICATION_EFFECT,
            "active_deck_binding": dict(_AUTHORIZED_ACTIVE_DECK_BINDING),
            "preregistration_binding": source_bindings[
                "postrun_aggregate_source_authorization"
            ],
            "registered_execution_preregistration_binding": source_bindings[
                "registered_execution_preregistration"
            ],
            "cases": manifest_cases,
        }
        aggregate_manifest_payload = self._write_canonical_readonly(
            self.aggregate_bundle / verifier.AGGREGATE_MANIFEST_NAME,
            aggregate_manifest,
        )
        self._freeze_tree(self.aggregate_bundle)

        aggregate_record = {
            "schema_version": 1,
            "record_type": verifier.AGGREGATE_RECEIPT_RECORD_TYPE,
            "evidence_class": verifier.AGGREGATE_EVIDENCE_CLASS,
            "qualification_effect": verifier.AGGREGATE_QUALIFICATION_EFFECT,
            "consumption_rule": verifier.CONSUMPTION_RULE,
            "publication_root_identity": self._publication_identity(),
            "aggregate_bundle": {
                "path": str(self.aggregate_bundle),
                "manifest_sha256": _sha256(aggregate_manifest_payload),
            },
            "aggregate_analysis": {
                "path": str(self.aggregate_analysis),
                "sha256": _sha256(self.aggregate_analysis.read_bytes()),
            },
            "source_bindings": source_bindings,
            "raw_cases": raw_cases,
        }
        aggregate_payload = self._write_canonical_readonly(
            self.aggregate_receipt,
            aggregate_record,
        )
        self.aggregate_binding = {
            "path": str(self.aggregate_receipt),
            "sha256": _sha256(aggregate_payload),
        }

        figures = self.packet_root / "figures"
        figures.mkdir(parents=True)
        payloads = {
            "PRESSURE_REVIEW_PACKET.md": b"engineering calibration review packet\n",
            "figures/terminal_mhd_pic_pressure_comparison.png": b"pressure png\n",
            "figures/terminal_profile_overlays.png": b"profile png\n",
            "pressure_review_metrics.json": _canonical_json_bytes(
                {
                    "schema_version": 1,
                    "record_type": verifier.REVIEW_METRICS_RECORD_TYPE,
                    "watermark": verifier.WATERMARK,
                    "qualification_effect": verifier.QUALIFICATION_EFFECT,
                    "aggregate_receipt": self.aggregate_binding,
                    "cases": [
                        {
                            "case_id": case_id,
                            "problem_ps_p0": float(problem_ps_p0),
                            "terminal_particle_count": 1,
                            "particle_efficiency": 1.0,
                            "zone_cycles_per_second": 1.0,
                            "particle_updates_per_second": 1.0,
                            "tracked_gpu_memory_high_water_bytes": 1,
                        }
                        for case_id, problem_ps_p0 in zip(
                            verifier.RAW_CASE_IDS,
                            verifier.RAW_CASE_PRESSURES,
                        )
                    ],
                }
            ),
        }
        for relative, payload in payloads.items():
            self._write_readonly(self.packet_root / relative, payload)
        inventory = {
            "schema_version": 1,
            "record_type": verifier.INVENTORY_RECORD_TYPE,
            "members": [
                {
                    "path": relative,
                    "sha256": _sha256(payloads[relative]),
                    "size": len(payloads[relative]),
                }
                for relative in sorted(payloads)
            ],
        }
        inventory_payload = self._write_canonical_readonly(
            self.packet_root / verifier.INVENTORY_NAME,
            inventory,
        )
        figures.chmod(0o555)
        self.packet_root.chmod(0o555)

        packet_record = {
            "schema_version": 1,
            "record_type": verifier.PACKET_RECEIPT_RECORD_TYPE,
            "watermark": verifier.WATERMARK,
            "qualification_effect": verifier.QUALIFICATION_EFFECT,
            "consumption_rule": verifier.CONSUMPTION_RULE,
            "publication_root_identity": self._publication_identity(),
            "aggregate_receipt": self.aggregate_binding,
            "packet_root": str(self.packet_root),
            "inventory_sha256": _sha256(inventory_payload),
            "source_bindings": source_bindings,
        }
        packet_payload = self._write_canonical_readonly(
            self.packet_receipt,
            packet_record,
        )
        self.packet_binding = {
            "path": str(self.packet_receipt),
            "sha256": _sha256(packet_payload),
        }
        self._replace_seal(self.aggregate_receipt)
        self._replace_seal(self.packet_receipt)

    def verify(self) -> dict[str, object]:
        return verifier.verify_pressure_review_packet_binding(
            self.packet_binding,
            self.aggregate_binding,
            authorized_pic_root=self.root,
        )

    def consume(
        self,
        receipt_path: str | os.PathLike[str] | None = None,
        aggregate_receipt_binding: object | None = None,
    ) -> dict[str, object]:
        return verifier.consume_published_pressure_pilot_review_packet(
            self.packet_receipt if receipt_path is None else receipt_path,
            aggregate_receipt_binding=(
                self.aggregate_binding
                if aggregate_receipt_binding is None
                else aggregate_receipt_binding
            ),
            authorized_pic_root=self.root,
        )

    def bind_worker_verified_source_archive(self) -> None:
        aggregate = json.loads(self.aggregate_receipt.read_text(encoding="utf-8"))
        source_bindings = aggregate["source_bindings"]
        closure_payload = _canonical_json_bytes(
            {
                "postrun_aggregate_source_authorization": source_bindings[
                    "postrun_aggregate_source_authorization"
                ],
                "reviewed_source_closure": source_bindings["reviewed_source_closure"],
            }
        )
        source_bindings["runtime_source_archive"] = {
            "execution_mode": "worker_extracted_git_archive_head_verified",
            "git_commit": "e" * 40,
            "archive_sha256": "f" * 64,
            "verified_source_closure_sha256": _sha256(closure_payload),
        }
        self.rewrite_receipt(
            self.aggregate_receipt,
            lambda value: value.__setitem__("source_bindings", source_bindings),
        )
        metrics_path = self.packet_root / "pressure_review_metrics.json"
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        metrics["aggregate_receipt"] = dict(self.aggregate_binding)
        self.rewrite_packet_member(
            "pressure_review_metrics.json",
            _canonical_json_bytes(metrics),
        )

        def bind_packet(value: dict[str, object]) -> None:
            value["aggregate_receipt"] = dict(self.aggregate_binding)
            value["source_bindings"] = source_bindings

        self.rewrite_receipt(self.packet_receipt, bind_packet)

    def rewrite_receipt(self, receipt: Path, mutate: object) -> None:
        receipt.chmod(0o644)
        value = json.loads(receipt.read_text(encoding="utf-8"))
        mutate(value)
        payload = _canonical_json_bytes(value)
        receipt.write_bytes(payload)
        receipt.chmod(0o444)
        if receipt == self.packet_receipt:
            self.packet_binding["sha256"] = _sha256(payload)
        elif receipt == self.aggregate_receipt:
            self.aggregate_binding["sha256"] = _sha256(payload)
        else:
            raise AssertionError("unknown receipt")
        self._replace_seal(receipt)

    def cross_bind_aggregate_receipt(self) -> None:
        self.rewrite_receipt(
            self.packet_receipt,
            lambda value: value.__setitem__(
                "aggregate_receipt",
                dict(self.aggregate_binding),
            ),
        )

    def rewrite_aggregate_manifest(self, mutate: object) -> None:
        manifest_path = self.aggregate_bundle / verifier.AGGREGATE_MANIFEST_NAME
        manifest_path.chmod(0o644)
        value = json.loads(manifest_path.read_text(encoding="utf-8"))
        mutate(value)
        payload = _canonical_json_bytes(value)
        manifest_path.write_bytes(payload)
        manifest_path.chmod(0o444)
        self.rewrite_receipt(
            self.aggregate_receipt,
            lambda receipt: receipt["aggregate_bundle"].__setitem__(
                "manifest_sha256",
                _sha256(payload),
            ),
        )
        self.cross_bind_aggregate_receipt()

    def rewrite_raw_descriptor(self, case_index: int, mutate: object) -> None:
        aggregate = json.loads(self.aggregate_receipt.read_text(encoding="utf-8"))
        raw_case = aggregate["raw_cases"][case_index]
        descriptor_path = Path(raw_case["artifact_dir"]) / raw_case["descriptor_path"]
        descriptor_path.chmod(0o644)
        value = json.loads(descriptor_path.read_text(encoding="utf-8"))
        mutate(value)
        payload = _canonical_json_bytes(value)
        descriptor_path.write_bytes(payload)
        descriptor_path.chmod(0o444)
        self.rewrite_receipt(
            self.aggregate_receipt,
            lambda receipt: receipt["raw_cases"][case_index].__setitem__(
                "descriptor_sha256",
                _sha256(payload),
            ),
        )
        self.cross_bind_aggregate_receipt()

    def rewrite_manifest_case(self, case_index: int, mutate: object) -> None:
        self.rewrite_aggregate_manifest(
            lambda manifest: mutate(manifest["cases"][case_index])
        )
        self.rewrite_raw_descriptor(
            case_index,
            lambda descriptor: mutate(descriptor["manifest_case"]),
        )

    def rewrite_raw_receipt(self, receipt: Path, payload: bytes) -> None:
        receipt.chmod(0o644)
        receipt.write_bytes(payload)
        receipt.chmod(0o444)
        if receipt == self.packet_receipt:
            self.packet_binding["sha256"] = _sha256(payload)
        elif receipt == self.aggregate_receipt:
            self.aggregate_binding["sha256"] = _sha256(payload)
        else:
            raise AssertionError("unknown receipt")
        self._replace_seal(receipt)

    def add_guard(self, receipt: Path) -> None:
        self._write_readonly(
            self.guard_path(receipt),
            b"receipt publication is not authoritative\n",
        )

    def remove_seal(self, receipt: Path) -> None:
        seal = self.seal_path(receipt)
        seal.chmod(0o644)
        seal.unlink()

    def tamper_seal_hash(self, receipt: Path) -> None:
        seal = self.seal_path(receipt)
        seal.chmod(0o644)
        value = json.loads(seal.read_text(encoding="utf-8"))
        value["receipt_sha256"] = "0" * 64
        seal.write_bytes(_canonical_json_bytes(value))
        seal.chmod(0o444)

    def add_packet_member(self, relative: str, payload: bytes) -> None:
        self.packet_root.chmod(0o755)
        self._write_readonly(self.packet_root / relative, payload)
        self.packet_root.chmod(0o555)

    def tamper_packet_member(self, relative: str, payload: bytes) -> None:
        self.tamper_readonly_file(self.packet_root / relative, payload)

    def rewrite_packet_member(self, relative: str, payload: bytes) -> None:
        self.tamper_packet_member(relative, payload)
        inventory_path = self.packet_root / verifier.INVENTORY_NAME
        inventory_path.chmod(0o644)
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        member = next(item for item in inventory["members"] if item["path"] == relative)
        member["sha256"] = _sha256(payload)
        member["size"] = len(payload)
        inventory_payload = _canonical_json_bytes(inventory)
        inventory_path.write_bytes(inventory_payload)
        inventory_path.chmod(0o444)
        self.rewrite_receipt(
            self.packet_receipt,
            lambda value: value.__setitem__(
                "inventory_sha256",
                _sha256(inventory_payload),
            ),
        )

    def tamper_readonly_file(self, path: Path, payload: bytes) -> None:
        path.chmod(0o644)
        path.write_bytes(payload)
        path.chmod(0o444)


@contextmanager
def _published_packet() -> Iterator[_PublishedPacket]:
    packet = _PublishedPacket()
    try:
        yield packet
    finally:
        packet.close()


class _SealedPressureGateAttestations:
    def __init__(self, packet: _PublishedPacket) -> None:
        self.packet = packet
        self.root = packet.root
        self.archive_root = self.root / verifier.PRESSURE_GATE_ATTESTATION_ROOT_NAME
        self.archive_root.mkdir()
        self.now = datetime(2026, 6, 5, 9, 14, 0, tzinfo=timezone.utc)
        self.operator_id = "operator.alpha"
        self.reviewer_id = "reviewer.alpha"
        self.recomputed_utc = "2026-06-05T09:10:00Z"
        self.reanalysis_sealed_utc = "2026-06-05T09:11:00Z"
        self.reviewed_utc = "2026-06-05T09:12:00Z"
        self.reviewer_sealed_utc = "2026-06-05T09:13:00Z"
        self.selected_case = {"case_id": "ps_p0_0p10", "problem_ps_p0": 0.1}
        self.git_commit = "a" * 40
        self.source_archive_sha256 = "b" * 64
        self.source_closure = [
            {"path": path, "sha256": f"{index + 1:064x}"}
            for index, path in enumerate(verifier.PRESSURE_REANALYSIS_SOURCE_PATHS)
        ]
        self.helper_source_closure = _json_clone(self.source_closure)
        aggregate = json.loads(packet.aggregate_receipt.read_text(encoding="utf-8"))
        self.manifest_sha256 = aggregate["aggregate_bundle"]["manifest_sha256"]
        self.analysis_sha256 = aggregate["aggregate_analysis"]["sha256"]
        self.result = {
            "packet_receipt_sha256": packet.packet_binding["sha256"],
            "aggregate_receipt_sha256": packet.aggregate_binding["sha256"],
            "manifest_sha256": self.manifest_sha256,
            "analysis_result_sha256": self.analysis_sha256,
            "status": "pass_engineering_calibration_only",
        }
        self.source_authorization = {
            "execution_mode": verifier.PRESSURE_REANALYSIS_EXECUTION_MODE,
            "git_commit": self.git_commit,
            "source_archive_sha256": self.source_archive_sha256,
            "source_closure_sha256": _source_closure_sha256(self.source_closure),
            "source_closure": _json_clone(self.source_closure),
            "historical_production_source_authorization": dict(
                verifier.AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION
            ),
        }
        reanalysis = {
            "schema_version": 1,
            "record_type": verifier.PRESSURE_REANALYSIS_RECORD_TYPE,
            "qualification_effect": verifier.PRESSURE_REANALYSIS_QUALIFICATION_EFFECT,
            "operator_id": self.operator_id,
            "recomputed_utc": self.recomputed_utc,
            "sealed_utc": self.reanalysis_sealed_utc,
            "operator_statement": verifier.PRESSURE_REANALYSIS_OPERATOR_STATEMENT,
            "evidence": {
                "published_pressure_pilot_receipt": dict(packet.aggregate_binding),
                "published_pressure_pilot_review_packet_receipt": dict(
                    packet.packet_binding
                ),
                "pilot_bundle_manifest_sha256": self.manifest_sha256,
                "aggregate_pilot_analysis_sha256": self.analysis_sha256,
            },
            "source_authorization": _json_clone(self.source_authorization),
            "result": dict(self.result),
        }
        self.reanalysis_path, self.reanalysis_binding = self._write_attestation(
            (
                "20260605T091100Z-q011-section54-pressure-reanalysis-"
                f"{self.operator_id}"
            ),
            reanalysis,
        )
        reviewer = {
            "schema_version": 1,
            "record_type": verifier.PRESSURE_REVIEWER_RECORD_TYPE,
            "qualification_effect": verifier.PRESSURE_REVIEWER_QUALIFICATION_EFFECT,
            "selection_method": "human_review_only",
            "reviewer_id": self.reviewer_id,
            "reviewed_utc": self.reviewed_utc,
            "sealed_utc": self.reviewer_sealed_utc,
            "rationale": "Selected the exact registered pressure case after review.",
            "reviewer_statement": verifier.PRESSURE_REVIEWER_STATEMENT,
            "published_pressure_pilot_receipt": dict(packet.aggregate_binding),
            "published_pressure_pilot_review_packet_receipt": dict(
                packet.packet_binding
            ),
            "authoritative_reanalysis_attestation": dict(self.reanalysis_binding),
            "selected_case": dict(self.selected_case),
        }
        self.reviewer_path, self.reviewer_binding = self._write_attestation(
            (
                "20260605T091300Z-q011-section54-pressure-selection-"
                f"{self.reviewer_id}"
            ),
            reviewer,
        )
        self.archive_root.chmod(0o555)

    def _write_attestation(
        self, directory_name: str, value: dict[str, object]
    ) -> tuple[Path, dict[str, str]]:
        directory = self.archive_root / directory_name
        directory.mkdir()
        path = directory / verifier.PRESSURE_GATE_ATTESTATION_FILENAME
        payload = _canonical_json_bytes(value)
        path.write_bytes(payload)
        path.chmod(0o400)
        directory.chmod(0o500)
        return path, {"path": str(path), "sha256": _sha256(payload)}

    def _attributes(self, kind: str) -> tuple[Path, dict[str, str]]:
        if kind == "reanalysis":
            return self.reanalysis_path, self.reanalysis_binding
        if kind == "reviewer":
            return self.reviewer_path, self.reviewer_binding
        raise AssertionError(f"unknown pressure-gate attestation kind: {kind}")

    def rewrite(
        self,
        kind: str,
        mutate: object,
        *,
        canonical: bool = True,
        rebind: bool = True,
    ) -> None:
        path, binding = self._attributes(kind)
        directory = path.parent
        directory.chmod(0o700)
        path.chmod(0o600)
        value = json.loads(path.read_text(encoding="utf-8"))
        mutate(value)
        payload = (
            _canonical_json_bytes(value)
            if canonical
            else json.dumps(value, sort_keys=True).encode("utf-8")
        )
        path.write_bytes(payload)
        path.chmod(0o400)
        directory.chmod(0o500)
        if rebind:
            binding["sha256"] = _sha256(payload)

    def rename_directory(self, kind: str, directory_name: str) -> None:
        path, binding = self._attributes(kind)
        self.archive_root.chmod(0o755)
        destination = self.archive_root / directory_name
        os.rename(path.parent, destination)
        self.archive_root.chmod(0o555)
        renamed = destination / verifier.PRESSURE_GATE_ATTESTATION_FILENAME
        binding["path"] = str(renamed)
        if kind == "reanalysis":
            self.reanalysis_path = renamed
        else:
            self.reviewer_path = renamed

    def add_member(self, kind: str) -> None:
        path, _ = self._attributes(kind)
        path.parent.chmod(0o700)
        extra = path.parent / "extra-member"
        extra.write_bytes(b"unexpected\n")
        extra.chmod(0o400)
        path.parent.chmod(0o500)

    def off_root_binding(self, kind: str) -> dict[str, str]:
        path, _ = self._attributes(kind)
        outside = self.root / f"off-root-{kind}-attestation.json"
        outside.write_bytes(path.read_bytes())
        outside.chmod(0o400)
        return {"path": str(outside), "sha256": _sha256(outside.read_bytes())}

    def consume_reanalysis(
        self, binding: object | None = None
    ) -> dict[str, object]:
        return verifier.consume_sealed_pressure_reanalysis_attestation(
            self.reanalysis_binding if binding is None else binding,
            aggregate_receipt_binding=self.packet.aggregate_binding,
            packet_receipt_binding=self.packet.packet_binding,
            pilot_bundle_manifest_sha256=self.manifest_sha256,
            aggregate_pilot_analysis_sha256=self.analysis_sha256,
            authorized_pic_root=self.root,
            expected_result=self.result,
            now=self.now,
        )

    def consume_reviewer(
        self,
        binding: object | None = None,
        *,
        reanalysis_verification: object | None = None,
    ) -> dict[str, object]:
        reanalysis = (
            self.consume_reanalysis()
            if reanalysis_verification is None
            else reanalysis_verification
        )
        return verifier.consume_sealed_pressure_reviewer_attestation(
            self.reviewer_binding if binding is None else binding,
            aggregate_receipt_binding=self.packet.aggregate_binding,
            packet_receipt_binding=self.packet.packet_binding,
            reanalysis_verification=reanalysis,
            selected_case=self.selected_case,
            authorized_pic_root=self.root,
            now=self.now,
        )


@contextmanager
def _sealed_pressure_gate(
    packet: _PublishedPacket,
) -> Iterator[_SealedPressureGateAttestations]:
    yield _SealedPressureGateAttestations(packet)


class Q011PressureReviewPacketVerifierTests(unittest.TestCase):
    def assert_rejected(self, packet: _PublishedPacket, pattern: str) -> None:
        with self.assertRaisesRegex(
            verifier.PressureReviewPacketVerificationError,
            pattern,
        ):
            packet.verify()

    def test_consumes_exact_sealed_pressure_gate_attestations_and_source_snapshot(
        self,
    ) -> None:
        with _published_packet() as packet, _sealed_pressure_gate(packet) as gate:
            reanalysis = gate.consume_reanalysis()
            self.assertEqual(reanalysis["binding"], gate.reanalysis_binding)
            self.assertEqual(reanalysis["operator_id"], gate.operator_id)
            self.assertEqual(reanalysis["result"], gate.result)
            self.assertEqual(
                reanalysis["source_authorization"],
                gate.source_authorization,
            )

            reviewer = gate.consume_reviewer(reanalysis_verification=reanalysis)
            self.assertEqual(reviewer["binding"], gate.reviewer_binding)
            self.assertEqual(reviewer["reviewer_id"], gate.reviewer_id)
            self.assertEqual(reviewer["selected_case"], gate.selected_case)
            self.assertEqual(
                reviewer["authoritative_reanalysis_attestation"],
                gate.reanalysis_binding,
            )

            verifier.validate_pressure_reanalysis_source_snapshot(
                reanalysis,
                git_commit=gate.git_commit,
                source_archive_sha256=gate.source_archive_sha256,
                helper_source_closure=gate.helper_source_closure,
                authorized_pic_root=gate.root,
            )

    def test_public_reanalysis_consumers_reopen_and_reject_forged_results(self) -> None:
        with _published_packet() as packet, _sealed_pressure_gate(packet) as gate:
            reanalysis = gate.consume_reanalysis()
            forged = _json_clone(reanalysis)
            forged["operator_id"] = "operator.forged"
            nonexistent = _json_clone(reanalysis)
            nonexistent["binding"] = {
                "path": str(
                    gate.archive_root
                    / (
                        "20260605T091100Z-q011-section54-pressure-reanalysis-"
                        "operator.missing"
                    )
                    / verifier.PRESSURE_GATE_ATTESTATION_FILENAME
                ),
                "sha256": "0" * 64,
            }
            altered_sealed_utc = _json_clone(reanalysis)
            altered_sealed_utc["sealed_utc"] = "2026-06-05T09:10:30Z"
            cases = (
                (
                    "forged dict",
                    forged,
                    "differs from re-consumed sealed attestation",
                ),
                ("nonexistent binding", nonexistent, "cannot be inspected"),
                (
                    "altered sealed_utc",
                    altered_sealed_utc,
                    "differs from re-consumed sealed attestation",
                ),
            )
            for name, candidate, pattern in cases:
                with self.subTest(consumer="reviewer", case=name), self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    gate.consume_reviewer(reanalysis_verification=candidate)
                with self.subTest(
                    consumer="source snapshot",
                    case=name,
                ), self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    verifier.validate_pressure_reanalysis_source_snapshot(
                        candidate,
                        git_commit=gate.git_commit,
                        source_archive_sha256=gate.source_archive_sha256,
                        helper_source_closure=gate.helper_source_closure,
                        authorized_pic_root=gate.root,
                    )
            with self.subTest(
                consumer="source snapshot",
                case="caller-supplied authorized root",
            ), self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "must be a direct publication child",
            ):
                verifier.validate_pressure_reanalysis_source_snapshot(
                    reanalysis,
                    git_commit=gate.git_commit,
                    source_archive_sha256=gate.source_archive_sha256,
                    helper_source_closure=gate.helper_source_closure,
                    authorized_pic_root=gate.root.parent,
                )

    def test_reanalysis_attestation_rejects_malformed_tampered_and_unsafe_trees(
        self,
    ) -> None:
        cases = (
            (
                "malformed",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value.pop("operator_statement"),
                ),
                "unexpected keys",
            ),
            (
                "tampered",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value.__setitem__("operator_id", "operator.tampered"),
                    rebind=False,
                ),
                "binding hash mismatch",
            ),
            (
                "wrong file mode",
                lambda gate: gate.reanalysis_path.chmod(0o444),
                "payload mode must be 0400",
            ),
            (
                "wrong directory mode",
                lambda gate: gate.reanalysis_path.parent.chmod(0o555),
                "directory mode must be 0500",
            ),
            (
                "extra member",
                lambda gate: gate.add_member("reanalysis"),
                "tree must contain only attestation.json",
            ),
            (
                "noncanonical JSON",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value,
                    canonical=False,
                ),
                "is not canonical JSON",
            ),
        )
        for name, prepare, pattern in cases:
            with self.subTest(name=name), _published_packet() as packet, _sealed_pressure_gate(
                packet
            ) as gate:
                prepare(gate)
                with self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    gate.consume_reanalysis()

        with self.subTest(name="off-root"), _published_packet() as packet, _sealed_pressure_gate(
            packet
        ) as gate:
            with self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "outside the fixed pressure-gate attestation layout",
            ):
                gate.consume_reanalysis(gate.off_root_binding("reanalysis"))

    def test_reanalysis_attestation_rejects_time_directory_mode_and_source_drift(
        self,
    ) -> None:
        def future(gate: _SealedPressureGateAttestations) -> None:
            gate.rewrite(
                "reanalysis",
                lambda value: value.__setitem__(
                    "sealed_utc",
                    "2026-06-05T09:20:00Z",
                ),
            )
            gate.rename_directory(
                "reanalysis",
                "20260605T092000Z-q011-section54-pressure-reanalysis-operator.alpha",
            )

        cases = (
            (
                "future",
                future,
                "timestamps are out of order",
            ),
            (
                "pre-gate",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value.__setitem__(
                        "recomputed_utc",
                        "2026-06-05T09:09:41Z",
                    ),
                ),
                "timestamps are out of order",
            ),
            (
                "timestamp order",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value.__setitem__(
                        "recomputed_utc",
                        "2026-06-05T09:12:00Z",
                    ),
                ),
                "timestamps are out of order",
            ),
            (
                "directory name",
                lambda gate: gate.rename_directory(
                    "reanalysis",
                    "20260605T091100Z-q011-section54-pressure-reanalysis-operator.wrong",
                ),
                "directory name differs",
            ),
            (
                "execution mode",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value["source_authorization"].__setitem__(
                        "execution_mode",
                        "direct_source_tree_reanalysis",
                    ),
                ),
                "execution_mode must equal",
            ),
            (
                "source closure path",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value["source_authorization"]["source_closure"][
                        0
                    ].__setitem__("path", "tst/publication/forged.py"),
                ),
                "source_closure\\[0\\].path drifted",
            ),
            (
                "source closure digest",
                lambda gate: gate.rewrite(
                    "reanalysis",
                    lambda value: value["source_authorization"].__setitem__(
                        "source_closure_sha256",
                        "0" * 64,
                    ),
                ),
                "source_closure_sha256 drifted",
            ),
        )
        for name, prepare, pattern in cases:
            with self.subTest(name=name), _published_packet() as packet, _sealed_pressure_gate(
                packet
            ) as gate:
                prepare(gate)
                with self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    gate.consume_reanalysis()

    def test_reviewer_attestation_rejects_time_identity_rationale_and_layout_drift(
        self,
    ) -> None:
        def future(gate: _SealedPressureGateAttestations) -> None:
            gate.rewrite(
                "reviewer",
                lambda value: value.__setitem__(
                    "sealed_utc",
                    "2026-06-05T09:20:00Z",
                ),
            )
            gate.rename_directory(
                "reviewer",
                "20260605T092000Z-q011-section54-pressure-selection-reviewer.alpha",
            )

        cases = (
            ("future", future, "timestamps are out of order"),
            (
                "pre-gate",
                lambda gate: gate.rewrite(
                    "reviewer",
                    lambda value: value.__setitem__(
                        "reviewed_utc",
                        "2026-06-05T09:09:41Z",
                    ),
                ),
                "timestamps are out of order",
            ),
            (
                "timestamp order",
                lambda gate: gate.rewrite(
                    "reviewer",
                    lambda value: value.__setitem__(
                        "reviewed_utc",
                        "2026-06-05T09:14:00Z",
                    ),
                ),
                "timestamps are out of order",
            ),
            (
                "directory name",
                lambda gate: gate.rename_directory(
                    "reviewer",
                    "20260605T091300Z-q011-section54-pressure-selection-reviewer.wrong",
                ),
                "directory name differs",
            ),
            (
                "reviewer ID",
                lambda gate: gate.rewrite(
                    "reviewer",
                    lambda value: value.__setitem__("reviewer_id", "Reviewer.Alpha"),
                ),
                "canonical lowercase reviewer/operator ID",
            ),
            (
                "short rationale",
                lambda gate: gate.rewrite(
                    "reviewer",
                    lambda value: value.__setitem__("rationale", "too short"),
                ),
                "trimmed and contain 16-1024 characters",
            ),
            (
                "multiline rationale",
                lambda gate: gate.rewrite(
                    "reviewer",
                    lambda value: value.__setitem__(
                        "rationale",
                        "Selected after review.\nInjected second line.",
                    ),
                ),
                "single-line printable ASCII",
            ),
        )
        for name, prepare, pattern in cases:
            with self.subTest(name=name), _published_packet() as packet, _sealed_pressure_gate(
                packet
            ) as gate:
                prepare(gate)
                with self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    gate.consume_reviewer()

        with self.subTest(name="off-root"), _published_packet() as packet, _sealed_pressure_gate(
            packet
        ) as gate:
            with self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "outside the fixed pressure-gate attestation layout",
            ):
                gate.consume_reviewer(gate.off_root_binding("reviewer"))

    def test_reanalysis_source_snapshot_rejects_frozen_snapshot_mismatch(self) -> None:
        with _published_packet() as packet, _sealed_pressure_gate(packet) as gate:
            reanalysis = gate.consume_reanalysis()
            mismatched_closure = _json_clone(gate.helper_source_closure)
            mismatched_closure[0]["sha256"] = "f" * 64
            missing_closure = _json_clone(gate.helper_source_closure)
            missing_closure.pop()
            cases = (
                ("Git commit", "c" * 40, gate.source_archive_sha256, gate.helper_source_closure),
                ("archive", gate.git_commit, "d" * 64, gate.helper_source_closure),
                ("source closure", gate.git_commit, gate.source_archive_sha256, mismatched_closure),
            )
            for name, git_commit, archive_sha256, closure in cases:
                with self.subTest(name=name), self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    "differs from the frozen qualifying-plan source snapshot",
                ):
                    verifier.validate_pressure_reanalysis_source_snapshot(
                        reanalysis,
                        git_commit=git_commit,
                        source_archive_sha256=archive_sha256,
                        helper_source_closure=closure,
                        authorized_pic_root=gate.root,
                    )
            with self.subTest(name="missing source"), self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "omits a reanalysis source",
            ):
                verifier.validate_pressure_reanalysis_source_snapshot(
                    reanalysis,
                    git_commit=gate.git_commit,
                    source_archive_sha256=gate.source_archive_sha256,
                    helper_source_closure=missing_closure,
                    authorized_pic_root=gate.root,
                )

    def test_accepts_exact_canonical_immutable_cross_bound_packet(self) -> None:
        with _published_packet() as packet:
            result = packet.verify()
            self.assertEqual(
                set(result),
                {
                    "packet_receipt",
                    "aggregate_receipt",
                    "packet_root",
                    "inventory_sha256",
                    "source_bindings",
                },
            )
            self.assertEqual(result["packet_receipt"], packet.packet_binding)
            self.assertEqual(result["aggregate_receipt"], packet.aggregate_binding)
            self.assertEqual(result["packet_root"], str(packet.packet_root))
            self.assertEqual(
                result["source_bindings"],
                json.loads(packet.packet_receipt.read_text(encoding="utf-8"))[
                    "source_bindings"
                ],
            )

    def test_stable_consumer_api_and_normalized_parsed_result(self) -> None:
        signature = inspect.signature(
            verifier.consume_published_pressure_pilot_review_packet
        )
        self.assertEqual(
            list(signature.parameters),
            ["receipt_path", "aggregate_receipt_binding", "authorized_pic_root"],
        )
        self.assertEqual(
            signature.parameters["receipt_path"].kind,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
        for name in ("aggregate_receipt_binding", "authorized_pic_root"):
            self.assertEqual(
                signature.parameters[name].kind,
                inspect.Parameter.KEYWORD_ONLY,
            )
        self.assertTrue(
            issubclass(verifier.PressureReviewPacketVerificationError, ValueError)
        )

        with _published_packet() as packet:
            noncanonical = f"{packet.publication}/./{packet.packet_receipt.name}"
            result = packet.consume(receipt_path=noncanonical)
            self.assertEqual(
                set(result),
                {
                    "receipt_binding",
                    "aggregate_receipt_binding",
                    "packet_receipt",
                    "aggregate_receipt",
                    "aggregate_bundle",
                    "aggregate_analysis",
                    "source_bindings",
                    "inventory",
                },
            )
            self.assertEqual(result["receipt_binding"], packet.packet_binding)
            self.assertEqual(
                result["aggregate_receipt_binding"],
                packet.aggregate_binding,
            )
            packet_record = json.loads(
                packet.packet_receipt.read_text(encoding="utf-8")
            )
            inventory = json.loads(
                (packet.packet_root / verifier.INVENTORY_NAME).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(result["packet_receipt"], packet_record)
            aggregate_record = json.loads(
                packet.aggregate_receipt.read_text(encoding="utf-8")
            )
            self.assertEqual(result["aggregate_receipt"], aggregate_record)
            self.assertEqual(result["aggregate_bundle"], aggregate_record["aggregate_bundle"])
            self.assertEqual(result["aggregate_analysis"], aggregate_record["aggregate_analysis"])
            self.assertEqual(result["source_bindings"], packet_record["source_bindings"])
            self.assertEqual(result["inventory"], inventory)
            result["source_bindings"]["reviewed_source_closure"][0]["sha256"] = "0" * 64
            self.assertNotEqual(
                result["source_bindings"],
                result["packet_receipt"]["source_bindings"],
            )
            result["aggregate_receipt"]["raw_cases"][0]["case_id"] = "detached-mutation"
            self.assertNotEqual(result["aggregate_receipt"], aggregate_record)

    def test_stable_consumer_raises_module_error_for_invalid_inputs(self) -> None:
        with _published_packet() as packet:
            drifted = dict(packet.aggregate_binding)
            drifted["sha256"] = "0" * 64
            with self.assertRaises(verifier.PressureReviewPacketVerificationError):
                packet.consume(aggregate_receipt_binding=drifted)
            with self.assertRaises(verifier.PressureReviewPacketVerificationError):
                packet.consume(receipt_path=packet.root / "outside.json")
            with self.assertRaises(verifier.PressureReviewPacketVerificationError):
                verifier.consume_published_pressure_pilot_review_packet(
                    object(),
                    aggregate_receipt_binding=packet.aggregate_binding,
                    authorized_pic_root=packet.root,
                )

    def test_rejects_empty_raw_cases(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_receipt(
                packet.aggregate_receipt,
                lambda value: value.__setitem__("raw_cases", []),
            )
            with self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "exact ordered four-case set",
            ):
                packet.consume()

    def test_rejects_raw_case_id_and_shape_drift(self) -> None:
        def forged_case_id(value: dict[str, object]) -> None:
            value["raw_cases"][0]["case_id"] = "ps_p0_forged"

        def missing_runtime_artifact(value: dict[str, object]) -> None:
            del value["raw_cases"][0]["runtime_artifacts"]["athena_stderr.txt"]

        for mutation, pattern in (
            (forged_case_id, "exact ordered four-case set"),
            (missing_runtime_artifact, "runtime_artifacts has unexpected members"),
        ):
            with self.subTest(pattern=pattern), _published_packet() as packet:
                packet.rewrite_receipt(packet.aggregate_receipt, mutation)
                with self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    pattern,
                ):
                    packet.consume()

    def test_rejects_unauthorized_active_deck_with_rebound_receipts(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_aggregate_manifest(
                lambda manifest: manifest.__setitem__(
                    "active_deck_binding",
                    {
                        "path": "inputs/publication/unauthorized.athinput",
                        "sha256": "0" * 64,
                    },
                )
            )
            self.assert_rejected(packet, "active deck is not the authorized")

    def test_rejects_injected_override_with_rebound_case_contract(self) -> None:
        def inject_extra_override(case: dict[str, object]) -> None:
            case["overrides"].insert(
                len(_AUTHORIZED_COMMON_OVERRIDES),
                "problem/unauthorized_injected_override=true",
            )

        with _published_packet() as packet:
            packet.rewrite_manifest_case(0, inject_extra_override)
            self.assert_rejected(packet, "overrides differs from the authorized")

    def test_rejects_runtime_profile_drift_with_rebound_descriptor(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_raw_descriptor(
                0,
                lambda descriptor: descriptor.__setitem__(
                    "runtime_profile",
                    "frontier_unauthorized_profile",
                ),
            )
            self.assert_rejected(packet, "runtime_profile")

    def test_rejects_rank_count_drift_with_rebound_descriptor(self) -> None:
        def drift_rank_count(descriptor: dict[str, object]) -> None:
            descriptor["parallel_ranks"] = 2
            descriptor["rank_gpu_bindings"].append(
                {
                    "host": "frontier00002",
                    "rank": 1,
                    "rocr_visible_device": _AUTHORIZED_ROCR_VISIBLE_DEVICE,
                }
            )

        with _published_packet() as packet:
            packet.rewrite_raw_descriptor(0, drift_rank_count)
            self.assert_rejected(packet, "parallel_ranks")

    def test_rejects_device_drift_with_rebound_descriptor(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_raw_descriptor(
                0,
                lambda descriptor: descriptor["rank_gpu_bindings"][0].__setitem__(
                    "rocr_visible_device",
                    1,
                ),
            )
            self.assert_rejected(packet, "rocr_visible_device")

    def test_production_root_requires_worker_verified_source_archive(self) -> None:
        with _published_packet() as packet, patch.object(
            verifier,
            "AUTHORIZED_PRODUCTION_PIC_ROOT",
            packet.root,
        ), patch.object(
            verifier,
            "AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING",
            dict(packet.packet_binding),
        ), patch.object(
            verifier,
            "AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING",
            dict(packet.aggregate_binding),
        ):
            with self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "production binding must be worker verified",
            ):
                packet.consume()
            packet.bind_worker_verified_source_archive()
            source_bindings = json.loads(
                packet.aggregate_receipt.read_text(encoding="utf-8")
            )["source_bindings"]
            with patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_SOURCE_BINDINGS",
                source_bindings,
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING",
                dict(packet.packet_binding),
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING",
                dict(packet.aggregate_binding),
            ):
                result = packet.consume()
            self.assertEqual(result["receipt_binding"], packet.packet_binding)

    def test_production_root_requires_exact_immutable_receipt_pair(self) -> None:
        with _published_packet() as packet:
            packet.bind_worker_verified_source_archive()
            source_bindings = json.loads(
                packet.aggregate_receipt.read_text(encoding="utf-8")
            )["source_bindings"]
            with patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_PIC_ROOT",
                packet.root,
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_SOURCE_BINDINGS",
                source_bindings,
            ):
                self.assert_rejected(
                    packet,
                    "not the authorized immutable production receipt",
                )

    def test_isolated_installed_context_has_no_repository_dependency(self) -> None:
        source = Path(verifier.__file__).read_bytes()
        imported: list[str] = []
        original_import = builtins.__import__

        def isolated_import(
            name: str,
            globals: object = None,
            locals: object = None,
            fromlist: object = (),
            level: int = 0,
        ) -> object:
            imported.append(name)
            if name == "tst" or name.startswith("tst.") or name.startswith(
                ("publish_q011_", "analyze_q011_", "render_q011_")
            ):
                raise ImportError(f"repository module is unavailable: {name}")
            return original_import(name, globals, locals, fromlist, level)

        installed_builtins = dict(vars(builtins))
        installed_builtins["__import__"] = isolated_import
        module = types.ModuleType("q011_pressure_review_packet_verifier")
        module.__file__ = "/installed/control_plane/q011_pressure_review_packet_verifier.py"
        module.__dict__["__builtins__"] = installed_builtins
        exec(
            compile(source, module.__file__, "exec", dont_inherit=True),
            module.__dict__,
        )
        with _published_packet() as packet:
            result = module.consume_published_pressure_pilot_review_packet(
                packet.packet_receipt,
                aggregate_receipt_binding=packet.aggregate_binding,
                authorized_pic_root=packet.root,
            )
        self.assertEqual(result["receipt_binding"], packet.packet_binding)
        self.assertFalse(
            any(
                name == "tst"
                or name.startswith("tst.")
                or name.startswith(("publish_q011_", "analyze_q011_", "render_q011_"))
                for name in imported
            )
        )

    def test_verifier_source_uses_no_mutation_or_publication_apis(self) -> None:
        tree = ast.parse(inspect.getsource(verifier))
        forbidden_calls = {
            "chmod",
            "fsync",
            "link",
            "mkdir",
            "remove",
            "rename",
            "replace",
            "rmdir",
            "symlink",
            "unlink",
            "write",
            "write_bytes",
            "write_text",
        }

        def is_datetime_parse_replace(node: ast.Call) -> bool:
            receiver = node.func.value
            return (
                node.func.attr == "replace"
                and isinstance(receiver, ast.Call)
                and isinstance(receiver.func, ast.Attribute)
                and receiver.func.attr == "strptime"
                and isinstance(receiver.func.value, ast.Name)
                and receiver.func.value.id == "datetime"
            )

        observed = {
            node.func.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and not is_datetime_parse_replace(node)
        }
        self.assertFalse(observed & forbidden_calls)
        forbidden_open_flags = {
            "O_APPEND",
            "O_CREAT",
            "O_EXCL",
            "O_RDWR",
            "O_TRUNC",
            "O_WRONLY",
        }
        observed_os_attributes = {
            node.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id == "os"
        }
        self.assertFalse(observed_os_attributes & forbidden_open_flags)
        imports = {
            alias.name.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        } | {
            (node.module or "").split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module != "__future__"
        }
        self.assertTrue(imports <= set(sys.stdlib_module_names))

    def test_rejects_packet_and_aggregate_guards(self) -> None:
        for receipt_name in ("packet_receipt", "aggregate_receipt"):
            with self.subTest(receipt=receipt_name), _published_packet() as packet:
                packet.add_guard(getattr(packet, receipt_name))
                self.assert_rejected(packet, "publication guard is present")

    def test_rejects_missing_packet_and_aggregate_success_seals(self) -> None:
        for receipt_name in ("packet_receipt", "aggregate_receipt"):
            with self.subTest(receipt=receipt_name), _published_packet() as packet:
                packet.remove_seal(getattr(packet, receipt_name))
                self.assert_rejected(packet, "success seal cannot be opened")

    def test_rejects_success_seal_that_does_not_bind_receipt_hash(self) -> None:
        with _published_packet() as packet:
            packet.tamper_seal_hash(packet.packet_receipt)
            self.assert_rejected(packet, "success seal does not bind")

    def test_rejects_aggregate_cross_binding_drift(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_receipt(
                packet.packet_receipt,
                lambda value: value["aggregate_receipt"].__setitem__("sha256", "0" * 64),
            )
            self.assert_rejected(packet, "does not bind the supplied aggregate receipt")

    def test_rejects_source_binding_drift_across_receipts(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_receipt(
                packet.packet_receipt,
                lambda value: value["source_bindings"][
                    "reviewed_source_closure"
                ][0].__setitem__("sha256", "0" * 64),
            )
            self.assert_rejected(packet, "source bindings differ")

    def test_production_rejects_cross_bound_but_unauthorized_source_tuple(self) -> None:
        with _published_packet() as packet:
            packet.bind_worker_verified_source_archive()
            original = json.loads(
                packet.aggregate_receipt.read_text(encoding="utf-8")
            )["source_bindings"]
            forged = json.loads(json.dumps(original))
            forged["registered_execution_preregistration"]["sha256"] = "0" * 64
            packet.rewrite_receipt(
                packet.aggregate_receipt,
                lambda value: value.__setitem__("source_bindings", forged),
            )

            def cross_bind_forged(value: dict[str, object]) -> None:
                value["aggregate_receipt"] = dict(packet.aggregate_binding)
                value["source_bindings"] = forged

            packet.rewrite_receipt(packet.packet_receipt, cross_bind_forged)
            with patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_PIC_ROOT",
                packet.root,
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_SOURCE_BINDINGS",
                original,
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING",
                dict(packet.packet_binding),
            ), patch.object(
                verifier,
                "AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING",
                dict(packet.aggregate_binding),
            ):
                self.assert_rejected(packet, "authorized immutable production source tuple")

    def test_rejects_cross_bound_malformed_source_bindings(self) -> None:
        with _published_packet() as packet:
            packet.rewrite_receipt(
                packet.aggregate_receipt,
                lambda value: value["source_bindings"][
                    "postrun_aggregate_source_authorization"
                ].__setitem__("path", "../escaped.json"),
            )

            def cross_bind_malformed_source(value: dict[str, object]) -> None:
                value["aggregate_receipt"] = dict(packet.aggregate_binding)
                value["source_bindings"][
                    "postrun_aggregate_source_authorization"
                ]["path"] = "../escaped.json"

            packet.rewrite_receipt(
                packet.packet_receipt,
                cross_bind_malformed_source,
            )
            self.assert_rejected(packet, "is unsafe")

    def test_rejects_aggregate_bundle_analysis_and_raw_tree_drift(self) -> None:
        with _published_packet() as packet:
            member = next(
                path
                for path in packet.aggregate_bundle.rglob("*")
                if path.is_file() and path.name != verifier.AGGREGATE_MANIFEST_NAME
            )
            packet.tamper_readonly_file(member, b"forged aggregate member\n")
            self.assert_rejected(packet, "aggregate bundle member hash drifted")
        with _published_packet() as packet:
            packet.tamper_readonly_file(packet.aggregate_analysis, b'{"status":"forged"}\n')
            self.assert_rejected(packet, "aggregate analysis hash")
        with _published_packet() as packet:
            aggregate = json.loads(packet.aggregate_receipt.read_text(encoding="utf-8"))
            raw_root = Path(aggregate["raw_cases"][0]["artifact_dir"])
            raw_member = next((raw_root / "output" / "bin").iterdir())
            packet.tamper_readonly_file(raw_member, b"forged raw source\n")
            self.assert_rejected(packet, "raw case .* member hash drifted")

    def test_rejects_packet_member_drift_and_extra_member(self) -> None:
        with _published_packet() as packet:
            packet.tamper_packet_member(
                "PRESSURE_REVIEW_PACKET.md",
                b"forged review packet\n",
            )
            self.assert_rejected(packet, "does not match the packet member")
        with _published_packet() as packet:
            packet.add_packet_member("undeclared.txt", b"undeclared\n")
            self.assert_rejected(packet, "closure is not exact")

    def test_rejects_coherently_rehashed_review_metrics_drift(self) -> None:
        with _published_packet() as packet:
            metrics_path = packet.packet_root / "pressure_review_metrics.json"
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            metrics["aggregate_receipt"]["sha256"] = "0" * 64
            packet.rewrite_packet_member(
                "pressure_review_metrics.json",
                _canonical_json_bytes(metrics),
            )
            self.assert_rejected(
                packet,
                "metrics does not bind the supplied aggregate receipt",
            )

        with _published_packet() as packet:
            metrics_path = packet.packet_root / "pressure_review_metrics.json"
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            metrics["cases"][0]["problem_ps_p0"] = 0.2
            packet.rewrite_packet_member(
                "pressure_review_metrics.json",
                _canonical_json_bytes(metrics),
            )
            self.assert_rejected(packet, "differs from registered pressure")

    def test_resource_limits_fail_closed(self) -> None:
        limits = (
            ("MAX_RETAINED_FILE_BYTES", 1, "retained-file size limit"),
            ("MAX_JSON_BYTES", 1, "JSON size limit"),
            ("MAX_DIRECTORY_ENTRIES", 1, "entry limit"),
            ("MAX_TREE_ENTRIES", 1, "tree-entry limit"),
            ("MAX_TREE_DEPTH", 1, "(?:path|tree)-depth limit"),
        )
        for name, value, pattern in limits:
            with self.subTest(limit=name), _published_packet() as packet, patch.object(
                verifier,
                name,
                value,
            ):
                self.assert_rejected(packet, pattern)

    def test_directory_entry_limit_does_not_use_unbounded_listdir(self) -> None:
        with _published_packet() as packet, patch.object(
            verifier.os,
            "listdir",
            side_effect=AssertionError("unbounded listdir must not be used"),
        ):
            self.assertEqual(packet.verify()["packet_receipt"], packet.packet_binding)

    def test_deeply_nested_json_is_normalized_to_module_error(self) -> None:
        with _published_packet() as packet:
            deeply_nested = b'{"nested":' + (b"[" * 2000) + b"0" + (b"]" * 2000) + b"}\n"
            packet.rewrite_raw_receipt(packet.packet_receipt, deeply_nested)
            self.assert_rejected(packet, "valid strict JSON")

    def test_rejects_noncanonical_or_writable_packet_receipt(self) -> None:
        with _published_packet() as packet:
            value = json.loads(packet.packet_receipt.read_text(encoding="utf-8"))
            packet.rewrite_raw_receipt(
                packet.packet_receipt,
                json.dumps(value, sort_keys=True).encode("utf-8"),
            )
            self.assert_rejected(packet, "is not canonical JSON")
        with _published_packet() as packet:
            packet.packet_receipt.chmod(0o644)
            self.assert_rejected(packet, "must be read-only")

    def test_nonblocking_open_rejects_fifo_special_file(self) -> None:
        if not getattr(os, "O_NONBLOCK", 0):
            self.skipTest("platform does not expose O_NONBLOCK")
        with _published_packet() as packet:
            packet.packet_receipt.chmod(0o644)
            packet.packet_receipt.unlink()
            os.mkfifo(packet.packet_receipt, 0o444)
            packet.packet_receipt.chmod(0o444)
            original_open = os.open
            observed = False

            def audited_open(
                path: str | bytes | os.PathLike[str] | os.PathLike[bytes],
                flags: int,
                *args: object,
                **kwargs: object,
            ) -> int:
                nonlocal observed
                if path == packet.packet_receipt.name and kwargs.get("dir_fd") is not None:
                    observed = True
                    self.assertTrue(flags & os.O_NONBLOCK)
                return original_open(path, flags, *args, **kwargs)

            with patch.object(verifier.os, "open", side_effect=audited_open):
                self.assert_rejected(packet, "must be a regular file")
            self.assertTrue(observed)

    def test_authorized_root_intermediate_transplant_is_rejected(self) -> None:
        with _published_packet() as packet:
            parent = packet.root.parent
            retained_parent = parent.with_name(parent.name + "-retained")
            original_validate_packet = verifier._validate_packet_receipt
            transplanted = False

            def transplant_parent(*args: object, **kwargs: object) -> tuple[dict[str, str], Path, str]:
                nonlocal transplanted
                result = original_validate_packet(*args, **kwargs)
                if not transplanted:
                    os.rename(parent, retained_parent)
                    shutil.copytree(retained_parent, parent)
                    transplanted = True
                return result

            try:
                with patch.object(
                    verifier,
                    "_validate_packet_receipt",
                    side_effect=transplant_parent,
                ), self.assertRaisesRegex(
                    verifier.PressureReviewPacketVerificationError,
                    "authorized PIC root component .* path identity changed",
                ):
                    packet.verify()
            finally:
                if retained_parent.exists():
                    if parent.exists():
                        _make_tree_writable(parent)
                        shutil.rmtree(parent)
                    os.rename(retained_parent, parent)
            self.assertTrue(transplanted)

    def test_rejects_packet_root_path_replacement_between_stable_reads(self) -> None:
        with _published_packet() as packet:
            original_scan = verifier._scan_packet_tree
            calls = 0

            def replace_after_first_scan(packet_root_fd: int) -> dict[str, object]:
                nonlocal calls
                result = original_scan(packet_root_fd)
                calls += 1
                if calls == 1:
                    retained = packet.publication / "review-packet-retained"
                    os.rename(packet.packet_root, retained)
                    shutil.copytree(retained, packet.packet_root)
                return result

            with patch.object(
                verifier,
                "_scan_packet_tree",
                side_effect=replace_after_first_scan,
            ):
                self.assert_rejected(
                    packet,
                    r"packet root (?:path identity or metadata )?changed during verification",
                )

    def test_temporary_publication_root_transplant_cannot_redirect_packet_root(self) -> None:
        with _published_packet() as packet:
            retained_publication = packet.root / "publication-retained"
            original_scan = verifier._scan_packet_tree
            original_validate_packet = verifier._validate_packet_receipt
            validation_calls = 0
            scan_calls = 0

            def transplant_publication(
                *args: object,
                **kwargs: object,
            ) -> tuple[dict[str, str], Path, str]:
                nonlocal validation_calls
                result = original_validate_packet(*args, **kwargs)
                validation_calls += 1
                os.rename(packet.publication, retained_publication)
                shutil.copytree(retained_publication, packet.publication)
                alternate_member = packet.packet_root / "PRESSURE_REVIEW_PACKET.md"
                alternate_member.chmod(0o644)
                alternate_member.write_bytes(b"transplanted forged packet\n")
                alternate_member.chmod(0o444)
                return result

            def restore_before_scan(packet_root_fd: int) -> dict[str, object]:
                nonlocal scan_calls
                scan_calls += 1
                if scan_calls == 1:
                    _make_tree_writable(packet.publication)
                    shutil.rmtree(packet.publication)
                    os.rename(retained_publication, packet.publication)
                return original_scan(packet_root_fd)

            with patch.object(
                verifier,
                "_validate_packet_receipt",
                side_effect=transplant_publication,
            ), patch.object(
                verifier,
                "_scan_packet_tree",
                side_effect=restore_before_scan,
            ), self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                "publication root changed during verification",
            ):
                verifier.consume_published_pressure_pilot_review_packet(
                    packet.packet_receipt,
                    aggregate_receipt_binding=packet.aggregate_binding,
                    authorized_pic_root=packet.root,
                )
            self.assertEqual(validation_calls, 1)
            self.assertEqual(scan_calls, 2)

    def test_rejects_packet_root_replacement_after_scans(self) -> None:
        with _published_packet() as packet:
            original_guard_check = verifier._require_guard_absent
            guard_checks = 0

            def replace_before_final_markers(
                publication_fd: int,
                receipt_name: str,
                label: str,
            ) -> None:
                nonlocal guard_checks
                guard_checks += 1
                if guard_checks == 3:
                    retained = packet.publication / "review-packet-retained"
                    os.rename(packet.packet_root, retained)
                    shutil.copytree(retained, packet.packet_root)
                original_guard_check(publication_fd, receipt_name, label)

            with patch.object(
                verifier,
                "_require_guard_absent",
                side_effect=replace_before_final_markers,
            ), self.assertRaisesRegex(
                verifier.PressureReviewPacketVerificationError,
                r"packet root (?:path identity or metadata )?changed during verification",
            ):
                verifier.consume_published_pressure_pilot_review_packet(
                    packet.packet_receipt,
                    aggregate_receipt_binding=packet.aggregate_binding,
                    authorized_pic_root=packet.root,
                )


if __name__ == "__main__":
    unittest.main()
