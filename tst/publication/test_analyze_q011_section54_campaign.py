#!/usr/bin/env python3
"""Focused adversarial tests for the Q-011 Section 5.4 campaign admission gate."""

from __future__ import annotations

from contextlib import contextmanager
import copy
import hashlib
import json
from pathlib import Path
import stat
import struct
import tempfile
from typing import Any, Callable, Iterator
import unittest
from unittest import mock

from tst.publication import analyze_q011_section54_campaign as campaign
from tst.publication import immutable_orion_tree


_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_TIMES = tuple(float(value) for value in range(0, 1300, 100))
_RECEIPT = {
    "schema_version": 1,
    "artifact_role": campaign.ARTIFACT_ROLE,
    "qualification_effect": "retained_qualifying_campaign_attempt",
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_FREEZE_ID = "03a7bd9a-7d4c-4e37-a12b-46de3817eff2"


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fnv1a64(payload: bytes) -> int:
    digest = 14695981039346656037
    for value in payload:
        digest ^= value
        digest = (digest * 1099511628211) & ((1 << 64) - 1)
    return digest


def _make_writable_tree(tree: Path) -> None:
    for path in [tree, *tree.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


def _elf_executable() -> bytes:
    return b"\x7fELF" + bytes((2, 1, 1)) + b"\0" * 57


def _mesh_bin(kind: str, time: float, cycle: int) -> bytes:
    field = "dens" if kind == "rho" else kind
    parameter_header = (
        "<mesh>\n"
        "nx1=1\nnx2=1\nnx3=1\nnghost=0\n"
        "x1min=0.0\nx1max=1.0\n"
        "x2min=0.0\nx2max=1.0\n"
        "x3min=0.0\nx3max=1.0\n"
        "<meshblock>\n"
        "nx1=1\nnx2=1\nnx3=1\n"
    ).encode("ascii")
    header = (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={time:.1f}\n".encode("ascii")
        + f"  cycle={cycle}\n".encode("ascii")
        + b"  size of location=8\n"
        + b"  size of variable=4\n"
        + b"  number of variables=1\n"
        + f"  variables:  {field}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
    )
    block = (
        struct.pack("<10i", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        + struct.pack("<6d", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        + struct.pack("<f", 1.0)
    )
    return header + block


def _particle_vtk(
    time: float,
    cycle: int,
    *,
    ptags: tuple[int, int] = (10, 11),
    nranks: int = 4,
) -> bytes:
    points = ((100.0, 1.0, 0.0), (200.0, 2.0, 0.0))
    integer_scalars = {
        "gid": (0, 1),
        "ptag": ptags,
        "species": (0, 0),
        "cr_source": (1, 1),
    }
    real_scalars = {
        "macro_weight": (1.0, 2.0),
        "birth_time": (45.0, 46.0),
        "deltaf_f0": (1.0, 1.0),
        "deltaf_weight": (0.0, 0.0),
    }
    velocities = ((3.0, 4.0, 0.0), (6.0, 8.0, 0.0))
    payload = bytearray(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= {nranks}  "
            f"cycle={cycle}  variables=prtcl_all\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n"
            f"POINTS {len(points)} float\n"
        ).encode("ascii")
    )
    payload.extend(struct.pack(">6f", *(value for point in points for value in point)))
    payload.extend(f"\n\nPOINT_DATA {len(points)}\n".encode("ascii"))
    for name, values in integer_scalars.items():
        payload.extend(f"\nSCALARS {name} int\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2i", *values))
    for name, values in real_scalars.items():
        payload.extend(f"\nSCALARS {name} float\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2f", *values))
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(
        struct.pack(">6f", *(value for velocity in velocities for value in velocity))
    )
    return bytes(payload)


def _restart_payload(cycle: int) -> bytes:
    return (
        b"<job>\nbasename=q011\n<time>\ncycle="
        + str(cycle).encode("ascii")
        + b"\n<par_end>\n"
        + struct.pack("<4I", cycle, 1, 2, 3)
    )


def _restart_marker(payload: bytes) -> bytes:
    return (
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={len(payload)}\n"
        f"fnv1a64={_fnv1a64(payload):016x}\n"
    ).encode("ascii")


def _stdout_telemetry(*, omit: frozenset[str] = frozenset()) -> bytes:
    lines = ["AthenaK retained stdout fixture"]
    for name in sorted(campaign._Q017_REQUIRED_NAMES - omit):
        value = 2.0 if name == "schema_version" else 1.0
        lines.append(f"q017.telemetry.{name}={value}")
    return ("\n".join(lines) + "\n").encode("ascii")


def _write_file(root: Path, relative: str, payload: bytes) -> dict[str, str]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {"path": relative, "sha256": _sha256(payload)}


def _product(
    root: Path,
    kind: str,
    relative: str,
    payload: bytes,
    snapshot_time: float | None,
) -> dict[str, Any]:
    binding = _write_file(root, relative, payload)
    return {
        "kind": kind,
        "path": binding["path"],
        "sha256": binding["sha256"],
        "snapshot_time": snapshot_time,
    }


def _frozen_candidate_path(suffix: str) -> str:
    return (
        campaign.ORION_BULK_ROOT / "clean_candidates" / _FREEZE_ID / suffix
    ).as_posix()


def _clean_candidate_manifest(
    executable: dict[str, str],
    deck: dict[str, str],
    analyzer: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": 4,
        "freeze_id": _FREEZE_ID,
        "created_utc": "2026-06-01T12:00:00Z",
        "prepared_artifacts": {
            "inventory_path": campaign.PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH,
            "inventory_sha256": "c" * 64,
            "paper_decks": [
                {"path": campaign.ACTIVE_DECK_SOURCE_PATH, "sha256": deck["sha256"]},
                {"path": "inputs/tests/pic_fixture.athinput", "sha256": "d" * 64},
            ],
            "analyzers": [
                {
                    "path": campaign.CAMPAIGN_ANALYZER_SOURCE_PATH,
                    "sha256": analyzer["sha256"],
                },
                {
                    "path": "tst/publication/analyze_q011_section54_outputs.py",
                    "sha256": "e" * 64,
                },
            ],
        },
        "source": {
            "archive_path": _frozen_candidate_path("source.tar"),
            "archive_sha256": "3" * 64,
            "commit_path": _frozen_candidate_path("source.commit"),
            "commit_sha256": "4" * 64,
            "source_bundle_sha256": "2" * 64,
            "git_commit": "1" * 40,
            "git_tree": "5" * 40,
            "worktree_status": "clean",
            "submodule_status": "absent",
            "submodules": [],
        },
        "build": {
            "profile_id": "hip-mpi-release-paper-pic",
            "profile_path": _frozen_candidate_path("build_profile.json"),
            "profile_sha256": "6" * 64,
            "profile_receipt_path": _frozen_candidate_path("profile_receipt.json"),
            "profile_receipt_sha256": "7" * 64,
            "source_archive_sha256": "3" * 64,
            "source_commit_sha256": "4" * 64,
            "source_bundle_sha256": "2" * 64,
            "toolchain": "Frontier fixture toolchain",
            "build_invocations_sha256": "8" * 64,
            "executable_path": _frozen_candidate_path("athena"),
            "executable_sha256": executable["sha256"],
        },
    }


def _append_restart_publication(
    root: Path, products: list[dict[str, Any]], time: float, cycle: int
) -> None:
    restart = _product(
        root, "restart", f"rst/q011.{cycle:05d}.rst", _restart_payload(cycle), time
    )
    products.append(restart)
    products.append(
        _product(
            root,
            "restart_complete",
            restart["path"] + ".complete",
            _restart_marker((root / restart["path"]).read_bytes()),
            time,
        )
    )
    manifest_payload = (
        json.dumps(
            {
                "schema": "ATHENAK_RESTART_MANIFEST_V1",
                "members": [
                    {
                        "path": restart["path"],
                        "size": len((root / restart["path"]).read_bytes()),
                        "fnv1a64": f"{_fnv1a64((root / restart['path']).read_bytes()):016x}",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    manifest = _product(
        root, "restart_manifest", restart["path"] + ".manifest", manifest_payload, time
    )
    products.append(manifest)
    products.append(
        _product(
            root,
            "restart_manifest_complete",
            manifest["path"] + ".complete",
            _restart_marker(manifest_payload),
            time,
        )
    )


def _manifest(root: Path) -> dict[str, Any]:
    executable = _write_file(root, "bindings/athena", _elf_executable())
    (root / executable["path"]).chmod(0o755)
    deck = _write_file(root, "bindings/section54.athinput", b"<job>\nbasename=q011\n")
    analyzer = _write_file(
        root,
        "bindings/analyze_q011_section54_campaign.py",
        campaign.ANALYZER_PATH.read_bytes(),
    )
    candidate_payload = (
        json.dumps(
            _clean_candidate_manifest(executable, deck, analyzer),
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    candidate_manifest = _write_file(
        root, "bindings/clean_candidate_manifest.json", candidate_payload
    )
    preregistration = _write_file(
        root,
        "bindings/q011_section54_preregistration.json",
        campaign.PREREGISTRATION_PATH.read_bytes(),
    )
    products = []
    for cycle, time in enumerate(_TIMES):
        suffix = f"{int(time):05d}"
        for kind in ("rho", "bmag", "prtcl_jx", "j2"):
            products.append(
                _product(
                    root,
                    kind,
                    f"bin/q011.{kind}.{suffix}.bin",
                    _mesh_bin(kind, time, cycle),
                    time,
                )
            )
        products.append(
            _product(
                root,
                "prtcl_all",
                f"pvtk/q011.prtcl_all.{suffix}.part.vtk",
                _particle_vtk(time, cycle),
                time,
            )
        )
        _append_restart_publication(root, products, time, cycle)
    products.append(_product(root, "stdout", "stdout.txt", _stdout_telemetry(), None))
    return {
        "schema_version": 1,
        "record_type": "q011_section54_campaign_run_manifest",
        "campaign_id": campaign.CAMPAIGN_ID,
        "qualification_scope": campaign.QUALIFICATION_SCOPE,
        "authorized_orion_campaign_root": str(root),
        "run_identity": {
            "variant": "three_level_amr_root_dx12_finest_dx3",
            "seed": 23050101,
            "physical_mode": "paper_mhd_pic_vl2_tsc",
            "attempt_id": "q011-amr-seed-23050101-attempt-001",
        },
        "candidate_binding": {
            "git_commit": "1" * 40,
            "source_bundle_sha256": "2" * 64,
            "clean_candidate_manifest": candidate_manifest,
        },
        "artifact_bindings": {
            "executable": executable,
            "deck": deck,
            "analyzer": analyzer,
            "preregistration": preregistration,
        },
        "products": products,
    }


def _find_product(manifest: dict[str, Any], kind: str, time: float | None) -> dict[str, Any]:
    return next(
        product
        for product in manifest["products"]
        if product["kind"] == kind and product["snapshot_time"] == time
    )


def _rewrite_product(root: Path, product: dict[str, Any], payload: bytes) -> None:
    (root / product["path"]).write_bytes(payload)
    product["sha256"] = _sha256(payload)


def _rewrite_clean_candidate(
    root: Path,
    manifest: dict[str, Any],
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    binding = manifest["candidate_binding"]["clean_candidate_manifest"]
    path = root / binding["path"]
    candidate = json.loads(path.read_text(encoding="utf-8"))
    mutate(candidate)
    payload = (json.dumps(candidate, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.write_bytes(payload)
    binding["sha256"] = _sha256(payload)


def _rewrite_bound_artifact(
    root: Path, manifest: dict[str, Any], name: str, payload: bytes
) -> None:
    binding = manifest["artifact_bindings"][name]
    (root / binding["path"]).write_bytes(payload)
    binding["sha256"] = _sha256(payload)


def _rewrite_restart_payload(
    root: Path, manifest: dict[str, Any], time: float, payload: bytes
) -> None:
    restart = _find_product(manifest, "restart", time)
    restart_marker = _find_product(manifest, "restart_complete", time)
    restart_manifest = _find_product(manifest, "restart_manifest", time)
    manifest_marker = _find_product(manifest, "restart_manifest_complete", time)
    _rewrite_product(root, restart, payload)
    _rewrite_product(root, restart_marker, _restart_marker(payload))
    publication = json.loads((root / restart_manifest["path"]).read_text(encoding="utf-8"))
    publication["members"][0]["size"] = len(payload)
    publication["members"][0]["fnv1a64"] = f"{_fnv1a64(payload):016x}"
    publication_payload = (
        json.dumps(publication, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    _rewrite_product(root, restart_manifest, publication_payload)
    _rewrite_product(root, manifest_marker, _restart_marker(publication_payload))


@contextmanager
def _frozen_fixture(
    mutate: Callable[[Path, dict[str, Any]], None] | None = None,
    *,
    receipt: dict[str, Any] | None = None,
    external_side_effect: Exception | None = None,
) -> Iterator[tuple[Path, Path, str]]:
    with tempfile.TemporaryDirectory() as directory:
        authorized_root = Path(directory)
        tree = authorized_root / "campaign"
        tree.mkdir()
        manifest = _manifest(tree)
        if mutate is not None:
            mutate(tree, manifest)
        (tree / campaign.MANIFEST_NAME).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        patch_kwargs = (
            {"side_effect": external_side_effect}
            if external_side_effect is not None
            else {
                "return_value": {
                    "freeze_id": _FREEZE_ID,
                    "referenced_manifest_sha256": manifest["candidate_binding"][
                        "clean_candidate_manifest"
                    ]["sha256"],
                    "retained_executable_sha256": manifest["artifact_bindings"][
                        "executable"
                    ]["sha256"],
                    "validated_submodule_count": 0,
                    "validation": "fixture_control_plane_closure",
                }
            }
        )
        try:
            frozen = immutable_orion_tree.freeze_tree(
                tree,
                _RECEIPT if receipt is None else receipt,
                authorized_root=authorized_root,
            )
            with mock.patch.object(
                campaign, "_validate_external_clean_candidate_closure", **patch_kwargs
            ):
                yield authorized_root, tree, frozen["inventory_sha256"]
        finally:
            _make_writable_tree(tree)


def _qualify(
    authorized_root: Path,
    tree: Path,
    inventory_sha256: str,
) -> dict[str, Any]:
    return campaign.qualify_campaign(
        tree,
        inventory_sha256,
        authorized_orion_root=authorized_root,
    )


class Q011Section54CampaignAdmissionTests(unittest.TestCase):
    def assert_rejected(self, result: dict[str, Any], code: str) -> None:
        self.assertFalse(result["admitted_for_follow_on_numerical_qualification"])
        self.assertIsNone(result["admission"])
        self.assertEqual(result["failure_reasons"][0]["code"], code)

    def test_valid_fixture_is_admitted_without_claiming_numerical_closure(self) -> None:
        with _frozen_fixture() as fixture:
            result = _qualify(*fixture)
        self.assertTrue(result["admitted_for_follow_on_numerical_qualification"])
        self.assertFalse(result["final_claim_closure"])
        self.assertEqual(result["failure_reasons"], [])
        admission = result["admission"]
        self.assertEqual(set(admission["endpoint_products"]), {"500.0", "1200.0"})
        self.assertEqual(len(admission["retained_snapshot_products"]), 13)
        self.assertEqual(len(admission["restart_publications"]), 13)
        self.assertEqual(
            admission["particle_endpoints"]["500.0"]["execution_header"]["nranks"], 4
        )
        self.assertEqual(
            admission["preregistration_binding"]["binding_scope"],
            "complete_retained_bytes_equal_invoked_frozen_policy",
        )
        self.assertEqual(
            admission["numerical_qualification_status"],
            "not_evaluated_by_artifact_admission_slice",
        )

    def test_path_escape_is_rejected_without_partial_success(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["artifact_bindings"]["deck"]["path"] = "../outside.athinput"

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "unsafe_relative_path")

    def test_fictional_text_executable_is_rejected_even_when_self_consistent(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_bound_artifact(root, manifest, "executable", b"fictional executable\n")
            digest = manifest["artifact_bindings"]["executable"]["sha256"]
            _rewrite_clean_candidate(
                root,
                manifest,
                lambda candidate: candidate["build"].__setitem__(
                    "executable_sha256", digest
                ),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_executable_elf")

    def test_external_clean_candidate_closure_failure_is_rejected(self) -> None:
        error = campaign.QualificationError(
            "fixture external candidate closure failed",
            code="external_clean_candidate_closure_error",
        )
        with _frozen_fixture(external_side_effect=error) as fixture:
            self.assert_rejected(
                _qualify(*fixture), "external_clean_candidate_closure_error"
            )

    def test_fictional_plain_mesh_bin_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(root, _find_product(manifest, "rho", 500.0), b"rho\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_mesh_bin")

    def test_valid_mesh_bin_with_wrong_observable_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "rho", 500.0)
            _rewrite_product(root, product, _mesh_bin("bmag", 500.0, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_mesh_bin")

    def test_duplicate_ptag_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            _rewrite_product(root, product, _particle_vtk(500.0, 5, ptags=(10, 10)))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_prtcl_all")

    def test_pvtk_embedded_time_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            _rewrite_product(root, product, _particle_vtk(999.0, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "snapshot_metadata_drift")

    def test_pvtk_manifest_time_must_directly_match_embedded_time(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "prtcl_all", 500.0)
            product["snapshot_time"] = 500.00000075
            _rewrite_product(root, product, _particle_vtk(499.99999925, 5))

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "snapshot_metadata_drift")

    def test_wrong_freeze_receipt_semantics_are_rejected(self) -> None:
        receipt = dict(_RECEIPT)
        receipt["qualification_effect"] = "fictional_effect"
        with _frozen_fixture(receipt=receipt) as fixture:
            self.assert_rejected(_qualify(*fixture), "freeze_receipt_semantics_drift")

    def test_cadence_omission_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            product = _find_product(manifest, "j2", 100.0)
            manifest["products"].remove(product)
            (root / product["path"]).unlink()

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "missing_snapshot_cadence")

    def test_duplicate_endpoint_is_rejected_as_ambiguous_cadence(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            manifest["products"].append(
                _product(
                    root,
                    "rho",
                    "bin/q011.rho.00500.duplicate.bin",
                    _mesh_bin("rho", 500.0, 5),
                    500.0,
                )
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "ambiguous_snapshot_cadence")

    def test_fictional_plain_restart_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_restart_payload(root, manifest, 500.0, b"restart\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_restart_payload")

    def test_restart_marker_corruption_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            marker = _find_product(manifest, "restart_complete", 500.0)
            _rewrite_product(root, marker, b"restart complete\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "restart_marker_corruption")

    def test_fictional_plain_stdout_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(root, _find_product(manifest, "stdout", None), b"completed\n")

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_q017_telemetry")

    def test_stdout_telemetry_omission_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite_product(
                root,
                _find_product(manifest, "stdout", None),
                _stdout_telemetry(omit=frozenset({"particles.total"})),
            )

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "invalid_q017_telemetry")

    def test_complete_preregistration_checksum_drift_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            payload = campaign.PREREGISTRATION_PATH.read_bytes() + b" "
            _rewrite_bound_artifact(root, manifest, "preregistration", payload)

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "hash_drift")

    def test_forged_outer_candidate_identity_is_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["candidate_binding"]["git_commit"] = "9" * 40

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "clean_candidate_binding_drift")

    def test_forged_self_consistent_deck_or_analyzer_is_rejected(self) -> None:
        for name in ("deck", "analyzer"):
            with self.subTest(name=name):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, name: str = name
                ) -> None:
                    _rewrite_bound_artifact(
                        root, manifest, name, f"forged {name}\n".encode("ascii")
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(
                        _qualify(*fixture), "clean_candidate_binding_drift"
                    )

    def test_clean_candidate_schema_aliases_are_rejected(self) -> None:
        for value in (True, 4.0, "4"):
            with self.subTest(value=value):
                def mutate(
                    root: Path, manifest: dict[str, Any], *, value: object = value
                ) -> None:
                    _rewrite_clean_candidate(
                        root,
                        manifest,
                        lambda candidate: candidate.__setitem__("schema_version", value),
                    )

                with _frozen_fixture(mutate) as fixture:
                    self.assert_rejected(_qualify(*fixture), "schema_type_error")

    def test_duplicate_prepared_path_is_rejected(self) -> None:
        def mutate(root: Path, manifest: dict[str, Any]) -> None:
            def duplicate(candidate: dict[str, Any]) -> None:
                decks = candidate["prepared_artifacts"]["paper_decks"]
                decks.append(copy.deepcopy(decks[0]))

            _rewrite_clean_candidate(root, manifest, duplicate)

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "clean_candidate_schema_error")

    def test_boolean_seed_alias_is_rejected(self) -> None:
        def mutate(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["run_identity"]["seed"] = True

        with _frozen_fixture(mutate) as fixture:
            self.assert_rejected(_qualify(*fixture), "schema_type_error")

    def test_weighted_spectrum_utility_uses_preregistered_fixed_bins(self) -> None:
        policy = json.loads(campaign.PREREGISTRATION_PATH.read_text(encoding="utf-8"))
        report = campaign.weighted_spectrum_from_chi(
            [0.5, 1.0, 2.0, 2048.0],
            [3.0, 4.0, 5.0, 6.0],
            policy,
        )
        self.assertEqual(report["underflow_count"], 1)
        self.assertEqual(report["overflow_count"], 1)
        self.assertEqual(report["underflow_weight"], 3.0)
        self.assertEqual(report["overflow_weight"], 6.0)
        self.assertEqual(report["macro_weight_in_bins"], 9.0)
        self.assertEqual(report["admitted_weight"], 18.0)
        self.assertEqual(report["overflow_macro_weight_fraction"], 1.0 / 3.0)
        self.assertFalse(report["overflow_gate"]["passed"])


if __name__ == "__main__":
    unittest.main()
