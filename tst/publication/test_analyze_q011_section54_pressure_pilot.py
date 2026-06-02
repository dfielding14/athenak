#!/usr/bin/env python3
"""Adversarial tests for the nonqualifying Q-011 pressure-pilot analyzer."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import struct
import tempfile
from typing import Any, Callable, Iterator
import unittest

import numpy as np

from tst.publication import analyze_q011_section54_pressure_pilot as pilot


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fnv1a64(payload: bytes) -> str:
    value = 14695981039346656037
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def _put(root: Path, relative: str, payload: bytes) -> dict[str, str]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {"path": relative, "sha256": _sha256(payload)}


def _rewrite(root: Path, binding: dict[str, str], payload: bytes) -> None:
    (root / binding["path"]).write_bytes(payload)
    binding["sha256"] = _sha256(payload)


def _parameters(ps_p0: str) -> bytes:
    return (
        "<mesh>\n"
        "nx1=100\n"
        "nx2=20\n"
        "nx3=1\n"
        "nghost=2\n"
        "x1min=0.0\n"
        "x1max=1200\n"
        "x2min=0.0\n"
        "x2max=240\n"
        "x3min=0.0\n"
        "x3max=1.0\n"
        "<meshblock>\n"
        "nx1=100\n"
        "nx2=20\n"
        "nx3=1\n"
        "<mesh_refinement>\n"
        "refinement=none\n"
        "num_levels=1\n"
        "<time>\n"
        "tlim=60\n"
        "nlim=4096\n"
        "ndiag=50\n"
        "<mhd>\n"
        "gamma=1.66666666667\n"
        "<particles>\n"
        "pic_physical_mode=paper_mhd_pic_vl2_tsc\n"
        "<problem>\n"
        f"ps_p0={ps_p0}\n"
        "ps_enable_curvature_amr=false\n"
        "ps_feedback_diag_dcycle=50\n"
        "<output1>\n"
        "variable=mhd_w_bcc\n"
        "id=mhd_w_bcc\n"
        "dt=15\n"
        "<output2>\n"
        "variable=mhd_bmag\n"
        "id=bmag\n"
        "dt=15\n"
        "<output3>\n"
        "variable=prtcl_jx\n"
        "id=prtcl_jx\n"
        "dt=15\n"
        "<output4>\n"
        "variable=mhd_j2\n"
        "id=j2\n"
        "dt=15\n"
        "<output5>\n"
        "variable=prtcl_all\n"
        "id=prtcl_all\n"
        "dt=15\n"
        "<output6>\n"
        "file_type=rst\n"
        "dt=15\n"
    ).encode("ascii")


def _binary(
    time: float,
    ps_p0: str,
    fields: tuple[str, ...],
    *,
    negative_pressure: bool = False,
) -> bytes:
    shape = (1, 20, 100)
    x = np.arange(100, dtype=np.float32)
    base = np.broadcast_to(x, shape)
    values = {
        "dens": 1.0 + 0.001 * base,
        "eint": np.full(shape, -1.0 if negative_pressure else 2.0, dtype=np.float32),
        "velx": np.full(shape, 30.0, dtype=np.float32),
        "vely": np.zeros(shape, dtype=np.float32),
        "velz": np.zeros(shape, dtype=np.float32),
        "bcc1": np.ones(shape, dtype=np.float32),
        "bcc2": np.full(shape, 0.5, dtype=np.float32),
        "bcc3": np.zeros(shape, dtype=np.float32),
        "bmag": np.full(shape, np.sqrt(1.25), dtype=np.float32),
        "prtcl_jx": np.full(shape, 0.25, dtype=np.float32),
        "j2": np.full(shape, 0.125, dtype=np.float32),
    }
    parameters = _parameters(ps_p0)
    header = (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={time}\n".encode("ascii")
        + b"  cycle=50\n"
        b"  size of location=8\n"
        b"  size of variable=4\n"
        + f"  number of variables={len(fields)}\n".encode("ascii")
        + b"  variables:  "
        + b"  ".join(name.encode("ascii") for name in fields)
        + b"  \n"
        + f"  header offset={len(parameters)}\n".encode("ascii")
        + parameters
    )
    block = (
        struct.pack("<10i", 0, 99, 0, 19, 0, 0, 0, 0, 0, 0)
        + struct.pack("<6d", 0.0, 1200.0, 0.0, 240.0, 0.0, 1.0)
        + np.concatenate([values[name].reshape(-1) for name in fields]).astype("<f4").tobytes()
    )
    return header + block


def _particle_vtk(time: float, *, source: tuple[int, int] = (1, 1),
                  birth_time: tuple[float, float] = (45.0, 46.0)) -> bytes:
    points = ((100.0, 1.0, 0.0), (200.0, 2.0, 0.0))
    integer_scalars = {
        "gid": (0, 1),
        "ptag": (10, 11),
        "species": (0, 0),
        "cr_source": source,
    }
    real_scalars = {
        "macro_weight": (1.0, 2.0),
        "birth_time": birth_time,
        "deltaf_f0": (1.0, 1.0),
        "deltaf_weight": (0.0, 0.0),
    }
    payload = bytearray(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {time}  nranks= 1  "
            "cycle=50  variables=prtcl_all\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n"
            "POINTS 2 float\n"
        ).encode("ascii")
    )
    payload.extend(struct.pack(">6f", *(value for point in points for value in point)))
    payload.extend(b"\n\nPOINT_DATA 2\n")
    for name, values in integer_scalars.items():
        payload.extend(f"\nSCALARS {name} int\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2i", *values))
    for name, values in real_scalars.items():
        payload.extend(f"\nSCALARS {name} float\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(struct.pack(">2f", *values))
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(struct.pack(">6f", 3.0, 4.0, 0.0, 6.0, 8.0, 0.0))
    return bytes(payload)


def _stdout() -> bytes:
    telemetry = {
        "schema_version": 2.0,
        "mpi.ranks": 1.0,
        "cycles": 100.0,
        "meshblocks.total": 1.0,
        "mesh.active_cells": 2000.0,
        "particles.total": 2.0,
        "throughput.zone_cycles_per_second": 1000.0,
        "throughput.particle_updates_per_second": 2000.0,
        "load.meshblock_efficiency": 1.0,
        "load.particle_efficiency": 1.0,
        "amr.enabled": 0.0,
        "particle_memory.resident_records.bytes_total": 448.0,
        "particle_memory.invalid_records": 0.0,
        (
            "particle_memory.athenak_owned_tracked_kokkos_views."
            "allocated_snapshot_bytes_total"
        ): 4096.0,
        (
            "particle_memory.athenak_owned_tracked_kokkos_views."
            "allocated_high_water_bytes_rank_max"
        ): 8192.0,
    }
    lines = [f"q017.telemetry.{name}={value:.17e}" for name, value in telemetry.items()]
    lines.extend([
        (
            "pic_parallel_shock feedback_diag_cfg: pusher=2 pic_feedback_mode=1 "
            "couple_moments_to_mhd=1 mom_feedback=1 eng_feedback=1 "
            "pic_enable_2d3v=1 use_delta_feedback=1"
        ),
        (
            "pic_parallel_shock feedback_diag: cycle=50 time=5.000000e+01 npart=2 "
            "ncell=2000 j_rms=(1,2,3) dpdt_rms=(4,5,6) dedt_rms=7"
        ),
        (
            "pic_parallel_shock source_transaction_diag: cycle=50 "
            "applied=(1,2,3,4,5) expected=(1,2,3,4,5)"
        ),
        "pic_parallel_shock removed_cr_sink: count=1 mass=1 momentum=(1,2,3) energy=4",
    ])
    return ("\n".join(lines) + "\n").encode("ascii")


def _marker(payload: bytes) -> bytes:
    return (
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={len(payload)}\n"
        f"fnv1a64={_fnv1a64(payload)}\n"
    ).encode("ascii")


def _restart(root: Path, case_id: str) -> dict[str, Any]:
    artifact_path = f"cases/{case_id}/rst/{case_id}.00004.rst"
    artifact = (f"terminal restart {case_id}\n").encode("ascii")
    artifact_binding = _put(root, artifact_path, artifact)
    artifact_marker = _put(root, artifact_path + ".complete", _marker(artifact))
    manifest_path = artifact_path + ".manifest"
    manifest = {
        "schema": "ATHENAK_RESTART_MANIFEST_V1",
        "members": [{
            "path": f"rst/{case_id}.00004.rst",
            "size": len(artifact),
            "fnv1a64": _fnv1a64(artifact),
        }],
    }
    manifest_payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    return {
        "time": 60.0,
        "manifest": _put(root, manifest_path, manifest_payload),
        "manifest_complete": _put(root, manifest_path + ".complete", _marker(manifest_payload)),
        "members": [{"artifact": artifact_binding, "complete": artifact_marker}],
    }


def _rank_shard_restart(
    root: Path,
    manifest: dict[str, Any],
    *,
    case_index: int = 0,
    ranks: tuple[int, ...] = (0, 1),
) -> None:
    case = manifest["cases"][case_index]
    case_id = case["case_id"]
    restart = case["terminal_restart"]
    for member in restart["members"]:
        (root / member["artifact"]["path"]).unlink()
        (root / member["complete"]["path"]).unlink()
    outer_members = []
    manifest_members = []
    for rank in ranks:
        declared = f"rst/rank_{rank:08d}/{case_id}.00004.rst"
        relative = f"cases/{case_id}/{declared}"
        payload = f"terminal restart {case_id} rank {rank}\n".encode("ascii")
        outer_members.append({
            "artifact": _put(root, relative, payload),
            "complete": _put(root, relative + ".complete", _marker(payload)),
        })
        manifest_members.append({
            "path": declared,
            "size": len(payload),
            "fnv1a64": _fnv1a64(payload),
        })
    manifest_payload = (
        json.dumps(
            {"schema": "ATHENAK_RESTART_MANIFEST_V1", "members": manifest_members},
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    _rewrite(root, restart["manifest"], manifest_payload)
    _rewrite(root, restart["manifest_complete"], _marker(manifest_payload))
    restart["members"] = outer_members


def _bundle(root: Path) -> dict[str, Any]:
    cases = []
    for case_id, ps_p0, argv in pilot._CASES:
        snapshots = []
        for index, time in enumerate(pilot._TIMES):
            paths = pilot._expected_snapshot_paths(case_id, index)
            snapshots.append({
                "time": time,
                "mhd_w_bcc": _put(root, paths["mhd_w_bcc"], _binary(time, argv, pilot._MHD_FIELDS)),
                "bmag": _put(root, paths["bmag"], _binary(time, argv, ("bmag",))),
                "prtcl_jx": _put(root, paths["prtcl_jx"], _binary(time, argv, ("prtcl_jx",))),
                "j2": _put(root, paths["j2"], _binary(time, argv, ("j2",))),
                "prtcl_all": _put(root, paths["prtcl_all"], _particle_vtk(time)),
            })
        cases.append({
            "case_id": case_id,
            "ps_p0": ps_p0,
            "overrides": [*pilot._FIXED_OVERRIDES, f"problem/ps_p0={argv}"],
            "snapshots": snapshots,
            "stdout": _put(root, f"cases/{case_id}/stdout.txt", _stdout()),
            "terminal_restart": _restart(root, case_id),
        })
    return {
        "schema_version": 1,
        "record_type": "q011_section54_pressure_pilot_bundle_manifest",
        "evidence_class": pilot.EVIDENCE_CLASS,
        "qualification_effect": pilot.QUALIFICATION_EFFECT,
        "active_deck_binding": {
            "path": pilot.ACTIVE_DECK_PATH.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(pilot.ACTIVE_DECK_PATH.read_bytes()),
        },
        "preregistration_binding": {
            "path": pilot.PREREGISTRATION_PATH.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(pilot.PREREGISTRATION_PATH.read_bytes()),
        },
        "cases": cases,
    }


def _write_manifest(root: Path, manifest: dict[str, Any]) -> str:
    payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    (root / pilot.MANIFEST_NAME).write_bytes(payload)
    return _sha256(payload)


@contextmanager
def _fixture(
    mutate: Callable[[Path, dict[str, Any]], None] | None = None,
) -> Iterator[tuple[Path, dict[str, Any], str]]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory).resolve() / "bundle"
        root.mkdir()
        manifest = _bundle(root)
        if mutate is not None:
            mutate(root, manifest)
        digest = _write_manifest(root, manifest)
        yield root, manifest, digest


def _analyze(fixture: tuple[Path, dict[str, Any], str]) -> dict[str, Any]:
    root, _manifest, digest = fixture
    return pilot.analyze_pressure_pilot_bundle(root, digest)


class Q011Section54PressurePilotTests(unittest.TestCase):
    def test_valid_bundle_emits_engineering_only_overlay_profiles(self) -> None:
        with _fixture() as fixture:
            result = _analyze(fixture)
        self.assertEqual(result["status"], "pass_engineering_calibration_only")
        self.assertEqual(result["evidence_class"], "engineering_calibration_only")
        self.assertFalse(result["scientific_evidence_eligible"])
        self.assertFalse(result["sun_bai_claim"])
        self.assertFalse(result["frontier_execution_authorized"])
        self.assertEqual(result["overlay_profile_fields"], ["rho", "p", "vx", "|B|"])
        self.assertEqual(len(result["overlay_profiles"]), 20)
        self.assertEqual(len(result["case_summaries"]), 4)
        profile = result["overlay_profiles"][0]
        self.assertEqual(len(profile["x1_c_over_omega_pi"]), 100)
        self.assertTrue(all(value >= 0.0 for value in profile["p_y_average"]))
        self.assertEqual(
            result["case_summaries"][0]["terminal_restart"],
            {"layout": "shared", "member_count": 1},
        )

    def test_exact_four_case_identity_and_override_contract_fail_closed(self) -> None:
        mutations = {
            "missing case": lambda _root, manifest: manifest["cases"].pop(),
            "pressure drift": lambda _root, manifest: manifest["cases"][0].__setitem__("ps_p0", 0.9),
            "override drift": lambda _root, manifest: manifest["cases"][0]["overrides"].__setitem__(0, "mesh/nx1=99"),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_raw_checksum_and_inventory_drift_fail_closed(self) -> None:
        def checksum(_root: Path, manifest: dict[str, Any]) -> None:
            manifest["cases"][0]["snapshots"][0]["bmag"]["sha256"] = "0" * 64

        def inventory(root: Path, _manifest: dict[str, Any]) -> None:
            _put(root, "cases/ps_p0_1p00/undeclared.txt", b"undeclared\n")

        for label, mutate in (("checksum", checksum), ("inventory", inventory)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_malformed_binary_and_internal_time_drift_fail_closed(self) -> None:
        def malformed(root: Path, manifest: dict[str, Any]) -> None:
            _rewrite(root, manifest["cases"][0]["snapshots"][0]["mhd_w_bcc"], b"not Athena\n")

        def time_drift(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["snapshots"][0]["mhd_w_bcc"]
            _rewrite(root, binding, _binary(1.0, "1.0", pilot._MHD_FIELDS))

        for label, mutate in (("malformed", malformed), ("time", time_drift)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_mhd_variable_and_pressure_drift_fail_closed(self) -> None:
        def variables(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["snapshots"][0]["mhd_w_bcc"]
            _rewrite(root, binding, _binary(0.0, "1.0", pilot._MHD_FIELDS[:-1]))

        def pressure(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["snapshots"][0]["mhd_w_bcc"]
            _rewrite(root, binding, _binary(0.0, "1.0", pilot._MHD_FIELDS, negative_pressure=True))

        for label, mutate in (("variables", variables), ("pressure", pressure)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_particle_provenance_and_terminal_startup_cohort_fail_closed(self) -> None:
        def provenance(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["snapshots"][0]["prtcl_all"]
            _rewrite(root, binding, _particle_vtk(0.0, source=(1, 2)))

        def cohort(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["snapshots"][-1]["prtcl_all"]
            _rewrite(root, binding, _particle_vtk(60.0, birth_time=(44.0, 46.0)))

        for label, mutate in (("provenance", provenance), ("cohort", cohort)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_q017_and_shock_stdout_requirements_fail_closed(self) -> None:
        def missing_telemetry(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["stdout"]
            payload = (root / binding["path"]).read_bytes()
            payload = b"\n".join(
                line for line in payload.splitlines()
                if not line.startswith(b"q017.telemetry.schema_version=")
            ) + b"\n"
            _rewrite(root, binding, payload)

        def missing_diagnostic(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["stdout"]
            payload = (root / binding["path"]).read_bytes()
            payload = b"\n".join(
                line for line in payload.splitlines()
                if not line.startswith(b"pic_parallel_shock source_transaction_diag:")
            ) + b"\n"
            _rewrite(root, binding, payload)

        for label, mutate in (("telemetry", missing_telemetry), ("diagnostic", missing_diagnostic)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_restart_marker_and_manifest_integrity_fail_closed(self) -> None:
        def marker(root: Path, manifest: dict[str, Any]) -> None:
            binding = manifest["cases"][0]["terminal_restart"]["members"][0]["complete"]
            _rewrite(root, binding, b"ATHENAK_RESTART_COMPLETE_V1\nsize=1\nfnv1a64=0000000000000000\n")

        def restart_manifest(root: Path, manifest: dict[str, Any]) -> None:
            restart = manifest["cases"][0]["terminal_restart"]
            path = root / restart["manifest"]["path"]
            decoded = json.loads(path.read_text(encoding="utf-8"))
            decoded["members"][0]["fnv1a64"] = "0" * 16
            payload = (json.dumps(decoded, indent=2, sort_keys=True) + "\n").encode("utf-8")
            _rewrite(root, restart["manifest"], payload)
            _rewrite(root, restart["manifest_complete"], _marker(payload))

        for label, mutate in (("marker", marker), ("manifest", restart_manifest)):
            with self.subTest(label=label):
                with _fixture(mutate) as fixture:
                    with self.assertRaises(pilot.PilotAnalysisError):
                        _analyze(fixture)

    def test_rank_sharded_restart_is_accepted_and_rank_gaps_fail_closed(self) -> None:
        with _fixture(lambda root, manifest: _rank_shard_restart(root, manifest)) as fixture:
            result = _analyze(fixture)
        self.assertEqual(
            result["case_summaries"][0]["terminal_restart"],
            {"layout": "rank_sharded", "member_count": 2},
        )

        with _fixture(
            lambda root, manifest: _rank_shard_restart(root, manifest, ranks=(0, 2))
        ) as fixture:
            with self.assertRaisesRegex(pilot.PilotAnalysisError, "not contiguous"):
                _analyze(fixture)

    def test_manifest_sha256_is_an_external_fail_closed_binding(self) -> None:
        with _fixture() as (root, _manifest, _digest):
            with self.assertRaisesRegex(pilot.PilotAnalysisError, "manifest SHA-256 drifted"):
                pilot.analyze_pressure_pilot_bundle(root, "0" * 64)

    def test_preregistration_is_explicitly_nonqualifying_and_non_authorizing(self) -> None:
        policy = json.loads(pilot.PREREGISTRATION_PATH.read_text(encoding="utf-8"))
        self.assertEqual(policy["classification"], "engineering_calibration_only")
        self.assertEqual(
            policy["qualification_effect"],
            "none_no_sun_bai_claim_no_execution_authorization",
        )
        self.assertFalse(policy["execution_policy"]["frontier_execution_authorized_by_this_record"])
        self.assertFalse(policy["execution_policy"]["scheduler_commands_authorized_by_this_record"])
        self.assertFalse(policy["execution_policy"]["scientific_evidence_eligible"])
        self.assertFalse(policy["execution_policy"]["sun_bai_claim"])
        self.assertEqual(
            [item["problem_ps_p0"] for item in policy["pilot_contract"]["cases"]],
            [1.0, 0.05, 0.1, 0.2],
        )


if __name__ == "__main__":
    unittest.main()
