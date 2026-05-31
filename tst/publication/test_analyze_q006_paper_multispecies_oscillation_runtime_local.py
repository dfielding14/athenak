#!/usr/bin/env python3
"""Focused tests for the bounded serial-host Q-006 runtime-local successor."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import stat
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

from tst.publication import analyze_q006_paper_multispecies_oscillation_runtime_local as q006
from tst.publication.pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q006_paper_multispecies_oscillation_runtime_local_successor_2026-05-31.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _initial_arrays() -> tuple[np.ndarray, ...]:
    density = np.ones((16, 4, 4, 4))
    mhd_velocity = np.zeros(density.shape + (3,))
    mhd_velocity[..., 1] = -0.3
    cell_volume = np.ones_like(density)
    count = 16 * 4 * 4 * 4 * 128
    particle_velocity = np.zeros((count, 3))
    particle_velocity[:, 1] = 0.1
    macro_weight = np.ones(count)
    species = np.arange(count, dtype=np.int64) % 2
    tags = np.arange(count, dtype=np.int64)
    return density, mhd_velocity, cell_volume, particle_velocity, macro_weight, species, tags


class Q006PaperMultispeciesOscillationRuntimeLocalTests(unittest.TestCase):
    def test_three_additive_runtime_local_decks_are_exact_isothermal_momentum_only(self) -> None:
        decks = q006.validate_decks()
        self.assertEqual(
            [deck["grid_setup"] for deck in decks],
            ["uniform", "smr", "audited_amr_runtime_local"],
        )
        self.assertTrue(all(not deck["qualifying_evidence"] for deck in decks))
        for path in q006.DECKS.values():
            text = path.read_text(encoding="utf-8")
            self.assertIn("eos             = isothermal", text)
            self.assertIn("couple_moments_energy_to_mhd       = false", text)
            self.assertNotIn("gamma       =", text)

    def test_additive_registration_preserves_deltaf_allowance_and_historical_hashes(self) -> None:
        registration = q006.validate_registration()
        self.assertTrue(registration["fresh_dispatch"])
        self.assertTrue(registration["restart_dispatch"])
        self.assertTrue(registration["deltaf_allowance_preserved"])
        self.assertTrue(registration["fullf_allowance_added"])
        self.assertTrue(q006.validate_historical_preparation_bindings()["preserved"])

    def test_static_descriptor_retains_only_bounded_nonclaims(self) -> None:
        descriptor = q006.static_descriptor()
        self.assertEqual(descriptor["qualification_effect"], "none")
        self.assertFalse(descriptor["qualifying_evidence"])
        self.assertFalse(descriptor["long_horizon_qualification_claimed"])
        self.assertFalse(descriptor["true_amr_policy_qualification_claimed"])
        self.assertFalse(descriptor["mpi_qualification_claimed"])
        self.assertFalse(descriptor["gpu_qualification_claimed"])
        self.assertFalse(descriptor["frontier_qualification_claimed"])
        self.assertFalse(descriptor["external_review_claimed"])
        contract = descriptor["analytical_contract"]
        self.assertEqual(contract["oscillation_omega"], 2.0)
        self.assertEqual(contract["oscillation_period"], np.pi)

    def test_initial_weighted_mechanics_recompute_exact_com_and_kinetic_energy(self) -> None:
        report = q006.analyze_snapshot_arrays(*_initial_arrays())
        self.assertFalse(report["qualifying_evidence"])
        self.assertEqual(report["particle_count"], 131072)
        self.assertEqual(report["species_mass_density"], [1.5, 1.5])
        np.testing.assert_allclose(report["volume_averaged_gas_velocity"], [0.0, -0.3, 0.0])
        np.testing.assert_allclose(report["total_momentum"], [0.0, 0.0, 0.0], atol=1.0e-10)
        self.assertAlmostEqual(report["total_kinetic_energy_density"], 0.06)

    def test_snapshot_recompute_rejects_nonfinite_nonunique_and_volume_drift(self) -> None:
        arrays = list(_initial_arrays())
        arrays[3] = arrays[3].copy()
        arrays[3][0, 0] = np.nan
        with self.assertRaisesRegex(q006.AuditError, "finite"):
            q006.analyze_snapshot_arrays(*arrays)
        arrays = list(_initial_arrays())
        arrays[6] = arrays[6].copy()
        arrays[6][1] = arrays[6][0]
        with self.assertRaisesRegex(q006.AuditError, "unique"):
            q006.analyze_snapshot_arrays(*arrays)
        arrays = list(_initial_arrays())
        arrays[2] = arrays[2] * 0.5
        with self.assertRaisesRegex(q006.AuditError, "physical volume"):
            q006.analyze_snapshot_arrays(*arrays)

    def test_raw_snapshot_extractor_hashes_inputs(self) -> None:
        arrays = _initial_arrays()
        particle = ParticleVTKData(
            points=np.zeros((arrays[3].shape[0], 3)),
            scalars={
                "macro_weight": arrays[4],
                "species": arrays[5],
                "ptag": arrays[6],
            },
            vectors={"vel": arrays[3]},
        )
        mhd = {
            "mb_data": {
                "dens": arrays[0],
                "velx": arrays[1][..., 0],
                "vely": arrays[1][..., 1],
                "velz": arrays[1][..., 2],
                "bcc1": np.zeros_like(arrays[0]),
                "bcc2": np.zeros_like(arrays[0]),
                "bcc3": np.ones_like(arrays[0]),
            },
            "mb_logical": np.column_stack((
                np.arange(16), np.zeros((16, 3), dtype=np.int64)
            )),
            "time": 0.0,
            "cycle": 0,
            "x1min": 0.0,
            "x1max": 16.0,
            "x2min": 0.0,
            "x2max": 8.0,
            "x3min": 0.0,
            "x3max": 8.0,
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mhd_path = root / "snapshot.bin"
            particle_path = root / "snapshot.vtk"
            mhd_path.write_bytes(b"mhd")
            particle_payload = (
                b"# vtk DataFile Version 2.0\n"
                b"# AthenaK particle data at time= 0  nranks= 1  cycle=0  "
                b"variables=prtcl_all\nBINARY\nparticles"
            )
            particle_path.write_bytes(particle_payload)
            with mock.patch.object(q006, "_read_mhd_binary", return_value=mhd), \
                 mock.patch.object(q006, "read_particle_vtk", return_value=particle):
                report = q006.extract_runtime_snapshot(mhd_path, particle_path)
        self.assertEqual(report["immutable_artifacts"]["mhd_w_bcc"]["sha256"],
                         hashlib.sha256(b"mhd").hexdigest())
        self.assertEqual(report["immutable_artifacts"]["prtcl_all"]["sha256"],
                         hashlib.sha256(particle_payload).hexdigest())
        self.assertEqual(
            report["pvtk_execution_metadata"],
            {"time": 0.0, "nranks": 1, "cycle": 0, "variables": "prtcl_all"},
        )

    def test_direct_cli_static_import_path_works_without_pythonpath(self) -> None:
        proc = subprocess.run(
            [sys.executable, str(q006.__file__), "static"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["campaign_id"], q006.CAMPAIGN_ID)

    def test_inventory_freeze_and_verify_are_orion_scoped(self) -> None:
        original_root = q006.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            (tree / "stdout.txt").write_text("bounded\n", encoding="utf-8")
            (tree / "bin").mkdir()
            (tree / "bin/snapshot.bin").write_bytes(b"snapshot")
            q006.ORION_BULK_ROOT = Path(directory)
            try:
                report = q006.freeze_tree(tree)
                self.assertEqual(report["inventoried_file_count"], 3)
                self.assertTrue(report["recursively_read_only"])
                verified = q006.verify_frozen_tree(tree, report["inventory_sha256"])
                self.assertEqual(verified["inventory_sha256"], report["inventory_sha256"])
                self.assertEqual(verified["writable_entries"], [])
            finally:
                for path in [tree, *tree.rglob("*")]:
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                q006.ORION_BULK_ROOT = original_root

    def test_strict_deck_parser_rejects_duplicate_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.athinput"
            path.write_text("<mesh>\nnx1 = 16\nnx1 = 16\n", encoding="utf-8")
            with self.assertRaisesRegex(q006.AuditError, "duplicate parameter"):
                q006.parse_athinput(path)

    def test_exact_deck_binding_rejects_undeclared_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "extra.athinput"
            path.write_text(
                q006.DECKS["uniform"].read_text(encoding="utf-8")
                + "\n<problem>\nextra = true\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(q006.AuditError, "deck SHA-256 drifted"):
                q006.validate_deck("uniform", path)

    def test_invocation_rejects_noninteger_schema_version(self) -> None:
        pinned = {
            "pinned_executable_root": "/retained/pin",
            "pinned_executable_root_inventory_sha256": "1" * 64,
            "pinned_executable_path": "/retained/pin/bin/athena",
            "pinned_executable_sha256": "2" * 64,
        }
        invocation = {
            "schema_version": 1,
            "launch_style": "direct_serial_host_execution",
            "mpi_used": False,
            "slurm_used": False,
            "frontier_used": False,
            "kronos_used": False,
            "pinned_executable_root": pinned["pinned_executable_root"],
            "pinned_executable_root_inventory_sha256":
                pinned["pinned_executable_root_inventory_sha256"],
            "executable_realpath": pinned["pinned_executable_path"],
            "executable_sha256": pinned["pinned_executable_sha256"],
        }
        for value in (True, 1.0, "1"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(q006.AuditError, "schema drifted"):
                    q006._validate_serial_invocation_identity(
                        {**invocation, "schema_version": value},
                        pinned,
                        label="invocation test",
                    )

    def test_retained_json_guards_reject_boolean_and_returncode_aliases(self) -> None:
        expected = {
            "clean_empty_directory_configure": True,
            "configure_returncode": 0,
            "mpi_used": False,
        }
        for replacement in (
            {"clean_empty_directory_configure": 1},
            {"clean_empty_directory_configure": 1.0},
            {"configure_returncode": False},
            {"configure_returncode": 0.0},
            {"mpi_used": 0},
            {"mpi_used": 0.0},
        ):
            with self.subTest(replacement=replacement):
                with self.assertRaisesRegex(q006.AuditError, "primitive type drifted"):
                    q006._require_exact_match({**expected, **replacement}, expected, "metadata")
        for value in (False, 0.0):
            with self.subTest(returncode=value):
                with self.assertRaisesRegex(q006.AuditError, "returncode drifted"):
                    q006._require_exact_int(value, 0, "returncode")

    def test_parser_returncode_sidecar_requires_canonical_text(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "returncode.txt"
            for expected in (0, 1):
                with self.subTest(expected=expected):
                    path.write_text(f"{expected}\n", encoding="utf-8")
                    q006._require_canonical_returncode_sidecar(
                        path, expected, "parser returncode"
                    )
            for alias in ("00\n", "01\n", "+0\n", " 0\n", "0 \n", "0", "0\n\n"):
                with self.subTest(alias=alias):
                    path.write_text(alias, encoding="utf-8")
                    with self.assertRaisesRegex(q006.AuditError, "sidecar drifted"):
                        q006._require_canonical_returncode_sidecar(
                            path, int(alias.strip()), "parser returncode"
                        )

    def test_runtime_command_sidecar_requires_canonical_text(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "command.txt"
            argv = ["/retained/bin/athena", "-i", "/retained/deck.athinput", "time/nlim=0"]
            path.write_text(" ".join(argv) + "\n", encoding="utf-8")
            q006._require_canonical_command_sidecar(path, argv, "runtime command")
            for alias in (
                "  ".join(argv) + "\n",
                " ".join(argv),
                "THIS CONTRADICTS THE JSON INVOCATION\n",
            ):
                with self.subTest(alias=alias):
                    path.write_text(alias, encoding="utf-8")
                    with self.assertRaisesRegex(q006.AuditError, "sidecar drifted"):
                        q006._require_canonical_command_sidecar(path, argv, "runtime command")

    def test_parser_command_json_sidecar_requires_canonical_serialization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "command.json"
            argv = ["/retained/bin/athena", "-i", "/retained/deck.athinput", "time/nlim=0"]
            path.write_text(json.dumps(argv, indent=2) + "\n", encoding="utf-8")
            q006._require_canonical_json_argv_sidecar(path, argv, "parser command sidecar")
            for alias in (
                json.dumps(argv) + "\n",
                json.dumps(argv, indent=4) + "\n",
                json.dumps(argv, indent=2),
            ):
                with self.subTest(alias=alias):
                    path.write_text(alias, encoding="utf-8")
                    with self.assertRaisesRegex(q006.AuditError, "drifted"):
                        q006._require_canonical_json_argv_sidecar(
                            path, argv, "parser command sidecar"
                        )

    def test_q006_freeze_receipt_requires_closed_shape_and_unique_json_keys(self) -> None:
        expected = {
            "schema_version": 1,
            "artifact_role": q006.ARTIFACT_ROLE,
            "qualification_effect": q006.QUALIFICATION_EFFECT,
            "inventory_excludes": q006.INVENTORY_NAME,
            "freeze_policy": "remove all owner, group and other write bits recursively",
        }
        with self.assertRaisesRegex(q006.AuditError, "object keys drifted"):
            q006._require_exact_match({**expected, "unexpected": "alias"}, expected, "receipt")
        with self.assertRaisesRegex(q006.AuditError, "duplicated"):
            q006._load_retained_json_text(
                '{"schema_version": 1, "schema_version": 1}',
                "Q-006 freeze receipt",
            )

    def test_pvtk_execution_header_rejects_spacing_and_line_ending_aliases(self) -> None:
        canonical = (
            b"# vtk DataFile Version 2.0\n"
            b"# AthenaK particle data at time= 0  nranks= 1  cycle=0  "
            b"variables=prtcl_all\nBINARY\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.vtk"
            path.write_bytes(canonical)
            self.assertEqual(q006._read_pvtk_execution_metadata(path)["cycle"], 0)
            for alias in (
                canonical.replace(b"\n", b"\r\n"),
                canonical.replace(b"time= 0", b"time=  0"),
                canonical.replace(b"nranks= 1", b"nranks=  1"),
                canonical.replace(b"nranks= 1", b"nranks= 01"),
                canonical.replace(b"cycle=0", b"cycle= 0"),
                canonical.replace(b"cycle=0", b"cycle=00"),
                canonical.replace(b"variables=prtcl_all", b"variables=prtcl_all "),
            ):
                with self.subTest(alias=alias):
                    path.write_bytes(alias)
                    with self.assertRaisesRegex(q006.AuditError, "header|variables"):
                        q006._read_pvtk_execution_metadata(path)

    def test_particle_vtk_reader_rejects_noncanonical_counts_and_line_endings(self) -> None:
        canonical = (
            b"# vtk DataFile Version 2.0\n"
            b"# fixture\n"
            b"BINARY\n"
            b"DATASET UNSTRUCTURED_GRID\n"
            b"\nPOINTS 1 float\n"
            + struct.pack(">3f", 0.0, 0.0, 0.0)
            + b"\n\nPOINT_DATA 1\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.vtk"
            path.write_bytes(canonical)
            self.assertEqual(read_particle_vtk(path).points.shape, (1, 3))
            for alias in (
                canonical.replace(b"\n", b"\r\n"),
                canonical.replace(b"POINTS 1", b"POINTS 01"),
                canonical.replace(b"POINTS 1", b"POINTS  1"),
                canonical.replace(b"POINT_DATA 1", b"POINT_DATA 01"),
                canonical.replace(b"POINT_DATA 1", b"POINT_DATA  1"),
            ):
                with self.subTest(alias=alias):
                    path.write_bytes(alias)
                    with self.assertRaises(ValueError):
                        read_particle_vtk(path)

    def test_readiness_sidecar_binds_only_new_q006_runtime_local_files(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-006")
        self.assertEqual(sidecar["qualification_effect"], "none")
        self.assertFalse(sidecar["claim_closure"])
        expected_paths = {
            "src/pgen/tests/q006_paper_multispecies_oscillation_runtime_local.cpp",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_smr_runtime_local.athinput",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_audited_amr_runtime_local.athinput",
            "tst/publication/immutable_orion_tree.py",
            "tst/publication/analyze_q006_paper_multispecies_oscillation_runtime_local.py",
            "tst/publication/test_analyze_q006_paper_multispecies_oscillation_runtime_local.py",
        }
        self.assertEqual(set(sidecar["artifact_bindings"]), expected_paths)
        for relative, expected in sidecar["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected)
        shared_paths = {
            "docs/source/engineering/pic_mhd_model_contract.md",
            "src/CMakeLists.txt",
            "src/particles/particles.cpp",
            "src/pgen/pgen.cpp",
            "src/pgen/pgen.hpp",
            "tst/publication/pvtk_particles.py",
            "tst/scripts/particles/pic_parser_contract_guards.py",
        }
        self.assertEqual(set(sidecar["shared_artifact_bindings"]), shared_paths)
        for relative, expected in sidecar["shared_artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected)
        runtime = sidecar["bounded_runtime_probe"]
        self.assertTrue(runtime["tree_freeze"]["recursively_read_only"])
        self.assertEqual(runtime["tree_freeze"]["writable_entries"], 0)
        for text in ("long-horizon", "MPI", "GPU", "Frontier", "external review"):
            self.assertTrue(any(text in item for item in sidecar["explicitly_not_claimed"]))


if __name__ == "__main__":
    unittest.main()
