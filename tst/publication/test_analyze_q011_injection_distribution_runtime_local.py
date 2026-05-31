#!/usr/bin/env python3
"""Focused tests for the bounded local Q-011 injection runtime audit."""

from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import stat
import tarfile
import tempfile
import unittest
from unittest import mock

import numpy as np

from tst.publication import analyze_q011_injection_distribution_runtime_local as q011
from tst.publication import immutable_orion_tree
from tst.publication.pvtk_particles import ParticleVTKData


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_injection_distribution_runtime_local_successor_v4_2026-05-31.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _make_writable_tree(tree: Path) -> None:
    for path in [tree, *tree.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


def _payload(count: int = 512) -> ParticleVTKData:
    indices = np.arange(count, dtype=np.float64)
    mu = -1.0 + 2.0 * (indices + 0.5) / count
    phi = 2.0 * np.pi * np.mod(indices * 0.6180339887498949, 1.0)
    st = np.sqrt(1.0 - mu * mu)
    direction = np.column_stack((mu, st * np.cos(phi), st * np.sin(phi)))
    deck = q011.validate_deck()
    relative_velocity = float(deck["injection_speed_relative_to_surface"]) * direction
    surface_speed = float(deck["shock_speed"])
    light_speed = float(deck["numerical_light_speed"])
    gamma_surface = 1.0 / np.sqrt(1.0 - (surface_speed / light_speed) ** 2)
    denominator = 1.0 + surface_speed * relative_velocity[:, 0] / light_speed**2
    velocity = relative_velocity.copy()
    velocity[:, 0] = (relative_velocity[:, 0] + surface_speed) / denominator
    velocity[:, 1] = relative_velocity[:, 1] / (gamma_surface * denominator)
    velocity[:, 2] = relative_velocity[:, 2] / (gamma_surface * denominator)
    points = np.column_stack((
        np.full(count, float(deck["clamped_surface_output_x1"])),
        np.mod(indices * 0.73, 16.0),
        np.zeros(count),
    ))
    return ParticleVTKData(
        points=points,
        scalars={
            "gid": np.zeros(count, dtype=np.int64),
            "species": np.zeros(count, dtype=np.int64),
            "ptag": np.arange(count, dtype=np.int64),
            "cr_source": np.ones(count, dtype=np.int64),
            "macro_weight": np.ones(count),
            "birth_time": np.zeros(count),
            "deltaf_f0": np.ones(count),
            "deltaf_weight": np.zeros(count),
        },
        vectors={"vel": velocity},
    )


def _emitted_payload(count: int = 512, timestep: float = 0.5) -> ParticleVTKData:
    tags = np.arange(count, dtype=np.int64)
    replay = q011._expected_cycle_one_payload(tags, timestep)
    return ParticleVTKData(
        points=replay["points"].astype(np.float32).astype(np.float64),
        scalars={
            "gid": np.zeros(count, dtype=np.int64),
            "species": np.zeros(count, dtype=np.int64),
            "ptag": tags,
            "cr_source": np.ones(count, dtype=np.int64),
            "macro_weight": np.ones(count),
            "birth_time": np.zeros(count),
            "deltaf_f0": np.ones(count),
            "deltaf_weight": np.zeros(count),
        },
        vectors={"vel": replay["velocity"].astype(np.float32).astype(np.float64)},
    )


class Q011InjectionDistributionRuntimeLocalTests(unittest.TestCase):
    def test_deck_and_committed_source_bind_reduced_nonqualifying_sampler(self) -> None:
        deck = q011.validate_deck()
        source = q011.validate_source_contract()
        self.assertEqual(deck["audit_role"], q011.ARTIFACT_ROLE)
        self.assertEqual(deck["qualification_effect"], "none")
        self.assertFalse(deck["qualifying_evidence"])
        self.assertEqual(deck["audit_method"], q011.AUDIT_METHOD)
        self.assertAlmostEqual(deck["shock_speed"], 10.00000000005)
        self.assertEqual(
            source["sampler"],
            "full_sphere_isotropic_monomomentum_relative_to_ideal_surface",
        )

    def test_payload_audits_surface_monoenergetic_shell_and_full_sphere_bounds(self) -> None:
        report = q011.analyze_particle_payload(_payload())
        self.assertFalse(report["qualifying_evidence"])
        self.assertEqual(report["particle_count"], 512)
        self.assertEqual(report["provenance"]["cr_source"], "shock_injected")
        self.assertTrue(
            report["shock_surface_placement"]["single_surface_coordinate"]
        )
        self.assertEqual(
            report["shock_surface_placement"]["maximum_absolute_x1_residual"],
            0.0,
        )
        sampler = report["monoenergetic_full_sphere_sampler"]
        self.assertLess(sampler["maximum_absolute_relative_speed_residual"], 1.0e-12)
        self.assertLess(
            sampler["maximum_absolute_relative_momentum_residual"], 1.0e-12
        )
        self.assertTrue(all(count > 0 for count in sampler["octant_counts"]))

    def test_cycle_one_payload_replays_transport_and_reflecting_wall(self) -> None:
        report = q011.analyze_particle_payload(
            _emitted_payload(),
            cycle_one_timestep=0.5,
        )
        replay = report["exact_cycle_one_emission_replay"]
        self.assertGreater(replay["reflected_particle_count"], 0)
        self.assertTrue(replay["serialized_float32_payload_exact"])
        self.assertEqual(replay["maximum_absolute_point_residual"], 0.0)
        self.assertEqual(replay["maximum_absolute_velocity_residual"], 0.0)
        self.assertEqual(
            report["shock_surface_placement"]["payload_state"],
            "emitted_after_one_push_with_boundary_handling",
        )

    def test_payload_rejects_planar_drift_surface_drift_and_wrong_provenance(self) -> None:
        payload = _payload()
        cases = []
        planar = ParticleVTKData(
            points=payload.points,
            scalars=payload.scalars,
            vectors={"vel": payload.vectors["vel"].copy()},
        )
        planar.vectors["vel"][:, 2] = 0.0
        cases.append(("monoenergetic shell", planar))
        shifted = ParticleVTKData(
            points=payload.points.copy(),
            scalars=payload.scalars,
            vectors=payload.vectors,
        )
        shifted.points[0, 0] += 1.0e-4
        cases.append(("shock surface", shifted))
        wrong_source = {**payload.scalars, "cr_source": payload.scalars["cr_source"].copy()}
        wrong_source["cr_source"][0] = 0
        cases.append((
            "shock_injected",
            ParticleVTKData(payload.points, wrong_source, payload.vectors),
        ))
        for expected, candidate in cases:
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(q011.AuditError, expected):
                    q011.analyze_particle_payload(candidate)

    def test_immutable_extractor_requires_verified_runtime_and_executable_trees(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            runtime_tree = base / "runtime"
            executable_tree = base / "executable"
            runtime_tree.mkdir()
            executable_tree.mkdir()
            pvtk = runtime_tree / "snapshot.part.vtk"
            executable = executable_tree / "athena"
            pvtk_payload = (
                b"# vtk DataFile Version 2.0\n"
                b"# AthenaK particle data at time= 0.5  nranks= 1  cycle=1  "
                b"variables=prtcl_all\nBINARY\nparticle artifact"
            )
            pvtk.write_bytes(pvtk_payload)
            executable.write_bytes(b"host executable")
            q011.ORION_BULK_ROOT = base
            try:
                runtime_freeze = q011.freeze_tree(runtime_tree)
                executable_freeze = q011.freeze_tree(executable_tree)
                with mock.patch.object(
                    q011,
                    "read_particle_vtk",
                    return_value=_emitted_payload(),
                ):
                    report = q011.extract_runtime_artifact(
                        pvtk,
                        artifact_root=runtime_tree,
                        artifact_inventory_sha256=runtime_freeze["inventory_sha256"],
                        executable=executable,
                        executable_root=executable_tree,
                        executable_root_inventory_sha256=executable_freeze["inventory_sha256"],
                        expected_executable_sha256=hashlib.sha256(
                            b"host executable"
                        ).hexdigest(),
                    )
            finally:
                _make_writable_tree(runtime_tree)
                _make_writable_tree(executable_tree)
                q011.ORION_BULK_ROOT = original_root
        self.assertEqual(
            report["immutable_runtime_artifact"]["sha256"],
            hashlib.sha256(pvtk_payload).hexdigest(),
        )
        self.assertEqual(
            report["host_executable"]["sha256"],
            hashlib.sha256(b"host executable").hexdigest(),
        )
        self.assertEqual(
            report["pvtk_execution_metadata"],
            {"time": 0.5, "nranks": 1, "cycle": 1, "variables": "prtcl_all"},
        )
        self.assertTrue(report["immutable_runtime_tree"]["recursively_read_only"])
        self.assertTrue(report["immutable_executable_tree"]["recursively_read_only"])

    def test_immutable_extractor_rejects_unfrozen_payload_and_executable_drift(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            runtime_tree = base / "runtime"
            executable_tree = base / "executable"
            outside = base / "outside"
            runtime_tree.mkdir()
            executable_tree.mkdir()
            outside.mkdir()
            pvtk = runtime_tree / "snapshot.part.vtk"
            executable = executable_tree / "athena"
            forged = outside / "forged.part.vtk"
            payload = (
                b"# vtk DataFile Version 2.0\n"
                b"# AthenaK particle data at time= 0.5  nranks= 1  cycle=1  "
                b"variables=prtcl_all\nBINARY\nparticle artifact"
            )
            pvtk.write_bytes(payload)
            forged.write_bytes(payload)
            executable.write_bytes(b"host executable")
            q011.ORION_BULK_ROOT = base
            try:
                runtime_freeze = q011.freeze_tree(runtime_tree)
                executable_freeze = q011.freeze_tree(executable_tree)
                arguments = {
                    "artifact_root": runtime_tree,
                    "artifact_inventory_sha256": runtime_freeze["inventory_sha256"],
                    "executable": executable,
                    "executable_root": executable_tree,
                    "executable_root_inventory_sha256": executable_freeze[
                        "inventory_sha256"
                    ],
                    "expected_executable_sha256": hashlib.sha256(
                        b"host executable"
                    ).hexdigest(),
                }
                with self.assertRaisesRegex(q011.AuditError, "must remain below"):
                    q011.extract_runtime_artifact(forged, **arguments)
                with self.assertRaisesRegex(q011.AuditError, "executable SHA-256 drifted"):
                    q011.extract_runtime_artifact(
                        pvtk,
                        **{
                            **arguments,
                            "expected_executable_sha256": "0" * 64,
                        },
                    )
            finally:
                _make_writable_tree(runtime_tree)
                _make_writable_tree(executable_tree)
                q011.ORION_BULK_ROOT = original_root

    def test_immutable_extractor_rejects_payload_mutation_after_initial_verification(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            runtime_tree = base / "runtime"
            executable_tree = base / "executable"
            runtime_tree.mkdir()
            executable_tree.mkdir()
            pvtk = runtime_tree / "snapshot.part.vtk"
            executable = executable_tree / "athena"
            pvtk.write_bytes(
                b"# vtk DataFile Version 2.0\n"
                b"# AthenaK particle data at time= 0.5  nranks= 1  cycle=1  "
                b"variables=prtcl_all\nBINARY\nparticle artifact"
            )
            executable.write_bytes(b"host executable")
            q011.ORION_BULK_ROOT = base
            try:
                runtime_freeze = q011.freeze_tree(runtime_tree)
                executable_freeze = q011.freeze_tree(executable_tree)
                original_verify = immutable_orion_tree._verify_frozen_tree_anchored
                verify_count = 0

                def mutate_runtime_payload(*args: object, **kwargs: object):
                    nonlocal verify_count
                    verified = original_verify(*args, **kwargs)
                    verify_count += 1
                    if verify_count == 1:
                        pvtk.chmod(pvtk.stat().st_mode | stat.S_IWUSR)
                        pvtk.write_bytes(b"substituted payload")
                        pvtk.chmod(pvtk.stat().st_mode & ~stat.S_IWUSR)
                    return verified

                with mock.patch.object(
                    immutable_orion_tree,
                    "_verify_frozen_tree_anchored",
                    side_effect=mutate_runtime_payload,
                ):
                    with self.assertRaisesRegex(
                        q011.AuditError, "(SHA-256|artifact hash) drifted"
                    ):
                        q011.extract_runtime_artifact(
                            pvtk,
                            artifact_root=runtime_tree,
                            artifact_inventory_sha256=runtime_freeze["inventory_sha256"],
                            executable=executable,
                            executable_root=executable_tree,
                            executable_root_inventory_sha256=executable_freeze[
                                "inventory_sha256"
                            ],
                            expected_executable_sha256=hashlib.sha256(
                                b"host executable"
                            ).hexdigest(),
                        )
            finally:
                _make_writable_tree(runtime_tree)
                _make_writable_tree(executable_tree)
                q011.ORION_BULK_ROOT = original_root

    def test_inventory_freeze_and_verify_are_orion_scoped(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            (tree / "stdout.txt").write_text("bounded\n", encoding="utf-8")
            (tree / "pvtk").mkdir()
            (tree / "pvtk/snapshot.vtk").write_bytes(b"snapshot")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                report = q011.freeze_tree(tree)
                self.assertEqual(report["inventoried_file_count"], 3)
                self.assertTrue(report["recursively_read_only"])
                self.assertEqual(report["writable_entries"], [])
                verified = q011.verify_frozen_tree(tree, report["inventory_sha256"])
                self.assertEqual(
                    verified["inventory_sha256"],
                    report["inventory_sha256"],
                )
                self.assertEqual(verified["inventoried_file_count"], 3)
                self.assertTrue(verified["recursively_read_only"])
            finally:
                _make_writable_tree(tree)
                q011.ORION_BULK_ROOT = original_root

    def test_payload_rejects_lossy_integer_provenance(self) -> None:
        payload = _payload()
        scalars = {**payload.scalars, "cr_source": payload.scalars["cr_source"].astype(float)}
        scalars["cr_source"][0] = 1.5
        with self.assertRaisesRegex(q011.AuditError, "exact integers"):
            q011.analyze_particle_payload(
                ParticleVTKData(payload.points, scalars, payload.vectors)
            )

    def test_payload_rejects_weight_and_deltaf_metadata_drift(self) -> None:
        payload = _payload()
        for name, value, expected in (
            ("gid", 1, "carrier gid zero"),
            ("macro_weight", 2.0, "macro weights must be one"),
            ("deltaf_f0", 2.0, "delta-f f0 values must be one"),
            ("deltaf_weight", 1.0, "delta-f weights must be zero"),
        ):
            with self.subTest(name=name):
                scalars = {**payload.scalars, name: payload.scalars[name].copy()}
                scalars[name][0] = value
                with self.assertRaisesRegex(q011.AuditError, expected):
                    q011.analyze_particle_payload(
                        ParticleVTKData(payload.points, scalars, payload.vectors)
                    )

    def test_cycle_one_payload_rejects_single_float32_ulp_drift(self) -> None:
        payload = _emitted_payload()
        points = payload.points.copy()
        points[0, 1] = np.nextafter(
            np.float32(points[0, 1]), np.float32(np.inf), dtype=np.float32
        )
        with self.assertRaisesRegex(q011.AuditError, "float32 points drifted"):
            q011.analyze_particle_payload(
                ParticleVTKData(points, payload.scalars, payload.vectors),
                cycle_one_timestep=0.5,
            )

    def test_runtime_header_rejects_wrong_cycle_rank_and_variable(self) -> None:
        mutations = (
            (b"nranks= 1", b"nranks= 2", "nranks=1"),
            (b"cycle=1", b"cycle=2", "cycle=1"),
            (b"variables=prtcl_all", b"variables=prtcl_all_extra", "variables=prtcl_all"),
        )
        for old, new, expected in mutations:
            with self.subTest(mutation=new):
                with tempfile.TemporaryDirectory() as directory:
                    pvtk = Path(directory) / "snapshot.part.vtk"
                    header = (
                        b"# vtk DataFile Version 2.0\n"
                        b"# AthenaK particle data at time= 0.5  nranks= 1  cycle=1  "
                        b"variables=prtcl_all\nBINARY\n"
                    )
                    pvtk.write_bytes(header.replace(old, new))
                    with self.assertRaisesRegex(q011.AuditError, expected):
                        q011._read_pvtk_execution_metadata(pvtk)

    def test_runtime_header_rejects_spacing_and_line_ending_aliases(self) -> None:
        canonical = (
            b"# vtk DataFile Version 2.0\n"
            b"# AthenaK particle data at time= 0.5  nranks= 1  cycle=1  "
            b"variables=prtcl_all\nBINARY\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            pvtk = Path(directory) / "snapshot.part.vtk"
            pvtk.write_bytes(canonical)
            self.assertEqual(q011._read_pvtk_execution_metadata(pvtk)["cycle"], 1)
            for alias in (
                canonical.replace(b"\n", b"\r\n"),
                canonical.replace(b"time= 0.5", b"time=  0.5"),
                canonical.replace(b"nranks= 1", b"nranks=  1"),
                canonical.replace(b"nranks= 1", b"nranks= 01"),
                canonical.replace(b"cycle=1", b"cycle= 1"),
                canonical.replace(b"cycle=1", b"cycle=01"),
                canonical.replace(b"variables=prtcl_all", b"variables=prtcl_all "),
            ):
                with self.subTest(alias=alias):
                    pvtk.write_bytes(alias)
                    with self.assertRaisesRegex(q011.AuditError, "header|variables"):
                        q011._read_pvtk_execution_metadata(pvtk)

    def test_freezer_rejects_reserved_symlink_before_outside_overwrite(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            outside = Path(directory) / "outside.txt"
            outside.write_text("preserve\n", encoding="utf-8")
            (tree / q011.FREEZE_RECEIPT_NAME).symlink_to(outside)
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                with self.assertRaisesRegex(q011.AuditError, "symlink"):
                    q011.freeze_tree(tree)
                self.assertEqual(outside.read_text(encoding="utf-8"), "preserve\n")
            finally:
                q011.ORION_BULK_ROOT = original_root

    def test_freezer_rejects_existing_inventory_before_receipt_write(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            (tree / q011.INVENTORY_NAME).write_text("preserve\n", encoding="utf-8")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                with self.assertRaisesRegex(q011.AuditError, "reserved metadata"):
                    q011.freeze_tree(tree)
                self.assertFalse((tree / q011.FREEZE_RECEIPT_NAME).exists())
                self.assertEqual(
                    (tree / q011.INVENTORY_NAME).read_text(encoding="utf-8"),
                    "preserve\n",
                )
            finally:
                q011.ORION_BULK_ROOT = original_root

    def test_freezer_rejects_line_separator_path_before_metadata_or_chmod(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            payload = tree / "payload\rname"
            payload.write_text("preserve\n", encoding="utf-8")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                with self.assertRaisesRegex(q011.AuditError, "line separators"):
                    q011.freeze_tree(tree)
                self.assertFalse((tree / q011.FREEZE_RECEIPT_NAME).exists())
                self.assertFalse((tree / q011.INVENTORY_NAME).exists())
                self.assertTrue(tree.stat().st_mode & stat.S_IWUSR)
                self.assertTrue(payload.stat().st_mode & stat.S_IWUSR)
            finally:
                q011.ORION_BULK_ROOT = original_root

    def test_frozen_verify_rejects_semantically_dishonest_receipt(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            (tree / "payload").write_text("bounded\n", encoding="utf-8")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                q011.freeze_tree(tree)
                _make_writable_tree(tree)
                receipt = tree / q011.FREEZE_RECEIPT_NAME
                payload = json.loads(receipt.read_text(encoding="utf-8"))
                payload["freeze_policy"] = "dishonest"
                receipt.write_text(
                    json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                inventory = tree / q011.INVENTORY_NAME
                lines = inventory.read_text(encoding="utf-8").splitlines()
                inventory.write_text(
                    "\n".join(
                        f"{_sha256(receipt)}  {q011.FREEZE_RECEIPT_NAME}"
                        if line.endswith(f"  {q011.FREEZE_RECEIPT_NAME}")
                        else line
                        for line in lines
                    )
                    + "\n",
                    encoding="utf-8",
                )
                for path in reversed([tree, *tree.rglob("*")]):
                    path.chmod(path.stat().st_mode & ~stat.S_IWUSR)
                with self.assertRaisesRegex(q011.AuditError, "freeze receipt freeze_policy"):
                    q011.verify_frozen_tree(tree, _sha256(inventory))
            finally:
                _make_writable_tree(tree)
                q011.ORION_BULK_ROOT = original_root

    def test_freezer_orders_rendered_relative_paths_lexically(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            (tree / "a-b").write_text("flat\n", encoding="utf-8")
            (tree / "a").mkdir()
            (tree / "a/b").write_text("nested\n", encoding="utf-8")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                q011.freeze_tree(tree)
                paths = [
                    line.partition("  ")[2]
                    for line in (tree / q011.INVENTORY_NAME).read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                self.assertEqual(paths, sorted(paths))
            finally:
                _make_writable_tree(tree)
                q011.ORION_BULK_ROOT = original_root

    def test_freezer_rejects_invalid_receipt_before_metadata_or_chmod(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            payload = tree / "payload"
            payload.write_text("preserve\n", encoding="utf-8")
            with self.assertRaisesRegex(q011.AuditError, "freeze receipt freeze_policy"):
                immutable_orion_tree.freeze_tree(
                    tree,
                    {
                        "schema_version": 1,
                        "artifact_role": "bounded",
                        "qualification_effect": "none",
                        "inventory_excludes": q011.INVENTORY_NAME,
                        "freeze_policy": "dishonest",
                    },
                    authorized_root=base,
                    error_type=q011.AuditError,
                    label="test invalid receipt",
                )
            self.assertFalse((tree / q011.FREEZE_RECEIPT_NAME).exists())
            self.assertFalse((tree / q011.INVENTORY_NAME).exists())
            self.assertTrue(tree.stat().st_mode & stat.S_IWUSR)
            self.assertTrue(payload.stat().st_mode & stat.S_IWUSR)

    def test_frozen_verify_rejects_alias_symlink_fifo_and_hardlink(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            q011.ORION_BULK_ROOT = base
            try:
                for entry_type in ("root_alias", "symlink", "fifo", "hardlink"):
                    with self.subTest(entry_type=entry_type):
                        tree = base / entry_type
                        tree.mkdir()
                        payload = tree / "payload"
                        payload.write_text("bounded\n", encoding="utf-8")
                        report = q011.freeze_tree(tree)
                        digest = report["inventory_sha256"]
                        if entry_type == "root_alias":
                            alias = base / "alias"
                            alias.symlink_to(tree, target_is_directory=True)
                            with self.assertRaisesRegex(q011.AuditError, "canonical"):
                                q011.verify_frozen_tree(alias, digest)
                            alias.unlink()
                        else:
                            tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                            if entry_type == "symlink":
                                (tree / "alias").symlink_to(payload)
                                expected = "symlink"
                            elif entry_type == "fifo":
                                os.mkfifo(tree / "fifo")
                                expected = "unsupported"
                            else:
                                os.link(payload, tree / "linked")
                                expected = "hard-linked"
                            tree.chmod(tree.stat().st_mode & ~stat.S_IWUSR)
                            with self.assertRaisesRegex(q011.AuditError, expected):
                                q011.verify_frozen_tree(tree, digest)
                        _make_writable_tree(tree)
            finally:
                q011.ORION_BULK_ROOT = original_root

    def test_freezer_rejects_alias_symlink_fifo_and_hardlink_before_mutation(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            q011.ORION_BULK_ROOT = base
            try:
                for entry_type in ("root_alias", "symlink", "fifo", "hardlink"):
                    with self.subTest(entry_type=entry_type):
                        tree = base / entry_type
                        tree.mkdir()
                        payload = tree / "payload"
                        payload.write_text("bounded\n", encoding="utf-8")
                        candidate = tree
                        if entry_type == "root_alias":
                            candidate = base / "alias"
                            candidate.symlink_to(tree, target_is_directory=True)
                            expected = "canonical"
                        elif entry_type == "symlink":
                            (tree / "alias").symlink_to(payload)
                            expected = "symlink"
                        elif entry_type == "fifo":
                            os.mkfifo(tree / "fifo")
                            expected = "unsupported"
                        else:
                            os.link(payload, tree / "linked")
                            expected = "hard-linked"
                        with self.assertRaisesRegex(q011.AuditError, expected):
                            q011.freeze_tree(candidate)
                        self.assertFalse((tree / q011.FREEZE_RECEIPT_NAME).exists())
                        self.assertFalse((tree / q011.INVENTORY_NAME).exists())
                        self.assertTrue(tree.stat().st_mode & stat.S_IWUSR)
                        self.assertTrue(payload.stat().st_mode & stat.S_IWUSR)
                        if entry_type == "root_alias":
                            candidate.unlink()
            finally:
                q011.ORION_BULK_ROOT = original_root

    def test_frozen_verify_rejects_writable_entry(self) -> None:
        original_root = q011.ORION_BULK_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "retained"
            tree.mkdir()
            payload = tree / "payload"
            payload.write_text("bounded\n", encoding="utf-8")
            q011.ORION_BULK_ROOT = Path(directory)
            try:
                report = q011.freeze_tree(tree)
                payload.chmod(payload.stat().st_mode | stat.S_IWUSR)
                with self.assertRaisesRegex(q011.AuditError, "writable"):
                    q011.verify_frozen_tree(tree, report["inventory_sha256"])
            finally:
                _make_writable_tree(tree)
                q011.ORION_BULK_ROOT = original_root

    def test_source_archive_validator_rejects_bytecode_cache_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "source.tar"
            with tarfile.open(archive_path, "w") as archive:
                payload = b"cache\n"
                member = tarfile.TarInfo("src/__pycache__/module.cpython-311.pyc")
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
            with self.assertRaisesRegex(q011.AuditError, "bytecode cache"):
                immutable_orion_tree.validate_source_archive(
                    archive_path,
                    _sha256(archive_path),
                    error_type=q011.AuditError,
                    label="test source archive",
                )

    def test_source_archive_validator_rejects_noncanonical_member_name(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "source.tar"
            with tarfile.open(archive_path, "w") as archive:
                payload = b"source\n"
                member = tarfile.TarInfo("src//module.py")
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
            with self.assertRaisesRegex(q011.AuditError, "unsafe"):
                immutable_orion_tree.validate_source_archive(
                    archive_path,
                    _sha256(archive_path),
                    error_type=q011.AuditError,
                    label="test source archive",
                )

    def test_executable_validator_rejects_non_elf_payload(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "athena"
            executable.write_bytes(b"not an executable\n")
            executable.chmod(0o555)
            with self.assertRaisesRegex(q011.AuditError, "not an ELF"):
                immutable_orion_tree.validate_executable_elf(
                    executable,
                    _sha256(executable),
                    error_type=q011.AuditError,
                    label="test executable",
                )

    def test_source_archive_dependencies_reject_unrelated_archive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "source.tar"
            with tarfile.open(archive_path, "w") as archive:
                payload = b"unrelated\n"
                member = tarfile.TarInfo("src/unrelated.cpp")
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
            with self.assertRaisesRegex(q011.AuditError, "absent"):
                immutable_orion_tree.validate_source_archive_dependencies(
                    archive_path,
                    {"src/required.cpp": hashlib.sha256(b"required\n").hexdigest()},
                    {"src/required.cpp"},
                    error_type=q011.AuditError,
                    label="test source archive dependencies",
                )

    def test_strict_deck_parser_rejects_duplicate_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.athinput"
            path.write_text("<mesh>\nnx1 = 8\nnx1 = 8\n", encoding="utf-8")
            with self.assertRaisesRegex(q011.AuditError, "duplicate parameter"):
                q011.parse_athinput(path)

    def test_exact_deck_binding_rejects_undeclared_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "extra.athinput"
            path.write_text(
                q011.DECK.read_text(encoding="utf-8") + "\n<problem>\nextra = true\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(q011.AuditError, "deck SHA-256 drifted"):
                q011.validate_deck(path)

    def test_source_binding_rejects_digest_drift_even_when_snippets_remain(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pic_parallel_shock.cpp"
            path.write_text(
                q011.SOURCE.read_text(encoding="utf-8") + "\n// digest drift\n",
                encoding="utf-8",
            )
            with mock.patch.object(q011, "SOURCE", path):
                with self.assertRaisesRegex(q011.AuditError, "source SHA-256 drifted"):
                    q011.validate_source_contract()

    def test_readiness_sidecar_binds_only_new_q011_files_and_nonclaim_gaps(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-011")
        self.assertEqual(sidecar["qualification_effect"], "none")
        self.assertFalse(sidecar["claim_closure"])
        expected_paths = {
            "inputs/tests/pic_q011_injection_distribution_runtime_local.athinput",
            "tst/publication/immutable_orion_tree.py",
            "tst/publication/analyze_q011_injection_distribution_runtime_local.py",
            "tst/publication/test_analyze_q011_injection_distribution_runtime_local.py",
        }
        self.assertEqual(set(sidecar["artifact_bindings"]), expected_paths)
        for relative, expected_sha256 in sidecar["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        gaps = sidecar["explicit_gaps"]
        for word in ("calibration", "AMR", "GPU", "external review"):
            self.assertTrue(any(word in gap for gap in gaps), word)
        runtime = sidecar["bounded_runtime_probe"]
        self.assertTrue(runtime["tree_freeze"]["recursively_read_only"])
        self.assertEqual(runtime["tree_freeze"]["writable_entries"], 0)


if __name__ == "__main__":
    unittest.main()
