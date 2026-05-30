#!/usr/bin/env python3
"""Regression tests for the structured Frontier F1 offline analyzers."""

from __future__ import annotations

import base64
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
import unittest

PUBLICATION_DIR = Path(__file__).resolve().parent
CONTROL_PLANE_DIR = PUBLICATION_DIR / "frontier_control_plane"
import sys
sys.path.insert(0, str(PUBLICATION_DIR))
sys.path.insert(0, str(CONTROL_PLANE_DIR))

from control_plane_common import PRODUCTION_RUNTIME_LOADED_MODULES
from control_plane_common import PRODUCTION_RUNTIME_MODULEFILES
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH
from frontier_f1_gpu_paper_coupling_analysis import parse_particle_vtk_bytes as parse_coupling_vtk
from frontier_f1_gpu_paper_coupling_analysis import registered_case_output_paths
from frontier_f1_gpu_paper_coupling_analysis import require_execution_context
from frontier_f1_gpu_relativistic_gyro_analysis import EXPECTED_BASELINE
from frontier_f1_gpu_relativistic_gyro_analysis import RUNTIME_ALLOWLIST_KEYS
from frontier_f1_gpu_relativistic_gyro_analysis import _structured_execution_context
from frontier_f1_gpu_relativistic_gyro_analysis import parse_particle_vtk_bytes as parse_gyro_vtk
from frontier_f1_gpu_relativistic_gyro_analysis import registered_particle_output
from frontier_f1_gpu_relativistic_gyro_analysis import require_linked_library
from frontier_f1_structured_artifacts import StructuredArtifactTree
from frontier_f1_structured_artifacts import REVIEWED_FRONTIER_MPICH_DIAGNOSTIC_SHA256
from frontier_f1_structured_artifacts import load_inventory
from frontier_f1_structured_artifacts import read_inventory_bytes
from frontier_f1_structured_artifacts import require_inventory_sha256
from frontier_f1_structured_artifacts import validate_frontier_mpich_diagnostic_stderr
from frontier_f1_structured_artifacts import write_result_exclusive


class FrontierF1StructuredAnalysisTests(unittest.TestCase):
    def test_mpich_stderr_accepts_only_exact_reviewed_transcript(self) -> None:
        diagnostic = base64.b64decode(
            (
                PUBLICATION_DIR
                / "readiness/frontier_mpich_diagnostic_stderr.base64"
            ).read_bytes().replace(b"\n", b""),
            validate=True,
        ).decode("utf-8")
        self.assertEqual(
            hashlib.sha256(diagnostic.encode("utf-8")).hexdigest(),
            REVIEWED_FRONTIER_MPICH_DIAGNOSTIC_SHA256,
        )
        validate_frontier_mpich_diagnostic_stderr(diagnostic)
        for forged in (
            diagnostic + "warning: runtime drift\n",
            diagnostic + "PE 0:   MALICIOUS_DIAGNOSTIC = runtime drift accepted\n",
            diagnostic.replace("CRAY MPICH version 8.1.31.9",
                               "CRAY MPICH version 8.1.31.9 arbitrary suffix"),
            diagnostic.replace("(CH4)\n", "(CH4) arbitrary suffix\n"),
            diagnostic + "PE 0:   MPICH_GPU_SUPPORT_ENABLED                      = 1\n",
            diagnostic.replace("MPICH_GPU_SUPPORT_ENABLED                      = 1",
                               "MPICH_GPU_SUPPORT_ENABLED                      = 0"),
            diagnostic.replace("MPI BUILD INFO : Wed", "MPI BUILD INFO : \nWed"),
            "",
        ):
            with self.subTest(forged=forged):
                with self.assertRaises(ValueError):
                    validate_frontier_mpich_diagnostic_stderr(forged)

    def _values(self) -> dict[str, str]:
        return {
            **{key: "<unset>" for key in RUNTIME_ALLOWLIST_KEYS},
            **EXPECTED_BASELINE,
            "LOADEDMODULES": ":".join(PRODUCTION_RUNTIME_LOADED_MODULES),
            "_LMFILES_": ":".join(PRODUCTION_RUNTIME_MODULEFILES),
            "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH,
        }

    def _allowlist(
        self,
        root: Path,
        name: str,
        *,
        values: dict[str, str] | None = None,
        keys: list[str] | None = None,
        suffix: str = "",
    ) -> Path:
        path = root / name
        selected = values or self._values()
        path.write_text(
            "".join(f"{key}={selected[key]}\n" for key in (keys or RUNTIME_ALLOWLIST_KEYS))
            + suffix,
            encoding="utf-8",
        )
        return path

    def _publish_inventory(self, root: Path) -> None:
        (root / "analysis").mkdir(mode=0o700)
        (root / "analysis").chmod(0o700)
        records = []
        for path in sorted(root.rglob("*")):
            if "analysis" in path.relative_to(root).parts or not path.is_file():
                continue
            data = path.read_bytes()
            records.append(
                {
                    "path": path.relative_to(root).as_posix(),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "size": len(data),
                }
            )
            path.chmod(0o444)
        for path in sorted(root.rglob("*"), reverse=True):
            if path.is_dir() and path != root / "analysis":
                path.chmod(0o555)
        inventory_path = root / "artifact_inventory.json"
        inventory_path.write_text(
            json.dumps({"schema_version": 1, "files": records}, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        inventory_path.chmod(0o444)
        root.chmod(0o555)

    def _particle_vtk(self, *, suffix: bytes = b"") -> bytes:
        return (
            b"# AthenaK particle data at time= 0.0 nranks=1 cycle=0\n"
            b"POINTS 0 float\n"
            b"\nSCALARS gid int\nLOOKUP_TABLE default\n"
            b"\nSCALARS ptag int\nLOOKUP_TABLE default\n"
            b"\nSCALARS species int\nLOOKUP_TABLE default\n"
            b"\nSCALARS deltaf_f0 float\nLOOKUP_TABLE default\n"
            b"\nSCALARS deltaf_weight float\nLOOKUP_TABLE default\n"
            b"\nVECTORS vel float\n"
            + suffix
        )

    def test_structured_contexts_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._allowlist(root, "f1-gyro.environment.allowlist.txt")
            self._allowlist(root, "f1-coupling-coeff0.environment.allowlist.txt")
            self._allowlist(root, "f1-coupling-coeff7.environment.allowlist.txt")
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                self.assertEqual(
                    _structured_execution_context(artifact_tree, inventory)[
                        "runtime_profile"
                    ],
                    "frontier_minimum_supported",
                )
                self.assertEqual(
                    require_execution_context(artifact_tree, inventory)["runtime_profile"],
                    "frontier_minimum_supported",
                )

    def test_gyro_rejects_missing_structured_allowlist_even_with_legacy_artifacts(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "snapshot_verification.txt").write_text(
                "immutable snapshot verified\n", encoding="utf-8"
            )
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                with self.assertRaisesRegex(ValueError, "omits required path"):
                    _structured_execution_context(artifact_tree, inventory)

    def test_coupling_rejects_missing_and_partial_structured_allowlists(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                with self.assertRaisesRegex(ValueError, "incomplete"):
                    require_execution_context(artifact_tree, inventory)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._allowlist(root, "f1-coupling-coeff0.environment.allowlist.txt")
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                with self.assertRaisesRegex(ValueError, "incomplete"):
                    require_execution_context(artifact_tree, inventory)

    def test_structured_allowlist_rejects_writable_reordered_and_duplicate_records(
        self,
    ) -> None:
        variants = ("writable", "reordered", "duplicate")
        for variant in variants:
            with self.subTest(variant=variant), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                path = self._allowlist(root, "f1-gyro.environment.allowlist.txt")
                if variant == "reordered":
                    keys = [*RUNTIME_ALLOWLIST_KEYS]
                    keys[0], keys[1] = keys[1], keys[0]
                    self._allowlist(root, path.name, keys=keys)
                elif variant == "duplicate":
                    self._allowlist(
                        root,
                        path.name,
                        suffix=f"{RUNTIME_ALLOWLIST_KEYS[0]}=forged\n",
                    )
                self._publish_inventory(root)
                if variant == "writable":
                    path.chmod(0o644)
                with StructuredArtifactTree(root) as artifact_tree:
                    with self.assertRaises(ValueError):
                        inventory = load_inventory(artifact_tree)
                        _structured_execution_context(artifact_tree, inventory)

    def test_structured_allowlist_rejects_module_provenance_drift(self) -> None:
        for key in ("LOADEDMODULES", "_LMFILES_", "MODULEPATH"):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                values = self._values()
                values[key] += ":/tmp/forged"
                self._allowlist(
                    root, "f1-gyro.environment.allowlist.txt", values=values
                )
                self._publish_inventory(root)
                with StructuredArtifactTree(root) as artifact_tree:
                    inventory = load_inventory(artifact_tree)
                    with self.assertRaisesRegex(ValueError, "reviewed provenance"):
                        _structured_execution_context(artifact_tree, inventory)

    def test_inventory_rejects_checksum_drift_and_boolean_schema_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self._allowlist(root, "f1-gyro.environment.allowlist.txt")
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                path.chmod(0o644)
                path.write_bytes(path.read_bytes() + b"forged=1\n")
                path.chmod(0o444)
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    read_inventory_bytes(artifact_tree, inventory, path.name)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory_path = root / "artifact_inventory.json"
            self._publish_inventory(root)
            root.chmod(0o755)
            inventory_path.chmod(0o644)
            inventory_path.write_text(
                '{"files": [], "schema_version": true}\n', encoding="utf-8"
            )
            inventory_path.chmod(0o444)
            root.chmod(0o555)
            with StructuredArtifactTree(root) as artifact_tree:
                with self.assertRaisesRegex(ValueError, "malformed"):
                    load_inventory(artifact_tree)

    def test_parent_inventory_binding_rejects_mismatched_digest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                with self.assertRaisesRegex(ValueError, "parent binding"):
                    require_inventory_sha256(artifact_tree, "0" * 64)

    def test_result_publication_is_read_only_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._publish_inventory(root)
            output = root / "analysis" / "analysis.json"
            with StructuredArtifactTree(root) as artifact_tree:
                load_inventory(artifact_tree)
                write_result_exclusive(
                    artifact_tree, "analysis/analysis.json", {"status": "pass"}
                )
                with self.assertRaises(FileExistsError):
                    write_result_exclusive(
                        artifact_tree, "analysis/analysis.json", {"status": "forged"}
                    )
                receipt = root / "analysis" / "offline_analysis_receipt.json"
                write_result_exclusive(
                    artifact_tree,
                    "analysis/offline_analysis_receipt.json",
                    {"status": "bound"},
                )
                self.assertEqual(stat.S_IMODE(receipt.stat().st_mode), 0o444)
            self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o444)

    def test_inventory_rejects_duplicate_keys_noncanonical_paths_and_unlisted_files(
        self,
    ) -> None:
        variants = ("duplicate-key", "noncanonical-path", "unlisted-file")
        for variant in variants:
            with self.subTest(variant=variant), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                nested = root / "a"
                nested.mkdir()
                (nested / "b").write_text("verified\n", encoding="utf-8")
                self._publish_inventory(root)
                root.chmod(0o755)
                inventory_path = root / "artifact_inventory.json"
                if variant == "duplicate-key":
                    inventory_path.chmod(0o644)
                    inventory_path.write_text(
                        '{"schema_version":1,"schema_version":1,"files":[]}\n',
                        encoding="utf-8",
                    )
                    inventory_path.chmod(0o444)
                elif variant == "noncanonical-path":
                    inventory_path.chmod(0o644)
                    inventory_path.write_text(
                        inventory_path.read_text(encoding="utf-8").replace(
                            '"a/b"', '"a//b"'
                        ),
                        encoding="utf-8",
                    )
                    inventory_path.chmod(0o444)
                else:
                    extra = root / "unlisted.txt"
                    extra.write_text("forged\n", encoding="utf-8")
                    extra.chmod(0o444)
                root.chmod(0o555)
                with StructuredArtifactTree(root) as artifact_tree:
                    with self.assertRaises(ValueError):
                        load_inventory(artifact_tree)

    def test_inventory_rejects_empty_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "empty").mkdir()
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                with self.assertRaisesRegex(ValueError, "empty directory"):
                    load_inventory(artifact_tree)

    def test_pinned_tree_rejects_root_substitution_and_analysis_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "root"
            replacement = base / "replacement"
            root.mkdir()
            replacement.mkdir()
            for directory in (root, replacement):
                self._allowlist(directory, "f1-gyro.environment.allowlist.txt")
                self._publish_inventory(directory)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                root.rename(base / "detached")
                replacement.rename(root)
                with self.assertRaisesRegex(ValueError, "root path changed"):
                    read_inventory_bytes(
                        artifact_tree,
                        inventory,
                        "f1-gyro.environment.allowlist.txt",
                    )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "root"
            outside = Path(temporary) / "outside"
            root.mkdir()
            outside.mkdir()
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                load_inventory(artifact_tree)
                root.chmod(0o755)
                (root / "analysis").rename(root / "analysis-detached")
                (root / "analysis").symlink_to(outside, target_is_directory=True)
                root.chmod(0o555)
                try:
                    with self.assertRaises(ValueError):
                        write_result_exclusive(
                            artifact_tree, "analysis/analysis.json", {"status": "forged"}
                        )
                finally:
                    root.chmod(0o755)
                    (root / "analysis").unlink()
                    (root / "analysis-detached").rename(root / "analysis")

    def test_pinned_tree_rejects_nested_directory_substitution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "root"
            nested = root / "nested"
            root.mkdir()
            nested.mkdir()
            artifact = nested / "artifact.txt"
            artifact.write_text("verified\n", encoding="utf-8")
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                inventory = load_inventory(artifact_tree)
                root.chmod(0o755)
                detached = root / "nested-detached"
                nested.rename(detached)
                detached.chmod(0o755)
                (detached / "artifact.txt").unlink()
                detached.rmdir()
                nested.mkdir()
                replacement = nested / "artifact.txt"
                replacement.write_text("verified\n", encoding="utf-8")
                replacement.chmod(0o444)
                nested.chmod(0o555)
                root.chmod(0o555)
                with self.assertRaisesRegex(ValueError, "directory identity changed"):
                    read_inventory_bytes(artifact_tree, inventory, "nested/artifact.txt")

    def test_pinned_tree_rejects_analysis_directory_substitution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "root"
            root.mkdir()
            self._publish_inventory(root)
            with StructuredArtifactTree(root) as artifact_tree:
                load_inventory(artifact_tree)
                root.chmod(0o755)
                (root / "analysis").rename(base / "analysis-detached")
                (root / "analysis").mkdir()
                root.chmod(0o555)
                with self.assertRaisesRegex(ValueError, "analysis directory changed"):
                    write_result_exclusive(
                        artifact_tree, "analysis/analysis.json", {"status": "forged"}
                    )

    def test_pinned_tree_rejects_open_analysis_directory_mode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._publish_inventory(root)
            (root / "analysis").chmod(0o755)
            with self.assertRaisesRegex(ValueError, "analysis directory changed"):
                with StructuredArtifactTree(root):
                    pass

    def test_linkage_matcher_rejects_library_prefix_lookalike(self) -> None:
        require_linked_library(
            "  libamdhip64.so.6 => /opt/rocm/lib/libamdhip64.so.6 (0x1)\n",
            "libamdhip64",
        )
        with self.assertRaises(ValueError):
            require_linked_library(
                "  libamdhip64evil.so => /tmp/libamdhip64evil.so (0x1)\n",
                "libamdhip64",
            )

    def test_gyro_requires_exact_registered_cycle_two_artifact(self) -> None:
        expected = "output/pvtk/f1_gpu_relativistic_gyro.prtcl_all.00002.part.vtk"
        self.assertEqual(registered_particle_output({expected: {}}), expected)
        with self.assertRaisesRegex(ValueError, "differs"):
            registered_particle_output(
                {
                    expected: {},
                    "output/pvtk/f1_gpu_relativistic_gyro.prtcl_all.zzzzz.part.vtk": {},
                }
            )

    def test_coupling_requires_exact_registered_output_set(self) -> None:
        expected = {}
        for label in ("coeff0", "coeff7"):
            basename = "f1_gpu_paper_coupling_" + label
            output_dir = f"output/{label}"
            for cycle in range(3):
                expected[f"{output_dir}/pvtk/{basename}.prtcl_all.{cycle:05d}.part.vtk"] = {}
            for file_id in ("mhd_bcc", "mhd_u_e", "mhd_u_m1", "mhd_u_m2", "mhd_u_m3"):
                for cycle in range(3):
                    expected[f"{output_dir}/bin/{basename}.{file_id}.{cycle:05d}.bin"] = {}
        registered_case_output_paths(expected, "coeff0")
        expected["output/coeff0/bin/unregistered.bin"] = {}
        with self.assertRaisesRegex(ValueError, "differs"):
            registered_case_output_paths(expected, "coeff0")
        del expected["output/coeff0/bin/unregistered.bin"]
        for path in ("output/coeff8/unregistered.bin", "output/unregistered.bin"):
            with self.subTest(path=path):
                expected[path] = {}
                with self.assertRaisesRegex(ValueError, "differs"):
                    registered_case_output_paths(expected, "coeff0")
                del expected[path]

    def test_particle_vtk_parsers_reject_suffix(self) -> None:
        for parser in (parse_gyro_vtk, parse_coupling_vtk):
            with self.subTest(parser=parser.__module__):
                parser(self._particle_vtk(), label="particles.vtk")
                with self.assertRaisesRegex(ValueError, "suffix"):
                    parser(self._particle_vtk(suffix=b"forged"), label="particles.vtk")


if __name__ == "__main__":
    unittest.main()
