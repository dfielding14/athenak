#!/usr/bin/env python3
"""Regression tests for the structured Frontier F2 runtime-metadata analyzer."""

from __future__ import annotations

import base64
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest

PUBLICATION_DIR = Path(__file__).resolve().parent
CONTROL_PLANE_DIR = PUBLICATION_DIR / "frontier_control_plane"
sys.path.insert(0, str(PUBLICATION_DIR))
sys.path.insert(0, str(CONTROL_PLANE_DIR))

from control_plane_common import PRODUCTION_RUNTIME_LOADED_MODULES
from control_plane_common import PRODUCTION_RUNTIME_MODULEFILES
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH
from frontier_f1_structured_artifacts import TRUSTED_PYTHON
from frontier_f2_multirank_runtime_metadata_analysis import EXPECTED_BASELINE
from frontier_f2_multirank_runtime_metadata_analysis import RUNTIME_ALLOWLIST_KEYS
from frontier_f2_multirank_runtime_metadata_analysis import analyze


class FrontierF2MultirankRuntimeMetadataAnalysisTests(unittest.TestCase):
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
        inventory_path = root / "artifact_inventory.json"
        inventory_path.write_text(
            json.dumps({"schema_version": 1, "files": records}, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        inventory_path.chmod(0o444)
        root.chmod(0o555)

    def _artifacts(self, root: Path) -> None:
        values = {
            **{key: "<unset>" for key in RUNTIME_ALLOWLIST_KEYS},
            **EXPECTED_BASELINE,
            "LOADEDMODULES": ":".join(PRODUCTION_RUNTIME_LOADED_MODULES),
            "_LMFILES_": ":".join(PRODUCTION_RUNTIME_MODULEFILES),
            "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH,
        }
        allowlist = root / "f2-runtime-metadata.environment.allowlist.txt"
        allowlist.write_text(
            "".join(f"{key}={values[key]}\n" for key in RUNTIME_ALLOWLIST_KEYS),
            encoding="utf-8",
        )
        preflight = "".join(
            "PIC trusted GPU launch: "
            f"rank={rank} host=frontier00001 ROCR_VISIBLE_DEVICES={rank} "
            "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
            for rank in range(8)
        )
        runtime = (
            "Number of parallel ranks = 8\n"
            "PIC runtime model: physical_mode=extended_mhd_pic state=momentum_p_over_m "
            "C=3 background=coupled feedback=coupled induction=ideal_mhd_only "
            "deposition=tsc deltaf=off deltaf_adapt=off deltaf_adapt_interval=0 "
            "expanding_box=off expansion_law=linear wave_damping=off nu_in=0 "
            "lb_cost_per_particle=0 max_cell_cross=2 theta_max=0.3 restart_schema=7\n"
        )
        (root / "athena_stdout.txt").write_text(preflight + runtime, encoding="utf-8")
        diagnostic = base64.b64decode(
            (
                PUBLICATION_DIR / "readiness/frontier_mpich_diagnostic_stderr.base64"
            ).read_bytes().replace(b"\n", b""),
            validate=True,
        )
        (root / "athena_stderr.txt").write_bytes(diagnostic)

    def test_multirank_runtime_metadata_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            self._publish_inventory(root)
            result = analyze(root)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["parallel_ranks"], 8)
            self.assertEqual(result["hosts"], ["frontier00001"])
            self.assertEqual(
                [binding["rocr_visible_device"] for binding in result["rank_gpu_bindings"]],
                list(range(8)),
            )

    def test_multirank_runtime_metadata_rejects_rank_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            stdout = root / "athena_stdout.txt"
            stdout.write_text(
                stdout.read_text(encoding="utf-8").replace("rank=7 ", "rank=6 ", 1),
                encoding="utf-8",
            )
            self._publish_inventory(root)
            with self.assertRaisesRegex(ValueError, "exactly ranks 0 through 7"):
                analyze(root)

    def test_multirank_runtime_metadata_rejects_gpu_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            stdout = root / "athena_stdout.txt"
            stdout.write_text(
                stdout.read_text(encoding="utf-8").replace(
                    "ROCR_VISIBLE_DEVICES=7", "ROCR_VISIBLE_DEVICES=0", 1
                ),
                encoding="utf-8",
            )
            self._publish_inventory(root)
            with self.assertRaisesRegex(ValueError, "exactly devices 0 through 7"):
                analyze(root)

    def test_multirank_runtime_metadata_rejects_multinode_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            stdout = root / "athena_stdout.txt"
            stdout.write_text(
                stdout.read_text(encoding="utf-8").replace(
                    "rank=7 host=frontier00001", "rank=7 host=frontier00002", 1
                ),
                encoding="utf-8",
            )
            self._publish_inventory(root)
            with self.assertRaisesRegex(ValueError, "one Frontier node"):
                analyze(root)

    def test_multirank_runtime_metadata_rejects_runtime_model_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            stdout = root / "athena_stdout.txt"
            stdout.write_text(
                stdout.read_text(encoding="utf-8").replace(
                    "background=coupled", "background=passive_mhd", 1
                ),
                encoding="utf-8",
            )
            self._publish_inventory(root)
            with self.assertRaisesRegex(ValueError, "missing expected tokens"):
                analyze(root)

    def test_trusted_runner_publishes_receipt_and_recomputes_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "artifacts"
            root.mkdir()
            self._artifacts(root)
            self._publish_inventory(root)
            snapshot = base / "snapshot"
            snapshot.mkdir()
            analyzer = snapshot / "000-frontier_f2_multirank_runtime_metadata_analysis.py"
            helper = snapshot / "frontier_f1_structured_artifacts.py"
            shutil.copyfile(Path(__file__).with_name(
                "frontier_f2_multirank_runtime_metadata_analysis.py"
            ), analyzer)
            shutil.copyfile(Path(__file__).with_name(
                "frontier_f1_structured_artifacts.py"
            ), helper)
            analyzer.chmod(0o444)
            helper.chmod(0o444)
            snapshot.chmod(0o555)
            subprocess.run(
                [TRUSTED_PYTHON, "-I", "-B", str(analyzer), "--artifact-dir", str(root)],
                check=True,
            )
            inventory = root / "artifact_inventory.json"
            result = root / "analysis" / "analysis.json"
            receipt = root / "analysis" / "offline_analysis_receipt.json"
            for path in (result, receipt):
                self.assertTrue(path.is_file())
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o444)
            descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            try:
                subprocess.run(
                    [
                        TRUSTED_PYTHON,
                        "-I",
                        "-B",
                        str(analyzer),
                        "--artifact-dir",
                        str(root),
                        "--artifact-dir-fd",
                        str(descriptor),
                        "--verify-artifact-inventory-sha256",
                        hashlib.sha256(inventory.read_bytes()).hexdigest(),
                        "--verify-result-sha256",
                        hashlib.sha256(result.read_bytes()).hexdigest(),
                    ],
                    check=True,
                    pass_fds=(descriptor,),
                )
            finally:
                os.close(descriptor)


if __name__ == "__main__":
    unittest.main()
