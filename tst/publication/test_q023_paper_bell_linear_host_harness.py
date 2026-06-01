#!/usr/bin/env python3
"""Host-contract and registration tests for the Q-023 Section 5.2 generator."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
HARNESS = REPO_ROOT / "tst/publication/q023_paper_bell_linear_host_harness.cpp"
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/q023_paper_bell_linear_source_local_implementation_successor_v3_2026-06-01.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inventory_sha256(root: Path) -> str:
    rows = [
        f"{_sha256(path)}  {path.relative_to(root).as_posix()}\n"
        for path in sorted(path for path in root.rglob("*") if path.is_file())
    ]
    return hashlib.sha256("".join(rows).encode("utf-8")).hexdigest()


def _basis(dimension: int) -> tuple[list[float], list[float], list[float]]:
    raw = [1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0]
    norm = math.sqrt(sum(value*value for value in raw))
    parallel = [value / norm for value in raw]
    if dimension == 1:
        transverse_a = [0.0, 1.0, 0.0]
    else:
        transverse_a = [-parallel[1], parallel[0], 0.0]
        transverse_norm = math.sqrt(sum(value*value for value in transverse_a))
        transverse_a = [value / transverse_norm for value in transverse_a]
    transverse_b = [
        parallel[1]*transverse_a[2] - parallel[2]*transverse_a[1],
        parallel[2]*transverse_a[0] - parallel[0]*transverse_a[2],
        parallel[0]*transverse_a[1] - parallel[1]*transverse_a[0],
    ]
    return parallel, transverse_a, transverse_b


class Q023PaperBellLinearHostHarnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary_directory = tempfile.TemporaryDirectory()
        binary = Path(cls._temporary_directory.name) / "q023-paper-bell-linear-host"
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(HARNESS),
                "-o",
                str(binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        cls.lines = subprocess.run(
            [str(binary)],
            cwd=REPO_ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary_directory.cleanup()

    def test_host_eigenmode_matches_section52_phase_relation(self) -> None:
        sample_lines = [line.split() for line in self.lines if line.startswith("sample ")]
        self.assertEqual(len(sample_lines), 27)
        for fields in sample_lines:
            dimension = int(fields[1])
            epsilon = float(fields[2])
            phase = float(fields[3])
            measured_b = [float(value) for value in fields[4:7]]
            measured_u = [float(value) for value in fields[7:10]]
            parallel, transverse_a, transverse_b = _basis(dimension)
            growth = math.sqrt(1.0 - epsilon*epsilon)
            expected_b = [
                parallel[axis]
                + 1.0e-6*(
                    math.cos(phase)*transverse_a[axis]
                    - math.sin(phase)*transverse_b[axis]
                )
                for axis in range(3)
            ]
            expected_u = [
                1.0e-6*(
                    (-epsilon*math.cos(phase) + growth*math.sin(phase))
                    * transverse_a[axis]
                    + (growth*math.cos(phase) + epsilon*math.sin(phase))
                    * transverse_b[axis]
                )
                for axis in range(3)
            ]
            for measured, expected in zip(measured_b, expected_b):
                self.assertAlmostEqual(measured, expected, places=14)
            for measured, expected in zip(measured_u, expected_u):
                self.assertAlmostEqual(measured, expected, places=14)

    def test_vector_potential_curl_matches_eigenmode(self) -> None:
        curl_lines = [line.split() for line in self.lines if line.startswith("curl ")]
        self.assertEqual(len(curl_lines), 9)
        for fields in curl_lines:
            for residual in fields[3:]:
                self.assertLess(abs(float(residual)), 3.0e-12)

    def test_builtin_registration_covers_fresh_and_restart_paths(self) -> None:
        pgen = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
        self.assertEqual(pgen.count('compare("q023_paper_bell_linear")'), 2)
        self.assertEqual(pgen.count("Q023PaperBellLinear(pin, false);"), 1)
        self.assertEqual(pgen.count("Q023PaperBellLinear(pin, true);"), 1)
        self.assertIn(
            "void Q023PaperBellLinear(ParameterInput *pin, const bool restart);",
            (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "pgen/tests/q023_paper_bell_linear.cpp",
            (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8"),
        )

    def test_source_local_sidecar_hashes_and_limits(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["campaign_id"], "Q023-PAPER-BELL-LINEAR")
        self.assertEqual(sidecar["qualification_effect"], "none")
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        self.assertEqual(sidecar["section52_qualification"], "not_claimed")
        self.assertEqual(
            sidecar["physical_1d_carrier_contract"]["semantics"],
            "transverse_invariant_thin_2d3v_carrier",
        )
        self.assertEqual(
            sidecar["physical_1d_carrier_contract"]["mesh"],
            {"nx1": 32, "nx2": 4, "nx3": 1},
        )
        smoke = sidecar["source_local_runtime_smoke"]
        self.assertEqual(
            smoke["qualification_effect"],
            "bounded_source_local_preparation_smoke_only",
        )
        self.assertEqual(smoke["retention"]["writable_entries"], 0)
        self.assertEqual(
            [item["carrier_shape"] for item in smoke["cycle_zero_initializations"]],
            [[1, 4, 32], [1, 32, 64], [32, 64, 128]],
        )
        self.assertTrue(
            all(item["result"] == "pass"
                for item in smoke["cycle_zero_initializations"])
        )
        self.assertEqual(smoke["restart_continuation"]["result"], "pass")
        self.assertEqual(smoke["restart_continuation"]["final_cycle"], 1)
        artifact_root = Path(smoke["artifact_root"])
        self.assertEqual(
            _inventory_sha256(artifact_root),
            smoke["retention"]["inventory_sha256"],
        )
        executable = smoke["debug_executable"]
        self.assertEqual(_sha256(Path(executable["path"])), executable["sha256"])
        for item in smoke["cycle_zero_initializations"]:
            self.assertEqual(
                _sha256(artifact_root / item["selected_raw_mhd_w_bcc_path"]),
                item["selected_raw_mhd_w_bcc_sha256"],
            )
            self.assertLess(
                item["source_local_velocity_magnetic_ratio_absolute_error"],
                smoke["source_local_velocity_magnetic_ratio_diagnostic_limit"],
            )
        restart = smoke["restart_continuation"]
        self.assertEqual(
            _sha256(artifact_root / restart["loaded_restart_path"]),
            restart["loaded_restart_sha256"],
        )
        self.assertEqual(
            _sha256(artifact_root / restart["selected_continuation_raw_mhd_w_bcc_path"]),
            restart["selected_continuation_raw_mhd_w_bcc_sha256"],
        )
        for artifact in sidecar["source_local_artifacts"]:
            self.assertEqual(_sha256(REPO_ROOT / artifact["path"]), artifact["sha256"])
        self.assertTrue(
            any("clean Frontier executable" in item
                for item in sidecar["remaining_open_dependencies"])
        )


if __name__ == "__main__":
    unittest.main()
