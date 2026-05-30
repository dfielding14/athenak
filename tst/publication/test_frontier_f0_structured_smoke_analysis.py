#!/usr/bin/env python3
"""Regression tests for the structured Frontier F0 offline analyzer."""

from __future__ import annotations

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
from frontier_f0_structured_smoke_analysis import EXPECTED_BASELINE
from frontier_f0_structured_smoke_analysis import RUNTIME_ALLOWLIST_KEYS
from frontier_f0_structured_smoke_analysis import _write_result, analyze


class FrontierF0StructuredSmokeAnalysisTests(unittest.TestCase):
    def _artifacts(self, root: Path, *, profile: str = "frontier_minimum_supported") -> None:
        values = {
            **{key: "<unset>" for key in RUNTIME_ALLOWLIST_KEYS},
            **EXPECTED_BASELINE,
            "PIC_FRONTIER_PROFILE": profile,
            "LOADEDMODULES": ":".join(PRODUCTION_RUNTIME_LOADED_MODULES),
            "_LMFILES_": ":".join(PRODUCTION_RUNTIME_MODULEFILES),
            "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH,
        }
        allowlist = root / "f0-parser.environment.allowlist.txt"
        allowlist.write_text(
            "".join(f"{key}={values[key]}\n" for key in RUNTIME_ALLOWLIST_KEYS),
            encoding="utf-8",
        )
        allowlist.chmod(0o400)
        (root / "athena_stdout.txt").write_text("AthenaK parser smoke\n", encoding="utf-8")
        (root / "athena_stderr.txt").write_text("", encoding="utf-8")

    def test_structured_smoke_artifacts_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            result = analyze(root)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["runtime_profile"], "frontier_minimum_supported")
            self.assertEqual(set(result["runtime_artifacts"]), {
                "f0-parser.environment.allowlist.txt",
                "athena_stdout.txt",
                "athena_stderr.txt",
            })

    def test_structured_smoke_rejects_profile_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root, profile="frontier_xnack1_experimental")
            with self.assertRaisesRegex(ValueError, "PIC_FRONTIER_PROFILE"):
                analyze(root)

    def test_structured_smoke_rejects_writable_allowlist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            (root / "f0-parser.environment.allowlist.txt").chmod(0o600)
            with self.assertRaisesRegex(ValueError, "read-only"):
                analyze(root)

    def test_structured_smoke_rejects_module_provenance_drift(self) -> None:
        variants = {
            "loaded-modules": {
                "LOADEDMODULES": ":".join(PRODUCTION_RUNTIME_LOADED_MODULES[1:]),
            },
            "modulefiles": {
                "_LMFILES_": ":".join(PRODUCTION_RUNTIME_MODULEFILES) + ":/tmp/forged.lua",
            },
            "modulepath-prefix": {
                "MODULEPATH": "/tmp/forged:" + PRODUCTION_RUNTIME_MODULEPATH,
            },
            "modulepath-suffix": {
                "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH + ":/tmp/forged",
            },
            "modulepath-missing-component": {
                "MODULEPATH": PRODUCTION_RUNTIME_MODULEPATH.split(":", 1)[1],
            },
        }
        for label, replacements in variants.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self._artifacts(root)
                allowlist = root / "f0-parser.environment.allowlist.txt"
                values = dict(
                    line.split("=", 1)
                    for line in allowlist.read_text(encoding="utf-8").splitlines()
                )
                values.update(replacements)
                allowlist.chmod(0o600)
                allowlist.write_text(
                    "".join(f"{key}={values[key]}\n" for key in RUNTIME_ALLOWLIST_KEYS),
                    encoding="utf-8",
                )
                allowlist.chmod(0o400)
                with self.assertRaisesRegex(ValueError, "reviewed provenance"):
                    analyze(root)

    def test_result_publication_is_read_only_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._artifacts(root)
            output = root / "analysis.json"
            _write_result(output, analyze(root))
            self.assertEqual(output.stat().st_mode & 0o777, 0o444)
            with self.assertRaises(FileExistsError):
                _write_result(output, analyze(root))


if __name__ == "__main__":
    unittest.main()
