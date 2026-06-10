"""Focused tests for the R14/R15 C2P and heat-flux stability launcher."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys


REPOSITORY = Path(__file__).resolve().parents[3]
LAUNCHER = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_fast_c2p_stable.py"


def load_launcher():
    name = "cgl_lf_stage_i_fast_c2p_stable_for_tests"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_preflight_default_is_restricted_to_r14_r15():
    launcher = load_launcher()
    calls: list[tuple[Path, tuple[str, ...]]] = []

    def validate(root: Path, cases: tuple[str, ...]) -> dict[str, object]:
        calls.append((root, cases))
        return {"cases": list(cases)}

    corrected = SimpleNamespace(validate_provenance=validate)
    launcher.restrict_validate_provenance(corrected)

    root = Path("/tmp/c2p-stable-campaign")
    assert corrected.validate_provenance(root) == {"cases": ["R14", "R15"]}
    assert calls == [(root, ("R14", "R15"))]
