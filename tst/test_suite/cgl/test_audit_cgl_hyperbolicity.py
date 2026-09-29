"""Focused regressions for the external retained-state hyperbolicity audit."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest


AUDIT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/"
    "analysis/audit_cgl_hyperbolicity.py"
)


def load_audit():
    name = "audit_cgl_hyperbolicity_tests"
    spec = importlib.util.spec_from_file_location(name, AUDIT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def audit():
    return load_audit()


def bhoriya_unexpanded_discriminant(
    density: float,
    p_parallel: float,
    p_perp: float,
    b_squared: float,
    b_parallel: float,
) -> float:
    """Return b^2 - 4ac from Bhoriya et al.'s unexpanded coefficients."""
    mu_squared = b_parallel * b_parallel / b_squared
    q_squared = (
        b_squared
        + 2.0 * p_perp
        + (2.0 * p_parallel - p_perp) * mu_squared
    )
    a = 2.0 * density
    b = -q_squared
    c = -(
        3.0 * mu_squared * mu_squared * p_parallel * p_parallel
        - mu_squared * p_perp * p_perp * (mu_squared - 1.0)
        - 3.0 * b_parallel * b_parallel * p_parallel
        + 3.0
        * mu_squared
        * p_parallel
        * p_perp
        * (mu_squared - 2.0)
    ) / (2.0 * density)
    return b * b - 4.0 * a * c


@pytest.mark.parametrize(
    (
        "density",
        "p_parallel",
        "p_perp",
        "b_squared",
        "b_parallel",
        "expected",
    ),
    [
        pytest.param(2.5, 3.0, 0.75, 4.0, 0.0, 30.25, id="perpendicular"),
        pytest.param(
            0.4,
            0.5,
            0.5,
            2.0,
            math.sqrt(2.0),
            0.25,
            id="isotropic-parallel",
        ),
        pytest.param(
            7.0,
            1.0,
            0.5,
            1.0,
            math.sqrt(0.65),
            1.083125,
            id="oblique-sign-changing",
        ),
    ],
)
def test_literature_correct_formula_matches_unexpanded_bhoriya_oracle(
    audit,
    density: float,
    p_parallel: float,
    p_perp: float,
    b_squared: float,
    b_parallel: float,
    expected: float,
) -> None:
    oracle = bhoriya_unexpanded_discriminant(
        density, p_parallel, p_perp, b_squared, b_parallel
    )
    _, observed = audit.discriminant(
        np.asarray([p_parallel]),
        np.asarray([p_perp]),
        np.asarray([b_squared]),
        np.asarray([b_parallel]),
        "literature-correct",
    )

    assert oracle == pytest.approx(expected)
    assert observed.item() == pytest.approx(oracle)


def test_literature_correct_formula_is_density_independent(audit) -> None:
    state = {
        "p_parallel": 1.3,
        "p_perp": 0.8,
        "b_squared": 2.2,
        "b_parallel": math.sqrt(0.9),
    }
    oracle_values = [
        bhoriya_unexpanded_discriminant(density=density, **state)
        for density in (0.125, 1.0, 64.0)
    ]
    _, observed = audit.discriminant(
        np.asarray([state["p_parallel"]]),
        np.asarray([state["p_perp"]]),
        np.asarray([state["b_squared"]]),
        np.asarray([state["b_parallel"]]),
        "literature-correct",
    )

    assert oracle_values[0] == pytest.approx(oracle_values[1])
    assert oracle_values[1] == pytest.approx(oracle_values[2])
    assert observed.item() == pytest.approx(oracle_values[0])


def test_legacy_and_literature_correct_formulas_diverge_on_known_fixture(audit) -> None:
    p_parallel = np.asarray([1.0])
    p_perp = np.asarray([0.5])
    b_squared = np.asarray([1.0])
    mu_squared = 0.65
    b_parallel = np.asarray([math.sqrt(mu_squared)])

    _, legacy = audit.discriminant(
        p_parallel, p_perp, b_squared, b_parallel, "qualified-legacy"
    )
    _, correct = audit.discriminant(
        p_parallel, p_perp, b_squared, b_parallel, "literature-correct"
    )
    expected_correction = (
        12.0 * mu_squared * mu_squared * 1.0 * (1.0 - 0.5)
    )

    assert legacy.item() == pytest.approx(-1.451875)
    assert correct.item() == pytest.approx(1.083125)
    assert correct.item() - legacy.item() == pytest.approx(expected_correction)
    assert legacy.item() < 0.0 < correct.item()


def test_expand_inputs_rejects_elf_bin_and_keeps_athena_snapshot(
    audit, tmp_path: Path
) -> None:
    elf = tmp_path / "athena-executable.bin"
    elf.write_bytes(b"\x7fELF\x02\x01\x01\x00not-an-athena-snapshot")
    snapshot = tmp_path / "retained-state.bin"
    snapshot.write_bytes(b"Athena binary output version=1.1\n")

    assert not audit.is_athena_binary(elf)
    assert audit.is_athena_binary(snapshot)
    assert audit.expand_inputs([str(tmp_path)]) == [snapshot.resolve()]

    with pytest.raises(audit.AuditError, match="matched no .bin files"):
        audit.expand_inputs([str(elf)])


@pytest.mark.parametrize(
    ("formula_id", "reference_fragment"),
    [
        pytest.param(
            "qualified-legacy",
            "qualified AthenaK commit 9e07542281e4e6d125582f253df3ad2e3b8b154d",
            id="legacy",
        ),
        pytest.param(
            "literature-correct",
            "Bhoriya et al. (2024), eigenvalue coefficient c",
            id="literature-correct",
        ),
    ],
)
def test_json_report_records_selected_formula_provenance(
    audit,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    formula_id: str,
    reference_fragment: str,
) -> None:
    helper = tmp_path / "bin_convert.py"
    helper.write_text("# fixture helper\n", encoding="utf-8")
    snapshot = tmp_path / "snapshot.bin"
    snapshot.write_bytes(b"Athena binary output version=1.1\n")
    args = SimpleNamespace(
        bin_convert=helper,
        inputs=[str(snapshot)],
        formula=formula_id,
        format="json",
        hash_inputs=False,
        bfloor=None,
        fail_on_negative=False,
    )
    monkeypatch.setattr(audit, "parse_args", lambda: args)
    monkeypatch.setattr(audit, "expand_inputs", lambda _inputs: [snapshot.resolve()])
    monkeypatch.setattr(audit, "import_bin_convert", lambda _path: object())
    monkeypatch.setattr(
        audit,
        "git_provenance",
        lambda _path: {"root": "/fixture/repository", "revision": "abc123"},
    )
    monkeypatch.setattr(
        audit,
        "scan_snapshot",
        lambda _key, _paths, _helper, _args: {"aggregate": {"negative": 0}},
    )

    assert audit.main() == 0
    report = json.loads(capsys.readouterr().out)
    provenance = report["provenance"]

    assert audit.FORMULA_IDS == ("qualified-legacy", "literature-correct")
    assert provenance["formula_id"] == formula_id
    assert reference_fragment in provenance["formula_reference"]
    assert provenance["script_path"] == str(AUDIT.resolve())
    assert provenance["script_version"] == audit.SCRIPT_VERSION
    assert provenance["writes_files"] is False
