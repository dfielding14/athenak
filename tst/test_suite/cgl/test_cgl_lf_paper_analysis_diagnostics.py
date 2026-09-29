from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[3]
ANALYZER = REPOSITORY / "scripts/analyze_cgl_lf_paper.py"


def load_analyzer():
    name = "cgl_lf_paper_analysis_diagnostics"
    spec = importlib.util.spec_from_file_location(name, ANALYZER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_pressure_balance_detects_exact_compensation():
    analyzer = load_analyzer()
    coordinate = np.linspace(0.0, 2.0 * math.pi, 64, endpoint=False)
    fluctuation = np.sin(coordinate)[None, None, :]
    record = analyzer.pressure_balance_diagnostics(
        2.0 + fluctuation,
        3.0 - fluctuation,
    )
    assert record["available"] is True
    assert math.isclose(record["correlation"], -1.0, abs_tol=1.0e-14)
    assert math.isclose(
        record["perpendicular_on_magnetic_regression_slope"],
        -1.0,
        abs_tol=1.0e-14,
    )
    assert math.isclose(
        record["normalized_residual_variance"], 0.0, abs_tol=1.0e-14
    )


def test_compressive_velocity_power_fraction_distinguishes_modes():
    analyzer = load_analyzer()
    shape = (8, 8, 8)
    lengths = (1.0, 1.0, 1.0)
    x = (np.arange(shape[2]) + 0.5) / shape[2]
    wave = np.broadcast_to(
        np.sin(2.0 * math.pi * x)[None, None, :], shape
    )
    zero = np.zeros(shape)
    dk = 2.0 * math.pi

    longitudinal = [
        wave,
        zero,
        zero,
    ]
    transverse = [
        zero,
        wave,
        zero,
    ]
    for velocity, expected in ((longitudinal, 1.0), (transverse, 0.0)):
        total = analyzer.shell_spectrum(velocity, lengths, dk)
        compressive = analyzer.compressive_velocity_spectrum(
            velocity, lengths, dk
        )
        result = analyzer.compressive_velocity_power_fraction(
            total, compressive
        )
        assert result["available"] is True
        assert math.isclose(result["fraction"], expected, abs_tol=1.0e-14)


def test_scalar_diagnostic_ensemble_preserves_descriptive_uncertainty():
    analyzer = load_analyzer()
    records = [
        {
            "available": True,
            "definition": "test",
            "scope": "descriptive",
            "fraction": 0.1,
        },
        {
            "available": True,
            "definition": "test",
            "scope": "descriptive",
            "fraction": 0.3,
        },
    ]
    result = analyzer.mean_scalar_diagnostics(records, [4.0, 10.0])
    assert result["available"] is True
    assert math.isclose(result["fraction_mean"], 0.2)
    assert result["uncertainty"]["fraction"]["available"] is True
