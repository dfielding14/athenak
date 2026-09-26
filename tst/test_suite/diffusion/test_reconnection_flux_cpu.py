"""Sign, orientation, and normalization input for the flux quadrature."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "vis/python"))
from reconnection_flux import sheet_flux, sheet_geometry  # noqa: E402


def test_flux_integral_and_reversed_points():
    x = np.linspace(-1, 1, 101)
    # Az = x^2, hence By = -2x. Linear By has exact trapezoidal quadrature.
    assert sheet_flux(x, -2*x, 0, 0.6) == pytest.approx(0.36)
    assert sheet_flux(x, -2*x, 0.6, 0) == pytest.approx(-0.36)
    assert sheet_flux(x, -2*x, 0, 1.01) == pytest.approx(1.01**2)
    with pytest.raises(ValueError, match="within"):
        sheet_flux(x, -2*x, 0, 2)


def test_layer_geometry_and_unbounded_length():
    x, y = np.linspace(-1, 1, 401), np.linspace(-0.2, 0.2, 401)
    bx = np.tanh(y[:, None]/0.05) * np.exp(-x[None, :]**2/0.3**2)
    by, rho = np.zeros_like(bx), np.full_like(bx, 4.0)
    delta, normalized, length, aspect = sheet_geometry(x, y, bx, by, rho, 0, 0, 1, 0.1)
    assert delta == pytest.approx(0.05, rel=0.001)
    assert normalized == pytest.approx(1, rel=0.001)
    assert length == pytest.approx(0.3*np.sqrt(np.log(2)), rel=0.001)
    assert aspect == pytest.approx(delta/length)
    bx[:] = np.tanh(y[:, None]/0.05)
    assert np.isnan(sheet_geometry(x, y, bx, by, rho, 0, 0, 1, 0.1)[2])
