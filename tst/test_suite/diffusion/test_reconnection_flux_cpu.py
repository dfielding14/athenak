"""Sign, orientation, and normalization input for the flux quadrature."""

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "vis/python"))
from reconnection_flux import sheet_flux, sheet_geometry  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "benchmarks/reconnection"))
import analyze as campaign  # noqa: E402
import compare as comparison  # noqa: E402


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


def test_campaign_boundary_budget_and_history_cadence(tmp_path, monkeypatch):
    case = tmp_path / "case"
    (case / "bin").mkdir(parents=True)
    xf, yf = np.linspace(-0.5, 0.5, 9), np.linspace(-0.25, 0.25, 5)
    x, y = (xf[1:]+xf[:-1])/2, (yf[1:]+yf[:-1])/2
    xx, yy = np.meshgrid(x, y)
    final = 2-2e-7  # Six-significant-digit binary time rounds up to 2.
    snapshots = {}
    for number, time in enumerate((0, final)):
        path = case / "bin" / f"test.state.{number:05d}.bin"
        path.touch()
        bx, by, bz = 2*yy, (1+time)*xx, np.ones_like(xx)
        fields = dict(bcc1=bx, bcc2=by, bcc3=bz, dens=bz,
                      mom1=0.2*bz, mom2=0*bz, mom3=0*bz,
                      ener=3+0.02+(bx*bx+by*by+bz*bz)/2)
        snapshots[path] = dict(Time=float(f"{time:.6g}"), NumCycles=number,
                               x1v=x, x2v=y, x1f=xf, x2f=yf,
                               **{key: value[None] for key, value in fields.items()})
    header = ["<problem>", "b0=1", "density=1", "<mhd>", "d_i=0.1",
              "gamma=1.6666666666666667", "ohmic_resistivity=0.01",
              "<mesh>", "x3min=-1", "x3max=1"]
    raw = dict(header=header, Nx1=8, Nx2=4, Nx3=1, n_mbs=1,
               nx1_mb=8, nx2_mb=4, nx3_mb=1, nx1_out_mb=8, nx2_out_mb=4,
               nx3_out_mb=1, mb_logical=np.array([[0, 0, 0, 0]]))
    monkeypatch.setattr(campaign, "read_binary", lambda path: raw)
    monkeypatch.setattr(campaign, "read_binary_as_athdf", lambda path, **kw: snapshots[path])
    first = campaign.snapshot(next(iter(snapshots)), 1e-6)
    assert first["psi_ref_minus_X"] == pytest.approx(-0.125)
    assert first["fixed_X_verified"] and first["O_candidates_midplane"] == 0
    assert first["outward_energy_flux_approx"] == pytest.approx(-0.005)
    assert first["internal_energy"] == pytest.approx(1.5)
    assert first["total_energy_snapshot"] == pytest.approx(
        first["internal_energy"]+first["magnetic_energy_bcc"]+first["kinetic_energy"])

    user = dict(time=np.array([0, 1, final]), q_max_edge=np.array([0.2, 2, 0.2]),
                frac_qstar=np.array([0, 0.3, 0]), frac_q1=np.array([0, 0.2, 0]))
    for key, value in dict(x_q=0.1, etaJ2_cc=0.04, heat_frac=0.5, x_etaJz=-0.02,
                           x_Ez=-0.25, ref_etaJz=-0.01, ref_Ez=-0.125).items():
        user[key] = np.full(3, value)
    mhd = dict(time=user["time"], **{"tot-E": np.full(3, 4.0)})
    for kind in ("user", "mhd"):
        (case / f"test.{kind}.hst").touch()
    monkeypatch.setattr(campaign, "hst", lambda path: user if ".user." in path else mhd)
    summary = campaign.analyze(case, tmp_path, 1e-6)
    assert summary["max_q"] == 2 and summary["max_q_occupancy"] == 0.2
    assert summary["end"] == final
    with (tmp_path / "case.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[-1]["time_binary"]) == 2
    assert float(rows[-1]["rate_minus_physical_EMF"]) == pytest.approx(0, abs=1e-14)
    assert float(rows[-1]["etaJ2_cc"]) == pytest.approx(0.02)  # Per unit depth.
    (case / "run.json").write_text(json.dumps(dict(status="complete",
        physics=dict(model="constant", cells_per_di=2), arguments=dict(cfl=0.4))))
    monkeypatch.setattr(comparison, "hst", lambda path: user)
    budget = comparison.flux_budget(case, tmp_path)
    assert budget["mean_flux_rate"] == pytest.approx(-0.125)
    assert budget["mean_rate_discrepancy"] == pytest.approx(0, abs=1e-14)
    # A genuinely uncovered endpoint must not be silently clamped to history.
    snapshots[next(reversed(snapshots))]["Time"] = 2.001
    with pytest.raises(ValueError, match="does not cover"):
        campaign.analyze(case, tmp_path, 1e-6)


def test_campaign_field_restriction_preserves_affine_fields(tmp_path, monkeypatch):
    states, cases = {}, []
    for n in (4, 8):
        case = tmp_path / str(n)
        (case / "bin").mkdir(parents=True)
        path = case / "bin" / "test.state.00001.bin"
        path.touch()
        coordinate = (np.arange(n)+0.5)/n
        x, y = np.meshgrid(coordinate, coordinate)
        states[path] = dict(Time=0.2, bcc1=x[None], bcc2=y[None], bcc3=(x+2*y)[None])
        cases.append(case)
    monkeypatch.setattr(comparison, "read_binary_as_athdf", lambda path, **kw: states[path])
    difference = comparison.field_difference(*cases)
    assert max(value for key, value in difference.items() if key.startswith("B")) < 1e-15
