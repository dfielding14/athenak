"""Passive analysis contracts using analytic fields and synthetic on-disk histories."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from test_cgl_lf_physics_benchmark import write_binary


REPOSITORY = Path(__file__).resolve().parents[3]
ANALYZER = REPOSITORY / "scripts/analyze_cgl_lf_physics_benchmark.py"
INPUT = REPOSITORY / "inputs/cgl_lf_paper/cgl_lf_physics_benchmark_beta10.athinput"


@pytest.fixture(scope="module")
def analyzer():
    name = "cgl_lf_physics_benchmark_passive_regression"
    spec = importlib.util.spec_from_file_location(name, ANALYZER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def passive_header(analyzer):
    header = analyzer.parse_input(INPUT.read_text().splitlines())
    header["mhd"].update(passive="true", iso_sound_speed=str(np.sqrt(5.0)))
    header["problem"]["passive_delta"] = "true"
    return header


def analytic_fields():
    shape = (8, 8, 16)
    phase = 2*np.pi*(np.arange(shape[-1])+.5)/shape[-1]
    wave = np.broadcast_to(np.cos(phase)[None, None, :], shape)
    return {"dens": 1-.05*wave, "eint": np.full(shape, 5.),
            "p_perp": 5+.25*wave, "velx": .5*wave,
            "vely": np.zeros(shape), "velz": np.zeros(shape),
            "bcc1": np.zeros(shape), "bcc2": np.zeros(shape),
            "bcc3": np.sqrt(2*(1+.25*wave))}


def test_passive_model_requires_actual_sound_speed_and_preserves_active_defaults(analyzer):
    active = analyzer.parse_input(INPUT.read_text().splitlines())
    assert "iso_sound_speed" not in analyzer.infer_model(active, {})
    header = passive_header(analyzer)
    model = analyzer.infer_model(header, {})
    assert model["iso_sound_speed"] == pytest.approx(np.sqrt(5.))
    with pytest.raises(ValueError, match="model.iso_sound_speed"):
        analyzer.infer_model(header, {"model": {"iso_sound_speed": 1.}})
    del header["mhd"]["iso_sound_speed"]
    with pytest.raises(ValueError, match="requires retained mhd/iso_sound_speed"):
        analyzer.infer_model(header, {})
    for value in ("nan", "inf", "0", "-1"):
        header["mhd"]["iso_sound_speed"] = value
        with pytest.raises(ValueError, match="finite and positive"):
            analyzer.infer_model(header, {})


def test_dynamic_pressure_and_mach_are_independent_of_passive_thermal_state(analyzer):
    model = analyzer.infer_model(passive_header(analyzer), {})
    fields = analytic_fields()
    scalars, spectra, _ = analyzer.snapshot_products(fields, (1, 1, 1), model, .1)
    assert scalars["Mach_isothermal"] == pytest.approx(np.sqrt(.125/5))
    # p_dyn = 5-.25*cos(kx) exactly compensates p_B=1+.25*cos(kx),
    # while passive p_perp=5+.25*cos(kx) correlates positively with p_B.
    assert scalars["dynamic_pressure_correlation"] == pytest.approx(-1.)
    assert scalars["dynamic_pressure_normalized_residual_variance"] < 1.e-26
    assert scalars["pressure_correlation"] == pytest.approx(1.)
    assert spectra["dynamic_isothermal_pressure"]["spectral_integral"] == pytest.approx(.25**2/2)
    for key in ("p_parallel", "p_perp", "magnetic_pressure"):
        assert key in spectra
    hotter = dict(fields, eint=4*fields["eint"], p_perp=4*fields["p_perp"])
    hot_scalars, hot_spectra, _ = analyzer.snapshot_products(hotter, (1, 1, 1), model, .1)
    assert hot_scalars["Mach_isothermal"] == scalars["Mach_isothermal"]
    assert hot_scalars["Mach_isotropic_proxy"] == pytest.approx(.5*scalars["Mach_isotropic_proxy"])
    np.testing.assert_array_equal(hot_spectra["dynamic_isothermal_pressure"]["power"],
                                  spectra["dynamic_isothermal_pressure"]["power"])


def test_passive_restart_encoding_initialization_does_not_mask_physics_changes(analyzer):
    header = passive_header(analyzer)
    initial = analyzer.comparable_model_sections(header)
    header["mhd"]["passive_restart_encoding"] = "1"
    assert analyzer.comparable_model_sections(header) == initial
    header["mhd"]["iso_sound_speed"] = "3"
    assert analyzer.comparable_model_sections(header) != initial
    header["mhd"]["passive_restart_encoding"] = "2"
    with pytest.raises(ValueError, match="unsupported retained passive_restart_encoding"):
        analyzer.comparable_model_sections(header)


def test_passive_ledger_uses_measured_kicks_without_interpreting_invariants_as_energy(analyzer):
    time = np.array([0., 1., 3.])
    user = {"time": time, "force_work": np.array([11., 13., 17.])}
    # Even an accidentally supplied tot-E must not create an active energy budget.
    mhd = {"time": time, "cgl-J": np.array([1.e8, 2.e8, -3.e8]),
           "cgl-A": np.array([7., 8., 9.]), "thermal-U": np.array([4., 6., 9.]),
           "tot-E": np.array([100., 101., 101.])}
    model = {"passive": True, "record_injected_work": True}
    result = analyzer.history_products(user, mhd, model, .5, 2.5, .5)
    assert result["forcing_work"]["actual_applied_work"] == pytest.approx(4.)
    assert result["forcing_work"]["actual_mean_total_power"] == pytest.approx(2.)
    assert result["forcing_work"]["requested_window_covered"]
    assert not result["energy_budget"]["available"]
    assert not result["energy_budget"]["applicable"]
    assert "conserved_energy_history" not in result
    assert "residual_E_minus_work" not in result["energy_budget"]
    assert result["thermal_energy_history"]["thermal-U"] == [4., 6., 9.]
    assert analyzer.energy_coverage_reasons(result, model) == []
    assert analyzer.common_energy_history(result) is None
    partial = analyzer.history_products(user, mhd, model, 0., 4., 1.)
    assert "covers only [0.0, 3.0]" in analyzer.energy_coverage_reasons(partial, model)[0]
    disabled = analyzer.history_products(user, mhd, dict(model, record_injected_work=False), 0., 3., 1.)
    assert "unavailable" in analyzer.energy_coverage_reasons(disabled, model)[0]


def test_passive_cli_exports_distinct_dynamics_and_thermal_diagnostics(analyzer, tmp_path):
    run = tmp_path / "passive_synthetic_not_solver_evidence"
    (run / "bin").mkdir(parents=True)
    header = passive_header(analyzer)
    header["job"]["basename"] = "synthetic"
    fields = analytic_fields()
    for axis, count in enumerate(fields["dens"].shape[::-1], 1):
        header["mesh"][f"nx{axis}"] = str(count)
        header["mesh"][f"x{axis}max"] = "1.0"
        header["meshblock"][f"nx{axis}"] = str(count)
    force = {"force1": fields["velx"], "force2": fields["velx"], "force3": fields["velz"]}
    for index, time in enumerate((0., 1., 2.)):
        if index:
            header["mhd"]["passive_restart_encoding"] = "1"
        write_binary(run / "bin" / f"synthetic.mhd_w_bcc.{index:05d}.bin", time, fields, header)
        write_binary(run / "bin" / f"synthetic.turb_force.{index:05d}.bin", time, force, header)
    (run / "synthetic.user.hst").write_text(
        "# Athena++ history data\n"
        "# [0]=time [1]=volume [2]=kinetic [3]=magnetic [4]=therm_cgl [5]=force_work\n"
        "0 1 .0625 1 7.5 0\n1 1 .0625 1 7.5 .32\n2 1 .0625 1 7.5 .64\n")
    (run / "synthetic.mhd.hst").write_text(
        "# Athena++ history data\n# [0]=time [1]=cgl-J [2]=cgl-A [3]=thermal-U\n"
        "0 99 44 7.5\n1 -77 101 7.5\n2 22 13 7.5\n")
    (run / "benchmark_metadata.json").write_text(json.dumps({
        "purpose": "Synthetic analyzer regression; no solver was run", "simulation": {}, "outputs": {}}))
    output = run / "analysis"
    proc = subprocess.run([sys.executable, str(ANALYZER), str(run), "--time-start", "0",
        "--time-end", "2", "--block-duration", ".5", "--output-dir", str(output)],
        cwd=run, capture_output=True, text=True, timeout=120)
    (run / "analysis.log").write_text(proc.stdout+proc.stderr)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    result = json.loads((output / "metrics.json").read_text())
    assert result["history"]["forcing_work"]["requested_window_covered"]
    assert not result["history"]["energy_budget"]["applicable"]
    assert not any("forcing/conserved-energy" in reason for reason in result["adequacy"]["reasons"])
    assert "Mach_isothermal" in result["scalars"]
    assert "dynamic_isothermal_pressure" in result["spectra"]
    assert result["scalars"]["dynamic_pressure_correlation"]["mean"] == pytest.approx(-1., abs=1.e-10)
    report = (output / "report.md").read_text()
    assert "single-run physics benchmark supplement" in report
    assert "cgl-J and cgl-A are invariant integrals, not energies" in report
    audit = json.loads((output / "figure-audit.json").read_text())
    assert all(record["ok"] for record in audit["figures"].values())
    for group in ("marginality", "pressure_balance", "gradients", "spectra_energy"):
        assert (output / f"{group}.png").stat().st_size > 100
        assert (output / f"{group}.pdf").stat().st_size > 100
