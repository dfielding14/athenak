"""Analytic regression checks for the retained-run physics benchmark analyzer."""

import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np
import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
ANALYZER = REPOSITORY / "scripts/analyze_cgl_lf_physics_benchmark.py"


@pytest.fixture(scope="module")
def analyzer():
    name = "cgl_lf_physics_benchmark_regression"
    spec = importlib.util.spec_from_file_location(name, ANALYZER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def uniform_fields(shape):
    return {
        "dens": np.ones(shape),
        "eint": np.full(shape, 5.0),
        "p_perp": np.full(shape, 5.0),
        "velx": np.zeros(shape),
        "vely": np.zeros(shape),
        "velz": np.zeros(shape),
        "bcc1": np.zeros(shape),
        "bcc2": np.zeros(shape),
        "bcc3": np.ones(shape),
    }


MODEL = {"B0": 1.0, "mirror_X": 1.0, "firehose_X": -2.0,
         "gamma_sound": 5.0 / 3.0}


def test_physical_shell_parseval_and_helmholtz_power(analyzer):
    shape = (8, 16, 8)
    lengths = (1.0, 2.0, 1.0)
    y = (np.arange(shape[1]) + 0.5) * lengths[1] / shape[1]
    transverse = np.broadcast_to(np.sin(np.pi * y)[None, :, None], shape)
    longitudinal = np.broadcast_to(np.cos(np.pi * y)[None, :, None], shape)
    bulk = np.full(shape, 3.0)
    fields = [transverse, longitudinal, bulk]
    rec = analyzer.spectrum(fields, lengths)
    assert rec["dk"] == pytest.approx(np.pi)
    assert rec["real_space_power"] == pytest.approx(1.0)
    assert rec["spectral_integral"] == pytest.approx(1.0)
    assert np.sum(rec["power"]) * rec["dk"] == pytest.approx(1.0)
    assert rec["power"][1] * rec["dk"] == pytest.approx(1.0)
    assert rec["parseval_relative_error"] < 1.0e-14
    with_bulk = analyzer.spectrum(fields, lengths, remove_mean=False)
    assert with_bulk["spectral_integral"] == pytest.approx(10.0)
    assert with_bulk["power"][0] * with_bulk["dk"] == pytest.approx(9.0)
    force = analyzer.force_products(dict(zip(analyzer.FORCE_FIELDS, fields)), lengths)
    assert force["total_acceleration_power"] == pytest.approx(1.0)
    assert force["compressive_acceleration_power"] == pytest.approx(0.5)
    assert force["solenoidal_acceleration_power"] == pytest.approx(0.5)
    assert force["solenoidal_fraction"] == pytest.approx(0.5)
    assert force["zero_mode_power"] == pytest.approx(9.0)


def test_density_weighted_kinetic_parseval_keeps_transformed_mean(analyzer):
    shape = (4, 4, 64)
    phase = 2 * np.pi * (np.arange(shape[2]) + 0.5) / shape[2]
    cosine = np.broadcast_to(np.cos(phase)[None, None, :], shape)
    sine = np.broadcast_to(np.sin(phase)[None, None, :], shape)
    fields = uniform_fields(shape)
    fields.update(dens=1.0 + 0.8 * cosine, velx=0.25 + sine,
                  vely=0.3 * cosine)
    _, spectra, _ = analyzer.snapshot_products(fields, (1, 1, 1), MODEL, 0.1)
    rho = fields["dens"]
    velocity = [fields[key] for key in ("velx", "vely", "velz")]
    # Independent integrated kinetic identity: E_fluct = E - |P|^2/(2M).
    total = 0.5 * sum(np.mean(rho * value**2) for value in velocity)
    bulk = 0.5 * sum(np.mean(rho * value)**2 for value in velocity) / np.mean(rho)
    assert spectra["kinetic"]["spectral_integral"] == pytest.approx(total - bulk,
                                                                      rel=1.0e-13)


def test_full_wavenumber_resolution_keeps_parallel_modes_and_rejects_high_kparallel(analyzer):
    shape = (32, 32, 32)
    phase = 2 * np.pi * (np.arange(32) + 0.5) / 32
    z = phase[:, None, None]
    x = phase[None, None, :]
    parallel = np.broadcast_to(np.sin(z), shape)
    unresolved = np.broadcast_to(np.sin(x + 6 * z), shape)
    low = analyzer.spectrum([parallel], (1, 1, 1))
    high = analyzer.spectrum([unresolved], (1, 1, 1))
    # The resolved band ends at mode radius 4. A pure parallel mode belongs
    # in perpendicular shell zero; a low-k_perp mode with k_parallel=6 does not.
    assert low["spectral_integral"] == pytest.approx(0.5)
    assert low["full_k_resolved_power"][0] * low["dk"] == pytest.approx(0.5)
    assert high["power"][1] * high["dk"] == pytest.approx(0.5)
    assert np.sum(high["full_k_resolved_power"]) * high["dk"] < 1.0e-27
    dc = analyzer.spectrum([np.ones(shape)], (1, 1, 1), remove_mean=False)
    assert dc["spectral_integral"] == pytest.approx(1.0)
    assert np.sum(dc["full_k_resolved_power"]) == 0.0


def test_local_strain_is_not_parallel_velocity_gradient(analyzer):
    shape = (4, 4, 128)
    phase = 2 * np.pi * (np.arange(shape[2]) + 0.5) / shape[2]
    wave = np.broadcast_to(np.sin(phase)[None, None, :], shape)
    fields = uniform_fields(shape)
    fields.update(vely=np.ones(shape), bcc1=np.ones(shape),
                  bcc2=wave, bcc3=np.zeros(shape))
    scalars, spectra, _ = analyzer.snapshot_products(fields, (1, 1, 1), MODEL, 0.1)
    assert scalars["S_parallel_rms"] == 0.0
    assert scalars["div_u_rms"] == 0.0
    assert spectra["S_parallel"]["spectral_integral"] == 0.0
    for name in ("perp_gradient_u_perp", "parallel_gradient_u_perp",
                 "perp_gradient_u_parallel"):
        assert scalars[name + "_rms"] == 0.0
        assert spectra[name]["spectral_integral"] == 0.0
    # B=(1,sin(kx),0) is divergence-free. Uniform u has zero strain, while
    # b.grad(u.b) = k*cos(kx)/(1+sin(kx)^2)^2 is nonzero due to curvature.
    analytic = 2 * np.pi * np.cos(phase) / (1 + np.sin(phase)**2)**2
    assert scalars["curvature_plus_discretization_rms"] == pytest.approx(
        np.sqrt(np.mean(analytic**2)), rel=0.003)
    assert spectra["b_grad_u_parallel"]["spectral_integral"] > 1.0


def test_pressure_balance_preserves_sign_and_residual(analyzer):
    wave = 0.3 * np.cos(2 * np.pi * (np.arange(64) + 0.5) / 64)
    magnetic = 2.0 + wave
    balanced = analyzer.pressure_balance(5.0 - wave, magnetic)
    aligned = analyzer.pressure_balance(5.0 + wave, magnetic)
    assert balanced["correlation"] == pytest.approx(-1.0)
    assert balanced["covariance"] == pytest.approx(-0.045)
    assert balanced["normalized_residual_variance"] < 1.0e-27
    assert aligned["correlation"] == pytest.approx(1.0)
    assert aligned["normalized_residual_variance"] == pytest.approx(2.0)
    empty = analyzer.pressure_balance(np.ones(8), np.ones(8))
    assert empty["correlation"] is None
    assert empty["normalized_residual_variance"] is None


def test_threshold_nearness_is_distinct_from_exceedance(analyzer):
    x = np.array([-2.5, -2.0625, -2.0, -1.9375, 0.0,
                  0.9375, 1.0, 1.0625, 1.5])
    shape = (4, 4, len(x))
    fields = uniform_fields(shape)
    fields["p_perp"] = 5.0 + 0.5 * np.broadcast_to(x, shape)
    scalars, _, distributions = analyzer.snapshot_products(
        fields, (1, 1, 1), MODEL, 0.125)
    np.testing.assert_array_equal(distributions["X"][0, 0], x)
    for name in ("mirror", "firehose"):
        assert scalars[name + "_strict"] == pytest.approx(2 / 9)
        assert scalars[name + "_inclusive"] == pytest.approx(3 / 9)
        assert scalars[name + "_near"] == pytest.approx(3 / 9)
        assert scalars[name + "_near_interior"] == pytest.approx(2 / 9)


def test_volume_pdf_and_irregular_time_weights(analyzer):
    values = np.array([0.25, 1.5, 2.5])
    edges = np.array([0.0, 1.0, 3.0])
    pdf = analyzer.volume_pdf(values, edges, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(pdf, [1 / 6, 5 / 12])
    assert np.sum(pdf * np.diff(edges)) == pytest.approx(1.0)
    rec = analyzer.temporal_summary([0.0, 1.0, 4.0], [0.0, 2.0, 8.0], 2.0)
    assert rec["mean"] == pytest.approx(4.0)
    np.testing.assert_allclose(rec["block_means"], [2.0, 6.0])
    assert rec["block_sd"] == pytest.approx(np.sqrt(8.0))
    with pytest.raises(ValueError, match="unique and increasing"):
        analyzer.time_weights([0.0, 0.0])


def write_history(path, times, values):
    lines = ["# Athena++ history data", "# [0]=time [1]=value"]
    lines += [f"{time:.17g} {value:.17g}" for time, value in zip(times, values)]
    path.write_text("\n".join(lines) + "\n")


def test_restart_history_discards_abandoned_future_and_deduplicates(analyzer, tmp_path):
    first, second = tmp_path / "first.hst", tmp_path / "second.hst"
    write_history(first, [0, 1, 2, 3], [0, 10, 20, 30])
    write_history(second, [1.5, 2.5, 2.5], [15, 25, 26])
    data, audit = analyzer.merge_histories([first, second])
    np.testing.assert_array_equal(data["time"], [0, 1, 1.5, 2.5])
    np.testing.assert_array_equal(data["value"], [0, 10, 15, 26])
    assert audit["available"]
    assert audit["files"][0]["sha256"] == hashlib.sha256(first.read_bytes()).hexdigest()
    # The same rollback can appear inside a single append-only restarted history.
    combined = tmp_path / "combined.hst"
    write_history(combined, [0, 1, 2, 3, 1.5, 2.5, 2.5],
                  [0, 10, 20, 30, 15, 25, 26])
    merged, _ = analyzer.merge_histories([combined])
    np.testing.assert_array_equal(merged["time"], data["time"])
    np.testing.assert_array_equal(merged["value"], data["value"])


def test_missing_energy_inputs_and_file_provenance_are_explicit(analyzer, tmp_path):
    data, audit = analyzer.merge_histories([])
    assert data == {} and audit["available"] is False and audit["reason"]
    result = analyzer.history_products({}, {}, {"record_injected_work": True}, 0, 1, 1)
    assert result["available"] is False
    assert result["energy_budget"]["available"] is False
    assert result["reason"]
    with pytest.raises(ValueError, match="retained file not found"):
        analyzer.paths_from_metadata(tmp_path, {"snapshots": ["absent.bin"]}, "snapshots")
    retained = tmp_path / "metadata.json"
    retained.write_bytes(b'{"revision":"synthetic-only"}\n')
    info = analyzer.retained_file(retained)
    assert info["sha256"] == hashlib.sha256(retained.read_bytes()).hexdigest()
    assert info["bytes"] == retained.stat().st_size
    assert Path(info["path"]) == retained.resolve()


def test_canonical_model_uses_amplitude_blend_and_rejects_conflicting_metadata(analyzer):
    path = REPOSITORY / "inputs/cgl_lf_paper/cgl_lf_physics_benchmark_beta10.athinput"
    header = analyzer.parse_input(path.read_text().splitlines())
    model = analyzer.infer_model(header, {})
    assert model["B0"] == 1.0
    assert model["mirror_X"] == 1.0
    assert model["firehose_X"] == -2.0
    assert model["expected_solenoidal_power_fraction"] == pytest.approx(0.5)
    with pytest.raises(ValueError, match="disagrees"):
        analyzer.infer_model(header, {"model": {"B0": 2.0}})
    for key in ("nu_coll", "limiter_nu_coll", "backup_limiters", "mirror_limiter",
                "firehose_limiter", "passive"):
        incomplete = {section: values.copy() for section, values in header.items()}
        del incomplete["mhd"][key]
        with pytest.raises(ValueError, match="required model fields"):
            analyzer.infer_model(incomplete, {})


def test_energy_balance_uses_measured_work_with_a_signed_residual(analyzer):
    user = {"time": np.array([0.0, 1.0, 2.0]),
            "volume": np.full(3, 2.0), "force_work": np.array([0.0, 0.3, 1.0])}
    mhd = {"time": np.array([0.0, 1.0, 2.0]),
           "tot-E": np.array([16.0, 16.2, 16.8])}
    model = {"record_injected_work": True, "dedt_per_volume": 1000.0,
             "passive": False}
    result = analyzer.history_products(user, mhd, model, 0, 2, 1)
    budget = result["energy_budget"]
    assert budget["available"]
    assert budget["actual_applied_work"] == pytest.approx(1.0)
    assert budget["actual_mean_total_power"] == pytest.approx(0.5)
    assert budget["conserved_total_energy_change"] == pytest.approx(0.8)
    assert budget["residual_E_minus_work"] == pytest.approx(-0.2)
    without_work = analyzer.history_products(
        {"time": user["time"], "volume": user["volume"]}, mhd, model, 0, 2, 1)
    assert without_work["energy_budget"]["available"] is False
    assert without_work["energy_budget"]["reason"]
    unknown_mode = {key: value for key, value in model.items() if key != "passive"}
    unknown = analyzer.history_products(user, mhd, unknown_mode, 0, 2, 1)
    assert unknown["energy_budget"]["available"] is False


def write_binary(path, time, fields, header, lengths=(1.0, 1.0, 1.0)):
    """One real v1.1 MeshBlock, following src/outputs/binary.cpp's disk layout."""
    shape = next(iter(fields.values())).shape
    parameters = "".join(
        f"<{block}>\n" + "".join(f"{key}={value}\n" for key, value in items.items())
        for block, items in header.items()).encode()
    names = list(fields)
    prefix = (
        "Athena binary output version=1.1\n"
        "  size of preheader=5\n"
        f"  time={time:.17g}\n  cycle={int(time * 100)}\n"
        "  size of location=8\n  size of variable=4\n"
        f"  number of variables={len(names)}\n"
        f"  variables: {' '.join(names)}\n"
        f"  header offset={len(parameters)}\n").encode()
    ng = int(header["mesh"]["nghost"])
    indices = [index for n in shape[::-1] for index in (ng, ng + n - 1)]
    block = struct.pack("<10i", *indices, 0, 0, 0, 0)
    block += struct.pack("<6d", *[value for length in lengths for value in (0, length)])
    block += b"".join(np.asarray(fields[name], dtype="<f4").tobytes() for name in names)
    path.write_bytes(prefix + parameters + block)


@pytest.fixture
def synthetic_run(analyzer, tmp_path):
    """Diagnostic data only: these fields were not evolved by the solver."""
    run = tmp_path / "synthetic_not_a_simulation"
    (run / "bin").mkdir(parents=True)
    canonical = REPOSITORY / "inputs/cgl_lf_paper/cgl_lf_physics_benchmark_beta10.athinput"
    header = analyzer.parse_input(canonical.read_text().splitlines())
    header["job"]["basename"] = "synthetic"
    for axis in (1, 2, 3):
        header["mesh"][f"nx{axis}"] = "8"
        header["mesh"][f"x{axis}max"] = "1.0"
        header["meshblock"][f"nx{axis}"] = "8"
    phase = 2 * np.pi * (np.arange(8) + 0.5) / 8
    wave = np.broadcast_to(np.cos(phase)[None, None, :], (8, 8, 8))
    user_rows, mhd_rows = [], []
    for index, time in enumerate((0.0, 1.0, 2.0)):
        fields = uniform_fields((8, 8, 8))
        fields.update(velx=(0.1 + 0.01 * time) * wave,
                      vely=0.1 * wave, bcc3=np.sqrt(1.0 + 0.1 * wave),
                      p_perp=5.0 - 0.05 * wave)
        force = {"force1": wave, "force2": wave, "force3": np.zeros_like(wave)}
        write_binary(run / "bin" / f"synthetic.mhd_w_bcc.{index:05d}.bin",
                     time, fields, header)
        write_binary(run / "bin" / f"synthetic.turb_force.{index:05d}.bin",
                     time, force, header)
        kinetic = 0.25 * ((0.1 + 0.01 * time)**2 + 0.1**2)
        energy = 8.0 + kinetic
        work = kinetic - 0.005
        user_rows.append(f"{time} 1 {kinetic:.17g} .5 7.5 {work:.17g}")
        mhd_rows.append(f"{time} {energy:.17g}")
    (run / "synthetic.user.hst").write_text(
        "# Athena++ history data\n"
        "# [0]=time [1]=volume [2]=kinetic [3]=magnetic [4]=therm_cgl [5]=force_work\n"
        + "\n".join(user_rows) + "\n")
    (run / "synthetic.mhd.hst").write_text(
        "# Athena++ history data\n# [0]=time [1]=tot-E\n"
        + "\n".join(mhd_rows) + "\n")
    metadata = {"purpose": "Synthetic analyzer regression; no solver was run",
                "simulation": {}, "outputs": {}}
    (run / "benchmark_metadata.json").write_text(json.dumps(metadata))
    return run, header


def run_cli(run, output, metadata=None, start=0, end=2):
    command = [sys.executable, str(ANALYZER), str(run), "--time-start", str(start),
               "--time-end", str(end), "--block-duration", "0.5",
               "--output-dir", str(output)]
    if metadata is not None:
        command += ["--metadata", str(metadata)]
    proc = subprocess.run(command, cwd=run, capture_output=True, text=True, timeout=90)
    (run / (output.name + "-stdout.log")).write_text(proc.stdout + proc.stderr)
    return proc


def test_real_binary_reader_and_cli_export_with_unknown_simulation(analyzer, synthetic_run):
    run, _ = synthetic_run
    output = run / "analysis-complete"
    proc = run_cli(run, output)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    result = json.loads((output / "metrics.json").read_text())
    assert result["sampling"]["times"] == [0.0, 1.0, 2.0]
    assert result["snapshots"][0]["info"]["shape_zyx"] == [8, 8, 8]
    assert result["snapshots"][0]["info"]["lengths_xyz"] == [1.0, 1.0, 1.0]
    assert result["provenance"]["simulation_revision"] is None
    assert "analysis_revision" in result["provenance"]
    assert result["adequacy"]["classification"] == "inconclusive"
    assert any("revision" in reason for reason in result["adequacy"]["reasons"])
    assert result["history"]["energy_budget"]["available"]
    assert abs(result["history"]["energy_budget"]["residual_E_minus_work"]) < 1.0e-14
    assert result["forcing_decomposition"]["available"]
    assert result["model"]["expected_solenoidal_power_fraction"] == pytest.approx(0.5)
    for name in ("marginality", "pressure_balance", "gradients", "spectra_energy"):
        for suffix in ("png", "pdf"):
            assert (output / f"{name}.{suffix}").stat().st_size > 100
    assert (output / "report.md").is_file()
    # Missing ledgers/force snapshots remain explicit, even when primitive
    # snapshots and plotting are otherwise usable.
    metadata = {"purpose": "Synthetic missing-data regression",
                "simulation": {}, "outputs": {
                    "user_history": [], "mhd_history": [], "force_snapshots": []}}
    missing = run / "missing_metadata.json"
    missing.write_text(json.dumps(metadata))
    absent_output = run / "analysis-missing"
    proc = run_cli(run, absent_output, missing)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    absent = json.loads((absent_output / "metrics.json").read_text())
    assert absent["history"]["energy_budget"]["available"] is False
    assert absent["forcing_decomposition"]["available"] is False
    assert absent["forcing_decomposition"]["reason"]
    assert absent["adequacy"]["classification"] == "inconclusive"


def test_real_binary_reader_rejects_missing_primitive_fields(analyzer, synthetic_run):
    run, header = synthetic_run
    fields = uniform_fields((8, 8, 8))
    del fields["p_perp"]
    path = run / "missing-pressure.bin"
    write_binary(path, 0.0, fields, header)
    with pytest.raises(ValueError, match="missing fields"):
        analyzer.read_uniform(path, analyzer.FIELDS)


def test_cli_explains_undefined_startup_metrics_without_shrinking_window(synthetic_run):
    run, header = synthetic_run
    write_binary(run / "bin" / "synthetic.mhd_w_bcc.00000.bin", 0.0,
                 uniform_fields((8, 8, 8)), header)
    zero_force = {name: np.zeros((8, 8, 8))
                  for name in ("force1", "force2", "force3")}
    write_binary(run / "bin" / "synthetic.turb_force.00000.bin", 0.0,
                 zero_force, header)
    output = run / "analysis-undefined-startup"
    proc = run_cli(run, output, start=0.5, end=2.0)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    result = json.loads((output / "metrics.json").read_text())
    assert result["requested_window"] == [0.5, 2.0]
    assert result["retained_window"] == [0.5, 2.0]
    assert result["sampling"]["times"] == [0.0, 1.0, 2.0]
    assert result["scalars"]["thermal_density"]["mean"] == pytest.approx(7.5)
    assert result["scalars"]["u_parallel_rms"]["mean"] == 0.0
    reasons = result["unavailable_metric_reasons"]
    names = ("scalars.pressure_correlation",
             "scalars.pressure_normalized_residual_variance",
             "forcing_decomposition.solenoidal_fraction")
    report = (output / "report.md").read_text()
    for name in names:
        assert reasons[name]["available"] is False
        assert reasons[name]["undefined_sample_times"] == [0.0]
        assert reasons[name]["undefined_bracketing_times"] == [0.0]
        assert reasons[name]["requested_window"] == [0.5, 2.0]
        assert reasons[name]["reason"]
        assert [row["time"] for row in reasons[name]["valid_samples"]] == [1.0, 2.0]
        assert name in report
    for snapshot in result["snapshots"][1:]:
        assert snapshot["scalars"]["pressure_correlation"] == pytest.approx(-1.0)
        assert snapshot["scalars"]["pressure_normalized_residual_variance"] < 1.0e-9
    force = result["forcing_decomposition"]
    assert force["snapshots"][0]["values"]["solenoidal_fraction"] is None
    for snapshot in force["snapshots"][1:]:
        assert snapshot["values"]["solenoidal_fraction"] == pytest.approx(0.5)
    assert "solenoidal_fraction" not in force["statistics"]
    # P(t) rises linearly from zero at t=0 to P at t=1, then stays at P.
    # Its integral on [.5,2] is (3/8+1)*P, so the mean is (11/12)*P.
    phase = 2 * np.pi * (np.arange(8) + 0.5) / 8
    stored_wave = np.cos(phase).astype(np.float32).astype(float)
    later_power = 2 * np.mean(stored_wave**2)
    assert force["statistics"]["total_acceleration_power"]["mean"] == pytest.approx(
        (11 / 12) * later_power, rel=1.0e-13)


def test_cli_preserves_reserved_projection_counter_without_activity_inference(
        analyzer, synthetic_run):
    run, _ = synthetic_run
    revision = "7a37710f6c224e24e7c7f364e7e0b812b3a9494c"
    path = run / "benchmark_metadata.json"
    metadata = json.loads(path.read_text())
    metadata.update(simulation={"revision": revision}, launch={"returncode": 0})
    path.write_text(json.dumps(metadata))
    history = run / "synthetic.mhd.hst"
    lines = history.read_text().splitlines()
    counters = ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos",
                "lf_hardbd", "lf_hwproj")
    lines[1] += " " + " ".join(f"[{index}]={name}"
                               for index, name in enumerate(counters, 2))
    # Deliberately synthetic nonzero reserved values must be retained verbatim;
    # they cannot establish that the audited executable measured projections.
    for index, (hardbd, reserved) in enumerate(zip((10, 15, 21), (7, 9, 14)), 2):
        lines[index] += f" 0 0 0 0 {hardbd} {reserved}"
    history.write_text("\n".join(lines) + "\n")
    output = run / "analysis-reserved-counter"
    proc = run_cli(run, output, start=0.5, end=2.0)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    result = json.loads((output / "metrics.json").read_text())
    health = result["simulation_integrity"]
    reserved = health["physical_stage_activity"]["lf_hwproj"]
    assert reserved["available"]
    assert reserved["full_retained_min"] == 7.0
    assert reserved["full_retained_max"] == 14.0
    assert reserved["last"] == 14.0
    assert reserved["window_increment"] == 6.0
    assert reserved["instrumentation"] == "reserved_uninstrumented"
    assert reserved["audited_simulation_revision"] == revision
    assert "never increments" in reserved["meaning"]
    assert "absence" in reserved["meaning"]
    assert health["classification"] == "consistent"
    assert health["reasons"] == []
    assert "restart-persistent" in health["LF_counter_lifecycle"]
    hardbd = health["physical_stage_activity"]["lf_hardbd"]
    assert hardbd["window_increment"] == 8.5
    assert "repeated cell-stage events" in hardbd["meaning"]
    definitions = result["definitions"]["LF_counters"]
    assert all(word in definitions for word in ("uninstrumented", "halo", "active", "restart"))
    report = (output / "report.md").read_text()
    assert "Projection-count limitation" in report
    assert "zero does not demonstrate absence of projections" in report
    assert "including refreshed halo cells" in report
    # A different or unknown revision must not inherit the audited label.
    data, _ = analyzer.merge_histories([history])
    unknown = analyzer.solver_health(
        [{"directory": run, "metadata": {"simulation": {"revision": "unknown"},
                                          "launch": {"returncode": 0}}}],
        {}, data, 0.5, 2.0)
    unknown_counter = unknown["physical_stage_activity"]["lf_hwproj"]
    assert unknown_counter["instrumentation"] == "not_verified_for_retained_revision"
    for name in ("full_retained_min", "full_retained_max", "last", "window_increment"):
        assert unknown_counter[name] == reserved[name]


@pytest.mark.parametrize("defect", ["time", "geometry", "model"])
def test_cli_rejects_unmatched_force_snapshot(synthetic_run, defect):
    run, header = synthetic_run
    force = {name: np.ones((8, 8, 8)) for name in ("force1", "force2", "force3")}
    time, lengths = 1.0, (1.0, 1.0, 1.0)
    if defect == "time":
        time = 1.125
    elif defect == "geometry":
        lengths = (2.0, 1.0, 1.0)
        header["mesh"]["x1max"] = "2.0"
    else:
        header["turb_driving"]["dedt"] = "0.64"
    write_binary(run / "bin" / "synthetic.turb_force.00001.bin", time, force,
                 header, lengths)
    proc = run_cli(run, run / "analysis-rejected")
    assert proc.returncode != 0, "unmatched forcing data silently accepted"
    assert any(word in (proc.stdout + proc.stderr).lower()
               for word in ("match", "geometry", "forcing", "force"))
