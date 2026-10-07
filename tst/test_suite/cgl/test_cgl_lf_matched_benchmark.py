"""Matched comparison contracts; synthetic fields are not solver evidence."""

import copy
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from test_cgl_lf_physics_benchmark import uniform_fields, write_binary


REPOSITORY = Path(__file__).resolve().parents[3]
COMPARISON = REPOSITORY / "scripts/compare_cgl_lf_physics_benchmark.py"
INPUT = REPOSITORY / "inputs/cgl_lf_paper/cgl_lf_physics_benchmark_matched_beta10.athinput"
LAUNCHER = REPOSITORY / "scripts/run_cgl_lf_matched_benchmark.py"


@pytest.fixture(scope="module")
def comparison():
    name = "cgl_lf_matched_benchmark_regression"
    spec = importlib.util.spec_from_file_location(name, COMPARISON)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def launcher():
    name = "cgl_lf_matched_launcher_regression"
    spec = importlib.util.spec_from_file_location(name, LAUNCHER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("interrupted", [True, False])
def test_launcher_histories_survive_interruption_and_finalize_exactly(
        launcher, comparison, monkeypatch, tmp_path, interrupted):
    """Exercise actual launch metadata writes and analyzer discovery without Slurm."""
    binary = tmp_path / "synthetic-athena"
    binary.write_text("Synthetic executable fixture; never executed.\n")
    manifest = tmp_path / "binary-manifest.json"
    manifest.write_text(json.dumps({"revision": "retained-simulation-revision",
                                    "binary_sha256": launcher.sha(binary)}))
    cache = tmp_path / "CMakeCache.txt"
    cache.write_text("PROBLEM:STRING=built_in_pgens\nKokkos_ENABLE_HIP:BOOL=ON\n"
                     "Athena_ENABLE_MPI:BOOL=ON\nAthena_SINGLE_PRECISION:BOOL=OFF\n")
    run = tmp_path / "run"
    monkeypatch.setattr(sys, "argv", [str(LAUNCHER), str(run), "--mode", "active",
        "--input", str(INPUT), "--executable", str(binary),
        "--build-manifest", str(manifest), "--build-cache", str(cache),
        "--job-id", "synthetic-no-job-submitted"])

    def git_result(command, **kwargs):
        assert command[:1] == ["git"]
        if "diff" in command:
            return ""
        assert command[-2:] == ["rev-parse", "HEAD"]
        return "launcher-checkout-revision\n"

    def simulated_slurm(command, *, cwd, stdout, stderr):
        assert command[0] == "srun" and cwd == run
        (run / "hst").mkdir()
        for kind in ("user", "mhd"):
            (run / "hst" / f"fixture.{kind}.hst").write_text(
                "# Athena++ history data\n# [1]=time [2]=value\n0 1\n1 2\n")
        if interrupted:
            raise KeyboardInterrupt("simulated launcher interruption")
        return launcher.subprocess.CompletedProcess(command, 0)

    with monkeypatch.context() as launch_tools:
        launch_tools.setattr(launcher.subprocess, "check_output", git_result)
        launch_tools.setattr(launcher.subprocess, "run", simulated_slurm)
        if interrupted:
            with pytest.raises(KeyboardInterrupt, match="simulated launcher interruption"):
                launcher.main()
        else:
            assert launcher.main() == 0
    if not interrupted:
        # A completed segment's exact inventory must not silently grow later.
        for kind in ("user", "mhd"):
            (run / "hst" / f"unretained.{kind}.hst").write_text(
                "# Athena++ history data\n# [1]=time [2]=value\n2 999\n")

    metadata_path = run / "benchmark_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    assert metadata["simulation"]["revision"] == "retained-simulation-revision"
    assert metadata["launch"]["source_checkout_revision"] == "launcher-checkout-revision"
    assert ("returncode" in metadata["launch"]) == (not interrupted)
    segments = comparison.single.load_segments(run, metadata, metadata_path)
    for kind in ("user_history", "mhd_history"):
        suffix = kind.removesuffix("_history")
        if interrupted:
            assert kind not in metadata["outputs"]
        else:
            assert metadata["outputs"][kind] == [f"hst/fixture.{suffix}.hst"]
        paths, boundaries, _ = comparison.single.segment_paths(
            segments, kind, f"**/*.{suffix}.hst")
        assert paths == [run / "hst" / f"fixture.{suffix}.hst"]
        history, provenance = comparison.single.merge_histories(paths, boundaries)
        assert provenance["available"]
        np.testing.assert_array_equal(history["time"], [0., 1.])
        np.testing.assert_array_equal(history["value"], [1., 2.])


@pytest.fixture
def headers(comparison):
    active = comparison.single.parse_input(INPUT.read_text().splitlines())
    passive = copy.deepcopy(active)
    passive["mhd"]["passive"] = "true"
    passive["problem"]["passive_delta"] = "true"
    return active, passive


def test_matching_allows_runtime_administration_and_numeric_spellings(comparison, headers):
    active, passive = headers
    passive["job"]["basename"] = "separate_passive_output"
    passive["output2"].update(file_number="99", dt="0.125")
    passive["time"].update(tlim="20", restart_time="4")
    passive["mhd"]["passive_restart_encoding"] = "1"
    passive["time"]["cfl_number"] = "3e-1"
    result = comparison.compare_parameters(active, passive)
    assert result["matched"]
    exemptions = {row["parameter"] for row in result["exempted_differences"]}
    assert {"job/basename", "output2/dt", "time/tlim", "time/restart_time",
            "mhd/passive", "problem/passive_delta"} <= exemptions
    assert "time/cfl_number" not in exemptions


@pytest.mark.parametrize("section,key,value", [
    ("mhd", "reconstruct", "plm"),
    ("mesh", "nghost", "2"),
    ("mesh", "nx1", "96"),
    ("mhd", "iso_sound_speed", "1"),
    ("time", "cfl_number", "0.300000000000001"),
    ("time", "sts_integrator", "rkc2"),
    ("turb_driving", "rseed", "314159"),
    ("mhd", "lf_k_parallel", "3"),
])
def test_matching_rejects_changed_physics_and_numerics(comparison, headers, section, key, value):
    active, passive = headers
    passive[section][key] = value
    with pytest.raises(ValueError, match="Unmatched effective parameters") as error:
        comparison.compare_parameters(active, passive)
    assert f"{section}/{key}" in str(error.value)


def test_matching_rejects_wrong_modes_and_conflicting_aliases(comparison, headers):
    active, passive = headers
    with pytest.raises(ValueError, match="must have mhd/passive"):
        comparison.compare_parameters(passive, active)
    with pytest.raises(ValueError, match="must have mhd/passive"):
        comparison.compare_parameters(active, active)
    passive["problem"]["passive_delta"] = "false"
    with pytest.raises(ValueError, match="passive_delta disagrees"):
        comparison.compare_parameters(active, passive)


def test_matching_requires_sound_speed_and_supported_restart_encoding(comparison, headers):
    active, passive = headers
    missing = copy.deepcopy(active)
    del missing["mhd"]["iso_sound_speed"]
    with pytest.raises(ValueError, match="iso_sound_speed"):
        comparison.compare_parameters(missing, passive)
    passive["mhd"]["passive_restart_encoding"] = "2"
    with pytest.raises(ValueError, match="unsupported passive_restart_encoding"):
        comparison.compare_parameters(active, passive)


def test_pdf_prescan_uses_true_union_and_endpoint_brackets(comparison, headers, tmp_path):
    scans, retained_samples = {}, {"B_over_B0": [], "X": []}
    # A single-cell tail in each real binary guards against per-run histogram
    # rebinning or sampled extrema. t=0 and t=2 bracket the requested interval;
    # the much larger t=3 values must not broaden the selected interval's bins.
    tails = {"active": [(1., -.5), (2., .75), (1., 0.), (16., 20.)],
             "passive": [(.5, -3.), (1., 0.), (4., 12.), (32., 30.)]}
    for mode, original in zip(("active", "passive"), headers):
        run = tmp_path / mode
        (run / "bin").mkdir(parents=True)
        header = copy.deepcopy(original)
        for axis in (1, 2, 3):
            header["mesh"].update({f"nx{axis}": "8", f"x{axis}max": "1"})
            header["meshblock"][f"nx{axis}"] = "8"
        (run / "benchmark_metadata.json").write_text(json.dumps({
            "classification": "synthetic plumbing only; no solver was run",
            "simulation": {}, "outputs": {}}))
        for index, (magnetic, anisotropy) in enumerate(tails[mode]):
            fields = uniform_fields((8, 8, 8))
            fields["bcc3"].flat[0] = magnetic
            fields["p_perp"].flat[0] = 5 + .5*anisotropy*magnetic**2
            write_binary(run / "bin" / f"synthetic.mhd_w_bcc.{index:05d}.bin",
                         float(index), fields, header)
            if index < 3:
                retained_samples["B_over_B0"].extend(fields["bcc3"].ravel())
                retained_samples["X"].extend(
                    (2*(fields["p_perp"]-fields["eint"])/fields["bcc3"]**2).ravel())
        scans[mode] = comparison.prescan(run, .25, 1.75, mode)
        assert [row["info"]["time"] for row in scans[mode]["snapshots"]] == [0., 1., 2.]
    edges = comparison.shared_edges(scans, 32)
    assert edges["observed_union_limits"] == {"B_over_B0": [.5, 4.], "X": [-3., 12.]}
    for name, samples in retained_samples.items():
        bins = np.asarray(edges["edges"][name])
        counts, _ = np.histogram(samples, bins)
        assert bins[0] < min(samples) and bins[-1] > max(samples)
        assert counts.sum() == len(samples)
        assert counts[0] and counts[-1]
        pdf = comparison.single.volume_pdf(np.asarray(samples), bins, 1.)
        np.testing.assert_allclose(pdf*np.diff(bins), counts/len(samples), atol=1.e-15)
        assert pdf@np.diff(bins) == pytest.approx(1., abs=1.e-14)


def test_block_contrast_matches_physical_time_despite_different_snapshot_cadence(comparison):
    active_times = np.array([0., .7, 2., 3.6, 4.])
    passive_times = np.array([0., 1.2, 2.8, 4.])
    # Linear trajectories have exact time means independently of cadence.
    active = comparison.single.temporal_summary(active_times, active_times, 2., 0., 4.)
    passive = comparison.single.temporal_summary(passive_times, .5*passive_times+.5, 2., 0., 4.)
    contrast = comparison.summary_difference(active, passive)
    assert contrast["active_minus_passive"] == pytest.approx(.5)
    assert contrast["active_over_passive"] == pytest.approx(4/3)
    assert contrast["block_difference_means"] == pytest.approx([0., 1.])
    assert contrast["block_difference_sd"] == pytest.approx(1/np.sqrt(2))
    displaced = dict(passive, block_starts=[.5, 2.5])
    unpaired = comparison.summary_difference(active, displaced)
    assert unpaired["block_difference_means"] == []
    assert unpaired["block_difference_sd"] is None
    assert unpaired["block_comparison_reason"]


def test_undefined_pressure_bins_are_not_silently_dropped_from_blocks(comparison):
    active = dict(mean=[2., None], block_means=[[1., None], [3., 2.]],
                  block_starts=[0., 2.], block_duration=2.)
    passive = dict(mean=[1., 0.], block_means=[[.5, 0.], [1.5, 0.]],
                   block_starts=[0., 2.], block_duration=2.)
    contrast = comparison.summary_difference(active, passive)
    assert contrast["active_minus_passive"] == [1., None]
    assert contrast["active_over_passive"] == [2., None]
    assert contrast["block_difference_min"] == [.5, None]
    assert contrast["block_difference_max"] == [1.5, None]
    assert contrast["block_difference_sd"][0] == pytest.approx(1/np.sqrt(2))
    assert contrast["block_difference_sd"][1] is None
    assert json.loads(json.dumps(contrast, allow_nan=False)) == contrast
    scalar = comparison.summary_difference(
        dict(mean=2., block_means=[1., 3.], block_starts=[0., 2.], block_duration=2.),
        dict(mean=0., block_means=[0., 0.], block_starts=[0., 2.], block_duration=2.))
    assert scalar["active_minus_passive"] == 2.
    assert scalar["active_over_passive"] is None
    json.dumps(scalar, allow_nan=False)


def test_failed_segment_provenance_remains_visible_after_successful_segment(comparison):
    data = {mode: {"provenance": {"retained_metadata": {"classification": "plumbing only"}},
                   "simulation_integrity": {"segments": [{"returncode": 0}]}}
            for mode in ("active", "passive")}
    assert "FAILED RUN" not in comparison.figure_status(data)
    data["passive"]["simulation_integrity"]["segments"].insert(0, {"returncode": 1})
    status = comparison.figure_status(data)
    assert "FAILED RUN: passive" in status and "PLUMBING ONLY" in status
    evidence = comparison.physical_evidence(data, {"failed_runs": ["passive"]})
    assert any("Failed simulation provenance: passive" in reason for reason in evidence["reasons"])
    json.dumps(evidence, allow_nan=False)


def test_healthy_plumbing_does_not_automatically_validate_physics(comparison):
    data = {mode: {"adequacy": {"reasons": []}} for mode in ("active", "passive")}
    evidence = comparison.physical_evidence(data, {"classification": "consistent", "failed_runs": []})
    assert evidence["classification"] == "inconclusive"
    assert not evidence["scientific_review_completed"]
    assert evidence["reasons"]
