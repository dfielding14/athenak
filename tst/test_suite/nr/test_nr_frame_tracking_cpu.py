"""Regression coverage for serial non-relativistic frame tracking."""

from pathlib import Path
import shutil
from subprocess import PIPE, run

import athena_read
import numpy as np
import pytest
import test_suite.testutils as testutils


HYDRO_INPUT = "inputs/frame_tracking_hydro.athinput"
MHD_INPUT = "inputs/frame_tracking_mhd.athinput"
THREE_D_INPUT = "inputs/frame_tracking_3d.athinput"
AMBIGUOUS_INPUT = "inputs/frame_tracking_invalid_ambiguous.athinput"
LEGACY_INPUT = "inputs/frame_tracking_legacy.athinput"


def remove_outputs(basename: str) -> None:
    for path in Path(".").glob(f"{basename}.*.hst"):
        path.unlink()
    for path in Path("rst").glob(f"{basename}.*.rst") if Path("rst").exists() else []:
        path.unlink()
    for path in Path("tab").glob(f"{basename}.*.tab") if Path("tab").exists() else []:
        path.unlink()


def input_with_frame_key(tmp_path: Path, key: str, value: str) -> Path:
    modified = Path(HYDRO_INPUT).read_text().replace(
        "<output1>", f"{key} = {value}\n\n<output1>", 1
    )
    path = tmp_path / f"frame_tracking_{key}.athinput"
    path.write_text(modified)
    return path


def test_hydro_structured_history() -> None:
    basename = "FrameTrackHydro"
    try:
        assert testutils.run(HYDRO_INPUT)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        required = {
            "ft_vf_x1",
            "ft_dx_x1",
            "ft_pos_x1",
            "ft_err_x1",
            "ft_dv_x1",
            "ft_vf_x2",
            "ft_weight",
            "ft_misses",
            "ft_limit",
            "ft_skip",
        }
        assert required.issubset(data)
        assert np.all(np.isfinite(data["ft_vf_x1"]))
        assert np.any(np.abs(data["ft_vf_x1"]) > 0.0)
        assert np.max(data["ft_misses"]) == 0.0
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


def test_mhd_boost_preserves_fields_and_updates_conserved_state() -> None:
    basename = "FrameTrackMHD"
    try:
        assert testutils.run(MHD_INPUT)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        assert np.all(np.isfinite(data["ft_vf_x1"]))
        assert np.any(np.abs(data["ft_vf_x1"]) > 0.0)
        assert np.max(data["ft_misses"]) == 0.0
        active = np.flatnonzero(np.abs(data["ft_dv_x1"]) > 0.0)
        assert active.size
        index = int(active[0])
        tabs = sorted(Path("tab").glob(f"{basename}.mhd_u_bcc.*.tab"))
        assert index > 0 and len(tabs) > index
        before = athena_read.tab(str(tabs[index - 1]))
        after = athena_read.tab(str(tabs[index]))
        dv = data["ft_dv_x1"][index]
        tol = 100.0 * np.finfo(float).eps
        np.testing.assert_allclose(after["dens"], before["dens"], rtol=tol, atol=tol)
        for field in ("bcc1", "bcc2", "bcc3"):
            np.testing.assert_allclose(after[field], before[field], rtol=tol, atol=tol)
        np.testing.assert_allclose(
            after["mom1"], before["mom1"] + before["dens"] * dv, rtol=tol, atol=tol
        )
        np.testing.assert_allclose(after["mom2"], before["mom2"], rtol=tol, atol=tol)
        np.testing.assert_allclose(after["mom3"], before["mom3"], rtol=tol, atol=tol)
        np.testing.assert_allclose(
            after["ener"],
            before["ener"] + before["mom1"] * dv + 0.5 * before["dens"] * dv**2,
            rtol=tol,
            atol=tol,
        )
    finally:
        remove_outputs(basename)
        testutils.cleanup()


def test_three_axis_history_has_complete_schema_and_actuation() -> None:
    basename = "FrameTrack3D"
    try:
        assert testutils.run(THREE_D_INPUT)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        expected = {
            f"ft_{quantity}_x{axis}"
            for axis in (1, 2, 3)
            for quantity in ("vf", "dx", "pos", "err", "dv")
        } | {"ft_weight", "ft_misses", "ft_recov", "ft_limit", "ft_skip"}
        assert expected.issubset(data)
        assert len(expected) == 20
        for axis in (1, 2, 3):
            assert np.any(np.abs(data[f"ft_dv_x{axis}"]) > 0.0)
        assert np.max(data["ft_misses"]) == 0.0
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


@pytest.mark.parametrize("args", [["coord/special_rel=true"], ["coord/general_rel=true"]])
def test_relativistic_configurations_are_rejected(args: list[str]) -> None:
    result = run(
        ["./athena", "-i", HYDRO_INPUT, *args],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "supports non-relativistic hydro/MHD only" in output


def test_dynamical_relativistic_configuration_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "frame_tracking_dyn_gr.athinput"
    path.write_text(
        Path(HYDRO_INPUT).read_text() + "\n<adm>\nframe_tracking_guard = true\n"
    )
    result = run(
        ["./athena", "-i", str(path)],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "supports non-relativistic hydro/MHD only" in output


def test_requested_missing_fluid_is_rejected() -> None:
    result = run(
        ["./athena", "-i", MHD_INPUT, "frame_tracking/tracked_fluid=hydro"],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    assert result.returncode != 0
    assert "does not select an available hydro or mhd fluid" in (
        result.stdout + result.stderr
    )


def test_ambiguous_dual_fluid_configuration_is_rejected() -> None:
    result = run(
        ["./athena", "-i", AMBIGUOUS_INPUT],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "set tracked_fluid=hydro or tracked_fluid=mhd explicitly" in output


@pytest.mark.parametrize(
    ("removed_key", "value", "canonical_key"),
    [
        ("enable", "true", "enabled"),
        ("target_variable", "density", "target"),
        ("axis", "x1", "axes"),
        ("ft_tau_avg", "0.1", "tau_avg"),
        ("ft_max_boost_change_rate", "0.1", "max_boost_change_rate"),
        ("ndiag", "1", "diagnostic_every"),
        ("z_target", "0.0", "x3_target"),
    ],
)
def test_removed_aliases_fail_with_canonical_replacement(
    tmp_path: Path, removed_key: str, value: str, canonical_key: str
) -> None:
    input_path = input_with_frame_key(tmp_path, removed_key, value)
    result = run(
        ["./athena", "-i", str(input_path)],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert f"Removed key <frame_tracking>/{removed_key}; use {canonical_key}." in output


def test_initialization_summary_records_interpreted_configuration() -> None:
    basename = "FrameTrackSummary"
    try:
        result = run(
            ["./athena", "-i", HYDRO_INPUT, f"job/basename={basename}"],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 0
        assert "FrameTracker configuration: fluid=hydro axes=x1,x2" in output
        assert "target=density" in output
        assert "slew=per_time state=new" in output
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


def test_legacy_frame_state_warns_and_remains_finite() -> None:
    basename = "FrameTrackLegacy"
    try:
        result = run(
            [
                "./athena",
                "-i",
                LEGACY_INPUT,
                f"job/basename={basename}",
            ],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        assert result.returncode == 0
        assert "restoring legacy frame state without complete controller history" in (
            result.stdout + result.stderr
        )
        assert "slew=per_time state=legacy_restart" in (result.stdout + result.stderr)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        assert np.all(np.isfinite(data["ft_vf_x1"]))
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()
