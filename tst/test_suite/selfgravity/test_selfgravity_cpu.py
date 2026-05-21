"""
Regression tests for multigrid self-gravity.

These tests exercise the periodic Jeans-wave pgen in both hydro and MHD mode.
They keep the mesh small but use multiple MeshBlocks so the multigrid pack and
boundary paths are active.
"""

from pathlib import Path
import re
import subprocess

import pytest

import test_suite.testutils as testutils


_DEFECT_RE = re.compile(r"Final defect norm = ([0-9.eE+-]+)")
_REPO_ROOT = Path(__file__).resolve().parents[3]


def _run_selfgravity(input_file, basename):
    command = [
        "./athena",
        "-i",
        str(_REPO_ROOT / input_file),
        f"job/basename={basename}",
        "gravity/show_defect=true",
    ]
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, check=False)
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "nan" not in combined.lower()
    matches = _DEFECT_RE.findall(combined)
    assert matches, combined
    assert float(matches[-1]) < 2.0e-8
    assert "k/k_J = 2" in combined
    return combined


def _remove_binary_outputs(basename):
    bin_dir = Path("bin")
    if not bin_dir.exists():
        return
    for path in bin_dir.glob(f"{basename}*.bin"):
        path.unlink(missing_ok=True)
    try:
        bin_dir.rmdir()
    except OSError:
        pass


@pytest.mark.parametrize(
    "input_file,basename",
    [
        ("inputs/tests/selfgravity.athinput", "selfgravity_hydro_cpu"),
        ("inputs/tests/selfgravity_mhd.athinput", "selfgravity_mhd_cpu"),
    ],
)
def test_jeans_wave_selfgravity_converges(input_file, basename):
    try:
        _run_selfgravity(input_file, basename)
    finally:
        _remove_binary_outputs(basename)
        testutils.cleanup()


def test_gravity_potential_binary_output():
    basename = "selfgravity_phi_cpu"
    try:
        _run_selfgravity("inputs/tests/selfgravity.athinput", basename)
        bin_outputs = list(Path("bin").glob(f"{basename}*.bin"))
        assert bin_outputs, "grav_phi binary output was not written"
        assert all(path.stat().st_size > 0 for path in bin_outputs)
    finally:
        _remove_binary_outputs(basename)
        testutils.cleanup()
