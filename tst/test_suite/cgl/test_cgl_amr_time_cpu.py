"""Physical-time uniform refinement crosses its boundary once and survives restart."""

import os
from pathlib import Path
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[3]
TEMPLATE = REPO / "inputs/tests/cgl_amr_primitive_current_churn_restart.athinput"
sys.path.insert(0, str(REPO / "vis/python"))
import bin_convert  # noqa: E402


def run(tmp_path, *args):
    tmp_path.mkdir(exist_ok=True)
    return subprocess.run(
        [str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()),
         *args, "-d", str(tmp_path)], cwd=tmp_path,
        capture_output=True, text=True, check=False,
    )


def input_file(tmp_path, start):
    template = TEMPLATE.read_text()
    path = tmp_path / "time.athinput"
    path.write_text(template.replace("method = user", f"method = time\ntime_start = {start}")
                    + "\n<amr_criterion2>\nmethod = slope\nvariable = mhd_w_d\n"
                    + "value_max = 1000\n")
    return path


def test_cgl_amr_time_boundary_and_restart(tmp_path):
    start = 0.0015
    deck = input_file(tmp_path, start)
    result = run(tmp_path, "-i", str(deck), "mesh/nx1=8", "mesh/nx2=8",
                 "mesh/nx3=8", "meshblock/nx3=8", "time/nlim=-1",
                 "mesh_refinement/max_nmb_per_rank=8",
                 f"time/tlim={start}", "output2/dcycle=1", "output3/dcycle=1")
    assert result.returncode == 0, result.stdout + result.stderr
    states = [bin_convert.read_binary(str(path))
              for path in sorted((tmp_path / "bin").glob("*.state.*.bin"))]
    assert all(state["n_mbs"] == 1 for state in states if state["time"] < start)
    assert states[-1]["time"] == start
    assert states[-1]["n_mbs"] == 8

    checkpoint = sorted((tmp_path / "rst").rglob("*.rst"))[-1]
    resumed = tmp_path / "resumed"
    result = run(resumed, "-r", str(checkpoint), f"time/tlim={2 * start}")
    assert result.returncode == 0, result.stdout + result.stderr
    states = [bin_convert.read_binary(str(path))
              for path in sorted((resumed / "bin").glob("*.state.*.bin"))]
    assert states and all(state["n_mbs"] == 8 for state in states)
    assert states[-1]["time"] == 2 * start


def test_cgl_amr_time_overlay_on_uniform_lf_checkpoint(tmp_path):
    template = TEMPLATE.read_text()
    template = template[:template.index("<mesh_refinement>")] + template[template.index("<time>"):]
    template = template.replace("integrator = rk2", "integrator = rk2\nsts_integrator = rkl2")
    template = template.replace("<mhd>", "<mhd>\ncgl_heat_flux = landau_fluid\n"
                                "cgl_heat_flux_integrator = sts\nlf_k_parallel = 32\n"
                                "cgl_lf_record_pressure_work = true\n")
    template = template.replace("dcycle = 3", "dt = 0.01")
    deck = tmp_path / "uniform.athinput"
    deck.write_text(template)
    result = run(tmp_path, "-i", str(deck), "mesh/nx1=8", "mesh/nx2=8",
                 "mesh/nx3=8", "meshblock/nx3=8", "time/nlim=1",
                 "output3/dcycle=1")
    assert result.returncode == 0, result.stdout + result.stderr
    checkpoint = sorted((tmp_path / "rst").rglob("*.rst"))[-1]
    with checkpoint.open("rb") as stream:
        header = stream.read(20000)
    assert b"cgl_amr_repair_restart_version" not in header
    assert b"cgl_lf_record_pressure_work" in header

    overlay = tmp_path / "amr.athinput"
    overlay.write_text("<mesh_refinement>\nrefinement = adaptive\nnum_levels = 2\n"
                       "ncycle_check = 1\nrefinement_interval = 1\n"
                       "prolong_primitives = true\nmax_nmb_per_rank = 8\n"
                       "<amr_criterion1>\nmethod = time\ntime_start = 0.01\n"
                       "<mhd>\ncgl_lf_record_pressure_work = false\n")
    refined = tmp_path / "refined"
    result = run(refined, "-r", str(checkpoint), "-i", str(overlay),
                 "time/nlim=-1", "time/tlim=0.01")
    assert result.returncode == 0, result.stdout + result.stderr
    states = [bin_convert.read_binary(str(path))
              for path in sorted((refined / "bin").glob("*.state.*.bin"))]
    assert all(state["n_mbs"] == 1 for state in states[:-1])
    assert any(state["time"] < 0.01 and state["cycle"] > 1 for state in states)
    assert states[-2]["time"] == states[-1]["time"] == 0.01
    final_cycle = states[-1]["cycle"]
    assert states[-2]["cycle"] == final_cycle
    assert states[-1]["n_mbs"] == 8
    assert len(list((refined / "rst").rglob("*.rst"))) == 1
    checkpoint = sorted((refined / "rst").rglob("*.rst"))[-1]
    with checkpoint.open("rb") as stream:
        assert b"cgl_amr_repair_restart_version" in stream.read(20000)
    resumed = tmp_path / "resumed"
    # Raising the level ceiling must not refine again until the later threshold.
    result = run(resumed, "-r", str(checkpoint), "time/tlim=0.02",
                 "amr_criterion1/time_start=0.02", "mesh_refinement/num_levels=3",
                 "mesh_refinement/max_nmb_per_rank=64")
    assert result.returncode == 0, result.stdout + result.stderr
    states = [bin_convert.read_binary(str(path))
              for path in sorted((resumed / "bin").glob("*.state.*.bin"))]
    assert all(state["n_mbs"] == 8 for state in states[:-1])
    assert any(0.01 < state["time"] < 0.02 for state in states)
    assert states[-2]["time"] == states[-1]["time"] == 0.02
    assert states[-2]["cycle"] == states[-1]["cycle"]
    assert states[-1]["n_mbs"] == 64
    assert len(list((resumed / "rst").rglob("*.rst"))) == 1
    checkpoint = sorted((resumed / "rst").rglob("*.rst"))[-1]
    final_cycle = states[-1]["cycle"]
    held = tmp_path / "held"
    result = run(held, "-r", str(checkpoint), f"time/nlim={final_cycle + 1}",
                 "time/tlim=0.03")
    assert result.returncode == 0, result.stdout + result.stderr
    state = bin_convert.read_binary(str(sorted((held / "bin").glob("*.state.*.bin"))[-1]))
    assert state["cycle"] == final_cycle + 1 and state["n_mbs"] == 64
    assert state["time"] > 0.02


@pytest.mark.parametrize("start", ["-1", "nan", "inf"])
def test_cgl_amr_time_rejects_invalid_start(tmp_path, start):
    result = run(tmp_path, "-i", str(input_file(tmp_path, start)), "time/nlim=0")
    assert result.returncode != 0
    assert "AMR time_start must be finite and nonnegative" in result.stdout + result.stderr
