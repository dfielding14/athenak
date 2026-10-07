"""Production passive wall decoding must agree across independent kernels."""

import os
from pathlib import Path
import shlex
import subprocess


REPOSITORY = Path(__file__).resolve().parents[3]


def test_passive_wall_survives_canonical_c2p(tmp_path):
    supplied = os.environ.get("ATHENAK_CGL_PASSIVE_WALL_BINARY")
    if supplied:
        executable = Path(supplied).resolve()
        assert executable.is_file()
    else:
        build = tmp_path / "build"
        commands = [
            ["cmake", "-S", str(REPOSITORY), "-B", str(build),
             "-DPROBLEM=unit_tests/cgl_passive_wall_decode_test", "-DAthena_ENABLE_MPI=OFF",
             "-DKokkos_ENABLE_HIP=OFF", "-DKokkos_ENABLE_SERIAL=ON", "-DCMAKE_BUILD_TYPE=Release"],
            ["cmake", "--build", str(build), "--target", "athena", "-j4"],
        ]
        for index, command in enumerate(commands):
            result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
            (tmp_path / f"build-{index}.log").write_text(result.stdout+result.stderr)
            assert result.returncode == 0, result.stdout+result.stderr
        executable = build / "src/athena"
    input_file = tmp_path / "passive-wall.athinput"
    input_file.write_text("""<job>
basename = passive_wall
<mesh>
nghost = 3
nx1 = 8
nx2 = 8
nx3 = 8
x1min = 0
x1max = 1
x2min = 0
x2max = 1
x3min = 0
x3max = 1
ix1_bc = periodic
ox1_bc = periodic
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = periodic
ox3_bc = periodic
<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8
<time>
evolution = dynamic
integrator = rk2
nlim = 0
tlim = 0.001
cfl_number = 0.3
<mhd>
eos = cgl
passive = true
iso_sound_speed = 2.23606797749979
reconstruct = ppm4
rsolver = hlle
fofc = false
mirror_limiter = true
firehose_limiter = true
mirror_threshold = 1
firehose_threshold = 2
backup_limiters = false
limiter_nu_coll = 1e10
nu_coll = 0
dfloor = 1e-12
pfloor = 1e-12
bfloor = 1e-10
<problem>
""")
    launcher = shlex.split(os.environ.get("ATHENAK_CGL_PASSIVE_WALL_LAUNCHER", ""))
    command = launcher + [str(executable), "-i", str(input_file)]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    (tmp_path / "passive-wall.log").write_text(result.stdout+result.stderr)
    assert result.returncode == 0, result.stdout+result.stderr
    assert "Passive production wall decode checks passed" in result.stdout
    old = os.environ.get("ATHENAK_CGL_PASSIVE_WALL_OLD_BINARY")
    if old:
        control = tmp_path / "unpatched-negative-control"
        control.mkdir()
        result = subprocess.run(launcher + [old, "-i", str(input_file)], cwd=control,
                                capture_output=True, text=True)
        (control / "passive-wall.log").write_text(result.stdout+result.stderr)
        assert result.returncode != 0, "unpatched independent-kernel decode unexpectedly passed"
        assert any("Passive wall decode regression failed: "+label in result.stdout
                   for label in ("captured wall survives canonical C2P",
                                 "encoded walls survive independent canonical C2P"))
