"""Focused regression for passive-CGL signal-speed, flux, and timestep isolation."""

import math
from pathlib import Path
import re
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[3]
NX1 = 16
CFL = 0.3
ISO_CS = 0.7
BX = math.sqrt(0.65)
BY = math.sqrt(0.35)
PRESSURE_STATES = ((1.0, 0.5), (1.5, 1.0))
FLOW_FIELDS = ("dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3")


def _run(command, cwd=None):
    result = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _input_text(
    ppar,
    pperp,
    nlim,
    passive=True,
    guide_b2=0.25,
    guide_b3=-0.15,
    field_amp=0.2,
):
    return f"""\
<job>
basename = passive_fast_path

<mesh>
nghost = 2
nx1 = {NX1}
x1min = 0.0
x1max = 1.0
ix1_bc = periodic
ox1_bc = periodic
nx2 = 1
x2min = 0.0
x2max = 1.0
ix2_bc = periodic
ox2_bc = periodic
nx3 = 1
x3min = 0.0
x3max = 1.0
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = {NX1}
nx2 = 1
nx3 = 1

<time>
evolution = dynamic
integrator = rk2
cfl_number = {CFL}
nlim = {nlim}
tlim = 10.0
ndiag = 1

<mhd>
eos = cgl
passive = {str(passive).lower()}
iso_sound_speed = {ISO_CS}
reconstruct = plm
rsolver = hlle
gamma = 1.666666666666667
dfloor = 1.0e-12
pfloor = 1.0e-12
tfloor = 1.0e-12
sfloor = 1.0e-12
bfloor = 1.0e-10

<problem>
pgen_name = divb_amr
rho0 = 1.0
ppar0 = {ppar}
pperp0 = {pperp}
vx0 = 0.16
vy0 = -0.04
vz0 = 0.02
guide_b1 = {BX}
guide_b2 = {guide_b2}
guide_b3 = {guide_b3}
field_amp = {field_amp}
field_k = 1.0
refine_levels = 0

<output1>
file_type = tab
variable = mhd_w_bcc
data_format = %24.16e
dcycle = 1
slice_x2 = 0.0
slice_x3 = 0.0
"""


def _diagnostic_timesteps(stdout):
    matches = re.findall(
        r"elapsed=\S+\s+cycle=(\d+)\s+time=(\S+)\s+dt=(\S+)", stdout
    )
    assert matches, stdout
    by_cycle = {}
    for cycle, _time, timestep in matches:
        by_cycle.setdefault(int(cycle), float(timestep))
    return by_cycle


def _read_final_table(run_dir):
    paths = sorted((run_dir / "tab").glob("passive_fast_path.mhd_w_bcc.*.tab"))
    assert paths, f"no table output in {run_dir}"
    lines = paths[-1].read_text().splitlines()
    labels = lines[1].lstrip("#").split()
    columns = {label: [] for label in labels}
    for line in lines[2:]:
        values = line.split()
        assert len(values) == len(labels), line
        for label, value in zip(labels, values):
            columns[label].append(float(value))
    return columns


def _isothermal_fast_speed(density, bx, by, bz):
    asq = ISO_CS * ISO_CS * density
    bperp2 = by * by + bz * bz
    qsq = bx * bx + bperp2 + asq
    tmp = bx * bx + bperp2 - asq
    return math.sqrt(0.5 * (qsq + math.sqrt(tmp * tmp + 4.0 * asq * bperp2))
                     / density)


def _active_cgl_fast_speed(density, ppar, pperp, bx, by, bz):
    bx2 = bx * bx
    b2 = bx2 + by * by + bz * bz
    mu2 = bx2 / b2
    qsq = b2 + 2.0 * pperp + (2.0 * ppar - pperp) * mu2
    disc = (
        qsq * qsq
        + 4.0 * pperp * pperp * (1.0 - mu2) * mu2
        - 12.0 * ppar * pperp * mu2 * (2.0 - mu2)
        + 12.0 * ppar * ppar * mu2 * mu2
        - 12.0 * bx2 * ppar
    )
    return math.sqrt(0.5 * (qsq + math.sqrt(abs(disc))) / density)


def _legacy_active_cgl_fast_speed(density, ppar, pperp, bx, by, bz):
    bx2 = bx * bx
    b2 = bx2 + by * by + bz * bz
    mu2 = bx2 / b2
    qsq = b2 + 2.0 * pperp + (2.0 * ppar - pperp) * mu2
    disc = (
        qsq * qsq
        + 4.0 * pperp * pperp * (1.0 - mu2) * mu2
        - 12.0 * ppar * pperp * mu2 * (2.0 - mu2)
        + 12.0 * ppar * pperp * mu2 * mu2
        - 12.0 * bx2 * ppar
    )
    return math.sqrt(0.5 * (qsq + math.sqrt(abs(disc))) / density)


def _initial_timestep(table, speed):
    limits = []
    dx = table["x1v"][1] - table["x1v"][0]
    for density, vx, bx, by, bz in zip(
        table["dens"], table["velx"], table["bcc1"], table["bcc2"], table["bcc3"]
    ):
        limits.append(dx / (abs(vx) + speed(density, bx, by, bz)))
    return CFL * min(limits)


def _assert_hard_bound_admissible(table, ppar, pperp):
    for bx, by, bz in zip(table["bcc1"], table["bcc2"], table["bcc3"]):
        b_squared = bx * bx + by * by + bz * bz
        anisotropy = pperp - ppar
        assert -1.5 * b_squared < anisotropy < b_squared


def test_passive_cgl_fast_overload_isolated_from_flow_and_timestep(tmp_path):
    unit_build_dir = tmp_path / "unit-build"
    _run([
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(unit_build_dir),
        "-G",
        "Unix Makefiles",
        "-DPROBLEM=unit_tests/cgl_passive_fast_path_test",
        "-DAthena_ENABLE_MPI=OFF",
        "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
    ])
    object_target = (
        "src/CMakeFiles/athena.dir/pgen/unit_tests/"
        "cgl_passive_fast_path_test.cpp.o"
    )
    _run(
        ["make", "-f", "src/CMakeFiles/athena.dir/build.make", object_target, "-j4"],
        cwd=unit_build_dir,
    )
    main_source = tmp_path / "main.cpp"
    main_source.write_text(
        "void RunCglPassiveFastPathChecks();\n"
        "int main() { RunCglPassiveFastPathChecks(); }\n"
    )
    compiler = None
    for line in (unit_build_dir / "CMakeCache.txt").read_text().splitlines():
        if line.startswith("CMAKE_CXX_COMPILER:"):
            compiler = line.split("=", 1)[1]
            break
    assert compiler is not None
    unit_executable = tmp_path / "cgl_passive_fast_path_test"
    _run([
        compiler,
        str(main_source),
        str(unit_build_dir / object_target),
        "-std=c++17",
        "-o",
        str(unit_executable),
    ])
    unit_result = _run([str(unit_executable)])
    assert "Passive CGL signal-speed and flux checks passed" in unit_result.stdout

    build_dir = tmp_path / "integration-build"
    _run([
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(build_dir),
        "-DPROBLEM=built_in_pgens",
        "-DAthena_ENABLE_MPI=OFF",
        "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
    ])
    _run(["cmake", "--build", str(build_dir), "--target", "athena", "-j4"])
    executable = build_dir / "src" / "athena"

    initial_dir = tmp_path / "initial"
    initial_dir.mkdir()
    initial_input = initial_dir / "athinput.passive_fast_path"
    initial_input.write_text(_input_text(*PRESSURE_STATES[0], nlim=0))
    initial_result = _run([str(executable), "-i", str(initial_input)], cwd=initial_dir)
    initial_table = _read_final_table(initial_dir)
    initial_timesteps = _diagnostic_timesteps(initial_result.stdout)
    for pressure_state in PRESSURE_STATES:
        _assert_hard_bound_admissible(initial_table, *pressure_state)

    results = []
    for index, (ppar, pperp) in enumerate(PRESSURE_STATES):
        run_dir = tmp_path / f"run-{index}"
        run_dir.mkdir()
        input_path = run_dir / "athinput.passive_fast_path"
        input_path.write_text(_input_text(ppar, pperp, nlim=2))
        result = _run([str(executable), "-i", str(input_path)], cwd=run_dir)
        results.append((_diagnostic_timesteps(result.stdout), _read_final_table(run_dir)))

    timesteps_a, table_a = results[0]
    timesteps_b, table_b = results[1]
    assert timesteps_a == timesteps_b

    expected_passive_dt = _initial_timestep(initial_table, _isothermal_fast_speed)
    assert math.isclose(initial_timesteps[0], expected_passive_dt, rel_tol=2.0e-6)
    assert math.isclose(timesteps_a[0], expected_passive_dt, rel_tol=2.0e-6)

    active_dts = []
    for ppar, pperp in PRESSURE_STATES:
        active_dts.append(_initial_timestep(
            initial_table,
            lambda density, bx, by, bz: _active_cgl_fast_speed(
                density, ppar, pperp, bx, by, bz
            )
        ))
    assert not math.isclose(active_dts[0], active_dts[1], rel_tol=0.05)
    assert all(
        not math.isclose(timesteps_a[0], active_dt, rel_tol=0.02)
        for active_dt in active_dts
    )

    active_dir = tmp_path / "active-oblique"
    active_dir.mkdir()
    active_input = active_dir / "athinput.active_fast_path"
    active_input.write_text(_input_text(
        *PRESSURE_STATES[0],
        nlim=0,
        passive=False,
        guide_b2=BY,
        guide_b3=0.0,
        field_amp=0.0,
    ))
    active_result = _run([str(executable), "-i", str(active_input)], cwd=active_dir)
    active_table = _read_final_table(active_dir)
    active_initial_dt = _diagnostic_timesteps(active_result.stdout)[0]
    expected_active_dt = _initial_timestep(
        active_table,
        lambda density, bx, by, bz: _active_cgl_fast_speed(
            density, PRESSURE_STATES[0][0], PRESSURE_STATES[0][1], bx, by, bz
        ),
    )
    legacy_active_dt = _initial_timestep(
        active_table,
        lambda density, bx, by, bz: _legacy_active_cgl_fast_speed(
            density, PRESSURE_STATES[0][0], PRESSURE_STATES[0][1], bx, by, bz
        ),
    )
    assert math.isclose(active_initial_dt, expected_active_dt, rel_tol=2.0e-6)
    assert not math.isclose(active_initial_dt, legacy_active_dt, rel_tol=5.0e-3)

    for field in FLOW_FIELDS:
        assert table_a[field] == table_b[field], field
    assert table_a["eint"] != table_b["eint"]
    assert table_a["p_perp"] != table_b["p_perp"]
