"""Replay maps use the production LF bound without changing the evolved state."""

import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


def executable():
    return Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()


def source(*, uniform=False, nu=0, forcing=False, nlim=3):
    text = passive_acceptance.input_text(
        passive=False, active_cgl=True, nx=16, mode=1 if uniform else 0,
        amp=0 if uniform else 0.1, thermal=0 if uniform else 0.2,
        nlim=nlim, tlim=1, forcing=forcing)
    return text + f"""
<time>
sts_integrator = rkl2
sts_safety = 0.3
<mhd>
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
cgl_lf_diagnostics = none
cgl_lf_strict_admissibility = true
lf_k_parallel = 6.283185307179586
lf_coefficient_mode = background
lf_c_parallel0 = 1
nu_coll = {nu}
<problem>
bx = 1
by = 0
bz = 0
<output2>
single_file_per_rank = false
"""


def run(directory, text=None, *, flags=(), restart=None, ranks=1, error=None):
    directory.mkdir()
    command = [str(executable()), "-d", str(directory)]
    if text is not None:
        deck = directory / "input.athinput"
        deck.write_text(text)
        command += ["-i", str(deck)]
    if restart is not None:
        command += ["-r", str(restart)]
    if ranks > 1:
        command = ["mpirun", "-np", str(ranks), *command]
    result = subprocess.run(command + list(flags), cwd=directory,
                            capture_output=True, text=True, timeout=120)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    if error is not None:
        assert result.returncode != 0 and error in result.stdout + result.stderr
        return None
    assert result.returncode == 0, result.stdout + result.stderr
    checkpoint = sorted((directory / "rst").glob("*.rst"))[-1]
    state = passive_acceptance.decode_restart(checkpoint, False, active_cgl=True)
    state["checkpoint"] = checkpoint
    state["history"] = passive_acceptance.history(
        directory / "passive_validation.mhd.hst")
    return state


def read_map(path):
    with path.open("rb") as stream:
        assert stream.readline() == b"# AthenaK LF stiffness v1\n"
        metadata = dict(item.split("=", 1) for item in stream.readline().decode().split())
        blocks = [list(map(int, stream.readline().split()[1:]))
                  for _ in range(int(metadata["nblocks"]))]
        assert stream.readline() == b"data\n"
        dtype = ("<" if metadata["endian"] == "little" else ">")
        dtype += f"f{metadata['real_bytes']}"
        values = np.frombuffer(stream.read(), dtype=dtype)
    shape = tuple(int(metadata[name]) for name in ("nblocks", "nx3", "nx2", "nx1"))
    return metadata, blocks, values.reshape(shape)


def assert_same_state(first, second):
    assert first["time"] == second["time"]
    assert first["dt"] == second["dt"]
    np.testing.assert_array_equal(first["payload"], second["payload"])
    assert first["history"].keys() == second["history"].keys()
    for name in first["history"]:
        np.testing.assert_array_equal(first["history"][name], second["history"][name])


def require_mpi():
    config = subprocess.run([str(executable()), "-c"], text=True,
                            capture_output=True, check=True).stdout
    if not re.search(r"MPI parallelism:\s+ON", config) or not shutil.which("mpirun"):
        pytest.skip("requires an MPI executable and mpirun")


@pytest.mark.parametrize("nu", [0, 10])
def test_stiffness_map_matches_uniform_analytic_bound(tmp_path, nu):
    directory = tmp_path / "uniform"
    run(directory, source(uniform=True, nu=nu),
        flags=("mhd/cgl_lf_stiffness_cycle=1",))
    paths = sorted(directory.glob("lf-stiffness.*.bin"))
    assert {path.name for path in paths} == {
        "lf-stiffness.1.pre.rank0.bin", "lf-stiffness.1.post.rank0.bin"}
    chi = 8 / (np.sqrt(8*np.pi)*2*np.pi + (3*np.pi-8)*nu)
    expected_dt = 1 / (2*16**2*chi)
    for path in paths:
        metadata, blocks, values = read_map(path)
        assert values.size == 16
        assert len(blocks) == 4
        np.testing.assert_allclose(values, expected_dt, rtol=3e-15, atol=0)
        np.testing.assert_allclose(float(metadata["raw_dt_rank"]), expected_dt,
                                   rtol=3e-15, atol=0)
        np.testing.assert_allclose(float(metadata["dt_sts"]), 0.3*expected_dt,
                                   rtol=3e-15, atol=0)
        ratio = float(metadata["dt_sweep"]) / float(metadata["dt_sts"])
        stages = max(3, int(np.ceil(0.5*(-1 + np.sqrt(9 + 16*ratio)))))
        stages += stages % 2 == 0
        assert int(metadata["nstages"]) == stages
        faces = np.genfromtxt(path.with_suffix(".tsv"), names=True, comments="#",
                              delimiter="\t", skip_header=2)
        assert len(faces) == 8
        np.testing.assert_array_equal(faces["valid"], 1)
        np.testing.assert_array_equal(faces["qpar_ratio"], 0)
        np.testing.assert_array_equal(faces["qperp_ratio"], 0)
        np.testing.assert_array_equal(faces["nu_actual"], nu)
        np.testing.assert_array_equal(faces["Brow_parallel"], 0)
        np.testing.assert_array_equal(faces["Brow_perp"], 0)
        np.testing.assert_allclose(2/np.maximum(
            faces["row_parallel"].reshape(4, 2).sum(axis=1),
            faces["row_perp"].reshape(4, 2).sum(axis=1)), expected_dt,
            rtol=3e-15, atol=0)


def test_logarithmic_bound_and_hot_faces_remain_representable(tmp_path):
    deck = Path(__file__).resolve().parents[3] / "inputs/unit_tests/cgl_lf_uniform_timestep.athinput"
    subprocess.run([str(executable()), "-i", str(deck), "time/nlim=1",
                    "mhd/cgl_lf_stiffness_cycle=0", "mhd/lf_c_parallel0=1e-200",
                    "mhd/lf_k_parallel=1", "mhd/nu_coll=1", "mesh/x1max=1.6e-199",
                    "problem/pperp0=1.5"], cwd=tmp_path, check=True, capture_output=True)
    # dx = c_parallel = 1e-200. The collisional chi and dx squared both underflow.
    expected_dt = (3*np.pi - 8)/16
    for phase in ("pre", "post"):
        path = tmp_path / f"lf-stiffness.0.{phase}.rank0.bin"
        _, _, values = read_map(path)
        np.testing.assert_allclose(values, expected_dt, rtol=3e-13)
        faces = np.genfromtxt(path.with_suffix(".tsv"), names=True, skip_header=2)
        assert np.all(faces["valid"] == 2)
        np.testing.assert_allclose(2/faces["row_parallel"].sum(), expected_dt, rtol=3e-13)


@pytest.mark.parametrize("ranks", [1, 2])
def test_stiffness_map_preserves_forced_active_state(tmp_path, ranks):
    if ranks > 1:
        require_mpi()
    text = source(forcing=True)
    baseline = run(tmp_path / "off", text, ranks=ranks)
    assert not list((tmp_path / "off").glob("lf-stiffness.*"))
    mapped = run(tmp_path / "on", text, ranks=ranks,
                 flags=("mhd/cgl_lf_stiffness_cycle=1",))
    assert_same_state(baseline, mapped)
    assert len(list((tmp_path / "on").glob("lf-stiffness.*.bin"))) == 2*ranks
    for phase in ("pre", "post"):
        entries = [read_map(path) for path in sorted(
            (tmp_path / "on").glob(f"lf-stiffness.1.{phase}.rank*.bin"))]
        assert sum(values.size for _, _, values in entries) == 16**3
        gids = [block[0] for _, blocks, _ in entries for block in blocks]
        assert sorted(gids) == list(range(8))
        for metadata, _, values in entries:
            assert np.isfinite(values).all() and (values > 0).all()
            assert float(metadata["raw_dt_rank"]) == values.min()


def test_stiffness_map_is_rank_layout_independent(tmp_path):
    require_mpi()
    maps = []
    for ranks in (1, 2):
        directory = tmp_path / f"rank{ranks}"
        run(directory, source(nlim=1), ranks=ranks,
            flags=("mhd/cgl_lf_stiffness_cycle=0",))
        by_gid = {}
        for path in directory.glob("lf-stiffness.0.pre.rank*.bin"):
            _, blocks, values = read_map(path)
            by_gid.update((block[0], values[index]) for index, block in enumerate(blocks))
        maps.append(by_gid)
    assert set(maps[0]) == set(maps[1]) == set(range(4))
    assert np.ptp(np.concatenate([values.ravel() for values in maps[0].values()])) > 0
    for gid in maps[0]:
        np.testing.assert_array_equal(maps[0][gid], maps[1][gid])


def test_stiffness_face_evidence_includes_magnetic_gradient(tmp_path):
    root = Path(__file__).resolve().parents[3]
    text = (root / "inputs/unit_tests/cgl_lf_quant_grad_b.athinput").read_text()
    text += """
<output1>
file_type = hst
data_format = %25.17e
dcycle = 1
<output2>
file_type = rst
single_file_per_rank = false
dcycle = 1000
"""
    directory = tmp_path / "grad_b"
    # This fixture has uniform temperatures and a spatially varying magnetic field.
    run(directory, text, flags=("job/basename=passive_validation",
                                "mhd/cgl_lf_stiffness_cycle=0"))
    table = np.genfromtxt(directory / "lf-stiffness.0.pre.rank0.tsv", names=True,
                          comments="#", delimiter="\t", skip_header=2)
    assert np.all(table["valid"] == 1)
    assert np.any(table["Brow_perp"] > 0)
    assert np.any(np.abs(table["grad_b_parallel"]) > 0)
    assert np.any(table["qperp_ratio"] > 0)
    for component in ("parallel", "perp"):
        np.testing.assert_allclose(table[f"Trow_{component}"] +
                                   table[f"Brow_{component}"],
                                   table[f"row_{component}"], rtol=3e-15, atol=0)


def test_stiffness_faces_report_real_scattering_even_with_zero_normal_field(tmp_path):
    # A 2D state varying only along B has conductive transverse faces with zero rows.
    text = source(uniform=True, nu=2, nlim=1) + """
<mesh>
nx2 = 4
<meshblock>
nx2 = 4
<problem>
thermal_amplitude = 1
<mhd>
limiter_nu_coll = 100
"""
    maps = []
    for limiter in (False, True):
        directory = tmp_path / str(limiter)
        run(directory, text + f"\n<mhd>\nmirror_limiter = {str(limiter).lower()}\n",
            flags=("mhd/cgl_lf_stiffness_cycle=0",))
        _, blocks, values = read_map(directory / "lf-stiffness.0.pre.rank0.bin")
        maps.append(values)
        faces = np.genfromtxt(directory / "lf-stiffness.0.pre.rank0.tsv", names=True,
                              comments="#", delimiter="\t", skip_header=2)
        expected_nu = 2 + 100*limiter*(faces["pperp_face"]-faces["ppar_face"] >= 0.5)
        np.testing.assert_array_equal(faces["nu_actual"], expected_nu)
        transverse = faces[faces["dir"] == 1]
        np.testing.assert_array_equal(transverse["bn_ct"], 0)
        np.testing.assert_array_equal(transverse["row_parallel"], 0)
        assert np.any(transverse["nu_actual"] > 2) == limiter
        for face in transverse:
            location = next(block[2] for block in blocks if block[0] == face["gid"])
            x = (4*location + face["i"] + 0.5)/16
            temperatures = 2 - (2/3)*np.cos(2*np.pi*(x + np.array([-1, 0, 1])/16))
            left, right = np.diff(temperatures)
            slope = 2*left*right/(left+right) if left*right > 0 else 0
            gradient = 16*slope
            expected = abs(gradient/face["ppar_face"] /
                           (2*np.pi + (3*np.pi-8)*face["nu_actual"]/np.sqrt(8*np.pi)))
            np.testing.assert_allclose(face["qpar_ratio"], expected, rtol=2e-12,
                                       atol=2e-15)
        assert np.any(transverse["qpar_ratio"] > 0)
    # The controller intentionally excludes limiter scattering.
    np.testing.assert_array_equal(maps[0], maps[1])


def test_stiffness_map_can_be_enabled_on_older_active_restart(tmp_path):
    seed = run(tmp_path / "seed", source(nlim=1))
    old_checkpoint = tmp_path / "older.rst"
    raw = seed["checkpoint"].read_bytes()
    header, payload = raw.split(b"<par_end>\n", 1)
    header = re.sub(rb"(?m)^cgl_lf_stiffness_cycle\s*=.*\n", b"", header)
    old_checkpoint.write_bytes(header + b"<par_end>\n" + payload)
    baseline = run(tmp_path / "off", restart=old_checkpoint, flags=("time/nlim=3",))
    mapped = run(tmp_path / "on", restart=old_checkpoint,
                 flags=("time/nlim=3", "mhd/cgl_lf_stiffness_cycle=1"))
    assert_same_state(baseline, mapped)
    assert len(list((tmp_path / "on").glob("lf-stiffness.1.*.bin"))) == 2


@pytest.mark.parametrize("cycle", ["-2", "1.5", "nan"])
def test_stiffness_map_rejects_invalid_cycle(tmp_path, cycle):
    run(tmp_path / "invalid", source(), flags=(f"mhd/cgl_lf_stiffness_cycle={cycle}",),
        error="cgl_lf_stiffness_cycle must be")
