"""One-time uniform restart refinement preserves physics and restart state."""

import os
from pathlib import Path
import re
import shutil
import struct
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance
from test_suite.cgl import test_cgl_lf_cooling


MODES = ("ideal", "isothermal", "cgl", "passive", "adiabatic")


def input_text(mode, *, blocks=1, nlim=0, nghost=4, shape=None,
               block_shape=(8, 8, 8)):
    shape = shape or (8*blocks, 8, 8)
    source = passive_acceptance.input_text(
        passive=mode in ("passive", "adiabatic"),
        adiabatic=mode == "adiabatic", ideal_mhd=mode == "ideal",
        active_cgl=mode == "cgl", nx=8, nlim=nlim, tlim=10,
        amp=0, thermal=0,
    )
    return source + f"""
<mesh>
nghost = {nghost}
nx1 = {shape[0]}
nx2 = {shape[1]}
nx3 = {shape[2]}
x1max = {shape[0]/8}
x2max = {shape[1]/8}
x3max = {shape[2]/8}
<meshblock>
nx1 = {block_shape[0]}
nx2 = {block_shape[1]}
nx3 = {block_shape[2]}
<mhd>
gamma = 1.6666666666666667
tfloor = 0
sfloor = 0
<problem>
p_gas0 = 2
<output2>
single_file_per_rank = false
dcycle = 1
"""


def mpi_available():
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    config = subprocess.run([str(executable), "-c"], capture_output=True,
                            text=True, check=True).stdout
    return re.search(r"MPI parallelism:\s+ON", config) and shutil.which("mpirun")


def run(directory, mode, *, source=None, restart=None, refine=False, reblock=False,
        flags=(), ranks=1, error=False):
    directory.mkdir()
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    command = [str(executable), "-d", str(directory)]
    if source is not None:
        deck = directory / "input.athinput"
        deck.write_text(source)
        command += ["-i", str(deck)]
    if restart is not None:
        command += ["-r", str(restart)]
    if refine:
        command += ["mesh/refine_on_restart=true"]
    if reblock:
        command += ["mesh/reblock_on_restart=true"]
    if ranks > 1:
        if not mpi_available():
            pytest.skip("requires an MPI executable and mpirun")
        command = ["mpirun", "-np", str(ranks), *command]
    result = subprocess.run(command + list(flags), cwd=directory,
                            capture_output=True, text=True, check=False, timeout=120)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output)
    if error:
        assert result.returncode != 0, output
        return output
    assert result.returncode == 0, output
    checkpoint = sorted((directory / "rst").glob("*.rst"))[-1]
    state = passive_acceptance.decode_restart(
        checkpoint, mode in ("passive", "adiabatic"),
        active_cgl=mode == "cgl", adiabatic=mode == "adiabatic",
        ideal_mhd=mode == "ideal")
    raw = checkpoint.read_bytes()
    header = raw.index(b"<par_end>\n") + len(b"<par_end>\n")
    nmb = len(state["u"])
    box = np.frombuffer(raw, dtype="=f8", count=9, offset=header + 8)
    state.update(path=checkpoint, text=raw[:header].decode(),
                 block_shape=tuple(reversed(state["active"].shape[-3:])),
                 shape=struct.unpack_from("=3i", raw, header + 8 + 72 + 4),
                 lengths=tuple(box[3:6] - box[:3]),
                 locations=np.frombuffer(raw, dtype="=i4", count=4*nmb,
                                         offset=header + 252).reshape(nmb, 4),
                 offset=len(raw) - state["payload"].nbytes,
                 metadata=raw[header + 252 + 20*nmb:
                              len(raw) - state["payload"].nbytes - 8])
    return state


def global_fields(state):
    ng = state["ng"]
    n1, n2, n3 = state["shape"]
    bx, by, bz = state["block_shape"]
    u = np.empty((state["u"].shape[1], n3, n2, n1))
    faces = [np.full((n3, n2, n1+1), np.nan),
             np.full((n3, n2+1, n1), np.nan),
             np.full((n3+1, n2, n1), np.nan)]
    for m, (x, y, z, _) in enumerate(state["locations"]):
        i, j, k = bx*x, by*y, bz*z
        u[:, k:k+bz, j:j+by, i:i+bx] = state["active"][m]
        for axis, field in enumerate(state["faces"]):
            sizes = [bz, by, bx]
            sizes[2-axis] += 1
            nk, nj, ni = sizes
            target = faces[axis][k:k+nk, j:j+nj, i:i+ni]
            values = field[m, ng:ng+nk, ng:ng+nj, ng:ng+ni]
            shared = np.isfinite(target)
            np.testing.assert_array_equal(target[shared], values[shared])
            target[:] = values
    assert np.isfinite(u).all() and all(np.isfinite(b).all() for b in faces)
    return u, faces


def curl_field(state):
    # Curl edge-centered A with the same differences used by discrete face divergence.
    nx, ny, nz = state["shape"]
    dx, dy, dz = np.array(state["lengths"])/state["shape"]
    xc, yc, zc = [(np.arange(n) + 0.5)/n for n in (nx, ny, nz)]
    xf, yf, zf = [np.arange(n+1) % n/n for n in (nx, ny, nz)]
    k = 2*np.pi
    az = 0.025*np.cos(k*zc[:, None, None])*np.sin(k*yf[None, :, None])*np.sin(k*xf)
    ax = 0.017*np.sin(k*zf[:, None, None])*np.sin(k*yf[None, :, None])*np.cos(k*xc)
    ay = 0.011*np.sin(k*zf[:, None, None])*np.cos(k*yc[None, :, None])*np.sin(k*xf)
    return [1 + np.diff(az, axis=1)/dy - np.diff(ay, axis=0)/dz,
            np.diff(ax, axis=0)/dz - np.diff(az, axis=2)/dx,
            1 + np.diff(ay, axis=2)/dx - np.diff(ax, axis=1)/dy]


def patch_checkpoint(state, mode, path, *, nonuniform=False, poison_ghosts=False):
    ng = state["ng"]
    nx, ny, nz = state["block_shape"]
    active_cells = (slice(ng, ng+nz), slice(ng, ng+ny), slice(ng, ng+nx))
    raw = bytearray(state["path"].read_bytes())
    payload = np.frombuffer(raw, dtype="=f8", offset=state["offset"])
    payload = payload.reshape(state["payload"].shape)
    u = payload[:, :state["u"][0].size].reshape(state["u"].shape)
    faces, offset = [], state["u"][0].size
    for field in state["faces"]:
        faces.append(payload[:, offset:offset+field[0].size].reshape(field.shape))
        offset += field[0].size
    if nonuniform:
        magnetic = curl_field(state)
    for m, (lx, ly, lz, _) in enumerate(state["locations"]):
        if nonuniform:
            x = (nx*lx + np.arange(nx) + 0.5)/state["shape"][0]
            y = (ny*ly + np.arange(ny) + 0.5)/state["shape"][1]
            z = (nz*lz + np.arange(nz) + 0.5)/state["shape"][2]
            zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
            for axis in range(3):
                sizes = [nz, ny, nx]
                sizes[2-axis] += 1
                nk, nj, ni = sizes
                faces[axis][m, ng:ng+nk, ng:ng+nj, ng:ng+ni] = magnetic[axis][
                    nz*lz:nz*lz+nk, ny*ly:ny*ly+nj, nx*lx:nx*lx+ni]
            bx = 0.5*(faces[0][(m, *active_cells)] +
                      faces[0][m, ng:ng+nz, ng:ng+ny, ng+1:ng+nx+1])
            by = 0.5*(faces[1][(m, *active_cells)] +
                      faces[1][m, ng:ng+nz, ng+1:ng+ny+1, ng:ng+nx])
            bz = 0.5*(faces[2][(m, *active_cells)] +
                      faces[2][m, ng+1:ng+nz+1, ng:ng+ny, ng:ng+nx])
            rho = 1 + 0.1*np.cos(2*np.pi*xx)
            velocity = np.array([0.1*np.sin(2*np.pi*yy),
                                 0.1*np.cos(2*np.pi*zz),
                                 0.1*np.sin(2*np.pi*xx)])
            pp, pt = 2 + 0.1*np.sin(2*np.pi*yy), 2 + 0.1*np.cos(2*np.pi*zz)
            energy = 0.5*rho*np.sum(velocity**2, axis=0) + 0.5*(bx*bx + by*by + bz*bz)
            active = u[(m, slice(None), *active_cells)]
            active[0], active[1:4] = rho, rho*velocity
            if mode == "ideal":
                active[4] = energy + 3
            elif mode != "isothermal":
                bm = np.sqrt(bx*bx + by*by + bz*bz)
                active[5] = rho*(np.log(pt/pp) + 2*np.log(rho) - 3*np.log(bm))
                active[4] = (energy + 0.5*pp + pt if mode == "cgl" else
                             rho*(np.log(pp) + 2*np.log(bm) - 3*np.log(rho)))
                if mode == "adiabatic":
                    active[6] = energy + 3
        if poison_ghosts:
            active = np.zeros(u[m].shape[1:], dtype=bool)
            active[active_cells] = True
            u[m, :, ~active] = np.nan
            for axis, field in enumerate(faces):
                active = np.zeros(field[m].shape, dtype=bool)
                sizes = [nz, ny, nx]
                sizes[2-axis] += 1
                nk, nj, ni = sizes
                active[ng:ng+nk, ng:ng+nj, ng:ng+ni] = True
                field[m, ~active] = np.nan
    path.write_bytes(raw)
    return path


@pytest.mark.parametrize("mode", MODES)
def test_refine_restart_constant_state_on_more_ranks(tmp_path, mode):
    initial = run(tmp_path / "initial", mode, source=input_text(mode))
    refined = run(tmp_path / "fine", mode, restart=initial["path"], refine=True, ranks=8)
    assert refined["shape"] == (16, 16, 16) and len(refined["u"]) == 8
    assert refined["u"].shape[2:] == initial["u"].shape[2:]
    assert (refined["time"], refined["cycle"]) == (initial["time"], initial["cycle"])
    assert refined["dt"] < initial["dt"]
    assert re.search(r"^refine_on_restart\s*=\s*(?:false|0)\s",
                     refined["text"], re.MULTILINE)
    u, faces = global_fields(refined)
    np.testing.assert_allclose(
        u, np.broadcast_to(initial["active"][0, :, :1, :1, :1], u.shape),
        rtol=2e-14, atol=2e-14)
    for first, second in zip(faces, initial["faces"]):
        ng = initial["ng"]
        np.testing.assert_allclose(first, second[0, ng, ng, ng], rtol=0, atol=2e-15)
    held = run(tmp_path / "held", mode, restart=refined["path"], ranks=2)
    np.testing.assert_array_equal(held["payload"], refined["payload"])
    second = run(tmp_path / "second", mode, restart=held["path"], refine=True)
    assert second["shape"] == (32, 32, 32) and len(second["u"]) == 64
    assert second["time"] == initial["time"] and second["cycle"] == initial["cycle"]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("nghost", [2, 4])
def test_refine_restart_nonuniform_conservation_ct_and_ghosts(tmp_path, mode, nghost):
    initial = run(tmp_path / "initial", mode,
                  source=input_text(mode, blocks=2, nghost=nghost),
                  flags=("mhd/fofc=false",))
    clean = patch_checkpoint(initial, mode, tmp_path / "nonuniform.rst", nonuniform=True)
    dirty = patch_checkpoint(initial, mode, tmp_path / "poisoned.rst",
                             nonuniform=True, poison_ghosts=True)
    refined = run(tmp_path / "fine", mode, restart=clean, refine=True)
    poisoned = run(tmp_path / "poisoned", mode, restart=dirty, refine=True, ranks=8)
    np.testing.assert_array_equal(poisoned["payload"], refined["payload"])
    # The source view still carries the initial payload, so decode its patched arrays.
    source = passive_acceptance.decode_restart(
        clean, mode in ("passive", "adiabatic"), active_cgl=mode == "cgl",
        adiabatic=mode == "adiabatic", ideal_mhd=mode == "ideal")
    variables = [0, 1, 2, 3] if mode == "isothermal" else [0, 1, 2, 3, 4]
    if mode == "cgl":
        assert np.isfinite(refined["active"][:, 5]).all()
    elif mode in ("passive", "adiabatic"):
        variables.append(5)
        if mode == "adiabatic":
            variables.append(6)
        assert np.min(refined["pp"]) > 0 and np.min(refined["pt"]) > 0
    np.testing.assert_allclose(
        np.mean(refined["active"][:, variables], axis=(0, 2, 3, 4)),
        np.mean(source["active"][:, variables], axis=(0, 2, 3, 4)),
        rtol=3e-13, atol=3e-13)
    _, b = global_fields(refined)
    n1, n2, n3 = refined["shape"]
    lx, ly, lz = refined["lengths"]
    divb = (n1/lx*np.diff(b[0], axis=2) + n2/ly*np.diff(b[1], axis=1) +
            n3/lz*np.diff(b[2], axis=0))
    assert np.max(np.abs(divb)) < 2e-13
    assert np.max(np.abs(np.diff(b[0], axis=2))) > 1e-3
    np.testing.assert_allclose([np.mean(field) for field in b], [1, 0, 1], atol=2e-15)


def test_refine_restart_passive_gas_matches_ideal_mhd(tmp_path):
    states = {}
    for mode in ("ideal", "adiabatic"):
        initial = run(tmp_path / f"initial_{mode}", mode,
                      source=input_text(mode, blocks=2))
        source = patch_checkpoint(initial, mode, tmp_path / f"{mode}.rst",
                                  nonuniform=True)
        states[mode] = run(tmp_path / f"fine_{mode}", mode, restart=source, refine=True)
    ideal, passive = states["ideal"], states["adiabatic"]
    np.testing.assert_array_equal(passive["active"][:, :4], ideal["active"][:, :4])
    np.testing.assert_array_equal(passive["active"][:, 6], ideal["active"][:, 4])
    for first, second in zip(passive["faces"], ideal["faces"]):
        np.testing.assert_array_equal(first, second)


@pytest.mark.parametrize("mode", ["cgl", "passive", "adiabatic"])
def test_refine_restart_preserves_pressure_floors(tmp_path, mode):
    floor = 1e-8
    initial = run(tmp_path / "initial", mode, source=input_text(mode),
                  flags=(f"mhd/pfloor={floor}", f"problem/p0={floor}",
                         f"problem/p_gas0={floor}"))
    refined = run(tmp_path / "fine", mode, restart=initial["path"], refine=True)
    u, b = refined["active"], refined["bcc"]
    if mode == "cgl":
        rho = u[:, 0]
        bsqr = np.sum(b*b, axis=1)
        thermal = u[:, 4] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/rho - 0.5*bsqr
        ratio = np.exp(u[:, 5]/rho - 2*np.log(rho) + 1.5*np.log(bsqr))
        pp = thermal/(0.5 + ratio)
        pt = ratio*pp
    else:
        pp, pt = refined["pp"], refined["pt"]
    assert np.min(pp) >= floor*(1 - 1e-10)
    assert np.min(pt) >= floor*(1 - 1e-10)
    if mode == "adiabatic":
        gas = (2/3)*(u[:, 6] - 0.5*np.sum(u[:, 1:4]**2, axis=1)/u[:, 0] -
                     0.5*np.sum(b*b, axis=1))
        assert np.min(gas) >= floor*(1 - 1e-8)


@pytest.mark.parametrize("mode", ["ideal", "cgl"])
def test_refine_restart_matches_native_uniform_amr_transfer(tmp_path, mode):
    initial = run(tmp_path / "initial", mode, source=input_text(mode, blocks=2))
    checkpoint = patch_checkpoint(initial, mode, tmp_path / "source.rst", nonuniform=True)
    coarse = run(tmp_path / "coarse", mode, restart=checkpoint, flags=("time/nlim=1",))
    native = run(tmp_path / "native", mode, restart=checkpoint, flags=("time/nlim=1",),
                 source="""
<mesh_refinement>
refinement = adaptive
num_levels = 2
ncycle_check = 1
refinement_interval = 1
prolong_primitives = true
max_nmb_per_rank = 16
<amr_criterion1>
method = time
time_start = 0
""")
    refined = run(tmp_path / "fine", mode, restart=coarse["path"], refine=True)
    assert native["time"] == refined["time"] and native["cycle"] == refined["cycle"]
    # Native AMR retains the coarse root dimensions; all its leaves are uniformly fine.
    native["shape"] = refined["shape"]
    native_u, native_b = global_fields(native)
    fine_u, fine_b = global_fields(refined)
    np.testing.assert_allclose(fine_u, native_u, rtol=5e-13, atol=5e-13)
    for first, second in zip(fine_b, native_b):
        np.testing.assert_array_equal(first, second)


def test_refine_restart_requires_checkpoint(tmp_path):
    message = run(tmp_path / "invalid", "ideal", source=input_text("ideal"),
                  refine=True, error=True)
    assert "refine_on_restart requires a shared checkpoint" in message


@pytest.mark.parametrize("flag,message", [
    ("meshblock/nx1=4", "unchanged MeshBlock dimensions"),
    ("mesh/nx1=16", "leave mesh dimensions unchanged"),
    ("mesh/x1max=2", "unchanged physical domain"),
])
def test_refine_restart_rejects_incompatible_geometry(tmp_path, flag, message):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal"))
    output = run(tmp_path / "invalid", "ideal", restart=initial["path"],
                 refine=True, flags=(flag,), error=True)
    assert message in output


def test_refine_restart_rejects_rank_local_checkpoint(tmp_path):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal"))
    directory = tmp_path / "rank_00000000"
    directory.mkdir()
    checkpoint = directory / initial["path"].name
    shutil.copyfile(initial["path"], checkpoint)
    output = run(tmp_path / "invalid", "ideal", restart=checkpoint,
                 refine=True, error=True)
    assert "refine_on_restart requires a shared checkpoint" in output


def test_refine_restart_requires_enough_final_blocks(tmp_path):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal"))
    output = run(tmp_path / "invalid", "ideal", restart=initial["path"],
                 refine=True, ranks=9, error=True)
    assert "Fewer MeshBlocks (nmb_total=8) than MPI ranks (nranks=9)" in output


def test_refine_restart_rejects_rectangular_cells(tmp_path):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal", blocks=2)
                  + "\n<mesh>\nx1max = 1\n")
    output = run(tmp_path / "invalid", "ideal", restart=initial["path"],
                 refine=True, error=True)
    assert "requires cubic cells for native magnetic prolongation" in output


@pytest.mark.parametrize("mode", ["cgl", "adiabatic"])
@pytest.mark.parametrize("operation", ["refine", "reblock", "both"])
def test_refine_restart_preserves_forcing_lf_and_cooling_clocks(
        tmp_path, mode, operation):
    source = test_cgl_lf_cooling.input_text(nlim=3, history="dcycle = 100")
    source += """
<mesh>
nx1 = 8
nx2 = 8
nx3 = 8
<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8
<mhd>
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
<time>
sts_integrator = rkl2
<turb_driving>
dedt = 0.2
"""
    if mode == "adiabatic":
        source += test_cgl_lf_cooling.adiabatic_settings()
    initial = run(tmp_path / "initial", mode, source=source)
    assert initial["time"] > 0 and initial["cycle"] == 3
    reblock = operation in ("reblock", "both")
    flags = tuple(f"meshblock/nx{i}=4" for i in (1, 2, 3)) if reblock else ()
    refined = run(tmp_path / "fine", mode, restart=initial["path"],
                  refine=operation in ("refine", "both"), reblock=reblock,
                  flags=flags, ranks=2)
    assert refined["time"] == initial["time"] and refined["cycle"] == initial["cycle"]
    assert refined["metadata"] == initial["metadata"]
    for key in ("cooling_alpha", "cooling_history_last_time",
                "cooling_history_last_cycle", "cooling_target_fraction",
                "cooling_sample_time"):
        pattern = rf"^{key}\s*=\s*(\S+)"
        before = re.search(pattern, initial["text"], re.MULTILINE)
        after = re.search(pattern, refined["text"], re.MULTILINE)
        assert before and after and before.group(1) == after.group(1), key
    continued = run(tmp_path / "continued", mode, restart=refined["path"],
                    flags=("time/nlim=4",))
    assert continued["cycle"] == 4 and continued["time"] > refined["time"]
