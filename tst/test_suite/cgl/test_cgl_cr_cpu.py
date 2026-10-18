"""Frozen CGL restart integration and magnetic-only Boris accuracy.

Run with ATHENAK_TEST_EXE pointing to a double-precision built-in-pgen executable.
"""

import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[3]
TEMPLATE = REPO / "inputs/tests/cgl_lf_turb_driving_amr.athinput"
OVERLAY = REPO / "inputs/particles/frozen_cgl.athinput"
sys.path.insert(0, str(REPO / "vis/python"))
import bin_convert  # noqa: E402
import athena_read  # noqa: E402


def run(directory, *args, ranks=1, error=None):
    directory.mkdir(parents=True, exist_ok=True)
    exe = str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve())
    command = [exe, *map(str, args), "-d", str(directory)]
    if ranks > 1:
        command = ["mpirun", "-np", str(ranks), *command]
    result = subprocess.run(command, cwd=directory, capture_output=True, text=True,
                            check=False, timeout=120)
    output = result.stdout + result.stderr
    if error is None:
        assert result.returncode == 0, output
    else:
        assert result.returncode != 0, output
        assert error in output
    return result


def latest(directory, kind, pattern):
    paths = sorted((directory / kind).rglob(pattern))
    assert paths, f"no {pattern} in {directory / kind}"
    return paths[-1]


def state(directory):
    return bin_convert.read_binary(str(latest(directory, "bin", "*.state.*.bin")))


def checkpoint(directory):
    return latest(directory, "rst", "*.rst")


def background(directory, *, uniform=False, two_d=False, ranks=1, along_x=False,
               passive=False):
    directory.mkdir(parents=True, exist_ok=True)
    text = TEMPLATE.read_text().replace("nx1 = 24", "nx1 = 16")
    text = text.replace("nx2 = 24", "nx2 = 8")
    if not two_d:
        text = text.replace("nx3 = 1", "nx3 = 8")
    text = text.replace("method = user", "method = time\ntime_start = 0.015")
    text = text.replace("variable = turb_force", "variable = mhd_w_bcc")
    text = text.replace("id = force", "id = state")
    if passive:
        text = text.replace("passive = false", "passive = true\niso_sound_speed = 1")
        text = text.replace("pgen_name = divb_amr", "pgen_name = cgl_lf_paper\n"
                            "passive_delta = true\nbeta0 = 2\nb0 = 1")
        # Passive J/A does not support AMR in the CGL solver.
        start, end = text.index("<mesh_refinement>"), text.index("<time>")
        text = text[:start] + text[end:]
    deck = directory / "background.athinput"
    deck.write_text(text)
    flags = ["time/nlim=0", "time/tlim=0.1"]
    if uniform:
        flags += ["problem/field_amp=0", "problem/guide_b2=0",
                  f"problem/guide_b1={2 if along_x else 0}",
                  f"problem/guide_b3={0 if along_x else 2}"]
    run(directory, "-i", deck, *flags, ranks=ranks)
    assert state(directory)["n_mbs"] == 2
    return checkpoint(directory)


def particles(directory):
    records = []
    times = []
    for rank_dir in sorted((directory / "prst").glob("rank_*")):
        path = sorted(rank_dir.glob("*.prst"))[-1]
        raw = np.fromfile(path, dtype=np.float64)
        assert len(raw) == 3 + 17 * int(raw[2]), path
        times.append(raw[0])
        records.append(raw[3:].reshape(-1, 17))
    assert records, f"no particle checkpoints in {directory}"
    assert len(set(times)) == 1
    rows = np.concatenate(records)
    assert np.all(np.isfinite(rows))
    assert len(np.unique(rows[:, 1])) == len(rows)
    return times[0], rows[np.argsort(rows[:, 1])]


def frozen(directory, restart, end, *flags, ranks=1, error=None):
    directory.mkdir(parents=True, exist_ok=True)
    overlay = directory / "particles.athinput"
    overlay.write_text(OVERLAY.read_text() + "\n<particles>\n"
                       "exchange_gyro_only_substeps = true\n"
                       "subcycle_per_particle_gyro = false\n<problem>\nB0z = 73\n"
                       "<output98>\nfile_type = trk\ndt = -1\nnparticles = 17\n")
    return run(directory, "-r", restart, "-i", overlay, f"time/tlim={end}",
               "particles/ppc=0.015625", "particles/check_consistency_mode=full",
               "problem/particle_velocity=uniform", "problem/v0x=1",
               "problem/v0y=0", "problem/v0z=0.2", *flags, ranks=ranks, error=error)


def check_tracks(directory, rows, ntrack):
    frames = []
    with next((directory / "trk").glob("*.trk")).open("rb") as stream:
        while line := stream.readline():
            if not line.strip():
                continue
            assert line.startswith(b"# AthenaK tracked particle data")
            assert int(re.search(rb"ntracked_prtcls=\s*(\d+)", line)[1]) == ntrack
            assert not stream.readline().strip()  # Exactly one text separator.
            data = np.frombuffer(stream.read(24 * ntrack), dtype=np.float32)
            assert data.size == 6 * ntrack
            frames.append(data.reshape(ntrack, 6))
    assert frames and all(np.isfinite(frame).all() for frame in frames)
    np.testing.assert_array_equal(rows[:ntrack, 1], np.arange(ntrack))
    np.testing.assert_array_equal(frames[-1], rows[:ntrack, 3:9].astype(np.float32))


def active_restart_arrays(path):
    # This fixture has two 8^3 blocks, six CGL variables, no scalars,
    # and a shared double-precision restart. Its fluid/face payload is the tail.
    nx, ny, nz = 12, 12, 12
    shapes = [(6, nz, ny, nx), (nz, ny, nx + 1),
              (nz, ny + 1, nx), (nz + 1, ny, nx)]
    sizes = [int(np.prod(shape)) for shape in shapes]
    payload = path.read_bytes()[-2 * sum(sizes) * 8:]
    blocks = np.frombuffer(payload, dtype=np.float64).reshape(2, sum(sizes))
    arrays = []
    offset = 0
    k = slice(2, 10)
    kface = slice(2, 11)
    slices = [(slice(None), k, slice(2, 10), slice(2, 10)),
              (k, slice(2, 10), slice(2, 11)),
              (k, slice(2, 11), slice(2, 10)),
              (kface, slice(2, 10), slice(2, 10))]
    for shape, size, active in zip(shapes, sizes, slices):
        array = blocks[:, offset:offset + size].reshape(2, *shape)
        arrays.append(array[(slice(None), *active)])
        offset += size
    return arrays


@pytest.mark.parametrize("passive", [False, True])
def test_frozen_cgl_preserves_background_and_disables_evolution(tmp_path, passive):
    source = tmp_path / "source"
    restart = background(source, passive=passive)
    loaded = tmp_path / "loaded"
    frozen(loaded, restart, 0.03, "time/nlim=0")
    held = tmp_path / "frozen"
    frozen(held, restart, 0.03)
    assert state(held)["time"] == pytest.approx(0.03)
    np.testing.assert_array_equal(state(source)["mb_logical"], state(held)["mb_logical"])
    for before, after in zip(active_restart_arrays(restart),
                             active_restart_arrays(checkpoint(held))):
        np.testing.assert_allclose(after, before, rtol=2e-14, atol=2e-14)
    for before, after in zip(active_restart_arrays(checkpoint(loaded)),
                             active_restart_arrays(checkpoint(held))):
        np.testing.assert_array_equal(after, before)
    _, rows = particles(held)
    assert len(rows) == 16
    assert np.max(np.linalg.norm(rows[:, 13:16], axis=1)) > 0.02

    # Without particles, LF/forcing and supported AMR must still run.
    evolved = tmp_path / "evolved"
    run(evolved, "-r", restart, "time/nlim=-1", "time/tlim=0.02")
    if passive:
        assert state(evolved)["n_mbs"] == state(source)["n_mbs"]
    else:
        assert state(evolved)["n_mbs"] > state(source)["n_mbs"]
    history = athena_read.hst(str(next(evolved.glob("*.mhd.hst"))))
    assert history["lf_nstage"][-1] > 0
    assert not list(evolved.rglob("*.prst"))


@pytest.mark.parametrize("interpolation", ["tsc", "trilinear"])
def test_frozen_cgl_gathers_loaded_field_and_gyro_phase_converges(tmp_path, interpolation):
    restart = background(tmp_path / "source", uniform=True)
    errors = []
    seed = (f"particles/interpolation={interpolation}",
            "problem/particle_position=fixed", "problem/particle_x1=0.25",
            "problem/particle_x2=0.5", "problem/particle_x3=0.5", "problem/B0z=73")
    for cfl in (0.08, 0.04, 0.02):
        directory = tmp_path / f"gyro_{cfl}"
        frozen(directory, restart, 0.1, f"particles/cfl_part={cfl}",
               *seed)
        time, rows = particles(directory)
        assert time == pytest.approx(0.1)
        np.testing.assert_allclose(rows[:, 10:13], [[0, 0, 2]] * len(rows), atol=1e-13)
        np.testing.assert_allclose(np.linalg.norm(rows[:, 6:9], axis=1),
                                   np.sqrt(1.04), rtol=2e-13)
        expected = np.array([np.cos(20 * time), -np.sin(20 * time), 0.2])
        errors.append(np.max(np.linalg.norm(rows[:, 6:9] - expected, axis=1)))
        if cfl == 0.04:
            fine_steps = rows
    assert errors[1] < 0.4 * errors[0], errors
    assert errors[2] < 0.4 * errors[1], errors
    # The strict eight-substep boundary must not be rounded up to nine.
    for per_particle in (False, True):
        subcycled = tmp_path / f"subcycled_{per_particle}"
        frozen(subcycled, restart, 0.1, "particles/cfl_part=0.04",
               "particles/subcycle=true",
               f"particles/exchange_gyro_only_substeps={str(not per_particle).lower()}",
               f"particles/subcycle_per_particle_gyro={str(per_particle).lower()}", *seed)
        _, rows = particles(subcycled)
        np.testing.assert_allclose(rows[:, 3:], fine_steps[:, 3:], rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("per_particle", [False, True])
def test_frozen_cgl_particle_checkpoint_continuation(tmp_path, per_particle):
    restart = background(tmp_path / "source", uniform=True)
    uninterrupted = tmp_path / "uninterrupted"
    partial = tmp_path / "partial"
    resumed = tmp_path / "resumed"
    flags = ("particles/cfl_part=0.04", "problem/particle_position=tag_random",
             "particles/nspecies=2", "particles/mass_log_spacing=2",
             f"particles/subcycle={str(per_particle).lower()}",
             f"particles/subcycle_per_particle_gyro={str(per_particle).lower()}",
             f"particles/exchange_gyro_only_substeps={str(not per_particle).lower()}",
             "output98/dt=0.01")
    frozen(uninterrupted, restart, 0.08, *flags)
    # Stop by cycle so checkpointing does not insert a shortened physical step.
    frozen(partial, restart, 0.08, *flags, f"time/nlim={2 if per_particle else 10}")
    assert particles(partial)[0] < 0.08
    particle_restart = latest(partial / "prst", "rank_00000000", "*.prst")
    # Particle array order is arbitrary after migration or a GPU packing kernel.
    packed = np.fromfile(particle_restart, dtype=np.float64)
    particle_restart = tmp_path / "shuffled.prst"
    np.concatenate((packed[:3], packed[3:].reshape(-1, 17)[::-1].ravel())).tofile(
        particle_restart)
    frozen(resumed, checkpoint(partial), 0.08, *flags,
           "particles/seed_on_restart=false", "problem/prtcl_rst_flag=1",
           f"problem/prtcl_res_file={particle_restart}")
    ref_time, reference = particles(uninterrupted)
    resumed_time, result = particles(resumed)
    assert ref_time == resumed_time == pytest.approx(0.08)
    np.testing.assert_array_equal(result[:, :3], reference[:, :3])
    np.testing.assert_allclose(result[:, 3:], reference[:, 3:], rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(np.bincount(result[:, 2].astype(int)), [16, 16])
    np.testing.assert_array_equal(result[:, 9], 0.1 * 2 ** result[:, 2])
    np.testing.assert_allclose(np.linalg.norm(result[:, 6:9], axis=1),
                               np.sqrt(1.04), rtol=2e-13)
    phase = 2 * ref_time / result[:, 9]
    expected = np.column_stack((np.cos(phase), -np.sin(phase), np.full(len(result), 0.2)))
    # Boris phase error is O(phase * gyro_fraction**2 / 12).
    np.testing.assert_allclose(result[:, 6:9], expected, rtol=0, atol=2.2e-4)
    check_tracks(uninterrupted, reference, 17)
    check_tracks(resumed, result, 17)

    packed[0] += 1.0
    packed.tofile(particle_restart)
    frozen(tmp_path / "mismatched", checkpoint(partial), 0.08, *flags,
           "particles/seed_on_restart=false", "problem/prtcl_rst_flag=1",
           f"problem/prtcl_res_file={particle_restart}",
           error="Particle and fluid checkpoints must have the same time")


def test_frozen_cgl_rejects_two_dimensions(tmp_path):
    restart = background(tmp_path / "source", uniform=True, two_d=True)
    frozen(tmp_path / "rejected", restart, 0.1, error="require 3D periodic")


def test_cgl_particles_require_frozen_restart(tmp_path):
    source = tmp_path / "source"
    restart = background(source, uniform=True)
    frozen(tmp_path / "live", restart, 0.1, "particles/frozen_mhd=false",
           error="CGL particles require")
    deck = tmp_path / "direct.athinput"
    deck.write_text((source / "background.athinput").read_text() + OVERLAY.read_text())
    run(tmp_path / "direct", "-i", deck, error="require a restart")
