"""Balanced rich tracks must sample the endpoint without changing particle dynamics."""

import itertools
import re

import numpy as np
import pytest

from test_suite.cgl.test_cgl_cr_cpu import (
    OVERLAY, background, frozen, particles, state,
)


FIELDS = "tag,track_tag,species,mass,x,y,z,vx,vy,vz,bx,by,bz,dx,dy,dz,dpar"


def frames(path):
    result = []
    with path.open("rb") as stream:
        while line := stream.readline():
            if not line.strip():
                continue
            assert line.startswith(b"# AthenaK tracked particle data at time=")
            lines = [line]
            while (line := stream.readline()).strip():
                lines.append(line)
            header = dict(re.findall(rb"(\w+)=\s*([^\s]+)", b" ".join(lines)))
            assert header[b"trk_format"] == b"rich_v1"
            assert header[b"fields"].decode() == FIELDS
            assert int(header[b"nfields"]) == 17
            count = int(header[b"record_count"])
            data = np.frombuffer(stream.read(4*17*count), dtype="<f4")
            assert data.size == 17*count
            data = data.reshape(count, 17)
            result.append((float(header[b"time"]), int(header[b"cycle"]),
                           data[np.argsort(data[:, 0])]))
    return result


def endpoint_field(background_state, positions, method):
    # The fixture has two adjacent x blocks; sample its periodic global B independently.
    order = np.argsort(background_state["mb_logical"][:, 0])
    field = np.array([
        np.concatenate([background_state["mb_data"][name][m] for m in order], axis=2)
        for name in ("bcc1", "bcc2", "bcc3")
    ])
    shape = np.array(field.shape[1:][::-1])
    values = []
    for position in positions:
        coordinate = position*shape
        if method == "trilinear":
            base = np.floor(coordinate - 0.5).astype(int)
            fraction = coordinate - 0.5 - base
            weights = np.array([1 - fraction, fraction]).T
        else:
            center = np.floor(coordinate).astype(int)
            fraction = coordinate - center
            base = center - 1
            weights = np.array([0.5*(1 - fraction)**2,
                                0.75 - (fraction - 0.5)**2, 0.5*fraction**2]).T
        value = np.zeros(3)
        for offset in itertools.product(range(weights.shape[1]), repeat=3):
            i, j, k = (base + offset) % shape
            value += field[:, k, j, i]*np.prod(weights[np.arange(3), offset])
        values.append(value)
    return np.array(values)


@pytest.mark.parametrize("interpolation", ["tsc", "trilinear"])
def test_rich_tracks_cover_bins_and_preserve_push(tmp_path, monkeypatch, interpolation):
    restart = background(tmp_path / "source")
    overlay = tmp_path / "rich.athinput"
    overlay.write_text(OVERLAY.read_text() +
                       "\n<output98>\nparticles_per_species = 5\n")
    monkeypatch.setattr("test_suite.cgl.test_cgl_cr_cpu.OVERLAY", overlay)
    flags = ("particles/nspecies=3", "particles/mass_log_spacing=2",
             f"particles/interpolation={interpolation}", "output98/dt=0.01")
    tracked, control = tmp_path / "tracked", tmp_path / "control"
    frozen(tracked, restart, 0.05, *flags)
    frozen(control, restart, 0.05, *flags, "output98/particles_per_species=0")
    time, rows = particles(tracked)
    np.testing.assert_array_equal(rows, particles(control)[1])
    saved = frames(next((tracked / "trk").rglob("*.trk")))
    assert len(saved) > 2
    assert len(set(cycle for _, cycle, _ in saved)) == len(saved)
    assert saved[0][0] > 0
    assert saved[-1][0] == time == pytest.approx(0.05)
    ids = np.concatenate([np.arange(5) + 16*s for s in range(3)])
    for _, _, data in saved:
        assert np.isfinite(data).all()
        np.testing.assert_array_equal(data[:, 0], np.arange(15))
        np.testing.assert_array_equal(data[:, 1], ids)
        np.testing.assert_array_equal(data[:, 2], np.repeat(np.arange(3), 5))
        np.testing.assert_array_equal(data[:, 3], rows[ids, 9].astype(np.float32))
        np.testing.assert_allclose(np.linalg.norm(data[:, 7:10], axis=1),
                                   np.sqrt(1.04), rtol=2e-7)
    data = saved[-1][2]
    np.testing.assert_array_equal(data[:, 4:10], rows[ids, 3:9].astype(np.float32))
    np.testing.assert_array_equal(data[:, 13:17], rows[ids, 13:17].astype(np.float32))
    expected = endpoint_field(state(tmp_path / "source"), rows[ids, 3:6], interpolation)
    np.testing.assert_allclose(data[:, 10:13], expected, rtol=0, atol=2e-7)
    assert np.max(np.abs(data[:, 10:13] - rows[ids, 10:13])) > 1e-6


def test_rich_tracks_reject_invalid_count(tmp_path, monkeypatch):
    restart = background(tmp_path / "source", uniform=True)
    overlay = tmp_path / "invalid.athinput"
    overlay.write_text(OVERLAY.read_text() +
                       "\n<output98>\nparticles_per_species = -1\n")
    monkeypatch.setattr("test_suite.cgl.test_cgl_cr_cpu.OVERLAY", overlay)
    frozen(tmp_path / "negative", restart, 0.01, "output98/dt=0.01",
           error="particles_per_species must be nonnegative")
    frozen(tmp_path / "too_many", restart, 0.01, "output98/dt=0.01",
           "output98/particles_per_species=17",
           error="particles_per_species <= available particles per bin")
