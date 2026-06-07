"""Helpers for exact upper-periodic-boundary Ito tracer tests."""

from pathlib import Path
import struct
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "tst/inputs/particles_ito_upper_boundary.athinput")
RESTART_HEADER = struct.Struct("@16s6iQ")
TARGETS = {
    0: (1.0, 0.75, 0.75),
    1: (0.75, 1.0, 0.75),
    2: (0.75, 0.75, 1.0),
    3: (1.0, 1.0, 1.0),
}
WRAPPED = {
    tag: tuple(0.0 if coordinate == 1.0 else coordinate for coordinate in position)
    for tag, position in TARGETS.items()
}


def latest_rank_zero_restart(run_dir):
    """Return the latest rank-zero restart and all matching rank files."""
    rank_zero = sorted(
        (run_dir / "rst/rank_00000000").glob("ito_upper_boundary.*.rst")
    )
    assert rank_zero
    selected = rank_zero[-1]
    rank_files = sorted((run_dir / "rst").glob(f"rank_*/{selected.name}"))
    assert rank_files
    return selected, rank_files


def _particle_section(path):
    data = bytearray(path.read_bytes())
    offset = data.rfind(b"ATHKPRTCLMC")
    assert offset >= 0
    values = RESTART_HEADER.unpack_from(data, offset)
    magic, version, enabled, nrdata, nidata, nlocal, nschedules, _ = values
    assert magic.rstrip(b"\0") == b"ATHKPRTCLMC"
    assert version == 4
    assert enabled == 2
    covariance_model = struct.unpack_from("@i", data, offset + RESTART_HEADER.size)[0]
    assert covariance_model == 0

    fixed_bytes = (
        RESTART_HEADER.size
        + struct.calcsize("@i")
        + 3 * nschedules * struct.calcsize("@i")
        + nidata * nlocal * struct.calcsize("@i")
        + nlocal * struct.calcsize("@Q")
    )
    real_count = nschedules + nrdata * nlocal
    assert real_count > 0
    real_size, remainder = divmod(len(data) - offset - fixed_bytes, real_count)
    assert remainder == 0
    assert real_size in (4, 8)

    real_base = (
        offset
        + RESTART_HEADER.size
        + struct.calcsize("@i")
        + nschedules * real_size
        + 3 * nschedules * struct.calcsize("@i")
    )
    int_base = real_base + nrdata * nlocal * real_size
    tag_base = int_base + nidata * nlocal * struct.calcsize("@i")
    tags = [
        struct.unpack_from("@Q", data, tag_base + p * struct.calcsize("@Q"))[0]
        for p in range(nlocal)
    ]
    gids = [
        struct.unpack_from("@i", data, int_base + p * struct.calcsize("@i"))[0]
        for p in range(nlocal)
    ]
    return data, real_base, real_size, nrdata, nlocal, tags, gids


def patch_exact_upper_coordinates(rank_files):
    """Place the four seeded particles exactly on upper faces and the corner."""
    initial_owner = {}
    initial_gid = {}
    patched = set()
    real_format = {4: "@f", 8: "@d"}
    for rank_file in rank_files:
        rank = int(rank_file.parent.name.removeprefix("rank_"))
        data, real_base, real_size, nrdata, nlocal, tags, gids = _particle_section(
            rank_file
        )
        assert nrdata >= 3
        for p, (tag, gid) in enumerate(zip(tags, gids)):
            assert tag in TARGETS
            initial_owner[tag] = rank
            initial_gid[tag] = gid
            patched.add(tag)
            for component, coordinate in zip((0, 2, 4), TARGETS[tag]):
                position_offset = real_base + (component * nlocal + p) * real_size
                struct.pack_into(
                    real_format[real_size], data, position_offset, coordinate
                )
        rank_file.write_bytes(data)

    assert patched == set(TARGETS)
    assert len(set(initial_owner.values())) == 1
    assert len(set(initial_gid.values())) == 1
    return initial_owner, initial_gid


def final_history(run_dir):
    """Return final-cycle particle state indexed by exact uint64 tag."""
    history = read_history(
        run_dir
        / "prtcl_thermo_history"
        / "ito_upper_boundary.prtcl_thermo_history.thp"
    )
    final_cycle = np.max(history["cycle"])
    keep = history["cycle"] == final_cycle
    state = {}
    for index in np.flatnonzero(keep):
        tag = int(history["tag"][index])
        state[tag] = {
            "gid": int(history["gid"][index]),
            "position": tuple(
                float(history[name][index]) for name in ("x1", "x2", "x3")
            ),
        }
    assert set(state) == set(TARGETS)
    return int(final_cycle), state


def assert_wrapped_state(state, initial_gid):
    """Check half-open coordinates and GID reassignment after exact-face wrapping."""
    for tag, expected in WRAPPED.items():
        np.testing.assert_array_equal(state[tag]["position"], expected)
        assert all(0.0 <= coordinate < 1.0 for coordinate in state[tag]["position"])
        assert state[tag]["gid"] != initial_gid[tag]
    assert len({entry["gid"] for entry in state.values()}) == len(state)


def restart_tag_owners(rank_files):
    """Return the rank owning each tag in a per-rank restart set."""
    owners = {}
    for rank_file in rank_files:
        rank = int(rank_file.parent.name.removeprefix("rank_"))
        _, _, _, _, _, tags, _ = _particle_section(rank_file)
        for tag in tags:
            assert tag not in owners
            owners[tag] = rank
    assert set(owners) == set(TARGETS)
    return owners
