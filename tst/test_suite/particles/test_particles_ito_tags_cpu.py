"""CPU regression tests for exact 64-bit Ito tracer tags."""

from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "tst/inputs/particles_ito_tags.athinput")
RUN_ROOT = Path("run_particles_ito_tags_cpu")
HISTORY_HEADER = struct.Struct("32siiii")
UINT64_MAX = 2**64 - 1


def _history_path(run_dir):
    return run_dir / "prtcl_thermo_history/ito_tags.prtcl_thermo_history.thp"


def _read_v3_history(run_dir):
    path = _history_path(run_dir)
    with path.open("rb") as handle:
        _, version, real_size, _, _ = HISTORY_HEADER.unpack(
            handle.read(HISTORY_HEADER.size)
        )
    assert version == 3
    assert real_size in (4, 8)
    data = read_history(path)
    assert data["tag"].dtype == np.uint64
    return data


def _read_vtk_tags(run_dir, count):
    paths = sorted((run_dir / "pvtk").glob("*.vtk"))
    assert paths
    payload = paths[-1].read_bytes()

    def read_word(name):
        marker = (
            f"SCALARS {name} unsigned_int\nLOOKUP_TABLE default\n".encode("ascii")
        )
        offset = payload.index(marker) + len(marker)
        return np.frombuffer(payload, dtype=">u4", count=count, offset=offset)

    low = read_word("ptag_low32").astype(np.uint64)
    high = read_word("ptag_high32").astype(np.uint64)
    return (high << np.uint64(32)) | low


def _state_at_last_cycle(run_dir):
    data = _read_v3_history(run_dir)
    keep = data["cycle"] == np.max(data["cycle"])
    order = np.argsort(data["tag"][keep])
    tags = data["tag"][keep][order]
    assert np.unique(tags).size == tags.size
    state = np.column_stack(
        (
            data["gid"][keep][order],
            data["x1"][keep][order],
            data["x2"][keep][order],
            data["x3"][keep][order],
        )
    )
    return tags, state


def _run(run_dir, next_tag, count, nlim=1):
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    return testutils.run(
        INPUT,
        [
            "-d",
            str(run_dir),
            f"particles/next_tracer_tag={next_tag}",
            f"tracer_seed1/count_per_event={count}",
            f"time/nlim={nlim}",
        ],
    )


def test_ito_uint64_tag_boundaries_in_thermo_history_v3():
    """Thermo history preserves tags across every integer-width boundary."""
    cases = (
        (2**31 - 2, 3),
        (2**32 - 2, 3),
        (2**63 - 1, 3),
        (UINT64_MAX - 1, 1),
    )
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        for index, (first_tag, count) in enumerate(cases):
            run_dir = RUN_ROOT / f"boundary_{index}"
            assert _run(run_dir, first_tag, count)
            data = _read_v3_history(run_dir)
            cycle_zero = data["tag"][data["cycle"] == 0]
            expected = np.arange(
                first_tag, first_tag + count, dtype=np.uint64
            )
            np.testing.assert_array_equal(np.sort(cycle_zero), expected)
            np.testing.assert_array_equal(np.sort(_read_vtk_tags(run_dir, count)),
                                          expected)
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)


def test_ito_uint64_tags_restart_round_trip_and_trajectory():
    """Restart preserves exact tags and reproduces the uninterrupted trajectory."""
    first_tag = 2**63 - 17
    count = 64
    split_dir = RUN_ROOT / "split"
    round_trip_dir = RUN_ROOT / "round_trip"
    resumed_dir = RUN_ROOT / "resumed"
    uninterrupted_dir = RUN_ROOT / "uninterrupted"
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        assert _run(split_dir, first_tag, count, nlim=1)
        restart = sorted(
            (split_dir / "rst/rank_00000000").glob("ito_tags.*.rst")
        )[-1]

        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(restart),
                "-d",
                str(round_trip_dir),
                "time/nlim=1",
            ]
        )
        split_tags, split_state = _state_at_last_cycle(split_dir)
        round_trip_tags, round_trip_state = _state_at_last_cycle(round_trip_dir)
        np.testing.assert_array_equal(round_trip_tags, split_tags)
        np.testing.assert_array_equal(round_trip_state, split_state)

        assert testutils.run_command(
            [
                "./athena",
                "-r",
                str(restart),
                "-d",
                str(resumed_dir),
                "time/nlim=2",
            ]
        )
        assert _run(uninterrupted_dir, first_tag, count, nlim=2)
        resumed_tags, resumed_state = _state_at_last_cycle(resumed_dir)
        uninterrupted_tags, uninterrupted_state = _state_at_last_cycle(
            uninterrupted_dir
        )
        np.testing.assert_array_equal(resumed_tags, uninterrupted_tags)
        np.testing.assert_array_equal(resumed_state, uninterrupted_state)
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)


def test_ito_uint64_tag_exhaustion_fails_closed():
    """Seeding at UINT64_MAX is rejected before tags can wrap to zero."""
    run_dir = RUN_ROOT / "exhausted"
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    try:
        run_dir.parent.mkdir(parents=True)
        result = subprocess.run(
            [
                "./athena",
                "-i",
                INPUT,
                "-d",
                str(run_dir),
                f"particles/next_tracer_tag={UINT64_MAX}",
                "tracer_seed1/count_per_event=1",
                "time/nlim=1",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "particle tag space exhausted" in result.stdout + result.stderr
        assert not list(run_dir.rglob("*.thp"))
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)
