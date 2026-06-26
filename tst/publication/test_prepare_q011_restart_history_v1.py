#!/usr/bin/env python3
"""Tests for Q011 restart-history continuation preparation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

try:
    from tst.publication import prepare_q011_restart_history_v1 as preparation
except ModuleNotFoundError:
    import prepare_q011_restart_history_v1 as preparation


_MESH_FORMAT = "<ii9d19i19iddii"
_BANNER = b"# Athena++ history data\n"
_HEADER = b"#  [1]=time      [2]=dt       [3]=mass\n"


def _row(time: float, mass: float, *, timestep: float = 0.1) -> bytes:
    return f" {time:.17e} {timestep:.17e} {mass:.17e}\n".encode("ascii")


def _history(*segments: list[tuple[float, float]]) -> bytes:
    payload = bytearray()
    for segment in segments:
        payload.extend(_BANNER)
        payload.extend(_HEADER)
        for time, mass in segment:
            payload.extend(_row(time, mass))
    return bytes(payload)


def _write_restart(
    path: Path,
    *,
    time: float = 2.0,
    timestep: float = 0.1,
    cycle: int = 20,
    meshblocks: int = 4,
    ranks: int = 2,
    sparse_size: int | None = None,
) -> bytes:
    mesh_header = struct.pack(
        _MESH_FORMAT,
        meshblocks,
        1,
        *([0.0] * 9),
        *([0] * 19),
        *([0] * 19),
        time,
        timestep,
        cycle,
        ranks,
    )
    prefix = b"<job>\nbasename=q011-history-fixture\n<par_end>\n" + mesh_header
    with path.open("wb") as stream:
        stream.write(prefix)
        if sparse_size is not None:
            if sparse_size <= len(prefix):
                raise AssertionError("sparse fixture must exceed its restart header")
            stream.seek(sparse_size - 1)
            stream.write(b"\0")
        else:
            stream.write(b"synthetic-restart-payload")
    return prefix


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class _GuardedReader:
    def __init__(self, stream, read_sizes: list[int]):
        self._stream = stream
        self._read_sizes = read_sizes

    def __enter__(self):
        self._stream.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return self._stream.__exit__(exc_type, exc_value, traceback)

    def read(self, size: int = -1) -> bytes:
        if size <= 0 or size > preparation._IO_CHUNK_SIZE:
            raise AssertionError(f"unbounded restart read requested: {size}")
        self._read_sizes.append(size)
        return self._stream.read(size)

    def __getattr__(self, name: str):
        return getattr(self._stream, name)


class PrepareRestartHistoryTests(unittest.TestCase):
    def _paths(self, root: Path) -> tuple[Path, Path, Path, Path]:
        return (
            root / "q011.rst",
            root / "q011.user.hst",
            root / "q011.user.hst.pre-continuation",
            root / "q011.history-preparation.json",
        )

    def test_prunes_deduplicates_archives_and_writes_canonical_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            restart, history, archive, receipt = self._paths(root)
            sparse_size = 8 * 1024 * 1024
            _write_restart(restart, sparse_size=sparse_size)
            original = _history(
                [(0.0, 10.0), (1.0, 11.0), (2.0, 12.0), (3.0, 13.0)],
                [(2.0, 12.0), (4.0, 14.0)],
            )
            history.write_bytes(original)

            real_open = preparation._open_binary
            read_sizes: list[int] = []

            def guarded_open(path: Path):
                stream = real_open(path)
                if Path(path) == restart:
                    return _GuardedReader(stream, read_sizes)
                return stream

            with mock.patch.object(preparation, "_open_binary", side_effect=guarded_open), mock.patch.object(
                Path,
                "read_bytes",
                side_effect=AssertionError("full-file helper must not be used"),
            ):
                record = preparation.prepare_restart_history(
                    restart=restart,
                    history=history,
                    archive=archive,
                    receipt=receipt,
                )

            expected = (
                _BANNER
                + b"#  [1]=time [2]=dt [3]=mass\n"
                + _row(0.0, 10.0)
                + _row(1.0, 11.0)
                + _row(2.0, 12.0)
            )
            self.assertEqual(history.read_bytes(), expected)
            self.assertEqual(archive.read_bytes(), original)
            self.assertTrue(read_sizes)
            self.assertLessEqual(max(read_sizes), preparation._IO_CHUNK_SIZE)

            receipt_payload = receipt.read_bytes()
            observed = json.loads(receipt_payload)
            self.assertEqual(observed, record)
            self.assertEqual(receipt_payload, preparation._canonical_json(observed))
            self.assertEqual(record["record_type"], preparation.RECORD_TYPE)
            self.assertEqual(record["restart"]["byte_count"], sparse_size)
            self.assertEqual(record["restart"]["cycle"], 20)
            self.assertEqual(record["restart"]["time"], 2.0)
            self.assertEqual(record["history"]["input_row_count"], 6)
            self.assertEqual(record["history"]["retained_row_count"], 3)
            self.assertEqual(record["history"]["dropped_row_count"], 2)
            self.assertEqual(record["history"]["duplicate_row_count"], 1)
            self.assertEqual(record["history"]["input_sha256"], _sha256(original))
            self.assertEqual(record["history"]["normalized_sha256"], _sha256(expected))
            self.assertEqual(record["archive"]["sha256"], _sha256(original))

    def test_tight_tolerance_collapses_restart_appended_time(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            restart, history, archive, receipt = self._paths(root)
            _write_restart(restart)
            nearby = math.nextafter(2.0, math.inf)
            original = _history(
                [(0.0, 10.0), (1.0, 11.0), (2.0, 12.0)],
                [(nearby, 12.0), (2.0 + 1.0e-8, 13.0)],
            )
            history.write_bytes(original)

            record = preparation.prepare_restart_history(
                restart=restart,
                history=history,
                archive=archive,
                receipt=receipt,
            )

            self.assertEqual(record["history"]["retained_row_count"], 3)
            self.assertEqual(record["history"]["duplicate_row_count"], 1)
            self.assertEqual(record["history"]["dropped_row_count"], 1)
            lines = history.read_text(encoding="ascii").splitlines()
            self.assertEqual(len(lines), 5)
            self.assertEqual([float(line.split()[0]) for line in lines[2:]], [0.0, 1.0, 2.0])

    def test_malformed_history_fails_without_mutation(self) -> None:
        cases = {
            "nonfinite": _BANNER + _HEADER + b"0.0 0.1 1e999\n",
            "row width": _BANNER + _HEADER + b"0.0 0.1\n",
            "column drift": _history([(0.0, 1.0)])
            + _BANNER
            + b"#  [1]=time [2]=dt [3]=energy\n"
            + _row(1.0, 2.0),
            "nonmonotonic": _history([(1.0, 1.0), (0.5, 2.0)]),
            "conflicting state": _history([(0.0, 1.0), (2.0, 2.0)], [(2.0, 3.0)]),
        }
        for name, original in cases.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                restart, history, archive, receipt = self._paths(root)
                _write_restart(restart)
                history.write_bytes(original)

                with self.assertRaises(preparation.PreparationError):
                    preparation.prepare_restart_history(
                        restart=restart,
                        history=history,
                        archive=archive,
                        receipt=receipt,
                    )

                self.assertEqual(history.read_bytes(), original)
                self.assertFalse(archive.exists())
                self.assertFalse(receipt.exists())
                self.assertEqual(list(root.glob(".*.tmp")), [])

    def test_invalid_restart_chronology_or_cardinality_fails_closed(self) -> None:
        invalid = {
            "negative cycle": {"cycle": -1},
            "zero timestep": {"timestep": 0.0},
            "nonfinite time": {"time": math.nan},
            "too many ranks": {"meshblocks": 2, "ranks": 3},
        }
        for name, overrides in invalid.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                restart, history, archive, receipt = self._paths(root)
                _write_restart(restart, **overrides)
                original = _history([(0.0, 1.0)])
                history.write_bytes(original)

                with self.assertRaisesRegex(
                    preparation.PreparationError, "chronology or cardinality"
                ):
                    preparation.prepare_restart_history(
                        restart=restart,
                        history=history,
                        archive=archive,
                        receipt=receipt,
                    )

                self.assertEqual(history.read_bytes(), original)
                self.assertFalse(archive.exists())
                self.assertFalse(receipt.exists())

    def test_publish_failures_restore_original_and_remove_outputs(self) -> None:
        for failing_destination in ("history", "receipt"):
            with self.subTest(destination=failing_destination), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                restart, history, archive, receipt = self._paths(root)
                _write_restart(restart)
                original = _history([(0.0, 1.0), (2.0, 2.0), (3.0, 3.0)])
                history.write_bytes(original)
                destination = history if failing_destination == "history" else receipt
                real_replace = preparation._atomic_replace

                def injected_failure(source: Path, target: Path) -> None:
                    if Path(target) == destination:
                        raise OSError("injected publication failure")
                    real_replace(source, target)

                with mock.patch.object(
                    preparation, "_atomic_replace", side_effect=injected_failure
                ), self.assertRaisesRegex(OSError, "injected publication failure"):
                    preparation.prepare_restart_history(
                        restart=restart,
                        history=history,
                        archive=archive,
                        receipt=receipt,
                    )

                self.assertEqual(history.read_bytes(), original)
                self.assertFalse(archive.exists())
                self.assertFalse(receipt.exists())
                self.assertEqual(list(root.glob(".*.tmp")), [])

    def test_cli_accepts_only_named_paths_and_emits_the_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            restart, history, archive, receipt = self._paths(root)
            _write_restart(restart, time=1.0, cycle=10)
            history.write_bytes(_history([(0.0, 1.0), (1.0, 2.0)]))
            script = Path(preparation.__file__).resolve()

            completed = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--restart",
                    str(restart),
                    "--history",
                    str(history),
                    "--archive",
                    str(archive),
                    "--receipt",
                    str(receipt),
                ],
                check=False,
                capture_output=True,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr.decode())
            self.assertEqual(completed.stdout, receipt.read_bytes())
            self.assertEqual(json.loads(completed.stdout)["restart"]["cycle"], 10)


if __name__ == "__main__":
    unittest.main()
