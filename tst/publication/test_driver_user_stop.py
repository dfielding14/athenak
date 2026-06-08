from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
DECK = REPO_ROOT / "inputs/tests/driver_user_stop.athinput"


def _run(tmp_path: Path, *, failure: bool) -> subprocess.CompletedProcess[str]:
    executable_dir = os.environ.get("ATHENA_DRIVER_USER_STOP_EXE_DIR")
    if not executable_dir:
        pytest.skip("ATHENA_DRIVER_USER_STOP_EXE_DIR is required")
    executable = Path(executable_dir) / "athena"
    if not executable.is_file():
        raise RuntimeError(f"driver user-stop executable is missing: {executable}")
    run_dir = tmp_path / ("failure" if failure else "success")
    run_dir.mkdir()
    return subprocess.run(
        [
            str(executable),
            "-i",
            str(DECK),
            "-d",
            str(run_dir),
            f"problem/test_stop_failure={'true' if failure else 'false'}",
        ],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


@pytest.mark.parametrize(("failure", "expected_exit"), [(False, 0), (True, 1)])
def test_orderly_user_stop_publishes_final_state_and_exit_disposition(
    tmp_path: Path, failure: bool, expected_exit: int
) -> None:
    result = _run(tmp_path, failure=failure)
    assert result.returncode == expected_exit, result.stdout
    assert "Terminating on user request" in result.stdout
    assert (
        f"user_stop_reason_code=17 user_stop_failure={'true' if failure else 'false'}"
        in result.stdout
    )
    assert "cycle=2" in result.stdout
    run_dir = tmp_path / ("failure" if failure else "success")
    restarts = sorted((run_dir / "rst").glob("driver_user_stop.*.rst"))
    tables = sorted((run_dir / "tab").glob("driver_user_stop.state.*.tab"))
    # Initialization writes sequence 00000; orderly finalization must add 00001.
    assert len(restarts) >= 2
    assert len(tables) >= 2
    assert all(path.stat().st_size > 0 for path in [*restarts, *tables])
