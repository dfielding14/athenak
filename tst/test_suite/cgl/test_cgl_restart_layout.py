"""CGL checkpoint slots retain their EOS and scalar interpretation on restart."""

import os
from pathlib import Path
import subprocess

import pytest

from test_suite.cgl.test_refine_on_restart import input_text


def run(directory, *, source=None, restart=None, flags=(), error=False):
    directory.mkdir()
    command = [str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()),
               "-d", str(directory)]
    if source is not None:
        deck = directory / "input.athinput"
        deck.write_text(source)
        command += ["-i", str(deck)]
    if restart is not None:
        command += ["-r", str(restart)]
    result = subprocess.run(command + list(flags), cwd=directory,
                            capture_output=True, text=True, check=False, timeout=30)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output)
    if error:
        assert result.returncode != 0 and "CGL restart layout" in output, output
        return None
    assert result.returncode == 0, output
    return sorted((directory / "rst").rglob("*.rst"))[-1]


@pytest.mark.parametrize("mode", ["cgl", "ideal"])
@pytest.mark.parametrize("transfer", ["shared", "rank_local", "refine", "reblock"])
def test_cgl_restart_rejects_equal_sized_eos_alias(tmp_path, mode, transfer):
    scalars = 0 if mode == "cgl" else 1
    source = input_text(mode) + f"\n<mhd>\nnscalars = {scalars}\n"
    if transfer == "rank_local":
        source += "\n<output2>\nsingle_file_per_rank = true\n"
    checkpoint = run(tmp_path / "initial", source=source)
    flags = [f"mhd/eos={'ideal' if mode == 'cgl' else 'cgl'}",
             f"mhd/nscalars={1-scalars}"]
    if transfer == "refine":
        flags += ["mesh/refine_on_restart=true"]
    elif transfer == "reblock":
        flags += ["mesh/reblock_on_restart=true", "meshblock/nx1=4"]
    run(tmp_path / "invalid", restart=checkpoint, flags=flags, error=True)


def test_cgl_restart_rejects_changed_scalar_count_through_overlay(tmp_path):
    checkpoint = run(tmp_path / "initial", source=input_text("cgl"))
    overlay = "\n<mhd>\nnscalars = 1\n"
    run(tmp_path / "invalid", source=overlay, restart=checkpoint, error=True)
