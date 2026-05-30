"""Regression tests for final-output policy and opt-in IO timing."""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
INPUT_FILE = "inputs/io_finalization_timing.athinput"


def _run_case(tmp_path: Path, *overrides: str):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", INPUT_FILE, "-d", str(run_dir), *overrides],
        check=True,
        capture_output=True,
        text=True,
    )
    bin_files = sorted((run_dir / "bin").glob("*.bin"))
    rst_files = sorted((run_dir / "rst").glob("*.rst"))
    return proc.stdout, bin_files, rst_files


def test_default_policy_retains_terminal_outputs(tmp_path):
    stdout, bin_files, rst_files = _run_case(tmp_path)
    assert "[output-io]" not in stdout
    assert len(bin_files) == 2
    assert len(rst_files) == 2


def test_restart_only_suppresses_terminal_diagnostic(tmp_path):
    _, bin_files, rst_files = _run_case(
        tmp_path, "time/final_output_policy=restart_only"
    )
    assert len(bin_files) == 1
    assert len(rst_files) == 2


def test_none_suppresses_all_terminal_outputs(tmp_path):
    _, bin_files, rst_files = _run_case(tmp_path, "time/final_output_policy=none")
    assert len(bin_files) == 1
    assert len(rst_files) == 1


def test_timing_is_opt_in(tmp_path):
    stdout, _, _ = _run_case(tmp_path, "time/output_timing=true")
    assert stdout.count("[output-io] event=initial ") == 2
    assert stdout.count("[output-io] event=final ") == 2
    assert (
        "event=initial block=output1 type=bin distribution=shared "
        "elapsed_max_s=" in stdout
    )
    assert (
        "event=initial block=output2 type=rst distribution=shared "
        "elapsed_max_s=" in stdout
    )


def test_invalid_final_output_policy_is_rejected(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "time/final_output_policy=invalid",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "final_output_policy = 'invalid' not implemented" in proc.stdout


def test_terminal_restart_resume_advances_counter_without_overwrite(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    subprocess.run(
        [
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    terminal_restart = run_dir / "rst" / "io_policy.00001.rst"
    saved_terminal_bytes = terminal_restart.read_bytes()
    subprocess.run(
        [
            "./athena",
            "-r",
            str(terminal_restart),
            "-d",
            str(run_dir),
            "time/tlim=0.02",
            "time/nlim=2",
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert terminal_restart.read_bytes() == saved_terminal_bytes
    assert (run_dir / "rst" / "io_policy.00002.rst").exists()


def test_origin_main_shared_restart_fixture_resumes(tmp_path):
    run_dir = tmp_path / "origin_main_shared_resume"
    run_dir.mkdir()
    restart = FIXTURES / "rst" / "shared" / "io_legacy_shared.00001.rst"
    subprocess.run(
        [
            "./athena",
            "-r",
            str(restart),
            "-d",
            str(run_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
