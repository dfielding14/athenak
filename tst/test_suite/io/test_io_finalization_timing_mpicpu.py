"""MPI regression for opt-in output timing reduction and final policy."""

from pathlib import Path
import subprocess


INPUT_FILE = "inputs/io_finalization_timing.athinput"


def test_mpi_timing_reports_one_rank_maximum_per_event(tmp_path: Path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            INPUT_FILE,
            "-d",
            str(run_dir),
            "mesh/nx1=32",
            "meshblock/nx1=16",
            "time/output_timing=true",
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert proc.stdout.count("[output-io] event=initial ") == 2
    assert proc.stdout.count("[output-io] event=final ") == 1
    assert (
        "event=initial block=output1 type=bin distribution=shared elapsed_max_s="
        in proc.stdout
    )
    assert (
        "event=final block=output2 type=rst distribution=shared elapsed_max_s="
        in proc.stdout
    )
    assert len(list((run_dir / "bin").glob("*.bin"))) == 1
    assert len(list((run_dir / "rst").glob("*.rst"))) == 2
