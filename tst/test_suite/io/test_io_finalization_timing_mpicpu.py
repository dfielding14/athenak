"""MPI regression for opt-in performance timing and final policy."""

from pathlib import Path
import subprocess


INPUT_FILE = "inputs/io_finalization_timing.athinput"


def test_mpi_timing_reports_one_rank_maximum(tmp_path: Path):
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
            "time/performance_timing=true",
            "time/final_output_policy=restart_only",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert proc.stdout.count("PERFORMANCE_TIMING synchronized=true") == 1
    assert proc.stdout.count("PERFORMANCE_REGION outputs ") == 1
    assert len(list((run_dir / "bin").glob("*.bin"))) == 1
    assert len(list((run_dir / "rst").glob("*.rst"))) == 2
