"""Multi-cycle grid regression for weak-upwind anisotropy transport (T-B4)."""

from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.xfail(
    strict=True,
    reason="Known T-B4 multi-cycle weak-field transport failure; blocks WO1 merge "
           "pending flux design",
)
@pytest.mark.parametrize("velocity", [10.0, -10.0])
def test_weak_field_contact_remains_bounded(tmp_path, velocity):
    result = subprocess.run(
        [str(Path("athena").resolve()), "-i",
         str(ROOT / "inputs/unit_tests/cgl_weak_field_transport.athinput"),
         f"problem/weak_field_velocity={velocity}"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert int(re.findall(r"cycle=(\d+)", result.stdout)[-1]) == 50
    snapshots = sorted((tmp_path / "tab").glob("*.tab"))
    cycles = set()
    for path in snapshots:
        header = path.read_text().splitlines()
        cycles.add(int(re.search(r"cycle=(\d+)", header[0])[1]))
        assert header[1].split()[-6:] == [
            "dens", "velx", "vely", "velz", "eint", "p_perp",
        ]
        # The 1D table has gid, i, x1v followed by those six primitive fields.
        data = np.loadtxt(path)
        assert data.shape == (128, 9)
        assert np.all(np.isfinite(data)), path
        assert np.min(data[:, [3, 7, 8]]) > 1.0e-10, path
        ratio = data[:, 8] / data[:, 7]
        assert np.min(ratio) >= 0.5, (path, np.min(ratio))
        assert np.max(ratio) <= 2.0, (path, np.max(ratio))
    assert cycles == set(range(51))
