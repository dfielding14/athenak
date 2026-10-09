"""Multi-cycle grid regression for weak-upwind anisotropy transport (T-B4)."""

from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.xfail(
    strict=True,
    reason="Accepted T-B4 sharp-contact limitation of A transport; see "
           "docs/cgl_lf_surgical_rebuild_work_order.md",
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


@pytest.mark.parametrize("velocity", [10.0, -10.0])
def test_weak_field_first_cycle(tmp_path, velocity):
    result = subprocess.run(
        [str(Path("athena").resolve()), "-i",
         str(ROOT / "inputs/unit_tests/cgl_weak_field_transport.athinput"),
         f"problem/weak_field_velocity={velocity}", "time/nlim=1"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert int(re.findall(r"cycle=(\d+)", result.stdout)[-1]) == 1
    snapshots = sorted((tmp_path / "tab").glob("*.tab"))
    cycles = {int(re.search(r"cycle=(\d+)", path.read_text())[1]) for path in snapshots}
    assert cycles == {0, 1}
    data = np.loadtxt(snapshots[-1])
    assert np.all(np.isfinite(data))
    assert np.min(data[:, [3, 7, 8]]) > 1.0e-10
    ratio = data[:, 8] / data[:, 7]
    assert np.min(ratio) >= 0.5
    assert np.max(ratio) <= 2.0


def test_weak_field_smooth_interface_converges(tmp_path):
    """Resolved A mixing converges; the archived smooth case is fully magnetized."""
    executable = str(Path("athena").resolve())
    source = (ROOT / "inputs/unit_tests/cgl_weak_field_transport.athinput").read_text()
    source = source.replace("variable = mhd_w", "variable = mhd_w_bcc")
    source = source.replace("data_format = %24.16e\ndcycle = 1",
                            "data_format = %24.16e\ndcycle = 1000000")
    source += (
        "\n<output3>\nfile_type = tab\nvariable = mhd_u\n"
        "data_format = %24.16e\ndcycle = 1000000\n"
    )
    errors = []
    for nx in (128, 256, 512):
        mirrored = []
        for velocity in (10.0, -10.0):
            case = tmp_path / f"nx{nx}_v{velocity:+g}"
            case.mkdir()
            input_path = case / "smooth.athinput"
            input_path.write_text(source)
            result = subprocess.run(
                [executable, "-i", str(input_path), f"mesh/nx1={nx}",
                 "meshblock/nx1=128", f"problem/weak_field_velocity={velocity}",
                 "problem/profile_width=0.04", "mhd/bfloor=1.0e-14",
                 "time/nlim=10000", "time/tlim=0.012"],
                cwd=case, capture_output=True, text=True, check=False,
            )
            (case / "run.log").write_text(result.stdout + result.stderr)
            assert result.returncode == 0, result.stdout + result.stderr
            snapshots = sorted((case / "tab").glob("*.mhd_w_bcc.*.tab"))
            assert len(snapshots) == 2
            time = float(re.search(r"time=\s*(\S+)", snapshots[-1].read_text())[1])
            assert abs(time - 0.012) < 1.0e-14
            data = np.loadtxt(snapshots[-1])
            data = data[np.argsort(data[:, 2])]
            assert data.shape == (nx, 12)
            assert np.all(np.isfinite(data))
            assert np.min(data[:, [3, 7, 8]]) > 1.0e-10
            sign = np.sign(velocity)
            exact_b = 1.0e-12 + (1.0 - 1.0e-12) * 0.5 * (
                1.0 + np.tanh(sign * (data[:, 2] - velocity*time - 0.5) / 0.04)
            )
            exact_p = 1.5 - 0.5*exact_b**2
            error = np.array([
                np.mean(np.abs(np.log(data[:, 8] / data[:, 7]))),
                np.mean(np.abs(data[:, 10] - exact_b)),
                np.mean(np.abs(data[:, 7] - exact_p)),
                np.mean(np.abs(data[:, 8] - exact_p)),
            ])
            if velocity > 0:
                errors.append(error)
            state = data[:, 3:].copy()
            state[:, 1] *= sign
            mirrored.append(state if velocity > 0 else state[::-1])
            conserved = sorted((case / "tab").glob("*.mhd_u.*.tab"))
            initial_u, final_u = np.loadtxt(conserved[0]), np.loadtxt(conserved[-1])
            # Equal boundary mass/momentum fluxes; net energy inflow is |v|/4.
            np.testing.assert_allclose(final_u[:, 3:7].mean(axis=0),
                                       initial_u[:, 3:7].mean(axis=0), atol=1e-10, rtol=0)
            assert abs(np.mean(final_u[:, 7] - initial_u[:, 7])
                       - 0.25*abs(velocity)*time) < 1.0e-9
            print(f"nx={nx} velocity={velocity:g} errors(logR,B,ppar,pperp)={error}")
        np.testing.assert_allclose(mirrored[0], mirrored[1], atol=2e-10, rtol=2e-10)
    errors = np.asarray(errors)
    assert np.all(errors > 0.0)
    # The archived fixed-width profile gives approximately second-order log-ratio errors.
    assert np.all(errors[1:, 0] < 0.4*errors[:-1, 0]), errors
    assert np.all(errors[1:, 1:] < 0.6*errors[:-1, 1:]), errors
