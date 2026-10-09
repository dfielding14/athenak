"""CGL startup guards using a supplied built_in_pgens executable.

Set ATHENAK_TEST_EXE to a built executable or launch adapter; otherwise use
./athena in the usual staged test directory. This module performs no build.
"""

import os
from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize(
    ("lf", "parameters", "error"),
    [
        pytest.param(True, (f"{name}={value}",), name, id=f"{name}-{value}")
        for name in ("dfloor", "pfloor", "bfloor", "sigma_max",
                     "lf_k_parallel", "lf_c_parallel0", "firehose_threshold",
                     "mirror_threshold")
        for value in ("0", "-1", "nan", "inf")
    ] + [
        pytest.param(True, (f"{name}={value}",), name, id=f"{name}-{value}")
        for name in ("tfloor", "sfloor", "nu_coll", "limiter_nu_coll",
                     "limiter_backup_nu")
        for value in ("-1", "nan", "inf")
    ] + [
        pytest.param(True, (f"{limiter}=true", f"limiter_nu_coll={value}"),
                     "limiter_nu_coll", id=f"{limiter}-rate-{value}")
        for limiter in ("mirror_limiter", "firehose_limiter")
        for value in ("-1", "nan", "inf")
    ] + [
        pytest.param(lf, ("passive=true", *sound_speed), "passive",
                     id=f"passive-lf-{lf}-sound-speed-{bool(sound_speed)}")
        for lf in (False, True)
        for sound_speed in ((), ("iso_sound_speed=1",))
    ] + [
        pytest.param(True, (f"{name}={value}",), name, id=f"{name}-{value}")
        for name in ("firehose_backup_factor", "mirror_backup_factor")
        for value in ("0.99", "nan", "inf")
    ] + [
        pytest.param(True, (f"{limiter}=true",),
                     ("Parameter name 'limiter_nu_coll' not found",),
                     id=f"{limiter}-missing-rate")
        for limiter in ("mirror_limiter", "firehose_limiter")
    ] + [
        pytest.param(True, ("limiter_hardwall=true",),
                     ("limiter_hardwall=true is no longer supported", "limiter_nu_coll=1e10"),
                     id="legacy-hardwall-rejected"),
        pytest.param(True, ("backup_limiters=true",),
                     ("<mhd>/backup_limiters requires",), id="backup-needs-limiter"),
        pytest.param(True, ("cgl_firehose_threshold=oblique", "firehose_threshold=2"),
                     ("<mhd>/firehose_threshold conflicts",), id="conflicting-legacy"),
        pytest.param(True, ("cgl_firehose_threshold=invalid",),
                     ("<mhd>/cgl_firehose_threshold must be",), id="invalid-legacy"),
        pytest.param(True, (), None, id="active-default"),
        pytest.param(False, ("passive=false",), None, id="active-no-lf"),
        pytest.param(True, ("tfloor=0", "sfloor=0", "nu_coll=0",
                           "limiter_nu_coll=0"), None, id="zero-floors-and-rates"),
        pytest.param(True, ("limiter_nu_coll=7",), None, id="disabled-rate"),
        pytest.param(True, ("nu_coll=3",), None, id="background-rate"),
        pytest.param(True, ("firehose_threshold=2", "mirror_threshold=1",
                           "firehose_backup_factor=1", "mirror_backup_factor=1",
                           "limiter_backup_nu=0"), None, id="numeric-wall-boundaries"),
        pytest.param(True, ("cgl_firehose_threshold=parallel", "firehose_threshold=2"),
                     None, id="matching-parallel"),
        pytest.param(True, ("cgl_firehose_threshold=oblique", "firehose_threshold=1.4"),
                     None, id="matching-oblique"),
        pytest.param(True, ("backup_limiters=true", "mirror_limiter=true",
                           "limiter_nu_coll=0"), None, id="backup-with-zero-soft-rate"),
        pytest.param(True, ("firehose_limiter=true", "limiter_nu_coll=0"), None,
                     id="firehose-zero-rate"),
        pytest.param(True, ("mirror_limiter=true", "limiter_nu_coll=3"), None,
                     id="mirror-positive-rate"),
    ],
)
def test_cgl_constructor_parameters(tmp_path, lf, parameters, error):
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).absolute()
    assert executable.is_file(), f"Supply a built_in_pgens executable: {executable}"
    source = REPO_ROOT / "inputs/tests/cgl_lf_decay.athinput"
    contents = source.read_text().split("<output1>", 1)[0]
    contents = contents.replace("passive = false\n", "")
    if not lf:
        contents = contents.replace("cgl_heat_flux = landau_fluid\n", "")
        contents = contents.replace("sts_integrator = rkl2", "sts_integrator = none")
    # File blocks can add parameters; command-line overrides require existing keys.
    contents += "\n<mhd>\nlf_coefficient_mode = background\nlf_c_parallel0 = 1\n"
    contents += "\n".join(parameters) + "\n"
    input_file = tmp_path / "startup.athinput"
    input_file.write_text(contents)
    result = subprocess.run(
        [str(executable), "-i", str(input_file), "time/nlim=0",
         "mesh/nx1=8", "meshblock/nx1=8"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    output = result.stdout + result.stderr
    (tmp_path / "startup.log").write_text(output)
    if error:
        assert result.returncode != 0, output
        assert "FATAL ERROR" in output, output
        if isinstance(error, tuple):
            for marker in error:
                assert marker in output, output
        elif error == "passive":
            assert f"<mhd>/{error}" in output, output
            assert "not supported" in output, output
        else:
            assert f"<mhd>/{error}" in output, output
            assert "must be finite" in output, output
    else:
        assert result.returncode == 0, output
        assert "Terminating on cycle limit" in output, output
        assert "cycle=0" in output, output
