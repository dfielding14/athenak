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
                     "lf_k_parallel", "lf_c_parallel0")
        for value in ("0", "-1", "nan", "inf")
    ] + [
        pytest.param(True, (f"{name}={value}",), name, id=f"{name}-{value}")
        for name in ("tfloor", "sfloor", "nu_coll", "limiter_nu_coll")
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
        pytest.param(True, (), None, id="active-default"),
        pytest.param(False, ("passive=false",), None, id="active-no-lf"),
        pytest.param(True, ("tfloor=0", "sfloor=0", "nu_coll=0",
                           "limiter_nu_coll=0"), None, id="zero-floors-and-rates"),
        pytest.param(True, ("limiter_nu_coll=7",), None, id="disabled-rate"),
        pytest.param(True, ("nu_coll=3",), None, id="background-rate"),
        pytest.param(True, ("mirror_limiter=true",), None,
                     id="mirror-default-rate"),
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
        assert f"<mhd>/{error}" in output, output
        if error == "passive":
            assert "not supported" in output, output
        else:
            assert "must be finite" in output, output
    else:
        assert result.returncode == 0, output
        assert "Terminating on cycle limit" in output, output
        assert "cycle=0" in output, output
