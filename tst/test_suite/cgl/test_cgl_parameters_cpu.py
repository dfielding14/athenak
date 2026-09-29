"""Construct the production CGL EOS in double and single precision."""

from pathlib import Path
import shlex
import subprocess

import pytest

from test_cgl_c2p_pressure_floor_cpu import REPO_ROOT, _cmake_cache_value, _run


@pytest.mark.parametrize("single_precision", [False, True])
def test_cgl_constructor_parameters(tmp_path, single_precision):
    build = tmp_path / "build"
    _run([
        "cmake", "-S", str(REPO_ROOT), "-B", str(build), "-G", "Unix Makefiles",
        "-DAthena_ENABLE_MPI=OFF", "-DKokkos_ENABLE_HIP=OFF",
        "-DKokkos_ENABLE_SERIAL=ON",
        f"-DAthena_SINGLE_PRECISION={'ON' if single_precision else 'OFF'}",
    ])
    sources = (
        "eos/cgl_mhd.cpp", "eos/eos.cpp", "parameter_input.cpp",
        "globals.cpp", "outputs/io_wrapper.cpp",
    )
    objects = [Path("src/CMakeFiles/athena.dir") / f"{s}.o" for s in sources]
    _run([
        "make", "-f", "src/CMakeFiles/athena.dir/build.make",
        *map(str, objects), "-j4",
    ], cwd=build)
    _run(["cmake", "--build", str(build), "--target", "kokkoscore", "-j4"])
    flags = {}
    for line in (build / "src/CMakeFiles/athena.dir/flags.make").read_text().splitlines():
        if " = " in line:
            key, value = line.split(" = ", 1)
            flags[key] = shlex.split(value)
    main = tmp_path / "main.cpp"
    main.write_text(
        '#include <iostream>\n#include "athena.hpp"\n#include "mesh/mesh.hpp"\n'
        '#include "eos/eos.hpp"\n#include "parameter_input.hpp"\n'
        'int main() { ParameterInput pin; pin.LoadFromStream(std::cin);\n'
        'CGLMHD eos(nullptr, &pin);\n'
        'std::cout << "CGL constructor accepted coll=" << eos.eos_data.coll\n'
        ' << " passive=" << eos.eos_data.passive << " mirror=" << eos.eos_data.mlim\n'
        ' << " firehose=" << eos.eos_data.flim << " backup=" << eos.eos_data.backup_lim\n'
        ' << " rate=" << eos.eos_data.lim_coll << std::endl; }\n'
    )
    executable = tmp_path / "cgl_constructor"
    _run([
        _cmake_cache_value(build, "CMAKE_CXX_COMPILER"),
        *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"], *flags["CXX_FLAGS"],
        str(main), *(str(build / obj) for obj in objects),
        str(build / "kokkos/core/src/libkokkoscore.a"), "-o", str(executable),
    ])
    for value, error in (
        ("0", "must be positive"),
        ("-1", "must be positive"),
        (None, "set <mhd>/bfloor explicitly" if single_precision else None),
        ("1e-20", "set <mhd>/bfloor explicitly" if single_precision else None),
        ("1e-10", None),
    ):
        parameters = "<mhd>\n" + (f"bfloor = {value}\n" if value is not None else "")
        result = subprocess.run(
            [str(executable)], input=parameters, capture_output=True,
            text=True, check=False,
        )
        output = result.stdout + result.stderr
        if error:
            assert result.returncode != 0, output
            assert error in output
        else:
            assert result.returncode == 0, output
            assert "CGL constructor accepted" in output
    for parameters, expected, error in (
        ("", "coll=0 passive=0 mirror=0 firehose=0 backup=0 rate=0", None),
        ("mirror_limiter = false\nfirehose_limiter = false\n",
         "coll=0 passive=0 mirror=0 firehose=0 backup=0 rate=0", None),
        ("passive = false\n", "passive=0", None),
        ("limiter_hardwall = false\n", "coll=0", None),
        ("limiter_hardwall = true\n", None,
         "limiter_hardwall=true is no longer supported"),
        ("mirror_limiter = true\nlimiter_nu_coll = 7\nlimiter_hardwall = true\n",
         None, "limiter_nu_coll=1e10"),
        ("passive = true\niso_sound_speed = 1\n", "passive=1", None),
        ("mirror_limiter = true\n", None, "limiter_nu_coll"),
        ("firehose_limiter = true\n", None, "limiter_nu_coll"),
        ("backup_limiters = true\n", None,
         "backup_limiters requires mirror_limiter or firehose_limiter"),
        ("mirror_limiter = false\nbackup_limiters = true\n", None,
         "backup_limiters requires mirror_limiter or firehose_limiter"),
        ("mirror_limiter = true\nlimiter_nu_coll = 3\nbackup_limiters = true\n",
         "coll=1 passive=0 mirror=1 firehose=0 backup=1 rate=3", None),
        ("firehose_limiter = true\nlimiter_nu_coll = 0\n",
         "coll=1 passive=0 mirror=0 firehose=1 backup=0 rate=0", None),
    ):
        result = subprocess.run(
            [str(executable)], input="<mhd>\nbfloor = 1e-10\n" + parameters,
            capture_output=True, text=True, check=False,
        )
        output = result.stdout + result.stderr
        if error:
            assert result.returncode != 0, output
            assert error in output
        else:
            assert result.returncode == 0, output
            assert expected in output
