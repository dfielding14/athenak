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
        'std::cout << "CGL constructor accepted" << std::endl; }\n'
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
