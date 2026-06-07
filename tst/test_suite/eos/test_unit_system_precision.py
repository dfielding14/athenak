"""Check primitive-solver unit conversions in both Real precisions."""

from pathlib import Path
import os
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
UNIT_DIR = ROOT / "src/eos/primitive-solver"

ATHENA_STUB = r"""
#ifndef ATHENA_HPP_
#define ATHENA_HPP_
#if SINGLE_PRECISION_ENABLED
using Real = float;
#else
using Real = double;
#endif
#define KOKKOS_INLINE_FUNCTION
#endif
"""

PROBE = r"""
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

#include "unit_system.hpp"

using Primitive::UnitSystem;

static_assert(std::is_same_v<UnitSystem::ConversionReal, double>);
static_assert(std::is_same_v<decltype(UnitSystem::energy), double>);
static_assert(std::is_same_v<
    decltype(std::declval<const UnitSystem&>().EnergyConversion(
        std::declval<const UnitSystem&>())), double>);

int main() {
  const std::array<UnitSystem, 5> systems{
      Primitive::MakeCGS(),
      Primitive::MakeGeometricKilometer(),
      Primitive::MakeGeometricSolar(),
      Primitive::MakeNuclear(),
      Primitive::MakeMKS(),
  };

  for (const auto& units : systems) {
    const std::array<double, 13> scales{
        units.c, units.G, units.kb, units.Msun, units.MeV,
        units.length, units.time, units.density, units.mass, units.energy,
        units.pressure, units.temperature, units.chemicalPotential,
    };
    for (const double scale : scales) {
      if (!std::isfinite(scale) || scale == 0.0) return 1;
    }
  }

  for (const auto& from : systems) {
    for (const auto& to : systems) {
      const std::array<double, 12> conversions{
          from.LengthConversion(to),
          from.TimeConversion(to),
          from.VelocityConversion(to),
          from.DensityConversion(to),
          from.MassConversion(to),
          from.MassDensityConversion(to),
          from.EnergyConversion(to),
          from.EnergyDensityConversion(to),
          from.EntropyConversion(to),
          from.PressureConversion(to),
          from.TemperatureConversion(to),
          from.ChemicalPotentialConversion(to),
      };
      for (const double conversion : conversions) {
        if (!std::isfinite(conversion) || conversion == 0.0) return 2;
      }
    }
  }

  const auto geometric_km = Primitive::MakeGeometricKilometer();
  const auto geometric_solar = Primitive::MakeGeometricSolar();
  if (geometric_km.energy == 0.0 || geometric_km.temperature == 0.0 ||
      geometric_solar.energy == 0.0) {
    return 3;
  }
  return 0;
}
"""


@pytest.mark.parametrize("single_precision", [False, True])
def test_unit_system_conversion_range(tmp_path, single_precision):
    """Compile the real implementation and exercise every built-in conversion."""
    compiler = shutil.which(os.environ.get("CXX", "c++"))
    if compiler is None:
        pytest.fail("A C++ compiler is required for the unit-system precision test")

    stub_dir = tmp_path / "include"
    stub_dir.mkdir()
    (stub_dir / "athena.hpp").write_text(ATHENA_STUB)
    probe = tmp_path / "unit_system_probe.cpp"
    probe.write_text(PROBE)
    executable = tmp_path / "unit_system_probe"

    command = [
        compiler,
        "-std=c++17",
        f"-DSINGLE_PRECISION_ENABLED={int(single_precision)}",
        "-I",
        str(stub_dir),
        "-I",
        str(UNIT_DIR),
        str(probe),
        str(UNIT_DIR / "unit_system.cpp"),
        "-o",
        str(executable),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)
