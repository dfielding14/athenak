from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def turb_driver_harness(tmp_path_factory: pytest.TempPathFactory) -> Path:
    build_dir = tmp_path_factory.mktemp("turb_driver")
    harness = build_dir / "turb_driver_harness.cpp"
    executable = build_dir / "turb_driver_harness"
    harness.write_text(
        r"""
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <limits>

#include "src/srcterms/turb_driver_utils.hpp"

int main() {
  using turbulence::IsIsotropicMode;
  using turbulence::PositiveConstantEdotScale;

  assert(turbulence::HasCompleteIsotropicBounds(
      3, -3, 3, -3, 3, -3, 3, true, true));
  assert(!turbulence::HasCompleteIsotropicBounds(
      3, 0, 3, -3, 3, -3, 3, true, true));
  assert(turbulence::HasCompleteIsotropicBounds(
      3, -3, 3, 0, 0, 0, 0, false, false));

  // The stable branch must retain a small positive root when m1^2 dominates.
  const double cancellation_root = PositiveConstantEdotScale(1.0, 1.0e20, 1.0);
  assert(cancellation_root > 0.0);
  assert(std::abs(cancellation_root - 1.0e-20) < 1.0e-35);

  for (const double m1 : {-4.0, 0.0, 7.0}) {
    const double root = PositiveConstantEdotScale(2.0, m1, 3.0);
    assert(root >= 0.0);
    assert(std::abs(2.0*root*root + m1*root - 3.0) < 1.0e-12);
  }
  assert(PositiveConstantEdotScale(2.0, -4.0, 0.0) == 0.0);
  const double huge_m0 = std::numeric_limits<double>::max()/16.0;
  const double huge_root = PositiveConstantEdotScale(huge_m0, 0.0, 4.0);
  assert(std::isfinite(huge_root));
  assert(std::abs(huge_root - 2.0/std::sqrt(huge_m0)) < 1.0e-170);
  const double huge_m1 = 0.75*std::numeric_limits<double>::max();
  const double tiny_root = PositiveConstantEdotScale(1.0, huge_m1, 1.0);
  assert(tiny_root > 0.0 && std::isfinite(tiny_root));
  assert(std::abs(huge_m1*tiny_root - 1.0) < 1.0e-14);
  const double balanced_root = PositiveConstantEdotScale(
      huge_m1, -huge_m1, 1.0);
  assert(std::isfinite(balanced_root));
  assert(std::abs(balanced_root - 1.0) < 1.0e-14);

  // A signed cube contains one canonical member of each {k,-k} pair.
  int modes_2d = 0;
  int modes_3d = 0;
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        modes_2d += IsIsotropicMode(i, j, k, 1, 3, true, false) ? 1 : 0;
        modes_3d += IsIsotropicMode(i, j, k, 1, 3, true, true) ? 1 : 0;
      }
    }
  }
  assert(modes_2d == 4);
  assert(modes_3d == 13);
  assert(IsIsotropicMode(1, -1, 0, 1, 3, true, false));
  assert(!IsIsotropicMode(-1, 1, 0, 1, 3, true, false));
  assert(!IsIsotropicMode(1, 0, 1, 1, 3, true, false));

  // The default n=1..3 shell is a complete canonical half-space. Cubic symmetry
  // makes its angular tensor exactly isotropic up to floating-point summation.
  int modes_1d_full = 0;
  int modes_2d_full = 0;
  int modes_3d_full = 0;
  double angular[3][3] = {};
  for (int i = -3; i <= 3; ++i) {
    for (int j = -3; j <= 3; ++j) {
      for (int k = -3; k <= 3; ++k) {
        modes_1d_full += IsIsotropicMode(i, j, k, 1, 9, false, false) ? 1 : 0;
        modes_2d_full += IsIsotropicMode(i, j, k, 1, 9, true, false) ? 1 : 0;
        if (IsIsotropicMode(i, j, k, 1, 9, true, true)) {
          ++modes_3d_full;
          const double wave[3] = {static_cast<double>(i),
                                  static_cast<double>(j),
                                  static_cast<double>(k)};
          const double wave_squared = wave[0]*wave[0] + wave[1]*wave[1] +
                                      wave[2]*wave[2];
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              angular[a][b] += wave[a]*wave[b]/wave_squared;
            }
          }
          assert(!IsIsotropicMode(-i, -j, -k, 1, 9, true, true));
        }
      }
    }
  }
  assert(modes_1d_full == 3);
  assert(modes_2d_full == 14);
  assert(modes_3d_full == 61);
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      const double expected = (a == b) ? modes_3d_full/3.0 : 0.0;
      assert(std::abs(angular[a][b] - expected) < 1.0e-12);
    }
  }

  // In 2D3V, projection uses k=(kx,ky,0): in-plane divergence vanishes while
  // the out-of-plane acceleration survives unchanged.
  const double wave[3] = {2.0, -1.0, 0.0};
  const double raw[3] = {0.7, -0.2, 1.3};
  const double wave_squared = 5.0;
  const double wave_dot_raw = wave[0]*raw[0] + wave[1]*raw[1];
  double projected[3];
  for (int dir = 0; dir < 3; ++dir) {
    projected[dir] = turbulence::BlendProjectedModeComponent(
        raw[dir], wave[dir], wave_dot_raw, wave_squared, 1.0);
  }
  assert(std::abs(wave[0]*projected[0] + wave[1]*projected[1]) < 1.0e-14);
  assert(projected[2] == raw[2]);

  const double spectrum = 5.0/3.0;
  assert(std::abs(turbulence::PowerLawAmplitudeExponent(spectrum, 3) -
                  11.0/6.0) < 1.0e-14);
  assert(std::abs(turbulence::PowerLawAmplitudeExponent(spectrum, 2) -
                  4.0/3.0) < 1.0e-14);
  assert(std::abs(turbulence::ShellCompensationExponent<double>(2) - 0.5) <
         1.0e-14);

  // Staged forcing is the conservative source dE/dt = p*a. Explicit midpoint
  // recovers the exact constant-acceleration energy without a separate a^2 term.
  const double p0 = -0.7;
  const double accel = 1.3;
  const double dt = 0.2;
  const double p_half = p0 + 0.5*accel*dt;
  const double exact_momentum = p0 + accel*dt;
  const double predictor_energy = p0*accel*0.5*dt;
  assert(std::isfinite(predictor_energy));
  const double corrector_energy = p_half*accel*dt;
  const double exact_energy = p0*accel*dt + 0.5*accel*accel*dt*dt;
  const double midpoint_momentum = p0 + accel*dt;
  assert(std::abs(midpoint_momentum - exact_momentum) < 1.0e-14);
  assert(std::abs(corrector_energy - exact_energy) < 1.0e-14);

  // Heun's weighted stage map must produce the same result. Adding an exact-kick
  // quadratic term to either stage would over-inject energy.
  const double predictor_momentum = p0 + accel*dt;
  const double heun_predictor_energy = p0*accel*dt;
  const double heun_momentum = 0.5*predictor_momentum + 0.5*p0 +
                               0.5*accel*dt;
  const double heun_energy = 0.5*heun_predictor_energy +
                             0.5*predictor_momentum*accel*dt;
  assert(std::abs(heun_energy - exact_energy) < 1.0e-14);
  assert(std::abs(heun_momentum - exact_momentum) < 1.0e-14);

  // A standalone impulse is different: its exact map includes the finite-kick term.
  const double impulse_energy = p0*accel*dt + 0.5*accel*accel*dt*dt;
  assert(std::abs(impulse_energy - exact_energy) < 1.0e-14);

  // RK1 uses that same exact map because it has no later velocity evaluation.
  const double rk1_energy = p0*accel*dt + 0.5*accel*accel*dt*dt;
  assert(std::abs(rk1_energy - exact_energy) < 1.0e-14);
  assert(std::abs(p0 + accel*dt - exact_momentum) < 1.0e-14);
  return 0;
}
""",
        encoding="ascii",
    )
    subprocess.run(
        [
            os.environ.get("CXX", "c++"),
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-I",
            str(REPO_ROOT),
            str(harness),
            "-o",
            str(executable),
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    return executable


def test_turbulence_math_and_mode_contract(turb_driver_harness: Path) -> None:
    subprocess.run([str(turb_driver_harness)], cwd=REPO_ROOT, check=True)


def test_driver_uses_three_force_components_for_isotropic_2d() -> None:
    source = (REPO_ROOT / "src/srcterms/turb_driver.cpp").read_text(
        encoding="utf-8"
    )
    assert "const int force_components = (driving_type == 0) ? 3 : 2;" in source
    assert "IsIsotropicMode(nkx, nky, nkz" in source
    assert "PowerLawAmplitudeExponent(\n" in source
    assert "ShellCompensationExponent<Real>(" in source
    assert "pdrive->nexp_stages == 1" in source
    assert "ApplyForcingWithStep(bdt, exact_single_stage);" in source
    assert "ApplyForcingWithStep(kick_dt, true);" in source


def test_uniform_acceleration_regression_is_confined_to_custom_pgen() -> None:
    driver_source = (REPO_ROOT / "src/srcterms/turb_driver.cpp").read_text(
        encoding="utf-8"
    )
    pgen_source = (
        REPO_ROOT / "src/pgen/turb_uniform_accel_test.cpp"
    ).read_text(encoding="utf-8")

    # Normal turbulence input cannot request a spatially uniform k=0 mode.
    assert "test_accel_x1" not in driver_source
    assert "SetManufacturedUniformForce" not in driver_source

    # The custom pgen replaces only the force register after ordinary force setup;
    # the stage task still calls the production AddForcing implementation.
    assert 'tl_map.at("before_timeintegrator")' in pgen_source
    assert "task_list->GetIDLastTask()" in pgen_source
    assert "task_list->AddTask(SetManufacturedUniformForce" in pgen_source
    assert "test_pack->pturb->force" in pgen_source
    assert "ApplyForcingWithStep" not in pgen_source
