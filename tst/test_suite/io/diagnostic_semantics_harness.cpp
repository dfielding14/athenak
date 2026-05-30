//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file diagnostic_semantics_harness.cpp
//! \brief Direct analytic regression harness for generic output diagnostics.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "outputs/diagnostic_semantics.hpp"

namespace {

[[noreturn]] void Fail(const char *message) {
  std::cerr << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Require(bool condition, const char *message) {
  if (!condition) Fail(message);
}

void RequireClose(double actual, double expected, const char *message) {
  if (std::fabs(actual - expected) > 1.0e-12) Fail(message);
}

void TestGeometry() {
  auto origin = output_diagnostics::BuildGeometry(0.0, 0.0, 0.0);
  Require(origin.valid, "Origin geometry must be valid.");
  RequireClose(origin.theta, 0.0, "Origin theta must be canonical zero.");
  RequireClose(origin.phi, 0.0, "Origin phi must be canonical zero.");
  RequireClose(origin.costheta, 1.0, "Origin costheta must be canonical one.");

  double above_one = std::nextafter(1.0, 2.0);
  RequireClose(output_diagnostics::ClampUnit(above_one), 1.0,
               "Positive floating overshoot was not clamped.");
  RequireClose(output_diagnostics::ClampUnit(-above_one), -1.0,
               "Negative floating overshoot was not clamped.");
  Require(!output_diagnostics::BuildGeometry(
              std::numeric_limits<double>::infinity(), 0.0, 0.0).valid,
          "Nonfinite geometry must be rejected.");
}

void TestFlow() {
  auto flow = output_diagnostics::BuildFlow(1.0, 0.0, 0.0, 2.0, 4.0, -6.0, 8.0);
  Require(flow.valid, "Finite positive-density flow must be valid.");
  RequireClose(flow.radial_velocity, 2.0, "Radial velocity is incorrect.");
  RequireClose(flow.theta_velocity, -4.0, "Theta velocity is incorrect.");
  RequireClose(flow.phi_velocity, -3.0, "Phi velocity is incorrect.");
  RequireClose(flow.radial_mass_flux, 4.0, "Radial mass flux is incorrect.");
  RequireClose(flow.vertical_mass_flux, 0.0,
               "Midplane vertical mass flux must be canonical zero.");

  auto axis = output_diagnostics::BuildFlow(0.0, 0.0, 1.0, 2.0, 4.0, -6.0, 8.0);
  Require(axis.valid, "Axis flow must be valid.");
  RequireClose(axis.phi_velocity, 0.0, "Axis phi velocity must be canonical zero.");
  RequireClose(axis.theta_velocity, 0.0, "Axis theta velocity must be canonical zero.");
  RequireClose(axis.radial_velocity, 4.0, "Axis radial velocity is incorrect.");

  Require(!output_diagnostics::BuildFlow(1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0).valid,
          "Zero density must be rejected.");
  Require(!output_diagnostics::BuildFlow(1.0, 0.0, 0.0, -1.0, 1.0, 2.0, 3.0).valid,
          "Negative density must be rejected.");
  Require(!output_diagnostics::BuildFlow(
              1.0, 0.0, 0.0, 1.0, std::numeric_limits<double>::quiet_NaN(),
              2.0, 3.0).valid,
          "Nonfinite momentum must be rejected.");
}

void TestEnergy() {
  auto hydro = output_diagnostics::BuildEnergyFlux(
      1.0, 0.0, 0.0, 2.0, 4.0, -6.0, 8.0, 41.5, 1.4, false,
      0.0, 0.0, 0.0, true);
  Require(hydro.valid, "Finite Hydro energy state must be valid.");
  RequireClose(hydro.kinetic_radial, 29.0*2.0,
               "Hydro kinetic radial flux is incorrect.");
  RequireClose(hydro.thermal_radial, 17.5*2.0,
               "Hydro thermal radial flux is incorrect.");
  RequireClose(hydro.total_radial, 46.5*2.0,
               "Hydro total radial flux is incorrect.");

  auto mhd = output_diagnostics::BuildEnergyFlux(
      1.0, 0.0, 0.0, 2.0, 4.0, -6.0, 8.0, 44.5, 1.4, true,
      1.0, 2.0, -1.0, true);
  Require(mhd.valid, "Finite MHD energy state must be valid.");
  RequireClose(mhd.magnetic_radial, 20.0,
               "MHD magnetic radial flux is incorrect.");
  RequireClose(mhd.total_radial, 113.0,
               "MHD total radial flux is incorrect.");

  auto kinetic_only = output_diagnostics::BuildEnergyFlux(
      1.0, 0.0, 0.0, 2.0, 4.0, -6.0, 8.0,
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(), false, 0.0, 0.0, 0.0, false);
  Require(kinetic_only.valid,
          "Kinetic-only flux must not require total energy or gamma.");
  RequireClose(kinetic_only.kinetic_radial, 58.0,
               "Kinetic-only radial flux is incorrect.");

  auto magnetic_only = output_diagnostics::BuildEnergyFlux(
      1.0, 0.0, 0.0, 2.0, 4.0, -6.0, 8.0,
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(), true, 1.0, 2.0, -1.0, false);
  Require(magnetic_only.valid,
          "Magnetic-only flux must not require total energy or gamma.");
  RequireClose(magnetic_only.magnetic_radial, 20.0,
               "Magnetic-only radial flux is incorrect.");
}

}  // namespace

int main() {
  TestGeometry();
  TestFlow();
  TestEnergy();
  return EXIT_SUCCESS;
}
