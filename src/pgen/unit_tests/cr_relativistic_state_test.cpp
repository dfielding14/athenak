//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_state_test.cpp
//! \brief Host/device parity and limit tests for passive relativistic CR state helpers.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "particles/relativistic_state.hpp"
#include "pgen/pgen.hpp"

namespace {

using particles::relativistic::GammaFromW;
using particles::relativistic::KineticEnergyFromW;
using particles::relativistic::MaxVelocityParsingBeta;
using particles::relativistic::MaxVelocityParsingGamma;
using particles::relativistic::MaxVelocityShadowGamma;
using particles::relativistic::ScaledNorm3;
using particles::relativistic::StateStatus;
using particles::relativistic::ValidateWState;
using particles::relativistic::Vector3;
using particles::relativistic::VelocityFromW;
using particles::relativistic::WFromVelocity;

enum class ProbeOperation : int {
  gamma_from_w = 0,
  velocity_from_w = 1,
  kinetic_energy_from_w = 2,
  w_from_velocity = 3,
  validate_w_state = 4
};

struct ProbeInput {
  ProbeOperation operation;
  Vector3 value;
  Real c_model;
};

struct ProbeResult {
  StateStatus status;
  Vector3 vector;
  Real scalar;
  Real gamma;
};

KOKKOS_INLINE_FUNCTION
ProbeResult EvaluateProbe(const ProbeInput &input) {
  ProbeResult result{StateStatus::success, {0.0, 0.0, 0.0}, 0.0, 0.0};
  if (input.operation == ProbeOperation::gamma_from_w) {
    result.status = GammaFromW(input.value, input.c_model, result.gamma);
  } else if (input.operation == ProbeOperation::velocity_from_w) {
    result.status = VelocityFromW(input.value, input.c_model, result.vector,
                                  result.gamma);
  } else if (input.operation == ProbeOperation::kinetic_energy_from_w) {
    result.status = KineticEnergyFromW(input.value, input.c_model, result.scalar);
  } else if (input.operation == ProbeOperation::w_from_velocity) {
    result.status = WFromVelocity(input.value, input.c_model, result.vector,
                                  result.gamma);
  } else {
    result.status = ValidateWState(input.value, input.c_model);
  }
  return result;
}

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "CR relativistic state helper test failed: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

bool NearlyEqual(const Real a, const Real b, const Real tolerance) {
  if (!std::isfinite(a) || !std::isfinite(b)) {return a == b;}
  Real scale = std::max(static_cast<Real>(1.0), std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= tolerance*scale;
}

void Require(const bool condition, const std::string &message) {
  if (!condition) {Fail(message);}
}

void RequireStatus(const StateStatus actual, const StateStatus expected,
                   const std::string &message) {
  Require(actual == expected, message);
}

void RequireVectorNear(const Vector3 &actual, const Vector3 &expected,
                       const Real tolerance, const std::string &message) {
  Require(NearlyEqual(actual.x, expected.x, tolerance) &&
          NearlyEqual(actual.y, expected.y, tolerance) &&
          NearlyEqual(actual.z, expected.z, tolerance), message);
}

void CheckHostDeviceParity(const std::vector<ProbeInput> &inputs, const Real tolerance) {
  int nprobe = static_cast<int>(inputs.size());
  DvceArray1D<ProbeInput> device_inputs("cr_rel_state_probe_inputs", nprobe);
  DvceArray1D<ProbeResult> device_results("cr_rel_state_probe_results", nprobe);
  auto host_inputs = Kokkos::create_mirror_view(device_inputs);
  for (int n=0; n<nprobe; ++n) {host_inputs(n) = inputs[n];}
  Kokkos::deep_copy(device_inputs, host_inputs);

  Kokkos::parallel_for("cr_rel_state_probe",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nprobe),
  KOKKOS_LAMBDA(const int n) {
    device_results(n) = EvaluateProbe(device_inputs(n));
  });
  auto host_results = Kokkos::create_mirror_view_and_copy(HostMemSpace(), device_results);

  for (int n=0; n<nprobe; ++n) {
    ProbeResult expected = EvaluateProbe(inputs[n]);
    ProbeResult actual = host_results(n);
    Require(actual.status == expected.status,
            "host/device status mismatch for probe " + std::to_string(n));
    RequireVectorNear(actual.vector, expected.vector, tolerance,
                      "host/device vector mismatch for probe " + std::to_string(n));
    Require(NearlyEqual(actual.scalar, expected.scalar, tolerance),
            "host/device scalar mismatch for probe " + std::to_string(n));
    Require(NearlyEqual(actual.gamma, expected.gamma, tolerance),
            "host/device gamma mismatch for probe " + std::to_string(n));
  }
}

void RunHelperTests() {
  const Real epsilon = std::numeric_limits<Real>::epsilon();
  const Real tolerance = 128.0*epsilon;
  const Real infinity = std::numeric_limits<Real>::infinity();
  const Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();
  const Real max_value = std::numeric_limits<Real>::max();
  const Real min_value = std::numeric_limits<Real>::min();
  const Real sqrt_max = std::sqrt(max_value);
  const Real parsing_beta_max = MaxVelocityParsingBeta();
  const Real parsing_gamma_max = MaxVelocityParsingGamma();
  const Real shadow_gamma_max = MaxVelocityShadowGamma();

  Real gamma = 0.0;
  Real energy = 0.0;
  Vector3 velocity{0.0, 0.0, 0.0};
  Vector3 w{0.0, 0.0, 0.0};

  RequireStatus(GammaFromW({0.0, 0.0, 0.0}, 1.0, gamma), StateStatus::success,
                "zero-momentum gamma status");
  Require(NearlyEqual(gamma, 1.0, tolerance), "zero-momentum gamma");
  RequireStatus(VelocityFromW({0.0, 0.0, 0.0}, 1.0, velocity, gamma),
                StateStatus::success, "zero-momentum velocity status");
  RequireVectorNear(velocity, {0.0, 0.0, 0.0}, tolerance, "zero-momentum velocity");
  RequireStatus(KineticEnergyFromW({0.0, 0.0, 0.0}, 1.0, energy),
                StateStatus::success, "zero-momentum kinetic-energy status");
  Require(energy == 0.0, "zero-momentum kinetic energy");

  Vector3 small_w{1.0e-4, -2.0e-4, 3.0e-4};
  Real small_w2 = small_w.x*small_w.x + small_w.y*small_w.y + small_w.z*small_w.z;
  RequireStatus(KineticEnergyFromW(small_w, 1.0, energy), StateStatus::success,
                "nonrelativistic kinetic-energy status");
  Require(NearlyEqual(energy, 0.5*small_w2, 128.0*std::sqrt(epsilon)),
          "nonrelativistic kinetic-energy limit");

  Vector3 moderate_w{3.0, 4.0, 0.0};
  RequireStatus(GammaFromW(moderate_w, 2.0, gamma), StateStatus::success,
                "moderate gamma status");
  Require(NearlyEqual(gamma, std::sqrt(7.25), tolerance), "moderate gamma");
  RequireStatus(VelocityFromW(moderate_w, 2.0, velocity, gamma), StateStatus::success,
                "moderate velocity status");
  Require(ScaledNorm3(velocity) < 2.0, "moderate reconstructed velocity is subluminal");

  Vector3 overflow_prone{0.75*sqrt_max, 0.75*sqrt_max, 0.0};
  Real overflow_norm = ScaledNorm3(overflow_prone);
  Require(std::isfinite(overflow_norm), "scaled norm survives naive-square overflow");
  RequireStatus(GammaFromW(overflow_prone, 1.0, gamma), StateStatus::success,
                "scaled gamma survives naive-square overflow");
  Require(std::isfinite(gamma), "scaled gamma is finite");
  RequireStatus(VelocityFromW(overflow_prone, 1.0, velocity, gamma),
                StateStatus::velocity_shadow_unrepresentable,
                "overflow-prone velocity shadow is rejected deterministically");

  Vector3 high_w{1000.0, 0.0, 0.0};
  RequireStatus(GammaFromW(high_w, 1.0, gamma), StateStatus::success,
                "high-momentum gamma status");
  Require(std::fabs(gamma/1000.0 - 1.0) < 1.0e-6,
          "high-momentum gamma asymptote");
  RequireStatus(KineticEnergyFromW(high_w, 1.0, energy), StateStatus::success,
                "high-momentum kinetic-energy status");
  Require(std::fabs(energy/1000.0 - 1.0) < 2.0e-3,
          "high-momentum kinetic-energy asymptote");
  RequireStatus(VelocityFromW(high_w, 1.0, velocity, gamma), StateStatus::success,
                "high-momentum velocity-shadow status");
  Require(ScaledNorm3(velocity) < 1.0, "high-momentum shadow is subluminal");
  RequireStatus(ValidateWState(high_w, 1.0), StateStatus::success,
                "high-momentum state validation");

  Real below_gamma = 0.9*shadow_gamma_max;
  Real above_gamma = 2.0*shadow_gamma_max;
  Vector3 below_shadow_limit{std::sqrt((below_gamma - 1.0)*(below_gamma + 1.0)),
                             0.0, 0.0};
  Vector3 above_shadow_limit{std::sqrt((above_gamma - 1.0)*(above_gamma + 1.0)),
                             0.0, 0.0};
  RequireStatus(VelocityFromW(below_shadow_limit, 1.0, velocity, gamma),
                StateStatus::success, "velocity shadow below gamma cap");
  Require(ScaledNorm3(velocity) < 1.0, "below-cap shadow is subluminal");
  RequireStatus(ValidateWState(below_shadow_limit, 1.0), StateStatus::success,
                "state validation below gamma cap");
  RequireStatus(VelocityFromW(above_shadow_limit, 1.0, velocity, gamma),
                StateStatus::velocity_shadow_unrepresentable,
                "velocity shadow above gamma cap");
  RequireStatus(ValidateWState(above_shadow_limit, 1.0),
                StateStatus::velocity_shadow_unrepresentable,
                "state validation above gamma cap");

  Real near_below_gamma = (1.0 - 128.0*epsilon)*shadow_gamma_max;
  Real near_above_gamma = (1.0 + 128.0*epsilon)*shadow_gamma_max;
  Vector3 near_below_shadow_limit{
    std::sqrt((near_below_gamma - 1.0)*(near_below_gamma + 1.0)), 0.0, 0.0};
  Vector3 near_above_shadow_limit{
    std::sqrt((near_above_gamma - 1.0)*(near_above_gamma + 1.0)), 0.0, 0.0};
  RequireStatus(VelocityFromW(near_below_shadow_limit, 1.0, velocity, gamma),
                StateStatus::success, "velocity shadow immediately below gamma cap");
  RequireStatus(VelocityFromW(near_above_shadow_limit, 1.0, velocity, gamma),
                StateStatus::velocity_shadow_unrepresentable,
                "velocity shadow immediately above gamma cap");

  Vector3 source_velocity{0.3, -0.4, 0.2};
  RequireStatus(WFromVelocity(source_velocity, 1.0, w, gamma), StateStatus::success,
                "velocity-to-momentum status");
  RequireStatus(VelocityFromW(w, 1.0, velocity, gamma), StateStatus::success,
                "momentum-to-velocity round-trip status");
  RequireVectorNear(velocity, source_velocity, tolerance, "velocity round trip");
  Require(ScaledNorm3(velocity) < 1.0, "round-trip reconstructed velocity is subluminal");

  Vector3 source_w{0.7, -0.2, 0.4};
  RequireStatus(VelocityFromW(source_w, 1.0, velocity, gamma), StateStatus::success,
                "momentum-to-velocity status");
  RequireStatus(WFromVelocity(velocity, 1.0, w, gamma), StateStatus::success,
                "velocity-to-momentum round-trip status");
  RequireVectorNear(w, source_w, tolerance, "momentum round trip");

  RequireStatus(WFromVelocity({0.99, 0.0, 0.0}, 1.0, w, gamma), StateStatus::success,
                "near-limit velocity conversion status");
  RequireStatus(VelocityFromW(w, 1.0, velocity, gamma), StateStatus::success,
                "near-limit momentum conversion status");
  Require(ScaledNorm3(velocity) < 1.0, "near-limit reconstructed velocity is subluminal");

  RequireStatus(GammaFromW({0.0, 0.0, 0.0}, 0.0, gamma),
                StateStatus::invalid_c_model, "zero c_model rejection");
  RequireStatus(GammaFromW({0.0, 0.0, 0.0}, -1.0, gamma),
                StateStatus::invalid_c_model, "negative c_model rejection");
  RequireStatus(GammaFromW({0.0, 0.0, 0.0}, quiet_nan, gamma),
                StateStatus::invalid_c_model, "NaN c_model rejection");
  RequireStatus(GammaFromW({0.0, 0.0, 0.0}, infinity, gamma),
                StateStatus::invalid_c_model, "infinite c_model rejection");
  RequireStatus(GammaFromW({quiet_nan, 0.0, 0.0}, 1.0, gamma),
                StateStatus::nonfinite_input, "NaN momentum rejection");
  RequireStatus(GammaFromW({infinity, 0.0, 0.0}, 1.0, gamma),
                StateStatus::nonfinite_input, "infinite momentum rejection");
  RequireStatus(WFromVelocity({quiet_nan, 0.0, 0.0}, 1.0, w, gamma),
                StateStatus::nonfinite_input, "NaN velocity rejection");
  RequireStatus(WFromVelocity({infinity, 0.0, 0.0}, 1.0, w, gamma),
                StateStatus::nonfinite_input, "infinite velocity rejection");
  RequireStatus(WFromVelocity({0.0, 0.0, 0.0}, 0.0, w, gamma),
                StateStatus::invalid_c_model, "conversion zero c_model rejection");
  RequireStatus(WFromVelocity({0.0, 0.0, 0.0}, -1.0, w, gamma),
                StateStatus::invalid_c_model, "conversion negative c_model rejection");
  RequireStatus(WFromVelocity({0.0, 0.0, 0.0}, quiet_nan, w, gamma),
                StateStatus::invalid_c_model, "conversion NaN c_model rejection");
  RequireStatus(WFromVelocity({0.0, 0.0, 0.0}, infinity, w, gamma),
                StateStatus::invalid_c_model, "conversion infinite c_model rejection");
  RequireStatus(WFromVelocity({1.0, 0.0, 0.0}, 1.0, w, gamma),
                StateStatus::nonsubluminal_velocity, "light-speed velocity rejection");
  RequireStatus(WFromVelocity({1.01, 0.0, 0.0}, 1.0, w, gamma),
                StateStatus::nonsubluminal_velocity, "superluminal velocity rejection");
  Vector3 below_parsing_limit{
    std::nextafter(parsing_beta_max, static_cast<Real>(0.0)), 0.0, 0.0};
  Vector3 at_parsing_limit{parsing_beta_max, 0.0, 0.0};
  Vector3 above_parsing_limit{
    std::nextafter(parsing_beta_max, static_cast<Real>(1.0)), 0.0, 0.0};
  RequireStatus(WFromVelocity(below_parsing_limit, 1.0, w, gamma),
                StateStatus::success, "velocity parsing below gamma cap");
  RequireStatus(WFromVelocity(at_parsing_limit, 1.0, w, gamma),
                StateStatus::velocity_shadow_unrepresentable,
                "velocity parsing at gamma cap");
  RequireStatus(WFromVelocity(above_parsing_limit, 1.0, w, gamma),
                StateStatus::velocity_shadow_unrepresentable,
                "velocity parsing above gamma cap");
  Real round_trip_gamma = 0.5*parsing_gamma_max;
  Vector3 high_gamma_round_trip{
    std::sqrt((round_trip_gamma - 1.0)*(round_trip_gamma + 1.0)), 0.0, 0.0};
  RequireStatus(VelocityFromW(high_gamma_round_trip, 1.0, velocity, gamma),
                StateStatus::success, "high-gamma output-shadow status");
  RequireStatus(WFromVelocity(velocity, 1.0, w, gamma),
                StateStatus::success, "high-gamma inverse parsing status");
  RequireVectorNear(w, high_gamma_round_trip, 64.0*std::sqrt(epsilon),
                    "high-gamma momentum round trip");
  RequireStatus(GammaFromW({max_value, 0.0, 0.0}, min_value, gamma),
                StateStatus::nonfinite_result, "nonfinite gamma-result rejection");

  std::vector<ProbeInput> probes = {
    {ProbeOperation::gamma_from_w, {0.0, 0.0, 0.0}, 1.0},
    {ProbeOperation::velocity_from_w, small_w, 1.0},
    {ProbeOperation::kinetic_energy_from_w, moderate_w, 2.0},
    {ProbeOperation::gamma_from_w, overflow_prone, 1.0},
    {ProbeOperation::velocity_from_w, overflow_prone, 1.0},
    {ProbeOperation::velocity_from_w, high_w, 1.0},
    {ProbeOperation::kinetic_energy_from_w, high_w, 1.0},
    {ProbeOperation::validate_w_state, below_shadow_limit, 1.0},
    {ProbeOperation::validate_w_state, above_shadow_limit, 1.0},
    {ProbeOperation::validate_w_state, near_below_shadow_limit, 1.0},
    {ProbeOperation::validate_w_state, near_above_shadow_limit, 1.0},
    {ProbeOperation::velocity_from_w, source_w, 1.0},
    {ProbeOperation::w_from_velocity, source_velocity, 1.0},
    {ProbeOperation::w_from_velocity, {0.99, 0.0, 0.0}, 1.0},
    {ProbeOperation::w_from_velocity, below_parsing_limit, 1.0},
    {ProbeOperation::w_from_velocity, at_parsing_limit, 1.0},
    {ProbeOperation::w_from_velocity, above_parsing_limit, 1.0},
    {ProbeOperation::validate_w_state, moderate_w, 2.0},
    {ProbeOperation::gamma_from_w, {max_value, 0.0, 0.0}, min_value},
    {ProbeOperation::gamma_from_w, {0.0, 0.0, 0.0}, 0.0},
    {ProbeOperation::gamma_from_w, {quiet_nan, 0.0, 0.0}, 1.0},
    {ProbeOperation::gamma_from_w, {infinity, 0.0, 0.0}, 1.0},
    {ProbeOperation::validate_w_state, {quiet_nan, 0.0, 0.0}, 1.0},
    {ProbeOperation::validate_w_state, {infinity, 0.0, 0.0}, 1.0},
    {ProbeOperation::w_from_velocity, {0.0, 0.0, 0.0}, 0.0},
    {ProbeOperation::w_from_velocity, {0.0, 0.0, 0.0}, -1.0},
    {ProbeOperation::w_from_velocity, {0.0, 0.0, 0.0}, quiet_nan},
    {ProbeOperation::w_from_velocity, {0.0, 0.0, 0.0}, infinity},
    {ProbeOperation::w_from_velocity, {quiet_nan, 0.0, 0.0}, 1.0},
    {ProbeOperation::w_from_velocity, {infinity, 0.0, 0.0}, 1.0},
    {ProbeOperation::w_from_velocity, {1.0, 0.0, 0.0}, 1.0},
    {ProbeOperation::w_from_velocity, {1.01, 0.0, 0.0}, 1.0},
  };
  CheckHostDeviceParity(probes, tolerance);
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void)pin;
  (void)restart;
  RunHelperTests();
  std::cout << "CR relativistic state helper tests passed" << std::endl;
}
