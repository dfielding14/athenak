//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_relativistic_pusher_kernel_test.cpp
//! \brief Direct host/device and analytical tests for the relativistic CR pusher kernel.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "particles/relativistic_pusher.hpp"
#include "pgen/pgen.hpp"

namespace {

using particles::relativistic::Add;
using particles::relativistic::Cross;
using particles::relativistic::HigueraCaryPush;
using particles::relativistic::IsFinite;
using particles::relativistic::MaxVelocityShadowGamma;
using particles::relativistic::PushResult;
using particles::relativistic::PushStatus;
using particles::relativistic::Scale;
using particles::relativistic::Vector3;

struct ProbeInput {
  Vector3 w;
  Vector3 cE;
  Vector3 b;
  Real alpha_s;
  Real c_model;
  Real dt;
};

struct ProbeResult {
  PushStatus status;
  PushResult push;
};

KOKKOS_INLINE_FUNCTION
ProbeResult EvaluateProbe(const ProbeInput &input) {
  ProbeResult result{
    PushStatus::invalid_input,
    {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0.0, false}
  };
  result.status = HigueraCaryPush(input.w, input.cE, input.b, input.alpha_s,
                                 input.c_model, input.dt, result.push);
  return result;
}

[[noreturn]] void Fail(const std::string &message) {
  std::cout << "CR relativistic pusher kernel test failed: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string FormatReal(const Real value) {
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
  return stream.str();
}

std::string FormatVector(const Vector3 &value) {
  return "{" + FormatReal(value.x) + ", " + FormatReal(value.y) + ", " +
         FormatReal(value.z) + "}";
}

std::string PushStatusName(const PushStatus status) {
  if (status == PushStatus::success) {return "success";}
  if (status == PushStatus::invalid_input) {return "invalid_input";}
  if (status == PushStatus::invalid_state) {return "invalid_state";}
  if (status == PushStatus::arithmetic_range) {return "arithmetic_range";}
  if (status == PushStatus::invalid_rotation) {return "invalid_rotation";}
  if (status == PushStatus::invalid_result) {return "invalid_result";}
  return "unknown(" + std::to_string(static_cast<int>(status)) + ")";
}

bool NearlyEqual(const Real actual, const Real expected, const Real tolerance) {
  if (!std::isfinite(actual) || !std::isfinite(expected)) {return actual == expected;}
  Real scale = std::max(static_cast<Real>(1.0),
                        std::max(std::fabs(actual), std::fabs(expected)));
  return std::fabs(actual - expected) <= tolerance*scale;
}

void Require(const bool condition, const std::string &message) {
  if (!condition) {Fail(message);}
}

void RequireNear(const Real actual, const Real expected, const Real tolerance,
                 const std::string &message) {
  Require(NearlyEqual(actual, expected, tolerance),
          message + ": actual=" + FormatReal(actual) +
          ", expected=" + FormatReal(expected) +
          ", tolerance=" + FormatReal(tolerance));
}

void RequireVectorNear(const Vector3 &actual, const Vector3 &expected,
                       const Real tolerance, const std::string &message) {
  Require(NearlyEqual(actual.x, expected.x, tolerance) &&
          NearlyEqual(actual.y, expected.y, tolerance) &&
          NearlyEqual(actual.z, expected.z, tolerance),
          message + ": actual=" + FormatVector(actual) +
          ", expected=" + FormatVector(expected) +
          ", tolerance=" + FormatReal(tolerance));
}

void RequireStatus(const PushStatus actual, const PushStatus expected,
                   const std::string &message) {
  Require(actual == expected, message + ": actual=" + PushStatusName(actual) +
          ", expected=" + PushStatusName(expected));
}

PushResult RequireSuccessfulPush(const ProbeInput &input, const std::string &message) {
  ProbeResult probe = EvaluateProbe(input);
  RequireStatus(probe.status, PushStatus::success, message);
  return probe.push;
}

Real AnalyticalGamma(const Vector3 &w, const Real c_model) {
  Real wx = w.x/c_model;
  Real wy = w.y/c_model;
  Real wz = w.z/c_model;
  return std::sqrt(1.0 + wx*wx + wy*wy + wz*wz);
}

Vector3 ClassicalBorisPush(const ProbeInput &input) {
  Vector3 epsilon = Scale(0.5*input.dt*input.alpha_s, input.cE);
  Vector3 t = Scale(0.5*input.dt*input.alpha_s, input.b);
  Vector3 w_minus = Add(input.w, epsilon);
  Real t2 = t.x*t.x + t.y*t.y + t.z*t.z;
  Vector3 s = Scale(2.0/(1.0 + t2), t);
  Vector3 w_prime = Add(w_minus, Cross(w_minus, t));
  return Add(Add(w_minus, Cross(w_prime, s)), epsilon);
}

void CheckHostDeviceParity(const std::vector<ProbeInput> &inputs,
                           const std::vector<std::string> &labels,
                           const Real tolerance) {
  Require(inputs.size() == labels.size(), "host/device probe label count");
  int nprobe = static_cast<int>(inputs.size());
  DvceArray1D<ProbeInput> device_inputs("cr_rel_push_probe_inputs", nprobe);
  DvceArray1D<ProbeResult> device_results("cr_rel_push_probe_results", nprobe);
  auto host_inputs = Kokkos::create_mirror_view(device_inputs);
  for (int n=0; n<nprobe; ++n) {host_inputs(n) = inputs[n];}
  Kokkos::deep_copy(device_inputs, host_inputs);

  Kokkos::parallel_for("cr_rel_push_probe",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nprobe),
  KOKKOS_LAMBDA(const int n) {
    device_results(n) = EvaluateProbe(device_inputs(n));
  });
  auto host_results = Kokkos::create_mirror_view_and_copy(HostMemSpace(), device_results);

  for (int n=0; n<nprobe; ++n) {
    ProbeResult expected = EvaluateProbe(inputs[n]);
    ProbeResult actual = host_results(n);
    std::string prefix = "host/device " + labels[n] + " probe " + std::to_string(n);
    RequireStatus(actual.status, expected.status, prefix + " status");
    if (actual.status != PushStatus::success) {continue;}
    RequireVectorNear(actual.push.w, expected.push.w, tolerance, prefix + " momentum");
    RequireVectorNear(actual.push.velocity, expected.push.velocity, tolerance,
                      prefix + " velocity");
    RequireNear(actual.push.gamma, expected.push.gamma, tolerance, prefix + " gamma");
    RequireNear(actual.push.work, expected.push.work, tolerance, prefix + " work");
    Require(actual.push.used_conjugate_root == expected.push.used_conjugate_root,
            prefix + " conjugate-root flag mismatch");
  }
}

void RunKernelTests() {
  const Real epsilon = std::numeric_limits<Real>::epsilon();
  const Real tolerance = 512.0*epsilon;
  const Real low_speed_tolerance = 1024.0*epsilon;
  const Real infinity = std::numeric_limits<Real>::infinity();
  const Real quiet_nan = std::numeric_limits<Real>::quiet_NaN();
  const Real max_value = std::numeric_limits<Real>::max();
  const Vector3 zero{0.0, 0.0, 0.0};

  ProbeInput zero_field{{0.7, -0.4, 0.2}, zero, zero, 1.25, 2.0, 0.75};
  PushResult zero_field_result = RequireSuccessfulPush(zero_field, "zero-field identity");
  Real zero_field_gamma = AnalyticalGamma(zero_field.w, zero_field.c_model);
  RequireVectorNear(zero_field_result.w, zero_field.w, tolerance,
                    "zero-field momentum identity");
  RequireVectorNear(zero_field_result.velocity,
                    Scale(1.0/zero_field_gamma, zero_field.w), tolerance,
                    "zero-field velocity identity");
  RequireNear(zero_field_result.gamma, zero_field_gamma, tolerance,
              "zero-field gamma identity");
  RequireNear(zero_field_result.work, 0.0, tolerance, "zero-field work identity");
  Require(!zero_field_result.used_conjugate_root,
          "zero-field push unexpectedly used conjugate root");

  ProbeInput electric_1d{{0.4, 0.0, 0.0}, {0.75, 0.0, 0.0}, zero, 1.25, 2.0, 0.4};
  PushResult electric_1d_result =
      RequireSuccessfulPush(electric_1d, "aligned uniform-electric-field push");
  Vector3 electric_1d_expected_w{
    electric_1d.w.x + electric_1d.alpha_s*electric_1d.dt*electric_1d.cE.x, 0.0, 0.0};
  Real electric_1d_expected_gamma =
      AnalyticalGamma(electric_1d_expected_w, electric_1d.c_model);
  Real electric_1d_expected_work = electric_1d.c_model*electric_1d.c_model*
      (electric_1d_expected_gamma - AnalyticalGamma(electric_1d.w, electric_1d.c_model));
  RequireVectorNear(electric_1d_result.w, electric_1d_expected_w, tolerance,
                    "aligned uniform-electric-field closed-form momentum");
  RequireNear(electric_1d_result.work, electric_1d_expected_work, tolerance,
              "aligned uniform-electric-field closed-form work");

  ProbeInput magnetic{{0.8, -0.4, 0.25}, zero, {0.2, -0.3, 1.1}, 0.7, 1.0, 0.35};
  PushResult magnetic_forward = RequireSuccessfulPush(magnetic, "magnetic-only push");
  RequireNear(magnetic_forward.gamma, AnalyticalGamma(magnetic.w, magnetic.c_model),
              tolerance, "magnetic-only gamma conservation");
  RequireNear(magnetic_forward.work, 0.0, tolerance, "magnetic-only zero work");
  ProbeInput magnetic_reverse{
    magnetic_forward.w, magnetic.cE, magnetic.b, magnetic.alpha_s,
    magnetic.c_model, -magnetic.dt};
  PushResult magnetic_reversed =
      RequireSuccessfulPush(magnetic_reverse, "magnetic-only reverse push");
  RequireVectorNear(magnetic_reversed.w, magnetic.w, tolerance,
                    "magnetic-only reversibility");
  RequireNear(magnetic_reversed.gamma, AnalyticalGamma(magnetic.w, magnetic.c_model),
              tolerance, "magnetic-only reversed gamma");

  ProbeInput positive_alpha{{0.5, 0.0, 0.0}, {0.4, 0.0, 0.0}, zero, 2.0, 1.0, 0.25};
  ProbeInput negative_alpha = positive_alpha;
  negative_alpha.alpha_s = -positive_alpha.alpha_s;
  PushResult positive_alpha_result =
      RequireSuccessfulPush(positive_alpha, "positive-alpha push");
  PushResult negative_alpha_result =
      RequireSuccessfulPush(negative_alpha, "negative-alpha push");
  RequireVectorNear(positive_alpha_result.w, {0.7, 0.0, 0.0}, tolerance,
                    "positive-alpha momentum direction");
  RequireVectorNear(negative_alpha_result.w, {0.3, 0.0, 0.0}, tolerance,
                    "negative-alpha momentum direction");
  Require(positive_alpha_result.work > 0.0, "positive-alpha push did not gain energy");
  Require(negative_alpha_result.work < 0.0, "negative-alpha push did not lose energy");

  ProbeInput signed_general{
    {0.25, -0.15, 0.05}, {0.3, 0.1, -0.2}, {0.2, -0.4, 0.5}, -0.8, 1.0, 0.3};
  PushResult signed_general_result =
      RequireSuccessfulPush(signed_general, "negative-alpha general push");
  ProbeInput flipped_fields{
    signed_general.w, Scale(-1.0, signed_general.cE), Scale(-1.0, signed_general.b),
    -signed_general.alpha_s, signed_general.c_model, signed_general.dt};
  PushResult flipped_fields_result =
      RequireSuccessfulPush(flipped_fields, "field-flipped positive-alpha general push");
  RequireVectorNear(signed_general_result.w, flipped_fields_result.w, tolerance,
                    "negative-alpha field-sign equivalence momentum");
  RequireNear(signed_general_result.work, flipped_fields_result.work, tolerance,
              "negative-alpha field-sign equivalence work");

  ProbeInput low_speed{
    {1.0e-5, -2.0e-5, 0.5e-5}, {2.0e-6, 1.0e-6, -1.0e-6},
    {0.4, -0.2, 0.3}, 0.8, 1.0, 0.25};
  PushResult low_speed_result = RequireSuccessfulPush(low_speed, "low-speed push");
  RequireVectorNear(low_speed_result.w, ClassicalBorisPush(low_speed),
                    low_speed_tolerance, "low-speed classical Boris limit");
  RequireVectorNear(low_speed_result.velocity, low_speed_result.w,
                    low_speed_tolerance, "low-speed velocity-momentum limit");

  Real high_gamma_scale = 0.25*MaxVelocityShadowGamma();
  ProbeInput high_gamma{
    {high_gamma_scale, static_cast<Real>(-0.25)*high_gamma_scale,
     static_cast<Real>(0.125)*high_gamma_scale},
    {0.01, 0.02, -0.01}, {0.1, -0.2, 0.3}, 0.75, 1.0, 0.125};
  PushResult high_gamma_result = RequireSuccessfulPush(high_gamma, "high-gamma push");
  Require(high_gamma_result.gamma > 100.0,
          "high-gamma probe did not exercise a relativistic state: gamma=" +
          FormatReal(high_gamma_result.gamma));
  Require(IsFinite(high_gamma_result.w) && IsFinite(high_gamma_result.velocity) &&
          std::isfinite(high_gamma_result.gamma) && std::isfinite(high_gamma_result.work),
          "high-gamma push produced a nonfinite result");

  ProbeInput conjugate_root{{0.2, 0.0, 0.0}, zero, {0.0, 0.0, 3.0}, 1.0, 1.0, 1.0};
  PushResult conjugate_root_result =
      RequireSuccessfulPush(conjugate_root, "conjugate-root push");
  Require(conjugate_root_result.used_conjugate_root,
          "conjugate-root probe did not enter the sigma < 0 branch");
  RequireVectorNear(conjugate_root_result.w,
                    {-0.07586942058399232, -0.18505088765053057, 0.0},
                    tolerance, "conjugate-root independently evaluated momentum");
  RequireNear(conjugate_root_result.gamma,
              AnalyticalGamma(conjugate_root.w, conjugate_root.c_model), tolerance,
              "conjugate-root magnetic gamma conservation");

  Real rejected_output_scale = 2.0*MaxVelocityShadowGamma();
  std::vector<ProbeInput> rejected = {
    {{0.1, 0.0, 0.0}, zero, zero, 0.0, 1.0, 0.1},
    {{0.1, 0.0, 0.0}, zero, zero, quiet_nan, 1.0, 0.1},
    {{quiet_nan, 0.0, 0.0}, zero, zero, 1.0, 1.0, 0.1},
    {{0.1, 0.0, 0.0}, {infinity, 0.0, 0.0}, zero, 1.0, 1.0, 0.1},
    {{0.1, 0.0, 0.0}, zero, {0.0, infinity, 0.0}, 1.0, 1.0, 0.1},
    {{0.1, 0.0, 0.0}, zero, zero, 1.0, 1.0, infinity},
    {{0.1, 0.0, 0.0}, zero, zero, 1.0, 0.0, 0.1},
    {{0.1, 0.0, 0.0}, {max_value, 0.0, 0.0}, zero, 4.0, 1.0, 1.0},
    {{0.1, 0.0, 0.0}, zero, {max_value, 0.0, 0.0}, 4.0, 1.0, 1.0},
    {zero, {rejected_output_scale, 0.0, 0.0}, zero, 1.0, 1.0, 1.0},
  };
  std::vector<PushStatus> rejected_statuses = {
    PushStatus::invalid_input,
    PushStatus::invalid_input,
    PushStatus::invalid_input,
    PushStatus::invalid_input,
    PushStatus::invalid_input,
    PushStatus::invalid_input,
    PushStatus::invalid_state,
    PushStatus::arithmetic_range,
    PushStatus::arithmetic_range,
    PushStatus::invalid_result,
  };
  std::vector<std::string> rejected_labels = {
    "zero alpha",
    "NaN alpha",
    "NaN momentum",
    "infinite electric field",
    "infinite magnetic field",
    "infinite timestep",
    "zero model light speed",
    "overflowing electric half-kick",
    "overflowing magnetic rotation",
    "unrepresentable output velocity shadow",
  };
  for (std::size_t n=0; n<rejected.size(); ++n) {
    RequireStatus(EvaluateProbe(rejected[n]).status, rejected_statuses[n],
                  rejected_labels[n] + " rejection");
  }

  std::vector<ProbeInput> parity_inputs = {
    zero_field,
    electric_1d,
    magnetic,
    magnetic_reverse,
    positive_alpha,
    negative_alpha,
    signed_general,
    flipped_fields,
    low_speed,
    high_gamma,
    conjugate_root,
  };
  std::vector<std::string> parity_labels = {
    "zero-field identity",
    "uniform-electric closed form",
    "magnetic-only forward",
    "magnetic-only reverse",
    "positive alpha",
    "negative alpha",
    "negative-alpha general",
    "field-flipped positive-alpha general",
    "low-speed limit",
    "high-gamma stability",
    "conjugate-root branch",
  };
  parity_inputs.insert(parity_inputs.end(), rejected.begin(), rejected.end());
  parity_labels.insert(parity_labels.end(), rejected_labels.begin(),
                       rejected_labels.end());
  CheckHostDeviceParity(parity_inputs, parity_labels, tolerance);
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void)pin;
  (void)restart;
  RunKernelTests();
  std::cout << "CR relativistic pusher kernel tests passed" << std::endl;
}
