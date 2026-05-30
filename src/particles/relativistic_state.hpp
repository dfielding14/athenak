#ifndef PARTICLES_RELATIVISTIC_STATE_HPP_
#define PARTICLES_RELATIVISTIC_STATE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file relativistic_state.hpp
//! \brief Device-capable helpers for passive relativistic CR tracer state.

#include <algorithm>
#include <limits>

#include "athena.hpp"

namespace particles {
namespace relativistic {

struct Vector3 {
  Real x, y, z;
};

enum class StateStatus : int {
  success = 0,
  invalid_c_model = 1,
  nonfinite_input = 2,
  nonsubluminal_velocity = 3,
  nonfinite_result = 4,
  velocity_shadow_unrepresentable = 5
};

KOKKOS_INLINE_FUNCTION
bool IsFinite(const Real value) {
  return Kokkos::isfinite(value);
}

KOKKOS_INLINE_FUNCTION
bool IsFinite(const Vector3 &value) {
  return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
}

KOKKOS_INLINE_FUNCTION
Real MaxVelocityShadowGamma() {
  Real epsilon = std::numeric_limits<Real>::epsilon();
  return 1.0/Kokkos::sqrt(8.0*epsilon);
}

KOKKOS_INLINE_FUNCTION
Real MaxVelocityParsingGamma() {
  Real epsilon = std::numeric_limits<Real>::epsilon();
  return 0.5/Kokkos::sqrt(Kokkos::sqrt(epsilon));
}

KOKKOS_INLINE_FUNCTION
Real MaxVelocityParsingBeta() {
  Real gamma = MaxVelocityParsingGamma();
  return Kokkos::sqrt((gamma - 1.0)*(gamma + 1.0))/gamma;
}

KOKKOS_INLINE_FUNCTION
Real ScaledNorm3(const Vector3 &value) {
  Real scale = Kokkos::max(Kokkos::fabs(value.x), Kokkos::fabs(value.y));
  scale = Kokkos::max(scale, Kokkos::fabs(value.z));
  if (scale == 0.0) {return 0.0;}
  Real sx = value.x/scale;
  Real sy = value.y/scale;
  Real sz = value.z/scale;
  return scale*Kokkos::sqrt(sx*sx + sy*sy + sz*sz);
}

KOKKOS_INLINE_FUNCTION
StateStatus GammaFromW(const Vector3 &w, const Real c_model, Real &gamma) {
  if (!IsFinite(c_model) || c_model <= 0.0) {
    return StateStatus::invalid_c_model;
  }
  if (!IsFinite(w)) {return StateStatus::nonfinite_input;}

  Real scale = Kokkos::max(c_model, Kokkos::fabs(w.x));
  scale = Kokkos::max(scale, Kokkos::fabs(w.y));
  scale = Kokkos::max(scale, Kokkos::fabs(w.z));
  Real cs = c_model/scale;
  Real wxs = w.x/scale;
  Real wys = w.y/scale;
  Real wzs = w.z/scale;
  gamma = (scale/c_model)*Kokkos::sqrt(cs*cs + wxs*wxs + wys*wys + wzs*wzs);
  if (!IsFinite(gamma) || gamma < 1.0) {
    return StateStatus::nonfinite_result;
  }
  return StateStatus::success;
}

KOKKOS_INLINE_FUNCTION
StateStatus VelocityFromW(const Vector3 &w, const Real c_model, Vector3 &velocity,
                          Real &gamma) {
  StateStatus status = GammaFromW(w, c_model, gamma);
  if (status != StateStatus::success) {return status;}
  if (gamma > MaxVelocityShadowGamma()) {
    return StateStatus::velocity_shadow_unrepresentable;
  }

  Real scale = Kokkos::max(c_model, Kokkos::fabs(w.x));
  scale = Kokkos::max(scale, Kokkos::fabs(w.y));
  scale = Kokkos::max(scale, Kokkos::fabs(w.z));
  Real cs = c_model/scale;
  Real wxs = w.x/scale;
  Real wys = w.y/scale;
  Real wzs = w.z/scale;
  Real scaled_four_norm = Kokkos::sqrt(cs*cs + wxs*wxs + wys*wys + wzs*wzs);
  Real velocity_scale = cs/scaled_four_norm;
  velocity = {w.x*velocity_scale, w.y*velocity_scale, w.z*velocity_scale};
  if (!IsFinite(velocity)) {return StateStatus::nonfinite_result;}

  Real speed = ScaledNorm3(velocity);
  if (!IsFinite(speed)) {return StateStatus::nonfinite_result;}
  if (!(speed < c_model)) {return StateStatus::nonsubluminal_velocity;}
  return StateStatus::success;
}

KOKKOS_INLINE_FUNCTION
StateStatus KineticEnergyFromW(const Vector3 &w, const Real c_model, Real &energy) {
  Real gamma = 0.0;
  StateStatus status = GammaFromW(w, c_model, gamma);
  if (status != StateStatus::success) {return status;}

  Real momentum = ScaledNorm3(w);
  energy = momentum*(momentum/(gamma + 1.0));
  if (!IsFinite(energy) || energy < 0.0) {
    return StateStatus::nonfinite_result;
  }
  return StateStatus::success;
}

KOKKOS_INLINE_FUNCTION
StateStatus WFromVelocity(const Vector3 &velocity, const Real c_model, Vector3 &w,
                          Real &gamma) {
  if (!IsFinite(c_model) || c_model <= 0.0) {
    return StateStatus::invalid_c_model;
  }
  if (!IsFinite(velocity)) {return StateStatus::nonfinite_input;}

  Real speed = ScaledNorm3(velocity);
  if (!IsFinite(speed)) {return StateStatus::nonfinite_result;}
  Real beta = speed/c_model;
  if (!IsFinite(beta)) {return StateStatus::nonfinite_result;}
  if (!(beta < 1.0)) {return StateStatus::nonsubluminal_velocity;}
  if (!(beta < MaxVelocityParsingBeta())) {
    return StateStatus::velocity_shadow_unrepresentable;
  }

  Real inv_gamma = Kokkos::sqrt((1.0 - beta)*(1.0 + beta));
  if (!IsFinite(inv_gamma) || inv_gamma <= 0.0) {
    return StateStatus::nonfinite_result;
  }
  gamma = 1.0/inv_gamma;
  if (!IsFinite(gamma) || gamma < 1.0) {
    return StateStatus::nonfinite_result;
  }
  if (gamma >= MaxVelocityParsingGamma()) {
    return StateStatus::velocity_shadow_unrepresentable;
  }
  w = {velocity.x/inv_gamma, velocity.y/inv_gamma, velocity.z/inv_gamma};
  if (!IsFinite(w)) {
    return StateStatus::nonfinite_result;
  }
  return StateStatus::success;
}

KOKKOS_INLINE_FUNCTION
StateStatus ValidateWState(const Vector3 &w, const Real c_model) {
  Real gamma = 0.0;
  Vector3 velocity{0.0, 0.0, 0.0};
  return VelocityFromW(w, c_model, velocity, gamma);
}

} // namespace relativistic
} // namespace particles

#endif // PARTICLES_RELATIVISTIC_STATE_HPP_
