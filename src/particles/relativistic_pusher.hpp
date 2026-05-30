#ifndef PARTICLES_RELATIVISTIC_PUSHER_HPP_
#define PARTICLES_RELATIVISTIC_PUSHER_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file relativistic_pusher.hpp
//! \brief Device-capable reference Higuera-Cary push for passive relativistic CRs.

#include <limits>

#include "athena.hpp"
#include "relativistic_state.hpp"

namespace particles {
namespace relativistic {

enum class PushStatus : int {
  success = 0,
  invalid_input = 1,
  invalid_state = 2,
  arithmetic_range = 3,
  invalid_rotation = 4,
  invalid_result = 5
};

struct PushResult {
  Vector3 w;
  Vector3 velocity;
  Real gamma;
  Real work;
  bool used_conjugate_root;
};

KOKKOS_INLINE_FUNCTION
Vector3 Cross(const Vector3 &a, const Vector3 &b) {
  return {a.y*b.z - a.z*b.y,
          a.z*b.x - a.x*b.z,
          a.x*b.y - a.y*b.x};
}

KOKKOS_INLINE_FUNCTION
Real Dot(const Vector3 &a, const Vector3 &b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

KOKKOS_INLINE_FUNCTION
Vector3 Add(const Vector3 &a, const Vector3 &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

KOKKOS_INLINE_FUNCTION
Vector3 Scale(const Real scale, const Vector3 &value) {
  return {scale*value.x, scale*value.y, scale*value.z};
}

KOKKOS_INLINE_FUNCTION
bool CheckedProduct(const Real a, const Real b, Real &product) {
  if (!IsFinite(a) || !IsFinite(b)) {return false;}
  if (a != 0.0 && Kokkos::fabs(b) > std::numeric_limits<Real>::max()/Kokkos::fabs(a)) {
    return false;
  }
  product = a*b;
  return IsFinite(product);
}

KOKKOS_INLINE_FUNCTION
bool CheckedScalarAdd(const Real a, const Real b, Real &sum) {
  if (!IsFinite(a) || !IsFinite(b)) {return false;}
  sum = a + b;
  return IsFinite(sum);
}

KOKKOS_INLINE_FUNCTION
bool CheckedScale(const Real scale, const Vector3 &value, Vector3 &result) {
  if (!IsFinite(scale) || !IsFinite(value)) {return false;}
  const Real norm = ScaledNorm3(value);
  if (!IsFinite(norm)) {return false;}
  if (norm != 0.0 &&
      Kokkos::fabs(scale) > std::numeric_limits<Real>::max()/norm) {
    return false;
  }
  result = Scale(scale, value);
  return IsFinite(result);
}

KOKKOS_INLINE_FUNCTION
bool CheckedAdd(const Vector3 &a, const Vector3 &b, Vector3 &result) {
  if (!IsFinite(a) || !IsFinite(b)) {return false;}
  result = Add(a, b);
  return IsFinite(result);
}

KOKKOS_INLINE_FUNCTION
bool CheckedDot(const Vector3 &a, const Vector3 &b, Real &result) {
  if (!IsFinite(a) || !IsFinite(b)) {return false;}
  const Real anorm = ScaledNorm3(a);
  const Real bnorm = ScaledNorm3(b);
  if (!IsFinite(anorm) || !IsFinite(bnorm)) {return false;}
  if (anorm == 0.0 || bnorm == 0.0) {
    result = 0.0;
    return true;
  }
  Real scale = 0.0;
  if (!CheckedProduct(anorm, bnorm, scale)) {return false;}
  const Real normalized = (a.x/anorm)*(b.x/bnorm) +
                          (a.y/anorm)*(b.y/bnorm) +
                          (a.z/anorm)*(b.z/bnorm);
  return CheckedProduct(scale, normalized, result);
}

KOKKOS_INLINE_FUNCTION
bool CheckedCross(const Vector3 &a, const Vector3 &b, Vector3 &result) {
  if (!IsFinite(a) || !IsFinite(b)) {return false;}
  const Real anorm = ScaledNorm3(a);
  const Real bnorm = ScaledNorm3(b);
  if (!IsFinite(anorm) || !IsFinite(bnorm)) {return false;}
  if (anorm == 0.0 || bnorm == 0.0) {
    result = {0.0, 0.0, 0.0};
    return true;
  }
  Real scale = 0.0;
  if (!CheckedProduct(anorm, bnorm, scale)) {return false;}
  return CheckedScale(scale, Cross(Scale(1.0/anorm, a), Scale(1.0/bnorm, b)),
                      result);
}

KOKKOS_INLINE_FUNCTION
Real SafeSquareRootLimit() {
  return 0.25*Kokkos::sqrt(std::numeric_limits<Real>::max());
}

KOKKOS_INLINE_FUNCTION
PushStatus HigueraCaryRotationGamma(const Vector3 &w_minus, const Vector3 &beta,
                                    const Real c_model, Real &gamma_rotation,
                                    bool &used_conjugate_root) {
  Real gamma_minus = 0.0;
  if (GammaFromW(w_minus, c_model, gamma_minus) != StateStatus::success) {
    return PushStatus::invalid_state;
  }
  if (!IsFinite(beta)) {return PushStatus::invalid_input;}

  Real beta_norm = ScaledNorm3(beta);
  if (!IsFinite(beta_norm) || gamma_minus > SafeSquareRootLimit() ||
      beta_norm > SafeSquareRootLimit()) {
    return PushStatus::arithmetic_range;
  }

  Real beta_dot_w = 0.0;
  if (!CheckedDot(beta, w_minus, beta_dot_w)) {
    return PushStatus::arithmetic_range;
  }
  Real beta_dot_u = beta_dot_w/c_model;
  if (!IsFinite(beta_dot_u) || Kokkos::fabs(beta_dot_u) > SafeSquareRootLimit()) {
    return PushStatus::arithmetic_range;
  }

  Real beta2 = beta_norm*beta_norm;
  Real sigma = (gamma_minus - beta_norm)*(gamma_minus + beta_norm);
  Real positive_term = beta2 + beta_dot_u*beta_dot_u;
  Real radical = Kokkos::hypot(sigma, 2.0*Kokkos::sqrt(positive_term));
  if (!IsFinite(sigma) || !IsFinite(positive_term) || !IsFinite(radical)) {
    return PushStatus::arithmetic_range;
  }

  Real gamma2 = 0.0;
  used_conjugate_root = (sigma < 0.0);
  if (!used_conjugate_root) {
    gamma2 = 0.5*(sigma + radical);
  } else {
    Real denominator = radical - sigma;
    if (!IsFinite(denominator) || denominator <= 0.0) {
      return PushStatus::invalid_rotation;
    }
    gamma2 = 2.0*positive_term/denominator;
  }
  if (!IsFinite(gamma2) || gamma2 <= 0.0) {
    return PushStatus::invalid_rotation;
  }
  gamma_rotation = Kokkos::sqrt(gamma2);
  if (!IsFinite(gamma_rotation) || gamma_rotation <= 0.0) {
    return PushStatus::invalid_rotation;
  }
  return PushStatus::success;
}

KOKKOS_INLINE_FUNCTION
PushStatus HigueraCaryPush(const Vector3 &w, const Vector3 &cE, const Vector3 &b,
                          const Real alpha_s, const Real c_model, const Real dt,
                          PushResult &result) {
  if (!IsFinite(w) || !IsFinite(cE) || !IsFinite(b) || !IsFinite(alpha_s) ||
      alpha_s == 0.0 || !IsFinite(dt)) {
    return PushStatus::invalid_input;
  }

  Real gamma_old = 0.0;
  Vector3 velocity_old{0.0, 0.0, 0.0};
  if (VelocityFromW(w, c_model, velocity_old, gamma_old) != StateStatus::success) {
    return PushStatus::invalid_state;
  }

  Real half_step_scale = 0.0;
  if (!CheckedProduct(0.5*dt, alpha_s, half_step_scale)) {
    return PushStatus::arithmetic_range;
  }
  Vector3 epsilon{0.0, 0.0, 0.0};
  Vector3 beta{0.0, 0.0, 0.0};
  Vector3 w_minus{0.0, 0.0, 0.0};
  if (!CheckedScale(half_step_scale, cE, epsilon) ||
      !CheckedScale(half_step_scale, b, beta) ||
      !CheckedAdd(w, epsilon, w_minus)) {
    return PushStatus::arithmetic_range;
  }

  Real gamma_rotation = 0.0;
  bool used_conjugate_root = false;
  PushStatus status = HigueraCaryRotationGamma(
      w_minus, beta, c_model, gamma_rotation, used_conjugate_root);
  if (status != PushStatus::success) {return status;}

  Vector3 t = Scale(1.0/gamma_rotation, beta);
  Real tnorm = ScaledNorm3(t);
  if (!IsFinite(t) || !IsFinite(tnorm) || tnorm > SafeSquareRootLimit()) {
    return PushStatus::arithmetic_range;
  }
  Real t2 = tnorm*tnorm;
  Vector3 s = Scale(2.0/(1.0 + t2), t);
  Vector3 w_cross_t{0.0, 0.0, 0.0};
  Vector3 w_prime{0.0, 0.0, 0.0};
  Vector3 w_prime_cross_s{0.0, 0.0, 0.0};
  Vector3 w_plus{0.0, 0.0, 0.0};
  Vector3 w_new{0.0, 0.0, 0.0};
  if (!IsFinite(s) || !CheckedCross(w_minus, t, w_cross_t) ||
      !CheckedAdd(w_minus, w_cross_t, w_prime) ||
      !CheckedCross(w_prime, s, w_prime_cross_s) ||
      !CheckedAdd(w_minus, w_prime_cross_s, w_plus) ||
      !CheckedAdd(w_plus, epsilon, w_new)) {
    return PushStatus::arithmetic_range;
  }

  Real gamma_new = 0.0;
  Vector3 velocity_new{0.0, 0.0, 0.0};
  if (VelocityFromW(w_new, c_model, velocity_new, gamma_new) !=
      StateStatus::success) {
    return PushStatus::invalid_result;
  }

  Real denominator = gamma_old + gamma_new;
  Vector3 w_sum{0.0, 0.0, 0.0};
  Real cE_dot_w_sum = 0.0;
  Real work_scale = 0.0;
  Real work_numerator = 0.0;
  if (!IsFinite(denominator) || denominator <= 0.0 ||
      !CheckedAdd(w, w_new, w_sum) || !CheckedDot(cE, w_sum, cE_dot_w_sum) ||
      !CheckedProduct(alpha_s, dt, work_scale) ||
      !CheckedProduct(work_scale, cE_dot_w_sum, work_numerator)) {
    return PushStatus::invalid_result;
  }
  Real work = work_numerator/denominator;
  if (!IsFinite(work)) {return PushStatus::invalid_result;}

  result = {w_new, velocity_new, gamma_new, work, used_conjugate_root};
  return PushStatus::success;
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
StateStatus StoreAuthoritativeWAndSynchronizeVelocityShadow(
    const RealView &pr, const int p, const Vector3 &w, const Real c_model) {
  Vector3 velocity{0.0, 0.0, 0.0};
  Real gamma = 0.0;
  StateStatus status = VelocityFromW(w, c_model, velocity, gamma);
  if (status == StateStatus::success) {
    pr(IPWX,p) = w.x;
    pr(IPWY,p) = w.y;
    pr(IPWZ,p) = w.z;
    pr(IPVX,p) = velocity.x;
    pr(IPVY,p) = velocity.y;
    pr(IPVZ,p) = velocity.z;
  }
  return status;
}

} // namespace relativistic
} // namespace particles

#endif // PARTICLES_RELATIVISTIC_PUSHER_HPP_
