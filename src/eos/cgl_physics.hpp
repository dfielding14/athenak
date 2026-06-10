#ifndef EOS_CGL_PHYSICS_HPP_
#define EOS_CGL_PHYSICS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_physics.hpp
//! \brief Shared CGL limiter and Landau-fluid closure predicates and constants.

#include <limits>

#include "athena.hpp"

namespace cgl {

constexpr Real kFirehoseObliqueThreshold = -0.7;
constexpr Real kFirehoseParallelThreshold = -1.0;
// Emergency numerical overshoot bound, distinct from either activation policy.
constexpr Real kFirehoseHardBound = -1.5;
constexpr Real kMirrorThreshold = 0.5;
constexpr Real kMirrorHardBound = 1.0;
constexpr Real kBackupCollisionRate = 1.0e10;
constexpr Real kSqrtTwoOverPi = 0.7978845608028654;
constexpr Real kSqrtEightOverPi = 1.5957691216057308;
constexpr Real kSqrtTwoPi = 2.5066282746310002;
constexpr Real kSqrtEightPi = 5.013256549262000;
constexpr Real kThreePiMinusEight = 1.4247779607693793;

KOKKOS_INLINE_FUNCTION
bool FirehoseLimiterActive(const Real paniso, const Real bsqr,
                           const Real firehose_threshold) {
  return paniso <= firehose_threshold*bsqr;
}

KOKKOS_INLINE_FUNCTION
bool FirehoseHardBoundViolated(const Real paniso, const Real bsqr) {
  return paniso <= kFirehoseHardBound*bsqr;
}

KOKKOS_INLINE_FUNCTION
bool MirrorLimiterActive(const Real paniso, const Real bsqr) {
  return paniso >= kMirrorThreshold*bsqr;
}

KOKKOS_INLINE_FUNCTION
bool MirrorHardBoundViolated(const Real paniso, const Real bsqr) {
  return paniso >= kMirrorHardBound*bsqr;
}

KOKKOS_INLINE_FUNCTION
bool ApplyHardwallLimiter(Real &ppar, Real &pperp, const Real bsqr,
                          const bool mirror, const bool firehose,
                          const Real firehose_threshold) {
  const Real paniso = pperp - ppar;
  Real limited_anisotropy = paniso;
  if (firehose && paniso < firehose_threshold*bsqr) {
    limited_anisotropy = firehose_threshold*bsqr;
  } else if (mirror && paniso > kMirrorThreshold*bsqr) {
    limited_anisotropy = kMirrorThreshold*bsqr;
  }
  if (limited_anisotropy == paniso) {
    return false;
  }

  // Scattering preserves internal energy while pinning Delta p to the bound.
  const Real piso = ONE_3RD*ppar + TWO_3RDS*pperp;
  pperp = piso + ONE_3RD*limited_anisotropy;
  ppar = piso - TWO_3RDS*limited_anisotropy;
  return true;
}

KOKKOS_INLINE_FUNCTION
Real LimiterCollisionRate(const Real ppar, const Real pperp, const Real bsqr,
                          const Real limiter_rate, const bool mirror,
                          const bool firehose, const Real firehose_threshold,
                          const bool backup) {
  const Real paniso = pperp - ppar;
  const Real rate = fmax(limiter_rate, static_cast<Real>(0.0));
  Real nu = 0.0;
  if (firehose) {
    if (backup && FirehoseHardBoundViolated(paniso, bsqr)) {
      nu = kBackupCollisionRate;
    } else if (FirehoseLimiterActive(paniso, bsqr, firehose_threshold)) {
      nu = rate;
    }
  }
  if (mirror) {
    if (backup && MirrorHardBoundViolated(paniso, bsqr)) {
      nu = fmax(nu, kBackupCollisionRate);
    } else if (MirrorLimiterActive(paniso, bsqr)) {
      nu = fmax(nu, rate);
    }
  }
  return nu;
}

KOKKOS_INLINE_FUNCTION
Real LimitedHeatFlux(const Real q, const Real qmax) {
  if (!(qmax > 0.0)) return 0.0;
  const Real qabs = fabs(q);
  // Keep each ratio at most unity so neither a product nor a sum can overflow.
  if (qabs <= qmax) return q/(1.0 + qabs/qmax);
  const Real magnitude = qmax/(1.0 + qmax/qabs);
  return (q < 0.0) ? -magnitude : magnitude;
}

KOKKOS_INLINE_FUNCTION
Real LimitedHeatFluxFromRatio(const Real ratio, const Real qmax) {
  if (!(qmax > 0.0)) return 0.0;
  const Real ratio_abs = fabs(ratio);
  if (ratio_abs <= 1.0) {
    return qmax*(ratio/(1.0 + ratio_abs));
  }
  const Real magnitude = qmax/(1.0 + 1.0/ratio_abs);
  return (ratio < 0.0) ? -magnitude : magnitude;
}

KOKKOS_INLINE_FUNCTION
Real PositiveProduct3(const Real a, const Real b, const Real c) {
  const Real inner = b*c;
  const Real direct = a*inner;
  const bool inner_is_usable =
      inner == 0.0 || fabs(inner) >= std::numeric_limits<Real>::min();
  if (Kokkos::isfinite(direct) &&
      inner_is_usable &&
      (direct > 0.0 || a == 0.0 || b == 0.0 || c == 0.0)) {
    return direct;
  }
  const Real alternate1 = (a*b)*c;
  if (Kokkos::isfinite(alternate1) && alternate1 > 0.0) {
    return alternate1;
  }
  const Real alternate2 = (a*c)*b;
  if (Kokkos::isfinite(alternate2) && alternate2 > 0.0) {
    return alternate2;
  }
  if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0) ||
      !Kokkos::isfinite(a) || !Kokkos::isfinite(b) ||
      !Kokkos::isfinite(c)) {
    return direct;
  }

  const Real log_product =
      Kokkos::log2(a) + Kokkos::log2(b) + Kokkos::log2(c);
  const Real max_log = Kokkos::log2(std::numeric_limits<Real>::max());
  if (log_product > max_log) {
    return std::numeric_limits<Real>::infinity();
  }
  const Real min_log =
      Kokkos::log2(std::numeric_limits<Real>::denorm_min());
  if (log_product < min_log) {
    return 0.0;
  }
  return Kokkos::exp2(log_product);
}

KOKKOS_INLINE_FUNCTION
Real ParallelHeatFluxRatio(const Real cparallel, const Real rho,
                           const Real ppar, const Real lf_k,
                           const Real nu, const Real grad_tpar) {
  if (nu == std::numeric_limits<Real>::infinity()) {
    return -0.0*grad_tpar;
  }
  const Real collision_over_speed =
      (kThreePiMinusEight*nu/kSqrtEightPi)/cparallel;
  const Real response = 1.0/(lf_k + collision_over_speed);
  const Real temperature_term = rho*grad_tpar/ppar;
  const Real ratio = -response*temperature_term;
  if (Kokkos::isfinite(ratio) && ratio != 0.0) {
    return ratio;
  }
  const Real reordered_temperature =
      (rho <= ppar) ? (rho/ppar)*grad_tpar
                    : rho*(grad_tpar/ppar);
  const Real reordered_ratio = -response*reordered_temperature;
  if (Kokkos::isfinite(reordered_ratio) && reordered_ratio != 0.0) {
    return reordered_ratio;
  }

  if (!Kokkos::isfinite(cparallel) || !Kokkos::isfinite(rho) ||
      !Kokkos::isfinite(ppar) || !Kokkos::isfinite(lf_k) ||
      !Kokkos::isfinite(nu) || !Kokkos::isfinite(grad_tpar) ||
      !(cparallel > 0.0) || !(rho > 0.0) || !(ppar > 0.0) ||
      !(lf_k > 0.0) || !(nu >= 0.0)) {
    return ratio;
  }
  if (grad_tpar == 0.0) {
    return -0.0*grad_tpar;
  }

  const Real negative_infinity = -std::numeric_limits<Real>::infinity();
  const Real log_k = Kokkos::log2(lf_k);
  const Real log_collision_over_speed =
      (nu > 0.0)
          ? Kokkos::log2(kThreePiMinusEight) + Kokkos::log2(nu) -
                Kokkos::log2(kSqrtEightPi) - Kokkos::log2(cparallel)
          : negative_infinity;
  const Real log_denom_scale = fmax(log_k, log_collision_over_speed);
  const Real scaled_denom =
      Kokkos::exp2(log_k - log_denom_scale) +
      Kokkos::exp2(log_collision_over_speed - log_denom_scale);
  const Real log_response =
      -log_denom_scale - Kokkos::log2(scaled_denom);
  const Real log_temperature =
      Kokkos::log2(rho) + Kokkos::log2(fabs(grad_tpar)) -
      Kokkos::log2(ppar);
  const Real log_ratio = log_response + log_temperature;
  const Real max_log = Kokkos::log2(std::numeric_limits<Real>::max());
  if (log_ratio > max_log) {
    const Real infinity = std::numeric_limits<Real>::infinity();
    return (grad_tpar < 0.0) ? infinity : -infinity;
  }
  const Real magnitude = Kokkos::exp2(log_ratio);
  return (grad_tpar < 0.0) ? magnitude : -magnitude;
}

KOKKOS_INLINE_FUNCTION
Real LimitedParallelHeatFlux(const Real cparallel, const Real rho,
                             const Real ppar, const Real lf_k,
                             const Real nu, const Real grad_tpar,
                             Real &signed_ratio) {
  signed_ratio = ParallelHeatFluxRatio(
      cparallel, rho, ppar, lf_k, nu, grad_tpar);
  const Real qmax = PositiveProduct3(cparallel, kSqrtEightOverPi, ppar);
  return LimitedHeatFluxFromRatio(signed_ratio, qmax);
}

KOKKOS_INLINE_FUNCTION
Real PerpendicularHeatFluxRatio(const Real cparallel, const Real rho,
                                const Real ppar, const Real pperp,
                                const Real bmag_inv, const Real lf_k,
                                const Real nu, const Real grad_tperp,
                                const Real grad_b) {
  if (nu == std::numeric_limits<Real>::infinity()) {
    return -0.0*grad_tperp;
  }
  const Real collision_over_speed =
      (nu/kSqrtTwoPi)/cparallel;
  const Real response = 1.0/(lf_k + collision_over_speed);
  const Real temperature_term = rho*grad_tperp/pperp;
  const Real magnetic_term =
      (1.0 - pperp/ppar)*grad_b*bmag_inv;
  const Real ratio = -response*(temperature_term - magnetic_term);
  if (Kokkos::isfinite(ratio) && ratio != 0.0) {
    return ratio;
  }
  const Real reordered_temperature =
      (rho <= pperp) ? (rho/pperp)*grad_tperp
                     : rho*(grad_tperp/pperp);
  const Real reordered_magnetic =
      ((ppar - pperp)/ppar)*(grad_b*bmag_inv);
  const Real reordered_ratio =
      -response*(reordered_temperature - reordered_magnetic);
  if (Kokkos::isfinite(reordered_ratio) && reordered_ratio != 0.0) {
    return reordered_ratio;
  }

  if (!Kokkos::isfinite(cparallel) || !Kokkos::isfinite(rho) ||
      !Kokkos::isfinite(ppar) || !Kokkos::isfinite(pperp) ||
      !Kokkos::isfinite(bmag_inv) || !Kokkos::isfinite(lf_k) ||
      !Kokkos::isfinite(nu) || !Kokkos::isfinite(grad_tperp) ||
      !Kokkos::isfinite(grad_b) || !(cparallel > 0.0) ||
      !(rho > 0.0) || !(ppar > 0.0) || !(pperp > 0.0) ||
      !(bmag_inv >= 0.0) || !(lf_k > 0.0) || !(nu >= 0.0)) {
    return ratio;
  }

  const Real negative_infinity = -std::numeric_limits<Real>::infinity();
  const Real log_k = Kokkos::log2(lf_k);
  const Real log_collision_over_speed =
      (nu > 0.0)
          ? Kokkos::log2(nu) - Kokkos::log2(kSqrtTwoPi) -
                Kokkos::log2(cparallel)
          : negative_infinity;
  const Real log_denom_scale = fmax(log_k, log_collision_over_speed);
  const Real scaled_denom =
      Kokkos::exp2(log_k - log_denom_scale) +
      Kokkos::exp2(log_collision_over_speed - log_denom_scale);
  const Real log_response =
      -log_denom_scale - Kokkos::log2(scaled_denom);

  const Real pressure_difference = ppar - pperp;
  const bool temperature_zero = (grad_tperp == 0.0);
  const bool magnetic_zero =
      (pressure_difference == 0.0 || grad_b == 0.0 || bmag_inv == 0.0);
  if (temperature_zero && magnetic_zero) {
    return -0.0*grad_tperp;
  }

  Real log_temperature = negative_infinity;
  Real temperature_sign = 0.0;
  if (!temperature_zero) {
    log_temperature = Kokkos::log2(rho) + Kokkos::log2(fabs(grad_tperp)) -
                      Kokkos::log2(pperp);
    temperature_sign = (grad_tperp < 0.0) ? -1.0 : 1.0;
  }

  Real log_magnetic = negative_infinity;
  Real magnetic_sign = 0.0;
  if (!magnetic_zero) {
    log_magnetic = Kokkos::log2(fabs(pressure_difference)) +
                   Kokkos::log2(fabs(grad_b)) + Kokkos::log2(bmag_inv) -
                   Kokkos::log2(ppar);
    const bool magnetic_product_negative =
        (pressure_difference < 0.0) != (grad_b < 0.0);
    magnetic_sign = magnetic_product_negative ? 1.0 : -1.0;
  }

  const Real log_term_scale = fmax(log_temperature, log_magnetic);
  const Real scaled_sum =
      temperature_sign*Kokkos::exp2(log_temperature - log_term_scale) +
      magnetic_sign*Kokkos::exp2(log_magnetic - log_term_scale);
  if (scaled_sum == 0.0) {
    return 0.0;
  }

  const Real log_ratio = log_response + log_term_scale +
                         Kokkos::log2(fabs(scaled_sum));
  const Real max_log = Kokkos::log2(std::numeric_limits<Real>::max());
  if (log_ratio > max_log) {
    const Real infinity = std::numeric_limits<Real>::infinity();
    return (scaled_sum > 0.0) ? -infinity : infinity;
  }
  const Real magnitude = Kokkos::exp2(log_ratio);
  return (scaled_sum > 0.0) ? -magnitude : magnitude;
}

KOKKOS_INLINE_FUNCTION
Real LimitedPerpendicularHeatFlux(const Real cparallel, const Real rho,
                                  const Real ppar, const Real pperp,
                                  const Real bmag_inv, const Real lf_k,
                                  const Real nu, const Real grad_tperp,
                                  const Real grad_b, Real &signed_ratio) {
  signed_ratio = PerpendicularHeatFluxRatio(
      cparallel, rho, ppar, pperp, bmag_inv, lf_k, nu, grad_tperp, grad_b);
  const Real qmax = PositiveProduct3(cparallel, kSqrtTwoOverPi, pperp);
  return LimitedHeatFluxFromRatio(signed_ratio, qmax);
}

} // namespace cgl

#endif // EOS_CGL_PHYSICS_HPP_
