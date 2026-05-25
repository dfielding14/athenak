#ifndef EOS_CGL_PHYSICS_HPP_
#define EOS_CGL_PHYSICS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_physics.hpp
//! \brief Shared CGL limiter and Landau-fluid closure predicates and constants.

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
  return (qmax > 0.0) ? q*qmax/(qmax + fabs(q)) : 0.0;
}

} // namespace cgl

#endif // EOS_CGL_PHYSICS_HPP_
