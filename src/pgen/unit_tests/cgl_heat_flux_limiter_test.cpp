//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_heat_flux_limiter_test.cpp
//! \brief Unit checks for overflow-safe CGL Landau-fluid heat-flux limiting.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "diffusion/cgl_landau_fluid_arithmetic.hpp"
#include "diffusion/sts_rkl2.hpp"
#include "eos/cgl_physics.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kTolerance =
    static_cast<Real>(64.0)*std::numeric_limits<Real>::epsilon();

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL heat-flux limiter test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireClose(const std::string &label, const Real got, const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(got) || std::abs(got - expected) > kTolerance*scale) {
    std::cout << "CGL heat-flux limiter test failed for " << label
              << ": got=" << got << ", expected=" << expected << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireRelativeClose(const std::string &label, const Real got,
                          const Real expected) {
  if (!std::isfinite(got) || got == 0.0 || expected == 0.0 ||
      std::abs((got - expected)/expected) > kTolerance) {
    std::cout << "CGL heat-flux limiter test failed for " << label
              << ": got=" << got << ", expected=" << expected << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireOverflowClose(const std::string &label, const Real got,
                          const Real expected) {
  constexpr Real overflow_tolerance =
      static_cast<Real>(4096.0)*std::numeric_limits<Real>::epsilon();
  if (!std::isfinite(got) || got == 0.0 || expected == 0.0 ||
      std::abs((got - expected)/expected) > overflow_tolerance) {
    std::cout << "CGL heat-flux limiter test failed for " << label
              << ": got=" << got << ", expected=" << expected << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void CheckModerateValues() {
  RequireClose("positive moderate value", cgl::LimitedHeatFlux(3.0, 2.0), 1.2);
  RequireClose("negative moderate value", cgl::LimitedHeatFlux(-3.0, 2.0), -1.2);
  Require("positive sign", cgl::LimitedHeatFlux(3.0, 2.0) > 0.0);
  Require("negative sign", cgl::LimitedHeatFlux(-3.0, 2.0) < 0.0);
  Require("zero input", cgl::LimitedHeatFlux(0.0, 2.0) == 0.0);
  Require("zero cap", cgl::LimitedHeatFlux(3.0, 0.0) == 0.0);
  Require("negative cap", cgl::LimitedHeatFlux(3.0, -2.0) == 0.0);
}

void CheckBackupLimiterPolicy() {
  Require("configured backup remains enabled",
          cgl::EffectiveBackupLimiter(true, false, false, true));
  Require("strict LF does not imply backup",
          !cgl::EffectiveBackupLimiter(false, true, true, true));
  Require("relaxed LF enables backup",
          cgl::EffectiveBackupLimiter(false, true, true, false));
  Require("relaxed LF without instability limiters does not enable backup",
          !cgl::EffectiveBackupLimiter(false, true, false, false));
  Require("unrelated relaxed configuration does not enable backup",
          !cgl::EffectiveBackupLimiter(false, false, true, false));

  constexpr Real ppar = 1.0;
  constexpr Real bsqr = 1.0;
  constexpr Real limiter_rate = 20.0;
  const bool relaxed_backup =
      cgl::EffectiveBackupLimiter(false, true, true, false);
  RequireClose(
      "relaxed mirror ordinary limiter rate",
      cgl::LimiterCollisionRate(
          ppar, 1.75, bsqr, limiter_rate, true, false,
          cgl::kFirehoseObliqueThreshold, relaxed_backup),
      limiter_rate);
  RequireClose(
      "relaxed mirror hard-bound backup rate",
      cgl::LimiterCollisionRate(
          ppar, 2.25, bsqr, limiter_rate, true, false,
          cgl::kFirehoseObliqueThreshold, relaxed_backup),
      cgl::kBackupCollisionRate);
  RequireClose(
      "relaxed firehose hard-bound backup rate",
      cgl::LimiterCollisionRate(
          3.0, 1.0, bsqr, limiter_rate, false, true,
          cgl::kFirehoseParallelThreshold, relaxed_backup),
      cgl::kBackupCollisionRate);
  RequireClose(
      "strict unconfigured hard bound retains finite limiter rate",
      cgl::LimiterCollisionRate(
          ppar, 2.25, bsqr, limiter_rate, true, false,
          cgl::kFirehoseObliqueThreshold, false),
      limiter_rate);
}

void CheckFiniteOverflowCases() {
  const Real maximum = std::numeric_limits<Real>::max();
  const Real half_maximum = maximum/static_cast<Real>(2.0);
  const Real equal_scales = cgl::LimitedHeatFlux(half_maximum, half_maximum);
  RequireClose("equal near-maximum scales", equal_scales,
               maximum/static_cast<Real>(4.0));

  const Real unequal_scales = cgl::LimitedHeatFlux(maximum, half_maximum);
  Require("unequal near-maximum scales finite", std::isfinite(unequal_scales));
  Require("unequal near-maximum scales positive", unequal_scales > 0.0);
  Require("unequal near-maximum scales capped", unequal_scales < half_maximum);

  RequireRelativeClose("maximum equal scales",
                       cgl::LimitedHeatFlux(maximum, maximum), half_maximum);
  RequireRelativeClose("negative maximum equal scales",
                       cgl::LimitedHeatFlux(-maximum, maximum), -half_maximum);

  const Real scale = maximum/static_cast<Real>(16.0);
  RequireRelativeClose("overflowing old numerator",
                       cgl::LimitedHeatFlux(3.0*scale, 2.0*scale), 1.2*scale);
  RequireRelativeClose("negative overflowing old numerator",
                       cgl::LimitedHeatFlux(-3.0*scale, 2.0*scale), -1.2*scale);

  const Real tiny_cap = std::numeric_limits<Real>::denorm_min();
  const Real capped = cgl::LimitedHeatFlux(maximum, tiny_cap);
  Require("maximum input with tiny cap finite", std::isfinite(capped));
  Require("maximum input with tiny cap", capped == tiny_cap);
  Require("tiny input with maximum cap",
          cgl::LimitedHeatFlux(tiny_cap, maximum) == tiny_cap);
}

void CheckInfiniteAsymptotes() {
  const Real infinity = std::numeric_limits<Real>::infinity();
  RequireClose("positive infinite input", cgl::LimitedHeatFlux(infinity, 7.0), 7.0);
  RequireClose("negative infinite input", cgl::LimitedHeatFlux(-infinity, 7.0), -7.0);
  RequireClose("infinite cap", cgl::LimitedHeatFlux(7.0, infinity), 7.0);
  Require("NaN input remains NaN",
          std::isnan(cgl::LimitedHeatFlux(
              std::numeric_limits<Real>::quiet_NaN(), 7.0)));
}

void CheckRatioLimiter() {
  const Real maximum = std::numeric_limits<Real>::max();
  const Real infinity = std::numeric_limits<Real>::infinity();
  const Real negative_zero = -static_cast<Real>(0.0);

  RequireClose("ratio moderate positive",
               cgl::LimitedHeatFluxFromRatio(1.5, 2.0), 1.2);
  RequireClose("ratio moderate negative",
               cgl::LimitedHeatFluxFromRatio(-1.5, 2.0), -1.2);
  RequireRelativeClose("ratio positive endpoint",
                       cgl::LimitedHeatFluxFromRatio(1.0, maximum),
                       maximum/static_cast<Real>(2.0));
  RequireRelativeClose("ratio negative endpoint",
                       cgl::LimitedHeatFluxFromRatio(-1.0, maximum),
                       -maximum/static_cast<Real>(2.0));
  Require("ratio positive infinity caps",
          cgl::LimitedHeatFluxFromRatio(infinity, 7.0) == 7.0);
  Require("ratio negative infinity caps",
          cgl::LimitedHeatFluxFromRatio(-infinity, 7.0) == -7.0);
  Require("ratio positive maximum caps",
          cgl::LimitedHeatFluxFromRatio(maximum, maximum) == maximum);
  Require("ratio negative maximum caps",
          cgl::LimitedHeatFluxFromRatio(-maximum, maximum) == -maximum);
  Require("ratio positive zero preserved",
          !std::signbit(cgl::LimitedHeatFluxFromRatio(0.0, 7.0)));
  Require("ratio negative zero preserved",
          std::signbit(cgl::LimitedHeatFluxFromRatio(negative_zero, 7.0)));
  Require("ratio NaN propagates",
          std::isnan(cgl::LimitedHeatFluxFromRatio(
              std::numeric_limits<Real>::quiet_NaN(), 7.0)));
}

void CheckPositiveProductEndpoints() {
  const Real maximum = std::numeric_limits<Real>::max();
  const Real tiny = std::numeric_limits<Real>::denorm_min();
  const Real inverse_maximum = static_cast<Real>(1.0)/maximum;

  RequireRelativeClose(
      "positive product overflowing inner pair",
      cgl::PositiveProduct3(inverse_maximum, cgl::kSqrtEightOverPi, maximum),
      cgl::kSqrtEightOverPi);
  RequireRelativeClose(
      "positive product tiny final factor",
      cgl::PositiveProduct3(maximum, cgl::kSqrtTwoOverPi, tiny),
      static_cast<Real>(
          static_cast<long double>(maximum)*
          static_cast<long double>(cgl::kSqrtTwoOverPi)*
          static_cast<long double>(tiny)));
  Require("positive product true overflow",
          std::isinf(cgl::PositiveProduct3(maximum, maximum, maximum)));
  Require("positive product true underflow",
          cgl::PositiveProduct3(tiny, tiny, tiny) == 0.0);
  RequireRelativeClose(
      "four-factor overflow cancellation",
      cgl::PositiveProduct4(static_cast<Real>(4.0)/maximum, maximum,
                            cgl::kSqrtEightOverPi, 1.0),
      static_cast<Real>(4.0)*cgl::kSqrtEightOverPi);
  RequireRelativeClose(
      "four-factor subnormal cancellation",
      cgl::PositiveProduct4(tiny, maximum, cgl::kSqrtTwoOverPi, 1.0),
      static_cast<Real>(
          static_cast<long double>(tiny)*
          static_cast<long double>(maximum)*
          static_cast<long double>(cgl::kSqrtTwoOverPi)));
}

Real ReferencePerpendicularHeatFlux(
    const Real cparallel, const Real rho, const Real ppar,
    const Real pperp, const Real bmag, const Real lf_k,
    const Real nu, const Real grad_tperp, const Real grad_b) {
  const Real denominator = cgl::kSqrtTwoPi*cparallel*lf_k + nu;
  const Real chi_perp = static_cast<Real>(2.0)*cparallel*cparallel/denominator;
  const Real q_unlimited = -chi_perp*(
      rho*grad_tperp -
      pperp*(static_cast<Real>(1.0) - pperp/ppar)*grad_b/bmag);
  const Real qmax = cgl::kSqrtTwoOverPi*cparallel*pperp;
  return cgl::LimitedHeatFlux(q_unlimited, qmax);
}

Real ReferenceParallelHeatFlux(
    const Real cparallel, const Real rho, const Real ppar,
    const Real lf_k, const Real nu, const Real grad_tpar) {
  const Real denominator =
      cgl::kSqrtEightPi*cparallel*lf_k +
      cgl::kThreePiMinusEight*nu;
  const Real chi_parallel =
      static_cast<Real>(8.0)*cparallel*cparallel/denominator;
  const Real q_unlimited = -chi_parallel*rho*grad_tpar;
  const Real qmax = cgl::kSqrtEightOverPi*cparallel*ppar;
  return cgl::LimitedHeatFlux(q_unlimited, qmax);
}

void CheckParallelClosureOrdinaryAgreement() {
  constexpr Real cparallel = 2.0;
  constexpr Real rho = 3.0;
  constexpr Real ppar = 5.0;
  constexpr Real lf_k = 0.7;
  constexpr Real nu = 0.3;
  constexpr Real grad_tpar = 0.4;
  Real ratio = 0.0;
  const Real got = cgl::LimitedParallelHeatFlux(
      cparallel, rho, ppar, lf_k, nu, grad_tpar, ratio);
  const Real expected = ReferenceParallelHeatFlux(
      cparallel, rho, ppar, lf_k, nu, grad_tpar);
  RequireRelativeClose("parallel ordinary agreement", got, expected);
  Require("parallel ordinary ratio finite", std::isfinite(ratio));
  Require("parallel ordinary sign preserved", got*ratio > 0.0);
}

void CheckParallelClosureOverflowEndpoints() {
  const Real maximum = std::numeric_limits<Real>::max();

  {
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        maximum, 1.0, 1.0, 1.0, 0.0, -tiny, ratio);
    const Real expected = static_cast<Real>(
        static_cast<long double>(maximum)*
        static_cast<long double>(cgl::kSqrtEightOverPi)*
        static_cast<long double>(tiny));
    RequireClose("parallel subnormal ratio", ratio, tiny);
    Require("parallel subnormal-ratio flux finite", std::isfinite(q));
    RequireRelativeClose("parallel subnormal-ratio flux", q, expected);
  }

  {
    const Real small_ratio = static_cast<Real>(4.0)/maximum;
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        maximum, 1.0, 1.0, 1.0, 0.0, -small_ratio, ratio);
    const Real expected =
        static_cast<Real>(4.0)*cgl::kSqrtEightOverPi/(1.0 + small_ratio);
    RequireRelativeClose("parallel overflowing cap small ratio", ratio, small_ratio);
    Require("parallel overflowing cap finite flux", std::isfinite(q));
    RequireRelativeClose("parallel overflowing cap representable flux", q, expected);
  }

  {
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        1.0, maximum/static_cast<Real>(4.0), 1.0, 1.0, 0.0,
        16.0, ratio);
    Require("parallel product overflow ratio negative",
            std::isinf(ratio) && ratio < 0.0);
    Require("parallel product overflow limited finite", std::isfinite(q));
    Require("parallel product overflow sign preserved", q < 0.0);
    RequireClose("parallel product overflow cap", q,
                 -cgl::kSqrtEightOverPi);
  }

  {
    const Real cparallel = static_cast<Real>(1.0)/maximum;
    const Real ppar = maximum/static_cast<Real>(4.0);
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        cparallel, 1.0, ppar, 1.0, 0.0, 1.0, ratio);
    Require("parallel coefficient-underflow ratio finite", std::isfinite(ratio));
    Require("parallel coefficient-underflow flux finite", std::isfinite(q));
    Require("parallel coefficient-underflow flux remains nonzero", q < 0.0);
  }

  {
    const Real cparallel = static_cast<Real>(1.0)/maximum;
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        cparallel, maximum, maximum, 1.0, 0.0, -1.0, ratio);
    RequireClose("parallel max pressure tiny speed ratio", ratio, 1.0);
    RequireRelativeClose("parallel max pressure tiny speed cap",
                         q, static_cast<Real>(0.5)*cgl::kSqrtEightOverPi);
  }

  {
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    Real ratio = 0.0;
    const Real q = cgl::LimitedParallelHeatFlux(
        maximum, tiny, tiny, 1.0, 0.0, -1.0, ratio);
    const Real qmax =
        cgl::PositiveProduct3(maximum, cgl::kSqrtEightOverPi, tiny);
    RequireClose("parallel max speed tiny pressure ratio", ratio, 1.0);
    RequireRelativeClose("parallel max speed tiny pressure cap",
                         q, static_cast<Real>(0.5)*qmax);
  }

  Real positive_ratio = 0.0;
  Real negative_ratio = 0.0;
  const Real positive_flux = cgl::LimitedParallelHeatFlux(
      1.0, 1.0, 1.0, 1.0, 0.0, -1.0, positive_ratio);
  const Real negative_flux = cgl::LimitedParallelHeatFlux(
      1.0, 1.0, 1.0, 1.0, 0.0, 1.0, negative_ratio);
  Require("parallel negative gradient positive flux",
          positive_ratio > 0.0 && positive_flux > 0.0);
  Require("parallel positive gradient negative flux",
          negative_ratio < 0.0 && negative_flux < 0.0);

  Real suppressed_ratio = 1.0;
  const Real suppressed_flux = cgl::LimitedParallelHeatFlux(
      1.0, 1.0, 1.0, 1.0, std::numeric_limits<Real>::infinity(),
      1.0, suppressed_ratio);
  Require("parallel infinite collision ratio zero", suppressed_ratio == 0.0);
  Require("parallel infinite collision flux zero", suppressed_flux == 0.0);
}

void CheckPerpendicularClosureOrdinaryAgreement() {
  constexpr Real cparallel = 2.0;
  constexpr Real rho = 3.0;
  constexpr Real ppar = 5.0;
  constexpr Real pperp = 4.0;
  constexpr Real bmag = 2.0;
  constexpr Real lf_k = 0.7;
  constexpr Real nu = 0.3;
  constexpr Real grad_tperp = 0.4;
  constexpr Real grad_b = -0.2;
  Real ratio = 0.0;
  const Real got = cgl::LimitedPerpendicularHeatFlux(
      cparallel, rho, ppar, pperp, static_cast<Real>(1.0)/bmag, lf_k, nu,
      grad_tperp, grad_b, ratio);
  const Real expected = ReferencePerpendicularHeatFlux(
      cparallel, rho, ppar, pperp, bmag, lf_k, nu,
      grad_tperp, grad_b);
  RequireRelativeClose("perpendicular ordinary agreement", got, expected);
  Require("perpendicular ordinary ratio finite", std::isfinite(ratio));
  Require("perpendicular ordinary sign preserved", got*ratio > 0.0);
}

void CheckPerpendicularClosureOverflowEndpoints() {
  const Real maximum = std::numeric_limits<Real>::max();

  {
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        maximum, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
        -tiny, 0.0, ratio);
    const Real expected = static_cast<Real>(
        static_cast<long double>(maximum)*
        static_cast<long double>(cgl::kSqrtTwoOverPi)*
        static_cast<long double>(tiny));
    RequireClose("perpendicular subnormal ratio", ratio, tiny);
    Require("perpendicular subnormal-ratio flux finite", std::isfinite(q));
    RequireRelativeClose("perpendicular subnormal-ratio flux", q, expected);
  }

  {
    const Real small_ratio = static_cast<Real>(4.0)/maximum;
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        maximum, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
        -small_ratio, 0.0, ratio);
    const Real expected =
        static_cast<Real>(4.0)*cgl::kSqrtTwoOverPi/(1.0 + small_ratio);
    RequireRelativeClose(
        "perpendicular overflowing cap small ratio", ratio, small_ratio);
    Require("perpendicular overflowing cap finite flux", std::isfinite(q));
    RequireRelativeClose(
        "perpendicular overflowing cap representable flux", q, expected);
  }

  {
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        1.0, maximum/static_cast<Real>(4.0), 1.0, 1.0, 1.0,
        1.0, 0.0, 16.0, 0.0, ratio);
    Require("temperature-product overflow ratio negative",
            std::isinf(ratio) && ratio < 0.0);
    Require("temperature-product overflow limited finite", std::isfinite(q));
    Require("temperature-product overflow sign preserved", q < 0.0);
    RequireClose("temperature-product overflow cap", q,
                 -cgl::kSqrtTwoOverPi);
  }

  {
    const Real pperp = maximum/static_cast<Real>(8.0);
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        1.0, 1.0, 1.0, pperp, 1.0, 1.0, 0.0,
        0.0, 16.0, ratio);
    const Real qmax = cgl::kSqrtTwoOverPi*pperp;
    Require("anisotropy-product overflow ratio negative",
            std::isinf(ratio) && ratio < 0.0);
    Require("anisotropy-product overflow limited finite", std::isfinite(q));
    Require("anisotropy-product overflow sign preserved", q < 0.0);
    RequireRelativeClose("anisotropy-product overflow cap", q, -qmax);
  }

  {
    const Real cparallel = static_cast<Real>(1.0)/maximum;
    const Real pperp = maximum/static_cast<Real>(4.0);
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        cparallel, 1.0, pperp, pperp, 1.0, 1.0, 0.0,
        1.0, 0.0, ratio);
    Require("coefficient-underflow ratio finite", std::isfinite(ratio));
    Require("coefficient-underflow flux finite", std::isfinite(q));
    Require("coefficient-underflow flux remains nonzero", q < 0.0);
  }

  {
    const Real cparallel = static_cast<Real>(1.0)/maximum;
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        cparallel, maximum, maximum, maximum, 1.0, 1.0, 0.0,
        -1.0, 0.0, ratio);
    RequireClose("perpendicular max pressure tiny speed ratio", ratio, 1.0);
    RequireRelativeClose("perpendicular max pressure tiny speed cap",
                         q, static_cast<Real>(0.5)*cgl::kSqrtTwoOverPi);
  }

  {
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    Real ratio = 0.0;
    const Real q = cgl::LimitedPerpendicularHeatFlux(
        maximum, tiny, tiny, tiny, 1.0, 1.0, 0.0,
        -1.0, 0.0, ratio);
    const Real qmax =
        cgl::PositiveProduct3(maximum, cgl::kSqrtTwoOverPi, tiny);
    RequireClose("perpendicular max speed tiny pressure ratio", ratio, 1.0);
    RequireRelativeClose("perpendicular max speed tiny pressure cap",
                         q, static_cast<Real>(0.5)*qmax);
  }
}

void CheckPerpendicularClosureCancellationSigns() {
  constexpr Real cparallel = 1.0;
  const Real ppar = std::numeric_limits<Real>::max()/static_cast<Real>(8.0);
  const Real pperp = static_cast<Real>(2.0)*ppar;
  const Real rho = pperp;
  constexpr Real bmag = 1.0;
  constexpr Real lf_k = 1.0;
  constexpr Real nu = 0.0;
  constexpr Real grad_tperp = 16.0;
  constexpr Real cancelling_grad_b = -16.0;

  Real ratio = 1.0;
  const Real cancelled = cgl::LimitedPerpendicularHeatFlux(
      cparallel, rho, ppar, pperp, static_cast<Real>(1.0)/bmag, lf_k, nu,
      grad_tperp, cancelling_grad_b, ratio);
  Require("perpendicular cancellation ratio zero", ratio == 0.0);
  Require("perpendicular cancellation flux zero", cancelled == 0.0);

  const Real grad_b_more_negative = std::nextafter(
      cancelling_grad_b, -std::numeric_limits<Real>::infinity());
  const Real grad_b_less_negative = std::nextafter(
      cancelling_grad_b, static_cast<Real>(0.0));
  Real positive_ratio = 0.0;
  Real negative_ratio = 0.0;
  const Real positive_flux = cgl::LimitedPerpendicularHeatFlux(
      cparallel, rho, ppar, pperp, static_cast<Real>(1.0)/bmag, lf_k, nu,
      grad_tperp, grad_b_more_negative, positive_ratio);
  const Real negative_flux = cgl::LimitedPerpendicularHeatFlux(
      cparallel, rho, ppar, pperp, static_cast<Real>(1.0)/bmag, lf_k, nu,
      grad_tperp, grad_b_less_negative, negative_ratio);
  Require("perpendicular cancellation positive perturbation",
          positive_ratio > 0.0 && positive_flux > 0.0);
  Require("perpendicular cancellation negative perturbation",
          negative_ratio < 0.0 && negative_flux < 0.0);

  Real suppressed_ratio = 1.0;
  const Real suppressed_flux = cgl::LimitedPerpendicularHeatFlux(
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      std::numeric_limits<Real>::infinity(), 1.0, 1.0,
      suppressed_ratio);
  Require("perpendicular infinite collision ratio zero", suppressed_ratio == 0.0);
  Require("perpendicular infinite collision flux zero", suppressed_flux == 0.0);

  Real zero_inverse_field_ratio = 1.0;
  const Real zero_inverse_field_flux = cgl::LimitedPerpendicularHeatFlux(
      1.0, 1.0, std::numeric_limits<Real>::denorm_min(),
      std::numeric_limits<Real>::max(), 0.0, 1.0, 0.0,
      0.0, std::numeric_limits<Real>::max(), zero_inverse_field_ratio);
  Require("perpendicular zero inverse field ratio zero",
          zero_inverse_field_ratio == 0.0);
  Require("perpendicular zero inverse field flux zero",
          zero_inverse_field_flux == 0.0);
}

void CheckWeightedFluxArithmetic() {
  constexpr Real cparallel = 2.0;
  constexpr Real ppar = 5.0;
  constexpr Real pperp = 4.0;
  constexpr Real qpar_ratio = 0.4;
  constexpr Real qperp_ratio = -1.7;
  constexpr Real bhdir = 0.6;
  constexpr Real bmag_inv = 0.5;
  constexpr Real dt_sweep = 0.1;
  constexpr Real rkl_weight = 0.03;

  const Real qpar = cgl::LimitedHeatFluxFromRatioAndScale(
      qpar_ratio, cparallel, cgl::kSqrtEightOverPi, ppar);
  const Real qperp = cgl::LimitedHeatFluxFromRatioAndScale(
      qperp_ratio, cparallel, cgl::kSqrtTwoOverPi, pperp);
  const Real stage_weight = dt_sweep*rkl_weight;

  auto weighted_qpar = cgl_lf::LimitedHeatFlux(
      qpar_ratio, cparallel, cgl::kSqrtEightOverPi, ppar);
  weighted_qpar = cgl_lf::Multiply(weighted_qpar, bhdir);
  weighted_qpar = cgl_lf::Multiply(weighted_qpar, dt_sweep);
  weighted_qpar = cgl_lf::Multiply(weighted_qpar, rkl_weight);
  auto weighted_qperp = cgl_lf::LimitedHeatFlux(
      qperp_ratio, cparallel, cgl::kSqrtTwoOverPi, pperp);
  weighted_qperp = cgl_lf::Multiply(weighted_qperp, bhdir);
  weighted_qperp = cgl_lf::Multiply(weighted_qperp, dt_sweep);
  weighted_qperp = cgl_lf::Multiply(weighted_qperp, rkl_weight);

  const Real energy = cgl_lf::Materialize(cgl_lf::Add(
      weighted_qperp,
      cgl_lf::Multiply(weighted_qpar, static_cast<Real>(0.5))));
  const Real moment = cgl_lf::Materialize(
      cgl_lf::Multiply(weighted_qperp, bmag_inv));
  RequireRelativeClose(
      "weighted energy moderate equivalence", energy,
      stage_weight*bhdir*(qperp + static_cast<Real>(0.5)*qpar));
  RequireRelativeClose(
      "weighted moment moderate equivalence", moment,
      stage_weight*bhdir*qperp*bmag_inv);
}

void CheckWeightedFluxOverflowEndpoints() {
  const Real maximum = std::numeric_limits<Real>::max();
  const Real inverse_maximum = static_cast<Real>(1.0)/maximum;

  auto physical_overflow = cgl_lf::LimitedHeatFlux(
      1.0, maximum, cgl::kSqrtEightOverPi, 4.0);
  Require("unweighted heat flux truly overflows",
          std::isinf(cgl_lf::Materialize(physical_overflow)));
  const Real weighted_flux = cgl_lf::Materialize(
      cgl_lf::Multiply(physical_overflow, inverse_maximum));
  RequireOverflowClose(
      "overflowing heat flux becomes representable after stage weight",
      weighted_flux, static_cast<Real>(2.0)*cgl::kSqrtEightOverPi);

  auto weighted_moment = cgl_lf::LimitedHeatFlux(
      1.0, maximum, cgl::kSqrtTwoOverPi, 1.0);
  Require("raw qperp over B truly overflows",
          std::isinf(cgl_lf::Materialize(weighted_moment)*maximum));
  weighted_moment = cgl_lf::Multiply(weighted_moment, inverse_maximum);
  weighted_moment = cgl_lf::Multiply(weighted_moment, inverse_maximum);
  weighted_moment = cgl_lf::Multiply(weighted_moment, maximum);
  RequireOverflowClose(
      "qperp over B avoids overflowing intermediate",
      cgl_lf::Materialize(weighted_moment),
      static_cast<Real>(0.5)*cgl::kSqrtTwoOverPi);

  const Real divf = maximum/static_cast<Real>(4.0);
  Require("raw dt times divF truly overflows",
          std::isinf(static_cast<Real>(8.0)*divf));
  auto weighted_rhs = cgl_lf::Multiply(
      cgl_lf::FromReal(divf), static_cast<Real>(8.0));
  weighted_rhs = cgl_lf::Multiply(weighted_rhs, static_cast<Real>(0.125));
  RequireOverflowClose(
      "weighted dt times divF avoids overflowing intermediate",
      cgl_lf::Materialize(weighted_rhs), divf);

  auto weighted_work = cgl_lf::LimitedHeatFlux(
      1.0, maximum, cgl::kSqrtTwoOverPi, 1.0);
  weighted_work = cgl_lf::Multiply(weighted_work, inverse_maximum);
  weighted_work = cgl_lf::Multiply(weighted_work, 4.0);
  RequireOverflowClose(
      "diagnostic work weights before physical power overflow",
      cgl_lf::Materialize(weighted_work),
      static_cast<Real>(2.0)*cgl::kSqrtTwoOverPi);

  const auto true_overflow =
      cgl_lf::Multiply(cgl_lf::FromReal(maximum), maximum);
  Require("true weighted overflow remains visible",
          std::isinf(cgl_lf::Materialize(true_overflow)));

  auto endpoint_round_trip =
      cgl_lf::Multiply(cgl_lf::FromReal(maximum), 2.0);
  endpoint_round_trip = cgl_lf::Multiply(endpoint_round_trip, 0.5);
  Require("near-maximum scaled round trip remains finite",
          cgl_lf::Materialize(endpoint_round_trip) == maximum);

  const auto positive = cgl_lf::Multiply(
      cgl_lf::FromReal(maximum), static_cast<Real>(4.0));
  const auto negative = cgl_lf::Multiply(
      cgl_lf::FromReal(-maximum), static_cast<Real>(4.0));
  Require("overflowing energy components cancel before materialization",
          cgl_lf::Materialize(cgl_lf::Add(positive, negative)) == 0.0);

  const Real minimum = std::numeric_limits<Real>::min();
  const Real below_minimum = std::nextafter(
      minimum, static_cast<Real>(0.0));
  Require(
      "close scaled cancellation retains subnormal residual",
      cgl_lf::Materialize(cgl_lf::Add(
          cgl_lf::FromReal(minimum), cgl_lf::FromReal(-below_minimum))) ==
          std::numeric_limits<Real>::denorm_min());
}

void CheckWeightedRKLCacheAlgebra() {
  constexpr Real first_rhs = 7.0;
  constexpr Real current_rhs = -3.0;
  for (const int nstages : {3, 5, 101}) {
    const Real first_weight = cgl_lf::FirstStageRKLWeight(nstages);
    const auto first_coeffs =
        parabolic::ComputeRKL2Coefficients(1, nstages);
    RequireRelativeClose(
        "first-stage RKL weight matches controller",
        first_weight, first_coeffs.muj_tilde);
    const Real cached_rhs = first_weight*first_rhs;
    for (int stage = 1; stage <= nstages; ++stage) {
      const auto coeffs =
          parabolic::ComputeRKL2Coefficients(stage, nstages);
      const Real old_terms =
          coeffs.gammaj_tilde*first_rhs +
          coeffs.muj_tilde*current_rhs;
      const Real new_terms =
          cgl_lf::CachedRHSCoefficient(
              coeffs.gammaj_tilde, nstages)*cached_rhs +
          coeffs.muj_tilde*current_rhs;
      RequireClose("weighted cached-RHS equivalence", new_terms, old_terms);
    }
  }

  RequireClose("single-stage first weight",
               cgl_lf::FirstStageRKLWeight(1), 1.0);
  RequireClose("single-stage cached coefficient",
               cgl_lf::CachedRHSCoefficient(0.0, 1), 0.0);
  RequireClose("single-stage weighted RHS", 1.0*current_rhs, current_rhs);
  constexpr Real synthetic_gamma = -0.375;
  RequireClose(
      "single-stage synthetic cached-RHS identity",
      cgl_lf::CachedRHSCoefficient(synthetic_gamma, 1)*first_rhs,
      synthetic_gamma*first_rhs);

  const Real maximum = std::numeric_limits<Real>::max();
  const Real state = static_cast<Real>(0.75)*maximum;
  Require("raw RKL state term overflows before cancellation",
          std::isinf(static_cast<Real>(1.5)*state));
  RequireOverflowClose(
      "scaled RKL state terms cancel to a finite state",
      cgl_lf::WeightedRKL2Update(
          1.5, state, -0.5, state, 0.0, 0.0, 0.0, 0.0, 0.0),
      state);
}

} // namespace

void RunCglHeatFluxLimiterChecks() {
  CheckModerateValues();
  CheckBackupLimiterPolicy();
  CheckFiniteOverflowCases();
  CheckInfiniteAsymptotes();
  CheckRatioLimiter();
  CheckPositiveProductEndpoints();
  CheckParallelClosureOrdinaryAgreement();
  CheckParallelClosureOverflowEndpoints();
  CheckPerpendicularClosureOrdinaryAgreement();
  CheckPerpendicularClosureOverflowEndpoints();
  CheckPerpendicularClosureCancellationSigns();
  CheckWeightedFluxArithmetic();
  CheckWeightedFluxOverflowEndpoints();
  CheckWeightedRKLCacheAlgebra();
  std::cout << "CGL heat-flux limiter checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglHeatFluxLimiterChecks();
}
