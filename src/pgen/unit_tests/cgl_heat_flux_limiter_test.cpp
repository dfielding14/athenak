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

void CheckModerateValues() {
  RequireClose("positive moderate value", cgl::LimitedHeatFlux(3.0, 2.0), 1.2);
  RequireClose("negative moderate value", cgl::LimitedHeatFlux(-3.0, 2.0), -1.2);
  Require("positive sign", cgl::LimitedHeatFlux(3.0, 2.0) > 0.0);
  Require("negative sign", cgl::LimitedHeatFlux(-3.0, 2.0) < 0.0);
  Require("zero input", cgl::LimitedHeatFlux(0.0, 2.0) == 0.0);
  Require("zero cap", cgl::LimitedHeatFlux(3.0, 0.0) == 0.0);
  Require("negative cap", cgl::LimitedHeatFlux(3.0, -2.0) == 0.0);
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

} // namespace

void RunCglHeatFluxLimiterChecks() {
  CheckModerateValues();
  CheckFiniteOverflowCases();
  CheckInfiniteAsymptotes();
  CheckRatioLimiter();
  CheckPositiveProductEndpoints();
  CheckParallelClosureOrdinaryAgreement();
  CheckParallelClosureOverflowEndpoints();
  CheckPerpendicularClosureOrdinaryAgreement();
  CheckPerpendicularClosureOverflowEndpoints();
  CheckPerpendicularClosureCancellationSigns();
  std::cout << "CGL heat-flux limiter checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglHeatFluxLimiterChecks();
}
