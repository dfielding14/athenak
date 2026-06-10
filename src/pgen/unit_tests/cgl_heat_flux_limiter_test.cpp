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

} // namespace

void RunCglHeatFluxLimiterChecks() {
  CheckModerateValues();
  CheckFiniteOverflowCases();
  CheckInfiniteAsymptotes();
  std::cout << "CGL heat-flux limiter checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglHeatFluxLimiterChecks();
}
