//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_fast_speed_test.cpp
//! \brief Unit checks for the CGL fast-magnetosonic characteristic speed.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr Real kTol = 2.0e-13;

void Fail(const std::string &label, const Real got, const Real expected) {
  std::cout << "CGL fast-speed test failed for " << label
            << ": got=" << got << ", expected=" << expected << std::endl;
  std::exit(EXIT_FAILURE);
}

void RequireClose(const std::string &label, const Real got, const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(got) || std::abs(got - expected) > kTol*scale) {
    Fail(label, got, expected);
  }
}

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL fast-speed test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

Real LiteratureDiscriminant(const Real ppar, const Real pperp, const Real bx,
                            const Real by, const Real bz) {
  const Real bx2 = bx*bx;
  const Real b2 = bx2 + by*by + bz*bz;
  const Real mu2 = bx2/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return qsq*qsq + 4.0*pperp*pperp*(1.0 - mu2)*mu2
       - 12.0*ppar*pperp*mu2*(2.0 - mu2)
       + 12.0*ppar*ppar*mu2*mu2 - 12.0*bx2*ppar;
}

Real LiteratureFastSpeed(const Real density, const Real ppar, const Real pperp,
                         const Real bx, const Real by, const Real bz) {
  const Real b2 = bx*bx + by*by + bz*bz;
  const Real mu2 = bx*bx/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return std::sqrt(0.5*(qsq + std::sqrt(LiteratureDiscriminant(
      ppar, pperp, bx, by, bz)))/density);
}

Real LegacyDiscriminant(const Real ppar, const Real pperp, const Real bx,
                        const Real by, const Real bz) {
  const Real bx2 = bx*bx;
  const Real b2 = bx2 + by*by + bz*bz;
  const Real mu2 = bx2/b2;
  const Real qsq = b2 + 2.0*pperp + (2.0*ppar - pperp)*mu2;
  return qsq*qsq + 4.0*pperp*pperp*(1.0 - mu2)*mu2
       - 12.0*ppar*pperp*mu2*(2.0 - mu2)
       + 12.0*ppar*pperp*mu2*mu2 - 12.0*bx2*ppar;
}

void CheckObliqueRegression(const EOS_Data &eos) {
  // This state is inside the CGL hyperbolic bounds, but the legacy
  // p_parallel*p_perp coefficient makes its fast-mode discriminant negative.
  const Real density = 1.0;
  const Real ppar = 1.0;
  const Real pperp = 0.5;
  const Real bx = std::sqrt(0.65);
  const Real by = std::sqrt(0.35);
  const Real bz = 0.0;
  const Real corrected = LiteratureDiscriminant(ppar, pperp, bx, by, bz);
  const Real legacy = LegacyDiscriminant(ppar, pperp, bx, by, bz);
  Require("oblique corrected discriminant is positive", corrected > 1.0);
  Require("oblique legacy discriminant is negative", legacy < -1.0);
  RequireClose("oblique corrected discriminant", corrected, 1.083125);
  RequireClose(
      "oblique fast speed",
      eos.IdealMHDFastSpeed(density, ppar, pperp, bx, by, bz, eos.bfloor),
      1.4169920456419176);
}

void CheckDirectionalLimits(const EOS_Data &eos) {
  const Real perpendicular = eos.IdealMHDFastSpeed(
      2.0, 1.2, 0.7, 0.0, 1.0, 0.0, eos.bfloor);
  RequireClose("perpendicular fast speed", perpendicular, std::sqrt(1.2));

  const Real density = 1.3;
  const Real ppar = 0.9;
  const Real pperp = 0.6;
  const Real bx = 0.8;
  const Real parallel = eos.IdealMHDFastSpeed(
      density, ppar, pperp, bx, 0.0, 0.0, eos.bfloor);
  RequireClose(
      "parallel fast speed", parallel,
      LiteratureFastSpeed(density, ppar, pperp, bx, 0.0, 0.0));
}

} // namespace

void RunCglFastSpeedChecks() {
  EOS_Data eos{};
  eos.gamma = 5.0/3.0;
  eos.bfloor = 1.0e-10;
  CheckObliqueRegression(eos);
  CheckDirectionalLimits(eos);
  std::cout << "CGL fast-speed checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglFastSpeedChecks();
}
