//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_c2p_pressure_floor_test.cpp
//! \brief Unit checks for generic CGL C2P pressure floors and corrected total energy.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {

constexpr bool kSinglePrecision = sizeof(Real) == sizeof(float);
constexpr Real kTol = kSinglePrecision ? 1.0e-5 : 1.0e-13;
constexpr Real kDensity = 1.7;
constexpr Real kVx = 0.11;
constexpr Real kVy = -0.07;
constexpr Real kVz = 0.19;
constexpr Real kBx = 0.91;
constexpr Real kBy = -0.28;
constexpr Real kBz = 0.37;
constexpr Real kPressureFloor = 1.0;

void RequireClose(const std::string &label, const Real got, const Real expected) {
  const Real error = std::abs(got - expected);
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(got) || !std::isfinite(expected) || error > kTol*scale) {
    std::cout << "CGL C2P pressure-floor test failed for " << label
              << ": got " << got << ", expected " << expected
              << ", relative error " << error/scale << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireRelativeClose(const std::string &label, const Real got,
                          const Real expected, const Real tolerance) {
  const Real error = std::abs(got - expected);
  const Real scale = std::abs(expected);
  if (!std::isfinite(got) || !std::isfinite(expected) || scale == 0.0 ||
      error > tolerance*scale) {
    std::cout << "CGL C2P pressure-floor test failed for " << label
              << ": got " << got << ", expected " << expected
              << ", relative error " << error/scale << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL C2P pressure-floor test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

EOS_Data MakeCglEOS() {
  EOS_Data eos{};
  eos.is_ideal = true;
  eos.is_cgl = true;
  eos.dfloor = 1.0e-12;
  eos.pfloor = kPressureFloor;
  eos.bfloor = 1.0e-12;
  return eos;
}

Real KineticEnergy() {
  return 0.5*kDensity*(SQR(kVx) + SQR(kVy) + SQR(kVz));
}

Real MagneticEnergy() {
  return 0.5*(SQR(kBx) + SQR(kBy) + SQR(kBz));
}

MHDCons1D MakeConserved(const Real p_parallel, const Real p_perp) {
  MHDCons1D u{};
  const Real bmag = std::sqrt(SQR(kBx) + SQR(kBy) + SQR(kBz));
  u.d = kDensity;
  u.mx = kDensity*kVx;
  u.my = kDensity*kVy;
  u.mz = kDensity*kVz;
  u.e = 0.5*p_parallel + p_perp + KineticEnergy() + MagneticEnergy();
  u.mu = CGLConservedAnisotropy(kDensity, p_parallel, p_perp, bmag);
  u.bx = kBx;
  u.by = kBy;
  u.bz = kBz;
  return u;
}

void CheckFloorCase(const std::string &label, const Real initial_parallel,
                    const Real initial_perp, const Real expected_parallel,
                    const Real expected_perp) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(initial_parallel, initial_perp);
  HydPrim1D w{};
  bool dfloor_used = false;
  bool efloor_used = false;
  bool tfloor_used = false;
  bool bfloor_used = false;

  SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);

  Require(label + " pressure floor flag", efloor_used);
  Require(label + " no density floor", !dfloor_used);
  Require(label + " no temperature floor", !tfloor_used);
  Require(label + " no magnetic floor", !bfloor_used);
  const Real relative_tolerance = 1000.0*std::numeric_limits<Real>::epsilon();
  RequireRelativeClose(label + " p_parallel", w.e, expected_parallel,
                       relative_tolerance);
  RequireRelativeClose(label + " p_perp", w.pp, expected_perp,
                       relative_tolerance);

  const Real expected_energy =
      0.5*expected_parallel + expected_perp + KineticEnergy() + MagneticEnergy();
  RequireClose(label + " expected total energy", u.e, expected_energy);
  RequireClose(label + " primitive/conserved energy consistency", u.e,
               0.5*w.e + w.pp + KineticEnergy() + MagneticEnergy());
  const Real bmag = std::sqrt(SQR(kBx) + SQR(kBy) + SQR(kBz));
  RequireClose(label + " conserved anisotropy", u.mu,
               CGLConservedAnisotropy(kDensity, expected_parallel, expected_perp, bmag));

  const MHDCons1D corrected = u;
  HydPrim1D recovered{};
  dfloor_used = false;
  efloor_used = false;
  tfloor_used = false;
  bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, recovered, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  Require(label + " repeated conversion floor is precision-compatible",
          !efloor_used || kSinglePrecision);
  RequireRelativeClose(label + " idempotent p_parallel", recovered.e,
                       expected_parallel, relative_tolerance);
  RequireRelativeClose(label + " idempotent p_perp", recovered.pp, expected_perp,
                       relative_tolerance);
  RequireClose(label + " idempotent total energy", u.e, corrected.e);
  RequireClose(label + " idempotent conserved anisotropy", u.mu, corrected.mu);
}

void CheckWideDynamicRangeRoundTrip() {
  const Real exponent =
      0.15*static_cast<Real>(std::numeric_limits<Real>::max_exponent10);
  const Real large = std::pow(static_cast<Real>(10.0), exponent);
  const Real small = 1.0/large;
  const Real density = large;
  const Real p_parallel = small;
  const Real p_perp = large;
  const Real bmag = small;
  const Real eint = 0.5*p_parallel + p_perp;
  const Real tolerance = 1000.0*std::numeric_limits<Real>::epsilon();

  const Real anisotropy =
      CGLConservedAnisotropy(density, p_parallel, p_perp, bmag);
  Require("wide-range conserved anisotropy is finite", std::isfinite(anisotropy));

  Real recovered_parallel = 0.0;
  Real recovered_perp = 0.0;
  CGLRecoverPressuresFromInternalEnergyAndAnisotropy(
      density, eint, anisotropy, bmag, recovered_parallel, recovered_perp);
  RequireRelativeClose("wide-range p_parallel", recovered_parallel, p_parallel,
                       tolerance);
  RequireRelativeClose("wide-range p_perp", recovered_perp, p_perp, tolerance);
  RequireRelativeClose("wide-range internal energy",
                       0.5*recovered_parallel + recovered_perp, eint, tolerance);
}

void CheckIsotropicExactRoundTrip() {
  constexpr Real density = 1.7;
  constexpr Real pressure = 2.0;
  constexpr Real bmag = 0.37;
  const Real anisotropy =
      CGLConservedAnisotropy(density, pressure, pressure, bmag);
  Real recovered_parallel = 0.0;
  Real recovered_perp = 0.0;
  CGLRecoverPressuresFromInternalEnergyAndAnisotropy(
      density, 1.5*pressure, anisotropy, bmag, recovered_parallel, recovered_perp);
  Require("isotropic p_parallel is exact", recovered_parallel == pressure);
  Require("isotropic p_perp is exact", recovered_perp == pressure);
}

void CheckExtremeLogRatioFloor(const std::string &label, const Real log_p_ratio,
                               const Real bmag, const Real expected_parallel,
                               const Real expected_perp) {
  EOS_Data eos = MakeCglEOS();
  eos.pfloor = 1.0e-12;
  eos.bfloor = 1.0e-12;

  MHDCons1D u{};
  u.d = 1.0;
  u.e = 3.0 + 0.5*SQR(bmag);
  u.mu = log_p_ratio - 3.0*std::log(bmag);
  u.bx = bmag;

  HydPrim1D w{};
  bool dfloor_used = false;
  bool efloor_used = false;
  bool tfloor_used = false;
  bool bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);

  Require(label + " pressure floor flag", efloor_used);
  Require(label + " no density floor", !dfloor_used);
  Require(label + " no temperature floor", !tfloor_used);
  Require(label + " no magnetic floor", !bfloor_used);
  RequireClose(label + " p_parallel", w.e, expected_parallel);
  RequireClose(label + " p_perp", w.pp, expected_perp);
  Require(label + " total energy is finite", std::isfinite(u.e));
  Require(label + " conserved anisotropy is finite", std::isfinite(u.mu));
  RequireClose(label + " conserved anisotropy consistency", u.mu,
               CGLConservedAnisotropy(u.d, w.e, w.pp, bmag));
  RequireClose(label + " primitive/conserved energy consistency", u.e,
               0.5*w.e + w.pp + 0.5*SQR(bmag));

  const MHDCons1D corrected = u;
  HydPrim1D recovered{};
  dfloor_used = false;
  efloor_used = false;
  tfloor_used = false;
  bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, recovered, dfloor_used, efloor_used, tfloor_used,
                   bfloor_used);
  Require(label + " repeated conversion floor is precision-compatible",
          !efloor_used || kSinglePrecision);
  RequireClose(label + " idempotent p_parallel", recovered.e, expected_parallel);
  RequireClose(label + " idempotent p_perp", recovered.pp, expected_perp);
  RequireClose(label + " idempotent total energy", u.e, corrected.e);
  RequireClose(label + " idempotent conserved anisotropy", u.mu, corrected.mu);
}

void CheckNegativeInternalEnergyFloor() {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(2.0, 2.0);
  u.e = KineticEnergy() + MagneticEnergy() - 1.0;

  HydPrim1D w{};
  bool dfloor_used = false;
  bool efloor_used = false;
  bool tfloor_used = false;
  bool bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);

  Require("negative-eint pressure floor flag", efloor_used);
  RequireClose("negative-eint p_parallel", w.e, kPressureFloor);
  RequireClose("negative-eint p_perp", w.pp, kPressureFloor);
  Require("negative-eint total energy is finite", std::isfinite(u.e));
  Require("negative-eint conserved anisotropy is finite", std::isfinite(u.mu));
  RequireClose("negative-eint corrected total energy", u.e,
               1.5*kPressureFloor + KineticEnergy() + MagneticEnergy());
}

} // namespace

void RunCglC2PPressureFloorChecks() {
  CheckFloorCase("parallel-only", 0.25, 2.0, kPressureFloor, 2.0);
  CheckFloorCase("perpendicular-only", 2.0, 0.25, 2.0, kPressureFloor);
  CheckFloorCase("both-pressure", 0.25, 0.5, kPressureFloor, kPressureFloor);
  CheckWideDynamicRangeRoundTrip();
  CheckIsotropicExactRoundTrip();
  CheckExtremeLogRatioFloor("large-positive-log-ratio", 1000.0, 1.0e-11,
                            1.0e-12, 3.0);
  CheckExtremeLogRatioFloor("large-negative-log-ratio", -1000.0, 1.0e-11,
                            6.0, 1.0e-12);
  CheckNegativeInternalEnergyFloor();

  std::cout << "CGL C2P pressure-floor checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglC2PPressureFloorChecks();
}
