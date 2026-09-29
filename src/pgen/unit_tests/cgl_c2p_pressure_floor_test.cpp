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
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
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
  eos.sigma_max = std::numeric_limits<Real>::max();
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

void CheckIdempotence(const std::string &label, MHDCons1D &u,
                      const HydPrim1D &first, const EOS_Data &eos) {
  const MHDCons1D corrected = u;
  HydPrim1D recovered{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, recovered, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  Require(label + " no repeated pressure floor", !efloor_used);
  Require(label + " bitwise conserved state",
          std::memcmp(&u, &corrected, sizeof(u)) == 0);
  Require(label + " bitwise primitive state",
          std::memcmp(&recovered, &first, sizeof(first)) == 0);
  const Real kinetic = 0.5*(1.0/u.d)*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  const Real magnetic = 0.5*(SQR(u.bx) + SQR(u.by) + SQR(u.bz));
  const Real tolerance = kSinglePrecision ? 32.0*std::numeric_limits<Real>::epsilon()
                                        : 1.0e-14;
  RequireRelativeClose(label + " internal energy identity", u.e - kinetic - magnetic,
                       0.5*first.e + first.pp, tolerance);
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

  CheckIdempotence(label, u, w, eos);
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

  CheckIdempotence(label, u, w, eos);
}

void CheckMagneticFloorCase(const Real bmag, const Real pressure) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(2.0, 1.0);
  u.bx = bmag;
  u.by = u.bz = 0.0;
  u.e = 1.5*pressure + KineticEnergy() + 0.5*SQR(bmag);
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  Require("magnetic floor flag", bfloor_used);
  Require("magnetic floor isotropy", w.e == w.pp);
  CheckIdempotence("magnetic floor", u, w, eos);
}

void CheckUnresolvedInternalEnergyFloor() {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u{};
  u.d = 1.0;
  u.mx = 1.0e12;
  u.bx = 1.0e8;
  u.e = Kokkos::nextafter(static_cast<Real>(0.5*SQR(u.mx) + 0.5*SQR(u.bx)),
                          static_cast<Real>(0.0));
  u.mu = CGLConservedAnisotropy(u.d, 1.0, 1.0, u.bx);
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  Require("unresolved internal energy uses floor", efloor_used);
  Require("representable floor pressures", w.e >= eos.pfloor && w.pp >= eos.pfloor);
  CheckIdempotence("unresolved internal energy", u, w, eos);
}

void CheckDensityFloorRatios() {
  EOS_Data eos = MakeCglEOS();
  eos.dfloor = 1.0e-4;
  eos.pfloor = 1.0e-12;
  const Real densities[] = {5.0e-5, 1.0e-5, 1.0e-6, 1.0e-8};
  const Real ratios[] = {1.0, 0.3, 4.0};
  const Real tolerance = kSinglePrecision ? 2.0e-5 : 1.0e-12;
  for (const Real density : densities) {
    for (const Real ratio : ratios) {
      MHDCons1D u{};
      u.d = density;
      u.mx = 0.2*density;
      u.bx = 1.0;
      u.e = 1.0e-3*(0.5 + ratio) + 0.5*SQR(u.mx)/density + 0.5;
      u.mu = CGLConservedAnisotropy(density, 1.0e-3, 1.0e-3*ratio, 1.0);
      const Real original_energy = u.e;
      HydPrim1D w{};
      bool dfloor_used = false, efloor_used = false;
      bool tfloor_used = false, bfloor_used = false;
      SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
      Require("density floor flag", dfloor_used);
      Require("density floor value", u.d == eos.dfloor && w.d == eos.dfloor);
      Require("density floor preserves total energy", u.e == original_energy);
      RequireRelativeClose("density floor preserves pressure ratio", w.pp/w.e,
                           ratio, tolerance);
      CheckIdempotence("density floor", u, w, eos);
    }
  }
}

void CheckInvalidAnisotropy() {
  EOS_Data eos = MakeCglEOS();
  eos.dfloor = 1.0e-4;
  eos.pfloor = 1.0e-12;
  const Real logarithms[] = {std::numeric_limits<Real>::quiet_NaN(),
                            std::numeric_limits<Real>::infinity(),
                            -std::numeric_limits<Real>::infinity(), 1000.0,
                            std::log(std::numeric_limits<Real>::max())
                                - static_cast<Real>(1.0)};
  for (const Real logarithm : logarithms) {
    MHDCons1D u{};
    u.d = 1.0e-20;
    u.bx = 1.0;
    u.e = 3.0e-3 + 0.5;
    u.mu = u.d*logarithm;
    HydPrim1D w{};
    bool dfloor_used = false, efloor_used = false;
    bool tfloor_used = false, bfloor_used = false;
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
    Require("invalid anisotropy recovery flagged", efloor_used);
    Require("invalid anisotropy recovers finite A", std::isfinite(u.mu));
    Require("invalid anisotropy recovers isotropy", w.e == w.pp);
    Require("invalid anisotropy recovers positive pressures", w.e >= eos.pfloor);
    CheckIdempotence("invalid anisotropy", u, w, eos);
  }
}

void CheckInvalidDensity() {
  const EOS_Data eos = MakeCglEOS();
  const Real densities[] = {0.0, -1.0, std::numeric_limits<Real>::quiet_NaN(),
                            std::numeric_limits<Real>::infinity()};
  for (const Real density : densities) {
    MHDCons1D u{};
    u.d = density;
    u.bx = 1.0;
    u.e = 3.5;
    u.mu = 1.0;
    HydPrim1D w{};
    bool dfloor_used = false, efloor_used = false;
    bool tfloor_used = false, bfloor_used = false;
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
    Require("invalid density recovery flagged", dfloor_used && efloor_used);
    Require("invalid density floored", u.d == eos.dfloor);
    Require("invalid density recovers finite A", std::isfinite(u.mu));
    Require("invalid density recovers isotropy", w.e == w.pp);
    CheckIdempotence("invalid density", u, w, eos);
  }
}

void CheckMagnetizationDensityFloor() {
  EOS_Data eos = MakeCglEOS();
  eos.pfloor = 1.0e-12;
  eos.sigma_max = 2.0;
  for (const bool magnetic_moment : {false, true}) {
    MHDCons1D u{};
    u.d = 0.25;
    u.mx = 0.1;
    u.bx = 2.0;
    u.e = 2.5 + 0.5*SQR(u.mx)/u.d + 0.5*SQR(u.bx);
    u.mu = magnetic_moment ? 1.0 : CGLConservedAnisotropy(u.d, 1.0, 2.0, u.bx);
    const Real original_energy = u.e;
    HydPrim1D w{};
    bool dfloor_used = false, efloor_used = false;
    bool tfloor_used = false, bfloor_used = false;
    if (magnetic_moment) {
      SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                         tfloor_used, bfloor_used);
    } else {
      SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
    }
    Require("magnetization density floor flag", dfloor_used);
    Require("magnetization ceiling", u.d == 2.0 && w.d == 2.0);
    Require("magnetization floor preserves total energy", u.e == original_energy);
    if (magnetic_moment) {
      RequireClose("magnetization floor preserves magnetic moment pressure", w.pp, 2.0);
      RequireClose("magnetization floor moment internal energy", u.e,
                   0.5*w.e + w.pp + 0.5*SQR(u.mx)/u.d + 0.5*SQR(u.bx));
    } else {
      RequireRelativeClose("magnetization floor preserves pressure ratio", w.pp/w.e,
                           2.0, kSinglePrecision ? 2.0e-5 : 1.0e-12);
      CheckIdempotence("magnetization floor", u, w, eos);
    }
  }
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

// The map is tested directly so roundoff cannot be hidden by C2P repair.
void CheckCollisionMap() {
  EOS_Data eos = MakeCglEOS();
  eos.firehose_threshold = 1.4;
  eos.mirror_threshold = 1.0;
  eos.firehose_backup_factor = 1.0;
  eos.mirror_backup_factor = 2.0;
  const Real mean_tolerance = kSinglePrecision ? 8.0*std::numeric_limits<Real>::epsilon()
                                              : 1.0e-14;
  for (Real nudt : {static_cast<Real>(0.0), static_cast<Real>(1.0),
                    static_cast<Real>(1.0e10)}) {
    eos.lim_coll = nudt;
    Real previous = -std::numeric_limits<Real>::max();
    for (int n = 0; n <= 80; ++n) {
      const Real delta = -2.0 + 0.05*n;
      MHDPrim1D w{};
      w.bx = 1.0;
      w.e = 10.0 - TWO_3RDS*delta;
      w.pp = 10.0 + ONE_3RD*delta;
      const Real initial = w.pp - w.e;
      const Real piso = ONE_3RD*w.e + TWO_3RDS*w.pp;
      const MHDPrim1D original = w;
      const Real threshold = initial < -0.7 ? -0.7 : 0.5;
      const Real expected = initial < -0.7 || initial > 0.5
          ? threshold + (initial - threshold)/(1.0 + nudt) : initial;
      SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
      RequireClose("limiter backward Euler", w.pp - w.e, expected);
      Require("monotone limiter", w.pp - w.e >= previous);
      previous = w.pp - w.e;
      RequireRelativeClose("rate mean pressure", ONE_3RD*w.e + TWO_3RDS*w.pp,
                           piso, mean_tolerance);
      if (nudt == 0.0) {
        Require("zero-rate identity", std::memcmp(&w, &original, sizeof(w)) == 0);
      }
      if (nudt == 1.0e10 && (initial < -0.7 || initial > 0.5)) {
        Require("stiff soft-threshold limit",
                std::abs(w.pp - w.e - threshold) <
                    (kSinglePrecision ? 4.0e-6 : 2.0e-10));
      }
      SingleCollWalls_CGLMHD(w, eos, true);
      const MHDPrim1D projected = w;
      SingleCollWalls_CGLMHD(w, eos, true);
      Require("bitwise idempotent walls", std::memcmp(&w, &projected, sizeof(w)) == 0);
    }
  }
  eos.lim_coll = 1.0;
  Real adjacent[2];
  for (int n = 0; n < 2; ++n) {
    const Real delta = n == 0 ? 0.99 : 1.01;
    MHDPrim1D w{};
    w.bx = 1.0;
    w.e = 3.0 - TWO_3RDS*delta;
    w.pp = 3.0 + ONE_3RD*delta;
    SingleCollRates_CGLMHD(w, eos, 1.0, true, false);
    adjacent[n] = w.pp - w.e;
  }
  RequireClose("no backup-band hysteresis", adjacent[1] - adjacent[0], 0.01);

  // Backup walls ignore soft-limiter flags; the fluid firehose wall ignores backup.
  for (bool backup : {false, true}) {
    eos.mlim = eos.flim = false;
    for (Real delta : {static_cast<Real>(-2.0), static_cast<Real>(2.0)}) {
      MHDPrim1D w{};
      w.bx = 1.0;
      w.e = 3.0 - TWO_3RDS*delta;
      w.pp = 3.0 + ONE_3RD*delta;
      SingleCollWalls_CGLMHD(w, eos, backup);
      const Real expected = backup ? (delta < 0.0 ? -0.7 : 1.0)
                                   : (delta < 0.0 ? -1.0 : delta);
      RequireClose("independent wall flags", w.pp - w.e, expected);
    }
  }

  std::mt19937_64 generator(81);
  std::uniform_real_distribution<double> exponent(-30.0, 30.0);
  for (int n = 0; n < 10000; ++n) {
    MHDPrim1D w{};
    w.e = std::exp(exponent(generator));
    w.pp = std::exp(exponent(generator));
    w.bx = std::exp(exponent(generator));
    const Real piso = ONE_3RD*w.e + TWO_3RDS*w.pp;
    eos.nu_coll = n%3 == 0 ? 0.0 : 0.4;
    eos.lim_coll = n%3 == 0 ? 0.0 : (n%3 == 1 ? 1.0 : 1.0e10);
    SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
    SingleCollWalls_CGLMHD(w, eos, n%2);
    Require("randomized positive pressures", w.e > 0.0 && w.pp > 0.0);
    Require("randomized wall admissibility",
            !cgl::HardBoundViolated(w.pp - w.e, SQR(w.bx), eos, n%2));
    RequireRelativeClose("randomized conserved mean pressure",
                         ONE_3RD*w.e + TWO_3RDS*w.pp, piso, mean_tolerance);
    const MHDPrim1D projected = w;
    SingleCollWalls_CGLMHD(w, eos, n%2);
    Require("randomized bitwise idempotent walls",
            std::memcmp(&w, &projected, sizeof(w)) == 0);
    w.d = std::exp(exponent(generator));
    const Real eint = 0.5*w.e + w.pp;
    const Real bmag = fmax(w.bx, eos.bfloor);
    const Real candidate = CGLConservedAnisotropy(w.d, w.e, w.pp, bmag);
    const Real admissible = CGLWallAdmissibleAnisotropy(w, eint, candidate, eos, n%2);
    Real recovered_parallel, recovered_perp;
    CGLRecoverPressuresFromInternalEnergyAndAnisotropy(
        w.d, eint, admissible, bmag, recovered_parallel, recovered_perp);
    Require("encoded wall admissibility",
            !cgl::HardBoundViolated(recovered_perp - recovered_parallel,
                                    SQR(w.bx), eos, n%2));
    RequireRelativeClose("encoded wall mean pressure",
                         ONE_3RD*recovered_parallel + TWO_3RDS*recovered_perp,
                         piso, mean_tolerance);
    Require("encoded wall fixed point",
            CGLWallAdmissibleAnisotropy(w, eint, admissible, eos, n%2) == admissible);
  }
}

} // namespace

void RunCglC2PPressureFloorChecks() {
  CheckFloorCase("parallel-only", 0.25, 2.0, kPressureFloor, 2.0);
  CheckFloorCase("perpendicular-only", 2.0, 0.25, 2.0, kPressureFloor);
  CheckFloorCase("both-pressure", 0.25, 0.5, kPressureFloor, kPressureFloor);
  for (Real pressure : {kPressureFloor, static_cast<Real>(4.0/3.0)}) {
    CheckMagneticFloorCase(0.0, pressure);
    CheckMagneticFloorCase(MakeCglEOS().bfloor, pressure);
  }
  CheckUnresolvedInternalEnergyFloor();
  CheckWideDynamicRangeRoundTrip();
  CheckIsotropicExactRoundTrip();
  CheckExtremeLogRatioFloor("large-positive-log-ratio", 1000.0, 1.0e-11,
                            2.0, 2.0);
  CheckExtremeLogRatioFloor("large-negative-log-ratio", -1000.0, 1.0e-11,
                            6.0, 1.0e-12);
  CheckDensityFloorRatios();
  CheckInvalidAnisotropy();
  CheckInvalidDensity();
  CheckMagnetizationDensityFloor();
  CheckNegativeInternalEnergyFloor();
  CheckCollisionMap();

  std::cout << "CGL C2P pressure-floor checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglC2PPressureFloorChecks();
}
