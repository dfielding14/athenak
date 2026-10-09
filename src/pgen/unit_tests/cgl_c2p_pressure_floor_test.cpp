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
#include <fstream>
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
  eos.sigma_max = std::numeric_limits<float>::max();
  eos.firehose_threshold = 2.0;
  eos.mirror_threshold = 1.0;
  eos.firehose_backup_factor = 1.0;
  eos.mirror_backup_factor = 2.0;
  eos.limiter_backup_nu = 1.0e10;
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
                      const HydPrim1D &first, const EOS_Data &eos,
                      const bool magnetic_moment) {
  const MHDCons1D corrected = u;
  HydPrim1D recovered{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, recovered, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, recovered, dfloor_used, efloor_used,
                     tfloor_used, bfloor_used);
  }
  Require(label + " no repeated density or pressure floor", !dfloor_used && !efloor_used);
  Require(label + " unchanged conserved fields",
          u.d == corrected.d && u.mx == corrected.mx && u.my == corrected.my &&
          u.mz == corrected.mz && u.e == corrected.e && u.mu == corrected.mu &&
          u.bx == corrected.bx && u.by == corrected.by && u.bz == corrected.bz);
  Require(label + " unchanged primitive fields",
          recovered.d == first.d && recovered.vx == first.vx &&
          recovered.vy == first.vy && recovered.vz == first.vz &&
          recovered.e == first.e && recovered.pp == first.pp);
  const Real kinetic = 0.5*(1.0/u.d)*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  const Real magnetic = 0.5*(SQR(u.bx) + SQR(u.by) + SQR(u.bz));
  const Real tolerance = kSinglePrecision ? 32.0*std::numeric_limits<Real>::epsilon()
                                        : 1.0e-14;
  RequireRelativeClose(label + " internal energy identity", u.e - kinetic - magnetic,
                       0.5*first.e + first.pp, tolerance);
}

void CheckDensityRepair(const std::string &label, MHDCons1D u, const EOS_Data &eos,
                         const Real expected_density, const Real initial_parallel,
                         const Real initial_perp, const bool magnetic_moment) {
  const MHDCons1D original = u;
  const Real bsqr = SQR(u.bx) + SQR(u.by) + SQR(u.bz);
  const Real bmag = std::sqrt(bsqr);
  const Real kinetic = 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/expected_density;
  const Real internal = original.e - kinetic - 0.5*bsqr;
  const Real ratio = initial_perp/initial_parallel;
  const Real expected_parallel = magnetic_moment ? 2.0*(internal - initial_perp)
                                                : internal/(0.5 + ratio);
  const Real expected_perp = magnetic_moment ? initial_perp : ratio*expected_parallel;
  const Real tolerance = kSinglePrecision ? 2.0e-5 : 1.0e-12;
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  }
  Require(label + " density floor flag", dfloor_used);
  Require(label + " no other floors", !efloor_used && !tfloor_used && !bfloor_used);
  Require(label + " density floor value",
          u.d == expected_density && w.d == expected_density);
  Require(label + " preserved total energy", u.e == original.e);
  Require(label + " preserved momentum",
          u.mx == original.mx && u.my == original.my && u.mz == original.mz);
  Require(label + " preserved magnetic field",
          u.bx == original.bx && u.by == original.by && u.bz == original.bz);
  RequireClose(label + " vx", w.vx, original.mx/expected_density);
  RequireClose(label + " vy", w.vy, original.my/expected_density);
  RequireClose(label + " vz", w.vz, original.mz/expected_density);
  RequireRelativeClose(label + " p_parallel", w.e, expected_parallel, tolerance);
  RequireRelativeClose(label + " p_perp", w.pp, expected_perp, tolerance);
  if (magnetic_moment) {
    Require(label + " preserved magnetic moment", u.mu == original.mu);
  } else {
    RequireRelativeClose(label + " preserved pressure ratio", w.pp/w.e, ratio, tolerance);
    const Real expected_anisotropy = expected_density*
        (std::log(ratio) + 2.0*std::log(expected_density) - 3.0*std::log(bmag));
    RequireClose(label + " repaired conserved anisotropy", u.mu, expected_anisotropy);
    Require(label + " changed conserved anisotropy", u.mu != original.mu);
  }
  CheckIdempotence(label, u, w, eos, magnetic_moment);
}

void CheckDensityFloorRatios() {
  EOS_Data eos = MakeCglEOS();
  eos.dfloor = 1.0e-4;
  eos.pfloor = 1.0e-12;
  const Real densities[] = {5.0e-5, 1.0e-5, 1.0e-6, 1.0e-8};
  const Real ratios[] = {1.0, 0.3, 4.0};
  for (const Real density : densities) {
    for (const Real ratio : ratios) {
      MHDCons1D u{};
      u.d = density;
      u.mx = 0.2*density;
      u.my = -0.1*density;
      u.mz = 0.3*density;
      u.bx = 1.0;
      u.e = 1.0e-3*(0.5 + ratio)
          + 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/density + 0.5;
      u.mu = density*(std::log(ratio) + 2.0*std::log(density));
      CheckDensityRepair("density floor", u, eos, eos.dfloor, 1.0e-3,
                          1.0e-3*ratio, false);
    }
  }
}

void CheckMagneticMomentDensityFloor() {
  EOS_Data eos = MakeCglEOS();
  eos.dfloor = 2.0*kDensity;
  eos.pfloor = 1.0e-12;
  MHDCons1D u = MakeConserved(1.0, 2.0);
  u.mu = 2.0/std::sqrt(SQR(kBx) + SQR(kBy) + SQR(kBz));
  CheckDensityRepair("moment density floor", u, eos, eos.dfloor, 1.0, 2.0, true);
}

void CheckNonpositiveDensityFloor() {
  EOS_Data eos = MakeCglEOS();
  eos.dfloor = 1.0;
  eos.pfloor = 1.0e-12;
  for (const Real density : {0.0, -1.0}) {
    MHDCons1D u{};
    u.d = density;
    u.bx = 1.0;
    u.e = 3.5;
    HydPrim1D w{};
    bool dfloor_used = false, efloor_used = false;
    bool tfloor_used = false, bfloor_used = false;
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
    Require("nonpositive density floor flag", dfloor_used);
    Require("nonpositive density no other floors",
            !efloor_used && !tfloor_used && !bfloor_used);
    Require("nonpositive density floor value", u.d == 1.0 && w.d == 1.0);
    RequireClose("nonpositive density p_parallel", w.e, 2.0);
    RequireClose("nonpositive density p_perp", w.pp, 2.0);
    Require("nonpositive density preserved energy", u.e == 3.5);
    Require("nonpositive density canonical anisotropy", u.mu == 0.0);
    CheckIdempotence("nonpositive density", u, w, eos, false);
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
    u.mu = magnetic_moment ? 1.0
        : u.d*(std::log(2.0) + 2.0*std::log(u.d) - 3.0*std::log(u.bx));
    CheckDensityRepair(magnetic_moment ? "moment magnetization floor"
                                      : "anisotropy magnetization floor",
                        u, eos, 2.0, 1.0, 2.0, magnetic_moment);
  }
}

void CheckFloorCase(const std::string &label, const Real initial_parallel,
                    const Real initial_perp, const Real expected_parallel,
                    const Real expected_perp, const bool magnetic_moment = false) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(initial_parallel, initial_perp);
  const Real bmag = std::sqrt(SQR(kBx) + SQR(kBy) + SQR(kBz));
  if (magnetic_moment) u.mu = initial_perp/bmag;
  HydPrim1D w{};
  bool dfloor_used = false;
  bool efloor_used = false;
  bool tfloor_used = false;
  bool bfloor_used = false;

  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  }

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
  if (magnetic_moment) {
    RequireClose(label + " magnetic moment", u.mu, expected_perp/bmag);
  } else {
    RequireClose(label + " conserved anisotropy", u.mu,
                 CGLConservedAnisotropy(
                     kDensity, expected_parallel, expected_perp, bmag));
  }
  CheckIdempotence(label, u, w, eos, magnetic_moment);
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

void CheckPressureRatioUnderflowBoundaries() {
  const int large_exponent = kSinglePrecision ? 80 : 600;
  const Real large = std::ldexp(static_cast<Real>(1.0), large_exponent);
  const int double_exponents[] = {1021, 1022, 1023, 1073, 1074, 1075, 1200};
  const int float_exponents[] = {125, 126, 127, 148, 149, 150, 160};
  const Real tolerance = 2048.0*std::numeric_limits<Real>::epsilon();
  auto check = [&](const int exponent, const bool pressure_floor) {
    const Real small = std::ldexp(static_cast<Real>(1.0), large_exponent - exponent);
    EOS_Data eos = MakeCglEOS();
    eos.pfloor = pressure_floor
        ? std::ldexp(static_cast<Real>(1.0), -large_exponent) : small/16.0;
    for (const bool perpendicular_large : {false, true}) {
      const std::string label = "ratio 2^-" + std::to_string(exponent)
          + (perpendicular_large ? " perpendicular-large" : " parallel-large")
          + (pressure_floor ? " pressure-floor" : " healthy");
      const Real initial_parallel = perpendicular_large ? small : large;
      const Real initial_perp = perpendicular_large ? large : small;
      MHDCons1D u{};
      u.d = u.bx = 1.0;
      u.mx = 0.125;
      u.my = -0.25;
      u.mz = 0.5;
      u.e = 0.5*initial_parallel + initial_perp
          + 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz)) + 0.5;
      u.mu = (perpendicular_large ? exponent : -exponent)*std::log(2.0);
      const MHDCons1D original = u;
      HydPrim1D w{};
      bool dfloor_used = false, efloor_used = false;
      bool tfloor_used = false, bfloor_used = false;
      SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
      Require(label + " pressure-floor flag", efloor_used == pressure_floor);
      Require(label + " no other floors", !dfloor_used && !tfloor_used && !bfloor_used);
      RequireRelativeClose(label + " p_parallel", w.e,
                           std::max(initial_parallel, eos.pfloor), tolerance);
      RequireRelativeClose(label + " p_perp", w.pp,
                           std::max(initial_perp, eos.pfloor), tolerance);
      Require(label + " admissible pressures", w.e >= eos.pfloor && w.pp >= eos.pfloor);
      Require(label + " preserved primary fields",
              u.d == original.d && u.mx == original.mx && u.my == original.my &&
              u.mz == original.mz && u.bx == original.bx && u.by == original.by &&
              u.bz == original.bz);
      if (!pressure_floor) {
        Require(label + " unchanged energy and anisotropy",
                u.e == original.e && u.mu == original.mu);
      }
      CheckIdempotence(label, u, w, eos, false);
    }
  };
  for (int i = 0; i < 7; ++i) {
    check(kSinglePrecision ? float_exponents[i] : double_exponents[i], false);
  }
  check(kSinglePrecision ? 180 : 1220, true);
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

  CheckIdempotence(label, u, w, eos, false);
}

void CheckNegativeInternalEnergyFloor(const bool magnetic_moment = false) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(2.0, 2.0);
  if (magnetic_moment) u.mu = 2.0/std::sqrt(SQR(kBx) + SQR(kBy) + SQR(kBz));
  u.e = KineticEnergy() + MagneticEnergy() - 1.0;

  HydPrim1D w{};
  bool dfloor_used = false;
  bool efloor_used = false;
  bool tfloor_used = false;
  bool bfloor_used = false;
  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  }

  Require("negative-eint pressure floor flag", efloor_used);
  RequireClose("negative-eint p_parallel", w.e, kPressureFloor);
  const Real expected_perp = magnetic_moment ? 2.0 : kPressureFloor;
  RequireClose("negative-eint p_perp", w.pp, expected_perp);
  Require("negative-eint total energy is finite", std::isfinite(u.e));
  Require("negative-eint conserved anisotropy is finite", std::isfinite(u.mu));
  RequireClose("negative-eint corrected total energy", u.e,
               0.5*kPressureFloor + expected_perp + KineticEnergy() + MagneticEnergy());
  CheckIdempotence("negative internal energy", u, w, eos, magnetic_moment);
}

void CheckMagneticFloorCase(const Real bmag, const Real pressure,
                            const bool magnetic_moment) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u = MakeConserved(2.0, 1.0);
  u.bx = bmag;
  u.by = u.bz = 0.0;
  if (magnetic_moment) u.mu = pressure/eos.bfloor;
  const Real expected_energy = 1.5*pressure + KineticEnergy() + 0.5*SQR(bmag);
  u.e = expected_energy;
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  }
  Require("magnetic floor flag", bfloor_used);
  Require("magnetic floor isotropy", w.e == w.pp);
  RequireClose("magnetic floor pressure", w.e, pressure);
  RequireClose("magnetic floor total energy", u.e, expected_energy);
  CheckIdempotence("magnetic floor", u, w, eos, magnetic_moment);
}

void CheckUnresolvedInternalEnergyFloor(const bool magnetic_moment) {
  const EOS_Data eos = MakeCglEOS();
  MHDCons1D u{};
  u.d = 1.0;
  u.mx = 1.0e12;
  u.bx = 1.0e8;
  const Real kinetic = 0.5*SQR(u.mx);
  const Real magnetic = 0.5*SQR(u.bx);
  u.e = std::nextafter(static_cast<Real>(kinetic + magnetic), static_cast<Real>(0.0));
  u.mu = magnetic_moment ? 1.0/u.bx : -3.0*std::log(u.bx);
  const MHDCons1D original = u;
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  if (magnetic_moment) {
    SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                     tfloor_used, bfloor_used);
  } else {
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
  }
  Require("unresolved internal energy uses pressure floor", efloor_used);
  Require("unresolved internal energy no other floors",
          !dfloor_used && !tfloor_used && !bfloor_used);
  Require("representable floor pressures", std::isfinite(w.e) && std::isfinite(w.pp)
          && w.e >= eos.pfloor && w.pp >= eos.pfloor);
  // A rounded K+M may already leave positive thermal energy when K and M are subtracted.
  Require("representable total energy", std::isfinite(u.e) && u.e > original.e);
  Require("representable anisotropy or moment", std::isfinite(u.mu));
  Require("unresolved floor preserved primary fields",
          u.d == original.d && u.mx == original.mx && u.my == original.my &&
          u.mz == original.mz && u.bx == original.bx && u.by == original.by &&
          u.bz == original.bz);
  CheckIdempotence("unresolved internal energy", u, w, eos, magnetic_moment);
}

void CheckMagneticMomentDivisionRounding() {
  const EOS_Data eos = MakeCglEOS();
  const Real bmag = kSinglePrecision ? 41.0 : 49.0;
  const Real rounded_moment = eos.pfloor/bmag;
  Require("moment division fixture rounds pressure down",
          rounded_moment*bmag < eos.pfloor);
  MHDCons1D u{};
  u.d = 1.0;
  u.bx = bmag;
  u.mu = 0.25/bmag;
  u.e = 3.0 + 0.5*SQR(bmag);
  HydPrim1D w{};
  bool dfloor_used = false, efloor_used = false;
  bool tfloor_used = false, bfloor_used = false;
  SingleC2P_CGLMHDFromMagneticMoment(u, eos, w, dfloor_used, efloor_used,
                                   tfloor_used, bfloor_used);
  Require("moment division pressure floor flag", efloor_used);
  Require("moment division no other floors",
          !dfloor_used && !tfloor_used && !bfloor_used);
  RequireClose("moment division p_parallel", w.e, 5.5);
  RequireClose("moment division p_perp", w.pp, eos.pfloor);
  RequireClose("moment division total energy", u.e, 3.75 + 0.5*SQR(bmag));
  Require("moment encoding clears pressure floor", u.mu > rounded_moment
          && u.mu*bmag >= eos.pfloor && w.pp >= eos.pfloor);
  CheckIdempotence("moment division", u, w, eos, true);
}

// Keep the physical oracle separate from C2P and the shared threshold helpers.
void CheckCollisionMap() {
  EOS_Data eos = MakeCglEOS();
  eos.firehose_threshold = 1.4;
  const Real mean_tolerance = kSinglePrecision ? 8.0*std::numeric_limits<Real>::epsilon()
                                              : 1.0e-14;
  for (const Real nudt : {Real(0.0), Real(1.0), Real(1.0e10)}) {
    eos.lim_coll = nudt;
    Real previous = -std::numeric_limits<Real>::max();
    for (int n = 0; n <= 80; ++n) {
      const Real delta = -2.0 + 0.05*n;
      MHDPrim1D w{};
      w.bx = 1.0;
      w.e = 10.0 - TWO_3RDS*delta;
      w.pp = 10.0 + ONE_3RD*delta;
      const MHDPrim1D original = w;
      const Real initial = w.pp - w.e;
      const Real piso = ONE_3RD*w.e + TWO_3RDS*w.pp;
      const Real threshold = initial < -0.7 ? -0.7 : 0.5;
      const Real expected = initial < -0.7 || initial > 0.5
          ? threshold + (initial - threshold)/(1.0 + nudt) : initial;
      SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
      RequireClose("limiter backward Euler", w.pp-w.e, expected);
      Require("monotone limiter", w.pp-w.e >= previous);
      previous = w.pp-w.e;
      RequireRelativeClose("rate mean pressure", ONE_3RD*w.e + TWO_3RDS*w.pp,
                           piso, mean_tolerance);
      if (nudt == 0.0) {
        Require("zero-rate identity", std::memcmp(&w, &original, sizeof(w)) == 0);
      }
      SingleCollWalls_CGLMHD(w, eos, true);
      const MHDPrim1D projected = w;
      Require("wall bounds", w.pp-w.e >= -0.7 && w.pp-w.e <= 1.0);
      SingleCollWalls_CGLMHD(w, eos, true);
      Require("bitwise idempotent walls", std::memcmp(&w, &projected, sizeof(w)) == 0);
    }
  }

  // Exact soft-threshold equality adds no relaxation.
  eos.firehose_threshold = 2.0;
  eos.nu_coll = 0.0;
  eos.lim_coll = 1.0e10;
  for (const bool firehose : {false, true}) {
    MHDPrim1D w{};
    w.bx = 1.0;
    w.e = firehose ? 2.0 : 1.0;
    w.pp = firehose ? 1.0 : 1.5;
    const MHDPrim1D original = w;
    SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
    Require("soft-threshold equality identity",
            std::memcmp(&w, &original, sizeof(w)) == 0);
  }
  eos.firehose_threshold = 1.4;

  // Background decay must precede the soft-threshold test and full-step map.
  for (const Real decay : {Real(1.0), Real(0.5), Real(0.25)}) {
    MHDPrim1D w{};
    w.bx = 1.0;
    w.e = 2.0;
    w.pp = 4.0;
    eos.nu_coll = -std::log(decay);
    eos.lim_coll = 1.0;
    const Real background = 2.0*decay;
    const Real expected = background > 0.5 ? 0.5 + (background-0.5)/2.0 : background;
    SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
    RequireClose("background then full limiter kick", w.pp-w.e, expected);
    RequireClose("background and limiter preserve energy", 0.5*w.e+w.pp, 5.0);
  }

  // Products of finite rates and timesteps may overflow in the stiff limit.
  for (const bool background : {false, true}) {
    MHDPrim1D w{};
    w.bx = 1.0;
    w.e = 2.0;
    w.pp = 4.0;
    eos.nu_coll = background ? std::numeric_limits<Real>::max() : 0.0;
    eos.lim_coll = std::numeric_limits<Real>::max();
    SingleCollRates_CGLMHD(w, eos, 2.0, true, true);
    RequireClose("overflowing rate product has finite limit", w.pp-w.e,
                 background ? 0.0 : 0.5);
    RequireClose("stiff rate preserves energy", 0.5*w.e+w.pp, 5.0);
  }

  // Backup walls ignore soft flags; the fluid firehose wall ignores backup.
  eos.mlim = eos.flim = false;
  for (const bool backup : {false, true}) {
    for (const Real delta : {Real(-2.0), Real(2.0)}) {
      MHDPrim1D w{};
      w.bx = 1.0;
      w.e = 3.0 - TWO_3RDS*delta;
      w.pp = 3.0 + ONE_3RD*delta;
      SingleCollWalls_CGLMHD(w, eos, backup);
      const Real expected = backup ? (delta < 0.0 ? -0.7 : 1.0)
                                   : (delta < 0.0 ? -1.0 : delta);
      RequireClose("independent wall flags", w.pp-w.e, expected);
    }
  }
  std::cout << "CGL collision-map checks passed" << std::endl;
}

// Exercise the real grid collision operator, then decode in a separate kernel.
// These doubles are the three cells that failed a strict HIP restart at t=10.
void CheckProductionWallEncoding(MeshBlockPack *pack, const std::string &filename) {
  Require("production wall fixture uses double", !kSinglePrecision);
  constexpr int ncase = 6;
  double cases[ncase][9] = {
    {1.0580936533619145, 0.4998698938877745, 0.162893050030097,
     -0.03657178263370929, 11.505868899072333, -0.1425980922790349,
     0.23714419782252238, 0.03332376975719438, 1.0067785045133788},
    {1.068962309818397, -0.2853674066747329, 1.1165715708919983,
     -0.47723398025244956, 12.294172720527772, -0.5135541903454163,
     0.4807561523963997, -0.5410872057353117, 0.9013372287919911},
    {1.045423160356516, -0.3245337512214358, 0.894759895491961,
     -0.4707708062214894, 11.851465659314885, -0.501664065376218,
     0.7615518523764249, -0.45361804497395064, 0.7156273634530386},
    // Isotropic interior, exact fluid wall, and a genuinely outside-wall state.
    {1., 0., 0., 0., 8., 0., 0., 0., 1.},
    {1., 0., 0., 0., 4., std::log(2.)-std::log(3.), 0., 0., 1.},
    {1., 0., 0., 0., 3.5, -std::log(4.), 0., 0., 1.}
  };
  const int nmb = pack->nmb_thispack;
  auto *peos = pack->pmhd->peos;
  const EOS_Data eos = peos->eos_data;
  const bool backup = eos.backup_lim;
  Require("active production wall fixture", eos.is_cgl && !eos.passive);
  Real expected_delta[ncase]{};
  Real expected_thermal[ncase]{};
  if (backup) {
    Require("asymmetric backup fixture policy",
            eos.firehose_threshold == 1.4 && eos.firehose_backup_factor == 1.1 &&
            eos.mirror_threshold == 0.7 && eos.mirror_backup_factor == 1.3 &&
            eos.nu_coll == 0.3 && eos.lim_coll == 7.0);
    for (int c=0; c<ncase; ++c) {
      const Real rho = 1.7 + 0.13*c;
      const Real mx = 0.31 - 0.06*c, my = -0.23 + 0.04*c, mz = 0.11 + 0.02*c;
      const Real bx = 0.91 - 0.03*c, by = -0.28 + 0.02*c, bz = 0.37 + 0.01*c;
      const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
      // Independent numeric bounds; do not reuse the production threshold helpers.
      const Real lower = 1.1*(-0.5*1.4*bsqr), upper = 1.3*(0.5*0.7*bsqr);
      const Real delta = c == 0 ? -1.4*bsqr : c == 1 ? 1.2*bsqr :
                         c == 2 ? lower : c == 4 ? upper : c == 5 ? 0.44*bsqr : 0.0;
      const Real thermal = 6.0 + 0.2*c;
      const Real ppar = (2.0/3.0)*(thermal - delta);
      const Real pperp = (2.0/3.0)*thermal + delta/3.0;
      cases[c][0] = rho; cases[c][1] = mx; cases[c][2] = my; cases[c][3] = mz;
      cases[c][4] = thermal + 0.5/rho*(SQR(mx) + SQR(my) + SQR(mz)) + 0.5*bsqr;
      cases[c][5] = rho*(std::log(pperp/ppar) + 2.0*std::log(rho)
                        - 3.0*std::log(std::sqrt(bsqr)));
      cases[c][6] = bx; cases[c][7] = by; cases[c][8] = bz;
      expected_delta[c] = std::min(upper, std::max(lower, delta));
      expected_thermal[c] = thermal;
      if (c == 5) {
        // This interior state sees background relaxation, then a finite soft kick.
        const Real after_background = delta*std::exp(-0.3*0.2);
        Require("finite-rate fixture activates soft mirror relaxation",
                after_background > 0.35*bsqr && delta < upper);
        expected_delta[c] = 0.35*bsqr + (after_background - 0.35*bsqr)/(1.0 + 7.0*0.2);
      }
    }
  }
  DvceArray5D<Real> cons("wall_cons", nmb, 6, 1, 1, ncase);
  DvceArray5D<Real> prim("wall_prim", nmb, 6, 1, 1, ncase);
  DvceArray5D<Real> decoded("wall_decoded", nmb, 6, 1, 1, ncase);
  DvceArray5D<Real> bcc("wall_bcc", nmb, 3, 1, 1, ncase);
  Kokkos::View<int*> hard_bad("wall_bad", nmb*ncase);
  auto copy_to_host = [](const DvceArray5D<Real> &view) {
    auto host = Kokkos::create_mirror(view);
    Kokkos::deep_copy(host, view);
    return host;
  };
  auto initial = Kokkos::create_mirror(cons);
  auto initial_b = Kokkos::create_mirror(bcc);
  for (int m=0; m<nmb; ++m) {
    for (int c=0; c<ncase; ++c) {
      for (int v=0; v<6; ++v) initial(m,v,0,0,c) = cases[c][v];
      for (int v=0; v<3; ++v) initial_b(m,v,0,0,c) = cases[c][v+6];
    }
  }
  Kokkos::deep_copy(cons, initial);
  Kokkos::deep_copy(bcc, initial_b);
  // The same decoder is launched independently before and after the real operator.
  auto decode = [&](DvceArray5D<Real> dest) {
    Kokkos::parallel_for("production_wall_independent_c2p", nmb*ncase,
    KOKKOS_LAMBDA(const int index) {
      const int m = index/ncase, c = index%ncase;
      MHDCons1D u{cons(m,IDN,0,0,c), cons(m,IM1,0,0,c), cons(m,IM2,0,0,c),
                  cons(m,IM3,0,0,c), cons(m,IEN,0,0,c), cons(m,IAN,0,0,c),
                  bcc(m,IBX,0,0,c), bcc(m,IBY,0,0,c), bcc(m,IBZ,0,0,c)};
      HydPrim1D w{};
      bool df=false, ef=false, tf=false, bf=false;
      SingleC2P_CGLMHD(u, eos, w, df, ef, tf, bf);
      if (df || ef || tf || bf) Kokkos::abort("Unexpected wall fixture floor");
      dest(m,IDN,0,0,c)=w.d; dest(m,IVX,0,0,c)=w.vx;
      dest(m,IVY,0,0,c)=w.vy; dest(m,IVZ,0,0,c)=w.vz;
      dest(m,IPR,0,0,c)=w.e; dest(m,IPP,0,0,c)=w.pp;
      const Real bsqr = SQR(u.bx) + SQR(u.by) + SQR(u.bz);
      const Real lower = backup ? 1.1*(-0.5*1.4*bsqr) : -bsqr;
      const Real upper = backup ? 1.3*(0.5*0.7*bsqr) :
                                 std::numeric_limits<Real>::max();
      hard_bad(index) = !Kokkos::isfinite(w.e) || !Kokkos::isfinite(w.pp) ||
                        !(w.e > 0.0) || !(w.pp > 0.0) ||
                        w.pp-w.e < lower || w.pp-w.e > upper;
    });
  };
  decode(prim);
  if (backup) {
    peos->Collisions(prim, bcc, cons, 0.2, CGLCollisionMode::full, 5, 5, 0, 0, 0, 0);
  }
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0, ncase-1, 0, 0, 0, 0);
  decode(decoded);
  auto corrected = copy_to_host(cons);
  auto cached = copy_to_host(prim);
  auto fresh = copy_to_host(decoded);
  auto final_b = copy_to_host(bcc);
  auto bad = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hard_bad);
  for (int m=0; m<nmb; ++m) {
    for (int c=0; c<ncase; ++c) {
      for (int v=0; v<5; ++v) {
        Require("wall conserves rho/momentum/total E bitwise",
                std::memcmp(&corrected(m,v,0,0,c), &initial(m,v,0,0,c),
                            sizeof(Real)) == 0);
      }
      Require("production wall survives independent C2P", bad(m*ncase+c) == 0);
      if (backup) {
        RequireClose("encoded backup/rate anisotropy",
                     fresh(m,IPP,0,0,c)-fresh(m,IPR,0,0,c), expected_delta[c]);
        RequireClose("encoded backup/rate parallel pressure", fresh(m,IPR,0,0,c),
                     (2.0/3.0)*(expected_thermal[c] - expected_delta[c]));
        RequireClose("encoded backup/rate perpendicular pressure", fresh(m,IPP,0,0,c),
                     (2.0/3.0)*expected_thermal[c] + expected_delta[c]/3.0);
      } else if (c != 3) {
        const Real bsqr = SQR(initial_b(m,IBX,0,0,c)) + SQR(initial_b(m,IBY,0,0,c)) +
                          SQR(initial_b(m,IBZ,0,0,c));
        RequireClose("encoded wall retains boundary anisotropy",
                     fresh(m,IPP,0,0,c)-fresh(m,IPR,0,0,c), -bsqr);
      }
      for (int v : {IPR, IPP}) {
        Require("production wall cached pressure is canonical",
                cached(m,v,0,0,c) == fresh(m,v,0,0,c));
      }
    }
    Require("interior conservative A stays untouched",
            std::memcmp(&corrected(m,IAN,0,0,3), &initial(m,IAN,0,0,3),
                        sizeof(Real)) == 0);
  }
  Require("wall leaves magnetic field bitwise unchanged",
          std::memcmp(initial_b.data(), final_b.data(),
                      final_b.size()*sizeof(Real)) == 0);
  // Serialize the actual corrected conserved state, restore into the device view,
  // and check the independent decoder again. No primitive cache is serialized.
  const std::size_t bytes = corrected.size()*sizeof(Real);
  {
    std::ofstream stream(filename, std::ios::binary);
    stream.write(reinterpret_cast<const char*>(corrected.data()), bytes);
    Require("write corrected wall state", stream.good());
  }
  auto restored = Kokkos::create_mirror(cons);
  {
    std::ifstream stream(filename, std::ios::binary);
    stream.read(reinterpret_cast<char*>(restored.data()), bytes);
    Require("read corrected wall state", stream.good());
    Require("wall state has exact serialized length", stream.peek() == EOF);
  }
  Require("serialized wall state bitwise identity",
          std::memcmp(corrected.data(), restored.data(), bytes) == 0);
  Kokkos::deep_copy(cons, restored);
  decode(prim);
  auto restarted = copy_to_host(prim);
  Require("serialized state independently recovers identical primitives",
          std::memcmp(fresh.data(), restarted.data(), fresh.size()*sizeof(Real)) == 0);
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0, ncase-1, 0, 0, 0, 0);
  auto repeated = copy_to_host(cons);
  Require("production wall encoded fixed point",
          std::memcmp(corrected.data(), repeated.data(), bytes) == 0);
  std::cout << "CGL production wall encoding/restart checks passed: "
            << ncase << " states per block" << std::endl;
  if (backup) {
    std::cout << "CGL asymmetric backup wall/full-rate checks passed" << std::endl;
  }

  // Reproduce the LF sweep exit. Conversion must refresh the pressure cache
  // from stored E and A before the following wall operation can skip a cell.
  peos->CGLAnisotropyToMagneticMoment(cons, bcc, 0, ncase-1, 0, 0, 0, 0);
  peos->CGLRefreshPrimFromMagneticMoment(cons, bcc, prim, 0, ncase-1, 0, 0, 0, 0);
  peos->CGLMagneticMomentToAnisotropy(cons, bcc, prim, 0, ncase-1, 0, 0, 0, 0);
  auto converted = copy_to_host(cons);
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0, ncase-1, 0, 0, 0, 0);
  decode(decoded);
  auto transition_cons = copy_to_host(cons);
  auto transition_cached = copy_to_host(prim);
  auto transition_fresh = copy_to_host(decoded);
  auto transition_bad = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hard_bad);
  for (int m=0; m<nmb; ++m) {
    for (int c=0; c<ncase; ++c) {
      const bool same_pressures =
          transition_cached(m,IPR,0,0,c) == transition_fresh(m,IPR,0,0,c) &&
          transition_cached(m,IPP,0,0,c) == transition_fresh(m,IPP,0,0,c);
      if (transition_bad(m*ncase+c) != 0 || !same_pressures) {
        const Real bsqr = SQR(initial_b(m,IBX,0,0,c)) + SQR(initial_b(m,IBY,0,0,c)) +
                          SQR(initial_b(m,IBZ,0,0,c));
        std::cout.precision(std::numeric_limits<Real>::max_digits10);
        std::cout << "LF conversion wall case " << c << ": cached delta="
                  << transition_cached(m,IPP,0,0,c)-transition_cached(m,IPR,0,0,c)
                  << ", stored delta="
                  << transition_fresh(m,IPP,0,0,c)-transition_fresh(m,IPR,0,0,c)
                  << ", fluid wall=" << -bsqr << std::endl;
      }
      Require("LF conversion wall survives independent C2P",
              transition_bad(m*ncase+c) == 0);
      Require("LF conversion wall cached pressure is canonical", same_pressures);
      for (int v=0; v<5; ++v) {
        Require("LF conversion wall conserves rho/momentum/total E bitwise",
                std::memcmp(&transition_cons(m,v,0,0,c), &converted(m,v,0,0,c),
                            sizeof(Real)) == 0);
      }
    }
  }
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0, ncase-1, 0, 0, 0, 0);
  auto transition_repeated = copy_to_host(cons);
  Require("LF conversion wall encoded fixed point",
          std::memcmp(transition_cons.data(), transition_repeated.data(), bytes) == 0);
  std::cout << "CGL LF conversion wall encoding checks passed: "
            << ncase << " states per block" << std::endl;
}

} // namespace

void RunCglC2PPressureFloorChecks() {
  CheckFloorCase("parallel-only", 0.25, 2.0, kPressureFloor, 2.0);
  CheckFloorCase("perpendicular-only", 2.0, 0.25, 2.0, kPressureFloor);
  CheckFloorCase("both-pressure", 0.25, 0.5, kPressureFloor, kPressureFloor);
  CheckFloorCase("moment parallel-only", 0.25, 2.0, kPressureFloor, 2.0, true);
  CheckFloorCase("moment perpendicular-only", 2.0, 0.25, 2.0, kPressureFloor, true);
  CheckFloorCase("moment both-pressure", 0.25, 0.5, kPressureFloor, kPressureFloor, true);
  CheckWideDynamicRangeRoundTrip();
  CheckPressureRatioUnderflowBoundaries();
  CheckIsotropicExactRoundTrip();
  CheckExtremeLogRatioFloor("large-positive-log-ratio", 1000.0, 1.0e-11,
                            1.0e-12, 3.0);
  CheckExtremeLogRatioFloor("large-negative-log-ratio", -1000.0, 1.0e-11,
                            6.0, 1.0e-12);
  CheckNegativeInternalEnergyFloor();
  CheckNegativeInternalEnergyFloor(true);
  for (const bool magnetic_moment : {false, true}) {
    for (const Real pressure : {kPressureFloor, static_cast<Real>(4.0/3.0)}) {
      CheckMagneticFloorCase(0.0, pressure, magnetic_moment);
      CheckMagneticFloorCase(MakeCglEOS().bfloor, pressure, magnetic_moment);
    }
    CheckUnresolvedInternalEnergyFloor(magnetic_moment);
  }
  CheckMagneticMomentDivisionRounding();
  CheckDensityFloorRatios();
  CheckMagneticMomentDensityFloor();
  CheckNonpositiveDensityFloor();
  CheckMagnetizationDensityFloor();
  CheckCollisionMap();

  std::cout << "CGL C2P pressure-floor checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglC2PPressureFloorChecks();
  if (!kSinglePrecision) {
    CheckProductionWallEncoding(pmy_mesh_->pmb_pack, "wall_state.bin");
  }
  // Leave a finite uniform state for zero-step application startup.
  auto *pmhd = pmy_mesh_->pmb_pack->pmhd;
  Kokkos::deep_copy(pmhd->u0, 0.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IDN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 1.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IEN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 2.0);
  Kokkos::deep_copy(pmhd->b0.x1f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x2f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x3f, 1.0);
}
