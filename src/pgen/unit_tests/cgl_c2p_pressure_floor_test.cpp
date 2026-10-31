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
#include "eos/cgl_passive.hpp"
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

  // A tiny inward kick must raise the smaller pressure even when the change
  // cannot be represented by subtracting two large anisotropies or pressures.
  const Real small_pressure = kSinglePrecision ? Real(9.51513186e-7)
                                               : Real(1.4305114746093706e-6);
  const Real large_pressure = kSinglePrecision ? Real(0.000243587376)
                                               : Real(51539607552.0);
  const Real epsilon = std::numeric_limits<Real>::epsilon();
  for (const bool reverse : {false, true}) {
    for (const bool soft : {false, true}) {
      for (const Real rate : {epsilon/8, epsilon/4, epsilon/2, epsilon}) {
        MHDPrim1D w{};
        w.e = reverse ? large_pressure : small_pressure;
        w.pp = reverse ? small_pressure : large_pressure;
        w.bx = kSinglePrecision ? Real(0.00039366493) : Real(0.00042286396669162041);
        eos.nu_coll = soft ? 0.0 : rate;
        eos.lim_coll = soft ? rate : 0.0;
        const long double ppar = w.e, pperp = w.pp;
        const long double delta = pperp-ppar;
        const long double threshold = (reverse ? -0.5L*eos.firehose_threshold
                                               : 0.5L*eos.mirror_threshold)*SQR(w.bx);
        const long double exchange = soft
            ? (threshold-delta)*(static_cast<long double>(rate)/(1.0L+rate))
            : delta*std::expm1(-static_cast<long double>(rate));
        const Real expected_ppar = static_cast<Real>(ppar-(2.0L/3.0L)*exchange);
        const Real expected_pperp = static_cast<Real>(pperp+exchange/3.0L);
        SingleCollRates_CGLMHD(w, eos, 1.0, true, true);
        Require("tiny relaxation preserves minimum pressure",
                w.e >= small_pressure && w.pp >= small_pressure);
        RequireRelativeClose("tiny-kick parallel pressure", w.e, expected_ppar, 64*epsilon);
        RequireRelativeClose("tiny-kick perpendicular pressure", w.pp, expected_pperp,
                             64*epsilon);
        RequireRelativeClose("tiny-kick thermal energy", 0.5*w.e+w.pp,
                             static_cast<Real>(0.5L*ppar+pperp), 8*epsilon);
      }
    }
  }

  std::cout << "CGL collision-map checks passed" << std::endl;
}

// Independent pressure-form oracle from Stephen's 0e6d5e034 mirror correction.
// The backup branch replaces the soft branch; the combined map has jumps.
void CheckStephenLimiterMap() {
  EOS_Data eos = MakeCglEOS();
  eos.firehose_threshold = 1.4;
  eos.firehose_backup_factor = 2.0;
  const Real dt = 0.125;
  const Real offset = 32.0*std::numeric_limits<Real>::epsilon();
  const Real states[][2] = {{1.4, 0.7}, {0.5, 1.0}, {2.0, 1.0}, {1.0, 2.0},
                            {2.0, 0.5}, {0.5, 2.0}};
  const Real rates[][2] = {{0.0, 0.0}, {0.0, 1.0e10}, {4.0, 0.0},
                           {4.0, 13.0}, {4.0, 1.0e10}, {1.0e10, 13.0}};
  int checks = 0;
  for (const bool backup : {false, true}) {
    eos.backup_lim = backup;
    for (int flags = 0; flags < 4; ++flags) {
      const bool mirror = (flags & 1) != 0, firehose = (flags & 2) != 0;
      for (const Real background : {Real(0.0), Real(4.0), Real(20.0)}) {
        eos.nu_coll = background;
        for (const auto &rates_case : rates) {
          eos.lim_coll = rates_case[0];
          eos.limiter_backup_nu = rates_case[1];
          for (const auto &state : states) {
            for (const int side : {-1, 0, 1}) {
              MHDPrim1D w{};
              w.d = w.bx = 1.0;
              w.e = state[0];
              w.pp = state[1] + side*offset;
              const MHDPrim1D original = w;
              const long double thermal = 0.5L*w.e + w.pp;
              const long double piso = (w.e + 2.0L*w.pp)/3.0L;
              const long double delta = (static_cast<long double>(w.pp)-w.e)
                  * std::exp(-static_cast<long double>(background)*dt);
              long double perp = piso + delta/3.0L, parallel = perp-delta;
              long double target = 0.0L, rate = 0.0L;
              if (firehose && delta <= -0.7L) {
                target = backup && delta <= -1.0L ? -1.0L : -0.7L;
                rate = backup && delta <= -1.0L ? rates_case[1] : rates_case[0];
              } else if (mirror && delta >= 0.5L) {
                target = backup && delta >= 1.0L ? 1.0L : 0.5L;
                rate = backup && delta >= 1.0L ? rates_case[1] : rates_case[0];
              }
              const long double nudt = rate*dt;
              const long double expected_perp =
                  (3.0L*perp + nudt*(2.0L*perp+parallel+target))/(3.0L+3.0L*nudt);
              const long double expected_parallel =
                  (3.0L*parallel + nudt*(2.0L*perp+parallel-2.0L*target))
                  /(3.0L+3.0L*nudt);
              SingleCollRates_CGLMHD(w, eos, dt, mirror, firehose);
              RequireClose("Stephen limiter parallel pressure", w.e, expected_parallel);
              RequireClose("Stephen limiter perpendicular pressure", w.pp, expected_perp);
              RequireClose("Stephen limiter thermal energy", 0.5*w.e+w.pp, thermal);
              Require("Stephen limiter positive pressures", w.e > 0.0 && w.pp > 0.0);
              Require("Stephen limiter preserves other primitives",
                      w.d == original.d && w.vx == original.vx &&
                      w.vy == original.vy && w.vz == original.vz &&
                      w.bx == original.bx && w.by == original.by && w.bz == original.bz);
              if (background == 0.0 && rate == 0.0L) {
                Require("Stephen zero-rate identity",
                        std::memcmp(&w, &original, sizeof(w)) == 0);
              }
              ++checks;
            }
          }
        }
      }
    }
  }
  // A configurable soft onset can lie beyond the fluid-capped backup target.
  // The backup still acts first, as it does in the LF scattering coefficient.
  eos.nu_coll = 0.0;
  eos.firehose_threshold = 4.0;
  eos.lim_coll = 64.0;
  eos.limiter_backup_nu = 8.0;
  MHDPrim1D capped{};
  capped.d = capped.bx = 1.0;
  capped.e = 2.0;
  capped.pp = 0.5;
  const MHDPrim1D original = capped;
  eos.backup_lim = false;
  SingleCollRates_CGLMHD(capped, eos, dt, false, true);
  Require("inactive soft limiter with backup disabled",
          std::memcmp(&capped, &original, sizeof(capped)) == 0);
  eos.backup_lim = true;
  SingleCollRates_CGLMHD(capped, eos, dt, false, true);
  RequireClose("capped backup parallel pressure", capped.e, 11.0L/6.0L);
  RequireClose("capped backup perpendicular pressure", capped.pp, 7.0L/12.0L);
  RequireClose("capped backup thermal energy", 0.5*capped.e+capped.pp, 1.5L);
  std::cout << "CGL Stephen two-tier limiter checks passed: " << checks << std::endl;
}

// Finite backups must survive conserved-state reload and LF representation changes.
void CheckProductionFiniteBackup(MeshBlockPack *pack, const std::string &filename) {
  auto *peos = pack->pmhd->peos;
  const EOS_Data saved_eos = peos->eos_data;
  const int nmb = pack->nmb_thispack;
  constexpr int ncase = 4;
  constexpr Real dt = 2.0e-4;
  DvceArray5D<Real> cons("finite_backup_cons", nmb, 6, 1, 1, ncase);
  DvceArray5D<Real> prim("finite_backup_prim", nmb, 6, 1, 1, ncase);
  DvceArray5D<Real> bcc("finite_backup_bcc", nmb, 3, 1, 1, ncase);
  auto copy_to_host = [](const DvceArray5D<Real> &view) {
    auto host = Kokkos::create_mirror(view);
    Kokkos::deep_copy(host, view);
    return host;
  };
  for (const bool passive : {false, true}) {
    EOS_Data eos = saved_eos;
    eos.passive = passive;
    eos.backup_lim = eos.mlim = eos.flim = true;
    eos.firehose_threshold = 1.4;
    eos.mirror_threshold = 1.0;
    eos.firehose_backup_factor = eos.mirror_backup_factor = 2.0;
    eos.nu_coll = 0.0;
    eos.lim_coll = 64.0;
    eos.limiter_backup_nu = 1.0e10;
    peos->eos_data = eos;
    auto initial = Kokkos::create_mirror(cons);
    auto initial_b = Kokkos::create_mirror(bcc);
    for (int m=0; m<nmb; ++m) for (int c=0; c<ncase; ++c) {
      const Real delta = c == 0 ? -1.2 : c == 1 ? 1.2 : c == 2 ? -0.8 : 0.8;
      const Real pp = 4.0 - TWO_3RDS*delta, pt = 4.0 + ONE_3RD*delta;
      initial(m,IDN,0,0,c) = 1.0;
      initial(m,IM1,0,0,c) = 0.1;
      initial(m,IM2,0,0,c) = -0.2;
      initial(m,IM3,0,0,c) = 0.3;
      const auto q = cgl::PassiveEncode(1.0, pp, pt, 1.0);
      initial(m,IEN,0,0,c) = passive ? q.j : 6.57;
      initial(m,IAN,0,0,c) = passive ? q.a : CGLConservedAnisotropy(1.0, pp, pt, 1.0);
      for (int v=0; v<3; ++v) initial_b(m,v,0,0,c) = v == IBZ ? 1.0 : 0.0;
    }
    Kokkos::deep_copy(cons, initial);
    Kokkos::deep_copy(bcc, initial_b);
    auto decode = [&]() {
      par_for("finite_backup_independent_c2p", DevExeSpace(), 0,nmb-1, 0,0, 0,0, 0,ncase-1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        MHDCons1D u{cons(m,IDN,k,j,i), cons(m,IM1,k,j,i), cons(m,IM2,k,j,i),
                    cons(m,IM3,k,j,i), cons(m,IEN,k,j,i), cons(m,IAN,k,j,i), 0., 0., 1.};
        HydPrim1D w{};
        if (passive) {
          w.d=u.d; w.vx=u.mx/u.d; w.vy=u.my/u.d; w.vz=u.mz/u.d;
          cgl::PassiveDecode(u.d, {u.e, u.mu}, 1.0, w.e, w.pp);
        } else {
          bool df=false, ef=false, tf=false, bf=false;
          SingleC2P_CGLMHD(u, eos, w, df, ef, tf, bf);
          if (df || ef || tf || bf) Kokkos::abort("Finite backup fixture needs a floor");
        }
        prim(m,IDN,k,j,i)=w.d; prim(m,IVX,k,j,i)=w.vx;
        prim(m,IVY,k,j,i)=w.vy; prim(m,IVZ,k,j,i)=w.vz;
        prim(m,IPR,k,j,i)=w.e; prim(m,IPP,k,j,i)=w.pp;
      });
    };
    decode();
    const auto before = copy_to_host(prim);
    peos->Collisions(prim, bcc, cons, dt, 0,ncase-1, 0,0, 0,0);
    const auto corrected = copy_to_host(cons);
    const auto cached = copy_to_host(prim);
    decode();
    const auto fresh = copy_to_host(prim);
    for (int m=0; m<nmb; ++m) for (int c=0; c<ncase; ++c) {
      const Real delta = before(m,IPP,0,0,c)-before(m,IPR,0,0,c);
      const Real target = c == 0 ? -1.0 : c == 1 ? 1.0 : c == 2 ? -0.7 : 0.5;
      const Real rate = c < 2 ? 1.0e10 : 64.0;
      RequireClose("production two-tier map", fresh(m,IPP,0,0,c)-fresh(m,IPR,0,0,c),
                   target+(delta-target)/(1.0+rate*dt));
      RequireClose("production collision thermal energy", 0.5*fresh(m,IPR,0,0,c)
                   +fresh(m,IPP,0,0,c), 0.5*before(m,IPR,0,0,c)+before(m,IPP,0,0,c));
      for (int v=IDN; v <= (passive ? IM3 : IEN); ++v) {
        Require("collision preserves primary conserved fields",
                initial(m,v,0,0,c) == corrected(m,v,0,0,c));
      }
      for (int v : {IPR, IPP}) {
        Require("collision pressure cache matches independent decoder",
                cached(m,v,0,0,c) == fresh(m,v,0,0,c));
      }
    }
    const std::size_t bytes = corrected.size()*sizeof(Real);
    {
      std::ofstream file(filename, std::ios::binary);
      file.write(reinterpret_cast<const char*>(corrected.data()), bytes);
      Require("write finite backup state", file.good());
    }
    auto restored = Kokkos::create_mirror(cons);
    {
      std::ifstream file(filename, std::ios::binary);
      file.read(reinterpret_cast<char*>(restored.data()), bytes);
      Require("read finite backup state", file.good() && file.peek() == EOF);
    }
    Kokkos::deep_copy(cons, restored);
    decode();
    const auto reloaded = copy_to_host(prim);
    Require("persistent collision reload preserves pressure bits",
            std::memcmp(fresh.data(), reloaded.data(), fresh.size()*sizeof(Real)) == 0);
    // Zero-time calls and LF storage conversions must not add a hard projection.
    peos->Collisions(prim, bcc, cons, 0.0, 0,ncase-1, 0,0, 0,0);
    peos->CGLAnisotropyToMagneticMoment(cons, bcc, 0,ncase-1, 0,0, 0,0);
    peos->CGLRefreshPrimFromMagneticMoment(cons, bcc, prim, 0,ncase-1, 0,0, 0,0);
    peos->CGLMagneticMomentToAnisotropy(cons, bcc, prim, 0,ncase-1, 0,0, 0,0);
    decode();
    const auto converted = copy_to_host(prim);
    for (int m=0; m<nmb; ++m) {
      Require("finite firehose residual survives", converted(m,IPP,0,0,0)
              -converted(m,IPR,0,0,0) < -1.0);
      Require("finite mirror residual survives", converted(m,IPP,0,0,1)
              -converted(m,IPR,0,0,1) > 1.0);
      for (int c=0; c<ncase; ++c) for (int v : {IPR, IPP}) {
        RequireClose("LF conversion retains finite relaxation", converted(m,v,0,0,c),
                     fresh(m,v,0,0,c));
      }
    }
  }
  peos->eos_data = saved_eos;
  std::cout << "CGL finite backup/reload checks passed: active and passive" << std::endl;
}

// A tiny anisotropy loss must increase the small pressure, even at a binade edge.
void CheckProductionTinyCollisionRates(MeshBlockPack *pack) {
  auto *peos = pack->pmhd->peos;
  Require("tiny-rate fixture has no LF policy override", pack->pmhd->pcgl_lf == nullptr);
  const EOS_Data saved_eos = peos->eos_data;
  const int nmb = pack->nmb_thispack;
  constexpr Real state[9] = {0.3, 0., 0., 0., 51539607552., 17.706158816790296,
                             0.00042286396669162041, 0., 0.};
  DvceArray5D<Real> cons("tiny_rate_cons", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> prim("tiny_rate_prim", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> decoded("tiny_rate_decoded", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> bcc("tiny_rate_bcc", nmb, 3, 1, 1, 1);
  auto initial = Kokkos::create_mirror(cons);
  auto initial_b = Kokkos::create_mirror(bcc);
  for (int m=0; m<nmb; ++m) {
    for (int v=0; v<6; ++v) initial(m,v,0,0,0) = state[v];
    for (int v=0; v<3; ++v) initial_b(m,v,0,0,0) = state[v+6];
  }
  auto copy_to_host = [](const DvceArray5D<Real> &view) {
    auto host = Kokkos::create_mirror(view);
    Kokkos::deep_copy(host, view);
    return host;
  };
  for (const bool soft : {false, true}) {
    EOS_Data eos = saved_eos;
    eos.backup_lim = false;
    eos.flim = false;
    eos.mlim = soft;
    eos.mirror_threshold = 1.0;
    eos.nu_coll = soft ? 0.0 : 0.5*std::numeric_limits<Real>::epsilon();
    eos.lim_coll = soft ? std::numeric_limits<Real>::epsilon() : 0.0;
    peos->eos_data = eos;
    Kokkos::deep_copy(cons, initial);
    Kokkos::deep_copy(bcc, initial_b);
    auto decode = [&](DvceArray5D<Real> dest) {
      Kokkos::parallel_for("tiny_rate_independent_c2p", nmb,
      KOKKOS_LAMBDA(const int m) {
        MHDCons1D u{cons(m,IDN,0,0,0), cons(m,IM1,0,0,0), cons(m,IM2,0,0,0),
                    cons(m,IM3,0,0,0), cons(m,IEN,0,0,0), cons(m,IAN,0,0,0),
                    bcc(m,IBX,0,0,0), bcc(m,IBY,0,0,0), bcc(m,IBZ,0,0,0)};
        HydPrim1D w{};
        bool df=false, ef=false, tf=false, bf=false;
        SingleC2P_CGLMHD(u, eos, w, df, ef, tf, bf);
        if (df || ef || tf || bf || u.d != cons(m,IDN,0,0,0) ||
            u.e != cons(m,IEN,0,0,0) || u.mu != cons(m,IAN,0,0,0)) {
          Kokkos::abort("Tiny-rate stored state requires C2P repair");
        }
        dest(m,IDN,0,0,0)=w.d; dest(m,IVX,0,0,0)=w.vx;
        dest(m,IVY,0,0,0)=w.vy; dest(m,IVZ,0,0,0)=w.vz;
        dest(m,IPR,0,0,0)=w.e; dest(m,IPP,0,0,0)=w.pp;
      });
    };
    decode(prim);
    const auto before = copy_to_host(prim);
    // Put the canonical input exactly on its floor to expose any pressure decrease.
    eos.pfloor = before(0,IPR,0,0,0);
    peos->eos_data = eos;
    peos->Collisions(prim, bcc, cons, 1.0, 0, 0, 0, 0, 0, 0);
    decode(decoded);
    const auto after = copy_to_host(cons);
    const auto cached = copy_to_host(prim);
    const auto fresh = copy_to_host(decoded);
    const auto final_b = copy_to_host(bcc);
    for (int m=0; m<nmb; ++m) {
      for (int v=0; v<5; ++v) {
        Require("tiny rates conserve rho/momentum/total E bitwise",
                std::memcmp(&after(m,v,0,0,0), &initial(m,v,0,0,0), sizeof(Real)) == 0);
      }
      for (int v=0; v<6; ++v) {
        Require("tiny-rate primitive cache is canonical",
                cached(m,v,0,0,0) == fresh(m,v,0,0,0));
      }
      const long double ppar = before(m,IPR,0,0,0), pperp = before(m,IPP,0,0,0);
      const long double delta = pperp - ppar;
      const long double nudt = soft ? eos.lim_coll : eos.nu_coll;
      const long double threshold = 0.5L*state[6]*state[6];
      const long double loss = soft ? (delta-threshold)*nudt/(1.0L+nudt)
                                    : delta*(-std::expm1(-nudt));
      Require("tiny rates increase the small pressure", fresh(m,IPR,0,0,0) > ppar);
      RequireRelativeClose("tiny-rate parallel pressure", fresh(m,IPR,0,0,0),
                           static_cast<Real>(ppar + (2.0L/3.0L)*loss), kTol);
      RequireRelativeClose("tiny-rate perpendicular pressure", fresh(m,IPP,0,0,0),
                           static_cast<Real>(pperp - loss/3.0L), kTol);
    }
    Require("tiny rates leave magnetic field bitwise unchanged",
            std::memcmp(initial_b.data(), final_b.data(), final_b.size()*sizeof(Real)) == 0);
  }
  peos->eos_data = saved_eos;
  std::cout << "CGL production tiny-rate cancellation checks passed: 2 states per block"
            << std::endl;
}

void CheckPassiveRatioRecovery() {
  EOS_Data eos{};
  eos.dfloor = 1.0e-12;
  eos.pfloor = 1.0e-12;
  eos.bfloor = 1.0e-10;
  eos.sigma_max = 1.0e30;
  const Real large_cap = kSinglePrecision ? Real(1.0e3) : Real(1.0e10);
  for (const Real cap : {Real(1.0), Real(100.0), large_cap}) {
    eos.passive_max_pressure_ratio = cap;
    for (const bool perpendicular_hot : {false, true}) {
      const Real pp = perpendicular_hot ? Real(1.0e-6) : Real(1.0e5);
      const Real pt = perpendicular_hot ? Real(1.0e5) : Real(1.0e-6);
      const Real energy = pt + 0.5*pp;
      MHDCons1D u{};
      u.d = 1.7;
      u.bx = 0.3;
      const auto q = cgl::PassiveEncode(u.d, pp, pt, u.bx);
      u.e = q.j;
      u.mu = q.a;
      HydPrim1D w{};
      for (int cycle = 0; cycle < 4; ++cycle) {
        bool d = false, f = false, b = false, ratio = false;
        cgl::PassiveC2P(u, eos, w, d, f, b, &ratio);
        Require("ratio recovery respects floors",
                w.e >= eos.pfloor && w.pp >= eos.pfloor);
        RequireClose("ratio recovery thermal energy", w.pp + 0.5*w.e, energy);
        Require("ordinary ratio event only on first recovery", ratio == (cycle == 0));
        const Real pressure_ratio = perpendicular_hot ? w.pp/w.e : w.e/w.pp;
        if (cap <= 100.0) {
          RequireClose("ratio recovery reaches requested boundary", pressure_ratio, cap);
        }
        const auto unchanged_e = u.e, unchanged_mu = u.mu;
        ratio = false;
        cgl::PassiveC2P(u, eos, w, d, f, b, &ratio);
        Require("ordinary ratio recovery is idempotent R=" + std::to_string(cap) +
                " perpendicular_hot=" + std::to_string(perpendicular_hot),
                !ratio && u.e == unchanged_e && u.mu == unchanged_mu);
        // Start both representations outside the cap, then cycle the repaired state.
        u.e = cycle == 0 ? energy : w.pp + 0.5*w.e;
        u.mu = (cycle == 0 ? pt : w.pp)/u.bx;
        const Real moment_energy = u.e;
        cgl::PassiveMomentC2P(u, eos, w, d, f, b, u.bx, false, &ratio);
        if (cycle == 0) Require("unbounded LF moment is repaired", ratio);
        Require("moment ratio repair retains U", u.e == moment_energy);
        Require("moment ratio recovery respects floors",
                w.e >= eos.pfloor && w.pp >= eos.pfloor);
        d = f = b = ratio = false;
        cgl::PassiveMomentC2P(u, eos, w, d, f, b, u.bx, false, &ratio);
        Require("moment ratio recovery is idempotent", !ratio && !f);
        const auto encoded = cgl::PassiveEncodeFloored(
            u.d, w.e, w.pp, u.bx, eos.pfloor, cap);
        u.e = encoded.j;
        u.mu = encoded.a;
      }
    }
  }
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
  CheckStephenLimiterMap();
  CheckPassiveRatioRecovery();

  std::cout << "CGL C2P pressure-floor checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglC2PPressureFloorChecks();
  if (!kSinglePrecision) {
    CheckProductionFiniteBackup(pmy_mesh_->pmb_pack, "collision_state.bin");
    CheckProductionTinyCollisionRates(pmy_mesh_->pmb_pack);
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
