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
  eos.sigma_max = std::numeric_limits<float>::max();
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

  std::cout << "CGL C2P pressure-floor checks passed" << std::endl;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglC2PPressureFloorChecks();
}
