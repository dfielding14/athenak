//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_amr_projection_test.cpp
//! \brief Unit checks for CGL AMR projection and nonfinite magnetic-field policy.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/cgl_amr_projection.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL AMR projection test failed for " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

EOS_Data TestEOS() {
  EOS_Data eos{};
  eos.is_ideal = true;
  eos.is_cgl = true;
  eos.dfloor = 1.0e-12;
  eos.pfloor = 1.0e-12;
  eos.bfloor = 1.0e-10;
  eos.sigma_max = std::numeric_limits<Real>::max();
  eos.firehose_threshold = 2.0;
  eos.mirror_threshold = 1.0;
  eos.firehose_backup_factor = 1.0;
  eos.mirror_backup_factor = 2.0;
  eos.limiter_backup_nu = 1.0e10;
  return eos;
}

cgl::amr::ProjectionReport Project(const Real rho, const Real vx, const Real bx,
                                   const Real by, const Real bz, const Real U,
                                   const Real delta, MHDPrim1D &w, HydCons1D &u) {
  return cgl::amr::ProjectUDeltaToCGL(
      rho, vx, 0.0, 0.0, bx, by, bz, U, delta, TestEOS(),
      cgl::amr::anisotropy, w, u);
}

void CheckFiniteProjection() {
  MHDPrim1D w{};
  HydCons1D u{};
  const auto clean = Project(1.0, 0.1, 0.7, 0.2, -0.15, 1.5, 0.1, w, u);
  Require("clean projection report", clean.repairs == cgl::amr::kNone);
  Require("clean primitive state", cgl::amr::Finite3(w.d, w.e, w.pp));
  Require("clean conserved state", cgl::amr::Finite3(u.d, u.e, u.mu));

  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const Real infinity = std::numeric_limits<Real>::infinity();
  const auto repaired = Project(nan, infinity, 0.7, 0.2, -0.15, nan, nan, w, u);
  Require("nonfinite thermodynamics reported",
          (repaired.repairs & cgl::amr::kNonfiniteThermo) != 0u);
  Require("density repair reported",
          (repaired.repairs & cgl::amr::kDensityFloor) != 0u);
  Require("energy repair reported",
          (repaired.repairs & cgl::amr::kInternalEnergyFloor) != 0u);
  Require("repaired primitive state", cgl::amr::Finite3(w.d, w.e, w.pp));
  Require("repaired conserved state", cgl::amr::Finite3(u.d, u.e, u.mu));

  MHDCons1D conserved{};
  conserved.d = -1.0;
  conserved.e = 1.5 + 0.5*(SQR(0.7) + SQR(0.2) + SQR(-0.15));
  conserved.mu = 1.0;
  conserved.bx = 0.7;
  conserved.by = 0.2;
  conserved.bz = -0.15;
  const auto subfloor_density = cgl::amr::ProjectConservedToCGL(
      conserved, TestEOS(), cgl::amr::magnetic_moment, w, u);
  Require("conserved density floor reported",
          (subfloor_density.repairs & cgl::amr::kDensityFloor) != 0u);
  Require("conserved density projection finite",
          cgl::amr::Finite3(u.d, u.e, u.mu));

  conserved.d = nan;
  const auto nonfinite_density = cgl::amr::ProjectConservedToCGL(
      conserved, TestEOS(), cgl::amr::magnetic_moment, w, u);
  Require("conserved nonfinite density reported",
          (nonfinite_density.repairs & cgl::amr::kNonfiniteThermo) != 0u);
  Require("conserved nonfinite density floor reported",
          (nonfinite_density.repairs & cgl::amr::kDensityFloor) != 0u);
  Require("conserved nonfinite density projection finite",
          cgl::amr::Finite3(u.d, u.e, u.mu));
  std::cout << "CGL AMR finite projection checks passed" << std::endl;
}

void CheckStoredProjection(const MHDPrim1D &w, const HydCons1D &u,
                           const cgl::amr::ProjectionReport &report,
                           const EOS_Data &eos, const cgl::amr::SlotRepresentation slot,
                           const Real expected_delta, const Real expected_internal) {
  MHDCons1D stored{u.d, u.mx, u.my, u.mz, u.e, u.mu, w.bx, w.by, w.bz};
  const MHDCons1D original = stored;
  const Real bsqr = SQR(w.bx) + SQR(w.by) + SQR(w.bz);
  const Real kinetic = 0.5*(1.0/u.d)*(SQR(u.mx) + SQR(u.my) + SQR(u.mz));
  const Real internal = u.e - kinetic - 0.5*bsqr;
  const Real tolerance = 32.0*std::numeric_limits<Real>::epsilon()*
                         fmax(1.0, fmax(fabs(expected_internal), bsqr));
  for (int repeat=0; repeat<2; ++repeat) {
    HydPrim1D decoded{};
    bool df=false, ef=false, tf=false, bf=false;
    if (slot == cgl::amr::anisotropy) {
      SingleC2P_CGLMHD(stored, eos, decoded, df, ef, tf, bf);
    } else {
      SingleC2P_CGLMHDFromMagneticMoment(stored, eos, decoded, df, ef, tf, bf);
    }
    Require("projected state needs no C2P repair", !df && !ef && !tf && !bf);
    Require("projected conserved state is a C2P fixed point",
            stored.d == original.d && stored.mx == original.mx &&
            stored.my == original.my && stored.mz == original.mz &&
            stored.e == original.e && stored.mu == original.mu &&
            stored.bx == original.bx && stored.by == original.by && stored.bz == original.bz);
    Require("projected pressure cache matches stored state", decoded.e == w.e &&
            decoded.pp == w.pp && report.p_parallel == w.e && report.p_perp == w.pp);
    const Real delta = decoded.pp-decoded.e;
    // With the nondefault firehose coefficient 1, backup raises its wall to -B^2/2.
    Require("independent projected wall bounds", delta >= -bsqr &&
            (!eos.backup_lim || (delta >= -0.5*bsqr && delta <= bsqr)));
    Require("projected wall retains expected anisotropy",
            std::abs(delta-expected_delta) <= tolerance);
    Require("stored energy recovers projected thermal energy",
            std::abs(0.5*decoded.e+decoded.pp-internal) <= tolerance &&
            std::abs(internal-expected_internal) <= tolerance);
    MHDPrim1D wall = w;
    SingleCollWalls_CGLMHD(wall, eos, eos.backup_lim);
    Require("projected pressures need no second wall", wall.e == w.e && wall.pp == w.pp);
  }
}

void CheckCanonicalWallProjection() {
  EOS_Data eos = TestEOS();
  eos.firehose_threshold = 1.0;
  eos.mlim = eos.flim = true;
  // rho, velocity, B, U, Delta: the final case exercises stored-momentum rounding.
  const Real cases[][9] = {
    {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, -2.0},
    {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, -1.0},
    {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, 2.0},
    {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, -0.75},
    {Real(1.3), Real(0.7), Real(-0.2), Real(0.9), Real(0.7), Real(0.2), Real(0.4),
     4.0, -2.0}
  };
  for (const bool backup : {false, true}) {
    eos.backup_lim = backup;
    for (const auto &q : cases) {
      MHDPrim1D w{};
      HydCons1D u{};
      const Real bsqr = SQR(q[4]) + SQR(q[5]) + SQR(q[6]);
      Real expected_delta = fmax(q[8], backup ? -0.5*bsqr : -bsqr);
      if (backup) expected_delta = fmin(expected_delta, bsqr);
      const auto report = cgl::amr::ProjectUDeltaToCGL(
          q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], eos,
          cgl::amr::anisotropy, w, u);
      Require("wall projection preserves density and momentum",
              u.d == q[0] && u.mx == q[0]*q[1] && u.my == q[0]*q[2] && u.mz == q[0]*q[3]);
      Require("wall projection preserves magnetic field",
              w.bx == q[4] && w.by == q[5] && w.bz == q[6]);
      CheckStoredProjection(w, u, report, eos, cgl::amr::anisotropy,
                            expected_delta, q[7]);
    }
  }

  // A dyadic interior state retains magnetic-moment storage during an LF sweep.
  eos.backup_lim = false;
  MHDPrim1D w{};
  HydCons1D u{};
  const auto moment = cgl::amr::ProjectUDeltaToCGL(
      1.0, 0.5, -0.25, 0.125, 1.0, 0.0, 0.0, 3.0, 0.0, eos,
      cgl::amr::magnetic_moment, w, u);
  Require("AMR retains moment representation", u.mu == 2.0 && u.mu == w.pp/w.bx);
  CheckStoredProjection(w, u, moment, eos, cgl::amr::magnetic_moment, 0.0, 3.0);
  const MHDCons1D stored{u.d, u.mx, u.my, u.mz, u.e, u.mu, w.bx, w.by, w.bz};
  MHDPrim1D repeated_w{};
  HydCons1D repeated_u{};
  cgl::amr::ProjectConservedToCGL(stored, eos, cgl::amr::magnetic_moment,
                                 repeated_w, repeated_u);
  Require("dyadic moment projection fixed point",
          repeated_u.d == u.d && repeated_u.mx == u.mx && repeated_u.my == u.my &&
          repeated_u.mz == u.mz && repeated_u.e == u.e && repeated_u.mu == u.mu &&
          repeated_w.e == w.e && repeated_w.pp == w.pp);
  std::cout << "CGL AMR canonical wall encoding checks passed" << std::endl;
}

void CheckRepresentableProjectionFloor() {
  // vx, Bx, floor, U, Delta: kinetic loss, magnetic loss, moderate-kinetic weak B.
  const Real weak_floor = Real(0.75) + Real(3.0)*std::numeric_limits<Real>::epsilon();
  const Real cases[][5] = {
    {1.0e12, 1.0, 1.0, 3.0, -1.0},
    {0.0, 1.0, Real(1.0e-20), Real(1.5e-20), 0.0},
    {4.0, 0.0, weak_floor, Real(1.5)*weak_floor, 0.0}
  };
  const char *labels[] = {"kinetic cancellation", "magnetic cancellation",
                          "weak-field moderate kinetic rounding"};
  for (int c=0; c<3; ++c) {
    const auto &q = cases[c];
    const std::string label = labels[c];
    EOS_Data eos = TestEOS();
    eos.pfloor = q[2];
    MHDPrim1D w{};
    HydCons1D u{};
    const auto report = cgl::amr::ProjectUDeltaToCGL(
        1.0, q[0], 0.0, 0.0, q[1], 0.0, 0.0, q[3], q[4], eos,
        cgl::amr::anisotropy, w, u);
    Require(label + " reports energy repair",
            (report.repairs & cgl::amr::kInternalEnergyFloor) != 0u);
    Require(label + " reports pressure repair",
            (report.repairs & (cgl::amr::kParallelPressureFloor |
                               cgl::amr::kPerpPressureFloor)) != 0u);
    Require("AMR energy rounding preserves primary fields",
            u.d == 1.0 && u.mx == q[0] && u.my == 0.0 && u.mz == 0.0 &&
            w.bx == q[1] && w.by == 0.0 && w.bz == 0.0);
    Require("AMR rounded floor is finite and wall-admissible",
            std::isfinite(u.e) && std::isfinite(u.mu) &&
            w.e >= eos.pfloor && w.pp >= eos.pfloor && w.pp-w.e >= -SQR(q[1]));
    Require("AMR report returns encoded pressures", report.p_parallel == w.e &&
            report.p_perp == w.pp && report.U == w.pp+0.5*w.e);
    MHDCons1D stored{u.d, u.mx, u.my, u.mz, u.e, u.mu, w.bx, w.by, w.bz};
    for (int repeat=0; repeat<2; ++repeat) {
      HydPrim1D decoded{};
      bool df=false, ef=false, tf=false, bf=false;
      SingleC2P_CGLMHD(stored, eos, decoded, df, ef, tf, bf);
      Require("AMR rounded floor needs no second repair", !df && !ef && !tf);
      Require("AMR rounded floor is stored consistently",
              stored.e == u.e && stored.mu == u.mu && decoded.e == w.e && decoded.pp == w.pp);
    }
  }
  std::cout << "CGL AMR representable pressure-floor checks passed" << std::endl;
}

void CheckRejectedMagneticField(const std::string &mode) {
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const Real infinity = std::numeric_limits<Real>::infinity();
  const Real maximum = std::numeric_limits<Real>::max();
  Real bx = 0.7;
  Real by = 0.2;
  Real bz = -0.15;
  if (mode == "nan-component") {
    bx = nan;
  } else if (mode == "infinite-component") {
    by = infinity;
  } else if (mode == "overflowed-magnitude") {
    bz = maximum;
  } else {
    Require("recognized rejection mode", false);
  }

  MHDPrim1D w{};
  HydCons1D u{};
  (void) Project(1.0, 0.1, bx, by, bz, 1.5, 0.1, w, u);
  std::cout << "CGL AMR projection accepted " << mode << std::endl;
}

} // namespace

void RunCglAMRProjectionStandaloneChecks(const char *mode) {
  Kokkos::initialize();
  {
    const std::string selected(mode);
    if (selected == "finite") {
      CheckFiniteProjection();
      CheckCanonicalWallProjection();
      CheckRepresentableProjectionFloor();
    } else {
      CheckRejectedMagneticField(selected);
    }
  }
  Kokkos::finalize();
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  CheckFiniteProjection();
  CheckCanonicalWallProjection();
  CheckRepresentableProjectionFloor();
}
