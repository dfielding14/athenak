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
  eos.dfloor = 1.0e-12;
  eos.pfloor = 1.0e-12;
  eos.bfloor = 1.0e-10;
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
  std::cout << "CGL AMR finite projection checks passed" << std::endl;
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
}
