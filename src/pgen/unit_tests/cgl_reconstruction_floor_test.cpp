//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_reconstruction_floor_test.cpp
//! \brief Exercise pressure floors through the production directional reconstruction.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"

namespace {

void CheckWrapper(const ReconstructionMethod method, const int direction,
                  const bool is_cgl, const bool apply_floors, const int orientation) {
  constexpr int nvars = 7;
  constexpr int width = 5;
  constexpr int center = 2;
  const Real value[nvars] = {0.5, -1.0, 2.0, -3.0, 1.0, 1.0, 0.125};
  const Real slope[nvars] = {0.125, 0.25, -0.125, 0.125, 0.25, -0.25, 0.03125};
  EOS_Data eos{};
  eos.is_ideal = true;
  eos.is_cgl = is_cgl;
  eos.gamma = 5.0/3.0;
  eos.dfloor = 0.5;
  eos.pfloor = 1.0;

  DvceArray5D<Real> q("reconstruction_stencil", 1, nvars, width, width, width);
  auto host_q = Kokkos::create_mirror(q);
  for (int n = 0; n < nvars; ++n) {
    for (int k = 0; k < width; ++k) {
      for (int j = 0; j < width; ++j) {
        for (int i = 0; i < width; ++i) {
          const int coordinate = direction == 0 ? i : (direction == 1 ? j : k);
          host_q(0,n,k,j,i) = value[n] + orientation*slope[n]*(coordinate - center);
        }
      }
    }
  }
  Kokkos::deep_copy(q, host_q);
  DvceArray2D<Real> faces("reconstruction_faces", 2, nvars);
  const size_t scratch_size = 2*ScrArray2D<Real>::shmem_size(nvars, width);
  par_for_outer("cgl_reconstruction_wrappers", DevExeSpace(), scratch_size, 0, 0, 0,
  KOKKOS_LAMBDA(TeamMember_t member, const int) {
    ScrArray2D<Real> ql(member.team_scratch(0), nvars, width);
    ScrArray2D<Real> qr(member.team_scratch(0), nvars, width);
    if (method == ReconstructionMethod::wenoz) {
      if (direction == 0) {
        WENOZX1(member, eos, apply_floors, 0, center, center, center, center, q, ql, qr);
      } else if (direction == 1) {
        WENOZX2(member, eos, apply_floors, 0, center, center, center, center, q, ql, qr);
      } else {
        WENOZX3(member, eos, apply_floors, 0, center, center, center, center, q, ql, qr);
      }
    } else {
      const bool extrema = method == ReconstructionMethod::ppmx;
      if (direction == 0) {
        PiecewiseParabolicX1(member, eos, extrema, apply_floors,
                             0, center, center, center, center, q, ql, qr);
      } else if (direction == 1) {
        PiecewiseParabolicX2(member, eos, extrema, apply_floors,
                             0, center, center, center, center, q, ql, qr);
      } else {
        PiecewiseParabolicX3(member, eos, extrema, apply_floors,
                             0, center, center, center, center, q, ql, qr);
      }
    }
    member.team_barrier();
    par_for_inner(member, 0, nvars-1, [&](const int n) {
      faces(0,n) = ql(n, direction == 0 ? center+1 : center);
      faces(1,n) = qr(n, center);
    });
  });

  const auto host_faces = Kokkos::create_mirror_view_and_copy(HostMemSpace(), faces);
  const char *method_name = method == ReconstructionMethod::wenoz ? "wenoz"
      : (method == ReconstructionMethod::ppmx ? "ppmx" : "ppm4");
  const Real tolerance = 64.0*std::numeric_limits<Real>::epsilon();
  for (int side = 0; side < 2; ++side) {
    for (int n = 0; n < nvars; ++n) {
      // Every method reconstructs an affine profile exactly before its optional floors.
      Real expected = value[n] + (side == 0 ? 0.5 : -0.5)*orientation*slope[n];
      if (apply_floors && method != ReconstructionMethod::ppm4) {
        if (n == IDN) expected = std::max(expected, Real(0.5));
        if (n == IEN) expected = std::max(expected, is_cgl ? Real(1.0) : Real(1.5));
        if (is_cgl && n == IPP) expected = std::max(expected, Real(1.0));
      }
      const Real got = host_faces(side,n);
      if (!std::isfinite(got) || std::abs(got - expected) >
          tolerance*std::max(Real(1.0), std::abs(expected))) {
        std::cout << "CGL reconstruction floor check failed: method=" << method_name
                  << " direction=" << direction << " cgl=" << is_cgl
                  << " floors=" << apply_floors << " orientation=" << orientation
                  << " side=" << side << " variable=" << n
                  << " got=" << got << " expected=" << expected << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
}

void RunCglReconstructionFloorChecks() {
  for (const auto method : {ReconstructionMethod::ppmx, ReconstructionMethod::wenoz,
                            ReconstructionMethod::ppm4}) {
    for (int direction = 0; direction < 3; ++direction) {
      for (const bool is_cgl : {false, true}) {
        for (const bool apply_floors : {false, true}) {
          for (const int orientation : {-1, 1}) {
            CheckWrapper(method, direction, is_cgl, apply_floors, orientation);
          }
        }
      }
    }
  }
  std::cout << "CGL directional reconstruction floor checks passed" << std::endl;
}

} // namespace

void RunCglReconstructionFloorStandaloneChecks() {
  Kokkos::initialize();
  RunCglReconstructionFloorChecks();
  Kokkos::finalize();
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  RunCglReconstructionFloorChecks();
  auto *pmhd = pmy_mesh_->pmb_pack->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl || pmhd->peos->eos_data.passive) {
    std::cout << "CGL reconstruction fixture requires active CGL MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // The isolated wrapper checks leave a uniform physical state for zero-step startup.
  Kokkos::deep_copy(pmhd->u0, 0.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IDN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 1.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IEN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 2.0);
  Kokkos::deep_copy(pmhd->b0.x1f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x2f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x3f, 1.0);
}
