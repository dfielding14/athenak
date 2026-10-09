//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_c2p_storage_test.cpp
//! \brief Check repaired CGL states stored by the production EOS device kernels.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

void Require(const std::string &label, const bool condition) {
  if (!condition) {
    std::cout << "CGL C2P storage test failed: " << label << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void RequireClose(const std::string &label, const Real got, const Real expected) {
  const Real tolerance = 128.0*std::numeric_limits<Real>::epsilon();
  Require(label, std::isfinite(got) &&
                 std::abs(got - expected) <= tolerance*std::max(Real(1.0),
                                                              std::abs(expected)));
}

auto Snapshot(const DvceArray5D<Real> &view) {
  // An independent allocation prevents Serial snapshots from aliasing live storage.
  auto host = Kokkos::create_mirror(view);
  Kokkos::deep_copy(host, view);
  return host;
}

void CheckDensityFloorStorage(MeshBlockPack *pack) {
  auto *peos = pack->pmhd->peos;
  const EOS_Data saved_eos = peos->eos_data;
  peos->eos_data.dfloor = 1.0;
  peos->eos_data.pfloor = 1.0;
  peos->eos_data.bfloor = 1.0e-10;
  peos->eos_data.sigma_max = std::numeric_limits<Real>::max();
  peos->eos_data.hardwall_lim = false;

  const int nmb = pack->nmb_thispack;
  DvceArray5D<Real> cons("storage_cons", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> prim("storage_prim", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> bcc("storage_bcc", nmb, 3, 1, 1, 1);
  DvceFaceFld4D<Real> b("storage_b", nmb, 1, 1, 1);
  Kokkos::deep_copy(b.x1f, 0.0);
  Kokkos::deep_copy(b.x2f, 0.0);
  Kokkos::deep_copy(b.x3f, 1.0);
  Kokkos::deep_copy(bcc, 0.0);
  Kokkos::deep_copy(Kokkos::subview(bcc, Kokkos::ALL, static_cast<int>(IBZ),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 1.0);

  const char *names[] = {"anisotropy", "magnetic moment", "moment refresh"};
  for (int mode = 0; mode < 3; ++mode) {
    const std::string label = std::string(names[mode]) + " density floor";
    auto initial = Kokkos::create_mirror(cons);
    Kokkos::deep_copy(initial, 0.0);
    // rho=1/2, mx=1/2, Bz=1, p_parallel=2, p_perp=4 give E=23/4.
    for (int m = 0; m < nmb; ++m) {
      initial(m,IDN,0,0,0) = 0.5;
      initial(m,IM1,0,0,0) = 0.5;
      initial(m,IEN,0,0,0) = 5.75;
      initial(m,IAN,0,0,0) = mode == 0 ? 0.5*std::log(0.5) : 4.0;
    }
    Kokkos::deep_copy(cons, initial);

    auto convert = [&]() {
      if (mode == 0) {
        peos->ConsToPrim(cons, b, prim, bcc, false, 0, 0, 0, 0, 0, 0);
      } else if (mode == 1) {
        peos->CGLMagneticMomentToPrim(cons, b, prim, bcc, 0, 0, 0, 0, 0, 0);
      } else {
        peos->CGLRefreshPrimFromMagneticMoment(cons, bcc, prim, 0, 0, 0, 0, 0, 0);
      }
    };
    convert();
    const auto corrected = Snapshot(cons);
    const auto recovered = Snapshot(prim);
    // At rho=1, fixed E and momentum give U=41/8; A repair preserves p_perp/p_parallel.
    const Real expected_parallel = mode == 0 ? 2.05 : 2.25;
    const Real expected_perp = mode == 0 ? 4.1 : 4.0;
    const Real expected_slot = mode == 0 ? std::log(2.0) : 4.0;
    for (int m = 0; m < nmb; ++m) {
      RequireClose(label + " stored density", corrected(m,IDN,0,0,0), 1.0);
      RequireClose(label + " stored extra slot", corrected(m,IAN,0,0,0), expected_slot);
      for (int v = IM1; v <= IEN; ++v) {
        Require(label + " unchanged momentum and energy",
                corrected(m,v,0,0,0) == initial(m,v,0,0,0));
      }
      RequireClose(label + " primitive density", recovered(m,IDN,0,0,0), 1.0);
      RequireClose(label + " velocity x", recovered(m,IVX,0,0,0), 0.5);
      RequireClose(label + " velocity y", recovered(m,IVY,0,0,0), 0.0);
      RequireClose(label + " velocity z", recovered(m,IVZ,0,0,0), 0.0);
      RequireClose(label + " parallel pressure", recovered(m,IPR,0,0,0),
                   expected_parallel);
      RequireClose(label + " perpendicular pressure", recovered(m,IPP,0,0,0),
                   expected_perp);
      RequireClose(label + " energy identity", corrected(m,IEN,0,0,0),
                   0.5*recovered(m,IPR,0,0,0) + recovered(m,IPP,0,0,0) + 0.625);
    }

    convert();
    const auto repeated_cons = Snapshot(cons);
    const auto repeated_prim = Snapshot(prim);
    Require(label + " exact stored fixed point",
            std::memcmp(corrected.data(), repeated_cons.data(),
                        corrected.size()*sizeof(Real)) == 0);
    Require(label + " exact primitive fixed point",
            std::memcmp(recovered.data(), repeated_prim.data(),
                        recovered.size()*sizeof(Real)) == 0);
  }
  peos->eos_data = saved_eos;
  std::cout << "CGL C2P density-floor storage checks passed" << std::endl;
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;
  auto *pack = pmy_mesh_->pmb_pack;
  auto *pmhd = pack->pmhd;
  Require("active CGL without scalars", pmhd != nullptr &&
      pmhd->peos->eos_data.is_cgl && !pmhd->peos->eos_data.passive &&
      pmhd->nscalars == 0);
  CheckDensityFloorStorage(pack);

  // Normal zero-step startup still needs a finite live state after the isolated fixture.
  Kokkos::deep_copy(pmhd->u0, 0.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IDN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 1.0);
  Kokkos::deep_copy(Kokkos::subview(pmhd->u0, Kokkos::ALL, static_cast<int>(IEN),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 8.0);
  Kokkos::deep_copy(pmhd->b0.x1f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x2f, 0.0);
  Kokkos::deep_copy(pmhd->b0.x3f, 1.0);
}
