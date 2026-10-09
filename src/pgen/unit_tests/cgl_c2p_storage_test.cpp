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

void CheckFloorStorage(MeshBlockPack *pack, const std::string &case_name,
                       const Real density, const Real initial_parallel,
                       const Real initial_perp, const Real density_floor,
                       const Real sigma_max, const Real bmag = 1.0,
                       const bool unresolved = false) {
  auto *peos = pack->pmhd->peos;
  const EOS_Data saved_eos = peos->eos_data;
  peos->eos_data.dfloor = density_floor;
  peos->eos_data.pfloor = 1.0;
  peos->eos_data.bfloor = 1.0e-10;
  peos->eos_data.sigma_max = sigma_max;
  peos->eos_data.hardwall_lim = false;

  const int nmb = pack->nmb_thispack;
  DvceArray5D<Real> cons("storage_cons", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> prim("storage_prim", nmb, 6, 1, 1, 1);
  DvceArray5D<Real> bcc("storage_bcc", nmb, 3, 1, 1, 1);
  DvceFaceFld4D<Real> b("storage_b", nmb, 1, 1, 1);
  Kokkos::deep_copy(b.x1f, 0.0);
  Kokkos::deep_copy(b.x2f, 0.0);
  Kokkos::deep_copy(b.x3f, bmag);
  Kokkos::deep_copy(bcc, 0.0);
  Kokkos::deep_copy(Kokkos::subview(bcc, Kokkos::ALL, static_cast<int>(IBZ),
      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), bmag);

  const Real expected_density =
      std::max(density, std::max(density_floor, SQR(bmag)/sigma_max));
  const Real momentum = unresolved ? std::ldexp(Real(1.0), 40) : 0.5;
  const Real initial_kinetic = 0.5*(1.0/density)*SQR(momentum);
  const Real magnetic = 0.5*SQR(bmag);
  const Real initial_energy = unresolved
      ? std::nextafter(initial_kinetic + magnetic, Real(0.0))
      : 0.5*initial_parallel + initial_perp + initial_kinetic + magnetic;
  const Real expected_kinetic = 0.5*(1.0/expected_density)*SQR(momentum);
  const char *names[] = {"anisotropy", "magnetic moment", "moment refresh"};
  for (int mode = 0; mode < 3; ++mode) {
    const std::string label = std::string(names[mode]) + " " + case_name;
    auto initial = Kokkos::create_mirror(cons);
    Kokkos::deep_copy(initial, 0.0);
    // Analytic energy and pressure ratios keep the storage oracle independent of C2P.
    for (int m = 0; m < nmb; ++m) {
      initial(m,IDN,0,0,0) = density;
      initial(m,IM1,0,0,0) = momentum;
      initial(m,IEN,0,0,0) = initial_energy;
      initial(m,IAN,0,0,0) = mode == 0
          ? density*std::log(initial_perp/initial_parallel*density*density/
                             (bmag*bmag*bmag)) : initial_perp/bmag;
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
    // Density repair releases kinetic energy at fixed total E and momentum.
    const Real thermal = initial_energy - expected_kinetic - magnetic;
    const Real ratio = initial_perp/initial_parallel;
    const Real unfloored_parallel =
        mode == 0 ? thermal/(0.5 + ratio) : 2.0*(thermal - initial_perp);
    const Real unfloored_perp = mode == 0 ? ratio*unfloored_parallel : initial_perp;
    const Real expected_parallel = std::max(Real(1.0), unfloored_parallel);
    const Real expected_perp = std::max(Real(1.0), unfloored_perp);
    const Real expected_slot = mode == 0
        ? expected_density*std::log(expected_perp/expected_parallel*
                                    expected_density*expected_density/(bmag*bmag*bmag))
        : expected_perp/bmag;
    const bool pressure_repair = unfloored_parallel < 1.0 || unfloored_perp < 1.0;
    const Real expected_energy = pressure_repair
        ? 0.5*expected_parallel + expected_perp + expected_kinetic + magnetic
        : initial_energy;
    auto &events = pack->pmesh->ecounter;
    const int density_events = events.neos_dfloor;
    const int pressure_events = events.neos_efloor;
    convert();
    Require(label + " expected density repairs",
            events.neos_dfloor == density_events + nmb*(expected_density > density));
    Require(label + " expected pressure repairs",
            events.neos_efloor == pressure_events + nmb*pressure_repair);
    const auto corrected = Snapshot(cons);
    const auto recovered = Snapshot(prim);
    for (int m = 0; m < nmb; ++m) {
      RequireClose(label + " stored density", corrected(m,IDN,0,0,0), expected_density);
      RequireClose(label + " stored extra slot", corrected(m,IAN,0,0,0), expected_slot);
      for (int v = IM1; v <= IM3; ++v) {
        Require(label + " unchanged momentum",
                corrected(m,v,0,0,0) == initial(m,v,0,0,0));
      }
      if (unresolved) {
        // The first floor E rounds to K; U=-1/2 forces the representability correction.
        Require(label + " representable stored energy", std::isfinite(corrected(m,IEN,0,0,0))
                && corrected(m,IEN,0,0,0) > initial_kinetic + magnetic);
      } else {
        RequireClose(label + " stored energy", corrected(m,IEN,0,0,0), expected_energy);
      }
      if (!pressure_repair) {
        Require(label + " unchanged energy", corrected(m,IEN,0,0,0) == initial_energy);
      }
      RequireClose(label + " primitive density", recovered(m,IDN,0,0,0), expected_density);
      RequireClose(label + " velocity x", recovered(m,IVX,0,0,0), momentum/expected_density);
      RequireClose(label + " velocity y", recovered(m,IVY,0,0,0), 0.0);
      RequireClose(label + " velocity z", recovered(m,IVZ,0,0,0), 0.0);
      Require(label + " finite pressures above floor",
              std::isfinite(recovered(m,IPR,0,0,0)) && recovered(m,IPR,0,0,0) >= 1.0 &&
              std::isfinite(recovered(m,IPP,0,0,0)) && recovered(m,IPP,0,0,0) >= 1.0);
      const Real stored_thermal = corrected(m,IEN,0,0,0) - expected_kinetic - magnetic;
      if (!unresolved) {
        // A representable E can shift pressures while preserving the floored ratio or moment.
        const Real floored_ratio = expected_perp/expected_parallel;
        const Real canonical_parallel = mode == 0
            ? stored_thermal/(0.5 + floored_ratio) : 2.0*(stored_thermal - expected_perp);
        const Real canonical_perp = mode == 0
            ? floored_ratio*canonical_parallel : expected_perp;
        RequireClose(label + " parallel pressure", recovered(m,IPR,0,0,0),
                     canonical_parallel);
        RequireClose(label + " perpendicular pressure", recovered(m,IPP,0,0,0),
                     canonical_perp);
      }
      const Real rounded_moment = Real(1.0)/bmag;
      if (mode != 0 && expected_perp == 1.0 && rounded_moment*bmag < 1.0) {
        Require(label + " upward moment encoding", corrected(m,IAN,0,0,0) > rounded_moment);
      }
      RequireClose(label + " internal energy identity", stored_thermal,
                   0.5*recovered(m,IPR,0,0,0) + recovered(m,IPP,0,0,0));
    }

    convert();
    Require(label + " no repeated density repair",
            events.neos_dfloor == density_events + nmb*(expected_density > density));
    Require(label + " no repeated pressure repair",
            events.neos_efloor == pressure_events + nmb*pressure_repair);
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
  std::cout << "CGL C2P " << case_name << " storage checks passed" << std::endl;
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
  const Real max_sigma = std::numeric_limits<Real>::max();
  CheckFloorStorage(pack, "density-floor", 0.5, 2.0, 4.0, 1.0, max_sigma);
  CheckFloorStorage(pack, "magnetization-floor", 0.5, 2.0, 4.0, 0.25, 1.0);
  CheckFloorStorage(pack, "parallel-pressure-floor", 2.0, 0.25, 2.0, 0.25, max_sigma);
  CheckFloorStorage(pack, "perpendicular-pressure-floor", 2.0, 2.0, 0.25, 0.25, max_sigma);
  CheckFloorStorage(pack, "both-pressure-floor", 2.0, 0.25, 0.5, 0.25, max_sigma);
  const Real rounding_b = sizeof(Real) == sizeof(float) ? 41.0 : 49.0;
  Require("moment division fixture rounds pressure down", (Real(1.0)/rounding_b)*rounding_b < 1.0);
  CheckFloorStorage(pack, "moment-division-rounding", 1.0, 5.5, 0.25, 0.25,
                    max_sigma, rounding_b);
  CheckFloorStorage(pack, "unresolved-internal-energy", 2.0, 1.0, 1.0, 0.25,
                    max_sigma, 1.0, true);

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
