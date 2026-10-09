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
#include "eos/ideal_c2p_mhd.hpp"
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
  const Real b_eff = std::max(bmag, peos->eos_data.bfloor);
  const bool low_b = bmag <= peos->eos_data.bfloor;

  const int nmb = pack->nmb_thispack;
  const int scalar = pack->pmhd->nmhd;
  DvceArray5D<Real> cons("storage_cons", nmb, scalar+1, 1, 1, 1);
  DvceArray5D<Real> prim("storage_prim", nmb, scalar+1, 1, 1, 1);
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
                             (b_eff*b_eff*b_eff)) : initial_perp/b_eff;
      initial(m,scalar,0,0,0) = 0.25;
    }
    Kokkos::deep_copy(cons, initial);
    Kokkos::deep_copy(prim, 0.0);
    Kokkos::deep_copy(Kokkos::subview(prim, Kokkos::ALL, scalar,
        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), 0.25/density);

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
    const Real unfloored_parallel = low_b ? 2.0*thermal/3.0
        : (mode == 0 ? thermal/(0.5 + ratio) : 2.0*(thermal - initial_perp));
    const Real unfloored_perp = low_b ? unfloored_parallel
        : (mode == 0 ? ratio*unfloored_parallel : initial_perp);
    const Real expected_parallel = std::max(Real(1.0), unfloored_parallel);
    const Real expected_perp = std::max(Real(1.0), unfloored_perp);
    const Real expected_slot = mode == 0
        ? expected_density*std::log(expected_perp/expected_parallel*
                                    expected_density*expected_density/(b_eff*b_eff*b_eff))
        : expected_perp/b_eff;
    const bool pressure_repair = unfloored_parallel < 1.0 || unfloored_perp < 1.0;
    const Real expected_energy = pressure_repair
        ? 0.5*expected_parallel + expected_perp + expected_kinetic + magnetic
        : initial_energy;
    auto &events = pack->pmesh->ecounter;
    const int density_events = events.neos_dfloor;
    const int pressure_events = events.neos_efloor;
    if (mode == 0) {
      const auto saved_fofc = pack->pmhd->fofc;
      pack->pmhd->fofc = DvceArray4D<bool>("storage_fofc", nmb, 1, 1, 1);
      Kokkos::deep_copy(pack->pmhd->fofc, false);
      const auto initial_prim = Snapshot(prim);
      const auto initial_bcc = Snapshot(bcc);
      const int fofc_events = events.nfofc;
      peos->ConsToPrim(cons, b, prim, bcc, true, 0, 0, 0, 0, 0, 0);
      const auto tested_cons = Snapshot(cons);
      const auto tested_prim = Snapshot(prim);
      const auto tested_bcc = Snapshot(bcc);
      const auto flags = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pack->pmhd->fofc);
      const bool expect_flag = expected_density > density || pressure_repair;
      Require(label + " test-floor conserved state unchanged",
              std::memcmp(initial.data(), tested_cons.data(), initial.size()*sizeof(Real)) == 0);
      Require(label + " test-floor primitives unchanged",
              std::memcmp(initial_prim.data(), tested_prim.data(),
                          initial_prim.size()*sizeof(Real)) == 0);
      Require(label + " test-floor magnetic field unchanged",
              std::memcmp(initial_bcc.data(), tested_bcc.data(),
                          initial_bcc.size()*sizeof(Real)) == 0);
      Require(label + " test-floor counters",
              events.nfofc == fofc_events + nmb*expect_flag &&
              events.neos_dfloor == density_events && events.neos_efloor == pressure_events);
      for (int m = 0; m < nmb; ++m) {
        Require(label + " test-floor cell flag", flags(m,0,0,0) == expect_flag);
      }
      pack->pmhd->fofc = saved_fofc;
    }
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
      if (low_b) {
        Require(label + " repaired low-B slot",
                corrected(m,IAN,0,0,0) != initial(m,IAN,0,0,0));
      }
      Require(label + " preserved scalar mass", corrected(m,scalar,0,0,0) == 0.25);
      RequireClose(label + " scalar fraction", recovered(m,scalar,0,0,0),
                   0.25/expected_density);
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
        const Real canonical_parallel = low_b ? 2.0*stored_thermal/3.0
            : (mode == 0 ? stored_thermal/(0.5 + floored_ratio)
                         : 2.0*(stored_thermal - expected_perp));
        const Real canonical_perp = low_b ? canonical_parallel
            : (mode == 0 ? floored_ratio*canonical_parallel : expected_perp);
        RequireClose(label + " parallel pressure", recovered(m,IPR,0,0,0),
                     canonical_parallel);
        RequireClose(label + " perpendicular pressure", recovered(m,IPP,0,0,0),
                     canonical_perp);
      }
      const Real rounded_moment = Real(1.0)/b_eff;
      if (!low_b && mode != 0 && expected_perp == 1.0 && rounded_moment*bmag < 1.0) {
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

void CheckExtremeRatiosOnDevice() {
  const int exponent = sizeof(Real) == sizeof(float) ? 80 : 600;
  const Real large = std::ldexp(Real(1.0), exponent);
  const Real small = std::ldexp(Real(1.0), -exponent);
  const Real below_floor = std::ldexp(Real(1.0), -exponent-20);
  EOS_Data base_eos{};
  base_eos.is_cgl = base_eos.is_ideal = true;
  base_eos.dfloor = base_eos.bfloor = 1.0e-12;
  base_eos.sigma_max = std::numeric_limits<Real>::max();
  DvceArray2D<Real> pressures("extreme_ratio_pressures", 4, 2);
  DvceArray2D<int> checks("extreme_ratio_checks", 4, 4);
  Kokkos::parallel_for("cgl_extreme_ratio_checks",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, 4), KOKKOS_LAMBDA(const int c) {
    const bool floor_case = c >= 2;
    const Real initial_small = floor_case ? below_floor : small;
    const Real parallel = c%2 == 0 ? initial_small : large;
    const Real perp = c%2 == 0 ? large : initial_small;
    EOS_Data eos = base_eos;
    eos.pfloor = floor_case ? small : small/16.0;
    MHDCons1D u{};
    u.d = u.bz = 1.0;
    u.e = 0.5*parallel + perp + 0.5;
    u.mu = log(perp) - log(parallel);
    const Real initial_energy = u.e, initial_anisotropy = u.mu;
    HydPrim1D w{};
    bool df=false, ef=false, tf=false, bf=false;
    SingleC2P_CGLMHD(u, eos, w, df, ef, tf, bf);
    pressures(c,0) = w.e;
    pressures(c,1) = w.pp;
    checks(c,0) = !df && ef == floor_case && !tf && !bf;
    checks(c,1) = floor_case || (u.e == initial_energy && u.mu == initial_anisotropy);
    checks(c,2) = Kokkos::isfinite(u.e) && Kokkos::isfinite(u.mu);
    const Real repaired_energy = u.e, repaired_anisotropy = u.mu;
    const Real repaired_parallel = w.e, repaired_perp = w.pp;
    df = ef = tf = bf = false;
    SingleC2P_CGLMHD(u, eos, w, df, ef, tf, bf);
    checks(c,3) = !df && !ef && !tf && !bf &&
        u.e == repaired_energy && u.mu == repaired_anisotropy &&
        w.e == repaired_parallel && w.pp == repaired_perp;
  });
  const auto host_p = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pressures);
  const auto host_checks = Kokkos::create_mirror_view_and_copy(HostMemSpace(), checks);
  const Real tolerance = 2048.0*std::numeric_limits<Real>::epsilon();
  for (int c = 0; c < 4; ++c) {
    const std::string label = "extreme ratio device case " + std::to_string(c);
    const Real floor = c >= 2 ? small : small/16.0;
    for (int p = 0; p < 2; ++p) {
      const Real expected = p == c%2 ? small : large;
      Require(label + " pressure bounds", std::isfinite(host_p(c,p)) && host_p(c,p) >= floor);
      Require(label + " analytic pressure", std::abs(host_p(c,p) - expected) <=
                                           tolerance*expected);
    }
    Require(label + " first floor flags", host_checks(c,0));
    Require(label + " healthy E/A unchanged", host_checks(c,1));
    Require(label + " finite stored E/A", host_checks(c,2));
    Require(label + " exact second conversion", host_checks(c,3));
  }
  std::cout << "CGL C2P extreme-ratio device checks passed" << std::endl;
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
  // Only the isolated views own a scalar; restore metadata before live mesh startup.
  const int saved_nscalars = pmhd->nscalars;
  pmhd->nscalars = 1;
  const Real max_sigma = std::numeric_limits<Real>::max();
  CheckFloorStorage(pack, "density-floor", 0.5, 2.0, 4.0, 1.0, max_sigma);
  CheckFloorStorage(pack, "magnetization-floor", 0.5, 2.0, 4.0, 0.25, 1.0);
  CheckFloorStorage(pack, "parallel-pressure-floor", 2.0, 0.25, 2.0, 0.25, max_sigma);
  CheckFloorStorage(pack, "perpendicular-pressure-floor", 2.0, 2.0, 0.25, 0.25, max_sigma);
  CheckFloorStorage(pack, "both-pressure-floor", 2.0, 0.25, 0.5, 0.25, max_sigma);
  CheckFloorStorage(pack, "low-magnetic-field", 2.0, 4.0, 1.0, 0.25, max_sigma, 0.0);
  const Real rounding_b = sizeof(Real) == sizeof(float) ? 41.0 : 49.0;
  Require("moment division fixture rounds pressure down", (Real(1.0)/rounding_b)*rounding_b < 1.0);
  CheckFloorStorage(pack, "moment-division-rounding", 1.0, 5.5, 0.25, 0.25,
                    max_sigma, rounding_b);
  CheckFloorStorage(pack, "unresolved-internal-energy", 2.0, 1.0, 1.0, 0.25,
                    max_sigma, 1.0, true);
  pmhd->nscalars = saved_nscalars;
  CheckExtremeRatiosOnDevice();

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
