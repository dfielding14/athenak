//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_lf_cooling.cpp
//! \brief CGL turbulence with history-controlled cooling at fixed pressure difference.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {

void Fail(const std::string &message) {
  std::cout << "### FATAL ERROR: cgl_lf_cooling: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

// SetReal truncates checkpoint parameters to six significant digits.
void SaveReal(ParameterInput *pin, const std::string &block,
              const std::string &name, const Real value) {
  std::ostringstream text;
  text << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
  pin->SetString(block, name, text.str());
}

struct Cooling {
  Mesh *mesh;
  ParameterInput *pin;
  UserHistoryFnPtr paper_history;
  std::string history_block;
  Real target, tau, alpha = 0.0;
  Real target_fraction = 1.0, sample_time = 0.0;
  Real history_dt = 0.0, last_time = -1.0;
  int history_dcycle = 0, last_cycle = -1;
  DvceArray1D<std::uint64_t> guard_hits;
  std::uint64_t last_guard_hits = 0, interval_guard_hits = 0;
  int guard_cycle = 0;
  Real interval_guard_fraction = 0.0;

  TaskStatus Source(Driver *pdrive, int stage) {
    if (alpha == 0.0) return TaskStatus::complete;
    auto *pmhd = mesh->pmb_pack->pmhd;
    pmhd->RequireCGLAnisotropyRepresentation("CGL cooling source");
    auto u = pmhd->u0;
    auto w = pmhd->w0;
    auto bcc = pmhd->bcc0;
    const auto eos = pmhd->peos->eos_data;
    constexpr Real cooling_max_ratio = 10.0;
    const Real age = fmax(mesh->time - sample_time, 0.0);
    const Real remaining = (1.0 - target_fraction)*std::exp(-age/tau);
    const Real excess_fraction = target_fraction == 0.0 ? 1.0 :
        remaining/(target_fraction + remaining);
    const Real fraction = excess_fraction*(-std::expm1(-mesh->dt/tau));
    const Real retained = target_fraction == 0.0 ? std::exp(-mesh->dt/tau) :
        (target_fraction + remaining*std::exp(-mesh->dt/tau))/
        (target_fraction + remaining);
    if (fraction == 0.0) return TaskStatus::complete;
    const auto &indcs = mesh->mb_indcs;
    const int ni = indcs.nx1 + 2*indcs.ng;
    const int nj = indcs.nx2 > 1 ? indcs.nx2 + 2*indcs.ng : 1;
    const int nk = indcs.nx3 > 1 ? indcs.nx3 + 2*indcs.ng : 1;
    const int is = indcs.is, ie = indcs.ie;
    const int js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    auto guard_hits_ = guard_hits;

    // Cooling periodic ghosts identically avoids another boundary exchange.
    par_for("cgl_cooling", DevExeSpace(), 0, mesh->pmb_pack->nmb_thispack - 1,
            0, nk - 1, 0, nj - 1, 0, ni - 1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      if (eos.passive_adiabatic) {
        const Real rho = u(m,IDN,k,j,i);
        const Real ek = (0.5/rho)*(SQR(u(m,IM1,k,j,i)) + SQR(u(m,IM2,k,j,i)) +
                                  SQR(u(m,IM3,k,j,i)));
        const Real em = 0.5*(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) +
                             SQR(bcc(m,IBZ,k,j,i)));
        const Real floor = fmax(eos.pfloor, fmax(rho*eos.tfloor,
            eos.sfloor*pow(rho, eos.gamma)))/(eos.gamma - 1.0);
        const Real old_energy = u(m,IEG,k,j,i);
        const Real thermal_minimum = fraction > 0.5 ?
            fmax(retained*w(m,IEG,k,j,i), floor) : floor;
        Real energy = fraction > 0.5 ? (thermal_minimum + ek) + em :
            fmax(old_energy - fraction*w(m,IEG,k,j,i), (floor + ek) + em);
        // Retain a representable thermal floor without adding energy.
        while (energy < old_energy && energy - ek - em < thermal_minimum) {
          energy = Kokkos::nextafter(energy, std::numeric_limits<Real>::infinity());
        }
        if (energy < old_energy) {
          u(m,IEG,k,j,i) = energy;
          w(m,IEG,k,j,i) = energy - ek - em;
        }
      }
      Real ppar = w(m,IPR,k,j,i);
      Real pperp = w(m,IPP,k,j,i);
      const Real delta = pperp - ppar;
      const Real piso = TWO_3RDS*pperp + ONE_3RD*ppar;
      // shortcut: cap cooling at a 10:1 pressure ratio; revisit if many cells hit it.
      const Real cooling_floor = fmax(eos.pfloor, fabs(delta)/(cooling_max_ratio - 1.0));
      const Real available = fmax(fmin(ppar, pperp) - cooling_floor, 0.0);
      const Real requested = fraction*piso;
      const Real floor_available = fmax(fmin(ppar, pperp) - eos.pfloor, 0.0);
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke &&
          fmin(requested, floor_available) > available) {
        Kokkos::atomic_increment(&guard_hits_(0));
      }
      const Real decrement = fmin(requested, available);
      if (decrement == 0.0) return;
      // Limit the common decrement, never floor the two pressures independently.
      if (delta >= 0.0) {
        ppar = fmax(fraction > 0.5 ? retained*piso - TWO_3RDS*delta :
                    ppar - decrement, cooling_floor);
        pperp = ppar + delta;
      } else {
        pperp = fmax(fraction > 0.5 ? retained*piso + ONE_3RD*delta :
                     pperp - decrement, cooling_floor);
        ppar = pperp - delta;
      }
      const Real bx = bcc(m,IBX,k,j,i);
      const Real by = bcc(m,IBY,k,j,i);
      const Real bz = bcc(m,IBZ,k,j,i);
      const Real bmag = fmax(sqrt(SQR(bx) + SQR(by) + SQR(bz)), eos.bfloor);
      if (eos.passive_adiabatic) {
        // Passive cooling cannot alter the independent MHD energy or its controller.
        const auto q = cgl::PassiveEncodeFloored(u(m,IDN,k,j,i), ppar, pperp,
                                                 bmag, eos.pfloor);
        u(m,IEN,k,j,i) = q.j;
        u(m,IAN,k,j,i) = q.a;
        cgl::PassiveDecode(u(m,IDN,k,j,i), q, bmag, ppar, pperp);
        w(m,IPR,k,j,i) = ppar;
        w(m,IPP,k,j,i) = pperp;
        return;
      }
      const Real old_energy = u(m,IEN,k,j,i);
      Real energy = old_energy - 1.5*decrement;
      if (fraction > 0.5) {
        const Real ek = (0.5/u(m,IDN,k,j,i))*(SQR(u(m,IM1,k,j,i)) +
            SQR(u(m,IM2,k,j,i)) + SQR(u(m,IM3,k,j,i)));
        const Real em = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
        const Real thermal = 0.5*ppar + pperp;
        energy = (thermal + ek) + em;
        // Near-complete cooling retains the small thermal remainder directly.
        while (energy < old_energy && energy - ek - em < thermal) {
          energy = Kokkos::nextafter(energy, std::numeric_limits<Real>::infinity());
        }
      }
      if (energy >= old_energy) return;
      Real anisotropy = CGLConservedAnisotropy(u(m,IDN,k,j,i), ppar, pperp, bmag);
      // Keep the pressure cache consistent with the stored energy and anisotropy.
      CGLRecoverPressuresFromTotalEnergyAndAnisotropy(
          u(m,IDN,k,j,i), u(m,IM1,k,j,i), u(m,IM2,k,j,i), u(m,IM3,k,j,i),
          energy, anisotropy, bx, by, bz, eos.bfloor, ppar, pperp);
      if (ppar < eos.pfloor || pperp < eos.pfloor) {
        // Tiny thermal remainders can be lost against kinetic/magnetic energy.
        // Reuse the EOS's upward rounding to keep the stored state above its floors.
        MHDCons1D candidate;
        candidate.d = u(m,IDN,k,j,i);
        candidate.mx = u(m,IM1,k,j,i);
        candidate.my = u(m,IM2,k,j,i);
        candidate.mz = u(m,IM3,k,j,i);
        candidate.e = energy;
        candidate.mu = anisotropy;
        candidate.bx = bx;
        candidate.by = by;
        candidate.bz = bz;
        HydPrim1D recovered;
        bool dfloor_used = false, efloor_used = false;
        bool tfloor_used = false, bfloor_used = false;
        SingleC2P_CGLMHD(candidate, eos, recovered, dfloor_used, efloor_used,
                         tfloor_used, bfloor_used);
        if (candidate.e >= old_energy) return;  // No change without actual cooling.
        energy = candidate.e;
        anisotropy = candidate.mu;
        ppar = recovered.e;
        pperp = recovered.pp;
      }
      u(m,IEN,k,j,i) = energy;
      u(m,IAN,k,j,i) = anisotropy;
      w(m,IPR,k,j,i) = ppar;
      w(m,IPP,k,j,i) = pperp;
    });
    // Guard the post-RK LF timestep when the passive cap is tighter than cooling's.
    if (eos.passive && eos.passive_max_pressure_ratio != 0.0 &&
        eos.passive_max_pressure_ratio < cooling_max_ratio && pmhd->has_cgl_lf_split) {
      pmhd->ConToPrim(pdrive, stage);
    }
    return TaskStatus::complete;
  }
};

Cooling cooling;

void CoolingHistory(HistoryData *pdata, Mesh *pm) {
  cooling.paper_history(pdata, pm);
  const bool due = cooling.last_cycle < 0 ||
      (cooling.history_dcycle > 0 && pm->ncycle%cooling.history_dcycle == 0) ||
      (cooling.history_dt > 0.0 &&
       static_cast<float>(pm->time) >=
       static_cast<float>(cooling.last_time + cooling.history_dt));
  // Finalize also calls history off cadence and can repeat the last sample.
  if (!pm->pmb_pack->frozen_mhd && due && pm->ncycle != cooling.last_cycle) {
    const auto eos = pm->pmb_pack->pmhd->peos->eos_data;
    const std::string thermal = eos.passive_adiabatic ? "therm_mhd" : "therm_cgl";
    const Real heat_capacity = eos.passive_adiabatic ? 1.0/(eos.gamma - 1.0) : 1.5;
    Real sums[2] = {};
    for (int n = 0; n < pdata->nhist; ++n) {
      if (pdata->label[n] == thermal) sums[0] = pdata->hdata[n];
      if (pdata->label[n] == "volume") sums[1] = pdata->hdata[n];
    }
#if MPI_PARALLEL_ENABLED
    // History's own later Reduce needs the original rank-local hdata untouched.
    MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (!std::isfinite(sums[0]) || sums[0] <= 0.0 ||
        !std::isfinite(sums[1]) || sums[1] <= 0.0) {
      Fail("history requires positive finite thermal energy and volume");
    }
    cooling.target_fraction = fmin(heat_capacity*cooling.target*sums[1]/sums[0], 1.0);
    cooling.sample_time = pm->time;
    cooling.alpha = (1.0 - cooling.target_fraction)/cooling.tau;
    cooling.last_time = cooling.last_cycle < 0 ? pm->time :
        cooling.last_time + cooling.history_dt;
    cooling.last_cycle = pm->ncycle;
    SaveReal(cooling.pin, "problem", "cooling_alpha", cooling.alpha);
    SaveReal(cooling.pin, "problem", "cooling_target_fraction", cooling.target_fraction);
    SaveReal(cooling.pin, "problem", "cooling_sample_time", cooling.sample_time);
    SaveReal(cooling.pin, "problem", "cooling_history_last_time", cooling.last_time);
    cooling.pin->SetInteger("problem", "cooling_history_last_cycle", cooling.last_cycle);
  }
  pdata->label[pdata->nhist] = "cool_alpha";
  pdata->hdata[pdata->nhist++] = global_variable::my_rank == 0 ? cooling.alpha : 0.0;
  if (pm->ncycle != cooling.guard_cycle) {
    std::uint64_t hits;
    Kokkos::deep_copy(hits, Kokkos::subview(cooling.guard_hits, 0));
    cooling.interval_guard_hits = hits - cooling.last_guard_hits;
    const auto &indcs = pm->mb_indcs;
    const Real updates = Real(pm->nmb_total)*indcs.nx1*indcs.nx2*indcs.nx3*
                         (pm->ncycle - cooling.guard_cycle);
    // Generic history sums local counts and fractions over all ranks.
    cooling.interval_guard_fraction = cooling.interval_guard_hits/updates;
    cooling.last_guard_hits = hits;
    cooling.guard_cycle = pm->ncycle;
  }
  pdata->label[pdata->nhist] = "cool_hits";
  pdata->hdata[pdata->nhist++] = cooling.interval_guard_hits;
  pdata->label[pdata->nhist] = "cool_frac";
  pdata->hdata[pdata->nhist++] = cooling.interval_guard_fraction;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \brief Reuse the paper setup, with a pgen-local split user source and history wrapper.

void ProblemGenerator::CGLLFCooling(ParameterInput *pin, const bool restart) {
  auto *pm = pmy_mesh_;
  auto *pmbp = pm->pmb_pack;
  auto *pmhd = pmbp->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl ||
      (pmhd->peos->eos_data.passive && !pmhd->peos->eos_data.passive_adiabatic) ||
      (pmbp->frozen_mhd && !restart) || pm->multilevel) {
    Fail("requires active CGL or adiabatic passive CGL on a uniform mesh");
  }
  for (int face = 0; face < 6; ++face) {
    if (pm->mesh_bcs[face] != BoundaryFlag::periodic) {
      Fail("requires periodic boundaries");
    }
  }
  if (user_srcs) Fail("uses a split user source; leave problem/user_srcs = false");
  CGLLFPaper(pin, restart);
  cooling = Cooling{};
  cooling.mesh = pm;
  cooling.pin = pin;
  cooling.paper_history = user_hist_func;
  cooling.guard_hits = DvceArray1D<std::uint64_t>("cgl_cooling_guard_hits", 1);
  cooling.guard_cycle = pm->ncycle;
  cooling.sample_time = pm->time;
  pgen_final_func = [](ParameterInput *, Mesh *) {
    cooling.guard_hits = DvceArray1D<std::uint64_t>();
  };
  cooling.target = pin->GetReal("problem", "cooling_piso");
  cooling.tau = pin->GetReal("problem", "cooling_tau");
  if (!std::isfinite(cooling.target) || cooling.target <= 0.0 ||
      !std::isfinite(cooling.tau) || cooling.tau <= 0.0) {
    Fail("cooling_piso and cooling_tau must be finite and positive");
  }
  for (const auto &block : pin->block) {
    const auto &name = block.block_name;
    if (name.compare(0, 6, "output") != 0 ||
        !pin->DoesParameterExist(name, "file_type") ||
        pin->GetString(name, "file_type") != "hst") continue;
    if (!cooling.history_block.empty()) Fail("requires exactly one history output");
    cooling.history_block = name;
  }
  if (cooling.history_block.empty()) Fail("requires a history output");
  const auto &hist = cooling.history_block;
  if (pin->DoesParameterExist(hist, "dcycle")) {
    cooling.history_dcycle = pin->GetInteger(hist, "dcycle");
    if (cooling.history_dcycle <= 0) Fail("history dcycle must be positive");
  } else {
    cooling.history_dt = pin->GetReal(hist, "dt");
    if (!std::isfinite(cooling.history_dt) || cooling.history_dt <= 0.0) {
      Fail("history dt must be finite and positive");
    }
  }
  if (restart && pin->DoesParameterExist("problem", "cooling_alpha")) {
    cooling.alpha = pin->GetReal("problem", "cooling_alpha");
    cooling.last_time = pin->GetReal("problem", "cooling_history_last_time");
    cooling.last_cycle = pin->GetInteger("problem", "cooling_history_last_cycle");
    if (!std::isfinite(cooling.alpha) || cooling.alpha < 0.0 ||
        !std::isfinite(cooling.last_time) || cooling.last_time < 0.0 ||
        static_cast<float>(cooling.last_time) > static_cast<float>(pm->time) ||
        cooling.last_cycle < 0 ||
        cooling.last_cycle > pm->ncycle) {
      Fail("invalid cooling controller restart state");
    }
    // Undo final-output-only clock advances before Outputs is constructed.
    SaveReal(pin, hist, "last_time", cooling.last_time);
    const bool has_fraction = pin->DoesParameterExist("problem", "cooling_target_fraction");
    const bool has_time = pin->DoesParameterExist("problem", "cooling_sample_time");
    if (has_fraction != has_time) Fail("incomplete sampled cooling restart state");
    if (has_fraction) {
      cooling.target_fraction = pin->GetReal("problem", "cooling_target_fraction");
      cooling.sample_time = pin->GetReal("problem", "cooling_sample_time");
      if (!std::isfinite(cooling.target_fraction) || cooling.target_fraction < 0.0 ||
          cooling.target_fraction > 1.0 || !std::isfinite(cooling.sample_time) ||
          cooling.sample_time < 0.0 || cooling.sample_time > pm->time) {
        Fail("invalid sampled cooling restart state");
      }
    } else if (!pmbp->frozen_mhd) {
      // Legacy coefficients lack an actual sample time; wait for the next sample.
      cooling.alpha = 0.0;
      SaveReal(pin, "problem", "cooling_alpha", cooling.alpha);
    }
  } else if (restart) {
    // An older uncooled checkpoint starts cooling at its next history sample.
    cooling.last_time = pin->GetOrAddReal(hist, "last_time", pm->time);
    cooling.last_cycle = pm->ncycle;
  }
  user_hist = true;
  SaveReal(pin, "problem", "cooling_target_fraction", cooling.target_fraction);
  SaveReal(pin, "problem", "cooling_sample_time", cooling.sample_time);
  pin->SetBoolean("problem", "user_hist", true);
  user_hist_func = CoolingHistory;
  if (pmbp->frozen_mhd) return;
  // An exact source inside every RK stage would corrupt the RK combination of A.
  // Split once before native collisions, which already refresh restriction/timesteps.
  TaskID none(0);
  const auto task = pmbp->tl_map["after_timeintegrator"]->InsertTask(
      &Cooling::Source, &cooling, none, pmhd->id.cglcoll);
  if (task == none) Fail("requires the standalone MHD collision task");
}

#ifdef CGL_COOLING_USER_PROBLEM
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  CGLLFCooling(pin, restart);
}
#endif
