//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_timed_amr.cpp
//! \brief Turbulence problem generator with a time-dependent AMR refinement hook.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "globals.hpp"

namespace {

Real t_refine = -1.0;

void TurbTimedAMRHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 0;
  return;
}

void TurbTimedAMRRefinementCondition(MeshBlockPack* pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  if (pmesh->pmr == nullptr) return;

  int nmb = pmbp->nmb_thispack;
  int mbs = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;

  static bool refinement_triggered = false;
  if (pmesh->time > t_refine && !refinement_triggered && global_variable::my_rank == 0) {
    std::cout << "### REFINEMENT TRIGGERED at t=" << pmesh->time
              << " (t_refine=" << t_refine << "), nmb=" << nmb << std::endl;
    refinement_triggered = true;
  }

  if (pmesh->time > t_refine) {
    for (int m = 0; m < nmb; ++m) {
      refine_flag.h_view(m + mbs) = 1;
    }
    refine_flag.template modify<HostMemSpace>();
    refine_flag.template sync<DevExeSpace>();
  }
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::TurbTimedAMR()
//! \brief Problem generator for turbulence with time-dependent AMR.

void ProblemGenerator::TurbTimedAMR(ParameterInput *pin, const bool restart) {
  t_refine = pin->GetOrAddReal("problem", "t_refine", 1.0e10);
  user_ref_func = TurbTimedAMRRefinementCondition;
  user_hist_func = TurbTimedAMRHistory;

  if (restart) {
    if (pmy_mesh_->pmb_pack->ppart != nullptr &&
        pin->GetOrAddBoolean("particles", "seed_on_restart", false)) {
      PartRandom(pin, true);
    }
    return;
  }

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Turbulence problem generator can only be run with Hydro and/or MHD, "
              << "but no <hydro> or <mhd> block in input file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int &is = indcs.is, &ie = indcs.ie;
  int &js = indcs.js, &je = indcs.je;
  int &ks = indcs.ks, &ke = indcs.ke;

  Real cs = pin->GetOrAddReal("eos", "iso_sound_speed", 1.0);
  Real beta = pin->GetOrAddReal("problem", "beta", 1.0);

  if (pmbp->phydro != nullptr) {
    Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = rho0*cs*cs/eos.gamma;

    par_for("pgen_turb_timed_hydro", DevExeSpace(), 0, (pmbp->nmb_thispack - 1),
    ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m, IDN, k, j, i) = rho0;
      u0(m, IM1, k, j, i) = 0.0;
      u0(m, IM2, k, j, i) = 0.0;
      u0(m, IM3, k, j, i) = 0.0;
      if (eos.is_ideal) {
        u0(m, IEN, k, j, i) = p0/gm1;
      }
    });
  }

  if (pmbp->pmhd != nullptr) {
    Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = rho0*cs*cs/eos.gamma;
    Real b0_mag = cs*std::sqrt(2.0*p0/beta);

    par_for("pgen_turb_timed_mhd", DevExeSpace(), 0, (pmbp->nmb_thispack - 1),
    ks, ke, js, je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m, IDN, k, j, i) = rho0;
      u0(m, IM1, k, j, i) = 0.0;
      u0(m, IM2, k, j, i) = 0.0;
      u0(m, IM3, k, j, i) = 0.0;

      b0.x1f(m, k, j, i) = 0.0;
      b0.x2f(m, k, j, i) = 0.0;
      b0.x3f(m, k, j, i) = b0_mag;
      if (i == ie) { b0.x1f(m, k, j, i + 1) = 0.0; }
      if (j == je) { b0.x2f(m, k, j + 1, i) = 0.0; }
      if (k == ke) { b0.x3f(m, k + 1, j, i) = b0_mag; }

      if (eos.is_ideal) {
        u0(m, IEN, k, j, i) = p0/gm1 + 0.5*b0_mag*b0_mag;
      }
    });
  }
}
