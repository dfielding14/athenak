//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_drag.cpp
//! \brief drag coupling between massive particles and hydro/MHD fluids

#include <algorithm>
#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "particles.hpp"

namespace particles {

//----------------------------------------------------------------------------------------
//! \fn Particles::ApplyDrag
//! \brief Apply particle drag and optional equal/opposite backreaction during a fluid
//! stage.  The default model is a stopping-time force; model=user dispatches to the pgen.

TaskStatus Particles::ApplyDrag(Driver *pdriver, int stage) {
  if (!drag_enabled) {return TaskStatus::complete;}

  Real beta_dt = (pdriver->beta[stage-1])*(pmy_pack->pmesh->dt);
  if (drag_model == DragParticlesModel::stopping_time) {
    ApplyStoppingTimeDrag(beta_dt);
  } else if (drag_model == DragParticlesModel::user) {
    if ((pmy_pack->pmesh->pgen == nullptr) ||
        (pmy_pack->pmesh->pgen->user_particle_drag_func == nullptr)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<drag_particles>/model=user requires a pgen-enrolled "
                << "particle drag callback" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    (pmy_pack->pmesh->pgen->user_particle_drag_func)(pmy_pack, beta_dt);
  }

  if ((pmy_pack->phydro != nullptr) && (pmy_pack->pmhd == nullptr)) {
    (void) pmy_pack->phydro->ConToPrim(pdriver, stage);
  } else if ((pmy_pack->pmhd != nullptr) && (pmy_pack->phydro == nullptr)) {
    (void) pmy_pack->pmhd->ConToPrim(pdriver, stage);
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn Particles::ApplyStoppingTimeDrag
//! \brief Standard drag law dv_p/dt = (v_g - v_p)/t_stop with user-selected
//! cell-to-particle interpolation and particle-to-cell deposition stencils.

void Particles::ApplyStoppingTimeDrag(const Real bdt) {
  if (nprtcl_thispack <= 0) {return;}

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  int npart = nprtcl_thispack;
  auto gids = pmy_pack->gids;
  int nmb = pmy_pack->nmb_thispack;

  auto &mbsize = pmy_pack->pmb->mb_size;
  auto pr = prtcl_rdata;
  auto pi = prtcl_idata;

  bool use_mhd = (pmy_pack->pmhd != nullptr);
  auto w0 = use_mhd ? pmy_pack->pmhd->w0 : pmy_pack->phydro->w0;
  auto u0 = use_mhd ? pmy_pack->pmhd->u0 : pmy_pack->phydro->u0;
  EOS_Data eos_data = use_mhd ? pmy_pack->pmhd->peos->eos_data
                              : pmy_pack->phydro->peos->eos_data;

  Real stopping_time = drag_stopping_time;
  Real default_mass = drag_particle_mass;
  bool backreaction = drag_backreaction;
  bool include_energy = drag_include_energy && eos_data.is_ideal;
  bool orbital_terms = drag_orbital_terms;
  DragParticlesCoupling interpolation = drag_interpolation;
  DragParticlesCoupling deposition = drag_deposition;
  Real omega0 = drag_omega0;
  Real qshear = drag_qshear;

  par_for("particle_drag_stopping_time", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    if ((m < 0) || (m >= nmb)) {return;}

    Real pmass = pr(IPM,p);
    if (pmass <= 0.0) {
      pmass = default_mass;
      pr(IPM,p) = pmass;
    }

    Real vx_old = pr(IPVX,p);
    Real vy_old = multi_d ? pr(IPVY,p) : 0.0;
    Real vz_old = three_d ? pr(IPVZ,p) : 0.0;

    Real vx_pred = vx_old;
    Real vy_pred = vy_old;
    Real vz_pred = vz_old;
    if (orbital_terms) {
      vx_pred += bdt*(2.0*omega0*vy_old);
      vy_pred += bdt*(-(2.0 - qshear)*omega0*vx_old);
    }

    DragParticleStencil interp_stencil;
    BuildDragParticleStencil(interpolation, pr(IPX,p), pr(IPY,p), pr(IPZ,p),
                             mbsize.d_view(m).x1min, mbsize.d_view(m).x2min,
                             mbsize.d_view(m).x3min, mbsize.d_view(m).dx1,
                             mbsize.d_view(m).dx2, mbsize.d_view(m).dx3, is, ie, js,
                             je, ks, ke, multi_d, three_d, interp_stencil);
    Real vg1, vg2, vg3;
    InterpolateDragVelocity(w0, m, interp_stencil, multi_d, three_d, vg1, vg2, vg3);
    Real drag_factor = 1.0 - exp(-bdt/stopping_time);

    Real dv1_drag = drag_factor*(vg1 - vx_pred);
    Real dv2_drag = multi_d ? drag_factor*(vg2 - vy_pred) : 0.0;
    Real dv3_drag = three_d ? drag_factor*(vg3 - vz_pred) : 0.0;

    Real vx_new = vx_pred + dv1_drag;
    Real vy_new = vy_pred + dv2_drag;
    Real vz_new = vz_pred + dv3_drag;

    pr(IPVX,p) = vx_new;
    if (multi_d) {pr(IPVY,p) = vy_new;}
    if (three_d) {pr(IPVZ,p) = vz_new;}

    if (backreaction) {
      Real vol = mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*mbsize.d_view(m).dx3;
      Real inv_vol = 1.0/vol;
      Real denergy = 0.0;
      if (include_energy) {
        Real ke_old = 0.5*pmass*(vx_pred*vx_pred + vy_pred*vy_pred + vz_pred*vz_pred);
        Real ke_new = 0.5*pmass*(vx_new*vx_new + vy_new*vy_new + vz_new*vz_new);
        denergy = -(ke_new - ke_old);
      }
      DragParticleStencil deposit_stencil;
      BuildDragParticleStencil(deposition, pr(IPX,p), pr(IPY,p), pr(IPZ,p),
                               mbsize.d_view(m).x1min, mbsize.d_view(m).x2min,
                               mbsize.d_view(m).x3min, mbsize.d_view(m).dx1,
                               mbsize.d_view(m).dx2, mbsize.d_view(m).dx3, is, ie, js,
                               je, ks, ke, multi_d, three_d, deposit_stencil);
      DepositDragBackreaction(u0, m, deposit_stencil, inv_vol, -pmass*dv1_drag,
                              -pmass*dv2_drag, -pmass*dv3_drag, denergy, multi_d,
                              three_d, include_energy);
    }
  });
}

} // namespace particles
