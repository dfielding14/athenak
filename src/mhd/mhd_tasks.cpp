//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_tasks.cpp
//! \brief functions that control MHD tasks stored in tasklists in MeshBlockPack

#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <cmath>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/resistivity.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "bvals/bvals.hpp"
#include "shearing_box/shearing_box.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "particles/particles.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void MHD::AssembleMHDTasks
//! \brief Adds mhd tasks to appropriate task lists used by time integrators.
//! Called by MeshBlockPack::AddPhysics() function directly after MHD constructor
//! See comments Hydro::AssembleHydroTasks() function for more details.

void MHD::AssembleMHDTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // assemble "before_timeintegrator" task list
  id.savest = tl["before_timeintegrator"]->AddTask(&MHD::SaveMHDState, this, none);

  // assemble "before_stagen" task list
  id.irecv = tl["before_stagen"]->AddTask(&MHD::InitRecv, this, none);

  // assemble "stagen" task list
  id.copyu     = tl["stagen"]->AddTask(&MHD::CopyCons, this, none);
  id.flux      = tl["stagen"]->AddTask(&MHD::Fluxes, this, id.copyu);
  id.sendf     = tl["stagen"]->AddTask(&MHD::SendFlux, this, id.flux);
  id.recvf     = tl["stagen"]->AddTask(&MHD::RecvFlux, this, id.sendf);
  id.rkupdt    = tl["stagen"]->AddTask(&MHD::RKUpdate, this, id.recvf);
  id.srctrms   = tl["stagen"]->AddTask(&MHD::MHDSrcTerms, this, id.rkupdt);
  id.picdampu  = tl["stagen"]->AddTask(&MHD::ApplyPICWaveDamping, this, id.srctrms);
  id.efld      = tl["stagen"]->AddTask(&MHD::CornerE, this, id.picdampu);
  id.efldsrc   = tl["stagen"]->AddTask(&MHD::EFieldSrc, this, id.efld);
  id.sende     = tl["stagen"]->AddTask(&MHD::SendE, this, id.efldsrc);
  id.recve     = tl["stagen"]->AddTask(&MHD::RecvE, this, id.sende);
  id.ct        = tl["stagen"]->AddTask(&MHD::CT, this, id.recve);
  id.expboxb   = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxB, this, id.ct);
  id.expboxu   = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxU, this, id.expboxb);
  id.expboxfb  = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxFeedback, this,
                                      id.expboxu);
  id.expboxdampu = tl["stagen"]->AddTask(&MHD::ApplyPICExpandingBoxWaveDamping, this,
                                        id.expboxfb);
  id.sendu_oa  = tl["stagen"]->AddTask(&MHD::SendU_OA, this, id.expboxdampu);
  id.recvu_oa  = tl["stagen"]->AddTask(&MHD::RecvU_OA, this, id.sendu_oa);
  id.restu     = tl["stagen"]->AddTask(&MHD::RestrictU, this, id.recvu_oa);
  id.sendu     = tl["stagen"]->AddTask(&MHD::SendU, this, id.restu);
  id.recvu     = tl["stagen"]->AddTask(&MHD::RecvU, this, id.sendu);
  id.sendu_shr = tl["stagen"]->AddTask(&MHD::SendU_Shr, this, id.recvu);
  id.recvu_shr = tl["stagen"]->AddTask(&MHD::RecvU_Shr, this, id.sendu_shr);
  id.sendb_oa  = tl["stagen"]->AddTask(&MHD::SendB_OA, this, id.recvu_shr);
  id.recvb_oa  = tl["stagen"]->AddTask(&MHD::RecvB_OA, this, id.sendb_oa);
  id.restb     = tl["stagen"]->AddTask(&MHD::RestrictB, this, id.recvb_oa);
  id.sendb     = tl["stagen"]->AddTask(&MHD::SendB, this, id.restb);
  id.recvb     = tl["stagen"]->AddTask(&MHD::RecvB, this, id.sendb);
  id.sendb_shr = tl["stagen"]->AddTask(&MHD::SendB_Shr, this, id.recvb);
  id.recvb_shr = tl["stagen"]->AddTask(&MHD::RecvB_Shr, this, id.sendb_shr);
  id.bcs       = tl["stagen"]->AddTask(&MHD::ApplyPhysicalBCs, this, id.recvb_shr);
  id.prol      = tl["stagen"]->AddTask(&MHD::Prolongate, this, id.bcs);
  id.c2p       = tl["stagen"]->AddTask(&MHD::ConToPrim, this, id.prol);
  id.newdt     = tl["stagen"]->AddTask(&MHD::NewTimeStep, this, id.c2p);

  // assemble "after_stagen" task list
  id.csend = tl["after_stagen"]->AddTask(&MHD::ClearSend, this, none);
  // although RecvFlux/U/E/B functions check that all recvs complete, add ClearRecv to
  // task list anyways to catch potential bugs in MPI communication logic
  id.crecv = tl["after_stagen"]->AddTask(&MHD::ClearRecv, this, id.csend);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::SaveMHDState
//! \brief Copy primitives and bcc before step to enable computation of time derivatives,
//! for example to compute jcon in GRMHD.

TaskStatus MHD::SaveMHDState(Driver *pdrive, int stage) {
  if (wbcc_saved) {
    Kokkos::deep_copy(DevExeSpace(), wsaved, w0);
    Kokkos::deep_copy(DevExeSpace(), bccsaved, bcc0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI), and
//! initialize all boundary receive status flags to waiting (with or without MPI).  Note
//! this must be done for communication of BOTH conserved (cell-centered) and
//! face-centered fields AND their fluxes (with SMR/AMR).

TaskStatus MHD::InitRecv(Driver *pdrive, int stage) {
  // post receives for U
  TaskStatus tstat = pbval_u->InitRecv(nmhd+nscalars);
  if (tstat != TaskStatus::complete) return tstat;
  // post receives for B
  tstat = pbval_b->InitRecv(3);
  if (tstat != TaskStatus::complete) return tstat;

  // with SMR/AMR post receives for fluxes of U, always post receives for fluxes of B
  // do not post receives for fluxes when stage < 0 (i.e. ICs)
  if (stage >= 0) {
    // with SMR/AMR, post receives for fluxes of U
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->InitFluxRecv(nmhd+nscalars);
      if (tstat != TaskStatus::complete) return tstat;
    }
    // post receives for fluxes of B, which are used even with uniform grids
    tstat = pbval_b->InitFluxRecv(3);
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with orbital advection post receives for U and B
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = porb_u->InitRecv();
    if (tstat != TaskStatus::complete) return tstat;
    tstat = porb_b->InitRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with shearing box boundaries caluclate x2-distance x1-boundaries have sheared and
  // with MPI post receives for U and B
  // only execute when (shearing box defined) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    Real qom = (psrc->qshear)*(psrc->omega0);
    Real time = pmy_pack->pmesh->time;
    if (stage == pdrive->nexp_stages) {
      time += pmy_pack->pmesh->dt;
    }
    tstat = psbox_u->InitRecv(qom, time);
    if (tstat != TaskStatus::complete) return tstat;
    tstat = psbox_b->InitRecv(qom, time);
    if (tstat != TaskStatus::complete) return tstat;
  }

  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::CopyCons
//! \brief Simple task list function that copies u0 --> u1, and b0 --> b1 in first stage

TaskStatus MHD::CopyCons(Driver *pdrive, int stage) {
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
    Kokkos::deep_copy(DevExeSpace(), b1.x1f, b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), b1.x2f, b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), b1.x3f, b0.x3f);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::Fluxes
//! \brief Wrapper task list function that calls everything necessary to compute fluxes
//! of conserved variables

TaskStatus MHD::Fluxes(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart != nullptr) && ppart->UsesExpandingBox()) {
    RefreshPICExpandingBoxPhysicalB(pmy_pack->pmesh->time);
  }

  // select which calculate_flux function to call based on rsolver_method
  if (rsolver_method == MHD_RSolver::advect) {
    CalculateFluxes<MHD_RSolver::advect>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::llf) {
    CalculateFluxes<MHD_RSolver::llf>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::hlle) {
    CalculateFluxes<MHD_RSolver::hlle>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::hlld) {
    CalculateFluxes<MHD_RSolver::hlld>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::llf_sr) {
    CalculateFluxes<MHD_RSolver::llf_sr>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::hlle_sr) {
    CalculateFluxes<MHD_RSolver::hlle_sr>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::llf_gr) {
    CalculateFluxes<MHD_RSolver::llf_gr>(pdrive, stage);
  } else if (rsolver_method == MHD_RSolver::hlle_gr) {
    CalculateFluxes<MHD_RSolver::hlle_gr>(pdrive, stage);
  }

  // Add viscous, resistive, heat-flux, etc fluxes
  if (pvisc != nullptr) {
    pvisc->IsotropicViscousFlux(w0, pvisc->nu_iso, peos->eos_data, uflx);
  }
  if ((presist != nullptr) && (peos->eos_data.is_ideal)) {
    presist->OhmicEnergyFlux(b0, uflx);
  }
  if (pcond != nullptr) {
    pcond->AddHeatFlux(w0, peos->eos_data, uflx);
  }

  // call FOFC if necessary
  if (use_fofc) {
    FOFC(pdrive, stage);
  } else if (pmy_pack->pcoord->is_general_relativistic) {
    if (pmy_pack->pcoord->coord_data.bh_excise) {
      FOFC(pdrive, stage);
    }
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::RefreshPICExpandingBoxPhysicalB
//! \brief Populate physical face fields from divergence-preserving comoving face fluxes.

void MHD::RefreshPICExpandingBoxPhysicalB(const Real time) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart == nullptr) || !ppart->UsesExpandingBox()) return;
  bphys_time = time;

  const auto geom = particles::PICExpandingBoxGeometryAt(
      ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
      ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3, time);
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto b1 = b0.x1f;
  auto b2 = b0.x2f;
  auto b3 = b0.x3f;
  auto bp1 = bphys.x1f;
  auto bp2 = bphys.x2f;
  auto bp3 = bphys.x3f;
  const int b1k = static_cast<int>(b1.extent(1)) - 1;
  const int b1j = static_cast<int>(b1.extent(2)) - 1;
  const int b1i = static_cast<int>(b1.extent(3)) - 1;
  const int b2k = static_cast<int>(b2.extent(1)) - 1;
  const int b2j = static_cast<int>(b2.extent(2)) - 1;
  const int b2i = static_cast<int>(b2.extent(3)) - 1;
  const int b3k = static_cast<int>(b3.extent(1)) - 1;
  const int b3j = static_cast<int>(b3.extent(2)) - 1;
  const int b3i = static_cast<int>(b3.extent(3)) - 1;
  par_for("pic_expanding_box_physical_b1", DevExeSpace(), 0, nmb1,
          0, b1k, 0, b1j, 0, b1i,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    bp1(m, k, j, i) = geom.inv_area1*b1(m, k, j, i);
  });
  par_for("pic_expanding_box_physical_b2", DevExeSpace(), 0, nmb1,
          0, b2k, 0, b2j, 0, b2i,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    bp2(m, k, j, i) = geom.inv_area2*b2(m, k, j, i);
  });
  par_for("pic_expanding_box_physical_b3", DevExeSpace(), 0, nmb1,
          0, b3k, 0, b3j, 0, b3i,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    bp3(m, k, j, i) = geom.inv_area3*b3(m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::SendFlux
//! \brief Wrapper task list function to pack/send restricted values of fluxes of
//! conserved variables at fine/coarse boundaries

TaskStatus MHD::SendFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // Only execute BoundaryValues function with SMR/SMR
  if (pmy_pack->pmesh->multilevel)  {
    tstat = pbval_u->PackAndSendFluxCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RecvFlux
//! \brief Wrapper task list function to recv/unpack restricted values of fluxes of
//! conserved variables at fine/coarse boundaries

TaskStatus MHD::RecvFlux(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // Only execute BoundaryValues function with SMR/SMR
  if (pmy_pack->pmesh->multilevel) {
    tstat = pbval_u->RecvAndUnpackFluxCC(uflx);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::MHDSrcTerms
//! \brief Wrapper task list function to apply source terms to conservative vars
//! Note source terms must be computed using only primitives (w0), as the conserved
//! variables (u0) have already been partially updated when this fn called.

TaskStatus MHD::MHDSrcTerms(Driver *pdrive, int stage) {
  Real beta_dt = (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt);
  auto *ppart = pmy_pack->ppart;
  if ((ppart != nullptr) &&
      (ppart->pic_background_mode == PICBackgroundMode::passive_mhd)) {
    return TaskStatus::complete;
  }

  // Add source terms for various physics
  if (psrc->const_accel)  psrc->ConstantAccel(w0, peos->eos_data, beta_dt, u0);
  if (psrc->ism_cooling)  psrc->ISMCooling(w0, peos->eos_data, beta_dt, u0);
  if (psrc->cgm_cooling)  psrc->CGMCooling(w0, peos->eos_data, beta_dt, u0);
  if (psrc->rel_cooling)  psrc->RelCooling(w0, peos->eos_data, beta_dt, u0);
  if (psrc->shearing_box) psrc->ShearingBox(w0, bcc0, peos->eos_data, beta_dt, u0);

  // Add coordinate source terms in GR.  Again, must be computed with only primitives.
  if (pmy_pack->pcoord->is_general_relativistic &&
      !pmy_pack->pcoord->is_dynamical_relativistic) {
    pmy_pack->pcoord->CoordSrcTerms(w0, bcc0, peos->eos_data, beta_dt, u0);
  } else if (pmy_pack->pcoord->is_dynamical_relativistic) {
    pmy_pack->pdyngr->AddCoordTerms(w0, bcc0, beta_dt, u0, pmy_pack->pmesh->mb_indcs.ng);
  }

  // Add user source terms
  if (pmy_pack->pmesh->pgen->user_srcs) {
    (pmy_pack->pmesh->pgen->user_srcs_func)(pmy_pack->pmesh, beta_dt);
  }

  // WS-H: optional particle-feedback source split for fluid momentum/energy.
  // Legacy coupled modes reuse cycle-fixed deposited rates across RK stages.
  // Paper VL2 uses predictor rho/J in stage 1 and the final particle impulse
  // deposited after the stage-2 midpoint kick.
  if ((ppart != nullptr) && ppart->couple_moments_to_mhd) {
    const bool add_mom = ppart->couple_moments_momentum_to_mhd;
    const bool add_eng = ppart->couple_moments_energy_to_mhd;
    const bool apply_feedback_here = (ppart->couple_fluid_feedback_order ==
                                      CoupledFluidFeedbackOrder::mhd_src_terms);
    const bool use_delta_feedback =
        ((ppart->pusher == ParticlesPusher::boris_lin) ||
         (ppart->pusher == ParticlesPusher::boris_tsc)) &&
        (ppart->pic_feedback_mode == PICFeedbackMode::coupled);
    if ((add_mom || add_eng) && apply_feedback_here && !ppart->UsesExpandingBox()) {
      if (add_eng && (nmhd <= IEN)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Particle energy feedback requires an ideal-MHD energy variable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is, ie = indcs.ie;
      int js = indcs.js, je = indcs.je;
      int ks = indcs.ks, ke = indcs.ke;
      int nmb1 = pmy_pack->nmb_thispack - 1;
      const Real mom_coef = ppart->couple_moments_momentum_coeff;
      const Real eng_coef = ppart->couple_moments_energy_coeff;
      const bool use_deltaf = ppart->UsesDeltaF();
      const bool paper_vl2_predictor =
          ppart->UsesPaperVL2Coupling() && (stage == 1);
      Real background_density_scale = 1.0;
      if (ppart->UsesExpandingBox()) {
        const auto geom = particles::PICExpandingBoxGeometryAt(
            ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
            ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3,
            pmy_pack->pmesh->time + pmy_pack->pmesh->dt);
        background_density_scale = geom.inv_a1*geom.inv_a2*geom.inv_a3;
      }
      const Real background_rho =
          background_density_scale*ppart->pic_deltaf_background_rho;
      const Real background_jx = ppart->pic_deltaf_background_jx;
      const Real background_jy = ppart->pic_deltaf_background_jy;
      const Real background_jz = ppart->pic_deltaf_background_jz;
      auto mom = ppart->moments;
      auto bcc = bcc0;
      auto w = w0;
      auto u = u0;

      par_for("prtcl_fluid_feedback_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (paper_vl2_predictor && !use_deltaf) {
          const Real rho = mom(m, particles::Particles::IMOM_RHO, k, j, i);
          const Real jx = mom(m, particles::Particles::IMOM_JX, k, j, i);
          const Real jy = mom(m, particles::Particles::IMOM_JY, k, j, i);
          const Real jz = mom(m, particles::Particles::IMOM_JZ, k, j, i);
          const Real bx = bcc(m, IBX, k, j, i);
          const Real by = bcc(m, IBY, k, j, i);
          const Real bz = bcc(m, IBZ, k, j, i);
          const Real cex = -(w(m, IVY, k, j, i)*bz - w(m, IVZ, k, j, i)*by);
          const Real cey = -(w(m, IVZ, k, j, i)*bx - w(m, IVX, k, j, i)*bz);
          const Real cez = -(w(m, IVX, k, j, i)*by - w(m, IVY, k, j, i)*bx);
          if (add_mom) {
            u(m, IM1, k, j, i) -= beta_dt*mom_coef*(rho*cex + jy*bz - jz*by);
            u(m, IM2, k, j, i) -= beta_dt*mom_coef*(rho*cey + jz*bx - jx*bz);
            u(m, IM3, k, j, i) -= beta_dt*mom_coef*(rho*cez + jx*by - jy*bx);
          }
          if (add_eng) {
            u(m, IEN, k, j, i) -= beta_dt*eng_coef*(jx*cex + jy*cey + jz*cez);
          }
        } else if (use_deltaf) {
          const Real rho = background_rho +
              mom(m, particles::Particles::IMOM_RHO, k, j, i);
          const Real jx = background_jx +
              mom(m, particles::Particles::IMOM_JX, k, j, i);
          const Real jy = background_jy +
              mom(m, particles::Particles::IMOM_JY, k, j, i);
          const Real jz = background_jz +
              mom(m, particles::Particles::IMOM_JZ, k, j, i);
          const Real bx = bcc(m, IBX, k, j, i);
          const Real by = bcc(m, IBY, k, j, i);
          const Real bz = bcc(m, IBZ, k, j, i);
          const Real cex = -(w(m, IVY, k, j, i)*bz - w(m, IVZ, k, j, i)*by);
          const Real cey = -(w(m, IVZ, k, j, i)*bx - w(m, IVX, k, j, i)*bz);
          const Real cez = -(w(m, IVX, k, j, i)*by - w(m, IVY, k, j, i)*bx);
          if (add_mom) {
            u(m, IM1, k, j, i) -= beta_dt*mom_coef*(rho*cex + jy*bz - jz*by);
            u(m, IM2, k, j, i) -= beta_dt*mom_coef*(rho*cey + jz*bx - jx*bz);
            u(m, IM3, k, j, i) -= beta_dt*mom_coef*(rho*cez + jx*by - jy*bx);
          }
          if (add_eng) {
            u(m, IEN, k, j, i) -= beta_dt*eng_coef*(jx*cex + jy*cey + jz*cez);
          }
        } else if (use_delta_feedback) {
          if (add_mom) {
            u(m, IM1, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPXDT, k, j, i);
            u(m, IM2, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPYDT, k, j, i);
            u(m, IM3, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPZDT, k, j, i);
          }
          if (add_eng) {
            u(m, IEN, k, j, i) -= beta_dt*eng_coef*
                                  mom(m, particles::Particles::IMOM_DEDT, k, j, i);
          }
        } else {
          Real jx = mom(m, particles::Particles::IMOM_JX, k, j, i);
          Real jy = mom(m, particles::Particles::IMOM_JY, k, j, i);
          Real jz = mom(m, particles::Particles::IMOM_JZ, k, j, i);

          Real bx = bcc(m, IBX, k, j, i);
          Real by = bcc(m, IBY, k, j, i);
          Real bz = bcc(m, IBZ, k, j, i);

          Real fx = -(jy*bz - jz*by);
          Real fy = -(jz*bx - jx*bz);
          Real fz = -(jx*by - jy*bx);

          if (add_mom) {
            u(m, IM1, k, j, i) += beta_dt*mom_coef*fx;
            u(m, IM2, k, j, i) += beta_dt*mom_coef*fy;
            u(m, IM3, k, j, i) += beta_dt*mom_coef*fz;
          }
          if (add_eng) {
            u(m, IEN, k, j, i) += beta_dt*eng_coef*(jx*bx + jy*by + jz*bz);
          }
        }
      });
    }
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICWaveDamping
//! \brief Apply reduced high-frequency ion-neutral friction to transverse ion momentum.

TaskStatus MHD::ApplyPICWaveDamping(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart != nullptr) && ppart->UsesExpandingBox()) {
    return TaskStatus::complete;
  }
  return ApplyPICWaveDampingMap(pdrive, stage);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICExpandingBoxWaveDamping
//! \brief Apply reduced wave damping after final expanding-box physical-frame sources.

TaskStatus MHD::ApplyPICExpandingBoxWaveDamping(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart == nullptr) || !ppart->UsesExpandingBox()) {
    return TaskStatus::complete;
  }
  return ApplyPICWaveDampingMap(pdrive, stage);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICWaveDampingMap
//! \brief Apply the exact reduced high-frequency ion-neutral friction map.

TaskStatus MHD::ApplyPICWaveDampingMap(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart == nullptr) || !ppart->UsesPICWaveDamping() ||
      (stage != pdrive->nexp_stages)) {
    return TaskStatus::complete;
  }
  const Real factor = exp(-ppart->pic_ion_neutral_collision_rate*
                          pmy_pack->pmesh->dt);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool ideal = peos->eos_data.is_ideal;
  auto u = u0;
  par_for("pic_ion_neutral_friction", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real den = u(m, IDN, k, j, i);
    const Real my0 = u(m, IM2, k, j, i);
    const Real mz0 = u(m, IM3, k, j, i);
    u(m, IM2, k, j, i) *= factor;
    u(m, IM3, k, j, i) *= factor;
    if (ideal) {
      const Real transverse_ekin_loss =
          0.5*(1.0 - factor*factor)*(my0*my0 + mz0*mz0)/den;
      u(m, IEN, k, j, i) -= transverse_ekin_loss;
    }
  });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICExpandingBoxU
//! \brief Rescale fluid conserved variables for a uniform anisotropic expanding box.

TaskStatus MHD::ApplyPICExpandingBoxU(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart == nullptr) || !ppart->UsesExpandingBox() ||
      (stage != pdrive->nexp_stages)) {
    return TaskStatus::complete;
  }
  // Apply the exact expansion flow once after the RK and CT updates. Treating
  // this split flow as a stage-local source map over-expands multistage
  // integrators.
  const Real bdt = pmy_pack->pmesh->dt;
  const Real t0 = pmy_pack->pmesh->time;
  const Real t1 = t0 + bdt;
  const auto geom0 = particles::PICExpandingBoxGeometryAt(
      ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
      ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3, t0);
  const auto geom1 = particles::PICExpandingBoxGeometryAt(
      ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
      ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3, t1);
  const Real r1 = geom0.a1/geom1.a1;
  const Real r2 = geom0.a2/geom1.a2;
  const Real r3 = geom0.a3/geom1.a3;
  const Real rvol = r1*r2*r3;
  const bool ideal = peos->eos_data.is_ideal;
  const Real gamma = peos->eos_data.gamma;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int scalar_start = nmhd;
  const int scalar_end = nmhd + nscalars;
  auto u = u0;
  auto b1 = b0.x1f;
  auto b2 = b0.x2f;
  auto b3 = b0.x3f;

  par_for("pic_expanding_box_u", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real den0 = u(m, IDN, k, j, i);
    const Real mx0 = u(m, IM1, k, j, i);
    const Real my0 = u(m, IM2, k, j, i);
    const Real mz0 = u(m, IM3, k, j, i);
    const Real den1 = rvol*den0;
    const Real mx1 = rvol*r1*mx0;
    const Real my1 = rvol*r2*my0;
    const Real mz1 = rvol*r3*mz0;
    if (ideal) {
      const Real bx = 0.5*(b1(m, k, j, i) + b1(m, k, j, i+1));
      const Real by = 0.5*(b2(m, k, j, i) + b2(m, k, j+1, i));
      const Real bz = 0.5*(b3(m, k, j, i) + b3(m, k+1, j, i));
      const Real ekin0 = 0.5*(mx0*mx0 + my0*my0 + mz0*mz0)/den0;
      const Real emag0 = 0.5*(bx*bx*geom0.inv_area1*geom0.inv_area1 +
                              by*by*geom0.inv_area2*geom0.inv_area2 +
                              bz*bz*geom0.inv_area3*geom0.inv_area3);
      const Real etherm0 = u(m, IEN, k, j, i) - ekin0 - emag0;
      const Real ekin1 = 0.5*(mx1*mx1 + my1*my1 + mz1*mz1)/den1;
      const Real emag1 = 0.5*(bx*bx*geom1.inv_area1*geom1.inv_area1 +
                              by*by*geom1.inv_area2*geom1.inv_area2 +
                              bz*bz*geom1.inv_area3*geom1.inv_area3);
      u(m, IEN, k, j, i) = pow(rvol, gamma)*etherm0 + ekin1 + emag1;
    }
    u(m, IDN, k, j, i) = den1;
    u(m, IM1, k, j, i) = mx1;
    u(m, IM2, k, j, i) = my1;
    u(m, IM3, k, j, i) = mz1;
    for (int n = scalar_start; n < scalar_end; ++n) {
      u(m, n, k, j, i) *= rvol;
    }
  });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICExpandingBoxFeedback
//! \brief Apply qualified feedback after mapping gas state to the final physical frame.

TaskStatus MHD::ApplyPICExpandingBoxFeedback(Driver *pdrive, int stage) {
  auto *ppart = pmy_pack->ppart;
  if ((ppart == nullptr) || !ppart->UsesExpandingBox() ||
      !ppart->couple_moments_to_mhd || (stage != pdrive->nexp_stages)) {
    return TaskStatus::complete;
  }
  const bool add_mom = ppart->couple_moments_momentum_to_mhd;
  const bool add_eng = ppart->couple_moments_energy_to_mhd;
  if (!(add_mom || add_eng)) return TaskStatus::complete;

  const Real dt = pmy_pack->pmesh->dt;
  const Real mom_coef = ppart->couple_moments_momentum_coeff;
  const Real eng_coef = ppart->couple_moments_energy_coeff;
  const bool use_deltaf = ppart->UsesDeltaF();
  const auto geom = particles::PICExpandingBoxGeometryAt(
      ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
      ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3,
      pmy_pack->pmesh->time + dt);
  const Real background_rho =
      geom.inv_a1*geom.inv_a2*geom.inv_a3*ppart->pic_deltaf_background_rho;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto mom = ppart->moments;
  auto u = u0;
  auto b1 = b0.x1f;
  auto b2 = b0.x2f;
  auto b3 = b0.x3f;

  par_for("prtcl_expanding_box_feedback", DevExeSpace(), 0, nmb1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (use_deltaf) {
      const Real den = u(m, IDN, k, j, i);
      const Real vx = u(m, IM1, k, j, i)/den;
      const Real vy = u(m, IM2, k, j, i)/den;
      const Real vz = u(m, IM3, k, j, i)/den;
      const Real bx = static_cast<Real>(0.5)*
          (b1(m, k, j, i) + b1(m, k, j, i+1))*geom.inv_area1;
      const Real by = static_cast<Real>(0.5)*
          (b2(m, k, j, i) + b2(m, k, j+1, i))*geom.inv_area2;
      const Real bz = static_cast<Real>(0.5)*
          (b3(m, k, j, i) + b3(m, k+1, j, i))*geom.inv_area3;
      const Real cex = -(vy*bz - vz*by);
      const Real cey = -(vz*bx - vx*bz);
      const Real cez = -(vx*by - vy*bx);
      const Real rho = background_rho +
          mom(m, particles::Particles::IMOM_RHO, k, j, i);
      const Real jx = mom(m, particles::Particles::IMOM_JX, k, j, i);
      const Real jy = mom(m, particles::Particles::IMOM_JY, k, j, i);
      const Real jz = mom(m, particles::Particles::IMOM_JZ, k, j, i);
      if (add_mom) {
        u(m, IM1, k, j, i) -= dt*mom_coef*(rho*cex + jy*bz - jz*by);
        u(m, IM2, k, j, i) -= dt*mom_coef*(rho*cey + jz*bx - jx*bz);
        u(m, IM3, k, j, i) -= dt*mom_coef*(rho*cez + jx*by - jy*bx);
      }
      if (add_eng) {
        u(m, IEN, k, j, i) -= dt*eng_coef*(jx*cex + jy*cey + jz*cez);
      }
    } else {
      if (add_mom) {
        u(m, IM1, k, j, i) -= dt*mom_coef*
                              mom(m, particles::Particles::IMOM_DPXDT, k, j, i);
        u(m, IM2, k, j, i) -= dt*mom_coef*
                              mom(m, particles::Particles::IMOM_DPYDT, k, j, i);
        u(m, IM3, k, j, i) -= dt*mom_coef*
                              mom(m, particles::Particles::IMOM_DPZDT, k, j, i);
      }
      if (add_eng) {
        u(m, IEN, k, j, i) -= dt*eng_coef*
                              mom(m, particles::Particles::IMOM_DEDT, k, j, i);
      }
    }
  });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPICExpandingBoxB
//! \brief Keep raw face fields as divergence-preserving comoving magnetic fluxes.

TaskStatus MHD::ApplyPICExpandingBoxB(Driver *pdrive, int stage) {
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::SendU_OA
//! \brief Wrapper task list function to pack/send data for orbital advection

TaskStatus MHD::SendU_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = porb_u->PackAndSendCC(u0);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::RecvU_OA
//! \brief Wrapper task list function to recv/unpack data for orbital advection

TaskStatus MHD::RecvU_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    Real qom = (psrc->qshear)*(psrc->omega0);
    tstat = porb_u->RecvAndUnpackCC(u0, recon_method, qom);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RestrictU
//! \brief Wrapper task list function to restrict conserved vars

TaskStatus MHD::RestrictU(Driver *pdrive, int stage) {
  // Only execute Mesh function with SMR/AMR
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(u0, coarse_u0);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::SendU
//! \brief Wrapper task list function to pack/send cell-centered conserved variables

TaskStatus MHD::SendU(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->PackAndSendCC(u0, coarse_u0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RecvU
//! \brief Wrapper task list function to receive/unpack cell-centered conserved variables

TaskStatus MHD::RecvU(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_u->RecvAndUnpackCC(u0, coarse_u0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::SendU_Shr
//! \brief Wrapper task list function to pack/send data for shearing box boundaries

TaskStatus MHD::SendU_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_u->PackAndSendCC(u0, recon_method);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::RecvU_Shr
//! \brief Wrapper task list function to recv/unpack data for shearing box boundaries
//! Orbital remap is performed in this step.

TaskStatus MHD::RecvU_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_u->RecvAndUnpackCC(u0);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::EFieldSrc
//! \brief Wrapper task list function to apply source terms to electric field

TaskStatus MHD::EFieldSrc(Driver *pdrive, int stage) {
  // only execute when (shearing box defined) AND (2D)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->two_d)) {
    psrc->SBoxEField(b0, efld);
  }

  // PR2: add deposited particle current to edge-centered electric fields.
  auto *ppart = pmy_pack->ppart;
  if ((ppart != nullptr) && ppart->AddsCRCurrentToCT()) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int is = indcs.is, ie = indcs.ie;
    int js = indcs.js, je = indcs.je;
    int ks = indcs.ks, ke = indcs.ke;
    int nmb1 = pmy_pack->nmb_thispack - 1;
    const Real jcoef = ppart->couple_j_to_efield_coeff;

    if (ppart->couple_j_to_efield_representation ==
        CoupledCurrentRepresentation::edge_staggered) {
      auto jx_e = ppart->j_edge_x1e;
      auto jy_e = ppart->j_edge_x2e;
      auto jz_e = ppart->j_edge_x3e;

      if (pmy_pack->pmesh->one_d) {
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_edge_1d", DevExeSpace(), 0, nmb1, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int i) {
          e2(m,ks  ,js  ,i) += jcoef*jy_e(m,ks  ,js,i);
          e2(m,ke+1,js  ,i) += jcoef*jy_e(m,ke+1,js,i);
          e3(m,ks  ,js  ,i) += jcoef*jz_e(m,ks,js  ,i);
          e3(m,ks  ,je+1,i) += jcoef*jz_e(m,ks,je+1,i);
        });
      } else if (pmy_pack->pmesh->two_d) {
        auto e1 = efld.x1e;
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_edge_2d_e1", DevExeSpace(), 0, nmb1, js, je+1,
                is, ie,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e1(m,ks  ,j,i) += jcoef*jx_e(m,ks  ,j,i);
          e1(m,ke+1,j,i) += jcoef*jx_e(m,ke+1,j,i);
        });
        par_for("prtcl_efldsrc_edge_2d_e2", DevExeSpace(), 0, nmb1, js, je,
                is, ie+1,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e2(m,ks  ,j,i) += jcoef*jy_e(m,ks  ,j,i);
          e2(m,ke+1,j,i) += jcoef*jy_e(m,ke+1,j,i);
        });
        par_for("prtcl_efldsrc_edge_2d_e3", DevExeSpace(), 0, nmb1, js, je+1,
                is, ie+1,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e3(m,ks  ,j,i) += jcoef*jz_e(m,ks,j,i);
        });
      } else {
        auto e1 = efld.x1e;
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_edge_3d_e1", DevExeSpace(), 0, nmb1, ks, ke+1,
                js, je+1, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e1(m,k,j,i) += jcoef*jx_e(m,k,j,i);
        });
        par_for("prtcl_efldsrc_edge_3d_e2", DevExeSpace(), 0, nmb1, ks, ke+1,
                js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e2(m,k,j,i) += jcoef*jy_e(m,k,j,i);
        });
        par_for("prtcl_efldsrc_edge_3d_e3", DevExeSpace(), 0, nmb1, ks, ke,
                js, je+1, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e3(m,k,j,i) += jcoef*jz_e(m,k,j,i);
        });
      }
    } else {
      auto mom = ppart->moments;
      if (pmy_pack->pmesh->one_d) {
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_1d", DevExeSpace(), 0, nmb1, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int i) {
          e2(m,ks  ,js  ,i) += jcoef*mom(m, particles::Particles::IMOM_JY, ks, js, i);
          e2(m,ke+1,js  ,i) += jcoef*mom(m, particles::Particles::IMOM_JY, ks, js, i);
          e3(m,ks  ,js  ,i) += jcoef*mom(m, particles::Particles::IMOM_JZ, ks, js, i);
          e3(m,ks  ,je+1,i) += jcoef*mom(m, particles::Particles::IMOM_JZ, ks, js, i);
        });
      } else if (pmy_pack->pmesh->two_d) {
        auto e1 = efld.x1e;
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_2d_e1", DevExeSpace(), 0, nmb1, js, je+1, is, ie,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e1(m,ks  ,j,i) += jcoef*mom(m, particles::Particles::IMOM_JX, ks, j, i);
          e1(m,ke+1,j,i) += jcoef*mom(m, particles::Particles::IMOM_JX, ks, j, i);
        });
        par_for("prtcl_efldsrc_2d_e2", DevExeSpace(), 0, nmb1, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e2(m,ks  ,j,i) += jcoef*mom(m, particles::Particles::IMOM_JY, ks, j, i);
          e2(m,ke+1,j,i) += jcoef*mom(m, particles::Particles::IMOM_JY, ks, j, i);
        });
        par_for("prtcl_efldsrc_2d_e3", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int j, const int i) {
          e3(m,ks  ,j,i) += jcoef*mom(m, particles::Particles::IMOM_JZ, ks, j, i);
        });
      } else {
        auto e1 = efld.x1e;
        auto e2 = efld.x2e;
        auto e3 = efld.x3e;
        par_for("prtcl_efldsrc_3d_e1", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1,
                is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e1(m,k,j,i) += jcoef*mom(m, particles::Particles::IMOM_JX, k, j, i);
        });
        par_for("prtcl_efldsrc_3d_e2", DevExeSpace(), 0, nmb1, ks, ke+1, js, je,
                is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e2(m,k,j,i) += jcoef*mom(m, particles::Particles::IMOM_JY, k, j, i);
        });
        par_for("prtcl_efldsrc_3d_e3", DevExeSpace(), 0, nmb1, ks, ke, js, je+1,
                is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          e3(m,k,j,i) += jcoef*mom(m, particles::Particles::IMOM_JZ, k, j, i);
        });
      }
    }
  }

  // PR2 step II: optional fluid-feedback ordering mode for parity experiments.
  if ((ppart != nullptr) && ppart->couple_moments_to_mhd &&
      (ppart->couple_fluid_feedback_order == CoupledFluidFeedbackOrder::efield_src)) {
    const bool add_mom = ppart->couple_moments_momentum_to_mhd;
    const bool add_eng = ppart->couple_moments_energy_to_mhd;
    const bool use_delta_feedback =
        ((ppart->pusher == ParticlesPusher::boris_lin) ||
         (ppart->pusher == ParticlesPusher::boris_tsc)) &&
        (ppart->pic_feedback_mode == PICFeedbackMode::coupled);
    if ((add_mom || add_eng) && !ppart->UsesExpandingBox()) {
      if (add_eng && (nmhd <= IEN)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Particle energy feedback requires an ideal-MHD energy variable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      Real beta_dt = (pdrive->beta[stage-1])*(pmy_pack->pmesh->dt);
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int is = indcs.is, ie = indcs.ie;
      int js = indcs.js, je = indcs.je;
      int ks = indcs.ks, ke = indcs.ke;
      int nmb1 = pmy_pack->nmb_thispack - 1;
      const Real mom_coef = ppart->couple_moments_momentum_coeff;
      const Real eng_coef = ppart->couple_moments_energy_coeff;
      const bool use_deltaf = ppart->UsesDeltaF();
      Real background_density_scale = 1.0;
      if (ppart->UsesExpandingBox()) {
        const auto geom = particles::PICExpandingBoxGeometryAt(
            ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
            ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3,
            pmy_pack->pmesh->time + pmy_pack->pmesh->dt);
        background_density_scale = geom.inv_a1*geom.inv_a2*geom.inv_a3;
      }
      const Real background_rho =
          background_density_scale*ppart->pic_deltaf_background_rho;
      const Real background_jx = ppart->pic_deltaf_background_jx;
      const Real background_jy = ppart->pic_deltaf_background_jy;
      const Real background_jz = ppart->pic_deltaf_background_jz;
      auto mom = ppart->moments;
      auto bcc = bcc0;
      auto w = w0;
      auto u = u0;

      par_for("prtcl_fluid_feedback_src_efield", DevExeSpace(), 0, nmb1,
              ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (use_deltaf) {
          const Real rho = background_rho +
              mom(m, particles::Particles::IMOM_RHO, k, j, i);
          const Real jx = background_jx +
              mom(m, particles::Particles::IMOM_JX, k, j, i);
          const Real jy = background_jy +
              mom(m, particles::Particles::IMOM_JY, k, j, i);
          const Real jz = background_jz +
              mom(m, particles::Particles::IMOM_JZ, k, j, i);
          const Real bx = bcc(m, IBX, k, j, i);
          const Real by = bcc(m, IBY, k, j, i);
          const Real bz = bcc(m, IBZ, k, j, i);
          const Real cex = -(w(m, IVY, k, j, i)*bz - w(m, IVZ, k, j, i)*by);
          const Real cey = -(w(m, IVZ, k, j, i)*bx - w(m, IVX, k, j, i)*bz);
          const Real cez = -(w(m, IVX, k, j, i)*by - w(m, IVY, k, j, i)*bx);
          if (add_mom) {
            u(m, IM1, k, j, i) -= beta_dt*mom_coef*(rho*cex + jy*bz - jz*by);
            u(m, IM2, k, j, i) -= beta_dt*mom_coef*(rho*cey + jz*bx - jx*bz);
            u(m, IM3, k, j, i) -= beta_dt*mom_coef*(rho*cez + jx*by - jy*bx);
          }
          if (add_eng) {
            u(m, IEN, k, j, i) -= beta_dt*eng_coef*(jx*cex + jy*cey + jz*cez);
          }
        } else if (use_delta_feedback) {
          if (add_mom) {
            u(m, IM1, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPXDT, k, j, i);
            u(m, IM2, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPYDT, k, j, i);
            u(m, IM3, k, j, i) -= beta_dt*mom_coef*
                                  mom(m, particles::Particles::IMOM_DPZDT, k, j, i);
          }
          if (add_eng) {
            u(m, IEN, k, j, i) -= beta_dt*eng_coef*
                                  mom(m, particles::Particles::IMOM_DEDT, k, j, i);
          }
        } else {
          Real jx = mom(m, particles::Particles::IMOM_JX, k, j, i);
          Real jy = mom(m, particles::Particles::IMOM_JY, k, j, i);
          Real jz = mom(m, particles::Particles::IMOM_JZ, k, j, i);

          Real bx = bcc(m, IBX, k, j, i);
          Real by = bcc(m, IBY, k, j, i);
          Real bz = bcc(m, IBZ, k, j, i);

          Real fx = -(jy*bz - jz*by);
          Real fy = -(jz*bx - jx*bz);
          Real fz = -(jx*by - jy*bx);

          if (add_mom) {
            u(m, IM1, k, j, i) += beta_dt*mom_coef*fx;
            u(m, IM2, k, j, i) += beta_dt*mom_coef*fy;
            u(m, IM3, k, j, i) += beta_dt*mom_coef*fz;
          }
          if (add_eng) {
            u(m, IEN, k, j, i) += beta_dt*eng_coef*(jx*bx + jy*by + jz*bz);
          }
        }
      });
    }
  }

  if ((ppart != nullptr) && ppart->UsesExpandingBox()) {
    const auto geom = particles::PICExpandingBoxGeometryAt(
        ppart->pic_expansion_law, ppart->pic_expansion_rate_x1,
        ppart->pic_expansion_rate_x2, ppart->pic_expansion_rate_x3,
        pmy_pack->pmesh->time);
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie;
    const int js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nmb1 = pmy_pack->nmb_thispack - 1;
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    par_for("pic_expanding_box_emf1", DevExeSpace(), 0, nmb1, ks, ke+1,
            js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      e1(m, k, j, i) *= geom.a1;
    });
    par_for("pic_expanding_box_emf2", DevExeSpace(), 0, nmb1, ks, ke+1,
            js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      e2(m, k, j, i) *= geom.a2;
    });
    par_for("pic_expanding_box_emf3", DevExeSpace(), 0, nmb1, ks, ke,
            js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      e3(m, k, j, i) *= geom.a3;
    });
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::SendE
//! \brief Wrapper task list function to pack/send fluxes of magnetic fields
//! (i.e. edge-centered electric field E) at MeshBlock boundaries. This is performed both
//! at MeshBlock boundaries at the same level (to keep magnetic flux in-sync on different
//! MeshBlocks), and at fine/coarse boundaries with SMR/AMR using restricted values of E.

TaskStatus MHD::SendE(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  tstat = pbval_b->PackAndSendFluxFC(efld);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RecvE
//! \brief Wrapper task list function to recv/unpack fluxes of magnetic fields
//! (i.e. edge-centered electric field E) at MeshBlock boundaries

TaskStatus MHD::RecvE(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  tstat = pbval_b->RecvAndUnpackFluxFC(efld);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::SendB_OA
//! \brief Wrapper task list function to pack/send data for orbital advection

TaskStatus MHD::SendB_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = porb_b->PackAndSendFC(b0);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::RecvB_OA
//! \brief Wrapper task list function to recv/unpack data for orbital advection

TaskStatus MHD::RecvB_OA(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    Real qom = (psrc->qshear)*(psrc->omega0);
    tstat = porb_b->RecvAndUnpackFC(b0, recon_method, qom);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::SendB
//! \brief Wrapper task list function to pack/send face-centered magnetic fields

TaskStatus MHD::SendB(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_b->PackAndSendFC(b0, coarse_b0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RecvB
//! \brief Wrapper task list function to recv/unpack face-centered magnetic fields

TaskStatus MHD::RecvB(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_b->RecvAndUnpackFC(b0, coarse_b0);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::SendB_Shr
//! \brief Wrapper task list function to pack/send data for shearing box boundaries

TaskStatus MHD::SendB_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_b->PackAndSendFC(b0, recon_method);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::RecvB_Shr
//! \brief Wrapper task list function to recv/unpack data for shearing box boundaries
//! Orbital remap is performed in this step.

TaskStatus MHD::RecvB_Shr(Driver *pdrive, int stage) {
  TaskStatus tstat = TaskStatus::complete;
  // only execute when (shearing box defined) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_b->RecvAndUnpackFC(b0);
  }
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ApplyPhysicalBCs
//! \brief Wrapper task list function to call funtions that set physical and user BCs

TaskStatus MHD::ApplyPhysicalBCs(Driver *pdrive, int stage) {
  // do not apply BCs if domain is strictly periodic
  if (pmy_pack->pmesh->strictly_periodic) return TaskStatus::complete;

  // physical BCs
  pbval_u->HydroBCs((pmy_pack), (pbval_u->u_in), u0);
  pbval_b->BFieldBCs((pmy_pack), (pbval_b->b_in), b0);

  // user BCs
  if (pmy_pack->pmesh->pgen->user_bcs) {
    (pmy_pack->pmesh->pgen->user_bcs_func)(pmy_pack->pmesh);
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList MHD::Prolongate
//! \brief Wrapper task list function to prolongate conserved (or primitive) variables
//! at fine/coarse bundaries with SMR/AMR

TaskStatus MHD::Prolongate(Driver *pdrive, int stage) {
  if (pmy_pack->pmesh->multilevel) {  // only prolongate with SMR/AMR
    pbval_u->FillCoarseInBndryCC(u0, coarse_u0);
    pbval_b->FillCoarseInBndryFC(b0, coarse_b0);
    if (pmy_pack->pmesh->pmr->prolong_prims) {
      pbval_u->ConsToPrimCoarseBndry(coarse_u0, coarse_b0, coarse_w0);
      pbval_u->ProlongateCC(w0, coarse_w0);
      pbval_b->ProlongateFC(b0, coarse_b0);
      pbval_u->PrimToConsFineBndry(w0, b0, u0);
    } else {
      pbval_u->ProlongateCC(u0, coarse_u0);
      pbval_b->ProlongateFC(b0, coarse_b0);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ConToPrim
//! \brief Wrapper task list function to call ConsToPrim over entire mesh (including gz)

TaskStatus MHD::ConToPrim(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  auto *ppart = pmy_pack->ppart;
  const bool expanding_box = ((ppart != nullptr) && ppart->UsesExpandingBox());
  Real state_time = pmy_pack->pmesh->time;
  if (expanding_box && (stage == pdrive->nexp_stages)) {
    state_time += pmy_pack->pmesh->dt;
  }
  if (expanding_box) RefreshPICExpandingBoxPhysicalB(state_time);
  const DvceFaceFld4D<Real> &bfc = expanding_box ? bphys : b0;
  peos->ConsToPrim(u0, bfc, w0, bcc0, false, 0, n1m1, 0, n2m1, 0, n3m1);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed. Used in
//! TaskList and in Driver::InitBoundaryValuesAndPrimitives()
//! If stage=(last stage):      clears U, B, Flx_U, Flx_B, U_OA, B_OA, U_Shr, BShr
//! If (last stage)>stage>=(0): clears U, B, Flx_U, Flx_B,             U_Shr, B_Shr
//! If stage=(-1):              clears U, B
//! If stage=(-4):              clears                                 U_Shr, B_Shr

TaskStatus MHD::ClearSend(Driver *pdrive, int stage) {
  TaskStatus tstat;
  if ((stage >= 0) || (stage == -1)) {
    // check sends of U complete
    TaskStatus tstat = pbval_u->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
    // check sends of B complete
    tstat = pbval_b->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with SMR/AMR check sends for fluxes of U complete.  Always check sends of E complete
  // do not check flux send for ICs (stage < 0)
  if (stage >= 0) {
    // with SMR/AMR check sends of restricted fluxes of U complete
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->ClearFluxSend();
      if (tstat != TaskStatus::complete) return tstat;
    }
    // check sends of restricted fluxes of B complete even for uniform grids
    tstat = pbval_b->ClearFluxSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with orbital advection check sends for U and B complete
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = porb_u->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
    tstat = porb_b->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with shearing box boundaries check sends of U and B complete
  // only execute when (shearing box defined) AND (stage>=0 or -4) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && ((stage >= 0) || (stage == -4)) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_u->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
    tstat = psbox_b->ClearSend();
    if (tstat != TaskStatus::complete) return tstat;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed. Used in
//! TaskList and in Driver::InitBoundaryValuesAndPrimitives()
//! If stage=(last stage):      clears U, B, Flx_U, Flx_B, U_OA, B_OA, U_Shr, BShr
//! If (last stage)>stage>=(0): clears U, B, Flx_U, Flx_B,             U_Shr, B_Shr
//! If stage=(-1):              clears U, B
//! If stage=(-4):              clears                                 U_Shr, B_Shr

TaskStatus MHD::ClearRecv(Driver *pdrive, int stage) {
  TaskStatus tstat;
  if ((stage >= 0) || (stage == -1)) {
    // check receives of U complete
    tstat = pbval_u->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
    // check receives of B complete
    tstat = pbval_b->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with SMR/AMR check recvs for fluxes of U complete.  Always check recvs of E complete
  // do not check flux receives when stage < 0 (i.e. ICs)
  if (stage >= 0) {
    // with SMR/AMR check receives of restricted fluxes of U complete
    if (pmy_pack->pmesh->multilevel) {
      tstat = pbval_u->ClearFluxRecv();
      if (tstat != TaskStatus::complete) return tstat;
    }
    // with SMR/AMR check receives of restricted fluxes of B complete
    tstat = pbval_b->ClearFluxRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with orbital advection check receives of U and B are complete
  // only execute when (shearing box defined) AND (last stage) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && (stage == pdrive->nexp_stages) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = porb_u->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
    tstat = porb_b->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  // with shearing box boundaries check receives of U and B complete
  // only execute when (shearing box defined) AND (stage>=0 or -4) AND (3D OR 2d_r_phi)
  if ((psrc->shearing_box) && ((stage >= 0) || (stage == -4)) &&
      (pmy_pack->pmesh->three_d || psrc->shearing_box_r_phi)) {
    tstat = psbox_u->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
    tstat = psbox_b->ClearRecv();
    if (tstat != TaskStatus::complete) return tstat;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::RestrictB
//! \brief Wrapper function that restricts face-centered variables (magnetic field)

TaskStatus MHD::RestrictB(Driver *pdrive, int stage) {
  // Only execute Mesh function with SMR/AMR
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictFC(b0, coarse_b0);
  }
  return TaskStatus::complete;
}

} // namespace mhd
