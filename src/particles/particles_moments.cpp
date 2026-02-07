//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_moments.cpp
//! \brief task wrappers and kernels for particle moment deposition/communication

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"

namespace particles {

namespace {
// PR2 coupling inserts wrappers into stagen; keep PR1 default in before_timeintegrator.
KOKKOS_INLINE_FUNCTION
bool RunMomentWrappersAtStage(const bool couple_to_mhd, const int stage) {
  if (couple_to_mhd) {
    return (stage == 1);
  }
  return (stage == 0);
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::SaveOldPositions()
//! \brief Store particle positions before pushing.

TaskStatus Particles::SaveOldPositions(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }

  int npart = nprtcl_thispack;
  if (npart <= 0) return TaskStatus::complete;

  if (x1_old.extent_int(0) < npart) {
    Kokkos::resize(x1_old, npart);
    Kokkos::resize(x2_old, npart);
    Kokkos::resize(x3_old, npart);
  }

  auto &pr = prtcl_rdata;
  auto &x1 = x1_old;
  auto &x2 = x2_old;
  auto &x3 = x3_old;

  par_for("save_old_positions", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    x1(p) = pr(IPX,p);
    x2(p) = pr(IPY,p);
    x3(p) = pr(IPZ,p);
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ZeroMoments()
//! \brief Zero particle moment arrays before local deposition.

TaskStatus Particles::ZeroMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }

  Kokkos::deep_copy(moments, static_cast<Real>(0.0));
  if (coarse_moments.size() > 0) {
    Kokkos::deep_copy(coarse_moments, static_cast<Real>(0.0));
  }
  if (couple_moments_to_mhd &&
      couple_j_to_efield_representation == CoupledCurrentRepresentation::edge_staggered &&
      (j_edge_x1e.size() > 0)) {
    Kokkos::deep_copy(j_edge_x1e, static_cast<Real>(0.0));
    Kokkos::deep_copy(j_edge_x2e, static_cast<Real>(0.0));
    Kokkos::deep_copy(j_edge_x3e, static_cast<Real>(0.0));
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::InitRecvMoments()
//! \brief Post non-blocking receives for deposited moments.

TaskStatus Particles::InitRecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  return pbval_mom->InitRecv(NMOM);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::DepositMoments()
//! \brief Deposit rho/J moments to mesh cells.

TaskStatus Particles::DepositMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }

  int npart = nprtcl_thispack;
  if (npart <= 0) return TaskStatus::complete;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ng = indcs.ng;
  const int i_min = is - ng;
  const int i_max = ie + ng;
  const int j_min = js - ng;
  const int j_max = je + ng;
  const int k_min = ks - ng;
  const int k_max = ke + ng;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int nmb = pmy_pack->nmb_thispack;
  const int gids = pmy_pack->gids;
  const int nspecies_local = nspecies;
  const Real qscale = deposit_qscale;

  auto &size = pmy_pack->pmb->mb_size;
  auto &pr = prtcl_rdata;
  auto &pi = prtcl_idata;
  auto &qspecies = species_charge;
  auto &mom = moments;

  par_for("deposit_moments", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) return;

    int sp = pi(PSP,p);
    if (sp < 0 || sp >= nspecies_local) return;

    Real q_macro = qscale*qspecies(sp);

    Real x = pr(IPX,p);
    Real y = pr(IPY,p);
    Real z = pr(IPZ,p);
    Real x1min = size.d_view(m).x1min;
    Real x2min = size.d_view(m).x2min;
    Real x3min = size.d_view(m).x3min;
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;

    int ip = static_cast<int>((x - x1min)/dx1) + is;
    int jp = static_cast<int>((y - x2min)/dx2) + js;
    if (ip < i_min) ip = i_min;
    if (ip > i_max) ip = i_max;
    if (jp < j_min) jp = j_min;
    if (jp > j_max) jp = j_max;

    int kp = ks;
    if (three_d) {
      kp = static_cast<int>((z - x3min)/dx3) + ks;
      if (kp < k_min) kp = k_min;
      if (kp > k_max) kp = k_max;
    }

    Kokkos::atomic_add(&mom(m, IMOM_RHO, kp, jp, ip), q_macro);
    Kokkos::atomic_add(&mom(m, IMOM_JX,  kp, jp, ip), q_macro*pr(IPVX,p));
    Kokkos::atomic_add(&mom(m, IMOM_JY,  kp, jp, ip), q_macro*pr(IPVY,p));
    Kokkos::atomic_add(&mom(m, IMOM_JZ,  kp, jp, ip), q_macro*pr(IPVZ,p));
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::RestrictMoments()
//! \brief Restrict deposited moments to coarse storage in multilevel runs.

TaskStatus Particles::RestrictMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  if (!(pmy_pack->pmesh->multilevel)) return TaskStatus::complete;
  if (coarse_moments.size() == 0) return TaskStatus::complete;

  pmy_pack->pmesh->pmr->RestrictCC(moments, coarse_moments);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::SendMoments()
//! \brief Pack and send deposited moments to neighbors.

TaskStatus Particles::SendMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  return pbval_mom->PackAndSendCC(moments, coarse_moments);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::RecvMoments()
//! \brief Receive and unpack deposited moments with additive accumulation.

TaskStatus Particles::RecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  return pbval_mom->RecvAndUnpackCC(moments, coarse_moments, CCRecvOp::accumulate);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearRecvMoments()
//! \brief Wait for moment receive completion.

TaskStatus Particles::ClearRecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  return pbval_mom->ClearRecv();
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearSendMoments()
//! \brief Wait for moment send completion.

TaskStatus Particles::ClearSendMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  return pbval_mom->ClearSend();
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ApplyMomentPhysicalBCs()
//! \brief Apply physical BCs to deposited rho/J moments at mesh boundaries.

TaskStatus Particles::ApplyMomentPhysicalBCs(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }

  MeshBoundaryValues::HydroBCs(pmy_pack, pbval_mom->u_in, moments);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ProlongateMoments()
//! \brief Fill coarse boundary state and prolongate to fine moment ghosts with AMR/SMR.

TaskStatus Particles::ProlongateMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  if (!(pmy_pack->pmesh->multilevel)) return TaskStatus::complete;
  if (coarse_moments.size() == 0) return TaskStatus::complete;

  pbval_mom->FillCoarseInBndryCC(moments, coarse_moments);
  pbval_mom->ProlongateCC(moments, coarse_moments);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ConvertCoupledCurrentRepresentation()
//! \brief Convert deposited cell-centered J to edge-centered J for coupled E updates.

TaskStatus Particles::ConvertCoupledCurrentRepresentation(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!couple_moments_to_mhd) return TaskStatus::complete;
  if (!RunMomentWrappersAtStage(couple_moments_to_mhd, stage)) {
    return TaskStatus::complete;
  }
  if (couple_j_to_efield_representation != CoupledCurrentRepresentation::edge_staggered) {
    return TaskStatus::complete;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto mom = moments;
  auto jx_e = j_edge_x1e;
  auto jy_e = j_edge_x2e;
  auto jz_e = j_edge_x3e;

  if (pmy_pack->pmesh->one_d) {
    par_for("convert_jy_cc_to_edge_1d", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int i) {
      Real jedge = mom(m, IMOM_JY, ks, js, i);
      jy_e(m,ks  ,js,i) = jedge;
      jy_e(m,ke+1,js,i) = jedge;
    });

    par_for("convert_jz_cc_to_edge_1d", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int i) {
      Real jedge = mom(m, IMOM_JZ, ks, js, i);
      jz_e(m,ks,js  ,i) = jedge;
      jz_e(m,ks,je+1,i) = jedge;
    });
    return TaskStatus::complete;
  }

  if (pmy_pack->pmesh->two_d) {
    par_for("convert_jx_cc_to_edge_2d", DevExeSpace(), 0, nmb1, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      Real jedge = 0.5*(mom(m, IMOM_JX, ks, j-1, i) + mom(m, IMOM_JX, ks, j, i));
      jx_e(m,ks  ,j,i) = jedge;
      jx_e(m,ke+1,j,i) = jedge;
    });

    par_for("convert_jy_cc_to_edge_2d", DevExeSpace(), 0, nmb1, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      Real jedge = 0.5*(mom(m, IMOM_JY, ks, j, i-1) + mom(m, IMOM_JY, ks, j, i));
      jy_e(m,ks  ,j,i) = jedge;
      jy_e(m,ke+1,j,i) = jedge;
    });

    par_for("convert_jz_cc_to_edge_2d", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      jz_e(m,ks,j,i) = 0.25*(mom(m, IMOM_JZ, ks, j-1, i-1) +
                             mom(m, IMOM_JZ, ks, j-1, i  ) +
                             mom(m, IMOM_JZ, ks, j,   i-1) +
                             mom(m, IMOM_JZ, ks, j,   i  ));
    });
    return TaskStatus::complete;
  }

  par_for("convert_jx_cc_to_edge_3d", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    jx_e(m,k,j,i) = 0.25*(mom(m, IMOM_JX, k-1, j-1, i) +
                          mom(m, IMOM_JX, k-1, j,   i) +
                          mom(m, IMOM_JX, k,   j-1, i) +
                          mom(m, IMOM_JX, k,   j,   i));
  });

  par_for("convert_jy_cc_to_edge_3d", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    jy_e(m,k,j,i) = 0.25*(mom(m, IMOM_JY, k-1, j, i-1) +
                          mom(m, IMOM_JY, k-1, j, i  ) +
                          mom(m, IMOM_JY, k,   j, i-1) +
                          mom(m, IMOM_JY, k,   j, i  ));
  });

  par_for("convert_jz_cc_to_edge_3d", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    jz_e(m,k,j,i) = 0.25*(mom(m, IMOM_JZ, k, j-1, i-1) +
                          mom(m, IMOM_JZ, k, j-1, i  ) +
                          mom(m, IMOM_JZ, k, j,   i-1) +
                          mom(m, IMOM_JZ, k, j,   i  ));
  });

  return TaskStatus::complete;
}

} // namespace particles
