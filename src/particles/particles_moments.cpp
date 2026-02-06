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
//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::SaveOldPositions()
//! \brief Store particle positions before pushing.

TaskStatus Particles::SaveOldPositions(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!deposit_moments) return TaskStatus::complete;

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
  (void)stage;
  if (!deposit_moments) return TaskStatus::complete;

  Kokkos::deep_copy(moments, static_cast<Real>(0.0));
  if (coarse_moments.size() > 0) {
    Kokkos::deep_copy(coarse_moments, static_cast<Real>(0.0));
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::InitRecvMoments()
//! \brief Post non-blocking receives for deposited moments.

TaskStatus Particles::InitRecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  return pbval_mom->InitRecv(NMOM);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::DepositMoments()
//! \brief Deposit rho/J moments to mesh cells.

TaskStatus Particles::DepositMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!deposit_moments) return TaskStatus::complete;

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
//! \fn TaskStatus Particles::SendMoments()
//! \brief Pack and send deposited moments to neighbors.

TaskStatus Particles::SendMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  return pbval_mom->PackAndSendCC(moments, coarse_moments);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::RecvMoments()
//! \brief Receive and unpack deposited moments with additive accumulation.

TaskStatus Particles::RecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  return pbval_mom->RecvAndUnpackCC(moments, coarse_moments, CCRecvOp::accumulate);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearRecvMoments()
//! \brief Wait for moment receive completion.

TaskStatus Particles::ClearRecvMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  return pbval_mom->ClearRecv();
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearSendMoments()
//! \brief Wait for moment send completion.

TaskStatus Particles::ClearSendMoments(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!(deposit_moments) || pbval_mom == nullptr) return TaskStatus::complete;
  return pbval_mom->ClearSend();
}

} // namespace particles
