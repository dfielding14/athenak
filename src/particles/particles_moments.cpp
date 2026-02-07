//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_moments.cpp
//! \brief task wrappers and kernels for particle moment deposition/communication

#include <cmath>

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

KOKKOS_INLINE_FUNCTION
bool UseDirectEdgeCurrentDeposit(const bool deposit_moments,
                                 const bool couple_to_mhd,
                                 const CoupledCurrentRepresentation repr,
                                 const CoupledCurrentDepositionMode mode) {
  return (deposit_moments &&
          couple_to_mhd &&
          (repr == CoupledCurrentRepresentation::edge_staggered) &&
          (mode == CoupledCurrentDepositionMode::direct_staggered));
}

KOKKOS_INLINE_FUNCTION
int ClampCellIndex(const int i, const int imin, const int imax) {
  if (i < imin) return imin;
  if (i > imax) return imax;
  return i;
}

KOKKOS_INLINE_FUNCTION
void PosToCellAndFrac(const Real x, const Real xmin, const Real dx, const int is,
                      const int imin, const int imax, int &icell, Real &frac) {
  Real xrel = (x - xmin)/dx;
  Real ifloor = floor(xrel);
  int iloc = static_cast<int>(ifloor) + is;
  frac = xrel - ifloor;
  if (frac < static_cast<Real>(0.0)) frac = static_cast<Real>(0.0);
  if (frac > static_cast<Real>(1.0)) frac = static_cast<Real>(1.0);
  icell = ClampCellIndex(iloc, imin, imax);
}

inline DvceEdgeFld4D<Real> MakeEdgeCurrentAlias(const DvceArray4D<Real> &x1e,
                                                const DvceArray4D<Real> &x2e,
                                                const DvceArray4D<Real> &x3e) {
  DvceEdgeFld4D<Real> edge_curr("prtcl_jedge_alias", 0, 0, 0, 0);
  edge_curr.x1e = x1e;
  edge_curr.x2e = x2e;
  edge_curr.x3e = x3e;
  return edge_curr;
}

inline bool RunEdgeCurrentWrappersAtStage(const bool deposit_moments,
                                          const bool couple_to_mhd,
                                          const CoupledCurrentRepresentation repr,
                                          const CoupledCurrentDepositionMode mode,
                                          const int stage) {
  if (!UseDirectEdgeCurrentDeposit(deposit_moments, couple_to_mhd, repr, mode)) {
    return false;
  }
  return RunMomentWrappersAtStage(couple_to_mhd, stage);
}

// PR4b scaffolding: trajectory-based direct deposition requires pre-push old
// positions from before_timeintegrator (stage 0), not stagen (stage 1).
KOKKOS_INLINE_FUNCTION
bool RunSaveOldPositionsAtStage(const bool couple_to_mhd,
                                const CoupledCurrentDepositionMode deposit_mode,
                                const int stage) {
  if (!couple_to_mhd) {
    return (stage == 0);
  }
  if (deposit_mode == CoupledCurrentDepositionMode::direct_staggered) {
    return (stage == 0);
  }
  return (stage == 1);
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::SaveOldPositions()
//! \brief Store particle positions before pushing.

TaskStatus Particles::SaveOldPositions(Driver *pdriver, int stage) {
  (void)pdriver;
  if (!deposit_moments) return TaskStatus::complete;
  if (!RunSaveOldPositionsAtStage(couple_moments_to_mhd, couple_j_deposition_mode,
                                  stage)) {
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

  if (!UseDirectEdgeCurrentDeposit(deposit_moments, couple_moments_to_mhd,
                                   couple_j_to_efield_representation,
                                   couple_j_deposition_mode)) {
    return TaskStatus::complete;
  }

  const Real inv_dt = static_cast<Real>(1.0) / pmy_pack->pmesh->dt;
  auto &x1prev = x1_old;
  auto &x2prev = x2_old;
  auto &x3prev = x3_old;
  auto &jx_e = j_edge_x1e;
  auto &jy_e = j_edge_x2e;
  auto &jz_e = j_edge_x3e;

  par_for("deposit_direct_edge_currents", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) return;

    int sp = pi(PSP,p);
    if (sp < 0 || sp >= nspecies_local) return;

    Real q_macro = qscale*qspecies(sp);
    Real coeff = q_macro*inv_dt;

    Real x1min = size.d_view(m).x1min;
    Real x2min = size.d_view(m).x2min;
    Real x3min = size.d_view(m).x3min;
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;

    int i0, i1, j0, j1, k0, k1;
    Real dxo, dxn, dyo, dyn, dzo, dzn;
    PosToCellAndFrac(x1prev(p), x1min, dx1, is, i_min, i_max, i0, dxo);
    PosToCellAndFrac(pr(IPX,p), x1min, dx1, is, i_min, i_max, i1, dxn);
    PosToCellAndFrac(x2prev(p), x2min, dx2, js, j_min, j_max, j0, dyo);
    PosToCellAndFrac(pr(IPY,p), x2min, dx2, js, j_min, j_max, j1, dyn);
    if (three_d) {
      PosToCellAndFrac(x3prev(p), x3min, dx3, ks, k_min, k_max, k0, dzo);
      PosToCellAndFrac(pr(IPZ,p), x3min, dx3, ks, k_min, k_max, k1, dzn);
    } else {
      k0 = ks;
      k1 = ks;
      dzo = static_cast<Real>(0.0);
      dzn = static_cast<Real>(0.0);
    }

    Real vy = pr(IPVY,p);
    Real vz = pr(IPVZ,p);

    Real dxp1 = ((i1 == i0) ? static_cast<Real>(0.5)*(dxn + dxo)
                            : static_cast<Real>(0.0));
    Real wx1_1 = static_cast<Real>(0.5)*(dxp1 + dxo +
                                         static_cast<Real>(i1 > i0));
    Real wx1_2 = static_cast<Real>(0.5)*(dxn + dxp1 +
      static_cast<Real>(static_cast<int>(i1 > i0) + i0 - i1));
    Real fx1_1 = (static_cast<Real>(i1 > i0) + dxp1 - dxo)*coeff;
    Real fx1_2 = (static_cast<Real>(i1 - i0 - static_cast<int>(i1 > i0)) +
                  dxn - dxp1)*coeff;

    if (pmy_pack->pmesh->one_d) {
      Real fx2_1 = static_cast<Real>(0.5)*vy*q_macro;
      Real fx2_2 = static_cast<Real>(0.5)*vy*q_macro;
      Real fx3_1 = static_cast<Real>(0.5)*vz*q_macro;
      Real fx3_2 = static_cast<Real>(0.5)*vz*q_macro;

      if (i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jy_e(m, ks, js, i0), fx2_1*(static_cast<Real>(1.0) - wx1_1));
        Kokkos::atomic_add(&jy_e(m, ke+1, js, i0),
                           fx2_1*(static_cast<Real>(1.0) - wx1_1));
        Kokkos::atomic_add(&jz_e(m, ks, js, i0), fx3_1*(static_cast<Real>(1.0) - wx1_1));
        Kokkos::atomic_add(&jz_e(m, ks, je+1, i0),
                           fx3_1*(static_cast<Real>(1.0) - wx1_1));
      }
      if ((i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jy_e(m, ks, js, i0+1), fx2_1*wx1_1);
        Kokkos::atomic_add(&jy_e(m, ke+1, js, i0+1), fx2_1*wx1_1);
        Kokkos::atomic_add(&jz_e(m, ks, js, i0+1), fx3_1*wx1_1);
        Kokkos::atomic_add(&jz_e(m, ks, je+1, i0+1), fx3_1*wx1_1);
      }
      if (i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jy_e(m, ks, js, i1), fx2_2*(static_cast<Real>(1.0) - wx1_2));
        Kokkos::atomic_add(&jy_e(m, ke+1, js, i1),
                           fx2_2*(static_cast<Real>(1.0) - wx1_2));
        Kokkos::atomic_add(&jz_e(m, ks, js, i1), fx3_2*(static_cast<Real>(1.0) - wx1_2));
        Kokkos::atomic_add(&jz_e(m, ks, je+1, i1),
                           fx3_2*(static_cast<Real>(1.0) - wx1_2));
      }
      if ((i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jy_e(m, ks, js, i1+1), fx2_2*wx1_2);
        Kokkos::atomic_add(&jy_e(m, ke+1, js, i1+1), fx2_2*wx1_2);
        Kokkos::atomic_add(&jz_e(m, ks, js, i1+1), fx3_2*wx1_2);
        Kokkos::atomic_add(&jz_e(m, ks, je+1, i1+1), fx3_2*wx1_2);
      }
      return;
    }

    Real dxp2 = ((j1 == j0) ? static_cast<Real>(0.5)*(dyn + dyo)
                            : static_cast<Real>(0.0));
    Real wx2_1 = static_cast<Real>(0.5)*(dxp2 + dyo +
                                         static_cast<Real>(j1 > j0));
    Real wx2_2 = static_cast<Real>(0.5)*(dyn + dxp2 +
      static_cast<Real>(static_cast<int>(j1 > j0) + j0 - j1));
    Real fx2_1 = (static_cast<Real>(j1 > j0) + dxp2 - dyo)*coeff;
    Real fx2_2 = (static_cast<Real>(j1 - j0 - static_cast<int>(j1 > j0)) +
                  dyn - dxp2)*coeff;

    if (pmy_pack->pmesh->two_d) {
      Real fx3_1 = static_cast<Real>(0.5)*vz*q_macro;
      Real fx3_2 = static_cast<Real>(0.5)*vz*q_macro;
      int kf = ks;
      int kf2 = ke + 1;

      if (j0 >= j_min && j0 <= j_max && i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jx_e(m, kf, j0, i0), fx1_1*(static_cast<Real>(1.0) - wx2_1));
        Kokkos::atomic_add(&jx_e(m, kf2, j0, i0), fx1_1*(static_cast<Real>(1.0) - wx2_1));
      }
      if ((j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
          i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jx_e(m, kf, j0+1, i0), fx1_1*wx2_1);
        Kokkos::atomic_add(&jx_e(m, kf2, j0+1, i0), fx1_1*wx2_1);
      }
      if (j1 >= j_min && j1 <= j_max && i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jx_e(m, kf, j1, i1), fx1_2*(static_cast<Real>(1.0) - wx2_2));
        Kokkos::atomic_add(&jx_e(m, kf2, j1, i1), fx1_2*(static_cast<Real>(1.0) - wx2_2));
      }
      if ((j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
          i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jx_e(m, kf, j1+1, i1), fx1_2*wx2_2);
        Kokkos::atomic_add(&jx_e(m, kf2, j1+1, i1), fx1_2*wx2_2);
      }

      if (j0 >= j_min && j0 <= j_max && i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jy_e(m, kf, j0, i0), fx2_1*(static_cast<Real>(1.0) - wx1_1));
        Kokkos::atomic_add(&jy_e(m, kf2, j0, i0), fx2_1*(static_cast<Real>(1.0) - wx1_1));
      }
      if (j0 >= j_min && j0 <= j_max && (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jy_e(m, kf, j0, i0+1), fx2_1*wx1_1);
        Kokkos::atomic_add(&jy_e(m, kf2, j0, i0+1), fx2_1*wx1_1);
      }
      if (j1 >= j_min && j1 <= j_max && i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jy_e(m, kf, j1, i1), fx2_2*(static_cast<Real>(1.0) - wx1_2));
        Kokkos::atomic_add(&jy_e(m, kf2, j1, i1), fx2_2*(static_cast<Real>(1.0) - wx1_2));
      }
      if (j1 >= j_min && j1 <= j_max && (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jy_e(m, kf, j1, i1+1), fx2_2*wx1_2);
        Kokkos::atomic_add(&jy_e(m, kf2, j1, i1+1), fx2_2*wx1_2);
      }

      if (j0 >= j_min && j0 <= j_max && i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jz_e(m, ks, j0, i0),
                           fx3_1*(static_cast<Real>(1.0) - wx1_1)*
                           (static_cast<Real>(1.0) - wx2_1));
      }
      if (j0 >= j_min && j0 <= j_max && (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jz_e(m, ks, j0, i0+1),
                           fx3_1*wx1_1*(static_cast<Real>(1.0) - wx2_1));
      }
      if ((j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
          i0 >= i_min && i0 <= i_max) {
        Kokkos::atomic_add(&jz_e(m, ks, j0+1, i0),
                           fx3_1*(static_cast<Real>(1.0) - wx1_1)*wx2_1);
      }
      if ((j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
          (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jz_e(m, ks, j0+1, i0+1), fx3_1*wx1_1*wx2_1);
      }
      if (j1 >= j_min && j1 <= j_max && i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jz_e(m, ks, j1, i1),
                           fx3_2*(static_cast<Real>(1.0) - wx1_2)*
                           (static_cast<Real>(1.0) - wx2_2));
      }
      if (j1 >= j_min && j1 <= j_max && (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jz_e(m, ks, j1, i1+1),
                           fx3_2*wx1_2*(static_cast<Real>(1.0) - wx2_2));
      }
      if ((j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
          i1 >= i_min && i1 <= i_max) {
        Kokkos::atomic_add(&jz_e(m, ks, j1+1, i1),
                           fx3_2*(static_cast<Real>(1.0) - wx1_2)*wx2_2);
      }
      if ((j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
          (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
        Kokkos::atomic_add(&jz_e(m, ks, j1+1, i1+1), fx3_2*wx1_2*wx2_2);
      }
      return;
    }

    Real dxp3 = ((k1 == k0) ? static_cast<Real>(0.5)*(dzn + dzo)
                            : static_cast<Real>(0.0));
    Real wx3_1 = static_cast<Real>(0.5)*(dxp3 + dzo +
                                         static_cast<Real>(k1 > k0));
    Real wx3_2 = static_cast<Real>(0.5)*(dzn + dxp3 +
      static_cast<Real>(static_cast<int>(k1 > k0) + k0 - k1));
    Real fx3_1 = (static_cast<Real>(k1 > k0) + dxp3 - dzo)*coeff;
    Real fx3_2 = (static_cast<Real>(k1 - k0 - static_cast<int>(k1 > k0)) +
                  dzn - dxp3)*coeff;

    if (k0 >= k_min && k0 <= (k_max + 1) &&
        j0 >= j_min && j0 <= (j_max + 1) &&
        i0 >= i_min && i0 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k0, j0, i0), fx1_1*(static_cast<Real>(1.0) - wx2_1)*
                         (static_cast<Real>(1.0) - wx3_1));
    }
    if (k0 >= k_min && k0 <= (k_max + 1) &&
        (j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
        i0 >= i_min && i0 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k0, j0+1, i0), fx1_1*wx2_1*
                         (static_cast<Real>(1.0) - wx3_1));
    }
    if ((k0 + 1) >= k_min && (k0 + 1) <= (k_max + 1) &&
        j0 >= j_min && j0 <= (j_max + 1) &&
        i0 >= i_min && i0 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k0+1, j0, i0), fx1_1*
                         (static_cast<Real>(1.0) - wx2_1)*wx3_1);
    }
    if ((k0 + 1) >= k_min && (k0 + 1) <= (k_max + 1) &&
        (j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
        i0 >= i_min && i0 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k0+1, j0+1, i0), fx1_1*wx2_1*wx3_1);
    }
    if (k1 >= k_min && k1 <= (k_max + 1) &&
        j1 >= j_min && j1 <= (j_max + 1) &&
        i1 >= i_min && i1 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k1, j1, i1), fx1_2*(static_cast<Real>(1.0) - wx2_2)*
                         (static_cast<Real>(1.0) - wx3_2));
    }
    if (k1 >= k_min && k1 <= (k_max + 1) &&
        (j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
        i1 >= i_min && i1 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k1, j1+1, i1), fx1_2*wx2_2*
                         (static_cast<Real>(1.0) - wx3_2));
    }
    if ((k1 + 1) >= k_min && (k1 + 1) <= (k_max + 1) &&
        j1 >= j_min && j1 <= (j_max + 1) &&
        i1 >= i_min && i1 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k1+1, j1, i1), fx1_2*
                         (static_cast<Real>(1.0) - wx2_2)*wx3_2);
    }
    if ((k1 + 1) >= k_min && (k1 + 1) <= (k_max + 1) &&
        (j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
        i1 >= i_min && i1 <= i_max) {
      Kokkos::atomic_add(&jx_e(m, k1+1, j1+1, i1), fx1_2*wx2_2*wx3_2);
    }

    if (k0 >= k_min && k0 <= (k_max + 1) &&
        j0 >= j_min && j0 <= j_max &&
        i0 >= i_min && i0 <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k0, j0, i0), fx2_1*(static_cast<Real>(1.0) - wx1_1)*
                         (static_cast<Real>(1.0) - wx3_1));
    }
    if (k0 >= k_min && k0 <= (k_max + 1) &&
        j0 >= j_min && j0 <= j_max &&
        (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k0, j0, i0+1), fx2_1*wx1_1*
                         (static_cast<Real>(1.0) - wx3_1));
    }
    if ((k0 + 1) >= k_min && (k0 + 1) <= (k_max + 1) &&
        j0 >= j_min && j0 <= j_max &&
        i0 >= i_min && i0 <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k0+1, j0, i0), fx2_1*
                         (static_cast<Real>(1.0) - wx1_1)*wx3_1);
    }
    if ((k0 + 1) >= k_min && (k0 + 1) <= (k_max + 1) &&
        j0 >= j_min && j0 <= j_max &&
        (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k0+1, j0, i0+1), fx2_1*wx1_1*wx3_1);
    }
    if (k1 >= k_min && k1 <= (k_max + 1) &&
        j1 >= j_min && j1 <= j_max &&
        i1 >= i_min && i1 <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k1, j1, i1), fx2_2*(static_cast<Real>(1.0) - wx1_2)*
                         (static_cast<Real>(1.0) - wx3_2));
    }
    if (k1 >= k_min && k1 <= (k_max + 1) &&
        j1 >= j_min && j1 <= j_max &&
        (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k1, j1, i1+1), fx2_2*wx1_2*
                         (static_cast<Real>(1.0) - wx3_2));
    }
    if ((k1 + 1) >= k_min && (k1 + 1) <= (k_max + 1) &&
        j1 >= j_min && j1 <= j_max &&
        i1 >= i_min && i1 <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k1+1, j1, i1), fx2_2*
                         (static_cast<Real>(1.0) - wx1_2)*wx3_2);
    }
    if ((k1 + 1) >= k_min && (k1 + 1) <= (k_max + 1) &&
        j1 >= j_min && j1 <= j_max &&
        (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jy_e(m, k1+1, j1, i1+1), fx2_2*wx1_2*wx3_2);
    }

    if (k0 >= k_min && k0 <= k_max &&
        j0 >= j_min && j0 <= (j_max + 1) &&
        i0 >= i_min && i0 <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k0, j0, i0), fx3_1*(static_cast<Real>(1.0) - wx1_1)*
                         (static_cast<Real>(1.0) - wx2_1));
    }
    if (k0 >= k_min && k0 <= k_max &&
        j0 >= j_min && j0 <= (j_max + 1) &&
        (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k0, j0, i0+1), fx3_1*wx1_1*
                         (static_cast<Real>(1.0) - wx2_1));
    }
    if (k0 >= k_min && k0 <= k_max &&
        (j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
        i0 >= i_min && i0 <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k0, j0+1, i0), fx3_1*
                         (static_cast<Real>(1.0) - wx1_1)*wx2_1);
    }
    if (k0 >= k_min && k0 <= k_max &&
        (j0 + 1) >= j_min && (j0 + 1) <= (j_max + 1) &&
        (i0 + 1) >= i_min && (i0 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k0, j0+1, i0+1), fx3_1*wx1_1*wx2_1);
    }
    if (k1 >= k_min && k1 <= k_max &&
        j1 >= j_min && j1 <= (j_max + 1) &&
        i1 >= i_min && i1 <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k1, j1, i1), fx3_2*(static_cast<Real>(1.0) - wx1_2)*
                         (static_cast<Real>(1.0) - wx2_2));
    }
    if (k1 >= k_min && k1 <= k_max &&
        j1 >= j_min && j1 <= (j_max + 1) &&
        (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k1, j1, i1+1), fx3_2*wx1_2*
                         (static_cast<Real>(1.0) - wx2_2));
    }
    if (k1 >= k_min && k1 <= k_max &&
        (j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
        i1 >= i_min && i1 <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k1, j1+1, i1), fx3_2*
                         (static_cast<Real>(1.0) - wx1_2)*wx2_2);
    }
    if (k1 >= k_min && k1 <= k_max &&
        (j1 + 1) >= j_min && (j1 + 1) <= (j_max + 1) &&
        (i1 + 1) >= i_min && (i1 + 1) <= (i_max + 1)) {
      Kokkos::atomic_add(&jz_e(m, k1, j1+1, i1+1), fx3_2*wx1_2*wx2_2);
    }
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
//! \fn TaskStatus Particles::InitRecvEdgeCurrents()
//! \brief Post receives for direct edge-current synchronization.

TaskStatus Particles::InitRecvEdgeCurrents(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }
  return pbval_jedge->InitFluxRecv(3);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::SendEdgeCurrents()
//! \brief Pack/send direct edge-currents to neighboring MeshBlocks.

TaskStatus Particles::SendEdgeCurrents(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }
  auto jedge = MakeEdgeCurrentAlias(j_edge_x1e, j_edge_x2e, j_edge_x3e);
  return pbval_jedge->PackAndSendFluxFC(jedge);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::RecvEdgeCurrents()
//! \brief Receive and accumulate direct edge-currents from neighboring MeshBlocks.

TaskStatus Particles::RecvEdgeCurrents(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }
  TaskStatus tstat = pbval_jedge->ClearFluxRecv();
  if (tstat != TaskStatus::complete) return tstat;

  if (pmy_pack->nmb_thispack <= 0) return TaskStatus::complete;
  auto jedge = MakeEdgeCurrentAlias(j_edge_x1e, j_edge_x2e, j_edge_x3e);
  DvceArray2D<int> nflx("prtcl_jedge_nflx", pmy_pack->nmb_thispack, 48);
  par_for("init_prtcl_jedge_nflx", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, 47,
  KOKKOS_LAMBDA(const int m, const int n) {
    nflx(m,n) = 1;
  });
  pbval_jedge->SumBoundaryFluxes(jedge, true, nflx);
  if (pmy_pack->pmesh->multilevel) {
    pbval_jedge->SumBoundaryFluxes(jedge, false, nflx);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearRecvEdgeCurrents()
//! \brief Verify completion of edge-current receives.

TaskStatus Particles::ClearRecvEdgeCurrents(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }
  return pbval_jedge->ClearFluxRecv();
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::ClearSendEdgeCurrents()
//! \brief Verify completion of edge-current sends.

TaskStatus Particles::ClearSendEdgeCurrents(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }
  return pbval_jedge->ClearFluxSend();
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
  if (couple_j_deposition_mode != CoupledCurrentDepositionMode::cc_convert) {
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
