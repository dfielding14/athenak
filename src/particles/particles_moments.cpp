//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_moments.cpp
//! \brief task wrappers and kernels for particle moment deposition/communication

#include <cmath>
#include <iostream>

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

KOKKOS_INLINE_FUNCTION
bool InBoundsInclusive(const int i, const int imin, const int imax) {
  return (i >= imin && i <= imax);
}

template <bool STAGGERED, int O>
KOKKOS_INLINE_FUNCTION
void ShapeOrder(const int i, const Real di, int &i_min, Real S[O + 1]) {
  static_assert((O == 1 || O == 2),
                "ShapeOrder currently supports only O=1 and O=2");
  if constexpr (O == 1) {
    if constexpr (!STAGGERED) {
      i_min = i;
      S[0] = static_cast<Real>(1.0) - di;
      S[1] = di;
    } else {
      if (di < static_cast<Real>(0.5)) {
        i_min = i - 1;
        S[0] = static_cast<Real>(0.5) - di;
        S[1] = static_cast<Real>(1.0) - S[0];
      } else {
        i_min = i;
        S[0] = static_cast<Real>(1.5) - di;
        S[1] = static_cast<Real>(1.0) - S[0];
      }
    }
  } else {
    if constexpr (!STAGGERED) {
      if (di < static_cast<Real>(0.5)) {
        i_min = i - 1;
        S[0] = static_cast<Real>(0.5) *
               (static_cast<Real>(0.5) - di) *
               (static_cast<Real>(0.5) - di);
        S[1] = static_cast<Real>(0.75) - di*di;
        S[2] = static_cast<Real>(1.0) - S[0] - S[1];
      } else {
        i_min = i;
        S[0] = static_cast<Real>(0.5) *
               (static_cast<Real>(1.5) - di) *
               (static_cast<Real>(1.5) - di);
        Real d1 = static_cast<Real>(1.0) - di;
        S[1] = static_cast<Real>(0.75) - d1*d1;
        S[2] = static_cast<Real>(1.0) - S[0] - S[1];
      }
    } else {
      i_min = i - 1;
      S[0] = static_cast<Real>(0.5) *
             (static_cast<Real>(1.0) - di) *
             (static_cast<Real>(1.0) - di);
      S[2] = static_cast<Real>(0.5)*di*di;
      S[1] = static_cast<Real>(1.0) - S[0] - S[2];
    }
  }
}

template <int O>
KOKKOS_INLINE_FUNCTION
void ShapeForDeposit(const int i_init, const Real di_init,
                     const int i_fin, const Real di_fin,
                     int &i_min, int &i_max,
                     Real iS[O + 2], Real fS[O + 2]) {
  static_assert((O == 1 || O == 2),
                "ShapeForDeposit currently supports only O=1 and O=2");
  int i_init_min = 0;
  int i_fin_min = 0;
  Real iS_local[O + 1];
  Real fS_local[O + 1];

  ShapeOrder<false, O>(i_init, di_init, i_init_min, iS_local);
  ShapeOrder<false, O>(i_fin, di_fin, i_fin_min, fS_local);

  if (i_init_min < i_fin_min) {
    i_min = i_init_min;
    i_max = i_min + O + 1;
    for (int j = 0; j < O + 1; ++j) iS[j] = iS_local[j];
    iS[O + 1] = static_cast<Real>(0.0);
    fS[0] = static_cast<Real>(0.0);
    for (int j = 0; j < O + 1; ++j) fS[j + 1] = fS_local[j];
  } else if (i_init_min > i_fin_min) {
    i_min = i_fin_min;
    i_max = i_min + O + 1;
    iS[0] = static_cast<Real>(0.0);
    for (int j = 0; j < O + 1; ++j) iS[j + 1] = iS_local[j];
    for (int j = 0; j < O + 1; ++j) fS[j] = fS_local[j];
    fS[O + 1] = static_cast<Real>(0.0);
  } else {
    i_min = i_init_min;
    i_max = i_min + O;
    for (int j = 0; j < O + 1; ++j) iS[j] = iS_local[j];
    iS[O + 1] = static_cast<Real>(0.0);
    for (int j = 0; j < O + 1; ++j) fS[j] = fS_local[j];
    fS[O + 1] = static_cast<Real>(0.0);
  }
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

  if (deposit_order == 1) {
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

  if (deposit_order == 2) {
    par_for("deposit_direct_edge_currents_o2", DevExeSpace(), 0, npart-1,
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

      Real iS_x1[4];
      Real fS_x1[4];
      int i1_min = 0;
      int i1_max = 0;
      ShapeForDeposit<2>(i0, dxo, i1, dxn, i1_min, i1_max, iS_x1, fS_x1);

      if (pmy_pack->pmesh->one_d) {
        Real wx1[4];
        Real wx23[4];
        for (int i = 0; i < 4; ++i) {
          wx1[i] = fS_x1[i] - iS_x1[i];
          wx23[i] = static_cast<Real>(0.5)*(fS_x1[i] + iS_x1[i]);
        }

        Real jx1_contrib[4];
        Real qdx1dt = coeff;
        Real qvy = q_macro*vy;
        Real qvz = q_macro*vz;
        jx1_contrib[0] = -qdx1dt*wx1[0];
        for (int i = 1; i < 4; ++i) {
          jx1_contrib[i] = jx1_contrib[i-1] - qdx1dt*wx1[i];
        }

        const int di_x1 = i1_max - i1_min;
        for (int i = 0; i <= di_x1; ++i) {
          int ii = i1_min + i;
          if (InBoundsInclusive(ii, i_min, i_max + 1)) {
            Real jy_val = qvy*wx23[i];
            Real jz_val = qvz*wx23[i];
            Kokkos::atomic_add(&jy_e(m, ks, js, ii), jy_val);
            Kokkos::atomic_add(&jy_e(m, ke+1, js, ii), jy_val);
            Kokkos::atomic_add(&jz_e(m, ks, js, ii), jz_val);
            Kokkos::atomic_add(&jz_e(m, ks, je+1, ii), jz_val);
          }
        }
        return;
      }

      Real iS_x2[4];
      Real fS_x2[4];
      int i2_min = 0;
      int i2_max = 0;
      ShapeForDeposit<2>(j0, dyo, j1, dyn, i2_min, i2_max, iS_x2, fS_x2);

      if (pmy_pack->pmesh->two_d) {
        Real wx1[4][4];
        Real wx2[4][4];
        Real wx3[4][4];
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            wx1[i][j] = static_cast<Real>(0.5)*(fS_x1[i] - iS_x1[i])*
                        (fS_x2[j] + iS_x2[j]);
            wx2[i][j] = static_cast<Real>(0.5)*(fS_x1[i] + iS_x1[i])*
                        (fS_x2[j] - iS_x2[j]);
            wx3[i][j] = static_cast<Real>(1.0/3.0)*
                        (fS_x2[j]*(static_cast<Real>(0.5)*iS_x1[i] + fS_x1[i]) +
                         iS_x2[j]*(static_cast<Real>(0.5)*fS_x1[i] + iS_x1[i]));
          }
        }

        Real jx1_contrib[4][4];
        Real jy1_contrib[4][4];
        Real qdx1dt = coeff;
        Real qdx2dt = coeff;
        Real qvz = q_macro*vz;

        for (int j = 0; j < 4; ++j) {
          jx1_contrib[0][j] = -qdx1dt*wx1[0][j];
        }
        for (int i = 1; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            jx1_contrib[i][j] = jx1_contrib[i-1][j] - qdx1dt*wx1[i][j];
          }
        }

        for (int i = 0; i < 4; ++i) {
          jy1_contrib[i][0] = -qdx2dt*wx2[i][0];
        }
        for (int j = 1; j < 4; ++j) {
          for (int i = 0; i < 4; ++i) {
            jy1_contrib[i][j] = jy1_contrib[i][j-1] - qdx2dt*wx2[i][j];
          }
        }

        int di_x1 = i1_max - i1_min;
        int di_x2 = i2_max - i2_min;
        int kf = ks;
        int kf2 = ke + 1;

        for (int i = 0; i < di_x1; ++i) {
          int ii = i1_min + i;
          if (!InBoundsInclusive(ii, i_min, i_max)) continue;
          for (int j = 0; j <= di_x2; ++j) {
            int jj = i2_min + j;
            if (!InBoundsInclusive(jj, j_min, j_max + 1)) continue;
            Real val = jx1_contrib[i][j];
            Kokkos::atomic_add(&jx_e(m, kf, jj, ii), val);
            Kokkos::atomic_add(&jx_e(m, kf2, jj, ii), val);
          }
        }

        for (int i = 0; i <= di_x1; ++i) {
          int ii = i1_min + i;
          if (!InBoundsInclusive(ii, i_min, i_max + 1)) continue;
          for (int j = 0; j < di_x2; ++j) {
            int jj = i2_min + j;
            if (!InBoundsInclusive(jj, j_min, j_max)) continue;
            Real val = jy1_contrib[i][j];
            Kokkos::atomic_add(&jy_e(m, kf, jj, ii), val);
            Kokkos::atomic_add(&jy_e(m, kf2, jj, ii), val);
          }
        }

        for (int i = 0; i <= di_x1; ++i) {
          int ii = i1_min + i;
          if (!InBoundsInclusive(ii, i_min, i_max + 1)) continue;
          for (int j = 0; j <= di_x2; ++j) {
            int jj = i2_min + j;
            if (!InBoundsInclusive(jj, j_min, j_max + 1)) continue;
            Kokkos::atomic_add(&jz_e(m, ks, jj, ii), qvz*wx3[i][j]);
          }
        }
        return;
      }

      Real iS_x3[4];
      Real fS_x3[4];
      int i3_min = 0;
      int i3_max = 0;
      ShapeForDeposit<2>(k0, dzo, k1, dzn, i3_min, i3_max, iS_x3, fS_x3);

      Real wx1[4][4][4];
      Real wx2[4][4][4];
      Real wx3[4][4][4];
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 4; ++k) {
            wx1[i][j][k] = static_cast<Real>(1.0/3.0)*
                           (fS_x1[i] - iS_x1[i])*
                           ((iS_x2[j]*iS_x3[k] + fS_x2[j]*fS_x3[k]) +
                            static_cast<Real>(0.5)*
                            (iS_x3[k]*fS_x2[j] + iS_x2[j]*fS_x3[k]));
            wx2[i][j][k] = static_cast<Real>(1.0/3.0)*
                           (fS_x2[j] - iS_x2[j])*
                           ((iS_x1[i]*iS_x3[k] + fS_x1[i]*fS_x3[k]) +
                            static_cast<Real>(0.5)*
                            (iS_x3[k]*fS_x1[i] + iS_x1[i]*fS_x3[k]));
            wx3[i][j][k] = static_cast<Real>(1.0/3.0)*
                           (fS_x3[k] - iS_x3[k])*
                           ((iS_x1[i]*iS_x2[j] + fS_x1[i]*fS_x2[j]) +
                            static_cast<Real>(0.5)*
                            (iS_x1[i]*fS_x2[j] + iS_x2[j]*fS_x1[i]));
          }
        }
      }

      Real jx1_contrib[4][4][4];
      Real jy1_contrib[4][4][4];
      Real jz1_contrib[4][4][4];
      Real qdx1dt = coeff;
      Real qdx2dt = coeff;
      Real qdx3dt = coeff;

      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
          jx1_contrib[0][j][k] = -qdx1dt*wx1[0][j][k];
        }
      }
      for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 4; ++k) {
            jx1_contrib[i][j][k] =
                jx1_contrib[i-1][j][k] - qdx1dt*wx1[i][j][k];
          }
        }
      }

      for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
          jy1_contrib[i][0][k] = -qdx2dt*wx2[i][0][k];
        }
      }
      for (int i = 0; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
          for (int k = 0; k < 4; ++k) {
            jy1_contrib[i][j][k] =
                jy1_contrib[i][j-1][k] - qdx2dt*wx2[i][j][k];
          }
        }
      }

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          jz1_contrib[i][j][0] = -qdx3dt*wx3[i][j][0];
        }
      }
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          for (int k = 1; k < 4; ++k) {
            jz1_contrib[i][j][k] =
                jz1_contrib[i][j][k-1] - qdx3dt*wx3[i][j][k];
          }
        }
      }

      int di_x1 = i1_max - i1_min;
      int di_x2 = i2_max - i2_min;
      int di_x3 = i3_max - i3_min;

      for (int i = 0; i < di_x1; ++i) {
        int ii = i1_min + i;
        if (!InBoundsInclusive(ii, i_min, i_max)) continue;
        for (int j = 0; j <= di_x2; ++j) {
          int jj = i2_min + j;
          if (!InBoundsInclusive(jj, j_min, j_max + 1)) continue;
          for (int k = 0; k <= di_x3; ++k) {
            int kk = i3_min + k;
            if (!InBoundsInclusive(kk, k_min, k_max + 1)) continue;
            Kokkos::atomic_add(&jx_e(m, kk, jj, ii), jx1_contrib[i][j][k]);
          }
        }
      }

      for (int i = 0; i <= di_x1; ++i) {
        int ii = i1_min + i;
        if (!InBoundsInclusive(ii, i_min, i_max + 1)) continue;
        for (int j = 0; j < di_x2; ++j) {
          int jj = i2_min + j;
          if (!InBoundsInclusive(jj, j_min, j_max)) continue;
          for (int k = 0; k <= di_x3; ++k) {
            int kk = i3_min + k;
            if (!InBoundsInclusive(kk, k_min, k_max + 1)) continue;
            Kokkos::atomic_add(&jy_e(m, kk, jj, ii), jy1_contrib[i][j][k]);
          }
        }
      }

      for (int i = 0; i <= di_x1; ++i) {
        int ii = i1_min + i;
        if (!InBoundsInclusive(ii, i_min, i_max + 1)) continue;
        for (int j = 0; j <= di_x2; ++j) {
          int jj = i2_min + j;
          if (!InBoundsInclusive(jj, j_min, j_max + 1)) continue;
          for (int k = 0; k < di_x3; ++k) {
            int kk = i3_min + k;
            if (!InBoundsInclusive(kk, k_min, k_max)) continue;
            Kokkos::atomic_add(&jz_e(m, kk, jj, ii), jz1_contrib[i][j][k]);
          }
        }
      }
    });
    return TaskStatus::complete;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Unsupported deposit_order=" << deposit_order
            << " for direct edge-current deposition (supported: 1, 2)"
            << std::endl;
  std::exit(EXIT_FAILURE);
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
//! \fn TaskStatus Particles::ApplyEdgeCurrentPhysicalBCs()
//! \brief Apply physical BCs to direct edge-current state before EFieldSrc.

TaskStatus Particles::ApplyEdgeCurrentPhysicalBCs(Driver *pdriver, int stage) {
  (void)pdriver;
  if (pbval_jedge == nullptr) return TaskStatus::complete;
  if (!RunEdgeCurrentWrappersAtStage(deposit_moments, couple_moments_to_mhd,
                                     couple_j_to_efield_representation,
                                     couple_j_deposition_mode, stage)) {
    return TaskStatus::complete;
  }

  auto &pm = pmy_pack->pmesh;
  if (pm->strictly_periodic) return TaskStatus::complete;

  auto &indcs = pm->mb_indcs;
  int ng = indcs.ng;
  if (ng <= 0) return TaskStatus::complete;

  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int nmb = pmy_pack->nmb_thispack;
  if (nmb <= 0) return TaskStatus::complete;

  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto jx_e = j_edge_x1e;
  auto jy_e = j_edge_x2e;
  auto jz_e = j_edge_x3e;

  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
    par_for("prtcl_jedge_bc_x1_jx", DevExeSpace(), 0, nmb - 1, 0, n3, 0, n2,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      BoundaryFlag bix = mb_bcs.d_view(m, BoundaryFace::inner_x1);
      switch (bix) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jx_e(m, k, j, is - i - 1) = -jx_e(m, k, j, is + i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jx_e(m, k, j, is - i - 1) = jx_e(m, k, j, is);
          }
          break;
        default:
          break;
      }

      BoundaryFlag box = mb_bcs.d_view(m, BoundaryFace::outer_x1);
      switch (box) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jx_e(m, k, j, ie + i + 1) = -jx_e(m, k, j, ie - i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jx_e(m, k, j, ie + i + 1) = jx_e(m, k, j, ie);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x1_jy", DevExeSpace(), 0, nmb - 1, 0, n3, 0, n2 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      BoundaryFlag bix = mb_bcs.d_view(m, BoundaryFace::inner_x1);
      switch (bix) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jy_e(m, k, j, is - i - 1) = jy_e(m, k, j, is + i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jy_e(m, k, j, is - i - 1) = jy_e(m, k, j, is);
          }
          break;
        default:
          break;
      }

      BoundaryFlag box = mb_bcs.d_view(m, BoundaryFace::outer_x1);
      switch (box) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jy_e(m, k, j, ie + i + 2) = jy_e(m, k, j, ie - i + 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jy_e(m, k, j, ie + i + 2) = jy_e(m, k, j, ie + 1);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x1_jz", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      BoundaryFlag bix = mb_bcs.d_view(m, BoundaryFace::inner_x1);
      switch (bix) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jz_e(m, k, j, is - i - 1) = jz_e(m, k, j, is + i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jz_e(m, k, j, is - i - 1) = jz_e(m, k, j, is);
          }
          break;
        default:
          break;
      }

      BoundaryFlag box = mb_bcs.d_view(m, BoundaryFace::outer_x1);
      switch (box) {
        case BoundaryFlag::reflect:
          for (int i = 0; i < ng; ++i) {
            jz_e(m, k, j, ie + i + 2) = jz_e(m, k, j, ie - i + 1);
          }
          break;
        case BoundaryFlag::outflow:
          for (int i = 0; i < ng; ++i) {
            jz_e(m, k, j, ie + i + 2) = jz_e(m, k, j, ie + 1);
          }
          break;
        default:
          break;
      }
    });
  }

  if (pm->one_d) return TaskStatus::complete;

  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    par_for("prtcl_jedge_bc_x2_jx", DevExeSpace(), 0, nmb - 1, 0, n3, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      BoundaryFlag biy = mb_bcs.d_view(m, BoundaryFace::inner_x2);
      switch (biy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jx_e(m, k, js - j - 1, i) = jx_e(m, k, js + j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jx_e(m, k, js - j - 1, i) = jx_e(m, k, js, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boy = mb_bcs.d_view(m, BoundaryFace::outer_x2);
      switch (boy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jx_e(m, k, je + j + 2, i) = jx_e(m, k, je - j + 1, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jx_e(m, k, je + j + 2, i) = jx_e(m, k, je + 1, i);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x2_jy", DevExeSpace(), 0, nmb - 1, 0, n3, 0, n1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      BoundaryFlag biy = mb_bcs.d_view(m, BoundaryFace::inner_x2);
      switch (biy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jy_e(m, k, js - j - 1, i) = -jy_e(m, k, js + j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jy_e(m, k, js - j - 1, i) = jy_e(m, k, js, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boy = mb_bcs.d_view(m, BoundaryFace::outer_x2);
      switch (boy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jy_e(m, k, je + j + 1, i) = -jy_e(m, k, je - j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jy_e(m, k, je + j + 1, i) = jy_e(m, k, je, i);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x2_jz", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      BoundaryFlag biy = mb_bcs.d_view(m, BoundaryFace::inner_x2);
      switch (biy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jz_e(m, k, js - j - 1, i) = jz_e(m, k, js + j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jz_e(m, k, js - j - 1, i) = jz_e(m, k, js, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boy = mb_bcs.d_view(m, BoundaryFace::outer_x2);
      switch (boy) {
        case BoundaryFlag::reflect:
          for (int j = 0; j < ng; ++j) {
            jz_e(m, k, je + j + 2, i) = jz_e(m, k, je - j + 1, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int j = 0; j < ng; ++j) {
            jz_e(m, k, je + j + 2, i) = jz_e(m, k, je + 1, i);
          }
          break;
        default:
          break;
      }
    });
  }

  if (pm->two_d) return TaskStatus::complete;

  if (pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic) {
    par_for("prtcl_jedge_bc_x3_jx", DevExeSpace(), 0, nmb - 1, 0, n2, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      BoundaryFlag biz = mb_bcs.d_view(m, BoundaryFace::inner_x3);
      switch (biz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jx_e(m, ks - k - 1, j, i) = jx_e(m, ks + k, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jx_e(m, ks - k - 1, j, i) = jx_e(m, ks, j, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boz = mb_bcs.d_view(m, BoundaryFace::outer_x3);
      switch (boz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jx_e(m, ke + k + 2, j, i) = jx_e(m, ke - k + 1, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jx_e(m, ke + k + 2, j, i) = jx_e(m, ke + 1, j, i);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x3_jy", DevExeSpace(), 0, nmb - 1, 0, n2 - 1, 0, n1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      BoundaryFlag biz = mb_bcs.d_view(m, BoundaryFace::inner_x3);
      switch (biz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jy_e(m, ks - k - 1, j, i) = jy_e(m, ks + k, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jy_e(m, ks - k - 1, j, i) = jy_e(m, ks, j, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boz = mb_bcs.d_view(m, BoundaryFace::outer_x3);
      switch (boz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jy_e(m, ke + k + 2, j, i) = jy_e(m, ke - k + 1, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jy_e(m, ke + k + 2, j, i) = jy_e(m, ke + 1, j, i);
          }
          break;
        default:
          break;
      }
    });

    par_for("prtcl_jedge_bc_x3_jz", DevExeSpace(), 0, nmb - 1, 0, n2, 0, n1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      BoundaryFlag biz = mb_bcs.d_view(m, BoundaryFace::inner_x3);
      switch (biz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jz_e(m, ks - k - 1, j, i) = -jz_e(m, ks + k, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jz_e(m, ks - k - 1, j, i) = jz_e(m, ks, j, i);
          }
          break;
        default:
          break;
      }

      BoundaryFlag boz = mb_bcs.d_view(m, BoundaryFace::outer_x3);
      switch (boz) {
        case BoundaryFlag::reflect:
          for (int k = 0; k < ng; ++k) {
            jz_e(m, ke + k + 1, j, i) = -jz_e(m, ke - k, j, i);
          }
          break;
        case BoundaryFlag::outflow:
          for (int k = 0; k < ng; ++k) {
            jz_e(m, ke + k + 1, j, i) = jz_e(m, ke, j, i);
          }
          break;
        default:
          break;
      }
    });
  }

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
