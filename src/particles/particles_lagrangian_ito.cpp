//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file particles_lagrangian_ito.cpp
//! \brief second-moment Ito tracers derived from Monte-Carlo mass-flux tracers

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "bvals/bvals.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mesh/mesh_refinement.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

KOKKOS_INLINE_FUNCTION
std::uint64_t ItoSplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

KOKKOS_INLINE_FUNCTION
Real ItoHashReal(std::uint64_t x) {
  return static_cast<Real>((ItoSplitMix64(x) >> 11) *
                           (1.0/9007199254740992.0));
}

KOKKOS_INLINE_FUNCTION
bool ItoResidualWithinTolerance(const Real residual, const Real diagonal_product,
                                const Real residual_correction,
                                const Real rel_tol) {
  Real diagonal_scale = sqrt(fabs(diagonal_product));
  Real residual_scale = fmax(diagonal_scale, residual_correction);
  return fabs(residual) <= rel_tol*residual_scale;
}

// Factor a symmetric positive-semidefinite 3x3 matrix Q = L L^T. This is the
// fixed-size form of pivoted Cholesky; keeping it unrolled matters in the
// per-particle hot path.
KOKKOS_INLINE_FUNCTION
bool ItoFactorCovariance(const Real q[6], Real l[3][3], const Real rel_tol) {
  Real a[3][3] = {
    {q[0], q[3], q[4]},
    {q[3], q[1], q[5]},
    {q[4], q[5], q[2]}
  };
  if (!Kokkos::isfinite(q[0]) || !Kokkos::isfinite(q[1]) ||
      !Kokkos::isfinite(q[2]) || !Kokkos::isfinite(q[3]) ||
      !Kokkos::isfinite(q[4]) || !Kokkos::isfinite(q[5])) {
    return false;
  }
  l[0][0] = 0.0;
  l[0][1] = 0.0;
  l[0][2] = 0.0;
  l[1][0] = 0.0;
  l[1][1] = 0.0;
  l[1][2] = 0.0;
  l[2][0] = 0.0;
  l[2][1] = 0.0;
  l[2][2] = 0.0;

  int pivot = 0;
  Real best_diagonal = a[0][0];
  if (a[1][1] > best_diagonal) {
    pivot = 1;
    best_diagonal = a[1][1];
  }
  if (a[2][2] > best_diagonal) {
    pivot = 2;
    best_diagonal = a[2][2];
  }
  int remaining1 = (pivot == 0) ? 1 : ((pivot == 1) ? 0 : 1);
  int remaining2 = (pivot == 2) ? 0 : 2;

  Real pivot_correction = 0.0;
  Real diagonal = a[pivot][pivot];
  Real pivot_scale = fmax(fabs(a[pivot][pivot]), pivot_correction);
  Real pivot_tol = rel_tol*pivot_scale;
  if (diagonal < -pivot_tol) return false;
  if (diagonal <= pivot_tol) {
    if (!ItoResidualWithinTolerance(a[0][0], a[0][0]*a[0][0], 0.0, rel_tol) ||
        !ItoResidualWithinTolerance(a[0][1], a[0][0]*a[1][1], 0.0, rel_tol) ||
        !ItoResidualWithinTolerance(a[0][2], a[0][0]*a[2][2], 0.0, rel_tol) ||
        !ItoResidualWithinTolerance(a[1][1], a[1][1]*a[1][1], 0.0, rel_tol) ||
        !ItoResidualWithinTolerance(a[1][2], a[1][1]*a[2][2], 0.0, rel_tol) ||
        !ItoResidualWithinTolerance(a[2][2], a[2][2]*a[2][2], 0.0, rel_tol)) {
      return false;
    }
    return true;
  }

  l[pivot][0] = sqrt(diagonal);
  l[remaining1][0] = a[remaining1][pivot]/l[pivot][0];
  l[remaining2][0] = a[remaining2][pivot]/l[pivot][0];
  if (!Kokkos::isfinite(l[remaining1][0]) ||
      !Kokkos::isfinite(l[remaining2][0])) {
    return false;
  }

  Real diagonal1 =
      a[remaining1][remaining1] - l[remaining1][0]*l[remaining1][0];
  Real diagonal2 =
      a[remaining2][remaining2] - l[remaining2][0]*l[remaining2][0];
  if (diagonal2 > diagonal1) {
    int swap = remaining1;
    remaining1 = remaining2;
    remaining2 = swap;
    Real swap_diagonal = diagonal1;
    diagonal1 = diagonal2;
    diagonal2 = swap_diagonal;
  }

  pivot = remaining1;
  pivot_correction = l[pivot][0]*l[pivot][0];
  diagonal = diagonal1;
  pivot_scale = fmax(fabs(a[pivot][pivot]), pivot_correction);
  pivot_tol = rel_tol*pivot_scale;
  if (diagonal < -pivot_tol) return false;
  if (diagonal <= pivot_tol) {
    Real cross_correction = l[remaining1][0]*l[remaining2][0];
    Real residual1 = diagonal1;
    Real residual12 = a[remaining1][remaining2] - cross_correction;
    Real residual2 = diagonal2;
    if (!ItoResidualWithinTolerance(
            residual1, a[remaining1][remaining1]*a[remaining1][remaining1],
            fabs(pivot_correction), rel_tol) ||
        !ItoResidualWithinTolerance(
            residual12, a[remaining1][remaining1]*a[remaining2][remaining2],
            fabs(cross_correction), rel_tol) ||
        !ItoResidualWithinTolerance(
            residual2, a[remaining2][remaining2]*a[remaining2][remaining2],
            fabs(l[remaining2][0]*l[remaining2][0]), rel_tol)) {
      return false;
    }
    return true;
  }

  l[pivot][1] = sqrt(diagonal);
  Real residual =
      a[remaining2][pivot] - l[remaining2][0]*l[pivot][0];
  l[remaining2][1] = residual/l[pivot][1];
  if (!Kokkos::isfinite(l[remaining2][1])) return false;

  pivot = remaining2;
  pivot_correction =
      l[pivot][0]*l[pivot][0] + l[pivot][1]*l[pivot][1];
  diagonal = a[pivot][pivot] - pivot_correction;
  pivot_scale = fmax(fabs(a[pivot][pivot]), pivot_correction);
  pivot_tol = rel_tol*pivot_scale;
  if (diagonal < -pivot_tol) return false;
  if (diagonal <= pivot_tol) {
    if (!ItoResidualWithinTolerance(
            diagonal, a[pivot][pivot]*a[pivot][pivot],
            fabs(l[pivot][0]*l[pivot][0]) + fabs(l[pivot][1]*l[pivot][1]),
            rel_tol)) {
      return false;
    }
    return true;
  }
  l[pivot][2] = sqrt(diagonal);
  return true;
}

void FatalIto(const std::string &msg) {
  std::cout << "### FATAL ERROR in particles_lagrangian_ito.cpp" << std::endl
            << msg << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

} // namespace

namespace particles {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::BuildItoCoefficients

TaskStatus Particles::BuildItoCoefficients(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  bool full_finite_step =
      (ito_covariance_model == ItoCovarianceModel::full_finite_step);
  Real dt = pmy_pack->pmesh->dt;
  if (dt <= 0.0) FatalIto("Ito-2 requires a positive fluid timestep");

  auto &coeff = ito_coeff;
  auto &invalid = ito_invalid;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &uold = (pmy_pack->phydro != nullptr) ? pmy_pack->phydro->u1 :
                                               pmy_pack->pmhd->u1;
  auto &saved = (pmy_pack->phydro != nullptr) ? pmy_pack->phydro->uflxidnsaved :
                                                pmy_pack->pmhd->uflxidnsaved;
  auto &flx1 = saved.x1f;
  auto &flx2 = saved.x2f;
  auto &flx3 = saved.x3f;
  constexpr Real prob_tol = (sizeof(Real) == sizeof(float)) ? 2.0e-5 : 2.0e-12;

  Kokkos::deep_copy(coeff, 0.0);
  Kokkos::deep_copy(invalid, 0);
  par_for("ito2_build_coeff", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real mass = uold(m,IDN,k,j,i);
    if (!(mass > 0.0) || !Kokkos::isfinite(mass)) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }

    Real raw1l = -flx1(m,k,j,i)/mass;
    Real raw1r =  flx1(m,k,j,i+1)/mass;
    Real raw2l = multi_d ? -flx2(m,k,j,i)/mass : 0.0;
    Real raw2r = multi_d ?  flx2(m,k,j+1,i)/mass : 0.0;
    Real raw3l = three_d ? -flx3(m,k,j,i)/mass : 0.0;
    Real raw3r = three_d ?  flx3(m,k+1,j,i)/mass : 0.0;
    if (!Kokkos::isfinite(raw1l) || !Kokkos::isfinite(raw1r) ||
        !Kokkos::isfinite(raw2l) || !Kokkos::isfinite(raw2r) ||
        !Kokkos::isfinite(raw3l) || !Kokkos::isfinite(raw3r)) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }

    Real p1l = fmax(raw1l, 0.0);
    Real p1r = fmax(raw1r, 0.0);
    Real p2l = fmax(raw2l, 0.0);
    Real p2r = fmax(raw2r, 0.0);
    Real p3l = fmax(raw3l, 0.0);
    Real p3r = fmax(raw3r, 0.0);
    Real outgoing = p1l + p1r + p2l + p2r + p3l + p3r;

    Real cp1 = p1r + p1l, cm1 = p1r - p1l;
    Real cp2 = p2r + p2l, cm2 = p2r - p2l;
    Real cp3 = p3r + p3l, cm3 = p3r - p3l;
    Real var1 = cp1 - cm1*cm1;
    Real var2 = cp2 - cm2*cm2;
    Real var3 = cp3 - cm3*cm3;
    if (!Kokkos::isfinite(outgoing) || outgoing > 1.0 + prob_tol ||
        var1 < -prob_tol || var2 < -prob_tol || var3 < -prob_tol) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
    }

    Real dx1 = mbsize.d_view(m).dx1;
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    Real mean1 = dx1*cm1;
    Real mean2 = multi_d ? dx2*cm2 : 0.0;
    Real mean3 = three_d ? dx3*cm3 : 0.0;
    coeff(m,ITO_M1,k,j,i) = mean1;
    coeff(m,ITO_M2,k,j,i) = mean2;
    coeff(m,ITO_M3,k,j,i) = mean3;
    coeff(m,ITO_Q11,k,j,i) = dx1*dx1*fmax(var1, 0.0);
    coeff(m,ITO_Q22,k,j,i) = multi_d ? dx2*dx2*fmax(var2, 0.0) : 0.0;
    coeff(m,ITO_Q33,k,j,i) = three_d ? dx3*dx3*fmax(var3, 0.0) : 0.0;
    if (full_finite_step) {
      coeff(m,ITO_Q12,k,j,i) = multi_d ? -mean1*mean2 : 0.0;
      coeff(m,ITO_Q13,k,j,i) = three_d ? -mean1*mean3 : 0.0;
      coeff(m,ITO_Q23,k,j,i) = three_d ? -mean2*mean3 : 0.0;
    }
  });

  HostArray1D<int> h_invalid("ito_invalid_host", 1);
  Kokkos::deep_copy(h_invalid, invalid);
  if (h_invalid(0) != 0) {
    FatalIto("invalid Monte-Carlo probabilities while constructing Ito-2 coefficients");
  }

  int nx1 = ie - is + 1;
  int nx2 = je - js + 1;
  int nx3 = ke - ks + 1;
  std::int64_t nmkji = static_cast<std::int64_t>(nmb)*nx3*nx2*nx1;
  std::int64_t nkji = static_cast<std::int64_t>(nx3)*nx2*nx1;
  std::int64_t nji = static_cast<std::int64_t>(nx2)*nx1;
  Real max_outgoing = 0.0;
  Kokkos::parallel_reduce("ito2_max_outgoing",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<std::int64_t>>(0, nmkji),
  KOKKOS_LAMBDA(const std::int64_t idx, Real &max_value) {
    int m = static_cast<int>(idx/nkji);
    std::int64_t remainder = idx - static_cast<std::int64_t>(m)*nkji;
    int kk = static_cast<int>(remainder/nji);
    int jj = static_cast<int>((remainder - static_cast<std::int64_t>(kk)*nji)/nx1);
    int ii = static_cast<int>(
        remainder - static_cast<std::int64_t>(kk)*nji -
        static_cast<std::int64_t>(jj)*nx1);
    int k = kk + ks;
    int j = jj + js;
    int i = ii + is;
    Real dx1 = mbsize.d_view(m).dx1;
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    Real mean1 = coeff(m,ITO_M1,k,j,i);
    Real mean2 = coeff(m,ITO_M2,k,j,i);
    Real mean3 = coeff(m,ITO_M3,k,j,i);
    Real outgoing = (coeff(m,ITO_Q11,k,j,i) + mean1*mean1)/(dx1*dx1);
    if (multi_d) {
      outgoing += (coeff(m,ITO_Q22,k,j,i) + mean2*mean2)/(dx2*dx2);
    }
    if (three_d) {
      outgoing += (coeff(m,ITO_Q33,k,j,i) + mean3*mean3)/(dx3*dx3);
    }
    max_value = fmax(max_value, outgoing);
  }, Kokkos::Max<Real>(max_outgoing));

  dtnew = (max_outgoing > 0.0) ?
          dt*ito_probability_target/max_outgoing :
          std::numeric_limits<Real>::max();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------

TaskStatus Particles::RestrictItoCoefficients(Driver *pdriver, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pmy_pack->pmesh->pmr->RestrictCC(ito_coeff, coarse_ito_coeff);
  }
  return TaskStatus::complete;
}

TaskStatus Particles::InitRecvItoCoefficients(Driver *pdriver, int stage) {
  return pbval_ito->InitRecv(ito_ncoeff);
}

TaskStatus Particles::SendItoCoefficients(Driver *pdriver, int stage) {
  return pbval_ito->PackAndSendCC(ito_coeff, coarse_ito_coeff);
}

TaskStatus Particles::RecvItoCoefficients(Driver *pdriver, int stage) {
  return pbval_ito->RecvAndUnpackCC(ito_coeff, coarse_ito_coeff);
}

TaskStatus Particles::ClearRecvItoCoefficients(Driver *pdriver, int stage) {
  return pbval_ito->ClearRecv();
}

TaskStatus Particles::ClearSendItoCoefficients(Driver *pdriver, int stage) {
  return pbval_ito->ClearSend();
}

TaskStatus Particles::ProlongateItoCoefficients(Driver *pdriver, int stage) {
  if (pmy_pack->pmesh->multilevel) {
    pbval_ito->FillCoarseInBndryCC(ito_coeff, coarse_ito_coeff);

    int nmb = pmy_pack->nmb_thispack;
    int nnghbr = pmy_pack->pmb->nnghbr;
    auto &nghbr = pmy_pack->pmb->nghbr;
    auto &mblev = pmy_pack->pmb->mb_lev;
    auto &rbuf = pbval_ito->recvbuf;
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    bool multi_d = pmy_pack->pmesh->multi_d;
    bool three_d = pmy_pack->pmesh->three_d;
    int ncoeff = ito_ncoeff;
    auto &coeff = ito_coeff;
    auto &coarse_coeff = coarse_ito_coeff;

    Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb*nnghbr, Kokkos::AUTO);
    Kokkos::parallel_for("ito2_prolong_central_coefficients", policy,
    KOKKOS_LAMBDA(TeamMember_t tmember) {
      int m = tmember.league_rank()/nnghbr;
      int n = tmember.league_rank() - m*nnghbr;
      if (nghbr.d_view(m,n).gid >= 0 &&
          nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        int il = rbuf[n].iprol[0].bis;
        int iu = rbuf[n].iprol[0].bie;
        int jl = rbuf[n].iprol[0].bjs;
        int ju = rbuf[n].iprol[0].bje;
        int kl = rbuf[n].iprol[0].bks;
        int ku = rbuf[n].iprol[0].bke;
        int ni = iu - il + 1;
        int nj = ju - jl + 1;
        int nkji = (ku - kl + 1)*nj*ni;
        int nji = nj*ni;

        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
        [&](const int idx) {
          int k = idx/nji + kl;
          int j = (idx % nji)/ni + jl;
          int i = idx % ni + il;
          int fi = (i - indcs.cis)*2 + indcs.is;
          int fj = (j - indcs.cjs)*2 + indcs.js;
          int fk = (k - indcs.cks)*2 + indcs.ks;
          int nchild2 = multi_d ? 2 : 1;
          int nchild3 = three_d ? 2 : 1;

          for (int dk=0; dk<nchild3; ++dk) {
            for (int dj=0; dj<nchild2; ++dj) {
              for (int di=0; di<2; ++di) {
                for (int v=0; v<ncoeff; ++v) {
                  coeff(m,v,fk+dk,fj+dj,fi+di) = coarse_coeff(m,v,k,j,i);
                }
              }
            }
          }
        });
      }
      tmember.team_barrier();
    });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::PushIto2

TaskStatus Particles::PushIto2(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  bool full_finite_step =
      (ito_covariance_model == ItoCovarianceModel::full_finite_step);
  int nmb = pmy_pack->nmb_thispack;
  int gids = pmy_pack->gids;
  int ncycle = pmy_pack->pmesh->ncycle;
  std::int64_t rseed = random_seed;

  auto &coeff = ito_coeff;
  auto &invalid = ito_invalid;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &pi = prtcl_idata;
  auto &ptag = prtcl_tag;
  auto &pr = prtcl_rdata;
  constexpr Real sqrt_three = 1.7320508075688772935;
  constexpr Real coeff_tol = (sizeof(Real) == sizeof(float)) ? 2.0e-5 : 2.0e-12;

  Kokkos::deep_copy(invalid, 0);
  par_for("ito2_push", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }

    Real gx = (pr(LMCX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1
              + is - 0.5;
    Real gy = multi_d ?
              (pr(LMCY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2
              + js - 0.5 : static_cast<Real>(js);
    Real gz = three_d ?
              (pr(LMCZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3
              + ks - 0.5 : static_cast<Real>(ks);
    int i0 = static_cast<int>(floor(gx));
    int j0 = multi_d ? static_cast<int>(floor(gy)) : js;
    int k0 = three_d ? static_cast<int>(floor(gz)) : ks;
    Real wx = gx - i0;
    Real wy = multi_d ? gy - j0 : 0.0;
    Real wz = three_d ? gz - k0 : 0.0;
    if (i0 < is-1 || i0+1 > ie+1 ||
        (multi_d && (j0 < js-1 || j0+1 > je+1)) ||
        (three_d && (k0 < ks-1 || k0+1 > ke+1))) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }

    Real mean[3] = {0.0, 0.0, 0.0};
    Real q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int nk = three_d ? 2 : 1;
    int nj = multi_d ? 2 : 1;
    int nq = full_finite_step ? 6 : 3;
    for (int dk=0; dk<nk; ++dk) {
      Real wk = three_d ? (dk == 0 ? 1.0 - wz : wz) : 1.0;
      for (int dj=0; dj<nj; ++dj) {
        Real wj = multi_d ? (dj == 0 ? 1.0 - wy : wy) : 1.0;
        for (int di=0; di<2; ++di) {
          Real wi = (di == 0 ? 1.0 - wx : wx);
          Real weight = wi*wj*wk;
          for (int n=0; n<3; ++n) {
            mean[n] += weight*coeff(m,ITO_M1+n,k0+dk,j0+dj,i0+di);
          }
          for (int n=0; n<nq; ++n) {
            q[n] += weight*coeff(m,ITO_Q11+n,k0+dk,j0+dj,i0+di);
          }
        }
      }
    }

    std::uint64_t base = ItoSplitMix64(ptag(p));
    base ^= ItoSplitMix64(static_cast<std::uint64_t>(ncycle) +
                          0x9e3779b97f4a7c15ULL);
    base ^= ItoSplitMix64(static_cast<std::uint64_t>(rseed) +
                          0xbf58476d1ce4e5b9ULL);
    Real xi1 = sqrt_three*(2.0*ItoHashReal(base ^ 0x243f6a8885a308d3ULL) - 1.0);
    Real xi2 = sqrt_three*(2.0*ItoHashReal(base ^ 0x13198a2e03707344ULL) - 1.0);
    Real xi3 = sqrt_three*(2.0*ItoHashReal(base ^ 0xa4093822299f31d0ULL) - 1.0);
    if (!Kokkos::isfinite(mean[0]) || !Kokkos::isfinite(mean[1]) ||
        !Kokkos::isfinite(mean[2]) || !Kokkos::isfinite(q[0]) ||
        !Kokkos::isfinite(q[1]) || !Kokkos::isfinite(q[2]) ||
        !Kokkos::isfinite(q[3]) || !Kokkos::isfinite(q[4]) ||
        !Kokkos::isfinite(q[5])) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }
    Real stochastic[3] = {0.0, 0.0, 0.0};
    if (full_finite_step) {
      Real l[3][3];
      if (!ItoFactorCovariance(q, l, coeff_tol)) {
        Kokkos::atomic_fetch_add(&invalid(0), 1);
        return;
      }
      Real xi[3] = {xi1, xi2, xi3};
      for (int row=0; row<3; ++row) {
        for (int column=0; column<3; ++column) {
          stochastic[row] += l[row][column]*xi[column];
        }
      }
    } else {
      if (q[0] < -coeff_tol || q[1] < -coeff_tol || q[2] < -coeff_tol) {
        Kokkos::atomic_fetch_add(&invalid(0), 1);
        return;
      }
      stochastic[0] = sqrt(fmax(q[0], 0.0))*xi1;
      stochastic[1] = sqrt(fmax(q[1], 0.0))*xi2;
      stochastic[2] = sqrt(fmax(q[2], 0.0))*xi3;
    }
    Real dx1 = mean[0] + stochastic[0];
    Real dx2 = multi_d ? mean[1] + stochastic[1] : 0.0;
    Real dx3 = three_d ? mean[2] + stochastic[2] : 0.0;
    Real lx = mbsize.d_view(m).x1max - mbsize.d_view(m).x1min;
    Real ly = mbsize.d_view(m).x2max - mbsize.d_view(m).x2min;
    Real lz = mbsize.d_view(m).x3max - mbsize.d_view(m).x3min;
    if (!Kokkos::isfinite(dx1) || !Kokkos::isfinite(dx2) || !Kokkos::isfinite(dx3) ||
        fabs(dx1) >= lx || (multi_d && fabs(dx2) >= ly) ||
        (three_d && fabs(dx3) >= lz)) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }

    pr(LMCX,p) += dx1;
    if (multi_d) pr(LMCY,p) += dx2;
    if (three_d) pr(LMCZ,p) += dx3;
    pi(PLASTMOVE,p) = 0;
    pi(PLASTLEVEL,p) = mblev.d_view(m);
  });

  HostArray1D<int> h_invalid("ito_push_invalid_host", 1);
  Kokkos::deep_copy(h_invalid, invalid);
  if (h_invalid(0) != 0) {
    FatalIto("Ito-2 particle push was invalid or crossed multiple MeshBlocks per axis");
  }
  return TaskStatus::complete;
}

} // namespace particles
