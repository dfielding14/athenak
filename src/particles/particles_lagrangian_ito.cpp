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
#include <string>

#include "athena.hpp"
#include "bvals/bvals.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mesh/mesh_refinement.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"

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

void FatalIto(const std::string &msg) {
  std::cout << "### FATAL ERROR in particles_lagrangian_ito.cpp" << std::endl
            << msg << std::endl;
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

    Real p1l = fmax(-flx1(m,k,j,i)/mass, 0.0);
    Real p1r = fmax( flx1(m,k,j,i+1)/mass, 0.0);
    Real p2l = multi_d ? fmax(-flx2(m,k,j,i)/mass, 0.0) : 0.0;
    Real p2r = multi_d ? fmax( flx2(m,k,j+1,i)/mass, 0.0) : 0.0;
    Real p3l = three_d ? fmax(-flx3(m,k,j,i)/mass, 0.0) : 0.0;
    Real p3r = three_d ? fmax( flx3(m,k+1,j,i)/mass, 0.0) : 0.0;
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
    coeff(m,ITO_U1,k,j,i) = dx1*cm1/dt;
    coeff(m,ITO_U2,k,j,i) = multi_d ? dx2*cm2/dt : 0.0;
    coeff(m,ITO_U3,k,j,i) = three_d ? dx3*cm3/dt : 0.0;
    coeff(m,ITO_KAPPA1,k,j,i) = 0.5*dx1*dx1*fmax(var1, 0.0)/dt;
    coeff(m,ITO_KAPPA2,k,j,i) = multi_d ?
                                0.5*dx2*dx2*fmax(var2, 0.0)/dt : 0.0;
    coeff(m,ITO_KAPPA3,k,j,i) = three_d ?
                                0.5*dx3*dx3*fmax(var3, 0.0)/dt : 0.0;
  });

  HostArray1D<int> h_invalid("ito_invalid_host", 1);
  Kokkos::deep_copy(h_invalid, invalid);
  if (h_invalid(0) != 0) {
    FatalIto("invalid Monte-Carlo probabilities while constructing Ito-2 coefficients");
  }
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
  return pbval_ito->InitRecv(ITO_NCOEFF);
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
    pbval_ito->ProlongateCC(ito_coeff, coarse_ito_coeff);
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
  int nmb = pmy_pack->nmb_thispack;
  int gids = pmy_pack->gids;
  Real dt = pmy_pack->pmesh->dt;
  int ncycle = pmy_pack->pmesh->ncycle;
  std::int64_t rseed = random_seed;

  auto &coeff = ito_coeff;
  auto &invalid = ito_invalid;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &pi = prtcl_idata;
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

    Real u[3] = {0.0, 0.0, 0.0};
    Real kappa[3] = {0.0, 0.0, 0.0};
    int nk = three_d ? 2 : 1;
    int nj = multi_d ? 2 : 1;
    for (int dk=0; dk<nk; ++dk) {
      Real wk = three_d ? (dk == 0 ? 1.0 - wz : wz) : 1.0;
      for (int dj=0; dj<nj; ++dj) {
        Real wj = multi_d ? (dj == 0 ? 1.0 - wy : wy) : 1.0;
        for (int di=0; di<2; ++di) {
          Real wi = (di == 0 ? 1.0 - wx : wx);
          Real weight = wi*wj*wk;
          for (int n=0; n<3; ++n) {
            u[n] += weight*coeff(m,ITO_U1+n,k0+dk,j0+dj,i0+di);
            kappa[n] += weight*coeff(m,ITO_KAPPA1+n,k0+dk,j0+dj,i0+di);
          }
        }
      }
    }

    std::uint64_t base = ItoSplitMix64(static_cast<std::uint64_t>(pi(PTAG,p)));
    base ^= ItoSplitMix64(static_cast<std::uint64_t>(ncycle) +
                          0x9e3779b97f4a7c15ULL);
    base ^= ItoSplitMix64(static_cast<std::uint64_t>(rseed) +
                          0xbf58476d1ce4e5b9ULL);
    Real xi1 = sqrt_three*(2.0*ItoHashReal(base ^ 0x243f6a8885a308d3ULL) - 1.0);
    Real xi2 = sqrt_three*(2.0*ItoHashReal(base ^ 0x13198a2e03707344ULL) - 1.0);
    Real xi3 = sqrt_three*(2.0*ItoHashReal(base ^ 0xa4093822299f31d0ULL) - 1.0);
    if (!Kokkos::isfinite(u[0]) || !Kokkos::isfinite(u[1]) ||
        !Kokkos::isfinite(u[2]) || !Kokkos::isfinite(kappa[0]) ||
        !Kokkos::isfinite(kappa[1]) || !Kokkos::isfinite(kappa[2]) ||
        kappa[0] < -coeff_tol || kappa[1] < -coeff_tol || kappa[2] < -coeff_tol) {
      Kokkos::atomic_fetch_add(&invalid(0), 1);
      return;
    }
    Real dx1 = u[0]*dt + sqrt(2.0*fmax(kappa[0], 0.0)*dt)*xi1;
    Real dx2 = multi_d ? u[1]*dt + sqrt(2.0*fmax(kappa[1], 0.0)*dt)*xi2 : 0.0;
    Real dx3 = three_d ? u[2]*dt + sqrt(2.0*fmax(kappa[2], 0.0)*dt)*xi3 : 0.0;
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
