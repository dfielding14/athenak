//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"
#include "relativistic_pusher.hpp"
#include "trilinear_gather.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace particles {
namespace {

struct ParticleBField {
  Real x, y, z;
};

KOKKOS_INLINE_FUNCTION
bool IsCellCenteredTrilinearCoordinateSupported(Real x, Real xmin, Real dx,
                                                int is, int nx, int ng) {
  if (!relativistic::IsFinite(x) || !relativistic::IsFinite(xmin) ||
      !relativistic::IsFinite(dx) || dx <= 0.0) {
    return false;
  }
  const Real lower = xmin + (0.5 - static_cast<Real>(is))*dx;
  const Real upper =
      xmin + (static_cast<Real>(nx + 2*ng - is) - 0.5)*dx;
  return relativistic::IsFinite(lower) && relativistic::IsFinite(upper) &&
         x >= lower && x < upper;
}

KOKKOS_INLINE_FUNCTION
int StepsForRatio(Real ratio) {
  if (ratio <= 1.0) {return 1;}
  if (!relativistic::IsFinite(ratio) ||
      ratio >= static_cast<Real>(std::numeric_limits<int>::max())) {
    return std::numeric_limits<int>::max();
  }
  int n = static_cast<int>(ratio);
  if (static_cast<Real>(n) < ratio) {++n;}
  return n;
}

struct LinLegacyGather {
  template<typename FaceField, typename CellField, typename BlockSize>
  KOKKOS_INLINE_FUNCTION
  static ParticleBField Gather(const FaceField &b0, const CellField &bcc,
                               const BlockSize &mbsize, int m, Real x1, Real x2,
                               Real x3, int is, int js, int ks, int nx1, int nx2,
                               int nx3, bool multi_d, bool three_d) {
    (void)bcc;
    int ip = static_cast<int>((x1 - mbsize.d_view(m).x1min)/
                              mbsize.d_view(m).dx1) + is;
    int jp = static_cast<int>((x2 - mbsize.d_view(m).x2min)/
                              mbsize.d_view(m).dx2) + js;
    int kp = static_cast<int>((x3 - mbsize.d_view(m).x3min)/
                              mbsize.d_view(m).dx3) + ks;

    Real x1v = LeftEdgeX(ip, nx1, mbsize.d_view(m).x1min,
                         mbsize.d_view(m).x1max);
    Real x2v = LeftEdgeX(jp, nx2, mbsize.d_view(m).x2min,
                         mbsize.d_view(m).x2max);
    Real x3v = LeftEdgeX(kp, nx3, mbsize.d_view(m).x3min,
                         mbsize.d_view(m).x3max);

    ParticleBField bf;
    bf.x = b0.x1f(m,kp,jp,ip) +
           (x1 - x1v)*(b0.x1f(m,kp,jp,ip+1) - b0.x1f(m,kp,jp,ip))/
           mbsize.d_view(m).dx1;
    bf.y = b0.x2f(m,kp,jp,ip);
    if (multi_d) {
      bf.y += (x2 - x2v)*(b0.x2f(m,kp,jp+1,ip) - b0.x2f(m,kp,jp,ip))/
              mbsize.d_view(m).dx2;
    }
    bf.z = b0.x3f(m,kp,jp,ip);
    if (three_d) {
      bf.z += (x3 - x3v)*(b0.x3f(m,kp+1,jp,ip) - b0.x3f(m,kp,jp,ip))/
              mbsize.d_view(m).dx3;
    }
    return bf;
  }
};

struct TrilinearGather {
  template<typename FaceField, typename CellField, typename BlockSize>
  KOKKOS_INLINE_FUNCTION
  static ParticleBField Gather(const FaceField &b0, const CellField &bcc,
                               const BlockSize &mbsize, int m, Real x1, Real x2,
                               Real x3, int is, int js, int ks, int nx1, int nx2,
                               int nx3, bool multi_d, bool three_d) {
    (void)b0;
    (void)nx1;
    (void)nx2;
    (void)nx3;
    ParticleVector3 bf = GatherCellCenteredB(
        bcc, mbsize, m, x1, x2, x3, is, js, ks, multi_d, three_d);
    return ParticleBField{bf.x, bf.y, bf.z};
  }
};

struct TSCGather {
  template<typename FaceField, typename CellField, typename BlockSize>
  KOKKOS_INLINE_FUNCTION
  static ParticleBField Gather(const FaceField &b0, const CellField &bcc,
                               const BlockSize &mbsize, int m, Real x1, Real x2,
                               Real x3, int is, int js, int ks, int nx1, int nx2,
                               int nx3, bool multi_d, bool three_d) {
    (void)b0;
    (void)nx1;
    (void)nx2;
    (void)nx3;
    Real wei1[3], wei2[3], wei3[3];

    Real a = (x1 - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 +
             static_cast<Real>(is);
    int ip = static_cast<int>(a);
    int i0 = ip - 1;
    Real d = a - static_cast<Real>(ip);
    wei1[0] = 0.5*(1.0 - d)*(1.0 - d);
    wei1[1] = 0.75 - (d - 0.5)*(d - 0.5);
    wei1[2] = 0.5*d*d;

    int j0 = js;
    wei2[0] = 1.0;
    wei2[1] = 0.0;
    wei2[2] = 0.0;
    if (multi_d) {
      a = (x2 - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 +
          static_cast<Real>(js);
      int jp = static_cast<int>(a);
      j0 = jp - 1;
      d = a - static_cast<Real>(jp);
      wei2[0] = 0.5*(1.0 - d)*(1.0 - d);
      wei2[1] = 0.75 - (d - 0.5)*(d - 0.5);
      wei2[2] = 0.5*d*d;
    }

    int k0 = ks;
    wei3[0] = 1.0;
    wei3[1] = 0.0;
    wei3[2] = 0.0;
    if (three_d) {
      a = (x3 - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 +
          static_cast<Real>(ks);
      int kp = static_cast<int>(a);
      k0 = kp - 1;
      d = a - static_cast<Real>(kp);
      wei3[0] = 0.5*(1.0 - d)*(1.0 - d);
      wei3[1] = 0.75 - (d - 0.5)*(d - 0.5);
      wei3[2] = 0.5*d*d;
    }

    ParticleBField bf{0.0, 0.0, 0.0};
    for (int k=0; k<3; ++k) {
      int kk = k0 + k;
      for (int j=0; j<3; ++j) {
        int jj = j0 + j;
        for (int i=0; i<3; ++i) {
          int ii = i0 + i;
          Real w = wei1[i]*wei2[j]*wei3[k];
          bf.x += w*bcc(m,IBX,kk,jj,ii);
          bf.y += w*bcc(m,IBY,kk,jj,ii);
          bf.z += w*bcc(m,IBZ,kk,jj,ii);
        }
      }
    }
    return bf;
  }
};

template<typename RealView>
KOKKOS_INLINE_FUNCTION
void FinishBorisPush(const RealView &pr, int p, Real dt, bool multi_d, bool three_d,
                     Real x1_old, Real x2_old, Real x3_old, Real x1, Real x2,
                     Real x3, const ParticleBField &bf) {
  pr(IPBX,p) = bf.x;
  pr(IPBY,p) = bf.y;
  pr(IPBZ,p) = bf.z;

  Real bmag = Kokkos::sqrt(bf.x*bf.x + bf.y*bf.y + bf.z*bf.z);
  if (pr(IPM,p) < 0.0) {
    if (bmag > 0.0) {
      pr(IPVX,p) = bf.x/bmag;
      pr(IPVY,p) = bf.y/bmag;
      pr(IPVZ,p) = bf.z/bmag;
    }
  } else {
    Real tx = 0.5*dt*bf.x/pr(IPM,p);
    Real ty = 0.5*dt*bf.y/pr(IPM,p);
    Real tz = 0.5*dt*bf.z/pr(IPM,p);
    Real tsq = tx*tx + ty*ty + tz*tz;

    Real v1 = pr(IPVX,p);
    Real v2 = pr(IPVY,p);
    Real v3 = pr(IPVZ,p);
    Real vp1 = v1 + v2*tz - v3*ty;
    Real vp2 = v2 + v3*tx - v1*tz;
    Real vp3 = v3 + v1*ty - v2*tx;
    Real sfac = 2.0/(1.0 + tsq);
    Real sx = tx*sfac;
    Real sy = ty*sfac;
    Real sz = tz*sfac;

    pr(IPVX,p) = v1 + vp2*sz - vp3*sy;
    pr(IPVY,p) = v2 + vp3*sx - vp1*sz;
    pr(IPVZ,p) = v3 + vp1*sy - vp2*sx;
  }

  pr(IPX,p) = x1 + 0.5*dt*pr(IPVX,p);
  if (multi_d) {pr(IPY,p) = x2 + 0.5*dt*pr(IPVY,p);}
  if (three_d) {pr(IPZ,p) = x3 + 0.5*dt*pr(IPVZ,p);}

  Real dx1 = pr(IPX,p) - x1_old;
  Real dx2 = pr(IPY,p) - x2_old;
  Real dx3 = pr(IPZ,p) - x3_old;
  pr(IPDX,p) += dx1;
  pr(IPDY,p) += dx2;
  pr(IPDZ,p) += dx3;
  if (bmag > 0.0) {
    pr(IPDB,p) += (dx1*bf.x + dx2*bf.y + dx3*bf.z)/bmag;
  }
}

template<typename GatherPolicy>
void RunBorisGather(const std::string &label, MeshBlockPack *pmy_pack,
                    DvceArray2D<Real> &pr, DvceArray2D<int> &pi, int npart,
                    int gids, int is, int js, int ks, int nx1, int nx2, int nx3,
                    bool multi_d, bool three_d, Real dt) {
  auto mbsize = pmy_pack->pmb->mb_size;
  auto b0 = pmy_pack->pmhd->b0;
  auto bcc = pmy_pack->pmhd->bcc0;

  par_for(label,DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;

    Real x1_old = pr(IPX,p);
    Real x2_old = pr(IPY,p);
    Real x3_old = pr(IPZ,p);

    Real x1 = x1_old + 0.5*dt*pr(IPVX,p);
    Real x2 = x2_old;
    Real x3 = x3_old;
    if (multi_d) {x2 += 0.5*dt*pr(IPVY,p);}
    if (three_d) {x3 += 0.5*dt*pr(IPVZ,p);}

    ParticleBField bf = GatherPolicy::Gather(b0, bcc, mbsize, m, x1, x2, x3,
                                             is, js, ks, nx1, nx2, nx3,
                                             multi_d, three_d);
    FinishBorisPush(pr, p, dt, multi_d, three_d, x1_old, x2_old, x3_old,
                    x1, x2, x3, bf);
  });
}

void RunDriftStep(MeshBlockPack *pmy_pack, DvceArray2D<Real> &pr, int npart,
                  bool multi_d, bool three_d, Real dt) {
  (void)pmy_pack;
  par_for("part_drift",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    Real x1_old = pr(IPX,p);
    Real x2_old = pr(IPY,p);
    Real x3_old = pr(IPZ,p);
    pr(IPX,p) += 0.5*dt*pr(IPVX,p);

    if (multi_d) {
      pr(IPY,p) += 0.5*dt*pr(IPVY,p);
    }

    if (three_d) {
      pr(IPZ,p) += 0.5*dt*pr(IPVZ,p);
    }

    pr(IPDX,p) += pr(IPX,p) - x1_old;
    if (multi_d) {pr(IPDY,p) += pr(IPY,p) - x2_old;}
    if (three_d) {pr(IPDZ,p) += pr(IPZ,p) - x3_old;}
  });
}

[[noreturn]] void FatalRelativisticPush(const int status) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Particle pusher 'relativistic_hc' failed with status "
            << status << std::endl;
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#endif
  std::exit(EXIT_FAILURE);
}

void RunRelativisticPrescribedStep(DvceArray2D<Real> &pr, int npart, Real c_model,
                                   Real dt) {
  DvceArray1D<int> failure("part_relativistic_push_failure", 1);
  Kokkos::deep_copy(failure, 0);
  par_for("part_relativistic_hc_prescribed",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    relativistic::Vector3 w_old{pr(IPWX,p), pr(IPWY,p), pr(IPWZ,p)};
    relativistic::Vector3 velocity_old{0.0, 0.0, 0.0};
    Real gamma_old = 0.0;
    auto state_status =
        relativistic::VelocityFromW(w_old, c_model, velocity_old, gamma_old);
    if (state_status != relativistic::StateStatus::success) {
      Kokkos::atomic_max(&failure(0), static_cast<int>(state_status));
      return;
    }

    relativistic::PushResult result;
    auto push_status = relativistic::HigueraCaryPush(
        w_old, {pr(IPCEX,p), pr(IPCEY,p), pr(IPCEZ,p)},
        {pr(IPBX,p), pr(IPBY,p), pr(IPBZ,p)}, pr(IPALPHA,p), c_model, dt, result);
    if (push_status != relativistic::PushStatus::success) {
      Kokkos::atomic_max(&failure(0), 100 + static_cast<int>(push_status));
      return;
    }

    Real half_dt = 0.0;
    relativistic::Vector3 x_old{pr(IPX,p), pr(IPY,p), pr(IPZ,p)};
    relativistic::Vector3 old_half_drift{0.0, 0.0, 0.0};
    relativistic::Vector3 x_mid{0.0, 0.0, 0.0};
    relativistic::Vector3 new_half_drift{0.0, 0.0, 0.0};
    relativistic::Vector3 x_new{0.0, 0.0, 0.0};
    relativistic::Vector3 negative_x_old{0.0, 0.0, 0.0};
    relativistic::Vector3 displacement{0.0, 0.0, 0.0};
    relativistic::Vector3 displacement_old{pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p)};
    relativistic::Vector3 displacement_new{0.0, 0.0, 0.0};
    Real work_new = 0.0;
    if (!relativistic::CheckedProduct(0.5, dt, half_dt) ||
        !relativistic::CheckedScale(half_dt, velocity_old, old_half_drift) ||
        !relativistic::CheckedAdd(x_old, old_half_drift, x_mid) ||
        !relativistic::CheckedScale(half_dt, result.velocity, new_half_drift) ||
        !relativistic::CheckedAdd(x_mid, new_half_drift, x_new) ||
        !relativistic::CheckedScale(-1.0, x_old, negative_x_old) ||
        !relativistic::CheckedAdd(x_new, negative_x_old, displacement) ||
        !relativistic::CheckedAdd(displacement_old, displacement, displacement_new) ||
        !relativistic::CheckedScalarAdd(pr(IPWORK,p), result.work, work_new)) {
      Kokkos::atomic_max(&failure(0), 200);
      return;
    }
    relativistic::Vector3 b{pr(IPBX,p), pr(IPBY,p), pr(IPBZ,p)};
    Real bmag = relativistic::ScaledNorm3(b);
    Real db_new = pr(IPDB,p);
    if (!relativistic::IsFinite(bmag) || !relativistic::IsFinite(db_new)) {
      Kokkos::atomic_max(&failure(0), 201);
      return;
    }
    if (bmag > 0.0) {
      Real displacement_dot_b = 0.0;
      Real db_increment = 0.0;
      if (!relativistic::CheckedDot(displacement, b, displacement_dot_b)) {
        Kokkos::atomic_max(&failure(0), 202);
        return;
      }
      db_increment = displacement_dot_b/bmag;
      if (!relativistic::CheckedScalarAdd(pr(IPDB,p), db_increment, db_new)) {
        Kokkos::atomic_max(&failure(0), 203);
        return;
      }
    }

    state_status =
        relativistic::StoreAuthoritativeWAndSynchronizeVelocityShadow(
            pr, p, result.w, c_model);
    if (state_status != relativistic::StateStatus::success) {
      Kokkos::atomic_max(&failure(0), 300 + static_cast<int>(state_status));
      return;
    }
    pr(IPWORK,p) = work_new;
    pr(IPX,p) = x_new.x;
    pr(IPY,p) = x_new.y;
    pr(IPZ,p) = x_new.z;
    pr(IPDX,p) = displacement_new.x;
    pr(IPDY,p) = displacement_new.y;
    pr(IPDZ,p) = displacement_new.z;
    pr(IPDB,p) = db_new;
  });
  auto failure_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), failure);
  if (failure_h(0) != 0) {FatalRelativisticPush(failure_h(0));}
}

void RunRelativisticMHDIdealStep(MeshBlockPack *pmy_pack, DvceArray2D<Real> &pr,
                                 DvceArray2D<int> &pi, int npart, int gids, int is,
                                 int js, int ks, bool multi_d, bool three_d,
                                 Real c_model, Real dt) {
  auto mbsize = pmy_pack->pmb->mb_size;
  auto w0 = pmy_pack->pmhd->w0;
  auto bcc0 = pmy_pack->pmhd->bcc0;
  int nmb = pmy_pack->nmb_thispack;
  int nx1 = pmy_pack->pmesh->mb_indcs.nx1;
  int nx2 = pmy_pack->pmesh->mb_indcs.nx2;
  int nx3 = pmy_pack->pmesh->mb_indcs.nx3;
  int ng = pmy_pack->pmesh->mb_indcs.ng;
  DvceArray1D<int> failure("part_relativistic_mhd_ideal_push_failure", 1);
  Kokkos::deep_copy(failure, 0);
  par_for("part_relativistic_hc_mhd_ideal",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    relativistic::Vector3 w_old{pr(IPWX,p), pr(IPWY,p), pr(IPWZ,p)};
    relativistic::Vector3 velocity_old{0.0, 0.0, 0.0};
    Real gamma_old = 0.0;
    auto state_status =
        relativistic::VelocityFromW(w_old, c_model, velocity_old, gamma_old);
    if (state_status != relativistic::StateStatus::success) {
      Kokkos::atomic_max(&failure(0), static_cast<int>(state_status));
      return;
    }

    Real half_dt = 0.0;
    relativistic::Vector3 x_old{pr(IPX,p), pr(IPY,p), pr(IPZ,p)};
    relativistic::Vector3 old_half_drift{0.0, 0.0, 0.0};
    relativistic::Vector3 x_mid{0.0, 0.0, 0.0};
    if (!relativistic::CheckedProduct(0.5, dt, half_dt) ||
        !relativistic::CheckedScale(half_dt, velocity_old, old_half_drift) ||
        !relativistic::CheckedAdd(x_old, old_half_drift, x_mid)) {
      Kokkos::atomic_max(&failure(0), 200);
      return;
    }

    int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      Kokkos::atomic_max(&failure(0), 207);
      return;
    }
    if (!IsCellCenteredTrilinearCoordinateSupported(
            x_mid.x, mbsize.d_view(m).x1min, mbsize.d_view(m).dx1, is, nx1,
            ng) ||
        (multi_d && !IsCellCenteredTrilinearCoordinateSupported(
                        x_mid.y, mbsize.d_view(m).x2min, mbsize.d_view(m).dx2,
                        js, nx2, ng)) ||
        (three_d && !IsCellCenteredTrilinearCoordinateSupported(
                       x_mid.z, mbsize.d_view(m).x3min, mbsize.d_view(m).dx3,
                       ks, nx3, ng))) {
      Kokkos::atomic_max(&failure(0), 208);
      return;
    }
    auto stencil = ConstructCellCenteredTrilinearStencil(
        mbsize, m, x_mid.x, x_mid.y, x_mid.z, is, js, ks, multi_d, three_d);
    if (stencil.i0 < 0 || stencil.i0 + 1 >= nx1 + 2*ng ||
        stencil.j0 < 0 || stencil.j0 + stencil.nj - 1 >= nx2 + 2*ng ||
        stencil.k0 < 0 || stencil.k0 + stencil.nk - 1 >= nx3 + 2*ng) {
      Kokkos::atomic_max(&failure(0), 208);
      return;
    }
    NewtonianCellCenteredFields fields;
    fields.u = stencil.Gather(w0, m, IVX, IVY, IVZ);
    fields.b = GatherCellCenteredB(stencil, bcc0, m);
    relativistic::Vector3 fluid_velocity{
        fields.u.x, fields.u.y, fields.u.z};
    relativistic::Vector3 b{fields.b.x, fields.b.y, fields.b.z};
    Real fluid_speed = relativistic::ScaledNorm3(fluid_velocity);
    if (!relativistic::IsFinite(fluid_velocity) || !relativistic::IsFinite(b) ||
        !relativistic::IsFinite(fluid_speed) || fluid_speed >= c_model) {
      Kokkos::atomic_max(&failure(0), 201);
      return;
    }
    relativistic::Vector3 uxB{0.0, 0.0, 0.0};
    relativistic::Vector3 cE{0.0, 0.0, 0.0};
    if (!relativistic::CheckedCross(fluid_velocity, b, uxB) ||
        !relativistic::CheckedScale(-1.0, uxB, cE)) {
      Kokkos::atomic_max(&failure(0), 202);
      return;
    }

    relativistic::PushResult result;
    auto push_status = relativistic::HigueraCaryPush(
        w_old, cE, b, pr(IPALPHA,p), c_model, dt, result);
    if (push_status != relativistic::PushStatus::success) {
      Kokkos::atomic_max(&failure(0), 100 + static_cast<int>(push_status));
      return;
    }

    relativistic::Vector3 new_half_drift{0.0, 0.0, 0.0};
    relativistic::Vector3 x_new{0.0, 0.0, 0.0};
    relativistic::Vector3 negative_x_old{0.0, 0.0, 0.0};
    relativistic::Vector3 displacement{0.0, 0.0, 0.0};
    relativistic::Vector3 displacement_old{pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p)};
    relativistic::Vector3 displacement_new{0.0, 0.0, 0.0};
    Real work_new = 0.0;
    if (!relativistic::CheckedScale(half_dt, result.velocity, new_half_drift) ||
        !relativistic::CheckedAdd(x_mid, new_half_drift, x_new) ||
        !relativistic::CheckedScale(-1.0, x_old, negative_x_old) ||
        !relativistic::CheckedAdd(x_new, negative_x_old, displacement) ||
        !relativistic::CheckedAdd(displacement_old, displacement, displacement_new) ||
        !relativistic::CheckedScalarAdd(pr(IPWORK,p), result.work, work_new)) {
      Kokkos::atomic_max(&failure(0), 203);
      return;
    }
    Real bmag = relativistic::ScaledNorm3(b);
    Real db_new = pr(IPDB,p);
    if (!relativistic::IsFinite(bmag) || !relativistic::IsFinite(db_new)) {
      Kokkos::atomic_max(&failure(0), 204);
      return;
    }
    if (bmag > 0.0) {
      Real displacement_dot_b = 0.0;
      Real db_increment = 0.0;
      if (!relativistic::CheckedDot(displacement, b, displacement_dot_b)) {
        Kokkos::atomic_max(&failure(0), 205);
        return;
      }
      db_increment = displacement_dot_b/bmag;
      if (!relativistic::CheckedScalarAdd(pr(IPDB,p), db_increment, db_new)) {
        Kokkos::atomic_max(&failure(0), 206);
        return;
      }
    }

    state_status =
        relativistic::StoreAuthoritativeWAndSynchronizeVelocityShadow(
            pr, p, result.w, c_model);
    if (state_status != relativistic::StateStatus::success) {
      Kokkos::atomic_max(&failure(0), 300 + static_cast<int>(state_status));
      return;
    }
    pr(IPBX,p) = b.x;
    pr(IPBY,p) = b.y;
    pr(IPBZ,p) = b.z;
    pr(IPCEX,p) = cE.x;
    pr(IPCEY,p) = cE.y;
    pr(IPCEZ,p) = cE.z;
    pr(IPWORK,p) = work_new;
    pr(IPX,p) = x_new.x;
    pr(IPY,p) = x_new.y;
    pr(IPZ,p) = x_new.z;
    pr(IPDX,p) = displacement_new.x;
    pr(IPDY,p) = displacement_new.y;
    pr(IPDZ,p) = displacement_new.z;
    pr(IPDB,p) = db_new;
  });
  auto failure_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), failure);
  if (failure_h(0) != 0) {FatalRelativisticPush(failure_h(0));}
}

struct RelativisticEnvelope {
  Real dx_min;
  Real block_length_min;
  Real b_max;
  Real e_max;
  Real alpha_max;
  Real w_min;
};

[[noreturn]] void FatalRelativisticBound(const std::string &message) {
  if (global_variable::my_rank == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << message << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#endif
  std::exit(EXIT_FAILURE);
}

RelativisticEnvelope ComputeRelativisticEnvelope(MeshBlockPack *pmy_pack,
                                                 DvceArray2D<Real> &pr,
                                                 int npart,
                                                 RelativisticFieldSource source,
                                                 Real alpha_s,
                                                 Real c_model) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const int nmb = pmy_pack->nmb_thispack;
  Real dx_min = std::numeric_limits<Real>::max();
  Real block_length_min = std::numeric_limits<Real>::max();
  for (int m=0; m<nmb; ++m) {
    dx_min = std::min(dx_min, mbsize.h_view(m).dx1);
    dx_min = std::min(dx_min, mbsize.h_view(m).dx2);
    dx_min = std::min(dx_min, mbsize.h_view(m).dx3);
    block_length_min = std::min(block_length_min,
        mbsize.h_view(m).x1max - mbsize.h_view(m).x1min);
    block_length_min = std::min(block_length_min,
        mbsize.h_view(m).x2max - mbsize.h_view(m).x2min);
    block_length_min = std::min(block_length_min,
        mbsize.h_view(m).x3max - mbsize.h_view(m).x3min);
  }

  Real u_max = 0.0;
  Real b_max = 0.0;
  Real e_max = 0.0;
  int invalid_fields = 0;
  if (source == RelativisticFieldSource::mhd_ideal) {
    auto w0 = pmy_pack->pmhd->w0;
    auto bcc0 = pmy_pack->pmhd->bcc0;
    const int ni = indcs.nx1;
    const int nj = indcs.nx2;
    const int nk = indcs.nx3;
    const int nmkji = nmb*nk*nj*ni;
    Kokkos::parallel_reduce("part_relativistic_grid_envelope",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, Real &local_u_max, Real &local_b_max,
                    int &local_invalid) {
        int i = indcs.is + idx % ni;
        int j = indcs.js + (idx/ni) % nj;
        int k = indcs.ks + (idx/(ni*nj)) % nk;
        int m = idx/(ni*nj*nk);
        relativistic::Vector3 u{w0(m,IVX,k,j,i), w0(m,IVY,k,j,i),
                                w0(m,IVZ,k,j,i)};
        relativistic::Vector3 b{bcc0(m,IBX,k,j,i), bcc0(m,IBY,k,j,i),
                                bcc0(m,IBZ,k,j,i)};
        Real umag = relativistic::ScaledNorm3(u);
        Real bmag = relativistic::ScaledNorm3(b);
        if (!relativistic::IsFinite(u) || !relativistic::IsFinite(b) ||
            !relativistic::IsFinite(umag) || !relativistic::IsFinite(bmag)) {
          local_invalid = 1;
          return;
        }
        local_u_max = Kokkos::max(local_u_max, umag);
        local_b_max = Kokkos::max(local_b_max, bmag);
      }, Kokkos::Max<Real>(u_max), Kokkos::Max<Real>(b_max),
         Kokkos::Max<int>(invalid_fields));
  }

  Real alpha_max = std::abs(alpha_s);
  Real w_min = std::numeric_limits<Real>::max();
  Real particle_b_max = 0.0;
  int invalid_particles = 0;
  if (npart > 0) {
    bool prescribed = (source == RelativisticFieldSource::prescribed_test);
    Kokkos::parallel_reduce("part_relativistic_particle_envelope",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, npart),
      KOKKOS_LAMBDA(const int p, Real &local_w_min, Real &local_b_max,
                    Real &local_e_max, Real &local_alpha_max,
                    int &local_invalid) {
        relativistic::Vector3 w{pr(IPWX,p), pr(IPWY,p), pr(IPWZ,p)};
        Real wmag = relativistic::ScaledNorm3(w);
        Real alpha = Kokkos::abs(pr(IPALPHA,p));
        if (!relativistic::IsFinite(w) || !relativistic::IsFinite(wmag) ||
            !relativistic::IsFinite(alpha) || alpha <= 0.0) {
          local_invalid = 1;
          return;
        }
        local_w_min = Kokkos::min(local_w_min, wmag);
        local_alpha_max = Kokkos::max(local_alpha_max, alpha);
        if (prescribed) {
          relativistic::Vector3 b{pr(IPBX,p), pr(IPBY,p), pr(IPBZ,p)};
          relativistic::Vector3 ce{pr(IPCEX,p), pr(IPCEY,p), pr(IPCEZ,p)};
          Real bmag = relativistic::ScaledNorm3(b);
          Real emag = relativistic::ScaledNorm3(ce);
          if (!relativistic::IsFinite(b) || !relativistic::IsFinite(ce) ||
              !relativistic::IsFinite(bmag) || !relativistic::IsFinite(emag)) {
            local_invalid = 1;
            return;
          }
          local_b_max = Kokkos::max(local_b_max, bmag);
          local_e_max = Kokkos::max(local_e_max, emag);
        }
      }, Kokkos::Min<Real>(w_min), Kokkos::Max<Real>(particle_b_max),
         Kokkos::Max<Real>(e_max), Kokkos::Max<Real>(alpha_max),
         Kokkos::Max<int>(invalid_particles));
  }
  b_max = std::max(b_max, particle_b_max);

  int particle_count = npart;
#if MPI_PARALLEL_ENABLED
  Real maxima[4] = {u_max, b_max, e_max, alpha_max};
  Real minima[3] = {dx_min, block_length_min, w_min};
  int invalid = std::max(invalid_fields, invalid_particles);
  MPI_Allreduce(MPI_IN_PLACE, maxima, 4, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, minima, 3, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &particle_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  u_max = maxima[0];
  b_max = maxima[1];
  e_max = maxima[2];
  alpha_max = maxima[3];
  dx_min = minima[0];
  block_length_min = minima[1];
  w_min = minima[2];
  invalid_fields = invalid;
  invalid_particles = invalid;
#endif
  if (particle_count == 0) {w_min = 0.0;}
  if (source == RelativisticFieldSource::mhd_ideal) {
    if (!relativistic::CheckedProduct(u_max, b_max, e_max)) {
      FatalRelativisticBound("relativistic MHD envelope overflowed Umax*Bmax");
    }
    if (u_max >= c_model) {
      FatalRelativisticBound(
          "Particle pusher 'relativistic_hc' failed with status 201: "
          "sampled MHD envelope has |u_fluid| >= c_model");
    }
  }
  if (invalid_fields > 0 || invalid_particles > 0 ||
      !std::isfinite(dx_min) || dx_min <= 0.0 ||
      !std::isfinite(block_length_min) || block_length_min <= 0.0 ||
      !std::isfinite(b_max) || b_max < 0.0 ||
      !std::isfinite(e_max) || e_max < 0.0 ||
      !std::isfinite(alpha_max) || alpha_max <= 0.0 ||
      !std::isfinite(w_min) || w_min < 0.0) {
    FatalRelativisticBound(
        "invalid relativistic particle timestep envelope:"
        " invalid_fields=" + std::to_string(invalid_fields) +
        " invalid_particles=" + std::to_string(invalid_particles) +
        " dx_min=" + std::to_string(dx_min) +
        " block_length_min=" + std::to_string(block_length_min) +
        " b_max=" + std::to_string(b_max) +
        " e_max=" + std::to_string(e_max) +
        " alpha_max=" + std::to_string(alpha_max) +
        " w_min=" + std::to_string(w_min) +
        " u_max=" + std::to_string(u_max) +
        " c_model=" + std::to_string(c_model));
  }
  return RelativisticEnvelope{
      dx_min, block_length_min, b_max, e_max, alpha_max, w_min};
}

Real RelativisticOuterTimestepBound(const RelativisticEnvelope &envelope,
                                   bool subcycle, int subcycle_max_steps,
                                   Real cell_fraction, Real block_fraction,
                                   Real gyro_angle_max, Real electric_kick_max,
                                   Real c_model) {
  Real alpha_e = 0.0;
  Real alpha_b = 0.0;
  Real outer_bound = 0.0;
  Real bound = std::min(cell_fraction*envelope.dx_min/c_model,
                        block_fraction*envelope.block_length_min/c_model);
  if (envelope.e_max > 0.0) {
    if (!relativistic::CheckedProduct(envelope.alpha_max, envelope.e_max,
                                      alpha_e)) {
      FatalRelativisticBound("relativistic outer electric bound overflowed");
    }
    bound = std::min(bound, electric_kick_max/alpha_e);
  }
  if (envelope.b_max > 0.0) {
    if (!relativistic::CheckedProduct(envelope.alpha_max, envelope.b_max,
                                      alpha_b)) {
      FatalRelativisticBound("relativistic outer gyro bound overflowed");
    }
    bound = std::min(bound, gyro_angle_max/alpha_b);
  }
  Real ncap = subcycle ? static_cast<Real>(subcycle_max_steps) : 1.0;
  if (!relativistic::CheckedProduct(ncap, bound, outer_bound) ||
      outer_bound <= 0.0) {
    FatalRelativisticBound("invalid relativistic outer particle timestep bound");
  }
  return outer_bound;
}

} // namespace

//----------------------------------------------------------------------------------------
// ComputeSubcycleSteps()
// Determine the global number of particle substeps required by active motion constraints.

int Particles::ComputeSubcycleSteps(Real dt, int &cell_steps, int &block_steps,
                                    int &gyro_steps) {
  cell_steps = 1;
  block_steps = 1;
  gyro_steps = 1;
  if (!subcycle) {return 1;}

  if (nprtcl_thispack > 0) {
    auto &pr = prtcl_rdata;
    auto &pi = prtcl_idata;
    auto &mbsize = pmy_pack->pmb->mb_size;
    int gids = pmy_pack->gids;
    int npart = nprtcl_thispack;
    bool multi_d = pmy_pack->pmesh->multi_d;
    bool three_d = pmy_pack->pmesh->three_d;
    Real abs_dt_half = 0.5*((dt >= 0.0) ? dt : -dt);
    Real abs_dt = ((dt >= 0.0) ? dt : -dt);
    Real cell_fraction = subcycle_cell_fraction;
    Real block_fraction = subcycle_meshblock_fraction;
    Real gyro_fraction = subcycle_gyro_fraction;
    bool boris_pusher = (pusher == ParticlesPusher::boris);

    Kokkos::parallel_reduce("particle_subcycle_steps",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, npart),
      KOKKOS_LAMBDA(const int p, int &cell_max, int &block_max, int &gyro_max) {
        int m = pi(PGID,p) - gids;
        Real vx = pr(IPVX,p);
        Real vy = multi_d ? pr(IPVY,p) : 0.0;
        Real vz = three_d ? pr(IPVZ,p) : 0.0;

        Real cell_ratio = abs_dt_half*Kokkos::abs(vx)/
                          (cell_fraction*mbsize.d_view(m).dx1);
        if (multi_d) {
          cell_ratio = Kokkos::max(cell_ratio, abs_dt_half*Kokkos::abs(vy)/
                                   (cell_fraction*mbsize.d_view(m).dx2));
        }
        if (three_d) {
          cell_ratio = Kokkos::max(cell_ratio, abs_dt_half*Kokkos::abs(vz)/
                                   (cell_fraction*mbsize.d_view(m).dx3));
        }
        cell_max = Kokkos::max(cell_max, StepsForRatio(cell_ratio));

        Real lx1 = mbsize.d_view(m).x1max - mbsize.d_view(m).x1min;
        Real block_ratio = abs_dt_half*Kokkos::abs(vx)/(block_fraction*lx1);
        if (multi_d) {
          Real lx2 = mbsize.d_view(m).x2max - mbsize.d_view(m).x2min;
          block_ratio = Kokkos::max(block_ratio, abs_dt_half*Kokkos::abs(vy)/
                                    (block_fraction*lx2));
        }
        if (three_d) {
          Real lx3 = mbsize.d_view(m).x3max - mbsize.d_view(m).x3min;
          block_ratio = Kokkos::max(block_ratio, abs_dt_half*Kokkos::abs(vz)/
                                    (block_fraction*lx3));
        }
        block_max = Kokkos::max(block_max, StepsForRatio(block_ratio));

        if (boris_pusher && pr(IPM,p) > 0.0) {
          Real bmag = Kokkos::sqrt(pr(IPBX,p)*pr(IPBX,p) + pr(IPBY,p)*pr(IPBY,p) +
                                   pr(IPBZ,p)*pr(IPBZ,p));
          Real gyro_ratio = abs_dt*bmag/(pr(IPM,p)*gyro_fraction);
          gyro_max = Kokkos::max(gyro_max, StepsForRatio(gyro_ratio));
        }
      }, Kokkos::Max<int>(cell_steps), Kokkos::Max<int>(block_steps),
      Kokkos::Max<int>(gyro_steps));
  }

#if MPI_PARALLEL_ENABLED
  int local_steps[3] = {cell_steps, block_steps, gyro_steps};
  int global_steps[3] = {1, 1, 1};
  MPI_Allreduce(local_steps, global_steps, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  cell_steps = global_steps[0];
  block_steps = global_steps[1];
  gyro_steps = global_steps[2];
#endif

  int nsub = std::max(cell_steps, std::max(block_steps, gyro_steps));
  if (nsub > subcycle_max_steps) {
    if (subcycle_strict) {
      if (global_variable::my_rank == 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle subcycling requires " << nsub
                  << " substeps, but <particles>/subcycle_max_steps = "
                  << subcycle_max_steps << std::endl
                  << "  cell_steps=" << cell_steps
                  << " block_steps=" << block_steps
                  << " gyro_steps=" << gyro_steps << std::endl;
      }
      std::exit(EXIT_FAILURE);
    }
    nsub = subcycle_max_steps;
  }
  return std::max(1, nsub);
}

//----------------------------------------------------------------------------------------
// ComputeRelativisticSubcycleSteps()
// Use one conservative MPI-global schedule for the full frozen-field outer step.

int Particles::ComputeRelativisticSubcycleSteps(Real dt, int &cell_steps,
                                                int &block_steps, int &gyro_steps,
                                                int &electric_steps) {
  cell_steps = 1;
  block_steps = 1;
  gyro_steps = 1;
  electric_steps = 1;
  if (!subcycle) {return 1;}

  const auto envelope = ComputeRelativisticEnvelope(
      pmy_pack, prtcl_rdata, nprtcl_thispack, relativistic_field_source,
      alpha_s, c_model);
  const Real abs_dt = std::abs(dt);
  Real c_dt = 0.0;
  if (!relativistic::CheckedProduct(c_model, abs_dt, c_dt)) {
    FatalRelativisticBound("relativistic crossing subcycle request overflowed");
  }
  cell_steps = StepsForRatio(
      c_dt/(subcycle_cell_fraction*envelope.dx_min));
  block_steps = StepsForRatio(
      c_dt/(subcycle_meshblock_fraction*envelope.block_length_min));
  Real alpha_e = 0.0;
  Real kick = 0.0;
  if (!relativistic::CheckedProduct(envelope.alpha_max, envelope.e_max, alpha_e) ||
      !relativistic::CheckedProduct(alpha_e, abs_dt, kick)) {
    FatalRelativisticBound("relativistic electric-kick subcycle request overflowed");
  }
  Real w_floor = std::max(0.0, envelope.w_min - kick);
  Real gamma_floor = std::hypot(1.0, w_floor/c_model);
  if (!std::isfinite(gamma_floor) || gamma_floor < 1.0) {
    FatalRelativisticBound("invalid conservative relativistic gamma floor");
  }
  if (envelope.b_max > 0.0) {
    Real alpha_b = 0.0;
    Real gyro_request = 0.0;
    if (!relativistic::CheckedProduct(
            envelope.alpha_max, envelope.b_max, alpha_b) ||
        !relativistic::CheckedProduct(alpha_b, abs_dt, gyro_request)) {
      FatalRelativisticBound("relativistic gyro subcycle request overflowed");
    }
    gyro_steps = StepsForRatio(
        gyro_request/(subcycle_gyro_fraction*gamma_floor));
  }
  if (envelope.e_max > 0.0) {
    electric_steps = StepsForRatio(kick/subcycle_electric_kick_max);
  }

  int nsub = std::max(std::max(cell_steps, block_steps),
                      std::max(gyro_steps, electric_steps));
  if (nsub > subcycle_max_steps) {
    FatalRelativisticBound(
        "Relativistic particle subcycling requires " + std::to_string(nsub) +
        " substeps, but <particles>/subcycle_max_steps = " +
        std::to_string(subcycle_max_steps) + "; cap clipping is unsupported");
  }
  return std::max(1, nsub);
}

//----------------------------------------------------------------------------------------
// LogSubcycle()
// Optional diagnostic reporting the active particle subcycle constraint.

void Particles::LogSubcycle(int nsub, int cell_steps, int block_steps, int gyro_steps) {
  if (!log_performance || global_variable::my_rank != 0) {return;}
  const char *constraint = "cell";
  int active = cell_steps;
  if (block_steps > active) {
    constraint = "meshblock";
    active = block_steps;
  }
  if (gyro_steps > active) {
    constraint = "gyro";
    active = gyro_steps;
  }
  std::cout << "Particle subcycling: steps=" << nsub
            << " active_constraint=" << constraint
            << " cell_steps=" << cell_steps
            << " meshblock_steps=" << block_steps
            << " gyro_steps=" << gyro_steps << std::endl;
}

//----------------------------------------------------------------------------------------
// LogRelativisticSubcycle()
// Optional diagnostic reporting the active acceleration-aware constraint.

void Particles::LogRelativisticSubcycle(int nsub, int cell_steps, int block_steps,
                                        int gyro_steps, int electric_steps) {
  if (!log_performance || global_variable::my_rank != 0) {return;}
  const char *constraint = "cell";
  int active = cell_steps;
  if (block_steps > active) {
    constraint = "meshblock";
    active = block_steps;
  }
  if (gyro_steps > active) {
    constraint = "gyro";
    active = gyro_steps;
  }
  if (electric_steps > active) {
    constraint = "electric_kick";
  }
  std::cout << "Particle relativistic subcycling: steps=" << nsub
            << " active_constraint=" << constraint
            << " cell_steps=" << cell_steps
            << " meshblock_steps=" << block_steps
            << " gyro_steps=" << gyro_steps
            << " electric_steps=" << electric_steps << std::endl;
}

//----------------------------------------------------------------------------------------
// MarkRelativisticTimestepBoundDirty()

void Particles::MarkRelativisticTimestepBoundDirty() {
  if (pusher == ParticlesPusher::relativistic_hc &&
      !relativistic_timestep_bound_override) {
    relativistic_timestep_bound_dirty = true;
  }
}

//----------------------------------------------------------------------------------------
// RefreshRelativisticTimestepBound()

void Particles::RefreshRelativisticTimestepBound() {
  if (pusher != ParticlesPusher::relativistic_hc) {return;}
  if (relativistic_timestep_bound_override) {return;}
#if MPI_PARALLEL_ENABLED
  int dirty = relativistic_timestep_bound_dirty ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &dirty, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (dirty == 0) {return;}
#else
  if (!relativistic_timestep_bound_dirty) {return;}
#endif
  const auto envelope = ComputeRelativisticEnvelope(
      pmy_pack, prtcl_rdata, nprtcl_thispack, relativistic_field_source,
      alpha_s, c_model);
  dtnew = RelativisticOuterTimestepBound(
      envelope, subcycle, subcycle_max_steps, subcycle_cell_fraction,
      subcycle_meshblock_fraction, subcycle_gyro_fraction,
      subcycle_electric_kick_max, c_model);
  relativistic_timestep_bound_dirty = false;
}

//----------------------------------------------------------------------------------------
// ExchangeAfterSubcycle()
// Update ownership between intermediate particle substeps.

void Particles::ExchangeAfterSubcycle() {
  pbval_part->SetNewPrtclGID();
#if MPI_PARALLEL_ENABLED
  pbval_part->CountSendsAndRecvs();
  pbval_part->InitPrtclRecv();
  pbval_part->PackAndSendPrtcls();
  while (pbval_part->RecvAndUnpackPrtcls() == TaskStatus::incomplete) {}
  pbval_part->ClearPrtclRecv();
  pbval_part->ClearPrtclSend();
#else
  pmy_pack->pmesh->UpdateParticleCounts();
  CheckConsistency("particle subcycle exchange");
#endif
}

//----------------------------------------------------------------------------------------
//! \fn  void Particles::ParticlesPush
//  \brief

TaskStatus Particles::Push(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &pr = prtcl_rdata;
  Real dt = pmy_pack->pmesh->dt;
  last_push_dt = dt;
  int gids = pmy_pack->gids;
  int cell_steps, block_steps, gyro_steps;
  int electric_steps = 1;
  int nsub = 1;
  if (pusher == ParticlesPusher::relativistic_hc) {
    nsub = ComputeRelativisticSubcycleSteps(
        dt, cell_steps, block_steps, gyro_steps, electric_steps);
    LogRelativisticSubcycle(
        nsub, cell_steps, block_steps, gyro_steps, electric_steps);
  } else {
    nsub = ComputeSubcycleSteps(dt, cell_steps, block_steps, gyro_steps);
    LogSubcycle(nsub, cell_steps, block_steps, gyro_steps);
  }
  last_subcycle_steps = nsub;
  Real dt_sub = dt/static_cast<Real>(nsub);

  for (int sub=0; sub<nsub; ++sub) {
    int npart = nprtcl_thispack;
    if (npart > 0) {
      switch (pusher) {
        case ParticlesPusher::drift:
          RunDriftStep(pmy_pack, prtcl_rdata, npart, multi_d, three_d, dt_sub);
          break;

        case ParticlesPusher::boris:
          switch (interpolation) {
            case ParticleInterpolation::lin_legacy:
              RunBorisGather<LinLegacyGather>("part_boris_lin", pmy_pack, prtcl_rdata,
                                              prtcl_idata, npart, gids, is, js, ks,
                                              nx1, nx2, nx3, multi_d, three_d, dt_sub);
              break;
            case ParticleInterpolation::trilinear:
              RunBorisGather<TrilinearGather>("part_boris_trilinear", pmy_pack,
                                              prtcl_rdata, prtcl_idata, npart, gids,
                                              is, js, ks, nx1, nx2, nx3, multi_d,
                                              three_d, dt_sub);
              break;
            case ParticleInterpolation::tsc:
              RunBorisGather<TSCGather>("part_boris_tsc", pmy_pack, prtcl_rdata,
                                        prtcl_idata, npart, gids, is, js, ks, nx1,
                                        nx2, nx3, multi_d, three_d, dt_sub);
              break;
          }
          break;

        case ParticlesPusher::relativistic_hc:
          switch (relativistic_field_source) {
            case RelativisticFieldSource::prescribed_test:
              RunRelativisticPrescribedStep(prtcl_rdata, npart, c_model, dt_sub);
              break;
            case RelativisticFieldSource::mhd_ideal:
              RunRelativisticMHDIdealStep(
                  pmy_pack, prtcl_rdata, prtcl_idata, npart, gids, is, js, ks,
                  multi_d, three_d, c_model, dt_sub);
              break;
          }
          break;

        default:
          break;
        }
    }
    if (sub < nsub - 1) {
      ExchangeAfterSubcycle();
      CheckMotionBounds("particle subcycle exchange");
    }
  }

  MarkRelativisticTimestepBoundDirty();
  CheckMotionBounds("particle push");
  LogPerformance("push", nprtcl_thispack, 0, 0, 0, 0);
  return TaskStatus::complete;
}
} // namespace particles
