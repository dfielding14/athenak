//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include <algorithm>
#include <cstdlib>
#include <iostream>
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
int StepsForRatio(Real ratio) {
  if (ratio <= 1.0) {return 1;}
  int n = static_cast<int>(ratio);
  if (static_cast<Real>(n) <= ratio) {++n;}
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
  int gids = pmy_pack->gids;
  int cell_steps, block_steps, gyro_steps;
  int nsub = ComputeSubcycleSteps(dt, cell_steps, block_steps, gyro_steps);
  LogSubcycle(nsub, cell_steps, block_steps, gyro_steps);
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
          RunRelativisticPrescribedStep(prtcl_rdata, npart, c_model, dt_sub);
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

  CheckMotionBounds("particle push");
  LogPerformance("push", nprtcl_thispack, 0, 0, 0, 0);
  return TaskStatus::complete;
}
} // namespace particles
