//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"
#include "units/units.hpp"

namespace particles {
namespace {

bool PushValidationEnabled() {
  static int enabled = []() -> int {
    const char *env = std::getenv("ATHENA_PIC_VALIDATE_PUSH");
    if (env == nullptr) return 0;
    return (std::atoi(env) != 0) ? 1 : 0;
  }();
  return (enabled != 0);
}

enum PushValidateCounter {
  kBadGID = 0,
  kBadSpecies = 1,
  kBadX = 2,
  kBadY = 3,
  kBadZ = 4,
  kBadVx = 5,
  kBadVy = 6,
  kBadVz = 7,
  kBadQOverM = 8,
  kNumPushValidateCounters = 9
};

KOKKOS_INLINE_FUNCTION
Real TSCWeight(const Real d) {
  const Real ad = fabs(d);
  if (ad < static_cast<Real>(0.5)) return static_cast<Real>(0.75) - ad*ad;
  if (ad < static_cast<Real>(1.5)) {
    const Real t = static_cast<Real>(1.5) - ad;
    return static_cast<Real>(0.5)*t*t;
  }
  return static_cast<Real>(0.0);
}

template <typename SizeViewType>
KOKKOS_INLINE_FUNCTION
void InterpolateLinearDevice(const SizeViewType &size_d, const DvceArray5D<Real> &bcc,
                             const DvceArray5D<Real> &w0, const bool use_no_mhd_bcc,
                             const int m, const Real x, const Real y, const Real z,
                             const int is, const int ie, const int js, const int je,
                             const int ks, const int ke, const int nx3,
                             const bool allow_2d3v, Real &Bx, Real &By,
                             Real &Bz, Real &Ux, Real &Uy, Real &Uz) {
  const bool use_bz_channel = (nx3 > 1) || allow_2d3v;
  const Real dx1 = size_d(m).dx1;
  const Real dx2 = size_d(m).dx2;
  const Real dx3 = (nx3 > 1) ? size_d(m).dx3 : static_cast<Real>(1.0);

  const Real fx = (x - size_d(m).x1min)/dx1;
  const Real fy = (y - size_d(m).x2min)/dx2;
  const Real fz = (nx3 > 1) ? (z - size_d(m).x3min)/dx3 : static_cast<Real>(0.0);

  int i0 = static_cast<int>(floor(fx)) + is;
  int j0 = static_cast<int>(floor(fy)) + js;
  int k0 = (nx3 > 1) ? static_cast<int>(floor(fz)) + ks : ks;

  i0 = (i0 < is) ? is : ((i0 > ie - 1) ? ie - 1 : i0);
  j0 = (j0 < js) ? js : ((j0 > je - 1) ? je - 1 : j0);
  k0 = (k0 < ks) ? ks : ((k0 > ke - 1) ? ke - 1 : k0);

  const int i1 = i0 + 1;
  const int j1 = j0 + 1;
  const int k1 = (nx3 > 1) ? k0 + 1 : k0;

  const Real wx = fx - floor(fx);
  const Real wy = fy - floor(fy);
  const Real wz = (nx3 > 1) ? fz - floor(fz) : static_cast<Real>(0.0);

  Bx = By = Bz = static_cast<Real>(0.0);
  Ux = Uy = Uz = static_cast<Real>(0.0);

  const int nk = (nx3 > 1) ? 2 : 1;
  for (int dk = 0; dk < nk; ++dk) {
    const Real wk = (nx3 > 1) ? ((dk == 0) ? (static_cast<Real>(1.0) - wz) : wz)
                               : static_cast<Real>(1.0);
    const int kk = (dk == 0) ? k0 : k1;
    for (int dj = 0; dj <= 1; ++dj) {
      const Real wj = (dj == 0) ? (static_cast<Real>(1.0) - wy) : wy;
      const int jj = (dj == 0) ? j0 : j1;
      for (int di = 0; di <= 1; ++di) {
        const Real wi = (di == 0) ? (static_cast<Real>(1.0) - wx) : wx;
        const int ii = (di == 0) ? i0 : i1;
        const Real w = wi*wj*wk;
        Bx += w*bcc(m, IBX, kk, jj, ii);
        By += w*bcc(m, IBY, kk, jj, ii);
        if (!use_no_mhd_bcc) {
          Ux += w*w0(m, IVX, kk, jj, ii);
          Uy += w*w0(m, IVY, kk, jj, ii);
          Uz += w*w0(m, IVZ, kk, jj, ii);
        }
        if (use_bz_channel) {
          Bz += w*bcc(m, IBZ, kk, jj, ii);
        }
      }
    }
  }
  if (!use_bz_channel) {
    Bz = static_cast<Real>(0.0);
    Uz = static_cast<Real>(0.0);
  }
}

template <typename SizeViewType>
KOKKOS_INLINE_FUNCTION
void InterpolateTSCDevice(const SizeViewType &size_d, const DvceArray5D<Real> &bcc,
                          const DvceArray5D<Real> &w0, const bool use_no_mhd_bcc,
                          const int m, const Real x, const Real y, const Real z,
                          const int is, const int ie, const int js, const int je,
                          const int ks, const int ke, const int nx3,
                          const bool allow_2d3v, Real &Bx, Real &By,
                          Real &Bz, Real &Ux, Real &Uy, Real &Uz) {
  const bool use_bz_channel = (nx3 > 1) || allow_2d3v;
  const Real dx1 = size_d(m).dx1;
  const Real dx2 = size_d(m).dx2;
  const Real dx3 = (nx3 > 1) ? size_d(m).dx3 : static_cast<Real>(1.0);

  const Real fx = (x - size_d(m).x1min)/dx1 - static_cast<Real>(0.5);
  const Real fy = (y - size_d(m).x2min)/dx2 - static_cast<Real>(0.5);
  const Real fz = (nx3 > 1)
                      ? (z - size_d(m).x3min)/dx3 - static_cast<Real>(0.5)
                      : static_cast<Real>(0.0);

  const int ic = static_cast<int>(floor(fx)) + is;
  const int jc = static_cast<int>(floor(fy)) + js;
  const int kc = (nx3 > 1) ? static_cast<int>(floor(fz)) + ks : ks;

  const Real di = fx - floor(fx);
  const Real dj = fy - floor(fy);
  const Real dk = (nx3 > 1) ? fz - floor(fz) : static_cast<Real>(0.0);

  const Real wx[3] = {TSCWeight(di + static_cast<Real>(1.0)),
                      TSCWeight(di),
                      TSCWeight(di - static_cast<Real>(1.0))};
  const Real wy[3] = {TSCWeight(dj + static_cast<Real>(1.0)),
                      TSCWeight(dj),
                      TSCWeight(dj - static_cast<Real>(1.0))};

  Real wz[3] = {static_cast<Real>(1.0), static_cast<Real>(0.0),
                static_cast<Real>(0.0)};
  int kz[3] = {kc, kc, kc};
  if (nx3 > 1) {
    wz[0] = TSCWeight(dk + static_cast<Real>(1.0));
    wz[1] = TSCWeight(dk);
    wz[2] = TSCWeight(dk - static_cast<Real>(1.0));
    kz[0] = kc - 1;
    kz[1] = kc;
    kz[2] = kc + 1;
  }

  int ix[3] = {ic - 1, ic, ic + 1};
  int iy[3] = {jc - 1, jc, jc + 1};
  for (int n = 0; n < 3; ++n) {
    ix[n] = (ix[n] < is) ? is : ((ix[n] > ie) ? ie : ix[n]);
    iy[n] = (iy[n] < js) ? js : ((iy[n] > je) ? je : iy[n]);
    if (nx3 > 1) {
      kz[n] = (kz[n] < ks) ? ks : ((kz[n] > ke) ? ke : kz[n]);
    }
  }

  Bx = By = Bz = static_cast<Real>(0.0);
  Ux = Uy = Uz = static_cast<Real>(0.0);

  const int nk = (nx3 > 1) ? 3 : 1;
  for (int kk = 0; kk < nk; ++kk) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        const Real w = wx[ii]*wy[jj]*wz[kk];
        Bx += w*bcc(m, IBX, kz[kk], iy[jj], ix[ii]);
        By += w*bcc(m, IBY, kz[kk], iy[jj], ix[ii]);
        if (!use_no_mhd_bcc) {
          Ux += w*w0(m, IVX, kz[kk], iy[jj], ix[ii]);
          Uy += w*w0(m, IVY, kz[kk], iy[jj], ix[ii]);
          Uz += w*w0(m, IVZ, kz[kk], iy[jj], ix[ii]);
        }
        if (use_bz_channel) {
          Bz += w*bcc(m, IBZ, kz[kk], iy[jj], ix[ii]);
        }
      }
    }
  }
  if (!use_bz_channel) {
    Bz = static_cast<Real>(0.0);
    Uz = static_cast<Real>(0.0);
  }
}

}  // namespace

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3, Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal, Real R200, Real rho_mean);

//----------------------------------------------------------------------------------------
//! \fn void Particles::ParticlesPush
//  \brief

TaskStatus Particles::Push(Driver *pdriver, int stage) {
  switch (pusher) {
  case ParticlesPusher::drift:
    return PushDrift(pdriver, stage);
  case ParticlesPusher::rk4_gravity:
    return PushStars(pdriver, stage);
  case ParticlesPusher::boris_lin:
  case ParticlesPusher::boris_tsc:
    return PushCosmicRays(pdriver, stage);
  default:
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushDrift
//! \brief Simple drift pusher updating positions using velocities

TaskStatus Particles::PushDrift(Driver *pdriver, int stage) {
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &pr = prtcl_rdata;
  auto dt_ = (pmy_pack->pmesh->dt);

  par_for(
      "part_update", DevExeSpace(), 0, (nprtcl_thispack - 1),
      KOKKOS_LAMBDA(const int p) {
        pr(IPX, p) += dt_ * pr(IPVX, p);

        if (multi_d) {
          pr(IPY, p) += dt_ * pr(IPVY, p);
        }

        if (three_d) {
          pr(IPZ, p) += dt_ * pr(IPVZ, p);
        }
      });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushStars
//! \brief RK4 gravity pusher for star particles

TaskStatus Particles::PushStars(Driver *pdriver, int stage) {
  auto &pr = prtcl_rdata;
  auto dt_ = (pmy_pack->pmesh->dt);

  Real G = pmy_pack->punit->grav_constant();
  Real r_s = r_scale;
  Real rho_s = rho_scale;
  Real m_g = m_gal;
  Real a_g = a_gal;
  Real z_g = z_gal;
  Real r_m = r_200;
  Real rho_m = rho_mean;
  Real h = par_grav_dx;

  par_for(
      "part_rk4_gravity", DevExeSpace(), 0, (nprtcl_thispack - 1),
      KOKKOS_LAMBDA(const int p) {
        // Current state
        Real x = pr(IPX, p);
        Real y = pr(IPY, p);
        Real z = pr(IPZ, p);
        Real vx = pr(IPVX, p);
        Real vy = pr(IPVY, p);
        Real vz = pr(IPVZ, p);

        auto compute_acc = [&](Real px, Real py, Real pz, Real &ax, Real &ay,
                               Real &az) {
          Real phi_xp =
              GravPot(px + h, py, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_xm =
              GravPot(px - h, py, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_yp =
              GravPot(px, py + h, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_ym =
              GravPot(px, py - h, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_zp =
              GravPot(px, py, pz + h, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_zm =
              GravPot(px, py, pz - h, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);

          ax = -(phi_xp - phi_xm) / (2.0 * h);
          ay = -(phi_yp - phi_ym) / (2.0 * h);
          az = -(phi_zp - phi_zm) / (2.0 * h);
        };

        // RK4 coefficients
        Real k1x, k1y, k1z, k1vx, k1vy, k1vz;
        Real k2x, k2y, k2z, k2vx, k2vy, k2vz;
        Real k3x, k3y, k3z, k3vx, k3vy, k3vz;
        Real k4x, k4y, k4z, k4vx, k4vy, k4vz;

        // K1: derivatives at current state
        k1x = vx;
        k1y = vy;
        k1z = vz;
        compute_acc(x, y, z, k1vx, k1vy, k1vz);

        // K2: derivatives at midpoint using K1
        Real x_mid1 = x + 0.5 * dt_ * k1x;
        Real y_mid1 = y + 0.5 * dt_ * k1y;
        Real z_mid1 = z + 0.5 * dt_ * k1z;
        Real vx_mid1 = vx + 0.5 * dt_ * k1vx;
        Real vy_mid1 = vy + 0.5 * dt_ * k1vy;
        Real vz_mid1 = vz + 0.5 * dt_ * k1vz;

        k2x = vx_mid1;
        k2y = vy_mid1;
        k2z = vz_mid1;
        compute_acc(x_mid1, y_mid1, z_mid1, k2vx, k2vy, k2vz);

        // K3: derivatives at midpoint using K2
        Real x_mid2 = x + 0.5 * dt_ * k2x;
        Real y_mid2 = y + 0.5 * dt_ * k2y;
        Real z_mid2 = z + 0.5 * dt_ * k2z;
        Real vx_mid2 = vx + 0.5 * dt_ * k2vx;
        Real vy_mid2 = vy + 0.5 * dt_ * k2vy;
        Real vz_mid2 = vz + 0.5 * dt_ * k2vz;

        k3x = vx_mid2;
        k3y = vy_mid2;
        k3z = vz_mid2;
        compute_acc(x_mid2, y_mid2, z_mid2, k3vx, k3vy, k3vz);

        // K4: derivatives at endpoint using K3
        Real x_end = x + dt_ * k3x;
        Real y_end = y + dt_ * k3y;
        Real z_end = z + dt_ * k3z;
        Real vx_end = vx + dt_ * k3vx;
        Real vy_end = vy + dt_ * k3vy;
        Real vz_end = vz + dt_ * k3vz;

        k4x = vx_end;
        k4y = vy_end;
        k4z = vz_end;
        compute_acc(x_end, y_end, z_end, k4vx, k4vy, k4vz);

        // Final RK4 update
        pr(IPX, p) = x + dt_ / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        pr(IPY, p) = y + dt_ / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
        pr(IPZ, p) = z + dt_ / 6.0 * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);

        pr(IPVX, p) = vx + dt_ / 6.0 * (k1vx + 2.0 * k2vx + 2.0 * k3vx + k4vx);
        pr(IPVY, p) = vy + dt_ / 6.0 * (k1vy + 2.0 * k2vy + 2.0 * k3vy + k4vy);
        pr(IPVZ, p) = vz + dt_ / 6.0 * (k1vz + 2.0 * k2vz + 2.0 * k3vz + k4vz);
      });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushCosmicRays
//  \brief Boris pusher for cosmic ray particles in electromagnetic fields

TaskStatus Particles::PushCosmicRays(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  const Real dt = pmy_pack->pmesh->dt;
  const int gids = pmy_pack->gids;
  const bool use_tsc = (pusher == ParticlesPusher::boris_tsc);

  const bool no_mhd_mode = (pic_background_mode == PICBackgroundMode::no_mhd);
  const bool has_field_carrier = (pmy_pack->pmhd != nullptr) || no_mhd_mode;
  if (!has_field_carrier) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Boris pushers require MHD fields, or "
              << "<particles>/pic_background_mode=no_mhd with "
              << "pic_no_mhd_bx/by/bz" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const bool use_no_mhd_bcc =
      (pic_background_mode == PICBackgroundMode::no_mhd) ||
      (pmy_pack->pmhd == nullptr);
  auto bcc = use_no_mhd_bcc ? pic_no_mhd_bcc0 : pmy_pack->pmhd->bcc0;
  auto w0 = (use_no_mhd_bcc || pmy_pack->pmhd == nullptr) ?
            DvceArray5D<Real>() : pmy_pack->pmhd->w0;

  // Create local copies to avoid GPU lambda capture warning
  const bool track_displacement_local = track_displacement;
  const int nx3_local = indcs.nx3;  // avoid capturing host ref
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto size_d = size.d_view;
  const bool allow_2d3v_local = (pic_enable_2d3v && (nx3_local == 1));
  const bool use_vz_component = (nx3_local > 1) || allow_2d3v_local;

  const Real qscale = deposit_qscale;
  auto mspecies = species_mass;
  const int nspecies_local = nspecies;
  const int nmb_local = pmy_pack->nmb_thispack;
  const Real inv_dt = (dt > 0.0) ? (1.0/dt) : 0.0;
  const bool expanding_box_local =
      (pic_expanding_box_mode == PICExpandingBoxMode::on);
  const Real exp_rate_x1_local = pic_expansion_rate_x1;
  const Real exp_rate_x2_local = pic_expansion_rate_x2;
  const Real exp_rate_x3_local = pic_expansion_rate_x3;

  if (PushValidationEnabled()) {
    Kokkos::View<int*> bad_counts("push_validate_bad_counts",
                                  kNumPushValidateCounters);
    Kokkos::deep_copy(bad_counts, 0);
    Kokkos::View<int> first_bad_idx("push_validate_first_bad_idx");
    Kokkos::deep_copy(first_bad_idx, -1);
    Kokkos::View<int*> first_bad_i("push_validate_first_bad_i", 2);
    Kokkos::deep_copy(first_bad_i, -1);
    Kokkos::View<Real*> first_bad_r("push_validate_first_bad_r", 7);
    Kokkos::deep_copy(first_bad_r, static_cast<Real>(0.0));

    const int npart_local = nprtcl_thispack;
    if (npart_local > 0) {
      par_for("push_validate_particles", DevExeSpace(), 0, npart_local - 1,
      KOKKOS_LAMBDA(const int p) {
        const int gid = pi(PGID, p);
        const int sp = pi(PSP, p);
        const Real x = pr(IPX, p);
        const Real y = pr(IPY, p);
        const Real z = pr(IPZ, p);
        const Real vx = pr(IPVX, p);
        const Real vy = pr(IPVY, p);
        const Real vz = pr(IPVZ, p);
        const Real qom = pr(IPM, p);

        bool bad = false;
        if (gid < gids || gid >= (gids + nmb_local)) {
          Kokkos::atomic_add(&bad_counts(kBadGID), 1);
          bad = true;
        }
        if (sp < 0 || sp >= nspecies_local) {
          Kokkos::atomic_add(&bad_counts(kBadSpecies), 1);
          bad = true;
        }
        if (!isfinite(x)) {
          Kokkos::atomic_add(&bad_counts(kBadX), 1);
          bad = true;
        }
        if (!isfinite(y)) {
          Kokkos::atomic_add(&bad_counts(kBadY), 1);
          bad = true;
        }
        if (!isfinite(z)) {
          Kokkos::atomic_add(&bad_counts(kBadZ), 1);
          bad = true;
        }
        if (!isfinite(vx)) {
          Kokkos::atomic_add(&bad_counts(kBadVx), 1);
          bad = true;
        }
        if (!isfinite(vy)) {
          Kokkos::atomic_add(&bad_counts(kBadVy), 1);
          bad = true;
        }
        if (!isfinite(vz)) {
          Kokkos::atomic_add(&bad_counts(kBadVz), 1);
          bad = true;
        }
        if (!isfinite(qom)) {
          Kokkos::atomic_add(&bad_counts(kBadQOverM), 1);
          bad = true;
        }

        if (bad) {
          const int prev = Kokkos::atomic_compare_exchange(&first_bad_idx(), -1, p);
          if (prev == -1) {
            first_bad_i(0) = gid;
            first_bad_i(1) = sp;
            first_bad_r(0) = x;
            first_bad_r(1) = y;
            first_bad_r(2) = z;
            first_bad_r(3) = vx;
            first_bad_r(4) = vy;
            first_bad_r(5) = vz;
            first_bad_r(6) = qom;
          }
        }
      });
    }
    Kokkos::fence();

    auto h_bad_counts = Kokkos::create_mirror_view_and_copy(HostMemSpace(), bad_counts);
    int local_counts[kNumPushValidateCounters];
    int global_counts[kNumPushValidateCounters];
    int local_bad_total = 0;
    for (int n = 0; n < kNumPushValidateCounters; ++n) {
      local_counts[n] = h_bad_counts(n);
      global_counts[n] = local_counts[n];
      local_bad_total += local_counts[n];
    }

#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(local_counts, global_counts, kNumPushValidateCounters,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    int global_bad_total = 0;
    for (int n = 0; n < kNumPushValidateCounters; ++n) {
      global_bad_total += global_counts[n];
    }

    if (global_bad_total > 0) {
      int h_first_bad_idx = -1;
      Kokkos::deep_copy(h_first_bad_idx, first_bad_idx);
      if (local_bad_total > 0) {
        auto h_first_bad_i =
            Kokkos::create_mirror_view_and_copy(HostMemSpace(), first_bad_i);
        auto h_first_bad_r =
            Kokkos::create_mirror_view_and_copy(HostMemSpace(), first_bad_r);
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "ATHENA_PIC_VALIDATE_PUSH detected invalid particle state."
                  << std::endl
                  << "rank=" << global_variable::my_rank
                  << " cycle=" << pmy_pack->pmesh->ncycle
                  << " time=" << pmy_pack->pmesh->time
                  << " local_counts(gid,species,x,y,z,vx,vy,vz,q_over_m)=("
                  << local_counts[kBadGID] << ","
                  << local_counts[kBadSpecies] << ","
                  << local_counts[kBadX] << ","
                  << local_counts[kBadY] << ","
                  << local_counts[kBadZ] << ","
                  << local_counts[kBadVx] << ","
                  << local_counts[kBadVy] << ","
                  << local_counts[kBadVz] << ","
                  << local_counts[kBadQOverM] << ")"
                  << std::endl;
        if (h_first_bad_idx >= 0) {
          std::cout << "first_bad_particle: p=" << h_first_bad_idx
                    << " gid=" << h_first_bad_i(0)
                    << " sp=" << h_first_bad_i(1)
                    << " x=(" << h_first_bad_r(0) << ","
                    << h_first_bad_r(1) << ","
                    << h_first_bad_r(2) << ")"
                    << " v=(" << h_first_bad_r(3) << ","
                    << h_first_bad_r(4) << ","
                    << h_first_bad_r(5) << ")"
                    << " q_over_m=" << h_first_bad_r(6)
                    << std::endl;
        }
      }
      if (global_variable::my_rank == 0) {
        std::cout << "push_validate global_counts(gid,species,x,y,z,vx,vy,vz,q_over_m)=("
                  << global_counts[kBadGID] << ","
                  << global_counts[kBadSpecies] << ","
                  << global_counts[kBadX] << ","
                  << global_counts[kBadY] << ","
                  << global_counts[kBadZ] << ","
                  << global_counts[kBadVx] << ","
                  << global_counts[kBadVy] << ","
                  << global_counts[kBadVz] << ","
                  << global_counts[kBadQOverM] << ")"
                  << std::endl;
      }
      std::exit(EXIT_FAILURE);
    }
  }

  // Midpoint E+B Boris pusher
  par_for(
      "push_cr", DevExeSpace(), 0, nprtcl_thispack - 1,
      KOKKOS_LAMBDA(const int p) {
        // Get particle properties
        int m = pi(PGID, p) - gids;
        if (m < 0 || m >= nmb_local) return;
        Real x = pr(IPX, p);
        Real y = pr(IPY, p);
        Real z = (nx3_local > 1) ? pr(IPZ, p) : 0.0;
        Real vx = pr(IPVX, p);
        Real vy = pr(IPVY, p);
        Real vz = use_vz_component ? pr(IPVZ, p) : 0.0;
        Real vx_old = vx;
        Real vy_old = vy;
        Real vz_old = vz;
        Real q_over_m = pr(IPM, p);
        int sp = pi(PSP, p);
        if (sp < 0 || sp >= nspecies_local) return;
        Real m_macro = qscale*mspecies(sp);

        if (expanding_box_local) {
          // Apply half-step expansion/compression source update around Boris.
          Real sx = exp(-0.5*dt*exp_rate_x1_local);
          Real sy = exp(-0.5*dt*exp_rate_x2_local);
          Real sz = exp(-0.5*dt*exp_rate_x3_local);
          vx *= sx;
          vy *= sy;
          if (use_vz_component) {
            vz *= sz;
          }
        }

        // Drift to midpoint with old velocity.
        Real dt_half = 0.5*dt;
        Real x_mid = x + dt_half*vx;
        Real y_mid = y + dt_half*vy;
        Real z_mid = z + dt_half*vz;

        // Interpolate midpoint B and fluid velocity for frozen-in cE = -u x B.
        Real Bx = 0.0, By = 0.0, Bz = 0.0;
        Real Ux = 0.0, Uy = 0.0, Uz = 0.0;
        if (use_tsc) {
          InterpolateTSCDevice(size_d, bcc, w0, use_no_mhd_bcc, m,
                               x_mid, y_mid, z_mid,
                               is, ie, js, je, ks, ke, nx3_local,
                               allow_2d3v_local, Bx, By, Bz, Ux, Uy, Uz);
        } else {
          InterpolateLinearDevice(size_d, bcc, w0, use_no_mhd_bcc, m,
                                  x_mid, y_mid, z_mid,
                                  is, ie, js, je, ks, ke, nx3_local,
                                  allow_2d3v_local, Bx, By, Bz, Ux, Uy, Uz);
        }
        if (!use_vz_component) {
          Bz = 0.0;
          Uz = 0.0;
        }

        Real cEx = -(Uy*Bz - Uz*By);
        Real cEy = -(Uz*Bx - Ux*Bz);
        Real cEz = -(Ux*By - Uy*Bx);
        if (!use_vz_component) cEz = 0.0;

        Real v_old_sq = vx_old*vx_old + vy_old*vy_old + vz_old*vz_old;

        // Half electric acceleration
        Real qdt_2m = q_over_m*dt_half;
        vx += qdt_2m*cEx;
        vy += qdt_2m*cEy;
        vz += qdt_2m*cEz;

        // Magnetic rotation
        Real tx = qdt_2m*Bx;
        Real ty = qdt_2m*By;
        Real tz = qdt_2m*Bz;
        Real t2 = tx*tx + ty*ty + tz*tz;
        Real sx = 2.0*tx/(1.0 + t2);
        Real sy = 2.0*ty/(1.0 + t2);
        Real sz = 2.0*tz/(1.0 + t2);

        Real vpx = vx + (vy*tz - vz*ty);
        Real vpy = vy + (vz*tx - vx*tz);
        Real vpz = vz + (vx*ty - vy*tx);

        vx += vpy*sz - vpz*sy;
        vy += vpz*sx - vpx*sz;
        vz += vpx*sy - vpy*sx;

        // Half electric acceleration
        vx += qdt_2m*cEx;
        vy += qdt_2m*cEy;
        vz += qdt_2m*cEz;

        if (expanding_box_local) {
          Real sx = exp(-0.5*dt*exp_rate_x1_local);
          Real sy = exp(-0.5*dt*exp_rate_x2_local);
          Real sz = exp(-0.5*dt*exp_rate_x3_local);
          vx *= sx;
          vy *= sy;
          if (use_vz_component) {
            vz *= sz;
          }
        }

        // Complete drift from midpoint to full-step position.
        pr(IPX, p) = x_mid + dt_half*vx;
        pr(IPY, p) = y_mid + dt_half*vy;
        pr(IPZ, p) = (nx3_local > 1) ? (z_mid + dt_half*vz) : 0.0;

        // Store updated velocity
        pr(IPVX, p) = vx;
        pr(IPVY, p) = vy;
        pr(IPVZ, p) = use_vz_component ? vz : 0.0;

        // Store sampled midpoint EM fields at particle.
        pr(IPBX, p) = Bx;
        pr(IPBY, p) = By;
        pr(IPBZ, p) = use_vz_component ? Bz : 0.0;
        pr(IPEX, p) = cEx;
        pr(IPEY, p) = cEy;
        pr(IPEZ, p) = use_vz_component ? cEz : 0.0;

        Real v_new_sq = vx*vx + vy*vy + vz*vz;
        pr(IPDPX, p) = m_macro*(vx - vx_old)*inv_dt;
        pr(IPDPY, p) = m_macro*(vy - vy_old)*inv_dt;
        pr(IPDPZ, p) = m_macro*(vz - vz_old)*inv_dt;
        pr(IPDE, p) = 0.5*m_macro*(v_new_sq - v_old_sq)*inv_dt;
        pr(IPEBDOT, p) = cEx*Bx + cEy*By + cEz*Bz;

        // Update displacement tracking if enabled
        if (track_displacement_local) {
          pr(IPDX, p) += dt*vx;
          pr(IPDY, p) += dt*vy;
          if (use_vz_component) {
            pr(IPDZ, p) += dt*vz;
          }
          // Parallel displacement
          Real B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
          if (B_mag > 0.0) {
            pr(IPDB, p) += dt*(vx*Bx + vy*By + vz*Bz)/B_mag;
          }
        }
      });

  return TaskStatus::complete;
}

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3, Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal, Real R200, Real rho_mean) {
  Real R = sqrt(x1 * x1 + x2 * x2);
  Real r2 = x1 * x1 + x2 * x2 + x3 * x3;
  Real r = sqrt(r2);

  // Avoid division by zero
  constexpr Real tiny = 1.0e-20;
  r = fmax(r, tiny);

  // NFW component
  Real x = r / r_s;
  Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log(1 + x) / x;

  // Miyamoto-Nagai model
  Real phi_MN =
      -G * M_gal / sqrt(R * R + SQR(sqrt(x3 * x3 + z_gal * z_gal) + a_gal));

  // Outer component
  Real term1 = (4.0 / 3.0) * pow(5 * R200, 1.5) * sqrt(r);
  Real term2 = (1.0 / 6.0) * r2;
  Real phi_Outer = 4 * M_PI * G * rho_mean * (term1 + term2);

  // Total potential
  Real phi = phi_NFW + phi_MN + phi_Outer;

  return phi;
}

} // namespace particles
