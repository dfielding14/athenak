//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"
#include "units/units.hpp"
#include "globals.hpp"

namespace particles {
KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3, Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal, Real R200, Real rho_mean);

KOKKOS_INLINE_FUNCTION
Real HashToUnitReal(int64_t seed) {
  uint64_t z = static_cast<uint64_t>(seed);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z = z ^ (z >> 31);
  return static_cast<Real>(z & 0x7FFFFFFFULL) / static_cast<Real>(0x80000000ULL);
}

KOKKOS_INLINE_FUNCTION
bool IsPhysicalInflowParticleBoundary(BoundaryFlag bc) {
  return (bc != BoundaryFlag::block && bc != BoundaryFlag::periodic &&
          bc != BoundaryFlag::shear_periodic && bc != BoundaryFlag::undef);
}

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
  case ParticlesPusher::lagrangian_mc:
    return PushLagrangianMC(pdriver, stage);
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
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  Real dt = pmy_pack->pmesh->dt;
  int gids = pmy_pack->gids;
  bool use_tsc = (pusher == ParticlesPusher::boris_tsc);

  // Access MHD fields if available
  bool has_mhd = (pmy_pack->pmhd != nullptr);
  if (!has_mhd) {
    std::cout << "### FATAL ERROR in Particles::PushCosmicRays: "
              << "Boris pushers require MHD fields" << std::endl;
    return TaskStatus::fail;
  }

  Particles *pt = this;
  
  // Create local copies to avoid GPU lambda capture warning
  bool track_displacement_local = track_displacement;
  const int nx3_local = indcs.nx3;  // avoid capturing host ref

  // Boris pusher algorithm
  par_for(
      "push_cr", DevExeSpace(), 0, nprtcl_thispack - 1,
      KOKKOS_LAMBDA(const int p) {
        // Get particle properties
        int m = pi(PGID, p) - gids;
        Real x = pr(IPX, p);
        Real y = pr(IPY, p);
        Real z = (indcs.nx3 > 1) ? pr(IPZ, p) : 0.0;
        Real vx = pr(IPVX, p);
        Real vy = pr(IPVY, p);
        Real vz = (indcs.nx3 > 1) ? pr(IPVZ, p) : 0.0;
        Real q_over_m = pr(IPM, p);

        // Interpolate B-field to particle position
        Real Bx = 0.0, By = 0.0, Bz = 0.0;
        if (use_tsc) {
          pt->InterpolateTSC(m, x, y, z, Bx, By, Bz);
        } else {
          pt->InterpolateLinear(m, x, y, z, Bx, By, Bz);
        }

        // Boris algorithm
        Real dt_half = 0.5 * dt;
        Real qdt_2m = q_over_m * dt_half;

        // No electric field for now
        Real Ex = 0.0, Ey = 0.0, Ez = 0.0;

        // Half electric acceleration
        vx += qdt_2m * Ex;
        vy += qdt_2m * Ey;
        vz += qdt_2m * Ez;

        // Magnetic rotation
        Real tx = qdt_2m * Bx;
        Real ty = qdt_2m * By;
        Real tz = qdt_2m * Bz;
        Real t2 = tx * tx + ty * ty + tz * tz;
        Real s = 2.0 / (1.0 + t2);

        Real vx1 = vx + vy * tz - vz * ty;
        Real vy1 = vy + vz * tx - vx * tz;
        Real vz1 = vz + vx * ty - vy * tx;

        vx += s * (vy1 * tz - vz1 * ty);
        vy += s * (vz1 * tx - vx1 * tz);
        vz += s * (vx1 * ty - vy1 * tx);

        // Half electric acceleration
        vx += qdt_2m * Ex;
        vy += qdt_2m * Ey;
        vz += qdt_2m * Ez;

        // Update position
        pr(IPX, p) += dt * vx;
        pr(IPY, p) += dt * vy;
        if (indcs.nx3 > 1) {
          pr(IPZ, p) += dt * vz;
        }

        // Store updated velocity
        pr(IPVX, p) = vx;
        pr(IPVY, p) = vy;
        if (indcs.nx3 > 1) {
          pr(IPVZ, p) = vz;
        }

        // Store B-field at particle
        pr(IPBX, p) = Bx;
        pr(IPBY, p) = By;
        if (nx3_local > 1) {
          pr(IPBZ, p) = Bz;
        }

        // Update displacement tracking if enabled
        if (track_displacement_local) {
          pr(IPDX, p) += dt * vx;
          pr(IPDY, p) += dt * vy;
          if (nx3_local > 1) {
            pr(IPDZ, p) += dt * vz;
          }
          // Parallel displacement
          Real B_mag = sqrt(Bx * Bx + By * By + Bz * Bz);
          if (B_mag > 0.0) {
            pr(IPDB, p) += dt * (vx * Bx + vy * By + vz * Bz) / B_mag;
          }
        }
      });

  return TaskStatus::complete;
}

KOKKOS_INLINE_FUNCTION
void Particles::InterpolateLinear(int m, Real x, Real y, Real z, Real &Bx, Real &By,
                                  Real &Bz) const {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto bcc = pmy_pack->pmhd->bcc0;

  Real dx1 = size.d_view(m).dx1;
  Real dx2 = size.d_view(m).dx2;
  Real dx3 = (indcs.nx3 > 1) ? size.d_view(m).dx3 : 1.0;

  Real fx = (x - size.d_view(m).x1min) / dx1;
  Real fy = (y - size.d_view(m).x2min) / dx2;
  Real fz = (indcs.nx3 > 1) ? (z - size.d_view(m).x3min) / dx3 : 0.0;

  int i0 = static_cast<int>(floor(fx)) + indcs.is;
  int j0 = static_cast<int>(floor(fy)) + indcs.js;
  int k0 = (indcs.nx3 > 1) ? static_cast<int>(floor(fz)) + indcs.ks : indcs.ks;

  i0 = (i0 < indcs.is) ? indcs.is : ((i0 > indcs.ie - 1) ? indcs.ie - 1 : i0);
  j0 = (j0 < indcs.js) ? indcs.js : ((j0 > indcs.je - 1) ? indcs.je - 1 : j0);
  k0 = (k0 < indcs.ks) ? indcs.ks : ((k0 > indcs.ke - 1) ? indcs.ke - 1 : k0);

  int i1 = i0 + 1;
  int j1 = j0 + 1;
  int k1 = (indcs.nx3 > 1) ? k0 + 1 : k0;

  Real wx = fx - floor(fx);
  Real wy = fy - floor(fy);
  Real wz = (indcs.nx3 > 1) ? fz - floor(fz) : 0.0;

  Bx = By = Bz = 0.0;

  for (int dk = 0; dk <= (indcs.nx3 > 1 ? 1 : 0); ++dk) {
    Real wk = (indcs.nx3 > 1) ? (dk == 0 ? (1.0 - wz) : wz) : 1.0;
    int kk = (dk == 0) ? k0 : k1;
    for (int dj = 0; dj <= 1; ++dj) {
      Real wj = (dj == 0) ? (1.0 - wy) : wy;
      int jj = (dj == 0) ? j0 : j1;
      for (int di = 0; di <= 1; ++di) {
        Real wi = (di == 0) ? (1.0 - wx) : wx;
        int ii = (di == 0) ? i0 : i1;
        Real w = wi * wj * wk;
        Bx += w * bcc(m, IBX, kk, jj, ii);
        By += w * bcc(m, IBY, kk, jj, ii);
        if (indcs.nx3 > 1) {
          Bz += w * bcc(m, IBZ, kk, jj, ii);
        }
      }
    }
  }
  if (indcs.nx3 == 1) {
    Bz = 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
void Particles::InterpolateTSC(int m, Real x, Real y, Real z, Real &Bx, Real &By,
                               Real &Bz) const {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto bcc = pmy_pack->pmhd->bcc0;

  Real dx1 = size.d_view(m).dx1;
  Real dx2 = size.d_view(m).dx2;
  Real dx3 = (indcs.nx3 > 1) ? size.d_view(m).dx3 : 1.0;

  Real fx = (x - size.d_view(m).x1min) / dx1 - 0.5;
  Real fy = (y - size.d_view(m).x2min) / dx2 - 0.5;
  Real fz = (indcs.nx3 > 1) ? (z - size.d_view(m).x3min) / dx3 - 0.5 : 0.0;

  int ic = static_cast<int>(floor(fx)) + indcs.is;
  int jc = static_cast<int>(floor(fy)) + indcs.js;
  int kc = (indcs.nx3 > 1) ? static_cast<int>(floor(fz)) + indcs.ks : indcs.ks;

  Real di = fx - floor(fx);
  Real dj = fy - floor(fy);
  Real dk = (indcs.nx3 > 1) ? fz - floor(fz) : 0.0;

  auto weight = [](Real d) {
    Real ad = fabs(d);
    if (ad < 0.5) return 0.75 - ad * ad;
    if (ad < 1.5) {
      Real t = 1.5 - ad;
      return 0.5 * t * t;
    }
    return 0.0;
  };

  Real wx[3] = {weight(di + 1.0), weight(di), weight(di - 1.0)};
  Real wy[3] = {weight(dj + 1.0), weight(dj), weight(dj - 1.0)};
  Real wz[3] = {0.0, 1.0, 0.0};
  int kz[3] = {kc, kc, kc};
  if (indcs.nx3 > 1) {
    wz[0] = weight(dk + 1.0);
    wz[1] = weight(dk);
    wz[2] = weight(dk - 1.0);
    kz[0] = kc - 1;
    kz[1] = kc;
    kz[2] = kc + 1;
  }

  int ix[3] = {ic - 1, ic, ic + 1};
  int iy[3] = {jc - 1, jc, jc + 1};
  for (int n = 0; n < 3; ++n) {
    ix[n] = (ix[n] < indcs.is) ? indcs.is : ((ix[n] > indcs.ie) ? indcs.ie : ix[n]);
    iy[n] = (iy[n] < indcs.js) ? indcs.js : ((iy[n] > indcs.je) ? indcs.je : iy[n]);
    if (indcs.nx3 > 1) {
      kz[n] = (kz[n] < indcs.ks) ? indcs.ks : ((kz[n] > indcs.ke) ? indcs.ke : kz[n]);
    }
  }

  Bx = By = Bz = 0.0;
  int nk = (indcs.nx3 > 1) ? 3 : 1;
  for (int kk = 0; kk < nk; ++kk) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        Real w = wx[ii] * wy[jj] * wz[kk];
        Bx += w * bcc(m, IBX, kz[kk], iy[jj], ix[ii]);
        By += w * bcc(m, IBY, kz[kk], iy[jj], ix[ii]);
        if (indcs.nx3 > 1) {
          Bz += w * bcc(m, IBZ, kz[kk], iy[jj], ix[ii]);
        }
      }
    }
  }
  if (indcs.nx3 == 1) {
    Bz = 0.0;
  }
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

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::PushLagrangianMC
//! \brief push with Lagrangian Monte Carlo method (Genel+ 2013, MNRAS.435.1426G)
//!        AMR corrections are applied after particle communication and mesh refinement.

TaskStatus Particles::PushLagrangianMC(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &gids = pmy_pack->gids;
  auto &mblev = pmy_pack->pmb->mb_lev;

  auto &u1_ = (pmy_pack->phydro != nullptr)?pmy_pack->phydro->u1:pmy_pack->pmhd->u1;
  auto &uflxidn_ = (pmy_pack->phydro != nullptr)?
                    pmy_pack->phydro->uflxidnsaved:pmy_pack->pmhd->uflxidnsaved;
  auto &flx1_ = uflxidn_.x1f;
  auto &flx2_ = uflxidn_.x2f;
  auto &flx3_ = uflxidn_.x3f;

  int ncycle = pmy_pack->pmesh->ncycle;
  int64_t rseed = random_seed;  // Capture random_seed to avoid implicit 'this' capture

  int nmb = pmy_pack->nmb_thispack;

  par_for("part_update",DevExeSpace(),0,(nprtcl_thispack-1),
  KOKKOS_LAMBDA(const int p) {
    if (pi(PLASTMOVE,p) >= 0) {
      // only update particles that are not frozen or marked for deletion

      int m = pi(PGID,p) - gids;

      // Bounds check: ensure m is valid for this rank
      if (m < 0 || m >= nmb) {
        pi(PLASTMOVE,p) = -1;  // Mark as frozen
        return;
      }

      int ip = (pr(LMCX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
      int jp = js;
      int kp = ks;

      if (multi_d) {
        jp = (pr(LMCY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
      }

      if (three_d) {
        kp = (pr(LMCZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;
      }

      // Minimal bounds check to prevent out-of-bounds flux array access
      int ie = is + indcs.nx1;
      int je = js + indcs.nx2;
      int ke = ks + indcs.nx3;

      if (ip < is || ip >= ie || jp < js || jp >= je || kp < ks || kp >= ke) {
        // Particle is outside active zone - skip this timestep
        return;
      }

      // get normalized fluxes based on local density
      Real mass = u1_(m,IDN,kp,jp,ip);

      // by convention, these values will be positive when there is outflow
      // with respect to the current particle's cell
      Real flx1_left = -flx1_(m,kp,jp,ip) / mass;
      Real flx1_right = flx1_(m,kp,jp,ip+1) / mass;
      Real flx2_left = (multi_d) ? -flx2_(m,kp,jp,ip) / mass : 0.;
      Real flx2_right = (multi_d) ? flx2_(m,kp,jp+1,ip) / mass : 0.;
      Real flx3_left = (three_d) ? -flx3_(m,kp,jp,ip) / mass : 0.;
      Real flx3_right = (three_d) ? flx3_(m,kp+1,jp,ip) / mass : 0.;

      flx1_left = flx1_left < 0 ? 0 : flx1_left;
      flx1_right = flx1_right < 0 ? 0 : flx1_right;
      flx2_left = flx2_left < 0 ? 0 : flx2_left;
      flx2_right = flx2_right < 0 ? 0 : flx2_right;
      flx3_left = flx3_left < 0 ? 0 : flx3_left;
      flx3_right = flx3_right < 0 ? 0 : flx3_right;

      // Deterministic random seed: tag * prime1 + ncycle * prime2 + input_seed
      // Using large primes to avoid correlations: 7919 and 104729
      int64_t det_seed = pi(PTAG,p) * 7919 + ncycle * 104729 + rseed;

      // Hash-based pseudo-random number generation (splitmix64 algorithm)
      // Fast, stateless, device-compatible alternative to Kokkos random pool
      uint64_t z = static_cast<uint64_t>(det_seed);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      z = z ^ (z >> 31);
      Real rand = static_cast<Real>(z & 0x7FFFFFFFULL) / static_cast<Real>(0x80000000ULL);

      // save refinement level of current zone
      pi(PLASTLEVEL,p) = mblev.d_view(m);

      // save parity of current zone stored as (i_isodd,j_isodd,k_isodd) * 8
      pi(PLASTMOVE,p) = 32 * (ip % 2) + 16 * (jp % 2) + 8 * (kp % 2);

      if (rand < flx1_left) {
        pr(LMCX,p) -= mbsize.d_view(m).dx1;
        pi(PLASTMOVE,p) += 1;
      } else if (rand < flx1_left + flx1_right) {
        pr(LMCX,p) += mbsize.d_view(m).dx1;
        pi(PLASTMOVE,p) += 2;
      } else if (multi_d && rand < flx1_left + flx1_right + flx2_left) {
        pr(LMCY,p) -= mbsize.d_view(m).dx2;
        pi(PLASTMOVE,p) += 3;
      } else if (multi_d && rand < flx1_left + flx1_right + flx2_left + flx2_right) {
        pr(LMCY,p) += mbsize.d_view(m).dx2;
        pi(PLASTMOVE,p) += 4;
      } else if (three_d && rand < flx1_left + flx1_right + flx2_left + flx2_right
                                + flx3_left) {
        pr(LMCZ,p) -= mbsize.d_view(m).dx3;
        pi(PLASTMOVE,p) += 5;
      } else if (three_d && rand < flx1_left + flx1_right + flx2_left + flx2_right
                                + flx3_left + flx3_right) {
        pr(LMCZ,p) += mbsize.d_view(m).dx3;
        pi(PLASTMOVE,p) += 6;
      }
    }
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::InjectLagrangianMCInflow
//! \brief add MC tracer particles through physical domain faces with inward mass flux

TaskStatus Particles::InjectLagrangianMCInflow(Driver *pdriver, int stage) {
  if (!lmc_inject_at_inflow || pusher != ParticlesPusher::lagrangian_mc ||
      lmc_mass_per_particle <= 0.0) {
    return TaskStatus::complete;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  bool multi_d = pmy_pack->pmesh->multi_d;
  bool three_d = pmy_pack->pmesh->three_d;
  int nmb = pmy_pack->nmb_thispack;

  int x1_face_size = nx2*nx3;
  int x2_face_size = nx1*nx3;
  int x3_face_size = nx1*nx2;
  int slots_per_mb = 2*x1_face_size;
  if (multi_d) { slots_per_mb += 2*x2_face_size; }
  if (three_d) { slots_per_mb += 2*x3_face_size; }
  int nslots = nmb*slots_per_mb;
  if (nslots <= 0) {
    return TaskStatus::complete;
  }

  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto &uflxidn_ = (pmy_pack->phydro != nullptr)?
                   pmy_pack->phydro->uflxidnsaved:pmy_pack->pmhd->uflxidnsaved;
  auto &flx1_ = uflxidn_.x1f;
  auto &flx2_ = uflxidn_.x2f;
  auto &flx3_ = uflxidn_.x3f;
  Real mass_per_particle = lmc_mass_per_particle;
  int64_t seed = lmc_inject_seed;
  int ncycle = pmy_pack->pmesh->ncycle;
  int gids = pmy_pack->gids;

  if (lmc_scheduled_tracking) {
    int max_injected_track_thisrank = -1;
    if (nprtcl_thispack > 0) {
      auto pi_for_track_scan = prtcl_idata;
      const int ntrack_initial = lmc_track_ninitial;
      Kokkos::parallel_reduce("lmc_track_slot_scan",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
          KOKKOS_LAMBDA(const int p, int &max_slot) {
            const int slot = pi_for_track_scan(PTRACK,p);
            if (slot >= ntrack_initial && slot > max_slot) {
              max_slot = slot;
            }
          }, Kokkos::Max<int>(max_injected_track_thisrank));
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &max_injected_track_thisrank, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
#endif
    if (max_injected_track_thisrank >= lmc_track_ninitial) {
      lmc_track_next_injected_slot =
          std::max(lmc_track_next_injected_slot,
                   max_injected_track_thisrank - lmc_track_ninitial + 1);
    }
  }

  DualArray2D<int> inject_count("lmc_inflow_count", nslots, 4);
  auto &count = inject_count;
  par_for("lmc_inflow_count", DevExeSpace(), 0, nslots-1,
  KOKKOS_LAMBDA(const int slot) {
    int m = slot / slots_per_mb;
    int rem = slot - m*slots_per_mb;
    int face = BoundaryFace::inner_x1;
    int i = is;
    int j = js;
    int k = ks;
    Real signed_inflow = 0.0;

    if (rem < x1_face_size) {
      int q = rem;
      k = ks + q/nx2;
      j = js + (q - (q/nx2)*nx2);
      i = is;
      face = BoundaryFace::inner_x1;
      signed_inflow = flx1_(m,k,j,is);
    } else if (rem < 2*x1_face_size) {
      int q = rem - x1_face_size;
      k = ks + q/nx2;
      j = js + (q - (q/nx2)*nx2);
      i = ie;
      face = BoundaryFace::outer_x1;
      signed_inflow = -flx1_(m,k,j,ie+1);
    } else if (multi_d && rem < 2*x1_face_size + x2_face_size) {
      int q = rem - 2*x1_face_size;
      k = ks + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      j = js;
      face = BoundaryFace::inner_x2;
      signed_inflow = flx2_(m,k,js,i);
    } else if (multi_d && rem < 2*x1_face_size + 2*x2_face_size) {
      int q = rem - 2*x1_face_size - x2_face_size;
      k = ks + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      j = je;
      face = BoundaryFace::outer_x2;
      signed_inflow = -flx2_(m,k,je+1,i);
    } else if (three_d && rem < 2*x1_face_size + 2*x2_face_size + x3_face_size) {
      int q = rem - 2*x1_face_size - 2*x2_face_size;
      j = js + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      k = ks;
      face = BoundaryFace::inner_x3;
      signed_inflow = flx3_(m,ks,j,i);
    } else {
      int q = rem - 2*x1_face_size - 2*x2_face_size - x3_face_size;
      j = js + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      k = ke;
      face = BoundaryFace::outer_x3;
      signed_inflow = -flx3_(m,ke+1,j,i);
    }

    int np = 0;
    if (IsPhysicalInflowParticleBoundary(mb_bcs.d_view(m,face))) {
      Real vol = mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*mbsize.d_view(m).dx3;
      Real expected = fmax(signed_inflow, 0.0)*vol/mass_per_particle;
      np = static_cast<int>(expected);
      Real frac = expected - static_cast<Real>(np);
      int64_t local_seed = seed + static_cast<int64_t>(gids + m)*999983LL +
                           static_cast<int64_t>(face)*13007LL +
                           static_cast<int64_t>(i)*7919LL +
                           static_cast<int64_t>(j)*104729LL +
                           static_cast<int64_t>(k)*524287LL +
                           static_cast<int64_t>(ncycle)*15485863LL;
      if (HashToUnitReal(local_seed) < frac) {
        np += 1;
      }
    }
    count.d_view(slot,0) = np;
    count.d_view(slot,1) = 0;
    count.d_view(slot,2) = 0;
    count.d_view(slot,3) = 0;
  });

  inject_count.template modify<DevExeSpace>();
  inject_count.template sync<HostMemSpace>();

  auto slot_face = [&](const int slot) {
    const int rem = slot - (slot/slots_per_mb)*slots_per_mb;
    if (rem < x1_face_size) {
      return BoundaryFace::inner_x1;
    } else if (rem < 2*x1_face_size) {
      return BoundaryFace::outer_x1;
    } else if (multi_d && rem < 2*x1_face_size + x2_face_size) {
      return BoundaryFace::inner_x2;
    } else if (multi_d && rem < 2*x1_face_size + 2*x2_face_size) {
      return BoundaryFace::outer_x2;
    } else if (three_d && rem < 2*x1_face_size + 2*x2_face_size + x3_face_size) {
      return BoundaryFace::inner_x3;
    }
    return BoundaryFace::outer_x3;
  };

  int nnew = 0;
  int remaining = lmc_max_inject_per_step;
  for (int slot=0; slot<nslots; ++slot) {
    int np = inject_count.h_view(slot,0);
    if (lmc_max_inject_per_step >= 0) {
      if (remaining <= 0) {
        np = 0;
      } else if (np > remaining) {
        np = remaining;
      }
      remaining -= np;
    }
    inject_count.h_view(slot,0) = np;
    inject_count.h_view(slot,1) = nnew;
    nnew += np;
  }

  int ntrack_new_total = 0;
  int ntrack_new_thisrank = 0;
  int first_track_slot =
      lmc_track_ninitial + lmc_track_next_injected_slot;
  if (lmc_scheduled_tracking) {
    const int ntrack_injected =
        lmc_track_nparticles_total - lmc_track_ninitial;
    if (ntrack_injected > 0 &&
        lmc_track_next_injected_slot < ntrack_injected) {
      const Real time = pmy_pack->pmesh->time;
      int due_total = 0;
      if (time >= lmc_track_inject_t_start) {
        if (lmc_track_inject_t_stop <= lmc_track_inject_t_start) {
          due_total = ntrack_injected;
        } else {
          const Real frac =
              (time - lmc_track_inject_t_start)/
              (lmc_track_inject_t_stop - lmc_track_inject_t_start);
          due_total = static_cast<int>(
              std::floor(frac*static_cast<Real>(ntrack_injected) + 0.5));
        }
      }
      due_total = std::max(0, std::min(ntrack_injected, due_total));
      const int due_remaining =
          std::max(0, due_total - lmc_track_next_injected_slot);
      if (due_remaining > 0) {
        int eligible_thisrank = 0;
        for (int slot=0; slot<nslots; ++slot) {
          if (slot_face(slot) == lmc_track_face) {
            eligible_thisrank += inject_count.h_view(slot,0);
          }
        }
        std::vector<int> eligible_eachrank(global_variable::nranks, 0);
        eligible_eachrank[global_variable::my_rank] = eligible_thisrank;
#if MPI_PARALLEL_ENABLED
        MPI_Allgather(&eligible_thisrank, 1, MPI_INT, eligible_eachrank.data(), 1,
                      MPI_INT, MPI_COMM_WORLD);
#endif
        int eligible_before = 0;
        int eligible_total = 0;
        for (int n=0; n<global_variable::nranks; ++n) {
          if (n < global_variable::my_rank) {
            eligible_before += eligible_eachrank[n];
          }
          eligible_total += eligible_eachrank[n];
        }
        ntrack_new_total = std::min(due_remaining, eligible_total);
        const int local_first_global = std::max(eligible_before, 0);
        const int local_last_global =
            std::min(eligible_before + eligible_thisrank, ntrack_new_total);
        ntrack_new_thisrank = std::max(0, local_last_global - local_first_global);
        const int skip_thisrank =
            std::max(0, local_first_global - eligible_before);
        first_track_slot =
            lmc_track_ninitial + lmc_track_next_injected_slot + local_first_global;

        int eligible_seen = 0;
        int track_prefix = 0;
        for (int slot=0; slot<nslots; ++slot) {
          int track_count = 0;
          if (slot_face(slot) == lmc_track_face) {
            const int np = inject_count.h_view(slot,0);
            const int begin = eligible_seen;
            const int end = eligible_seen + np;
            const int assign_begin = std::max(begin, skip_thisrank);
            const int assign_end =
                std::min(end, skip_thisrank + ntrack_new_thisrank);
            track_count = std::max(0, assign_end - assign_begin);
            eligible_seen = end;
          }
          inject_count.h_view(slot,2) = track_count;
          inject_count.h_view(slot,3) = track_prefix;
          track_prefix += track_count;
        }
      }
    }
  }

  inject_count.template modify<HostMemSpace>();
  inject_count.template sync<DevExeSpace>();

  if (lmc_initial_next_tag <= 0) {
    int max_tag_thisrank = -1;
    if (nprtcl_thispack > 0) {
      auto old_pi_for_tags = prtcl_idata;
      Kokkos::parallel_reduce("lmc_initial_tag_scan",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
          KOKKOS_LAMBDA(const int p, int &max_tag) {
            if (old_pi_for_tags(PTAG,p) > max_tag) {
              max_tag = old_pi_for_tags(PTAG,p);
            }
          }, Kokkos::Max<int>(max_tag_thisrank));
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &max_tag_thisrank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    lmc_initial_next_tag = max_tag_thisrank + 1;
    next_tag = std::max(next_tag, lmc_initial_next_tag);
  }

  std::vector<int> nnew_eachrank(global_variable::nranks, 0);
  nnew_eachrank[global_variable::my_rank] = nnew;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&nnew, 1, MPI_INT, nnew_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#endif
  int tagstart = next_tag;
  int nnew_total = 0;
  for (int n=0; n<global_variable::nranks; ++n) {
    if (n < global_variable::my_rank) {
      tagstart += nnew_eachrank[n];
    }
    nnew_total += nnew_eachrank[n];
  }
  next_tag += nnew_total;
  lmc_track_next_injected_slot += ntrack_new_total;

  if (nnew_total <= 0) {
    return TaskStatus::complete;
  }
  if (nnew <= 0) {
    pmy_pack->pmesh->CountParticles();
    return TaskStatus::complete;
  }

  int old_npart = nprtcl_thispack;
  int new_npart = old_npart + nnew;
  DvceArray2D<Real> new_pr("lmc_inflow_rdata", nrdata, new_npart);
  DvceArray2D<int> new_pi("lmc_inflow_idata", nidata, new_npart);
  auto old_pr = prtcl_rdata;
  auto old_pi = prtcl_idata;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto counts = inject_count;

  if (old_npart > 0) {
    par_for("lmc_inflow_copy", DevExeSpace(), 0, old_npart-1,
    KOKKOS_LAMBDA(const int p) {
      for (int n=0; n<nidata; ++n) {
        new_pi(n,p) = old_pi(n,p);
      }
      for (int n=0; n<nrdata; ++n) {
        new_pr(n,p) = old_pr(n,p);
      }
    });
  }

  par_for("lmc_inflow_fill", DevExeSpace(), 0, nslots-1,
  KOKKOS_LAMBDA(const int slot) {
    int np = counts.d_view(slot,0);
    if (np <= 0) {
      return;
    }

    int m = slot / slots_per_mb;
    int rem = slot - m*slots_per_mb;
    int i = is;
    int j = js;
    int k = ks;

    if (rem < x1_face_size) {
      int q = rem;
      k = ks + q/nx2;
      j = js + (q - (q/nx2)*nx2);
      i = is;
    } else if (rem < 2*x1_face_size) {
      int q = rem - x1_face_size;
      k = ks + q/nx2;
      j = js + (q - (q/nx2)*nx2);
      i = ie;
    } else if (multi_d && rem < 2*x1_face_size + x2_face_size) {
      int q = rem - 2*x1_face_size;
      k = ks + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      j = js;
    } else if (multi_d && rem < 2*x1_face_size + 2*x2_face_size) {
      int q = rem - 2*x1_face_size - x2_face_size;
      k = ks + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      j = je;
    } else if (three_d && rem < 2*x1_face_size + 2*x2_face_size + x3_face_size) {
      int q = rem - 2*x1_face_size - 2*x2_face_size;
      j = js + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      k = ks;
    } else {
      int q = rem - 2*x1_face_size - 2*x2_face_size - x3_face_size;
      j = js + q/nx1;
      i = is + (q - (q/nx1)*nx1);
      k = ke;
    }

    Real x = mbsize.d_view(m).x1min +
             (static_cast<Real>(i - is) + 0.5)*mbsize.d_view(m).dx1;
    Real y = mbsize.d_view(m).x2min +
             (static_cast<Real>(j - js) + 0.5)*mbsize.d_view(m).dx2;
    Real z = mbsize.d_view(m).x3min +
             (static_cast<Real>(k - ks) + 0.5)*mbsize.d_view(m).dx3;
    int first = old_npart + counts.d_view(slot,1);
    for (int n=0; n<np; ++n) {
      int p = first + n;
      new_pi(PGID,p) = gids + m;
      new_pi(PTAG,p) = tagstart + counts.d_view(slot,1) + n;
      new_pi(PLASTMOVE,p) = 0;
      new_pi(PLASTLEVEL,p) = mblev.d_view(m);
      new_pi(PTRACK,p) = -1;
      if (n < counts.d_view(slot,2)) {
        new_pi(PTRACK,p) = first_track_slot + counts.d_view(slot,3) + n;
      }
      new_pr(LMCX,p) = x;
      new_pr(LMCY,p) = y;
      new_pr(LMCZ,p) = z;
    }
  });
  Kokkos::fence();

  prtcl_rdata = new_pr;
  prtcl_idata = new_pi;
  nprtcl_thispack = new_npart;
  pmy_pack->pmesh->CountParticles();

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::AdjustMeshRefinement
//! \brief update locations of particles that enter meshblocks with new refinement levels

TaskStatus Particles::AdjustMeshRefinement(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &gids = pmy_pack->gids;
  auto &mblev = pmy_pack->pmb->mb_lev;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;

  auto &uflxidn_ = (pmy_pack->phydro != nullptr)?
                   pmy_pack->phydro->uflxidnsaved:pmy_pack->pmhd->uflxidnsaved;
  auto &flx1_ = uflxidn_.x1f;
  auto &flx2_ = uflxidn_.x2f;
  auto &flx3_ = uflxidn_.x3f;

  int ncycle = pmy_pack->pmesh->ncycle;
  int64_t rseed = random_seed;  // Capture random_seed to avoid implicit 'this' capture
  int nmb = pmy_pack->nmb_thispack;

  par_for("particle_meshshift",DevExeSpace(),0,(nprtcl_thispack-1),
  KOKKOS_LAMBDA(const int p) {
    if (pi(PLASTMOVE,p) >= 0) {
      // only update particles that are not frozen or marked for deletion

      int m = pi(PGID,p) - gids;

      // Bounds check: ensure m is valid for this rank
      if (m < 0 || m >= nmb) {
        pi(PLASTMOVE,p) = -1;  // Mark as frozen
        return;
      }

      int level = mblev.d_view(m);

      int lastlevel = pi(PLASTLEVEL,p);
      int lastmove = pi(PLASTMOVE,p);

      // oddness of the last cell that the particle lived in
      int i_parity = lastmove / 32;
      int j_parity = (lastmove % 32) / 16;
      int k_parity = (lastmove % 16) / 8;

      // direction of last move:
      //   1 -> "left" x1 face was chosen
      //   2 -> "right" x1 face was chosen
      //   3 -> "left" x2 face was chosen
      //   4 -> "right" x2 face was chosen
      //   5 -> "left" x3 face was chosen
      //   6 -> "right" x3 face was chosen
      lastmove = lastmove % 8;

      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = multi_d ? mbsize.d_view(m).dx2 : 0.;
      Real dx3 = three_d ? mbsize.d_view(m).dx3 : 0.;

      if (level > lastlevel) {
        // this is a higher refinement level, i.e., the zones are smaller now

        if (lastmove == 1) {
          // came from zone to right (dx--)
          pr(LMCX,p) += dx1/2;

          pr(LMCY,p) -= dx2/2;
          pr(LMCZ,p) -= dx3/2;
        } else if (lastmove == 2) {
          // came from zone to left (dx++)
          pr(LMCX,p) -= dx1/2;

          pr(LMCY,p) -= dx2/2;
          pr(LMCZ,p) -= dx3/2;
        } else if (multi_d && lastmove == 3) {
          // came from zone above (dy--)
          pr(LMCY,p) += dx2/2;

          pr(LMCX,p) -= dx1/2;
          pr(LMCZ,p) -= dx3/2;
        } else if (multi_d && lastmove == 4) {
          // came from zone below (dy++)
          pr(LMCY,p) -= dx2/2;

          pr(LMCX,p) -= dx1/2;
          pr(LMCZ,p) -= dx3/2;
        } else if (three_d && lastmove == 5) {
          // came from zone in front (dz--)
          pr(LMCZ,p) += dx3/2;

          pr(LMCX,p) -= dx1/2;
          pr(LMCY,p) -= dx2/2;
        } else if (three_d && lastmove == 6) {
          // came from zone behind (dz++)
          pr(LMCZ,p) -= dx3/2;

          pr(LMCX,p) -= dx1/2;
          pr(LMCY,p) -= dx2/2;
        }

        int ip = (pr(LMCX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
        int jp = js;
        int kp = ks;

        if (multi_d) {
          jp = (pr(LMCY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
        }

        if (three_d) {
          kp = (pr(LMCZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;
        }

        // get fluxes into the four zones that the particle could have ended up in
        Real flx1 = 0.;
        Real flx2 = 0.;
        Real flx3 = 0.;
        Real flx4 = 0.;

        if (lastmove == 1) {
          // came from zone to the right
          flx1 = -flx1_(m,kp,jp,ip+1);
          flx2 = (multi_d) ? -flx1_(m,kp,jp+1,ip+1) : 0.;
          flx3 = (three_d) ? -flx1_(m,kp+1,jp,ip+1) : 0.;
          flx4 = (multi_d && three_d) ? -flx1_(m,kp+1,jp+1,ip+1) : 0.;
        } else if (lastmove == 2) {
          // came from zone to the left
          flx1 = flx1_(m,kp,jp,ip);
          flx2 = (multi_d) ? flx1_(m,kp,jp+1,ip) : 0.;
          flx3 = (three_d) ? flx1_(m,kp+1,jp,ip) : 0.;
          flx4 = (multi_d && three_d) ? flx1_(m,kp+1,jp+1,ip) : 0.;
        } else if (lastmove == 3) {
          // came from zone above. is at least multi_d
          flx1 = -flx2_(m,kp,jp+1,ip);
          flx2 = -flx2_(m,kp,jp+1,ip+1);
          flx3 = (three_d) ? -flx2_(m,kp+1,jp+1,ip) : 0.;
          flx4 = (three_d) ? -flx2_(m,kp+1,jp+1,ip+1) : 0.;
        } else if (lastmove == 4) {
          // came from zone below. is at least multi_d
          flx1 = flx2_(m,kp,jp,ip);
          flx2 = flx2_(m,kp,jp,ip+1);
          flx3 = (three_d) ? flx2_(m,kp+1,jp,ip) : 0.;
          flx4 = (three_d) ? flx2_(m,kp+1,jp,ip+1) : 0.;
        } else if (lastmove == 5) {
          // came from zone in front. is three_d
          flx1 = -flx3_(m,kp+1,jp,ip);
          flx2 = -flx3_(m,kp+1,jp,ip+1);
          flx3 = -flx3_(m,kp+1,jp+1,ip);
          flx4 = -flx3_(m,kp+1,jp+1,ip+1);
        } else if (lastmove == 6) {
          // came from zone behind. is three_d
          flx1 = flx3_(m,kp,jp,ip);
          flx2 = flx3_(m,kp,jp,ip+1);
          flx3 = flx3_(m,kp,jp+1,ip);
          flx4 = flx3_(m,kp,jp+1,ip+1);
        }

        flx1 = (flx1 < 0) ? 0. : flx1;
        flx2 = (flx2 < 0) ? 0. : flx2;
        flx3 = (flx3 < 0) ? 0. : flx3;
        flx4 = (flx4 < 0) ? 0. : flx4;

        Real flx_total = flx1 + flx2 + flx3 + flx4;
        flx_total = (flx_total > 0) ? flx_total : 1.e-10;

        flx1 /= flx_total;
        flx2 /= flx_total;
        flx3 /= flx_total;
        flx4 /= flx_total;

        // Deterministic random seed: tag * prime1 + ncycle * prime2 + input_seed
        // Using large primes to avoid correlations: 7919 and 104729
        // Add 1 to differentiate from the seed used in PushLagrangianMC
        int64_t det_seed = pi(PTAG,p) * 7919 + ncycle * 104729 + rseed + 1;

        // Hash-based pseudo-random number generation (splitmix64 algorithm)
        // Fast, stateless, device-compatible alternative to Kokkos random pool
        uint64_t z = static_cast<uint64_t>(det_seed);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        Real rand = static_cast<Real>(z & 0x7FFFFFFFULL) / static_cast<Real>(0x80000000ULL);

        int target_zone = 4;
        if (rand < flx1) {
          target_zone = 1;
        } else if (rand < flx1 + flx2) {
          target_zone = 2;
        } else if (rand < flx1 + flx2 + flx3) {
          target_zone = 3;
        }

        if (lastmove == 1 || lastmove == 2) {
          if (target_zone == 2) {
            pr(LMCY,p) += mbsize.d_view(m).dx2;
          } else if (target_zone == 3) {
            pr(LMCZ,p) += mbsize.d_view(m).dx3;
          } else if (target_zone == 4) {
            pr(LMCY,p) += mbsize.d_view(m).dx2;
            pr(LMCZ,p) += mbsize.d_view(m).dx3;
          }
        } else if (lastmove == 3 || lastmove == 4) {
          if (target_zone == 2) {
            pr(LMCX,p) += mbsize.d_view(m).dx1;
          } else if (target_zone == 3) {
            pr(LMCZ,p) += mbsize.d_view(m).dx3;
          } else if (target_zone == 4) {
            pr(LMCX,p) += mbsize.d_view(m).dx1;
            pr(LMCZ,p) += mbsize.d_view(m).dx3;
          }
        } else if (lastmove == 5 || lastmove == 6) {
          if (target_zone == 2) {
            pr(LMCX,p) += mbsize.d_view(m).dx1;
          } else if (target_zone == 3) {
            pr(LMCY,p) += mbsize.d_view(m).dx2;
          } else if (target_zone == 4) {
            pr(LMCX,p) += mbsize.d_view(m).dx1;
            pr(LMCY,p) += mbsize.d_view(m).dx2;
          }
        }

      } else if (level < lastlevel) {
        // this is a lower refinement level, i.e., the zones are larger now,
        // there's nothing special to do other than to move the particle to
        // the center of the new zone

        if (i_parity) {
          pr(LMCX,p) -= mbsize.d_view(m).dx1/4;
        } else {
          pr(LMCX,p) += mbsize.d_view(m).dx1/4;
        }
        if (multi_d) {
          if (j_parity) {
            pr(LMCY,p) -= mbsize.d_view(m).dx2/4;
          } else {
            pr(LMCY,p) += mbsize.d_view(m).dx2/4;
          }
        }
        if (three_d) {
          if (k_parity) {
            pr(LMCZ,p) -= mbsize.d_view(m).dx3/4;
          } else {
            pr(LMCZ,p) += mbsize.d_view(m).dx3/4;
          }
        }

        if (lastmove == 1) {
          // came from zone to right (dx--)
          pr(LMCX,p) -= mbsize.d_view(m).dx1/2;
        } else if (lastmove == 2) {
          // came from zone to left (dx++)
          pr(LMCX,p) += mbsize.d_view(m).dx1/2;
        } else if (lastmove == 3) {
          // came from zone above (dy--)
          pr(LMCY,p) -= mbsize.d_view(m).dx2/2;
        } else if (lastmove == 4) {
          // came from zone below (dy++)
          pr(LMCY,p) += mbsize.d_view(m).dx2/2;
        } else if (lastmove == 5) {
          // came from zone in front (dz--)
          pr(LMCZ,p) -= mbsize.d_view(m).dx3/2;
        } else if (lastmove == 6) {
          // came from zone behind (dz++)
          pr(LMCZ,p) += mbsize.d_view(m).dx3/2;
        }
      }
    }
  });

  return TaskStatus::complete;
}

} // namespace particles
