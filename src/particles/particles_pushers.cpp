//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "particles.hpp"
#include "units/units.hpp"
#include "mhd/mhd.hpp"

namespace particles {
KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real R200, Real rho_mean);

//----------------------------------------------------------------------------------------
//! \fn void Particles::ParticlesPush
//  \brief

TaskStatus Particles::Push(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto dt_ = (pmy_pack->pmesh->dt);
  auto gids = pmy_pack->gids;

  switch (pusher) {
    case ParticlesPusher::drift: {
      par_for("part_update",DevExeSpace(),0,(nprtcl_thispack-1),
      KOKKOS_LAMBDA(const int p) {
        int m = pi(PGID,p) - gids;
        int ip = (pr(IPX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
        pr(IPX,p) += dt_*pr(IPVX,p);

        if (multi_d) {
          int jp = (pr(IPY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
          pr(IPY,p) += dt_*pr(IPVY,p);
        }

        if (three_d) {
          int kp = (pr(IPZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;
          pr(IPZ,p) += dt_*pr(IPVZ,p);
        }
      });
    break;}
    case ParticlesPusher::rk4_gravity:{
      Real G = pmy_pack->punit->grav_constant();
      Real r_s = r_scale;
      Real rho_s = rho_scale;
      Real m_g = m_gal;
      Real a_g = a_gal;
      Real z_g = z_gal;
      Real r_m = r_200;
      Real rho_m = rho_mean;
      Real h = par_grav_dx;

      par_for("part_rk4_gravity", DevExeSpace(), 0, (nprtcl_thispack-1),
      KOKKOS_LAMBDA(const int p) {
        // Current state
        Real x  = pr(IPX, p);
        Real y  = pr(IPY, p);
        Real z  = pr(IPZ, p);
        Real vx = pr(IPVX, p);
        Real vy = pr(IPVY, p);
        Real vz = pr(IPVZ, p);

	auto compute_acc = [&](Real px, Real py, Real pz, Real& ax, Real& ay, Real& az) {
          Real phi_xp = GravPot(px+h, py, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_xm = GravPot(px-h, py, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_yp = GravPot(px, py+h, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_ym = GravPot(px, py-h, pz, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_zp = GravPot(px, py, pz+h, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);
          Real phi_zm = GravPot(px, py, pz-h, G, r_s, rho_s, m_g, a_g, z_g, r_m, rho_m);

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
        Real x_mid1  =  x + 0.5*dt_*k1x;
        Real y_mid1  =  y + 0.5*dt_*k1y;
        Real z_mid1  =  z + 0.5*dt_*k1z;
        Real vx_mid1 = vx + 0.5*dt_*k1vx;
        Real vy_mid1 = vy + 0.5*dt_*k1vy;
        Real vz_mid1 = vz + 0.5*dt_*k1vz;
        
        k2x = vx_mid1;
        k2y = vy_mid1;
        k2z = vz_mid1;
        compute_acc(x_mid1, y_mid1, z_mid1, k2vx, k2vy, k2vz);
        
        // K3: derivatives at midpoint using K2
        Real x_mid2  = x  + 0.5*dt_*k2x;
        Real y_mid2  = y  + 0.5*dt_*k2y;
        Real z_mid2  = z  + 0.5*dt_*k2z;
        Real vx_mid2 = vx + 0.5*dt_*k2vx;
        Real vy_mid2 = vy + 0.5*dt_*k2vy;
        Real vz_mid2 = vz + 0.5*dt_*k2vz;
        
        k3x = vx_mid2;
        k3y = vy_mid2;
        k3z = vz_mid2;
        compute_acc(x_mid2, y_mid2, z_mid2, k3vx, k3vy, k3vz);
        
        // K4: derivatives at endpoint using K3
        Real x_end  = x  + dt_*k3x;
        Real y_end  = y  + dt_*k3y;
        Real z_end  = z  + dt_*k3z;
        Real vx_end = vx + dt_*k3vx;
        Real vy_end = vy + dt_*k3vy;
        Real vz_end = vz + dt_*k3vz;
        
        k4x = vx_end;
        k4y = vy_end;
        k4z = vz_end;
        compute_acc(x_end, y_end, z_end, k4vx, k4vy, k4vz);
        
        // Final RK4 update
        pr(IPX, p) = x + dt_/6.0 * (k1x + 2.0*k2x + 2.0*k3x + k4x);
        pr(IPY, p) = y + dt_/6.0 * (k1y + 2.0*k2y + 2.0*k3y + k4y);
        pr(IPZ, p) = z + dt_/6.0 * (k1z + 2.0*k2z + 2.0*k3z + k4z);
        
        pr(IPVX, p) = vx + dt_/6.0 * (k1vx + 2.0*k2vx + 2.0*k3vx + k4vx);
        pr(IPVY, p) = vy + dt_/6.0 * (k1vy + 2.0*k2vy + 2.0*k3vy + k4vy);
        pr(IPVZ, p) = vz + dt_/6.0 * (k1vz + 2.0*k2vz + 2.0*k3vz + k4vz);
      });
    break;}
    case ParticlesPusher::boris_lin:
    case ParticlesPusher::boris_tsc:
      return PushCosmicRays(pdriver, stage);
  default:
    break;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushCosmicRays
//  \brief Boris pusher for cosmic ray particles in electromagnetic fields

TaskStatus Particles::PushCosmicRays(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  Real dt = pmy_pack->pmesh->dt;
  int gids = pmy_pack->gids;
  bool use_tsc = (pusher == ParticlesPusher::boris_tsc);
  
  // Access MHD fields if available
  bool has_mhd = (pmy_pack->pmhd != nullptr);
  DvceArray5D<Real> bcc;
  if (has_mhd) {
    bcc = pmy_pack->pmhd->bcc0;
  }
  
  // Boris pusher algorithm
  par_for("push_cr", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    // Get particle properties
    int m = pi(PGID,p) - gids;
    Real x = pr(IPX,p);
    Real y = pr(IPY,p);
    Real z = (indcs.nx3 > 1) ? pr(IPZ,p) : 0.0;
    Real vx = pr(IPVX,p);
    Real vy = pr(IPVY,p);
    Real vz = (indcs.nx3 > 1) ? pr(IPVZ,p) : 0.0;
    Real q_over_m = pr(IPM,p);
    
    // Interpolate B-field to particle position
    Real Bx = 0.0, By = 0.0, Bz = 0.0;
    if (has_mhd) {
      // Simple linear interpolation for now
      // TODO: Add TSC interpolation option
      Real dx1 = size.d_view(m).dx1;
      Real dx2 = size.d_view(m).dx2;
      Real dx3 = (indcs.nx3 > 1) ? size.d_view(m).dx3 : 1.0;
      
      int i = static_cast<int>((x - size.d_view(m).x1min) / dx1) + indcs.is;
      int j = static_cast<int>((y - size.d_view(m).x2min) / dx2) + indcs.js;
      int k = (indcs.nx3 > 1) ? 
              static_cast<int>((z - size.d_view(m).x3min) / dx3) + indcs.ks : indcs.ks;
      
      // Ensure indices are within bounds
      i = (i < indcs.is) ? indcs.is : ((i > indcs.ie) ? indcs.ie : i);
      j = (j < indcs.js) ? indcs.js : ((j > indcs.je) ? indcs.je : j);
      k = (k < indcs.ks) ? indcs.ks : ((k > indcs.ke) ? indcs.ke : k);
      
      // Simple nearest-neighbor for now
      Bx = bcc(m,IBX,k,j,i);
      By = bcc(m,IBY,k,j,i);
      if (indcs.nx3 > 1) {
        Bz = bcc(m,IBZ,k,j,i);
      }
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
    Real t2 = tx*tx + ty*ty + tz*tz;
    Real s = 2.0 / (1.0 + t2);
    
    Real vx1 = vx + vy*tz - vz*ty;
    Real vy1 = vy + vz*tx - vx*tz;
    Real vz1 = vz + vx*ty - vy*tx;
    
    vx += s*(vy1*tz - vz1*ty);
    vy += s*(vz1*tx - vx1*tz);
    vz += s*(vx1*ty - vy1*tx);
    
    // Half electric acceleration
    vx += qdt_2m * Ex;
    vy += qdt_2m * Ey;
    vz += qdt_2m * Ez;
    
    // Update position
    pr(IPX,p) += dt * vx;
    pr(IPY,p) += dt * vy;
    if (indcs.nx3 > 1) {
      pr(IPZ,p) += dt * vz;
    }
    
    // Store updated velocity
    pr(IPVX,p) = vx;
    pr(IPVY,p) = vy;
    if (indcs.nx3 > 1) {
      pr(IPVZ,p) = vz;
    }
    
    // Store B-field at particle
    pr(IPBX,p) = Bx;
    pr(IPBY,p) = By;
    if (indcs.nx3 > 1) {
      pr(IPBZ,p) = Bz;
    }
    
    // Update displacement tracking if enabled
    if (track_displacement) {
      pr(IPDX,p) += dt * vx;
      pr(IPDY,p) += dt * vy;
      if (indcs.nx3 > 1) {
        pr(IPDZ,p) += dt * vz;
      }
      // Parallel displacement
      Real B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
      if (B_mag > 0.0) {
        pr(IPDB,p) += dt * (vx*Bx + vy*By + vz*Bz) / B_mag;
      }
    }
  });
  
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushStars
//  \brief Separate method for star particle dynamics

TaskStatus Particles::PushStars(Driver *pdriver, int stage) {
  // For now, just call the RK4 gravity pusher code
  // This is a placeholder for future star-specific dynamics
  return Push(pdriver, stage);
}

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real R200, Real rho_mean) {
  Real R = sqrt(x1*x1 + x2*x2);
  Real r2 = x1*x1 + x2*x2 + x3*x3;
  Real r = sqrt(r2);

  // Avoid division by zero
  constexpr Real tiny = 1.0e-20;
  r = fmax(r, tiny);

  // NFW component
  Real x = r / r_s;
  Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log(1 + x) / x;

  // Miyamoto-Nagai model
  Real phi_MN = -G * M_gal / sqrt(R*R + SQR(sqrt(x3*x3 + z_gal*z_gal) + a_gal));

  // Outer component
  Real term1 = (4.0 / 3.0) * pow(5 * R200, 1.5) * sqrt(r);
  Real term2 = (1.0 / 6.0) * r2;
  Real phi_Outer = 4 * M_PI * G * rho_mean * (term1 + term2);

  // Total potential
  Real phi = phi_NFW + phi_MN + phi_Outer;

  return phi;
}

} // namespace particles
