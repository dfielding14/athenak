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
#include "utils/galaxy_potential.hpp"

namespace particles {

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
      Real c_out = GalaxyPotOuterCoef(r_200);
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
          auto phi = [&](Real qx, Real qy, Real qz) {
            return GalaxyPot(qx, qy, qz, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
          };
          Real phi_xp = phi(px+h, py, pz);
          Real phi_xm = phi(px-h, py, pz);
          Real phi_yp = phi(px, py+h, pz);
          Real phi_ym = phi(px, py-h, pz);
          Real phi_zp = phi(px, py, pz+h);
          Real phi_zm = phi(px, py, pz-h);

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
  default:
    break;
  }

  return TaskStatus::complete;
}

} // namespace particles
