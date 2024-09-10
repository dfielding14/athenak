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
#include "mhd/mhd.hpp"
#include "coordinates/cell_locations.hpp"
namespace particles {
//----------------------------------------------------------------------------------------
//! \fn  void Particles::ParticlesPush
//  \brief

TaskStatus Particles::Push(Driver *pdriver, int stage) {

  auto &npart = nprtcl_thispack;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &b0_ = pmy_pack->pmhd->b0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int ie = indcs.ie;
  const int je = indcs.je;
  const int ke = indcs.ke;
  const int ngh = indcs.ng;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto gids = pmy_pack->gids;
  const bool &multi_d = pmy_pack->pmesh->multi_d;
  const bool &three_d = pmy_pack->pmesh->three_d;
  auto dt_ = (pmy_pack->pmesh->dt);

  switch (pusher) {
    case ParticlesPusher::drift:{

      par_for("part_update",DevExeSpace(),0,(npart-1),
      KOKKOS_LAMBDA(const int p) {
        int m = pi(PGID,p) - gids;
        int ip = (pr(IPX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
        pr(IPX,p) += 0.5*dt_*pr(IPVX,p);

        if (multi_d) {
          int jp = (pr(IPY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
          pr(IPY,p) += 0.5*dt_*pr(IPVY,p);
        }

        if (three_d) {
          int kp = (pr(IPZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;
          pr(IPZ,p) += 0.5*dt_*pr(IPVZ,p);
        }
      });
    break;}
    case ParticlesPusher::boris:{
      auto &b0_ = pmy_pack->pmhd->b0;

      // First half-step in space
      par_for("part_boris",DevExeSpace(),0,(npart-1),
      KOKKOS_LAMBDA(const int p) {

        Real x[3]; //Half-step increment.
        // Get metric components at new location x1,x2,x3
        Real x1p = pr(IPX,p);
        Real x2p = pr(IPY,p);
        Real x3p = pr(IPZ,p);
        
        x[0] = pr(IPX,p) + dt_/(2.0) * pr(IPVX,p);
        if (multi_d) { x[1] = pr(IPY,p) + dt_/(2.0)*pr(IPVY,p);}
        if (three_d) { x[2] = pr(IPZ,p) + dt_/(2.0)*pr(IPVZ,p);}

        Real uB[3]; //Evolution of the velocity due to the magnetic field.

        int m = pi(PGID,p) - gids;
        int ip = (x[0] - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
        int jp = (x[1] - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
        int kp = (x[2] - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;

        Real &x1min = mbsize.d_view(m).x1min;
        Real &x2min = mbsize.d_view(m).x2min;
        Real &x3min = mbsize.d_view(m).x3min;
        Real &x1max = mbsize.d_view(m).x1max;
        Real &x2max = mbsize.d_view(m).x2max;
        Real &x3max = mbsize.d_view(m).x3max;


        Real x1v,x2v,x3v;
        x1v = LeftEdgeX(ip, indcs.nx1, x1min, x1max);
        x2v = LeftEdgeX(jp, indcs.nx2, x2min, x2max);
        x3v = LeftEdgeX(kp, indcs.nx3, x3min, x3max);
        Real Dx,Dy,Dz;
        Dx = (x1max - x1min)/indcs.nx1;
        Dy = (x2max - x2min)/indcs.nx2;
        Dz = (x3max - x3min)/indcs.nx3;

        // Interpolate Magnetic Field at new particle location x1, x2, x3
        // Store it in an array for convenience 
        Real B[3] = {0.0, 0.0, 0.0};
        B[0] = b0_.x1f(m, kp, jp, ip) + (x[0] - x1v)*(b0_.x1f(m, kp, jp, ip+1) - b0_.x1f(m, kp, jp, ip))/Dx;
        B[1] = b0_.x2f(m, kp, jp, ip) + (x[1] - x2v)*(b0_.x2f(m, kp, jp+1, ip) - b0_.x2f(m, kp, jp, ip))/Dy;
        B[2] = b0_.x3f(m, kp, jp, ip) + (x[2] - x3v)*(b0_.x3f(m, kp+1, jp, ip) - b0_.x3f(m, kp, jp, ip))/Dz;

        Real mod_t_sqr = 0.0;
        Real t[3];
        for (int i1 = 0; i1 < 3; ++i1 ){
          t[i1] = B[i1]/pr(IPM,p)  /(2.0)*dt_;
          mod_t_sqr += SQR(t[i1]);
        }

        Real up[3] = {pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p)};
        // Save the vector product of u and t 
        Real vec_ut[3] = {
        up[1]*t[2] - up[2]*t[1],
        up[2]*t[0] - up[0]*t[2],
        up[0]*t[1] - up[1]*t[0],
        };
        uB[0] = up[0] + 2.0/(1.0+mod_t_sqr)*( (up[1] + vec_ut[1])*t[2] - (up[2] + vec_ut[2])*t[1] );

        if (multi_d) { uB[1] = up[1] + 2.0/(1.0+mod_t_sqr)*( (up[2] + vec_ut[2])*t[0] - (up[0] + vec_ut[0])*t[2] ); }
        if (three_d) { uB[2] = up[2] + 2.0/(1.0+mod_t_sqr)*( (up[0] + vec_ut[0])*t[1] - (up[1] + vec_ut[1])*t[0] ); }

        // Finally update velocity in local space
        pr(IPVX,p) = uB[0];
        pr(IPVY,p) = uB[1];
        pr(IPVZ,p) = uB[2];

        pr(IPX,p) = x[0] + dt_/(2.0)*  pr(IPVX,p);
        if (multi_d) { pr(IPY,p) = x[1] + dt_/(2.0) * pr(IPVY,p);}
        if (three_d) { pr(IPZ,p) = x[2] + dt_/(2.0) * pr(IPVZ,p);}

	pr(IPBX,p) = B[0];
	pr(IPBY,p) = B[1];
	pr(IPBZ,p) = B[2];
	
        Real Dx1 = pr(IPX,p)-x1p;
        Real Dx2 = pr(IPY,p)-x2p;
	Real Dx3 = pr(IPZ,p)-x3p;

        pr(IPDX,p) += Dx1;
	pr(IPDY,p) += Dx2;
	pr(IPDZ,p) += Dx3;
	
	Real Bmag = std::sqrt( SQR(B[0]) + SQR(B[1]) + SQR(B[2]) );

	pr(IPDB,p) += (Dx1 * B[0] + Dx2 * B[1] + Dx3 * B[2] ) / Bmag;	

      });
    break;}
  default:
    break;
  }

  return TaskStatus::complete;
}
} // namespace particles

