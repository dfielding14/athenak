//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn  void Particles::ParticlesPush
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
  int npart = nprtcl_thispack;

  switch (pusher) {
    case ParticlesPusher::drift:

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

    break;
    case ParticlesPusher::boris_lin:
      {
        auto &b0 = pmy_pack->pmhd->b0;
        par_for("part_boris_lin",DevExeSpace(),0,(npart-1),
        KOKKOS_LAMBDA(const int p) {
          int m = pi(PGID,p) - gids;

          Real x1_old = pr(IPX,p);
          Real x2_old = pr(IPY,p);
          Real x3_old = pr(IPZ,p);

          Real x1 = pr(IPX,p) + 0.5*dt_*pr(IPVX,p);
          Real x2 = pr(IPY,p);
          Real x3 = pr(IPZ,p);
          if (multi_d) {x2 += 0.5*dt_*pr(IPVY,p);}
          if (three_d) {x3 += 0.5*dt_*pr(IPVZ,p);}

          int ip = static_cast<int>((x1 - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1) + is;
          int jp = static_cast<int>((x2 - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2) + js;
          int kp = static_cast<int>((x3 - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3) + ks;

          Real x1v = LeftEdgeX(ip, indcs.nx1, mbsize.d_view(m).x1min,
                               mbsize.d_view(m).x1max);
          Real x2v = LeftEdgeX(jp, indcs.nx2, mbsize.d_view(m).x2min,
                               mbsize.d_view(m).x2max);
          Real x3v = LeftEdgeX(kp, indcs.nx3, mbsize.d_view(m).x3min,
                               mbsize.d_view(m).x3max);

          Real Bx = b0.x1f(m,kp,jp,ip) +
                    (x1 - x1v)*(b0.x1f(m,kp,jp,ip+1) - b0.x1f(m,kp,jp,ip))/
                    mbsize.d_view(m).dx1;
          Real By = b0.x2f(m,kp,jp,ip);
          if (multi_d) {
            By += (x2 - x2v)*(b0.x2f(m,kp,jp+1,ip) - b0.x2f(m,kp,jp,ip))/
                  mbsize.d_view(m).dx2;
          }
          Real Bz = b0.x3f(m,kp,jp,ip);
          if (three_d) {
            Bz += (x3 - x3v)*(b0.x3f(m,kp+1,jp,ip) - b0.x3f(m,kp,jp,ip))/
                  mbsize.d_view(m).dx3;
          }

          pr(IPBX,p) = Bx;
          pr(IPBY,p) = By;
          pr(IPBZ,p) = Bz;

          Real tx = 0.5*dt_*Bx/pr(IPM,p);
          Real ty = 0.5*dt_*By/pr(IPM,p);
          Real tz = 0.5*dt_*Bz/pr(IPM,p);
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

          pr(IPX,p) = x1 + 0.5*dt_*pr(IPVX,p);
          if (multi_d) {pr(IPY,p) = x2 + 0.5*dt_*pr(IPVY,p);}
          if (three_d) {pr(IPZ,p) = x3 + 0.5*dt_*pr(IPVZ,p);}

          Real dx1 = pr(IPX,p) - x1_old;
          Real dx2 = pr(IPY,p) - x2_old;
          Real dx3 = pr(IPZ,p) - x3_old;
          pr(IPDX,p) += dx1;
          pr(IPDY,p) += dx2;
          pr(IPDZ,p) += dx3;

          Real bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
          if (bmag > 0.0) {
            pr(IPDB,p) += (dx1*Bx + dx2*By + dx3*Bz)/bmag;
          }
        });
        break;
      }
    case ParticlesPusher::boris_tsc:
      {
        auto &bcc = pmy_pack->pmhd->bcc0;
        par_for("part_boris_tsc",DevExeSpace(),0,(npart-1),
        KOKKOS_LAMBDA(const int p) {
          int m = pi(PGID,p) - gids;

          Real x1_old = pr(IPX,p);
          Real x2_old = pr(IPY,p);
          Real x3_old = pr(IPZ,p);

          Real isgr = static_cast<Real>(is);
          Real jsgr = static_cast<Real>(js);
          Real ksgr = static_cast<Real>(ks);
          Real wei1[3], wei2[3], wei3[3];

          pr(IPX,p) += 0.5*dt_*pr(IPVX,p);
          Real a = (pr(IPX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + isgr;
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
            pr(IPY,p) += 0.5*dt_*pr(IPVY,p);
            a = (pr(IPY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + jsgr;
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
            pr(IPZ,p) += 0.5*dt_*pr(IPVZ,p);
            a = (pr(IPZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ksgr;
            int kp = static_cast<int>(a);
            k0 = kp - 1;
            d = a - static_cast<Real>(kp);
            wei3[0] = 0.5*(1.0 - d)*(1.0 - d);
            wei3[1] = 0.75 - (d - 0.5)*(d - 0.5);
            wei3[2] = 0.5*d*d;
          }

          Real Bx = 0.0;
          Real By = 0.0;
          Real Bz = 0.0;
          for (int k=0; k<3; ++k) {
            int kk = k0 + k;
            for (int j=0; j<3; ++j) {
              int jj = j0 + j;
              for (int i=0; i<3; ++i) {
                int ii = i0 + i;
                Real w = wei1[i]*wei2[j]*wei3[k];
                Bx += w*bcc(m,IBX,kk,jj,ii);
                By += w*bcc(m,IBY,kk,jj,ii);
                Bz += w*bcc(m,IBZ,kk,jj,ii);
              }
            }
          }

          pr(IPBX,p) = Bx;
          pr(IPBY,p) = By;
          pr(IPBZ,p) = Bz;

          Real bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
          if (pr(IPM,p) < 0.0) {
            if (bmag > 0.0) {
              pr(IPVX,p) = Bx/bmag;
              pr(IPVY,p) = By/bmag;
              pr(IPVZ,p) = Bz/bmag;
            }
          } else {
            Real tx = 0.5*dt_*Bx/pr(IPM,p);
            Real ty = 0.5*dt_*By/pr(IPM,p);
            Real tz = 0.5*dt_*Bz/pr(IPM,p);

            Real v1 = pr(IPVX,p);
            Real v2 = pr(IPVY,p);
            Real v3 = pr(IPVZ,p);
            Real tsq = tx*tx + ty*ty + tz*tz;
            Real sfac = 2.0/(1.0 + tsq);
            Real sx = tx*sfac;
            Real sy = ty*sfac;
            Real sz = tz*sfac;

            Real vp1 = v1 + v2*tz - v3*ty;
            Real vp2 = v2 + v3*tx - v1*tz;
            Real vp3 = v3 + v1*ty - v2*tx;

            pr(IPVX,p) = v1 + vp2*sz - vp3*sy;
            pr(IPVY,p) = v2 + vp3*sx - vp1*sz;
            pr(IPVZ,p) = v3 + vp1*sy - vp2*sx;
          }

          pr(IPX,p) += 0.5*dt_*pr(IPVX,p);
          if (multi_d) {pr(IPY,p) += 0.5*dt_*pr(IPVY,p);}
          if (three_d) {pr(IPZ,p) += 0.5*dt_*pr(IPVZ,p);}

          Real dx1 = pr(IPX,p) - x1_old;
          Real dx2 = pr(IPY,p) - x2_old;
          Real dx3 = pr(IPZ,p) - x3_old;
          pr(IPDX,p) += dx1;
          pr(IPDY,p) += dx2;
          pr(IPDZ,p) += dx3;
          if (bmag > 0.0) {
            pr(IPDB,p) += (dx1*Bx + dx2*By + dx3*Bz)/bmag;
          }
        });
        break;
      }
  default:
    break;
  }

  return TaskStatus::complete;
}
} // namespace particles
