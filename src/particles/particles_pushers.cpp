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
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "particles/field_interpolation.hpp"
#include "particles.hpp"
#include "units/units.hpp"

namespace particles {
namespace {

KOKKOS_INLINE_FUNCTION
Real ParticleBoundaryRoundoff() {
#if SINGLE_PRECISION_ENABLED
  return static_cast<Real>(1.1920928955078125e-7);
#else
  return static_cast<Real>(2.2204460492503131e-16);
#endif
}

KOKKOS_INLINE_FUNCTION
Real InteriorUpperBoundary(const Real lower, const Real upper, const Real dx) {
  const Real scale = fmax(static_cast<Real>(1.0),
                          fmax(fabs(lower), fmax(fabs(upper), fabs(dx))));
  const Real tol = static_cast<Real>(64.0)*ParticleBoundaryRoundoff()*scale;
  const Real interior = upper - tol;
  return (interior > lower) ? interior : static_cast<Real>(0.5)*(lower + upper);
}

KOKKOS_INLINE_FUNCTION
void ReflectCoordinateAtPhysicalBoundary(const BoundaryFlag inner_flag,
                                         const BoundaryFlag outer_flag,
                                         const Real lower, const Real upper,
                                         const Real dx, Real &x, Real &v) {
  if (inner_flag == BoundaryFlag::reflect && x < lower) {
    x = static_cast<Real>(2.0)*lower - x;
    v = -v;
    if (x < lower) {
      x = lower;
    }
  }
  if (outer_flag == BoundaryFlag::reflect && x >= upper) {
    x = static_cast<Real>(2.0)*upper - x;
    v = -v;
    if (x >= upper) {
      x = InteriorUpperBoundary(lower, upper, dx);
    }
  }
}

template <typename SizeView, typename BoundaryView>
KOKKOS_INLINE_FUNCTION
void ApplyReflectiveParticleBCs(const int m, const SizeView size,
                                const BoundaryView mb_bcs,
                                const bool multi_d, const bool three_d,
                                Real &x, Real &y, Real &z,
                                Real &vx, Real &vy, Real &vz) {
  const auto block_size = size.d_view(m);
  ReflectCoordinateAtPhysicalBoundary(
      mb_bcs.d_view(m, BoundaryFace::inner_x1),
      mb_bcs.d_view(m, BoundaryFace::outer_x1),
      block_size.x1min, block_size.x1max, block_size.dx1, x, vx);
  if (multi_d) {
    ReflectCoordinateAtPhysicalBoundary(
        mb_bcs.d_view(m, BoundaryFace::inner_x2),
        mb_bcs.d_view(m, BoundaryFace::outer_x2),
        block_size.x2min, block_size.x2max, block_size.dx2, y, vy);
  }
  if (three_d) {
    ReflectCoordinateAtPhysicalBoundary(
        mb_bcs.d_view(m, BoundaryFace::inner_x3),
        mb_bcs.d_view(m, BoundaryFace::outer_x3),
        block_size.x3min, block_size.x3max, block_size.dx3, z, vz);
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
  if (UsesPaperVL2Coupling() && stage == 0) {
    return TaskStatus::complete;
  }
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus status = TaskStatus::fail;
  switch (pusher) {
  case ParticlesPusher::drift:
    status = PushDrift(pdriver, stage);
    break;
  case ParticlesPusher::rk4_gravity:
    status = PushStars(pdriver, stage);
    break;
  case ParticlesPusher::boris_lin:
  case ParticlesPusher::boris_tsc:
    status = PushCosmicRays(pdriver, stage);
    break;
  default:
    break;
  }
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::push, q017_timer.seconds());
  return status;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushDrift
//! \brief Simple drift pusher updating positions using velocities

TaskStatus Particles::PushDrift(Driver *pdriver, int stage) {
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int gids = pmy_pack->gids;
  const int nmb = pmy_pack->nmb_thispack;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto dt_ = (pmy_pack->pmesh->dt);
  size.template sync<DevExeSpace>();
  mb_bcs.template sync<DevExeSpace>();
  auto size_view = size;
  auto mb_bcs_view = mb_bcs;

  par_for(
      "part_update", DevExeSpace(), 0, (nprtcl_thispack - 1),
      KOKKOS_LAMBDA(const int p) {
        const int m = pi(PGID, p) - gids;
        pr(IPX, p) += dt_ * pr(IPVX, p);

        if (multi_d) {
          pr(IPY, p) += dt_ * pr(IPVY, p);
        }

        if (three_d) {
          pr(IPZ, p) += dt_ * pr(IPVZ, p);
        }
        if (m >= 0 && m < nmb) {
          Real x = pr(IPX, p);
          Real y = pr(IPY, p);
          Real z = pr(IPZ, p);
          Real vx = pr(IPVX, p);
          Real vy = pr(IPVY, p);
          Real vz = pr(IPVZ, p);
          ApplyReflectiveParticleBCs(m, size_view, mb_bcs_view, multi_d, three_d,
                                     x, y, z, vx, vy, vz);
          pr(IPX, p) = x;
          pr(IPY, p) = y;
          pr(IPZ, p) = z;
          pr(IPVX, p) = vx;
          pr(IPVY, p) = vy;
          pr(IPVZ, p) = vz;
        }
      });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushStars
//! \brief RK4 gravity pusher for star particles

TaskStatus Particles::PushStars(Driver *pdriver, int stage) {
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int gids = pmy_pack->gids;
  const int nmb = pmy_pack->nmb_thispack;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto dt_ = (pmy_pack->pmesh->dt);
  size.template sync<DevExeSpace>();
  mb_bcs.template sync<DevExeSpace>();
  auto size_view = size;
  auto mb_bcs_view = mb_bcs;

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
        const int m = pi(PGID, p) - gids;
        if (m >= 0 && m < nmb) {
          Real xr = pr(IPX, p);
          Real yr = pr(IPY, p);
          Real zr = pr(IPZ, p);
          Real vxr = pr(IPVX, p);
          Real vyr = pr(IPVY, p);
          Real vzr = pr(IPVZ, p);
          ApplyReflectiveParticleBCs(m, size_view, mb_bcs_view, multi_d, three_d,
                                     xr, yr, zr, vxr, vyr, vzr);
          pr(IPX, p) = xr;
          pr(IPY, p) = yr;
          pr(IPZ, p) = zr;
          pr(IPVX, p) = vxr;
          pr(IPVY, p) = vyr;
          pr(IPVZ, p) = vzr;
        }
      });

  return TaskStatus::complete;
}

template <typename SizeView>
KOKKOS_INLINE_FUNCTION
void InterpolateLinearFields(const RegionIndcs indcs, const SizeView size,
                             const DvceArray5D<Real> bcc,
                             const DvceArray5D<Real> w0,
                             const bool use_mhd_fluid_velocity, const int m,
                             Real x, Real y, Real z, Real &Bx, Real &By,
                             Real &Bz, Real &Ux, Real &Uy, Real &Uz,
                             bool allow_2d3v) {
  const bool three_d = (indcs.nx3 > 1);
  const bool use_bz_channel = three_d || allow_2d3v;
  Real dx1 = size.d_view(m).dx1;
  Real dx2 = size.d_view(m).dx2;
  Real dx3 = three_d ? size.d_view(m).dx3 : 1.0;

  Real fx = (x - size.d_view(m).x1min) / dx1;
  Real fy = (y - size.d_view(m).x2min) / dx2;
  Real fz = three_d ? (z - size.d_view(m).x3min) / dx3 : 0.0;

  int i0 = static_cast<int>(floor(fx)) + indcs.is;
  int j0 = static_cast<int>(floor(fy)) + indcs.js;
  int k0 = three_d ? static_cast<int>(floor(fz)) + indcs.ks : indcs.ks;

  i0 = (i0 < indcs.is) ? indcs.is : ((i0 > indcs.ie - 1) ? indcs.ie - 1 : i0);
  j0 = (j0 < indcs.js) ? indcs.js : ((j0 > indcs.je - 1) ? indcs.je - 1 : j0);
  k0 = (k0 < indcs.ks) ? indcs.ks : ((k0 > indcs.ke - 1) ? indcs.ke - 1 : k0);

  int i1 = i0 + 1;
  int j1 = j0 + 1;
  int k1 = three_d ? k0 + 1 : k0;

  Real wx = fx - floor(fx);
  Real wy = fy - floor(fy);
  Real wz = three_d ? fz - floor(fz) : 0.0;

  Bx = By = Bz = 0.0;
  Ux = Uy = Uz = 0.0;

  for (int dk = 0; dk <= (three_d ? 1 : 0); ++dk) {
    Real wk = three_d ? (dk == 0 ? (1.0 - wz) : wz) : 1.0;
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
        if (use_mhd_fluid_velocity) {
          Ux += w*w0(m, IVX, kk, jj, ii);
          Uy += w*w0(m, IVY, kk, jj, ii);
          Uz += w*w0(m, IVZ, kk, jj, ii);
        }
        if (use_bz_channel) {
          Bz += w * bcc(m, IBZ, kk, jj, ii);
        }
      }
    }
  }
  if (!use_bz_channel) {
    Bz = 0.0;
    Uz = 0.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushPaperCosmicRaysVL2
//! \brief Sun & Bai VL2 split pusher for paper MHD-PIC coupling.

TaskStatus Particles::PushPaperCosmicRaysVL2(Driver *pdriver, int stage) {
  (void)pdriver;
  if (stage != 1 && stage != 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Paper MHD-PIC VL2 particle push requires stage 1 or 2"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // Stage 1 deposits predictor rho/J at x_ini without changing particle
  // momenta. Stage 2 applies the full-step Boris kick at x_mid before the
  // impulse deposit. Position updates are a separate post-deposit task.
  if (stage == 1) return TaskStatus::complete;
  if (pmy_pack->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Paper MHD-PIC VL2 particle push requires active MHD fields"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const RegionIndcs indcs = pmy_pack->pmesh->mb_indcs;
  const Real dt = pmy_pack->pmesh->dt;
  const Real dt_half = static_cast<Real>(0.5)*dt;
  const Real inv_dt = (dt > 0.0) ? (static_cast<Real>(1.0)/dt) :
                                   static_cast<Real>(0.0);
  const int gids = pmy_pack->gids;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx3 = indcs.nx3;
  const bool allow_2d3v = (pic_enable_2d3v && (nx3 == 1));
  const bool use_vz_component = (nx3 > 1) || allow_2d3v;
  const Real qscale = deposit_qscale;
  const int nspecies_local = nspecies;
  const Real light_speed = pic_cr_light_speed;
  const bool deltaf_local = UsesDeltaF();
  const bool adaptive_deltaf_local = UsesAdaptiveDeltaF();
  const PICDeltaFBackground deltaf_background_local = pic_deltaf_background;
  const Real deltaf_p0_local = pic_deltaf_p0;
  const Real deltaf_kappa_local = pic_deltaf_kappa;
  const Real deltaf_adaptive_xi_local = pic_deltaf_adaptive_xi;
  const Real deltaf_adaptive_p0_local = pic_deltaf_adaptive_p0;
  const Real deltaf_drift_x1_local = pic_deltaf_drift_x1;
  const Real deltaf_drift_x2_local = pic_deltaf_drift_x2;
  const Real deltaf_drift_x3_local = pic_deltaf_drift_x3;
  const Real deltaf_aniso_x1_local = pic_deltaf_aniso_x1;
  const Real deltaf_aniso_x2_local = pic_deltaf_aniso_x2;
  const Real deltaf_aniso_x3_local = pic_deltaf_aniso_x3;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &size = pmy_pack->pmb->mb_size;
  auto bcc = pmy_pack->pmhd->bcc0;
  auto w0 = pmy_pack->pmhd->w0;
  auto mspecies = species_mass;
  size.template sync<DevExeSpace>();
  auto size_view = size;

  par_for(
      "push_cr_paper_vl2_midpoint_kick",
      DevExeSpace(), 0, nprtcl_thispack - 1,
      KOKKOS_LAMBDA(const int p) {
        const int m = pi(PGID, p) - gids;
        if (m < 0 || m >= nmb) return;
        Real x = pr(IPX, p);
        Real y = pr(IPY, p);
        Real z = (nx3 > 1) ? pr(IPZ, p) : static_cast<Real>(0.0);
        Real state_x = pr(IPVX, p);
        Real state_y = pr(IPVY, p);
        Real state_z = use_vz_component ? pr(IPVZ, p) : static_cast<Real>(0.0);
        Real vx, vy, vz;
        CRVelocityFromState(true, light_speed, state_x, state_y, state_z,
                            vx, vy, vz);

        Real Bx = 0.0, By = 0.0, Bz = 0.0;
        Real Ux = 0.0, Uy = 0.0, Uz = 0.0;
        InterpolateTSCFields(indcs, size_view, bcc, w0, true, m, x, y, z,
                             Bx, By, Bz, Ux, Uy, Uz, allow_2d3v);
        if (!use_vz_component) {
          Bz = 0.0;
          Uz = 0.0;
        }
        const Real cEx = -(Uy*Bz - Uz*By);
        const Real cEy = -(Uz*Bx - Ux*Bz);
        const Real cEz = use_vz_component ? -(Ux*By - Uy*Bx) :
                                           static_cast<Real>(0.0);
        const Real state_x_before = state_x;
        const Real state_y_before = state_y;
        const Real state_z_before = state_z;
        const Real qdt_2m = pr(IPM, p)*dt_half;

        state_x += qdt_2m*cEx;
        state_y += qdt_2m*cEy;
        state_z += qdt_2m*cEz;
        const Real inv_gamma_minus =
            static_cast<Real>(1.0)/
            CRLorentzFactor(state_x, state_y, state_z, light_speed);
        const Real tx = qdt_2m*Bx*inv_gamma_minus;
        const Real ty = qdt_2m*By*inv_gamma_minus;
        const Real tz = qdt_2m*Bz*inv_gamma_minus;
        const Real t2 = tx*tx + ty*ty + tz*tz;
        const Real rot_x = static_cast<Real>(2.0)*tx/(static_cast<Real>(1.0) + t2);
        const Real rot_y = static_cast<Real>(2.0)*ty/(static_cast<Real>(1.0) + t2);
        const Real rot_z = static_cast<Real>(2.0)*tz/(static_cast<Real>(1.0) + t2);
        const Real state_px = state_x + (state_y*tz - state_z*ty);
        const Real state_py = state_y + (state_z*tx - state_x*tz);
        const Real state_pz = state_z + (state_x*ty - state_y*tx);
        state_x += state_py*rot_z - state_pz*rot_y;
        state_y += state_pz*rot_x - state_px*rot_z;
        state_z += state_px*rot_y - state_py*rot_x;
        state_x += qdt_2m*cEx;
        state_y += qdt_2m*cEy;
        state_z += qdt_2m*cEz;

        const Real energy_before =
            CRKineticEnergy(true, light_speed, state_x_before, state_y_before,
                            state_z_before);
        const Real energy_after =
            CRKineticEnergy(true, light_speed, state_x, state_y, state_z);
        const int sp = pi(PSP, p);
        if (sp < 0 || sp >= nspecies_local) return;
        Real weight = pr(IPWT, p);
        if (weight <= static_cast<Real>(0.0)) weight = static_cast<Real>(1.0);
        const Real macro_mass = qscale*weight*mspecies(sp);
        pr(IPDPX, p) = macro_mass*(state_x - state_x_before)*inv_dt;
        pr(IPDPY, p) = macro_mass*(state_y - state_y_before)*inv_dt;
        pr(IPDPZ, p) = macro_mass*(state_z - state_z_before)*inv_dt;
        pr(IPDE, p) = macro_mass*(energy_after - energy_before)*inv_dt;
        pr(IPEBDOT, p) = cEx*Bx + cEy*By + cEz*Bz;
        pr(IPBX, p) = Bx;
        pr(IPBY, p) = By;
        pr(IPBZ, p) = use_vz_component ? Bz : static_cast<Real>(0.0);
        pr(IPEX, p) = cEx;
        pr(IPEY, p) = cEy;
        pr(IPEZ, p) = use_vz_component ? cEz : static_cast<Real>(0.0);

        pr(IPVX, p) = state_x;
        pr(IPVY, p) = state_y;
        pr(IPVZ, p) = use_vz_component ? state_z : static_cast<Real>(0.0);
        if (deltaf_local) {
          const Real f0_current = adaptive_deltaf_local ?
              PICAdaptiveDeltaFBackgroundValue(
                  deltaf_kappa_local, deltaf_p0_local, deltaf_adaptive_p0_local,
                  deltaf_adaptive_xi_local, static_cast<Real>(1.0),
                  static_cast<Real>(1.0), static_cast<Real>(1.0),
                  state_x, state_y, state_z) :
              PICDeltaFBackgroundValue(
                  deltaf_background_local, deltaf_p0_local, deltaf_kappa_local,
                  deltaf_drift_x1_local, deltaf_drift_x2_local,
                  deltaf_drift_x3_local, deltaf_aniso_x1_local,
                  deltaf_aniso_x2_local, deltaf_aniso_x3_local,
                  static_cast<Real>(1.0), static_cast<Real>(1.0),
                  static_cast<Real>(1.0), state_x, state_y, state_z);
          pr(IPDFWT, p) = static_cast<Real>(1.0) - f0_current/pr(IPF0, p);
        }
      });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::DriftPaperCosmicRaysHalfStep
//! \brief Drift paper-mode CRs after depositing predictor or impulse moments.

TaskStatus Particles::DriftPaperCosmicRaysHalfStep(Driver *pdriver, int stage) {
  (void)pdriver;
  if (stage != 1 && stage != 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Paper MHD-PIC VL2 half drift requires stage 1 or 2"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const RegionIndcs indcs = pmy_pack->pmesh->mb_indcs;
  const Real dt_half = static_cast<Real>(0.5)*pmy_pack->pmesh->dt;
  const int gids = pmy_pack->gids;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx3 = indcs.nx3;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool allow_2d3v = (pic_enable_2d3v && (nx3 == 1));
  const bool use_vz_component = (nx3 > 1) || allow_2d3v;
  const bool track_displacement_local = track_displacement;
  const Real light_speed = pic_cr_light_speed;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto bcc = pmy_pack->pmhd->bcc0;
  auto w0 = pmy_pack->pmhd->w0;
  size.template sync<DevExeSpace>();
  mb_bcs.template sync<DevExeSpace>();
  auto size_view = size;
  auto mb_bcs_view = mb_bcs;

  par_for("drift_cr_paper_vl2_half_step", DevExeSpace(), 0, nprtcl_thispack - 1,
      KOKKOS_LAMBDA(const int p) {
        const int m = pi(PGID, p) - gids;
        if (m < 0 || m >= nmb) return;
        Real x = pr(IPX, p);
        Real y = pr(IPY, p);
        Real z = (nx3 > 1) ? pr(IPZ, p) : static_cast<Real>(0.0);
        Real state_x = pr(IPVX, p);
        Real state_y = pr(IPVY, p);
        Real state_z = use_vz_component ? pr(IPVZ, p) : static_cast<Real>(0.0);
        Real vx, vy, vz;
        CRVelocityFromState(true, light_speed, state_x, state_y, state_z,
                            vx, vy, vz);
        if (track_displacement_local) {
          Real Bx = pr(IPBX, p);
          Real By = pr(IPBY, p);
          Real Bz = use_vz_component ? pr(IPBZ, p) : static_cast<Real>(0.0);
          if (stage == 1) {
            Real Ux = 0.0;
            Real Uy = 0.0;
            Real Uz = 0.0;
            InterpolateTSCFields(indcs, size_view, bcc, w0, true, m, x, y, z,
                                 Bx, By, Bz, Ux, Uy, Uz, allow_2d3v);
          }
          const Real bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
          if (bmag > static_cast<Real>(0.0)) {
            pr(IPDB, p) += dt_half*(vx*Bx + vy*By + vz*Bz)/bmag;
          }
        }
        x += dt_half*vx;
        y += dt_half*vy;
        if (nx3 > 1) z += dt_half*vz;
        ApplyReflectiveParticleBCs(m, size_view, mb_bcs_view, multi_d, three_d,
                                   x, y, z, state_x, state_y, state_z);
        pr(IPX, p) = x;
        pr(IPY, p) = y;
        pr(IPZ, p) = z;
        pr(IPVX, p) = state_x;
        pr(IPVY, p) = state_y;
        pr(IPVZ, p) = use_vz_component ? state_z : static_cast<Real>(0.0);
        if (track_displacement_local) {
          pr(IPDX, p) += dt_half*vx;
          pr(IPDY, p) += dt_half*vy;
          if (use_vz_component) pr(IPDZ, p) += dt_half*vz;
        }
      });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::PushCosmicRays
//  \brief Boris pusher for cosmic ray particles in electromagnetic fields

TaskStatus Particles::PushCosmicRays(Driver *pdriver, int stage) {
  if (UsesPaperVL2Coupling()) {
    return PushPaperCosmicRaysVL2(pdriver, stage);
  }
  const RegionIndcs indcs = pmy_pack->pmesh->mb_indcs;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto &size = pmy_pack->pmb->mb_size;
  Real dt = pmy_pack->pmesh->dt;
  int gids = pmy_pack->gids;
  bool use_tsc = (pusher == ParticlesPusher::boris_tsc);

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

  // Create local copies to avoid GPU lambda capture warning
  bool track_displacement_local = track_displacement;
  const bool multi_d_local = pmy_pack->pmesh->multi_d;
  const int nx3_local = indcs.nx3;  // avoid capturing host ref
  const bool allow_2d3v_local = (pic_enable_2d3v && (nx3_local == 1));
  const bool use_vz_component = (nx3_local > 1) || allow_2d3v_local;
  const bool use_mhd_fluid_velocity =
      (!no_mhd_mode && (pmy_pack->pmhd != nullptr));
  auto bcc = use_mhd_fluid_velocity ? pmy_pack->pmhd->bcc0 : pic_no_mhd_bcc0;
  auto w0 = use_mhd_fluid_velocity ? pmy_pack->pmhd->w0 : DvceArray5D<Real>();
  size.template sync<DevExeSpace>();
  auto size_view = size;

  const Real qscale = deposit_qscale;
  auto mspecies = species_mass;
  const int nspecies_local = nspecies;
  const int nmb_local = pmy_pack->nmb_thispack;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  mb_bcs.template sync<DevExeSpace>();
  auto mb_bcs_view = mb_bcs;
  const Real inv_dt = (dt > 0.0) ? (1.0/dt) : 0.0;
  const bool expanding_box_local =
      (pic_expanding_box_mode == PICExpandingBoxMode::on);
  const PICExpansionLaw expansion_law_local = pic_expansion_law;
  const Real exp_rate_x1_local = pic_expansion_rate_x1;
  const Real exp_rate_x2_local = pic_expansion_rate_x2;
  const Real exp_rate_x3_local = pic_expansion_rate_x3;
  const bool momentum_state_local = UsesRelativisticCRState();
  const Real light_speed_local = pic_cr_light_speed;
  const bool deltaf_local = UsesDeltaF();
  const bool adaptive_deltaf_local = UsesAdaptiveDeltaF();
  const PICDeltaFBackground deltaf_background_local = pic_deltaf_background;
  const Real deltaf_p0_local = pic_deltaf_p0;
  const Real deltaf_kappa_local = pic_deltaf_kappa;
  const Real deltaf_adaptive_xi_local = pic_deltaf_adaptive_xi;
  const Real deltaf_adaptive_p0_local = pic_deltaf_adaptive_p0;
  const Real deltaf_drift_x1_local = pic_deltaf_drift_x1;
  const Real deltaf_drift_x2_local = pic_deltaf_drift_x2;
  const Real deltaf_drift_x3_local = pic_deltaf_drift_x3;
  const Real deltaf_aniso_x1_local = pic_deltaf_aniso_x1;
  const Real deltaf_aniso_x2_local = pic_deltaf_aniso_x2;
  const Real deltaf_aniso_x3_local = pic_deltaf_aniso_x3;
  const Real time_start = pmy_pack->pmesh->time;
  const Real time_mid = time_start + 0.5*dt;
  const Real time_end = time_start + dt;
  const Real a1_start = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x1_local, time_start) : 1.0;
  const Real a2_start = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x2_local, time_start) : 1.0;
  const Real a3_start = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x3_local, time_start) : 1.0;
  const Real a1_mid = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x1_local, time_mid) : 1.0;
  const Real a2_mid = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x2_local, time_mid) : 1.0;
  const Real a3_mid = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x3_local, time_mid) : 1.0;
  const Real a1_end = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x1_local, time_end) : 1.0;
  const Real a2_end = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x2_local, time_end) : 1.0;
  const Real a3_end = expanding_box_local ? PICScaleFactor(
      expansion_law_local, exp_rate_x3_local, time_end) : 1.0;

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
        Real state_x = pr(IPVX, p);
        Real state_y = pr(IPVY, p);
        Real state_z = use_vz_component ? pr(IPVZ, p) : 0.0;
        Real vx, vy, vz;
        CRVelocityFromState(momentum_state_local, light_speed_local,
                            state_x, state_y, state_z, vx, vy, vz);
        Real q_over_m = pr(IPM, p);
        int sp = pi(PSP, p);
        if (sp < 0 || sp >= nspecies_local) return;
        Real weight = pr(IPWT, p);
        if (weight <= 0.0) weight = 1.0;
        Real m_macro = qscale*weight*mspecies(sp);

        if (expanding_box_local) {
          state_x *= a1_start/a1_mid;
          state_y *= a2_start/a2_mid;
          if (use_vz_component) {
            state_z *= a3_start/a3_mid;
          }
          CRVelocityFromState(momentum_state_local, light_speed_local,
                              state_x, state_y, state_z, vx, vy, vz);
        }
        const Real state_x_before_em = state_x;
        const Real state_y_before_em = state_y;
        const Real state_z_before_em = state_z;

        // Drift to midpoint with old velocity.
        Real dt_half = 0.5*dt;
        Real x_mid = x + dt_half*vx/a1_mid;
        Real y_mid = y + dt_half*vy/a2_mid;
        Real z_mid = z + dt_half*vz/a3_mid;

        // Interpolate midpoint B and fluid velocity for frozen-in cE = -u x B.
        Real Bx = 0.0, By = 0.0, Bz = 0.0;
        Real Ux = 0.0, Uy = 0.0, Uz = 0.0;
        if (use_tsc) {
          InterpolateTSCFields(indcs, size_view, bcc, w0, use_mhd_fluid_velocity,
                               m, x_mid, y_mid, z_mid, Bx, By, Bz, Ux, Uy, Uz,
                               allow_2d3v_local);
        } else {
          InterpolateLinearFields(indcs, size_view, bcc, w0,
                                  use_mhd_fluid_velocity, m, x_mid, y_mid,
                                  z_mid, Bx, By, Bz, Ux, Uy, Uz,
                                  allow_2d3v_local);
        }
        if (!use_vz_component) {
          Bz = 0.0;
          Uz = 0.0;
        }
        if (expanding_box_local && no_mhd_mode) {
          const Real volume_mid = a1_mid*a2_mid*a3_mid;
          Bx *= a1_mid/volume_mid;
          By *= a2_mid/volume_mid;
          Bz *= a3_mid/volume_mid;
        }

        Real cEx = -(Uy*Bz - Uz*By);
        Real cEy = -(Uz*Bx - Ux*Bz);
        Real cEz = -(Ux*By - Uy*Bx);
        if (!use_vz_component) cEz = 0.0;

        // Half electric acceleration
        Real qdt_2m = q_over_m*dt_half;
        state_x += qdt_2m*cEx;
        state_y += qdt_2m*cEy;
        state_z += qdt_2m*cEz;

        // Magnetic rotation
        const Real inv_gamma_minus = momentum_state_local ?
            1.0/CRLorentzFactor(state_x, state_y, state_z, light_speed_local) : 1.0;
        Real tx = qdt_2m*Bx*inv_gamma_minus;
        Real ty = qdt_2m*By*inv_gamma_minus;
        Real tz = qdt_2m*Bz*inv_gamma_minus;
        Real t2 = tx*tx + ty*ty + tz*tz;
        Real rot_x = 2.0*tx/(1.0 + t2);
        Real rot_y = 2.0*ty/(1.0 + t2);
        Real rot_z = 2.0*tz/(1.0 + t2);

        Real state_px = state_x + (state_y*tz - state_z*ty);
        Real state_py = state_y + (state_z*tx - state_x*tz);
        Real state_pz = state_z + (state_x*ty - state_y*tx);

        state_x += state_py*rot_z - state_pz*rot_y;
        state_y += state_pz*rot_x - state_px*rot_z;
        state_z += state_px*rot_y - state_py*rot_x;

        // Half electric acceleration
        state_x += qdt_2m*cEx;
        state_y += qdt_2m*cEy;
        state_z += qdt_2m*cEz;
        Real feedback_state_x_before = state_x_before_em;
        Real feedback_state_y_before = state_y_before_em;
        Real feedback_state_z_before = state_z_before_em;

        if (expanding_box_local) {
          state_x *= a1_mid/a1_end;
          state_y *= a2_mid/a2_end;
          feedback_state_x_before *= a1_mid/a1_end;
          feedback_state_y_before *= a2_mid/a2_end;
          if (use_vz_component) {
            state_z *= a3_mid/a3_end;
            feedback_state_z_before *= a3_mid/a3_end;
          }
        }
        const Real feedback_energy_before =
            CRKineticEnergy(momentum_state_local, light_speed_local,
                            feedback_state_x_before, feedback_state_y_before,
                            feedback_state_z_before);
        const Real feedback_energy_after =
            CRKineticEnergy(momentum_state_local, light_speed_local,
                            state_x, state_y, state_z);
        CRVelocityFromState(momentum_state_local, light_speed_local,
                            state_x, state_y, state_z, vx, vy, vz);

        // Complete drift from midpoint to full-step position, then apply wall
        // reflection. Feedback diagnostics below remain the EM-push delta only.
        Real x_new = x_mid + dt_half*vx/a1_mid;
        Real y_new = y_mid + dt_half*vy/a2_mid;
        Real z_new = (nx3_local > 1) ? (z_mid + dt_half*vz/a3_mid) : 0.0;
        ApplyReflectiveParticleBCs(m, size_view, mb_bcs_view, multi_d_local,
                                   (nx3_local > 1), x_new, y_new, z_new,
                                   state_x, state_y, state_z);
        pr(IPX, p) = x_new;
        pr(IPY, p) = y_new;
        pr(IPZ, p) = z_new;

        // Store velocity in engineering mode and p/m in physical PIC modes.
        pr(IPVX, p) = state_x;
        pr(IPVY, p) = state_y;
        pr(IPVZ, p) = use_vz_component ? state_z : 0.0;

        // Store sampled midpoint EM fields at particle.
        pr(IPBX, p) = Bx;
        pr(IPBY, p) = By;
        pr(IPBZ, p) = use_vz_component ? Bz : 0.0;
        pr(IPEX, p) = cEx;
        pr(IPEY, p) = cEy;
        pr(IPEZ, p) = use_vz_component ? cEz : 0.0;

        pr(IPDPX, p) = m_macro*(state_x - feedback_state_x_before)*inv_dt;
        pr(IPDPY, p) = m_macro*(state_y - feedback_state_y_before)*inv_dt;
        pr(IPDPZ, p) = m_macro*(state_z - feedback_state_z_before)*inv_dt;
        pr(IPDE, p) = m_macro*(feedback_energy_after - feedback_energy_before)*inv_dt;
        pr(IPEBDOT, p) = cEx*Bx + cEy*By + cEz*Bz;
        if (deltaf_local) {
          const Real f0_current = adaptive_deltaf_local ?
              PICAdaptiveDeltaFBackgroundValue(
                  deltaf_kappa_local, deltaf_p0_local, deltaf_adaptive_p0_local,
                  deltaf_adaptive_xi_local, a1_end, a2_end, a3_end,
                  state_x, state_y, state_z) :
              PICDeltaFBackgroundValue(
                  deltaf_background_local, deltaf_p0_local, deltaf_kappa_local,
                  deltaf_drift_x1_local, deltaf_drift_x2_local,
                  deltaf_drift_x3_local, deltaf_aniso_x1_local,
                  deltaf_aniso_x2_local, deltaf_aniso_x3_local,
                  a1_end, a2_end, a3_end, state_x, state_y, state_z);
          pr(IPDFWT, p) = 1.0 - f0_current/pr(IPF0, p);
        }

        // Update displacement tracking if enabled
        if (track_displacement_local) {
          CRVelocityFromState(momentum_state_local, light_speed_local,
                              state_x, state_y, state_z, vx, vy, vz);
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
