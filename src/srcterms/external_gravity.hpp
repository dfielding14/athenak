#ifndef SRCTERMS_EXTERNAL_GRAVITY_HPP_
#define SRCTERMS_EXTERNAL_GRAVITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file external_gravity.hpp
//! \brief Standard external gravitational potentials for fluid source terms.

#include <cmath>

#include "athena.hpp"

class MeshBlockPack;
class ParameterInput;

namespace external_gravity {

enum Model {
  none = 0,
  uniform = 1,
  cartesian_harmonic = 2,
  spherical_harmonic = 3,
  point_mass = 4,
  plummer = 5,
  nfw = 6,
  logarithmic_halo = 7,
  miyamoto_nagai = 8,
  outer_cgm = 9,
  gotham = 10,
  user = 11
};

enum BoundaryVelocityMode {
  copy_velocity = 0,
  zero_normal_velocity = 1,
  no_inflow_velocity = 2,
  reflect_normal_velocity = 3
};

struct ExternalGravityData {
  int model = Model::none;
  Real grav_constant = 1.0;

  Real mesh_x1min = 0.0;
  Real mesh_x1max = 1.0;
  Real mesh_x2min = 0.0;
  Real mesh_x2max = 1.0;
  Real mesh_x3min = 0.0;
  Real mesh_x3max = 1.0;

  Real x1_origin = 0.0;
  Real x2_origin = 0.0;
  Real x3_origin = 0.0;

  Real g1 = 0.0;
  Real g2 = 0.0;
  Real g3 = 0.0;

  Real omega1 = 0.0;
  Real omega2 = 0.0;
  Real omega3 = 0.0;
  Real omega = 0.0;

  Real mass = 1.0;
  Real softening = 0.0;

  Real r_scale = 1.0;
  Real rho_scale = 1.0;

  Real v0 = 1.0;
  Real core_radius = 1.0;
  Real flattening_q = 1.0;

  Real disk_mass = 1.0;
  Real disk_a = 1.0;
  Real disk_b = 0.1;

  Real rho_mean = 0.0;
  Real r_outer = 1.0;

  int taper_gravity = 0;
  Real taper_width_x1_inner = 0.0;
  Real taper_width_x1_outer = 0.0;
  Real taper_width_x2_inner = 0.0;
  Real taper_width_x2_outer = 0.0;
  Real taper_width_x3_inner = 0.0;
  Real taper_width_x3_outer = 0.0;

  Real user_fd_step = 1.0e-5;

  int boundary_velocity = BoundaryVelocityMode::no_inflow_velocity;
  Real boundary_sound_speed = -1.0;
  Real boundary_density_floor = -1.0;
  Real boundary_pressure_floor = -1.0;
  Real boundary_max_exponent = 60.0;
};

ExternalGravityData ParseInput(MeshBlockPack *pp, ParameterInput *pin);

KOKKOS_FUNCTION
Real UserExternalGravityPotential(Real x1, Real x2, Real x3);

KOKKOS_INLINE_FUNCTION
Real RadiusRatioLog1p(const Real r, const Real r_scale) {
  const Real x = r/fmax(r_scale, 1.0e-20);
  if (x < 1.0e-8) {
    return 1.0 - 0.5*x + ONE_3RD*x*x;
  }
  return log1p(x)/x;
}

KOKKOS_INLINE_FUNCTION
Real NFWMassFactor(const Real x) {
  if (x < 1.0e-6) {
    return 0.5*x*x - (2.0/3.0)*x*x*x + 0.75*x*x*x*x;
  }
  return log1p(x) - x/(1.0 + x);
}

KOKKOS_INLINE_FUNCTION
Real Potential(const ExternalGravityData &grav, const Real x1, const Real x2,
               const Real x3) {
  const Real x = x1 - grav.x1_origin;
  const Real y = x2 - grav.x2_origin;
  const Real z = x3 - grav.x3_origin;
  const Real r_cyl2 = fma(x, x, y*y);
  const Real r_cyl = sqrt(r_cyl2);
  const Real r2 = fma(z, z, r_cyl2);
  const Real r = sqrt(fmax(r2, 1.0e-20));
  const Real gconst = grav.grav_constant;

  switch (grav.model) {
    case Model::uniform:
      return -(grav.g1*x + grav.g2*y + grav.g3*z);

    case Model::cartesian_harmonic:
      return 0.5*(SQR(grav.omega1)*SQR(x) + SQR(grav.omega2)*SQR(y)
                  + SQR(grav.omega3)*SQR(z));

    case Model::spherical_harmonic:
      return 0.5*SQR(grav.omega)*r2;

    case Model::point_mass:
    case Model::plummer:
      return -gconst*grav.mass/sqrt(r2 + SQR(grav.softening));

    case Model::nfw:
      return -4.0*M_PI*gconst*grav.rho_scale*SQR(grav.r_scale)*
             RadiusRatioLog1p(r, grav.r_scale);

    case Model::logarithmic_halo:
      return 0.5*SQR(grav.v0)*log(SQR(grav.core_radius) + r_cyl2
                                  + SQR(z/grav.flattening_q));

    case Model::miyamoto_nagai:
      return -gconst*grav.disk_mass/
             sqrt(r_cyl2 + SQR(grav.disk_a + sqrt(SQR(z) + SQR(grav.disk_b))));

    case Model::outer_cgm: {
      const Real c_outer = (4.0/3.0)*pow(grav.r_outer, 1.5);
      return 4.0*M_PI*gconst*grav.rho_mean*(c_outer*sqrt(r) + (1.0/6.0)*r2);
    }

    case Model::gotham: {
      const Real phi_nfw = -4.0*M_PI*gconst*grav.rho_scale*SQR(grav.r_scale)*
                           RadiusRatioLog1p(r, grav.r_scale);
      const Real phi_disk = -gconst*grav.disk_mass/
          sqrt(r_cyl2 + SQR(grav.disk_a + sqrt(SQR(z) + SQR(grav.disk_b))));
      const Real c_outer = (4.0/3.0)*pow(grav.r_outer, 1.5);
      const Real phi_outer = 4.0*M_PI*gconst*grav.rho_mean*
                             (c_outer*sqrt(r) + (1.0/6.0)*r2);
      return phi_nfw + phi_disk + phi_outer;
    }

    case Model::user:
      return UserExternalGravityPotential(x1, x2, x3);

    default:
      return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
Real SmoothStep01(const Real s) {
  const Real t = fmin(1.0, fmax(0.0, s));
  return t*t*(3.0 - 2.0*t);
}

KOKKOS_INLINE_FUNCTION
Real BoundaryTaper(const Real x, const Real xmin, const Real xmax,
                   const Real w_inner, const Real w_outer) {
  Real taper = 1.0;
  if (w_inner > 0.0) {
    taper *= SmoothStep01((x - xmin)/w_inner);
  }
  if (w_outer > 0.0) {
    taper *= SmoothStep01((xmax - x)/w_outer);
  }
  return taper;
}

KOKKOS_INLINE_FUNCTION
Real GravityTaper(const ExternalGravityData &grav, const Real x1, const Real x2,
                  const Real x3) {
  if (grav.taper_gravity == 0) return 1.0;
  return BoundaryTaper(x1, grav.mesh_x1min, grav.mesh_x1max,
                       grav.taper_width_x1_inner, grav.taper_width_x1_outer)*
         BoundaryTaper(x2, grav.mesh_x2min, grav.mesh_x2max,
                       grav.taper_width_x2_inner, grav.taper_width_x2_outer)*
         BoundaryTaper(x3, grav.mesh_x3min, grav.mesh_x3max,
                       grav.taper_width_x3_inner, grav.taper_width_x3_outer);
}

KOKKOS_INLINE_FUNCTION
void FiniteDifferenceAcceleration(const ExternalGravityData &grav, const Real x1,
                                  const Real x2, const Real x3,
                                  Real &a1, Real &a2, Real &a3) {
  const Real h = fmax(grav.user_fd_step, 1.0e-12);
  a1 = -(Potential(grav, x1 + h, x2, x3) - Potential(grav, x1 - h, x2, x3))/(2.0*h);
  a2 = -(Potential(grav, x1, x2 + h, x3) - Potential(grav, x1, x2 - h, x3))/(2.0*h);
  a3 = -(Potential(grav, x1, x2, x3 + h) - Potential(grav, x1, x2, x3 - h))/(2.0*h);
}

KOKKOS_INLINE_FUNCTION
void Acceleration(const ExternalGravityData &grav, const Real x1, const Real x2,
                  const Real x3, Real &a1, Real &a2, Real &a3) {
  const Real x = x1 - grav.x1_origin;
  const Real y = x2 - grav.x2_origin;
  const Real z = x3 - grav.x3_origin;
  const Real r_cyl2 = fma(x, x, y*y);
  const Real r2 = fma(z, z, r_cyl2);
  const Real r = sqrt(r2);
  const Real gconst = grav.grav_constant;

  a1 = 0.0;
  a2 = 0.0;
  a3 = 0.0;

  switch (grav.model) {
    case Model::uniform:
      a1 = grav.g1;
      a2 = grav.g2;
      a3 = grav.g3;
      break;

    case Model::cartesian_harmonic:
      a1 = -SQR(grav.omega1)*x;
      a2 = -SQR(grav.omega2)*y;
      a3 = -SQR(grav.omega3)*z;
      break;

    case Model::spherical_harmonic:
      a1 = -SQR(grav.omega)*x;
      a2 = -SQR(grav.omega)*y;
      a3 = -SQR(grav.omega)*z;
      break;

    case Model::point_mass:
    case Model::plummer: {
      const Real denom = pow(r2 + SQR(grav.softening), 1.5);
      const Real coef = (denom > 0.0) ? -gconst*grav.mass/denom : 0.0;
      a1 = coef*x;
      a2 = coef*y;
      a3 = coef*z;
      break;
    }

    case Model::nfw: {
      if (r > 0.0) {
        const Real rs = fmax(grav.r_scale, 1.0e-20);
        const Real mf = NFWMassFactor(r/rs);
        const Real coef = -4.0*M_PI*gconst*grav.rho_scale*rs*rs*rs*mf/(r*r*r);
        a1 = coef*x;
        a2 = coef*y;
        a3 = coef*z;
      }
      break;
    }

    case Model::logarithmic_halo: {
      const Real q = fmax(grav.flattening_q, 1.0e-20);
      const Real denom = SQR(grav.core_radius) + r_cyl2 + SQR(z/q);
      const Real v02 = SQR(grav.v0);
      a1 = -v02*x/denom;
      a2 = -v02*y/denom;
      a3 = -v02*z/(q*q*denom);
      break;
    }

    case Model::miyamoto_nagai: {
      const Real zzb = sqrt(SQR(z) + SQR(grav.disk_b));
      const Real azb = grav.disk_a + zzb;
      const Real denom = pow(r_cyl2 + SQR(azb), 1.5);
      const Real coef = (denom > 0.0) ? -gconst*grav.disk_mass/denom : 0.0;
      a1 = coef*x;
      a2 = coef*y;
      a3 = (zzb > 0.0) ? coef*azb*z/zzb : 0.0;
      break;
    }

    case Model::outer_cgm: {
      if (r > 0.0) {
        const Real c_outer = (4.0/3.0)*pow(grav.r_outer, 1.5);
        const Real dphi_dr = 4.0*M_PI*gconst*grav.rho_mean*
                             (0.5*c_outer/sqrt(r) + ONE_3RD*r);
        const Real coef = -dphi_dr/r;
        a1 = coef*x;
        a2 = coef*y;
        a3 = coef*z;
      }
      break;
    }

    case Model::gotham: {
      Real nfw_a1 = 0.0, nfw_a2 = 0.0, nfw_a3 = 0.0;
      if (r > 0.0) {
        const Real rs = fmax(grav.r_scale, 1.0e-20);
        const Real mf = NFWMassFactor(r/rs);
        const Real coef = -4.0*M_PI*gconst*grav.rho_scale*rs*rs*rs*mf/(r*r*r);
        nfw_a1 = coef*x;
        nfw_a2 = coef*y;
        nfw_a3 = coef*z;
      }

      const Real zzb = sqrt(SQR(z) + SQR(grav.disk_b));
      const Real azb = grav.disk_a + zzb;
      const Real disk_denom = pow(r_cyl2 + SQR(azb), 1.5);
      const Real disk_coef = (disk_denom > 0.0) ? -gconst*grav.disk_mass/disk_denom : 0.0;

      Real outer_a1 = 0.0, outer_a2 = 0.0, outer_a3 = 0.0;
      if (r > 0.0) {
        const Real c_outer = (4.0/3.0)*pow(grav.r_outer, 1.5);
        const Real dphi_dr = 4.0*M_PI*gconst*grav.rho_mean*
                             (0.5*c_outer/sqrt(r) + ONE_3RD*r);
        const Real outer_coef = -dphi_dr/r;
        outer_a1 = outer_coef*x;
        outer_a2 = outer_coef*y;
        outer_a3 = outer_coef*z;
      }

      a1 = nfw_a1 + disk_coef*x + outer_a1;
      a2 = nfw_a2 + disk_coef*y + outer_a2;
      a3 = nfw_a3 + ((zzb > 0.0) ? disk_coef*azb*z/zzb : 0.0) + outer_a3;
      break;
    }

    case Model::user:
      FiniteDifferenceAcceleration(grav, x1, x2, x3, a1, a2, a3);
      break;

    default:
      break;
  }

  const Real taper = GravityTaper(grav, x1, x2, x3);
  a1 *= taper;
  a2 *= taper;
  a3 *= taper;
}

} // namespace external_gravity

#endif  // SRCTERMS_EXTERNAL_GRAVITY_HPP_
