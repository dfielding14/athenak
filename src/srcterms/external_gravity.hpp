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

struct ExternalGravityData {
  int model = Model::none;
  Real grav_constant = 1.0;

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

} // namespace external_gravity

#endif  // SRCTERMS_EXTERNAL_GRAVITY_HPP_
