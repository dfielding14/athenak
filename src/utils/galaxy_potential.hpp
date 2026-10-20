#ifndef UTILS_GALAXY_POTENTIAL_HPP_
#define UTILS_GALAXY_POTENTIAL_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file galaxy_potential.hpp
//! \brief Static GOTHAM galaxy potential: NFW halo + Miyamoto-Nagai disk + outer term.
//! The one copy used by both the gas (pgen/gotham.cpp) and the rk4_gravity star-particle
//! pusher, so the two feel the same gravity. The legacy pgen/cgm_cooling_flow_full.cpp
//! keeps its own copy.

#include <cmath>

#include "athena.hpp"

//! \fn Real GalaxyPotOuterCoef(const Real r_200)
//! \brief c_outer = (4/3)(5 r_200)^(3/2): the outer term of GalaxyPot is the potential of
//! the density rho_mean [1 + (r/(5 r_200))^(-3/2)].
KOKKOS_INLINE_FUNCTION
Real GalaxyPotOuterCoef(const Real r_200) {
  return (4.0/3.0) * pow(5 * r_200, 1.5);
}

//! \fn Real GalaxyPot()
//! \brief Potential at (x1, x2, x3) in code units, with c_outer from GalaxyPotOuterCoef.
KOKKOS_INLINE_FUNCTION
Real GalaxyPot(Real x1, Real x2, Real x3,
               Real G, Real r_s, Real rho_s,
               Real M_gal, Real a_gal, Real z_gal,
               Real c_outer, Real rho_mean) {
  const Real R2 = fma(x1, x1 , x2*x2);
  const Real r2 = fma(x3 , x3 , R2);
  const Real r  = sqrt(fmax(r2, 1e-20));

  // NFW component
  Real x = r / r_s;
  Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log1p(x) / x;

  // Miyamoto-Nagai model
  Real phi_MN = -G * M_gal / sqrt(R2 + SQR(sqrt(fma(x3 , x3 , z_gal*z_gal)) + a_gal));

  // Outer component
  Real phi_Outer = 4 * M_PI * G * rho_mean * (c_outer * sqrt(r) + (1.0/6.0) * r2);

  // Total potential
  return phi_NFW + phi_MN + phi_Outer;
}

#endif // UTILS_GALAXY_POTENTIAL_HPP_
