#ifndef SRCTERMS_CGMCOOLING_HPP_
#define SRCTERMS_CGMCOOLING_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgmcooling.hpp
//! \brief function to implement CGM cooling

// Athena++ headers
#include "athena.hpp"
#include "cooling_tables.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real CGMCoolFn()
//! \brief WiersmaCooling at redshift z = 0 taken from Wiersma et al (2009)

KOKKOS_INLINE_FUNCTION
Real CGMCoolFn(Real temp, Real nH, Real Z) {
  
  // Real logt = log10(temp);

  Real He2Habundance = 10**-1.07 * (0.71553 + 0.28447*Z); // Asplund+09, Groves+04
  Real X = (1 - 0.014*Z2Zsun) / (1.+4.*He2Habundance);
  Real Y = 4.*He2Habundance * X;
  Real iHe = find_he_bin(Y);

  // linear interpolation of cooling tables
  // int ipps  = static_cast<int>(25.0*logt) - 103;
  // ipps = (ipps < 100)? ipps : 100;
  //ipps = (ipps > 0 )? ipps : 0;
  //Real x0 = 4.12 + 0.04*static_cast<Real>(ipps);
  //Real dx = logt - x0;
  //Real logcool = (lhd[ipps+1]*dx - lhd[ipps]*(dx - 0.04))*25.0;
  //return pow(10.0,logcool);
  return 42.0;
}

KOKKOS_INLINE_FUNCTION
int find_he_bin(Real val) {
    // Check if val is outside the bounds of the He_mf_bins array
    if (val < He_mf_bins[0]) {
        return 0;  // Return the first index
    }
    if (val > He_mf_bins[He_mf_bins_TOTAL_SIZE - 1]) {
        return He_mf_bins_TOTAL_SIZE - 1;  // Return the last index
    }

    // Calculate the step size (dx) between consecutive elements
    Real dx = He_mf_bins[1] - He_mf_bins[0];

    // Calculate the bin index and return it
    // We assume that the spacing in He_mf_bins is regular
    return static_cast<int>((val + (dx / 2) - He_mf_bins[0]) / dx);
}

#endif // SRCTERMS_CGMCOOLING_HPP_
