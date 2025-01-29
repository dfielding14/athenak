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

// Forward declare helper functions
//KOKKOS_INLINE_FUNCTION size_t find_he_bin(Real val);

//----------------------------------------------------------------------------------------
//! \fn Real CGMCoolFn()
//! \brief WiersmaCooling at redshift z = 0 taken from Wiersma et al (2009)
// Inputs are Temperature in K, Metallicity in Zsun, and density in cm^-3
KOKKOS_INLINE_FUNCTION
Real CGMCoolFn(Real temp, Real nH, Real Z) {
  return 1.0;
  /*
  // If out of range, no cooling
  if (temp < Tbins(0) or temp > Tbins(Tbins_TOTAL_SIZE - 1)) return 0.;
  if (nH < nHbins(0)  or nH > nHbins(nHbins_TOTAL_SIZE - 1)) return 0.;
  
  
  // Compute Helium Mass Fraction
  Real He2Habundance = pow(10,-1.07) * (0.71553 + 0.28447*Z); // Asplund+09, Groves+04
  Real X = (1 - 0.014*Z) / (1.+4.*He2Habundance);
  Real Y = 4.*He2Habundance * X;
  size_t iHe = find_he_bin(Y);
   

  // Convert input values to log space
  Real log_temp = log10(temp);
  Real log_density = log10(nH);

  // Locate indices in Tbins and nHbins
  size_t i = 0, j = 0;
  while (i < Tbins_TOTAL_SIZE - 1 && log10(Tbins(i + 1)) < log_temp) ++i;
  while (j < nHbins_TOTAL_SIZE - 1 && log10(nHbins(j + 1)) < log_density) ++j;

  
  // Logarithms of the Tbins and nHbins bounding the point
  Real log_T0 = log10(Tbins(i));
  Real log_T1 = log10(Tbins(i + 1));
  Real log_nH0 = log10(nHbins(j));
  Real log_nH1 = log10(nHbins(j + 1));

  // Compute weights for bilinear interpolation
  Real t = (log_temp - log_T0) / (log_T1 - log_T0);
  Real u = (log_density - log_nH0) / (log_nH1 - log_nH0);

  // Corner values in log space from the H_He cooling grid
  Real log_C00 = log10(H_He_Cooling(iHe, i, j));
  Real log_C10 = log10(H_He_Cooling(iHe, i + 1, j));
  Real log_C01 = log10(H_He_Cooling(iHe, i, j + 1));
  Real log_C11 = log10(H_He_Cooling(iHe, i + 1, j + 1));

  // Corner values in log space from the Metal cooling grid
  Real log_M00 = log10(Metal_Cooling(i, j));
  Real log_M10 = log10(Metal_Cooling(i + 1, j));
  Real log_M01 = log10(Metal_Cooling(i, j + 1));
  Real log_M11 = log10(Metal_Cooling(i + 1, j + 1));

  // Bilinear interpolation in log space
  Real log_interpolated_prim_cooling =
      (1 - t) * (1 - u) * log_C00 +
      t * (1 - u) * log_C10 +
      (1 - t) * u * log_C01 +
      t * u * log_C11;
  Real log_interpolated_metal_cooling =
      (1 - t) * (1 - u) * log_M00 +
      t * (1 - u) * log_M10 +
      (1 - t) * u * log_M01 +
      t * u * log_M11;

  // Convert back to linear space and return
  Real prim_cooling = pow(10.0, log_interpolated_prim_cooling);
  Real metal_cooling = pow(10.0, log_interpolated_metal_cooling);

  return prim_cooling + Z*metal_cooling;
  
  return 1.0;
  */
}

/*
KOKKOS_INLINE_FUNCTION
size_t find_he_bin(Real val) {
    // Check if val is outside the bounds of the He_mf_bins array
    if (val < He_mf_bins(0)) {
        return 0;  // Return the first index
    }
    if (val > He_mf_bins(He_mf_bins_TOTAL_SIZE - 1)) {
        return He_mf_bins_TOTAL_SIZE - 1;  // Return the last index
    }

    // Calculate the step size (dx) between consecutive elements
    Real dx = He_mf_bins(1) - He_mf_bins(0);

    // Calculate the bin index and return it
    // We assume that the spacing in He_mf_bins is regular
    return static_cast<size_t>((val + (dx / 2) - He_mf_bins(0)) / dx);
}
*/

#endif // SRCTERMS_CGMCOOLING_HPP_
