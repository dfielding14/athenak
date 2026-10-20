#ifndef SRCTERMS_CGM_COOLING_TERMS_HPP_
#define SRCTERMS_CGM_COOLING_TERMS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgm_cooling_terms.hpp
//! \brief The two halves of a CGM net-rate evaluation (cgm_cooling::NetRate): the terms
//! that depend only on temperature and those that depend only on the cell. Kept apart
//! from cgm_cooling.hpp so that srcterms.hpp can hold tables of TempTerms without
//! including the cooling tables.

#include "athena.hpp"

namespace cgm_cooling {

//! Terms that depend only on temperature (cgm_cooling::TempPart)
struct TempTerms {
  Real temp;          // temperature they were evaluated at [K]
  int iT;             // cooling-table bin holding temp (clamped to the table)
  Real t;             // position in that bin, linear in log T
  Real in_table;      // 1 for Tfloor <= temp <= Tceil, else 0
  Real below_table;   // 1 for temp < Tfloor (low-T fit), else 0
  Real prim_CIE;      // CIE H+He Lambda [erg cm^3/s]
  Real metal_CIE;     // CIE metal Lambda at solar metallicity [erg cm^3/s]
  Real lowT;          // Koyama & Inutsuka (2002) Lambda at solar metallicity [erg cm^3/s]
  Real heat_damp;     // disk-heating suppression: (1e4/T)^8 above 1e4 K, else 1
  Real neutral_frac;  // neutral fraction in the shielding estimate
};

//! Terms that depend only on the cell: density, metallicity, position and width
//! (cgm_cooling::CellPart)
struct CellTerms {
  Real X_rho;   // hydrogen mass density X*rho [code]
  Real nH;      // hydrogen number density [cm^-3]
  Real Z;       // metallicity [solar]
  int jN;       // PIE-table nH bin holding nH (clamped to the table)
  Real u;       // position in that bin, linear in log nH
  Real in_nH;   // 1 for nH >= nHfloor, else 0
  Real gamma0;  // disk heating before heat_damp (0 if heating is off)
  Real dx_cgs;  // cell width [cm], the shielding length
};

}  // namespace cgm_cooling

#endif  // SRCTERMS_CGM_COOLING_TERMS_HPP_
