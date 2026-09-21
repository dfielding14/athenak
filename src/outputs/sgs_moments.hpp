//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
#ifndef OUTPUTS_SGS_MOMENTS_HPP_
#define OUTPUTS_SGS_MOMENTS_HPP_

#include "athena.hpp"

// Raw moments in output order. Padded lanes in grouped reductions contribute zero.
KOKKOS_INLINE_FUNCTION
Real HydroSGS3DMoment(int n, const MHDCons1D &s) {
  if (n >= 10) { return 0.0; }
  if (n == 0) { return s.d; }
  if (n < 4) { return (n == 1) ? s.mx : ((n == 2) ? s.my : s.mz); }
  int t = n - 4;
  int a = (t < 3) ? 0 : ((t < 5) ? 1 : 2);
  int b = (t < 3) ? t : ((t < 5) ? t - 2 : 2);
  Real ma = (a == 0) ? s.mx : ((a == 1) ? s.my : s.mz);
  Real mb = (b == 0) ? s.mx : ((b == 1) ? s.my : s.mz);
  return ma*mb/s.d;
}

KOKKOS_INLINE_FUNCTION
Real MHDSGSMoment(int n, const MHDCons1D &s) {
  if (n >= 59) { return 0.0; }
  if (n == 0) { return s.d; }
  if (n < 4) { return (n == 1) ? s.mx : ((n == 2) ? s.my : s.mz); }
  if (n == 4) { return s.e; }  // Conserved total energy, including kinetic and magnetic.
  if (n < 8) { return (n == 5) ? s.bx : ((n == 6) ? s.by : s.bz); }
  if (n < 14) { return HydroSGS3DMoment(n - 4, s); }
  if (n < 20) {
    int t = n - 14;
    int a = (t < 3) ? 0 : ((t < 5) ? 1 : 2);
    int b = (t < 3) ? t : ((t < 5) ? t - 2 : 2);
    Real ba = (a == 0) ? s.bx : ((a == 1) ? s.by : s.bz);
    Real bb = (b == 0) ? s.bx : ((b == 1) ? s.by : s.bz);
    return ba*bb;
  }
  if (n >= 29 && n < 32) {
    Real ma = (n == 29) ? s.mx : ((n == 30) ? s.my : s.mz);
    return ma*s.e/s.d;
  }
  int t = (n < 29) ? n - 20 : ((n < 41) ? n - 32 : ((n < 50) ? n - 41 : n - 50));
  int a = t/3, b = t%3;
  Real ma = (a == 0) ? s.mx : ((a == 1) ? s.my : s.mz);
  Real bb = (b == 0) ? s.bx : ((b == 1) ? s.by : s.bz);
  if (n < 29) { return ma*bb/s.d; }
  if (n < 41) {
    Real mb = (b == 0) ? s.mx : ((b == 1) ? s.my : s.mz);
    return ma*mb*mb/s.d/s.d;
  }
  if (n < 50) { return ma*bb*bb/s.d; }
  Real ba = (a == 0) ? s.bx : ((a == 1) ? s.by : s.bz);
  return ma*ba*bb/s.d;
}

#endif  // OUTPUTS_SGS_MOMENTS_HPP_
