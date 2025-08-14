#ifndef SRCTERMS_TURB_DRIVER_HPP_
#define SRCTERMS_TURB_DRIVER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver.hpp
//  \brief defines turbulence driver class, which implements data and functions for
//  randomly forced turbulence which evolves via an Ornstein-Uhlenbeck stochastic process

#include <memory>
#include <cmath>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/random.hpp"

//----------------------------------------------------------------------------------------
//! \class TurbulenceDriver

// add enum and forward decl before class
enum class TurbBasis { Cartesian = 0, SphericalFB = 1 };

class TurbulenceDriver {
 public:
  TurbulenceDriver(MeshBlockPack *pp, ParameterInput *pin);
  ~TurbulenceDriver();

  DvceArray5D<Real> force, force_tmp1, force_tmp2;  // arrays used for turb forcing
  RNG_State rstate;                    // random state

  DualArray2D<Real> aka, akb; //to store amplitude coefficients
  DualArray1D<Real> kx_mode, ky_mode, kz_mode;
  DvceArray3D<Real> xcos, xsin, ycos, ysin, zcos, zsin;

  /* ===================== Spherical Fourier–Bessel basis data ===================== */
  TurbBasis basis_type_ = TurbBasis::Cartesian;  // basis switch

  // Quantum numbers for each mode when using SFB basis
  DualArray1D<int> l_mode, m_mode, n_mode;
  DualArray1D<Real> kln_mode;  // radial wavenumber = x_ln / r0_turb
  DualArray1D<Real> xln_root;  // spherical bessel roots

  // User-supplied limits (only used if basis_type_ == SphericalFB)
  int lmax = 0;   // maximum l
  int nmax = 0;   // maximum radial index n
  Real r0_turb = -1.0; // outer radius of driven region
  
  // For efficient mode selection
  Real kmin_sfb = 0.0;  // minimum wavenumber
  Real kmax_sfb = 0.0;  // maximum wavenumber

  // Pre-computed jl*Ylm per mode and cell (real and imaginary parts)
  DvceArray4D<Real> sfb_basis_real; // (mode,nk,nj,ni)
  DvceArray4D<Real> sfb_basis_imag;
  
  // Vector spherical harmonics components for proper divergence-free projection
  DvceArray6D<Real> sfb_vector_basis_real; // (nmb,mode,dir,nk,nj,ni) 
  DvceArray6D<Real> sfb_vector_basis_imag;

  // -----------------------------------------------------------------------------
  // Helper functions – defined in turb_driver.cpp but declared here so that they
  // can be called from device lambdas using KOKKOS_INLINE_FUNCTION wrappers.

  KOKKOS_INLINE_FUNCTION static Real sph_bessel_jl(int l, Real x);
  KOKKOS_INLINE_FUNCTION static Real legendre_Plm(int l, int m, Real x);
  KOKKOS_INLINE_FUNCTION static Real ylm_norm(int l, int m);
  KOKKOS_INLINE_FUNCTION static Real sph_bessel_jl_prime(int l, Real x);

  // Pre-compute lookup tables (host side)
  void BuildCartesianModeTable(ParameterInput *pin);
  void BuildSFBModeTable(ParameterInput *pin);
  void ComputeSphericalBesselRoots();
  static Real FindSphericalBesselRoot(int l, int n);

  // -----------------------------------------------------------------------------

  // parameters of driving
  int nlow, nhigh, spect_form;
  int mode_count;
  int rseed;  // random seed for turbulence driving
  Real kpeak;
  Real tcorr, dedt, tdriv_duration, tdriv_start;
  Real expo, exp_prl, exp_prp;
  int driving_type, turb_flag;
  int min_kz, max_kz, min_kx, max_kx, min_ky, max_ky;
  Real sol_fraction; // To store fraction of energy in solenoidal modes
  Real dt_turb_update,dt_turb_thresh;
  // Real t_last_update;
  int n_turb_updates_yet;

  // drive with constant edot or constant acceleration
  bool constant_edot;

  // spatially varying driving
  Real x_turb_scale_height, y_turb_scale_height, z_turb_scale_height;
  Real x_turb_center, y_turb_center, z_turb_center;


  // functions
  void IncludeInitializeModesTask(std::shared_ptr<TaskList> tl, TaskID start);
  void IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start);
  TaskStatus InitializeModes(Driver *pdrive, int stage);
  TaskStatus UpdateForcing(Driver *pdrive, int stage);
  TaskStatus AddForcing(Driver *pdrive, int stage);
  void Initialize();
  void ResizeArrays(int new_nmb);  // Resize arrays when mesh changes

 private:
  bool first_time = true;   // flag to enable initialization on first call
  MeshBlockPack *pmy_pack;  // ptr to MeshBlockPack containing this TurbulenceDriver
  int current_nmb = 0;      // current number of MeshBlocks for array sizing
};

// ========================== Inline math utilities =============================
KOKKOS_INLINE_FUNCTION
static Real sph_bessel_jl_rec(int l, Real x) {
  if (l == 0) return (x != 0.0 ? sin(x)/x : 1.0);
  if (l == 1) return (x != 0.0 ? sin(x)/(x*x) - cos(x)/x : 0.0);
  Real jlm2 = sin(x)/x;                         // j_0(x)
  Real jlm1 = sin(x)/(x*x) - cos(x)/x;          // j_1(x)
  for (int n=2; n<=l; ++n) {
    Real jl = ((2.0*n-1.0)/x)*jlm1 - jlm2;
    jlm2 = jlm1;
    jlm1 = jl;
  }
  return jlm1;
}

KOKKOS_INLINE_FUNCTION
Real TurbulenceDriver::sph_bessel_jl(int l, Real x) {
  return sph_bessel_jl_rec(l, x);
}

// Derivative of spherical Bessel function
KOKKOS_INLINE_FUNCTION
Real TurbulenceDriver::sph_bessel_jl_prime(int l, Real x) {
  if (x == 0.0) {
    return (l == 1) ? 1.0/3.0 : 0.0;
  }
  // Use recurrence: j'_l(x) = j_{l-1}(x) - (l+1)/x * j_l(x)
  if (l == 0) {
    return -sph_bessel_jl(1, x);
  }
  return sph_bessel_jl(l-1, x) - (l+1.0)/x * sph_bessel_jl(l, x);
}

// Associated Legendre polynomial P_l^m using forward recurrence
KOKKOS_INLINE_FUNCTION
Real TurbulenceDriver::legendre_Plm(int l, int m, Real x) {
  if (m < 0) m = -m; // use |m|; caller responsible for Condon-Shortley phase
  // compute P_m^m
  Real pmm = 1.0;
  if (m > 0) {
    Real somx2 = sqrt((1.0 - x)*(1.0 + x));
    Real fact = 1.0;
    for (int i=1; i<=m; ++i) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m) return pmm;
  // compute P_{m+1}^m
  Real pmmp1 = x * (2.0*m + 1.0) * pmm;
  if (l == m+1) return pmmp1;
  // upward recurrence
  Real pll = 0.0;
  for (int ll = m + 2; ll <= l; ++ll) {
    pll = ((2.0*ll - 1.0)*x*pmmp1 - (ll + m - 1.0)*pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

// Normalisation factor of spherical harmonic Y_l^m (real/imag parts share same norm)
KOKKOS_INLINE_FUNCTION
Real TurbulenceDriver::ylm_norm(int l, int m) {
  if (m < 0) m = -m;
  auto factorial = [](int n)->Real{
    Real f=1.0; for(int i=2;i<=n;++i) f*=i; return f; };
  Real num = (2.0*l + 1.0) * factorial(l - m);
  Real den = 4.0*M_PI * factorial(l + m);
  return sqrt(num / den);
}

#endif  // SRCTERMS_TURB_DRIVER_HPP_
