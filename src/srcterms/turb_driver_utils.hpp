#ifndef SRCTERMS_TURB_DRIVER_UTILS_HPP_
#define SRCTERMS_TURB_DRIVER_UTILS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver_utils.hpp
//! \brief Small host-side helpers used by the stochastic turbulence driver.

#include <cmath>

namespace turbulence {

// Keep exactly one member of each {k, -k} pair.  A complex amplitude for this
// representative already supplies both members when the real force field is formed.
inline bool IsCanonicalMode(const int nkx, const int nky, const int nkz) {
  if (nkx != 0) return nkx > 0;
  if (nky != 0) return nky > 0;
  return nkz > 0;
}

inline bool IsIsotropicMode(const int nkx, const int nky, const int nkz,
                            const int nlow_sqr, const int nhigh_sqr,
                            const bool multi_d, const bool three_d) {
  if (nkx == 0 && nky == 0 && nkz == 0) return false;
  if ((!multi_d && nky != 0) || (!three_d && nkz != 0)) return false;
  if (!IsCanonicalMode(nkx, nky, nkz)) return false;

  int nsqr = nkx*nkx;
  if (multi_d) nsqr += nky*nky;
  if (three_d) nsqr += nkz*nkz;
  return nsqr >= nlow_sqr && nsqr <= nhigh_sqr;
}

inline int SpatialDimension(const bool multi_d, const bool three_d) {
  return three_d ? 3 : (multi_d ? 2 : 1);
}

inline bool HasCompleteIsotropicBounds(
    const int nhigh, const int min_kx, const int max_kx,
    const int min_ky, const int max_ky, const int min_kz, const int max_kz,
    const bool multi_d, const bool three_d) {
  if (min_kx != -nhigh || max_kx != nhigh) return false;
  if (multi_d && (min_ky != -nhigh || max_ky != nhigh)) return false;
  if (three_d && (min_kz != -nhigh || max_kz != nhigh)) return false;
  return true;
}

template <typename T>
inline T PowerLawAmplitudeExponent(const T spectrum_exponent,
                                   const int spatial_dimension) {
  return static_cast<T>(0.5)*
      (spectrum_exponent + static_cast<T>(spatial_dimension - 1));
}

template <typename T>
inline T ShellCompensationExponent(const int spatial_dimension) {
  return static_cast<T>(0.5)*static_cast<T>(spatial_dimension - 1);
}

template <typename T>
inline T BlendProjectedModeComponent(const T amplitude, const T wave_component,
                                     const T wave_dot_amplitude,
                                     const T wave_squared,
                                     const T solenoidal_fraction) {
  const T compressive = wave_component*wave_dot_amplitude/wave_squared;
  const T solenoidal = amplitude - compressive;
  return solenoidal_fraction*solenoidal +
         (static_cast<T>(1.0) - solenoidal_fraction)*compressive;
}

// Solve m0*s^2 + m1*s = dedt for the non-negative root without subtractive
// cancellation when m1 is positive.
template <typename T>
inline T PositiveConstantEdotScale(const T m0, const T m1, const T dedt) {
  if (dedt == static_cast<T>(0.0)) return static_cast<T>(0.0);
  const T disc = std::hypot(
      m1, static_cast<T>(2.0)*std::sqrt(m0)*std::sqrt(dedt));
  if (m1 >= static_cast<T>(0.0)) {
    return dedt/(static_cast<T>(0.5)*m1 + static_cast<T>(0.5)*disc);
  }
  return (-static_cast<T>(0.5)*m1)/m0 +
         (static_cast<T>(0.5)*disc)/m0;
}

}  // namespace turbulence

#endif  // SRCTERMS_TURB_DRIVER_UTILS_HPP_
