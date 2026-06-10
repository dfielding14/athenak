#ifndef DIFFUSION_CGL_LANDAU_FLUID_ARITHMETIC_HPP_
#define DIFFUSION_CGL_LANDAU_FLUID_ARITHMETIC_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid_arithmetic.hpp
//! \brief Overflow-safe weighted arithmetic for the CGL Landau-fluid split update.

#include <limits>

#include "athena.hpp"

namespace cgl_lf {

struct ScaledValue {
  Real mantissa = 0.0;
  int exponent = 0;
  bool infinite = false;
  bool invalid = false;
};

KOKKOS_INLINE_FUNCTION
ScaledValue FromReal(const Real value) {
  ScaledValue result;
  if (Kokkos::isnan(value)) {
    result.mantissa = value;
    result.invalid = true;
    return result;
  }
  if (value == 0.0) {
    result.mantissa = value;
    return result;
  }
  if (Kokkos::isinf(value)) {
    result.mantissa = copysign(static_cast<Real>(1.0), value);
    result.infinite = true;
    return result;
  }
  result.mantissa = frexp(value, &result.exponent);
  return result;
}

KOKKOS_INLINE_FUNCTION
ScaledValue Multiply(ScaledValue value, const Real factor) {
  if (value.invalid || Kokkos::isnan(factor)) {
    return FromReal(std::numeric_limits<Real>::quiet_NaN());
  }
  if (factor == 0.0) {
    if (value.infinite) {
      return FromReal(std::numeric_limits<Real>::quiet_NaN());
    }
    const bool negative =
        Kokkos::signbit(value.mantissa) != Kokkos::signbit(factor);
    return FromReal(copysign(static_cast<Real>(0.0),
                             negative ? static_cast<Real>(-1.0)
                                      : static_cast<Real>(1.0)));
  }
  if (value.mantissa == 0.0) {
    if (Kokkos::isinf(factor)) {
      return FromReal(std::numeric_limits<Real>::quiet_NaN());
    }
    const bool negative =
        Kokkos::signbit(value.mantissa) != Kokkos::signbit(factor);
    value.mantissa = copysign(
        static_cast<Real>(0.0),
        negative ? static_cast<Real>(-1.0) : static_cast<Real>(1.0));
    return value;
  }
  if (Kokkos::isinf(factor)) {
    value.mantissa = copysign(
        static_cast<Real>(1.0), value.mantissa*factor);
    value.infinite = true;
    return value;
  }

  int factor_exponent = 0;
  const Real factor_mantissa = frexp(factor, &factor_exponent);
  value.mantissa *= factor_mantissa;
  value.exponent += factor_exponent;
  int product_exponent = 0;
  value.mantissa = frexp(value.mantissa, &product_exponent);
  value.exponent += product_exponent;
  return value;
}

KOKKOS_INLINE_FUNCTION
ScaledValue Add(const ScaledValue &a, const ScaledValue &b) {
  if (a.invalid || b.invalid) {
    return FromReal(std::numeric_limits<Real>::quiet_NaN());
  }
  if (a.mantissa == 0.0) return b;
  if (b.mantissa == 0.0) return a;

  if (a.infinite || b.infinite) {
    if (a.infinite && b.infinite &&
        Kokkos::signbit(a.mantissa) != Kokkos::signbit(b.mantissa)) {
      return FromReal(std::numeric_limits<Real>::quiet_NaN());
    }
    return a.infinite ? a : b;
  }

  const int common_exponent =
      (a.exponent > b.exponent) ? a.exponent : b.exponent;
  const Real aligned_a =
      scalbn(a.mantissa, a.exponent - common_exponent);
  const Real aligned_b =
      scalbn(b.mantissa, b.exponent - common_exponent);
  const Real sum = aligned_a + aligned_b;
  if (sum == 0.0) {
    return FromReal(sum);
  }

  ScaledValue result;
  int sum_exponent = 0;
  result.mantissa = frexp(sum, &sum_exponent);
  result.exponent = common_exponent + sum_exponent;
  return result;
}

KOKKOS_INLINE_FUNCTION
Real Materialize(const ScaledValue &value) {
  if (value.invalid) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  if (value.infinite) {
    return copysign(std::numeric_limits<Real>::infinity(), value.mantissa);
  }
  return scalbn(value.mantissa, value.exponent);
}

KOKKOS_INLINE_FUNCTION
Real WeightedRKL2Update(const Real muj, const Real state1,
                        const Real nuj, const Real state2,
                        const Real initial_coefficient, const Real initial_state,
                        const Real cached_coefficient, const Real cached_rhs,
                        const Real current_rhs) {
  ScaledValue value = Multiply(FromReal(state1), muj);
  value = Add(value, Multiply(FromReal(state2), nuj));
  value = Add(
      value, Multiply(FromReal(initial_state), initial_coefficient));
  value = Add(value, Multiply(FromReal(cached_rhs), cached_coefficient));
  value = Add(value, FromReal(current_rhs));
  return Materialize(value);
}

KOKKOS_INLINE_FUNCTION
Real LimitedRatio(const Real ratio) {
  const Real ratio_abs = fabs(ratio);
  if (ratio_abs <= 1.0) {
    return ratio/(1.0 + ratio_abs);
  }
  const Real magnitude = 1.0/(1.0 + 1.0/ratio_abs);
  return copysign(magnitude, ratio);
}

KOKKOS_INLINE_FUNCTION
ScaledValue LimitedHeatFlux(const Real ratio, const Real cparallel,
                            const Real normalization, const Real pressure) {
  ScaledValue value = FromReal(LimitedRatio(ratio));
  value = Multiply(value, cparallel);
  value = Multiply(value, normalization);
  return Multiply(value, pressure);
}

inline Real FirstStageRKLWeight(const int nstages) {
  if (nstages == 1) {
    return 1.0;
  }
  const Real stages = static_cast<Real>(nstages);
  return (static_cast<Real>(4.0)/static_cast<Real>(3.0))/
         (stages*stages + stages - static_cast<Real>(2.0));
}

inline Real CachedRHSCoefficient(const Real gamma_tilde, const int nstages) {
  return gamma_tilde/FirstStageRKLWeight(nstages);
}

} // namespace cgl_lf

#endif // DIFFUSION_CGL_LANDAU_FLUID_ARITHMETIC_HPP_
