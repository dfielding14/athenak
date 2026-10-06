#ifndef DIFFUSION_LIMITERS_HPP_
#define DIFFUSION_LIMITERS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \brief Shared van Leer slope means for anisotropic thermal diffusion.
#include <limits>

#include "athena.hpp"

KOKKOS_INLINE_FUNCTION
Real VanLeerLimiter(const Real a, const Real b) {
  const Real product = a*b;
  if (product > 0.0 && Kokkos::isfinite(product) &&
      product >= std::numeric_limits<Real>::min()) {
    // Preserve the original arithmetic when its intermediates are representable.
    const Real numerator = 2.0*a*b;
    const Real sum = a+b;
    if (Kokkos::isfinite(numerator) && Kokkos::isfinite(sum)) {
      return 2.0*a*b/(a+b);
    }
  }
  if (!((a > 0.0 && b > 0.0) || (a < 0.0 && b < 0.0))) return 0.0;
  if (!Kokkos::isfinite(a) || !Kokkos::isfinite(b)) return 2.0*a*b/(a+b);
  // Same-sign harmonic mean without products or sums of large slopes.
  const Real smaller = fmin(fabs(a), fabs(b));
  const Real larger = fmax(fabs(a), fabs(b));
  return copysign(smaller/(0.5 + 0.5*(smaller/larger)), a);
}

KOKKOS_INLINE_FUNCTION
Real VL4Limiter(const Real a, const Real b, const Real c, const Real d) {
  return VanLeerLimiter(VanLeerLimiter(a,b),VanLeerLimiter(c,d));
}

#endif // DIFFUSION_LIMITERS_HPP_
