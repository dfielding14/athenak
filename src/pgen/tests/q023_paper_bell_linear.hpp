//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q023_paper_bell_linear.hpp
//! \brief Shared Q-023 Sun and Bai Section 5.2 Bell seed-carrier math.

#ifndef PGEN_TESTS_Q023_PAPER_BELL_LINEAR_HPP_
#define PGEN_TESTS_Q023_PAPER_BELL_LINEAR_HPP_

#include <cmath>

#if !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#include "athena.hpp"
#endif

#if defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#define Q023_INLINE inline
#else
#define Q023_INLINE KOKKOS_INLINE_FUNCTION
#endif

namespace q023_paper_bell_linear {

struct Vector3 {
  double x1;
  double x2;
  double x3;
};

struct Basis {
  Vector3 parallel;
  Vector3 transverse_a;
  Vector3 transverse_b;
};

struct ModeParameters {
  int dimension;
  double epsilon;
  double amplitude;
  double rho;
  double pressure;
  double b_g;
  double u_a;
  double k0;
};

struct ModeSample {
  Vector3 magnetic;
  Vector3 velocity;
};

Q023_INLINE Vector3 Add(const Vector3 a, const Vector3 b) {
  return {a.x1 + b.x1, a.x2 + b.x2, a.x3 + b.x3};
}

Q023_INLINE Vector3 Scale(const Vector3 value, const double factor) {
  return {factor*value.x1, factor*value.x2, factor*value.x3};
}

Q023_INLINE double Dot(const Vector3 a, const Vector3 b) {
  return a.x1*b.x1 + a.x2*b.x2 + a.x3*b.x3;
}

Q023_INLINE Vector3 Cross(const Vector3 a, const Vector3 b) {
  return {
    a.x2*b.x3 - a.x3*b.x2,
    a.x3*b.x1 - a.x1*b.x3,
    a.x1*b.x2 - a.x2*b.x1
  };
}

Q023_INLINE Vector3 Normalize(const Vector3 value) {
  const double inv_norm = 1.0/std::sqrt(Dot(value, value));
  return Scale(value, inv_norm);
}

Q023_INLINE Basis ModeBasis(const int dimension) {
  const Vector3 parallel = Normalize({
    1.0,
    (dimension >= 2) ? 2.0 : 0.0,
    (dimension >= 3) ? 4.0 : 0.0
  });
  const Vector3 transverse_a = (dimension == 1) ?
      Vector3{0.0, 1.0, 0.0} :
      Normalize(Vector3{-parallel.x2, parallel.x1, 0.0});
  return {parallel, transverse_a, Cross(parallel, transverse_a)};
}

Q023_INLINE double NormalizedGrowthRate(const double epsilon) {
  return std::sqrt(1.0 - epsilon*epsilon);
}

Q023_INLINE ModeSample EigenmodeAtPhase(const ModeParameters parameters,
                                        const double phase) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const double growth = NormalizedGrowthRate(parameters.epsilon);
  const Vector3 magnetic = Add(
      Scale(basis.parallel, parameters.b_g),
      Add(Scale(basis.transverse_a, parameters.amplitude*cs),
          Scale(basis.transverse_b, -parameters.amplitude*sn)));
  const double velocity_scale = parameters.amplitude/std::sqrt(parameters.rho);
  const Vector3 velocity = Add(
      Scale(basis.transverse_a,
            velocity_scale*(-parameters.epsilon*cs + growth*sn)),
      Scale(basis.transverse_b,
            velocity_scale*(growth*cs + parameters.epsilon*sn)));
  return {magnetic, velocity};
}

Q023_INLINE Vector3 VectorPotentialAt(const ModeParameters parameters,
                                     const Vector3 position) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double phase = parameters.k0*Dot(basis.parallel, position);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const Vector3 guide = Scale(Cross(basis.parallel, position),
                              0.5*parameters.b_g);
  const Vector3 perturbation = Add(
      Scale(basis.transverse_a, parameters.amplitude*cs/parameters.k0),
      Scale(basis.transverse_b, -parameters.amplitude*sn/parameters.k0));
  return Add(guide, perturbation);
}

}  // namespace q023_paper_bell_linear

#undef Q023_INLINE

#endif  // PGEN_TESTS_Q023_PAPER_BELL_LINEAR_HPP_
