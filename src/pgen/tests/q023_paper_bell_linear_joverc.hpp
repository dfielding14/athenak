//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q023_paper_bell_linear_joverc.hpp
//! \brief Corrected Bai et al. Section 5.2 unstable Bell eigenmode carrier.

#ifndef PGEN_TESTS_Q023_PAPER_BELL_LINEAR_JOVERC_HPP_
#define PGEN_TESTS_Q023_PAPER_BELL_LINEAR_JOVERC_HPP_

#include "q023_paper_bell_linear.hpp"

#if defined(Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT)
#define Q023_JOVERC_INLINE inline
#else
#define Q023_JOVERC_INLINE KOKKOS_INLINE_FUNCTION
#endif

namespace q023_paper_bell_linear_joverc {

using q023_paper_bell_linear::Add;
using q023_paper_bell_linear::Basis;
using q023_paper_bell_linear::Cross;
using q023_paper_bell_linear::Dot;
using q023_paper_bell_linear::ModeBasis;
using q023_paper_bell_linear::ModeParameters;
using q023_paper_bell_linear::ModeSample;
using q023_paper_bell_linear::Scale;
using q023_paper_bell_linear::Vector3;

Q023_JOVERC_INLINE ModeSample UnstableEigenmodeAtPhase(
    const ModeParameters parameters, const double phase) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const double growth = std::sqrt(1.0 - parameters.epsilon*parameters.epsilon);
  const Vector3 magnetic = Add(
      Scale(basis.parallel, parameters.b_g),
      Add(Scale(basis.transverse_a, parameters.amplitude*cs),
          Scale(basis.transverse_b, parameters.amplitude*sn)));
  const double velocity_scale = parameters.amplitude/std::sqrt(parameters.rho);
  const Vector3 velocity = Add(
      Scale(basis.transverse_a,
            velocity_scale*(-parameters.epsilon*cs + growth*sn)),
      Scale(basis.transverse_b,
            velocity_scale*(-growth*cs - parameters.epsilon*sn)));
  return {magnetic, velocity};
}

Q023_JOVERC_INLINE Vector3 UnstableVectorPotentialAt(
    const ModeParameters parameters, const Vector3 position) {
  const Basis basis = ModeBasis(parameters.dimension);
  const double phase = parameters.k0*Dot(basis.parallel, position);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const Vector3 guide = Scale(Cross(basis.parallel, position),
                              0.5*parameters.b_g);
  const Vector3 perturbation = Add(
      Scale(basis.transverse_a, -parameters.amplitude*cs/parameters.k0),
      Scale(basis.transverse_b, -parameters.amplitude*sn/parameters.k0));
  return Add(guide, perturbation);
}

}  // namespace q023_paper_bell_linear_joverc

#undef Q023_JOVERC_INLINE

#endif  // PGEN_TESTS_Q023_PAPER_BELL_LINEAR_JOVERC_HPP_
