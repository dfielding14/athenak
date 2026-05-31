//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#define Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#include "../../src/pgen/tests/q023_paper_bell_linear.hpp"

#include <iomanip>
#include <iostream>

namespace {

using q023_paper_bell_linear::Add;
using q023_paper_bell_linear::ModeParameters;
using q023_paper_bell_linear::Scale;
using q023_paper_bell_linear::Vector3;
using q023_paper_bell_linear::VectorPotentialAt;

Vector3 CurlAt(const ModeParameters parameters, const Vector3 position,
               const double spacing) {
  const Vector3 dx1 = {spacing, 0.0, 0.0};
  const Vector3 dx2 = {0.0, spacing, 0.0};
  const Vector3 dx3 = {0.0, 0.0, spacing};
  const Vector3 a_x1p = VectorPotentialAt(parameters, Add(position, dx1));
  const Vector3 a_x1m = VectorPotentialAt(parameters, Add(position, Scale(dx1, -1.0)));
  const Vector3 a_x2p = VectorPotentialAt(parameters, Add(position, dx2));
  const Vector3 a_x2m = VectorPotentialAt(parameters, Add(position, Scale(dx2, -1.0)));
  const Vector3 a_x3p = VectorPotentialAt(parameters, Add(position, dx3));
  const Vector3 a_x3m = VectorPotentialAt(parameters, Add(position, Scale(dx3, -1.0)));
  const double inv_width = 0.5/spacing;
  return {
    ((a_x2p.x3 - a_x2m.x3) - (a_x3p.x2 - a_x3m.x2))*inv_width,
    ((a_x3p.x1 - a_x3m.x1) - (a_x1p.x3 - a_x1m.x3))*inv_width,
    ((a_x1p.x2 - a_x1m.x2) - (a_x2p.x1 - a_x2m.x1))*inv_width
  };
}

}  // namespace

int main() {
  using q023_paper_bell_linear::Dot;
  using q023_paper_bell_linear::EigenmodeAtPhase;
  using q023_paper_bell_linear::ModeBasis;

  std::cout << std::setprecision(17);
  for (const int dimension : {1, 2, 3}) {
    for (const double epsilon : {0.1, 0.4, 0.8}) {
      const ModeParameters parameters = {
        dimension, epsilon, 1.0e-6, 1.0, 1.0, 1.0, 1.0, 2.0*std::acos(-1.0)
      };
      for (const double phase : {0.0, 0.37, 1.7}) {
        const auto sample = EigenmodeAtPhase(parameters, phase);
        std::cout << "sample " << dimension << ' ' << epsilon << ' ' << phase << ' '
                  << sample.magnetic.x1 << ' ' << sample.magnetic.x2 << ' '
                  << sample.magnetic.x3 << ' ' << sample.velocity.x1 << ' '
                  << sample.velocity.x2 << ' ' << sample.velocity.x3 << '\n';
      }
      const Vector3 position = {0.23, 0.41, 0.67};
      const auto basis = ModeBasis(dimension);
      const auto expected = EigenmodeAtPhase(
          parameters, parameters.k0*Dot(basis.parallel, position)).magnetic;
      const auto measured = CurlAt(parameters, position, 1.0e-5);
      std::cout << "curl " << dimension << ' ' << epsilon << ' '
                << measured.x1 - expected.x1 << ' '
                << measured.x2 - expected.x2 << ' '
                << measured.x3 - expected.x3 << '\n';
    }
  }
}
