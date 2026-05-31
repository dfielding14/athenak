//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#define Q029_HALL_BELL_LINEAR_HOST_CONTRACT 1
#include "../../src/pgen/tests/q029_hall_bell_linear.cpp"

#include <iomanip>
#include <iostream>

int main() {
  using q029_hall_bell_linear::ChiH;
  using q029_hall_bell_linear::CurrentDensity;
  using q029_hall_bell_linear::EigenmodeAtPhase;
  using q029_hall_bell_linear::IsPositivePreparedChiH;
  using q029_hall_bell_linear::ModeBasis;
  using q029_hall_bell_linear::ModeParameters;
  using q029_hall_bell_linear::NormalizedGrowthRate;
  using q029_hall_bell_linear::VectorPotentialAt;

  std::cout << std::setprecision(17);
  for (const double ppc : {4.0, 16.0}) {
    std::cout << "current_density " << ppc << ' '
              << CurrentDensity(ppc, 0.25, 2.0, 10.0) << '\n';
  }
  for (const double alpha_h : {0.125, 0.5, 2.0}) {
    std::cout << "chi_h " << alpha_h << ' '
              << ChiH(alpha_h, 4.0, 2.0, 1.0) << '\n';
  }
  for (const double chi_h : {
      -1.0, 0.0, 0.25, 0.2500000000005, 0.250000000002, 0.5, 1.0, 2.0}) {
    std::cout << "positive_prepared " << chi_h << ' '
              << IsPositivePreparedChiH(chi_h) << '\n';
  }
  for (const int dimension : {1, 2, 3}) {
    const ModeParameters parameters = {
      dimension, 0.4, 1.0e-6, 1.0, 1.0, 1.0, 1.0, 2.0*std::acos(-1.0)
    };
    const auto basis = ModeBasis(dimension);
    const auto sample = EigenmodeAtPhase(parameters, 0.37);
    const auto vector_potential = VectorPotentialAt(parameters, {0.23, 0.41, 0.67});
    std::cout << "q023_carrier " << dimension << ' '
              << NormalizedGrowthRate(parameters.epsilon) << ' '
              << basis.parallel.x1 << ' ' << basis.parallel.x2 << ' '
              << basis.parallel.x3 << ' ' << sample.magnetic.x1 << ' '
              << sample.magnetic.x2 << ' ' << sample.magnetic.x3 << ' '
              << sample.velocity.x1 << ' ' << sample.velocity.x2 << ' '
              << sample.velocity.x3 << ' ' << vector_potential.x1 << ' '
              << vector_potential.x2 << ' ' << vector_potential.x3 << '\n';
  }
}
