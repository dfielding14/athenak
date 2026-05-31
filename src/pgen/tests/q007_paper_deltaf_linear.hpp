//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q007_paper_deltaf_linear.hpp
//! \brief Shared Q-007 Sun and Bai true-delta-f linear preparation math.

#ifndef PGEN_TESTS_Q007_PAPER_DELTAF_LINEAR_HPP_
#define PGEN_TESTS_Q007_PAPER_DELTAF_LINEAR_HPP_

#include <cmath>

namespace q007_paper_deltaf_linear {

constexpr double kPi = 3.141592653589793238462643383279502884;

inline double KappaNormalization(const double number_density, const double p0,
                                 const double kappa) {
  return number_density*std::tgamma(kappa + 1.0)/
      (std::pow(kPi*kappa*p0*p0, 1.5)*std::tgamma(kappa - 0.5));
}

inline double IsotropicKappaDistribution(const double number_density,
                                         const double p, const double p0,
                                         const double kappa) {
  const double shape = 1.0 + p*p/(kappa*p0*p0);
  return KappaNormalization(number_density, p0, kappa)*
      std::pow(shape, -(kappa + 1.0));
}

inline double AnisotropicKappaDistribution(const double number_density,
                                           const double px, const double py,
                                           const double pz, const double p0,
                                           const double kappa,
                                           const double xi) {
  const double xi2 = xi*xi;
  const double shape = 1.0 + (px*px + xi2*(py*py + pz*pz))/(kappa*p0*p0);
  return xi2*KappaNormalization(number_density, p0, kappa)*
      std::pow(shape, -(kappa + 1.0));
}

inline double AthenaKTransverseAnisotropyScale(const double xi) {
  return 1.0/xi;
}

inline double ResonantWavenumber(const double mass, const double omega0,
                                 const double p0) {
  return mass*omega0/p0;
}

inline double ResonantWavelength(const double mass, const double omega0,
                                 const double p0) {
  return 2.0*kPi/ResonantWavenumber(mass, omega0, p0);
}

inline double GyroresonanceQ2(const double k, const double mass,
                             const double omega0, const double p0,
                             const double kappa) {
  const double ratio = mass*omega0/(std::abs(k)*p0);
  const double prefactor = std::sqrt(kPi)/std::pow(kappa, 1.5)*
      std::tgamma(kappa + 1.0)/std::tgamma(kappa - 0.5);
  return prefactor*ratio*
      std::pow(1.0 + ratio*ratio/kappa, -kappa);
}

inline double CRSILowDensityGrowthRate(const double k, const double mass,
                                       const double omega0, const double p0,
                                       const double kappa,
                                       const double mncr_over_rho0,
                                       const double vd_over_ua) {
  const double direction = (k >= 0.0) ? 1.0 : -1.0;
  return -0.5*mncr_over_rho0*omega0*
      (1.0 - vd_over_ua*direction)*
      GyroresonanceQ2(k, mass, omega0, p0, kappa);
}

inline double CRPAILowDensityGrowthRate(const double k, const double mass,
                                        const double omega0, const double p0,
                                        const double kappa,
                                        const double mncr_over_rho0,
                                        const double ua, const double xi,
                                        const int signed_branch) {
  const double xi2 = xi*xi;
  const double branch = static_cast<double>(signed_branch);
  return -0.5*mncr_over_rho0*omega0*
      (1.0 + branch*(1.0 - xi2)/xi2*omega0/(std::abs(k)*ua))*
      GyroresonanceQ2(k, mass, omega0, p0, kappa)/xi2;
}

}  // namespace q007_paper_deltaf_linear

#endif  // PGEN_TESTS_Q007_PAPER_DELTAF_LINEAR_HPP_
