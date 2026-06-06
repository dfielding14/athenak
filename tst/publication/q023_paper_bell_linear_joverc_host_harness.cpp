// Host-only contract harness for the corrected Q023 J/c Bell predecessor.

#include <iomanip>
#include <iostream>

#define Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT 1
#include "../../src/pgen/tests/q023_paper_bell_linear_joverc.cpp"

int main() {
  using namespace q023_paper_bell_linear_joverc;
  std::cout << std::setprecision(17);
  for (const double ppc : {0.5, 1.0, 4.0, 4.5}) {
    std::cout << "ppc " << ppc << " " << HasPositiveIntegralPPC(ppc) << "\n";
  }
  const double charge = 4.0*M_PI*1.0e-6;
  const double mass = 2.0;
  const double speed = 2.5;
  for (const double volume : {1.0/128.0, 1.0/819.2, 4.588798980430252e-5}) {
    const double qscale = RequiredDepositQScale(1.0, charge, speed, volume,
                                                1.0, 2.0*M_PI);
    std::cout << "closure " << volume << " " << qscale << " "
              << DepositedJOverC(1.0, qscale, charge, speed, volume) << " "
              << HasRequiredDepositedJOverC(1.0, qscale, charge, speed, volume,
                                            1.0, 2.0*M_PI)
              << "\n";
  }
  std::cout << "species "
            << HasRequiredSpeciesChargeOverMass(mass, charge, 2.0*M_PI*1.0e-6)
            << " "
            << HasRequiredSpeciesChargeOverMass(mass, charge, 1.0e-6)
            << "\n";
  std::cout << "velocity_match "
            << HasMatchingVelocityVectors(2.5, 0.0, 0.0, 2.5, 0.0, 0.0) << " "
            << HasMatchingVelocityVectors(2.5, 0.0, 0.0, 0.0, 0.0, 0.0) << " "
            << HasMatchingVelocityVectors(2.5, 0.0, 0.0, 2.0, 0.0, 0.0) << "\n";
}
