//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_transform_check.cpp
//! \brief Unit-style checks for CGL conserved anisotropy and magnetic-moment transforms.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {

void RequireClose(const char *label, const Real got, const Real expected, const Real tol) {
  const Real err = std::abs(got - expected);
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (err > tol*scale) {
    std::cout << "CGL transform check failed for " << label
              << ": got " << got << ", expected " << expected
              << ", relative error " << err/scale << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void) pin;
  if (restart) return;

  const Real rho = 1.7;
  const Real vx = 0.11;
  const Real vy = -0.07;
  const Real vz = 0.19;
  const Real mx = rho*vx;
  const Real my = rho*vy;
  const Real mz = rho*vz;
  const Real p_parallel = 0.83;
  const Real p_perp = 1.24;
  const Real bx = 0.91;
  const Real by = -0.28;
  const Real bz = 0.37;
  const Real bmag = std::sqrt(SQR(bx) + SQR(by) + SQR(bz));
  const Real bfloor = 1.0e-20;
  const Real pfloor = 1.0e-20;
  const Real kinetic = 0.5*rho*(SQR(vx) + SQR(vy) + SQR(vz));
  const Real magnetic = 0.5*SQR(bmag);
  const Real total_energy = 0.5*p_parallel + p_perp + kinetic + magnetic;

  const Real anisotropy = CGLConservedAnisotropy(rho, p_parallel, p_perp, bmag);
  Real recovered_parallel = 0.0, recovered_perp = 0.0;
  CGLRecoverPressuresFromTotalEnergyAndAnisotropy(rho, mx, my, mz, total_energy,
                                                  anisotropy, bx, by, bz, bfloor,
                                                  recovered_parallel, recovered_perp);
  RequireClose("p_parallel", recovered_parallel, p_parallel, 1.0e-13);
  RequireClose("p_perp", recovered_perp, p_perp, 1.0e-13);

  const Real magnetic_moment =
      CGLConservedAnisotropyToMagneticMoment(rho, mx, my, mz, total_energy,
                                             anisotropy, bx, by, bz, bfloor, pfloor);
  Real moment_energy = total_energy;
  Real roundtrip_anisotropy = 0.0;
  CGLMagneticMomentToConservedAnisotropy(rho, mx, my, mz, bx, by, bz, bfloor, pfloor,
                                         magnetic_moment, moment_energy,
                                         roundtrip_anisotropy);
  RequireClose("A round trip", roundtrip_anisotropy, anisotropy, 1.0e-13);
  RequireClose("energy round trip", moment_energy, total_energy, 1.0e-13);

  std::cout << "CGL transform round trips passed; magnetic_moment="
            << magnetic_moment << std::endl;
}
