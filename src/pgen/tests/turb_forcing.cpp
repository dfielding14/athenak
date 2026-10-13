//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_forcing.cpp
//! \brief Uniform periodic states and direct conservative forcing-kick checks.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/turb_driver.hpp"

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  auto *pack = pmy_mesh_->pmb_pack;
  const auto &in = pmy_mesh_->mb_indcs;
  auto *turb = pack->pturb;
  const bool check_kick = pin->GetOrAddBoolean("problem", "check_kick", false);
  DvceArray5D<Real> states[2], prims[2];
  bool ideal[2], magnetic[2];
  Real gamma[2];
  int count = 0;
  if (pack->phydro != nullptr) {
    states[count] = pack->phydro->u0;
    prims[count] = pack->phydro->w0;
    ideal[count] = pack->phydro->peos->eos_data.is_ideal;
    magnetic[count] = false;
    gamma[count++] = pack->phydro->peos->eos_data.gamma;
  }
  if (pack->pmhd != nullptr) {
    states[count] = pack->pmhd->u0;
    prims[count] = pack->pmhd->w0;
    ideal[count] = pack->pmhd->peos->eos_data.is_ideal;
    magnetic[count] = true;
    gamma[count++] = pack->pmhd->peos->eos_data.gamma;
    Kokkos::deep_copy(pack->pmhd->b0.x1f, 0.0);
    Kokkos::deep_copy(pack->pmhd->b0.x2f, 0.0);
    Kokkos::deep_copy(pack->pmhd->b0.x3f, 1.0);
  }
  if (count == 0 ||
      (check_kick && turb == nullptr)) {
    std::cerr << "### FATAL ERROR: turb_forcing requires a fluid and forcing for checks"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for (int fluid = 0; fluid < count; ++fluid) {
    auto u = Kokkos::create_mirror_view(states[fluid]);
    Kokkos::deep_copy(u, 0.0);
    if (check_kick && u.extent(1) <= IEN) {
      std::cerr << "### FATAL ERROR: isothermal kick checks require nscalars = 1\n";
      std::exit(EXIT_FAILURE);
    }
    for (int m = 0; m < pack->nmb_thispack; ++m) {
      for (int k = in.ks; k <= in.ke; ++k) {
        for (int j = in.js; j <= in.je; ++j) {
          for (int i = in.is; i <= in.ie; ++i) {
            const Real rho = check_kick ? 1.3 + fluid + 0.01*(i-in.is) : 1.0;
            u(m,IDN,k,j,i) = rho;
            Real kinetic = 0.0;
            for (int d = 0; d < 3; ++d) {
              const Real v = check_kick ? 0.4 + fluid - 0.3*d : 0.0;
              u(m,IM1+d,k,j,i) = rho*v;
              kinetic += 0.5*rho*v*v;
            }
            if (ideal[fluid]) {
              u(m,IEN,k,j,i) = (check_kick ? 11.0 : 1.0/(gamma[fluid]-1.0)) +
                              kinetic + (magnetic[fluid] ? 0.5 : 0.0);
            } else if (check_kick) {
              // The passive scalar at IEN detects energy writes to isothermal fluids.
              u(m,IEN,k,j,i) = 17.0;
            }
          }
        }
      }
    }
    Kokkos::deep_copy(states[fluid], u);
    if (check_kick) Kokkos::deep_copy(prims[fluid], -3.0);
  }
  if (!check_kick) return;

  Kokkos::deep_copy(turb->force, 0.7);
  turb->ApplyImpulse(0.125);
  Real max_error = 0.0;
  for (int fluid = 0; fluid < count; ++fluid) {
    auto u = Kokkos::create_mirror_view_and_copy(HostMemSpace(), states[fluid]);
    for (int m = 0; m < pack->nmb_thispack; ++m) {
      for (int k = in.ks; k <= in.ke; ++k) {
        for (int j = in.js; j <= in.je; ++j) {
          for (int i = in.is; i <= in.ie; ++i) {
            const Real rho = 1.3 + fluid + 0.01*(i-in.is);
            Real kinetic = 0.0;
            if (!std::isfinite(u(m,IDN,k,j,i))) {
              std::cerr << "### FATAL ERROR: nonfinite forcing density\n";
              std::exit(EXIT_FAILURE);
            }
            max_error = fmax(max_error, std::abs(u(m,IDN,k,j,i)-rho));
            for (int d = 0; d < 3; ++d) {
              const Real mom = u(m,IM1+d,k,j,i);
              if (!std::isfinite(mom)) {
                std::cerr << "### FATAL ERROR: nonfinite forcing momentum\n";
                std::exit(EXIT_FAILURE);
              }
              const Real expected = rho*(0.4+fluid-0.3*d+0.7*0.125);
              max_error = fmax(max_error, std::abs(mom-expected));
              kinetic += 0.5*mom*mom/rho;
            }
            const Real thermal = ideal[fluid] ? u(m,IEN,k,j,i)-kinetic-
                                 (magnetic[fluid] ? 0.5 : 0.0) : u(m,IEN,k,j,i);
            if (!std::isfinite(thermal)) {
              std::cerr << "### FATAL ERROR: nonfinite forcing energy/scalar\n";
              std::exit(EXIT_FAILURE);
            }
            max_error = fmax(max_error, std::abs(thermal-(ideal[fluid] ? 11.0 : 17.0)));
          }
        }
      }
    }
  }
  if (max_error > 2.0e-14 || !std::isfinite(max_error)) {
    std::cerr << "### FATAL ERROR: forcing kick changed thermal energy/scalar or impulse;"
              << " max error=" << max_error << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "turb_forcing: " << count << " fluid(s); max kick error=" << max_error
            << std::endl;
}
