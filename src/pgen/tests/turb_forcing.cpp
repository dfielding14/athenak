//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_forcing.cpp
//! \brief Direct modal-process and conservative forcing-kick regression checks.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen/pgen.hpp"

namespace {
void Require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "### FATAL ERROR: turb_forcing: " << message << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
}  // namespace

void ProblemGenerator::TurbForcing(ParameterInput *pin, const bool restart) {
  Require(!restart, "direct checks require a fresh run");
  Turb(pin, restart);
  Mesh *pm = pmy_mesh_;
  auto *pack = pm->pmb_pack;
  auto *turb = pack->pturb;
  Require(turb != nullptr, "requires turb_driving");
  const auto metadata = turb->RestartMetadata();
  Require(metadata.tcorr > 1.0e-6, "OU check requires finite correlation time");
  turb->InitializeModes(nullptr, 0);
  Require(turb->n_turb_updates_yet == 1, "first OU innovation missing");
  auto previous_real = Kokkos::create_mirror(turb->mode_amp_real.d_view);
  Kokkos::deep_copy(previous_real, turb->mode_amp_real.d_view);
  auto previous_imag = Kokkos::create_mirror(turb->mode_amp_imag.d_view);
  Kokkos::deep_copy(previous_imag, turb->mode_amp_imag.d_view);
  for (int n = 0; n < turb->mode_count; ++n) {
    Real power = 0.0;
    for (int d = 0; d < 3; ++d) {
      power += SQR(previous_real(d,n)) + SQR(previous_imag(d,n));
    }
    Require(std::isfinite(power) && power > 0.0,
            "selected mode has zero or nonfinite normalization");
  }
  pm->time = 0.5 * metadata.dt_update;
  turb->InitializeModes(nullptr, 0);
  Require(turb->n_turb_updates_yet == 1, "OU state advanced between update boundaries");
  for (int n = 0; n < turb->mode_count; ++n) {
    for (int d = 0; d < 3; ++d) {
      Require(turb->mode_amp_real.h_view(d,n) == previous_real(d,n) &&
              turb->mode_amp_imag.h_view(d,n) == previous_imag(d,n),
              "held OU coefficient changed");
    }
  }
  pm->time = metadata.dt_update;
  turb->InitializeModes(nullptr, 0);
  Require(turb->n_turb_updates_yet == 2, "OU update boundary missed");
  const Real f = std::exp(-metadata.dt_update / metadata.tcorr);
  const Real g = std::sqrt(1.0 - f*f);
  for (int n = 0; n < turb->mode_count; ++n) {
    for (int d = 0; d < 3; ++d) {
      const Real er = f * previous_real(d,n) + g * turb->mode_noise_real.h_view(d,n);
      const Real ei = f * previous_imag(d,n) + g * turb->mode_noise_imag.h_view(d,n);
      Require(std::abs(turb->mode_amp_real.h_view(d,n) - er) < 1.0e-14 &&
              std::abs(turb->mode_amp_imag.h_view(d,n) - ei) < 1.0e-14,
              "OU recurrence does not use tcorr and dt_update");
    }
  }
  pm->time = 0.0;

  // Deliberately inconsistent primitives detect accidental use of stale velocities.
  // Isothermal inputs allocate one passive scalar at IEN to catch unguarded writes.
  auto force = Kokkos::create_mirror_view(turb->force);
  const auto &in = pm->mb_indcs;
  for (int m = 0; m < pack->nmb_thispack; ++m) {
    for (int k = in.ks; k <= in.ke; ++k) {
      for (int j = in.js; j <= in.je; ++j) {
        for (int i = in.is; i <= in.ie; ++i) {
          force(m,0,k,j,i) = 0.7;
          force(m,1,k,j,i) = -0.3;
          force(m,2,k,j,i) = 0.2;
        }
      }
    }
  }
  Kokkos::deep_copy(turb->force, force);
  DvceArray5D<Real> states[2], prims[2];
  bool ideal[2] = {false, false};
  int count = 0;
  if (pack->phydro != nullptr) {
    states[count] = pack->phydro->u0;
    prims[count] = pack->phydro->w0;
    ideal[count++] = pack->phydro->peos->eos_data.is_ideal;
  }
  if (pack->pmhd != nullptr) {
    states[count] = pack->pmhd->u0;
    prims[count] = pack->pmhd->w0;
    ideal[count++] = pack->pmhd->peos->eos_data.is_ideal;
  }
  auto saved0 = Kokkos::create_mirror(states[0]);
  Kokkos::deep_copy(saved0, states[0]);
  auto saved1 = Kokkos::create_mirror(states[count-1]);
  Kokkos::deep_copy(saved1, states[count-1]);
  for (int fluid = 0; fluid < count; ++fluid) {
    auto u = Kokkos::create_mirror_view_and_copy(HostMemSpace(), states[fluid]);
    for (int m = 0; m < pack->nmb_thispack; ++m) {
      for (int k = in.ks; k <= in.ke; ++k) {
        for (int j = in.js; j <= in.je; ++j) {
          for (int i = in.is; i <= in.ie; ++i) {
            const Real rho = 1.3 + 1.7*fluid + 0.01*(i-in.is);
            u(m,IDN,k,j,i) = rho;
            u(m,IM1,k,j,i) = rho*(0.4 + fluid);
            u(m,IM2,k,j,i) = rho*(-0.2 - 0.3*fluid);
            u(m,IM3,k,j,i) = rho*(0.6 - 0.7*fluid);
            if (ideal[fluid]) {
              u(m,IEN,k,j,i) = 11.0 + (SQR(u(m,IM1,k,j,i)) +
                  SQR(u(m,IM2,k,j,i)) + SQR(u(m,IM3,k,j,i))) / (2*rho);
            } else {
              u(m,IEN,k,j,i) = 17.0;
            }
          }
        }
      }
    }
    Kokkos::deep_copy(states[fluid], u);
    Kokkos::deep_copy(prims[fluid], -3.0);
  }
  turb->ApplyForcingWithStep(0.125);
  Real max_error = 0.0;
  for (int fluid = 0; fluid < count; ++fluid) {
    auto u = Kokkos::create_mirror_view_and_copy(HostMemSpace(), states[fluid]);
    for (int m = 0; m < pack->nmb_thispack; ++m) {
      for (int k = in.ks; k <= in.ke; ++k) {
        for (int j = in.js; j <= in.je; ++j) {
          for (int i = in.is; i <= in.ie; ++i) {
            const Real rho = 1.3 + 1.7*fluid + 0.01*(i-in.is);
            const Real expected[] = {rho*(0.4+fluid+0.7*0.125),
                rho*(-0.2-0.3*fluid-0.3*0.125), rho*(0.6-0.7*fluid+0.2*0.125)};
            for (int d = 0; d < 3; ++d) {
              max_error = fmax(max_error, std::abs(u(m,IM1+d,k,j,i)-expected[d]));
            }
            const Real thermal = ideal[fluid] ? u(m,IEN,k,j,i) -
                (SQR(u(m,IM1,k,j,i)) + SQR(u(m,IM2,k,j,i)) +
                 SQR(u(m,IM3,k,j,i))) / (2*rho) : u(m,IEN,k,j,i);
            max_error = fmax(max_error, std::abs(thermal-(ideal[fluid] ? 11.0 : 17.0)));
          }
        }
      }
    }
  }
  Require(max_error < 2.0e-14,
          "kick changed thermal energy/scalar or wrong fluid impulse");
  Kokkos::deep_copy(states[0], saved0);
  if (count == 2) Kokkos::deep_copy(states[1], saved1);
  std::cout << "turb_forcing: " << turb->mode_count << " nonzero modes; " << count
            << " fluid(s); max kick error=" << max_error << std::endl;
}
