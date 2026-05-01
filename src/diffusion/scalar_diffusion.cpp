//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_diffusion.cpp
//! \brief Implements functions for ScalarDiffusion class. This implements isotropic
//! diffusion for passive scalars, in which the diffusive flux is proportional to the
//! negative local scalar gradient: F = -kappa * rho * grad(s)

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "scalar_diffusion.hpp"

namespace {

parabolic::ParabolicIntegratorMode ParseScalarDiffusivityIntegrator(std::string block,
    ParameterInput *pin) {
  std::string integrator = pin->GetOrAddString(block, "scalar_diffusivity_integrator",
                                               "explicit");
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  } else if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<" << block << ">/scalar_diffusivity_integrator = '" << integrator
            << "' must be 'explicit' or 'sts'" << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \brief ScalarDiffusion constructor

ScalarDiffusion::ScalarDiffusion(std::string block, MeshBlockPack *pp,
                                 ParameterInput *pin) :
  pmy_pack(pp) {
  // Read scalar diffusivity coefficient
  kappa = pin->GetReal(block, "scalar_diffusivity");
  mode = ParseScalarDiffusivityIntegrator(block, pin);

  // Compute diffusion timestep constraint
  dtnew = std::numeric_limits<float>::max();
  if (kappa <= 0.0) {
    // Treat non-positive diffusivity as "disabled" (but allow the parameter to exist).
    return;
  }
  auto size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pp->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pp->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  for (int m = 0; m < (pp->nmb_thispack); ++m) {
    dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx1)/kappa);
    if (pp->pmesh->multi_d) {
      dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx2)/kappa);
    }
    if (pp->pmesh->three_d) {
      dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx3)/kappa);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief ScalarDiffusion destructor

ScalarDiffusion::~ScalarDiffusion() {
}

//----------------------------------------------------------------------------------------
//! \fn void ScalarDiffusion::AddScalarDiffusionFlux()
//! \brief Adds scalar diffusion fluxes to face-centered fluxes of conserved variables.
//! The diffusive flux for each scalar n is: F_n = -kappa * rho_face * grad(s_n)
//! where s_n is the primitive scalar (concentration) at index nhydro+n in w0.

void ScalarDiffusion::AddScalarDiffusionFlux(const DvceArray5D<Real> &w0,
    const int nhydro, const int nscalars, DvceFaceFld5D<Real> &flx) {
  // Return early if no scalars
  if (nscalars == 0 || kappa <= 0.0) return;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  Real kappa_ = kappa;
  int nhydro_ = nhydro;
  int nscalars_ = nscalars;

  //--------------------------------------------------------------------------------------
  // fluxes in x1-direction

  auto &flx1 = flx.x1f;

  par_for("scalar_diff1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Face-centered density (arithmetic average)
    Real rho_face = 0.5 * (w0(m,IDN,k,j,i) + w0(m,IDN,k,j,i-1));

    // Add diffusion flux for each scalar
    for (int n = 0; n < nscalars_; ++n) {
      int ns = nhydro_ + n;  // index in w0/flx arrays

      // Scalar gradient at face: ds/dx
      Real ds_dx = (w0(m,ns,k,j,i) - w0(m,ns,k,j,i-1)) / size.d_view(m).dx1;

      // Diffusive flux: F = -kappa * rho * grad(s)
      flx1(m,ns,k,j,i) -= kappa_ * rho_face * ds_dx;
    }
  });
  if (pmy_pack->pmesh->one_d) {return;}

  //--------------------------------------------------------------------------------------
  // fluxes in x2-direction

  auto &flx2 = flx.x2f;

  par_for("scalar_diff2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Face-centered density
    Real rho_face = 0.5 * (w0(m,IDN,k,j,i) + w0(m,IDN,k,j-1,i));

    // Add diffusion flux for each scalar
    for (int n = 0; n < nscalars_; ++n) {
      int ns = nhydro_ + n;

      // Scalar gradient at face: ds/dy
      Real ds_dy = (w0(m,ns,k,j,i) - w0(m,ns,k,j-1,i)) / size.d_view(m).dx2;

      // Diffusive flux: F = -kappa * rho * grad(s)
      flx2(m,ns,k,j,i) -= kappa_ * rho_face * ds_dy;
    }
  });
  if (pmy_pack->pmesh->two_d) {return;}

  //--------------------------------------------------------------------------------------
  // fluxes in x3-direction

  auto &flx3 = flx.x3f;

  par_for("scalar_diff3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Face-centered density
    Real rho_face = 0.5 * (w0(m,IDN,k,j,i) + w0(m,IDN,k-1,j,i));

    // Add diffusion flux for each scalar
    for (int n = 0; n < nscalars_; ++n) {
      int ns = nhydro_ + n;

      // Scalar gradient at face: ds/dz
      Real ds_dz = (w0(m,ns,k,j,i) - w0(m,ns,k-1,j,i)) / size.d_view(m).dx3;

      // Diffusive flux: F = -kappa * rho * grad(s)
      flx3(m,ns,k,j,i) -= kappa_ * rho_face * ds_dz;
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ScalarDiffusion::NewTimeStep()
//! \brief Compute new time step for scalar diffusion.
//! The diffusion timestep constraint is: dt < fac * dx^2 / kappa
//! where fac = 1/(2*ndim) for explicit stability.

void ScalarDiffusion::NewTimeStep(const DvceArray5D<Real> &w0,
    const int nhydro, const int nscalars) {
  // Return early if no scalars or diffusion disabled
  (void)w0;
  (void)nhydro;
  if (nscalars == 0 || kappa <= 0.0) {
    dtnew = std::numeric_limits<float>::max();
    return;
  }

  auto &size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }

  // Find minimum dx^2/kappa across all meshblocks
  dtnew = std::numeric_limits<float>::max();
  for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
    dtnew = std::min(dtnew, fac * SQR(size.h_view(m).dx1) / kappa);
    if (pmy_pack->pmesh->multi_d) {
      dtnew = std::min(dtnew, fac * SQR(size.h_view(m).dx2) / kappa);
    }
    if (pmy_pack->pmesh->three_d) {
      dtnew = std::min(dtnew, fac * SQR(size.h_view(m).dx3) / kappa);
    }
  }

  return;
}
