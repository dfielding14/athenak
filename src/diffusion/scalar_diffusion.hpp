#ifndef DIFFUSION_SCALAR_DIFFUSION_HPP_
#define DIFFUSION_SCALAR_DIFFUSION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_diffusion.hpp
//! \brief Contains data and functions that implement scalar diffusion for passive scalars.
//! The diffusive flux is F = -kappa * rho * grad(s), where s is the scalar concentration.

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"

//----------------------------------------------------------------------------------------
//! \class ScalarDiffusion
//! \brief data and functions that implement scalar diffusion in Hydro

class ScalarDiffusion {
 public:
  ScalarDiffusion(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~ScalarDiffusion();

  // data
  Real dtnew;         // diffusion timestep constraint
  Real kappa;         // scalar diffusivity coefficient D

  // function to add scalar diffusion fluxes
  void AddScalarDiffusionFlux(const DvceArray5D<Real> &w, const int nhydro,
                              const int nscalars, DvceFaceFld5D<Real> &f);

  // function to compute diffusion timestep
  void NewTimeStep(const DvceArray5D<Real> &w, const int nhydro, const int nscalars);

 private:
  MeshBlockPack* pmy_pack;
};
#endif // DIFFUSION_SCALAR_DIFFUSION_HPP_
