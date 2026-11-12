#ifndef DIFFUSION_VISCOSITY_HPP_
#define DIFFUSION_VISCOSITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file viscosity.hpp
//  \brief Contains data and functions that implement various formulations for
//  viscosity, including Smagorinsky eddy viscosity for 2D Newtonian Hydro.

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \class Viscosity
//  \brief data and functions that implement viscosity in Hydro and MHD

class Viscosity {
 public:
  Viscosity(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~Viscosity();

  // data
  Real dtnew;
  Real nu_iso;     // coefficient of isotropic kinematic shear viscosity

  // function to add viscous fluxes to Hydro and/or MHD fluxes
  void IsotropicViscousFlux(const DvceArray5D<Real> &w, const Real nu,
                            const EOS_Data &eos, DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w);

 private:
  MeshBlockPack* pmy_pack;
  Real smagorinsky_cs = 0.0;
};

#endif // DIFFUSION_VISCOSITY_HPP_
