#ifndef DIFFUSION_RESISTIVITY_HPP_
#define DIFFUSION_RESISTIVITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.hpp
//  \brief Contains data and functions that implement various non-ideal MHD (resistive)
//  processes, such as Ohmic diffusion. TODO(@user): add ambipolar diffusion, Hall effect

#include "athena.hpp"
#include "diffusion/current_limited_resistivity.hpp"
#include "diffusion/sts_types.hpp"
#include "parameter_input.hpp"
#include "mesh/meshblock.hpp"

//----------------------------------------------------------------------------------------
//! \class Resistivity
//  \brief data and functions that implement various resistive physics

class Resistivity {
 public:
  Resistivity(MeshBlockPack *pp, ParameterInput *pin);
  ~Resistivity();

  // data
  Real dtnew;
  Real eta_ohm;
  bool current_limited = false;
  bool b_rec_jump = false;
  Real b_rec_radius = 0.0, b_rec_floor = 0.0;
  ::current_limited::Parameters current_limited_params{};
  parabolic::DiffusionSelection mode = parabolic::DiffusionSelection::explicit_only;

  // functions to add resistive E-Field and energy flux
  void OhmicEField(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void OhmicEnergyFlux(const DvceFaceFld4D<Real> &b, DvceFaceFld5D<Real> &flx);
  void NewTimeStep();
  void PrepareBRec();

 private:
  void CurrentLimitedEField(const DvceFaceFld4D<Real> &b, DvceEdgeFld4D<Real> &efld);
  void CurrentLimitedEnergyFlux(const DvceFaceFld4D<Real> &b, DvceFaceFld5D<Real> &flx);
  MeshBlockPack* pmy_pack;
};

#endif // DIFFUSION_RESISTIVITY_HPP_
