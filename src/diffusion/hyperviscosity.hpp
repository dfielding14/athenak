#ifndef DIFFUSION_HYPERVISCOSITY_HPP_
#define DIFFUSION_HYPERVISCOSITY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hyperviscosity.hpp
//! \brief Fourth-derivative numerical viscosity for Newtonian Hydro and MHD.

#include <string>

#include "athena.hpp"
#include "diffusion/sts_types.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \class HyperViscosity
//! \brief Conservative velocity-biharmonic damping with explicit or STS integration.

class HyperViscosity {
 public:
  HyperViscosity(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~HyperViscosity() = default;

  Real dtnew;
  Real nu4;  // fourth-derivative kinematic coefficient, dimensions L^4/T
  parabolic::ParabolicIntegratorMode mode =
      parabolic::ParabolicIntegratorMode::explicit_mode;

  void AddHyperViscousFlux(const DvceArray5D<Real> &w, const EOS_Data &eos,
                           DvceFaceFld5D<Real> &f);
  void NewTimeStep();

 private:
  MeshBlockPack *pmy_pack;
  DvceArray5D<Real> lapv;  // cached cell-centered Laplacians of velocity components
};

#endif // DIFFUSION_HYPERVISCOSITY_HPP_
