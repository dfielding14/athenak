#ifndef DIFFUSION_CGL_LANDAU_FLUID_HPP_
#define DIFFUSION_CGL_LANDAU_FLUID_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.hpp
//! \brief CGL Landau-fluid heat-flux closure advanced by MHD-owned STS.

#include "athena.hpp"
#include "diffusion/sts_types.hpp"

class ParameterInput;

class CGLLandauFluid {
 public:
  CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin);

  Real dtnew;
  Real lf_k_parallel;
  bool lf_coeff_local;
  Real lf_c_parallel0;
  parabolic::ParabolicIntegratorMode mode;

  void AddHeatFluxes(const DvceArray5D<Real> &w, const DvceArray5D<Real> &bcc,
                     const EOS_Data &eos, DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos);

 private:
  MeshBlockPack *pmy_pack;
  DvceArray4D<Real> tpar_, tperp_, bmag_;
};

#endif // DIFFUSION_CGL_LANDAU_FLUID_HPP_
