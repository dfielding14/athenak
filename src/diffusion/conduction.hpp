#ifndef DIFFUSION_CONDUCTION_HPP_
#define DIFFUSION_CONDUCTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.hpp
//! \brief Contains data and functions that implement various formulations for conduction.
//  Currently only isotropic conduction implemented

#include <string>

#include "athena.hpp"
#include "diffusion/sts_types.hpp"
#include "parameter_input.hpp"

//----------------------------------------------------------------------------------------
//! \class Conduction
//! \brief data and functions that implement thermal conduction in Hydro and MHD

class Conduction {
 public:
  Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~Conduction();

  // data
  Real dtnew;
  std::string iso_cond_type = "none";       // "constant", "spitzer", or "spitzer_limited"
  std::string cgl_heat_flux_type = "none";  // "landau_fluid" for CGL LF heat flux
  Real kappa_iso = 0.0;                      // isotropic thermal conductivity
  Real kappa_iso_limit = 0.0;                // limit to isotropic thermal conductivity
  Real lf_k_parallel = 0.0;                  // CGL LF |k_parallel| from the closure
  bool lf_coeff_local = true;                // true: local c_parallel; false: background
  Real lf_c_parallel0 = 0.0;                 // background c_parallel for LF coefficients
  parabolic::ParabolicIntegratorMode mode = parabolic::ParabolicIntegratorMode::explicit_mode;

  // functions
  bool IsCGLLandauFluidHeatFlux() const { return cgl_heat_flux_type == "landau_fluid"; }
  void AddHeatFluxes(const DvceArray5D<Real> &w, const EOS_Data &eos,
                     DvceFaceFld5D<Real> &f);
  void AddCGLLandauFluidHeatFluxes(const DvceArray5D<Real> &w,
                                   const DvceArray5D<Real> &bcc,
                                   const EOS_Data &eos,
                                   DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                       DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);

 private:
  MeshBlockPack* pmy_pack;
};
#endif // DIFFUSION_CONDUCTION_HPP_
