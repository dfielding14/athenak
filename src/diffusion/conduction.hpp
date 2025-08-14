#ifndef DIFFUSION_CONDUCTION_HPP_
#define DIFFUSION_CONDUCTION_HPP_

#include <string>

#include "athena.hpp"
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
  Real kappa;         // thermal conductivity
  bool tdep_kappa;    // temperature-dependent conductivity
  Real kappa_ceiling; // ceiling of thermal conductivity
  bool sat_hflux;     // saturtion of heat flux

  // function to add heat fluxes to Hydro and/or MHD fluxes
  void AddHeatFlux(const DvceArray5D<Real> &w, const EOS_Data &eos,
                   DvceFaceFld5D<Real> &f);
  void IsotropicHeatFlux(const DvceArray5D<Real> &w, const EOS_Data &eos,
                         DvceFaceFld5D<Real> &f);
  void TempDependentHeatFlux(const DvceArray5D<Real> &w, const EOS_Data &eos,
                             DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);

 private:
  MeshBlockPack* pmy_pack;
};
#endif // DIFFUSION_CONDUCTION_HPP_
