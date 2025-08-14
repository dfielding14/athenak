#ifndef DIFFUSION_RESISTIVITY_HPP_
#define DIFFUSION_RESISTIVITY_HPP_

#include "athena.hpp"
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

  // functions to add resistive E-Field and energy flux
  void OhmicEField(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void OhmicEnergyFlux(const DvceFaceFld4D<Real> &b, DvceFaceFld5D<Real> &flx);

 private:
  MeshBlockPack* pmy_pack;
};

#endif // DIFFUSION_RESISTIVITY_HPP_
