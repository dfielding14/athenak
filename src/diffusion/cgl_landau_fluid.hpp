#ifndef DIFFUSION_CGL_LANDAU_FLUID_HPP_
#define DIFFUSION_CGL_LANDAU_FLUID_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_landau_fluid.hpp
//! \brief CGL Landau-fluid heat-flux closure advanced by MHD-owned split sweeps.

#include <cstdint>

#include "athena.hpp"
#include "diffusion/sts_types.hpp"

class ParameterInput;

struct CGLLFDiagnostics {
  std::uint64_t nstage = 0;
  std::uint64_t dfloor = 0;
  std::uint64_t pfloor = 0;
  std::uint64_t nonfinite = 0;
  std::uint64_t nonpositive = 0;
  std::uint64_t mirror = 0;
  std::uint64_t firehose = 0;
  std::uint64_t hard_bound = 0;
};

class CGLLandauFluid {
 public:
  CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin);

  Real dtnew;
  Real lf_k_parallel;
  bool lf_coeff_local;
  Real lf_c_parallel0;
  bool strict_admissibility;
  parabolic::ParabolicIntegratorMode mode;
  CGLLFDiagnostics diagnostics;

  void AddHeatFluxes(const DvceArray5D<Real> &w, const DvceArray5D<Real> &bcc,
                     const EOS_Data &eos, DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos);
  void RecordAdmissibility(const DvceArray5D<Real> &u, const DvceArray5D<Real> &w,
                           const DvceArray5D<Real> &bcc, const EOS_Data &eos,
                           int dfloor_delta, int pfloor_delta);

 private:
  MeshBlockPack *pmy_pack;
  DvceArray4D<Real> tpar_, tperp_, bmag_;
};

#endif // DIFFUSION_CGL_LANDAU_FLUID_HPP_
