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
#include "diffusion/sts_rkl2.hpp"
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
  std::uint64_t qfaces = 0;
  std::uint64_t qpar_cap = 0;
  std::uint64_t qpar_cap10 = 0;
  std::uint64_t qperp_cap = 0;
  std::uint64_t qperp_cap10 = 0;
  Real qpar_work = 0.0;
  Real qperp_work = 0.0;
  Real pressure_work = 0.0;
  Real anisotropic_pressure_work = 0.0;
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
  void AdvanceHeatFluxWorkDiagnostics(Real dt_sweep,
                                      const parabolic::RKL2Coefficients &coeffs,
                                      int stage, int nstages);
  void AdvancePressureWorkDiagnostics(Real beta_dt, Real gam0, Real gam1, int stage,
                                      Real pressure_power, Real anisotropic_power);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos);
  void RecordAdmissibility(const DvceArray5D<Real> &u, const DvceArray5D<Real> &w,
                           const DvceArray5D<Real> &bcc, const EOS_Data &eos,
                           int dfloor_delta, int pfloor_delta);

 private:
  void AccumulateHeatFluxDiagnostics(const array_sum::GlobalSum &stats);

  MeshBlockPack *pmy_pack;
  DvceArray4D<Real> tpar_, tperp_, bmag_;
  Real stage_qpar_power_ = 0.0;
  Real stage_qperp_power_ = 0.0;
  Real sweep_qpar_work_ = 0.0;
  Real sweep_qperp_work_ = 0.0;
  Real sweep_qpar_work1_ = 0.0;
  Real sweep_qperp_work1_ = 0.0;
  Real sweep_qpar_work2_ = 0.0;
  Real sweep_qperp_work2_ = 0.0;
  Real sweep_qpar_rhs_ = 0.0;
  Real sweep_qperp_rhs_ = 0.0;
  Real pressure_work_cycle_start_ = 0.0;
  Real anisotropic_pressure_work_cycle_start_ = 0.0;
};

#endif // DIFFUSION_CGL_LANDAU_FLUID_HPP_
