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
  std::uint64_t hardwall_projection = 0;
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

enum class CGLLFDiagnosticsMode {
  full,
  none
};

enum class CGLLFArithmeticMode {
  safe,
  fast
};

enum class CGLLFSTSFluxMode {
  weighted,
  physical
};

enum class CGLLFProfileBucket {
  heat_flux_total = 0,
  heat_flux_precompute,
  heat_flux_flux1,
  heat_flux_flux1_gradients,
  heat_flux_flux1_face_state,
  heat_flux_flux1_closure,
  heat_flux_flux2,
  heat_flux_flux2_gradients,
  heat_flux_flux2_face_state,
  heat_flux_flux2_closure,
  heat_flux_flux3,
  heat_flux_flux3_gradients,
  heat_flux_flux3_face_state,
  heat_flux_flux3_closure,
  heat_flux_work_diagnostics,
  timestep_reduction,
  sweep_begin_conversion,
  sts_clear_flux,
  sts_update_copies,
  sts_update_kernel,
  primitive_refresh,
  admissibility,
  sweep_end_conversion,
  post_sweep_collisions,
  parabolic_init_recv,
  parabolic_send_flux,
  parabolic_recv_flux,
  parabolic_restrict_u,
  parabolic_send_u,
  parabolic_recv_u,
  parabolic_physical_bcs,
  parabolic_prolongate,
  count
};

constexpr int kCGLLFProfileBucketCount =
    static_cast<int>(CGLLFProfileBucket::count);

class CGLLandauFluid {
 public:
  CGLLandauFluid(MeshBlockPack *pp, ParameterInput *pin);

  Real dtnew;
  Real lf_k_parallel;
  bool lf_coeff_local;
  Real lf_c_parallel0;
  bool strict_admissibility;
  bool effective_backup_limiter;
  CGLLFDiagnosticsMode diagnostics_mode;
  CGLLFArithmeticMode arithmetic_mode;
  CGLLFSTSFluxMode sts_flux_mode;
  parabolic::ParabolicIntegratorMode mode;
  CGLLFDiagnostics diagnostics;

  void AddHeatFluxes(const DvceArray5D<Real> &w, const DvceArray5D<Real> &bcc,
                     const EOS_Data &eos, Real dt_sweep, Real rkl_weight,
                     DvceFaceFld5D<Real> &f);
  void AdvanceHeatFluxWorkDiagnostics(const parabolic::RKL2Coefficients &coeffs,
                                      int stage, int nstages);
  void AdvancePressureWorkDiagnostics(Real beta_dt, Real gam0, Real gam1, int stage,
                                      Real pressure_power, Real anisotropic_power);
  void ResetHeatFluxDiagnostics();
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos);
  void RecordAdmissibility(const DvceArray5D<Real> &u, const DvceArray5D<Real> &w,
                           const DvceArray5D<Real> &bcc, const EOS_Data &eos,
                           int dfloor_delta, int pfloor_delta,
                           const char *sweep_name, int stage, int nstages);
  bool ProfileEnabled() const {return profile_enabled_;}
  void AddProfileTime(CGLLFProfileBucket bucket, Real seconds);
  void ReportProfile(const char *context) const;
  bool UsesWeightedSTSFlux() const {
    return sts_flux_mode == CGLLFSTSFluxMode::weighted;
  }

 private:
  void AccumulateHeatFluxDiagnostics(const array_sum::GlobalSum &stats);

  MeshBlockPack *pmy_pack;
  bool profile_enabled_ = false;
  bool profile_detail_enabled_ = false;
  int profile_last_nstages_ = 0;
  int profile_max_nstages_ = 0;
  Real profile_seconds_[kCGLLFProfileBucketCount] = {};
  std::uint64_t profile_counts_[kCGLLFProfileBucketCount] = {};
  Real profile_detail_sink_ = 0.0;
  DvceArray4D<Real> tpar_, tperp_, bmag_;
  Real stage_qpar_work_ = 0.0;
  Real stage_qperp_work_ = 0.0;
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

class CGLLFProfileRegion {
 public:
  CGLLFProfileRegion(CGLLandauFluid *profile, CGLLFProfileBucket bucket);
  ~CGLLFProfileRegion();
  CGLLFProfileRegion(const CGLLFProfileRegion&) = delete;
  CGLLFProfileRegion& operator=(const CGLLFProfileRegion&) = delete;

 private:
  CGLLandauFluid *profile_;
  CGLLFProfileBucket bucket_;
  Kokkos::Timer timer_;
  bool active_;
};

#endif // DIFFUSION_CGL_LANDAU_FLUID_HPP_
