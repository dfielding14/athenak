#ifndef DRIVER_DRIVER_HPP_
#define DRIVER_DRIVER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file driver.hpp
//  \brief definitions for Driver class
//
// Note ProblemGenerator object is stored in Driver and is called in Initialize(). If the
// pgen class contains analysis routines that are run at end of execution, they can be
// called in Finalize().

#include <ctime>
#include <limits>
#include <memory>
#include <string>

#include "diffusion/sts_types.hpp"
#include "diffusion/sts_rkl2.hpp"
#include "parameter_input.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

//----------------------------------------------------------------------------------------
//! \class Driver

class Driver {
 public:
  enum class STSSweep {none, pre, post};

  struct STSController {
    bool enabled = false;
    bool explicit_split = false;
    parabolic::STSIntegrator integrator = parabolic::STSIntegrator::none;
    STSSweep sweep = STSSweep::none;
    Real dt_cycle = 0.0;
    Real dt_sweep = 0.0;
    Real dt_parabolic_min = std::numeric_limits<float>::max();
    int nstages = 0;
    int current_stage = 0;
    // Internal LF chunks share one representation conversion and collision hook.
    bool first_chunk = true, last_chunk = true;
    parabolic::RKL2Coefficients coeffs;
  };

  Driver(ParameterInput *pin, Mesh *pmesh, Real wtlim, Kokkos::Timer* ptimer);
  ~Driver() = default;

  // data
  TimeEvolution time_evolution;
  DvceArray6D<Real> impl_src;  // stiff source terms used in ImEx integrators
  STSController sts;

  // folowing data only relevant for runs involving time evolution
  Real tlim;      // stopping time
  int nlim;       // cycle-limit
  int ndiag;      // cycles between output of diagnostic information
  // variables for various SSP and ImEx RK integrators
  std::string integrator;          // integrator name (rk1, rk2, rk3)
  int nimp_stages;                 // number of implicit stages (ImEx only)
  int nexp_stages;                 // number of explicit stages (both SSP-RK and ImEx)
  Real gam0[4], gam1[4], beta[4];  // weights and fractional timestep per explicit stage
  Real delta[4];                   // weights for updating the intermediate stage (u1)
  Real a_twid[4][4], a_impl;       // matrix elements for implicit stages in ImEx
  Real cfl_limit;                  // maximum CFL number for integrator
  Real gamma;                      // gamma value for the IMEX_new integrator
  Kokkos::Timer* pwall_clock_;     // timer for tracking the wall clock
  Real wall_time;

  // functions
  void ExecuteTaskList(Mesh *pm, std::string tl, int stage);
  void Initialize(Mesh *pmesh, ParameterInput *pin, Outputs *pout, bool rflag);
  void Execute(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
  void Finalize(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
  void InitBoundaryValuesAndPrimitives(Mesh *pm);

 private:
  Kokkos::Timer run_time_;      // generalized timer for cpu/gpu/etc
  std::uint64_t nmb_updated_;   // running total of MB updated during run
  std::uint64_t npart_updated_; // running total of particles updated during run
  std::uint64_t last_diag_nmb_updated_; // global MB updates at previous diagnostic
  int last_diag_cycle_;         // completed cycles at previous diagnostic
  double last_diag_time_;       // rank-zero wall seconds at previous diagnostic
  float lb_efficiency_;         // measure of how efficient was load balancing
  // Optional, collisionless LF-only transaction; restart files never store debt.
  bool merge_sts_requested_ = false, merge_sts_enabled_ = false;
  Real pending_sts_half_ = 0.0, pending_sts_cycle_ = 0.0, pending_sts_time_ = 0.0;
  Real merge_deferred_ = 0.0, merge_consumed_ = 0.0, merge_flushed_ = 0.0;
  Real merge_snapshot_seconds_ = 0.0;
  std::uint64_t merge_accepted_ = 0, merge_rejected_ = 0;
  std::uint64_t merge_rejected_cfl_ = 0, merge_rejected_admissibility_ = 0;
  std::uint64_t merge_accepted_stages_ = 0, merge_rejected_stages_ = 0;
  DvceArray5D<Real> merge_u_backup_, merge_w_backup_;
  Real cgl_lf_max_chunk_ratio_ = 0.0;
  std::uint64_t lf_chunk_sweeps_ = 0, lf_chunks_ = 0, lf_chunk_rhs_ = 0;
  void ConfigureMergedSTS(Mesh *pm);
  bool CanDeferSTSPost(Mesh *pm, Outputs *pout) const;
  void RunMergedSTSSweep(Mesh *pm, Real duration, Real old_cycle);
  bool TryMergedSTSPre(Mesh *pm);
  void FlushPendingSTS(Mesh *pm, bool select_next_dt=false);
  void ResetSTSController();
  void ValidateSTSConfiguration(Mesh *pm);
  void RefreshSTSCycleState(Mesh *pm);
  void BeginSTSSweep(Mesh *pm, STSSweep sweep);
  void RunSTSSweep(Mesh *pm, STSSweep sweep);
  void SetSTSStage(int stage);
  void EndSTSSweep();
  void OutputCycleDiagnostics(Mesh *pm);
  Real UpdateWallClock();
};
#endif // DRIVER_DRIVER_HPP_
