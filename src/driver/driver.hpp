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

#include <array>
#include <cstdint>
#include <ctime>
#include <memory>
#include <string>

#include "parameter_input.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

// Coarse performance regions used by opt-in, synchronization-based diagnostics.
// Regions may be nested (for example cooling inside stage_tasks), so their times
// are not intended to sum to the total runtime.
enum class PerformanceRegion : int {
  before_integrator = 0,
  before_stage,
  stage_tasks,
  after_stage,
  after_integrator,
  outputs,
  amr,
  timestep,
  trml_cooling,
  tracer_save_flux,
  frame_total,
  frame_control,
  frame_boundary,
  frame_timestep,
  particle_push,
  particle_comm,
  particle_adjust,
  particle_seed,
  count
};

//----------------------------------------------------------------------------------------
//! \class Driver

class Driver {
 public:
  Driver(ParameterInput *pin, Mesh *pmesh, Real wtlim, Kokkos::Timer* ptimer);
  ~Driver() = default;

  // data
  TimeEvolution time_evolution;
  DvceArray6D<Real> impl_src;  // stiff source terms used in ImEx integrators

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
  bool performance_timing = false;

  // functions
  void ExecuteTaskList(Mesh *pm, std::string tl, int stage);
  void Initialize(Mesh *pmesh, ParameterInput *pin, Outputs *pout, bool rflag);
  void Execute(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
  void Finalize(Mesh *pmesh, ParameterInput *pin, Outputs *pout);
  void InitBoundaryValuesAndPrimitives(Mesh *pm);
  void StartPerformanceRegion(PerformanceRegion region);
  void StopPerformanceRegion(PerformanceRegion region);

 private:
  Kokkos::Timer run_time_;      // generalized timer for cpu/gpu/etc
  std::uint64_t nmb_updated_;   // running total of MB updated during run
  std::uint64_t npart_updated_; // running total of particles updated during run
  float lb_efficiency_;         // measure of how efficient was load balancing
  static constexpr int nperformance_regions_ =
      static_cast<int>(PerformanceRegion::count);
  std::array<Kokkos::Timer, nperformance_regions_> performance_timers_;
  std::array<double, nperformance_regions_> performance_seconds_{};
  std::array<std::uint64_t, nperformance_regions_> performance_calls_{};
  void ResetPerformanceTiming();
  void ReportPerformanceTiming(double run_seconds);
  void OutputCycleDiagnostics(Mesh *pm);
  Real UpdateWallClock();
};
#endif // DRIVER_DRIVER_HPP_
