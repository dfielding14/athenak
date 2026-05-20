#ifndef PGEN_PGEN_HPP_
#define PGEN_PGEN_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.hpp
//  \brief definitions for ProblemGenerator class

#include <functional>
#include <memory>
#include <vector>

#include "geodesic-grid/spherical_grid.hpp"
#include "parameter_input.hpp"

struct EOS_Data;

namespace cooling {
struct RuntimeData;
struct TimestepData;
}

using ProblemFinalizeFnPtr = void (*)(ParameterInput *pin, Mesh *pm);
using UserBoundaryFnPtr = void (*)(Mesh* pm);
using UserSrctermFnPtr = void (*)(Mesh* pm, const Real bdt);
using UserTimestepFnPtr = void (*)(Mesh* pm);
using UserRefinementFnPtr = void (*)(MeshBlockPack* pmbp);
using UserHistoryFnPtr = void (*)(HistoryData *pdata, Mesh *pm);
using UserCoolingSourceFnPtr = void (*)(MeshBlockPack *pmbp,
    const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
    const cooling::RuntimeData &runtime, const Real bdt,
    DvceArray5D<Real> &u0, Real &gross_energy, Real &net_energy);
using UserCoolingTimeStepFnPtr = void (*)(MeshBlockPack *pmbp,
    const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
    const cooling::RuntimeData &runtime, const cooling::TimestepData &timestep,
    Real &dtnew);

//----------------------------------------------------------------------------------------
//! \class ProblemGenerator

class ProblemGenerator {
 public:
  // constructor for new problems
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh);
  // constructor for restarts
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh, IOWrapper resfile,
                   bool single_file_per_rank=false);
  ~ProblemGenerator() = default;

  // true if user BCs are specified on any face
  bool user_bcs;

  // true if user srcterms are specified
  bool user_srcs;

  // true if user timestep constraints are specified
  bool user_dt;

  // storage for user timestep constraints
  Real dtnew = 1e10;

  // true if user history outputs are specified
  bool user_hist;

  // vector of SphericalGrid objects for analysis
  std::vector<std::unique_ptr<SphericalGrid>> spherical_grids;

  // function pointer for final work after main loop (e.g. compute errors).  Called by
  // Driver::Finalize()
  ProblemFinalizeFnPtr pgen_final_func=nullptr;
  // function pointer for user-enrolled BCs.  Called in ApplyPhysicalBCs in task list
  UserBoundaryFnPtr user_bcs_func=nullptr;
  UserSrctermFnPtr user_srcs_func=nullptr;
  UserTimestepFnPtr user_time_step_func=nullptr;
  UserRefinementFnPtr user_ref_func=nullptr;
  UserHistoryFnPtr user_hist_func=nullptr;
  // Optional hooks used when <cooling>/cooling_model or heating_model is "user".
  // UserProblem() should set these to host launcher functions that run any needed
  // device kernels and report domain-integrated gross/net cooling energies.
  UserCoolingSourceFnPtr user_cooling_func=nullptr;
  UserCoolingTimeStepFnPtr user_cooling_timestep_func=nullptr;

  // predefined problem generator functions (default test suite)
  void CallProblemGenerator(ParameterInput *pin, bool is_restart);
  void Advection(ParameterInput *pin, const bool restart);
  void AlfvenWave(ParameterInput *pin, const bool restart);
  void BondiAccretion(ParameterInput *pin, const bool restart);
  void CShock(ParameterInput *pin, const bool restart);
  void DivBAMR(ParameterInput *pin, const bool restart);
  void Diffusion(ParameterInput *pin, const bool restart);
  void LinearWave(ParameterInput *pin, const bool restart);
  void LWImplode(ParameterInput *pin, const bool restart);
  void Monopole(ParameterInput *pin, const bool restart);
  void MRI3d(ParameterInput *pin, const bool restart);
  void OrszagTang(ParameterInput *pin, const bool restart);
  void ShockTube(ParameterInput *pin, const bool restart);
  void Shwave(ParameterInput *pin, const bool restart);
  void SphericalCollapse(ParameterInput *pin, const bool restart);
  void RadiationLinearWave(ParameterInput *pin, const bool restart);
  void RadiationBeam(ParameterInput *pin, const bool restart);
  void Z4cBoostedPuncture(ParameterInput *pin, const bool restart);
  void Z4cLinearWave(ParameterInput *pin, const bool restart);

  // Generic error output function (using difference u0-u1)
  void OutputErrors(ParameterInput *pin, Mesh *pm);

  // template for user-specified problem generator
  void UserProblem(ParameterInput *pin, const bool restart);

 private:
  bool single_file_per_rank; // for restart file naming
  Mesh* pmy_mesh_;
};

#endif // PGEN_PGEN_HPP_
