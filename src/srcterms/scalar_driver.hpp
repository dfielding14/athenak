#ifndef SRCTERMS_SCALAR_DRIVER_HPP_
#define SRCTERMS_SCALAR_DRIVER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_driver.hpp
//  \brief defines stochastic Ornstein-Uhlenbeck forcing for a passive scalar

#include <memory>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/random.hpp"

enum class ScalarForcingNormalization { variance_rate, source_rms };
enum class ScalarForcingSpectrum { parabolic, power_law };

struct ScalarForcingRestartMetadata {
  int version;
  int mode_count;
  int n_updates;
  int scalar_index;
  int nlow;
  int nhigh;
  int normalization;
  int spectrum;
  Real tcorr;
  Real dt_update;
  Real variance_rate;
  Real source_rms;
  Real npeak;
  Real expo;
};

//----------------------------------------------------------------------------------------
//! \class ScalarForcingDriver

class ScalarForcingDriver {
 public:
  ScalarForcingDriver(MeshBlockPack *pp, ParameterInput *pin);
  ~ScalarForcingDriver();

  DvceArray5D<Real> force;
  RNG_State rstate;
  DualArray1D<Real> mode_amp_real;
  DualArray1D<Real> mode_amp_imag;
  DualArray1D<Real> mode_noise_real;
  DualArray1D<Real> mode_noise_imag;
  DualArray1D<Real> mode_weight;
  DualArray1D<Real> kx_mode;
  DualArray1D<Real> ky_mode;
  DualArray1D<Real> kz_mode;
  DvceArray3D<Real> xcos, xsin, ycos, ysin, zcos, zsin;

  int mode_count;
  int n_updates_yet;
  int scalar_index;

  void IncludeInitializeModesTask(std::shared_ptr<TaskList> tl, TaskID start);
  void IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start);
  TaskStatus InitializeModes(Driver *pdrive, int stage);
  TaskStatus EnsureBasisSize(Driver *pdrive, int stage);
  TaskStatus UpdateForcing(Driver *pdrive, int stage);
  TaskStatus AddForcing(Driver *pdrive, int stage);

  ScalarForcingRestartMetadata RestartMetadata() const;
  void ValidateRestartMetadata(const ScalarForcingRestartMetadata &metadata) const;

 private:
  void Initialize();
  void BuildBasis();
  void RenderForce();
  void ApplyForcingWithStep(Real bdt);

  MeshBlockPack *pmy_pack;

  int nlow, nhigh;
  int rseed;
  int dimension;
  Real tcorr, dt_update;
  Real variance_rate, source_rms;
  Real npeak, expo;
  ScalarForcingNormalization normalization;
  ScalarForcingSpectrum spectrum;

  Real lx, ly, lz;
  int current_nmb_;
  int last_nmb_created_;
  int last_nmb_deleted_;
};

#endif  // SRCTERMS_SCALAR_DRIVER_HPP_
