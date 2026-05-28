#ifndef SRCTERMS_TURB_DRIVER_HPP_
#define SRCTERMS_TURB_DRIVER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver.hpp
//  \brief defines the stochastic turbulence forcing driver

#include <memory>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/random.hpp"

enum class TurbNormalization { edot, accel_rms };
enum class TurbLocalization { none, include, exclude };
enum class TurbSpectrum { parabolic, power_law };

// Native restart records are only intended for restarts from the same executable
// precision and feature version, consistent with the existing AthenaK restart format.
struct TurbulenceRestartMetadata {
  int version;
  int mode_count;
  int n_updates;
  int nlow;
  int nhigh;
  int driving_type;
  int min_kx;
  int max_kx;
  int min_ky;
  int max_ky;
  int min_kz;
  int max_kz;
  int use_npeak;
  int turb_flag;
  int tile_nx;
  int tile_ny;
  int tile_nz;
  int normalization;
  int localization;
  int spectrum;
  Real tcorr;
  Real dt_update;
  Real dedt;
  Real accel_rms;
  Real sol_fraction;
  Real kpeak;
  Real npeak;
  Real expo;
  Real exp_prp;
  Real exp_prl;
  Real tdriv_duration;
  Real tdriv_start;
  Real sigma_x1;
  Real sigma_x2;
  Real sigma_x3;
  Real center_x1;
  Real center_x2;
  Real center_x3;
};

//----------------------------------------------------------------------------------------
//! \class TurbulenceDriver

class TurbulenceDriver {
 public:
  TurbulenceDriver(MeshBlockPack* pp, ParameterInput* pin);
  ~TurbulenceDriver();

  DvceArray5D<Real> force;                  // normalized acceleration on the mesh
  RNG_State rstate;                         // RNG state for modal OU innovations
  DualArray2D<Real> mode_amp_real;          // authoritative real modal OU state
  DualArray2D<Real> mode_amp_imag;          // authoritative imaginary modal OU state
  DualArray2D<Real> mode_noise_real;        // current modal innovation
  DualArray2D<Real> mode_noise_imag;        // current modal innovation
  DualArray1D<Real> kx_mode, ky_mode, kz_mode;
  DvceArray3D<Real> xcos, xsin, ycos, ysin, zcos, zsin;

  int mode_count;
  int n_turb_updates_yet;

  void IncludeInitializeModesTask(std::shared_ptr<TaskList> tl, TaskID start);
  void IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start);
  TaskStatus InitializeModes(Driver* pdrive, int stage);
  TaskStatus EnsureBasisSize(Driver* pdrive, int stage);
  TaskStatus UpdateForcing(Driver* pdrive, int stage);
  TaskStatus AddForcing(Driver* pdrive, int stage);
  void ApplyForcingWithStep(Real bdt);

  TurbulenceRestartMetadata RestartMetadata() const;
  void ValidateRestartMetadata(const TurbulenceRestartMetadata& metadata) const;

 private:
  void Initialize();
  void BuildBasis();
  void RenderForce();

  MeshBlockPack* pmy_pack;  // MeshBlockPack containing this driver

  int nlow, nhigh;
  int rseed;
  Real kpeak, npeak;
  bool use_npeak;
  Real tcorr, dedt, tdriv_duration, tdriv_start;
  Real expo, exp_prl, exp_prp;
  int driving_type, turb_flag;
  int min_kz, max_kz, min_kx, max_kx, min_ky, max_ky;
  Real sol_fraction;
  Real dt_update;
  TurbNormalization normalization;
  Real accel_rms;
  TurbSpectrum spectrum;

  TurbLocalization localization;
  Real sigma_x1, sigma_x2, sigma_x3;
  Real center_x1, center_x2, center_x3;

  int tile_nx, tile_ny, tile_nz;
  Real tile_lx, tile_ly, tile_lz;
  Real inv_tile_lx, inv_tile_ly, inv_tile_lz;
  Real domain_x1min, domain_x2min, domain_x3min;

  int current_nmb_;
  int last_nmb_created_;
  int last_nmb_deleted_;
};

#endif  // SRCTERMS_TURB_DRIVER_HPP_
