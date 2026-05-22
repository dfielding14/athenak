#ifndef SRCTERMS_INITIAL_PERTURBATIONS_HPP_
#define SRCTERMS_INITIAL_PERTURBATIONS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file initial_perturbations.hpp
//! \brief One-time Fourier perturbations applied to initial conditions.

#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/random.hpp"

class InitialPerturbations {
 public:
  InitialPerturbations(MeshBlockPack *pp, ParameterInput *pin,
                       const std::string &block_name = "initial_perturbations");

  bool AnyEnabled() const {
    return DensityEnabled() || VelocityEnabled() || MagneticEnabled();
  }
  void Apply();

 private:
  bool DensityEnabled() const {
    return perturb_density && density_rms > 0.0;
  }
  bool VelocityEnabled() const {
    return perturb_velocity && velocity_rms > 0.0;
  }
  bool MagneticEnabled() const {
    return perturb_magnetic && magnetic_rms > 0.0;
  }

  void InitializeModes();
  void GenerateAmplitudes(DualArray2D<Real> &amp_real, DualArray2D<Real> &amp_imag,
                          const int ncomp, const bool project_velocity,
                          const bool vector_potential);
  void BuildDensityField();
  void BuildVelocityField();
  void BuildVectorPotential();
  void CurlVectorPotential();
  void ApplyDensityField();
  void ApplyVelocityField();
  void ApplyMagneticField();
  Real CellScalarRMS(const DvceArray4D<Real> &field, const bool remove_mean);
  Real CellVectorRMS(const DvceArray5D<Real> &field, const bool remove_mean);
  Real FaceCenteredBRMS(const DvceFaceFld4D<Real> &field);
  void Validate() const;

  MeshBlockPack *pmy_pack;
  const std::string block_name_;

  bool perturb_density, perturb_velocity, perturb_magnetic;
  bool density_fractional;
  bool remove_density_mean, remove_velocity_mean;
  Real density_rms, velocity_rms, magnetic_rms;

  int nlow, nhigh;
  int min_kx, max_kx, min_ky, max_ky, min_kz, max_kz;
  int mode_count;
  Real spectral_slope;
  Real f_solenoidal;
  int rseed;
  RNG_State rstate;

  int localization_mode;  // 0=none, 1=include Gaussian region, 2=exclude Gaussian region
  Real x1_center, x2_center, x3_center;
  Real x1_scale, x2_scale, x3_scale;

  Real domain_x1min, domain_x2min, domain_x3min;
  Real lx, ly, lz;

  DualArray1D<Real> kx_mode, ky_mode, kz_mode;
  DualArray2D<Real> rho_amp_real, rho_amp_imag;
  DualArray2D<Real> vel_amp_real, vel_amp_imag;
  DualArray2D<Real> avec_amp_real, avec_amp_imag;

  DvceArray4D<Real> rho_pert;
  DvceArray5D<Real> vel_pert;
  DvceEdgeFld4D<Real> avec;
  DvceFaceFld4D<Real> b_pert;
};

void ApplyInitialPerturbations(Mesh *pm, ParameterInput *pin);

#endif  // SRCTERMS_INITIAL_PERTURBATIONS_HPP_
