//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file TRML_frame_tracking.cpp
//! \brief Minimal turbulent radiative mixing-layer setup demonstrating FrameTracker use.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "parameter_input.hpp"
#include "pgen.hpp"
#include "srcterms/frame_tracker.hpp"

namespace {

constexpr Real kPi = 3.141592653589793238462643383279502884;

struct TRMLFrameTrackingData {
  Real gm1 = 2.0/3.0;
  Real rho_hot = 1.0;
  Real pressure = 1.0;
  Real density_contrast = 100.0;
  Real shear_velocity = 1.0;
  Real drift_velocity_x1 = 0.0;
  Real drift_velocity_x2 = 0.0;
  Real drift_velocity_x3 = 0.0;
  Real z_interface = 0.0;
  Real smoothing_thickness = 0.02;
  Real perturbation_amplitude = 0.01;
  int perturbation_mode_x1 = 1;
  int perturbation_mode_x2 = 1;
  Real x1min = -0.5;
  Real x1max = 0.5;
  Real x2min = -0.5;
  Real x2max = 0.5;

  bool cooling_enabled = true;
  Real t_cool_start = 0.0;
  Real t_cool_0 = 1.0;
  Real cfl_cool = 0.25;
  Real T_peak_over_T_cold = 3.16227766017;
  Real T_cutoff_over_T_hot = 1.0;
  Real beta_lo = 1.0;
  Real beta_hi = 1.0;
  Real T_floor = 0.0;
  Real T_cold = 0.01;
  Real T_hot = 1.0;
  Real T_peak = 0.0316227766017;
  Real T_cutoff = 1.0;
};

TRMLFrameTrackingData trml_data;

void FatalTRMLInput(const std::string &message) {
  std::cout << "### FATAL ERROR in TRML_frame_tracking input" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real CoolingShape(const Real temperature, const TRMLFrameTrackingData &data) {
  if (!data.cooling_enabled || temperature <= data.T_floor ||
      temperature >= data.T_cutoff) {
    return 0.0;
  }
  const Real T = fmax(temperature, 1.0e-30);
  if (T < data.T_peak) {
    return pow(T/data.T_peak, data.beta_lo);
  }
  return pow(data.T_peak/T, data.beta_hi);
}

KOKKOS_INLINE_FUNCTION
void SetHydroState(const DvceArray5D<Real> &u0, const int m, const int k,
                   const int j, const int i, const Real density,
                   const Real pressure, const Real vx, const Real vy,
                   const Real vz, const Real gm1) {
  u0(m, IDN, k, j, i) = density;
  u0(m, IM1, k, j, i) = density*vx;
  u0(m, IM2, k, j, i) = density*vy;
  u0(m, IM3, k, j, i) = density*vz;
  u0(m, IEN, k, j, i) = pressure/gm1 +
      0.5*density*(SQR(vx) + SQR(vy) + SQR(vz));
}

KOKKOS_INLINE_FUNCTION
void EvaluateTRMLState(const Real x1, const Real x2, const Real x3,
                       const TRMLFrameTrackingData &data,
                       Real &density, Real &pressure,
                       Real &vx, Real &vy, Real &vz, Real &cold_fraction) {
  const Real width = fmax(data.smoothing_thickness, 1.0e-30);
  const Real profile = tanh((x3 - data.z_interface)/width);
  cold_fraction = 0.5*(1.0 - profile);
  density = data.rho_hot*(1.0 + cold_fraction*(data.density_contrast - 1.0));
  pressure = data.pressure;

  vx = data.drift_velocity_x1 + 0.5*data.shear_velocity*profile;
  vy = data.drift_velocity_x2;
  vz = data.drift_velocity_x3;

  const Real Lx1 = data.x1max - data.x1min;
  const Real Lx2 = data.x2max - data.x2min;
  if (data.perturbation_amplitude != 0.0 && Lx1 > 0.0 && Lx2 > 0.0) {
    const Real phase1 = 2.0*kPi*data.perturbation_mode_x1*
                        (x1 - data.x1min)/Lx1;
    const Real phase2 = 2.0*kPi*data.perturbation_mode_x2*
                        (x2 - data.x2min)/Lx2;
    const Real envelope = exp(-SQR((x3 - data.z_interface)/(4.0*width)));
    vz += data.perturbation_amplitude*data.shear_velocity*
          sin(phase1)*sin(phase2)*envelope;
  }
}

void TRMLZBoundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const TRMLFrameTrackingData data = trml_data;

  Real frame_x1 = 0.0;
  Real frame_x2 = 0.0;
  Real frame_x3 = 0.0;
  Real frame_v1 = 0.0;
  Real frame_v2 = 0.0;
  Real frame_v3 = 0.0;
  if (pmbp->pframe_tracker != nullptr) {
    frame_x1 = pmbp->pframe_tracker->FrameDisplacement(0);
    frame_x2 = pmbp->pframe_tracker->FrameDisplacement(1);
    frame_x3 = pmbp->pframe_tracker->FrameDisplacement(2);
    frame_v1 = pmbp->pframe_tracker->FrameVelocity(0);
    frame_v2 = pmbp->pframe_tracker->FrameVelocity(1);
    frame_v3 = pmbp->pframe_tracker->FrameVelocity(2);
  }

  par_for("trml_inner_x3", DevExeSpace(), 0, nmb1, 0, ng-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int kg, const int j, const int i) {
    if (mb_bcs.d_view(m, BoundaryFace::inner_x3) != BoundaryFlag::user) {
      return;
    }
    const RegionSize block_size = size.d_view(m);
    const int kb = ks - kg - 1;
    const Real x1 = CellCenterX(i - is, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, block_size.x2min, block_size.x2max);
    const Real x3 = CellCenterX(kb - ks, indcs.nx3, block_size.x3min, block_size.x3max);
    Real density, pressure, vx, vy, vz, cold_fraction;
    EvaluateTRMLState(x1 + frame_x1, x2 + frame_x2, x3 + frame_x3, data,
                      density, pressure, vx, vy, vz, cold_fraction);
    SetHydroState(u0, m, kb, j, i, density, pressure, vx - frame_v1,
                  vy - frame_v2, vz - frame_v3, data.gm1);
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      u0(m, n, kb, j, i) = density*cold_fraction;
    }
  });

  par_for("trml_outer_x3", DevExeSpace(), 0, nmb1, 0, ng-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int kg, const int j, const int i) {
    if (mb_bcs.d_view(m, BoundaryFace::outer_x3) != BoundaryFlag::user) {
      return;
    }
    const RegionSize block_size = size.d_view(m);
    const int kb = indcs.ke + kg + 1;
    const Real x1 = CellCenterX(i - is, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, block_size.x2min, block_size.x2max);
    const Real x3 = CellCenterX(kb - ks, indcs.nx3, block_size.x3min, block_size.x3max);
    Real density, pressure, vx, vy, vz, cold_fraction;
    EvaluateTRMLState(x1 + frame_x1, x2 + frame_x2, x3 + frame_x3, data,
                      density, pressure, vx, vy, vz, cold_fraction);
    SetHydroState(u0, m, kb, j, i, density, pressure, vx - frame_v1,
                  vy - frame_v2, vz - frame_v3, data.gm1);
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      u0(m, n, kb, j, i) = density*cold_fraction;
    }
  });
}

void TRMLCoolingSource(Mesh *pm, const Real bdt) {
  if (!trml_data.cooling_enabled || pm->time < trml_data.t_cool_start) {
    return;
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  const TRMLFrameTrackingData data = trml_data;

  par_for("trml_cooling", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real density = w0(m, IDN, k, j, i);
    if (density <= 0.0) {
      return;
    }
    const Real eint = u0(m, IEN, k, j, i) -
        0.5*(SQR(u0(m, IM1, k, j, i)) + SQR(u0(m, IM2, k, j, i)) +
             SQR(u0(m, IM3, k, j, i)))/density;
    const Real temperature = data.gm1*eint/density;
    const Real shape = CoolingShape(temperature, data);
    if (shape <= 0.0) {
      return;
    }
    const Real floor_eint = density*data.T_floor/data.gm1;
    const Real max_de = fmax(0.0, eint - floor_eint);
    const Real requested_de = bdt*eint*shape/data.t_cool_0;
    u0(m, IEN, k, j, i) -= fmin(max_de, requested_de);
  });
}

void TRMLCoolingTimeStep(Mesh *pm) {
  if (!trml_data.cooling_enabled || pm->time < trml_data.t_cool_start) {
    pm->pgen->dtnew = std::numeric_limits<float>::max();
    return;
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  const TRMLFrameTrackingData data = trml_data;
  Real dtnew = static_cast<Real>(std::numeric_limits<float>::max());

  Kokkos::parallel_reduce("trml_cooling_newdt",
  Kokkos::RangePolicy<>(DevExeSpace(), 0,
                        pmbp->nmb_thispack*indcs.nx3*indcs.nx2*indcs.nx1),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    const int nkji = indcs.nx3*indcs.nx2*indcs.nx1;
    const int nji = indcs.nx2*indcs.nx1;
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji + ks;
    int j = (idx - m*nkji - (k - ks)*nji)/indcs.nx1 + js;
    int i = (idx - m*nkji - (k - ks)*nji - (j - js)*indcs.nx1) + is;

    const Real density = w0(m, IDN, k, j, i);
    if (density <= 0.0) {
      return;
    }
    const Real eint = u0(m, IEN, k, j, i) -
        0.5*(SQR(u0(m, IM1, k, j, i)) + SQR(u0(m, IM2, k, j, i)) +
             SQR(u0(m, IM3, k, j, i)))/density;
    const Real temperature = data.gm1*eint/density;
    const Real shape = CoolingShape(temperature, data);
    if (shape > 0.0) {
      min_dt = fmin(min_dt, data.cfl_cool*data.t_cool_0/shape);
    }
  }, Kokkos::Min<Real>(dtnew));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &dtnew, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  pm->pgen->dtnew = dtnew;
}

void ReadTRMLParameters(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) {
    FatalTRMLInput("TRML_frame_tracking requires a <hydro> block.");
  }
  if (!pmbp->phydro->peos->eos_data.is_ideal) {
    FatalTRMLInput("TRML_frame_tracking requires an ideal-gas hydro EOS.");
  }

  TRMLFrameTrackingData data;
  data.gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  data.rho_hot = pin->GetOrAddReal("problem", "rho_0", 1.0);
  data.pressure = pin->GetOrAddReal("problem", "pgas_0", 1.0);
  data.density_contrast = pin->GetOrAddReal("problem", "density_contrast", 100.0);
  data.shear_velocity = pin->GetOrAddReal("problem", "velocity", 1.0);
  data.drift_velocity_x1 = pin->GetOrAddReal("problem", "drift_velocity_x1", 0.0);
  data.drift_velocity_x2 = pin->GetOrAddReal("problem", "drift_velocity_x2", 0.0);
  data.drift_velocity_x3 = pin->GetOrAddReal("problem", "drift_velocity_x3", 0.1);
  data.z_interface = pin->GetOrAddReal("problem", "z_interface", 0.0);
  data.smoothing_thickness = pin->GetOrAddReal("problem", "smoothing_thickness", 0.02);
  data.perturbation_amplitude =
      pin->GetOrAddReal("problem", "perturbation_amplitude", 0.01);
  data.perturbation_mode_x1 = pin->GetOrAddInteger("problem", "perturbation_mode_x1", 1);
  data.perturbation_mode_x2 = pin->GetOrAddInteger("problem", "perturbation_mode_x2", 1);
  data.x1min = pm->mesh_size.x1min;
  data.x1max = pm->mesh_size.x1max;
  data.x2min = pm->mesh_size.x2min;
  data.x2max = pm->mesh_size.x2max;
  data.cooling_enabled = pin->GetOrAddBoolean("problem", "cooling_enabled", true);
  data.t_cool_start = pin->GetOrAddReal("problem", "t_cool_start", 0.0);
  data.t_cool_0 = pin->GetOrAddReal("problem", "t_cool_0", 1.0);
  data.cfl_cool = pin->GetOrAddReal("problem", "cfl_cool", 0.25);
  data.T_peak_over_T_cold = pin->GetOrAddReal("problem", "T_peak_over_T_cold",
                                              3.16227766017);
  data.T_cutoff_over_T_hot = pin->GetOrAddReal("problem", "T_cutoff_over_T_hot", 1.0);
  data.beta_lo = pin->GetOrAddReal("problem", "beta_lo", 1.0);
  data.beta_hi = pin->GetOrAddReal("problem", "beta_hi", 1.0);

  if (data.rho_hot <= 0.0 || data.pressure <= 0.0 || data.density_contrast <= 1.0) {
    FatalTRMLInput("Require rho_0 > 0, pgas_0 > 0, and density_contrast > 1.");
  }
  if (data.smoothing_thickness <= 0.0) {
    FatalTRMLInput("Require problem/smoothing_thickness > 0.");
  }
  if (data.cooling_enabled && (data.t_cool_0 <= 0.0 || data.cfl_cool <= 0.0)) {
    FatalTRMLInput("Cooling requires t_cool_0 > 0 and cfl_cool > 0.");
  }

  data.T_hot = data.pressure/data.rho_hot;
  data.T_cold = data.T_hot/data.density_contrast;
  data.T_peak = data.T_peak_over_T_cold*data.T_cold;
  data.T_cutoff = data.T_cutoff_over_T_hot*data.T_hot;
  data.T_floor = pin->GetOrAddReal("problem", "T_floor", 0.5*data.T_cold);
  if (data.T_peak <= 0.0 || data.T_cutoff <= data.T_floor) {
    FatalTRMLInput("Require positive T_peak and T_cutoff > T_floor.");
  }

  trml_data = data;

  if (global_variable::my_rank == 0) {
    std::cout << "TRML_frame_tracking: rho_hot=" << data.rho_hot
              << ", rho_cold=" << data.rho_hot*data.density_contrast
              << ", T_hot=" << data.T_hot
              << ", T_cold=" << data.T_cold
              << ", cooling=" << (data.cooling_enabled ? "on" : "off")
              << std::endl;
  }
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem()
//! \brief Initialize a pressure-balanced TRML and enroll moving-frame x3 boundaries.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  ReadTRMLParameters(pin, pmy_mesh_);

  user_bcs_func = TRMLZBoundary;
  if (trml_data.cooling_enabled) {
    user_srcs = true;
    user_srcs_func = TRMLCoolingSource;
    user_dt = true;
    user_time_step_func = TRMLCoolingTimeStep;
  }

  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &u0 = pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const TRMLFrameTrackingData data = trml_data;

  par_for("pgen_trml_frame_tracking", DevExeSpace(), 0, pmbp->nmb_thispack-1,
  ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const RegionSize block_size = size.d_view(m);
    const Real x1 = CellCenterX(i - is, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, block_size.x2min, block_size.x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, block_size.x3min, block_size.x3max);
    Real density, pressure, vx, vy, vz, cold_fraction;
    EvaluateTRMLState(x1, x2, x3, data, density, pressure, vx, vy, vz,
                      cold_fraction);
    SetHydroState(u0, m, k, j, i, density, pressure, vx, vy, vz, data.gm1);
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      u0(m, n, k, j, i) = density*cold_fraction;
    }
  });
}
