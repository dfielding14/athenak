//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file simple_TRML_extended.cpp
//! \brief Extended experimental radiatively cooling turbulent mixing layer.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream> // cout
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "driver/driver.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "srcterms/frame_tracker.hpp"

Real glob_rho_cold;
Real glob_rho_hot;
Real glob_pres;
Real glob_t_cool_min;
Real glob_T_cutoff_over_T_cold;
Real glob_T_ci_over_T_cold;
Real glob_T_ih_over_T_cold;
Real glob_beta;
Real glob_custom_min_timestep;
Real glob_velocity;
Real glob_fixed_frame_velocity_x3;
Real glob_cooling_ramp_time;
Real glob_cooling_ramp_min_factor;
Real glob_cooling_ramp_factor;
bool glob_zero_gradient_vx;
int glob_lower_x3_user_bc;
int glob_lower_x3_mass_balance_pressure_mode;
int glob_lower_x3_mass_balance_tangential_mode;
int glob_lower_x3_mass_balance_velocity_frame;
Real glob_lower_x3_mass_balance_vz;
Real glob_lower_x3_mass_balance_max_abs_vz;
Real glob_shear_vel_thresh;
Real glob_vy_vel_thresh;

constexpr int kLowerX3UserBCReflect = 0;
constexpr int kLowerX3UserBCOutflow = 1;
constexpr int kLowerX3UserBCMassBalance = 2;
constexpr int kMassBalancePressureFixed = 0;
constexpr int kMassBalancePressureTotal = 1;
constexpr int kMassBalanceTangentialZeroGradient = 0;
constexpr int kMassBalanceTangentialReservoir = 1;
constexpr int kMassBalanceVelocityGrid = 0;
constexpr int kMassBalanceVelocityLab = 1;

// These are rank-local history diagnostics. AthenaK performs the global sum
// when the history file is written, so CoolingSrc does not need an MPI
// collective at every stage.
Real glob_cooling_rate;
Real glob_cooling_rk_u0;
Real glob_cooling_rk_u1;

void CoolingSrc(Mesh *pm, Real bdt, Driver *pdrive, int stage);
void CoolingTimestep(Mesh *pm);
void TRMLZBoundary(Mesh *pm);
void HistoryOutput(HistoryData *pdata, Mesh *pm);

void FatalSimpleTRMLInput(const std::string &message) {
  std::cout << "### FATAL ERROR in simple_TRML input" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real shiftedtanh(Real x, Real low, Real high) {
  return (high - low) / 2 * tanh(x) + (high + low) / 2;
}

Real CoolingRampFactorAtTime(const Real time) {
  if (glob_cooling_ramp_time <= 0.0) {
    return 1.0;
  }
  Real s = time/glob_cooling_ramp_time;
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;
  const Real smoothstep = s*s*(3.0 - 2.0*s);
  return glob_cooling_ramp_min_factor +
         (1.0 - glob_cooling_ramp_min_factor)*smoothstep;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // Grab meshblock pointers
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;

  // Indices
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;

  glob_rho_cold = pin->GetReal("problem", "rho_cold");
  glob_rho_hot = pin->GetReal("problem", "rho_hot");
  glob_pres = pin->GetReal("problem", "pres");
  glob_t_cool_min = pin->GetReal("problem", "t_cool_min");
  glob_T_cutoff_over_T_cold = pin->GetReal("problem", "T_cutoff_over_T_cold");
  glob_T_ci_over_T_cold = pin->GetReal("problem", "T_ci_over_T_cold");
  glob_T_ih_over_T_cold = pin->GetReal("problem", "T_ih_over_T_cold");
  glob_beta = pin->GetReal("problem", "beta");
  glob_velocity = pin->GetReal("problem", "velocity");
  if (pin->DoesParameterExist("problem", "initial_profile")) {
    const std::string initial_profile_name =
        pin->GetString("problem", "initial_profile");
    if (initial_profile_name != "tanh" &&
        initial_profile_name != "smooth_tanh") {
      FatalSimpleTRMLInput(
          "<problem>/initial_profile='" + initial_profile_name +
          "' is no longer supported by simple_TRML. This pgen now only uses "
          "the smooth tanh TRML IC; put 1D-front/table ICs in a separate pgen.");
    }
  }
  if (pin->DoesParameterExist("problem", "front_j") ||
      pin->DoesParameterExist("problem", "front_step_x3")) {
    FatalSimpleTRMLInput(
        "<problem>/front_j and <problem>/front_step_x3 have been removed from "
        "simple_TRML. Use <problem>/fixed_frame_velocity_x3 and, if needed, "
        "<problem>/lower_x3_mass_balance_vz instead.");
  }
  glob_fixed_frame_velocity_x3 =
      pin->GetOrAddReal("problem", "fixed_frame_velocity_x3", 0.0);
  glob_cooling_ramp_time =
      pin->GetOrAddReal("problem", "cooling_ramp_time", -1.0);
  glob_cooling_ramp_min_factor =
      pin->GetOrAddReal("problem", "cooling_ramp_min_factor", 0.0);
  if (glob_cooling_ramp_min_factor < 0.0 ||
      glob_cooling_ramp_min_factor > 1.0) {
    FatalSimpleTRMLInput(
        "Require 0 <= <problem>/cooling_ramp_min_factor <= 1.");
  }
  glob_zero_gradient_vx =
      pin->GetOrAddBoolean("problem", "zero_gradient_vx", true);
  std::string lower_x3_user_bc =
      pin->GetOrAddString("problem", "lower_x3_user_bc", "reflect");
  if (lower_x3_user_bc == "reflect" || lower_x3_user_bc == "reflecting" ||
      lower_x3_user_bc == "cold_reflect") {
    glob_lower_x3_user_bc = kLowerX3UserBCReflect;
  } else if (lower_x3_user_bc == "outflow") {
    glob_lower_x3_user_bc = kLowerX3UserBCOutflow;
  } else if (lower_x3_user_bc == "mass_balance" ||
             lower_x3_user_bc == "mass-balanced" ||
             lower_x3_user_bc == "balanced_outflow") {
    glob_lower_x3_user_bc = kLowerX3UserBCMassBalance;
  } else {
    FatalSimpleTRMLInput("Invalid <problem>/lower_x3_user_bc='" +
                         lower_x3_user_bc +
                         "'; expected 'reflect', 'outflow', or 'mass_balance'.");
  }
  std::string mass_balance_pressure_mode =
      pin->GetOrAddString("problem", "lower_x3_mass_balance_pressure_mode",
                          "fixed");
  if (mass_balance_pressure_mode == "fixed" ||
      mass_balance_pressure_mode == "initial" ||
      mass_balance_pressure_mode == "thermal") {
    glob_lower_x3_mass_balance_pressure_mode = kMassBalancePressureFixed;
  } else if (mass_balance_pressure_mode == "total_pressure" ||
             mass_balance_pressure_mode == "ram_pressure" ||
             mass_balance_pressure_mode == "total") {
    glob_lower_x3_mass_balance_pressure_mode = kMassBalancePressureTotal;
  } else {
    FatalSimpleTRMLInput("Invalid <problem>/lower_x3_mass_balance_pressure_mode='" +
                         mass_balance_pressure_mode +
                         "'; expected 'fixed' or 'total_pressure'.");
  }
  std::string mass_balance_tangential_mode =
      pin->GetOrAddString("problem", "lower_x3_mass_balance_tangential_mode",
                          "zero_gradient");
  if (mass_balance_tangential_mode == "zero_gradient" ||
      mass_balance_tangential_mode == "copy") {
    glob_lower_x3_mass_balance_tangential_mode =
        kMassBalanceTangentialZeroGradient;
  } else if (mass_balance_tangential_mode == "reservoir" ||
             mass_balance_tangential_mode == "fixed") {
    glob_lower_x3_mass_balance_tangential_mode =
        kMassBalanceTangentialReservoir;
  } else {
    FatalSimpleTRMLInput("Invalid <problem>/lower_x3_mass_balance_tangential_mode='" +
                         mass_balance_tangential_mode +
                         "'; expected 'zero_gradient' or 'reservoir'.");
  }
  glob_lower_x3_mass_balance_vz =
      pin->GetOrAddReal("problem", "lower_x3_mass_balance_vz", 0.0);
  glob_lower_x3_mass_balance_max_abs_vz =
      pin->GetOrAddReal("problem", "lower_x3_mass_balance_max_abs_vz", -1.0);
  std::string mass_balance_velocity_frame =
      pin->GetOrAddString("problem",
                          "lower_x3_mass_balance_velocity_frame",
                          "grid");
  if (mass_balance_velocity_frame == "grid" ||
      mass_balance_velocity_frame == "frame" ||
      mass_balance_velocity_frame == "moving_frame" ||
      mass_balance_velocity_frame == "comoving") {
    glob_lower_x3_mass_balance_velocity_frame = kMassBalanceVelocityGrid;
  } else if (mass_balance_velocity_frame == "lab" ||
             mass_balance_velocity_frame == "laboratory" ||
             mass_balance_velocity_frame == "inertial" ||
             mass_balance_velocity_frame == "physical") {
    glob_lower_x3_mass_balance_velocity_frame = kMassBalanceVelocityLab;
  } else {
    FatalSimpleTRMLInput(
        "Invalid <problem>/lower_x3_mass_balance_velocity_frame='" +
        mass_balance_velocity_frame + "'; expected 'grid' or 'lab'.");
  }

  glob_shear_vel_thresh =
      pin->GetOrAddReal("problem", "hist_shear_vel_frac", 0.45) * glob_velocity;
  glob_vy_vel_thresh =
      pin->GetOrAddReal("problem", "hist_vy_vel_frac", 0.08) * glob_velocity;

  glob_custom_min_timestep =
      pin->GetOrAddReal("problem", "custom_min_timestep", -1.0);

  // The history value is zero until the first complete timestep. The two RK
  // registers are reset again at stage 1 of every timestep.
  glob_cooling_rate = 0.0;
  glob_cooling_rk_u0 = 0.0;
  glob_cooling_rk_u1 = 0.0;
  glob_cooling_ramp_factor = CoolingRampFactorAtTime(pmy_mesh_->time);

  Real rho_cold = glob_rho_cold;
  Real rho_hot = glob_rho_hot;
  Real pres = glob_pres;
  Real t_cool_min = glob_t_cool_min;
  Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
  Real beta = glob_beta;
  Real sharpness = pin->GetReal("problem", "phase_sharpness");
  Real velocity = glob_velocity;
  Real fixed_frame_velocity_x3 = glob_fixed_frame_velocity_x3;

  Real x1min = pin->GetReal("mesh", "x1min");
  Real x1max = pin->GetReal("mesh", "x1max");
  Real x2min = pin->GetReal("mesh", "x2min");
  Real x2max = pin->GetReal("mesh", "x2max");
  Real x3min = pin->GetReal("mesh", "x3min");
  Real x3max = pin->GetReal("mesh", "x3max");

  uint max_init_p = pin->GetInteger("problem", "max_init_perturb_log2freq");
  uint min_init_p = pin->GetInteger("problem", "min_init_perturb_log2freq");
  Real init_p_sharp = pin->GetReal("problem", "init_perturb_sharpness");
  Real init_p_vel_frac = pin->GetReal("problem", "init_perturb_vel_frac");

  // Use the stage-aware callback so the cooling history follows the selected RK
  // method. Existing pgens can continue using user_srcs_func without receiving
  // stage data.
  user_stage_srcs_func = CoolingSrc;
  user_time_step_func = CoolingTimestep;
  user_bcs_func = TRMLZBoundary;
  user_hist_func = HistoryOutput;

  // for restarting
  if (restart)
    return;

  if (pmbp->phydro != nullptr) {
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    int nhydro = pmbp->phydro->nhydro;
    int nscalars = pmbp->phydro->nscalars;

    // Set initial conditions
    par_for(
        "TRML_problem", DevExeSpace(), 0, (pmbp->nmb_thispack - 1), ks, ke, js,
        je, is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;

          Real coordx = CellCenterX(i - is, indcs.nx1, x1min, x1max);
          Real coordy = CellCenterX(j - js, indcs.nx2, x2min, x2max);
          Real coordz = CellCenterX(k - ks, indcs.nx3, x3min, x3max);
          Real dens = shiftedtanh(coordz * sharpness, rho_cold, rho_hot);
          Real cold_fraction = shiftedtanh(coordz * sharpness, 1.0, 0.0);
          Real vx = shiftedtanh(coordz * sharpness, -0.5 * velocity,
                                0.5 * velocity);
          // The original lab-frame IC has no vertical flow.  A positive
          // fixed_frame_velocity_x3 means the grid moves in +x3, so stored
          // velocities use v_grid = v_lab - V_frame.
          Real vz = -fixed_frame_velocity_x3;
          u0(m, IDN, k, j, i) = dens;
          u0(m, IM1, k, j, i) = dens*vx;
          u0(m, IM2, k, j, i) = 0.0;
          u0(m, IM3, k, j, i) = dens*vz;

          for (int i0 = min_init_p; i0 <= max_init_p; i0++) {
            for (int j0 = min_init_p; j0 <= max_init_p; j0++) {
              u0(m, IM3, k, j, i) +=
                  velocity * init_p_vel_frac *
                  sin(2 * 3.141592653589 * coordy * pow(2, j0)) *
                  sin(2 * 3.141592653589 * coordx * pow(2, i0)) *
                  u0(m, IDN, k, j, i) *
                  exp(-coordz * coordz * init_p_sharp / 2);
            }
          }

          if (eos.is_ideal) {
            u0(m, IEN, k, j, i) = pres / gm1 + 0.5 *
                                                   (SQR(u0(m, IM1, k, j, i)) +
                                                    SQR(u0(m, IM2, k, j, i)) +
                                                    SQR(u0(m, IM3, k, j, i))) /
                                                   u0(m, IDN, k, j, i);
          }
          for (int n = nhydro; n < nhydro + nscalars; ++n) {
            u0(m, n, k, j, i) = dens * cold_fraction;
          }
        });
  }

  return;
}

void CoolingSrc(Mesh *pm, Real bdt, Driver *pdrive, int stage) {
  MeshBlockPack *pmbp = pm->pmb_pack;

  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack) * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;

  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  Real rho_cold = glob_rho_cold;
  Real pres = glob_pres;
  Real t_cool_min = glob_t_cool_min;
  Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
  Real beta = glob_beta;

  Real T_cold = pres / rho_cold;
  Real cooling_ramp_factor = CoolingRampFactorAtTime(pm->time);
  glob_cooling_ramp_factor = cooling_ramp_factor;
  Real cooling_bdt = bdt*cooling_ramp_factor;
  // Sum the positive amount of energy removed on this rank by this source call.
  // deltaE already contains the stage's beta*dt, so it must not be multiplied
  // by beta again.
  Real cooling_energy_this_stage;
  pdrive->StartPerformanceRegion(PerformanceRegion::trml_cooling);
  Kokkos::parallel_reduce(
      "TRML_cooling_src", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &stage_cooling) {
        int m = (idx) / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        int i = (idx - m * nkji - k * nji - j * nx1) + is;
        k += ks;
        j += js;

        Real dens = w0(m, IDN, k, j, i);
        // For ideal hydro, w0(IEN) is already primitive internal-energy
        // density.
        Real eint = w0(m, IEN, k, j, i);
        Real temp = eint / dens * gm1;

        Real temp_at_t_cool_min =
            (beta > 0) ? T_cold : (T_cutoff_over_T_cold * T_cold);
        // This is the exact solution for the cooling
        Real deltaE;
        if (beta != 0.0) {
          Real partofcooling = 1 - pow(temp / temp_at_t_cool_min, -beta) *
                                   cooling_bdt / t_cool_min * beta;
          if (partofcooling > 0)
            deltaE = eint * pow(partofcooling, 1 / beta) - eint;
          else
            deltaE = dens * T_cold / gm1 - eint;
        } else {
          deltaE = eint * exp(-cooling_bdt / t_cool_min) - eint;
        }

        const Real eint_floor = dens*T_cold/gm1;
        if ((temp / T_cold) <= 1.0 || (temp / T_cold) > T_cutoff_over_T_cold) {
          deltaE = 0.0;
        } else {
          deltaE = fmax(deltaE, eint_floor - eint);
        }

        Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
        stage_cooling += -deltaE * vol;

        u0(m, IEN, k, j, i) += deltaE;
      },
      Kokkos::Sum<Real>(cooling_energy_this_stage));
  pdrive->StopPerformanceRegion(PerformanceRegion::trml_cooling);

  // Mirror Hydro::CopyCons for the scalar cooling diagnostic. u1 stores the
  // initial state for RK1/RK2/RK3, while AthenaK's low-storage RK4 also updates
  // u1 at later stages with the driver's delta coefficient.
  if (stage == 1) {
    glob_cooling_rk_u0 = 0.0;
    glob_cooling_rk_u1 = glob_cooling_rk_u0;
  } else if (pdrive->integrator == "rk4") {
    glob_cooling_rk_u1 += pdrive->delta[stage - 1] * glob_cooling_rk_u0;
  }

  // Mirror Hydro::RKUpdate, then add the cooling from this source call exactly
  // where the hydro energy update adds deltaE. This automatically gives the
  // correct retained cooling for RK1, RK2, RK3, and the driver's two-register
  // RK4 implementation.
  glob_cooling_rk_u0 = pdrive->gam0[stage - 1] * glob_cooling_rk_u0 +
                       pdrive->gam1[stage - 1] * glob_cooling_rk_u1 +
                       cooling_energy_this_stage;

  // Publish only a completed-timestep value. This remains rank-local until
  // HistoryOutput returns it on every rank and AthenaK's history writer
  // performs its normal MPI sum.
  if (stage == pdrive->nexp_stages) {
    glob_cooling_rate = glob_cooling_rk_u0 / pm->dt;
  }
}

void CoolingTimestep(Mesh *pm) {
  if (glob_custom_min_timestep > 0.0)
    pm->pgen->dtnew = glob_custom_min_timestep;
  else
    pm->pgen->dtnew = glob_t_cool_min;
}

void TRMLZBoundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &mb_bcs = pmbp->pmb->mb_bcs;

  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;

  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  auto gm1 = eos.gamma - 1.0;
  int nhydro = pmbp->phydro->nhydro;
  int nscalars = pmbp->phydro->nscalars;

  Real rho_cold = glob_rho_cold;
  Real rho_hot = glob_rho_hot;
  Real velocity = glob_velocity;
  Real pres = glob_pres;
  Real fixed_frame_velocity_x3 = glob_fixed_frame_velocity_x3;
  bool zero_gradient_vx = glob_zero_gradient_vx;
  int lower_x3_user_bc = glob_lower_x3_user_bc;
  int mass_balance_pressure_mode = glob_lower_x3_mass_balance_pressure_mode;
  int mass_balance_tangential_mode = glob_lower_x3_mass_balance_tangential_mode;
  int mass_balance_velocity_frame = glob_lower_x3_mass_balance_velocity_frame;
  Real lower_x3_mass_balance_vz = glob_lower_x3_mass_balance_vz;

  // Frame tracking is restricted to x3.  The x1 reservoir shear velocities are
  // therefore unchanged; the optional zero-gradient vx condition instead copies
  // the adjacent velocity.
  bool frame_tracking = (pmbp->pframe_tracker != nullptr);
  Real frame_v2 = frame_tracking ? pmbp->pframe_tracker->FrameVelocity(1) : 0.0;
  Real frame_v3 = fixed_frame_velocity_x3 +
                  (frame_tracking ? pmbp->pframe_tracker->FrameVelocity(2) : 0.0);
  Real lower_mass_balance_vz_grid = 0.0;
  Real lower_mass_balance_vz_lab = frame_v3;
  Real mass_balance_top_mdot_average = 0.0;
  Real mass_balance_area = 0.0;
  bool have_live_top_mdot_average = false;
  if (lower_x3_user_bc == kLowerX3UserBCMassBalance) {
    if (frame_tracking &&
        pmbp->pframe_tracker->TopMdotAverageWeightTime() > 0.0) {
      mass_balance_top_mdot_average = pmbp->pframe_tracker->LastTopMdotAverage();
      mass_balance_area = pmbp->pframe_tracker->TopMdotArea();
      have_live_top_mdot_average = true;
    }
    if (mass_balance_area <= 0.0) {
      const RegionSize &mesh_size = pm->mesh_size;
      mass_balance_area = (mesh_size.x1max - mesh_size.x1min)*
                          (mesh_size.x2max - mesh_size.x2min);
    }
    if (rho_cold > 0.0 && mass_balance_area > 0.0) {
      // Mdot_top = integral rho*v3*dA at upper x3, positive outward.
      // Mdot_bot = -integral rho*v3*dA at lower x3, positive outward.
      // Setting v3_bot = Mdot_top/(rho_cold*A) gives Mdot_bot = -Mdot_top,
      // measured in the velocity frame selected by
      // lower_x3_mass_balance_velocity_frame.
      Real lower_mass_balance_vz_target = lower_x3_mass_balance_vz;
      if (have_live_top_mdot_average) {
        lower_mass_balance_vz_target =
            mass_balance_top_mdot_average/(rho_cold*mass_balance_area);
      }
      if (glob_lower_x3_mass_balance_max_abs_vz > 0.0 &&
          std::fabs(lower_mass_balance_vz_target) >
          glob_lower_x3_mass_balance_max_abs_vz) {
        lower_mass_balance_vz_target =
            std::copysign(glob_lower_x3_mass_balance_max_abs_vz,
                          lower_mass_balance_vz_target);
      }
      if (mass_balance_velocity_frame == kMassBalanceVelocityLab) {
        lower_mass_balance_vz_lab = lower_mass_balance_vz_target;
        lower_mass_balance_vz_grid = lower_mass_balance_vz_lab - frame_v3;
      } else {
        lower_mass_balance_vz_grid = lower_mass_balance_vz_target;
        lower_mass_balance_vz_lab = lower_mass_balance_vz_grid + frame_v3;
      }
    }
  }
  Real lower_mass_balance_pressure = pres;
  if (mass_balance_pressure_mode == kMassBalancePressureTotal) {
    // Thermal pressure support for the cold reservoir against the hot inflow's
    // ram pressure in the cold-gas frame. Stored velocities are grid-frame
    // velocities, so v_lab = v_grid + v_frame.  When the mass-balance target
    // is declared in the lab frame, the top-Mdot-derived hot velocity is
    // already lab-frame; otherwise it is a grid-frame velocity and must be
    // transformed before taking the physical relative velocity.
    Real hot_flux_vz = 0.0;
    if (rho_hot > 0.0 && mass_balance_area > 0.0) {
      hot_flux_vz = mass_balance_top_mdot_average/(rho_hot*mass_balance_area);
    }
    const Real hot_lab_vz = have_live_top_mdot_average ?
        ((mass_balance_velocity_frame == kMassBalanceVelocityLab) ?
         hot_flux_vz : hot_flux_vz + frame_v3) : 0.0;
    const Real relative_hot_cold_vz = hot_lab_vz - lower_mass_balance_vz_lab;
    lower_mass_balance_pressure =
        pres + rho_hot*SQR(relative_hot_cold_vz);
    lower_mass_balance_pressure = std::max(lower_mass_balance_pressure,
                                           eos.pfloor);
  }

  par_for(
      "TRML_boundary", DevExeSpace(), 0, (pmbp->nmb_thispack - 1), js, je, is,
      ie, KOKKOS_LAMBDA(int m, int j, int i) {
        if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
          for (int k = 0; k < ng; k++) {
            int ghost_inner_k = ks - k - 1;
            if (lower_x3_user_bc == kLowerX3UserBCOutflow) {
              for (int n = 0; n < nhydro + nscalars; ++n) {
                u0(m, n, ghost_inner_k, j, i) = u0(m, n, ks, j, i);
              }
              continue;
            }
            if (lower_x3_user_bc == kLowerX3UserBCMassBalance) {
              u0(m, IDN, ghost_inner_k, j, i) = rho_cold;
              if (mass_balance_tangential_mode ==
                  kMassBalanceTangentialZeroGradient) {
                const Real active_rho = u0(m, IDN, ks, j, i);
                const Real inv_active_rho =
                    (active_rho > 0.0) ? 1.0/active_rho : 0.0;
                u0(m, IM1, ghost_inner_k, j, i) =
                    rho_cold*u0(m, IM1, ks, j, i)*inv_active_rho;
                u0(m, IM2, ghost_inner_k, j, i) =
                    rho_cold*u0(m, IM2, ks, j, i)*inv_active_rho;
              } else {
                u0(m, IM1, ghost_inner_k, j, i) =
                    rho_cold*(-0.5*velocity);
                u0(m, IM2, ghost_inner_k, j, i) =
                    frame_tracking ? -rho_cold*frame_v2 : 0.0;
              }
              u0(m, IM3, ghost_inner_k, j, i) =
                  rho_cold*lower_mass_balance_vz_grid;
              u0(m, IEN, ghost_inner_k, j, i) =
                  lower_mass_balance_pressure/gm1 +
                  0.5*(SQR(u0(m, IM1, ghost_inner_k, j, i)) +
                       SQR(u0(m, IM2, ghost_inner_k, j, i)) +
                       SQR(u0(m, IM3, ghost_inner_k, j, i))) /
                  u0(m, IDN, ghost_inner_k, j, i);
              for (int n = nhydro; n < nhydro + nscalars; ++n) {
                u0(m, n, ghost_inner_k, j, i) = rho_cold;
              }
              continue;
            }
            u0(m, IDN, ghost_inner_k, j, i) = rho_cold;
            if (zero_gradient_vx) {
              u0(m, IM1, ghost_inner_k, j, i) =
                  rho_cold * u0(m, IM1, ks, j, i) / u0(m, IDN, ks, j, i);
            } else {
              u0(m, IM1, ghost_inner_k, j, i) = rho_cold * (-0.5 * velocity);
            }
            if (frame_tracking) {
              u0(m, IM2, ghost_inner_k, j, i) = -rho_cold * frame_v2;
            } else {
              u0(m, IM2, ghost_inner_k, j, i) =
                  rho_cold * u0(m, IM2, ks, j, i) / u0(m, IDN, ks, j, i);
            }
            // Fixed cold reservoir wall: reflect the normal velocity in the
            // stored/grid frame.  This makes the wall stationary in the
            // computational frame whether the frame velocity is zero, fixed,
            // or adaptively tracked.
            u0(m, IM3, ghost_inner_k, j, i) =
                -rho_cold * u0(m, IM3, ks + k, j, i) / u0(m, IDN, ks + k, j, i);
            u0(m, IEN, ghost_inner_k, j, i) =
                pres / gm1 + 0.5 *
                                 (SQR(u0(m, IM1, ghost_inner_k, j, i)) +
                                  SQR(u0(m, IM2, ghost_inner_k, j, i)) +
                                  SQR(u0(m, IM3, ghost_inner_k, j, i))) /
                                 u0(m, IDN, ghost_inner_k, j, i);
            for (int n = nhydro; n < nhydro + nscalars; ++n) {
              u0(m, n, ghost_inner_k, j, i) = rho_cold;
            }
          }
        }

        if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user) {
          for (int k = 0; k < ng; k++) {
            int ghost_outer_k = ke + k + 1;
            u0(m, IDN, ghost_outer_k, j, i) = rho_hot;
            // Hot reservoir with zero-gradient transverse/normal flow:
            // enforce rho, P, and the x1 shear velocity, while copying v_y and
            // v_z from the adjacent active zone. Frame tracking is x3-only, so
            // it does not transform vx.
            u0(m, IM1, ghost_outer_k, j, i) =
                rho_hot * (0.5 * velocity);
            const Real active_rho = u0(m, IDN, ke, j, i);
            const Real inv_active_rho =
                (active_rho > 0.0) ? 1.0/active_rho : 0.0;
            u0(m, IM2, ghost_outer_k, j, i) =
                rho_hot * u0(m, IM2, ke, j, i) * inv_active_rho;
            u0(m, IM3, ghost_outer_k, j, i) =
                rho_hot * u0(m, IM3, ke, j, i) * inv_active_rho;
            u0(m, IEN, ghost_outer_k, j, i) =
                pres / gm1 + 0.5 *
                                 (SQR(u0(m, IM1, ghost_outer_k, j, i)) +
                                  SQR(u0(m, IM2, ghost_outer_k, j, i)) +
                                  SQR(u0(m, IM3, ghost_outer_k, j, i))) /
                                 u0(m, IDN, ghost_outer_k, j, i);
            for (int n = nhydro; n < nhydro + nscalars; ++n) {
              u0(m, n, ghost_outer_k, j, i) = 0.0;
            }
          }
        }
      });
}

#define TRMLHISTVARS 29

namespace Kokkos { // reduction identity must be defined in Kokkos namespace
template <>
struct reduction_identity<array_sum::array_type<Real, TRMLHISTVARS>> {
  KOKKOS_FORCEINLINE_FUNCTION static array_sum::array_type<Real, TRMLHISTVARS>
  sum() {
    return array_sum::array_type<Real, TRMLHISTVARS>();
  }
};
} // namespace Kokkos

void HistoryOutput(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack) * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;
  Real lxmax3 = pm->mesh_size.x3max;
  Real lxmin3 = pm->mesh_size.x3min;

  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  auto gm1 = eos.gamma - 1.0;

  Real velocity = glob_velocity;
  Real pres = glob_pres;
  Real rho_cold = glob_rho_cold;
  Real T_cold = pres / rho_cold;
  Real T_ci_over_T_cold = glob_T_ci_over_T_cold;
  Real T_ih_over_T_cold = glob_T_ih_over_T_cold;
  Real shear_vel_thresh = glob_shear_vel_thresh;
  Real vy_vel_thresh = glob_vy_vel_thresh;
  bool frame_tracking = (pmbp->pframe_tracker != nullptr);
  Real frame_v2_hist = frame_tracking ? pmbp->pframe_tracker->FrameVelocity(1) : 0.0;
  Real frame_v3_hist = glob_fixed_frame_velocity_x3 +
                       (frame_tracking ? pmbp->pframe_tracker->FrameVelocity(2) : 0.0);

  pdata->nhist = TRMLHISTVARS;
  int sum_min_index = 1;
  int sum_max_index = 17;
  pdata->label[0] = "cooling_rate ";

  pdata->label[1] = "M_flux_top "; // Mass flux out of the top
  pdata->label[2] = "M_flux_bot "; // Mass flux out of the bottom

  pdata->label[3] = "E_flux_top "; // Energy flux out of the top
  pdata->label[4] = "E_flux_bot "; // Energy flux out of the bottom

  pdata->label[5] = "vysq_c "; // vy^2 for cold gas
  pdata->label[6] = "vysq_i "; //         intermediate gas
  pdata->label[7] = "vysq_h "; //         hot gas

  pdata->label[8] = "kex_c ";  // x vel kinetic energy for cold gas
  pdata->label[9] = "kex_i ";  //                         intermediate gas
  pdata->label[10] = "kex_h "; //                        hot gas

  pdata->label[11] = "momx_c "; // x momentum for cold gas
  pdata->label[12] = "momx_i "; //               intermediate gas
  pdata->label[13] = "momx_h "; //               hot gas

  pdata->label[14] = "V_c "; // Volume for cold gas
  pdata->label[15] = "V_i "; //           intermediate gas
  pdata->label[16] = "V_h "; //           hot gas

  pdata->label[17] = "zavg_i "; // Average z pos for intermediate gas

  pdata->label[18] = "zmin_shear ";
  pdata->label[19] = "zmax_shear ";
  pdata->label[20] = "zmin_intermediate ";
  pdata->label[21] = "zmax_intermediate ";
  pdata->label[22] = "zmin_vy ";
  pdata->label[23] = "zmax_vy ";
  pdata->label[24] = "Mtop_lab ";
  pdata->label[25] = "Mbot_lab ";
  pdata->label[26] = "Etop_lab ";
  pdata->label[27] = "Ebot_lab ";
  pdata->label[28] = "cool_ramp ";

  array_sum::array_type<Real, TRMLHISTVARS> history_sum;
  Real min_z_intermediate = std::numeric_limits<Real>::max();
  Real max_z_intermediate = -std::numeric_limits<Real>::max();

  Real min_z_shear = std::numeric_limits<Real>::max();
  Real max_z_shear = -std::numeric_limits<Real>::max();

  Real min_z_vely = std::numeric_limits<Real>::max();
  Real max_z_vely = -std::numeric_limits<Real>::max();
  int cnt_z_intermediate = 0;
  int cnt_z_shear = 0;
  int cnt_z_vely = 0;

  Kokkos::parallel_reduce(
      "TRML_Hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx,
                    array_sum::array_type<Real, TRMLHISTVARS> &mb_sum,
                    Real &min_z_intermediate_, Real &max_z_intermediate_,
                    Real &min_z_shear_, Real &max_z_shear_, Real &min_z_vely_,
                    Real &max_z_vely_, int &cnt_z_intermediate_,
                    int &cnt_z_shear_, int &cnt_z_vely_) {
        int m = (idx) / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        int i = (idx - m * nkji - k * nji - j * nx1) + is;
        k += ks;
        j += js;

        Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);
        Real dz = size.d_view(m).dx3;
        Real dA = size.d_view(m).dx1 * size.d_view(m).dx2;

        Real dens = w0(m, IDN, k, j, i);
        Real velx = w0(m, IVX, k, j, i);
        Real vely = w0(m, IVY, k, j, i);
        Real velz = w0(m, IVZ, k, j, i);
        Real velx_lab = velx;
        Real vely_lab = vely + frame_v2_hist;
        Real velz_lab = velz + frame_v3_hist;

        Real temp = 0.0;
        Real eint = 0.0;
        temp = w0(m, IEN, k, j, i) / dens * gm1;
        eint = w0(m, IEN, k, j, i);
        Real ekin = 0.5 * w0(m, IDN, k, j, i) *
                    (SQR(w0(m, IVX, k, j, i)) + SQR(w0(m, IVY, k, j, i)) +
                     SQR(w0(m, IVZ, k, j, i)));
        Real etot = eint + ekin;
        Real ekin_lab = 0.5*dens*
                        (SQR(velx_lab) + SQR(vely_lab) + SQR(velz_lab));
        Real etot_lab = eint + ekin_lab;
        Real work = (gm1)*eint;

        if (fabs(x3v - lxmax3) < dz) {
          mb_sum.the_array[1] += dens * dA * velz;
          mb_sum.the_array[3] += (etot + work) * dA * velz;
          mb_sum.the_array[24] += dens * dA * velz_lab;
          mb_sum.the_array[26] += (etot_lab + work) * dA * velz_lab;
        }

        if (fabs(x3v - lxmin3) < dz) {
          mb_sum.the_array[2] -= dens * dA * velz;
          mb_sum.the_array[4] -= (etot + work) * dA * velz;
          mb_sum.the_array[25] -= dens * dA * velz_lab;
          mb_sum.the_array[27] -= (etot_lab + work) * dA * velz_lab;
        }

        int temp_bin;
        if (temp / T_cold < T_ci_over_T_cold)
          temp_bin = 0;
        else if (temp / T_cold >= T_ci_over_T_cold &&
                 temp / T_cold < T_ih_over_T_cold)
          temp_bin = 1;
        else
          temp_bin = 2;

        // Split data into temp bins
        mb_sum.the_array[5 + temp_bin] += SQR(vely) * vol;
        mb_sum.the_array[8 + temp_bin] += 0.5 * dens * SQR(velx) * vol;
        mb_sum.the_array[11 + temp_bin] += dens * velx * vol;
        mb_sum.the_array[14 + temp_bin] += vol;
        if (temp_bin == 1) {
          mb_sum.the_array[17] += x3v * vol;
          min_z_intermediate_ = fmin(x3v, min_z_intermediate_);
          max_z_intermediate_ = fmax(x3v, max_z_intermediate_);
          cnt_z_intermediate_ += 1;
        }

        if (fabs(velx) < shear_vel_thresh) {
          min_z_shear_ = fmin(x3v, min_z_shear_);
          max_z_shear_ = fmax(x3v, max_z_shear_);
          cnt_z_shear_ += 1;
        }

        if (fabs(vely) > vy_vel_thresh) {
          min_z_vely_ = fmin(x3v, min_z_vely_);
          max_z_vely_ = fmax(x3v, max_z_vely_);
          cnt_z_vely_ += 1;
        }
      },
      Kokkos::Sum<array_sum::array_type<Real, TRMLHISTVARS>>(history_sum),
      Kokkos::Min<Real>(min_z_intermediate),
      Kokkos::Max<Real>(max_z_intermediate), Kokkos::Min<Real>(min_z_shear),
      Kokkos::Max<Real>(max_z_shear), Kokkos::Min<Real>(min_z_vely),
      Kokkos::Max<Real>(max_z_vely), Kokkos::Sum<int>(cnt_z_intermediate),
      Kokkos::Sum<int>(cnt_z_shear), Kokkos::Sum<int>(cnt_z_vely));

#if MPI_PARALLEL_ENABLED
  Real m_min[3] = {min_z_intermediate, min_z_shear, min_z_vely};
  Real m_max[3] = {max_z_intermediate, max_z_shear, max_z_vely};
  Real gm_min[3];
  Real gm_max[3];
  int loc_counts[3] = {cnt_z_intermediate, cnt_z_shear, cnt_z_vely};
  int glob_counts[3] = {0, 0, 0};
  // MPI_Allreduce(MPI_IN_PLACE, &dtnew, 1, MPI_ATHENA_REAL, MPI_MIN,
  // MPI_COMM_WORLD);
  MPI_Allreduce(m_min, gm_min, 3, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(m_max, gm_max, 3, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(loc_counts, glob_counts, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  min_z_intermediate = gm_min[0];
  min_z_shear = gm_min[1];
  min_z_vely = gm_min[2];
  max_z_intermediate = gm_max[0];
  max_z_shear = gm_max[1];
  max_z_vely = gm_max[2];
  cnt_z_intermediate = glob_counts[0];
  cnt_z_shear = glob_counts[1];
  cnt_z_vely = glob_counts[2];
#endif
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  if (cnt_z_intermediate <= 0) {
    min_z_intermediate = nan;
    max_z_intermediate = nan;
  }
  if (cnt_z_shear <= 0) {
    min_z_shear = nan;
    max_z_shear = nan;
  }
  if (cnt_z_vely <= 0) {
    min_z_vely = nan;
    max_z_vely = nan;
  }

  // store data into hdata array
  for (int n = sum_min_index; n <= sum_max_index; n++) {
    pdata->hdata[n] = history_sum.the_array[n];
  }
  for (int n = 24; n <= 27; n++) {
    pdata->hdata[n] = history_sum.the_array[n];
  }
  // The cooling rate is rank-local here. Return it on every rank so AthenaK's
  // normal history MPI_Reduce produces the domain-wide rate without
  // synchronizing every stage.
  pdata->hdata[0] = glob_cooling_rate;

  // The history variables added undergo an MPI_Reduce before writing the
  // output. Since we have already reduced the following variables, we shall
  // only store them in the master processor. For quantities that are already
  // calculated in other parts of the code or which are maxima/minima, we store
  // their values in the master processer, since the current history function
  // can only perform a global sum
  if (global_variable::my_rank == 0) {
    pdata->hdata[18] = min_z_shear;
    pdata->hdata[19] = max_z_shear;
    pdata->hdata[20] = min_z_intermediate;
    pdata->hdata[21] = max_z_intermediate;
    pdata->hdata[22] = min_z_vely;
    pdata->hdata[23] = max_z_vely;
    pdata->hdata[28] = CoolingRampFactorAtTime(pm->time);
  } else {
    pdata->hdata[18] = 0.0;
    pdata->hdata[19] = 0.0;
    pdata->hdata[20] = 0.0;
    pdata->hdata[21] = 0.0;
    pdata->hdata[22] = 0.0;
    pdata->hdata[23] = 0.0;
    pdata->hdata[28] = 0.0;
  }
}
