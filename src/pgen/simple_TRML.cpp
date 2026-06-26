//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file simple_TRML.cpp
//! \brief Simple radiatively cooling turbulent mixing layer.

#include <cmath>
#include <iostream> // cout
#include <limits>

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
bool glob_zero_gradient_vx;
Real glob_shear_vel_thresh;
Real glob_vy_vel_thresh;

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

KOKKOS_INLINE_FUNCTION
Real shiftedtanh(Real x, Real low, Real high) {
  return (high - low) / 2 * tanh(x) + (high + low) / 2;
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
  glob_zero_gradient_vx =
      pin->GetOrAddBoolean("problem", "zero_gradient_vx", true);

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

  Real rho_cold = glob_rho_cold;
  Real rho_hot = glob_rho_hot;
  Real pres = glob_pres;
  Real t_cool_min = glob_t_cool_min;
  Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
  Real beta = glob_beta;
  Real sharpness = pin->GetReal("problem", "phase_sharpness");
  Real velocity = glob_velocity;

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
          u0(m, IDN, k, j, i) = dens;
          u0(m, IM1, k, j, i) =
              dens *
              shiftedtanh(coordz * sharpness, -0.5 * velocity, 0.5 * velocity);
          u0(m, IM2, k, j, i) = 0.0;
          u0(m, IM3, k, j, i) = 0.0;

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
  // Sum the positive amount of energy removed on this rank by this source call.
  // deltaE already contains the stage's beta*dt, so it must not be multiplied
  // by beta again.
  Real cooling_energy_this_stage;
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
          Real partofcooling = 1 - pow(temp / temp_at_t_cool_min, -beta) * bdt /
                                       t_cool_min * beta;
          if (partofcooling > 0)
            deltaE = eint * pow(partofcooling, 1 / beta) - eint;
          else
            deltaE = dens * T_cold / gm1 - eint;
        } else {
          deltaE = eint * exp(-bdt / t_cool_min) - eint;
        }

        if ((temp / T_cold) < 1 || (temp / T_cold) > T_cutoff_over_T_cold)
          deltaE = 0.0;

        Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
        stage_cooling += -deltaE * vol;

        u0(m, IEN, k, j, i) += deltaE;
      },
      Kokkos::Sum<Real>(cooling_energy_this_stage));

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
  bool zero_gradient_vx = glob_zero_gradient_vx;

  // The fixed reservoir states are uniform, so frame displacement does not change
  // them. Fixed velocities are transformed into the tracked grid frame; the optional
  // zero-gradient vx condition instead copies the adjacent grid-frame velocity.
  bool frame_tracking = (pmbp->pframe_tracker != nullptr);
  Real frame_v1 = frame_tracking ? pmbp->pframe_tracker->FrameVelocity(0) : 0.0;
  Real frame_v2 = frame_tracking ? pmbp->pframe_tracker->FrameVelocity(1) : 0.0;
  Real frame_v3 = frame_tracking ? pmbp->pframe_tracker->FrameVelocity(2) : 0.0;

  par_for(
      "TRML_boundary", DevExeSpace(), 0, (pmbp->nmb_thispack - 1), js, je, is,
      ie, KOKKOS_LAMBDA(int m, int j, int i) {
        if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
          for (int k = 0; k < ng; k++) {
            int ghost_inner_k = ks - k - 1;
            u0(m, IDN, ghost_inner_k, j, i) = rho_cold;
            if (zero_gradient_vx) {
              u0(m, IM1, ghost_inner_k, j, i) =
                  rho_cold * u0(m, IM1, ks, j, i) / u0(m, IDN, ks, j, i);
            } else if (frame_tracking) {
              u0(m, IM1, ghost_inner_k, j, i) =
                  rho_cold * (-0.5 * velocity - frame_v1);
            } else {
              u0(m, IM1, ghost_inner_k, j, i) = rho_cold * (-0.5 * velocity);
            }
            if (frame_tracking) {
              u0(m, IM2, ghost_inner_k, j, i) = -rho_cold * frame_v2;
              u0(m, IM3, ghost_inner_k, j, i) = -rho_cold * frame_v3;
            } else {
              u0(m, IM2, ghost_inner_k, j, i) =
                  rho_cold * u0(m, IM2, ks, j, i) / u0(m, IDN, ks, j, i);
              u0(m, IM3, ghost_inner_k, j, i) =
                  rho_cold * u0(m, IM3, ks, j, i) / u0(m, IDN, ks, j, i);
            }
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
            if (zero_gradient_vx) {
              u0(m, IM1, ghost_outer_k, j, i) =
                  rho_hot * u0(m, IM1, ke, j, i) / u0(m, IDN, ke, j, i);
            } else if (frame_tracking) {
              u0(m, IM1, ghost_outer_k, j, i) =
                  rho_hot * (0.5 * velocity - frame_v1);
            } else {
              u0(m, IM1, ghost_outer_k, j, i) = rho_hot * (0.5 * velocity);
            }
            if (frame_tracking) {
              u0(m, IM2, ghost_outer_k, j, i) = -rho_hot * frame_v2;
              u0(m, IM3, ghost_outer_k, j, i) = -rho_hot * frame_v3;
            } else {
              u0(m, IM2, ghost_outer_k, j, i) =
                  rho_hot * u0(m, IM2, ke, j, i) / u0(m, IDN, ke, j, i);
              u0(m, IM3, ghost_outer_k, j, i) =
                  rho_hot * u0(m, IM3, ke, j, i) / u0(m, IDN, ke, j, i);
            }
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

#define TRMLHISTVARS 24

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

        Real temp = 0.0;
        Real eint = 0.0;
        temp = w0(m, IEN, k, j, i) / dens * gm1;
        eint = w0(m, IEN, k, j, i);
        Real ekin = 0.5 * w0(m, IDN, k, j, i) *
                    (SQR(w0(m, IVX, k, j, i)) + SQR(w0(m, IVY, k, j, i)) +
                     SQR(w0(m, IVZ, k, j, i)));
        Real etot = eint + ekin;
        Real work = (gm1)*eint;

        if (fabs(x3v - lxmax3) < dz) {
          mb_sum.the_array[1] += dens * dA * velz;
          mb_sum.the_array[3] += (etot + work) * dA * velz;
        }

        if (fabs(x3v - lxmin3) < dz) {
          mb_sum.the_array[2] -= dens * dA * velz;
          mb_sum.the_array[4] -= (etot + work) * dA * velz;
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
  } else {
    pdata->hdata[18] = 0.0;
    pdata->hdata[19] = 0.0;
    pdata->hdata[20] = 0.0;
    pdata->hdata[21] = 0.0;
    pdata->hdata[22] = 0.0;
    pdata->hdata[23] = 0.0;
  }
}
