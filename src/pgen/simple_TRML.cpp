//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file simple_TRML.cpp
//! \brief Minimal radiatively cooling turbulent mixing layer.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "driver/driver.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "parameter_input.hpp"
#include "pgen.hpp"

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
Real glob_initial_vz;
Real glob_cooling_ramp_time;
Real glob_cooling_ramp_min_factor;
Real glob_shear_vel_thresh;
Real glob_vy_vel_thresh;

// Rank-local cooling-history registers. AthenaK performs the global sum when
// the history file is written.
Real glob_cooling_rate;
Real glob_cooling_rk_u0;
Real glob_cooling_rk_u1;

void CoolingSrc(Mesh *pm, Real bdt, Driver *pdrive, int stage);
void CoolingTimestep(Mesh *pm);
void TRMLOuterX3Boundary(Mesh *pm);
void HistoryOutput(HistoryData *pdata, Mesh *pm);

void FatalSimpleTRMLInput(const std::string &message) {
  std::cout << "### FATAL ERROR in simple_TRML input" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real shiftedtanh(Real x, Real low, Real high) {
  return (high - low)/2.0*tanh(x) + (high + low)/2.0;
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
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;

  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;

  // This pgen supplies only an outer-x3 user boundary. AthenaK still permits
  // ix3_bc=user syntactically, so reject it before leaving inner ghosts
  // unfilled.
  if (pin->GetString("mesh", "ix3_bc") == "user") {
    FatalSimpleTRMLInput(
        "<mesh>/ix3_bc=user is unsupported. Use an AthenaK native inner-x3 "
        "boundary such as 'outflow' or 'reflecting'.");
  }

  glob_rho_cold = pin->GetReal("problem", "rho_cold");
  glob_rho_hot = pin->GetReal("problem", "rho_hot");
  glob_pres = pin->GetReal("problem", "pres");
  glob_t_cool_min = pin->GetReal("problem", "t_cool_min");
  glob_T_cutoff_over_T_cold =
      pin->GetReal("problem", "T_cutoff_over_T_cold");
  glob_T_ci_over_T_cold = pin->GetReal("problem", "T_ci_over_T_cold");
  glob_T_ih_over_T_cold = pin->GetReal("problem", "T_ih_over_T_cold");
  glob_beta = pin->GetReal("problem", "beta");
  glob_velocity = pin->GetReal("problem", "velocity");
  glob_initial_vz = pin->GetOrAddReal("problem", "initial_vz", 0.0);
  glob_cooling_ramp_time =
      pin->GetOrAddReal("problem", "cooling_ramp_time", -1.0);
  glob_cooling_ramp_min_factor =
      pin->GetOrAddReal("problem", "cooling_ramp_min_factor", 0.0);
  if (glob_cooling_ramp_min_factor < 0.0 ||
      glob_cooling_ramp_min_factor > 1.0) {
    FatalSimpleTRMLInput(
        "Require 0 <= <problem>/cooling_ramp_min_factor <= 1.");
  }

  glob_shear_vel_thresh =
      pin->GetOrAddReal("problem", "hist_shear_vel_frac", 0.45)*glob_velocity;
  glob_vy_vel_thresh =
      pin->GetOrAddReal("problem", "hist_vy_vel_frac", 0.08)*glob_velocity;
  glob_custom_min_timestep =
      pin->GetOrAddReal("problem", "custom_min_timestep", -1.0);

  glob_cooling_rate = 0.0;
  glob_cooling_rk_u0 = 0.0;
  glob_cooling_rk_u1 = 0.0;

  Real rho_cold = glob_rho_cold;
  Real rho_hot = glob_rho_hot;
  Real pres = glob_pres;
  Real velocity = glob_velocity;
  Real initial_vz = glob_initial_vz;
  Real sharpness = pin->GetReal("problem", "phase_sharpness");
  uint max_init_p = pin->GetInteger("problem", "max_init_perturb_log2freq");
  uint min_init_p = pin->GetInteger("problem", "min_init_perturb_log2freq");
  Real init_p_sharp = pin->GetReal("problem", "init_perturb_sharpness");
  Real init_p_vel_frac = pin->GetReal("problem", "init_perturb_vel_frac");

  user_stage_srcs_func = CoolingSrc;
  user_time_step_func = CoolingTimestep;
  user_bcs_func = TRMLOuterX3Boundary;
  user_hist_func = HistoryOutput;

  if (restart) return;

  if (pmbp->phydro != nullptr) {
    auto &u0 = pmbp->phydro->u0;
    Real gm1 = eos.gamma - 1.0;
    int nhydro = pmbp->phydro->nhydro;
    int nscalars = pmbp->phydro->nscalars;

    par_for(
        "TRML_problem", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js,
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
          Real dens = shiftedtanh(coordz*sharpness, rho_cold, rho_hot);
          Real cold_fraction = shiftedtanh(coordz*sharpness, 1.0, 0.0);
          Real vx = shiftedtanh(coordz*sharpness, -0.5*velocity,
                                0.5*velocity);

          u0(m, IDN, k, j, i) = dens;
          u0(m, IM1, k, j, i) = dens*vx;
          u0(m, IM2, k, j, i) = 0.0;
          // initial_vz is the stored velocity itself: no sign conversion.
          u0(m, IM3, k, j, i) = dens*initial_vz;

          for (int i0 = min_init_p; i0 <= max_init_p; ++i0) {
            for (int j0 = min_init_p; j0 <= max_init_p; ++j0) {
              u0(m, IM3, k, j, i) +=
                  velocity*init_p_vel_frac*
                  sin(2.0*3.141592653589*coordy*pow(2, j0))*
                  sin(2.0*3.141592653589*coordx*pow(2, i0))*dens*
                  exp(-coordz*coordz*init_p_sharp/2.0);
            }
          }

          if (eos.is_ideal) {
            u0(m, IEN, k, j, i) =
                pres/gm1 +
                0.5*(SQR(u0(m, IM1, k, j, i)) +
                     SQR(u0(m, IM2, k, j, i)) +
                     SQR(u0(m, IM3, k, j, i)))/dens;
          }
          for (int n = nhydro; n < nhydro + nscalars; ++n) {
            u0(m, n, k, j, i) = dens*cold_fraction;
          }
        });
  }
}

void CoolingSrc(Mesh *pm, Real bdt, Driver *pdrive, int stage) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  Real rho_cold = glob_rho_cold;
  Real t_cool_min = glob_t_cool_min;
  Real T_cutoff_over_T_cold = glob_T_cutoff_over_T_cold;
  Real beta = glob_beta;
  Real T_cold = glob_pres/rho_cold;
  Real cooling_bdt = bdt*CoolingRampFactorAtTime(pm->time);

  Real cooling_energy_this_stage;
  Kokkos::parallel_reduce(
      "TRML_cooling_src", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &stage_cooling) {
        int m = idx/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/nx1;
        int i = (idx - m*nkji - k*nji - j*nx1) + is;
        k += ks;
        j += js;

        Real dens = w0(m, IDN, k, j, i);
        Real eint = w0(m, IEN, k, j, i);
        Real temp = eint/dens*gm1;
        Real temp_at_t_cool_min =
            (beta > 0.0) ? T_cold : T_cutoff_over_T_cold*T_cold;

        Real deltaE;
        if (beta != 0.0) {
          Real partofcooling =
              1.0 - pow(temp/temp_at_t_cool_min, -beta)*
                        cooling_bdt/t_cool_min*beta;
          if (partofcooling > 0.0) {
            deltaE = eint*pow(partofcooling, 1.0/beta) - eint;
          } else {
            deltaE = dens*T_cold/gm1 - eint;
          }
        } else {
          deltaE = eint*exp(-cooling_bdt/t_cool_min) - eint;
        }

        const Real eint_floor = dens*T_cold/gm1;
        if (temp <= T_cold || temp > T_cutoff_over_T_cold*T_cold) {
          deltaE = 0.0;
        } else {
          deltaE = fmax(deltaE, eint_floor - eint);
        }

        Real vol = size.d_view(m).dx1*size.d_view(m).dx2*
                   size.d_view(m).dx3;
        stage_cooling += -deltaE*vol;
        u0(m, IEN, k, j, i) += deltaE;
      },
      Kokkos::Sum<Real>(cooling_energy_this_stage));

  if (stage == 1) {
    glob_cooling_rk_u0 = 0.0;
    glob_cooling_rk_u1 = 0.0;
  } else if (pdrive->integrator == "rk4") {
    glob_cooling_rk_u1 += pdrive->delta[stage - 1]*glob_cooling_rk_u0;
  }

  glob_cooling_rk_u0 =
      pdrive->gam0[stage - 1]*glob_cooling_rk_u0 +
      pdrive->gam1[stage - 1]*glob_cooling_rk_u1 +
      cooling_energy_this_stage;

  if (stage == pdrive->nexp_stages) {
    glob_cooling_rate = glob_cooling_rk_u0/pm->dt;
  }
}

void CoolingTimestep(Mesh *pm) {
  if (glob_custom_min_timestep > 0.0) {
    pm->pgen->dtnew = glob_custom_min_timestep;
  } else {
    pm->pgen->dtnew = glob_t_cool_min;
  }
}

void TRMLOuterX3Boundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ke = indcs.ke;

  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  int nhydro = pmbp->phydro->nhydro;
  int nscalars = pmbp->phydro->nscalars;
  Real rho_hot = glob_rho_hot;
  Real velocity = glob_velocity;
  Real pres = glob_pres;

  par_for(
      "TRML_outer_x3_boundary", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
      js, je, is, ie, KOKKOS_LAMBDA(int m, int j, int i) {
        if (mb_bcs.d_view(m, BoundaryFace::outer_x3) != BoundaryFlag::user) {
          return;
        }
        const Real active_rho = u0(m, IDN, ke, j, i);
        const Real inv_active_rho =
            (active_rho > 0.0) ? 1.0/active_rho : 0.0;
        const Real vy = u0(m, IM2, ke, j, i)*inv_active_rho;
        const Real vz = u0(m, IM3, ke, j, i)*inv_active_rho;

        for (int k = 0; k < ng; ++k) {
          int ghost_k = ke + k + 1;
          u0(m, IDN, ghost_k, j, i) = rho_hot;
          u0(m, IM1, ghost_k, j, i) = rho_hot*(0.5*velocity);
          u0(m, IM2, ghost_k, j, i) = rho_hot*vy;
          u0(m, IM3, ghost_k, j, i) = rho_hot*vz;
          u0(m, IEN, ghost_k, j, i) =
              pres/gm1 + 0.5*rho_hot*
              (SQR(0.5*velocity) + SQR(vy) + SQR(vz));
          for (int n = nhydro; n < nhydro + nscalars; ++n) {
            u0(m, n, ghost_k, j, i) = 0.0;
          }
        }
      });
}

#define TRMLHISTVARS 25

namespace Kokkos {
template <>
struct reduction_identity<array_sum::array_type<Real, TRMLHISTVARS>> {
  KOKKOS_FORCEINLINE_FUNCTION
  static array_sum::array_type<Real, TRMLHISTVARS> sum() {
    return array_sum::array_type<Real, TRMLHISTVARS>();
  }
};
} // namespace Kokkos

void HistoryOutput(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  Real lxmax3 = pm->mesh_size.x3max;
  Real lxmin3 = pm->mesh_size.x3min;

  auto &w0 = pmbp->phydro->w0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  Real T_cold = glob_pres/glob_rho_cold;
  Real T_ci_over_T_cold = glob_T_ci_over_T_cold;
  Real T_ih_over_T_cold = glob_T_ih_over_T_cold;
  Real shear_vel_thresh = glob_shear_vel_thresh;
  Real vy_vel_thresh = glob_vy_vel_thresh;

  pdata->nhist = TRMLHISTVARS;
  pdata->label[0] = "cooling_rate ";
  pdata->label[1] = "M_flux_top ";
  pdata->label[2] = "M_flux_bot ";
  pdata->label[3] = "E_flux_top ";
  pdata->label[4] = "E_flux_bot ";
  pdata->label[5] = "vysq_c ";
  pdata->label[6] = "vysq_i ";
  pdata->label[7] = "vysq_h ";
  pdata->label[8] = "kex_c ";
  pdata->label[9] = "kex_i ";
  pdata->label[10] = "kex_h ";
  pdata->label[11] = "momx_c ";
  pdata->label[12] = "momx_i ";
  pdata->label[13] = "momx_h ";
  pdata->label[14] = "V_c ";
  pdata->label[15] = "V_i ";
  pdata->label[16] = "V_h ";
  pdata->label[17] = "zavg_i ";
  pdata->label[18] = "zmin_shear ";
  pdata->label[19] = "zmax_shear ";
  pdata->label[20] = "zmin_intermediate ";
  pdata->label[21] = "zmax_intermediate ";
  pdata->label[22] = "zmin_vy ";
  pdata->label[23] = "zmax_vy ";
  pdata->label[24] = "cool_ramp ";

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
      KOKKOS_LAMBDA(
          const int &idx,
          array_sum::array_type<Real, TRMLHISTVARS> &mb_sum,
          Real &min_z_intermediate_, Real &max_z_intermediate_,
          Real &min_z_shear_, Real &max_z_shear_, Real &min_z_vely_,
          Real &max_z_vely_, int &cnt_z_intermediate_, int &cnt_z_shear_,
          int &cnt_z_vely_) {
        int m = idx/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/nx1;
        int i = (idx - m*nkji - k*nji - j*nx1) + is;
        k += ks;
        j += js;

        Real vol = size.d_view(m).dx1*size.d_view(m).dx2*
                   size.d_view(m).dx3;
        Real x3v = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
        Real dz = size.d_view(m).dx3;
        Real dA = size.d_view(m).dx1*size.d_view(m).dx2;

        Real dens = w0(m, IDN, k, j, i);
        Real velx = w0(m, IVX, k, j, i);
        Real vely = w0(m, IVY, k, j, i);
        Real velz = w0(m, IVZ, k, j, i);
        Real eint = w0(m, IEN, k, j, i);
        Real temp = eint/dens*gm1;
        Real ekin = 0.5*dens*(SQR(velx) + SQR(vely) + SQR(velz));
        Real etot = eint + ekin;
        Real pressure = gm1*eint;

        if (fabs(x3v - lxmax3) < dz) {
          mb_sum.the_array[1] += dens*dA*velz;
          mb_sum.the_array[3] += (etot + pressure)*dA*velz;
        }
        if (fabs(x3v - lxmin3) < dz) {
          mb_sum.the_array[2] -= dens*dA*velz;
          mb_sum.the_array[4] -= (etot + pressure)*dA*velz;
        }

        int temp_bin;
        if (temp/T_cold < T_ci_over_T_cold) {
          temp_bin = 0;
        } else if (temp/T_cold < T_ih_over_T_cold) {
          temp_bin = 1;
        } else {
          temp_bin = 2;
        }

        mb_sum.the_array[5 + temp_bin] += SQR(vely)*vol;
        mb_sum.the_array[8 + temp_bin] += 0.5*dens*SQR(velx)*vol;
        mb_sum.the_array[11 + temp_bin] += dens*velx*vol;
        mb_sum.the_array[14 + temp_bin] += vol;
        if (temp_bin == 1) {
          mb_sum.the_array[17] += x3v*vol;
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
  Real local_min[3] = {min_z_intermediate, min_z_shear, min_z_vely};
  Real local_max[3] = {max_z_intermediate, max_z_shear, max_z_vely};
  Real global_min[3];
  Real global_max[3];
  int local_counts[3] = {cnt_z_intermediate, cnt_z_shear, cnt_z_vely};
  int global_counts[3] = {0, 0, 0};
  MPI_Allreduce(local_min, global_min, 3, MPI_ATHENA_REAL, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(local_max, global_max, 3, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(local_counts, global_counts, 3, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  min_z_intermediate = global_min[0];
  min_z_shear = global_min[1];
  min_z_vely = global_min[2];
  max_z_intermediate = global_max[0];
  max_z_shear = global_max[1];
  max_z_vely = global_max[2];
  cnt_z_intermediate = global_counts[0];
  cnt_z_shear = global_counts[1];
  cnt_z_vely = global_counts[2];
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

  for (int n = 1; n <= 17; ++n) {
    pdata->hdata[n] = history_sum.the_array[n];
  }
  pdata->hdata[0] = glob_cooling_rate;

  if (global_variable::my_rank == 0) {
    pdata->hdata[18] = min_z_shear;
    pdata->hdata[19] = max_z_shear;
    pdata->hdata[20] = min_z_intermediate;
    pdata->hdata[21] = max_z_intermediate;
    pdata->hdata[22] = min_z_vely;
    pdata->hdata[23] = max_z_vely;
    pdata->hdata[24] = CoolingRampFactorAtTime(pm->time);
  } else {
    for (int n = 18; n <= 24; ++n) {
      pdata->hdata[n] = 0.0;
    }
  }
}
