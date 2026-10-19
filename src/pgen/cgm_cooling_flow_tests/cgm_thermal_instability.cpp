//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgm_thermal_instability.cpp
//! \brief Thermal-instability test of the CGM cooling (<hydro>/cgm_cooling = true).
//! A uniform medium at (nH0, T0) is made an exact equilibrium of the code's own net rate
//! by setting the heating amplitude, then seeded with an isobaric density perturbation:
//! a single mode (mode_x, mode_y) and/or white noise. The predicted isobaric growth rate
//! of that rate function is printed; history adds the mode amplitude and the cold gas.
//! With balance = false the input hrate is kept (a cooling medium). churn_speed > 0 makes
//! the mesh refine a band in x1 moving at that speed and derefine elsewhere (AMR churn).

#include <cmath>
#include <iostream>
#include <random>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "mesh/mesh_refinement.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/cgm_cooling.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"

namespace {
Real rho0, t_cold, kx, ky;  // background density, cold threshold [code], wavenumbers
Real churn_v, churn_w;      // speed and half width of the refined band
void TIHistory(HistoryData *pdata, Mesh *pm);
void ChurnRefinement(MeshBlockPack *pmbp);
}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  user_hist_func = TIHistory;
  auto *psrc = (pmbp->phydro != nullptr) ? pmbp->phydro->psrc : nullptr;
  if (psrc == nullptr || !psrc->cgm_cooling) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "cgm_thermal_instability needs <hydro> with cgm_cooling = true"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &units = pmbp->punit;
  const Real temp_unit = units->temperature_cgs();
  const Real nH_unit = units->density_cgs()/units->atomic_mass_unit_cgs;
  const Real X = 0.75;  // as in cgm_cooling::CalculateRates
  const Real gamma = pmbp->phydro->peos->eos_data.gamma, gm1 = gamma - 1.0;
  rho0 = pin->GetReal("problem", "nH0")/(X*nH_unit);
  const Real temp0 = pin->GetReal("problem", "T0")/temp_unit;  // code T = P/rho
  const Real zmet = pin->GetOrAddReal("problem", "metallicity", 0.02/3.0);
  t_cold = pin->GetOrAddReal("problem", "T_cold", 1.0e3)/temp_unit;
  auto &size = pmy_mesh_->mesh_size;
  kx = 2.0*M_PI*pin->GetOrAddInteger("problem", "mode_x", 0)/(size.x1max - size.x1min);
  ky = 2.0*M_PI*pin->GetOrAddInteger("problem", "mode_y", 0)/(size.x2max - size.x2min);

  // Net cooling rate per volume [code] at (rho, T) at the box centre, evaluated with the
  // same host tables and arguments CGMCooling uses. It is linear in hrate.
  const Real xc = 0.5*(size.x1min + size.x1max), yc = 0.5*(size.x2min + size.x2max);
  const Real zc = 0.5*(size.x3min + size.x3max);
  const Real dx1 = (size.x1max - size.x1min)/pmy_mesh_->mesh_indcs.nx1;
  auto net = [&](Real rho, Real temp) {
    return cgm_cooling::CalculateRates(rho, rho*temp/gm1, zmet, xc, yc, zc, dx1, gm1,
        temp_unit, nH_unit, units->pressure_cgs()/units->time_cgs()/nH_unit/nH_unit,
        units->pressure_cgs()/units->time_cgs()/nH_unit, units->length_cgs(),
        psrc->hrate, psrc->hscale_norm, psrc->hscale_height, psrc->hscale_radius,
        psrc->T_max, std::pow(10.0, Tbins_ARR[0]),
        std::pow(10.0, Tbins_ARR[Tbins_DIM_0-1]),
        std::pow(10.0, nHbins_ARR[0]), std::pow(10.0, nHbins_ARR[nHbins_DIM_0-1]),
        psrc->Tbins.h_view, psrc->nHbins.h_view, psrc->Metal_Cooling.h_view,
        psrc->H_He_Cooling.h_view, psrc->Metal_Cooling_CIE.h_view,
        psrc->H_He_Cooling_CIE.h_view).net_cooling_rate;
  };
  const Real hrate_in = psrc->hrate;
  const bool balance = pin->GetOrAddBoolean("problem", "balance", true);
  psrc->hrate = 0.0;
  const Real r_cool = net(rho0, temp0);
  psrc->hrate = 1.0;
  psrc->hrate = balance ? r_cool/(r_cool - net(rho0, temp0)) : hrate_in;
  // isobaric growth rate: -(gamma-1)/(gamma rho T) [T dL/dT|rho - rho dL/drho|T]
  const Real d = 1.0e-4;
  const Real dldt = (net(rho0, temp0*(1+d)) - net(rho0, temp0*(1-d)))/(2*d);
  const Real dldr = (net(rho0*(1+d), temp0) - net(rho0*(1-d), temp0))/(2*d);
  const Real t_cool = rho0*temp0/gm1/r_cool;
  const Real sigma = -gm1/(gamma*rho0*temp0)*(dldt - dldr);
  if (global_variable::my_rank == 0) {
    std::cout << "TI: hrate = " << psrc->hrate << "  t_cool = " << t_cool
              << "  sigma_isobaric = " << sigma << "  sigma*t_cool = " << sigma*t_cool
              << std::endl;
  }
  churn_v = pin->GetOrAddReal("problem", "churn_speed", 0.0);
  churn_w = pin->GetOrAddReal("problem", "churn_width", 0.0);
  if (churn_v > 0.0) user_ref_func = ChurnRefinement;
  if (restart) return;

  const Real amp = pin->GetOrAddReal("problem", "amp", 1.0e-3);
  const Real noise = pin->GetOrAddReal("problem", "noise", 0.0);
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nscal = pmbp->phydro->nscalars, nhyd = pmbp->phydro->nhydro;
  auto u0 = Kokkos::create_mirror_view(pmbp->phydro->u0);
  auto &mbsize = pmbp->pmb->mb_size.h_view;
  std::mt19937 rng(pin->GetOrAddInteger("problem", "seed", 1) + global_variable::my_rank);
  std::uniform_real_distribution<Real> uni(-1.0, 1.0);
  const Real e0 = rho0*temp0/gm1;  // isobaric: uniform pressure
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    for (int k=ks; k<ks+nx3; ++k)
    for (int j=js; j<js+nx2; ++j)
    for (int i=is; i<is+nx1; ++i) {
      Real x = CellCenterX(i-is, nx1, mbsize(m).x1min, mbsize(m).x1max);
      Real y = CellCenterX(j-js, nx2, mbsize(m).x2min, mbsize(m).x2max);
      Real rho = rho0*(1.0 + amp*std::cos(kx*x + ky*y) + noise*uni(rng));
      u0(m,IDN,k,j,i) = rho;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = e0;
      for (int n=0; n<nscal; ++n) u0(m,nhyd+n,k,j,i) = (n == 0) ? zmet*rho : 0.0;
    }
  }
  Kokkos::deep_copy(pmbp->phydro->u0, u0);
}

namespace {
//! history: mode amplitude (cos, sin parts of 2<(rho-rho0) e^{-ik.x}>/rho0) and cold gas
void TIHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 4;
  pdata->label[0] = "amp_cos";
  pdata->label[1] = "amp_sin";
  pdata->label[2] = "M_cold";
  pdata->label[3] = "V_cold";
  auto &w0 = pm->pmb_pack->phydro->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1, js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const Real gm1 = pm->pmb_pack->phydro->peos->eos_data.gamma - 1.0;
  const Real vol_tot = (pm->mesh_size.x1max - pm->mesh_size.x1min)*
                       (pm->mesh_size.x2max - pm->mesh_size.x2min)*
                       (pm->mesh_size.x3max - pm->mesh_size.x3min);
  const Real r0 = rho0, tc = t_cold, kx_ = kx, ky_ = ky;
  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("TIHist",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, pm->pmb_pack->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real rho = w0(m,IDN,k,j,i);
    bool cold = (gm1*w0(m,IEN,k,j,i)/rho < tc);
    array_sum::GlobalSum h;
    h.the_array[0] = 2.0*(rho - r0)*cos(kx_*x + ky_*y)*vol/(r0*vol_tot);
    h.the_array[1] = 2.0*(rho - r0)*sin(kx_*x + ky_*y)*vol/(r0*vol_tot);
    h.the_array[2] = cold ? rho*vol : 0.0;
    h.the_array[3] = cold ? vol : 0.0;
    for (int n=4; n<NHISTORY_VARIABLES; ++n) h.the_array[n] = 0.0;
    mb_sum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  for (int n=0; n<4; ++n) pdata->hdata[n] = sum.the_array[n];
}

//! refine MeshBlocks whose centre is within churn_w of a band moving along x1 at churn_v
void ChurnRefinement(MeshBlockPack *pmbp) {
  Mesh *pm = pmbp->pmesh;
  const int mbs = pm->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pm->pmr->refine_flag;
  auto &size = pmbp->pmb->mb_size;
  const Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min, w = churn_w;
  const Real xb = pm->mesh_size.x1min + std::fmod(churn_v*pm->time, lx);
  Kokkos::parallel_for("ChurnRef", Kokkos::RangePolicy<>(DevExeSpace(), 0,
      pmbp->nmb_thispack), KOKKOS_LAMBDA(const int m) {
    Real d = fabs(0.5*(size.d_view(m).x1min + size.d_view(m).x1max) - xb);
    refine_flag.d_view(m+mbs) = (fmin(d, lx - d) < w) ? 1 : -1;
  });
}
}  // namespace
