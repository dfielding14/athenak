//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cgl_lf_paper.cpp
//! \brief CGL-Landau-fluid paper-validation problem generator.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/conduction.hpp"
#include "outputs/outputs.hpp"
#include "pgen.hpp"

namespace {

enum class CglLfPaperMode { turbulence=0, np_mode=1, fast_wave=2,
                            oblique_iaw=3, linear_wave_scan=4 };

void CglLfPaperHistory(HistoryData *pdata, Mesh *pm);

[[noreturn]] void FatalInput(const char *msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

CglLfPaperMode ParseMode(const std::string &mode) {
  if (mode == "turbulence") return CglLfPaperMode::turbulence;
  if (mode == "np_mode") return CglLfPaperMode::np_mode;
  if (mode == "fast_wave") return CglLfPaperMode::fast_wave;
  if (mode == "oblique_iaw") return CglLfPaperMode::oblique_iaw;
  if (mode == "linear_wave_scan") return CglLfPaperMode::linear_wave_scan;
  FatalInput("<problem>/mode must be turbulence, np_mode, fast_wave, oblique_iaw, "
             "or linear_wave_scan");
}

void ValidateForcingMode(const std::string &forcing) {
  if (forcing == "alfvenic" || forcing == "random" || forcing == "sonic_corr") return;
  FatalInput("<problem>/forcing_mode must be alfvenic, random, or sonic_corr");
}

KOKKOS_INLINE_FUNCTION
Real PositiveFloor(const Real x, const Real floor) {
  return (x > floor) ? x : floor;
}

KOKKOS_INLINE_FUNCTION
Real CglLfPaperPhase(const Real x1, const Real x3, const Real x1min,
                     const Real x3min, const Real kx, const Real kz) {
  return kx*(x1 - x1min) + kz*(x3 - x3min);
}

KOKKOS_INLINE_FUNCTION
void CglLfPaperBField(const int mode, const Real x1, const Real x3,
                      const Real x1min, const Real x3min, const Real dx1,
                      const Real dx3, const Real b0, const Real amp,
                      const Real kx, const Real kz, const Real kmag,
                      Real &bx, Real &by, Real &bz) {
  bx = 0.0;
  by = 0.0;
  bz = b0;

  if (mode == static_cast<int>(CglLfPaperMode::fast_wave)) {
    const Real ph = CglLfPaperPhase(x1, x3, x1min, x3min, kx, 0.0);
    bz = b0*(1.0 + amp*sin(ph));
    return;
  }

  if (mode == static_cast<int>(CglLfPaperMode::linear_wave_scan)) {
    const Real ph = CglLfPaperPhase(x1, x3, x1min, x3min, kx, kz);
    by = amp*b0*sin(ph);
    return;
  }

  if (mode == static_cast<int>(CglLfPaperMode::np_mode) ||
      mode == static_cast<int>(CglLfPaperMode::oblique_iaw)) {
    const Real ph = CglLfPaperPhase(x1, x3, x1min, x3min, kx, kz);
    const Real field_amp = (mode == static_cast<int>(CglLfPaperMode::np_mode)) ?
                           amp : 0.1*amp;
    const Real ay_amp = (kmag > 0.0) ? field_amp*b0/kmag : 0.0;
    const Real fd1 = (dx1 > 0.0) ? 2.0*sin(0.5*kx*dx1)/dx1 : kx;
    const Real fd3 = (dx3 > 0.0) ? 2.0*sin(0.5*kz*dx3)/dx3 : kz;
    bx = -ay_amp*fd3*cos(ph);
    bz =  b0 + ay_amp*fd1*cos(ph);
  }
}

KOKKOS_INLINE_FUNCTION
Real CglLfPaperLimiterNu(const Real ppar, const Real pperp, const Real bx,
                         const Real by, const Real bz, const Real lim_coll,
                         const bool mlim, const bool flim, const bool backup_lim) {
  const Real paniso = pperp - ppar;
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  const Real limiter_nu = fmax(lim_coll, static_cast<Real>(0.0));
  const Real backup_nu = static_cast<Real>(1.0e10);
  Real nu_eff = 0.0;

  if (flim && backup_lim) {
    if ((paniso <= static_cast<Real>(-0.7)*bsqr) && (paniso > -bsqr)) {
      nu_eff = fmax(nu_eff, limiter_nu);
    } else if (paniso <= -bsqr) {
      nu_eff = fmax(nu_eff, backup_nu);
    }
  } else if (flim && paniso <= static_cast<Real>(-0.7)*bsqr) {
    nu_eff = fmax(nu_eff, limiter_nu);
  }

  if (mlim && backup_lim) {
    if ((paniso >= static_cast<Real>(0.5)*bsqr) && (paniso < bsqr)) {
      nu_eff = fmax(nu_eff, limiter_nu);
    } else if (paniso >= bsqr) {
      nu_eff = fmax(nu_eff, backup_nu);
    }
  } else if (mlim && paniso >= static_cast<Real>(0.5)*bsqr) {
    nu_eff = fmax(nu_eff, limiter_nu);
  }

  return nu_eff;
}

KOKKOS_INLINE_FUNCTION
Real CglLfPaperChiPerp(const Real cpar, const Real lf_k_parallel, const Real nu_eff) {
  const Real denom = static_cast<Real>(2.5066282746310002)*cpar*lf_k_parallel + nu_eff;
  return (denom > 0.0) ? static_cast<Real>(2.0)*SQR(cpar)/denom : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real CglLfPaperChiParallel(const Real cpar, const Real lf_k_parallel,
                           const Real nu_eff) {
  const Real denom = static_cast<Real>(5.013256549262000)*cpar*lf_k_parallel
                   + static_cast<Real>(1.4247779607693793)*nu_eff;
  return (denom > 0.0) ? static_cast<Real>(8.0)*SQR(cpar)/denom : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real CglLfPaperCappedFlux(const Real q_unlimited, const Real q_max) {
  return (q_max > 0.0) ? q_unlimited*q_max/(q_max + fabs(q_unlimited)) : 0.0;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief CGL-LF validation initial data for turbulence and Majeski/Squire wave analogues.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_hist_func = CglLfPaperHistory;
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->pmhd == nullptr || pmbp->phydro != nullptr) {
    FatalInput("cgl_lf_paper requires an <mhd> block and no <hydro> block");
  }
  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  if (!eos.is_cgl) {
    FatalInput("cgl_lf_paper requires <mhd>/eos = cgl");
  }

  const std::string mode_name = pin->GetOrAddString("problem", "mode", "turbulence");
  const int mode_id = static_cast<int>(ParseMode(mode_name));
  const std::string forcing_mode =
      pin->GetOrAddString("problem", "forcing_mode", "alfvenic");
  ValidateForcingMode(forcing_mode);

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real b0 = pin->GetOrAddReal("problem", "B0", 1.0);
  const Real beta0 = pin->GetOrAddReal("problem", "beta0", 10.0);
  const Real p_beta = 0.5*beta0*SQR(b0);
  const Real ppar0 = pin->GetOrAddReal("problem", "p_parallel0", p_beta);
  const Real pperp0 = pin->GetOrAddReal("problem", "p_perp0", p_beta);
  const Real amp = pin->DoesParameterExist("problem", "alpha") ?
                   pin->GetReal("problem", "alpha") :
                   pin->GetOrAddReal("problem", "amp", 1.0e-3);
  const Real velocity_amp = pin->GetOrAddReal("problem", "velocity_amp", amp);

  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  const Real x3min = pmy_mesh_->mesh_size.x3min;
  const Real x3max = pmy_mesh_->mesh_size.x3max;
  const Real lx1 = x1max - x1min;
  const Real lx3 = x3max - x3min;
  const Real kpar_mode = pin->GetOrAddReal("problem", "k_parallel_mode", 1.0);
  const Real kperp_over_kpar = pin->GetOrAddReal("problem", "k_perp_over_k_parallel", 4.0);
  const Real default_kperp_mode = kperp_over_kpar*kpar_mode*lx1/lx3;
  const Real kperp_mode = pin->GetOrAddReal("problem", "k_perp_mode", default_kperp_mode);
  const Real kx = 2.0*M_PI*kperp_mode/lx1;
  const Real kz = 2.0*M_PI*kpar_mode/lx3;
  const Real kmag = std::sqrt(SQR(kx) + SQR(kz));

  if (rho0 <= 0.0 || b0 <= 0.0 || beta0 <= 0.0 || ppar0 <= 0.0 || pperp0 <= 0.0) {
    FatalInput("rho0, B0, beta0, p_parallel0, and p_perp0 must be positive");
  }

  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &b = pmbp->pmhd->b0;
  auto &bcc = pmbp->pmhd->bcc0;
  const Real pfloor = eos.pfloor;
  const Real dfloor = eos.dfloor;
  const Real vA0 = b0/std::sqrt(rho0);

  par_for("cgl_lf_paper_prim", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
  ks, ke, js, je, is, ie, KOKKOS_LAMBDA(const int m, const int k, const int j,
                                        const int i) {
    const Real x1v = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real dx1 = size.d_view(m).dx1;
    const Real dx3 = size.d_view(m).dx3;

    Real bx, by, bz;
    CglLfPaperBField(mode_id, x1v, x3v, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                     kmag, bx, by, bz);

    Real rho = rho0;
    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    Real ppar = ppar0;
    Real pperp = pperp0;
    const Real ph = CglLfPaperPhase(x1v, x3v, x1min, x3min, kx, kz);

    if (mode_id == static_cast<int>(CglLfPaperMode::np_mode)) {
      const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
      const Real pbal = PositiveFloor(pperp0 + 0.5*(SQR(b0) - bsqr), pfloor);
      ppar = pbal;
      pperp = pbal;
    } else if (mode_id == static_cast<int>(CglLfPaperMode::fast_wave)) {
      const Real comp = amp*sin(ph);
      rho = PositiveFloor(rho0*(1.0 + comp), dfloor);
      vx = velocity_amp*vA0*sin(ph);
      ppar = PositiveFloor(ppar0*(1.0 + comp), pfloor);
      pperp = PositiveFloor(pperp0*(1.0 + comp), pfloor);
    } else if (mode_id == static_cast<int>(CglLfPaperMode::oblique_iaw)) {
      const Real comp = amp*sin(ph);
      rho = PositiveFloor(rho0*(1.0 + comp), dfloor);
      const Real inv_k = (kmag > 0.0) ? 1.0/kmag : 0.0;
      vx = velocity_amp*sqrt(ppar0/rho0)*kx*inv_k*sin(ph);
      vz = velocity_amp*sqrt(ppar0/rho0)*kz*inv_k*sin(ph);
      ppar = PositiveFloor(ppar0*(1.0 + comp), pfloor);
      pperp = PositiveFloor(pperp0*(1.0 + 0.25*comp), pfloor);
    } else if (mode_id == static_cast<int>(CglLfPaperMode::linear_wave_scan)) {
      vy = velocity_amp*vA0*sin(ph);
    }

    w0(m,IDN,k,j,i) = rho;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = vy;
    w0(m,IVZ,k,j,i) = vz;
    w0(m,IPR,k,j,i) = ppar;
    w0(m,IPP,k,j,i) = pperp;
  });

  par_for("cgl_lf_paper_b", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
  ks, ke, js, je, is, ie, KOKKOS_LAMBDA(const int m, const int k, const int j,
                                        const int i) {
    const Real x1v = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Real dx1 = size.d_view(m).dx1;
    const Real dx3 = size.d_view(m).dx3;
    Real bx, by, bz;

    CglLfPaperBField(mode_id, x1f, x3v, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                     kmag, bx, by, bz);
    b.x1f(m,k,j,i) = bx;
    if (i == ie) {
      const Real x1r = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
      CglLfPaperBField(mode_id, x1r, x3v, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                       kmag, bx, by, bz);
      b.x1f(m,k,j,i+1) = bx;
    }

    CglLfPaperBField(mode_id, x1v, x3v, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                     kmag, bx, by, bz);
    b.x2f(m,k,j,i) = by;
    if (j == je) b.x2f(m,k,j+1,i) = by;

    CglLfPaperBField(mode_id, x1v, x3f, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                     kmag, bx, by, bz);
    b.x3f(m,k,j,i) = bz;
    if (k == ke) {
      const Real x3r = LeftEdgeX(k+1-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
      CglLfPaperBField(mode_id, x1v, x3r, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                       kmag, bx, by, bz);
      b.x3f(m,k+1,j,i) = bz;
    }

    CglLfPaperBField(mode_id, x1v, x3v, x1min, x3min, dx1, dx3, b0, amp, kx, kz,
                     kmag, bx, by, bz);
    bcc(m,IBX,k,j,i) = bx;
    bcc(m,IBY,k,j,i) = by;
    bcc(m,IBZ,k,j,i) = bz;
  });

  auto &u0 = pmbp->pmhd->u0;
  pmbp->pmhd->peos->PrimToCons(w0, bcc, u0, is, ie, js, je, ks, ke);
}

namespace {

void CglLfPaperHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 20;
  pdata->label[0] = "vol";
  pdata->label[1] = "mass";
  pdata->label[2] = "ekin";
  pdata->label[3] = "emag";
  pdata->label[4] = "ecgl";
  pdata->label[5] = "b2";
  pdata->label[6] = "b4";
  pdata->label[7] = "beta";
  pdata->label[8] = "dp";
  pdata->label[9] = "absdp";
  pdata->label[10] = "mir";
  pdata->label[11] = "fire";
  pdata->label[12] = "bmir";
  pdata->label[13] = "bfire";
  pdata->label[14] = "nueff";
  pdata->label[15] = "qpar1";
  pdata->label[16] = "qperp1";
  pdata->label[17] = "qpar10";
  pdata->label[18] = "qperp10";
  pdata->label[19] = "hfpow";

  auto &w = pm->pmb_pack->pmhd->w0;
  auto &bcc = pm->pmb_pack->pmhd->bcc0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  EOS_Data eos = pm->pmb_pack->pmhd->peos->eos_data;
  Conduction *pcond = pm->pmb_pack->pmhd->pcond;

  const bool has_lf = (pcond != nullptr) && pcond->IsCGLLandauFluidHeatFlux();
  const Real lf_k = has_lf ? pcond->lf_k_parallel : 0.0;
  const bool lf_local = has_lf ? pcond->lf_coeff_local : true;
  const Real lf_cpar0 = has_lf ? pcond->lf_c_parallel0 : 0.0;
  const Real dfloor = eos.dfloor;
  const Real pfloor = eos.pfloor;
  const Real tfloor = eos.tfloor;
  const Real bfloor = eos.bfloor;
  const Real nu_coll = eos.nu_coll;
  const Real lim_coll = eos.lim_coll;
  const bool mlim = eos.mlim;
  const bool flim = eos.flim;
  const bool backup = eos.backup_lim;

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;

  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("cgl_lf_paper_history",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/nkji;
    const int k0 = (idx - m*nkji)/nji;
    const int j0 = (idx - m*nkji - k0*nji)/nx1;
    const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
    const int k = k0 + ks;
    const int j = j0 + js;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    const Real rho = PositiveFloor(w(m,IDN,k,j,i), dfloor);
    const Real vx = w(m,IVX,k,j,i);
    const Real vy = w(m,IVY,k,j,i);
    const Real vz = w(m,IVZ,k,j,i);
    const Real ppar = PositiveFloor(w(m,IPR,k,j,i), pfloor);
    const Real pperp = PositiveFloor(w(m,IPP,k,j,i), pfloor);
    const Real bx = bcc(m,IBX,k,j,i);
    const Real by = bcc(m,IBY,k,j,i);
    const Real bz = bcc(m,IBZ,k,j,i);
    const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
    const Real bmag = sqrt(bsqr);
    const Real p_iso = (ppar + 2.0*pperp)/3.0;
    const Real paniso = pperp - ppar;
    const Real beta = (bsqr > SQR(bfloor)) ? 2.0*p_iso/bsqr : 0.0;

    const bool mir = mlim && ((backup && paniso >= 0.5*bsqr && paniso < bsqr) ||
                              (!backup && paniso >= 0.5*bsqr));
    const bool fire = flim && ((backup && paniso <= -0.7*bsqr && paniso > -bsqr) ||
                               (!backup && paniso <= -0.7*bsqr));
    const bool bmir = mlim && backup && paniso >= bsqr;
    const bool bfire = flim && backup && paniso <= -bsqr;
    const Real nu_eff = fmax(nu_coll, static_cast<Real>(0.0)) +
        CglLfPaperLimiterNu(ppar, pperp, bx, by, bz, lim_coll, mlim, flim, backup);

    array_sum::GlobalSum hvars;
    hvars.the_array[0] = vol;
    hvars.the_array[1] = rho*vol;
    hvars.the_array[2] = 0.5*rho*(SQR(vx) + SQR(vy) + SQR(vz))*vol;
    hvars.the_array[3] = 0.5*bsqr*vol;
    hvars.the_array[4] = (0.5*ppar + pperp)*vol;
    hvars.the_array[5] = bsqr*vol;
    hvars.the_array[6] = SQR(bsqr)*vol;
    hvars.the_array[7] = beta*vol;
    hvars.the_array[8] = paniso*vol;
    hvars.the_array[9] = fabs(paniso)*vol;
    hvars.the_array[10] = mir ? vol : 0.0;
    hvars.the_array[11] = fire ? vol : 0.0;
    hvars.the_array[12] = bmir ? vol : 0.0;
    hvars.the_array[13] = bfire ? vol : 0.0;
    hvars.the_array[14] = nu_eff*vol;

    Real qpar_ratio = 0.0;
    Real qperp_ratio = 0.0;
    Real hfpow = 0.0;
    if (has_lf && lf_k > 0.0 && bmag > bfloor) {
      const Real tpar_ip = w(m,IPR,k,j,i+1)/PositiveFloor(w(m,IDN,k,j,i+1), dfloor);
      const Real tpar_im = w(m,IPR,k,j,i-1)/PositiveFloor(w(m,IDN,k,j,i-1), dfloor);
      const Real tperp_ip = w(m,IPP,k,j,i+1)/PositiveFloor(w(m,IDN,k,j,i+1), dfloor);
      const Real tperp_im = w(m,IPP,k,j,i-1)/PositiveFloor(w(m,IDN,k,j,i-1), dfloor);
      Real grad_tpar_x = 0.5*(tpar_ip - tpar_im)/size.d_view(m).dx1;
      Real grad_tperp_x = 0.5*(tperp_ip - tperp_im)/size.d_view(m).dx1;
      Real grad_b_x = 0.5*(sqrt(SQR(bcc(m,IBX,k,j,i+1)) +
                                      SQR(bcc(m,IBY,k,j,i+1)) +
                                      SQR(bcc(m,IBZ,k,j,i+1))) -
                             sqrt(SQR(bcc(m,IBX,k,j,i-1)) +
                                      SQR(bcc(m,IBY,k,j,i-1)) +
                                      SQR(bcc(m,IBZ,k,j,i-1))))/size.d_view(m).dx1;
      Real grad_tpar_y = 0.0;
      Real grad_tperp_y = 0.0;
      Real grad_b_y = 0.0;
      Real grad_tpar_z = 0.0;
      Real grad_tperp_z = 0.0;
      Real grad_b_z = 0.0;
      if (multi_d) {
        const Real tpar_jp = w(m,IPR,k,j+1,i)/PositiveFloor(w(m,IDN,k,j+1,i), dfloor);
        const Real tpar_jm = w(m,IPR,k,j-1,i)/PositiveFloor(w(m,IDN,k,j-1,i), dfloor);
        const Real tperp_jp = w(m,IPP,k,j+1,i)/PositiveFloor(w(m,IDN,k,j+1,i), dfloor);
        const Real tperp_jm = w(m,IPP,k,j-1,i)/PositiveFloor(w(m,IDN,k,j-1,i), dfloor);
        grad_tpar_y = 0.5*(tpar_jp - tpar_jm)/size.d_view(m).dx2;
        grad_tperp_y = 0.5*(tperp_jp - tperp_jm)/size.d_view(m).dx2;
        grad_b_y = 0.5*(sqrt(SQR(bcc(m,IBX,k,j+1,i)) +
                                   SQR(bcc(m,IBY,k,j+1,i)) +
                                   SQR(bcc(m,IBZ,k,j+1,i))) -
                          sqrt(SQR(bcc(m,IBX,k,j-1,i)) +
                                   SQR(bcc(m,IBY,k,j-1,i)) +
                                   SQR(bcc(m,IBZ,k,j-1,i))))/size.d_view(m).dx2;
      }
      if (three_d) {
        const Real tpar_kp = w(m,IPR,k+1,j,i)/PositiveFloor(w(m,IDN,k+1,j,i), dfloor);
        const Real tpar_km = w(m,IPR,k-1,j,i)/PositiveFloor(w(m,IDN,k-1,j,i), dfloor);
        const Real tperp_kp = w(m,IPP,k+1,j,i)/PositiveFloor(w(m,IDN,k+1,j,i), dfloor);
        const Real tperp_km = w(m,IPP,k-1,j,i)/PositiveFloor(w(m,IDN,k-1,j,i), dfloor);
        grad_tpar_z = 0.5*(tpar_kp - tpar_km)/size.d_view(m).dx3;
        grad_tperp_z = 0.5*(tperp_kp - tperp_km)/size.d_view(m).dx3;
        grad_b_z = 0.5*(sqrt(SQR(bcc(m,IBX,k+1,j,i)) +
                                   SQR(bcc(m,IBY,k+1,j,i)) +
                                   SQR(bcc(m,IBZ,k+1,j,i))) -
                          sqrt(SQR(bcc(m,IBX,k-1,j,i)) +
                                   SQR(bcc(m,IBY,k-1,j,i)) +
                                   SQR(bcc(m,IBZ,k-1,j,i))))/size.d_view(m).dx3;
      }

      const Real bhx = bx/bmag;
      const Real bhy = by/bmag;
      const Real bhz = bz/bmag;
      const Real gradpar_tpar = bhx*grad_tpar_x + bhy*grad_tpar_y + bhz*grad_tpar_z;
      const Real gradpar_tperp = bhx*grad_tperp_x + bhy*grad_tperp_y + bhz*grad_tperp_z;
      const Real gradpar_b = bhx*grad_b_x + bhy*grad_b_y + bhz*grad_b_z;
      const Real cpar = lf_local ? sqrt(fmax(ppar/rho, tfloor)) : lf_cpar0;
      const Real chi_par = CglLfPaperChiParallel(cpar, lf_k, nu_eff);
      const Real chi_perp = CglLfPaperChiPerp(cpar, lf_k, nu_eff);
      const Real qpar_l = -chi_par*rho*gradpar_tpar;
      const Real qperp_l = -chi_perp*(rho*gradpar_tperp -
          pperp*(1.0 - pperp/ppar)*gradpar_b/bmag);
      const Real qpar_max = static_cast<Real>(1.5957691216057308)*cpar*ppar;
      const Real qperp_max = static_cast<Real>(0.7978845608028654)*cpar*pperp;
      qpar_ratio = (qpar_max > 0.0) ? fabs(qpar_l)/qpar_max : 0.0;
      qperp_ratio = (qperp_max > 0.0) ? fabs(qperp_l)/qperp_max : 0.0;
      const Real qpar = CglLfPaperCappedFlux(qpar_l, qpar_max);
      const Real qperp = CglLfPaperCappedFlux(qperp_l, qperp_max);
      hfpow = -(qpar*gradpar_tpar + qperp*gradpar_tperp);
    }
    hvars.the_array[15] = (qpar_ratio > 1.0) ? vol : 0.0;
    hvars.the_array[16] = (qperp_ratio > 1.0) ? vol : 0.0;
    hvars.the_array[17] = (qpar_ratio > 10.0) ? vol : 0.0;
    hvars.the_array[18] = (qperp_ratio > 10.0) ? vol : 0.0;
    hvars.the_array[19] = hfpow*vol;
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
}

} // namespace
