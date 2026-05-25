//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_lf_paper.cpp
//! \brief Reduced/current-API initializer for MKS24 CGL-LF forced-turbulence runs.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_physics.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/turb_driver.hpp"

namespace {

void Fail(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void CGLLFPaperHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 22;
  pdata->label[0] = "volume";
  pdata->label[1] = "p_parallel";
  pdata->label[2] = "p_perp";
  pdata->label[3] = "force_prp2";
  pdata->label[4] = "force_prl2";
  pdata->label[5] = "vel_prp2";
  pdata->label[6] = "mass";
  pdata->label[7] = "kinetic";
  pdata->label[8] = "magnetic";
  pdata->label[9] = "therm_cgl";
  pdata->label[10] = "b2";
  pdata->label[11] = "b4";
  pdata->label[12] = "delta_p";
  pdata->label[13] = "abs_dp";
  pdata->label[14] = "beta";
  pdata->label[15] = "mirror_vol";
  pdata->label[16] = "fire_vol";
  pdata->label[17] = "hard_vol";
  pdata->label[18] = "nu_eff";
  pdata->label[19] = "force_pwr";
  pdata->label[20] = "vel_prl2";
  pdata->label[21] = "force_work";

  auto *pmbp = pm->pmb_pack;
  auto w0 = pmbp->pmhd->w0;
  auto bcc0 = pmbp->pmhd->bcc0;
  const auto eos = pmbp->pmhd->peos->eos_data;
  DvceArray5D<Real> force;
  const bool has_force = (pmbp->pturb != nullptr);
  if (has_force) force = pmbp->pturb->force;
  auto size = pmbp->pmb->mb_size;
  const auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("cgl_lf_paper_hist", Kokkos::RangePolicy<>(DevExeSpace(),
  0, nmkji), KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &result) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    const Real rho = w0(m,IDN,k,j,i);
    const Real vx = w0(m,IVX,k,j,i);
    const Real vy = w0(m,IVY,k,j,i);
    const Real vz = w0(m,IVZ,k,j,i);
    const Real ppar = w0(m,IPR,k,j,i);
    const Real pperp = w0(m,IPP,k,j,i);
    const Real bx = bcc0(m,IBX,k,j,i);
    const Real by = bcc0(m,IBY,k,j,i);
    const Real bz = bcc0(m,IBZ,k,j,i);
    const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
    const Real paniso = pperp - ppar;
    const Real piso = ONE_3RD*ppar + TWO_3RDS*pperp;
    const Real beta = 2.0*piso/fmax(bsqr, SQR(eos.bfloor));
    const Real nu_eff = eos.nu_coll +
        cgl::LimiterCollisionRate(ppar, pperp, bsqr, eos.lim_coll, eos.mlim,
                                  eos.flim, eos.firehose_threshold, eos.backup_lim);
    array_sum::GlobalSum cell;
    cell.the_array[0] = vol;
    cell.the_array[1] = vol*ppar;
    cell.the_array[2] = vol*pperp;
    cell.the_array[3] = has_force ?
        vol*(SQR(force(m,0,k,j,i)) + SQR(force(m,1,k,j,i))) : 0.0;
    cell.the_array[4] = has_force ? vol*SQR(force(m,2,k,j,i)) : 0.0;
    cell.the_array[5] = vol*(SQR(vx) + SQR(vy));
    cell.the_array[6] = vol*rho;
    cell.the_array[7] = 0.5*vol*rho*(SQR(vx) + SQR(vy) + SQR(vz));
    cell.the_array[8] = 0.5*vol*bsqr;
    cell.the_array[9] = vol*(pperp + 0.5*ppar);
    cell.the_array[10] = vol*bsqr;
    cell.the_array[11] = vol*SQR(bsqr);
    cell.the_array[12] = vol*paniso;
    cell.the_array[13] = vol*fabs(paniso);
    cell.the_array[14] = vol*beta;
    cell.the_array[15] = cgl::MirrorLimiterActive(paniso, bsqr) ? vol : 0.0;
    cell.the_array[16] =
        cgl::FirehoseLimiterActive(paniso, bsqr, eos.firehose_threshold) ? vol : 0.0;
    cell.the_array[17] =
        (cgl::MirrorHardBoundViolated(paniso, bsqr) ||
         cgl::FirehoseHardBoundViolated(paniso, bsqr)) ? vol : 0.0;
    cell.the_array[18] = vol*nu_eff;
    cell.the_array[19] = has_force ?
        vol*rho*(vx*force(m,0,k,j,i) + vy*force(m,1,k,j,i) +
                 vz*force(m,2,k,j,i)) : 0.0;
    cell.the_array[20] = vol*SQR(vz);
    for (int n = 21; n < NREDUCTION_VARIABLES; ++n) {
      cell.the_array[n] = 0.0;
    }
    result += cell;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  Kokkos::fence();

  for (int n = 0; n < pdata->nhist; ++n) {
    pdata->hdata[n] = sum.the_array[n];
  }
  pdata->hdata[21] =
      (has_force && pmbp->pturb->record_injected_work &&
       global_variable::my_rank == 0) ? pmbp->pturb->injected_work : 0.0;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::CGLLFPaper()
//! \brief Initialize the active-pressure MKS24 forced-turbulence configuration.

void ProblemGenerator::CGLLFPaper(ParameterInput *pin, const bool restart) {
  user_hist_func = CGLLFPaperHistory;

  auto *pmbp = pmy_mesh_->pmb_pack;
  auto *pmhd = pmbp->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl || pmhd->pcgl_lf == nullptr) {
    Fail("cgl_lf_paper requires <mhd>/eos = cgl and cgl_heat_flux = landau_fluid");
  }

  const std::string mode =
      pin->GetOrAddString("problem", "paper_mode", "turbulence");
  if (mode != "turbulence") {
    Fail("<problem>/paper_mode currently supports only turbulence");
  }
  const bool passive_delta =
      pin->GetOrAddBoolean("problem", "passive_delta", false);
  if (passive_delta != pmhd->peos->eos_data.passive) {
    Fail("<problem>/passive_delta must match <mhd>/passive for paper smoke cases");
  }
  if (restart) return;

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real b0z = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real beta0 = pin->GetOrAddReal("problem", "beta0", 10.0);
  if (rho0 <= 0.0 || beta0 <= 0.0 || b0z == 0.0) {
    Fail("cgl_lf_paper requires rho0 > 0, beta0 > 0, and b0 != 0");
  }
  const Real p0 = 0.5*beta0*SQR(b0z);
  const Real ppar0 = pin->GetOrAddReal("problem", "p_parallel0", p0);
  const Real pperp0 = pin->GetOrAddReal("problem", "p_perp0", p0);
  if (ppar0 <= 0.0 || pperp0 <= 0.0) {
    Fail("cgl_lf_paper requires positive initial parallel and perpendicular pressures");
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto w0 = pmhd->w0;
  auto bcc0 = pmhd->bcc0;
  auto b0 = pmhd->b0;

  par_for("cgl_lf_paper_prim", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w0(m,IDN,k,j,i) = rho0;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    w0(m,IPR,k,j,i) = ppar0;
    w0(m,IPP,k,j,i) = pperp0;
    bcc0(m,IBX,k,j,i) = 0.0;
    bcc0(m,IBY,k,j,i) = 0.0;
    bcc0(m,IBZ,k,j,i) = b0z;
  });
  par_for("cgl_lf_paper_b1", DevExeSpace(), 0, nmb - 1, ks, ke, js, je,
          is, ie + 1, KOKKOS_LAMBDA(const int m, const int k, const int j,
                                    const int i) {
    b0.x1f(m,k,j,i) = 0.0;
  });
  par_for("cgl_lf_paper_b2", DevExeSpace(), 0, nmb - 1, ks, ke, js, je + 1,
          is, ie, KOKKOS_LAMBDA(const int m, const int k, const int j,
                                const int i) {
    b0.x2f(m,k,j,i) = 0.0;
  });
  par_for("cgl_lf_paper_b3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je,
          is, ie, KOKKOS_LAMBDA(const int m, const int k, const int j,
                                const int i) {
    b0.x3f(m,k,j,i) = b0z;
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
}
