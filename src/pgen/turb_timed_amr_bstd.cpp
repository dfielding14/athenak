#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "globals.hpp"

// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);

// User-defined refinement condition
void RefinementCondition(MeshBlockPack* pmbp);

// Static variables to store refinement parameters
static Real t_refine = -1.0;
static Real bstd_refine = 0.5;
static Real bstd_derefine = -1.0;
static Real bstd_mean_floor = 1.0e-20;

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for turbulence with B-std AMR criterion

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // Read refinement controls (t_refine < 0 means active immediately)
  t_refine = pin->GetOrAddReal("problem","t_refine",-1.0);
  bstd_refine = pin->GetOrAddReal("problem","bstd_refine",0.5);
  bstd_derefine = pin->GetOrAddReal("problem","bstd_derefine",-1.0);
  bstd_mean_floor = pin->GetOrAddReal("problem","bstd_mean_floor",1.0e-20);

  // Register the refinement condition function
  user_ref_func = RefinementCondition;
  // enroll user history function
  user_hist_func = TurbulentHistory;

  if (restart) {
    return;
  }

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Read problem parameters
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    Real rho0 = pin->GetOrAddReal("problem","rho0",1.0);
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = rho0*cs*cs/eos.gamma;

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = rho0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {
    Real rho0 = pin->GetOrAddReal("problem","rho0",1.0);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = rho0*cs*cs/eos.gamma;
    Real B0 = cs*std::sqrt(2.0*p0/beta);


    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = rho0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      // initialize B
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = B0;
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0;}

      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*B0*B0 +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RefinementCondition()
//  \brief Refinement based on the per-block B std/mean ratio

void RefinementCondition(MeshBlockPack* pmbp) {
  Mesh *pmesh       = pmbp->pmesh;
  int nmb           = pmbp->nmb_thispack;
  int mbs           = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;

  if (pmbp->pmhd == nullptr) {
    return;
  }

  auto &indcs = pmesh->mb_indcs;
  int &is = indcs.is, nx1 = indcs.nx1;
  int &js = indcs.js, nx2 = indcs.nx2;
  int &ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3 * nx2 * nx1;
  const int nji  = nx2 * nx1;

  auto &bcc = pmbp->pmhd->bcc0;

  Real time = pmesh->time;
  Real t_ref = t_refine;
  Real refine_thresh = bstd_refine;
  Real derefine_thresh = bstd_derefine;
  Real mean_floor = bstd_mean_floor;

  if (t_ref >= 0.0 && time < t_ref) {
    return;
  }

  par_for_outer("UserRefineCondBstd",DevExeSpace(), 0, 0, 0, (nmb-1),
  KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {
    Real sum_b = 0.0;
    Real sum_b2 = 0.0;

    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
    [=](const int idx, Real &local_sum) {
      int k = (idx)/nji;
      int j = (idx - k*nji)/nx1;
      int i = (idx - k*nji - j*nx1) + is;
      j += js;
      k += ks;

      Real bx = bcc(m,IBX,k,j,i);
      Real by = bcc(m,IBY,k,j,i);
      Real bz = bcc(m,IBZ,k,j,i);
      Real bmag = sqrt(bx*bx + by*by + bz*bz);
      local_sum += bmag;
    },Kokkos::Sum<Real>(sum_b));

    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
    [=](const int idx, Real &local_sum) {
      int k = (idx)/nji;
      int j = (idx - k*nji)/nx1;
      int i = (idx - k*nji - j*nx1) + is;
      j += js;
      k += ks;

      Real bx = bcc(m,IBX,k,j,i);
      Real by = bcc(m,IBY,k,j,i);
      Real bz = bcc(m,IBZ,k,j,i);
      Real bmag = sqrt(bx*bx + by*by + bz*bz);
      local_sum += bmag*bmag;
    },Kokkos::Sum<Real>(sum_b2));

    Real inv_ncells = 1.0/static_cast<Real>(nkji);
    Real bmean = sum_b * inv_ncells;
    Real bvar = sum_b2 * inv_ncells - bmean*bmean;
    if (bvar < 0.0) bvar = 0.0;
    Real bstd = sqrt(bvar);
    Real denom = bmean + mean_floor;
    Real ratio = (denom > 0.0) ? (bstd/denom) : 0.0;

    if (refine_thresh > 0.0 && ratio > refine_thresh) {
      refine_flag.d_view(m+mbs) = 1;
    } else if (derefine_thresh >= 0.0 && ratio < derefine_thresh) {
      refine_flag.d_view(m+mbs) = -1;
    }
  });

  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

//----------------------------------------------------------------------------------------
// Function for computing history variables (same as turb.cpp)
// 0-2 = <Bx>, <By>, <Bz>
// 3 = <B^2>
// 4 = <B^4>
// 5 = <(d_j B_i)(d_j B_i)>
// 6 = <(B_j d_j B_i)(B_k d_k B_i)>
// 7 = <|BxJ|^2>
// 8 = <|B.J|^2>
// 9 = <U^2>
// 10 = <(d_j U_i)(d_j U_i)>

void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
  // Only compute if MHD is enabled
  if (pm->pmb_pack->pmhd == nullptr) {
    pdata->nhist = 0;
    return;
  }

  pdata->nhist = 12;
  pdata->label[0] = "Bx";
  pdata->label[1] = "By";
  pdata->label[2] = "Bz";
  pdata->label[3] = "B^2";
  pdata->label[4] = "B^4";
  pdata->label[5] = "dB^2";
  pdata->label[6] = "BdB^2";
  pdata->label[7] = "|BxJ|^2";
  pdata->label[8] = "|B.J|^2";
  pdata->label[9] = "U^2";
  pdata->label[10] = "dU";
  pdata->label[11] = "|divB|";

  // capture class variables for kernel
  auto &bcc = pm->pmb_pack->pmhd->bcc0;
  auto &b = pm->pmb_pack->pmhd->b0;
  auto &w0_ = pm->pmb_pack->pmhd->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;
  bool &multi_d = pm->multi_d;
  bool &three_d = pm->three_d;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Real total_vol = 0.0;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real dx_squared = size.d_view(m).dx1 * size.d_view(m).dx1;

    // MHD conserved variables:
    array_sum::GlobalSum hvars;

    // calculate mean B
    hvars.the_array[0] = bcc(m,IBX,k,j,i);
    hvars.the_array[1] = bcc(m,IBY,k,j,i);
    hvars.the_array[2] = bcc(m,IBZ,k,j,i);

    // 0 = < B^2 >
    Real B_mag_sq = bcc(m,IBX,k,j,i)*bcc(m,IBX,k,j,i)
                  + bcc(m,IBY,k,j,i)*bcc(m,IBY,k,j,i)
                  + bcc(m,IBZ,k,j,i)*bcc(m,IBZ,k,j,i);
    hvars.the_array[3] = B_mag_sq*vol;
    // 0 = < B^4 >
    Real B_fourth = B_mag_sq*B_mag_sq;
    hvars.the_array[4] = B_fourth*vol;
    // 1 = < (d_j B_i)(d_j B_i) >
    hvars.the_array[5] = (
      ((b.x1f(m,k,j,i+1)-b.x1f(m,k,j,i))*(b.x1f(m,k,j,i+1)-b.x1f(m,k,j,i))
     + (b.x2f(m,k,j+1,i)-b.x2f(m,k,j,i))*(b.x2f(m,k,j+1,i)-b.x2f(m,k,j,i))
     + (b.x3f(m,k+1,j,i)-b.x3f(m,k,j,i))*(b.x3f(m,k+1,j,i)-b.x3f(m,k,j,i))
     + 0.25*(bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))
           *(bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))
     + 0.25*(bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))
           *(bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))
     + 0.25*(bcc(m,IBY,k,j,i+1)-bcc(m,IBY,k,j,i-1))
           *(bcc(m,IBY,k,j,i+1)-bcc(m,IBY,k,j,i-1))
     + 0.25*(bcc(m,IBY,k+1,j,i)-bcc(m,IBY,k-1,j,i))
           *(bcc(m,IBY,k+1,j,i)-bcc(m,IBY,k-1,j,i))
     + 0.25*(bcc(m,IBZ,k,j,i+1)-bcc(m,IBZ,k,j,i-1))
           *(bcc(m,IBZ,k,j,i+1)-bcc(m,IBZ,k,j,i-1))
     + 0.25*(bcc(m,IBZ,k,j+1,i)-bcc(m,IBZ,i,j-1,i))
           *(bcc(m,IBZ,k,j+1,i)-bcc(m,IBZ,i,j-1,i)))
       / dx_squared)*vol;
    // 2 = < (B_j d_j B_i)(B_k d_k B_i) >
    Real bdb1 = bcc(m,IBX,k,j,i)*(b.x1f(m,k,j,i+1)-b.x1f(m,k,j,i))
                +0.5*bcc(m,IBY,k,j,i)*(bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))
                +0.5*bcc(m,IBZ,k,j,i)*(bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i));
    Real bdb2 = bcc(m,IBY,k,j,i)*(b.x2f(m,k,j+1,i)-b.x2f(m,k,j,i))
                +0.5*bcc(m,IBZ,k,j,i)*(bcc(m,IBY,k+1,j,i)-bcc(m,IBY,k-1,j,i))
                +0.5*bcc(m,IBX,k,j,i)*(bcc(m,IBY,k,j,i+1)-bcc(m,IBY,k,j,i-1));
    Real bdb3 = bcc(m,IBZ,k,j,i)*(b.x3f(m,k+1,j,i)-b.x3f(m,k,j,i))
                +0.5*bcc(m,IBX,k,j,i)*(bcc(m,IBZ,k,j,i+1)-bcc(m,IBZ,k,j,i-1))
                +0.5*bcc(m,IBY,k,j,i)*(bcc(m,IBZ,k,j+1,i)-bcc(m,IBZ,k,j-1,i));
    hvars.the_array[6] = ((bdb1*bdb1 + bdb2*bdb2 + bdb3*bdb3) / dx_squared)*vol;
    // 3 = < |BxJ|^2 >
    Real Jx = 0.5*(bcc(m,IBZ,k,j+1,i)-bcc(m,IBZ,k,j-1,i))
             -0.5*(bcc(m,IBY,k+1,j,i)-bcc(m,IBY,k-1,j,i));
    Real Jy = 0.5*(bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))
             -0.5*(bcc(m,IBZ,k,j,i+1)-bcc(m,IBZ,k,j,i-1));
    Real Jz = 0.5*(bcc(m,IBY,k,j,i+1)-bcc(m,IBY,k,j,i-1))
             -0.5*(bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i));
    hvars.the_array[7] =((
       (bcc(m,IBY,k,j,i)*Jz - bcc(m,IBZ,k,j,i)*Jy)
      *(bcc(m,IBY,k,j,i)*Jz - bcc(m,IBZ,k,j,i)*Jy)
      +(bcc(m,IBZ,k,j,i)*Jx - bcc(m,IBX,k,j,i)*Jz)
      *(bcc(m,IBZ,k,j,i)*Jx - bcc(m,IBX,k,j,i)*Jz)
      +(bcc(m,IBX,k,j,i)*Jy - bcc(m,IBY,k,j,i)*Jx)
      *(bcc(m,IBX,k,j,i)*Jy - bcc(m,IBY,k,j,i)*Jx))
                    / dx_squared)*vol;
    // 4 = < |B.J|^2 >
    hvars.the_array[8] = (
      ((bcc(m,IBX,k,j,i)*Jx + bcc(m,IBY,k,j,i)*Jy + bcc(m,IBZ,k,j,i)*Jz)
      *(bcc(m,IBX,k,j,i)*Jx + bcc(m,IBY,k,j,i)*Jy + bcc(m,IBZ,k,j,i)*Jz)
                          )/dx_squared)*vol;
    // 5 = < U^2 >
    hvars.the_array[9] += ((w0_(m,IVX,k,j,i)*w0_(m,IVX,k,j,i))
                        + (w0_(m,IVY,k,j,i)*w0_(m,IVY,k,j,i))
                        + (w0_(m,IVZ,k,j,i)*w0_(m,IVZ,k,j,i)))*vol;
    // 6 = < (d_j U_i)(d_j U_i) >
    hvars.the_array[10] +=
    (((0.25*(w0_(m,IVX,k,j,i+1)-w0_(m,IVX,k,j,i-1))
           *(w0_(m,IVX,k,j,i+1)-w0_(m,IVX,k,j,i-1))
     + 0.25*(w0_(m,IVY,k,j+1,i)-w0_(m,IVY,k,j-1,i))
           *(w0_(m,IVY,k,j+1,i)-w0_(m,IVY,k,j-1,i))
     + 0.25*(w0_(m,IVZ,k+1,j,i)-w0_(m,IVZ,k-1,j,i))
           *(w0_(m,IVZ,k+1,j,i)-w0_(m,IVZ,k-1,j,i))
     + 0.25*(w0_(m,IVX,k,j+1,i)-w0_(m,IVX,k,j-1,i))
           *(w0_(m,IVX,k,j+1,i)-w0_(m,IVX,k,j-1,i))
     + 0.25*(w0_(m,IVX,k+1,j,i)-w0_(m,IVX,k-1,j,i))
           *(w0_(m,IVX,k+1,j,i)-w0_(m,IVX,k-1,j,i))
     + 0.25*(w0_(m,IVY,k,j,i+1)-w0_(m,IVY,k,j,i-1))
           *(w0_(m,IVY,k,j,i+1)-w0_(m,IVY,k,j,i-1))
     + 0.25*(w0_(m,IVY,k+1,j,i)-w0_(m,IVY,k-1,j,i))
           *(w0_(m,IVY,k+1,j,i)-w0_(m,IVY,k-1,j,i))
     + 0.25*(w0_(m,IVZ,k,j,i+1)-w0_(m,IVZ,k,j,i-1))
           *(w0_(m,IVZ,k,j,i+1)-w0_(m,IVZ,k,j,i-1))
     + 0.25*(w0_(m,IVZ,k,j+1,i)-w0_(m,IVZ,k,j-1,i))
           *(w0_(m,IVZ,k,j+1,i)-w0_(m,IVZ,k,j-1,i))))
     / dx_squared)*vol;

    // Compute divergence of B
    Real divb = (b.x1f(m,k,j,i+1) - b.x1f(m,k,j,i))/size.d_view(m).dx1;
    if (multi_d) {
      divb += (b.x2f(m,k,j+1,i) - b.x2f(m,k,j,i))/size.d_view(m).dx2;
    }
    if (three_d) {
      divb += (b.x3f(m,k+1,j,i) - b.x3f(m,k,j,i))/size.d_view(m).dx3;
    }

    // Sum of divergence of B
    hvars.the_array[11] += fabs(divb)*vol;
    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  // Kokkos::fence();

  // store data in HistoryData array
  for (int n=0; n<nhist_; n++) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}
