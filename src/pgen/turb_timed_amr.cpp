#include <cmath>
#include <iostream>

#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);

// Store eddy_time as static variable for use in RefinementCondition
static Real eddy_time = 1.0;

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for turbulence with time-dependent AMR

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // enroll user history function
  user_hist_func = TurbulentHistory;

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Read problem parameters
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);
  eddy_time = pin->GetOrAddReal("problem","eddy_time",1.0);

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = 1.0/eos.gamma;

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_n;
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
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    Real B0 = cs*std::sqrt(2.0*d_i/beta);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = 1.0/eos.gamma;

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = 1.0;
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
//! \fn int RefinementCondition()
//  \brief Time-dependent refinement condition for turbulence

int RefinementCondition(MeshBlockPack* pmbp) {
  // Get current simulation time
  Real time = pmbp->pmesh->time;

  // Define refinement schedule (in units of eddy_time)
  // Stage 1: t = 2.0 -> refine to level 1
  // Stage 2: t = 3.0 -> refine to level 2
  // Stage 3: t = 3.5 -> refine to level 3
  // Stage 4: t = 3.75 -> refine to level 4
  // Stage 5: t = 3.875 -> refine to level 5
  Real refine_times[] = {2.0, 3.0, 3.5, 3.75, 3.875};

  // Get current maximum refinement level in the mesh
  // We'll check the first meshblock's level and assume uniform refinement
  int current_level = 0;
  if (pmbp->nmb_thispack > 0) {
    current_level = pmbp->pmb->mb_lev.h_view(0);
  }

  // Determine which refinement stage we're at
  int target_level = 0;
  for (int i = 0; i < 5; i++) {
    if (time >= refine_times[i] * eddy_time - 0.001) {  // Small tolerance
      target_level = i + 1;
    }
  }

  // If we need to refine to reach the target level
  if (current_level < target_level) {
    // Return 1 to trigger refinement for all blocks
    return 1;
  }

  // No refinement needed
  return 0;
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

  pdata->nhist = 11;
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

  // capture class variables for kernel
  auto &bcc = pm->pmb_pack->pmhd->bcc0;
  auto &b = pm->pmb_pack->pmhd->b0;
  auto &w0_ = pm->pmb_pack->pmhd->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

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

    // Compute cell-centered B
    Real bccx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
    Real bccy = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
    Real bccz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
    Real b2 = SQR(bccx) + SQR(bccy) + SQR(bccz);

    // Compute velocity
    Real ux = w0_(m,IVX,k,j,i);
    Real uy = w0_(m,IVY,k,j,i);
    Real uz = w0_(m,IVZ,k,j,i);
    Real u2 = SQR(ux) + SQR(uy) + SQR(uz);

    // Basic quantities
    mb_sum.the_array[0] += bccx*vol;
    mb_sum.the_array[1] += bccy*vol;
    mb_sum.the_array[2] += bccz*vol;
    mb_sum.the_array[3] += b2*vol;
    mb_sum.the_array[4] += b2*b2*vol;
    mb_sum.the_array[9] += u2*vol;

    // Placeholder for gradient terms (would need proper implementation)
    mb_sum.the_array[5] += 0.0;  // dB^2
    mb_sum.the_array[6] += 0.0;  // BdB^2
    mb_sum.the_array[7] += 0.0;  // |BxJ|^2
    mb_sum.the_array[8] += 0.0;  // |B.J|^2
    mb_sum.the_array[10] += 0.0; // dU^2
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  // store data in HistoryData array
  for (int n=0; n<nhist_; n++) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}
