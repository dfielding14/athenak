//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
#include <cmath> // For cos and M_PI
#include <cstdlib> // For exit and EXIT_FAILURE

#include "athena.hpp"
#include "utils/random.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart)
//  \brief Problem Generator for Dynamo with Zero Net Magnetic Field

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // Enroll user history function
  user_hist_func = TurbulentHistory;

  // Exit early for restarts
  if (restart) return;

  // Fetch necessary parameters
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);
  Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
  Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
  Real B0 = cs * std::sqrt(2.0 * d_i / beta);
  Real p0 = 1.0 / (pin->GetOrAddReal("eos", "gamma", 5.0/3.0) - 1.0);

  // Initialize MHD variables
  // if (pmbp->pmhd != nullptr) {
  //   auto &u0 = pmbp->pmhd->u0;
  //   auto &b0 = pmbp->pmhd->b0;
  //   EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  //   Real gm1 = eos.gamma - 1.0;

  //   // Zero-Net Magnetic Field Initialization Parameters
  //   Real rseed = pin->GetReal("problem", "rseed", -1.0);
  //   bool first_time_ = true;
  //   int64_t iran = static_cast<int64_t>(rseed);
  //   int nrange = pin->GetOrAddInteger("problem", "nrange", 2);
  //   Real total_amp = 0.0;
  //   Real k_par = 2.0 * M_PI;
  //   Real b0_amp = pin->GetReal("problem", "b0_amp", B0); // Desired RMS magnetic field strength

  //   // Parallel loop to initialize magnetic fields with zero net flux
  //   par_for("pgen_turb", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  //   KOKKOS_LAMBDA(int m, int k, int j, int i) {
  //     // Initialize magnetic field components to zero
  //     b0.x1f(m, k, j, i) = 0.0;
  //     b0.x2f(m, k, j, i) = 0.0;
  //     b0.x3f(m, k, j, i) = 0.0;

  //     // Calculate positions
  //     Real x1 = pcoord->x1v(i);
  //     Real x2 = pcoord->x2v(j);
  //     Real x3 = pcoord->x3v(k);

  //     // Loop over wave numbers to superimpose cosine waves
  //     for (int nkx = -nrange; nkx <= nrange; nkx++) {
  //       for (int nky = -nrange; nky <= nrange; nky++) {
  //         for (int nkz = -nrange; nkz <= nrange; nkz++) {
  //           Real nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);

  //           if (nsqr > 0 && nsqr <= SQR(nrange)) {
  //             // Avoid double-counting by only initializing positive wave numbers
  //             if (nkx < 0 || (nkx == 0 && nky < 0) || (nkx == 0 && nky == 0 && nkz < 0)) {
  //               continue;
  //             }

  //             // Random amplitudes and phases
  //             Real ampx = (nkx != 0) ? 0.0 : RanGaussian(&iran);
  //             Real ampy = (nky != 0) ? 0.0 : RanGaussian(&iran);
  //             Real ampz = (nkz != 0) ? 0.0 : RanGaussian(&iran);

  //             Real ph1 = 2.0 * M_PI * Ran2(&iran);
  //             Real ph2 = 2.0 * M_PI * Ran2(&iran);
  //             Real ph3 = 2.0 * M_PI * Ran2(&iran);

  //             // Accumulate total amplitude for normalization
  //             total_amp += SQR(ampx) + SQR(ampy) + SQR(ampz);

  //             // Superimpose cosine waves
  //             b0.x1f(m, k, j, i) += ampx * cos(k_par * (nky * x2 + nkz * x3) + ph1);
  //             b0.x2f(m, k, j, i) += ampy * cos(k_par * (nkx * x1 + nkz * x3) + ph2);
  //             b0.x3f(m, k, j, i) += ampz * cos(k_par * (nkx * x1 + nky * x2) + ph3);
  //           }
  //         }
  //       }
  //     }

  //     // Apply boundary conditions for face-centered fields
  //     if (i == ie) {
  //       b0.x1f(m, k, j, i + 1) = b0.x1f(m, k, j, i);
  //     }
  //     if (j == je) {
  //       b0.x2f(m, k, je + 1, i) = b0.x2f(m, k, j, i);
  //     }
  //     if (k == ke) {
  //       b0.x3f(m, ke + 1, j, i) = b0.x3f(m, k, j, i);
  //     }

  //     // Update conserved variables if using ideal EOS
  //     if (eos.is_ideal) {
  //       u0(m, IEN, k, j, i) = p0 / gm1 + 0.5 * (SQR(u0(m, IM1, k, j, i)) +
  //                                               SQR(u0(m, IM2, k, j, i)) +
  //                                               SQR(u0(m, IM3, k, j, i))) / u0(m, IDN, k, j, i);
  //     }
  //   });

  //   // Normalize the magnetic field to achieve zero net flux
  //   Real Brms = b0_amp / sqrt(0.5 * total_amp);

  //   // Parallel loop to normalize all face-centered magnetic fields
  //   par_for("normalize_b0", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  //   KOKKOS_LAMBDA(int m, int k, int j, int i) {
  //     b0.x1f(m, k, j, i) *= Brms;
  //     b0.x2f(m, k, j, i) *= Brms;
  //     b0.x3f(m, k, j, i) *= Brms;

  //     // Update conserved variables after normalization
  //     if (eos.is_ideal) {
  //       u0(m, IEN, k, j, i) += 0.5 * (SQR(b0.x1f(m, k, j, i)) +
  //                                                SQR(b0.x2f(m, k, j, i)) +
  //                                                SQR(b0.x3f(m, k, j, i)));
  //     }
  //   });
  // }


  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {
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
      Real x1 = pcoord->x1f(i);
      b0.x3f(m,k,j,i) = B0 * std::tanh(16.0*x1);
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0;}

      if (eos.is_ideal) {
        Real x1v = pcoord->x1v(i);
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*SQR(B0 * std::tanh(16.0*x1v)) +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }


  return;
}

//----------------------------------------------------------------------------------------
// Function for computing history variables
// 0 = < B^4 >
// 1 = < (d_j B_i)(d_j B_i) >
// 2 = < (B_j d_j B_i)(B_k d_k B_i) >
// 3 = < |BxJ|^2 >
// 4 = < |B.J|^2 >
// 5 = < U^2 >
// 6 = < (d_j U_i)(d_j U_i) >
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
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

  // capture class variabels for kernel
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
    hvars.the_array[0] = bcc(m,IBX,k,j,i)*vol;
    hvars.the_array[1] = bcc(m,IBY,k,j,i)*vol;
    hvars.the_array[2] = bcc(m,IBZ,k,j,i)*vol;

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

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
  return;
}

</rewritten_chunk>