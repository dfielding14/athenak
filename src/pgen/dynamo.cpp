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

  // Things to do only for sims starting from ICs
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with MHD, but no "
       << "<mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Fetch necessary parameters
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);
  Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
  Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
  Real B0 = cs * std::sqrt(2.0 * d_i / beta);
  Real p0 = 1.0 / (pin->GetOrAddReal("eos", "gamma", 5.0/3.0) - 1.0);

  auto &size = pmbp->pmb->mb_size;

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
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      b0.x3f(m,k,j,i) = B0 * std::tanh(16.0*x1v);
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0 * std::tanh(16.0*x1v);}

      if (eos.is_ideal) {
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
// < Bx >
// < By >
// < Bz >
// < B^2 >
// < B^4 >
// < (d_j B_i)(d_j B_i) >
// < (B_j d_j B_i)(B_k d_k B_i) >
// < |BxJ|^2 >
// < |B.J|^2 >
// < U^2 >
// < (d_j U_i)(d_j U_i) >
// < |J|^2 >
// < Pi:dv >
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 14;
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
  pdata->label[10] = "|dU|^2";
  pdata->label[11] = "|J|^2";
  pdata->label[12] = "Pi:dv";

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
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;
    Real dx1_squared = dx1 * dx1;
    Real dx2_squared = dx2 * dx2;
    Real dx3_squared = dx3 * dx3;

    // gradient components
    // face-centered
    Real db1f_dx1 = (b.x1f(m,k,j,i+1)-b.x1f(m,k,j,i));
    Real db1f_dx2 = (b.x1f(m,k,j+1,i)-b.x1f(m,k,j,i));
    Real db1f_dx3 = (b.x1f(m,k+1,j,i)-b.x1f(m,k,j,i));
    Real db2f_dx1 = (b.x2f(m,k,j,i+1)-b.x2f(m,k,j,i));
    Real db2f_dx2 = (b.x2f(m,k,j+1,i)-b.x2f(m,k,j,i));
    Real db2f_dx3 = (b.x2f(m,k+1,j,i)-b.x2f(m,k,j,i));
    Real db3f_dx1 = (b.x3f(m,k,j,i+1)-b.x3f(m,k,j,i));
    Real db3f_dx2 = (b.x3f(m,k,j+1,i)-b.x3f(m,k,j,i));
    Real db3f_dx3 = (b.x3f(m,k+1,j,i)-b.x3f(m,k,j,i));
    // cell-centered
    Real db1c_dx1 = (bcc(m,IBX,k,j,i+1)-bcc(m,IBX,k,j,i-1))*0.5;
    Real db1c_dx2 = (bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))*0.5;
    Real db1c_dx3 = (bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))*0.5;
    Real db2c_dx1 = (bcc(m,IBX,k,j,i+1)-bcc(m,IBX,k,j,i-1))*0.5;
    Real db2c_dx2 = (bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))*0.5;
    Real db2c_dx3 = (bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))*0.5;
    Real db3c_dx1 = (bcc(m,IBX,k,j,i+1)-bcc(m,IBX,k,j,i-1))*0.5;
    Real db3c_dx2 = (bcc(m,IBX,k,j+1,i)-bcc(m,IBX,k,j-1,i))*0.5;
    Real db3c_dx3 = (bcc(m,IBX,k+1,j,i)-bcc(m,IBX,k-1,j,i))*0.5;
    // gradient terms
    Real db1_dx1 = (4/3.)*db1f_dx1 - (1/6.)*db1c_dx1;
    Real db2_dx1 = (4/3.)*db2f_dx1 - (1/6.)*db2c_dx1;
    Real db3_dx1 = (4/3.)*db3f_dx1 - (1/6.)*db3c_dx1;
    Real db1_dx2 = (4/3.)*db1f_dx2 - (1/6.)*db1c_dx2;
    Real db2_dx2 = (4/3.)*db2f_dx2 - (1/6.)*db2c_dx2;
    Real db3_dx2 = (4/3.)*db3f_dx2 - (1/6.)*db3c_dx2;
    Real db1_dx3 = (4/3.)*db1f_dx3 - (1/6.)*db1c_dx3;
    Real db2_dx3 = (4/3.)*db2f_dx3 - (1/6.)*db2c_dx3;
    Real db3_dx3 = (4/3.)*db3f_dx3 - (1/6.)*db3c_dx3;
    // central values
    Real b1 = bcc(m,IBX,k,j,i);
    Real b2 = bcc(m,IBY,k,j,i);
    Real b3 = bcc(m,IBZ,k,j,i);

    // MHD conserved variables:
    array_sum::GlobalSum hvars;

    // calculate mean B
    hvars.the_array[0] += b1*vol;
    hvars.the_array[1] += b2*vol;
    hvars.the_array[2] += b3*vol;

    // < B^2 >
    Real B_mag_sq = b1*b1 + b2*b2 + b3*b3;
    hvars.the_array[3] += B_mag_sq*vol;

    // < B^4 >
    hvars.the_array[4] += B_mag_sq*B_mag_sq*vol;

    // < (d_j B_i)(d_j B_i) >
    hvars.the_array[5] += ((SQR(db1_dx1) + SQR(db2_dx1) + SQR(db3_dx1))/dx1_squared
                          +(SQR(db1_dx2) + SQR(db2_dx2) + SQR(db3_dx2))/dx2_squared
                          +(SQR(db1_dx3) + SQR(db2_dx3) + SQR(db3_dx3))/dx3_squared)*vol;

    // < (B_j d_j B_i)(B_k d_k B_i) >
    Real bdb1 = b1*db1_dx1 + b2*db1_dx2 + b3*db1_dx3;
    Real bdb2 = b1*db2_dx1 + b2*db2_dx2 + b3*db2_dx3;
    Real bdb3 = b1*db3_dx1 + b2*db3_dx2 + b3*db3_dx3;
    hvars.the_array[6] += ( SQR(bdb1)/dx1_squared
                          + SQR(bdb2)/dx2_squared
                          + SQR(bdb3)/dx3_squared)*vol;

    // < |BxJ|^2 >
    Real Jx = (db3_dx2/dx2 - db2_dx3/dx3);
    Real Jy = (db1_dx3/dx3 - db3_dx1/dx1);
    Real Jz = (db2_dx1/dx1 - db1_dx2/dx2);
    hvars.the_array[7] += ( SQR(b2*Jz - b3*Jy)
                           +SQR(b3*Jx - b1*Jz)
                           +SQR(b1*Jy - b2*Jx))*vol;

    // < |B.J|^2 >
    hvars.the_array[8] = SQR(b1*Jx + b2*Jy + b3*Jz)*vol;

    // < U^2 >
    hvars.the_array[9] += ((w0_(m,IVX,k,j,i)*w0_(m,IVX,k,j,i))
                         + (w0_(m,IVY,k,j,i)*w0_(m,IVY,k,j,i))
                         + (w0_(m,IVZ,k,j,i)*w0_(m,IVZ,k,j,i)))*vol;
    // < (d_j U_i)(d_j U_i) >
    Real dvx_dx1 = (w0_(m,IVX,k,j,i+1)-w0_(m,IVX,k,j,i-1))*0.5;
    Real dvx_dx2 = (w0_(m,IVX,k,j+1,i)-w0_(m,IVX,k,j-1,i))*0.5;
    Real dvx_dx3 = (w0_(m,IVX,k+1,j,i)-w0_(m,IVX,k-1,j,i))*0.5;
    Real dvy_dx1 = (w0_(m,IVY,k,j,i+1)-w0_(m,IVY,k,j,i-1))*0.5;
    Real dvy_dx2 = (w0_(m,IVY,k,j+1,i)-w0_(m,IVY,k,j-1,i))*0.5;
    Real dvy_dx3 = (w0_(m,IVY,k+1,j,i)-w0_(m,IVY,k-1,j,i))*0.5;
    Real dvz_dx1 = (w0_(m,IVZ,k,j,i+1)-w0_(m,IVZ,k,j,i-1))*0.5;
    Real dvz_dx2 = (w0_(m,IVZ,k,j+1,i)-w0_(m,IVZ,k,j-1,i))*0.5;
    Real dvz_dx3 = (w0_(m,IVZ,k+1,j,i)-w0_(m,IVZ,k-1,j,i))*0.5;
    hvars.the_array[10] +=((SQR(dvx_dx1) + SQR(dvy_dx1) + SQR(dvz_dx1))/dx1_squared
                          +(SQR(dvx_dx2) + SQR(dvy_dx2) + SQR(dvz_dx2))/dx2_squared
                          +(SQR(dvx_dx3) + SQR(dvy_dx3) + SQR(dvz_dx3))/dx3_squared)*vol;

    // < |J|^2 >
    hvars.the_array[11] += SQR(Jx) + SQR(Jy) + SQR(Jz)*vol;

    // <Pi_ij d_i v_j>
    Real Pi_xx = (2*dvx_dx1 - (2/3.)*(dvx_dx1 + dvy_dx2 + dvz_dx3));
    Real Pi_yy = (2*dvy_dx2 - (2/3.)*(dvx_dx1 + dvy_dx2 + dvz_dx3));
    Real Pi_zz = (2*dvz_dx3 - (2/3.)*(dvx_dx1 + dvy_dx2 + dvz_dx3));
    Real Pi_xy = (dvx_dx2 + dvy_dx1);
    Real Pi_xz = (dvx_dx3 + dvz_dx1);
    Real Pi_yz = (dvy_dx3 + dvz_dx2);
    hvars.the_array[12] += (Pi_xx*dvx_dx1 + Pi_yy*dvy_dx2 + Pi_zz*dvz_dx3
                           + Pi_xy*(dvx_dx2 + dvy_dx1)
                           + Pi_xz*(dvx_dx3 + dvz_dx1)
                           + Pi_yz*(dvy_dx3 + dvz_dx2))*vol*w0_(m,IDN,k,j,i);
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
