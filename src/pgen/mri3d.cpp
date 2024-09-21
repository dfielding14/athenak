//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mri3d.cpp
//! \brief Problem generator for 3D MRI.
//!
//! PURPOSE:  Problem generator for 3D MRI. Based on the initial conditions described in
//! "Local Three-dimensional Magnetohydrodynamic Simulations of Accretion Disks" by
//! Hawley, Gammie & Balbus, or HGB.  AthanK version based on pgen/hgb.cpp in Athena++.
//!
//! Several different field configurations are possible:
//! - ifield = 1 - Bz=B0 sin(nwx*kx*x1) field with zero-net-flux [default] (nwx input)
//! - ifield = 2 - uniform Bz
//! - ifield = 3 - uniform By
//! Random perturbations to the pressure are added in the initial conditions to seed MRI
//!
//! REFERENCE: Hawley, J. F., Gammie, C.F. & Balbus, S. A., ApJ 440, 742-763 (1995).

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "shearing_box/shearing_box.hpp"
#include "pgen.hpp"

#include <Kokkos_Random.hpp>

// User-defined history function
void MRIHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::_()
//  \brief

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // enroll user history function
  user_hist_func = MRIHistory;
  if (restart) return;

  // First, do some error checks
  if (!(pmy_mesh_->three_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "mri3d problem generator only works in 2D (nx3=1)" << std::endl;
    exit(EXIT_FAILURE);
  }
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd != nullptr) {
    if (pmbp->pmhd->psrc == nullptr) {
      std::cout <<"### FATAL ERROR in "<< __FILE__ <<" at line " <<__LINE__ << std::endl
                << "Shearing box source terms not enabled for mri3d problem" << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "mri3d problem generator only works in mhd" << std::endl;
    exit(EXIT_FAILURE);
  }

  // initialize problem variables
  Real amp   = pin->GetReal("problem","amp");
  Real beta  = pin->GetReal("problem","beta");
  int nwx    = pin->GetOrAddInteger("problem","nwx",1);
  int ifield = pin->GetOrAddInteger("problem","ifield",1);

  // background density, pressure, and magnetic field
  Real d0 = 1.0;
  Real p0 = 1.0;
  Real binit = std::sqrt(2.0*p0/beta);

  Real x1size = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real kx = 2.0*(M_PI/x1size)*(static_cast<Real>(nwx));

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;

  // Initialize magnetic field first, so entire arrays are initialized before adding
  // magnetic energy to conserved variables in next loop.  For 3D shearing box
  // B1=Bx, B2=By, B3=Bz
  // ifield = 1 - Bz=binit sin(kx*xav1) field with zero-net-flux [default]
  // ifield = 2 - uniform Bz
  // ifield = 3 - uniform By
  auto b0 = pmbp->pmhd->b0;
  par_for("mri3d", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    if (ifield == 1) {
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = binit*sin(kx*x1v);
      if (i==ie) b0.x1f(m,k,j,ie+1) = 0.0;
      if (j==je) b0.x2f(m,k,je+1,i) = 0.0;
      if (k==ke) b0.x3f(m,ke+1,j,i) = binit*sin(kx*x1v);
    } else if (ifield == 2) {
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = binit;
      if (i==ie) b0.x1f(m,k,j,ie+1) = 0.0;
      if (j==je) b0.x2f(m,k,je+1,i) = 0.0;
      if (k==ke) b0.x3f(m,ke+1,j,i) = binit;
    } else if (ifield == 3) {
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = binit;
      b0.x3f(m,k,j,i) = 0.0;
      if (i==ie) b0.x1f(m,k,j,ie+1) = 0.0;
      if (j==je) b0.x2f(m,k,je+1,i) = binit;
      if (k==ke) b0.x3f(m,ke+1,j,i) = 0.0;
    }
  });

  // Initialize conserved variables
  // Only sets up random perturbations in pressure to seed MRI
  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  auto u0 = pmbp->pmhd->u0;
  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
  par_for("mri3d-u", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real rd = d0;
    Real rp = p0;
    auto rand_gen = rand_pool64.get_state();  // get random number state this thread
    Real rval = amp*(rand_gen.frand() - 0.5);
    if (eos.is_ideal) {
      rp = p0*(1.0 + rval);;
    } else {
      rd = d0*(1.0 + rval);
    }
    u0(m,IDN,k,j,i) = rd;
    rval = amp*(rand_gen.frand() - 0.5);
    u0(m,IM1,k,j,i) = d0*rval;
    rval = amp*(rand_gen.frand() - 0.5);
    u0(m,IM2,k,j,i) = d0*rval;
    rval = amp*(rand_gen.frand() - 0.5);
    u0(m,IM3,k,j,i) = d0*rval;
    if (eos.is_ideal) {
      u0(m,IEN,k,j,i) = rp/gm1 + 0.5*SQR(0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i))) +
                                 0.5*SQR(0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i)));
    }
    rand_pool64.free_state(rand_gen);  // free state for use by other threads
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MRIHistory()
//  \brief Compute and store MRI history data.  Adds Reynolds and Maxwell stress and net
//  magnetic flux to usual list of MHD history variables

void MRIHistory(HistoryData *pdata, Mesh *pm) {
  auto &eos_data = pm->pmb_pack->pmhd->peos->eos_data;
  int &nmhd_ = pm->pmb_pack->pmhd->nmhd;

  // set number of and names of history variables for mhd
  if (eos_data.is_ideal) {
    pdata->nhist = 16;
  } else {
    pdata->nhist = 15;
  }
  pdata->label[IDN] = "mass";
  pdata->label[IM1] = "1-mom";
  pdata->label[IM2] = "2-mom";
  pdata->label[IM3] = "3-mom";
  if (eos_data.is_ideal) {
    pdata->label[IEN] = "tot-E";
  }
  pdata->label[nmhd_  ] = "1-KE";
  pdata->label[nmhd_+1] = "2-KE";
  pdata->label[nmhd_+2] = "3-KE";
  pdata->label[nmhd_+3] = "1-ME";
  pdata->label[nmhd_+4] = "2-ME";
  pdata->label[nmhd_+5] = "3-ME";
  pdata->label[nmhd_+6] = "1-bcc";
  pdata->label[nmhd_+7] = "2-bcc";
  pdata->label[nmhd_+8] = "3-bcc";
  pdata->label[nmhd_+9] = "dVxVy";
  pdata->label[nmhd_+10] = "dBxBy";

  // capture class variabels for kernel
  auto &u0_ = pm->pmb_pack->pmhd->u0;
  auto &bx1f = pm->pmb_pack->pmhd->b0.x1f;
  auto &bx2f = pm->pmb_pack->pmhd->b0.x2f;
  auto &bx3f = pm->pmb_pack->pmhd->b0.x3f;
  auto &bcc = pm->pmb_pack->pmhd->bcc0;
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

    // MHD conserved variables:
    array_sum::GlobalSum hvars;
    hvars.the_array[IDN] = vol*u0_(m,IDN,k,j,i);
    hvars.the_array[IM1] = vol*u0_(m,IM1,k,j,i);
    hvars.the_array[IM2] = vol*u0_(m,IM2,k,j,i);
    hvars.the_array[IM3] = vol*u0_(m,IM3,k,j,i);
    if (eos_data.is_ideal) {
      hvars.the_array[IEN] = vol*u0_(m,IEN,k,j,i);
    }

    // MHD KE
    hvars.the_array[nmhd_  ] = vol*0.5*SQR(u0_(m,IM1,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nmhd_+1] = vol*0.5*SQR(u0_(m,IM2,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nmhd_+2] = vol*0.5*SQR(u0_(m,IM3,k,j,i))/u0_(m,IDN,k,j,i);

    // MHD ME
    hvars.the_array[nmhd_+3] = vol*0.25*(SQR(bx1f(m,k,j,i+1)) + SQR(bx1f(m,k,j,i)));
    hvars.the_array[nmhd_+4] = vol*0.25*(SQR(bx2f(m,k,j+1,i)) + SQR(bx2f(m,k,j,i)));
    hvars.the_array[nmhd_+5] = vol*0.25*(SQR(bx3f(m,k+1,j,i)) + SQR(bx3f(m,k,j,i)));

    // net B fluxes
    hvars.the_array[nmhd_+6] = vol*bcc(m,IBX,k,j,i);
    hvars.the_array[nmhd_+7] = vol*bcc(m,IBY,k,j,i);
    hvars.the_array[nmhd_+8] = vol*bcc(m,IBZ,k,j,i);

    // Reynolds and Maxwell stresses
    hvars.the_array[nmhd_+9] = vol*u0_(m,IM1,k,j,i)*u0_(m,IM2,k,j,i)/u0_(m,IDN,k,j,i);
    hvars.the_array[nmhd_+10] = -vol*bcc(m,IBX,k,j,i)*bcc(m,IBY,k,j,i);

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
  return;
}