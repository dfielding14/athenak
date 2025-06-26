//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sfb_turb.cpp
//  \brief Problem generator for turbulence with Spherical Fourier-Bessel driving
#include <iostream> // cout

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
void SFBTurbulentHistory(HistoryData *pdata, Mesh *pm);


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::SFBTurb_()
//  \brief Problem Generator for spherical Fourier-Bessel turbulence

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "SFB Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // enroll user history function
  user_hist_func = SFBTurbulentHistory;

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Read problem parameters
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  Real beta = pin->GetOrAddReal("problem","beta",1.0);
  Real rho0 = pin->GetOrAddReal("problem","rho0",1.0);
  Real p0 = pin->GetOrAddReal("problem","p0",1.0);
  
  // Radial profile parameters for initial density/pressure
  Real r_core = pin->GetOrAddReal("problem","r_core",0.1);
  Real r_outer = pin->GetOrAddReal("problem","r_outer",1.0);
  Real rho_outer = pin->GetOrAddReal("problem","rho_outer",0.1);
  
  // Center of the sphere
  Real x_center = pin->GetOrAddReal("problem","x_center",0.0);
  Real y_center = pin->GetOrAddReal("problem","y_center",0.0);
  Real z_center = pin->GetOrAddReal("problem","z_center",0.0);

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;

    // Set initial conditions with radial profile
    par_for("pgen_sfb_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = pmbp->pmb->mb_size.d_view(m).x1min;
      Real &x1max = pmbp->pmb->mb_size.d_view(m).x1max;
      Real &x2min = pmbp->pmb->mb_size.d_view(m).x2min;
      Real &x2max = pmbp->pmb->mb_size.d_view(m).x2max;
      Real &x3min = pmbp->pmb->mb_size.d_view(m).x3min;
      Real &x3max = pmbp->pmb->mb_size.d_view(m).x3max;
      
      Real x = CellCenterX(i-is, indcs.nx1, x1min, x1max) - x_center;
      Real y = CellCenterX(j-js, indcs.nx2, x2min, x2max) - y_center;
      Real z = CellCenterX(k-ks, indcs.nx3, x3min, x3max) - z_center;
      Real r = sqrt(x*x + y*y + z*z);
      
      // Smooth density profile
      Real density;
      if (r < r_core) {
        density = rho0;
      } else if (r < r_outer) {
        // Smooth transition
        Real t = (r - r_core) / (r_outer - r_core);
        density = rho0 * (1.0 - t) + rho_outer * t;
      } else {
        density = rho_outer;
      }
      
      u0(m,IDN,k,j,i) = density;
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
    Real B0 = cs*std::sqrt(2.0*rho0/beta);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;

    // Set initial conditions with radial profile
    par_for("pgen_sfb_turb_mhd", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = pmbp->pmb->mb_size.d_view(m).x1min;
      Real &x1max = pmbp->pmb->mb_size.d_view(m).x1max;
      Real &x2min = pmbp->pmb->mb_size.d_view(m).x2min;
      Real &x2max = pmbp->pmb->mb_size.d_view(m).x2max;
      Real &x3min = pmbp->pmb->mb_size.d_view(m).x3min;
      Real &x3max = pmbp->pmb->mb_size.d_view(m).x3max;
      
      Real x = CellCenterX(i-is, indcs.nx1, x1min, x1max) - x_center;
      Real y = CellCenterX(j-js, indcs.nx2, x2min, x2max) - y_center;
      Real z = CellCenterX(k-ks, indcs.nx3, x3min, x3max) - z_center;
      Real r = sqrt(x*x + y*y + z*z);
      
      // Smooth density profile
      Real density;
      if (r < r_core) {
        density = rho0;
      } else if (r < r_outer) {
        // Smooth transition
        Real t = (r - r_core) / (r_outer - r_core);
        density = rho0 * (1.0 - t) + rho_outer * t;
      } else {
        density = rho_outer;
      }
      
      u0(m,IDN,k,j,i) = density;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      // initialize B - weak uniform field
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
//! \fn void SFBTurbulentHistory()
//  \brief User-defined history function for SFB turbulence

void SFBTurbulentHistory(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is; int ie = indcs.ie;
  int js = indcs.js; int je = indcs.je;
  int ks = indcs.ks; int ke = indcs.ke;
  
  Real vol = 0.0, dEk = 0.0, Ek = 0.0, ME = 0.0;
  Real drho_rms = 0.0, mean_rho = 0.0;

  if (pmbp->phydro != nullptr) {
    auto &u0 = pmbp->phydro->u0;
    const int nmkji = (pmbp->nmb_thispack)*indcs.nx3*indcs.nx2*indcs.nx1;
    const int nkji = indcs.nx3*indcs.nx2*indcs.nx1;
    const int nji  = indcs.nx2*indcs.nx1;
    
    Kokkos::parallel_reduce("SFB_turb_hist", Kokkos::RangePolicy<>(DevExeSpace(),0,nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &lvol, Real &lEk, Real &lmrho, Real &ldrho2) {
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/indcs.nx1;
      int i = (idx - m*nkji - k*nji - j*indcs.nx1) + is;
      k += ks; j += js;
      
      Real dx1 = pmbp->pmb->mb_size.d_view(m).dx1;
      Real dx2 = pmbp->pmb->mb_size.d_view(m).dx2;
      Real dx3 = pmbp->pmb->mb_size.d_view(m).dx3;
      Real dvol = dx1*dx2*dx3;
      
      lvol += dvol;
      Real rho = u0(m,IDN,k,j,i);
      Real vx = u0(m,IM1,k,j,i)/rho;
      Real vy = u0(m,IM2,k,j,i)/rho;
      Real vz = u0(m,IM3,k,j,i)/rho;
      Real ke = 0.5*rho*(vx*vx + vy*vy + vz*vz);
      
      lEk += ke*dvol;
      lmrho += rho*dvol;
      ldrho2 += rho*rho*dvol;
    }, Kokkos::Sum<Real>(vol), Kokkos::Sum<Real>(Ek), 
       Kokkos::Sum<Real>(mean_rho), Kokkos::Sum<Real>(drho_rms));
  }

#ifdef MPI_PARALLEL
  Real mpi_buf[4], mpi_recv[4];
  mpi_buf[0] = vol;
  mpi_buf[1] = Ek;
  mpi_buf[2] = mean_rho;
  mpi_buf[3] = drho_rms;
  MPI_Allreduce(mpi_buf, mpi_recv, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  vol = mpi_recv[0];
  Ek = mpi_recv[1];
  mean_rho = mpi_recv[2];
  drho_rms = mpi_recv[3];
#endif

  // Calculate RMS values
  mean_rho /= vol;
  drho_rms = sqrt(drho_rms/vol - mean_rho*mean_rho);
  Real Mach = sqrt(2.0*Ek/vol)/1.0; // Assuming cs = 1.0

  // Write history output
  if (pdata != nullptr) {
    pdata->hdata[0] = pm->time;
    pdata->hdata[1] = pm->dt;
    pdata->hdata[2] = Ek/vol;
    pdata->hdata[3] = Mach;
    pdata->hdata[4] = drho_rms/mean_rho;
    pdata->nhist = 5;
    
    if (global_variable::my_rank == 0) {
      pdata->label[0] = "time";
      pdata->label[1] = "dt";
      pdata->label[2] = "Ek/V";
      pdata->label[3] = "Mach";
      pdata->label[4] = "drho_rms/rho";
    }
  }
  
  return;
}