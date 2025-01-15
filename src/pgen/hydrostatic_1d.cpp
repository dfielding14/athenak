//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mass_removal_test.cpp
//  \brief Problem generator for testing mass removal
#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

void UserSource(Mesh* pm, const Real bdt);
void GravitySource(Mesh* pm, const Real bdt);
void Static(Mesh* pm);

Real scale_height, g_acc;

//----------------------------------------------------------------------------------------
//  \brief Problem Generator for mass removal

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Enroll user functions 
  user_srcs_func = UserSource;
  user_bcs_func = Static;

  if (restart) return;

  // Read in constants
  g_acc = pin->GetOrAddReal("problem", "grav_acc",1.0);
  scale_height = 1/g_acc;

  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;

  // Initialize Hydro variables -------------------------------
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  Real &H = scale_height;

  // Set initial conditions
  par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real z = CellCenterX(k-ks, nx3, x3min, x3max);
    Real rho = exp(-z/H);

    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = rho/gm1 +
       0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
       SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
  });

  return;
}

//===========================================================================//
//                              Source Terms                                 //
//===========================================================================//

void UserSource(Mesh* pm, const Real bdt) {
  GravitySource(pm, bdt);
}

void GravitySource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  Real &g = g_acc;

  par_for("gravity_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real src = bdt*g*w0(m,IDN,k,j,i);
	
    u0(m,IM3,k,j,i) -= src;
    u0(m,IEN,k,j,i) -= src*w0(m,IVZ,k,j,i);
  });

  return;
}

void Static(Mesh* pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;

  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  Real &H = scale_height;

  par_for("static_x3", DevExeSpace(),0,nmb1,0,(ng-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
      Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
      Real rho = exp(-x3v/H);
      u0(m,IDN,k,j,i) = rho;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = rho/gm1;
    }
    
    if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
      Real x3v = CellCenterX((ke+k+1)-ks, indcs.nx3, x3min, x3max);
      Real rho = exp(-x3v/H);
      u0(m,IDN,ke+k+1,j,i) = rho;
      u0(m,IM1,ke+k+1,j,i) = 0.0;
      u0(m,IM2,ke+k+1,j,i) = 0.0;
      u0(m,IM3,ke+k+1,j,i) = 0.0;
      u0(m,IEN,ke+k+1,j,i) = rho/gm1;
    }

  });
}
