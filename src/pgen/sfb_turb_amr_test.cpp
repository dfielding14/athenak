//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sfb_turb_amr_test.cpp
//! \brief Test problem generator for SFB turbulence driving with AMR
//! Sets up a uniform medium with static refinement to test SFB turbulence continuity

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Sets up SFB turbulence test with AMR

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Read problem parameters
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real pgas0 = pin->GetOrAddReal("problem", "pgas0", 1.0);
  Real vx0 = pin->GetOrAddReal("problem", "vx0", 0.0);
  Real vy0 = pin->GetOrAddReal("problem", "vy0", 0.0);
  Real vz0 = pin->GetOrAddReal("problem", "vz0", 0.0);
  
  // Initialize Hydro variables
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  
  Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  
  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  
  par_for("pgen_sfb_turb_amr", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IDN,k,j,i) = rho0;
    u0(m,IM1,k,j,i) = rho0*vx0;
    u0(m,IM2,k,j,i) = rho0*vy0;
    u0(m,IM3,k,j,i) = rho0*vz0;
    Real ekin = 0.5*rho0*(vx0*vx0 + vy0*vy0 + vz0*vz0);
    u0(m,IEN,k,j,i) = pgas0/gm1 + ekin;
    
    // Also set primitives
    w0(m,IDN,k,j,i) = rho0;
    w0(m,IVX,k,j,i) = vx0;
    w0(m,IVY,k,j,i) = vy0;
    w0(m,IVZ,k,j,i) = vz0;
    w0(m,IPR,k,j,i) = pgas0;
  });
  
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int RefinementCondition()
//  \brief refinement condition for SFB turbulence AMR test - using static refinement only

int RefinementCondition(MeshBlockPack* pmbp) {
  // Using static refinement from input file, no dynamic refinement
  return 0;
}