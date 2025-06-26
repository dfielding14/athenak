//========================================================================================
// AthenaK turbulence driving with moving AMR refinement test
// Copyright(C) 2025 AthenaK code team
// Licensed under the 3-clause BSD License
//========================================================================================
//! \file turb_amr_wave_test.cpp
//! \brief Problem generator for testing turbulence driving with dynamically moving
//!        refinement regions (traveling wave) to diagnose derefinement issues

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
#include "driver/driver.hpp"
#include "pgen/pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator()
//! \brief Sets up turbulence driving test with traveling wave refinement

void ProblemGenerator(ParameterInput *pin, const bool restart) {
  if (restart) return;

  // Read problem parameters
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real pgas0 = pin->GetOrAddReal("problem", "pgas0", 1.0);
  Real vx0 = pin->GetOrAddReal("problem", "vx0", 0.0);
  Real vy0 = pin->GetOrAddReal("problem", "vy0", 0.0);
  Real vz0 = pin->GetOrAddReal("problem", "vz0", 0.0);

  // Read wave parameters for moving refinement
  Real wave_speed = pin->GetOrAddReal("problem", "wave_speed", 0.5);
  Real wave_width = pin->GetOrAddReal("problem", "wave_width", 0.3);
  
  // Get mesh and coordinates
  auto &mesh = pmy_mesh_;
  auto &indcs = mesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;

  // Capture variables for kernel
  Real gam = peos->eos_data.gamma;
  
  // Initialize uniform flow
  auto &w0 = pmb->w0;
  par_for("turb_amr_wave_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m,IDN,k,j,i) = rho0;
    w0(m,IVX,k,j,i) = vx0;
    w0(m,IVY,k,j,i) = vy0;
    w0(m,IVZ,k,j,i) = vz0;
    
    // Set pressure/internal energy
    Real eint = pgas0/(gam - 1.0);
    w0(m,IEN,k,j,i) = eint/rho0;
  });

  // Convert to conserved variables
  if (phydro != nullptr) {
    phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void UserWorkAfterLoop()
//! \brief Save turbulence diagnostics for analysis

void UserWorkAfterLoop(Mesh *pmesh, ParameterInput *pin, SimTime &tm) {
  // This function can be used to output diagnostics about the turbulence
  // driver array sizes and force field continuity
  if (Globals::my_rank == 0) {
    std::cout << "# Turbulence AMR wave test completed" << std::endl;
    std::cout << "# Final time: " << tm.time << std::endl;
    std::cout << "# Final cycle: " << tm.ncycle << std::endl;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real RefinementCondition()
//! \brief Traveling wave refinement condition for testing AMR derefinement

Real RefinementCondition(MeshBlockPack *pmbp) {
  auto &pmb = pmbp->pmb_pack;
  auto &indcs = pmbp->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;

  // Get current time
  Real time = pmbp->pmesh->time;
  
  // Wave parameters - should match those in problem generator
  Real wave_speed = 0.5;  // Wave propagation speed
  Real wave_width = 0.3;  // Width of refined region
  Real wave_period = 4.0; // Period for wave to cross domain
  
  // Calculate current wave center position (moves in x1 direction)
  Real x1_center = -1.0 + fmod(2.0 + wave_speed * time, 2.0);
  
  // Check if any cell in this MeshBlock is within the wave
  Real threshold = 0.0;
  
  par_reduce(loop_pattern_mdrange_tag, "wave_refine", DevExeSpace(),
    0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i, Real &result) {
      // Get cell center position
      Real &x1min = pmb(m)->mb_size.x1min;
      Real &x1max = pmb(m)->mb_size.x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      
      // Calculate distance from wave center
      Real dist = fabs(x1v - x1_center);
      
      // Handle periodic boundary
      if (dist > 1.0) dist = 2.0 - dist;
      
      // Set refinement level based on distance
      Real refine_level = 0.0;
      if (dist < wave_width) {
        refine_level = 1.0;  // Refine to next level
      }
      
      // Take maximum over all cells
      result = fmax(result, refine_level);
    },
    Kokkos::Max<Real>(threshold));

  return threshold;
}