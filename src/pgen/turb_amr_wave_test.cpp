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

#include <Kokkos_Core.hpp>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "mesh/mesh_refinement.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/turb_driver.hpp"
#include "driver/driver.hpp"
#include "pgen.hpp"

// Static variables to store wave parameters (accessible to RefinementCondition)
static Real wave_speed_g = 0.5;
static Real wave_width_g = 0.15;
static Real domain_size_g = 2.0;

void RefinementCondition(MeshBlockPack* pmbp);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Sets up turbulence driving test with traveling wave refinement

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {

  user_ref_func  = RefinementCondition;

  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Read problem parameters
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real pgas0 = pin->GetOrAddReal("problem", "pgas0", 1.0);
  Real vx0 = pin->GetOrAddReal("problem", "vx0", 0.0);
  Real vy0 = pin->GetOrAddReal("problem", "vy0", 0.0);
  Real vz0 = pin->GetOrAddReal("problem", "vz0", 0.0);
  
  // Read wave parameters and store in static variables
  wave_speed_g = pin->GetOrAddReal("problem", "wave_speed", 0.5);
  wave_width_g = pin->GetOrAddReal("problem", "wave_width", 0.15);
  
  // Get domain size from mesh
  domain_size_g = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;

  // Initialize Hydro variables
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  
  Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  
  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  
  // Initialize uniform flow
  par_for("turb_amr_wave_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m,IDN,k,j,i) = rho0;
    w0(m,IVX,k,j,i) = vx0;
    w0(m,IVY,k,j,i) = vy0;
    w0(m,IVZ,k,j,i) = vz0;
    
    // Set pressure/internal energy
    w0(m,IEN,k,j,i) = pgas0/rho0/gm1;
  });

  // Convert to conserved variables
  pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RefinementCondition()
//! \brief Refinement condition based on traveling Gaussian wave
//!        Creates a narrow refined region that moves through the domain

void RefinementCondition(MeshBlockPack* pmbp) {
  // Get mesh parameters
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &nmb = pmbp->nmb_thispack;
  int mbs = pmbp->gids;  // starting global ID for this pack
  
  // Get current time from mesh
  Real time = pmbp->pmesh->time;
  
  // Print diagnostic info
  static int call_count = 0;
  call_count++;
  
  // Track turbulence force statistics
  Real force_mag_refined = 0.0;
  Real force_mag_unrefined = 0.0;
  int ncells_refined = 0;
  int ncells_unrefined = 0;
  
  // Check if turbulence driver exists
  bool has_turb = (pmbp->pturb != nullptr);
  
  // Get max refinement level early for debug output
  int max_level = pmbp->pmesh->max_level;
  
  if (global_variable::my_rank == 0 && call_count % 5 == 1) {
    std::cout << "### RefinementCondition called at cycle " << pmbp->pmesh->ncycle
              << ", time=" << time 
              << ", nmb=" << nmb << ", mbs=" << mbs 
              << ", has_turb=" << has_turb 
              << ", max_level=" << max_level << std::endl;
  }
  
  // Use wave parameters from input file
  Real wave_speed = wave_speed_g;
  Real L = domain_size_g;
  Real wave_width = wave_width_g;
  
  // Get the refinement flag array from MeshRefinement
  auto &refine_flag_ = pmbp->pmesh->pmr->refine_flag;
  
  // Capture variables for lambda
  int call_count_local = call_count;
  
  // Check each MeshBlock
  par_for("check_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    // Get MeshBlock bounds
    Real x1min = size.d_view(m).x1min;
    Real x1max = size.d_view(m).x1max;
    Real x1_center = 0.5 * (x1min + x1max);
    
    // Calculate wave center position (with periodic wrapping)
    // Position wave at x = -0.625 to align with MeshBlock center for 8x8x8 grid
    Real wave_center = -0.625 + fmod(wave_speed * time, L);
    if (wave_center > 1.0) wave_center -= L;
    
    // Calculate distance from wave center with periodic boundary consideration
    Real distance = fabs(x1_center - wave_center);
    // Handle periodic wrap-around
    if (distance > L/2.0) {
      distance = L - distance;
    }
    
    // Gaussian profile for refinement
    Real gaussian_value = exp(-distance*distance/(2.0*wave_width*wave_width));
    
    // Get current refinement level
    int current_level = pmbp->pmesh->lloc_eachmb[m+mbs].level;
    
    // Improved refinement logic with smoother transitions and hysteresis
    int flag = 0;
    
    // Define thresholds with hysteresis to prevent oscillations
    // ULTRA-aggressive thresholds to force level 3 refinement
    Real refine_threshold = 0.5;   // Lower threshold to refine more easily
    Real derefine_threshold = 0.001; // Extremely low threshold for derefinement
    
    // Force refinement at wave center by using very small thresholds
    if (current_level == 3) {
      refine_threshold = 0.5;    // Refine if gaussian > 0.5
    } else if (current_level == 4) {
      refine_threshold = 0.005;  // Refine if gaussian > 0.005 (adjusted for actual values)
    } else if (current_level == 5) {
      refine_threshold = 0.0001; // Refine if gaussian > 0.0001
    } else if (current_level == 6) {
      refine_threshold = 0.000001; // Refine if gaussian > 0.000001
    } else {
      refine_threshold = 0.0000001; // Force refinement for any detectable signal
    }
    
    // Apply refinement logic
    if (gaussian_value > refine_threshold && current_level < max_level) {
      flag = 1;   // refine in the wave peak
    } else if (gaussian_value < derefine_threshold && current_level > 0) {
      flag = -1;  // derefine far from wave
    } else {
      flag = 0;   // no change in transition region (hysteresis band)
    }
    
    // Debug output for first few blocks - show more frequently to catch refinement
    if (m < 8 && call_count_local <= 50 && call_count_local % 2 == 1) {
      printf("MB %d: level=%d/%d, x1_center=%f, wave_center=%f, distance=%f, gauss=%f, flag=%d\n", 
             m, current_level, max_level, x1_center, wave_center, distance, gaussian_value, flag);
    }
    
    // Set refinement flag
    refine_flag_.d_view(m+mbs) = flag;
  });
  
  // Calculate turbulence force statistics if turbulence driver is active
  if (has_turb && pmbp->pturb->force.extent(0) > 0) {
    auto force = pmbp->pturb->force;
    int is = indcs.is, ie = indcs.ie;
    int js = indcs.js, je = indcs.je;
    int ks = indcs.ks, ke = indcs.ke;
    
    // Calculate force magnitudes in refined vs unrefined regions
    const int nkji = indcs.nx3 * indcs.nx2 * indcs.nx1;
    const int nmkji = nmb * nkji;
    const int nji = indcs.nx2 * indcs.nx1;
    
    Kokkos::parallel_reduce("turb_force_stats", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &lforce_ref, Real &lforce_unref, 
                  int &lcells_ref, int &lcells_unref) {
      int m = idx / nkji;
      int k = (idx - m*nkji) / nji;
      int j = (idx - m*nkji - k*nji) / indcs.nx1;
      int i = (idx - m*nkji - k*nji - j*indcs.nx1) + is;
      k += ks;
      j += js;
      // Get cell center position
      Real x1 = CellCenterX(i, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      
      // Calculate wave position and distance
      Real wave_center_local = -0.625 + fmod(wave_speed * time, L);
      if (wave_center_local > 1.0) wave_center_local -= L;
      Real distance_local = fabs(x1 - wave_center_local);
      if (distance_local > L/2.0) distance_local = L - distance_local;
      
      // Calculate Gaussian value
      Real gaussian_value_local = exp(-distance_local*distance_local/(2.0*wave_width*wave_width));
      
      // Calculate force magnitude
      Real fmag_local = 0.0;
      if (m < force.extent(0)) {  // Check bounds
        fmag_local = sqrt(SQR(force(m,0,k,j,i)) + SQR(force(m,1,k,j,i)) + SQR(force(m,2,k,j,i)));
      }
      
      // Accumulate statistics based on location
      if (gaussian_value_local > 0.5) {  // In refined region
        lforce_ref += fmag_local;
        lcells_ref += 1;
      } else {  // In unrefined region
        lforce_unref += fmag_local;
        lcells_unref += 1;
      }
    }, Kokkos::Sum<Real>(force_mag_refined), Kokkos::Sum<Real>(force_mag_unrefined),
       Kokkos::Sum<int>(ncells_refined), Kokkos::Sum<int>(ncells_unrefined));
    
    // Print force statistics periodically
    if (global_variable::my_rank == 0 && call_count % 10 == 1) {
      Real avg_force_refined = (ncells_refined > 0) ? force_mag_refined/ncells_refined : 0.0;
      Real avg_force_unrefined = (ncells_unrefined > 0) ? force_mag_unrefined/ncells_unrefined : 0.0;
      
      std::cout << "### Turbulence Force Statistics at time=" << time << std::endl;
      std::cout << "    Refined region: avg_force=" << avg_force_refined 
                << " (ncells=" << ncells_refined << ")" << std::endl;
      std::cout << "    Unrefined region: avg_force=" << avg_force_unrefined 
                << " (ncells=" << ncells_unrefined << ")" << std::endl;
      std::cout << "    Force ratio (refined/unrefined): " 
                << ((avg_force_unrefined > 0) ? avg_force_refined/avg_force_unrefined : 0.0) 
                << std::endl;
    }
  }
  
  return;
}

