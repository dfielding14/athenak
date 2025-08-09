//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_mhd_amr_wave.cpp
//! \brief Problem generator for testing turbulence driving with MHD and dynamically 
//!        moving refinement regions (traveling wave) to test div(B) preservation
//!
//! PURPOSE: This test problem is designed to verify that div(B) remains zero when
//! refinement boundaries move through the domain. It combines:
//! - Turbulence driving (to create complex field structures)
//! - MHD with initial magnetic field
//! - A traveling wave refinement pattern that continuously creates and destroys
//!   fine/coarse boundaries
//!
//! The moving refinement pattern ensures that blocks are repeatedly refined and
//! derefined, testing whether the div(B) fix properly handles dynamically changing
//! mesh structure. This is a stringent test of the AMR+MHD implementation.
//!
//! Test created for issue #595: https://github.com/IAS-Astrophysics/athenak/issues/595

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/turb_driver.hpp"
#include "driver/driver.hpp"
#include "pgen.hpp"

// Global variables for traveling wave refinement
static Real wave_speed;
static Real wave_width;
static Real wave_amplitude;
static int wave_direction;  // 0=x, 1=y, 2=z

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Sets up turbulence driving test with MHD and traveling wave refinement

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Read problem parameters for initial conditions
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real pgas0 = pin->GetOrAddReal("problem", "pgas0", 1.0);
  Real vx0 = pin->GetOrAddReal("problem", "vx0", 0.0);
  Real vy0 = pin->GetOrAddReal("problem", "vy0", 0.0);
  Real vz0 = pin->GetOrAddReal("problem", "vz0", 0.0);
  
  // Magnetic field parameters
  Real b0 = pin->GetOrAddReal("problem", "b0", 0.1);      // Field strength
  Real beta = pin->GetOrAddReal("problem", "beta", 100.0); // Plasma beta
  int field_config = pin->GetOrAddInteger("problem", "field_config", 0); 
  // field_config: 0=uniform-x, 1=uniform-y, 2=uniform-z, 3=diagonal, 4=helical
  
  // Calculate field strength from plasma beta if specified
  if (beta > 0.0 && beta < 1.0e10) {
    b0 = std::sqrt(2.0 * pgas0 / beta);
  }
  
  // Traveling wave refinement parameters
  wave_speed = pin->GetOrAddReal("problem", "wave_speed", 0.1);
  wave_width = pin->GetOrAddReal("problem", "wave_width", 0.2);
  wave_amplitude = pin->GetOrAddReal("problem", "wave_amplitude", 0.0);
  wave_direction = pin->GetOrAddInteger("problem", "wave_direction", 0);

  // Initialize MHD variables
  auto &u0 = pmbp->pmhd->u0;
  auto &b0_ = pmbp->pmhd->b0;
  auto &w0 = pmbp->pmhd->w0;
  
  Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  
  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  
  // Initialize uniform flow with optional density perturbation
  par_for("turb_mhd_wave_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    
    // Add optional density wave for visualization
    Real rho = rho0;
    if (wave_amplitude > 0.0) {
      Real phase = 0.0;
      if (wave_direction == 0) phase = 2.0 * M_PI * x1v;
      else if (wave_direction == 1) phase = 2.0 * M_PI * x2v;
      else if (wave_direction == 2) phase = 2.0 * M_PI * x3v;
      rho = rho0 * (1.0 + wave_amplitude * sin(phase));
    }
    
    w0(m,IDN,k,j,i) = rho;
    w0(m,IVX,k,j,i) = vx0;
    w0(m,IVY,k,j,i) = vy0;
    w0(m,IVZ,k,j,i) = vz0;
    w0(m,IPR,k,j,i) = pgas0;
    
    // Set conserved variables (without magnetic energy yet)
    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = rho * vx0;
    u0(m,IM2,k,j,i) = rho * vy0;
    u0(m,IM3,k,j,i) = rho * vz0;
    Real ekin = 0.5 * rho * (vx0*vx0 + vy0*vy0 + vz0*vz0);
    u0(m,IEN,k,j,i) = pgas0/gm1 + ekin;
  });
  
  // Set magnetic field configuration
  // field_config: 0=uniform-x, 1=uniform-y, 2=uniform-z, 3=diagonal, 4=helical
  if (field_config == 0) {
    // Uniform field in x-direction
    par_for("set_bx", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x1f(m,k,j,i) = b0;
    });
    par_for("set_by", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x2f(m,k,j,i) = 0.0;
    });
    par_for("set_bz", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x3f(m,k,j,i) = 0.0;
    });
    
  } else if (field_config == 1) {
    // Uniform field in y-direction
    par_for("set_bx", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x1f(m,k,j,i) = 0.0;
    });
    par_for("set_by", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x2f(m,k,j,i) = b0;
    });
    par_for("set_bz", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x3f(m,k,j,i) = 0.0;
    });
    
  } else if (field_config == 2) {
    // Uniform field in z-direction
    par_for("set_bx", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x1f(m,k,j,i) = 0.0;
    });
    par_for("set_by", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x2f(m,k,j,i) = 0.0;
    });
    par_for("set_bz", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x3f(m,k,j,i) = b0;
    });
    
  } else if (field_config == 3) {
    // Diagonal field (tests all components)
    Real bcomp = b0 / std::sqrt(3.0);
    par_for("set_bx", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x1f(m,k,j,i) = bcomp;
    });
    par_for("set_by", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x2f(m,k,j,i) = bcomp;
    });
    par_for("set_bz", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x3f(m,k,j,i) = bcomp;
    });
    
  } else if (field_config == 4) {
    // Helical field (more complex test)
    par_for("set_bx_helical", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1f = CellCenterX(i-is, nx1+1, x1min, x1max); // face position
      
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      
      b0_.x1f(m,k,j,i) = b0 * cos(2.0 * M_PI * x2v);
    });
    par_for("set_by_helical", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real x2f = CellCenterX(j-js, nx2+1, x2min, x2max); // face position
      
      b0_.x2f(m,k,j,i) = b0 * sin(2.0 * M_PI * x1v);
    });
    par_for("set_bz_helical", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0_.x3f(m,k,j,i) = b0 * 0.1; // Small z-component
    });
  }
  
  // Add magnetic energy to total energy
  par_for("add_emag", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Calculate cell-centered B from face-centered values
    Real bx = 0.5 * (b0_.x1f(m,k,j,i) + b0_.x1f(m,k,j,i+1));
    Real by = 0.5 * (b0_.x2f(m,k,j,i) + b0_.x2f(m,k,j+1,i));
    Real bz = 0.5 * (b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j,i));
    Real emag = 0.5 * (bx*bx + by*by + bz*bz);
    u0(m,IEN,k,j,i) += emag;
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn int RefinementCondition()
//! \brief Refinement condition for traveling wave test
//!
//! Creates a moving band of refinement that travels through the domain.
//! This continuously creates and destroys fine/coarse boundaries, testing
//! whether div(B) is properly maintained during dynamic refinement.

int RefinementCondition(MeshBlockPack* pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &nmb = pmbp->nmb_thispack;
  
  // Get current time from mesh
  Real time = pmbp->pmesh->time;
  
  // Calculate current position of refinement wave
  Real wave_pos = wave_speed * time;
  
  // Get domain bounds (assumes symmetric domain [-L,L])
  Real domain_min = -0.5;
  Real domain_max = 0.5;
  Real domain_length = domain_max - domain_min;
  
  // Wrap wave position periodically
  wave_pos = std::fmod(wave_pos, domain_length);
  if (wave_pos < 0.0) wave_pos += domain_length;
  wave_pos += domain_min;
  
  auto ref_flag = Kokkos::View<int*>("ref_flag", nmb);
  Real wave_min = wave_pos - 0.5 * wave_width;
  Real wave_max = wave_pos + 0.5 * wave_width;
  
  par_for("check_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    
    // Check if MeshBlock overlaps with traveling wave refinement region
    bool refine = false;
    
    if (wave_direction == 0) {
      // Wave travels in x-direction
      if (wave_min <= domain_max && wave_max >= domain_min) {
        // Normal case: wave is within domain
        if (x1max > wave_min && x1min < wave_max) {
          refine = true;
        }
      } else {
        // Handle periodic wrapping
        Real wrapped_min = wave_min;
        Real wrapped_max = wave_max;
        if (wrapped_min < domain_min) {
          wrapped_min += domain_length;
          if (x1min < wave_max || x1max > wrapped_min) {
            refine = true;
          }
        } else if (wrapped_max > domain_max) {
          wrapped_max -= domain_length;
          if (x1max > wave_min || x1min < wrapped_max) {
            refine = true;
          }
        }
      }
    } else if (wave_direction == 1) {
      // Wave travels in y-direction
      if (wave_min <= domain_max && wave_max >= domain_min) {
        if (x2max > wave_min && x2min < wave_max) {
          refine = true;
        }
      } else {
        // Handle periodic wrapping
        Real wrapped_min = wave_min;
        Real wrapped_max = wave_max;
        if (wrapped_min < domain_min) {
          wrapped_min += domain_length;
          if (x2min < wave_max || x2max > wrapped_min) {
            refine = true;
          }
        } else if (wrapped_max > domain_max) {
          wrapped_max -= domain_length;
          if (x2max > wave_min || x2min < wrapped_max) {
            refine = true;
          }
        }
      }
    } else if (wave_direction == 2) {
      // Wave travels in z-direction
      if (wave_min <= domain_max && wave_max >= domain_min) {
        if (x3max > wave_min && x3min < wave_max) {
          refine = true;
        }
      } else {
        // Handle periodic wrapping
        Real wrapped_min = wave_min;
        Real wrapped_max = wave_max;
        if (wrapped_min < domain_min) {
          wrapped_min += domain_length;
          if (x3min < wave_max || x3max > wrapped_min) {
            refine = true;
          }
        } else if (wrapped_max > domain_max) {
          wrapped_max -= domain_length;
          if (x3max > wave_min || x3min < wrapped_max) {
            refine = true;
          }
        }
      }
    }
    
    ref_flag(m) = refine ? 1 : 0;
  });
  
  // Copy results back to host
  auto h_ref_flag = Kokkos::create_mirror_view(ref_flag);
  Kokkos::deep_copy(h_ref_flag, ref_flag);
  
  // Find maximum refinement flag
  int max_flag = 0;
  for (int m = 0; m < nmb; ++m) {
    max_flag = std::max(max_flag, h_ref_flag(m));
  }
  
  return max_flag;
}

