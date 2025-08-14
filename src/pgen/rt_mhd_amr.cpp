//! and uses AMR to refine the interface region. It can be used to verify that the
//! magnetic field remains divergence-free across refinement boundaries.
//!
//! The problem was created specifically to test the fix for issue #595:
//! https://github.com/IAS-Astrophysics/athenak/issues/595
//!
//! REFERENCE: Based on the RT instability test described in 
//! Stone et al. (2008), ApJS, 178, 137 with added magnetic field

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
#include "pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Sets up Rayleigh-Taylor instability with magnetic field

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  // Read problem parameters
  Real amp = pin->GetOrAddReal("problem","amp",0.01);
  Real drat = pin->GetOrAddReal("problem","drat",2.0);
  Real b0 = pin->GetOrAddReal("problem","b0",0.1);
  int ipert = pin->GetOrAddInteger("problem","ipert",1);
  
  // Initialize MHD variables
  auto &u0 = pmbp->pmhd->u0;
  auto &b0_ = pmbp->pmhd->b0;
  auto &w0 = pmbp->pmhd->w0;
  
  Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  Real grav = 0.1;
  
  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  
  par_for("pgen_rt_mhd", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
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
    
    // Set density and pressure based on height
    Real den, pres;
    if (x2v > 0.0) {
      den = drat;
      pres = 1.0 + grav*(drat*x2v);
    } else {
      den = 1.0;
      pres = 1.0 + grav*(drat*x2v);
    }
    
    // Add perturbation
    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    
    if (ipert == 1) {
      // Single mode perturbation
      vy = amp * sin(2.0*M_PI*x1v) * cos(M_PI*x2v);
    } else if (ipert == 2) {
      // Multi-mode perturbation
      vy = amp * (sin(2.0*M_PI*x1v) + 0.5*sin(4.0*M_PI*x1v)) * cos(M_PI*x2v);
      if (indcs.nx3 > 1) {
        vy += amp * 0.25 * sin(2.0*M_PI*x3v) * cos(M_PI*x2v);
      }
    }
    
    u0(m,IDN,k,j,i) = den;
    u0(m,IM1,k,j,i) = den*vx;
    u0(m,IM2,k,j,i) = den*vy;
    u0(m,IM3,k,j,i) = den*vz;
    
    Real ekin = 0.5*den*(vx*vx + vy*vy + vz*vz);
    u0(m,IEN,k,j,i) = pres/gm1 + ekin;
    
    // Also set primitives
    w0(m,IDN,k,j,i) = den;
    w0(m,IVX,k,j,i) = vx;
    w0(m,IVY,k,j,i) = vy;
    w0(m,IVZ,k,j,i) = vz;
    w0(m,IPR,k,j,i) = pres;
  });
  
  // Set magnetic field (uniform horizontal field)
  par_for("pgen_rt_b", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0_.x1f(m,k,j,i) = b0;
  });
  
  par_for("pgen_rt_b", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0_.x2f(m,k,j,i) = 0.0;
  });
  
  par_for("pgen_rt_b", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0_.x3f(m,k,j,i) = 0.0;
  });
  
  // Add magnetic energy to total energy
  par_for("pgen_rt_emag", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IEN,k,j,i) += 0.5*b0*b0;
  });
  
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int RefinementCondition()
//  \brief refinement condition for RT instability
//  Refine at the interface and regions with high density gradients

int RefinementCondition(MeshBlockPack* pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is, &ie = indcs.ie;
  int &js = indcs.js, &je = indcs.je;
  int &ks = indcs.ks, &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  
  auto &u0 = pmbp->pmhd->u0;
  
  auto ref_flag = Kokkos::View<int*>("ref_flag", nmb);
  
  par_for("check_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    
    // Refine near the interface (y=0)
    if (x2min < 0.1 && x2max > -0.1) {
      ref_flag(m) = 1;
    } else {
      ref_flag(m) = 0;
    }
    
    // Also check for large density gradients
    Real max_drho = 0.0;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je-1; j++) {
        for (int i=is; i<=ie; i++) {
          Real drho = fabs(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j,i));
          max_drho = fmax(max_drho, drho);
        }
      }
    }
    
    if (max_drho > 0.5) {
      ref_flag(m) = 1;
    }
  });
  
  // Copy device view to host
  auto ref_flag_h = Kokkos::create_mirror_view(ref_flag);
  Kokkos::deep_copy(ref_flag_h, ref_flag);
  
  // Check if any cell requires refinement
  int refine = 0;
  for (int m=0; m<nmb; m++) {
    if (ref_flag_h(m) > 0) refine = 1;
  }
  
  return refine;
}

