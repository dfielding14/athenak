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
//  \brief Sets up turbulence test with AMR

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
  
  par_for("pgen_turb_amr", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
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
//  \brief refinement condition for turbulence AMR test

int RefinementCondition(MeshBlockPack* pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is, &ie = indcs.ie;
  int &js = indcs.js, &je = indcs.je;
  int &ks = indcs.ks, &ke = indcs.ke;
  int &nmb = pmbp->nmb_thispack;
  
  // Define refinement region parameters
  Real x_ref_min = -0.25;
  Real x_ref_max = 0.25;
  Real y_ref_min = -0.25;
  Real y_ref_max = 0.25;
  Real z_ref_min = -0.25;
  Real z_ref_max = 0.25;
  
  auto ref_flag = Kokkos::View<int*>("ref_flag", nmb);
  
  par_for("check_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    
    // Check if MeshBlock overlaps with refinement region
    bool refine = false;
    if (x1max > x_ref_min && x1min < x_ref_max &&
        x2max > y_ref_min && x2min < y_ref_max &&
        x3max > z_ref_min && x3min < z_ref_max) {
      refine = true;
    }
    
    if (refine) {
      ref_flag(m) = 1;  // flag for refinement
    } else {
      ref_flag(m) = 0;  // no refinement
    }
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