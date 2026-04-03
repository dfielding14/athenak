//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_sphere_advection.cpp
//  \brief Problem generator for passive scalar sphere advection with uniform flow.
//
//  Initializes a sphere of scalar concentration and advects it with a constant velocity
//  field pointing along (1,1,1). Designed for scalar_only runs with diffusion enabled.

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen.hpp"
#include "globals.hpp"

namespace {
void NoOpSource(Mesh *pm, const Real bdt) {
  (void)pm;
  (void)bdt;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar sphere advection.
//
//  Initializes density, uniform velocity, and passive scalar(s) in a sphere.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Scalar sphere advection requires Hydro, but no <hydro> block "
              << "in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  int nscalars = pmbp->phydro->nscalars;
  if (nscalars == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Scalar sphere advection requires nscalars > 0" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Enroll a no-op source term if user_srcs is enabled in the input.
  if (pin->GetOrAddBoolean("problem", "user_srcs", false)) {
    user_srcs_func = NoOpSource;
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  int nhydro = pmbp->phydro->nhydro;
  int nmb = pmbp->nmb_thispack;

  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  const bool is_ideal = eos.is_ideal;
  Real gm1 = eos.gamma - 1.0;

  // Domain extents for defaults
  Real x1min = pmy_mesh_->mesh_size.x1min;
  Real x1max = pmy_mesh_->mesh_size.x1max;
  Real x2min = pmy_mesh_->mesh_size.x2min;
  Real x2max = pmy_mesh_->mesh_size.x2max;
  Real x3min = pmy_mesh_->mesh_size.x3min;
  Real x3max = pmy_mesh_->mesh_size.x3max;
  Real lx = x1max - x1min;
  Real ly = x2max - x2min;
  Real lz = x3max - x3min;
  Real min_len = lx;
  if (pmy_mesh_->multi_d && ly < min_len) min_len = ly;
  if (pmy_mesh_->three_d && lz < min_len) min_len = lz;

  // Input parameters
  Real den = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real temp = pin->GetOrAddReal("problem", "temp0", 1.0);
  Real sphere_radius = pin->GetOrAddReal("problem", "sphere_radius", 0.1*min_len);
  Real sphere_x1 = pin->GetOrAddReal("problem", "sphere_x1", 0.5*(x1min + x1max));
  Real sphere_x2 = pin->GetOrAddReal("problem", "sphere_x2", 0.5*(x2min + x2max));
  Real sphere_x3 = pin->GetOrAddReal("problem", "sphere_x3", 0.5*(x3min + x3max));
  Real scalar_in = 1.0;
  if (pin->DoesParameterExist("problem", "scalar_in")) {
    scalar_in = pin->GetReal("problem", "scalar_in");
  } else {
    scalar_in = pin->GetOrAddReal("problem", "scalar_init", 1.0);
  }
  Real scalar_out = pin->GetOrAddReal("problem", "scalar_out", 0.0);
  Real flow_speed = pin->GetOrAddReal("problem", "flow_speed", 1.0);

  Real vcomp = flow_speed / std::sqrt(3.0);
  Real mom1 = den * vcomp;
  Real mom2 = den * vcomp;
  Real mom3 = den * vcomp;
  Real ke = 0.5 * den * (vcomp*vcomp + vcomp*vcomp + vcomp*vcomp);
  Real eint = (is_ideal ? den * temp / gm1 : 0.0);
  Real radius2 = sphere_radius * sphere_radius;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar sphere advection problem:" << std::endl
              << "  sphere_center = (" << sphere_x1 << ", " << sphere_x2
              << ", " << sphere_x3 << ")" << std::endl
              << "  sphere_radius = " << sphere_radius << std::endl
              << "  scalar_in/out = " << scalar_in << " / " << scalar_out << std::endl
              << "  flow_speed = " << flow_speed
              << " (v = " << vcomp << " along (1,1,1))" << std::endl
              << "  nscalars = " << nscalars << std::endl;
  }

  par_for("pgen_scalar_sphere", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min_m = size.d_view(m).x1min;
    Real &x1max_m = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, x1min_m, x1max_m);

    Real &x2min_m = size.d_view(m).x2min;
    Real &x2max_m = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, nx2, x2min_m, x2max_m);

    Real &x3min_m = size.d_view(m).x3min;
    Real &x3max_m = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, nx3, x3min_m, x3max_m);

    Real dx = x1v - sphere_x1;
    Real dy = x2v - sphere_x2;
    Real dz = x3v - sphere_x3;
    Real r2 = dx*dx + dy*dy + dz*dz;
    Real scalar = (r2 <= radius2) ? scalar_in : scalar_out;

    u0(m,IDN,k,j,i) = den;
    u0(m,IM1,k,j,i) = mom1;
    u0(m,IM2,k,j,i) = mom2;
    u0(m,IM3,k,j,i) = mom3;
    if (is_ideal) {
      u0(m,IEN,k,j,i) = eint + ke;
    }
    for (int ns = 0; ns < nscalars; ++ns) {
      u0(m, nhydro+ns, k, j, i) = den * scalar;
    }
  });

  return;
}
