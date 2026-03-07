//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_mixing_test.cpp
//  \brief Passive-scalar advection-diffusion test with a top-hat blob in uniform flow.
//
//  This user pgen initializes a single passive scalar as a circular/spherical blob
//  embedded in a uniform-density medium with constant diagonal velocity.  The reference
//  verification path uses scalar_only=true so the velocity field is frozen and the scalar
//  evolves under advection plus optional scalar diffusion with a known periodic solution.

#include <algorithm>
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

Real ActiveBoxMinLength(Mesh *pm) {
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  Real min_len = lx;
  if (pm->multi_d) min_len = std::min(min_len, ly);
  if (pm->three_d) min_len = std::min(min_len, lz);
  return min_len;
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem generator for passive-scalar advection-diffusion verification.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "scalar_mixing_test requires a <hydro> block." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (pmy_mesh_->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "scalar_mixing_test only supports 2D or 3D problems." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int nscalars = pmbp->phydro->nscalars;
  if (nscalars < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "scalar_mixing_test requires nscalars >= 1." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  auto &mesh_size = pmy_mesh_->mesh_size;
  auto &u0 = pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;

  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  int nmb = pmbp->nmb_thispack;
  int nhydro = pmbp->phydro->nhydro;

  Real x1min = mesh_size.x1min;
  Real x1max = mesh_size.x1max;
  Real x2min = mesh_size.x2min;
  Real x2max = mesh_size.x2max;
  Real x3min = mesh_size.x3min;
  Real x3max = mesh_size.x3max;

  Real lx = x1max - x1min;
  Real ly = x2max - x2min;
  Real lz = x3max - x3min;

  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real temp0 = pin->GetOrAddReal("problem", "temp0", 1.0);
  Real scalar_inside = pin->GetOrAddReal("problem", "scalar_inside", 1.0);
  Real scalar_outside = pin->GetOrAddReal("problem", "scalar_outside", 0.0);
  Real crossing_time = pin->GetOrAddReal("problem", "crossing_time", 1.0);

  if (crossing_time <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "problem/crossing_time must be positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real blob_radius = pin->GetOrAddReal("problem", "blob_radius",
                                       0.125 * ActiveBoxMinLength(pmy_mesh_));
  if (blob_radius <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "problem/blob_radius must be positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real center_x1 = pin->GetOrAddReal("problem", "blob_center_x1", x1min + 0.25*lx);
  Real center_x2 = pin->GetOrAddReal("problem", "blob_center_x2", x2min + 0.25*ly);
  Real center_x3 = pin->GetOrAddReal("problem", "blob_center_x3", x3min + 0.25*lz);

  Real vx0 = lx / crossing_time;
  Real vy0 = pmy_mesh_->multi_d ? ly / crossing_time : 0.0;
  Real vz0 = pmy_mesh_->three_d ? lz / crossing_time : 0.0;
  Real gm1 = eos.is_ideal ? (eos.gamma - 1.0) : 0.0;
  bool three_d = pmy_mesh_->three_d;
  int active_dims = three_d ? 3 : 2;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar mixing advection-diffusion test:" << std::endl
              << "  active_dims = " << active_dims << std::endl
              << "  rho0 = " << rho0 << ", temp0 = " << temp0 << std::endl
              << "  scalar_inside = " << scalar_inside
              << ", scalar_outside = " << scalar_outside << std::endl
              << "  blob_radius = " << blob_radius << std::endl
              << "  blob_center = (" << center_x1 << ", " << center_x2 << ", "
              << (three_d ? center_x3 : 0.0) << ")" << std::endl
              << "  crossing_time = " << crossing_time << std::endl
              << "  velocity = (" << vx0 << ", " << vy0 << ", " << vz0 << ")" << std::endl;
  }

  par_for("scalar_mixing_test_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &mb_x1min = size.d_view(m).x1min;
    Real &mb_x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, mb_x1min, mb_x1max);

    Real &mb_x2min = size.d_view(m).x2min;
    Real &mb_x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, nx2, mb_x2min, mb_x2max);

    Real &mb_x3min = size.d_view(m).x3min;
    Real &mb_x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, nx3, mb_x3min, mb_x3max);

    Real dx1 = x1v - center_x1;
    if (dx1 > 0.5*lx) dx1 -= lx;
    if (dx1 < -0.5*lx) dx1 += lx;

    Real dx2 = x2v - center_x2;
    if (dx2 > 0.5*ly) dx2 -= ly;
    if (dx2 < -0.5*ly) dx2 += ly;

    Real dx3 = 0.0;
    if (three_d) {
      dx3 = x3v - center_x3;
      if (dx3 > 0.5*lz) dx3 -= lz;
      if (dx3 < -0.5*lz) dx3 += lz;
    }

    Real r2 = dx1*dx1 + dx2*dx2 + dx3*dx3;
    Real scalar0 = (r2 <= blob_radius*blob_radius) ? scalar_inside : scalar_outside;

    u0(m,IDN,k,j,i) = rho0;
    u0(m,IM1,k,j,i) = rho0 * vx0;
    u0(m,IM2,k,j,i) = rho0 * vy0;
    u0(m,IM3,k,j,i) = three_d ? rho0 * vz0 : 0.0;
    if (eos.is_ideal) {
      Real ekin = 0.5*rho0*(vx0*vx0 + vy0*vy0 + vz0*vz0);
      u0(m,IEN,k,j,i) = rho0 * temp0 / gm1 + ekin;
    }

    u0(m,nhydro,k,j,i) = rho0 * scalar0;
    for (int n = 1; n < nscalars; ++n) {
      u0(m,nhydro+n,k,j,i) = rho0 * scalar_outside;
    }
  });
}
