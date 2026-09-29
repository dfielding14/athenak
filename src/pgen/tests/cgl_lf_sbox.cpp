//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_lf_sbox.cpp
//! \brief Deterministic CGL Landau-fluid shearing-box regression initial conditions.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "shearing_box/shearing_box.hpp"

void ProblemGenerator::CGLLFSbox(ParameterInput *pin, const bool restart) {
  auto *pmbp = pmy_mesh_->pmb_pack;
  auto *pmhd = pmbp->pmhd;
  if (pmhd == nullptr || !pmhd->peos->eos_data.is_cgl || pmhd->pcgl_lf == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "cgl_lf_sbox requires CGL MHD with Landau-fluid heat flux."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!pmy_mesh_->three_d || pmhd->psbox_u == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "cgl_lf_sbox requires a three-dimensional shearing box."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (restart) {
    return;
  }

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real ppar0 = pin->GetOrAddReal("problem", "ppar0", 1.0);
  const Real pperp0 = pin->GetOrAddReal("problem", "pperp0", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.02);
  const Real bx0 = pin->GetOrAddReal("problem", "bx0", 0.3);
  const Real by0 = pin->GetOrAddReal("problem", "by0", 0.2);
  const Real bz0 = pin->GetOrAddReal("problem", "bz0", 0.1);
  const int nwx = pin->GetOrAddInteger("problem", "nwx", 1);
  const int nwy = pin->GetOrAddInteger("problem", "nwy", 1);
  const int nwz = pin->GetOrAddInteger("problem", "nwz", 1);
  if (rho0 <= 0.0 || ppar0 <= 0.0 || pperp0 <= 0.0 ||
      amp < 0.0 || amp >= 1.0 || (SQR(bx0) + SQR(by0) + SQR(bz0)) <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "cgl_lf_sbox received invalid initial-state parameters."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const Real xmin = pmy_mesh_->mesh_size.x1min;
  const Real ymin = pmy_mesh_->mesh_size.x2min;
  const Real zmin = pmy_mesh_->mesh_size.x3min;
  const Real kx = 2.0*M_PI*static_cast<Real>(nwx)/
                  (pmy_mesh_->mesh_size.x1max - xmin);
  const Real ky = 2.0*M_PI*static_cast<Real>(nwy)/
                  (pmy_mesh_->mesh_size.x2max - ymin);
  const Real kz = 2.0*M_PI*static_cast<Real>(nwz)/
                  (pmy_mesh_->mesh_size.x3max - zmin);

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto w0 = pmhd->w0;
  auto bcc0 = pmhd->bcc0;
  auto b0 = pmhd->b0;

  par_for("cgl_lf_sbox_prim", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                              size.d_view(m).x1max);
    const Real y = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                              size.d_view(m).x2max);
    const Real z = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
    const Real phase = kx*(x - xmin) + ky*(y - ymin) + kz*(z - zmin);
    w0(m,IDN,k,j,i) = rho0;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    w0(m,IPR,k,j,i) = ppar0*(1.0 + amp*sin(phase));
    w0(m,IPP,k,j,i) = pperp0*(1.0 + 0.5*amp*cos(phase));
    bcc0(m,IBX,k,j,i) = bx0;
    bcc0(m,IBY,k,j,i) = by0;
    bcc0(m,IBZ,k,j,i) = bz0;
  });

  par_for("cgl_lf_sbox_b1", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    b0.x1f(m,k,j,i) = bx0;
  });
  par_for("cgl_lf_sbox_b2", DevExeSpace(), 0, nmb - 1, ks, ke, js, je + 1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    b0.x2f(m,k,j,i) = by0;
  });
  par_for("cgl_lf_sbox_b3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    b0.x3f(m,k,j,i) = bz0;
  });

  pmhd->peos->PrimToCons(w0, bcc0, pmhd->u0, is, ie, js, je, ks, ke);
}
