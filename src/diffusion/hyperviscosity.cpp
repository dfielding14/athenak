//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hyperviscosity.cpp
//! \brief Conservative fourth-derivative velocity damping for Hydro and MHD.

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "hyperviscosity.hpp"

namespace {

parabolic::ParabolicIntegratorMode ParseHyperViscosityIntegrator(
    const std::string &block, ParameterInput *pin) {
  const std::string integrator =
      pin->GetOrAddString(block, "hyperviscosity_integrator", "explicit");
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  }
  if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<" << block << ">/hyperviscosity_integrator = '" << integrator
            << "' must be 'explicit' or 'sts'" << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------

HyperViscosity::HyperViscosity(std::string block, MeshBlockPack *pp,
                               ParameterInput *pin) :
    pmy_pack(pp),
    lapv("hypervisc_lapv", 1, 1, 1, 1, 1) {
  nu4 = pin->GetReal(block, "hyperviscosity");
  if (nu4 < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/hyperviscosity must be non-negative" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mode = ParseHyperViscosityIntegrator(block, pin);
  dtnew = std::numeric_limits<float>::max();
  if (nu4 == 0.0) {
    return;
  }

  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/hyperviscosity is not supported with SMR or AMR"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_pack->pcoord->is_special_relativistic ||
      pmy_pack->pcoord->is_general_relativistic ||
      pmy_pack->pcoord->is_dynamical_relativistic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/hyperviscosity is supported only for Newtonian "
              << "Hydro and MHD" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (indcs.ng < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/hyperviscosity requires <mesh>/nghost >= 2"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const int nmb = std::max(pmy_pack->nmb_thispack, pmy_pack->pmesh->nmb_maxperrank);
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  Kokkos::realloc(lapv, nmb, 3, ncells3, ncells2, ncells1);
}

//----------------------------------------------------------------------------------------
//! \brief Compute the explicit stability limit for -nu4*(discrete Laplacian)^2.

void HyperViscosity::NewTimeStep() {
  dtnew = std::numeric_limits<float>::max();
  if (nu4 <= 0.0) {
    return;
  }

  auto size = pmy_pack->pmb->mb_size;
  for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
    Real inv_dx2_sum = 1.0/SQR(size.h_view(m).dx1);
    if (pmy_pack->pmesh->multi_d) {
      inv_dx2_sum += 1.0/SQR(size.h_view(m).dx2);
    }
    if (pmy_pack->pmesh->three_d) {
      inv_dx2_sum += 1.0/SQR(size.h_view(m).dx3);
    }
    dtnew = std::min(dtnew, 1.0/(8.0*nu4*SQR(inv_dx2_sum)));
  }
}

//----------------------------------------------------------------------------------------
//! \brief Add conservative fourth-derivative momentum and ideal-gas energy fluxes.

void HyperViscosity::AddHyperViscousFlux(const DvceArray5D<Real> &w0,
                                         const EOS_Data &eos,
                                         DvceFaceFld5D<Real> &flx) {
  if (nu4 <= 0.0) {
    return;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int ilo = is - 1, ihi = ie + 1;
  const int jlo = multi_d ? js - 1 : js;
  const int jhi = multi_d ? je + 1 : je;
  const int klo = three_d ? ks - 1 : ks;
  const int khi = three_d ? ke + 1 : ke;
  const Real nu4_ = nu4;
  auto size = pmy_pack->pmb->mb_size;
  auto lapv_ = lapv;

  // One cached Laplacian evaluation supplies all subsequent face-gradient fluxes.
  par_for("hypervisc_lapv", DevExeSpace(), 0, nmb1, 0, 2, klo, khi, jlo, jhi,
          ilo, ihi, KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
    const int nv = IVX + n;
    const Real idx1 = 1.0/size.d_view(m).dx1;
    Real lap = (w0(m,nv,k,j,i+1) - 2.0*w0(m,nv,k,j,i) +
                w0(m,nv,k,j,i-1))*idx1*idx1;
    if (multi_d) {
      const Real idx2 = 1.0/size.d_view(m).dx2;
      lap += (w0(m,nv,k,j+1,i) - 2.0*w0(m,nv,k,j,i) +
              w0(m,nv,k,j-1,i))*idx2*idx2;
    }
    if (three_d) {
      const Real idx3 = 1.0/size.d_view(m).dx3;
      lap += (w0(m,nv,k+1,j,i) - 2.0*w0(m,nv,k,j,i) +
              w0(m,nv,k-1,j,i))*idx3*idx3;
    }
    lapv_(m,n,k,j,i) = lap;
  });

  auto flx1 = flx.x1f;
  par_for("hypervisc_x1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real coeff = 0.5*nu4_*(w0(m,IDN,k,j,i) + w0(m,IDN,k,j,i-1))/
                       size.d_view(m).dx1;
    const Real fm1 = coeff*(lapv_(m,0,k,j,i) - lapv_(m,0,k,j,i-1));
    const Real fm2 = coeff*(lapv_(m,1,k,j,i) - lapv_(m,1,k,j,i-1));
    const Real fm3 = coeff*(lapv_(m,2,k,j,i) - lapv_(m,2,k,j,i-1));
    flx1(m,IVX,k,j,i) += fm1;
    flx1(m,IVY,k,j,i) += fm2;
    flx1(m,IVZ,k,j,i) += fm3;
    if (eos.is_ideal) {
      flx1(m,IEN,k,j,i) +=
          0.5*((w0(m,IVX,k,j,i) + w0(m,IVX,k,j,i-1))*fm1 +
               (w0(m,IVY,k,j,i) + w0(m,IVY,k,j,i-1))*fm2 +
               (w0(m,IVZ,k,j,i) + w0(m,IVZ,k,j,i-1))*fm3);
    }
  });
  if (!multi_d) {
    return;
  }

  auto flx2 = flx.x2f;
  par_for("hypervisc_x2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real coeff = 0.5*nu4_*(w0(m,IDN,k,j,i) + w0(m,IDN,k,j-1,i))/
                       size.d_view(m).dx2;
    const Real fm1 = coeff*(lapv_(m,0,k,j,i) - lapv_(m,0,k,j-1,i));
    const Real fm2 = coeff*(lapv_(m,1,k,j,i) - lapv_(m,1,k,j-1,i));
    const Real fm3 = coeff*(lapv_(m,2,k,j,i) - lapv_(m,2,k,j-1,i));
    flx2(m,IVX,k,j,i) += fm1;
    flx2(m,IVY,k,j,i) += fm2;
    flx2(m,IVZ,k,j,i) += fm3;
    if (eos.is_ideal) {
      flx2(m,IEN,k,j,i) +=
          0.5*((w0(m,IVX,k,j,i) + w0(m,IVX,k,j-1,i))*fm1 +
               (w0(m,IVY,k,j,i) + w0(m,IVY,k,j-1,i))*fm2 +
               (w0(m,IVZ,k,j,i) + w0(m,IVZ,k,j-1,i))*fm3);
    }
  });
  if (!three_d) {
    return;
  }

  auto flx3 = flx.x3f;
  par_for("hypervisc_x3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real coeff = 0.5*nu4_*(w0(m,IDN,k,j,i) + w0(m,IDN,k-1,j,i))/
                       size.d_view(m).dx3;
    const Real fm1 = coeff*(lapv_(m,0,k,j,i) - lapv_(m,0,k-1,j,i));
    const Real fm2 = coeff*(lapv_(m,1,k,j,i) - lapv_(m,1,k-1,j,i));
    const Real fm3 = coeff*(lapv_(m,2,k,j,i) - lapv_(m,2,k-1,j,i));
    flx3(m,IVX,k,j,i) += fm1;
    flx3(m,IVY,k,j,i) += fm2;
    flx3(m,IVZ,k,j,i) += fm3;
    if (eos.is_ideal) {
      flx3(m,IEN,k,j,i) +=
          0.5*((w0(m,IVX,k,j,i) + w0(m,IVX,k-1,j,i))*fm1 +
               (w0(m,IVY,k,j,i) + w0(m,IVY,k-1,j,i))*fm2 +
               (w0(m,IVZ,k,j,i) + w0(m,IVZ,k-1,j,i))*fm3);
    }
  });
}
