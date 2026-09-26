//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.cpp
//  \brief Implements functions for Resistivity class.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "resistivity.hpp"
#include "current_density.hpp"

namespace {

[[noreturn]] void ResistivityError(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << ": " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

parabolic::DiffusionSelection ParseResistivityIntegrator(ParameterInput *pin) {
  std::string integrator = pin->GetOrAddString("mhd", "resistivity_integrator",
                                               "explicit");
  if (integrator == "explicit") {
    return parabolic::DiffusionSelection::explicit_only;
  }
  if (integrator == "sts") {
    return parabolic::DiffusionSelection::sts_only;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<mhd>/resistivity_integrator = '" << integrator
            << "' must be 'explicit' or 'sts'" << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls Resistivity base class constructor

Resistivity::Resistivity(MeshBlockPack *pp, ParameterInput *pin) :
  pmy_pack(pp) {
  // Read parameters for Ohmic diffusion (if any)
  eta_ohm = pin->GetReal("mhd","ohmic_resistivity");
  mode = ParseResistivityIntegrator(pin);
  std::string model = pin->GetOrAddString("mhd", "resistivity_model", "constant");
  if (model != "constant" && model != "current_limited") {
    ResistivityError("resistivity_model must be constant or current_limited");
  }
  current_limited = (model == "current_limited");
  if (current_limited) {
    auto &p = current_limited_params;
    p.eta0 = eta_ohm;
    p.eta_max = pin->GetReal("mhd", "eta_max");
    if (!std::isfinite(p.eta0) || !std::isfinite(p.eta_max) ||
        p.eta0 <= 0.0 || p.eta_max < p.eta0) {
      ResistivityError("current_limited requires finite 0 < ohmic_resistivity <= eta_max");
    }
    std::string method = pin->GetOrAddString("mhd", "b_rec_method", "constant");
    if (method != "constant" && method != "jump") {
      ResistivityError("b_rec_method must be constant or jump");
    }
    b_rec_jump = (method == "jump");
    p.d_i = pin->GetReal("mhd", "d_i");
    if (b_rec_jump) {
      b_rec_radius = pin->GetReal("mhd", "b_rec_radius");
      b_rec_floor = pin->GetReal("mhd", "b_rec_floor");
      if (!std::isfinite(b_rec_radius) || !std::isfinite(b_rec_floor) ||
          b_rec_radius <= 0.0 || b_rec_floor <= 0.0) {
        ResistivityError("jump requires finite positive b_rec_radius and b_rec_floor");
      }
      p.b_rec = b_rec_floor;
      auto *pm = pmy_pack->pmesh;
      const auto &indcs = pm->mb_indcs;
      const Real dx[3] = {pm->mesh_size.dx1, pm->mesh_size.dx2, pm->mesh_size.dx3};
      const int nx[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
      const int ndim = pm->three_d ? 3 : (pm->multi_d ? 2 : 1);
      for (int a=0; a<ndim; ++a) {
        const Real finest_dx = std::ldexp(dx[a], pm->root_level-pm->max_level);
        const Real reach = std::ceil(b_rec_radius/finest_dx);
        // An edge also uses an adjacent coefficient ghost cell. Validate against
        // future allowed AMR levels, not only the blocks initially present.
        if (!std::isfinite(reach) || reach > indcs.ng-1) {
          ResistivityError("jump requires nghost >= ceil(b_rec_radius/dx_finest) + 1");
        }
        if (nx[a] < (pm->multilevel ? 2 : 1)*indcs.ng) {
          ResistivityError("jump requires meshblock width >= nghost (2*nghost with SMR/AMR)");
        }
      }
    } else {
      p.b_rec = pin->GetReal("mhd", "b_rec");
    }
    if (!std::isfinite(p.d_i) || !std::isfinite(p.b_rec) ||
        p.d_i <= 0.0 || p.b_rec <= 0.0) {
      ResistivityError("current_limited requires finite positive d_i and b_rec");
    }
    // The MHD constructor has built its EOS, but pmy_pack->pmhd is not assigned yet.
    const Real dfloor = pin->GetReal("mhd", "dfloor");
    if (!std::isfinite(dfloor) || dfloor <= 0.0) {
      ResistivityError("current_limited requires a finite positive mhd/dfloor");
    }
    auto *coord = pmy_pack->pcoord;
    if (coord->is_special_relativistic || coord->is_general_relativistic ||
        coord->is_dynamical_relativistic) {
      ResistivityError("current_limited resistivity requires Newtonian MHD");
    }
    if (pin->DoesParameterExist("mhd", "eta_ad") &&
        pin->GetReal("mhd", "eta_ad") != 0.0) {
      ResistivityError("current_limited cannot be combined with eta_ad");
    }
    const Real root_ratio = std::sqrt(p.eta0)/std::sqrt(p.eta_max);
    p.q_star = 1.0-root_ratio;
    p.eta_star = std::sqrt(p.eta0)*std::sqrt(p.eta_max);
    pin->SetString("mhd", "resistivity_species", "ion");
  }
  if (mode == parabolic::DiffusionSelection::sts_only &&
      (!std::isfinite(eta_ohm) || eta_ohm <= 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "STS requires positive constant Ohmic resistivity" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  NewTimeStep();
  if (current_limited && global_variable::my_rank == 0) {
    const auto &p = current_limited_params;
    Real max_dx = 0.0;
    const auto &size = pmy_pack->pmb->mb_size;
    for (int m=0; m<pmy_pack->nmb_thispack; ++m) {
      max_dx = std::max(max_dx, size.h_view(m).dx1);
      if (pmy_pack->pmesh->multi_d) max_dx = std::max(max_dx, size.h_view(m).dx2);
      if (pmy_pack->pmesh->three_d) max_dx = std::max(max_dx, size.h_view(m).dx3);
    }
    std::cout << "Resistivity: current_limited, species=ion, b_rec_method="
              << (b_rec_jump ? "jump (experimental, frozen per cycle)" : "constant")
              << "\n  q_star=" << p.q_star << " eta_star=" << p.eta_star
              << "\n  rank-0 parabolic dt bound (before CFL)=" << dtnew
              << " min cells per d_i(rho=1)=" << p.d_i/max_dx << std::endl;
    if (p.d_i/max_dx < 4.0) {
      std::cout << "### WARNING: fewer than 4 cells per d_i(rho=1) on rank 0"
                << std::endl;
    }
    if (b_rec_jump) {
      const auto &indcs = pmy_pack->pmesh->mb_indcs;
      const int nx[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
      const int ndim = pmy_pack->pmesh->three_d ? 3 : (pmy_pack->pmesh->multi_d ? 2 : 1);
      Real field_ratio = 1.0, cache_cells = 1.0;
      for (int a=0; a<ndim; ++a) {
        field_ratio *= static_cast<Real>(nx[a]+2*indcs.ng)/(nx[a]+4);
        cache_cells *= nx[a]+2;
      }
      std::cout << "  jump radius=" << b_rec_radius << " floor=" << b_rec_floor
                << " nghost=" << indcs.ng
                << "\n  cell-array volume relative to nghost=2: " << field_ratio
                << "; cached B_rec cells per block: " << cache_cells << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Rebuild the experimental jump estimate from valid cell-centered field ghosts.
//! Called after initialization/regridding and before a complete split evolution cycle.

void Resistivity::PrepareBRec() {
  if (!b_rec_jump) return;
  auto *pm = pmy_pack->pmesh;
  const auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int n1 = indcs.nx1+2;
  const int n2 = pm->multi_d ? indcs.nx2+2 : 1;
  const int n3 = pm->three_d ? indcs.nx3+2 : 1;
  auto &cache = current_limited_params.b_rec_cells;
  if (cache.extent_int(0) < nmb) {
    Kokkos::realloc(cache, std::max(nmb, pm->nmb_maxperrank), n3, n2, n1);
  }
  const int offset = indcs.ng-1;
  current_limited_params.b_rec_offset = offset;
  auto bcc = pmy_pack->pmhd->bcc0;
  auto size = pmy_pack->pmb->mb_size;
  const Real radius = b_rec_radius, floor = b_rec_floor;
  const bool multi_d = pm->multi_d, three_d = pm->three_d;
  par_for("reconnecting_field_jump", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    cache(m,k,j,i) = ::current_limited::JumpEstimate(bcc, size.d_view(m), radius, floor,
        multi_d, three_d, m, k+(three_d ? offset : 0), j+(multi_d ? offset : 0), i+offset);
  });
}

//----------------------------------------------------------------------------------------
//! \brief Refresh the resistive stability bound after evolution or mesh refinement.

void Resistivity::NewTimeStep() {
  dtnew = std::numeric_limits<float>::max();
  if (eta_ohm == 0.0) return;
  const Real eta_bound = current_limited ? current_limited_params.eta_max : eta_ohm;
  auto size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  for (int m=0; m<(pmy_pack->nmb_thispack); ++m) {
    dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx1)/eta_bound);
    if (pmy_pack->pmesh->multi_d) {dtnew = std::min(dtnew,fac*SQR(size.h_view(m).dx2)/eta_bound);}
    if (pmy_pack->pmesh->three_d) {dtnew = std::min(dtnew,fac*SQR(size.h_view(m).dx3)/eta_bound);}
  }
}

//----------------------------------------------------------------------------------------
// Resistivity destructor

Resistivity::~Resistivity() {
}

//----------------------------------------------------------------------------------------
//! \fn OhmicEField()
//  \brief Adds electric field from Ohmic resistivity to corner-centered electric field
//  Using Ohm's Law to compute the electric field:  E + (v x B) = \eta J, then
//    E_{inductive} = - (v x B)  [computed in the MHD Riemann solver]
//    E_{resistive} = \eta J     [computed in this function]

void Resistivity::OhmicEField(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld) {
  if (current_limited && current_limited_params.eta_max != eta_ohm) {
    CurrentLimitedEField(b0, efld);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int nmb1 = pmy_pack->nmb_thispack - 1;

  //---- 1-D problem:
  //  copy face-centered E-fields to edges and return.
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;
    auto eta_o = eta_ohm;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm1", DevExeSpace(), scr_size, scr_level, 0, nmb1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(member, m, ks, js, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields
      par_for_inner(member, is, ie+1, [&](const int i) {
        e2(m,ks,  js  ,i) += eta_o*j2(i);
        e2(m,ke+1,js  ,i) += eta_o*j2(i);
        e3(m,ks  ,js  ,i) += eta_o*j3(i);
        e3(m,ks  ,je+1,i) += eta_o*j3(i);
      });
    });
    return;
  }

  //---- 2-D problem:
  if (pmy_pack->pmesh->two_d) {
    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;
    auto eta_o = eta_ohm;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm2", DevExeSpace(), scr_size, scr_level, 0, nmb1, js, je+1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(member, m, ks, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields
      par_for_inner(member, is, ie+1, [&](const int i) {
        e1(m,ks,  j,i) += eta_o*j1(i);
        e1(m,ke+1,j,i) += eta_o*j1(i);
        e2(m,ks,  j,i) += eta_o*j2(i);
        e2(m,ke+1,j,i) += eta_o*j2(i);
        e3(m,ks  ,j,i) += eta_o*j3(i);
      });
    });
    return;
  }

  //---- 3-D problem:

  // capture class variables for the kernels
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto eta_o = eta_ohm;

  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

  par_for_outer("ohm3", DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke+1, js, je+1,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

    CurrentDensity(member, m, k, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

    // Add E_{resistive} = \eta J to corner-centered electric fields
    par_for_inner(member, is, ie+1, [&](const int i) {
      e1(m,k,j,i) += eta_o*j1(i);
      e2(m,k,j,i) += eta_o*j2(i);
      e3(m,k,j,i) += eta_o*j3(i);
    });
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn OhmicEnergyFlux()
//  \brief Adds Poynting flux from Ohmic resistivity to energy flux
//  Total energy equation is dE/dt = - Div(F) where F = (E X B) = \eta (J X B)


void Resistivity::OhmicEnergyFlux(const DvceFaceFld4D<Real> &b,
                                  DvceFaceFld5D<Real> &flx) {
  if (current_limited && current_limited_params.eta_max != eta_ohm) {
    CurrentLimitedEnergyFlux(b, flx);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  Real qa = 0.25*eta_ohm;

  //------------------------------
  // energy fluxes in x1-direction

  auto &flx1 = flx.x1f;
  par_for("ohm_heat1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j2k   = -(b.x3f(m,k  ,j,i) - b.x3f(m,k  ,j,i-1))/size.d_view(m).dx1;
    Real j2kp1 = -(b.x3f(m,k+1,j,i) - b.x3f(m,k+1,j,i-1))/size.d_view(m).dx1;

    Real j3j   = (b.x2f(m,k,j  ,i) - b.x2f(m,k,j  ,i-1))/size.d_view(m).dx1;
    Real j3jp1 = (b.x2f(m,k,j+1,i) - b.x2f(m,k,j+1,i-1))/size.d_view(m).dx1;

    if (multi_d) {
      j3j   -= (b.x1f(m,k,j  ,i) - b.x1f(m,k,j-1,i))/size.d_view(m).dx2;
      j3jp1 -= (b.x1f(m,k,j+1,i) - b.x1f(m,k,j  ,i))/size.d_view(m).dx2;
    }
    if (three_d) {
      j2k   += (b.x1f(m,k  ,j,i) - b.x1f(m,k-1,j,i))/size.d_view(m).dx3;
      j2kp1 += (b.x1f(m,k+1,j,i) - b.x1f(m,k  ,j,i))/size.d_view(m).dx3;
    }

    // flx1 = (E X B)_{1} =  ((\eta J) X B)_{1} = \eta (J2*B3 - J3*B2)
    flx1(m,IEN,k,j,i) += qa*(j2k  *(b.x3f(m,k  ,j  ,i) + b.x3f(m,k  ,j  ,i-1)) +
                             j2kp1*(b.x3f(m,k+1,j  ,i) + b.x3f(m,k+1,j  ,i-1)) -
                             j3j  *(b.x2f(m,k  ,j  ,i) + b.x2f(m,k  ,j  ,i-1)) -
                             j3jp1*(b.x2f(m,k  ,j+1,i) + b.x2f(m,k  ,j+1,i-1)));
  });
  if (pmy_pack->pmesh->one_d) {return;}

  //------------------------------
  // energy fluxes in x2-direction

  auto &flx2 = flx.x2f;
  par_for("ohm_heat2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j1k   = (b.x3f(m,k  ,j,i) - b.x3f(m,k  ,j-1,i))/size.d_view(m).dx2;
    Real j1kp1 = (b.x3f(m,k+1,j,i) - b.x3f(m,k+1,j-1,i))/size.d_view(m).dx2;

    Real j3i   = (b.x2f(m,k,j,i  ) - b.x2f(m,k,j  ,i-1))/size.d_view(m).dx1
               - (b.x1f(m,k,j,i  ) - b.x1f(m,k,j-1,i  ))/size.d_view(m).dx2;
    Real j3ip1 = (b.x2f(m,k,j,i+1) - b.x2f(m,k,j  ,i  ))/size.d_view(m).dx1
               - (b.x1f(m,k,j,i+1) - b.x1f(m,k,j-1,i+1))/size.d_view(m).dx2;

    if (three_d) {
      j1k   -= (b.x2f(m,k  ,j,i) - b.x2f(m,k-1,j,i))/size.d_view(m).dx3;
      j1kp1 -= (b.x2f(m,k+1,j,i) - b.x2f(m,k  ,j,i))/size.d_view(m).dx3;
    }

    // E2 = \eta (J X B)_{2} = \eta (J3*B1 - J1*B3)
    flx2(m,IEN,k,j,i) += qa*(j3i  *(b.x1f(m,k  ,j,i  ) + b.x1f(m,k  ,j-1,i  )) +
                             j3ip1*(b.x1f(m,k  ,j,i+1) + b.x1f(m,k  ,j-1,i+1)) -
                             j1k  *(b.x3f(m,k  ,j,i  ) + b.x3f(m,k  ,j-1,i  )) -
                             j1kp1*(b.x3f(m,k+1,j,i  ) + b.x3f(m,k+1,j-1,i  )));
  });
  if (pmy_pack->pmesh->two_d) {return;}

  //------------------------------
  // energy fluxes in x3-direction

  auto &flx3 = flx.x3f;
  par_for("ohm_heat3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j1j   = (b.x3f(m,k,j  ,i) - b.x3f(m,k  ,j-1,i))/size.d_view(m).dx2
               - (b.x2f(m,k,j  ,i) - b.x2f(m,k-1,j  ,i))/size.d_view(m).dx3;
    Real j1jp1 = (b.x3f(m,k,j+1,i) - b.x3f(m,k  ,j  ,i))/size.d_view(m).dx2
               - (b.x2f(m,k,j+1,i) - b.x2f(m,k-1,j+1,i))/size.d_view(m).dx3;

    Real j2i   = -(b.x3f(m,k,j,i  ) - b.x3f(m,k  ,j,i-1))/size.d_view(m).dx1
                + (b.x1f(m,k,j,i  ) - b.x1f(m,k-1,j,i  ))/size.d_view(m).dx3;
    Real j2ip1 = -(b.x3f(m,k,j,i+1) - b.x3f(m,k  ,j,i  ))/size.d_view(m).dx1
                + (b.x1f(m,k,j,i+1) - b.x1f(m,k-1,j,i+1))/size.d_view(m).dx3;

    // E2 = \eta (J X B)_{2} = \eta (J1*B2 - J2*B1)
    flx3(m,IEN,k,j,i) += qa*(j1j  *(b.x2f(m,k,j  ,i  ) + b.x2f(m,k-1,j  ,i  )) +
                             j1jp1*(b.x2f(m,k,j+1,i  ) + b.x2f(m,k-1,j+1,i  )) -
                             j2i  *(b.x1f(m,k,j  ,i  ) + b.x1f(m,k-1,j  ,i  )) -
                             j2ip1*(b.x1f(m,k,j  ,i+1) + b.x1f(m,k-1,j  ,i+1)));
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \brief Add current-dependent eta J on each component's native CT edges.

void Resistivity::CurrentLimitedEField(const DvceFaceFld4D<Real> &b,
                                       DvceEdgeFld4D<Real> &efld) {
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto w = pmy_pack->pmhd->w0;
  auto size = pmy_pack->pmb->mb_size;
  const auto p = current_limited_params;
  const Real dfloor = pmy_pack->pmhd->peos->eos_data.dfloor;

  for (int c=0; c<3; ++c) {
    if (c == 0 && !multi_d) continue;
    auto e = (c == 0) ? efld.x1e : ((c == 1) ? efld.x2e : efld.x3e);
    const int iu = indcs.ie + (c != 0);
    const int ju = indcs.je + (c != 1);
    const int ku = indcs.ke + (c != 2);
    par_for("current_limited_emf", DevExeSpace(), 0, nmb-1,
            indcs.ks, ku, indcs.js, ju, indcs.is, iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const auto state = ::current_limited::EdgeState(
          b, w, size.d_view(m), p, dfloor, multi_d, three_d, c, m, k, j, i);
      const Real current = (c == 0) ? state.j1 : ((c == 1) ? state.j2 : state.j3);
      e(m,k,j,i) += state.eta*current;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \brief Resistive Poynting flux using the same edge coefficients as the CT update.
//! Each face averages its four bounding edge products, as in OhmicEnergyFlux.

void Resistivity::CurrentLimitedEnergyFlux(const DvceFaceFld4D<Real> &b,
                                           DvceFaceFld5D<Real> &flx) {
  const auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto w = pmy_pack->pmhd->w0;
  auto size = pmy_pack->pmb->mb_size;
  const auto p = current_limited_params;
  const Real dfloor = pmy_pack->pmhd->peos->eos_data.dfloor;
  const int ndim = three_d ? 3 : (multi_d ? 2 : 1);

  for (int a=0; a<ndim; ++a) {
    const int c1 = (a+1)%3, c2 = (a+2)%3;
    auto flux = (a == 0) ? flx.x1f : ((a == 1) ? flx.x2f : flx.x3f);
    auto b1 = (c1 == 0) ? b.x1f : ((c1 == 1) ? b.x2f : b.x3f);
    auto b2 = (c2 == 0) ? b.x1f : ((c2 == 1) ? b.x2f : b.x3f);
    const int iu = indcs.ie + (a == 0);
    const int ju = indcs.je + (a == 1);
    const int ku = indcs.ke + (a == 2);
    par_for("current_limited_energy_flux", DevExeSpace(), 0, nmb-1,
            indcs.ks, ku, indcs.js, ju, indcs.is, iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real poynting = 0.0;
      for (int side=0; side<2; ++side) {
        int edge1[3] = {i, j, k}, edge2[3] = {i, j, k};
        edge1[c2] += side;
        edge2[c1] += side;
        const auto s1 = ::current_limited::EdgeState(
            b, w, size.d_view(m), p, dfloor, multi_d, three_d, c1, m,
            edge1[2], edge1[1], edge1[0]);
        const auto s2 = ::current_limited::EdgeState(
            b, w, size.d_view(m), p, dfloor, multi_d, three_d, c2, m,
            edge2[2], edge2[1], edge2[0]);
        const Real j1 = (c1 == 0) ? s1.j1 : ((c1 == 1) ? s1.j2 : s1.j3);
        const Real j2 = (c2 == 0) ? s2.j1 : ((c2 == 1) ? s2.j2 : s2.j3);
        Real bc2 = b2(m,edge1[2],edge1[1],edge1[0]);
        Real bc1 = b1(m,edge2[2],edge2[1],edge2[0]);
        --edge1[a];
        --edge2[a];
        bc2 += b2(m,edge1[2],edge1[1],edge1[0]);
        bc1 += b1(m,edge2[2],edge2[1],edge2[0]);
        poynting += s1.eta*j1*bc2 - s2.eta*j2*bc1;
      }
      flux(m,IEN,k,j,i) += 0.25*poynting;
    });
  }
}
