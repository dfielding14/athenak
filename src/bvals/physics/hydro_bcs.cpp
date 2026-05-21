//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_bcs.cpp
//  \brief

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/external_gravity.hpp"

namespace {

bool UsesHydrostaticGravityBoundary(Mesh *pm) {
  for (int n=0; n<6; ++n) {
    if (pm->mesh_bcs[n] == BoundaryFlag::hydrostatic_gravity) return true;
  }
  return false;
}

KOKKOS_INLINE_FUNCTION
Real ClampExp(const Real exponent, const Real max_exponent) {
  return exp(fmin(max_exponent, fmax(-max_exponent, exponent)));
}

KOKKOS_INLINE_FUNCTION
void ApplyVelocityMode(const external_gravity::ExternalGravityData &grav,
                       const int normal, const int side,
                       Real &v1, Real &v2, Real &v3) {
  Real vn = (normal == 1) ? v1 : ((normal == 2) ? v2 : v3);

  if (grav.boundary_velocity ==
      external_gravity::BoundaryVelocityMode::zero_normal_velocity) {
    vn = 0.0;
  } else if (grav.boundary_velocity ==
             external_gravity::BoundaryVelocityMode::no_inflow_velocity) {
    if ((side < 0 && vn > 0.0) || (side > 0 && vn < 0.0)) {
      vn = 0.0;
    }
  } else if (grav.boundary_velocity ==
             external_gravity::BoundaryVelocityMode::reflect_normal_velocity) {
    vn = -vn;
  }

  if (normal == 1) {
    v1 = vn;
  } else if (normal == 2) {
    v2 = vn;
  } else {
    v3 = vn;
  }
}

KOKKOS_INLINE_FUNCTION
void FillHydrostaticGravityCell(DvceArray5D<Real> u0,
                                const external_gravity::ExternalGravityData &grav,
                                const EOS_Data &eos, const int nvar, const int nhydro,
                                const int normal, const int side,
                                const int m, const int kg, const int jg, const int ig,
                                const int ka, const int ja, const int ia,
                                const Real x1g, const Real x2g, const Real x3g,
                                const Real x1a, const Real x2a, const Real x3a) {
  const Real rho_floor = (grav.boundary_density_floor > 0.0) ?
                         grav.boundary_density_floor : eos.dfloor;
  const Real pressure_floor = (grav.boundary_pressure_floor > 0.0) ?
                              grav.boundary_pressure_floor : eos.pfloor;

  const Real rho_a = fmax(u0(m, IDN, ka, ja, ia), rho_floor);
  Real v1 = u0(m, IM1, ka, ja, ia)/rho_a;
  Real v2 = u0(m, IM2, ka, ja, ia)/rho_a;
  Real v3 = u0(m, IM3, ka, ja, ia)/rho_a;
  ApplyVelocityMode(grav, normal, side, v1, v2, v3);

  Real cs2;
  if (grav.boundary_sound_speed > 0.0) {
    cs2 = SQR(grav.boundary_sound_speed);
  } else if (eos.is_ideal) {
    const Real ekin_a = 0.5*(SQR(u0(m, IM1, ka, ja, ia)) +
                             SQR(u0(m, IM2, ka, ja, ia)) +
                             SQR(u0(m, IM3, ka, ja, ia)))/rho_a;
    const Real eint_a = fmax(u0(m, IEN, ka, ja, ia) - ekin_a,
                             pressure_floor/(eos.gamma - 1.0));
    cs2 = fmax((eos.gamma - 1.0)*eint_a/rho_a, 1.0e-20);
  } else {
    cs2 = fmax(SQR(eos.iso_cs), 1.0e-20);
  }

  const Real phi_g = external_gravity::Potential(grav, x1g, x2g, x3g);
  const Real phi_a = external_gravity::Potential(grav, x1a, x2a, x3a);
  const Real rho_g = fmax(rho_a*ClampExp(-(phi_g - phi_a)/cs2,
                                         grav.boundary_max_exponent),
                          rho_floor);
  const Real p_g = fmax(cs2*rho_g, pressure_floor);
  const Real ekin_g = 0.5*rho_g*(SQR(v1) + SQR(v2) + SQR(v3));

  u0(m, IDN, kg, jg, ig) = rho_g;
  u0(m, IM1, kg, jg, ig) = rho_g*v1;
  u0(m, IM2, kg, jg, ig) = rho_g*v2;
  u0(m, IM3, kg, jg, ig) = rho_g*v3;
  if (eos.is_ideal) {
    u0(m, IEN, kg, jg, ig) = p_g/(eos.gamma - 1.0) + ekin_g;
  }

  for (int n=nhydro; n<nvar; ++n) {
    u0(m, n, kg, jg, ig) = rho_g*u0(m, n, ka, ja, ia)/rho_a;
  }
}

} // namespace

//----------------------------------------------------------------------------------------
//! \!fn void BoundaryValues::HydroBCs()
//! \brief Apply physical boundary conditions for all Hydro variables at faces of MB which
//! are at the edge of the computational domain

void MeshBoundaryValues::HydroBCs(MeshBlockPack *ppack, DualArray2D<Real> u_in,
                                  DvceArray5D<Real> u0) {
  // loop over all MeshBlocks in this MeshBlockPack
  auto &pm = ppack->pmesh;
  auto &indcs = ppack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  auto &mb_bcs = ppack->pmb->mb_bcs;

  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int nvar = u0.extent_int(1);  // TODO(@user): 2nd index from L of in array must be NVAR
  int nmb = ppack->nmb_thispack;
  auto &size = ppack->pmb->mb_size;

  const bool use_hydrostatic_gravity = UsesHydrostaticGravityBoundary(pm);
  external_gravity::ExternalGravityData grav;
  EOS_Data eos;
  int nhydro = 0;
  if (use_hydrostatic_gravity) {
    if (ppack->phydro == nullptr || ppack->phydro->psrc == nullptr ||
        !ppack->phydro->psrc->external_gravity) {
      std::cout << "### FATAL ERROR in hydro boundary conditions" << std::endl
                << "hydrostatic_gravity boundaries require hydro with an enabled "
                << "<external_gravity> block" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    grav = ppack->phydro->psrc->external_gravity_data;
    eos = ppack->phydro->peos->eos_data;
    nhydro = ppack->phydro->nhydro;
  }

  // only apply BCs if not (periodic) or (shear_periodic)
  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic &&
      pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::shear_periodic) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    par_for("hydrobc_x1", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int n, int k, int j) {
      // apply physical boundaries to inner_x1
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,is-i-1) = -u0(m,n,k,j,is+i);
            } else {
              u0(m,n,k,j,is-i-1) =  u0(m,n,k,j,is+i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = u0(m,n,k,j,is);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = u_in.d_view(n,BoundaryFace::inner_x1);
          }
          break;
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,is-i-1) = fmin(0.0,u0(m,n,k,j,is));
            } else {
              u0(m,n  ,k,j,is-i-1) = u0(m,n,k,j,is);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = 0.0;
          }
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x1
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {  // reflect 1-velocity
              u0(m,n,k,j,ie+i+1) = -u0(m,n,k,j,ie-i);
            } else {
              u0(m,n,k,j,ie+i+1) =  u0(m,n,k,j,ie-i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = u0(m,n,k,j,ie);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = u_in.d_view(n,BoundaryFace::outer_x1);
          }
          break;
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,ie+i+1) = fmax(0.0,u0(m,n,k,j,ie));
            } else {
              u0(m,n  ,k,j,ie+i+1) = u0(m,n,k,j,ie);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = 0.0;
          }
          break;
        default:
          break;
      }
    });
  }

  if (use_hydrostatic_gravity &&
      (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::hydrostatic_gravity ||
       pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::hydrostatic_gravity)) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    int js = indcs.js;
    int ks = indcs.ks;
    par_for("hydrobc_grav_x1", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
      const Real x1min = size.d_view(m).x1min;
      const Real x1max = size.d_view(m).x1max;
      const Real x2min = size.d_view(m).x2min;
      const Real x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min;
      const Real x3max = size.d_view(m).x3max;
      const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);
      const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);

      if (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::hydrostatic_gravity) {
        const Real x1a = CellCenterX(0, nx1, x1min, x1max);
        for (int i=0; i<ng; ++i) {
          const int ig = is - i - 1;
          const Real x1g = CellCenterX(ig - is, nx1, x1min, x1max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 1, -1,
                                     m, k, j, ig, k, j, is,
                                     x1g, x2v, x3v, x1a, x2v, x3v);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::hydrostatic_gravity) {
        const Real x1a = CellCenterX(ie - is, nx1, x1min, x1max);
        for (int i=0; i<ng; ++i) {
          const int ig = ie + i + 1;
          const Real x1g = CellCenterX(ig - is, nx1, x1min, x1max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 1, 1,
                                     m, k, j, ig, k, j, ie,
                                     x1g, x2v, x3v, x1a, x2v, x3v);
        }
      }
    });
  }

  if (pm->one_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    int &js = indcs.js;
    int &je = indcs.je;
    par_for("hydrobc_x2", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int n, int k, int i) {
      // apply physical boundaries to inner_x2
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {  // reflect 2-velocity
              u0(m,n,k,js-j-1,i) = -u0(m,n,k,js+j,i);
            } else {
              u0(m,n,k,js-j-1,i) =  u0(m,n,k,js+j,i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = u0(m,n,k,js,i);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = u_in.d_view(n,BoundaryFace::inner_x2);
          }
          break;
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {
              u0(m,n,k,js-j-1,i) = fmin(0.0,u0(m,n,k,js,i));
            } else {
              u0(m,n,k,js-j-1,i) = u0(m,n,k,js,i);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = 0.0;
          }
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x2
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {  // reflect 2-velocity
              u0(m,n,k,je+j+1,i) = -u0(m,n,k,je-j,i);
            } else {
              u0(m,n,k,je+j+1,i) =  u0(m,n,k,je-j,i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = u0(m,n,k,je,i);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = u_in.d_view(n,BoundaryFace::outer_x2);
          }
          break;
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {
              u0(m,n,k,je+j+1,i) = fmax(0.0,u0(m,n,k,je,i));
            } else {
              u0(m,n,k,je+j+1,i) = u0(m,n,k,je,i);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = 0.0;
          }
          break;
        default:
          break;
      }
    });
  }

  if (use_hydrostatic_gravity &&
      (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::hydrostatic_gravity ||
       pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::hydrostatic_gravity)) {
    int &js = indcs.js;
    int &je = indcs.je;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    int is = indcs.is;
    int ks = indcs.ks;
    par_for("hydrobc_grav_x2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int i) {
      const Real x1min = size.d_view(m).x1min;
      const Real x1max = size.d_view(m).x1max;
      const Real x2min = size.d_view(m).x2min;
      const Real x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min;
      const Real x3max = size.d_view(m).x3max;
      const Real x1v = CellCenterX(i - is, nx1, x1min, x1max);
      const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);

      if (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::hydrostatic_gravity) {
        const Real x2a = CellCenterX(0, nx2, x2min, x2max);
        for (int j=0; j<ng; ++j) {
          const int jg = js - j - 1;
          const Real x2g = CellCenterX(jg - js, nx2, x2min, x2max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 2, -1,
                                     m, k, jg, i, k, js, i,
                                     x1v, x2g, x3v, x1v, x2a, x3v);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::hydrostatic_gravity) {
        const Real x2a = CellCenterX(je - js, nx2, x2min, x2max);
        for (int j=0; j<ng; ++j) {
          const int jg = je + j + 1;
          const Real x2g = CellCenterX(jg - js, nx2, x2min, x2max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 2, 1,
                                     m, k, jg, i, k, je, i,
                                     x1v, x2g, x3v, x1v, x2a, x3v);
        }
      }
    });
  }
  if (pm->two_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic) return;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  par_for("hydrobc_x3", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int n, int j, int i) {
    // apply physical boundaries to inner_x3
    switch (mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {  // reflect 3-velocity
            u0(m,n,ks-k-1,j,i) = -u0(m,n,ks+k,j,i);
          } else {
            u0(m,n,ks-k-1,j,i) =  u0(m,n,ks+k,j,i);
          }
        }
        break;
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = u0(m,n,ks,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = u_in.d_view(n,BoundaryFace::inner_x3);
        }
        break;
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {
            u0(m,n,ks-k-1,j,i) = fmin(0.0,u0(m,n,ks,j,i));
          } else {
            u0(m,n,ks-k-1,j,i) = u0(m,n,ks,j,i);
          }
        }
        break;
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = 0.0;
        }
        break;
      default:
        break;
    }

    // apply physical boundaries to outer_x3
    switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {  // reflect 3-velocity
            u0(m,n,ke+k+1,j,i) = -u0(m,n,ke-k,j,i);
          } else {
            u0(m,n,ke+k+1,j,i) =  u0(m,n,ke-k,j,i);
          }
        }
        break;
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = u0(m,n,ke,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = u_in.d_view(n,BoundaryFace::outer_x3);
        }
        break;
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {
            u0(m,n,ke+k+1,j,i) = fmax(0.0,u0(m,n,ke,j,i));
          } else {
            u0(m,n,ke+k+1,j,i) = u0(m,n,ke,j,i);
          }
        }
        break;
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = 0.0;
        }
        break;
      default:
        break;
    }
  });

  if (use_hydrostatic_gravity &&
      (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::hydrostatic_gravity ||
       pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::hydrostatic_gravity)) {
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    int is = indcs.is;
    int js = indcs.js;
    par_for("hydrobc_grav_x3", DevExeSpace(), 0,(nmb-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int j, int i) {
      const Real x1min = size.d_view(m).x1min;
      const Real x1max = size.d_view(m).x1max;
      const Real x2min = size.d_view(m).x2min;
      const Real x2max = size.d_view(m).x2max;
      const Real x3min = size.d_view(m).x3min;
      const Real x3max = size.d_view(m).x3max;
      const Real x1v = CellCenterX(i - is, nx1, x1min, x1max);
      const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);

      if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::hydrostatic_gravity) {
        const Real x3a = CellCenterX(0, nx3, x3min, x3max);
        for (int k=0; k<ng; ++k) {
          const int kg = ks - k - 1;
          const Real x3g = CellCenterX(kg - ks, nx3, x3min, x3max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 3, -1,
                                     m, kg, j, i, ks, j, i,
                                     x1v, x2v, x3g, x1v, x2v, x3a);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::hydrostatic_gravity) {
        const Real x3a = CellCenterX(ke - ks, nx3, x3min, x3max);
        for (int k=0; k<ng; ++k) {
          const int kg = ke + k + 1;
          const Real x3g = CellCenterX(kg - ks, nx3, x3min, x3max);
          FillHydrostaticGravityCell(u0, grav, eos, nvar, nhydro, 3, 1,
                                     m, kg, j, i, ke, j, i,
                                     x1v, x2v, x3g, x1v, x2v, x3a);
        }
      }
    });
  }

  return;
}
