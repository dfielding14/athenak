//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file viscosity.cpp
//  \brief Implements functions for Viscosity class. This includes isotropic shear
//  viscosity in a Newtonian fluid (in which stress is proportional to shear).
//  Viscosity may be added to Hydro and/or MHD independently.

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <string> // string

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "viscosity.hpp"

namespace {

// Evaluate the strain norm at a face using the same gradients as the stress.
KOKKOS_INLINE_FUNCTION
Real SmagorinskyNu(const DvceArray5D<Real> &w, int m, int k, int j, int i,
                   const RegionSize &size, Real cs, bool x1_face) {
  const int it = x1_face ? 0 : 1;
  const int jt = 1 - it;
  const int i_lo = i - jt;
  const int j_lo = j - it;
  const Real dn = x1_face ? size.dx1 : size.dx2;
  const Real dt = x1_face ? size.dx2 : size.dx1;
  const int vn = x1_face ? IVX : IVY;
  const int vt = x1_face ? IVY : IVX;
  const auto normal = [&](int v) {
    return (w(m,v,k,j,i) - w(m,v,k,j_lo,i_lo))/dn;
  };
  const auto transverse = [&](int v) {
    return ((w(m,v,k,j+jt,i+it) + w(m,v,k,j_lo+jt,i_lo+it)) -
            (w(m,v,k,j-jt,i-it) + w(m,v,k,j_lo-jt,i_lo-it)))/(4.0*dt);
  };
  const Real dvn_dn = normal(vn), dvt_dt = transverse(vt);
  const Real shear = normal(vt) + transverse(vn);
  const Real dvz_dn = normal(IVZ), dvz_dt = transverse(IVZ);
  const Real strain_sq = 2.0*(SQR(dvn_dn) + SQR(dvt_dt)) +
                         SQR(shear) + SQR(dvz_dn) + SQR(dvz_dt);
  return SQR(cs)*size.dx1*size.dx2*sqrt(strain_sq);
}

}  // namespace

//----------------------------------------------------------------------------------------
// ctor:
// Note first argument passes string ("hydro" or "mhd") denoting in which class this
// object is being constructed, and therefore which <block> in the input file from which
// the parameters are read.

Viscosity::Viscosity(std::string block, MeshBlockPack *pp,
                     ParameterInput *pin) :
  pmy_pack(pp) {
  nu_iso = pin->DoesParameterExist(block,"viscosity") ?
           pin->GetReal(block,"viscosity") : 0.0;
  if (block == "hydro" && pin->DoesParameterExist(block,"smagorinsky_cs")) {
    smagorinsky_cs = pin->GetReal(block,"smagorinsky_cs");
    if (!std::isfinite(smagorinsky_cs) || smagorinsky_cs < 0.0) {
      std::cout << "### FATAL ERROR: <hydro>/smagorinsky_cs must be finite and "
                << "non-negative." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (smagorinsky_cs > 0.0 &&
        (!pp->pmesh->two_d || pp->pcoord->is_special_relativistic ||
         pp->pcoord->is_general_relativistic || pp->pcoord->is_dynamical_relativistic)) {
      std::cout << "### FATAL ERROR: Smagorinsky viscosity requires 2D Newtonian Hydro."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (smagorinsky_cs > 0.0 && (!std::isfinite(nu_iso) || nu_iso < 0.0)) {
      std::cout << "### FATAL ERROR: <hydro>/viscosity must be finite and non-negative "
                << "with Smagorinsky viscosity." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // viscous timestep on MeshBlock(s) in this pack
  dtnew = std::numeric_limits<float>::max();
  if (nu_iso == 0.0 || smagorinsky_cs > 0.0) {return;}
  auto size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pp->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pp->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  for (int m=0; m<(pp->nmb_thispack); ++m) {
    dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx1)/nu_iso);
    if (pp->pmesh->multi_d) {dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx2)/nu_iso);}
    if (pp->pmesh->three_d) {dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx3)/nu_iso);}
  }
}

//----------------------------------------------------------------------------------------
// Viscosity destructor

Viscosity::~Viscosity() {
}

//----------------------------------------------------------------------------------------
//  \brief Update the explicit diffusion timestep using the current face viscosities.

void Viscosity::NewTimeStep(const DvceArray5D<Real> &w) {
  if (smagorinsky_cs == 0.0) {return;}

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ks = indcs.ks, ie = indcs.ie, je = indcs.je;
  const Real cs = smagorinsky_cs, nu = nu_iso;
  auto size = pmy_pack->pmb->mb_size;
  using Policy = Kokkos::MDRangePolicy<DevExeSpace, Kokkos::Rank<3>>;
  Real max_rate;
  Kokkos::parallel_reduce("SmagorinskyTimeStep",
      Policy({0, indcs.js, indcs.is}, {pmy_pack->nmb_thispack, je+2, ie+2}),
      KOKKOS_LAMBDA(int m, int j, int i, Real &max_local) {
    const auto &block_size = size.d_view(m);
    Real nu_sgs = 0.0;
    if (j <= je) {
      nu_sgs = SmagorinskyNu(w, m, ks, j, i, block_size, cs, true);
    }
    if (i <= ie) {
      nu_sgs = fmax(nu_sgs, SmagorinskyNu(w, m, ks, j, i, block_size, cs, false));
    }
    const Real h = fmin(block_size.dx1, block_size.dx2);
    max_local = fmax(max_local, (nu + 2.0*nu_sgs)/SQR(h));
  }, Kokkos::Max<Real>(max_rate));

  // The tangent viscosity doubles in pure shear, and normal stress adds a 4/3 factor.
  dtnew = (max_rate > 0.0) ? 0.1875/max_rate : std::numeric_limits<float>::max();
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsoViscousFlux
//  \brief Adds viscous fluxes to face-centered fluxes of conserved variables

void Viscosity::IsotropicViscousFlux(const DvceArray5D<Real> &w0, const Real nu_iso,
  const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  const Real cs = smagorinsky_cs;
  if (nu_iso == 0.0 && cs == 0.0) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  //--------------------------------------------------------------------------------------
  // fluxes in x1-direction

  int scr_level = 0;
  size_t scr_size = (ScrArray1D<Real>::shmem_size(ncells1)) * 3;
  auto flx1 = flx.x1f;

  par_for_outer("visc1",DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> fvx(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvy(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvz(member.team_scratch(scr_level), ncells1);

    // Add [2(dVx/dx)-(2/3)dVx/dx, dVy/dx, dVz/dx]
    par_for_inner(member, is, ie+1, [&](const int i) {
      fvx(i) = 4.0*(w0(m,IVX,k,j,i) - w0(m,IVX,k,j,i-1))/(3.0*size.d_view(m).dx1);
      fvy(i) =     (w0(m,IVY,k,j,i) - w0(m,IVY,k,j,i-1))/size.d_view(m).dx1;
      fvz(i) =     (w0(m,IVZ,k,j,i) - w0(m,IVZ,k,j,i-1))/size.d_view(m).dx1;
    });

    // In 2D/3D Add [(-2/3)dVy/dy, dVx/dy, 0]
    if (multi_d) {
      par_for_inner(member, is, ie+1, [&](const int i) {
        fvx(i) -= ((w0(m,IVY,k,j+1,i) + w0(m,IVY,k,j+1,i-1)) -
                   (w0(m,IVY,k,j-1,i) + w0(m,IVY,k,j-1,i-1)))/(6.0*size.d_view(m).dx2);
        fvy(i) += ((w0(m,IVX,k,j+1,i) + w0(m,IVX,k,j+1,i-1)) -
                   (w0(m,IVX,k,j-1,i) + w0(m,IVX,k,j-1,i-1)))/(4.0*size.d_view(m).dx2);
      });
    }

    // In 3D Add [(-2/3)dVz/dz, 0,  dVx/dz]
    if (three_d) {
      par_for_inner(member, is, ie+1, [&](const int i) {
        fvx(i) -= ((w0(m,IVZ,k+1,j,i) + w0(m,IVZ,k+1,j,i-1)) -
                   (w0(m,IVZ,k-1,j,i) + w0(m,IVZ,k-1,j,i-1)))/(6.0*size.d_view(m).dx3);
        fvz(i) += ((w0(m,IVX,k+1,j,i) + w0(m,IVX,k+1,j,i-1)) -
                   (w0(m,IVX,k-1,j,i) + w0(m,IVX,k-1,j,i-1)))/(4.0*size.d_view(m).dx3);
      });
    }

    // Sum viscous fluxes into fluxes of conserved variables; including energy fluxes
    par_for_inner(member, is, ie+1, [&](const int i) {
      Real nu = nu_iso;
      if (cs > 0.0) {nu += SmagorinskyNu(w0, m, k, j, i, size.d_view(m), cs, true);}
      Real nud = 0.5*nu*(w0(m,IDN,k,j,i) + w0(m,IDN,k,j,i-1));
      flx1(m,IVX,k,j,i) -= nud*fvx(i);
      flx1(m,IVY,k,j,i) -= nud*fvy(i);
      flx1(m,IVZ,k,j,i) -= nud*fvz(i);
      if (eos.is_ideal) {
        flx1(m,IEN,k,j,i) -= 0.5*nud*((w0(m,IVX,k,j,i-1) + w0(m,IVX,k,j,i))*fvx(i) +
                                      (w0(m,IVY,k,j,i-1) + w0(m,IVY,k,j,i))*fvy(i) +
                                      (w0(m,IVZ,k,j,i-1) + w0(m,IVZ,k,j,i))*fvz(i));
      }
    });
  });
  if (pmy_pack->pmesh->one_d) {return;}

  //--------------------------------------------------------------------------------------
  // fluxes in x2-direction

  auto flx2 = flx.x2f;

  par_for_outer("visc2",DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke, js, je+1,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> fvx(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvy(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvz(member.team_scratch(scr_level), ncells1);

    // Add [(dVx/dy+dVy/dx), 2(dVy/dy)-(2/3)(dVx/dx+dVy/dy), dVz/dy]
    par_for_inner(member, is, ie, [&](const int i) {
      fvx(i) = (w0(m,IVX,k,j,i  ) - w0(m,IVX,k,j-1,i  ))/size.d_view(m).dx2 +
              ((w0(m,IVY,k,j,i+1) + w0(m,IVY,k,j-1,i+1)) -
               (w0(m,IVY,k,j,i-1) + w0(m,IVY,k,j-1,i-1)))/(4.0*size.d_view(m).dx1);
      fvy(i) = (w0(m,IVY,k,j,i) - w0(m,IVY,k,j-1,i))*4.0/(3.0*size.d_view(m).dx2) -
              ((w0(m,IVX,k,j,i+1) + w0(m,IVX,k,j-1,i+1)) -
               (w0(m,IVX,k,j,i-1) + w0(m,IVX,k,j-1,i-1)))/(6.0*size.d_view(m).dx1);
      fvz(i) = (w0(m,IVZ,k,j,i  ) - w0(m,IVZ,k,j-1,i  ))/size.d_view(m).dx2;
    });

    // In 3D Add [0, (-2/3)dVz/dz, dVy/dz]
    if (three_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        fvy(i) -= ((w0(m,IVZ,k+1,j,i) + w0(m,IVZ,k+1,j-1,i)) -
                   (w0(m,IVZ,k-1,j,i) + w0(m,IVZ,k-1,j-1,i)))/(6.0*size.d_view(m).dx3);
        fvz(i) += ((w0(m,IVY,k+1,j,i) + w0(m,IVY,k+1,j-1,i)) -
                   (w0(m,IVY,k-1,j,i) + w0(m,IVY,k-1,j-1,i)))/(4.0*size.d_view(m).dx3);
      });
    }

    // Sum viscous fluxes into fluxes of conserved variables; including energy fluxes
    par_for_inner(member, is, ie, [&](const int i) {
      Real nu = nu_iso;
      if (cs > 0.0) {nu += SmagorinskyNu(w0, m, k, j, i, size.d_view(m), cs, false);}
      Real nud = 0.5*nu*(w0(m,IDN,k,j,i) + w0(m,IDN,k,j-1,i));
      flx2(m,IVX,k,j,i) -= nud*fvx(i);
      flx2(m,IVY,k,j,i) -= nud*fvy(i);
      flx2(m,IVZ,k,j,i) -= nud*fvz(i);
      if (eos.is_ideal) {
        flx2(m,IEN,k,j,i) -= 0.5*nud*((w0(m,IVX,k,j-1,i) + w0(m,IVX,k,j,i))*fvx(i) +
                                      (w0(m,IVY,k,j-1,i) + w0(m,IVY,k,j,i))*fvy(i) +
                                      (w0(m,IVZ,k,j-1,i) + w0(m,IVZ,k,j,i))*fvz(i));
      }
    });
  });
  if (pmy_pack->pmesh->two_d) {return;}

  //--------------------------------------------------------------------------------------
  // fluxes in x3-direction

  auto flx3 = flx.x3f;

  par_for_outer("visc3",DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke+1, js, je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> fvx(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvy(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> fvz(member.team_scratch(scr_level), ncells1);

    // Add [(dVx/dz+dVz/dx), (dVy/dz+dVz/dy), 2(dVz/dz)-(2/3)(dVx/dx+dVy/dy+dVz/dz)]
    par_for_inner(member, is, ie, [&](const int i) {
      fvx(i) = (w0(m,IVX,k,j,i  ) - w0(m,IVX,k-1,j,i  ))/size.d_view(m).dx3 +
              ((w0(m,IVZ,k,j,i+1) + w0(m,IVZ,k-1,j,i+1)) -
               (w0(m,IVZ,k,j,i-1) + w0(m,IVZ,k-1,j,i-1)))/(4.0*size.d_view(m).dx1);
      fvy(i) = (w0(m,IVY,k,j,i  ) - w0(m,IVY,k-1,j,i  ))/size.d_view(m).dx3 +
              ((w0(m,IVZ,k,j+1,i) + w0(m,IVZ,k-1,j+1,i)) -
               (w0(m,IVZ,k,j-1,i) + w0(m,IVZ,k-1,j-1,i)))/(4.0*size.d_view(m).dx2);
      fvz(i) = (w0(m,IVZ,k,j,i) - w0(m,IVZ,k-1,j,i))*4.0/(3.0*size.d_view(m).dx3) -
              ((w0(m,IVX,k,j,i+1) + w0(m,IVX,k-1,j,i+1)) -
               (w0(m,IVX,k,j,i-1) + w0(m,IVX,k-1,j,i-1)))/(6.0*size.d_view(m).dx1) -
              ((w0(m,IVY,k,j+1,i) + w0(m,IVY,k-1,j+1,i)) -
               (w0(m,IVY,k,j-1,i) + w0(m,IVY,k-1,j-1,i)))/(6.0*size.d_view(m).dx2);
    });

    // Sum viscous fluxes into fluxes of conserved variables; including energy fluxes
    par_for_inner(member, is, ie, [&](const int i) {
      Real nud = 0.5*nu_iso*(w0(m,IDN,k,j,i) + w0(m,IDN,k-1,j,i));
      flx3(m,IVX,k,j,i) -= nud*fvx(i);
      flx3(m,IVY,k,j,i) -= nud*fvy(i);
      flx3(m,IVZ,k,j,i) -= nud*fvz(i);
      if (eos.is_ideal) {
        flx3(m,IEN,k,j,i) -= 0.5*nud*((w0(m,IVX,k-1,j,i) + w0(m,IVX,k,j,i))*fvx(i) +
                                      (w0(m,IVY,k-1,j,i) + w0(m,IVY,k,j,i))*fvy(i) +
                                      (w0(m,IVZ,k-1,j,i) + w0(m,IVZ,k,j,i))*fvz(i));
      }
    });
  });

  return;
}
