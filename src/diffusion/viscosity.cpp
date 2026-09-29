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
#include <cstdlib>
#include <limits>
#include <iostream>
#include <string> // string

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "viscosity.hpp"

namespace {

parabolic::ParabolicIntegratorMode ParseViscosityIntegrator(std::string block,
    ParameterInput *pin) {
  std::string integrator = pin->GetOrAddString(block, "viscosity_integrator", "explicit");
  if (integrator == "explicit") {
    return parabolic::ParabolicIntegratorMode::explicit_mode;
  } else if (integrator == "sts") {
    return parabolic::ParabolicIntegratorMode::sts;
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<" << block << ">/viscosity_integrator = '" << integrator
            << "' must be 'explicit' or 'sts'" << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real PowerLawNu(Real temp, Real nu_ref, Real temp_ref, Real exponent,
                Real floor, Real ceiling) {
  const Real ratio = fmax(temp/temp_ref, static_cast<Real>(1.0e-30));
  return fmin(fmax(nu_ref*pow(ratio, exponent), floor), ceiling);
}

} // namespace

//----------------------------------------------------------------------------------------
// ctor:
// Note first argument passes string ("hydro" or "mhd") denoting in which class this
// object is being constructed, and therefore which <block> in the input file from which
// the parameters are read.

Viscosity::Viscosity(std::string block, MeshBlockPack *pp,
                     ParameterInput *pin) :
  pmy_pack(pp) {
  // Read coefficient of isotropic kinematic shear viscosity (must be present)
  nu_iso = pin->GetReal(block,"viscosity");
  tdep_nu = pin->GetOrAddBoolean(block, "tdep_viscosity", false);
  nu_tref = pin->GetOrAddReal(block, "viscosity_tref", 1.0);
  nu_exponent = pin->GetOrAddReal(block, "viscosity_exponent", 0.0);
  nu_floor = pin->GetOrAddReal(block, "viscosity_floor", 0.0);
  nu_ceiling = pin->GetOrAddReal(
      block, "viscosity_ceiling", static_cast<Real>(std::numeric_limits<float>::max()));
  if (nu_iso < 0.0 || nu_tref <= 0.0 || nu_floor < 0.0 || nu_ceiling < nu_floor) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Invalid viscosity coefficient, reference temperature, or power-law "
              << "saturation bounds in <" << block << ">" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mode = ParseViscosityIntegrator(block, pin);
  dtnew = std::numeric_limits<float>::max();
}

//----------------------------------------------------------------------------------------
// Viscosity destructor

Viscosity::~Viscosity() {
}

//----------------------------------------------------------------------------------------
//! \fn void Viscosity::NewTimeStep()
//! \brief Compute new time step for viscosity.

void Viscosity::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  dtnew = std::numeric_limits<float>::max();
  auto size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  if (!tdep_nu) {
    if (nu_iso <= 0.0) return;
    for (int m=0; m<(pmy_pack->nmb_thispack); ++m) {
      dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx1)/nu_iso);
      if (pmy_pack->pmesh->multi_d) {
        dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx2)/nu_iso);
      }
      if (pmy_pack->pmesh->three_d) {
        dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx3)/nu_iso);
      }
    }
    return;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool use_e = eos_data.use_e;
  const Real gm1 = eos_data.gamma - 1.0;
  const Real nu_ref = nu_iso;
  const Real temp_ref = nu_tref;
  const Real exponent = nu_exponent;
  const Real floor = nu_floor;
  const Real ceiling = nu_ceiling;

  Kokkos::parallel_reduce("visc_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &min_dt) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k-ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k-ks)*nji - (j-js)*nx1 + is;
    const Real temp = use_e ? gm1*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                            : w0(m,ITM,k,j,i);
    const Real nu = PowerLawNu(temp, nu_ref, temp_ref, exponent, floor, ceiling);
    if (nu > 0.0) {
      min_dt = fmin(min_dt, fac*SQR(size.d_view(m).dx1)/nu);
      if (multi_d) min_dt = fmin(min_dt, fac*SQR(size.d_view(m).dx2)/nu);
      if (three_d) min_dt = fmin(min_dt, fac*SQR(size.d_view(m).dx3)/nu);
    }
  }, Kokkos::Min<Real>(dtnew));
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsoViscousFlux
//  \brief Adds viscous fluxes to face-centered fluxes of conserved variables

void Viscosity::IsotropicViscousFlux(const DvceArray5D<Real> &w0, const Real nu_iso,
  const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  const bool variable_nu = tdep_nu;
  const bool use_e = eos.use_e;
  const Real gm1 = eos.gamma - 1.0;
  const Real nu_ref = nu_iso;
  const Real temp_ref = nu_tref;
  const Real exponent = nu_exponent;
  const Real floor = nu_floor;
  const Real ceiling = nu_ceiling;

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
      Real nu = nu_ref;
      if (variable_nu) {
        const Real tl = use_e ? gm1*w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)
                              : w0(m,ITM,k,j,i-1);
        const Real tr = use_e ? gm1*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                              : w0(m,ITM,k,j,i);
        nu = 0.5*(PowerLawNu(tl, nu_ref, temp_ref, exponent, floor, ceiling) +
                  PowerLawNu(tr, nu_ref, temp_ref, exponent, floor, ceiling));
      }
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
      Real nu = nu_ref;
      if (variable_nu) {
        const Real tl = use_e ? gm1*w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)
                              : w0(m,ITM,k,j-1,i);
        const Real tr = use_e ? gm1*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                              : w0(m,ITM,k,j,i);
        nu = 0.5*(PowerLawNu(tl, nu_ref, temp_ref, exponent, floor, ceiling) +
                  PowerLawNu(tr, nu_ref, temp_ref, exponent, floor, ceiling));
      }
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
      Real nu = nu_ref;
      if (variable_nu) {
        const Real tl = use_e ? gm1*w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)
                              : w0(m,ITM,k-1,j,i);
        const Real tr = use_e ? gm1*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                              : w0(m,ITM,k,j,i);
        nu = 0.5*(PowerLawNu(tl, nu_ref, temp_ref, exponent, floor, ceiling) +
                  PowerLawNu(tr, nu_ref, temp_ref, exponent, floor, ceiling));
      }
      Real nud = 0.5*nu*(w0(m,IDN,k,j,i) + w0(m,IDN,k-1,j,i));
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
