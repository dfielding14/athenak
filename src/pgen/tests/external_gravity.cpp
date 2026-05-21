//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file external_gravity.cpp
//! \brief Example problem generators for fixed external gravitational potentials.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/external_gravity.hpp"
#include "srcterms/srcterms.hpp"

namespace {

void RequireHydroExternalGravity(MeshBlockPack *pmbp, const char *pgen_name) {
  if (pmbp->phydro == nullptr || pmbp->pmhd != nullptr) {
    std::cout << "### FATAL ERROR in " << pgen_name << std::endl
              << "This external-gravity example pgen currently initializes hydro only."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmbp->phydro->psrc == nullptr || !pmbp->phydro->psrc->external_gravity) {
    std::cout << "### FATAL ERROR in " << pgen_name << std::endl
              << "This pgen requires an enabled <external_gravity> block." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

KOKKOS_INLINE_FUNCTION
void StoreHydroState(DvceArray5D<Real> u0, const EOS_Data &eos, const int nhydro,
                     const int nscalars, const int m, const int k, const int j,
                     const int i, const Real rho, const Real pressure,
                     const Real v1, const Real v2, const Real v3) {
  u0(m, IDN, k, j, i) = rho;
  u0(m, IM1, k, j, i) = rho*v1;
  u0(m, IM2, k, j, i) = rho*v2;
  u0(m, IM3, k, j, i) = rho*v3;
  if (eos.is_ideal) {
    u0(m, IEN, k, j, i) = pressure/(eos.gamma - 1.0) +
                           0.5*rho*(SQR(v1) + SQR(v2) + SQR(v3));
  }
  for (int n=nhydro; n<(nhydro + nscalars); ++n) {
    u0(m, n, k, j, i) = 0.0;
  }
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::ExternalGravityHydrostatic()
//! \brief Isothermal hydrostatic atmosphere in an arbitrary configured potential.

void ProblemGenerator::ExternalGravityHydrostatic(ParameterInput *pin,
                                                 const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  RequireHydroExternalGravity(pmbp, "external_gravity_hydrostatic");

  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;

  const auto grav = pmbp->phydro->psrc->external_gravity_data;
  const auto eos = pmbp->phydro->peos->eos_data;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real pressure0 = pin->GetOrAddReal("problem", "pressure0", 1.0);
  const Real x1_ref = pin->GetOrAddReal("problem", "x1_ref", grav.x1_origin);
  const Real x2_ref = pin->GetOrAddReal("problem", "x2_ref", grav.x2_origin);
  const Real x3_ref = pin->GetOrAddReal("problem", "x3_ref", grav.x3_origin);
  const Real v1 = pin->GetOrAddReal("problem", "v1", 0.0);
  const Real v2 = pin->GetOrAddReal("problem", "v2", 0.0);
  const Real v3 = pin->GetOrAddReal("problem", "v3", 0.0);

  Real cs2 = 1.0;
  if (pin->DoesParameterExist("problem", "sound_speed")) {
    cs2 = SQR(pin->GetReal("problem", "sound_speed"));
  } else if (eos.is_ideal) {
    cs2 = pressure0/rho0;
  } else {
    cs2 = SQR(eos.iso_cs);
  }
  if (cs2 <= 0.0) {
    std::cout << "### FATAL ERROR in external_gravity_hydrostatic" << std::endl
              << "The hydrostatic sound speed squared must be positive." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const Real phi_ref = external_gravity::Potential(grav, x1_ref, x2_ref, x3_ref);

  par_for("pgen_extgrav_hydrostatic", DevExeSpace(), 0,(nmb-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1v = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    const Real phi = external_gravity::Potential(grav, x1v, x2v, x3v);
    const Real rho = fmax(rho0*exp(-(phi - phi_ref)/cs2), eos.dfloor);
    const Real pressure = fmax(cs2*rho, eos.pfloor);
    StoreHydroState(u0, eos, nhydro, nscalars, m, k, j, i, rho, pressure, v1, v2, v3);
  });
}

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::ExternalGravityDisk()
//! \brief Rotating hydro disk initialized from the configured external acceleration.

void ProblemGenerator::ExternalGravityDisk(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  RequireHydroExternalGravity(pmbp, "external_gravity_disk");

  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;

  const auto grav = pmbp->phydro->psrc->external_gravity_data;
  const auto eos = pmbp->phydro->peos->eos_data;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;

  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real pressure0 = pin->GetOrAddReal("problem", "pressure0", 0.01);
  const Real r_scale = pin->GetOrAddReal("problem", "r_scale", 1.0);
  const Real z_scale = pin->GetOrAddReal("problem", "z_scale", 0.1);
  const Real density_floor = pin->GetOrAddReal("problem", "density_floor", eos.dfloor);
  const bool rotate = pin->GetOrAddBoolean("problem", "rotate", true);

  par_for("pgen_extgrav_disk", DevExeSpace(), 0,(nmb-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1v = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    const Real x = x1v - grav.x1_origin;
    const Real y = x2v - grav.x2_origin;
    const Real z = x3v - grav.x3_origin;
    const Real r_cyl = sqrt(SQR(x) + SQR(y));
    const Real rho = fmax(rho0*exp(-r_cyl/r_scale)*exp(-fabs(z)/z_scale), density_floor);
    const Real pressure = fmax(pressure0*rho/rho0, eos.pfloor);

    Real v1 = 0.0, v2 = 0.0, v3 = 0.0;
    if (rotate && r_cyl > 0.0) {
      Real a1, a2, a3;
      external_gravity::Acceleration(grav, x1v, x2v, x3v, a1, a2, a3);
      const Real a_r = (a1*x + a2*y)/r_cyl;
      const Real vphi = sqrt(fmax(-r_cyl*a_r, 0.0));
      v1 = -vphi*y/r_cyl;
      v2 =  vphi*x/r_cyl;
    }

    StoreHydroState(u0, eos, nhydro, nscalars, m, k, j, i, rho, pressure, v1, v2, v3);
  });
}
