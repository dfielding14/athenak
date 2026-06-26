//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q019_nonlinear_bell_saturation_engineering.cpp
//! \brief Compact finite-rigidity nonlinear Bell engineering initial conditions.

#include <algorithm>
#include <cmath>

#if !defined(Q019_NONLINEAR_BELL_SATURATION_ENGINEERING_HOST_CONTRACT)
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"
#endif

#if defined(Q019_NONLINEAR_BELL_SATURATION_ENGINEERING_HOST_CONTRACT) && \
    !defined(Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT)
#define Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT 1
#define Q019_ENGINEERING_UNDEF_Q019_HOST_CONTRACT 1
#endif
#include "q019_physics_first_nonlinear_bell_successor_v2.hpp"
#if defined(Q019_ENGINEERING_UNDEF_Q019_HOST_CONTRACT)
#undef Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT
#undef Q019_ENGINEERING_UNDEF_Q019_HOST_CONTRACT
#endif

namespace q019_nonlinear_bell_saturation_engineering {

using q019_physics_first_nonlinear_bell_successor_v2::AxisAlignedEigenmodeVelocityAt;
using q019_physics_first_nonlinear_bell_successor_v2::AxisAlignedVectorPotentialAt;
using q019_physics_first_nonlinear_bell_successor_v2::SeedParameters;
using q019_physics_first_nonlinear_bell_successor_v2::SeedTopology;
using q019_physics_first_nonlinear_bell_successor_v2::Vector3;

inline double RootCellVolume(const double x1_extent, const int nx1,
                             const double x2_extent, const int nx2,
                             const double x3_extent, const int nx3) {
  return (x1_extent/static_cast<double>(nx1))*
      (x2_extent/static_cast<double>(nx2))*
      (x3_extent/static_cast<double>(nx3));
}

inline double DepositedJOverC(const double ppc, const double qscale,
                              const double charge, const double stream_speed,
                              const double cell_volume) {
  return ppc*qscale*charge*stream_speed/cell_volume;
}

inline bool NearlyEqual(const double measured, const double expected) {
  const double scale = std::max(1.0, std::abs(expected));
  return std::isfinite(measured) && std::isfinite(expected) &&
      std::abs(measured - expected) <= 1.0e-12*scale;
}

}  // namespace q019_nonlinear_bell_saturation_engineering

#if !defined(Q019_NONLINEAR_BELL_SATURATION_ENGINEERING_HOST_CONTRACT)
namespace {

[[noreturn]] void Q019EngineeringFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q019EngineeringRequireString(ParameterInput *pin, const std::string &block,
                                  const std::string &name,
                                  const std::string &expected) {
  if (pin->GetString(block, name) != expected) {
    Q019EngineeringFatal("<" + block + ">/" + name + " must be " + expected);
  }
}

void Q019EngineeringRequireBoolean(ParameterInput *pin, const std::string &block,
                                   const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q019EngineeringFatal("<" + block + ">/" + name +
                         " does not match the coupled paper-TSC contract");
  }
}

void Q019EngineeringRequirePositive(const std::string &name, const Real value) {
  if (!std::isfinite(value) || value <= 0.0) {
    Q019EngineeringFatal(name + " must be finite and positive");
  }
}

void Q019EngineeringRequireClose(const std::string &name, const Real measured,
                                 const Real expected) {
  if (!q019_nonlinear_bell_saturation_engineering::NearlyEqual(
          measured, expected)) {
    Q019EngineeringFatal(name + " does not match the nonlinear Bell normalization");
  }
}

}  // namespace

void ProblemGenerator::Q019NonlinearBellSaturationEngineering(
    ParameterInput *pin, const bool restart) {
  using q019_nonlinear_bell_saturation_engineering::AxisAlignedEigenmodeVelocityAt;
  using q019_nonlinear_bell_saturation_engineering::AxisAlignedVectorPotentialAt;
  using q019_nonlinear_bell_saturation_engineering::DepositedJOverC;
  using q019_nonlinear_bell_saturation_engineering::RootCellVolume;
  using q019_nonlinear_bell_saturation_engineering::SeedParameters;
  using q019_nonlinear_bell_saturation_engineering::SeedTopology;
  using q019_nonlinear_bell_saturation_engineering::Vector3;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q019EngineeringFatal(
        "q019_nonlinear_bell_saturation_engineering requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q019EngineeringFatal(
        "q019_nonlinear_bell_saturation_engineering requires an ideal EOS");
  }
  if (pmy_mesh_->one_d || pmy_mesh_->two_d || pmy_mesh_->multilevel ||
      pmy_mesh_->mb_indcs.nx1 <= 1 || pmy_mesh_->mb_indcs.nx2 <= 1 ||
      pmy_mesh_->mb_indcs.nx3 <= 1) {
    Q019EngineeringFatal(
        "q019_nonlinear_bell_saturation_engineering requires a uniform 3D mesh");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q019EngineeringFatal(
        "q019_nonlinear_bell_saturation_engineering requires periodic boundaries");
  }

  const std::string block = "q019_nonlinear_bell_saturation_engineering";
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real epsilon = pin->GetReal(block, "epsilon");
  const Real omega = pin->GetReal(block, "omega");
  const Real eigenmode_amplitude = pin->GetReal(block, "eigenmode_amplitude");
  const Real broadband_amplitude = pin->GetReal(block, "broadband_amplitude");
  const int field_seed = pin->GetInteger(block, "field_seed");
  Q019EngineeringRequirePositive("rho", rho);
  Q019EngineeringRequirePositive("pressure", pressure);
  Q019EngineeringRequirePositive("b_g", b_g);
  Q019EngineeringRequirePositive("u_a", u_a);
  Q019EngineeringRequirePositive("wavelength", wavelength);
  Q019EngineeringRequirePositive("k0", k0);
  Q019EngineeringRequirePositive("epsilon", epsilon);
  Q019EngineeringRequirePositive("omega", omega);
  Q019EngineeringRequirePositive("eigenmode_amplitude", eigenmode_amplitude);
  Q019EngineeringRequirePositive("broadband_amplitude", broadband_amplitude);
  if (epsilon >= 1.0 || field_seed <= 0) {
    Q019EngineeringFatal("epsilon must be less than one and field_seed must be positive");
  }
  Q019EngineeringRequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q019EngineeringRequireClose("k0", k0, 2.0*M_PI/wavelength);

  Q019EngineeringRequireString(pin, "mhd", "eos", "ideal");
  Q019EngineeringRequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q019EngineeringRequireString(pin, "particles", "pusher", "boris_tsc");
  Q019EngineeringRequireString(pin, "particles", "cr_distribution", "center");
  Q019EngineeringRequireString(pin, "particles", "pic_physical_mode",
                               "paper_mhd_pic_vl2_tsc");
  Q019EngineeringRequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q019EngineeringRequireBoolean(pin, "particles", "deposit_moments", true);
  Q019EngineeringRequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q019EngineeringRequireBoolean(pin, "particles",
                                "couple_moments_momentum_to_mhd", true);
  Q019EngineeringRequireBoolean(pin, "particles",
                                "couple_moments_energy_to_mhd", true);
  if (pin->GetInteger("particles", "deposit_order") != 2 ||
      pin->GetInteger("particles", "nspecies") != 1) {
    Q019EngineeringFatal("coupled paper TSC requires deposit_order=2 and nspecies=1");
  }

  const Real ppc = pin->GetReal("particles", "ppc");
  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real cr_vx = pin->GetReal("particles", "cr_vx0");
  const Real cr_vy = pin->GetReal("particles", "cr_vy0");
  const Real cr_vz = pin->GetReal("particles", "cr_vz0");
  const Real mass = pin->GetReal("species0", "mass");
  const Real charge = pin->GetReal("species0", "charge");
  const Real species_vx = pin->GetReal("species0", "vx0");
  const Real species_vy = pin->GetReal("species0", "vy0");
  const Real species_vz = pin->GetReal("species0", "vz0");
  for (const auto value : {ppc, qscale, mass, charge}) {
    Q019EngineeringRequirePositive("particle normalization value", value);
  }
  if (std::floor(ppc) != ppc) {
    Q019EngineeringFatal("<particles>/ppc must be a positive integer");
  }
  const Real vstream = u_a/epsilon;
  Q019EngineeringRequireClose("cr_vx0", cr_vx, vstream);
  Q019EngineeringRequireClose("cr_vy0", cr_vy, 0.0);
  Q019EngineeringRequireClose("cr_vz0", cr_vz, 0.0);
  Q019EngineeringRequireClose("species0/vx0", species_vx, vstream);
  Q019EngineeringRequireClose("species0/vy0", species_vy, 0.0);
  Q019EngineeringRequireClose("species0/vz0", species_vz, 0.0);
  Q019EngineeringRequireClose("epsilon", epsilon, u_a/vstream);
  Q019EngineeringRequireClose("species0 charge-to-mass", charge/mass, omega/b_g);

  const Real x1_extent = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const Real x2_extent = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  const Real x3_extent = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  const Real cell_volume = RootCellVolume(
      x1_extent, pmy_mesh_->mesh_indcs.nx1,
      x2_extent, pmy_mesh_->mesh_indcs.nx2,
      x3_extent, pmy_mesh_->mesh_indcs.nx3);
  Q019EngineeringRequirePositive("root-cell volume", cell_volume);
  Q019EngineeringRequireClose(
      "root-cell volume", cell_volume,
      pmy_mesh_->mesh_size.dx1*pmy_mesh_->mesh_size.dx2*pmy_mesh_->mesh_size.dx3);
  Q019EngineeringRequireClose(
      "PPC*qscale*q*vstream/Vcell",
      DepositedJOverC(ppc, qscale, charge, vstream, cell_volume), 2.0*b_g*k0);
  for (const auto extent : {x1_extent, x2_extent, x3_extent}) {
    Q019EngineeringRequireClose("domain extent/wavelength", extent/wavelength,
                                std::round(extent/wavelength));
  }

  if (restart) return;

  const SeedParameters parameters = {
    3, field_seed, b_g, k0, epsilon, rho, eigenmode_amplitude,
    broadband_amplitude, 0.0, x1_extent, x2_extent, x3_extent,
    SeedTopology::shared_spectrum_only
  };
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;

  DvceArray4D<Real> a1("q019_engineering_a1", nmb, indcs.nx3 + 2*indcs.ng,
                       indcs.nx2 + 2*indcs.ng, indcs.nx1 + 2*indcs.ng);
  DvceArray4D<Real> a2("q019_engineering_a2", nmb, indcs.nx3 + 2*indcs.ng,
                       indcs.nx2 + 2*indcs.ng, indcs.nx1 + 2*indcs.ng);
  DvceArray4D<Real> a3("q019_engineering_a3", nmb, indcs.nx3 + 2*indcs.ng,
                       indcs.nx2 + 2*indcs.ng, indcs.nx1 + 2*indcs.ng);

  par_for("pgen_q019_engineering_vector_potential", DevExeSpace(),
          0, nmb - 1, ks, ke + 1, js, je + 1, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1v = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i - is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j - js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k - ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Vector3 av1 = AxisAlignedVectorPotentialAt(parameters, {x1v, x2f, x3f});
    const Vector3 av2 = AxisAlignedVectorPotentialAt(parameters, {x1f, x2v, x3f});
    const Vector3 av3 = AxisAlignedVectorPotentialAt(parameters, {x1f, x2f, x3v});
    a1(m, k, j, i) = av1.x1;
    a2(m, k, j, i) = av2.x2;
    a3(m, k, j, i) = av3.x3;
  });

  par_for("pgen_q019_engineering_ct_field", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m, k, j, i) = (a3(m, k, j + 1, i) - a3(m, k, j, i))/dx2 -
                          (a2(m, k + 1, j, i) - a2(m, k, j, i))/dx3;
    b0.x2f(m, k, j, i) = (a1(m, k + 1, j, i) - a1(m, k, j, i))/dx3 -
                          (a3(m, k, j, i + 1) - a3(m, k, j, i))/dx1;
    b0.x3f(m, k, j, i) = (a2(m, k, j, i + 1) - a2(m, k, j, i))/dx1 -
                          (a1(m, k, j + 1, i) - a1(m, k, j, i))/dx2;
    if (i == ie) {
      b0.x1f(m, k, j, i + 1) =
          (a3(m, k, j + 1, i + 1) - a3(m, k, j, i + 1))/dx2 -
          (a2(m, k + 1, j, i + 1) - a2(m, k, j, i + 1))/dx3;
    }
    if (j == je) {
      b0.x2f(m, k, j + 1, i) =
          (a1(m, k + 1, j + 1, i) - a1(m, k, j + 1, i))/dx3 -
          (a3(m, k, j + 1, i + 1) - a3(m, k, j + 1, i))/dx1;
    }
    if (k == ke) {
      b0.x3f(m, k + 1, j, i) =
          (a2(m, k + 1, j, i + 1) - a2(m, k + 1, j, i))/dx1 -
          (a1(m, k + 1, j + 1, i) - a1(m, k + 1, j, i))/dx2;
    }
  });

  par_for("pgen_q019_engineering_primitives", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const Vector3 velocity = AxisAlignedEigenmodeVelocityAt(parameters, {x1, x2, x3});
    w0(m, IDN, k, j, i) = rho;
    w0(m, IVX, k, j, i) = velocity.x1;
    w0(m, IVY, k, j, i) = velocity.x2;
    w0(m, IVZ, k, j, i) = velocity.x3;
    w0(m, IEN, k, j, i) = pressure;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
#endif
