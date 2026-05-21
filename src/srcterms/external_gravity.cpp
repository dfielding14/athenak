//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file external_gravity.cpp
//! \brief Input parsing and source-term application for external gravity.

#include "external_gravity.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "srcterms.hpp"
#include "eos/eos.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "units/units.hpp"

namespace {

constexpr const char *kBlock = "external_gravity";
constexpr const char *kBoundaryBlock = "external_gravity_boundary";

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

Real GetRealAny(ParameterInput *pin, const std::string &primary, const Real def,
                std::initializer_list<const char *> aliases = {}) {
  if (pin->DoesParameterExist(kBlock, primary)) return pin->GetReal(kBlock, primary);
  for (const char *alias : aliases) {
    if (pin->DoesParameterExist(kBlock, alias)) return pin->GetReal(kBlock, alias);
  }
  return pin->GetOrAddReal(kBlock, primary, def);
}

void FatalExternalGravityError(const std::string &msg) {
  std::cout << "### FATAL ERROR in external gravity: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void RequirePositive(const std::string &name, const Real val) {
  if (val <= 0.0) {
    FatalExternalGravityError("<external_gravity>/" + name + " must be positive");
  }
}

void RequireNonNegative(const std::string &name, const Real val) {
  if (val < 0.0) {
    FatalExternalGravityError("<external_gravity>/" + name + " must be non-negative");
  }
}

int ParseModel(const std::string &model) {
  const std::string m = Lower(model);
  if (m == "uniform" || m == "constant" || m == "constant_acceleration" ||
      m == "cartesian_linear" || m == "linear") {
    return external_gravity::Model::uniform;
  }
  if (m == "cartesian_harmonic" || m == "harmonic_cartesian") {
    return external_gravity::Model::cartesian_harmonic;
  }
  if (m == "spherical_harmonic" || m == "harmonic" || m == "solid_body") {
    return external_gravity::Model::spherical_harmonic;
  }
  if (m == "point_mass" || m == "kepler" || m == "keplerian") {
    return external_gravity::Model::point_mass;
  }
  if (m == "plummer") {
    return external_gravity::Model::plummer;
  }
  if (m == "nfw") {
    return external_gravity::Model::nfw;
  }
  if (m == "logarithmic" || m == "logarithmic_halo" || m == "log_halo") {
    return external_gravity::Model::logarithmic_halo;
  }
  if (m == "miyamoto_nagai" || m == "miyamoto-nagai" || m == "mn_disk") {
    return external_gravity::Model::miyamoto_nagai;
  }
  if (m == "outer" || m == "outer_cgm" || m == "uniform_background") {
    return external_gravity::Model::outer_cgm;
  }
  if (m == "gotham" || m == "galaxy" || m == "nfw_miyamoto_nagai_outer") {
    return external_gravity::Model::gotham;
  }
  if (m == "user" || m == "custom") {
    return external_gravity::Model::user;
  }

  FatalExternalGravityError("unknown <external_gravity>/model = '" + model + "'");
  return external_gravity::Model::none;
}

int ParseBoundaryVelocity(const std::string &mode) {
  const std::string m = Lower(mode);
  if (m == "copy" || m == "copy_velocity") {
    return external_gravity::BoundaryVelocityMode::copy_velocity;
  }
  if (m == "zero_normal" || m == "zero_normal_velocity") {
    return external_gravity::BoundaryVelocityMode::zero_normal_velocity;
  }
  if (m == "no_inflow" || m == "diode") {
    return external_gravity::BoundaryVelocityMode::no_inflow_velocity;
  }
  if (m == "reflect_normal" || m == "reflect") {
    return external_gravity::BoundaryVelocityMode::reflect_normal_velocity;
  }

  FatalExternalGravityError("unknown <external_gravity_boundary>/velocity = '" +
                            mode + "'");
  return external_gravity::BoundaryVelocityMode::no_inflow_velocity;
}

} // namespace

namespace external_gravity {

ExternalGravityData ParseInput(MeshBlockPack *pp, ParameterInput *pin) {
  ExternalGravityData grav;
  if (pin == nullptr || !pin->DoesBlockExist(kBlock)) {
    return grav;
  }

  const bool enabled = pin->GetOrAddBoolean(kBlock, "enabled", true);
  if (!enabled) {
    return grav;
  }

  grav.model = ParseModel(pin->GetOrAddString(kBlock, "model", "uniform"));
  const Real default_g = (pp != nullptr && pp->punit != nullptr) ?
                         pp->punit->grav_constant() : 1.0;
  grav.grav_constant = GetRealAny(pin, "G", default_g, {"grav_constant"});

  if (pp != nullptr && pp->pmesh != nullptr) {
    grav.mesh_x1min = pp->pmesh->mesh_size.x1min;
    grav.mesh_x1max = pp->pmesh->mesh_size.x1max;
    grav.mesh_x2min = pp->pmesh->mesh_size.x2min;
    grav.mesh_x2max = pp->pmesh->mesh_size.x2max;
    grav.mesh_x3min = pp->pmesh->mesh_size.x3min;
    grav.mesh_x3max = pp->pmesh->mesh_size.x3max;
  }

  grav.x1_origin = GetRealAny(pin, "x1_origin", 0.0, {"x0", "x_origin"});
  grav.x2_origin = GetRealAny(pin, "x2_origin", 0.0, {"y0", "y_origin"});
  grav.x3_origin = GetRealAny(pin, "x3_origin", 0.0, {"z0", "z_origin"});

  switch (grav.model) {
    case Model::uniform:
      grav.g1 = GetRealAny(pin, "g1", 0.0, {"gx", "accel_x1"});
      grav.g2 = GetRealAny(pin, "g2", 0.0, {"gy", "accel_x2"});
      grav.g3 = GetRealAny(pin, "g3", 0.0, {"gz", "accel_x3"});
      break;

    case Model::cartesian_harmonic: {
      const Real omega_default = GetRealAny(pin, "omega", 0.0);
      grav.omega1 = GetRealAny(pin, "omega1", omega_default, {"omega_x1"});
      grav.omega2 = GetRealAny(pin, "omega2", omega_default, {"omega_x2"});
      grav.omega3 = GetRealAny(pin, "omega3", omega_default, {"omega_x3"});
      break;
    }

    case Model::spherical_harmonic:
      grav.omega = GetRealAny(pin, "omega", 0.0);
      break;

    case Model::point_mass:
    case Model::plummer:
      grav.mass = GetRealAny(pin, "mass", 1.0, {"M"});
      grav.softening = GetRealAny(pin, "softening", 0.0, {"r_soft", "eps"});
      break;

    case Model::nfw:
      grav.r_scale = GetRealAny(pin, "r_scale", 1.0, {"rs", "scale_radius"});
      grav.rho_scale = GetRealAny(pin, "rho_scale", 1.0, {"rho_s"});
      break;

    case Model::logarithmic_halo:
      grav.v0 = GetRealAny(pin, "v0", 1.0, {"vc", "v_c"});
      grav.core_radius = GetRealAny(pin, "core_radius", 1.0, {"r_core", "rc"});
      grav.flattening_q = GetRealAny(pin, "q", 1.0, {"flattening", "q_phi"});
      break;

    case Model::miyamoto_nagai:
      grav.disk_mass = GetRealAny(pin, "disk_mass", 1.0, {"mass_gal"});
      grav.disk_a = GetRealAny(pin, "disk_a", 1.0, {"scale_gal", "a_gal"});
      grav.disk_b = GetRealAny(pin, "disk_b", 0.1, {"z_gal", "b_gal", "scale_height"});
      break;

    case Model::outer_cgm: {
      grav.rho_mean = GetRealAny(pin, "rho_mean", 0.0);
      const Real r200 = GetRealAny(pin, "r_200", 1.0);
      const Real outer_factor = GetRealAny(pin, "outer_radius_factor", 5.0);
      grav.r_outer = GetRealAny(pin, "r_outer", outer_factor*r200);
      break;
    }

    case Model::gotham: {
      grav.r_scale = GetRealAny(pin, "r_scale", 1.0, {"rs", "scale_radius"});
      grav.rho_scale = GetRealAny(pin, "rho_scale", 1.0, {"rho_s"});
      grav.disk_mass = GetRealAny(pin, "disk_mass", 1.0, {"mass_gal"});
      grav.disk_a = GetRealAny(pin, "disk_a", 1.0, {"scale_gal", "a_gal"});
      grav.disk_b = GetRealAny(pin, "disk_b", 0.1, {"z_gal", "b_gal", "scale_height"});
      grav.rho_mean = GetRealAny(pin, "rho_mean", 0.0);
      const Real r200 = GetRealAny(pin, "r_200", grav.r_scale);
      const Real outer_factor = GetRealAny(pin, "outer_radius_factor", 5.0);
      grav.r_outer = GetRealAny(pin, "r_outer", outer_factor*r200);
      break;
    }

    case Model::user:
      grav.user_fd_step = GetRealAny(pin, "user_fd_step", 1.0e-5, {"fd_step"});
      break;

    default:
      break;
  }

  grav.taper_gravity = pin->GetOrAddBoolean(kBlock, "taper_gravity", false) ? 1 : 0;
  if (grav.taper_gravity != 0) {
    const Real taper_width = GetRealAny(pin, "taper_width", 0.0);
    grav.taper_width_x1_inner = GetRealAny(pin, "taper_width_x1_inner", taper_width);
    grav.taper_width_x1_outer = GetRealAny(pin, "taper_width_x1_outer", taper_width);
    grav.taper_width_x2_inner = GetRealAny(pin, "taper_width_x2_inner", taper_width);
    grav.taper_width_x2_outer = GetRealAny(pin, "taper_width_x2_outer", taper_width);
    grav.taper_width_x3_inner = GetRealAny(pin, "taper_width_x3_inner", taper_width);
    grav.taper_width_x3_outer = GetRealAny(pin, "taper_width_x3_outer", taper_width);
  }

  if (pin->DoesBlockExist(kBoundaryBlock)) {
    grav.boundary_velocity = ParseBoundaryVelocity(
        pin->GetOrAddString(kBoundaryBlock, "velocity", "no_inflow"));
    grav.boundary_sound_speed =
        pin->GetOrAddReal(kBoundaryBlock, "sound_speed", -1.0);
    grav.boundary_density_floor =
        pin->GetOrAddReal(kBoundaryBlock, "density_floor", -1.0);
    grav.boundary_pressure_floor =
        pin->GetOrAddReal(kBoundaryBlock, "pressure_floor", -1.0);
    grav.boundary_max_exponent =
        pin->GetOrAddReal(kBoundaryBlock, "max_exponent", 60.0);
  }

  RequireNonNegative("softening", grav.softening);
  RequireNonNegative("user_fd_step", grav.user_fd_step);
  if (grav.taper_gravity != 0) {
    RequireNonNegative("taper_width_x1_inner", grav.taper_width_x1_inner);
    RequireNonNegative("taper_width_x1_outer", grav.taper_width_x1_outer);
    RequireNonNegative("taper_width_x2_inner", grav.taper_width_x2_inner);
    RequireNonNegative("taper_width_x2_outer", grav.taper_width_x2_outer);
    RequireNonNegative("taper_width_x3_inner", grav.taper_width_x3_inner);
    RequireNonNegative("taper_width_x3_outer", grav.taper_width_x3_outer);
  }
  if (grav.model == Model::point_mass || grav.model == Model::plummer) {
    RequirePositive("mass", grav.mass);
  }
  if (grav.model == Model::nfw || grav.model == Model::gotham) {
    RequirePositive("r_scale", grav.r_scale);
    RequirePositive("rho_scale", grav.rho_scale);
  }
  if (grav.model == Model::logarithmic_halo) {
    RequirePositive("v0", grav.v0);
    RequirePositive("core_radius", grav.core_radius);
    RequirePositive("q", grav.flattening_q);
  }
  if (grav.model == Model::miyamoto_nagai || grav.model == Model::gotham) {
    RequirePositive("disk_mass", grav.disk_mass);
    RequirePositive("disk_a", grav.disk_a);
    RequirePositive("disk_b", grav.disk_b);
  }
  if (grav.model == Model::outer_cgm || grav.model == Model::gotham) {
    RequirePositive("r_outer", grav.r_outer);
  }
  if (grav.boundary_sound_speed == 0.0) {
    FatalExternalGravityError("<external_gravity_boundary>/sound_speed must be positive "
                              "when specified");
  }
  if (grav.boundary_density_floor == 0.0 || grav.boundary_pressure_floor == 0.0) {
    FatalExternalGravityError("<external_gravity_boundary> floors must be positive "
                              "when specified");
  }
  if (grav.boundary_max_exponent <= 0.0) {
    FatalExternalGravityError("<external_gravity_boundary>/max_exponent must be "
                              "positive");
  }

  return grav;
}

} // namespace external_gravity

//----------------------------------------------------------------------------------------
//! \fn SourceTerms::ExternalGravity()
//! \brief Add momentum and ideal-gas energy source terms from an external potential.

void SourceTerms::ExternalGravity(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                                  const Real bdt, DvceArray5D<Real> &u0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;

  external_gravity::ExternalGravityData grav = external_gravity_data;

  par_for("external_gravity", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min;
    const Real x1max = size.d_view(m).x1max;
    const Real x2min = size.d_view(m).x2min;
    const Real x2max = size.d_view(m).x2max;
    const Real x3min = size.d_view(m).x3min;
    const Real x3max = size.d_view(m).x3max;

    const Real x1v = CellCenterX(i - is, nx1, x1min, x1max);
    const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);
    const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);

    Real a1, a2, a3;
    external_gravity::Acceleration(grav, x1v, x2v, x3v, a1, a2, a3);

    const Real density = w0(m, IDN, k, j, i);
    const Real src1 = bdt*density*a1;
    const Real src2 = bdt*density*a2;
    const Real src3 = bdt*density*a3;

    u0(m, IM1, k, j, i) += src1;
    u0(m, IM2, k, j, i) += src2;
    u0(m, IM3, k, j, i) += src3;
    if (eos_data.is_ideal) {
      u0(m, IEN, k, j, i) += src1*w0(m, IVX, k, j, i)
                            + src2*w0(m, IVY, k, j, i)
                            + src3*w0(m, IVZ, k, j, i);
    }
  });
}
