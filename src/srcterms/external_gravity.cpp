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

  grav.x1_origin = GetRealAny(pin, "x1_origin", 0.0, {"x0", "x_origin"});
  grav.x2_origin = GetRealAny(pin, "x2_origin", 0.0, {"y0", "y_origin"});
  grav.x3_origin = GetRealAny(pin, "x3_origin", 0.0, {"z0", "z_origin"});

  grav.g1 = GetRealAny(pin, "g1", 0.0, {"gx", "accel_x1"});
  grav.g2 = GetRealAny(pin, "g2", 0.0, {"gy", "accel_x2"});
  grav.g3 = GetRealAny(pin, "g3", 0.0, {"gz", "accel_x3"});

  const Real omega_default = GetRealAny(pin, "omega", 0.0);
  grav.omega = omega_default;
  grav.omega1 = GetRealAny(pin, "omega1", omega_default, {"omega_x1"});
  grav.omega2 = GetRealAny(pin, "omega2", omega_default, {"omega_x2"});
  grav.omega3 = GetRealAny(pin, "omega3", omega_default, {"omega_x3"});

  grav.mass = GetRealAny(pin, "mass", 1.0, {"M"});
  grav.softening = GetRealAny(pin, "softening", 0.0, {"r_soft", "eps"});

  grav.r_scale = GetRealAny(pin, "r_scale", 1.0, {"rs", "scale_radius"});
  grav.rho_scale = GetRealAny(pin, "rho_scale", 1.0, {"rho_s"});

  grav.v0 = GetRealAny(pin, "v0", 1.0, {"vc", "v_c"});
  grav.core_radius = GetRealAny(pin, "core_radius", 1.0, {"r_core", "rc"});
  grav.flattening_q = GetRealAny(pin, "q", 1.0, {"flattening", "q_phi"});

  grav.disk_mass = GetRealAny(pin, "disk_mass", grav.mass, {"mass_gal"});
  grav.disk_a = GetRealAny(pin, "disk_a", 1.0, {"scale_gal", "a_gal"});
  grav.disk_b = GetRealAny(pin, "disk_b", 0.1, {"z_gal", "b_gal", "scale_height"});

  grav.rho_mean = GetRealAny(pin, "rho_mean", 0.0);
  const Real r200 = GetRealAny(pin, "r_200", grav.r_scale);
  const Real outer_factor = GetRealAny(pin, "outer_radius_factor", 5.0);
  grav.r_outer = GetRealAny(pin, "r_outer", outer_factor*r200);

  RequireNonNegative("softening", grav.softening);
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
    const Real x1l = LeftEdgeX(i - is, nx1, x1min, x1max);
    const Real x1r = LeftEdgeX(i + 1 - is, nx1, x1min, x1max);

    const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);
    const Real x2l = LeftEdgeX(j - js, nx2, x2min, x2max);
    const Real x2r = LeftEdgeX(j + 1 - js, nx2, x2min, x2max);

    const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);
    const Real x3l = LeftEdgeX(k - ks, nx3, x3min, x3max);
    const Real x3r = LeftEdgeX(k + 1 - ks, nx3, x3min, x3max);

    const Real phi1l = external_gravity::Potential(grav, x1l, x2v, x3v);
    const Real phi1r = external_gravity::Potential(grav, x1r, x2v, x3v);
    const Real phi2l = external_gravity::Potential(grav, x1v, x2l, x3v);
    const Real phi2r = external_gravity::Potential(grav, x1v, x2r, x3v);
    const Real phi3l = external_gravity::Potential(grav, x1v, x2v, x3l);
    const Real phi3r = external_gravity::Potential(grav, x1v, x2v, x3r);

    const Real tiny = 1.0e-20;
    const Real a1 = -(phi1r - phi1l)/fmax(x1r - x1l, tiny);
    const Real a2 = -(phi2r - phi2l)/fmax(x2r - x2l, tiny);
    const Real a3 = -(phi3r - phi3l)/fmax(x3r - x3l, tiny);

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
