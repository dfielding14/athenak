//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cloud_crushing.cpp
//! \brief Cold-cloud crushing test with ISM cooling and user-selectable inflow boundary.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "parameter_input.hpp"
#include "pgen.hpp"
#include "srcterms/frame_tracker.hpp"
#include "srcterms/ismcooling.hpp"
#include "units/units.hpp"

namespace {

enum CloudCrushingBoundaryMode {
  kBoundarySedovTaylor = 0,
  kBoundaryConstant = 1
};

struct PhaseState {
  Real temp_cgs = 0.0;
  Real number_density_cgs = 0.0;
  Real density_code = 0.0;
  Real pressure_code = 0.0;
};

struct CloudCrushingData {
  Real pressure_over_k = 0.0;
  Real hrate = 0.0;
  Real n_unit = 1.0;
  Real pressure_unit = 1.0;
  Real velocity_unit = 1.0;
  Real time_unit = 1.0;
  Real mu = 1.0;
  Real gm1 = 2.0/3.0;

  PhaseState cold;
  PhaseState warm;

  Real cloud_radius = 1.0;
  Real cloud_smoothing_width = 0.0;
  Real cloud_center_x1 = 0.0;
  Real cloud_center_x2 = 0.0;
  Real cloud_center_x3 = 0.0;
  int cloud_tracer_scalar_index = -1;

  int boundary_mode = kBoundarySedovTaylor;
  Real constant_density_code = 0.0;
  Real constant_pressure_code = 0.0;
  Real constant_vx_code = 0.0;
  Real constant_vy_code = 0.0;
  Real constant_vz_code = 0.0;
  Real boundary_timestep_factor = 0.15;

  Real sedov_energy_cgs = 1.0e51;
  Real sedov_beta = 1.15167;
  Real sedov_start_time = 0.0;
  Real sedov_origin_x1 = -30.0;
  Real sedov_origin_x2 = 0.0;
  Real sedov_origin_x3 = 0.0;
  Real sedov_radius_at_start = 30.0;
  Real sedov_age_at_start_cgs = 0.0;
  Real ambient_mass_density_cgs = 1.0;
};

CloudCrushingData cloud_crushing;

void FatalCloudCrushingInput(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string NormalizeToken(std::string input) {
  input.erase(std::remove_if(input.begin(), input.end(),
                             [](unsigned char c) { return std::isspace(c); }),
              input.end());
  std::transform(input.begin(), input.end(), input.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  for (char &c : input) {
    if (c == '-') {
      c = '_';
    }
  }
  return input;
}

int ParseBoundaryMode(const std::string &mode_raw) {
  const std::string mode = NormalizeToken(mode_raw);
  if (mode == "sedov" || mode == "sedov_taylor" || mode == "tvns") {
    return kBoundarySedovTaylor;
  }
  if (mode == "constant" || mode == "const") {
    return kBoundaryConstant;
  }
  FatalCloudCrushingInput("Invalid problem/inner_x1_boundary '" + mode_raw +
                          "'. Expected sedov or constant.");
  return kBoundarySedovTaylor;
}

const char *BoundaryModeName(const int mode) {
  return (mode == kBoundaryConstant) ? "constant" : "sedov";
}

bool ReadProblemRealAny(ParameterInput *pin, const std::vector<std::string> &names,
                        Real &value) {
  for (const std::string &name : names) {
    if (pin->DoesParameterExist("problem", name)) {
      value = pin->GetReal("problem", name);
      return true;
    }
  }
  return false;
}

Real EquilibriumFunction(const Real temp, const Real pressure_over_k, const Real hrate) {
  return pressure_over_k*ISMCoolFn(temp)/temp - hrate;
}

Real BisectEquilibriumRoot(const Real temp_lo, const Real temp_hi,
                           const Real pressure_over_k, const Real hrate) {
  Real lo = temp_lo;
  Real hi = temp_hi;
  Real flo = EquilibriumFunction(lo, pressure_over_k, hrate);
  for (int iter = 0; iter < 80; ++iter) {
    const Real mid = std::sqrt(lo*hi);
    const Real fmid = EquilibriumFunction(mid, pressure_over_k, hrate);
    if (flo*fmid <= 0.0) {
      hi = mid;
    } else {
      lo = mid;
      flo = fmid;
    }
  }
  return std::sqrt(lo*hi);
}

void FindStableISMEquilibria(const Real pressure_over_k, const Real hrate,
                             const Real temp_min, const Real temp_max,
                             Real &temp_cold, Real &temp_warm) {
  std::vector<Real> stable_roots;
  constexpr int nbins = 4096;
  Real prev_t = temp_min;
  Real prev_f = EquilibriumFunction(prev_t, pressure_over_k, hrate);

  for (int n = 1; n <= nbins; ++n) {
    const Real frac = static_cast<Real>(n)/static_cast<Real>(nbins);
    const Real temp = temp_min*std::pow(temp_max/temp_min, frac);
    const Real f = EquilibriumFunction(temp, pressure_over_k, hrate);
    if (prev_f == 0.0 || f == 0.0 || prev_f*f < 0.0) {
      const Real root = BisectEquilibriumRoot(prev_t, temp, pressure_over_k, hrate);
      const bool positive_slope_crossing = (prev_f < 0.0 && f > 0.0);
      if (positive_slope_crossing) {
        stable_roots.push_back(root);
      }
    }
    prev_t = temp;
    prev_f = f;
  }

  if (stable_roots.size() < 2) {
    FatalCloudCrushingInput("Could not find both cold and warm stable ISM equilibria. "
                            "Adjust hydro_srcterms/hrate, problem/pressure_over_k, "
                            "or the temperature search interval.");
  }
  temp_cold = stable_roots.front();
  temp_warm = stable_roots.back();
}

PhaseState MakePhaseState(const Real temp_cgs, const CloudCrushingData &data) {
  PhaseState state;
  state.temp_cgs = temp_cgs;
  state.number_density_cgs = data.pressure_over_k/temp_cgs;
  state.density_code = state.number_density_cgs/data.n_unit;
  state.pressure_code = (data.pressure_over_k*units::Units::k_boltzmann_cgs)/
                        data.pressure_unit;
  return state;
}

void ReadCloudCrushingParameters(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) {
    FatalCloudCrushingInput("cloud_crushing requires a <hydro> block.");
  }
  if (pmbp->punit == nullptr) {
    FatalCloudCrushingInput("cloud_crushing with ISM cooling requires a <units> block.");
  }
  if (!pin->DoesBlockExist("hydro_srcterms") ||
      !pin->GetOrAddBoolean("hydro_srcterms", "ism_cooling", false)) {
    FatalCloudCrushingInput("cloud_crushing requires hydro_srcterms/ism_cooling = true.");
  }
  if (!pmbp->phydro->peos->eos_data.is_ideal) {
    FatalCloudCrushingInput("cloud_crushing requires an ideal-gas hydro EOS.");
  }

  CloudCrushingData data;
  data.pressure_over_k =
      pin->GetOrAddReal("problem", "pressure_over_k", 3162.277660168379);
  data.hrate = pin->GetReal("hydro_srcterms", "hrate");
  data.n_unit = pmbp->punit->density_cgs()/pmbp->punit->mu()/
                units::Units::atomic_mass_unit_cgs;
  data.pressure_unit = pmbp->punit->pressure_cgs();
  data.velocity_unit = pmbp->punit->velocity_cgs();
  data.time_unit = pmbp->punit->time_cgs();
  data.mu = pmbp->punit->mu();
  data.gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;

  const Real temp_min = pin->GetOrAddReal("problem", "equilibrium_temp_min", 1.0);
  const Real temp_max = pin->GetOrAddReal("problem", "equilibrium_temp_max", 1.0e7);
  Real temp_cold = 0.0;
  Real temp_warm = 0.0;
  FindStableISMEquilibria(data.pressure_over_k, data.hrate, temp_min, temp_max,
                          temp_cold, temp_warm);
  data.cold = MakePhaseState(temp_cold, data);
  data.warm = MakePhaseState(temp_warm, data);

  const RegionSize &mesh_size = pm->mesh_size;
  data.cloud_radius = pin->GetOrAddReal("problem", "cloud_radius", 1.0);
  data.cloud_smoothing_width = pin->GetOrAddReal(
      "problem", "cloud_smoothing_width", 0.1*data.cloud_radius);
  data.cloud_center_x1 = pin->GetOrAddReal("problem", "cloud_center_x1",
                                           0.5*(mesh_size.x1min + mesh_size.x1max));
  data.cloud_center_x2 = pin->GetOrAddReal("problem", "cloud_center_x2", 0.0);
  data.cloud_center_x3 = pin->GetOrAddReal("problem", "cloud_center_x3", 0.0);
  data.cloud_tracer_scalar_index =
      pin->GetOrAddInteger("problem", "cloud_tracer_scalar_index", -1);
  if (data.cloud_tracer_scalar_index < -1 ||
      data.cloud_tracer_scalar_index >= pmbp->phydro->nscalars) {
    FatalCloudCrushingInput(
        "problem/cloud_tracer_scalar_index must be -1 or index an existing hydro "
        "passive scalar.");
  }

  std::string boundary_mode = "sedov";
  if (pin->DoesParameterExist("problem", "inner_x1_boundary")) {
    boundary_mode = pin->GetString("problem", "inner_x1_boundary");
  } else if (pin->DoesParameterExist("problem", "inner_x1_bc")) {
    boundary_mode = pin->GetString("problem", "inner_x1_bc");
  } else if (pin->DoesParameterExist("problem", "boundary_mode")) {
    boundary_mode = pin->GetString("problem", "boundary_mode");
  }
  data.boundary_mode = ParseBoundaryMode(boundary_mode);
  if (pin->DoesParameterExist("problem", "constant_inner_x1") &&
      pin->GetBoolean("problem", "constant_inner_x1")) {
    data.boundary_mode = kBoundaryConstant;
  }

  data.constant_density_code = data.warm.density_code;
  data.constant_pressure_code = data.warm.pressure_code;
  data.constant_vx_code = 0.0;
  data.constant_vy_code = 0.0;
  data.constant_vz_code = 0.0;
  (void) ReadProblemRealAny(pin, {"constant_density", "constant_boundary_density",
                                  "inner_x1_density"},
                            data.constant_density_code);
  (void) ReadProblemRealAny(pin, {"constant_pressure", "constant_boundary_pressure",
                                  "inner_x1_pressure"},
                            data.constant_pressure_code);
  (void) ReadProblemRealAny(pin, {"constant_vx", "constant_v1", "constant_boundary_vx",
                                  "constant_boundary_v1", "inner_x1_vx", "inner_x1_v1"},
                            data.constant_vx_code);
  (void) ReadProblemRealAny(pin, {"constant_vy", "constant_v2", "constant_boundary_vy",
                                  "constant_boundary_v2", "inner_x1_vy", "inner_x1_v2"},
                            data.constant_vy_code);
  (void) ReadProblemRealAny(pin, {"constant_vz", "constant_v3", "constant_boundary_vz",
                                  "constant_boundary_v3", "inner_x1_vz", "inner_x1_v3"},
                            data.constant_vz_code);
  Real constant_density_cgs = 0.0;
  if (ReadProblemRealAny(pin, {"constant_density_cgs", "constant_boundary_density_cgs",
                               "inner_x1_density_cgs"},
                         constant_density_cgs)) {
    data.constant_density_code = constant_density_cgs/pmbp->punit->density_cgs();
  }
  Real constant_pressure_cgs = 0.0;
  if (ReadProblemRealAny(pin, {"constant_pressure_cgs", "constant_boundary_pressure_cgs",
                               "inner_x1_pressure_cgs"},
                         constant_pressure_cgs)) {
    data.constant_pressure_code = constant_pressure_cgs/data.pressure_unit;
  }
  Real constant_vx_cgs = 0.0;
  if (ReadProblemRealAny(pin, {"constant_vx_cgs", "constant_v1_cgs",
                               "constant_boundary_vx_cgs", "constant_boundary_v1_cgs",
                               "inner_x1_vx_cgs", "inner_x1_v1_cgs"},
                         constant_vx_cgs)) {
    data.constant_vx_code = constant_vx_cgs/data.velocity_unit;
  }
  Real constant_vy_cgs = 0.0;
  if (ReadProblemRealAny(pin, {"constant_vy_cgs", "constant_v2_cgs",
                               "constant_boundary_vy_cgs", "constant_boundary_v2_cgs",
                               "inner_x1_vy_cgs", "inner_x1_v2_cgs"},
                         constant_vy_cgs)) {
    data.constant_vy_code = constant_vy_cgs/data.velocity_unit;
  }
  Real constant_vz_cgs = 0.0;
  if (ReadProblemRealAny(pin, {"constant_vz_cgs", "constant_v3_cgs",
                               "constant_boundary_vz_cgs", "constant_boundary_v3_cgs",
                               "inner_x1_vz_cgs", "inner_x1_v3_cgs"},
                         constant_vz_cgs)) {
    data.constant_vz_code = constant_vz_cgs/data.velocity_unit;
  }
  if (data.boundary_mode == kBoundaryConstant &&
      (data.constant_density_code <= 0.0 || data.constant_pressure_code <= 0.0)) {
    FatalCloudCrushingInput("Constant inner_x1 boundary requires positive density "
                            "and pressure.");
  }
  data.boundary_timestep_factor = pin->GetOrAddReal(
      "problem", "boundary_timestep_factor", data.boundary_timestep_factor);
  if (data.boundary_timestep_factor <= 0.0) {
    FatalCloudCrushingInput("problem/boundary_timestep_factor must be positive.");
  }

  const Real origin_distance =
      pin->GetOrAddReal("problem", "sedov_origin_distance", 30.0);
  data.sedov_origin_x1 = pin->GetOrAddReal("problem", "sedov_origin_x1",
                                           mesh_size.x1min - origin_distance);
  data.sedov_origin_x2 = pin->GetOrAddReal("problem", "sedov_origin_x2", 0.0);
  data.sedov_origin_x3 = pin->GetOrAddReal("problem", "sedov_origin_x3", 0.0);
  data.sedov_energy_cgs = pin->GetOrAddReal("problem", "sedov_energy_cgs", 1.0e51);
  data.sedov_beta = pin->GetOrAddReal("problem", "sedov_beta", 1.15167);
  data.sedov_start_time = pin->GetOrAddReal("problem", "sedov_start_time", 0.0);
  data.sedov_radius_at_start = pin->GetOrAddReal("problem", "sedov_radius_at_start",
                                                 origin_distance);
  data.ambient_mass_density_cgs = data.warm.number_density_cgs*data.mu*
                                  units::Units::atomic_mass_unit_cgs;
  const Real gamma = data.gm1 + 1.0;
  if (data.boundary_mode == kBoundarySedovTaylor) {
    if (data.sedov_energy_cgs <= 0.0 || data.sedov_beta <= 0.0 ||
        data.sedov_radius_at_start <= 0.0 || data.ambient_mass_density_cgs <= 0.0) {
      FatalCloudCrushingInput("Sedov-Taylor parameters must be positive.");
    }
    if (gamma <= 1.0 || gamma >= 2.0) {
      FatalCloudCrushingInput("The analytic TVNS Sedov profile used by cloud_crushing "
                              "requires 1 < hydro/gamma < 2.");
    }
    const Real radius_start_cgs = data.sedov_radius_at_start*pmbp->punit->length_cgs();
    data.sedov_age_at_start_cgs =
        std::sqrt(std::pow(radius_start_cgs/data.sedov_beta, 5)*
                  data.ambient_mass_density_cgs/data.sedov_energy_cgs);
  }

  cloud_crushing = data;

  if (global_variable::my_rank == 0) {
    std::cout << "CloudCrushing ISM equilibrium: P/k_B=" << data.pressure_over_k
              << " K cm^-3, hrate=" << data.hrate << " erg s^-1" << std::endl;
    std::cout << " cold: T=" << data.cold.temp_cgs
              << " K, n=" << data.cold.number_density_cgs << " cm^-3" << std::endl;
    std::cout << " warm: T=" << data.warm.temp_cgs
              << " K, n=" << data.warm.number_density_cgs << " cm^-3" << std::endl;
    std::cout << " inner_x1 boundary mode=" << BoundaryModeName(data.boundary_mode)
              << std::endl;
    if (data.boundary_mode == kBoundaryConstant) {
      std::cout << " Constant inner_x1 boundary: density="
                << data.constant_density_code
                << ", pressure=" << data.constant_pressure_code
                << ", velocity=(" << data.constant_vx_code << ", "
                << data.constant_vy_code << ", " << data.constant_vz_code
                << ") in code units" << std::endl;
    } else {
      std::cout << " Sedov boundary: E=" << data.sedov_energy_cgs
                << " erg, R_start=" << data.sedov_radius_at_start
                << ", age_start=" << data.sedov_age_at_start_cgs/data.time_unit
                << " code time, full TVNS interior profile" << std::endl;
    }
    std::cout << " inner_x1 boundary timestep factor="
              << data.boundary_timestep_factor << std::endl;
    if (data.cloud_tracer_scalar_index >= 0) {
      std::cout << " original-cloud tracer scalar index="
                << data.cloud_tracer_scalar_index << std::endl;
    }
  }
}

KOKKOS_INLINE_FUNCTION
void SetHydroState(const DvceArray5D<Real> &u0, const int m, const int k, const int j,
                   const int i, const Real density, const Real pressure,
                   const Real vx, const Real vy, const Real vz, const Real gm1) {
  u0(m, IDN, k, j, i) = density;
  u0(m, IM1, k, j, i) = density*vx;
  u0(m, IM2, k, j, i) = density*vy;
  u0(m, IM3, k, j, i) = density*vz;
  u0(m, IEN, k, j, i) = pressure/gm1 +
      0.5*density*(SQR(vx) + SQR(vy) + SQR(vz));
}

KOKKOS_INLINE_FUNCTION
Real TvnsLogXi5FromLogEta(const Real log_eta, const Real gamma) {
  const Real nu1 = -(13.0*SQR(gamma) - 7.0*gamma + 12.0)/
                   ((3.0*gamma - 1.0)*(2.0*gamma + 1.0));
  const Real nu2 = 5.0*(gamma - 1.0)/(2.0*gamma + 1.0);
  const Real eta = exp(log_eta);
  const Real similarity_v = (1.0 + eta)/gamma;
  const Real a = 0.5*(gamma + 1.0)*similarity_v;
  const Real b = ((gamma + 1.0)/(7.0 - gamma))*
                 (5.0 - (3.0*gamma - 1.0)*similarity_v);
  const Real c = ((gamma + 1.0)/(gamma - 1.0))*eta;
  return -2.0*log(a) + nu1*log(b) + nu2*log(c);
}

KOKKOS_INLINE_FUNCTION
void TvnsProfile(const Real xi, const Real gamma, Real &density_ratio,
                 Real &velocity_shape, Real &sound_speed_shape) {
  const Real xi_floor = 1.0e-12;
  const Real eta_shock = (gamma - 1.0)/(gamma + 1.0);
  const Real log_eta_shock = log(eta_shock);
  Real log_eta = log_eta_shock;

  if (xi < 1.0) {
    const Real xi_limited = (xi > xi_floor) ? xi : xi_floor;
    const Real target = 5.0*log(xi_limited);
    Real lo = -700.0;
    Real hi = log_eta_shock;
    for (int iter = 0; iter < 80; ++iter) {
      const Real mid = 0.5*(lo + hi);
      if (TvnsLogXi5FromLogEta(mid, gamma) < target) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    log_eta = 0.5*(lo + hi);
  }
  const Real eta = exp(log_eta);
  const Real similarity_v = (1.0 + eta)/gamma;

  const Real nu1 = -(13.0*SQR(gamma) - 7.0*gamma + 12.0)/
                   ((3.0*gamma - 1.0)*(2.0*gamma + 1.0));
  const Real nu3 = 3.0/(2.0*gamma + 1.0);
  const Real nu4 = -nu1/(2.0 - gamma);
  const Real nu5 = -2.0/(2.0 - gamma);
  const Real compression = (gamma + 1.0)/(gamma - 1.0);
  const Real a = ((gamma + 1.0)/(gamma - 1.0))*eta;
  const Real b = ((gamma + 1.0)/(7.0 - gamma))*
                 (5.0 - (3.0*gamma - 1.0)*similarity_v);
  const Real c = ((gamma + 1.0)/(gamma - 1.0))*(1.0 - similarity_v);
  density_ratio = compression*pow(a, nu3)*pow(b, nu4)*pow(c, nu5);
  velocity_shape = similarity_v;
  sound_speed_shape = gamma*(gamma - 1.0)*(1.0 - similarity_v)*
                      SQR(similarity_v)/(2.0*eta);
}

void SedovBoundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const CloudCrushingData data = cloud_crushing;

  if (data.boundary_mode == kBoundaryConstant) {
    par_for("constant_inner_x1", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, ng-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x1) != BoundaryFlag::user) {
        return;
      }
      const int ib = is - i - 1;
      SetHydroState(u0, m, k, j, ib, data.constant_density_code,
                    data.constant_pressure_code, data.constant_vx_code,
                    data.constant_vy_code, data.constant_vz_code, data.gm1);
      for (int n = nhydro; n < nhydro + nscalars; ++n) {
        u0(m, n, k, j, ib) = 0.0;
      }
    });
    return;
  }

  Real sedov_age_cgs = data.sedov_age_at_start_cgs;
  if (pm->time > data.sedov_start_time) {
    sedov_age_cgs += (pm->time - data.sedov_start_time)*data.time_unit;
  }
  const Real shock_radius_cgs = data.sedov_beta*
      std::pow(data.sedov_energy_cgs*SQR(sedov_age_cgs)/data.ambient_mass_density_cgs,
               0.2);
  const Real shock_radius_code = shock_radius_cgs/pm->pmb_pack->punit->length_cgs();
  const Real shock_speed_code = (0.4*shock_radius_cgs/sedov_age_cgs)/data.velocity_unit;
  const Real gamma = data.gm1 + 1.0;
  Real frame_x1 = 0.0;
  Real frame_x2 = 0.0;
  Real frame_x3 = 0.0;
  Real frame_v1 = 0.0;
  Real frame_v2 = 0.0;
  Real frame_v3 = 0.0;
  if (pmbp->pframe_tracker != nullptr) {
    frame_x1 = pmbp->pframe_tracker->FrameDisplacement(0);
    frame_x2 = pmbp->pframe_tracker->FrameDisplacement(1);
    frame_x3 = pmbp->pframe_tracker->FrameDisplacement(2);
    frame_v1 = pmbp->pframe_tracker->FrameVelocity(0);
    frame_v2 = pmbp->pframe_tracker->FrameVelocity(1);
    frame_v3 = pmbp->pframe_tracker->FrameVelocity(2);
  }

  par_for("sedov_inner_x1", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (mb_bcs.d_view(m, BoundaryFace::inner_x1) != BoundaryFlag::user) {
      return;
    }
    const RegionSize block_size = size.d_view(m);
    const int ib = is - i - 1;
    const Real x1 = CellCenterX(ib - is, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, block_size.x2min, block_size.x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, block_size.x3min, block_size.x3max);
    const Real dx1 = x1 + frame_x1 - data.sedov_origin_x1;
    const Real dx2 = x2 + frame_x2 - data.sedov_origin_x2;
    const Real dx3 = x3 + frame_x3 - data.sedov_origin_x3;
    const Real radius = sqrt(SQR(dx1) + SQR(dx2) + SQR(dx3));

    if (radius <= shock_radius_code) {
      const Real xi = radius/shock_radius_code;
      const Real xi_profile = (xi > 1.0e-12) ? xi : 1.0e-12;
      Real density_ratio = 0.0;
      Real velocity_shape = 0.0;
      Real sound_speed_shape = 0.0;
      TvnsProfile(xi, gamma, density_ratio, velocity_shape, sound_speed_shape);
      const Real density = density_ratio*data.warm.density_code;
      const Real radial_speed = shock_speed_code*xi*velocity_shape;
      const Real pressure =
          density*SQR(shock_speed_code*xi_profile)*sound_speed_shape/gamma;
      Real vx = -frame_v1;
      Real vy = -frame_v2;
      Real vz = -frame_v3;
      if (radius > 0.0) {
        vx += radial_speed*dx1/radius;
        vy += radial_speed*dx2/radius;
        vz += radial_speed*dx3/radius;
      }
      SetHydroState(u0, m, k, j, ib, density, pressure, vx, vy, vz, data.gm1);
    } else {
      SetHydroState(u0, m, k, j, ib, data.warm.density_code, data.warm.pressure_code,
                    -frame_v1, -frame_v2, -frame_v3, data.gm1);
    }
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      u0(m, n, k, j, ib) = 0.0;
    }
  });
}

void CloudCrushingTimeStep(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const CloudCrushingData data = cloud_crushing;
  Real min_dx1 = static_cast<Real>(std::numeric_limits<float>::max());
  auto &size = pmbp->pmb->mb_size;
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    min_dx1 = std::min(min_dx1, size.h_view(m).dx1);
  }

  const Real gamma = data.gm1 + 1.0;
  Real signal_speed = 0.0;
  if (data.boundary_mode == kBoundaryConstant) {
    const Real cs = std::sqrt(gamma*data.constant_pressure_code/
                              data.constant_density_code);
    signal_speed = std::fabs(data.constant_vx_code) + cs;
  } else {
    Real sedov_age_cgs = data.sedov_age_at_start_cgs;
    if (pm->time > data.sedov_start_time) {
      sedov_age_cgs += (pm->time - data.sedov_start_time)*data.time_unit;
    }
    const Real shock_radius_cgs = data.sedov_beta*
        std::pow(data.sedov_energy_cgs*SQR(sedov_age_cgs)/data.ambient_mass_density_cgs,
                 0.2);
    const Real shock_speed_code = (0.4*shock_radius_cgs/sedov_age_cgs)/
                                  data.velocity_unit;
    signal_speed = 2.0*shock_speed_code;
    if (pmbp->pframe_tracker != nullptr) {
      signal_speed += std::fabs(pmbp->pframe_tracker->FrameVelocity(0));
    }
  }

  Real dtnew_loc = static_cast<Real>(std::numeric_limits<float>::max());
  if (min_dx1 > 0.0 && signal_speed > 0.0) {
    dtnew_loc = data.boundary_timestep_factor*min_dx1/signal_speed;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &dtnew_loc, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  pm->pgen->dtnew = dtnew_loc;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem()
//! \brief Initialize a pressure-balanced cold cloud in a warm ISM and enroll Sedov BCs.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  ReadCloudCrushingParameters(pin, pmy_mesh_);
  user_bcs_func = SedovBoundary;
  user_dt = true;
  user_time_step_func = CloudCrushingTimeStep;

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  const CloudCrushingData data = cloud_crushing;

  par_for("pgen_cloud_crushing", DevExeSpace(), 0, pmbp->nmb_thispack-1,
  ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const RegionSize block_size = size.d_view(m);
    const Real x1 = CellCenterX(i - is, indcs.nx1, block_size.x1min, block_size.x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, block_size.x2min, block_size.x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, block_size.x3min, block_size.x3max);
    const Real radius = sqrt(SQR(x1 - data.cloud_center_x1) +
                             SQR(x2 - data.cloud_center_x2) +
                             SQR(x3 - data.cloud_center_x3));
    Real cold_fraction = (radius < data.cloud_radius) ? 1.0 : 0.0;
    if (data.cloud_smoothing_width > 0.0) {
      cold_fraction = 0.5*(1.0 - tanh((radius - data.cloud_radius)/
                                      data.cloud_smoothing_width));
    }
    const Real density = data.warm.density_code +
                         cold_fraction*(data.cold.density_code - data.warm.density_code);
    SetHydroState(u0, m, k, j, i, density, data.warm.pressure_code,
                  0.0, 0.0, 0.0, data.gm1);
    for (int n = nhydro; n < nhydro + nscalars; ++n) {
      u0(m, n, k, j, i) = 0.0;
    }
    if (data.cloud_tracer_scalar_index >= 0) {
      u0(m, nhydro + data.cloud_tracer_scalar_index, k, j, i) =
          density*cold_fraction;
    }
  });
}
