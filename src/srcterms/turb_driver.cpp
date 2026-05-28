//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver.cpp
//  \brief implementation of functions in TurbulenceDriver

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "ion-neutral/ion-neutral.hpp"
#include "driver/driver.hpp"
#include "utils/random.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_hyd.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "turb_driver.hpp"
#include "globals.hpp"

namespace {

void FatalTurbulenceError(const std::string& message) {
  std::cout << "### FATAL ERROR in turbulence driver: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

TurbulenceDriver::TurbulenceDriver(MeshBlockPack* pp, ParameterInput* pin)
    : pmy_pack(pp),
      force("force", 1, 1, 1, 1, 1),
      mode_amp_real("mode_amp_real", 1, 1),
      mode_amp_imag("mode_amp_imag", 1, 1),
      mode_noise_real("mode_noise_real", 1, 1),
      mode_noise_imag("mode_noise_imag", 1, 1),
      kx_mode("kx_mode", 1),
      ky_mode("ky_mode", 1),
      kz_mode("kz_mode", 1),
      xcos("xcos", 1, 1, 1),
      xsin("xsin", 1, 1, 1),
      ycos("ycos", 1, 1, 1),
      ysin("ysin", 1, 1, 1),
      zcos("zcos", 1, 1, 1),
      zsin("zsin", 1, 1, 1) {
  // Allocate up to the AMR capacity, matching the evolved fluid arrays.
  int nmb = pmy_pack->nmb_thispack;
  int nmb_alloc = std::max(nmb, pmy_pack->pmesh->nmb_maxperrank);
  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int ncells1 = indcs.nx1 + 2 * (indcs.ng);
  int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2 * (indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2 * (indcs.ng)) : 1;

  // Initialize AMR tracking variables
  current_nmb_ = nmb;
  Mesh* pm = pmy_pack->pmesh;
  if (pm->adaptive && pm->pmr != nullptr) {
    last_nmb_created_ = pm->pmr->nmb_created;
    last_nmb_deleted_ = pm->pmr->nmb_deleted;
  } else {
    last_nmb_created_ = 0;
    last_nmb_deleted_ = 0;
  }

  Kokkos::realloc(force, nmb_alloc, 3, ncells3, ncells2, ncells1);

  const std::string block_name = "turb_driving";
  // Default values are written to restart parameter dumps as text. Use that
  // same representation from the first cycle so restart comparisons are exact.
  auto get_serialized_real = [pin, &block_name](const char* name, Real default_value) {
    pin->GetOrAddReal(block_name, name, default_value);
    return pin->GetReal(block_name, name);
  };
  const std::pair<const char*, const char*> removed_parameters[] = {
      {"constant_edot", "normalization = edot or normalization = accel_rms"},
      {"x_turb_scale_height", "sigma_x1"},
      {"y_turb_scale_height", "sigma_x2"},
      {"z_turb_scale_height", "sigma_x3"},
      {"x_turb_center", "center_x1"},
      {"y_turb_center", "center_x2"},
      {"z_turb_center", "center_x3"},
      {"tile_factor", "tile_nx, tile_ny, and tile_nz"},
      {"tile_driving", "tile_nx, tile_ny, and tile_nz"},
      {"dt_turb_update", "dt_update"},
      {"spect_form", "spectrum"}
  };
  for (const auto& parameter : removed_parameters) {
    if (pin->DoesParameterExist(block_name, parameter.first)) {
      FatalTurbulenceError("removed parameter '" + std::string(parameter.first) +
                           "'; use '" + std::string(parameter.second) + "'");
    }
  }

  // range of modes included, corresponding to kmin and kmax
  nlow = pin->GetOrAddInteger(block_name, "nlow", 1);
  nhigh = pin->GetOrAddInteger(block_name, "nhigh", 3);
  if (nlow < 1 || nhigh < nlow) {
    FatalTurbulenceError("nlow and nhigh must satisfy 1 <= nlow <= nhigh");
  }
  // Peak of power when spectral form is parabolic. Interpret npeak in
  // tile-local x1 mode units once the tile dimensions have been established.
  use_npeak = pin->DoesParameterExist(block_name, "npeak");
  if (use_npeak) {
    npeak = pin->GetReal(block_name, "npeak");
    kpeak = 0.0;
  } else {
    npeak = 0.0;
    kpeak = get_serialized_real("kpeak", 4.0 * M_PI);
  }
  std::string spectrum_name = pin->GetOrAddString(block_name, "spectrum", "parabolic");
  if (spectrum_name == "parabolic") {
    spectrum = TurbSpectrum::parabolic;
  } else if (spectrum_name == "power_law") {
    spectrum = TurbSpectrum::power_law;
  } else {
    FatalTurbulenceError("spectrum must be parabolic or power_law");
  }
  // driving type - 0 for 3D isotropic, 1 for planar (xy) driving
  driving_type = pin->GetOrAddInteger(block_name, "driving_type", 0);
  if (driving_type != 0 && driving_type != 1) {
    FatalTurbulenceError("driving_type must be zero or one");
  }
  // min kz zero should be 0 for including kz modes and 1 for not including
  min_kz = pin->GetOrAddInteger(block_name, "min_kz", 0);
  max_kz = pin->GetOrAddInteger(block_name, "max_kz", nhigh);
  min_kx = pin->GetOrAddInteger(block_name, "min_kx", 0);
  max_kx = pin->GetOrAddInteger(block_name, "max_kx", nhigh);
  min_ky = pin->GetOrAddInteger(block_name, "min_ky", 0);
  max_ky = pin->GetOrAddInteger(block_name, "max_ky", nhigh);
  // power-law exponent for isotropic driving
  expo = get_serialized_real("expo", 5.0 / 3.0);
  exp_prp = get_serialized_real("exp_prp", 5.0 / 3.0);
  exp_prl = get_serialized_real("exp_prl", 0.0);
  // correlation time
  tcorr = get_serialized_real("tcorr", 0.0);
  if (tcorr < 0.0) {
    FatalTurbulenceError("tcorr must not be negative");
  }
  // update time for the turbulence driver
  dt_update = get_serialized_real("dt_update", 0.01);
  // To store fraction of energy in solenoidal modes
  sol_fraction = get_serialized_real("sol_fraction", 1.0);
  if (dt_update <= 0.0) {
    FatalTurbulenceError("dt_update must be greater than zero");
  }
  if (sol_fraction < 0.0 || sol_fraction > 1.0) {
    FatalTurbulenceError("sol_fraction must lie between zero and one");
  }

  // random seed for turbulence driving
  // Non-negative values give reproducible sequences; negative values fall back
  // to the internal default (seed = 1).
  rseed = pin->GetOrAddInteger(block_name, "rseed", -1);

  std::string normalization_name =
      pin->GetOrAddString(block_name, "normalization", "edot");
  bool has_dedt = pin->DoesParameterExist(block_name, "dedt");
  bool has_accel_rms = pin->DoesParameterExist(block_name, "accel_rms");
  if (normalization_name == "edot") {
    normalization = TurbNormalization::edot;
    if (!has_dedt) {
      FatalTurbulenceError("normalization = edot requires dedt");
    }
    if (has_accel_rms) {
      FatalTurbulenceError("accel_rms is not used with normalization = edot");
    }
    dedt = pin->GetReal(block_name, "dedt");
    accel_rms = 0.0;
    if (dedt < 0.0) {
      FatalTurbulenceError("dedt must not be negative");
    }
  } else if (normalization_name == "accel_rms") {
    normalization = TurbNormalization::accel_rms;
    if (!has_accel_rms) {
      FatalTurbulenceError("normalization = accel_rms requires accel_rms");
    }
    if (has_dedt) {
      FatalTurbulenceError("dedt is not used with normalization = accel_rms");
    }
    dedt = 0.0;
    accel_rms = pin->GetReal(block_name, "accel_rms");
    if (accel_rms < 0.0) {
      FatalTurbulenceError("accel_rms must not be negative");
    }
  } else {
    FatalTurbulenceError("normalization must be edot or accel_rms");
  }

  sigma_x1 = get_serialized_real("sigma_x1", -1.0);
  sigma_x2 = get_serialized_real("sigma_x2", -1.0);
  sigma_x3 = get_serialized_real("sigma_x3", -1.0);
  center_x1 = get_serialized_real("center_x1", 0.0);
  center_x2 = get_serialized_real("center_x2", 0.0);
  center_x3 = get_serialized_real("center_x3", 0.0);
  std::string localization_name =
      pin->GetOrAddString(block_name, "localization", "none");
  if (localization_name == "none") {
    localization = TurbLocalization::none;
  } else if (localization_name == "include") {
    localization = TurbLocalization::include;
  } else if (localization_name == "exclude") {
    localization = TurbLocalization::exclude;
  } else {
    FatalTurbulenceError("localization must be none, include, or exclude");
  }
  bool has_envelope = (sigma_x1 > 0.0 || sigma_x2 > 0.0 || sigma_x3 > 0.0);
  if (localization == TurbLocalization::none && has_envelope) {
    FatalTurbulenceError("sigma_x* requires localization = include or exclude");
  } else if (localization != TurbLocalization::none && !has_envelope) {
    FatalTurbulenceError("localization requires a positive sigma_x* value");
  }

  tile_nx = pin->GetOrAddInteger(block_name, "tile_nx", 1);
  tile_ny = pin->GetOrAddInteger(block_name, "tile_ny", 1);
  tile_nz = pin->GetOrAddInteger(block_name, "tile_nz", 1);

  domain_x1min = pm->mesh_size.x1min;
  domain_x2min = pm->mesh_size.x2min;
  domain_x3min = pm->mesh_size.x3min;

  auto& mesh_indcs_root = pm->mesh_indcs;
  if (tile_nx < 1 || tile_ny < 1 || tile_nz < 1) {
    FatalTurbulenceError("tile counts must be greater than or equal to one");
  }

  if (mesh_indcs_root.nx1 % tile_nx != 0) {
    FatalTurbulenceError("tile_nx must evenly divide nx1");
  }
  if (mesh_indcs_root.nx2 <= 1) {
    if (tile_ny != 1) {
      FatalTurbulenceError("tile_ny must be one for a one-dimensional x2 grid");
    }
  } else if (mesh_indcs_root.nx2 % tile_ny != 0) {
    FatalTurbulenceError("tile_ny must evenly divide nx2");
  }
  if (mesh_indcs_root.nx3 <= 1) {
    if (tile_nz != 1) {
      FatalTurbulenceError("tile_nz must be one for a one-dimensional x3 grid");
    }
  } else if (mesh_indcs_root.nx3 % tile_nz != 0) {
    FatalTurbulenceError("tile_nz must evenly divide nx3");
  }

  Real lx_global = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly_global = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz_global = pm->mesh_size.x3max - pm->mesh_size.x3min;

  tile_lx = lx_global / static_cast<Real>(tile_nx);
  tile_ly = (tile_ny > 0) ? ly_global / static_cast<Real>(tile_ny) : ly_global;
  tile_lz = (tile_nz > 0) ? lz_global / static_cast<Real>(tile_nz) : lz_global;
  if (use_npeak) {
    kpeak = npeak * 2.0 * M_PI / tile_lx;
  }
  if (spectrum == TurbSpectrum::parabolic) {
    if (nhigh == nlow) {
      FatalTurbulenceError("parabolic spectrum requires nhigh greater than nlow");
    }
    if (use_npeak && (npeak < nlow || npeak > nhigh)) {
      FatalTurbulenceError("npeak must lie between nlow and nhigh");
    }
  }

  inv_tile_lx = (tile_lx > 0.0) ? 1.0 / tile_lx : 0.0;
  inv_tile_ly = (tile_ly > 0.0) ? 1.0 / tile_ly : 0.0;
  inv_tile_lz = (tile_lz > 0.0) ? 1.0 / tile_lz : 0.0;

  // decaying/constant energy injection - 1 for decaying, 2 continuously driven
  turb_flag = pin->GetOrAddInteger(block_name, "turb_flag", 2);
  if (turb_flag == 1) {
    tdriv_duration =
        get_serialized_real("tdriv_duration",
                            tcorr);  // If not specified, drive for one correlation time
  } else {
    tdriv_duration = static_cast<Real>(
        std::numeric_limits<float>::max());  // For constantly stirred turbulence, set
                                             // this to float max
  }
  tdriv_start = get_serialized_real(
      "tdriv_start", 0.);  // If not specified, start driving at t=0
  if (global_variable::my_rank == 0) {
    std::cout << "Initialising turbulence driving module" << std::endl
              << " dedt = " << dedt << " tcorr = " << tcorr
              << " dt_update = " << dt_update << std::endl;
  }
  n_turb_updates_yet = 0;

  Real nlow_sqr = nlow * nlow;
  Real nhigh_sqr = nhigh * nhigh;

  mode_count = 0;

  // Count Cartesian modes
  int nkx, nky, nkz;
  Real nsqr;
  for (nkx = min_kx; nkx <= max_kx; nkx++) {
    for (nky = min_ky; nky <= max_ky; nky++) {
      for (nkz = min_kz; nkz <= max_kz; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        nsqr = 0.0;
        bool flag_prl = true;
        if (driving_type == 0) {
          nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);
        } else if (driving_type == 1) {
          nsqr = SQR(nkx) + SQR(nky);
          Real nprlsqr = SQR(nkz);
          if (nprlsqr >= nlow_sqr && nprlsqr <= nhigh_sqr) {
            flag_prl = true;
          } else {
            flag_prl = false;
          }
        }
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr && flag_prl) {
          mode_count++;
        }
      }
    }
  }

  if (mode_count == 0) {
    std::cout << "ERROR: mode_count is 0! Check turbulence driving parameters."
              << std::endl;
    std::cout << "  nlow=" << nlow << ", nhigh=" << nhigh << std::endl;
    std::cout << "  driving_type=" << driving_type << std::endl;
    exit(EXIT_FAILURE);
  }

  Kokkos::realloc(mode_amp_real, 3, mode_count);
  Kokkos::realloc(mode_amp_imag, 3, mode_count);
  Kokkos::realloc(mode_noise_real, 3, mode_count);
  Kokkos::realloc(mode_noise_imag, 3, mode_count);

  // Allocate Cartesian mode arrays
  Kokkos::realloc(kx_mode, mode_count);
  Kokkos::realloc(ky_mode, mode_count);
  Kokkos::realloc(kz_mode, mode_count);

  Kokkos::realloc(xcos, nmb_alloc, mode_count, ncells1);
  Kokkos::realloc(xsin, nmb_alloc, mode_count, ncells1);
  Kokkos::realloc(ycos, nmb_alloc, mode_count, ncells2);
  Kokkos::realloc(ysin, nmb_alloc, mode_count, ncells2);
  Kokkos::realloc(zcos, nmb_alloc, mode_count, ncells3);
  Kokkos::realloc(zsin, nmb_alloc, mode_count, ncells3);

  Initialize();
}

//----------------------------------------------------------------------------------------
// destructor

TurbulenceDriver::~TurbulenceDriver() {}

//----------------------------------------------------------------------------------------
//! \fn  noid Initialize
//  \brief Function to initialize the driver

void TurbulenceDriver::Initialize() {
  int nmb = pmy_pack->nmb_thispack;
  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int ncells1 = indcs.nx1 + 2 * (indcs.ng);
  int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2 * (indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2 * (indcs.ng)) : 1;

  auto force_ = force;
  par_for(
      "force_init_pgen", DevExeSpace(), 0, nmb - 1, 0, 2, 0, ncells3 - 1, 0, ncells2 - 1,
      0, ncells1 - 1, KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
        force_(m, n, k, j, i) = 0.0;
      });
  for (int dir = 0; dir < 3; ++dir) {
    for (int n = 0; n < mode_count; ++n) {
      mode_amp_real.h_view(dir, n) = 0.0;
      mode_amp_imag.h_view(dir, n) = 0.0;
      mode_noise_real.h_view(dir, n) = 0.0;
      mode_noise_imag.h_view(dir, n) = 0.0;
    }
  }
  mode_amp_real.template modify<HostMemSpace>();
  mode_amp_real.template sync<DevExeSpace>();
  mode_amp_imag.template modify<HostMemSpace>();
  mode_amp_imag.template sync<DevExeSpace>();

  // Initialize RNG state for the Ornstein-Uhlenbeck forcing. Use a negative
  // idum so that RanSt() takes the initialization branch on first use.
  if (rseed >= 0) {
    // Non-negative seeds give reproducible sequences; treat 0 as 1.
    int seed = (rseed > 0) ? rseed : 1;
    rstate.idum = -static_cast<decltype(rstate.idum)>(seed);
  } else {
    // Negative rseed falls back to internal default seed = 1.
    rstate.idum = -1;
  }
  rstate.iset = 0;

  auto kx_mode_ = kx_mode;
  auto ky_mode_ = ky_mode;
  auto kz_mode_ = kz_mode;

  // Cartesian plane-wave precomputations
  Real dkx, dky, dkz, kx, ky, kz;
  Real lx = tile_lx;
  Real ly = tile_ly;
  Real lz = tile_lz;
  dkx = (lx > 0.0) ? 2.0 * M_PI / lx : 0.0;
  dky = (ly > 0.0) ? 2.0 * M_PI / ly : 0.0;
  dkz = (lz > 0.0) ? 2.0 * M_PI / lz : 0.0;

  int nmode = 0;
  int nkx, nky, nkz;
  Real nsqr;
  Real nlow_sqr = nlow * nlow;
  Real nhigh_sqr = nhigh * nhigh;
  for (nkx = min_kx; nkx <= max_kx; nkx++) {
    for (nky = min_ky; nky <= max_ky; nky++) {
      for (nkz = min_kz; nkz <= max_kz; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        nsqr = 0.0;
        bool flag_prl = true;
        if (driving_type == 0) {
          nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);
        } else if (driving_type == 1) {
          nsqr = SQR(nkx) + SQR(nky);
          Real nprlsqr = SQR(nkz);
          if (nprlsqr >= nlow_sqr && nprlsqr <= nhigh_sqr) {
            flag_prl = true;
          } else {
            flag_prl = false;
          }
        }
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr && flag_prl) {
          kx = dkx * nkx;
          ky = dky * nky;
          kz = dkz * nkz;
          kx_mode_.h_view(nmode) = kx;
          ky_mode_.h_view(nmode) = ky;
          kz_mode_.h_view(nmode) = kz;
          nmode++;
        }
      }
    }
  }

  kx_mode_.template modify<HostMemSpace>();
  kx_mode_.template sync<DevExeSpace>();
  ky_mode_.template modify<HostMemSpace>();
  ky_mode_.template sync<DevExeSpace>();
  kz_mode_.template modify<HostMemSpace>();
  kz_mode_.template sync<DevExeSpace>();

  BuildBasis();
}

//----------------------------------------------------------------------------------------
//! \fn BuildBasis()
// \brief Render the geometry-dependent trigonometric basis for current MeshBlocks.

void TurbulenceDriver::BuildBasis() {
  int nmb = pmy_pack->nmb_thispack;
  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2 * indcs.ng) : 1;
  int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2 * indcs.ng) : 1;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;

  auto kx_mode_ = kx_mode;
  auto ky_mode_ = ky_mode;
  auto kz_mode_ = kz_mode;
  auto xcos_ = xcos;
  auto xsin_ = xsin;
  auto ycos_ = ycos;
  auto ysin_ = ysin;
  auto zcos_ = zcos;
  auto zsin_ = zsin;

  auto size_view = pmy_pack->pmb->mb_size;
  size_view.template modify<HostMemSpace>();
  size_view.template sync<DevExeSpace>();
  const int drivingtype = driving_type;
  const bool tile_enabled = (tile_nx > 1 || tile_ny > 1 || tile_nz > 1);
  const int tile_nx_local = tile_nx;
  const int tile_ny_local = tile_ny;
  const int tile_nz_local = tile_nz;
  const Real tile_lx_local = tile_lx;
  const Real tile_ly_local = tile_ly;
  const Real tile_lz_local = tile_lz;
  const Real inv_tile_lx_local = inv_tile_lx;
  const Real inv_tile_ly_local = inv_tile_ly;
  const Real inv_tile_lz_local = inv_tile_lz;
  const Real domain_x1min_local = domain_x1min;
  const Real domain_x2min_local = domain_x2min;
  const Real domain_x3min_local = domain_x3min;

  par_for(
      "turb_basis_x", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, is, ie,
      KOKKOS_LAMBDA(int m, int n, int i) {
        Real x1v = CellCenterX(i - is, nx1, size_view.d_view(m).x1min,
                               size_view.d_view(m).x1max);
        Real arg = x1v;
        if (tile_enabled && tile_nx_local > 1) {
          Real rel = x1v - domain_x1min_local;
          int tile_i = static_cast<int>(floor(rel * inv_tile_lx_local));
          tile_i =
              (tile_i < 0) ? 0 : ((tile_i >= tile_nx_local) ? tile_nx_local - 1 : tile_i);
          arg = x1v - (domain_x1min_local + tile_i * tile_lx_local);
        }
        xsin_(m, n, i) = sin(kx_mode_.d_view(n) * arg);
        xcos_(m, n, i) = cos(kx_mode_.d_view(n) * arg);
      });

  par_for(
      "turb_basis_y", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, js, je,
      KOKKOS_LAMBDA(int m, int n, int j) {
        Real x2v = CellCenterX(j - js, nx2, size_view.d_view(m).x2min,
                               size_view.d_view(m).x2max);
        Real arg = x2v;
        if (tile_enabled && tile_ny_local > 1) {
          Real rel = x2v - domain_x2min_local;
          int tile_j = static_cast<int>(floor(rel * inv_tile_ly_local));
          tile_j =
              (tile_j < 0) ? 0 : ((tile_j >= tile_ny_local) ? tile_ny_local - 1 : tile_j);
          arg = x2v - (domain_x2min_local + tile_j * tile_ly_local);
        }
        ysin_(m, n, j) = sin(ky_mode_.d_view(n) * arg);
        ycos_(m, n, j) = cos(ky_mode_.d_view(n) * arg);
        if (ncells2 == 1) {
          ysin_(m, n, j) = 0.0;
          ycos_(m, n, j) = 1.0;
        }
      });

  par_for(
      "turb_basis_z", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, ks, ke,
      KOKKOS_LAMBDA(int m, int n, int k) {
        Real x3v = CellCenterX(k - ks, nx3, size_view.d_view(m).x3min,
                               size_view.d_view(m).x3max);
        Real arg = x3v;
        if (tile_enabled && tile_nz_local > 1) {
          Real rel = x3v - domain_x3min_local;
          int tile_k = static_cast<int>(floor(rel * inv_tile_lz_local));
          tile_k =
              (tile_k < 0) ? 0 : ((tile_k >= tile_nz_local) ? tile_nz_local - 1 : tile_k);
          arg = x3v - (domain_x3min_local + tile_k * tile_lz_local);
        }
        zsin_(m, n, k) = sin(kz_mode_.d_view(n) * arg);
        zcos_(m, n, k) = cos(kz_mode_.d_view(n) * arg);
        if (ncells3 == 1 || drivingtype == 1) {
          zsin_(m, n, k) = 0.0;
          zcos_(m, n, k) = 1.0;
        }
      });
}

//----------------------------------------------------------------------------------------
//! \fn  void IncludeModeEvolutionTasks
//  \brief Includes task in the operator split task list that constructs new modes with
//  random amplitudes and phases that can be used to evolve the force via an O-U process
//  Called by MeshBlockPack::AddPhysics() function

void TurbulenceDriver::IncludeInitializeModesTask(std::shared_ptr<TaskList> tl,
                                                  TaskID start) {
  //  We check for mesh changes, then initialize modes and update the forcing
  auto id_resize = tl->AddTask(&TurbulenceDriver::EnsureBasisSize, this, start);
  auto id_init = tl->AddTask(&TurbulenceDriver::InitializeModes, this, id_resize);
  auto id_add = tl->AddTask(&TurbulenceDriver::UpdateForcing, this, id_init);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void IncludeForcingTasks
//  \brief includes task in the stage_run task list for adding random forcing to fluid
//  as an explicit source terms in each stage of integrator
//  Called by MeshBlockPack::AddPhysics() function

void TurbulenceDriver::IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start) {
  // These must be inserted after update task, but before the source terms
  // We apply the forcing in each step of the time integration,
  // note that we do not update the forcing in each RK stage
  if (pmy_pack->pionn == nullptr) {
    if (pmy_pack->phydro != nullptr) {
      auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                               pmy_pack->phydro->id.rkupdt, pmy_pack->phydro->id.srctrms);
    }
    if (pmy_pack->pmhd != nullptr) {
      auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                               pmy_pack->pmhd->id.rkupdt, pmy_pack->pmhd->id.srctrms);
    }
  } else {
    auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                             pmy_pack->pionn->id.n_rkupdt, pmy_pack->pionn->id.n_flux);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn InitializeModes()
// \brief Evolve the modal OU state to the update index required at the current time.

TaskStatus TurbulenceDriver::InitializeModes(Driver* pdrive, int stage) {
  if (pmy_pack == nullptr) {
    return TaskStatus::complete;
  }

  Mesh* pm = pmy_pack->pmesh;
  if (pm == nullptr) {
    return TaskStatus::complete;
  }

  Real current_time = pm->time;
  if (current_time < tdriv_start) return TaskStatus::complete;
  Real t_since_start = current_time - tdriv_start;
  int n_turb_updates_reqd = static_cast<int>(t_since_start / dt_update) + 1;

  int nlow_sqr = SQR(nlow);
  int nhigh_sqr = SQR(nhigh);

  auto mode_amp_real_ = mode_amp_real;
  auto mode_amp_imag_ = mode_amp_imag;
  auto mode_noise_real_ = mode_noise_real;
  auto mode_noise_imag_ = mode_noise_imag;

  Real dkx, dky, dkz, kx, ky, kz;
  Real lx = tile_lx;
  Real ly = tile_ly;
  Real lz = tile_lz;
  dkx = (lx > 0.0) ? 2.0 * M_PI / lx : 0.0;
  dky = (ly > 0.0) ? 2.0 * M_PI / ly : 0.0;
  dkz = (lz > 0.0) ? 2.0 * M_PI / lz : 0.0;

  Real& ex = expo;
  Real& ex_prp = exp_prp;
  Real& ex_prl = exp_prl;
  Real norm, kprl, kprp, kiso;
  Real khigh = nhigh * fmax(fmax(dkx, dky), dkz);
  Real klow = nlow * fmin(fmin(dkx, dky), dkz);
  Real parab_prefact = 0.0;
  if (spectrum == TurbSpectrum::parabolic) {
    parab_prefact = -4.0 / pow(khigh - klow, 2.0);
  }
  Real& k_peak = kpeak;

  // Now compute new force using new random amplitudes and phases
  // Advance modal state only when another configured OU update boundary is reached.

  if ((t_since_start < tdriv_duration) ||
      turb_flag != 1) {  // Update forcing if continuous or t<tdriv_duration
    for (int i_turb_update = n_turb_updates_yet; i_turb_update < n_turb_updates_reqd;
         i_turb_update++) {
      int no_dir = 3;
      int nmode = 0;

      // Cartesian mode generation
      int nkx, nky, nkz, nsqr;

      for (nkx = min_kx; nkx <= max_kx; nkx++) {
        for (nky = min_ky; nky <= max_ky; nky++) {
          for (nkz = min_kz; nkz <= max_kz; nkz++) {
            if (nkx == 0 && nky == 0 && nkz == 0) continue;
            norm = 0.0;
            nsqr = 0.0;
            bool flag_prl = true;
            if (driving_type == 0) {
              nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);
            } else if (driving_type == 1) {
              nsqr = SQR(nkx) + SQR(nky);
              Real nprlsqr = SQR(nkz);
              if (nprlsqr >= nlow_sqr && nprlsqr <= nhigh_sqr) {
                flag_prl = true;
              } else {
                flag_prl = false;
              }
            }
            if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr && flag_prl) {
              kx = dkx * nkx;
              ky = dky * nky;
              kz = dkz * nkz;

              Real k[3] = {kx, ky, kz};
              // Always define kiso; used below for the solenoidal/compressive split
              kiso = sqrt(SQR(kx) + SQR(ky) + SQR(kz));

              // Generate Fourier amplitudes

              if (driving_type == 0) {
                if (kiso > 1e-16) {
                  if (spectrum == TurbSpectrum::power_law) {
                    norm = 1.0 / pow(kiso, (ex + 2.0) / 2.0);  // power-law driving
                  } else if (spectrum == TurbSpectrum::parabolic) {
                    norm = fabs(parab_prefact * pow(kiso - k_peak, 2.0) +
                                1.0);  // parabola in k-space
                    norm = pow(norm, 0.5) * pow(k_peak / kiso, (no_dir - 1) / 2.0);
                  } else {
                    norm = 0.0;
                  }
                } else {
                  norm = 0.0;
                }
              } else if (driving_type == 1) {
                no_dir = 2;
                kprl = sqrt(SQR(kx));
                kprp = sqrt(SQR(ky) + SQR(kz));
                if (kprl > 1e-16 && kprp > 1e-16) {
                  if (spectrum == TurbSpectrum::power_law) {
                    norm =
                        1.0 / pow(kprp, (ex_prp + 1.0) / 2.0) / pow(kprl, ex_prl / 2.0);
                  } else if (spectrum == TurbSpectrum::parabolic) {
                    norm = fabs(parab_prefact * pow(kprp - k_peak, 2.0) +
                                1.0);  // parabola in kperp-space
                    norm = pow(norm, 0.5) * pow(k_peak / kprp, (no_dir - 1) / 2.0);
                  }
                } else {
                  norm = 0.0;
                }
              }
              // Generate complex Fourier amplitudes for this mode:
              //   amp_real_dir (real part) and amp_imag_dir (imaginary part),
              // scaled by norm. Also accumulate k·Re(A) and k·Im(A) to construct
              // solenoidal/compressive projections below.
              Real k_dot_amp_imag = 0.0;
              Real k_dot_amp_real = 0.0;

              for (int dir = 0; dir < 3; dir++) {
                mode_noise_real_.h_view(dir, nmode) = 0.0;
                mode_noise_imag_.h_view(dir, nmode) = 0.0;
              }
              for (int dir = 0; dir < no_dir; dir++) {
                Real amp_real_dir = norm * RanGaussianSt(&(rstate));
                Real amp_imag_dir = norm * RanGaussianSt(&(rstate));
                mode_noise_real_.h_view(dir, nmode) = amp_real_dir;
                mode_noise_imag_.h_view(dir, nmode) = amp_imag_dir;

                k_dot_amp_imag += k[dir] * amp_imag_dir;  // k·Im(A)
                k_dot_amp_real += k[dir] * amp_real_dir;  // k·Re(A)
              }

              // Now decompose into solenoidal/compressive modes.
              if (norm > 0.) {
                for (int dir = 0; dir < no_dir; dir++) {
                  // Compressible (longitudinal) projections:
                  //   A_div = k (k·Re(A)) / |k|^2,  B_div = k (k·Im(A)) / |k|^2
                  Real A_div = k[dir] * k_dot_amp_real / SQR(kiso);
                  Real B_div = k[dir] * k_dot_amp_imag / SQR(kiso);

                  // Solenoidal parts (divergence-free):
                  //   A_sol = A - A_div,  B_sol = B - B_div
                  Real A_sol = mode_noise_real_.h_view(dir, nmode) - A_div;
                  Real B_sol = mode_noise_imag_.h_view(dir, nmode) - B_div;

                  // Blend in amplitude-space:
                  //   sol_fraction = 1.0 -> purely solenoidal,
                  //   sol_fraction = 0.0 -> purely compressive.
                  mode_noise_real_.h_view(dir, nmode) =
                      sol_fraction * A_sol + (1.0 - sol_fraction) * A_div;
                  mode_noise_imag_.h_view(dir, nmode) =
                      sol_fraction * B_sol + (1.0 - sol_fraction) * B_div;
                }
              }

              nmode++;
            }
          }
        }
      }

      mode_noise_real_.template modify<HostMemSpace>();
      mode_noise_imag_.template modify<HostMemSpace>();

      Real fcorr, gcorr;
      if ((tcorr <= 1e-6) || (i_turb_update == 0)) {  // use whitenoise
        fcorr = 0.0;
        gcorr = 1.0;
      } else {
        fcorr = std::exp(-dt_update / tcorr);
        gcorr = std::sqrt(1.0 - fcorr * fcorr);
      }
      // Modal coefficients are the OU state. Rendering them on another mesh therefore
      // cannot change the random process or consume additional random numbers.
      for (int dir = 0; dir < 3; ++dir) {
        for (int n = 0; n < mode_count; ++n) {
          mode_amp_real_.h_view(dir, n) =
              fcorr * mode_amp_real_.h_view(dir, n) +
              gcorr * mode_noise_real_.h_view(dir, n);
          mode_amp_imag_.h_view(dir, n) =
              fcorr * mode_amp_imag_.h_view(dir, n) +
              gcorr * mode_noise_imag_.h_view(dir, n);
        }
      }
      mode_amp_real_.template modify<HostMemSpace>();
      mode_amp_real_.template sync<DevExeSpace>();
      mode_amp_imag_.template modify<HostMemSpace>();
      mode_amp_imag_.template sync<DevExeSpace>();
    }  // end of for loop over i_turb_update
  }
  n_turb_updates_yet = n_turb_updates_reqd;
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn RenderForce
//  \brief render the authoritative modal OU state on the current MeshBlock geometry

void TurbulenceDriver::RenderForce() {
  auto& indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  auto force_ = force;
  auto mode_amp_real_ = mode_amp_real;
  auto mode_amp_imag_ = mode_amp_imag;
  auto xcos_ = xcos;
  auto xsin_ = xsin;
  auto ycos_ = ycos;
  auto ysin_ = ysin;
  auto zcos_ = zcos;
  auto zsin_ = zsin;
  const int mode_count_ = mode_count;

  mode_amp_real_.template sync<DevExeSpace>();
  mode_amp_imag_.template sync<DevExeSpace>();
  par_for(
      "turb_force_zero", DevExeSpace(), 0, nmb - 1, 0, 2, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int dir, int k, int j, int i) {
        force_(m, dir, k, j, i) = 0.0;
      });
  par_for(
      "turb_force_render", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        for (int n = 0; n < mode_count_; ++n) {
          Real forc_real =
              (xcos_(m, n, i) * ycos_(m, n, j) - xsin_(m, n, i) * ysin_(m, n, j)) *
                  zcos_(m, n, k) -
              (xsin_(m, n, i) * ycos_(m, n, j) + xcos_(m, n, i) * ysin_(m, n, j)) *
                  zsin_(m, n, k);
          Real forc_imag =
              (ycos_(m, n, j) * zsin_(m, n, k) + ysin_(m, n, j) * zcos_(m, n, k)) *
                  xcos_(m, n, i) +
              (ycos_(m, n, j) * zcos_(m, n, k) - ysin_(m, n, j) * zsin_(m, n, k)) *
                  xsin_(m, n, i);
          for (int dir = 0; dir < 3; ++dir) {
            force_(m, dir, k, j, i) += mode_amp_real_.d_view(dir, n) * forc_real -
                                        mode_amp_imag_.d_view(dir, n) * forc_imag;
          }
        }
      });
}

//----------------------------------------------------------------------------------------
//! \fn UpdateForcing
//  \brief render, localize, remove net acceleration, and normalize one forcing field
//

TaskStatus TurbulenceDriver::UpdateForcing(Driver* pdrive, int stage) {
  if (pmy_pack == nullptr) {
    return TaskStatus::complete;
  }

  Mesh* pm = pmy_pack->pmesh;
  if (pm == nullptr) {
    return TaskStatus::complete;
  }

  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  int& nx1 = indcs.nx1;
  int& nx2 = indcs.nx2;
  int& nx3 = indcs.nx3;

  Real dt = pm->dt;
  Real current_time = pm->time;
  Real t_since_start = current_time - tdriv_start;

  DvceArray5D<Real> u0, u0_;
  if (pmy_pack->phydro != nullptr) u0 = (pmy_pack->phydro->u0);
  if (pmy_pack->pmhd != nullptr) u0 = (pmy_pack->pmhd->u0);
  bool flag_twofl = false;
  if (pmy_pack->pionn != nullptr) {
    u0 = (pmy_pack->phydro->u0);
    u0_ = (pmy_pack->pmhd->u0);
    flag_twofl = true;
  }

  auto force_ = force;

  const int nmkji = nmb * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;
  // Copy the DualView handle by value, sync device, and use this in all kernels
  auto mb_size = pmy_pack->pmb->mb_size;
  mb_size.template modify<HostMemSpace>();
  mb_size.template sync<DevExeSpace>();

  const Real sigma_x1_ = sigma_x1;
  const Real sigma_x2_ = sigma_x2;
  const Real sigma_x3_ = sigma_x3;
  const Real center_x1_ = center_x1;
  const Real center_x2_ = center_x2;
  const Real center_x3_ = center_x3;
  const TurbLocalization localization_ = localization;

  if ((pm->ncycle >= 1) && (current_time >= tdriv_start) &&
      ((t_since_start < tdriv_duration) || turb_flag != 1)) {
    RenderForce();

    if (localization_ != TurbLocalization::none) {
      par_for(
          "force_localization", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
            Real exponent = 0.0;
            if (sigma_x1_ > 0.0) {
              Real x1v = CellCenterX(i - is, nx1, mb_size.d_view(m).x1min,
                                     mb_size.d_view(m).x1max);
              exponent += SQR(x1v - center_x1_) / (2.0 * SQR(sigma_x1_));
            }
            if (sigma_x2_ > 0.0) {
              Real x2v = CellCenterX(j - js, nx2, mb_size.d_view(m).x2min,
                                     mb_size.d_view(m).x2max);
              exponent += SQR(x2v - center_x2_) / (2.0 * SQR(sigma_x2_));
            }
            if (sigma_x3_ > 0.0) {
              Real x3v = CellCenterX(k - ks, nx3, mb_size.d_view(m).x3min,
                                     mb_size.d_view(m).x3max);
              exponent += SQR(x3v - center_x3_) / (2.0 * SQR(sigma_x3_));
            }
            Real gaussian = std::exp(-exponent);
            Real weight = (localization_ == TurbLocalization::include) ?
                          gaussian : (1.0 - gaussian);
            force_(m, 0, k, j, i) *= weight;
            force_(m, 1, k, j, i) *= weight;
            force_(m, 2, k, j, i) *= weight;
          });
    }

    Real t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
    Kokkos::parallel_reduce(
        "net_mom_1", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
        KOKKOS_LAMBDA(const int& idx, Real& sum_t0, Real& sum_t1, Real& sum_t2,
                      Real& sum_t3) {
          // compute n,k,j,i indices of thread
          int m = (idx) / nkji;
          int k = (idx - m * nkji) / nji;
          int j = (idx - m * nkji - k * nji) / nx1;
          int i = (idx - m * nkji - k * nji - j * nx1) + is;
          k += ks;
          j += js;
          Real vol =
              mb_size.d_view(m).dx1 * mb_size.d_view(m).dx2 * mb_size.d_view(m).dx3;
          Real den = u0(m, IDN, k, j, i);
          if (flag_twofl) {
            den += u0_(m, IDN, k, j, i);
          }
          sum_t0 += den * vol;
          sum_t1 += den * force_(m, 0, k, j, i) * vol;
          sum_t2 += den * force_(m, 1, k, j, i) * vol;
          sum_t3 += den * force_(m, 2, k, j, i) * vol;
        },
        Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1), Kokkos::Sum<Real>(t2),
        Kokkos::Sum<Real>(t3));

#if MPI_PARALLEL_ENABLED
    Real m[4], gm[4];
    m[0] = t0;
    m[1] = t1;
    m[2] = t2;
    m[3] = t3;
    MPI_Allreduce(m, gm, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    t0 = gm[0];
    t1 = gm[1];
    t2 = gm[2];
    t3 = gm[3];
#endif

    if (t0 <= 0.0) {
      FatalTurbulenceError("mass integral is not positive while normalizing forcing");
    }
    par_for(
        "force_remove_net_mom", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          force_(m, 0, k, j, i) -= t1 / t0;
          force_(m, 1, k, j, i) -= t2 / t0;
          force_(m, 2, k, j, i) -= t3 / t0;
        });

    t0 = 0.0;
    t1 = 0.0;
    Real totvol = 0.0;
    bool normalize_edot = (normalization == TurbNormalization::edot);
    Kokkos::parallel_reduce(
        "net_mom_2", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
        KOKKOS_LAMBDA(const int& idx, Real& sum_t0, Real& sum_t1, Real& totvol_) {
          // compute n,k,j,i indices of thread
          int m = (idx) / nkji;
          int k = (idx - m * nkji) / nji;
          int j = (idx - m * nkji - k * nji) / nx1;
          int i = (idx - m * nkji - k * nji - j * nx1) + is;
          k += ks;
          j += js;
          Real vol =
              mb_size.d_view(m).dx1 * mb_size.d_view(m).dx2 * mb_size.d_view(m).dx3;

          Real den = u0(m, IDN, k, j, i);
          Real mom1 = u0(m, IM1, k, j, i);
          Real mom2 = u0(m, IM2, k, j, i);
          Real mom3 = u0(m, IM3, k, j, i);
          if (flag_twofl) {
            den += u0_(m, IDN, k, j, i);
            mom1 += u0_(m, IM1, k, j, i);
            mom2 += u0_(m, IM2, k, j, i);
            mom3 += u0_(m, IM3, k, j, i);
          }
          Real a1 = force_(m, 0, k, j, i);
          Real a2 = force_(m, 1, k, j, i);
          Real a3 = force_(m, 2, k, j, i);

          if (normalize_edot) {
            sum_t0 += den * 0.5 * (a1 * a1 + a2 * a2 + a3 * a3) * dt * vol;
            sum_t1 += (mom1 * a1 + mom2 * a2 + mom3 * a3) * vol;
          } else {
            sum_t0 += (a1 * a1 + a2 * a2 + a3 * a3) * vol;
          }
          totvol_ += vol;
        },
        Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1), Kokkos::Sum<Real>(totvol));

#if MPI_PARALLEL_ENABLED
    m[0] = t0;
    m[1] = t1;
    m[2] = totvol;
    MPI_Allreduce(m, gm, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    t0 = gm[0];
    t1 = gm[1];
    totvol = gm[2];
#endif

    if (totvol <= 0.0) {
      FatalTurbulenceError("volume integral is not positive while normalizing forcing");
    }
    Real m0 = t0 / totvol;
    Real m1 = t1 / totvol;

    Real s = 0.0;
    if (normalization == TurbNormalization::edot) {
      // Solve m0*s^2 + m1*s = dedt using its non-negative root.
      if (m0 > 1.0e-30) {
        s = (-m1 + sqrt(m1 * m1 + 4.0 * m0 * dedt)) / (2.0 * m0);
      } else if (dedt > 0.0) {
        FatalTurbulenceError("cannot inject non-zero dedt with a zero forcing field");
      }
    } else {
      // Match the volume-weighted RMS acceleration independently of AMR layout.
      if (m0 > 1.0e-30) {
        s = accel_rms / sqrt(m0);
      } else if (accel_rms > 0.0) {
        FatalTurbulenceError(
            "cannot impose non-zero accel_rms with a zero forcing field");
      }
    }
    par_for(
        "force_norm", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          force_(m, 0, k, j, i) *= s;
          force_(m, 1, k, j, i) *= s;
          force_(m, 2, k, j, i) *= s;
        });
  } else {  // set force to zero
    par_for(
        "force_zero", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          force_(m, 0, k, j, i) = 0.0;
          force_(m, 1, k, j, i) = 0.0;
          force_(m, 2, k, j, i) = 0.0;
        });
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn apply forcing

//
// @brief Adds forcing in the turbulence driver.
//
// This function applies forcing in the turbulence driver based on the provided
// driver and stage. It updates the conserved variables with the applied forces
// and handles both relativistic and non-relativistic cases. Additionally, it
// supports two-fluid and magnetohydrodynamic (MHD) scenarios.
//
// @param pdrive Pointer to the driver object.
// @param stage The current stage of the driver.
// @return TaskStatus indicating the completion status of the task.
//
// The function performs the following main steps:
// 1. Applies forcing to the conserved variables using a parallel loop.
// 2. Handles relativistic transformations if required.
//

void TurbulenceDriver::ApplyForcingWithStep(Real bdt) {
  if (pmy_pack == nullptr) {
    return;
  }

  Mesh* pm = pmy_pack->pmesh;
  if (pm == nullptr) {
    return;
  }

  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  int& nx1 = indcs.nx1;
  int& nx2 = indcs.nx2;
  int& nx3 = indcs.nx3;

  Real current_time = pm->time;
  Real t_since_start = current_time - tdriv_start;

  EquationOfState* peos;

  DvceArray5D<Real> u0, u0_;
  DvceArray5D<Real> w0, w0_;
  DvceFaceFld4D<Real>* bcc0;
  if (pmy_pack->phydro != nullptr) u0 = (pmy_pack->phydro->u0);
  if (pmy_pack->phydro != nullptr) w0 = (pmy_pack->phydro->w0);
  if (pmy_pack->phydro != nullptr) peos = (pmy_pack->phydro->peos);
  if (pmy_pack->pmhd != nullptr) u0 = (pmy_pack->pmhd->u0);
  if (pmy_pack->pmhd != nullptr) w0 = (pmy_pack->pmhd->w0);
  if (pmy_pack->pmhd != nullptr) bcc0 = &(pmy_pack->pmhd->b0);
  if (pmy_pack->pmhd != nullptr) peos = pmy_pack->pmhd->peos;
  bool flag_twofl = false;
  if (pmy_pack->pionn != nullptr) {
    u0 = (pmy_pack->phydro->u0);
    u0_ = (pmy_pack->pmhd->u0);
    w0 = (pmy_pack->phydro->w0);
    w0_ = (pmy_pack->pmhd->w0);
    flag_twofl = true;
  }

  bool flag_relativistic = pmy_pack->pcoord->is_special_relativistic;

  auto force_ = force;
  const int nmkji = nmb * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;

  auto eos = peos->eos_data;  // copy-by-value (POD expected)

  if ((current_time >= tdriv_start) &&
      ((t_since_start < tdriv_duration) || turb_flag != 1)) {
    par_for(
        "push", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          Real a1 = force_(m, 0, k, j, i);
          Real a2 = force_(m, 1, k, j, i);
          Real a3 = force_(m, 2, k, j, i);

          Real den = w0(m, IDN, k, j, i);
          auto& ux = w0(m, IVX, k, j, i);
          auto& uy = w0(m, IVY, k, j, i);
          auto& uz = w0(m, IVZ, k, j, i);

          Real Fv = (a1 * ux + a2 * uy + a3 * uz);
          if (flag_relativistic) {
            // Compute Lorentz factor
            Real ut = 1. + ux * ux + uy * uy + uz * uz;
            ut = sqrt(ut);
            den /= ut;
            Fv = (a1 * ux + a2 * uy + a3 * uz) / ut;
          }
          u0(m, IM1, k, j, i) += den * a1 * bdt;
          u0(m, IM2, k, j, i) += den * a2 * bdt;
          u0(m, IM3, k, j, i) += den * a3 * bdt;
          if (eos.is_ideal) {
            u0(m, IEN, k, j, i) +=
                (Fv + 0.5 * (a1 * a1 + a2 * a2 + a3 * a3) * bdt) * den * bdt;
            // u0(m,IEN,k,j,i) += Fv*den*bdt;
          }

          if (flag_twofl) {
            den = u0_(m, IDN, k, j, i);
            u0_(m, IM1, k, j, i) += den * a1 * bdt;
            u0_(m, IM2, k, j, i) += den * a2 * bdt;
            u0_(m, IM3, k, j, i) += den * a3 * bdt;
            u0_(m, IEN, k, j, i) +=
                (Fv + 0.5 * (a1 * a1 + a2 * a2 + a3 * a3) * bdt) * den * bdt;
            // u0_(m,IEN,k,j,i) += Fv*den*bdt;
          }
        });

    // Relativistic case will require a Lorentz transformation
    if (flag_relativistic) {
      if (pmy_pack->pmhd != nullptr) {
        auto& b = *bcc0;

        par_for(
            "net_mom_4", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
            KOKKOS_LAMBDA(int m, int k, int j, int i) {
              // load single state conserved variables
              MHDCons1D u;
              u.d = u0(m, IDN, k, j, i);
              u.mx = u0(m, IM1, k, j, i);
              u.my = u0(m, IM2, k, j, i);
              u.mz = u0(m, IM3, k, j, i);
              u.e = u0(m, IEN, k, j, i);

              u.bx = 0.5 * (b.x1f(m, k, j, i) + b.x1f(m, k, j, i + 1));
              u.by = 0.5 * (b.x2f(m, k, j, i) + b.x2f(m, k, j + 1, i));
              u.bz = 0.5 * (b.x3f(m, k, j, i) + b.x3f(m, k + 1, j, i));

              // Compute (S^i S_i) (eqn C2)
              Real s2 = SQR(u.mx) + SQR(u.my) + SQR(u.mz);
              Real b2 = SQR(u.bx) + SQR(u.by) + SQR(u.bz);
              Real rpar = (u.bx * u.mx + u.by * u.my + u.bz * u.mz) / u.d;

              // call c2p function
              // (inline function in ideal_c2p_mhd.hpp file)
              HydPrim1D w;
              bool dfloor_used = false, efloor_used = false;
              // bool vceiling_used = false;
              bool c2p_failure = false;
              int iter_used = 0;
              SingleC2P_IdealSRMHD(u, eos, s2, b2, rpar, w, dfloor_used, efloor_used,
                                   c2p_failure, iter_used);
              // apply velocity ceiling if necessary
              Real lor = sqrt(1.0 + SQR(w.vx) + SQR(w.vy) + SQR(w.vz));
              if (lor > eos.gamma_max) {
                // vceiling_used = true;
                Real factor = sqrt((SQR(eos.gamma_max) - 1.0) / (SQR(lor) - 1.0));
                w.vx *= factor;
                w.vy *= factor;
                w.vz *= factor;
              }

              // Temporarily store primitives in conserved state
              u0(m, IDN, k, j, i) = w.d;
              u0(m, IM1, k, j, i) = w.vx;
              u0(m, IM2, k, j, i) = w.vy;
              u0(m, IM3, k, j, i) = w.vz;
              u0(m, IEN, k, j, i) = w.e;
            });
      } else {
        auto eos = peos->eos_data;  // copy-by-value (POD expected)

        par_for(
            "net_mom_4", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
            KOKKOS_LAMBDA(int m, int k, int j, int i) {
              u0(m, IEN, k, j, i) = fmin(u0(m, IEN, k, j, i), 40. * u0(m, IDN, k, j, i));

              // load single state conserved variables
              HydCons1D u;
              u.d = u0(m, IDN, k, j, i);
              u.mx = u0(m, IM1, k, j, i);
              u.my = u0(m, IM2, k, j, i);
              u.mz = u0(m, IM3, k, j, i);
              u.e = u0(m, IEN, k, j, i);

              // Compute (S^i S_i) (eqn C2)
              Real s2 = SQR(u.mx) + SQR(u.my) + SQR(u.mz);

              // call c2p function
              // (inline function in ideal_c2p_mhd.hpp file)
              HydPrim1D w;
              bool dfloor_used = false, efloor_used = false;
              // bool vceiling_used = false;
              bool c2p_failure = false;
              int iter_used = 0;
              SingleC2P_IdealSRHyd(u, eos, s2, w, dfloor_used, efloor_used, c2p_failure,
                                   iter_used);
              // apply velocity ceiling if necessary
              Real lor = sqrt(1.0 + SQR(w.vx) + SQR(w.vy) + SQR(w.vz));
              if (lor > eos.gamma_max) {
                // vceiling_used = true;
                Real factor = sqrt((SQR(eos.gamma_max) - 1.0) / (SQR(lor) - 1.0));
                w.vx *= factor;
                w.vy *= factor;
                w.vz *= factor;
              }

              u0(m, IDN, k, j, i) = w.d;
              u0(m, IM1, k, j, i) = w.vx;
              u0(m, IM2, k, j, i) = w.vy;
              u0(m, IM3, k, j, i) = w.vz;
              u0(m, IEN, k, j, i) = w.e;
            });
      }

      // remove net momentum
      Real t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
      Kokkos::parallel_reduce(
          "net_mom_3", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
          KOKKOS_LAMBDA(const int& idx, Real& sum_t0, Real& sum_t1, Real& sum_t2,
                        Real& sum_t3) {
            // compute n,k,j,i indices of thread
            int m = (idx) / nkji;
            int k = (idx - m * nkji) / nji;
            int j = (idx - m * nkji - k * nji) / nx1;
            int i = (idx - m * nkji - k * nji - j * nx1) + is;
            k += ks;
            j += js;

            Real u_t = sqrt(1. + u0(m, IVX, k, j, i) * u0(m, IVX, k, j, i) +
                            u0(m, IVY, k, j, i) * u0(m, IVY, k, j, i) +
                            u0(m, IVZ, k, j, i) * u0(m, IVZ, k, j, i));

            Real den = u0(m, IDN, k, j, i) * u_t;
            Real mom1 = den * u0(m, IVX, k, j, i);
            Real mom2 = den * u0(m, IVY, k, j, i);
            Real mom3 = den * u0(m, IVZ, k, j, i);

            sum_t0 += den;
            sum_t1 += mom1;
            sum_t2 += mom2;
            sum_t3 += mom3;
          },
          Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1), Kokkos::Sum<Real>(t2),
          Kokkos::Sum<Real>(t3));

#if MPI_PARALLEL_ENABLED
      Real m[4], gm[4];
      m[0] = t0;
      m[1] = t1;
      m[2] = t2;
      m[3] = t3;
      MPI_Allreduce(m, gm, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      t0 = gm[0];
      t1 = gm[1];
      t2 = gm[2];
      t3 = gm[3];
#endif

      // Compute average velocity
      Real uA_x = t1 / t0;
      Real uA_y = t2 / t0;
      Real uA_z = t3 / t0;

      Real uA_0 = sqrt(1. + uA_x * uA_x + uA_y * uA_y + uA_z * uA_z);
      Real betaA = sqrt(uA_x * uA_x + uA_y * uA_y + uA_z * uA_z) / uA_0;

      Real vx = uA_x / uA_0;
      Real vy = uA_y / uA_0;
      Real vz = uA_z / uA_0;

      if (pmy_pack->pmhd != nullptr) {
        auto b = *bcc0;             // copy handle by value
        auto eos = peos->eos_data;  // copy-by-value (POD expected)

        par_for(
            "net_mom_4", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
            KOKKOS_LAMBDA(int m, int k, int j, int i) {
              u0(m, IEN, k, j, i) = fmin(u0(m, IEN, k, j, i), 40. * u0(m, IDN, k, j, i));

              // load single state conserved variables
              MHDPrim1D u;
              u.d = u0(m, IDN, k, j, i);
              u.vx = u0(m, IM1, k, j, i);
              u.vy = u0(m, IM2, k, j, i);
              u.vz = u0(m, IM3, k, j, i);
              u.e = u0(m, IEN, k, j, i);

              u.bx = 0.5 * (b.x1f(m, k, j, i) + b.x1f(m, k, j, i + 1));
              u.by = 0.5 * (b.x2f(m, k, j, i) + b.x2f(m, k, j + 1, i));
              u.bz = 0.5 * (b.x3f(m, k, j, i) + b.x3f(m, k + 1, j, i));

              HydCons1D u_out;
              SingleP2C_IdealSRMHD(u, eos.gamma, u_out);

              Real en = u_out.d + u_out.e;
              Real sx = u_out.mx;
              Real sy = u_out.my;
              Real sz = u_out.mz;

              Real dens = u_out.d;

              auto& w = u;

              Real lorentz = sqrt(1. + w.vx * w.vx + w.vy * w.vy + w.vz * w.vz);
              Real beta = sqrt(w.vx * w.vx + w.vy * w.vy + w.vz * w.vz) / lorentz;

              u0(m, IDN, k, j, i) = dens;  // *uA_0*(1.-beta*betaA);

              // Does not require knowledge of v
              u0(m, IEN, k, j, i) = uA_0 * en - uA_0 * (sx * vx + sy * vy + sz * vz);
              u0(m, IEN, k, j, i) -= u0(m, IDN, k, j, i);

              u0(m, IM1, k, j, i) =
                  sx + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vx;
              u0(m, IM2, k, j, i) =
                  sy + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vy;
              u0(m, IM3, k, j, i) =
                  sz + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vz;

              u0(m, IM1, k, j, i) -= uA_0 * en * vx;
              u0(m, IM2, k, j, i) -= uA_0 * en * vy;
              u0(m, IM3, k, j, i) -= uA_0 * en * vz;
            });
      } else {
        auto eos = peos->eos_data;  // copy-by-value (POD expected)

        par_for(
            "net_mom_4", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
            KOKKOS_LAMBDA(int m, int k, int j, int i) {
              u0(m, IEN, k, j, i) = fmin(u0(m, IEN, k, j, i), 40. * u0(m, IDN, k, j, i));

              // load single state conserved variables
              HydPrim1D u;
              u.d = u0(m, IDN, k, j, i);
              u.vx = u0(m, IM1, k, j, i);
              u.vy = u0(m, IM2, k, j, i);
              u.vz = u0(m, IM3, k, j, i);
              u.e = u0(m, IEN, k, j, i);

              HydCons1D u_out;
              SingleP2C_IdealSRHyd(u, eos.gamma, u_out);

              Real en = u_out.d + u_out.e;
              Real sx = u_out.mx;
              Real sy = u_out.my;
              Real sz = u_out.mz;

              Real dens = u_out.d;

              auto& w = u;

              Real lorentz = sqrt(1. + w.vx * w.vx + w.vy * w.vy + w.vz * w.vz);
              Real beta = sqrt(w.vx * w.vx + w.vy * w.vy + w.vz * w.vz) / lorentz;

              u0(m, IDN, k, j, i) = dens;  //*uA_0*(1.-beta*betaA);

              // Does not require knowledge of v
              u0(m, IEN, k, j, i) = uA_0 * en - uA_0 * (sx * vx + sy * vy + sz * vz);
              u0(m, IEN, k, j, i) -= u0(m, IDN, k, j, i);
              u0(m, IM1, k, j, i) =
                  sx + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vx;
              u0(m, IM2, k, j, i) =
                  sy + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vy;
              u0(m, IM3, k, j, i) =
                  sz + (uA_0 - 1.) / (betaA * betaA) * (sx * vx + sy * vy + sz * vz) * vz;
              u0(m, IM1, k, j, i) -= uA_0 * en * vx;
              u0(m, IM2, k, j, i) -= uA_0 * en * vy;
              u0(m, IM3, k, j, i) -= uA_0 * en * vz;
            });
      }
    }  // end relativistic case
  }
  return;
}

TaskStatus TurbulenceDriver::AddForcing(Driver* pdrive, int stage) {
  if (pmy_pack == nullptr) {
    return TaskStatus::complete;
  }

  Mesh* pm = pmy_pack->pmesh;
  if (pm == nullptr) {
    return TaskStatus::complete;
  }

  Real dt = pm->dt;
  Real bdt = dt;
  if (pdrive != nullptr && stage > 0) {
    bdt = (pdrive->beta[stage - 1]) * dt;
  }

  ApplyForcingWithStep(bdt);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn EnsureBasisSize()
// \brief Detect mesh/AMR changes and ensure forcing basis and arrays match current mesh.
//        Recomputes basis for all blocks when change is detected.

TaskStatus TurbulenceDriver::EnsureBasisSize(Driver* pdrive, int stage) {
  if (pmy_pack == nullptr) {
    return TaskStatus::complete;
  }

  Mesh* pm = pmy_pack->pmesh;
  if (pm == nullptr) return TaskStatus::complete;

  // Update cached domain offsets in case AMR has modified the root-grid geometry.
  domain_x1min = pm->mesh_size.x1min;
  domain_x2min = pm->mesh_size.x2min;
  domain_x3min = pm->mesh_size.x3min;

  // --- change detection (idempotent) ---
  int nmb = pmy_pack->nmb_thispack;
  int nmb_alloc = std::max(nmb, pm->nmb_maxperrank);
  bool needs_resize = false;
  if (nmb != current_nmb_) needs_resize = true;
  if (pm->adaptive && pm->pmr != nullptr) {
    if (pm->pmr->nmb_created != last_nmb_created_ ||
        pm->pmr->nmb_deleted != last_nmb_deleted_) {
      needs_resize = true;
    }
  }
  if (!needs_resize) return TaskStatus::complete;

  // --- resize/rebuild path ---

  auto& indcs = pmy_pack->pmesh->mb_indcs;
  int ncells1 = indcs.nx1 + 2 * (indcs.ng);
  int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2 * (indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2 * (indcs.ng)) : 1;

  // Retain the AMR-capacity allocation used by the fluid state arrays.
  if (force.extent(0) < nmb_alloc) {
    Kokkos::realloc(force, nmb_alloc, 3, ncells3, ncells2, ncells1);
    Kokkos::realloc(xcos, nmb_alloc, mode_count, ncells1);
    Kokkos::realloc(xsin, nmb_alloc, mode_count, ncells1);
    Kokkos::realloc(ycos, nmb_alloc, mode_count, ncells2);
    Kokkos::realloc(ysin, nmb_alloc, mode_count, ncells2);
    Kokkos::realloc(zcos, nmb_alloc, mode_count, ncells3);
    Kokkos::realloc(zsin, nmb_alloc, mode_count, ncells3);
  }

  // Modal coefficients are not modified here. Regridding only changes geometry-
  // dependent Fourier basis values; UpdateForcing renders the same OU state.

  if (pmy_pack->pmb == nullptr) {
    return TaskStatus::complete;
  }

  BuildBasis();

  // Update tracking variables
  current_nmb_ = nmb;

  if (pm->adaptive && pm->pmr != nullptr) {
    last_nmb_created_ = pm->pmr->nmb_created;
    last_nmb_deleted_ = pm->pmr->nmb_deleted;
  }

  return TaskStatus::complete;
}

TurbulenceRestartMetadata TurbulenceDriver::RestartMetadata() const {
  TurbulenceRestartMetadata metadata{};
  metadata.version = 1;
  metadata.mode_count = mode_count;
  metadata.n_updates = n_turb_updates_yet;
  metadata.nlow = nlow;
  metadata.nhigh = nhigh;
  metadata.driving_type = driving_type;
  metadata.min_kx = min_kx;
  metadata.max_kx = max_kx;
  metadata.min_ky = min_ky;
  metadata.max_ky = max_ky;
  metadata.min_kz = min_kz;
  metadata.max_kz = max_kz;
  metadata.use_npeak = static_cast<int>(use_npeak);
  metadata.turb_flag = turb_flag;
  metadata.tile_nx = tile_nx;
  metadata.tile_ny = tile_ny;
  metadata.tile_nz = tile_nz;
  metadata.normalization = static_cast<int>(normalization);
  metadata.localization = static_cast<int>(localization);
  metadata.spectrum = static_cast<int>(spectrum);
  metadata.tcorr = tcorr;
  metadata.dt_update = dt_update;
  metadata.dedt = dedt;
  metadata.accel_rms = accel_rms;
  metadata.sol_fraction = sol_fraction;
  metadata.kpeak = kpeak;
  metadata.npeak = npeak;
  metadata.expo = expo;
  metadata.exp_prp = exp_prp;
  metadata.exp_prl = exp_prl;
  metadata.tdriv_duration = tdriv_duration;
  metadata.tdriv_start = tdriv_start;
  metadata.sigma_x1 = sigma_x1;
  metadata.sigma_x2 = sigma_x2;
  metadata.sigma_x3 = sigma_x3;
  metadata.center_x1 = center_x1;
  metadata.center_x2 = center_x2;
  metadata.center_x3 = center_x3;
  return metadata;
}

void TurbulenceDriver::ValidateRestartMetadata(
    const TurbulenceRestartMetadata& metadata) const {
  TurbulenceRestartMetadata expected = RestartMetadata();
  auto check = [](bool mismatch, const char* key) {
    if (mismatch) {
      FatalTurbulenceError("restart turbulence configuration differs for '" +
                           std::string(key) + "'");
    }
  };
  check(metadata.version != expected.version, "version");
  check(metadata.mode_count != expected.mode_count, "mode_count");
  check(metadata.nlow != expected.nlow, "nlow");
  check(metadata.nhigh != expected.nhigh, "nhigh");
  check(metadata.driving_type != expected.driving_type, "driving_type");
  check(metadata.min_kx != expected.min_kx, "min_kx");
  check(metadata.max_kx != expected.max_kx, "max_kx");
  check(metadata.min_ky != expected.min_ky, "min_ky");
  check(metadata.max_ky != expected.max_ky, "max_ky");
  check(metadata.min_kz != expected.min_kz, "min_kz");
  check(metadata.max_kz != expected.max_kz, "max_kz");
  check(metadata.use_npeak != expected.use_npeak, "npeak selection");
  check(metadata.turb_flag != expected.turb_flag, "turb_flag");
  check(metadata.tile_nx != expected.tile_nx, "tile_nx");
  check(metadata.tile_ny != expected.tile_ny, "tile_ny");
  check(metadata.tile_nz != expected.tile_nz, "tile_nz");
  check(metadata.normalization != expected.normalization, "normalization");
  check(metadata.localization != expected.localization, "localization");
  check(metadata.spectrum != expected.spectrum, "spectrum");
  check(metadata.tcorr != expected.tcorr, "tcorr");
  check(metadata.dt_update != expected.dt_update, "dt_update");
  check(metadata.dedt != expected.dedt, "dedt");
  check(metadata.accel_rms != expected.accel_rms, "accel_rms");
  check(metadata.sol_fraction != expected.sol_fraction, "sol_fraction");
  check(metadata.npeak != expected.npeak, "npeak");
  check(metadata.kpeak != expected.kpeak, "kpeak");
  check(metadata.expo != expected.expo, "expo");
  check(metadata.exp_prp != expected.exp_prp, "exp_prp");
  check(metadata.exp_prl != expected.exp_prl, "exp_prl");
  check(metadata.tdriv_duration != expected.tdriv_duration, "tdriv_duration");
  check(metadata.tdriv_start != expected.tdriv_start, "tdriv_start");
  check(metadata.sigma_x1 != expected.sigma_x1, "sigma_x1");
  check(metadata.sigma_x2 != expected.sigma_x2, "sigma_x2");
  check(metadata.sigma_x3 != expected.sigma_x3, "sigma_x3");
  check(metadata.center_x1 != expected.center_x1, "center_x1");
  check(metadata.center_x2 != expected.center_x2, "center_x2");
  check(metadata.center_x3 != expected.center_x3, "center_x3");
}
