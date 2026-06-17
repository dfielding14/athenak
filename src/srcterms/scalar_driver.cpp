//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_driver.cpp
//  \brief stochastic Ornstein-Uhlenbeck forcing for a passive scalar

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "scalar_driver.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

[[noreturn]] void FatalScalarForcingError(const std::string &message) {
  std::cout << "### FATAL ERROR in scalar forcing driver: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

bool IsCanonicalMode(int nx, int ny, int nz) {
  if (nx != 0) return nx > 0;
  if (ny != 0) return ny > 0;
  return nz > 0;
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor

ScalarForcingDriver::ScalarForcingDriver(MeshBlockPack *pp, ParameterInput *pin)
    : force("scalar_force", 1, 1, 1, 1, 1),
      mode_amp_real("scalar_mode_amp_real", 1),
      mode_amp_imag("scalar_mode_amp_imag", 1),
      mode_noise_real("scalar_mode_noise_real", 1),
      mode_noise_imag("scalar_mode_noise_imag", 1),
      mode_weight("scalar_mode_weight", 1),
      kx_mode("scalar_kx_mode", 1),
      ky_mode("scalar_ky_mode", 1),
      kz_mode("scalar_kz_mode", 1),
      xcos("scalar_xcos", 1, 1, 1),
      xsin("scalar_xsin", 1, 1, 1),
      ycos("scalar_ycos", 1, 1, 1),
      ysin("scalar_ysin", 1, 1, 1),
      zcos("scalar_zcos", 1, 1, 1),
      zsin("scalar_zsin", 1, 1, 1),
      pmy_pack(pp) {
  if (pmy_pack == nullptr || pmy_pack->pmesh == nullptr) {
    FatalScalarForcingError("driver requires a valid MeshBlockPack");
  }
  if (pmy_pack->phydro == nullptr || pmy_pack->pmhd != nullptr) {
    FatalScalarForcingError("driver currently requires single-fluid hydro");
  }
  if (!(pmy_pack->phydro->scalar_only)) {
    FatalScalarForcingError("driver requires <hydro>/scalar_only = true");
  }
  if (pmy_pack->phydro->nscalars <= 0) {
    FatalScalarForcingError("driver requires at least one passive scalar");
  }

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int nmb_alloc = std::max(nmb, pm->nmb_maxperrank);
  int ncells1 = indcs.nx1 + 2 * indcs.ng;
  int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2 * indcs.ng : 1;
  int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2 * indcs.ng : 1;
  Kokkos::realloc(force, nmb_alloc, 1, ncells3, ncells2, ncells1);

  current_nmb_ = nmb;
  if (pm->adaptive && pm->pmr != nullptr) {
    last_nmb_created_ = pm->pmr->nmb_created;
    last_nmb_deleted_ = pm->pmr->nmb_deleted;
  } else {
    last_nmb_created_ = 0;
    last_nmb_deleted_ = 0;
  }

  const std::string block_name = "scalar_driving";
  auto get_serialized_real = [pin, &block_name](const char *name, Real default_value) {
    pin->GetOrAddReal(block_name, name, default_value);
    return pin->GetReal(block_name, name);
  };

  scalar_index = pin->GetOrAddInteger(block_name, "scalar_index", 0);
  if (scalar_index < 0 || scalar_index >= pmy_pack->phydro->nscalars) {
    FatalScalarForcingError("scalar_index is outside the available passive scalars");
  }

  nlow = pin->GetOrAddInteger(block_name, "nlow", 1);
  nhigh = pin->GetOrAddInteger(block_name, "nhigh", 3);
  if (nlow < 1 || nhigh < nlow) {
    FatalScalarForcingError("nlow and nhigh must satisfy 1 <= nlow <= nhigh");
  }

  std::string spectrum_name =
      pin->GetOrAddString(block_name, "spectrum", "parabolic");
  if (spectrum_name == "parabolic") {
    spectrum = ScalarForcingSpectrum::parabolic;
  } else if (spectrum_name == "power_law") {
    spectrum = ScalarForcingSpectrum::power_law;
  } else {
    FatalScalarForcingError("spectrum must be parabolic or power_law");
  }
  npeak = get_serialized_real("npeak", 0.5 * (nlow + nhigh));
  expo = get_serialized_real("expo", 5.0 / 3.0);
  if (spectrum == ScalarForcingSpectrum::parabolic) {
    if (nhigh == nlow) {
      FatalScalarForcingError("parabolic spectrum requires nhigh greater than nlow");
    }
    if (npeak < nlow || npeak > nhigh) {
      FatalScalarForcingError("npeak must lie between nlow and nhigh");
    }
  }

  tcorr = get_serialized_real("tcorr", 0.0);
  dt_update = get_serialized_real("dt_update", 0.01);
  if (tcorr < 0.0) {
    FatalScalarForcingError("tcorr must not be negative");
  }
  if (dt_update <= 0.0) {
    FatalScalarForcingError("dt_update must be greater than zero");
  }
  rseed = pin->GetOrAddInteger(block_name, "rseed", -1);

  std::string normalization_name =
      pin->GetOrAddString(block_name, "normalization", "variance_rate");
  bool has_variance_rate = pin->DoesParameterExist(block_name, "variance_rate");
  bool has_source_rms = pin->DoesParameterExist(block_name, "source_rms");
  if (normalization_name == "variance_rate") {
    normalization = ScalarForcingNormalization::variance_rate;
    if (!has_variance_rate) {
      FatalScalarForcingError(
          "normalization = variance_rate requires variance_rate");
    }
    if (has_source_rms) {
      FatalScalarForcingError(
          "source_rms is not used with normalization = variance_rate");
    }
    variance_rate = pin->GetReal(block_name, "variance_rate");
    source_rms = 0.0;
    if (variance_rate < 0.0) {
      FatalScalarForcingError("variance_rate must not be negative");
    }
  } else if (normalization_name == "source_rms") {
    normalization = ScalarForcingNormalization::source_rms;
    if (!has_source_rms) {
      FatalScalarForcingError("normalization = source_rms requires source_rms");
    }
    if (has_variance_rate) {
      FatalScalarForcingError(
          "variance_rate is not used with normalization = source_rms");
    }
    source_rms = pin->GetReal(block_name, "source_rms");
    variance_rate = 0.0;
    if (source_rms < 0.0) {
      FatalScalarForcingError("source_rms must not be negative");
    }
  } else {
    FatalScalarForcingError(
        "normalization must be variance_rate or source_rms");
  }

  dimension = 1 + static_cast<int>(pm->multi_d) + static_cast<int>(pm->three_d);
  lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  if (lx <= 0.0 || (pm->multi_d && ly <= 0.0) || (pm->three_d && lz <= 0.0)) {
    FatalScalarForcingError("active mesh dimensions must have positive length");
  }

  const int ny_min = pm->multi_d ? -nhigh : 0;
  const int ny_max = pm->multi_d ? nhigh : 0;
  const int nz_min = pm->three_d ? -nhigh : 0;
  const int nz_max = pm->three_d ? nhigh : 0;
  const int nlow2 = nlow * nlow;
  const int nhigh2 = nhigh * nhigh;
  mode_count = 0;
  for (int nx = -nhigh; nx <= nhigh; ++nx) {
    for (int ny = ny_min; ny <= ny_max; ++ny) {
      for (int nz = nz_min; nz <= nz_max; ++nz) {
        if (!IsCanonicalMode(nx, ny, nz)) continue;
        int nsqr = nx * nx + ny * ny + nz * nz;
        if (nsqr >= nlow2 && nsqr <= nhigh2) ++mode_count;
      }
    }
  }
  if (mode_count == 0) {
    FatalScalarForcingError("selected wavenumber shell contains no modes");
  }

  Kokkos::realloc(mode_amp_real, mode_count);
  Kokkos::realloc(mode_amp_imag, mode_count);
  Kokkos::realloc(mode_noise_real, mode_count);
  Kokkos::realloc(mode_noise_imag, mode_count);
  Kokkos::realloc(mode_weight, mode_count);
  Kokkos::realloc(kx_mode, mode_count);
  Kokkos::realloc(ky_mode, mode_count);
  Kokkos::realloc(kz_mode, mode_count);
  Kokkos::realloc(xcos, nmb_alloc, mode_count, ncells1);
  Kokkos::realloc(xsin, nmb_alloc, mode_count, ncells1);
  Kokkos::realloc(ycos, nmb_alloc, mode_count, ncells2);
  Kokkos::realloc(ysin, nmb_alloc, mode_count, ncells2);
  Kokkos::realloc(zcos, nmb_alloc, mode_count, ncells3);
  Kokkos::realloc(zsin, nmb_alloc, mode_count, ncells3);

  n_updates_yet = 0;
  Initialize();

  if (global_variable::my_rank == 0) {
    std::cout << "Initialising scalar forcing module" << std::endl
              << " scalar_index = " << scalar_index
              << " modes = " << mode_count
              << " tcorr = " << tcorr
              << " dt_update = " << dt_update << std::endl;
  }
}

ScalarForcingDriver::~ScalarForcingDriver() {}

//----------------------------------------------------------------------------------------
//! \brief initialize modal data, RNG state, and geometry-dependent basis

void ScalarForcingDriver::Initialize() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int ncells1 = indcs.nx1 + 2 * indcs.ng;
  int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2 * indcs.ng : 1;
  int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2 * indcs.ng : 1;
  auto force_ = force;
  par_for(
      "scalar_force_init", DevExeSpace(), 0, nmb - 1, 0, 0, 0, ncells3 - 1,
      0, ncells2 - 1, 0, ncells1 - 1,
      KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
        force_(m, n, k, j, i) = 0.0;
      });

  if (rseed >= 0) {
    int seed = (rseed > 0) ? rseed : 1;
    rstate.idum = -static_cast<decltype(rstate.idum)>(seed);
  } else {
    rstate.idum = -1;
  }
  rstate.iset = 0;

  const Real dkx = 2.0 * M_PI / lx;
  const Real dky = pmy_pack->pmesh->multi_d ? 2.0 * M_PI / ly : 0.0;
  const Real dkz = pmy_pack->pmesh->three_d ? 2.0 * M_PI / lz : 0.0;
  const int ny_min = pmy_pack->pmesh->multi_d ? -nhigh : 0;
  const int ny_max = pmy_pack->pmesh->multi_d ? nhigh : 0;
  const int nz_min = pmy_pack->pmesh->three_d ? -nhigh : 0;
  const int nz_max = pmy_pack->pmesh->three_d ? nhigh : 0;
  const int nlow2 = nlow * nlow;
  const int nhigh2 = nhigh * nhigh;

  int nmode = 0;
  for (int nx = -nhigh; nx <= nhigh; ++nx) {
    for (int ny = ny_min; ny <= ny_max; ++ny) {
      for (int nz = nz_min; nz <= nz_max; ++nz) {
        if (!IsCanonicalMode(nx, ny, nz)) continue;
        int nsqr = nx * nx + ny * ny + nz * nz;
        if (nsqr < nlow2 || nsqr > nhigh2) continue;

        Real nmag = std::sqrt(static_cast<Real>(nsqr));
        Real weight = 0.0;
        if (spectrum == ScalarForcingSpectrum::power_law) {
          weight = std::pow(nmag, -0.5 * (expo + dimension - 1.0));
        } else {
          Real width = static_cast<Real>(nhigh - nlow);
          Real shell_power =
              std::fabs(1.0 - 4.0 * SQR((nmag - npeak) / width));
          weight = std::sqrt(shell_power) *
                   std::pow(npeak / nmag, 0.5 * (dimension - 1.0));
        }

        mode_amp_real.h_view(nmode) = 0.0;
        mode_amp_imag.h_view(nmode) = 0.0;
        mode_noise_real.h_view(nmode) = 0.0;
        mode_noise_imag.h_view(nmode) = 0.0;
        mode_weight.h_view(nmode) = weight;
        kx_mode.h_view(nmode) = dkx * nx;
        ky_mode.h_view(nmode) = dky * ny;
        kz_mode.h_view(nmode) = dkz * nz;
        ++nmode;
      }
    }
  }

  mode_amp_real.template modify<HostMemSpace>();
  mode_amp_real.template sync<DevExeSpace>();
  mode_amp_imag.template modify<HostMemSpace>();
  mode_amp_imag.template sync<DevExeSpace>();
  mode_weight.template modify<HostMemSpace>();
  mode_weight.template sync<DevExeSpace>();
  kx_mode.template modify<HostMemSpace>();
  kx_mode.template sync<DevExeSpace>();
  ky_mode.template modify<HostMemSpace>();
  ky_mode.template sync<DevExeSpace>();
  kz_mode.template modify<HostMemSpace>();
  kz_mode.template sync<DevExeSpace>();

  BuildBasis();
}

//----------------------------------------------------------------------------------------
//! \brief build separable trigonometric basis arrays on active cells

void ScalarForcingDriver::BuildBasis() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;

  auto size = pmy_pack->pmb->mb_size;
  size.template modify<HostMemSpace>();
  size.template sync<DevExeSpace>();
  auto kx_mode_ = kx_mode;
  auto ky_mode_ = ky_mode;
  auto kz_mode_ = kz_mode;
  auto xcos_ = xcos;
  auto xsin_ = xsin;
  auto ycos_ = ycos;
  auto ysin_ = ysin;
  auto zcos_ = zcos;
  auto zsin_ = zsin;

  par_for(
      "scalar_basis_x", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, is, ie,
      KOKKOS_LAMBDA(int m, int n, int i) {
        Real x = CellCenterX(i - is, nx1, size.d_view(m).x1min,
                             size.d_view(m).x1max);
        xcos_(m, n, i) = cos(kx_mode_.d_view(n) * x);
        xsin_(m, n, i) = sin(kx_mode_.d_view(n) * x);
      });
  par_for(
      "scalar_basis_y", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, js, je,
      KOKKOS_LAMBDA(int m, int n, int j) {
        Real y = CellCenterX(j - js, nx2, size.d_view(m).x2min,
                             size.d_view(m).x2max);
        ycos_(m, n, j) = cos(ky_mode_.d_view(n) * y);
        ysin_(m, n, j) = sin(ky_mode_.d_view(n) * y);
      });
  par_for(
      "scalar_basis_z", DevExeSpace(), 0, nmb - 1, 0, mode_count - 1, ks, ke,
      KOKKOS_LAMBDA(int m, int n, int k) {
        Real z = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                             size.d_view(m).x3max);
        zcos_(m, n, k) = cos(kz_mode_.d_view(n) * z);
        zsin_(m, n, k) = sin(kz_mode_.d_view(n) * z);
      });
}

//----------------------------------------------------------------------------------------
//! \brief add tasks that update modal state and render the scalar source

void ScalarForcingDriver::IncludeInitializeModesTask(std::shared_ptr<TaskList> tl,
                                                     TaskID start) {
  auto id_resize = tl->AddTask(&ScalarForcingDriver::EnsureBasisSize, this, start);
  auto id_modes = tl->AddTask(&ScalarForcingDriver::InitializeModes, this, id_resize);
  (void) tl->AddTask(&ScalarForcingDriver::UpdateForcing, this, id_modes);
}

//----------------------------------------------------------------------------------------
//! \brief insert the scalar source after the conservative RK update

void ScalarForcingDriver::IncludeAddForcingTask(std::shared_ptr<TaskList> tl,
                                                TaskID start) {
  (void) start;
  (void) tl->InsertTask(&ScalarForcingDriver::AddForcing, this,
                        pmy_pack->phydro->id.rkupdt,
                        pmy_pack->phydro->id.srctrms);
}

//----------------------------------------------------------------------------------------
//! \brief evolve modal OU coefficients to the update index required at current time

TaskStatus ScalarForcingDriver::InitializeModes(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  Real current_time = std::max(pmy_pack->pmesh->time, static_cast<Real>(0.0));
  int updates_required = static_cast<int>(current_time / dt_update) + 1;

  for (int update = n_updates_yet; update < updates_required; ++update) {
    for (int n = 0; n < mode_count; ++n) {
      Real weight = mode_weight.h_view(n);
      mode_noise_real.h_view(n) = weight * RanGaussianSt(&rstate);
      mode_noise_imag.h_view(n) = weight * RanGaussianSt(&rstate);
    }

    Real fcorr = 0.0;
    Real gcorr = 1.0;
    if (tcorr > 1.0e-6 && update > 0) {
      fcorr = std::exp(-dt_update / tcorr);
      gcorr = std::sqrt(1.0 - fcorr * fcorr);
    }
    for (int n = 0; n < mode_count; ++n) {
      mode_amp_real.h_view(n) =
          fcorr * mode_amp_real.h_view(n) + gcorr * mode_noise_real.h_view(n);
      mode_amp_imag.h_view(n) =
          fcorr * mode_amp_imag.h_view(n) + gcorr * mode_noise_imag.h_view(n);
    }
  }

  if (updates_required > n_updates_yet) {
    mode_amp_real.template modify<HostMemSpace>();
    mode_amp_real.template sync<DevExeSpace>();
    mode_amp_imag.template modify<HostMemSpace>();
    mode_amp_imag.template sync<DevExeSpace>();
  }
  n_updates_yet = updates_required;
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \brief render the authoritative modal state on the current mesh

void ScalarForcingDriver::RenderForce() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  auto force_ = force;
  auto amp_real = mode_amp_real;
  auto amp_imag = mode_amp_imag;
  auto xcos_ = xcos;
  auto xsin_ = xsin;
  auto ycos_ = ycos;
  auto ysin_ = ysin;
  auto zcos_ = zcos;
  auto zsin_ = zsin;
  const int mode_count_ = mode_count;

  mode_amp_real.template sync<DevExeSpace>();
  mode_amp_imag.template sync<DevExeSpace>();
  par_for(
      "scalar_force_render", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real value = 0.0;
        for (int n = 0; n < mode_count_; ++n) {
          Real phase_real =
              (xcos_(m, n, i) * ycos_(m, n, j) -
               xsin_(m, n, i) * ysin_(m, n, j)) *
                  zcos_(m, n, k) -
              (xsin_(m, n, i) * ycos_(m, n, j) +
               xcos_(m, n, i) * ysin_(m, n, j)) *
                  zsin_(m, n, k);
          Real phase_imag =
              (ycos_(m, n, j) * zsin_(m, n, k) +
               ysin_(m, n, j) * zcos_(m, n, k)) *
                  xcos_(m, n, i) +
              (ycos_(m, n, j) * zcos_(m, n, k) -
               ysin_(m, n, j) * zsin_(m, n, k)) *
                  xsin_(m, n, i);
          value += amp_real.d_view(n) * phase_real -
                   amp_imag.d_view(n) * phase_imag;
        }
        force_(m, 0, k, j, i) = value;
      });
}

//----------------------------------------------------------------------------------------
//! \brief remove the mass-weighted mean and normalize the rendered source

TaskStatus ScalarForcingDriver::UpdateForcing(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  RenderForce();

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nmkji = nmb * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;
  const int scalar_component = pmy_pack->phydro->nhydro + scalar_index;
  auto u0 = pmy_pack->phydro->u0;
  auto force_ = force;
  auto size = pmy_pack->pmb->mb_size;
  size.template modify<HostMemSpace>();
  size.template sync<DevExeSpace>();

  Real mass = 0.0;
  Real source_integral = 0.0;
  Kokkos::parallel_reduce(
      "scalar_force_mean", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &sum_mass, Real &sum_source) {
        int m = idx / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        int i = (idx - m * nkji - k * nji - j * nx1) + is;
        k += ks;
        j += js;
        Real vol =
            size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
        Real rho = u0(m, IDN, k, j, i);
        sum_mass += rho * vol;
        sum_source += rho * force_(m, 0, k, j, i) * vol;
      },
      Kokkos::Sum<Real>(mass), Kokkos::Sum<Real>(source_integral));

#if MPI_PARALLEL_ENABLED
  Real local_mean_data[2] = {mass, source_integral};
  Real global_mean_data[2] = {0.0, 0.0};
  MPI_Allreduce(local_mean_data, global_mean_data, 2, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  mass = global_mean_data[0];
  source_integral = global_mean_data[1];
#endif

  if (mass <= 0.0) {
    FatalScalarForcingError("mass integral is not positive");
  }
  Real source_mean = source_integral / mass;
  par_for(
      "scalar_force_project_mean", DevExeSpace(), 0, nmb - 1, ks,
      ks + nx3 - 1, js, js + nx2 - 1, is, is + nx1 - 1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        force_(m, 0, k, j, i) -= source_mean;
      });

  Real scalar_integral = 0.0;
  Real projected_source_integral = 0.0;
  Real scalar_source_integral = 0.0;
  Real source_square_integral = 0.0;
  Kokkos::parallel_reduce(
      "scalar_force_norm", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &sum_scalar, Real &sum_source,
                    Real &sum_scalar_source, Real &sum_source_square) {
        int m = idx / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        int i = (idx - m * nkji - k * nji - j * nx1) + is;
        k += ks;
        j += js;
        Real vol =
            size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
        Real rho = u0(m, IDN, k, j, i);
        Real scalar = u0(m, scalar_component, k, j, i);
        Real source = force_(m, 0, k, j, i);
        sum_scalar += scalar * vol;
        sum_source += rho * source * vol;
        sum_scalar_source += scalar * source * vol;
        sum_source_square += rho * source * source * vol;
      },
      Kokkos::Sum<Real>(scalar_integral),
      Kokkos::Sum<Real>(projected_source_integral),
      Kokkos::Sum<Real>(scalar_source_integral),
      Kokkos::Sum<Real>(source_square_integral));

#if MPI_PARALLEL_ENABLED
  Real local_norm_data[4] = {scalar_integral, projected_source_integral,
                             scalar_source_integral, source_square_integral};
  Real global_norm_data[4] = {0.0, 0.0, 0.0, 0.0};
  MPI_Allreduce(local_norm_data, global_norm_data, 4, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  scalar_integral = global_norm_data[0];
  projected_source_integral = global_norm_data[1];
  scalar_source_integral = global_norm_data[2];
  source_square_integral = global_norm_data[3];
#endif

  Real mean_scalar = scalar_integral / mass;
  Real mean_source = projected_source_integral / mass;
  Real covariance =
      scalar_source_integral / mass - mean_scalar * mean_source;
  Real mean_source_square = source_square_integral / mass;
  Real scale = 0.0;

  if (normalization == ScalarForcingNormalization::source_rms) {
    if (mean_source_square > 1.0e-30) {
      scale = source_rms / std::sqrt(mean_source_square);
    } else if (source_rms > 0.0) {
      FatalScalarForcingError(
          "cannot impose non-zero source_rms with a zero forcing field");
    }
  } else {
    Real dt = pmy_pack->pmesh->dt;
    if (dt <= 0.0) {
      FatalScalarForcingError(
          "variance_rate normalization requires a positive timestep");
    }
    if (variance_rate > 0.0) {
      if (mean_source_square <= 1.0e-30) {
        FatalScalarForcingError(
            "cannot inject non-zero variance_rate with a zero forcing field");
      }
      Real discriminant =
          std::sqrt(covariance * covariance +
                    2.0 * dt * mean_source_square * variance_rate);
      if (covariance >= 0.0) {
        scale = 2.0 * variance_rate / (covariance + discriminant);
      } else {
        scale = (-covariance + discriminant) /
                (dt * mean_source_square);
      }
    }
  }

  par_for(
      "scalar_force_scale", DevExeSpace(), 0, nmb - 1, ks, ks + nx3 - 1,
      js, js + nx2 - 1, is, is + nx1 - 1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        force_(m, 0, k, j, i) *= scale;
      });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \brief add the scalar source for one RK stage

void ScalarForcingDriver::ApplyForcingWithStep(Real bdt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int scalar_component = pmy_pack->phydro->nhydro + scalar_index;
  auto u0 = pmy_pack->phydro->u0;
  auto force_ = force;

  par_for(
      "scalar_force_apply", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real rho = u0(m, IDN, k, j, i);
        u0(m, scalar_component, k, j, i) +=
            rho * force_(m, 0, k, j, i) * bdt;
      });
}

TaskStatus ScalarForcingDriver::AddForcing(Driver *pdrive, int stage) {
  Real bdt = pmy_pack->pmesh->dt;
  if (pdrive != nullptr && stage > 0) {
    bdt = pdrive->beta[stage - 1] * pmy_pack->pmesh->dt;
  }
  ApplyForcingWithStep(bdt);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \brief rebuild geometry-dependent arrays after AMR topology changes

TaskStatus ScalarForcingDriver::EnsureBasisSize(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  Mesh *pm = pmy_pack->pmesh;
  int nmb = pmy_pack->nmb_thispack;
  bool mesh_changed = (nmb != current_nmb_);
  if (pm->adaptive && pm->pmr != nullptr) {
    mesh_changed = mesh_changed ||
                   pm->pmr->nmb_created != last_nmb_created_ ||
                   pm->pmr->nmb_deleted != last_nmb_deleted_;
  }
  if (!mesh_changed) return TaskStatus::complete;

  auto &indcs = pm->mb_indcs;
  int nmb_alloc = std::max(nmb, pm->nmb_maxperrank);
  int ncells1 = indcs.nx1 + 2 * indcs.ng;
  int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2 * indcs.ng : 1;
  int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2 * indcs.ng : 1;
  if (force.extent(0) < static_cast<std::size_t>(nmb_alloc)) {
    Kokkos::realloc(force, nmb_alloc, 1, ncells3, ncells2, ncells1);
    Kokkos::realloc(xcos, nmb_alloc, mode_count, ncells1);
    Kokkos::realloc(xsin, nmb_alloc, mode_count, ncells1);
    Kokkos::realloc(ycos, nmb_alloc, mode_count, ncells2);
    Kokkos::realloc(ysin, nmb_alloc, mode_count, ncells2);
    Kokkos::realloc(zcos, nmb_alloc, mode_count, ncells3);
    Kokkos::realloc(zsin, nmb_alloc, mode_count, ncells3);
  }
  BuildBasis();
  current_nmb_ = nmb;
  if (pm->adaptive && pm->pmr != nullptr) {
    last_nmb_created_ = pm->pmr->nmb_created;
    last_nmb_deleted_ = pm->pmr->nmb_deleted;
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \brief return the configuration and dynamic update count stored in restart files

ScalarForcingRestartMetadata ScalarForcingDriver::RestartMetadata() const {
  ScalarForcingRestartMetadata metadata{};
  metadata.version = 1;
  metadata.mode_count = mode_count;
  metadata.n_updates = n_updates_yet;
  metadata.scalar_index = scalar_index;
  metadata.nlow = nlow;
  metadata.nhigh = nhigh;
  metadata.normalization = static_cast<int>(normalization);
  metadata.spectrum = static_cast<int>(spectrum);
  metadata.tcorr = tcorr;
  metadata.dt_update = dt_update;
  metadata.variance_rate = variance_rate;
  metadata.source_rms = source_rms;
  metadata.npeak = npeak;
  metadata.expo = expo;
  return metadata;
}

void ScalarForcingDriver::ValidateRestartMetadata(
    const ScalarForcingRestartMetadata &metadata) const {
  ScalarForcingRestartMetadata expected = RestartMetadata();
  auto check = [](bool mismatch, const char *key) {
    if (mismatch) {
      FatalScalarForcingError(
          "restart scalar-forcing configuration differs for '" +
          std::string(key) + "'");
    }
  };
  check(metadata.version != expected.version, "version");
  check(metadata.mode_count != expected.mode_count, "mode_count");
  check(metadata.scalar_index != expected.scalar_index, "scalar_index");
  check(metadata.nlow != expected.nlow, "nlow");
  check(metadata.nhigh != expected.nhigh, "nhigh");
  check(metadata.normalization != expected.normalization, "normalization");
  check(metadata.spectrum != expected.spectrum, "spectrum");
  check(metadata.tcorr != expected.tcorr, "tcorr");
  check(metadata.dt_update != expected.dt_update, "dt_update");
  check(metadata.variance_rate != expected.variance_rate, "variance_rate");
  check(metadata.source_rms != expected.source_rms, "source_rms");
  check(metadata.npeak != expected.npeak, "npeak");
  check(metadata.expo != expected.expo, "expo");
}
