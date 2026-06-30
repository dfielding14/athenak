//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file initial_perturbations.cpp
//! \brief One-time Fourier perturbations applied to initial conditions.

#include "initial_perturbations.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

bool ListContains(std::string list, const std::string &needle) {
  list = Lower(list);
  for (char &c : list) {
    if (c == ',' || c == ';') c = ' ';
  }
  std::stringstream ss(list);
  std::string item;
  while (ss >> item) {
    if (item == needle) return true;
    if (needle == "density" && item == "rho") return true;
    if (needle == "magnetic" &&
        (item == "b" || item == "field" || item == "magnetic_field")) {
      return true;
    }
    if (needle == "velocity" && (item == "v" || item == "vel")) return true;
  }
  return false;
}

bool IsCanonicalMode(const int nkx, const int nky, const int nkz) {
  if (nkx != 0) return nkx > 0;
  if (nky != 0) return nky > 0;
  return nkz > 0;
}

Real GetRealAny(ParameterInput *pin, const std::string &block,
                const std::string &primary, const Real def,
                std::initializer_list<const char *> aliases = {}) {
  if (pin->DoesParameterExist(block, primary)) return pin->GetReal(block, primary);
  for (const char *alias : aliases) {
    if (pin->DoesParameterExist(block, alias)) return pin->GetReal(block, alias);
  }
  return pin->GetOrAddReal(block, primary, def);
}

int GetIntegerAny(ParameterInput *pin, const std::string &block,
                  const std::string &primary, const int def,
                  std::initializer_list<const char *> aliases = {}) {
  if (pin->DoesParameterExist(block, primary)) return pin->GetInteger(block, primary);
  for (const char *alias : aliases) {
    if (pin->DoesParameterExist(block, alias)) return pin->GetInteger(block, alias);
  }
  return pin->GetOrAddInteger(block, primary, def);
}

bool GetBooleanAny(ParameterInput *pin, const std::string &block,
                   const std::string &primary, const bool def,
                   std::initializer_list<const char *> aliases = {}) {
  if (pin->DoesParameterExist(block, primary)) return pin->GetBoolean(block, primary);
  for (const char *alias : aliases) {
    if (pin->DoesParameterExist(block, alias)) return pin->GetBoolean(block, alias);
  }
  return pin->GetOrAddBoolean(block, primary, def);
}

void FatalInitialPerturbationError(const std::string &msg) {
  std::cout << "### FATAL ERROR in initial perturbations: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void ReduceRealSum(Real *vals, const int n) {
#if MPI_PARALLEL_ENABLED
  Real global_vals[8];
  if (n > 8) {
    FatalInitialPerturbationError("internal MPI reduction buffer is too small");
  }
  MPI_Allreduce(vals, global_vals, n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int i = 0; i < n; ++i) vals[i] = global_vals[i];
#else
  (void) vals;
  (void) n;
#endif
}

KOKKOS_INLINE_FUNCTION
Real PerturbationWindow(const Real x1, const Real x2, const Real x3,
                        const int mode, const Real c1, const Real c2, const Real c3,
                        const Real s1, const Real s2, const Real s3,
                        const bool multi_d, const bool three_d) {
  if (mode == 0) return 1.0;

  bool active = false;
  Real arg = 0.0;
  if (s1 > 0.0) {
    arg += SQR((x1 - c1)/s1);
    active = true;
  }
  if (multi_d && s2 > 0.0) {
    arg += SQR((x2 - c2)/s2);
    active = true;
  }
  if (three_d && s3 > 0.0) {
    arg += SQR((x3 - c3)/s3);
    active = true;
  }
  if (!active) return 1.0;

  const Real w = exp(-0.5*arg);
  return (mode == 2) ? (1.0 - w) : w;
}

template <typename KView, typename AmpView>
KOKKOS_INLINE_FUNCTION
Real FourierComponent(const KView &kx, const KView &ky, const KView &kz,
                      const AmpView &amp_real, const AmpView &amp_imag,
                      const int comp, const int mode_count,
                      const Real x1, const Real x2, const Real x3,
                      const Real x1_origin, const Real x2_origin, const Real x3_origin) {
  Real val = 0.0;
  const Real xx1 = x1 - x1_origin;
  const Real xx2 = x2 - x2_origin;
  const Real xx3 = x3 - x3_origin;
  for (int n = 0; n < mode_count; ++n) {
    const Real phase = kx(n)*xx1 + ky(n)*xx2 + kz(n)*xx3;
    val += amp_real(comp,n)*cos(phase) - amp_imag(comp,n)*sin(phase);
  }
  return val;
}

}  // namespace

InitialPerturbations::InitialPerturbations(MeshBlockPack *pp, ParameterInput *pin,
                                           const std::string &block_name) :
    pmy_pack(pp),
    block_name_(block_name),
    perturb_density(false),
    perturb_velocity(false),
    perturb_magnetic(false),
    density_fractional(true),
    remove_density_mean(false),
    remove_velocity_mean(false),
    density_rms(0.0),
    velocity_rms(0.0),
    magnetic_rms(0.0),
    nlow(1),
    nhigh(3),
    min_kx(0),
    max_kx(3),
    min_ky(0),
    max_ky(3),
    min_kz(0),
    max_kz(3),
    mode_count(0),
    spectral_slope(5.0/3.0),
    f_solenoidal(1.0),
    rseed(-1),
    localization_mode(0),
    x1_center(0.0),
    x2_center(0.0),
    x3_center(0.0),
    x1_scale(-1.0),
    x2_scale(-1.0),
    x3_scale(-1.0),
    domain_x1min(0.0),
    domain_x2min(0.0),
    domain_x3min(0.0),
    lx(1.0),
    ly(1.0),
    lz(1.0),
    kx_mode("initial_pert_kx", 1),
    ky_mode("initial_pert_ky", 1),
    kz_mode("initial_pert_kz", 1),
    rho_amp_real("initial_pert_rho_amp_real", 1, 1),
    rho_amp_imag("initial_pert_rho_amp_imag", 1, 1),
    vel_amp_real("initial_pert_vel_amp_real", 1, 1),
    vel_amp_imag("initial_pert_vel_amp_imag", 1, 1),
    avec_amp_real("initial_pert_avec_amp_real", 1, 1),
    avec_amp_imag("initial_pert_avec_amp_imag", 1, 1),
    rho_pert("initial_pert_rho", 1, 1, 1, 1),
    vel_pert("initial_pert_vel", 1, 1, 1, 1, 1),
    avec("initial_pert_avec", 1, 1, 1, 1),
    b_pert("initial_pert_b", 1, 1, 1, 1) {
  if (pmy_pack == nullptr || pin == nullptr) {
    return;
  }

  Mesh *pm = pmy_pack->pmesh;
  domain_x1min = pm->mesh_size.x1min;
  domain_x2min = pm->mesh_size.x2min;
  domain_x3min = pm->mesh_size.x3min;
  lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  lz = pm->mesh_size.x3max - pm->mesh_size.x3min;

  const Real x1_mid = 0.5*(pm->mesh_size.x1min + pm->mesh_size.x1max);
  const Real x2_mid = 0.5*(pm->mesh_size.x2min + pm->mesh_size.x2max);
  const Real x3_mid = 0.5*(pm->mesh_size.x3min + pm->mesh_size.x3max);

  std::string variables = pin->GetOrAddString(block_name_, "variables", "");
  perturb_density = GetBooleanAny(pin, block_name_, "perturb_density", false,
                                  {"density", "perturb_rho"});
  perturb_velocity = GetBooleanAny(pin, block_name_, "perturb_velocity", false,
                                   {"velocity", "perturb_v"});
  perturb_magnetic = GetBooleanAny(pin, block_name_, "perturb_magnetic", false,
                                   {"magnetic", "perturb_b", "perturb_B"});

  perturb_density = perturb_density || ListContains(variables, "density");
  perturb_velocity = perturb_velocity || ListContains(variables, "velocity");
  perturb_magnetic = perturb_magnetic || ListContains(variables, "magnetic");

  density_rms = GetRealAny(pin, block_name_, "density_rms", 0.0,
                           {"rho_rms", "drho_rms"});
  velocity_rms = GetRealAny(pin, block_name_, "velocity_rms", 0.0,
                            {"v_rms", "vel_rms"});
  magnetic_rms = GetRealAny(pin, block_name_, "magnetic_rms", 0.0,
                            {"b_rms", "B_rms"});
  if (density_rms > 0.0) perturb_density = true;
  if (velocity_rms > 0.0) perturb_velocity = true;
  if (magnetic_rms > 0.0) perturb_magnetic = true;

  density_fractional = GetBooleanAny(pin, block_name_, "density_fractional", true,
                                     {"rho_fractional"});
  remove_density_mean = GetBooleanAny(pin, block_name_, "remove_density_mean", false,
                                      {"rho_remove_mean"});
  remove_velocity_mean = GetBooleanAny(pin, block_name_, "remove_velocity_mean", false,
                                       {"vel_remove_mean"});

  nlow = GetIntegerAny(pin, block_name_, "nlow", 1, {"kmin"});
  nhigh = GetIntegerAny(pin, block_name_, "nhigh", 3, {"kmax"});
  min_kx = GetIntegerAny(pin, block_name_, "min_kx", -nhigh);
  max_kx = GetIntegerAny(pin, block_name_, "max_kx", nhigh);
  min_ky = GetIntegerAny(pin, block_name_, "min_ky", -nhigh);
  max_ky = GetIntegerAny(pin, block_name_, "max_ky", nhigh);
  min_kz = GetIntegerAny(pin, block_name_, "min_kz", -nhigh);
  max_kz = GetIntegerAny(pin, block_name_, "max_kz", nhigh);
  spectral_slope = GetRealAny(pin, block_name_, "spectral_slope", 5.0/3.0,
                              {"slope", "expo"});
  f_solenoidal = GetRealAny(pin, block_name_, "f_solenoidal", 1.0,
                            {"sol_fraction"});
  rseed = GetIntegerAny(pin, block_name_, "rseed", -1);

  std::string loc = Lower(pin->GetOrAddString(block_name_, "localization", "none"));
  if (loc == "include" || loc == "region" || loc == "add" || loc == "only") {
    localization_mode = 1;
  } else if (loc == "exclude" || loc == "avoid" || loc == "outside") {
    localization_mode = 2;
  } else if (loc == "none" || loc == "off" || loc == "false") {
    localization_mode = 0;
  } else {
    FatalInitialPerturbationError("localization must be none, include, or exclude");
  }

  x1_center = GetRealAny(pin, block_name_, "x1_center", x1_mid,
                         {"x_center", "x_turb_center"});
  x2_center = GetRealAny(pin, block_name_, "x2_center", x2_mid,
                         {"y_center", "y_turb_center"});
  x3_center = GetRealAny(pin, block_name_, "x3_center", x3_mid,
                         {"z_center", "z_turb_center"});
  x1_scale = GetRealAny(pin, block_name_, "x1_scale", -1.0,
                        {"x_scale", "x_turb_scale_height"});
  x2_scale = GetRealAny(pin, block_name_, "x2_scale", -1.0,
                        {"y_scale", "y_turb_scale_height"});
  x3_scale = GetRealAny(pin, block_name_, "x3_scale", -1.0,
                        {"z_scale", "z_turb_scale_height"});
  if (localization_mode == 0 && (x1_scale > 0.0 || x2_scale > 0.0 || x3_scale > 0.0)) {
    localization_mode = 1;
  }

  if (pm->one_d) {
    min_ky = 0;
    max_ky = 0;
    min_kz = 0;
    max_kz = 0;
  } else if (pm->two_d) {
    min_kz = 0;
    max_kz = 0;
  }

  Validate();
  if (!AnyEnabled()) {
    return;
  }
  InitializeModes();

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;

  if (DensityEnabled()) {
    Kokkos::realloc(rho_pert, nmb, ncells3, ncells2, ncells1);
  }
  if (VelocityEnabled()) {
    Kokkos::realloc(vel_pert, nmb, 3, ncells3, ncells2, ncells1);
  }
  if (MagneticEnabled()) {
    Kokkos::realloc(avec.x1e, nmb, ncells3 + 1, ncells2 + 1, ncells1);
    Kokkos::realloc(avec.x2e, nmb, ncells3 + 1, ncells2, ncells1 + 1);
    Kokkos::realloc(avec.x3e, nmb, ncells3, ncells2 + 1, ncells1 + 1);
    Kokkos::realloc(b_pert.x1f, nmb, ncells3, ncells2, ncells1 + 1);
    Kokkos::realloc(b_pert.x2f, nmb, ncells3, ncells2 + 1, ncells1);
    Kokkos::realloc(b_pert.x3f, nmb, ncells3 + 1, ncells2, ncells1);
  }

  const int seed = (rseed > 0) ? rseed : 1;
  rstate.idum = -static_cast<int64_t>(seed);
  rstate.iset = 0;
}

void InitialPerturbations::Validate() const {
  if (!std::isfinite(density_rms) || !std::isfinite(velocity_rms) ||
      !std::isfinite(magnetic_rms)) {
    FatalInitialPerturbationError("rms amplitudes must be finite");
  }
  if (density_rms < 0.0 || velocity_rms < 0.0 || magnetic_rms < 0.0) {
    FatalInitialPerturbationError("rms amplitudes must be non-negative");
  }
  if (!AnyEnabled()) return;
  if (pmy_pack->pcoord->is_special_relativistic ||
      pmy_pack->pcoord->is_general_relativistic) {
    FatalInitialPerturbationError("non-relativistic hydro/MHD is required");
  }
  if (pmy_pack->phydro == nullptr && pmy_pack->pmhd == nullptr) {
    FatalInitialPerturbationError("requires an active <hydro> or <mhd> block");
  }
  if (pmy_pack->pionn != nullptr) {
    FatalInitialPerturbationError("ion-neutral initial perturbations are not "
                                  "implemented");
  }
  if (MagneticEnabled() && pmy_pack->pmhd == nullptr) {
    FatalInitialPerturbationError("magnetic perturbations require an active <mhd> block");
  }
  if (nlow < 0 || nhigh < nlow) {
    FatalInitialPerturbationError("require 0 <= nlow <= nhigh");
  }
  if (min_kx > max_kx || min_ky > max_ky || min_kz > max_kz) {
    FatalInitialPerturbationError("each min_k* must be <= the corresponding max_k*");
  }
  if (!std::isfinite(spectral_slope)) {
    FatalInitialPerturbationError("spectral_slope must be finite");
  }
  if (VelocityEnabled() && (!std::isfinite(f_solenoidal) ||
                           f_solenoidal < 0.0 || f_solenoidal > 1.0)) {
    FatalInitialPerturbationError("f_solenoidal must lie in [0, 1]");
  }
  if (localization_mode != 0) {
    const bool active_scale = x1_scale > 0.0 ||
        (pmy_pack->pmesh->multi_d && x2_scale > 0.0) ||
        (pmy_pack->pmesh->three_d && x3_scale > 0.0);
    if (!active_scale) {
      FatalInitialPerturbationError("localization requires a positive scale in an active "
                                    "mesh dimension");
    }
  }
}

void InitialPerturbations::InitializeModes() {
  if (!AnyEnabled()) return;

  const Real dkx = (lx > 0.0) ? 2.0*M_PI/lx : 0.0;
  const Real dky = (ly > 0.0) ? 2.0*M_PI/ly : 0.0;
  const Real dkz = (lz > 0.0) ? 2.0*M_PI/lz : 0.0;
  const int nlow_sqr = nlow*nlow;
  const int nhigh_sqr = nhigh*nhigh;

  mode_count = 0;
  for (int nkx = min_kx; nkx <= max_kx; ++nkx) {
    for (int nky = min_ky; nky <= max_ky; ++nky) {
      for (int nkz = min_kz; nkz <= max_kz; ++nkz) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        if (!IsCanonicalMode(nkx, nky, nkz)) continue;
        const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) ++mode_count;
      }
    }
  }
  if (mode_count == 0) {
    FatalInitialPerturbationError("mode_count is zero; check nlow/nhigh and mode bounds");
  }

  Kokkos::realloc(kx_mode, mode_count);
  Kokkos::realloc(ky_mode, mode_count);
  Kokkos::realloc(kz_mode, mode_count);

  int nmode = 0;
  for (int nkx = min_kx; nkx <= max_kx; ++nkx) {
    for (int nky = min_ky; nky <= max_ky; ++nky) {
      for (int nkz = min_kz; nkz <= max_kz; ++nkz) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        if (!IsCanonicalMode(nkx, nky, nkz)) continue;
        const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) {
          kx_mode.h_view(nmode) = dkx*static_cast<Real>(nkx);
          ky_mode.h_view(nmode) = dky*static_cast<Real>(nky);
          kz_mode.h_view(nmode) = dkz*static_cast<Real>(nkz);
          ++nmode;
        }
      }
    }
  }
  kx_mode.template modify<HostMemSpace>();
  ky_mode.template modify<HostMemSpace>();
  kz_mode.template modify<HostMemSpace>();
  kx_mode.template sync<DevExeSpace>();
  ky_mode.template sync<DevExeSpace>();
  kz_mode.template sync<DevExeSpace>();
}

void InitialPerturbations::GenerateAmplitudes(DualArray2D<Real> &amp_real,
                                              DualArray2D<Real> &amp_imag,
                                              const int ncomp,
                                              const bool project_velocity,
                                              const bool vector_potential) {
  Kokkos::realloc(amp_real, ncomp, mode_count);
  Kokkos::realloc(amp_imag, ncomp, mode_count);

  for (int n = 0; n < mode_count; ++n) {
    const Real kx = kx_mode.h_view(n);
    const Real ky = ky_mode.h_view(n);
    const Real kz = kz_mode.h_view(n);
    const Real k[3] = {kx, ky, kz};
    const Real kiso = std::sqrt(kx*kx + ky*ky + kz*kz);
    Real norm = 0.0;
    if (kiso > 1.0e-16) {
      norm = 1.0/std::pow(kiso, 0.5*(spectral_slope + 2.0));
      if (vector_potential) {
        norm /= kiso;
      }
    }

    Real dot_real = 0.0;
    Real dot_imag = 0.0;
    for (int c = 0; c < ncomp; ++c) {
      const Real ar = norm*RanGaussianSt(&rstate);
      const Real ai = norm*RanGaussianSt(&rstate);
      amp_real.h_view(c,n) = ar;
      amp_imag.h_view(c,n) = ai;
      if (project_velocity && c < 3) {
        dot_real += k[c]*ar;
        dot_imag += k[c]*ai;
      }
    }

    if (project_velocity && kiso > 1.0e-16) {
      for (int c = 0; c < ncomp; ++c) {
        const Real div_real = k[c]*dot_real/(kiso*kiso);
        const Real div_imag = k[c]*dot_imag/(kiso*kiso);
        const Real sol_real = amp_real.h_view(c,n) - div_real;
        const Real sol_imag = amp_imag.h_view(c,n) - div_imag;
        amp_real.h_view(c,n) = f_solenoidal*sol_real + (1.0 - f_solenoidal)*div_real;
        amp_imag.h_view(c,n) = f_solenoidal*sol_imag + (1.0 - f_solenoidal)*div_imag;
      }
    }
  }

  amp_real.template modify<HostMemSpace>();
  amp_imag.template modify<HostMemSpace>();
  amp_real.template sync<DevExeSpace>();
  amp_imag.template sync<DevExeSpace>();
}

void InitialPerturbations::BuildDensityField() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;

  auto field = rho_pert;
  auto kx = kx_mode.d_view;
  auto ky = ky_mode.d_view;
  auto kz = kz_mode.d_view;
  auto ar = rho_amp_real.d_view;
  auto ai = rho_amp_imag.d_view;
  auto mb_size = pmy_pack->pmb->mb_size;

  const int nmode = mode_count;
  const int loc_mode = localization_mode;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const Real c1 = x1_center, c2 = x2_center, c3 = x3_center;
  const Real s1 = x1_scale, s2 = x2_scale, s3 = x3_scale;
  const Real o1 = domain_x1min, o2 = domain_x2min, o3 = domain_x3min;

  par_for("initial_density_perturbation", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, nx1, mb_size.d_view(m).x1min,
                                mb_size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, nx2, mb_size.d_view(m).x2min,
                                mb_size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, nx3, mb_size.d_view(m).x3min,
                                mb_size.d_view(m).x3max);
    const Real win = PerturbationWindow(x1, x2, x3, loc_mode, c1, c2, c3,
                                        s1, s2, s3, multi_d, three_d);
    field(m,k,j,i) = win*FourierComponent(kx, ky, kz, ar, ai, 0, nmode,
                                          x1, x2, x3, o1, o2, o3);
  });
}

void InitialPerturbations::BuildVelocityField() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;

  auto field = vel_pert;
  auto kx = kx_mode.d_view;
  auto ky = ky_mode.d_view;
  auto kz = kz_mode.d_view;
  auto ar = vel_amp_real.d_view;
  auto ai = vel_amp_imag.d_view;
  auto mb_size = pmy_pack->pmb->mb_size;

  const int nmode = mode_count;
  const int loc_mode = localization_mode;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const Real c1 = x1_center, c2 = x2_center, c3 = x3_center;
  const Real s1 = x1_scale, s2 = x2_scale, s3 = x3_scale;
  const Real o1 = domain_x1min, o2 = domain_x2min, o3 = domain_x3min;

  par_for("initial_velocity_perturbation", DevExeSpace(),
          0, nmb1, 0, 2, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int c, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, nx1, mb_size.d_view(m).x1min,
                                mb_size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, nx2, mb_size.d_view(m).x2min,
                                mb_size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, nx3, mb_size.d_view(m).x3min,
                                mb_size.d_view(m).x3max);
    const Real win = PerturbationWindow(x1, x2, x3, loc_mode, c1, c2, c3,
                                        s1, s2, s3, multi_d, three_d);
    field(m,c,k,j,i) = win*FourierComponent(kx, ky, kz, ar, ai, c, nmode,
                                            x1, x2, x3, o1, o2, o3);
  });
}

void InitialPerturbations::BuildVectorPotential() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;

  auto a1 = avec.x1e;
  auto a2 = avec.x2e;
  auto a3 = avec.x3e;
  auto kx = kx_mode.d_view;
  auto ky = ky_mode.d_view;
  auto kz = kz_mode.d_view;
  auto ar = avec_amp_real.d_view;
  auto ai = avec_amp_imag.d_view;
  auto mb_size = pmy_pack->pmb->mb_size;

  const int nmode = mode_count;
  const int loc_mode = localization_mode;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const Real c1 = x1_center, c2 = x2_center, c3 = x3_center;
  const Real s1 = x1_scale, s2 = x2_scale, s3 = x3_scale;
  const Real o1 = domain_x1min, o2 = domain_x2min, o3 = domain_x3min;

  par_for("initial_vector_potential_a1", DevExeSpace(),
          0, nmb1, 0, ncells3, 0, ncells2, 0, ncells1 - 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, nx1, mb_size.d_view(m).x1min,
                                mb_size.d_view(m).x1max);
    const Real x2 = LeftEdgeX(j - js, nx2, mb_size.d_view(m).x2min,
                              mb_size.d_view(m).x2max);
    const Real x3 = LeftEdgeX(k - ks, nx3, mb_size.d_view(m).x3min,
                              mb_size.d_view(m).x3max);
    const Real win = PerturbationWindow(x1, x2, x3, loc_mode, c1, c2, c3,
                                        s1, s2, s3, multi_d, three_d);
    a1(m,k,j,i) = win*FourierComponent(kx, ky, kz, ar, ai, 0, nmode,
                                       x1, x2, x3, o1, o2, o3);
  });

  par_for("initial_vector_potential_a2", DevExeSpace(),
          0, nmb1, 0, ncells3, 0, ncells2 - 1, 0, ncells1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = LeftEdgeX(i - is, nx1, mb_size.d_view(m).x1min,
                              mb_size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, nx2, mb_size.d_view(m).x2min,
                                mb_size.d_view(m).x2max);
    const Real x3 = LeftEdgeX(k - ks, nx3, mb_size.d_view(m).x3min,
                              mb_size.d_view(m).x3max);
    const Real win = PerturbationWindow(x1, x2, x3, loc_mode, c1, c2, c3,
                                        s1, s2, s3, multi_d, three_d);
    a2(m,k,j,i) = win*FourierComponent(kx, ky, kz, ar, ai, 1, nmode,
                                       x1, x2, x3, o1, o2, o3);
  });

  par_for("initial_vector_potential_a3", DevExeSpace(),
          0, nmb1, 0, ncells3 - 1, 0, ncells2, 0, ncells1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = LeftEdgeX(i - is, nx1, mb_size.d_view(m).x1min,
                              mb_size.d_view(m).x1max);
    const Real x2 = LeftEdgeX(j - js, nx2, mb_size.d_view(m).x2min,
                              mb_size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, nx3, mb_size.d_view(m).x3min,
                                mb_size.d_view(m).x3max);
    const Real win = PerturbationWindow(x1, x2, x3, loc_mode, c1, c2, c3,
                                        s1, s2, s3, multi_d, three_d);
    a3(m,k,j,i) = win*FourierComponent(kx, ky, kz, ar, ai, 2, nmode,
                                       x1, x2, x3, o1, o2, o3);
  });
}

void InitialPerturbations::CurlVectorPotential() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;

  auto a1 = avec.x1e;
  auto a2 = avec.x2e;
  auto a3 = avec.x3e;
  auto b1 = b_pert.x1f;
  auto b2 = b_pert.x2f;
  auto b3 = b_pert.x3f;
  auto mb_size = pmy_pack->pmb->mb_size;

  par_for("initial_bpert_x1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real curl = 0.0;
    if (multi_d) {
      curl += (a3(m,k,j+1,i) - a3(m,k,j,i))/mb_size.d_view(m).dx2;
    }
    if (three_d) {
      curl -= (a2(m,k+1,j,i) - a2(m,k,j,i))/mb_size.d_view(m).dx3;
    }
    b1(m,k,j,i) = curl;
  });

  par_for("initial_bpert_x2", DevExeSpace(), 0, nmb1, ks, ke, js, je + 1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real curl = -(a3(m,k,j,i+1) - a3(m,k,j,i))/mb_size.d_view(m).dx1;
    if (three_d) {
      curl += (a1(m,k+1,j,i) - a1(m,k,j,i))/mb_size.d_view(m).dx3;
    }
    b2(m,k,j,i) = curl;
  });

  par_for("initial_bpert_x3", DevExeSpace(), 0, nmb1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real curl = (a2(m,k,j,i+1) - a2(m,k,j,i))/mb_size.d_view(m).dx1;
    if (multi_d) {
      curl -= (a1(m,k,j+1,i) - a1(m,k,j,i))/mb_size.d_view(m).dx2;
    }
    b3(m,k,j,i) = curl;
  });
}

Real InitialPerturbations::CellScalarRMS(const DvceArray4D<Real> &field,
                                         const bool remove_mean) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  auto mb_size = pmy_pack->pmb->mb_size;

  Real sum = 0.0;
  Real vol = 0.0;
  if (remove_mean) {
    Kokkos::parallel_reduce("initial_scalar_mean",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int idx, Real &sum_, Real &vol_) {
      const int m = idx/nkji;
      const int k0 = (idx - m*nkji)/nji;
      const int j0 = (idx - m*nkji - k0*nji)/nx1;
      const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
      const int j = j0 + js;
      const int k = k0 + ks;
      const Real dv = mb_size.d_view(m).dx1*mb_size.d_view(m).dx2*
                      mb_size.d_view(m).dx3;
      sum_ += field(m,k,j,i)*dv;
      vol_ += dv;
    }, Kokkos::Sum<Real>(sum), Kokkos::Sum<Real>(vol));
    Real vals[2] = {sum, vol};
    ReduceRealSum(vals, 2);
    const Real mean = (vals[1] > 0.0) ? vals[0]/vals[1] : 0.0;
    auto field_mut = field;
    par_for("initial_scalar_subtract_mean", DevExeSpace(),
            0, nmb - 1, ks, ks + nx3 - 1, js, js + nx2 - 1, is, is + nx1 - 1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      field_mut(m,k,j,i) -= mean;
    });
  }

  sum = 0.0;
  vol = 0.0;
  Kokkos::parallel_reduce("initial_scalar_rms",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &sum_, Real &vol_) {
    const int m = idx/nkji;
    const int k0 = (idx - m*nkji)/nji;
    const int j0 = (idx - m*nkji - k0*nji)/nx1;
    const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real dv = mb_size.d_view(m).dx1*mb_size.d_view(m).dx2*
                    mb_size.d_view(m).dx3;
    sum_ += SQR(field(m,k,j,i))*dv;
    vol_ += dv;
  }, Kokkos::Sum<Real>(sum), Kokkos::Sum<Real>(vol));
  Real vals[2] = {sum, vol};
  ReduceRealSum(vals, 2);
  return (vals[1] > 0.0 && vals[0] > 0.0) ? std::sqrt(vals[0]/vals[1]) : 0.0;
}

Real InitialPerturbations::CellVectorRMS(const DvceArray5D<Real> &field,
                                         const bool remove_mean) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  auto mb_size = pmy_pack->pmb->mb_size;

  if (remove_mean) {
    Real s0 = 0.0, s1 = 0.0, s2 = 0.0, vol = 0.0;
    Kokkos::parallel_reduce("initial_vector_mean",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int idx, Real &sx, Real &sy, Real &sz, Real &vol_) {
      const int m = idx/nkji;
      const int k0 = (idx - m*nkji)/nji;
      const int j0 = (idx - m*nkji - k0*nji)/nx1;
      const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
      const int j = j0 + js;
      const int k = k0 + ks;
      const Real dv = mb_size.d_view(m).dx1*mb_size.d_view(m).dx2*
                      mb_size.d_view(m).dx3;
      sx += field(m,0,k,j,i)*dv;
      sy += field(m,1,k,j,i)*dv;
      sz += field(m,2,k,j,i)*dv;
      vol_ += dv;
    }, Kokkos::Sum<Real>(s0), Kokkos::Sum<Real>(s1),
       Kokkos::Sum<Real>(s2), Kokkos::Sum<Real>(vol));
    Real vals[4] = {s0, s1, s2, vol};
    ReduceRealSum(vals, 4);
    const Real mx = (vals[3] > 0.0) ? vals[0]/vals[3] : 0.0;
    const Real my = (vals[3] > 0.0) ? vals[1]/vals[3] : 0.0;
    const Real mz = (vals[3] > 0.0) ? vals[2]/vals[3] : 0.0;
    auto field_mut = field;
    par_for("initial_vector_subtract_mean", DevExeSpace(),
            0, nmb - 1, ks, ks + nx3 - 1, js, js + nx2 - 1, is, is + nx1 - 1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      field_mut(m,0,k,j,i) -= mx;
      field_mut(m,1,k,j,i) -= my;
      field_mut(m,2,k,j,i) -= mz;
    });
  }

  Real sum = 0.0;
  Real vol = 0.0;
  Kokkos::parallel_reduce("initial_vector_rms",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &sum_, Real &vol_) {
    const int m = idx/nkji;
    const int k0 = (idx - m*nkji)/nji;
    const int j0 = (idx - m*nkji - k0*nji)/nx1;
    const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real dv = mb_size.d_view(m).dx1*mb_size.d_view(m).dx2*
                    mb_size.d_view(m).dx3;
    sum_ += (SQR(field(m,0,k,j,i)) + SQR(field(m,1,k,j,i)) +
             SQR(field(m,2,k,j,i)))*dv;
    vol_ += dv;
  }, Kokkos::Sum<Real>(sum), Kokkos::Sum<Real>(vol));
  Real vals[2] = {sum, vol};
  ReduceRealSum(vals, 2);
  return (vals[1] > 0.0 && vals[0] > 0.0) ? std::sqrt(vals[0]/vals[1]) : 0.0;
}

Real InitialPerturbations::FaceCenteredBRMS(const DvceFaceFld4D<Real> &field) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  auto b1 = field.x1f;
  auto b2 = field.x2f;
  auto b3 = field.x3f;
  auto mb_size = pmy_pack->pmb->mb_size;

  Real sum = 0.0;
  Real vol = 0.0;
  Kokkos::parallel_reduce("initial_b_rms", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &sum_, Real &vol_) {
    const int m = idx/nkji;
    const int k0 = (idx - m*nkji)/nji;
    const int j0 = (idx - m*nkji - k0*nji)/nx1;
    const int i = (idx - m*nkji - k0*nji - j0*nx1) + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real dv = mb_size.d_view(m).dx1*mb_size.d_view(m).dx2*mb_size.d_view(m).dx3;
    const Real dbx = 0.5*(b1(m,k,j,i) + b1(m,k,j,i+1));
    const Real dby = 0.5*(b2(m,k,j,i) + b2(m,k,j+1,i));
    const Real dbz = 0.5*(b3(m,k,j,i) + b3(m,k+1,j,i));
    sum_ += (SQR(dbx) + SQR(dby) + SQR(dbz))*dv;
    vol_ += dv;
  }, Kokkos::Sum<Real>(sum), Kokkos::Sum<Real>(vol));
  Real vals[2] = {sum, vol};
  ReduceRealSum(vals, 2);
  return (vals[1] > 0.0 && vals[0] > 0.0) ? std::sqrt(vals[0]/vals[1]) : 0.0;
}

void InitialPerturbations::ApplyDensityField() {
  BuildDensityField();
  const Real measured_rms = CellScalarRMS(rho_pert, remove_density_mean);
  if (measured_rms <= 0.0 || density_rms <= 0.0) return;
  const Real scale = density_rms/measured_rms;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool fractional = density_fractional;

  DvceArray5D<Real> u0;
  int nfluid = 0;
  int nscalars = 0;
  EOS_Data eos;
  bool is_mhd = false;
  DvceFaceFld4D<Real> b0("unused_b0", 1, 1, 1, 1);
  if (pmy_pack->pmhd != nullptr) {
    u0 = pmy_pack->pmhd->u0;
    nfluid = pmy_pack->pmhd->nmhd;
    nscalars = pmy_pack->pmhd->nscalars;
    eos = pmy_pack->pmhd->peos->eos_data;
    is_mhd = true;
    b0 = pmy_pack->pmhd->b0;
  } else {
    u0 = pmy_pack->phydro->u0;
    nfluid = pmy_pack->phydro->nhydro;
    nscalars = pmy_pack->phydro->nscalars;
    eos = pmy_pack->phydro->peos->eos_data;
  }

  auto field = rho_pert;
  const Real dfloor = eos.dfloor;
  const bool ideal = eos.is_ideal;
  const int nvar = nfluid + nscalars;

  par_for("apply_initial_density_perturbation", DevExeSpace(),
          0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real old_rho = u0(m,IDN,k,j,i);
    Real new_rho = fractional ? old_rho*(1.0 + scale*field(m,k,j,i))
                              : old_rho + scale*field(m,k,j,i);
    new_rho = fmax(new_rho, dfloor);
    const Real ratio = (old_rho > 0.0) ? new_rho/old_rho : 1.0;

    const Real old_m1 = u0(m,IM1,k,j,i);
    const Real old_m2 = u0(m,IM2,k,j,i);
    const Real old_m3 = u0(m,IM3,k,j,i);
    const Real old_ke = 0.5*(SQR(old_m1) + SQR(old_m2) + SQR(old_m3))/old_rho;
    Real old_me = 0.0;
    if (is_mhd) {
      const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      old_me = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    }

    u0(m,IDN,k,j,i) = new_rho;
    u0(m,IM1,k,j,i) = ratio*old_m1;
    u0(m,IM2,k,j,i) = ratio*old_m2;
    u0(m,IM3,k,j,i) = ratio*old_m3;
    for (int n = nfluid; n < nvar; ++n) {
      u0(m,n,k,j,i) *= ratio;
    }
    if (ideal) {
      Real gas_e = u0(m,IEN,k,j,i) - old_ke - old_me;
      gas_e = fmax(gas_e, 0.0);
      const Real new_ke = ratio*old_ke;
      u0(m,IEN,k,j,i) = ratio*gas_e + new_ke + old_me;
    }
  });
}

void InitialPerturbations::ApplyVelocityField() {
  BuildVelocityField();
  const Real measured_rms = CellVectorRMS(vel_pert, remove_velocity_mean);
  if (measured_rms <= 0.0 || velocity_rms <= 0.0) return;
  const Real scale = velocity_rms/measured_rms;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;

  DvceArray5D<Real> u0;
  EOS_Data eos;
  bool is_mhd = false;
  DvceFaceFld4D<Real> b0("unused_b0", 1, 1, 1, 1);
  if (pmy_pack->pmhd != nullptr) {
    u0 = pmy_pack->pmhd->u0;
    eos = pmy_pack->pmhd->peos->eos_data;
    is_mhd = true;
    b0 = pmy_pack->pmhd->b0;
  } else {
    u0 = pmy_pack->phydro->u0;
    eos = pmy_pack->phydro->peos->eos_data;
  }

  auto field = vel_pert;
  const bool ideal = eos.is_ideal;

  par_for("apply_initial_velocity_perturbation", DevExeSpace(),
          0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real rho = u0(m,IDN,k,j,i);
    const Real old_m1 = u0(m,IM1,k,j,i);
    const Real old_m2 = u0(m,IM2,k,j,i);
    const Real old_m3 = u0(m,IM3,k,j,i);
    const Real old_ke = 0.5*(SQR(old_m1) + SQR(old_m2) + SQR(old_m3))/rho;
    Real old_me = 0.0;
    if (is_mhd) {
      const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      old_me = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    }

    const Real v1 = old_m1/rho + scale*field(m,0,k,j,i);
    const Real v2 = old_m2/rho + scale*field(m,1,k,j,i);
    const Real v3 = old_m3/rho + scale*field(m,2,k,j,i);
    u0(m,IM1,k,j,i) = rho*v1;
    u0(m,IM2,k,j,i) = rho*v2;
    u0(m,IM3,k,j,i) = rho*v3;
    if (ideal) {
      Real gas_e = u0(m,IEN,k,j,i) - old_ke - old_me;
      gas_e = fmax(gas_e, 0.0);
      const Real new_ke = 0.5*rho*(SQR(v1) + SQR(v2) + SQR(v3));
      u0(m,IEN,k,j,i) = gas_e + new_ke + old_me;
    }
  });
}

void InitialPerturbations::ApplyMagneticField() {
  BuildVectorPotential();
  CurlVectorPotential();
  const Real measured_rms = FaceCenteredBRMS(b_pert);
  if (measured_rms <= 0.0 || magnetic_rms <= 0.0) return;
  const Real scale = magnetic_rms/measured_rms;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;

  auto b0 = pmy_pack->pmhd->b0;
  auto bcc0 = pmy_pack->pmhd->bcc0;
  auto u0 = pmy_pack->pmhd->u0;
  auto db1 = b_pert.x1f;
  auto db2 = b_pert.x2f;
  auto db3 = b_pert.x3f;
  const bool ideal = pmy_pack->pmhd->peos->eos_data.is_ideal;

  if (ideal) {
    par_for("initial_magnetic_energy_update", DevExeSpace(),
            0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
      const Real by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      const Real dbx = 0.5*scale*(db1(m,k,j,i) + db1(m,k,j,i+1));
      const Real dby = 0.5*scale*(db2(m,k,j,i) + db2(m,k,j+1,i));
      const Real dbz = 0.5*scale*(db3(m,k,j,i) + db3(m,k+1,j,i));
      const Real old_me = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
      const Real new_me = 0.5*(SQR(bx + dbx) + SQR(by + dby) + SQR(bz + dbz));
      u0(m,IEN,k,j,i) += new_me - old_me;
    });
  }

  par_for("apply_initial_bx_perturbation", DevExeSpace(),
          0, nmb1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) += scale*db1(m,k,j,i);
  });
  par_for("apply_initial_by_perturbation", DevExeSpace(),
          0, nmb1, ks, ke, js, je + 1, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x2f(m,k,j,i) += scale*db2(m,k,j,i);
  });
  par_for("apply_initial_bz_perturbation", DevExeSpace(),
          0, nmb1, ks, ke + 1, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x3f(m,k,j,i) += scale*db3(m,k,j,i);
  });

  par_for("initial_bcc_after_perturbation", DevExeSpace(),
          0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
  });
}

void InitialPerturbations::Apply() {
  if (!AnyEnabled()) return;

  if (DensityEnabled()) {
    GenerateAmplitudes(rho_amp_real, rho_amp_imag, 1, false, false);
    ApplyDensityField();
  }
  if (VelocityEnabled()) {
    GenerateAmplitudes(vel_amp_real, vel_amp_imag, 3, true, false);
    ApplyVelocityField();
  }
  if (MagneticEnabled()) {
    GenerateAmplitudes(avec_amp_real, avec_amp_imag, 3, false, true);
    ApplyMagneticField();
  }
  Kokkos::fence();
}

void ApplyInitialPerturbations(Mesh *pm, ParameterInput *pin) {
  if (pm == nullptr || pin == nullptr || pm->pmb_pack == nullptr) {
    return;
  }

  const bool has_plural = pin->DoesBlockExist("initial_perturbations");
  const bool has_singular = pin->DoesBlockExist("initial_perturbation");
  if (has_plural && has_singular) {
    FatalInitialPerturbationError("use only one of <initial_perturbations> or "
                                  "<initial_perturbation>");
  }
  if (has_plural) {
    InitialPerturbations perturbations(pm->pmb_pack, pin, "initial_perturbations");
    perturbations.Apply();
  } else if (has_singular) {
    InitialPerturbations perturbation(pm->pmb_pack, pin, "initial_perturbation");
    perturbation.Apply();
  }
}
