//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_mixing.cpp
//  \brief Problem generator for scalar mixing in frozen turbulent velocity field
//
//  Passive scalars with selectable source terms (mean gradient, reaction, sponge).
//  This simulates mixing of background gradients and/or reactive reservoirs by turbulence.
//
//  When scalar_only=true in <hydro>, velocity is frozen and only the scalars evolve.

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen.hpp"
#include "globals.hpp"
#include "utils/random.hpp"

#include <Kokkos_Random.hpp>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

// Function declarations
void ScalarForcingSource(Mesh* pm, const Real bdt);
void ScalarMixingHistory(HistoryData *pdata, Mesh *pm);

namespace {

constexpr Real kTiny = 1.0e-16;

enum class VelocityMethod { Projection, StreamFunction };

struct TurbulenceConfig {
  VelocityMethod method = VelocityMethod::Projection;
  Real v_rms = 1.0;
  int nlow = 1;
  int nhigh = 4;
  Real expo = 5.0/3.0;
  Real sol_frac = 1.0;
  int rseed = 12345;
  Real k_crit = 16.0;

  int nx1 = 1;
  int nx2 = 1;
  int nx3 = 1;
  int global_nx1 = 1;
  int global_nx2 = 1;
  int global_nx3 = 1;

  bool mesh_multi_d = false;
  bool mesh_three_d = false;
  bool projection_true_2d = false;
  bool stream_2d = false;
  bool active_v3 = false;
  bool divfree_scalar_flux = false;
  bool dump_generator_diagnostics = false;
  bool stream_allow_poor_fit = false;

  Real lx = 1.0;
  Real ly = 1.0;
  Real lz = 1.0;
  Real x1min = 0.0;
  Real x1max = 1.0;
  Real x2min = 0.0;
  Real x2max = 1.0;
  Real x3min = 0.0;
  Real x3max = 1.0;
  Real dkx = 0.0;
  Real dky = 0.0;
  Real dkz = 0.0;
  Real k_crit_mag = 0.0;
  std::string basename;
};

struct ModeCatalog {
  int total_modes = 0;
  std::vector<int> nkx, nky, nkz, shell;
  std::vector<Real> prob, boost;
  std::vector<Real> kx, ky, kz, kiso;

  int KeptModes() const { return static_cast<int>(nkx.size()); }
};

struct VectorCoefficients {
  std::vector<Real> aka0, aka1, aka2;
  std::vector<Real> akb0, akb1, akb2;
};

struct ScalarCoefficients {
  std::vector<Real> aka, akb;
};

struct SignedMode {
  int nkx, nky, nkz;
  int shell;
  Real kx, ky, kz;
  Real weight;
};

struct SignedScalarCoeffMode {
  int nkx, nky, nkz;
  int shell;
  Real kx, ky, kz;
  std::complex<Real> coeff;
};

struct AccumulatedVelocityMode {
  int shell = 0;
  std::array<std::complex<Real>, 3> coeff = {
      std::complex<Real>(0.0, 0.0), std::complex<Real>(0.0, 0.0),
      std::complex<Real>(0.0, 0.0)};
};

struct ShellFitEval {
  std::vector<Real> pred;
  std::vector<Real> grad1;
  std::vector<Real> grad2;
  Real loss = 0.0;
};

enum class ClebschSeedModel { LocalTriad, CamilloSymmetric, Flat };

enum class ClebschSector {
  LocalLocal = 0,
  LowHigh = 1,
  HighLow = 2,
  HighHighCancel = 3
};

constexpr int kClebschSectorCount = 4;

struct ShellFitStartSummary {
  ClebschSeedModel seed_model = ClebschSeedModel::LocalTriad;
  Real initial_loss = 0.0;
  Real final_loss = 0.0;
  Real fitted_slope = 0.0;
  Real leakage_fraction = 0.0;
  bool converged = false;
  bool chosen = false;
};

struct ShellFitResult {
  std::vector<Real> shell_amp1;
  std::vector<Real> shell_amp2;
  std::vector<Real> shell_power1;
  std::vector<Real> shell_power2;
  std::vector<Real> predicted_shell_energy;
  std::vector<ShellFitStartSummary> fit_starts;
  Real initial_loss = 0.0;
  Real final_loss = 0.0;
  Real leakage_fraction = 0.0;
  Real fitted_slope = 0.0;
  bool converged = false;
  ClebschSeedModel chosen_seed_model = ClebschSeedModel::LocalTriad;
};

struct RealizedSpectrumSummary {
  std::vector<Real> shell_energy;
  Real loss = 0.0;
  Real leakage_fraction = 0.0;
  Real fitted_slope = 0.0;
};

struct ClebschRealization {
  ScalarCoefficients phi1_coeffs;
  ScalarCoefficients phi2_coeffs;
  RealizedSpectrumSummary summary;
  int trial_index = 0;
  int trial_count = 1;
};

struct ClebschVelocityModes {
  VectorCoefficients coeffs;
  std::vector<Real> shell_energy;
};

struct ClebschSectorResponse {
  std::array<std::vector<Real>, kClebschSectorCount> raw_shell_energy;
  std::array<std::vector<Real>, kClebschSectorCount> fitted_shell_energy;
  std::array<Real, kClebschSectorCount> raw_total = {0.0, 0.0, 0.0, 0.0};
  std::array<Real, kClebschSectorCount> fitted_total = {0.0, 0.0, 0.0, 0.0};
};

struct ClebschQualityGate {
  Real deterministic_slope_tolerance = 0.5;
  Real realized_slope_tolerance = 0.5;
  Real leakage_tolerance = 0.05;
  bool deterministic_converged = false;
  bool deterministic_slope_ok = false;
  bool deterministic_leakage_ok = false;
  bool deterministic_pass = false;
  bool realized_slope_ok = false;
  bool realized_leakage_ok = false;
  bool realized_pass = false;
  bool pass = false;
  bool allow_poor_fit = false;
  bool forced_continue = false;
};

struct VelocityStats {
  Real vmean1 = 0.0;
  Real vmean2 = 0.0;
  Real vmean3 = 0.0;
  Real scale = 1.0;
};

enum ScalarForcingMode {
  kScalarNone = 0,
  kScalarMeanGradient = 1,
  kScalarMassConservingAC = 2,
  kScalarSponge = 3,
  kScalarSpongeReaction = 4
};

enum ReactionType {
  kReactionKpp = 1,
  kReactionBistable = 2
};

struct ScalarForcingConfig {
  int mode = kScalarNone;
  Real G = 1.0;
  Real tau = 1.0;
  Real theta1 = 1.0e-6;
  Real theta2 = 0.5;
  Real theta3 = 1.0;
  Real a = 0.5;
  Real chiL0 = 1.0;
  Real chiR0 = 1.0;
  Real thetaL = 1.0e-6;
  Real thetaR = 1.0;
  Real xL0 = 0.0;
  Real xL1 = 0.0;
  Real xR0 = 0.0;
  Real xR1 = 0.0;
  Real mask_w = 0.0;
  int reaction_type = kReactionBistable;
  Real theta_floor = 1.0e-6;
  bool clip_reaction_interval = true;
};

std::vector<ScalarForcingConfig> scalar_forcing_cfg;

bool UseTrue2DVelocity(ParameterInput *pin, Mesh *pm) {
  if (!pm->two_d) return false;
  if (!(pin->DoesParameterExist("mesh", "ix3_bc") &&
        pin->DoesParameterExist("mesh", "ox3_bc"))) {
    return false;
  }
  return (pin->GetString("mesh", "ix3_bc") == "reflect" &&
          pin->GetString("mesh", "ox3_bc") == "reflect");
}

int GetTurbulenceSeed(ParameterInput *pin) {
  if (pin->DoesParameterExist("problem", "turb_rseed")) {
    return pin->GetInteger("problem", "turb_rseed");
  }
  return pin->GetOrAddInteger("problem", "turb_seed", 12345);
}

int MixSeed(int seed, int salt) {
  long long mixed = static_cast<long long>(seed) + 104729LL * salt;
  mixed %= 2147483629LL;
  if (mixed <= 0) mixed += 2147483629LL;
  return static_cast<int>(mixed);
}

[[noreturn]] void FatalProblemSetup(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

int ShellFromIndexSqr(int nsqr) {
  if (nsqr <= 0) return 0;
  return static_cast<int>(std::ceil(std::sqrt(static_cast<Real>(nsqr)) - 1.0e-12));
}

Real TargetShellEnergy(int shell, Real expo) {
  return 1.0 / std::pow(static_cast<Real>(shell), expo);
}

using WallClock = std::chrono::steady_clock;
using WallTimePoint = WallClock::time_point;

Real ElapsedWallSeconds(const WallTimePoint &start_time) {
  return std::chrono::duration<Real>(WallClock::now() - start_time).count();
}

Real GlobalMaxReal(Real local_value) {
#if MPI_PARALLEL_ENABLED
  Real global_value = 0.0;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  return global_value;
#else
  return local_value;
#endif
}

const char *VelocityMethodName(VelocityMethod method) {
  return (method == VelocityMethod::Projection) ? "projection" : "stream_function";
}

const char *ClebschSeedModelName(ClebschSeedModel model) {
  switch (model) {
    case ClebschSeedModel::LocalTriad:
      return "local_triad";
    case ClebschSeedModel::CamilloSymmetric:
      return "camillo_symmetric";
    case ClebschSeedModel::Flat:
      return "flat";
  }
  return "unknown";
}

const char *ClebschSectorName(ClebschSector sector) {
  switch (sector) {
    case ClebschSector::LocalLocal:
      return "local_local";
    case ClebschSector::LowHigh:
      return "low_high";
    case ClebschSector::HighLow:
      return "high_low";
    case ClebschSector::HighHighCancel:
      return "high_high_cancel";
  }
  return "unknown";
}

Real UniformCellCenter(int index, int ncell, Real xmin, Real xmax) {
  if (ncell <= 0) return 0.5 * (xmin + xmax);
  return xmin + (static_cast<Real>(index) + 0.5) * (xmax - xmin) / static_cast<Real>(ncell);
}

std::vector<Real> ComputeScalarShellEnergy(const ModeCatalog &catalog,
                                           const ScalarCoefficients &coeffs,
                                           int nshell_max) {
  std::vector<Real> shell_energy(nshell_max + 1, 0.0);
  for (int n = 0; n < catalog.KeptModes(); ++n) {
    const int shell = catalog.shell[n];
    if (shell < 0 || shell > nshell_max) continue;
    shell_energy[shell] += 0.5 * (coeffs.aka[n] * coeffs.aka[n] +
                                  coeffs.akb[n] * coeffs.akb[n]);
  }
  return shell_energy;
}

std::vector<Real> ComputeVectorShellEnergy(const ModeCatalog &catalog,
                                           const VectorCoefficients &coeffs,
                                           int nshell_max) {
  std::vector<Real> shell_energy(nshell_max + 1, 0.0);
  for (int n = 0; n < catalog.KeptModes(); ++n) {
    const int shell = catalog.shell[n];
    if (shell < 0 || shell > nshell_max) continue;
    shell_energy[shell] += 0.5 * (
        coeffs.aka0[n] * coeffs.aka0[n] + coeffs.akb0[n] * coeffs.akb0[n] +
        coeffs.aka1[n] * coeffs.aka1[n] + coeffs.akb1[n] * coeffs.akb1[n] +
        coeffs.aka2[n] * coeffs.aka2[n] + coeffs.akb2[n] * coeffs.akb2[n]);
  }
  return shell_energy;
}

std::vector<Real> SampleScalarFieldXY(const TurbulenceConfig &cfg,
                                      const ModeCatalog &catalog,
                                      const ScalarCoefficients &coeffs,
                                      int nx, int ny, Real x3_slice) {
  std::vector<Real> field(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0.0);
  for (int j = 0; j < ny; ++j) {
    const Real x2 = UniformCellCenter(j, ny, cfg.x2min, cfg.x2max);
    for (int i = 0; i < nx; ++i) {
      const Real x1 = UniformCellCenter(i, nx, cfg.x1min, cfg.x1max);
      Real value = 0.0;
      for (int n = 0; n < catalog.KeptModes(); ++n) {
        const Real phase = catalog.kx[n] * x1 + catalog.ky[n] * x2 + catalog.kz[n] * x3_slice;
        value += coeffs.aka[n] * std::cos(phase) - coeffs.akb[n] * std::sin(phase);
      }
      field[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + i] = value;
    }
  }
  return field;
}

std::string JsonEscape(const std::string &value) {
  std::string escaped;
  escaped.reserve(value.size() + 8);
  for (char ch : value) {
    switch (ch) {
      case '\\':
        escaped += "\\\\";
        break;
      case '"':
        escaped += "\\\"";
        break;
      case '\n':
        escaped += "\\n";
        break;
      case '\r':
        escaped += "\\r";
        break;
      case '\t':
        escaped += "\\t";
        break;
      default:
        escaped += ch;
        break;
    }
  }
  return escaped;
}

void WriteIndent(std::ostream &os, int indent) {
  os << std::string(indent, ' ');
}

template <typename T>
void WriteJsonNumericArray(std::ostream &os, const std::vector<T> &values) {
  os << "[";
  for (std::size_t n = 0; n < values.size(); ++n) {
    if (n > 0) os << ", ";
    os << values[n];
  }
  os << "]";
}

void WriteJsonModeCatalog(std::ostream &os, const ModeCatalog &catalog, int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  WriteIndent(os, indent + 2);
  os << "\"total_modes\": " << catalog.total_modes << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"kept_modes\": " << catalog.KeptModes() << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"nkx\": ";
  WriteJsonNumericArray(os, catalog.nkx);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"nky\": ";
  WriteJsonNumericArray(os, catalog.nky);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"nkz\": ";
  WriteJsonNumericArray(os, catalog.nkz);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"shell\": ";
  WriteJsonNumericArray(os, catalog.shell);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"prob\": ";
  WriteJsonNumericArray(os, catalog.prob);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"boost\": ";
  WriteJsonNumericArray(os, catalog.boost);
  os << "\n";
  WriteIndent(os, indent);
  os << "}";
}

void WriteJsonScalarCoefficients(std::ostream &os, const ScalarCoefficients &coeffs,
                                 int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  WriteIndent(os, indent + 2);
  os << "\"aka\": ";
  WriteJsonNumericArray(os, coeffs.aka);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"akb\": ";
  WriteJsonNumericArray(os, coeffs.akb);
  os << "\n";
  WriteIndent(os, indent);
  os << "}";
}

void WriteJsonVectorCoefficients(std::ostream &os, const VectorCoefficients &coeffs,
                                 int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  WriteIndent(os, indent + 2);
  os << "\"aka0\": ";
  WriteJsonNumericArray(os, coeffs.aka0);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"aka1\": ";
  WriteJsonNumericArray(os, coeffs.aka1);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"aka2\": ";
  WriteJsonNumericArray(os, coeffs.aka2);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"akb0\": ";
  WriteJsonNumericArray(os, coeffs.akb0);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"akb1\": ";
  WriteJsonNumericArray(os, coeffs.akb1);
  os << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"akb2\": ";
  WriteJsonNumericArray(os, coeffs.akb2);
  os << "\n";
  WriteIndent(os, indent);
  os << "}";
}

void WriteJsonField2D(std::ostream &os, const std::vector<Real> &field,
                      int nx, int ny, Real x3_slice, int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  WriteIndent(os, indent + 2);
  os << "\"shape\": [" << ny << ", " << nx << "],\n";
  WriteIndent(os, indent + 2);
  os << "\"x3_slice\": " << x3_slice << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"data\": ";
  WriteJsonNumericArray(os, field);
  os << "\n";
  WriteIndent(os, indent);
  os << "}";
}

void WriteJsonClebschFitStarts(std::ostream &os,
                               const std::vector<ShellFitStartSummary> &starts,
                               int indent) {
  WriteIndent(os, indent);
  os << "[\n";
  for (std::size_t idx = 0; idx < starts.size(); ++idx) {
    const auto &start = starts[idx];
    WriteIndent(os, indent + 2);
    os << "{\n";
    WriteIndent(os, indent + 4);
    os << "\"seed_model\": \"" << ClebschSeedModelName(start.seed_model) << "\",\n";
    WriteIndent(os, indent + 4);
    os << "\"initial_loss\": " << start.initial_loss << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"final_loss\": " << start.final_loss << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"fitted_slope\": " << start.fitted_slope << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"estimated_out_of_band_leakage\": " << start.leakage_fraction << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"converged\": " << (start.converged ? "true" : "false") << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"chosen\": " << (start.chosen ? "true" : "false") << "\n";
    WriteIndent(os, indent + 2);
    os << "}";
    if (idx + 1 < starts.size()) os << ",";
    os << "\n";
  }
  WriteIndent(os, indent);
  os << "]";
}

void WriteJsonClebschSectorBreakdown(std::ostream &os,
                                     const std::array<std::vector<Real>, kClebschSectorCount> &shells,
                                     const std::array<Real, kClebschSectorCount> &totals,
                                     int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  for (int sector_idx = 0; sector_idx < kClebschSectorCount; ++sector_idx) {
    const ClebschSector sector = static_cast<ClebschSector>(sector_idx);
    WriteIndent(os, indent + 2);
    os << "\"" << ClebschSectorName(sector) << "\": {\n";
    WriteIndent(os, indent + 4);
    os << "\"shell_energy\": ";
    WriteJsonNumericArray(os, shells[sector_idx]);
    os << ",\n";
    WriteIndent(os, indent + 4);
    os << "\"total\": " << totals[sector_idx] << "\n";
    WriteIndent(os, indent + 2);
    os << "}";
    if (sector_idx + 1 < kClebschSectorCount) os << ",";
    os << "\n";
  }
  WriteIndent(os, indent);
  os << "}";
}

void WriteDiagnosticsMetadata(std::ostream &os, const TurbulenceConfig &cfg,
                              const std::string &construction,
                              const VelocityStats &stats,
                              Real init_wall_seconds, int indent) {
  WriteIndent(os, indent);
  os << "{\n";
  WriteIndent(os, indent + 2);
  os << "\"basename\": \"" << JsonEscape(cfg.basename) << "\",\n";
  WriteIndent(os, indent + 2);
  os << "\"velocity_method\": \"" << VelocityMethodName(cfg.method) << "\",\n";
  WriteIndent(os, indent + 2);
  os << "\"construction\": \"" << construction << "\",\n";
  WriteIndent(os, indent + 2);
  os << "\"active_velocity_components\": " << (cfg.active_v3 ? 3 : 2) << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"projection_true_2d\": " << (cfg.projection_true_2d ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"stream_2d\": " << (cfg.stream_2d ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"divfree_scalar_flux\": " << (cfg.divfree_scalar_flux ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"init_wall_seconds\": " << init_wall_seconds << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"mesh\": {\n";
  WriteIndent(os, indent + 4);
  os << "\"global_nx1\": " << cfg.global_nx1 << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"global_nx2\": " << cfg.global_nx2 << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"global_nx3\": " << cfg.global_nx3 << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x1min\": " << cfg.x1min << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x1max\": " << cfg.x1max << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x2min\": " << cfg.x2min << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x2max\": " << cfg.x2max << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x3min\": " << cfg.x3min << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"x3max\": " << cfg.x3max << "\n";
  WriteIndent(os, indent + 2);
  os << "},\n";
  WriteIndent(os, indent + 2);
  os << "\"target\": {\n";
  WriteIndent(os, indent + 4);
  os << "\"v_rms\": " << cfg.v_rms << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"nlow\": " << cfg.nlow << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"nhigh\": " << cfg.nhigh << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"expo\": " << cfg.expo << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"sol_frac\": " << cfg.sol_frac << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"rseed\": " << cfg.rseed << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"k_crit\": " << cfg.k_crit << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"k_crit_mag\": " << cfg.k_crit_mag << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"turb_stream_allow_poor_fit\": "
     << (cfg.stream_allow_poor_fit ? "true" : "false") << "\n";
  WriteIndent(os, indent + 2);
  os << "},\n";
  WriteIndent(os, indent + 2);
  os << "\"normalization\": {\n";
  WriteIndent(os, indent + 4);
  os << "\"vmean\": [" << stats.vmean1 << ", " << stats.vmean2 << ", "
     << stats.vmean3 << "],\n";
  WriteIndent(os, indent + 4);
  os << "\"scale\": " << stats.scale << "\n";
  WriteIndent(os, indent + 2);
  os << "}\n";
  WriteIndent(os, indent);
  os << "}";
}

std::string DiagnosticsPath(const TurbulenceConfig &cfg) {
  return cfg.basename + ".turb_init_diag.json";
}

void WriteProjectionDiagnostics(const TurbulenceConfig &cfg,
                                const ModeCatalog &catalog,
                                const VectorCoefficients &coeffs,
                                const VelocityStats &stats,
                                Real init_wall_seconds) {
  if (global_variable::my_rank != 0 || !cfg.dump_generator_diagnostics) return;

  std::ofstream os(DiagnosticsPath(cfg));
  os << std::setprecision(17);
  const std::vector<Real> velocity_shell =
      ComputeVectorShellEnergy(catalog, coeffs, cfg.nhigh);

  os << "{\n";
  WriteIndent(os, 2);
  os << "\"metadata\": ";
  WriteDiagnosticsMetadata(os, cfg,
                           cfg.projection_true_2d ? "projection_2d" : "projection_3d",
                           stats, init_wall_seconds, 2);
  os << ",\n";
  WriteIndent(os, 2);
  os << "\"catalogs\": {\n";
  WriteIndent(os, 4);
  os << "\"velocity\": ";
  WriteJsonModeCatalog(os, catalog, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"coefficients\": {\n";
  WriteIndent(os, 4);
  os << "\"velocity\": ";
  WriteJsonVectorCoefficients(os, coeffs, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"spectra\": {\n";
  WriteIndent(os, 4);
  os << "\"velocity_shell_energy\": ";
  WriteJsonNumericArray(os, velocity_shell);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"slices\": {}\n";
  os << "}\n";
}

void WriteStream2DDiagnostics(const TurbulenceConfig &cfg,
                              const ModeCatalog &catalog,
                              const ScalarCoefficients &psi_coeffs,
                              const VectorCoefficients &vel_coeffs,
                              const VelocityStats &stats,
                              Real init_wall_seconds) {
  if (global_variable::my_rank != 0 || !cfg.dump_generator_diagnostics) return;

  const std::vector<Real> psi_field =
      SampleScalarFieldXY(cfg, catalog, psi_coeffs, cfg.global_nx1, cfg.global_nx2,
                          UniformCellCenter(0, 1, cfg.x3min, cfg.x3max));
  const std::vector<Real> psi_shell =
      ComputeScalarShellEnergy(catalog, psi_coeffs, cfg.nhigh);
  const std::vector<Real> velocity_shell =
      ComputeVectorShellEnergy(catalog, vel_coeffs, cfg.nhigh);

  std::ofstream os(DiagnosticsPath(cfg));
  os << std::setprecision(17);
  os << "{\n";
  WriteIndent(os, 2);
  os << "\"metadata\": ";
  WriteDiagnosticsMetadata(os, cfg, "stream_2d", stats, init_wall_seconds, 2);
  os << ",\n";
  WriteIndent(os, 2);
  os << "\"catalogs\": {\n";
  WriteIndent(os, 4);
  os << "\"psi\": ";
  WriteJsonModeCatalog(os, catalog, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"coefficients\": {\n";
  WriteIndent(os, 4);
  os << "\"psi\": ";
  WriteJsonScalarCoefficients(os, psi_coeffs, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"velocity\": ";
  WriteJsonVectorCoefficients(os, vel_coeffs, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"spectra\": {\n";
  WriteIndent(os, 4);
  os << "\"psi_shell_energy\": ";
  WriteJsonNumericArray(os, psi_shell);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"velocity_shell_energy\": ";
  WriteJsonNumericArray(os, velocity_shell);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"slices\": {\n";
  WriteIndent(os, 4);
  os << "\"psi_xy\": ";
  WriteJsonField2D(os, psi_field, cfg.global_nx1, cfg.global_nx2,
                   UniformCellCenter(0, 1, cfg.x3min, cfg.x3max), 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "}\n";
  os << "}\n";
}

void WriteStream3DDiagnostics(const TurbulenceConfig &cfg,
                              const ModeCatalog &catalog1,
                              const ScalarCoefficients &phi1_coeffs,
                              const ModeCatalog &catalog2,
                              const ScalarCoefficients &phi2_coeffs,
                              const ModeCatalog &velocity_catalog,
                              const ShellFitResult &fit,
                              const ClebschRealization &realization,
                              const ClebschVelocityModes &velocity_modes,
                              const ClebschSectorResponse &sector_response,
                              const ClebschQualityGate &quality_gate,
                              int phi_nhigh,
                              const VelocityStats &stats,
                              Real init_wall_seconds) {
  if (global_variable::my_rank != 0 || !cfg.dump_generator_diagnostics) return;

  const int slice_k = cfg.global_nx3 / 2;
  const Real x3_slice = UniformCellCenter(slice_k, cfg.global_nx3, cfg.x3min, cfg.x3max);
  const std::vector<Real> phi1_xy =
      SampleScalarFieldXY(cfg, catalog1, phi1_coeffs, cfg.global_nx1, cfg.global_nx2, x3_slice);
  const std::vector<Real> phi2_xy =
      SampleScalarFieldXY(cfg, catalog2, phi2_coeffs, cfg.global_nx1, cfg.global_nx2, x3_slice);
  const std::vector<Real> phi1_shell =
      ComputeScalarShellEnergy(catalog1, phi1_coeffs, phi_nhigh);
  const std::vector<Real> phi2_shell =
      ComputeScalarShellEnergy(catalog2, phi2_coeffs, phi_nhigh);

  std::ofstream os(DiagnosticsPath(cfg));
  os << std::setprecision(17);
  os << "{\n";
  WriteIndent(os, 2);
  os << "\"metadata\": ";
  WriteDiagnosticsMetadata(os, cfg, "stream_3d", stats, init_wall_seconds, 2);
  os << ",\n";
  WriteIndent(os, 2);
  os << "\"catalogs\": {\n";
  WriteIndent(os, 4);
  os << "\"phi1\": ";
  WriteJsonModeCatalog(os, catalog1, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"phi2\": ";
  WriteJsonModeCatalog(os, catalog2, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"velocity\": ";
  WriteJsonModeCatalog(os, velocity_catalog, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"coefficients\": {\n";
  WriteIndent(os, 4);
  os << "\"phi1\": ";
  WriteJsonScalarCoefficients(os, phi1_coeffs, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"phi2\": ";
  WriteJsonScalarCoefficients(os, phi2_coeffs, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"spectra\": {\n";
  WriteIndent(os, 4);
  os << "\"phi_shell_max\": " << phi_nhigh << ",\n";
  WriteIndent(os, 4);
  os << "\"phi_shell_max_used\": " << phi_nhigh << ",\n";
  WriteIndent(os, 4);
  os << "\"phi1_shell_energy\": ";
  WriteJsonNumericArray(os, phi1_shell);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"phi2_shell_energy\": ";
  WriteJsonNumericArray(os, phi2_shell);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"shell_weights_phi1\": ";
  WriteJsonNumericArray(os, fit.shell_amp1);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"shell_weights_phi2\": ";
  WriteJsonNumericArray(os, fit.shell_amp2);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"predicted_velocity_shell_energy\": ";
  WriteJsonNumericArray(os, fit.predicted_shell_energy);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_retained_velocity_shell_energy\": ";
  WriteJsonNumericArray(os, realization.summary.shell_energy);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_retained_velocity_shell_energy_raw\": ";
  WriteJsonNumericArray(os, velocity_modes.shell_energy);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"fit\": {\n";
  WriteIndent(os, 4);
  os << "\"initial_loss\": " << fit.initial_loss << ",\n";
  WriteIndent(os, 4);
  os << "\"final_loss\": " << fit.final_loss << ",\n";
  WriteIndent(os, 4);
  os << "\"fitted_slope\": " << fit.fitted_slope << ",\n";
  WriteIndent(os, 4);
  os << "\"estimated_out_of_band_leakage\": " << fit.leakage_fraction << ",\n";
  WriteIndent(os, 4);
  os << "\"converged\": " << (fit.converged ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"seed_model_chosen\": \"" << ClebschSeedModelName(fit.chosen_seed_model)
     << "\",\n";
  WriteIndent(os, 4);
  os << "\"fit_starts\": ";
  WriteJsonClebschFitStarts(os, fit.fit_starts, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"trial_count\": " << realization.trial_count << ",\n";
  WriteIndent(os, 4);
  os << "\"quality_gate_pass\": " << (quality_gate.pass ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"quality_gate_forced_continue\": "
     << (quality_gate.forced_continue ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"chosen_trial\": " << realization.trial_index << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_loss\": " << realization.summary.loss << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_slope\": " << realization.summary.fitted_slope << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_out_of_band_leakage\": " << realization.summary.leakage_fraction << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"sector_response\": {\n";
  WriteIndent(os, 4);
  os << "\"raw_tensor\": ";
  WriteJsonClebschSectorBreakdown(os, sector_response.raw_shell_energy,
                                  sector_response.raw_total, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"fitted_response\": ";
  WriteJsonClebschSectorBreakdown(os, sector_response.fitted_shell_energy,
                                  sector_response.fitted_total, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"quality_gate\": {\n";
  WriteIndent(os, 4);
  os << "\"deterministic_slope_tolerance\": "
     << quality_gate.deterministic_slope_tolerance << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_slope_tolerance\": "
     << quality_gate.realized_slope_tolerance << ",\n";
  WriteIndent(os, 4);
  os << "\"leakage_tolerance\": " << quality_gate.leakage_tolerance << ",\n";
  WriteIndent(os, 4);
  os << "\"deterministic_converged\": "
     << (quality_gate.deterministic_converged ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"deterministic_slope_ok\": "
     << (quality_gate.deterministic_slope_ok ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"deterministic_leakage_ok\": "
     << (quality_gate.deterministic_leakage_ok ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"deterministic_pass\": "
     << (quality_gate.deterministic_pass ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_slope_ok\": "
     << (quality_gate.realized_slope_ok ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_leakage_ok\": "
     << (quality_gate.realized_leakage_ok ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"realized_pass\": "
     << (quality_gate.realized_pass ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"allow_poor_fit\": "
     << (quality_gate.allow_poor_fit ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"forced_continue\": "
     << (quality_gate.forced_continue ? "true" : "false") << ",\n";
  WriteIndent(os, 4);
  os << "\"pass\": " << (quality_gate.pass ? "true" : "false") << "\n";
  WriteIndent(os, 2);
  os << "},\n";
  WriteIndent(os, 2);
  os << "\"slices\": {\n";
  WriteIndent(os, 4);
  os << "\"phi1_xy_mid\": ";
  WriteJsonField2D(os, phi1_xy, cfg.global_nx1, cfg.global_nx2, x3_slice, 4);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"phi2_xy_mid\": ";
  WriteJsonField2D(os, phi2_xy, cfg.global_nx1, cfg.global_nx2, x3_slice, 4);
  os << "\n";
  WriteIndent(os, 2);
  os << "}\n";
  os << "}\n";
}

KOKKOS_INLINE_FUNCTION
Real SmoothTopHat(Real x, Real x0, Real x1, Real w) {
  if (w <= 0.0) {
    return (x >= x0 && x <= x1) ? 1.0 : 0.0;
  }
  return 0.5 * (tanh((x - x0) / w) - tanh((x - x1) / w));
}

KOKKOS_INLINE_FUNCTION
Real BistableReaction(Real theta, Real tau, Real t1, Real t2, Real t3) {
  const Real inv_tau = (tau != 0.0) ? 1.0 / tau : 0.0;
  return -inv_tau * (theta - t1) * (theta - t2) * (theta - t3);
}

KOKKOS_INLINE_FUNCTION
Real KppReaction(Real theta, Real tau, Real t1, Real t3) {
  const Real denom = t3 - t1;
  const Real inv = (tau != 0.0 && denom != 0.0) ? 1.0 / (tau * denom) : 0.0;
  return inv * (theta - t1) * (t3 - theta);
}

TurbulenceConfig ReadTurbulenceConfig(ParameterInput *pin, Mesh *pm,
                                      bool projection_true_2d) {
  TurbulenceConfig cfg;

  cfg.method = pin->GetOrAddBoolean("problem", "turb_use_stream_function", false)
                   ? VelocityMethod::StreamFunction
                   : VelocityMethod::Projection;
  cfg.v_rms = pin->GetOrAddReal("problem", "turb_v_rms", 1.0);
  cfg.nlow = pin->GetOrAddInteger("problem", "turb_nlow", 1);
  cfg.nhigh = pin->GetOrAddInteger("problem", "turb_nhigh", 4);
  cfg.expo = pin->GetOrAddReal("problem", "turb_expo", 5.0/3.0);
  cfg.sol_frac = pin->GetOrAddReal("problem", "turb_sol_frac", 1.0);
  cfg.rseed = GetTurbulenceSeed(pin);
  cfg.k_crit = pin->GetOrAddReal("problem", "turb_k_crit", 16.0);
  cfg.divfree_scalar_flux = pin->GetOrAddBoolean("problem", "divfree_scalar_flux", false);
  cfg.dump_generator_diagnostics =
      pin->GetOrAddBoolean("problem", "turb_dump_generator_diagnostics", false);
  cfg.stream_allow_poor_fit =
      pin->GetOrAddBoolean("problem", "turb_stream_allow_poor_fit", false);
  cfg.basename = pin->GetString("job", "basename");

  cfg.nx1 = pm->mb_indcs.nx1;
  cfg.nx2 = pm->mb_indcs.nx2;
  cfg.nx3 = pm->mb_indcs.nx3;
  cfg.mesh_multi_d = (cfg.nx2 > 1);
  cfg.mesh_three_d = (cfg.nx3 > 1);
  cfg.projection_true_2d = projection_true_2d;

  if (cfg.nlow < 1) {
    FatalProblemSetup("<problem>/turb_nlow must be >= 1.");
  }
  if (cfg.nhigh < cfg.nlow) {
    FatalProblemSetup("<problem>/turb_nhigh must be >= turb_nlow.");
  }
  if (cfg.k_crit <= 0.0) {
    FatalProblemSetup("<problem>/turb_k_crit must be > 0.");
  }

  cfg.stream_2d = (cfg.method == VelocityMethod::StreamFunction && cfg.nx3 == 1);
  if (cfg.method == VelocityMethod::StreamFunction) {
    if (cfg.nx2 == 1) {
      FatalProblemSetup("Stream-function initialization is only supported in 2D/3D.");
    }
    cfg.active_v3 = !cfg.stream_2d;
  } else {
    cfg.active_v3 = !projection_true_2d;
  }

  cfg.lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  cfg.ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  cfg.lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  cfg.x1min = pm->mesh_size.x1min;
  cfg.x1max = pm->mesh_size.x1max;
  cfg.x2min = pm->mesh_size.x2min;
  cfg.x2max = pm->mesh_size.x2max;
  cfg.x3min = pm->mesh_size.x3min;
  cfg.x3max = pm->mesh_size.x3max;
  cfg.global_nx1 = std::max(1, static_cast<int>(std::llround(cfg.lx / pm->mesh_size.dx1)));
  cfg.global_nx2 = (cfg.mesh_multi_d && pm->mesh_size.dx2 > 0.0)
                       ? std::max(1, static_cast<int>(std::llround(cfg.ly / pm->mesh_size.dx2)))
                       : 1;
  cfg.global_nx3 = (cfg.mesh_three_d && pm->mesh_size.dx3 > 0.0)
                       ? std::max(1, static_cast<int>(std::llround(cfg.lz / pm->mesh_size.dx3)))
                       : 1;
  cfg.dkx = 2.0*M_PI/cfg.lx;
  cfg.dky = (cfg.nx2 > 1) ? 2.0*M_PI/cfg.ly : 0.0;
  cfg.dkz = (cfg.nx3 > 1) ? 2.0*M_PI/cfg.lz : 0.0;
  cfg.k_crit_mag = cfg.k_crit * cfg.dkx;

  return cfg;
}

Real ModeAcceptanceProbability(const TurbulenceConfig &cfg, Real kiso) {
  if (kiso < cfg.k_crit_mag) return 1.0;
  return (cfg.k_crit_mag * cfg.k_crit_mag) / (kiso * kiso);
}

bool IsHalfSpaceMode(int nkx, int nky, int nkz) {
  if (nkx == 0) {
    if (nky < 0) return false;
    if (nky == 0 && nkz <= 0) return false;
  }
  return !(nkx == 0 && nky == 0 && nkz == 0);
}

ModeCatalog BuildModeCatalog(const TurbulenceConfig &cfg, int nlow, int nhigh,
                             RNG_State *rstate) {
  ModeCatalog catalog;
  const int nlow_sqr = nlow * nlow;
  const int nhigh_sqr = nhigh * nhigh;

  for (int nkx = 0; nkx <= nhigh; ++nkx) {
    for (int nky = (cfg.mesh_multi_d ? -nhigh : 0);
         nky <= (cfg.mesh_multi_d ? nhigh : 0); ++nky) {
      for (int nkz = (cfg.mesh_three_d ? -nhigh : 0);
           nkz <= (cfg.mesh_three_d ? nhigh : 0); ++nkz) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;

        // Keep half of the kx=0 plane to avoid double-counting the implicit conjugates.
        if (nkx == 0) {
          if (nky < 0) continue;
          if (nky == 0 && nkz <= 0) continue;
        }

        const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr < nlow_sqr || nsqr > nhigh_sqr) continue;

        ++catalog.total_modes;
        const Real kx = cfg.dkx * nkx;
        const Real ky = cfg.dky * nky;
        const Real kz = cfg.dkz * nkz;
        const Real kiso = std::sqrt(kx*kx + ky*ky + kz*kz);
        const Real p_accept = ModeAcceptanceProbability(cfg, kiso);

        if (RanSt(rstate) >= p_accept) continue;

        catalog.nkx.push_back(nkx);
        catalog.nky.push_back(nky);
        catalog.nkz.push_back(nkz);
        catalog.shell.push_back(ShellFromIndexSqr(nsqr));
        catalog.prob.push_back(p_accept);
        catalog.boost.push_back(1.0/std::sqrt(p_accept));
        catalog.kx.push_back(kx);
        catalog.ky.push_back(ky);
        catalog.kz.push_back(kz);
        catalog.kiso.push_back(kiso);
      }
    }
  }

  return catalog;
}

ModeCatalog BuildModeCatalog(const TurbulenceConfig &cfg, RNG_State *rstate) {
  return BuildModeCatalog(cfg, cfg.nlow, cfg.nhigh, rstate);
}

Real ProjectionAmplitudeNorm(const TurbulenceConfig &cfg, Real kiso) {
  const Real spectral_ndim = cfg.active_v3 ? 3.0 : 2.0;
  const Real norm_exp = 0.5 * (cfg.expo + spectral_ndim - 1.0);
  return (kiso > kTiny) ? 1.0/std::pow(kiso, norm_exp) : 0.0;
}

Real Stream2DPsiAmplitudeNorm(const TurbulenceConfig &cfg, Real kiso) {
  return (kiso > kTiny) ? 1.0/std::pow(kiso, 0.5 * (cfg.expo + 3.0)) : 0.0;
}

VectorCoefficients GenerateProjectionCoefficients(const TurbulenceConfig &cfg,
                                                  const ModeCatalog &catalog,
                                                  RNG_State *rstate) {
  VectorCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka0.resize(nmodes);
  coeffs.aka1.resize(nmodes);
  coeffs.aka2.resize(nmodes);
  coeffs.akb0.resize(nmodes);
  coeffs.akb1.resize(nmodes);
  coeffs.akb2.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real kiso = catalog.kiso[n];
    const Real norm = ProjectionAmplitudeNorm(cfg, kiso) * catalog.boost[n];
    const Real k_vec[3] = {catalog.kx[n], catalog.ky[n], catalog.kz[n]};
    Real a[3], b[3];

    for (int dir = 0; dir < 3; ++dir) {
      if (dir == 2 && !cfg.active_v3) {
        a[dir] = 0.0;
        b[dir] = 0.0;
      } else {
        a[dir] = norm * RanGaussianSt(rstate);
        b[dir] = norm * RanGaussianSt(rstate);
      }
    }

    const Real kiso_sqr = kiso * kiso;
    const Real k_dot_a = k_vec[0]*a[0] + k_vec[1]*a[1] + k_vec[2]*a[2];
    const Real k_dot_b = k_vec[0]*b[0] + k_vec[1]*b[1] + k_vec[2]*b[2];

    for (int dir = 0; dir < 3; ++dir) {
      const Real a_sol = a[dir] - k_vec[dir] * k_dot_a / kiso_sqr;
      const Real b_sol = b[dir] - k_vec[dir] * k_dot_b / kiso_sqr;
      const Real a_div = k_vec[dir] * k_dot_a / kiso_sqr;
      const Real b_div = k_vec[dir] * k_dot_b / kiso_sqr;
      a[dir] = cfg.sol_frac * a_sol + (1.0 - cfg.sol_frac) * a_div;
      b[dir] = cfg.sol_frac * b_sol + (1.0 - cfg.sol_frac) * b_div;
    }

    coeffs.aka0[n] = a[0];
    coeffs.aka1[n] = a[1];
    coeffs.aka2[n] = a[2];
    coeffs.akb0[n] = b[0];
    coeffs.akb1[n] = b[1];
    coeffs.akb2[n] = b[2];
  }

  return coeffs;
}

ScalarCoefficients GenerateStream2DPsiCoefficients(const TurbulenceConfig &cfg,
                                                   const ModeCatalog &catalog,
                                                   RNG_State *rstate) {
  ScalarCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka.resize(nmodes);
  coeffs.akb.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real norm = Stream2DPsiAmplitudeNorm(cfg, catalog.kiso[n]) * catalog.boost[n];
    coeffs.aka[n] = norm * RanGaussianSt(rstate);
    coeffs.akb[n] = norm * RanGaussianSt(rstate);
  }

  return coeffs;
}

ScalarCoefficients GenerateShellWeightedScalarCoefficients(const ModeCatalog &catalog,
                                                           const std::vector<Real> &shell_amp,
                                                           RNG_State *rstate) {
  ScalarCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka.resize(nmodes);
  coeffs.akb.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const int shell = catalog.shell[n];
    const Real norm = shell_amp[shell] * catalog.boost[n];
    coeffs.aka[n] = norm * RanGaussianSt(rstate);
    coeffs.akb[n] = norm * RanGaussianSt(rstate);
  }

  return coeffs;
}

VectorCoefficients ConvertPsiToVelocityCoefficients(const ModeCatalog &catalog,
                                                    const ScalarCoefficients &psi_coeffs) {
  VectorCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka0.resize(nmodes);
  coeffs.aka1.resize(nmodes);
  coeffs.aka2.assign(nmodes, 0.0);
  coeffs.akb0.resize(nmodes);
  coeffs.akb1.resize(nmodes);
  coeffs.akb2.assign(nmodes, 0.0);

  for (int n = 0; n < nmodes; ++n) {
    const Real kx = catalog.kx[n];
    const Real ky = catalog.ky[n];
    const Real a = psi_coeffs.aka[n];
    const Real b = psi_coeffs.akb[n];

    // psi = a cos(k·x) - b sin(k·x)
    // u = (dpsi/dy, -dpsi/dx, 0) stays in the same cosine/sine basis.
    coeffs.aka0[n] = -ky * b;
    coeffs.akb0[n] =  ky * a;
    coeffs.aka1[n] =  kx * b;
    coeffs.akb1[n] = -kx * a;
  }

  return coeffs;
}

VectorCoefficients ConvertVelocityToVectorPotentialCoefficients(
    const ModeCatalog &catalog, const VectorCoefficients &vel_coeffs) {
  VectorCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka0.resize(nmodes);
  coeffs.aka1.resize(nmodes);
  coeffs.aka2.resize(nmodes);
  coeffs.akb0.resize(nmodes);
  coeffs.akb1.resize(nmodes);
  coeffs.akb2.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real kx = catalog.kx[n];
    const Real ky = catalog.ky[n];
    const Real kz = catalog.kz[n];
    const Real k2 = kx*kx + ky*ky + kz*kz;
    if (k2 <= kTiny) {
      coeffs.aka0[n] = 0.0;
      coeffs.aka1[n] = 0.0;
      coeffs.aka2[n] = 0.0;
      coeffs.akb0[n] = 0.0;
      coeffs.akb1[n] = 0.0;
      coeffs.akb2[n] = 0.0;
      continue;
    }

    coeffs.aka0[n] = -(ky*vel_coeffs.akb2[n] - kz*vel_coeffs.akb1[n]) / k2;
    coeffs.aka1[n] = -(kz*vel_coeffs.akb0[n] - kx*vel_coeffs.akb2[n]) / k2;
    coeffs.aka2[n] = -(kx*vel_coeffs.akb1[n] - ky*vel_coeffs.akb0[n]) / k2;
    coeffs.akb0[n] =  (ky*vel_coeffs.aka2[n] - kz*vel_coeffs.aka1[n]) / k2;
    coeffs.akb1[n] =  (kz*vel_coeffs.aka0[n] - kx*vel_coeffs.aka2[n]) / k2;
    coeffs.akb2[n] =  (kx*vel_coeffs.aka1[n] - ky*vel_coeffs.aka0[n]) / k2;
  }

  return coeffs;
}

KOKKOS_INLINE_FUNCTION
void AccumulateVelocityFromVectorPotential(
    const DvceArray1D<Real> &d_kx, const DvceArray1D<Real> &d_ky,
    const DvceArray1D<Real> &d_kz, const DvceArray1D<Real> &d_aA0,
    const DvceArray1D<Real> &d_aA1, const DvceArray1D<Real> &d_aA2,
    const DvceArray1D<Real> &d_bA0, const DvceArray1D<Real> &d_bA1,
    const DvceArray1D<Real> &d_bA2, int mode_count, Real x1v, Real x2v, Real x3v,
    Real dx1, Real dx2, Real dx3, bool multi_d, bool three_d,
    Real &vx, Real &vy, Real &vz) {
  for (int n = 0; n < mode_count; ++n) {
    const Real kx = d_kx(n);
    const Real ky = d_ky(n);
    const Real kz = d_kz(n);

    const Real kx_eff = (dx1 > 0.0) ? (2.0/dx1) * sin(0.5*kx*dx1) : 0.0;
    const Real ky_eff = (multi_d && dx2 > 0.0) ? (2.0/dx2) * sin(0.5*ky*dx2) : 0.0;
    const Real kz_eff = (three_d && dx3 > 0.0) ? (2.0/dx3) * sin(0.5*kz*dx3) : 0.0;
    const Real keff2 = kx_eff*kx_eff + ky_eff*ky_eff + kz_eff*kz_eff;
    if (keff2 <= kTiny) continue;

    Real aAx = d_aA0(n);
    Real aAy = d_aA1(n);
    Real aAz = d_aA2(n);
    Real bAx = d_bA0(n);
    Real bAy = d_bA1(n);
    Real bAz = d_bA2(n);

    const Real inv_keff2 = 1.0 / keff2;
    const Real k_dot_a = kx_eff*aAx + ky_eff*aAy + kz_eff*aAz;
    const Real k_dot_b = kx_eff*bAx + ky_eff*bAy + kz_eff*bAz;
    aAx -= kx_eff * k_dot_a * inv_keff2;
    aAy -= ky_eff * k_dot_a * inv_keff2;
    aAz -= kz_eff * k_dot_a * inv_keff2;
    bAx -= kx_eff * k_dot_b * inv_keff2;
    bAy -= ky_eff * k_dot_b * inv_keff2;
    bAz -= kz_eff * k_dot_b * inv_keff2;

    const Real aVx = -(ky_eff*bAz - kz_eff*bAy);
    const Real aVy = -(kz_eff*bAx - kx_eff*bAz);
    const Real aVz = -(kx_eff*bAy - ky_eff*bAx);
    const Real bVx =  (ky_eff*aAz - kz_eff*aAy);
    const Real bVy =  (kz_eff*aAx - kx_eff*aAz);
    const Real bVz =  (kx_eff*aAy - ky_eff*aAx);

    const Real phase = kx*x1v + ky*x2v + kz*x3v;
    const Real cosk = cos(phase);
    const Real sink = sin(phase);

    vx += aVx*cosk - bVx*sink;
    vy += aVy*cosk - bVy*sink;
    vz += aVz*cosk - bVz*sink;
  }
}

void SynthesizeVelocityFromVectorModes(MeshBlockPack *pmbp, const TurbulenceConfig &cfg,
                                       Real den, const ModeCatalog &catalog,
                                       const VectorCoefficients &coeffs) {
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;

  const int nmodes = catalog.KeptModes();
  if (nmodes == 0) return;

  DvceArray1D<Real> d_kx("scalar_mix_kx", nmodes);
  DvceArray1D<Real> d_ky("scalar_mix_ky", nmodes);
  DvceArray1D<Real> d_kz("scalar_mix_kz", nmodes);
  DvceArray1D<Real> d_aka0("scalar_mix_aka0", nmodes);
  DvceArray1D<Real> d_aka1("scalar_mix_aka1", nmodes);
  DvceArray1D<Real> d_aka2("scalar_mix_aka2", nmodes);
  DvceArray1D<Real> d_akb0("scalar_mix_akb0", nmodes);
  DvceArray1D<Real> d_akb1("scalar_mix_akb1", nmodes);
  DvceArray1D<Real> d_akb2("scalar_mix_akb2", nmodes);

  auto h_kx = Kokkos::create_mirror_view(d_kx);
  auto h_ky = Kokkos::create_mirror_view(d_ky);
  auto h_kz = Kokkos::create_mirror_view(d_kz);
  auto h_aka0 = Kokkos::create_mirror_view(d_aka0);
  auto h_aka1 = Kokkos::create_mirror_view(d_aka1);
  auto h_aka2 = Kokkos::create_mirror_view(d_aka2);
  auto h_akb0 = Kokkos::create_mirror_view(d_akb0);
  auto h_akb1 = Kokkos::create_mirror_view(d_akb1);
  auto h_akb2 = Kokkos::create_mirror_view(d_akb2);

  for (int n = 0; n < nmodes; ++n) {
    h_kx(n) = catalog.kx[n];
    h_ky(n) = catalog.ky[n];
    h_kz(n) = catalog.kz[n];
    h_aka0(n) = coeffs.aka0[n];
    h_aka1(n) = coeffs.aka1[n];
    h_aka2(n) = coeffs.aka2[n];
    h_akb0(n) = coeffs.akb0[n];
    h_akb1(n) = coeffs.akb1[n];
    h_akb2(n) = coeffs.akb2[n];
  }

  Kokkos::deep_copy(d_kx, h_kx);
  Kokkos::deep_copy(d_ky, h_ky);
  Kokkos::deep_copy(d_kz, h_kz);
  Kokkos::deep_copy(d_aka0, h_aka0);
  Kokkos::deep_copy(d_aka1, h_aka1);
  Kokkos::deep_copy(d_aka2, h_aka2);
  Kokkos::deep_copy(d_akb0, h_akb0);
  Kokkos::deep_copy(d_akb1, h_akb1);
  Kokkos::deep_copy(d_akb2, h_akb2);

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  const bool active_v3 = cfg.active_v3;

  par_for("scalar_mix_turb_vel_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    for (int n = 0; n < nmodes; ++n) {
      const Real phase = d_kx(n)*x1v + d_ky(n)*x2v + d_kz(n)*x3v;
      const Real cosk = cos(phase);
      const Real sink = sin(phase);
      vx += d_aka0(n)*cosk - d_akb0(n)*sink;
      vy += d_aka1(n)*cosk - d_akb1(n)*sink;
      vz += d_aka2(n)*cosk - d_akb2(n)*sink;
    }

    u0(m,IM1,k,j,i) = den * vx;
    u0(m,IM2,k,j,i) = den * vy;
    u0(m,IM3,k,j,i) = active_v3 ? den * vz : 0.0;
  });
}

VelocityStats NormalizeVelocityField(MeshBlockPack *pmbp, const TurbulenceConfig &cfg) {
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  const bool active_v3 = cfg.active_v3;

  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  Real sum_mass = 0.0;
  Real sum_mom1 = 0.0;
  Real sum_mom2 = 0.0;
  Real sum_mom3 = 0.0;

  Kokkos::parallel_reduce("scalar_mix_turb_vel_mean",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &s_mass, Real &s_mom1, Real &s_mom2, Real &s_mom3) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    const Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    const Real rho = u0(m,IDN,k,j,i);
    s_mass += rho * vol;
    s_mom1 += u0(m,IM1,k,j,i) * vol;
    s_mom2 += u0(m,IM2,k,j,i) * vol;
    s_mom3 += u0(m,IM3,k,j,i) * vol;
  }, Kokkos::Sum<Real>(sum_mass), Kokkos::Sum<Real>(sum_mom1),
     Kokkos::Sum<Real>(sum_mom2), Kokkos::Sum<Real>(sum_mom3));

#if MPI_PARALLEL_ENABLED
  Real local_mean[4] = {sum_mass, sum_mom1, sum_mom2, sum_mom3};
  Real global_mean[4];
  MPI_Allreduce(local_mean, global_mean, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  sum_mass = global_mean[0];
  sum_mom1 = global_mean[1];
  sum_mom2 = global_mean[2];
  sum_mom3 = global_mean[3];
#endif

  const Real vmean1 = sum_mom1 / sum_mass;
  const Real vmean2 = sum_mom2 / sum_mass;
  const Real vmean3 = active_v3 ? sum_mom3 / sum_mass : 0.0;

  par_for("scalar_mix_turb_vel_submean", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real rho = u0(m,IDN,k,j,i);
    u0(m,IM1,k,j,i) -= rho * vmean1;
    u0(m,IM2,k,j,i) -= rho * vmean2;
    u0(m,IM3,k,j,i) -= rho * vmean3;
  });

  Real sum_v1sqr = 0.0;
  Real sum_v2sqr = 0.0;
  Real sum_v3sqr = 0.0;
  Real sum_vol = 0.0;

  Kokkos::parallel_reduce("scalar_mix_turb_vel_rms",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &s_v1sqr, Real &s_v2sqr, Real &s_v3sqr, Real &s_vol) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    const Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    const Real rho = u0(m,IDN,k,j,i);
    const Real v1 = u0(m,IM1,k,j,i)/rho;
    const Real v2 = u0(m,IM2,k,j,i)/rho;
    const Real v3 = active_v3 ? u0(m,IM3,k,j,i)/rho : 0.0;

    s_v1sqr += v1*v1 * vol;
    s_v2sqr += v2*v2 * vol;
    s_v3sqr += v3*v3 * vol;
    s_vol += vol;
  }, Kokkos::Sum<Real>(sum_v1sqr), Kokkos::Sum<Real>(sum_v2sqr),
     Kokkos::Sum<Real>(sum_v3sqr), Kokkos::Sum<Real>(sum_vol));

#if MPI_PARALLEL_ENABLED
  Real local_rms[4] = {sum_v1sqr, sum_v2sqr, sum_v3sqr, sum_vol};
  Real global_rms[4];
  MPI_Allreduce(local_rms, global_rms, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  sum_v1sqr = global_rms[0];
  sum_v2sqr = global_rms[1];
  sum_v3sqr = global_rms[2];
  sum_vol = global_rms[3];
#endif

  const Real v1_rms = std::sqrt(sum_v1sqr / sum_vol);
  const Real v2_rms = std::sqrt(sum_v2sqr / sum_vol);
  const Real v3_rms = std::sqrt(sum_v3sqr / sum_vol);
  const Real current_vrms = std::sqrt(v1_rms*v1_rms + v2_rms*v2_rms + v3_rms*v3_rms);
  const Real scale = (current_vrms > kTiny) ? cfg.v_rms / current_vrms : 0.0;

  if (global_variable::my_rank == 0) {
    std::cout << "  vmean = (" << vmean1 << ", " << vmean2 << ", " << vmean3 << ")"
              << std::endl
              << "  component RMS before: (" << v1_rms << ", " << v2_rms << ", "
              << v3_rms << "), total = " << current_vrms << std::endl
              << "  uniform scale factor = " << scale
              << " (target v_rms = " << cfg.v_rms << ")" << std::endl
              << "  component RMS after: (" << v1_rms*scale << ", " << v2_rms*scale
              << ", " << v3_rms*scale << ")" << std::endl;
  }

  par_for("scalar_mix_turb_vel_normalize", DevExeSpace(), 0, nmb-1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IM1,k,j,i) *= scale;
    u0(m,IM2,k,j,i) *= scale;
    u0(m,IM3,k,j,i) = active_v3 ? u0(m,IM3,k,j,i) * scale : 0.0;
  });

  return {vmean1, vmean2, vmean3, scale};
}

void PopulateScalarFaceVelocitiesFromVectorPotential(
    MeshBlockPack *pmbp, const TurbulenceConfig &cfg, const ModeCatalog &catalog,
    const VectorCoefficients &apot_coeffs, const VelocityStats &stats) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int ng = indcs.ng;
  const int nmb = pmbp->nmb_thispack;
  const bool multi_d = cfg.mesh_multi_d;
  const bool three_d = cfg.mesh_three_d;
  const int nmodes = catalog.KeptModes();

  if (nmodes == 0) return;

  auto &hydro = *pmbp->phydro;
  const int ncells1 = nx1 + 2*ng;
  const int ncells2 = multi_d ? nx2 + 2*ng : 1;
  const int ncells3 = three_d ? nx3 + 2*ng : 1;
  if (!hydro.scalar_vface) {
    hydro.scalar_vface = std::make_unique<DvceFaceFld4D<Real>>(
        "scalar_vface", nmb, ncells3, ncells2, ncells1);
  }
  hydro.use_scalar_face_velocity = true;

  DvceArray1D<Real> d_kx("scalar_face_kx", nmodes);
  DvceArray1D<Real> d_ky("scalar_face_ky", nmodes);
  DvceArray1D<Real> d_kz("scalar_face_kz", nmodes);
  DvceArray1D<Real> d_aA0("scalar_face_aA0", nmodes);
  DvceArray1D<Real> d_aA1("scalar_face_aA1", nmodes);
  DvceArray1D<Real> d_aA2("scalar_face_aA2", nmodes);
  DvceArray1D<Real> d_bA0("scalar_face_bA0", nmodes);
  DvceArray1D<Real> d_bA1("scalar_face_bA1", nmodes);
  DvceArray1D<Real> d_bA2("scalar_face_bA2", nmodes);

  auto h_kx = Kokkos::create_mirror_view(d_kx);
  auto h_ky = Kokkos::create_mirror_view(d_ky);
  auto h_kz = Kokkos::create_mirror_view(d_kz);
  auto h_aA0 = Kokkos::create_mirror_view(d_aA0);
  auto h_aA1 = Kokkos::create_mirror_view(d_aA1);
  auto h_aA2 = Kokkos::create_mirror_view(d_aA2);
  auto h_bA0 = Kokkos::create_mirror_view(d_bA0);
  auto h_bA1 = Kokkos::create_mirror_view(d_bA1);
  auto h_bA2 = Kokkos::create_mirror_view(d_bA2);

  for (int n = 0; n < nmodes; ++n) {
    h_kx(n) = catalog.kx[n];
    h_ky(n) = catalog.ky[n];
    h_kz(n) = catalog.kz[n];
    h_aA0(n) = apot_coeffs.aka0[n];
    h_aA1(n) = apot_coeffs.aka1[n];
    h_aA2(n) = apot_coeffs.aka2[n];
    h_bA0(n) = apot_coeffs.akb0[n];
    h_bA1(n) = apot_coeffs.akb1[n];
    h_bA2(n) = apot_coeffs.akb2[n];
  }

  Kokkos::deep_copy(d_kx, h_kx);
  Kokkos::deep_copy(d_ky, h_ky);
  Kokkos::deep_copy(d_kz, h_kz);
  Kokkos::deep_copy(d_aA0, h_aA0);
  Kokkos::deep_copy(d_aA1, h_aA1);
  Kokkos::deep_copy(d_aA2, h_aA2);
  Kokkos::deep_copy(d_bA0, h_bA0);
  Kokkos::deep_copy(d_bA1, h_bA1);
  Kokkos::deep_copy(d_bA2, h_bA2);

  auto &size = pmbp->pmb->mb_size;
  auto sface_x1f = hydro.scalar_vface->x1f;
  auto sface_x2f = hydro.scalar_vface->x2f;
  auto sface_x3f = hydro.scalar_vface->x3f;

  int il = is - ng;
  int iu = ie + ng + 1;
  int jl = multi_d ? js - ng : js;
  int ju = multi_d ? je + ng : je;
  int kl = three_d ? ks - ng : ks;
  int ku = three_d ? ke + ng : ke;
  par_for("scalar_face_x1", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    const Real x1v = LeftEdgeX(i-is, nx1, x1min, x1max);

    Real x2v = 0.0;
    if (multi_d) {
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      x2v = CellCenterX(j-js, nx2, x2min, x2max);
    }

    Real x3v = 0.0;
    if (three_d) {
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    }

    Real vx = 0.0;
    Real vy = 0.0;
    Real vz = 0.0;
    AccumulateVelocityFromVectorPotential(
        d_kx, d_ky, d_kz, d_aA0, d_aA1, d_aA2, d_bA0, d_bA1, d_bA2, nmodes,
        x1v, x2v, x3v, size.d_view(m).dx1, size.d_view(m).dx2, size.d_view(m).dx3,
        multi_d, three_d, vx, vy, vz);
    sface_x1f(m,k,j,i) = (vx - stats.vmean1) * stats.scale;
  });

  if (multi_d) {
    il = is - ng;
    iu = ie + ng;
    jl = js - ng;
    ju = je + ng + 1;
    kl = three_d ? ks - ng : ks;
    ku = three_d ? ke + ng : ke;
    par_for("scalar_face_x2", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      const Real x2v = LeftEdgeX(j-js, nx2, x2min, x2max);

      Real x3v = 0.0;
      if (three_d) {
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }

      Real vx = 0.0;
      Real vy = 0.0;
      Real vz = 0.0;
      AccumulateVelocityFromVectorPotential(
          d_kx, d_ky, d_kz, d_aA0, d_aA1, d_aA2, d_bA0, d_bA1, d_bA2, nmodes,
          x1v, x2v, x3v, size.d_view(m).dx1, size.d_view(m).dx2, size.d_view(m).dx3,
          multi_d, three_d, vx, vy, vz);
      sface_x2f(m,k,j,i) = (vy - stats.vmean2) * stats.scale;
    });
  }

  if (three_d) {
    il = is - ng;
    iu = ie + ng;
    jl = js - ng;
    ju = je + ng;
    kl = ks - ng;
    ku = ke + ng + 1;
    par_for("scalar_face_x3", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      if (!cfg.active_v3) {
        sface_x3f(m,k,j,i) = 0.0;
        return;
      }

      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      const Real x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);

      Real vx = 0.0;
      Real vy = 0.0;
      Real vz = 0.0;
      AccumulateVelocityFromVectorPotential(
          d_kx, d_ky, d_kz, d_aA0, d_aA1, d_aA2, d_bA0, d_bA1, d_bA2, nmodes,
          x1v, x2v, x3v, size.d_view(m).dx1, size.d_view(m).dx2, size.d_view(m).dx3,
          multi_d, three_d, vx, vy, vz);
      sface_x3f(m,k,j,i) = (vz - stats.vmean3) * stats.scale;
    });
  }
}

std::vector<int> CountModesPerShell(const ModeCatalog &catalog, int nhigh) {
  std::vector<int> counts(nhigh + 1, 0);
  for (int shell : catalog.shell) {
    if (shell >= 1 && shell <= nhigh) ++counts[shell];
  }
  return counts;
}

void EnsureShellCoverage(const std::vector<int> &counts, int shell_lo, int shell_hi,
                         const std::string &label) {
  for (int shell = shell_lo; shell <= shell_hi; ++shell) {
    if (counts[shell] == 0) {
      FatalProblemSetup("Stream-function mode has no accepted modes in shell "
                        + std::to_string(shell) + " for " + label + ".");
    }
  }
}

std::vector<SignedMode> BuildSignedModes(const ModeCatalog &catalog) {
  std::vector<SignedMode> signed_modes;
  signed_modes.reserve(2 * catalog.KeptModes());

  for (int n = 0; n < catalog.KeptModes(); ++n) {
    const Real variance_weight = 0.5 * catalog.boost[n] * catalog.boost[n];
    signed_modes.push_back({catalog.nkx[n], catalog.nky[n], catalog.nkz[n], catalog.shell[n],
                            catalog.kx[n], catalog.ky[n], catalog.kz[n], variance_weight});
    signed_modes.push_back({-catalog.nkx[n], -catalog.nky[n], -catalog.nkz[n],
                            catalog.shell[n], -catalog.kx[n], -catalog.ky[n],
                            -catalog.kz[n], variance_weight});
  }

  return signed_modes;
}

std::vector<SignedScalarCoeffMode> BuildSignedScalarCoeffModes(
    const ModeCatalog &catalog, const ScalarCoefficients &coeffs) {
  std::vector<SignedScalarCoeffMode> signed_modes;
  signed_modes.reserve(2 * catalog.KeptModes());

  for (int n = 0; n < catalog.KeptModes(); ++n) {
    const std::complex<Real> ck(0.5 * coeffs.aka[n], 0.5 * coeffs.akb[n]);
    signed_modes.push_back({catalog.nkx[n], catalog.nky[n], catalog.nkz[n], catalog.shell[n],
                            catalog.kx[n], catalog.ky[n], catalog.kz[n], ck});
    signed_modes.push_back({-catalog.nkx[n], -catalog.nky[n], -catalog.nkz[n],
                            catalog.shell[n], -catalog.kx[n], -catalog.ky[n],
                            -catalog.kz[n], std::conj(ck)});
  }

  return signed_modes;
}

Real FitLogSlope(const std::vector<Real> &spectrum, int nlow, int nhigh);
Real LeakageFraction(const std::vector<Real> &spectrum, const TurbulenceConfig &cfg,
                     int nkout);

Real SpectrumShapeLoss(const std::vector<Real> &spectrum, const TurbulenceConfig &cfg,
                       int nkout) {
  const Real leak_ref_low = TargetShellEnergy(std::max(cfg.nlow, 1), cfg.expo);
  const Real leak_ref_high = TargetShellEnergy(cfg.nhigh, cfg.expo);
  const Real smooth_weight = 1.0e-2;
  Real loss = 0.0;

  for (int shell = 1; shell <= nkout; ++shell) {
    const Real pred = std::max(spectrum[shell], kTiny);
    if (shell >= cfg.nlow && shell <= cfg.nhigh) {
      const Real target = TargetShellEnergy(shell, cfg.expo);
      const Real resid = std::log(pred) - std::log(target);
      loss += resid * resid;
    } else {
      const bool low_leak = (shell < cfg.nlow);
      const Real leak_weight = low_leak ? 4.0 : 6.0;
      const Real leak_ref = low_leak ? leak_ref_low : leak_ref_high;
      const Real scaled = pred / leak_ref;
      loss += leak_weight * scaled * scaled;
    }
  }

  for (int shell = 2; shell <= nkout; ++shell) {
    const Real prev = std::max(spectrum[shell - 1], kTiny);
    const Real curr = std::max(spectrum[shell], kTiny);
    const Real diff = std::log(curr) - std::log(prev);
    loss += smooth_weight * diff * diff;
  }

  return loss;
}

RealizedSpectrumSummary SummarizeRealizedSpectrum(std::vector<Real> shell_energy,
                                                 const TurbulenceConfig &cfg,
                                                 int nkout) {
  RealizedSpectrumSummary summary;
  summary.shell_energy = std::move(shell_energy);

  Real pred_mean = 0.0;
  Real target_mean = 0.0;
  int nband = 0;
  for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
    pred_mean += summary.shell_energy[shell];
    target_mean += TargetShellEnergy(shell, cfg.expo);
    ++nband;
  }
  if (pred_mean > kTiny && nband > 0) {
    const Real scale = target_mean / pred_mean;
    for (int shell = 1; shell <= nkout; ++shell) summary.shell_energy[shell] *= scale;
  }

  summary.loss = SpectrumShapeLoss(summary.shell_energy, cfg, nkout);
  summary.leakage_fraction = LeakageFraction(summary.shell_energy, cfg, nkout);
  summary.fitted_slope = FitLogSlope(summary.shell_energy, cfg.nlow, cfg.nhigh);
  return summary;
}

std::uint64_t EncodeSignedModeKey(int nkx, int nky, int nkz, int offset, int span) {
  const std::uint64_t x = static_cast<std::uint64_t>(nkx + offset);
  const std::uint64_t y = static_cast<std::uint64_t>(nky + offset);
  const std::uint64_t z = static_cast<std::uint64_t>(nkz + offset);
  const std::uint64_t stride = static_cast<std::uint64_t>(span);
  return (x * stride + y) * stride + z;
}

bool SignedModeInBounds(int nkx, int nky, int nkz, int max_n) {
  return (nkx >= -max_n && nkx <= max_n &&
          nky >= -max_n && nky <= max_n &&
          nkz >= -max_n && nkz <= max_n);
}

std::unordered_map<std::uint64_t, int> BuildSignedModeLookup(
    const std::vector<SignedMode> &signed_modes, int max_n) {
  const int offset = max_n;
  const int span = 2 * offset + 1;
  std::unordered_map<std::uint64_t, int> lookup;
  lookup.reserve(2 * signed_modes.size());
  for (int n = 0; n < static_cast<int>(signed_modes.size()); ++n) {
    lookup.emplace(EncodeSignedModeKey(signed_modes[n].nkx, signed_modes[n].nky,
                                       signed_modes[n].nkz, offset, span), n);
  }
  return lookup;
}

std::unordered_map<std::uint64_t, int> BuildSignedModeLookup(
    const std::vector<SignedScalarCoeffMode> &signed_modes, int max_n) {
  const int offset = max_n;
  const int span = 2 * offset + 1;
  std::unordered_map<std::uint64_t, int> lookup;
  lookup.reserve(2 * signed_modes.size());
  for (int n = 0; n < static_cast<int>(signed_modes.size()); ++n) {
    lookup.emplace(EncodeSignedModeKey(signed_modes[n].nkx, signed_modes[n].nky,
                                       signed_modes[n].nkz, offset, span), n);
  }
  return lookup;
}

ClebschVelocityModes BuildClebschVelocityModes(
    const ModeCatalog &velocity_catalog, const ModeCatalog &catalog1,
    const ScalarCoefficients &phi1_coeffs, const ModeCatalog &catalog2,
    const ScalarCoefficients &phi2_coeffs, int npot_shells, int nkout) {
  ClebschVelocityModes result;
  const int nmodes = velocity_catalog.KeptModes();
  result.coeffs.aka0.assign(nmodes, 0.0);
  result.coeffs.aka1.assign(nmodes, 0.0);
  result.coeffs.aka2.assign(nmodes, 0.0);
  result.coeffs.akb0.assign(nmodes, 0.0);
  result.coeffs.akb1.assign(nmodes, 0.0);
  result.coeffs.akb2.assign(nmodes, 0.0);
  result.shell_energy.assign(nkout + 1, 0.0);

  const auto signed1 = BuildSignedScalarCoeffModes(catalog1, phi1_coeffs);
  const auto signed2 = BuildSignedScalarCoeffModes(catalog2, phi2_coeffs);
  const auto lookup2 = BuildSignedModeLookup(signed2, npot_shells);
  const int offset = npot_shells;
  const int span = 2 * offset + 1;

  for (int n = 0; n < nmodes; ++n) {
    const int rx = velocity_catalog.nkx[n];
    const int ry = velocity_catalog.nky[n];
    const int rz = velocity_catalog.nkz[n];
    std::array<std::complex<Real>, 3> accum = {
        std::complex<Real>(0.0, 0.0), std::complex<Real>(0.0, 0.0),
        std::complex<Real>(0.0, 0.0)};

    for (const auto &mode_p : signed1) {
      const int qx = rx - mode_p.nkx;
      const int qy = ry - mode_p.nky;
      const int qz = rz - mode_p.nkz;
      if (!SignedModeInBounds(qx, qy, qz, npot_shells)) continue;

      const auto it = lookup2.find(EncodeSignedModeKey(qx, qy, qz, offset, span));
      if (it == lookup2.end()) continue;
      const auto &mode_q = signed2[it->second];

      const Real cx = mode_p.ky*mode_q.kz - mode_p.kz*mode_q.ky;
      const Real cy = mode_p.kz*mode_q.kx - mode_p.kx*mode_q.kz;
      const Real cz = mode_p.kx*mode_q.ky - mode_p.ky*mode_q.kx;
      if ((cx*cx + cy*cy + cz*cz) <= kTiny) continue;

      const std::complex<Real> amp = -(mode_p.coeff * mode_q.coeff);
      accum[0] += amp * cx;
      accum[1] += amp * cy;
      accum[2] += amp * cz;
    }

    const std::complex<Real> ux = velocity_catalog.boost[n] * accum[0];
    const std::complex<Real> uy = velocity_catalog.boost[n] * accum[1];
    const std::complex<Real> uz = velocity_catalog.boost[n] * accum[2];
    result.coeffs.aka0[n] = 2.0 * ux.real();
    result.coeffs.akb0[n] = 2.0 * ux.imag();
    result.coeffs.aka1[n] = 2.0 * uy.real();
    result.coeffs.akb1[n] = 2.0 * uy.imag();
    result.coeffs.aka2[n] = 2.0 * uz.real();
    result.coeffs.akb2[n] = 2.0 * uz.imag();
    result.shell_energy[velocity_catalog.shell[n]] +=
        2.0 * (std::norm(ux) + std::norm(uy) + std::norm(uz));
  }

  return result;
}

RealizedSpectrumSummary EvaluateClebschRealizationSpectrum(
    const ModeCatalog &velocity_catalog, const ModeCatalog &catalog1,
    const ScalarCoefficients &phi1_coeffs, const ModeCatalog &catalog2,
    const ScalarCoefficients &phi2_coeffs, const TurbulenceConfig &cfg,
    int npot_shells) {
  const ClebschVelocityModes velocity_modes =
      BuildClebschVelocityModes(velocity_catalog, catalog1, phi1_coeffs,
                                catalog2, phi2_coeffs, npot_shells, cfg.nhigh);
  return SummarizeRealizedSpectrum(velocity_modes.shell_energy, cfg, cfg.nhigh);
}

int ClebschRealizationTrials(const TurbulenceConfig &cfg) {
  if (cfg.k_crit >= cfg.nhigh) return 1;
  if (cfg.k_crit >= 8.0) return 2;
  if (cfg.k_crit >= 4.0) return 4;
  return 8;
}

ClebschRealization SelectClebschRealization(const TurbulenceConfig &cfg,
                                            const ModeCatalog &velocity_catalog,
                                            const ModeCatalog &catalog1,
                                            const ModeCatalog &catalog2,
                                            const std::vector<Real> &shell_amp1,
                                            const std::vector<Real> &shell_amp2,
                                            int npot_shells) {
  ClebschRealization best;
  best.trial_count = ClebschRealizationTrials(cfg);
  bool have_best = false;

  for (int trial = 0; trial < best.trial_count; ++trial) {
    RNG_State rstate1;
    RNG_State rstate2;
    rstate1.idum = -MixSeed(cfg.rseed, 1000 + 2*trial);
    rstate2.idum = -MixSeed(cfg.rseed, 1001 + 2*trial);

    ClebschRealization candidate;
    candidate.phi1_coeffs =
        GenerateShellWeightedScalarCoefficients(catalog1, shell_amp1, &rstate1);
    candidate.phi2_coeffs =
        GenerateShellWeightedScalarCoefficients(catalog2, shell_amp2, &rstate2);
    candidate.summary = EvaluateClebschRealizationSpectrum(
        velocity_catalog, catalog1, candidate.phi1_coeffs, catalog2,
        candidate.phi2_coeffs, cfg, npot_shells);
    candidate.trial_index = trial;
    candidate.trial_count = best.trial_count;

    if (!have_best || candidate.summary.loss < best.summary.loss) {
      best = std::move(candidate);
      have_best = true;
    }
  }

  return best;
}

int TensorIndex(int kout, int s, int t, int npot_shells) {
  const int ns = npot_shells + 1;
  return (kout * ns + s) * ns + t;
}

int Stream3DPotentialNHigh(const TurbulenceConfig &cfg) {
  int max_supported_shell = cfg.global_nx1 / 2;
  if (cfg.mesh_multi_d) {
    max_supported_shell = std::min(max_supported_shell, cfg.global_nx2 / 2);
  }
  if (cfg.mesh_three_d) {
    max_supported_shell = std::min(max_supported_shell, cfg.global_nx3 / 2);
  }
  return std::max(1, std::min(cfg.nhigh, max_supported_shell));
}

std::vector<Real> BuildClebschShellResponseTensor(const ModeCatalog &velocity_catalog,
                                                  const ModeCatalog &catalog1,
                                                  const ModeCatalog &catalog2,
                                                  int npot_shells, int nkout) {
  std::vector<Real> tensor((nkout + 1) * (npot_shells + 1) * (npot_shells + 1), 0.0);

  // For phi_j(x) = Re[sum_p c_j(p) exp(i p·x)], the Clebsch field satisfies
  // u_hat(r) = -sum_{p+q=r} (p x q) c_1(p) c_2(q). Averaging over independent
  // Gaussian coefficients gives the shell response
  // E_u(K) = sum_{s,t} T[K,s,t] P_s P_t, where P_s is the scalar-power weight
  // in shell s. In 3D stream mode we evaluate that response on the retained
  // output velocity catalog directly, including the output importance-sampling
  // boost, so the shell fit targets the exact spectral representation that is
  // later synthesized into the grid.
  const auto signed1 = BuildSignedModes(catalog1);
  const auto signed2 = BuildSignedModes(catalog2);
  const auto lookup2 = BuildSignedModeLookup(signed2, npot_shells);
  const int offset = npot_shells;
  const int span = 2 * offset + 1;

  for (int n = 0; n < velocity_catalog.KeptModes(); ++n) {
    const int rx = velocity_catalog.nkx[n];
    const int ry = velocity_catalog.nky[n];
    const int rz = velocity_catalog.nkz[n];
    const int out_shell = velocity_catalog.shell[n];
    const Real out_weight = 2.0 * velocity_catalog.boost[n] * velocity_catalog.boost[n];

    for (const auto &mode_p : signed1) {
      const int qx = rx - mode_p.nkx;
      const int qy = ry - mode_p.nky;
      const int qz = rz - mode_p.nkz;
      if (!SignedModeInBounds(qx, qy, qz, npot_shells)) continue;

      const auto it = lookup2.find(EncodeSignedModeKey(qx, qy, qz, offset, span));
      if (it == lookup2.end()) continue;
      const auto &mode_q = signed2[it->second];

      const Real cx = mode_p.ky*mode_q.kz - mode_p.kz*mode_q.ky;
      const Real cy = mode_p.kz*mode_q.kx - mode_p.kx*mode_q.kz;
      const Real cz = mode_p.kx*mode_q.ky - mode_p.ky*mode_q.kx;
      const Real geom = cx*cx + cy*cy + cz*cz;
      if (geom <= kTiny) continue;

      tensor[TensorIndex(out_shell, mode_p.shell, mode_q.shell, npot_shells)] +=
          out_weight * mode_p.weight * mode_q.weight * geom;
    }
  }

  return tensor;
}

void SymmetrizeClebschTensor(std::vector<Real> *tensor, int npot_shells, int nkout) {
  for (int kout = 1; kout <= nkout; ++kout) {
    for (int s = 1; s <= npot_shells; ++s) {
      for (int t = s + 1; t <= npot_shells; ++t) {
        const int idx_st = TensorIndex(kout, s, t, npot_shells);
        const int idx_ts = TensorIndex(kout, t, s, npot_shells);
        const Real avg = 0.5 * ((*tensor)[idx_st] + (*tensor)[idx_ts]);
        (*tensor)[idx_st] = avg;
        (*tensor)[idx_ts] = avg;
      }
    }
  }
}

void EnsureReachableTargetBand(const std::vector<Real> &tensor, const TurbulenceConfig &cfg,
                               int npot_shells) {
  for (int kout = cfg.nlow; kout <= cfg.nhigh; ++kout) {
    Real shell_sum = 0.0;
    for (int s = 1; s <= npot_shells; ++s) {
      for (int t = 1; t <= npot_shells; ++t) {
        shell_sum += tensor[TensorIndex(kout, s, t, npot_shells)];
      }
    }
    if (shell_sum <= kTiny) {
      FatalProblemSetup("3D stream-function shell response has zero support in target shell "
                        + std::to_string(kout) + ".");
    }
  }
}

ClebschSector ClassifyClebschSector(int kout, int s, int t) {
  if (s > kout && t > kout) return ClebschSector::HighHighCancel;
  const int half_shell = (kout + 1) / 2;
  if (s < half_shell && t >= half_shell) return ClebschSector::LowHigh;
  if (t < half_shell && s >= half_shell) return ClebschSector::HighLow;
  return ClebschSector::LocalLocal;
}

ClebschSectorResponse ComputeClebschSectorResponse(
    const std::vector<Real> &tensor, int npot_shells, int nkout,
    const std::vector<Real> *shell_power1 = nullptr,
    const std::vector<Real> *shell_power2 = nullptr) {
  ClebschSectorResponse response;
  for (int sector_idx = 0; sector_idx < kClebschSectorCount; ++sector_idx) {
    response.raw_shell_energy[sector_idx].assign(nkout + 1, 0.0);
    response.fitted_shell_energy[sector_idx].assign(nkout + 1, 0.0);
  }

  for (int kout = 1; kout <= nkout; ++kout) {
    for (int s = 1; s <= npot_shells; ++s) {
      for (int t = 1; t <= npot_shells; ++t) {
        const Real raw_value = tensor[TensorIndex(kout, s, t, npot_shells)];
        if (raw_value <= kTiny) continue;
        const int sector_idx = static_cast<int>(ClassifyClebschSector(kout, s, t));
        response.raw_shell_energy[sector_idx][kout] += raw_value;
        response.raw_total[sector_idx] += raw_value;

        if (shell_power1 != nullptr && shell_power2 != nullptr) {
          const Real fitted_value = raw_value * (*shell_power1)[s] * (*shell_power2)[t];
          response.fitted_shell_energy[sector_idx][kout] += fitted_value;
          response.fitted_total[sector_idx] += fitted_value;
        }
      }
    }
  }

  return response;
}

ShellFitEval EvaluateClebschShellFit(const std::vector<Real> &tensor,
                                     const TurbulenceConfig &cfg,
                                     int npot_shells,
                                     const std::vector<Real> &shell_power1,
                                     const std::vector<Real> &shell_power2) {
  const int nkout = cfg.nhigh;
  const int ns = npot_shells + 1;
  const Real leak_ref_low = TargetShellEnergy(std::max(cfg.nlow, 1), cfg.expo);
  const Real leak_ref_high = TargetShellEnergy(cfg.nhigh, cfg.expo);
  const Real smooth_weight = 1.0e-2;

  ShellFitEval eval;
  eval.pred.assign(nkout + 1, 0.0);
  eval.grad1.assign(npot_shells + 1, 0.0);
  eval.grad2.assign(npot_shells + 1, 0.0);

  std::vector<Real> kernel1((nkout + 1) * ns, 0.0);
  std::vector<Real> kernel2((nkout + 1) * ns, 0.0);
  for (int kout = 1; kout <= nkout; ++kout) {
    for (int s = 1; s <= npot_shells; ++s) {
      Real accum = 0.0;
      for (int t = 1; t <= npot_shells; ++t) {
        accum += tensor[TensorIndex(kout, s, t, npot_shells)] * shell_power2[t];
      }
      kernel1[kout * ns + s] = accum;
      eval.pred[kout] += shell_power1[s] * accum;
    }
    for (int t = 1; t <= npot_shells; ++t) {
      Real accum = 0.0;
      for (int s = 1; s <= npot_shells; ++s) {
        accum += tensor[TensorIndex(kout, s, t, npot_shells)] * shell_power1[s];
      }
      kernel2[kout * ns + t] = accum;
    }
  }

  for (int kout = 1; kout <= nkout; ++kout) {
    const Real pred = std::max(eval.pred[kout], kTiny);
    Real dloss_dpred = 0.0;

    if (kout >= cfg.nlow && kout <= cfg.nhigh) {
      const Real target = TargetShellEnergy(kout, cfg.expo);
      const Real resid = std::log(pred) - std::log(target);
      eval.loss += resid * resid;
      dloss_dpred = 2.0 * resid / pred;
    } else {
      const bool low_leak = (kout < cfg.nlow);
      const Real leak_weight = low_leak ? 4.0 : 6.0;
      const Real leak_ref = low_leak ? leak_ref_low : leak_ref_high;
      const Real scaled = pred / leak_ref;
      eval.loss += leak_weight * scaled * scaled;
      dloss_dpred = 2.0 * leak_weight * scaled / leak_ref;
    }

    for (int s = 1; s <= npot_shells; ++s) {
      eval.grad1[s] += dloss_dpred * kernel1[kout * ns + s];
      eval.grad2[s] += dloss_dpred * kernel2[kout * ns + s];
    }
  }

  for (int s = 2; s <= npot_shells; ++s) {
    const Real prev = std::max(shell_power1[s-1], kTiny);
    const Real curr = std::max(shell_power1[s], kTiny);
    const Real diff = std::log(curr) - std::log(prev);
    eval.loss += smooth_weight * diff * diff;
    eval.grad1[s] += 2.0 * smooth_weight * diff / curr;
    eval.grad1[s-1] -= 2.0 * smooth_weight * diff / prev;
  }
  for (int s = 2; s <= npot_shells; ++s) {
    const Real prev = std::max(shell_power2[s-1], kTiny);
    const Real curr = std::max(shell_power2[s], kTiny);
    const Real diff = std::log(curr) - std::log(prev);
    eval.loss += smooth_weight * diff * diff;
    eval.grad2[s] += 2.0 * smooth_weight * diff / curr;
    eval.grad2[s-1] -= 2.0 * smooth_weight * diff / prev;
  }

  return eval;
}

Real FitLogSlope(const std::vector<Real> &spectrum, int nlow, int nhigh) {
  Real sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  int npts = 0;
  for (int shell = nlow; shell <= nhigh; ++shell) {
    if (spectrum[shell] <= kTiny) continue;
    const Real x = std::log(static_cast<Real>(shell));
    const Real y = std::log(spectrum[shell]);
    sx += x;
    sy += y;
    sxx += x*x;
    sxy += x*y;
    ++npts;
  }
  if (npts < 2) return 0.0;
  const Real denom = npts * sxx - sx*sx;
  if (std::abs(denom) <= kTiny) return 0.0;
  return (npts * sxy - sx * sy) / denom;
}

Real LeakageFraction(const std::vector<Real> &spectrum, const TurbulenceConfig &cfg,
                     int nkout) {
  Real in_band = 0.0;
  Real out_band = 0.0;
  for (int shell = 1; shell <= nkout; ++shell) {
    if (shell >= cfg.nlow && shell <= cfg.nhigh) {
      in_band += spectrum[shell];
    } else {
      out_band += spectrum[shell];
    }
  }
  const Real total = in_band + out_band;
  return (total > kTiny) ? out_band / total : 0.0;
}

void InitializeClebschSeedPowers(ClebschSeedModel seed_model,
                                 const TurbulenceConfig &cfg, int npot_shells,
                                 std::vector<Real> *power1,
                                 std::vector<Real> *power2) {
  power1->assign(npot_shells + 1, 0.0);
  power2->assign(npot_shells + 1, 0.0);

  if (seed_model == ClebschSeedModel::Flat) {
    for (int shell = 1; shell <= npot_shells; ++shell) {
      (*power1)[shell] = 1.0;
      (*power2)[shell] = 1.0;
    }
    return;
  }

  const Real beta = (seed_model == ClebschSeedModel::LocalTriad)
                        ? 0.25 * (cfg.expo + 8.0)
                        : 0.5 * (cfg.expo + 4.0);
  const Real delta = (seed_model == ClebschSeedModel::LocalTriad) ? 0.5 : 0.0;
  for (int shell = 1; shell <= npot_shells; ++shell) {
    (*power1)[shell] = std::pow(static_cast<Real>(shell), -2.0 * (beta - delta));
    (*power2)[shell] = std::pow(static_cast<Real>(shell), -2.0 * (beta + delta));
  }
}

bool ClebschFitBeats(const ShellFitResult &candidate, const ShellFitResult &current_best,
                     const TurbulenceConfig &cfg) {
  constexpr Real eps = 1.0e-12;
  if (candidate.final_loss < current_best.final_loss - eps) return true;
  if (candidate.final_loss > current_best.final_loss + eps) return false;

  if (candidate.leakage_fraction < current_best.leakage_fraction - eps) return true;
  if (candidate.leakage_fraction > current_best.leakage_fraction + eps) return false;

  const Real candidate_slope_err = std::abs(candidate.fitted_slope + cfg.expo);
  const Real current_slope_err = std::abs(current_best.fitted_slope + cfg.expo);
  return candidate_slope_err < current_slope_err;
}

ShellFitResult FitClebschShellWeightsFromStart(
    const std::vector<Real> &tensor, const TurbulenceConfig &cfg, int npot_shells,
    const std::vector<Real> &start_power1, const std::vector<Real> &start_power2,
    ClebschSeedModel seed_model) {
  ShellFitResult result;
  result.shell_amp1.assign(npot_shells + 1, 0.0);
  result.shell_amp2.assign(npot_shells + 1, 0.0);
  result.shell_power1.assign(npot_shells + 1, 0.0);
  result.shell_power2.assign(npot_shells + 1, 0.0);
  result.chosen_seed_model = seed_model;

  std::vector<Real> shell_power1 = start_power1;
  std::vector<Real> shell_power2 = start_power2;
  const int ref_shell = std::min(std::max(cfg.nlow, 1), npot_shells);
  const Real ref = shell_power1[ref_shell];
  if (ref > kTiny) {
    for (int shell = 1; shell <= npot_shells; ++shell) {
      shell_power1[shell] /= ref;
      shell_power2[shell] /= ref;
    }
  }

  auto balance_reference_shell = [&](std::vector<Real> *power1, std::vector<Real> *power2) {
    const Real ref1 = std::max((*power1)[ref_shell], kTiny);
    const Real ref2 = std::max((*power2)[ref_shell], kTiny);
    const Real scale = std::sqrt(ref2 / ref1);
    for (int shell = 1; shell <= npot_shells; ++shell) {
      (*power1)[shell] *= scale;
      (*power2)[shell] /= scale;
    }
  };
  balance_reference_shell(&shell_power1, &shell_power2);

  ShellFitEval current = EvaluateClebschShellFit(tensor, cfg, npot_shells,
                                                 shell_power1, shell_power2);
  result.initial_loss = current.loss;
  std::vector<Real> best_power1 = shell_power1;
  std::vector<Real> best_power2 = shell_power2;
  ShellFitEval best_eval = current;
  auto rescale_band_mean = [&](std::vector<Real> *power1, std::vector<Real> *power2) {
    ShellFitEval scaled = EvaluateClebschShellFit(tensor, cfg, npot_shells, *power1, *power2);
    Real pred_mean = 0.0;
    Real target_mean = 0.0;
    for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
      pred_mean += scaled.pred[shell];
      target_mean += TargetShellEnergy(shell, cfg.expo);
    }
    if (pred_mean <= kTiny) return;
    const Real scale_sqrt = std::sqrt(target_mean / pred_mean);
    for (int shell = 1; shell <= npot_shells; ++shell) {
      (*power1)[shell] *= scale_sqrt;
      (*power2)[shell] *= scale_sqrt;
    }
  };
  auto multiplicative_sweep = [&](std::vector<Real> *power1, std::vector<Real> *power2) {
    const int nkout = cfg.nhigh;
    const int ns = npot_shells + 1;
    std::vector<Real> pred(nkout + 1, 0.0);
    std::vector<Real> kernel1((nkout + 1) * ns, 0.0);
    std::vector<Real> kernel2((nkout + 1) * ns, 0.0);

    for (int kout = cfg.nlow; kout <= cfg.nhigh; ++kout) {
      for (int s = 1; s <= npot_shells; ++s) {
        Real accum = 0.0;
        for (int t = 1; t <= npot_shells; ++t) {
          accum += tensor[TensorIndex(kout, s, t, npot_shells)] * (*power2)[t];
        }
        kernel1[kout * ns + s] = accum;
        pred[kout] += (*power1)[s] * accum;
      }
      for (int t = 1; t <= npot_shells; ++t) {
        Real accum = 0.0;
        for (int s = 1; s <= npot_shells; ++s) {
          accum += tensor[TensorIndex(kout, s, t, npot_shells)] * (*power1)[s];
        }
        kernel2[kout * ns + t] = accum;
      }
    }

    std::vector<Real> next1 = *power1;
    std::vector<Real> next2 = *power2;
    for (int shell = 1; shell <= npot_shells; ++shell) {
      Real num1 = 0.0;
      Real den1 = 0.0;
      Real num2 = 0.0;
      Real den2 = 0.0;
      for (int kout = cfg.nlow; kout <= cfg.nhigh; ++kout) {
        const Real target = TargetShellEnergy(kout, cfg.expo);
        const Real inv_pred = 1.0 / std::max(pred[kout], kTiny);
        num1 += kernel1[kout * ns + shell] * target * inv_pred;
        den1 += kernel1[kout * ns + shell];
        num2 += kernel2[kout * ns + shell] * target * inv_pred;
        den2 += kernel2[kout * ns + shell];
      }
      if (den1 > kTiny) {
        next1[shell] *= num1 / den1;
      }
      if (den2 > kTiny) {
        next2[shell] *= num2 / den2;
      }
      next1[shell] = std::min(1.0e30, std::max(1.0e-30, next1[shell]));
      next2[shell] = std::min(1.0e30, std::max(1.0e-30, next2[shell]));
    }

    for (int shell = 2; shell < npot_shells; ++shell) {
      const Real smooth1 = std::exp((std::log(next1[shell - 1]) + std::log(next1[shell]) +
                                     std::log(next1[shell + 1])) / 3.0);
      const Real smooth2 = std::exp((std::log(next2[shell - 1]) + std::log(next2[shell]) +
                                     std::log(next2[shell + 1])) / 3.0);
      next1[shell] = std::pow(next1[shell], 0.8) * std::pow(smooth1, 0.2);
      next2[shell] = std::pow(next2[shell], 0.8) * std::pow(smooth2, 0.2);
    }

    *power1 = std::move(next1);
    *power2 = std::move(next2);
    balance_reference_shell(power1, power2);
    rescale_band_mean(power1, power2);
  };

  int stalled_sweeps = 0;
  for (int iter = 0; iter < 400; ++iter) {
    multiplicative_sweep(&shell_power1, &shell_power2);
    current = EvaluateClebschShellFit(tensor, cfg, npot_shells, shell_power1, shell_power2);
    if (current.loss < best_eval.loss) {
      best_power1 = shell_power1;
      best_power2 = shell_power2;
      best_eval = current;
      stalled_sweeps = 0;
    } else {
      ++stalled_sweeps;
      if (stalled_sweeps >= 50) break;
    }
  }

  shell_power1 = best_power1;
  shell_power2 = best_power2;
  current = best_eval;
  Real step_scale = 0.1;

  for (int iter = 0; iter < 250; ++iter) {
    std::vector<Real> grad_log1(npot_shells + 1, 0.0);
    std::vector<Real> grad_log2(npot_shells + 1, 0.0);
    Real max_grad = 0.0;
    for (int shell = 1; shell <= npot_shells; ++shell) {
      grad_log1[shell] = shell_power1[shell] * current.grad1[shell];
      grad_log2[shell] = shell_power2[shell] * current.grad2[shell];
      max_grad = std::max(max_grad, std::abs(grad_log1[shell]));
      max_grad = std::max(max_grad, std::abs(grad_log2[shell]));
    }
    if (max_grad < 1.0e-10) break;

    Real step = step_scale / max_grad;
    bool accepted = false;
    for (int trial = 0; trial < 10; ++trial) {
      std::vector<Real> trial_power1(npot_shells + 1, 0.0);
      std::vector<Real> trial_power2(npot_shells + 1, 0.0);
      for (int shell = 1; shell <= npot_shells; ++shell) {
        Real log_power1 = std::log(std::max(shell_power1[shell], kTiny));
        log_power1 -= step * grad_log1[shell];
        log_power1 = std::min(30.0, std::max(-30.0, log_power1));
        trial_power1[shell] = std::exp(log_power1);

        Real log_power2 = std::log(std::max(shell_power2[shell], kTiny));
        log_power2 -= step * grad_log2[shell];
        log_power2 = std::min(30.0, std::max(-30.0, log_power2));
        trial_power2[shell] = std::exp(log_power2);
      }
      balance_reference_shell(&trial_power1, &trial_power2);

      ShellFitEval trial_eval = EvaluateClebschShellFit(tensor, cfg, npot_shells,
                                                        trial_power1, trial_power2);
      if (trial_eval.loss < current.loss) {
        shell_power1 = std::move(trial_power1);
        shell_power2 = std::move(trial_power2);
        current = std::move(trial_eval);
        if (current.loss < best_eval.loss) {
          best_power1 = shell_power1;
          best_power2 = shell_power2;
          best_eval = current;
        }
        step_scale = std::min(step_scale * 1.2, 2.0);
        accepted = true;
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      step_scale *= 0.5;
      if (step_scale < 1.0e-6) break;
    }
  }

  shell_power1 = best_power1;
  shell_power2 = best_power2;
  current = best_eval;

  Real pred_mean = 0.0;
  Real target_mean = 0.0;
  int nband = 0;
  for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
    pred_mean += current.pred[shell];
    target_mean += TargetShellEnergy(shell, cfg.expo);
    ++nband;
  }
  if (pred_mean > kTiny && nband > 0) {
    const Real scale = target_mean / pred_mean;
    const Real scale_sqrt = std::sqrt(scale);
    for (int shell = 1; shell <= npot_shells; ++shell) {
      shell_power1[shell] *= scale_sqrt;
      shell_power2[shell] *= scale_sqrt;
    }
    current = EvaluateClebschShellFit(tensor, cfg, npot_shells,
                                      shell_power1, shell_power2);
  }

  for (int shell = 1; shell <= npot_shells; ++shell) {
    result.shell_amp1[shell] = std::sqrt(std::max(shell_power1[shell], 0.0));
    result.shell_amp2[shell] = std::sqrt(std::max(shell_power2[shell], 0.0));
  }
  result.shell_power1 = shell_power1;
  result.shell_power2 = shell_power2;
  result.predicted_shell_energy = current.pred;
  result.final_loss = current.loss;
  result.leakage_fraction = LeakageFraction(current.pred, cfg, cfg.nhigh);
  result.fitted_slope = FitLogSlope(current.pred, cfg.nlow, cfg.nhigh);
  result.converged = (result.final_loss <= result.initial_loss);

  return result;
}

ShellFitResult FitClebschShellWeights(const std::vector<Real> &tensor,
                                      const TurbulenceConfig &cfg,
                                      int npot_shells) {
  const std::array<ClebschSeedModel, 3> seed_models = {
      ClebschSeedModel::LocalTriad,
      ClebschSeedModel::CamilloSymmetric,
      ClebschSeedModel::Flat};
  std::vector<ShellFitStartSummary> fit_starts;
  fit_starts.reserve(seed_models.size());

  ShellFitResult best_result;
  bool have_best = false;
  std::size_t best_start_index = 0;
  for (std::size_t seed_idx = 0; seed_idx < seed_models.size(); ++seed_idx) {
    std::vector<Real> start_power1;
    std::vector<Real> start_power2;
    InitializeClebschSeedPowers(seed_models[seed_idx], cfg, npot_shells,
                                &start_power1, &start_power2);
    ShellFitResult candidate = FitClebschShellWeightsFromStart(
        tensor, cfg, npot_shells, start_power1, start_power2, seed_models[seed_idx]);
    fit_starts.push_back({seed_models[seed_idx], candidate.initial_loss, candidate.final_loss,
                          candidate.fitted_slope, candidate.leakage_fraction,
                          candidate.converged, false});
    if (!have_best || ClebschFitBeats(candidate, best_result, cfg)) {
      best_result = std::move(candidate);
      best_start_index = seed_idx;
      have_best = true;
    }
  }

  fit_starts[best_start_index].chosen = true;
  best_result.fit_starts = std::move(fit_starts);
  return best_result;
}

ClebschQualityGate EvaluateClebschQualityGate(const ShellFitResult &fit,
                                              const RealizedSpectrumSummary &realization,
                                              const TurbulenceConfig &cfg) {
  ClebschQualityGate gate;
  gate.allow_poor_fit = cfg.stream_allow_poor_fit;
  gate.deterministic_converged = fit.converged;
  gate.deterministic_slope_ok =
      (std::abs(fit.fitted_slope + cfg.expo) <= gate.deterministic_slope_tolerance);
  gate.deterministic_leakage_ok = (fit.leakage_fraction <= gate.leakage_tolerance);
  gate.deterministic_pass = gate.deterministic_converged &&
                            gate.deterministic_slope_ok &&
                            gate.deterministic_leakage_ok;
  gate.realized_slope_ok =
      (std::abs(realization.fitted_slope + cfg.expo) <= gate.realized_slope_tolerance);
  gate.realized_leakage_ok = (realization.leakage_fraction <= gate.leakage_tolerance);
  gate.realized_pass = gate.realized_slope_ok && gate.realized_leakage_ok;
  gate.pass = gate.deterministic_pass && gate.realized_pass;
  gate.forced_continue = gate.allow_poor_fit && !gate.pass;
  return gate;
}

void LogModeSummary(const TurbulenceConfig &cfg, const ModeCatalog &catalog,
                    const std::string &label) {
  std::cout << "  " << label << ": total_modes_in_shell=" << catalog.total_modes
            << ", kept_via_importance_sampling=" << catalog.KeptModes() << std::endl;
  if (cfg.method == VelocityMethod::StreamFunction && cfg.stream_2d) {
    std::cout << "  " << label << ": stream-function coefficients preserve "
              << "E_u(k) ~ k^(-expo) through |psi_k| ~ k^{-(expo+3)/2}" << std::endl;
  }
}

void InitTurbulentVelocity(MeshBlockPack *pmbp, ParameterInput *pin, Real den,
                           bool projection_true_2d) {
  const WallTimePoint init_start_time = WallClock::now();
  Real recorded_init_wall_seconds = -1.0;
  auto log_init_time = [&init_start_time, &recorded_init_wall_seconds]() {
    const Real wall_seconds = (recorded_init_wall_seconds > 0.0)
                                  ? recorded_init_wall_seconds
                                  : GlobalMaxReal(ElapsedWallSeconds(init_start_time));
    if (global_variable::my_rank == 0) {
      std::cout << "  velocity_init_wall_seconds=" << wall_seconds << std::endl;
    }
  };

  const TurbulenceConfig cfg = ReadTurbulenceConfig(pin, pmbp->pmesh, projection_true_2d);
  auto &hydro = *pmbp->phydro;
  hydro.use_scalar_face_velocity = false;
  hydro.scalar_vface.reset();
  const bool face_requested = cfg.divfree_scalar_flux;
  const bool face_scalar_ok = hydro.scalar_only && (hydro.nscalars > 0);

  if (face_requested && !hydro.scalar_only && global_variable::my_rank == 0) {
    std::cout << "  NOTE: divfree_scalar_flux ignored (scalar_only=false)." << std::endl;
  } else if (face_requested && hydro.nscalars == 0 && global_variable::my_rank == 0) {
    std::cout << "  NOTE: divfree_scalar_flux ignored (nscalars=0)." << std::endl;
  }

  if (cfg.method == VelocityMethod::Projection) {
    RNG_State rstate;
    rstate.idum = -cfg.rseed;
    const ModeCatalog catalog = BuildModeCatalog(cfg, &rstate);

    if (catalog.KeptModes() == 0) {
      std::cout << "### WARNING: No turbulent modes kept in range [" << cfg.nlow << ", "
                << cfg.nhigh << "]. Using zero velocity." << std::endl;
      log_init_time();
      return;
    }

    if (global_variable::my_rank == 0) {
      std::cout << "Initializing turbulent velocity field:" << std::endl
                << "  velocity_method=projection" << std::endl
                << "  v_rms=" << cfg.v_rms << ", nlow=" << cfg.nlow
                << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo
                << ", sol_frac=" << cfg.sol_frac << std::endl
                << "  true_2d=" << (cfg.projection_true_2d ? "true" : "false")
                << ", active_velocity_components=" << (cfg.active_v3 ? 3 : 2)
                << std::endl
                << "  k_crit=" << cfg.k_crit
                << " (wavenumber=" << cfg.k_crit_mag << ")" << std::endl;
      LogModeSummary(cfg, catalog, "velocity");
    }

    const VectorCoefficients coeffs = GenerateProjectionCoefficients(cfg, catalog, &rstate);
    SynthesizeVelocityFromVectorModes(pmbp, cfg, den, catalog, coeffs);
    const VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
    bool use_face_velocity = face_requested && face_scalar_ok;
    if (use_face_velocity && cfg.sol_frac < 1.0) {
      use_face_velocity = false;
      if (global_variable::my_rank == 0) {
        std::cout << "  NOTE: divfree_scalar_flux ignored when turb_sol_frac < 1 in "
                  << "projection mode." << std::endl;
      }
    }
    if (use_face_velocity) {
      if (global_variable::my_rank == 0) {
        std::cout << "  using divergence-free face velocity for scalar fluxes" << std::endl;
      }
      const VectorCoefficients apot_coeffs =
          ConvertVelocityToVectorPotentialCoefficients(catalog, coeffs);
      PopulateScalarFaceVelocitiesFromVectorPotential(pmbp, cfg, catalog, apot_coeffs, stats);
    }
    const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
    recorded_init_wall_seconds = init_wall_seconds;
    WriteProjectionDiagnostics(cfg, catalog, coeffs, stats, init_wall_seconds);
    log_init_time();
    return;
  }

  if (global_variable::my_rank == 0) {
    std::cout << "Initializing turbulent velocity field:" << std::endl
              << "  velocity_method=stream_function" << std::endl
              << "  v_rms=" << cfg.v_rms << ", nlow=" << cfg.nlow
              << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo
              << ", sol_frac=" << cfg.sol_frac << std::endl
              << "  active_velocity_components=" << (cfg.active_v3 ? 3 : 2) << std::endl
              << "  k_crit=" << cfg.k_crit
              << " (wavenumber=" << cfg.k_crit_mag << ")" << std::endl
              << "  turb_sol_frac is ignored in stream-function mode." << std::endl;
  }

  if (cfg.stream_2d) {
    RNG_State rstate;
    rstate.idum = -cfg.rseed;
    const ModeCatalog catalog = BuildModeCatalog(cfg, &rstate);

    if (catalog.KeptModes() == 0) {
      std::cout << "### WARNING: No turbulent modes kept in range [" << cfg.nlow << ", "
                << cfg.nhigh << "]. Using zero velocity." << std::endl;
      log_init_time();
      return;
    }

    if (global_variable::my_rank == 0) {
      std::cout << "  stream_geometry=2D (nx3=1 forces v3=0)" << std::endl;
      LogModeSummary(cfg, catalog, "psi");
    }

    const ScalarCoefficients psi_coeffs = GenerateStream2DPsiCoefficients(cfg, catalog,
                                                                          &rstate);
    const VectorCoefficients vel_coeffs = ConvertPsiToVelocityCoefficients(catalog, psi_coeffs);
    SynthesizeVelocityFromVectorModes(pmbp, cfg, den, catalog, vel_coeffs);
    const VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
    if (face_requested && face_scalar_ok) {
      if (global_variable::my_rank == 0) {
        std::cout << "  using divergence-free face velocity for scalar fluxes" << std::endl;
      }
      const VectorCoefficients apot_coeffs =
          ConvertVelocityToVectorPotentialCoefficients(catalog, vel_coeffs);
      PopulateScalarFaceVelocitiesFromVectorPotential(pmbp, cfg, catalog, apot_coeffs, stats);
    }
    const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
    recorded_init_wall_seconds = init_wall_seconds;
    WriteStream2DDiagnostics(cfg, catalog, psi_coeffs, vel_coeffs, stats, init_wall_seconds);
    log_init_time();
    return;
  }

  RNG_State rstate1;
  RNG_State rstate2;
  RNG_State rstatev;
  rstate1.idum = -MixSeed(cfg.rseed, 1);
  rstate2.idum = -MixSeed(cfg.rseed, 2);
  rstatev.idum = -MixSeed(cfg.rseed, 3);

  const int phi_nhigh = Stream3DPotentialNHigh(cfg);
  if (phi_nhigh < cfg.nhigh) {
    FatalProblemSetup("3D stream-function mode requires potential support through "
                      "<problem>/turb_nhigh on the current mesh. Increase the grid size "
                      "or lower <problem>/turb_nhigh.");
  }

  const ModeCatalog catalog1 = BuildModeCatalog(cfg, 1, phi_nhigh, &rstate1);
  const ModeCatalog catalog2 = BuildModeCatalog(cfg, 1, phi_nhigh, &rstate2);
  const ModeCatalog velocity_catalog = BuildModeCatalog(cfg, &rstatev);
  if (catalog1.KeptModes() == 0 || catalog2.KeptModes() == 0) {
    std::cout << "### WARNING: No turbulent modes kept in range [" << cfg.nlow << ", "
              << cfg.nhigh << "] for 3D stream-function initialization. Using zero velocity."
              << std::endl;
    log_init_time();
    return;
  }
  if (velocity_catalog.KeptModes() == 0) {
    std::cout << "### WARNING: No output velocity modes kept in range [" << cfg.nlow << ", "
              << cfg.nhigh << "] for 3D stream-function initialization. Using zero velocity."
              << std::endl;
    log_init_time();
    return;
  }

  const std::vector<int> shell_counts1 = CountModesPerShell(catalog1, phi_nhigh);
  const std::vector<int> shell_counts2 = CountModesPerShell(catalog2, phi_nhigh);
  const std::vector<int> velocity_shell_counts = CountModesPerShell(velocity_catalog, cfg.nhigh);
  EnsureShellCoverage(shell_counts1, 1, phi_nhigh, "phi_1");
  EnsureShellCoverage(shell_counts2, 1, phi_nhigh, "phi_2");
  EnsureShellCoverage(velocity_shell_counts, cfg.nlow, cfg.nhigh, "velocity");

  std::vector<Real> tensor = BuildClebschShellResponseTensor(
      velocity_catalog, catalog1, catalog2, phi_nhigh, cfg.nhigh);
  SymmetrizeClebschTensor(&tensor, phi_nhigh, cfg.nhigh);
  EnsureReachableTargetBand(tensor, cfg, phi_nhigh);
  const ShellFitResult fit = FitClebschShellWeights(tensor, cfg, phi_nhigh);
  const ClebschSectorResponse sector_response =
      ComputeClebschSectorResponse(tensor, phi_nhigh, cfg.nhigh,
                                   &fit.shell_power1, &fit.shell_power2);

  if (global_variable::my_rank == 0) {
    std::cout << "  stream_geometry=3D (Clebsch form v = grad(phi1) x grad(phi2))"
              << std::endl;
    std::cout << "  phi_shell_max=" << phi_nhigh
              << " (target velocity shells extend to " << cfg.nhigh << ")" << std::endl;
    LogModeSummary(cfg, catalog1, "phi_1");
    LogModeSummary(cfg, catalog2, "phi_2");
    LogModeSummary(cfg, velocity_catalog, "velocity");
    std::cout << "  3D stream shell fit:"
              << " initial_loss=" << fit.initial_loss
              << ", final_loss=" << fit.final_loss
              << ", fitted_in_band_slope=" << fit.fitted_slope
              << ", estimated_out_of_band_leakage=" << fit.leakage_fraction
              << ", seed_model=" << ClebschSeedModelName(fit.chosen_seed_model)
              << std::endl;
    for (const auto &start : fit.fit_starts) {
      std::cout << "    fit_start[" << ClebschSeedModelName(start.seed_model) << "]:"
                << " initial_loss=" << start.initial_loss
                << ", final_loss=" << start.final_loss
                << ", fitted_slope=" << start.fitted_slope
                << ", estimated_out_of_band_leakage=" << start.leakage_fraction
                << (start.chosen ? " <- chosen" : "") << std::endl;
    }
    if (!fit.converged || std::abs(fit.fitted_slope + cfg.expo) > 0.5 ||
        fit.leakage_fraction > 0.05) {
      std::cout << "  WARNING: 3D stream-function shell fit only matches the requested "
                << "velocity spectrum approximately for this band/seed." << std::endl;
    }
  }

  const ClebschRealization realization =
      SelectClebschRealization(cfg, velocity_catalog, catalog1, catalog2,
                               fit.shell_amp1, fit.shell_amp2, phi_nhigh);
  if (global_variable::my_rank == 0) {
    std::cout << "  3D stream realization selection:"
              << " trials=" << realization.trial_count
              << ", chosen_trial=" << realization.trial_index
              << ", realized_in_band_slope=" << realization.summary.fitted_slope
              << ", realized_out_of_band_leakage="
              << realization.summary.leakage_fraction
              << ", realized_loss=" << realization.summary.loss << std::endl;
  }

  const ClebschQualityGate quality_gate =
      EvaluateClebschQualityGate(fit, realization.summary, cfg);
  if (global_variable::my_rank == 0) {
    std::cout << "  3D stream quality gate:"
              << " deterministic_pass="
              << (quality_gate.deterministic_pass ? "true" : "false")
              << ", realized_pass=" << (quality_gate.realized_pass ? "true" : "false")
              << ", overall_pass=" << (quality_gate.pass ? "true" : "false")
              << std::endl;
  }
  if (!quality_gate.pass) {
    const std::string gate_msg =
        "3D stream-function fit failed quality gate: deterministic_converged="
        + std::string(quality_gate.deterministic_converged ? "true" : "false")
        + ", deterministic_slope=" + std::to_string(fit.fitted_slope)
        + ", deterministic_leakage=" + std::to_string(fit.leakage_fraction)
        + ", realized_slope=" + std::to_string(realization.summary.fitted_slope)
        + ", realized_leakage=" + std::to_string(realization.summary.leakage_fraction)
        + ". Set <problem>/turb_stream_allow_poor_fit=true to continue anyway.";
    if (!cfg.stream_allow_poor_fit) {
      FatalProblemSetup(gate_msg);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "  WARNING: forcing 3D stream-function initialization despite failed "
                   "quality gate because turb_stream_allow_poor_fit=true."
                << std::endl;
      std::cout << "  WARNING: " << gate_msg << std::endl;
    }
  }

  const ClebschVelocityModes velocity_modes =
      BuildClebschVelocityModes(velocity_catalog, catalog1, realization.phi1_coeffs,
                                catalog2, realization.phi2_coeffs, phi_nhigh, cfg.nhigh);
  SynthesizeVelocityFromVectorModes(pmbp, cfg, den, velocity_catalog, velocity_modes.coeffs);
  const VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
  if (face_requested && face_scalar_ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "  using divergence-free face velocity for scalar fluxes" << std::endl;
    }
    const VectorCoefficients apot_coeffs =
        ConvertVelocityToVectorPotentialCoefficients(velocity_catalog, velocity_modes.coeffs);
    PopulateScalarFaceVelocitiesFromVectorPotential(pmbp, cfg, velocity_catalog,
                                                    apot_coeffs, stats);
  }
  const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
  recorded_init_wall_seconds = init_wall_seconds;
  WriteStream3DDiagnostics(cfg, catalog1, realization.phi1_coeffs, catalog2,
                           realization.phi2_coeffs, velocity_catalog, fit,
                           realization, velocity_modes, sector_response,
                           quality_gate, phi_nhigh, stats, init_wall_seconds);
  log_init_time();
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar mixing
//
//  Initializes density, turbulent velocity, and passive scalars.
//  Scalars can be initialized from a shared x1 step profile or per-scalar uniform offsets.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  const int nscalars = pmbp->phydro->nscalars;
  const int nhydro = pmbp->phydro->nhydro;
  const bool true_2d = UseTrue2DVelocity(pin, pmy_mesh_);

  if (nscalars == 0) {
    FatalProblemSetup("Scalar mixing requires nscalars > 0.");
  }

  const Real legacy_G = pin->GetOrAddReal("problem", "mean_gradient", 1.0);
  const Real scalar_init_default = pin->GetOrAddReal("problem", "scalar_init", 0.5);
  const bool scalar_use_x1_step =
      pin->GetOrAddBoolean("problem", "scalar_use_x1_step", false);
  const Real scalar_step_x1 = pin->GetOrAddReal("problem", "scalar_step_x1", 0.0);
  const Real scalar_left = pin->GetOrAddReal("problem", "scalar_left", 0.0);
  const Real scalar_right = pin->GetOrAddReal("problem", "scalar_right", 1.0);
  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  const Real lx = x1max - x1min;

  scalar_forcing_cfg.resize(nscalars);
  std::vector<Real> scalar_init(nscalars, scalar_init_default);
  const Real scalar_noise_default = pin->GetOrAddReal("problem", "scalar_noise", 0.0);
  std::vector<Real> scalar_noise(nscalars, scalar_noise_default);
  bool clipped_inits = false;
  bool clipped_noise = false;
  bool noise_enabled = (scalar_noise_default > 0.0);

  for (int ns = 0; ns < nscalars; ++ns) {
    ScalarForcingConfig cfg;
    const std::string prefix = "scalar" + std::to_string(ns) + "_";
    const int default_mode = (ns == 0) ? kScalarMeanGradient : kScalarNone;
    const Real default_G = (ns == 0) ? legacy_G : 1.0;

    cfg.mode = pin->GetOrAddInteger("problem", prefix + "mode", default_mode);
    cfg.G = pin->GetOrAddReal("problem", prefix + "mean_gradient", default_G);
    cfg.tau = pin->GetOrAddReal("problem", prefix + "tau", 1.0);
    cfg.theta1 = pin->GetOrAddReal("problem", prefix + "theta1", 1.0e-6);
    cfg.theta3 = pin->GetOrAddReal("problem", prefix + "theta3", 1.0);
    cfg.a = pin->GetOrAddReal("problem", prefix + "a", 0.5);
    if (pin->DoesParameterExist("problem", prefix + "theta2")) {
      cfg.theta2 = pin->GetReal("problem", prefix + "theta2");
    } else {
      cfg.theta2 = cfg.theta1 + cfg.a * (cfg.theta3 - cfg.theta1);
    }
    cfg.chiL0 = pin->GetOrAddReal("problem", prefix + "chiL0", 1.0);
    cfg.chiR0 = pin->GetOrAddReal("problem", prefix + "chiR0", 1.0);
    cfg.thetaL = pin->GetOrAddReal("problem", prefix + "thetaL", cfg.theta1);
    cfg.thetaR = pin->GetOrAddReal("problem", prefix + "thetaR", cfg.theta3);
    cfg.xL0 = pin->GetOrAddReal("problem", prefix + "xL0", x1min + 0.20 * lx);
    cfg.xL1 = pin->GetOrAddReal("problem", prefix + "xL1", x1min + 0.30 * lx);
    cfg.xR0 = pin->GetOrAddReal("problem", prefix + "xR0", x1min + 0.70 * lx);
    cfg.xR1 = pin->GetOrAddReal("problem", prefix + "xR1", x1min + 0.80 * lx);
    cfg.mask_w = pin->GetOrAddReal("problem", prefix + "mask_w", 0.02 * lx);
    cfg.reaction_type = pin->GetOrAddInteger("problem", prefix + "reaction_type",
                                             kReactionBistable);
    cfg.theta_floor = pin->GetOrAddReal("problem", prefix + "theta_floor", cfg.theta1);
    cfg.clip_reaction_interval =
        pin->GetOrAddBoolean("problem", prefix + "clip_reaction_interval", true);

    scalar_forcing_cfg[ns] = cfg;

    Real init = pin->GetOrAddReal("problem", prefix + "init", scalar_init_default);
    if (init < cfg.theta_floor) {
      init = cfg.theta_floor;
      clipped_inits = true;
    }
    scalar_init[ns] = init;

    Real noise = pin->GetOrAddReal("problem", prefix + "noise", scalar_noise_default);
    if (noise < 0.0) {
      noise = 0.0;
      clipped_noise = true;
    }
    scalar_noise[ns] = noise;
    if (noise > 0.0) {
      noise_enabled = true;
    }
  }

  user_srcs_func = ScalarForcingSource;
  user_hist_func = ScalarMixingHistory;

  if (restart) {
    return;
  }

  // Capture variables for kernel
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  const int nx1 = indcs.nx1;
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  const Real gm1 = eos.gamma - 1.0;

  // Read initial conditions
  const Real den = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real temp = pin->GetOrAddReal("problem", "temp0", 1.0);
  const int nmb = pmbp->nmb_thispack;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar mixing problem:" << std::endl
              << "  scalar_ic = "
              << (scalar_use_x1_step ? "x1_step" : "uniform") << std::endl
              << "  true_2d = " << (true_2d ? "true" : "false") << std::endl
              << "  nscalars = " << nscalars << std::endl
              << "  scalar_init (default) = " << scalar_init_default << std::endl
              << "  scalar0_mode = " << scalar_forcing_cfg[0].mode << std::endl
              << "  scalar0_mean_gradient = " << scalar_forcing_cfg[0].G << std::endl;
    if (scalar_use_x1_step) {
      std::cout << "  scalar_left = " << scalar_left << std::endl
                << "  scalar_right = " << scalar_right << std::endl
                << "  scalar_step_x1 = " << scalar_step_x1 << std::endl;
    }
    if (clipped_inits) {
      std::cout << "  NOTE: scalar_init clipped to theta_floor for positivity." << std::endl;
    }
    if (noise_enabled) {
      std::cout << "  scalar_noise default = " << scalar_noise_default << std::endl;
    }
    if (clipped_noise) {
      std::cout << "  NOTE: negative scalar_noise values clipped to 0." << std::endl;
    }
  }

  // Initialize density and energy
  par_for("scalar_mix_pgen_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IDN,k,j,i) = den;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = den * temp / gm1;
  });

  // Initialize turbulent velocity field
  InitTurbulentVelocity(pmbp, pin, den, true_2d);

  // Update energy with kinetic energy
  par_for("scalar_mix_pgen_energy", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real ke_cell = 0.5 * (SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
                                SQR(u0(m,IM3,k,j,i))) / u0(m,IDN,k,j,i);
    u0(m,IEN,k,j,i) += ke_cell;
  });

  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
  for (int ns = 0; ns < nscalars; ++ns) {
    const Real init = scalar_init[ns];
    const Real noise = scalar_noise[ns];
    const Real theta_floor = scalar_forcing_cfg[ns].theta_floor;
    const std::string label = "scalar_mix_pgen_scalar_" + std::to_string(ns);
    par_for(label, DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real rho = u0(m, IDN, k, j, i);
      Real theta = init;
      if (scalar_use_x1_step) {
        Real &x1min_ = size.d_view(m).x1min;
        Real &x1max_ = size.d_view(m).x1max;
        const Real x1v = CellCenterX(i-is, nx1, x1min_, x1max_);
        theta = (x1v < scalar_step_x1) ? scalar_left : scalar_right;
      }
      if (noise > 0.0) {
        auto rand_gen = rand_pool64.get_state();
        const Real r = 2.0*static_cast<Real>(rand_gen.frand()) - 1.0;
        rand_pool64.free_state(rand_gen);
        theta += noise * r;
      }
      if (theta < theta_floor) {
        theta = theta_floor;
      }
      u0(m, nhydro + ns, k, j, i) = rho * theta;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn ScalarForcingSource
//  \brief Apply per-scalar source terms (mean gradient, reaction, and/or sponge).

void ScalarForcingSource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto &size = pmbp->pmb->mb_size;
  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;

  if (nscalars == 0 || scalar_forcing_cfg.size() < static_cast<size_t>(nscalars)) {
    return;
  }

  const int nmkji = pmbp->nmb_thispack * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;

  for (int ns = 0; ns < nscalars; ++ns) {
    const ScalarForcingConfig cfg = scalar_forcing_cfg[ns];
    if (cfg.mode == kScalarNone) {
      continue;
    }

    if (cfg.mode == kScalarMeanGradient) {
      const std::string label = "scalar_mean_grad_" + std::to_string(ns);
      par_for(label, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real density = w0(m, IDN, k, j, i);
        const Real vx = w0(m, IVX, k, j, i);
        u0(m, nhydro + ns, k, j, i) += density * cfg.G * vx * bdt;
      });
      continue;
    }

    if (cfg.mode == kScalarMassConservingAC) {
      Real num_local = 0.0;
      Real den_local = 0.0;
      const std::string reduce_label = "scalar_mc_ac_reduce_" + std::to_string(ns);
      Kokkos::parallel_reduce(reduce_label, Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, Real &s_num, Real &s_den) {
        const int m = idx / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        int i = (idx - m * nkji - k * nji - j * nx1) + is;
        k += ks;
        j += js;

        const Real rho = w0(m, IDN, k, j, i);
        const Real theta = w0(m, nhydro + ns, k, j, i);
        const Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
        const Real reaction = BistableReaction(theta, cfg.tau, cfg.theta1, cfg.theta2,
                                               cfg.theta3);
        s_num += rho * reaction * vol;
        s_den += rho * vol;
      }, Kokkos::Sum<Real>(num_local), Kokkos::Sum<Real>(den_local));

#if MPI_PARALLEL_ENABLED
      Real local_sum[2] = {num_local, den_local};
      Real global_sum[2] = {0.0, 0.0};
      MPI_Allreduce(local_sum, global_sum, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      num_local = global_sum[0];
      den_local = global_sum[1];
#endif
      const Real rbar = (den_local > 0.0) ? num_local / den_local : 0.0;
      const std::string update_label = "scalar_mc_ac_update_" + std::to_string(ns);
      par_for(update_label, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real rho = w0(m, IDN, k, j, i);
        const Real theta = w0(m, nhydro + ns, k, j, i);
        const Real reaction = BistableReaction(theta, cfg.tau, cfg.theta1, cfg.theta2,
                                               cfg.theta3);
        u0(m, nhydro + ns, k, j, i) += rho * (reaction - rbar) * bdt;
      });
      continue;
    }

    if (cfg.mode == kScalarSponge || cfg.mode == kScalarSpongeReaction) {
      const std::string label = (cfg.mode == kScalarSponge)
                                    ? "scalar_sponge_" + std::to_string(ns)
                                    : "scalar_sponge_react_" + std::to_string(ns);
      par_for(label, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real rho = w0(m, IDN, k, j, i);
        const Real theta = w0(m, nhydro + ns, k, j, i);
        Real &x1min_ = size.d_view(m).x1min;
        Real &x1max_ = size.d_view(m).x1max;
        const Real x = CellCenterX(i - is, nx1, x1min_, x1max_);

        const Real chiL = cfg.chiL0 * SmoothTopHat(x, cfg.xL0, cfg.xL1, cfg.mask_w);
        const Real chiR = cfg.chiR0 * SmoothTopHat(x, cfg.xR0, cfg.xR1, cfg.mask_w);
        const Real sponge = -chiL * (theta - cfg.thetaL) - chiR * (theta - cfg.thetaR);

        Real reaction = 0.0;
        if (cfg.mode == kScalarSpongeReaction) {
          if (cfg.reaction_type == kReactionKpp) {
            Real theta_eff = theta;
            if (cfg.clip_reaction_interval) {
              if (theta_eff < cfg.theta1) theta_eff = cfg.theta1;
              if (theta_eff > cfg.theta3) theta_eff = cfg.theta3;
            }
            reaction = KppReaction(theta_eff, cfg.tau, cfg.theta1, cfg.theta3);
          } else {
            reaction = BistableReaction(theta, cfg.tau, cfg.theta1, cfg.theta2, cfg.theta3);
          }
        }

        u0(m, nhydro + ns, k, j, i) += rho * (sponge + reaction) * bdt;
      });
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn ScalarMixingHistory
//  \brief User history outputs: volume-integrated |grad theta_s|^2 for each scalar.

void ScalarMixingHistory(HistoryData *pdata, Mesh *pm) {
  auto *phydro = pm->pmb_pack->phydro;
  if (phydro == nullptr || phydro->nscalars == 0) {
    pdata->nhist = 0;
    for (int n = 0; n < NHISTORY_VARIABLES; ++n) {
      pdata->hdata[n] = 0.0;
    }
    return;
  }

  if (phydro->nscalars > NHISTORY_VARIABLES) {
    FatalProblemSetup("ScalarMixingHistory requires NHISTORY_VARIABLES >= nscalars.");
  }

  pdata->nhist = phydro->nscalars;
  for (int n = 0; n < pdata->nhist; ++n) {
    pdata->label[n] = "gradth2_s" + std::to_string(n);
  }

  auto &w0 = phydro->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  const int nhist = pdata->nhist;
  const int nhydro = phydro->nhydro;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int nx1 = indcs.nx1;
  const int js = indcs.js;
  const int nx2 = indcs.nx2;
  const int ks = indcs.ks;
  const int nx3 = indcs.nx3;
  const int nmkji = pm->pmb_pack->nmb_thispack * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("scalar_hist_grad2",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx / nkji;
    int k = (idx - m*nkji) / nji;
    int j = (idx - m*nkji - k*nji) / nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    const Real vol = dx1 * dx2 * dx3;

    array_sum::GlobalSum hvars;
    for (int ns = 0; ns < nhist; ++ns) {
      const Real dtheta_dx =
          (w0(m, nhydro + ns, k, j, i+1) - w0(m, nhydro + ns, k, j, i-1)) / (2.0 * dx1);
      Real dtheta_dy = 0.0;
      Real dtheta_dz = 0.0;
      if (multi_d) {
        dtheta_dy =
            (w0(m, nhydro + ns, k, j+1, i) - w0(m, nhydro + ns, k, j-1, i)) / (2.0 * dx2);
      }
      if (three_d) {
        dtheta_dz =
            (w0(m, nhydro + ns, k+1, j, i) - w0(m, nhydro + ns, k-1, j, i)) / (2.0 * dx3);
      }
      const Real grad2 = dtheta_dx*dtheta_dx + dtheta_dy*dtheta_dy + dtheta_dz*dtheta_dz;
      hvars.the_array[ns] = grad2 * vol;
    }

    for (int n = nhist; n < NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  for (int n = 0; n < pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
}
