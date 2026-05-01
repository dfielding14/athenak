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
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
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

enum class VelocityMethod { Projection, Stream2D, Clebsch };

struct TurbulenceConfig {
  VelocityMethod method = VelocityMethod::Projection;
  Real v_rms = 1.0;
  int nlow = 1;
  int nhigh = 4;
  Real expo = 5.0/3.0;
  Real alpha = 1.0/3.0;
  Real phi_slope = 11.0/3.0;
  Real velocity_slope = 5.0/3.0;
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
  bool clebsch_3d = false;
  bool active_v3 = false;
  bool dump_generator_diagnostics = false;
  bool legacy_stream_bool_used = false;
  bool clebsch_alpha_from_expo = false;

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

struct VelocityStats {
  Real vmean1 = 0.0;
  Real vmean2 = 0.0;
  Real vmean3 = 0.0;
  Real scale = 1.0;
  bool scalar_face_velocity = false;
  Real face_divergence_ratio = 0.0;
  Real face_divergence_linf = 0.0;
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

void RejectRemovedProblemParameter(ParameterInput *pin, const std::string &name,
                                   const std::string &replacement) {
  if (!pin->DoesParameterExist("problem", name)) return;
  FatalProblemSetup("<problem>/" + name + " was removed from scalar_mixing. "
                    + replacement);
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
  switch (method) {
    case VelocityMethod::Projection:
      return "projection";
    case VelocityMethod::Stream2D:
      return "stream_2d";
    case VelocityMethod::Clebsch:
      return "clebsch";
  }
  return "unknown";
}

Real VelocitySlopeFromAlpha(Real alpha) {
  return 2.0 * alpha + 1.0;
}

Real ClebschPhiSlopeFromAlpha(Real alpha) {
  return (alpha >= 0.0) ? (3.0 + 2.0 * alpha) : (3.0 + alpha);
}

Real UniformCellCenter(int index, int ncell, Real xmin, Real xmax) {
  if (ncell <= 0) return 0.5 * (xmin + xmax);
  return xmin + (static_cast<Real>(index) + 0.5) * (xmax - xmin) / static_cast<Real>(ncell);
}

KOKKOS_INLINE_FUNCTION
Real X1SineScalarProfile(Real x1, Real x1min, Real lx) {
  if (lx <= 0.0) return 0.0;
  const Real phase = 2.0 * M_PI * (x1 - x1min) / lx;
  return 0.5 * (1.0 - cos(phase));
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
  os << "\"clebsch_3d\": " << (cfg.clebsch_3d ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"legacy_stream_bool_used\": "
     << (cfg.legacy_stream_bool_used ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"clebsch_alpha_from_expo\": "
     << (cfg.clebsch_alpha_from_expo ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"scalar_face_velocity\": "
     << (stats.scalar_face_velocity ? "true" : "false") << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"face_divergence_ratio\": " << stats.face_divergence_ratio << ",\n";
  WriteIndent(os, indent + 2);
  os << "\"face_divergence_linf\": " << stats.face_divergence_linf << ",\n";
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
  os << "\"alpha\": " << cfg.alpha << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"phi_slope\": " << cfg.phi_slope << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"velocity_slope\": " << cfg.velocity_slope << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"rseed\": " << cfg.rseed << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"k_crit\": " << cfg.k_crit << ",\n";
  WriteIndent(os, indent + 4);
  os << "\"k_crit_mag\": " << cfg.k_crit_mag << "\n";
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

void WriteClebsch3DDiagnostics(const TurbulenceConfig &cfg,
                               const ModeCatalog &catalog,
                               const ScalarCoefficients &phi1_coeffs,
                               const ScalarCoefficients &phi2_coeffs,
                               const VelocityStats &stats,
                               Real init_wall_seconds) {
  if (global_variable::my_rank != 0 || !cfg.dump_generator_diagnostics) return;

  const int slice_k = cfg.global_nx3 / 2;
  const Real x3_slice = UniformCellCenter(slice_k, cfg.global_nx3, cfg.x3min, cfg.x3max);
  const std::vector<Real> phi1_xy =
      SampleScalarFieldXY(cfg, catalog, phi1_coeffs, cfg.global_nx1, cfg.global_nx2, x3_slice);
  const std::vector<Real> phi2_xy =
      SampleScalarFieldXY(cfg, catalog, phi2_coeffs, cfg.global_nx1, cfg.global_nx2, x3_slice);
  const std::vector<Real> phi1_shell =
      ComputeScalarShellEnergy(catalog, phi1_coeffs, cfg.nhigh);
  const std::vector<Real> phi2_shell =
      ComputeScalarShellEnergy(catalog, phi2_coeffs, cfg.nhigh);

  std::ofstream os(DiagnosticsPath(cfg));
  os << std::setprecision(17);
  os << "{\n";
  WriteIndent(os, 2);
  os << "\"metadata\": ";
  WriteDiagnosticsMetadata(os, cfg, "clebsch_3d", stats, init_wall_seconds, 2);
  os << ",\n";
  WriteIndent(os, 2);
  os << "\"catalogs\": {\n";
  WriteIndent(os, 4);
  os << "\"phi\": ";
  WriteJsonModeCatalog(os, catalog, 4);
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
  os << "\"phi1_shell_energy\": ";
  WriteJsonNumericArray(os, phi1_shell);
  os << ",\n";
  WriteIndent(os, 4);
  os << "\"phi2_shell_energy\": ";
  WriteJsonNumericArray(os, phi2_shell);
  os << "\n";
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

  RejectRemovedProblemParameter(pin, "turb_sol_frac",
                                "Projection mode is now always solenoidal.");
  RejectRemovedProblemParameter(pin, "divfree_scalar_flux",
                                "scalar_mixing now always uses the divergence-free "
                                "scalar face velocity.");

  cfg.v_rms = pin->GetOrAddReal("problem", "turb_v_rms", 1.0);
  cfg.nlow = pin->GetOrAddInteger("problem", "turb_nlow", 1);
  cfg.nhigh = pin->GetOrAddInteger("problem", "turb_nhigh", 4);
  cfg.expo = pin->GetOrAddReal("problem", "turb_expo", 5.0/3.0);
  cfg.rseed = GetTurbulenceSeed(pin);
  cfg.k_crit = pin->GetOrAddReal("problem", "turb_k_crit", 16.0);
  cfg.dump_generator_diagnostics =
      pin->GetOrAddBoolean("problem", "turb_dump_generator_diagnostics", false);
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

  if (pin->DoesParameterExist("problem", "turb_velocity_method")) {
    const std::string method_name =
        pin->GetOrAddString("problem", "turb_velocity_method", "projection");
    if (method_name == "projection") {
      cfg.method = VelocityMethod::Projection;
    } else if (method_name == "stream_2d") {
      cfg.method = VelocityMethod::Stream2D;
    } else if (method_name == "clebsch") {
      cfg.method = VelocityMethod::Clebsch;
    } else {
      FatalProblemSetup("Unknown <problem>/turb_velocity_method=" + method_name
                        + ". Expected projection, stream_2d, or clebsch.");
    }
  } else {
    const bool legacy_stream =
        pin->GetOrAddBoolean("problem", "turb_use_stream_function", false);
    if (legacy_stream) {
      cfg.legacy_stream_bool_used = true;
      cfg.method = (cfg.nx3 == 1) ? VelocityMethod::Stream2D : VelocityMethod::Clebsch;
    } else {
      cfg.method = VelocityMethod::Projection;
    }
  }

  cfg.stream_2d = (cfg.method == VelocityMethod::Stream2D);
  cfg.clebsch_3d = (cfg.method == VelocityMethod::Clebsch);
  if (cfg.method == VelocityMethod::Stream2D) {
    if (cfg.nx2 == 1 || cfg.nx3 != 1) {
      FatalProblemSetup("stream_2d initialization requires a true 2D mesh with nx3=1.");
    }
    cfg.active_v3 = false;
  } else if (cfg.method == VelocityMethod::Clebsch) {
    if (cfg.nx2 == 1 || cfg.nx3 == 1) {
      FatalProblemSetup("clebsch initialization requires a 3D mesh.");
    }
    cfg.active_v3 = true;
  } else {
    cfg.active_v3 = !projection_true_2d;
  }

  cfg.alpha = 0.5 * (cfg.expo - 1.0);
  if (cfg.clebsch_3d) {
    if (pin->DoesParameterExist("problem", "turb_alpha")) {
      cfg.alpha = pin->GetReal("problem", "turb_alpha");
    } else {
      cfg.clebsch_alpha_from_expo = true;
    }
    cfg.velocity_slope = VelocitySlopeFromAlpha(cfg.alpha);
    cfg.phi_slope = ClebschPhiSlopeFromAlpha(cfg.alpha);
    cfg.expo = cfg.velocity_slope;
  } else {
    cfg.velocity_slope = cfg.expo;
    cfg.phi_slope = 0.0;
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
      a[dir] -= k_vec[dir] * k_dot_a / kiso_sqr;
      b[dir] -= k_vec[dir] * k_dot_b / kiso_sqr;
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

ScalarCoefficients GenerateShellPowerLawScalarCoefficients(const ModeCatalog &catalog,
                                                           int shell_lo, int shell_hi,
                                                           Real slope,
                                                           RNG_State *rstate) {
  ScalarCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka.resize(nmodes);
  coeffs.akb.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real boost = catalog.boost[n];
    coeffs.aka[n] = boost * RanGaussianSt(rstate);
    coeffs.akb[n] = boost * RanGaussianSt(rstate);
  }

  const std::vector<Real> raw_shell = ComputeScalarShellEnergy(catalog, coeffs, shell_hi);
  std::vector<Real> shell_scale(shell_hi + 1, 0.0);
  for (int shell = shell_lo; shell <= shell_hi; ++shell) {
    if (raw_shell[shell] <= kTiny) {
      FatalProblemSetup("Cannot normalize Clebsch scalar shell "
                        + std::to_string(shell) + " because the retained shell energy is zero.");
    }
    shell_scale[shell] = std::sqrt(TargetShellEnergy(shell, slope) / raw_shell[shell]);
  }

  for (int n = 0; n < nmodes; ++n) {
    const Real scale = shell_scale[catalog.shell[n]];
    coeffs.aka[n] *= scale;
    coeffs.akb[n] *= scale;
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

VectorCoefficients ConvertScalarToGradientCoefficients(const ModeCatalog &catalog,
                                                       const ScalarCoefficients &phi_coeffs) {
  VectorCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka0.resize(nmodes);
  coeffs.aka1.resize(nmodes);
  coeffs.aka2.resize(nmodes);
  coeffs.akb0.resize(nmodes);
  coeffs.akb1.resize(nmodes);
  coeffs.akb2.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real a = phi_coeffs.aka[n];
    const Real b = phi_coeffs.akb[n];
    coeffs.aka0[n] = -catalog.kx[n] * b;
    coeffs.aka1[n] = -catalog.ky[n] * b;
    coeffs.aka2[n] = -catalog.kz[n] * b;
    coeffs.akb0[n] =  catalog.kx[n] * a;
    coeffs.akb1[n] =  catalog.ky[n] * a;
    coeffs.akb2[n] =  catalog.kz[n] * a;
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

void SynthesizeClebschVelocityFromGradientModes(
    MeshBlockPack *pmbp, const TurbulenceConfig &cfg, Real den,
    const ModeCatalog &catalog, const VectorCoefficients &grad1_coeffs,
    const VectorCoefficients &grad2_coeffs) {
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
  DvceArray1D<Real> d_g1a0("scalar_mix_g1a0", nmodes);
  DvceArray1D<Real> d_g1a1("scalar_mix_g1a1", nmodes);
  DvceArray1D<Real> d_g1a2("scalar_mix_g1a2", nmodes);
  DvceArray1D<Real> d_g1b0("scalar_mix_g1b0", nmodes);
  DvceArray1D<Real> d_g1b1("scalar_mix_g1b1", nmodes);
  DvceArray1D<Real> d_g1b2("scalar_mix_g1b2", nmodes);
  DvceArray1D<Real> d_g2a0("scalar_mix_g2a0", nmodes);
  DvceArray1D<Real> d_g2a1("scalar_mix_g2a1", nmodes);
  DvceArray1D<Real> d_g2a2("scalar_mix_g2a2", nmodes);
  DvceArray1D<Real> d_g2b0("scalar_mix_g2b0", nmodes);
  DvceArray1D<Real> d_g2b1("scalar_mix_g2b1", nmodes);
  DvceArray1D<Real> d_g2b2("scalar_mix_g2b2", nmodes);

  auto h_kx = Kokkos::create_mirror_view(d_kx);
  auto h_ky = Kokkos::create_mirror_view(d_ky);
  auto h_kz = Kokkos::create_mirror_view(d_kz);
  auto h_g1a0 = Kokkos::create_mirror_view(d_g1a0);
  auto h_g1a1 = Kokkos::create_mirror_view(d_g1a1);
  auto h_g1a2 = Kokkos::create_mirror_view(d_g1a2);
  auto h_g1b0 = Kokkos::create_mirror_view(d_g1b0);
  auto h_g1b1 = Kokkos::create_mirror_view(d_g1b1);
  auto h_g1b2 = Kokkos::create_mirror_view(d_g1b2);
  auto h_g2a0 = Kokkos::create_mirror_view(d_g2a0);
  auto h_g2a1 = Kokkos::create_mirror_view(d_g2a1);
  auto h_g2a2 = Kokkos::create_mirror_view(d_g2a2);
  auto h_g2b0 = Kokkos::create_mirror_view(d_g2b0);
  auto h_g2b1 = Kokkos::create_mirror_view(d_g2b1);
  auto h_g2b2 = Kokkos::create_mirror_view(d_g2b2);

  for (int n = 0; n < nmodes; ++n) {
    h_kx(n) = catalog.kx[n];
    h_ky(n) = catalog.ky[n];
    h_kz(n) = catalog.kz[n];
    h_g1a0(n) = grad1_coeffs.aka0[n];
    h_g1a1(n) = grad1_coeffs.aka1[n];
    h_g1a2(n) = grad1_coeffs.aka2[n];
    h_g1b0(n) = grad1_coeffs.akb0[n];
    h_g1b1(n) = grad1_coeffs.akb1[n];
    h_g1b2(n) = grad1_coeffs.akb2[n];
    h_g2a0(n) = grad2_coeffs.aka0[n];
    h_g2a1(n) = grad2_coeffs.aka1[n];
    h_g2a2(n) = grad2_coeffs.aka2[n];
    h_g2b0(n) = grad2_coeffs.akb0[n];
    h_g2b1(n) = grad2_coeffs.akb1[n];
    h_g2b2(n) = grad2_coeffs.akb2[n];
  }

  Kokkos::deep_copy(d_kx, h_kx);
  Kokkos::deep_copy(d_ky, h_ky);
  Kokkos::deep_copy(d_kz, h_kz);
  Kokkos::deep_copy(d_g1a0, h_g1a0);
  Kokkos::deep_copy(d_g1a1, h_g1a1);
  Kokkos::deep_copy(d_g1a2, h_g1a2);
  Kokkos::deep_copy(d_g1b0, h_g1b0);
  Kokkos::deep_copy(d_g1b1, h_g1b1);
  Kokkos::deep_copy(d_g1b2, h_g1b2);
  Kokkos::deep_copy(d_g2a0, h_g2a0);
  Kokkos::deep_copy(d_g2a1, h_g2a1);
  Kokkos::deep_copy(d_g2a2, h_g2a2);
  Kokkos::deep_copy(d_g2b0, h_g2b0);
  Kokkos::deep_copy(d_g2b1, h_g2b1);
  Kokkos::deep_copy(d_g2b2, h_g2b2);

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  const bool active_v3 = cfg.active_v3;

  par_for("scalar_mix_clebsch_vel_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
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

    Real g1x = 0.0;
    Real g1y = 0.0;
    Real g1z = 0.0;
    Real g2x = 0.0;
    Real g2y = 0.0;
    Real g2z = 0.0;
    for (int n = 0; n < nmodes; ++n) {
      const Real phase = d_kx(n)*x1v + d_ky(n)*x2v + d_kz(n)*x3v;
      const Real cosk = cos(phase);
      const Real sink = sin(phase);
      g1x += d_g1a0(n)*cosk - d_g1b0(n)*sink;
      g1y += d_g1a1(n)*cosk - d_g1b1(n)*sink;
      g1z += d_g1a2(n)*cosk - d_g1b2(n)*sink;
      g2x += d_g2a0(n)*cosk - d_g2b0(n)*sink;
      g2y += d_g2a1(n)*cosk - d_g2b1(n)*sink;
      g2z += d_g2a2(n)*cosk - d_g2b2(n)*sink;
    }

    const Real vx = g1y*g2z - g1z*g2y;
    const Real vy = g1z*g2x - g1x*g2z;
    const Real vz = g1x*g2y - g1y*g2x;

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
    const VectorCoefficients &apot_coeffs) {
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
  const int nmodes = catalog.KeptModes();

  if (nmodes == 0) return;

  auto &hydro = *pmbp->phydro;
  const bool multi_d = cfg.mesh_multi_d;
  const bool three_d = cfg.mesh_three_d;

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
    sface_x1f(m,k,j,i) = vx;
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
      sface_x2f(m,k,j,i) = vy;
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
      sface_x3f(m,k,j,i) = vz;
    });
  }
}

void EnsureScalarFaceVelocityField(MeshBlockPack *pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &hydro = *pmbp->phydro;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = pmbp->pmesh->multi_d ? indcs.nx2 + 2*indcs.ng : 1;
  const int ncells3 = pmbp->pmesh->three_d ? indcs.nx3 + 2*indcs.ng : 1;
  if (!hydro.scalar_vface) {
    hydro.scalar_vface = std::make_unique<DvceFaceFld4D<Real>>(
        "scalar_vface", pmbp->nmb_thispack, ncells3, ncells2, ncells1);
  }
  hydro.use_scalar_face_velocity = true;
}

void ZeroScalarFaceVelocities(MeshBlockPack *pmbp) {
  EnsureScalarFaceVelocityField(pmbp);
  auto &hydro = *pmbp->phydro;
  Kokkos::deep_copy(hydro.scalar_vface->x1f, 0.0);
  Kokkos::deep_copy(hydro.scalar_vface->x2f, 0.0);
  Kokkos::deep_copy(hydro.scalar_vface->x3f, 0.0);
}

void AverageScalarFaceVelocitiesToCellCenters(MeshBlockPack *pmbp,
                                              const TurbulenceConfig &cfg) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;

  auto &hydro = *pmbp->phydro;
  auto sface_x1f = hydro.scalar_vface->x1f;
  auto sface_x2f = hydro.scalar_vface->x2f;
  auto sface_x3f = hydro.scalar_vface->x3f;
  auto &u0 = hydro.u0;
  const bool multi_d = cfg.mesh_multi_d;
  const bool use_v3 = cfg.mesh_three_d && cfg.active_v3;

  par_for("scalar_face_to_cell_center", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real rho = u0(m, IDN, k, j, i);
    const Real v1 = 0.5 * (sface_x1f(m, k, j, i) + sface_x1f(m, k, j, i + 1));
    const Real v2 = multi_d
                        ? 0.5 * (sface_x2f(m, k, j, i) + sface_x2f(m, k, j + 1, i))
                        : 0.0;
    const Real v3 = use_v3
                        ? 0.5 * (sface_x3f(m, k, j, i) + sface_x3f(m, k + 1, j, i))
                        : 0.0;
    u0(m, IM1, k, j, i) = rho * v1;
    u0(m, IM2, k, j, i) = rho * v2;
    u0(m, IM3, k, j, i) = rho * v3;
  });
}

void ApplyVelocityStatsToScalarFaceVelocities(MeshBlockPack *pmbp,
                                              const TurbulenceConfig &cfg,
                                              const VelocityStats &stats) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ng = indcs.ng;
  const int nmb = pmbp->nmb_thispack;

  auto &hydro = *pmbp->phydro;
  auto sface_x1f = hydro.scalar_vface->x1f;
  auto sface_x2f = hydro.scalar_vface->x2f;
  auto sface_x3f = hydro.scalar_vface->x3f;
  const bool multi_d = cfg.mesh_multi_d;
  const bool use_v3 = cfg.mesh_three_d && cfg.active_v3;
  const bool three_d = cfg.mesh_three_d;

  const int x1_kl = three_d ? ks - ng : ks;
  const int x1_ku = three_d ? ke + ng : ke;
  const int x1_jl = multi_d ? js - ng : js;
  const int x1_ju = multi_d ? je + ng : je;
  par_for("scalar_face_apply_stats_x1", DevExeSpace(), 0, nmb-1,
          x1_kl, x1_ku, x1_jl, x1_ju, is - ng, ie + ng + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    sface_x1f(m, k, j, i) = (sface_x1f(m, k, j, i) - stats.vmean1) * stats.scale;
  });

  if (multi_d) {
    const int x2_kl = three_d ? ks - ng : ks;
    const int x2_ku = three_d ? ke + ng : ke;
    par_for("scalar_face_apply_stats_x2", DevExeSpace(), 0, nmb-1,
            x2_kl, x2_ku, js - ng, je + ng + 1, is - ng, ie + ng,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      sface_x2f(m, k, j, i) = (sface_x2f(m, k, j, i) - stats.vmean2) * stats.scale;
    });
  } else {
    Kokkos::deep_copy(hydro.scalar_vface->x2f, 0.0);
  }

  if (three_d) {
    par_for("scalar_face_apply_stats_x3", DevExeSpace(), 0, nmb-1,
            ks - ng, ke + ng + 1, js - ng, je + ng, is - ng, ie + ng,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      sface_x3f(m, k, j, i) = use_v3
                                  ? (sface_x3f(m, k, j, i) - stats.vmean3) * stats.scale
                                  : 0.0;
    });
  } else {
    Kokkos::deep_copy(hydro.scalar_vface->x3f, 0.0);
  }
}

void UpdateScalarFaceVelocityDiagnostics(MeshBlockPack *pmbp,
                                         const TurbulenceConfig &cfg,
                                         VelocityStats &stats) {
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
  const int nmb = pmbp->nmb_thispack;
  const int nmkji = nmb * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;

  auto &hydro = *pmbp->phydro;
  auto sface_x1f = hydro.scalar_vface->x1f;
  auto sface_x2f = hydro.scalar_vface->x2f;
  auto sface_x3f = hydro.scalar_vface->x3f;
  auto &size = pmbp->pmb->mb_size;
  const bool multi_d = cfg.mesh_multi_d;
  const bool three_d = cfg.mesh_three_d;

  Real sum_div2 = 0.0;
  Real sum_grad2 = 0.0;
  Real max_abs_div = 0.0;
  Kokkos::parallel_reduce("scalar_face_divergence_norm",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &s_div2, Real &s_grad2, Real &s_maxdiv) {
    const int m = idx / nkji;
    int k = (idx - m * nkji) / nji;
    int j = (idx - m * nkji - k * nji) / nx1;
    int i = (idx - m * nkji - k * nji - j * nx1) + is;
    k += ks;
    j += js;

    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    const Real d1 = (sface_x1f(m, k, j, i + 1) - sface_x1f(m, k, j, i)) / dx1;
    const Real d2 = multi_d
                        ? (sface_x2f(m, k, j + 1, i) - sface_x2f(m, k, j, i)) / dx2
                        : 0.0;
    const Real d3 = three_d
                        ? (sface_x3f(m, k + 1, j, i) - sface_x3f(m, k, j, i)) / dx3
                        : 0.0;
    const Real div = d1 + d2 + d3;
    s_div2 += div * div;
    s_grad2 += d1 * d1 + d2 * d2 + d3 * d3;
    const Real abs_div = fabs(div);
    if (abs_div > s_maxdiv) s_maxdiv = abs_div;
  }, Kokkos::Sum<Real>(sum_div2), Kokkos::Sum<Real>(sum_grad2),
     Kokkos::Max<Real>(max_abs_div));

#if MPI_PARALLEL_ENABLED
  Real local_sum[2] = {sum_div2, sum_grad2};
  Real global_sum[2];
  MPI_Allreduce(local_sum, global_sum, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  sum_div2 = global_sum[0];
  sum_grad2 = global_sum[1];
  Real global_max_div = 0.0;
  MPI_Allreduce(&max_abs_div, &global_max_div, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  max_abs_div = global_max_div;
#endif

  stats.scalar_face_velocity = true;
  stats.face_divergence_ratio = (sum_grad2 > kTiny) ? std::sqrt(sum_div2 / sum_grad2) : 0.0;
  stats.face_divergence_linf = max_abs_div;
}

KOKKOS_INLINE_FUNCTION
void AccumulateScalarFieldAndGradient(const DvceArray1D<Real> &d_kx,
                                      const DvceArray1D<Real> &d_ky,
                                      const DvceArray1D<Real> &d_kz,
                                      const DvceArray1D<Real> &d_aka,
                                      const DvceArray1D<Real> &d_akb,
                                      int mode_count, Real x1v, Real x2v, Real x3v,
                                      Real &phi, Real &d1phi, Real &d2phi, Real &d3phi) {
  for (int n = 0; n < mode_count; ++n) {
    const Real phase = d_kx(n) * x1v + d_ky(n) * x2v + d_kz(n) * x3v;
    const Real cosk = cos(phase);
    const Real sink = sin(phase);
    const Real a = d_aka(n);
    const Real b = d_akb(n);
    const Real sin_term = a * sink + b * cosk;
    phi += a * cosk - b * sink;
    d1phi -= d_kx(n) * sin_term;
    d2phi -= d_ky(n) * sin_term;
    d3phi -= d_kz(n) * sin_term;
  }
}

void PopulateScalarFaceVelocitiesFromClebschPotential(
    MeshBlockPack *pmbp, const TurbulenceConfig &cfg, const ModeCatalog &catalog,
    const ScalarCoefficients &phi1_coeffs, const ScalarCoefficients &phi2_coeffs) {
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
  const int nmodes = catalog.KeptModes();

  if (nmodes == 0) return;

  DvceArray1D<Real> d_kx("scalar_clebsch_kx", nmodes);
  DvceArray1D<Real> d_ky("scalar_clebsch_ky", nmodes);
  DvceArray1D<Real> d_kz("scalar_clebsch_kz", nmodes);
  DvceArray1D<Real> d_phi1a("scalar_clebsch_phi1a", nmodes);
  DvceArray1D<Real> d_phi1b("scalar_clebsch_phi1b", nmodes);
  DvceArray1D<Real> d_phi2a("scalar_clebsch_phi2a", nmodes);
  DvceArray1D<Real> d_phi2b("scalar_clebsch_phi2b", nmodes);

  auto h_kx = Kokkos::create_mirror_view(d_kx);
  auto h_ky = Kokkos::create_mirror_view(d_ky);
  auto h_kz = Kokkos::create_mirror_view(d_kz);
  auto h_phi1a = Kokkos::create_mirror_view(d_phi1a);
  auto h_phi1b = Kokkos::create_mirror_view(d_phi1b);
  auto h_phi2a = Kokkos::create_mirror_view(d_phi2a);
  auto h_phi2b = Kokkos::create_mirror_view(d_phi2b);

  for (int n = 0; n < nmodes; ++n) {
    h_kx(n) = catalog.kx[n];
    h_ky(n) = catalog.ky[n];
    h_kz(n) = catalog.kz[n];
    h_phi1a(n) = phi1_coeffs.aka[n];
    h_phi1b(n) = phi1_coeffs.akb[n];
    h_phi2a(n) = phi2_coeffs.aka[n];
    h_phi2b(n) = phi2_coeffs.akb[n];
  }

  Kokkos::deep_copy(d_kx, h_kx);
  Kokkos::deep_copy(d_ky, h_ky);
  Kokkos::deep_copy(d_kz, h_kz);
  Kokkos::deep_copy(d_phi1a, h_phi1a);
  Kokkos::deep_copy(d_phi1b, h_phi1b);
  Kokkos::deep_copy(d_phi2a, h_phi2a);
  Kokkos::deep_copy(d_phi2b, h_phi2b);

  const int ncells1 = nx1 + 2 * ng;
  const int ncells2 = nx2 + 2 * ng;
  const int ncells3 = nx3 + 2 * ng;
  DvceEdgeFld4D<Real> apot("scalar_clebsch_apot", nmb, ncells3, ncells2, ncells1);
  auto apot_x1e = apot.x1e;
  auto apot_x2e = apot.x2e;
  auto apot_x3e = apot.x3e;

  auto &size = pmbp->pmb->mb_size;
  auto &hydro = *pmbp->phydro;
  auto sface_x1f = hydro.scalar_vface->x1f;
  auto sface_x2f = hydro.scalar_vface->x2f;
  auto sface_x3f = hydro.scalar_vface->x3f;

  par_for("scalar_clebsch_apot_x1", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng + 1, js - ng, je + ng + 1, is - ng, ie + ng,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    const Real x1v = CellCenterX(i - is, nx1, x1min, x1max);
    const Real x2v = LeftEdgeX(j - js, nx2, x2min, x2max);
    const Real x3v = LeftEdgeX(k - ks, nx3, x3min, x3max);
    Real phi1 = 0.0, d1phi1 = 0.0, d2phi1 = 0.0, d3phi1 = 0.0;
    Real phi2 = 0.0, d1phi2 = 0.0, d2phi2 = 0.0, d3phi2 = 0.0;
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi1a, d_phi1b, nmodes,
                                     x1v, x2v, x3v, phi1, d1phi1, d2phi1, d3phi1);
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi2a, d_phi2b, nmodes,
                                     x1v, x2v, x3v, phi2, d1phi2, d2phi2, d3phi2);
    apot_x1e(m, k, j, i) = 0.5 * (phi1 * d1phi2 - phi2 * d1phi1);
  });

  par_for("scalar_clebsch_apot_x2", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng + 1, js - ng, je + ng, is - ng, ie + ng + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    const Real x1v = LeftEdgeX(i - is, nx1, x1min, x1max);
    const Real x2v = CellCenterX(j - js, nx2, x2min, x2max);
    const Real x3v = LeftEdgeX(k - ks, nx3, x3min, x3max);
    Real phi1 = 0.0, d1phi1 = 0.0, d2phi1 = 0.0, d3phi1 = 0.0;
    Real phi2 = 0.0, d1phi2 = 0.0, d2phi2 = 0.0, d3phi2 = 0.0;
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi1a, d_phi1b, nmodes,
                                     x1v, x2v, x3v, phi1, d1phi1, d2phi1, d3phi1);
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi2a, d_phi2b, nmodes,
                                     x1v, x2v, x3v, phi2, d1phi2, d2phi2, d3phi2);
    apot_x2e(m, k, j, i) = 0.5 * (phi1 * d2phi2 - phi2 * d2phi1);
  });

  par_for("scalar_clebsch_apot_x3", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng, js - ng, je + ng + 1, is - ng, ie + ng + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    const Real x1v = LeftEdgeX(i - is, nx1, x1min, x1max);
    const Real x2v = LeftEdgeX(j - js, nx2, x2min, x2max);
    const Real x3v = CellCenterX(k - ks, nx3, x3min, x3max);
    Real phi1 = 0.0, d1phi1 = 0.0, d2phi1 = 0.0, d3phi1 = 0.0;
    Real phi2 = 0.0, d1phi2 = 0.0, d2phi2 = 0.0, d3phi2 = 0.0;
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi1a, d_phi1b, nmodes,
                                     x1v, x2v, x3v, phi1, d1phi1, d2phi1, d3phi1);
    AccumulateScalarFieldAndGradient(d_kx, d_ky, d_kz, d_phi2a, d_phi2b, nmodes,
                                     x1v, x2v, x3v, phi2, d1phi2, d2phi2, d3phi2);
    apot_x3e(m, k, j, i) = 0.5 * (phi1 * d3phi2 - phi2 * d3phi1);
  });

  par_for("scalar_clebsch_face_x1", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng, js - ng, je + ng, is - ng, ie + ng + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    sface_x1f(m, k, j, i) =
        (apot_x3e(m, k, j + 1, i) - apot_x3e(m, k, j, i)) / size.d_view(m).dx2
        - (apot_x2e(m, k + 1, j, i) - apot_x2e(m, k, j, i)) / size.d_view(m).dx3;
  });

  par_for("scalar_clebsch_face_x2", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng, js - ng, je + ng + 1, is - ng, ie + ng,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    sface_x2f(m, k, j, i) =
        (apot_x1e(m, k + 1, j, i) - apot_x1e(m, k, j, i)) / size.d_view(m).dx3
        - (apot_x3e(m, k, j, i + 1) - apot_x3e(m, k, j, i)) / size.d_view(m).dx1;
  });

  par_for("scalar_clebsch_face_x3", DevExeSpace(), 0, nmb-1,
          ks - ng, ke + ng + 1, js - ng, je + ng, is - ng, ie + ng,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    sface_x3f(m, k, j, i) =
        (apot_x2e(m, k, j, i + 1) - apot_x2e(m, k, j, i)) / size.d_view(m).dx1
        - (apot_x1e(m, k, j + 1, i) - apot_x1e(m, k, j, i)) / size.d_view(m).dx2;
  });
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
      FatalProblemSetup("Velocity generator has no accepted modes in shell "
                        + std::to_string(shell) + " for " + label + ".");
    }
  }
}

int MaxSupportedShell(const TurbulenceConfig &cfg) {
  int max_supported_shell = cfg.global_nx1 / 2;
  if (cfg.mesh_multi_d) {
    max_supported_shell = std::min(max_supported_shell, cfg.global_nx2 / 2);
  }
  if (cfg.mesh_three_d) {
    max_supported_shell = std::min(max_supported_shell, cfg.global_nx3 / 2);
  }
  return std::max(1, std::min(cfg.nhigh, max_supported_shell));
}

void LogModeSummary(const TurbulenceConfig &cfg, const ModeCatalog &catalog,
                    const std::string &label) {
  std::cout << "  " << label << ": total_modes_in_shell=" << catalog.total_modes
            << ", kept_via_importance_sampling=" << catalog.KeptModes() << std::endl;
  if (cfg.method == VelocityMethod::Stream2D && cfg.stream_2d) {
    std::cout << "  " << label << ": stream_2d coefficients preserve "
              << "E_u(k) ~ k^(-expo) through |psi_k| ~ k^{-(expo+3)/2}" << std::endl;
  }
}

void InitTurbulentVelocity(MeshBlockPack *pmbp, ParameterInput *pin, Real den,
                           bool projection_true_2d) {
  (void)den;
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
  if (!hydro.scalar_only) {
    FatalProblemSetup("scalar_mixing now requires <hydro>/scalar_only = true because "
                      "the problem always uses the divergence-free scalar face velocity.");
  }
  if (cfg.legacy_stream_bool_used && global_variable::my_rank == 0) {
    std::cout << "  NOTE: turb_use_stream_function is deprecated; prefer "
              << "turb_velocity_method=stream_2d or turb_velocity_method=clebsch."
              << std::endl;
  }

  ZeroScalarFaceVelocities(pmbp);

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
                << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo << std::endl
                << "  true_2d=" << (cfg.projection_true_2d ? "true" : "false")
                << ", active_velocity_components=" << (cfg.active_v3 ? 3 : 2)
                << std::endl
                << "  k_crit=" << cfg.k_crit
                << " (wavenumber=" << cfg.k_crit_mag << ")" << std::endl;
      LogModeSummary(cfg, catalog, "velocity");
    }

    const VectorCoefficients coeffs = GenerateProjectionCoefficients(cfg, catalog, &rstate);
    const VectorCoefficients apot_coeffs =
        ConvertVelocityToVectorPotentialCoefficients(catalog, coeffs);
    PopulateScalarFaceVelocitiesFromVectorPotential(pmbp, cfg, catalog, apot_coeffs);
    AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
    VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
    ApplyVelocityStatsToScalarFaceVelocities(pmbp, cfg, stats);
    AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
    UpdateScalarFaceVelocityDiagnostics(pmbp, cfg, stats);
    const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
    recorded_init_wall_seconds = init_wall_seconds;
    WriteProjectionDiagnostics(cfg, catalog, coeffs, stats, init_wall_seconds);
    log_init_time();
    return;
  }

  if (global_variable::my_rank == 0) {
    std::cout << "Initializing turbulent velocity field:" << std::endl
              << "  velocity_method=" << VelocityMethodName(cfg.method) << std::endl
              << "  v_rms=" << cfg.v_rms << ", nlow=" << cfg.nlow
              << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo << std::endl
              << "  active_velocity_components=" << (cfg.active_v3 ? 3 : 2) << std::endl
              << "  k_crit=" << cfg.k_crit
              << " (wavenumber=" << cfg.k_crit_mag << ")" << std::endl;
    if (cfg.clebsch_3d) {
      std::cout << "  alpha=" << cfg.alpha
                << ", velocity_slope=" << cfg.velocity_slope
                << ", phi_slope=" << cfg.phi_slope << std::endl;
      if (cfg.clebsch_alpha_from_expo) {
        std::cout << "  NOTE: clebsch now prefers turb_alpha; using alpha=(turb_expo-1)/2 "
                  << "for compatibility." << std::endl;
      }
    }
  }

  if (cfg.method == VelocityMethod::Stream2D) {
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
    const VectorCoefficients apot_coeffs =
        ConvertVelocityToVectorPotentialCoefficients(catalog, vel_coeffs);
    PopulateScalarFaceVelocitiesFromVectorPotential(pmbp, cfg, catalog, apot_coeffs);
    AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
    VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
    ApplyVelocityStatsToScalarFaceVelocities(pmbp, cfg, stats);
    AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
    UpdateScalarFaceVelocityDiagnostics(pmbp, cfg, stats);
    const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
    recorded_init_wall_seconds = init_wall_seconds;
    WriteStream2DDiagnostics(cfg, catalog, psi_coeffs, vel_coeffs, stats, init_wall_seconds);
    log_init_time();
    return;
  }

  RNG_State catalog_state;
  catalog_state.idum = -MixSeed(cfg.rseed, 1);

  if (MaxSupportedShell(cfg) < cfg.nhigh) {
    FatalProblemSetup("clebsch mode requires turb_nhigh to lie within the resolved Fourier band.");
  }

  const ModeCatalog catalog = BuildModeCatalog(cfg, cfg.nlow, cfg.nhigh, &catalog_state);
  if (catalog.KeptModes() == 0) {
    std::cout << "### WARNING: No turbulent modes kept in range [" << cfg.nlow << ", "
              << cfg.nhigh << "] for clebsch initialization. Using zero velocity."
              << std::endl;
    log_init_time();
    return;
  }
  const std::vector<int> shell_counts = CountModesPerShell(catalog, cfg.nhigh);
  EnsureShellCoverage(shell_counts, cfg.nlow, cfg.nhigh, "phi");

  RNG_State coeff_state1;
  RNG_State coeff_state2;
  coeff_state1.idum = -MixSeed(cfg.rseed, 101);
  coeff_state2.idum = -MixSeed(cfg.rseed, 102);
  const ScalarCoefficients phi1_coeffs =
      GenerateShellPowerLawScalarCoefficients(catalog, cfg.nlow, cfg.nhigh,
                                              cfg.phi_slope, &coeff_state1);
  const ScalarCoefficients phi2_coeffs =
      GenerateShellPowerLawScalarCoefficients(catalog, cfg.nlow, cfg.nhigh,
                                              cfg.phi_slope, &coeff_state2);

  if (global_variable::my_rank == 0) {
    std::cout << "  stream_geometry=3D (Clebsch form v = grad(phi1) x grad(phi2))"
              << std::endl;
    LogModeSummary(cfg, catalog, "phi(shared)");
    std::cout << "  clebsch scalar spectra:"
              << " phi_target_slope=" << cfg.phi_slope
              << ", target_velocity_slope=" << cfg.velocity_slope
              << ", realized velocity spectrum is measured from hydro output FFTs"
              << std::endl;
  }
  PopulateScalarFaceVelocitiesFromClebschPotential(pmbp, cfg, catalog,
                                                   phi1_coeffs, phi2_coeffs);
  AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
  VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
  ApplyVelocityStatsToScalarFaceVelocities(pmbp, cfg, stats);
  AverageScalarFaceVelocitiesToCellCenters(pmbp, cfg);
  UpdateScalarFaceVelocityDiagnostics(pmbp, cfg, stats);
  const Real init_wall_seconds = GlobalMaxReal(ElapsedWallSeconds(init_start_time));
  recorded_init_wall_seconds = init_wall_seconds;
  WriteClebsch3DDiagnostics(cfg, catalog, phi1_coeffs, phi2_coeffs,
                            stats, init_wall_seconds);
  log_init_time();
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar mixing
//
//  Initializes density, turbulent velocity, and passive scalars.
//  Scalars start from a shared box-scale x1 sine profile.

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

  RejectRemovedProblemParameter(pin, "scalar_init",
                                "The scalar base state is now fixed to the box-scale "
                                "x1 sine profile.");
  RejectRemovedProblemParameter(pin, "scalar_use_x1_step",
                                "The scalar base state is now fixed to the box-scale "
                                "x1 sine profile.");
  RejectRemovedProblemParameter(pin, "scalar_step_x1",
                                "The scalar base state is now fixed to the box-scale "
                                "x1 sine profile.");
  RejectRemovedProblemParameter(pin, "scalar_left",
                                "The scalar base state is now fixed to the box-scale "
                                "x1 sine profile.");
  RejectRemovedProblemParameter(pin, "scalar_right",
                                "The scalar base state is now fixed to the box-scale "
                                "x1 sine profile.");

  const Real legacy_G = pin->GetOrAddReal("problem", "mean_gradient", 1.0);
  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  const Real lx = x1max - x1min;

  scalar_forcing_cfg.resize(nscalars);
  const Real scalar_noise_default = pin->GetOrAddReal("problem", "scalar_noise", 0.0);
  std::vector<Real> scalar_noise(nscalars, scalar_noise_default);
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

    RejectRemovedProblemParameter(pin, prefix + "init",
                                  "The scalar base state is now fixed to the box-scale "
                                  "x1 sine profile.");

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
              << "  scalar_ic = x1_sine" << std::endl
              << "  true_2d = " << (true_2d ? "true" : "false") << std::endl
              << "  nscalars = " << nscalars << std::endl
              << "  scalar_sine_period = " << lx << std::endl
              << "  scalar0_mode = " << scalar_forcing_cfg[0].mode << std::endl
              << "  scalar0_mean_gradient = " << scalar_forcing_cfg[0].G << std::endl;
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
    const Real noise = scalar_noise[ns];
    const Real theta_floor = scalar_forcing_cfg[ns].theta_floor;
    const std::string label = "scalar_mix_pgen_scalar_" + std::to_string(ns);
    par_for(label, DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real rho = u0(m, IDN, k, j, i);
      Real &x1min_ = size.d_view(m).x1min;
      const Real x1v = CellCenterX(i-is, nx1, x1min_, size.d_view(m).x1max);
      Real theta = X1SineScalarProfile(x1v, x1min, lx);
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
