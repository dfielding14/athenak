//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_mixing_perfect_powerlaw.cpp
//  \brief Scalar mixing with exact centered-shell velocity power-law initialization
//
//  Passive scalars with selectable source terms (mean gradient, reaction, sponge).
//  This simulates mixing of background gradients and/or reactive reservoirs by turbulence.
//
//  When scalar_only=true in <hydro>, velocity is frozen and only the scalars evolve.

#include <algorithm>
#include <cmath>
#include <cstdlib>
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
constexpr Real kPi = 3.141592653589793238462643383279502884L;
constexpr Real kDefaultScalar0WavelengthsPerBox = 1.0;
constexpr Real kDefaultScalar1WavelengthsPerBox = 2.0;

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
  std::string spectrum_contract = "exact_shell";

  int nx1 = 1;
  int nx2 = 1;
  int nx3 = 1;

  bool mesh_multi_d = false;
  bool mesh_three_d = false;
  bool projection_true_2d = false;
  bool stream_2d = false;
  bool active_v3 = false;
  bool divfree_scalar_flux = false;

  Real lx = 1.0;
  Real ly = 1.0;
  Real lz = 1.0;
  Real dkx = 0.0;
  Real dky = 0.0;
  Real dkz = 0.0;
  Real k_crit_mag = 0.0;

  bool exact_shell_contract = true;
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

struct CandidateMode {
  int nkx = 0;
  int nky = 0;
  int nkz = 0;
  int shell = 0;
  Real kx = 0.0;
  Real ky = 0.0;
  Real kz = 0.0;
  Real kiso = 0.0;
  Real p_accept = 1.0;
};

struct SignedMode {
  int nkx, nky, nkz;
  Real kx, ky, kz;
  Real kiso;
  Real variance_weight;
};

struct ClebschPairSamples {
  std::vector<int> idx1, idx2, shell;
  std::vector<Real> geom;

  int Count() const { return static_cast<int>(shell.size()); }
};

struct Stream3DBetaCalibration {
  Real asymptotic_beta = 0.0;
  Real asymptotic_slope = 0.0;
  Real calibrated_beta = 0.0;
  Real calibrated_slope = 0.0;
  int sample_count = 0;
  bool bracketed = false;
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
  return static_cast<int>(std::floor(std::sqrt(static_cast<Real>(nsqr)) + 0.5 - 1.0e-12));
}

Real FitLogSlope(const std::vector<Real> &spectrum, int nlow, int nhigh) {
  Real sx = 0.0;
  Real sy = 0.0;
  Real sxx = 0.0;
  Real sxy = 0.0;
  int npts = 0;
  for (int shell = nlow; shell <= nhigh; ++shell) {
    if (shell < 1 || shell >= static_cast<int>(spectrum.size())) continue;
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

KOKKOS_INLINE_FUNCTION
Real NormalizeBoxCoordinate(Real x, Real xmin, Real xmax) {
  const Real length = xmax - xmin;
  return (fabs(length) > kTiny) ? (x - xmin) / length : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real SineScalarIC1D(Real x1v,
                    Real x1min, Real x1max,
                    Real wavelengths_per_box) {
  const Real mode = 2.0 * kPi * wavelengths_per_box;
  const Real xi1 = NormalizeBoxCoordinate(x1v, x1min, x1max);
  return 0.5 * (1.0 + sin(mode * xi1));
}

KOKKOS_INLINE_FUNCTION
Real SineScalarIC(Real x1v, Real x2v, Real x3v,
                  Real x1min, Real x1max,
                  Real x2min, Real x2max,
                  Real x3min, Real x3max,
                  Real wavelengths_per_box,
                  bool multi_d, bool three_d) {
  const Real mode = 2.0 * kPi * wavelengths_per_box;
  const Real xi1 = NormalizeBoxCoordinate(x1v, x1min, x1max);
  Real pattern = sin(mode * xi1);
  if (multi_d) {
    const Real xi2 = NormalizeBoxCoordinate(x2v, x2min, x2max);
    pattern *= sin(mode * xi2);
  }
  if (three_d) {
    const Real xi3 = NormalizeBoxCoordinate(x3v, x3min, x3max);
    pattern *= sin(mode * xi3);
  }
  return 0.5 * (1.0 + pattern);
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
  cfg.spectrum_contract =
      pin->GetOrAddString("problem", "turb_spectrum_contract", "exact_shell");
  cfg.divfree_scalar_flux = pin->GetOrAddBoolean("problem", "divfree_scalar_flux", false);

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
  cfg.dkx = 2.0*M_PI/cfg.lx;
  cfg.dky = (cfg.nx2 > 1) ? 2.0*M_PI/cfg.ly : 0.0;
  cfg.dkz = (cfg.nx3 > 1) ? 2.0*M_PI/cfg.lz : 0.0;
  cfg.k_crit_mag = cfg.k_crit * cfg.dkx;
  cfg.exact_shell_contract = (cfg.spectrum_contract == "exact_shell");

  if (!(cfg.spectrum_contract == "exact_shell" || cfg.spectrum_contract == "statistical")) {
    FatalProblemSetup("<problem>/turb_spectrum_contract must be 'exact_shell' or 'statistical'.");
  }

  return cfg;
}

Real ModeAcceptanceProbability(const TurbulenceConfig &cfg, Real kiso) {
  if (kiso < cfg.k_crit_mag) return 1.0;
  return (cfg.k_crit_mag * cfg.k_crit_mag) / (kiso * kiso);
}

bool SupportsExactShellContract(const TurbulenceConfig &cfg) {
  if (cfg.method == VelocityMethod::StreamFunction) return cfg.stream_2d;
  if (cfg.method == VelocityMethod::Projection) {
    return cfg.mesh_three_d && cfg.active_v3 && std::abs(cfg.sol_frac - 1.0) <= 1.0e-12;
  }
  return false;
}

void ValidateExactShellContract(const TurbulenceConfig &cfg) {
  if (!cfg.exact_shell_contract) return;
  if (SupportsExactShellContract(cfg)) return;

  if (cfg.method == VelocityMethod::StreamFunction) {
    FatalProblemSetup("exact_shell is only supported for 2D stream-function initialization.");
  }
  FatalProblemSetup(
      "exact_shell is only supported for 3D projection initialization with turb_sol_frac=1.");
}

void AppendCandidateMode(const CandidateMode &candidate, ModeCatalog &catalog) {
  catalog.nkx.push_back(candidate.nkx);
  catalog.nky.push_back(candidate.nky);
  catalog.nkz.push_back(candidate.nkz);
  catalog.shell.push_back(candidate.shell);
  catalog.prob.push_back(candidate.p_accept);
  catalog.boost.push_back(1.0);
  catalog.kx.push_back(candidate.kx);
  catalog.ky.push_back(candidate.ky);
  catalog.kz.push_back(candidate.kz);
  catalog.kiso.push_back(candidate.kiso);
}

ModeCatalog BuildExactShellModeCatalog(const TurbulenceConfig &cfg, RNG_State *rstate) {
  ModeCatalog catalog;
  std::vector<std::vector<CandidateMode>> shell_candidates(cfg.nhigh + 1);

  for (int nkx = 0; nkx <= cfg.nhigh; ++nkx) {
    for (int nky = (cfg.mesh_multi_d ? -cfg.nhigh : 0);
         nky <= (cfg.mesh_multi_d ? cfg.nhigh : 0); ++nky) {
      for (int nkz = (cfg.mesh_three_d ? -cfg.nhigh : 0);
           nkz <= (cfg.mesh_three_d ? cfg.nhigh : 0); ++nkz) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;

        // Keep half of the kx=0 plane to avoid double-counting the implicit conjugates.
        if (nkx == 0) {
          if (nky < 0) continue;
          if (nky == 0 && nkz <= 0) continue;
        }

        const int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        const int shell = ShellFromIndexSqr(nsqr);
        if (shell < cfg.nlow || shell > cfg.nhigh) continue;

        const CandidateMode candidate{
            nkx, nky, nkz, shell, cfg.dkx * nkx, cfg.dky * nky, cfg.dkz * nkz,
            std::sqrt((cfg.dkx * nkx)*(cfg.dkx * nkx) +
                      (cfg.dky * nky)*(cfg.dky * nky) +
                      (cfg.dkz * nkz)*(cfg.dkz * nkz)),
            1.0};
        shell_candidates[shell].push_back(candidate);
      }
    }
  }

  for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
    auto &candidates = shell_candidates[shell];
    if (candidates.empty()) {
      FatalProblemSetup("exact_shell requires at least one Fourier mode in every in-band shell.");
    }
    catalog.total_modes += static_cast<int>(candidates.size());
    for (auto &candidate : candidates) {
      candidate.p_accept = ModeAcceptanceProbability(cfg, candidate.kiso);
    }

    Real expected_count = 0.0;
    for (const auto &candidate : candidates) {
      expected_count += candidate.p_accept;
    }
    const int quota = std::min(static_cast<int>(candidates.size()),
                               std::max(1, static_cast<int>(std::llround(expected_count))));

    if (quota >= static_cast<int>(candidates.size())) {
      std::sort(candidates.begin(), candidates.end(),
                [](const CandidateMode &lhs, const CandidateMode &rhs) {
                  if (lhs.nkx != rhs.nkx) return lhs.nkx < rhs.nkx;
                  if (lhs.nky != rhs.nky) return lhs.nky < rhs.nky;
                  return lhs.nkz < rhs.nkz;
                });
      for (const auto &candidate : candidates) {
        AppendCandidateMode(candidate, catalog);
      }
      continue;
    }

    std::vector<std::pair<Real, int>> weighted_keys;
    weighted_keys.reserve(candidates.size());
    for (int idx = 0; idx < static_cast<int>(candidates.size()); ++idx) {
      const Real u = std::max(RanSt(rstate), 1.0e-12);
      weighted_keys.emplace_back(std::log(u) / candidates[idx].p_accept, idx);
    }
    std::sort(weighted_keys.begin(), weighted_keys.end(),
              [](const std::pair<Real, int> &lhs, const std::pair<Real, int> &rhs) {
                return lhs.first > rhs.first;
              });

    std::vector<CandidateMode> selected;
    selected.reserve(quota);
    for (int idx = 0; idx < quota; ++idx) {
      selected.push_back(candidates[weighted_keys[idx].second]);
    }
    std::sort(selected.begin(), selected.end(),
              [](const CandidateMode &lhs, const CandidateMode &rhs) {
                if (lhs.nkx != rhs.nkx) return lhs.nkx < rhs.nkx;
                if (lhs.nky != rhs.nky) return lhs.nky < rhs.nky;
                return lhs.nkz < rhs.nkz;
              });
    for (const auto &candidate : selected) {
      AppendCandidateMode(candidate, catalog);
    }
  }

  return catalog;
}

ModeCatalog BuildModeCatalog(const TurbulenceConfig &cfg, RNG_State *rstate) {
  if (cfg.exact_shell_contract) {
    return BuildExactShellModeCatalog(cfg, rstate);
  }

  ModeCatalog catalog;
  const int nlow_sqr = cfg.nlow * cfg.nlow;
  const int nhigh_sqr = cfg.nhigh * cfg.nhigh;

  for (int nkx = 0; nkx <= cfg.nhigh; ++nkx) {
    for (int nky = (cfg.mesh_multi_d ? -cfg.nhigh : 0);
         nky <= (cfg.mesh_multi_d ? cfg.nhigh : 0); ++nky) {
      for (int nkz = (cfg.mesh_three_d ? -cfg.nhigh : 0);
           nkz <= (cfg.mesh_three_d ? cfg.nhigh : 0); ++nkz) {
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

  if (cfg.exact_shell_contract) {
    std::vector<std::vector<int>> shell_members(cfg.nhigh + 1);
    for (int n = 0; n < nmodes; ++n) {
      shell_members[catalog.shell[n]].push_back(n);
    }

    for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
      const auto &members = shell_members[shell];
      if (members.empty()) continue;

      std::vector<Real> weights(members.size(), 0.0);
      Real weight_sum = 0.0;
      for (int idx = 0; idx < static_cast<int>(members.size()); ++idx) {
        const int n = members[idx];
        const Real weight = 1.0 / std::pow(catalog.kiso[n], cfg.expo + 2.0);
        weights[idx] = weight;
        weight_sum += weight;

        const Real k_vec[3] = {catalog.kx[n], catalog.ky[n], catalog.kz[n]};
        Real a[3] = {0.0, 0.0, 0.0};
        Real b[3] = {0.0, 0.0, 0.0};
        Real mode_energy = 0.0;

        for (int attempt = 0; attempt < 64; ++attempt) {
          for (int dir = 0; dir < 3; ++dir) {
            a[dir] = RanGaussianSt(rstate);
            b[dir] = RanGaussianSt(rstate);
          }
          const Real kiso_sqr = catalog.kiso[n] * catalog.kiso[n];
          const Real k_dot_a = k_vec[0]*a[0] + k_vec[1]*a[1] + k_vec[2]*a[2];
          const Real k_dot_b = k_vec[0]*b[0] + k_vec[1]*b[1] + k_vec[2]*b[2];
          for (int dir = 0; dir < 3; ++dir) {
            a[dir] -= k_vec[dir] * k_dot_a / kiso_sqr;
            b[dir] -= k_vec[dir] * k_dot_b / kiso_sqr;
          }
          mode_energy = 0.5 * (a[0]*a[0] + a[1]*a[1] + a[2]*a[2] +
                               b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
          if (mode_energy > kTiny) break;
        }
        if (mode_energy <= kTiny) {
          FatalProblemSetup("Failed to draw a non-degenerate solenoidal projection mode.");
        }

        const Real inv_norm = 1.0 / std::sqrt(mode_energy);
        coeffs.aka0[n] = a[0] * inv_norm;
        coeffs.aka1[n] = a[1] * inv_norm;
        coeffs.aka2[n] = a[2] * inv_norm;
        coeffs.akb0[n] = b[0] * inv_norm;
        coeffs.akb1[n] = b[1] * inv_norm;
        coeffs.akb2[n] = b[2] * inv_norm;
      }

      const Real shell_target = std::pow(static_cast<Real>(shell), -cfg.expo);
      for (int idx = 0; idx < static_cast<int>(members.size()); ++idx) {
        const int n = members[idx];
        const Real mode_target = shell_target * weights[idx] / weight_sum;
        const Real scale = std::sqrt(mode_target);
        coeffs.aka0[n] *= scale;
        coeffs.aka1[n] *= scale;
        coeffs.aka2[n] *= scale;
        coeffs.akb0[n] *= scale;
        coeffs.akb1[n] *= scale;
        coeffs.akb2[n] *= scale;
      }
    }

    return coeffs;
  }

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

  if (cfg.exact_shell_contract) {
    std::vector<std::vector<int>> shell_members(cfg.nhigh + 1);
    for (int n = 0; n < nmodes; ++n) {
      shell_members[catalog.shell[n]].push_back(n);
    }

    for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
      const auto &members = shell_members[shell];
      if (members.empty()) continue;

      std::vector<Real> weights(members.size(), 0.0);
      Real weight_sum = 0.0;
      for (int idx = 0; idx < static_cast<int>(members.size()); ++idx) {
        const int n = members[idx];
        const Real weight = 1.0 / std::pow(catalog.kiso[n], cfg.expo + 1.0);
        weights[idx] = weight;
        weight_sum += weight;
      }

      const Real shell_target = std::pow(static_cast<Real>(shell), -cfg.expo);
      for (int idx = 0; idx < static_cast<int>(members.size()); ++idx) {
        const int n = members[idx];
        const Real mode_target = shell_target * weights[idx] / weight_sum;
        const Real amplitude = std::sqrt(2.0 * mode_target) / catalog.kiso[n];
        const Real phase = 2.0 * kPi * RanSt(rstate);
        coeffs.aka[n] = amplitude * std::cos(phase);
        coeffs.akb[n] = amplitude * std::sin(phase);
      }
    }

    return coeffs;
  }

  for (int n = 0; n < nmodes; ++n) {
    const Real norm = Stream2DPsiAmplitudeNorm(cfg, catalog.kiso[n]) * catalog.boost[n];
    coeffs.aka[n] = norm * RanGaussianSt(rstate);
    coeffs.akb[n] = norm * RanGaussianSt(rstate);
  }

  return coeffs;
}

ScalarCoefficients GenerateStream3DPhiCoefficients(const ModeCatalog &catalog,
                                                   RNG_State *rstate,
                                                   Real beta) {
  ScalarCoefficients coeffs;
  const int nmodes = catalog.KeptModes();
  coeffs.aka.resize(nmodes);
  coeffs.akb.resize(nmodes);

  for (int n = 0; n < nmodes; ++n) {
    const Real kiso = catalog.kiso[n];
    const Real norm = (kiso > kTiny)
                          ? catalog.boost[n] / std::pow(kiso, beta)
                          : 0.0;
    coeffs.aka[n] = norm * RanGaussianSt(rstate);
    coeffs.akb[n] = norm * RanGaussianSt(rstate);
  }

  return coeffs;
}

std::vector<SignedMode> BuildSignedModes(const ModeCatalog &catalog) {
  std::vector<SignedMode> signed_modes;
  signed_modes.reserve(2 * catalog.KeptModes());

  for (int n = 0; n < catalog.KeptModes(); ++n) {
    const Real variance_weight = 0.5 * catalog.boost[n] * catalog.boost[n];
    signed_modes.push_back({catalog.nkx[n], catalog.nky[n], catalog.nkz[n],
                            catalog.kx[n], catalog.ky[n], catalog.kz[n],
                            catalog.kiso[n], variance_weight});
    signed_modes.push_back({-catalog.nkx[n], -catalog.nky[n], -catalog.nkz[n],
                            -catalog.kx[n], -catalog.ky[n], -catalog.kz[n],
                            catalog.kiso[n], variance_weight});
  }

  return signed_modes;
}

ClebschPairSamples SampleClebschPairResponses(const std::vector<SignedMode> &signed1,
                                              const std::vector<SignedMode> &signed2,
                                              int nhigh, int nsamples,
                                              RNG_State *rstate) {
  ClebschPairSamples samples;
  samples.idx1.reserve(nsamples);
  samples.idx2.reserve(nsamples);
  samples.shell.reserve(nsamples);
  samples.geom.reserve(nsamples);

  if (signed1.empty() || signed2.empty() || nsamples <= 0) return samples;

  const int n1 = static_cast<int>(signed1.size());
  const int n2 = static_cast<int>(signed2.size());
  const int max_shell = 2 * nhigh;

  while (samples.Count() < nsamples) {
    const int i1 = std::min(static_cast<int>(RanSt(rstate) * n1), n1 - 1);
    const int i2 = std::min(static_cast<int>(RanSt(rstate) * n2), n2 - 1);
    const auto &mode_p = signed1[i1];
    const auto &mode_q = signed2[i2];

    const int rx = mode_p.nkx + mode_q.nkx;
    const int ry = mode_p.nky + mode_q.nky;
    const int rz = mode_p.nkz + mode_q.nkz;
    if (rx == 0 && ry == 0 && rz == 0) continue;

    const int out_shell = ShellFromIndexSqr(rx*rx + ry*ry + rz*rz);
    if (out_shell < 1 || out_shell > max_shell) continue;

    const Real cx = mode_p.ky*mode_q.kz - mode_p.kz*mode_q.ky;
    const Real cy = mode_p.kz*mode_q.kx - mode_p.kx*mode_q.kz;
    const Real cz = mode_p.kx*mode_q.ky - mode_p.ky*mode_q.kx;
    const Real geom = cx*cx + cy*cy + cz*cz;
    if (geom <= kTiny) continue;

    samples.idx1.push_back(i1);
    samples.idx2.push_back(i2);
    samples.shell.push_back(out_shell);
    samples.geom.push_back(geom);
  }

  return samples;
}

std::vector<Real> EstimateClebschShellSpectrum(const std::vector<SignedMode> &signed1,
                                               const std::vector<SignedMode> &signed2,
                                               const ClebschPairSamples &samples,
                                               int nhigh, Real beta) {
  std::vector<Real> power1(signed1.size(), 0.0);
  std::vector<Real> power2(signed2.size(), 0.0);
  const Real two_beta = 2.0 * beta;

  for (int n = 0; n < static_cast<int>(signed1.size()); ++n) {
    power1[n] = signed1[n].variance_weight / std::pow(signed1[n].kiso, two_beta);
  }
  for (int n = 0; n < static_cast<int>(signed2.size()); ++n) {
    power2[n] = signed2[n].variance_weight / std::pow(signed2[n].kiso, two_beta);
  }

  std::vector<Real> spectrum(2 * nhigh + 1, 0.0);
  for (int sample = 0; sample < samples.Count(); ++sample) {
    spectrum[samples.shell[sample]] +=
        samples.geom[sample] * power1[samples.idx1[sample]] * power2[samples.idx2[sample]];
  }

  return spectrum;
}

Stream3DBetaCalibration CalibrateStream3DBeta(const TurbulenceConfig &cfg,
                                              const ModeCatalog &catalog1,
                                              const ModeCatalog &catalog2) {
  Stream3DBetaCalibration result;
  result.asymptotic_beta = 0.25 * (cfg.expo + 8.0);
  result.calibrated_beta = result.asymptotic_beta;

  const auto signed1 = BuildSignedModes(catalog1);
  const auto signed2 = BuildSignedModes(catalog2);
  if (signed1.empty() || signed2.empty()) return result;

  RNG_State sample_state;
  sample_state.idum = -MixSeed(cfg.rseed, 91);
  const int requested_samples = std::max(262144, 4096 * cfg.nhigh);
  const ClebschPairSamples samples =
      SampleClebschPairResponses(signed1, signed2, cfg.nhigh, requested_samples, &sample_state);
  result.sample_count = samples.Count();
  if (samples.Count() == 0) return result;

  const Real target_slope = -cfg.expo;
  Real best_beta = result.asymptotic_beta;
  Real best_slope = 0.0;
  Real best_err = 1.0e99;

  auto evaluate_beta = [&](Real beta) {
    const std::vector<Real> spectrum =
        EstimateClebschShellSpectrum(signed1, signed2, samples, cfg.nhigh, beta);
    const Real slope = FitLogSlope(spectrum, cfg.nlow, cfg.nhigh);
    const Real err = std::abs(slope - target_slope);
    if (err < best_err) {
      best_err = err;
      best_beta = beta;
      best_slope = slope;
    }
    return slope;
  };

  result.asymptotic_slope = evaluate_beta(result.asymptotic_beta);
  result.calibrated_slope = result.asymptotic_slope;

  if (best_err <= 0.05) {
    result.calibrated_beta = best_beta;
    result.calibrated_slope = best_slope;
    return result;
  }

  Real beta_lo = std::max(0.25, result.asymptotic_beta - 1.0);
  Real beta_hi = result.asymptotic_beta + 1.0;
  Real slope_lo = evaluate_beta(beta_lo);
  Real slope_hi = evaluate_beta(beta_hi);

  auto is_bracketed = [&](Real slo, Real shi) {
    return ((slo - target_slope) * (shi - target_slope) <= 0.0);
  };

  int expand_iter = 0;
  while (!is_bracketed(slope_lo, slope_hi) && expand_iter < 8) {
    if (slope_lo > target_slope && slope_hi > target_slope) {
      beta_lo = beta_hi;
      slope_lo = slope_hi;
      beta_hi += 1.0;
      slope_hi = evaluate_beta(beta_hi);
    } else if (slope_lo < target_slope && slope_hi < target_slope) {
      if (beta_lo <= 0.25) break;
      beta_hi = beta_lo;
      slope_hi = slope_lo;
      beta_lo = std::max(0.25, beta_lo - 1.0);
      slope_lo = evaluate_beta(beta_lo);
    } else {
      break;
    }
    ++expand_iter;
  }

  result.bracketed = is_bracketed(slope_lo, slope_hi);
  if (result.bracketed) {
    for (int iter = 0; iter < 12; ++iter) {
      const Real beta_mid = 0.5 * (beta_lo + beta_hi);
      const Real slope_mid = evaluate_beta(beta_mid);
      if ((slope_lo - target_slope) * (slope_mid - target_slope) <= 0.0) {
        beta_hi = beta_mid;
        slope_hi = slope_mid;
      } else {
        beta_lo = beta_mid;
        slope_lo = slope_mid;
      }
    }
  }

  result.calibrated_beta = best_beta;
  result.calibrated_slope = best_slope;
  return result;
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
  MPI_Allreduce(local_mean, global_mean, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
  MPI_Allreduce(local_rms, global_rms, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

void PopulateScalarFaceVelocitiesFromClebsch(
    MeshBlockPack *pmbp, const TurbulenceConfig &cfg,
    const ModeCatalog &catalog1, const ScalarCoefficients &phi1_coeffs,
    const ModeCatalog &catalog2, const ScalarCoefficients &phi2_coeffs,
    const VelocityStats &stats) {
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

  const int nmodes1 = catalog1.KeptModes();
  const int nmodes2 = catalog2.KeptModes();
  if (nmodes1 == 0 || nmodes2 == 0) return;

  auto &hydro = *pmbp->phydro;
  const int ncells1 = nx1 + 2*ng;
  const int ncells2 = multi_d ? nx2 + 2*ng : 1;
  const int ncells3 = three_d ? nx3 + 2*ng : 1;
  if (!hydro.scalar_vface) {
    hydro.scalar_vface = std::make_unique<DvceFaceFld4D<Real>>(
        "scalar_vface", nmb, ncells3, ncells2, ncells1);
  }
  hydro.use_scalar_face_velocity = true;

  DvceArray1D<Real> d_kx1("scalar_face_stream_kx1", nmodes1);
  DvceArray1D<Real> d_ky1("scalar_face_stream_ky1", nmodes1);
  DvceArray1D<Real> d_kz1("scalar_face_stream_kz1", nmodes1);
  DvceArray1D<Real> d_a1("scalar_face_stream_a1", nmodes1);
  DvceArray1D<Real> d_b1("scalar_face_stream_b1", nmodes1);
  DvceArray1D<Real> d_kx2("scalar_face_stream_kx2", nmodes2);
  DvceArray1D<Real> d_ky2("scalar_face_stream_ky2", nmodes2);
  DvceArray1D<Real> d_kz2("scalar_face_stream_kz2", nmodes2);
  DvceArray1D<Real> d_a2("scalar_face_stream_a2", nmodes2);
  DvceArray1D<Real> d_b2("scalar_face_stream_b2", nmodes2);

  auto h_kx1 = Kokkos::create_mirror_view(d_kx1);
  auto h_ky1 = Kokkos::create_mirror_view(d_ky1);
  auto h_kz1 = Kokkos::create_mirror_view(d_kz1);
  auto h_a1 = Kokkos::create_mirror_view(d_a1);
  auto h_b1 = Kokkos::create_mirror_view(d_b1);
  auto h_kx2 = Kokkos::create_mirror_view(d_kx2);
  auto h_ky2 = Kokkos::create_mirror_view(d_ky2);
  auto h_kz2 = Kokkos::create_mirror_view(d_kz2);
  auto h_a2 = Kokkos::create_mirror_view(d_a2);
  auto h_b2 = Kokkos::create_mirror_view(d_b2);

  for (int n = 0; n < nmodes1; ++n) {
    h_kx1(n) = catalog1.kx[n];
    h_ky1(n) = catalog1.ky[n];
    h_kz1(n) = catalog1.kz[n];
    h_a1(n) = phi1_coeffs.aka[n];
    h_b1(n) = phi1_coeffs.akb[n];
  }
  for (int n = 0; n < nmodes2; ++n) {
    h_kx2(n) = catalog2.kx[n];
    h_ky2(n) = catalog2.ky[n];
    h_kz2(n) = catalog2.kz[n];
    h_a2(n) = phi2_coeffs.aka[n];
    h_b2(n) = phi2_coeffs.akb[n];
  }

  Kokkos::deep_copy(d_kx1, h_kx1);
  Kokkos::deep_copy(d_ky1, h_ky1);
  Kokkos::deep_copy(d_kz1, h_kz1);
  Kokkos::deep_copy(d_a1, h_a1);
  Kokkos::deep_copy(d_b1, h_b1);
  Kokkos::deep_copy(d_kx2, h_kx2);
  Kokkos::deep_copy(d_ky2, h_ky2);
  Kokkos::deep_copy(d_kz2, h_kz2);
  Kokkos::deep_copy(d_a2, h_a2);
  Kokkos::deep_copy(d_b2, h_b2);

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
  par_for("scalar_face_stream_x1", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
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

    Real g1x = 0.0, g1y = 0.0, g1z = 0.0;
    Real g2x = 0.0, g2y = 0.0, g2z = 0.0;
    for (int n = 0; n < nmodes1; ++n) {
      const Real phase = d_kx1(n)*x1v + d_ky1(n)*x2v + d_kz1(n)*x3v;
      const Real basis = d_a1(n)*sin(phase) + d_b1(n)*cos(phase);
      g1x -= d_kx1(n) * basis;
      g1y -= d_ky1(n) * basis;
      g1z -= d_kz1(n) * basis;
    }
    for (int n = 0; n < nmodes2; ++n) {
      const Real phase = d_kx2(n)*x1v + d_ky2(n)*x2v + d_kz2(n)*x3v;
      const Real basis = d_a2(n)*sin(phase) + d_b2(n)*cos(phase);
      g2x -= d_kx2(n) * basis;
      g2y -= d_ky2(n) * basis;
      g2z -= d_kz2(n) * basis;
    }

    const Real vx = g1y*g2z - g1z*g2y;
    sface_x1f(m,k,j,i) = (vx - stats.vmean1) * stats.scale;
  });

  if (multi_d) {
    il = is - ng;
    iu = ie + ng;
    jl = js - ng;
    ju = je + ng + 1;
    kl = three_d ? ks - ng : ks;
    ku = three_d ? ke + ng : ke;
    par_for("scalar_face_stream_x2", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
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

      Real g1x = 0.0, g1y = 0.0, g1z = 0.0;
      Real g2x = 0.0, g2y = 0.0, g2z = 0.0;
      for (int n = 0; n < nmodes1; ++n) {
        const Real phase = d_kx1(n)*x1v + d_ky1(n)*x2v + d_kz1(n)*x3v;
        const Real basis = d_a1(n)*sin(phase) + d_b1(n)*cos(phase);
        g1x -= d_kx1(n) * basis;
        g1y -= d_ky1(n) * basis;
        g1z -= d_kz1(n) * basis;
      }
      for (int n = 0; n < nmodes2; ++n) {
        const Real phase = d_kx2(n)*x1v + d_ky2(n)*x2v + d_kz2(n)*x3v;
        const Real basis = d_a2(n)*sin(phase) + d_b2(n)*cos(phase);
        g2x -= d_kx2(n) * basis;
        g2y -= d_ky2(n) * basis;
        g2z -= d_kz2(n) * basis;
      }

      const Real vy = g1z*g2x - g1x*g2z;
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
    par_for("scalar_face_stream_x3", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      const Real x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);

      Real g1x = 0.0, g1y = 0.0, g1z = 0.0;
      Real g2x = 0.0, g2y = 0.0, g2z = 0.0;
      for (int n = 0; n < nmodes1; ++n) {
        const Real phase = d_kx1(n)*x1v + d_ky1(n)*x2v + d_kz1(n)*x3v;
        const Real basis = d_a1(n)*sin(phase) + d_b1(n)*cos(phase);
        g1x -= d_kx1(n) * basis;
        g1y -= d_ky1(n) * basis;
        g1z -= d_kz1(n) * basis;
      }
      for (int n = 0; n < nmodes2; ++n) {
        const Real phase = d_kx2(n)*x1v + d_ky2(n)*x2v + d_kz2(n)*x3v;
        const Real basis = d_a2(n)*sin(phase) + d_b2(n)*cos(phase);
        g2x -= d_kx2(n) * basis;
        g2y -= d_ky2(n) * basis;
        g2z -= d_kz2(n) * basis;
      }

      const Real vz = g1x*g2y - g1y*g2x;
      sface_x3f(m,k,j,i) = (vz - stats.vmean3) * stats.scale;
    });
  }
}

void SynthesizeClebschVelocity(MeshBlockPack *pmbp, const TurbulenceConfig &cfg,
                               Real den, const ModeCatalog &catalog1,
                               const ScalarCoefficients &phi1_coeffs,
                               const ModeCatalog &catalog2,
                               const ScalarCoefficients &phi2_coeffs) {
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

  const int nmodes1 = catalog1.KeptModes();
  const int nmodes2 = catalog2.KeptModes();
  if (nmodes1 == 0 || nmodes2 == 0) return;

  DvceArray1D<Real> d_kx1("scalar_mix_stream_kx1", nmodes1);
  DvceArray1D<Real> d_ky1("scalar_mix_stream_ky1", nmodes1);
  DvceArray1D<Real> d_kz1("scalar_mix_stream_kz1", nmodes1);
  DvceArray1D<Real> d_a1("scalar_mix_stream_a1", nmodes1);
  DvceArray1D<Real> d_b1("scalar_mix_stream_b1", nmodes1);
  DvceArray1D<Real> d_kx2("scalar_mix_stream_kx2", nmodes2);
  DvceArray1D<Real> d_ky2("scalar_mix_stream_ky2", nmodes2);
  DvceArray1D<Real> d_kz2("scalar_mix_stream_kz2", nmodes2);
  DvceArray1D<Real> d_a2("scalar_mix_stream_a2", nmodes2);
  DvceArray1D<Real> d_b2("scalar_mix_stream_b2", nmodes2);

  auto h_kx1 = Kokkos::create_mirror_view(d_kx1);
  auto h_ky1 = Kokkos::create_mirror_view(d_ky1);
  auto h_kz1 = Kokkos::create_mirror_view(d_kz1);
  auto h_a1 = Kokkos::create_mirror_view(d_a1);
  auto h_b1 = Kokkos::create_mirror_view(d_b1);
  auto h_kx2 = Kokkos::create_mirror_view(d_kx2);
  auto h_ky2 = Kokkos::create_mirror_view(d_ky2);
  auto h_kz2 = Kokkos::create_mirror_view(d_kz2);
  auto h_a2 = Kokkos::create_mirror_view(d_a2);
  auto h_b2 = Kokkos::create_mirror_view(d_b2);

  for (int n = 0; n < nmodes1; ++n) {
    h_kx1(n) = catalog1.kx[n];
    h_ky1(n) = catalog1.ky[n];
    h_kz1(n) = catalog1.kz[n];
    h_a1(n) = phi1_coeffs.aka[n];
    h_b1(n) = phi1_coeffs.akb[n];
  }
  for (int n = 0; n < nmodes2; ++n) {
    h_kx2(n) = catalog2.kx[n];
    h_ky2(n) = catalog2.ky[n];
    h_kz2(n) = catalog2.kz[n];
    h_a2(n) = phi2_coeffs.aka[n];
    h_b2(n) = phi2_coeffs.akb[n];
  }

  Kokkos::deep_copy(d_kx1, h_kx1);
  Kokkos::deep_copy(d_ky1, h_ky1);
  Kokkos::deep_copy(d_kz1, h_kz1);
  Kokkos::deep_copy(d_a1, h_a1);
  Kokkos::deep_copy(d_b1, h_b1);
  Kokkos::deep_copy(d_kx2, h_kx2);
  Kokkos::deep_copy(d_ky2, h_ky2);
  Kokkos::deep_copy(d_kz2, h_kz2);
  Kokkos::deep_copy(d_a2, h_a2);
  Kokkos::deep_copy(d_b2, h_b2);

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;

  par_for("scalar_mix_stream3d_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
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

    Real g1x = 0.0, g1y = 0.0, g1z = 0.0;
    Real g2x = 0.0, g2y = 0.0, g2z = 0.0;

    for (int n = 0; n < nmodes1; ++n) {
      const Real phase = d_kx1(n)*x1v + d_ky1(n)*x2v + d_kz1(n)*x3v;
      const Real basis = d_a1(n)*sin(phase) + d_b1(n)*cos(phase);
      g1x -= d_kx1(n) * basis;
      g1y -= d_ky1(n) * basis;
      g1z -= d_kz1(n) * basis;
    }
    for (int n = 0; n < nmodes2; ++n) {
      const Real phase = d_kx2(n)*x1v + d_ky2(n)*x2v + d_kz2(n)*x3v;
      const Real basis = d_a2(n)*sin(phase) + d_b2(n)*cos(phase);
      g2x -= d_kx2(n) * basis;
      g2y -= d_ky2(n) * basis;
      g2z -= d_kz2(n) * basis;
    }

    const Real vx = g1y*g2z - g1z*g2y;
    const Real vy = g1z*g2x - g1x*g2z;
    const Real vz = g1x*g2y - g1y*g2x;

    u0(m,IM1,k,j,i) = den * vx;
    u0(m,IM2,k,j,i) = den * vy;
    u0(m,IM3,k,j,i) = den * vz;
  });
}

void LogModeSummary(const TurbulenceConfig &cfg, const ModeCatalog &catalog,
                    const std::string &label) {
  std::cout << "  " << label << ": total_modes_in_shell=" << catalog.total_modes;
  if (cfg.exact_shell_contract) {
    std::cout << ", kept_via_shell_quota=" << catalog.KeptModes();
  } else {
    std::cout << ", kept_via_importance_sampling=" << catalog.KeptModes();
  }
  std::cout << std::endl;
  if (cfg.method == VelocityMethod::StreamFunction && cfg.stream_2d) {
    std::cout << "  " << label << ": stream-function coefficients preserve "
              << "E_u(k) ~ k^(-expo) through |psi_k| ~ k^{-(expo+3)/2}" << std::endl;
  }
  if (cfg.exact_shell_contract) {
    std::cout << "  " << label << ": exact_shell enforces raw centered-shell energy "
              << "E(k_shell) ~ k_shell^(-expo)" << std::endl;
  }
}

void InitTurbulentVelocity(MeshBlockPack *pmbp, ParameterInput *pin, Real den,
                           bool projection_true_2d) {
  const TurbulenceConfig cfg = ReadTurbulenceConfig(pin, pmbp->pmesh, projection_true_2d);
  ValidateExactShellContract(cfg);
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
      return;
    }

    if (global_variable::my_rank == 0) {
      std::cout << "Initializing turbulent velocity field:" << std::endl
                << "  velocity_method=projection" << std::endl
                << "  v_rms=" << cfg.v_rms << ", nlow=" << cfg.nlow
                << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo
                << ", sol_frac=" << cfg.sol_frac << std::endl
                << "  spectrum_contract=" << cfg.spectrum_contract << std::endl
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
    return;
  }

  if (global_variable::my_rank == 0) {
    std::cout << "Initializing turbulent velocity field:" << std::endl
              << "  velocity_method=stream_function" << std::endl
              << "  v_rms=" << cfg.v_rms << ", nlow=" << cfg.nlow
              << ", nhigh=" << cfg.nhigh << ", expo=" << cfg.expo
              << ", sol_frac=" << cfg.sol_frac << std::endl
              << "  spectrum_contract=" << cfg.spectrum_contract << std::endl
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
    return;
  }

  if (cfg.exact_shell_contract) {
    FatalProblemSetup("exact_shell does not support 3D stream-function initialization.");
  }

  RNG_State rstate1;
  RNG_State rstate2;
  rstate1.idum = -MixSeed(cfg.rseed, 1);
  rstate2.idum = -MixSeed(cfg.rseed, 2);

  const ModeCatalog catalog1 = BuildModeCatalog(cfg, &rstate1);
  const ModeCatalog catalog2 = BuildModeCatalog(cfg, &rstate2);
  if (catalog1.KeptModes() == 0 || catalog2.KeptModes() == 0) {
    std::cout << "### WARNING: No turbulent modes kept in range [" << cfg.nlow << ", "
              << cfg.nhigh << "] for 3D stream-function initialization. Using zero velocity."
              << std::endl;
    return;
  }

  const Stream3DBetaCalibration beta_fit = CalibrateStream3DBeta(cfg, catalog1, catalog2);
  if (global_variable::my_rank == 0) {
    std::cout << "  stream_geometry=3D direct-sampled Clebsch "
                 "(v = grad(phi1) x grad(phi2))"
              << std::endl;
    LogModeSummary(cfg, catalog1, "phi_1");
    LogModeSummary(cfg, catalog2, "phi_2");
    std::cout << "  phi_amplitude_norm = k^(-beta), asymptotic_beta="
              << beta_fit.asymptotic_beta
              << ", calibrated_beta=" << beta_fit.calibrated_beta
              << ", estimated_in_band_slope(asymptotic)=" << beta_fit.asymptotic_slope
              << ", estimated_in_band_slope(calibrated)=" << beta_fit.calibrated_slope
              << ", response_samples=" << beta_fit.sample_count << std::endl;
  }

  const ScalarCoefficients phi1_coeffs =
      GenerateStream3DPhiCoefficients(catalog1, &rstate1, beta_fit.calibrated_beta);
  const ScalarCoefficients phi2_coeffs =
      GenerateStream3DPhiCoefficients(catalog2, &rstate2, beta_fit.calibrated_beta);
  SynthesizeClebschVelocity(pmbp, cfg, den, catalog1, phi1_coeffs, catalog2, phi2_coeffs);
  const VelocityStats stats = NormalizeVelocityField(pmbp, cfg);
  if (face_requested && face_scalar_ok) {
    if (global_variable::my_rank == 0) {
      std::cout << "  using divergence-free face velocity for scalar fluxes" << std::endl;
    }
    PopulateScalarFaceVelocitiesFromClebsch(
        pmbp, cfg, catalog1, phi1_coeffs, catalog2, phi2_coeffs, stats);
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar mixing
//
//  Initializes density, turbulent velocity, and passive scalars.
//  Scalar 0 uses the existing uniform/step initialization path, while scalar 1
//  is initialized to a smooth sinusoidal product if it is present.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  const int nscalars = pmbp->phydro->nscalars;
  const int nhydro = pmbp->phydro->nhydro;
  const bool true_2d = UseTrue2DVelocity(pin, pmy_mesh_);
  const bool multi_d = pmy_mesh_->multi_d;
  const bool three_d = pmy_mesh_->three_d;
  const bool has_sine_scalar = (nscalars > 1);

  if (nscalars == 0) {
    FatalProblemSetup("Scalar mixing requires nscalars > 0.");
  }

  const Real legacy_G = pin->GetOrAddReal("problem", "mean_gradient", 1.0);
  const Real scalar_init_default = pin->GetOrAddReal("problem", "scalar_init", 0.5);
  const bool scalar0_use_sine_x1 =
      pin->GetOrAddBoolean("problem", "scalar0_use_sine_x1", true);
  const Real scalar0_wavelengths_per_box =
      pin->GetOrAddReal("problem", "scalar0_wavelengths_per_box",
                        kDefaultScalar0WavelengthsPerBox);
  if (scalar0_wavelengths_per_box <= 0.0) {
    FatalProblemSetup("problem/scalar0_wavelengths_per_box must be > 0.");
  }
  const bool scalar_use_x1_step =
      pin->GetOrAddBoolean("problem", "scalar_use_x1_step", false);
  const Real scalar_step_x1 = pin->GetOrAddReal("problem", "scalar_step_x1", 0.0);
  const Real scalar_left = pin->GetOrAddReal("problem", "scalar_left", 0.0);
  const Real scalar_right = pin->GetOrAddReal("problem", "scalar_right", 1.0);
  const Real scalar1_wavelengths_per_box =
      pin->GetOrAddReal("problem", "scalar1_wavelengths_per_box",
                        kDefaultScalar1WavelengthsPerBox);
  if (scalar1_wavelengths_per_box <= 0.0) {
    FatalProblemSetup("problem/scalar1_wavelengths_per_box must be > 0.");
  }
  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  const Real x2min = pmy_mesh_->mesh_size.x2min;
  const Real x2max = pmy_mesh_->mesh_size.x2max;
  const Real x3min = pmy_mesh_->mesh_size.x3min;
  const Real x3max = pmy_mesh_->mesh_size.x3max;
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
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  const Real gm1 = eos.gamma - 1.0;

  // Read initial conditions
  const Real den = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real temp = pin->GetOrAddReal("problem", "temp0", 1.0);
  const int nmb = pmbp->nmb_thispack;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar mixing problem:" << std::endl
              << "  scalar0_ic = "
              << (scalar0_use_sine_x1 ? "sine_x1_global"
                                      : (scalar_use_x1_step ? "x1_step" : "uniform"))
              << std::endl
              << "  true_2d = " << (true_2d ? "true" : "false") << std::endl
              << "  nscalars = " << nscalars << std::endl
              << "  scalar_init (default) = " << scalar_init_default << std::endl
              << "  scalar0_wavelengths_per_box = " << scalar0_wavelengths_per_box
              << std::endl
              << "  scalar1_ic = "
              << (has_sine_scalar
                      ? (three_d ? "sine_product_3d"
                                 : (multi_d ? "sine_product_2d"
                                            : "sine_x1"))
                      : "disabled (set nscalars >= 2)")
              << std::endl
              << "  scalar1_wavelengths_per_box = " << scalar1_wavelengths_per_box
              << std::endl
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
    const bool sine_scalar = (ns == 1);
    const std::string label = "scalar_mix_pgen_scalar_" + std::to_string(ns);
    par_for(label, DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real rho = u0(m, IDN, k, j, i);
      Real theta = init;
      if (ns == 0 && scalar0_use_sine_x1) {
        Real &mb_x1min = size.d_view(m).x1min;
        Real &mb_x1max = size.d_view(m).x1max;
        const Real x1v = CellCenterX(i-is, nx1, mb_x1min, mb_x1max);
        theta = SineScalarIC1D(x1v, x1min, x1max, scalar0_wavelengths_per_box);
      } else if (sine_scalar) {
        Real &mb_x1min = size.d_view(m).x1min;
        Real &mb_x1max = size.d_view(m).x1max;
        const Real x1v = CellCenterX(i-is, nx1, mb_x1min, mb_x1max);
        Real x2v = 0.0;
        Real x3v = 0.0;
        if (multi_d) {
          Real &mb_x2min = size.d_view(m).x2min;
          Real &mb_x2max = size.d_view(m).x2max;
          x2v = CellCenterX(j-js, nx2, mb_x2min, mb_x2max);
        }
        if (three_d) {
          Real &mb_x3min = size.d_view(m).x3min;
          Real &mb_x3max = size.d_view(m).x3max;
          x3v = CellCenterX(k-ks, nx3, mb_x3min, mb_x3max);
        }
        theta = SineScalarIC(x1v, x2v, x3v,
                             x1min, x1max,
                             x2min, x2max,
                             x3min, x3max,
                             scalar1_wavelengths_per_box,
                             multi_d, three_d);
      } else if (scalar_use_x1_step) {
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
      MPI_Allreduce(local_sum, global_sum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
