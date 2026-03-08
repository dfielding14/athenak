//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_mixing.cpp
//  \brief Problem generator for scalar mixing in frozen turbulent velocity field
//
//  Single passive scalar with mean gradient forcing: ds/dt = G * v_x
//  This simulates mixing of a background scalar gradient by turbulence.
//
//  When scalar_only=true in <hydro>, velocity is frozen and only the scalar evolves.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
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
void MeanGradientForcing(Mesh* pm, const Real bdt);

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

  bool mesh_multi_d = false;
  bool mesh_three_d = false;
  bool projection_true_2d = false;
  bool stream_2d = false;
  bool active_v3 = false;

  Real lx = 1.0;
  Real ly = 1.0;
  Real lz = 1.0;
  Real dkx = 0.0;
  Real dky = 0.0;
  Real dkz = 0.0;
  Real k_crit_mag = 0.0;
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

struct ShellFitEval {
  std::vector<Real> pred;
  std::vector<Real> grad;
  Real loss = 0.0;
};

struct ShellFitResult {
  std::vector<Real> shell_amp;
  std::vector<Real> predicted_shell_energy;
  Real initial_loss = 0.0;
  Real final_loss = 0.0;
  Real leakage_fraction = 0.0;
  Real fitted_slope = 0.0;
  bool converged = false;
};

// Store mean gradient parameter (read in ProblemGenerator, used in source term)
Real mean_gradient_G = 1.0;

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

  return cfg;
}

Real ModeAcceptanceProbability(const TurbulenceConfig &cfg, Real kiso) {
  if (kiso < cfg.k_crit_mag) return 1.0;
  return (cfg.k_crit_mag * cfg.k_crit_mag) / (kiso * kiso);
}

ModeCatalog BuildModeCatalog(const TurbulenceConfig &cfg, RNG_State *rstate) {
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

void NormalizeVelocityField(MeshBlockPack *pmbp, const TurbulenceConfig &cfg) {
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
}

std::vector<int> CountModesPerShell(const ModeCatalog &catalog, int nhigh) {
  std::vector<int> counts(nhigh + 1, 0);
  for (int shell : catalog.shell) {
    if (shell >= 1 && shell <= nhigh) ++counts[shell];
  }
  return counts;
}

void EnsureShellCoverage(const std::vector<int> &counts, const TurbulenceConfig &cfg,
                         const std::string &label) {
  for (int shell = cfg.nlow; shell <= cfg.nhigh; ++shell) {
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
    signed_modes.push_back({-catalog.nkx[n], -catalog.nky[n], -catalog.nkz[n], catalog.shell[n],
                            -catalog.kx[n], -catalog.ky[n], -catalog.kz[n], variance_weight});
  }

  return signed_modes;
}

int TensorIndex(int kout, int s, int t, int nhigh) {
  const int ns = nhigh + 1;
  return (kout * ns + s) * ns + t;
}

std::vector<Real> BuildClebschShellResponseTensor(const TurbulenceConfig &cfg,
                                                  const ModeCatalog &catalog1,
                                                  const ModeCatalog &catalog2) {
  const int nhigh = cfg.nhigh;
  const int nkout = 2 * nhigh;
  std::vector<Real> tensor((nkout + 1) * (nhigh + 1) * (nhigh + 1), 0.0);

  // For phi_j(x) = Re[sum_p c_j(p) exp(i p·x)], the Clebsch field satisfies
  // u_hat(r) = -sum_{p+q=r} (p x q) c_1(p) c_2(q). Averaging over independent
  // Gaussian coefficients gives the shell response
  // E_u(K) = sum_{s,t} T[K,s,t] P_s P_t, where P_s is the scalar-power weight
  // in shell s. Local-triad counting then implies E_u(k) ~ k^(8-4 beta) when
  // |phi_k| ~ k^(-beta), which motivates beta ~= (expo + 8)/4.
  const auto signed1 = BuildSignedModes(catalog1);
  const auto signed2 = BuildSignedModes(catalog2);

  for (const auto &mode_p : signed1) {
    for (const auto &mode_q : signed2) {
      const int rx = mode_p.nkx + mode_q.nkx;
      const int ry = mode_p.nky + mode_q.nky;
      const int rz = mode_p.nkz + mode_q.nkz;
      if (rx == 0 && ry == 0 && rz == 0) continue;

      const int out_shell = ShellFromIndexSqr(rx*rx + ry*ry + rz*rz);
      if (out_shell < 1 || out_shell > nkout) continue;

      const Real cx = mode_p.ky*mode_q.kz - mode_p.kz*mode_q.ky;
      const Real cy = mode_p.kz*mode_q.kx - mode_p.kx*mode_q.kz;
      const Real cz = mode_p.kx*mode_q.ky - mode_p.ky*mode_q.kx;
      const Real geom = cx*cx + cy*cy + cz*cz;
      if (geom <= kTiny) continue;

      tensor[TensorIndex(out_shell, mode_p.shell, mode_q.shell, nhigh)] +=
          mode_p.weight * mode_q.weight * geom;
    }
  }

  return tensor;
}

void SymmetrizeClebschTensor(std::vector<Real> *tensor, int nhigh) {
  const int nkout = 2 * nhigh;
  for (int kout = 1; kout <= nkout; ++kout) {
    for (int s = 1; s <= nhigh; ++s) {
      for (int t = s + 1; t <= nhigh; ++t) {
        const int idx_st = TensorIndex(kout, s, t, nhigh);
        const int idx_ts = TensorIndex(kout, t, s, nhigh);
        const Real avg = 0.5 * ((*tensor)[idx_st] + (*tensor)[idx_ts]);
        (*tensor)[idx_st] = avg;
        (*tensor)[idx_ts] = avg;
      }
    }
  }
}

void EnsureReachableTargetBand(const std::vector<Real> &tensor, const TurbulenceConfig &cfg) {
  for (int kout = cfg.nlow; kout <= cfg.nhigh; ++kout) {
    Real shell_sum = 0.0;
    for (int s = 1; s <= cfg.nhigh; ++s) {
      for (int t = 1; t <= cfg.nhigh; ++t) {
        shell_sum += tensor[TensorIndex(kout, s, t, cfg.nhigh)];
      }
    }
    if (shell_sum <= kTiny) {
      FatalProblemSetup("3D stream-function shell response has zero support in target shell "
                        + std::to_string(kout) + ".");
    }
  }
}

ShellFitEval EvaluateClebschShellFit(const std::vector<Real> &tensor,
                                     const TurbulenceConfig &cfg,
                                     const std::vector<Real> &shell_power) {
  const int nhigh = cfg.nhigh;
  const int nkout = 2 * nhigh;
  const int ns = nhigh + 1;
  const Real leak_ref = TargetShellEnergy(std::max(cfg.nlow, 1), cfg.expo);
  const Real smooth_weight = 1.0e-2;

  ShellFitEval eval;
  eval.pred.assign(nkout + 1, 0.0);
  eval.grad.assign(nhigh + 1, 0.0);

  std::vector<Real> kernel((nkout + 1) * ns, 0.0);
  for (int kout = 1; kout <= nkout; ++kout) {
    for (int s = 1; s <= nhigh; ++s) {
      Real accum = 0.0;
      for (int t = 1; t <= nhigh; ++t) {
        accum += tensor[TensorIndex(kout, s, t, nhigh)] * shell_power[t];
      }
      kernel[kout * ns + s] = accum;
      eval.pred[kout] += shell_power[s] * accum;
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
      const Real leak_weight = (kout < cfg.nlow) ? 4.0 : 2.0;
      const Real scaled = pred / leak_ref;
      eval.loss += leak_weight * scaled * scaled;
      dloss_dpred = 2.0 * leak_weight * scaled / leak_ref;
    }

    for (int s = 1; s <= nhigh; ++s) {
      eval.grad[s] += dloss_dpred * 2.0 * kernel[kout * ns + s];
    }
  }

  for (int s = 2; s <= nhigh; ++s) {
    const Real prev = std::max(shell_power[s-1], kTiny);
    const Real curr = std::max(shell_power[s], kTiny);
    const Real diff = std::log(curr) - std::log(prev);
    eval.loss += smooth_weight * diff * diff;
    eval.grad[s] += 2.0 * smooth_weight * diff / curr;
    eval.grad[s-1] -= 2.0 * smooth_weight * diff / prev;
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

Real LeakageFraction(const std::vector<Real> &spectrum, const TurbulenceConfig &cfg) {
  Real in_band = 0.0;
  Real out_band = 0.0;
  for (int shell = 1; shell <= 2*cfg.nhigh; ++shell) {
    if (shell >= cfg.nlow && shell <= cfg.nhigh) {
      in_band += spectrum[shell];
    } else {
      out_band += spectrum[shell];
    }
  }
  const Real total = in_band + out_band;
  return (total > kTiny) ? out_band / total : 0.0;
}

ShellFitResult FitClebschShellWeights(const std::vector<Real> &tensor,
                                      const TurbulenceConfig &cfg) {
  ShellFitResult result;
  result.shell_amp.assign(cfg.nhigh + 1, 0.0);

  std::vector<Real> shell_power(cfg.nhigh + 1, 0.0);
  const Real beta_asym = 0.25 * (cfg.expo + 8.0);
  for (int shell = 1; shell <= cfg.nhigh; ++shell) {
    shell_power[shell] = std::pow(static_cast<Real>(shell), -2.0 * beta_asym);
  }
  const Real ref = shell_power[cfg.nlow];
  if (ref > kTiny) {
    for (int shell = 1; shell <= cfg.nhigh; ++shell) shell_power[shell] /= ref;
  }

  ShellFitEval current = EvaluateClebschShellFit(tensor, cfg, shell_power);
  result.initial_loss = current.loss;
  std::vector<Real> best_power = shell_power;
  ShellFitEval best_eval = current;
  Real step_scale = 0.25;

  for (int iter = 0; iter < 250; ++iter) {
    std::vector<Real> grad_log(cfg.nhigh + 1, 0.0);
    Real max_grad = 0.0;
    for (int shell = 1; shell <= cfg.nhigh; ++shell) {
      grad_log[shell] = shell_power[shell] * current.grad[shell];
      max_grad = std::max(max_grad, std::abs(grad_log[shell]));
    }
    if (max_grad < 1.0e-10) break;

    Real step = step_scale / max_grad;
    bool accepted = false;
    for (int trial = 0; trial < 10; ++trial) {
      std::vector<Real> trial_power(cfg.nhigh + 1, 0.0);
      for (int shell = 1; shell <= cfg.nhigh; ++shell) {
        Real log_power = std::log(std::max(shell_power[shell], kTiny));
        log_power -= step * grad_log[shell];
        log_power = std::min(30.0, std::max(-30.0, log_power));
        trial_power[shell] = std::exp(log_power);
      }

      ShellFitEval trial_eval = EvaluateClebschShellFit(tensor, cfg, trial_power);
      if (trial_eval.loss < current.loss) {
        shell_power = std::move(trial_power);
        current = std::move(trial_eval);
        if (current.loss < best_eval.loss) {
          best_power = shell_power;
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

  shell_power = best_power;
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
    const Real scale = std::sqrt(target_mean / pred_mean);
    for (int shell = 1; shell <= cfg.nhigh; ++shell) shell_power[shell] *= scale;
    current = EvaluateClebschShellFit(tensor, cfg, shell_power);
  }

  for (int shell = 1; shell <= cfg.nhigh; ++shell) {
    result.shell_amp[shell] = std::sqrt(std::max(shell_power[shell], 0.0));
  }
  result.predicted_shell_energy = current.pred;
  result.final_loss = current.loss;
  result.leakage_fraction = LeakageFraction(current.pred, cfg);
  result.fitted_slope = FitLogSlope(current.pred, cfg.nlow, cfg.nhigh);
  result.converged = (result.final_loss <= result.initial_loss);

  return result;
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
  std::cout << "  " << label << ": total_modes_in_shell=" << catalog.total_modes
            << ", kept_via_importance_sampling=" << catalog.KeptModes() << std::endl;
  if (cfg.method == VelocityMethod::StreamFunction && cfg.stream_2d) {
    std::cout << "  " << label << ": stream-function coefficients preserve "
              << "E_u(k) ~ k^(-expo) through |psi_k| ~ k^{-(expo+3)/2}" << std::endl;
  }
}

void InitTurbulentVelocity(MeshBlockPack *pmbp, ParameterInput *pin, Real den,
                           bool projection_true_2d) {
  const TurbulenceConfig cfg = ReadTurbulenceConfig(pin, pmbp->pmesh, projection_true_2d);

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
                << "  true_2d=" << (cfg.projection_true_2d ? "true" : "false")
                << ", active_velocity_components=" << (cfg.active_v3 ? 3 : 2)
                << std::endl
                << "  k_crit=" << cfg.k_crit
                << " (wavenumber=" << cfg.k_crit_mag << ")" << std::endl;
      LogModeSummary(cfg, catalog, "velocity");
    }

    const VectorCoefficients coeffs = GenerateProjectionCoefficients(cfg, catalog, &rstate);
    SynthesizeVelocityFromVectorModes(pmbp, cfg, den, catalog, coeffs);
    NormalizeVelocityField(pmbp, cfg);
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
    NormalizeVelocityField(pmbp, cfg);
    return;
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

  const std::vector<int> shell_counts1 = CountModesPerShell(catalog1, cfg.nhigh);
  const std::vector<int> shell_counts2 = CountModesPerShell(catalog2, cfg.nhigh);
  EnsureShellCoverage(shell_counts1, cfg, "phi_1");
  EnsureShellCoverage(shell_counts2, cfg, "phi_2");

  std::vector<Real> tensor = BuildClebschShellResponseTensor(cfg, catalog1, catalog2);
  SymmetrizeClebschTensor(&tensor, cfg.nhigh);
  EnsureReachableTargetBand(tensor, cfg);
  const ShellFitResult fit = FitClebschShellWeights(tensor, cfg);

  if (global_variable::my_rank == 0) {
    std::cout << "  stream_geometry=3D (Clebsch form v = grad(phi1) x grad(phi2))"
              << std::endl;
    LogModeSummary(cfg, catalog1, "phi_1");
    LogModeSummary(cfg, catalog2, "phi_2");
    std::cout << "  3D stream shell fit:"
              << " initial_loss=" << fit.initial_loss
              << ", final_loss=" << fit.final_loss
              << ", fitted_in_band_slope=" << fit.fitted_slope
              << ", estimated_out_of_band_leakage=" << fit.leakage_fraction
              << std::endl;
    if (!fit.converged || std::abs(fit.fitted_slope + cfg.expo) > 0.5 ||
        fit.leakage_fraction > 0.05) {
      std::cout << "  WARNING: 3D stream-function shell fit only matches the requested "
                << "velocity spectrum approximately for this band/seed." << std::endl;
    }
  }

  const ScalarCoefficients phi1_coeffs =
      GenerateShellWeightedScalarCoefficients(catalog1, fit.shell_amp, &rstate1);
  const ScalarCoefficients phi2_coeffs =
      GenerateShellWeightedScalarCoefficients(catalog2, fit.shell_amp, &rstate2);
  SynthesizeClebschVelocity(pmbp, cfg, den, catalog1, phi1_coeffs, catalog2, phi2_coeffs);
  NormalizeVelocityField(pmbp, cfg);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar mixing
//
//  Initializes density, turbulent velocity, and a single passive scalar.
//  Scalar is initialized to scalar_init (default 0.5) to allow symmetric fluctuations.
//  Mean gradient forcing ds/dt = G*v_x is applied via UserSourceTerm.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  const int nscalars = pmbp->phydro->nscalars;
  const int nhydro = pmbp->phydro->nhydro;
  const bool true_2d = UseTrue2DVelocity(pin, pmy_mesh_);

  // Read mean gradient parameter
  mean_gradient_G = pin->GetOrAddReal("problem", "mean_gradient", 1.0);

  // Enroll user source term
  user_srcs_func = MeanGradientForcing;

  if (restart) return;

  // Capture variables for kernel
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  const auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  const Real gm1 = eos.gamma - 1.0;

  // Read initial conditions
  const Real den = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real temp = pin->GetOrAddReal("problem", "temp0", 1.0);
  const Real scalar_init = pin->GetOrAddReal("problem", "scalar_init", 0.5);

  const int nmb = pmbp->nmb_thispack;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar mixing problem:" << std::endl
              << "  mean_gradient G = " << mean_gradient_G << std::endl
              << "  scalar_init = " << scalar_init << std::endl
              << "  true_2d = " << (true_2d ? "true" : "false") << std::endl
              << "  nscalars = " << nscalars << std::endl;
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

  // Initialize passive scalar to offset value (allows symmetric fluctuations)
  par_for("scalar_mix_pgen_scalar", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real rho = u0(m, IDN, k, j, i);
    for (int ns = 0; ns < nscalars; ++ns) {
      u0(m, nhydro+ns, k, j, i) = rho * scalar_init;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn MeanGradientForcing
//  \brief Apply mean gradient forcing to passive scalar: ds/dt = G * v_x
//
//  This simulates the effect of a background mean scalar gradient in the x-direction.
//  Turbulence mixes this gradient, creating scalar fluctuations correlated with velocity.

void MeanGradientForcing(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  const int nhydro = pmbp->phydro->nhydro;
  const Real G = mean_gradient_G;

  par_for("scalar_mix_mean_grad_forcing", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real density = w0(m, IDN, k, j, i);
    const Real vx = w0(m, IVX, k, j, i);
    u0(m, nhydro, k, j, i) += density * G * vx * bdt;
  });
}
