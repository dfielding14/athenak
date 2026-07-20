//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q033_crpai_transport_runtime_local.cpp
//! \brief Bounded source-local CRPAI transport-runtime diagnostic carrier.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

[[noreturn]] void Q033Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q033RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(measured) ||
      std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q033Fatal(label + " does not match the Q-033 bounded local runtime contract");
  }
}

void Q033RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q033Fatal("<" + block + ">/" + name +
              " does not match the Q-033 bounded local runtime contract");
  }
}

void Q033RequireInteger(ParameterInput *pin, const std::string &block,
                        const std::string &name, const int expected) {
  if (pin->GetInteger(block, name) != expected) {
    Q033Fatal("<" + block + ">/" + name +
              " does not match the Q-033 bounded local runtime contract");
  }
}

void Q033RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q033Fatal("<" + block + ">/" + name +
              " does not match the Q-033 bounded local runtime contract");
  }
}

void Q033RequireReal(ParameterInput *pin, const std::string &block,
                     const std::string &name, const Real expected) {
  Q033RequireClose("<" + block + ">/" + name, pin->GetReal(block, name), expected);
}

KOKKOS_INLINE_FUNCTION
std::uint64_t Q033SplitMix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30))*0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27))*0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

KOKKOS_INLINE_FUNCTION
Real Q033Uniform01(const int seed, const int item, const int component) {
  std::uint64_t key = static_cast<std::uint64_t>(seed);
  key ^= (static_cast<std::uint64_t>(item + 1))*0xbf58476d1ce4e5b9ULL;
  key ^= (static_cast<std::uint64_t>(component + 1))*0xd2b74407b1ce6e93ULL;
  const std::uint64_t bits = Q033SplitMix64(key);
  return static_cast<Real>(bits >> 11)*
         static_cast<Real>(1.0/9007199254740992.0);
}

KOKKOS_INLINE_FUNCTION
Real Q033WavePhase(const int seed, const int mode) {
  return static_cast<Real>(2.0*M_PI)*Q033Uniform01(seed, mode, 0);
}

KOKKOS_INLINE_FUNCTION
Real Q033WaveAmplitude(const Real amplitude, const int mode) {
  return amplitude/static_cast<Real>(mode*mode);
}

KOKKOS_INLINE_FUNCTION
Real Q033WaveBy(const Real x1, const Real length, const Real amplitude,
                const int seed) {
  Real value = 0.0;
  for (int mode = 1; mode <= 3; ++mode) {
    const Real phase = static_cast<Real>(2.0*M_PI*mode)*x1/length +
                       Q033WavePhase(seed, mode);
    value += Q033WaveAmplitude(amplitude, mode)*std::cos(phase);
  }
  return value;
}

KOKKOS_INLINE_FUNCTION
Real Q033WaveBz(const Real x1, const Real length, const Real amplitude,
                const int seed) {
  Real value = 0.0;
  for (int mode = 1; mode <= 3; ++mode) {
    const Real phase = static_cast<Real>(2.0*M_PI*mode)*x1/length +
                       Q033WavePhase(seed, mode);
    const Real handedness = (mode % 2 == 0) ? -1.0 : 1.0;
    value += handedness*Q033WaveAmplitude(amplitude, mode)*std::sin(phase);
  }
  return value;
}

}  // namespace

void ProblemGenerator::Q033CRPAITransportRuntimeLocal(ParameterInput *pin,
                                                       const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q033Fatal("q033_crpai_transport_runtime_local requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q033Fatal("q033_crpai_transport_runtime_local requires ideal-MHD EOS");
  }
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic ||
      pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
      pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
      pin->DoesBlockExist("z4c") || pin->DoesBlockExist("turb_driving") ||
      pin->DoesBlockExist("initial_turb")) {
    Q033Fatal("q033_crpai_transport_runtime_local requires the exact Newtonian "
              "single-fluid MHD source task path");
  }
  if (pmy_mesh_->multilevel) {
    Q033Fatal("q033_crpai_transport_runtime_local rejects AMR/SMR");
  }
  if (!pmy_mesh_->two_d || pmy_mesh_->one_d || pmy_mesh_->three_d ||
      pmy_mesh_->mesh_indcs.nx1 != 32 || pmy_mesh_->mesh_indcs.nx2 != 4 ||
      pmy_mesh_->mesh_indcs.nx3 != 1 || pmy_mesh_->mb_indcs.nx1 != 32 ||
      pmy_mesh_->mb_indcs.nx2 != 4 || pmy_mesh_->mb_indcs.nx3 != 1 ||
      pmbp->nmb_thispack != 1) {
    Q033Fatal("q033_crpai_transport_runtime_local requires the exact serial "
              "32x4x1 thin 2D3V one-meshblock carrier");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q033Fatal("q033_crpai_transport_runtime_local requires periodic boundaries");
  }

  Q033RequireString(pin, "time", "integrator", "rk1");
  Q033RequireReal(pin, "time", "cfl_number", 0.1);
  Q033RequireString(pin, "mhd", "eos", "ideal");
  Q033RequireString(pin, "mhd", "reconstruct", "plm");
  Q033RequireString(pin, "mhd", "rsolver", "llf");
  Q033RequireReal(pin, "mhd", "gamma", 1.66666666667);

  Q033RequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q033RequireReal(pin, "particles", "ppc", 4.0);
  Q033RequireString(pin, "particles", "pusher", "boris_tsc");
  Q033RequireInteger(pin, "particles", "nspecies", 1);
  Q033RequireString(pin, "particles", "cr_distribution", "center");
  Q033RequireBoolean(pin, "particles", "deposit_moments", true);
  Q033RequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q033RequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q033RequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", true);
  Q033RequireString(pin, "particles", "pic_background_mode", "coupled");
  Q033RequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q033RequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q033RequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q033RequireReal(pin, "particles", "pic_cr_light_speed", 300.0);
  Q033RequireString(pin, "particles", "pic_cr_initial_state", "momentum");
  Q033RequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q033RequireString(pin, "particles", "pic_wave_damping_mode",
                    "ion_neutral_friction");
  Q033RequireReal(pin, "particles", "pic_ion_neutral_collision_rate", 1.0e-4);
  Q033RequireString(pin, "particles", "pic_deltaf_mode", "physical");
  Q033RequireString(pin, "particles", "pic_deltaf_f0", "kappa_aniso");
  Q033RequireString(pin, "particles", "pic_deltaf_adapt_mode",
                    "global_bikappa_moments_experimental");
  Q033RequireString(pin, "particles", "pic_expanding_box_mode", "on");
  Q033RequireString(pin, "particles", "pic_expansion_law", "exponential");
  Q033RequireReal(pin, "particles", "pic_expansion_rate_x1", 0.0);
  Q033RequireReal(pin, "particles", "pic_expansion_rate_x2", 1.0e-5);
  Q033RequireReal(pin, "particles", "pic_expansion_rate_x3", 1.0e-5);
  Q033RequireBoolean(pin, "particles", "track_displacement", true);
  Q033RequireReal(pin, "species0", "mass", 1.0);
  Q033RequireReal(pin, "species0", "charge", 1.0);

  const std::string block = "q033_crpai_transport_runtime_local";
  Q033RequireString(pin, block, "campaign_id", "Q033-CRPAI-TRANSPORT-RUNTIME-LOCAL");
  Q033RequireString(pin, block, "deck_role",
                    "bounded_source_local_runtime_diagnostic_only_"
                    "not_transport_calibration");
  Q033RequireString(pin, block, "qualification_effect", "none");
  Q033RequireString(pin, block, "frontier_authorization", "not_bound");
  Q033RequireString(pin, block, "q022_closure", "not_claimed");
  Q033RequireString(pin, block, "physical_calibration", "not_claimed");
  Q033RequireString(pin, block, "external_review", "not_claimed");
  Q033RequireString(pin, block, "momentum_distribution",
                    "deterministic_antipodal_bounded_prolate_lattice");
  Q033RequireString(pin, block, "wave_spectrum",
                    "deterministic_seeded_three_mode_transverse_ct_carrier");
  Q033RequireInteger(pin, block, "wave_mode_count", 3);
  const int wave_seed = pin->GetInteger(block, "wave_seed");
  const int momentum_seed = pin->GetInteger(block, "momentum_seed");
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real guide_field = pin->GetReal(block, "guide_field");
  const Real wave_amplitude = pin->GetReal(block, "wave_amplitude");
  const Real momentum_p0 = pin->GetReal(block, "momentum_p0");
  const Real momentum_xi = pin->GetReal(block, "momentum_xi");
  if (wave_seed != 330033 || momentum_seed != 330034 || rho <= 0.0 ||
      pressure <= 0.0 || guide_field <= 0.0 || wave_amplitude <= 0.0 ||
      momentum_p0 <= 0.0 || momentum_xi <= 1.0) {
    Q033Fatal("q033_crpai_transport_runtime_local diagnostic seed contract is invalid");
  }

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  const Real length = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;

  par_for("pgen_q033_local_ct_field", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real by = Q033WaveBy(x1, length, wave_amplitude, wave_seed);
    const Real bz = Q033WaveBz(x1, length, wave_amplitude, wave_seed);
    b0.x1f(m, k, j, i) = guide_field;
    b0.x2f(m, k, j, i) = by;
    b0.x3f(m, k, j, i) = bz;
    if (i == ie) b0.x1f(m, k, j, i + 1) = guide_field;
    if (j == je) b0.x2f(m, k, j + 1, i) = by;
    if (k == ke) b0.x3f(m, k + 1, j, i) = bz;
  });

  par_for("pgen_q033_local_primitives", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m, IDN, k, j, i) = rho;
    w0(m, IVX, k, j, i) = 0.0;
    w0(m, IVY, k, j, i) = 0.0;
    w0(m, IVZ, k, j, i) = 0.0;
    w0(m, IEN, k, j, i) = pressure;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);

  auto *ppart = pmbp->ppart;
  ppart->pic_deltaf_adaptive_xi = momentum_xi;
  ppart->pic_deltaf_adaptive_p0 = momentum_p0;
  auto &pi = ppart->prtcl_idata;
  auto &pr = ppart->prtcl_rdata;
  const Real kappa = ppart->pic_deltaf_kappa;
  par_for("pgen_q033_local_antipodal_prolate_momenta", DevExeSpace(),
          0, ppart->nprtcl_thispack - 1,
  KOKKOS_LAMBDA(const int p) {
    const int pair = pi(PTAG, p)/2;
    const Real sign = ((pi(PTAG, p) % 2) == 0) ? 1.0 : -1.0;
    const Real u_parallel = Q033Uniform01(momentum_seed, pair, 0);
    const Real u_perp = Q033Uniform01(momentum_seed, pair, 1);
    const Real phi = static_cast<Real>(2.0*M_PI)*
                     Q033Uniform01(momentum_seed, pair, 2);
    const Real parallel = sign*momentum_p0*momentum_xi*(0.25 + 0.75*u_parallel);
    const Real perp = sign*momentum_p0*std::sqrt(u_perp/momentum_xi);
    pr(IPVX, p) = parallel;
    pr(IPVY, p) = perp*std::cos(phi);
    pr(IPVZ, p) = perp*std::sin(phi);
    pr(IPF0, p) = particles::PICAdaptiveDeltaFBackgroundValue(
        kappa, momentum_p0, momentum_p0, momentum_xi, 1.0, 1.0, 1.0,
        pr(IPVX, p), pr(IPVY, p), pr(IPVZ, p));
    pr(IPDFWT, p) = 0.0;
  });
  Kokkos::fence();
}
