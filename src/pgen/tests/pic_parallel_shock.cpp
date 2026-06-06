//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pic_parallel_shock.cpp
//! \brief Parallel-shock MHD-PIC benchmark problem generator.
//!
//! Minimal AthenaK adaptation of the Section 5.4 parallel-shock setup in:
//!   Sun & Bai (2022), "The MHD-PIC Module in Athena++"
//! using:
//! - reflecting wall at inner x1
//! - inflow at outer x1
//! - background field B0 parallel to x
//! - eta-based shock-surface CR injection with conservative gas subtraction
//! - AMR refinement based on density/pressure curvature metrics.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "outputs/restart_utils.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

struct ShockCell {
  int m, k, j, i;
  int gid;
  Real x1c, x2c, x3c;
  Real dx1, dx2, dx3;
  Real area;
  Real vol;
};

struct GlobalShockCell {
  int owner_rank;
  int local_cell_index;
  int gid;
  Real x1c, x2c, x3c;
  Real dx1, dx2, dx3;
  Real area;
};

struct InjectedParticle {
  std::int64_t tag;
  int gid;
  int m, k, j, i;
  Real vol;
  Real x1, x2, x3;
  Real vx, vy, vz;
};

struct GasDelta {
  int m, k, j, i;
  Real vol;
  Real dm;
  Real dmx;
  Real dmy;
  Real dmz;
  Real de;
};

// Runtime controls populated by ProblemGenerator::PICParallelShock().
Real ps_rho0 = 1.0;
Real ps_p0 = 1.0;
Real ps_u0 = 10.0;
Real ps_b0 = 1.0;
Real ps_eta = 1.0e-3;
Real ps_vinj_over_u0 = std::sqrt(10.0);
Real ps_shock_speed = 0.0;
Real ps_xshock0 = 0.0;
Real ps_inject_half_width_cells = 0.5;
Real ps_inject_t_start = 0.0;
Real ps_inject_t_stop = 1.0e99;
Real ps_remove_birth_time_before = -1.0;
Real ps_seed_noise_amp = 0.0;
int ps_seed_noise_seed = 1234;
std::array<Real, 4> ps_seed_noise_phase_by = {0.0, 0.0, 0.0, 0.0};
std::array<Real, 4> ps_seed_noise_phase_bz = {0.0, 0.0, 0.0, 0.0};
enum class PSShockSpeedModel { finite_mach, ideal_surface };
PSShockSpeedModel ps_shock_speed_model = PSShockSpeedModel::finite_mach;
Real ps_refine_curv = 1.0;
Real ps_derefine_curv = 0.1;
Real ps_rho_floor_frac = 1.0e-6;
Real ps_p_floor_frac = 1.0e-8;
bool ps_enable_injection = true;
bool ps_enable_subtraction = true;
bool ps_enable_curvature_amr = true;
bool ps_enable_frame_tracking = false;
bool ps_use_2d3v = false;
enum class PSFrameMode { velocity, recenter };
PSFrameMode ps_frame_mode = PSFrameMode::velocity;
Real ps_frame_t_start = 0.0;
Real ps_frame_t_ramp = 0.0;
Real ps_frame_vfrac = 1.0;
Real ps_frame_dv_max = 1.0e99;
bool ps_frame_apply_to_particles = true;
bool ps_frame_apply_to_inflow = true;
bool ps_frame_require_uniform = true;
int ps_frame_diag_dcycle = 200;
int ps_feedback_diag_dcycle = 0;
Real ps_recenter_x_target = 2.0;
Real ps_recenter_x_trigger = 3.0;
Real ps_recenter_dx1 = 0.0;
int ps_recenter_shift_cells = 2;
Real ps_recenter_vshock_model = -1.0;
int ps_inject_species = 0;
int ps_inject_seed = 1234;
Real ps_particle_mass = 1.0;
Real ps_particle_charge = 1.0;
Real ps_particle_q_over_m = 1.0;
Real ps_particle_macro_mass = 1.0;
Real ps_mass_reservoir_global = 0.0;
int ps_injection_transaction_cycle = std::numeric_limits<int>::min();
std::vector<GasDelta> ps_injection_transaction_gas_deltas;
std::array<Real, 5> ps_injection_transaction_expected_global = {};
std::array<Real, 5> ps_injection_transaction_applied_local = {};
std::array<Real, 5> ps_injection_transaction_abs_global = {};
Real ps_injection_transaction_terms_global = 1.0;
bool ps_test_source_transaction_terms_override = false;
Real ps_test_source_transaction_terms = 1.0;
Real ps_injected_cr_count_global = 0.0;
Real ps_injected_cr_mass_global = 0.0;
Real ps_injected_cr_momentum_x1_global = 0.0;
Real ps_injected_cr_momentum_x2_global = 0.0;
Real ps_injected_cr_momentum_x3_global = 0.0;
Real ps_injected_cr_energy_global = 0.0;
bool ps_removed_excluded_early_cohort = false;
Real ps_removed_cr_count_global = 0.0;
Real ps_removed_cr_mass_global = 0.0;
Real ps_removed_cr_momentum_x1_global = 0.0;
Real ps_removed_cr_momentum_x2_global = 0.0;
Real ps_removed_cr_momentum_x3_global = 0.0;
Real ps_removed_cr_energy_global = 0.0;
bool ps_escape_ledger_complete = true;
int ps_escape_audit_calls = 0;
Real ps_escape_last_audit_time = 0.0;
Real ps_escaped_injected_cr_count_global = 0.0;
Real ps_escaped_injected_cr_mass_global = 0.0;
Real ps_escaped_injected_cr_momentum_x1_global = 0.0;
Real ps_escaped_injected_cr_momentum_x2_global = 0.0;
Real ps_escaped_injected_cr_momentum_x3_global = 0.0;
Real ps_escaped_injected_cr_energy_global = 0.0;
Real ps_escaped_initial_cr_count_global = 0.0;
std::int64_t ps_particle_population_audit_calls = 0;
bool ps_cr_ledger_complete = true;
bool ps_tag_seeded = false;
bool ps_tag_progression_validated = false;
std::int64_t ps_injection_tag_floor = 0;
std::int64_t ps_next_tag = 0;
ParameterInput *ps_pin = nullptr;

void HashParallelShockRestartBytes(std::uint64_t &hash, const void *data,
                                   const std::size_t size) {
  constexpr std::uint64_t fnv_prime = 1099511628211ULL;
  const auto *bytes = static_cast<const unsigned char *>(data);
  for (std::size_t n = 0; n < size; ++n) {
    hash ^= static_cast<std::uint64_t>(bytes[n]);
    hash *= fnv_prime;
  }
}

template <typename T>
void HashParallelShockRestartControl(std::uint64_t &hash, const char *name,
                                     const T &value) {
  while (*name != '\0') {
    HashParallelShockRestartBytes(hash, name, 1);
    ++name;
  }
  constexpr unsigned char separator = 0;
  HashParallelShockRestartBytes(hash, &separator, sizeof(separator));
  HashParallelShockRestartBytes(hash, &value, sizeof(value));
}

std::string ParallelShockRestartControlFingerprint() {
  constexpr std::uint64_t fnv_offset_basis = 14695981039346656037ULL;
  std::uint64_t hash = fnv_offset_basis;
  constexpr char schema[] = "athenak_pic_parallel_shock_restart_controls_v1";
  HashParallelShockRestartBytes(hash, schema, sizeof(schema));

  // Diagnostics cadence and initial-only seed noise do not alter continuation.
  HashParallelShockRestartControl(hash, "ps_rho0", ps_rho0);
  HashParallelShockRestartControl(hash, "ps_p0", ps_p0);
  HashParallelShockRestartControl(hash, "ps_u0", ps_u0);
  HashParallelShockRestartControl(hash, "ps_b0", ps_b0);
  HashParallelShockRestartControl(hash, "ps_eta", ps_eta);
  HashParallelShockRestartControl(hash, "ps_vinj_over_u0", ps_vinj_over_u0);
  HashParallelShockRestartControl(hash, "ps_inject_half_width_cells",
                                  ps_inject_half_width_cells);
  HashParallelShockRestartControl(hash, "ps_inject_t_start", ps_inject_t_start);
  HashParallelShockRestartControl(hash, "ps_inject_t_stop", ps_inject_t_stop);
  HashParallelShockRestartControl(hash, "ps_remove_birth_time_before",
                                  ps_remove_birth_time_before);
  HashParallelShockRestartControl(hash, "ps_shock_speed_model",
                                  static_cast<int>(ps_shock_speed_model));
  HashParallelShockRestartControl(hash, "ps_refine_curv", ps_refine_curv);
  HashParallelShockRestartControl(hash, "ps_derefine_curv", ps_derefine_curv);
  HashParallelShockRestartControl(hash, "ps_rho_floor_frac", ps_rho_floor_frac);
  HashParallelShockRestartControl(hash, "ps_p_floor_frac", ps_p_floor_frac);
  HashParallelShockRestartControl(hash, "ps_enable_injection",
                                  static_cast<int>(ps_enable_injection));
  HashParallelShockRestartControl(hash, "ps_enable_subtraction",
                                  static_cast<int>(ps_enable_subtraction));
  HashParallelShockRestartControl(hash, "ps_enable_curvature_amr",
                                  static_cast<int>(ps_enable_curvature_amr));
  HashParallelShockRestartControl(hash, "ps_test_source_transaction_terms_override",
                                  static_cast<int>(
                                      ps_test_source_transaction_terms_override));
  HashParallelShockRestartControl(hash, "ps_test_source_transaction_terms",
                                  ps_test_source_transaction_terms);
  HashParallelShockRestartControl(hash, "ps_inject_species", ps_inject_species);
  HashParallelShockRestartControl(hash, "ps_inject_seed", ps_inject_seed);
  HashParallelShockRestartControl(hash, "ps_particle_mass", ps_particle_mass);
  HashParallelShockRestartControl(hash, "ps_particle_charge", ps_particle_charge);
  HashParallelShockRestartControl(hash, "ps_particle_q_over_m",
                                  ps_particle_q_over_m);
  HashParallelShockRestartControl(hash, "ps_particle_macro_mass",
                                  ps_particle_macro_mass);

  HashParallelShockRestartControl(hash, "ps_enable_frame_tracking",
                                  static_cast<int>(ps_enable_frame_tracking));
  HashParallelShockRestartControl(hash, "ps_frame_mode",
                                  static_cast<int>(ps_frame_mode));
  HashParallelShockRestartControl(hash, "ps_frame_t_start", ps_frame_t_start);
  HashParallelShockRestartControl(hash, "ps_frame_t_ramp", ps_frame_t_ramp);
  HashParallelShockRestartControl(hash, "ps_frame_vfrac", ps_frame_vfrac);
  HashParallelShockRestartControl(hash, "ps_frame_dv_max", ps_frame_dv_max);
  HashParallelShockRestartControl(hash, "ps_frame_apply_to_particles",
                                  static_cast<int>(ps_frame_apply_to_particles));
  HashParallelShockRestartControl(hash, "ps_frame_apply_to_inflow",
                                  static_cast<int>(ps_frame_apply_to_inflow));
  HashParallelShockRestartControl(hash, "ps_frame_require_uniform",
                                  static_cast<int>(ps_frame_require_uniform));
  HashParallelShockRestartControl(hash, "ps_recenter_x_target",
                                  ps_recenter_x_target);
  HashParallelShockRestartControl(hash, "ps_recenter_x_trigger",
                                  ps_recenter_x_trigger);
  HashParallelShockRestartControl(hash, "ps_recenter_dx1", ps_recenter_dx1);
  HashParallelShockRestartControl(hash, "ps_recenter_shift_cells",
                                  ps_recenter_shift_cells);
  HashParallelShockRestartControl(hash, "ps_recenter_vshock_model",
                                  ps_recenter_vshock_model);

  HashParallelShockRestartControl(hash, "ps_shock_speed", ps_shock_speed);
  HashParallelShockRestartControl(hash, "ps_xshock0", ps_xshock0);
  HashParallelShockRestartControl(hash, "ps_use_2d3v",
                                  static_cast<int>(ps_use_2d3v));

  std::ostringstream fingerprint;
  fingerprint << "v1:" << std::hex << std::setfill('0') << std::setw(16) << hash;
  return fingerprint.str();
}

void ValidateAndStoreParallelShockRestartControls(ParameterInput *pin,
                                                  const bool restart) {
  constexpr char block[] = "problem";
  constexpr char parameter[] = "ps_restart_control_fingerprint";
  const std::string current = ParallelShockRestartControlFingerprint();
  if (restart) {
    if (!pin->DoesParameterExist(block, parameter)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart metadata is missing the "
                << "continuation-control fingerprint." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    const std::string checkpointed = pin->GetString(block, parameter);
    if (checkpointed != current) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart continuation-control fingerprint "
                << "mismatch; injection/frame controls must not be overridden."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  pin->SetString(block, parameter, current);
}

inline bool FrameModeVelocity() {
  return ps_frame_mode == PSFrameMode::velocity;
}

inline bool FrameModeRecenter() {
  return ps_frame_mode == PSFrameMode::recenter;
}

inline Real FrameRampFactor(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  if (!ps_enable_frame_tracking) return 0.0;
  if (t <= ps_frame_t_start) return 0.0;
  if (ps_frame_t_ramp <= 0.0) return 1.0;

  const Real s = (t - ps_frame_t_start)/ps_frame_t_ramp;
  if (s >= 1.0) return 1.0;
  if (s <= 0.0) return 0.0;
  return 0.5*(1.0 - std::cos(M_PI*s));
}

inline Real FrameVelocityOffset(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  const Real vfollow = ps_frame_vfrac*ps_shock_speed;
  return -vfollow*FrameRampFactor(t);
}

inline Real FrameOffsetDisplacement(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  if (!ps_enable_frame_tracking) return 0.0;
  if (t <= ps_frame_t_start) return 0.0;

  const Real vfollow = ps_frame_vfrac*ps_shock_speed;
  if (ps_frame_t_ramp <= 0.0) {
    return -vfollow*(t - ps_frame_t_start);
  }

  const Real tr = ps_frame_t_ramp;
  const Real dt = t - ps_frame_t_start;
  if (dt >= tr) {
    const Real iramp = 0.5*tr;
    const Real itail = dt - tr;
    return -vfollow*(iramp + itail);
  }

  const Real s = dt/tr;
  const Real iramp = tr*(0.5*s - std::sin(M_PI*s)/(2.0*M_PI));
  return -vfollow*iramp;
}

inline Real ShockSurfaceUnshiftedX1(const Real t) {
  if (FrameModeRecenter()) {
    Real vmodel = ps_shock_speed;
    if (ps_recenter_vshock_model > 0.0) vmodel = ps_recenter_vshock_model;
    return ps_xshock0 + vmodel*t;
  }
  return ps_xshock0 + ps_shock_speed*t + FrameOffsetDisplacement(t);
}

inline int RecenterEventCount(const Real t) {
  if (!ps_enable_frame_tracking) return 0;
  if (!FrameModeRecenter()) return 0;
  if (ps_recenter_shift_cells < 1 || ps_recenter_dx1 <= 0.0) return 0;

  const Real xu = ShockSurfaceUnshiftedX1(t);
  const Real dx_shift = ps_recenter_shift_cells*ps_recenter_dx1;
  const Real trigger = ps_recenter_x_trigger;
  const Real target = ps_recenter_x_target;
  if (xu <= trigger) return 0;
  const Real overshoot = xu - target;
  if (overshoot <= 0.0) return 0;
  return static_cast<int>(std::ceil(overshoot/dx_shift));
}

inline Real RecenterDisplacement(const Real t) {
  if (!ps_enable_frame_tracking) return 0.0;
  if (!FrameModeRecenter()) return 0.0;
  const int nevents = RecenterEventCount(t);
  return static_cast<Real>(nevents*ps_recenter_shift_cells)*ps_recenter_dx1;
}

inline Real ShockSurfaceModelX1(const Real t) {
  return ShockSurfaceUnshiftedX1(t) - RecenterDisplacement(t);
}

inline std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline Real UniformFromUint64(std::uint64_t x) {
  constexpr Real inv = 1.0/static_cast<Real>(std::numeric_limits<std::uint64_t>::max());
  return static_cast<Real>(x)*inv;
}

inline void ConfigureSeedNoisePhases() {
  std::uint64_t state = static_cast<std::uint64_t>(ps_seed_noise_seed);
  for (int n = 0; n < 4; ++n) {
    state = SplitMix64(state + 0x9e3779b97f4a7c15ULL);
    ps_seed_noise_phase_by[n] = 2.0*M_PI*UniformFromUint64(state);
    state = SplitMix64(state + 0xbf58476d1ce4e5b9ULL);
    ps_seed_noise_phase_bz[n] = 2.0*M_PI*UniformFromUint64(state);
  }
}

inline Real EstimateShockSpeed(const Real gamma, const Real rho0, const Real p0,
                               const Real u0) {
  // Finite-Mach hydrodynamic estimate for a piston-driven shock reflected off
  // the inner wall. For B parallel to shock normal (Bx-only), magnetic terms
  // do not change the 1D compression ratio in this minimal setup.
  const Real cs2 = gamma*p0/rho0;
  if (cs2 <= 0.0) return 0.0;
  const Real ms2 = SQR(u0)/cs2;
  if (ms2 <= 1.0) return 0.0;
  const Real r = ((gamma + 1.0)*ms2)/((gamma - 1.0)*ms2 + 2.0);
  if (r <= 1.0) return 0.0;
  return u0/(r - 1.0);
}

inline Real IdealSurfaceShockSpeed(const Real gamma, const Real u0) {
  // Sun & Bai Section 5.4 injects at x = u_sh' t with
  // u_sh' = (Gamma - 1) u0 / 2 in the reflecting-wall frame.
  if (gamma <= 1.0 || u0 <= 0.0) return 0.0;
  return 0.5*(gamma - 1.0)*u0;
}

void EncodeCRStateFromVelocity(const particles::Particles *ppart,
                               const Real vx, const Real vy, const Real vz,
                               Real &state_x, Real &state_y, Real &state_z) {
  state_x = vx;
  state_y = vy;
  state_z = vz;
  if (!ppart->UsesRelativisticCRState()) return;

  const Real light_speed = ppart->pic_cr_light_speed;
  const Real v2 = vx*vx + vy*vy + vz*vz;
  if (v2 >= light_speed*light_speed) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injected velocity must satisfy |v| < "
              << "<particles>/pic_cr_light_speed in a momentum-state mode."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const Real gamma = 1.0/std::sqrt(1.0 - v2/(light_speed*light_speed));
  state_x *= gamma;
  state_y *= gamma;
  state_z *= gamma;
}

Real VelocityMagnitudeFromMomentumMagnitude(const particles::Particles *ppart,
                                            const Real momentum) {
  if (!ppart->UsesRelativisticCRState()) return momentum;
  const Real light_speed = ppart->pic_cr_light_speed;
  return momentum/std::sqrt(1.0 + SQR(momentum/light_speed));
}

void BoostRelativeVelocityFromSurface(const particles::Particles *ppart,
                                      const Real surface_vx,
                                      const Real relative_vx,
                                      const Real relative_vy,
                                      const Real relative_vz,
                                      Real &vx, Real &vy, Real &vz) {
  vx = surface_vx + relative_vx;
  vy = relative_vy;
  vz = relative_vz;
  if (!ppart->UsesRelativisticCRState()) return;

  const Real light_speed = ppart->pic_cr_light_speed;
  if (std::abs(surface_vx) >= light_speed) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock shock-surface speed must satisfy |v| < "
              << "<particles>/pic_cr_light_speed in a momentum-state mode."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const Real beta2 = SQR(surface_vx/light_speed);
  const Real gamma_surface = 1.0/std::sqrt(1.0 - beta2);
  const Real denominator = 1.0 + surface_vx*relative_vx/SQR(light_speed);
  if (!(denominator > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock velocity boost has a non-positive "
              << "denominator." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  vx = (surface_vx + relative_vx)/denominator;
  vy = relative_vy/(gamma_surface*denominator);
  vz = relative_vz/(gamma_surface*denominator);
}

void ApplyFrameShiftToFluid(Mesh *pm, const Real dvx) {
  if (dvx == 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto *pmhd = pmbp->pmhd;
  auto &u0 = pmhd->u0;
  auto &w0 = pmhd->w0;
  auto &b0 = pmhd->b0;
  auto &bcc0 = pmhd->bcc0;

  pmhd->peos->ConsToPrim(u0, b0, w0, bcc0, false, is, ie, js, je, ks, ke);
  par_for("ps_frame_shift_w0", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w0(m, IVX, k, j, i) += dvx;
  });
  pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}

void ApplyFrameShiftToParticles(Mesh *pm, const Real dvx) {
  if (dvx == 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;

  auto *ppart = pmbp->ppart;
  if (ppart->nprtcl_thispack <= 0) return;
  if (ppart->UsesRelativisticCRState()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle frame shifting is not yet defined "
              << "for momentum-state modes; disable "
              << "<problem>/ps_frame_apply_to_particles." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  auto &prtcl_rdata = ppart->prtcl_rdata;
  par_for("ps_frame_shift_particles", DevExeSpace(), 0, ppart->nprtcl_thispack - 1,
  KOKKOS_LAMBDA(const int p) {
    prtcl_rdata(IPVX, p) += dvx;
  });
}

void FatalParticleMigrationError(const char *message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl
            << "pic_parallel_shock recenter particle migration failed: "
            << message << std::endl;
  restart_utils::AbortOnFatalError();
}

void CopyPackedParticleDataToDevice(particles::Particles *ppart,
                                    const HostArray2D<int> &h_pi_source,
                                    const HostArray2D<Real> &h_pr_source, const int npart,
                                    const int destination_offset) {
  if (npart <= 0) return;

  const int ni = ppart->nidata;
  const int nr = ppart->nrdata;
  HostArray2D<int> h_pi_packed("ps_pi_packed", ni, npart);
  HostArray2D<Real> h_pr_packed("ps_pr_packed", nr, npart);
  for (int p = 0; p < npart; ++p) {
    for (int q = 0; q < ni; ++q) h_pi_packed(q, p) = h_pi_source(q, p);
    for (int q = 0; q < nr; ++q) h_pr_packed(q, p) = h_pr_source(q, p);
  }

  auto d_pi_packed = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_pi_packed);
  auto d_pr_packed = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_pr_packed);
  auto prtcl_idata = ppart->prtcl_idata;
  auto prtcl_rdata = ppart->prtcl_rdata;
  par_for("ps_copy_packed_particles", DevExeSpace(), 0, npart - 1,
  KOKKOS_LAMBDA(const int p) {
    for (int q = 0; q < ni; ++q) {
      prtcl_idata(q, destination_offset + p) = d_pi_packed(q, p);
    }
    for (int q = 0; q < nr; ++q) {
      prtcl_rdata(q, destination_offset + p) = d_pr_packed(q, p);
    }
  });
  Kokkos::fence();
}

Real ParallelShockLedgerTolerance(const Real lhs, const Real rhs,
                                  const Real accumulated_terms = 1.0) {
  const Real scale = std::max({std::abs(lhs), std::abs(rhs),
                               std::abs(ps_particle_macro_mass),
                               static_cast<Real>(1.0)});
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real relative_bound = std::max(accumulated_terms, static_cast<Real>(1.0))*eps;
  const Real summation_bound =
      (relative_bound < 0.5) ? relative_bound/(1.0 - relative_bound) : 1.0;
  return (8.0*summation_bound + 64.0*eps)*scale;
}

bool ParallelShockLedgerValuesAgree(const Real lhs, const Real rhs,
                                    const Real accumulated_terms = 1.0) {
  return std::abs(lhs - rhs) <=
      ParallelShockLedgerTolerance(lhs, rhs, accumulated_terms);
}

bool ParallelShockSourceTransactionValuesAgree(const Real lhs, const Real rhs,
                                               const Real absolute_contributions,
                                               const Real accumulated_terms) {
  if (!std::isfinite(lhs) || !std::isfinite(rhs) ||
      !std::isfinite(absolute_contributions) || absolute_contributions < 0.0 ||
      !std::isfinite(accumulated_terms) || accumulated_terms < 0.0 ||
      !std::isfinite(ps_particle_macro_mass)) {
    return false;
  }
  const Real scale = std::max({std::abs(lhs), std::abs(rhs),
                               std::abs(absolute_contributions),
                               std::abs(ps_particle_macro_mass),
                               static_cast<Real>(1.0)});
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real relative_bound =
      std::max(accumulated_terms, static_cast<Real>(1.0))*eps;
  if (!std::isfinite(scale) || !std::isfinite(relative_bound) ||
      relative_bound >= 0.5) {
    return false;
  }
  const Real summation_bound = relative_bound/(1.0 - relative_bound);
  const Real tolerance = (8.0*summation_bound + 64.0*eps)*scale;
  return std::isfinite(tolerance) && std::abs(lhs - rhs) <= tolerance;
}

bool ParallelShockLedgerValueExceeds(const Real lhs, const Real rhs) {
  return lhs > rhs + ParallelShockLedgerTolerance(lhs, rhs);
}

void RejectDuplicateParallelShockRestartLedgers(ParameterInput *pin,
                                                const bool restart) {
  if (!restart || pin == nullptr) return;
  constexpr std::array<const char *, 19> cr_ledger_fields = {
    "ps_cr_ledger_schema", "ps_cr_ledger_complete", "ps_mass_reservoir_global",
    "ps_injected_cr_count_global", "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global", "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global", "ps_injected_cr_energy_global",
    "ps_removed_excluded_early_cohort", "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global", "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global", "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global", "ps_tag_seeded", "ps_injection_tag_floor",
    "ps_next_tag"
  };
  constexpr std::array<const char *, 11> escape_ledger_fields = {
    "ps_escape_ledger_schema", "ps_escape_ledger_complete",
    "ps_escape_audit_calls", "ps_escape_last_audit_time",
    "ps_escaped_injected_cr_count_global", "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global", "ps_escaped_initial_cr_count_global"
  };
  int duplicate_local = 0;
  for (const char *field : cr_ledger_fields) {
    duplicate_local +=
        pin->ParameterWasDuplicatedInLoadedInput("problem", field) ? 1 : 0;
  }
  for (const char *field : escape_ledger_fields) {
    duplicate_local +=
        pin->ParameterWasDuplicatedInLoadedInput("problem", field) ? 1 : 0;
  }
#if MPI_PARALLEL_ENABLED
  int duplicate_global = 0;
  MPI_Allreduce(&duplicate_local, &duplicate_global, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#else
  const int duplicate_global = duplicate_local;
#endif
  if (duplicate_global != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock restart contains duplicate CR or escape "
              << "ledger metadata." << std::endl;
    restart_utils::AbortOnFatalError();
  }
}

bool PaperVL2EscapeStageChronologyIsValid(const particles::Particles *ppart,
                                          const Mesh *pm, const int stage,
                                          const Real audit_time) {
  if (!ppart->UsesPaperVL2Coupling()) return true;
  if (stage != 1 && stage != 2) return false;
  if (pm->ncycle < 0 ||
      pm->ncycle > (std::numeric_limits<int>::max() - stage + 1)/2) {
    return false;
  }
  const int expected_calls_before_stage = 2*pm->ncycle + stage - 1;
  const Real expected_last_time = (stage == 1) ? pm->time : pm->time + 0.5*pm->dt;
  const Real expected_audit_time =
      (stage == 1) ? pm->time + 0.5*pm->dt : pm->time + pm->dt;
  return ps_escape_audit_calls == expected_calls_before_stage &&
      ParallelShockLedgerValuesAgree(ps_escape_last_audit_time,
                                      expected_last_time) &&
      ParallelShockLedgerValuesAgree(audit_time, expected_audit_time);
}

void ValidatePaperVL2CommittedEscapeChronology(const particles::Particles *ppart,
                                               const int committed_cycle,
                                               const Real committed_time,
                                               const char *context,
                                               const bool collective = true) {
  if (ppart == nullptr || !ppart->UsesPaperVL2Coupling()) return;
  int invalid_local = 0;
  if (committed_cycle < 0 ||
      committed_cycle > std::numeric_limits<int>::max()/2 ||
      !std::isfinite(committed_time) || committed_time < 0.0) {
    invalid_local = 1;
  } else {
    const int expected_calls = 2*committed_cycle;
    const Real expected_last_time = (expected_calls == 0) ? 0.0 : committed_time;
    if (ps_escape_audit_calls != expected_calls ||
        !ParallelShockLedgerValuesAgree(ps_escape_last_audit_time,
                                        expected_last_time)) {
      invalid_local = 1;
    }
  }
#if MPI_PARALLEL_ENABLED
  int invalid_global = invalid_local;
  if (collective) {
    MPI_Allreduce(&invalid_local, &invalid_global, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
  }
#else
  const int invalid_global = invalid_local;
#endif
  if (invalid_global != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock " << context
              << " paper-VL2 escape-audit chronology is invalid." << std::endl;
    restart_utils::AbortOnFatalError();
  }
}

void ValidateParallelShockRuntimeLedger(const char *context, const Real current_time) {
  const auto invalid_nonnegative_ledger = [](const Real value) {
    return !std::isfinite(value) || value < 0.0;
  };
  const auto invalid_count_ledger = [&](const Real value) {
    return invalid_nonnegative_ledger(value) || std::floor(value) != value;
  };
  const bool invalid_removed_before_sink =
      !ps_removed_excluded_early_cohort &&
      (ps_removed_cr_count_global != 0.0 || ps_removed_cr_mass_global != 0.0 ||
       ps_removed_cr_momentum_x1_global != 0.0 ||
       ps_removed_cr_momentum_x2_global != 0.0 ||
       ps_removed_cr_momentum_x3_global != 0.0 ||
       ps_removed_cr_energy_global != 0.0);
  const bool invalid_sink_completion =
      ps_removed_excluded_early_cohort &&
      (ps_remove_birth_time_before < 0.0 ||
       current_time < ps_remove_birth_time_before);
  const Real expected_injected_mass =
      ps_injected_cr_count_global*ps_particle_macro_mass;
  const Real expected_removed_mass =
      ps_removed_cr_count_global*ps_particle_macro_mass;
  const Real expected_escaped_mass =
      ps_escaped_injected_cr_count_global*ps_particle_macro_mass;
  const auto max_tag = static_cast<std::int64_t>(std::numeric_limits<int>::max());
  const bool injected_count_fits_tag =
      std::isfinite(ps_injected_cr_count_global) &&
      ps_injected_cr_count_global >= 0.0 &&
      std::floor(ps_injected_cr_count_global) == ps_injected_cr_count_global &&
      ps_injected_cr_count_global <= static_cast<Real>(max_tag);
  const auto injected_tag_count = injected_count_fits_tag ?
      static_cast<std::int64_t>(ps_injected_cr_count_global) : 0;
  const bool invalid_tag_sum =
      !injected_count_fits_tag || ps_injection_tag_floor < 0 ||
      ps_injection_tag_floor > max_tag ||
      injected_tag_count > max_tag - ps_injection_tag_floor;
  const bool invalid_tag_window =
      invalid_tag_sum || ps_next_tag < 0 || ps_next_tag > max_tag ||
      (!ps_tag_seeded &&
       (ps_injection_tag_floor != 0 || ps_next_tag != 0 ||
        ps_injected_cr_count_global != 0.0)) ||
      (ps_tag_seeded &&
       (ps_injection_tag_floor + injected_tag_count != ps_next_tag));
  const bool invalid_escape_before_audit =
      ps_escape_audit_calls == 0 &&
      (ps_escape_last_audit_time != 0.0 ||
       ps_escaped_injected_cr_count_global != 0.0 ||
       ps_escaped_injected_cr_mass_global != 0.0 ||
       ps_escaped_injected_cr_momentum_x1_global != 0.0 ||
       ps_escaped_injected_cr_momentum_x2_global != 0.0 ||
       ps_escaped_injected_cr_momentum_x3_global != 0.0 ||
       ps_escaped_injected_cr_energy_global != 0.0 ||
       ps_escaped_initial_cr_count_global != 0.0);
  const bool invalid_escape_audit_time =
      ps_escape_audit_calls < 0 || !std::isfinite(ps_escape_last_audit_time) ||
      ps_escape_last_audit_time < 0.0 ||
      ps_escape_last_audit_time >
          current_time + ParallelShockLedgerTolerance(ps_escape_last_audit_time,
                                                       current_time);
  const bool invalid_empty_injected_ledger =
      ps_injected_cr_count_global == 0.0 &&
      (ps_injected_cr_momentum_x1_global != 0.0 ||
       ps_injected_cr_momentum_x2_global != 0.0 ||
       ps_injected_cr_momentum_x3_global != 0.0 ||
       ps_injected_cr_energy_global != 0.0);
  const bool invalid_empty_removed_ledger =
      ps_removed_cr_count_global == 0.0 &&
      (ps_removed_cr_momentum_x1_global != 0.0 ||
       ps_removed_cr_momentum_x2_global != 0.0 ||
       ps_removed_cr_momentum_x3_global != 0.0 ||
       ps_removed_cr_energy_global != 0.0);
  const bool invalid_empty_escape_ledger =
      ps_escaped_injected_cr_count_global == 0.0 &&
      (ps_escaped_injected_cr_momentum_x1_global != 0.0 ||
       ps_escaped_injected_cr_momentum_x2_global != 0.0 ||
       ps_escaped_injected_cr_momentum_x3_global != 0.0 ||
       ps_escaped_injected_cr_energy_global != 0.0);
  if (!std::isfinite(ps_particle_macro_mass) || ps_particle_macro_mass <= 0.0 ||
      !std::isfinite(ps_mass_reservoir_global) ||
      ps_mass_reservoir_global < 0.0 ||
      ps_mass_reservoir_global >= ps_particle_macro_mass ||
      invalid_count_ledger(ps_injected_cr_count_global) ||
      invalid_nonnegative_ledger(ps_injected_cr_mass_global) ||
      !std::isfinite(ps_injected_cr_momentum_x1_global) ||
      !std::isfinite(ps_injected_cr_momentum_x2_global) ||
      !std::isfinite(ps_injected_cr_momentum_x3_global) ||
      invalid_nonnegative_ledger(ps_injected_cr_energy_global) ||
      invalid_count_ledger(ps_removed_cr_count_global) ||
      invalid_nonnegative_ledger(ps_removed_cr_mass_global) ||
      !std::isfinite(ps_removed_cr_momentum_x1_global) ||
      !std::isfinite(ps_removed_cr_momentum_x2_global) ||
      !std::isfinite(ps_removed_cr_momentum_x3_global) ||
      invalid_nonnegative_ledger(ps_removed_cr_energy_global) ||
      !ps_escape_ledger_complete ||
      invalid_count_ledger(ps_escaped_injected_cr_count_global) ||
      invalid_nonnegative_ledger(ps_escaped_injected_cr_mass_global) ||
      !std::isfinite(ps_escaped_injected_cr_momentum_x1_global) ||
      !std::isfinite(ps_escaped_injected_cr_momentum_x2_global) ||
      !std::isfinite(ps_escaped_injected_cr_momentum_x3_global) ||
      invalid_nonnegative_ledger(ps_escaped_injected_cr_energy_global) ||
      invalid_count_ledger(ps_escaped_initial_cr_count_global) ||
      ps_escaped_initial_cr_count_global != 0.0 ||
      invalid_escape_before_audit || invalid_escape_audit_time ||
      invalid_empty_injected_ledger || invalid_empty_removed_ledger ||
      invalid_empty_escape_ledger ||
      invalid_tag_window ||
      !std::isfinite(expected_injected_mass) ||
      !std::isfinite(expected_removed_mass) ||
      !std::isfinite(expected_escaped_mass) ||
      !ParallelShockLedgerValuesAgree(ps_injected_cr_mass_global,
                                      expected_injected_mass) ||
      !ParallelShockLedgerValuesAgree(ps_removed_cr_mass_global,
                                      expected_removed_mass) ||
      !ParallelShockLedgerValuesAgree(ps_escaped_injected_cr_mass_global,
                                      expected_escaped_mass) ||
      ps_removed_cr_count_global + ps_escaped_injected_cr_count_global >
          ps_injected_cr_count_global ||
      ParallelShockLedgerValueExceeds(
          ps_removed_cr_mass_global + ps_escaped_injected_cr_mass_global,
                                      ps_injected_cr_mass_global) ||
      invalid_removed_before_sink || invalid_sink_completion) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock " << context
              << " CR ledger numeric metadata is invalid." << std::endl;
    restart_utils::AbortOnFatalError();
  }
}

void StoreRuntimeStateForRestart(const Real current_time) {
  if (ps_pin == nullptr) return;
  ValidateParallelShockRuntimeLedger("runtime", current_time);
  ps_pin->SetReal("problem", "ps_mass_reservoir_global", ps_mass_reservoir_global);
  ps_pin->SetReal("problem", "ps_injected_cr_count_global",
                  ps_injected_cr_count_global);
  ps_pin->SetReal("problem", "ps_injected_cr_mass_global",
                  ps_injected_cr_mass_global);
  ps_pin->SetReal("problem", "ps_injected_cr_momentum_x1_global",
                  ps_injected_cr_momentum_x1_global);
  ps_pin->SetReal("problem", "ps_injected_cr_momentum_x2_global",
                  ps_injected_cr_momentum_x2_global);
  ps_pin->SetReal("problem", "ps_injected_cr_momentum_x3_global",
                  ps_injected_cr_momentum_x3_global);
  ps_pin->SetReal("problem", "ps_injected_cr_energy_global",
                  ps_injected_cr_energy_global);
  ps_pin->SetBoolean("problem", "ps_removed_excluded_early_cohort",
                     ps_removed_excluded_early_cohort);
  ps_pin->SetReal("problem", "ps_removed_cr_count_global",
                  ps_removed_cr_count_global);
  ps_pin->SetReal("problem", "ps_removed_cr_mass_global",
                  ps_removed_cr_mass_global);
  ps_pin->SetReal("problem", "ps_removed_cr_momentum_x1_global",
                  ps_removed_cr_momentum_x1_global);
  ps_pin->SetReal("problem", "ps_removed_cr_momentum_x2_global",
                  ps_removed_cr_momentum_x2_global);
  ps_pin->SetReal("problem", "ps_removed_cr_momentum_x3_global",
                  ps_removed_cr_momentum_x3_global);
  ps_pin->SetReal("problem", "ps_removed_cr_energy_global",
                  ps_removed_cr_energy_global);
  ps_pin->SetInteger("problem", "ps_escape_ledger_schema", 1);
  ps_pin->SetBoolean("problem", "ps_escape_ledger_complete",
                     ps_escape_ledger_complete);
  ps_pin->SetInteger("problem", "ps_escape_audit_calls", ps_escape_audit_calls);
  ps_pin->SetReal("problem", "ps_escape_last_audit_time",
                  ps_escape_last_audit_time);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_count_global",
                  ps_escaped_injected_cr_count_global);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_mass_global",
                  ps_escaped_injected_cr_mass_global);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_momentum_x1_global",
                  ps_escaped_injected_cr_momentum_x1_global);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_momentum_x2_global",
                  ps_escaped_injected_cr_momentum_x2_global);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_momentum_x3_global",
                  ps_escaped_injected_cr_momentum_x3_global);
  ps_pin->SetReal("problem", "ps_escaped_injected_cr_energy_global",
                  ps_escaped_injected_cr_energy_global);
  ps_pin->SetReal("problem", "ps_escaped_initial_cr_count_global",
                  ps_escaped_initial_cr_count_global);
  ps_pin->SetInteger("problem", "ps_cr_ledger_schema", 3);
  ps_pin->SetBoolean("problem", "ps_cr_ledger_complete", ps_cr_ledger_complete);
  ps_pin->SetBoolean("problem", "ps_tag_seeded", ps_tag_seeded);
  if (ps_next_tag < 0 ||
      ps_next_tag > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock next particle tag cannot be stored in "
              << "restart metadata." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ps_pin->SetInteger("problem", "ps_injection_tag_floor",
                     static_cast<int>(ps_injection_tag_floor));
  ps_pin->SetInteger("problem", "ps_next_tag",
                     static_cast<int>(ps_next_tag));
}

void ObserveParallelShockParticleDestruction(particles::Particles *ppart, Mesh *pm,
                                             const int stage) {
  if (ppart == nullptr || pm == nullptr || ppart->pbval_part == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle-destruction observer received an "
              << "invalid runtime object." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  Real local[8] = {};
  Real audit_time = pm->time + pm->dt;
  if (!std::isfinite(pm->time) || !std::isfinite(pm->dt) || pm->dt <= 0.0) {
    local[7] += 1.0;
  } else if (ppart->UsesPaperVL2Coupling()) {
    if (stage == 1) {
      audit_time = pm->time + 0.5*pm->dt;
    } else if (stage == 2) {
      audit_time = pm->time + pm->dt;
    } else {
      local[7] += 1.0;
    }
  }
  if (!std::isfinite(audit_time) ||
      audit_time + ParallelShockLedgerTolerance(audit_time,
                                                 ps_escape_last_audit_time) <
          ps_escape_last_audit_time) {
    local[7] += 1.0;
  }
  if (!PaperVL2EscapeStageChronologyIsValid(ppart, pm, stage, audit_time)) {
    local[7] += 1.0;
  }

  const int npart = ppart->nprtcl_thispack;
  const int ndestroy = ppart->pbval_part->nprtcl_destroy;
  const auto &destroylist = ppart->pbval_part->destroylist.h_view;
  std::vector<int> destroyed_indices;
  if (npart < 0 || ndestroy < 0 || ndestroy > npart ||
      static_cast<std::size_t>(ndestroy) > destroylist.extent(0)) {
    local[7] += 1.0;
  } else {
    destroyed_indices.reserve(static_cast<std::size_t>(ndestroy));
    for (int n = 0; n < ndestroy; ++n) {
      const ParticleLocationData entry = destroylist(n);
      const int p = entry.prtcl_indx;
      if (p < 0 || p >= npart) {
        local[7] += 1.0;
        continue;
      }
      destroyed_indices.push_back(p);
      if (entry.destruction_reason !=
              static_cast<int>(ParticleDestructionReason::physical_boundary) ||
          entry.physical_boundary_mask != particle_boundary_outer_x1) {
        local[7] += 1.0;
      }
    }
    std::sort(destroyed_indices.begin(), destroyed_indices.end());
    if (std::adjacent_find(destroyed_indices.begin(), destroyed_indices.end()) !=
        destroyed_indices.end()) {
      local[7] += 1.0;
    }
  }
  if (local[7] == 0.0 && ndestroy > 0) {
    auto &pr = ppart->prtcl_rdata;
    auto &pi = ppart->prtcl_idata;
    auto &destroylist_d = ppart->pbval_part->destroylist;
    const int nrdata = ppart->nrdata;
    const int inject_species = ps_inject_species;
    const Real q_over_m = ps_particle_q_over_m;
    const Real macro_mass = ps_particle_macro_mass;
    const bool momentum_state = ppart->UsesRelativisticCRState();
    const Real light_speed = ppart->pic_cr_light_speed;
    Kokkos::parallel_reduce(
        "ps_accumulate_particle_escape_sink",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, ndestroy),
        KOKKOS_LAMBDA(const int n, Real &count, Real &mass, Real &momentum_x1,
                      Real &momentum_x2, Real &momentum_x3, Real &energy_sum,
                      Real &initial_count, Real &invalid) {
          const int p = destroylist_d.d_view(n).prtcl_indx;
          bool finite_payload = true;
          for (int q = 0; q < nrdata; ++q) {
            finite_payload = finite_payload && isfinite(pr(q, p));
          }
          if (!finite_payload) {
            invalid += 1.0;
            return;
          }
          const int source = pi(PCRSOURCE, p);
          if (source == static_cast<int>(CRParticleSource::initial)) {
            invalid += 1.0;
            return;
          }
          if (source != static_cast<int>(CRParticleSource::shock_injected) ||
              pi(PSP, p) != inject_species || pi(PTAG, p) < 0 ||
              pr(IPM, p) != q_over_m || pr(IPWT, p) != 1.0) {
            invalid += 1.0;
            return;
          }
          const Real state_x = pr(IPVX, p);
          const Real state_y = pr(IPVY, p);
          const Real state_z = pr(IPVZ, p);
          const Real energy = particles::CRKineticEnergy(
              momentum_state, light_speed, state_x, state_y, state_z);
          if (!isfinite(energy) || energy < 0.0) {
            invalid += 1.0;
            return;
          }
          count += 1.0;
          mass += macro_mass;
          momentum_x1 += macro_mass*state_x;
          momentum_x2 += macro_mass*state_y;
          momentum_x3 += macro_mass*state_z;
          energy_sum += macro_mass*energy;
        },
        Kokkos::Sum<Real>(local[0]),
        Kokkos::Sum<Real>(local[1]),
        Kokkos::Sum<Real>(local[2]),
        Kokkos::Sum<Real>(local[3]),
        Kokkos::Sum<Real>(local[4]),
        Kokkos::Sum<Real>(local[5]),
        Kokkos::Sum<Real>(local[6]),
        Kokkos::Sum<Real>(local[7]));
  }

#if MPI_PARALLEL_ENABLED
  Real global[8] = {};
  MPI_Allreduce(local, global, 8, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#else
  Real *global = local;
#endif
  if (global[7] != 0.0 ||
      !ParallelShockLedgerValuesAgree(global[1],
                                      global[0]*ps_particle_macro_mass,
                                      global[0])) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock found an unaccounted or invalid particle "
              << "destruction event." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_escape_audit_calls == std::numeric_limits<int>::max()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle escape audit-call counter overflow."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ++ps_escape_audit_calls;
  ps_escape_last_audit_time = audit_time;
  ps_escaped_injected_cr_count_global += global[0];
  ps_escaped_injected_cr_mass_global =
      ps_escaped_injected_cr_count_global*ps_particle_macro_mass;
  ps_escaped_injected_cr_momentum_x1_global += global[2];
  ps_escaped_injected_cr_momentum_x2_global += global[3];
  ps_escaped_injected_cr_momentum_x3_global += global[4];
  ps_escaped_injected_cr_energy_global += global[5];
  ps_escaped_initial_cr_count_global += global[6];
  StoreRuntimeStateForRestart(audit_time);
  if (global_variable::my_rank == 0 && (global[0] > 0.0 || global[6] > 0.0)) {
    std::cout << std::setprecision(17)
              << "pic_parallel_shock outer_x1_escape_sink: time=" << audit_time
              << " injected_count=" << ps_escaped_injected_cr_count_global
              << " injected_mass=" << ps_escaped_injected_cr_mass_global
              << " injected_momentum=("
              << ps_escaped_injected_cr_momentum_x1_global << ","
              << ps_escaped_injected_cr_momentum_x2_global << ","
              << ps_escaped_injected_cr_momentum_x3_global << ")"
              << " injected_energy=" << ps_escaped_injected_cr_energy_global
              << " initial_count=" << ps_escaped_initial_cr_count_global
              << std::endl;
  }
}

void ValidateParallelShockParticlePopulation(Mesh *pm) {
  if (pm == nullptr || pm->pmb_pack == nullptr || pm->pmb_pack->ppart == nullptr) {
    return;
  }
  if (ps_particle_population_audit_calls ==
      std::numeric_limits<std::int64_t>::max()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle-population audit counter overflow."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ++ps_particle_population_audit_calls;
  auto *ppart = pm->pmb_pack->ppart;
  auto &pi = ppart->prtcl_idata;
  auto &pr = ppart->prtcl_rdata;
  const int npart = ppart->nprtcl_thispack;
  const int inject_species = ps_inject_species;
  const Real q_over_m = ps_particle_q_over_m;
  Real active_injected_local = 0.0;
  Real invalid_local = 0.0;
  Kokkos::parallel_reduce(
      "ps_validate_active_injected_population",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, npart),
      KOKKOS_LAMBDA(const int p, Real &active_injected, Real &invalid) {
        const int source = pi(PCRSOURCE, p);
        if (source == static_cast<int>(CRParticleSource::initial)) return;
        if (source != static_cast<int>(CRParticleSource::shock_injected) ||
            pi(PSP, p) != inject_species || pi(PTAG, p) < 0 ||
            pr(IPM, p) != q_over_m || pr(IPWT, p) != 1.0) {
          invalid += 1.0;
          return;
        }
        active_injected += 1.0;
      },
      Kokkos::Sum<Real>(active_injected_local),
      Kokkos::Sum<Real>(invalid_local));
#if MPI_PARALLEL_ENABLED
  Real local[2] = {active_injected_local, invalid_local};
  Real global[2] = {};
  MPI_Allreduce(local, global, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  const Real active_injected_global = global[0];
  const Real invalid_global = global[1];
#else
  const Real active_injected_global = active_injected_local;
  const Real invalid_global = invalid_local;
#endif
  const Real accounted =
      active_injected_global + ps_removed_cr_count_global +
      ps_escaped_injected_cr_count_global;
  if (invalid_global != 0.0 ||
      !ParallelShockLedgerValuesAgree(accounted, ps_injected_cr_count_global,
                                      ps_injected_cr_count_global)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock active + startup-removed + escaped injected "
              << "particle count does not match the injection ledger." << std::endl;
    restart_utils::AbortOnFatalError();
  }
}

void MigrateParticlesAfterHostEdit(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;
  auto *ppart = pmbp->ppart;
  if (ppart->pbval_part == nullptr) {
    pm->CountParticles();
    return;
  }

  if (ppart->NewGID(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("NewGID did not complete");
  }
  if (ppart->SendCnt(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("SendCnt did not complete");
  }
  if (ppart->InitRecv(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("InitRecv did not complete");
  }
  if (ppart->SendP(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("SendP did not complete");
  }

  int nwait = 0;
  while (ppart->RecvP(nullptr, 0) != TaskStatus::complete) {
    ++nwait;
    if (nwait > 1000000) {
      FatalParticleMigrationError("RecvP did not complete");
    }
  }
  if (ppart->ClearRecv(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("ClearRecv did not complete");
  }
  if (ppart->ClearSend(nullptr, 0) != TaskStatus::complete) {
    FatalParticleMigrationError("ClearSend did not complete");
  }
  pm->CountParticles();
}

void ValidateLocalParticleOwnership(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;
  auto *ppart = pmbp->ppart;
  const int np = ppart->nprtcl_thispack;
  if (np <= 0) return;

  pmbp->pmb->mb_size.template sync<HostMemSpace>();
  auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
  auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_idata);
  const int gids = pmbp->gids;
  for (int p = 0; p < np; ++p) {
    const int gid = h_pi(PGID, p);
    const int m = gid - gids;
    if (m < 0 || m >= pmbp->nmb_thispack) {
      FatalParticleMigrationError("particle gid is not local after recenter");
    }
    const auto &size = pmbp->pmb->mb_size.h_view(m);
    const Real tol1 = 16.0*std::numeric_limits<Real>::epsilon()*
        std::max(static_cast<Real>(1.0), size.x1max - size.x1min);
    const Real tol2 = 16.0*std::numeric_limits<Real>::epsilon()*
        std::max(static_cast<Real>(1.0), size.x2max - size.x2min);
    const Real tol3 = 16.0*std::numeric_limits<Real>::epsilon()*
        std::max(static_cast<Real>(1.0), size.x3max - size.x3min);
    const Real x1 = h_pr(IPX, p);
    const Real x2 = h_pr(IPY, p);
    const Real x3 = h_pr(IPZ, p);
    if (x1 < size.x1min - tol1 || x1 > size.x1max + tol1 ||
        x2 < size.x2min - tol2 || x2 > size.x2max + tol2 ||
        x3 < size.x3min - tol3 || x3 > size.x3max + tol3) {
      FatalParticleMigrationError("particle position lies outside its local gid block");
    }
  }
}

void UpdateOuterInflowState(Mesh *pm, const Real frame_vx) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;
  auto &u_in = pmbp->pmhd->pbval_u->u_in;
  auto &b_in = pmbp->pmhd->pbval_b->b_in;
  const Real gamma = pmbp->pmhd->peos->eos_data.gamma;
  const Real gm1 = gamma - 1.0;
  const Real vin = -ps_u0 + frame_vx;
  const Real ein = ps_p0/gm1 + 0.5*ps_rho0*SQR(vin) + 0.5*SQR(ps_b0);

  u_in.h_view(IDN, BoundaryFace::outer_x1) = ps_rho0;
  u_in.h_view(IM1, BoundaryFace::outer_x1) = ps_rho0*vin;
  u_in.h_view(IM2, BoundaryFace::outer_x1) = 0.0;
  u_in.h_view(IM3, BoundaryFace::outer_x1) = 0.0;
  u_in.h_view(IEN, BoundaryFace::outer_x1) = ein;
  b_in.h_view(IBX, BoundaryFace::outer_x1) = ps_b0;
  b_in.h_view(IBY, BoundaryFace::outer_x1) = 0.0;
  b_in.h_view(IBZ, BoundaryFace::outer_x1) = 0.0;
  u_in.template modify<HostMemSpace>();
  u_in.template sync<DevExeSpace>();
  b_in.template modify<HostMemSpace>();
  b_in.template sync<DevExeSpace>();
}

void ApplyRecenteringShiftToParticles(Mesh *pm, const Real xshift) {
  if (xshift <= 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;

  auto *ppart = pmbp->ppart;
  const int np_old = ppart->nprtcl_thispack;
  if (np_old <= 0) return;

  const int nr = ppart->nrdata;
  const int ni = ppart->nidata;
  auto h_pr_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_rdata);
  auto h_pi_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_idata);
  HostArray2D<Real> h_pr_new("ps_pr_recent", nr, np_old);
  HostArray2D<int> h_pi_new("ps_pi_recent", ni, np_old);

  const Real xmin = pm->mesh_size.x1min;
  const Real xmax = pm->mesh_size.x1max;
  int np_new = 0;
  for (int p = 0; p < np_old; ++p) {
    const Real xnew = h_pr_old(IPX, p) - xshift;
    if (!(xnew > xmin && xnew < xmax)) continue;

    for (int q = 0; q < nr; ++q) h_pr_new(q, np_new) = h_pr_old(q, p);
    for (int q = 0; q < ni; ++q) h_pi_new(q, np_new) = h_pi_old(q, p);
    h_pr_new(IPX, np_new) = xnew;
    ++np_new;
  }

  Kokkos::resize(ppart->prtcl_rdata, nr, np_new);
  Kokkos::resize(ppart->prtcl_idata, ni, np_new);
  CopyPackedParticleDataToDevice(ppart, h_pi_new, h_pr_new, np_new, 0);
  ppart->nprtcl_thispack = np_new;
}

void RemoveExcludedEarlyInjectedParticles(Mesh *pm) {
  if (ps_removed_excluded_early_cohort || ps_remove_birth_time_before < 0.0 ||
      pm->time < ps_remove_birth_time_before) {
    return;
  }
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;

  auto *ppart = pmbp->ppart;
  const int np_old = ppart->nprtcl_thispack;
  const int nr = ppart->nrdata;
  const int ni = ppart->nidata;
  auto h_pr_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_rdata);
  auto h_pi_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_idata);
  HostArray2D<Real> h_pr_new("ps_pr_filtered", nr, np_old);
  HostArray2D<int> h_pi_new("ps_pi_filtered", ni, np_old);

  int np_new = 0;
  Real removed_local[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int p = 0; p < np_old; ++p) {
    const bool remove =
        h_pi_old(PCRSOURCE, p) == static_cast<int>(CRParticleSource::shock_injected) &&
        h_pr_old(IPT_BIRTH, p) < ps_remove_birth_time_before;
    if (remove) {
      const Real mass = ps_particle_macro_mass*h_pr_old(IPWT, p);
      removed_local[0] += 1.0;
      removed_local[1] += mass;
      removed_local[2] += mass*h_pr_old(IPVX, p);
      removed_local[3] += mass*h_pr_old(IPVY, p);
      removed_local[4] += mass*h_pr_old(IPVZ, p);
      removed_local[5] += mass*particles::CRKineticEnergy(
          ppart->UsesRelativisticCRState(), ppart->pic_cr_light_speed,
          h_pr_old(IPVX, p), h_pr_old(IPVY, p), h_pr_old(IPVZ, p));
      continue;
    }
    for (int q = 0; q < nr; ++q) h_pr_new(q, np_new) = h_pr_old(q, p);
    for (int q = 0; q < ni; ++q) h_pi_new(q, np_new) = h_pi_old(q, p);
    ++np_new;
  }

  Kokkos::resize(ppart->prtcl_rdata, nr, np_new);
  Kokkos::resize(ppart->prtcl_idata, ni, np_new);
  CopyPackedParticleDataToDevice(ppart, h_pi_new, h_pr_new, np_new, 0);
  ppart->nprtcl_thispack = np_new;
#if MPI_PARALLEL_ENABLED
  Real removed_global[6] = {};
  MPI_Allreduce(removed_local, removed_global, 6, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#else
  Real *removed_global = removed_local;
#endif
  const Real expected_removed_mass = removed_global[0]*ps_particle_macro_mass;
  if (!ParallelShockLedgerValuesAgree(removed_global[1], expected_removed_mass,
                                      removed_global[0])) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock removed-particle mass accounting is inconsistent."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ps_removed_cr_count_global += removed_global[0];
  ps_removed_cr_mass_global =
      ps_removed_cr_count_global*ps_particle_macro_mass;
  ps_removed_cr_momentum_x1_global += removed_global[2];
  ps_removed_cr_momentum_x2_global += removed_global[3];
  ps_removed_cr_momentum_x3_global += removed_global[4];
  ps_removed_cr_energy_global += removed_global[5];
  ps_removed_excluded_early_cohort = true;
  StoreRuntimeStateForRestart(pm->time);
  pm->CountParticles();
  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(17)
              << "pic_parallel_shock removed_cr_sink: count="
              << ps_removed_cr_count_global
              << " mass=" << ps_removed_cr_mass_global
              << " momentum=(" << ps_removed_cr_momentum_x1_global << ","
              << ps_removed_cr_momentum_x2_global << ","
              << ps_removed_cr_momentum_x3_global << ")"
              << " energy=" << ps_removed_cr_energy_global << std::endl;
  }
}

void ApplyRecenteringShift(Mesh *pm, const int nshift) {
  if (nshift <= 0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ng = indcs.ng;
  if (nshift > ng) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock recenter shift must be <= nghost."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  auto *pmhd = pmbp->pmhd;
  auto h_u0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->u0);
  auto h_x1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x1f);
  auto h_x2f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x2f);
  auto h_x3f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x3f);
  const int nvar = static_cast<int>(pmhd->u0.extent_int(1));
  const Real vin = -ps_u0 + FrameVelocityOffset(pm->time + pm->dt);
  const Real gm1 = pmhd->peos->eos_data.gamma - 1.0;
  const Real ein = ps_p0/gm1 + 0.5*ps_rho0*SQR(vin) + 0.5*SQR(ps_b0);

  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          if (src <= ie + ng) {
            for (int n = 0; n < nvar; ++n) {
              h_u0(m, n, k, j, i) = h_u0(m, n, k, j, src);
            }
          } else {
            for (int n = 0; n < nvar; ++n) {
              h_u0(m, n, k, j, i) = h_u0(m, n, k, j, ie);
            }
            h_u0(m, IDN, k, j, i) = ps_rho0;
            h_u0(m, IM1, k, j, i) = ps_rho0*vin;
            h_u0(m, IM2, k, j, i) = 0.0;
            h_u0(m, IM3, k, j, i) = 0.0;
            h_u0(m, IEN, k, j, i) = ein;
          }
        }
      }
    }
  }

  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie + 1; ++i) {
          const int src = i + nshift;
          h_x1f(m, k, j, i) = (src <= ie + 1 + ng) ?
              h_x1f(m, k, j, src) : ps_b0;
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je + 1; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          h_x2f(m, k, j, i) = (src <= ie + ng) ? h_x2f(m, k, j, src) : 0.0;
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          h_x3f(m, k, j, i) = (src <= ie + ng) ? h_x3f(m, k, j, src) : 0.0;
        }
      }
    }
  }

  Kokkos::deep_copy(pmhd->u0, h_u0);
  Kokkos::deep_copy(pmhd->b0.x1f, h_x1f);
  Kokkos::deep_copy(pmhd->b0.x2f, h_x2f);
  Kokkos::deep_copy(pmhd->b0.x3f, h_x3f);

  auto &u0 = pmhd->u0;
  auto &b0 = pmhd->b0;
  auto &w0 = pmhd->w0;
  auto &bcc0 = pmhd->bcc0;
  pmhd->peos->ConsToPrim(u0, b0, w0, bcc0, false, is, ie, js, je, ks, ke);

  if (ps_frame_apply_to_particles) {
    const Real xshift = static_cast<Real>(nshift)*ps_recenter_dx1;
    ApplyRecenteringShiftToParticles(pm, xshift);
    MigrateParticlesAfterHostEdit(pm);
    ValidateLocalParticleOwnership(pm);
  }
}

inline Real TaggedUniform01(const std::int64_t tag, const std::uint64_t stream) {
  const std::uint64_t tag_bits = static_cast<std::uint64_t>(tag);
  const std::uint64_t seed_bits = static_cast<std::uint64_t>(ps_inject_seed);
  return UniformFromUint64(SplitMix64(
      seed_bits ^ (tag_bits + 1ULL)*0x9e3779b97f4a7c15ULL ^
      (stream + 1ULL)*0xbf58476d1ce4e5b9ULL));
}

inline Real ClampInsideDomain(const Real x, const Real xmin, const Real xmax) {
  const Real span = xmax - xmin;
  const Real eps = std::max(static_cast<Real>(1.0e-12)*span,
                            static_cast<Real>(1.0e-14));
  return std::min(std::max(x, xmin + eps), xmax - eps);
}

void SeedNextTag(particles::Particles *ppart, const Real current_time) {
  if (ps_tag_seeded && ps_tag_progression_validated) return;

  int local_max = -1;
  int local_baseline_max = -1;
  bool local_payload_invalid = false;
  std::vector<int> local_tags;
  int npart = ppart->nprtcl_thispack;
  local_tags.reserve(npart);
  if (npart > 0) {
    auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
    auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_idata);
    for (int p = 0; p < npart; ++p) {
      const int tag = h_pi(PTAG, p);
      const int source = h_pi(PCRSOURCE, p);
      local_payload_invalid = local_payload_invalid || tag < 0;
      local_tags.push_back(tag);
      local_max = std::max(local_max, tag);
      for (int q = 0; q < ppart->nrdata; ++q) {
        local_payload_invalid = local_payload_invalid || !std::isfinite(h_pr(q, p));
      }
      if (source == static_cast<int>(CRParticleSource::initial)) {
        local_baseline_max = std::max(local_baseline_max, tag);
        local_payload_invalid =
            local_payload_invalid || (ps_tag_seeded && tag >= ps_injection_tag_floor);
      } else if (source == static_cast<int>(CRParticleSource::shock_injected)) {
        local_payload_invalid =
            local_payload_invalid || !ps_tag_seeded ||
            tag < ps_injection_tag_floor || tag >= ps_next_tag ||
            h_pi(PSP, p) != ps_inject_species ||
            h_pr(IPM, p) != ps_particle_q_over_m || h_pr(IPWT, p) != 1.0 ||
            h_pr(IPT_BIRTH, p) < ps_inject_t_start ||
            h_pr(IPT_BIRTH, p) > ps_inject_t_stop ||
            h_pr(IPT_BIRTH, p) > current_time ||
            (ps_removed_excluded_early_cohort &&
             h_pr(IPT_BIRTH, p) < ps_remove_birth_time_before);
      } else {
        local_payload_invalid = true;
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  int payload_invalid = local_payload_invalid ? 1 : 0;
  int global_payload_invalid = 0;
  MPI_Allreduce(&payload_invalid, &global_payload_invalid, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (global_payload_invalid != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle provenance payload is invalid."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  const int nranks = global_variable::nranks;
  int global_max = -1;
  int global_baseline_max = -1;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_baseline_max, &global_baseline_max, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  std::vector<int> send_counts(nranks, 0);
  std::vector<int> recv_counts(nranks, 0);
  for (const int tag : local_tags) ++send_counts[tag % nranks];
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  std::vector<int> send_displs(nranks, 0);
  std::vector<int> recv_displs(nranks, 0);
  int send_total = 0;
  int recv_total = 0;
  for (int rank = 0; rank < nranks; ++rank) {
    if (send_counts[rank] > std::numeric_limits<int>::max() - send_total ||
        recv_counts[rank] > std::numeric_limits<int>::max() - recv_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock particle tag audit exceeds MPI count limits."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    send_displs[rank] = send_total;
    recv_displs[rank] = recv_total;
    send_total += send_counts[rank];
    recv_total += recv_counts[rank];
  }
  const std::int64_t audited_tag_count = ps_tag_seeded ? ps_next_tag :
      static_cast<std::int64_t>(global_max) + 1;
  const std::int64_t max_owner_tags =
      (audited_tag_count + static_cast<std::int64_t>(nranks) - 1) / nranks;
  int local_owner_overflow = recv_total > max_owner_tags ? 1 : 0;
  int global_owner_overflow = 0;
  MPI_Allreduce(&local_owner_overflow, &global_owner_overflow, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (global_owner_overflow != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle tag audit exceeds owner allocation bounds."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  std::vector<int> send_tags(send_total, 0);
  std::vector<int> recv_tags(recv_total, 0);
  std::vector<int> send_cursor = send_displs;
  for (const int tag : local_tags) {
    const int owner = tag % nranks;
    send_tags[send_cursor[owner]++] = tag;
  }
  MPI_Alltoallv(send_tags.data(), send_counts.data(), send_displs.data(), MPI_INT,
                recv_tags.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                MPI_COMM_WORLD);
  std::sort(recv_tags.begin(), recv_tags.end());
  int local_duplicate =
      std::adjacent_find(recv_tags.begin(), recv_tags.end()) != recv_tags.end() ? 1 : 0;
  int global_duplicate = 0;
  MPI_Allreduce(&local_duplicate, &global_duplicate, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (global_duplicate != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle tags are not globally unique."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
#else
  if (local_payload_invalid) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle provenance payload is invalid."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  std::sort(local_tags.begin(), local_tags.end());
  if (std::adjacent_find(local_tags.begin(), local_tags.end()) != local_tags.end()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock particle tags are not globally unique."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  int global_max = local_max;
  int global_baseline_max = local_baseline_max;
#endif

  const std::int64_t start = static_cast<std::int64_t>(global_max) + 1;
  if (start > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle tag range exhausted before pic_parallel_shock injection."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_tag_seeded) {
    const std::int64_t baseline_start =
        static_cast<std::int64_t>(global_baseline_max) + 1;
    if (ps_next_tag < start || ps_injection_tag_floor < baseline_start) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock persisted CR tag progression is invalid."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ps_tag_progression_validated = true;
    return;
  }
  ps_injection_tag_floor = start;
  ps_next_tag = start;
  ps_tag_seeded = true;
  ps_tag_progression_validated = true;
  StoreRuntimeStateForRestart(current_time);
}

void ApplyParallelShockGasSubtraction(Mesh *pm, const Real stage_weight) {
  if (!ps_enable_subtraction) return;
  if (!(stage_weight > 0.0) || !std::isfinite(stage_weight)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock gas-subtraction stage weight must be finite "
              << "and positive." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;
  auto *pmhd = pmbp->pmhd;
  const int nsub = static_cast<int>(ps_injection_transaction_gas_deltas.size());
  HostArray1D<int> h_m("ps_m", nsub);
  HostArray1D<int> h_k("ps_k", nsub);
  HostArray1D<int> h_j("ps_j", nsub);
  HostArray1D<int> h_i("ps_i", nsub);
  HostArray1D<Real> h_dm("ps_dm", nsub);
  HostArray1D<Real> h_dmx("ps_dmx", nsub);
  HostArray1D<Real> h_dmy("ps_dmy", nsub);
  HostArray1D<Real> h_dmz("ps_dmz", nsub);
  HostArray1D<Real> h_de("ps_de", nsub);
  for (int n = 0; n < nsub; ++n) {
    const GasDelta &d = ps_injection_transaction_gas_deltas[n];
    h_m(n) = d.m;
    h_k(n) = d.k;
    h_j(n) = d.j;
    h_i(n) = d.i;
    h_dm(n) = stage_weight*d.dm;
    h_dmx(n) = stage_weight*d.dmx;
    h_dmy(n) = stage_weight*d.dmy;
    h_dmz(n) = stage_weight*d.dmz;
    h_de(n) = stage_weight*d.de;
  }

  auto d_m = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_m);
  auto d_k = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_k);
  auto d_j = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_j);
  auto d_i = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_i);
  auto d_dm = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dm);
  auto d_dmx = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmx);
  auto d_dmy = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmy);
  auto d_dmz = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmz);
  auto d_de = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_de);

  auto &u0 = pmhd->u0;
  auto &b0 = pmhd->b0;
  const Real rho_floor = ps_rho_floor_frac*ps_rho0;
  const Real p_floor = ps_p_floor_frac*ps_p0;
  const Real gm1 = pmhd->peos->eos_data.gamma - 1.0;
  int local_density_floor_clips = 0;
  int local_pressure_floor_clips = 0;
  int local_nonfinite_cells = 0;
  Kokkos::parallel_reduce(
    "ps_gas_subtract_validate", Kokkos::RangePolicy<>(DevExeSpace(), 0, nsub),
  KOKKOS_LAMBDA(const int n, int &density_clips, int &pressure_clips,
                int &nonfinite_cells) {
    const int m = d_m(n);
    const int k = d_k(n);
    const int j = d_j(n);
    const int i = d_i(n);
    const Real dm = d_dm(n);
    if (!isfinite(dm) || !isfinite(d_dmx(n)) || !isfinite(d_dmy(n)) ||
        !isfinite(d_dmz(n)) || !isfinite(d_de(n))) {
      ++nonfinite_cells;
      return;
    }
    if (dm <= 0.0) return;
    const Real rho = u0(m, IDN, k, j, i) - dm;
    if (!isfinite(rho)) {
      ++nonfinite_cells;
      return;
    }
    if (rho < rho_floor) {
      ++density_clips;
      return;
    }
    const Real mx = u0(m, IM1, k, j, i) - d_dmx(n);
    const Real my = u0(m, IM2, k, j, i) - d_dmy(n);
    const Real mz = u0(m, IM3, k, j, i) - d_dmz(n);
    const Real energy = u0(m, IEN, k, j, i) - d_de(n);
    const Real bx = 0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    const Real by = 0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    const Real bz = 0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
    const Real kin = 0.5*(SQR(mx) + SQR(my) + SQR(mz))/rho;
    const Real efloor = p_floor/gm1 + kin + 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    if (!isfinite(mx) || !isfinite(my) || !isfinite(mz) || !isfinite(energy) ||
        !isfinite(bx) || !isfinite(by) || !isfinite(bz) || !isfinite(efloor)) {
      ++nonfinite_cells;
      return;
    }
    if (energy < efloor) ++pressure_clips;
  }, Kokkos::Sum<int>(local_density_floor_clips),
     Kokkos::Sum<int>(local_pressure_floor_clips),
     Kokkos::Sum<int>(local_nonfinite_cells));
#if MPI_PARALLEL_ENABLED
  int global_density_floor_clips = 0;
  int global_pressure_floor_clips = 0;
  int global_nonfinite_cells = 0;
  MPI_Allreduce(&local_density_floor_clips, &global_density_floor_clips, 1,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_pressure_floor_clips, &global_pressure_floor_clips, 1,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_nonfinite_cells, &global_nonfinite_cells, 1,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  int global_density_floor_clips = local_density_floor_clips;
  int global_pressure_floor_clips = local_pressure_floor_clips;
  int global_nonfinite_cells = local_nonfinite_cells;
#endif
  if (global_density_floor_clips > 0 || global_pressure_floor_clips > 0 ||
      global_nonfinite_cells > 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock gas subtraction would violate a fluid floor: "
              << "density_floor_cells=" << global_density_floor_clips
              << " pressure_floor_cells=" << global_pressure_floor_clips
              << " nonfinite_cells=" << global_nonfinite_cells
              << ". The process is stopping before clipping or checkpoint publication."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (nsub <= 0) return;
  par_for("ps_gas_subtract", DevExeSpace(), 0, nsub - 1,
  KOKKOS_LAMBDA(const int n) {
    const int m = d_m(n);
    const int k = d_k(n);
    const int j = d_j(n);
    const int i = d_i(n);
    const Real dm = d_dm(n);
    if (dm <= 0.0) return;

    u0(m, IDN, k, j, i) -= dm;
    u0(m, IM1, k, j, i) -= d_dmx(n);
    u0(m, IM2, k, j, i) -= d_dmy(n);
    u0(m, IM3, k, j, i) -= d_dmz(n);
    u0(m, IEN, k, j, i) -= d_de(n);
  });
  std::array<Real, 5> stage_delta = {};
  for (const GasDelta &d : ps_injection_transaction_gas_deltas) {
    stage_delta[0] += d.dm*d.vol;
    stage_delta[1] += d.dmx*d.vol;
    stage_delta[2] += d.dmy*d.vol;
    stage_delta[3] += d.dmz*d.vol;
    stage_delta[4] += d.de*d.vol;
  }
  const bool paper_vl2 =
      (pm->pmb_pack != nullptr) && (pm->pmb_pack->ppart != nullptr) &&
      pm->pmb_pack->ppart->UsesPaperVL2Coupling();
  for (int n=0; n<5; ++n) {
    if (paper_vl2) {
      // VL2 stage 2 resets to the saved cycle-start state before replaying the
      // full source transaction. The stage-1 half-step is predictor-only.
      ps_injection_transaction_applied_local[n] = stage_weight*stage_delta[n];
    } else {
      // Qualified legacy shock injection uses SSPRK1/2/3, whose source-only
      // low-storage recurrence is S <- beta*(S + Delta).
      ps_injection_transaction_applied_local[n] =
          stage_weight*(ps_injection_transaction_applied_local[n] + stage_delta[n]);
    }
  }
}

void PrepareParallelShockInjectionTransaction(Mesh *pm) {
  if (ps_injection_transaction_cycle == pm->ncycle) return;
  ps_injection_transaction_cycle = pm->ncycle;
  ps_injection_transaction_gas_deltas.clear();
  ps_injection_transaction_expected_global.fill(0.0);
  ps_injection_transaction_applied_local.fill(0.0);
  ps_injection_transaction_abs_global.fill(0.0);
  ps_injection_transaction_terms_global = 1.0;
  if (!ps_enable_injection) return;
  if (pm->time < ps_inject_t_start || pm->time > ps_inject_t_stop) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr || pmbp->ppart == nullptr) return;
  if (!(pm->dt > 0.0) || !std::isfinite(pm->dt)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injection timestep must be finite and positive."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  // Particle creation and reservoir consumption are irreversible.  Commit them
  // once per physical cycle, then replay only the matching RK-weighted fluid
  // subtraction on later stages.
  auto *ppart = pmbp->ppart;
  SeedNextTag(ppart, pm->time);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const bool three_d = pm->three_d;

  auto &mb_size = pmbp->pmb->mb_size;
  auto &mb_gid = pmbp->pmb->mb_gid;
  mb_size.template sync<HostMemSpace>();
  mb_gid.template sync<HostMemSpace>();

  const Real xshock = ShockSurfaceModelX1(pm->time);
  std::vector<ShockCell> cells;
  cells.reserve(static_cast<std::size_t>(pmbp->nmb_thispack)*indcs.nx2*indcs.nx3);

  Real area_total = 0.0;
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    const Real x1min = mb_size.h_view(m).x1min;
    const Real x2min = mb_size.h_view(m).x2min;
    const Real x3min = mb_size.h_view(m).x3min;
    const Real dx1 = mb_size.h_view(m).dx1;
    const Real dx2 = mb_size.h_view(m).dx2;
    const Real dx3 = mb_size.h_view(m).dx3;
    const Real half_width = ps_inject_half_width_cells*dx1;

    for (int k = ks; k <= ke; ++k) {
      Real x3c = x3min + (static_cast<Real>(k - ks) + 0.5)*dx3;
      if (!three_d) x3c = 0.0;
      for (int j = js; j <= je; ++j) {
        const Real x2c = x2min + (static_cast<Real>(j - js) + 0.5)*dx2;
        for (int i = is; i <= ie; ++i) {
          const Real x1c = x1min + (static_cast<Real>(i - is) + 0.5)*dx1;
          if (!(xshock >= x1c - half_width && xshock < x1c + half_width)) continue;

          ShockCell cell;
          cell.m = m;
          cell.k = k;
          cell.j = j;
          cell.i = i;
          cell.gid = mb_gid.h_view(m);
          cell.x1c = x1c;
          cell.x2c = x2c;
          cell.x3c = x3c;
          cell.dx1 = dx1;
          cell.dx2 = dx2;
          cell.dx3 = dx3;
          cell.area = dx2*dx3;
          cell.vol = dx1*dx2*dx3;
          area_total += cell.area;
          cells.push_back(cell);
        }
      }
    }
  }

  Real reduced_area_global = area_total;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&area_total, &reduced_area_global, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  std::vector<GlobalShockCell> local_global_cells;
  local_global_cells.reserve(cells.size());
  for (std::size_t n = 0; n < cells.size(); ++n) {
    GlobalShockCell cell;
    cell.owner_rank = global_variable::my_rank;
    cell.local_cell_index = static_cast<int>(n);
    cell.gid = cells[n].gid;
    cell.x1c = cells[n].x1c;
    cell.x2c = cells[n].x2c;
    cell.x3c = cells[n].x3c;
    cell.dx1 = cells[n].dx1;
    cell.dx2 = cells[n].dx2;
    cell.dx3 = cells[n].dx3;
    cell.area = cells[n].area;
    local_global_cells.push_back(cell);
  }
  std::vector<GlobalShockCell> global_cells;
#if MPI_PARALLEL_ENABLED
  const int local_cell_bytes =
      static_cast<int>(local_global_cells.size()*sizeof(GlobalShockCell));
  std::vector<int> cell_bytes_eachrank(global_variable::nranks, 0);
  MPI_Allgather(&local_cell_bytes, 1, MPI_INT, cell_bytes_eachrank.data(), 1,
                MPI_INT, MPI_COMM_WORLD);
  std::vector<int> cell_byte_offsets(global_variable::nranks, 0);
  int global_cell_bytes = 0;
  for (int r = 0; r < global_variable::nranks; ++r) {
    cell_byte_offsets[r] = global_cell_bytes;
    global_cell_bytes += cell_bytes_eachrank[r];
  }
  global_cells.resize(static_cast<std::size_t>(global_cell_bytes)/
                      sizeof(GlobalShockCell));
  MPI_Allgatherv(local_global_cells.data(), local_cell_bytes, MPI_BYTE,
                 global_cells.data(), cell_bytes_eachrank.data(),
                 cell_byte_offsets.data(), MPI_BYTE, MPI_COMM_WORLD);
#else
  global_cells = local_global_cells;
#endif
  std::sort(global_cells.begin(), global_cells.end(),
            [](const GlobalShockCell &lhs, const GlobalShockCell &rhs) {
    return std::tie(lhs.x3c, lhs.x2c, lhs.x1c, lhs.dx3, lhs.dx2, lhs.dx1,
                    lhs.gid, lhs.owner_rank, lhs.local_cell_index) <
           std::tie(rhs.x3c, rhs.x2c, rhs.x1c, rhs.dx3, rhs.dx2, rhs.dx1,
                    rhs.gid, rhs.owner_rank, rhs.local_cell_index);
  });
  std::vector<Real> global_area_prefix(global_cells.size(), 0.0);
  Real global_running_area = 0.0;
  for (std::size_t n = 0; n < global_cells.size(); ++n) {
    global_running_area += global_cells[n].area;
    global_area_prefix[n] = global_running_area;
  }
  const Real area_tolerance =
      64.0*std::numeric_limits<Real>::epsilon()*
      std::max(static_cast<Real>(1.0), std::abs(reduced_area_global));
  if (std::abs(global_running_area - reduced_area_global) > area_tolerance) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock gathered shock-surface area does not match "
              << "the global reduction." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  int ninj_global = 0;
  Real reservoir_after = ps_mass_reservoir_global;
  if (global_running_area > 0.0) {
    const Real sweep_speed = ps_u0 + ps_shock_speed;
    const Real swept_mass = ps_eta*ps_rho0*sweep_speed*pm->dt*global_running_area;
    const Real mass_budget = ps_mass_reservoir_global + swept_mass;
    if (mass_budget > 0.0) {
      const Real ninj_real = std::floor(mass_budget/ps_particle_macro_mass);
      if (!std::isfinite(ninj_real) || ninj_real >
          static_cast<Real>(std::numeric_limits<int>::max())) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "pic_parallel_shock injected particle count exceeds int range."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
      ninj_global = static_cast<int>(ninj_real);
      reservoir_after = mass_budget -
          static_cast<Real>(ninj_global)*ps_particle_macro_mass;
    }
  }
  if (ninj_global <= 0) {
    ps_mass_reservoir_global = reservoir_after;
    StoreRuntimeStateForRestart(pm->time);
    return;
  }

  const std::int64_t tag_base = ps_next_tag;
  if (tag_base >
      static_cast<std::int64_t>(std::numeric_limits<int>::max()) -
      static_cast<std::int64_t>(ninj_global)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle tag range exhausted during pic_parallel_shock injection."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  std::map<std::tuple<int, int, int, int>, GasDelta> gas_deltas;
  std::vector<InjectedParticle> injected;
  injected.reserve(static_cast<std::size_t>(
      std::ceil(static_cast<Real>(ninj_global)*area_total/global_running_area)));
  Real injected_local[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Real injected_abs_local[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  const Real frame_vx = FrameVelocityOffset(pm->time);
  const Real surface_vx = ps_shock_speed + frame_vx;
  const Real pinj = ps_vinj_over_u0*ps_u0;
  const Real vinj = VelocityMagnitudeFromMomentumMagnitude(ppart, pinj);
  const auto &mesh_size = pm->mesh_size;
  const Real x1min = mesh_size.x1min;
  const Real x1max = mesh_size.x1max;
  const Real x2min = mesh_size.x2min;
  const Real x2max = mesh_size.x2max;
  const Real x3min = mesh_size.x3min;
  const Real x3max = mesh_size.x3max;
  for (int n = 0; n < ninj_global; ++n) {
    const std::int64_t tag = tag_base + static_cast<std::int64_t>(n);
    const Real draw = TaggedUniform01(tag, 0)*global_running_area;
    auto it = std::lower_bound(global_area_prefix.begin(), global_area_prefix.end(),
                               draw);
    std::size_t idx =
        static_cast<std::size_t>(std::distance(global_area_prefix.begin(), it));
    if (idx >= global_cells.size()) idx = global_cells.size() - 1;
    const GlobalShockCell &global_cell = global_cells[idx];
    if (global_cell.owner_rank != global_variable::my_rank) continue;
    const ShockCell &cell = cells[global_cell.local_cell_index];

    Real dirx = 0.0;
    Real diry = 0.0;
    Real dirz = 0.0;
    if (three_d || ps_use_2d3v) {
      // Section 5.4 requires isotropy relative to the ideal shock surface.
      const Real mu = 2.0*TaggedUniform01(tag, 1) - 1.0;
      const Real phi = 2.0*M_PI*TaggedUniform01(tag, 2);
      const Real st = std::sqrt(std::max(static_cast<Real>(0.0), 1.0 - mu*mu));
      dirx = mu;
      diry = st*std::cos(phi);
      dirz = st*std::sin(phi);
    } else {
      const Real phi = 2.0*M_PI*TaggedUniform01(tag, 1);
      dirx = std::cos(phi);
      diry = std::sin(phi);
      dirz = 0.0;
    }

    InjectedParticle part;
    part.tag = tag;
    part.gid = cell.gid;
    part.m = cell.m;
    part.k = cell.k;
    part.j = cell.j;
    part.i = cell.i;
    part.vol = cell.vol;
    part.x1 = xshock;
    part.x2 = cell.x2c + (TaggedUniform01(tag, 3) - 0.5)*cell.dx2;
    part.x3 = three_d ?
        (cell.x3c + (TaggedUniform01(tag, 4) - 0.5)*cell.dx3) : 0.0;
    part.x1 = ClampInsideDomain(part.x1, x1min, x1max);
    part.x2 = ClampInsideDomain(part.x2, x2min, x2max);
    if (three_d) part.x3 = ClampInsideDomain(part.x3, x3min, x3max);
    BoostRelativeVelocityFromSurface(ppart, surface_vx, vinj*dirx, vinj*diry,
                                     vinj*dirz, part.vx, part.vy, part.vz);
    injected.push_back(part);

    Real state_x, state_y, state_z;
    EncodeCRStateFromVelocity(ppart, part.vx, part.vy, part.vz,
                              state_x, state_y, state_z);
    const Real particle_energy = particles::CRKineticEnergy(
        ppart->UsesRelativisticCRState(), ppart->pic_cr_light_speed,
        state_x, state_y, state_z);
    injected_local[0] += 1.0;
    injected_local[1] += ps_particle_macro_mass;
    injected_local[2] += ps_particle_macro_mass*state_x;
    injected_local[3] += ps_particle_macro_mass*state_y;
    injected_local[4] += ps_particle_macro_mass*state_z;
    injected_local[5] += ps_particle_macro_mass*particle_energy;
    injected_abs_local[0] += ps_particle_macro_mass;
    injected_abs_local[1] += std::abs(ps_particle_macro_mass*state_x);
    injected_abs_local[2] += std::abs(ps_particle_macro_mass*state_y);
    injected_abs_local[3] += std::abs(ps_particle_macro_mass*state_z);
    injected_abs_local[4] += std::abs(ps_particle_macro_mass*particle_energy);
    if (!ps_enable_subtraction) continue;
    const auto key = std::make_tuple(cell.m, cell.k, cell.j, cell.i);
    auto dit = gas_deltas.find(key);
    if (dit == gas_deltas.end()) {
      GasDelta delta;
      delta.m = cell.m;
      delta.k = cell.k;
      delta.j = cell.j;
      delta.i = cell.i;
      delta.vol = cell.vol;
      delta.dm = 0.0;
      delta.dmx = 0.0;
      delta.dmy = 0.0;
      delta.dmz = 0.0;
      delta.de = 0.0;
      dit = gas_deltas.emplace(key, delta).first;
    }
    const Real inv_vol = 1.0/part.vol;
    const Real mass_rho = ps_particle_macro_mass*inv_vol;
    dit->second.dm += mass_rho;
    dit->second.dmx += mass_rho*state_x;
    dit->second.dmy += mass_rho*state_y;
    dit->second.dmz += mass_rho*state_z;
    dit->second.de += mass_rho*particle_energy;
  }

  const int old_npart = ppart->nprtcl_thispack;
  const int ninject = static_cast<int>(injected.size());
  const int new_npart = old_npart + ninject;
  if (ps_enable_subtraction) {
    for (const auto &kv : gas_deltas) {
      ps_injection_transaction_gas_deltas.push_back(kv.second);
    }
  }
#if MPI_PARALLEL_ENABLED
  Real injected_global[6] = {};
  Real injected_abs_global[5] = {};
  MPI_Allreduce(injected_local, injected_global, 6, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(injected_abs_local, injected_abs_global, 5, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#else
  Real *injected_global = injected_local;
  Real *injected_abs_global = injected_abs_local;
#endif
  if (std::abs(injected_global[0] - static_cast<Real>(ninj_global)) > 0.5) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock global injected-particle accounting "
              << "does not match the budget." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const Real expected_injected_mass = injected_global[0]*ps_particle_macro_mass;
  if (!ParallelShockLedgerValuesAgree(injected_global[1], expected_injected_mass,
                                      injected_global[0])) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injected-particle mass accounting is inconsistent."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_enable_subtraction) {
    ps_injection_transaction_terms_global =
        std::max(static_cast<Real>(1.0), injected_global[0]);
    for (int n=0; n<5; ++n) {
      ps_injection_transaction_expected_global[n] = injected_global[n + 1];
      ps_injection_transaction_abs_global[n] = injected_abs_global[n];
    }
  }
  ps_mass_reservoir_global = reservoir_after;
  ps_next_tag = tag_base + static_cast<std::int64_t>(ninj_global);
  ps_injected_cr_count_global += injected_global[0];
  ps_injected_cr_mass_global =
      ps_injected_cr_count_global*ps_particle_macro_mass;
  ps_injected_cr_momentum_x1_global += injected_global[2];
  ps_injected_cr_momentum_x2_global += injected_global[3];
  ps_injected_cr_momentum_x3_global += injected_global[4];
  ps_injected_cr_energy_global += injected_global[5];
  StoreRuntimeStateForRestart(pm->time);
  if (ninject <= 0) {
    pm->CountParticles();
    return;
  }

  HostArray2D<int> h_pi_new("ps_pi_new", ppart->nidata, ninject);
  HostArray2D<Real> h_pr_new("ps_pr_new", ppart->nrdata, ninject);
  for (int n = 0; n < ninject; ++n) {
    for (int q = 0; q < ppart->nidata; ++q) h_pi_new(q, n) = 0;
    for (int q = 0; q < ppart->nrdata; ++q) h_pr_new(q, n) = 0.0;
  }

  for (int n = 0; n < ninject; ++n) {
    h_pi_new(PGID, n) = injected[n].gid;
    h_pi_new(PSP, n) = ps_inject_species;
    h_pi_new(PTAG, n) = static_cast<int>(injected[n].tag);
    h_pi_new(PCRSOURCE, n) = static_cast<int>(CRParticleSource::shock_injected);

    h_pr_new(IPX, n) = injected[n].x1;
    h_pr_new(IPY, n) = injected[n].x2;
    h_pr_new(IPZ, n) = injected[n].x3;
    EncodeCRStateFromVelocity(ppart, injected[n].vx, injected[n].vy, injected[n].vz,
                              h_pr_new(IPVX, n), h_pr_new(IPVY, n), h_pr_new(IPVZ, n));
    h_pr_new(IPM, n) = ps_particle_q_over_m;
    h_pr_new(IPWT, n) = 1.0;
    const Real a1 = particles::PICScaleFactor(
        ppart->pic_expansion_law, ppart->pic_expansion_rate_x1, pm->time);
    const Real a2 = particles::PICScaleFactor(
        ppart->pic_expansion_law, ppart->pic_expansion_rate_x2, pm->time);
    const Real a3 = particles::PICScaleFactor(
        ppart->pic_expansion_law, ppart->pic_expansion_rate_x3, pm->time);
    h_pr_new(IPF0, n) = particles::PICDeltaFBackgroundValue(
        ppart->pic_deltaf_background, ppart->pic_deltaf_p0, ppart->pic_deltaf_kappa,
        ppart->pic_deltaf_drift_x1, ppart->pic_deltaf_drift_x2,
        ppart->pic_deltaf_drift_x3, ppart->pic_deltaf_aniso_x1,
        ppart->pic_deltaf_aniso_x2, ppart->pic_deltaf_aniso_x3, a1, a2, a3,
        h_pr_new(IPVX, n), h_pr_new(IPVY, n), h_pr_new(IPVZ, n));
    h_pr_new(IPDFWT, n) = 0.0;
    h_pr_new(IPT_BIRTH, n) = pm->time;
  }

  Kokkos::resize(ppart->prtcl_idata, ppart->nidata, new_npart);
  Kokkos::resize(ppart->prtcl_rdata, ppart->nrdata, new_npart);
  CopyPackedParticleDataToDevice(ppart, h_pi_new, h_pr_new, ninject, old_npart);
  ppart->nprtcl_thispack = new_npart;
  pm->CountParticles();
}

void ParallelShockSource(Mesh *pm, const Real bdt) {
  if (!ps_enable_injection || bdt <= 0.0) return;
  if (pm->time < ps_inject_t_start || pm->time > ps_inject_t_stop) return;
  if (ps_injection_transaction_cycle != pm->ncycle) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injection transaction was not prepared before "
              << "the RK source stage." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ApplyParallelShockGasSubtraction(pm, bdt/pm->dt);
}

void ValidateParallelShockInjectionTransaction(Mesh *pm) {
  if (!ps_enable_injection || !ps_enable_subtraction) return;
#if MPI_PARALLEL_ENABLED
  std::array<Real, 5> applied_global = {};
  Real delta_terms_global = 0.0;
  const Real delta_terms_local =
      static_cast<Real>(ps_injection_transaction_gas_deltas.size());
  MPI_Allreduce(ps_injection_transaction_applied_local.data(), applied_global.data(),
                5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&delta_terms_local, &delta_terms_global, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#else
  const std::array<Real, 5> &applied_global = ps_injection_transaction_applied_local;
  const Real delta_terms_global =
      static_cast<Real>(ps_injection_transaction_gas_deltas.size());
#endif
  // Count particle aggregation, cell aggregation and worst-case qualified
  // three-stage SSPRK shadow arithmetic, plus a conservative rank reduction depth.
  Real terms = ps_injection_transaction_terms_global +
      8.0*delta_terms_global + 4.0*global_variable::nranks;
  if (ps_test_source_transaction_terms_override) {
    terms = ps_test_source_transaction_terms;
  }
  for (int n=0; n<5; ++n) {
    if (!ParallelShockSourceTransactionValuesAgree(
            applied_global[n], ps_injection_transaction_expected_global[n],
            ps_injection_transaction_abs_global[n], terms)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock applied gas-subtraction transaction does not "
                << "match the injected-particle ledger." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  if (global_variable::my_rank == 0 && ps_feedback_diag_dcycle > 0 &&
      (pm->ncycle % ps_feedback_diag_dcycle) == 0) {
    std::cout << "pic_parallel_shock source_transaction_diag: cycle=" << pm->ncycle
              << " applied=(" << applied_global[0] << "," << applied_global[1]
              << "," << applied_global[2] << "," << applied_global[3] << ","
              << applied_global[4] << ") expected=("
              << ps_injection_transaction_expected_global[0] << ","
              << ps_injection_transaction_expected_global[1] << ","
              << ps_injection_transaction_expected_global[2] << ","
              << ps_injection_transaction_expected_global[3] << ","
              << ps_injection_transaction_expected_global[4] << ")" << std::endl;
  }
}

void ParallelShockRefinement(MeshBlockPack *pmbp) {
  if (!ps_enable_curvature_amr) return;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;
  if (!pmbp->pmesh->multilevel) return;
  if (!pmbp->pmesh->multi_d) return;

  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &w0 = pmbp->pmhd->w0;
  auto &refine_flag = pmbp->pmesh->pmr->refine_flag;
  const int mbs = pmbp->pmesh->gids_eachrank[global_variable::my_rank];
  const Real refine_thresh = ps_refine_curv;
  const Real deref_thresh = ps_derefine_curv;
  const Real tiny = 1.0e-30;

  par_for("ps_curvature_amr", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
  KOKKOS_LAMBDA(const int m) {
    Real max_grho = 0.0;
    Real max_gprs = 0.0;
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real rho = fmax(w0(m, IDN, k, j, i), tiny);
          Real prs = fmax(w0(m, IEN, k, j, i), tiny);

          Real gx_rho = fabs((w0(m, IDN, k, j, i + 1) +
                              w0(m, IDN, k, j, i - 1))/rho - 2.0);
          Real gy_rho = fabs((w0(m, IDN, k, j + 1, i) +
                              w0(m, IDN, k, j - 1, i))/rho - 2.0);
          Real gx_prs = fabs((w0(m, IEN, k, j, i + 1) +
                              w0(m, IEN, k, j, i - 1))/prs - 2.0);
          Real gy_prs = fabs((w0(m, IEN, k, j + 1, i) +
                              w0(m, IEN, k, j - 1, i))/prs - 2.0);

          Real grho = gx_rho + gy_rho;
          Real gprs = gx_prs + gy_prs;
          max_grho = fmax(max_grho, grho);
          max_gprs = fmax(max_gprs, gprs);
        }
      }
    }

    if (fmax(max_grho, max_gprs) > refine_thresh) {
      refine_flag.d_view(m + mbs) = 1;
    } else if (max_grho < deref_thresh && max_gprs < deref_thresh) {
      refine_flag.d_view(m + mbs) = -1;
    }
  });
}

void MaybePrintFeedbackDiagnostics(Mesh *pm) {
  if (pm == nullptr) return;
  if (ps_feedback_diag_dcycle < 1) return;
  if ((pm->ncycle % ps_feedback_diag_dcycle) != 0) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;
  auto *ppart = pmbp->ppart;
  if (!ppart->deposit_moments) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  constexpr int nfield = 8;
  const int idx[nfield] = {
      particles::Particles::IMOM_JX,
      particles::Particles::IMOM_JY,
      particles::Particles::IMOM_JZ,
      particles::Particles::IMOM_DPXDT,
      particles::Particles::IMOM_DPYDT,
      particles::Particles::IMOM_DPZDT,
      particles::Particles::IMOM_DEDT,
      particles::Particles::IMOM_EBDOT,
  };

  Real sum_abs_local[nfield] = {0.0};
  Real sum_sq_local[nfield] = {0.0};
  Real max_abs_local[nfield] = {0.0};
  std::int64_t ncell_local = 0;
  constexpr int nprt_field = 8;
  const int pr_idx[nprt_field] = {
      IPDPX, IPDPY, IPDPZ, IPDE, IPEBDOT, IPBX, IPBY, IPBZ};
  Real pr_sum_abs_local[nprt_field] = {0.0};
  Real pr_sum_sq_local[nprt_field] = {0.0};
  Real pr_max_abs_local[nprt_field] = {0.0};

  auto h_mom = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->moments);
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          for (int n = 0; n < nfield; ++n) {
            const Real v = h_mom(m, idx[n], k, j, i);
            const Real av = std::abs(v);
            sum_abs_local[n] += av;
            sum_sq_local[n] += v*v;
            if (av > max_abs_local[n]) max_abs_local[n] = av;
          }
          ncell_local += 1;
        }
      }
    }
  }
  auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
  for (int p = 0; p < ppart->nprtcl_thispack; ++p) {
    for (int n = 0; n < nprt_field; ++n) {
      const Real v = h_pr(pr_idx[n], p);
      const Real av = std::abs(v);
      pr_sum_abs_local[n] += av;
      pr_sum_sq_local[n] += v*v;
      if (av > pr_max_abs_local[n]) pr_max_abs_local[n] = av;
    }
  }

  std::int64_t npart_local = static_cast<std::int64_t>(ppart->nprtcl_thispack);
  std::int64_t npart_global = npart_local;
  std::int64_t ncell_global = ncell_local;
  Real sum_abs_global[nfield] = {0.0};
  Real sum_sq_global[nfield] = {0.0};
  Real max_abs_global[nfield] = {0.0};
  Real pr_sum_abs_global[nprt_field] = {0.0};
  Real pr_sum_sq_global[nprt_field] = {0.0};
  Real pr_max_abs_global[nprt_field] = {0.0};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(sum_abs_local, sum_abs_global, nfield, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(sum_sq_local, sum_sq_global, nfield, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(max_abs_local, max_abs_global, nfield, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(pr_sum_abs_local, pr_sum_abs_global, nprt_field, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pr_sum_sq_local, pr_sum_sq_global, nprt_field, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pr_max_abs_local, pr_max_abs_global, nprt_field, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&ncell_local, &ncell_global, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&npart_local, &npart_global, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
#else
  for (int n = 0; n < nfield; ++n) {
    sum_abs_global[n] = sum_abs_local[n];
    sum_sq_global[n] = sum_sq_local[n];
    max_abs_global[n] = max_abs_local[n];
  }
  for (int n = 0; n < nprt_field; ++n) {
    pr_sum_abs_global[n] = pr_sum_abs_local[n];
    pr_sum_sq_global[n] = pr_sum_sq_local[n];
    pr_max_abs_global[n] = pr_max_abs_local[n];
  }
#endif

  if (global_variable::my_rank != 0) return;
  const Real inv_ncell = (ncell_global > 0) ?
      1.0/static_cast<Real>(ncell_global) : 0.0;
  std::array<Real, nfield> rms{};
  for (int n = 0; n < nfield; ++n) {
    rms[n] = std::sqrt(std::max(sum_sq_global[n]*inv_ncell, 0.0));
  }
  const Real inv_npart = (npart_global > 0) ?
      1.0/static_cast<Real>(npart_global) : 0.0;
  std::array<Real, nprt_field> pr_rms{};
  for (int n = 0; n < nprt_field; ++n) {
    pr_rms[n] = std::sqrt(std::max(pr_sum_sq_global[n]*inv_npart, 0.0));
  }

  static bool printed_header = false;
  if (!printed_header) {
    const bool is_boris =
        (ppart->pusher == ParticlesPusher::boris_lin) ||
        (ppart->pusher == ParticlesPusher::boris_tsc);
    const bool use_delta_feedback =
        (is_boris &&
         (ppart->pic_feedback_mode == PICFeedbackMode::coupled) &&
         ppart->couple_moments_to_mhd &&
         (ppart->couple_moments_momentum_to_mhd || ppart->couple_moments_energy_to_mhd));
    std::cout << "pic_parallel_shock feedback_diag_cfg: "
              << "pusher=" << static_cast<int>(ppart->pusher)
              << " pic_feedback_mode=" << static_cast<int>(ppart->pic_feedback_mode)
              << " couple_moments_to_mhd=" << (ppart->couple_moments_to_mhd ? 1 : 0)
              << " mom_feedback=" << (ppart->couple_moments_momentum_to_mhd ? 1 : 0)
              << " eng_feedback=" << (ppart->couple_moments_energy_to_mhd ? 1 : 0)
              << " pic_enable_2d3v=" << (ppart->pic_enable_2d3v ? 1 : 0)
              << " use_delta_feedback=" << (use_delta_feedback ? 1 : 0)
              << std::endl;
    printed_header = true;
  }

  const Real step_mom_l1 = pm->dt*ppart->couple_moments_momentum_coeff*
      (sum_abs_global[3] + sum_abs_global[4] + sum_abs_global[5]);
  const Real step_eng_l1 = pm->dt*ppart->couple_moments_energy_coeff*
      sum_abs_global[6];
  auto old_flags = std::cout.flags();
  auto old_prec = std::cout.precision();
  std::cout << std::scientific << std::setprecision(12);
  std::cout << "pic_parallel_shock feedback_diag: cycle=" << pm->ncycle
            << " time=" << pm->time
            << " npart=" << npart_global
            << " ncell=" << ncell_global
            << " j_rms=(" << rms[0] << "," << rms[1] << "," << rms[2] << ")"
            << " dpdt_rms=(" << rms[3] << "," << rms[4] << "," << rms[5] << ")"
            << " dedt_rms=" << rms[6]
            << " ebdot_rms=" << rms[7]
            << " dpdt_l1=(" << sum_abs_global[3] << ","
            << sum_abs_global[4] << "," << sum_abs_global[5] << ")"
            << " dedt_l1=" << sum_abs_global[6]
            << " step_mom_l1=" << step_mom_l1
            << " step_eng_l1=" << step_eng_l1
            << " max_abs_dpdt=(" << max_abs_global[3] << ","
            << max_abs_global[4] << "," << max_abs_global[5] << ")"
            << " pr_dpdt_rms=(" << pr_rms[0] << "," << pr_rms[1] << ","
            << pr_rms[2] << ")"
            << " pr_dedt_rms=" << pr_rms[3]
            << " pr_ebdot_rms=" << pr_rms[4]
            << " pr_b_rms=(" << pr_rms[5] << "," << pr_rms[6] << ","
            << pr_rms[7] << ")"
            << " pr_dpdt_l1=(" << pr_sum_abs_global[0] << ","
            << pr_sum_abs_global[1] << "," << pr_sum_abs_global[2] << ")"
            << " pr_max_abs_dpdt=(" << pr_max_abs_global[0] << ","
            << pr_max_abs_global[1] << "," << pr_max_abs_global[2] << ")"
            << std::endl;
  std::cout.flags(old_flags);
  std::cout.precision(old_prec);
}

void CompleteParallelShockCycle(Mesh *pm) {
  auto *ppart = (pm != nullptr && pm->pmb_pack != nullptr) ?
      pm->pmb_pack->ppart : nullptr;
  ValidatePaperVL2CommittedEscapeChronology(ppart, pm->ncycle + 1,
                                            pm->time + pm->dt,
                                            "cycle completion", false);
  StoreRuntimeStateForRestart(pm->time + pm->dt);
  MaybePrintFeedbackDiagnostics(pm);
}

void ParallelShockCheckpoint(ParameterInput *pin, Mesh *pm) {
  (void)pin;
  if (pm == nullptr || pm->pmb_pack == nullptr) return;
  ValidatePaperVL2CommittedEscapeChronology(pm->pmb_pack->ppart, pm->ncycle,
                                            pm->time, "checkpoint");
  ValidateParallelShockParticlePopulation(pm);
  StoreRuntimeStateForRestart(pm->time);
}

void ParallelShockFinalize(ParameterInput *pin, Mesh *pm) {
  (void)pin;
  if (pm == nullptr || pm->pmb_pack == nullptr) return;
  ValidatePaperVL2CommittedEscapeChronology(pm->pmb_pack->ppart, pm->ncycle,
                                            pm->time, "run end");
  ValidateParallelShockParticlePopulation(pm);
  StoreRuntimeStateForRestart(pm->time);
  if (global_variable::my_rank == 0) {
    std::cout << "pic_parallel_shock escape_accounting_telemetry:"
              << " population_audit_calls=" << ps_particle_population_audit_calls
              << " destruction_audit_calls=" << ps_escape_audit_calls
              << " population_audit_policy=checkpoint_restart_run_end"
              << " production_pilot_required=1" << std::endl;
  }
}

void ParallelShockWorkInLoop(Mesh *pm) {
  if (pm == nullptr || pm->dt <= 0.0) return;
  if (ps_frame_diag_dcycle < 1) ps_frame_diag_dcycle = 1;
  ValidateParallelShockInjectionTransaction(pm);

  if (!ps_enable_frame_tracking) {
    CompleteParallelShockCycle(pm);
    return;
  }

  if (FrameModeVelocity()) {
    const Real v_now = FrameVelocityOffset(pm->time);
    const Real v_next_target = FrameVelocityOffset(pm->time + pm->dt);
    Real dv = v_next_target - v_now;
    const bool clipped = (ps_frame_dv_max > 0.0 && std::abs(dv) > ps_frame_dv_max);
    if (ps_frame_dv_max > 0.0 && std::abs(dv) > ps_frame_dv_max) {
      dv = std::copysign(ps_frame_dv_max, dv);
    }
    if (std::abs(dv) <= 1.0e-15) {
      CompleteParallelShockCycle(pm);
      return;
    }

    ApplyFrameShiftToFluid(pm, dv);
    if (ps_frame_apply_to_particles) {
      ApplyFrameShiftToParticles(pm, dv);
    }
    if (ps_frame_apply_to_inflow) {
      UpdateOuterInflowState(pm, v_next_target);
    }

    if (global_variable::my_rank == 0 &&
        ((pm->ncycle % ps_frame_diag_dcycle) == 0)) {
      const Real x_model = ShockSurfaceModelX1(pm->time + pm->dt);
      std::cout << "pic_parallel_shock frame(velocity): cycle=" << pm->ncycle
                << " time=" << pm->time
                << " v_now=" << v_now
                << " v_target=" << v_next_target
                << " dv=" << dv
                << " clipped=" << (clipped ? 1 : 0)
                << " xshock_model=" << x_model << std::endl;
    }
    CompleteParallelShockCycle(pm);
    return;
  }

  if (FrameModeRecenter()) {
    const int ev_now = RecenterEventCount(pm->time);
    const int ev_next = RecenterEventCount(pm->time + pm->dt);
    const int dev = ev_next - ev_now;
    int nshift = 0;
    if (dev > 0) {
      for (int n = 0; n < dev; ++n) {
        ApplyRecenteringShift(pm, ps_recenter_shift_cells);
      }
      nshift = dev*ps_recenter_shift_cells;
    }

    if (global_variable::my_rank == 0 &&
        ((pm->ncycle % ps_frame_diag_dcycle) == 0 || dev > 0)) {
      const Real x_model = ShockSurfaceModelX1(pm->time + pm->dt);
      const Real xshift = static_cast<Real>(nshift)*ps_recenter_dx1;
      std::cout << "pic_parallel_shock frame(recenter): cycle=" << pm->ncycle
                << " time=" << pm->time
                << " ev_now=" << ev_now
                << " ev_next=" << ev_next
                << " dev=" << dev
                << " nshift=" << nshift
                << " xshift=" << xshift
                << " xshock_model=" << x_model << std::endl;
    }
  }
  CompleteParallelShockCycle(pm);
}

void ParallelShockWorkBeforeLoop(Mesh *pm) {
  if (pm == nullptr || pm->dt <= 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp != nullptr && pmbp->ppart != nullptr) {
    SeedNextTag(pmbp->ppart, pm->time);
  }
  RemoveExcludedEarlyInjectedParticles(pm);
  PrepareParallelShockInjectionTransaction(pm);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::PICParallelShock()
//! \brief Parallel-shock MHD-PIC benchmark setup.

void ProblemGenerator::PICParallelShock(ParameterInput *pin, const bool restart) {
  ps_pin = pin;
  ps_particle_population_audit_calls = 0;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires an active <mhd> block." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (pmbp->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires an active <particles> block."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (pmbp->phydro != nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock is MHD-only; do not enable <hydro>."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ideal MHD EOS." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::reflect) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects <mesh>/ix1_bc=reflect." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const BoundaryFlag ox1_bc = pmy_mesh_->mesh_bcs[BoundaryFace::outer_x1];
  if (!(ox1_bc == BoundaryFlag::inflow || ox1_bc == BoundaryFlag::outflow)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects <mesh>/ox1_bc=inflow or outflow."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
      pmy_mesh_->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects periodic y boundaries." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  // Parse runtime controls.
  ps_rho0 = pin->GetOrAddReal("problem", "ps_rho0", 1.0);
  ps_p0 = pin->GetOrAddReal("problem", "ps_p0", 1.0);
  ps_u0 = pin->GetOrAddReal("problem", "ps_u0", 30.0);
  ps_b0 = pin->GetOrAddReal("problem", "ps_b0", 1.0);
  ps_eta = pin->GetOrAddReal("problem", "ps_eta", 1.0e-3);
  ps_vinj_over_u0 = pin->GetOrAddReal("problem", "ps_vinj_over_u0", std::sqrt(10.0));
  ps_inject_half_width_cells = pin->GetOrAddReal(
      "problem", "ps_inject_half_width_cells", 0.5);
  ps_inject_t_start = pin->GetOrAddReal("problem", "ps_inject_t_start", 0.0);
  ps_inject_t_stop = pin->GetOrAddReal("problem", "ps_inject_t_stop", 1.0e99);
  ps_remove_birth_time_before = pin->GetOrAddReal(
      "problem", "ps_remove_birth_time_before", -1.0);
  std::string shock_speed_model = pin->GetOrAddString(
      "problem", "ps_shock_speed_model", "finite_mach");
  ps_seed_noise_amp = pin->GetOrAddReal("problem", "ps_seed_noise_amp", 0.0);
  ps_seed_noise_seed = pin->GetOrAddInteger("problem", "ps_seed_noise_seed", 1234);
  ps_refine_curv = pin->GetOrAddReal("problem", "ps_refine_curv", 1.0);
  ps_derefine_curv = pin->GetOrAddReal("problem", "ps_derefine_curv", 0.1);
  ps_rho_floor_frac = pin->GetOrAddReal("problem", "ps_rho_floor_frac", 1.0e-6);
  ps_p_floor_frac = pin->GetOrAddReal("problem", "ps_p_floor_frac", 1.0e-8);
  ps_enable_injection = pin->GetOrAddBoolean("problem", "ps_enable_injection", true);
  ps_enable_subtraction = pin->GetOrAddBoolean(
      "problem", "ps_enable_gas_subtraction", true);
  ps_enable_curvature_amr = pin->GetOrAddBoolean(
      "problem", "ps_enable_curvature_amr", true);
  ps_test_source_transaction_terms_override = pin->GetOrAddBoolean(
      "problem", "ps_test_source_transaction_terms_override", false);
  ps_test_source_transaction_terms = pin->GetOrAddReal(
      "problem", "ps_test_source_transaction_terms", 1.0);
  ps_enable_frame_tracking = pin->GetOrAddBoolean(
      "problem", "ps_enable_frame_tracking", false);
  std::string frame_mode = pin->GetOrAddString("problem", "ps_frame_mode",
                                                "velocity");
  ps_frame_t_start = pin->GetOrAddReal("problem", "ps_frame_t_start", 0.0);
  ps_frame_t_ramp = pin->GetOrAddReal("problem", "ps_frame_t_ramp", 0.0);
  ps_frame_vfrac = pin->GetOrAddReal("problem", "ps_frame_vfrac", 1.0);
  ps_frame_dv_max = pin->GetOrAddReal("problem", "ps_frame_dv_max", 1.0e99);
  ps_frame_apply_to_particles = pin->GetOrAddBoolean(
      "problem", "ps_frame_apply_to_particles", true);
  ps_frame_apply_to_inflow = pin->GetOrAddBoolean(
      "problem", "ps_frame_apply_to_inflow", true);
  ps_frame_require_uniform = pin->GetOrAddBoolean(
      "problem", "ps_frame_require_uniform", true);
  ps_frame_diag_dcycle = pin->GetOrAddInteger(
      "problem", "ps_frame_diag_dcycle", 200);
  ps_feedback_diag_dcycle = pin->GetOrAddInteger(
      "problem", "ps_feedback_diag_dcycle", 0);
  ps_recenter_x_target = pin->GetOrAddReal("problem", "ps_recenter_x_target", 2.0);
  ps_recenter_x_trigger = pin->GetOrAddReal("problem", "ps_recenter_x_trigger", 3.0);
  ps_recenter_shift_cells = pin->GetOrAddInteger("problem",
                                                  "ps_recenter_shift_cells", 2);
  ps_recenter_vshock_model = pin->GetOrAddReal("problem",
                                                "ps_recenter_vshock_model", -1.0);
  ps_inject_species = pin->GetOrAddInteger("problem", "ps_inject_species", 0);
  ps_inject_seed = pin->GetOrAddInteger("problem", "ps_inject_seed", 1234);

  const std::string integrator = pin->GetString("time", "integrator");
  if (ps_enable_injection && integrator != "rk1" && integrator != "rk2"
      && integrator != "rk3") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injection is qualified only with time/integrator="
              << "rk1, rk2, or rk3; received '" << integrator << "'." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  if (!std::isfinite(ps_rho0) || !std::isfinite(ps_p0) || !std::isfinite(ps_u0) ||
      !std::isfinite(ps_b0) || !std::isfinite(ps_eta) ||
      ps_rho0 <= 0.0 || ps_p0 <= 0.0 || ps_u0 <= 0.0 || ps_b0 <= 0.0 || ps_eta < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock parameters ps_rho0/ps_p0/ps_u0/ps_b0 must be > 0 "
              << "and ps_eta must be >= 0." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!std::isfinite(ps_rho_floor_frac) || !std::isfinite(ps_p_floor_frac) ||
      ps_rho_floor_frac < 0.0 || ps_p_floor_frac < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock floor fractions must be finite and non-negative."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_enable_subtraction &&
      pmbp->ppart->pic_background_mode != PICBackgroundMode::coupled) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock gas subtraction requires "
              << "<particles>/pic_background_mode=coupled." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!std::isfinite(ps_vinj_over_u0) || ps_vinj_over_u0 <= 0.0 ||
      !std::isfinite(ps_inject_half_width_cells) ||
      !std::isfinite(ps_inject_t_start) || !std::isfinite(ps_inject_t_stop) ||
      ps_inject_t_stop < ps_inject_t_start ||
      !std::isfinite(ps_remove_birth_time_before)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock injection controls must be finite, "
              << "ps_vinj_over_u0 must be positive, and the injection interval "
              << "must be ordered." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!std::isfinite(ps_seed_noise_amp) || ps_seed_noise_amp < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_seed_noise_amp >= 0."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (std::abs(ps_inject_half_width_cells - 0.5) > 1.0e-15) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_inject_half_width_cells = 0.5 "
              << "for a unique shock-surface carrier cell." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (shock_speed_model == "finite_mach") {
    ps_shock_speed_model = PSShockSpeedModel::finite_mach;
  } else if (shock_speed_model == "ideal_surface") {
    ps_shock_speed_model = PSShockSpeedModel::ideal_surface;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "ps_shock_speed_model must be 'finite_mach' or "
              << "'ideal_surface'." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!std::isfinite(ps_frame_t_start) || !std::isfinite(ps_frame_t_ramp) ||
      !std::isfinite(ps_frame_vfrac) || !std::isfinite(ps_frame_dv_max) ||
      ps_frame_t_ramp < 0.0 || ps_frame_vfrac < 0.0 || ps_frame_dv_max < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock frame controls require ps_frame_t_ramp >= 0, "
              << "ps_frame_vfrac >= 0, and ps_frame_dv_max >= 0." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_frame_diag_dcycle < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_frame_diag_dcycle >= 1."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_feedback_diag_dcycle < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_feedback_diag_dcycle >= 0."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (frame_mode == "velocity") {
    ps_frame_mode = PSFrameMode::velocity;
  } else if (frame_mode == "recenter") {
    ps_frame_mode = PSFrameMode::recenter;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "ps_frame_mode must be 'velocity' or 'recenter'." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  const int nspecies = pmbp->ppart->nspecies;
  if (ps_inject_species < 0 || ps_inject_species >= nspecies) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "ps_inject_species is out of range for configured species."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  std::string species_block = "species" + std::to_string(ps_inject_species);
  ps_particle_mass = pin->GetOrAddReal(species_block, "mass", 1.0);
  ps_particle_charge = pin->GetOrAddReal(species_block, "charge", 1.0);
  if (!std::isfinite(ps_particle_mass) || !std::isfinite(ps_particle_charge) ||
      ps_particle_mass <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Injected species mass must be > 0." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ps_particle_q_over_m = ps_particle_charge/ps_particle_mass;
  Real qscale = pin->GetOrAddReal("particles", "deposit_qscale", 1.0);
  ps_particle_macro_mass = qscale*ps_particle_mass;
  if (!std::isfinite(ps_particle_macro_mass) || ps_particle_macro_mass <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Computed injected macro-mass must be > 0." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  // Select the source-local surface model in the reflecting-wall frame.
  Real gamma = pmbp->pmhd->peos->eos_data.gamma;
  ps_shock_speed = (ps_shock_speed_model == PSShockSpeedModel::ideal_surface) ?
      IdealSurfaceShockSpeed(gamma, ps_u0) :
      EstimateShockSpeed(gamma, ps_rho0, ps_p0, ps_u0);
  if (!std::isfinite(ps_shock_speed) || ps_shock_speed <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires supersonic inflow with positive "
              << "finite-Mach shock-speed estimate." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ps_xshock0 = pmy_mesh_->mesh_size.x1min;
  ps_recenter_dx1 = pmy_mesh_->mesh_size.dx1;
  ps_use_2d3v = (pmy_mesh_->two_d && pmbp->ppart->pic_enable_2d3v);
  ValidateAndStoreParallelShockRestartControls(pin, restart);
  RejectDuplicateParallelShockRestartLedgers(pin, restart);
  const bool has_ledger_schema =
      restart && pin->DoesParameterExist("problem", "ps_cr_ledger_schema");
  if (has_ledger_schema) {
    const int ledger_schema = pin->GetInteger("problem", "ps_cr_ledger_schema");
    const std::array<const char *, 18> ledger_fields = {
      "ps_cr_ledger_complete", "ps_mass_reservoir_global",
      "ps_injected_cr_count_global", "ps_injected_cr_mass_global",
      "ps_injected_cr_momentum_x1_global", "ps_injected_cr_momentum_x2_global",
      "ps_injected_cr_momentum_x3_global", "ps_injected_cr_energy_global",
      "ps_removed_excluded_early_cohort",
      "ps_removed_cr_count_global", "ps_removed_cr_mass_global",
      "ps_removed_cr_momentum_x1_global", "ps_removed_cr_momentum_x2_global",
      "ps_removed_cr_momentum_x3_global", "ps_removed_cr_energy_global",
      "ps_tag_seeded", "ps_injection_tag_floor", "ps_next_tag"
    };
    bool ledger_fields_complete = true;
    for (const char *field : ledger_fields) {
      ledger_fields_complete =
          ledger_fields_complete && pin->DoesParameterExist("problem", field);
    }
    if (ledger_schema != 3 || !ledger_fields_complete) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart has unsupported or incomplete CR "
                << "ledger metadata." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ps_cr_ledger_complete = pin->GetBoolean("problem", "ps_cr_ledger_complete");
  } else {
    if (restart) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart requires complete schema-3 CR ledger "
                << "metadata." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ps_cr_ledger_complete = true;
  }
  if (!ps_cr_ledger_complete) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock restart CR ledger is explicitly incomplete."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ps_mass_reservoir_global = restart ?
      pin->GetOrAddReal("problem", "ps_mass_reservoir_global", 0.0) : 0.0;
  ps_injected_cr_count_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_count_global", 0.0) : 0.0;
  ps_injected_cr_mass_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_mass_global", 0.0) : 0.0;
  ps_injected_cr_momentum_x1_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_momentum_x1_global", 0.0) :
      0.0;
  ps_injected_cr_momentum_x2_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_momentum_x2_global", 0.0) :
      0.0;
  ps_injected_cr_momentum_x3_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_momentum_x3_global", 0.0) :
      0.0;
  ps_injected_cr_energy_global = restart ?
      pin->GetOrAddReal("problem", "ps_injected_cr_energy_global", 0.0) : 0.0;
  ps_injection_transaction_cycle = std::numeric_limits<int>::min();
  ps_injection_transaction_gas_deltas.clear();
  ps_removed_excluded_early_cohort = restart ?
      pin->GetOrAddBoolean("problem", "ps_removed_excluded_early_cohort", false) :
      false;
  ps_removed_cr_count_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_count_global", 0.0) : 0.0;
  ps_removed_cr_mass_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_mass_global", 0.0) : 0.0;
  ps_removed_cr_momentum_x1_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_momentum_x1_global", 0.0) : 0.0;
  ps_removed_cr_momentum_x2_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_momentum_x2_global", 0.0) : 0.0;
  ps_removed_cr_momentum_x3_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_momentum_x3_global", 0.0) : 0.0;
  ps_removed_cr_energy_global = restart ?
      pin->GetOrAddReal("problem", "ps_removed_cr_energy_global", 0.0) : 0.0;
  const bool has_escape_ledger_schema =
      restart && pin->DoesParameterExist("problem", "ps_escape_ledger_schema");
  if (restart) {
    const std::array<const char *, 10> escape_ledger_fields = {
      "ps_escape_ledger_complete", "ps_escape_audit_calls",
      "ps_escape_last_audit_time", "ps_escaped_injected_cr_count_global",
      "ps_escaped_injected_cr_mass_global",
      "ps_escaped_injected_cr_momentum_x1_global",
      "ps_escaped_injected_cr_momentum_x2_global",
      "ps_escaped_injected_cr_momentum_x3_global",
      "ps_escaped_injected_cr_energy_global",
      "ps_escaped_initial_cr_count_global"
    };
    bool escape_ledger_fields_complete = has_escape_ledger_schema;
    for (const char *field : escape_ledger_fields) {
      escape_ledger_fields_complete =
          escape_ledger_fields_complete && pin->DoesParameterExist("problem", field);
    }
    if (!escape_ledger_fields_complete ||
        pin->GetInteger("problem", "ps_escape_ledger_schema") != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart requires complete schema-1 particle "
                << "escape ledger metadata." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  ps_escape_ledger_complete = restart ?
      pin->GetBoolean("problem", "ps_escape_ledger_complete") : true;
  ps_escape_audit_calls = restart ?
      pin->GetInteger("problem", "ps_escape_audit_calls") : 0;
  ps_escape_last_audit_time = restart ?
      pin->GetReal("problem", "ps_escape_last_audit_time") : 0.0;
  ps_escaped_injected_cr_count_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_count_global") : 0.0;
  ps_escaped_injected_cr_mass_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_mass_global") : 0.0;
  ps_escaped_injected_cr_momentum_x1_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_momentum_x1_global") : 0.0;
  ps_escaped_injected_cr_momentum_x2_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_momentum_x2_global") : 0.0;
  ps_escaped_injected_cr_momentum_x3_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_momentum_x3_global") : 0.0;
  ps_escaped_injected_cr_energy_global = restart ?
      pin->GetReal("problem", "ps_escaped_injected_cr_energy_global") : 0.0;
  ps_escaped_initial_cr_count_global = restart ?
      pin->GetReal("problem", "ps_escaped_initial_cr_count_global") : 0.0;
  ps_tag_seeded = has_ledger_schema ?
      pin->GetBoolean("problem", "ps_tag_seeded") : false;
  ps_tag_progression_validated = false;
  ps_injection_tag_floor = 0;
  if (has_ledger_schema ||
      (restart && pin->DoesParameterExist("problem", "ps_injection_tag_floor"))) {
    const int stored_tag_floor =
        pin->GetInteger("problem", "ps_injection_tag_floor");
    if (stored_tag_floor < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart metadata has negative "
                << "ps_injection_tag_floor." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ps_injection_tag_floor = static_cast<std::int64_t>(stored_tag_floor);
  }
  ps_next_tag = 0;
  if (has_ledger_schema ||
      (restart && pin->DoesParameterExist("problem", "ps_next_tag"))) {
    const int stored_next_tag = pin->GetInteger("problem", "ps_next_tag");
    if (stored_next_tag < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock restart metadata has negative ps_next_tag."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ps_next_tag = static_cast<std::int64_t>(stored_next_tag);
    if (!has_ledger_schema) {
      ps_tag_seeded = true;
    }
  }
  ValidateParallelShockRuntimeLedger(restart ? "restart" : "runtime",
                                     pmy_mesh_->time);
  ValidatePaperVL2CommittedEscapeChronology(pmbp->ppart, pmy_mesh_->ncycle,
                                            pmy_mesh_->time,
                                            restart ? "restart" : "initial state");
  if (restart && pmbp->ppart != nullptr) {
    SeedNextTag(pmbp->ppart, pmy_mesh_->time);
    ValidateParallelShockParticlePopulation(pmy_mesh_);
  }
  StoreRuntimeStateForRestart(pmy_mesh_->time);
  ConfigureSeedNoisePhases();

  if (ps_enable_frame_tracking && ps_frame_require_uniform &&
      (pmy_mesh_->adaptive || pmy_mesh_->multilevel)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock frame tracking currently supports uniform grids "
              << "only (no AMR/SMR)." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ps_enable_frame_tracking && FrameModeRecenter()) {
    if (ps_frame_apply_to_particles && pmbp->ppart->UsesPaperVL2Coupling()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock paper-mode recentering with particle shifts "
                << "is not qualified because its coordinate-removal sink is not "
                << "included in the physical-boundary escape ledger." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (!(ps_recenter_x_target > pmy_mesh_->mesh_size.x1min &&
          ps_recenter_x_target < pmy_mesh_->mesh_size.x1max)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_target must lie strictly inside the x1 domain."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (!(ps_recenter_x_trigger > pmy_mesh_->mesh_size.x1min)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_trigger must be greater than x1min."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (!(ps_recenter_x_trigger < pmy_mesh_->mesh_size.x1max)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_trigger must be less than x1max."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (ps_recenter_dx1 <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock recentering requires positive mesh dx1."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (ps_recenter_shift_cells < 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_shift_cells must be >= 1 in recenter mode."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (ps_recenter_shift_cells > pmy_mesh_->mb_indcs.ng) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_shift_cells must be <= nghost in recenter mode."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    const Real vmodel = ps_recenter_vshock_model;
    const bool use_default = std::abs(vmodel + 1.0) < 1.0e-12;
    const bool use_override = vmodel > 0.0;
    if (!(use_default || use_override)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_vshock_model must be > 0 or left at default -1."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }

  // Enroll benchmark callbacks.
  user_srcs = true;
  user_srcs_func = ParallelShockSource;
  user_ref_func = ParallelShockRefinement;
  user_work_before_loop_func = ParallelShockWorkBeforeLoop;
  user_work_in_loop = true;
  user_work_in_loop_func = ParallelShockWorkInLoop;
  pgen_checkpoint_func = ParallelShockCheckpoint;
  pgen_final_func = ParallelShockFinalize;
  pmbp->ppart->particle_destruction_observer =
      ObserveParallelShockParticleDestruction;

  // The inflow reservoir is not stored in restart files, so rebuild it before
  // Driver::Initialize() fills ghost zones on both new and restarted runs.
  UpdateOuterInflowState(pmy_mesh_, FrameVelocityOffset(pmy_mesh_->time));

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mb_size;
  Real gm1 = gamma - 1.0;
  const bool add_seed_noise = (ps_seed_noise_amp > 0.0);
  const Real seed_noise_amp = ps_seed_noise_amp*ps_b0;
  const Real x1min_global = pmy_mesh_->mesh_size.x1min;
  const Real lx_global = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const Real seed_noise_norm = 1.0/1.875;
  const Real ph_by_1 = ps_seed_noise_phase_by[0];
  const Real ph_by_2 = ps_seed_noise_phase_by[1];
  const Real ph_by_3 = ps_seed_noise_phase_by[2];
  const Real ph_by_4 = ps_seed_noise_phase_by[3];
  const Real ph_bz_1 = ps_seed_noise_phase_bz[0];
  const Real ph_bz_2 = ps_seed_noise_phase_bz[1];
  const Real ph_bz_3 = ps_seed_noise_phase_bz[2];
  const Real ph_bz_4 = ps_seed_noise_phase_bz[3];
  const Real upstream_rho0 = ps_rho0;
  const Real upstream_p0 = ps_p0;
  const Real upstream_u0 = ps_u0;
  const Real upstream_b0 = ps_b0;

  // Set uniform upstream state (flow toward the reflecting wall).
  par_for("pgen_pic_parallel_shock", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    u0(m, IDN, k, j, i) = upstream_rho0;
    u0(m, IM1, k, j, i) = -upstream_rho0*upstream_u0;
    u0(m, IM2, k, j, i) = 0.0;
    u0(m, IM3, k, j, i) = 0.0;

    Real dby = 0.0;
    Real dbz = 0.0;
    if (add_seed_noise) {
      const Real x1c = size.d_view(m).x1min +
          (static_cast<Real>(i - is) + 0.5)*size.d_view(m).dx1;
      const Real xfrac = (x1c - x1min_global)/lx_global;
      const Real xphase = 2.0*M_PI*xfrac;

      const Real nby = sin(1.0*xphase + ph_by_1) +
                       0.5*sin(2.0*xphase + ph_by_2) +
                       0.25*sin(4.0*xphase + ph_by_3) +
                       0.125*sin(8.0*xphase + ph_by_4);
      const Real nbz = sin(1.0*xphase + ph_bz_1) +
                       0.5*sin(2.0*xphase + ph_bz_2) +
                       0.25*sin(4.0*xphase + ph_bz_3) +
                       0.125*sin(8.0*xphase + ph_bz_4);
      dby = seed_noise_amp*seed_noise_norm*nby;
      dbz = seed_noise_amp*seed_noise_norm*nbz;
    }

    b0.x1f(m, k, j, i) = upstream_b0;
    b0.x2f(m, k, j, i) = dby;
    b0.x3f(m, k, j, i) = dbz;
    if (i == ie) b0.x1f(m, k, j, i + 1) = upstream_b0;
    if (j == je) b0.x2f(m, k, j + 1, i) = dby;
    if (k == ke) b0.x3f(m, k + 1, j, i) = dbz;
  });

  par_for("pgen_pic_parallel_shock_e", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real bx = 0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    Real by = 0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    Real bz = 0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
    Real ekin = 0.5*(SQR(u0(m, IM1, k, j, i)) +
                     SQR(u0(m, IM2, k, j, i)) +
                     SQR(u0(m, IM3, k, j, i)))/u0(m, IDN, k, j, i);
    Real emag = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    u0(m, IEN, k, j, i) = upstream_p0/gm1 + ekin + emag;
  });

  // Keep size data synced for host-side injection cell selection.
  size.template modify<HostMemSpace>();
  size.template sync<DevExeSpace>();
}
