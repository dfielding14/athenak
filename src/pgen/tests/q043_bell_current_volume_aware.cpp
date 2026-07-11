//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q043_bell_current_volume_aware.cpp
//! \brief Volume-aware Bell current oracle and corrected Section 5.2 eigenmode.

#include <algorithm>
#include <cmath>

#if !defined(Q043_BELL_CURRENT_VOLUME_AWARE_HOST_CONTRACT)
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#endif

#if defined(Q043_BELL_CURRENT_VOLUME_AWARE_HOST_CONTRACT) && \
    !defined(Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT)
#define Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT 1
#define Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT 1
#endif
#if defined(Q043_BELL_CURRENT_VOLUME_AWARE_HOST_CONTRACT) && \
    !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#define Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#define Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#endif
#include "q023_paper_bell_linear_joverc.hpp"
#if defined(Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#undef Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#undef Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#endif
#if defined(Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT)
#undef Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT
#undef Q043_VOLUME_AWARE_UNDEF_Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT
#endif

namespace q043_bell_current_volume_aware {

using q023_paper_bell_linear::Basis;
using q023_paper_bell_linear::ModeBasis;
using q023_paper_bell_linear::ModeParameters;
using q023_paper_bell_linear::Vector3;
using q023_paper_bell_linear_joverc::UnstableEigenmodeAtPhase;
using q023_paper_bell_linear_joverc::UnstableVectorPotentialAt;

inline double RootCellVolume(const double x1_extent, const int root_nx1,
                             const double x2_extent, const int root_nx2,
                             const double x3_extent, const int root_nx3) {
  return (x1_extent/static_cast<double>(root_nx1))*
      (x2_extent/static_cast<double>(root_nx2))*
      (x3_extent/static_cast<double>(root_nx3));
}

enum class SourceMode {
  uniform_current_oracle,
  corrected_linear_eigenmode
};

inline bool SourceModeAmplitudeIsValid(const SourceMode mode, const double amplitude) {
  if (!std::isfinite(amplitude)) return false;
  if (mode == SourceMode::uniform_current_oracle) return amplitude == 0.0;
  return amplitude > 0.0 && amplitude <= 1.0e-2;
}

inline bool SourceModeSpeciesMassIsValid(const SourceMode mode,
                                         const double species_mass) {
  if (!std::isfinite(species_mass) || species_mass <= 0.0) return false;
  if (mode == SourceMode::uniform_current_oracle) return true;
  return std::abs(species_mass - 1.0) <= 1.0e-12;
}

inline bool HasRequiredSpeciesChargeOverMass(const double species_mass,
                                             const double species_charge,
                                             const double expected_q_over_mc) {
  if (!std::isfinite(species_mass) || species_mass <= 0.0 ||
      !std::isfinite(species_charge) || !std::isfinite(expected_q_over_mc)) {
    return false;
  }
  const double measured_q_over_mc = species_charge/species_mass;
  const double scale = std::max(1.0, std::abs(expected_q_over_mc));
  return std::abs(measured_q_over_mc - expected_q_over_mc) <= 1.0e-12*scale;
}

inline bool PositiveIntegralPPCIsValid(const double ppc) {
  return std::isfinite(ppc) && ppc >= 1.0 && std::floor(ppc) == ppc;
}

inline double DepositedJOverC(const double ppc, const double qscale,
                              const double species_charge, const double stream_speed,
                              const double root_cell_volume) {
  return ppc*qscale*species_charge*stream_speed/root_cell_volume;
}

inline double RequiredDepositedJOverC(const double b_g, const double k0) {
  return 2.0*b_g*k0;
}

inline double RequiredDepositQScale(const double ppc, const double species_charge,
                                    const double stream_speed,
                                    const double root_cell_volume,
                                    const double b_g, const double k0) {
  return RequiredDepositedJOverC(b_g, k0)*root_cell_volume/
      (ppc*species_charge*stream_speed);
}

inline bool HasRequiredDepositedJOverC(const double ppc, const double qscale,
                                       const double species_charge,
                                       const double stream_speed,
                                       const double root_cell_volume,
                                       const double b_g, const double k0) {
  const double measured =
      DepositedJOverC(ppc, qscale, species_charge, stream_speed, root_cell_volume);
  const double expected = RequiredDepositedJOverC(b_g, k0);
  const double scale = std::max(1.0, std::abs(expected));
  return std::isfinite(measured) && std::isfinite(expected) &&
      std::abs(measured - expected) <= 1.0e-12*scale;
}

}  // namespace q043_bell_current_volume_aware

#if !defined(Q043_BELL_CURRENT_VOLUME_AWARE_HOST_CONTRACT)
namespace {

[[noreturn]] void Q043VolumeAwareFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q043VolumeAwareRequireClose(const std::string &label, const Real measured,
                            const Real expected) {
  if (!std::isfinite(measured) || !std::isfinite(expected)) {
    Q043VolumeAwareFatal(label + " must be finite for the corrected Q-043 "
                   "Section 5.2 contract");
  }
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q043VolumeAwareFatal(label + " does not match the corrected Q-043 Section 5.2 "
                   "contract");
  }
}

void Q043VolumeAwareRequireFinite(const std::string &label, const Real value) {
  if (!std::isfinite(value)) {
    Q043VolumeAwareFatal(label + " must be finite for the corrected Q-043 "
                   "Section 5.2 contract");
  }
}

void Q043VolumeAwareRequireString(ParameterInput *pin, const std::string &block,
                             const std::string &name,
                             const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q043VolumeAwareFatal("<" + block + ">/" + name +
                    " does not match the corrected Q-043 Section 5.2 contract");
  }
}

void Q043VolumeAwareRequireBoolean(ParameterInput *pin, const std::string &block,
                              const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q043VolumeAwareFatal("<" + block + ">/" + name +
                    " does not match the corrected Q-043 Section 5.2 contract");
  }
}

}  // namespace

void ProblemGenerator::Q043BellCurrentVolumeAware(ParameterInput *pin,
                                                  const bool restart) {
  using q043_bell_current_volume_aware::Basis;
  using q043_bell_current_volume_aware::HasRequiredSpeciesChargeOverMass;
  using q043_bell_current_volume_aware::HasRequiredDepositedJOverC;
  using q043_bell_current_volume_aware::ModeBasis;
  using q043_bell_current_volume_aware::ModeParameters;
  using q043_bell_current_volume_aware::PositiveIntegralPPCIsValid;
  using q043_bell_current_volume_aware::SourceMode;
  using q043_bell_current_volume_aware::SourceModeAmplitudeIsValid;
  using q043_bell_current_volume_aware::SourceModeSpeciesMassIsValid;
  using q043_bell_current_volume_aware::UnstableEigenmodeAtPhase;
  using q043_bell_current_volume_aware::UnstableVectorPotentialAt;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware requires ideal-MHD EOS");
  }
  if (pmy_mesh_->multilevel) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware rejects AMR/SMR");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware requires periodic boundaries");
  }

  const std::string block = "q043_bell_current_volume_aware";
  const int dimension = pin->GetInteger(block, "dimension");
  const int mesh_dimension = pmy_mesh_->one_d ? 1 : (pmy_mesh_->two_d ? 2 : 3);
  const bool thin_2d3v_1d_carrier =
      dimension == 1 && mesh_dimension == 2 &&
      pmy_mesh_->mb_indcs.nx2 == 4 && pmy_mesh_->mb_indcs.nx3 == 1;
  const bool native_dimension =
      dimension >= 2 && dimension <= 3 && dimension == mesh_dimension;
  if (!(thin_2d3v_1d_carrier || native_dimension)) {
    Q043VolumeAwareFatal("<q043_bell_current_volume_aware>/dimension does not match "
                    "the required native mesh or exact thin 2D3V carrier");
  }

  Q043VolumeAwareRequireString(pin, block, "campaign_id",
                          "Q043-BELL-CURRENT-VOLUME-AWARE");
  Q043VolumeAwareRequireString(pin, block, "supersedes_campaign_id",
                          "Q023-PAPER-BELL-LINEAR");
  Q043VolumeAwareRequireString(pin, block, "current_normalization",
                          "deposited_j_over_c_equals_2_b_g_k0");
  Q043VolumeAwareRequireString(
      pin, block, "deposition_measure",
      "ppc_times_deposit_qscale_times_species_charge_times_v_cr_over_root_cell_volume");
  Q043VolumeAwareRequireString(
      pin, block, "root_cell_volume",
      "global_root_mesh_extents_over_global_root_mesh_counts");
  Q043VolumeAwareRequireString(
      pin, block, "deposit_qscale_semantics",
      "root_cell_macro_charge_volume_aware");
  Q043VolumeAwareRequireString(
      pin, block, "fixed_qscale",
      "forbidden_across_dimensions_and_resolutions");
  Q043VolumeAwareRequireString(pin, block, "epsilon_grid", "0.1,0.2,0.4,0.6,0.8");
  Q043VolumeAwareRequireString(pin, block, "timestep",
                          "open_clean_candidate_timestep_freeze");
  const std::string source_mode_name = pin->GetString(block, "source_mode");
  const bool uniform_current_oracle =
      source_mode_name.compare("uniform_current_oracle") == 0;
  const bool corrected_linear_eigenmode =
      source_mode_name.compare("corrected_linear_eigenmode") == 0;
  if (!uniform_current_oracle && !corrected_linear_eigenmode) {
    Q043VolumeAwareFatal("<q043_bell_current_volume_aware>/source_mode must be "
                         "uniform_current_oracle or corrected_linear_eigenmode");
  }
  const SourceMode source_mode = uniform_current_oracle
      ? SourceMode::uniform_current_oracle
      : SourceMode::corrected_linear_eigenmode;
  Q043VolumeAwareRequireString(
      pin, block, "q_over_mc_representation",
      uniform_current_oracle
          ? "species_charge_over_species_mass_equals_omega_over_b_g"
          : "species_charge_equals_q_over_mc_only_with_explicit_species_mass_equals_1");
  Q043VolumeAwareRequireString(
      pin, block, "initial_eigenmode",
      uniform_current_oracle ? "uniform_zero_perturbation_parallel_stream"
                             : "section52_positive_current_unstable_eigenmode");

  Q043VolumeAwareRequireString(pin, "time", "evolution", "dynamic");
  Q043VolumeAwareRequireString(pin, "time", "integrator", "rk2");
  Q043VolumeAwareRequireString(pin, "mhd", "eos", "ideal");
  Q043VolumeAwareRequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q043VolumeAwareRequireString(pin, "particles", "pusher", "boris_tsc");
  Q043VolumeAwareRequireString(pin, "particles", "cr_distribution", "center");
  if (pin->GetInteger("particles", "nspecies") != 1) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware requires exactly one species");
  }
  Q043VolumeAwareRequireBoolean(pin, "particles", "deposit_moments", true);
  const std::string physical_mode =
      pin->GetString("particles", "pic_physical_mode");
  const bool paper_vl2_tsc =
      physical_mode.compare("paper_mhd_pic_vl2_tsc") == 0;
  if (!paper_vl2_tsc && physical_mode.compare("paper_mhd_pic") != 0) {
    Q043VolumeAwareFatal("<particles>/pic_physical_mode does not match the corrected "
                    "Q-043 Section 5.2 contract");
  }
  if (pin->GetInteger("particles", "deposit_order") != (paper_vl2_tsc ? 2 : 1)) {
    Q043VolumeAwareFatal("<particles>/deposit_order does not match the corrected "
                    "Q-043 Section 5.2 physical mode");
  }
  Q043VolumeAwareRequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q043VolumeAwareRequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q043VolumeAwareRequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", true);
  Q043VolumeAwareRequireString(pin, "particles", "couple_j_to_efield_representation",
                          "cell_centered");
  Q043VolumeAwareRequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");
  Q043VolumeAwareRequireString(pin, "particles", "couple_fluid_feedback_order",
                          "mhd_src_terms");
  Q043VolumeAwareRequireString(pin, "particles", "pic_background_mode", "coupled");
  Q043VolumeAwareRequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q043VolumeAwareRequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q043VolumeAwareRequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q043VolumeAwareRequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  const std::string hall_mode = pin->GetString("particles", "pic_cr_hall_mode");
  if (hall_mode.compare("off") != 0 && hall_mode.compare("full") != 0) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware supports only the "
                         "physical CR-Hall off and full modes");
  }
  Q043VolumeAwareRequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q043VolumeAwareRequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q043VolumeAwareRequireString(pin, "particles", "pic_expanding_box_mode", "off");

  const Real epsilon_default = pin->GetReal(block, "epsilon_default");
  const Real epsilon = pin->GetOrAddReal(block, "epsilon", epsilon_default);
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real amplitude = pin->GetReal(block, "amplitude");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real seed_wavenumber = pin->GetOrAddReal(block, "seed_wavenumber", k0);
  const Real omega = pin->GetReal(block, "omega");
  const Real c_over_v_cr = pin->GetReal(block, "c_over_v_cr");
  Q043VolumeAwareRequireFinite("epsilon_default", epsilon_default);
  Q043VolumeAwareRequireFinite("epsilon", epsilon);
  Q043VolumeAwareRequireFinite("rho", rho);
  Q043VolumeAwareRequireFinite("pressure", pressure);
  Q043VolumeAwareRequireFinite("amplitude", amplitude);
  Q043VolumeAwareRequireFinite("b_g", b_g);
  Q043VolumeAwareRequireFinite("u_a", u_a);
  Q043VolumeAwareRequireFinite("wavelength", wavelength);
  Q043VolumeAwareRequireFinite("k0", k0);
  Q043VolumeAwareRequireFinite("seed_wavenumber", seed_wavenumber);
  Q043VolumeAwareRequireFinite("omega", omega);
  Q043VolumeAwareRequireFinite("c_over_v_cr", c_over_v_cr);
  if (!(epsilon > 0.0 && epsilon < 1.0) || rho <= 0.0 || pressure <= 0.0 ||
      !SourceModeAmplitudeIsValid(source_mode, amplitude) || b_g <= 0.0 ||
      wavelength <= 0.0 || k0 <= 0.0 || seed_wavenumber <= 0.0 ||
      omega <= 0.0 || c_over_v_cr <= 1.0) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware physical normalization "
                    "contract is invalid");
  }
  Q043VolumeAwareRequireClose("rho", rho, 1.0);
  Q043VolumeAwareRequireClose("b_g", b_g, 1.0);
  Q043VolumeAwareRequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q043VolumeAwareRequireClose("wavelength", wavelength, 1.0);
  Q043VolumeAwareRequireClose("k0", k0, 2.0*M_PI/wavelength);
  Q043VolumeAwareRequireClose("omega", omega, 1.0e-6*k0*u_a);

  const Basis basis = ModeBasis(dimension);
  const Real cr_vx = pin->GetReal("particles", "cr_vx0");
  const Real cr_vy = pin->GetReal("particles", "cr_vy0");
  const Real cr_vz = pin->GetReal("particles", "cr_vz0");
  const Real v_cr = std::sqrt(cr_vx*cr_vx + cr_vy*cr_vy + cr_vz*cr_vz);
  const Real light_speed = pin->GetReal("particles", "pic_cr_light_speed");
  const Real ppc = pin->GetReal("particles", "ppc");
  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real species_mass = pin->GetReal("species0", "mass");
  const Real species_charge = pin->GetReal("species0", "charge");
  const Real species_vx = pin->GetReal("species0", "vx0");
  const Real species_vy = pin->GetReal("species0", "vy0");
  const Real species_vz = pin->GetReal("species0", "vz0");
  const Real species_v_cr =
      std::sqrt(species_vx*species_vx + species_vy*species_vy + species_vz*species_vz);
  const Real root_cell_volume = q043_bell_current_volume_aware::RootCellVolume(
      pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min,
      pmy_mesh_->mesh_indcs.nx1,
      pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min,
      pmy_mesh_->mesh_indcs.nx2,
      pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min,
      pmy_mesh_->mesh_indcs.nx3);
  Q043VolumeAwareRequireFinite("cr_vx0", cr_vx);
  Q043VolumeAwareRequireFinite("cr_vy0", cr_vy);
  Q043VolumeAwareRequireFinite("cr_vz0", cr_vz);
  Q043VolumeAwareRequireFinite("v_cr", v_cr);
  Q043VolumeAwareRequireFinite("pic_cr_light_speed", light_speed);
  Q043VolumeAwareRequireFinite("ppc", ppc);
  Q043VolumeAwareRequireFinite("deposit_qscale", qscale);
  Q043VolumeAwareRequireFinite("species mass", species_mass);
  Q043VolumeAwareRequireFinite("species charge", species_charge);
  Q043VolumeAwareRequireFinite("species0 vx0", species_vx);
  Q043VolumeAwareRequireFinite("species0 vy0", species_vy);
  Q043VolumeAwareRequireFinite("species0 vz0", species_vz);
  Q043VolumeAwareRequireFinite("species0 stream speed", species_v_cr);
  Q043VolumeAwareRequireFinite("root cell volume", root_cell_volume);
  if (v_cr <= 0.0 || !PositiveIntegralPPCIsValid(ppc) || qscale <= 0.0 ||
      species_mass <= 0.0 ||
      species_charge <= 0.0 || species_v_cr <= 0.0 || root_cell_volume <= 0.0) {
    Q043VolumeAwareFatal("q043_bell_current_volume_aware particle normalization "
                    "contract is invalid; PPC must be a positive integer because "
                    "AthenaK realizes a discrete global particle count");
  }
  Q043VolumeAwareRequireClose(
      "root cell volume", root_cell_volume,
      pmy_mesh_->mesh_size.dx1*pmy_mesh_->mesh_size.dx2*pmy_mesh_->mesh_size.dx3);
  Q043VolumeAwareRequireClose("epsilon", epsilon, u_a/v_cr);
  Q043VolumeAwareRequireClose("cr_vx0", cr_vx, v_cr*basis.parallel.x1);
  Q043VolumeAwareRequireClose("cr_vy0", cr_vy, v_cr*basis.parallel.x2);
  Q043VolumeAwareRequireClose("cr_vz0", cr_vz, v_cr*basis.parallel.x3);
  Q043VolumeAwareRequireClose("species0 vx0", species_vx, cr_vx);
  Q043VolumeAwareRequireClose("species0 vy0", species_vy, cr_vy);
  Q043VolumeAwareRequireClose("species0 vz0", species_vz, cr_vz);
  Q043VolumeAwareRequireClose("pic_cr_light_speed", light_speed, c_over_v_cr*v_cr);
  if (!SourceModeSpeciesMassIsValid(source_mode, species_mass)) {
    Q043VolumeAwareFatal("corrected_linear_eigenmode requires species_mass=1; "
                         "uniform_current_oracle accepts any positive species_mass");
  }
  if (!HasRequiredSpeciesChargeOverMass(species_mass, species_charge, omega/b_g)) {
    Q043VolumeAwareFatal("species_charge/species_mass must equal omega/b_g; "
                         "species mass must not enter the deposited-current closure");
  }
  if (!HasRequiredDepositedJOverC(
          ppc, qscale, species_charge, species_v_cr, root_cell_volume, b_g, k0)) {
    Q043VolumeAwareFatal("PPC*deposit_qscale*species_charge*v_CR/V_root_cell must equal "
                    "2*b_g*k0; artificial light speed and MeshBlock decomposition "
                    "must not enter the required deposited current");
  }

  if (restart) return;

  const ModeParameters parameters = {
    dimension, epsilon, amplitude, rho, pressure, b_g, u_a, seed_wavenumber
  };
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;

  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 2;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 2;
  DvceArray4D<Real> a1;
  DvceArray4D<Real> a2;
  DvceArray4D<Real> a3;
  Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

  par_for("pgen_q043_volume_aware_bell_vector_potential", DevExeSpace(),
          0, nmb - 1, ks, ke + 1, js, je + 1, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const int nx1 = indcs.nx1;
    const int nx2 = indcs.nx2;
    const int nx3 = indcs.nx3;
    const Real x1v = CellCenterX(i - is, nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i - is, nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j - js, nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k - ks, nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    a1(m, k, j, i) = UnstableVectorPotentialAt(parameters, {x1v, x2f, x3f}).x1;
    a2(m, k, j, i) = UnstableVectorPotentialAt(parameters, {x1f, x2v, x3f}).x2;
    a3(m, k, j, i) = UnstableVectorPotentialAt(parameters, {x1f, x2f, x3v}).x3;
  });

  par_for("pgen_q043_volume_aware_bell_ct_field", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m, k, j, i) = (a3(m, k, j + 1, i) - a3(m, k, j, i))/dx2 -
                          (a2(m, k + 1, j, i) - a2(m, k, j, i))/dx3;
    b0.x2f(m, k, j, i) = (a1(m, k + 1, j, i) - a1(m, k, j, i))/dx3 -
                          (a3(m, k, j, i + 1) - a3(m, k, j, i))/dx1;
    b0.x3f(m, k, j, i) = (a2(m, k, j, i + 1) - a2(m, k, j, i))/dx1 -
                          (a1(m, k, j + 1, i) - a1(m, k, j, i))/dx2;
    if (i == ie) {
      b0.x1f(m, k, j, i + 1) =
          (a3(m, k, j + 1, i + 1) - a3(m, k, j, i + 1))/dx2 -
          (a2(m, k + 1, j, i + 1) - a2(m, k, j, i + 1))/dx3;
    }
    if (j == je) {
      b0.x2f(m, k, j + 1, i) =
          (a1(m, k + 1, j + 1, i) - a1(m, k, j + 1, i))/dx3 -
          (a3(m, k, j + 1, i + 1) - a3(m, k, j + 1, i))/dx1;
    }
    if (k == ke) {
      b0.x3f(m, k + 1, j, i) =
          (a2(m, k + 1, j, i + 1) - a2(m, k + 1, j, i))/dx1 -
          (a1(m, k + 1, j + 1, i) - a1(m, k + 1, j, i))/dx2;
    }
  });

  par_for("pgen_q043_volume_aware_bell_primitives", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const Real phase = parameters.k0*
        (basis.parallel.x1*x1 + basis.parallel.x2*x2 + basis.parallel.x3*x3);
    const auto sample = UnstableEigenmodeAtPhase(parameters, phase);
    w0(m, IDN, k, j, i) = parameters.rho;
    w0(m, IVX, k, j, i) = sample.velocity.x1;
    w0(m, IVY, k, j, i) = sample.velocity.x2;
    w0(m, IVZ, k, j, i) = sample.velocity.x3;
    w0(m, IEN, k, j, i) = parameters.pressure;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
#endif
