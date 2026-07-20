//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q023_paper_bell_linear_joverc.cpp
//! \brief Corrected volume-aware Section 5.2 Bell linear predecessor.

#include <algorithm>
#include <cmath>

#if !defined(Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT)
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

#if defined(Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT) && \
    !defined(Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#define Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#define Q023_JOVERC_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT 1
#endif
#include "q023_paper_bell_linear_joverc.hpp"
#if defined(Q023_JOVERC_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT)
#undef Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#undef Q023_JOVERC_UNDEF_Q023_PAPER_BELL_LINEAR_HOST_CONTRACT
#endif

namespace q023_paper_bell_linear_joverc {

using q023_paper_bell_linear::Basis;
using q023_paper_bell_linear::ModeBasis;
using q023_paper_bell_linear::ModeParameters;
using q023_paper_bell_linear::Vector3;

inline double RootCellVolume(const double x1_extent, const int root_nx1,
                             const double x2_extent, const int root_nx2,
                             const double x3_extent, const int root_nx3) {
  return (x1_extent/static_cast<double>(root_nx1))*
      (x2_extent/static_cast<double>(root_nx2))*
      (x3_extent/static_cast<double>(root_nx3));
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

inline bool HasPositiveIntegralPPC(const double ppc) {
  return std::isfinite(ppc) && ppc > 0.0 && std::floor(ppc) == ppc;
}

inline bool HasMatchingVelocityVectors(const double global_vx, const double global_vy,
                                       const double global_vz, const double species_vx,
                                       const double species_vy,
                                       const double species_vz) {
  const double global[3] = {global_vx, global_vy, global_vz};
  const double species[3] = {species_vx, species_vy, species_vz};
  for (int component = 0; component < 3; ++component) {
    if (!std::isfinite(global[component]) || !std::isfinite(species[component])) {
      return false;
    }
    const double scale = std::max(1.0, std::abs(species[component]));
    if (std::abs(global[component] - species[component]) > 1.0e-12*scale) {
      return false;
    }
  }
  return true;
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

}  // namespace q023_paper_bell_linear_joverc

#if !defined(Q023_PAPER_BELL_LINEAR_JOVERC_HOST_CONTRACT)
namespace {

[[noreturn]] void Q023JOverCFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q023JOverCRequireClose(const std::string &label, const Real measured,
                            const Real expected) {
  if (!std::isfinite(measured) || !std::isfinite(expected)) {
    Q023JOverCFatal(label + " must be finite for the corrected Q-023-JOVERC "
                   "Section 5.2 contract");
  }
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q023JOverCFatal(label + " does not match the corrected Q-023-JOVERC Section 5.2 "
                   "contract");
  }
}

void Q023JOverCRequireFinite(const std::string &label, const Real value) {
  if (!std::isfinite(value)) {
    Q023JOverCFatal(label + " must be finite for the corrected Q-023-JOVERC "
                   "Section 5.2 contract");
  }
}

void Q023JOverCRequireString(ParameterInput *pin, const std::string &block,
                             const std::string &name,
                             const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q023JOverCFatal("<" + block + ">/" + name +
                    " does not match the corrected Q-023-JOVERC Section 5.2 contract");
  }
}

void Q023JOverCRequireBoolean(ParameterInput *pin, const std::string &block,
                              const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q023JOverCFatal("<" + block + ">/" + name +
                    " does not match the corrected Q-023-JOVERC Section 5.2 contract");
  }
}

}  // namespace

void ProblemGenerator::Q023PaperBellLinearJOverC(ParameterInput *pin,
                                                  const bool restart) {
  using q023_paper_bell_linear_joverc::Basis;
  using q023_paper_bell_linear_joverc::HasRequiredSpeciesChargeOverMass;
  using q023_paper_bell_linear_joverc::HasRequiredDepositedJOverC;
  using q023_paper_bell_linear_joverc::HasMatchingVelocityVectors;
  using q023_paper_bell_linear_joverc::HasPositiveIntegralPPC;
  using q023_paper_bell_linear_joverc::ModeBasis;
  using q023_paper_bell_linear_joverc::ModeParameters;
  using q023_paper_bell_linear_joverc::UnstableEigenmodeAtPhase;
  using q023_paper_bell_linear_joverc::UnstableVectorPotentialAt;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc requires ideal-MHD EOS");
  }
  if (pmy_mesh_->multilevel) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc rejects AMR/SMR");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc requires periodic boundaries");
  }

  const std::string block = "q023_paper_bell_linear_joverc";
  const int dimension = pin->GetInteger(block, "dimension");
  const int mesh_dimension = pmy_mesh_->one_d ? 1 : (pmy_mesh_->two_d ? 2 : 3);
  const bool thin_2d3v_1d_carrier =
      dimension == 1 && mesh_dimension == 2 &&
      pmy_mesh_->mb_indcs.nx2 == 4 && pmy_mesh_->mb_indcs.nx3 == 1;
  const bool native_dimension =
      dimension >= 2 && dimension <= 3 && dimension == mesh_dimension;
  if (!(thin_2d3v_1d_carrier || native_dimension)) {
    Q023JOverCFatal("<q023_paper_bell_linear_joverc>/dimension does not match "
                    "the required native mesh or exact thin 2D3V carrier");
  }

  Q023JOverCRequireString(pin, block, "campaign_id",
                          "Q023-PAPER-BELL-LINEAR-JOVERC");
  Q023JOverCRequireString(pin, block, "supersedes_campaign_id",
                          "Q023-PAPER-BELL-LINEAR");
  Q023JOverCRequireString(pin, block, "foundational_current_campaign_id",
                          "Q043-BELL-CURRENT-VOLUME-AWARE");
  Q023JOverCRequireString(pin, block, "foundational_raw_oracle_id",
                          "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE");
  Q023JOverCRequireString(
      pin, block, "foundational_binding_status",
      "source_deck_raw_oracle_sequence_bound_registered_admission_hardening_pending");
  Q023JOverCRequireString(
      pin, block, "foundational_current_lineage_commit",
      "883e4679ae797392f25f86e46b952feef753e6aa");
  Q023JOverCRequireString(
      pin, block, "foundational_current_source_sha256",
      "dd44fff4fe39bf1640dad6019ded29ce7d434b50a0c9831a8749383c6241f3e3");
  Q023JOverCRequireString(
      pin, block, "foundational_raw_oracle_analyzer_sha256",
      "f0a8db7fb47798d4bdc9ec742cef09e5ae0c69ff3ccb2f25d065602b23a26bed");
  Q023JOverCRequireString(
      pin, block, "foundational_raw_oracle_deck_manifest_sha256",
      "77d6ec862314d308c41490dde1137f6c384d1da24b6b72ad43e357545cff98b0");
  Q023JOverCRequireString(
      pin, block, "foundational_supersession_sha256",
      "2f83afc03e83038db168fed7d40bdc85bc7bca1eccbba52d18ede0d319dce342");
  Q023JOverCRequireString(
      pin, block, "foundational_registered_admission_binding_status",
      "pending_hardened_successor_digest_and_schema");
  Q023JOverCRequireString(
      pin, block, "foundational_registered_admission_binding_permitted",
      "false_until_hardened_successor_lands");
  Q023JOverCRequireString(pin, block, "foundational_raw_oracle_case_count", "132");
  Q023JOverCRequireString(
      pin, block, "foundational_raw_oracle_registered_pass",
      "required_before_linear_qualification");
  Q023JOverCRequireString(
      pin, block, "foundational_registered_execution_admission",
      "required_before_any_q023_execution_or_qualification");
  Q023JOverCRequireString(
      pin, block, "dependency_effect",
      "non_authorizing_downstream_predecessor_dependency_only");
  Q023JOverCRequireString(
      pin, block, "deck_role",
      "non_authorizing_q023_joverc_linear_predecessor");
  Q023JOverCRequireString(pin, block, "current_normalization",
                          "deposited_j_over_c_equals_2_b_g_k0");
  Q023JOverCRequireString(
      pin, block, "deposition_measure",
      "ppc_times_deposit_qscale_times_species_charge_times_species0_actual_v_cr_"
      "over_root_cell_volume");
  Q023JOverCRequireString(
      pin, block, "particle_velocity_semantics",
      "species0_vx0_vy0_vz0_are_actual_and_must_match_particles_cr_velocity");
  Q023JOverCRequireString(
      pin, block, "root_cell_volume",
      "global_root_mesh_extents_over_global_root_mesh_counts");
  Q023JOverCRequireString(
      pin, block, "deposit_qscale_semantics",
      "root_cell_macro_charge_volume_aware");
  Q023JOverCRequireString(
      pin, block, "fixed_qscale",
      "forbidden_across_dimensions_and_resolutions");
  Q023JOverCRequireString(
      pin, block, "decomposition_contract",
      "dimension_appropriate_multidirectional_matrix_required");
  Q023JOverCRequireString(pin, block, "raw_output_layout", "shared_mpi_io");
  Q023JOverCRequireString(pin, block, "epsilon_grid", "0.1,0.2,0.4,0.6,0.8");
  Q023JOverCRequireString(pin, block, "timestep",
                          "open_clean_candidate_timestep_freeze");
  Q023JOverCRequireString(pin, block, "source_mode", "corrected_linear_eigenmode");
  Q023JOverCRequireString(
      pin, block, "q_over_mc_representation",
      "species_charge_over_species_mass_equals_omega_over_b_g");
  Q023JOverCRequireString(
      pin, block, "initial_eigenmode", "section52_right_polarized_eigenmode");
  Q023JOverCRequireString(
      pin, block, "qualification_effect",
      "none_non_authorizing_predecessor_preparation_only");
  Q023JOverCRequireBoolean(pin, block, "launch_authorized", false);
  Q023JOverCRequireBoolean(pin, block, "qualification_eligible", false);

  Q023JOverCRequireString(pin, "time", "evolution", "dynamic");
  Q023JOverCRequireString(pin, "time", "integrator", "vl2");
  Q023JOverCRequireString(pin, "mhd", "eos", "ideal");
  Q023JOverCRequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q023JOverCRequireString(pin, "particles", "pusher", "boris_tsc");
  Q023JOverCRequireString(pin, "particles", "cr_distribution", "center");
  if (pin->GetInteger("particles", "nspecies") != 1) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc requires exactly one species");
  }
  Q023JOverCRequireBoolean(pin, "particles", "deposit_moments", true);
  if (pin->GetInteger("particles", "deposit_order") != 2) {
    Q023JOverCFatal("Q-023-JOVERC requires TSC moment deposition");
  }
  Q023JOverCRequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q023JOverCRequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q023JOverCRequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", true);
  Q023JOverCRequireClose(
      "<particles>/couple_j_to_efield_coeff",
      pin->GetReal("particles", "couple_j_to_efield_coeff"), 1.0);
  Q023JOverCRequireClose(
      "<particles>/couple_moments_momentum_coeff",
      pin->GetReal("particles", "couple_moments_momentum_coeff"), 1.0);
  Q023JOverCRequireClose(
      "<particles>/couple_moments_energy_coeff",
      pin->GetReal("particles", "couple_moments_energy_coeff"), 1.0);
  Q023JOverCRequireString(pin, "particles", "couple_j_to_efield_representation",
                          "cell_centered");
  Q023JOverCRequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");
  Q023JOverCRequireString(pin, "particles", "couple_fluid_feedback_order",
                          "mhd_src_terms");
  Q023JOverCRequireString(pin, "particles", "pic_background_mode", "coupled");
  Q023JOverCRequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q023JOverCRequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q023JOverCRequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q023JOverCRequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  Q023JOverCRequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q023JOverCRequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q023JOverCRequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q023JOverCRequireString(pin, "particles", "pic_expanding_box_mode", "off");

  const Real epsilon_default = pin->GetReal(block, "epsilon_default");
  const Real epsilon = pin->GetOrAddReal(block, "epsilon", epsilon_default);
  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real amplitude = pin->GetReal(block, "amplitude");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real omega = pin->GetReal(block, "omega");
  const Real c_over_v_cr = pin->GetReal(block, "c_over_v_cr");
  Q023JOverCRequireFinite("epsilon_default", epsilon_default);
  Q023JOverCRequireFinite("epsilon", epsilon);
  Q023JOverCRequireFinite("rho", rho);
  Q023JOverCRequireFinite("pressure", pressure);
  Q023JOverCRequireFinite("amplitude", amplitude);
  Q023JOverCRequireFinite("b_g", b_g);
  Q023JOverCRequireFinite("u_a", u_a);
  Q023JOverCRequireFinite("wavelength", wavelength);
  Q023JOverCRequireFinite("k0", k0);
  Q023JOverCRequireFinite("omega", omega);
  Q023JOverCRequireFinite("c_over_v_cr", c_over_v_cr);
  if (!(epsilon > 0.0 && epsilon < 1.0) || rho <= 0.0 || pressure <= 0.0 ||
      amplitude <= 0.0 || amplitude > 1.0e-2 || b_g <= 0.0 ||
      wavelength <= 0.0 || k0 <= 0.0 || omega <= 0.0 || c_over_v_cr <= 1.0) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc physical normalization "
                    "contract is invalid");
  }
  Q023JOverCRequireClose("rho", rho, 1.0);
  Q023JOverCRequireClose("b_g", b_g, 1.0);
  Q023JOverCRequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q023JOverCRequireClose("wavelength", wavelength, 1.0);
  Q023JOverCRequireClose("k0", k0, 2.0*M_PI/wavelength);
  Q023JOverCRequireClose("omega", omega, 1.0e-6*k0*u_a);

  const Basis basis = ModeBasis(dimension);
  const Real cr_vx = pin->GetReal("particles", "cr_vx0");
  const Real cr_vy = pin->GetReal("particles", "cr_vy0");
  const Real cr_vz = pin->GetReal("particles", "cr_vz0");
  const Real species_vx = pin->GetReal("species0", "vx0");
  const Real species_vy = pin->GetReal("species0", "vy0");
  const Real species_vz = pin->GetReal("species0", "vz0");
  const Real v_cr = std::sqrt(
      species_vx*species_vx + species_vy*species_vy + species_vz*species_vz);
  const Real light_speed = pin->GetReal("particles", "pic_cr_light_speed");
  const Real ppc = pin->GetReal("particles", "ppc");
  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real species_mass = pin->GetReal("species0", "mass");
  const Real species_charge = pin->GetReal("species0", "charge");
  const Real root_cell_volume = q023_paper_bell_linear_joverc::RootCellVolume(
      pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min,
      pmy_mesh_->mesh_indcs.nx1,
      pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min,
      pmy_mesh_->mesh_indcs.nx2,
      pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min,
      pmy_mesh_->mesh_indcs.nx3);
  Q023JOverCRequireFinite("cr_vx0", cr_vx);
  Q023JOverCRequireFinite("cr_vy0", cr_vy);
  Q023JOverCRequireFinite("cr_vz0", cr_vz);
  Q023JOverCRequireFinite("species0/vx0", species_vx);
  Q023JOverCRequireFinite("species0/vy0", species_vy);
  Q023JOverCRequireFinite("species0/vz0", species_vz);
  Q023JOverCRequireFinite("v_cr", v_cr);
  Q023JOverCRequireFinite("pic_cr_light_speed", light_speed);
  Q023JOverCRequireFinite("ppc", ppc);
  Q023JOverCRequireFinite("deposit_qscale", qscale);
  Q023JOverCRequireFinite("species mass", species_mass);
  Q023JOverCRequireFinite("species charge", species_charge);
  Q023JOverCRequireFinite("root cell volume", root_cell_volume);
  if (v_cr <= 0.0) {
    Q023JOverCFatal("<species0>/vx0,vy0,vz0 are the actual particle initialization "
                    "velocity and must define a nonzero finite stream");
  }
  if (!HasPositiveIntegralPPC(ppc)) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc particle normalization "
                    "contract is invalid; PPC must be a positive integer because "
                    "the deposited-current closure uses the discrete particle count");
  }
  if (qscale <= 0.0 || species_mass <= 0.0 || species_charge <= 0.0 ||
      root_cell_volume <= 0.0) {
    Q023JOverCFatal("q023_paper_bell_linear_joverc particle normalization "
                    "contract requires positive qscale, species mass and charge, "
                    "and root-cell volume");
  }
  Q023JOverCRequireClose(
      "root cell volume", root_cell_volume,
      pmy_mesh_->mesh_size.dx1*pmy_mesh_->mesh_size.dx2*pmy_mesh_->mesh_size.dx3);
  if (!HasMatchingVelocityVectors(
          cr_vx, cr_vy, cr_vz, species_vx, species_vy, species_vz)) {
    Q023JOverCFatal("<species0>/vx0,vy0,vz0 are the actual particle initialization "
                    "velocity and must match <particles>/cr_vx0,cr_vy0,cr_vz0");
  }
  Q023JOverCRequireClose("epsilon", epsilon, u_a/v_cr);
  Q023JOverCRequireClose("species0/vx0", species_vx, v_cr*basis.parallel.x1);
  Q023JOverCRequireClose("species0/vy0", species_vy, v_cr*basis.parallel.x2);
  Q023JOverCRequireClose("species0/vz0", species_vz, v_cr*basis.parallel.x3);
  Q023JOverCRequireClose("pic_cr_light_speed", light_speed, c_over_v_cr*v_cr);
  if (!HasRequiredSpeciesChargeOverMass(species_mass, species_charge, omega/b_g)) {
    Q023JOverCFatal("species_charge/species_mass must equal omega/b_g; "
                         "species mass must not enter the deposited-current closure");
  }
  if (!HasRequiredDepositedJOverC(
          ppc, qscale, species_charge, v_cr, root_cell_volume, b_g, k0)) {
    Q023JOverCFatal("PPC*deposit_qscale*species_charge*v_CR/V_root_cell must equal "
                    "2*b_g*k0; artificial light speed and MeshBlock decomposition "
                    "must not enter the required deposited current");
  }

  if (restart) return;

  const ModeParameters parameters = {
    dimension, epsilon, amplitude, rho, pressure, b_g, u_a, k0
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

  par_for("pgen_q023_joverc_bell_vector_potential", DevExeSpace(),
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

  par_for("pgen_q023_joverc_bell_ct_field", DevExeSpace(),
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

  par_for("pgen_q023_joverc_bell_primitives", DevExeSpace(),
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
