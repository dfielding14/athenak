//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q007_paper_deltaf_linear.cpp
//! \brief Bounded source-local Sun and Bai true-delta-f linear preparations.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

#include "q007_paper_deltaf_linear.hpp"

namespace {

struct Q007ModeContract {
  std::string block;
  std::string campaign_id;
  std::string section_anchor;
  std::string deltaf_background;
  std::string case_role;
  Real source_local_x1max;
  Real paper_dx;
  Real light_speed;
  Real p0;
  Real kappa;
  Real xi;
  Real gas_vx;
  Real seed_amplitude;
};

[[noreturn]] void Q007Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q007RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(measured) ||
      std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q007Fatal(label + " does not match the Q-007 preparation contract");
  }
}

void Q007RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q007Fatal("<" + block + ">/" + name +
              " does not match the Q-007 preparation contract");
  }
}

void Q007RequireInteger(ParameterInput *pin, const std::string &block,
                        const std::string &name, const int expected) {
  if (pin->GetInteger(block, name) != expected) {
    Q007Fatal("<" + block + ">/" + name +
              " does not match the Q-007 preparation contract");
  }
}

void Q007RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q007Fatal("<" + block + ">/" + name +
              " does not match the Q-007 preparation contract");
  }
}

void Q007RequireReal(ParameterInput *pin, const std::string &block,
                     const std::string &name, const Real expected) {
  Q007RequireClose("<" + block + ">/" + name, pin->GetReal(block, name), expected);
}

void Q007ValidateCommon(ParameterInput *pin, const Q007ModeContract &mode) {
  Q007RequireString(pin, "time", "evolution", "dynamic");
  Q007RequireString(pin, "time", "integrator", "rk2");
  Q007RequireReal(pin, "time", "cfl_number", 0.1);
  Q007RequireInteger(pin, "time", "nlim", 0);
  Q007RequireReal(pin, "time", "tlim", 0.0);

  Q007RequireString(pin, "mhd", "eos", "isothermal");
  Q007RequireReal(pin, "mhd", "iso_sound_speed", 1.0);
  Q007RequireString(pin, "mhd", "reconstruct", "plm");
  Q007RequireString(pin, "mhd", "rsolver", "llf");

  Q007RequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q007RequireReal(pin, "particles", "ppc", 0.0625);
  Q007RequireString(pin, "particles", "pusher", "boris_tsc");
  Q007RequireInteger(pin, "particles", "nspecies", 8);
  Q007RequireString(pin, "particles", "cr_distribution", "center");
  Q007RequireBoolean(pin, "particles", "deposit_moments", true);
  Q007RequireInteger(pin, "particles", "deposit_order", 1);
  Q007RequireReal(pin, "particles", "deposit_qscale", 1.0);
  Q007RequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q007RequireReal(pin, "particles", "couple_j_to_efield_coeff", 1.0);
  Q007RequireString(pin, "particles", "couple_j_to_efield_representation",
                    "cell_centered");
  Q007RequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");
  Q007RequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q007RequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", false);
  Q007RequireString(pin, "particles", "couple_fluid_feedback_order",
                    "mhd_src_terms");
  Q007RequireString(pin, "particles", "pic_physical_mode", "paper_mhd_pic");
  Q007RequireString(pin, "particles", "pic_background_mode", "coupled");
  Q007RequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q007RequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q007RequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q007RequireReal(pin, "particles", "pic_cr_light_speed", mode.light_speed);
  Q007RequireString(pin, "particles", "pic_cr_initial_state", "momentum");
  Q007RequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q007RequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q007RequireInteger(pin, "particles", "pic_max_cell_cross", 1);
  Q007RequireReal(pin, "particles", "pic_theta_max", 0.3);
  Q007RequireString(pin, "particles", "pic_deltaf_mode", "physical");
  Q007RequireString(pin, "particles", "pic_deltaf_f0", mode.deltaf_background);
  Q007RequireReal(pin, "particles", "pic_deltaf_p0", mode.p0);
  Q007RequireReal(pin, "particles", "pic_deltaf_kappa", mode.kappa);
  Q007RequireReal(pin, "particles", "pic_deltaf_drift_x1", 0.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_drift_x2", 0.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_drift_x3", 0.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_background_rho", 1.0e-4);
  Q007RequireReal(pin, "particles", "pic_deltaf_background_jx", 0.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_background_jy", 0.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_background_jz", 0.0);
  Q007RequireString(pin, "particles", "pic_deltaf_adapt_mode", "off");
  Q007RequireString(pin, "particles", "pic_expanding_box_mode", "off");

  const Real transverse_anisotropy =
      static_cast<Real>(q007_paper_deltaf_linear::AthenaKTransverseAnisotropyScale(
          mode.xi));
  Q007RequireReal(pin, "particles", "pic_deltaf_aniso_x1", 1.0);
  Q007RequireReal(pin, "particles", "pic_deltaf_aniso_x2", transverse_anisotropy);
  Q007RequireReal(pin, "particles", "pic_deltaf_aniso_x3", transverse_anisotropy);

  for (int species = 0; species < 8; ++species) {
    const std::string block = "species" + std::to_string(species);
    Q007RequireReal(pin, block, "mass", 1.0);
    Q007RequireReal(pin, block, "charge", 1.0);
    Q007RequireReal(pin, block, "vx0", 0.0);
    Q007RequireReal(pin, block, "vy0", 0.0);
    Q007RequireReal(pin, block, "vz0", 0.0);
  }

  Q007RequireString(pin, mode.block, "campaign_id", mode.campaign_id);
  Q007RequireString(pin, mode.block, "section_anchor", mode.section_anchor);
  Q007RequireString(
      pin, mode.block, "deck_role",
      "cycle_zero_true_deltaf_parser_and_analytic_mapping_preparation_only");
  Q007RequireString(pin, mode.block, "qualification_effect", "none");
  Q007RequireString(pin, mode.block, "frontier_authorization", "not_bound");
  Q007RequireString(pin, mode.block, "runtime_evolution",
                    "blocked_cycle_zero_only");
  Q007RequireString(pin, mode.block, "physical_loading",
                    "blocked_missing_paper_log_bin_weighted_loader");
  Q007RequireString(pin, mode.block, "initial_wave_spectrum",
                    "blocked_missing_random_phase_four_branch_loader");
  Q007RequireString(pin, mode.block, "theory_runtime_comparison",
                    "blocked_missing_runtime_artifacts_and_independent_review");
  Q007RequireInteger(pin, mode.block, "paper_momentum_bin_count", 8);
  Q007RequireInteger(pin, mode.block, "paper_particles_per_cell_per_bin", 256);
  Q007RequireInteger(pin, mode.block, "source_local_placeholder_particles_total", 8);
  Q007RequireReal(pin, mode.block, "paper_rho0", 1.0);
  Q007RequireReal(pin, mode.block, "paper_b0", 1.0);
  Q007RequireReal(pin, mode.block, "paper_ua", 1.0);
  Q007RequireReal(pin, mode.block, "paper_mncr_over_rho0", 1.0e-4);
  Q007RequireReal(pin, mode.block, "paper_domain_x1", 96000.0);
  Q007RequireReal(pin, mode.block, "paper_dx", mode.paper_dx);
  Q007RequireReal(pin, mode.block, "source_local_x1", mode.source_local_x1max);
  Q007RequireReal(pin, mode.block, "paper_light_speed", mode.light_speed);
  Q007RequireReal(pin, mode.block, "paper_p0", mode.p0);
  Q007RequireReal(pin, mode.block, "paper_kappa", mode.kappa);
  Q007RequireReal(pin, mode.block, "paper_xi", mode.xi);
  Q007RequireReal(pin, mode.block, "paper_seed_amplitude", mode.seed_amplitude);
  Q007RequireString(pin, mode.block, "case_role", mode.case_role);
  if (mode.deltaf_background.compare("kappa_iso") == 0) {
    Q007RequireReal(pin, mode.block, "paper_vd", 2.0);
    Q007RequireReal(pin, mode.block, "source_local_gas_vx", -2.0);
    Q007RequireString(pin, mode.block, "handedness_mapping",
                      "paper_both_forward_polarizations_static_mapping_only");
  } else {
    Q007RequireReal(pin, mode.block, "source_local_gas_vx", 0.0);
    Q007RequireString(pin, mode.block, "handedness_mapping",
                      "blocked_pending_manuscript_text_caption_review");
  }
}

void Q007Prepare(Mesh *pmesh, ParameterInput *pin, const bool restart,
                 const Q007ModeContract &mode) {
  MeshBlockPack *pmbp = pmesh->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q007Fatal("Q-007 true-delta-f preparation requires MHD and particles");
  }
  if (pmbp->pmhd->peos->eos_data.is_ideal) {
    Q007Fatal("Q-007 true-delta-f preparation requires exact isothermal MHD");
  }
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic ||
      pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
      pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
      pin->DoesBlockExist("z4c") || pin->DoesBlockExist("turb_driving") ||
      pin->DoesBlockExist("initial_turb")) {
    Q007Fatal("Q-007 true-delta-f preparation requires the exact Newtonian "
              "single-fluid isothermal-MHD task path");
  }
  if (pmesh->multilevel || pmesh->adaptive) {
    Q007Fatal("Q-007 true-delta-f preparation rejects AMR and SMR");
  }
  if (!pmesh->two_d || pmesh->one_d || pmesh->three_d ||
      pmesh->mesh_indcs.nx1 != 32 || pmesh->mesh_indcs.nx2 != 4 ||
      pmesh->mesh_indcs.nx3 != 1 || pmesh->mb_indcs.nx1 != 32 ||
      pmesh->mb_indcs.nx2 != 4 || pmesh->mb_indcs.nx3 != 1 ||
      pmbp->nmb_thispack != 1 || pmbp->ppart->nprtcl_thispack != 8) {
    Q007Fatal("Q-007 true-delta-f preparation requires its exact serial "
              "32x4x1 thin-2D3V one-MeshBlock eight-placeholder carrier");
  }
  if (!pmesh->strictly_periodic) {
    Q007Fatal("Q-007 true-delta-f preparation requires periodic boundaries");
  }
  Q007RequireClose("global x1min", pmesh->mesh_size.x1min, 0.0);
  Q007RequireClose("global x1max", pmesh->mesh_size.x1max,
                   mode.source_local_x1max);
  Q007RequireClose("global x2min", pmesh->mesh_size.x2min, 0.0);
  Q007RequireClose("global x2max", pmesh->mesh_size.x2max, 4.0);
  Q007RequireClose("global x3min", pmesh->mesh_size.x3min, 0.0);
  Q007RequireClose("global x3max", pmesh->mesh_size.x3max, 1.0);

  Q007ValidateCommon(pin, mode);
  if (restart) return;

  auto &indcs = pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  const Real gas_vx = mode.gas_vx;

  par_for("pgen_q007_uniform_isothermal_paper_state", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m, k, j, i) = 1.0;
    b0.x2f(m, k, j, i) = 0.0;
    b0.x3f(m, k, j, i) = 0.0;
    if (i == ie) b0.x1f(m, k, j, i + 1) = 1.0;
    if (j == je) b0.x2f(m, k, j + 1, i) = 0.0;
    if (k == ke) b0.x3f(m, k + 1, j, i) = 0.0;
    w0(m, IDN, k, j, i) = 1.0;
    w0(m, IVX, k, j, i) = gas_vx;
    w0(m, IVY, k, j, i) = 0.0;
    w0(m, IVZ, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) = 1.0;
    bcc0(m, IBY, k, j, i) = 0.0;
    bcc0(m, IBZ, k, j, i) = 0.0;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}

}  // namespace

void ProblemGenerator::Q007PaperCRSILinearPreparation(ParameterInput *pin,
                                                       const bool restart) {
  const Q007ModeContract mode = {
    "q007_paper_crsi_linear_preparation",
    "Q007-PAPER-CRSI-LINEAR-PREPARATION",
    "sun_bai_2023_section_5_5_1_crsi",
    "kappa_iso",
    "crsi_isotropic_kappa_gas_drift_minus_vd",
    320.0, 10.0, 300.0, 300.0, 1.25, 1.0, -2.0, 1.0e-3
  };
  Q007Prepare(pmy_mesh_, pin, restart, mode);
}

void ProblemGenerator::Q007PaperCRPAILinearPreparation(ParameterInput *pin,
                                                        const bool restart) {
  const std::string block = "q007_paper_crpai_linear_preparation";
  const Real xi = pin->GetReal(block, "paper_xi");
  Q007ModeContract mode = {
    block,
    "Q007-PAPER-CRPAI-LINEAR-PREPARATION",
    "sun_bai_2023_section_5_5_2_crpai",
    "kappa_aniso",
    "",
    640.0, 20.0, 3.0e4, 300.0, 1.75, xi, 0.0, 1.0e-3
  };
  if (std::abs(xi - static_cast<Real>(0.99)) < static_cast<Real>(1.0e-12)) {
    mode.case_role = "crpai_prolate_xi_0p99_signed_branch_mapping_only";
  } else if (std::abs(xi - static_cast<Real>(1.01)) <
             static_cast<Real>(1.0e-12)) {
    mode.case_role = "crpai_oblate_xi_1p01_signed_branch_mapping_only";
  } else {
    Q007Fatal("<q007_paper_crpai_linear_preparation>/paper_xi must be "
              "exactly 0.99 or 1.01");
  }
  Q007Prepare(pmy_mesh_, pin, restart, mode);
}
