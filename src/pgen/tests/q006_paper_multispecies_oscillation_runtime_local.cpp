//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q006_paper_multispecies_oscillation_runtime_local.cpp
//! \brief Bounded serial-host exact-isothermal Q-006 Section 5.3 mechanics successor.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

#include "src/athena.hpp"
#include "src/globals.hpp"
#include "src/parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

int q006_runtime_local_amr_seed = 60053;

[[noreturn]] void Q006RuntimeLocalFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q006RuntimeLocalRequireClose(const std::string &label, const Real measured,
                                  const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(measured) ||
      std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q006RuntimeLocalFatal(
        label + " does not match the bounded Q-006 runtime-local contract");
  }
}

void Q006RuntimeLocalRequireString(ParameterInput *pin, const std::string &block,
                                   const std::string &name,
                                   const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q006RuntimeLocalFatal("<" + block + ">/" + name +
                          " does not match the bounded Q-006 runtime-local contract");
  }
}

void Q006RuntimeLocalRequireInteger(ParameterInput *pin, const std::string &block,
                                    const std::string &name, const int expected) {
  if (pin->GetInteger(block, name) != expected) {
    Q006RuntimeLocalFatal("<" + block + ">/" + name +
                          " does not match the bounded Q-006 runtime-local contract");
  }
}

void Q006RuntimeLocalRequireBoolean(ParameterInput *pin, const std::string &block,
                                    const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q006RuntimeLocalFatal("<" + block + ">/" + name +
                          " does not match the bounded Q-006 runtime-local contract");
  }
}

void Q006RuntimeLocalRequireReal(ParameterInput *pin, const std::string &block,
                                const std::string &name, const Real expected) {
  Q006RuntimeLocalRequireClose("<" + block + ">/" + name,
                              pin->GetReal(block, name), expected);
}

std::uint64_t Q006RuntimeLocalSplitMix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30))*0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27))*0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

Real Q006RuntimeLocalAMRUniform01(const int seed, const int cycle, const int gid) {
  std::uint64_t key = static_cast<std::uint64_t>(seed);
  key ^= (static_cast<std::uint64_t>(cycle + 1))*0xbf58476d1ce4e5b9ULL;
  key ^= (static_cast<std::uint64_t>(gid + 1))*0xd2b74407b1ce6e93ULL;
  const std::uint64_t bits = Q006RuntimeLocalSplitMix64(key);
  return static_cast<Real>(bits >> 11)*static_cast<Real>(1.0/9007199254740992.0);
}

void Q006RuntimeLocalAuditedAMRRefinement(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  if (!pmesh->adaptive || pmesh->pmr == nullptr) {
    Q006RuntimeLocalFatal("Q-006 runtime-local audited AMR callback requires AMR");
  }

  auto &refine_flag = pmesh->pmr->refine_flag;
  const int nmb = pmbp->nmb_thispack;
  const int mbs = pmesh->gids_eachrank[global_variable::my_rank];
  const int cycle = pmesh->ncycle;
  const Real refine_probability = static_cast<Real>(0.10);
  const Real derefine_probability = static_cast<Real>(0.60);
  for (int m = 0; m < nmb; ++m) {
    const Real draw = Q006RuntimeLocalAMRUniform01(
        q006_runtime_local_amr_seed, cycle, mbs + m);
    refine_flag.h_view(mbs + m) =
        (draw < refine_probability) ? 1 :
        ((draw < refine_probability + derefine_probability) ? -1 : 0);
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

}  // namespace

void ProblemGenerator::Q006PaperMultispeciesOscillationRuntimeLocal(
    ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (global_variable::nranks != 1) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor is bounded to one serial rank");
  }
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor requires MHD and particles");
  }
  if (pmbp->pmhd->peos->eos_data.is_ideal) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor requires exact isothermal MHD");
  }
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic ||
      pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
      pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
      pin->DoesBlockExist("z4c") || pin->DoesBlockExist("turb_driving") ||
      pin->DoesBlockExist("initial_turb")) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor requires the Newtonian "
                          "single-fluid isothermal-MHD task path");
  }
  if (!pmy_mesh_->three_d || pmy_mesh_->one_d || pmy_mesh_->two_d ||
      pmy_mesh_->mesh_indcs.nx1 != 16 || pmy_mesh_->mesh_indcs.nx2 != 8 ||
      pmy_mesh_->mesh_indcs.nx3 != 8 || pmy_mesh_->mb_indcs.nx1 != 4 ||
      pmy_mesh_->mb_indcs.nx2 != 4 || pmy_mesh_->mb_indcs.nx3 != 4) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor requires its exact 16x8x8 "
                          "periodic 3D carrier with 4x4x4 MeshBlocks");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q006RuntimeLocalFatal("Q-006 runtime-local successor requires periodic boundaries");
  }
  Q006RuntimeLocalRequireClose("global x1min", pmy_mesh_->mesh_size.x1min, 0.0);
  Q006RuntimeLocalRequireClose("global x1max", pmy_mesh_->mesh_size.x1max, 16.0);
  Q006RuntimeLocalRequireClose("global x2min", pmy_mesh_->mesh_size.x2min, 0.0);
  Q006RuntimeLocalRequireClose("global x2max", pmy_mesh_->mesh_size.x2max, 8.0);
  Q006RuntimeLocalRequireClose("global x3min", pmy_mesh_->mesh_size.x3min, 0.0);
  Q006RuntimeLocalRequireClose("global x3max", pmy_mesh_->mesh_size.x3max, 8.0);

  Q006RuntimeLocalRequireString(pin, "time", "integrator", "rk2");
  Q006RuntimeLocalRequireReal(pin, "time", "cfl_number", 0.1);
  Q006RuntimeLocalRequireString(pin, "mhd", "eos", "isothermal");
  Q006RuntimeLocalRequireReal(pin, "mhd", "iso_sound_speed", 1.0);
  Q006RuntimeLocalRequireString(pin, "mhd", "reconstruct", "plm");
  Q006RuntimeLocalRequireString(pin, "mhd", "rsolver", "llf");

  Q006RuntimeLocalRequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q006RuntimeLocalRequireReal(pin, "particles", "ppc", 128.0);
  Q006RuntimeLocalRequireString(pin, "particles", "pusher", "boris_tsc");
  Q006RuntimeLocalRequireInteger(pin, "particles", "nspecies", 2);
  Q006RuntimeLocalRequireString(pin, "particles", "cr_distribution", "center");
  Q006RuntimeLocalRequireBoolean(pin, "particles", "deposit_moments", true);
  Q006RuntimeLocalRequireInteger(pin, "particles", "deposit_order", 2);
  Q006RuntimeLocalRequireReal(pin, "particles", "deposit_qscale", 0.0234375);
  Q006RuntimeLocalRequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q006RuntimeLocalRequireReal(pin, "particles", "couple_j_to_efield_coeff", 1.0);
  Q006RuntimeLocalRequireString(pin, "particles", "couple_j_to_efield_representation",
                                "cell_centered");
  Q006RuntimeLocalRequireString(pin, "particles", "couple_j_deposition_mode",
                                "cc_convert");
  Q006RuntimeLocalRequireBoolean(pin, "particles",
                                 "couple_moments_momentum_to_mhd", true);
  Q006RuntimeLocalRequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", false);
  Q006RuntimeLocalRequireString(pin, "particles", "couple_fluid_feedback_order",
                                "mhd_src_terms");
  Q006RuntimeLocalRequireReal(pin, "particles", "couple_moments_momentum_coeff", 1.0);
  Q006RuntimeLocalRequireReal(pin, "particles", "couple_moments_energy_coeff", 0.0);
  Q006RuntimeLocalRequireString(pin, "particles", "pic_physical_mode", "paper_mhd_pic");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_background_mode", "coupled");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q006RuntimeLocalRequireReal(pin, "particles", "pic_cr_light_speed", 1000.0);
  Q006RuntimeLocalRequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q006RuntimeLocalRequireInteger(pin, "particles", "pic_max_cell_cross", 1);
  Q006RuntimeLocalRequireReal(pin, "particles", "pic_theta_max", 0.1);
  Q006RuntimeLocalRequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q006RuntimeLocalRequireString(pin, "particles", "pic_expanding_box_mode", "off");
  Q006RuntimeLocalRequireReal(pin, "species0", "mass", 1.0);
  Q006RuntimeLocalRequireReal(pin, "species0", "charge", -1.0);
  Q006RuntimeLocalRequireReal(pin, "species0", "vx0", 0.0);
  Q006RuntimeLocalRequireReal(pin, "species0", "vy0", 0.1);
  Q006RuntimeLocalRequireReal(pin, "species0", "vz0", 0.0);
  Q006RuntimeLocalRequireReal(pin, "species1", "mass", 1.0);
  Q006RuntimeLocalRequireReal(pin, "species1", "charge", 1.0);
  Q006RuntimeLocalRequireReal(pin, "species1", "vx0", 0.0);
  Q006RuntimeLocalRequireReal(pin, "species1", "vy0", 0.1);
  Q006RuntimeLocalRequireReal(pin, "species1", "vz0", 0.0);

  const std::string block = "q006_paper_multispecies_oscillation_runtime_local";
  Q006RuntimeLocalRequireString(pin, block, "campaign_id",
                                "Q006-PAPER-MULTISPECIES-OSCILLATION-RUNTIME-LOCAL");
  Q006RuntimeLocalRequireString(pin, block, "qualification_effect", "none");
  Q006RuntimeLocalRequireString(pin, block, "frontier_authorization", "not_bound");
  Q006RuntimeLocalRequireString(pin, block, "runtime_scope",
                                "bounded_serial_host_mechanics_only");
  Q006RuntimeLocalRequireString(pin, block, "long_horizon_qualification", "not_claimed");
  Q006RuntimeLocalRequireString(pin, block, "true_amr_policy_qualification",
                                "not_claimed");
  Q006RuntimeLocalRequireString(pin, block, "mpi_qualification", "not_claimed");
  Q006RuntimeLocalRequireString(pin, block, "gpu_qualification", "not_claimed");
  Q006RuntimeLocalRequireString(pin, block, "external_review", "not_claimed");
  Q006RuntimeLocalRequireString(pin, block, "runtime_eos_contract",
                                "exact_isothermal_cs_1_fullf_momentum_only");
  Q006RuntimeLocalRequireString(
      pin, block, "ppc_semantics",
      "aggregate_128_round_robin_yields_64_per_species_per_cell");
  Q006RuntimeLocalRequireReal(pin, block, "rho", 1.0);
  Q006RuntimeLocalRequireReal(pin, block, "paper_cs", 1.0);
  Q006RuntimeLocalRequireReal(pin, block, "b_g", 1.0);
  Q006RuntimeLocalRequireReal(pin, block, "omega", 1.0);
  Q006RuntimeLocalRequireReal(pin, block, "species_mass_density_ratio", 1.5);
  Q006RuntimeLocalRequireReal(pin, block, "gas_uy", -0.3);
  Q006RuntimeLocalRequireReal(pin, block, "species_vy", 0.1);
  Q006RuntimeLocalRequireReal(pin, block, "artificial_c", 1000.0);
  Q006RuntimeLocalRequireReal(pin, block, "dt_omega_target", 0.1);

  const std::string grid_setup = pin->GetString(block, "grid_setup");
  if (grid_setup.compare("uniform") == 0) {
    Q006RuntimeLocalRequireString(pin, block, "deck_role",
                                  "bounded_serial_host_uniform_mechanics_only");
    Q006RuntimeLocalRequireString(pin, block, "amr_policy", "not_applicable_uniform");
    if (pmy_mesh_->multilevel || pmy_mesh_->adaptive) {
      Q006RuntimeLocalFatal(
          "Q-006 runtime-local uniform carrier rejects mesh refinement");
    }
  } else if (grid_setup.compare("smr") == 0) {
    Q006RuntimeLocalRequireString(pin, block, "deck_role",
                                  "bounded_serial_host_smr_mechanics_only");
    Q006RuntimeLocalRequireString(pin, block, "amr_policy",
                                  "paper_static_one_eighth_region");
    if (!pmy_mesh_->multilevel || pmy_mesh_->adaptive) {
      Q006RuntimeLocalFatal("Q-006 runtime-local SMR carrier requires static refinement");
    }
  } else if (grid_setup.compare("audited_amr_runtime_local") == 0) {
    Q006RuntimeLocalRequireString(pin, block, "deck_role",
                                  "bounded_serial_host_audited_amr_mechanics_only");
    Q006RuntimeLocalRequireString(
        pin, block, "amr_policy",
        "deterministic_audited_randomized_runtime_local_10_refine_60_derefine_"
        "not_qualified");
    Q006RuntimeLocalRequireInteger(pin, block, "amr_seed", 60053);
    if (!pmy_mesh_->multilevel || !pmy_mesh_->adaptive) {
      Q006RuntimeLocalFatal("Q-006 runtime-local audited AMR carrier requires AMR");
    }
    q006_runtime_local_amr_seed = pin->GetInteger(block, "amr_seed");
    user_ref_func = Q006RuntimeLocalAuditedAMRRefinement;
  } else {
    Q006RuntimeLocalFatal("<q006_paper_multispecies_oscillation_runtime_local>/"
                          "grid_setup is unsupported");
  }

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
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
  const Real rho = pin->GetReal(block, "rho");
  const Real gas_uy = pin->GetReal(block, "gas_uy");
  const Real b_g = pin->GetReal(block, "b_g");

  par_for("pgen_q006_runtime_local_uniform_isothermal_state", DevExeSpace(),
          0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m, k, j, i) = 0.0;
    b0.x2f(m, k, j, i) = 0.0;
    b0.x3f(m, k, j, i) = b_g;
    if (i == ie) b0.x1f(m, k, j, i + 1) = 0.0;
    if (j == je) b0.x2f(m, k, j + 1, i) = 0.0;
    if (k == ke) b0.x3f(m, k + 1, j, i) = b_g;
    w0(m, IDN, k, j, i) = rho;
    w0(m, IVX, k, j, i) = 0.0;
    w0(m, IVY, k, j, i) = gas_uy;
    w0(m, IVZ, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) = 0.0;
    bcc0(m, IBY, k, j, i) = 0.0;
    bcc0(m, IBZ, k, j, i) = b_g;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
