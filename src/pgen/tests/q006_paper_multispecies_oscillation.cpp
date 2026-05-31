//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q006_paper_multispecies_oscillation.cpp
//! \brief Source-local Q-006 Sun and Bai Section 5.3 oscillation preparation.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

int q006_amr_seed = 60053;

[[noreturn]] void Q006Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q006RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (!std::isfinite(measured) ||
      std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q006Fatal(label + " does not match the Q-006 Section 5.3 preparation contract");
  }
}

void Q006RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q006Fatal("<" + block + ">/" + name +
              " does not match the Q-006 Section 5.3 preparation contract");
  }
}

void Q006RequireInteger(ParameterInput *pin, const std::string &block,
                        const std::string &name, const int expected) {
  if (pin->GetInteger(block, name) != expected) {
    Q006Fatal("<" + block + ">/" + name +
              " does not match the Q-006 Section 5.3 preparation contract");
  }
}

void Q006RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q006Fatal("<" + block + ">/" + name +
              " does not match the Q-006 Section 5.3 preparation contract");
  }
}

void Q006RequireReal(ParameterInput *pin, const std::string &block,
                     const std::string &name, const Real expected) {
  Q006RequireClose("<" + block + ">/" + name, pin->GetReal(block, name), expected);
}

std::uint64_t Q006SplitMix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30))*0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27))*0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

Real Q006AuditedAMRUniform01(const int seed, const int cycle, const int gid) {
  std::uint64_t key = static_cast<std::uint64_t>(seed);
  key ^= (static_cast<std::uint64_t>(cycle + 1))*0xbf58476d1ce4e5b9ULL;
  key ^= (static_cast<std::uint64_t>(gid + 1))*0xd2b74407b1ce6e93ULL;
  const std::uint64_t bits = Q006SplitMix64(key);
  return static_cast<Real>(bits >> 11)*
         static_cast<Real>(1.0/9007199254740992.0);
}

void Q006AuditedAMRRefinement(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  if (!pmesh->adaptive || pmesh->pmr == nullptr) {
    Q006Fatal("Q-006 audited AMR callback requires adaptive mesh refinement");
  }

  auto &refine_flag = pmesh->pmr->refine_flag;
  const int nmb = pmbp->nmb_thispack;
  const int mbs = pmesh->gids_eachrank[global_variable::my_rank];
  const int cycle = pmesh->ncycle;
  const Real refine_probability = static_cast<Real>(0.10);
  const Real derefine_probability = static_cast<Real>(0.60);
  for (int m = 0; m < nmb; ++m) {
    const Real draw = Q006AuditedAMRUniform01(q006_amr_seed, cycle, mbs + m);
    refine_flag.h_view(mbs + m) =
        (draw < refine_probability) ? 1 :
        ((draw < refine_probability + derefine_probability) ? -1 : 0);
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

}  // namespace

void ProblemGenerator::Q006PaperMultispeciesOscillation(ParameterInput *pin,
                                                         const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q006Fatal("q006_paper_multispecies_oscillation requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q006Fatal("q006_paper_multispecies_oscillation uses the ideal-MHD energy-feedback "
              "compatibility preparation because exact Section 5.3 isothermal "
              "coupled runtime is not supported");
  }
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic ||
      pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
      pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
      pin->DoesBlockExist("z4c") || pin->DoesBlockExist("turb_driving") ||
      pin->DoesBlockExist("initial_turb")) {
    Q006Fatal("q006_paper_multispecies_oscillation requires the exact Newtonian "
              "single-fluid MHD source task path");
  }
  if (!pmy_mesh_->three_d || pmy_mesh_->one_d || pmy_mesh_->two_d ||
      pmy_mesh_->mesh_indcs.nx1 != 16 || pmy_mesh_->mesh_indcs.nx2 != 8 ||
      pmy_mesh_->mesh_indcs.nx3 != 8 || pmy_mesh_->mb_indcs.nx1 != 4 ||
      pmy_mesh_->mb_indcs.nx2 != 4 || pmy_mesh_->mb_indcs.nx3 != 4) {
    Q006Fatal("q006_paper_multispecies_oscillation requires its exact 16x8x8 "
              "periodic 3D carrier with 4x4x4 MeshBlocks");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q006Fatal("q006_paper_multispecies_oscillation requires periodic boundaries");
  }
  Q006RequireClose("global x1min", pmy_mesh_->mesh_size.x1min, 0.0);
  Q006RequireClose("global x1max", pmy_mesh_->mesh_size.x1max, 16.0);
  Q006RequireClose("global x2min", pmy_mesh_->mesh_size.x2min, 0.0);
  Q006RequireClose("global x2max", pmy_mesh_->mesh_size.x2max, 8.0);
  Q006RequireClose("global x3min", pmy_mesh_->mesh_size.x3min, 0.0);
  Q006RequireClose("global x3max", pmy_mesh_->mesh_size.x3max, 8.0);

  Q006RequireString(pin, "time", "integrator", "rk2");
  Q006RequireReal(pin, "time", "cfl_number", 0.1);
  Q006RequireString(pin, "mhd", "eos", "ideal");
  Q006RequireString(pin, "mhd", "reconstruct", "plm");
  Q006RequireString(pin, "mhd", "rsolver", "llf");
  Q006RequireReal(pin, "mhd", "gamma", 1.66666666667);

  Q006RequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q006RequireReal(pin, "particles", "ppc", 128.0);
  Q006RequireString(pin, "particles", "pusher", "boris_tsc");
  Q006RequireInteger(pin, "particles", "nspecies", 2);
  Q006RequireString(pin, "particles", "cr_distribution", "center");
  Q006RequireBoolean(pin, "particles", "deposit_moments", true);
  Q006RequireInteger(pin, "particles", "deposit_order", 1);
  Q006RequireReal(pin, "particles", "deposit_qscale", 0.0234375);
  Q006RequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q006RequireReal(pin, "particles", "couple_j_to_efield_coeff", 1.0);
  Q006RequireString(pin, "particles", "couple_j_to_efield_representation",
                    "cell_centered");
  Q006RequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");
  Q006RequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q006RequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", true);
  Q006RequireString(pin, "particles", "couple_fluid_feedback_order",
                    "mhd_src_terms");
  Q006RequireReal(pin, "particles", "couple_moments_momentum_coeff", 1.0);
  Q006RequireReal(pin, "particles", "couple_moments_energy_coeff", 1.0);
  Q006RequireString(pin, "particles", "pic_physical_mode", "paper_mhd_pic");
  Q006RequireString(pin, "particles", "pic_background_mode", "coupled");
  Q006RequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q006RequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q006RequireReal(pin, "particles", "pic_cr_light_speed", 1000.0);
  Q006RequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  Q006RequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q006RequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q006RequireInteger(pin, "particles", "pic_max_cell_cross", 1);
  Q006RequireReal(pin, "particles", "pic_theta_max", 0.1);
  Q006RequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q006RequireString(pin, "particles", "pic_expanding_box_mode", "off");
  Q006RequireReal(pin, "species0", "mass", 1.0);
  Q006RequireReal(pin, "species0", "charge", -1.0);
  Q006RequireReal(pin, "species0", "vx0", 0.0);
  Q006RequireReal(pin, "species0", "vy0", 0.1);
  Q006RequireReal(pin, "species0", "vz0", 0.0);
  Q006RequireReal(pin, "species1", "mass", 1.0);
  Q006RequireReal(pin, "species1", "charge", 1.0);
  Q006RequireReal(pin, "species1", "vx0", 0.0);
  Q006RequireReal(pin, "species1", "vy0", 0.1);
  Q006RequireReal(pin, "species1", "vz0", 0.0);

  const std::string block = "q006_paper_multispecies_oscillation";
  Q006RequireString(pin, block, "campaign_id", "Q006-PAPER-MULTISPECIES-OSCILLATION");
  Q006RequireString(pin, block, "qualification_effect", "none");
  Q006RequireString(pin, block, "frontier_authorization", "not_bound");
  Q006RequireString(pin, block, "section53_runtime_evidence", "not_claimed");
  Q006RequireString(pin, block, "true_amr_policy_qualification", "not_claimed");
  Q006RequireString(pin, block, "mpi_qualification", "not_claimed");
  Q006RequireString(pin, block, "gpu_qualification", "not_claimed");
  Q006RequireString(pin, block, "external_review", "not_claimed");
  Q006RequireString(pin, block, "paper_eos_target", "isothermal_cs_1");
  Q006RequireString(pin, block, "runtime_eos_compatibility",
                    "ideal_energy_feedback_preparation_only");
  Q006RequireString(pin, block, "ppc_semantics",
                    "aggregate_128_round_robin_yields_64_per_species_per_cell");
  Q006RequireReal(pin, block, "rho", 1.0);
  Q006RequireReal(pin, block, "paper_cs", 1.0);
  Q006RequireReal(pin, block, "compatibility_pressure", 0.6);
  Q006RequireReal(pin, block, "b_g", 1.0);
  Q006RequireReal(pin, block, "omega", 1.0);
  Q006RequireReal(pin, block, "species_mass_density_ratio", 1.5);
  Q006RequireReal(pin, block, "gas_uy", -0.3);
  Q006RequireReal(pin, block, "species_vy", 0.1);
  Q006RequireReal(pin, block, "artificial_c", 1000.0);
  Q006RequireReal(pin, block, "dt_omega_target", 0.1);

  const std::string grid_setup = pin->GetString(block, "grid_setup");
  if (grid_setup.compare("uniform") == 0) {
    Q006RequireString(pin, block, "deck_role",
                      "source_local_uniform_deck_freeze_preparation_only_not_authorized");
    Q006RequireString(pin, block, "amr_policy", "not_applicable_uniform");
    if (pmy_mesh_->multilevel || pmy_mesh_->adaptive) {
      Q006Fatal("Q-006 uniform preparation rejects mesh refinement");
    }
  } else if (grid_setup.compare("smr") == 0) {
    Q006RequireString(pin, block, "deck_role",
                      "source_local_smr_deck_freeze_preparation_only_not_authorized");
    Q006RequireString(pin, block, "amr_policy", "paper_static_one_eighth_region");
    if (!pmy_mesh_->multilevel || pmy_mesh_->adaptive) {
      Q006Fatal("Q-006 SMR preparation requires static mesh refinement");
    }
  } else if (grid_setup.compare("audited_amr_preparation") == 0) {
    Q006RequireString(pin, block, "deck_role",
                      "source_local_audited_amr_deck_freeze_preparation_only_"
                      "not_authorized");
    Q006RequireString(
        pin, block, "amr_policy",
        "deterministic_audited_randomized_preparation_10_refine_60_derefine_"
        "not_qualified");
    Q006RequireInteger(pin, block, "amr_seed", 60053);
    if (!pmy_mesh_->multilevel || !pmy_mesh_->adaptive) {
      Q006Fatal("Q-006 audited AMR preparation requires adaptive mesh refinement");
    }
    q006_amr_seed = pin->GetInteger(block, "amr_seed");
    user_ref_func = Q006AuditedAMRRefinement;
  } else {
    Q006Fatal("<q006_paper_multispecies_oscillation>/grid_setup is unsupported");
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
  const Real pressure = pin->GetReal(block, "compatibility_pressure");
  const Real gas_uy = pin->GetReal(block, "gas_uy");
  const Real b_g = pin->GetReal(block, "b_g");

  par_for("pgen_q006_uniform_paper_state", DevExeSpace(),
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
    w0(m, IEN, k, j, i) = pressure;
    bcc0(m, IBX, k, j, i) = 0.0;
    bcc0(m, IBY, k, j, i) = 0.0;
    bcc0(m, IBZ, k, j, i) = b_g;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}
