//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q032_reduced_static_neutral_local.cpp
//! \brief Bounded local mechanics carrier for the reduced static-neutral damping map.

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
#include "pgen/pgen.hpp"

namespace {

[[noreturn]] void Q032Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Q032RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q032Fatal(label + " does not match the Q-032 bounded local mechanics contract");
  }
}

void Q032RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q032Fatal("<" + block + ">/" + name +
              " does not match the Q-032 bounded local mechanics contract");
  }
}

void Q032RequireInteger(ParameterInput *pin, const std::string &block,
                        const std::string &name, const int expected) {
  if (pin->GetInteger(block, name) != expected) {
    Q032Fatal("<" + block + ">/" + name +
              " does not match the Q-032 bounded local mechanics contract");
  }
}

void Q032RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q032Fatal("<" + block + ">/" + name +
              " does not match the Q-032 bounded local mechanics contract");
  }
}

void Q032RequireReal(ParameterInput *pin, const std::string &block,
                     const std::string &name, const Real expected) {
  Q032RequireClose("<" + block + ">/" + name, pin->GetReal(block, name), expected);
}

}  // namespace

void ProblemGenerator::Q032ReducedStaticNeutralLocal(ParameterInput *pin,
                                                      const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q032Fatal("q032_reduced_static_neutral_local requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q032Fatal("q032_reduced_static_neutral_local requires ideal-MHD EOS");
  }
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic ||
      pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
      pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
      pin->DoesBlockExist("z4c")) {
    Q032Fatal("q032_reduced_static_neutral_local requires the exact Newtonian "
              "single-fluid MHD source task path");
  }
  if (pmy_mesh_->multilevel) {
    Q032Fatal("q032_reduced_static_neutral_local rejects AMR/SMR");
  }
  if (!pmy_mesh_->two_d || pmy_mesh_->one_d || pmy_mesh_->three_d ||
      pmy_mesh_->mesh_indcs.nx1 != 32 || pmy_mesh_->mesh_indcs.nx2 != 4 ||
      pmy_mesh_->mesh_indcs.nx3 != 1 || pmy_mesh_->mb_indcs.nx1 != 32 ||
      pmy_mesh_->mb_indcs.nx2 != 4 || pmy_mesh_->mb_indcs.nx3 != 1) {
    Q032Fatal("q032_reduced_static_neutral_local requires the exact thin "
              "32x4x1 2D3V x1-parallel carrier");
  }
  if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic ||
      pmy_mesh_->mesh_bcs[BoundaryFace::outer_x1] != BoundaryFlag::periodic ||
      pmy_mesh_->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
      pmy_mesh_->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic) {
    Q032Fatal("q032_reduced_static_neutral_local requires periodic active "
              "carrier boundaries");
  }

  Q032RequireString(pin, "time", "integrator", "rk1");
  Q032RequireReal(pin, "time", "cfl_number", 0.1);
  Q032RequireString(pin, "mhd", "eos", "ideal");
  Q032RequireString(pin, "mhd", "reconstruct", "plm");
  Q032RequireString(pin, "mhd", "rsolver", "llf");
  Q032RequireReal(pin, "mhd", "gamma", 1.66666666667);

  Q032RequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q032RequireReal(pin, "particles", "ppc", 0.0);
  Q032RequireString(pin, "particles", "pusher", "boris_tsc");
  Q032RequireInteger(pin, "particles", "nspecies", 1);
  Q032RequireString(pin, "particles", "cr_distribution", "center");
  Q032RequireBoolean(pin, "particles", "deposit_moments", false);
  Q032RequireBoolean(pin, "particles", "couple_moments_to_mhd", false);
  Q032RequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", false);
  Q032RequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", false);
  Q032RequireReal(pin, "particles", "cr_vx0", 0.0);
  Q032RequireReal(pin, "particles", "cr_vy0", 0.0);
  Q032RequireReal(pin, "particles", "cr_vz0", 0.0);
  Q032RequireString(pin, "particles", "pic_physical_mode", "extended_mhd_pic");
  Q032RequireString(pin, "particles", "pic_background_mode", "coupled");
  Q032RequireString(pin, "particles", "pic_feedback_mode", "test_particle");
  Q032RequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q032RequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q032RequireReal(pin, "particles", "pic_cr_light_speed", 3.0);
  Q032RequireString(pin, "particles", "pic_cr_initial_state", "momentum");
  Q032RequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q032RequireString(pin, "particles", "pic_wave_damping_mode",
                    "ion_neutral_friction");
  Q032RequireReal(pin, "particles", "pic_ion_neutral_collision_rate", 0.7);
  Q032RequireInteger(pin, "particles", "pic_max_cell_cross", 2);
  Q032RequireReal(pin, "particles", "pic_theta_max", 0.3);
  Q032RequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q032RequireString(pin, "particles", "pic_expanding_box_mode", "off");
  Q032RequireReal(pin, "species0", "mass", 1.0);
  Q032RequireReal(pin, "species0", "charge", 0.0);

  Q032RequireInteger(pin, "problem", "wave_flag", 1);
  Q032RequireReal(pin, "problem", "amp", 0.01);
  Q032RequireReal(pin, "problem", "vflow", 0.2);
  Q032RequireBoolean(pin, "problem", "along_x1", true);
  Q032RequireBoolean(pin, "problem", "along_x2", false);
  Q032RequireBoolean(pin, "problem", "along_x3", false);

  const std::string block = "q032_reduced_static_neutral_local";
  Q032RequireString(pin, block, "campaign_id", "Q032-REDUCED-STATIC-NEUTRAL-LOCAL");
  Q032RequireString(pin, block, "deck_role",
                    "bounded_source_local_mechanics_only_not_plotnikov_evidence");
  Q032RequireString(pin, block, "qualification_effect", "none");
  Q032RequireString(pin, block, "frontier_authorization", "not_bound");
  Q032RequireString(pin, block, "plotnikov_qualification", "not_claimed");
  Q032RequireString(pin, block, "physical_scope",
                    "reduced_static_neutral_high_frequency_transverse_friction_only");
  Q032RequireString(pin, block, "carrier_semantics",
                    "transverse_invariant_thin_2d3v_x1_parallel_alfven");
  Q032RequireString(pin, block, "nu_in_parameter",
                    "particles/pic_ion_neutral_collision_rate");
  Q032RequireReal(pin, block, "nu_in", 0.7);
  Q032RequireString(pin, block, "attenuation_factor", "exp(-nu_in*dt)");

  // Keep the MHD initialization and restart-derived boundary-state behavior in
  // the reviewed built-in linear-wave implementation.
  LinearWave(pin, restart);
}
