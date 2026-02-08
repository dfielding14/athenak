//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include <vector>
#include <Kokkos_Random.hpp>
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"
#include "units/units.hpp"
#include "utils/sn_scheduler.hpp"

namespace particles {
namespace {

const char *BoundaryFaceName(BoundaryFace face) {
  switch (face) {
    case BoundaryFace::inner_x1:
      return "mesh/ix1_bc";
    case BoundaryFace::outer_x1:
      return "mesh/ox1_bc";
    case BoundaryFace::inner_x2:
      return "mesh/ix2_bc";
    case BoundaryFace::outer_x2:
      return "mesh/ox2_bc";
    case BoundaryFace::inner_x3:
      return "mesh/ix3_bc";
    case BoundaryFace::outer_x3:
      return "mesh/ox3_bc";
    default:
      return "mesh/?_bc";
  }
}

const char *BoundaryFlagName(BoundaryFlag flag) {
  switch (flag) {
    case BoundaryFlag::reflect:
      return "reflect";
    case BoundaryFlag::inflow:
      return "inflow";
    case BoundaryFlag::outflow:
      return "outflow";
    case BoundaryFlag::diode:
      return "diode";
    case BoundaryFlag::user:
      return "user";
    case BoundaryFlag::periodic:
      return "periodic";
    case BoundaryFlag::shear_periodic:
      return "shear_periodic";
    case BoundaryFlag::vacuum:
      return "vacuum";
    case BoundaryFlag::block:
      return "block";
    default:
      return "undef";
  }
}

bool MomentBoundaryFlagSupported(BoundaryFlag flag) {
  return ((flag == BoundaryFlag::periodic) ||
          (flag == BoundaryFlag::outflow) ||
          (flag == BoundaryFlag::reflect));
}

void ValidateMomentBoundaryPolicy(Mesh *pmesh) {
  const std::array<BoundaryFace, 6> faces = {
      BoundaryFace::inner_x1, BoundaryFace::outer_x1,
      BoundaryFace::inner_x2, BoundaryFace::outer_x2,
      BoundaryFace::inner_x3, BoundaryFace::outer_x3};

  for (const auto face : faces) {
    BoundaryFlag flag = pmesh->mesh_bcs[face];
    if (!MomentBoundaryFlagSupported(flag)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/deposit_moments=true does not support "
                << BoundaryFaceName(face) << "=" << BoundaryFlagName(flag)
                << " in PR3c (supported: periodic/outflow/reflect)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // initialize vector for star particles
  std::vector<std::array<Real, 9>> particle_list;

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
      // read number of particles per cell, and calculate number of particles this pack
      Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
      // compute number of particles as real number, since ppc can be < 1
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
      Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
      nprtcl_thispack = static_cast<int>(r_npart); // then cast to integer
    } else if (ptype.compare("star") == 0) {
      particle_type = ParticleType::star;
      nprtcl_thispack = 0; // initialize to zero

      // Load particles from file
      std::string particle_file = pin->GetString("particles","star_particle_file");
      std::ifstream infile(particle_file);
      if (!infile) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Unable to open particle file: " << particle_file
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // Read particle positions from file
      int nmb = pmy_pack->nmb_thispack;
      auto &size = pmy_pack->pmb->mb_size;

      // Skip header comments
      while (infile.peek() == '#') {
        infile.ignore(1000, '\n');
      }

      std::array<Real, 8> p;
      while (infile >> p[0] >> p[1] >> p[2] >> p[3] >>
             p[4] >> p[5] >> p[6] >> p[7]) {
        for (int m=0; m<nmb; ++m) {
          // Loop over particles to see which are in this meshblock
          if (p[0] > size.h_view(m).x1min && p[0] <= size.h_view(m).x1max &&
              p[1] > size.h_view(m).x2min && p[1] <= size.h_view(m).x2max &&
              p[2] > size.h_view(m).x3min && p[2] <= size.h_view(m).x3max) {
            // Add particle to the list if it is within the mesh block bounds
            std::array<Real, 9> new_particle = {p[0], p[1], p[2],
                                                p[3], p[4], p[5],
                                                p[6], p[7], static_cast<Real>(m)};
            particle_list.push_back(new_particle);
          }
        }
      }
      infile.close();
      nprtcl_thispack = particle_list.size();

      // Print number of particles loaded
      std::cout << "Loaded " << nprtcl_thispack
                << " star particles from file " << particle_file << std::endl;
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("rk4_gravity") == 0) {
      pusher = ParticlesPusher::rk4_gravity;
      // load gravity constants
      r_scale   = pin->GetReal("potential", "r_scale");
      rho_scale = pin->GetReal("potential", "rho_scale");
      m_gal     = pin->GetReal("potential", "mass_gal");
      a_gal     = pin->GetReal("potential", "scale_gal");
      z_gal     = pin->GetReal("potential", "z_gal");
      r_200     = pin->GetReal("potential", "r_200");
      rho_mean  = pin->GetReal("potential", "rho_mean");
      par_grav_dx = pin->GetOrAddReal("particles", "grav_dx", 1e-6);
    } else if (ppush.compare("boris_lin") == 0) {
      pusher = ParticlesPusher::boris_lin;
    } else if (ppush.compare("boris_tsc") == 0) {
      pusher = ParticlesPusher::boris_tsc;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particles only work in 2D/3D, but 1D problem initialized" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // stars must be 3D
  if (particle_type == ParticleType::star && !pmy_pack->pmesh->three_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Star particles only work in 3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        // Determine number of real data slots based on dimensionality and options
        track_displacement = pin->GetOrAddBoolean("particles","track_displacement",false);
        if (track_displacement) {
          nrdata = IPDB + 1;  // includes displacement indices
        } else if (pmy_pack->pmesh->three_d) {
          nrdata = IPBZ + 1;  // up to Bz
        } else {
          nrdata = IPBY + 1;  // up to By in 2D
        }

        nidata = 3;  // PGID, PTAG, PSP (species)
        break;
      }
    case ParticleType::star:
      {
        nrdata = 9;
        nidata = 3;
        break;
      }
    default:
      break;
  }

  // PR1 deposition controls
  deposit_moments = pin->GetOrAddBoolean("particles", "deposit_moments", false);
  deposit_order = pin->GetOrAddInteger("particles", "deposit_order", 1);
  deposit_qscale = pin->GetOrAddReal("particles", "deposit_qscale", 1.0);
  couple_moments_to_mhd = pin->GetOrAddBoolean("particles",
                                               "couple_moments_to_mhd", false);
  couple_j_to_efield_coeff = pin->GetOrAddReal("particles",
                                                "couple_j_to_efield_coeff", 1.0);
  std::string j_repr = pin->GetOrAddString("particles",
                                           "couple_j_to_efield_representation",
                                           "cell_centered");
  if (j_repr.compare("cell_centered") == 0) {
    couple_j_to_efield_representation = CoupledCurrentRepresentation::cell_centered;
  } else if (j_repr.compare("edge_staggered") == 0) {
    couple_j_to_efield_representation = CoupledCurrentRepresentation::edge_staggered;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/couple_j_to_efield_representation: "
              << j_repr << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string j_deposit_default =
      (couple_j_to_efield_representation == CoupledCurrentRepresentation::edge_staggered)
          ? "direct_staggered"
          : "cc_convert";
  std::string j_deposit_mode = pin->GetOrAddString(
      "particles", "couple_j_deposition_mode", j_deposit_default);
  if (j_deposit_mode.compare("cc_convert") == 0) {
    couple_j_deposition_mode = CoupledCurrentDepositionMode::cc_convert;
  } else if (j_deposit_mode.compare("direct_staggered") == 0) {
    couple_j_deposition_mode = CoupledCurrentDepositionMode::direct_staggered;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/couple_j_deposition_mode: "
              << j_deposit_mode << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string feedback_order = pin->GetOrAddString("particles",
                                                    "couple_fluid_feedback_order",
                                                    "mhd_src_terms");
  if (feedback_order.compare("mhd_src_terms") == 0) {
    couple_fluid_feedback_order = CoupledFluidFeedbackOrder::mhd_src_terms;
  } else if (feedback_order.compare("efield_src") == 0) {
    couple_fluid_feedback_order = CoupledFluidFeedbackOrder::efield_src;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/couple_fluid_feedback_order: "
              << feedback_order << std::endl;
    std::exit(EXIT_FAILURE);
  }
  couple_moments_momentum_to_mhd = pin->GetOrAddBoolean(
      "particles", "couple_moments_momentum_to_mhd", false);
  couple_moments_energy_to_mhd = pin->GetOrAddBoolean(
      "particles", "couple_moments_energy_to_mhd", false);
  couple_moments_momentum_coeff = pin->GetOrAddReal(
      "particles", "couple_moments_momentum_coeff", 1.0);
  couple_moments_energy_coeff = pin->GetOrAddReal(
      "particles", "couple_moments_energy_coeff", 1.0);
  cr_vx0 = pin->GetOrAddReal("particles", "cr_vx0", 0.0);
  cr_vy0 = pin->GetOrAddReal("particles", "cr_vy0", 0.0);
  cr_vz0 = pin->GetOrAddReal("particles", "cr_vz0", 0.0);

  // Staged PR5+ PIC runtime controls (parse + validation only at this step)
  std::string pic_background_mode_str = pin->GetOrAddString(
      "particles", "pic_background_mode", "coupled");
  if (pic_background_mode_str.compare("coupled") == 0) {
    pic_background_mode = PICBackgroundMode::coupled;
  } else if (pic_background_mode_str.compare("passive_mhd") == 0) {
    pic_background_mode = PICBackgroundMode::passive_mhd;
  } else if (pic_background_mode_str.compare("no_mhd") == 0) {
    pic_background_mode = PICBackgroundMode::no_mhd;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_background_mode: "
              << pic_background_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const char *pic_feedback_mode_default =
      ((pic_background_mode == PICBackgroundMode::passive_mhd) ||
       (pic_background_mode == PICBackgroundMode::no_mhd)) ?
      "test_particle" : "coupled";
  std::string pic_feedback_mode_str = pin->GetOrAddString(
      "particles", "pic_feedback_mode", pic_feedback_mode_default);
  if (pic_feedback_mode_str.compare("coupled") == 0) {
    pic_feedback_mode = PICFeedbackMode::coupled;
  } else if (pic_feedback_mode_str.compare("test_particle") == 0) {
    pic_feedback_mode = PICFeedbackMode::test_particle;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_feedback_mode: "
              << pic_feedback_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string pic_interp_scheme_str = pin->GetOrAddString(
      "particles", "pic_interp_scheme", "tsc");
  if (pic_interp_scheme_str.compare("tsc") == 0) {
    pic_interp_scheme = PICInterpolationScheme::tsc;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_interp_scheme: "
              << pic_interp_scheme_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_cr_light_speed = pin->GetOrAddReal("particles", "pic_cr_light_speed", 1.0);
  if (pic_cr_light_speed <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_cr_light_speed must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_max_cell_cross = pin->GetOrAddInteger("particles", "pic_max_cell_cross", 2);
  if (pic_max_cell_cross <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_max_cell_cross must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_theta_max = pin->GetOrAddReal("particles", "pic_theta_max", 0.3);
  if (pic_theta_max <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_theta_max must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string pic_deltaf_mode_str = pin->GetOrAddString(
      "particles", "pic_deltaf_mode", "off");
  if (pic_deltaf_mode_str.compare("off") == 0) {
    pic_deltaf_mode = PICDeltaFMode::off;
  } else if (pic_deltaf_mode_str.compare("on") == 0) {
    pic_deltaf_mode = PICDeltaFMode::on;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_deltaf_mode: "
              << pic_deltaf_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_deltaf_f0 = pin->GetOrAddString("particles", "pic_deltaf_f0", "");
  if ((pic_deltaf_mode == PICDeltaFMode::on) && pic_deltaf_f0.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_deltaf_mode=on requires "
              << "<particles>/pic_deltaf_f0 to define background f0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_sort_interval = pin->GetOrAddInteger("particles", "pic_sort_interval", 0);
  if (pic_sort_interval < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_sort_interval must be >= 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string pic_intermediate_arrays_mode_str = pin->GetOrAddString(
      "particles", "pic_intermediate_arrays", "auto");
  if (pic_intermediate_arrays_mode_str.compare("auto") == 0) {
    pic_intermediate_arrays_mode = PICIntermediateArraysMode::auto_mode;
  } else if (pic_intermediate_arrays_mode_str.compare("off") == 0) {
    pic_intermediate_arrays_mode = PICIntermediateArraysMode::off;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_intermediate_arrays: "
              << pic_intermediate_arrays_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string pic_expanding_box_mode_str = pin->GetOrAddString(
      "particles", "pic_expanding_box_mode", "off");
  if (pic_expanding_box_mode_str.compare("off") == 0) {
    pic_expanding_box_mode = PICExpandingBoxMode::off;
  } else if (pic_expanding_box_mode_str.compare("on") == 0) {
    pic_expanding_box_mode = PICExpandingBoxMode::on;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_expanding_box_mode: "
              << pic_expanding_box_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_expansion_rate_x1 = pin->GetOrAddReal("particles", "pic_expansion_rate_x1", 0.0);
  pic_expansion_rate_x2 = pin->GetOrAddReal("particles", "pic_expansion_rate_x2", 0.0);
  pic_expansion_rate_x3 = pin->GetOrAddReal("particles", "pic_expansion_rate_x3", 0.0);
  pic_no_mhd_bx = pin->GetOrAddReal("particles", "pic_no_mhd_bx", 0.0);
  pic_no_mhd_by = pin->GetOrAddReal("particles", "pic_no_mhd_by", 0.0);
  pic_no_mhd_bz = pin->GetOrAddReal("particles", "pic_no_mhd_bz", 0.0);
  if ((pic_expanding_box_mode == PICExpandingBoxMode::off) &&
      ((pic_expansion_rate_x1 != 0.0) ||
       (pic_expansion_rate_x2 != 0.0) ||
       (pic_expansion_rate_x3 != 0.0))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_expansion_rate_x1/x2/x3 require "
              << "<particles>/pic_expanding_box_mode=on" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if ((pic_feedback_mode == PICFeedbackMode::test_particle) &&
      (couple_moments_to_mhd || couple_moments_momentum_to_mhd ||
       couple_moments_energy_to_mhd)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_feedback_mode=test_particle does not support "
              << "particle-to-MHD coupling toggles "
              << "(couple_moments_to_mhd and momentum/energy feedback)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (pic_background_mode == PICBackgroundMode::passive_mhd) {
    if (!(pin->DoesBlockExist("mhd"))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_background_mode=passive_mhd requires an active "
                << "<mhd> block" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pic_feedback_mode != PICFeedbackMode::test_particle) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_background_mode=passive_mhd requires "
                << "<particles>/pic_feedback_mode=test_particle unless coupled "
                << "feedback is explicitly implemented for this mode"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pic_background_mode == PICBackgroundMode::no_mhd) {
    if (pic_feedback_mode != PICFeedbackMode::test_particle) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_background_mode=no_mhd requires "
                << "<particles>/pic_feedback_mode=test_particle" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  const bool use_fluid_feedback = (couple_moments_momentum_to_mhd ||
                                   couple_moments_energy_to_mhd);

  // PR1 runtime scope guard for moment deposition
  if (deposit_moments) {
    if (particle_type != ParticleType::cosmic_ray) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<particles>/deposit_moments=true requires "
                << "<particles>/particle_type=cosmic_ray" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const bool direct_edge_mode =
        (couple_moments_to_mhd &&
         couple_j_deposition_mode ==
         CoupledCurrentDepositionMode::direct_staggered);
    const bool valid_direct_order = (deposit_order == 1 || deposit_order == 2);
    if ((!direct_edge_mode && deposit_order != 1) ||
        (direct_edge_mode && !valid_direct_order)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<particles>/deposit_order=" << deposit_order
                << " is not supported (only deposit_order=1, or "
                << "deposit_order={1,2} with coupled "
                << "couple_j_deposition_mode=direct_staggered)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ValidateMomentBoundaryPolicy(pmy_pack->pmesh);
    if (pin->DoesBlockExist("shearing_box")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<particles>/deposit_moments=true does not support "
                << "<shearing_box> in PR1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // PR2 runtime scope guard for current coupling
  if (use_fluid_feedback && !couple_moments_to_mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Fluid feedback coupling requires "
              << "<particles>/couple_moments_to_mhd=true" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (couple_moments_to_mhd) {
    if (!deposit_moments) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true requires "
                << "<particles>/deposit_moments=true" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(pin->DoesBlockExist("mhd"))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true requires an active "
                << "<mhd> block in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesBlockExist("radiation")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true does not support "
                << "radiation+MHD compositions in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("hydro")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true does not support "
                << "ion-neutral/two-fluid compositions in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->DoesBlockExist("adm") || pin->DoesBlockExist("z4c")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true does not support "
                << "numerical relativity task paths in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (couple_moments_energy_to_mhd) {
      std::string mhd_eos = pin->GetString("mhd", "eos");
      if (mhd_eos == "isothermal") {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<particles>/couple_moments_energy_to_mhd=true requires "
                  << "<mhd>/eos=ideal" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if ((couple_j_deposition_mode == CoupledCurrentDepositionMode::direct_staggered) &&
        (couple_j_to_efield_representation !=
         CoupledCurrentRepresentation::edge_staggered)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_j_deposition_mode=direct_staggered requires "
                << "<particles>/couple_j_to_efield_representation=edge_staggered"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (couple_j_to_efield_representation == CoupledCurrentRepresentation::edge_staggered
        && (pmy_pack->pcoord->is_special_relativistic ||
            pmy_pack->pcoord->is_general_relativistic ||
            pmy_pack->pcoord->is_dynamical_relativistic)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_j_to_efield_representation=edge_staggered "
                << "requires non-relativistic Cartesian MHD in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pic_background_mode == PICBackgroundMode::no_mhd && pin->DoesBlockExist("mhd")) {
    std::cout << "### WARNING in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_background_mode=no_mhd is active while an <mhd> "
              << "block is present; no_mhd particle field carrier will be used"
              << std::endl;
  }
  if (use_fluid_feedback &&
      (pmy_pack->pcoord->is_special_relativistic ||
       pmy_pack->pcoord->is_general_relativistic)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Fluid momentum/energy feedback is limited to non-relativistic "
              << "MHD in PR2" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);

  if (deposit_moments) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(moments, nmb, NMOM, ncells3, ncells2, ncells1);
    if (ppack->pmesh->multilevel) {
      int n_ccells1 = indcs.cnx1 + 2*(indcs.ng);
      int n_ccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
      int n_ccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(coarse_moments, nmb, NMOM, n_ccells3, n_ccells2, n_ccells1);
      Kokkos::deep_copy(coarse_moments, static_cast<Real>(0.0));
    }

    Kokkos::realloc(x1_old, nprtcl_thispack);
    Kokkos::realloc(x2_old, nprtcl_thispack);
    Kokkos::realloc(x3_old, nprtcl_thispack);

    pbval_mom = new MeshBoundaryValuesCC(ppack, pin, false, CCCommMode::synchronize);
    pbval_mom->InitializeBuffers(NMOM);

    if (couple_moments_to_mhd &&
        (couple_j_to_efield_representation ==
         CoupledCurrentRepresentation::edge_staggered)) {
      Kokkos::realloc(j_edge_x1e, nmb, ncells3+1, ncells2+1, ncells1);
      Kokkos::realloc(j_edge_x2e, nmb, ncells3+1, ncells2, ncells1+1);
      Kokkos::realloc(j_edge_x3e, nmb, ncells3, ncells2+1, ncells1+1);
      Kokkos::deep_copy(j_edge_x1e, static_cast<Real>(0.0));
      Kokkos::deep_copy(j_edge_x2e, static_cast<Real>(0.0));
      Kokkos::deep_copy(j_edge_x3e, static_cast<Real>(0.0));
      if (couple_j_deposition_mode == CoupledCurrentDepositionMode::direct_staggered) {
        pbval_jedge = new MeshBoundaryValuesFC(ppack, pin);
        pbval_jedge->InitializeBuffers(3);
      }
    }
  }

  if (pic_background_mode == PICBackgroundMode::no_mhd) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb = ppack->nmb_thispack;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(pic_no_mhd_bcc0, nmb, NMAG, ncells3, ncells2, ncells1);
    auto bcc = pic_no_mhd_bcc0;
    const Real bx = pic_no_mhd_bx;
    const Real by = pic_no_mhd_by;
    const Real bz = pic_no_mhd_bz;
    par_for("pic_no_mhd_bcc_init", DevExeSpace(), 0, nmb - 1,
            0, ncells3 - 1, 0, ncells2 - 1, 0, ncells1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      bcc(m, IBX, k, j, i) = bx;
      bcc(m, IBY, k, j, i) = by;
      bcc(m, IBZ, k, j, i) = bz;
    });
  }

  // Initialize particles based on type
  if (particle_type == ParticleType::cosmic_ray) {
    InitializeCosmicRays(pin);
  } else if (particle_type == ParticleType::star) {
    InitializeStars(particle_list);
  }
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
  delete pbval_part;
  if (pbval_mom != nullptr) {
    delete pbval_mom;
  }
  if (pbval_jedge != nullptr) {
    delete pbval_jedge;
  }
}

//----------------------------------------------------------------------------------------
// InitializeCosmicRays()
// Initializes cosmic ray particles with species support

void Particles::InitializeCosmicRays(ParameterInput *pin) {
  // Read number of species
  nspecies = pin->GetOrAddInteger("particles","nspecies",1);

  // Allocate species arrays
  Kokkos::realloc(species_mass, nspecies);
  Kokkos::realloc(species_charge, nspecies);

  // Read species properties
  auto h_mass = Kokkos::create_mirror_view(species_mass);
  auto h_charge = Kokkos::create_mirror_view(species_charge);

  for (int s=0; s<nspecies; ++s) {
    std::string block = "species" + std::to_string(s);
    h_mass(s) = pin->GetOrAddReal(block,"mass",1.0);
    h_charge(s) = pin->GetOrAddReal(block,"charge",1.0);
  }

  Kokkos::deep_copy(species_mass, h_mass);
  Kokkos::deep_copy(species_charge, h_charge);

  // Initialize particle positions and velocities
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  int nmb = pmy_pack->nmb_thispack;
  const int gids_local = pmy_pack->gids;

  // Determine distribution type for initial positions
  std::string dist = pin->GetOrAddString("particles", "cr_distribution", "center");
  const bool random_dist = (dist.compare("random") == 0);
  const bool track_disp_local = track_displacement;  // avoid capturing 'this'
  const int nx3_local = indcs.nx3;                   // avoid capturing host ref

  // Avoid division by zero if there are fewer particles than meshblocks
  int particles_per_mb = std::max(1, (nprtcl_thispack + nmb - 1) / nmb);

  // Random number pool for optional random placement
  Kokkos::Random_XorShift64_Pool<DevExeSpace> rand_pool64(pmy_pack->gids);

  // Create local copies to avoid GPU lambda capture warnings
  int nspecies_local = nspecies;
  auto species_charge_local = species_charge;
  auto species_mass_local = species_mass;
  const Real cr_vx0_local = cr_vx0;
  const Real cr_vy0_local = cr_vy0;
  const Real cr_vz0_local = cr_vz0;
  // Make sure geometry is available on device
  size.template sync<DevExeSpace>();
  auto size_view = size;
  // Capture handle by value for use as size_view.d_view(...) on device.

  // Simple uniform distribution for now
  par_for("init_cr", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    // Determine which meshblock this particle belongs to
    int m = p / particles_per_mb;
    if (m >= nmb) m = nmb - 1;

    // Set GID and species
    pi(PGID,p) = gids_local + m;
    pi(PTAG,p) = p;
    pi(PSP,p) = p % nspecies_local;  // Round-robin species assignment

    // Choose position within the mesh block
    Real rx = 0.5, ry = 0.5, rz = 0.5;
    if (random_dist) {
      auto rand_gen = rand_pool64.get_state();
      rx = rand_gen.drand();
      ry = rand_gen.drand();
      if (nx3_local > 1) {
        rz = rand_gen.drand();
      }
      rand_pool64.free_state(rand_gen);
    }

    const Real x1min = size_view.d_view(m).x1min;
    const Real x1max = size_view.d_view(m).x1max;
    const Real x2min = size_view.d_view(m).x2min;
    const Real x2max = size_view.d_view(m).x2max;
    const Real x3min = size_view.d_view(m).x3min;
    const Real x3max = size_view.d_view(m).x3max;

    pr(IPX,p) = x1min + rx*(x1max - x1min);
    pr(IPY,p) = x2min + ry*(x2max - x2min);
    pr(IPZ,p) = (nx3_local > 1) ? (x3min + rz*(x3max - x3min)) : 0.0;

    // Initialize velocity with deterministic CR controls
    pr(IPVX,p) = cr_vx0_local;
    pr(IPVY,p) = cr_vy0_local;
    pr(IPVZ,p) = cr_vz0_local;

    // Set mass/charge ratio
    int species = pi(PSP,p);
    pr(IPM,p) = species_charge_local(species) / species_mass_local(species);

    // Initialize B-field components to zero
    pr(IPBX,p) = 0.0;
    pr(IPBY,p) = 0.0;
    if (nx3_local > 1) {
      pr(IPBZ,p) = 0.0;
    }

    // Initialize displacement tracking if enabled
    if (track_disp_local) {
      pr(IPDX,p) = 0.0;
      pr(IPDY,p) = 0.0;
      pr(IPDZ,p) = 0.0;
      pr(IPDB,p) = 0.0;
    }
  });

  // Set timestep
  auto &dx = size.h_view(0);
  dtnew = std::min(dx.dx1, dx.dx2);
  if (indcs.nx3 > 1) {
    dtnew = std::min(dtnew, dx.dx3);
  }
}

//----------------------------------------------------------------------------------------
// InitializeStars()
// Initializes star particles from pre-loaded particle list

void Particles::InitializeStars(std::vector<std::array<Real, 9>> &particle_list) {
  auto &size = pmy_pack->pmb->mb_size;
  const int &gids = pmy_pack->gids;

  // Copy to device-accessible arrays
  // First create and populate a host view
  HostArray2D<Real> host_pos("host_positions", 9, nprtcl_thispack);
  for (size_t i = 0; i < nprtcl_thispack; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      host_pos(j, i) = particle_list[i][j];
    }
  }

  // Then create the device view and copy data
  auto pos_data = Kokkos::create_mirror_view_and_copy(DevExeSpace(), host_pos);

  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  Real unit_time = pmy_pack->punit->time_cgs();

  // Initialize particles
  par_for("star_par", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    int m = static_cast<int>(pos_data(8, p));
    pi(PGID,p) = gids + m;
    pi(NSN,p) = 0;  // track number of SNe for star particle
    pr(IPX,p)  = pos_data(0, p);
    pr(IPY,p)  = pos_data(1, p);
    pr(IPZ,p)  = pos_data(2, p);
    pr(IPVX,p) = pos_data(3, p);
    pr(IPVY,p) = pos_data(4, p);
    pr(IPVZ,p) = pos_data(5, p);
    pr(IPT_CREATE, p) = pos_data(6, p);  // creation time of star particle
    pr(IPMASS, p)     = pos_data(7, p);  // mass of star particle
    pr(IPT_NEXT_SN,p) = GetNthSNTime(pr(IPMASS,p), pr(IPT_CREATE,p), unit_time, 0);
  });

  dtnew = std::min(size.h_view(0).dx1, size.h_view(0).dx2);
  dtnew = std::min(dtnew, size.h_view(0).dx3);
}

//----------------------------------------------------------------------------------------
// CreatePaticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1];
    }

    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = tagstart + p;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = myrank + nranks*p;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assinment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace particles
