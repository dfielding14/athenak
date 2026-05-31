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
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <vector>
#include <cstdlib>
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "mhd/mhd.hpp"
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
          (flag == BoundaryFlag::inflow) ||
          (flag == BoundaryFlag::outflow) ||
          (flag == BoundaryFlag::reflect));
}

bool DirectEdgeCurrentBoundaryFlagSupported(BoundaryFlag flag) {
  return ((flag == BoundaryFlag::periodic) ||
          (flag == BoundaryFlag::outflow) ||
          (flag == BoundaryFlag::reflect));
}

int ParticlesInGlobalMeshBlock(const Real particles_per_mb, const int gid) {
  const Real left = std::floor(particles_per_mb*static_cast<Real>(gid));
  const Real right = std::floor(particles_per_mb*static_cast<Real>(gid + 1));
  const Real nblock = std::max(static_cast<Real>(0.0), right - left);
  if (nblock > static_cast<Real>(std::numeric_limits<int>::max())) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/ppc creates too many particles for one MeshBlock"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return static_cast<int>(nblock);
}

int CountCosmicRayParticlesThisPack(MeshBlockPack *ppack, const Real ppc) {
  auto &indcs = ppack->pmesh->mb_indcs;
  const int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
  const Real particles_per_mb = ppc*static_cast<Real>(ncells);
  int npart = 0;
  for (int m = 0; m < ppack->nmb_thispack; ++m) {
    const int gid = ppack->gids + m;
    const int nblock = ParticlesInGlobalMeshBlock(particles_per_mb, gid);
    if (nblock > std::numeric_limits<int>::max() - npart) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/ppc creates too many particles for an int counter"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    npart += nblock;
  }
  return npart;
}

KOKKOS_INLINE_FUNCTION
std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30))*0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27))*0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

KOKKOS_INLINE_FUNCTION
Real DeterministicUniform01(const int seed, const int gid, const int pinmb,
                            const int component) {
  std::uint64_t key = static_cast<std::uint64_t>(seed);
  key ^= (static_cast<std::uint64_t>(gid) + 0x9e3779b97f4a7c15ULL)*
         0xbf58476d1ce4e5b9ULL;
  key ^= (static_cast<std::uint64_t>(pinmb + 1) + 0x94d049bb133111ebULL)*
         0x9e3779b97f4a7c15ULL;
  key ^= (static_cast<std::uint64_t>(component + 1))*0xd2b74407b1ce6e93ULL;
  const std::uint64_t bits = SplitMix64(key);
  return static_cast<Real>(bits >> 11)*
         static_cast<Real>(1.0/9007199254740992.0);
}

void ValidateMomentBoundaryPolicy(Mesh *pmesh) {
  std::array<BoundaryFace, 6> faces = {
      BoundaryFace::inner_x1, BoundaryFace::outer_x1,
      BoundaryFace::inner_x2, BoundaryFace::outer_x2,
      BoundaryFace::inner_x3, BoundaryFace::outer_x3};
  int nfaces = 2;
  if (pmesh->multi_d) nfaces += 2;
  if (pmesh->three_d) nfaces += 2;

  for (int n = 0; n < nfaces; ++n) {
    const auto face = faces[n];
    BoundaryFlag flag = pmesh->mesh_bcs[face];
    if (!MomentBoundaryFlagSupported(flag)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/deposit_moments=true does not support "
                << BoundaryFaceName(face) << "=" << BoundaryFlagName(flag)
                << " in PR3c (supported: periodic/inflow/outflow/reflect)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

void ValidateDirectEdgeCurrentBoundaryPolicy(Mesh *pmesh) {
  std::array<BoundaryFace, 6> faces = {
      BoundaryFace::inner_x1, BoundaryFace::outer_x1,
      BoundaryFace::inner_x2, BoundaryFace::outer_x2,
      BoundaryFace::inner_x3, BoundaryFace::outer_x3};
  int nfaces = 2;
  if (pmesh->multi_d) nfaces += 2;
  if (pmesh->three_d) nfaces += 2;

  for (int n = 0; n < nfaces; ++n) {
    const auto face = faces[n];
    BoundaryFlag flag = pmesh->mesh_bcs[face];
    if (!DirectEdgeCurrentBoundaryFlagSupported(flag)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_j_deposition_mode=direct_staggered does not "
                << "support " << BoundaryFaceName(face) << "="
                << BoundaryFlagName(flag)
                << " in PR4 (supported: periodic/outflow/reflect; use "
                << "couple_j_deposition_mode=cc_convert for inflow boundaries)"
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
      if (ppc < 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<particles>/ppc must be >= 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      nprtcl_thispack = CountCosmicRayParticlesThisPack(pmy_pack, ppc);
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
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Unsupported value for <particles>/particle_type: "
                << ptype << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (particle_type == ParticleType::star && pmy_pack->punit == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/particle_type=star requires a <units> block."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("rk4_gravity") == 0) {
      pusher = ParticlesPusher::rk4_gravity;
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

  if (particle_type == ParticleType::star &&
      (pusher == ParticlesPusher::boris_lin ||
       pusher == ParticlesPusher::boris_tsc)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Boris pushers are incompatible with star particles; use "
              << "<particles>/particle_type=cosmic_ray for boris_lin/boris_tsc, "
              << "or use drift/rk4_gravity for star particles."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (particle_type == ParticleType::cosmic_ray &&
      pusher == ParticlesPusher::rk4_gravity) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pusher=rk4_gravity requires "
              << "<particles>/particle_type=star; cosmic-ray particles support "
              << "drift, boris_lin, or boris_tsc."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pusher == ParticlesPusher::rk4_gravity) {
    // Gravity constants are only meaningful for the star-particle RK4 pusher.
    r_scale   = pin->GetReal("potential", "r_scale");
    rho_scale = pin->GetReal("potential", "rho_scale");
    m_gal     = pin->GetReal("potential", "mass_gal");
    a_gal     = pin->GetReal("potential", "scale_gal");
    z_gal     = pin->GetReal("potential", "z_gal");
    r_200     = pin->GetReal("potential", "r_200");
    rho_mean  = pin->GetReal("potential", "rho_mean");
    par_grav_dx = pin->GetOrAddReal("particles", "grav_dx", 1e-6);
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
        // CR payload now always carries sampled B/cE and per-step feedback deltas.
        track_displacement = pin->GetOrAddBoolean("particles","track_displacement",false);
        nrdata = IPT_BIRTH + 1;

        nidata = 4;  // PGID, PTAG, PSP (species), PCRSOURCE
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

  // A publication run selects one coherent physical model. The engineering
  // default preserves the historical independent-toggle interface.
  std::string pic_physical_mode_str = pin->GetOrAddString(
      "particles", "pic_physical_mode", "engineering");
  if (pic_physical_mode_str.compare("engineering") == 0) {
    pic_physical_mode = PICPhysicalMode::engineering;
  } else if (pic_physical_mode_str.compare("paper_test_particle") == 0) {
    pic_physical_mode = PICPhysicalMode::paper_test_particle;
  } else if (pic_physical_mode_str.compare("paper_mhd_pic") == 0) {
    pic_physical_mode = PICPhysicalMode::paper_mhd_pic;
  } else if (pic_physical_mode_str.compare("extended_mhd_pic") == 0) {
    pic_physical_mode = PICPhysicalMode::extended_mhd_pic;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_physical_mode: "
              << pic_physical_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool paper_coupled_mode =
      (pic_physical_mode == PICPhysicalMode::paper_mhd_pic);

  // PR1 deposition controls
  deposit_moments = pin->GetOrAddBoolean("particles", "deposit_moments",
                                         paper_coupled_mode);
  deposit_order = pin->GetOrAddInteger("particles", "deposit_order", 1);
  deposit_qscale = pin->GetOrAddReal("particles", "deposit_qscale", 1.0);
  couple_moments_to_mhd = pin->GetOrAddBoolean("particles",
                                               "couple_moments_to_mhd",
                                               paper_coupled_mode);
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
      "particles", "couple_moments_momentum_to_mhd", paper_coupled_mode);
  couple_moments_energy_to_mhd = pin->GetOrAddBoolean(
      "particles", "couple_moments_energy_to_mhd", paper_coupled_mode);
  couple_moments_momentum_coeff = pin->GetOrAddReal(
      "particles", "couple_moments_momentum_coeff", 1.0);
  couple_moments_energy_coeff = pin->GetOrAddReal(
      "particles", "couple_moments_energy_coeff", 1.0);
  cr_vx0 = pin->GetOrAddReal("particles", "cr_vx0", 0.0);
  cr_vy0 = pin->GetOrAddReal("particles", "cr_vy0", 0.0);
  cr_vz0 = pin->GetOrAddReal("particles", "cr_vz0", 0.0);

  // Staged PR5+ PIC runtime controls (parse + validation only at this step)
  std::string pic_cr_hall_mode_str = pin->GetOrAddString(
      "particles", "pic_cr_hall_mode", "off");
  if (pic_cr_hall_mode_str.compare("off") == 0) {
    pic_cr_hall_mode = PICCRHallMode::off;
  } else if (pic_cr_hall_mode_str.compare("current_to_ct_experimental") == 0) {
    pic_cr_hall_mode = PICCRHallMode::current_to_ct_experimental;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_cr_hall_mode: "
              << pic_cr_hall_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string pic_wave_damping_mode_str = pin->GetOrAddString(
      "particles", "pic_wave_damping_mode", "off");
  if (pic_wave_damping_mode_str.compare("off") == 0) {
    pic_wave_damping_mode = PICWaveDampingMode::off;
  } else if (pic_wave_damping_mode_str.compare("ion_neutral_friction") == 0) {
    pic_wave_damping_mode = PICWaveDampingMode::ion_neutral_friction;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_wave_damping_mode: "
              << pic_wave_damping_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }
  pic_ion_neutral_collision_rate = pin->GetOrAddReal(
      "particles", "pic_ion_neutral_collision_rate", 0.0);
  if (pic_ion_neutral_collision_rate < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_ion_neutral_collision_rate must be >= 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  pic_enable_2d3v = pin->GetOrAddBoolean("particles", "pic_enable_2d3v", false);
  if ((pusher == ParticlesPusher::boris_lin || pusher == ParticlesPusher::boris_tsc) &&
      pmy_pack->pmesh->two_d && !pic_enable_2d3v) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Boris pushers in 2D require <particles>/pic_enable_2d3v=true; "
              << "the reduced 2D/2V Lorentz-force path is not implemented."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_cr_light_speed = pin->GetOrAddReal("particles", "pic_cr_light_speed", 1.0);
  if (pic_cr_light_speed <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_cr_light_speed must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const char *pic_cr_initial_state_default =
      (pic_physical_mode == PICPhysicalMode::engineering) ? "velocity" : "momentum";
  std::string pic_cr_initial_state_str = pin->GetOrAddString(
      "particles", "pic_cr_initial_state", pic_cr_initial_state_default);
  if (pic_cr_initial_state_str.compare("velocity") == 0) {
    pic_cr_initial_state = PICCRInitialState::velocity;
  } else if (pic_cr_initial_state_str.compare("momentum") == 0) {
    pic_cr_initial_state = PICCRInitialState::momentum;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_cr_initial_state: "
              << pic_cr_initial_state_str << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode == PICPhysicalMode::engineering) &&
      (pic_cr_initial_state != PICCRInitialState::velocity)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_physical_mode=engineering requires "
              << "<particles>/pic_cr_initial_state=velocity" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode == PICPhysicalMode::engineering) &&
      (std::abs(pic_cr_light_speed - 1.0) >
       16.0*std::numeric_limits<Real>::epsilon())) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_cr_light_speed is reserved in "
              << "<particles>/pic_physical_mode=engineering; select a momentum-state "
              << "physical mode to change particle mechanics" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_max_cell_cross = pin->GetOrAddInteger("particles", "pic_max_cell_cross", 2);
  if (pic_max_cell_cross <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_max_cell_cross must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  {
    const auto &indcs = pmy_pack->pmesh->mb_indcs;
    int min_mb_cells = indcs.nx1;
    if (pmy_pack->pmesh->multi_d) {
      min_mb_cells = std::min(min_mb_cells, indcs.nx2);
    }
    if (pmy_pack->pmesh->three_d) {
      min_mb_cells = std::min(min_mb_cells, indcs.nx3);
    }
    if (pic_max_cell_cross > min_mb_cells) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_max_cell_cross must not exceed the smallest "
                << "active MeshBlock dimension because particle exchange only "
                << "supports nearest-neighbor MeshBlock crossings."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
  } else if (pic_deltaf_mode_str.compare("quiet_start") == 0) {
    pic_deltaf_mode = PICDeltaFMode::quiet_start;
  } else if (pic_deltaf_mode_str.compare("on") == 0) {
    pic_deltaf_mode = (pic_physical_mode == PICPhysicalMode::engineering) ?
        PICDeltaFMode::quiet_start : PICDeltaFMode::physical;
  } else if (pic_deltaf_mode_str.compare("physical") == 0) {
    pic_deltaf_mode = PICDeltaFMode::physical;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_deltaf_mode: "
              << pic_deltaf_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_deltaf_f0 = pin->GetOrAddString("particles", "pic_deltaf_f0", "");
  if ((pic_deltaf_mode != PICDeltaFMode::off) && pic_deltaf_f0.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_deltaf_mode requires "
              << "<particles>/pic_deltaf_f0 to select a background distribution"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_deltaf_mode != PICDeltaFMode::off) &&
      !(pic_deltaf_f0.compare("uniform") == 0 ||
        pic_deltaf_f0.compare("uniform_quiet") == 0 ||
        pic_deltaf_f0.compare("kappa_iso") == 0 ||
        pic_deltaf_f0.compare("kappa_drift") == 0 ||
        pic_deltaf_f0.compare("kappa_aniso") == 0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_deltaf_f0: "
              << pic_deltaf_f0 << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pic_deltaf_f0.compare("kappa_iso") == 0) {
    pic_deltaf_background = PICDeltaFBackground::kappa_iso;
  } else if (pic_deltaf_f0.compare("kappa_drift") == 0) {
    pic_deltaf_background = PICDeltaFBackground::kappa_drift;
  } else if (pic_deltaf_f0.compare("kappa_aniso") == 0) {
    pic_deltaf_background = PICDeltaFBackground::kappa_aniso;
  } else {
    pic_deltaf_background = PICDeltaFBackground::uniform;
  }
  pic_deltaf_p0 = pin->GetOrAddReal("particles", "pic_deltaf_p0", 1.0);
  pic_deltaf_kappa = pin->GetOrAddReal("particles", "pic_deltaf_kappa", 1.25);
  pic_deltaf_drift_x1 = pin->GetOrAddReal("particles", "pic_deltaf_drift_x1", 0.0);
  pic_deltaf_drift_x2 = pin->GetOrAddReal("particles", "pic_deltaf_drift_x2", 0.0);
  pic_deltaf_drift_x3 = pin->GetOrAddReal("particles", "pic_deltaf_drift_x3", 0.0);
  pic_deltaf_aniso_x1 = pin->GetOrAddReal("particles", "pic_deltaf_aniso_x1", 1.0);
  pic_deltaf_aniso_x2 = pin->GetOrAddReal("particles", "pic_deltaf_aniso_x2", 1.0);
  pic_deltaf_aniso_x3 = pin->GetOrAddReal("particles", "pic_deltaf_aniso_x3", 1.0);
  pic_deltaf_background_rho = pin->GetOrAddReal(
      "particles", "pic_deltaf_background_rho", 0.0);
  pic_deltaf_background_jx = pin->GetOrAddReal(
      "particles", "pic_deltaf_background_jx", 0.0);
  pic_deltaf_background_jy = pin->GetOrAddReal(
      "particles", "pic_deltaf_background_jy", 0.0);
  pic_deltaf_background_jz = pin->GetOrAddReal(
      "particles", "pic_deltaf_background_jz", 0.0);
  std::string pic_deltaf_adapt_mode_str = pin->GetOrAddString(
      "particles", "pic_deltaf_adapt_mode", "off");
  if (pic_deltaf_adapt_mode_str.compare("off") == 0) {
    pic_deltaf_adapt_mode = PICDeltaFAdaptMode::off;
  } else if (pic_deltaf_adapt_mode_str.compare(
                 "global_bikappa_moments_experimental") == 0) {
    pic_deltaf_adapt_mode =
        PICDeltaFAdaptMode::global_bikappa_moments_experimental;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_deltaf_adapt_mode: "
              << pic_deltaf_adapt_mode_str << std::endl;
    std::exit(EXIT_FAILURE);
  }
  pic_deltaf_adapt_interval = pin->GetOrAddReal(
      "particles", "pic_deltaf_adapt_interval", 0.0);
  if (pic_deltaf_adapt_interval < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_deltaf_adapt_interval must be >= 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  pic_deltaf_adaptive_p0 = pic_deltaf_p0;
  if (UsesDeltaF() &&
      (pic_deltaf_p0 <= 0.0 || pic_deltaf_kappa <= 0.0 ||
       pic_deltaf_aniso_x1 <= 0.0 || pic_deltaf_aniso_x2 <= 0.0 ||
       pic_deltaf_aniso_x3 <= 0.0 || pic_deltaf_background_rho < 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Physical delta-f requires positive p0, kappa, anisotropy "
              << "scales and non-negative background charge density" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_sort_interval = pin->GetOrAddInteger("particles", "pic_sort_interval", 0);
  if (pic_sort_interval < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_sort_interval must be >= 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_load_balance_cost_per_particle = pin->GetOrAddReal(
      "particles", "pic_load_balance_cost_per_particle", 0.0);
  if (pic_load_balance_cost_per_particle < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_load_balance_cost_per_particle must be >= 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_q017_sync_kernel_timers = pin->GetOrAddBoolean(
      "particles", "pic_q017_sync_kernel_timers", false);

  pic_random_seed = pin->GetOrAddInteger("particles", "pic_random_seed", 0);
  if (pic_random_seed < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_random_seed must be >= 0" << std::endl;
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

  std::string pic_expansion_law_str = pin->GetOrAddString(
      "particles", "pic_expansion_law", "linear");
  if (pic_expansion_law_str.compare("linear") == 0) {
    pic_expansion_law = PICExpansionLaw::linear;
  } else if (pic_expansion_law_str.compare("reciprocal_linear") == 0) {
    pic_expansion_law = PICExpansionLaw::reciprocal_linear;
  } else if (pic_expansion_law_str.compare("exponential") == 0) {
    pic_expansion_law = PICExpansionLaw::exponential;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/pic_expansion_law: "
              << pic_expansion_law_str << std::endl;
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
  if (pic_expanding_box_mode == PICExpandingBoxMode::on) {
    const Real tstart = pmy_pack->pmesh->time;
    const Real tlim = pin->GetOrAddReal("time", "tlim", tstart);
    const std::array<Real, 2> times = {tstart, tlim};
    const std::array<Real, 3> rates = {
        pic_expansion_rate_x1, pic_expansion_rate_x2, pic_expansion_rate_x3};
    for (const Real time : times) {
      for (const Real rate : rates) {
        const Real scale = PICScaleFactor(pic_expansion_law, rate, time);
        if (!std::isfinite(scale) || scale <= 0.0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "<particles>/pic_expansion_law and rate must produce finite, "
                    << "positive scale factors through <time>/tlim" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  if (UsesExpandingBox() && pin->DoesBlockExist("shearing_box")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_expanding_box_mode=on does not support "
              << "<shearing_box> or orbital advection" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool uses_expanding_mhd = (UsesExpandingBox() && pin->DoesBlockExist("mhd"));
  if (uses_expanding_mhd) {
    auto reject_expanding_mhd = [](const char *reason) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_expanding_box_mode=on with active MHD "
                << "does not support " << reason << std::endl;
      std::exit(EXIT_FAILURE);
    };
    if (!pmy_pack->pmesh->strictly_periodic) {
      reject_expanding_mhd("non-periodic boundaries");
    }
    if (pmy_pack->pmesh->multilevel) {
      reject_expanding_mhd("SMR or AMR meshes");
    }
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic ||
        pmy_pack->pcoord->is_dynamical_relativistic) {
      reject_expanding_mhd("relativistic coordinate paths");
    }
    if (pin->GetString("mhd", "eos").compare("ideal") != 0) {
      reject_expanding_mhd("<mhd>/eos values other than ideal");
    }
    if (pin->GetString("time", "evolution").compare("dynamic") != 0) {
      reject_expanding_mhd("<time>/evolution values other than dynamic");
    }
    if ((particle_type != ParticleType::cosmic_ray) ||
        ((pusher != ParticlesPusher::boris_lin) &&
         (pusher != ParticlesPusher::boris_tsc))) {
      reject_expanding_mhd(
          "particle types or pushers other than cosmic-ray Boris pushers");
    }
    if (pin->GetReal("particles", "ppc") > 0.0) {
      const Real tstart = pmy_pack->pmesh->time;
      const Real tlim = pin->GetReal("time", "tlim");
      const auto geom_start = PICExpandingBoxGeometryAt(
          pic_expansion_law, pic_expansion_rate_x1, pic_expansion_rate_x2,
          pic_expansion_rate_x3, tstart);
      const auto geom_end = PICExpandingBoxGeometryAt(
          pic_expansion_law, pic_expansion_rate_x1, pic_expansion_rate_x2,
          pic_expansion_rate_x3, tlim);
      if ((geom_end.a1 < geom_start.a1) || (geom_end.a2 < geom_start.a2) ||
          (geom_end.a3 < geom_start.a3)) {
        reject_expanding_mhd("contracting scale factors with active particles");
      }
    }
    if ((pmy_pack->pmhd != nullptr) && pmy_pack->pmhd->use_fofc) {
      reject_expanding_mhd("<mhd>/fofc=true");
    }
    if (pin->DoesParameterExist("mhd", "viscosity") ||
        pin->DoesParameterExist("mhd", "ohmic_resistivity") ||
        pin->DoesParameterExist("mhd", "conductivity") ||
        pin->DoesParameterExist("mhd", "tdep_conductivity")) {
      reject_expanding_mhd("MHD viscosity, resistivity, or conductivity");
    }
    if (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
        pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("adm") ||
        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("turb_driving") ||
        pin->DoesBlockExist("initial_turb")) {
      reject_expanding_mhd("coupled fluid, radiation, relativity, or turbulence blocks");
    }
    if (couple_moments_to_mhd &&
        ((couple_j_to_efield_representation ==
          CoupledCurrentRepresentation::edge_staggered) ||
         (couple_j_deposition_mode ==
          CoupledCurrentDepositionMode::direct_staggered))) {
      reject_expanding_mhd("staggered or direct CR-current deposition");
    }
    if (UsesDeltaF() &&
        ((pic_deltaf_background_jx != 0.0) ||
         (pic_deltaf_background_jy != 0.0) ||
         (pic_deltaf_background_jz != 0.0))) {
      reject_expanding_mhd("nonzero delta-f background current");
    }
    const bool uses_coupled_feedback =
        (couple_moments_to_mhd || couple_moments_momentum_to_mhd ||
         couple_moments_energy_to_mhd);
    const bool common_coupled_feedback =
        (deposit_moments && couple_moments_to_mhd &&
         couple_moments_momentum_to_mhd && couple_moments_energy_to_mhd &&
         (pic_background_mode == PICBackgroundMode::coupled) &&
         (pic_feedback_mode == PICFeedbackMode::coupled) &&
         (couple_fluid_feedback_order == CoupledFluidFeedbackOrder::mhd_src_terms) &&
         (couple_j_to_efield_representation ==
          CoupledCurrentRepresentation::cell_centered) &&
         (couple_j_deposition_mode == CoupledCurrentDepositionMode::cc_convert) &&
         (pic_cr_hall_mode == PICCRHallMode::off));
    const bool qualified_conservative_feedback =
        (common_coupled_feedback && !UsesDeltaF());
    const bool experimental_adaptive_deltaf_feedback =
        (common_coupled_feedback && UsesAdaptiveDeltaF() &&
         (pic_physical_mode == PICPhysicalMode::extended_mhd_pic));
    const bool admitted_coupled_feedback =
        (qualified_conservative_feedback || experimental_adaptive_deltaf_feedback);
    if (uses_coupled_feedback && !admitted_coupled_feedback) {
      reject_expanding_mhd("particle feedback outside the admitted cell-centered "
                           "source splits");
    }
    if (UsesPICWaveDamping() && !admitted_coupled_feedback) {
      reject_expanding_mhd("reduced ion-neutral friction without an admitted "
                           "cell-centered source split");
    }
    if (pic_background_mode == PICBackgroundMode::no_mhd) {
      reject_expanding_mhd("<particles>/pic_background_mode=no_mhd with an active "
                           "<mhd> block");
    }
    if (pin->DoesParameterExist("problem", "user_hist") &&
        pin->GetBoolean("problem", "user_hist")) {
      reject_expanding_mhd("user-defined history callbacks");
    }
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
    if (direct_edge_mode) {
      ValidateDirectEdgeCurrentBoundaryPolicy(pmy_pack->pmesh);
    }
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
  if ((pic_physical_mode != PICPhysicalMode::engineering) &&
      (particle_type != ParticleType::cosmic_ray)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_physical_mode=" << pic_physical_mode_str
              << " requires <particles>/particle_type=cosmic_ray" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode != PICPhysicalMode::engineering) &&
      (pusher != ParticlesPusher::boris_tsc)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_physical_mode=" << pic_physical_mode_str
              << " requires <particles>/pusher=boris_tsc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode != PICPhysicalMode::extended_mhd_pic) &&
      (pic_cr_hall_mode != PICCRHallMode::off)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_cr_hall_mode=" << pic_cr_hall_mode_str
              << " requires <particles>/pic_physical_mode=extended_mhd_pic"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode != PICPhysicalMode::extended_mhd_pic) &&
      (pic_wave_damping_mode != PICWaveDampingMode::off)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_wave_damping_mode=" << pic_wave_damping_mode_str
              << " requires <particles>/pic_physical_mode=extended_mhd_pic"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_physical_mode != PICPhysicalMode::extended_mhd_pic) &&
      UsesAdaptiveDeltaF()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_deltaf_adapt_mode=" << pic_deltaf_adapt_mode_str
              << " requires <particles>/pic_physical_mode=extended_mhd_pic"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (UsesAdaptiveDeltaF() &&
      (!UsesDeltaF() ||
       (pic_deltaf_background != PICDeltaFBackground::kappa_aniso) ||
       (pic_expanding_box_mode != PICExpandingBoxMode::on) ||
       (pic_deltaf_adapt_interval <= 0.0) ||
       (pic_deltaf_kappa <= 1.0) ||
       (pic_deltaf_drift_x1 != 0.0) || (pic_deltaf_drift_x2 != 0.0) ||
       (pic_deltaf_drift_x3 != 0.0) ||
       (pic_deltaf_aniso_x1 != 1.0) || (pic_deltaf_aniso_x2 != 1.0) ||
       (pic_deltaf_aniso_x3 != 1.0))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_deltaf_adapt_mode="
              << "global_bikappa_moments_experimental requires physical "
              << "kappa_aniso delta-f, expanding-box mode, positive adapt interval, "
              << "kappa > 1, zero drift, and unit configured anisotropy scales"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction) &&
      (pic_ion_neutral_collision_rate <= 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_wave_damping_mode=ion_neutral_friction requires "
              << "<particles>/pic_ion_neutral_collision_rate > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction) &&
      ((pic_background_mode != PICBackgroundMode::coupled) ||
       !(pin->DoesBlockExist("mhd")))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_wave_damping_mode=ion_neutral_friction requires "
              << "an active coupled <mhd> background" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction) {
    auto reject_wave_damping_composition = [](const char *reason) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_wave_damping_mode=ion_neutral_friction requires "
                << "the Newtonian single-fluid MHD source task path; does not support "
                << reason << std::endl;
      std::exit(EXIT_FAILURE);
    };
    if (pin->DoesBlockExist("radiation")) {
      reject_wave_damping_composition("<radiation> task paths");
    }
    if (pin->DoesBlockExist("ion-neutral")) {
      reject_wave_damping_composition("<ion-neutral> alternate task lists");
    }
    if (pin->DoesBlockExist("hydro")) {
      reject_wave_damping_composition("<hydro> compositions");
    }
    if (pin->DoesBlockExist("adm") || pin->DoesBlockExist("z4c")) {
      reject_wave_damping_composition("dynamical-GR <adm>/<z4c> task paths");
    }
    if (pmy_pack->pcoord->is_special_relativistic) {
      reject_wave_damping_composition("<coord>/special_rel=true");
    }
    if (pmy_pack->pcoord->is_general_relativistic) {
      reject_wave_damping_composition("<coord>/general_rel=true");
    }
    if (pmy_pack->pcoord->is_dynamical_relativistic) {
      reject_wave_damping_composition("dynamical-relativistic coordinates");
    }
  }
  if (pic_physical_mode == PICPhysicalMode::paper_test_particle) {
    if ((pic_feedback_mode != PICFeedbackMode::test_particle) ||
        couple_moments_to_mhd || couple_moments_momentum_to_mhd ||
        couple_moments_energy_to_mhd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_physical_mode=paper_test_particle requires "
                << "test-particle feedback with all particle-to-MHD coupling "
                << "toggles disabled" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (pic_physical_mode == PICPhysicalMode::paper_mhd_pic) {
    const bool exact_isothermal_deltaf_paper_feedback =
        (pin->GetString("mhd", "eos").compare("isothermal") == 0) &&
        UsesDeltaF() && !couple_moments_energy_to_mhd;
    if ((pic_background_mode != PICBackgroundMode::coupled) ||
        (pic_feedback_mode != PICFeedbackMode::coupled) ||
        !deposit_moments || !couple_moments_to_mhd ||
        !couple_moments_momentum_to_mhd ||
        (!couple_moments_energy_to_mhd &&
         !exact_isothermal_deltaf_paper_feedback)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_physical_mode=paper_mhd_pic requires coupled "
                << "MHD background, coupled feedback, moment deposition, and "
                << "conservative momentum feedback. Energy feedback is required "
                << "for ideal MHD; exact isothermal paper delta-f uses momentum-only "
                << "feedback with energy feedback disabled." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((couple_j_to_efield_representation ==
         CoupledCurrentRepresentation::edge_staggered) ||
        (couple_j_deposition_mode ==
         CoupledCurrentDepositionMode::direct_staggered)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/pic_physical_mode=paper_mhd_pic rejects "
                << "direct-current CT induction options; use the ideal-MHD "
                << "paper induction path" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if ((pic_cr_hall_mode == PICCRHallMode::current_to_ct_experimental) &&
      !couple_moments_to_mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_cr_hall_mode=current_to_ct_experimental "
              << "requires <particles>/couple_moments_to_mhd=true" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if ((global_variable::my_rank == 0) &&
      (particle_type == ParticleType::cosmic_ray)) {
    const char *state_name = UsesRelativisticCRState() ? "momentum_p_over_m" :
                                                        "velocity";
    const char *induction_name = AddsCRCurrentToCT() ? "cr_current_to_ct" :
                                                       "ideal_mhd_only";
    std::cout << "PIC runtime model: physical_mode=" << pic_physical_mode_str
              << " state=" << state_name
              << " C=" << pic_cr_light_speed
              << " background=" << pic_background_mode_str
              << " feedback=" << pic_feedback_mode_str
              << " induction=" << induction_name
              << " deposition=tsc"
              << " deltaf=" << pic_deltaf_mode_str
              << " deltaf_adapt=" << pic_deltaf_adapt_mode_str
              << " deltaf_adapt_interval=" << pic_deltaf_adapt_interval
              << " expanding_box=" << pic_expanding_box_mode_str
              << " expansion_law=" << pic_expansion_law_str
              << " wave_damping=" << pic_wave_damping_mode_str
              << " nu_in=" << pic_ion_neutral_collision_rate
              << " lb_cost_per_particle=" << pic_load_balance_cost_per_particle
              << " max_cell_cross=" << pic_max_cell_cross
              << " theta_max=" << pic_theta_max
              << " restart_schema=" << PIC_RESTART_SCHEMA_VERSION << std::endl;
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
    if (!(pmy_pack->pmesh->strictly_periodic) &&
        (std::getenv("ATHENA_SKIP_MOM_INFLOW_ZERO") == nullptr)) {
      // Use zero-valued inflow moments unless a problem callback overwrites them.
      // This keeps coupled inflow boundaries well-defined.
      Kokkos::deep_copy(pbval_mom->u_in.d_view, static_cast<Real>(0.0));
      pbval_mom->u_in.template modify<DevExeSpace>();
      pbval_mom->u_in.template sync<HostMemSpace>();
    }

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
    int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
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
//! \fn void Particles::UpdateAfterAMR
//! \brief Refresh the retained pack pointer and validate particle-owned AMR capacity.

void Particles::UpdateAfterAMR(MeshBlockPack *new_pp) {
  pmy_pack = new_pp;
  const int nmb = pmy_pack->nmb_thispack;
  auto require_meshblock_capacity = [nmb](const int allocated, const char *name) {
    if (allocated < nmb) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle AMR refresh found insufficient " << name
                << " MeshBlock capacity: allocated=" << allocated
                << " required=" << nmb << std::endl;
      std::exit(EXIT_FAILURE);
    }
  };

  if (deposit_moments) {
    require_meshblock_capacity(moments.extent_int(0), "moment");
    if (pmy_pack->pmesh->multilevel) {
      require_meshblock_capacity(coarse_moments.extent_int(0), "coarse-moment");
    }
    if (x1_old.extent_int(0) < nprtcl_thispack) {
      Kokkos::resize(x1_old, nprtcl_thispack);
      Kokkos::resize(x2_old, nprtcl_thispack);
      Kokkos::resize(x3_old, nprtcl_thispack);
    }
    if (j_edge_x1e.extent_int(0) > 0) {
      require_meshblock_capacity(j_edge_x1e.extent_int(0), "edge-current");
    }
  }
  if (pic_background_mode == PICBackgroundMode::no_mhd) {
    require_meshblock_capacity(pic_no_mhd_bcc0.extent_int(0), "no-MHD-field");
  }
  // Boundary helpers retain this Particles object or the stable MeshBlockPack pointer.
  // Their kernels resolve the reconstructed MeshBlock and neighbor state dynamically.
}

//----------------------------------------------------------------------------------------
// InitializeCosmicRays()
// Initializes cosmic ray particles with species support

void Particles::InitializeCosmicRays(ParameterInput *pin) {
  // Read number of species
  nspecies = pin->GetOrAddInteger("particles","nspecies",1);
  if (nspecies <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/nspecies must be > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Allocate species arrays
  Kokkos::realloc(species_mass, nspecies);
  Kokkos::realloc(species_charge, nspecies);
  Kokkos::realloc(species_vx0, nspecies);
  Kokkos::realloc(species_vy0, nspecies);
  Kokkos::realloc(species_vz0, nspecies);

  // Read species properties
  auto h_mass = Kokkos::create_mirror_view(species_mass);
  auto h_charge = Kokkos::create_mirror_view(species_charge);
  auto h_vx0 = Kokkos::create_mirror_view(species_vx0);
  auto h_vy0 = Kokkos::create_mirror_view(species_vy0);
  auto h_vz0 = Kokkos::create_mirror_view(species_vz0);

  for (int s=0; s<nspecies; ++s) {
    std::string block = "species" + std::to_string(s);
    h_mass(s) = pin->GetOrAddReal(block,"mass",1.0);
    h_charge(s) = pin->GetOrAddReal(block,"charge",1.0);
    h_vx0(s) = pin->GetOrAddReal(block, "vx0", cr_vx0);
    h_vy0(s) = pin->GetOrAddReal(block, "vy0", cr_vy0);
    h_vz0(s) = pin->GetOrAddReal(block, "vz0", cr_vz0);
    if (h_mass(s) <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/mass must be > 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (UsesRelativisticCRState() &&
        (pic_cr_initial_state == PICCRInitialState::velocity)) {
      const Real v2 = h_vx0(s)*h_vx0(s) + h_vy0(s)*h_vy0(s) + h_vz0(s)*h_vz0(s);
      if (v2 >= pic_cr_light_speed*pic_cr_light_speed) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<" << block << ">/vx0,vy0,vz0 must define |v| < "
                  << "<particles>/pic_cr_light_speed when "
                  << "<particles>/pic_cr_initial_state=velocity" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  Kokkos::deep_copy(species_mass, h_mass);
  Kokkos::deep_copy(species_charge, h_charge);
  Kokkos::deep_copy(species_vx0, h_vx0);
  Kokkos::deep_copy(species_vy0, h_vy0);
  Kokkos::deep_copy(species_vz0, h_vz0);

  // Initialize particle positions and velocities
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  int nmb = pmy_pack->nmb_thispack;
  const int gids_local = pmy_pack->gids;

  // Determine distribution type for initial positions
  std::string dist = pin->GetOrAddString("particles", "cr_distribution", "center");
  const bool center_dist = (dist.compare("center") == 0);
  const bool random_dist = (dist.compare("random") == 0);
  if (!center_dist && !random_dist) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported value for <particles>/cr_distribution: "
              << dist << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool deltaf_quiet_start =
      (random_dist && pic_deltaf_mode != PICDeltaFMode::off);
  const bool track_disp_local = track_displacement;  // avoid capturing 'this'
  const int nx1_local = indcs.nx1;                   // avoid capturing host refs
  const int nx2_local = indcs.nx2;
  const int nx3_local = indcs.nx3;                   // avoid capturing host ref

  const Real ppc = pin->GetOrAddReal("particles", "ppc", 1.0);
  const int ncells = nx1_local*nx2_local*nx3_local;
  const Real particles_per_mb_real = ppc*static_cast<Real>(ncells);
  DvceArray1D<int> block_offsets("cr_block_offsets", nmb + 1);
  HostArray1D<int> h_block_offsets("cr_block_offsets_host", nmb + 1);
  h_block_offsets(0) = 0;
  for (int m = 0; m < nmb; ++m) {
    const int gid = gids_local + m;
    h_block_offsets(m + 1) =
        h_block_offsets(m) + ParticlesInGlobalMeshBlock(particles_per_mb_real, gid);
  }
  if (h_block_offsets(nmb) != nprtcl_thispack) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Internal PIC particle-count mismatch during initialization."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Kokkos::deep_copy(block_offsets, h_block_offsets);

  // Create local copies to avoid GPU lambda capture warnings
  int nspecies_local = nspecies;
  auto species_charge_local = species_charge;
  auto species_mass_local = species_mass;
  auto species_vx0_local = species_vx0;
  auto species_vy0_local = species_vy0;
  auto species_vz0_local = species_vz0;
  const bool deltaf_quiet_start_local = deltaf_quiet_start;
  const int pic_random_seed_local = pic_random_seed;
  const bool momentum_state_local = UsesRelativisticCRState();
  const bool initialize_from_velocity =
      (pic_cr_initial_state == PICCRInitialState::velocity);
  const Real light_speed_local = pic_cr_light_speed;
  const PICDeltaFBackground deltaf_background_local = pic_deltaf_background;
  const bool adaptive_deltaf_local = UsesAdaptiveDeltaF();
  const Real deltaf_p0_local = pic_deltaf_p0;
  const Real deltaf_kappa_local = pic_deltaf_kappa;
  const Real deltaf_drift_x1_local = pic_deltaf_drift_x1;
  const Real deltaf_drift_x2_local = pic_deltaf_drift_x2;
  const Real deltaf_drift_x3_local = pic_deltaf_drift_x3;
  const Real deltaf_aniso_x1_local = pic_deltaf_aniso_x1;
  const Real deltaf_aniso_x2_local = pic_deltaf_aniso_x2;
  const Real deltaf_aniso_x3_local = pic_deltaf_aniso_x3;
  const Real birth_time_local = pmy_pack->pmesh->time;
  // Make sure geometry is available on device
  size.template sync<DevExeSpace>();
  auto size_view = size;
  // Capture handle by value for use as size_view.d_view(...) on device.
  const Real root_cell_vol =
      pmy_pack->pmesh->mesh_size.dx1*
      pmy_pack->pmesh->mesh_size.dx2*
      pmy_pack->pmesh->mesh_size.dx3;

  // Simple uniform distribution for now
  par_for("init_cr", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(const int p) {
    // Determine which meshblock this particle belongs to
    int m = 0;
    while ((m + 1) < nmb && p >= block_offsets(m + 1)) {
      ++m;
    }
    const int pinmb = p - block_offsets(m);
    const int position_index = pinmb/nspecies_local;
    const int gid = gids_local + m;

    // Set GID and species
    pi(PGID,p) = gid;
    pi(PTAG,p) = p;
    pi(PSP,p) = pinmb % nspecies_local;  // Round-robin species assignment per block
    pi(PCRSOURCE,p) = static_cast<int>(CRParticleSource::initial);

    // Choose position within the mesh block
    Real rx = 0.5, ry = 0.5, rz = 0.5;
    if (center_dist) {
      const int particles_this_block = block_offsets(m + 1) - block_offsets(m);
      const int positions_this_block =
          (particles_this_block + nspecies_local - 1)/nspecies_local;
      if (positions_this_block > 0 && ncells > 0) {
        int cell_linear;
        if (positions_this_block <= ncells) {
          cell_linear = static_cast<int>(
              (static_cast<Real>(position_index) + 0.5)*
              static_cast<Real>(ncells)/static_cast<Real>(positions_this_block));
          if (cell_linear >= ncells) cell_linear = ncells - 1;
        } else {
          cell_linear = position_index % ncells;
        }

        const int ci = cell_linear % nx1_local;
        int rem = cell_linear/nx1_local;
        const int cj = rem % nx2_local;
        const int ck = rem/nx2_local;
        rx = (static_cast<Real>(ci) + 0.5)/static_cast<Real>(nx1_local);
        ry = (static_cast<Real>(cj) + 0.5)/static_cast<Real>(nx2_local);
        rz = (nx3_local > 1) ?
            (static_cast<Real>(ck) + 0.5)/static_cast<Real>(nx3_local) : 0.5;
      }
    } else if (random_dist) {
      if (deltaf_quiet_start_local) {
        // Low-discrepancy placement for staged quiet-start noise control.
        int n2 = position_index + 1;
        Real w2 = 0.5;
        rx = 0.0;
        while (n2 > 0) {
          rx += w2*static_cast<Real>(n2 % 2);
          n2 /= 2;
          w2 *= 0.5;
        }

        int n3 = position_index + 1;
        Real w3 = 1.0/3.0;
        ry = 0.0;
        while (n3 > 0) {
          ry += w3*static_cast<Real>(n3 % 3);
          n3 /= 3;
          w3 *= 1.0/3.0;
        }

        if (nx3_local > 1) {
          int n5 = position_index + 1;
          Real w5 = 0.2;
          rz = 0.0;
          while (n5 > 0) {
            rz += w5*static_cast<Real>(n5 % 5);
            n5 /= 5;
            w5 *= 0.2;
          }
        }
      } else {
        rx = DeterministicUniform01(pic_random_seed_local, gid, position_index, 0);
        ry = DeterministicUniform01(pic_random_seed_local, gid, position_index, 1);
        if (nx3_local > 1) {
          rz = DeterministicUniform01(pic_random_seed_local, gid, position_index, 2);
        }
      }
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

    // Set species-dependent initialization values
    int species = pi(PSP,p);
    Real state_x = species_vx0_local(species);
    Real state_y = species_vy0_local(species);
    Real state_z = species_vz0_local(species);
    if (momentum_state_local && initialize_from_velocity) {
      const Real v2 = state_x*state_x + state_y*state_y + state_z*state_z;
      const Real gamma = 1.0/sqrt(1.0 - v2/(light_speed_local*light_speed_local));
      state_x *= gamma;
      state_y *= gamma;
      state_z *= gamma;
    }
    pr(IPVX,p) = state_x;
    pr(IPVY,p) = state_y;
    pr(IPVZ,p) = state_z;
    pr(IPM,p) = species_charge_local(species) / species_mass_local(species);

    // Initialize B-field components to zero
    pr(IPBX,p) = 0.0;
    pr(IPBY,p) = 0.0;
    pr(IPBZ,p) = 0.0;
    pr(IPEX,p) = 0.0;
    pr(IPEY,p) = 0.0;
    pr(IPEZ,p) = 0.0;
    pr(IPDPX,p) = 0.0;
    pr(IPDPY,p) = 0.0;
    pr(IPDPZ,p) = 0.0;
    pr(IPDE,p) = 0.0;
    pr(IPEBDOT,p) = 0.0;
    pr(IPWT,p) = (size_view.d_view(m).dx1*
                  size_view.d_view(m).dx2*
                  size_view.d_view(m).dx3)/root_cell_vol;
    if (adaptive_deltaf_local) {
      pr(IPF0,p) = PICAdaptiveDeltaFBackgroundValue(
          deltaf_kappa_local, deltaf_p0_local, deltaf_p0_local, 1.0,
          1.0, 1.0, 1.0, state_x, state_y, state_z);
    } else {
      pr(IPF0,p) = PICDeltaFBackgroundValue(
          deltaf_background_local, deltaf_p0_local, deltaf_kappa_local,
          deltaf_drift_x1_local, deltaf_drift_x2_local, deltaf_drift_x3_local,
          deltaf_aniso_x1_local, deltaf_aniso_x2_local, deltaf_aniso_x3_local,
          1.0, 1.0, 1.0, state_x, state_y, state_z);
    }
    pr(IPDFWT,p) = 0.0;
    pr(IPT_BIRTH,p) = birth_time_local;

    // Initialize displacement tracking if enabled
    if (track_disp_local) {
      pr(IPDX,p) = 0.0;
      pr(IPDY,p) = 0.0;
      pr(IPDZ,p) = 0.0;
      pr(IPDB,p) = 0.0;
    }
  });
  Kokkos::fence();

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
// NewTimeStep()
// Calculates the particle timestep limit.

void Particles::NewTimeStep() {
  const Real max_dt = std::numeric_limits<Real>::max();
  if (particle_type != ParticleType::cosmic_ray) {
    auto &size = pmy_pack->pmb->mb_size;
    size.template sync<HostMemSpace>();
    Real dt_part = max_dt;
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      Real dx_min = size.h_view(m).dx1;
      if (pmy_pack->pmesh->multi_d) {
        dx_min = std::min(dx_min, size.h_view(m).dx2);
      }
      if (pmy_pack->pmesh->three_d) {
        dx_min = std::min(dx_min, size.h_view(m).dx3);
      }
      dt_part = std::min(dt_part, dx_min);
    }
    dtnew = dt_part;
    return;
  }

  if (nprtcl_thispack <= 0) {
    dtnew = max_dt;
    return;
  }

  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int gids = pmy_pack->gids;
  const int nmb = pmy_pack->nmb_thispack;
  const Real max_cell_cross = static_cast<Real>(pic_max_cell_cross);
  const bool momentum_state_local = UsesRelativisticCRState();
  const Real light_speed_local = pic_cr_light_speed;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;

  size.template sync<DevExeSpace>();
  auto size_view = size;

  Real dt_part = max_dt;
  Kokkos::parallel_reduce(
      "ParticlesNewTimeStep", Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
      KOKKOS_LAMBDA(const int &p, Real &min_dt) {
        int m = pi(PGID, p) - gids;
        if (m < 0 || m >= nmb) return;

        Real candidate = max_dt;
        Real vx, vy, vz;
        CRVelocityFromState(momentum_state_local, light_speed_local,
                            pr(IPVX, p), pr(IPVY, p), pr(IPVZ, p), vx, vy, vz);
        vx = fabs(vx);
        vy = fabs(vy);
        vz = fabs(vz);

        if (vx > 0.0) {
          candidate = fmin(candidate, max_cell_cross*size_view.d_view(m).dx1/vx);
        }
        if (multi_d && vy > 0.0) {
          candidate = fmin(candidate, max_cell_cross*size_view.d_view(m).dx2/vy);
        }
        if (three_d && vz > 0.0) {
          candidate = fmin(candidate, max_cell_cross*size_view.d_view(m).dx3/vz);
        }
        min_dt = fmin(min_dt, candidate);
      },
      Kokkos::Min<Real>(dt_part));

  if (pusher == ParticlesPusher::boris_lin || pusher == ParticlesPusher::boris_tsc) {
    Real qom_max = 0.0;
    Kokkos::parallel_reduce(
        "ParticlesNewTimeStepQom",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
        KOKKOS_LAMBDA(const int &p, Real &max_qom) {
          max_qom = fmax(max_qom, fabs(pr(IPM, p)));
        },
        Kokkos::Max<Real>(qom_max));

    DvceArray5D<Real> bcc;
    if (pic_background_mode == PICBackgroundMode::no_mhd) {
      bcc = pic_no_mhd_bcc0;
    } else if (pmy_pack->pmhd != nullptr) {
      bcc = pmy_pack->pmhd->bcc0;
    }

    if (bcc.size() > 0 && qom_max > 0.0) {
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      const int is = indcs.is;
      const int js = indcs.js;
      const int ks = indcs.ks;
      const int nx1 = indcs.nx1;
      const int nx2 = indcs.nx2;
      const int nx3 = indcs.nx3;
      const int nkji = nx3*nx2*nx1;
      const int nji = nx2*nx1;
      const int nmkji = nmb*nkji;
      Real b2_max = 0.0;

      Kokkos::parallel_reduce(
          "ParticlesNewTimeStepBmax",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
          KOKKOS_LAMBDA(const int &idx, Real &max_b2) {
            const int m = idx/nkji;
            int rem = idx - m*nkji;
            int k = rem/nji;
            rem -= k*nji;
            int j = rem/nx1;
            int i = rem - j*nx1;
            k += ks;
            j += js;
            i += is;

            const Real bx = bcc(m, IBX, k, j, i);
            const Real by = bcc(m, IBY, k, j, i);
            const Real bz = bcc(m, IBZ, k, j, i);
            max_b2 = fmax(max_b2, bx*bx + by*by + bz*bz);
          },
          Kokkos::Max<Real>(b2_max));

      if (b2_max > 0.0) {
        const Real omega_max = qom_max*std::sqrt(b2_max);
        dt_part = std::min(dt_part, pic_theta_max/omega_max);
      }
    }
  }

  dtnew = dt_part;
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
