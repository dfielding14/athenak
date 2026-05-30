//========================================================================================
// Athena++ astrophysical MHD code, Kokkos version
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_random.cpp
//! \brief Problem generator that initializes random particle positions and velocities.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

#include "globals.hpp"
#include "parameter_input.hpp"
#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "particles/relativistic_pusher.hpp"
#include "pgen.hpp"

#include <Kokkos_Random.hpp>

namespace {

struct PartRandomBField {
  Real x, y, z;
};

enum PartRandomBProfile {
  kUniformField=0,
  kLinearCrossField=1,
  kSinusoidalDivFreeField=2,
  kMirrorField=3,
  kGradBField=4,
  kTurbulentField=5
};

constexpr Real kTwoPi = 6.283185307179586476925286766559;

KOKKOS_INLINE_FUNCTION
Real WrapUnit(Real x) {
  while (x < 0.0) {x += 1.0;}
  while (x >= 1.0) {x -= 1.0;}
  return x;
}

KOKKOS_INLINE_FUNCTION
std::uint64_t HashParticleState(std::uint64_t value) {
  value += 0x9E3779B97F4A7C15ULL;
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
  return value ^ (value >> 31);
}

KOKKOS_INLINE_FUNCTION
Real TagRandomUnit(std::uint64_t seed, int species, int tag, int stream) {
  std::uint64_t value = seed ^ (static_cast<std::uint64_t>(species) << 32);
  value ^= static_cast<std::uint64_t>(static_cast<unsigned int>(tag));
  value ^= static_cast<std::uint64_t>(stream + 1) * 0xD2B74407B1CE6E93ULL;
  return static_cast<Real>(HashParticleState(value) >> 11)/
         static_cast<Real>(9007199254740992.0);
}

KOKKOS_INLINE_FUNCTION
Real CenterUnit(Real xmin, Real xmax, Real mesh_min, Real mesh_max) {
  return (0.5*(xmin + xmax) - mesh_min)/(mesh_max - mesh_min);
}

KOKKOS_INLINE_FUNCTION
Real HalfWidthUnit(Real xmin, Real xmax, Real mesh_min, Real mesh_max) {
  return 0.5*(xmax - xmin)/(mesh_max - mesh_min);
}

KOKKOS_INLINE_FUNCTION
Real PeriodicDistance(Real a, Real b) {
  Real d = fabs(a - b);
  return fmin(d, 1.0 - d);
}

KOKKOS_INLINE_FUNCTION
bool OverlapsMovingBox(const RegionSize &sz, const RegionSize &mesh_size,
                       Real c1, Real c2, Real c3, bool multi_d, bool three_d) {
  bool overlaps = PeriodicDistance(
      CenterUnit(sz.x1min, sz.x1max, mesh_size.x1min, mesh_size.x1max), c1)
      <= HalfWidthUnit(sz.x1min, sz.x1max, mesh_size.x1min, mesh_size.x1max) + 0.17;
  if (multi_d) {
    overlaps = overlaps && PeriodicDistance(
        CenterUnit(sz.x2min, sz.x2max, mesh_size.x2min, mesh_size.x2max), c2)
        <= HalfWidthUnit(sz.x2min, sz.x2max, mesh_size.x2min, mesh_size.x2max) + 0.13;
  }
  if (three_d) {
    overlaps = overlaps && PeriodicDistance(
        CenterUnit(sz.x3min, sz.x3max, mesh_size.x3min, mesh_size.x3max), c3)
        <= HalfWidthUnit(sz.x3min, sz.x3max, mesh_size.x3min, mesh_size.x3max) + 0.11;
  }
  return overlaps;
}

KOKKOS_INLINE_FUNCTION
PartRandomBField PartRandomLinearCrossField(Real x1, Real x2, Real b0x, Real b0y,
                                            Real b0z, Real bgrad) {
  return PartRandomBField{b0x + bgrad*x2, b0y + bgrad*x1, b0z};
}

KOKKOS_INLINE_FUNCTION
PartRandomBField PartRandomProfileField(int profile, Real x1, Real x2, Real x3,
                                        Real b0x, Real b0y, Real b0z,
                                        Real bgrad, Real bamp, Real kwave) {
  PartRandomBField bf{b0x, b0y, b0z};
  if (profile == kLinearCrossField) {
    bf = PartRandomLinearCrossField(x1, x2, b0x, b0y, b0z, bgrad);
  } else if (profile == kSinusoidalDivFreeField) {
    Real k = kTwoPi*kwave;
    bf.x += bamp*k*std::sin(k*x1)*(std::cos(k*x2) - std::cos(k*x3));
    bf.y += bamp*k*std::sin(k*x2)*(std::cos(k*x3) - std::cos(k*x1));
    bf.z += bamp*k*std::sin(k*x3)*(std::cos(k*x1) - std::cos(k*x2));
  } else if (profile == kMirrorField) {
    bf.x += -bgrad*x1*x3;
    bf.y += -bgrad*x2*x3;
    bf.z += bgrad*x3*x3;
  } else if (profile == kGradBField) {
    bf.z += bgrad*x1;
  } else if (profile == kTurbulentField) {
    Real k1 = kTwoPi*kwave;
    Real k2 = 2.0*k1;
    bf.x += bamp*k1*std::sin(k1*x1)*(std::cos(k1*x2) - std::cos(k1*x3));
    bf.y += bamp*k1*std::sin(k1*x2)*(std::cos(k1*x3) - std::cos(k1*x1));
    bf.z += bamp*k1*std::sin(k1*x3)*(std::cos(k1*x1) - std::cos(k1*x2));
    bf.x += 0.35*bamp*k2*std::sin(k2*(x1 + 0.13))*
            (std::cos(k2*(x2 - 0.07)) - std::cos(k2*(x3 + 0.11)));
    bf.y += 0.35*bamp*k2*std::sin(k2*(x2 - 0.07))*
            (std::cos(k2*(x3 + 0.11)) - std::cos(k2*(x1 + 0.13)));
    bf.z += 0.35*bamp*k2*std::sin(k2*(x3 + 0.11))*
            (std::cos(k2*(x1 + 0.13)) - std::cos(k2*(x2 - 0.07)));
  }
  return bf;
}

void PartRandomRefinementCondition(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &mblev = pmbp->pmb->mb_lev;
  auto &mb_size = pmbp->pmb->mb_size;
  const int nmb = pmbp->nmb_thispack;
  const int mbs = pmesh->gids_eachrank[global_variable::my_rank];
  const int root_level = pmesh->root_level;
  const int target_level = pmesh->max_level;
  const bool multi_d = pmesh->multi_d;
  const bool three_d = pmesh->three_d;
  const RegionSize mesh_size = pmesh->mesh_size;
  const Real phase = static_cast<Real>(pmesh->ncycle);

  par_for("part_random_refinement", DevExeSpace(), 0, nmb-1,
  KOKKOS_LAMBDA(int m) {
    bool refine_region = OverlapsMovingBox(
        mb_size.d_view(m), mesh_size, WrapUnit(0.21 + 0.19*phase),
        WrapUnit(0.37 + 0.11*phase), WrapUnit(0.53 + 0.07*phase),
        multi_d, three_d);
    refine_region = refine_region || OverlapsMovingBox(
        mb_size.d_view(m), mesh_size, WrapUnit(0.74 - 0.13*phase),
        WrapUnit(0.63 + 0.17*phase), WrapUnit(0.28 - 0.09*phase),
        multi_d, three_d);
    const int level = mblev.d_view(m);
    if (refine_region && level < target_level) {
      refine_flag.d_view(m + mbs) = 1;
    } else if ((!refine_region) && level > root_level) {
      refine_flag.d_view(m + mbs) = -1;
    }
  });

  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::PartRandom()
//! \brief Problem Generator for random particle positions/velocities

void ProblemGenerator::PartRandom(ParameterInput *pin, const bool restart) {
  std::string particle_refinement =
      pin->GetOrAddString("problem","particle_refinement","none");
  if (particle_refinement.compare("moving_boxes") == 0) {
    user_ref_func = PartRandomRefinementCondition;
  } else if (particle_refinement.compare("none") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown particle_refinement = '"
              << particle_refinement << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Random particles test requires <particles> block in input file"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  int prtcl_rst_flag = pmbp->ppart->prtcl_rst_flag;
  if (restart && !prtcl_rst_flag) return;

  // capture variables for the kernel
  auto &mbsize = pmbp->pmb->mb_size;
  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  auto &npart = pmbp->ppart->nprtcl_thispack;
  auto &npart_spec = pmbp->ppart->nprtcl_perspec_thispack;
  auto &nspecies = pmbp->ppart->nspecies;
  auto gids = pmbp->gids;
  auto gide = pmbp->gide;
  auto nmb_thispack = pmbp->nmb_thispack;

  if (prtcl_rst_flag) {
    std::string prst_fname = pin->GetString("problem","prtcl_res_file");
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", global_variable::my_rank);
    std::size_t rank_pos = prst_fname.find("rank_00000000");
    if (rank_pos != std::string::npos) {
      prst_fname.replace(rank_pos, sizeof("rank_00000000") - 1, rank_dir);
    }

    FILE* pfile = std::fopen(prst_fname.c_str(),"rb");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' could not be opened" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Real gen_data[3];
    if (std::fread(gen_data, sizeof(Real), 3, pfile) != 3) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' has an incomplete header" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    pmbp->ppart->dtnew = gen_data[1];
    pmy_mesh_->time = gen_data[0];
    HostArray1D<Real> host_part_data("host_part_data", 17*npart);
    if (std::fread(host_part_data.data(), sizeof(Real), 17*npart, pfile) !=
        static_cast<std::size_t>(17*npart)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle restart file '" << prst_fname
                << "' has incomplete particle data" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fclose(pfile);

    DvceArray1D<Real> part_data("part_data", 17*npart);
    Kokkos::deep_copy(part_data, host_part_data);
    par_for("part_restart",DevExeSpace(),0,(npart-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PGID,p) = static_cast<int>(part_data(17*p+0));
      pi(PTAG,p) = static_cast<int>(part_data(17*p+1));
      pi(PSP,p) = static_cast<int>(part_data(17*p+2));
      pr(IPX,p) = part_data(17*p+3);
      pr(IPY,p) = part_data(17*p+4);
      pr(IPZ,p) = part_data(17*p+5);
      pr(IPVX,p) = part_data(17*p+6);
      pr(IPVY,p) = part_data(17*p+7);
      pr(IPVZ,p) = part_data(17*p+8);
      pr(IPM,p) = part_data(17*p+9);
      pr(IPBX,p) = part_data(17*p+10);
      pr(IPBY,p) = part_data(17*p+11);
      pr(IPBZ,p) = part_data(17*p+12);
      pr(IPDX,p) = part_data(17*p+13);
      pr(IPDY,p) = part_data(17*p+14);
      pr(IPDZ,p) = part_data(17*p+15);
      pr(IPDB,p) = part_data(17*p+16);
    });
    if (restart) {
      pmbp->ppart->CheckConsistency("particle restart load");
      return;
    }
  } else {
    Real min_mass = pin->GetOrAddReal("particles","min_mass",1.0);
    Real mass_log_spacing = pin->GetOrAddReal("particles","mass_log_spacing",1.0);
    bool has_B0x = pin->DoesParameterExist("problem","B0x");
    bool has_B0y = pin->DoesParameterExist("problem","B0y");
    bool has_B0z = pin->DoesParameterExist("problem","B0z");
    Real B0x = pin->GetOrAddReal("problem","B0x",0.0);
    Real B0y = pin->GetOrAddReal("problem","B0y",0.0);
    Real B0z = pin->GetOrAddReal("problem","B0z",1.0);
    bool relativistic_hc =
        (pmbp->ppart->pusher == ParticlesPusher::relativistic_hc);
    bool initial_state_is_w = false;
    Real w0x = 0.0, w0y = 0.0, w0z = 0.0;
    Real cE0x = 0.0, cE0y = 0.0, cE0z = 0.0;
    Real c_model = pmbp->ppart->c_model;
    Real alpha_s = pmbp->ppart->alpha_s;
    if (relativistic_hc) {
      if (!has_B0x || !has_B0y || !has_B0z ||
          !(pin->DoesParameterExist("problem","relativistic_initial_state")) ||
          !(pin->DoesParameterExist("problem","cE0x")) ||
          !(pin->DoesParameterExist("problem","cE0y")) ||
          !(pin->DoesParameterExist("problem","cE0z"))) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' prescribed-test "
                  << "initialization requires explicit B0x, B0y, B0z, "
                  << "relativistic_initial_state, and cE0x, cE0y, cE0z" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string initial_state =
          pin->GetString("problem","relativistic_initial_state");
      initial_state_is_w = (initial_state.compare("w") == 0);
      if (!initial_state_is_w && initial_state.compare("velocity") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Unknown relativistic_initial_state = '"
                  << initial_state << "'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      bool has_w = pin->DoesParameterExist("problem","w0x") ||
                   pin->DoesParameterExist("problem","w0y") ||
                   pin->DoesParameterExist("problem","w0z");
      if (initial_state_is_w) {
        if (!(pin->DoesParameterExist("problem","w0x")) ||
            !(pin->DoesParameterExist("problem","w0y")) ||
            !(pin->DoesParameterExist("problem","w0z")) ||
            pin->DoesParameterExist("problem","v0x") ||
            pin->DoesParameterExist("problem","v0y") ||
            pin->DoesParameterExist("problem","v0z") ||
            pin->DoesParameterExist("problem","v0")) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Particle pusher 'relativistic_hc' w "
                    << "initialization requires exactly w0x, w0y, w0z and no "
                    << "velocity initialization values" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        w0x = pin->GetReal("problem","w0x");
        w0y = pin->GetReal("problem","w0y");
        w0z = pin->GetReal("problem","w0z");
      } else if (has_w) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' velocity "
                  << "initialization rejects w0x, w0y, w0z" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cE0x = pin->GetReal("problem","cE0x");
      cE0y = pin->GetReal("problem","cE0y");
      cE0z = pin->GetReal("problem","cE0z");
      if (!std::isfinite(B0x) || !std::isfinite(B0y) || !std::isfinite(B0z) ||
          !std::isfinite(w0x) || !std::isfinite(w0y) || !std::isfinite(w0z) ||
          !std::isfinite(cE0x) || !std::isfinite(cE0y) || !std::isfinite(cE0z)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' prescribed-test "
                  << "fields and initial state must be finite" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    std::string particle_position =
        pin->GetOrAddString("problem","particle_position","random");
    std::string particle_velocity =
        pin->GetOrAddString("problem","particle_velocity","random");
    bool has_v0x = pin->DoesParameterExist("problem","v0x");
    bool has_v0y = pin->DoesParameterExist("problem","v0y");
    bool has_v0z = pin->DoesParameterExist("problem","v0z");
    bool has_v0 = pin->DoesParameterExist("problem","v0");
    if (relativistic_hc && !initial_state_is_w &&
        (particle_velocity.compare("uniform") != 0 ||
         !has_v0x || !has_v0y || !has_v0z || has_v0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher 'relativistic_hc' velocity "
                << "initialization requires particle_velocity = uniform, explicit "
                << "v0x, v0y, v0z, and no scalar v0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (particle_position.compare("random") != 0 &&
        particle_position.compare("tag_random") != 0 &&
        particle_position.compare("center") != 0 &&
        particle_position.compare("meshblock_center") != 0 &&
        particle_position.compare("fixed") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Unknown particle_position = '"
                << particle_position << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (particle_velocity.compare("random") != 0 &&
        particle_velocity.compare("uniform") != 0 &&
        particle_velocity.compare("isotropic") != 0 &&
        particle_velocity.compare("isotropic_tag_random") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Unknown particle_velocity = '"
                << particle_velocity << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    bool center_particles = (particle_position.compare("center") == 0);
    bool meshblock_center_particles =
        (particle_position.compare("meshblock_center") == 0);
    bool block_center_particles = center_particles || meshblock_center_particles;
    bool fixed_particles = (particle_position.compare("fixed") == 0);
    bool tag_random_particles = (particle_position.compare("tag_random") == 0);
    Real fixed_x1 = pin->GetOrAddReal("problem","particle_x1",0.0);
    Real fixed_x2 = pin->GetOrAddReal("problem","particle_x2",0.0);
    Real fixed_x3 = pin->GetOrAddReal("problem","particle_x3",0.0);
    if (fixed_particles) {
      const RegionSize &mesh_size = pmy_mesh_->mesh_size;
      if (fixed_x1 < mesh_size.x1min || fixed_x1 >= mesh_size.x1max ||
          fixed_x2 < mesh_size.x2min || fixed_x2 >= mesh_size.x2max ||
          fixed_x3 < mesh_size.x3min || fixed_x3 >= mesh_size.x3max) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Fixed particle position is outside the mesh"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    bool uniform_velocity = (particle_velocity.compare("uniform") == 0);
    bool isotropic_velocity = (particle_velocity.compare("isotropic") == 0);
    bool tag_isotropic_velocity =
        (particle_velocity.compare("isotropic_tag_random") == 0);
    Real v0x = pin->GetOrAddReal("problem","v0x",1.0);
    Real v0y = pin->GetOrAddReal("problem","v0y",0.0);
    Real v0z = pin->GetOrAddReal("problem","v0z",0.0);
    Real v0 = pin->GetOrAddReal("problem","v0",1.0);
    std::uint64_t particle_seed = static_cast<std::uint64_t>(
        pin->GetOrAddInteger("problem","particle_seed",0));
    RegionSize mesh_size = pmy_mesh_->mesh_size;

    // initialize particles
    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
    DvceArray1D<int> relativistic_status("part_relativistic_init_status", 1);
    Kokkos::deep_copy(relativistic_status, 0);
    par_for("part_update",DevExeSpace(),0,(npart-1),
    KOKKOS_LAMBDA(const int p) {
      auto rand_gen = rand_pool64.get_state();  // get random number state this thread
      int spec = (npart_spec > 0) ? p/npart_spec : 0;
      spec = (spec < 0) ? 0 : spec;
      spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
      int p_in_spec = (npart_spec > 0) ? p - spec*npart_spec : p;
      int tag = pi(PTAG,p);

      // choose parent MeshBlock randomly unless a deterministic layout is requested
      int m = static_cast<int>(rand_gen.frand()*(gide - gids + 1.0));
      if ((meshblock_center_particles || tag_random_particles) && nmb_thispack > 0) {
        m = p_in_spec % nmb_thispack;
      }
      m = (m < 0) ? 0 : m;
      m = (m > nmb_thispack - 1) ? nmb_thispack - 1 : m;
      pi(PGID,p) = gids + m;
      pi(PSP,p) = spec;

      Real rand = rand_gen.frand();
      pr(IPX,p) = tag_random_particles ?
                  mesh_size.x1min + TagRandomUnit(particle_seed, spec, tag, 0)*
                  (mesh_size.x1max - mesh_size.x1min) :
                  fixed_particles ? fixed_x1 :
                  block_center_particles ?
                  0.5*(mbsize.d_view(m).x1min + mbsize.d_view(m).x1max) :
                  (1. - rand)*mbsize.d_view(m).x1min + rand*mbsize.d_view(m).x1max;
      if (!tag_random_particles) {
        pr(IPX,p) = fmin(pr(IPX,p),mbsize.d_view(m).x1max);
        pr(IPX,p) = fmax(pr(IPX,p),mbsize.d_view(m).x1min);
      }

      rand = rand_gen.frand();
      pr(IPY,p) = tag_random_particles ?
                  mesh_size.x2min + TagRandomUnit(particle_seed, spec, tag, 1)*
                  (mesh_size.x2max - mesh_size.x2min) :
                  fixed_particles ? fixed_x2 :
                  block_center_particles ?
                  0.5*(mbsize.d_view(m).x2min + mbsize.d_view(m).x2max) :
                  (1. - rand)*mbsize.d_view(m).x2min + rand*mbsize.d_view(m).x2max;
      if (!tag_random_particles) {
        pr(IPY,p) = fmin(pr(IPY,p),mbsize.d_view(m).x2max);
        pr(IPY,p) = fmax(pr(IPY,p),mbsize.d_view(m).x2min);
      }

      rand = rand_gen.frand();
      pr(IPZ,p) = tag_random_particles ?
                  mesh_size.x3min + TagRandomUnit(particle_seed, spec, tag, 2)*
                  (mesh_size.x3max - mesh_size.x3min) :
                  fixed_particles ? fixed_x3 :
                  block_center_particles ?
                  0.5*(mbsize.d_view(m).x3min + mbsize.d_view(m).x3max) :
                  (1. - rand)*mbsize.d_view(m).x3min + rand*mbsize.d_view(m).x3max;
      if (!tag_random_particles) {
        pr(IPZ,p) = fmin(pr(IPZ,p),mbsize.d_view(m).x3max);
        pr(IPZ,p) = fmax(pr(IPZ,p),mbsize.d_view(m).x3min);
      }

      if (uniform_velocity) {
        pr(IPVX,p) = v0x;
        pr(IPVY,p) = v0y;
        pr(IPVZ,p) = v0z;
      } else if (isotropic_velocity || tag_isotropic_velocity) {
        Real mu = tag_isotropic_velocity ?
                  2.0*TagRandomUnit(particle_seed, spec, tag, 3) - 1.0 :
                  2.0*rand_gen.frand() - 1.0;
        Real phi = tag_isotropic_velocity ?
                   kTwoPi*TagRandomUnit(particle_seed, spec, tag, 4) :
                   kTwoPi*rand_gen.frand();
        Real sintheta = sqrt(fmax(0.0, 1.0 - mu*mu));
        pr(IPVX,p) = v0*sintheta*cos(phi);
        pr(IPVY,p) = v0*sintheta*sin(phi);
        pr(IPVZ,p) = v0*mu;
      } else {
        pr(IPVX,p) = 2.0*(rand_gen.frand() - 0.5);
        pr(IPVY,p) = 2.0*(rand_gen.frand() - 0.5);
        pr(IPVZ,p) = 2.0*(rand_gen.frand() - 0.5);
      }
      pr(IPM,p) = min_mass*pow(mass_log_spacing, spec);
      pr(IPBX,p) = B0x;
      pr(IPBY,p) = B0y;
      pr(IPBZ,p) = B0z;
      pr(IPDX,p) = 0.0;
      pr(IPDY,p) = 0.0;
      pr(IPDZ,p) = 0.0;
      pr(IPDB,p) = 0.0;
      if (relativistic_hc) {
        particles::relativistic::Vector3 initial_w{w0x, w0y, w0z};
        particles::relativistic::StateStatus status =
            particles::relativistic::StateStatus::success;
        if (!initial_state_is_w) {
          Real gamma = 0.0;
          status = particles::relativistic::WFromVelocity(
              {pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p)}, c_model, initial_w, gamma);
        } else {
          status = particles::relativistic::ValidateWState(initial_w, c_model);
        }
        pr(IPCEX,p) = cE0x;
        pr(IPCEY,p) = cE0y;
        pr(IPCEZ,p) = cE0z;
        pr(IPWORK,p) = 0.0;
        pr(IPALPHA,p) = alpha_s;
        if (status == particles::relativistic::StateStatus::success) {
          status =
              particles::relativistic::StoreAuthoritativeWAndSynchronizeVelocityShadow(
                  pr, p, initial_w, c_model);
        }
        if (status != particles::relativistic::StateStatus::success) {
          Kokkos::atomic_max(&relativistic_status(0), static_cast<int>(status));
        }
      }

      rand_pool64.free_state(rand_gen);  // free state for use by other threads
    });
    auto relativistic_status_h =
        Kokkos::create_mirror_view_and_copy(HostMemSpace(), relativistic_status);
    if (relativistic_status_h(0) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher 'relativistic_hc' could not "
                << "initialize a finite representable authoritative w state"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  Real B0x = pin->GetOrAddReal("problem","B0x",0.0);
  Real B0y = pin->GetOrAddReal("problem","B0y",0.0);
  Real B0z = pin->GetOrAddReal("problem","B0z",1.0);
  std::string b_profile = pin->GetOrAddString("problem","B_profile","uniform");
  bool gather_diagnostic_only =
      pin->GetOrAddBoolean("problem","relativistic_gather_diagnostic_only",false);
  if (pmbp->ppart != nullptr &&
      pmbp->ppart->pusher == ParticlesPusher::relativistic_hc &&
      !gather_diagnostic_only &&
      b_profile.compare("uniform") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle pusher 'relativistic_hc' prescribed-test "
              << "execution requires <problem>/B_profile = uniform" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  int b_profile_id = kUniformField;
  if (b_profile.compare("uniform") == 0) {
    b_profile_id = kUniformField;
  } else if (b_profile.compare("linear_cross") == 0) {
    b_profile_id = kLinearCrossField;
  } else if (b_profile.compare("sinusoidal_divb_free") == 0) {
    b_profile_id = kSinusoidalDivFreeField;
  } else if (b_profile.compare("mirror") == 0) {
    b_profile_id = kMirrorField;
  } else if (b_profile.compare("gradb") == 0) {
    b_profile_id = kGradBField;
  } else if (b_profile.compare("turbulent") == 0) {
    b_profile_id = kTurbulentField;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown B_profile = '" << b_profile << "'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real Bgrad = pin->GetOrAddReal("problem","Bgrad",1.0);
  Real Bamp = pin->GetOrAddReal("problem","Bamp",0.05);
  Real Bwave = pin->GetOrAddReal("problem","Bwave_number",1.0);
  Real cfl_part = pin->GetOrAddReal("particles","cfl_part",0.05);

  if (pmbp->pmhd != nullptr) {
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = 1.0/eos.gamma;
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    auto &indcs = pmy_mesh_->mb_indcs;
    int &is = indcs.is, &ie = indcs.ie;
    int &js = indcs.js, &je = indcs.je;
    int &ks = indcs.ks, &ke = indcs.ke;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    auto &mb_size = pmbp->pmb->mb_size;
    bool is_ideal = eos.is_ideal;
    int b_profile_id_ = b_profile_id;
    par_for("pgen_mhd", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real x1c = CellCenterX(i - is, nx1, mb_size.d_view(m).x1min,
                             mb_size.d_view(m).x1max);
      Real x2c = CellCenterX(j - js, nx2, mb_size.d_view(m).x2min,
                             mb_size.d_view(m).x2max);
      Real x3c = CellCenterX(k - ks, nx3, mb_size.d_view(m).x3min,
                             mb_size.d_view(m).x3max);
      Real x1f = LeftEdgeX(i - is, nx1, mb_size.d_view(m).x1min,
                           mb_size.d_view(m).x1max);
      Real x2f = LeftEdgeX(j - js, nx2, mb_size.d_view(m).x2min,
                           mb_size.d_view(m).x2max);
      Real x3f = LeftEdgeX(k - ks, nx3, mb_size.d_view(m).x3min,
                           mb_size.d_view(m).x3max);
      PartRandomBField bf = PartRandomProfileField(
          b_profile_id_, x1c, x2c, x3c, B0x, B0y, B0z, Bgrad, Bamp, Bwave);
      u0(m,IDN,k,j,i) = 1.0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      bcc0(m,IBX,k,j,i) = bf.x;
      bcc0(m,IBY,k,j,i) = bf.y;
      bcc0(m,IBZ,k,j,i) = bf.z;
      b0.x1f(m,k,j,i) = PartRandomProfileField(
          b_profile_id_, x1f, x2c, x3c, B0x, B0y, B0z, Bgrad, Bamp, Bwave).x;
      b0.x2f(m,k,j,i) = PartRandomProfileField(
          b_profile_id_, x1c, x2f, x3c, B0x, B0y, B0z, Bgrad, Bamp, Bwave).y;
      b0.x3f(m,k,j,i) = PartRandomProfileField(
          b_profile_id_, x1c, x2c, x3f, B0x, B0y, B0z, Bgrad, Bamp, Bwave).z;
      if (i==ie) {
        Real x1fr = LeftEdgeX(i + 1 - is, nx1, mb_size.d_view(m).x1min,
                              mb_size.d_view(m).x1max);
        b0.x1f(m,k,j,i+1) = PartRandomProfileField(
            b_profile_id_, x1fr, x2c, x3c, B0x, B0y, B0z, Bgrad, Bamp, Bwave).x;
      }
      if (j==je) {
        Real x2fr = LeftEdgeX(j + 1 - js, nx2, mb_size.d_view(m).x2min,
                              mb_size.d_view(m).x2max);
        b0.x2f(m,k,j+1,i) = PartRandomProfileField(
            b_profile_id_, x1c, x2fr, x3c, B0x, B0y, B0z, Bgrad, Bamp, Bwave).y;
      }
      if (k==ke) {
        Real x3fr = LeftEdgeX(k + 1 - ks, nx3, mb_size.d_view(m).x3min,
                              mb_size.d_view(m).x3max);
        b0.x3f(m,k+1,j,i) = PartRandomProfileField(
            b_profile_id_, x1c, x2c, x3fr, B0x, B0y, B0z, Bgrad, Bamp, Bwave).z;
      }
      if (is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*(bf.x*bf.x + bf.y*bf.y + bf.z*bf.z);
      }
    });
  }

  // set timestep (which will remain constant for entire run
  // Assumes uniform mesh (no SMR or AMR)
  // Assumes velocities normalized to one, so dt=min(dx)
  Real &dtnew_ = pmbp->ppart->dtnew;
  dtnew_ = std::min(mbsize.h_view(0).dx1, mbsize.h_view(0).dx2);
  dtnew_ = std::min(dtnew_, mbsize.h_view(0).dx3);
  Real min_mass = pin->GetOrAddReal("particles","min_mass",1.0);
  Real bmag = std::sqrt(B0x*B0x + B0y*B0y + B0z*B0z);
  if (bmag > 0.0) {
    dtnew_ = std::min(dtnew_, cfl_part*min_mass/bmag);
  }
  std::string particle_position =
      pin->GetOrAddString("problem","particle_position","random");
  if ((!restart) && particle_position.compare("tag_random") == 0) {
    pmbp->ppart->RemapAfterAMR();
  }

  return;
}

#if PART_RANDOM_USER_PROBLEM_ENABLED
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  PartRandom(pin, restart);
}
#endif
