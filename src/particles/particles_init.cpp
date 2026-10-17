//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_init.cpp
//! \brief Seed or restore particles without modifying the loaded fluid snapshot.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles.hpp"

namespace {

void ParticleInputError(const std::string &message) {
  std::cout << "### FATAL ERROR: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

// The donor's tag-based generator makes seeding independent of execution order.
KOKKOS_INLINE_FUNCTION
Real TagRandomUnit(std::uint64_t seed, int species, int tag, int stream) {
  std::uint64_t value = seed ^ (static_cast<std::uint64_t>(species) << 32);
  value ^= static_cast<std::uint64_t>(static_cast<unsigned int>(tag));
  value ^= static_cast<std::uint64_t>(stream + 1) * 0xD2B74407B1CE6E93ULL;
  value += 0x9E3779B97F4A7C15ULL;
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
  value ^= value >> 31;
  return static_cast<Real>(value >> 11)/static_cast<Real>(9007199254740992.0);
}

} // namespace

namespace particles {

void Particles::InitializeFromSnapshot(ParameterInput *pin) {
  if (pmy_pack->pmesh->nprtcl_total == 0) {
    ParticleInputError("Frozen tracers require at least one particle globally");
  }
  auto &pr = prtcl_rdata;
  auto &pi = prtcl_idata;
  const int npart = nprtcl_thispack;
  if (prtcl_rst_flag) {
    std::string fname = pin->GetString("problem", "prtcl_res_file");
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", global_variable::my_rank);
    std::size_t pos = fname.find("rank_00000000");
    if (pos != std::string::npos) {fname.replace(pos, 13, rank_dir);}
    FILE *file = std::fopen(fname.c_str(), "rb");
    if (file == nullptr) {ParticleInputError("Cannot open particle restart " + fname);}
    Real header[3];
    if (std::fread(header, sizeof(Real), 3, file) != 3) {
      ParticleInputError("Incomplete particle restart header");
    }
    if (!std::isfinite(header[0]) || !std::isfinite(header[1]) || header[1] <= 0.0 ||
        header[2] != static_cast<Real>(npart)) {
      ParticleInputError("Invalid particle restart header");
    }
    if (header[0] != pmy_pack->pmesh->time) {
      ParticleInputError("Particle and fluid checkpoints must have the same time");
    }
    HostArray2D<Real> hr("particle_restart_real", nrdata, npart);
    HostArray2D<int> hi("particle_restart_int", nidata, npart);
    for (int p=0; p<npart; ++p) {
      Real record[17];
      if (std::fread(record, sizeof(Real), 17, file) != 17) {
        ParticleInputError("Incomplete particle restart data");
      }
      for (int n=0; n<17; ++n) {
        if (!std::isfinite(record[n])) {ParticleInputError("Nonfinite particle restart");}
      }
      for (int n=0; n<3; ++n) {
        if (record[n] < 0.0 || static_cast<double>(record[n]) >
            std::numeric_limits<int>::max() || std::trunc(record[n]) != record[n]) {
          ParticleInputError("Invalid integer in particle restart");
        }
        hi(n,p) = static_cast<int>(record[n]);
      }
      if (hi(PGID,p) >= pmy_pack->pmesh->nmb_total || hi(PSP,p) >= nspecies ||
          record[9] <= 0.0) {
        ParticleInputError("Invalid particle restart owner, species, or gyro parameter");
      }
      const int fields[14] = {IPX, IPY, IPZ, IPVX, IPVY, IPVZ, IPM,
                              IPBX, IPBY, IPBZ, IPDX, IPDY, IPDZ, IPDB};
      for (int n=0; n<14; ++n) {hr(fields[n],p) = record[n+3];}
    }
    if (std::fgetc(file) != EOF) {ParticleInputError("Unexpected particle restart data");}
    std::fclose(file);
    Kokkos::deep_copy(pr, hr);
    Kokkos::deep_copy(pi, hi);
  } else {
    if (!pin->GetOrAddBoolean("particles", "seed_on_restart", false)) {
      ParticleInputError("Frozen tracers require seed_on_restart=true or prtcl_rst_flag=1");
    }
    std::string position = pin->GetOrAddString("problem", "particle_position", "tag_random");
    std::string velocity = pin->GetOrAddString("problem", "particle_velocity",
                                              "isotropic_tag_random");
    if (position != "tag_random" && position != "fixed") {
      ParticleInputError("particle_position must be tag_random or fixed");
    }
    if (velocity != "isotropic_tag_random" && velocity != "uniform") {
      ParticleInputError("particle_velocity must be isotropic_tag_random or uniform");
    }
    Real mass = pin->GetOrAddReal("particles", "min_mass", 1.0);
    Real spacing = pin->GetOrAddReal("particles", "mass_log_spacing", 1.0);
    Real vx = pin->GetOrAddReal("problem", "v0x", 1.0);
    Real vy = pin->GetOrAddReal("problem", "v0y", 0.0);
    Real vz = pin->GetOrAddReal("problem", "v0z", 0.0);
    Real speed = std::sqrt(vx*vx + vy*vy + vz*vz);
    if (!std::isfinite(speed) || speed <= 0.0 || !std::isfinite(mass) || mass <= 0.0 ||
        !std::isfinite(spacing) || spacing <= 0.0) {
      ParticleInputError("Particle speed, min_mass, and mass_log_spacing must be positive");
    }
    for (int s=0; s<nspecies; ++s) {
      Real species_mass = mass*std::pow(spacing, s);
      if (!std::isfinite(species_mass) || species_mass <= 0.0) {
        ParticleInputError("Particle species gyro parameter overflow/underflow");
      }
    }
    Real x = pin->GetOrAddReal("problem", "particle_x1", 0.0);
    Real y = pin->GetOrAddReal("problem", "particle_x2", 0.0);
    Real z = pin->GetOrAddReal("problem", "particle_x3", 0.0);
    auto mesh = pmy_pack->pmesh->mesh_size;
    bool fixed = (position == "fixed");
    if (fixed && !(x >= mesh.x1min && x < mesh.x1max && y >= mesh.x2min &&
                   y < mesh.x2max && z >= mesh.x3min && z < mesh.x3max)) {
      ParticleInputError("Fixed particle position must lie inside the periodic domain");
    }
    bool isotropic = (velocity == "isotropic_tag_random");
    std::uint64_t seed = pin->GetOrAddInteger("problem", "particle_seed", 0);
    int nper = nprtcl_perspec_thispack;
    int gids = pmy_pack->gids;
    Kokkos::deep_copy(pr, 0.0);
    par_for("seed_snapshot_particles", DevExeSpace(), 0, npart-1,
    KOKKOS_LAMBDA(const int p) {
      int spec = p/nper;
      int tag = pi(PTAG,p);
      pi(PSP,p) = spec;
      pi(PGID,p) = gids;
      pr(IPX,p) = fixed ? x : mesh.x1min +
          (mesh.x1max - mesh.x1min)*TagRandomUnit(seed, spec, tag, 0);
      pr(IPY,p) = fixed ? y : mesh.x2min +
          (mesh.x2max - mesh.x2min)*TagRandomUnit(seed, spec, tag, 1);
      pr(IPZ,p) = fixed ? z : mesh.x3min +
          (mesh.x3max - mesh.x3min)*TagRandomUnit(seed, spec, tag, 2);
      if (isotropic) {
        Real mu = 2.0*TagRandomUnit(seed, spec, tag, 3) - 1.0;
        Real phi = 6.2831853071795864769*TagRandomUnit(seed, spec, tag, 4);
        Real sintheta = Kokkos::sqrt(Kokkos::fmax(0.0, 1.0 - mu*mu));
        pr(IPVX,p) = speed*sintheta*Kokkos::cos(phi);
        pr(IPVY,p) = speed*sintheta*Kokkos::sin(phi);
        pr(IPVZ,p) = speed*mu;
      } else {
        pr(IPVX,p) = vx;
        pr(IPVY,p) = vy;
        pr(IPVZ,p) = vz;
      }
      pr(IPM,p) = mass*Kokkos::pow(spacing, spec);
    });
  }
  if (check_consistency) {SetConsistencyReference();}
  RemapAfterAMR();
  CheckConsistency("frozen snapshot particle initialization");
}

} // namespace particles
