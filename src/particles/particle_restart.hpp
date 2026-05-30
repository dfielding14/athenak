#ifndef PARTICLES_PARTICLE_RESTART_HPP_
#define PARTICLES_PARTICLE_RESTART_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file particle_restart.hpp
//! \brief Shared header format for sidecar particle restart files.

#include "athena.hpp"

namespace particles {

static constexpr int kParticleRestartMagicSize = 16;
static constexpr int kParticleRestartBasenameSize = 128;
static constexpr int kParticleRestartVersion = 2;
static constexpr int kParticleRestartEndianMarker = 0x01020304;
static constexpr const char *kParticleRestartMagic = "ATHKPRTCLRST1";

struct ParticleRestartHeader {
  char magic[kParticleRestartMagicSize];
  int version;
  int header_size;
  int endian_marker;
  int real_size;
  int int_size;
  int ncycle;
  int total_particles;
  int nrdata;
  int nidata;
  int max_tag;
  int basename_length;
  Real time;
  Real star_mass_formed_total;
  Real star_mass_accreted_total;
  char basename[kParticleRestartBasenameSize];
};

} // namespace particles

#endif // PARTICLES_PARTICLE_RESTART_HPP_
