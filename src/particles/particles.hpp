#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.hpp
//  \brief definitions for Particles class

#include <map>
#include <memory>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {drift, rk4_gravity, leap_frog, lagrangian_tracer, lagrangian_mc,
                            boris_lin, boris_tsc};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray, star};

//----------------------------------------------------------------------------------------
//! \struct ParticlesTaskIDs
//  \brief container to hold TaskIDs of all particles tasks

struct ParticlesTaskIDs {
  TaskID push;
  TaskID newgid;
  TaskID count;
  TaskID irecv;
  TaskID sendp;
  TaskID recvp;
  TaskID csend;
  TaskID crecv;
};

namespace particles {

//----------------------------------------------------------------------------------------
//! \class Particles

class Particles {
  friend class ParticlesBoundaryValues;
 public:
  Particles(MeshBlockPack *ppack, ParameterInput *pin);
  ~Particles();

  // data
  ParticleType particle_type;
  int nprtcl_thispack;             // number of particles this MeshBlockPack
  int nrdata, nidata;
  DvceArray2D<Real> prtcl_rdata;   // real number properties each particle (x,v,etc.)
  DvceArray2D<int>  prtcl_idata;   // integer properties each particle (gid, tag, etc.)
  Real dtnew;

  ParticlesPusher pusher;
  
  // Cosmic ray specific
  int nspecies;                    // number of CR species
  bool track_displacement;         // enable displacement tracking
  DvceArray1D<Real> species_mass;  // mass per species
  DvceArray1D<Real> species_charge;// charge per species
 
  // Constants for rk4_gravity pusher
  Real r_scale;
  Real rho_scale;
  Real m_gal;
  Real a_gal;
  Real z_gal;
  Real r_200;
  Real rho_mean;
  Real par_grav_dx;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;

  // container to hold names of TaskIDs
  ParticlesTaskIDs id;

  // functions...
  void CreateParticleTags(ParameterInput *pin);
  void AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus Push(Driver *pdriver, int stage);
  TaskStatus NewGID(Driver *pdriver, int stage);
  TaskStatus SendCnt(Driver *pdriver, int stage);
  TaskStatus InitRecv(Driver *pdriver, int stage);
  TaskStatus SendP(Driver *pdriver, int stage);
  TaskStatus RecvP(Driver *pdriver, int stage);
  TaskStatus ClearSend(Driver *pdriver, int stage);
  TaskStatus ClearRecv(Driver *pdriver, int stage);
  
  // Cosmic ray specific methods
  void InitializeCosmicRays(ParameterInput *pin);
  TaskStatus PushCosmicRays(Driver *pdriver, int stage);
  TaskStatus PushStars(Driver *pdriver, int stage);
  
  // Field interpolation methods
  KOKKOS_INLINE_FUNCTION
  void InterpolateLinear(int m, Real x, Real y, Real z,
                        Real &Bx, Real &By, Real &Bz) const;
  KOKKOS_INLINE_FUNCTION
  void InterpolateTSC(int m, Real x, Real y, Real z,
                     Real &Bx, Real &By, Real &Bz) const;

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
