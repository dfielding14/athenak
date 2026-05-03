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
#include "bvals/bvals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"

// forward declarations

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {
  drift,
  rk4_gravity,
  leap_frog,
  lagrangian_tracer,
  lagrangian_mc,
  boris_lin,
  boris_tsc
};

// constants that enumerate ParticleTypes
enum class ParticleType { cosmic_ray, star, lagrangian_mc };

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
  TaskID inject; // Lagrangian MC inflow particle injection
  TaskID mradj;  // Mesh refinement adjustment for lagrangian_mc (AMR)
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
  int nprtcl_thispack; // number of particles this MeshBlockPack
  int nrdata, nidata;
  DvceArray2D<Real>
      prtcl_rdata; // real number properties each particle (x,v,etc.)
  DvceArray2D<int>
      prtcl_idata; // integer properties each particle (gid, tag, etc.)
  Real dtnew;

  ParticlesPusher pusher;

  // Lagrangian MC tracer particle specific
  int64_t random_seed;  // Base seed for deterministic RNG
  Real min_radius;      // Particles within this radius won't be updated (-1 = disabled)
  bool lmc_inject_at_inflow;  // Add new MC tracers through physical inflow faces
  int64_t lmc_inject_seed;    // Base seed for deterministic inflow injection sampling
  int lmc_max_inject_per_step;  // Optional per-rank cap; negative disables the cap
  int next_tag;          // Next globally unique tag for newly injected particles
  Real lmc_mass_per_particle;  // Mass represented by one MC tracer

  // Cosmic ray specific
  int nspecies;                     // number of CR species
  bool track_displacement;          // enable displacement tracking
  DvceArray1D<Real> species_mass;   // mass per species
  DvceArray1D<Real> species_charge; // charge per species

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
  void InitializeStars(ParameterInput *pin);
  TaskStatus PushDrift(Driver *pdriver, int stage);
  TaskStatus PushCosmicRays(Driver *pdriver, int stage);
  TaskStatus PushStars(Driver *pdriver, int stage);

  // Lagrangian MC tracer particle methods
  void ReallocateParticles(int new_count);
  TaskStatus PushLagrangianMC(Driver *pdriver, int stage);
  TaskStatus InjectLagrangianMCInflow(Driver *pdriver, int stage);
  TaskStatus AdjustMeshRefinement(Driver *pdriver, int stage);

  // Field interpolation methods
  KOKKOS_INLINE_FUNCTION
  void InterpolateLinear(int m, Real x, Real y, Real z, Real &Bx, Real &By,
                         Real &Bz) const;
  KOKKOS_INLINE_FUNCTION
  void InterpolateTSC(int m, Real x, Real y, Real z, Real &Bx, Real &By,
                      Real &Bz) const;

private:
  MeshBlockPack *pmy_pack; // ptr to MeshBlockPack containing this Particles
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
