#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.hpp
//  \brief definitions for Particles class

#include <algorithm>
#include <map>
#include <memory>
#include <cstdint>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations
class IOWrapper;

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {drift, leap_frog, lagrangian_tracer, lagrangian_mc};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray, lagrangian_mc};

enum class TracerSeedWeight {mass, volume};
enum class TracerSeedRegion {all, box, sphere, slab};
enum class TracerSeedTarget {none, density, temperature, pressure, entropy, scalar};

struct TracerSeedSchedule {
  int id = 0;
  std::string block_name;
  Real start_time = 0.0;
  Real end_time = 0.0;
  Real cadence = 0.0;
  Real next_time = 0.0;
  int event_index = 0;
  bool complete = false;
  int count_per_event = 0;
  int seed = 0;
  TracerSeedWeight weight = TracerSeedWeight::mass;
  TracerSeedRegion region = TracerSeedRegion::all;
  Real x1min = 0.0, x1max = 0.0;
  Real x2min = 0.0, x2max = 0.0;
  Real x3min = 0.0, x3max = 0.0;
  Real center1 = 0.0, center2 = 0.0, center3 = 0.0, radius = 0.0;
  int slab_axis = 1;
  Real slab_min = 0.0, slab_max = 0.0;
  TracerSeedTarget target = TracerSeedTarget::none;
  int scalar_index = -1;
  bool has_target_min = false;
  bool has_target_max = false;
  Real target_min = 0.0;
  Real target_max = 0.0;
};

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
  TaskID mradj;
  TaskID seed;
  TaskID sample;
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
//  DvceArray1D<int>  prtcl_gid;     // GID of MeshBlock containing each par
//  DvceArray2D<Real> prtcl_pos;     // positions
//  DvceArray2D<Real> prtcl_vel;     // velocities
  DvceArray2D<Real> prtcl_rdata;   // real number properties each particle (x,v,etc.)
  DvceArray2D<int>  prtcl_idata;   // integer properties each particle (gid, tag, etc.)
  Real dtnew;
  std::int64_t random_seed = 0;
  std::int64_t next_tracer_tag = 0;

  ParticlesPusher pusher;

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
  TaskStatus AdjustMeshRefinement(Driver *pdriver, int stage);
  TaskStatus SeedDueTracers(Driver *pdriver, int stage);
  TaskStatus SampleThermodynamicsTask(Driver *pdriver, int stage);
  TaskStatus PushLagrangianMC(Driver *pdriver, int stage);
  void SeedInitialTracers();
  void SampleThermodynamics();
  void RemapAfterMeshRefinement();
  void WriteRestartData(IOWrapper &resfile, bool single_file_per_rank);
  void ReadRestartData(IOWrapper &resfile, bool single_file_per_rank);
  int GetLagrangianMCScalarCount() const {return std::max(0, nrdata - LMC_SCALAR0);}
  bool IsLagrangianMC() const {return particle_type == ParticleType::lagrangian_mc;}

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles
  std::vector<TracerSeedSchedule> seed_schedules_;

  void ParseTracerSeedSchedules(ParameterInput *pin);
  void SeedTracersAtTime(Real event_time, bool initial_only);
  void AppendParticles(const HostArray2D<Real> &new_rdata,
                       const HostArray2D<int> &new_idata, int nnew);
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
