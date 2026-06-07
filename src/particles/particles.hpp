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
#include "particles/tracer_fields.hpp"

// forward declarations
class IOWrapper;

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {drift, leap_frog, lagrangian_tracer, lagrangian_mc, ito2};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray, lagrangian_mc, lagrangian_ito};

enum ItoCoefficientIndex {
  ITO_M1=0, ITO_M2=1, ITO_M3=2,
  ITO_Q11=3, ITO_Q22=4, ITO_Q33=5,
  ITO_Q12=6, ITO_Q13=7, ITO_Q23=8,
  ITO_NCOEFF=9
};

enum class TracerSeedWeight {mass, volume};
enum class TracerSeedRegion {all, box, sphere, slab};

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
  bool has_target = false;
  particles::TracerField target_field;
  bool has_target_min = false;
  bool has_target_max = false;
  Real target_min = 0.0;
  Real target_max = 0.0;
};

//----------------------------------------------------------------------------------------
//! \struct ParticlesTaskIDs
//  \brief container to hold TaskIDs of all particles tasks

struct ParticlesTaskIDs {
  TaskID ito_build;
  TaskID ito_restrict;
  TaskID ito_irecv;
  TaskID ito_send;
  TaskID ito_recv;
  TaskID ito_crecv;
  TaskID ito_csend;
  TaskID ito_prolong;
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
  DvceArray2D<int>  prtcl_idata;   // small integer properties (gid, seed id, etc.)
  DvceArray1D<std::uint64_t> prtcl_tag;  // globally unique particle tags
  Real dtnew;
  Real ito_probability_target = 0.99;
  std::int64_t random_seed = 0;
  std::uint64_t next_tracer_tag = 0;
  DvceArray5D<Real> ito_coeff;
  DvceArray5D<Real> coarse_ito_coeff;
  DvceArray1D<int> ito_invalid;

  ParticlesPusher pusher;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;
  MeshBoundaryValuesCC *pbval_ito = nullptr;

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
  TaskStatus PushLagrangianMC(Driver *pdriver, int stage);
  TaskStatus BuildItoCoefficients(Driver *pdriver, int stage);
  TaskStatus RestrictItoCoefficients(Driver *pdriver, int stage);
  TaskStatus InitRecvItoCoefficients(Driver *pdriver, int stage);
  TaskStatus SendItoCoefficients(Driver *pdriver, int stage);
  TaskStatus RecvItoCoefficients(Driver *pdriver, int stage);
  TaskStatus ClearRecvItoCoefficients(Driver *pdriver, int stage);
  TaskStatus ClearSendItoCoefficients(Driver *pdriver, int stage);
  TaskStatus ProlongateItoCoefficients(Driver *pdriver, int stage);
  TaskStatus PushIto2(Driver *pdriver, int stage);
  void SeedInitialTracers();
  void RemapAfterMeshRefinement();
  void WriteRestartData(IOWrapper &resfile, bool single_file_per_rank);
  void ReadRestartData(IOWrapper &resfile, bool single_file_per_rank);
  int GetLagrangianMCScalarCount() const;
  bool IsLagrangianMC() const {return particle_type == ParticleType::lagrangian_mc;}
  bool IsIto2() const {return particle_type == ParticleType::lagrangian_ito;}
  bool IsFluxTracer() const {return IsLagrangianMC() || IsIto2();}

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles
  std::vector<TracerSeedSchedule> seed_schedules_;

  void CheckMassFloorCompatibility();
  void ParseTracerSeedSchedules(ParameterInput *pin);
  void SeedTracersAtTime(Real event_time, bool initial_only);
  void AppendParticles(const HostArray2D<Real> &new_rdata,
                       const HostArray2D<int> &new_idata,
                       const HostArray1D<std::uint64_t> &new_tags, int nnew);
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
