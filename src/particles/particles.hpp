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
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations

// constants that enumerate ParticlesPusher options
enum class ParticlesPusher {drift, leap_frog, lagrangian_tracer, lagrangian_mc,
                            boris, relativistic_hc};

// constants that enumerate particle field-gather interpolation options
enum class ParticleInterpolation {lin_legacy, trilinear, tsc};

// constants that enumerate relativistic CR field-source options
enum class RelativisticFieldSource {prescribed_test};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray};

// constants that enumerate particle consistency check levels
enum class ParticlesConsistencyMode {none, counts, local, full};

// constants that enumerate AMR remap implementations
enum class ParticlesAMRRemapMode {device_table, host_tree};

// constants that enumerate particle MPI metadata exchange implementations
enum class ParticlesExchangeMode {alltoall_counts, allgather};

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
  int nprtcl_perspec_thispack;     // number of particles per species in this pack
  int nrdata, nidata;
  int nspecies;
//  DvceArray1D<int>  prtcl_gid;     // GID of MeshBlock containing each par
//  DvceArray2D<Real> prtcl_pos;     // positions
//  DvceArray2D<Real> prtcl_vel;     // velocities
  DvceArray2D<Real> prtcl_rdata;   // real number properties each particle (x,v,etc.)
  DvceArray2D<int>  prtcl_idata;   // integer properties each particle (gid, tag, etc.)
  Real dtnew;
  int is_dynamic;
  int prtcl_rst_flag;
  bool check_consistency;
  bool check_motion_bounds;
  bool log_performance;
  bool validate_amr_lookup;
  bool subcycle;
  bool subcycle_strict;
  int amr_lookup_max_cells;
  int subcycle_max_steps;
  Real subcycle_cell_fraction;
  Real subcycle_meshblock_fraction;
  Real subcycle_gyro_fraction;
  Real c_model;
  Real alpha_s;
  ParticlesConsistencyMode consistency_mode;
  ParticlesAMRRemapMode amr_remap_mode;
  ParticlesExchangeMode exchange_mode;
  bool reference_counts_set;
  int reference_nprtcl_total;
  std::vector<int> reference_nprtcl_eachspec;

  ParticlesPusher pusher;
  ParticleInterpolation interpolation;
  RelativisticFieldSource relativistic_field_source;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;

  // container to hold names of TaskIDs
  ParticlesTaskIDs id;

  // functions...
  void CreateParticleTags(ParameterInput *pin);
  void SetConsistencyReference();
  void CheckConsistency(const std::string &label, bool check_rank=true);
  void CheckMotionBounds(const std::string &label);
  void LogPerformance(const std::string &label, int64_t npushed, int64_t nsent,
                      int64_t nrecv, int64_t nremapped, int64_t ndiag,
                      int64_t nmessages=0, int64_t nbytes=0);
  int ComputeSubcycleSteps(Real dt, int &cell_steps, int &block_steps,
                           int &gyro_steps);
  void LogSubcycle(int nsub, int cell_steps, int block_steps, int gyro_steps);
  void ExchangeAfterSubcycle();
  void RemapAfterAMR();
  void AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus Push(Driver *pdriver, int stage);
  TaskStatus NewGID(Driver *pdriver, int stage);
  TaskStatus SendCnt(Driver *pdriver, int stage);
  TaskStatus InitRecv(Driver *pdriver, int stage);
  TaskStatus SendP(Driver *pdriver, int stage);
  TaskStatus RecvP(Driver *pdriver, int stage);
  TaskStatus ClearSend(Driver *pdriver, int stage);
  TaskStatus ClearRecv(Driver *pdriver, int stage);

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
