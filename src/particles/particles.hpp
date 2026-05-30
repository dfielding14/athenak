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
enum class ParticlesPusher {drift, gravity, leap_frog, lagrangian_tracer, lagrangian_mc};

// constants that enumerate ParticleTypes
enum class ParticleType {cosmic_ray, star};
enum class StarGravityForceMethod {direct, tree};
enum class StarGravityIntegrator {kdk, rk4};
enum class StarGravityExternalMode {none, constant, point_mass, user};
enum class StarGravityTimestepMode {fixed, acceleration, pair_orbit};
enum class StarGravityPeriodicMode {domain, minimum_image, error};

//----------------------------------------------------------------------------------------
//! \struct ParticlesTaskIDs
//  \brief container to hold TaskIDs of all particles tasks

struct ParticlesTaskIDs {
  TaskID form;
  TaskID accrete;
  TaskID push;
  TaskID newgid;
  TaskID count;
  TaskID irecv;
  TaskID sendp;
  TaskID recvp;
  TaskID csend;
  TaskID crecv;
  TaskID gravity_finish;
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

  ParticlesPusher pusher;

  // Star-particle controls
  bool star_init_from_file = false;
  bool star_init_from_restart = false;
  bool star_formation_enabled = false;
  bool star_accretion_enabled = false;
  bool star_remove_gas_on_formation = true;
  int star_formation_interval = 1;
  int star_formation_max_per_cycle = -1;
  int star_accretion_radius_cells = 0;
  int next_star_tag = 0;
  Real star_formation_density_threshold = 0.0;
  Real star_formation_particle_mass = 0.0;
  Real star_formation_density_floor = 0.0;
  Real star_accretion_rate = 0.0;
  Real star_accretion_max_fraction = 0.25;
  Real star_accretion_density_floor = 0.0;
  Real star_mass_formed_total = 0.0;
  Real star_mass_accreted_total = 0.0;

  // Star-particle gravity controls and diagnostics
  bool star_gravity_enabled = false;
  bool star_gravity_self_enabled = true;
  bool star_gravity_exact_diagnostics = true;
  bool star_gravity_conserve_accretion_momentum = true;
  StarGravityForceMethod star_gravity_force_method = StarGravityForceMethod::direct;
  StarGravityIntegrator star_gravity_integrator = StarGravityIntegrator::kdk;
  StarGravityExternalMode star_gravity_external_mode = StarGravityExternalMode::none;
  StarGravityTimestepMode star_gravity_timestep_mode = StarGravityTimestepMode::fixed;
  StarGravityPeriodicMode star_gravity_periodic_mode = StarGravityPeriodicMode::domain;
  Real star_gravity_constant = 0.0;
  Real star_gravity_softening = 0.0;
  Real star_gravity_timestep_eta = 0.02;
  Real star_gravity_max_timestep = -1.0;
  Real star_gravity_tree_theta = 0.7;
  Real star_gravity_base_timestep = 0.0;
  Real star_gravity_external_ax = 0.0;
  Real star_gravity_external_ay = 0.0;
  Real star_gravity_external_az = 0.0;
  Real star_gravity_external_mass = 0.0;
  Real star_gravity_external_x = 0.0;
  Real star_gravity_external_y = 0.0;
  Real star_gravity_external_z = 0.0;
  Real star_gravity_external_softening = 0.0;
  Real star_gravity_diag_dt = 0.0;
  Real star_gravity_diag_ekin = 0.0;
  Real star_gravity_diag_epot = 0.0;
  Real star_gravity_diag_px = 0.0;
  Real star_gravity_diag_py = 0.0;
  Real star_gravity_diag_pz = 0.0;
  Real star_gravity_diag_lx = 0.0;
  Real star_gravity_diag_ly = 0.0;
  Real star_gravity_diag_lz = 0.0;
  Real star_gravity_diag_com_x = 0.0;
  Real star_gravity_diag_com_y = 0.0;
  Real star_gravity_diag_com_z = 0.0;
  Real star_gravity_diag_amax = 0.0;
  Real star_gravity_diag_rmin = 0.0;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;

  // container to hold names of TaskIDs
  ParticlesTaskIDs id;

  // functions...
  void CreateParticleTags(ParameterInput *pin);
  void AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus FormStars(Driver *pdriver, int stage);
  TaskStatus AccreteStars(Driver *pdriver, int stage);
  TaskStatus Push(Driver *pdriver, int stage);
  TaskStatus FinishGravity(Driver *pdriver, int stage);
  TaskStatus NewGID(Driver *pdriver, int stage);
  TaskStatus SendCnt(Driver *pdriver, int stage);
  TaskStatus InitRecv(Driver *pdriver, int stage);
  TaskStatus SendP(Driver *pdriver, int stage);
  TaskStatus RecvP(Driver *pdriver, int stage);
  TaskStatus ClearSend(Driver *pdriver, int stage);
  TaskStatus ClearRecv(Driver *pdriver, int stage);

 private:
  MeshBlockPack* pmy_pack;  // ptr to MeshBlockPack containing this Particles

  void ParseStarGravity(ParameterInput *pin);
  void LoadStarsFromFile(ParameterInput *pin);
  void LoadStarsFromRestart(ParameterInput *pin);
  int FindLocalMeshBlockForPosition(Real x, Real y, Real z) const;
  void RefreshNextStarTag();
  void UpdateParticleCounts();
  void GravityKDKDrift();
  void GravityKDKFinalKick();
  void GravityRK4Step();
  void UpdateGravityDiagnosticsAndTimestep();
  void ComputeStarGravityAcceleration(const std::vector<Real> &x,
                                      const std::vector<Real> &y,
                                      const std::vector<Real> &z,
                                      const std::vector<Real> &mass,
                                      const std::vector<int> &tag,
                                      Real eval_time,
                                      std::vector<Real> &ax,
                                      std::vector<Real> &ay,
                                      std::vector<Real> &az,
                                      Real *min_pair=nullptr,
                                      Real *potential_energy=nullptr);
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
