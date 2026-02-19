#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.hpp
//  \brief definitions for Particles class

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

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
enum class ParticleType { cosmic_ray, star };

// constants for PR2 current representation used by coupled E-field source updates
enum class CoupledCurrentRepresentation { cell_centered, edge_staggered };
enum class CoupledCurrentDepositionMode { cc_convert, direct_staggered };
enum class CoupledFluidFeedbackOrder { mhd_src_terms, efield_src };

// constants for staged PIC runtime controls used by PR5+ test-suite expansion
enum class PICBackgroundMode { coupled, passive_mhd, no_mhd };
enum class PICFeedbackMode { coupled, test_particle };
enum class PICInterpolationScheme { tsc };
enum class PICDeltaFMode { off, on };
enum class PICIntermediateArraysMode { auto_mode, off };
enum class PICExpandingBoxMode { off, on };

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
  TaskID save_old;
  TaskID zero_mom;
  TaskID irecv_mom;
  TaskID dep_mom;
  TaskID rest_mom;
  TaskID send_mom;
  TaskID recv_mom;
  TaskID crecv_mom;
  TaskID csend_mom;
  TaskID bcs_mom;
  TaskID prol_mom;
  TaskID irecv_jedge;
  TaskID send_jedge;
  TaskID recv_jedge;
  TaskID crecv_jedge;
  TaskID csend_jedge;
  TaskID bcs_jedge;
  TaskID convert_j_edge;
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

  // Cosmic ray specific
  int nspecies;                     // number of CR species
  bool track_displacement;          // enable displacement tracking
  DvceArray1D<Real> species_mass;   // mass per species
  DvceArray1D<Real> species_charge; // charge per species
  DvceArray1D<Real> species_vx0;    // optional per-species vx initializer
  DvceArray1D<Real> species_vy0;    // optional per-species vy initializer
  DvceArray1D<Real> species_vz0;    // optional per-species vz initializer
  bool deposit_moments = false;     // enable particle moment deposition
  int deposit_order = 1;            // deposition shape order
  Real deposit_qscale = 1.0;        // scaling of particle macro-charge
  bool couple_moments_to_mhd = false;  // PR2 opt-in current coupling to MHD
  Real couple_j_to_efield_coeff = 1.0; // PR2 current-to-E coupling coefficient
  CoupledCurrentRepresentation couple_j_to_efield_representation =
      CoupledCurrentRepresentation::cell_centered;
  CoupledCurrentDepositionMode couple_j_deposition_mode =
      CoupledCurrentDepositionMode::cc_convert;
  CoupledFluidFeedbackOrder couple_fluid_feedback_order =
      CoupledFluidFeedbackOrder::mhd_src_terms;
  bool couple_moments_momentum_to_mhd = false;  // PR2 opt-in momentum feedback
  bool couple_moments_energy_to_mhd = false;    // PR2 opt-in energy feedback
  Real couple_moments_momentum_coeff = 1.0;     // momentum feedback coefficient
  Real couple_moments_energy_coeff = 1.0;       // energy feedback coefficient
  PICBackgroundMode pic_background_mode = PICBackgroundMode::coupled;
  PICFeedbackMode pic_feedback_mode = PICFeedbackMode::coupled;
  PICInterpolationScheme pic_interp_scheme = PICInterpolationScheme::tsc;
  bool pic_enable_2d3v = false;    // keep vz/Bz channels active when nx3==1
  PICDeltaFMode pic_deltaf_mode = PICDeltaFMode::off;
  PICIntermediateArraysMode pic_intermediate_arrays_mode =
      PICIntermediateArraysMode::auto_mode;
  PICExpandingBoxMode pic_expanding_box_mode = PICExpandingBoxMode::off;
  Real pic_cr_light_speed = 1.0;  // artificial particle light speed (staged use)
  int pic_max_cell_cross = 2;     // staged guard control (cells per step)
  Real pic_theta_max = 0.3;       // staged guard control (gyro angle per step)
  int pic_sort_interval = 0;      // staged sorting cadence (0 disables re-sorting)
  Real pic_expansion_rate_x1 = 0.0;
  Real pic_expansion_rate_x2 = 0.0;
  Real pic_expansion_rate_x3 = 0.0;
  Real pic_no_mhd_bx = 0.0;
  Real pic_no_mhd_by = 0.0;
  Real pic_no_mhd_bz = 0.0;
  std::string pic_deltaf_f0 = "";
  DvceArray5D<Real> pic_no_mhd_bcc0;
  Real cr_vx0 = 0.0;                // deterministic CR vx initialization
  Real cr_vy0 = 0.0;                // deterministic CR vy initialization
  Real cr_vz0 = 0.0;                // deterministic CR vz initialization
  static constexpr int NMOM = 9;
  static constexpr int IMOM_RHO = 0;
  static constexpr int IMOM_JX  = 1;
  static constexpr int IMOM_JY  = 2;
  static constexpr int IMOM_JZ  = 3;
  static constexpr int IMOM_DPXDT = 4;
  static constexpr int IMOM_DPYDT = 5;
  static constexpr int IMOM_DPZDT = 6;
  static constexpr int IMOM_DEDT = 7;
  static constexpr int IMOM_EBDOT = 8;
  DvceArray5D<Real> moments;
  DvceArray5D<Real> coarse_moments;
  DvceArray4D<Real> j_edge_x1e, j_edge_x2e, j_edge_x3e;
  DvceArray1D<Real> x1_old, x2_old, x3_old;

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
  MeshBoundaryValuesCC *pbval_mom = nullptr;
  MeshBoundaryValuesFC *pbval_jedge = nullptr;

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
  TaskStatus SaveOldPositions(Driver *pdriver, int stage);
  TaskStatus ZeroMoments(Driver *pdriver, int stage);
  TaskStatus InitRecvMoments(Driver *pdriver, int stage);
  TaskStatus DepositMoments(Driver *pdriver, int stage);
  TaskStatus RestrictMoments(Driver *pdriver, int stage);
  TaskStatus SendMoments(Driver *pdriver, int stage);
  TaskStatus RecvMoments(Driver *pdriver, int stage);
  TaskStatus ClearRecvMoments(Driver *pdriver, int stage);
  TaskStatus ClearSendMoments(Driver *pdriver, int stage);
  TaskStatus ApplyMomentPhysicalBCs(Driver *pdriver, int stage);
  TaskStatus ProlongateMoments(Driver *pdriver, int stage);
  TaskStatus InitRecvEdgeCurrents(Driver *pdriver, int stage);
  TaskStatus SendEdgeCurrents(Driver *pdriver, int stage);
  TaskStatus RecvEdgeCurrents(Driver *pdriver, int stage);
  TaskStatus ClearRecvEdgeCurrents(Driver *pdriver, int stage);
  TaskStatus ClearSendEdgeCurrents(Driver *pdriver, int stage);
  TaskStatus ApplyEdgeCurrentPhysicalBCs(Driver *pdriver, int stage);
  TaskStatus ConvertCoupledCurrentRepresentation(Driver *pdriver, int stage);

  // Cosmic ray specific methods
  void InitializeCosmicRays(ParameterInput *pin);
  void InitializeStars(std::vector<std::array<Real, 9>> &particle_list);
  TaskStatus PushDrift(Driver *pdriver, int stage);
  TaskStatus PushCosmicRays(Driver *pdriver, int stage);
  TaskStatus PushStars(Driver *pdriver, int stage);

 private:
  MeshBlockPack *pmy_pack; // ptr to MeshBlockPack containing this Particles
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
