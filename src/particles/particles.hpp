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
#include <cmath>
#include <cstddef>
#include <cstdint>
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
enum class PICPhysicalMode { engineering = 0, paper_test_particle = 1,
                             paper_mhd_pic = 2, extended_mhd_pic = 3,
                             paper_mhd_pic_vl2_tsc = 4 };
enum class PICCRHallMode { off, current_to_ct_experimental };
enum class PICWaveDampingMode { off, ion_neutral_friction };
enum class PICCRInitialState { velocity, momentum };
enum class PICInterpolationScheme { tsc };
enum class PICDeltaFMode { off, quiet_start, physical };
enum class PICDeltaFBackground { uniform, kappa_iso, kappa_drift, kappa_aniso };
enum class PICDeltaFAdaptMode { off, global_bikappa_moments_experimental };
enum class PICIntermediateArraysMode { auto_mode, off };
enum class PICExpandingBoxMode { off, on };
enum class PICExpansionLaw { linear, reciprocal_linear, exponential };
enum class CRParticleSource { initial = 0, shock_injected = 1 };
enum class Q017ParticleTimer { adaptive_deltaf, push, deposition, migration, count };

//----------------------------------------------------------------------------------------
//! \struct ParticlesTaskIDs
//  \brief container to hold TaskIDs of all particles tasks

struct ParticlesTaskIDs {
  TaskID adapt_deltaf;
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

class Particles;
using ParticleDestructionObserverFnPtr = void (*)(Particles*, Mesh*, int);

KOKKOS_INLINE_FUNCTION
Real CRLorentzFactor(const Real ux, const Real uy, const Real uz,
                     const Real light_speed) {
  return sqrt(1.0 + (ux*ux + uy*uy + uz*uz)/(light_speed*light_speed));
}

KOKKOS_INLINE_FUNCTION
Real CRKineticEnergy(const bool momentum_state, const Real light_speed,
                    const Real sx, const Real sy, const Real sz) {
  if (momentum_state) {
    const Real state_squared = sx*sx + sy*sy + sz*sz;
    return state_squared/(CRLorentzFactor(sx, sy, sz, light_speed) + 1.0);
  }
  return 0.5*(sx*sx + sy*sy + sz*sz);
}

KOKKOS_INLINE_FUNCTION
void CRVelocityFromState(const bool momentum_state, const Real light_speed,
                         const Real sx, const Real sy, const Real sz,
                         Real &vx, Real &vy, Real &vz) {
  const Real inv_gamma = momentum_state ?
      1.0/CRLorentzFactor(sx, sy, sz, light_speed) : 1.0;
  vx = sx*inv_gamma;
  vy = sy*inv_gamma;
  vz = sz*inv_gamma;
}

KOKKOS_INLINE_FUNCTION
Real PICScaleFactor(const PICExpansionLaw law, const Real initial_rate,
                    const Real time) {
  if (law == PICExpansionLaw::linear) {
    return 1.0 + initial_rate*time;
  }
  if (law == PICExpansionLaw::reciprocal_linear) {
    return 1.0/(1.0 + initial_rate*time);
  }
  return exp(initial_rate*time);
}

struct PICExpandingBoxGeometry {
  Real a1, a2, a3;
  Real inv_a1, inv_a2, inv_a3;
  Real inv_area1, inv_area2, inv_area3;
};

KOKKOS_INLINE_FUNCTION
PICExpandingBoxGeometry PICExpandingBoxGeometryAt(
    const PICExpansionLaw law, const Real rate_x1, const Real rate_x2,
    const Real rate_x3, const Real time) {
  const Real a1 = PICScaleFactor(law, rate_x1, time);
  const Real a2 = PICScaleFactor(law, rate_x2, time);
  const Real a3 = PICScaleFactor(law, rate_x3, time);
  return {a1, a2, a3, 1.0/a1, 1.0/a2, 1.0/a3,
          1.0/(a2*a3), 1.0/(a1*a3), 1.0/(a1*a2)};
}

KOKKOS_INLINE_FUNCTION
Real PICDeltaFBackgroundValue(
    const PICDeltaFBackground background, const Real p0, const Real kappa,
    const Real drift_x1, const Real drift_x2, const Real drift_x3,
    const Real aniso_x1, const Real aniso_x2, const Real aniso_x3,
    const Real scale_x1, const Real scale_x2, const Real scale_x3,
    const Real sx, const Real sy, const Real sz) {
  if (background == PICDeltaFBackground::uniform) return 1.0;
  const Real px = (scale_x1*sx - drift_x1)/aniso_x1;
  const Real py = (scale_x2*sy - drift_x2)/aniso_x2;
  const Real pz = (scale_x3*sz - drift_x3)/aniso_x3;
  const Real shape = 1.0 + (px*px + py*py + pz*pz)/(kappa*p0*p0);
  return fmax(pow(shape, -(kappa + 1.0)), static_cast<Real>(1.0e-30));
}

KOKKOS_INLINE_FUNCTION
Real PICAdaptiveDeltaFBackgroundValue(
    const Real kappa, const Real reference_p0, const Real fitted_p0,
    const Real xi, const Real scale_x1, const Real scale_x2,
    const Real scale_x3, const Real sx, const Real sy, const Real sz) {
  const Real xi2 = xi*xi;
  const Real xi4 = xi2*xi2;
  const Real volume = scale_x1*scale_x2*scale_x3;
  const Real shape =
      1.0 + (xi4*sx*sx + xi2*(sy*sy + sz*sz))/(kappa*fitted_p0*fitted_p0);
  const Real normalization =
      (xi4/volume)*pow(reference_p0/fitted_p0, static_cast<Real>(3.0));
  return fmax(normalization*pow(shape, -(kappa + 1.0)),
              static_cast<Real>(1.0e-30));
}

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
  PICPhysicalMode pic_physical_mode = PICPhysicalMode::engineering;
  PICCRHallMode pic_cr_hall_mode = PICCRHallMode::off;
  PICWaveDampingMode pic_wave_damping_mode = PICWaveDampingMode::off;
  PICCRInitialState pic_cr_initial_state = PICCRInitialState::velocity;
  PICInterpolationScheme pic_interp_scheme = PICInterpolationScheme::tsc;
  bool pic_enable_2d3v = false;    // keep vz/Bz channels active when nx3==1
  PICDeltaFMode pic_deltaf_mode = PICDeltaFMode::off;
  PICDeltaFBackground pic_deltaf_background = PICDeltaFBackground::uniform;
  PICDeltaFAdaptMode pic_deltaf_adapt_mode = PICDeltaFAdaptMode::off;
  PICIntermediateArraysMode pic_intermediate_arrays_mode =
      PICIntermediateArraysMode::auto_mode;
  PICExpandingBoxMode pic_expanding_box_mode = PICExpandingBoxMode::off;
  PICExpansionLaw pic_expansion_law = PICExpansionLaw::linear;
  Real pic_cr_light_speed = 1.0;  // artificial CR light speed in momentum-state modes
  Real pic_ion_neutral_collision_rate = 0.0; // reduced high-frequency IN damping rate
  int pic_max_cell_cross = 2;     // particle cell-crossing timestep limit
  Real pic_theta_max = 0.3;       // Boris gyro-angle timestep limit
  int pic_sort_interval = 0;      // staged sorting cadence (0 disables re-sorting)
  int pic_random_seed = 0;        // deterministic seed for random CR placement
  Real pic_load_balance_cost_per_particle = 0.0; // optional AMR balancing cost weight
  bool pic_q017_sync_kernel_timers = false; // opt-in fences for device elapsed timing
  Real pic_expansion_rate_x1 = 0.0;
  Real pic_expansion_rate_x2 = 0.0;
  Real pic_expansion_rate_x3 = 0.0;
  Real pic_deltaf_p0 = 1.0;
  Real pic_deltaf_kappa = 1.25;
  Real pic_deltaf_adapt_interval = 0.0;
  Real pic_deltaf_adaptive_xi = 1.0;
  Real pic_deltaf_adaptive_p0 = 1.0;
  std::int64_t pic_deltaf_adapt_last_bucket = -1;
  Real pic_deltaf_drift_x1 = 0.0;
  Real pic_deltaf_drift_x2 = 0.0;
  Real pic_deltaf_drift_x3 = 0.0;
  Real pic_deltaf_aniso_x1 = 1.0;
  Real pic_deltaf_aniso_x2 = 1.0;
  Real pic_deltaf_aniso_x3 = 1.0;
  Real pic_deltaf_background_rho = 0.0;
  Real pic_deltaf_background_jx = 0.0;
  Real pic_deltaf_background_jy = 0.0;
  Real pic_deltaf_background_jz = 0.0;
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
  DvceArray1D<PaperSmoothMomentRecord> paper_smooth_mom_records;
  DvceArray4D<Real> j_edge_x1e, j_edge_x2e, j_edge_x3e;
  DvceArray1D<Real> x1_old, x2_old, x3_old;

  // Constants for rk4_gravity pusher
  Real r_scale = 0.0;
  Real rho_scale = 0.0;
  Real m_gal = 0.0;
  Real a_gal = 0.0;
  Real z_gal = 0.0;
  Real r_200 = 0.0;
  Real rho_mean = 0.0;
  Real par_grav_dx = 1.0e-6;

  // Boundary communication buffers and functions for particles
  ParticlesBoundaryValues *pbval_part;
  MeshBoundaryValuesCC *pbval_mom = nullptr;
  MeshBoundaryValuesFC *pbval_jedge = nullptr;
  PaperSmoothMomentRecordTransport *paper_smooth_mom_transport = nullptr;
  ParticleDestructionObserverFnPtr particle_destruction_observer = nullptr;

  // container to hold names of TaskIDs
  ParticlesTaskIDs id;

  // functions...
  void CreateParticleTags(ParameterInput *pin);
  void UpdateAfterAMR(MeshBlockPack *new_pp);
  void AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  TaskStatus Push(Driver *pdriver, int stage);
  TaskStatus AdaptDeltaF(Driver *pdriver, int stage);
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
  TaskStatus DepositPaperSmoothMoments(Driver *pdriver, int stage);
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
  TaskStatus PushPaperCosmicRaysVL2(Driver *pdriver, int stage);
  TaskStatus DriftPaperCosmicRaysHalfStep(Driver *pdriver, int stage);
  TaskStatus PushStars(Driver *pdriver, int stage);
  void NewTimeStep();
  void Q017Fence() const {
    if (pic_q017_sync_kernel_timers) Kokkos::fence();
  }
  void AccumulateQ017Timer(Q017ParticleTimer timer, double seconds) {
    const int n = static_cast<int>(timer);
    q017_particle_time_[n] += seconds;
    q017_particle_calls_[n]++;
  }
  std::uint64_t Q017DirectViewAllocationBytes() const;
  std::uint64_t Q017OwnedKokkosViewAllocationBytes() const;
  void ObserveQ017OwnedKokkosViewAllocationBytes(std::uint64_t transient_bytes=0);
  std::uint64_t Q017PaperSmoothHostAllocationBytes() const;
  void ObserveQ017PaperSmoothHostAllocationBytes(std::uint64_t transient_bytes=0);
  void OutputQ017Telemetry() const;
  bool UsesRelativisticCRState() const {
    return pic_physical_mode != PICPhysicalMode::engineering;
  }
  bool UsesDeltaF() const {
    return pic_deltaf_mode == PICDeltaFMode::physical;
  }
  bool UsesAdaptiveDeltaF() const {
    return pic_deltaf_adapt_mode ==
           PICDeltaFAdaptMode::global_bikappa_moments_experimental;
  }
  bool UsesExpandingBox() const {
    return pic_expanding_box_mode == PICExpandingBoxMode::on;
  }
  bool UsesPaperVL2Coupling() const {
    return pic_physical_mode == PICPhysicalMode::paper_mhd_pic_vl2_tsc;
  }
  bool UsesPICWaveDamping() const {
    return pic_wave_damping_mode == PICWaveDampingMode::ion_neutral_friction;
  }
  static constexpr int PIC_RESTART_SCHEMA_VERSION = 7;
  static constexpr int NPIC_RESTART_MODEL_INTS = 31;
  static constexpr int NPIC_RESTART_CONFIG_REALS = 34;
  static constexpr int NPIC_RESTART_MODEL_REALS = 37;
  std::uint64_t RestartSpeciesConfigHash() const {
    std::uint64_t hash = 14695981039346656037ULL;
    auto hash_bytes = [&hash](const auto &value) {
      const auto *bytes = reinterpret_cast<const unsigned char *>(&value);
      for (std::size_t n=0; n<sizeof(value); ++n) {
        hash ^= bytes[n];
        hash *= 1099511628211ULL;
      }
    };
    const int restart_nspecies =
        (particle_type == ParticleType::cosmic_ray) ? nspecies : 0;
    hash_bytes(restart_nspecies);
    if (particle_type == ParticleType::cosmic_ray) {
      auto h_mass = Kokkos::create_mirror_view_and_copy(HostMemSpace(), species_mass);
      auto h_charge = Kokkos::create_mirror_view_and_copy(HostMemSpace(), species_charge);
      for (int s=0; s<nspecies; ++s) {
        hash_bytes(h_mass(s));
        hash_bytes(h_charge(s));
      }
    }
    return hash;
  }
  void FillRestartModelMetadata(
      std::array<int, NPIC_RESTART_MODEL_INTS> &model_ints,
      std::array<Real, NPIC_RESTART_MODEL_REALS> &model_reals) const {
    const std::uint64_t species_hash = RestartSpeciesConfigHash();
    model_ints = {static_cast<int>(pic_deltaf_mode),
                  static_cast<int>(pic_deltaf_background),
                  static_cast<int>(pic_expanding_box_mode),
                  static_cast<int>(pic_expansion_law),
                  static_cast<int>(pic_wave_damping_mode),
                  static_cast<int>(pic_deltaf_adapt_mode),
                  static_cast<int>(particle_type),
                  static_cast<int>(pusher),
                  (particle_type == ParticleType::cosmic_ray) ? nspecies : 0,
                  (particle_type == ParticleType::cosmic_ray &&
                   track_displacement) ? 1 : 0,
                  deposit_moments ? 1 : 0,
                  deposit_order,
                  couple_moments_to_mhd ? 1 : 0,
                  static_cast<int>(couple_j_to_efield_representation),
                  static_cast<int>(couple_j_deposition_mode),
                  static_cast<int>(couple_fluid_feedback_order),
                  couple_moments_momentum_to_mhd ? 1 : 0,
                  couple_moments_energy_to_mhd ? 1 : 0,
                  static_cast<int>(pic_background_mode),
                  static_cast<int>(pic_feedback_mode),
                  static_cast<int>(pic_cr_hall_mode),
                  static_cast<int>(pic_cr_initial_state),
                  static_cast<int>(pic_interp_scheme),
                  pic_enable_2d3v ? 1 : 0,
                  static_cast<int>(pic_intermediate_arrays_mode),
                  pic_max_cell_cross,
                  pic_sort_interval,
                  pic_random_seed,
                  static_cast<int>(species_hash & 0x3fffffULL),
                  static_cast<int>((species_hash >> 22) & 0x3fffffULL),
                  static_cast<int>((species_hash >> 44) & 0xfffffULL)};
    model_reals = {pic_expansion_rate_x1, pic_expansion_rate_x2,
                   pic_expansion_rate_x3, pic_deltaf_p0, pic_deltaf_kappa,
                   pic_deltaf_drift_x1, pic_deltaf_drift_x2, pic_deltaf_drift_x3,
                   pic_deltaf_aniso_x1, pic_deltaf_aniso_x2, pic_deltaf_aniso_x3,
                   pic_deltaf_background_rho, pic_deltaf_background_jx,
                   pic_deltaf_background_jy, pic_deltaf_background_jz,
                   pic_no_mhd_bx, pic_no_mhd_by, pic_no_mhd_bz,
                   pic_ion_neutral_collision_rate, pic_deltaf_adapt_interval,
                   deposit_qscale, couple_j_to_efield_coeff,
                   couple_moments_momentum_coeff, couple_moments_energy_coeff,
                   pic_theta_max, pic_load_balance_cost_per_particle,
                   r_scale, rho_scale, m_gal, a_gal, z_gal, r_200, rho_mean,
                   par_grav_dx,
                   pic_deltaf_adaptive_xi, pic_deltaf_adaptive_p0,
                   static_cast<Real>(pic_deltaf_adapt_last_bucket)};
  }
  bool MatchesRestartModelMetadata(
      const std::array<int, NPIC_RESTART_MODEL_INTS> &model_ints,
      const std::array<Real, NPIC_RESTART_MODEL_REALS> &model_reals) const {
    std::array<int, NPIC_RESTART_MODEL_INTS> expected_ints;
    std::array<Real, NPIC_RESTART_MODEL_REALS> expected_reals;
    FillRestartModelMetadata(expected_ints, expected_reals);
    if (model_ints != expected_ints) return false;
    for (int n=0; n<NPIC_RESTART_CONFIG_REALS; ++n) {
      if (model_reals[n] != expected_reals[n]) return false;
    }
    return true;
  }
  bool RestoreRestartModelState(
      const std::array<Real, NPIC_RESTART_MODEL_REALS> &model_reals) {
    if (!UsesAdaptiveDeltaF()) return true;
    constexpr int xi_index = NPIC_RESTART_CONFIG_REALS;
    constexpr int p0_index = NPIC_RESTART_CONFIG_REALS + 1;
    constexpr int bucket_index = NPIC_RESTART_CONFIG_REALS + 2;
    const Real xi = model_reals[xi_index];
    const Real p0 = model_reals[p0_index];
    const Real bucket_real = model_reals[bucket_index];
    const std::int64_t bucket = static_cast<std::int64_t>(bucket_real);
    if (!std::isfinite(xi) || xi <= 0.0 ||
        !std::isfinite(p0) || p0 <= 0.0 ||
        !std::isfinite(bucket_real) || static_cast<Real>(bucket) != bucket_real ||
        bucket < -1) {
      return false;
    }
    if (pic_deltaf_adapt_last_bucket >= 0 &&
        (pic_deltaf_adaptive_xi != xi || pic_deltaf_adaptive_p0 != p0 ||
         pic_deltaf_adapt_last_bucket != bucket)) {
      return false;
    }
    pic_deltaf_adaptive_xi = xi;
    pic_deltaf_adaptive_p0 = p0;
    pic_deltaf_adapt_last_bucket = bucket;
    return true;
  }
  bool AddsCRCurrentToCT() const {
    return couple_moments_to_mhd &&
           ((pic_physical_mode == PICPhysicalMode::engineering) ||
            ((pic_physical_mode == PICPhysicalMode::extended_mhd_pic) &&
             (pic_cr_hall_mode == PICCRHallMode::current_to_ct_experimental)));
  }

 private:
  MeshBlockPack *pmy_pack; // ptr to MeshBlockPack containing this Particles
  static constexpr int nq017_particle_timers =
      static_cast<int>(Q017ParticleTimer::count);
  std::array<double, nq017_particle_timers> q017_particle_time_{};
  std::array<std::uint64_t, nq017_particle_timers> q017_particle_calls_{};
  std::uint64_t q017_owned_kokkos_view_high_water_bytes_ = 0;
  std::uint64_t q017_paper_smooth_host_high_water_bytes_ = 0;
};

} // namespace particles
#endif // PARTICLES_PARTICLES_HPP_
