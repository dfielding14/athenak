//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q019_physics_first_nonlinear_bell_successor_v2.hpp
//! \brief Shared math for the physics-first nonlinear Bell successor.

#ifndef PGEN_TESTS_Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HPP_
#define PGEN_TESTS_Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HPP_

#include <cmath>
#include <cstdint>

#if !defined(Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT)
#include "athena.hpp"
#endif

#if defined(Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT)
#define Q019_NLB_INLINE inline
#else
#define Q019_NLB_INLINE KOKKOS_INLINE_FUNCTION
#endif

namespace q019_physics_first_nonlinear_bell_successor_v2 {

enum class Branch {
  high_rigidity_current_retention_candidate,
  finite_rigidity_self_consistent,
  finite_rigidity_early_time_predecessor
};

enum class SeedTopology {
  shared_spectrum_only,
  shared_spectrum_plus_separated_long_modes
};

enum class FiniteSamplingMode {
  not_applicable,
  cell_centered_nested_haar_octahedral_packets,
  independent_position_noise_seeded
};

struct Vector3 {
  double x1;
  double x2;
  double x3;
};

struct IntegerVector3 {
  int x1;
  int x2;
  int x3;
};

struct ComplexValue {
  double real;
  double imag;
};

enum class BoxEdgeDiagnosticStatus {
  not_sampled = 0,
  valid = 1,
  cadence_skipped = 2,
  no_global_cells = 3,
  zero_fluctuation = 4,
  numerical_unavailable = 5
};

struct SeedParameters {
  int dimension;
  int field_seed;
  double b_g;
  double k0;
  double epsilon;
  double rho;
  double eigenmode_amplitude;
  double broadband_amplitude;
  double separated_long_mode_amplitude;
  double x1_extent;
  double x2_extent;
  double x3_extent;
  SeedTopology seed_topology;
};

Q019_NLB_INLINE Vector3 Add(const Vector3 a, const Vector3 b) {
  return {a.x1 + b.x1, a.x2 + b.x2, a.x3 + b.x3};
}

Q019_NLB_INLINE Vector3 Subtract(const Vector3 a, const Vector3 b) {
  return {a.x1 - b.x1, a.x2 - b.x2, a.x3 - b.x3};
}

Q019_NLB_INLINE Vector3 Scale(const Vector3 value, const double factor) {
  return {factor*value.x1, factor*value.x2, factor*value.x3};
}

constexpr int kBoxEdgeUniqueModeCount2D = 7;
constexpr int kBoxEdgeUniqueModeCount3D = 23;
constexpr int kBoxEdgeMaximumUniqueModeCount = kBoxEdgeUniqueModeCount3D;
constexpr int kBoxEdgeModesPerReduction = 5;

Q019_NLB_INLINE int BoxEdgeUniqueModeCount(const int dimension) {
  return dimension == 2 ? kBoxEdgeUniqueModeCount2D :
      dimension == 3 ? kBoxEdgeUniqueModeCount3D : 0;
}

Q019_NLB_INLINE int BoxEdgeReductionGroupCount(const int dimension) {
  const int modes = BoxEdgeUniqueModeCount(dimension);
  return (modes + kBoxEdgeModesPerReduction - 1)/kBoxEdgeModesPerReduction;
}

Q019_NLB_INLINE IntegerVector3 BoxEdgeUniqueModeAt(
    const int index, const int dimension) {
  if (dimension == 2) {
    if (index == 0) return {0, 1, 0};
    if (index >= 1 && index <= 3) return {1, index - 2, 0};
    if (index >= 4 && index <= 6) return {2, index - 5, 0};
    return {0, 0, 0};
  }
  if (index == 0) return {0, 0, 1};
  if (index >= 1 && index <= 3) return {0, 1, index - 2};
  if (index >= 4 && index <= 12) {
    const int offset = index - 4;
    return {1, offset/3 - 1, offset%3 - 1};
  }
  if (index >= 13 && index <= 21) {
    const int offset = index - 13;
    return {2, offset/3 - 1, offset%3 - 1};
  }
  if (index == 22) return {3, 0, 0};
  return {0, 0, 0};
}

Q019_NLB_INLINE bool BoxEdgeFrozenAspectRatioIsValid(
    const int dimension, const double x1_extent, const double x2_extent,
    const double x3_extent) {
  if (!std::isfinite(x1_extent) || !std::isfinite(x2_extent) ||
      !std::isfinite(x3_extent) ||
      !(x1_extent > 0.0 && x2_extent > 0.0 && x3_extent > 0.0)) {
    return false;
  }
  const double tolerance = 1.0e-12;
  const bool x1_ratio =
      std::fabs(x1_extent/x2_extent - 2.0) <= tolerance;
  const bool transverse_ratio =
      dimension == 2 || (dimension == 3 &&
                         std::fabs(x2_extent/x3_extent - 1.0) <= tolerance);
  return (dimension == 2 || dimension == 3) && x1_ratio && transverse_ratio;
}

Q019_NLB_INLINE double BoxEdgePhysicalLowKCutoffSquaredWithoutTwoPi(
    const int dimension, const double x1_extent, const double x2_extent,
    const double x3_extent) {
  if (!BoxEdgeFrozenAspectRatioIsValid(
          dimension, x1_extent, x2_extent, x3_extent)) {
    return -1.0;
  }
  double minimum_extent = x1_extent < x2_extent ? x1_extent : x2_extent;
  if (dimension == 3 && x3_extent < minimum_extent) minimum_extent = x3_extent;
  return static_cast<double>(dimension)/(minimum_extent*minimum_extent);
}

Q019_NLB_INLINE double BoxEdgeModeWaveNumberSquaredWithoutTwoPi(
    const IntegerVector3 mode, const int dimension, const double x1_extent,
    const double x2_extent, const double x3_extent) {
  if (dimension != 2 && dimension != 3) return -1.0;
  const double k1 = static_cast<double>(mode.x1)/x1_extent;
  const double k2 = static_cast<double>(mode.x2)/x2_extent;
  const double k3 = dimension == 3 ?
      static_cast<double>(mode.x3)/x3_extent : 0.0;
  return k1*k1 + k2*k2 + k3*k3;
}

Q019_NLB_INLINE bool BoxEdgePhysicalLowKModeIsSelected(
    const IntegerVector3 mode, const int dimension, const double x1_extent,
    const double x2_extent, const double x3_extent) {
  if (mode.x1 == 0 && mode.x2 == 0 && mode.x3 == 0) return false;
  if (dimension == 2 && mode.x3 != 0) return false;
  const double cutoff_squared = BoxEdgePhysicalLowKCutoffSquaredWithoutTwoPi(
      dimension, x1_extent, x2_extent, x3_extent);
  const double mode_squared = BoxEdgeModeWaveNumberSquaredWithoutTwoPi(
      mode, dimension, x1_extent, x2_extent, x3_extent);
  if (!(cutoff_squared > 0.0 && mode_squared >= 0.0)) return false;
  return mode_squared <= cutoff_squared*(1.0 + 64.0*2.2204460492503131e-16);
}

Q019_NLB_INLINE bool BoxEdgeStatusCodeIsValid(const int status) {
  return status >= static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled) &&
      status <= static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable);
}

Q019_NLB_INLINE int BoxEdgeStatusBit(const int status) {
  return BoxEdgeStatusCodeIsValid(status) ? (1 << status) : 0;
}

Q019_NLB_INLINE bool BoxEdgeStatusRecordsValidMetric(const int status) {
  return status == static_cast<int>(BoxEdgeDiagnosticStatus::valid);
}

Q019_NLB_INLINE ComplexValue ComplexMultiply(
    const ComplexValue lhs, const ComplexValue rhs) {
  return {
    lhs.real*rhs.real - lhs.imag*rhs.imag,
    lhs.real*rhs.imag + lhs.imag*rhs.real
  };
}

Q019_NLB_INLINE ComplexValue ComplexConjugate(const ComplexValue value) {
  return {value.real, -value.imag};
}

Q019_NLB_INLINE ComplexValue SelectUnitComplexPower(
    const ComplexValue positive, const int exponent) {
  const ComplexValue factor =
      exponent >= 0 ? positive : ComplexConjugate(positive);
  ComplexValue value = {1.0, 0.0};
  const int count = exponent >= 0 ? exponent : -exponent;
  for (int n = 0; n < count; ++n) value = ComplexMultiply(value, factor);
  return value;
}

Q019_NLB_INLINE ComplexValue BoxEdgeModePhase(
    const ComplexValue x1_phase, const ComplexValue x2_phase,
    const ComplexValue x3_phase, const IntegerVector3 mode) {
  return ComplexMultiply(
      ComplexMultiply(SelectUnitComplexPower(x1_phase, mode.x1),
                      SelectUnitComplexPower(x2_phase, mode.x2)),
      SelectUnitComplexPower(x3_phase, mode.x3));
}

Q019_NLB_INLINE double BoxEdgePowerFraction(
    const double cell_count, const double transverse_fluctuation_sum,
    const double unique_mode_power_sum) {
  if (!(cell_count > 0.0 && transverse_fluctuation_sum > 0.0 &&
        unique_mode_power_sum >= 0.0)) {
    return -1.0;
  }
  return 2.0*unique_mode_power_sum/(cell_count*transverse_fluctuation_sum);
}

Q019_NLB_INLINE double BoxEdgeNextNominalTime(
    const int completed_slots, const double monitor_dt) {
  return static_cast<double>(completed_slots + 1)*monitor_dt;
}

Q019_NLB_INLINE double BoxEdgeLastNominalTime(
    const int completed_slots, const double monitor_dt) {
  return completed_slots > 0
      ? static_cast<double>(completed_slots)*monitor_dt : -1.0;
}

Q019_NLB_INLINE int BoxEdgeCrossedSlotCount(
    const int completed_slots, const double monitor_dt,
    const double completed_time) {
  if (completed_slots < 0 || !(monitor_dt > 0.0) || !std::isfinite(monitor_dt) ||
      !std::isfinite(completed_time) || completed_time < 0.0) {
    return -1;
  }
  const double next_nominal_time = BoxEdgeNextNominalTime(
      completed_slots, monitor_dt);
  double scale = std::fabs(completed_time);
  if (std::fabs(next_nominal_time) > scale) scale = std::fabs(next_nominal_time);
  if (scale < 1.0) scale = 1.0;
  const double tolerance = 64.0*2.2204460492503131e-16*scale;
  if (completed_time + tolerance < next_nominal_time) return 0;
  const double last_crossed_slot =
      std::floor((completed_time + tolerance)/monitor_dt);
  if (last_crossed_slot > 2.0e9) return -1;
  const int crossed_slots =
      static_cast<int>(last_crossed_slot) - completed_slots;
  return crossed_slots >= 0 ? crossed_slots : -1;
}

Q019_NLB_INLINE bool BoxEdgeFirstCrossingChronologyIsValid(
    const double prior_completed_time, const double nominal_time,
    const double completed_time) {
  if (!std::isfinite(prior_completed_time) || !std::isfinite(nominal_time) ||
      !std::isfinite(completed_time) || prior_completed_time < 0.0 ||
      nominal_time < 0.0 || completed_time < 0.0 ||
      !(prior_completed_time < completed_time)) {
    return false;
  }
  double scale = std::fabs(completed_time);
  if (std::fabs(nominal_time) > scale) scale = std::fabs(nominal_time);
  if (scale < 1.0) scale = 1.0;
  const double tolerance = 64.0*2.2204460492503131e-16*scale;
  return prior_completed_time < nominal_time &&
      nominal_time <= completed_time + tolerance;
}

Q019_NLB_INLINE double RootCellVolume(const double x1_extent, const int root_nx1,
                                      const double x2_extent, const int root_nx2,
                                      const double x3_extent, const int root_nx3) {
  return (x1_extent/static_cast<double>(root_nx1))*
      (x2_extent/static_cast<double>(root_nx2))*
      (x3_extent/static_cast<double>(root_nx3));
}

Q019_NLB_INLINE double RequiredDepositedJOverC(const double b_g, const double k0) {
  return 2.0*b_g*k0;
}

Q019_NLB_INLINE double DepositedJOverC(const double ppc, const double qscale,
                                       const double species_charge,
                                       const double guide_parallel_stream_speed,
                                       const double root_cell_volume) {
  return ppc*qscale*species_charge*guide_parallel_stream_speed/root_cell_volume;
}

Q019_NLB_INLINE double CRMassDensityOverRho0(const double ppc, const double qscale,
                                             const double root_cell_volume,
                                             const double rho0) {
  return ppc*qscale/(root_cell_volume*rho0);
}

Q019_NLB_INLINE double CRNumberDensity(const double rho_cr,
                                       const double species_mass) {
  return rho_cr/species_mass;
}

Q019_NLB_INLINE double ChargeDensityRatio(
    const double rho_cr_over_rho0, const double species_q_over_mc,
    const double background_q_over_mc) {
  const double cr_to_background_charge_density =
      rho_cr_over_rho0*species_q_over_mc/background_q_over_mc;
  return cr_to_background_charge_density/
      (1.0 + cr_to_background_charge_density);
}

Q019_NLB_INLINE double HallParameter(
    const double charge_density_ratio, const double relative_drift,
    const double alfven_speed) {
  return charge_density_ratio*relative_drift/alfven_speed;
}

Q019_NLB_INLINE double BackgroundIonGyrofrequency(
    const double background_q_over_mc, const double b_g) {
  return background_q_over_mc*b_g;
}

Q019_NLB_INLINE double BackgroundIonInertialLength(
    const double alfven_speed, const double background_ion_gyrofrequency) {
  return alfven_speed/background_ion_gyrofrequency;
}

Q019_NLB_INLINE double HallParameterFromCurrent(
    const double j_over_c, const double rho0, const double background_q_over_mc,
    const double alfven_speed, const double charge_density_ratio) {
  return j_over_c/(rho0*background_q_over_mc*alfven_speed)*
      (1.0 - charge_density_ratio);
}

Q019_NLB_INLINE double BaiHallLinearFactor(const double hall_parameter) {
  return 1.0 + 0.25*hall_parameter*hall_parameter;
}

Q019_NLB_INLINE double BaiHallGrowthRateFractionalShift(const double hall_parameter) {
  return 1.0/std::sqrt(BaiHallLinearFactor(hall_parameter)) - 1.0;
}

Q019_NLB_INLINE double BaiHallWavenumberFractionalShift(const double hall_parameter) {
  return 1.0/BaiHallLinearFactor(hall_parameter) - 1.0;
}

Q019_NLB_INLINE double BaiHallGrowthRateReductionFactor(const double hall_parameter) {
  return 1.0/std::sqrt(BaiHallLinearFactor(hall_parameter));
}

Q019_NLB_INLINE double BaiHallWavenumberReductionFactor(const double hall_parameter) {
  return 1.0/BaiHallLinearFactor(hall_parameter);
}

Q019_NLB_INLINE double CRInertiaParameter(const double rho_cr_over_rho0) {
  return rho_cr_over_rho0;
}

Q019_NLB_INLINE double CRMomentumLoadingParameter(
    const double rho_cr_over_rho0, const double relative_drift,
    const double alfven_speed) {
  return rho_cr_over_rho0*relative_drift/alfven_speed;
}

Q019_NLB_INLINE double CRRMSSpeedKineticLoadingProxyToBackgroundMagneticEnergy(
    const double rho_cr_over_rho0, const double rho0, const double total_speed,
    const double light_speed, const double b_g) {
  const double gamma = 1.0/std::sqrt(
      1.0 - total_speed*total_speed/(light_speed*light_speed));
  return rho_cr_over_rho0*rho0*light_speed*light_speed*(gamma - 1.0)/
      (0.5*b_g*b_g);
}

Q019_NLB_INLINE double FeedbackForceParameter(
    const double j_over_c, const double b_g, const double rho0,
    const double alfven_speed, const double k0) {
  return j_over_c*b_g/(rho0*alfven_speed*alfven_speed*k0);
}

Q019_NLB_INLINE bool LoadingAccountingIsFiniteAndPositive(
    const double inertia_parameter, const double momentum_loading_parameter) {
  return std::isfinite(inertia_parameter) &&
      std::isfinite(momentum_loading_parameter) && inertia_parameter > 0.0 &&
      momentum_loading_parameter > 0.0;
}

Q019_NLB_INLINE Vector3 DepositedSpeciesSumJOverC(
    const double ppc, const int nspecies, const double qscale,
    const Vector3 summed_species_charge_times_velocity,
    const double root_cell_volume) {
  const double scale = ppc*qscale/
      (static_cast<double>(nspecies)*root_cell_volume);
  return Scale(summed_species_charge_times_velocity, scale);
}

Q019_NLB_INLINE bool PositiveIntegralPPCIsValid(const double ppc) {
  return std::isfinite(ppc) && ppc >= 1.0 && std::floor(ppc) == ppc;
}

Q019_NLB_INLINE bool FinitePPCIsSupported(const double ppc) {
  return ppc == 24.0 || ppc == 48.0 || ppc == 96.0;
}

Q019_NLB_INLINE bool HighPPCIsSupported(const double ppc) {
  return ppc == 1.0 || ppc == 8.0 || ppc == 32.0;
}

Q019_NLB_INLINE bool FiniteSamplingModeIsValid(
    const Branch branch, const FiniteSamplingMode mode) {
  if (branch == Branch::high_rigidity_current_retention_candidate) {
    return mode == FiniteSamplingMode::not_applicable;
  }
  return mode == FiniteSamplingMode::cell_centered_nested_haar_octahedral_packets ||
      mode == FiniteSamplingMode::independent_position_noise_seeded;
}

Q019_NLB_INLINE bool OctahedralShellPacketLayoutIsValid(
    const int ppc, const int cells_per_meshblock, const int meshblocks,
    const int particles_this_pack) {
  return ppc > 0 && cells_per_meshblock > 0 && meshblocks > 0 &&
      ppc % 6 == 0 &&
      particles_this_pack == ppc*cells_per_meshblock*meshblocks &&
      (ppc*cells_per_meshblock) % ppc == 0;
}

Q019_NLB_INLINE bool BranchSamplingIsValid(const Branch branch, const double ppc,
                                            const int nspecies,
                                            const bool random_distribution) {
  if (!PositiveIntegralPPCIsValid(ppc)) return false;
  if (branch == Branch::high_rigidity_current_retention_candidate) {
    return HighPPCIsSupported(ppc) && nspecies == 1;
  }
  return FinitePPCIsSupported(ppc) && nspecies == 1 && random_distribution;
}

Q019_NLB_INLINE bool BranchRigidityIsValid(const Branch branch,
                                            const double k0_rg0) {
  if (!std::isfinite(k0_rg0) || k0_rg0 <= 0.0) return false;
  if (branch == Branch::high_rigidity_current_retention_candidate) {
    return k0_rg0 >= 128.0;
  }
  return k0_rg0 >= 4.0 && k0_rg0 <= 16.0;
}

Q019_NLB_INLINE bool ResolutionEnvelopeIsValid(
    const double initial_rg0_over_dx, const double minimum_evolving_rl_over_dx) {
  return std::isfinite(initial_rg0_over_dx) &&
      std::isfinite(minimum_evolving_rl_over_dx) &&
      minimum_evolving_rl_over_dx > 0.0 &&
      initial_rg0_over_dx >= minimum_evolving_rl_over_dx;
}

Q019_NLB_INLINE double SeedPhase(const int field_seed, const int mode) {
  const int reduced = (field_seed*104729 + mode*13007 + 7919) % 1000003;
  return 2.0*M_PI*static_cast<double>(reduced)/1000003.0;
}

Q019_NLB_INLINE std::uint64_t Mix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value ^= value >> 30;
  value *= 0xbf58476d1ce4e5b9ULL;
  value ^= value >> 27;
  value *= 0x94d049bb133111ebULL;
  value ^= value >> 31;
  return value;
}

Q019_NLB_INLINE std::uint64_t GlobalCellLinearId(
    const int global_i, const int global_j, const int global_k,
    const int root_nx1, const int root_nx2) {
  return (static_cast<std::uint64_t>(global_k)*static_cast<std::uint64_t>(root_nx2) +
          static_cast<std::uint64_t>(global_j))*static_cast<std::uint64_t>(root_nx1) +
      static_cast<std::uint64_t>(global_i);
}

Q019_NLB_INLINE std::uint64_t NestedPacketIdentity(
    const std::uint64_t global_cell_id, const int packet_within_cell,
    const int field_seed, const int particle_seed) {
  std::uint64_t value = Mix64(static_cast<std::uint64_t>(field_seed));
  value = Mix64(value ^ static_cast<std::uint64_t>(particle_seed));
  value = Mix64(value ^ global_cell_id);
  value = Mix64(value ^ static_cast<std::uint64_t>(packet_within_cell));
  return value;
}

Q019_NLB_INLINE double UniformOpen01(
    const std::uint64_t packet_identity, const int component) {
  const std::uint64_t bits = Mix64(
      packet_identity ^ static_cast<std::uint64_t>(component));
  return (static_cast<double>(bits >> 11) + 0.5)/
      9007199254740992.0;
}

Q019_NLB_INLINE Vector3 HaarRotatedOctahedralDirection(
    const std::uint64_t global_cell_id, const int packet_within_cell,
    const int direction, const int field_seed, const int particle_seed) {
  const std::uint64_t identity = NestedPacketIdentity(
      global_cell_id, packet_within_cell, field_seed, particle_seed);
  const double u1 = UniformOpen01(identity, 0);
  const double u2 = UniformOpen01(identity, 1);
  const double u3 = UniformOpen01(identity, 2);
  const double root_one_minus_u1 = std::sqrt(1.0 - u1);
  const double root_u1 = std::sqrt(u1);
  const double phase2 = 2.0*M_PI*u2;
  const double phase3 = 2.0*M_PI*u3;
  const double qx = root_one_minus_u1*std::sin(phase2);
  const double qy = root_one_minus_u1*std::cos(phase2);
  const double qz = root_u1*std::sin(phase3);
  const double qw = root_u1*std::cos(phase3);
  const Vector3 columns[3] = {
    {1.0 - 2.0*(qy*qy + qz*qz), 2.0*(qx*qy + qw*qz),
     2.0*(qx*qz - qw*qy)},
    {2.0*(qx*qy - qw*qz), 1.0 - 2.0*(qx*qx + qz*qz),
     2.0*(qy*qz + qw*qx)},
    {2.0*(qx*qz + qw*qy), 2.0*(qy*qz - qw*qx),
     1.0 - 2.0*(qx*qx + qy*qy)}
  };
  const int axis = direction/2;
  const double sign = (direction % 2 == 0) ? 1.0 : -1.0;
  return Scale(columns[axis], sign);
}

Q019_NLB_INLINE Vector3 NestedHaarOctahedralShellVelocityAtSample(
    const int ppc, const int sample, const std::uint64_t global_cell_id,
    const int field_seed, const int particle_seed,
    const double guide_parallel_drift, const double shell_speed) {
  const int direction = (sample % ppc) % 6;
  const int packet_within_cell = (sample % ppc)/6;
  const Vector3 shell = Scale(
      HaarRotatedOctahedralDirection(
          global_cell_id, packet_within_cell, direction, field_seed, particle_seed),
      shell_speed);
  return {guide_parallel_drift + shell.x1, shell.x2, shell.x3};
}

Q019_NLB_INLINE bool NoSubIonCellScaleEnvelopeIsValid(
    const double minimum_active_dx, const double background_ion_inertial_length) {
  return std::isfinite(minimum_active_dx) &&
      std::isfinite(background_ion_inertial_length) &&
      background_ion_inertial_length > 0.0 &&
      minimum_active_dx > background_ion_inertial_length;
}

Q019_NLB_INLINE double CharacteristicShellGyroradiusOverDx(
    const double shell_momentum_per_mass, const double q_over_mc,
    const double maximum_sampled_b_times_dx) {
  if (!(std::isfinite(shell_momentum_per_mass) && shell_momentum_per_mass > 0.0 &&
        std::isfinite(q_over_mc) && q_over_mc != 0.0 &&
        std::isfinite(maximum_sampled_b_times_dx) &&
        maximum_sampled_b_times_dx > 0.0)) {
    return 0.0;
  }
  return shell_momentum_per_mass/
      (std::fabs(q_over_mc)*maximum_sampled_b_times_dx);
}

Q019_NLB_INLINE Vector3 SharedSpectrumVectorPotentialAt(
    const SeedParameters parameters, const Vector3 position) {
  const double phase = parameters.k0*position.x1 + SeedPhase(parameters.field_seed, 0);
  Vector3 potential = {
    0.0,
    parameters.eigenmode_amplitude*std::sin(phase)/parameters.k0,
    parameters.eigenmode_amplitude*std::cos(phase)/parameters.k0
  };
  potential.x3 += parameters.broadband_amplitude/parameters.k0*
      (0.5*std::sin(parameters.k0*position.x1 + parameters.k0*position.x2 +
                    SeedPhase(parameters.field_seed, 1)) +
       0.35*std::sin(2.0*parameters.k0*position.x1 -
                     2.0*parameters.k0*position.x2 +
                     SeedPhase(parameters.field_seed, 2)));
  if (parameters.dimension == 3) {
    potential.x2 += 0.4*parameters.broadband_amplitude/parameters.k0*
        std::sin(parameters.k0*position.x1 + parameters.k0*position.x3 +
                 SeedPhase(parameters.field_seed, 3));
    potential.x1 += 0.3*parameters.broadband_amplitude/parameters.k0*
        std::sin(parameters.k0*position.x2 + parameters.k0*position.x3 +
                 SeedPhase(parameters.field_seed, 4));
  }
  return potential;
}

Q019_NLB_INLINE Vector3 SeparatedLongModeVectorPotentialAt(
    const SeedParameters parameters, const Vector3 position) {
  if (parameters.seed_topology !=
      SeedTopology::shared_spectrum_plus_separated_long_modes) {
    return {0.0, 0.0, 0.0};
  }
  const double k1 = 2.0*M_PI/parameters.x1_extent;
  const double k2 = 2.0*M_PI/parameters.x2_extent;
  const double k3 = 2.0*M_PI/parameters.x3_extent;
  Vector3 potential = {
    0.0,
    0.8*parameters.separated_long_mode_amplitude/k1*
        std::cos(k1*position.x1 + SeedPhase(parameters.field_seed, 6)),
    parameters.separated_long_mode_amplitude/k1*
        std::sin(k1*position.x1 + SeedPhase(parameters.field_seed, 5)) +
        0.6*parameters.separated_long_mode_amplitude/k2*
        std::cos(k2*position.x2 + SeedPhase(parameters.field_seed, 7))
  };
  if (parameters.dimension == 3) {
    potential.x1 += 0.7*parameters.separated_long_mode_amplitude/k3*
        std::sin(k3*position.x3 + SeedPhase(parameters.field_seed, 8));
    potential.x2 += 0.5*parameters.separated_long_mode_amplitude/k3*
        std::cos(k3*position.x3 + SeedPhase(parameters.field_seed, 9));
  }
  return potential;
}

Q019_NLB_INLINE Vector3 AxisAlignedPerturbationVectorPotentialAt(
    const SeedParameters parameters, const Vector3 position) {
  return Add(SharedSpectrumVectorPotentialAt(parameters, position),
             SeparatedLongModeVectorPotentialAt(parameters, position));
}

Q019_NLB_INLINE Vector3 AxisAlignedVectorPotentialAt(
    const SeedParameters parameters, const Vector3 position) {
  Vector3 potential = AxisAlignedPerturbationVectorPotentialAt(parameters, position);
  potential.x2 -= 0.5*parameters.b_g*position.x3;
  potential.x3 += 0.5*parameters.b_g*position.x2;
  return potential;
}

Q019_NLB_INLINE Vector3 AxisAlignedMagneticFieldAt(
    const SeedParameters parameters, const Vector3 position) {
  const double phase0 =
      parameters.k0*position.x1 + SeedPhase(parameters.field_seed, 0);
  const double phase1 = parameters.k0*(position.x1 + position.x2) +
      SeedPhase(parameters.field_seed, 1);
  const double phase2 = 2.0*parameters.k0*(position.x1 - position.x2) +
      SeedPhase(parameters.field_seed, 2);
  Vector3 field = {
    parameters.b_g + parameters.broadband_amplitude*
        (0.5*std::cos(phase1) - 0.7*std::cos(phase2)),
    parameters.eigenmode_amplitude*std::sin(phase0) -
        parameters.broadband_amplitude*
            (0.5*std::cos(phase1) + 0.7*std::cos(phase2)),
    parameters.eigenmode_amplitude*std::cos(phase0)
  };
  if (parameters.dimension == 3) {
    const double phase3 = parameters.k0*(position.x1 + position.x3) +
        SeedPhase(parameters.field_seed, 3);
    const double phase4 = parameters.k0*(position.x2 + position.x3) +
        SeedPhase(parameters.field_seed, 4);
    field.x1 -= 0.4*parameters.broadband_amplitude*std::cos(phase3);
    field.x2 += 0.3*parameters.broadband_amplitude*std::cos(phase4);
    field.x3 += 0.4*parameters.broadband_amplitude*std::cos(phase3) -
        0.3*parameters.broadband_amplitude*std::cos(phase4);
  }
  if (parameters.seed_topology ==
      SeedTopology::shared_spectrum_plus_separated_long_modes) {
    const double k1 = 2.0*M_PI/parameters.x1_extent;
    const double k2 = 2.0*M_PI/parameters.x2_extent;
    field.x1 -= 0.6*parameters.separated_long_mode_amplitude*
        std::sin(k2*position.x2 + SeedPhase(parameters.field_seed, 7));
    field.x2 -= parameters.separated_long_mode_amplitude*
        std::cos(k1*position.x1 + SeedPhase(parameters.field_seed, 5));
    field.x3 -= 0.8*parameters.separated_long_mode_amplitude*
        std::sin(k1*position.x1 + SeedPhase(parameters.field_seed, 6));
    if (parameters.dimension == 3) {
      const double k3 = 2.0*M_PI/parameters.x3_extent;
      field.x1 += 0.5*parameters.separated_long_mode_amplitude*
          std::sin(k3*position.x3 + SeedPhase(parameters.field_seed, 9));
      field.x2 += 0.7*parameters.separated_long_mode_amplitude*
          std::cos(k3*position.x3 + SeedPhase(parameters.field_seed, 8));
    }
  }
  return field;
}

Q019_NLB_INLINE double NominalRigidityLength(const double momentum_per_mass,
                                              const double omega) {
  return momentum_per_mass/omega;
}

Q019_NLB_INLINE Vector3 AxisAlignedEigenmodeVelocityAt(
    const SeedParameters parameters, const Vector3 position) {
  const double phase = parameters.k0*position.x1 + SeedPhase(parameters.field_seed, 0);
  const double cs = std::cos(phase);
  const double sn = std::sin(phase);
  const double growth = std::sqrt(1.0 - parameters.epsilon*parameters.epsilon);
  const double scale = parameters.eigenmode_amplitude/std::sqrt(parameters.rho);
  return {0.0, scale*(-parameters.epsilon*cs + growth*sn),
          scale*(growth*cs + parameters.epsilon*sn)};
}

}  // namespace q019_physics_first_nonlinear_bell_successor_v2

#undef Q019_NLB_INLINE

#endif  // PGEN_TESTS_Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HPP_
