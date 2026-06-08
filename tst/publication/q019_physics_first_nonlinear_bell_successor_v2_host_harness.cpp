//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#define Q019_PHYSICS_FIRST_NONLINEAR_BELL_SUCCESSOR_V2_HOST_CONTRACT 1
#include "../../src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace q019 = q019_physics_first_nonlinear_bell_successor_v2;

namespace {

struct SpectrumResult {
  double fraction;
  double fluctuation_sum;
  double unique_mode_power_sum;
};

bool SameMode(const q019::IntegerVector3 lhs, const q019::IntegerVector3 rhs) {
  return lhs.x1 == rhs.x1 && lhs.x2 == rhs.x2 && lhs.x3 == rhs.x3;
}

bool IsUniqueConjugateRepresentative(const q019::IntegerVector3 mode) {
  if (mode.x1 != 0) return mode.x1 > 0;
  if (mode.x2 != 0) return mode.x2 > 0;
  return mode.x3 > 0;
}

bool PhysicalLowKEnumerationIsComplete(
    const int dimension, const std::array<double, 3> &extents) {
  std::vector<q019::IntegerVector3> expected;
  for (int n1 = -6; n1 <= 6; ++n1) {
    for (int n2 = -6; n2 <= 6; ++n2) {
      for (int n3 = dimension == 3 ? -6 : 0; n3 <= (dimension == 3 ? 6 : 0); ++n3) {
        const q019::IntegerVector3 mode = {n1, n2, n3};
        if (IsUniqueConjugateRepresentative(mode) &&
            q019::BoxEdgePhysicalLowKModeIsSelected(
                mode, dimension, extents[0], extents[1], extents[2])) {
          expected.push_back(mode);
        }
      }
    }
  }
  if (expected.size() !=
      static_cast<std::size_t>(q019::BoxEdgeUniqueModeCount(dimension))) {
    return false;
  }
  for (int n = 0; n < q019::BoxEdgeUniqueModeCount(dimension); ++n) {
    const q019::IntegerVector3 actual = q019::BoxEdgeUniqueModeAt(n, dimension);
    if (std::count_if(
            expected.begin(), expected.end(),
            [actual](const q019::IntegerVector3 candidate) {
              return SameMode(actual, candidate);
            }) != 1) {
      return false;
    }
  }
  return true;
}

SpectrumResult EvaluatePureTransverseMode(
    const int dimension, const std::array<double, 3> &extents,
    const std::array<int, 3> &cells, const q019::IntegerVector3 pure_mode,
    const int x1_partitions) {
  const double two_pi = 2.0*std::acos(-1.0);
  const int total_cells = cells[0]*cells[1]*cells[2];
  double mean_b2 = 0.0;
  double mean_b3 = 0.0;
  for (int partition = 0; partition < x1_partitions; ++partition) {
    double partition_b2 = 0.0;
    double partition_b3 = 0.0;
    const int first_i = partition*cells[0]/x1_partitions;
    const int last_i = (partition + 1)*cells[0]/x1_partitions;
    for (int k = 0; k < cells[2]; ++k) {
      for (int j = 0; j < cells[1]; ++j) {
        for (int i = first_i; i < last_i; ++i) {
          const double x1 = (static_cast<double>(i) + 0.5)*extents[0]/cells[0];
          const double x2 = (static_cast<double>(j) + 0.5)*extents[1]/cells[1];
          const double x3 = (static_cast<double>(k) + 0.5)*extents[2]/cells[2];
          const double theta = two_pi*(
              pure_mode.x1*x1/extents[0] + pure_mode.x2*x2/extents[1] +
              (dimension == 3 ? pure_mode.x3*x3/extents[2] : 0.0));
          partition_b2 += std::cos(theta);
          partition_b3 += std::sin(theta);
        }
      }
    }
    mean_b2 += partition_b2;
    mean_b3 += partition_b3;
  }
  mean_b2 /= static_cast<double>(total_cells);
  mean_b3 /= static_cast<double>(total_cells);

  std::vector<std::array<double, 4>> amplitudes(
      q019::BoxEdgeUniqueModeCount(dimension), {0.0, 0.0, 0.0, 0.0});
  double fluctuation_sum = 0.0;
  for (int partition = 0; partition < x1_partitions; ++partition) {
    std::vector<std::array<double, 4>> partition_amplitudes(
        q019::BoxEdgeUniqueModeCount(dimension), {0.0, 0.0, 0.0, 0.0});
    double partition_fluctuation = 0.0;
    const int first_i = partition*cells[0]/x1_partitions;
    const int last_i = (partition + 1)*cells[0]/x1_partitions;
    for (int k = 0; k < cells[2]; ++k) {
      for (int j = 0; j < cells[1]; ++j) {
        for (int i = first_i; i < last_i; ++i) {
          const double x1 = (static_cast<double>(i) + 0.5)*extents[0]/cells[0];
          const double x2 = (static_cast<double>(j) + 0.5)*extents[1]/cells[1];
          const double x3 = (static_cast<double>(k) + 0.5)*extents[2]/cells[2];
          const double pure_theta = two_pi*(
              pure_mode.x1*x1/extents[0] + pure_mode.x2*x2/extents[1] +
              (dimension == 3 ? pure_mode.x3*x3/extents[2] : 0.0));
          const double db2 = std::cos(pure_theta) - mean_b2;
          const double db3 = std::sin(pure_theta) - mean_b3;
          partition_fluctuation += db2*db2 + db3*db3;
          const q019::ComplexValue x1_phase = {
            std::cos(two_pi*x1/extents[0]), -std::sin(two_pi*x1/extents[0])
          };
          const q019::ComplexValue x2_phase = {
            std::cos(two_pi*x2/extents[1]), -std::sin(two_pi*x2/extents[1])
          };
          const q019::ComplexValue x3_phase = dimension == 3 ?
              q019::ComplexValue{
                std::cos(two_pi*x3/extents[2]), -std::sin(two_pi*x3/extents[2])
              } : q019::ComplexValue{1.0, 0.0};
          for (int n = 0; n < q019::BoxEdgeUniqueModeCount(dimension); ++n) {
            const q019::ComplexValue phase = q019::BoxEdgeModePhase(
                x1_phase, x2_phase, x3_phase,
                q019::BoxEdgeUniqueModeAt(n, dimension));
            partition_amplitudes[n][0] += db2*phase.real;
            partition_amplitudes[n][1] += db2*phase.imag;
            partition_amplitudes[n][2] += db3*phase.real;
            partition_amplitudes[n][3] += db3*phase.imag;
          }
        }
      }
    }
    fluctuation_sum += partition_fluctuation;
    for (int n = 0; n < q019::BoxEdgeUniqueModeCount(dimension); ++n) {
      for (int component = 0; component < 4; ++component) {
        amplitudes[n][component] += partition_amplitudes[n][component];
      }
    }
  }
  double unique_mode_power_sum = 0.0;
  for (const auto &mode : amplitudes) {
    for (const double component : mode) {
      unique_mode_power_sum += component*component;
    }
  }
  return {
    q019::BoxEdgePowerFraction(
        static_cast<double>(total_cells), fluctuation_sum, unique_mode_power_sum),
    fluctuation_sum,
    unique_mode_power_sum
  };
}

}  // namespace

int main() {
  using q019::Add;
  using q019::BackgroundIonGyrofrequency;
  using q019::BackgroundIonInertialLength;
  using q019::BaiHallGrowthRateReductionFactor;
  using q019::BaiHallGrowthRateFractionalShift;
  using q019::BaiHallLinearFactor;
  using q019::BaiHallWavenumberReductionFactor;
  using q019::BaiHallWavenumberFractionalShift;
  using q019::BoxEdgeCrossedSlotCount;
  using q019::BoxEdgeDiagnosticStatus;
  using q019::BoxEdgeFirstCrossingChronologyIsValid;
  using q019::BoxEdgeLastNominalTime;
  using q019::BoxEdgeNextNominalTime;
  using q019::BoxEdgePhysicalLowKModeIsSelected;
  using q019::BoxEdgePowerFraction;
  using q019::BoxEdgeStatusBit;
  using q019::BoxEdgeStatusRecordsValidMetric;
  using q019::BoxEdgeUniqueModeCount;
  using q019::BoxEdgeUniqueModeAt;
  using q019::CharacteristicShellGyroradiusOverDx;
  using q019::ChargeDensityRatio;
  using q019::DepositedJOverC;
  using q019::GlobalCellLinearId;
  using q019::HallParameter;
  using q019::HallParameterFromCurrent;
  using q019::NestedHaarOctahedralShellVelocityAtSample;
  using q019::NoSubIonCellScaleEnvelopeIsValid;
  using q019::OctahedralShellPacketLayoutIsValid;
  using q019::ResolutionEnvelopeIsValid;
  using q019::RootCellVolume;
  using q019::Vector3;

  const double pi = std::acos(-1.0);
  const double k0 = 2.0*pi;
  const double qom = 10000.0;
  const double rho_cr = 3.0e-6;
  const double drift = 2.0*k0/(rho_cr*qom);
  const double shell = 4.0*qom/k0;
  const double volume = RootCellVolume(8.0, 1024, 4.0, 512, 1.0, 1);
  const double qscale = rho_cr*volume/48.0;
  const double omega_i = BackgroundIonGyrofrequency(qom, 1.0);
  const double di = BackgroundIonInertialLength(1.0, omega_i);
  const double exact_r = ChargeDensityRatio(rho_cr, qom, qom);
  const double lambda = HallParameterFromCurrent(2.0*k0, 1.0, qom, 1.0, exact_r);
  std::cout << std::setprecision(17);
  std::cout << "mapping " << di << ' ' << k0*di << ' ' << lambda << ' '
            << BaiHallLinearFactor(lambda) << ' '
            << BaiHallGrowthRateFractionalShift(lambda) << ' '
            << BaiHallWavenumberFractionalShift(lambda) << ' '
            << exact_r << ' '
            << HallParameter(exact_r, drift, 1.0) << ' '
            << BaiHallGrowthRateReductionFactor(lambda) << ' '
            << BaiHallWavenumberReductionFactor(lambda) << '\n';
  std::cout << "current "
            << DepositedJOverC(48.0, qscale, qom, drift, volume) << '\n';
  std::cout << "scale "
            << NoSubIonCellScaleEnvelopeIsValid(8.0/1536.0, di) << ' '
            << NoSubIonCellScaleEnvelopeIsValid(0.5*di, di) << ' '
            << (8.0/1536.0)/di << '\n';
  std::cout << "packet "
            << OctahedralShellPacketLayoutIsValid(48, 1024, 4, 48*1024*4) << ' '
            << OctahedralShellPacketLayoutIsValid(32, 1024, 4, 32*1024*4) << '\n';
  std::cout << "resolution "
            << ResolutionEnvelopeIsValid(10.2, 8.0) << ' '
            << ResolutionEnvelopeIsValid(7.9, 8.0) << ' '
            << CharacteristicShellGyroradiusOverDx(shell, qom, 1.0/16.0) << '\n';
  std::cout << "box_edge_selection "
            << BoxEdgeUniqueModeCount(2) << ' ' << BoxEdgeUniqueModeCount(3) << ' '
            << BoxEdgePhysicalLowKModeIsSelected({2, 0, 0}, 2, 8.0, 4.0, 1.0)
            << ' '
            << BoxEdgePhysicalLowKModeIsSelected({3, 0, 0}, 2, 8.0, 4.0, 1.0)
            << ' '
            << BoxEdgePhysicalLowKModeIsSelected({3, 0, 0}, 3, 16.0, 8.0, 8.0)
            << ' '
            << BoxEdgePhysicalLowKModeIsSelected({4, 0, 0}, 3, 16.0, 8.0, 8.0)
            << ' ' << PhysicalLowKEnumerationIsComplete(2, {8.0, 4.0, 1.0})
            << ' ' << PhysicalLowKEnumerationIsComplete(3, {16.0, 8.0, 8.0})
            << '\n';
  const SpectrumResult selected_2d = EvaluatePureTransverseMode(
      2, {8.0, 4.0, 1.0}, {32, 16, 1}, {2, 0, 0}, 1);
  const SpectrumResult unselected_2d = EvaluatePureTransverseMode(
      2, {8.0, 4.0, 1.0}, {32, 16, 1}, {3, 0, 0}, 1);
  const SpectrumResult selected_3d = EvaluatePureTransverseMode(
      3, {16.0, 8.0, 8.0}, {16, 8, 8}, {2, 1, -1}, 1);
  const SpectrumResult selected_3d_partitioned = EvaluatePureTransverseMode(
      3, {16.0, 8.0, 8.0}, {16, 8, 8}, {2, 1, -1}, 4);
  const SpectrumResult unselected_3d = EvaluatePureTransverseMode(
      3, {16.0, 8.0, 8.0}, {16, 8, 8}, {4, 0, 0}, 1);
  std::cout << "box_edge_numeric "
            << selected_2d.fraction << ' ' << unselected_2d.fraction << ' '
            << selected_3d.fraction << ' ' << unselected_3d.fraction << ' '
            << std::abs(selected_3d.fraction - selected_3d_partitioned.fraction)
            << ' ' << selected_3d.fluctuation_sum << ' '
            << selected_3d.unique_mode_power_sum << '\n';
  std::cout << "box_edge_schedule "
            << BoxEdgePowerFraction(8.0, 0.0, 16.0) << ' '
            << BoxEdgeCrossedSlotCount(1, 0.1, 0.35) << ' '
            << BoxEdgeFirstCrossingChronologyIsValid(0.09, 0.1, 0.11) << ' '
            << BoxEdgeFirstCrossingChronologyIsValid(0.1, 0.1, 0.11) << ' '
            << BoxEdgeLastNominalTime(4, 0.1) << ' '
            << BoxEdgeNextNominalTime(4, 0.1) << ' '
            << BoxEdgeStatusRecordsValidMetric(
                   static_cast<int>(BoxEdgeDiagnosticStatus::valid)) << ' '
            << BoxEdgeStatusRecordsValidMetric(
                   static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation)) << ' '
            << BoxEdgeStatusBit(
                   static_cast<int>(BoxEdgeDiagnosticStatus::cadence_skipped)) << ' '
            << BoxEdgeStatusBit(
                   static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation)) << '\n';
  std::cout << "box_edge_exact_cadence "
            << BoxEdgeCrossedSlotCount(0, 0.1, 0.1) << ' '
            << BoxEdgeCrossedSlotCount(1, 0.1, 0.2) << ' '
            << BoxEdgeCrossedSlotCount(2, 0.1, 0.3) << ' '
            << BoxEdgeCrossedSlotCount(3, 0.1, 0.4) << ' '
            << BoxEdgeCrossedSlotCount(0, 0.1, 0.099) << ' '
            << BoxEdgeLastNominalTime(4, 0.1) << ' '
            << BoxEdgeNextNominalTime(4, 0.1) << '\n';

  const std::uint64_t global_cell = GlobalCellLinearId(23, 11, 0, 64, 32);
  bool nested = true;
  for (int sample = 0; sample < 24; ++sample) {
    const Vector3 v24 = NestedHaarOctahedralShellVelocityAtSample(
        24, sample, global_cell, 19001, 29001, drift, shell);
    const Vector3 v48 = NestedHaarOctahedralShellVelocityAtSample(
        48, sample, global_cell, 19001, 29001, drift, shell);
    const Vector3 v96 = NestedHaarOctahedralShellVelocityAtSample(
        96, sample, global_cell, 19001, 29001, drift, shell);
    nested = nested && v24.x1 == v48.x1 && v24.x2 == v48.x2 && v24.x3 == v48.x3 &&
        v24.x1 == v96.x1 && v24.x2 == v96.x2 && v24.x3 == v96.x3;
  }
  const std::uint64_t global_cell_layout_a =
      GlobalCellLinearId(1*16 + 7, 1*8 + 3, 0, 64, 32);
  const std::uint64_t global_cell_layout_b =
      GlobalCellLinearId(2*8 + 7, 2*4 + 3, 0, 64, 32);
  const Vector3 decomposition_a = NestedHaarOctahedralShellVelocityAtSample(
      48, 31, global_cell_layout_a, 19001, 29001, drift, shell);
  const Vector3 decomposition_b = NestedHaarOctahedralShellVelocityAtSample(
      48, 31, global_cell_layout_b, 19001, 29001, drift, shell);
  const double decomposition_difference = std::max(
      {std::abs(decomposition_a.x1 - decomposition_b.x1),
       std::abs(decomposition_a.x2 - decomposition_b.x2),
       std::abs(decomposition_a.x3 - decomposition_b.x3)});
  std::cout << "nested " << nested << '\n';
  std::cout << "decomposition " << (global_cell_layout_a == global_cell_layout_b) << ' '
            << decomposition_difference << ' ' << global_cell_layout_a << ' '
            << global_cell_layout_b << '\n';

  for (const int ppc : {24, 48, 96}) {
    Vector3 sum = {0.0, 0.0, 0.0};
    double second[3][3] = {{0.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0}};
    double maximum_packet_mean_error = 0.0;
    double maximum_packet_second_error = 0.0;
    for (int sample = 0; sample < ppc; ++sample) {
      const Vector3 velocity =
          NestedHaarOctahedralShellVelocityAtSample(
              ppc, sample, global_cell, 19001, 29001, drift, shell);
      sum = Add(sum, velocity);
      const double centered[3] = {velocity.x1 - drift, velocity.x2, velocity.x3};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          second[i][j] += centered[i]*centered[j]/static_cast<double>(ppc);
        }
      }
    }
    for (int packet = 0; packet < ppc/6; ++packet) {
      double packet_sum[3] = {0.0, 0.0, 0.0};
      double packet_second[3][3] = {{0.0, 0.0, 0.0},
                                     {0.0, 0.0, 0.0},
                                     {0.0, 0.0, 0.0}};
      for (int direction = 0; direction < 6; ++direction) {
        const Vector3 velocity = NestedHaarOctahedralShellVelocityAtSample(
            ppc, packet*6 + direction, global_cell, 19001, 29001, drift, shell);
        const double centered[3] = {velocity.x1 - drift, velocity.x2, velocity.x3};
        for (int i = 0; i < 3; ++i) {
          packet_sum[i] += centered[i]/6.0;
          for (int j = 0; j < 3; ++j) {
            packet_second[i][j] += centered[i]*centered[j]/6.0;
          }
        }
      }
      for (int i = 0; i < 3; ++i) {
        maximum_packet_mean_error =
            std::max(maximum_packet_mean_error, std::abs(packet_sum[i]));
        for (int j = 0; j < 3; ++j) {
          const double expected = i == j ? shell*shell/3.0 : 0.0;
          maximum_packet_second_error = std::max(
              maximum_packet_second_error, std::abs(packet_second[i][j] - expected));
        }
      }
    }
    std::cout << "shell " << ppc << ' '
              << sum.x1/static_cast<double>(ppc) - drift << ' '
              << sum.x2/static_cast<double>(ppc) << ' '
              << sum.x3/static_cast<double>(ppc);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double expected = i == j ? shell*shell/3.0 : 0.0;
        std::cout << ' ' << second[i][j] - expected;
      }
    }
    std::cout << '\n';
    std::cout << "packet_moments " << ppc << ' ' << maximum_packet_mean_error << ' '
              << maximum_packet_second_error << '\n';

    double fourth[3][3][3][3] = {};
    constexpr int cells = 4096;
    for (int cell = 0; cell < cells; ++cell) {
      for (int sample = 0; sample < ppc; ++sample) {
        const Vector3 velocity = NestedHaarOctahedralShellVelocityAtSample(
            ppc, sample, static_cast<std::uint64_t>(cell), 19001, 29001, 0.0, 1.0);
        const double value[3] = {velocity.x1, velocity.x2, velocity.x3};
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
              for (int l = 0; l < 3; ++l) {
                fourth[i][j][k][l] += value[i]*value[j]*value[k]*value[l]/
                    static_cast<double>(cells*ppc);
              }
            }
          }
        }
      }
    }
    double squared_error = 0.0;
    double maximum_error = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          for (int l = 0; l < 3; ++l) {
            const double expected =
                ((i == j && k == l) + (i == k && j == l) + (i == l && j == k))/15.0;
            const double error = fourth[i][j][k][l] - expected;
            squared_error += error*error;
            maximum_error = std::max(maximum_error, std::abs(error));
          }
        }
      }
    }
    std::cout << "fourth " << ppc << ' ' << std::sqrt(squared_error/81.0) << ' '
              << maximum_error << '\n';
  }
}
