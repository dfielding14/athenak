//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file q019_physics_first_nonlinear_bell_successor_v2.cpp
//! \brief Corrected source-local nonlinear Bell production successor.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "outputs/restart_utils.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"
#include "q019_physics_first_nonlinear_bell_successor_v2.hpp"

namespace {

constexpr const char *kQ019RuntimeControllerBlock =
    "q019_nonlinear_bell_runtime_controller_v1";
constexpr int kQ019RuntimeControllerSchema = 1;
constexpr int kQ019ResolutionStopReason = 1901;
constexpr int kQ019BoxEdgeStopReason = 1902;
constexpr int kQ019PilotCycleStopReason = 1903;
constexpr int kQ019DiagnosticFailureReason = 1991;
constexpr std::array<const char *, 15> kQ019RuntimeControllerImmutableParameters = {
  "schema",
  "authority",
  "contract_id",
  "cadence_status",
  "stop_disposition",
  "source_case_id",
  "monitor_dt",
  "box_edge_monitor_enabled",
  "resolution_monitor_enabled",
  "diagnostic_failure_stop_armed",
  "resolution_stop_armed",
  "resolution_stop_B_over_B0",
  "box_edge_stop_armed",
  "box_edge_stop_ppm",
  "pilot_cycle_limit"
};
constexpr std::array<const char *, 11> kQ019RuntimeControllerMutableParameters = {
  "runtime_resolution_samples",
  "runtime_resolution_last_cycle",
  "runtime_resolution_last_time",
  "runtime_resolution_last_B_over_B0",
  "runtime_resolution_max_B_over_B0",
  "runtime_controller_triggered",
  "runtime_controller_trigger_failure",
  "runtime_controller_trigger_reason",
  "runtime_controller_trigger_cycle",
  "runtime_controller_trigger_time",
  "runtime_controller_trigger_metric"
};

using q019_physics_first_nonlinear_bell_successor_v2::Add;
using q019_physics_first_nonlinear_bell_successor_v2::
    AxisAlignedEigenmodeVelocityAt;
using q019_physics_first_nonlinear_bell_successor_v2::
    AxisAlignedMagneticFieldAt;
using q019_physics_first_nonlinear_bell_successor_v2::
    AxisAlignedVectorPotentialAt;
using q019_physics_first_nonlinear_bell_successor_v2::BackgroundIonGyrofrequency;
using q019_physics_first_nonlinear_bell_successor_v2::BackgroundIonInertialLength;
using q019_physics_first_nonlinear_bell_successor_v2::BaiHallGrowthRateReductionFactor;
using q019_physics_first_nonlinear_bell_successor_v2::BaiHallGrowthRateFractionalShift;
using q019_physics_first_nonlinear_bell_successor_v2::BaiHallLinearFactor;
using q019_physics_first_nonlinear_bell_successor_v2::BaiHallWavenumberReductionFactor;
using q019_physics_first_nonlinear_bell_successor_v2::BaiHallWavenumberFractionalShift;
using q019_physics_first_nonlinear_bell_successor_v2::Branch;
using q019_physics_first_nonlinear_bell_successor_v2::BranchRigidityIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::BranchSamplingIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeCrossedSlotCount;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeDiagnosticStatus;
using q019_physics_first_nonlinear_bell_successor_v2::
    BoxEdgeFirstCrossingChronologyIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeFrozenAspectRatioIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeModePhase;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeLastNominalTime;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeNextNominalTime;
using q019_physics_first_nonlinear_bell_successor_v2::
    BoxEdgePhysicalLowKModeIsSelected;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgePowerFraction;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeReductionGroupCount;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeStatusBit;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeStatusCodeIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeStatusRecordsValidMetric;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeUniqueModeCount;
using q019_physics_first_nonlinear_bell_successor_v2::BoxEdgeUniqueModeAt;
using q019_physics_first_nonlinear_bell_successor_v2::ChargeDensityRatio;
using q019_physics_first_nonlinear_bell_successor_v2::ComplexValue;
using q019_physics_first_nonlinear_bell_successor_v2::CRInertiaParameter;
using q019_physics_first_nonlinear_bell_successor_v2::
    CRRMSSpeedKineticLoadingProxyToBackgroundMagneticEnergy;
using q019_physics_first_nonlinear_bell_successor_v2::CRMassDensityOverRho0;
using q019_physics_first_nonlinear_bell_successor_v2::CRMomentumLoadingParameter;
using q019_physics_first_nonlinear_bell_successor_v2::CRNumberDensity;
using q019_physics_first_nonlinear_bell_successor_v2::DepositedJOverC;
using q019_physics_first_nonlinear_bell_successor_v2::DepositedSpeciesSumJOverC;
using q019_physics_first_nonlinear_bell_successor_v2::FeedbackForceParameter;
using q019_physics_first_nonlinear_bell_successor_v2::FiniteSamplingMode;
using q019_physics_first_nonlinear_bell_successor_v2::FiniteSamplingModeIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::GlobalCellLinearId;
using q019_physics_first_nonlinear_bell_successor_v2::HallParameter;
using q019_physics_first_nonlinear_bell_successor_v2::HallParameterFromCurrent;
using q019_physics_first_nonlinear_bell_successor_v2::IntegerVector3;
using q019_physics_first_nonlinear_bell_successor_v2::
    LoadingAccountingIsFiniteAndPositive;
using q019_physics_first_nonlinear_bell_successor_v2::
    NoSubIonCellScaleEnvelopeIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::NominalRigidityLength;
using q019_physics_first_nonlinear_bell_successor_v2::OctahedralShellPacketLayoutIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::RequiredDepositedJOverC;
using q019_physics_first_nonlinear_bell_successor_v2::ResolutionEnvelopeIsValid;
using q019_physics_first_nonlinear_bell_successor_v2::RootCellVolume;
using q019_physics_first_nonlinear_bell_successor_v2::
    NestedHaarOctahedralShellVelocityAtSample;
using q019_physics_first_nonlinear_bell_successor_v2::SeedParameters;
using q019_physics_first_nonlinear_bell_successor_v2::SeedTopology;
using q019_physics_first_nonlinear_bell_successor_v2::Vector3;
using q019_physics_first_nonlinear_bell_successor_v2::kBoxEdgeModesPerReduction;
using q019_physics_first_nonlinear_bell_successor_v2::kBoxEdgeMaximumUniqueModeCount;

[[noreturn]] void Q019NonlinearFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  restart_utils::AbortOnFatalError();
}

void Q019RequireClose(const std::string &label, const Real measured,
                      const Real expected) {
  if (!std::isfinite(measured) || !std::isfinite(expected)) {
    Q019NonlinearFatal(label + " must be finite");
  }
  const Real scale = std::max(static_cast<Real>(1.0), std::abs(expected));
  if (std::abs(measured - expected) > static_cast<Real>(1.0e-12)*scale) {
    Q019NonlinearFatal(label + " does not match the corrected nonlinear Bell contract");
  }
}

void Q019RequireMomentClose(const std::string &label, const Real measured,
                            const Real expected, const Real physical_scale) {
  if (!std::isfinite(measured) || !std::isfinite(expected) ||
      !std::isfinite(physical_scale) || physical_scale < 0.0) {
    Q019NonlinearFatal(label + " moment comparison must be finite");
  }
  const Real scale = std::max(
      {static_cast<Real>(1.0), std::abs(expected), physical_scale});
  const Real tolerance = static_cast<Real>(512.0)*
      std::numeric_limits<Real>::epsilon()*scale;
  if (std::abs(measured - expected) > tolerance) {
    std::ostringstream message;
    message << label << " does not match the nested isotropic-shell contract: measured="
            << measured << " expected=" << expected << " tolerance=" << tolerance;
    Q019NonlinearFatal(message.str());
  }
}

void Q019RequireString(ParameterInput *pin, const std::string &block,
                       const std::string &name, const std::string &expected) {
  if (pin->GetString(block, name).compare(expected) != 0) {
    Q019NonlinearFatal("<" + block + ">/" + name +
                       " does not match the corrected nonlinear Bell contract");
  }
}

void Q019RequireBoolean(ParameterInput *pin, const std::string &block,
                        const std::string &name, const bool expected) {
  if (pin->GetBoolean(block, name) != expected) {
    Q019NonlinearFatal("<" + block + ">/" + name +
                       " does not match the corrected nonlinear Bell contract");
  }
}

bool Q019MutableRuntimeBookkeepingParameter(
    const std::string &block_name, const std::string &parameter_name,
    const std::string &q019_block) {
  if (block_name.rfind("output", 0) == 0 &&
      (parameter_name == "file_number" || parameter_name == "last_time")) {
    return true;
  }
  if (block_name != q019_block) return false;
  return parameter_name == "runtime_box_edge_monitor_next_nominal_time" ||
      parameter_name == "runtime_box_edge_monitor_last_nominal_time" ||
      parameter_name == "runtime_box_edge_monitor_last_cycle" ||
      parameter_name == "runtime_box_edge_monitor_completed_slots" ||
      parameter_name == "runtime_box_edge_monitor_valid_samples" ||
      parameter_name == "runtime_box_edge_monitor_skipped_slots" ||
      parameter_name == "runtime_box_edge_monitor_last_status" ||
      parameter_name == "runtime_box_edge_monitor_status_mask" ||
      parameter_name == "runtime_box_edge_monitor_last_prior_time" ||
      parameter_name == "runtime_box_edge_monitor_last_time" ||
      parameter_name == "runtime_box_edge_monitor_last_power_fraction" ||
      parameter_name == "runtime_box_edge_monitor_max_power_fraction" ||
      parameter_name == "runtime_box_edge_monitor_last_fluctuation_mean";
}

std::string Q019MatrixIdentityPayload(ParameterInput *pin,
                                      const std::string &q019_block) {
  std::vector<std::tuple<std::string, std::string, std::string>> entries;
  for (const auto &input_block : pin->block) {
    if (input_block.block_name == "comment" ||
        input_block.block_name == kQ019RuntimeControllerBlock) continue;
    for (const auto &input_line : input_block.line) {
      if (input_block.block_name == q019_block &&
          input_line.param_name == "matrix_identity_fingerprint") {
        continue;
      }
      if (Q019MutableRuntimeBookkeepingParameter(
              input_block.block_name, input_line.param_name, q019_block)) {
        continue;
      }
      entries.emplace_back(input_block.block_name, input_line.param_name,
                           input_line.param_value);
    }
  }
  std::sort(entries.begin(), entries.end());
  std::ostringstream payload;
  for (const auto &[block_name, parameter_name, value] : entries) {
    payload << block_name << '/' << parameter_name << '=' << value << '\n';
  }
  return payload.str();
}

std::uint32_t Q019RotateRight(const std::uint32_t value, const int shift) {
  return (value >> shift) | (value << (32 - shift));
}

std::string Q019MatrixIdentityFingerprint(const std::string &payload) {
  constexpr std::array<std::uint32_t, 64> round_constants = {
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
    0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
    0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
    0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
    0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
    0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
    0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
    0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
    0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
    0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
    0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
    0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
    0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
    0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
    0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
    0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U
  };
  std::array<std::uint32_t, 8> state = {
    0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
    0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U
  };
  std::vector<unsigned char> bytes(payload.begin(), payload.end());
  const std::uint64_t payload_bits = static_cast<std::uint64_t>(bytes.size())*8U;
  bytes.push_back(0x80U);
  while (bytes.size() % 64U != 56U) bytes.push_back(0U);
  for (int shift = 56; shift >= 0; shift -= 8) {
    bytes.push_back(static_cast<unsigned char>((payload_bits >> shift) & 0xffU));
  }
  for (std::size_t offset = 0; offset < bytes.size(); offset += 64U) {
    std::array<std::uint32_t, 64> words = {};
    for (int n = 0; n < 16; ++n) {
      const std::size_t index = offset + static_cast<std::size_t>(4*n);
      words[n] = (static_cast<std::uint32_t>(bytes[index]) << 24) |
          (static_cast<std::uint32_t>(bytes[index + 1]) << 16) |
          (static_cast<std::uint32_t>(bytes[index + 2]) << 8) |
          static_cast<std::uint32_t>(bytes[index + 3]);
    }
    for (int n = 16; n < 64; ++n) {
      const std::uint32_t s0 = Q019RotateRight(words[n - 15], 7) ^
          Q019RotateRight(words[n - 15], 18) ^ (words[n - 15] >> 3);
      const std::uint32_t s1 = Q019RotateRight(words[n - 2], 17) ^
          Q019RotateRight(words[n - 2], 19) ^ (words[n - 2] >> 10);
      words[n] = words[n - 16] + s0 + words[n - 7] + s1;
    }
    std::uint32_t a = state[0];
    std::uint32_t b = state[1];
    std::uint32_t c = state[2];
    std::uint32_t d = state[3];
    std::uint32_t e = state[4];
    std::uint32_t f = state[5];
    std::uint32_t g = state[6];
    std::uint32_t h = state[7];
    for (int n = 0; n < 64; ++n) {
      const std::uint32_t sum1 = Q019RotateRight(e, 6) ^
          Q019RotateRight(e, 11) ^ Q019RotateRight(e, 25);
      const std::uint32_t choose = (e & f) ^ ((~e) & g);
      const std::uint32_t temporary1 =
          h + sum1 + choose + round_constants[n] + words[n];
      const std::uint32_t sum0 = Q019RotateRight(a, 2) ^
          Q019RotateRight(a, 13) ^ Q019RotateRight(a, 22);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temporary2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temporary1;
      d = c;
      c = b;
      b = a;
      a = temporary1 + temporary2;
    }
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
  }
  std::ostringstream fingerprint;
  fingerprint << std::hex << std::setfill('0');
  for (const std::uint32_t word : state) fingerprint << std::setw(8) << word;
  return fingerprint.str();
}

std::string Q019RuntimeControllerIdentityPayload(ParameterInput *pin) {
  std::vector<std::pair<std::string, std::string>> entries;
  std::array<bool, kQ019RuntimeControllerImmutableParameters.size()>
      immutable_seen = {};
  std::array<bool, kQ019RuntimeControllerMutableParameters.size()>
      mutable_seen = {};
  int fingerprint_count = 0;
  for (const auto &input_block : pin->block) {
    if (input_block.block_name != kQ019RuntimeControllerBlock) continue;
    for (const auto &input_line : input_block.line) {
      if (input_line.param_name == "controller_identity_fingerprint") {
        ++fingerprint_count;
        continue;
      }
      const auto immutable = std::find(
          kQ019RuntimeControllerImmutableParameters.begin(),
          kQ019RuntimeControllerImmutableParameters.end(),
          input_line.param_name);
      if (immutable != kQ019RuntimeControllerImmutableParameters.end()) {
        const auto index = static_cast<std::size_t>(
            immutable - kQ019RuntimeControllerImmutableParameters.begin());
        if (immutable_seen[index]) {
          Q019NonlinearFatal("Q019 runtime controller parameter is duplicated");
        }
        immutable_seen[index] = true;
        entries.emplace_back(input_line.param_name, input_line.param_value);
        continue;
      }
      const auto mutable_parameter = std::find(
          kQ019RuntimeControllerMutableParameters.begin(),
          kQ019RuntimeControllerMutableParameters.end(),
          input_line.param_name);
      if (mutable_parameter == kQ019RuntimeControllerMutableParameters.end()) {
        Q019NonlinearFatal("Q019 runtime controller contains an unknown parameter");
      }
      const auto index = static_cast<std::size_t>(
          mutable_parameter - kQ019RuntimeControllerMutableParameters.begin());
      if (mutable_seen[index]) {
        Q019NonlinearFatal("Q019 runtime controller parameter is duplicated");
      }
      mutable_seen[index] = true;
    }
  }
  if (fingerprint_count != 1 ||
      std::find(immutable_seen.begin(), immutable_seen.end(), false) !=
          immutable_seen.end()) {
    Q019NonlinearFatal("Q019 runtime controller immutable inventory drifted");
  }
  const int mutable_count = static_cast<int>(
      std::count(mutable_seen.begin(), mutable_seen.end(), true));
  if (mutable_count != 0 &&
      mutable_count != static_cast<int>(mutable_seen.size())) {
    Q019NonlinearFatal("Q019 runtime controller mutable inventory drifted");
  }
  std::sort(entries.begin(), entries.end());
  std::ostringstream payload;
  for (const auto &[parameter_name, value] : entries) {
    payload << kQ019RuntimeControllerBlock << '/' << parameter_name
            << '=' << value << '\n';
  }
  return payload.str();
}

struct Q019CanonicalMatrixIdentity {
  const char *case_id;
  const char *fingerprint;
};

constexpr Q019CanonicalMatrixIdentity q019_canonical_matrix_identities[] = {
  {"q019-q023-carrier-s1-onset-s0",
   "e0028fddc9f769bc4096773196528be1df3a35083386fc148be47640f4a71a4b"},
  {"q019-q023-carrier-s1-onset-s1",
   "32c355d803a53a999e28c0957dfbfbbf7eb9b0d8d4416ef7175ae9e49ffc82b1"},
  {"q019-q023-carrier-s2-window-fiducial-s0",
   "97a0b9c85a2e603a1dd92a3857d390c6cc070d6ba6c66a92d2619fe961b09629"},
  {"q019-q023-carrier-s2-window-fiducial-s1",
   "bf2c7a617ef5291fa9b4d9ec4c8bb29beec7459f48f53291f67114deb384605a"},
  {"q019-q023-carrier-s2-resolution-coarse-s0",
   "47eb08f359e93c8b6927069d7b87cf1635a216c503594792b50d229d8a4d2f52"},
  {"q019-q023-carrier-s2-resolution-fine-s0",
   "f580ffd28b6fbc96e5d886ae9ced50f3c2b73d290caf9ef86b994fba29eba0a8"},
  {"q019-q023-carrier-s2-cfl-small-s0",
   "5d5a1bc259debe274e00c8cfc3a514956359174b011d961b81047c44490ab710"},
  {"q019-q023-carrier-s2-cell-cross-one-s0",
   "7b7baa9429f813f8c5f1f2f15821969d2ff8a691616486cf894ffa3e953e2cde"},
  {"q019-q023-carrier-s2-riemann-hlld-s0",
   "e79a061ca125612a63ffcaf511e0cff7b33952dadd825198366655e2acff418c"},
  {"q019-q023-carrier-s2-reconstruct-wenoz-s0",
   "ca916cf8cd79059473acac2835f7f29b794b1db526ac11410ac3505efe9f9f59"},
  {"q019-q023-carrier-s2-qom1em4-s0",
   "03f66330d3ab1da4f06246a8582644553e9474bd5893fc24a0980b6db062f6bc"},
  {"q019-q023-carrier-s2-qom1em5-s0",
   "f6927649d2efe6281bd74003f10de86da3a343b7c90c443b698f50b96a1b2f14"},
  {"q019-q023-carrier-s3-3d-fiducial-s0",
   "023503a56d13b4f9f432435ac5e9e90637428a1282703237d5f1394c46ba0be7"},
  {"q019-q023-carrier-s3-3d-large-s0",
   "53c4b368e999a85aecc32c8eab356726c8dd602ad236951f6c02adb262bef4d7"},
  {"q019-hr-current-retention-s0",
   "799cb9c639ac432481b97b74b807283a9369ae94f91d020850d47561058c527e"},
  {"q019-hr-current-retention-s1",
   "db20573806afee27d9017aacf0def2ee91ad215d83172647a3beb1320eb16105"},
  {"q019-hr-current-retention-s2",
   "1117ae78944dc5eec328aa25c98a88bb45758ea7ea5380afc0232d206737a5ea"},
  {"q019-hr-fiducial-ppc8-centered-s0",
   "c7d312be0a20ee31688d0ba112d3e1cc48deaa6670f0939a9b85ffff359e5d3e"},
  {"q019-hr-fiducial-ppc32-centered-s0",
   "c990d8b56c6538988426039c36971a5a2a573f6e0bf73b8d1bf8a7f5b784f741"},
  {"q019-hr-fiducial-resolution-coarse-s0",
   "3ccdec0d0d207e9fb56b4696f1ec0c53b8d81ca9906b26b87606c747bcdb5aab"},
  {"q019-hr-fiducial-resolution-fine-s0",
   "e1f860e805a5ee6e4e817c61419f3c8305ea11bbdb3823ce527ebdec686f2bd8"},
  {"q019-hr-fiducial-particle-step-small-s0",
   "69b3aecc3e71914a122069e844959ad00c7453bd5b250232e83dfadb49533b5c"},
  {"q019-hr-fiducial-stochastic-position-s0",
   "a1cd6d1aed0844e7eb29cb11b3a742c315c14bf9094e37b658d57ec0ca0ac5bf"},
  {"q019-hr-fiducial-stochastic-position-s1",
   "78b5341ab6987262d8587d78446f24f2993f83f0a03fa78fe3c03b82ffb31f11"},
  {"q019-hr-fiducial-stochastic-position-s2",
   "b111215d22dc5a9bce3c5ed008d35685908d01a53731b81e77eda52f42d62edf"},
  {"q019-fr-grid-k4-rho1em06-s0",
   "d6292f5792ac29a8428a0c5b0a49eccf9d54eef454b006e6514f1e8b5a7b8a9e"},
  {"q019-fr-grid-k4-rho1em06-s1",
   "ed35d6e4d903665648805266f1c47f4cb5aec742f3c6312e5fe1b2091cd83ddd"},
  {"q019-fr-grid-k4-rho1em06-s2",
   "b44f225520174d9ee35c21933ecbd96a2f1d0d06969172fb9429eb22577ff709"},
  {"q019-fr-grid-k4-rho3em06-s0",
   "6b01ad9d66d36391e43e135988bc309e11247d5d9941682762fbd27745e67254"},
  {"q019-fr-grid-k4-rho3em06-s1",
   "0dc6892c65529a4670a15e20371b2532c374f3cb804b5ccce53deaea767c469f"},
  {"q019-fr-grid-k4-rho3em06-s2",
   "12aef7a88a2cadc410b345f5398a05d6ac3a07644b20ccab84391caea4700e52"},
  {"q019-fr-grid-k4-rho1em05-s0",
   "62606f69d9d6171d3d4c8e99ff94ccba8e3424feef0f209769853f3ccb4f85e5"},
  {"q019-fr-grid-k4-rho1em05-s1",
   "614d8bcb11048427926ab44bfea5d087aed42773442e26ea45d26d63f3685a83"},
  {"q019-fr-grid-k4-rho1em05-s2",
   "1146bba3c26825569fa6426681e956dc244413fb304dfe203be54ac5a8b4f18c"},
  {"q019-fr-grid-k8-rho1em06-s0",
   "d3762e89224da830cdf7872ee9b93fcde705afeb59df5c3b5ebf5ec58a79b2f7"},
  {"q019-fr-grid-k8-rho1em06-s1",
   "670666a31c3894f4a2e248192123eefa35b1cf3a8b68557a8b21c4af49e4fdc5"},
  {"q019-fr-grid-k8-rho1em06-s2",
   "0611feeed0225976cec8dafddbb4e9a6c71ffe53876e6292fe1c5ef50e84a2e2"},
  {"q019-fr-grid-k8-rho3em06-s0",
   "fdc59740a5a391c5e97d54e500fb3ba67c963447ce61d2814c8516c03fd37b30"},
  {"q019-fr-grid-k8-rho3em06-s1",
   "bd9f039566b13d4330e1ce50bb55a969cb2ac55d36fd8e336365e3fb4e82e15b"},
  {"q019-fr-grid-k8-rho3em06-s2",
   "dcbc569b3e72b4f5088136cae33c205692967152087062d3af6b4f72563c6bff"},
  {"q019-fr-grid-k8-rho1em05-s0",
   "e25135a6709149f0c8384749cdc11436cb2d3e49df1889bf81a0ae84cd045e71"},
  {"q019-fr-grid-k8-rho1em05-s1",
   "4bc9e145720f1a29ddd534a0f08481ca1459102f99c189f918f7393eeb2747e6"},
  {"q019-fr-grid-k8-rho1em05-s2",
   "1a0b73d814c6f07a9271c323809260f5ddf01fb849b9e0eabd97eceb148cbf22"},
  {"q019-fr-grid-k16-rho1em06-s0",
   "7be9247672d74f3804a01d4792cad4299a2e93d94050d1096e883a0663ac8590"},
  {"q019-fr-grid-k16-rho1em06-s1",
   "1f0d8560bc048399939eeeaded913635326bdb527eb400adf5fd8e1421a8f744"},
  {"q019-fr-grid-k16-rho1em06-s2",
   "8ac21753ed96b45bd80a1039f7f2c5b0e2adbecaf78059512dd908fc621029f0"},
  {"q019-fr-grid-k16-rho3em06-s0",
   "5f5ed92ecad913493d08764fd057a17ff2283250c4962097f0fe555a8172580c"},
  {"q019-fr-grid-k16-rho3em06-s1",
   "4c8513b715e0e8c18ac55dfa3f2444cb649c17432f5b9599a5e5634e69cef1d6"},
  {"q019-fr-grid-k16-rho3em06-s2",
   "fdc7df931f80029cbcf261ed11cb2ea936b2d5af2e6ac39a61f59754cd7cdb58"},
  {"q019-fr-grid-k16-rho1em05-s0",
   "44f9c2b15e9eed74d968384e9777f430e91216f41b3bd8a413106a1e9af9a17e"},
  {"q019-fr-grid-k16-rho1em05-s1",
   "b58262d06ec9dacd0ae074a09ce41af4071e3ad364765fa3b833bcb2e2cddf64"},
  {"q019-fr-grid-k16-rho1em05-s2",
   "2a778c3aafdf77dbea617a4b448b0a2700f12d52c84cedf7f98f50507a9333fe"},
  {"q019-fr-rigidity-isolation-k4-s0",
   "a752617ab4d07bc4b31d5de53639e2107d6d7affad5dfd4001e2a9c84d7030da"},
  {"q019-fr-rigidity-isolation-k8-s0",
   "0a911e8557691df7777f09bd49b4a37c5c55ef78bc9ff12ac3abce4971c0466b"},
  {"q019-fr-rigidity-isolation-k16-s0",
   "656d36ed9878b4e77056646bffa5a2fd8707ad552dc95ed839bf79c56c638118"},
  {"q019-fr-fiducial-ppc24-s0",
   "ef639c7a3e761f1a4945e4e14acc1fa1638124649d6b54f42049a3c3e1c1ab9b"},
  {"q019-fr-fiducial-ppc96-s0",
   "813018c931f5358ffb46c9c262fb1831fa1f35d8e47774a8b78e1cc33ed3a280"},
  {"q019-fr-fiducial-resolution-coarse-s0",
   "6e607bc1932cf6bad712f5881492429f993ab2feedfc87a29dd299255b450034"},
  {"q019-fr-fiducial-resolution-fine-s0",
   "5b94c5a55feea5a5606aafea0d913aa0c1ba27372cce9f882b6a2e16063f340d"},
  {"q019-fr-fiducial-particle-step-small-s0",
   "5204e252307f33794365360278b99de77f2cedb9a194405d7223f23572f23c03"},
  {"q019-fr-fiducial-noise-seeded-s0",
   "0c1dc586e5e189153576c6a696b71e10cb84ef4dbd3945848dc86495140eb341"},
  {"q019-fr-fiducial-riemann-hlld-s0",
   "3e38179885f367a406f479086f4864af8823260fa05c735b87b63c063b4e37c3"},
  {"q019-fr-fiducial-reconstruct-wenoz-s0",
   "ea9fed9133a533f8ce0dbb61d3f52e0518d344462a23f3abec888209b0267c14"},
  {"q019-fr-runtime-initializer-ppc24-s0",
   "71f73a64e2139159816d93f9021350df555618b966e69c0af7bc2da2d51c3fa6"},
  {"q019-fr-runtime-initializer-ppc96-s0",
   "e87a28bfb1fd14a2669a1d480978a6a6e8f7a7c1a8ceb61698b6e95f41cd6b6e"},
  {"q019-fr-3d-onset-small-s0",
   "57ba341bd55e14531afb0f2ca7abd9b11b910eee61a48e947bcbd349f2e277eb"},
  {"q019-fr-3d-onset-large-s0",
   "d2ff1c3bdedfc4927df49bb816c5acb8108f1d073931e4a0263122778cbfc105"},
  {"q019-fr-3d-onset-small-s1",
   "77c30c6173b61dfba27eb973c8185a353da2681c087a4f691f19ba4d345b5e11"},
  {"q019-fr-3d-onset-large-s1",
   "33fab9e7370205adde79ba44b56241e623d81dee127c72c140e2f499ecc459a6"},
  {"q019-fr-3d-onset-small-s2",
   "55bf3d61fca7a935274e355a9c71b91557d84e487095bec085ddca3e2df48840"},
  {"q019-fr-3d-onset-large-s2",
   "e63a208d2488f294ef7bd5ce0be86e675764ac26d6b027055ba085ea4f0df77d"},
  {"q019-fr-3d-onset-small-ppc48-s0",
   "b05065fd57f3bf77084404d2c6e8bbfccfccbaf612831e086f19e15ce3fda760"},
  {"q019-fr-3d-onset-small-resolution-coarse-s0",
   "b52c18e7550e26fdd5bf8b63ad6d7aa7563bd0043c3d782f444a785f05867f42"},
  {"q019-fr-3d-onset-small-resolution-fine-s0",
   "a8a30a52e707cb340dbdf66ede0ea96b92330ff0eb98dde4d407b7c58c454a03"},
  {"q019-fr-3d-onset-small-particle-step-small-s0",
   "99b17e031ebd803026351ee1b7eb3ddfe2c90f7fa7734ef285048681efce0e16"},
  {"q019-fr-3d-onset-large-long-mode-sensitivity-s0",
   "dbb37db10440bb8ed4e0c59b19dc8edae464e6f8d7a26c5a117a9144cc7304a3"},
  {"q019-fr-predecessor-k4-rho1em06-s0",
   "d45231c0c415a1abae767014956ee2076fdda5c98258818cd116e1efab7f6d95"},
  {"q019-fr-predecessor-k4-rho3em06-s0",
   "9a84c2412d8eefced3e0edb500e025f78cb05182db0fed926f906b7b9d006e13"},
  {"q019-fr-predecessor-k4-rho1em05-s0",
   "c4460fa6a8f1994996dda31eb426b39bb98815a3bdb2130d03c11bb90e15ae13"},
  {"q019-fr-predecessor-k8-rho1em06-s0",
   "3a0828f624a6ba09eaf6c991ab9af6dc9f87eab2951c112347de2ff55a609491"},
  {"q019-fr-predecessor-k8-rho3em06-s0",
   "2ba513470763c608a716c35e03553013aac8b8733aa1b92863195fd23a9a2cd8"},
  {"q019-fr-predecessor-k8-rho1em05-s0",
   "5fdf7839c8501bef907a942b697998743ab02b69ef3570da82cc243c57ab3471"},
  {"q019-fr-predecessor-k16-rho1em06-s0",
   "9d7e529dd3bc00575b3f49f0483fd6c7bc17c29001358a59ba7d16cc559ef5ef"},
  {"q019-fr-predecessor-k16-rho3em06-s0",
   "84a3da12481fd0043f6f38008fa243527dc639ba146c81a43731f6aeb484e400"},
  {"q019-fr-predecessor-k16-rho1em05-s0",
   "a7d8394a57fba2ce059bad445a6f1487cedfe3c64e20ca526c50f3a06e3e7e55"},
  {"q019-fr-predecessor-fiducial-ppc24-s0",
   "85eb47817a7a6c8eb178ddb8099a31ee942b2b2c097a4d15441f3b5774a2c8de"},
  {"q019-fr-predecessor-fiducial-ppc48-s0",
   "cf87e50926ec556e038467c2014dd033b5092f8ff02dffdfb1f1c122a8441b82"},
  {"q019-fr-predecessor-fiducial-resolution-coarse-s0",
   "48e395a9d29d987258774a64a200f72211719c738925ea02796ce7ddfb56a259"},
  {"q019-fr-predecessor-fiducial-resolution-fiducial-s0",
   "d50d4585c450cde4bf5444ac699eaede02737cf86a34a2370ff409d6a8e12f27"},
  {"q019-fr-predecessor-fiducial-particle-step-small-s0",
   "84cb50cc0d875e8a77d3c86e3ee8239be35e7e23cc39388b7d20125dac072d17"},
  {"q019-fr-predecessor-fiducial-noise-seeded-s0",
   "626eee951bdebca0e081085ca0ff8d2410f37d9442723ad6b86b61145aab5a6e"},
};

const char *Q019CanonicalMatrixFingerprint(const std::string &case_id) {
  for (const auto &identity : q019_canonical_matrix_identities) {
    if (case_id == identity.case_id) return identity.fingerprint;
  }
  return nullptr;
}

bool Q019RoleMatchesBranchAndDimension(const std::string &role, const Branch branch,
                                       const int dimension) {
  if (branch == Branch::high_rigidity_current_retention_candidate) {
    return dimension == 2 &&
        (role == "high_rigidity_current_retention_seed_ensemble" ||
         role == "high_rigidity_staged_convergence_control");
  }
  if (branch == Branch::high_rigidity_q023_carrier_candidate) {
    if (dimension == 3) {
      return role == "q023_carrier_3d_window_pilot" ||
          role == "q023_carrier_3d_box_control";
    }
    return dimension == 2 &&
        (role == "q023_carrier_linear_onset_pilot" ||
         role == "q023_carrier_nonlinear_window_pilot" ||
         role == "q023_carrier_numerical_control" ||
         role == "q023_carrier_rigidity_control");
  }
  if (branch == Branch::finite_rigidity_early_time_predecessor) {
    return dimension == 2 &&
        (role == "finite_rigidity_early_time_physics_reference_grid" ||
         role == "finite_rigidity_early_time_convergence_control");
  }
  if (dimension == 3) {
    return role == "finite_rigidity_3d_nonlinear_onset_box_small" ||
        role == "finite_rigidity_3d_nonlinear_onset_box_large" ||
        role == "finite_rigidity_3d_onset_matched_convergence_control" ||
        role == "finite_rigidity_3d_spectral_sensitivity_control";
  }
  return dimension == 2 &&
      (role == "finite_rigidity_density_drift_coupled_response_ensemble" ||
       role == "finite_rigidity_isolation_fixed_cr_distribution_control" ||
       role == "finite_rigidity_numerical_or_noise_control" ||
       role == "finite_rigidity_source_local_initializer_regression");
}

int q019_runtime_dimension = 0;
ParameterInput *q019_runtime_pin = nullptr;
bool q019_runtime_box_edge_monitor_enabled = false;
bool q019_runtime_box_edge_stop_armed = false;
int q019_runtime_box_edge_stop_ppm = -1;
int q019_runtime_box_edge_last_cycle = 0;
int q019_runtime_box_edge_completed_slots = 0;
int q019_runtime_box_edge_valid_samples = 0;
int q019_runtime_box_edge_skipped_slots = 0;
int q019_runtime_box_edge_last_status =
    static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled);
int q019_runtime_box_edge_status_mask = 0;
Real q019_runtime_box_edge_monitor_dt = 0.0;
Real q019_runtime_box_edge_next_nominal_time = 0.0;
Real q019_runtime_box_edge_last_nominal_time = -1.0;
Real q019_runtime_box_edge_last_prior_time = 0.0;
Real q019_runtime_box_edge_last_time = 0.0;
Real q019_runtime_box_edge_last_power_fraction = -1.0;
Real q019_runtime_box_edge_max_power_fraction = -1.0;
Real q019_runtime_box_edge_last_fluctuation_mean = -1.0;
bool q019_runtime_controller_enabled = false;
bool q019_runtime_controller_box_edge_monitor_enabled = false;
bool q019_runtime_resolution_monitor_enabled = false;
bool q019_runtime_diagnostic_failure_stop_armed = false;
bool q019_runtime_resolution_stop_armed = false;
bool q019_runtime_controller_box_edge_stop_armed = false;
bool q019_runtime_controller_triggered = false;
bool q019_runtime_controller_trigger_failure = false;
int q019_runtime_controller_trigger_reason = 0;
int q019_runtime_controller_trigger_cycle = 0;
int q019_runtime_resolution_samples = 0;
int q019_runtime_resolution_last_cycle = 0;
int q019_runtime_controller_box_edge_stop_ppm = -1;
int q019_runtime_controller_pilot_cycle_limit = -1;
Real q019_runtime_background_b = 0.0;
Real q019_runtime_controller_monitor_dt = 0.0;
Real q019_runtime_resolution_stop_b_over_b0 = -1.0;
Real q019_runtime_resolution_last_time = 0.0;
Real q019_runtime_resolution_last_b_over_b0 = -1.0;
Real q019_runtime_resolution_max_b_over_b0 = -1.0;
Real q019_runtime_controller_trigger_time = 0.0;
Real q019_runtime_controller_trigger_metric = -1.0;

void Q019StoreRuntimeControllerState() {
  if (!q019_runtime_controller_enabled || q019_runtime_pin == nullptr) return;
  const std::string block = kQ019RuntimeControllerBlock;
  q019_runtime_pin->SetInteger(block, "runtime_resolution_samples",
                               q019_runtime_resolution_samples);
  q019_runtime_pin->SetInteger(block, "runtime_resolution_last_cycle",
                               q019_runtime_resolution_last_cycle);
  q019_runtime_pin->SetReal(block, "runtime_resolution_last_time",
                            q019_runtime_resolution_last_time);
  q019_runtime_pin->SetReal(block, "runtime_resolution_last_B_over_B0",
                            q019_runtime_resolution_last_b_over_b0);
  q019_runtime_pin->SetReal(block, "runtime_resolution_max_B_over_B0",
                            q019_runtime_resolution_max_b_over_b0);
  q019_runtime_pin->SetBoolean(block, "runtime_controller_triggered",
                               q019_runtime_controller_triggered);
  q019_runtime_pin->SetBoolean(block, "runtime_controller_trigger_failure",
                               q019_runtime_controller_trigger_failure);
  q019_runtime_pin->SetInteger(block, "runtime_controller_trigger_reason",
                               q019_runtime_controller_trigger_reason);
  q019_runtime_pin->SetInteger(block, "runtime_controller_trigger_cycle",
                               q019_runtime_controller_trigger_cycle);
  q019_runtime_pin->SetReal(block, "runtime_controller_trigger_time",
                            q019_runtime_controller_trigger_time);
  q019_runtime_pin->SetReal(block, "runtime_controller_trigger_metric",
                            q019_runtime_controller_trigger_metric);
}

void Q019RequestRuntimeControllerStop(Mesh *pm, const int reason,
                                      const bool failure, const Real metric) {
  if (pm == nullptr || pm->pgen == nullptr) {
    Q019NonlinearFatal("Q019 runtime controller cannot address the driver");
  }
  const int completed_cycle = pm->ncycle + 1;
  const Real completed_time = pm->time + pm->dt;
  if (q019_runtime_controller_triggered) {
    if (q019_runtime_controller_trigger_reason != reason ||
        q019_runtime_controller_trigger_failure != failure) {
      Q019NonlinearFatal("Q019 runtime controllers requested conflicting stops");
    }
    return;
  }
  q019_runtime_controller_triggered = true;
  q019_runtime_controller_trigger_failure = failure;
  q019_runtime_controller_trigger_reason = reason;
  q019_runtime_controller_trigger_cycle = completed_cycle;
  q019_runtime_controller_trigger_time = completed_time;
  q019_runtime_controller_trigger_metric = metric;
  Q019StoreRuntimeControllerState();
  pm->pgen->RequestUserStop(reason, failure);
}

void Q019StoreBoxEdgeMonitorState() {
  if (q019_runtime_pin == nullptr) {
    Q019NonlinearFatal("Q019 box-edge monitor has no restart-state ParameterInput");
  }
  const std::string block = "q019_physics_first_nonlinear_bell_successor_v2";
  q019_runtime_pin->SetInteger(block, "runtime_box_edge_monitor_schema", 2);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_last_cycle",
      q019_runtime_box_edge_last_cycle);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_completed_slots",
      q019_runtime_box_edge_completed_slots);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_valid_samples",
      q019_runtime_box_edge_valid_samples);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_skipped_slots",
      q019_runtime_box_edge_skipped_slots);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_last_status",
      q019_runtime_box_edge_last_status);
  q019_runtime_pin->SetInteger(
      block, "runtime_box_edge_monitor_status_mask",
      q019_runtime_box_edge_status_mask);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_next_nominal_time",
      q019_runtime_box_edge_next_nominal_time);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_last_nominal_time",
      q019_runtime_box_edge_last_nominal_time);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_last_prior_time",
      q019_runtime_box_edge_last_prior_time);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_last_time",
      q019_runtime_box_edge_last_time);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_last_power_fraction",
      q019_runtime_box_edge_last_power_fraction);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_max_power_fraction",
      q019_runtime_box_edge_max_power_fraction);
  q019_runtime_pin->SetReal(
      block, "runtime_box_edge_monitor_last_fluctuation_mean",
      q019_runtime_box_edge_last_fluctuation_mean);
}

const char *Q019BoxEdgeStatusName(const int status) {
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled)) {
    return "not_sampled";
  }
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::valid)) return "valid";
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::cadence_skipped)) {
    return "cadence_skipped";
  }
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::no_global_cells)) {
    return "no_global_cells";
  }
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation)) {
    return "zero_fluctuation";
  }
  if (status == static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable)) {
    return "numerical_unavailable";
  }
  return "invalid_status";
}

void Q019RecordBoxEdgeDiagnostic(
    const int status, const int crossed_slots, const int completed_cycle,
    const Real prior_completed_time, const Real completed_time,
    const Real power_fraction, const Real fluctuation_mean) {
  if (!BoxEdgeStatusCodeIsValid(status) || crossed_slots <= 0) return;
  q019_runtime_box_edge_completed_slots += crossed_slots;
  q019_runtime_box_edge_skipped_slots += std::max(0, crossed_slots - 1);
  q019_runtime_box_edge_last_cycle = completed_cycle;
  q019_runtime_box_edge_last_prior_time = prior_completed_time;
  q019_runtime_box_edge_last_time = completed_time;
  q019_runtime_box_edge_last_status = status;
  q019_runtime_box_edge_status_mask |= BoxEdgeStatusBit(status);
  q019_runtime_box_edge_last_nominal_time = static_cast<Real>(
      BoxEdgeLastNominalTime(q019_runtime_box_edge_completed_slots,
                             q019_runtime_box_edge_monitor_dt));
  q019_runtime_box_edge_next_nominal_time = static_cast<Real>(
      BoxEdgeNextNominalTime(q019_runtime_box_edge_completed_slots,
                             q019_runtime_box_edge_monitor_dt));
  if (BoxEdgeStatusRecordsValidMetric(status)) {
    q019_runtime_box_edge_last_power_fraction = power_fraction;
    q019_runtime_box_edge_max_power_fraction =
        std::max(q019_runtime_box_edge_max_power_fraction, power_fraction);
    q019_runtime_box_edge_last_fluctuation_mean = fluctuation_mean;
    ++q019_runtime_box_edge_valid_samples;
  } else {
    q019_runtime_box_edge_last_power_fraction = -1.0;
    q019_runtime_box_edge_last_fluctuation_mean =
        status == static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation) ?
        0.0 : -1.0;
  }
  Q019StoreBoxEdgeMonitorState();
  if (status != static_cast<int>(BoxEdgeDiagnosticStatus::valid) &&
      global_variable::my_rank == 0) {
    std::cout << "Q019_BOX_EDGE_DIAGNOSTIC_UNAVAILABLE status="
              << Q019BoxEdgeStatusName(status)
              << " crossed_slots=" << crossed_slots
              << " completed_cycle=" << completed_cycle
              << " prior_completed_time=" << prior_completed_time
              << " completed_time=" << completed_time << std::endl;
  }
}

void Q019RuntimeBoxEdgeMonitor(Mesh *pm) {
  if (!q019_runtime_box_edge_monitor_enabled || pm == nullptr ||
      pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr) {
    return;
  }
  const int completed_cycle = pm->ncycle + 1;
  const Real prior_completed_time = pm->time;
  const Real completed_time = pm->time + pm->dt;
  const int crossed_slots = BoxEdgeCrossedSlotCount(
      q019_runtime_box_edge_completed_slots, q019_runtime_box_edge_monitor_dt,
      completed_time);
  if (crossed_slots == 0) return;
  if (crossed_slots < 0) {
    q019_runtime_box_edge_last_status =
        static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable);
    q019_runtime_box_edge_status_mask |=
        BoxEdgeStatusBit(q019_runtime_box_edge_last_status);
    Q019StoreBoxEdgeMonitorState();
    if (global_variable::my_rank == 0) {
      std::cout << "Q019_BOX_EDGE_DIAGNOSTIC_UNAVAILABLE status="
                << Q019BoxEdgeStatusName(q019_runtime_box_edge_last_status)
                << " reason=invalid_crossed_slot_count" << std::endl;
    }
    return;
  }
  if (crossed_slots > 1) {
    Q019RecordBoxEdgeDiagnostic(
        static_cast<int>(BoxEdgeDiagnosticStatus::cadence_skipped),
        crossed_slots, completed_cycle, prior_completed_time, completed_time,
        -1.0, -1.0);
    return;
  }
  if (!BoxEdgeFirstCrossingChronologyIsValid(
          prior_completed_time, q019_runtime_box_edge_next_nominal_time,
          completed_time)) {
    Q019RecordBoxEdgeDiagnostic(
        static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable),
        crossed_slots, completed_cycle, prior_completed_time, completed_time,
        -1.0, -1.0);
    return;
  }

  auto *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &bcc = pmbp->pmhd->bcc0;
  auto &size = pmbp->pmb->mb_size;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const int nmkji = nmb*nkji;
  array_sum::GlobalSum transverse_sums;
  Kokkos::parallel_reduce(
      "q019_box_edge_transverse_sums",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
    array_sum::GlobalSum values;
    values.the_array[0] = bcc(m, IBY, k, j, i);
    values.the_array[1] = bcc(m, IBZ, k, j, i);
    values.the_array[2] = 1.0;
    sum += values;
  }, Kokkos::Sum<array_sum::GlobalSum>(transverse_sums));
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, transverse_sums.the_array, 3, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
#endif
  const Real global_cells = transverse_sums.the_array[2];
  if (global_cells <= 0) {
    Q019RecordBoxEdgeDiagnostic(
        static_cast<int>(BoxEdgeDiagnosticStatus::no_global_cells),
        crossed_slots, completed_cycle, prior_completed_time, completed_time,
        -1.0, -1.0);
    return;
  }
  const Real mean_b2 = transverse_sums.the_array[0]/global_cells;
  const Real mean_b3 = transverse_sums.the_array[1]/global_cells;
  const int dimension = q019_runtime_dimension;
  const int unique_mode_count = BoxEdgeUniqueModeCount(dimension);
  const int reduction_group_count = BoxEdgeReductionGroupCount(dimension);
  const int last_group_first_mode =
      (reduction_group_count - 1)*kBoxEdgeModesPerReduction;
  const int fluctuation_reduction_index =
      4*(unique_mode_count - last_group_first_mode);
  Real mode_amplitudes[kBoxEdgeMaximumUniqueModeCount][4] = {};
  Real fluctuation_sum = 0.0;
  const Real x1min = pm->mesh_size.x1min;
  const Real x2min = pm->mesh_size.x2min;
  const Real x3min = pm->mesh_size.x3min;
  const Real x1extent = pm->mesh_size.x1max - x1min;
  const Real x2extent = pm->mesh_size.x2max - x2min;
  const Real x3extent = pm->mesh_size.x3max - x3min;
  const Real two_pi = 2.0*std::acos(-1.0);
  for (int n = 0; n < unique_mode_count; ++n) {
    if (!BoxEdgePhysicalLowKModeIsSelected(
            BoxEdgeUniqueModeAt(n, dimension), dimension,
            x1extent, x2extent, x3extent)) {
      Q019RecordBoxEdgeDiagnostic(
          static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable),
          crossed_slots, completed_cycle, prior_completed_time, completed_time,
          -1.0, -1.0);
      return;
    }
  }
  size.template sync<DevExeSpace>();
  auto size_view = size;
  for (int group = 0; group < reduction_group_count; ++group) {
    const int first_mode = group*kBoxEdgeModesPerReduction;
    array_sum::GlobalSum group_sums;
    Kokkos::parallel_reduce(
        "q019_box_edge_grouped_mode_amplitudes",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
      const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
      const Real x1 = CellCenterX(i - is, nx1, size_view.d_view(m).x1min,
                                  size_view.d_view(m).x1max);
      const Real x2 = CellCenterX(j - js, nx2, size_view.d_view(m).x2min,
                                  size_view.d_view(m).x2max);
      const Real x3 = CellCenterX(k - ks, nx3, size_view.d_view(m).x3min,
                                  size_view.d_view(m).x3max);
      const Real theta1 = two_pi*(x1 - x1min)/x1extent;
      const Real theta2 = two_pi*(x2 - x2min)/x2extent;
      const Real theta3 = dimension == 3 ?
          two_pi*(x3 - x3min)/x3extent : 0.0;
      const ComplexValue x1_phase = {cos(theta1), -sin(theta1)};
      const ComplexValue x2_phase = {cos(theta2), -sin(theta2)};
      const ComplexValue x3_phase = {cos(theta3), -sin(theta3)};
      const Real db2 = bcc(m, IBY, k, j, i) - mean_b2;
      const Real db3 = bcc(m, IBZ, k, j, i) - mean_b3;
      array_sum::GlobalSum values;
      for (int slot = 0; slot < kBoxEdgeModesPerReduction; ++slot) {
        const int n = first_mode + slot;
        if (n >= unique_mode_count) continue;
        const IntegerVector3 mode = BoxEdgeUniqueModeAt(n, dimension);
        const ComplexValue phase =
            BoxEdgeModePhase(x1_phase, x2_phase, x3_phase, mode);
        const int offset = 4*slot;
        values.the_array[offset    ] = db2*phase.real;
        values.the_array[offset + 1] = db2*phase.imag;
        values.the_array[offset + 2] = db3*phase.real;
        values.the_array[offset + 3] = db3*phase.imag;
      }
      if (group == reduction_group_count - 1) {
        values.the_array[fluctuation_reduction_index] = db2*db2 + db3*db3;
      }
      sum += values;
    }, Kokkos::Sum<array_sum::GlobalSum>(group_sums));
    Kokkos::fence();
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, group_sums.the_array, NREDUCTION_VARIABLES,
                  MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    for (int slot = 0; slot < kBoxEdgeModesPerReduction; ++slot) {
      const int n = first_mode + slot;
      if (n >= unique_mode_count) continue;
      for (int component = 0; component < 4; ++component) {
        mode_amplitudes[n][component] =
            group_sums.the_array[4*slot + component];
      }
    }
    if (group == reduction_group_count - 1) {
      fluctuation_sum = group_sums.the_array[fluctuation_reduction_index];
    }
  }
  Real unique_mode_power_sum = 0.0;
  for (int n = 0; n < unique_mode_count; ++n) {
    for (int component = 0; component < 4; ++component) {
      unique_mode_power_sum +=
          mode_amplitudes[n][component]*mode_amplitudes[n][component];
    }
  }
  if (fluctuation_sum == 0.0) {
    Q019RecordBoxEdgeDiagnostic(
        static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation),
        crossed_slots, completed_cycle, prior_completed_time, completed_time,
        -1.0, 0.0);
    return;
  }
  Real fraction = BoxEdgePowerFraction(
      global_cells, fluctuation_sum, unique_mode_power_sum);
  const Real tolerance = static_cast<Real>(4096.0)*
      std::numeric_limits<Real>::epsilon();
  if (!std::isfinite(fraction) || fraction < -tolerance ||
      fraction > 1.0 + tolerance || !std::isfinite(fluctuation_sum) ||
      !std::isfinite(unique_mode_power_sum)) {
    Q019RecordBoxEdgeDiagnostic(
        static_cast<int>(BoxEdgeDiagnosticStatus::numerical_unavailable),
        crossed_slots, completed_cycle, prior_completed_time, completed_time,
        -1.0, -1.0);
    return;
  }
  fraction = std::max(static_cast<Real>(0.0),
                      std::min(static_cast<Real>(1.0), fraction));
  Q019RecordBoxEdgeDiagnostic(
      static_cast<int>(BoxEdgeDiagnosticStatus::valid), crossed_slots,
      completed_cycle, prior_completed_time, completed_time, fraction,
      fluctuation_sum/global_cells);
}

void Q019RuntimeResolutionMonitor(Mesh *pm) {
  if (!q019_runtime_controller_enabled ||
      !q019_runtime_resolution_monitor_enabled || pm == nullptr ||
      pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr) {
    return;
  }
  const int completed_cycle = pm->ncycle + 1;
  if (q019_runtime_box_edge_last_cycle != completed_cycle) return;

  auto *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &bcc = pmbp->pmhd->bcc0;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  Real maximum_b2 = 0.0;
  Kokkos::parallel_reduce(
      "q019_runtime_resolution_maximum_b2",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &value) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
    const Real b1 = bcc(m, IBX, k, j, i);
    const Real b2 = bcc(m, IBY, k, j, i);
    const Real b3 = bcc(m, IBZ, k, j, i);
    value = fmax(value, b1*b1 + b2*b2 + b3*b3);
  }, Kokkos::Max<Real>(maximum_b2));
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &maximum_b2, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
#endif
  const Real sampled_b_over_b0 =
      std::sqrt(maximum_b2)/q019_runtime_background_b;
  if (!std::isfinite(sampled_b_over_b0) || sampled_b_over_b0 <= 0.0) {
    if (q019_runtime_diagnostic_failure_stop_armed) {
      Q019RequestRuntimeControllerStop(
          pm, kQ019DiagnosticFailureReason, true, -1.0);
    }
    return;
  }
  ++q019_runtime_resolution_samples;
  q019_runtime_resolution_last_cycle = completed_cycle;
  q019_runtime_resolution_last_time = pm->time + pm->dt;
  q019_runtime_resolution_last_b_over_b0 = sampled_b_over_b0;
  q019_runtime_resolution_max_b_over_b0 = std::max(
      q019_runtime_resolution_max_b_over_b0, sampled_b_over_b0);
  Q019StoreRuntimeControllerState();
  if (q019_runtime_resolution_stop_armed &&
      sampled_b_over_b0 >= q019_runtime_resolution_stop_b_over_b0) {
    Q019RequestRuntimeControllerStop(
        pm, kQ019ResolutionStopReason, false, sampled_b_over_b0);
  }
}

void Q019RuntimeDiagnostics(Mesh *pm) {
  const int completed_cycle = pm == nullptr ? 0 : pm->ncycle + 1;
  const int slots_before = q019_runtime_box_edge_completed_slots;
  const int expected_crossings =
      q019_runtime_controller_box_edge_monitor_enabled && pm != nullptr ?
      BoxEdgeCrossedSlotCount(slots_before, q019_runtime_box_edge_monitor_dt,
                              pm->time + pm->dt) : 0;
  Q019RuntimeBoxEdgeMonitor(pm);
  if (q019_runtime_diagnostic_failure_stop_armed &&
      q019_runtime_controller_box_edge_monitor_enabled &&
      expected_crossings != 0 &&
      q019_runtime_box_edge_last_cycle != completed_cycle) {
    Q019RequestRuntimeControllerStop(
        pm, kQ019DiagnosticFailureReason, true,
        static_cast<Real>(q019_runtime_box_edge_last_status));
    return;
  }
  Q019RuntimeResolutionMonitor(pm);
  if (!q019_runtime_controller_enabled ||
      q019_runtime_controller_triggered || pm == nullptr) {
    return;
  }
  if (q019_runtime_controller_box_edge_stop_armed &&
      q019_runtime_box_edge_last_cycle == completed_cycle) {
    if (q019_runtime_box_edge_last_status ==
        static_cast<int>(BoxEdgeDiagnosticStatus::valid)) {
      const Real threshold = static_cast<Real>(
          q019_runtime_controller_box_edge_stop_ppm)*1.0e-6;
      if (q019_runtime_box_edge_last_power_fraction >= threshold) {
        Q019RequestRuntimeControllerStop(
            pm, kQ019BoxEdgeStopReason, false,
            q019_runtime_box_edge_last_power_fraction);
      }
    } else {
      Q019RequestRuntimeControllerStop(
          pm, kQ019DiagnosticFailureReason, true,
          static_cast<Real>(q019_runtime_box_edge_last_status));
    }
  }
  if (!q019_runtime_controller_triggered &&
      q019_runtime_controller_pilot_cycle_limit > 0 &&
      completed_cycle >= q019_runtime_controller_pilot_cycle_limit) {
    Q019RequestRuntimeControllerStop(
        pm, kQ019PilotCycleStopReason, false,
        static_cast<Real>(completed_cycle));
  }
}

void Q019BoxEdgeHistory(HistoryData *pdata, Mesh *) {
  pdata->nhist = 12;
  if (q019_runtime_controller_enabled) {
    const char *labels[12] = {
      "boxedge", "boxmax", "boxstat", "boxmask",
      "resB", "resBmax", "rescyc", "ressamp",
      "stopwhy", "stopcyc", "stoptime", "stopval"
    };
    for (int n = 0; n < pdata->nhist; ++n) {
      pdata->label[n] = labels[n];
      pdata->hdata[n] = 0.0;
    }
    if (global_variable::my_rank != 0) return;
    pdata->hdata[0] = q019_runtime_box_edge_last_power_fraction;
    pdata->hdata[1] = q019_runtime_box_edge_max_power_fraction;
    pdata->hdata[2] = static_cast<Real>(q019_runtime_box_edge_last_status);
    pdata->hdata[3] = static_cast<Real>(q019_runtime_box_edge_status_mask);
    pdata->hdata[4] = q019_runtime_resolution_last_b_over_b0;
    pdata->hdata[5] = q019_runtime_resolution_max_b_over_b0;
    pdata->hdata[6] = static_cast<Real>(q019_runtime_resolution_last_cycle);
    pdata->hdata[7] = static_cast<Real>(q019_runtime_resolution_samples);
    pdata->hdata[8] = static_cast<Real>(q019_runtime_controller_trigger_reason);
    pdata->hdata[9] = static_cast<Real>(q019_runtime_controller_trigger_cycle);
    pdata->hdata[10] = q019_runtime_controller_trigger_time;
    pdata->hdata[11] = q019_runtime_controller_trigger_metric;
    return;
  }
  pdata->label[0] = "boxedge";
  pdata->label[1] = "boxmax";
  pdata->label[2] = "boxcyc";
  pdata->label[3] = "boxprior";
  pdata->label[4] = "boxtime";
  pdata->label[5] = "boxfluc";
  pdata->label[6] = "boxslots";
  pdata->label[7] = "boxvalid";
  pdata->label[8] = "boxskip";
  pdata->label[9] = "boxstat";
  pdata->label[10] = "boxmask";
  pdata->label[11] = "boxnext";
  for (int n = 0; n < pdata->nhist; ++n) pdata->hdata[n] = 0.0;
  if (global_variable::my_rank != 0) return;
  pdata->hdata[0] = q019_runtime_box_edge_last_power_fraction;
  pdata->hdata[1] = q019_runtime_box_edge_max_power_fraction;
  pdata->hdata[2] = static_cast<Real>(q019_runtime_box_edge_last_cycle);
  pdata->hdata[3] = q019_runtime_box_edge_last_prior_time;
  pdata->hdata[4] = q019_runtime_box_edge_last_time;
  pdata->hdata[5] = q019_runtime_box_edge_last_fluctuation_mean;
  pdata->hdata[6] = static_cast<Real>(q019_runtime_box_edge_completed_slots);
  pdata->hdata[7] = static_cast<Real>(q019_runtime_box_edge_valid_samples);
  pdata->hdata[8] = static_cast<Real>(q019_runtime_box_edge_skipped_slots);
  pdata->hdata[9] = static_cast<Real>(q019_runtime_box_edge_last_status);
  pdata->hdata[10] = static_cast<Real>(q019_runtime_box_edge_status_mask);
  pdata->hdata[11] = q019_runtime_box_edge_next_nominal_time;
}
void Q019FinalEvidenceStatus(ParameterInput *, Mesh *) {
  if (global_variable::my_rank != 0) return;
  std::cout << "Q019_BOX_EDGE_MONITOR_COMPLETED_SLOTS="
            << q019_runtime_box_edge_completed_slots << std::endl;
  std::cout << "Q019_BOX_EDGE_MONITOR_VALID_SAMPLES="
            << q019_runtime_box_edge_valid_samples << std::endl;
  std::cout << "Q019_BOX_EDGE_MONITOR_SKIPPED_SLOTS="
            << q019_runtime_box_edge_skipped_slots << std::endl;
  std::cout << "Q019_BOX_EDGE_MONITOR_LAST_STATUS="
            << Q019BoxEdgeStatusName(q019_runtime_box_edge_last_status) << std::endl;
  std::cout << "Q019_BOX_EDGE_MONITOR_STATUS_MASK="
            << q019_runtime_box_edge_status_mask << std::endl;
  std::cout << "Q019_BOX_EDGE_LAST_POWER_FRACTION="
            << q019_runtime_box_edge_last_power_fraction << std::endl;
  std::cout << "Q019_BOX_EDGE_MAX_POWER_FRACTION="
            << q019_runtime_box_edge_max_power_fraction << std::endl;
  std::cout << "Q019_RUNTIME_RESOLUTION_SAMPLES="
            << q019_runtime_resolution_samples << std::endl;
  std::cout << "Q019_RUNTIME_RESOLUTION_MAX_B_OVER_B0="
            << q019_runtime_resolution_max_b_over_b0 << std::endl;
  std::cout << "Q019_RUNTIME_CONTROLLER_TRIGGERED="
            << (q019_runtime_controller_triggered ? "true" : "false") << std::endl;
  std::cout << "Q019_RUNTIME_CONTROLLER_TRIGGER_REASON="
            << q019_runtime_controller_trigger_reason << std::endl;
  std::cout << "Q019_FINAL_EVIDENCE_STATUS=completed_not_acceptance_eligible"
            << std::endl;
  std::cout << "Q019_SATURATION_EVIDENCE_ELIGIBLE=false" << std::endl;
}

Branch ParseBranch(const std::string &name) {
  if (name.compare("high_rigidity_current_retention_candidate") == 0) {
    return Branch::high_rigidity_current_retention_candidate;
  }
  if (name.compare("high_rigidity_q023_carrier_candidate") == 0) {
    return Branch::high_rigidity_q023_carrier_candidate;
  }
  if (name.compare("finite_rigidity_self_consistent") == 0) {
    return Branch::finite_rigidity_self_consistent;
  }
  if (name.compare("finite_rigidity_early_time_predecessor") == 0) {
    return Branch::finite_rigidity_early_time_predecessor;
  }
  Q019NonlinearFatal("<q019_physics_first_nonlinear_bell_successor_v2>/branch "
                     "must select one explicit nonlinear Bell branch");
}

SeedTopology ParseSeedTopology(const std::string &name) {
  if (name.compare("shared_spectrum_only") == 0) {
    return SeedTopology::shared_spectrum_only;
  }
  if (name.compare("shared_spectrum_plus_separated_long_modes") == 0) {
    return SeedTopology::shared_spectrum_plus_separated_long_modes;
  }
  Q019NonlinearFatal("Q019 seed topology is invalid");
}

FiniteSamplingMode ParseFiniteSamplingMode(const std::string &name) {
  if (name.compare("not_applicable") == 0) {
    return FiniteSamplingMode::not_applicable;
  }
  if (name.compare("cell_centered_nested_haar_octahedral_packets") == 0) {
    return FiniteSamplingMode::cell_centered_nested_haar_octahedral_packets;
  }
  if (name.compare("independent_position_noise_seeded") == 0) {
    return FiniteSamplingMode::independent_position_noise_seeded;
  }
  Q019NonlinearFatal("Q019 finite sampling mode is invalid");
}

}  // namespace

void ProblemGenerator::Q019PhysicsFirstNonlinearBellSuccessorV2(
    ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->ppart == nullptr) {
    Q019NonlinearFatal("corrected nonlinear Bell successor requires MHD and particles");
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    Q019NonlinearFatal("corrected nonlinear Bell successor requires ideal-MHD EOS");
  }
  if (pmy_mesh_->multilevel) {
    Q019NonlinearFatal("corrected nonlinear Bell successor rejects AMR/SMR");
  }
  if (!pmy_mesh_->strictly_periodic) {
    Q019NonlinearFatal("corrected nonlinear Bell successor requires periodic boundaries");
  }
  const int mesh_dimension = pmy_mesh_->two_d ? 2 : (pmy_mesh_->three_d ? 3 : 1);
  if (mesh_dimension < 2) {
    Q019NonlinearFatal("corrected nonlinear Bell successor requires a 2D3V or 3D mesh");
  }

  const std::string block = "q019_physics_first_nonlinear_bell_successor_v2";
  const int dimension = pin->GetInteger(block, "dimension");
  if (dimension != mesh_dimension) {
    Q019NonlinearFatal("<" + block + ">/dimension must match the active mesh");
  }
  const std::string case_id = pin->GetString(block, "case_id");
  const std::string campaign_id = pin->GetString(block, "campaign_id");
  const std::string branch_name = pin->GetString(block, "branch");
  const std::string role = pin->GetString(block, "role");
  const std::string seed_topology_name = pin->GetString(block, "seed_topology");
  const std::string finite_sampling_mode_name =
      pin->GetString(block, "finite_sampling_mode");
  const std::string high_sampling_mode = pin->GetString(block, "high_sampling_mode");
  const int pic_max_cell_cross = pin->GetInteger("particles", "pic_max_cell_cross");
  const std::string box_pair_id = pin->GetString(block, "box_pair_id");
  const bool saturation_candidate = pin->GetBoolean(block, "saturation_candidate");
  const bool spectral_sensitivity_control =
      pin->GetBoolean(block, "spectral_sensitivity_control");
  if (case_id.rfind("q019-", 0) != 0) {
    Q019NonlinearFatal("Q019 case_id must use the q019- namespace");
  }
  Q019RequireString(pin, "job", "basename", case_id);
  const Branch branch = ParseBranch(branch_name);
  const bool high_rigidity =
      branch == Branch::high_rigidity_current_retention_candidate ||
      branch == Branch::high_rigidity_q023_carrier_candidate;
  const bool q023_carrier =
      branch == Branch::high_rigidity_q023_carrier_candidate;
  if (!Q019RoleMatchesBranchAndDimension(role, branch, dimension)) {
    Q019NonlinearFatal("Q019 role does not match its branch and dimension");
  }
  const SeedTopology seed_topology = ParseSeedTopology(seed_topology_name);
  const FiniteSamplingMode finite_sampling_mode = ParseFiniteSamplingMode(
      finite_sampling_mode_name);
  if (!FiniteSamplingModeIsValid(branch, finite_sampling_mode)) {
    Q019NonlinearFatal("finite sampling mode leaked across Q019 branches");
  }
  const bool quiet_isotropic_shell_packet =
      finite_sampling_mode ==
      FiniteSamplingMode::cell_centered_nested_haar_octahedral_packets;
  const bool high_stochastic_sampling =
      high_sampling_mode.compare("seeded_random_cold_beam_noise_control") == 0;
  if ((high_rigidity &&
       high_sampling_mode.compare("cell_centered_cold_beam") != 0 &&
       !high_stochastic_sampling) ||
      (!high_rigidity && high_sampling_mode.compare("not_applicable") != 0)) {
    Q019NonlinearFatal("high-rigidity sampling mode leaked across branches");
  }
  Q019RequireString(
      pin, block, "high_position_noise_reference_case_id",
      high_stochastic_sampling ? "q019-hr-fiducial-ppc8-centered-s0" :
      "not_applicable");
  const bool finite_predecessor =
      branch == Branch::finite_rigidity_early_time_predecessor;
  const std::string expected_campaign_id =
      q023_carrier ?
      "Q019-HR-JOVERC-Q023-CARRIER-FIXED-CURRENT-LIKE-NOHALL-V1"
                   : high_rigidity ?
      "Q019-HR-SIMILARITY-MAPPED-HALL-OMISSION-CANDIDATE-NONSHOCK-V4"
                    : (finite_predecessor ?
                       "Q019-FR-ISOTROPIC-SHELL-EARLY-TIME-PREDECESSOR-V4" :
                       "Q019-FR-ISOTROPIC-SHELL-ONSET-CANDIDATE-NONSHOCK-V4");
  Q019RequireString(pin, block, "campaign_id", expected_campaign_id);
  if (campaign_id != expected_campaign_id) {
    Q019NonlinearFatal("Q019 campaign identity drifted");
  }
  if (saturation_candidate) {
    Q019NonlinearFatal("Q019 current matrix cannot contain saturation candidates");
  }
  const bool expected_spectral_sensitivity =
      role == "finite_rigidity_3d_spectral_sensitivity_control";
  if (spectral_sensitivity_control != expected_spectral_sensitivity ||
      spectral_sensitivity_control !=
          (seed_topology == SeedTopology::shared_spectrum_plus_separated_long_modes)) {
    Q019NonlinearFatal("Q019 spectral-sensitivity identity drifted");
  }
  const bool onset_box_role =
      role == "finite_rigidity_3d_nonlinear_onset_box_small" ||
      role == "finite_rigidity_3d_nonlinear_onset_box_large" ||
      role == "finite_rigidity_3d_onset_matched_convergence_control" ||
      role == "finite_rigidity_3d_spectral_sensitivity_control";
  if ((onset_box_role && box_pair_id.rfind("finite-3d-onset-s", 0) != 0) ||
      (!onset_box_role && box_pair_id != "not_applicable")) {
    Q019NonlinearFatal("Q019 box-pair identity drifted");
  }
  Q019RequireString(pin, block, "deck_role",
                    "source_local_physics_first_candidate_not_authorized");
  Q019RequireString(
      pin, block, "source_lineage",
      q023_carrier ?
      "hardened_q043_then_registered_q023_joverc_then_q019_q023_carrier_v1" :
      "hardened_control_plane_q043_then_corrected_q023_then_q019_v4");
  Q019RequireString(
      pin, block, "matrix_scope",
      q023_carrier ? "q023_carrier_excluded_physical_pilot_design_v1" :
      "physics_first_preproduction_design_successor_v4");
  Q019RequireString(
      pin, block, "domain_time_status",
      q023_carrier ?
      "versioned_redesign_after_finite_rigidity_resource_overrun_execution_prohibited" :
      "engineering_candidate_pending_excluded_resource_and_window_pilots");
  Q019RequireString(
      pin, block, "applicability_scope",
      q023_carrier ?
      "q023_low_gyrofrequency_large_inertia_external_current_surrogate_"
      "periodic_ideal_mhd_no_finite_rigidity_or_physical_cr_density_claim" :
      "similarity_scaled_equal_cr_background_qom_mhd_scale_and_R_explicit_"
      "applicability_review_pending_not_strong_shock_mapping");
  Q019RequireString(
      pin, block, "energy_loading_regime",
      q023_carrier ?
      "large_inertial_carrier_reservoir_explicitly_reported_not_physical_cr_density" :
      "isotropic_shell_plus_drift_energy_accounted_periodic_system");
  Q019RequireString(
      pin, block, "energy_loading_gate_status",
      q023_carrier ?
      "external_current_surrogate_only_requires_measured_current_momentum_energy_"
      "invariance_and_conservation" :
      "unqualified_requires_coupled_response_grid_conservation_and_saturation_review");
  Q019RequireString(
      pin, block, "finite_distribution",
      high_rigidity ? "not_applicable_high_rigidity_cold_beam"
                    : "single_species_equal_weight_isotropic_shell_plus_guide_drift");
  Q019RequireString(
      pin, block, "finite_position_sampling",
      high_rigidity ? "not_applicable_high_rigidity_branch"
                    : (quiet_isotropic_shell_packet ?
                       "one_cell_centered_collocated_nested_haar_octahedral_packets_"
                       "per_cell" :
                       "seeded_random_unique_positions_noise_sensitivity"));
  Q019RequireString(
      pin, block, "finite_gyrophase_position_coupling",
      high_rigidity ? "not_applicable_high_rigidity_cold_beam"
                    : (quiet_isotropic_shell_packet ?
                       "nested_haar_rotated_antipodal_packets_exact_zero_rest_mean_"
                       "isotropic_second_moment_and_local_tsc_current_quiet_start" :
                       "nested_haar_rotated_octahedral_packets_independent_positions"));
  Q019RequireString(
      pin, block, "high_position_sampling",
      high_rigidity ? high_sampling_mode : "not_applicable");
  Q019RequireString(
      pin, block, "finite_packet_grouping_contract",
      high_rigidity ? "not_applicable" :
      (quiet_isotropic_shell_packet ?
       "ppc_divisible_by_six_nested_packets_bound_to_global_root_cell_identity_and_"
       "within_cell_packet_index_and_field_particle_seeds_with_collocated_packet_"
       "spatial_current_quiet_path_proved" :
       "ppc_divisible_by_six_nested_packets_bound_to_global_root_cell_identity_and_"
       "within_cell_packet_index_and_field_particle_seeds_velocity_sequence_only_"
       "independent_positions_not_spatially_current_quiet"));
  Q019RequireString(
      pin, block, "finite_initial_rho_jx_noise_pair_gate",
      high_rigidity ? "not_applicable" :
      "required_quiet_isotropic_shell_vs_independent_positions_first_snapshot_rho_jx");
  Q019RequireString(
      pin, block, "nested_spectrum_contract",
      "box_convergence_uses_identical_shared_spectrum_and_separate_long_mode_control");
  Q019RequireString(
      pin, block, "evolving_resolution_stop_contract",
      "postprocessing_fail_closed_gate_required_runtime_controller_removed_pending_"
      "excluded_pilot_benchmark_and_freeze");
  Q019RequireString(pin, block, "high_rigidity_validity_gate", "measured_not_assumed");
  Q019RequireString(
      pin, block, "high_rigidity_convergence_envelope",
      "staged_2d_ppc_resolution_timestep_and_stochastic_position_controls_required_"
      "before_3d_nonlinear_use");
  Q019RequireString(
      pin, block, "finite_rigidity_predecessor_gate",
      "measured_complex_mode_response_and_convergence_required");
  Q019RequireString(
      pin, block, "current_agreement_gate",
      "deposited_grid_vs_reconstructed_particle_required");
  Q019RequireString(
      pin, block, "morphology_and_relative_drift_diagnostics", "required");
  Q019RequireString(pin, block, "axis_alignment", "B0_parallel_JCR_parallel_x1");
  Q019RequireString(pin, block, "current_normalization",
                    "deposited_j_over_c_equals_2_b_g_k0");
  Q019RequireString(pin, block, "deposit_qscale_semantics",
                    "root_cell_macro_mass_volume_aware");
  Q019RequireString(pin, block, "deposited_prtcl_rho_semantics",
                    "single_species_charge_density");
  Q019RequireString(
      pin, block, "deposited_cr_mass_density_derivation",
      "prtcl_rho_times_species_mass_over_species_charge");
  Q019RequireBoolean(pin, block,
                     "species_mass_and_charge_bound_in_immutable_payload", true);
  Q019RequireString(pin, block, "qualification_effect",
                    "none_source_local_design_only");
  Q019RequireString(pin, block, "numeric_tolerance_status",
                    "machine_scale_moment_tolerance_and_pilot_science_thresholds_unset");
  Q019RequireString(
      pin, block, "no_hall_applicability",
      "bounded_linear_correction_candidate_pending_review_not_accepted");
  Q019RequireString(
      pin, block, "nonlinear_no_hall_applicability_status",
      "fail_closed_pending_registered_local_diagnostics_thresholds_and_external_review");
  Q019RequireString(
      pin, block, "no_hall_mapping_basis",
      q023_carrier ?
      "q023_cr_qom_separate_from_similarity_mapped_background_ion_qom_with_exact_"
      "Bai_R_k0di_Lambda_accounting_and_no_hall_claim_withheld" :
      "equal_similarity_scaled_cr_background_qom_with_exact_Bai_R_k0di_Lambda_"
      "resolved_scale_and_signed_Bai_linear_reduction");
  Q019RequireBoolean(pin, block, "no_hall_applicability_accepted", false);
  Q019RequireBoolean(pin, block,
                     "evolving_local_hall_applicability_diagnostics_complete", false);
  Q019RequireBoolean(pin, block, "bounded_hall_omission_candidate", true);
  Q019RequireBoolean(pin, block, "bounded_hall_omission_review_complete", false);
  Q019RequireBoolean(pin, block, "mhd_resolved_scale_applicability_accepted", false);
  Q019RequireBoolean(pin, block, "R_much_less_than_one_applicability_accepted", false);
  Q019RequireBoolean(pin, block, "no_subion_cell_scale_envelope_satisfied", true);
  Q019RequireBoolean(
      pin, block, "species_q_over_mc_matches_background", !q023_carrier);
  Q019RequireBoolean(pin, block,
                     "positive_finite_loading_accounting_satisfied", true);
  Q019RequireBoolean(pin, block, "resolution_envelope_satisfied_by_design", true);
  Q019RequireBoolean(pin, block, "strong_shock_applicability_authorized", false);
  Q019RequireBoolean(pin, block, "energy_loading_grid_review_complete", false);
  Q019RequireBoolean(pin, block, "universal_saturation_inference_authorized", false);
  Q019RequireBoolean(pin, block, "runtime_box_edge_monitor_installed", true);
  Q019RequireBoolean(pin, block, "runtime_box_edge_monitor_enabled", false);
  Q019RequireBoolean(pin, block, "runtime_box_edge_stop_boundary_frozen", false);
  Q019RequireBoolean(pin, block, "runtime_box_edge_stop_armed", false);
  Q019RequireBoolean(pin, block, "runtime_box_edge_monitor_passive", true);
  Q019RequireString(
      pin, block, "runtime_box_edge_monitor_selection_rule",
      "physical_wave_number_ball_abs_k_le_sqrt_dimension_times_2pi_over_"
      "shortest_active_extent");
  Q019RequireString(
      pin, block, "runtime_box_edge_monitor_geometry",
      "frozen_x1_to_transverse_aspect_ratio_2");
  Q019RequireString(
      pin, block, "runtime_box_edge_monitor_diagnostic_authority",
      "none_quarantined_excluded_pilot_only");
  Q019RequireBoolean(pin, block, "runtime_dominant_scale_stop_controller_installed",
                     false);
  Q019RequireBoolean(pin, block, "raw_production_authorized", false);
  Q019RequireBoolean(pin, block, "nonlinear_saturation_claim_authorized", false);
  Q019RequireBoolean(pin, block, "isolated_variable_claim_authorized", false);
  Q019RequireBoolean(pin, block,
                     "finite_packet_grouping_against_actual_initializer_proved",
                     quiet_isotropic_shell_packet);
  Q019RequireBoolean(pin, block,
                     "finite_density_parallel_current_quiet_by_construction",
                     quiet_isotropic_shell_packet);
  Q019RequireBoolean(pin, block,
                     "finite_density_parallel_current_quiet_measured_and_accepted",
                     false);
  Q019RequireBoolean(pin, block, "launch_authorized", false);
  Q019RequireBoolean(pin, block, "policy_authorized", false);
  Q019RequireBoolean(pin, block, "claim_authorized", false);
  Q019RequireBoolean(pin, block, "physical_pilot_authorized", false);
  Q019RequireString(
      pin, block, "physical_pilot_gate",
      q023_carrier ?
      "explicit_independent_q043_then_registered_q023_then_versioned_resource_"
      "redesign_then_excluded_q023_carrier_pilots" :
      "explicit_independent_q043_then_q023_then_finite_predecessor_source_"
      "compatibility_then_excluded_pilots");
  Q019RequireString(
      pin, block, "q043_independent_raw_cycle_one_oracle_id",
      "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one");
  Q019RequireBoolean(pin, block, "q043_independent_raw_cycle_one_oracle_bound", false);
  Q019RequireString(
      pin, block, "q023_independent_linear_predecessor_id",
      "Q023-PAPER-BELL-LINEAR-JOVERC_after_passed_Q043_oracle");
  Q019RequireBoolean(pin, block, "q023_independent_linear_predecessor_bound", false);
  Q019RequireBoolean(pin, block, "independent_prerequisites_complete", false);
  Q019RequireBoolean(pin, block, "nonlinear_execution_prerequisites_passed", false);
  Q019RequireBoolean(pin, "problem", "user_work_in_loop", true);
  Q019RequireBoolean(pin, "problem", "user_hist", true);
  Q019RequireString(pin, "mesh_refinement", "refinement", "none");
  if (pmy_mesh_->multilevel) {
    Q019NonlinearFatal(
        "Q019 box-edge Fourier monitor requires an unrefined uniform mesh");
  }

  Q019RequireString(pin, "time", "evolution", "dynamic");
  Q019RequireString(pin, "time", "integrator", "vl2");
  if (pin->GetInteger("time", "nlim") == 0 ||
      !(pin->GetReal("time", "tlim") > 0.0)) {
    Q019NonlinearFatal("Q019 candidate decks must genuinely advance");
  }
  Q019RequireString(pin, "mhd", "eos", "ideal");
  Q019RequireString(pin, "particles", "particle_type", "cosmic_ray");
  Q019RequireString(pin, "particles", "pusher", "boris_tsc");
  Q019RequireBoolean(pin, "particles", "deposit_moments", true);
  Q019RequireBoolean(pin, "particles", "couple_moments_to_mhd", true);
  Q019RequireBoolean(pin, "particles", "couple_moments_momentum_to_mhd", true);
  Q019RequireBoolean(pin, "particles", "couple_moments_energy_to_mhd", true);
  Q019RequireString(pin, "particles", "couple_j_to_efield_representation",
                    "cell_centered");
  Q019RequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");
  Q019RequireString(pin, "particles", "couple_fluid_feedback_order", "mhd_src_terms");
  Q019RequireString(pin, "particles", "pic_background_mode", "coupled");
  Q019RequireString(pin, "particles", "pic_feedback_mode", "coupled");
  Q019RequireString(pin, "particles", "pic_interp_scheme", "tsc");
  Q019RequireBoolean(pin, "particles", "pic_enable_2d3v", true);
  Q019RequireString(pin, "particles", "pic_cr_initial_state", "velocity");
  Q019RequireString(pin, "particles", "pic_cr_hall_mode", "off");
  Q019RequireString(pin, "particles", "pic_wave_damping_mode", "off");
  Q019RequireString(pin, "particles", "pic_deltaf_mode", "off");
  Q019RequireString(pin, "particles", "pic_expanding_box_mode", "off");
  if (pin->GetInteger("particles", "deposit_order") != 2) {
    Q019NonlinearFatal("corrected nonlinear Bell successor requires TSC deposition");
  }
  if (pic_max_cell_cross != pin->GetInteger(block, "pic_max_cell_cross")) {
    Q019NonlinearFatal("Q019 particle cell-crossing timestep control drifted");
  }
  Q019RequireClose("pic_theta_max", pin->GetReal("particles", "pic_theta_max"),
                   pin->GetReal(block, "pic_theta_max"));
  Q019RequireClose("immutable Q019 pic_theta_max",
                   pin->GetReal(block, "pic_theta_max"), 0.3);

  const std::string distribution = pin->GetString("particles", "cr_distribution");
  const bool random_distribution = distribution.compare("random") == 0;
  const bool center_distribution = distribution.compare("center") == 0;
  if (!random_distribution && !center_distribution) {
    Q019NonlinearFatal("corrected nonlinear Bell successor distribution is invalid");
  }
  const Real ppc = pin->GetReal("particles", "ppc");
  const int nspecies = pin->GetInteger("particles", "nspecies");
  if (!BranchSamplingIsValid(branch, ppc, nspecies, random_distribution)) {
    Q019NonlinearFatal("branch sampling contract failed; high rigidity requires PPC in "
                       "{1,8,32}, and finite rigidity requires one species, seeded "
                       "random positions, and PPC in {24,48,96}");
  }
  if (high_rigidity &&
      ((high_stochastic_sampling && (!random_distribution || ppc < 8.0)) ||
       (!high_stochastic_sampling && !center_distribution))) {
    Q019NonlinearFatal("high-rigidity position/noise sampling contract failed");
  }
  const int integral_ppc = static_cast<int>(ppc);
  const int expected_spatial_anchors =
      high_rigidity ? 0 : (quiet_isotropic_shell_packet ? 1 : integral_ppc);
  if (pin->GetInteger(
          block, "finite_packet_independent_spatial_anchors_per_cell_expected") !=
      expected_spatial_anchors) {
    Q019NonlinearFatal("finite packet spatial-anchor count drifted");
  }

  const Real rho = pin->GetReal(block, "rho");
  const Real pressure = pin->GetReal(block, "pressure");
  const Real b_g = pin->GetReal(block, "b_g");
  const Real u_a = pin->GetReal(block, "u_a");
  const Real wavelength = pin->GetReal(block, "wavelength");
  const Real k0 = pin->GetReal(block, "k0");
  const Real epsilon = pin->GetReal(block, "epsilon");
  const Real omega = pin->GetReal(block, "omega");
  const Real k0_rg0 = pin->GetReal(block, "k0_rg0");
  const Real rho_cr_over_rho0 = pin->GetReal(block, "rho_cr_over_rho0");
  const Real species_mass = pin->GetReal(block, "species_mass");
  const Real configured_n_cr = pin->GetReal(block, "n_cr");
  const Real eigenmode_amplitude = pin->GetReal(block, "eigenmode_amplitude");
  const Real broadband_amplitude = pin->GetReal(block, "broadband_amplitude");
  const Real separated_long_mode_amplitude =
      pin->GetReal(block, "separated_long_mode_amplitude");
  const Real finite_shell_speed = pin->GetReal(block, "finite_shell_speed");
  const int field_seed = pin->GetInteger(block, "field_seed");
  const int particle_seed = pin->GetInteger("particles", "pic_random_seed");
  const std::string identity_payload = Q019MatrixIdentityPayload(pin, block);
  const std::string identity_fingerprint =
      Q019MatrixIdentityFingerprint(identity_payload);
  const std::string configured_identity_fingerprint =
      pin->GetString(block, "matrix_identity_fingerprint");
  if (configured_identity_fingerprint != identity_fingerprint) {
    Q019NonlinearFatal(
        "immutable runtime deck semantics checksum mismatch: configured=" +
        configured_identity_fingerprint + " measured=" + identity_fingerprint);
  }
  const char *canonical_identity_fingerprint =
      Q019CanonicalMatrixFingerprint(case_id);
  if (canonical_identity_fingerprint == nullptr ||
      identity_fingerprint != canonical_identity_fingerprint) {
    Q019NonlinearFatal("Q019 runtime identity is not an admitted canonical matrix row");
  }
  Q019RequireString(
      pin, block, "runtime_identity_checksum_status",
      "sha256_cryptographic_immutable_runtime_semantics_bound_by_compiled_registry");
  Q019RequireBoolean(pin, block, "external_sha256_execution_receipt_required", true);
  Q019RequireBoolean(pin, block, "external_sha256_execution_receipt_bound", false);
  if (!(rho > 0.0 && pressure > 0.0 && b_g > 0.0 && u_a > 0.0 &&
        wavelength > 0.0 && k0 > 0.0 && epsilon > 0.0 && epsilon < 1.0 &&
        eigenmode_amplitude > 0.0 && eigenmode_amplitude <= 1.0e-2 &&
        broadband_amplitude >= 0.0 && broadband_amplitude <= 1.0e-2 &&
        separated_long_mode_amplitude >= 0.0 &&
        separated_long_mode_amplitude <= 1.0e-2 && rho_cr_over_rho0 > 0.0 &&
        species_mass == 1.0 &&
        field_seed > 0 && particle_seed >= 0)) {
    Q019NonlinearFatal("corrected nonlinear Bell physical or seed contract is invalid");
  }
  if (!BranchRigidityIsValid(branch, k0_rg0)) {
    Q019NonlinearFatal("branch rigidity separator failed");
  }
  if ((high_rigidity &&
       (finite_shell_speed != 0.0 ||
        (high_stochastic_sampling ? particle_seed <= 0 : particle_seed != 0))) ||
      (!high_rigidity && (finite_shell_speed <= 0.0 || particle_seed <= 0))) {
    Q019NonlinearFatal("finite-rigidity distribution parameters leaked across branches");
  }
  if ((seed_topology == SeedTopology::shared_spectrum_only &&
       separated_long_mode_amplitude != 0.0) ||
      (seed_topology == SeedTopology::shared_spectrum_plus_separated_long_modes &&
       separated_long_mode_amplitude <= 0.0)) {
    Q019NonlinearFatal("separated long-mode amplitude does not match seed topology");
  }
  if (!(rho > 0.0 && pressure > 0.0 && b_g > 0.0 && u_a > 0.0)) {
    Q019NonlinearFatal("Q019 background state must be finite and positive");
  }
  Q019RequireClose("pressure", pressure, 1.0);
  Q019RequireClose("u_a", u_a, b_g/std::sqrt(rho));
  Q019RequireClose("k0", k0, 2.0*M_PI/wavelength);
  const Real guide_parallel_stream_speed = u_a/epsilon;
  if (!(guide_parallel_stream_speed > 0.0)) {
    Q019NonlinearFatal("guide-parallel stream speed is invalid");
  }

  const Real x1_extent = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const Real x2_extent = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  const Real x3_extent = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  for (const auto value : {x1_extent/wavelength, x2_extent/wavelength,
                           x3_extent/wavelength}) {
    Q019RequireClose("domain extent in seed wavelengths", value, std::round(value));
  }
  if (!BoxEdgeFrozenAspectRatioIsValid(
          dimension, x1_extent, x2_extent, x3_extent)) {
    Q019NonlinearFatal(
        "corrected nonlinear Bell box-edge monitor requires frozen 2:1 "
        "x1-to-transverse geometry");
  }
  const int expected_box_edge_modes = BoxEdgeUniqueModeCount(dimension);
  if (pin->GetInteger(block, "runtime_box_edge_monitor_unique_mode_count") !=
          expected_box_edge_modes ||
      pin->GetInteger(block, "runtime_box_edge_monitor_cell_passes_per_sample") !=
          1 + BoxEdgeReductionGroupCount(dimension) ||
      pin->GetInteger(block, "runtime_box_edge_monitor_global_reductions_per_sample") !=
          1 + BoxEdgeReductionGroupCount(dimension)) {
    Q019NonlinearFatal("Q019 physical-low-k monitor cost contract drifted");
  }
  for (int n = 0; n < expected_box_edge_modes; ++n) {
    if (!BoxEdgePhysicalLowKModeIsSelected(
            BoxEdgeUniqueModeAt(n, dimension), dimension,
            x1_extent, x2_extent, x3_extent)) {
      Q019NonlinearFatal("Q019 physical-low-k mode enumeration drifted");
    }
  }
  const Real max_active_dx = dimension == 3 ?
      std::max({pmy_mesh_->mesh_size.dx1, pmy_mesh_->mesh_size.dx2,
                pmy_mesh_->mesh_size.dx3}) :
      std::max(pmy_mesh_->mesh_size.dx1, pmy_mesh_->mesh_size.dx2);
  const Real rigidity_momentum_per_mass =
      high_rigidity ? guide_parallel_stream_speed : finite_shell_speed;
  const Real nominal_rg0 =
      NominalRigidityLength(rigidity_momentum_per_mass, omega);
  Q019RequireClose("k0*r_g0", k0*nominal_rg0, k0_rg0);
  Q019RequireClose("nominal r_g0", pin->GetReal(block, "nominal_rg0"), nominal_rg0);
  Q019RequireClose("initial nominal r_g0/dx",
                   pin->GetReal(block, "initial_nominal_rg0_over_max_active_dx"),
                   nominal_rg0/max_active_dx);
  const Real minimum_evolving_rl_over_dx =
      pin->GetReal(block, "minimum_characteristic_shell_rl_over_dx");
  Q019RequireClose("immutable Q019 characteristic shell rL/dx floor",
                   minimum_evolving_rl_over_dx, 8.0);
  const Real required_initial_rg0_over_dx =
      pin->GetReal(block, "required_initial_rg0_over_max_active_dx");
  const bool three_d_onset_family = !high_rigidity && dimension == 3;
  const Real preregistered_onset_bperp_rms_over_b0 =
      pin->GetReal(block, "preregistered_nonlinear_onset_Bperp_rms_over_B0");
  const Real resolution_design_maximum_sampled_b_over_b0 =
      pin->GetReal(block, "resolution_design_maximum_sampled_B_over_B0");
  Q019RequireClose(
      "preregistered nonlinear-onset Bperp_rms/B0",
      preregistered_onset_bperp_rms_over_b0, three_d_onset_family ? 1.0 : 0.0);
  Q019RequireClose(
      "resolution-design maximum sampled B/B0",
      resolution_design_maximum_sampled_b_over_b0, three_d_onset_family ? 2.0 : 0.0);
  Q019RequireBoolean(
      pin, block, "common_nonlinear_onset_resolution_reachable", false);
  Q019RequireString(
      pin, block, "common_nonlinear_onset_reachability_status",
      three_d_onset_family ?
      "pilot_pending_not_inferred_from_resolution_design_envelope" :
      "not_applicable");
  Q019RequireBoolean(
      pin, block, "coarse_resolution_stop_or_intermittency_limitation",
      role == "finite_rigidity_3d_onset_matched_convergence_control" &&
      case_id == "q019-fr-3d-onset-small-resolution-coarse-s0");
  Q019RequireBoolean(pin, block, "runtime_resolution_stop_controller_installed", false);
  Q019RequireBoolean(pin, block, "runtime_resolution_guard_pilot_qualified", false);
  Q019RequireBoolean(
      pin, block, "postprocessing_resolution_gate_required", !high_rigidity);
  Q019RequireString(
      pin, block, "runtime_resolution_stop_controller_status",
      "not_installed_unbounded_per_cycle_scan_removed_pending_excluded_pilot_"
      "benchmark_and_freeze");
  const Real expected_required_initial_rg0_over_dx =
      high_rigidity ? 0.0 : minimum_evolving_rl_over_dx;
  Q019RequireClose("required initial r_g0/dx",
                   required_initial_rg0_over_dx,
                   expected_required_initial_rg0_over_dx);
  if (!high_rigidity &&
      !ResolutionEnvelopeIsValid(nominal_rg0/max_active_dx,
                                 required_initial_rg0_over_dx)) {
    Q019NonlinearFatal("finite-rigidity resolution stop/acceptance envelope failed");
  }
  const Real maximum_sampled_b_over_b0_before_resolution_stop =
      high_rigidity ? 0.0 :
      (nominal_rg0/max_active_dx)/minimum_evolving_rl_over_dx;
  Q019RequireClose(
      "maximum sampled B/B0 before characteristic resolution stop",
      pin->GetReal(block, "maximum_sampled_B_over_B0_before_resolution_stop"),
      maximum_sampled_b_over_b0_before_resolution_stop);
  Q019RequireClose("cr_vx0", pin->GetReal("particles", "cr_vx0"),
                   guide_parallel_stream_speed);
  Q019RequireClose("cr_vy0", pin->GetReal("particles", "cr_vy0"), 0.0);
  Q019RequireClose("cr_vz0", pin->GetReal("particles", "cr_vz0"), 0.0);
  const Real light_speed = pin->GetReal("particles", "pic_cr_light_speed");
  if (light_speed <= guide_parallel_stream_speed + finite_shell_speed) {
    Q019NonlinearFatal("artificial CR light speed must exceed every initial CR speed");
  }
  Q019RequireClose(
      "CR RMS-speed kinetic-loading proxy/background magnetic energy",
      CRRMSSpeedKineticLoadingProxyToBackgroundMagneticEnergy(
          rho_cr_over_rho0, rho,
          std::sqrt(guide_parallel_stream_speed*guide_parallel_stream_speed +
                    finite_shell_speed*finite_shell_speed),
          light_speed, b_g),
      pin->GetReal(
          block,
          "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy"));

  const Real species_charge = pin->GetReal("species0", "charge");
  Q019RequireClose("species0/charge", species_charge, omega/b_g);
  Vector3 summed_species_charge_times_velocity = {0.0, 0.0, 0.0};
  for (int species = 0; species < nspecies; ++species) {
    const std::string species_block = "species" + std::to_string(species);
    Q019RequireClose(species_block + "/mass", pin->GetReal(species_block, "mass"), 1.0);
    Q019RequireClose(species_block + "/charge", pin->GetReal(species_block, "charge"),
                     species_charge);
    Q019RequireClose(species_block + "/vx0", pin->GetReal(species_block, "vx0"),
                     guide_parallel_stream_speed);
    Q019RequireClose(species_block + "/vy0", pin->GetReal(species_block, "vy0"), 0.0);
    Q019RequireClose(species_block + "/vz0", pin->GetReal(species_block, "vz0"), 0.0);
    summed_species_charge_times_velocity = Add(
        summed_species_charge_times_velocity,
        {species_charge*pin->GetReal(species_block, "vx0"),
         species_charge*pin->GetReal(species_block, "vy0"),
         species_charge*pin->GetReal(species_block, "vz0")});
  }

  const Real qscale = pin->GetReal("particles", "deposit_qscale");
  const Real root_cell_volume = RootCellVolume(
      x1_extent, pmy_mesh_->mesh_indcs.nx1, x2_extent, pmy_mesh_->mesh_indcs.nx2,
      x3_extent, pmy_mesh_->mesh_indcs.nx3);
  Q019RequireClose("root cell volume", root_cell_volume,
                   pmy_mesh_->mesh_size.dx1*pmy_mesh_->mesh_size.dx2*
                   pmy_mesh_->mesh_size.dx3);
  Q019RequireClose(
      "PPC*deposit_qscale*species_charge*v_CR/V_root_cell",
      DepositedJOverC(ppc, qscale, species_charge, guide_parallel_stream_speed,
                      root_cell_volume),
      RequiredDepositedJOverC(b_g, k0));
  const Real measured_rho_cr_over_rho0 =
      CRMassDensityOverRho0(ppc, qscale, root_cell_volume, rho);
  Q019RequireClose("derived rho_CR/rho0", measured_rho_cr_over_rho0,
                   rho_cr_over_rho0);
  Q019RequireClose("derived n_CR", CRNumberDensity(rho_cr_over_rho0*rho, species_mass),
                   configured_n_cr);
  const Real background_q_over_mc_reference =
      pin->GetReal(block, "background_q_over_mc_reference");
  if (q023_carrier) {
    if (!(species_charge < background_q_over_mc_reference)) {
      Q019NonlinearFatal(
          "Q019 Q023 carrier requires CR q/(mc) below the background-ion mapping");
    }
  } else {
    Q019RequireClose("equal similarity-scaled ion CR/background q/(mc)",
                     species_charge, background_q_over_mc_reference);
  }
  Q019RequireClose("declared similarity-scaled ion q/(mc)",
                   background_q_over_mc_reference, 10000.0);
  const Real background_ion_gyrofrequency =
      BackgroundIonGyrofrequency(background_q_over_mc_reference, b_g);
  const Real background_ion_inertial_length =
      BackgroundIonInertialLength(u_a, background_ion_gyrofrequency);
  Q019RequireClose("background-ion gyrofrequency", background_ion_gyrofrequency,
                   pin->GetReal(block, "background_ion_gyrofrequency"));
  Q019RequireClose("background-ion inertial length", background_ion_inertial_length,
                   pin->GetReal(block, "background_ion_inertial_length"));
  Q019RequireClose("k0*d_i", k0*background_ion_inertial_length,
                   pin->GetReal(block, "k0_background_ion_inertial_length"));
  const Real minimum_active_dx = dimension == 3 ?
      std::min({pmy_mesh_->mesh_size.dx1, pmy_mesh_->mesh_size.dx2,
                pmy_mesh_->mesh_size.dx3}) :
      std::min(pmy_mesh_->mesh_size.dx1, pmy_mesh_->mesh_size.dx2);
  Q019RequireClose("minimum active dx/d_i",
                   minimum_active_dx/background_ion_inertial_length,
                   pin->GetReal(block, "minimum_active_dx_over_background_di"));
  if (!NoSubIonCellScaleEnvelopeIsValid(minimum_active_dx,
                                        background_ion_inertial_length)) {
    Q019NonlinearFatal("Q019 mapping resolves sub-ion cell scales and is prohibited");
  }
  const Real charge_density_ratio = ChargeDensityRatio(
      rho_cr_over_rho0, species_charge, background_q_over_mc_reference);
  const Real hall_parameter =
      HallParameter(charge_density_ratio, guide_parallel_stream_speed, u_a);
  const Real hall_parameter_from_current = HallParameterFromCurrent(
      RequiredDepositedJOverC(b_g, k0), rho, background_q_over_mc_reference, u_a,
      charge_density_ratio);
  Q019RequireClose("exact Bai charge-density ratio R", charge_density_ratio,
                   pin->GetReal(block, "charge_density_ratio_equal_background_qom"));
  Q019RequireClose("exact Bai Hall parameter Lambda", hall_parameter,
                   pin->GetReal(block, "hall_parameter_equal_background_qom"));
  Q019RequireClose("Hall parameter from current", hall_parameter_from_current,
                   pin->GetReal(block, "hall_parameter_from_current"));
  Q019RequireClose("Hall parameter identity", hall_parameter,
                   hall_parameter_from_current);
  Q019RequireClose("Lambda_Hall/(2 k0 d_i)",
                   hall_parameter_from_current/
                       (2.0*k0*background_ion_inertial_length),
                   pin->GetReal(block, "hall_parameter_over_twice_k0_di"));
  Q019RequireClose("background q/(mc) at Lambda=1 order-unity reference",
                   pin->GetReal(block,
                                "background_q_over_mc_at_lambda_equal_one_reference"),
                   RequiredDepositedJOverC(b_g, k0)/(rho*u_a)*
                       (1.0 - charge_density_ratio));
  const Real hall_order_unity_reference =
      pin->GetReal(block, "hall_order_unity_reference");
  Q019RequireClose("Hall order-unity reference",
                   hall_order_unity_reference, 1.0);
  Q019RequireClose("Hall order-unity reference margin",
                   hall_order_unity_reference/hall_parameter_from_current,
                   pin->GetReal(block, "hall_order_unity_reference_margin"));
  Q019RequireClose("Bai Hall linear factor",
                   BaiHallLinearFactor(hall_parameter_from_current),
                   pin->GetReal(block, "bai_hall_linear_factor"));
  Q019RequireClose("Bai Hall growth-rate fractional shift",
                   BaiHallGrowthRateFractionalShift(hall_parameter_from_current),
                   pin->GetReal(block, "bai_hall_growth_rate_fractional_shift"));
  Q019RequireClose("Bai Hall wavenumber fractional shift",
                   BaiHallWavenumberFractionalShift(hall_parameter_from_current),
                   pin->GetReal(block, "bai_hall_wavenumber_fractional_shift"));
  Q019RequireClose("Bai Hall growth-rate reduction factor",
                   BaiHallGrowthRateReductionFactor(hall_parameter_from_current),
                   pin->GetReal(block, "bai_hall_growth_rate_reduction_factor"));
  Q019RequireClose("Bai Hall wavenumber reduction factor",
                   BaiHallWavenumberReductionFactor(hall_parameter_from_current),
                   pin->GetReal(block, "bai_hall_wavenumber_reduction_factor"));
  const Real maximum_charge_density_ratio =
      q023_carrier ? static_cast<Real>(1.0e-3) : static_cast<Real>(1.0e-5);
  if (!(charge_density_ratio > 0.0 &&
        charge_density_ratio <= maximum_charge_density_ratio &&
        hall_parameter_from_current > 0.0 &&
        hall_parameter_from_current < hall_order_unity_reference)) {
    Q019NonlinearFatal(
        "Q019 explicit R/Lambda similarity-mapping design envelope failed");
  }
  const Real inertia_parameter = CRInertiaParameter(rho_cr_over_rho0);
  const Real momentum_loading = CRMomentumLoadingParameter(
      rho_cr_over_rho0, guide_parallel_stream_speed, u_a);
  Q019RequireClose("CR inertia parameter", inertia_parameter,
                   pin->GetReal(block, "cr_inertia_parameter"));
  Q019RequireClose("CR momentum loading", momentum_loading,
                   pin->GetReal(block, "cr_momentum_loading_parameter"));
  Q019RequireClose("feedback force parameter",
                   FeedbackForceParameter(RequiredDepositedJOverC(b_g, k0), b_g,
                                          rho, u_a, k0),
                   pin->GetReal(block, "feedback_force_parameter"));
  if (!LoadingAccountingIsFiniteAndPositive(inertia_parameter, momentum_loading)) {
    Q019NonlinearFatal("CR inertia or momentum-loading accounting is invalid");
  }
  const Vector3 configured_species_sum_j_over_c = DepositedSpeciesSumJOverC(
      ppc, nspecies, qscale, summed_species_charge_times_velocity, root_cell_volume);
  Q019RequireClose("species-summed guide-parallel deposited J_CR/c",
                   configured_species_sum_j_over_c.x1,
                   RequiredDepositedJOverC(b_g, k0));
  Q019RequireClose("species-summed transverse deposited Jy/c",
                   configured_species_sum_j_over_c.x2, 0.0);
  Q019RequireClose("species-summed transverse deposited Jz/c",
                   configured_species_sum_j_over_c.x3, 0.0);
  if (!high_rigidity) {
    const Real first_moment_scale = finite_shell_speed;
    const Real second_moment_scale = finite_shell_speed*finite_shell_speed;
    const std::uint64_t validation_global_cell_id = 0;
    for (int packet = 0; packet < integral_ppc/6; ++packet) {
      Vector3 centered_sum = {0.0, 0.0, 0.0};
      Real centered_second[3][3] = {};
      for (int direction = 0; direction < 6; ++direction) {
        const int sample = packet*6 + direction;
        const Vector3 velocity = NestedHaarOctahedralShellVelocityAtSample(
            integral_ppc, sample, validation_global_cell_id, field_seed, particle_seed,
            guide_parallel_stream_speed, finite_shell_speed);
        const Vector3 centered = {
          velocity.x1 - guide_parallel_stream_speed, velocity.x2, velocity.x3};
        centered_sum = Add(centered_sum, centered);
        const Real values[3] = {
          static_cast<Real>(centered.x1),
          static_cast<Real>(centered.x2),
          static_cast<Real>(centered.x3)
        };
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            centered_second[a][b] += values[a]*values[b]/6.0;
          }
        }
      }
      Q019RequireMomentClose(
          "finite-shell packet centered mean x1", centered_sum.x1/6.0, 0.0,
          first_moment_scale);
      Q019RequireMomentClose(
          "finite-shell packet centered mean x2", centered_sum.x2/6.0, 0.0,
          first_moment_scale);
      Q019RequireMomentClose(
          "finite-shell packet centered mean x3", centered_sum.x3/6.0, 0.0,
          first_moment_scale);
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          Q019RequireMomentClose(
              "finite-shell packet centered isotropic second moment",
              centered_second[a][b],
              a == b ? second_moment_scale/3.0 : 0.0, second_moment_scale);
        }
      }
    }
  }

  Q019RequireClose("nonrelativistic characteristic shell p_iso/m",
                   pin->GetReal(block, "characteristic_shell_p_iso_over_m"),
                   finite_shell_speed);
  q019_runtime_dimension = dimension;
  q019_runtime_pin = pin;
  q019_runtime_background_b = b_g;
  q019_runtime_controller_enabled = pin->DoesParameterExist(
      kQ019RuntimeControllerBlock, "schema");
  q019_runtime_controller_box_edge_monitor_enabled = false;
  q019_runtime_resolution_monitor_enabled = false;
  q019_runtime_diagnostic_failure_stop_armed = false;
  q019_runtime_resolution_stop_armed = false;
  q019_runtime_controller_box_edge_stop_armed = false;
  q019_runtime_controller_triggered = false;
  q019_runtime_controller_trigger_failure = false;
  q019_runtime_controller_trigger_reason = 0;
  q019_runtime_controller_trigger_cycle = 0;
  q019_runtime_resolution_samples = 0;
  q019_runtime_resolution_last_cycle = 0;
  q019_runtime_controller_box_edge_stop_ppm = -1;
  q019_runtime_controller_pilot_cycle_limit = -1;
  q019_runtime_controller_monitor_dt = 0.0;
  q019_runtime_resolution_stop_b_over_b0 = -1.0;
  q019_runtime_resolution_last_time = 0.0;
  q019_runtime_resolution_last_b_over_b0 = -1.0;
  q019_runtime_resolution_max_b_over_b0 = -1.0;
  q019_runtime_controller_trigger_time = 0.0;
  q019_runtime_controller_trigger_metric = -1.0;
  if (q019_runtime_controller_enabled) {
    const std::string controller = kQ019RuntimeControllerBlock;
    if (pin->GetInteger(controller, "schema") != kQ019RuntimeControllerSchema) {
      Q019NonlinearFatal("Q019 runtime controller schema drifted");
    }
    const std::string controller_authority =
        pin->GetString(controller, "authority");
    const bool runtime_regression =
        controller_authority == "runtime_regression_only";
    const bool excluded_pilot =
        controller_authority == "excluded_pilot_only";
    if (!runtime_regression && !excluded_pilot) {
      Q019NonlinearFatal("Q019 runtime controller authority is invalid");
    }
    Q019RequireString(
        pin, controller, "contract_id",
        runtime_regression ? "q019-runtime-controller-v1-regression" :
                             "q019-runtime-controller-v1-excluded-pilot");
    Q019RequireString(
        pin, controller, "cadence_status",
        runtime_regression ? "accelerated_runtime_regression" :
                             "excluded_pilot_candidate_not_frozen");
    Q019RequireString(pin, controller, "stop_disposition",
                      "accepted_guard_stop_or_diagnostic_failure");
    Q019RequireString(pin, controller, "source_case_id", case_id);
    const std::string controller_identity = Q019MatrixIdentityFingerprint(
        Q019RuntimeControllerIdentityPayload(pin));
    if (pin->GetString(controller, "controller_identity_fingerprint") !=
        controller_identity) {
      Q019NonlinearFatal("Q019 runtime controller identity fingerprint drifted");
    }
    q019_runtime_resolution_monitor_enabled =
        pin->GetBoolean(controller, "resolution_monitor_enabled");
    q019_runtime_controller_box_edge_monitor_enabled =
        pin->GetBoolean(controller, "box_edge_monitor_enabled");
    q019_runtime_diagnostic_failure_stop_armed =
        pin->GetBoolean(controller, "diagnostic_failure_stop_armed");
    q019_runtime_resolution_stop_armed =
        pin->GetBoolean(controller, "resolution_stop_armed");
    q019_runtime_controller_box_edge_stop_armed =
        pin->GetBoolean(controller, "box_edge_stop_armed");
    q019_runtime_resolution_stop_b_over_b0 =
        pin->GetReal(controller, "resolution_stop_B_over_B0");
    q019_runtime_controller_box_edge_stop_ppm =
        pin->GetInteger(controller, "box_edge_stop_ppm");
    q019_runtime_controller_monitor_dt =
        pin->GetReal(controller, "monitor_dt");
    q019_runtime_controller_pilot_cycle_limit =
        pin->GetInteger(controller, "pilot_cycle_limit");
    if (q019_runtime_resolution_monitor_enabled !=
            q019_runtime_controller_box_edge_monitor_enabled ||
        (q019_runtime_resolution_stop_armed &&
         !q019_runtime_resolution_monitor_enabled) ||
        (q019_runtime_controller_box_edge_stop_armed &&
         !q019_runtime_controller_box_edge_monitor_enabled) ||
        (q019_runtime_diagnostic_failure_stop_armed &&
         !q019_runtime_resolution_monitor_enabled) ||
        ((q019_runtime_resolution_stop_armed ||
          q019_runtime_controller_box_edge_stop_armed) &&
         !q019_runtime_diagnostic_failure_stop_armed) ||
        !(std::isfinite(q019_runtime_controller_monitor_dt) &&
          q019_runtime_controller_monitor_dt > 0.0) ||
        (q019_runtime_resolution_stop_armed &&
         !(std::isfinite(q019_runtime_resolution_stop_b_over_b0) &&
           q019_runtime_resolution_stop_b_over_b0 > 1.0)) ||
        (!q019_runtime_resolution_stop_armed &&
         q019_runtime_resolution_stop_b_over_b0 != -1.0) ||
        (q019_runtime_controller_box_edge_stop_armed &&
         (q019_runtime_controller_box_edge_stop_ppm <= 0 ||
          q019_runtime_controller_box_edge_stop_ppm > 1000000)) ||
        (!q019_runtime_controller_box_edge_stop_armed &&
         q019_runtime_controller_box_edge_stop_ppm != -1) ||
        q019_runtime_controller_pilot_cycle_limit == 0 ||
        q019_runtime_controller_pilot_cycle_limit < -1 ||
        (runtime_regression &&
         case_id != "q019-fr-runtime-initializer-ppc24-s0")) {
      Q019NonlinearFatal("Q019 runtime controller configuration is invalid");
    }
    if (q019_runtime_resolution_stop_armed) {
      Q019RequireClose(
          "Q019 excluded-pilot resolution stop threshold",
          q019_runtime_resolution_stop_b_over_b0,
          pin->GetReal(block, "maximum_sampled_B_over_B0_before_resolution_stop"));
    }
    q019_runtime_resolution_samples = pin->GetOrAddInteger(
        controller, "runtime_resolution_samples", 0);
    q019_runtime_resolution_last_cycle = pin->GetOrAddInteger(
        controller, "runtime_resolution_last_cycle", 0);
    q019_runtime_resolution_last_time = pin->GetOrAddReal(
        controller, "runtime_resolution_last_time", 0.0);
    q019_runtime_resolution_last_b_over_b0 = pin->GetOrAddReal(
        controller, "runtime_resolution_last_B_over_B0", -1.0);
    q019_runtime_resolution_max_b_over_b0 = pin->GetOrAddReal(
        controller, "runtime_resolution_max_B_over_B0", -1.0);
    q019_runtime_controller_triggered = pin->GetOrAddBoolean(
        controller, "runtime_controller_triggered", false);
    q019_runtime_controller_trigger_failure = pin->GetOrAddBoolean(
        controller, "runtime_controller_trigger_failure", false);
    q019_runtime_controller_trigger_reason = pin->GetOrAddInteger(
        controller, "runtime_controller_trigger_reason", 0);
    q019_runtime_controller_trigger_cycle = pin->GetOrAddInteger(
        controller, "runtime_controller_trigger_cycle", 0);
    q019_runtime_controller_trigger_time = pin->GetOrAddReal(
        controller, "runtime_controller_trigger_time", 0.0);
    q019_runtime_controller_trigger_metric = pin->GetOrAddReal(
        controller, "runtime_controller_trigger_metric", -1.0);
    const bool empty_resolution_state =
        q019_runtime_resolution_samples == 0 &&
        q019_runtime_resolution_last_cycle == 0 &&
        q019_runtime_resolution_last_time == 0.0 &&
        q019_runtime_resolution_last_b_over_b0 == -1.0 &&
        q019_runtime_resolution_max_b_over_b0 == -1.0;
    const Real restart_time_tolerance =
        static_cast<Real>(512.0)*std::numeric_limits<Real>::epsilon()*
        std::max(static_cast<Real>(1.0), std::abs(pmy_mesh_->time));
    const bool sampled_resolution_state =
        q019_runtime_resolution_samples > 0 &&
        q019_runtime_resolution_last_cycle > 0 &&
        q019_runtime_resolution_last_cycle <= pmy_mesh_->ncycle &&
        std::isfinite(q019_runtime_resolution_last_time) &&
        q019_runtime_resolution_last_time >= 0.0 &&
        q019_runtime_resolution_last_time <=
            pmy_mesh_->time + restart_time_tolerance &&
        std::isfinite(q019_runtime_resolution_last_b_over_b0) &&
        q019_runtime_resolution_last_b_over_b0 > 0.0 &&
        std::isfinite(q019_runtime_resolution_max_b_over_b0) &&
        q019_runtime_resolution_max_b_over_b0 >=
            q019_runtime_resolution_last_b_over_b0;
    const bool empty_trigger_state =
        !q019_runtime_controller_triggered &&
        !q019_runtime_controller_trigger_failure &&
        q019_runtime_controller_trigger_reason == 0 &&
        q019_runtime_controller_trigger_cycle == 0 &&
        q019_runtime_controller_trigger_time == 0.0 &&
        q019_runtime_controller_trigger_metric == -1.0;
    const bool trigger_reason_is_failure =
        q019_runtime_controller_trigger_reason == kQ019DiagnosticFailureReason;
    const bool trigger_reason_is_guard =
        q019_runtime_controller_trigger_reason == kQ019ResolutionStopReason ||
        q019_runtime_controller_trigger_reason == kQ019BoxEdgeStopReason ||
        q019_runtime_controller_trigger_reason == kQ019PilotCycleStopReason;
    const bool observed_trigger_state =
        q019_runtime_controller_triggered &&
        (trigger_reason_is_failure || trigger_reason_is_guard) &&
        q019_runtime_controller_trigger_failure == trigger_reason_is_failure &&
        q019_runtime_controller_trigger_cycle > 0 &&
        q019_runtime_controller_trigger_cycle <= pmy_mesh_->ncycle &&
        std::isfinite(q019_runtime_controller_trigger_time) &&
        q019_runtime_controller_trigger_time >= 0.0 &&
        q019_runtime_controller_trigger_time <=
            pmy_mesh_->time + restart_time_tolerance &&
        std::isfinite(q019_runtime_controller_trigger_metric);
    const bool resolution_state =
        empty_resolution_state || sampled_resolution_state;
    const bool trigger_state = empty_trigger_state || observed_trigger_state;
    if ((restart && !(resolution_state && trigger_state)) ||
        (!restart && !(empty_resolution_state && empty_trigger_state))) {
      Q019NonlinearFatal("Q019 runtime controller restart state drifted");
    }
  }
  q019_runtime_box_edge_monitor_enabled =
      pin->GetBoolean(block, "runtime_box_edge_monitor_enabled") ||
      q019_runtime_controller_box_edge_monitor_enabled;
  q019_runtime_box_edge_stop_armed =
      pin->GetBoolean(block, "runtime_box_edge_stop_armed");
  q019_runtime_box_edge_monitor_dt =
      q019_runtime_controller_enabled ? q019_runtime_controller_monitor_dt :
      pin->GetReal(block, "runtime_box_edge_monitor_dt");
  q019_runtime_box_edge_stop_ppm =
      pin->GetInteger(block, "runtime_box_edge_stop_ppm");
  if (!(std::isfinite(q019_runtime_box_edge_monitor_dt) &&
        q019_runtime_box_edge_monitor_dt > 0.0)) {
    Q019NonlinearFatal("Q019 box-edge monitor cadence is invalid");
  }
  if (!q019_runtime_controller_enabled) {
    Q019RequireClose("immutable Q019 box-edge monitor cadence",
                     q019_runtime_box_edge_monitor_dt, 0.1);
  }
  if (q019_runtime_box_edge_stop_armed ||
      q019_runtime_box_edge_stop_ppm != -1) {
    Q019NonlinearFatal(
      "Q019 source-local matrix cannot arm an unfrozen box-edge stop threshold");
  }
  if (restart) {
    if (pin->GetInteger(block, "runtime_box_edge_monitor_schema") != 2) {
      Q019NonlinearFatal("Q019 box-edge monitor restart schema drifted");
    }
    q019_runtime_box_edge_next_nominal_time =
        pin->GetReal(block, "runtime_box_edge_monitor_next_nominal_time");
    q019_runtime_box_edge_last_nominal_time =
        pin->GetReal(block, "runtime_box_edge_monitor_last_nominal_time");
    q019_runtime_box_edge_last_cycle =
        pin->GetInteger(block, "runtime_box_edge_monitor_last_cycle");
    q019_runtime_box_edge_completed_slots =
        pin->GetInteger(block, "runtime_box_edge_monitor_completed_slots");
    q019_runtime_box_edge_valid_samples =
        pin->GetInteger(block, "runtime_box_edge_monitor_valid_samples");
    q019_runtime_box_edge_skipped_slots =
        pin->GetInteger(block, "runtime_box_edge_monitor_skipped_slots");
    q019_runtime_box_edge_last_status =
        pin->GetInteger(block, "runtime_box_edge_monitor_last_status");
    q019_runtime_box_edge_status_mask =
        pin->GetInteger(block, "runtime_box_edge_monitor_status_mask");
    q019_runtime_box_edge_last_prior_time =
        pin->GetReal(block, "runtime_box_edge_monitor_last_prior_time");
    q019_runtime_box_edge_last_time =
        pin->GetReal(block, "runtime_box_edge_monitor_last_time");
    q019_runtime_box_edge_last_power_fraction =
        pin->GetReal(block, "runtime_box_edge_monitor_last_power_fraction");
    q019_runtime_box_edge_max_power_fraction =
        pin->GetReal(block, "runtime_box_edge_monitor_max_power_fraction");
    q019_runtime_box_edge_last_fluctuation_mean =
        pin->GetReal(block, "runtime_box_edge_monitor_last_fluctuation_mean");
    const Real chronology_tolerance = static_cast<Real>(128.0)*
        std::numeric_limits<Real>::epsilon()*
        std::max(static_cast<Real>(1.0), std::abs(pmy_mesh_->time));
    const Real expected_next_nominal_time = static_cast<Real>(
        BoxEdgeNextNominalTime(q019_runtime_box_edge_completed_slots,
                               q019_runtime_box_edge_monitor_dt));
    const Real expected_last_nominal_time = static_cast<Real>(
        BoxEdgeLastNominalTime(q019_runtime_box_edge_completed_slots,
                               q019_runtime_box_edge_monitor_dt));
    const bool schedule_state =
        q019_runtime_box_edge_completed_slots >= 0 &&
        q019_runtime_box_edge_valid_samples >= 0 &&
        q019_runtime_box_edge_valid_samples <=
            q019_runtime_box_edge_completed_slots &&
        q019_runtime_box_edge_skipped_slots >= 0 &&
        q019_runtime_box_edge_skipped_slots <=
            q019_runtime_box_edge_completed_slots -
                q019_runtime_box_edge_valid_samples &&
        BoxEdgeStatusCodeIsValid(q019_runtime_box_edge_last_status) &&
        q019_runtime_box_edge_status_mask >= 0 &&
        std::isfinite(q019_runtime_box_edge_next_nominal_time) &&
        std::isfinite(q019_runtime_box_edge_last_nominal_time) &&
        std::abs(q019_runtime_box_edge_next_nominal_time -
                 expected_next_nominal_time) <= chronology_tolerance &&
        std::abs(q019_runtime_box_edge_last_nominal_time -
                 expected_last_nominal_time) <= chronology_tolerance;
    const bool unobserved_payload =
        q019_runtime_box_edge_completed_slots == 0 &&
        q019_runtime_box_edge_valid_samples == 0 &&
        q019_runtime_box_edge_skipped_slots == 0 &&
        q019_runtime_box_edge_last_status ==
            static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled) &&
        q019_runtime_box_edge_status_mask == 0 &&
        q019_runtime_box_edge_last_cycle == 0 &&
        q019_runtime_box_edge_last_nominal_time == -1.0 &&
        q019_runtime_box_edge_last_prior_time == 0.0 &&
        q019_runtime_box_edge_last_time == 0.0 &&
        q019_runtime_box_edge_last_power_fraction == -1.0 &&
        q019_runtime_box_edge_max_power_fraction == -1.0 &&
        q019_runtime_box_edge_last_fluctuation_mean == -1.0;
    const bool unobserved_state =
        unobserved_payload &&
        (!q019_runtime_box_edge_monitor_enabled ||
         pmy_mesh_->time + chronology_tolerance <
             q019_runtime_box_edge_next_nominal_time);
    const bool maximum_state =
        (q019_runtime_box_edge_valid_samples == 0 &&
         q019_runtime_box_edge_max_power_fraction == -1.0) ||
        (q019_runtime_box_edge_valid_samples > 0 &&
         std::isfinite(q019_runtime_box_edge_max_power_fraction) &&
         q019_runtime_box_edge_max_power_fraction >= 0.0 &&
         q019_runtime_box_edge_max_power_fraction <= 1.0);
    const bool valid_metric_state =
        q019_runtime_box_edge_last_status ==
            static_cast<int>(BoxEdgeDiagnosticStatus::valid) &&
        q019_runtime_box_edge_valid_samples > 0 &&
        std::isfinite(q019_runtime_box_edge_last_power_fraction) &&
        q019_runtime_box_edge_last_power_fraction >= 0.0 &&
        q019_runtime_box_edge_last_power_fraction <= 1.0 &&
        q019_runtime_box_edge_max_power_fraction >=
            q019_runtime_box_edge_last_power_fraction &&
        std::isfinite(q019_runtime_box_edge_last_fluctuation_mean) &&
        q019_runtime_box_edge_last_fluctuation_mean > 0.0;
    const bool zero_fluctuation_state =
        q019_runtime_box_edge_last_status ==
            static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation) &&
        q019_runtime_box_edge_last_power_fraction == -1.0 &&
        q019_runtime_box_edge_last_fluctuation_mean == 0.0;
    const bool other_unavailable_state =
        q019_runtime_box_edge_last_status !=
            static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled) &&
        q019_runtime_box_edge_last_status !=
            static_cast<int>(BoxEdgeDiagnosticStatus::valid) &&
        q019_runtime_box_edge_last_status !=
            static_cast<int>(BoxEdgeDiagnosticStatus::zero_fluctuation) &&
        q019_runtime_box_edge_last_power_fraction == -1.0 &&
        q019_runtime_box_edge_last_fluctuation_mean == -1.0;
    const bool observed_state =
        q019_runtime_box_edge_completed_slots > 0 &&
        q019_runtime_box_edge_last_cycle > 0 &&
        q019_runtime_box_edge_last_cycle <= pmy_mesh_->ncycle &&
        (q019_runtime_box_edge_status_mask &
         BoxEdgeStatusBit(q019_runtime_box_edge_last_status)) != 0 &&
        (q019_runtime_box_edge_last_status !=
             static_cast<int>(BoxEdgeDiagnosticStatus::cadence_skipped) ||
         q019_runtime_box_edge_skipped_slots > 0) &&
        std::isfinite(q019_runtime_box_edge_last_prior_time) &&
        std::isfinite(q019_runtime_box_edge_last_time) &&
        BoxEdgeFirstCrossingChronologyIsValid(
            q019_runtime_box_edge_last_prior_time,
            q019_runtime_box_edge_last_nominal_time,
            q019_runtime_box_edge_last_time) &&
        q019_runtime_box_edge_last_time <=
            pmy_mesh_->time + chronology_tolerance &&
        q019_runtime_box_edge_next_nominal_time >
            pmy_mesh_->time &&
        maximum_state &&
        (valid_metric_state || zero_fluctuation_state ||
         other_unavailable_state);
    if (!(schedule_state && (unobserved_state || observed_state))) {
      Q019NonlinearFatal("Q019 box-edge monitor restart chronology drifted");
    }
  } else {
    if (!q019_runtime_controller_enabled &&
        (pin->GetInteger(block, "runtime_box_edge_monitor_schema") != 2 ||
        pin->GetInteger(block, "runtime_box_edge_monitor_last_cycle") != 0 ||
        pin->GetInteger(block, "runtime_box_edge_monitor_completed_slots") != 0 ||
        pin->GetInteger(block, "runtime_box_edge_monitor_valid_samples") != 0 ||
        pin->GetInteger(block, "runtime_box_edge_monitor_skipped_slots") != 0 ||
        pin->GetInteger(block, "runtime_box_edge_monitor_last_status") !=
            static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled) ||
        pin->GetInteger(block, "runtime_box_edge_monitor_status_mask") != 0)) {
      Q019NonlinearFatal("Q019 box-edge monitor initial integer state drifted");
    }
    if (!q019_runtime_controller_enabled) {
      Q019RequireClose(
          "initial box-edge next nominal time",
          pin->GetReal(block, "runtime_box_edge_monitor_next_nominal_time"),
          q019_runtime_box_edge_monitor_dt);
      Q019RequireClose(
          "initial box-edge last nominal time",
          pin->GetReal(block, "runtime_box_edge_monitor_last_nominal_time"), -1.0);
      Q019RequireClose(
          "initial box-edge last prior time",
          pin->GetReal(block, "runtime_box_edge_monitor_last_prior_time"), 0.0);
      Q019RequireClose("initial box-edge last time",
                       pin->GetReal(block, "runtime_box_edge_monitor_last_time"), 0.0);
      Q019RequireClose(
          "initial box-edge last power fraction",
          pin->GetReal(block, "runtime_box_edge_monitor_last_power_fraction"), -1.0);
      Q019RequireClose(
          "initial box-edge maximum power fraction",
          pin->GetReal(block, "runtime_box_edge_monitor_max_power_fraction"), -1.0);
      Q019RequireClose(
          "initial box-edge fluctuation mean",
          pin->GetReal(block, "runtime_box_edge_monitor_last_fluctuation_mean"), -1.0);
    }
    q019_runtime_box_edge_next_nominal_time =
        q019_runtime_box_edge_monitor_dt;
    q019_runtime_box_edge_last_nominal_time = -1.0;
    q019_runtime_box_edge_last_cycle = 0;
    q019_runtime_box_edge_completed_slots = 0;
    q019_runtime_box_edge_valid_samples = 0;
    q019_runtime_box_edge_skipped_slots = 0;
    q019_runtime_box_edge_last_status =
        static_cast<int>(BoxEdgeDiagnosticStatus::not_sampled);
    q019_runtime_box_edge_status_mask = 0;
    q019_runtime_box_edge_last_prior_time = 0.0;
    q019_runtime_box_edge_last_time = 0.0;
    q019_runtime_box_edge_last_power_fraction = -1.0;
    q019_runtime_box_edge_max_power_fraction = -1.0;
    q019_runtime_box_edge_last_fluctuation_mean = -1.0;
  }
  Q019StoreBoxEdgeMonitorState();
  Q019StoreRuntimeControllerState();
  user_work_in_loop = true;
  user_work_in_loop_func = Q019RuntimeDiagnostics;
  user_hist = true;
  user_hist_func = Q019BoxEdgeHistory;
  pgen_final_func = Q019FinalEvidenceStatus;

  if (restart && q019_runtime_controller_triggered) {
    RequestUserStop(q019_runtime_controller_trigger_reason,
                    q019_runtime_controller_trigger_failure);
  }

  if (restart) return;

  const SeedParameters parameters = {
    dimension, field_seed, b_g, k0, epsilon, rho, eigenmode_amplitude,
    broadband_amplitude, separated_long_mode_amplitude, x1_extent, x2_extent,
    x3_extent, seed_topology
  };
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;

  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = indcs.nx2 + 2*indcs.ng;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 2;
  DvceArray4D<Real> a1;
  DvceArray4D<Real> a2;
  DvceArray4D<Real> a3;
  Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

  par_for("pgen_q019_corrected_nonlinear_bell_vector_potential", DevExeSpace(),
          0, nmb - 1, ks, ke + 1, js, je + 1, is, ie + 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1v = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i - is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x2v = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j - js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real x3v = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k - ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Vector3 av1 = AxisAlignedVectorPotentialAt(parameters, {x1v, x2f, x3f});
    const Vector3 av2 = AxisAlignedVectorPotentialAt(parameters, {x1f, x2v, x3f});
    const Vector3 av3 = AxisAlignedVectorPotentialAt(parameters, {x1f, x2f, x3v});
    a1(m, k, j, i) = av1.x1;
    a2(m, k, j, i) = av2.x2;
    a3(m, k, j, i) = av3.x3;
  });

  par_for("pgen_q019_corrected_nonlinear_bell_ct_field", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m, k, j, i) = (a3(m, k, j + 1, i) - a3(m, k, j, i))/dx2 -
                          (a2(m, k + 1, j, i) - a2(m, k, j, i))/dx3;
    b0.x2f(m, k, j, i) = (a1(m, k + 1, j, i) - a1(m, k, j, i))/dx3 -
                          (a3(m, k, j, i + 1) - a3(m, k, j, i))/dx1;
    b0.x3f(m, k, j, i) = (a2(m, k, j, i + 1) - a2(m, k, j, i))/dx1 -
                          (a1(m, k, j + 1, i) - a1(m, k, j, i))/dx2;
    if (i == ie) {
      b0.x1f(m, k, j, i + 1) =
          (a3(m, k, j + 1, i + 1) - a3(m, k, j, i + 1))/dx2 -
          (a2(m, k + 1, j, i + 1) - a2(m, k, j, i + 1))/dx3;
    }
    if (j == je) {
      b0.x2f(m, k, j + 1, i) =
          (a1(m, k + 1, j + 1, i) - a1(m, k, j + 1, i))/dx3 -
          (a3(m, k, j + 1, i + 1) - a3(m, k, j + 1, i))/dx1;
    }
    if (k == ke) {
      b0.x3f(m, k + 1, j, i) =
          (a2(m, k + 1, j, i + 1) - a2(m, k + 1, j, i))/dx1 -
          (a1(m, k + 1, j + 1, i) - a1(m, k + 1, j, i))/dx2;
    }
  });

  par_for("pgen_q019_corrected_nonlinear_bell_primitives", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1 = CellCenterX(i - is, indcs.nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x2 = CellCenterX(j - js, indcs.nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x3 = CellCenterX(k - ks, indcs.nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const Vector3 velocity = AxisAlignedEigenmodeVelocityAt(parameters, {x1, x2, x3});
    w0(m, IDN, k, j, i) = rho;
    w0(m, IVX, k, j, i) = velocity.x1;
    w0(m, IVY, k, j, i) = velocity.x2;
    w0(m, IVZ, k, j, i) = velocity.x3;
    w0(m, IEN, k, j, i) = pressure;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);

  auto *ppart = pmbp->ppart;
  if (!high_rigidity) {
    const int finite_ppc = integral_ppc;
    const int nx1_local = indcs.nx1;
    const int nx2_local = indcs.nx2;
    const int nx3_local = indcs.nx3;
    const int cells_per_meshblock = nx1_local*nx2_local*nx3_local;
    const int particles_per_meshblock = finite_ppc*cells_per_meshblock;
    const int root_nx1 = pmy_mesh_->mesh_indcs.nx1;
    const int root_nx2 = pmy_mesh_->mesh_indcs.nx2;
    const int root_nx3 = pmy_mesh_->mesh_indcs.nx3;
    if (!OctahedralShellPacketLayoutIsValid(
            finite_ppc, cells_per_meshblock, nmb, ppart->nprtcl_thispack)) {
      Q019NonlinearFatal("finite-rigidity particle count must be exactly PPC times "
                         "MeshBlock cells with PPC divisible by six for every block");
    }
    HostArray1D<std::uint64_t> h_global_cell_base(
        "q019_global_root_cell_base_host", nmb);
    for (int m = 0; m < nmb; ++m) {
      const LogicalLocation &location = pmy_mesh_->lloc_eachmb[pmbp->gids + m];
      if (location.level != pmy_mesh_->root_level) {
        Q019NonlinearFatal("Q019 nested packet identity requires root-level MeshBlocks");
      }
      const int global_i = location.lx1*nx1_local;
      const int global_j = location.lx2*nx2_local;
      const int global_k = location.lx3*nx3_local;
      if (global_i < 0 || global_i + nx1_local > root_nx1 ||
          global_j < 0 || global_j + nx2_local > root_nx2 ||
          global_k < 0 || global_k + nx3_local > root_nx3) {
        Q019NonlinearFatal("Q019 global root-cell identity is outside the root mesh");
      }
      h_global_cell_base(m) = GlobalCellLinearId(
          global_i, global_j, global_k, root_nx1, root_nx2);
    }
    auto global_cell_base =
        Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_global_cell_base);
    auto &pi = ppart->prtcl_idata;
    auto &pr = ppart->prtcl_rdata;
    const bool momentum_state = ppart->UsesMomentumState();
    const Real cr_light_speed = ppart->pic_cr_light_speed;
    const int gids = pmbp->gids;
    int initializer_order_mismatches = 0;
    Kokkos::parallel_reduce(
        "pgen_q019_verify_initializer_block_order",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, ppart->nprtcl_thispack),
    KOKKOS_LAMBDA(const int p, int &mismatches) {
      const int expected_gid = gids + p/particles_per_meshblock;
      if (pi(PGID, p) != expected_gid) ++mismatches;
    }, initializer_order_mismatches);
    if (initializer_order_mismatches != 0) {
      Q019NonlinearFatal("finite-rigidity initializer block ordering does not support "
                         "the declared p modulo PPC isotropic-shell packet grouping");
    }
    if (quiet_isotropic_shell_packet) {
      par_for("pgen_q019_finite_cell_centered_isotropic_shell_packet", DevExeSpace(),
              0, ppart->nprtcl_thispack - 1,
      KOKKOS_LAMBDA(const int p) {
        const int m = p/particles_per_meshblock;
        const int cell_linear = (p % particles_per_meshblock)/finite_ppc;
        const int ci = cell_linear % nx1_local;
        const int remainder = cell_linear/nx1_local;
        const int cj = remainder % nx2_local;
        const int ck = remainder/nx2_local;
        pr(IPX, p) = size.d_view(m).x1min +
            (static_cast<Real>(ci) + 0.5)*size.d_view(m).dx1;
        pr(IPY, p) = size.d_view(m).x2min +
            (static_cast<Real>(cj) + 0.5)*size.d_view(m).dx2;
        pr(IPZ, p) = nx3_local > 1 ?
            size.d_view(m).x3min + (static_cast<Real>(ck) + 0.5)*size.d_view(m).dx3 :
            0.0;
      });
      Kokkos::fence();
    }
    par_for("pgen_q019_finite_nested_haar_octahedral_shell", DevExeSpace(),
            0, ppart->nprtcl_thispack - 1,
    KOKKOS_LAMBDA(const int p) {
      const int m = p/particles_per_meshblock;
      const int particle_in_meshblock = p % particles_per_meshblock;
      const int cell_linear = particle_in_meshblock/finite_ppc;
      const int sample_within_cell = particle_in_meshblock % finite_ppc;
      const int ci = cell_linear % nx1_local;
      const int remainder = cell_linear/nx1_local;
      const int cj = remainder % nx2_local;
      const int ck = remainder/nx2_local;
      const std::uint64_t global_cell_id = global_cell_base(m) +
          static_cast<std::uint64_t>(ck)*static_cast<std::uint64_t>(root_nx2)*
              static_cast<std::uint64_t>(root_nx1) +
          static_cast<std::uint64_t>(cj)*static_cast<std::uint64_t>(root_nx1) +
          static_cast<std::uint64_t>(ci);
      const Vector3 velocity = NestedHaarOctahedralShellVelocityAtSample(
          finite_ppc, sample_within_cell, global_cell_id, field_seed, particle_seed,
          guide_parallel_stream_speed, finite_shell_speed);
      const Real v2 = velocity.x1*velocity.x1 + velocity.x2*velocity.x2 +
                      velocity.x3*velocity.x3;
      const Real gamma = momentum_state ?
          1.0/std::sqrt(1.0 - v2/(cr_light_speed*cr_light_speed)) : 1.0;
      pr(IPVX, p) = gamma*velocity.x1;
      pr(IPVY, p) = gamma*velocity.x2;
      pr(IPVZ, p) = gamma*velocity.x3;
    });
    Kokkos::fence();
  }

  auto &particle_idata = ppart->prtcl_idata;
  auto &particle_rdata = ppart->prtcl_rdata;
  const int gids = pmbp->gids;
  par_for("pgen_q019_initialize_particle_magnetic_field_cache", DevExeSpace(),
          0, ppart->nprtcl_thispack - 1,
  KOKKOS_LAMBDA(const int p) {
    const int m = particle_idata(PGID, p) - gids;
    if (m < 0 || m >= nmb) return;
    const Vector3 field = AxisAlignedMagneticFieldAt(
        parameters,
        {particle_rdata(IPX, p), particle_rdata(IPY, p),
         dimension == 3 ? particle_rdata(IPZ, p) : 0.0});
    particle_rdata(IPBX, p) = field.x1;
    particle_rdata(IPBY, p) = field.x2;
    particle_rdata(IPBZ, p) = field.x3;
  });
  Kokkos::fence();
}
