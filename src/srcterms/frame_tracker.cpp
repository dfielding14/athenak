//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file frame_tracker.cpp
//! \brief Implements shared frame tracking for hydro/MHD problem generators.

#include "frame_tracker.hpp"

#include <math.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

enum FrameTrackingMode {
  kFTVelocity = 0,
  kFTPosition = 1,
  kFTPD = 2,
  kFTTopMdot = 3,
  kFTTopMdotAuto = 4,
  kFTTopMdotWindow = 5
};

enum FrameTrackingPositionSignal {
  kFTPosCentroid = 0,
  kFTPosBandMidpoint = 1,
  kFTPosBlend = 2
};

enum FrameTrackingVelocitySignal {
  kFTVelMaterialMean = 0,
  kFTVelPositionRate = 1
};

enum FrameTrackingBoostChangeMode {
  kFTBoostPerApply = 0,
  kFTBoostPerTime = 1
};

enum FrameTrackingTargetKind {
  kFTTargetDensity = 0,
  kFTTargetTemperature = 1,
  kFTTargetPressure = 2,
  kFTTargetEntropy = 3,
  kFTTargetScalar = 4,
  kFTTargetVelocity1 = 5,
  kFTTargetVelocity2 = 6,
  kFTTargetVelocity3 = 7,
  kFTTargetSpeed = 8,
  kFTTargetInternalEnergy = 9
};

enum FrameTrackingWeightMode {
  kFTWeightMass = 0,
  kFTWeightVolume = 1,
  kFTWeightTarget = 2,
  kFTWeightTracerMass = 3
};

constexpr int kFrameTrackingStateVersion = 3;

const std::vector<std::pair<std::string, std::string>> &RemovedAliases() {
  static const std::vector<std::pair<std::string, std::string>> aliases = {
      {"enable", "enabled"},
      {"ft_apply_every", "apply_every"},
      {"t_start", "start_time"},
      {"t_frame_tracking_start", "start_time"},
      {"ndiag", "diagnostic_every"},
      {"target_variable", "target"},
      {"derived_variable", "target"},
      {"min", "target_min"},
      {"max", "target_max"},
      {"ft_temp_lo", "target_min"},
      {"ft_temp_hi", "target_max"},
      {"center", "target_center"},
      {"reacquire_scale", "target_scale"},
      {"ft_mode", "mode"},
      {"ft_position_signal", "position_signal"},
      {"ft_max_boost_change_mode", "max_boost_change_mode"},
      {"ft_position_blend", "position_blend"},
      {"ft_tau_avg", "tau_avg"},
      {"ft_tau_relax", "tau_relax"},
      {"ft_tau_vel", "tau_vel"},
      {"ft_tau_int", "tau_int"},
      {"ft_int_max_abs", "int_max_abs"},
      {"ft_int_leak_tau", "int_leak_tau"},
      {"ft_int_unsat_only", "int_unsat_only"},
      {"ft_weight_floor", "weight_floor"},
      {"ft_min_global_weight", "min_global_weight"},
      {"ft_max_abs_boost", "max_abs_boost"},
      {"ft_max_boost_change", "max_boost_change"},
      {"ft_max_boost_change_rate", "max_boost_change_rate"},
      {"ft_reacquire_expand_factor", "reacquire_expand_factor"},
      {"ft_reacquire_max_expand", "reacquire_max_expand"},
      {"ft_reacquire_recover_updates", "reacquire_recover_updates"},
      {"ft_reacquire_leak_on_miss", "reacquire_leak_on_miss"},
      {"ft_boundary_guard_enable", "boundary_guard_enable"},
      {"ft_boundary_guard_cells", "boundary_guard_cells"},
      {"ft_boundary_guard_min_scale", "boundary_guard_min_scale"},
      {"axis", "axes"},
      {"directions", "axes"},
      {"z_target", "x3_target"},
      {"ft_z_target", "x3_target"},
      {"z_lower_limit", "x3_lower_limit"},
      {"ft_lowzlim", "x3_lower_limit"},
      {"z_upper_limit", "x3_upper_limit"},
      {"ft_uppzlim", "x3_upper_limit"},
      {"ft_use_lowzlim", "x3_lower_limit"},
      {"ft_use_uppzlim", "x3_upper_limit"},
      {"x1_target_position", "x1_target"},
      {"x2_target_position", "x2_target"},
      {"x3_target_position", "x3_target"},
      {"target_x1", "x1_target"},
      {"target_x2", "x2_target"},
      {"target_x3", "x3_target"},
      {"x1_min_limit", "x1_lower_limit"},
      {"x2_min_limit", "x2_lower_limit"},
      {"x3_min_limit", "x3_lower_limit"},
      {"x1_low_limit", "x1_lower_limit"},
      {"x2_low_limit", "x2_lower_limit"},
      {"x3_low_limit", "x3_lower_limit"},
      {"x1_max_limit", "x1_upper_limit"},
      {"x2_max_limit", "x2_upper_limit"},
      {"x3_max_limit", "x3_upper_limit"},
      {"x1_upp_limit", "x1_upper_limit"},
      {"x2_upp_limit", "x2_upper_limit"},
      {"x3_upp_limit", "x3_upper_limit"}
  };
  return aliases;
}

std::string TrimCopy(const std::string &input) {
  const auto first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const auto last = input.find_last_not_of(" \t\r\n");
  std::string trimmed = input.substr(first, last - first + 1);
  if (trimmed.size() >= 2) {
    const char front = trimmed.front();
    const char back = trimmed.back();
    if ((front == '"' && back == '"') || (front == '\'' && back == '\'')) {
      trimmed = trimmed.substr(1, trimmed.size() - 2);
      const auto inner_first = trimmed.find_first_not_of(" \t\r\n");
      if (inner_first == std::string::npos) {
        return "";
      }
      const auto inner_last = trimmed.find_last_not_of(" \t\r\n");
      trimmed = trimmed.substr(inner_first, inner_last - inner_first + 1);
    }
  }
  return trimmed;
}

std::string ToLowerCopy(std::string input) {
  std::transform(input.begin(), input.end(), input.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return input;
}

std::string NormalizeToken(const std::string &input) {
  std::string out = ToLowerCopy(TrimCopy(input));
  for (char &c : out) {
    if (c == '-' || c == ' ') {
      c = '_';
    }
  }
  return out;
}

void FatalFrameTrackingInput(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

bool ReadRealAny(ParameterInput *pin, const std::string &block,
                 const std::vector<std::string> &names, Real &value) {
  for (const auto &name : names) {
    if (pin->DoesParameterExist(block, name)) {
      value = pin->GetReal(block, name);
      return true;
    }
  }
  return false;
}

bool ReadIntegerAny(ParameterInput *pin, const std::string &block,
                    const std::vector<std::string> &names, int &value) {
  for (const auto &name : names) {
    if (pin->DoesParameterExist(block, name)) {
      value = pin->GetInteger(block, name);
      return true;
    }
  }
  return false;
}

bool ReadBooleanAny(ParameterInput *pin, const std::string &block,
                    const std::vector<std::string> &names, bool &value) {
  for (const auto &name : names) {
    if (pin->DoesParameterExist(block, name)) {
      value = pin->GetBoolean(block, name);
      return true;
    }
  }
  return false;
}

bool ReadStringAny(ParameterInput *pin, const std::string &block,
                   const std::vector<std::string> &names, std::string &value) {
  for (const auto &name : names) {
    if (pin->DoesParameterExist(block, name)) {
      value = pin->GetString(block, name);
      return true;
    }
  }
  return false;
}

int ParseFrameTrackingMode(const std::string &mode_raw) {
  const std::string mode = NormalizeToken(mode_raw);
  if (mode == "velocity") return kFTVelocity;
  if (mode == "position") return kFTPosition;
  if (mode == "pd" || mode == "position_velocity") return kFTPD;
  if (mode == "top_mdot" || mode == "mdot_top" ||
      mode == "delayed_top_mdot" || mode == "vtop") {
    return kFTTopMdot;
  }
  if (mode == "top_mdot_auto" || mode == "auto_top_mdot" ||
      mode == "settled_top_mdot") {
    return kFTTopMdotAuto;
  }
  if (mode == "top_mdot_window" || mode == "window_top_mdot" ||
      mode == "running_top_mdot" || mode == "finite_window_top_mdot") {
    return kFTTopMdotWindow;
  }
  FatalFrameTrackingInput("Invalid <frame_tracking>/mode '" + mode_raw +
                          "'. Expected one of: velocity, position, pd, "
                          "top_mdot, top_mdot_auto, top_mdot_window.");
  return kFTPD;
}

Real ParseMdotVelocitySign(const std::string &sign_raw) {
  const std::string sign = NormalizeToken(sign_raw);
  if (sign == "normal" || sign == "same" || sign == "positive" ||
      sign == "plus" || sign == "+1") {
    return 1.0;
  }
  if (sign == "flipped" || sign == "flip" || sign == "opposite" ||
      sign == "negative" || sign == "minus" || sign == "-1") {
    return -1.0;
  }
  FatalFrameTrackingInput("Invalid <frame_tracking>/mdot_velocity_sign '" +
                          sign_raw + "'. Expected normal or flipped.");
  return 1.0;
}

Real SmoothStep01(Real x) {
  x = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), x));
  return x*x*(static_cast<Real>(3.0) - static_cast<Real>(2.0)*x);
}

int ParseFrameTrackingPositionSignal(const std::string &signal_raw) {
  const std::string signal = NormalizeToken(signal_raw);
  if (signal == "centroid") return kFTPosCentroid;
  if (signal == "band_midpoint" || signal == "range_midpoint") return kFTPosBandMidpoint;
  if (signal == "blend") return kFTPosBlend;
  FatalFrameTrackingInput("Invalid <frame_tracking>/position_signal '" + signal_raw +
                          "'. Expected one of: centroid, band_midpoint, blend.");
  return kFTPosBlend;
}

int ParseFrameTrackingVelocitySignal(const std::string &signal_raw) {
  const std::string signal = NormalizeToken(signal_raw);
  if (signal == "material_mean" || signal == "mean_velocity") {
    return kFTVelMaterialMean;
  }
  if (signal == "position_rate" || signal == "position_derivative") {
    return kFTVelPositionRate;
  }
  FatalFrameTrackingInput("Invalid <frame_tracking>/velocity_signal '" + signal_raw +
                          "'. Expected one of: material_mean, position_rate.");
  return kFTVelMaterialMean;
}

int ParseFrameTrackingBoostChangeMode(const std::string &mode_raw) {
  const std::string mode = NormalizeToken(mode_raw);
  if (mode == "per_apply" || mode == "per_update") return kFTBoostPerApply;
  if (mode == "per_time") return kFTBoostPerTime;
  FatalFrameTrackingInput("Invalid <frame_tracking>/max_boost_change_mode '" + mode_raw +
                          "'. Expected one of: per_apply, per_time.");
  return kFTBoostPerApply;
}

int ParseWeightMode(const std::string &mode_raw) {
  const std::string mode = NormalizeToken(mode_raw);
  if (mode == "mass" || mode == "density") return kFTWeightMass;
  if (mode == "volume") return kFTWeightVolume;
  if (mode == "target" || mode == "value") return kFTWeightTarget;
  if (mode == "tracer_mass") return kFTWeightTracerMass;
  FatalFrameTrackingInput("Invalid <frame_tracking>/weight '" + mode_raw +
                          "'. Expected one of: mass, volume, target, tracer_mass.");
  return kFTWeightMass;
}

bool ParseScalarTargetSuffix(const std::string &target, int &scalar_index) {
  if (target.compare(0, 6, "scalar") != 0) {
    return false;
  }
  std::string suffix = target.substr(6);
  if (!suffix.empty() && suffix.front() == '_') {
    suffix.erase(suffix.begin());
  }
  if (suffix.empty()) {
    return true;
  }
  for (char c : suffix) {
    if (!std::isdigit(static_cast<unsigned char>(c))) {
      return false;
    }
  }
  scalar_index = std::atoi(suffix.c_str());
  return true;
}

int ParseTargetKind(const std::string &target_raw, int &scalar_index) {
  const std::string target = NormalizeToken(target_raw);
  if (target == "density" || target == "dens" || target == "rho") return kFTTargetDensity;
  if (target == "temperature" || target == "temp" || target == "t") {
    return kFTTargetTemperature;
  }
  if (target == "pressure" || target == "pres" || target == "p") return kFTTargetPressure;
  if (target == "entropy" || target == "adiabat") return kFTTargetEntropy;
  if (target == "internal_energy" || target == "eint") return kFTTargetInternalEnergy;
  if (target == "v1" || target == "vx" || target == "velx" ||
      target == "velocity1" || target == "velocity_x1" || target == "x1_velocity") {
    return kFTTargetVelocity1;
  }
  if (target == "v2" || target == "vy" || target == "vely" ||
      target == "velocity2" || target == "velocity_x2" || target == "x2_velocity") {
    return kFTTargetVelocity2;
  }
  if (target == "v3" || target == "vz" || target == "velz" ||
      target == "velocity3" || target == "velocity_x3" || target == "x3_velocity") {
    return kFTTargetVelocity3;
  }
  if (target == "speed" || target == "velocity" || target == "velocity_magnitude") {
    return kFTTargetSpeed;
  }
  if (ParseScalarTargetSuffix(target, scalar_index)) {
    return kFTTargetScalar;
  }
  FatalFrameTrackingInput("Invalid <frame_tracking>/target '" + target_raw +
                          "'. Expected density, temperature, pressure, entropy, "
                          "scalar/scalarN, v1/v2/v3, speed, or internal_energy.");
  return kFTTargetDensity;
}

const char *AxisName(const int axis) {
  if (axis == 0) return "x1";
  if (axis == 1) return "x2";
  return "x3";
}

KOKKOS_INLINE_FUNCTION
int AxisCellCount(const int axis, const RegionIndcs &indcs) {
  if (axis == 0) return indcs.nx1;
  if (axis == 1) return indcs.nx2;
  return indcs.nx3;
}

Real AxisGlobalMin(const int axis, const Mesh *pm) {
  if (axis == 0) return pm->mesh_size.x1min;
  if (axis == 1) return pm->mesh_size.x2min;
  return pm->mesh_size.x3min;
}

Real AxisGlobalMax(const int axis, const Mesh *pm) {
  if (axis == 0) return pm->mesh_size.x1max;
  if (axis == 1) return pm->mesh_size.x2max;
  return pm->mesh_size.x3max;
}

KOKKOS_INLINE_FUNCTION
Real FrameTrackingPressure(const DvceArray5D<Real> &w0, const EOS_Data &eos,
                           const int m, const int k, const int j, const int i) {
  const Real dens = w0(m, IDN, k, j, i);
  if (eos.use_e) {
    return eos.IdealGasPressure(w0(m, IEN, k, j, i));
  }
  if (eos.use_t) {
    return dens*w0(m, ITM, k, j, i);
  }
  return dens*SQR(eos.iso_cs);
}

KOKKOS_INLINE_FUNCTION
Real FrameTrackingTemperature(const DvceArray5D<Real> &w0, const EOS_Data &eos,
                              const int m, const int k, const int j, const int i) {
  const Real dens = w0(m, IDN, k, j, i);
  if (dens <= 0.0) {
    return 0.0;
  }
  if (eos.use_e) {
    return eos.IdealGasPressure(w0(m, IEN, k, j, i))/dens;
  }
  if (eos.use_t) {
    return w0(m, ITM, k, j, i);
  }
  return SQR(eos.iso_cs);
}

KOKKOS_INLINE_FUNCTION
Real FrameTrackingTargetValue(const int target_kind, const int scalar_var,
                              const DvceArray5D<Real> &w0, const EOS_Data &eos,
                              const int m, const int k, const int j, const int i) {
  const Real dens = w0(m, IDN, k, j, i);
  if (target_kind == kFTTargetDensity) {
    return dens;
  }
  if (target_kind == kFTTargetTemperature) {
    return FrameTrackingTemperature(w0, eos, m, k, j, i);
  }
  if (target_kind == kFTTargetPressure) {
    return FrameTrackingPressure(w0, eos, m, k, j, i);
  }
  if (target_kind == kFTTargetEntropy) {
    const Real pressure = FrameTrackingPressure(w0, eos, m, k, j, i);
    return (dens > 0.0 && eos.gamma > 0.0) ? pressure/pow(dens, eos.gamma) : 0.0;
  }
  if (target_kind == kFTTargetScalar) {
    return w0(m, scalar_var, k, j, i);
  }
  if (target_kind == kFTTargetVelocity1) {
    return w0(m, IVX, k, j, i);
  }
  if (target_kind == kFTTargetVelocity2) {
    return w0(m, IVY, k, j, i);
  }
  if (target_kind == kFTTargetVelocity3) {
    return w0(m, IVZ, k, j, i);
  }
  if (target_kind == kFTTargetSpeed) {
    return sqrt(SQR(w0(m, IVX, k, j, i)) +
                SQR(w0(m, IVY, k, j, i)) +
                SQR(w0(m, IVZ, k, j, i)));
  }
  if (target_kind == kFTTargetInternalEnergy) {
    if (eos.use_e) {
      return w0(m, IEN, k, j, i);
    }
    if (eos.use_t && eos.gamma > 1.0) {
      return dens*w0(m, ITM, k, j, i)/(eos.gamma - 1.0);
    }
    return 0.0;
  }
  return dens;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn FrameTracker::FrameTracker()
//! \brief Parse <frame_tracking> and initialize shared tracking state.

FrameTracker::FrameTracker(MeshBlockPack *pp, ParameterInput *pin,
                           const std::string &block_name) :
    pmy_pack(pp),
    block_name_(block_name) {
  RejectRemovedAliases(pin);
  enabled_ = pin->GetOrAddBoolean(block_name_, "enabled", true);
  apply_every_ = pin->GetOrAddInteger(block_name_, "apply_every", 1);
  start_time_ = pin->GetOrAddReal(block_name_, "start_time", 0.0);
  diagnostic_every_ = pin->GetOrAddInteger(block_name_, "diagnostic_every", -1);
  verbose_ = pin->GetOrAddBoolean(block_name_, "verbose", false);

  target_name_ = pin->GetOrAddString(block_name_, "target", "density");
  scalar_index_ = pin->GetOrAddInteger(block_name_, "scalar_index", 0);
  target_kind_ = ParseTargetKind(target_name_, scalar_index_);
  weight_mode_ = ParseWeightMode(pin->GetOrAddString(block_name_, "weight", "mass"));

  const Real huge = std::numeric_limits<Real>::max()/static_cast<Real>(4.0);
  target_min_ = -huge;
  target_max_ = huge;
  (void) ReadRealAny(pin, block_name_, {"target_min"}, target_min_);
  (void) ReadRealAny(pin, block_name_, {"target_max"}, target_max_);

  if (!ReadRealAny(pin, block_name_, {"target_center"}, target_center_)) {
    if (target_min_ > 0.0 && target_max_ > 0.0 && target_max_ < huge) {
      target_center_ = std::sqrt(target_min_*target_max_);
    } else if (target_min_ > -huge && target_max_ < huge) {
      target_center_ = 0.5*(target_min_ + target_max_);
    } else {
      target_center_ = 0.0;
    }
  }
  target_min_active_ = target_min_;
  target_max_active_ = target_max_;

  std::string target_scale = "auto";
  (void) ReadStringAny(pin, block_name_, {"target_scale"}, target_scale);
  target_scale = NormalizeToken(target_scale);
  if (target_scale == "log" || target_scale == "logarithmic") {
    use_log_reacquire_ = true;
  } else if (target_scale == "linear") {
    use_log_reacquire_ = false;
  } else if (target_scale == "auto") {
    use_log_reacquire_ =
        (target_min_ > 0.0 && target_max_ > 0.0 && target_center_ > 0.0);
  } else {
    FatalFrameTrackingInput("Invalid <frame_tracking>/target_scale '" + target_scale +
                            "'. Expected auto, linear, or log.");
  }

  mode_name_ = pin->GetOrAddString(block_name_, "mode", "pd");
  mode_ = ParseFrameTrackingMode(mode_name_);
  position_signal_name_ =
      pin->GetOrAddString(block_name_, "position_signal", "blend");
  position_signal_ = ParseFrameTrackingPositionSignal(position_signal_name_);
  velocity_signal_name_ =
      pin->GetOrAddString(block_name_, "velocity_signal", "material_mean");
  velocity_signal_ = ParseFrameTrackingVelocitySignal(velocity_signal_name_);
  boost_change_mode_name_ =
      pin->GetOrAddString(block_name_, "max_boost_change_mode", "per_apply");
  boost_change_mode_ = ParseFrameTrackingBoostChangeMode(
      boost_change_mode_name_);
  position_blend_ = pin->GetOrAddReal(block_name_, "position_blend", 0.7);
  tau_avg_ = pin->GetOrAddReal(block_name_, "tau_avg", 1.0);
  tau_relax_ = pin->GetOrAddReal(block_name_, "tau_relax", 1.0);
  tau_vel_ = pin->GetOrAddReal(block_name_, "tau_vel", 1.0);
  tau_int_ = pin->GetOrAddReal(block_name_, "tau_int", 0.0);
  int_max_abs_ = pin->GetOrAddReal(block_name_, "int_max_abs", 0.0);
  int_leak_tau_ = pin->GetOrAddReal(block_name_, "int_leak_tau", 1.0);
  int_unsat_only_ = pin->GetOrAddBoolean(block_name_, "int_unsat_only", true);
  weight_floor_ = pin->GetOrAddReal(block_name_, "weight_floor", 0.0);
  min_global_weight_ = pin->GetOrAddReal(block_name_, "min_global_weight", 0.0);
  max_abs_boost_ = pin->GetOrAddReal(block_name_, "max_abs_boost", 0.0);
  max_boost_change_ = pin->GetOrAddReal(block_name_, "max_boost_change", 0.0);
  max_boost_change_rate_ = pin->GetOrAddReal(block_name_, "max_boost_change_rate", -1.0);
  mdot_average_start_time_ =
      pin->GetOrAddReal(block_name_, "mdot_average_start_time", 0.0);
  mdot_average_end_time_ =
      pin->GetOrAddReal(block_name_, "mdot_average_end_time",
                        (start_time_ > 0.0) ? start_time_ : -1.0);
  mdot_ramp_time_ = pin->GetOrAddReal(block_name_, "mdot_ramp_time", 0.0);
  mdot_window_start_time_ =
      pin->DoesParameterExist(block_name_, "mdot_start_time") ?
      pin->GetReal(block_name_, "mdot_start_time") :
      pin->GetOrAddReal(block_name_, "mdot_window_start_time", 0.5);
  mdot_window_time_ =
      pin->GetOrAddReal(block_name_, "mdot_window_time", 1.0);
  mdot_density_mode_name_ =
      pin->GetOrAddString(block_name_, "mdot_density_mode", "fixed");
  mdot_velocity_sign_name_ =
      pin->GetOrAddString(block_name_, "mdot_velocity_sign", "normal");
  mdot_velocity_sign_ = ParseMdotVelocitySign(mdot_velocity_sign_name_);
  mdot_velocity_factor_ =
      pin->GetOrAddReal(block_name_, "mdot_velocity_factor", 1.0);
  mdot_flux_frame_name_ =
      pin->GetOrAddString(block_name_, "mdot_flux_frame", "grid");
  const std::string mdot_flux_frame = NormalizeToken(mdot_flux_frame_name_);
  if (mdot_flux_frame == "lab" || mdot_flux_frame == "laboratory" ||
      mdot_flux_frame == "inertial" || mdot_flux_frame == "physical") {
    mdot_flux_lab_frame_ = true;
  } else if (mdot_flux_frame == "grid" || mdot_flux_frame == "frame" ||
             mdot_flux_frame == "moving_frame" ||
             mdot_flux_frame == "comoving") {
    mdot_flux_lab_frame_ = false;
  } else {
    FatalFrameTrackingInput("Invalid <frame_tracking>/mdot_flux_frame '" +
                            mdot_flux_frame_name_ + "'. Expected grid or lab.");
  }
  if (pin->DoesParameterExist(block_name_, "mdot_density")) {
    mdot_density_ = pin->GetReal(block_name_, "mdot_density");
  } else if (pin->DoesBlockExist("problem") &&
             pin->DoesParameterExist("problem", "rho_cold")) {
    mdot_density_ = pin->GetReal("problem", "rho_cold");
  } else if (pin->DoesParameterExist(block_name_, "mdot_chi")) {
    Real rho_hot = 1.0;
    if (pin->DoesBlockExist("problem") &&
        pin->DoesParameterExist("problem", "rho_hot")) {
      rho_hot = pin->GetReal("problem", "rho_hot");
    }
    mdot_density_ = pin->GetReal(block_name_, "mdot_chi")*rho_hot;
  } else {
    mdot_density_ = 1.0;
  }
  const std::string mdot_density_mode = NormalizeToken(mdot_density_mode_name_);
  if (mdot_density_mode == "fixed" || mdot_density_mode == "explicit" ||
      mdot_density_mode == "rho_cold" || mdot_density_mode == "cold") {
    // Keep the density chosen above. This preserves the pre-existing top_mdot
    // and top_mdot_auto defaults unless the user explicitly selects another mode.
  } else if (mdot_density_mode == "geometric_mean" ||
             mdot_density_mode == "sqrt_hot_cold" ||
             mdot_density_mode == "interface") {
    if (!pin->DoesBlockExist("problem") ||
        !pin->DoesParameterExist("problem", "rho_hot") ||
        !pin->DoesParameterExist("problem", "rho_cold")) {
      FatalFrameTrackingInput("mdot_density_mode=geometric_mean requires "
                              "<problem>/rho_hot and <problem>/rho_cold.");
    }
    const Real rho_hot = pin->GetReal("problem", "rho_hot");
    const Real rho_cold = pin->GetReal("problem", "rho_cold");
    if (rho_hot <= 0.0 || rho_cold <= 0.0) {
      FatalFrameTrackingInput("mdot_density_mode=geometric_mean requires "
                              "positive rho_hot and rho_cold.");
    }
    mdot_density_ = std::sqrt(rho_hot*rho_cold);
  } else if (mdot_density_mode == "hot" || mdot_density_mode == "rho_hot") {
    if (!pin->DoesBlockExist("problem") ||
        !pin->DoesParameterExist("problem", "rho_hot")) {
      FatalFrameTrackingInput("mdot_density_mode=hot requires <problem>/rho_hot.");
    }
    mdot_density_ = pin->GetReal("problem", "rho_hot");
  } else {
    FatalFrameTrackingInput("Invalid <frame_tracking>/mdot_density_mode '" +
                            mdot_density_mode_name_ + "'. Expected fixed, "
                            "geometric_mean, rho_cold, or rho_hot.");
  }
  mdot_area_ = pin->GetOrAddReal(block_name_, "mdot_area", -1.0);
  if (mdot_area_ <= 0.0 && pmy_pack != nullptr && pmy_pack->pmesh != nullptr) {
    const RegionSize &mesh_size = pmy_pack->pmesh->mesh_size;
    mdot_area_ = (mesh_size.x1max - mesh_size.x1min)*
                 (mesh_size.x2max - mesh_size.x2min);
  }
  mdot_auto_detect_start_time_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_detect_start_time", 0.5);
  mdot_auto_fast_tau_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_fast_tau", 0.5);
  mdot_auto_slow_tau_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_slow_tau", 1.5);
  mdot_auto_rel_tol_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_rel_tol", 0.075);
  mdot_auto_abs_floor_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_abs_floor", 1.0e-12);
  mdot_auto_hold_time_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_hold_time", 0.5);
  mdot_auto_ramp_fraction_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_ramp_fraction", 0.25);
  mdot_auto_ramp_min_time_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_ramp_min_time", 0.25);
  mdot_auto_ramp_max_time_ =
      pin->GetOrAddReal(block_name_, "mdot_auto_ramp_max_time", 5.0);
  reacquire_expand_factor_ =
      pin->GetOrAddReal(block_name_, "reacquire_expand_factor", 1.0);
  reacquire_max_expand_ = pin->GetOrAddReal(block_name_, "reacquire_max_expand", 1.0);
  reacquire_recover_updates_ =
      pin->GetOrAddInteger(block_name_, "reacquire_recover_updates", 1);
  reacquire_leak_on_miss_ =
      pin->GetOrAddBoolean(block_name_, "reacquire_leak_on_miss", true);
  boundary_guard_enable_ =
      pin->GetOrAddBoolean(block_name_, "boundary_guard_enable", false);
  boundary_guard_cells_ = pin->GetOrAddInteger(block_name_, "boundary_guard_cells", 8);
  boundary_guard_min_scale_ =
      pin->GetOrAddReal(block_name_, "boundary_guard_min_scale", 0.15);

  ParseTrackedFluid(pin);
  ParseAxes(pin);
  ParseAxisControls(pin);
  ValidateConfiguration();
  RestoreFrameState(pin);
  PrintConfigurationSummary();
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::RestoreFrameState()
//! \brief Restore moving-frame state stored in restart input headers.

void FrameTracker::RestoreFrameState(ParameterInput *pin) {
  if (pin == nullptr) {
    return;
  }
  const bool has_state_version = pin->DoesParameterExist(block_name_, "state_version");
  const int state_version =
      has_state_version ? pin->GetInteger(block_name_, "state_version") : 0;
  if (state_version > kFrameTrackingStateVersion) {
    FatalFrameTrackingInput("<frame_tracking>/state_version=" +
                            std::to_string(state_version) +
                            " is newer than this FrameTracker implementation.");
  }
  const bool has_legacy_state =
      pin->DoesParameterExist(block_name_, "frame_velocity_x1") ||
      pin->DoesParameterExist(block_name_, "frame_velocity_x2") ||
      pin->DoesParameterExist(block_name_, "frame_velocity_x3") ||
      pin->DoesParameterExist(block_name_, "frame_displacement_x1") ||
      pin->DoesParameterExist(block_name_, "frame_displacement_x2") ||
      pin->DoesParameterExist(block_name_, "frame_displacement_x3");

  for (int axis = 0; axis < 3; ++axis) {
    const std::string name = AxisName(axis);
    (void) ReadRealAny(pin, block_name_, {"frame_velocity_" + name},
                       frame_velocity_[axis]);
    (void) ReadRealAny(pin, block_name_, {"frame_displacement_" + name},
                       frame_displacement_[axis]);
    axes_[axis].cumulative_boost = -frame_velocity_[axis];
    if (state_version >= 1) {
      AxisState &state = axes_[axis];
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_boost"},
                         state.last_boost);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_cumulative_boost"},
                         state.cumulative_boost);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_mean_x"},
                         state.last_mean_x);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_mean_v"},
                         state.last_mean_v);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_global_weight"},
                         state.last_global_weight);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_filtered_x"},
                         state.last_filtered_x);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_filtered_v"},
                         state.last_filtered_v);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_x_err"},
                         state.last_x_err);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_vel_cmd_pre"},
                         state.last_vel_cmd_pre);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_vel_cmd_post"},
                         state.last_vel_cmd_post);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_x_centroid"},
                         state.last_x_centroid);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_x_midpoint"},
                         state.last_x_midpoint);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_last_x_ctrl"},
                         state.last_x_ctrl);
      (void) ReadRealAny(pin, block_name_, {"state_" + name + "_i_term"}, state.i_term);
      (void) ReadBooleanAny(pin, block_name_, {"state_" + name + "_last_skip_flag"},
                            state.last_skip_flag);
      (void) ReadBooleanAny(pin, block_name_, {"state_" + name + "_last_slew_limited"},
                            state.last_slew_limited);
      (void) ReadBooleanAny(pin, block_name_, {"state_" + name + "_filter_initialized"},
                            state.filter_initialized);
    }
  }
  if (state_version >= 1) {
    restored_state_kind_ = 1;
    (void) ReadRealAny(pin, block_name_, {"state_last_apply_time"}, last_apply_time_);
    (void) ReadIntegerAny(pin, block_name_, {"state_miss_streak"}, miss_streak_);
    (void) ReadIntegerAny(pin, block_name_, {"state_recover_streak"}, recover_streak_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_integral"}, mdot_integral_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_weight_time"}, mdot_weight_time_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_last_top"}, mdot_last_top_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_last_average"},
                       mdot_last_average_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_last_vtop"}, mdot_last_vtop_);
    (void) ReadRealAny(pin, block_name_, {"state_mdot_last_ramp_factor"},
                       mdot_last_ramp_factor_);
    if (state_version >= 3) {
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_fast"},
                         mdot_auto_fast_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_slow"},
                         mdot_auto_slow_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_rel"},
                         mdot_auto_rel_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_hold_elapsed"},
                         mdot_auto_hold_elapsed_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_settle_time"},
                         mdot_auto_settle_time_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_locked_mdot"},
                         mdot_auto_locked_mdot_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_locked_vframe"},
                         mdot_auto_locked_vframe_);
      (void) ReadRealAny(pin, block_name_,
                         {"state_mdot_auto_frame_velocity_at_lock"},
                         mdot_auto_frame_velocity_at_lock_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_ramp_duration"},
                         mdot_auto_ramp_duration_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_target_vframe"},
                         mdot_auto_target_vframe_);
      (void) ReadRealAny(pin, block_name_, {"state_mdot_auto_ramp_factor"},
                         mdot_auto_ramp_factor_);
      (void) ReadBooleanAny(pin, block_name_,
                            {"state_mdot_auto_filter_initialized"},
                            mdot_auto_filter_initialized_);
      (void) ReadBooleanAny(pin, block_name_, {"state_mdot_auto_locked"},
                            mdot_auto_locked_);
    }
  } else if (has_legacy_state) {
    restored_state_kind_ = 2;
    if (global_variable::my_rank == 0) {
      std::cout << "FrameTracker warning: restoring legacy frame state without complete "
                << "controller history; continuation is finite but not restart-exact."
                << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::StoreStateInParameterInput()
//! \brief Store moving-frame state in the input header for restart files.

void FrameTracker::StoreStateInParameterInput(ParameterInput *pin) const {
  if (pin == nullptr) {
    return;
  }
  pin->SetInteger(block_name_, "state_version", kFrameTrackingStateVersion);
  pin->SetRealPrecise(block_name_, "state_last_apply_time", last_apply_time_);
  pin->SetInteger(block_name_, "state_miss_streak", miss_streak_);
  pin->SetInteger(block_name_, "state_recover_streak", recover_streak_);
  pin->SetRealPrecise(block_name_, "state_mdot_integral", mdot_integral_);
  pin->SetRealPrecise(block_name_, "state_mdot_weight_time", mdot_weight_time_);
  pin->SetRealPrecise(block_name_, "state_mdot_last_top", mdot_last_top_);
  pin->SetRealPrecise(block_name_, "state_mdot_last_average", mdot_last_average_);
  pin->SetRealPrecise(block_name_, "state_mdot_last_vtop", mdot_last_vtop_);
  pin->SetRealPrecise(block_name_, "state_mdot_last_ramp_factor",
                      mdot_last_ramp_factor_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_fast", mdot_auto_fast_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_slow", mdot_auto_slow_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_rel", mdot_auto_rel_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_hold_elapsed",
                      mdot_auto_hold_elapsed_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_settle_time",
                      mdot_auto_settle_time_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_locked_mdot",
                      mdot_auto_locked_mdot_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_locked_vframe",
                      mdot_auto_locked_vframe_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_frame_velocity_at_lock",
                      mdot_auto_frame_velocity_at_lock_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_ramp_duration",
                      mdot_auto_ramp_duration_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_target_vframe",
                      mdot_auto_target_vframe_);
  pin->SetRealPrecise(block_name_, "state_mdot_auto_ramp_factor",
                      mdot_auto_ramp_factor_);
  pin->SetBoolean(block_name_, "state_mdot_auto_filter_initialized",
                  mdot_auto_filter_initialized_);
  pin->SetBoolean(block_name_, "state_mdot_auto_locked", mdot_auto_locked_);
  for (int axis = 0; axis < 3; ++axis) {
    const std::string name = AxisName(axis);
    const AxisState &state = axes_[axis];
    pin->SetRealPrecise(block_name_, "frame_velocity_" + name, frame_velocity_[axis]);
    pin->SetRealPrecise(block_name_, "frame_displacement_" + name,
                        frame_displacement_[axis]);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_boost",
                        state.last_boost);
    pin->SetRealPrecise(block_name_, "state_" + name + "_cumulative_boost",
                        state.cumulative_boost);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_mean_x",
                        state.last_mean_x);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_mean_v",
                        state.last_mean_v);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_global_weight",
                        state.last_global_weight);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_filtered_x",
                        state.last_filtered_x);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_filtered_v",
                        state.last_filtered_v);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_x_err", state.last_x_err);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_vel_cmd_pre",
                        state.last_vel_cmd_pre);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_vel_cmd_post",
                        state.last_vel_cmd_post);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_x_centroid",
                        state.last_x_centroid);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_x_midpoint",
                        state.last_x_midpoint);
    pin->SetRealPrecise(block_name_, "state_" + name + "_last_x_ctrl", state.last_x_ctrl);
    pin->SetRealPrecise(block_name_, "state_" + name + "_i_term", state.i_term);
    pin->SetBoolean(block_name_, "state_" + name + "_last_skip_flag",
                    state.last_skip_flag);
    pin->SetBoolean(block_name_, "state_" + name + "_last_slew_limited",
                    state.last_slew_limited);
    pin->SetBoolean(block_name_, "state_" + name + "_filter_initialized",
                    state.filter_initialized);
  }
}

//----------------------------------------------------------------------------------------
//! \fn int FrameTracker::FillHistoryData()
//! \brief Copy controller state into a compact machine-readable history record.

int FrameTracker::FillHistoryData(std::string labels[], Real values[],
                                  const int max_values) const {
  int required = 5;
  for (const auto &state : axes_) {
    if (state.active) {
      required += 5;
    }
  }
  if (mode_ == kFTTopMdot) {
    required += 4;
  } else if (mode_ == kFTTopMdotAuto) {
    required += 11;
  } else if (mode_ == kFTTopMdotWindow) {
    required += 7;
  }
  if (required > max_values) {
    FatalFrameTrackingInput("FrameTracker history data exceed configured output "
                            "capacity.");
  }

  int n = 0;
  const AxisState *shared_state = nullptr;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    const AxisState &state = axes_[axis];
    if (shared_state == nullptr) {
      shared_state = &state;
    }
    const std::string name = AxisName(axis);
    labels[n] = "ft_vf_" + name;
    values[n++] = frame_velocity_[axis];
    labels[n] = "ft_dx_" + name;
    values[n++] = frame_displacement_[axis];
    labels[n] = "ft_pos_" + name;
    values[n++] = state.last_x_ctrl;
    labels[n] = "ft_err_" + name;
    values[n++] = state.last_x_err;
    labels[n] = "ft_dv_" + name;
    values[n++] = state.last_boost;
  }
  labels[n] = "ft_weight";
  values[n++] = (shared_state == nullptr) ? 0.0 : shared_state->last_global_weight;
  labels[n] = "ft_misses";
  values[n++] = static_cast<Real>(miss_streak_);
  labels[n] = "ft_recov";
  values[n++] = static_cast<Real>(recover_streak_);
  bool limited = false;
  bool skipped = false;
  for (const auto &state : axes_) {
    limited = limited || (state.active && state.last_slew_limited);
    skipped = skipped || (state.active && state.last_skip_flag);
  }
  labels[n] = "ft_limit";
  values[n++] = limited ? 1.0 : 0.0;
  labels[n] = "ft_skip";
  values[n++] = skipped ? 1.0 : 0.0;
  if (mode_ == kFTTopMdot) {
    labels[n] = "ft_mdot_top";
    values[n++] = mdot_last_top_;
    labels[n] = "ft_mdot_avg";
    values[n++] = mdot_last_average_;
    labels[n] = "ft_vtop";
    values[n++] = mdot_last_vtop_;
    labels[n] = "ft_ramp";
    values[n++] = mdot_last_ramp_factor_;
  } else if (mode_ == kFTTopMdotAuto) {
    labels[n] = "ft_mtop";
    values[n++] = mdot_last_top_;
    labels[n] = "ft_mfast";
    values[n++] = mdot_auto_fast_;
    labels[n] = "ft_mslow";
    values[n++] = mdot_auto_slow_;
    labels[n] = "ft_mrel";
    values[n++] = mdot_auto_rel_;
    labels[n] = "ft_hold";
    values[n++] = mdot_auto_hold_elapsed_;
    labels[n] = "ft_settle";
    values[n++] = mdot_auto_locked_ ? 1.0 : 0.0;
    labels[n] = "ft_tset";
    values[n++] = mdot_auto_settle_time_;
    labels[n] = "ft_mlock";
    values[n++] = mdot_auto_locked_mdot_;
    labels[n] = "ft_vlock";
    values[n++] = mdot_auto_locked_vframe_;
    labels[n] = "ft_vtgt";
    values[n++] = mdot_auto_target_vframe_;
    labels[n] = "ft_ramp";
    values[n++] = mdot_auto_ramp_factor_;
  } else if (mode_ == kFTTopMdotWindow) {
    labels[n] = "ft_mtop";
    values[n++] = mdot_last_top_;
    labels[n] = "ft_mavg";
    values[n++] = mdot_last_average_;
    labels[n] = "ft_vtgt";
    values[n++] = mdot_window_target_vframe_;
    labels[n] = "ft_wtime";
    values[n++] = mdot_weight_time_;
    labels[n] = "ft_wlo";
    values[n++] = mdot_window_low_time_;
    labels[n] = "ft_rho";
    values[n++] = mdot_density_;
    labels[n] = "ft_fac";
    values[n++] = mdot_velocity_factor_;
  }
  return n;
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::ParseTrackedFluid()
//! \brief Resolve which fluid supplies target samples and receives frame boosts.

void FrameTracker::ParseTrackedFluid(ParameterInput *pin) {
  const bool have_hydro = (pmy_pack != nullptr && pmy_pack->phydro != nullptr);
  const bool have_mhd = (pmy_pack != nullptr && pmy_pack->pmhd != nullptr);
  std::string requested;
  if (ReadStringAny(pin, block_name_, {"tracked_fluid"}, requested)) {
    requested = NormalizeToken(requested);
    if (requested == "hydro" && have_hydro) {
      track_mhd_ = false;
      return;
    }
    if (requested == "mhd" && have_mhd) {
      track_mhd_ = true;
      return;
    }
    FatalFrameTrackingInput("<frame_tracking>/tracked_fluid='" + requested +
                            "' does not select an available hydro or mhd fluid.");
  }
  if (have_hydro && have_mhd) {
    FatalFrameTrackingInput("<frame_tracking> found both <hydro> and <mhd>; set "
                            "tracked_fluid=hydro or tracked_fluid=mhd explicitly.");
  }
  track_mhd_ = have_mhd;
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::RejectRemovedAliases()
//! \brief Reject former parameter spellings so new runs use one stable interface.

void FrameTracker::RejectRemovedAliases(ParameterInput *pin) const {
  if (pin == nullptr) {
    return;
  }
  for (const auto &alias : RemovedAliases()) {
    if (pin->DoesParameterExist(block_name_, alias.first)) {
      FatalFrameTrackingInput("Removed key <" + block_name_ + ">/" + alias.first +
                              "; use " + alias.second + ".");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::PrintConfigurationSummary()
//! \brief Emit one rank-0 startup record of the interpreted active configuration.

void FrameTracker::PrintConfigurationSummary() const {
  if (!enabled_ || global_variable::my_rank != 0) {
    return;
  }
  const char *state_name =
      (restored_state_kind_ == 1) ? "versioned_restart" :
      (restored_state_kind_ == 2) ? "legacy_restart" : "new";
  std::cout << "FrameTracker configuration: fluid=" << (track_mhd_ ? "mhd" : "hydro")
            << " axes=";
  bool first = true;
  for (int axis = 0; axis < 3; ++axis) {
    if (axes_[axis].active) {
      std::cout << (first ? "" : ",") << AxisName(axis);
      first = false;
    }
  }
  std::cout << " target=" << NormalizeToken(target_name_)
            << " range=[" << std::setprecision(16) << target_min_ << ","
            << target_max_ << "] weight="
            << ((weight_mode_ == kFTWeightMass) ? "mass" :
                (weight_mode_ == kFTWeightVolume) ? "volume" :
                (weight_mode_ == kFTWeightTarget) ? "target" : "tracer_mass")
            << " mode=" << NormalizeToken(mode_name_)
            << " position_signal=" << NormalizeToken(position_signal_name_)
            << " velocity_signal=" << NormalizeToken(velocity_signal_name_)
            << " slew=" << NormalizeToken(boost_change_mode_name_)
            << " state=" << state_name << std::endl;
  if (mode_ == kFTTopMdot) {
    std::cout << "FrameTracker top_mdot: average_window=["
              << mdot_average_start_time_ << "," << mdot_average_end_time_
              << "] ramp_time=" << mdot_ramp_time_
              << " density=" << mdot_density_
              << " area=" << mdot_area_
              << " flux_frame=" << NormalizeToken(mdot_flux_frame_name_)
              << std::endl;
    std::cout << "FrameTracker top_mdot sign: Mdot_top=integral rho*v3*dA "
              << "at upper x3 with v3 measured in the "
              << NormalizeToken(mdot_flux_frame_name_)
              << " frame; positive is +x3/outward. "
              << "v_frame_target=Mdot_top/(rho_cold*A); "
              << "fluid_boost=-Delta v_frame." << std::endl;
  } else if (mode_ == kFTTopMdotAuto) {
    std::cout << "FrameTracker top_mdot_auto: detect_start="
              << mdot_auto_detect_start_time_
              << " fast_tau=" << mdot_auto_fast_tau_
              << " slow_tau=" << mdot_auto_slow_tau_
              << " rel_tol=" << mdot_auto_rel_tol_
              << " hold_time=" << mdot_auto_hold_time_
              << " ramp_fraction=" << mdot_auto_ramp_fraction_
              << " ramp_bounds=[" << mdot_auto_ramp_min_time_ << ","
              << mdot_auto_ramp_max_time_ << "] density=" << mdot_density_
              << " area=" << mdot_area_
              << " flux_frame=" << NormalizeToken(mdot_flux_frame_name_)
              << std::endl;
    std::cout << "FrameTracker top_mdot_auto sign: Mdot_top=integral "
              << "rho*v3*dA at upper x3 with v3 measured in the "
              << NormalizeToken(mdot_flux_frame_name_)
              << " frame; positive is +x3/outward. "
              << "v_frame_target=Mdot_top/(rho_cold*A); "
              << "fluid_boost=-Delta v_frame." << std::endl;
  } else if (mode_ == kFTTopMdotWindow) {
    std::cout << "FrameTracker top_mdot_window: start="
              << mdot_window_start_time_
              << " window=" << mdot_window_time_
              << " density_mode=" << NormalizeToken(mdot_density_mode_name_)
              << " density=" << mdot_density_
              << " area=" << mdot_area_
              << " velocity_sign=" << NormalizeToken(mdot_velocity_sign_name_)
              << " velocity_factor=" << mdot_velocity_factor_
              << " flux_frame=" << NormalizeToken(mdot_flux_frame_name_)
              << std::endl;
    std::cout << "FrameTracker top_mdot_window sign: Mdot_top=integral "
              << "rho*v3*dA at upper x3 with v3 measured in the "
              << NormalizeToken(mdot_flux_frame_name_)
              << " frame; positive is +x3/outward. "
              << "v_frame_target=sign*<Mdot_top>_window/(rho_frame*A); "
              << "sign=" << mdot_velocity_sign_
              << " and fluid_boost=-Delta v_frame." << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::ParseAxes()
//! \brief Parse active frame-tracking directions.

void FrameTracker::ParseAxes(ParameterInput *pin) {
  for (auto &axis : axes_) {
    axis.active = false;
  }

  const Mesh *pm = pmy_pack->pmesh;
  std::string axes_raw;
  if (!ReadStringAny(pin, block_name_, {"axes"}, axes_raw)) {
    if (pm->mesh_indcs.nx3 > 1) {
      axes_[2].active = true;
    } else if (pm->mesh_indcs.nx2 > 1) {
      axes_[1].active = true;
    } else {
      axes_[0].active = true;
    }
  } else {
    std::stringstream ss(axes_raw);
    std::string token;
    while (std::getline(ss, token, ',')) {
      token = NormalizeToken(token);
      if (token.empty()) {
        continue;
      }
      if (token == "all" || token == "xyz" || token == "x1x2x3") {
        axes_[0].active = true;
        axes_[1].active = true;
        axes_[2].active = true;
      } else if (token == "x1" || token == "1") {
        axes_[0].active = true;
      } else if (token == "x2" || token == "2") {
        axes_[1].active = true;
      } else if (token == "x3" || token == "3" || token == "z") {
        axes_[2].active = true;
      } else {
        FatalFrameTrackingInput("Invalid <frame_tracking>/axes token '" + token +
                                "'. Use x1, x2, x3, or a comma-separated list.");
      }
    }
  }

  bool active = false;
  for (int axis = 0; axis < 3; ++axis) {
    bool axis_flag = false;
    if (ReadBooleanAny(pin, block_name_, {"track_" + std::string(AxisName(axis))},
                       axis_flag)) {
      axes_[axis].active = axis_flag;
    }
    active = active || axes_[axis].active;
  }
  if (!active) {
    FatalFrameTrackingInput("<frame_tracking> must activate at least one axis.");
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::ParseAxisControls()
//! \brief Parse per-axis target positions and optional predicted-centroid limits.

void FrameTracker::ParseAxisControls(ParameterInput *pin) {
  const Mesh *pm = pmy_pack->pmesh;
  for (int axis = 0; axis < 3; ++axis) {
    const std::string name = AxisName(axis);
    axes_[axis].target_position = 0.5*(AxisGlobalMin(axis, pm) + AxisGlobalMax(axis, pm));
    axes_[axis].lower_limit = AxisGlobalMin(axis, pm);
    axes_[axis].upper_limit = AxisGlobalMax(axis, pm);

    std::vector<std::string> target_keys = {name + "_target"};
    (void) ReadRealAny(pin, block_name_, target_keys, axes_[axis].target_position);

    std::vector<std::string> lower_keys = {name + "_lower_limit"};
    std::vector<std::string> upper_keys = {name + "_upper_limit"};
    if (ReadRealAny(pin, block_name_, lower_keys, axes_[axis].lower_limit)) {
      axes_[axis].use_lower_limit = true;
    }
    if (ReadRealAny(pin, block_name_, upper_keys, axes_[axis].upper_limit)) {
      axes_[axis].use_upper_limit = true;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::ValidateConfiguration()
//! \brief Validate parsed frame-tracking settings against the active fluid.

void FrameTracker::ValidateConfiguration() {
  if (pmy_pack == nullptr || pmy_pack->pmesh == nullptr) {
    FatalFrameTrackingInput("FrameTracker requires a valid MeshBlockPack.");
  }
  if (pmy_pack->phydro == nullptr && pmy_pack->pmhd == nullptr) {
    FatalFrameTrackingInput("<frame_tracking> requires a <hydro> or <mhd> block.");
  }
  if (pmy_pack->pcoord != nullptr &&
      (pmy_pack->pcoord->is_special_relativistic ||
       pmy_pack->pcoord->is_general_relativistic ||
       pmy_pack->pcoord->is_dynamical_relativistic)) {
    FatalFrameTrackingInput("<frame_tracking> supports non-relativistic hydro/MHD only; "
                            "SR, GR, and dynamical-GR coordinates are unsupported.");
  }

  if (apply_every_ < 1) {
    FatalFrameTrackingInput("Require <frame_tracking>/apply_every >= 1.");
  }
  if (tau_avg_ <= 0.0 || tau_relax_ <= 0.0 || tau_vel_ <= 0.0) {
    FatalFrameTrackingInput("Frame-tracking tau_avg, tau_relax, and tau_vel "
                            "must be positive.");
  }
  if (target_max_ < target_min_) {
    FatalFrameTrackingInput("Require <frame_tracking>/target_max >= target_min.");
  }
  if (weight_floor_ < 0.0 || min_global_weight_ < 0.0 || max_abs_boost_ < 0.0 ||
      max_boost_change_ < 0.0 || int_max_abs_ < 0.0 || mdot_ramp_time_ < 0.0 ||
      mdot_auto_detect_start_time_ < 0.0 || mdot_auto_abs_floor_ < 0.0 ||
      mdot_auto_hold_time_ < 0.0 || mdot_auto_ramp_fraction_ < 0.0 ||
      mdot_auto_ramp_min_time_ < 0.0 || mdot_auto_ramp_max_time_ < 0.0 ||
      mdot_window_start_time_ < 0.0 || mdot_window_time_ <= 0.0 ||
      mdot_velocity_factor_ <= 0.0) {
    FatalFrameTrackingInput("Frame-tracking floors and boost/integral limits "
                            "must be non-negative; mdot_window_time and "
                            "mdot_velocity_factor must be positive.");
  }
  if (mdot_auto_fast_tau_ <= 0.0 || mdot_auto_slow_tau_ <= 0.0 ||
      mdot_auto_rel_tol_ < 0.0) {
    FatalFrameTrackingInput("Require positive mdot_auto_fast_tau and "
                            "mdot_auto_slow_tau, and non-negative "
                            "mdot_auto_rel_tol.");
  }
  if (mdot_auto_ramp_max_time_ < mdot_auto_ramp_min_time_) {
    FatalFrameTrackingInput("Require mdot_auto_ramp_max_time >= "
                            "mdot_auto_ramp_min_time.");
  }
  if (boost_change_mode_ == kFTBoostPerTime && max_boost_change_rate_ <= 0.0) {
    FatalFrameTrackingInput("Require max_boost_change_rate > 0 when "
                            "max_boost_change_mode=per_time.");
  }
  if (tau_int_ > 0.0 && int_leak_tau_ <= 0.0) {
    FatalFrameTrackingInput("Require int_leak_tau > 0 when tau_int > 0.");
  }
  if (position_blend_ < 0.0 || position_blend_ > 1.0) {
    FatalFrameTrackingInput("Require 0 <= position_blend <= 1.");
  }
  if (reacquire_expand_factor_ < 1.0 || reacquire_max_expand_ < 1.0 ||
      reacquire_recover_updates_ < 1) {
    FatalFrameTrackingInput("Require reacquire_expand_factor >= 1, "
                            "reacquire_max_expand >= 1, and "
                            "reacquire_recover_updates >= 1.");
  }
  if (boundary_guard_cells_ < 1 || boundary_guard_min_scale_ < 0.0 ||
      boundary_guard_min_scale_ > 1.0) {
    FatalFrameTrackingInput("Require boundary_guard_cells >= 1 and "
                            "0 <= boundary_guard_min_scale <= 1.");
  }
  if (mode_ == kFTTopMdot || mode_ == kFTTopMdotAuto ||
      mode_ == kFTTopMdotWindow) {
    if (!axes_[2].active || axes_[0].active || axes_[1].active) {
      FatalFrameTrackingInput("<frame_tracking>/mode=top_mdot, top_mdot_auto, "
                              "or top_mdot_window supports only axes=x3.");
    }
    if (mdot_density_ <= 0.0 || mdot_area_ <= 0.0) {
      FatalFrameTrackingInput("Require positive mdot_density and mdot_area for "
                              "<frame_tracking>/mode=top_mdot, top_mdot_auto, "
                              "or top_mdot_window.");
    }
    if (mode_ == kFTTopMdot &&
        mdot_average_end_time_ >= mdot_average_start_time_ &&
        mdot_average_end_time_ <= 0.0) {
      FatalFrameTrackingInput("Require mdot_average_end_time > 0, or set it below "
                              "mdot_average_start_time for an open-ended average.");
    }
  }

  const RegionIndcs &mesh_indcs = pmy_pack->pmesh->mesh_indcs;
  for (int axis = 0; axis < 3; ++axis) {
    if (axes_[axis].active && AxisCellCount(axis, mesh_indcs) <= 1 && axis != 0) {
      FatalFrameTrackingInput(std::string("<frame_tracking> requested inactive axis ") +
                              AxisName(axis) + ".");
    }
    if ((axes_[axis].use_lower_limit || axes_[axis].use_upper_limit) &&
        axes_[axis].lower_limit > axes_[axis].upper_limit) {
      FatalFrameTrackingInput(std::string("Require ") + AxisName(axis) +
                              "_lower_limit <= " + AxisName(axis) + "_upper_limit.");
    }
  }

  const bool is_mhd = track_mhd_;
  const int nscalars = is_mhd ? pmy_pack->pmhd->nscalars :
                                pmy_pack->phydro->nscalars;
  const EOS_Data &eos = is_mhd ? pmy_pack->pmhd->peos->eos_data :
                                 pmy_pack->phydro->peos->eos_data;
  if (target_kind_ == kFTTargetScalar &&
      (scalar_index_ < 0 || scalar_index_ >= nscalars)) {
    FatalFrameTrackingInput("<frame_tracking> scalar target requested scalar_index=" +
                            std::to_string(scalar_index_) +
                            ", but the active fluid has " + std::to_string(nscalars) +
                            " passive scalar(s).");
  }
  if (weight_mode_ == kFTWeightTracerMass && target_kind_ != kFTTargetScalar) {
    FatalFrameTrackingInput("<frame_tracking>/weight=tracer_mass requires "
                            "target=scalar or target=scalarN.");
  }
  if (target_kind_ == kFTTargetEntropy && !eos.is_ideal) {
    FatalFrameTrackingInput("<frame_tracking>/target=entropy requires an ideal-gas EOS.");
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskID FrameTracker::IncludeFrameTrackingTask()
//! \brief Add frame tracking to the after-time-integrator task list.

TaskID FrameTracker::IncludeFrameTrackingTask(std::shared_ptr<TaskList> tl,
                                              TaskID start) {
  if (tl == nullptr) {
    return TaskID(0);
  }
  TaskID dep = start;
  return tl->AddTask(&FrameTracker::Apply, this, dep);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus FrameTracker::Apply()
//! \brief Task-list wrapper for scheduled frame tracking.

TaskStatus FrameTracker::Apply(Driver *pdrive, int stage) {
  if (!enabled_ || pmy_pack == nullptr || pmy_pack->pmesh == nullptr) {
    return TaskStatus::complete;
  }
  Mesh *pm = pmy_pack->pmesh;
  const bool displacement_changed = AdvanceFrameDisplacement(pm->dt);
  bool boost_changed = false;
  if (((mode_ == kFTTopMdot) || (mode_ == kFTTopMdotAuto) ||
       pm->time > start_time_) &&
      pm->ncycle % apply_every_ == 0) {
    boost_changed = ApplyTracking();
  }

  if ((displacement_changed || boost_changed) && !pm->strictly_periodic) {
    // User and physical ghost states may depend on the lab-frame origin and velocity.
    // Refresh them before output and before the next cycle consumes ghost-zone fluxes.
    if (track_mhd_) {
      (void) pmy_pack->pmhd->ApplyPhysicalBCs(pdrive, stage);
      (void) pmy_pack->pmhd->ConToPrimGhostZones(pdrive, stage);
    } else {
      (void) pmy_pack->phydro->ApplyPhysicalBCs(pdrive, stage);
      (void) pmy_pack->phydro->ConToPrimGhostZones(pdrive, stage);
    }
  }

  if (boost_changed) {
    // The controller changes characteristic velocities after the usual fluid dt task.
    // Refresh dtnew from the boosted state so restarts follow the same next timestep.
    if (track_mhd_) {
      (void) pmy_pack->pmhd->NewTimeStep(pdrive, pdrive->nexp_stages);
    } else {
      (void) pmy_pack->phydro->NewTimeStep(pdrive, pdrive->nexp_stages);
    }
  }
  return TaskStatus::complete;
}

//! \fn void FrameTracker::SetActiveTargetRange()
//! \brief Expand the tracked target interval while reacquiring missed material.

void FrameTracker::SetActiveTargetRange() {
  const bool reacquire_enabled =
      (reacquire_expand_factor_ > 1.0) && (reacquire_max_expand_ > 1.0);
  Real expand = 1.0;
  if (reacquire_enabled && miss_streak_ > 0) {
    expand = std::pow(reacquire_expand_factor_, static_cast<Real>(miss_streak_));
    expand = std::min(expand, reacquire_max_expand_);
  }

  if (use_log_reacquire_) {
    const Real safe_lo = std::max(target_min_, static_cast<Real>(1.0e-30));
    const Real safe_hi = std::max(target_max_, safe_lo*(1.0 + 1.0e-8));
    const Real safe_center = std::max(target_center_, static_cast<Real>(1.0e-30));
    const Real log_center = std::log10(safe_center);
    const Real log_lo = std::log10(safe_lo);
    const Real log_hi = std::log10(safe_hi);
    const Real dlog_lo = std::max(static_cast<Real>(0.0), log_center - log_lo);
    const Real dlog_hi = std::max(static_cast<Real>(0.0), log_hi - log_center);
    target_min_active_ = std::pow(static_cast<Real>(10.0),
                                  log_center - dlog_lo*expand);
    target_max_active_ = std::pow(static_cast<Real>(10.0),
                                  log_center + dlog_hi*expand);
  } else {
    target_min_active_ = target_center_ - (target_center_ - target_min_)*expand;
    target_max_active_ = target_center_ + (target_max_ - target_center_)*expand;
  }
  if (target_max_active_ < target_min_active_) {
    std::swap(target_min_active_, target_max_active_);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::SampleMoments()
//! \brief Compute global weighted moments for all active axes in one mesh traversal.

bool FrameTracker::SampleMoments(std::array<MomentSample, 3> &samples) const {
  for (auto &sample : samples) {
    sample = MomentSample();
  }
  MeshBlockPack *pmbp = pmy_pack;
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  int nmb = pmbp->nmb_thispack;
  if (nmb <= 0) {
    return false;
  }

  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  const bool is_mhd = track_mhd_;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const EOS_Data eos = is_mhd ? pmbp->pmhd->peos->eos_data :
                                pmbp->phydro->peos->eos_data;
  const int nfluid = is_mhd ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  const int scalar_var = nfluid + scalar_index_;
  auto &size = pmbp->pmb->mb_size;
  const int target_kind = target_kind_;
  const int weight_mode = weight_mode_;
  const Real target_min = target_min_active_;
  const Real target_max = target_max_active_;
  const bool active_x1 = axes_[0].active;
  const bool active_x2 = axes_[1].active;
  const bool active_x3 = axes_[2].active;
  Real rank_weight = 0.0;
  Real rank_weighted_x1 = 0.0, rank_weighted_x2 = 0.0, rank_weighted_x3 = 0.0;
  Real rank_weighted_v1 = 0.0, rank_weighted_v2 = 0.0, rank_weighted_v3 = 0.0;
  Real rank_xmin1 = FLT_MAX, rank_xmin2 = FLT_MAX, rank_xmin3 = FLT_MAX;
  Real rank_xmax1 = -FLT_MAX, rank_xmax2 = -FLT_MAX, rank_xmax3 = -FLT_MAX;
  int rank_count = 0;

  Kokkos::parallel_reduce("frame_tracking_moments",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &weight_,
                Real &weighted_x1_, Real &weighted_x2_, Real &weighted_x3_,
                Real &weighted_v1_, Real &weighted_v2_, Real &weighted_v3_,
                Real &xmin1_, Real &xmin2_, Real &xmin3_,
                Real &xmax1_, Real &xmax2_, Real &xmax3_, int &count_) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    const Real dens = w0(m, IDN, k, j, i);
    if (dens <= 0.0) {
      return;
    }

    const Real value = FrameTrackingTargetValue(target_kind, scalar_var, w0, eos,
                                                m, k, j, i);
    if (value < target_min || value > target_max) {
      return;
    }

    const RegionSize block_size = size.d_view(m);
    const Real dvol = block_size.dx1*block_size.dx2*block_size.dx3;
    Real weight = dens*dvol;
    if (weight_mode == kFTWeightVolume) {
      weight = dvol;
    } else if (weight_mode == kFTWeightTarget) {
      weight = fabs(value)*dvol;
    } else if (weight_mode == kFTWeightTracerMass) {
      weight = dens*fmax(value, 0.0)*dvol;
    }
    if (weight <= 0.0) {
      return;
    }

    weight_ += weight;
    count_ += 1;
    if (active_x1) {
      const Real x = CellCenterX(i - is, nx1, block_size.x1min, block_size.x1max);
      weighted_x1_ += weight*x;
      weighted_v1_ += weight*w0(m, IVX, k, j, i);
      xmin1_ = fmin(xmin1_, x);
      xmax1_ = fmax(xmax1_, x);
    }
    if (active_x2) {
      const Real x = CellCenterX(j - js, nx2, block_size.x2min, block_size.x2max);
      weighted_x2_ += weight*x;
      weighted_v2_ += weight*w0(m, IVY, k, j, i);
      xmin2_ = fmin(xmin2_, x);
      xmax2_ = fmax(xmax2_, x);
    }
    if (active_x3) {
      const Real x = CellCenterX(k - ks, nx3, block_size.x3min, block_size.x3max);
      weighted_x3_ += weight*x;
      weighted_v3_ += weight*w0(m, IVZ, k, j, i);
      xmin3_ = fmin(xmin3_, x);
      xmax3_ = fmax(xmax3_, x);
    }
  }, Kokkos::Sum<Real>(rank_weight),
     Kokkos::Sum<Real>(rank_weighted_x1), Kokkos::Sum<Real>(rank_weighted_x2),
     Kokkos::Sum<Real>(rank_weighted_x3), Kokkos::Sum<Real>(rank_weighted_v1),
     Kokkos::Sum<Real>(rank_weighted_v2), Kokkos::Sum<Real>(rank_weighted_v3),
     Kokkos::Min<Real>(rank_xmin1), Kokkos::Min<Real>(rank_xmin2),
     Kokkos::Min<Real>(rank_xmin3), Kokkos::Max<Real>(rank_xmax1),
     Kokkos::Max<Real>(rank_xmax2), Kokkos::Max<Real>(rank_xmax3),
     Kokkos::Sum<int>(rank_count));

  Real global_sums[7] = {rank_weight, rank_weighted_x1, rank_weighted_x2,
                         rank_weighted_x3, rank_weighted_v1, rank_weighted_v2,
                         rank_weighted_v3};
  Real global_xmin[3] = {rank_xmin1, rank_xmin2, rank_xmin3};
  Real global_xmax[3] = {rank_xmax1, rank_xmax2, rank_xmax3};
  int global_count = rank_count;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, global_sums, 7, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, global_xmin, 3, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, global_xmax, 3, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (global_sums[0] <= 0.0 || global_count <= 0) {
    return false;
  }
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    samples[axis].global_weight = global_sums[0];
    samples[axis].mean_x = global_sums[1 + axis]/global_sums[0];
    samples[axis].mean_v = global_sums[4 + axis]/global_sums[0];
    samples[axis].x_min = global_xmin[axis];
    samples[axis].x_max = global_xmax[axis];
    samples[axis].valid = true;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn Real FrameTracker::SampleTopMassFlux()
//! \brief Compute the signed mass flux through the upper x3 active-cell layer.

Real FrameTracker::SampleTopMassFlux(const bool lab_frame) const {
  MeshBlockPack *pmbp = pmy_pack;
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ke = indcs.ke;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nmb = pmbp->nmb_thispack;
  if (nmb <= 0) {
    return 0.0;
  }

  const bool is_mhd = track_mhd_;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &size = pmbp->pmb->mb_size;
  const Real mesh_x3max = pm->mesh_size.x3max;
  const int nmji = nmb*nx2*nx1;
  const int nji = nx2*nx1;
  const Real frame_v3 = lab_frame ? frame_velocity_[2] : 0.0;
  Real rank_flux = 0.0;

  Kokkos::parallel_reduce("frame_tracking_top_mdot",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmji),
  KOKKOS_LAMBDA(const int idx, Real &flux_) {
    const int m = idx/nji;
    const int j = ((idx - m*nji)/nx1) + js;
    const int i = (idx - m*nji - (j - js)*nx1) + is;
    const RegionSize block_size = size.d_view(m);
    if (fabs(block_size.x3max - mesh_x3max) >
        static_cast<Real>(0.5)*block_size.dx3) {
      return;
    }
    // Signed top flux convention:
    //   Mdot_top = integral rho*v3*dA in the upper active-cell layer.
    // Positive is +x3/outward through the top boundary; negative is -x3/inward.
    // Stored velocities are grid-frame velocities.  If lab_frame is requested,
    // add the frame velocity so this samples the physical/lab mass flux.
    const Real dA = block_size.dx1*block_size.dx2;
    flux_ += w0(m, IDN, ke, j, i)*(w0(m, IVZ, ke, j, i) + frame_v3)*dA;
  }, Kokkos::Sum<Real>(rank_flux));

  Real global_flux = rank_flux;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &global_flux, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  return global_flux;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::ApplyX3FrameVelocityTarget()
//! \brief Move the x3 frame velocity toward a target lab-frame velocity.

bool FrameTracker::ApplyX3FrameVelocityTarget(const Real target_frame_velocity,
                                              const Real dt_boost) {
  Mesh *pm = pmy_pack->pmesh;
  AxisState &state = axes_[2];

  // Stored fluid velocities are frame-relative. If the lab-frame coordinate frame
  // changes by Delta v_frame, every stored fluid velocity receives the opposite
  // Galilean boost:
  //   fluid_boost = -Delta v_frame.
  const Real requested_frame_delta = target_frame_velocity - frame_velocity_[2];
  Real frame_delta = requested_frame_delta;
  state.last_slew_limited = false;

  Real dv_allowed = 0.0;
  if (boost_change_mode_ == kFTBoostPerApply) {
    dv_allowed = max_boost_change_;
  } else if (boost_change_mode_ == kFTBoostPerTime && max_boost_change_rate_ > 0.0) {
    dv_allowed = max_boost_change_rate_*dt_boost;
  }
  if (dv_allowed > 0.0 && std::fabs(frame_delta) > dv_allowed) {
    frame_delta = std::copysign(dv_allowed, frame_delta);
    state.last_slew_limited = true;
  }

  Real boost = -frame_delta;
  if (max_abs_boost_ > 0.0 && std::fabs(boost) > max_abs_boost_) {
    boost = std::copysign(max_abs_boost_, boost);
    frame_delta = -boost;
    state.last_slew_limited = true;
  }

  state.last_vel_cmd_pre = -requested_frame_delta;
  state.last_vel_cmd_post = boost;
  state.last_boost = boost;
  state.cumulative_boost += boost;
  frame_velocity_[2] = -state.cumulative_boost;
  state.last_x_err = frame_velocity_[2] - target_frame_velocity;
  state.last_skip_flag = (boost == 0.0);

  const bool is_mhd = track_mhd_;
  auto &u0 = is_mhd ? pmy_pack->pmhd->u0 : pmy_pack->phydro->u0;
  auto &w0 = is_mhd ? pmy_pack->pmhd->w0 : pmy_pack->phydro->w0;
  const EOS_Data eos = is_mhd ? pmy_pack->pmhd->peos->eos_data :
                                pmy_pack->phydro->peos->eos_data;
  const bool use_energy = eos.is_ideal;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool state_changed = (boost != 0.0);
  if (state_changed && nmb1 >= 0) {
    par_for("frame_tracking_x3_target_boost", DevExeSpace(), 0, nmb1, 0, n3-1,
            0, n2-1, 0, n1-1, KOKKOS_LAMBDA(const int m, const int k,
                                             const int j, const int i) {
      const Real rho = w0(m, IDN, k, j, i);
      if (rho <= 0.0) {
        return;
      }
      const Real mom_old3 = u0(m, IM3, k, j, i);
      w0(m, IVZ, k, j, i) += boost;
      u0(m, IM3, k, j, i) += rho*boost;
      if (use_energy) {
        u0(m, IEN, k, j, i) += mom_old3*boost +
                                static_cast<Real>(0.5)*rho*SQR(boost);
      }
    });
  }
  return state_changed;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::ApplyTopMdotTracking()
//! \brief Drive the x3 frame velocity from the time-averaged top mass flux only.

bool FrameTracker::ApplyTopMdotTracking(const Real dt_boost) {
  Mesh *pm = pmy_pack->pmesh;
  if (!axes_[2].active) {
    return false;
  }

  const Real top_mdot = SampleTopMassFlux(mdot_flux_lab_frame_);
  mdot_last_top_ = top_mdot;

  const Real interval_end = pm->time;
  Real interval_start = (last_apply_time_ >= 0.0) ? last_apply_time_ :
      (interval_end - dt_boost);
  interval_start = std::max(static_cast<Real>(0.0), interval_start);
  const Real avg_end = (mdot_average_end_time_ >= mdot_average_start_time_) ?
      mdot_average_end_time_ : interval_end;
  const Real overlap_start = std::max(interval_start, mdot_average_start_time_);
  const Real overlap_end = std::min(interval_end, avg_end);
  const Real overlap = std::max(static_cast<Real>(0.0), overlap_end - overlap_start);
  if (overlap > 0.0) {
    mdot_integral_ += top_mdot*overlap;
    mdot_weight_time_ += overlap;
  }

  mdot_last_average_ =
      (mdot_weight_time_ > 0.0) ? (mdot_integral_/mdot_weight_time_) : 0.0;
  const Real raw_vtop = (mdot_area_ > 0.0 && mdot_density_ > 0.0) ?
      (mdot_last_average_/(mdot_area_*mdot_density_)) : 0.0;

  if (pm->time <= start_time_) {
    mdot_last_ramp_factor_ = 0.0;
  } else if (mdot_ramp_time_ > 0.0) {
    mdot_last_ramp_factor_ = SmoothStep01((pm->time - start_time_)/mdot_ramp_time_);
  } else {
    mdot_last_ramp_factor_ = 1.0;
  }
  mdot_last_vtop_ = mdot_last_ramp_factor_*raw_vtop;

  AxisState &state = axes_[2];
  state.last_global_weight = mdot_weight_time_;
  state.last_mean_x = 0.0;
  state.last_mean_v = mdot_last_average_;
  state.last_filtered_x = 0.0;
  state.last_filtered_v = mdot_last_vtop_;
  state.last_x_centroid = 0.0;
  state.last_x_midpoint = 0.0;
  state.last_x_ctrl = 0.0;
  state.last_x_err = frame_velocity_[2] - mdot_last_vtop_;
  state.last_skip_flag = (pm->time <= start_time_ || mdot_weight_time_ <= 0.0);
  state.last_slew_limited = false;

  if (state.last_skip_flag) {
    state.last_vel_cmd_pre = 0.0;
    state.last_vel_cmd_post = 0.0;
    state.last_boost = 0.0;
    last_apply_time_ = pm->time;
    PrintTopMdotMessage(top_mdot, mdot_last_average_, mdot_last_vtop_,
                        mdot_last_ramp_factor_, 0.0);
    return false;
  }

  const bool state_changed = ApplyX3FrameVelocityTarget(mdot_last_vtop_, dt_boost);
  state.last_skip_flag = false;  // Preserve the existing fixed-window top_mdot semantics.
  last_apply_time_ = pm->time;
  PrintTopMdotMessage(top_mdot, mdot_last_average_, mdot_last_vtop_,
                      mdot_last_ramp_factor_, state.last_boost);
  return state_changed;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::ApplyTopMdotAutoTracking()
//! \brief Auto-lock the top-Mdot frame velocity after flux EMAs settle.

bool FrameTracker::ApplyTopMdotAutoTracking(const Real dt_boost) {
  Mesh *pm = pmy_pack->pmesh;
  if (!axes_[2].active) {
    return false;
  }

  const Real top_mdot = SampleTopMassFlux(mdot_flux_lab_frame_);
  mdot_last_top_ = top_mdot;

  AxisState &state = axes_[2];
  state.last_global_weight = mdot_weight_time_;
  state.last_mean_x = 0.0;
  state.last_mean_v = mdot_auto_slow_;
  state.last_filtered_x = 0.0;
  state.last_filtered_v = mdot_auto_target_vframe_;
  state.last_x_centroid = 0.0;
  state.last_x_midpoint = 0.0;
  state.last_x_ctrl = 0.0;
  state.last_x_err = frame_velocity_[2] - mdot_auto_target_vframe_;
  state.last_slew_limited = false;

  if (pm->time < mdot_auto_detect_start_time_) {
    mdot_auto_ramp_factor_ = 0.0;
    mdot_auto_target_vframe_ = frame_velocity_[2];
    mdot_last_average_ = mdot_auto_slow_;
    mdot_last_vtop_ = mdot_auto_target_vframe_;
    mdot_last_ramp_factor_ = mdot_auto_ramp_factor_;
    state.last_filtered_v = mdot_auto_target_vframe_;
    state.last_x_err = frame_velocity_[2] - mdot_auto_target_vframe_;
    state.last_skip_flag = true;
    state.last_vel_cmd_pre = 0.0;
    state.last_vel_cmd_post = 0.0;
    state.last_boost = 0.0;
    last_apply_time_ = pm->time;
    PrintTopMdotAutoMessage(0.0);
    return false;
  }

  if (!mdot_auto_locked_) {
    if (!mdot_auto_filter_initialized_) {
      mdot_auto_fast_ = top_mdot;
      mdot_auto_slow_ = top_mdot;
      mdot_auto_rel_ = 0.0;
      mdot_auto_hold_elapsed_ = 0.0;
      mdot_auto_filter_initialized_ = true;
      mdot_auto_ramp_factor_ = 0.0;
      mdot_auto_target_vframe_ = frame_velocity_[2];
      mdot_last_average_ = mdot_auto_slow_;
      mdot_last_vtop_ = mdot_auto_target_vframe_;
      mdot_last_ramp_factor_ = mdot_auto_ramp_factor_;
      state.last_skip_flag = true;
      state.last_vel_cmd_pre = 0.0;
      state.last_vel_cmd_post = 0.0;
      state.last_boost = 0.0;
      state.last_global_weight = mdot_weight_time_;
      state.last_mean_v = mdot_auto_slow_;
      state.last_filtered_v = mdot_auto_target_vframe_;
      state.last_x_err = frame_velocity_[2] - mdot_auto_target_vframe_;
      last_apply_time_ = pm->time;
      PrintTopMdotAutoMessage(0.0);
      return false;
    }

    Real alpha_fast = 1.0 - std::exp(-dt_boost/mdot_auto_fast_tau_);
    Real alpha_slow = 1.0 - std::exp(-dt_boost/mdot_auto_slow_tau_);
    alpha_fast = std::max(static_cast<Real>(0.0),
                          std::min(static_cast<Real>(1.0), alpha_fast));
    alpha_slow = std::max(static_cast<Real>(0.0),
                          std::min(static_cast<Real>(1.0), alpha_slow));
    mdot_auto_fast_ += alpha_fast*(top_mdot - mdot_auto_fast_);
    mdot_auto_slow_ += alpha_slow*(top_mdot - mdot_auto_slow_);
    mdot_weight_time_ += dt_boost;
    const Real rel_denom = std::max(std::fabs(mdot_auto_slow_),
                                    mdot_auto_abs_floor_);
    mdot_auto_rel_ = std::fabs(mdot_auto_fast_ - mdot_auto_slow_)/rel_denom;

    if (mdot_auto_rel_ <= mdot_auto_rel_tol_) {
      mdot_auto_hold_elapsed_ += dt_boost;
    } else {
      mdot_auto_hold_elapsed_ = 0.0;
    }

    if (mdot_auto_hold_elapsed_ >= mdot_auto_hold_time_) {
      mdot_auto_locked_ = true;
      mdot_auto_settle_time_ = pm->time;
      mdot_auto_locked_mdot_ = mdot_auto_slow_;
      mdot_auto_locked_vframe_ = mdot_auto_locked_mdot_/(mdot_density_*mdot_area_);
      mdot_auto_frame_velocity_at_lock_ = frame_velocity_[2];
      mdot_auto_ramp_duration_ =
          mdot_auto_ramp_fraction_*std::max(static_cast<Real>(0.0),
                                            mdot_auto_settle_time_);
      mdot_auto_ramp_duration_ =
          std::max(mdot_auto_ramp_min_time_,
                   std::min(mdot_auto_ramp_max_time_, mdot_auto_ramp_duration_));
    }
  }

  bool state_changed = false;
  if (mdot_auto_locked_) {
    if (mdot_auto_ramp_duration_ > 0.0) {
      mdot_auto_ramp_factor_ =
          SmoothStep01((pm->time - mdot_auto_settle_time_)/
                       mdot_auto_ramp_duration_);
    } else {
      mdot_auto_ramp_factor_ = 1.0;
    }
    mdot_auto_target_vframe_ =
        mdot_auto_frame_velocity_at_lock_ +
        mdot_auto_ramp_factor_*
        (mdot_auto_locked_vframe_ - mdot_auto_frame_velocity_at_lock_);
    state.last_filtered_v = mdot_auto_target_vframe_;
    state.last_mean_v = mdot_auto_slow_;
    state.last_global_weight = mdot_weight_time_;
    state_changed = ApplyX3FrameVelocityTarget(mdot_auto_target_vframe_, dt_boost);
    state.last_skip_flag = !state_changed;
  } else {
    mdot_auto_ramp_factor_ = 0.0;
    mdot_auto_target_vframe_ = frame_velocity_[2];
    state.last_filtered_v = mdot_auto_target_vframe_;
    state.last_mean_v = mdot_auto_slow_;
    state.last_global_weight = mdot_weight_time_;
    state.last_x_err = frame_velocity_[2] - mdot_auto_target_vframe_;
    state.last_skip_flag = true;
    state.last_vel_cmd_pre = 0.0;
    state.last_vel_cmd_post = 0.0;
    state.last_boost = 0.0;
  }

  mdot_last_average_ = mdot_auto_slow_;
  mdot_last_vtop_ = mdot_auto_target_vframe_;
  mdot_last_ramp_factor_ = mdot_auto_ramp_factor_;
  last_apply_time_ = pm->time;
  PrintTopMdotAutoMessage(state.last_boost);
  return state_changed;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::ApplyTopMdotWindowTracking()
//! \brief Drive x3 frame velocity from a finite-memory signed top-Mdot average.

bool FrameTracker::ApplyTopMdotWindowTracking(const Real dt_boost) {
  Mesh *pm = pmy_pack->pmesh;
  if (!axes_[2].active) {
    return false;
  }

  const Real top_mdot = SampleTopMassFlux(mdot_flux_lab_frame_);
  mdot_last_top_ = top_mdot;

  AxisState &state = axes_[2];
  state.last_mean_x = 0.0;
  state.last_filtered_x = 0.0;
  state.last_x_centroid = 0.0;
  state.last_x_midpoint = 0.0;
  state.last_x_ctrl = 0.0;
  state.last_slew_limited = false;

  if (pm->time < mdot_window_start_time_) {
    mdot_integral_ = 0.0;
    mdot_weight_time_ = 0.0;
    mdot_last_average_ = 0.0;
    mdot_window_low_time_ = mdot_window_start_time_;
    mdot_window_target_vframe_ = frame_velocity_[2];
    mdot_last_vtop_ = mdot_window_target_vframe_;
    mdot_last_ramp_factor_ = 0.0;
    state.last_global_weight = 0.0;
    state.last_mean_v = mdot_last_average_;
    state.last_filtered_v = mdot_window_target_vframe_;
    state.last_x_err = frame_velocity_[2] - mdot_window_target_vframe_;
    state.last_skip_flag = true;
    state.last_vel_cmd_pre = 0.0;
    state.last_vel_cmd_post = 0.0;
    state.last_boost = 0.0;
    last_apply_time_ = pm->time;
    PrintTopMdotWindowMessage(0.0);
    return false;
  }

  if (!mdot_window_times_.empty() && pm->time <= mdot_window_times_.back()) {
    mdot_window_times_.back() = pm->time;
    mdot_window_values_.back() = top_mdot;
  } else {
    mdot_window_times_.push_back(pm->time);
    mdot_window_values_.push_back(top_mdot);
  }

  const Real window_low =
      std::max(mdot_window_start_time_, pm->time - mdot_window_time_);
  mdot_window_low_time_ = window_low;
  while (mdot_window_times_.size() > 2 &&
         mdot_window_times_[1] <= window_low) {
    mdot_window_times_.pop_front();
    mdot_window_values_.pop_front();
  }

  auto interpolate_mdot = [&](const Real t_sample) {
    if (mdot_window_times_.empty()) {
      return top_mdot;
    }
    if (t_sample <= mdot_window_times_.front()) {
      return mdot_window_values_.front();
    }
    for (std::size_t n = 1; n < mdot_window_times_.size(); ++n) {
      if (t_sample <= mdot_window_times_[n]) {
        const Real t0 = mdot_window_times_[n - 1];
        const Real t1 = mdot_window_times_[n];
        const Real y0 = mdot_window_values_[n - 1];
        const Real y1 = mdot_window_values_[n];
        if (t1 <= t0) {
          return y1;
        }
        const Real f = (t_sample - t0)/(t1 - t0);
        return y0 + f*(y1 - y0);
      }
    }
    return mdot_window_values_.back();
  };

  const Real duration = std::max(static_cast<Real>(0.0), pm->time - window_low);
  Real integral = 0.0;
  if (duration > 0.0) {
    Real prev_time = window_low;
    Real prev_mdot = interpolate_mdot(window_low);
    for (std::size_t n = 0; n < mdot_window_times_.size(); ++n) {
      const Real sample_time = mdot_window_times_[n];
      if (sample_time <= window_low) {
        continue;
      }
      if (sample_time >= pm->time) {
        break;
      }
      const Real sample_mdot = mdot_window_values_[n];
      integral += static_cast<Real>(0.5)*(prev_mdot + sample_mdot)*
                  (sample_time - prev_time);
      prev_time = sample_time;
      prev_mdot = sample_mdot;
    }
    integral += static_cast<Real>(0.5)*(prev_mdot + top_mdot)*
                (pm->time - prev_time);
  }

  mdot_integral_ = integral;
  mdot_weight_time_ = duration;
  mdot_last_average_ = (duration > 0.0) ? (integral/duration) : 0.0;
  mdot_window_target_vframe_ =
      mdot_velocity_factor_*mdot_velocity_sign_*mdot_last_average_/
      (mdot_density_*mdot_area_);
  mdot_last_vtop_ = mdot_window_target_vframe_;
  mdot_last_ramp_factor_ = (duration > 0.0) ? 1.0 : 0.0;

  state.last_global_weight = mdot_weight_time_;
  state.last_mean_v = mdot_last_average_;
  state.last_filtered_v = mdot_window_target_vframe_;
  state.last_x_err = frame_velocity_[2] - mdot_window_target_vframe_;
  state.last_skip_flag = (duration <= 0.0);
  if (state.last_skip_flag) {
    state.last_vel_cmd_pre = 0.0;
    state.last_vel_cmd_post = 0.0;
    state.last_boost = 0.0;
    last_apply_time_ = pm->time;
    PrintTopMdotWindowMessage(0.0);
    return false;
  }

  const bool state_changed =
      ApplyX3FrameVelocityTarget(mdot_window_target_vframe_, dt_boost);
  state.last_skip_flag = !state_changed;
  last_apply_time_ = pm->time;
  PrintTopMdotWindowMessage(state.last_boost);
  return state_changed;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::ApplyTracking()
//! \brief Apply the configured frame-tracking boost and report whether the
//! fluid state changed.

bool FrameTracker::ApplyTracking() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  if (nmb1 < 0) {
    return false;
  }

  const bool is_mhd = track_mhd_;
  auto &u0 = is_mhd ? pmy_pack->pmhd->u0 : pmy_pack->phydro->u0;
  auto &w0 = is_mhd ? pmy_pack->pmhd->w0 : pmy_pack->phydro->w0;
  const EOS_Data eos = is_mhd ? pmy_pack->pmhd->peos->eos_data :
                                pmy_pack->phydro->peos->eos_data;

  Real dt_boost = pm->dt*static_cast<Real>(apply_every_);
  if (last_apply_time_ >= 0.0) {
    dt_boost = pm->time - last_apply_time_;
  }
  if (dt_boost <= 0.0) {
    dt_boost = std::max(pm->dt, static_cast<Real>(1.0e-20));
  }

  if (mode_ == kFTTopMdot) {
    return ApplyTopMdotTracking(dt_boost);
  } else if (mode_ == kFTTopMdotAuto) {
    return ApplyTopMdotAutoTracking(dt_boost);
  } else if (mode_ == kFTTopMdotWindow) {
    return ApplyTopMdotWindowTracking(dt_boost);
  }

  SetActiveTargetRange();

  std::array<MomentSample, 3> samples;
  bool any_active = false;
  for (int axis = 0; axis < 3; ++axis) {
    any_active = any_active || axes_[axis].active;
  }
  if (!any_active) {
    return false;
  }
  const bool have_sample = SampleMoments(samples);
  Real global_weight = 0.0;
  if (have_sample) {
    for (int axis = 0; axis < 3; ++axis) {
      if (axes_[axis].active) {
        global_weight = samples[axis].global_weight;
        break;
      }
    }
  }

  const bool low_weight_floor = have_sample && (global_weight <= weight_floor_);
  const bool low_weight_observable =
      have_sample && (min_global_weight_ > 0.0) && (global_weight < min_global_weight_);

  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    AxisState &state = axes_[axis];
    state.last_global_weight = have_sample ? samples[axis].global_weight : 0.0;
    state.last_mean_x = samples[axis].mean_x;
    state.last_mean_v = samples[axis].mean_v;
    state.last_x_centroid = samples[axis].mean_x;
    state.last_x_midpoint = 0.5*(samples[axis].x_min + samples[axis].x_max);
    state.last_x_ctrl = samples[axis].mean_x;
  }

  if (!have_sample || low_weight_floor || low_weight_observable) {
    miss_streak_ += 1;
    recover_streak_ = 0;
    SetActiveTargetRange();
    for (int axis = 0; axis < 3; ++axis) {
      if (!axes_[axis].active) {
        continue;
      }
      AxisState &state = axes_[axis];
      state.last_skip_flag = true;
      state.last_slew_limited = false;
      state.last_vel_cmd_pre = 0.0;
      state.last_vel_cmd_post = 0.0;
      if (reacquire_leak_on_miss_ && tau_int_ > 0.0 && int_leak_tau_ > 0.0) {
        state.i_term *= std::exp(-dt_boost/int_leak_tau_);
      }
    }
    last_apply_time_ = pm->time;
    PrintSkipMessage(have_sample, low_weight_floor, global_weight);
    return false;
  }

  if (miss_streak_ > 0) {
    recover_streak_ += 1;
    if (recover_streak_ >= reacquire_recover_updates_) {
      miss_streak_ = std::max(0, miss_streak_ - 1);
      recover_streak_ = 0;
    }
  } else {
    recover_streak_ = 0;
  }
  SetActiveTargetRange();

  bool primed_filter = false;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    AxisState &state = axes_[axis];
    const Real midpoint = 0.5*(samples[axis].x_min + samples[axis].x_max);
    Real x_ctrl = samples[axis].mean_x;
    if (position_signal_ == kFTPosBandMidpoint) {
      x_ctrl = midpoint;
    } else if (position_signal_ == kFTPosBlend) {
      x_ctrl = (1.0 - position_blend_)*samples[axis].mean_x + position_blend_*midpoint;
    }

    state.last_global_weight = samples[axis].global_weight;
    state.last_mean_x = samples[axis].mean_x;
    state.last_mean_v = samples[axis].mean_v;
    state.last_x_centroid = samples[axis].mean_x;
    state.last_x_midpoint = midpoint;
    state.last_x_ctrl = x_ctrl;
    state.last_skip_flag = false;

    if (!state.filter_initialized) {
      state.last_filtered_x = x_ctrl;
      state.last_filtered_v = samples[axis].mean_v;
      state.filter_initialized = true;
      state.last_x_err = state.last_filtered_x - state.target_position;
      state.i_term = 0.0;
      state.last_vel_cmd_pre = 0.0;
      state.last_vel_cmd_post = 0.0;
      state.last_slew_limited = false;
      state.last_skip_flag = true;
      primed_filter = true;
    } else {
      Real alpha = 1.0 - std::exp(-dt_boost/tau_avg_);
      alpha = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), alpha));
      const Real previous_filtered_x = state.last_filtered_x;
      state.last_filtered_x += alpha*(x_ctrl - state.last_filtered_x);
      if (velocity_signal_ == kFTVelPositionRate) {
        state.last_filtered_v =
            (state.last_filtered_x - previous_filtered_x)/dt_boost;
      } else {
        state.last_filtered_v += alpha*(samples[axis].mean_v - state.last_filtered_v);
      }
    }
  }

  if (primed_filter) {
    last_apply_time_ = pm->time;
    PrintPrimeMessage(samples);
    return false;
  }

  std::array<Real, 3> v_p = {0.0, 0.0, 0.0};
  std::array<Real, 3> v_d = {0.0, 0.0, 0.0};
  std::array<Real, 3> v_i = {0.0, 0.0, 0.0};
  std::array<Real, 3> cmd_pre = {0.0, 0.0, 0.0};
  std::array<Real, 3> boost = {0.0, 0.0, 0.0};

  const bool use_position_feedback = (mode_ != kFTVelocity);
  const bool use_integral = use_position_feedback && (tau_int_ > 0.0) &&
                            (int_max_abs_ > 0.0);
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    AxisState &state = axes_[axis];
    const Real x_err = state.last_filtered_x - state.target_position;
    state.last_x_err = x_err;

    if (use_position_feedback) {
      const Real position_weight = 1.0 - std::exp(-dt_boost/tau_relax_);
      v_p[axis] = -position_weight*x_err/tau_relax_;
    }
    if (mode_ == kFTVelocity) {
      v_d[axis] = -state.last_filtered_v;
    } else if (mode_ == kFTPD) {
      const Real vel_weight = 1.0 - std::exp(-dt_boost/tau_vel_);
      v_d[axis] = -vel_weight*state.last_filtered_v;
    }

    if (!use_integral) {
      state.i_term = 0.0;
    } else {
      state.i_term *= std::exp(-dt_boost/int_leak_tau_);
      const Real vel_pre_no_accum = v_p[axis] + v_d[axis] + state.i_term;
      bool allow_integrate = true;
      if (int_unsat_only_ && max_abs_boost_ > 0.0 &&
          std::fabs(vel_pre_no_accum) >= max_abs_boost_) {
        allow_integrate = false;
      }
      if (allow_integrate) {
        state.i_term += -(dt_boost/tau_int_)*x_err;
      }
      if (std::fabs(state.i_term) > int_max_abs_) {
        state.i_term = std::copysign(int_max_abs_, state.i_term);
      }
    }
    if (use_integral) {
      const Real position_weight = 1.0 - std::exp(-dt_boost/tau_relax_);
      v_i[axis] = position_weight*state.i_term;
    } else {
      v_i[axis] = 0.0;
    }
    cmd_pre[axis] = v_p[axis] + v_d[axis] + v_i[axis];

    Real cmd = cmd_pre[axis];
    if (boundary_guard_enable_ && boundary_guard_cells_ > 0) {
      const Real domain_dx =
          (AxisGlobalMax(axis, pm) - AxisGlobalMin(axis, pm)) /
          std::max(static_cast<Real>(1.0),
                   static_cast<Real>(AxisCellCount(axis, pm->mesh_indcs)));
      const Real guard_thickness = static_cast<Real>(boundary_guard_cells_)*domain_dx;
      if (guard_thickness > 0.0) {
        const Real dist_to_lower = state.last_filtered_x - AxisGlobalMin(axis, pm);
        const Real dist_to_upper = AxisGlobalMax(axis, pm) - state.last_filtered_x;
        const Real dist_to_edge = std::max(static_cast<Real>(0.0),
                                           std::min(dist_to_lower, dist_to_upper));
        if (dist_to_edge < guard_thickness) {
          Real frac = dist_to_edge/guard_thickness;
          frac = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), frac));
          Real guard_scale = boundary_guard_min_scale_ +
                             (1.0 - boundary_guard_min_scale_)*frac;
          guard_scale = std::max(boundary_guard_min_scale_,
                                 std::min(static_cast<Real>(1.0), guard_scale));
          cmd *= guard_scale;
        }
      }
    }

    if (max_abs_boost_ > 0.0 && std::fabs(cmd) > max_abs_boost_) {
      cmd = std::copysign(max_abs_boost_, cmd);
    }

    if (state.use_upper_limit || state.use_lower_limit) {
      const Real x_pred = state.last_filtered_x + cmd*dt_boost;
      if (state.use_upper_limit && x_pred > state.upper_limit) {
        cmd = (state.upper_limit - state.last_filtered_x)/dt_boost;
      }
      if (state.use_lower_limit && x_pred < state.lower_limit) {
        cmd = (state.lower_limit - state.last_filtered_x)/dt_boost;
      }
    }

    state.last_slew_limited = false;
    Real dv_allowed = 0.0;
    if (boost_change_mode_ == kFTBoostPerApply) {
      dv_allowed = max_boost_change_;
    } else if (boost_change_mode_ == kFTBoostPerTime && max_boost_change_rate_ > 0.0) {
      dv_allowed = max_boost_change_rate_*dt_boost;
    }
    if (dv_allowed > 0.0) {
      const Real delta_boost = cmd - state.last_boost;
      if (std::fabs(delta_boost) > dv_allowed) {
        cmd = state.last_boost + std::copysign(dv_allowed, delta_boost);
        state.last_slew_limited = true;
      }
    }
    if (max_abs_boost_ > 0.0 && std::fabs(cmd) > max_abs_boost_) {
      cmd = std::copysign(max_abs_boost_, cmd);
    }

    boost[axis] = cmd;
    state.last_vel_cmd_pre = cmd_pre[axis];
    state.last_vel_cmd_post = boost[axis];
    state.last_boost = boost[axis];
    state.cumulative_boost += boost[axis];
    frame_velocity_[axis] = -state.cumulative_boost;
  }
  last_apply_time_ = pm->time;

  const Real dv1 = boost[0];
  const Real dv2 = boost[1];
  const Real dv3 = boost[2];
  const Real dv2_total = SQR(dv1) + SQR(dv2) + SQR(dv3);
  const bool use_energy = eos.is_ideal;
  const bool state_changed = (dv1 != 0.0) || (dv2 != 0.0) || (dv3 != 0.0);
  if (state_changed) {
    par_for("frame_tracking_boost", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real rho = w0(m, IDN, k, j, i);
      if (rho <= 0.0) {
        return;
      }
      const Real mom_old1 = u0(m, IM1, k, j, i);
      const Real mom_old2 = u0(m, IM2, k, j, i);
      const Real mom_old3 = u0(m, IM3, k, j, i);
      w0(m, IVX, k, j, i) += dv1;
      w0(m, IVY, k, j, i) += dv2;
      w0(m, IVZ, k, j, i) += dv3;
      u0(m, IM1, k, j, i) += rho*dv1;
      u0(m, IM2, k, j, i) += rho*dv2;
      u0(m, IM3, k, j, i) += rho*dv3;
      if (use_energy) {
        u0(m, IEN, k, j, i) += mom_old1*dv1 + mom_old2*dv2 + mom_old3*dv3 +
                                0.5*rho*dv2_total;
      }
    });
  }

  PrintApplyMessage(samples, v_p, v_d, v_i, cmd_pre, boost);
  return state_changed;
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::AdvanceFrameDisplacement()
//! \brief Advance the lab-frame origin of the tracked grid frame.

bool FrameTracker::AdvanceFrameDisplacement(const Real dt) {
  if (dt <= 0.0) {
    return false;
  }
  bool changed = false;
  for (int axis = 0; axis < 3; ++axis) {
    const Real dx = frame_velocity_[axis]*dt;
    frame_displacement_[axis] += dx;
    changed = changed || (dx != 0.0);
  }
  return changed;
}

//! \fn void FrameTracker::PrintSkipMessage()
//! \brief Print diagnostic information for skipped updates.

void FrameTracker::PrintSkipMessage(const bool have_sample, const bool low_weight_floor,
                                    const Real global_weight) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  if (!have_sample) {
    std::cout << "FrameTracker skipped: no tracked material in selected target range"
              << std::endl;
  } else if (low_weight_floor) {
    std::cout << "FrameTracker skipped: global_weight <= weight_floor (" << global_weight
              << " <= " << weight_floor_ << ")" << std::endl;
  } else {
    std::cout << "FrameTracker skipped: global_weight < min_global_weight ("
              << global_weight << " < " << min_global_weight_ << ")" << std::endl;
  }
  std::cout << " ft_miss_streak=" << miss_streak_ << std::endl;
  std::cout << " target_min_active=" << target_min_active_ << std::endl;
  std::cout << " target_max_active=" << target_max_active_ << std::endl;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    std::cout << " " << AxisName(axis) << "_frame_velocity="
              << frame_velocity_[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_frame_displacement="
              << frame_displacement_[axis] << std::endl;
  }
}

void FrameTracker::PrintPrimeMessage(const std::array<MomentSample, 3> &samples) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  std::cout << "FrameTracker primed filter state; skipping first actuation" << std::endl;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    const AxisState &state = axes_[axis];
    std::cout << " " << AxisName(axis) << "_global_weight="
              << samples[axis].global_weight << std::endl;
    std::cout << " " << AxisName(axis) << "_centroid="
              << samples[axis].mean_x << std::endl;
    std::cout << " " << AxisName(axis) << "_mean_velocity="
              << samples[axis].mean_v << std::endl;
    std::cout << " " << AxisName(axis) << "_midpoint="
              << state.last_x_midpoint << std::endl;
    std::cout << " " << AxisName(axis) << "_control_position="
              << state.last_x_ctrl << std::endl;
    std::cout << " " << AxisName(axis) << "_filtered_position="
              << state.last_filtered_x << std::endl;
    std::cout << " " << AxisName(axis) << "_filtered_velocity="
              << state.last_filtered_v << std::endl;
    std::cout << " " << AxisName(axis) << "_position_error="
              << state.last_x_err << std::endl;
    std::cout << " " << AxisName(axis) << "_frame_velocity="
              << frame_velocity_[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_frame_displacement="
              << frame_displacement_[axis] << std::endl;
  }
  std::cout << " target_min_active=" << target_min_active_ << std::endl;
  std::cout << " target_max_active=" << target_max_active_ << std::endl;
}

void FrameTracker::PrintApplyMessage(const std::array<MomentSample, 3> &samples,
                                     const std::array<Real, 3> &v_p,
                                     const std::array<Real, 3> &v_d,
                                     const std::array<Real, 3> &v_i,
                                     const std::array<Real, 3> &cmd_pre,
                                     const std::array<Real, 3> &boost) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  std::cout << "FrameTracker applied boost" << std::endl;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    const AxisState &state = axes_[axis];
    std::cout << " " << AxisName(axis) << "_global_weight="
              << samples[axis].global_weight << std::endl;
    std::cout << " " << AxisName(axis) << "_centroid="
              << samples[axis].mean_x << std::endl;
    std::cout << " " << AxisName(axis) << "_mean_velocity="
              << samples[axis].mean_v << std::endl;
    std::cout << " " << AxisName(axis) << "_midpoint="
              << state.last_x_midpoint << std::endl;
    std::cout << " " << AxisName(axis) << "_control_position="
              << state.last_x_ctrl << std::endl;
    std::cout << " " << AxisName(axis) << "_filtered_position="
              << state.last_filtered_x << std::endl;
    std::cout << " " << AxisName(axis) << "_filtered_velocity="
              << state.last_filtered_v << std::endl;
    std::cout << " " << AxisName(axis) << "_position_error="
              << state.last_x_err << std::endl;
    std::cout << " " << AxisName(axis) << "_v_p=" << v_p[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_v_d=" << v_d[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_v_i=" << v_i[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_cmd_pre=" << cmd_pre[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_boost=" << boost[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_frame_velocity="
              << frame_velocity_[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_frame_displacement="
              << frame_displacement_[axis] << std::endl;
    std::cout << " " << AxisName(axis) << "_slew_limited="
              << state.last_slew_limited << std::endl;
  }
  std::cout << " target_min_active=" << target_min_active_ << std::endl;
  std::cout << " target_max_active=" << target_max_active_ << std::endl;
  std::cout << " ft_miss_streak=" << miss_streak_ << std::endl;
}

void FrameTracker::PrintTopMdotMessage(const Real top_mdot, const Real avg_mdot,
                                       const Real target_frame_velocity,
                                       const Real ramp_factor,
                                       const Real boost) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  const AxisState &state = axes_[2];
  std::cout << "FrameTracker top_mdot update" << std::endl;
  std::cout << " sign_contract=Mdot_top_integral_rho_v3_dA_upper_x3;"
            << "positive_plus_x3_outward;"
            << "v_frame_target=Mdot_top/(rho_cold*A);"
            << "fluid_boost=-Delta_v_frame" << std::endl;
  std::cout << " x3_mdot_top=" << top_mdot << std::endl;
  std::cout << " x3_mdot_average=" << avg_mdot << std::endl;
  std::cout << " x3_mdot_average_time=" << mdot_weight_time_ << std::endl;
  std::cout << " x3_mdot_area=" << mdot_area_ << std::endl;
  std::cout << " x3_mdot_density=" << mdot_density_ << std::endl;
  std::cout << " x3_mdot_ramp_factor=" << ramp_factor << std::endl;
  std::cout << " x3_target_frame_velocity=" << target_frame_velocity << std::endl;
  std::cout << " x3_boost=" << boost << std::endl;
  std::cout << " x3_frame_velocity=" << frame_velocity_[2] << std::endl;
  std::cout << " x3_frame_displacement=" << frame_displacement_[2] << std::endl;
  std::cout << " x3_velocity_error=" << state.last_x_err << std::endl;
  std::cout << " x3_slew_limited=" << state.last_slew_limited << std::endl;
  std::cout << " x3_skip=" << state.last_skip_flag << std::endl;
}

void FrameTracker::PrintTopMdotAutoMessage(const Real boost) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  const AxisState &state = axes_[2];
  std::cout << "FrameTracker top_mdot_auto update" << std::endl;
  std::cout << " sign_contract=Mdot_top_integral_rho_v3_dA_upper_x3;"
            << "positive_plus_x3_outward;"
            << "v_frame_target=Mdot_top/(rho_cold*A);"
            << "fluid_boost=-Delta_v_frame" << std::endl;
  std::cout << " x3_mdot_top=" << mdot_last_top_ << std::endl;
  std::cout << " x3_mdot_fast=" << mdot_auto_fast_ << std::endl;
  std::cout << " x3_mdot_slow=" << mdot_auto_slow_ << std::endl;
  std::cout << " x3_mdot_rel=" << mdot_auto_rel_ << std::endl;
  std::cout << " x3_mdot_hold_elapsed=" << mdot_auto_hold_elapsed_ << std::endl;
  std::cout << " x3_mdot_hold_required=" << mdot_auto_hold_time_ << std::endl;
  std::cout << " x3_mdot_settled=" << mdot_auto_locked_ << std::endl;
  std::cout << " x3_mdot_settle_time=" << mdot_auto_settle_time_ << std::endl;
  std::cout << " x3_mdot_locked=" << mdot_auto_locked_mdot_ << std::endl;
  std::cout << " x3_mdot_area=" << mdot_area_ << std::endl;
  std::cout << " x3_mdot_density=" << mdot_density_ << std::endl;
  std::cout << " x3_locked_frame_velocity=" << mdot_auto_locked_vframe_
            << std::endl;
  std::cout << " x3_target_frame_velocity=" << mdot_auto_target_vframe_
            << std::endl;
  std::cout << " x3_ramp_duration=" << mdot_auto_ramp_duration_ << std::endl;
  std::cout << " x3_ramp_factor=" << mdot_auto_ramp_factor_ << std::endl;
  std::cout << " x3_frame_delta_applied=" << -boost << std::endl;
  std::cout << " x3_fluid_boost=" << boost << std::endl;
  std::cout << " x3_frame_velocity=" << frame_velocity_[2] << std::endl;
  std::cout << " x3_frame_displacement=" << frame_displacement_[2] << std::endl;
  std::cout << " x3_velocity_error=" << state.last_x_err << std::endl;
  std::cout << " x3_slew_limited=" << state.last_slew_limited << std::endl;
  std::cout << " x3_skip=" << state.last_skip_flag << std::endl;
}

void FrameTracker::PrintTopMdotWindowMessage(const Real boost) const {
  if (global_variable::my_rank != 0 || diagnostic_every_ <= 0 ||
      pmy_pack->pmesh->ncycle % diagnostic_every_ != 0) {
    return;
  }
  const AxisState &state = axes_[2];
  std::cout << "FrameTracker top_mdot_window update" << std::endl;
  std::cout << " sign_contract=Mdot_top_integral_rho_v3_dA_upper_x3;"
            << "v3_frame=" << NormalizeToken(mdot_flux_frame_name_) << ";"
            << "positive_plus_x3_outward;"
            << "v_frame_target=sign*window_avg_Mdot/(rho_frame*A);"
            << "fluid_boost=-Delta_v_frame" << std::endl;
  std::cout << " x3_mdot_top=" << mdot_last_top_ << std::endl;
  std::cout << " x3_mdot_window_average=" << mdot_last_average_ << std::endl;
  std::cout << " x3_mdot_window_start_time=" << mdot_window_start_time_ << std::endl;
  std::cout << " x3_mdot_window_time=" << mdot_window_time_ << std::endl;
  std::cout << " x3_mdot_window_low_time=" << mdot_window_low_time_ << std::endl;
  std::cout << " x3_mdot_window_weight_time=" << mdot_weight_time_ << std::endl;
  std::cout << " x3_mdot_area=" << mdot_area_ << std::endl;
  std::cout << " x3_mdot_density_mode=" << NormalizeToken(mdot_density_mode_name_)
            << std::endl;
  std::cout << " x3_mdot_density=" << mdot_density_ << std::endl;
  std::cout << " x3_mdot_velocity_sign=" << mdot_velocity_sign_ << std::endl;
  std::cout << " x3_mdot_velocity_factor=" << mdot_velocity_factor_ << std::endl;
  std::cout << " x3_mdot_flux_frame=" << NormalizeToken(mdot_flux_frame_name_)
            << std::endl;
  std::cout << " x3_target_frame_velocity=" << mdot_window_target_vframe_
            << std::endl;
  std::cout << " x3_frame_delta_applied=" << -boost << std::endl;
  std::cout << " x3_fluid_boost=" << boost << std::endl;
  std::cout << " x3_frame_velocity=" << frame_velocity_[2] << std::endl;
  std::cout << " x3_frame_displacement=" << frame_displacement_[2] << std::endl;
  std::cout << " x3_velocity_error=" << state.last_x_err << std::endl;
  std::cout << " x3_slew_limited=" << state.last_slew_limited << std::endl;
  std::cout << " x3_skip=" << state.last_skip_flag << std::endl;
}
