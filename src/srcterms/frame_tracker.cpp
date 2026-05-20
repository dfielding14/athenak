//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file frame_tracker.cpp
//! \brief Implements shared frame tracking for hydro/MHD problem generators.

#include "frame_tracker.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
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
  kFTPD = 2
};

enum FrameTrackingPositionSignal {
  kFTPosCentroid = 0,
  kFTPosBandMidpoint = 1,
  kFTPosBlend = 2
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
  kFTWeightTarget = 2
};

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

Real GetOrAddRealAny(ParameterInput *pin, const std::string &block,
                     const std::vector<std::string> &names, const Real default_value) {
  Real value = default_value;
  if (ReadRealAny(pin, block, names, value)) {
    return value;
  }
  return pin->GetOrAddReal(block, names.front(), default_value);
}

int GetOrAddIntegerAny(ParameterInput *pin, const std::string &block,
                       const std::vector<std::string> &names, const int default_value) {
  int value = default_value;
  if (ReadIntegerAny(pin, block, names, value)) {
    return value;
  }
  return pin->GetOrAddInteger(block, names.front(), default_value);
}

bool GetOrAddBooleanAny(ParameterInput *pin, const std::string &block,
                        const std::vector<std::string> &names, const bool default_value) {
  bool value = default_value;
  if (ReadBooleanAny(pin, block, names, value)) {
    return value;
  }
  return pin->GetOrAddBoolean(block, names.front(), default_value);
}

std::string GetOrAddStringAny(ParameterInput *pin, const std::string &block,
                              const std::vector<std::string> &names,
                              const std::string &default_value) {
  std::string value = default_value;
  if (ReadStringAny(pin, block, names, value)) {
    return value;
  }
  return pin->GetOrAddString(block, names.front(), default_value);
}

int ParseFrameTrackingMode(const std::string &mode_raw) {
  const std::string mode = NormalizeToken(mode_raw);
  if (mode == "velocity") return kFTVelocity;
  if (mode == "position") return kFTPosition;
  if (mode == "pd" || mode == "position_velocity") return kFTPD;
  FatalFrameTrackingInput("Invalid <frame_tracking>/mode '" + mode_raw +
                          "'. Expected one of: velocity, position, pd.");
  return kFTPD;
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
  FatalFrameTrackingInput("Invalid <frame_tracking>/weight '" + mode_raw +
                          "'. Expected one of: mass, volume, target.");
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
int CoordIndex(const int axis, const int i, const int j, const int k,
               const RegionIndcs &indcs) {
  if (axis == 0) return i - indcs.is;
  if (axis == 1) return j - indcs.js;
  return k - indcs.ks;
}

KOKKOS_INLINE_FUNCTION
int AxisCellCount(const int axis, const RegionIndcs &indcs) {
  if (axis == 0) return indcs.nx1;
  if (axis == 1) return indcs.nx2;
  return indcs.nx3;
}

KOKKOS_INLINE_FUNCTION
Real AxisMin(const int axis, const RegionSize &size) {
  if (axis == 0) return size.x1min;
  if (axis == 1) return size.x2min;
  return size.x3min;
}

KOKKOS_INLINE_FUNCTION
Real AxisMax(const int axis, const RegionSize &size) {
  if (axis == 0) return size.x1max;
  if (axis == 1) return size.x2max;
  return size.x3max;
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
  enabled_ = GetOrAddBooleanAny(pin, block_name_, {"enabled", "enable"}, true);
  apply_every_ = GetOrAddIntegerAny(pin, block_name_,
                                    {"apply_every", "ft_apply_every"}, 1);
  start_time_ = GetOrAddRealAny(pin, block_name_,
                                {"start_time", "t_start", "t_frame_tracking_start"}, 0.0);
  diagnostic_every_ = GetOrAddIntegerAny(pin, block_name_,
                                         {"diagnostic_every", "ndiag"}, -1);
  verbose_ = GetOrAddBooleanAny(pin, block_name_, {"verbose"}, false);

  std::string target = GetOrAddStringAny(pin, block_name_,
                                         {"target", "target_variable"}, "density");
  if (NormalizeToken(target) == "derived") {
    target = GetOrAddStringAny(pin, block_name_, {"derived_variable"}, "temperature");
  }
  scalar_index_ = GetOrAddIntegerAny(pin, block_name_, {"scalar_index"}, 0);
  target_kind_ = ParseTargetKind(target, scalar_index_);
  weight_mode_ = ParseWeightMode(GetOrAddStringAny(pin, block_name_, {"weight"}, "mass"));

  const Real huge = std::numeric_limits<Real>::max()/static_cast<Real>(4.0);
  target_min_ = -huge;
  target_max_ = huge;
  if (target_kind_ == kFTTargetTemperature) {
    (void) ReadRealAny(pin, block_name_, {"target_min", "min", "ft_temp_lo"}, target_min_);
    (void) ReadRealAny(pin, block_name_, {"target_max", "max", "ft_temp_hi"}, target_max_);
  } else {
    (void) ReadRealAny(pin, block_name_, {"target_min", "min"}, target_min_);
    (void) ReadRealAny(pin, block_name_, {"target_max", "max"}, target_max_);
  }

  if (!ReadRealAny(pin, block_name_, {"target_center", "center"}, target_center_)) {
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
  (void) ReadStringAny(pin, block_name_, {"target_scale", "reacquire_scale"}, target_scale);
  target_scale = NormalizeToken(target_scale);
  if (target_scale == "log" || target_scale == "logarithmic") {
    use_log_reacquire_ = true;
  } else if (target_scale == "linear") {
    use_log_reacquire_ = false;
  } else if (target_scale == "auto") {
    use_log_reacquire_ = (target_min_ > 0.0 && target_max_ > 0.0 && target_center_ > 0.0);
  } else {
    FatalFrameTrackingInput("Invalid <frame_tracking>/target_scale '" + target_scale +
                            "'. Expected auto, linear, or log.");
  }

  mode_ = ParseFrameTrackingMode(GetOrAddStringAny(pin, block_name_,
                                                   {"mode", "ft_mode"}, "pd"));
  position_signal_ = ParseFrameTrackingPositionSignal(
      GetOrAddStringAny(pin, block_name_,
                        {"position_signal", "ft_position_signal"}, "blend"));
  boost_change_mode_ = ParseFrameTrackingBoostChangeMode(
      GetOrAddStringAny(pin, block_name_,
                        {"max_boost_change_mode", "ft_max_boost_change_mode"},
                        "per_apply"));
  position_blend_ = GetOrAddRealAny(pin, block_name_,
                                    {"position_blend", "ft_position_blend"}, 0.7);
  tau_avg_ = GetOrAddRealAny(pin, block_name_, {"tau_avg", "ft_tau_avg"}, 1.0);
  tau_relax_ = GetOrAddRealAny(pin, block_name_, {"tau_relax", "ft_tau_relax"}, 1.0);
  tau_vel_ = GetOrAddRealAny(pin, block_name_, {"tau_vel", "ft_tau_vel"}, 1.0);
  tau_int_ = GetOrAddRealAny(pin, block_name_, {"tau_int", "ft_tau_int"}, 0.0);
  int_max_abs_ = GetOrAddRealAny(pin, block_name_,
                                 {"int_max_abs", "ft_int_max_abs"}, 0.0);
  int_leak_tau_ = GetOrAddRealAny(pin, block_name_,
                                  {"int_leak_tau", "ft_int_leak_tau"}, 1.0);
  int_unsat_only_ = GetOrAddBooleanAny(pin, block_name_,
                                       {"int_unsat_only", "ft_int_unsat_only"}, true);
  weight_floor_ = GetOrAddRealAny(pin, block_name_,
                                  {"weight_floor", "ft_weight_floor"}, 0.0);
  min_global_weight_ = GetOrAddRealAny(pin, block_name_,
                                       {"min_global_weight", "ft_min_global_weight"}, 0.0);
  max_abs_boost_ = GetOrAddRealAny(pin, block_name_,
                                   {"max_abs_boost", "ft_max_abs_boost"}, 0.0);
  max_boost_change_ = GetOrAddRealAny(pin, block_name_,
                                      {"max_boost_change", "ft_max_boost_change"}, 0.0);
  max_boost_change_rate_ = GetOrAddRealAny(
      pin, block_name_, {"max_boost_change_rate", "ft_max_boost_change_rate"}, -1.0);
  reacquire_expand_factor_ = GetOrAddRealAny(
      pin, block_name_, {"reacquire_expand_factor", "ft_reacquire_expand_factor"}, 1.0);
  reacquire_max_expand_ = GetOrAddRealAny(
      pin, block_name_, {"reacquire_max_expand", "ft_reacquire_max_expand"}, 1.0);
  reacquire_recover_updates_ = GetOrAddIntegerAny(
      pin, block_name_, {"reacquire_recover_updates", "ft_reacquire_recover_updates"}, 1);
  reacquire_leak_on_miss_ = GetOrAddBooleanAny(
      pin, block_name_, {"reacquire_leak_on_miss", "ft_reacquire_leak_on_miss"}, true);
  boundary_guard_enable_ = GetOrAddBooleanAny(
      pin, block_name_, {"boundary_guard_enable", "ft_boundary_guard_enable"}, false);
  boundary_guard_cells_ = GetOrAddIntegerAny(
      pin, block_name_, {"boundary_guard_cells", "ft_boundary_guard_cells"}, 8);
  boundary_guard_min_scale_ = GetOrAddRealAny(
      pin, block_name_, {"boundary_guard_min_scale", "ft_boundary_guard_min_scale"}, 0.15);

  ParseAxes(pin);
  ParseAxisControls(pin);
  ValidateConfiguration();
  RestoreFrameState(pin);

  if (global_variable::my_rank == 0 && verbose_) {
    std::cout << "FrameTracker enabled from <" << block_name_ << "> with target='"
              << target << "', target_min=" << std::setprecision(16) << target_min_
              << ", target_max=" << target_max_ << ", axes=";
    bool first = true;
    for (int axis = 0; axis < 3; ++axis) {
      if (axes_[axis].active) {
        std::cout << (first ? "" : ",") << AxisName(axis);
        first = false;
      }
    }
    std::cout << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::RestoreFrameState()
//! \brief Restore moving-frame state stored in restart input headers.

void FrameTracker::RestoreFrameState(ParameterInput *pin) {
  if (pin == nullptr) {
    return;
  }
  for (int axis = 0; axis < 3; ++axis) {
    const std::string name = AxisName(axis);
    (void) ReadRealAny(pin, block_name_,
                       {"frame_velocity_" + name, name + "_frame_velocity"},
                       frame_velocity_[axis]);
    (void) ReadRealAny(pin, block_name_,
                       {"frame_displacement_" + name, name + "_frame_displacement"},
                       frame_displacement_[axis]);
    axes_[axis].cumulative_boost = -frame_velocity_[axis];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::StoreStateInParameterInput()
//! \brief Store moving-frame state in the input header for restart files.

void FrameTracker::StoreStateInParameterInput(ParameterInput *pin) const {
  if (pin == nullptr) {
    return;
  }
  for (int axis = 0; axis < 3; ++axis) {
    const std::string name = AxisName(axis);
    pin->SetReal(block_name_, "frame_velocity_" + name, frame_velocity_[axis]);
    pin->SetReal(block_name_, "frame_displacement_" + name, frame_displacement_[axis]);
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
  if (!ReadStringAny(pin, block_name_, {"axes", "axis", "directions"}, axes_raw)) {
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

    std::vector<std::string> target_keys = {
        name + "_target", name + "_target_position", "target_" + name};
    if (axis == 2) {
      target_keys.push_back("z_target");
      target_keys.push_back("ft_z_target");
    }
    (void) ReadRealAny(pin, block_name_, target_keys, axes_[axis].target_position);

    std::vector<std::string> lower_keys = {
        name + "_lower_limit", name + "_min_limit", name + "_low_limit"};
    std::vector<std::string> upper_keys = {
        name + "_upper_limit", name + "_max_limit", name + "_upp_limit"};
    if (axis == 2) {
      lower_keys.push_back("z_lower_limit");
      lower_keys.push_back("ft_lowzlim");
      upper_keys.push_back("z_upper_limit");
      upper_keys.push_back("ft_uppzlim");
    }

    bool use_lower_flag = false;
    bool use_upper_flag = false;
    if (axis == 2) {
      (void) ReadBooleanAny(pin, block_name_, {"ft_use_lowzlim"}, use_lower_flag);
      (void) ReadBooleanAny(pin, block_name_, {"ft_use_uppzlim"}, use_upper_flag);
    }
    if (ReadRealAny(pin, block_name_, lower_keys, axes_[axis].lower_limit)) {
      axes_[axis].use_lower_limit = true;
    } else {
      axes_[axis].use_lower_limit = use_lower_flag;
    }
    if (ReadRealAny(pin, block_name_, upper_keys, axes_[axis].upper_limit)) {
      axes_[axis].use_upper_limit = true;
    } else {
      axes_[axis].use_upper_limit = use_upper_flag;
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

  if (apply_every_ < 1) {
    FatalFrameTrackingInput("Require <frame_tracking>/apply_every >= 1.");
  }
  if (tau_avg_ <= 0.0 || tau_relax_ <= 0.0 || tau_vel_ <= 0.0) {
    FatalFrameTrackingInput("Frame-tracking tau_avg, tau_relax, and tau_vel must be positive.");
  }
  if (target_max_ < target_min_) {
    FatalFrameTrackingInput("Require <frame_tracking>/target_max >= target_min.");
  }
  if (weight_floor_ < 0.0 || min_global_weight_ < 0.0 || max_abs_boost_ < 0.0 ||
      max_boost_change_ < 0.0 || int_max_abs_ < 0.0) {
    FatalFrameTrackingInput("Frame-tracking floors and boost/integral limits must be non-negative.");
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

  const bool is_mhd = (pmy_pack->pmhd != nullptr);
  const int nscalars = is_mhd ? pmy_pack->pmhd->nscalars : pmy_pack->phydro->nscalars;
  const EOS_Data &eos = is_mhd ? pmy_pack->pmhd->peos->eos_data :
                                 pmy_pack->phydro->peos->eos_data;
  if (target_kind_ == kFTTargetScalar && (scalar_index_ < 0 || scalar_index_ >= nscalars)) {
    FatalFrameTrackingInput("<frame_tracking> scalar target requested scalar_index=" +
                            std::to_string(scalar_index_) + ", but the active fluid has " +
                            std::to_string(nscalars) + " passive scalar(s).");
  }
  if (target_kind_ == kFTTargetEntropy && !eos.is_ideal) {
    FatalFrameTrackingInput("<frame_tracking>/target=entropy requires an ideal-gas EOS.");
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::IncludeFrameTrackingTask()
//! \brief Add frame tracking to the after-time-integrator task list.

void FrameTracker::IncludeFrameTrackingTask(std::shared_ptr<TaskList> tl, TaskID start) {
  if (tl == nullptr) {
    return;
  }
  TaskID dep = start;
  (void) tl->AddTask(&FrameTracker::Apply, this, dep);
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus FrameTracker::Apply()
//! \brief Task-list wrapper for scheduled frame tracking.

TaskStatus FrameTracker::Apply(Driver *pdrive, int stage) {
  (void) pdrive;
  (void) stage;
  if (!enabled_ || pmy_pack == nullptr || pmy_pack->pmesh == nullptr) {
    return TaskStatus::complete;
  }
  Mesh *pm = pmy_pack->pmesh;
  AdvanceFrameDisplacement(pm->dt);
  if (pm->time <= start_time_ || pm->ncycle % apply_every_ != 0) {
    return TaskStatus::complete;
  }
  ApplyTracking();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
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
    target_min_active_ = std::pow(static_cast<Real>(10.0), log_center - dlog_lo*expand);
    target_max_active_ = std::pow(static_cast<Real>(10.0), log_center + dlog_hi*expand);
  } else {
    target_min_active_ = target_center_ - (target_center_ - target_min_)*expand;
    target_max_active_ = target_center_ + (target_max_ - target_center_)*expand;
  }
  if (target_max_active_ < target_min_active_) {
    std::swap(target_min_active_, target_max_active_);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool FrameTracker::SampleAxisMoments()
//! \brief Compute global weighted moments for the selected target material along one axis.

bool FrameTracker::SampleAxisMoments(const int axis, MomentSample &sample) const {
  sample = MomentSample();
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
  const bool is_mhd = (pmbp->pmhd != nullptr);
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
  Real rank_weight = 0.0;
  Real rank_weighted_x = 0.0;
  Real rank_weighted_v = 0.0;
  Real rank_xmin = FLT_MAX;
  Real rank_xmax = -FLT_MAX;
  int rank_count = 0;

  Kokkos::parallel_reduce("frame_tracking_moments",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &weight_, Real &weighted_x_, Real &weighted_v_,
                Real &xmin_, Real &xmax_, int &count_) {
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

    const Real value = FrameTrackingTargetValue(target_kind, scalar_var, w0, eos, m, k, j, i);
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
    }
    if (weight <= 0.0) {
      return;
    }

    const Real x = CellCenterX(CoordIndex(axis, i, j, k, indcs),
                               AxisCellCount(axis, indcs),
                               AxisMin(axis, block_size),
                               AxisMax(axis, block_size));
    const int velocity_index = IVX + axis;
    weight_ += weight;
    weighted_x_ += weight*x;
    weighted_v_ += weight*w0(m, velocity_index, k, j, i);
    xmin_ = fmin(xmin_, x);
    xmax_ = fmax(xmax_, x);
    count_ += 1;
  }, Kokkos::Sum<Real>(rank_weight),
     Kokkos::Sum<Real>(rank_weighted_x),
     Kokkos::Sum<Real>(rank_weighted_v),
     Kokkos::Min<Real>(rank_xmin),
     Kokkos::Max<Real>(rank_xmax),
     Kokkos::Sum<int>(rank_count));

  Real global_sums[3] = {rank_weight, rank_weighted_x, rank_weighted_v};
  Real global_xmin = rank_xmin;
  Real global_xmax = rank_xmax;
  int global_count = rank_count;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, global_sums, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_xmin, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_xmax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  sample.global_weight = global_sums[0];
  if (sample.global_weight <= 0.0 || global_count <= 0) {
    return false;
  }
  sample.mean_x = global_sums[1]/sample.global_weight;
  sample.mean_v = global_sums[2]/sample.global_weight;
  sample.x_min = global_xmin;
  sample.x_max = global_xmax;
  sample.valid = true;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::ApplyTracking()
//! \brief Apply the configured frame-tracking boost.

void FrameTracker::ApplyTracking() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  if (nmb1 < 0) {
    return;
  }

  const bool is_mhd = (pmy_pack->pmhd != nullptr);
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

  SetActiveTargetRange();

  std::array<MomentSample, 3> samples;
  bool have_sample = true;
  bool any_active = false;
  Real global_weight = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    if (!axes_[axis].active) {
      continue;
    }
    any_active = true;
    const bool axis_valid = SampleAxisMoments(axis, samples[axis]);
    have_sample = have_sample && axis_valid;
    if (axis_valid) {
      global_weight = samples[axis].global_weight;
    }
  }
  if (!any_active) {
    return;
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
    return;
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
      state.last_filtered_x += alpha*(x_ctrl - state.last_filtered_x);
      state.last_filtered_v += alpha*(samples[axis].mean_v - state.last_filtered_v);
    }
  }

  if (primed_filter) {
    last_apply_time_ = pm->time;
    PrintPrimeMessage(samples);
    return;
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
      v_p[axis] = -x_err/tau_relax_;
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
    v_i[axis] = use_integral ? state.i_term : 0.0;
    cmd_pre[axis] = v_p[axis] + v_d[axis] + v_i[axis];

    Real cmd = cmd_pre[axis];
    if (boundary_guard_enable_ && boundary_guard_cells_ > 0) {
      const Real domain_dx = (AxisGlobalMax(axis, pm) - AxisGlobalMin(axis, pm)) /
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

  PrintApplyMessage(samples, v_p, v_d, v_i, cmd_pre, boost);
}

//----------------------------------------------------------------------------------------
//! \fn void FrameTracker::AdvanceFrameDisplacement()
//! \brief Advance the lab-frame origin of the tracked grid frame.

void FrameTracker::AdvanceFrameDisplacement(const Real dt) {
  if (dt <= 0.0) {
    return;
  }
  for (int axis = 0; axis < 3; ++axis) {
    frame_displacement_[axis] += frame_velocity_[axis]*dt;
  }
}

//----------------------------------------------------------------------------------------
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
    std::cout << "FrameTracker skipped: global_weight < min_global_weight (" << global_weight
              << " < " << min_global_weight_ << ")" << std::endl;
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
