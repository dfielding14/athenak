#ifndef SRCTERMS_FRAME_TRACKER_HPP_
#define SRCTERMS_FRAME_TRACKER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file frame_tracker.hpp
//! \brief Defines a shared frame-tracking source-term helper.

#include <array>
#include <deque>
#include <memory>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"

class Driver;

//----------------------------------------------------------------------------------------
//! \class FrameTracker
//! \brief Keeps selected material near target positions with Galilean boosts.

class FrameTracker {
 public:
  FrameTracker(MeshBlockPack *pp, ParameterInput *pin,
               const std::string &block_name = "frame_tracking");
  ~FrameTracker() = default;

  TaskID IncludeFrameTrackingTask(std::shared_ptr<TaskList> tl, TaskID start);
  TaskStatus Apply(Driver *pdrive, int stage);
  bool ApplyTracking();

  void UpdateMeshBlockPack(MeshBlockPack *new_pp) { pmy_pack = new_pp; }
  Real FrameVelocity(const int axis) const {
    return (axis >= 0 && axis < 3) ? frame_velocity_[axis] : 0.0;
  }
  Real FrameDisplacement(const int axis) const {
    return (axis >= 0 && axis < 3) ? frame_displacement_[axis] : 0.0;
  }
  std::array<Real, 3> FrameVelocity() const { return frame_velocity_; }
  std::array<Real, 3> FrameDisplacement() const { return frame_displacement_; }
  Real LastTopMdot() const { return mdot_last_top_; }
  Real LastTopMdotAverage() const { return mdot_last_average_; }
  Real TopMdotAverageWeightTime() const { return mdot_weight_time_; }
  Real TopMdotArea() const { return mdot_area_; }
  int FillHistoryData(std::string labels[], Real values[], const int max_values) const;
  void StoreStateInParameterInput(ParameterInput *pin) const;

 private:
  struct AxisState {
    bool active = false;
    Real target_position = 0.0;
    bool use_lower_limit = false;
    Real lower_limit = 0.0;
    bool use_upper_limit = false;
    Real upper_limit = 0.0;

    Real last_boost = 0.0;
    Real cumulative_boost = 0.0;
    Real last_mean_x = 0.0;
    Real last_mean_v = 0.0;
    Real last_global_weight = 0.0;
    Real last_filtered_x = 0.0;
    Real last_filtered_v = 0.0;
    Real last_x_err = 0.0;
    Real last_vel_cmd_pre = 0.0;
    Real last_vel_cmd_post = 0.0;
    Real last_x_centroid = 0.0;
    Real last_x_midpoint = 0.0;
    Real last_x_ctrl = 0.0;
    Real i_term = 0.0;
    bool last_skip_flag = false;
    bool last_slew_limited = false;
    bool filter_initialized = false;
  };

  struct MomentSample {
    Real global_weight = 0.0;
    Real mean_x = 0.0;
    Real mean_v = 0.0;
    Real x_min = 0.0;
    Real x_max = 0.0;
    bool valid = false;
  };

  void ParseAxes(ParameterInput *pin);
  void ParseAxisControls(ParameterInput *pin);
  void ParseTrackedFluid(ParameterInput *pin);
  void RejectRemovedAliases(ParameterInput *pin) const;
  void PrintConfigurationSummary() const;
  void ValidateConfiguration();
  void SetActiveTargetRange();
  bool SampleMoments(std::array<MomentSample, 3> &samples) const;
  Real SampleTopMassFlux(const bool lab_frame) const;
  bool ApplyTopMdotTracking(const Real dt_boost);
  bool ApplyTopMdotAutoTracking(const Real dt_boost);
  bool ApplyTopMdotWindowTracking(const Real dt_boost);
  bool ApplyX3FrameVelocityTarget(const Real target_frame_velocity,
                                  const Real dt_boost);
  void PrintSkipMessage(const bool have_sample, const bool low_weight_floor,
                        const Real global_weight) const;
  void PrintPrimeMessage(const std::array<MomentSample, 3> &samples) const;
  void PrintApplyMessage(const std::array<MomentSample, 3> &samples,
                         const std::array<Real, 3> &v_p,
                         const std::array<Real, 3> &v_d,
                         const std::array<Real, 3> &v_i,
                         const std::array<Real, 3> &cmd_pre,
                         const std::array<Real, 3> &boost) const;
  void PrintTopMdotMessage(const Real top_mdot, const Real avg_mdot,
                           const Real target_frame_velocity,
                           const Real ramp_factor, const Real boost) const;
  void PrintTopMdotAutoMessage(const Real boost) const;
  void PrintTopMdotWindowMessage(const Real boost) const;
  bool AdvanceFrameDisplacement(const Real dt);
  void RestoreFrameState(ParameterInput *pin);

  MeshBlockPack *pmy_pack = nullptr;
  std::string block_name_;

  bool enabled_ = true;
  std::string target_name_ = "density";
  std::string mode_name_ = "pd";
  std::string boost_change_mode_name_ = "per_apply";
  int apply_every_ = 1;
  Real start_time_ = 0.0;
  int diagnostic_every_ = -1;
  bool verbose_ = false;

  int target_kind_ = 0;
  int scalar_index_ = 0;
  int weight_mode_ = 0;
  bool track_mhd_ = false;
  Real target_min_ = 0.0;
  Real target_max_ = 0.0;
  Real target_center_ = 0.0;
  Real target_min_active_ = 0.0;
  Real target_max_active_ = 0.0;
  bool use_log_reacquire_ = false;

  int mode_ = 0;
  int position_signal_ = 0;
  std::string position_signal_name_ = "blend";
  int velocity_signal_ = 0;
  std::string velocity_signal_name_ = "material_mean";
  int boost_change_mode_ = 0;
  Real position_blend_ = 0.0;
  Real tau_avg_ = 1.0;
  Real tau_relax_ = 1.0;
  Real tau_vel_ = 1.0;
  Real tau_int_ = 0.0;
  Real int_max_abs_ = 0.0;
  Real int_leak_tau_ = 1.0;
  bool int_unsat_only_ = true;
  Real weight_floor_ = 0.0;
  Real min_global_weight_ = 0.0;
  Real max_abs_boost_ = 0.0;
  Real max_boost_change_ = 0.0;
  Real max_boost_change_rate_ = -1.0;
  Real mdot_average_start_time_ = 0.0;
  Real mdot_average_end_time_ = 0.0;
  Real mdot_ramp_time_ = 0.0;
  Real mdot_density_ = 1.0;
  Real mdot_area_ = -1.0;
  std::string mdot_density_mode_name_ = "fixed";
  std::string mdot_velocity_sign_name_ = "normal";
  std::string mdot_flux_frame_name_ = "grid";
  bool mdot_flux_lab_frame_ = false;
  Real mdot_velocity_sign_ = 1.0;
  Real mdot_velocity_factor_ = 1.0;
  Real mdot_integral_ = 0.0;
  Real mdot_weight_time_ = 0.0;
  Real mdot_last_top_ = 0.0;
  Real mdot_last_average_ = 0.0;
  Real mdot_last_vtop_ = 0.0;
  Real mdot_last_ramp_factor_ = 0.0;
  Real mdot_window_start_time_ = 0.5;
  Real mdot_window_time_ = 1.0;
  Real mdot_window_low_time_ = 0.0;
  Real mdot_window_target_vframe_ = 0.0;
  std::deque<Real> mdot_window_times_;
  std::deque<Real> mdot_window_values_;
  Real mdot_auto_detect_start_time_ = 0.5;
  Real mdot_auto_fast_tau_ = 0.5;
  Real mdot_auto_slow_tau_ = 1.5;
  Real mdot_auto_rel_tol_ = 0.075;
  Real mdot_auto_abs_floor_ = 1.0e-12;
  Real mdot_auto_hold_time_ = 0.5;
  Real mdot_auto_ramp_fraction_ = 0.25;
  Real mdot_auto_ramp_min_time_ = 0.25;
  Real mdot_auto_ramp_max_time_ = 5.0;
  Real mdot_auto_fast_ = 0.0;
  Real mdot_auto_slow_ = 0.0;
  Real mdot_auto_rel_ = 0.0;
  Real mdot_auto_hold_elapsed_ = 0.0;
  Real mdot_auto_settle_time_ = -1.0;
  Real mdot_auto_locked_mdot_ = 0.0;
  Real mdot_auto_locked_vframe_ = 0.0;
  Real mdot_auto_frame_velocity_at_lock_ = 0.0;
  Real mdot_auto_ramp_duration_ = 0.0;
  Real mdot_auto_target_vframe_ = 0.0;
  Real mdot_auto_ramp_factor_ = 0.0;
  bool mdot_auto_filter_initialized_ = false;
  bool mdot_auto_locked_ = false;
  Real reacquire_expand_factor_ = 1.0;
  Real reacquire_max_expand_ = 1.0;
  int reacquire_recover_updates_ = 1;
  bool reacquire_leak_on_miss_ = true;
  bool boundary_guard_enable_ = false;
  int boundary_guard_cells_ = 8;
  Real boundary_guard_min_scale_ = 0.15;

  Real last_apply_time_ = -1.0;
  int miss_streak_ = 0;
  int recover_streak_ = 0;
  int restored_state_kind_ = 0;  // 0: new, 1: versioned restart, 2: legacy restart
  std::array<AxisState, 3> axes_;
  std::array<Real, 3> frame_velocity_ = {0.0, 0.0, 0.0};
  std::array<Real, 3> frame_displacement_ = {0.0, 0.0, 0.0};
};

#endif  // SRCTERMS_FRAME_TRACKER_HPP_
