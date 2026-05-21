#ifndef SRCTERMS_COOLING_HPP_
#define SRCTERMS_COOLING_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cooling.hpp
//! \brief General radiative cooling and heating source term selected by <cooling>.

#include <array>
#include <string>

#include "athena.hpp"

class MeshBlockPack;
class ParameterInput;
struct EOS_Data;

namespace cooling {

constexpr int MAX_TABLE_AXES = 3;
constexpr int MAX_POWERLAW_PIECES = 8;

enum class UnitSystem { code, cgs };
enum class DensityKind { mass_density, number_density, hydrogen_number_density };
enum class AxisKind { temperature, density, scalar };
enum class ValueScale { linear, log10 };
enum class BoundsBehavior { zero, clamp, fatal };
enum class ModelKind { none, constant, ism, cgm, table, powerlaw, piecewise_powerlaw, user };
enum class TableValueKind { lambda, gamma, modifier };

struct CompositionData {
  Real temperature_mu = 1.0;
  Real mu = 1.0;
  Real hydrogen_mass_fraction = 1.0;
  bool has_temperature_mu = false;
  bool has_mu = false;
  bool has_hydrogen_mass_fraction = false;
};

struct AxisData {
  AxisKind kind = AxisKind::temperature;
  UnitSystem units = UnitSystem::code;
  ValueScale scale = ValueScale::log10;
  DensityKind density_kind = DensityKind::mass_density;
  int scalar_index = 0;
  Real xmin = 0.0;
  Real xmax = 1.0;
  int n = 1;
};

struct TableData {
  bool enabled = false;
  int ndim = 0;
  TableValueKind value_kind = TableValueKind::lambda;
  BoundsBehavior bounds = BoundsBehavior::zero;
  UnitSystem value_units = UnitSystem::code;
  ValueScale value_scale = ValueScale::linear;
  AxisData axes[MAX_TABLE_AXES];
  int n0 = 1, n1 = 1, n2 = 1;
  DvceArray1D<Real> values;
};

struct PowerLawData {
  bool enabled = false;
  bool piecewise = false;
  AxisKind axis = AxisKind::temperature;
  UnitSystem axis_units = UnitSystem::code;
  DensityKind density_kind = DensityKind::mass_density;
  int scalar_index = 0;
  UnitSystem value_units = UnitSystem::code;
  BoundsBehavior bounds = BoundsBehavior::zero;
  Real reference_axis = 1.0;
  Real reference_value = 0.0;
  Real slope = 0.0;
  int nbreaks = 0;
  Real breaks[MAX_POWERLAW_PIECES - 1] = {};
  Real slopes[MAX_POWERLAW_PIECES] = {};
};

struct CoolingCellState {
  Real rho_code = 0.0;
  Real rho_cgs = 0.0;
  Real temp_code = 0.0;
  Real temp_cgs = 0.0;
  Real eint_code = 0.0;
  Real scalar0 = 0.0;
  int nfluid = 0;
  int nscalars = 0;
};

struct CoolingRates {
  Real lambda = 0.0;
  Real gamma = 0.0;
  UnitSystem lambda_units = UnitSystem::code;
  UnitSystem gamma_units = UnitSystem::code;
};

struct RuntimeData {
  UnitSystem default_units = UnitSystem::code;
  DensityKind cooling_density = DensityKind::mass_density;
  DensityKind heating_density = DensityKind::mass_density;
  CompositionData composition;
  Real density_cgs = 1.0;
  Real length_cgs = 1.0;
  Real pressure_cgs = 1.0;
  Real time_cgs = 1.0;
  Real temp_cgs = 1.0;
  Real edot_cgs_to_code = 1.0;
};

struct TimestepData {
  Real factor = 1.0;
  bool use_temperature_bounds = false;
  UnitSystem temperature_units = UnitSystem::code;
  Real temperature_min = 0.0;
  Real temperature_max = 0.0;
  bool use_density_bounds = false;
  UnitSystem density_units = UnitSystem::code;
  DensityKind density_kind = DensityKind::mass_density;
  Real density_min = 0.0;
  Real density_max = 0.0;
};

struct ModifierData {
  bool enabled = false;
  ModelKind model = ModelKind::none;
  BoundsBehavior bounds = BoundsBehavior::zero;
  Real constant = 1.0;
  TableData table;
  PowerLawData powerlaw;
};

struct ShieldingData {
  bool enabled = false;
  bool apply_to_heating = false;
  DensityKind density_kind = DensityKind::hydrogen_number_density;
  Real cross_section_cgs = 1.0e-17;
  Real transition_temperature_cgs = 8.0e3;
  Real transition_width_cgs = 1.5e3;
  Real length_code = 0.0;
};

class GeneralCooling {
 public:
  GeneralCooling(MeshBlockPack *pp, ParameterInput *pin);
  ~GeneralCooling() = default;

  bool Enabled() const { return enabled_; }
  bool TimestepEnabled() const { return enabled_ && timestep_enabled_; }
  bool HistoryEnabled() const { return enabled_ && history_enabled_; }

  void Apply(const DvceArray5D<Real> &w0, const EOS_Data &eos_data, const Real bdt,
             const Real history_bdt, DvceArray5D<Real> &u0);
  void NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data, Real &dtnew);

  int AddHistoryLabels(std::string *labels, int start, int max_labels) const;
  int AddHistoryData(Real *hdata, int start, int max_data, Real current_time);

 private:
  MeshBlockPack *pmy_pack_;
  bool enabled_ = false;
  bool timestep_enabled_ = false;
  bool history_enabled_ = false;
  bool history_gross_ = true;
  bool history_net_ = true;
  Real timestep_factor_ = 1.0;
  bool timestep_use_temperature_bounds_ = false;
  UnitSystem timestep_temperature_units_ = UnitSystem::code;
  Real timestep_temperature_min_ = 0.0;
  Real timestep_temperature_max_ = 0.0;
  bool timestep_use_density_bounds_ = false;
  UnitSystem timestep_density_units_ = UnitSystem::code;
  DensityKind timestep_density_kind_ = DensityKind::mass_density;
  Real timestep_density_min_ = 0.0;
  Real timestep_density_max_ = 0.0;

  ModelKind cooling_model_ = ModelKind::none;
  ModelKind heating_model_ = ModelKind::none;
  BoundsBehavior cooling_bounds_ = BoundsBehavior::zero;
  BoundsBehavior heating_bounds_ = BoundsBehavior::zero;
  RuntimeData runtime_;
  TableData cooling_table_;
  TableData heating_table_;
  TableData cgm_pie_table_;
  TableData cgm_cie_table_;
  PowerLawData cooling_powerlaw_;
  PowerLawData heating_powerlaw_;
  ModifierData cooling_modifier_;
  ModifierData heating_modifier_;
  ShieldingData shielding_;
  Real constant_heating_rate_ = 0.0;

  Real accumulated_gross_energy_ = 0.0;
  Real accumulated_net_energy_ = 0.0;
  Real last_history_time_ = 0.0;
};

} // namespace cooling

#endif // SRCTERMS_COOLING_HPP_
