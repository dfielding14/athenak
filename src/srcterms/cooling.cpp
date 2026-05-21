//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cooling.cpp
//! \brief General cooling source term selected by <cooling>.

#include "cooling.hpp"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "cooling_hooks.hpp"
#include "hydro/hydro.hpp"
#include "ismcooling.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "units/units.hpp"

namespace cooling {
namespace {

void FatalCoolingInput(const std::string &message) {
  std::cout << "### FATAL ERROR in <cooling> input" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

bool Has(ParameterInput *pin, const std::string &name) {
  return pin->DoesParameterExist("cooling", name);
}

UnitSystem ParseUnitSystem(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "code") return UnitSystem::code;
  if (v == "cgs") return UnitSystem::cgs;
  FatalCoolingInput(context + " must be 'code' or 'cgs'.");
  return UnitSystem::code;
}

DensityKind ParseDensityKind(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "mass_density") return DensityKind::mass_density;
  if (v == "number_density") return DensityKind::number_density;
  if (v == "hydrogen_number_density") return DensityKind::hydrogen_number_density;
  FatalCoolingInput(context + " has unknown density kind '" + value + "'.");
  return DensityKind::mass_density;
}

AxisKind ParseAxisKind(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "temperature") return AxisKind::temperature;
  if (v == "density") return AxisKind::density;
  if (v == "scalar" || v == "metallicity") return AxisKind::scalar;
  FatalCoolingInput(context + " has unknown axis kind '" + value + "'.");
  return AxisKind::temperature;
}

ValueScale ParseValueScale(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "linear") return ValueScale::linear;
  if (v == "log10") return ValueScale::log10;
  FatalCoolingInput(context + " must be 'linear' or 'log10'.");
  return ValueScale::linear;
}

BoundsBehavior ParseBounds(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "zero") return BoundsBehavior::zero;
  if (v == "clamp") return BoundsBehavior::clamp;
  if (v == "fatal") return BoundsBehavior::fatal;
  FatalCoolingInput(context + " must be 'zero', 'clamp', or 'fatal'.");
  return BoundsBehavior::zero;
}

TableValueKind ParseTableValueKind(const std::string &value,
                                   const std::string &context) {
  std::string v = Lower(value);
  if (v == "lambda") return TableValueKind::lambda;
  if (v == "gamma") return TableValueKind::gamma;
  if (v == "modifier") return TableValueKind::modifier;
  FatalCoolingInput(context + " must be 'Lambda', 'Gamma', or 'modifier'.");
  return TableValueKind::lambda;
}

ModelKind ParseCoolingModel(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "none") return ModelKind::none;
  if (v == "ism") return ModelKind::ism;
  if (v == "cgm") return ModelKind::cgm;
  if (v == "table") return ModelKind::table;
  if (v == "powerlaw") return ModelKind::powerlaw;
  if (v == "piecewise_powerlaw") return ModelKind::piecewise_powerlaw;
  if (v == "user") return ModelKind::user;
  FatalCoolingInput(context + " has unknown model '" + value + "'.");
  return ModelKind::none;
}

ModelKind ParseModifierModel(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "none") return ModelKind::none;
  if (v == "constant") return ModelKind::constant;
  if (v == "table") return ModelKind::table;
  if (v == "powerlaw") return ModelKind::powerlaw;
  if (v == "piecewise_powerlaw") return ModelKind::piecewise_powerlaw;
  FatalCoolingInput(context + " must be 'none', 'constant', 'table', 'powerlaw', "
                    "or 'piecewise_powerlaw'.");
  return ModelKind::none;
}

ModelKind ParseHeatingModel(const std::string &value, const std::string &context) {
  std::string v = Lower(value);
  if (v == "none") return ModelKind::none;
  if (v == "constant") return ModelKind::constant;
  if (v == "table") return ModelKind::table;
  if (v == "powerlaw") return ModelKind::powerlaw;
  if (v == "piecewise_powerlaw") return ModelKind::piecewise_powerlaw;
  if (v == "user") return ModelKind::user;
  FatalCoolingInput(context + " has unknown model '" + value + "'.");
  return ModelKind::none;
}

std::vector<Real> ParseRealList(const std::string &value, const std::string &context) {
  std::vector<Real> values;
  std::string normalized = value;
  std::replace(normalized.begin(), normalized.end(), ',', ' ');
  std::istringstream stream(normalized);
  Real next = 0.0;
  while (stream >> next) values.push_back(next);
  if (values.empty() && !value.empty()) {
    FatalCoolingInput("Could not parse real list for " + context + ".");
  }
  return values;
}

Real DensityValue(DensityKind kind, UnitSystem units, Real rho_code,
                  const RuntimeData &runtime) {
  if (units == UnitSystem::cgs) {
    const Real rho_cgs = rho_code*runtime.density_cgs;
    if (kind == DensityKind::mass_density) return rho_cgs;
    if (kind == DensityKind::number_density) {
      return rho_cgs/(runtime.composition.mu*units::Units::atomic_mass_unit_cgs);
    }
    return rho_cgs*runtime.composition.hydrogen_mass_fraction/
           units::Units::atomic_mass_unit_cgs;
  }

  if (kind == DensityKind::mass_density) return rho_code;
  if (kind == DensityKind::number_density) return rho_code/runtime.composition.mu;
  return rho_code*runtime.composition.hydrogen_mass_fraction;
}

Real RawAxisValue(const AxisData &axis, const CoolingCellState &state,
                  const DvceArray5D<Real> &w0, const RuntimeData &runtime,
                  int m, int k, int j, int i) {
  Real value = 0.0;
  if (axis.kind == AxisKind::temperature) {
    value = (axis.units == UnitSystem::cgs) ? state.temp_cgs : state.temp_code;
  } else if (axis.kind == AxisKind::density) {
    value = DensityValue(axis.density_kind, axis.units, state.rho_code, runtime);
  } else {
    if (axis.scalar_index < state.nscalars) {
      value = w0(m, state.nfluid + axis.scalar_index, k, j, i);
    }
  }
  return value;
}

Real AxisValue(const AxisData &axis, const CoolingCellState &state,
               const DvceArray5D<Real> &w0, const RuntimeData &runtime,
               int m, int k, int j, int i) {
  Real value = RawAxisValue(axis, state, w0, runtime, m, k, j, i);
  if (axis.scale == ValueScale::log10) value = log10(value);
  return value;
}

Real TransformValue(Real value, ValueScale scale) {
  return (scale == ValueScale::log10) ? pow(10.0, value) : value;
}

Real CoolingEdotCode(Real lambda, UnitSystem units, DensityKind density_kind,
                     Real rho_code, const RuntimeData &runtime) {
  const Real q = DensityValue(density_kind, units, rho_code, runtime);
  Real edot = q*q*lambda;
  if (units == UnitSystem::cgs) edot *= runtime.edot_cgs_to_code;
  return edot;
}

Real HeatingEdotCode(Real gamma, UnitSystem units, DensityKind density_kind,
                     Real rho_code, const RuntimeData &runtime) {
  const Real q = DensityValue(density_kind, units, rho_code, runtime);
  Real edot = q*gamma;
  if (units == UnitSystem::cgs) edot *= runtime.edot_cgs_to_code;
  return edot;
}

Real EvaluateTable(const TableData &table, const CoolingCellState &state,
                   const DvceArray5D<Real> &w0, const RuntimeData &runtime,
                   int m, int k, int j, int i, Real &bad_bounds) {
  if (!table.enabled || table.ndim <= 0) return 0.0;

  Real coord[MAX_TABLE_AXES] = {0.0, 0.0, 0.0};
  int idx0[MAX_TABLE_AXES] = {0, 0, 0};
  Real frac[MAX_TABLE_AXES] = {0.0, 0.0, 0.0};

  for (int a = 0; a < table.ndim; ++a) {
    const AxisData &axis = table.axes[a];
    Real raw_coord = RawAxisValue(axis, state, w0, runtime, m, k, j, i);
    if (axis.scale == ValueScale::log10) {
      if (raw_coord <= 0.0) {
        if (table.bounds == BoundsBehavior::fatal) bad_bounds += 1.0;
        if (table.bounds == BoundsBehavior::zero ||
            table.bounds == BoundsBehavior::fatal) {
          return 0.0;
        }
        coord[a] = axis.xmin;
      } else {
        coord[a] = log10(raw_coord);
      }
    } else {
      coord[a] = raw_coord;
    }
    const bool low = coord[a] < axis.xmin;
    const bool high = coord[a] > axis.xmax;
    if (low || high) {
      if (table.bounds == BoundsBehavior::fatal) bad_bounds += 1.0;
      if (table.bounds == BoundsBehavior::zero ||
          table.bounds == BoundsBehavior::fatal) return 0.0;
      coord[a] = fmin(fmax(coord[a], axis.xmin), axis.xmax);
    }
    if (axis.n <= 1) {
      idx0[a] = 0;
      frac[a] = 0.0;
    } else {
      const Real x = (coord[a] - axis.xmin)/(axis.xmax - axis.xmin)*(axis.n - 1);
      int lo = static_cast<int>(floor(x));
      lo = (lo < 0) ? 0 : lo;
      lo = (lo > axis.n - 2) ? axis.n - 2 : lo;
      idx0[a] = lo;
      frac[a] = x - static_cast<Real>(lo);
    }
  }

  const int n1 = table.n1;
  const int n2 = table.n2;
  auto values = table.values;
  Real interpolated = 0.0;
  if (table.ndim == 1) {
    const int i0 = idx0[0];
    const Real t = frac[0];
    interpolated = (1.0 - t)*values(i0) + t*values(i0 + 1);
  } else if (table.ndim == 2) {
    const int i0 = idx0[0];
    const int j0 = idx0[1];
    const Real t = frac[0];
    const Real u = frac[1];
    const int p00 = i0*n1 + j0;
    const int p10 = (i0 + 1)*n1 + j0;
    const int p01 = i0*n1 + (j0 + 1);
    const int p11 = (i0 + 1)*n1 + (j0 + 1);
    interpolated = (1.0 - t)*(1.0 - u)*values(p00) +
                   t*(1.0 - u)*values(p10) +
                   (1.0 - t)*u*values(p01) +
                   t*u*values(p11);
  } else {
    const int i0 = idx0[0];
    const int j0 = idx0[1];
    const int k0 = idx0[2];
    const Real t = frac[0];
    const Real u = frac[1];
    const Real v = frac[2];
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        for (int dk = 0; dk <= 1; ++dk) {
          const Real wi = di ? t : (1.0 - t);
          const Real wj = dj ? u : (1.0 - u);
          const Real wk = dk ? v : (1.0 - v);
          const int p = ((i0 + di)*n1 + (j0 + dj))*n2 + (k0 + dk);
          interpolated += wi*wj*wk*values(p);
        }
      }
    }
  }
  return TransformValue(interpolated, table.value_scale);
}

Real EvaluatePowerLaw(const PowerLawData &powerlaw, const CoolingCellState &state,
                      const DvceArray5D<Real> &w0, const RuntimeData &runtime,
                      int m, int k, int j, int i, Real &bad_bounds) {
  if (!powerlaw.enabled) return 0.0;
  AxisData axis;
  axis.kind = powerlaw.axis;
  axis.units = powerlaw.axis_units;
  axis.scale = ValueScale::linear;
  axis.density_kind = powerlaw.density_kind;
  axis.scalar_index = powerlaw.scalar_index;
  Real x = AxisValue(axis, state, w0, runtime, m, k, j, i);
  if (x <= 0.0 || powerlaw.reference_axis <= 0.0) {
    if (powerlaw.bounds == BoundsBehavior::fatal) bad_bounds += 1.0;
    return 0.0;
  }
  if (!powerlaw.piecewise) {
    return powerlaw.reference_value*pow(x/powerlaw.reference_axis, powerlaw.slope);
  }

  Real value = powerlaw.reference_value;
  Real x_anchor = powerlaw.reference_axis;
  int piece = 0;
  while (piece < powerlaw.nbreaks && x_anchor >= powerlaw.breaks[piece]) ++piece;
  if (x >= x_anchor) {
    while (piece < powerlaw.nbreaks && x > powerlaw.breaks[piece]) {
      value *= pow(powerlaw.breaks[piece]/x_anchor, powerlaw.slopes[piece]);
      x_anchor = powerlaw.breaks[piece];
      ++piece;
    }
  } else {
    while (piece > 0 && x < powerlaw.breaks[piece - 1]) {
      value *= pow(powerlaw.breaks[piece - 1]/x_anchor, powerlaw.slopes[piece]);
      x_anchor = powerlaw.breaks[piece - 1];
      --piece;
    }
  }
  return value*pow(x/x_anchor, powerlaw.slopes[piece]);
}

void ParsePowerLaw(ParameterInput *pin, const std::string &prefix, ModelKind model,
                   UnitSystem default_units, BoundsBehavior bounds,
                   PowerLawData &powerlaw) {
  if (model != ModelKind::powerlaw && model != ModelKind::piecewise_powerlaw) return;
  powerlaw.enabled = true;
  powerlaw.piecewise = (model == ModelKind::piecewise_powerlaw);
  powerlaw.axis = ParseAxisKind(pin->GetOrAddString("cooling", prefix + "_powerlaw_axis",
                                                    "temperature"),
                                prefix + "_powerlaw_axis");
  powerlaw.axis_units = ParseUnitSystem(pin->GetOrAddString("cooling",
                                      prefix + "_powerlaw_axis_units",
                                      default_units == UnitSystem::cgs ? "cgs" : "code"),
                                      prefix + "_powerlaw_axis_units");
  powerlaw.value_units = ParseUnitSystem(pin->GetOrAddString("cooling",
                                      prefix + "_powerlaw_value_units",
                                      default_units == UnitSystem::cgs ? "cgs" : "code"),
                                      prefix + "_powerlaw_value_units");
  powerlaw.bounds = bounds;
  if (powerlaw.axis == AxisKind::density) {
    powerlaw.density_kind = ParseDensityKind(pin->GetOrAddString("cooling",
                                      prefix + "_powerlaw_density_kind",
                                      "mass_density"),
                                      prefix + "_powerlaw_density_kind");
  }
  if (powerlaw.axis == AxisKind::scalar) {
    powerlaw.scalar_index = pin->GetOrAddInteger("cooling",
                                                 prefix + "_powerlaw_scalar_index", 0);
  }
  powerlaw.reference_axis = pin->GetReal("cooling", prefix + "_reference_axis");
  powerlaw.reference_value = pin->GetReal("cooling", prefix + "_reference_value");
  powerlaw.slope = pin->GetOrAddReal("cooling", prefix + "_slope", 0.0);
  if (powerlaw.piecewise) {
    const std::vector<Real> breaks = ParseRealList(pin->GetString("cooling",
        prefix + "_breaks"), prefix + "_breaks");
    const std::vector<Real> slopes = ParseRealList(pin->GetString("cooling",
        prefix + "_slopes"), prefix + "_slopes");
    if (breaks.size() + 1 != slopes.size()) {
      FatalCoolingInput(prefix + "_slopes must contain exactly one more entry than " +
                        prefix + "_breaks.");
    }
    if (slopes.size() > MAX_POWERLAW_PIECES) {
      FatalCoolingInput(prefix + " piecewise power law exceeds MAX_POWERLAW_PIECES.");
    }
    powerlaw.nbreaks = static_cast<int>(breaks.size());
    for (int n = 0; n < powerlaw.nbreaks; ++n) {
      if (n > 0 && breaks[n] <= breaks[n - 1]) {
        FatalCoolingInput(prefix + "_breaks must be strictly increasing.");
      }
      powerlaw.breaks[n] = breaks[n];
    }
    for (std::size_t n = 0; n < slopes.size(); ++n) powerlaw.slopes[n] = slopes[n];
  }
}

void ParseAxisExtra(const std::string &token, AxisData &axis) {
  const std::size_t eq = token.find('=');
  if (eq == std::string::npos) {
    if (axis.kind == AxisKind::density) {
      axis.density_kind = ParseDensityKind(token, "table density axis");
    }
    return;
  }
  const std::string key = token.substr(0, eq);
  const std::string value = token.substr(eq + 1);
  if (key == "scalar_index") {
    axis.scalar_index = std::stoi(value);
  } else if (key == "density_kind") {
    axis.density_kind = ParseDensityKind(value, "table density axis");
  }
}

bool IsValueScaleToken(const std::string &value) {
  const std::string v = Lower(value);
  return (v == "linear" || v == "log10");
}

void LoadTable(const std::string &path, BoundsBehavior default_bounds,
               const std::string &expected_value_kind, TableData &table) {
  std::ifstream in(path);
  if (!in) FatalCoolingInput("Could not open cooling table '" + path + "'.");

  std::string token;
  int version = 0;
  in >> token >> version;
  if (token != "ATHENAK_COOLING_TABLE" || version != 1) {
    FatalCoolingInput("Cooling table '" + path + "' does not start with "
                      "ATHENAK_COOLING_TABLE 1.");
  }
  table.bounds = default_bounds;
  const TableValueKind expected_kind =
      ParseTableValueKind(expected_value_kind, "expected table value_kind");
  bool saw_data = false;
  bool saw_value_kind = false;
  bool axis_seen[MAX_TABLE_AXES] = {false, false, false};
  while (in >> token) {
    if (token.empty()) continue;
    if (token[0] == '#') {
      std::string rest;
      std::getline(in, rest);
      continue;
    }
    if (token == "data") {
      saw_data = true;
      break;
    } else if (token == "ndim") {
      in >> table.ndim;
      if (table.ndim < 1 || table.ndim > MAX_TABLE_AXES) {
        FatalCoolingInput("Cooling table '" + path + "' has invalid ndim.");
      }
    } else if (token == "value_kind") {
      std::string value_kind;
      in >> value_kind;
      table.value_kind = ParseTableValueKind(value_kind, "table value_kind");
      saw_value_kind = true;
      if (table.value_kind != expected_kind) {
        FatalCoolingInput("Cooling table '" + path + "' has value_kind '" +
                          value_kind + "' but this use requires value_kind " +
                          expected_value_kind + ".");
      }
    } else if (token == "value_units") {
      std::string units;
      in >> units;
      table.value_units = ParseUnitSystem(units, "table value_units");
    } else if (token == "value_scale") {
      std::string scale;
      in >> scale;
      table.value_scale = ParseValueScale(scale, "table value_scale");
    } else if (token == "bounds") {
      std::string bounds;
      in >> bounds;
      table.bounds = ParseBounds(bounds, "table bounds");
    } else if (token == "data_order") {
      std::string order;
      in >> order;
      if (order != "c_row_major_last_axis_fastest") {
        FatalCoolingInput("Cooling table '" + path + "' uses unsupported data_order.");
      }
    } else if (token.rfind("axis", 0) == 0) {
      const int axis_index = std::stoi(token.substr(4));
      if (axis_index < 0 || axis_index >= MAX_TABLE_AXES) {
        FatalCoolingInput("Cooling table '" + path + "' has invalid axis index.");
      }
      AxisData axis;
      std::string rest;
      std::getline(in, rest);
      std::istringstream axis_line(rest);
      std::vector<std::string> fields;
      std::string field;
      while (axis_line >> field) fields.push_back(field);
      if (fields.size() < 5) {
        FatalCoolingInput("Cooling table '" + path + "' has incomplete axis line.");
      }
      axis.kind = ParseAxisKind(fields[0], "table axis kind");
      axis.units = ParseUnitSystem(fields[1], "table axis units");
      int next = 2;
      axis.scale = ValueScale::log10;
      if (IsValueScaleToken(fields[next])) {
        axis.scale = ParseValueScale(fields[next], "table axis scale");
        ++next;
      }
      if (static_cast<int>(fields.size()) < next + 3) {
        FatalCoolingInput("Cooling table '" + path + "' has incomplete axis extent.");
      }
      axis.xmin = std::stod(fields[next++]);
      axis.xmax = std::stod(fields[next++]);
      axis.n = std::stoi(fields[next++]);
      for (; next < static_cast<int>(fields.size()); ++next) {
        ParseAxisExtra(fields[next], axis);
      }
      if (axis.n < 1 || axis.xmax <= axis.xmin) {
        FatalCoolingInput("Cooling table '" + path + "' has invalid axis extent.");
      }
      if (axis.n < 2) {
        FatalCoolingInput("Cooling table '" + path + "' axes must have at least "
                          "two points.");
      }
      table.axes[axis_index] = axis;
      axis_seen[axis_index] = true;
    } else {
      FatalCoolingInput("Cooling table '" + path + "' has unknown header token '" +
                        token + "'.");
    }
  }
  if (!saw_data || table.ndim < 1 || !saw_value_kind) {
    FatalCoolingInput("Cooling table '" + path +
                      "' is missing data, ndim, or value_kind.");
  }
  for (int a = 0; a < table.ndim; ++a) {
    if (!axis_seen[a]) {
      FatalCoolingInput("Cooling table '" + path + "' is missing axis" +
                        std::to_string(a) + ".");
    }
  }
  for (int a = table.ndim; a < MAX_TABLE_AXES; ++a) {
    if (axis_seen[a]) {
      FatalCoolingInput("Cooling table '" + path + "' defines axis" +
                        std::to_string(a) + " beyond ndim.");
    }
  }

  table.n0 = table.axes[0].n;
  table.n1 = (table.ndim > 1) ? table.axes[1].n : 1;
  table.n2 = (table.ndim > 2) ? table.axes[2].n : 1;
  const int total = table.n0*table.n1*table.n2;
  std::vector<Real> host_values;
  host_values.reserve(total);
  Real value = 0.0;
  while (in >> value) host_values.push_back(value);
  if (static_cast<int>(host_values.size()) != total) {
    FatalCoolingInput("Cooling table '" + path + "' has " +
                      std::to_string(host_values.size()) + " values but expected " +
                      std::to_string(total) + ".");
  }
  table.values = DvceArray1D<Real>("cooling_table_values", total);
  auto host = Kokkos::create_mirror_view(table.values);
  for (int n = 0; n < total; ++n) {
    if (!std::isfinite(host_values[n])) {
      FatalCoolingInput("Cooling table '" + path + "' contains a non-finite value.");
    }
    host(n) = host_values[n];
  }
  Kokkos::deep_copy(table.values, host);
  table.enabled = true;
}

bool UsesCgsAxis(const TableData &table) {
  if (table.value_units == UnitSystem::cgs) return true;
  for (int a = 0; a < table.ndim; ++a) {
    if (table.axes[a].units == UnitSystem::cgs) return true;
  }
  return false;
}

bool UsesCgsTemperatureAxis(const TableData &table) {
  for (int a = 0; a < table.ndim; ++a) {
    if (table.axes[a].kind == AxisKind::temperature &&
        table.axes[a].units == UnitSystem::cgs) {
      return true;
    }
  }
  return false;
}

bool UsesCgsPowerLaw(const PowerLawData &powerlaw) {
  return powerlaw.enabled &&
         (powerlaw.axis_units == UnitSystem::cgs || powerlaw.value_units == UnitSystem::cgs);
}

bool UsesCgsTemperaturePowerLaw(const PowerLawData &powerlaw) {
  return powerlaw.enabled && powerlaw.axis == AxisKind::temperature &&
         powerlaw.axis_units == UnitSystem::cgs;
}

void ValidateDensityComposition(DensityKind density_kind, UnitSystem units,
                                const CompositionData &composition,
                                const std::string &context) {
  if (units != UnitSystem::cgs) return;
  if (density_kind == DensityKind::number_density && !composition.has_mu) {
    FatalCoolingInput(context + " uses cgs number_density and therefore requires "
                      "<cooling>/mu.");
  }
  if (density_kind == DensityKind::hydrogen_number_density &&
      !composition.has_hydrogen_mass_fraction) {
    FatalCoolingInput(context + " uses cgs hydrogen_number_density and therefore "
                      "requires <cooling>/hydrogen_mass_fraction or <cooling>/mu_H.");
  }
}

void ValidateTableComposition(const TableData &table, const CompositionData &composition,
                              const std::string &context) {
  for (int a = 0; a < table.ndim; ++a) {
    if (table.axes[a].kind == AxisKind::density) {
      ValidateDensityComposition(table.axes[a].density_kind, table.axes[a].units,
                                 composition, context + " density axis");
    }
  }
}

void ParseModifier(ParameterInput *pin, const std::string &prefix,
                   UnitSystem default_units, ModifierData &modifier) {
  modifier.model = ParseModifierModel(pin->GetOrAddString("cooling",
      prefix + "_modifier_model", "none"), prefix + "_modifier_model");
  if (modifier.model == ModelKind::none) return;
  modifier.enabled = true;
  modifier.bounds = ParseBounds(pin->GetOrAddString("cooling",
      prefix + "_modifier_bounds", "zero"), prefix + "_modifier_bounds");
  if (modifier.model == ModelKind::constant) {
    modifier.constant = pin->GetReal("cooling", prefix + "_modifier");
  } else if (modifier.model == ModelKind::table) {
    LoadTable(pin->GetString("cooling", prefix + "_modifier_table"),
              modifier.bounds, "modifier", modifier.table);
  } else {
    ParsePowerLaw(pin, prefix + "_modifier", modifier.model, default_units,
                  modifier.bounds, modifier.powerlaw);
  }
}

Real EvaluateModifier(const ModifierData &modifier, const CoolingCellState &state,
                      const DvceArray5D<Real> &w0, const RuntimeData &runtime,
                      int m, int k, int j, int i, Real &bad_bounds) {
  if (!modifier.enabled) return 1.0;
  if (modifier.model == ModelKind::constant) return modifier.constant;
  if (modifier.model == ModelKind::table) {
    return EvaluateTable(modifier.table, state, w0, runtime, m, k, j, i, bad_bounds);
  }
  return EvaluatePowerLaw(modifier.powerlaw, state, w0, runtime, m, k, j, i,
                          bad_bounds);
}

bool UsesCgsModifier(const ModifierData &modifier) {
  if (!modifier.enabled) return false;
  if (modifier.model == ModelKind::table) return UsesCgsAxis(modifier.table);
  return UsesCgsPowerLaw(modifier.powerlaw);
}

bool UsesCgsTemperatureModifier(const ModifierData &modifier) {
  if (!modifier.enabled) return false;
  if (modifier.model == ModelKind::table) return UsesCgsTemperatureAxis(modifier.table);
  return UsesCgsTemperaturePowerLaw(modifier.powerlaw);
}

void ValidateModifierComposition(const ModifierData &modifier,
                                 const CompositionData &composition,
                                 const std::string &context) {
  if (!modifier.enabled) return;
  if (modifier.model == ModelKind::table) {
    ValidateTableComposition(modifier.table, composition, context + " table");
  } else if ((modifier.model == ModelKind::powerlaw ||
              modifier.model == ModelKind::piecewise_powerlaw) &&
             modifier.powerlaw.axis == AxisKind::density) {
    ValidateDensityComposition(modifier.powerlaw.density_kind,
        modifier.powerlaw.axis_units, composition, context + " density axis");
  }
}

Real EvaluateShieldingFraction(const ShieldingData &shielding,
                               const CoolingCellState &state,
                               const RuntimeData &runtime, Real dx1_code) {
  if (!shielding.enabled) return 1.0;
  const Real n_density = DensityValue(shielding.density_kind, UnitSystem::cgs,
                                      state.rho_code, runtime);
  const Real neutral_frac =
      1.0 - 0.5*(1.0 + tanh((state.temp_cgs -
                             shielding.transition_temperature_cgs)/
                            shielding.transition_width_cgs));
  const Real length_code = (shielding.length_code > 0.0) ?
                           shielding.length_code : dx1_code;
  const Real tau = neutral_frac*n_density*shielding.cross_section_cgs*
                   length_code*runtime.length_cgs;
  return exp(-tau);
}

} // namespace

GeneralCooling::GeneralCooling(MeshBlockPack *pp, ParameterInput *pin) :
    pmy_pack_(pp) {
  enabled_ = pin->GetOrAddBoolean("cooling", "enabled", false);
  if (!enabled_) return;

  runtime_.default_units = ParseUnitSystem(pin->GetOrAddString("cooling", "units", "code"),
                                           "units");
  cooling_model_ = ParseCoolingModel(pin->GetOrAddString("cooling", "cooling_model",
                                                         "none"),
                                     "cooling_model");
  heating_model_ = ParseHeatingModel(pin->GetOrAddString("cooling", "heating_model",
                                                         "none"),
                                     "heating_model");
  if (cooling_model_ == ModelKind::none && heating_model_ == ModelKind::none) {
    FatalCoolingInput("enabled=true requires at least one cooling or heating model.");
  }

  cooling_bounds_ = ParseBounds(pin->GetOrAddString("cooling", "cooling_bounds", "zero"),
                                "cooling_bounds");
  heating_bounds_ = ParseBounds(pin->GetOrAddString("cooling", "heating_bounds", "zero"),
                                "heating_bounds");
  const char *default_cooling_density = "mass_density";
  if (cooling_model_ == ModelKind::ism) default_cooling_density = "number_density";
  if (cooling_model_ == ModelKind::cgm) {
    default_cooling_density = "hydrogen_number_density";
  }
  runtime_.cooling_density = ParseDensityKind(pin->GetOrAddString("cooling",
      "cooling_density", default_cooling_density), "cooling_density");
  runtime_.heating_density = ParseDensityKind(pin->GetOrAddString("cooling",
      "heating_density", "mass_density"), "heating_density");

  runtime_.composition.has_mu = Has(pin, "mu");
  if (runtime_.composition.has_mu) {
    runtime_.composition.mu = pin->GetReal("cooling", "mu");
  }
  runtime_.composition.has_temperature_mu = Has(pin, "temperature_mu");
  if (runtime_.composition.has_temperature_mu) {
    runtime_.composition.temperature_mu = pin->GetReal("cooling", "temperature_mu");
  } else if (runtime_.composition.has_mu) {
    runtime_.composition.temperature_mu = runtime_.composition.mu;
    runtime_.composition.has_temperature_mu = true;
  }
  if (Has(pin, "hydrogen_mass_fraction")) {
    runtime_.composition.hydrogen_mass_fraction =
        pin->GetReal("cooling", "hydrogen_mass_fraction");
    runtime_.composition.has_hydrogen_mass_fraction = true;
  }
  if (Has(pin, "mu_H")) {
    const Real mu_h = pin->GetReal("cooling", "mu_H");
    const Real x_from_mu_h = 1.0/mu_h;
    if (runtime_.composition.has_hydrogen_mass_fraction) {
      const Real rel = std::abs(runtime_.composition.hydrogen_mass_fraction -
                                x_from_mu_h)/x_from_mu_h;
      if (rel > 1.0e-10) {
        FatalCoolingInput("hydrogen_mass_fraction and mu_H are inconsistent.");
      }
    } else {
      runtime_.composition.hydrogen_mass_fraction = x_from_mu_h;
      runtime_.composition.has_hydrogen_mass_fraction = true;
    }
  }

  if (pmy_pack_->punit != nullptr) {
    runtime_.density_cgs = pmy_pack_->punit->density_cgs();
    runtime_.length_cgs = pmy_pack_->punit->length_cgs();
    runtime_.pressure_cgs = pmy_pack_->punit->pressure_cgs();
    runtime_.time_cgs = pmy_pack_->punit->time_cgs();
    const Real v_cgs = pmy_pack_->punit->velocity_cgs();
    runtime_.temp_cgs = v_cgs*v_cgs*runtime_.composition.temperature_mu*
                        units::Units::atomic_mass_unit_cgs/
                        units::Units::k_boltzmann_cgs;
    runtime_.edot_cgs_to_code = runtime_.time_cgs/runtime_.pressure_cgs;
  }

  timestep_enabled_ = pin->GetOrAddBoolean("cooling", "timestep", false);
  timestep_factor_ = pin->GetOrAddReal("cooling", "timestep_factor", 1.0);
  if (timestep_factor_ <= 0.0) FatalCoolingInput("timestep_factor must be positive.");
  timestep_use_temperature_bounds_ = pin->GetOrAddBoolean("cooling",
      "timestep_use_temperature_bounds", false);
  timestep_temperature_units_ = ParseUnitSystem(pin->GetOrAddString("cooling",
      "timestep_temperature_units", "code"), "timestep_temperature_units");
  if (timestep_use_temperature_bounds_) {
    timestep_temperature_min_ = pin->GetReal("cooling", "timestep_temperature_min");
    timestep_temperature_max_ = pin->GetReal("cooling", "timestep_temperature_max");
    if (timestep_temperature_max_ <= timestep_temperature_min_) {
      FatalCoolingInput("timestep_temperature_max must exceed timestep_temperature_min.");
    }
  }
  timestep_use_density_bounds_ = pin->GetOrAddBoolean("cooling",
      "timestep_use_density_bounds", false);
  timestep_density_units_ = ParseUnitSystem(pin->GetOrAddString("cooling",
      "timestep_density_units", "code"), "timestep_density_units");
  timestep_density_kind_ = ParseDensityKind(pin->GetOrAddString("cooling",
      "timestep_density_kind", "mass_density"), "timestep_density_kind");
  if (timestep_use_density_bounds_) {
    timestep_density_min_ = pin->GetReal("cooling", "timestep_density_min");
    timestep_density_max_ = pin->GetReal("cooling", "timestep_density_max");
    if (timestep_density_max_ <= timestep_density_min_) {
      FatalCoolingInput("timestep_density_max must exceed timestep_density_min.");
    }
  }

  history_enabled_ = pin->GetOrAddBoolean("cooling", "history", false);
  history_gross_ = pin->GetOrAddBoolean("cooling", "history_gross", true);
  history_net_ = pin->GetOrAddBoolean("cooling", "history_net", true);

  if (cooling_model_ == ModelKind::table) {
    LoadTable(pin->GetString("cooling", "cooling_table"), cooling_bounds_, "lambda",
              cooling_table_);
  }
  if (cooling_model_ == ModelKind::cgm) {
    LoadTable(pin->GetString("cooling", "cgm_pie_table"), cooling_bounds_, "lambda",
              cgm_pie_table_);
    LoadTable(pin->GetString("cooling", "cgm_cie_table"), cooling_bounds_, "lambda",
              cgm_cie_table_);
    if (cgm_pie_table_.value_units != cgm_cie_table_.value_units) {
      FatalCoolingInput("cgm_pie_table and cgm_cie_table must use the same value_units.");
    }
  }
  if (heating_model_ == ModelKind::table) {
    LoadTable(pin->GetString("cooling", "heating_table"), heating_bounds_, "gamma",
              heating_table_);
  }
  ParsePowerLaw(pin, "cooling", cooling_model_, runtime_.default_units, cooling_bounds_,
                cooling_powerlaw_);
  ParsePowerLaw(pin, "heating", heating_model_, runtime_.default_units, heating_bounds_,
                heating_powerlaw_);
  ParseModifier(pin, "cooling", runtime_.default_units, cooling_modifier_);
  ParseModifier(pin, "heating", runtime_.default_units, heating_modifier_);
  if (heating_model_ == ModelKind::constant) {
    constant_heating_rate_ = pin->GetReal("cooling", "heating_rate");
  }
  std::string shielding_model = pin->GetOrAddString("cooling", "shielding_model",
      cooling_model_ == ModelKind::cgm ? "cgm_approx" : "none");
  shielding_model = Lower(shielding_model);
  if (shielding_model == "cgm_approx") {
    shielding_.enabled = true;
    shielding_.apply_to_heating = pin->GetOrAddBoolean("cooling",
        "shielding_apply_to_heating", true);
    shielding_.density_kind = ParseDensityKind(pin->GetOrAddString("cooling",
        "shielding_density", "hydrogen_number_density"), "shielding_density");
    shielding_.cross_section_cgs = pin->GetOrAddReal("cooling",
        "shielding_cross_section_cgs", 1.0e-17);
    shielding_.transition_temperature_cgs = pin->GetOrAddReal("cooling",
        "shielding_transition_temperature_cgs", 8.0e3);
    shielding_.transition_width_cgs = pin->GetOrAddReal("cooling",
        "shielding_transition_width_cgs", 1.5e3);
    shielding_.length_code = pin->GetOrAddReal("cooling", "shielding_length_code", 0.0);
    if (shielding_.cross_section_cgs < 0.0 ||
        shielding_.transition_width_cgs <= 0.0 ||
        shielding_.length_code < 0.0) {
      FatalCoolingInput("shielding cross section, width, and length must be physical.");
    }
  } else if (shielding_model != "none") {
    FatalCoolingInput("shielding_model must be 'none' or 'cgm_approx'.");
  }
  if (history_enabled_) {
    last_history_time_ = pmy_pack_->pmesh->time;
  }

  bool needs_cgs = (runtime_.default_units == UnitSystem::cgs ||
                    cooling_model_ == ModelKind::ism ||
                    UsesCgsAxis(cooling_table_) || UsesCgsAxis(heating_table_) ||
                    UsesCgsAxis(cgm_pie_table_) || UsesCgsAxis(cgm_cie_table_) ||
                    UsesCgsPowerLaw(cooling_powerlaw_) ||
                    UsesCgsPowerLaw(heating_powerlaw_) ||
                    UsesCgsModifier(cooling_modifier_) ||
                    UsesCgsModifier(heating_modifier_) ||
                    shielding_.enabled ||
                    (timestep_use_temperature_bounds_ &&
                     timestep_temperature_units_ == UnitSystem::cgs) ||
                    (timestep_use_density_bounds_ &&
                     timestep_density_units_ == UnitSystem::cgs));
  if (needs_cgs && pmy_pack_->punit == nullptr) {
    FatalCoolingInput("cgs cooling/heating evaluation requires a <units> block.");
  }
  bool needs_cgs_temperature =
      (cooling_model_ == ModelKind::ism ||
       UsesCgsTemperatureAxis(cooling_table_) ||
       UsesCgsTemperatureAxis(heating_table_) ||
       UsesCgsTemperatureAxis(cgm_pie_table_) ||
       UsesCgsTemperatureAxis(cgm_cie_table_) ||
       UsesCgsTemperaturePowerLaw(cooling_powerlaw_) ||
       UsesCgsTemperaturePowerLaw(heating_powerlaw_) ||
       UsesCgsTemperatureModifier(cooling_modifier_) ||
       UsesCgsTemperatureModifier(heating_modifier_) ||
       shielding_.enabled ||
       (timestep_use_temperature_bounds_ &&
        timestep_temperature_units_ == UnitSystem::cgs));
  if (needs_cgs_temperature && !runtime_.composition.has_temperature_mu) {
    FatalCoolingInput("cgs temperature evaluation requires <cooling>/temperature_mu "
                      "or <cooling>/mu.");
  }
  UnitSystem cooling_density_units = runtime_.default_units;
  if (cooling_model_ == ModelKind::ism) cooling_density_units = UnitSystem::cgs;
  if (cooling_model_ == ModelKind::table) cooling_density_units = cooling_table_.value_units;
  if (cooling_model_ == ModelKind::cgm) cooling_density_units = cgm_pie_table_.value_units;
  if (cooling_model_ == ModelKind::powerlaw ||
      cooling_model_ == ModelKind::piecewise_powerlaw) {
    cooling_density_units = cooling_powerlaw_.value_units;
  }
  UnitSystem heating_density_units = runtime_.default_units;
  if (heating_model_ == ModelKind::table) heating_density_units = heating_table_.value_units;
  if (heating_model_ == ModelKind::powerlaw ||
      heating_model_ == ModelKind::piecewise_powerlaw) {
    heating_density_units = heating_powerlaw_.value_units;
  }
  ValidateDensityComposition(runtime_.cooling_density, cooling_density_units,
                             runtime_.composition, "cooling_density");
  ValidateDensityComposition(runtime_.heating_density, heating_density_units,
                             runtime_.composition, "heating_density");
  if (timestep_use_density_bounds_) {
    ValidateDensityComposition(timestep_density_kind_, timestep_density_units_,
                               runtime_.composition, "timestep_density_kind");
  }
  ValidateTableComposition(cooling_table_, runtime_.composition, "cooling table");
  ValidateTableComposition(heating_table_, runtime_.composition, "heating table");
  ValidateTableComposition(cgm_pie_table_, runtime_.composition, "cgm PIE table");
  ValidateTableComposition(cgm_cie_table_, runtime_.composition, "cgm CIE table");
  ValidateModifierComposition(cooling_modifier_, runtime_.composition,
                              "cooling modifier");
  ValidateModifierComposition(heating_modifier_, runtime_.composition,
                              "heating modifier");
  if (shielding_.enabled) {
    ValidateDensityComposition(shielding_.density_kind, UnitSystem::cgs,
                               runtime_.composition, "shielding_density");
  }
}

void GeneralCooling::Apply(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                           const Real bdt, const Real history_bdt,
                           DvceArray5D<Real> &u0) {
  if (!enabled_) return;
  if (!eos_data.is_ideal) {
    FatalCoolingInput("cooling source terms require an ideal-gas energy equation.");
  }
  if (cooling_model_ == ModelKind::user || heating_model_ == ModelKind::user) {
    if (pmy_pack_->pmesh->pgen == nullptr ||
        pmy_pack_->pmesh->pgen->user_cooling_func == nullptr) {
      FatalCoolingInput("cooling_model=user or heating_model=user requires "
                        "UserProblem() to set user_cooling_func.");
    }
    Real gross_energy = 0.0;
    Real net_energy = 0.0;
    pmy_pack_->pmesh->pgen->user_cooling_func(pmy_pack_, w0, eos_data, runtime_, bdt,
                                              history_bdt, u0, gross_energy,
                                              net_energy);
    if (history_enabled_) {
      accumulated_gross_energy_ += gross_energy;
      accumulated_net_energy_ += net_energy;
    }
    return;
  }
  auto &indcs = pmy_pack_->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack_->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  int nfluid = 0, nscalars = 0;
  if (pmy_pack_->phydro != nullptr) {
    nfluid = pmy_pack_->phydro->nhydro;
    nscalars = pmy_pack_->phydro->nscalars;
  } else if (pmy_pack_->pmhd != nullptr) {
    nfluid = pmy_pack_->pmhd->nmhd;
    nscalars = pmy_pack_->pmhd->nscalars;
  }

  const Real gm1 = eos_data.gamma - 1.0;
  const Real use_e = eos_data.use_e;
  const RuntimeData runtime = runtime_;
  const ModelKind cooling_model = cooling_model_;
  const ModelKind heating_model = heating_model_;
  const TableData cooling_table = cooling_table_;
  const TableData heating_table = heating_table_;
  const TableData cgm_pie_table = cgm_pie_table_;
  const TableData cgm_cie_table = cgm_cie_table_;
  const PowerLawData cooling_powerlaw = cooling_powerlaw_;
  const PowerLawData heating_powerlaw = heating_powerlaw_;
  const ModifierData cooling_modifier = cooling_modifier_;
  const ModifierData heating_modifier = heating_modifier_;
  const ShieldingData shielding = shielding_;
  const Real constant_heating_rate = constant_heating_rate_;
  const Real history_bdt_local = history_bdt;
  auto &size = pmy_pack_->pmb->mb_size;

  Real gross_energy = 0.0, net_energy = 0.0, bad_bounds = 0.0;
  Kokkos::parallel_reduce("general_cooling_source",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &gross_sum, Real &net_sum, Real &bad_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    CoolingCellState state;
    state.rho_code = w0(m,IDN,k,j,i);
    state.rho_cgs = state.rho_code*runtime.density_cgs;
    state.nfluid = nfluid;
    state.nscalars = nscalars;
    if (use_e) {
      state.temp_code = w0(m,IEN,k,j,i)/state.rho_code*gm1;
      state.eint_code = w0(m,IEN,k,j,i);
    } else {
      state.temp_code = w0(m,ITM,k,j,i);
      state.eint_code = w0(m,ITM,k,j,i)*state.rho_code/gm1;
    }
    state.temp_cgs = state.temp_code*runtime.temp_cgs;

    Real local_bad = 0.0;
    Real lambda = 0.0;
    UnitSystem lambda_units = runtime.default_units;
    Real shielding_frac = 1.0;
    if (shielding.enabled) {
      shielding_frac = EvaluateShieldingFraction(shielding, state, runtime,
                                                 size.d_view(m).dx1);
    }

    if (cooling_model == ModelKind::ism) {
      lambda = ISMCoolFn(state.temp_cgs);
      lambda_units = UnitSystem::cgs;
    } else if (cooling_model == ModelKind::cgm) {
      const Real lambda_pie = EvaluateTable(cgm_pie_table, state, w0, runtime,
                                            m, k, j, i, local_bad);
      const Real lambda_cie = EvaluateTable(cgm_cie_table, state, w0, runtime,
                                            m, k, j, i, local_bad);
      lambda = (1.0 - shielding_frac)*lambda_cie + shielding_frac*lambda_pie;
      lambda_units = cgm_pie_table.value_units;
    } else if (cooling_model == ModelKind::table) {
      lambda = EvaluateTable(cooling_table, state, w0, runtime, m, k, j, i, local_bad);
      lambda_units = cooling_table.value_units;
    } else if (cooling_model == ModelKind::powerlaw ||
               cooling_model == ModelKind::piecewise_powerlaw) {
      lambda = EvaluatePowerLaw(cooling_powerlaw, state, w0, runtime,
                                m, k, j, i, local_bad);
      lambda_units = cooling_powerlaw.value_units;
    }
    lambda *= EvaluateModifier(cooling_modifier, state, w0, runtime, m, k, j, i,
                               local_bad);

    Real gamma_heat = 0.0;
    UnitSystem gamma_units = runtime.default_units;
    if (heating_model == ModelKind::constant) {
      gamma_heat = constant_heating_rate;
      gamma_units = runtime.default_units;
    } else if (heating_model == ModelKind::table) {
      gamma_heat = EvaluateTable(heating_table, state, w0, runtime, m, k, j, i,
                                 local_bad);
      gamma_units = heating_table.value_units;
    } else if (heating_model == ModelKind::powerlaw ||
               heating_model == ModelKind::piecewise_powerlaw) {
      gamma_heat = EvaluatePowerLaw(heating_powerlaw, state, w0, runtime,
                                    m, k, j, i, local_bad);
      gamma_units = heating_powerlaw.value_units;
    }
    gamma_heat *= EvaluateModifier(heating_modifier, state, w0, runtime, m, k, j, i,
                                   local_bad);
    if (shielding.enabled && shielding.apply_to_heating) {
      gamma_heat *= (1.0 - shielding_frac);
    }

    const Real edot_cool = CoolingEdotCode(lambda, lambda_units,
                                           runtime.cooling_density,
                                           state.rho_code, runtime);
    const Real edot_heat = HeatingEdotCode(gamma_heat, gamma_units,
                                           runtime.heating_density,
                                           state.rho_code, runtime);
    const Real edot_net = edot_cool - edot_heat;
    const Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    u0(m,IEN,k,j,i) -= bdt*edot_net;
    gross_sum += edot_cool*history_bdt_local*dvol;
    net_sum += edot_net*history_bdt_local*dvol;
    bad_sum += local_bad;
  }, Kokkos::Sum<Real>(gross_energy), Kokkos::Sum<Real>(net_energy),
     Kokkos::Sum<Real>(bad_bounds));

  if (bad_bounds > 0.0) {
    FatalCoolingInput("fatal cooling/heating bounds were exceeded in at least one cell.");
  }
  if (history_enabled_) {
    accumulated_gross_energy_ += gross_energy;
    accumulated_net_energy_ += net_energy;
  }
}

void GeneralCooling::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data,
                                 Real &dtnew) {
  if (!TimestepEnabled()) return;
  if (!eos_data.is_ideal) {
    FatalCoolingInput("cooling timestep requires an ideal-gas energy equation.");
  }
  if (cooling_model_ == ModelKind::user || heating_model_ == ModelKind::user) {
    if (pmy_pack_->pmesh->pgen == nullptr ||
        pmy_pack_->pmesh->pgen->user_cooling_timestep_func == nullptr) {
      FatalCoolingInput("cooling timestep with cooling_model=user or heating_model=user "
                        "requires UserProblem() to set user_cooling_timestep_func.");
    }
    TimestepData timestep;
    timestep.factor = timestep_factor_;
    timestep.use_temperature_bounds = timestep_use_temperature_bounds_;
    timestep.temperature_units = timestep_temperature_units_;
    timestep.temperature_min = timestep_temperature_min_;
    timestep.temperature_max = timestep_temperature_max_;
    timestep.use_density_bounds = timestep_use_density_bounds_;
    timestep.density_units = timestep_density_units_;
    timestep.density_kind = timestep_density_kind_;
    timestep.density_min = timestep_density_min_;
    timestep.density_max = timestep_density_max_;
    pmy_pack_->pmesh->pgen->user_cooling_timestep_func(pmy_pack_, w0, eos_data,
                                                       runtime_, timestep,
                                                       dtnew);
    return;
  }
  auto &indcs = pmy_pack_->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack_->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  int nfluid = 0, nscalars = 0;
  if (pmy_pack_->phydro != nullptr) {
    nfluid = pmy_pack_->phydro->nhydro;
    nscalars = pmy_pack_->phydro->nscalars;
  } else if (pmy_pack_->pmhd != nullptr) {
    nfluid = pmy_pack_->pmhd->nmhd;
    nscalars = pmy_pack_->pmhd->nscalars;
  }

  const Real gm1 = eos_data.gamma - 1.0;
  const Real use_e = eos_data.use_e;
  const RuntimeData runtime = runtime_;
  const ModelKind cooling_model = cooling_model_;
  const ModelKind heating_model = heating_model_;
  const TableData cooling_table = cooling_table_;
  const TableData heating_table = heating_table_;
  const TableData cgm_pie_table = cgm_pie_table_;
  const TableData cgm_cie_table = cgm_cie_table_;
  const PowerLawData cooling_powerlaw = cooling_powerlaw_;
  const PowerLawData heating_powerlaw = heating_powerlaw_;
  const ModifierData cooling_modifier = cooling_modifier_;
  const ModifierData heating_modifier = heating_modifier_;
  const ShieldingData shielding = shielding_;
  const Real constant_heating_rate = constant_heating_rate_;
  const Real timestep_factor = timestep_factor_;
  const bool use_temp_bounds = timestep_use_temperature_bounds_;
  const UnitSystem temp_bounds_units = timestep_temperature_units_;
  const Real temp_min = timestep_temperature_min_;
  const Real temp_max = timestep_temperature_max_;
  const bool use_density_bounds = timestep_use_density_bounds_;
  const UnitSystem density_bounds_units = timestep_density_units_;
  const DensityKind density_bounds_kind = timestep_density_kind_;
  const Real density_min = timestep_density_min_;
  const Real density_max = timestep_density_max_;
  auto &size = pmy_pack_->pmb->mb_size;

  Real bad_bounds = 0.0;
  Kokkos::parallel_reduce("general_cooling_newdt",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt, Real &bad_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    CoolingCellState state;
    state.rho_code = w0(m,IDN,k,j,i);
    state.rho_cgs = state.rho_code*runtime.density_cgs;
    state.nfluid = nfluid;
    state.nscalars = nscalars;
    if (use_e) {
      state.temp_code = w0(m,IEN,k,j,i)/state.rho_code*gm1;
      state.eint_code = w0(m,IEN,k,j,i);
    } else {
      state.temp_code = w0(m,ITM,k,j,i);
      state.eint_code = w0(m,ITM,k,j,i)*state.rho_code/gm1;
    }
    state.temp_cgs = state.temp_code*runtime.temp_cgs;

    if (use_temp_bounds) {
      const Real dt_temp = (temp_bounds_units == UnitSystem::cgs) ?
                           state.temp_cgs : state.temp_code;
      if (dt_temp < temp_min || dt_temp > temp_max) return;
    }
    if (use_density_bounds) {
      const Real dt_density = DensityValue(density_bounds_kind, density_bounds_units,
                                           state.rho_code, runtime);
      if (dt_density < density_min || dt_density > density_max) return;
    }

    Real local_bad = 0.0;
    Real lambda = 0.0;
    UnitSystem lambda_units = runtime.default_units;
    Real shielding_frac = 1.0;
    if (shielding.enabled) {
      shielding_frac = EvaluateShieldingFraction(shielding, state, runtime,
                                                 size.d_view(m).dx1);
    }

    if (cooling_model == ModelKind::ism) {
      lambda = ISMCoolFn(state.temp_cgs);
      lambda_units = UnitSystem::cgs;
    } else if (cooling_model == ModelKind::cgm) {
      const Real lambda_pie = EvaluateTable(cgm_pie_table, state, w0, runtime,
                                            m, k, j, i, local_bad);
      const Real lambda_cie = EvaluateTable(cgm_cie_table, state, w0, runtime,
                                            m, k, j, i, local_bad);
      lambda = (1.0 - shielding_frac)*lambda_cie + shielding_frac*lambda_pie;
      lambda_units = cgm_pie_table.value_units;
    } else if (cooling_model == ModelKind::table) {
      lambda = EvaluateTable(cooling_table, state, w0, runtime, m, k, j, i, local_bad);
      lambda_units = cooling_table.value_units;
    } else if (cooling_model == ModelKind::powerlaw ||
               cooling_model == ModelKind::piecewise_powerlaw) {
      lambda = EvaluatePowerLaw(cooling_powerlaw, state, w0, runtime,
                                m, k, j, i, local_bad);
      lambda_units = cooling_powerlaw.value_units;
    }
    lambda *= EvaluateModifier(cooling_modifier, state, w0, runtime, m, k, j, i,
                               local_bad);
    Real gamma_heat = 0.0;
    UnitSystem gamma_units = runtime.default_units;
    if (heating_model == ModelKind::constant) {
      gamma_heat = constant_heating_rate;
      gamma_units = runtime.default_units;
    } else if (heating_model == ModelKind::table) {
      gamma_heat = EvaluateTable(heating_table, state, w0, runtime, m, k, j, i,
                                 local_bad);
      gamma_units = heating_table.value_units;
    } else if (heating_model == ModelKind::powerlaw ||
               heating_model == ModelKind::piecewise_powerlaw) {
      gamma_heat = EvaluatePowerLaw(heating_powerlaw, state, w0, runtime,
                                    m, k, j, i, local_bad);
      gamma_units = heating_powerlaw.value_units;
    }
    gamma_heat *= EvaluateModifier(heating_modifier, state, w0, runtime, m, k, j, i,
                                   local_bad);
    if (shielding.enabled && shielding.apply_to_heating) {
      gamma_heat *= (1.0 - shielding_frac);
    }
    const Real edot_cool = CoolingEdotCode(lambda, lambda_units,
                                           runtime.cooling_density,
                                           state.rho_code, runtime);
    const Real edot_heat = HeatingEdotCode(gamma_heat, gamma_units,
                                           runtime.heating_density,
                                           state.rho_code, runtime);
    const Real edot_net = edot_cool - edot_heat;
    const Real rate = fabs(edot_net) + FLT_MIN;
    min_dt = fmin(min_dt, timestep_factor*state.eint_code/rate);
    bad_sum += local_bad;
  }, Kokkos::Min<Real>(dtnew), Kokkos::Sum<Real>(bad_bounds));

  if (bad_bounds > 0.0) {
    FatalCoolingInput("fatal cooling/heating bounds were exceeded while computing dt.");
  }
}

int GeneralCooling::AddHistoryLabels(std::string *labels, int start, int max_labels) const {
  if (!HistoryEnabled()) return 0;
  int n = 0;
  if (history_gross_) {
    if (start + n >= max_labels) FatalCoolingInput("Not enough history columns.");
    labels[start + n++] = "cool_gross";
  }
  if (history_net_) {
    if (start + n >= max_labels) FatalCoolingInput("Not enough history columns.");
    labels[start + n++] = "cool_net";
  }
  return n;
}

int GeneralCooling::AddHistoryData(Real *hdata, int start, int max_data,
                                   Real current_time) {
  if (!HistoryEnabled()) return 0;
  Real elapsed = current_time - last_history_time_;
  const bool valid_elapsed = elapsed > 0.0;
  int n = 0;
  if (history_gross_) {
    if (start + n >= max_data) FatalCoolingInput("Not enough history data columns.");
    hdata[start + n++] = valid_elapsed ? accumulated_gross_energy_/elapsed : 0.0;
  }
  if (history_net_) {
    if (start + n >= max_data) FatalCoolingInput("Not enough history data columns.");
    hdata[start + n++] = valid_elapsed ? accumulated_net_energy_/elapsed : 0.0;
  }
  accumulated_gross_energy_ = 0.0;
  accumulated_net_energy_ = 0.0;
  last_history_time_ = current_time;
  return n;
}

} // namespace cooling
