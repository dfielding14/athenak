//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file tracer_fields.cpp
//! \brief parser and host-side evaluators for lagrangian_mc tracer fields

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "particles/tracer_fields.hpp"

namespace {

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

void FatalTracerField(const std::string &context, const std::string &msg) {
  std::cout << "### FATAL ERROR in tracer_fields.cpp" << std::endl
            << context << ": " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

bool StartsWith(const std::string &value, const std::string &prefix) {
  return value.compare(0, prefix.size(), prefix) == 0;
}

bool IsMHDOnly(particles::TracerFieldKind kind) {
  using particles::TracerFieldKind;
  return (kind == TracerFieldKind::b1 ||
          kind == TracerFieldKind::b2 ||
          kind == TracerFieldKind::b3 ||
          kind == TracerFieldKind::bmag ||
          kind == TracerFieldKind::magnetic_pressure ||
          kind == TracerFieldKind::beta ||
          kind == TracerFieldKind::alfven_speed);
}

void CheckAvailability(const particles::TracerField &field, bool has_mhd, int nscalars,
                       const std::string &context) {
  if (IsMHDOnly(field.kind) && !has_mhd) {
    FatalTracerField(context, "field '" + field.name + "' requires MHD");
  }
  if (field.kind == particles::TracerFieldKind::scalar &&
      (field.scalar_index < 0 || field.scalar_index >= nscalars)) {
    FatalTracerField(context, "field '" + field.name + "' is out of range for nscalars=" +
                              std::to_string(nscalars));
  }
}

void AppendUnique(std::vector<particles::TracerField> &fields,
                  const particles::TracerField &field) {
  for (const auto &existing : fields) {
    if (existing.name == field.name) return;
  }
  fields.push_back(field);
}

} // namespace

namespace particles {

TracerField ParseTracerFieldName(const std::string &field_name, bool has_mhd,
                                 int nscalars, const std::string &context) {
  const std::string name = Lower(field_name);
  TracerField field;

  if (name == "rho" || name == "density") {
    field = {TracerFieldKind::density, -1, "density"};
  } else if (name == "p" || name == "pressure") {
    field = {TracerFieldKind::pressure, -1, "pressure"};
  } else if (name == "t" || name == "temperature") {
    field = {TracerFieldKind::temperature, -1, "temperature"};
  } else if (name == "s" || name == "entropy" || name == "specific_entropy") {
    field = {TracerFieldKind::entropy, -1, "entropy"};
  } else if (name == "eint" || name == "internal_energy") {
    field = {TracerFieldKind::internal_energy, -1, "internal_energy"};
  } else if (name == "v1" || name == "vx") {
    field = {TracerFieldKind::v1, -1, "v1"};
  } else if (name == "v2" || name == "vy") {
    field = {TracerFieldKind::v2, -1, "v2"};
  } else if (name == "v3" || name == "vz") {
    field = {TracerFieldKind::v3, -1, "v3"};
  } else if (name == "vmag" || name == "velocity_magnitude") {
    field = {TracerFieldKind::vmag, -1, "vmag"};
  } else if (name == "sound_speed" || name == "cs") {
    field = {TracerFieldKind::sound_speed, -1, "sound_speed"};
  } else if (name == "mach" || name == "mach_number") {
    field = {TracerFieldKind::mach, -1, "mach"};
  } else if (name == "b1" || name == "bx") {
    field = {TracerFieldKind::b1, -1, "b1"};
  } else if (name == "b2" || name == "by") {
    field = {TracerFieldKind::b2, -1, "b2"};
  } else if (name == "b3" || name == "bz") {
    field = {TracerFieldKind::b3, -1, "b3"};
  } else if (name == "bmag" || name == "magnetic_field_magnitude") {
    field = {TracerFieldKind::bmag, -1, "bmag"};
  } else if (name == "magnetic_pressure" || name == "pmag") {
    field = {TracerFieldKind::magnetic_pressure, -1, "magnetic_pressure"};
  } else if (name == "beta" || name == "plasma_beta") {
    field = {TracerFieldKind::beta, -1, "beta"};
  } else if (name == "alfven_speed" || name == "va") {
    field = {TracerFieldKind::alfven_speed, -1, "alfven_speed"};
  } else if (StartsWith(name, "scalar")) {
    const std::string suffix = name.substr(6);
    if (suffix.empty()) {
      FatalTracerField(context, "scalar fields must be named scalarN");
    }
    try {
      field = {TracerFieldKind::scalar, std::stoi(suffix), name};
    } catch (...) {
      FatalTracerField(context, "scalar fields must be named scalarN");
    }
  } else {
    FatalTracerField(context, "unknown field '" + field_name + "'");
  }

  CheckAvailability(field, has_mhd, nscalars, context);
  return field;
}

std::vector<TracerField> ParseTracerFieldList(const std::string &field_list,
                                              bool has_mhd, int nscalars,
                                              const std::string &context) {
  std::string normalized = Lower(field_list);
  std::replace(normalized.begin(), normalized.end(), ',', ' ');
  std::stringstream stream(normalized);
  std::string token;
  std::vector<TracerField> fields;

  while (stream >> token) {
    if (token == "default") {
      AppendUnique(fields, ParseTracerFieldName("density", has_mhd, nscalars, context));
      AppendUnique(fields, ParseTracerFieldName("pressure", has_mhd, nscalars, context));
      AppendUnique(fields, ParseTracerFieldName("temperature", has_mhd, nscalars,
                                                context));
      AppendUnique(fields, ParseTracerFieldName("entropy", has_mhd, nscalars, context));
      AppendUnique(fields, ParseTracerFieldName("internal_energy", has_mhd, nscalars,
                                                context));
      AppendUnique(fields, ParseTracerFieldName("v1", has_mhd, nscalars, context));
      AppendUnique(fields, ParseTracerFieldName("v2", has_mhd, nscalars, context));
      AppendUnique(fields, ParseTracerFieldName("v3", has_mhd, nscalars, context));
    } else {
      AppendUnique(fields, ParseTracerFieldName(token, has_mhd, nscalars, context));
    }
  }

  if (fields.empty()) {
    FatalTracerField(context, "field list is empty");
  }
  return fields;
}

std::vector<std::string> TracerFieldNames(const std::vector<TracerField> &fields) {
  std::vector<std::string> names;
  names.reserve(fields.size());
  for (const auto &field : fields) names.push_back(field.name);
  return names;
}

Real EvaluateTracerFieldHost(const TracerField &field, const HostArray5D<Real> &w0,
                             const HostArray5D<Real> &bcc, bool has_mhd,
                             Real gamma, int nfluid, int m, int k, int j, int i) {
  const Real rho = w0(m,IDN,k,j,i);
  const Real eint = w0(m,IEN,k,j,i);
  const Real press = (gamma - 1.0)*eint;
  const Real v1 = w0(m,IVX,k,j,i);
  const Real v2 = w0(m,IVY,k,j,i);
  const Real v3 = w0(m,IVZ,k,j,i);
  const Real vsq = SQR(v1) + SQR(v2) + SQR(v3);
  const Real cs2 = std::max(gamma*press/rho, 0.0);
  const Real cs = std::sqrt(cs2);

  switch (field.kind) {
    case TracerFieldKind::density:
      return rho;
    case TracerFieldKind::pressure:
      return press;
    case TracerFieldKind::temperature:
      return press/rho;
    case TracerFieldKind::entropy:
      return std::log(press/std::pow(rho, gamma));
    case TracerFieldKind::internal_energy:
      return eint;
    case TracerFieldKind::v1:
      return v1;
    case TracerFieldKind::v2:
      return v2;
    case TracerFieldKind::v3:
      return v3;
    case TracerFieldKind::vmag:
      return std::sqrt(vsq);
    case TracerFieldKind::sound_speed:
      return cs;
    case TracerFieldKind::mach:
      return (cs > 0.0) ? std::sqrt(vsq)/cs : 0.0;
    case TracerFieldKind::scalar:
      return w0(m,nfluid+field.scalar_index,k,j,i);
    case TracerFieldKind::b1:
      return has_mhd ? bcc(m,IBX,k,j,i) : 0.0;
    case TracerFieldKind::b2:
      return has_mhd ? bcc(m,IBY,k,j,i) : 0.0;
    case TracerFieldKind::b3:
      return has_mhd ? bcc(m,IBZ,k,j,i) : 0.0;
    case TracerFieldKind::bmag:
    case TracerFieldKind::magnetic_pressure:
    case TracerFieldKind::beta:
    case TracerFieldKind::alfven_speed:
      {
        const Real b1 = bcc(m,IBX,k,j,i);
        const Real b2 = bcc(m,IBY,k,j,i);
        const Real b3 = bcc(m,IBZ,k,j,i);
        const Real b2tot = SQR(b1) + SQR(b2) + SQR(b3);
        const Real bmag = std::sqrt(b2tot);
        const Real pmag = 0.5*b2tot;
        if (field.kind == TracerFieldKind::bmag) return bmag;
        if (field.kind == TracerFieldKind::magnetic_pressure) return pmag;
        if (field.kind == TracerFieldKind::beta) {
          return (pmag > 0.0) ? press/pmag : std::numeric_limits<Real>::max();
        }
        return (rho > 0.0) ? bmag/std::sqrt(rho) : 0.0;
      }
  }
  return 0.0;
}

} // namespace particles
