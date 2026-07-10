#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cfloat>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_set>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// The standalone reducer uses the exact CGM cooling tables shipped with the
// GOTHAM AthenaK branch. AthenaK spells the table scalar type as Real; this
// one-off reducer uses double throughout its derived state.
using Real = double;
#include "../../src/srcterms/cooling_tables.hpp"

namespace fs = std::filesystem;

namespace {

constexpr int kMaxDims = 4;
constexpr int kExpectedVars = 8;
constexpr int kRecordMetadataBytes = 10 * sizeof(std::int32_t) + 6 * sizeof(double);
constexpr double kCoolingHydrogenMassFraction = 0.75;
constexpr double kCoolingSolarMetallicity = 0.02;
constexpr double kAtomicMassUnitCgs = 1.67262192369e-24;
constexpr double kMyrCgs = 3.15576e13;

enum class Scale : int { linear = 0, log = 1, symlog = 2 };

enum class Variable : int {
  radius,
  costheta,
  abs_costheta,
  density,
  temperature,
  internal_energy,
  velocity_r,
  velocity_theta,
  velocity_phi,
  scalar_0,
  scalar_1,
  scalar_2,
  speed,
  mach,
  radial_mach,
  entropy_proxy,
  angular_momentum_z,
  angular_momentum_z_kpc_km_s,
  temperature_kelvin,
  hydrogen_number_density,
  pressure_over_kb,
  velocity_r_km_s,
  velocity_theta_km_s,
  velocity_phi_km_s,
  absolute_velocity_r_km_s,
  sound_speed_km_s,
  absolute_radial_mach,
  absolute_z,
  cylindrical_radius,
  dust_to_total_metal,
  small_grain_fraction,
  cooling_rate_erg_s_cm3,
  cooling_time_myr,
};

enum class Weight : int {
  volume,
  mass,
  mdot,
  mdot_out,
  mdot_in,
  edot,
  edot_out,
  edot_in,
  edot_kin,
  edot_th,
  mdot_in_abs,
  edot_in_abs,
  edot_kin_out,
  edot_kin_in_abs,
  edot_th_out,
  edot_th_in_abs,
  edot_cool,
  ram_out,
  vertical_mdot_out,
  vertical_mdot_in_abs,
  vertical_edot_out,
  vertical_ram_out,
  metal_mdot_out,
  metal_mdot_in_abs,
  total_metal_mdot_out,
  total_metal_mdot_in_abs,
  dust_mdot_out,
  dust_mdot_in_abs,
  total_metal_mass,
  dust_mass,
};

struct DeviceAxis {
  int variable = 0;
  int scale = 0;
  int nbin = 0;
  int stride = 0;
  double minimum = 0.0;
  double maximum = 0.0;
  double transformed_minimum = 0.0;
  double inverse_step = 0.0;
  double linthresh = 1.0;
};

struct DeviceProduct {
  std::uint64_t offset = 0;
  int total_bins = 0;
  int ndim = 0;
  int weight = 0;
  DeviceAxis axes[kMaxDims];
};

struct Axis {
  Variable variable;
  Scale scale;
  int nbin;
  double minimum;
  double maximum;
  double linthresh = 1.0;
};

struct Product {
  std::string id;
  std::vector<Axis> axes;
  Weight weight;
  DeviceProduct device;
};

struct FieldMap {
  int density = -1;
  int velocity_x = -1;
  int velocity_y = -1;
  int velocity_z = -1;
  int internal_energy = -1;
  int scalar_0 = -1;
  int scalar_1 = -1;
  int scalar_2 = -1;
};

struct HeaderInfo {
  double time = 0.0;
  double gamma = 0.0;
  double velocity_km_s = 0.0;
  double temperature_kelvin_per_ratio = 0.0;
  double hydrogen_number_density_per_density = 0.0;
  double pressure_over_kb_per_internal_energy = 0.0;
  double length_cgs = 0.0;
  double time_cgs = 0.0;
  double pressure_cgs = 0.0;
  double hrate = 0.0;
  double hscale_norm = 0.0;
  double hscale_height = 0.0;
  double hscale_radius = 0.0;
  int cgm_cooling = 0;
  double domain_volume = 0.0;
  std::array<double, 3> domain_min{};
  std::array<double, 3> domain_max{};
  std::array<int, 3> root_blocks{};
  std::uint64_t cycle = 0;
  std::uint64_t header_digest = 0;
  std::uint64_t data_offset = 0;
  std::uint64_t record_size = 0;
  std::uint64_t blocks = 0;
  int location_size = 0;
  int variable_size = 0;
  int nvars = 0;
  int nx1 = 0;
  int nx2 = 0;
  int nx3 = 0;
  FieldMap fields;
  std::vector<std::string> variables;
};

struct HeaderWire {
  double time;
  double gamma;
  double velocity_km_s;
  double temperature_kelvin_per_ratio;
  double hydrogen_number_density_per_density;
  double pressure_over_kb_per_internal_energy;
  double length_cgs;
  double time_cgs;
  double pressure_cgs;
  double hrate;
  double hscale_norm;
  double hscale_height;
  double hscale_radius;
  int cgm_cooling;
  double domain_volume;
  std::array<double, 3> domain_min;
  std::array<double, 3> domain_max;
  std::array<int, 3> root_blocks;
  std::uint64_t cycle;
  std::uint64_t header_digest;
  std::uint64_t record_size;
  int location_size;
  int variable_size;
  int nvars;
  int nx1;
  int nx2;
  int nx3;
  FieldMap fields;
};

struct Options {
  fs::path input_dir;
  fs::path output_dir;
  std::string sequence;
  std::string output_number;
  std::string basename = "gotham";
  std::string products = "all";
  std::vector<std::string> only_products;
  int chunk_blocks = 32;
  std::uint64_t max_shards = 0;
  std::uint64_t expected_shards = 0;
  std::uint64_t max_blocks_per_shard = 0;
  double output_time = std::numeric_limits<double>::quiet_NaN();
  bool dry_run = false;
  bool list_products = false;
  bool help = false;
};

struct BlockKey {
  std::int32_t lx1;
  std::int32_t lx2;
  std::int32_t lx3;
  std::int32_t level;

  bool operator==(const BlockKey &other) const {
    return lx1 == other.lx1 && lx2 == other.lx2 && lx3 == other.lx3 &&
           level == other.level;
  }
};

struct BlockKeyHash {
  std::size_t operator()(const BlockKey &key) const {
    std::uint64_t hash = 1469598103934665603ULL;
    const auto *bytes = reinterpret_cast<const unsigned char *>(&key);
    for (std::size_t i = 0; i < sizeof(key); ++i) {
      hash ^= bytes[i];
      hash *= 1099511628211ULL;
    }
    return static_cast<std::size_t>(hash);
  }
};

std::string trim(std::string value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string value_after_equals(const std::string &line, const std::string &context) {
  const auto pos = line.find('=');
  if (pos == std::string::npos) {
    throw std::runtime_error("Missing '=' while reading " + context + ": " + line);
  }
  return trim(line.substr(pos + 1));
}

std::string require_line(std::ifstream &stream, const std::string &context) {
  std::string line;
  if (!std::getline(stream, line)) {
    throw std::runtime_error("Unexpected EOF while reading " + context);
  }
  return line;
}

std::uint64_t fnv1a(const std::string &value) {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const unsigned char byte : value) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  return hash;
}

std::map<std::string, std::map<std::string, std::string>>
parse_input_deck(const std::string &text) {
  std::map<std::string, std::map<std::string, std::string>> result;
  std::string section;
  std::istringstream lines(text);
  std::string raw;
  while (std::getline(lines, raw)) {
    const auto comment = raw.find('#');
    if (comment != std::string::npos) raw.resize(comment);
    const std::string line = trim(raw);
    if (line.empty()) continue;
    if (line.front() == '<' && line.back() == '>') {
      section = trim(line.substr(1, line.size() - 2));
      continue;
    }
    const auto equal = line.find('=');
    if (equal != std::string::npos && !section.empty()) {
      result[section][trim(line.substr(0, equal))] = trim(line.substr(equal + 1));
    }
  }
  return result;
}

std::string input_value(
    const std::map<std::string, std::map<std::string, std::string>> &deck,
    const std::string &section, const std::string &key) {
  const auto section_it = deck.find(section);
  if (section_it == deck.end()) {
    throw std::runtime_error("Embedded input is missing <" + section + ">");
  }
  const auto key_it = section_it->second.find(key);
  if (key_it == section_it->second.end()) {
    throw std::runtime_error("Embedded input is missing <" + section + ">/" + key);
  }
  return key_it->second;
}

std::string input_value_or(
    const std::map<std::string, std::map<std::string, std::string>> &deck,
    const std::string &section, const std::string &key, const std::string &fallback) {
  const auto section_it = deck.find(section);
  if (section_it == deck.end()) return fallback;
  const auto key_it = section_it->second.find(key);
  return key_it == section_it->second.end() ? fallback : key_it->second;
}

bool input_bool_or(const std::map<std::string, std::map<std::string, std::string>> &deck,
                   const std::string &section, const std::string &key, bool fallback) {
  std::string value = input_value_or(deck, section, key, fallback ? "true" : "false");
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (value == "true" || value == "1" || value == "yes") return true;
  if (value == "false" || value == "0" || value == "no") return false;
  throw std::runtime_error("Embedded input has invalid boolean <" + section + ">/" + key +
                           ": " + value);
}

FieldMap make_field_map(const std::vector<std::string> &variables) {
  FieldMap fields;
  for (int index = 0; index < static_cast<int>(variables.size()); ++index) {
    const auto &name = variables[index];
    if (name == "dens") fields.density = index;
    if (name == "velx") fields.velocity_x = index;
    if (name == "vely") fields.velocity_y = index;
    if (name == "velz") fields.velocity_z = index;
    if (name == "eint") fields.internal_energy = index;
    if (name == "s_00") fields.scalar_0 = index;
    if (name == "s_01") fields.scalar_1 = index;
    if (name == "s_02") fields.scalar_2 = index;
  }
  if (fields.density < 0 || fields.velocity_x < 0 || fields.velocity_y < 0 ||
      fields.velocity_z < 0 || fields.internal_energy < 0 || fields.scalar_0 < 0 ||
      fields.scalar_1 < 0 || fields.scalar_2 < 0) {
    throw std::runtime_error(
        "The shard does not contain the exact GOTHAM hydro_w fields "
        "(dens velx vely velz eint s_00 s_01 s_02)");
  }
  return fields;
}

HeaderInfo read_header(const fs::path &path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("Unable to open shard " + path.string());

  HeaderInfo info;
  const std::string magic = require_line(stream, "Athena binary magic");
  if (magic.find("Athena binary output version=1.1") != 0) {
    throw std::runtime_error("Unsupported Athena binary format in " + path.string());
  }

  const int preheader_size =
      std::stoi(value_after_equals(require_line(stream, "preheader size"), "preheader size"));
  std::map<std::string, std::string> preheader;
  for (int i = 0; i < preheader_size - 1; ++i) {
    const std::string line = require_line(stream, "preheader");
    const auto equal = line.find('=');
    if (equal == std::string::npos) {
      throw std::runtime_error("Malformed preheader line in " + path.string());
    }
    preheader[trim(line.substr(0, equal))] = trim(line.substr(equal + 1));
  }

  info.time = std::stod(preheader.at("time"));
  info.cycle = std::stoull(preheader.at("cycle"));
  info.location_size = std::stoi(preheader.at("size of location"));
  info.variable_size = std::stoi(preheader.at("size of variable"));
  info.nvars = std::stoi(
      value_after_equals(require_line(stream, "number of variables"), "number of variables"));

  const std::string variables_line = require_line(stream, "variables");
  const auto colon = variables_line.find(':');
  if (colon == std::string::npos) {
    throw std::runtime_error("Malformed variables line in " + path.string());
  }
  std::istringstream variable_stream(variables_line.substr(colon + 1));
  std::string variable;
  while (variable_stream >> variable) info.variables.push_back(variable);
  if (static_cast<int>(info.variables.size()) != info.nvars) {
    throw std::runtime_error("Variable count mismatch in " + path.string());
  }
  info.fields = make_field_map(info.variables);

  const std::uint64_t embedded_size = std::stoull(
      value_after_equals(require_line(stream, "header offset"), "header offset"));
  std::string embedded(embedded_size, '\0');
  stream.read(embedded.data(), static_cast<std::streamsize>(embedded_size));
  if (static_cast<std::uint64_t>(stream.gcount()) != embedded_size) {
    throw std::runtime_error("Truncated embedded input in " + path.string());
  }
  info.data_offset = static_cast<std::uint64_t>(stream.tellg());
  info.header_digest = fnv1a(variables_line + embedded);

  const auto deck = parse_input_deck(embedded);
  info.nx1 = std::stoi(input_value(deck, "meshblock", "nx1"));
  info.nx2 = std::stoi(input_value(deck, "meshblock", "nx2"));
  info.nx3 = std::stoi(input_value(deck, "meshblock", "nx3"));
  info.gamma = std::stod(input_value(deck, "hydro", "gamma"));
  const double x1min = std::stod(input_value(deck, "mesh", "x1min"));
  const double x1max = std::stod(input_value(deck, "mesh", "x1max"));
  const double x2min = std::stod(input_value(deck, "mesh", "x2min"));
  const double x2max = std::stod(input_value(deck, "mesh", "x2max"));
  const double x3min = std::stod(input_value(deck, "mesh", "x3min"));
  const double x3max = std::stod(input_value(deck, "mesh", "x3max"));
  const int mesh_nx1 = std::stoi(input_value(deck, "mesh", "nx1"));
  const int mesh_nx2 = std::stoi(input_value(deck, "mesh", "nx2"));
  const int mesh_nx3 = std::stoi(input_value(deck, "mesh", "nx3"));
  if (mesh_nx1 <= 0 || mesh_nx2 <= 0 || mesh_nx3 <= 0 ||
      mesh_nx1 % info.nx1 != 0 || mesh_nx2 % info.nx2 != 0 ||
      mesh_nx3 % info.nx3 != 0) {
    throw std::runtime_error("Embedded root mesh is not divisible by MeshBlocks in " +
                             path.string());
  }
  info.domain_min = {x1min, x2min, x3min};
  info.domain_max = {x1max, x2max, x3max};
  info.root_blocks = {mesh_nx1 / info.nx1, mesh_nx2 / info.nx2,
                      mesh_nx3 / info.nx3};
  info.domain_volume = (x1max - x1min) * (x2max - x2min) * (x3max - x3min);
  const double length_cgs = std::stod(input_value(deck, "units", "length_cgs"));
  const double mass_cgs = std::stod(input_value(deck, "units", "mass_cgs"));
  const double time_cgs = std::stod(input_value(deck, "units", "time_cgs"));
  const double mu = std::stod(input_value(deck, "units", "mu"));
  constexpr double kBoltzmannCgs = 1.380649e-16;
  constexpr double protonMassCgs = 1.67262192369e-24;
  const double density_cgs = mass_cgs / (length_cgs * length_cgs * length_cgs);
  const double pressure_cgs = mass_cgs / (length_cgs * time_cgs * time_cgs);
  info.velocity_km_s = length_cgs / time_cgs / 1.0e5;
  info.temperature_kelvin_per_ratio =
      (info.gamma - 1.0) * pressure_cgs / density_cgs * mu * protonMassCgs /
      kBoltzmannCgs;
  // The CGM source term uses n_H = X rho / m_p with X=0.75.
  info.hydrogen_number_density_per_density =
      kCoolingHydrogenMassFraction * density_cgs / kAtomicMassUnitCgs;
  info.pressure_over_kb_per_internal_energy =
      (info.gamma - 1.0) * pressure_cgs / kBoltzmannCgs;
  info.length_cgs = length_cgs;
  info.time_cgs = time_cgs;
  info.pressure_cgs = pressure_cgs;
  info.cgm_cooling = input_bool_or(deck, "hydro", "cgm_cooling", false) ? 1 : 0;
  info.hrate = std::stod(input_value_or(deck, "hydro", "hrate", "0.0"));
  info.hscale_norm = std::stod(input_value_or(deck, "hydro", "hscale_norm", "0.0"));
  info.hscale_height = std::stod(input_value_or(deck, "hydro", "hscale_height", "0.0"));
  info.hscale_radius = std::stod(input_value_or(deck, "hydro", "hscale_radius", "0.0"));

  const std::uint64_t cells =
      static_cast<std::uint64_t>(info.nx1) * info.nx2 * info.nx3;
  info.record_size = 10 * sizeof(std::int32_t) + 6 * info.location_size +
                     cells * info.nvars * info.variable_size;
  const std::uint64_t file_size = fs::file_size(path);
  if (file_size < info.data_offset ||
      (file_size - info.data_offset) % info.record_size != 0) {
    throw std::runtime_error("Shard data size is inconsistent with fixed records: " +
                             path.string());
  }
  info.blocks = (file_size - info.data_offset) / info.record_size;
  return info;
}

HeaderWire to_wire(const HeaderInfo &info) {
  return HeaderWire{info.time,
                    info.gamma,
                    info.velocity_km_s,
                    info.temperature_kelvin_per_ratio,
                    info.hydrogen_number_density_per_density,
                    info.pressure_over_kb_per_internal_energy,
                    info.length_cgs,
                    info.time_cgs,
                    info.pressure_cgs,
                    info.hrate,
                    info.hscale_norm,
                    info.hscale_height,
                    info.hscale_radius,
                    info.cgm_cooling,
                    info.domain_volume,
                    info.domain_min,
                    info.domain_max,
                    info.root_blocks,
                    info.cycle,
                    info.header_digest,
                    info.record_size,
                    info.location_size,
                    info.variable_size,
                    info.nvars,
                    info.nx1,
                    info.nx2,
                    info.nx3,
                    info.fields};
}

void validate_header(const HeaderInfo &info, const HeaderWire &reference,
                     const fs::path &path) {
  const auto fail = [&path](const std::string &what) {
    throw std::runtime_error("Shard metadata mismatch for " + what + ": " + path.string());
  };
  if (info.time != reference.time) fail("time");
  if (info.gamma != reference.gamma) fail("gamma");
  if (info.velocity_km_s != reference.velocity_km_s) fail("velocity unit");
  if (info.temperature_kelvin_per_ratio != reference.temperature_kelvin_per_ratio) {
    fail("temperature unit");
  }
  if (info.hydrogen_number_density_per_density !=
      reference.hydrogen_number_density_per_density) {
    fail("hydrogen number density unit");
  }
  if (info.pressure_over_kb_per_internal_energy !=
      reference.pressure_over_kb_per_internal_energy) {
    fail("pressure unit");
  }
  if (info.length_cgs != reference.length_cgs) fail("length unit");
  if (info.time_cgs != reference.time_cgs) fail("time unit");
  if (info.pressure_cgs != reference.pressure_cgs) fail("pressure unit");
  if (info.hrate != reference.hrate) fail("CGM heating rate");
  if (info.hscale_norm != reference.hscale_norm) fail("CGM heating normalization");
  if (info.hscale_height != reference.hscale_height) fail("CGM heating scale height");
  if (info.hscale_radius != reference.hscale_radius) fail("CGM heating scale radius");
  if (info.cgm_cooling != reference.cgm_cooling) fail("CGM cooling enablement");
  if (info.domain_volume != reference.domain_volume) fail("domain volume");
  if (info.domain_min != reference.domain_min || info.domain_max != reference.domain_max) {
    fail("domain bounds");
  }
  if (info.root_blocks != reference.root_blocks) fail("root MeshBlock counts");
  if (info.cycle != reference.cycle) fail("cycle");
  if (info.header_digest != reference.header_digest) fail("complete header digest");
  if (info.record_size != reference.record_size) fail("record size");
  if (info.location_size != reference.location_size) fail("location size");
  if (info.variable_size != reference.variable_size) fail("variable size");
  if (info.nvars != reference.nvars) fail("number of variables");
  if (info.nx1 != reference.nx1 || info.nx2 != reference.nx2 ||
      info.nx3 != reference.nx3) {
    fail("MeshBlock shape");
  }
  if (std::memcmp(&info.fields, &reference.fields, sizeof(FieldMap)) != 0) {
    fail("field ordering");
  }
}

void validate_gotham_contract(const HeaderWire &header) {
  if (header.location_size != static_cast<int>(sizeof(double)) ||
      header.variable_size != static_cast<int>(sizeof(float)) ||
      header.nvars != kExpectedVars || header.nx1 <= 0 || header.nx2 <= 0 ||
      header.nx3 <= 0 || header.root_blocks[0] <= 0 || header.root_blocks[1] <= 0 ||
      header.root_blocks[2] <= 0 || !std::isfinite(header.time_cgs) ||
      !(header.time_cgs > 0.0) ||
      header.record_size !=
          static_cast<std::uint64_t>(kRecordMetadataBytes +
                                     kExpectedVars * header.nx1 * header.nx2 * header.nx3 *
                                         sizeof(float))) {
    throw std::runtime_error(
        "This one-off executable only accepts the five GOTHAM full-cube layouts: "
        "double geometry, float fields, and 8 hydro_w variables");
  }
}

std::string variable_name(Variable variable) {
  switch (variable) {
    case Variable::radius: return "coord_r";
    case Variable::costheta: return "coord_costheta";
    case Variable::abs_costheta: return "coord_abscostheta";
    case Variable::density: return "hydro_w_d";
    case Variable::temperature: return "temperature";
    case Variable::internal_energy: return "hydro_w_e";
    case Variable::velocity_r: return "vel_sph_r";
    case Variable::velocity_theta: return "vel_sph_theta";
    case Variable::velocity_phi: return "vel_sph_phi";
    case Variable::scalar_0: return "hydro_w_s_00";
    case Variable::scalar_1: return "hydro_w_s_01";
    case Variable::scalar_2: return "hydro_w_s_02";
    case Variable::speed: return "speed";
    case Variable::mach: return "mach";
    case Variable::radial_mach: return "radial_mach";
    case Variable::entropy_proxy: return "entropy_proxy";
    case Variable::angular_momentum_z: return "specific_angular_momentum_z";
    case Variable::angular_momentum_z_kpc_km_s:
      return "specific_angular_momentum_z_kpc_km_s";
    case Variable::temperature_kelvin: return "temperature_kelvin";
    case Variable::hydrogen_number_density: return "hydrogen_number_density_cm3";
    case Variable::pressure_over_kb: return "pressure_over_kb_k_cm3";
    case Variable::velocity_r_km_s: return "velocity_r_km_s";
    case Variable::velocity_theta_km_s: return "velocity_theta_km_s";
    case Variable::velocity_phi_km_s: return "velocity_phi_km_s";
    case Variable::absolute_velocity_r_km_s: return "absolute_velocity_r_km_s";
    case Variable::sound_speed_km_s: return "sound_speed_km_s";
    case Variable::absolute_radial_mach: return "absolute_radial_mach";
    case Variable::absolute_z: return "absolute_z";
    case Variable::cylindrical_radius: return "cylindrical_radius";
    case Variable::dust_to_total_metal: return "dust_to_total_metal";
    case Variable::small_grain_fraction: return "small_grain_fraction";
    case Variable::cooling_rate_erg_s_cm3: return "cooling_rate_erg_s_cm3";
    case Variable::cooling_time_myr: return "cooling_time_myr";
  }
  throw std::runtime_error("Unknown variable");
}

std::string scale_name(Scale scale) {
  switch (scale) {
    case Scale::linear: return "linear";
    case Scale::log: return "log";
    case Scale::symlog: return "symlog";
  }
  throw std::runtime_error("Unknown scale");
}

std::string weight_name(Weight weight) {
  switch (weight) {
    case Weight::volume: return "volume";
    case Weight::mass: return "mass";
    case Weight::mdot: return "mdot_sph";
    case Weight::mdot_out: return "mdot_sph_out";
    case Weight::mdot_in: return "mdot_sph_in";
    case Weight::edot: return "edot_sph";
    case Weight::edot_out: return "edot_sph_out";
    case Weight::edot_in: return "edot_sph_in";
    case Weight::edot_kin: return "edot_sph_kin";
    case Weight::edot_th: return "edot_sph_th";
    case Weight::mdot_in_abs: return "mdot_sph_in_abs";
    case Weight::edot_in_abs: return "edot_sph_in_abs";
    case Weight::edot_kin_out: return "edot_sph_kin_out";
    case Weight::edot_kin_in_abs: return "edot_sph_kin_in_abs";
    case Weight::edot_th_out: return "edot_sph_th_out";
    case Weight::edot_th_in_abs: return "edot_sph_th_in_abs";
    case Weight::edot_cool: return "edot_cool";
    case Weight::ram_out: return "radial_ram_out";
    case Weight::vertical_mdot_out: return "vertical_mdot_out";
    case Weight::vertical_mdot_in_abs: return "vertical_mdot_in_abs";
    case Weight::vertical_edot_out: return "vertical_edot_out";
    case Weight::vertical_ram_out: return "vertical_ram_out";
    case Weight::metal_mdot_out: return "gas_metal_mdot_out";
    case Weight::metal_mdot_in_abs: return "gas_metal_mdot_in_abs";
    case Weight::total_metal_mdot_out: return "total_metal_mdot_out";
    case Weight::total_metal_mdot_in_abs: return "total_metal_mdot_in_abs";
    case Weight::dust_mdot_out: return "dust_mdot_out";
    case Weight::dust_mdot_in_abs: return "dust_mdot_in_abs";
    case Weight::total_metal_mass: return "total_metal_mass";
    case Weight::dust_mass: return "dust_mass";
  }
  throw std::runtime_error("Unknown weight");
}

double symlog_forward_host(double value, double linthresh) {
  const double absolute = std::abs(value);
  const double transformed =
      absolute <= linthresh ? absolute / linthresh
                            : 1.0 + std::log10(absolute / linthresh);
  return std::copysign(transformed, value);
}

double symlog_inverse_host(double value, double linthresh) {
  const double absolute = std::abs(value);
  const double physical =
      absolute <= 1.0 ? absolute * linthresh
                      : linthresh * std::pow(10.0, absolute - 1.0);
  return std::copysign(physical, value);
}

double transform_host(double value, Scale scale, double linthresh) {
  if (scale == Scale::log) return std::log10(value);
  if (scale == Scale::symlog) return symlog_forward_host(value, linthresh);
  return value;
}

Axis linear(Variable variable, int nbin, double minimum, double maximum) {
  return Axis{variable, Scale::linear, nbin, minimum, maximum, 1.0};
}

Axis log_axis(Variable variable, int nbin, double minimum, double maximum) {
  return Axis{variable, Scale::log, nbin, minimum, maximum, 1.0};
}

Axis symlog_axis(Variable variable, int nbin, double minimum, double maximum,
                 double linthresh) {
  return Axis{variable, Scale::symlog, nbin, minimum, maximum, linthresh};
}

void append_product(std::vector<Product> &products, const std::string &id,
                    std::initializer_list<Axis> axes, Weight weight) {
  Product product;
  product.id = id;
  product.axes.assign(axes.begin(), axes.end());
  product.weight = weight;
  products.push_back(std::move(product));
}

std::vector<Product> make_products(
    const std::string &selection,
    const std::vector<std::string> &only_products = {}) {
  if (selection != "original" && selection != "science" && selection != "all") {
    throw std::runtime_error("--products must be original, science, or all");
  }
  const bool original = selection == "original" || selection == "all";
  const bool science = selection == "science" || selection == "all";
  std::vector<Product> products;

  const Axis r32_wide = log_axis(Variable::radius, 32, 1.0e-1, 190.0);
  const Axis r128_wide = log_axis(Variable::radius, 128, 1.0e-1, 190.0);
  const Axis r32_core = log_axis(Variable::radius, 32, 1.0e-1, 20.0);
  const Axis r4_core = log_axis(Variable::radius, 4, 1.0e-1, 20.0);
  const Axis absct16 = linear(Variable::abs_costheta, 16, 0.0, 1.0);
  const Axis absct4 = linear(Variable::abs_costheta, 4, 0.0, 1.0);
  const Axis cost16 = linear(Variable::costheta, 16, 0.0, 1.0);
  const Axis temp128 = log_axis(Variable::temperature, 128, 1.0e-3, 7.0e4);
  const Axis rho128 = log_axis(Variable::density, 128, 1.0e-5, 1.0e5);
  const Axis eint128 = log_axis(Variable::internal_energy, 128, 5.0e-5, 5.0e3);
  const Axis vr256 = linear(Variable::velocity_r, 256, -100.0, 500.0);
  const Axis vt128 = linear(Variable::velocity_theta, 128, -25.0, 25.0);
  const Axis vp128 = linear(Variable::velocity_phi, 128, -25.0, 25.0);
  const Axis vr128_positive = log_axis(Variable::velocity_r, 128, 1.0e-1, 1.0e3);
  const Axis vr64_positive = log_axis(Variable::velocity_r, 64, 1.0e-1, 1.0e3);
  const Axis s0_128 = log_axis(Variable::scalar_0, 128, 1.0e-5, 0.3);
  const Axis s1_128 = log_axis(Variable::scalar_1, 128, 1.0e-5, 1.0);
  const Axis s2_128 = log_axis(Variable::scalar_2, 128, 1.0e-5, 1.0);

  if (original) {
    append_product(products, "output23",
                   {r32_wide, absct16, temp128, rho128}, Weight::mass);
    append_product(products, "output24", {r128_wide, absct16, eint128}, Weight::volume);
    append_product(products, "output25", {r128_wide, absct16, vr256}, Weight::mass);
    append_product(products, "output26", {r128_wide, absct16, vp128}, Weight::mass);
    append_product(products, "output27", {r128_wide, absct16, vt128}, Weight::mass);
    append_product(products, "output28", {r128_wide, absct16}, Weight::mdot);
    append_product(products, "output29", {absct16, r128_wide}, Weight::edot);
    append_product(products, "output30",
                   {r32_core, absct16, temp128, vr128_positive}, Weight::mdot);
    append_product(products, "output31",
                   {cost16, r32_core, temp128, vr128_positive}, Weight::edot);
    append_product(products, "output32",
                   {r32_core, absct16, s0_128, vr64_positive}, Weight::mass);
    append_product(products, "output33", {r32_core, absct16, s1_128}, Weight::mass);
    append_product(products, "output34", {r32_core, absct16, s2_128}, Weight::mass);
    append_product(products, "output35", {r4_core, absct4, s1_128, s0_128}, Weight::mass);
    append_product(products, "output36", {r4_core, absct4, s2_128, s0_128}, Weight::mass);
  }

  if (science) {
    const Axis r_full = log_axis(Variable::radius, 361, 6.4e-2, 262.144);
    const Axis abs_mu = linear(Variable::abs_costheta, 64, 0.0, 1.0);
    const Axis signed_mu = linear(Variable::costheta, 128, -1.0, 1.0);
    const Axis r_phase = log_axis(Variable::radius, 64, 6.4e-2, 128.0);
    const Axis abs_mu_phase = linear(Variable::abs_costheta, 16, 0.0, 1.0);
    const Axis z_full = log_axis(Variable::absolute_z, 192, 2.0e-3, 210.0);
    const Axis rcyl_full =
        log_axis(Variable::cylindrical_radius, 192, 2.0e-3, 290.0);
    const Axis pressure_kb =
        log_axis(Variable::pressure_over_kb, 256, 1.0e-6, 1.0e10);
    const Axis nh = log_axis(Variable::hydrogen_number_density, 192, 1.0e-10, 1.0e6);
    const Axis temp_k = log_axis(Variable::temperature_kelvin, 192, 1.0, 1.0e9);
    const Axis vr = symlog_axis(Variable::velocity_r_km_s, 192, -1.0e4, 1.0e4, 5.0);
    const Axis vtangential =
        symlog_axis(Variable::velocity_theta_km_s, 192, -2.5e3, 2.5e3, 5.0);
    const Axis vphi =
        symlog_axis(Variable::velocity_phi_km_s, 192, -2.5e3, 2.5e3, 5.0);
    const Axis mach =
        symlog_axis(Variable::mach, 128, 0.0, 1.0e3, 1.0e-3);
    const Axis radial_mach =
        symlog_axis(Variable::radial_mach, 192, -1.0e3, 1.0e3, 1.0e-2);
    const Axis abs_radial_mach =
        symlog_axis(Variable::absolute_radial_mach, 128, 0.0, 1.0e3, 1.0e-2);
    const Axis abs_vr =
        symlog_axis(Variable::absolute_velocity_r_km_s, 128, 0.0, 1.0e4, 1.0);
    const Axis sound_speed =
        log_axis(Variable::sound_speed_km_s, 128, 1.0e-1, 1.0e4);
    const Axis entropy =
        log_axis(Variable::entropy_proxy, 160, 1.0e-8, 1.0e9);
    const Axis jz = symlog_axis(Variable::angular_momentum_z_kpc_km_s, 192,
                               -1.0e6, 1.0e6, 10.0);
    const Axis metallicity =
        symlog_axis(Variable::scalar_0, 96, 0.0, 0.3, 1.0e-8);
    const Axis dust_small =
        symlog_axis(Variable::scalar_1, 128, 0.0, 1.0, 1.0e-8);
    const Axis dust_large =
        symlog_axis(Variable::scalar_2, 128, 0.0, 1.0, 1.0e-8);
    const Axis nh_phase =
        log_axis(Variable::hydrogen_number_density, 96, 1.0e-10, 1.0e6);
    const Axis temp_phase =
        log_axis(Variable::temperature_kelvin, 96, 1.0, 1.0e9);
    const Axis metallicity_phase =
        symlog_axis(Variable::scalar_0, 96, 0.0, 0.3, 1.0e-8);
    const Axis dust_small_phase =
        symlog_axis(Variable::scalar_1, 96, 0.0, 1.0, 1.0e-8);
    const Axis dust_large_phase =
        symlog_axis(Variable::scalar_2, 96, 0.0, 1.0, 1.0e-8);
    const Axis dust_to_metal =
        linear(Variable::dust_to_total_metal, 96, 0.0, 1.0);
    const Axis small_fraction =
        linear(Variable::small_grain_fraction, 96, 0.0, 1.0);
    const Axis cooling_rate = symlog_axis(Variable::cooling_rate_erg_s_cm3, 192,
                                          -1.0e-14, 1.0e-14, 1.0e-30);
    const Axis cooling_time =
        symlog_axis(Variable::cooling_time_myr, 192, -1.0e8, 1.0e8, 1.0e-1);
    const Axis temp_transport =
        log_axis(Variable::temperature_kelvin, 80, 1.0, 1.0e9);
    const Axis vr_transport =
        symlog_axis(Variable::velocity_r_km_s, 80, -1.0e4, 1.0e4, 5.0);

    // Radial state conditioned on folded polar angle.
    append_product(products, "science_r_theta_pressure_volume",
                   {r_full, abs_mu, pressure_kb}, Weight::volume);
    append_product(products, "science_r_theta_density_volume",
                   {r_full, abs_mu, nh}, Weight::volume);
    append_product(products, "science_r_theta_temperature_volume",
                   {r_full, abs_mu, temp_k}, Weight::volume);
    append_product(products, "science_r_theta_temperature_edot_cool",
                   {r_full, abs_mu, temp_k}, Weight::edot_cool);
    append_product(products, "science_r_theta_vr_volume",
                   {r_full, abs_mu, vr}, Weight::volume);
    append_product(products, "science_r_theta_vtheta_mass",
                   {r_full, abs_mu, vtangential}, Weight::mass);
    append_product(products, "science_r_theta_vphi_mass",
                   {r_full, abs_mu, vphi}, Weight::mass);
    append_product(products, "science_r_theta_mach_mass",
                   {r_full, abs_mu, mach}, Weight::mass);
    append_product(products, "science_r_theta_radial_mach_mass",
                   {r_full, abs_mu, radial_mach}, Weight::mass);
    append_product(products, "science_r_theta_entropy_mass",
                   {r_full, abs_mu, entropy}, Weight::mass);
    append_product(products, "science_r_theta_jz_mass",
                   {r_full, abs_mu, jz}, Weight::mass);
    append_product(products, "science_r_theta_metallicity_mass",
                   {r_full, abs_mu, metallicity}, Weight::mass);
    append_product(products, "science_r_theta_small_dust_mass",
                   {r_full, abs_mu, dust_small}, Weight::mass);
    append_product(products, "science_r_theta_large_dust_mass",
                   {r_full, abs_mu, dust_large}, Weight::mass);
    append_product(products, "science_r_theta_edot_cool_volume",
                   {r_full, abs_mu, cooling_rate}, Weight::volume);
    append_product(products, "science_r_theta_tcool_volume",
                   {r_full, abs_mu, cooling_time}, Weight::volume);

    // Physical phase pairs retain coarse radius and folded angle.
    append_product(products, "science_phase_density_temperature_volume",
                   {r_phase, abs_mu_phase, nh_phase, temp_phase}, Weight::volume);
    append_product(products, "science_phase_density_temperature_mass",
                   {r_phase, abs_mu_phase, nh_phase, temp_phase}, Weight::mass);
    append_product(products, "science_phase_temperature_metallicity_mass",
                   {r_phase, abs_mu_phase, temp_phase, metallicity_phase},
                   Weight::mass);
    append_product(products, "science_phase_temperature_small_dust_mass",
                   {r_phase, abs_mu_phase, temp_phase, dust_small_phase},
                   Weight::mass);
    append_product(products, "science_phase_temperature_large_dust_mass",
                   {r_phase, abs_mu_phase, temp_phase, dust_large_phase},
                   Weight::mass);
    append_product(products, "science_phase_metallicity_small_dust_mass",
                   {r_phase, abs_mu_phase, metallicity_phase, dust_small_phase},
                   Weight::mass);
    append_product(products, "science_phase_metallicity_large_dust_mass",
                   {r_phase, abs_mu_phase, metallicity_phase, dust_large_phase},
                   Weight::mass);
    append_product(products, "science_phase_dust_survival",
                   {r_phase, abs_mu_phase, temp_phase, dust_to_metal},
                   Weight::total_metal_mass);
    append_product(products, "science_phase_grain_processing",
                   {r_phase, abs_mu_phase, temp_phase, small_fraction},
                   Weight::dust_mass);
    append_product(products, "science_phase_density_temperature_edot_cool",
                   {r_phase, abs_mu_phase, nh_phase, temp_phase},
                   Weight::edot_cool);

    // Signed opening-angle and thermal structure of gross outflow.
    append_product(products, "science_phase_opening_angle_mdot_out",
                   {r_full, signed_mu, temp_k}, Weight::mdot_out);
    append_product(products, "science_phase_opening_angle_edot_out",
                   {r_full, signed_mu, temp_k}, Weight::edot_out);
    append_product(products, "science_transport_geometry_mdot_in_abs",
                   {r_full, signed_mu}, Weight::mdot_in_abs);
    append_product(products, "science_transport_geometry_edot_in_abs",
                   {r_full, signed_mu}, Weight::edot_in_abs);
    append_product(products, "science_transport_geometry_edot_kin_out",
                   {r_full, signed_mu}, Weight::edot_kin_out);
    append_product(products, "science_transport_geometry_edot_kin_in_abs",
                   {r_full, signed_mu}, Weight::edot_kin_in_abs);
    append_product(products, "science_transport_geometry_edot_th_out",
                   {r_full, signed_mu}, Weight::edot_th_out);
    append_product(products, "science_transport_geometry_edot_th_in_abs",
                   {r_full, signed_mu}, Weight::edot_th_in_abs);
    append_product(products, "science_transport_geometry_ram_out",
                   {r_full, signed_mu}, Weight::ram_out);

    // Four-dimensional thermokinematic transport.
    const auto append_thermokinematic =
        [&](const std::string &id, Weight weight) {
          append_product(products, id,
                         {r_phase, abs_mu_phase, temp_transport, vr_transport},
                         weight);
        };
    append_thermokinematic("science_transport_thermokinematic_mdot_net",
                           Weight::mdot);
    append_thermokinematic("science_transport_thermokinematic_mdot_out",
                           Weight::mdot_out);
    append_thermokinematic("science_transport_thermokinematic_mdot_in_abs",
                           Weight::mdot_in_abs);
    append_thermokinematic("science_transport_thermokinematic_edot_net",
                           Weight::edot);
    append_thermokinematic("science_transport_thermokinematic_edot_out",
                           Weight::edot_out);
    append_thermokinematic("science_transport_thermokinematic_edot_in_abs",
                           Weight::edot_in_abs);
    append_thermokinematic("science_transport_thermokinematic_edot_kin_net",
                           Weight::edot_kin);
    append_thermokinematic("science_transport_thermokinematic_edot_kin_out",
                           Weight::edot_kin_out);
    append_thermokinematic("science_transport_thermokinematic_edot_kin_in_abs",
                           Weight::edot_kin_in_abs);
    append_thermokinematic("science_transport_thermokinematic_edot_th_net",
                           Weight::edot_th);
    append_thermokinematic("science_transport_thermokinematic_edot_th_out",
                           Weight::edot_th_out);
    append_thermokinematic("science_transport_thermokinematic_edot_th_in_abs",
                           Weight::edot_th_in_abs);
    append_thermokinematic("science_transport_thermokinematic_edot_cool",
                           Weight::edot_cool);

    // Composition transport by thermal phase.
    append_product(products, "science_gas_metal_mdot_out",
                   {r_full, temp_k}, Weight::metal_mdot_out);
    append_product(products, "science_gas_metal_mdot_in_abs",
                   {r_full, temp_k}, Weight::metal_mdot_in_abs);
    append_product(products, "science_total_metal_mdot_out",
                   {r_full, temp_k}, Weight::total_metal_mdot_out);
    append_product(products, "science_total_metal_mdot_in_abs",
                   {r_full, temp_k}, Weight::total_metal_mdot_in_abs);
    append_product(products, "science_dust_mdot_out",
                   {r_full, temp_k}, Weight::dust_mdot_out);
    append_product(products, "science_dust_mdot_in_abs",
                   {r_full, temp_k}, Weight::dust_mdot_in_abs);

    // Radial dynamics.
    append_product(products, "science_radial_mach_mass",
                   {r_full, abs_radial_mach}, Weight::mass);
    append_product(products, "science_radial_mach_mdot_out",
                   {r_full, abs_radial_mach}, Weight::mdot_out);
    append_product(products, "science_radial_mach_mdot_in_abs",
                   {r_full, abs_radial_mach}, Weight::mdot_in_abs);
    append_product(products, "science_vr_sound_speed_mdot_out",
                   {abs_vr, sound_speed}, Weight::mdot_out);
    append_product(products, "science_vr_sound_speed_mdot_in_abs",
                   {abs_vr, sound_speed}, Weight::mdot_in_abs);
    append_product(products, "science_vr_sound_speed_edot_out",
                   {abs_vr, sound_speed}, Weight::edot_out);
    append_product(products, "science_vr_sound_speed_ram_out",
                   {abs_vr, sound_speed}, Weight::ram_out);

    // Vertical launch and breakout.
    append_product(products, "science_vertical_geometry_mdot_out",
                   {rcyl_full, z_full}, Weight::vertical_mdot_out);
    append_product(products, "science_vertical_geometry_mdot_in_abs",
                   {rcyl_full, z_full}, Weight::vertical_mdot_in_abs);
    append_product(products, "science_vertical_geometry_edot_out",
                   {rcyl_full, z_full}, Weight::vertical_edot_out);
    append_product(products, "science_vertical_geometry_ram_out",
                   {rcyl_full, z_full}, Weight::vertical_ram_out);
    append_product(products, "science_vertical_phase_mdot_out",
                   {z_full, temp_k}, Weight::vertical_mdot_out);
    append_product(products, "science_vertical_phase_edot_out",
                   {z_full, temp_k}, Weight::vertical_edot_out);

    // Angular-momentum transport.
    append_product(products, "science_angular_momentum_mdot_out",
                   {r_full, jz}, Weight::mdot_out);
    append_product(products, "science_angular_momentum_mdot_in_abs",
                   {r_full, jz}, Weight::mdot_in_abs);
  }

  if (!only_products.empty()) {
    std::vector<Product> filtered;
    filtered.reserve(only_products.size());
    for (const std::string &requested : only_products) {
      const auto match = std::find_if(
          products.begin(), products.end(),
          [&requested](const Product &product) { return product.id == requested; });
      if (match == products.end()) {
        throw std::runtime_error(
            "--only-product is not in the selected product set: " + requested);
      }
      if (std::any_of(
              filtered.begin(), filtered.end(),
              [&requested](const Product &product) {
                return product.id == requested;
              })) {
        throw std::runtime_error("Duplicate --only-product: " + requested);
      }
      filtered.push_back(*match);
    }
    products = std::move(filtered);
  }

  std::uint64_t offset = 0;
  for (auto &product : products) {
    if (product.axes.empty() || product.axes.size() > kMaxDims) {
      throw std::runtime_error("Invalid dimensionality for " + product.id);
    }
    product.device.offset = offset;
    product.device.ndim = static_cast<int>(product.axes.size());
    product.device.weight = static_cast<int>(product.weight);
    std::uint64_t stride = 1;
    for (int d = product.device.ndim - 1; d >= 0; --d) {
      const Axis &axis = product.axes[d];
      if (axis.nbin <= 0 || !(axis.maximum > axis.minimum)) {
        throw std::runtime_error("Invalid axis for " + product.id);
      }
      if (axis.scale == Scale::log && !(axis.minimum > 0.0)) {
        throw std::runtime_error("Log axis has a non-positive minimum for " + product.id);
      }
      DeviceAxis device_axis;
      device_axis.variable = static_cast<int>(axis.variable);
      device_axis.scale = static_cast<int>(axis.scale);
      device_axis.nbin = axis.nbin;
      device_axis.stride = static_cast<int>(stride);
      device_axis.minimum = axis.minimum;
      device_axis.maximum = axis.maximum;
      device_axis.linthresh = axis.linthresh;
      device_axis.transformed_minimum =
          transform_host(axis.minimum, axis.scale, axis.linthresh);
      const double transformed_max =
          transform_host(axis.maximum, axis.scale, axis.linthresh);
      device_axis.inverse_step =
          static_cast<double>(axis.nbin) /
          (transformed_max - device_axis.transformed_minimum);
      product.device.axes[d] = device_axis;
      stride *= static_cast<std::uint64_t>(axis.nbin + 2);
    }
    if (stride > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Product is too large for MPI reduction: " + product.id);
    }
    product.device.total_bins = static_cast<int>(stride);
    offset += stride;
  }
  return products;
}

bool needs_cooling_values(const std::vector<Product> &products) {
  for (const Product &product : products) {
    if (product.weight == Weight::edot_cool) return true;
    for (const Axis &axis : product.axes) {
      if (axis.variable == Variable::cooling_rate_erg_s_cm3 ||
          axis.variable == Variable::cooling_time_myr) {
        return true;
      }
    }
  }
  return false;
}

std::uint64_t total_bins(const std::vector<Product> &products) {
  if (products.empty()) return 0;
  const auto &last = products.back().device;
  return last.offset + static_cast<std::uint64_t>(last.total_bins);
}

KOKKOS_INLINE_FUNCTION
double symlog_forward_device(double value, double linthresh) {
  const double absolute = fabs(value);
  const double transformed =
      absolute <= linthresh ? absolute / linthresh
                            : 1.0 + log10(absolute / linthresh);
  return copysign(transformed, value);
}

KOKKOS_INLINE_FUNCTION
double transform_device(double value, int scale, double linthresh) {
  if (scale == static_cast<int>(Scale::log)) return log10(value);
  if (scale == static_cast<int>(Scale::symlog)) {
    return symlog_forward_device(value, linthresh);
  }
  return value;
}

KOKKOS_INLINE_FUNCTION
int bin_index(double value, const DeviceAxis &axis) {
  if (value < axis.minimum) return 0;
  if (!(value < axis.maximum)) return axis.nbin + 1;
  const double transformed = transform_device(value, axis.scale, axis.linthresh);
  const double position =
      (transformed - axis.transformed_minimum) * axis.inverse_step;
  if (position < 0.0) return 0;
  if (position >= static_cast<double>(axis.nbin)) return axis.nbin;
  return static_cast<int>(position) + 1;
}

struct CellValues {
  double radius;
  double costheta;
  double abs_costheta;
  double density;
  double temperature;
  double internal_energy;
  double velocity_r;
  double velocity_theta;
  double velocity_phi;
  double scalar_0;
  double scalar_1;
  double scalar_2;
  double speed;
  double mach;
  double radial_mach;
  double entropy_proxy;
  double angular_momentum_z;
  double angular_momentum_z_kpc_km_s;
  double temperature_kelvin;
  double hydrogen_number_density;
  double pressure_over_kb;
  double velocity_r_km_s;
  double velocity_theta_km_s;
  double velocity_phi_km_s;
  double absolute_velocity_r_km_s;
  double sound_speed_km_s;
  double absolute_radial_mach;
  double absolute_z;
  double cylindrical_radius;
  double dust_to_total_metal;
  double small_grain_fraction;
  double cooling_rate_erg_s_cm3;
  double cooling_time_myr;
  double cooling_luminosity_erg_s;
  double vertical_velocity_outward;
  double total_metal_fraction;
  double dust_fraction;
  double enthalpy_plus_ke;
  double edot;
  double edot_kin;
  double edot_th;

  KOKKOS_INLINE_FUNCTION
  double value(int variable) const {
    switch (static_cast<Variable>(variable)) {
      case Variable::radius: return radius;
      case Variable::costheta: return costheta;
      case Variable::abs_costheta: return abs_costheta;
      case Variable::density: return density;
      case Variable::temperature: return temperature;
      case Variable::internal_energy: return internal_energy;
      case Variable::velocity_r: return velocity_r;
      case Variable::velocity_theta: return velocity_theta;
      case Variable::velocity_phi: return velocity_phi;
      case Variable::scalar_0: return scalar_0;
      case Variable::scalar_1: return scalar_1;
      case Variable::scalar_2: return scalar_2;
      case Variable::speed: return speed;
      case Variable::mach: return mach;
      case Variable::radial_mach: return radial_mach;
      case Variable::entropy_proxy: return entropy_proxy;
      case Variable::angular_momentum_z: return angular_momentum_z;
      case Variable::angular_momentum_z_kpc_km_s:
        return angular_momentum_z_kpc_km_s;
      case Variable::temperature_kelvin: return temperature_kelvin;
      case Variable::hydrogen_number_density: return hydrogen_number_density;
      case Variable::pressure_over_kb: return pressure_over_kb;
      case Variable::velocity_r_km_s: return velocity_r_km_s;
      case Variable::velocity_theta_km_s: return velocity_theta_km_s;
      case Variable::velocity_phi_km_s: return velocity_phi_km_s;
      case Variable::absolute_velocity_r_km_s: return absolute_velocity_r_km_s;
      case Variable::sound_speed_km_s: return sound_speed_km_s;
      case Variable::absolute_radial_mach: return absolute_radial_mach;
      case Variable::absolute_z: return absolute_z;
      case Variable::cylindrical_radius: return cylindrical_radius;
      case Variable::dust_to_total_metal: return dust_to_total_metal;
      case Variable::small_grain_fraction: return small_grain_fraction;
      case Variable::cooling_rate_erg_s_cm3: return cooling_rate_erg_s_cm3;
      case Variable::cooling_time_myr: return cooling_time_myr;
    }
    return 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  double weighted_value(int weight, double volume) const {
    switch (static_cast<Weight>(weight)) {
      case Weight::volume: return volume;
      case Weight::mass: return volume * density;
      case Weight::mdot: return volume * density * velocity_r;
      case Weight::mdot_out: return volume * density * fmax(velocity_r, 0.0);
      case Weight::mdot_in: return volume * density * fmin(velocity_r, 0.0);
      case Weight::edot: return volume * edot;
      case Weight::edot_out: return volume * (velocity_r > 0.0 ? edot : 0.0);
      case Weight::edot_in: return volume * (velocity_r < 0.0 ? edot : 0.0);
      case Weight::edot_kin: return volume * edot_kin;
      case Weight::edot_th: return volume * edot_th;
      case Weight::mdot_in_abs: return volume * density * fmax(-velocity_r, 0.0);
      case Weight::edot_in_abs: return volume * (velocity_r < 0.0 ? -edot : 0.0);
      case Weight::edot_kin_out: return volume * (velocity_r > 0.0 ? edot_kin : 0.0);
      case Weight::edot_kin_in_abs:
        return volume * (velocity_r < 0.0 ? -edot_kin : 0.0);
      case Weight::edot_th_out: return volume * (velocity_r > 0.0 ? edot_th : 0.0);
      case Weight::edot_th_in_abs:
        return volume * (velocity_r < 0.0 ? -edot_th : 0.0);
      case Weight::edot_cool: return cooling_luminosity_erg_s;
      case Weight::ram_out:
        return volume * density * fmax(velocity_r, 0.0) * fmax(velocity_r, 0.0);
      case Weight::vertical_mdot_out:
        return volume * density * fmax(vertical_velocity_outward, 0.0);
      case Weight::vertical_mdot_in_abs:
        return volume * density * fmax(-vertical_velocity_outward, 0.0);
      case Weight::vertical_edot_out:
        return volume * enthalpy_plus_ke * fmax(vertical_velocity_outward, 0.0);
      case Weight::vertical_ram_out:
        return volume * density * fmax(vertical_velocity_outward, 0.0) *
               fmax(vertical_velocity_outward, 0.0);
      case Weight::metal_mdot_out:
        return volume * density * fmax(velocity_r, 0.0) * scalar_0;
      case Weight::metal_mdot_in_abs:
        return volume * density * fmax(-velocity_r, 0.0) * scalar_0;
      case Weight::total_metal_mdot_out:
        return volume * density * fmax(velocity_r, 0.0) * total_metal_fraction;
      case Weight::total_metal_mdot_in_abs:
        return volume * density * fmax(-velocity_r, 0.0) * total_metal_fraction;
      case Weight::dust_mdot_out:
        return volume * density * fmax(velocity_r, 0.0) * dust_fraction;
      case Weight::dust_mdot_in_abs:
        return volume * density * fmax(-velocity_r, 0.0) * dust_fraction;
      case Weight::total_metal_mass: return volume * density * total_metal_fraction;
      case Weight::dust_mass: return volume * density * dust_fraction;
    }
    return 0.0;
  }
};

void read_exact_at(int fd, void *buffer, std::uint64_t bytes, std::uint64_t offset,
                   const fs::path &path) {
  auto *destination = static_cast<unsigned char *>(buffer);
  std::uint64_t completed = 0;
  while (completed < bytes) {
    const ssize_t count =
        pread(fd, destination + completed, static_cast<size_t>(bytes - completed),
              static_cast<off_t>(offset + completed));
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) {
      throw std::runtime_error("Short pread from " + path.string() + ": " +
                               std::strerror(errno));
    }
    completed += static_cast<std::uint64_t>(count);
  }
}

struct CoolingTables {
  Kokkos::View<double *> temperature_bins;
  Kokkos::View<double *> hydrogen_bins;
  Kokkos::View<double *> metal_pie;
  Kokkos::View<double *> hhe_pie;
  Kokkos::View<double *> metal_cie;
  Kokkos::View<double *> hhe_cie;
};

CoolingTables make_cooling_tables() {
  CoolingTables tables;
  tables.temperature_bins = Kokkos::View<double *>("cooling_temperature_bins",
                                                   Tbins_TOTAL_SIZE);
  tables.hydrogen_bins =
      Kokkos::View<double *>("cooling_hydrogen_bins", nHbins_TOTAL_SIZE);
  tables.metal_pie = Kokkos::View<double *>("cooling_metal_pie",
                                            Metal_Cooling_TOTAL_SIZE);
  tables.hhe_pie =
      Kokkos::View<double *>("cooling_hhe_pie", H_He_Cooling_TOTAL_SIZE);
  tables.metal_cie = Kokkos::View<double *>("cooling_metal_cie",
                                            Metal_Cooling_CIE_TOTAL_SIZE);
  tables.hhe_cie = Kokkos::View<double *>("cooling_hhe_cie",
                                          H_He_Cooling_CIE_TOTAL_SIZE);

  auto temperature_bins = Kokkos::create_mirror_view(tables.temperature_bins);
  auto hydrogen_bins = Kokkos::create_mirror_view(tables.hydrogen_bins);
  auto metal_pie = Kokkos::create_mirror_view(tables.metal_pie);
  auto hhe_pie = Kokkos::create_mirror_view(tables.hhe_pie);
  auto metal_cie = Kokkos::create_mirror_view(tables.metal_cie);
  auto hhe_cie = Kokkos::create_mirror_view(tables.hhe_cie);
  for (int i = 0; i < Tbins_TOTAL_SIZE; ++i) temperature_bins(i) = Tbins_ARR[i];
  for (int i = 0; i < nHbins_TOTAL_SIZE; ++i) hydrogen_bins(i) = nHbins_ARR[i];
  for (int i = 0; i < Metal_Cooling_TOTAL_SIZE; ++i) metal_pie(i) = Metal_Cooling_ARR[i];
  for (int i = 0; i < H_He_Cooling_TOTAL_SIZE; ++i) hhe_pie(i) = H_He_Cooling_ARR[i];
  for (int i = 0; i < Metal_Cooling_CIE_TOTAL_SIZE; ++i) {
    metal_cie(i) = Metal_Cooling_CIE_ARR[i];
  }
  for (int i = 0; i < H_He_Cooling_CIE_TOTAL_SIZE; ++i) {
    hhe_cie(i) = H_He_Cooling_CIE_ARR[i];
  }
  Kokkos::deep_copy(tables.temperature_bins, temperature_bins);
  Kokkos::deep_copy(tables.hydrogen_bins, hydrogen_bins);
  Kokkos::deep_copy(tables.metal_pie, metal_pie);
  Kokkos::deep_copy(tables.hhe_pie, hhe_pie);
  Kokkos::deep_copy(tables.metal_cie, metal_cie);
  Kokkos::deep_copy(tables.hhe_cie, hhe_cie);
  return tables;
}

struct CoolingLambda {
  double cie;
  double pie;
};

KOKKOS_INLINE_FUNCTION
CoolingLambda cooling_lambdas_cgs(const CoolingTables &tables, double temperature, double nH,
                                  double metallicity) {
  const double log_temperature = log10(temperature);
  const double log_nH = log10(nH);
  const double temperature_floor = pow(10.0, tables.temperature_bins(0));
  const double temperature_ceiling =
      pow(10.0, tables.temperature_bins(Tbins_TOTAL_SIZE - 1));
  const double hydrogen_floor = pow(10.0, tables.hydrogen_bins(0));
  const double hydrogen_ceiling = pow(10.0, tables.hydrogen_bins(nHbins_TOTAL_SIZE - 1));
  const double low_temperature = temperature < temperature_floor ? 1.0 : 0.0;
  const double temperature_in_table =
      temperature >= temperature_floor && temperature <= temperature_ceiling ? 1.0 : 0.0;
  const double hydrogen_in_table =
      nH >= hydrogen_floor && nH <= hydrogen_ceiling ? 1.0 : 0.0;
  const double pie_in_table = temperature_in_table * hydrogen_in_table;

  int temperature_index = 0;
  int hydrogen_index = 0;
  while (temperature_index < Tbins_TOTAL_SIZE - 2 &&
         tables.temperature_bins(temperature_index + 1) < log_temperature) {
    ++temperature_index;
  }
  while (hydrogen_index < nHbins_TOTAL_SIZE - 2 &&
         tables.hydrogen_bins(hydrogen_index + 1) < log_nH) {
    ++hydrogen_index;
  }

  const double log_temperature0 = tables.temperature_bins(temperature_index);
  const double log_temperature1 = tables.temperature_bins(temperature_index + 1);
  const double temperature_weight =
      (log_temperature - log_temperature0) / (log_temperature1 - log_temperature0);
  const double inverse_temperature_weight = 1.0 - temperature_weight;
  const double log_hydrogen0 = tables.hydrogen_bins(hydrogen_index);
  const double log_hydrogen1 = tables.hydrogen_bins(hydrogen_index + 1);
  const double hydrogen_weight =
      (log_nH - log_hydrogen0) / (log_hydrogen1 - log_hydrogen0);
  const double inverse_hydrogen_weight = 1.0 - hydrogen_weight;
  const int pie00 = temperature_index * nHbins_TOTAL_SIZE + hydrogen_index;
  const int pie10 = (temperature_index + 1) * nHbins_TOTAL_SIZE + hydrogen_index;
  const int pie01 = temperature_index * nHbins_TOTAL_SIZE + hydrogen_index + 1;
  const int pie11 = (temperature_index + 1) * nHbins_TOTAL_SIZE + hydrogen_index + 1;

  const double primordial_pie =
      inverse_temperature_weight * inverse_hydrogen_weight * tables.hhe_pie(pie00) +
      temperature_weight * inverse_hydrogen_weight * tables.hhe_pie(pie10) +
      inverse_temperature_weight * hydrogen_weight * tables.hhe_pie(pie01) +
      temperature_weight * hydrogen_weight * tables.hhe_pie(pie11);
  const double metal_pie =
      inverse_temperature_weight * inverse_hydrogen_weight * tables.metal_pie(pie00) +
      temperature_weight * inverse_hydrogen_weight * tables.metal_pie(pie10) +
      inverse_temperature_weight * hydrogen_weight * tables.metal_pie(pie01) +
      temperature_weight * hydrogen_weight * tables.metal_pie(pie11);
  const double lambda_pie = pie_in_table * fma(metallicity, metal_pie, primordial_pie);

  const double primordial_cie = fma(
      temperature_weight, tables.hhe_cie(temperature_index + 1) -
                              tables.hhe_cie(temperature_index),
      tables.hhe_cie(temperature_index));
  const double metal_cie = fma(
      temperature_weight, tables.metal_cie(temperature_index + 1) -
                              tables.metal_cie(temperature_index),
      tables.metal_cie(temperature_index));
  const double lambda_cie_table = fma(metallicity, metal_cie, primordial_cie);
  const double low_temperature_e1 = exp(-1.184e5 / (temperature + 1.0e3));
  const double low_temperature_e2 = exp(-92.0 / temperature);
  const double lambda_low_temperature =
      metallicity *
      (2.0e-19 * low_temperature_e1 +
       2.8e-28 * sqrt(temperature) * low_temperature_e2);
  const double lambda_cie =
      temperature_in_table * lambda_cie_table +
      (1.0 - temperature_in_table) * low_temperature * lambda_low_temperature;

  return CoolingLambda{lambda_cie, lambda_pie};
}

KOKKOS_INLINE_FUNCTION
double net_cooling_rate_cgs(const CoolingTables &tables, const HeaderWire &header,
                            double x, double y, double z, double rho, double temperature,
                            double metallicity, double dx) {
  const double nH = rho * header.hydrogen_number_density_per_density;
  const CoolingLambda lambdas = cooling_lambdas_cgs(tables, temperature, nH, metallicity);

  const double cylindrical_radius_squared = fma(x, x, y * y);
  const double cylindrical_radius = sqrt(cylindrical_radius_squared);
  const double horizontal_falloff = exp(-cylindrical_radius / header.hscale_radius);
  const double vertical_scale_squared =
      header.hscale_height * header.hscale_height *
      (1.0 + cylindrical_radius_squared /
                 (header.hscale_radius * header.hscale_radius));
  const double vertical_falloff = exp(-(z * z) / vertical_scale_squared);
  double gamma_heating = header.hrate * header.hscale_norm *
                         header.hydrogen_number_density_per_density *
                         horizontal_falloff * vertical_falloff;
  const double hot = temperature > 1.0e4 ? 1.0 : 0.0;
  const double inverse_ratio = 1.0e4 / temperature;
  const double damping_factor =
      hot * inverse_ratio * inverse_ratio * inverse_ratio * inverse_ratio *
          inverse_ratio * inverse_ratio * inverse_ratio * inverse_ratio +
      (1.0 - hot);
  gamma_heating *= damping_factor;
  const double neutral_fraction =
      1.0 - 0.5 * (1.0 + tanh((temperature - 8.0e3) / 1.5e3));
  const double tau = neutral_fraction * nH * 1.0e-17 * dx * header.length_cgs;
  const double pie_fraction = exp(-tau);
  const double lambda_cooling =
      (1.0 - pie_fraction) * lambdas.cie + pie_fraction * lambdas.pie;
  return nH * (nH * lambda_cooling - gamma_heating);
}

void validate_chunk_metadata(const unsigned char *raw, int blocks,
                             const HeaderWire &header, const fs::path &path,
                             std::uint64_t first_block,
                             std::vector<BlockKey> &block_keys,
                             double &block_volume) {
  const auto left_edge = [](std::int64_t index, std::int64_t count, double minimum,
                            double maximum) {
    const double fraction = static_cast<double>(index) / static_cast<double>(count);
    return (fraction * maximum - fraction * minimum) -
           (0.5 * maximum - 0.5 * minimum) + (0.5 * minimum + 0.5 * maximum);
  };
  const auto geometry_matches = [](double actual, double expected, double minimum,
                                   double maximum) {
    const double scale =
        std::max({1.0, std::abs(minimum), std::abs(maximum), std::abs(expected)});
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    return std::abs(actual - expected) <= tolerance;
  };
  for (int block = 0; block < blocks; ++block) {
    const unsigned char *record = raw + block * header.record_size;
    std::int32_t indices[10];
    double geometry[6];
    std::memcpy(indices, record, sizeof(indices));
    std::memcpy(geometry, record + sizeof(indices), sizeof(geometry));
    const auto fail = [&path, first_block, block](const std::string &reason) {
      throw std::runtime_error("Invalid MeshBlock " +
                               std::to_string(first_block + block) + " in " +
                               path.string() + ": " + reason);
    };
    if (indices[1] - indices[0] + 1 != header.nx1 ||
        indices[3] - indices[2] + 1 != header.nx2 ||
        indices[5] - indices[4] + 1 != header.nx3) {
      fail("output index shape does not match embedded MeshBlock shape");
    }
    if (indices[6] < 0 || indices[7] < 0 || indices[8] < 0 || indices[9] < 0) {
      fail("negative logical location or refinement level");
    }
    for (const double value : geometry) {
      if (!std::isfinite(value)) fail("non-finite geometry");
    }
    if (!(geometry[1] > geometry[0]) || !(geometry[3] > geometry[2]) ||
        !(geometry[5] > geometry[4])) {
      fail("non-positive physical extent");
    }
    const int level = indices[9];
    if (level > 30) fail("refinement level is too large for logical-key validation");
    const std::int64_t level_factor = std::int64_t{1} << level;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const std::int64_t logical = indices[6 + dimension];
      const std::int64_t logical_count =
          static_cast<std::int64_t>(header.root_blocks[dimension]) * level_factor;
      if (logical >= logical_count) {
        fail("logical location lies outside the embedded root mesh");
      }
      const double expected_min =
          logical == 0
              ? header.domain_min[dimension]
              : left_edge(logical, logical_count, header.domain_min[dimension],
                          header.domain_max[dimension]);
      const double expected_max =
          logical == logical_count - 1
              ? header.domain_max[dimension]
              : left_edge(logical + 1, logical_count, header.domain_min[dimension],
                          header.domain_max[dimension]);
      if (!geometry_matches(geometry[2 * dimension], expected_min,
                            header.domain_min[dimension], header.domain_max[dimension]) ||
          !geometry_matches(geometry[2 * dimension + 1], expected_max,
                            header.domain_min[dimension], header.domain_max[dimension])) {
        fail("physical geometry does not match logical location and refinement level");
      }
    }
    block_keys.push_back(BlockKey{indices[6], indices[7], indices[8], indices[9]});
    block_volume += (geometry[1] - geometry[0]) * (geometry[3] - geometry[2]) *
                    (geometry[5] - geometry[4]);
  }
}

std::vector<fs::path> find_shards(const fs::path &input_dir,
                                  const std::string &sequence) {
  if (!fs::is_directory(input_dir)) {
    throw std::runtime_error("Input directory does not exist: " + input_dir.string());
  }
  const std::string filename = "gotham.hydro_w." + sequence + ".bin";
  std::vector<fs::path> paths;
  const fs::path direct = input_dir / filename;
  if (fs::is_regular_file(direct)) paths.push_back(direct);
  for (const auto &entry : fs::directory_iterator(input_dir)) {
    if (!entry.is_directory()) continue;
    const std::string name = entry.path().filename().string();
    if (name.rfind("node_", 0) != 0 && name.rfind("rank_", 0) != 0) continue;
    const fs::path candidate = entry.path() / filename;
    if (fs::is_regular_file(candidate)) paths.push_back(candidate);
  }
  std::sort(paths.begin(), paths.end());
  if (paths.empty()) {
    throw std::runtime_error("No shards found for " + filename + " beneath " +
                             input_dir.string());
  }
  return paths;
}

void validate_shard_manifest(const std::vector<fs::path> &paths,
                             std::uint64_t expected_shards) {
  if (expected_shards == 0) return;
  if (paths.size() != expected_shards) {
    throw std::runtime_error("Expected " + std::to_string(expected_shards) +
                             " shards, found " + std::to_string(paths.size()));
  }
  for (std::size_t index = 0; index < paths.size(); ++index) {
    std::ostringstream expected;
    expected << "node_" << std::setw(8) << std::setfill('0') << index;
    if (paths[index].parent_path().filename() != expected.str()) {
      throw std::runtime_error("Expected exact shard directory " + expected.str() +
                               ", found " +
                               paths[index].parent_path().filename().string());
    }
  }
}

std::vector<fs::path> broadcast_paths(std::vector<fs::path> paths, int rank) {
  std::string packed;
  if (rank == 0) {
    for (const auto &path : paths) {
      packed += path.string();
      packed.push_back('\n');
    }
  }
  std::uint64_t size = packed.size();
  MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  if (rank != 0) packed.resize(size);
  if (size > 0) {
    if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Packed shard manifest is too large");
    }
    MPI_Bcast(packed.data(), static_cast<int>(size), MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  if (rank != 0) {
    paths.clear();
    std::istringstream stream(packed);
    std::string line;
    while (std::getline(stream, line)) {
      if (!line.empty()) paths.emplace_back(line);
    }
  }
  return paths;
}

void process_chunk(const Kokkos::View<unsigned char *> &raw, int blocks,
                   const HeaderWire &header,
                   const Kokkos::View<DeviceProduct *> &products,
                   const Kokkos::View<double *> &histogram,
                   const Kokkos::View<std::uint64_t *> &invalid_cells,
                   const CoolingTables &cooling_tables, bool calculate_cooling) {
  const int product_count = static_cast<int>(products.extent(0));
  const int nx1 = header.nx1;
  const int nx2 = header.nx2;
  const int nx3 = header.nx3;
  const int cells = nx1 * nx2 * nx3;
  const std::uint64_t record_size = header.record_size;
  const FieldMap fields = header.fields;
  const double gamma = header.gamma;
  const double velocity_km_s = header.velocity_km_s;
  const double temperature_kelvin_per_ratio = header.temperature_kelvin_per_ratio;
  const double hydrogen_number_density_per_density =
      header.hydrogen_number_density_per_density;
  const double pressure_over_kb_per_internal_energy =
      header.pressure_over_kb_per_internal_energy;
  const double time_cgs = header.time_cgs;
  const double pressure_cgs = header.pressure_cgs;
  const int total_cells = blocks * cells;

  Kokkos::parallel_for(
      "gotham_pdf_fused_chunk", Kokkos::RangePolicy<>(0, total_cells),
      KOKKOS_LAMBDA(const int linear) {
        const int block = linear / cells;
        const int cell = linear - block * cells;
        const int i = cell % nx1;
        const int j = (cell / nx1) % nx2;
        const int k = cell / (nx1 * nx2);

        const unsigned char *record = raw.data() + block * record_size;
        const double *geometry =
            reinterpret_cast<const double *>(record + 10 * sizeof(std::int32_t));
        const float *data =
            reinterpret_cast<const float *>(record + kRecordMetadataBytes);
        const double x =
            geometry[0] + (static_cast<double>(i) + 0.5) *
                              (geometry[1] - geometry[0]) / nx1;
        const double y =
            geometry[2] + (static_cast<double>(j) + 0.5) *
                              (geometry[3] - geometry[2]) / nx2;
        const double z =
            geometry[4] + (static_cast<double>(k) + 0.5) *
                              (geometry[5] - geometry[4]) / nx3;
        const double volume = ((geometry[1] - geometry[0]) / nx1) *
                              ((geometry[3] - geometry[2]) / nx2) *
                              ((geometry[5] - geometry[4]) / nx3);
        const double rho = data[fields.density * cells + cell];
        const double vx = data[fields.velocity_x * cells + cell];
        const double vy = data[fields.velocity_y * cells + cell];
        const double vz = data[fields.velocity_z * cells + cell];
        const double eint = data[fields.internal_energy * cells + cell];
        const double scalar0 = data[fields.scalar_0 * cells + cell];
        const double scalar1 = data[fields.scalar_1 * cells + cell];
        const double scalar2 = data[fields.scalar_2 * cells + cell];
        if (!Kokkos::isfinite(rho) || !Kokkos::isfinite(vx) ||
            !Kokkos::isfinite(vy) || !Kokkos::isfinite(vz) ||
            !Kokkos::isfinite(eint) || !Kokkos::isfinite(scalar0) ||
            !Kokkos::isfinite(scalar1) || !Kokkos::isfinite(scalar2) ||
            !(rho > 0.0) || !(eint > 0.0)) {
          Kokkos::atomic_inc(&invalid_cells(0));
          return;
        }

        const double cylindrical_radius = sqrt(x * x + y * y);
        const double radius = sqrt(cylindrical_radius * cylindrical_radius + z * z);
        const double costheta = radius > 0.0 ? z / radius : 1.0;
        const double radial_velocity =
            radius > 0.0 ? (vx * x + vy * y + vz * z) / radius : 0.0;
        const double theta_velocity =
            radius > 0.0 && cylindrical_radius > 0.0
                ? z * (vx * x + vy * y) / (radius * cylindrical_radius) -
                      vz * cylindrical_radius / radius
                : 0.0;
        const double phi_velocity =
            cylindrical_radius > 0.0 ? (-vx * y + vy * x) / cylindrical_radius : 0.0;
        const double velocity_squared = vx * vx + vy * vy + vz * vz;
        const double speed = sqrt(velocity_squared);
        const double temperature = eint / rho;
        const double sound_speed = sqrt(gamma * (gamma - 1.0) * eint / rho);
        const double edot_kin = 0.5 * rho * velocity_squared * radial_velocity;
        const double edot_th = gamma * eint * radial_velocity;
        const double enthalpy_plus_ke = 0.5 * rho * velocity_squared + gamma * eint;
        const double dust_fraction = scalar1 + scalar2;
        const double total_metal_fraction = scalar0 + dust_fraction;
        const double vertical_velocity_outward = z >= 0.0 ? vz : -vz;
        double cooling_rate_erg_s_cm3 = 0.0;
        double cooling_time_myr = 0.0;
        double cooling_luminosity_erg_s = 0.0;
        if (calculate_cooling) {
          // Match AthenaK's derived cooling_time output, which reports the
          // physical net cooling/heating denominator before the separate
          // high-temperature ceiling safeguard is applied.
          cooling_rate_erg_s_cm3 =
              net_cooling_rate_cgs(cooling_tables, header, x, y, z, rho,
                                   temperature * temperature_kelvin_per_ratio,
                                   scalar0 / kCoolingSolarMetallicity,
                                   (geometry[1] - geometry[0]) / nx1);
          // AthenaK adds FLT_MIN to the code-unit source rate. In physical
          // units that term is FLT_MIN * pressure_cgs / time_cgs.
          cooling_time_myr =
              eint * pressure_cgs /
              (cooling_rate_erg_s_cm3 + FLT_MIN * pressure_cgs / time_cgs) /
              kMyrCgs;
          cooling_luminosity_erg_s =
              cooling_rate_erg_s_cm3 * volume * header.length_cgs * header.length_cgs *
              header.length_cgs;
        }

        CellValues values;
        values.radius = radius;
        values.costheta = costheta;
        values.abs_costheta = fabs(costheta);
        values.density = rho;
        values.temperature = temperature;
        values.internal_energy = eint;
        values.velocity_r = radial_velocity;
        values.velocity_theta = theta_velocity;
        values.velocity_phi = phi_velocity;
        values.scalar_0 = scalar0;
        values.scalar_1 = scalar1;
        values.scalar_2 = scalar2;
        values.speed = speed;
        values.mach = speed / sound_speed;
        values.radial_mach = radial_velocity / sound_speed;
        values.entropy_proxy = (gamma - 1.0) * eint / pow(rho, gamma);
        values.angular_momentum_z = x * vy - y * vx;
        values.angular_momentum_z_kpc_km_s =
            values.angular_momentum_z * velocity_km_s;
        values.temperature_kelvin = temperature * temperature_kelvin_per_ratio;
        values.hydrogen_number_density =
            rho * hydrogen_number_density_per_density;
        values.pressure_over_kb =
            eint * pressure_over_kb_per_internal_energy;
        values.velocity_r_km_s = radial_velocity * velocity_km_s;
        values.velocity_theta_km_s = theta_velocity * velocity_km_s;
        values.velocity_phi_km_s = phi_velocity * velocity_km_s;
        values.absolute_velocity_r_km_s = fabs(radial_velocity) * velocity_km_s;
        values.sound_speed_km_s = sound_speed * velocity_km_s;
        values.absolute_radial_mach = fabs(radial_velocity / sound_speed);
        values.absolute_z = fabs(z);
        values.cylindrical_radius = cylindrical_radius;
        values.dust_to_total_metal =
            total_metal_fraction > 0.0 ? dust_fraction / total_metal_fraction : 0.0;
        values.small_grain_fraction =
            dust_fraction > 0.0 ? scalar1 / dust_fraction : 0.0;
        values.cooling_rate_erg_s_cm3 = cooling_rate_erg_s_cm3;
        values.cooling_time_myr = cooling_time_myr;
        values.cooling_luminosity_erg_s = cooling_luminosity_erg_s;
        values.vertical_velocity_outward = vertical_velocity_outward;
        values.total_metal_fraction = total_metal_fraction;
        values.dust_fraction = dust_fraction;
        values.enthalpy_plus_ke = enthalpy_plus_ke;
        values.edot_kin = edot_kin;
        values.edot_th = edot_th;
        values.edot = edot_kin + edot_th;

        for (int product_index = 0; product_index < product_count; ++product_index) {
          const DeviceProduct product = products(product_index);
          int flat = 0;
          for (int dimension = 0; dimension < product.ndim; ++dimension) {
            const DeviceAxis axis = product.axes[dimension];
            flat += bin_index(values.value(axis.variable), axis) * axis.stride;
          }
          const double weight = values.weighted_value(product.weight, volume);
          Kokkos::atomic_add(&histogram(product.offset + flat), weight);
        }
      });
}

struct LocalStats {
  std::uint64_t shards = 0;
  std::uint64_t blocks = 0;
  std::uint64_t cells = 0;
  std::uint64_t bytes = 0;
  double block_volume = 0.0;
};

LocalStats process_shards(const std::vector<fs::path> &paths, const Options &options,
                          int rank, int nranks, const HeaderWire &reference,
                          const Kokkos::View<DeviceProduct *> &device_products,
                          const Kokkos::View<double *> &histogram,
                          std::vector<BlockKey> &block_keys,
                          const CoolingTables &cooling_tables, bool calculate_cooling) {
  LocalStats stats;
  const std::uint64_t capacity_bytes =
      static_cast<std::uint64_t>(options.chunk_blocks) * reference.record_size;
  Kokkos::View<unsigned char *, Kokkos::HostSpace> host_raw("host_raw", capacity_bytes);
  Kokkos::View<unsigned char *> device_raw("device_raw", capacity_bytes);
  Kokkos::View<std::uint64_t *> invalid_cells("invalid_cells", 1);
  Kokkos::deep_copy(invalid_cells, static_cast<std::uint64_t>(0));

  for (std::size_t shard_index = rank; shard_index < paths.size();
       shard_index += static_cast<std::size_t>(nranks)) {
    const fs::path &path = paths[shard_index];
    const HeaderInfo info = read_header(path);
    validate_header(info, reference, path);
    const std::uint64_t blocks =
        options.max_blocks_per_shard == 0
            ? info.blocks
            : std::min(info.blocks, options.max_blocks_per_shard);
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
      throw std::runtime_error("Unable to open shard with POSIX I/O: " + path.string());
    }
    try {
      for (std::uint64_t first = 0; first < blocks;
           first += static_cast<std::uint64_t>(options.chunk_blocks)) {
        const int chunk_blocks = static_cast<int>(
            std::min<std::uint64_t>(options.chunk_blocks, blocks - first));
        const std::uint64_t chunk_bytes =
            static_cast<std::uint64_t>(chunk_blocks) * info.record_size;
        read_exact_at(fd, host_raw.data(), chunk_bytes,
                      info.data_offset + first * info.record_size, path);
        validate_chunk_metadata(host_raw.data(), chunk_blocks, reference, path, first,
                                block_keys, stats.block_volume);
        auto host_chunk = Kokkos::subview(
            host_raw, std::make_pair(static_cast<std::uint64_t>(0), chunk_bytes));
        auto device_chunk = Kokkos::subview(
            device_raw, std::make_pair(static_cast<std::uint64_t>(0), chunk_bytes));
        Kokkos::deep_copy(device_chunk, host_chunk);
        process_chunk(device_raw, chunk_blocks, reference, device_products, histogram,
                      invalid_cells, cooling_tables, calculate_cooling);
        Kokkos::fence();
        stats.blocks += static_cast<std::uint64_t>(chunk_blocks);
        stats.cells += static_cast<std::uint64_t>(chunk_blocks) * reference.nx1 *
                       reference.nx2 * reference.nx3;
        stats.bytes += chunk_bytes;
      }
    } catch (...) {
      close(fd);
      throw;
    }
    close(fd);
    ++stats.shards;
  }
  auto host_invalid = Kokkos::create_mirror_view(invalid_cells);
  Kokkos::deep_copy(host_invalid, invalid_cells);
  if (host_invalid(0) != 0) {
    throw std::runtime_error("Encountered " + std::to_string(host_invalid(0)) +
                             " cells with non-finite or non-positive hydro state");
  }
  return stats;
}

void validate_global_mesh(const std::vector<BlockKey> &local_keys,
                          const LocalStats &local_stats,
                          const HeaderWire &header, int rank, int nranks) {
  if (local_keys.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() /
                                                   sizeof(BlockKey))) {
    throw std::runtime_error("Local AMR key list is too large for MPI_Gatherv");
  }
  const int local_bytes = static_cast<int>(local_keys.size() * sizeof(BlockKey));
  std::vector<int> byte_counts(rank == 0 ? nranks : 0);
  MPI_Gather(&local_bytes, 1, MPI_INT, rank == 0 ? byte_counts.data() : nullptr, 1,
             MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displacements;
  std::vector<BlockKey> global_keys;
  if (rank == 0) {
    displacements.resize(nranks);
    int total_bytes = 0;
    for (int r = 0; r < nranks; ++r) {
      displacements[r] = total_bytes;
      if (byte_counts[r] > std::numeric_limits<int>::max() - total_bytes) {
        throw std::runtime_error("Global AMR key list is too large for MPI_Gatherv");
      }
      total_bytes += byte_counts[r];
    }
    global_keys.resize(static_cast<std::size_t>(total_bytes) / sizeof(BlockKey));
  }
  MPI_Gatherv(local_keys.data(), local_bytes, MPI_BYTE,
              rank == 0 ? global_keys.data() : nullptr,
              rank == 0 ? byte_counts.data() : nullptr,
              rank == 0 ? displacements.data() : nullptr, MPI_BYTE, 0, MPI_COMM_WORLD);

  double global_volume = 0.0;
  MPI_Reduce(&local_stats.block_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  int valid = 1;
  if (rank == 0) {
    try {
      std::unordered_set<BlockKey, BlockKeyHash> leaves;
      leaves.reserve(global_keys.size() * 2);
      for (const BlockKey &key : global_keys) {
        if (!leaves.insert(key).second) {
          throw std::runtime_error("Duplicate AMR leaf key detected");
        }
      }
      for (const BlockKey &key : global_keys) {
        for (int ancestor_level = key.level - 1; ancestor_level >= 0;
             --ancestor_level) {
          const int shift = key.level - ancestor_level;
          const BlockKey ancestor{key.lx1 >> shift, key.lx2 >> shift,
                                  key.lx3 >> shift, ancestor_level};
          if (leaves.find(ancestor) != leaves.end()) {
            throw std::runtime_error("Ancestor/descendant AMR leaf overlap detected");
          }
        }
      }
      const double relative_volume_error =
          std::abs(global_volume - header.domain_volume) / header.domain_volume;
      if (relative_volume_error > 2.0e-9) {
        std::ostringstream message;
        message << "AMR leaves do not close the domain volume: got "
                << std::setprecision(17) << global_volume << ", expected "
                << header.domain_volume << ", relative error " << relative_volume_error;
        throw std::runtime_error(message.str());
      }
      std::cout << "validated AMR leaves=" << global_keys.size()
                << " domain_volume=" << std::setprecision(17) << global_volume << "\n";
    } catch (const std::exception &error) {
      valid = 0;
      std::cerr << "global AMR validation failed: " << error.what() << "\n";
    }
  }
  MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!valid) throw std::runtime_error("Global AMR validation failed");
}

std::vector<double> axis_edges(const Axis &axis) {
  std::vector<double> edges(axis.nbin + 1);
  const double transformed_min = transform_host(axis.minimum, axis.scale, axis.linthresh);
  const double transformed_max = transform_host(axis.maximum, axis.scale, axis.linthresh);
  for (int i = 0; i <= axis.nbin; ++i) {
    const double transformed =
        transformed_min + (transformed_max - transformed_min) *
                              static_cast<double>(i) / axis.nbin;
    if (axis.scale == Scale::log) {
      edges[i] = std::pow(10.0, transformed);
    } else if (axis.scale == Scale::symlog) {
      edges[i] = symlog_inverse_host(transformed, axis.linthresh);
    } else {
      edges[i] = transformed;
    }
  }
  return edges;
}

void write_product(const fs::path &output_dir, const std::string &basename,
                   const std::string &number, double time, const Product &product,
                   const std::vector<double> &values) {
  const fs::path directory = output_dir / product.id;
  fs::create_directories(directory);
  const fs::path header_path = directory / (basename + ".header.pdf");
  const fs::path data_path = directory / (basename + "." + number + ".pdf");
  const fs::path partial_header = header_path.string() + ".partial";
  const fs::path partial_data = data_path.string() + ".partial";
  fs::remove(partial_header);
  fs::remove(partial_data);

  std::ofstream header(partial_header);
  if (!header) throw std::runtime_error("Unable to create " + partial_header.string());
  header << "# AthenaK N-D PDF Output reconstructed from GOTHAM full-volume shards\n";
  header << "format = dense\n";
  header << "distribution = global_rebuild\n";
  header << "ndim = " << product.axes.size() << "\n";
  if (product.weight == Weight::volume || product.weight == Weight::mass) {
    header << "weight = " << weight_name(product.weight) << "\n";
  } else {
    header << "weight = variable\n";
    header << "weight_variable = " << weight_name(product.weight) << "\n";
  }
  header << "total_bins = " << product.device.total_bins << "\n\n";
  header << std::scientific << std::setprecision(15);
  for (int d = 0; d < static_cast<int>(product.axes.size()); ++d) {
    const Axis &axis = product.axes[d];
    const DeviceAxis &device_axis = product.device.axes[d];
    header << "# Dimension " << d + 1 << "\n";
    header << "variable_" << d + 1 << " = " << variable_name(axis.variable) << "\n";
    header << "nbin" << d + 1 << " = " << axis.nbin << "\n";
    header << "bin" << d + 1 << "_min = " << axis.minimum << "\n";
    header << "bin" << d + 1 << "_max = " << axis.maximum << "\n";
    header << "scale" << d + 1 << " = " << scale_name(axis.scale) << "\n";
    if (axis.scale == Scale::symlog) {
      header << "linthresh" << d + 1 << " = " << axis.linthresh << "\n";
    }
    header << "logscale" << d + 1 << " = "
           << (axis.scale == Scale::log ? "true" : "false") << "\n";
    header << "stride" << d + 1 << " = " << device_axis.stride << "\n\n";
  }
  header << "# Bin edges (nbin+1 values per dimension)\n";
  for (int d = 0; d < static_cast<int>(product.axes.size()); ++d) {
    header << "bin_edges_" << d + 1 << " =";
    for (const double edge : axis_edges(product.axes[d])) header << " " << edge;
    header << "\n";
  }
  header.close();
  if (!header) throw std::runtime_error("Failed while writing " + partial_header.string());

  std::ofstream data(partial_data, std::ios::binary);
  if (!data) throw std::runtime_error("Unable to create " + partial_data.string());
  data.write(reinterpret_cast<const char *>(&time), sizeof(time));
  data.write(reinterpret_cast<const char *>(values.data()),
             static_cast<std::streamsize>(values.size() * sizeof(double)));
  data.close();
  if (!data) throw std::runtime_error("Failed while writing " + partial_data.string());
  fs::rename(partial_header, header_path);
  fs::rename(partial_data, data_path);
}

void write_manifest(const fs::path &output_dir, const Options &options,
                    const HeaderWire &header, const std::vector<fs::path> &paths,
                    const std::vector<Product> &products, const LocalStats &global_stats,
                    double elapsed, const std::vector<double> &product_sums) {
  fs::create_directories(output_dir);
  const fs::path path = output_dir / "rebuild_manifest.json";
  const fs::path partial = path.string() + ".partial";
  fs::remove(partial);
  std::ofstream output(partial);
  if (!output) throw std::runtime_error("Unable to create " + partial.string());
  output << std::setprecision(17);
  output << "{\n";
  output << "  \"source_input_dir\": \"" << options.input_dir.string() << "\",\n";
  output << "  \"source_sequence\": \"" << options.sequence << "\",\n";
  output << "  \"output_number\": \"" << options.output_number << "\",\n";
  output << "  \"source_time\": " << header.time << ",\n";
  output << "  \"output_time\": "
         << (std::isfinite(options.output_time) ? options.output_time : header.time)
         << ",\n";
  output << "  \"source_cycle\": " << header.cycle << ",\n";
  output << "  \"source_header_fnv1a64\": " << header.header_digest << ",\n";
  output << "  \"gamma\": " << header.gamma << ",\n";
  output << "  \"product_set\": \"" << options.products << "\",\n";
  output << "  \"product_filter\": [";
  for (std::size_t i = 0; i < options.only_products.size(); ++i) {
    if (i > 0) output << ", ";
    output << "\"" << options.only_products[i] << "\"";
  }
  output << "],\n";
  output << "  \"geometry_to_logical_key_validated\": "
         << ((options.max_shards == 0 && options.max_blocks_per_shard == 0) ? "true"
                                                                            : "false")
         << ",\n";
  output << "  \"domain_min\": [" << header.domain_min[0] << ", "
         << header.domain_min[1] << ", " << header.domain_min[2] << "],\n";
  output << "  \"domain_max\": [" << header.domain_max[0] << ", "
         << header.domain_max[1] << ", " << header.domain_max[2] << "],\n";
  output << "  \"root_meshblocks\": [" << header.root_blocks[0] << ", "
         << header.root_blocks[1] << ", " << header.root_blocks[2] << "],\n";
  output << "  \"expected_domain_volume\": " << header.domain_volume << ",\n";
  output << "  \"summed_leaf_volume\": " << global_stats.block_volume << ",\n";
  output << "  \"shards_available\": " << paths.size() << ",\n";
  output << "  \"shards_processed\": " << global_stats.shards << ",\n";
  output << "  \"meshblocks_processed\": " << global_stats.blocks << ",\n";
  output << "  \"cells_processed\": " << global_stats.cells << ",\n";
  output << "  \"payload_bytes_read\": " << global_stats.bytes << ",\n";
  const std::uint64_t histogram_bytes = total_bins(products) * sizeof(double);
  const std::uint64_t raw_buffer_bytes =
      static_cast<std::uint64_t>(options.chunk_blocks) * header.record_size;
  output << "  \"histogram_bytes_per_rank\": " << histogram_bytes << ",\n";
  output << "  \"raw_buffer_bytes_per_rank\": " << raw_buffer_bytes << ",\n";
  output << "  \"planned_peak_buffer_bytes_per_rank\": "
         << 2 * (histogram_bytes + raw_buffer_bytes) << ",\n";
  output << "  \"elapsed_seconds\": " << elapsed << ",\n";
  output << "  \"products\": [\n";
  for (std::size_t i = 0; i < products.size(); ++i) {
    output << "    {\"id\": \"" << products[i].id << "\", \"weight\": \""
           << weight_name(products[i].weight) << "\", \"total_bins\": "
           << products[i].device.total_bins << ", \"sum\": " << product_sums[i] << "}";
    output << (i + 1 == products.size() ? "\n" : ",\n");
  }
  output << "  ]\n";
  output << "}\n";
  output.close();
  if (!output) throw std::runtime_error("Failed while writing " + partial.string());
  fs::rename(partial, path);
}

void reduce_and_write(const Options &options, const HeaderWire &header,
                      const std::vector<fs::path> &paths,
                      const std::vector<Product> &products,
                      const Kokkos::View<double *> &histogram, const LocalStats &local_stats,
                      double start_time, int rank) {
  auto host_histogram = Kokkos::create_mirror_view(histogram);
  Kokkos::deep_copy(host_histogram, histogram);
  Kokkos::fence();

  LocalStats global_stats;
  MPI_Reduce(&local_stats.shards, &global_stats.shards, 4, MPI_UNSIGNED_LONG_LONG,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_stats.block_volume, &global_stats.block_volume, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);

  std::vector<double> product_sums(products.size(), 0.0);
  if (rank == 0) fs::create_directories(options.output_dir);
  MPI_Barrier(MPI_COMM_WORLD);
  for (std::size_t index = 0; index < products.size(); ++index) {
    const Product &product = products[index];
    const double *send = host_histogram.data() + product.device.offset;
    std::vector<double> global_values;
    if (rank == 0) global_values.resize(product.device.total_bins);
    MPI_Reduce(send, rank == 0 ? global_values.data() : nullptr,
               product.device.total_bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      double sum = 0.0;
      for (const double value : global_values) sum += value;
      product_sums[index] = sum;
      const double output_time =
          std::isfinite(options.output_time) ? options.output_time : header.time;
      write_product(options.output_dir, options.basename, options.output_number,
                    output_time, product, global_values);
      std::cout << "wrote " << product.id << " bins=" << product.device.total_bins
                << " sum=" << std::setprecision(12) << sum << "\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  const double elapsed = MPI_Wtime() - start_time;
  if (rank == 0) {
    write_manifest(options.output_dir, options, header, paths, products, global_stats,
                   elapsed, product_sums);
    std::cout << "processed shards=" << global_stats.shards
              << " blocks=" << global_stats.blocks
              << " cells=" << global_stats.cells
              << " payload_GiB=" << static_cast<double>(global_stats.bytes) /
                                        (1024.0 * 1024.0 * 1024.0)
              << " elapsed_s=" << elapsed << "\n";
  }
}

std::uint64_t parse_u64(const std::string &value, const std::string &option) {
  std::size_t parsed = 0;
  const std::uint64_t result = std::stoull(value, &parsed);
  if (parsed != value.size()) throw std::runtime_error("Invalid value for " + option);
  return result;
}

Options parse_options(int argc, char **argv) {
  Options options;
  const auto next = [&argc, &argv](int &index, const std::string &option) {
    if (index + 1 >= argc) throw std::runtime_error("Missing value for " + option);
    return std::string(argv[++index]);
  };
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--input-dir") {
      options.input_dir = next(i, argument);
    } else if (argument == "--output-dir") {
      options.output_dir = next(i, argument);
    } else if (argument == "--sequence") {
      options.sequence = next(i, argument);
    } else if (argument == "--output-number") {
      options.output_number = next(i, argument);
    } else if (argument == "--basename") {
      options.basename = next(i, argument);
    } else if (argument == "--products") {
      options.products = next(i, argument);
    } else if (argument == "--only-product") {
      options.only_products.push_back(next(i, argument));
    } else if (argument == "--chunk-blocks") {
      options.chunk_blocks = static_cast<int>(parse_u64(next(i, argument), argument));
    } else if (argument == "--max-shards") {
      options.max_shards = parse_u64(next(i, argument), argument);
    } else if (argument == "--expected-shards") {
      options.expected_shards = parse_u64(next(i, argument), argument);
    } else if (argument == "--max-blocks-per-shard") {
      options.max_blocks_per_shard = parse_u64(next(i, argument), argument);
    } else if (argument == "--output-time") {
      options.output_time = std::stod(next(i, argument));
    } else if (argument == "--dry-run") {
      options.dry_run = true;
    } else if (argument == "--list-products") {
      options.list_products = true;
    } else if (argument == "--help" || argument == "-h") {
      options.help = true;
    } else {
      throw std::runtime_error("Unknown option: " + argument);
    }
  }
  if (options.output_number.empty()) options.output_number = options.sequence;
  if (options.chunk_blocks <= 0) throw std::runtime_error("--chunk-blocks must be positive");
  return options;
}

void print_help() {
  std::cout
      << "Usage: gotham_pdf_rebuild --input-dir DIR --sequence XXXXX --output-dir DIR "
         "[options]\n"
      << "Options:\n"
      << "  --products original|science|all   Product set (default: all)\n"
      << "  --only-product ID                 Restrict to one selected product; repeatable\n"
      << "  --output-number XXXXX             PDF payload number (default: source sequence)\n"
      << "  --output-time VALUE               Exact original PDF payload time\n"
      << "  --chunk-blocks N                  Blocks per host-to-device transfer (default: 32)\n"
      << "  --expected-shards N               Require exact contiguous node_ shard manifest\n"
      << "  --max-shards N                    Process only the first N shards\n"
      << "  --max-blocks-per-shard N          Process only the first N blocks of each shard\n"
      << "  --dry-run                         Validate layout and print the execution plan\n"
      << "  --list-products                   Print selected product definitions and exit\n";
}

void print_products(const std::vector<Product> &products) {
  for (const auto &product : products) {
    std::cout << product.id << " [";
    for (std::size_t d = 0; d < product.axes.size(); ++d) {
      if (d > 0) std::cout << ", ";
      std::cout << variable_name(product.axes[d].variable) << ":"
                << product.axes[d].nbin << ":" << scale_name(product.axes[d].scale);
    }
    std::cout << "] weight=" << weight_name(product.weight)
              << " total_bins=" << product.device.total_bins << "\n";
  }
  std::cout << "total_bins=" << total_bins(products)
            << " histogram_GiB="
            << static_cast<double>(total_bins(products) * sizeof(double)) /
                   (1024.0 * 1024.0 * 1024.0)
            << "\n";
}

int run(int argc, char **argv, int rank, int nranks) {
  Options options = parse_options(argc, argv);
  if (options.help) {
    if (rank == 0) print_help();
    return 0;
  }
  const std::vector<Product> products =
      make_products(options.products, options.only_products);
  if (options.list_products) {
    if (rank == 0) print_products(products);
    return 0;
  }
  if (options.input_dir.empty() || options.sequence.empty() || options.output_dir.empty()) {
    throw std::runtime_error("--input-dir, --sequence, and --output-dir are required");
  }

  std::vector<fs::path> paths;
  if (rank == 0) {
    paths = find_shards(options.input_dir, options.sequence);
    if (options.max_shards > 0 && paths.size() > options.max_shards) {
      paths.resize(options.max_shards);
    }
    validate_shard_manifest(paths, options.expected_shards);
  }
  paths = broadcast_paths(std::move(paths), rank);

  HeaderWire reference{};
  if (rank == 0) reference = to_wire(read_header(paths.front()));
  MPI_Bcast(&reference, sizeof(reference), MPI_BYTE, 0, MPI_COMM_WORLD);
  validate_gotham_contract(reference);
  const bool calculate_cooling = needs_cooling_values(products);
  if (calculate_cooling && !reference.cgm_cooling) {
    throw std::runtime_error(
        "Cooling products require the embedded input deck to enable cgm_cooling");
  }

  if (rank == 0) {
    std::cout << "source_time=" << std::setprecision(17) << reference.time
              << " cycle=" << reference.cycle << " shards=" << paths.size()
              << " ranks=" << nranks << " chunk_blocks=" << options.chunk_blocks << "\n";
    print_products(products);
  }
  if (options.dry_run) return 0;

  Kokkos::View<DeviceProduct *, Kokkos::HostSpace> host_products("host_products",
                                                                  products.size());
  for (std::size_t i = 0; i < products.size(); ++i) host_products(i) = products[i].device;
  Kokkos::View<DeviceProduct *> device_products("device_products", products.size());
  Kokkos::deep_copy(device_products, host_products);
  Kokkos::View<double *> histogram("histogram", total_bins(products));
  Kokkos::deep_copy(histogram, 0.0);
  CoolingTables cooling_tables;
  if (calculate_cooling) cooling_tables = make_cooling_tables();
  Kokkos::fence();

  const double start_time = MPI_Wtime();
  std::vector<BlockKey> block_keys;
  const LocalStats local = process_shards(paths, options, rank, nranks, reference,
                                          device_products, histogram, block_keys,
                                          cooling_tables, calculate_cooling);
  if (options.max_shards == 0 && options.max_blocks_per_shard == 0) {
    validate_global_mesh(block_keys, local, reference, rank, nranks);
  } else if (rank == 0) {
    std::cout << "skipping global AMR closure for intentionally limited input\n";
  }
  reduce_and_write(options, reference, paths, products, histogram, local, start_time, rank);
  return 0;
}

}  // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rank = 0;
  int nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  int status = 0;
  try {
    status = run(argc, argv, rank, nranks);
  } catch (const std::exception &error) {
    std::cerr << "rank " << rank << " fatal: " << error.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
    status = 1;
  }
  Kokkos::finalize();
  MPI_Finalize();
  return status;
}
