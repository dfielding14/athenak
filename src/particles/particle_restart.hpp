#ifndef PARTICLES_PARTICLE_RESTART_HPP_
#define PARTICLES_PARTICLE_RESTART_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_restart.hpp
//! \brief canonical typed-v2 particle restart codec for relativistic CR tracers

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace particles {
namespace restart {

constexpr std::size_t kV2HeaderBytes = 144;
constexpr std::uint16_t kV2SchemaMajor = 2;
constexpr std::uint16_t kV2SchemaMinor = 0;
constexpr std::uint32_t kV2FlagsRelativisticHC = 1U;
constexpr std::uint32_t kV2EndianMarker = 0x01020304U;
constexpr std::uint32_t kV2IntegerFields = 3U;
constexpr std::uint32_t kV2RealFields = 22U;
constexpr std::uint32_t kV2RecordBytes =
    kV2IntegerFields*sizeof(std::int32_t) + kV2RealFields*sizeof(double);
constexpr std::uint32_t kV2RejectRankCountChange = 1U;
constexpr std::size_t kV2HeaderChecksumOffset = 120;
constexpr std::array<unsigned char, 8> kV2Magic =
    {'A', 'K', 'P', 'R', 'S', 'T', '2', '\0'};

struct V2Header {
  std::uint32_t flags = kV2FlagsRelativisticHC;
  std::uint32_t saved_nranks = 0;
  std::uint32_t saved_rank = 0;
  std::uint64_t local_count = 0;
  std::uint64_t global_count = 0;
  std::int64_t mesh_cycle = 0;
  double mesh_time = 0.0;
  double mesh_dt = 0.0;
  double particle_dtnew = 0.0;
  std::uint64_t checkpoint_id = 0;
  std::uint64_t payload_bytes = 0;
  std::uint64_t payload_checksum = 0;
  std::uint64_t header_checksum = 0;
  std::uint64_t config_fingerprint = 0;
};

struct ManifestShard {
  std::uint32_t rank = 0;
  std::string relative_path;
  std::uint64_t byte_count = 0;
  std::uint64_t local_count = 0;
  std::uint64_t payload_checksum = 0;
  std::uint64_t header_checksum = 0;
};

struct Manifest {
  std::uint64_t checkpoint_id = 0;
  std::uint32_t saved_nranks = 0;
  std::uint64_t global_count = 0;
  std::int64_t mesh_cycle = 0;
  double mesh_time = 0.0;
  double mesh_dt = 0.0;
  double particle_dtnew = 0.0;
  std::uint64_t config_fingerprint = 0;
  std::uint64_t mesh_byte_count = 0;
  std::uint64_t mesh_checksum = 0;
  std::uint64_t mesh_topology_hash = 0;
  std::string mesh_restart;
  std::map<std::uint32_t, ManifestShard> shards;
};

struct MeshWitness {
  std::uint64_t checkpoint_id = 0;
  std::uint32_t saved_nranks = 0;
  std::int64_t mesh_cycle = 0;
  double mesh_time = 0.0;
  double mesh_dt = 0.0;
  std::uint64_t mesh_byte_count = 0;
  std::uint64_t mesh_checksum = 0;
  std::uint64_t mesh_topology_hash = 0;
};

struct V2Shard {
  V2Header header;
  std::vector<std::int32_t> idata;
  std::vector<double> rdata;
  std::string shard_path;
  std::string manifest_path;
};

template <typename RealType>
struct LegacyV1Shard {
  RealType mesh_time = 0.0;
  RealType mesh_dt = 0.0;
  std::vector<RealType> rdata;
  std::string shard_path;
};

inline std::uint64_t UpdateFNV1a(std::uint64_t hash, const void *data,
                                std::size_t nbyte) {
  const unsigned char *bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i=0; i<nbyte; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= 1099511628211ULL;
  }
  return hash;
}

inline std::uint64_t FNV1a(const void *data, std::size_t nbyte) {
  return UpdateFNV1a(1469598103934665603ULL, data, nbyte);
}

inline void Require(bool condition, const std::string &message) {
  if (!condition) {throw std::runtime_error(message);}
}

inline std::uint64_t CheckedMultiply(std::uint64_t lhs, std::uint64_t rhs,
                                     const std::string &label) {
  if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max()/lhs) {
    throw std::runtime_error(label + " overflows uint64");
  }
  return lhs*rhs;
}

inline std::uint64_t ParseU64(const std::string &text, const std::string &label,
                              std::uint64_t maximum =
                                  std::numeric_limits<std::uint64_t>::max()) {
  Require(!text.empty() && text.front() != '-',
          label + " is not an unsigned integer");
  std::size_t used = 0;
  std::uint64_t value = 0;
  try {
    value = std::stoull(text, &used);
  } catch (const std::exception &) {
    throw std::runtime_error(label + " is not an unsigned integer");
  }
  Require(used == text.size() && value <= maximum,
          label + " is not a bounded unsigned integer");
  return value;
}

inline std::int64_t ParseI64(const std::string &text, const std::string &label) {
  std::size_t used = 0;
  std::int64_t value = 0;
  try {
    value = std::stoll(text, &used);
  } catch (const std::exception &) {
    throw std::runtime_error(label + " is not a signed integer");
  }
  Require(used == text.size(), label + " is not a complete signed integer");
  return value;
}

inline double ParseF64(const std::string &text, const std::string &label) {
  std::size_t used = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &used);
  } catch (const std::exception &) {
    throw std::runtime_error(label + " is not a floating-point value");
  }
  Require(used == text.size() && std::isfinite(value),
          label + " is not a finite floating-point value");
  return value;
}

inline void PutU16(std::vector<unsigned char> &bytes, std::size_t offset,
                   std::uint16_t value) {
  bytes.at(offset) = static_cast<unsigned char>(value & 0xffU);
  bytes.at(offset + 1) = static_cast<unsigned char>((value >> 8) & 0xffU);
}

inline void PutU32(std::vector<unsigned char> &bytes, std::size_t offset,
                   std::uint32_t value) {
  for (int i=0; i<4; ++i) {
    bytes.at(offset + i) = static_cast<unsigned char>((value >> (8*i)) & 0xffU);
  }
}

inline void PutU64(std::vector<unsigned char> &bytes, std::size_t offset,
                   std::uint64_t value) {
  for (int i=0; i<8; ++i) {
    bytes.at(offset + i) = static_cast<unsigned char>((value >> (8*i)) & 0xffU);
  }
}

inline void PutI32(std::vector<unsigned char> &bytes, std::size_t offset,
                   std::int32_t value) {
  PutU32(bytes, offset, static_cast<std::uint32_t>(value));
}

inline void PutI64(std::vector<unsigned char> &bytes, std::size_t offset,
                   std::int64_t value) {
  PutU64(bytes, offset, static_cast<std::uint64_t>(value));
}

inline void PutF64(std::vector<unsigned char> &bytes, std::size_t offset,
                   double value) {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value), "f64 requires 64-bit double");
  std::memcpy(&bits, &value, sizeof(bits));
  PutU64(bytes, offset, bits);
}

inline std::uint16_t GetU16(const std::vector<unsigned char> &bytes,
                            std::size_t offset) {
  Require(offset + 2 <= bytes.size(), "typed-v2 u16 read exceeds payload");
  return static_cast<std::uint16_t>(bytes[offset]) |
         static_cast<std::uint16_t>(bytes[offset + 1]) << 8;
}

inline std::uint32_t GetU32(const std::vector<unsigned char> &bytes,
                            std::size_t offset) {
  Require(offset + 4 <= bytes.size(), "typed-v2 u32 read exceeds payload");
  std::uint32_t value = 0;
  for (int i=0; i<4; ++i) {
    value |= static_cast<std::uint32_t>(bytes[offset + i]) << (8*i);
  }
  return value;
}

inline std::uint64_t GetU64(const std::vector<unsigned char> &bytes,
                            std::size_t offset) {
  Require(offset + 8 <= bytes.size(), "typed-v2 u64 read exceeds payload");
  std::uint64_t value = 0;
  for (int i=0; i<8; ++i) {
    value |= static_cast<std::uint64_t>(bytes[offset + i]) << (8*i);
  }
  return value;
}

inline std::int32_t GetI32(const std::vector<unsigned char> &bytes,
                           std::size_t offset) {
  return static_cast<std::int32_t>(GetU32(bytes, offset));
}

inline std::int64_t GetI64(const std::vector<unsigned char> &bytes,
                           std::size_t offset) {
  return static_cast<std::int64_t>(GetU64(bytes, offset));
}

inline double GetF64(const std::vector<unsigned char> &bytes, std::size_t offset) {
  std::uint64_t bits = GetU64(bytes, offset);
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

inline std::uint64_t ComputeCheckpointID(std::int64_t cycle, double time,
                                         double mesh_dt, int file_number,
                                         int nranks) {
  std::vector<unsigned char> bytes(32, 0);
  PutI64(bytes, 0, cycle);
  PutF64(bytes, 8, time);
  PutF64(bytes, 16, mesh_dt);
  PutI32(bytes, 24, file_number);
  PutI32(bytes, 28, nranks);
  return FNV1a(bytes.data(), bytes.size());
}

inline std::uint64_t ComputeRelativisticConfigFingerprint(
    double c_model, double alpha_s, std::uint32_t field_source,
    std::uint32_t temporal_sampling, bool subcycle, bool subcycle_strict,
    int subcycle_max_steps, double subcycle_cell_fraction,
    double subcycle_meshblock_fraction, double subcycle_gyro_fraction,
    double subcycle_electric_kick_max) {
  std::vector<unsigned char> bytes(68, 0);
  PutF64(bytes, 0, c_model);
  PutF64(bytes, 8, alpha_s);
  PutU32(bytes, 16, field_source);
  PutU32(bytes, 20, temporal_sampling);
  PutU32(bytes, 24, subcycle ? 1U : 0U);
  PutU32(bytes, 28, subcycle_strict ? 1U : 0U);
  PutI32(bytes, 32, subcycle_max_steps);
  PutF64(bytes, 36, subcycle_cell_fraction);
  PutF64(bytes, 44, subcycle_meshblock_fraction);
  PutF64(bytes, 52, subcycle_gyro_fraction);
  PutF64(bytes, 60, subcycle_electric_kick_max);
  return FNV1a(bytes.data(), bytes.size());
}

inline std::uint64_t ComputeTopologyHash(const int *rank_eachmb, int nmb_total) {
  Require(nmb_total >= 0, "mesh topology hash rejects negative MeshBlock count");
  std::vector<unsigned char> bytes(
      sizeof(std::uint32_t) + static_cast<std::size_t>(nmb_total)*sizeof(std::int32_t),
      0);
  PutU32(bytes, 0, static_cast<std::uint32_t>(nmb_total));
  for (int gid=0; gid<nmb_total; ++gid) {
    PutI32(bytes, sizeof(std::uint32_t) + gid*sizeof(std::int32_t), rank_eachmb[gid]);
  }
  return FNV1a(bytes.data(), bytes.size());
}

inline std::vector<unsigned char> EncodeHeader(V2Header header) {
  std::vector<unsigned char> bytes(kV2HeaderBytes, 0);
  std::copy(kV2Magic.begin(), kV2Magic.end(), bytes.begin());
  PutU16(bytes, 8, kV2SchemaMajor);
  PutU16(bytes, 10, kV2SchemaMinor);
  PutU32(bytes, 12, static_cast<std::uint32_t>(kV2HeaderBytes));
  PutU32(bytes, 16, header.flags);
  PutU32(bytes, 20, kV2EndianMarker);
  PutU32(bytes, 24, kV2IntegerFields);
  PutU32(bytes, 28, kV2RealFields);
  PutU32(bytes, 32, kV2RecordBytes);
  PutU32(bytes, 36, header.saved_nranks);
  PutU32(bytes, 40, header.saved_rank);
  PutU32(bytes, 44, kV2RejectRankCountChange);
  PutU64(bytes, 48, header.local_count);
  PutU64(bytes, 56, header.global_count);
  PutI64(bytes, 64, header.mesh_cycle);
  PutF64(bytes, 72, header.mesh_time);
  PutF64(bytes, 80, header.mesh_dt);
  PutF64(bytes, 88, header.particle_dtnew);
  PutU64(bytes, 96, header.checkpoint_id);
  PutU64(bytes, 104, header.payload_bytes);
  PutU64(bytes, 112, header.payload_checksum);
  PutU64(bytes, kV2HeaderChecksumOffset, 0);
  PutU32(bytes, 128, 1U);
  PutU64(bytes, 132, header.config_fingerprint);
  header.header_checksum = FNV1a(bytes.data(), bytes.size());
  PutU64(bytes, kV2HeaderChecksumOffset, header.header_checksum);
  return bytes;
}

inline V2Header DecodeHeader(const std::vector<unsigned char> &bytes) {
  Require(bytes.size() >= kV2HeaderBytes, "typed-v2 restart has truncated header");
  Require(std::equal(kV2Magic.begin(), kV2Magic.end(), bytes.begin()),
          "typed-v2 restart has invalid magic");
  Require(GetU16(bytes, 8) == kV2SchemaMajor &&
          GetU16(bytes, 10) == kV2SchemaMinor,
          "typed-v2 restart schema version is unsupported");
  Require(GetU32(bytes, 12) == kV2HeaderBytes,
          "typed-v2 restart header width is unsupported");
  Require(GetU32(bytes, 16) == kV2FlagsRelativisticHC,
          "typed-v2 restart mode flags are unsupported");
  Require(GetU32(bytes, 20) == kV2EndianMarker,
          "typed-v2 restart endian marker is invalid");
  Require(GetU32(bytes, 24) == kV2IntegerFields &&
          GetU32(bytes, 28) == kV2RealFields &&
          GetU32(bytes, 32) == kV2RecordBytes,
          "typed-v2 restart record layout is unsupported");
  Require(GetU32(bytes, 44) == kV2RejectRankCountChange,
          "typed-v2 restart topology policy is unsupported");
  Require(GetU32(bytes, 128) == 1U,
          "typed-v2 restart requires a manifest");
  std::vector<unsigned char> checksum_bytes(bytes.begin(),
                                            bytes.begin() + kV2HeaderBytes);
  std::uint64_t stored_header_checksum = GetU64(bytes, kV2HeaderChecksumOffset);
  PutU64(checksum_bytes, kV2HeaderChecksumOffset, 0);
  Require(FNV1a(checksum_bytes.data(), checksum_bytes.size()) ==
          stored_header_checksum, "typed-v2 restart header checksum mismatch");

  V2Header header;
  header.flags = GetU32(bytes, 16);
  header.saved_nranks = GetU32(bytes, 36);
  header.saved_rank = GetU32(bytes, 40);
  header.local_count = GetU64(bytes, 48);
  header.global_count = GetU64(bytes, 56);
  header.mesh_cycle = GetI64(bytes, 64);
  header.mesh_time = GetF64(bytes, 72);
  header.mesh_dt = GetF64(bytes, 80);
  header.particle_dtnew = GetF64(bytes, 88);
  header.checkpoint_id = GetU64(bytes, 96);
  header.payload_bytes = GetU64(bytes, 104);
  header.payload_checksum = GetU64(bytes, 112);
  header.header_checksum = stored_header_checksum;
  header.config_fingerprint = GetU64(bytes, 132);
  Require(std::isfinite(header.mesh_time) && std::isfinite(header.mesh_dt) &&
          std::isfinite(header.particle_dtnew),
          "typed-v2 restart header contains non-finite timestep metadata");
  Require(header.saved_nranks > 0 && header.saved_rank < header.saved_nranks,
          "typed-v2 restart has invalid saved rank topology");
  Require(header.local_count <=
          static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
          "typed-v2 restart local particle count exceeds runtime int range");
  Require(header.payload_bytes ==
          CheckedMultiply(header.local_count, kV2RecordBytes,
                          "typed-v2 payload byte count"),
          "typed-v2 restart payload byte count is inconsistent");
  return header;
}

inline std::vector<unsigned char> EncodePayload(const std::vector<std::int32_t> &idata,
                                                const std::vector<double> &rdata) {
  Require(idata.size() % kV2IntegerFields == 0,
          "typed-v2 writer integer record count is inconsistent");
  std::size_t count = idata.size()/kV2IntegerFields;
  Require(rdata.size() == count*kV2RealFields,
          "typed-v2 writer real record count is inconsistent");
  std::vector<unsigned char> payload(count*kV2RecordBytes, 0);
  for (std::size_t p=0; p<count; ++p) {
    std::size_t offset = p*kV2RecordBytes;
    for (std::size_t n=0; n<kV2IntegerFields; ++n) {
      PutI32(payload, offset + n*sizeof(std::int32_t),
             idata[p*kV2IntegerFields + n]);
    }
    offset += kV2IntegerFields*sizeof(std::int32_t);
    for (std::size_t n=0; n<kV2RealFields; ++n) {
      double value = rdata[p*kV2RealFields + n];
      Require(std::isfinite(value), "typed-v2 writer rejects non-finite particle data");
      PutF64(payload, offset + n*sizeof(double), value);
    }
  }
  return payload;
}

inline std::vector<unsigned char> EncodeShard(V2Header &header,
                                              const std::vector<std::int32_t> &idata,
                                              const std::vector<double> &rdata) {
  std::vector<unsigned char> payload = EncodePayload(idata, rdata);
  header.payload_bytes = payload.size();
  header.payload_checksum = FNV1a(payload.data(), payload.size());
  std::vector<unsigned char> bytes = EncodeHeader(header);
  header.header_checksum = GetU64(bytes, kV2HeaderChecksumOffset);
  bytes.insert(bytes.end(), payload.begin(), payload.end());
  return bytes;
}

inline std::string ReplaceRankDirectory(std::string path, int rank) {
  char rank_dir[20];
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", rank);
  std::size_t rank_pos = path.find("rank_00000000");
  Require(rank_pos != std::string::npos,
          "typed-v2 restart input must name the rank_00000000 shard");
  path.replace(rank_pos, sizeof("rank_00000000") - 1, rank_dir);
  return path;
}

inline std::string ManifestPathFromRankZeroShard(const std::string &path) {
  std::size_t rank_pos = path.find("rank_00000000/");
  Require(rank_pos != std::string::npos,
          "typed-v2 restart input must name a rank_00000000 shard");
  return path.substr(0, rank_pos) +
         path.substr(rank_pos + sizeof("rank_00000000/") - 1) + ".manifest";
}

inline bool IsSafeRelativeShardPath(const std::string &path) {
  return !path.empty() && path.front() != '/' && path.find("..") == std::string::npos &&
         path.find('\\') == std::string::npos;
}

inline std::map<std::string, std::string> ParseKeyValueManifest(
    const std::string &path, std::vector<ManifestShard> &shards) {
  std::ifstream input(path);
  Require(input.good(), "typed-v2 restart manifest '" + path + "' could not be opened");
  std::map<std::string, std::string> values;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {continue;}
    std::istringstream stream(line);
    std::string key;
    stream >> key;
    if (key == "shard") {
      ManifestShard shard;
      std::string rank, byte_count, local_count, payload_checksum, header_checksum;
      Require(static_cast<bool>(stream >> rank >> shard.relative_path >>
                  byte_count >> local_count >> payload_checksum >> header_checksum),
              "typed-v2 restart manifest has malformed shard row");
      std::string extra;
      Require(!(stream >> extra), "typed-v2 restart manifest shard row has extra data");
      shard.rank = static_cast<std::uint32_t>(
          ParseU64(rank, "typed-v2 restart manifest shard rank",
                   std::numeric_limits<std::uint32_t>::max()));
      shard.byte_count = ParseU64(byte_count,
                                  "typed-v2 restart manifest shard byte count");
      shard.local_count = ParseU64(local_count,
                                   "typed-v2 restart manifest shard local count");
      shard.payload_checksum = ParseU64(
          payload_checksum, "typed-v2 restart manifest shard payload checksum");
      shard.header_checksum = ParseU64(
          header_checksum, "typed-v2 restart manifest shard header checksum");
      shards.push_back(shard);
    } else {
      std::string value, extra;
      Require(static_cast<bool>(stream >> value) && !(stream >> extra),
              "typed-v2 restart manifest has malformed key '" + key + "'");
      Require(values.emplace(key, value).second,
              "typed-v2 restart manifest repeats key '" + key + "'");
    }
  }
  return values;
}

inline Manifest ReadManifest(const std::string &path) {
  std::vector<ManifestShard> shard_rows;
  std::map<std::string, std::string> values = ParseKeyValueManifest(path, shard_rows);
  auto Value = [&values](const std::string &key) {
    auto it = values.find(key);
    Require(it != values.end(), "typed-v2 restart manifest is missing '" + key + "'");
    return it->second;
  };
  Require(Value("magic") == "AKPRST-MANIFEST",
          "typed-v2 restart manifest has invalid magic");
  Require(Value("version") == "2.0",
          "typed-v2 restart manifest schema version is unsupported");
  Require(Value("topology_policy") == "reject_rank_count_change",
          "typed-v2 restart manifest topology policy is unsupported");
  Require(Value("paired_mesh_checkpoint") == "required",
          "typed-v2 restart manifest must require paired mesh checkpoint metadata");
  Manifest manifest;
  manifest.checkpoint_id = ParseU64(Value("checkpoint_id"), "manifest checkpoint_id");
  manifest.saved_nranks = static_cast<std::uint32_t>(
      ParseU64(Value("saved_nranks"), "manifest saved_nranks",
               std::numeric_limits<std::uint32_t>::max()));
  manifest.global_count = ParseU64(Value("global_count"), "manifest global_count");
  manifest.mesh_cycle = ParseI64(Value("mesh_cycle"), "manifest mesh_cycle");
  manifest.mesh_time = ParseF64(Value("mesh_time"), "manifest mesh_time");
  manifest.mesh_dt = ParseF64(Value("mesh_dt"), "manifest mesh_dt");
  manifest.particle_dtnew = ParseF64(Value("particle_dtnew"), "manifest particle_dtnew");
  manifest.config_fingerprint =
      ParseU64(Value("config_fingerprint"), "manifest config_fingerprint");
  manifest.mesh_byte_count =
      ParseU64(Value("mesh_byte_count"), "manifest mesh_byte_count");
  manifest.mesh_checksum = ParseU64(Value("mesh_checksum"), "manifest mesh_checksum");
  manifest.mesh_topology_hash =
      ParseU64(Value("mesh_topology_hash"), "manifest mesh_topology_hash");
  manifest.mesh_restart = Value("mesh_restart");
  Require(std::isfinite(manifest.mesh_time) && std::isfinite(manifest.mesh_dt) &&
          std::isfinite(manifest.particle_dtnew),
          "typed-v2 restart manifest contains non-finite timestep metadata");
  Require(manifest.saved_nranks > 0 && shard_rows.size() == manifest.saved_nranks,
          "typed-v2 restart manifest shard coverage is incomplete");
  Require(IsSafeRelativeShardPath(manifest.mesh_restart),
          "typed-v2 restart manifest mesh checkpoint path is unsafe");
  std::uint64_t shard_count_sum = 0;
  for (const ManifestShard &shard : shard_rows) {
    Require(shard.rank < manifest.saved_nranks,
            "typed-v2 restart manifest shard rank is out of range");
    Require(IsSafeRelativeShardPath(shard.relative_path),
            "typed-v2 restart manifest shard path is unsafe");
    Require(manifest.shards.emplace(shard.rank, shard).second,
            "typed-v2 restart manifest repeats a shard rank");
    Require(shard.local_count <=
            std::numeric_limits<std::uint64_t>::max() - shard_count_sum,
            "typed-v2 restart manifest shard count sum overflows");
    shard_count_sum += shard.local_count;
  }
  Require(shard_count_sum == manifest.global_count,
          "typed-v2 restart manifest shard counts do not sum to global count");
  return manifest;
}

inline std::vector<unsigned char> ReadBytes(const std::string &path) {
  std::ifstream input(path, std::ios::binary);
  Require(input.good(), "typed-v2 restart shard '" + path + "' could not be opened");
  input.seekg(0, std::ios::end);
  std::streamoff end = input.tellg();
  Require(end >= 0, "typed-v2 restart shard size could not be read");
  input.seekg(0, std::ios::beg);
  std::vector<unsigned char> bytes(static_cast<std::size_t>(end));
  if (!bytes.empty()) {
    input.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
    Require(input.good(), "typed-v2 restart shard read was incomplete");
  }
  return bytes;
}

inline int FileNumber(const std::string &path);

inline MeshWitness ReadMeshWitness(const std::string &mesh_restart) {
  std::string witness_path = mesh_restart + ".rmeta";
  std::vector<ManifestShard> unused_rows;
  std::map<std::string, std::string> values =
      ParseKeyValueManifest(witness_path, unused_rows);
  Require(unused_rows.empty(), "mesh restart witness contains unexpected shard rows");
  auto Value = [&values](const std::string &key) {
    auto it = values.find(key);
    Require(it != values.end(), "mesh restart witness is missing '" + key + "'");
    return it->second;
  };
  Require(Value("magic") == "AKRST-WITNESS", "mesh restart witness has invalid magic");
  Require(Value("version") == "1", "mesh restart witness schema version is unsupported");
  MeshWitness witness;
  witness.checkpoint_id = ParseU64(Value("checkpoint_id"), "mesh witness checkpoint_id");
  witness.saved_nranks = static_cast<std::uint32_t>(
      ParseU64(Value("saved_nranks"), "mesh witness saved_nranks",
               std::numeric_limits<std::uint32_t>::max()));
  witness.mesh_cycle = ParseI64(Value("mesh_cycle"), "mesh witness mesh_cycle");
  witness.mesh_time = ParseF64(Value("mesh_time"), "mesh witness mesh_time");
  witness.mesh_dt = ParseF64(Value("mesh_dt"), "mesh witness mesh_dt");
  witness.mesh_byte_count =
      ParseU64(Value("mesh_byte_count"), "mesh witness mesh_byte_count");
  witness.mesh_checksum =
      ParseU64(Value("mesh_checksum"), "mesh witness mesh_checksum");
  witness.mesh_topology_hash =
      ParseU64(Value("mesh_topology_hash"), "mesh witness mesh_topology_hash");
  Require(std::isfinite(witness.mesh_time) && std::isfinite(witness.mesh_dt),
          "mesh restart witness contains non-finite timestep metadata");
  Require(witness.saved_nranks > 0 &&
          witness.saved_nranks <=
              static_cast<std::uint32_t>(std::numeric_limits<int>::max()),
          "mesh restart witness saved rank count exceeds runtime int range");
  std::vector<unsigned char> bytes = ReadBytes(mesh_restart);
  Require(bytes.size() == witness.mesh_byte_count,
          "mesh restart witness byte count does not match mesh checkpoint");
  Require(FNV1a(bytes.data(), bytes.size()) == witness.mesh_checksum,
          "mesh restart witness checksum does not match mesh checkpoint");
  Require(witness.checkpoint_id == ComputeCheckpointID(
              witness.mesh_cycle, witness.mesh_time, witness.mesh_dt,
              FileNumber(mesh_restart), witness.saved_nranks),
          "mesh restart witness checkpoint ID is inconsistent");
  return witness;
}

inline std::string Basename(const std::string &path) {
  std::size_t slash = path.rfind('/');
  return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

inline int FileNumber(const std::string &path) {
  std::string basename = Basename(path);
  std::size_t suffix = basename.rfind('.');
  Require(suffix != std::string::npos, "restart filename is missing suffix");
  std::size_t number_dot = basename.rfind('.', suffix - 1);
  Require(number_dot != std::string::npos && suffix == number_dot + 6,
          "restart filename is missing five-digit checkpoint number");
  for (std::size_t n=number_dot + 1; n<suffix; ++n) {
    Require(basename[n] >= '0' && basename[n] <= '9',
            "restart filename checkpoint number is malformed");
  }
  return std::stoi(basename.substr(number_dot + 1, 5));
}

inline std::string ReplaceLegacyRankDirectory(std::string path, int rank) {
  char rank_dir[20];
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", rank);
  std::size_t rank_pos = path.find("rank_00000000");
  if (rank_pos != std::string::npos) {
    path.replace(rank_pos, sizeof("rank_00000000") - 1, rank_dir);
  }
  return path;
}

template <typename RealType>
inline LegacyV1Shard<RealType> ReadLegacyV1Shard(const std::string &rank_zero_path,
                                                 int current_rank) {
  LegacyV1Shard<RealType> shard;
  shard.shard_path = ReplaceLegacyRankDirectory(rank_zero_path, current_rank);
  std::vector<unsigned char> bytes = ReadBytes(shard.shard_path);
  Require(bytes.size() >= 3*sizeof(RealType),
          "legacy particle restart has an incomplete header");
  RealType header[3];
  std::memcpy(header, bytes.data(), sizeof(header));
  Require(std::isfinite(header[0]) && std::isfinite(header[1]) &&
          std::isfinite(header[2]) && header[2] >= 0.0 &&
          header[2] <= static_cast<RealType>(std::numeric_limits<int>::max()) &&
          std::floor(header[2]) == header[2],
          "legacy particle restart has invalid header values");
  std::uint64_t count = static_cast<std::uint64_t>(header[2]);
  std::uint64_t expected_bytes =
      CheckedMultiply(3U + CheckedMultiply(17U, count,
                                          "legacy particle restart record count"),
                      sizeof(RealType), "legacy particle restart byte count");
  Require(expected_bytes == bytes.size(),
          "legacy particle restart byte count is inconsistent");
  shard.mesh_time = header[0];
  shard.mesh_dt = header[1];
  shard.rdata.resize(static_cast<std::size_t>(17U*count));
  if (!shard.rdata.empty()) {
    std::memcpy(shard.rdata.data(), bytes.data() + 3*sizeof(RealType),
                shard.rdata.size()*sizeof(RealType));
  }
  for (std::uint64_t p=0; p<count; ++p) {
    for (int n=0; n<17; ++n) {
      Require(std::isfinite(shard.rdata[17*p + n]),
              "legacy particle restart contains non-finite particle data");
    }
    for (int n=0; n<3; ++n) {
      RealType value = shard.rdata[17*p + n];
      Require(value >= static_cast<RealType>(std::numeric_limits<int>::min()) &&
              value <= static_cast<RealType>(std::numeric_limits<int>::max()) &&
              std::floor(value) == value,
              "legacy particle restart contains an invalid integer identifier");
    }
  }

  std::ifstream meta_probe(shard.shard_path + ".pmeta");
  if (meta_probe.good()) {
    meta_probe.close();
    std::vector<ManifestShard> unused_rows;
    std::map<std::string, std::string> values =
        ParseKeyValueManifest(shard.shard_path + ".pmeta", unused_rows);
    Require(unused_rows.empty(),
            "legacy particle restart metadata contains unexpected shard rows");
    auto Value = [&values](const std::string &key) {
      auto it = values.find(key);
      Require(it != values.end(),
              "legacy particle restart metadata is missing '" + key + "'");
      return it->second;
    };
    Require(Value("magic") == "AKPRST",
            "legacy particle restart metadata has invalid magic");
    Require(Value("version") == "1",
            "legacy particle restart metadata schema version is unsupported");
    Require(Value("binary_payload") == "legacy_prst",
            "legacy particle restart metadata payload kind is unsupported");
    Require(ParseU64(Value("real_size"), "legacy metadata real_size") == sizeof(RealType),
            "legacy particle restart metadata Real size is incompatible");
    Require(ParseU64(Value("record_real_count"),
                     "legacy metadata record_real_count") == 17U,
            "legacy particle restart metadata record width is incompatible");
    Require(ParseU64(Value("particle_count"), "legacy metadata particle_count") == count,
            "legacy particle restart metadata particle count is inconsistent");
    Require(ParseU64(Value("byte_count"), "legacy metadata byte_count") == bytes.size(),
            "legacy particle restart metadata byte count is inconsistent");
    Require(ParseU64(Value("checksum_fnv1a64"), "legacy metadata checksum_fnv1a64") ==
            FNV1a(bytes.data(), bytes.size()),
            "legacy particle restart metadata checksum mismatch");
  }
  return shard;
}

inline V2Shard ReadV2Shard(const std::string &rank_zero_path,
                           const std::string &paired_mesh_restart,
                           int current_rank, int current_nranks,
                           std::uint64_t current_topology_hash) {
  V2Shard shard;
  int particle_file_number = FileNumber(rank_zero_path);
  int mesh_file_number = FileNumber(paired_mesh_restart);
  Require(particle_file_number == mesh_file_number,
          "typed-v2 restart particle and mesh checkpoint numbers do not match");
  shard.manifest_path = ManifestPathFromRankZeroShard(rank_zero_path);
  Manifest manifest = ReadManifest(shard.manifest_path);
  Require(Basename(paired_mesh_restart) == Basename(manifest.mesh_restart),
          "typed-v2 restart particle manifest does not match paired mesh checkpoint");
  MeshWitness witness = ReadMeshWitness(paired_mesh_restart);
  Require(witness.checkpoint_id == manifest.checkpoint_id &&
          witness.saved_nranks == manifest.saved_nranks &&
          witness.mesh_cycle == manifest.mesh_cycle &&
          witness.mesh_time == manifest.mesh_time &&
          witness.mesh_dt == manifest.mesh_dt &&
          witness.mesh_byte_count == manifest.mesh_byte_count &&
          witness.mesh_checksum == manifest.mesh_checksum &&
          witness.mesh_topology_hash == manifest.mesh_topology_hash,
          "typed-v2 restart particle manifest does not match mesh witness");
  Require(current_nranks >= 0 &&
          manifest.saved_nranks == static_cast<std::uint32_t>(current_nranks),
          "typed-v2 restart rejects MPI rank-count change");
  Require(witness.mesh_topology_hash == current_topology_hash,
          "typed-v2 restart rejects changed rank-to-MeshBlock topology");
  Require(current_rank >= 0 &&
          static_cast<std::uint32_t>(current_rank) < manifest.saved_nranks,
          "typed-v2 restart current rank is out of range");
  const ManifestShard &manifest_shard =
      manifest.shards.at(static_cast<std::uint32_t>(current_rank));
  shard.shard_path = ReplaceRankDirectory(rank_zero_path, current_rank);
  std::size_t slash = shard.shard_path.rfind('/');
  std::size_t rank_slash = shard.shard_path.rfind('/', slash - 1);
  Require(slash != std::string::npos && rank_slash != std::string::npos,
          "typed-v2 restart shard path is malformed");
  std::string relative = shard.shard_path.substr(rank_slash + 1);
  Require(relative == manifest_shard.relative_path,
          "typed-v2 restart manifest shard path does not match requested shard");
  std::vector<unsigned char> bytes = ReadBytes(shard.shard_path);
  Require(bytes.size() == manifest_shard.byte_count,
          "typed-v2 restart shard byte count does not match manifest");
  shard.header = DecodeHeader(bytes);
  Require(shard.header.saved_rank == static_cast<std::uint32_t>(current_rank) &&
          shard.header.saved_nranks == manifest.saved_nranks &&
          shard.header.local_count == manifest_shard.local_count &&
          shard.header.global_count == manifest.global_count &&
          shard.header.mesh_cycle == manifest.mesh_cycle &&
          shard.header.mesh_time == manifest.mesh_time &&
          shard.header.mesh_dt == manifest.mesh_dt &&
          shard.header.particle_dtnew == manifest.particle_dtnew &&
          shard.header.checkpoint_id == manifest.checkpoint_id &&
          shard.header.config_fingerprint == manifest.config_fingerprint &&
          shard.header.payload_checksum == manifest_shard.payload_checksum &&
          shard.header.header_checksum == manifest_shard.header_checksum,
          "typed-v2 restart shard header does not match manifest");
  Require(shard.header.checkpoint_id == ComputeCheckpointID(
              shard.header.mesh_cycle, shard.header.mesh_time, shard.header.mesh_dt,
              mesh_file_number, current_nranks),
          "typed-v2 restart checkpoint ID does not match paired mesh checkpoint");
  Require(bytes.size() == kV2HeaderBytes + shard.header.payload_bytes,
          "typed-v2 restart shard byte count is inconsistent");
  const unsigned char *payload = bytes.data() + kV2HeaderBytes;
  Require(FNV1a(payload, shard.header.payload_bytes) == shard.header.payload_checksum,
          "typed-v2 restart payload checksum mismatch");
  shard.idata.resize(shard.header.local_count*kV2IntegerFields);
  shard.rdata.resize(shard.header.local_count*kV2RealFields);
  for (std::size_t p=0; p<shard.header.local_count; ++p) {
    std::size_t offset = kV2HeaderBytes + p*kV2RecordBytes;
    for (std::size_t n=0; n<kV2IntegerFields; ++n) {
      shard.idata[p*kV2IntegerFields + n] =
          GetI32(bytes, offset + n*sizeof(std::int32_t));
    }
    offset += kV2IntegerFields*sizeof(std::int32_t);
    for (std::size_t n=0; n<kV2RealFields; ++n) {
      double value = GetF64(bytes, offset + n*sizeof(double));
      Require(std::isfinite(value), "typed-v2 restart contains non-finite particle data");
      shard.rdata[p*kV2RealFields + n] = value;
    }
  }
  return shard;
}

} // namespace restart
} // namespace particles

#endif // PARTICLES_PARTICLE_RESTART_HPP_
