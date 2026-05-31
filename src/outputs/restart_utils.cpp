//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_utils.cpp
//! \brief Atomic-publication and integrity helpers for restart artifacts.

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "restart_utils.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace restart_utils {

[[noreturn]] void AbortOnFatalError() {
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {
    MPI_Finalized(&finalized);
  }
  if (initialized != 0 && finalized == 0) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#endif
  std::exit(EXIT_FAILURE);
}

namespace {

constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr const char *kMarkerMagic = "ATHENAK_RESTART_COMPLETE_V1";
constexpr const char *kManifestSchema = "ATHENAK_RESTART_MANIFEST_V1";

namespace fs = std::filesystem;

struct ManifestMember {
  std::string declared_path;
  fs::path resolved_path;
  FileDigest digest;
};

struct RestartManifest {
  std::vector<ManifestMember> members;
};

[[noreturn]] void Fail(const std::string &message) {
  std::cerr << "### FATAL ERROR in " << __FILE__ << ": " << message << std::endl;
  AbortOnFatalError();
}

std::string DirectoryName(const std::string &path) {
  const std::size_t slash = path.rfind('/');
  return (slash == std::string::npos) ? "." : path.substr(0, slash);
}

void SyncDirectory(const std::string &path) {
  const std::string dir = DirectoryName(path);
  const int fd = open(dir.c_str(), O_RDONLY | O_DIRECTORY);
  if (fd < 0) {
    Fail("unable to open restart directory '" + dir + "' for sync: " +
         std::strerror(errno));
  }
  if (fsync(fd) != 0) {
    const std::string error = std::strerror(errno);
    close(fd);
    Fail("unable to sync restart directory '" + dir + "': " + error);
  }
  if (close(fd) != 0) {
    Fail("unable to close restart directory '" + dir + "': " +
         std::strerror(errno));
  }
}

void WriteAtomicText(const std::string &path, const std::string &contents) {
  const std::string partial = path + ".partial";
  FILE *file = std::fopen(partial.c_str(), "wb");
  if (file == nullptr) {
    Fail("unable to open restart metadata partial '" + partial + "': " +
         std::strerror(errno));
  }
  const std::size_t written = std::fwrite(contents.data(), 1, contents.size(), file);
  if (written != contents.size() || std::fflush(file) != 0 || fsync(fileno(file)) != 0) {
    const std::string error = std::strerror(errno);
    std::fclose(file);
    Fail("unable to write restart metadata partial '" + partial + "': " + error);
  }
  if (std::fclose(file) != 0) {
    Fail("unable to close restart metadata partial '" + partial + "': " +
         std::strerror(errno));
  }
  PublishRestartArtifact(partial, path);
}

std::string HexDigest(const std::uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

bool EndsWith(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool DigestsMatch(const FileDigest &left, const FileDigest &right) {
  return left.size == right.size && left.fnv1a64 == right.fnv1a64;
}

bool ParseUnsigned(const std::string &text, const int base,
                   std::uint64_t &value) {
  if (text.empty()) return false;
  for (const char next : text) {
    int digit = -1;
    if (next >= '0' && next <= '9') {
      digit = next - '0';
    } else if (next >= 'a' && next <= 'f') {
      digit = next - 'a' + 10;
    } else if (next >= 'A' && next <= 'F') {
      digit = next - 'A' + 10;
    }
    if (digit < 0 || digit >= base) return false;
  }
  std::size_t consumed = 0;
  try {
    value = std::stoull(text, &consumed, base);
  } catch (...) {
    return false;
  }
  return consumed == text.size();
}

int HexDigit(const char value) {
  if (value >= '0' && value <= '9') return value - '0';
  if (value >= 'a' && value <= 'f') return value - 'a' + 10;
  if (value >= 'A' && value <= 'F') return value - 'A' + 10;
  return -1;
}

std::string JsonEscape(const std::string &value) {
  std::ostringstream escaped;
  for (const unsigned char next : value) {
    switch (next) {
      case '"': escaped << "\\\""; break;
      case '\\': escaped << "\\\\"; break;
      case '\b': escaped << "\\b"; break;
      case '\f': escaped << "\\f"; break;
      case '\n': escaped << "\\n"; break;
      case '\r': escaped << "\\r"; break;
      case '\t': escaped << "\\t"; break;
      default:
        if (next < 0x20) {
          escaped << "\\u00" << std::hex << std::setfill('0') << std::setw(2)
                  << static_cast<int>(next);
        } else {
          escaped << next;
        }
    }
  }
  return escaped.str();
}

class ManifestParser {
 public:
  explicit ManifestParser(const std::string &text) : text_(text) {}

  bool Parse(RestartManifest &manifest, std::string &error) {
    bool have_schema = false;
    bool have_members = false;
    if (!Expect('{')) return ReturnError(error);
    while (true) {
      if (Take('}')) break;
      std::string key;
      if (!ParseString(key) || !Expect(':')) return ReturnError(error);
      if (key == "schema") {
        std::string schema;
        if (have_schema) return Fail("duplicate schema field", error);
        if (!ParseString(schema)) return ReturnError(error);
        if (schema != kManifestSchema) {
          return Fail("unsupported schema '" + schema + "'", error);
        }
        have_schema = true;
      } else if (key == "members") {
        if (have_members) return Fail("duplicate members field", error);
        if (!ParseMembers(manifest.members)) return ReturnError(error);
        have_members = true;
      } else {
        return Fail("unknown top-level field '" + key + "'", error);
      }
      if (Take('}')) break;
      if (!Expect(',')) return ReturnError(error);
    }
    SkipWhitespace();
    if (position_ != text_.size()) {
      return Fail("trailing content", error);
    }
    if (!have_schema || !have_members) {
      return Fail("schema or members field is missing", error);
    }
    if (manifest.members.empty()) {
      return Fail("members list is empty", error);
    }
    return true;
  }

 private:
  void SkipWhitespace() {
    while (position_ < text_.size() &&
           std::isspace(static_cast<unsigned char>(text_[position_]))) {
      ++position_;
    }
  }

  bool Take(const char token) {
    SkipWhitespace();
    if (position_ >= text_.size() || text_[position_] != token) return false;
    ++position_;
    return true;
  }

  bool Expect(const char token) {
    if (Take(token)) return true;
    std::ostringstream message;
    message << "expected '" << token << "'";
    return SetError(message.str());
  }

  bool ParseString(std::string &value) {
    SkipWhitespace();
    if (position_ >= text_.size() || text_[position_] != '"') {
      return SetError("expected string");
    }
    ++position_;
    while (position_ < text_.size()) {
      const char next = text_[position_++];
      if (next == '"') return true;
      if (static_cast<unsigned char>(next) < 0x20) {
        return SetError("unescaped control byte in string");
      }
      if (next != '\\') {
        value.push_back(next);
        continue;
      }
      if (position_ >= text_.size()) return SetError("unterminated escape");
      const char escaped = text_[position_++];
      switch (escaped) {
        case '"': value.push_back('"'); break;
        case '\\': value.push_back('\\'); break;
        case '/': value.push_back('/'); break;
        case 'b': value.push_back('\b'); break;
        case 'f': value.push_back('\f'); break;
        case 'n': value.push_back('\n'); break;
        case 'r': value.push_back('\r'); break;
        case 't': value.push_back('\t'); break;
        case 'u': {
          if (position_ + 4 > text_.size()) {
            return SetError("unterminated unicode escape");
          }
          int decoded = 0;
          for (int i = 0; i < 4; ++i) {
            const int digit = HexDigit(text_[position_++]);
            if (digit < 0) return SetError("invalid unicode escape");
            decoded = 16*decoded + digit;
          }
          if (decoded > 0x7f) return SetError("non-ASCII unicode escape");
          value.push_back(static_cast<char>(decoded));
          break;
        }
        default: return SetError("unsupported string escape");
      }
    }
    return SetError("unterminated string");
  }

  bool ParseUInt64(std::uint64_t &value) {
    SkipWhitespace();
    const std::size_t begin = position_;
    while (position_ < text_.size() &&
           std::isdigit(static_cast<unsigned char>(text_[position_]))) {
      ++position_;
    }
    if (begin == position_) return SetError("expected unsigned integer");
    if (!ParseUnsigned(text_.substr(begin, position_ - begin), 10, value)) {
      return SetError("invalid unsigned integer");
    }
    return true;
  }

  bool ParseDigest(std::uint64_t &value) {
    std::string digest;
    if (!ParseString(digest)) return false;
    if (digest.size() != 16 || !ParseUnsigned(digest, 16, value)) {
      return SetError("fnv1a64 must contain exactly 16 hexadecimal digits");
    }
    return true;
  }

  bool ParseMember(ManifestMember &member) {
    bool have_path = false;
    bool have_size = false;
    bool have_digest = false;
    if (!Expect('{')) return false;
    while (true) {
      if (Take('}')) break;
      std::string key;
      if (!ParseString(key) || !Expect(':')) return false;
      if (key == "path") {
        if (have_path) return SetError("duplicate member path field");
        if (!ParseString(member.declared_path)) return false;
        have_path = true;
      } else if (key == "size") {
        if (have_size) return SetError("duplicate member size field");
        if (!ParseUInt64(member.digest.size)) return false;
        have_size = true;
      } else if (key == "fnv1a64") {
        if (have_digest) return SetError("duplicate member fnv1a64 field");
        if (!ParseDigest(member.digest.fnv1a64)) return false;
        have_digest = true;
      } else {
        return SetError("unknown member field '" + key + "'");
      }
      if (Take('}')) break;
      if (!Expect(',')) return false;
    }
    if (!have_path || !have_size || !have_digest) {
      return SetError("member path, size, or fnv1a64 field is missing");
    }
    return true;
  }

  bool ParseMembers(std::vector<ManifestMember> &members) {
    if (!Expect('[')) return false;
    while (true) {
      if (Take(']')) break;
      ManifestMember member;
      if (!ParseMember(member)) return false;
      members.push_back(member);
      if (Take(']')) break;
      if (!Expect(',')) return false;
    }
    return true;
  }

  bool SetError(const std::string &message) {
    std::ostringstream out;
    out << message << " at byte " << position_;
    error_ = out.str();
    return false;
  }

  bool Fail(const std::string &message, std::string &error) {
    SetError(message);
    return ReturnError(error);
  }

  bool ReturnError(std::string &error) {
    error = error_;
    return false;
  }

  const std::string &text_;
  std::size_t position_ = 0;
  std::string error_;
};

bool IsRankDirectory(const fs::path &path, int &rank) {
  const std::string name = path.filename().string();
  if (name.size() != 13 || name.compare(0, 5, "rank_") != 0) return false;
  const std::string digits = name.substr(5);
  for (const char digit : digits) {
    if (!std::isdigit(static_cast<unsigned char>(digit))) return false;
  }
  std::uint64_t parsed = 0;
  if (!ParseUnsigned(digits, 10, parsed) ||
      parsed > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    return false;
  }
  rank = static_cast<int>(parsed);
  return true;
}

bool ResolveManifestMembers(const std::string &manifest_path,
                            RestartManifest &manifest,
                            std::string &error) {
  const fs::path normalized_manifest = fs::path(manifest_path).lexically_normal();
  const fs::path rst_dir = normalized_manifest.parent_path();
  if (rst_dir.filename() != "rst") {
    error = "restart manifest is not located in an rst directory: " + manifest_path;
    return false;
  }
  const fs::path root = rst_dir.parent_path();
  std::set<std::string> unique_paths;
  for (auto &member : manifest.members) {
    const fs::path declared(member.declared_path);
    const fs::path normalized = declared.lexically_normal();
    if (declared.empty() || declared.is_absolute() ||
        normalized.generic_string() != declared.generic_string()) {
      error = "restart manifest member path is not canonical: " +
              member.declared_path;
      return false;
    }
    for (const auto &component : normalized) {
      if (component == "..") {
        error = "restart manifest member path escapes publication root: " +
                member.declared_path;
        return false;
      }
    }
    if (normalized.begin() == normalized.end() || *normalized.begin() != "rst") {
      error = "restart manifest member path is outside rst/: " + member.declared_path;
      return false;
    }
    member.resolved_path = (root / normalized).lexically_normal();
    if (!unique_paths.insert(member.resolved_path.generic_string()).second) {
      error = "restart manifest contains duplicate member path: " +
              member.declared_path;
      return false;
    }
  }
  return true;
}

bool ValidateManifestLayout(const std::string &manifest_path,
                            const RestartManifest &manifest,
                            std::string &error) {
  const fs::path manifest_file(manifest_path);
  std::string payload_name = manifest_file.filename().string();
  if (!EndsWith(payload_name, ".manifest")) {
    error = "restart manifest path lacks .manifest suffix: " + manifest_path;
    return false;
  }
  payload_name.resize(payload_name.size() - std::string(".manifest").size());

  bool ranked = false;
  std::set<int> ranks;
  for (const auto &member : manifest.members) {
    if (member.resolved_path.filename() != payload_name) {
      error = "restart manifest member filename does not match publication: " +
              member.declared_path;
      return false;
    }
    int rank = -1;
    const bool member_ranked = IsRankDirectory(member.resolved_path.parent_path(), rank);
    if (&member != &manifest.members.front() && member_ranked != ranked) {
      error = "restart manifest mixes shared and per-rank members: " + manifest_path;
      return false;
    }
    ranked = member_ranked;
    if (ranked) ranks.insert(rank);
  }
  if (!ranked && manifest.members.size() != 1) {
    error = "shared restart manifest must contain exactly one member: " + manifest_path;
    return false;
  }
  if (ranked) {
    if (ranks.size() != manifest.members.size()) {
      error = "restart manifest contains duplicate rank members: " + manifest_path;
      return false;
    }
    int expected = 0;
    for (const int rank : ranks) {
      if (rank != expected++) {
        error = "restart manifest rank members are not contiguous: " + manifest_path;
        return false;
      }
    }
  }
  return true;
}

bool ParseRestartManifest(const std::string &manifest_path,
                          RestartManifest &manifest,
                          std::string &error) {
  std::ifstream file(manifest_path);
  if (!file) {
    error = "completed restart artifact is missing: " + manifest_path;
    return false;
  }
  std::ostringstream contents;
  contents << file.rdbuf();
  const std::string manifest_text = contents.str();
  ManifestParser parser(manifest_text);
  std::string parser_error;
  if (!parser.Parse(manifest, parser_error)) {
    error = "restart manifest is malformed: " + manifest_path + ": " + parser_error;
    return false;
  }
  return ResolveManifestMembers(manifest_path, manifest, error) &&
         ValidateManifestLayout(manifest_path, manifest, error);
}

fs::path ManifestPathForMember(const fs::path &artifact_path) {
  const fs::path normalized = artifact_path.lexically_normal();
  const fs::path parent = normalized.parent_path();
  int rank = -1;
  if (IsRankDirectory(parent, rank)) {
    return parent.parent_path() / (normalized.filename().string() + ".manifest");
  }
  return fs::path(normalized.string() + ".manifest");
}

bool VerifyCompletedArtifact(const std::string &artifact_path,
                             FileDigest &actual,
                             std::string &error) {
  std::ifstream artifact(artifact_path, std::ios::binary);
  if (!artifact) {
    error = "completed restart artifact is missing: " + artifact_path;
    return false;
  }
  artifact.close();
  const std::string marker_path = artifact_path + ".complete";
  std::ifstream marker(marker_path);
  if (!marker) {
    error = "restart completion marker is missing: " + marker_path;
    return false;
  }
  std::string magic;
  std::string size_line;
  std::string hash_line;
  std::string trailing;
  if (!std::getline(marker, magic) || !std::getline(marker, size_line) ||
      !std::getline(marker, hash_line) || std::getline(marker, trailing) ||
      magic != kMarkerMagic || size_line.rfind("size=", 0) != 0 ||
      hash_line.rfind("fnv1a64=", 0) != 0) {
    error = "restart completion marker is malformed: " + marker_path;
    return false;
  }
  FileDigest expected;
  if (!ParseUnsigned(size_line.substr(5), 10, expected.size) ||
      !ParseUnsigned(hash_line.substr(8), 16, expected.fnv1a64)) {
    error = "restart completion marker has invalid values: " + marker_path;
    return false;
  }
  actual = ComputeFileDigest(artifact_path);
  if (!DigestsMatch(actual, expected)) {
    error = "restart checksum mismatch for completed artifact: " + artifact_path;
    return false;
  }
  return true;
}

}  // namespace

FileDigest ComputeFileDigest(const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    Fail("unable to open restart artifact '" + path + "' for checksum");
  }
  FileDigest digest;
  digest.fnv1a64 = kFnvOffset;
  char buffer[1024 * 1024];
  while (file) {
    file.read(buffer, sizeof(buffer));
    const std::streamsize count = file.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      digest.fnv1a64 ^= static_cast<unsigned char>(buffer[i]);
      digest.fnv1a64 *= kFnvPrime;
    }
    digest.size += static_cast<std::uint64_t>(count);
  }
  if (!file.eof()) {
    Fail("unable to read restart artifact '" + path + "' for checksum");
  }
  return digest;
}

void PublishRestartArtifact(const std::string &partial_path,
                            const std::string &final_path) {
  if (std::rename(partial_path.c_str(), final_path.c_str()) != 0) {
    Fail("unable to publish restart artifact '" + final_path + "': " +
         std::strerror(errno));
  }
  SyncDirectory(final_path);
}

void WriteCompletionMarker(const std::string &artifact_path,
                           const FileDigest &digest) {
  std::ostringstream marker;
  marker << kMarkerMagic << "\n"
         << "size=" << digest.size << "\n"
         << "fnv1a64=" << HexDigest(digest.fnv1a64) << "\n";
  WriteAtomicText(artifact_path + ".complete", marker.str());
}

void WriteRestartManifest(
    const std::string &manifest_path,
    const std::vector<std::pair<std::string, FileDigest>> &members) {
  std::ostringstream manifest;
  manifest << "{\n  \"schema\": \"ATHENAK_RESTART_MANIFEST_V1\",\n"
           << "  \"members\": [\n";
  for (std::size_t i = 0; i < members.size(); ++i) {
    manifest << "    {\"path\": \"" << JsonEscape(members[i].first) << "\", \"size\": "
             << members[i].second.size << ", \"fnv1a64\": \""
             << HexDigest(members[i].second.fnv1a64) << "\"}"
             << ((i + 1 == members.size()) ? "\n" : ",\n");
  }
  manifest << "  ]\n}\n";
  WriteAtomicText(manifest_path, manifest.str());
  WriteCompletionMarker(manifest_path, ComputeFileDigest(manifest_path));
}

bool VerifyRestartArtifact(const std::string &artifact_path, std::string &error) {
  FileDigest actual;
  if (!VerifyCompletedArtifact(artifact_path, actual, error)) return false;

  if (EndsWith(artifact_path, ".manifest")) {
    RestartManifest manifest;
    return ParseRestartManifest(artifact_path, manifest, error);
  }

  const fs::path normalized_artifact = fs::path(artifact_path).lexically_normal();
  const std::string manifest_path =
      ManifestPathForMember(normalized_artifact).generic_string();
  FileDigest manifest_digest;
  if (!VerifyCompletedArtifact(manifest_path, manifest_digest, error)) return false;

  RestartManifest manifest;
  if (!ParseRestartManifest(manifest_path, manifest, error)) return false;
  for (const auto &member : manifest.members) {
    if (member.resolved_path == normalized_artifact) {
      if (!DigestsMatch(member.digest, actual)) {
        error = "restart manifest digest mismatch for member: " + artifact_path;
        return false;
      }
      return true;
    }
  }
  error = "restart manifest does not bind requested artifact: " + artifact_path;
  return false;
}

}  // namespace restart_utils
