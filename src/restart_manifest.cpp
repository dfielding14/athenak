//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_manifest.cpp
//! \brief Strict parser and native reader for transactional node restart manifests.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "outputs/io_wrapper.hpp"
#include "restart_manifest.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr const char *kNodeRestartMagic = "AthenaK node restart manifest version=1";
constexpr const char *kNodeRestartPrefix = "AthenaK node restart manifest version=";
constexpr const char *kPayloadSuffix = ".payload.rst";
constexpr int kMaxNodeRestartPayloads = 1024 * 1024;
constexpr std::size_t kMaxNodeRestartSegments = 1024 * 1024;
constexpr std::uintmax_t kMaxNodeRestartManifestBytes = 64ULL*1024ULL*1024ULL;
constexpr std::size_t kMaxManifestSignatureBytes = 256;
constexpr std::size_t kMaxManifestScalarRecordBytes = 256;
constexpr std::size_t kMaxManifestPayloadRecordBytes = 4096;
constexpr std::size_t kMaxManifestSegmentRecordBytes = 256;
constexpr std::size_t kMaxManifestTrailingRecordBytes = 256;
constexpr std::size_t kMaxGeneratedPayloadPathBytes = 1024;

struct NodeRestartSpan {
  int node;
  int local_block_start;
  int payload_block_start;
  int count;
};

enum class BoundedLineResult { line, end_of_file, limit_exceeded, read_error };

bool EndsWith(const std::string &text, const std::string &suffix) {
  return text.size() >= suffix.size() &&
      text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool IsNodeDirectory(const std::string &directory) {
  return directory.size() == 13 && directory.rfind("node_", 0) == 0 &&
      std::all_of(directory.begin() + 5, directory.end(),
                  [](char ch) { return ch >= '0' && ch <= '9'; });
}

bool IsGeneratedPayloadLeaf(const std::string &leaf) {
  std::string suffix(kPayloadSuffix);
  if (!EndsWith(leaf, suffix)) return false;
  std::size_t suffix_begin = leaf.size() - suffix.size();
  std::size_t generation_begin = leaf.rfind(".g", suffix_begin);
  if (generation_begin == std::string::npos || generation_begin == 0 ||
      generation_begin + 2 == suffix_begin) {
    return false;
  }
  return std::all_of(leaf.begin() + generation_begin + 2, leaf.begin() + suffix_begin,
                     [](char ch) { return ch >= '0' && ch <= '9'; });
}

bool IsGeneratedPayloadArtifactLeaf(std::string leaf) {
  constexpr const char *temporary_suffix = ".tmp";
  if (EndsWith(leaf, temporary_suffix)) {
    leaf.resize(leaf.size() - std::string(temporary_suffix).size());
  }
  return IsGeneratedPayloadLeaf(leaf);
}

bool HasGeneratedPayloadContract(const std::filesystem::path &path) {
  return IsNodeDirectory(path.parent_path().filename().string()) &&
      IsGeneratedPayloadArtifactLeaf(path.filename().string());
}

bool IsWithinDirectory(const std::filesystem::path &directory,
                       const std::filesystem::path &path) {
  auto directory_part = directory.begin();
  auto path_part = path.begin();
  while (directory_part != directory.end() && path_part != path.end()) {
    if (*directory_part != *path_part) return false;
    ++directory_part;
    ++path_part;
  }
  return directory_part == directory.end() && path_part != path.end();
}

std::string ParentDirectory(const std::string &path) {
  std::size_t slash = path.rfind('/');
  if (slash == std::string::npos) return ".";
  if (slash == 0) return "/";
  return path.substr(0, slash);
}

std::string LeafName(const std::string &path) {
  std::size_t slash = path.rfind('/');
  return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

std::string JoinPath(const std::string &directory, const std::string &path) {
  return (directory == "/") ? directory + path : directory + "/" + path;
}

std::uint64_t CheckedMultiply(std::uint64_t left, std::uint64_t right,
                              const std::string &context) {
  if (left != 0 && right > std::numeric_limits<std::uint64_t>::max()/left) {
    FailNodeRestart(context + " overflows.");
  }
  return left*right;
}

std::uint64_t CheckedAdd(std::uint64_t left, std::uint64_t right,
                         const std::string &context) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    FailNodeRestart(context + " overflows.");
  }
  return left + right;
}

std::size_t CheckedSize(std::uint64_t value, const std::string &context) {
  if (value > std::numeric_limits<std::size_t>::max()) {
    FailNodeRestart(context + " does not fit in memory.");
  }
  return static_cast<std::size_t>(value);
}

std::uint64_t ParseUnsignedToken(const std::string &token, const std::string &context) {
  if (token.empty() ||
      !std::all_of(token.begin(), token.end(),
                   [](char ch) { return ch >= '0' && ch <= '9'; })) {
    FailNodeRestart("invalid " + context + ".");
  }
  std::size_t parsed = 0;
  std::uint64_t value = 0;
  try {
    value = std::stoull(token, &parsed);
  } catch (const std::exception &) {
    FailNodeRestart("invalid " + context + ".");
  }
  if (parsed != token.size()) {
    FailNodeRestart("invalid " + context + ".");
  }
  return value;
}

int ParseIntToken(const std::string &token, const std::string &context) {
  std::uint64_t value = ParseUnsignedToken(token, context);
  if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    FailNodeRestart(context + " exceeds INT_MAX.");
  }
  return static_cast<int>(value);
}

BoundedLineResult ReadBoundedLine(std::ifstream *input, std::string *line,
                                  std::size_t max_bytes) {
  line->clear();
  while (true) {
    int next = input->get();
    if (next == std::char_traits<char>::eof()) {
      if (input->bad()) return BoundedLineResult::read_error;
      return line->empty() ? BoundedLineResult::end_of_file : BoundedLineResult::line;
    }
    char ch = static_cast<char>(next);
    if (ch == '\n') return BoundedLineResult::line;
    if (line->size() >= max_bytes) return BoundedLineResult::limit_exceeded;
    line->push_back(ch);
  }
}

bool ReadOptionalLine(std::ifstream *input, std::string *line, std::size_t max_bytes,
                      const std::string &context) {
  BoundedLineResult result = ReadBoundedLine(input, line, max_bytes);
  if (result == BoundedLineResult::limit_exceeded) {
    FailNodeRestart(context + " exceeds the " + std::to_string(max_bytes) +
                    "-byte limit.");
  }
  if (result == BoundedLineResult::read_error) {
    FailNodeRestart("could not read " + context + ".");
  }
  return result == BoundedLineResult::line;
}

std::string ReadRequiredLine(std::ifstream *input, const std::string &context,
                             std::size_t max_bytes) {
  std::string line;
  if (!ReadOptionalLine(input, &line, max_bytes, context)) {
    FailNodeRestart("missing " + context + ".");
  }
  return line;
}

std::uint64_t ReadUnsignedRecord(std::ifstream *input, const std::string &key) {
  std::string line = ReadRequiredLine(
      input, "'" + key + "' record", kMaxManifestScalarRecordBytes);
  std::string prefix = key + "=";
  if (line.rfind(prefix, 0) != 0) {
    FailNodeRestart("expected '" + key + "' record, found '" + line + "'.");
  }
  return ParseUnsignedToken(line.substr(prefix.size()), "'" + key + "' value");
}

int ReadIntRecord(std::ifstream *input, const std::string &key) {
  std::uint64_t value = ReadUnsignedRecord(input, key);
  if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    FailNodeRestart("'" + key + "' value exceeds INT_MAX.");
  }
  return static_cast<int>(value);
}

std::vector<std::string> SplitRecord(const std::string &line) {
  std::istringstream input(line);
  std::vector<std::string> fields;
  std::string field;
  while (input >> field) fields.push_back(field);
  return fields;
}

std::string ManifestPayloadPrefix(const std::string &manifest_path) {
  std::string leaf = LeafName(manifest_path);
  constexpr const char *suffix = ".rst";
  if (leaf.size() <= 4 || leaf.compare(leaf.size() - 4, 4, suffix) != 0 ||
      leaf.find('/') != std::string::npos || leaf.find('\\') != std::string::npos) {
    FailNodeRestart("manifest filename does not follow the restart leaf contract.");
  }
  return leaf.substr(0, leaf.size() - 4) + ".g";
}

std::string ValidatePayloadPath(const NodeRestartPayload &payload,
                                const std::string &payload_prefix) {
  if (payload.path.empty() || payload.path[0] == '/' ||
      payload.path.find('\\') != std::string::npos) {
    FailNodeRestart("payload path is not a relative node-shard path.");
  }
  std::size_t slash = payload.path.find('/');
  if (slash == std::string::npos || slash == 0 || slash + 1 >= payload.path.size() ||
      payload.path.find('/', slash + 1) != std::string::npos) {
    FailNodeRestart("payload path must contain exactly one node-directory component.");
  }
  std::string directory = payload.path.substr(0, slash);
  std::string leaf = payload.path.substr(slash + 1);
  char expected_directory[32];
  std::snprintf(expected_directory, sizeof(expected_directory), "node_%08d",
                payload.node);
  if (directory != expected_directory || directory == "." || directory == ".." ||
      leaf == "." || leaf == "..") {
    FailNodeRestart("payload path does not match its declared node directory.");
  }
  std::string suffix(kPayloadSuffix);
  if (leaf.rfind(payload_prefix, 0) != 0 ||
      leaf.size() <= payload_prefix.size() + suffix.size() ||
      !EndsWith(leaf, suffix)) {
    FailNodeRestart("payload leaf does not match the generated restart contract.");
  }
  std::string generation = leaf.substr(
      payload_prefix.size(), leaf.size() - payload_prefix.size() - suffix.size());
  if (generation.empty() ||
      !std::all_of(generation.begin(), generation.end(),
                   [](char ch) { return ch >= '0' && ch <= '9'; })) {
    FailNodeRestart("payload leaf has an invalid generation token.");
  }
  return leaf;
}

std::uint64_t ExpectedPayloadBytes(std::uint64_t header_size, std::uint64_t data_size,
                                   int blocks) {
  return CheckedAdd(header_size,
                    CheckedMultiply(data_size, static_cast<std::uint64_t>(blocks),
                                    "payload byte count"),
                    "payload byte count");
}

void ValidatePayloadHeaders(const std::vector<NodeRestartPayload> &payloads,
                            std::uint64_t header_size) {
  constexpr std::size_t kCompareBytes = 1024*1024;
  std::ifstream canonical(payloads[0].full_path, std::ios::binary);
  if (!canonical.good()) {
    FailNodeRestart("canonical payload header could not be opened.");
  }
  for (std::size_t index = 1; index < payloads.size(); ++index) {
    std::ifstream duplicate(payloads[index].full_path, std::ios::binary);
    if (!duplicate.good()) {
      FailNodeRestart("replicated payload header could not be opened.");
    }
    canonical.clear();
    canonical.seekg(0, std::ios::beg);
    std::uint64_t remaining = header_size;
    std::vector<char> canonical_data(kCompareBytes);
    std::vector<char> duplicate_data(kCompareBytes);
    while (remaining > 0) {
      std::size_t bytes = static_cast<std::size_t>(
          std::min<std::uint64_t>(remaining, kCompareBytes));
      canonical.read(canonical_data.data(), static_cast<std::streamsize>(bytes));
      duplicate.read(duplicate_data.data(), static_cast<std::streamsize>(bytes));
      if (canonical.gcount() != static_cast<std::streamsize>(bytes) ||
          duplicate.gcount() != static_cast<std::streamsize>(bytes) ||
          !std::equal(canonical_data.begin(), canonical_data.begin() + bytes,
                      duplicate_data.begin())) {
        FailNodeRestart("replicated payload header does not match canonical payload 0.");
      }
      remaining -= bytes;
    }
  }
}

std::vector<NodeRestartSpan> RouteLocalSpans(
    const std::vector<NodeRestartSegment> &segments, int gid_start, int count) {
  std::vector<NodeRestartSpan> spans;
  int gid_end = gid_start + count;
  int covered = 0;
  for (const auto &segment : segments) {
    int segment_end = segment.gid_start + segment.count;
    int overlap_start = std::max(gid_start, segment.gid_start);
    int overlap_end = std::min(gid_end, segment_end);
    if (overlap_start >= overlap_end) continue;
    NodeRestartSpan span{
      segment.node,
      overlap_start - gid_start,
      segment.payload_block_start + overlap_start - segment.gid_start,
      overlap_end - overlap_start
    };
    if (!spans.empty() && spans.back().node == span.node &&
        spans.back().local_block_start + spans.back().count == span.local_block_start &&
        spans.back().payload_block_start + spans.back().count ==
            span.payload_block_start) {
      spans.back().count += span.count;
    } else {
      spans.push_back(span);
    }
    covered += span.count;
  }
  if (covered != count) {
    FailNodeRestart("validated segments do not route the local MeshBlock range.");
  }
  return spans;
}

}  // namespace

[[noreturn]] void FailNodeRestart(const std::string &message) {
  std::cerr << "### FATAL ERROR while reading node restart manifest: "
            << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(EXIT_FAILURE);
}

void CheckNodeRestartPayloadMarker(IOWrapper &input, bool single_file_per_rank,
                                   bool expected) {
  IOWrapperSizeT position = input.GetPosition(single_file_per_rank);
  char marker[kNodeRestartPayloadMarkerSize];
  std::size_t read = input.Read_bytes(marker, 1, kNodeRestartPayloadMarkerSize,
                                      single_file_per_rank);
  bool present = read == kNodeRestartPayloadMarkerSize &&
      std::memcmp(marker, kNodeRestartPayloadMarker,
                  kNodeRestartPayloadMarkerSize) == 0;
  if (!present && input.Seek(position, single_file_per_rank) != 0) {
    FailNodeRestart("restart payload marker probe could not restore the file position.");
  }
  if (expected && !present) {
    FailNodeRestart("node restart payload marker is absent or invalid.");
  }
  if (!expected && present) {
    FailNodeRestart("node payload paths are not supported restart entry points; "
                    "use the public manifest path.");
  }
}

bool NodeRestartManifest::LooksLikeManifest(const std::string &path) {
  std::ifstream input(path);
  std::string first_line;
  if (!input.good()) return false;
  BoundedLineResult result =
      ReadBoundedLine(&input, &first_line, kMaxManifestSignatureBytes);
  if (result == BoundedLineResult::limit_exceeded &&
      first_line.rfind(kNodeRestartPrefix, 0) == 0) {
    FailNodeRestart("manifest signature exceeds the " +
                    std::to_string(kMaxManifestSignatureBytes) + "-byte limit.");
  }
  return result == BoundedLineResult::line &&
      first_line.rfind(kNodeRestartPrefix, 0) == 0;
}

bool NodeRestartManifest::IsPayloadPath(const std::string &path) {
  std::filesystem::path normalized = std::filesystem::path(path).lexically_normal();
  if (HasGeneratedPayloadContract(normalized)) return true;
  std::error_code error;
  std::filesystem::path canonical = std::filesystem::canonical(normalized, error);
  return !error && HasGeneratedPayloadContract(canonical);
}

NodeRestartManifest NodeRestartManifest::Load(const std::string &path) {
  NodeRestartManifest manifest;
  manifest.manifest_path_ = path;
  std::error_code error;
  std::uintmax_t manifest_bytes = std::filesystem::file_size(path, error);
  if (error) {
    FailNodeRestart("manifest file size could not be determined.");
  }
  if (manifest_bytes > kMaxNodeRestartManifestBytes) {
    FailNodeRestart("node restart manifest exceeds the " +
                    std::to_string(kMaxNodeRestartManifestBytes) + "-byte limit.");
  }
  std::ifstream input(path);
  std::string line = ReadRequiredLine(
      &input, "manifest signature", kMaxManifestSignatureBytes);
  if (line != kNodeRestartMagic) {
    FailNodeRestart("unsupported manifest signature in '" + path + "'.");
  }
  if (ReadRequiredLine(&input, "'complete' record", kMaxManifestScalarRecordBytes) !=
      "complete=1") {
    FailNodeRestart("invalid completion record.");
  }
  int payload_count = ReadIntRecord(&input, "payload_count");
  if (payload_count <= 0 || payload_count > kMaxNodeRestartPayloads) {
    FailNodeRestart("'payload_count' value must be between 1 and 1048576.");
  }
  manifest.nmb_total_ = ReadIntRecord(&input, "nmb_total");
  if (manifest.nmb_total_ <= 0) {
    FailNodeRestart("'nmb_total' value must be positive.");
  }
  manifest.header_size_ = ReadUnsignedRecord(&input, "header_size");
  if (manifest.header_size_ == 0) {
    FailNodeRestart("'header_size' value must be positive.");
  }
  manifest.data_size_ = ReadUnsignedRecord(&input, "data_size");

  manifest.payloads_.reserve(static_cast<std::size_t>(payload_count));
  for (int expected_node = 0; expected_node < payload_count; ++expected_node) {
    line = ReadRequiredLine(
        &input, "payload inventory record", kMaxManifestPayloadRecordBytes);
    std::vector<std::string> fields = SplitRecord(line);
    if (fields.size() != 5 || fields[0] != "payload") {
      FailNodeRestart("malformed or missing payload inventory record.");
    }
    if (fields[4].size() > kMaxGeneratedPayloadPathBytes) {
      FailNodeRestart("payload inventory record generated path exceeds the " +
                      std::to_string(kMaxGeneratedPayloadPathBytes) + "-byte limit.");
    }
    NodeRestartPayload payload{
      ParseIntToken(fields[1], "payload node"),
      ParseIntToken(fields[2], "payload block count"),
      ParseUnsignedToken(fields[3], "payload byte count"),
      fields[4],
      ""
    };
    if (payload.node != expected_node) {
      FailNodeRestart("payload inventory must be ordered by contiguous node id.");
    }
    if (payload.blocks <= 0) {
      FailNodeRestart("payload block count must be positive.");
    }
    manifest.payloads_.push_back(payload);
  }

  bool saw_end = false;
  while (ReadOptionalLine(
      &input, &line, kMaxManifestSegmentRecordBytes, "segment inventory record")) {
    if (line == "end") {
      saw_end = true;
      break;
    }
    std::vector<std::string> fields = SplitRecord(line);
    if (fields.size() != 5 || fields[0] != "segment") {
      FailNodeRestart("unrecognized or malformed inventory record '" + line + "'.");
    }
    std::size_t segment_limit = std::min(
        kMaxNodeRestartSegments, static_cast<std::size_t>(manifest.nmb_total_));
    if (manifest.segments_.size() >= segment_limit) {
      FailNodeRestart("node restart manifest contains too many segment records.");
    }
    NodeRestartSegment segment{
      ParseIntToken(fields[1], "segment node"),
      ParseIntToken(fields[2], "segment gid start"),
      ParseIntToken(fields[3], "segment block count"),
      ParseIntToken(fields[4], "segment payload block start")
    };
    if (segment.count <= 0) {
      FailNodeRestart("segment block count must be positive.");
    }
    manifest.segments_.push_back(segment);
  }
  if (!saw_end) {
    FailNodeRestart("missing manifest terminator.");
  }
  if (ReadOptionalLine(
      &input, &line, kMaxManifestTrailingRecordBytes, "trailing record")) {
    FailNodeRestart("records found after the manifest terminator.");
  }

  std::string directory = ParentDirectory(path);
  std::filesystem::path canonical_directory =
      std::filesystem::canonical(directory, error);
  if (error) {
    FailNodeRestart("manifest directory could not be canonicalized.");
  }
  std::string payload_prefix = ManifestPayloadPrefix(path);
  std::string generation_leaf;
  for (auto &payload : manifest.payloads_) {
    std::string leaf = ValidatePayloadPath(payload, payload_prefix);
    if (generation_leaf.empty()) {
      generation_leaf = leaf;
    } else if (generation_leaf != leaf) {
      FailNodeRestart("node payload inventory refers to mixed generations.");
    }
    payload.full_path = JoinPath(directory, payload.path);
  }

  std::vector<int> mapped_blocks(manifest.payloads_.size(), 0);
  std::vector<int> next_payload_block(manifest.payloads_.size(), 0);
  int next_gid = 0;
  for (const auto &segment : manifest.segments_) {
    if (segment.node < 0 ||
        segment.node >= static_cast<int>(manifest.payloads_.size())) {
      FailNodeRestart("segment refers to an invalid node.");
    }
    if (segment.gid_start < next_gid) {
      FailNodeRestart("payload segments overlap in the declared MeshBlock range.");
    }
    if (segment.gid_start > next_gid) {
      FailNodeRestart("payload segments leave a gap in the declared MeshBlock range.");
    }
    if (segment.payload_block_start != next_payload_block[segment.node]) {
      FailNodeRestart("segment node-local payload offset is inconsistent.");
    }
    if (segment.count > manifest.nmb_total_ - next_gid ||
        segment.count > manifest.payloads_[segment.node].blocks -
                        segment.payload_block_start) {
      FailNodeRestart("segment node-local payload count is inconsistent.");
    }
    mapped_blocks[segment.node] += segment.count;
    next_payload_block[segment.node] += segment.count;
    next_gid += segment.count;
  }
  if (next_gid != manifest.nmb_total_) {
    FailNodeRestart("payload segments do not cover the declared MeshBlock range.");
  }

  for (std::size_t index = 0; index < manifest.payloads_.size(); ++index) {
    auto &payload = manifest.payloads_[index];
    if (mapped_blocks[index] != payload.blocks ||
        next_payload_block[index] != payload.blocks ||
        payload.bytes != ExpectedPayloadBytes(
            manifest.header_size_, manifest.data_size_, payload.blocks)) {
      FailNodeRestart("payload block count or byte count is inconsistent.");
    }
    std::filesystem::path canonical_payload =
        std::filesystem::canonical(payload.full_path, error);
    if (error || !IsWithinDirectory(canonical_directory, canonical_payload)) {
      FailNodeRestart(error ? "payload '" + payload.path + "' is absent or incomplete."
                            : "payload '" + payload.path +
                                  "' escapes the manifest directory through a symlink.");
    }
    payload.full_path = canonical_payload.string();
    std::ifstream payload_input(payload.full_path, std::ios::binary | std::ios::ate);
    std::streampos position = payload_input.tellg();
    std::uint64_t bytes = payload_input.good() && position >= 0
        ? static_cast<std::uint64_t>(position) : 0;
    if (bytes != payload.bytes) {
      FailNodeRestart("payload '" + payload.path + "' is absent or incomplete.");
    }
  }
  ValidatePayloadHeaders(manifest.payloads_, manifest.header_size_);
  return manifest;
}

void NodeRestartManifest::LoadLocalBlocks(int gid_start, int count,
                                          std::uint64_t data_size,
                                          std::vector<char> *blocks) const {
  if (data_size != data_size_) {
    FailNodeRestart("payload per-block byte count does not match the restart header.");
  }
  if (gid_start < 0 || count < 0 || gid_start > nmb_total_ - count) {
    FailNodeRestart("local MeshBlock range is outside the manifest inventory.");
  }
  std::uint64_t block_bytes = CheckedMultiply(
      data_size, static_cast<std::uint64_t>(count), "local restart buffer byte count");
  blocks->resize(CheckedSize(block_bytes, "local restart buffer byte count"));
  std::vector<NodeRestartSpan> spans = RouteLocalSpans(segments_, gid_start, count);

#if MPI_PARALLEL_ENABLED
  global_variable::InitializeNodeCommunicator();
#endif
  for (std::size_t node = 0; node < payloads_.size(); ++node) {
    std::vector<const NodeRestartSpan *> node_spans;
    for (const auto &span : spans) {
      if (span.node == static_cast<int>(node)) node_spans.push_back(&span);
    }

#if MPI_PARALLEL_ENABLED
    int max_spans = global_variable::NodeMax(static_cast<int>(node_spans.size()));
    if (max_spans == 0) continue;
    IOWrapper payload;
    payload.SetCommunicator(global_variable::node_comm);
    payload.Open(payloads_[node].full_path.c_str(), IOWrapper::FileMode::read);
    for (int index = 0; index < max_spans; ++index) {
      std::uint64_t bytes = 0;
      std::uint64_t source = header_size_;
      char *destination = nullptr;
      if (index < static_cast<int>(node_spans.size())) {
        const NodeRestartSpan &span = *node_spans[index];
        bytes = CheckedMultiply(data_size, static_cast<std::uint64_t>(span.count),
                                "node restart read byte count");
        source = CheckedAdd(
            header_size_,
            CheckedMultiply(data_size,
                            static_cast<std::uint64_t>(span.payload_block_start),
                            "node restart source offset"),
            "node restart source offset");
        std::uint64_t destination_offset = CheckedMultiply(
            data_size, static_cast<std::uint64_t>(span.local_block_start),
            "node restart destination offset");
        destination = (bytes == 0) ? nullptr :
            blocks->data() + CheckedSize(destination_offset,
                                         "node restart destination offset");
      }
      if (payload.Read_bytes_at_all(destination, 1, bytes, source) != bytes) {
        FailNodeRestart("payload data read was incomplete.");
      }
    }
    if (payload.Close() != MPI_SUCCESS) {
      FailNodeRestart("payload file could not be closed.");
    }
#else
    if (node_spans.empty()) continue;
    IOWrapper payload;
    payload.Open(payloads_[node].full_path.c_str(), IOWrapper::FileMode::read);
    for (const NodeRestartSpan *span : node_spans) {
      std::uint64_t bytes = CheckedMultiply(
          data_size, static_cast<std::uint64_t>(span->count),
          "node restart read byte count");
      std::uint64_t source = CheckedAdd(
          header_size_,
          CheckedMultiply(data_size,
                          static_cast<std::uint64_t>(span->payload_block_start),
                          "node restart source offset"),
          "node restart source offset");
      std::uint64_t destination_offset = CheckedMultiply(
          data_size, static_cast<std::uint64_t>(span->local_block_start),
          "node restart destination offset");
      char *destination = (bytes == 0) ? nullptr :
          blocks->data() + CheckedSize(destination_offset,
                                       "node restart destination offset");
      if (payload.Read_bytes_at(destination, 1, bytes, source) != bytes) {
        FailNodeRestart("payload data read was incomplete.");
      }
    }
    if (payload.Close() != 0) {
      FailNodeRestart("payload file could not be closed.");
    }
#endif
  }
}
