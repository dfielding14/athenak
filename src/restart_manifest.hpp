#ifndef RESTART_MANIFEST_HPP_
#define RESTART_MANIFEST_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_manifest.hpp
//! \brief Strict parser and native reader for transactional node restart manifests.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class IOWrapper;

constexpr char kNodeRestartPayloadMarker[] = "AthenaK node restart payload version=1\n";
constexpr std::size_t kNodeRestartPayloadMarkerSize =
    sizeof(kNodeRestartPayloadMarker) - 1;

struct NodeRestartPayload {
  int node;
  int blocks;
  std::uint64_t bytes;
  std::string path;
  std::string full_path;
};

struct NodeRestartSegment {
  int node;
  int gid_start;
  int count;
  int payload_block_start;
};

class NodeRestartManifest {
 public:
  static bool LooksLikeManifest(const std::string &path);
  static bool IsPayloadPath(const std::string &path);
  static NodeRestartManifest Load(const std::string &path);

  const std::string &CanonicalPayloadPath() const { return payloads_[0].full_path; }
  std::uint64_t HeaderSize() const { return header_size_; }
  std::uint64_t DataSize() const { return data_size_; }
  int NumMeshBlocks() const { return nmb_total_; }

  void LoadLocalBlocks(int gid_start, int count, std::uint64_t data_size,
                       std::vector<char> *blocks) const;

 private:
  std::string manifest_path_;
  int nmb_total_;
  std::uint64_t header_size_;
  std::uint64_t data_size_;
  std::vector<NodeRestartPayload> payloads_;
  std::vector<NodeRestartSegment> segments_;
};

[[noreturn]] void FailNodeRestart(const std::string &message);
void CheckNodeRestartPayloadMarker(IOWrapper &input, bool single_file_per_rank,
                                   bool expected);

#endif  // RESTART_MANIFEST_HPP_
