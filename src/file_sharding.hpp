#ifndef FILE_SHARDING_HPP_
#define FILE_SHARDING_HPP_
//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file file_sharding.hpp
//! \brief Common naming and mode helpers for files distributed across MPI ranks or nodes.

#include <cstdio>
#include <string>

enum class FileShardMode {
  shared,
  rank,
  node
};

inline bool IsRankSharded(FileShardMode mode) {
  return mode == FileShardMode::rank;
}

inline bool IsNodeSharded(FileShardMode mode) {
  return mode == FileShardMode::node;
}

inline bool IsSharded(FileShardMode mode) {
  return mode != FileShardMode::shared;
}

inline bool UsesIndependentFileIO(FileShardMode mode) {
  return mode == FileShardMode::rank;
}

inline const char *ShardDistributionName(FileShardMode mode) {
  if (mode == FileShardMode::rank) return "rank";
  if (mode == FileShardMode::node) return "node";
  return "shared";
}

inline std::string ShardDirectoryName(FileShardMode mode, int world_rank, int node_id) {
  char name[32];
  if (mode == FileShardMode::rank) {
    std::snprintf(name, sizeof(name), "rank_%08d", world_rank);
    return std::string(name);
  }
  if (mode == FileShardMode::node) {
    std::snprintf(name, sizeof(name), "node_%08d", node_id);
    return std::string(name);
  }
  return std::string();
}

#endif // FILE_SHARDING_HPP_
