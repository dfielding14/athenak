#ifndef OUTPUTS_RESTART_UTILS_HPP_
#define OUTPUTS_RESTART_UTILS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_utils.hpp
//! \brief Atomic-publication and integrity helpers for restart artifacts.

#include <cstdint>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace restart_utils {

inline constexpr std::uint64_t kMeshMetadataMagic = 0x4154484b4d455348ULL;
inline constexpr int kMeshMetadataVersion = 1;

// Fatal failures must terminate the whole MPI world. A rank-local exit can
// strand peers in the next collective.
[[noreturn]] void AbortOnFatalError();

struct FileDigest {
  std::uint64_t size = 0;
  std::uint64_t fnv1a64 = 0;
};

FileDigest ComputeFileDigest(const std::string &path);
void PublishRestartArtifact(const std::string &partial_path,
                            const std::string &final_path);
void WriteCompletionMarker(const std::string &artifact_path,
                           const FileDigest &digest);
void WriteRestartManifest(
    const std::string &manifest_path,
    const std::vector<std::pair<std::string, FileDigest>> &members);
// Manifest paths validate member layout. Payload paths additionally require a
// matching manifest member with the payload's size and checksum.
bool VerifyRestartArtifact(const std::string &artifact_path, std::string &error);
bool VerifyRestartManifestMemberCount(const std::string &manifest_path,
                                      std::size_t expected_members,
                                      std::string &error);

}  // namespace restart_utils

#endif  // OUTPUTS_RESTART_UTILS_HPP_
