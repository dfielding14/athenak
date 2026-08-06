//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file output_file_utils.hpp
//! \brief Shared filesystem and sequence helpers for output writers

#ifndef OUTPUTS_OUTPUT_FILE_UTILS_HPP_
#define OUTPUTS_OUTPUT_FILE_UTILS_HPP_

#include <sys/stat.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

namespace output_file_utils {

using FailureHandler = void (*)(const std::string &);
constexpr std::size_t kDefaultMaxWriterAllocationBytes = 512ULL*1024ULL*1024ULL;

inline std::size_t CheckedSizeAdd(std::size_t left, std::size_t right,
                                  const std::string &context,
                                  FailureHandler fail) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    fail(context + " size overflow.");
  }
  return left + right;
}

inline std::size_t CheckedSizeProduct(std::size_t left, std::size_t right,
                                      const std::string &context,
                                      FailureHandler fail) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max()/left) {
    fail(context + " size overflow.");
  }
  return left*right;
}

inline void RequireAllocationBudget(std::size_t bytes, std::size_t limit,
                                    const std::string &context,
                                    FailureHandler fail) {
  if (bytes > limit) {
    fail(context + " requires " + std::to_string(bytes) +
         " bytes, exceeding max_writer_allocation_bytes=" +
         std::to_string(limit) + ".");
  }
}

inline void EnsureDirectory(const std::string &path, mode_t mode,
                            const std::string &context, FailureHandler fail) {
  if (mkdir(path.c_str(), mode) == 0) {
    return;
  }

  int mkdir_errno = errno;
  if (mkdir_errno == EEXIST) {
    struct stat path_status;
    if (stat(path.c_str(), &path_status) != 0) {
      int stat_errno = errno;
      fail(context + " could not inspect output directory '" + path + "': " +
           std::strerror(stat_errno) + ".");
      return;
    }
    if (S_ISDIR(path_status.st_mode)) {
      return;
    }
    fail(context + " output path '" + path + "' exists but is not a directory.");
    return;
  }

  fail(context + " could not create output directory '" + path + "': " +
       std::strerror(mkdir_errno) + ".");
}

inline std::string FormatSequence(int sequence, const std::string &context,
                                  FailureHandler fail) {
  if (sequence < 0 || sequence >= std::numeric_limits<int>::max()) {
    fail(context + " file number " + std::to_string(sequence) +
         " is outside the publishable range.");
    return "";
  }

  std::string formatted = std::to_string(sequence);
  if (formatted.size() < 5) {
    formatted.insert(0, 5 - formatted.size(), '0');
  }
  return formatted;
}

inline int AdvanceFileNumber(int file_number, const std::string &context,
                             FailureHandler fail) {
  if (file_number >= std::numeric_limits<int>::max()) {
    fail(context + " file number cannot be advanced beyond int range.");
    return file_number;
  }
  return file_number + 1;
}

inline std::string TemporaryPath(const std::string &final_path) {
  return final_path + ".tmp";
}

template <typename Floating>
inline std::string FormatRoundTripScientific(Floating value) {
  std::ostringstream formatted;
  formatted.imbue(std::locale::classic());
  formatted << std::scientific
            << std::setprecision(std::numeric_limits<Floating>::max_digits10 - 1)
            << value;
  return formatted.str();
}

template <typename Floating>
inline std::string FormatSphericalSliceRadius(Floating radius) {
  return "r_" + FormatRoundTripScientific(radius);
}

inline bool TryDiscardOwnedPath(const std::string &path,
                                std::string *failure_message) {
  if (path.empty() || std::remove(path.c_str()) == 0 || errno == ENOENT) {
    return true;
  }
  int remove_errno = errno;
  *failure_message = "Could not discard owned output path '" + path + "': " +
      std::strerror(remove_errno) + ".";
  return false;
}

inline void DiscardOwnedPath(const std::string &path) {
  std::string failure_message;
  if (!TryDiscardOwnedPath(path, &failure_message)) {
    std::cerr << failure_message << std::endl;
  }
}

inline bool TryPublishTemporaryFile(const std::string &temporary_path,
                                    const std::string &final_path,
                                    const std::string &context,
                                    std::string *failure_message) {
  if (temporary_path == final_path) {
    *failure_message = context + " temporary and final output paths must differ.";
    return false;
  }
  if (std::rename(temporary_path.c_str(), final_path.c_str()) != 0) {
    int rename_errno = errno;
    DiscardOwnedPath(temporary_path);
    *failure_message = context + " could not publish output file '" + final_path +
        "': " + std::strerror(rename_errno) + ".";
    return false;
  }
  return true;
}

inline void PublishTemporaryFile(const std::string &temporary_path,
                                 const std::string &final_path,
                                 const std::string &context, FailureHandler fail) {
  std::string failure_message;
  if (!TryPublishTemporaryFile(temporary_path, final_path, context, &failure_message)) {
    fail(failure_message);
  }
}

inline std::string LexicallyNormalTarget(const std::string &target) {
  return std::filesystem::path(target).lexically_normal().generic_string();
}

inline bool IsSafePathComponent(const std::string &component) {
  return !component.empty() && component != "." && component != ".." &&
         component.find('/') == std::string::npos &&
         component.find('\\') == std::string::npos &&
         component.find('\0') == std::string::npos;
}

inline void ValidatePathComponent(const std::string &component,
                                  const std::string &context, FailureHandler fail) {
  if (!IsSafePathComponent(component)) {
    fail(context + " contains an unsafe path component '" + component + "'.");
  }
}

}  // namespace output_file_utils

#endif  // OUTPUTS_OUTPUT_FILE_UTILS_HPP_
