//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rcp02_output_file_utils_harness.cpp
//! \brief Direct regression harness for shared output-file utility helpers.

#include <sys/stat.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "outputs/output_file_utils.hpp"

namespace {

std::string failure_message;  // NOLINT(runtime/string)

[[noreturn]] void Fail(const std::string &message) {
  std::cerr << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string &message) {
  if (!condition) {
    Fail(message);
  }
}

void CaptureFailure(const std::string &message) {
  Require(failure_message.empty(), "Failure handler was called more than once.");
  failure_message = message;
}

void ResetFailure() {
  failure_message.clear();
}

std::string JoinPath(const std::string &directory, const std::string &leaf) {
  return directory + "/" + leaf;
}

void WriteFile(const std::string &path, const std::string &contents) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  Require(output.good(), "Could not create harness file '" + path + "'.");
  output << contents;
  output.close();
  Require(output.good(), "Could not write harness file '" + path + "'.");
}

std::string ReadFile(const std::string &path) {
  std::ifstream input(path, std::ios::binary);
  Require(input.good(), "Could not read harness file '" + path + "'.");
  std::ostringstream contents;
  contents << input.rdbuf();
  Require(input.good() || input.eof(), "Could not finish reading harness file.");
  return contents.str();
}

bool PathExists(const std::string &path) {
  struct stat status;
  return stat(path.c_str(), &status) == 0;
}

void RunExistingDirectoryAccepted(const std::string &root) {
  output_file_utils::EnsureDirectory(root, 0775, "harness", CaptureFailure);
  Require(failure_message.empty(), "Existing directory was rejected.");
}

void RunConflictingNonDirectoryRejected(const std::string &root) {
  std::string path = JoinPath(root, "conflict");
  WriteFile(path, "not a directory");
  output_file_utils::EnsureDirectory(path, 0775, "harness", CaptureFailure);
  Require(failure_message == "harness output path '" + path +
                                 "' exists but is not a directory.",
          "Conflicting non-directory did not report the expected failure.");
}

void RunSequenceFormatting() {
  Require(output_file_utils::FormatSequence(99999, "harness", CaptureFailure) ==
              "99999",
          "99999 was not formatted correctly.");
  Require(output_file_utils::FormatSequence(100000, "harness", CaptureFailure) ==
              "100000",
          "100000 was not formatted correctly.");
  Require(output_file_utils::FormatSequence(100001, "harness", CaptureFailure) ==
              "100001",
          "100001 was not formatted correctly.");
  Require(failure_message.empty(), "Valid sequence formatting reported a failure.");
}

void RunSequenceRejections() {
  Require(output_file_utils::FormatSequence(-1, "harness", CaptureFailure).empty(),
          "Negative sequence unexpectedly formatted.");
  Require(failure_message ==
              "harness file number -1 is outside the publishable range.",
          "Negative sequence did not report the expected failure.");
  ResetFailure();

  int maximum = std::numeric_limits<int>::max();
  Require(output_file_utils::FormatSequence(maximum, "harness", CaptureFailure).empty(),
          "INT_MAX sequence unexpectedly formatted.");
  Require(failure_message == "harness file number " + std::to_string(maximum) +
                                 " is outside the publishable range.",
          "INT_MAX sequence did not report the expected failure.");
}

void RunCheckedAdvanceBoundary() {
  int maximum = std::numeric_limits<int>::max();
  Require(output_file_utils::AdvanceFileNumber(maximum - 1, "harness",
                                                CaptureFailure) == maximum,
          "INT_MAX - 1 did not advance to INT_MAX.");
  Require(failure_message.empty(), "Advancing INT_MAX - 1 reported a failure.");

  Require(output_file_utils::AdvanceFileNumber(maximum, "harness",
                                                CaptureFailure) == maximum,
          "Rejected INT_MAX advance did not preserve its input.");
  Require(failure_message ==
              "harness file number cannot be advanced beyond int range.",
          "INT_MAX advance did not report the expected failure.");
}

void RunStaleTemporaryReplacement(const std::string &root) {
  std::string final_path = JoinPath(root, "published");
  std::string temporary_path = output_file_utils::TemporaryPath(final_path);
  WriteFile(final_path, "stale final");
  WriteFile(temporary_path, "fresh temporary");

  output_file_utils::PublishTemporaryFile(temporary_path, final_path, "harness",
                                          CaptureFailure);
  Require(failure_message.empty(), "Stale temporary replacement reported a failure.");
  Require(ReadFile(final_path) == "fresh temporary",
          "Published output did not replace the stale final file.");
  Require(!PathExists(temporary_path), "Published temporary file was not removed.");
}

void RunFailedPublicationCleanup(const std::string &root) {
  std::string temporary_path = JoinPath(root, "orphan.tmp");
  std::string final_path = JoinPath(JoinPath(root, "missing"), "published");
  WriteFile(temporary_path, "discard me");

  output_file_utils::PublishTemporaryFile(temporary_path, final_path, "harness",
                                          CaptureFailure);
  Require(failure_message.find("harness could not publish output file '" + final_path +
                               "': ") == 0,
          "Failed publication did not report the expected failure.");
  Require(!PathExists(temporary_path), "Failed publication left its temporary file.");
}

void RunFailedOwnedPathDiscard(const std::string &root) {
  std::string owned_path = JoinPath(root, "nonempty");
  Require(mkdir(owned_path.c_str(), 0775) == 0,
          "Could not create owned-path discard directory.");
  WriteFile(JoinPath(owned_path, "child"), "retain me");
  std::string discard_failure;
  Require(!output_file_utils::TryDiscardOwnedPath(owned_path, &discard_failure),
          "Nonempty owned directory was unexpectedly discarded.");
  Require(discard_failure.find("Could not discard owned output path '" + owned_path +
                               "': ") == 0,
          "Failed owned-path discard did not report the expected failure.");
  Require(PathExists(owned_path), "Failed owned-path discard removed its target.");
}

void RunLexicalTargetNormalization() {
  Require(output_file_utils::LexicallyNormalTarget("bin/../cart/output.{SEQ}.bin") ==
              "cart/output.{SEQ}.bin",
          "Lexical target normalization did not collapse a parent component.");
  Require(output_file_utils::LexicallyNormalTarget("rst/{NODE}/./output.rst") ==
              "rst/{NODE}/output.rst",
          "Lexical target normalization did not collapse a current component.");
}

void RunAdjacentRadiusTokens() {
  double lower = 0.25;
  double upper = std::nextafter(lower, std::numeric_limits<double>::infinity());
  Require(output_file_utils::FormatSphericalSliceRadius(lower) ==
              "r_2.5000000000000000e-01",
          "Canonical spherical-slice radius token changed unexpectedly.");
  Require(output_file_utils::FormatSphericalSliceRadius(lower) !=
              output_file_utils::FormatSphericalSliceRadius(upper),
          "Adjacent representable spherical-slice radii collapsed to one token.");
}

}  // namespace

int main(int argc, char *argv[]) {
  Require(argc == 3, "Usage: rcp02_output_file_utils_harness MODE ROOT");
  std::string mode = argv[1];
  std::string root = argv[2];
  if (mode == "existing_directory") {
    RunExistingDirectoryAccepted(root);
  } else if (mode == "conflicting_non_directory") {
    RunConflictingNonDirectoryRejected(root);
  } else if (mode == "sequence_formatting") {
    RunSequenceFormatting();
  } else if (mode == "sequence_rejections") {
    RunSequenceRejections();
  } else if (mode == "advance_boundary") {
    RunCheckedAdvanceBoundary();
  } else if (mode == "stale_temporary_replacement") {
    RunStaleTemporaryReplacement(root);
  } else if (mode == "failed_publication_cleanup") {
    RunFailedPublicationCleanup(root);
  } else if (mode == "failed_owned_path_discard") {
    RunFailedOwnedPathDiscard(root);
  } else if (mode == "lexical_target_normalization") {
    RunLexicalTargetNormalization();
  } else if (mode == "adjacent_radius_tokens") {
    RunAdjacentRadiusTokens();
  } else {
    Fail("Unknown harness mode.");
  }
  return EXIT_SUCCESS;
}
