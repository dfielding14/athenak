//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cbin_layout_harness.cpp
//! \brief Direct regression harness for checked coarsened-binary layout arithmetic.

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "outputs/coarsened_binary_layout.hpp"

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

void RunValid() {
  auto layout = coarsened_binary_layout::KernelLayout::Build(
      8, 8, 8, 2, 4, CaptureFailure);
  Require(failure_message.empty(), "Valid layout reported a failure.");
  Require(layout.coarsened_nout1 == 4 && layout.coarsened_nout2 == 4 &&
              layout.coarsened_nout3 == 4,
          "Valid layout retained incorrect coarsened extents.");
  Require(layout.coarsened_cells == 64, "Valid layout retained incorrect cell count.");
  Require(layout.coarsen_factor_cubed == 8,
          "Valid layout retained incorrect factor cube.");
  Require(layout.coarsened_plane == 16,
          "Valid layout retained incorrect coarsened plane.");
  Require(layout.coarsen_factor == 2 && layout.coarsen_factor_squared == 4 &&
              layout.coarsen_factor_cubed_kernel == 8,
          "Valid layout retained incorrect kernel factors.");
  Require(layout.coarsen_iterations == 512,
          "Valid layout retained incorrect coarsening range.");
  Require(layout.normalize_iterations == 256,
          "Valid layout retained incorrect normalization range.");
}

void RunInvalidDivisibility() {
  coarsened_binary_layout::KernelLayout::Build(8, 7, 8, 2, 1, CaptureFailure);
  Require(failure_message ==
              "Coarsened-binary output extents must be divisible by coarsen_factor.",
          "Invalid divisibility did not report the expected failure.");
}

void RunFactorCubeOverflow() {
  constexpr int factor = 3000000;
  coarsened_binary_layout::KernelLayout::Build(
      factor, factor, factor, factor, 1, CaptureFailure);
  Require(failure_message == "coarsened-binary factor cube overflows.",
          "Factor-cube overflow did not report the expected failure.");
}

void RunCoarseningRangeOverflow() {
  int even_maximum = std::numeric_limits<int>::max() - 1;
  coarsened_binary_layout::KernelLayout::Build(
      even_maximum, even_maximum, 4, 2, 1, CaptureFailure);
  Require(failure_message ==
              "coarsened-binary coarsening range exceeds int64 range.",
          "Coarsening-range overflow did not report the expected failure.");
}

void RunNormalizationRangeOverflow() {
  int maximum = std::numeric_limits<int>::max();
  coarsened_binary_layout::KernelLayout::Build(
      maximum, maximum, 1, 1, 4, CaptureFailure);
  Require(failure_message ==
              "coarsened-binary normalization range exceeds int64 range.",
          "Normalization-range overflow did not report the expected failure.");
}

void RunAllocationOverflow() {
  auto layout = coarsened_binary_layout::KernelLayout::Build(
      8, 8, 8, 2, 1, CaptureFailure);
  Require(failure_message.empty(), "Valid allocation layout reported a failure.");
  coarsened_binary_layout::CheckedAllocationElements(
      std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), layout,
      CaptureFailure);
  Require(failure_message == "coarsened-binary allocation elements overflows.",
          "Allocation overflow did not report the expected failure.");
}

void RunAllocationByteOverflow() {
  coarsened_binary_layout::CheckedAllocationBytes(
      std::numeric_limits<std::size_t>::max(), 2, CaptureFailure);
  Require(failure_message == "coarsened-binary allocation bytes overflows.",
          "Allocation-byte overflow did not report the expected failure.");
}

void RunWideKernelIntermediates() {
  auto plane_layout = coarsened_binary_layout::KernelLayout::Build(
      50000, 50000, 1, 1, 1, CaptureFailure);
  Require(failure_message.empty(), "Wide-plane layout reported a failure.");
  Require(plane_layout.coarsened_plane == 2500000000LL,
          "Wide-plane layout truncated its 64-bit plane stride.");

  auto factor_layout = coarsened_binary_layout::KernelLayout::Build(
      65536, 65536, 65536, 65536, 1, CaptureFailure);
  Require(failure_message.empty(), "Wide-factor layout reported a failure.");
  Require(factor_layout.coarsen_factor_squared == 4294967296LL,
          "Wide-factor layout truncated its 64-bit factor square.");
  Require(factor_layout.coarsen_factor_cubed_kernel == 281474976710656LL,
          "Wide-factor layout truncated its 64-bit factor cube.");
}

}  // namespace

int main(int argc, char *argv[]) {
  Require(argc == 2, "Usage: cbin_layout_harness MODE");
  std::string mode = argv[1];
  if (mode == "valid") {
    RunValid();
  } else if (mode == "invalid_divisibility") {
    RunInvalidDivisibility();
  } else if (mode == "factor_cube_overflow") {
    RunFactorCubeOverflow();
  } else if (mode == "coarsening_range_overflow") {
    RunCoarseningRangeOverflow();
  } else if (mode == "normalization_range_overflow") {
    RunNormalizationRangeOverflow();
  } else if (mode == "allocation_overflow") {
    RunAllocationOverflow();
  } else if (mode == "allocation_byte_overflow") {
    RunAllocationByteOverflow();
  } else if (mode == "wide_kernel_intermediates") {
    RunWideKernelIntermediates();
  } else {
    Fail("Unknown harness mode.");
  }
  return EXIT_SUCCESS;
}
