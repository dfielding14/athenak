//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coarsened_binary_layout.hpp
//! \brief Checked arithmetic for coarsened-binary kernel and allocation extents.

#ifndef OUTPUTS_COARSENED_BINARY_LAYOUT_HPP_
#define OUTPUTS_COARSENED_BINARY_LAYOUT_HPP_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace coarsened_binary_layout {

using FailureHandler = void (*)(const std::string &);
using Size = std::size_t;

inline Size CheckedCount(int value, FailureHandler fail, const std::string &context) {
  if (value < 0) {
    fail(context + " is negative.");
  }
  return static_cast<Size>(value);
}

inline Size CheckedMultiply(Size left, Size right, FailureHandler fail,
                            const std::string &context) {
  if (left != 0 && right > std::numeric_limits<Size>::max()/left) {
    fail(context + " overflows.");
  }
  return left*right;
}

inline std::int64_t CheckedKernelRange(Size value, FailureHandler fail,
                                       const std::string &context) {
  if (value > static_cast<Size>(std::numeric_limits<std::int64_t>::max())) {
    fail(context + " exceeds int64 range.");
  }
  return static_cast<std::int64_t>(value);
}

struct KernelLayout {
  int coarsened_nout1;
  int coarsened_nout2;
  int coarsened_nout3;
  Size coarsened_cells;
  Size coarsen_factor_cubed;
  std::int64_t coarsened_plane;
  std::int64_t coarsen_factor;
  std::int64_t coarsen_factor_squared;
  std::int64_t coarsen_factor_cubed_kernel;
  std::int64_t coarsen_iterations;
  std::int64_t normalize_iterations;

  static KernelLayout Build(int nout1, int nout2, int nout3, int coarsen_factor,
                            int number_of_moments, FailureHandler fail) {
    if (nout1 <= 0 || nout2 <= 0 || nout3 <= 0) {
      fail("coarsened-binary emitted extents must be positive.");
    }
    if (coarsen_factor <= 0 || number_of_moments <= 0) {
      fail("coarsened-binary factor and moment count must be positive.");
    }
    if (nout1 % coarsen_factor != 0 || nout2 % coarsen_factor != 0 ||
        nout3 % coarsen_factor != 0) {
      fail("Coarsened-binary output extents must be divisible by coarsen_factor.");
    }
    int cnout1 = nout1/coarsen_factor;
    int cnout2 = nout2/coarsen_factor;
    int cnout3 = nout3/coarsen_factor;
    Size plane = CheckedMultiply(
        CheckedCount(cnout1, fail, "coarsened-binary x1 extent"),
        CheckedCount(cnout2, fail, "coarsened-binary x2 extent"), fail,
        "coarsened-binary coarsened plane");
    Size cells = CheckedMultiply(
        plane,
        CheckedCount(cnout3, fail, "coarsened-binary x3 extent"), fail,
        "coarsened-binary coarsened-cell count");
    Size factor = CheckedCount(coarsen_factor, fail, "coarsened-binary factor");
    Size factor_squared =
        CheckedMultiply(factor, factor, fail, "coarsened-binary factor square");
    Size factor_cubed =
        CheckedMultiply(factor_squared, factor, fail, "coarsened-binary factor cube");
    Size coarsen_range = CheckedMultiply(
        cells, factor_cubed, fail, "coarsened-binary coarsening range");
    Size normalize_range = CheckedMultiply(
        CheckedCount(number_of_moments, fail, "coarsened-binary moment count"),
        cells, fail, "coarsened-binary normalization range");
    return {cnout1, cnout2, cnout3, cells, factor_cubed,
            CheckedKernelRange(plane, fail, "coarsened-binary coarsened plane"),
            CheckedKernelRange(factor, fail, "coarsened-binary factor"),
            CheckedKernelRange(factor_squared, fail, "coarsened-binary factor square"),
            CheckedKernelRange(factor_cubed, fail, "coarsened-binary factor cube"),
            CheckedKernelRange(coarsen_range, fail,
                               "coarsened-binary coarsening range"),
            CheckedKernelRange(normalize_range, fail,
                               "coarsened-binary normalization range")};
  }
};

inline Size CheckedAllocationElements(int variables, int meshblocks,
                                      const KernelLayout &layout,
                                      FailureHandler fail) {
  return CheckedMultiply(
      CheckedMultiply(CheckedCount(variables, fail, "coarsened-binary variable count"),
                      CheckedCount(meshblocks, fail, "coarsened-binary MeshBlock count"),
                      fail, "coarsened-binary allocation elements"),
      layout.coarsened_cells, fail, "coarsened-binary allocation elements");
}

inline Size CheckedAllocationBytes(Size elements, Size element_bytes,
                                   FailureHandler fail) {
  return CheckedMultiply(elements, element_bytes, fail,
                         "coarsened-binary allocation bytes");
}

}  // namespace coarsened_binary_layout

#endif  // OUTPUTS_COARSENED_BINARY_LAYOUT_HPP_
