//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_layout_harness.cpp
//! \brief Direct regression harness for checked restart-layout arithmetic.

#include <cstdlib>

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "restart_layout.hpp"

namespace {

using restart_layout::Size;

[[noreturn]] void Fail(const std::string &message) {
  throw std::runtime_error(message);
}

void Require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

template <typename Callable>
void RequireFailure(const std::string &expected, Callable callback) {
  try {
    callback();
  } catch (const std::runtime_error &error) {
    Require(std::string(error.what()).find(expected) != std::string::npos,
            "Unexpected failure message: " + std::string(error.what()));
    return;
  }
  throw std::runtime_error("Expected checked restart-layout failure.");
}

void RunLargeFieldWithoutAllocation() {
  const Size above_int_max =
      static_cast<Size>(std::numeric_limits<int>::max()) + 1;
  restart_layout::PayloadInput input{
    above_int_max, 1, 1, 1, 0, 0, 0, 0, 0, sizeof(double)
  };
  restart_layout::PayloadLayout layout =
      restart_layout::PayloadLayout::Build(input, Fail);
  Size expected = above_int_max*sizeof(double);
  Require(layout.hydro_bytes == expected, "Large hydro field byte count mismatch.");
  Require(layout.data_bytes == expected, "Large payload byte count mismatch.");
  Require(restart_layout::CheckedSizeT(above_int_max, Fail,
                                       "large subview count") == above_int_max,
          "Large subview count narrowed before IO.");
  Require(restart_layout::CheckedOffset(17, layout.data_bytes, 3, Fail,
                                        "restart test offset") ==
              17 + 3*expected,
          "CheckedOffset returned an unexpected byte offset.");
  Require(restart_layout::CheckedPayloadBytes(17, layout.data_bytes, 3, Fail,
                                              "restart test payload") ==
              17 + 3*expected,
          "CheckedPayloadBytes returned an unexpected byte count.");
}

void RunCheckedSubtractUnderflow() {
  RequireFailure("subtract underflows.", []() {
    restart_layout::CheckedSubtract(1, 2, Fail, "subtract");
  });
}

void RunManifestBudgetOverflow() {
  RequireFailure("node restart manifest exceeds the 8-byte limit.", []() {
    restart_layout::ManifestBudget budget{0, 8};
    budget.Add(8, Fail);
    budget.Add(1, Fail);
  });
}

void RunCheckedMultiplyOverflow() {
  RequireFailure("multiply overflows.", []() {
    restart_layout::CheckedMultiply(std::numeric_limits<Size>::max(), 2, Fail,
                                    "multiply");
  });
}

void RunCheckedAddOverflow() {
  RequireFailure("add overflows.", []() {
    restart_layout::CheckedAdd(std::numeric_limits<Size>::max(), 1, Fail, "add");
  });
}

void RunCheckedOffsetOverflow() {
  RequireFailure("offset overflows.", []() {
    restart_layout::CheckedOffset(1, std::numeric_limits<Size>::max(), 2, Fail,
                                  "offset");
  });
}

void RunCheckedPayloadBytesOverflow() {
  RequireFailure("payload overflows.", []() {
    restart_layout::CheckedPayloadBytes(1, std::numeric_limits<Size>::max(), 2, Fail,
                                        "payload");
  });
}

}  // namespace

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: restart_layout_harness MODE" << std::endl;
    return EXIT_FAILURE;
  }
  try {
    std::string mode = argv[1];
    if (mode == "large_field") {
      RunLargeFieldWithoutAllocation();
    } else if (mode == "multiply_overflow") {
      RunCheckedMultiplyOverflow();
    } else if (mode == "add_overflow") {
      RunCheckedAddOverflow();
    } else if (mode == "offset_overflow") {
      RunCheckedOffsetOverflow();
    } else if (mode == "payload_bytes_overflow") {
      RunCheckedPayloadBytesOverflow();
    } else if (mode == "subtract_underflow") {
      RunCheckedSubtractUnderflow();
    } else if (mode == "manifest_budget_overflow") {
      RunManifestBudgetOverflow();
    } else {
      throw std::runtime_error("Unknown harness mode.");
    }
  } catch (const std::runtime_error &error) {
    std::cerr << error.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
