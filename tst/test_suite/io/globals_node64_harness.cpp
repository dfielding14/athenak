//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file globals_node64_harness.cpp
//! \brief MPI regression harness for 64-bit node-local sum and prefix helpers.

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>

#include <mpi.h>

#include "globals.hpp"

namespace {

[[noreturn]] void Fail(const char *message) {
  std::cerr << message << std::endl;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  std::exit(EXIT_FAILURE);
}

}  // namespace

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_variable::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_variable::nranks);
  if (global_variable::nranks != 2) {
    Fail("globals_node64_harness requires exactly two MPI ranks.");
  }

  const std::uint64_t local =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()) +
      1ULL + static_cast<std::uint64_t>(global_variable::my_rank);
  const std::uint64_t prefix = global_variable::NodePrefixSum64(local);
  const std::uint64_t total = global_variable::NodeSum64(local);
  const std::uint64_t expected_prefix =
      global_variable::my_rank == 0
          ? 0
          : static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1ULL;
  const std::uint64_t expected_total =
      2ULL*static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 3ULL;
  if (prefix != expected_prefix) Fail("64-bit node prefix was truncated.");
  if (total != expected_total) Fail("64-bit node sum was truncated.");

  global_variable::FinalizeNodeCommunicator();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
