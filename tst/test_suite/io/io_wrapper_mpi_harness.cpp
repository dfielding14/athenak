//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper_mpi_harness.cpp
//! \brief Direct MPI regression harness for IOWrapper defensive paths.

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "outputs/io_wrapper.hpp"

namespace {

[[noreturn]] void Fail(const std::string& message) {
  std::cerr << message << std::endl;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string& message) {
  if (!condition) {
    Fail(message);
  }
}

void CheckMpi(int mpi_error, const char* context) {
  if (mpi_error != MPI_SUCCESS) {
    Fail(std::string(context) + " failed.");
  }
}

void PrepareEmptyFile(const char* filename) {
  IOWrapper file;
  file.Open(filename, IOWrapper::FileMode::write);
  file.Close();
}

void RunCollectiveParticipation(const char* filename, int rank) {
  setenv("ATHENAK_TEST_MAX_MPI_BYTES", "3", 1);
  std::string payload;
  IOWrapperSizeT offset = 0;
  if (rank == 0) {
    payload = "abcdefg";
  } else if (rank == 1) {
    payload = "XY";
    offset = 32;
  }

  IOWrapper output;
  output.Open(filename, IOWrapper::FileMode::write);
  std::size_t written = output.Write_any_type_at_all(
      payload.data(), payload.size(), offset, "byte");
  Require(written == payload.size(), "Collective byte write count mismatch.");
  output.Close();

  std::vector<char> buffer(payload.size());
  char dummy = '\0';
  void* read_buffer = buffer.empty() ? static_cast<void*>(&dummy)
                                     : static_cast<void*>(buffer.data());
  IOWrapper input;
  input.Open(filename, IOWrapper::FileMode::read);
  std::size_t read =
      input.Read_bytes_at_all(read_buffer, 1, payload.size(), offset);
  Require(read == payload.size(), "Collective byte read count mismatch.");
  input.Close();
  Require(std::string(buffer.begin(), buffer.end()) == payload,
          "Collective byte read data mismatch.");

  CheckMpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
  if (rank == 0) {
    std::cout << "collective participation ok" << std::endl;
  }
}

void RunOffsetRangeWrite(const char* filename) {
  char bytes[2] = {};
  IOWrapper output;
  output.Open(filename, IOWrapper::FileMode::write);
  output.Write_any_type_at_all(
      bytes, 2,
      static_cast<IOWrapperSizeT>(std::numeric_limits<MPI_Offset>::max()),
      "byte");
  Fail("Expected positioned write range rejection.");
}

void RunOffsetRangeRead(const char* filename) {
  char bytes[2] = {};
  PrepareEmptyFile(filename);
  IOWrapper input;
  input.Open(filename, IOWrapper::FileMode::read);
  input.Read_bytes_at_all(
      bytes, 1, 2,
      static_cast<IOWrapperSizeT>(std::numeric_limits<MPI_Offset>::max()));
  Fail("Expected positioned read range rejection.");
}

void RunByteCountOverflow(const char* filename) {
  double value = 0.0;
  IOWrapper output;
  output.Open(filename, IOWrapper::FileMode::write);
  output.Write_any_type_at_all(
      &value, std::numeric_limits<IOWrapperSizeT>::max(), 0, "double");
  Fail("Expected byte count overflow rejection.");
}

void RunChunkLimitDisagreement(const char* filename, int rank) {
  setenv("ATHENAK_TEST_MAX_MPI_BYTES", rank == 0 ? "3" : "4", 1);
  char value = 'x';
  IOWrapper output;
  output.Open(filename, IOWrapper::FileMode::write);
  output.Write_any_type_at_all(&value, 1, rank, "byte");
  Fail("Expected chunk limit disagreement rejection.");
}

}  // namespace

int main(int argc, char* argv[]) {
  CheckMpi(MPI_Init(&argc, &argv), "MPI_Init");
  int rank = 0;
  CheckMpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "MPI_Comm_rank");
  Require(argc == 3, "Usage: io_wrapper_mpi_harness MODE FILE");

  std::string mode = argv[1];
  if (mode == "collective") {
    RunCollectiveParticipation(argv[2], rank);
  } else if (mode == "offset_range_write") {
    RunOffsetRangeWrite(argv[2]);
  } else if (mode == "offset_range_read") {
    RunOffsetRangeRead(argv[2]);
  } else if (mode == "byte_count_overflow") {
    RunByteCountOverflow(argv[2]);
  } else if (mode == "chunk_limit_disagreement") {
    RunChunkLimitDisagreement(argv[2], rank);
  } else {
    Fail("Unknown harness mode.");
  }

  CheckMpi(MPI_Finalize(), "MPI_Finalize");
  return EXIT_SUCCESS;
}
