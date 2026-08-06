//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper_serial_harness.cpp
//! \brief Direct serial regression harness for IOWrapper seek and tell failures.

#include <sys/types.h>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "outputs/io_wrapper.hpp"

namespace {

[[noreturn]] void Fail(const std::string &message) {
  std::cerr << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string &message) {
  if (!condition) {
    Fail(message);
  }
}

void PrepareEmptyFile(const char *filename) {
  IOWrapper file;
  file.Open(filename, IOWrapper::FileMode::write);
  Require(file.Close() == 0, "Unable to close temporary harness file.");
}

void RunSeekFailure() {
  IOWrapper input;
  input.Open("/dev/stdin", IOWrapper::FileMode::read);
  Require(input.Seek(0) != 0, "Seek unexpectedly accepted a non-seekable stream.");
  Require(input.Close() == 0, "Unable to close stdin harness stream.");
  std::cout << "serial seek failure observed" << std::endl;
}

void RunMaximumOffset(const char *filename) {
  PrepareEmptyFile(filename);
  IOWrapper input;
  input.Open(filename, IOWrapper::FileMode::read);
  input.Seek(std::numeric_limits<IOWrapperSizeT>::max());
  Fail("Expected serial maximum-offset rejection.");
}

void RunNegativePosition() {
  IOWrapper input;
  input.Open("/dev/stdin", IOWrapper::FileMode::read);
  input.GetPosition();
  Fail("Expected serial position rejection.");
}

void RunPositionedReadSeekFailure() {
  char byte = '\0';
  IOWrapper input;
  input.Open("/dev/stdin", IOWrapper::FileMode::read);
  input.Read_bytes_at(&byte, 1, 1, 0);
  Fail("Expected positioned serial read rejection.");
}

void RunTerminalOffsetOverflow(const char *filename) {
  PrepareEmptyFile(filename);
  char bytes[2] = {'\0', '\0'};
  IOWrapper input;
  input.Open(filename, IOWrapper::FileMode::read);
  input.Read_bytes_at(bytes, 1, 2,
                      static_cast<IOWrapperSizeT>(std::numeric_limits<off_t>::max()));
  Fail("Expected terminal serial offset rejection.");
}

}  // namespace

int main(int argc, char *argv[]) {
  Require(argc == 3, "Usage: io_wrapper_serial_harness MODE FILE");
  std::string mode = argv[1];
  if (mode == "seek_stdin") {
    RunSeekFailure();
  } else if (mode == "maximum_offset") {
    RunMaximumOffset(argv[2]);
  } else if (mode == "negative_position") {
    RunNegativePosition();
  } else if (mode == "positioned_read_seek_failure") {
    RunPositionedReadSeekFailure();
  } else if (mode == "terminal_offset_overflow") {
    RunTerminalOffsetOverflow(argv[2]);
  } else {
    Fail("Unknown harness mode.");
  }
  return EXIT_SUCCESS;
}
