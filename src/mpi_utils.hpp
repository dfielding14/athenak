#ifndef MPI_UTILS_HPP_
#define MPI_UTILS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mpi_utils.hpp
//! \brief Small helpers for rendering fatal MPI errors and aborting all ranks.

#include <cstdlib>
#include <iostream>
#include <string>

#include "config.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace mpi_utils {

using FatalCleanupHook = void (*)();

inline FatalCleanupHook &FatalCleanupHookSlot() {
  static FatalCleanupHook hook = nullptr;
  return hook;
}

inline void SetFatalCleanupHook(FatalCleanupHook hook) {
  FatalCleanupHookSlot() = hook;
}

inline void RunFatalCleanupHook() {
  FatalCleanupHook hook = FatalCleanupHookSlot();
  if (hook != nullptr) hook();
}

[[noreturn]] inline void AbortWorld(const std::string& message) {
  RunFatalCleanupHook();
  std::cerr << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

#if MPI_PARALLEL_ENABLED
inline std::string MpiErrorString(int error) {
  char message[MPI_MAX_ERROR_STRING] = {};
  int message_length = 0;
  if (MPI_Error_string(error, message, &message_length) == MPI_SUCCESS &&
      message_length > 0) {
    return std::string(message, message_length);
  }
  return "MPI error code " + std::to_string(error);
}

inline void CheckMpi(int error, const char* context) {
  if (error != MPI_SUCCESS) {
    AbortWorld(std::string(context) + " failed with MPI error: " +
               MpiErrorString(error));
  }
}
#endif

}  // namespace mpi_utils

#endif  // MPI_UTILS_HPP_
