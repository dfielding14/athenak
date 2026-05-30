#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file globals.hpp
//  \brief namespace containing external global variables

#include <cstdint>

#include "config.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace global_variable {
extern int my_rank, nranks;
// The shared-memory communicator is initialized lazily when node-sharded IO is
// requested and remains valid until FinalizeNodeCommunicator().
extern int node_rank, node_size, node_id, nnodes;
extern bool node_comm_initialized;
#if MPI_PARALLEL_ENABLED
extern MPI_Comm node_comm;
#endif

// First initialization is collective over MPI_COMM_WORLD. After initialization,
// the reduction helpers are collective over the shared-memory node communicator.
void InitializeNodeCommunicator();
void FinalizeNodeCommunicator();
int NodePrefixSum(int local_count);
int NodeSum(int local_count);
std::uint64_t NodePrefixSum64(std::uint64_t local_count);
std::uint64_t NodeSum64(std::uint64_t local_count);
int NodeMin(int local_count);
int NodeMax(int local_count);
}  // namespace global_variable

#endif // GLOBALS_HPP_
