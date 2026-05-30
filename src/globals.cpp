//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file globals.cpp
//  \brief namespace containing global variables.
//
// Yes, we all know global variables should NEVER be used, but in fact they are ideal for,
// e.g., global constants that are set once and never changed.  To prevent name collisions
// global variables are wrapped in their own namespace.

#include "athena.hpp"
#include "globals.hpp"

namespace global_variable {
int my_rank;   // MPI rank of this process; set at start of main();
int nranks;    // total number of MPI ranks; set at start of main();
int node_rank = 0;
int node_size = 1;
int node_id = 0;
int nnodes = 1;
bool node_comm_initialized = false;
#if MPI_PARALLEL_ENABLED
MPI_Comm node_comm = MPI_COMM_NULL;
#endif

void InitializeNodeCommunicator() {
  if (node_comm_initialized) return;
#if MPI_PARALLEL_ENABLED
  if (MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, my_rank,
                          MPI_INFO_NULL, &node_comm) != MPI_SUCCESS ||
      MPI_Comm_rank(node_comm, &node_rank) != MPI_SUCCESS ||
      MPI_Comm_size(node_comm, &node_size) != MPI_SUCCESS) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int is_leader = (node_rank == 0) ? 1 : 0;
  int prefix = 0;
  MPI_Exscan(&is_leader, &prefix, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (my_rank == 0) prefix = 0;
  if (node_rank == 0) node_id = prefix;
  MPI_Bcast(&node_id, 1, MPI_INT, 0, node_comm);
  MPI_Allreduce(&is_leader, &nnodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  node_comm_initialized = true;
}

void FinalizeNodeCommunicator() {
#if MPI_PARALLEL_ENABLED
  if (node_comm_initialized && node_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&node_comm);
    node_comm = MPI_COMM_NULL;
  }
#endif
  node_comm_initialized = false;
}

int NodePrefixSum(int local_count) {
#if MPI_PARALLEL_ENABLED
  InitializeNodeCommunicator();
  int prefix = 0;
  MPI_Exscan(&local_count, &prefix, 1, MPI_INT, MPI_SUM, node_comm);
  if (node_rank == 0) prefix = 0;
  return prefix;
#else
  (void)local_count;
  return 0;
#endif
}

int NodeSum(int local_count) {
#if MPI_PARALLEL_ENABLED
  InitializeNodeCommunicator();
  int total = 0;
  MPI_Allreduce(&local_count, &total, 1, MPI_INT, MPI_SUM, node_comm);
  return total;
#else
  return local_count;
#endif
}

int NodeMin(int local_count) {
#if MPI_PARALLEL_ENABLED
  InitializeNodeCommunicator();
  int minimum = 0;
  MPI_Allreduce(&local_count, &minimum, 1, MPI_INT, MPI_MIN, node_comm);
  return minimum;
#else
  return local_count;
#endif
}

int NodeMax(int local_count) {
#if MPI_PARALLEL_ENABLED
  InitializeNodeCommunicator();
  int maximum = 0;
  MPI_Allreduce(&local_count, &maximum, 1, MPI_INT, MPI_MAX, node_comm);
  return maximum;
#else
  return local_count;
#endif
}
}  // namespace global_variable
