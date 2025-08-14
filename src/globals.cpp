// global variables are wrapped in their own namespace.

#include "athena.hpp"
#include "globals.hpp"

namespace global_variable {
int my_rank;   // MPI rank of this process; set at start of main();
int nranks;    // total number of MPI ranks; set at start of main();
}
