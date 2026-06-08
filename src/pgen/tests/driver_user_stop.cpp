//========================================================================================
// AthenaK test-only problem generator for orderly driver user-stop regression.
//========================================================================================

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

int stop_cycle = 0;
int stop_reason_code = 0;
bool stop_failure = false;

void RequestConfiguredStop(Mesh *pm) {
  if (pm != nullptr && pm->pgen != nullptr && pm->ncycle + 1 == stop_cycle) {
    pm->pgen->RequestUserStop(stop_reason_code, stop_failure);
  }
}

}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  stop_cycle = pin->GetInteger("problem", "test_stop_cycle");
  stop_reason_code = pin->GetInteger("problem", "test_stop_reason_code");
  stop_failure = pin->GetBoolean("problem", "test_stop_failure");
  if (stop_cycle <= 0 || stop_reason_code <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Driver user-stop regression requires positive cycle and reason."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Advection(pin, restart);
  user_work_in_loop = true;
  user_work_in_loop_func = RequestConfiguredStop;
}
