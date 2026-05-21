//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//! \brief implementation of functions in class Gravity

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../driver/driver.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "gravity.hpp"
#include "mg_gravity.hpp"
#include "../multigrid/multigrid.hpp"

namespace gravity { // NOLINT (build/namespace)
namespace {
void FatalGravityInput(const std::string &message) {
  std::cout << "### FATAL ERROR in Gravity input validation" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

bool IsPowerOfTwo(int n) {
  return (n > 0) && ((n & (n - 1)) == 0);
}

void ValidateGravityInput(MeshBlockPack *pmbp, ParameterInput *pin) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  if (indcs.nx2 <= 1 || indcs.nx3 <= 1) {
    FatalGravityInput("Self-gravity currently requires a 3D mesh.");
  }
  if (indcs.nx1 != indcs.nx2 || indcs.nx1 != indcs.nx3) {
    FatalGravityInput("Self-gravity requires logically cubic MeshBlocks: "
                      "meshblock/nx1 == meshblock/nx2 == meshblock/nx3.");
  }
  if (!IsPowerOfTwo(indcs.nx1)) {
    FatalGravityInput("Self-gravity requires MeshBlock sizes that are powers of two.");
  }
  if (pin->GetOrAddReal("gravity", "four_pi_G", -1.0) < 0.0) {
    FatalGravityInput("Set gravity/four_pi_G to a non-negative value.");
  }
  const bool has_threshold = pin->DoesParameterExist("gravity", "threshold");
  const bool has_niteration = pin->DoesParameterExist("gravity", "niteration");
  const Real threshold = pin->GetOrAddReal("gravity", "threshold", -1.0);
  const int niteration = pin->GetOrAddInteger("gravity", "niteration", -1);
  if (!has_threshold && !has_niteration) {
    FatalGravityInput("Set either gravity/threshold or gravity/niteration.");
  }
  if (threshold < 0.0 && niteration <= 0) {
    FatalGravityInput("Fixed-iteration self-gravity requires gravity/niteration > 0 "
                      "when gravity/threshold < 0.");
  }
  if (pin->GetOrAddInteger("gravity", "npresmooth", 1) < 0 ||
      pin->GetOrAddInteger("gravity", "npostsmooth", 1) < 0) {
    FatalGravityInput("gravity/npresmooth and gravity/npostsmooth must be non-negative.");
  }
  if (pin->GetOrAddInteger("gravity", "fmg_ncycle", 1) <= 0) {
    FatalGravityInput("gravity/fmg_ncycle must be positive.");
  }
  if (pin->GetOrAddInteger("gravity", "mg_nghost", 1) <= 0) {
    FatalGravityInput("gravity/mg_nghost must be positive.");
  }
  if (pin->GetOrAddReal("gravity", "omega", 1.15) <= 0.0) {
    FatalGravityInput("gravity/omega must be positive.");
  }
  std::string mg_bc = pin->GetOrAddString("gravity", "mg_bc", "none");
  if (mg_bc != "none" && mg_bc != "zerofixed" && mg_bc != "zerograd" &&
      mg_bc != "multipole") {
    FatalGravityInput("Unknown gravity/mg_bc = '" + mg_bc +
                      "'. Valid choices are none, zerofixed, zerograd, multipole.");
  }
  if (mg_bc == "multipole") {
    const int mporder = pin->GetOrAddInteger("gravity", "mporder", 4);
    if (mporder != 2 && mporder != 4) {
      FatalGravityInput("gravity/mporder must be 2 or 4 for multipole boundaries.");
    }
    const bool auto_mporigin = pin->GetOrAddBoolean("gravity", "auto_mporigin", true);
    const bool nodipole = pin->GetOrAddBoolean("gravity", "nodipole", false);
    if (auto_mporigin && nodipole) {
      FatalGravityInput("gravity/auto_mporigin and gravity/nodipole cannot both be true.");
    }
    if (!auto_mporigin &&
        (!pin->DoesParameterExist("gravity", "mporigin_x1") ||
         !pin->DoesParameterExist("gravity", "mporigin_x2") ||
         !pin->DoesParameterExist("gravity", "mporigin_x3"))) {
      FatalGravityInput("Set gravity/mporigin_x1, gravity/mporigin_x2, and "
                        "gravity/mporigin_x3 when gravity/auto_mporigin=false.");
    }
  }
}
} // namespace

//! constructor, initializes data structures and parameters
//-------------------------------------------------------------------------------------
//! \fn Gravity::Gravity(MeshBlockPack *pmbp, ParameterInput *pin)
//! \brief Gravity constructor
Gravity::Gravity(MeshBlockPack *pmbp, ParameterInput *pin):
    pmy_pack(pmbp),
    phi("phi",1,1,1,1,1),
    four_pi_G(-1.0) {

    ValidateGravityInput(pmbp, pin);
    four_pi_G = pin->GetOrAddReal("gravity", "four_pi_G",-1.0);

    // create multigrid driver/solver
    // The driver allocates multigrid instances for root level and meshblock levels
    pmgd = new MGGravityDriver(pmbp, pin);


    // Enroll CellCenteredBoundaryVariable object
    //gbvar.bvar_index = pmb->pbval->bvars.size();
    //pmb->pbval->bvars.push_back(&gbvar);
    //pmb->pbval->pgbvar = &gbvar;
    int nmb = pmy_pack->nmb_thispack;
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(phi, nmb, 1, ncells3, ncells2, ncells1);
}

//----------------------------------------------------------------------------------------
//! \fn Gravity::~Gravity()
//! \brief Gravity destructor
Gravity::~Gravity() {
    delete pmgd;
}

void Gravity::SolveStage(Driver *pdriver, int stage, Real dt) {
  if (pmgd != nullptr) {
    pmgd->Solve(pdriver, stage, dt);
  }
}
} // namespace gravity
