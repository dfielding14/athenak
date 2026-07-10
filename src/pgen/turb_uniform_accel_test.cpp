//========================================================================================
// AthenaK test-only manufactured uniform-acceleration problem generator.
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_uniform_accel_test.cpp
//! \brief Exercises the production turbulence source update with a uniform force.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/turb_driver.hpp"
#include "tasklist/task_list.hpp"

namespace {

MeshBlockPack *test_pack = nullptr;
Real test_accel_x1 = 0.0;
Real test_accel_x2 = 0.0;
Real test_accel_x3 = 0.0;

// This task is enrolled only by this custom problem-generator build.  It runs after
// the ordinary before-time-integrator tasks, replacing the generated finite-k force
// register without adding a k=0 option to the production turbulence driver.
TaskStatus SetManufacturedUniformForce(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (test_pack == nullptr || test_pack->pturb == nullptr) {
    return TaskStatus::fail;
  }

  auto &indcs = test_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = test_pack->nmb_thispack;
  auto force = test_pack->pturb->force;
  const Real a1 = test_accel_x1;
  const Real a2 = test_accel_x2;
  const Real a3 = test_accel_x3;

  par_for("turb_uniform_accel_test_force", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    force(m, 0, k, j, i) = a1;
    force(m, 1, k, j, i) = a2;
    force(m, 2, k, j, i) = a3;
  });
  return TaskStatus::complete;
}

}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr || pmbp->pturb == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "turb_uniform_accel_test requires <mhd> and <turb_driving> blocks"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "turb_uniform_accel_test requires an ideal MHD equation of state"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  test_pack = pmbp;
  test_accel_x1 = pin->GetReal("problem", "test_accel_x1");
  test_accel_x2 = pin->GetReal("problem", "test_accel_x2");
  test_accel_x3 = pin->GetReal("problem", "test_accel_x3");
  if (!std::isfinite(test_accel_x1) || !std::isfinite(test_accel_x2) ||
      !std::isfinite(test_accel_x3)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Manufactured acceleration components must be finite" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Add the override after all physics modules have assembled their ordinary
  // before-time-integrator tasks.  On a restart the task must be re-enrolled, too.
  auto task_list = pmbp->tl_map.at("before_timeintegrator");
  TaskID dependency = task_list->GetIDLastTask();
  task_list->AddTask(SetManufacturedUniformForce, dependency);
  if (restart) return;

  const Real rho0 = pin->GetReal("problem", "test_rho0");
  const Real pressure0 = pin->GetReal("problem", "test_pressure0");
  const Real vx0 = pin->GetReal("problem", "test_vx0");
  const Real vy0 = pin->GetReal("problem", "test_vy0");
  const Real vz0 = pin->GetReal("problem", "test_vz0");
  if (!std::isfinite(rho0) || rho0 <= 0.0 ||
      !std::isfinite(pressure0) || pressure0 <= 0.0 ||
      !std::isfinite(vx0) || !std::isfinite(vy0) || !std::isfinite(vz0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Manufactured density/pressure must be positive and all initial "
              << "velocities must be finite" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = pmbp->pmhd->u0;
  auto b0 = pmbp->pmhd->b0;
  const Real gamma = pmbp->pmhd->peos->eos_data.gamma;
  const Real energy0 = pressure0/(gamma - 1.0) +
      0.5*rho0*(vx0*vx0 + vy0*vy0 + vz0*vz0);

  par_for("turb_uniform_accel_test_init", DevExeSpace(),
          0, nmb - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m, IDN, k, j, i) = rho0;
    u0(m, IM1, k, j, i) = rho0*vx0;
    u0(m, IM2, k, j, i) = rho0*vy0;
    u0(m, IM3, k, j, i) = rho0*vz0;
    u0(m, IEN, k, j, i) = energy0;

    b0.x1f(m, k, j, i) = 0.0;
    b0.x2f(m, k, j, i) = 0.0;
    b0.x3f(m, k, j, i) = 0.0;
    if (i == ie) b0.x1f(m, k, j, i + 1) = 0.0;
    if (j == je) b0.x2f(m, k, j + 1, i) = 0.0;
    if (k == ke) b0.x3f(m, k + 1, j, i) = 0.0;
  });
}
