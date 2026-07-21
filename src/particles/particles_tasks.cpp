//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_tasks.cpp
//! \brief functions that control Particles tasks stored in tasklists in MeshBlockPack

#include <map>
#include <memory>
#include <string>
#include <iostream>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn  void Particles::AssembleTasks
//! \brief Adds hydro tasks to appropriate task lists used by time integrators.
//! Called by MeshBlockPack::AddPhysics() function directly after Hydro constructor.

void Particles::AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl,
                              TaskID after_integrator_start) {
  TaskID none(0);

  std::shared_ptr<TaskList> list = tl["before_timeintegrator"];
  TaskID push_dependency = none;
  if (pusher == ParticlesPusher::lagrangian_mc) {
    list = tl["after_timeintegrator"];
    // A frame change is instantaneous: lagrangian_mc tracers have positions but no
    // velocity state to boost. Move them afterward with the just-completed step's saved
    // mass fluxes; the boosted fluid fluxes then govern their next step.
    push_dependency = after_integrator_start;
  }

  id.push   = list->AddTask(&Particles::Push, this, push_dependency);
  id.newgid = list->AddTask(&Particles::NewGID, this, id.push);
  id.count  = list->AddTask(&Particles::SendCnt, this, id.newgid);
  id.irecv  = list->AddTask(&Particles::InitRecv, this, id.count);
  id.sendp  = list->AddTask(&Particles::SendP, this, id.irecv);
  id.recvp  = list->AddTask(&Particles::RecvP, this, id.sendp);
  id.crecv  = list->AddTask(&Particles::ClearRecv, this, id.recvp);
  id.csend  = list->AddTask(&Particles::ClearSend, this, id.crecv);
  if (pusher == ParticlesPusher::lagrangian_mc) {
    id.mradj  = list->AddTask(&Particles::AdjustMeshRefinement, this, id.csend);
    id.seed   = list->AddTask(&Particles::SeedDueTracers, this, id.mradj);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::NewGID
//! \brief Wrapper task list function to set new GID for particles that move between
//! MeshBlocks.

TaskStatus Particles::NewGID(Driver *pdrive, int stage) {
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->SetNewPrtclGID();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendCnt
//! \brief Wrapper task list function to set share number of particles communicated with
//! MPI between all ranks

TaskStatus Particles::SendCnt(Driver *pdrive, int stage) {
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->CountSendsAndRecvs();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI).

TaskStatus Particles::InitRecv(Driver *pdrive, int stage) {
  // post receives for particles
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->InitPrtclRecv();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendP()
//! \brief Wrapper task list function to pack/send particles

TaskStatus Particles::SendP(Driver *pdrive, int stage) {
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->PackAndSendPrtcls();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::RecvP
//! \brief Wrapper task list function to receive/unpack particles

TaskStatus Particles::RecvP(Driver *pdrive, int stage) {
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->RecvAndUnpackPrtcls();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}


//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed.

TaskStatus Particles::ClearSend(Driver *pdrive, int stage) {
  // check sends of particles complete
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->ClearPrtclSend();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed.

TaskStatus Particles::ClearRecv(Driver *pdrive, int stage) {
  // check receives of particles complete
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_comm);
  TaskStatus tstat = pbval_part->ClearPrtclRecv();
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_comm);
  return tstat;
}

//----------------------------------------------------------------------------------------

TaskStatus Particles::SeedDueTracers(Driver *pdrive, int stage) {
  pdrive->StartPerformanceRegion(PerformanceRegion::particle_seed);
  if (particle_type == ParticleType::lagrangian_mc) {
    SeedTracersAtTime(pmy_pack->pmesh->time + pmy_pack->pmesh->dt, false);
  }
  pdrive->StopPerformanceRegion(PerformanceRegion::particle_seed);
  return TaskStatus::complete;
}

} // namespace particles
