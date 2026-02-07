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
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn  void Particles::AssembleHydroTasks
//! \brief Adds hydro tasks to appropriate task lists used by time integrators.
//! Called by MeshBlockPack::AddPhysics() function directly after Hydro constructor.

void Particles::AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);

  // particle integration done in "before_timeintegrator" task list
  id.save_old = tl["before_timeintegrator"]->AddTask(&Particles::SaveOldPositions,
                                                      this, none);
  id.push = tl["before_timeintegrator"]->AddTask(&Particles::Push, this, id.save_old);
  id.zero_mom = tl["before_timeintegrator"]->AddTask(&Particles::ZeroMoments, this,
                                                      id.push);
  id.irecv_mom = tl["before_timeintegrator"]->AddTask(&Particles::InitRecvMoments,
                                                       this, id.zero_mom);
  id.dep_mom = tl["before_timeintegrator"]->AddTask(&Particles::DepositMoments, this,
                                                     id.irecv_mom);
  id.send_mom = tl["before_timeintegrator"]->AddTask(&Particles::SendMoments, this,
                                                      id.dep_mom);
  id.recv_mom = tl["before_timeintegrator"]->AddTask(&Particles::RecvMoments, this,
                                                      id.send_mom);
  id.crecv_mom = tl["before_timeintegrator"]->AddTask(&Particles::ClearRecvMoments,
                                                       this, id.recv_mom);
  id.csend_mom = tl["before_timeintegrator"]->AddTask(&Particles::ClearSendMoments,
                                                       this, id.crecv_mom);

  // WS-H: in coupled mode, move particle migration/communication after field updates.
  auto comm_tl = (couple_moments_to_mhd ? tl["after_timeintegrator"] :
                                          tl["before_timeintegrator"]);
  TaskID comm_dep = (couple_moments_to_mhd ? none : id.csend_mom);
  id.newgid = comm_tl->AddTask(&Particles::NewGID, this, comm_dep);
  id.count  = comm_tl->AddTask(&Particles::SendCnt, this, id.newgid);
  id.irecv  = comm_tl->AddTask(&Particles::InitRecv, this, id.count);
  id.sendp  = comm_tl->AddTask(&Particles::SendP, this, id.irecv);
  id.recvp  = comm_tl->AddTask(&Particles::RecvP, this, id.sendp);
  id.crecv  = comm_tl->AddTask(&Particles::ClearRecv, this, id.recvp);
  id.csend  = comm_tl->AddTask(&Particles::ClearSend, this, id.crecv);

  // PR2: optionally insert moment deposition wrappers immediately before MHD::EFieldSrc.
  if (couple_moments_to_mhd) {
    auto *pmhd = pmy_pack->pmhd;
    auto stagen_tl = tl["stagen"];
    if ((pmhd == nullptr) || (stagen_tl == nullptr)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true requires a valid MHD "
                << "stagen tasklist insertion point" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    const bool use_fluid_feedback = (couple_moments_momentum_to_mhd ||
                                     couple_moments_energy_to_mhd);
    TaskID insert_dep = (use_fluid_feedback ? pmhd->id.rkupdt : pmhd->id.efld);
    TaskID insert_loc = (use_fluid_feedback ? pmhd->id.srctrms : pmhd->id.efldsrc);
    const char *insert_name = (use_fluid_feedback ? "MHD::MHDSrcTerms" :
                                                    "MHD::EFieldSrc");
    if ((insert_dep == TaskID(0)) || (insert_loc == TaskID(0))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<particles>/couple_moments_to_mhd=true is only supported in the "
                << "single-fluid MHD stagen task path in PR2" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    TaskID sid = stagen_tl->InsertTask(&Particles::SaveOldPositions, this,
                                       insert_dep, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::SaveOldPositions before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::ZeroMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::ZeroMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::InitRecvMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::InitRecvMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::DepositMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::DepositMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::SendMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::SendMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::RecvMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::RecvMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::ClearRecvMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::ClearRecvMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::ClearSendMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::ClearSendMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::NewGID
//! \brief Wrapper task list function to set new GID for particles that move between
//! MeshBlocks.

TaskStatus Particles::NewGID(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_part->SetNewPrtclGID();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendCnt
//! \brief Wrapper task list function to set share number of partciles communicated with
//! MPI between all ranks

TaskStatus Particles::SendCnt(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_part->CountSendsAndRecvs();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI).

TaskStatus Particles::InitRecv(Driver *pdrive, int stage) {
  // post receives for particles
  TaskStatus tstat = pbval_part->InitPrtclRecv();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendP()
//! \brief Wrapper task list function to pack/send particles

TaskStatus Particles::SendP(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_part->PackAndSendPrtcls();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::RecvP
//! \brief Wrapper task list function to receive/unpack particles

TaskStatus Particles::RecvP(Driver *pdrive, int stage) {
  TaskStatus tstat = pbval_part->RecvAndUnpackPrtcls();
  return tstat;
}


//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed.

TaskStatus Particles::ClearSend(Driver *pdrive, int stage) {
  // check sends of particles complete
  TaskStatus tstat = pbval_part->ClearPrtclSend();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed.

TaskStatus Particles::ClearRecv(Driver *pdrive, int stage) {
  // check receives of particles complete
  TaskStatus tstat = pbval_part->ClearPrtclRecv();
  return tstat;
}

} // namespace particles
