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
  id.rest_mom = none;
  id.bcs_mom = none;
  id.prol_mom = none;
  id.irecv_jedge = none;
  id.send_jedge = none;
  id.recv_jedge = none;
  id.crecv_jedge = none;
  id.csend_jedge = none;
  id.bcs_jedge = none;
  id.convert_j_edge = none;

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
  id.rest_mom = tl["before_timeintegrator"]->AddTask(&Particles::RestrictMoments, this,
                                                      id.dep_mom);
  id.send_mom = tl["before_timeintegrator"]->AddTask(&Particles::SendMoments, this,
                                                      id.rest_mom);
  id.recv_mom = tl["before_timeintegrator"]->AddTask(&Particles::RecvMoments, this,
                                                      id.send_mom);
  id.crecv_mom = tl["before_timeintegrator"]->AddTask(&Particles::ClearRecvMoments,
                                                       this, id.recv_mom);
  id.csend_mom = tl["before_timeintegrator"]->AddTask(&Particles::ClearSendMoments,
                                                       this, id.crecv_mom);
  id.bcs_mom = tl["before_timeintegrator"]->AddTask(&Particles::ApplyMomentPhysicalBCs,
                                                     this, id.csend_mom);
  id.prol_mom = tl["before_timeintegrator"]->AddTask(&Particles::ProlongateMoments, this,
                                                      id.bcs_mom);

  // WS-H: in coupled mode, move particle migration/communication after field updates.
  auto comm_tl = (couple_moments_to_mhd ? tl["after_timeintegrator"] :
                                          tl["before_timeintegrator"]);
  TaskID comm_dep = (couple_moments_to_mhd ? none : id.prol_mom);
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
    const bool feedback_in_mhd_src =
        use_fluid_feedback &&
        (couple_fluid_feedback_order == CoupledFluidFeedbackOrder::mhd_src_terms);
    TaskID insert_dep = (feedback_in_mhd_src ? pmhd->id.rkupdt : pmhd->id.efld);
    TaskID insert_loc = (feedback_in_mhd_src ? pmhd->id.srctrms : pmhd->id.efldsrc);
    const char *insert_name = (feedback_in_mhd_src ? "MHD::MHDSrcTerms" :
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

    sid = stagen_tl->InsertTask(&Particles::RestrictMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::RestrictMoments before "
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

    sid = stagen_tl->InsertTask(&Particles::ApplyMomentPhysicalBCs, this, sid,
                                insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to insert Particles::ApplyMomentPhysicalBCs before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    sid = stagen_tl->InsertTask(&Particles::ProlongateMoments, this, sid, insert_loc);
    if (sid == TaskID(0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed to insert Particles::ProlongateMoments before "
                << insert_name << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // PR4 direct mode: synchronize deposited edge-currents before EFieldSrc.
    if ((couple_j_to_efield_representation ==
         CoupledCurrentRepresentation::edge_staggered) &&
        (couple_j_deposition_mode ==
         CoupledCurrentDepositionMode::direct_staggered)) {
      TaskID direct_dep = sid | pmhd->id.efld;
      TaskID direct_loc = pmhd->id.efldsrc;
      if ((direct_dep == TaskID(0)) || (direct_loc == TaskID(0))) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to locate MHD::EFieldSrc insertion point for "
                  << "direct edge-current synchronization" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      sid = stagen_tl->InsertTask(&Particles::InitRecvEdgeCurrents, this,
                                  direct_dep, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::InitRecvEdgeCurrents before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.irecv_jedge = sid;

      sid = stagen_tl->InsertTask(&Particles::SendEdgeCurrents, this, sid, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::SendEdgeCurrents before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.send_jedge = sid;

      sid = stagen_tl->InsertTask(&Particles::RecvEdgeCurrents, this, sid, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::RecvEdgeCurrents before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.recv_jedge = sid;

      sid = stagen_tl->InsertTask(&Particles::ClearRecvEdgeCurrents, this,
                                  sid, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::ClearRecvEdgeCurrents before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.crecv_jedge = sid;

      sid = stagen_tl->InsertTask(&Particles::ClearSendEdgeCurrents, this,
                                  sid, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::ClearSendEdgeCurrents before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.csend_jedge = sid;

      sid = stagen_tl->InsertTask(&Particles::ApplyEdgeCurrentPhysicalBCs, this,
                                  sid, direct_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::ApplyEdgeCurrentPhysicalBCs before "
                  << "MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.bcs_jedge = sid;
    }

    // Step I: convert cell-centered deposited J into edge representation
    // before EFieldSrc.
    if ((couple_j_to_efield_representation ==
         CoupledCurrentRepresentation::edge_staggered) &&
        (couple_j_deposition_mode == CoupledCurrentDepositionMode::cc_convert)) {
      // Keep conversion after both CornerE and deposited-moment wrappers.
      TaskID conv_dep = sid | pmhd->id.efld;
      TaskID conv_loc = pmhd->id.efldsrc;
      if ((conv_dep == TaskID(0)) || (conv_loc == TaskID(0))) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to locate MHD::EFieldSrc insertion point for "
                  << "Particles::ConvertCoupledCurrentRepresentation" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      id.convert_j_edge = stagen_tl->InsertTask(
          &Particles::ConvertCoupledCurrentRepresentation, this, conv_dep, conv_loc);
      if (id.convert_j_edge == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert Particles::ConvertCoupledCurrentRepresentation "
                  << "before MHD::EFieldSrc" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::NewGID
//! \brief Wrapper task list function to set new GID for particles that move between
//! MeshBlocks.

TaskStatus Particles::NewGID(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  TaskStatus tstat = pbval_part->SetNewPrtclGID();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendCnt
//! \brief Wrapper task list function to set share number of partciles communicated with
//! MPI between all ranks

TaskStatus Particles::SendCnt(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  TaskStatus tstat = pbval_part->CountSendsAndRecvs();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI).

TaskStatus Particles::InitRecv(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  // post receives for particles
  TaskStatus tstat = pbval_part->InitPrtclRecv();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendP()
//! \brief Wrapper task list function to pack/send particles

TaskStatus Particles::SendP(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  TaskStatus tstat = pbval_part->PackAndSendPrtcls();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::RecvP
//! \brief Wrapper task list function to receive/unpack particles

TaskStatus Particles::RecvP(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  TaskStatus tstat = pbval_part->RecvAndUnpackPrtcls();
  return tstat;
}


//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed.

TaskStatus Particles::ClearSend(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  // check sends of particles complete
  TaskStatus tstat = pbval_part->ClearPrtclSend();
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed.

TaskStatus Particles::ClearRecv(Driver *pdrive, int stage) {
  (void)pdrive;
  (void)stage;
  // check receives of particles complete
  TaskStatus tstat = pbval_part->ClearPrtclRecv();
  return tstat;
}

} // namespace particles
