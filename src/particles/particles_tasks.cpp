//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles_tasks.cpp
//! \brief functions that control Particles tasks stored in tasklists in MeshBlockPack

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <vector>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "mesh/mesh_refinement.hpp"
#include "bvals/bvals.hpp"
#include "mhd/mhd.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn TaskStatus Particles::AdaptDeltaF
//! \brief Fit the global adaptive bi-kappa reference state before a particle push.

TaskStatus Particles::AdaptDeltaF(Driver *pdriver, int stage) {
  (void)pdriver;
  (void)stage;
  if (!UsesAdaptiveDeltaF()) return TaskStatus::complete;

  constexpr Real pi = static_cast<Real>(3.141592653589793238462643383279502884L);
  const Real time = pmy_pack->pmesh->time;
  const Real bucket_roundoff =
      static_cast<Real>(64.0)*std::numeric_limits<Real>::epsilon()*
      std::max(static_cast<Real>(1.0), std::abs(time/pic_deltaf_adapt_interval));
  const std::int64_t bucket = static_cast<std::int64_t>(
      std::floor(time/pic_deltaf_adapt_interval + bucket_roundoff));
  if (bucket == pic_deltaf_adapt_last_bucket) return TaskStatus::complete;
  Q017Fence();
  Kokkos::Timer q017_timer;

  auto &pr = prtcl_rdata;
  Real total_weight = 0.0;
  Real perpendicular_moment = 0.0;
  Real parallel_moment = 0.0;
  Kokkos::parallel_reduce(
      "ParticlesAdaptDeltaFFirstMoments",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
      KOKKOS_LAMBDA(const int &p, Real &weight_sum, Real &perpendicular_sum,
                    Real &parallel_sum) {
        const Real weight = pr(IPWT, p);
        if (weight <= 0.0) return;
        weight_sum += weight;
        perpendicular_sum += weight*sqrt(pr(IPVY, p)*pr(IPVY, p) +
                                         pr(IPVZ, p)*pr(IPVZ, p));
        parallel_sum += weight*fabs(pr(IPVX, p));
      },
      Kokkos::Sum<Real>(total_weight),
      Kokkos::Sum<Real>(perpendicular_moment),
      Kokkos::Sum<Real>(parallel_moment));

#if MPI_PARALLEL_ENABLED
  Real local_first_moments[3] = {
      total_weight, perpendicular_moment, parallel_moment};
  Real global_first_moments[3];
  MPI_Allreduce(local_first_moments, global_first_moments, 3, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  total_weight = global_first_moments[0];
  perpendicular_moment = global_first_moments[1];
  parallel_moment = global_first_moments[2];
#endif

  if (!std::isfinite(total_weight) || total_weight <= 0.0 ||
      !std::isfinite(perpendicular_moment) || perpendicular_moment <= 0.0 ||
      !std::isfinite(parallel_moment) || parallel_moment <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Adaptive delta-f global bi-kappa fit requires finite, positive "
              << "particle weight and non-degenerate parallel/perpendicular moments"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const Real xi = (static_cast<Real>(2.0)/pi)*
                  perpendicular_moment/parallel_moment;
  const Real xi2 = xi*xi;
  const Real xi4 = xi2*xi2;
  Real shape_moment = 0.0;
  Kokkos::parallel_reduce(
      "ParticlesAdaptDeltaFShapeMoment",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nprtcl_thispack),
      KOKKOS_LAMBDA(const int &p, Real &shape_sum) {
        const Real weight = pr(IPWT, p);
        if (weight <= 0.0) return;
        shape_sum += weight*sqrt(xi4*pr(IPVX, p)*pr(IPVX, p) +
                                 xi2*(pr(IPVY, p)*pr(IPVY, p) +
                                      pr(IPVZ, p)*pr(IPVZ, p)));
      },
      Kokkos::Sum<Real>(shape_moment));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &shape_moment, 1, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
#endif

  const Real kappa = pic_deltaf_kappa;
  const Real prefactor =
      sqrt(pi*kappa)*(kappa - static_cast<Real>(1.0))*std::tgamma(kappa - 0.5)/
      (static_cast<Real>(2.0)*std::tgamma(kappa + static_cast<Real>(1.0)));
  const Real p0 = prefactor*shape_moment/total_weight;
  if (!std::isfinite(xi) || xi <= 0.0 ||
      !std::isfinite(p0) || p0 <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Adaptive delta-f global bi-kappa fit produced an invalid state"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  pic_deltaf_adaptive_xi = xi;
  pic_deltaf_adaptive_p0 = p0;
  pic_deltaf_adapt_last_bucket = bucket;
  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(17)
              << "PIC adaptive delta-f fit: time=" << time
              << " bucket=" << bucket << " xi=" << xi << " p0=" << p0
              << std::endl;
  }
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::adaptive_deltaf, q017_timer.seconds());
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn  void Particles::AssembleHydroTasks
//! \brief Adds hydro tasks to appropriate task lists used by time integrators.
//! Called by MeshBlockPack::AddPhysics() function directly after Hydro constructor.

void Particles::AssembleTasks(std::map<std::string, std::shared_ptr<TaskList>> tl) {
  TaskID none(0);
  const bool paper_vl2 = UsesPaperVL2Coupling();
  if (paper_vl2 && UsesDeltaF()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "<particles>/pic_physical_mode=paper_mhd_pic staged VL2 coupling "
              << "currently supports full-f only; delta-f staged source semantics "
              << "are not implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  id.adapt_deltaf = tl["before_timeintegrator"]->AddTask(&Particles::AdaptDeltaF,
                                                         this, none);
  id.save_old = tl["before_timeintegrator"]->AddTask(&Particles::SaveOldPositions,
                                                      this, id.adapt_deltaf);
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

  // Paper VL2 applies particle boundary handling at the midpoint and endpoint.
  // Legacy coupled modes migrate once after the full time integrator.
  auto comm_tl = (paper_vl2 ? tl["after_stagen"] :
                  (couple_moments_to_mhd ? tl["after_timeintegrator"] :
                                           tl["before_timeintegrator"]));
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

    TaskID sid = insert_dep;
    if (paper_vl2) {
      sid = stagen_tl->InsertTask(&Particles::Push, this, sid, insert_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Failed to insert staged Particles::Push before "
                  << insert_name << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    sid = stagen_tl->InsertTask(&Particles::SaveOldPositions, this, sid, insert_loc);
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

    if (paper_vl2) {
      sid = stagen_tl->InsertTask(&Particles::DriftPaperCosmicRaysHalfStep,
                                  this, sid, insert_loc);
      if (sid == TaskID(0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to insert staged particle half drift before "
                  << insert_name << std::endl;
        std::exit(EXIT_FAILURE);
      }
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
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->SetNewPrtclGID();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendCnt
//! \brief Wrapper task list function to set share number of partciles communicated with
//! MPI between all ranks

TaskStatus Particles::SendCnt(Driver *pdrive, int stage) {
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->CountSendsAndRecvs();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::InitRecv
//! \brief Wrapper task list function to post non-blocking receives (with MPI).

TaskStatus Particles::InitRecv(Driver *pdrive, int stage) {
  // post receives for particles
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->InitPrtclRecv();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::SendP()
//! \brief Wrapper task list function to pack/send particles

TaskStatus Particles::SendP(Driver *pdrive, int stage) {
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->PackAndSendPrtcls();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::RecvP
//! \brief Wrapper task list function to receive/unpack particles

TaskStatus Particles::RecvP(Driver *pdrive, int stage) {
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->RecvAndUnpackPrtcls();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}


//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearSend
//! \brief Wrapper task list function that checks all MPI sends have completed.

TaskStatus Particles::ClearSend(Driver *pdrive, int stage) {
  // check sends of particles complete
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->ClearPrtclSend();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn TaskList Particles::ClearRecv
//! \brief Wrapper task list function that checks all MPI receives have completed.

TaskStatus Particles::ClearRecv(Driver *pdrive, int stage) {
  // check receives of particles complete
  Q017Fence();
  Kokkos::Timer q017_timer;
  TaskStatus tstat = pbval_part->ClearPrtclRecv();
  ObserveQ017OwnedKokkosViewAllocationBytes();
  Q017Fence();
  AccumulateQ017Timer(Q017ParticleTimer::migration, q017_timer.seconds());
  return tstat;
}

//----------------------------------------------------------------------------------------
//! \fn std::uint64_t Particles::Q017DirectViewAllocationBytes()
//! \brief Return the device-addressable allocation span for direct particle views.

std::uint64_t Particles::Q017DirectViewAllocationBytes() const {
  return static_cast<std::uint64_t>(prtcl_rdata.span())*sizeof(Real) +
         static_cast<std::uint64_t>(prtcl_idata.span())*sizeof(int) +
         static_cast<std::uint64_t>(species_mass.span())*sizeof(Real) +
         static_cast<std::uint64_t>(species_charge.span())*sizeof(Real) +
         static_cast<std::uint64_t>(species_vx0.span())*sizeof(Real) +
         static_cast<std::uint64_t>(species_vy0.span())*sizeof(Real) +
         static_cast<std::uint64_t>(species_vz0.span())*sizeof(Real) +
         static_cast<std::uint64_t>(moments.span())*sizeof(Real) +
         static_cast<std::uint64_t>(coarse_moments.span())*sizeof(Real) +
         static_cast<std::uint64_t>(paper_smooth_mom_records.span())*
             sizeof(PaperSmoothMomentRecord) +
         static_cast<std::uint64_t>(j_edge_x1e.span())*sizeof(Real) +
         static_cast<std::uint64_t>(j_edge_x2e.span())*sizeof(Real) +
         static_cast<std::uint64_t>(j_edge_x3e.span())*sizeof(Real) +
         static_cast<std::uint64_t>(x1_old.span())*sizeof(Real) +
         static_cast<std::uint64_t>(x2_old.span())*sizeof(Real) +
         static_cast<std::uint64_t>(x3_old.span())*sizeof(Real) +
         static_cast<std::uint64_t>(pic_no_mhd_bcc0.span())*sizeof(Real);
}

std::uint64_t Particles::Q017OwnedKokkosViewAllocationBytes() const {
  std::uint64_t bytes = Q017DirectViewAllocationBytes();
  if (pbval_part != nullptr) bytes += pbval_part->Q017OwnedKokkosViewAllocationBytes();
  if (pbval_mom != nullptr) bytes += pbval_mom->Q017OwnedKokkosViewAllocationBytes();
  if (pbval_jedge != nullptr) bytes += pbval_jedge->Q017OwnedKokkosViewAllocationBytes();
  if (pmy_pack->pmesh->pmr != nullptr) {
    bytes += pmy_pack->pmesh->pmr->Q017OwnedKokkosViewAllocationBytes();
  }
  return bytes;
}

void Particles::ObserveQ017OwnedKokkosViewAllocationBytes(std::uint64_t transient_bytes) {
  q017_owned_kokkos_view_high_water_bytes_ =
      std::max(q017_owned_kokkos_view_high_water_bytes_,
               Q017OwnedKokkosViewAllocationBytes() + transient_bytes);
}

std::uint64_t Particles::Q017PaperSmoothHostAllocationBytes() const {
  return (paper_smooth_mom_transport == nullptr) ? 0 :
         paper_smooth_mom_transport->AllocationBytes();
}

void Particles::ObserveQ017PaperSmoothHostAllocationBytes(std::uint64_t transient_bytes) {
  q017_paper_smooth_host_high_water_bytes_ =
      std::max(q017_paper_smooth_host_high_water_bytes_,
               Q017PaperSmoothHostAllocationBytes() + transient_bytes);
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::OutputQ017Telemetry()
//! \brief Emit final-only particle timers and AthenaK-owned allocation telemetry.

void Particles::OutputQ017Telemetry() const {
  constexpr const char* timer_names[nq017_particle_timers] = {
    "adaptive_deltaf", "push", "deposition", "migration"
  };
  const int nranks = global_variable::nranks;
  auto *pm = pmy_pack->pmesh;
  const bool cosmic_ray = (particle_type == ParticleType::cosmic_ray);
  const int nsp = cosmic_ray ? nspecies : 0;
  const int nlevels = std::max(1, pm->max_level - pm->root_level + 1);
  const std::uint64_t record_bytes =
      static_cast<std::uint64_t>(nrdata)*sizeof(Real) +
      static_cast<std::uint64_t>(nidata)*sizeof(int);
  const std::uint64_t resident_bytes =
      static_cast<std::uint64_t>(nprtcl_thispack)*record_bytes;
  const std::uint64_t allocated_bytes = Q017DirectViewAllocationBytes();
  const std::uint64_t owned_allocated_bytes = Q017OwnedKokkosViewAllocationBytes();
  const std::uint64_t owned_high_water_bytes =
      std::max(owned_allocated_bytes, q017_owned_kokkos_view_high_water_bytes_);
  const std::uint64_t paper_smooth_host_allocated_bytes =
      Q017PaperSmoothHostAllocationBytes();
  const std::uint64_t paper_smooth_host_high_water_bytes =
      std::max(paper_smooth_host_allocated_bytes,
               q017_paper_smooth_host_high_water_bytes_);

  std::vector<std::uint64_t> species_counts(nsp, 0);
  std::vector<std::uint64_t> level_counts(nlevels, 0);
  std::vector<std::uint64_t> species_level_counts(nsp*nlevels, 0);
  std::uint64_t invalid_records = 0;
  auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  for (int p=0; p<nprtcl_thispack; ++p) {
    const int gid = h_pi(PGID, p);
    if (gid < 0 || gid >= pm->nmb_total) {
      invalid_records++;
      continue;
    }
    const int logical_level = pm->lloc_eachmb[gid].level;
    const int level_offset = logical_level - pm->root_level;
    if (level_offset < 0 || level_offset >= nlevels) {
      invalid_records++;
      continue;
    }
    level_counts[level_offset]++;
    if (cosmic_ray) {
      const int sp = h_pi(PSP, p);
      if (sp < 0 || sp >= nsp) {
        invalid_records++;
        continue;
      }
      species_counts[sp]++;
      species_level_counts[sp*nlevels + level_offset]++;
    }
  }

  std::array<double, nq017_particle_timers> sum_times = q017_particle_time_;
  std::array<double, nq017_particle_timers> max_times = q017_particle_time_;
  std::array<std::uint64_t, nq017_particle_timers> max_calls = q017_particle_calls_;
  std::uint64_t total_resident_bytes = resident_bytes;
  std::uint64_t total_allocated_bytes = allocated_bytes;
  std::uint64_t max_allocated_bytes = allocated_bytes;
  std::uint64_t total_owned_allocated_bytes = owned_allocated_bytes;
  std::uint64_t max_owned_allocated_bytes = owned_allocated_bytes;
  std::uint64_t sum_owned_high_water_bytes = owned_high_water_bytes;
  std::uint64_t max_owned_high_water_bytes = owned_high_water_bytes;
  std::uint64_t total_paper_smooth_host_allocated_bytes =
      paper_smooth_host_allocated_bytes;
  std::uint64_t max_paper_smooth_host_allocated_bytes =
      paper_smooth_host_allocated_bytes;
  std::uint64_t sum_paper_smooth_host_high_water_bytes =
      paper_smooth_host_high_water_bytes;
  std::uint64_t max_paper_smooth_host_high_water_bytes =
      paper_smooth_host_high_water_bytes;
  std::uint64_t total_invalid_records = invalid_records;
  std::vector<std::uint64_t> global_species_counts = species_counts;
  std::vector<std::uint64_t> global_level_counts = level_counts;
  std::vector<std::uint64_t> global_species_level_counts = species_level_counts;
#if MPI_PARALLEL_ENABLED
  MPI_Reduce(q017_particle_time_.data(), sum_times.data(), nq017_particle_timers,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(q017_particle_time_.data(), max_times.data(), nq017_particle_timers,
             MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(q017_particle_calls_.data(), max_calls.data(), nq017_particle_timers,
             MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&resident_bytes, &total_resident_bytes, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&allocated_bytes, &total_allocated_bytes, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&allocated_bytes, &max_allocated_bytes, 1, MPI_UINT64_T, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&owned_allocated_bytes, &total_owned_allocated_bytes, 1, MPI_UINT64_T,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&owned_allocated_bytes, &max_owned_allocated_bytes, 1, MPI_UINT64_T,
             MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&owned_high_water_bytes, &sum_owned_high_water_bytes, 1, MPI_UINT64_T,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&owned_high_water_bytes, &max_owned_high_water_bytes, 1, MPI_UINT64_T,
             MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&paper_smooth_host_allocated_bytes,
             &total_paper_smooth_host_allocated_bytes, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&paper_smooth_host_allocated_bytes,
             &max_paper_smooth_host_allocated_bytes, 1, MPI_UINT64_T, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&paper_smooth_host_high_water_bytes,
             &sum_paper_smooth_host_high_water_bytes, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&paper_smooth_host_high_water_bytes,
             &max_paper_smooth_host_high_water_bytes, 1, MPI_UINT64_T, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&invalid_records, &total_invalid_records, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  if (nsp > 0) {
    MPI_Reduce(species_counts.data(), global_species_counts.data(), nsp, MPI_UINT64_T,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(species_level_counts.data(), global_species_level_counts.data(),
               nsp*nlevels, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  MPI_Reduce(level_counts.data(), global_level_counts.data(), nlevels, MPI_UINT64_T,
             MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank != 0) return;

  std::cout << std::scientific << std::setprecision(17);
  auto print_scalar = [](const std::string &name, double value) {
    std::cout << "q017.telemetry." << name << "=" << value << std::endl;
  };
  for (int n=0; n<nq017_particle_timers; ++n) {
    const std::string prefix = std::string("timer.particle.") + timer_names[n];
    print_scalar(prefix + ".seconds_rank_max", max_times[n]);
    print_scalar(prefix + ".seconds_rank_mean", sum_times[n]/nranks);
    print_scalar(prefix + ".calls_rank_max", max_calls[n]);
  }
  print_scalar("particle_memory.sync_kernel_timers", pic_q017_sync_kernel_timers);
  print_scalar("particle_memory.record_bytes", record_bytes);
  print_scalar("particle_memory.root_level", pm->root_level);
  print_scalar("particle_memory.max_level", pm->max_level);
  print_scalar("particle_memory.resident_records.bytes_total", total_resident_bytes);
  print_scalar("particle_memory.direct_views.allocated_snapshot_bytes_total",
               total_allocated_bytes);
  print_scalar("particle_memory.direct_views.allocated_snapshot_bytes_rank_max",
               max_allocated_bytes);
  print_scalar("particle_memory.athenak_owned_tracked_kokkos_views."
               "allocated_snapshot_bytes_total",
               total_owned_allocated_bytes);
  print_scalar("particle_memory.athenak_owned_tracked_kokkos_views."
               "allocated_snapshot_bytes_rank_max",
               max_owned_allocated_bytes);
  print_scalar("particle_memory.athenak_owned_tracked_kokkos_views."
               "allocated_high_water_bytes_rank_sum",
               sum_owned_high_water_bytes);
  print_scalar("particle_memory.athenak_owned_tracked_kokkos_views."
               "allocated_high_water_bytes_rank_max",
               max_owned_high_water_bytes);
  print_scalar("particle_memory.paper_smooth_host_transport."
               "allocated_snapshot_bytes_total",
               total_paper_smooth_host_allocated_bytes);
  print_scalar("particle_memory.paper_smooth_host_transport."
               "allocated_snapshot_bytes_rank_max",
               max_paper_smooth_host_allocated_bytes);
  print_scalar("particle_memory.paper_smooth_host_transport."
               "allocated_high_water_bytes_rank_sum",
               sum_paper_smooth_host_high_water_bytes);
  print_scalar("particle_memory.paper_smooth_host_transport."
               "allocated_high_water_bytes_rank_max",
               max_paper_smooth_host_high_water_bytes);
  print_scalar("particle_memory.invalid_records", total_invalid_records);
  for (int sp=0; sp<nsp; ++sp) {
    const std::string prefix = "particle_memory.species." + std::to_string(sp);
    print_scalar(prefix + ".count", global_species_counts[sp]);
    print_scalar(prefix + ".resident_bytes", global_species_counts[sp]*record_bytes);
  }
  for (int level=0; level<nlevels; ++level) {
    const std::string prefix = "particle_memory.level." +
                               std::to_string(pm->root_level + level);
    print_scalar(prefix + ".count", global_level_counts[level]);
    print_scalar(prefix + ".resident_bytes", global_level_counts[level]*record_bytes);
  }
  for (int sp=0; sp<nsp; ++sp) {
    for (int level=0; level<nlevels; ++level) {
      const std::string prefix = "particle_memory.species." + std::to_string(sp) +
                                 ".level." + std::to_string(pm->root_level + level);
      const std::uint64_t count = global_species_level_counts[sp*nlevels + level];
      print_scalar(prefix + ".count", count);
      print_scalar(prefix + ".resident_bytes", count*record_bytes);
    }
  }
}

} // namespace particles
