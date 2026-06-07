//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_part.cpp
//! \brief

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/nghbr_index.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "bvals.hpp"

namespace {

[[noreturn]] void FatalParticleCountOverflow(const std::string &context) {
  std::cout << "### FATAL ERROR in bvals_part.cpp" << std::endl
            << "particle MPI count overflow while " << context << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

int CheckedMpiIntAdd(int lhs, int rhs, const std::string &context) {
  if (lhs < 0 || rhs < 0 || lhs > std::numeric_limits<int>::max() - rhs) {
    FatalParticleCountOverflow(context);
  }
  return lhs + rhs;
}

int CheckedMpiIntProduct(int lhs, int rhs, const std::string &context) {
  if (lhs < 0 || rhs < 0 ||
      (rhs != 0 && lhs > std::numeric_limits<int>::max()/rhs)) {
    FatalParticleCountOverflow(context);
  }
  return lhs*rhs;
}

int CheckedSizeToMpiInt(std::size_t value, const std::string &context) {
  if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    FatalParticleCountOverflow(context);
  }
  return static_cast<int>(value);
}

} // namespace

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::UpdateGID()
//! \brief Updates GID of particles that cross boundary of their parent MeshBlock.  If
//! the new GID is on a different rank, then store in sendlist_buf DvceArray: (1) index of
//! particle in prtcl array, (2) destination GID, and (3) destination rank.

KOKKOS_INLINE_FUNCTION
void UpdateGID(int &newgid, NeighborBlock nghbr, int myrank, int *pcounter,
               DualArray1D<ParticleLocationData> slist, int p) {
  newgid = nghbr.gid;
#if MPI_PARALLEL_ENABLED
  if (nghbr.rank != myrank) {
    int index = Kokkos::atomic_fetch_add(pcounter,1);
    slist.d_view(index).prtcl_indx = p;
    slist.d_view(index).dest_gid   = nghbr.gid;
    slist.d_view(index).dest_rank  = nghbr.rank;
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::SetNewGID()
//! \brief

TaskStatus ParticlesBoundaryValues::SetNewPrtclGID() {
  // create local references for variables in kernel
  auto gids = pmy_part->pmy_pack->gids;
  auto &pr = pmy_part->prtcl_rdata;
  auto &pi = pmy_part->prtcl_idata;
  int npart = pmy_part->nprtcl_thispack;
  auto &mbsize = pmy_part->pmy_pack->pmb->mb_size;
  auto &mblev = pmy_part->pmy_pack->pmb->mb_lev;
  auto &meshsize = pmy_part->pmy_pack->pmesh->mesh_size;
  auto myrank = global_variable::my_rank;
  auto &nghbr = pmy_part->pmy_pack->pmb->nghbr;
  auto &psendl = sendlist;
  auto &scount = send_count;
  Kokkos::deep_copy(scount, 0);
  int *pcounter = scount.data();
  bool &multi_d = pmy_part->pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_part->pmy_pack->pmesh->three_d;

  Kokkos::realloc(sendlist, npart);
  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    int mylevel = mblev.d_view(m);
    Real x1 = pr(IPX,p);
    Real x2 = pr(IPY,p);
    Real x3 = pr(IPZ,p);

    // length of MeshBlock in each direction
    Real lx = (mbsize.d_view(m).x1max - mbsize.d_view(m).x1min);
    Real ly = (mbsize.d_view(m).x2max - mbsize.d_view(m).x2min);
    Real lz = (mbsize.d_view(m).x3max - mbsize.d_view(m).x3min);

    // integer offset of particle relative to center of MeshBlock (-1,0,+1)
    int ix = static_cast<int>((x1 - mbsize.d_view(m).x1min + lx)/lx) - 1;
    int iy = static_cast<int>((x2 - mbsize.d_view(m).x2min + ly)/ly) - 1;
    int iz = static_cast<int>((x3 - mbsize.d_view(m).x3min + lz)/lz) - 1;

    // sublock indices for faces and edges with S/AMR
    int fx = (x1 < 0.5*(mbsize.d_view(m).x1min + mbsize.d_view(m).x1max))? 0 : 1;
    int fy = (x2 < 0.5*(mbsize.d_view(m).x2min + mbsize.d_view(m).x2max))? 0 : 1;
    int fz = (x3 < 0.5*(mbsize.d_view(m).x3min + mbsize.d_view(m).x3max))? 0 : 1;
    fy = multi_d ? fy : 0;
    fz = three_d ? fz : 0;

    // only update particle GID if it has crossed MeshBlock boundary
    if ((abs(ix) + abs(iy) + abs(iz)) != 0) {
      if (iz == 0) {
        if (iy == 0) {
          // x1 face
          int indx = NeighborIndex(ix,0,0,0,0);           // neighbor at same level
          if (nghbr.d_view(m,indx).lev > mylevel) {       // neighbor at finer level
            indx = NeighborIndex(ix,0,0,fy,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}  // neighbor at coarser level
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        } else if (ix == 0) {
          // x2 face
          int indx = NeighborIndex(0,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,0,fx,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        } else {
          // x1x2 edge
          int indx = NeighborIndex(ix,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,iy,0,fz,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        }
      } else if (iy == 0) {
        if (ix == 0) {
          // x3 face
          int indx = NeighborIndex(0,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,0,iz,fx,fy);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        } else {
          // x3x1 edge
          int indx = NeighborIndex(ix,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,0,iz,fy,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        }
      } else {
        if (ix == 0) {
          // x2x3 edge
          int indx = NeighborIndex(0,iy,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,iz,fx,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        } else {
          // corners
          int indx = NeighborIndex(ix,iy,iz,0,0);
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, pcounter, psendl, p);
        }
      }

      // reset x,y,z positions if particle crosses Mesh boundary using periodic BCs
      if (x1 < meshsize.x1min) {
        pr(IPX,p) += (meshsize.x1max - meshsize.x1min);
      } else if (x1 >= meshsize.x1max) {
        pr(IPX,p) -= (meshsize.x1max - meshsize.x1min);
      }
      if (x2 < meshsize.x2min) {
        pr(IPY,p) += (meshsize.x2max - meshsize.x2min);
      } else if (x2 >= meshsize.x2max) {
        pr(IPY,p) -= (meshsize.x2max - meshsize.x2min);
      }
      if (x3 < meshsize.x3min) {
        pr(IPZ,p) += (meshsize.x3max - meshsize.x3min);
      } else if (x3 >= meshsize.x3max) {
        pr(IPZ,p) -= (meshsize.x3max - meshsize.x3min);
      }
    }
  });
  HostArray1D<int> h_send_count("particle_send_count_host", 1);
  Kokkos::deep_copy(h_send_count, scount);
  nprtcl_send = h_send_count(0);
  Kokkos::resize(sendlist, nprtcl_send);
  // sync sendlist device array with host
  sendlist.template modify<DevExeSpace>();
  sendlist.template sync<HostMemSpace>();

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::CountSendsAndRecvs()
//! \brief

TaskStatus ParticlesBoundaryValues::CountSendsAndRecvs() {
#if MPI_PARALLEL_ENABLED
  // Sort sendlist on host by destrank.
  namespace KE = Kokkos::Experimental;
  std::sort(KE::begin(sendlist.h_view), KE::end(sendlist.h_view), SortByRank);
  // sync sendlist host array with device.  This results in sorted array on device
  sendlist.template modify<HostMemSpace>();
  sendlist.template sync<DevExeSpace>();

  // load STL::vector of ParticleMessageData with <sendrank, recvrank, nprtcls> for sends
  // from this rank. Length will be nsends; initially this length is unknown
  sends_thisrank.clear();
  if (nprtcl_send > 0) {
    int &myrank = global_variable::my_rank;
    int rank = sendlist.h_view(0).dest_rank;
    int nprtcl = 1;

    for (int n=1; n<nprtcl_send; ++n) {
      if (sendlist.h_view(n).dest_rank == rank) {
        ++nprtcl;
      } else {
        sends_thisrank.emplace_back(ParticleMessageData(myrank,rank,nprtcl));
        rank = sendlist.h_view(n).dest_rank;
        nprtcl = 1;
      }
    }
    sends_thisrank.emplace_back(ParticleMessageData(myrank,rank,nprtcl));
  }
  nsends = CheckedSizeToMpiInt(sends_thisrank.size(), "counting destination ranks");

  // Share number of ranks to send to amongst all ranks
  nsends_eachrank[global_variable::my_rank] = nsends;
  MPI_Allgather(&nsends, 1, MPI_INT, nsends_eachrank.data(), 1, MPI_INT, mpi_comm_part);

  // Now share ParticleMessageData amongst all ranks
  // First create vector of starting indices in full vector
  std::vector<int> nsends_displ;
  nsends_displ.resize(global_variable::nranks);
  nsends_displ[0] = 0;
  for (int n=1; n<(global_variable::nranks); ++n) {
    nsends_displ[n] = CheckedMpiIntAdd(
        nsends_displ[n-1], nsends_eachrank[n-1],
        "building particle-message Allgatherv displacements");
  }
  int nsends_allranks = CheckedMpiIntAdd(
      nsends_displ[global_variable::nranks - 1],
      nsends_eachrank[global_variable::nranks - 1],
      "counting particle messages across ranks");
  // Load ParticleMessageData on this rank into full vector
  sends_allranks.resize(nsends_allranks, ParticleMessageData(0,0,0));
  for (int n=0; n<nsends_eachrank[global_variable::my_rank]; ++n) {
    sends_allranks[n + nsends_displ[global_variable::my_rank]] = sends_thisrank[n];
  }

  // Share tuples using MPI derived data type for tuple of 3*int
  MPI_Datatype mpi_ituple;
  MPI_Type_contiguous(3, MPI_INT, &mpi_ituple);
  MPI_Type_commit(&mpi_ituple);
  MPI_Allgatherv(MPI_IN_PLACE, nsends_eachrank[global_variable::my_rank],
                   mpi_ituple, sends_allranks.data(), nsends_eachrank.data(),
                   nsends_displ.data(), mpi_ituple, mpi_comm_part);
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::InitPrtclRecv()
//! \brief

TaskStatus ParticlesBoundaryValues::InitPrtclRecv() {
#if MPI_PARALLEL_ENABLED
  // load STL::vector of ParticleMessageData with <sendrank,recvrank,nprtcl_recv> for
  // receives // on this rank. Length will be nrecvs, initially this length is unknown
  recvs_thisrank.clear();

  int nsends_allranks = CheckedSizeToMpiInt(
      sends_allranks.size(), "reading gathered particle messages");
  for (int n=0; n<nsends_allranks; ++n) {
    if (sends_allranks[n].recvrank == global_variable::my_rank) {
      recvs_thisrank.emplace_back(sends_allranks[n]);
    }
  }
  nrecvs = CheckedSizeToMpiInt(recvs_thisrank.size(), "counting source ranks");

  // Figure out how many particles will be received from all ranks
  nprtcl_recv=0;
  for (int n=0; n<nrecvs; ++n) {
    nprtcl_recv = CheckedMpiIntAdd(
        nprtcl_recv, recvs_thisrank[n].nprtcls,
        "counting particles received by a rank");
  }

  // Allocate receive buffer
  int recv_real_total = CheckedMpiIntProduct(
      pmy_part->nrdata, nprtcl_recv, "allocating received particle real data");
  int recv_int_total = CheckedMpiIntProduct(
      pmy_part->nidata, nprtcl_recv, "allocating received particle integer data");
  Kokkos::realloc(prtcl_rrecvbuf, recv_real_total);
  Kokkos::realloc(prtcl_irecvbuf, recv_int_total);
  Kokkos::realloc(prtcl_trecvbuf, nprtcl_recv);

  // Post non-blocking receives
  bool no_errors=true;
  rrecv_req.clear();
  irecv_req.clear();
  trecv_req.clear();
  for (int n=0; n<nrecvs; ++n) {
    rrecv_req.emplace_back(MPI_REQUEST_NULL);
    irecv_req.emplace_back(MPI_REQUEST_NULL);
    trecv_req.emplace_back(MPI_REQUEST_NULL);
  }

  // Init receives for Reals
  int data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    // calculate amount of data to be passed, get pointer to variables
    int data_size = CheckedMpiIntProduct(
        pmy_part->nrdata, recvs_thisrank[n].nprtcls,
        "forming a particle real-data receive count");
    int data_end = CheckedMpiIntAdd(
        data_start, data_size, "forming a particle real-data receive displacement");
    auto recv_ptr = Kokkos::subview(prtcl_rrecvbuf,
                                    std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 0; // 0 for Reals, 1 for ints

    // Post non-blocking receive
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                         mpi_comm_part, &(rrecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start = data_end;
  }
  // Init receives for 64-bit particle tags.
  data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    int data_size = recvs_thisrank[n].nprtcls;
    int data_end = CheckedMpiIntAdd(
        data_start, data_size, "forming a particle tag receive displacement");
    auto recv_ptr = Kokkos::subview(prtcl_trecvbuf, std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 2;
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_UINT64_T, drank, tag,
                         mpi_comm_part, &(trecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start = data_end;
  }
  // Init receives for ints
  data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    // calculate amount of data to be passed, get pointer to variables
    int data_size = CheckedMpiIntProduct(
        pmy_part->nidata, recvs_thisrank[n].nprtcls,
        "forming a particle integer-data receive count");
    int data_end = CheckedMpiIntAdd(
        data_start, data_size, "forming a particle integer-data receive displacement");
    auto recv_ptr = Kokkos::subview(prtcl_irecvbuf,
                                    std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 1; // 0 for Reals, 1 for ints

    // Post non-blocking receive
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_INT, drank, tag,
                         mpi_comm_part, &(irecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start = data_end;
  }

  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::PackAndSendPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::PackAndSendPrtcls() {
#if MPI_PARALLEL_ENABLED
  // Figure out how many particles will be sent from this ranks
  nprtcl_send=0;
  for (int n=0; n<nsends; ++n) {
    nprtcl_send = CheckedMpiIntAdd(
        nprtcl_send, sends_thisrank[n].nprtcls,
        "counting particles sent by a rank");
  }

  bool no_errors=true;
  if (nprtcl_send > 0) {
    // Allocate send buffer
    int send_real_total = CheckedMpiIntProduct(
        pmy_part->nrdata, nprtcl_send, "allocating sent particle real data");
    int send_int_total = CheckedMpiIntProduct(
        pmy_part->nidata, nprtcl_send, "allocating sent particle integer data");
    Kokkos::realloc(prtcl_rsendbuf, send_real_total);
    Kokkos::realloc(prtcl_isendbuf, send_int_total);
    Kokkos::realloc(prtcl_tsendbuf, nprtcl_send);

    // sendlist on device is already sorted by destrank in CountSendAndRecvs()
    // Use sendlist on device to load particles into send buffer ordered by dest_rank
    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto &pr = pmy_part->prtcl_rdata;
    auto &pi = pmy_part->prtcl_idata;
    auto &ptag = pmy_part->prtcl_tag;
    auto &rsendbuf = prtcl_rsendbuf;
    auto &isendbuf = prtcl_isendbuf;
    auto &tsendbuf = prtcl_tsendbuf;
    par_for("ppack",DevExeSpace(),0,(nprtcl_send-1), KOKKOS_LAMBDA(const int n) {
      int p = sendlist.d_view(n).prtcl_indx;
      for (int i=0; i<nidata; ++i) {
        isendbuf(nidata*n + i) = pi(i,p);
      }
      for (int i=0; i<nrdata; ++i) {
        rsendbuf(nrdata*n + i) = pr(i,p);
      }
      tsendbuf(n) = ptag(p);
    });

    // Post non-blocking sends
    Kokkos::fence();
    rsend_req.clear();
    isend_req.clear();
    tsend_req.clear();
    for (int n=0; n<nsends; ++n) {
      rsend_req.emplace_back(MPI_REQUEST_NULL);
      isend_req.emplace_back(MPI_REQUEST_NULL);
      tsend_req.emplace_back(MPI_REQUEST_NULL);
    }

    // Send Reals
    int data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = CheckedMpiIntProduct(
          nrdata, sends_thisrank[n].nprtcls,
          "forming a particle real-data send count");
      int data_end = CheckedMpiIntAdd(
          data_start, data_size, "forming a particle real-data send displacement");
      auto send_ptr = Kokkos::subview(prtcl_rsendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 0; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                           mpi_comm_part, &(rsend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start = data_end;
    }
    // Send 64-bit tags.
    data_start=0;
    for (int n=0; n<nsends; ++n) {
      int data_size = sends_thisrank[n].nprtcls;
      int data_end = CheckedMpiIntAdd(
          data_start, data_size, "forming a particle tag send displacement");
      auto send_ptr = Kokkos::subview(prtcl_tsendbuf,
                                      std::make_pair(data_start, data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 2;
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_UINT64_T, drank, tag,
                           mpi_comm_part, &(tsend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start = data_end;
    }
    // Send ints
    data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = CheckedMpiIntProduct(
          nidata, sends_thisrank[n].nprtcls,
          "forming a particle integer-data send count");
      int data_end = CheckedMpiIntAdd(
          data_start, data_size, "forming a particle integer-data send displacement");
      auto send_ptr = Kokkos::subview(prtcl_isendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 1; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_INT, drank, tag,
                           mpi_comm_part, &(isend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start = data_end;
    }
  }

  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::RecvAndUnpackPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::RecvAndUnpackPrtcls() {
#if MPI_PARALLEL_ENABLED
  // Sort sendlist on host by index in particle array
  namespace KE = Kokkos::Experimental;
  std::sort(KE::begin(sendlist.h_view), KE::end(sendlist.h_view), SortByIndex);
  // sync sendlist host array with device.  This results in sorted array on device
  sendlist.template modify<HostMemSpace>();
  sendlist.template sync<DevExeSpace>();

  // increase size of particle arrays if needed
  int particles_before_removal = CheckedMpiIntAdd(
      pmy_part->nprtcl_thispack, nprtcl_recv,
      "resizing particle arrays after boundary exchange");
  if (nprtcl_send > particles_before_removal) {
    FatalParticleCountOverflow("resizing particle arrays after boundary exchange");
  }
  int new_npart = particles_before_removal - nprtcl_send;
  if (nprtcl_recv > nprtcl_send) {
    Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
    Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
    Kokkos::resize(pmy_part->prtcl_tag, new_npart);
  }

  // check that particle communications have all completed
  bool bflag = false;
  bool no_errors=true;
  for (int n=0; n<nrecvs; ++n) {
    int test;
    int ierr = MPI_Test(&(rrecv_req[n]), &test, MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (!(static_cast<bool>(test))) {
      bflag = true;
    }
    ierr = MPI_Test(&(irecv_req[n]), &test, MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (!(static_cast<bool>(test))) {
      bflag = true;
    }
    ierr = MPI_Test(&(trecv_req[n]), &test, MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (!(static_cast<bool>(test))) {
      bflag = true;
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // exit if particle communications have not completed
  if (bflag) {return TaskStatus::incomplete;}

  // unpack particles into positions of sent particles
  if (nprtcl_recv > 0) {
    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto &pr = pmy_part->prtcl_rdata;
    auto &pi = pmy_part->prtcl_idata;
    auto &ptag = pmy_part->prtcl_tag;
    auto &rrecvbuf = prtcl_rrecvbuf;
    auto &irecvbuf = prtcl_irecvbuf;
    auto &trecvbuf = prtcl_trecvbuf;
    int &npart = pmy_part->nprtcl_thispack;
    par_for("punpack",DevExeSpace(),0,(nprtcl_recv-1), KOKKOS_LAMBDA(const int n) {
      int p;
      if (n < nprtcl_send) {
        p = sendlist.d_view(n).prtcl_indx; // place particles in holes created by sends
      } else {
        p = npart + (n - nprtcl_send);     // place particle at end of arrays
      }
      for (int i=0; i<nidata; ++i) {
        pi(i,p) = irecvbuf(nidata*n + i);
      }
      for (int i=0; i<nrdata; ++i) {
        pr(i,p) = rrecvbuf(nrdata*n + i);
      }
      ptag(p) = trecvbuf(n);
    });
  }

  // At this point have filled npart_recv holes in particle arrays from sends
  // If (nprtcl_recv < nprtcl_send), have to move particles from end of arrays to fill
  // remaining holes
  int nremain = nprtcl_send - nprtcl_recv;
  if (nremain > 0) {
    int &npart = pmy_part->nprtcl_thispack;
    int i_last_hole = nprtcl_send-1;
    int i_next_hole = nprtcl_recv;
    for (int n=1; n<=nremain; ++n) {
      int nend = npart-n;
      if (nend > sendlist.h_view(i_last_hole).prtcl_indx) {
        // copy particle from end into hole
        int next_hole = sendlist.h_view(i_next_hole).prtcl_indx;
        auto rdest = Kokkos::subview(pmy_part->prtcl_rdata, Kokkos::ALL, next_hole);
        auto rsrc  = Kokkos::subview(pmy_part->prtcl_rdata, Kokkos::ALL, nend);
        Kokkos::deep_copy(rdest, rsrc);
        auto idest = Kokkos::subview(pmy_part->prtcl_idata, Kokkos::ALL, next_hole);
        auto isrc  = Kokkos::subview(pmy_part->prtcl_idata, Kokkos::ALL, nend);
        Kokkos::deep_copy(idest, isrc);
        auto tdest = Kokkos::subview(pmy_part->prtcl_tag, next_hole);
        auto tsrc = Kokkos::subview(pmy_part->prtcl_tag, nend);
        Kokkos::deep_copy(tdest, tsrc);
        i_next_hole += 1;
      } else {
        // this index contains a hole, so do nothing except find new index of last hole
        i_last_hole -= 1;
      }
    }

    // shrink size of particle data arrays
    Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
    Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
    Kokkos::resize(pmy_part->prtcl_tag, new_npart);
  }

  // Update nparticles_thisrank.  Update cost array (use npart_thismb[nmb]?)
  pmy_part->nprtcl_thispack = new_npart;
#endif
  pmy_part->pmy_pack->pmesh->UpdateParticleCounts();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::ClearPrtclSend()
//! \brief

TaskStatus ParticlesBoundaryValues::ClearPrtclSend() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  // wait for all non-blocking sends for vars to finish before continuing
  for (int n=0; n<nsends; ++n) {
    int ierr = MPI_Wait(&(rsend_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(isend_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(tsend_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rsend_req.clear();
  isend_req.clear();
  tsend_req.clear();
#endif
  nsends=0;
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::ClearPrtclRecv()
//! \brief

TaskStatus ParticlesBoundaryValues::ClearPrtclRecv() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  // wait for all non-blocking receives to finish before continuing
  for (int n=0; n<nrecvs; ++n) {
    int ierr = MPI_Wait(&(rrecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(irecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(trecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rrecv_req.clear();
  irecv_req.clear();
  trecv_req.clear();
#endif
  nrecvs=0;
  return TaskStatus::complete;
}

} // namespace particles
