//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_part.cpp
//! \brief

#include <cstdlib>
#include <iostream>
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
//! \fn void ParticlesBoundaryValues::MarkForDestruction()
//! \brief Marks particles that leave a non-periodic mesh boundary for removal. The list
//! uses ParticleLocationData so it can be sorted and compacted with MPI send holes.

KOKKOS_INLINE_FUNCTION
void MarkForDestruction(int *pcounter, DualArray1D<ParticleLocationData> dlist, int p) {
  int index = Kokkos::atomic_fetch_add(pcounter, 1);
  dlist.d_view(index).prtcl_indx = p;
  dlist.d_view(index).dest_gid = 0;
  dlist.d_view(index).dest_rank = 0;
  return;
}

KOKKOS_INLINE_FUNCTION
bool IsPeriodicParticleBoundary(BoundaryFlag flag) {
  return (flag == BoundaryFlag::periodic || flag == BoundaryFlag::shear_periodic);
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
  auto &pdestroyl = destroylist;
  bool &multi_d = pmy_part->pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_part->pmy_pack->pmesh->three_d;
  auto &mb_bcs = pmy_part->pmy_pack->pmb->mb_bcs;

  Kokkos::View<int> atom_count("atom_count");
  Kokkos::View<int> atom_destroy_count("atom_destroy_count");
  Kokkos::deep_copy(atom_count, 0);
  Kokkos::deep_copy(atom_destroy_count, 0);

  Kokkos::realloc(sendlist, npart);
  Kokkos::realloc(destroylist, npart);
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
      bool nonperiodic_mesh_exit =
          (x1 < meshsize.x1min &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::inner_x1))) ||
          (x1 > meshsize.x1max &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::outer_x1))) ||
          (x2 < meshsize.x2min &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::inner_x2))) ||
          (x2 > meshsize.x2max &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::outer_x2))) ||
          (x3 < meshsize.x3min &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::inner_x3))) ||
          (x3 > meshsize.x3max &&
           !IsPeriodicParticleBoundary(mb_bcs.d_view(m,BoundaryFace::outer_x3)));

      if (nonperiodic_mesh_exit) {
        MarkForDestruction(&atom_destroy_count(), pdestroyl, p);
      } else if (iz == 0) {
        if (iy == 0) {
          // x1 face
          int indx = NeighborIndex(ix,0,0,0,0);           // neighbor at same level
          if (nghbr.d_view(m,indx).lev > mylevel) {       // neighbor at finer level
            indx = NeighborIndex(ix,0,0,fy,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}  // neighbor at coarser level
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        } else if (ix == 0) {
          // x2 face
          int indx = NeighborIndex(0,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,0,fx,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        } else {
          // x1x2 edge
          int indx = NeighborIndex(ix,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,iy,0,fz,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        }
      } else if (iy == 0) {
        if (ix == 0) {
          // x3 face
          int indx = NeighborIndex(0,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,0,iz,fx,fy);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        } else {
          // x3x1 edge
          int indx = NeighborIndex(ix,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,0,iz,fy,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        }
      } else {
        if (ix == 0) {
          // x2x3 edge
          int indx = NeighborIndex(0,iy,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,iz,fx,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        } else {
          // corners
          int indx = NeighborIndex(ix,iy,iz,0,0);
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, &atom_count(), psendl, p);
        }
      }

      // reset x,y,z positions if particle crosses Mesh boundary using periodic BCs
      if (!nonperiodic_mesh_exit && x1 < meshsize.x1min) {
        pr(IPX,p) += (meshsize.x1max - meshsize.x1min);
      } else if (!nonperiodic_mesh_exit && x1 > meshsize.x1max) {
        pr(IPX,p) -= (meshsize.x1max - meshsize.x1min);
      }
      if (!nonperiodic_mesh_exit && x2 < meshsize.x2min) {
        pr(IPY,p) += (meshsize.x2max - meshsize.x2min);
      } else if (!nonperiodic_mesh_exit && x2 > meshsize.x2max) {
        pr(IPY,p) -= (meshsize.x2max - meshsize.x2min);
      }
      if (!nonperiodic_mesh_exit && x3 < meshsize.x3min) {
        pr(IPZ,p) += (meshsize.x3max - meshsize.x3min);
      } else if (!nonperiodic_mesh_exit && x3 > meshsize.x3max) {
        pr(IPZ,p) -= (meshsize.x3max - meshsize.x3min);
      }
    }
  });
  Kokkos::fence();
  Kokkos::deep_copy(nprtcl_send, atom_count);
  Kokkos::deep_copy(nprtcl_destroy, atom_destroy_count);
  Kokkos::resize(sendlist, nprtcl_send);
  Kokkos::resize(destroylist, nprtcl_destroy);
  // sync sendlist device array with host
  sendlist.template modify<DevExeSpace>();
  sendlist.template sync<HostMemSpace>();
  destroylist.template modify<DevExeSpace>();
  destroylist.template sync<HostMemSpace>();

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
  nsends = sends_thisrank.size();

  // Share number of ranks to send to amongst all ranks
  nsends_eachrank[global_variable::my_rank] = nsends;
  MPI_Allgather(&nsends, 1, MPI_INT, nsends_eachrank.data(), 1, MPI_INT, mpi_comm_part);

  // Now share ParticleMessageData amongst all ranks
  // First create vector of starting indices in full vector
  std::vector<int> nsends_displ;
  nsends_displ.resize(global_variable::nranks);
  nsends_displ[0] = 0;
  for (int n=1; n<(global_variable::nranks); ++n) {
    nsends_displ[n] = nsends_displ[n-1] + nsends_eachrank[n-1];
  }
  int nsends_allranks = nsends_displ[global_variable::nranks - 1] +
                        nsends_eachrank[global_variable::nranks - 1];
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
  MPI_Type_free(&mpi_ituple);
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::InitPrtclRecv()
//! \brief

TaskStatus ParticlesBoundaryValues::InitPrtclRecv() {
  nprtcl_recv = 0;
#if MPI_PARALLEL_ENABLED
  // load STL::vector of ParticleMessageData with <sendrank,recvrank,nprtcl_recv> for
  // receives // on this rank. Length will be nrecvs, initially this length is unknown
  recvs_thisrank.clear();

  int nsends_allranks = sends_allranks.size();
  for (int n=0; n<nsends_allranks; ++n) {
    if (sends_allranks[n].recvrank == global_variable::my_rank) {
      recvs_thisrank.emplace_back(sends_allranks[n]);
    }
  }
  nrecvs = recvs_thisrank.size();

  // Figure out how many particles will be received from all ranks
  for (int n=0; n<nrecvs; ++n) {
    nprtcl_recv += recvs_thisrank[n].nprtcls;
  }

  // Post non-blocking receives
  bool no_errors=true;
  rrecv_req.clear();
  irecv_req.clear();
  if (nprtcl_recv > 0) {
    // Allocate receive buffer
    Kokkos::realloc(prtcl_rrecvbuf, (pmy_part->nrdata)*nprtcl_recv);
    Kokkos::realloc(prtcl_irecvbuf, (pmy_part->nidata)*nprtcl_recv);

    for (int n=0; n<nrecvs; ++n) {
      rrecv_req.emplace_back(MPI_REQUEST_NULL);
      irecv_req.emplace_back(MPI_REQUEST_NULL);
    }

    // Init receives for Reals
    int data_start=0;
    for (int n=0; n<nrecvs; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = (pmy_part->nrdata)*(recvs_thisrank[n].nprtcls);
      int data_end = data_start + data_size;
      auto recv_ptr = Kokkos::subview(prtcl_rrecvbuf,
                                      std::make_pair(data_start, data_end));
      int drank = recvs_thisrank[n].sendrank;
      int tag = 0; // 0 for Reals, 1 for ints

      // Post non-blocking receive
      int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                           mpi_comm_part, &(rrecv_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
    }
    // Init receives for ints
    data_start=0;
    for (int n=0; n<nrecvs; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = (pmy_part->nidata)*(recvs_thisrank[n].nprtcls);
      int data_end = data_start + data_size;
      auto recv_ptr = Kokkos::subview(prtcl_irecvbuf,
                                      std::make_pair(data_start, data_end));
      int drank = recvs_thisrank[n].sendrank;
      int tag = 1; // 0 for Reals, 1 for ints

      // Post non-blocking receive
      int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_INT, drank, tag,
                           mpi_comm_part, &(irecv_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
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
//! \fn void ParticlesBoundaryValues::PackAndSendPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::PackAndSendPrtcls() {
#if MPI_PARALLEL_ENABLED
  // Figure out how many particles will be sent from this ranks
  nprtcl_send=0;
  for (int n=0; n<nsends; ++n) {
    nprtcl_send += sends_thisrank[n].nprtcls;
  }

  bool no_errors=true;
  if (nprtcl_send > 0) {
    // Allocate send buffer
    Kokkos::realloc(prtcl_rsendbuf, (pmy_part->nrdata)*nprtcl_send);
    Kokkos::realloc(prtcl_isendbuf, (pmy_part->nidata)*nprtcl_send);

    // sendlist on device is already sorted by destrank in CountSendAndRecvs()
    // Use sendlist on device to load particles into send buffer ordered by dest_rank
    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto &pr = pmy_part->prtcl_rdata;
    auto &pi = pmy_part->prtcl_idata;
    auto &rsendbuf = prtcl_rsendbuf;
    auto &isendbuf = prtcl_isendbuf;
    par_for("ppack",DevExeSpace(),0,(nprtcl_send-1), KOKKOS_LAMBDA(const int n) {
      int p = sendlist.d_view(n).prtcl_indx;
      for (int i=0; i<nidata; ++i) {
        isendbuf(nidata*n + i) = pi(i,p);
      }
      for (int i=0; i<nrdata; ++i) {
        rsendbuf(nrdata*n + i) = pr(i,p);
      }
    });

    // Post non-blocking sends
    Kokkos::fence();
    rsend_req.clear();
    isend_req.clear();
    for (int n=0; n<nsends; ++n) {
      rsend_req.emplace_back(MPI_REQUEST_NULL);
      isend_req.emplace_back(MPI_REQUEST_NULL);
    }

    // Send Reals
    int data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = nrdata*(sends_thisrank[n].nprtcls);
      int data_end = data_start + data_size;
      auto send_ptr = Kokkos::subview(prtcl_rsendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 0; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                           mpi_comm_part, &(rsend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
    }
    // Send ints
    data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = nidata*(sends_thisrank[n].nprtcls);
      int data_end = data_start + data_size;
      auto send_ptr = Kokkos::subview(prtcl_isendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 1; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_INT, drank, tag,
                           mpi_comm_part, &(isend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
    }
  }

  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in posting non-blocking sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::RecvAndUnpackPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::RecvAndUnpackPrtcls() {
  namespace KE = Kokkos::Experimental;
  if (nprtcl_send > 0) {
    std::sort(KE::begin(sendlist.h_view), KE::end(sendlist.h_view), SortByIndex);
    sendlist.template modify<HostMemSpace>();
  }
  if (nprtcl_destroy > 0) {
    std::sort(KE::begin(destroylist.h_view), KE::end(destroylist.h_view), SortByIndex);
    destroylist.template modify<HostMemSpace>();
  }

  int local_change = nprtcl_send + nprtcl_recv + nprtcl_destroy;
  int global_change = local_change;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&local_change, &global_change, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // wait for particle communications to complete before global count updates
  bool no_errors=true;
  for (int n=0; n<nrecvs; ++n) {
    int ierr = MPI_Wait(&(rrecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(irecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in receiving particles"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif

  if (global_change == 0) {
    return TaskStatus::complete;
  }

  int old_npart = pmy_part->nprtcl_thispack;
  int new_npart = old_npart + nprtcl_recv - nprtcl_send - nprtcl_destroy;
  if (new_npart < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle boundary exchange produced a negative "
              << "particle count" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto old_r = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_part->prtcl_rdata);
  auto old_i = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_part->prtcl_idata);
  std::vector<char> remove_particle(old_npart, 0);
  int nmarked = 0;
  auto mark_removed = [&](int pindx) {
    if (pindx < 0 || pindx >= old_npart || remove_particle[pindx] != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Invalid particle removal index " << pindx
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    remove_particle[pindx] = 1;
    ++nmarked;
  };
  for (int n=0; n<nprtcl_send; ++n) {
    mark_removed(sendlist.h_view(n).prtcl_indx);
  }
  for (int n=0; n<nprtcl_destroy; ++n) {
    mark_removed(destroylist.h_view(n).prtcl_indx);
  }

  Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
  Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
  auto new_r = Kokkos::create_mirror_view(pmy_part->prtcl_rdata);
  auto new_i = Kokkos::create_mirror_view(pmy_part->prtcl_idata);

  int next_particle = 0;
  for (int p=0; p<old_npart; ++p) {
    if (remove_particle[p] == 0) {
      for (int i=0; i<pmy_part->nidata; ++i) {
        new_i(i,next_particle) = old_i(i,p);
      }
      for (int i=0; i<pmy_part->nrdata; ++i) {
        new_r(i,next_particle) = old_r(i,p);
      }
      ++next_particle;
    }
  }

#if MPI_PARALLEL_ENABLED
  if (nprtcl_recv > 0) {
    auto recv_r = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rrecvbuf);
    auto recv_i = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_irecvbuf);
    for (int n=0; n<nprtcl_recv; ++n) {
      for (int i=0; i<pmy_part->nidata; ++i) {
        new_i(i,next_particle) = recv_i((pmy_part->nidata)*n + i);
      }
      for (int i=0; i<pmy_part->nrdata; ++i) {
        new_r(i,next_particle) = recv_r((pmy_part->nrdata)*n + i);
      }
      ++next_particle;
    }
  }
#endif

  if (next_particle != new_npart || nmarked != nprtcl_send + nprtcl_destroy) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle compaction count mismatch" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Kokkos::deep_copy(pmy_part->prtcl_idata, new_i);
  Kokkos::deep_copy(pmy_part->prtcl_rdata, new_r);
  pmy_part->nprtcl_thispack = new_npart;
  pmy_part->UpdateParticleCounts();
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
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rsend_req.clear();
  isend_req.clear();
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
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rrecv_req.clear();
  irecv_req.clear();
#endif
  nrecvs=0;
  return TaskStatus::complete;
}

} // namespace particles
