//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_part.cpp
//! \brief

#include <cstdlib>
#include <cstdint>
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

namespace {

KOKKOS_INLINE_FUNCTION
void WrapParticleCoordinate(Real &x, Real xmin, Real xmax) {
  Real len = xmax - xmin;
  if (len <= 0.0) {return;}
  while (x < xmin) {x += len;}
  while (x >= xmax) {x -= len;}
}

KOKKOS_INLINE_FUNCTION
int ParticleLocationIndex(Real x, Real xmin, Real xmax, int nloc) {
  Real frac = (x - xmin)/(xmax - xmin);
  std::int64_t ix = static_cast<std::int64_t>(frac*static_cast<Real>(nloc));
  ix = (ix < 0) ? 0 : ix;
  ix = (ix > nloc - 1) ? nloc - 1 : ix;
  return static_cast<int>(ix);
}

bool BuildParticleLookupTable(Mesh *pm, int max_entries, bool track_level_transitions,
                              DvceArray1D<int> &gid_table, DvceArray1D<int> &rank_table,
                              DvceArray1D<int> &level_table, int &nloc1, int &nloc2,
                              int &nloc3) {
  int lev_offset = pm->max_level - pm->root_level;
  nloc1 = pm->nmb_rootx1 << lev_offset;
  nloc2 = pm->multi_d ? (pm->nmb_rootx2 << lev_offset) : 1;
  nloc3 = pm->three_d ? (pm->nmb_rootx3 << lev_offset) : 1;

  int64_t ntable = static_cast<int64_t>(nloc1)*static_cast<int64_t>(nloc2)*
                   static_cast<int64_t>(nloc3);
  if (ntable > static_cast<int64_t>(max_entries)) {return false;}

  HostArray1D<int> gid_h("particle_lookup_gid_h", ntable);
  HostArray1D<int> rank_h("particle_lookup_rank_h", ntable);
  Kokkos::deep_copy(gid_h, -1);
  Kokkos::deep_copy(rank_h, -1);
  HostArray1D<int> level_h;
  if (track_level_transitions) {
    level_h = HostArray1D<int>("particle_lookup_level_h", ntable);
    Kokkos::deep_copy(level_h, -1);
  }

  for (int gid=0; gid<pm->nmb_total; ++gid) {
    LogicalLocation &lloc = pm->lloc_eachmb[gid];
    int shift = pm->max_level - lloc.level;
    int ilo = lloc.lx1 << shift;
    int ihi = (lloc.lx1 + 1) << shift;
    int jlo = pm->multi_d ? (lloc.lx2 << shift) : 0;
    int jhi = pm->multi_d ? ((lloc.lx2 + 1) << shift) : 1;
    int klo = pm->three_d ? (lloc.lx3 << shift) : 0;
    int khi = pm->three_d ? ((lloc.lx3 + 1) << shift) : 1;

    for (int k=klo; k<khi; ++k) {
      for (int j=jlo; j<jhi; ++j) {
        for (int i=ilo; i<ihi; ++i) {
          std::int64_t idx = static_cast<std::int64_t>(i) +
                             static_cast<std::int64_t>(nloc1)*
                             (static_cast<std::int64_t>(j) +
                              static_cast<std::int64_t>(nloc2)*
                              static_cast<std::int64_t>(k));
          gid_h(idx) = gid;
          rank_h(idx) = pm->rank_eachmb[gid];
          if (track_level_transitions) {level_h(idx) = lloc.level;}
        }
      }
    }
  }

  for (std::int64_t idx=0; idx<ntable; ++idx) {
    if (gid_h(idx) < 0 || rank_h(idx) < 0 ||
        (track_level_transitions && level_h(idx) < 0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle lookup table has an unfilled logical slot"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  gid_table = DvceArray1D<int>("particle_lookup_gid", ntable);
  rank_table = DvceArray1D<int>("particle_lookup_rank", ntable);
  Kokkos::deep_copy(gid_table, gid_h);
  Kokkos::deep_copy(rank_table, rank_h);
  if (track_level_transitions) {
    level_table = DvceArray1D<int>("particle_lookup_level", ntable);
    Kokkos::deep_copy(level_table, level_h);
  }
  return true;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::UpdateGID()
//! \brief Updates GID of particles that cross boundary of their parent MeshBlock.  If
//! the new GID is on a different rank, then store in sendlist_buf DvceArray: (1) index of
//! particle in prtcl array, (2) destination GID, and (3) destination rank.

KOKKOS_INLINE_FUNCTION
void UpdateGID(int &newgid, NeighborBlock nghbr, int myrank, DvceArray1D<int> counter,
               DualArray1D<ParticleLocationData> slist, int p) {
  newgid = nghbr.gid;
#if MPI_PARALLEL_ENABLED
  if (nghbr.rank != myrank) {
    int index = Kokkos::atomic_fetch_add(&counter(0),1);
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
  Mesh *pm = pmy_part->pmy_pack->pmesh;
  auto gids = pmy_part->pmy_pack->gids;
  auto &pr = pmy_part->prtcl_rdata;
  auto &pi = pmy_part->prtcl_idata;
  int npart = pmy_part->nprtcl_thispack;
  bool track_level_transitions = pmy_part->log_performance;
  int64_t nlevel_transitions = 0;
  nprtcl_send = 0;
  if (npart <= 0) {
    if (pm->multilevel) {
      pmy_part->LogPerformance(
          "particle_multilevel_lookup", 0, 0, 0, 0, nlevel_transitions);
    }
    return TaskStatus::complete;
  }

  if (pm->multilevel) {
    if (pmy_part->amr_remap_mode == ParticlesAMRRemapMode::device_table) {
      DvceArray1D<int> gid_table;
      DvceArray1D<int> rank_table;
      DvceArray1D<int> level_table;
      int nloc1, nloc2, nloc3;
      bool built_table = BuildParticleLookupTable(
          pm, pmy_part->amr_lookup_max_cells, track_level_transitions, gid_table,
          rank_table, level_table, nloc1, nloc2, nloc3);
      if (built_table) {
        Kokkos::realloc(sendlist, std::max(1,npart));
        auto &slist = sendlist;
        DvceArray1D<int> atom_count("particle_lookup_send_count", 3);
        Kokkos::deep_copy(atom_count, 0);
        int myrank = global_variable::my_rank;
        auto &mblev = pmy_part->pmy_pack->pmb->mb_lev;
        int nmb_pack = mblev.d_view.extent_int(0);
        bool multi_d = pm->multi_d;
        bool three_d = pm->three_d;
        Real x1min = pm->mesh_size.x1min;
        Real x1max = pm->mesh_size.x1max;
        Real x2min = pm->mesh_size.x2min;
        Real x2max = pm->mesh_size.x2max;
        Real x3min = pm->mesh_size.x3min;
        Real x3max = pm->mesh_size.x3max;

        par_for("part_update_lookup",DevExeSpace(),0,(npart-1),
        KOKKOS_LAMBDA(const int p) {
          Real x1 = pr(IPX,p);
          Real x2 = pr(IPY,p);
          Real x3 = pr(IPZ,p);
          WrapParticleCoordinate(x1, x1min, x1max);
          if (multi_d) {WrapParticleCoordinate(x2, x2min, x2max);}
          if (three_d) {WrapParticleCoordinate(x3, x3min, x3max);}

          int ix = ParticleLocationIndex(x1, x1min, x1max, nloc1);
          int iy = multi_d ? ParticleLocationIndex(x2, x2min, x2max, nloc2) : 0;
          int iz = three_d ? ParticleLocationIndex(x3, x3min, x3max, nloc3) : 0;
          std::int64_t idx = static_cast<std::int64_t>(ix) +
                             static_cast<std::int64_t>(nloc1)*
                             (static_cast<std::int64_t>(iy) +
                              static_cast<std::int64_t>(nloc2)*
                              static_cast<std::int64_t>(iz));
          int dest_gid = gid_table(idx);
          int dest_rank = rank_table(idx);
          if (dest_gid < 0 || dest_rank < 0) {
            Kokkos::atomic_fetch_add(&atom_count(1), 1);
            return;
          }

          if (track_level_transitions) {
            int source_pack_index = pi(PGID,p) - gids;
            if (source_pack_index < 0 || source_pack_index >= nmb_pack) {
              Kokkos::atomic_fetch_add(&atom_count(1), 1);
              return;
            }
            int source_level = mblev.d_view(source_pack_index);
            if (level_table(idx) != source_level) {
              Kokkos::atomic_fetch_add(&atom_count(2), 1);
            }
          }
          pr(IPX,p) = x1;
          pr(IPY,p) = x2;
          pr(IPZ,p) = x3;
          pi(PGID,p) = dest_gid;
#if MPI_PARALLEL_ENABLED
          if (dest_rank != myrank) {
            int index = Kokkos::atomic_fetch_add(&atom_count(0), 1);
            slist.d_view(index).prtcl_indx = p;
            slist.d_view(index).dest_gid = dest_gid;
            slist.d_view(index).dest_rank = dest_rank;
          }
#endif
        });
        HostArray1D<int> count_h("particle_lookup_send_count_h", 3);
        Kokkos::deep_copy(count_h, atom_count);
        if (count_h(1) != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << count_h(1)
                    << " particles failed multilevel lookup validation"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        nlevel_transitions = count_h(2);
        nprtcl_send = count_h(0);
        Kokkos::resize(sendlist, nprtcl_send);
        sendlist.template modify<DevExeSpace>();
        sendlist.template sync<HostMemSpace>();
        pmy_part->LogPerformance(
            "particle_multilevel_lookup", 0, 0, 0, 0, nlevel_transitions);
        return TaskStatus::complete;
      }
    }

    auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr);
    auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pi);
    std::vector<ParticleLocationData> host_sendlist;
    host_sendlist.reserve(npart);
    int myrank = global_variable::my_rank;
    for (int p=0; p<npart; ++p) {
      WrapParticleCoordinate(pr_h(IPX,p), pm->mesh_size.x1min, pm->mesh_size.x1max);
      if (pm->multi_d) {
        WrapParticleCoordinate(pr_h(IPY,p), pm->mesh_size.x2min, pm->mesh_size.x2max);
      }
      if (pm->three_d) {
        WrapParticleCoordinate(pr_h(IPZ,p), pm->mesh_size.x3min, pm->mesh_size.x3max);
      }
      int dest_gid = pm->FindMeshBlockContainingPosition(pr_h(IPX,p), pr_h(IPY,p),
                                                         pr_h(IPZ,p));
      if (dest_gid < 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle host lookup failed for particle " << p
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (track_level_transitions) {
        int source_gid = pi_h(PGID,p);
        if (source_gid < 0 || source_gid >= pm->nmb_total) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Particle source GID is outside the mesh hierarchy"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (pm->lloc_eachmb[dest_gid].level != pm->lloc_eachmb[source_gid].level) {
          ++nlevel_transitions;
        }
      }
      pi_h(PGID,p) = dest_gid;
#if MPI_PARALLEL_ENABLED
      int dest_rank = pm->rank_eachmb[dest_gid];
      if (dest_rank != myrank) {
        ParticleLocationData loc;
        loc.prtcl_indx = p;
        loc.dest_gid = dest_gid;
        loc.dest_rank = dest_rank;
        host_sendlist.push_back(loc);
      }
#endif
    }
    Kokkos::deep_copy(pr, pr_h);
    Kokkos::deep_copy(pi, pi_h);
    nprtcl_send = static_cast<int>(host_sendlist.size());
    Kokkos::realloc(sendlist, std::max(1,nprtcl_send));
    for (int n=0; n<nprtcl_send; ++n) {
      sendlist.h_view(n) = host_sendlist[n];
    }
    Kokkos::resize(sendlist, nprtcl_send);
    sendlist.template modify<HostMemSpace>();
    sendlist.template sync<DevExeSpace>();
    pmy_part->LogPerformance(
        "particle_multilevel_lookup", 0, 0, 0, 0, nlevel_transitions);
    return TaskStatus::complete;
  }

  auto &mbsize = pmy_part->pmy_pack->pmb->mb_size;
  auto &mblev = pmy_part->pmy_pack->pmb->mb_lev;
  auto &meshsize = pm->mesh_size;
  auto myrank = global_variable::my_rank;
  auto &nghbr = pmy_part->pmy_pack->pmb->nghbr;
  auto &psendl = sendlist;
  DvceArray1D<int> atom_count("particle_send_count",1);
  Kokkos::deep_copy(atom_count, 0);
  bool &multi_d = pmy_part->pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_part->pmy_pack->pmesh->three_d;

  Kokkos::realloc(sendlist, std::max(1,npart));
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
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        } else if (ix == 0) {
          // x2 face
          int indx = NeighborIndex(0,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,0,fx,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        } else {
          // x1x2 edge
          int indx = NeighborIndex(ix,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,iy,0,fz,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        }
      } else if (iy == 0) {
        if (ix == 0) {
          // x3 face
          int indx = NeighborIndex(0,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,0,iz,fx,fy);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        } else {
          // x3x1 edge
          int indx = NeighborIndex(ix,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,0,iz,fy,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        }
      } else {
        if (ix == 0) {
          // x2x3 edge
          int indx = NeighborIndex(0,iy,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,iz,fx,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
        } else {
          // corners
          int indx = NeighborIndex(ix,iy,iz,0,0);
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank, atom_count, psendl, p);
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
  int counter = 0;
  auto atom_count_h = Kokkos::create_mirror_view(atom_count);
  Kokkos::deep_copy(atom_count_h, atom_count);
  counter = atom_count_h(0);
  nprtcl_send = counter;
  Kokkos::resize(sendlist, nprtcl_send);
  // sync sendlist device array with host
  sendlist.template modify<DevExeSpace>();
  sendlist.template sync<HostMemSpace>();

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::EnsureBufferCapacity()
//! \brief Grow particle MPI buffers only when the current capacity is too small.

void ParticlesBoundaryValues::EnsureBufferCapacity(int nsend, int nrecv) {
#if MPI_PARALLEL_ENABLED
  int rsend_needed = std::max(1, (pmy_part->nrdata)*nsend);
  int isend_needed = std::max(1, (pmy_part->nidata)*nsend);
  int rrecv_needed = std::max(1, (pmy_part->nrdata)*nrecv);
  int irecv_needed = std::max(1, (pmy_part->nidata)*nrecv);
  if (rsend_needed > rsend_capacity) {
    Kokkos::realloc(prtcl_rsendbuf, rsend_needed);
    rsend_capacity = rsend_needed;
  }
  if (isend_needed > isend_capacity) {
    Kokkos::realloc(prtcl_isendbuf, isend_needed);
    isend_capacity = isend_needed;
  }
  if (rrecv_needed > rrecv_capacity) {
    Kokkos::realloc(prtcl_rrecvbuf, rrecv_needed);
    rrecv_capacity = rrecv_needed;
  }
  if (irecv_needed > irecv_capacity) {
    Kokkos::realloc(prtcl_irecvbuf, irecv_needed);
    irecv_capacity = irecv_needed;
  }
#endif
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

  if (pmy_part->exchange_mode == ParticlesExchangeMode::alltoall_counts) {
    std::fill(send_counts_eachrank.begin(), send_counts_eachrank.end(), 0);
    std::fill(recv_counts_eachrank.begin(), recv_counts_eachrank.end(), 0);
    for (int n=0; n<nsends; ++n) {
      send_counts_eachrank[sends_thisrank[n].recvrank] = sends_thisrank[n].nprtcls;
    }
    MPI_Alltoall(send_counts_eachrank.data(), 1, MPI_INT,
                 recv_counts_eachrank.data(), 1, MPI_INT, mpi_comm_part);

    recvs_thisrank.clear();
    int myrank = global_variable::my_rank;
    for (int rank=0; rank<global_variable::nranks; ++rank) {
      if (recv_counts_eachrank[rank] > 0) {
        recvs_thisrank.emplace_back(ParticleMessageData(rank, myrank,
                                                        recv_counts_eachrank[rank]));
      }
    }
    nrecvs = recvs_thisrank.size();
    sends_allranks.clear();
    return TaskStatus::complete;
  }

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
  nprtcl_recv=0;
#if MPI_PARALLEL_ENABLED
  if (pmy_part->exchange_mode == ParticlesExchangeMode::allgather) {
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
  }

  // Figure out how many particles will be received from all ranks
  for (int n=0; n<nrecvs; ++n) {
    nprtcl_recv += recvs_thisrank[n].nprtcls;
  }

  EnsureBufferCapacity(nprtcl_send, nprtcl_recv);

  // Post non-blocking receives
  bool no_errors=true;
  rrecv_req.clear();
  irecv_req.clear();
  for (int n=0; n<nrecvs; ++n) {
    rrecv_req.emplace_back(MPI_REQUEST_NULL);
    irecv_req.emplace_back(MPI_REQUEST_NULL);
  }

  // Init receives for Reals
  int data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    // calculate amount of data to be passed, get pointer to variables
    int data_size = (pmy_part->nrdata)*(recvs_thisrank[n].nprtcls);
    int data_end = data_start + (pmy_part->nrdata)*recvs_thisrank[n].nprtcls;
    auto recv_ptr = Kokkos::subview(prtcl_rrecvbuf, std::make_pair(data_start, data_end));
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
    int data_end = data_start + (pmy_part->nidata)*recvs_thisrank[n].nprtcls;
    auto recv_ptr = Kokkos::subview(prtcl_irecvbuf, std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 1; // 0 for Reals, 1 for ints

    // Post non-blocking receive
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_INT, drank, tag,
                         mpi_comm_part, &(irecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start += data_size;
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
  EnsureBufferCapacity(nprtcl_send, nprtcl_recv);
  if (nprtcl_send > 0) {
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
      int data_end = data_start + nrdata*sends_thisrank[n].nprtcls;
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
      int data_end = data_start + nidata*sends_thisrank[n].nprtcls;
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
  int old_npart = pmy_part->nprtcl_thispack;
  int new_npart = old_npart + (nprtcl_recv - nprtcl_send);
  if (nprtcl_recv > nprtcl_send) {
    Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
    Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
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
    auto &rrecvbuf = prtcl_rrecvbuf;
    auto &irecvbuf = prtcl_irecvbuf;
    int npart = old_npart;
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
    });
  }

  // At this point have filled npart_recv holes in particle arrays from sends
  // If (nprtcl_recv < nprtcl_send), have to move particles from end of arrays to fill
  // remaining holes
  int nremain = nprtcl_send - nprtcl_recv;
  if (nremain > 0) {
    DvceArray1D<int> drop("particle_drop_mask", old_npart);
    Kokkos::deep_copy(drop, 0);
    auto &slist = sendlist;
    par_for("particle_mark_drops",DevExeSpace(),0,(nremain-1),
    KOKKOS_LAMBDA(const int n) {
      int send_index = nprtcl_recv + n;
      drop(slist.d_view(send_index).prtcl_indx) = 1;
    });

    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto old_pr = pmy_part->prtcl_rdata;
    auto old_pi = pmy_part->prtcl_idata;
    DvceArray2D<Real> new_pr("particle_compact_rdata", nrdata, new_npart);
    DvceArray2D<int> new_pi("particle_compact_idata", nidata, new_npart);
    Kokkos::parallel_scan("particle_compact",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, old_npart),
      KOKKOS_LAMBDA(const int p, int &offset, const bool final) {
        if (drop(p) == 0) {
          if (final) {
            for (int i=0; i<nidata; ++i) {
              new_pi(i,offset) = old_pi(i,p);
            }
            for (int i=0; i<nrdata; ++i) {
              new_pr(i,offset) = old_pr(i,p);
            }
          }
          offset += 1;
        }
      });
    pmy_part->prtcl_idata = new_pi;
    pmy_part->prtcl_rdata = new_pr;
  }

  // Update nparticles_thisrank.  Update cost array (use npart_thismb[nmb]?)
  pmy_part->nprtcl_thispack = new_npart;
  pmy_part->pmy_pack->pmesh->UpdateParticleCounts();
  int64_t bytes_sent = static_cast<int64_t>(nprtcl_send)*
                       (static_cast<int64_t>(pmy_part->nrdata)*sizeof(Real) +
                        static_cast<int64_t>(pmy_part->nidata)*sizeof(int));
  pmy_part->LogPerformance("particle_mpi_exchange", 0, nprtcl_send, nprtcl_recv,
                           0, 0, static_cast<int64_t>(2*nsends), bytes_sent);
  pmy_part->CheckConsistency("particle MPI exchange");
#endif
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
