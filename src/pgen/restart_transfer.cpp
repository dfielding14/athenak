//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_transfer.cpp
//! \brief Transfer source restart patches to a reblocked or refined uniform mesh.

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/restart_prolongation.hpp"
#include "mhd/mhd.hpp"
#include "pgen/restart_transfer.hpp"

namespace {

struct Segment {
  int block, source, destination, count;
};

struct Patch {
  int gid, field, block;
  int source[3], destination[3], count[3];
};

[[noreturn]] void Fail(const char *message) {
  std::cerr << "### FATAL ERROR: restart transfer: " << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

int MPIInt(std::size_t count) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    Fail("MPI staging exceeds the supported count");
  }
  return static_cast<int>(count);
}

std::vector<int> Offsets(const std::vector<int> &counts) {
  std::vector<int> offsets(counts.size() + 1, 0);
  for (std::size_t n = 0; n < counts.size(); ++n) {
    offsets[n+1] = MPIInt(static_cast<std::size_t>(offsets[n]) + counts[n]);
  }
  return offsets;
}

std::vector<Segment> SplitAxis(int start, int length, int cells, int blocks,
                               int source_ng, int destination_start) {
  std::vector<Segment> segments;
  const int total = cells*blocks;
  for (int n = 0; n < length;) {
    const int position = ((start + n)%total + total)%total;
    const int local = position%cells;
    const int count = std::min(length - n, cells - local);
    segments.push_back({position/cells, local + source_ng, n + destination_start, count});
    n += count;
  }
  return segments;
}

int Values(const Patch &patch, int nvar) {
  const int variables = patch.field == 0 ? nvar : 1;
  return MPIInt(static_cast<std::size_t>(variables)*patch.count[0]*patch.count[1]*
                patch.count[2]);
}

}  // namespace

void ReadRestartMHD(Mesh *pm, ParameterInput *pin, IOWrapper &file,
                    IOWrapperSizeT offset, IOWrapperSizeT block_bytes) {
  const auto *pack = pm->pmb_pack;
  if (pack->pmhd == nullptr || pack->phydro != nullptr || pack->prad != nullptr ||
      pack->pz4c != nullptr || pack->padm != nullptr || pack->pdyngr != nullptr ||
      pack->pnr != nullptr || pack->pionn != nullptr || pack->ppart != nullptr ||
      pack->pcoord->is_special_relativistic || pack->pcoord->is_general_relativistic ||
      pack->pcoord->is_dynamical_relativistic) {
    Fail("requires standalone Newtonian MHD without particles or coupled physics");
  }
  const auto &source = *pm->restart_transfer;
  const auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const auto &source_indcs = source.source_mb_indcs;
  const int source_ng = source_indcs.ng;
  const int nx = source_indcs.nx1, ny = source_indcs.nx2, nz = source_indcs.nx3;
  const int factor = source.refine ? 2 : 1;
  const int halo = source.refine ? 1 : 0;
  const int cx = indcs.nx1/factor, cy = indcs.nx2/factor, cz = indcs.nx3/factor;
  const int nvar = pm->pmb_pack->pmhd->nmhd + pm->pmb_pack->pmhd->nscalars;
  const int nlocal = pm->pmb_pack->nmb_thispack;
  const int nsource = static_cast<int>(source.source_gid_by_location.size());
  const int rank = global_variable::my_rank, nranks = global_variable::nranks;
  const int cells[3] = {nx, ny, nz};
  const int patch_cells[3] = {cx, cy, cz};
  const int roots[3] = {source.nmb_rootx1, source.nmb_rootx2, source.nmb_rootx3};

  // A shared face belongs to the block on its positive side, including periodic faces.
  std::vector<std::vector<Patch>> requests(nranks);
  for (int m = 0; m < nlocal; ++m) {
    const auto &location = pm->lloc_eachmb[pm->gids_eachrank[rank] + m];
    const int position[3] = {location.lx1, location.lx2, location.lx3};
    for (int field = 0; field < 4; ++field) {
      std::array<std::vector<Segment>, 3> segments;
      for (int axis = 0; axis < 3; ++axis) {
        segments[axis] = SplitAxis(position[axis]*patch_cells[axis] - halo,
            patch_cells[axis] + 2*halo + (field == axis + 1), cells[axis],
            roots[axis], source_ng, ng - halo);
      }
      for (const auto &z : segments[2]) {
        for (const auto &y : segments[1]) {
          for (const auto &x : segments[0]) {
            const int index = (z.block*roots[1] + y.block)*roots[0] + x.block;
            const int gid = source.source_gid_by_location[index];
            requests[gid%nranks].push_back({gid, field, m,
                {x.source, y.source, z.source},
                {x.destination, y.destination, z.destination},
                {x.count, y.count, z.count}});
          }
        }
      }
    }
  }

  std::vector<int> request_send(nranks), request_recv(nranks);
  for (int n = 0; n < nranks; ++n) request_send[n] = MPIInt(requests[n].size());
#if MPI_PARALLEL_ENABLED
  MPI_Alltoall(request_send.data(), 1, MPI_INT, request_recv.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
#else
  request_recv = request_send;
#endif
  const auto send_offset = Offsets(request_send), recv_offset = Offsets(request_recv);
  std::vector<Patch> outgoing(send_offset.back()), incoming(recv_offset.back());
  for (int n = 0; n < nranks; ++n) {
    std::copy(requests[n].begin(), requests[n].end(), outgoing.begin() + send_offset[n]);
  }
#if MPI_PARALLEL_ENABLED
  std::vector<int> send_bytes(nranks), recv_bytes(nranks);
  for (int n = 0; n < nranks; ++n) {
    send_bytes[n] = MPIInt(static_cast<std::size_t>(request_send[n])*sizeof(Patch));
    recv_bytes[n] = MPIInt(static_cast<std::size_t>(request_recv[n])*sizeof(Patch));
  }
  const auto send_byte_offset = Offsets(send_bytes);
  const auto recv_byte_offset = Offsets(recv_bytes);
  MPI_Alltoallv(outgoing.data(), send_bytes.data(), send_byte_offset.data(), MPI_BYTE,
                incoming.data(), recv_bytes.data(), recv_byte_offset.data(), MPI_BYTE,
                MPI_COMM_WORLD);
#else
  incoming = outgoing;
#endif

  std::vector<int> value_send(nranks), value_recv(nranks);
  std::vector<int> patch_offset(incoming.size());
  for (int n = 0; n < nranks; ++n) {
    for (int q = recv_offset[n]; q < recv_offset[n+1]; ++q) {
      value_send[n] = MPIInt(static_cast<std::size_t>(value_send[n]) +
                            Values(incoming[q], nvar));
    }
    for (int q = send_offset[n]; q < send_offset[n+1]; ++q) {
      value_recv[n] = MPIInt(static_cast<std::size_t>(value_recv[n]) +
                            Values(outgoing[q], nvar));
    }
  }
  const auto value_send_offset = Offsets(value_send);
  const auto value_recv_offset = Offsets(value_recv);
  for (int n = 0; n < nranks; ++n) {
    int position = value_send_offset[n];
    for (int q = recv_offset[n]; q < recv_offset[n+1]; ++q) {
      patch_offset[q] = position;
      position += Values(incoming[q], nvar);
    }
  }
  std::vector<Real> send_values(value_send_offset.back());
  std::vector<Real> recv_values(value_recv_offset.back());
  const int sx = nx + 2*source_ng, sy = ny + 2*source_ng, sz = nz + 2*source_ng;
  const std::size_t field_offset[4] = {0, static_cast<std::size_t>(nvar)*sz*sy*sx,
      static_cast<std::size_t>(nvar)*sz*sy*sx + static_cast<std::size_t>(sz)*sy*(sx+1),
      static_cast<std::size_t>(nvar)*sz*sy*sx + static_cast<std::size_t>(sz)*sy*(sx+1) +
          static_cast<std::size_t>(sz)*(sy+1)*sx};
  std::vector<Real> parent;
  const int parent_values = MPIInt(block_bytes/sizeof(Real));
  if (rank < nsource) {
    parent.resize(parent_values);
  }
  for (int gid = rank; gid < nsource; gid += nranks) {
    if (file.Read_Reals_at(parent.data(), parent.size(), offset + gid*block_bytes) !=
        parent.size()) {
      Fail("could not read a complete source MeshBlock");
    }
    for (std::size_t q = 0; q < incoming.size(); ++q) {
      const auto &patch = incoming[q];
      if (patch.gid != gid) continue;
      int position = patch_offset[q];
      const int ni = sx + (patch.field == 1), nj = sy + (patch.field == 2);
      const int nk = sz + (patch.field == 3);
      for (int v = 0; v < (patch.field == 0 ? nvar : 1); ++v) {
        for (int k = 0; k < patch.count[2]; ++k) {
          for (int j = 0; j < patch.count[1]; ++j) {
            const auto start = field_offset[patch.field] +
                ((static_cast<std::size_t>(v)*nk + patch.source[2] + k)*nj +
                 patch.source[1] + j)*ni + patch.source[0];
            std::copy_n(parent.begin() + start, patch.count[0],
                        send_values.begin() + position);
            position += patch.count[0];
          }
        }
      }
    }
  }
  parent.clear();
  parent.shrink_to_fit();
#if MPI_PARALLEL_ENABLED
  MPI_Alltoallv(send_values.data(), value_send.data(), value_send_offset.data(),
                MPI_ATHENA_REAL, recv_values.data(), value_recv.data(),
                value_recv_offset.data(), MPI_ATHENA_REAL, MPI_COMM_WORLD);
#else
  recv_values = send_values;
#endif
  send_values.clear();
  send_values.shrink_to_fit();

  const int hx = cx + 2*ng, hy = cy + 2*ng, hz = cz + 2*ng;
  HostArray5D<Real> staged_cons("restart-transfer-cons", nlocal, nvar, hz, hy, hx);
  HostFaceFld4D<Real> staged_b("restart-transfer-b", nlocal, hz, hy, hx);
  for (int n = 0; n < nranks; ++n) {
    int position = value_recv_offset[n];
    for (int q = send_offset[n]; q < send_offset[n+1]; ++q) {
      const auto &patch = outgoing[q];
      for (int v = 0; v < (patch.field == 0 ? nvar : 1); ++v) {
        for (int k = 0; k < patch.count[2]; ++k) {
          for (int j = 0; j < patch.count[1]; ++j) {
            for (int i = 0; i < patch.count[0]; ++i) {
              const int di = patch.destination[0] + i, dj = patch.destination[1] + j;
              const int dk = patch.destination[2] + k;
              const Real value = recv_values[position++];
              if (patch.field == 0) staged_cons(patch.block,v,dk,dj,di) = value;
              if (patch.field == 1) staged_b.x1f(patch.block,dk,dj,di) = value;
              if (patch.field == 2) staged_b.x2f(patch.block,dk,dj,di) = value;
              if (patch.field == 3) staged_b.x3f(patch.block,dk,dj,di) = value;
            }
          }
        }
      }
    }
  }
  if (source.refine) {
    DvceArray5D<Real> device_cons("restart-transfer-cons", nlocal, nvar, hz, hy, hx);
    DvceFaceFld4D<Real> device_b("restart-transfer-b", nlocal, hz, hy, hx);
    Kokkos::deep_copy(device_cons, staged_cons);
    Kokkos::deep_copy(device_b.x1f, staged_b.x1f);
    Kokkos::deep_copy(device_b.x2f, staged_b.x2f);
    Kokkos::deep_copy(device_b.x3f, staged_b.x3f);
    ProlongRestartMHD(pm, pin, device_cons, device_b);
  } else {
    Kokkos::deep_copy(pack->pmhd->u0, staged_cons);
    Kokkos::deep_copy(pack->pmhd->b0.x1f, staged_b.x1f);
    Kokkos::deep_copy(pack->pmhd->b0.x2f, staged_b.x2f);
    Kokkos::deep_copy(pack->pmhd->b0.x3f, staged_b.x3f);
  }
}
