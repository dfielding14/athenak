//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cgl_lf_stiffness.cpp
//! \brief Rank-local LF timestep maps for selected replay sweeps.

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "diffusion/cgl_landau_fluid.hpp"

void CGLLandauFluid::WriteStiffnessMap(const char *phase, Real dt_sweep,
                                      Real dt_sts, int nstages) {
  auto *pm = pmy_pack->pmesh;
  auto *pmhd = pmy_pack->pmhd;
  const auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int cells_per_block = nx1*nx2*nx3;
  const int ndim = pm->three_d ? 3 : (pm->multi_d ? 2 : 1);
  DvceArray1D<Real> cell_dt("cgl_lf_stiffness_dt", nmb*cells_per_block);
  DvceArray1D<int> worst("cgl_lf_stiffness_worst", nmb);
  DvceArray3D<Real> faces("cgl_lf_stiffness_faces", nmb, 2*ndim+1, face_columns);
  const Real raw_dt = CalculateTimeStep<TimestepSample::cells>(
      pmhd->w0, pmhd->bcc0, pmhd->b0, pmhd->peos->eos_data, cell_dt, worst, faces);
  const auto dt_host = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cell_dt);
  auto worst_host = Kokkos::create_mirror_view(worst);
  for (int m = 0; m < nmb; ++m) {
    int index = 0;
    for (int n = 1; n < cells_per_block; ++n) {
      if (dt_host(m*cells_per_block+n) < dt_host(m*cells_per_block+index)) index = n;
    }
    worst_host(m) = index;
  }
  Kokkos::deep_copy(worst, worst_host);
  CalculateTimeStep<TimestepSample::faces>(
      pmhd->w0, pmhd->bcc0, pmhd->b0, pmhd->peos->eos_data, cell_dt, worst, faces);
  const auto faces_host = Kokkos::create_mirror_view_and_copy(HostMemSpace(), faces);

  const std::uint16_t endian_test = 1;
  const char *endian = *reinterpret_cast<const unsigned char*>(&endian_test) ?
      "little" : "big";
  std::ostringstream metadata;
  metadata << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "cycle=" << pm->ncycle << " phase=" << phase
           << " rank=" << global_variable::my_rank
           << " nranks=" << global_variable::nranks
           << " nblocks=" << nmb << " total_blocks=" << pm->nmb_total
           << " nx1=" << nx1 << " nx2=" << nx2 << " nx3=" << nx3
           << " real_bytes=" << sizeof(Real) << " endian=" << endian
           << " root_level=" << pm->root_level
           << " root_blocks1=" << pm->nmb_rootx1 << " root_blocks2=" << pm->nmb_rootx2
           << " root_blocks3=" << pm->nmb_rootx3
           << " x1min=" << pm->mesh_size.x1min << " x1max=" << pm->mesh_size.x1max
           << " x2min=" << pm->mesh_size.x2min << " x2max=" << pm->mesh_size.x2max
           << " x3min=" << pm->mesh_size.x3min << " x3max=" << pm->mesh_size.x3max
           << " time=" << pm->time << " dt_outer=" << pm->dt
           << " dt_sweep=" << dt_sweep << " sts_safety=" << pm->sts_safety
           << " dt_sts=" << dt_sts << " nstages=" << nstages
           << " raw_dt_rank=" << raw_dt << '\n';
  const std::string stem = "lf-stiffness." + std::to_string(pm->ncycle) + "." +
      phase + ".rank" + std::to_string(global_variable::my_rank);
  const auto fail = [&](const std::string &path) {
    std::cerr << "### FATAL ERROR: cannot write LF stiffness map " << path << std::endl;
    std::exit(EXIT_FAILURE);
  };
  std::ofstream binary(stem + ".bin", std::ios::binary);
  binary << "# AthenaK LF stiffness v1\n" << metadata.str();
  for (int m = 0; m < nmb; ++m) {
    const int gid = pmy_pack->pmb->mb_gid.h_view(m);
    const auto &loc = pm->lloc_eachmb[gid];
    binary << "block " << gid << ' ' << loc.level << ' ' << loc.lx1 << ' '
           << loc.lx2 << ' ' << loc.lx3 << '\n';
  }
  binary << "data\n";
  binary.write(reinterpret_cast<const char*>(dt_host.data()),
               static_cast<std::streamsize>(dt_host.size()*sizeof(Real)));
  binary.close();
  if (!binary) fail(stem + ".bin");

  std::ofstream table(stem + ".tsv");
  table << "# Worst owned cell in each block, with all incident faces.\n"
        << "# " << metadata.str()
        << "rank\tgid\ti\tj\tk\tdir\tside\traw_dt\trho\tppar\tpperp\tbx\tby\tbz"
        << "\tvalid\trho_face\tppar_face\tpperp_face\tbbar\tbn_ct\tgrad_b_parallel"
        << "\tnu_actual\tqpar_ratio\tqperp_ratio\tTrow_parallel\tBrow_parallel"
        << "\tTrow_perp\tBrow_perp\trow_parallel\trow_perp\n"
        << std::setprecision(std::numeric_limits<Real>::max_digits10);
  for (int m = 0; m < nmb; ++m) {
    const int index = worst_host(m);
    const int i = index%nx1, j = (index/nx1)%nx2, k = index/(nx1*nx2);
    for (int side = 0; side < 2*ndim; ++side) {
      table << global_variable::my_rank << '\t' << pmy_pack->pmb->mb_gid.h_view(m)
            << '\t' << i << '\t' << j << '\t' << k << '\t' << side/2
            << '\t' << (side%2 ? 1 : -1) << '\t' << dt_host(m*cells_per_block+index);
      for (int n = 0; n < 6; ++n) table << '\t' << faces_host(m,2*ndim,n);
      for (int n = 0; n < face_columns; ++n) table << '\t' << faces_host(m,side,n);
      table << '\n';
    }
  }
  table.close();
  if (!table) fail(stem + ".tsv");
}
