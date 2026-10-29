//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_prolongation.cpp
//! \brief Use native mesh-transfer operators for one uniform restart refinement.

#include <cstdint>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/restart_prolongation.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

void ProlongRestartMHD(Mesh *pm, ParameterInput *pin,
                       DvceArray5D<Real> &coarse_cons, DvceFaceFld4D<Real> &coarse_b) {
  auto *pmhd = pm->pmb_pack->pmhd;
  auto old_cons = pmhd->coarse_u0;
  auto old_prim = pmhd->coarse_w0;
  auto old_b = pmhd->coarse_b0;
  pmhd->coarse_u0 = coarse_cons;
  pmhd->coarse_b0 = coarse_b;

  // Only the transfer uses this adapter; the evolved mesh remains single-level.
  MeshRefinement transfer(pm, pin);
  transfer.new_nmb_eachrank = pm->nmb_eachrank;
  transfer.new_gids_eachrank = pm->gids_eachrank;
  DualArray1D<int> identity("restart_refine_identity", pm->nmb_total);
  for (int m = 0; m < pm->nmb_total; ++m) identity.h_view(m) = m;
  identity.template modify<HostMemSpace>();
  identity.template sync<DevExeSpace>();
  Kokkos::deep_copy(transfer.refine_flag.d_view, 1);

  transfer.RefineCC(identity, pmhd->u0, coarse_cons);
  transfer.RefineFC(identity, pmhd->b0, coarse_b);
  Kokkos::fence();
  const auto eos = pmhd->peos->eos_data;
  if (eos.is_cgl && !eos.passive) {
    Kokkos::realloc(pmhd->coarse_w0, coarse_cons.extent_int(0),
                   coarse_cons.extent_int(1), coarse_cons.extent_int(2),
                   coarse_cons.extent_int(3), coarse_cons.extent_int(4));
    transfer.RefineCGLMHDPrimitives(identity, pmhd);
    auto counters = transfer.ExportCGLAMRRepairCounters();
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, counters.data(), static_cast<int>(counters.size()),
                  MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "Restart refinement CGL projection events: cells="
                << counters[MeshRefinement::cgl_amr_cells_repaired_index]
                << " energy=" << counters[MeshRefinement::cgl_amr_energy_repairs_index]
                << " anisotropy="
                << counters[MeshRefinement::cgl_amr_anisotropy_repairs_index]
                << " slope=" << counters[MeshRefinement::cgl_amr_slope_repairs_index]
                << std::endl;
    }
  }
  // Passive J/A and the independent gas energy retain conservative CC transfer.
  Kokkos::fence();
  pmhd->coarse_u0 = old_cons;
  pmhd->coarse_w0 = old_prim;
  pmhd->coarse_b0 = old_b;
#if MPI_PARALLEL_ENABLED
  MPI_Comm_free(&transfer.amr_comm);
#endif
}
