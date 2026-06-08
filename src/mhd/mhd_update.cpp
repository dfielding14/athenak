//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_update.cpp
//! \brief Performs explicit update of MHD conserved variables (u0) for each stage of the
//! SSP RK integrators (e.g. RK1, RK2, RK3) implemented in AthenaK, using weighted average
//! and partial time update of flux divergence. Source terms are added in the
//! MHDSrcTerms() function.

#include <cstdlib>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \brief Fail with a global count and first-cell report for a nonfinite CGL RK state.

void MHD::DiagnoseNonfiniteRKState(int stage, const char *phase,
                                  DvceArray5D<Real> state) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nmkji = pmy_pack->nmb_thispack * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;
  auto state_ = state;

  int bad_density = 0;
  int bad_momentum1 = 0;
  int bad_momentum2 = 0;
  int bad_momentum3 = 0;
  int bad_energy = 0;
  int bad_anisotropy = 0;
  Kokkos::parallel_reduce(
      "cgl_nonfinite_rk_state_audit",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, int& density_count,
                    int& momentum1_count, int& momentum2_count,
                    int& momentum3_count, int& energy_count,
                    int& anisotropy_count) {
        const int m = idx / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        const int i = idx - m * nkji - k * nji - j * nx1 + is;
        k += ks;
        j += js;
        if (!Kokkos::isfinite(state_(m, IDN, k, j, i))) ++density_count;
        if (!Kokkos::isfinite(state_(m, IM1, k, j, i))) ++momentum1_count;
        if (!Kokkos::isfinite(state_(m, IM2, k, j, i))) ++momentum2_count;
        if (!Kokkos::isfinite(state_(m, IM3, k, j, i))) ++momentum3_count;
        if (!Kokkos::isfinite(state_(m, IEN, k, j, i))) ++energy_count;
        if (!Kokkos::isfinite(state_(m, IAN, k, j, i))) ++anisotropy_count;
      },
      Kokkos::Sum<int>(bad_density),
      Kokkos::Sum<int>(bad_momentum1),
      Kokkos::Sum<int>(bad_momentum2),
      Kokkos::Sum<int>(bad_momentum3),
      Kokkos::Sum<int>(bad_energy),
      Kokkos::Sum<int>(bad_anisotropy));

  int bad_local[6] = {
      bad_density, bad_momentum1, bad_momentum2,
      bad_momentum3, bad_energy, bad_anisotropy};
  int bad_global[6] = {
      bad_density, bad_momentum1, bad_momentum2,
      bad_momentum3, bad_energy, bad_anisotropy};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(bad_local, bad_global, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  const int total_bad = bad_global[0] + bad_global[1] + bad_global[2] +
                        bad_global[3] + bad_global[4] + bad_global[5];
  if (total_bad == 0) return;

  int first_bad = std::numeric_limits<int>::max();
  Kokkos::parallel_reduce(
      "cgl_nonfinite_rk_state_first_cell",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int idx, int& minimum) {
        const int m = idx / nkji;
        int k = (idx - m * nkji) / nji;
        int j = (idx - m * nkji - k * nji) / nx1;
        const int i = idx - m * nkji - k * nji - j * nx1 + is;
        k += ks;
        j += js;
        const bool bad =
            !Kokkos::isfinite(state_(m, IDN, k, j, i)) ||
            !Kokkos::isfinite(state_(m, IM1, k, j, i)) ||
            !Kokkos::isfinite(state_(m, IM2, k, j, i)) ||
            !Kokkos::isfinite(state_(m, IM3, k, j, i)) ||
            !Kokkos::isfinite(state_(m, IEN, k, j, i)) ||
            !Kokkos::isfinite(state_(m, IAN, k, j, i));
        if (bad && idx < minimum) minimum = idx;
      },
      Kokkos::Min<int>(first_bad));

  int first_rank =
      (first_bad < nmkji) ? global_variable::my_rank : global_variable::nranks;
#if MPI_PARALLEL_ENABLED
  int global_first_rank = global_variable::nranks;
  MPI_Allreduce(&first_rank, &global_first_rank, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  first_rank = global_first_rank;
#endif
  if (global_variable::my_rank == first_rank && first_bad < nmkji) {
    const int m = first_bad / nkji;
    int k = (first_bad - m * nkji) / nji;
    int j = (first_bad - m * nkji - k * nji) / nx1;
    const int i = first_bad - m * nkji - k * nji - j * nx1 + is;
    k += ks;
    j += js;
    auto cell = Kokkos::subview(state_, m, Kokkos::ALL, k, j, i);
    auto host_cell = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cell);
    std::cout.precision(std::numeric_limits<Real>::max_digits10);
    std::cout << "CGL RK state first nonfinite cell: phase=" << phase
              << " rank=" << first_rank
              << " m=" << m << " k=" << k << " j=" << j << " i=" << i
              << " state={density:" << host_cell(IDN)
              << ",momentum1:" << host_cell(IM1)
              << ",momentum2:" << host_cell(IM2)
              << ",momentum3:" << host_cell(IM3)
              << ",energy:" << host_cell(IEN)
              << ",anisotropy:" << host_cell(IAN) << "}" << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    std::cout << "### FATAL ERROR in CGL RK state diagnostic: phase=" << phase
              << " time=" << pmy_pack->pmesh->time
              << " cycle=" << pmy_pack->pmesh->ncycle
              << " stage=" << stage
              << " nonfinite_counts={density:" << bad_global[0]
              << ",momentum1:" << bad_global[1]
              << ",momentum2:" << bad_global[2]
              << ",momentum3:" << bad_global[3]
              << ",energy:" << bad_global[4]
              << ",anisotropy:" << bad_global[5] << "}" << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------
//! \fn  void MHD::Update
//  \brief Explicit RK update including flux divergence terms

TaskStatus MHD::RKUpdate(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  Real &gam0 = pdriver->gam0[stage-1];
  Real &gam1 = pdriver->gam1[stage-1];
  Real beta_dt = (pdriver->beta[stage-1])*(pmy_pack->pmesh->dt);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  int nv1 = nmhd + nscalars - 1;
  auto u0_ = u0;
  auto u1_ = u1;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto &mbsize = pmy_pack->pmb->mb_size;

  if (diagnose_nonfinite_rk_update) {
    DiagnoseNonfiniteRKState(stage, "pre-u0", u0_);
    DiagnoseNonfiniteRKState(stage, "pre-u1", u1_);
  }

  // hierarchical parallel loop that updates conserved variables to intermediate step
  // using weights and fractional time step appropriate to stages of time-integrator used
  // Vector inner loop used for good performance on cpus
  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1);

  par_for_outer("mhd_update",DevExeSpace(),scr_size,scr_level,0,nmb1,0,nv1,ks,ke,js,je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int n, const int k, const int j) {
    ScrArray1D<Real> divf(member.team_scratch(scr_level), ncells1);

    // compute dF1/dx1
    par_for_inner(member, is, ie, [&](const int i) {
      divf(i) = (flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i))/mbsize.d_view(m).dx1;
    });
    member.team_barrier();

    // Add dF2/dx2
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (multi_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += (flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i))/mbsize.d_view(m).dx2;
      });
      member.team_barrier();
    }

    // Add dF3/dx3
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (three_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += (flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i))/mbsize.d_view(m).dx3;
      });
      member.team_barrier();
    }

    par_for_inner(member, is, ie, [&](const int i) {
      u0_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - beta_dt*divf(i);
    });
  });

  if (diagnose_nonfinite_rk_update) {
    const int nx1 = indcs.nx1;
    const int nx2 = indcs.nx2;
    const int nx3 = indcs.nx3;
    const int nmkji = (nmb1 + 1) * nx3 * nx2 * nx1;
    const int nkji = nx3 * nx2 * nx1;
    const int nji = nx2 * nx1;
    int bad_density = 0;
    int bad_momentum1 = 0;
    int bad_momentum2 = 0;
    int bad_momentum3 = 0;
    int bad_energy = 0;
    int bad_anisotropy = 0;
    Kokkos::parallel_reduce(
        "cgl_nonfinite_rk_update_audit",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
        KOKKOS_LAMBDA(const int idx, int& density_count,
                      int& momentum1_count, int& momentum2_count,
                      int& momentum3_count, int& energy_count,
                      int& anisotropy_count) {
          const int m = idx / nkji;
          int k = (idx - m * nkji) / nji;
          int j = (idx - m * nkji - k * nji) / nx1;
          const int i = idx - m * nkji - k * nji - j * nx1 + is;
          k += ks;
          j += js;
          if (!Kokkos::isfinite(u0_(m, IDN, k, j, i))) ++density_count;
          if (!Kokkos::isfinite(u0_(m, IM1, k, j, i))) ++momentum1_count;
          if (!Kokkos::isfinite(u0_(m, IM2, k, j, i))) ++momentum2_count;
          if (!Kokkos::isfinite(u0_(m, IM3, k, j, i))) ++momentum3_count;
          if (!Kokkos::isfinite(u0_(m, IEN, k, j, i))) ++energy_count;
          if (!Kokkos::isfinite(u0_(m, IAN, k, j, i))) ++anisotropy_count;
        },
        Kokkos::Sum<int>(bad_density),
        Kokkos::Sum<int>(bad_momentum1),
        Kokkos::Sum<int>(bad_momentum2),
        Kokkos::Sum<int>(bad_momentum3),
        Kokkos::Sum<int>(bad_energy),
        Kokkos::Sum<int>(bad_anisotropy));

    int bad_local[6] = {
        bad_density, bad_momentum1, bad_momentum2,
        bad_momentum3, bad_energy, bad_anisotropy};
    int bad_global[6] = {
        bad_density, bad_momentum1, bad_momentum2,
        bad_momentum3, bad_energy, bad_anisotropy};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(bad_local, bad_global, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    const int total_bad = bad_global[0] + bad_global[1] + bad_global[2] +
                          bad_global[3] + bad_global[4] + bad_global[5];
    if (total_bad > 0) {
      int first_bad = std::numeric_limits<int>::max();
      Kokkos::parallel_reduce(
          "cgl_nonfinite_rk_update_first_cell",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
          KOKKOS_LAMBDA(const int idx, int& minimum) {
            const int m = idx / nkji;
            int k = (idx - m * nkji) / nji;
            int j = (idx - m * nkji - k * nji) / nx1;
            const int i = idx - m * nkji - k * nji - j * nx1 + is;
            k += ks;
            j += js;
            const bool bad =
                !Kokkos::isfinite(u0_(m, IDN, k, j, i)) ||
                !Kokkos::isfinite(u0_(m, IM1, k, j, i)) ||
                !Kokkos::isfinite(u0_(m, IM2, k, j, i)) ||
                !Kokkos::isfinite(u0_(m, IM3, k, j, i)) ||
                !Kokkos::isfinite(u0_(m, IEN, k, j, i)) ||
                !Kokkos::isfinite(u0_(m, IAN, k, j, i));
            if (bad && idx < minimum) minimum = idx;
          },
          Kokkos::Min<int>(first_bad));

      int first_rank =
          (first_bad < nmkji) ? global_variable::my_rank : global_variable::nranks;
#if MPI_PARALLEL_ENABLED
      int global_first_rank = global_variable::nranks;
      MPI_Allreduce(&first_rank, &global_first_rank, 1, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD);
      first_rank = global_first_rank;
#endif
      if (global_variable::my_rank == first_rank && first_bad < nmkji) {
        const int m = first_bad / nkji;
        int k = (first_bad - m * nkji) / nji;
        int j = (first_bad - m * nkji - k * nji) / nx1;
        const int i = first_bad - m * nkji - k * nji - j * nx1 + is;
        k += ks;
        j += js;
        auto cell = Kokkos::subview(u0_, m, Kokkos::ALL, k, j, i);
        auto host_cell =
            Kokkos::create_mirror_view_and_copy(HostMemSpace(), cell);
        std::cout.precision(std::numeric_limits<Real>::max_digits10);
        std::cout << "CGL RK update first nonfinite cell: phase=post-u0 rank="
                  << first_rank << " m=" << m << " k=" << k
                  << " j=" << j << " i=" << i
                  << " state={density:" << host_cell(IDN)
                  << ",momentum1:" << host_cell(IM1)
                  << ",momentum2:" << host_cell(IM2)
                  << ",momentum3:" << host_cell(IM3)
                  << ",energy:" << host_cell(IEN)
                  << ",anisotropy:" << host_cell(IAN) << "}" << std::endl;
      }
#if MPI_PARALLEL_ENABLED
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      if (global_variable::my_rank == 0) {
        std::cout << "### FATAL ERROR in CGL RK update diagnostic: "
                  << "phase=post-u0 time=" << pmy_pack->pmesh->time
                  << " cycle=" << pmy_pack->pmesh->ncycle
                  << " stage=" << stage
                  << " nonfinite_counts={density:" << bad_global[0]
                  << ",momentum1:" << bad_global[1]
                  << ",momentum2:" << bad_global[2]
                  << ",momentum3:" << bad_global[3]
                  << ",energy:" << bad_global[4]
                  << ",anisotropy:" << bad_global[5] << "}" << std::endl;
      }
#if MPI_PARALLEL_ENABLED
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      std::exit(EXIT_FAILURE);
    }
  }
  return TaskStatus::complete;
}
} // namespace mhd
