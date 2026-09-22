//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coarsened_binary.cpp
//! \brief writes output data in binary format, which simply consists of each MeshBlock
//! written contiguously in order of "gid" in binary format.

#include <sys/stat.h>  // mkdir

#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include <utility>
#include <algorithm>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "athena.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"
#include "sgs_moments.hpp"

namespace {

int ActiveCoarsenFactor(int extent, int coarsen_factor) {
  return (extent > 1) ? coarsen_factor : 1;
}

[[noreturn]] void FatalCoarsenedBinaryError(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

// Form moments while filtering: no full-resolution SGS arrays or global atomic adds.
template <bool is_mhd, int powers = 1>
void CoarsenSGS(Mesh *pm, int factor, HostArray5D<Real> &output) {
  auto &indcs = pm->mb_indcs;
  int nx = output.extent_int(4), ny = output.extent_int(3);
  int nz = output.extent_int(2), nvars = output.extent_int(0);
  int nmb = pm->pmb_pack->nmb_thispack;
  int fx = ActiveCoarsenFactor(indcs.nx1, factor);
  int fy = ActiveCoarsenFactor(indcs.nx2, factor);
  int fz = ActiveCoarsenFactor(indcs.nx3, factor);
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int samples = fx * fy * fz;
  // Keep live sums small and expose enough teams even for whole-block filters.
  constexpr int group_size = 8, chunk_cells = 1024;
  int groups = (nvars + group_size - 1) / group_size;
  int chunks = (samples + chunk_cells - 1) / chunk_cells;
  int cells = nmb * nz * ny * nx;
  DvceArray5D<Real> u, bcc;
  if constexpr (is_mhd) {
    u = pm->pmb_pack->pmhd->u0;
    bcc = pm->pmb_pack->pmhd->bcc0;
  } else {
    u = pm->pmb_pack->phydro->u0;
  }
  DvceArray5D<Real> coarse(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                          "sgs_coarse"), nvars, nmb, nz, ny, nx);
  DvceArray3D<Real> partial;
  if (chunks > 1) {
    partial = DvceArray3D<Real>(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               "sgs_partial"), nvars, cells, chunks);
  }
  Kokkos::parallel_for("sgs_3d_moments",
    Kokkos::TeamPolicy<DevExeSpace>(cells * groups * chunks, Kokkos::AUTO),
  KOKKOS_LAMBDA(const TeamMember_t &team) {
    // Neighboring teams reuse the same fine-state chunk across moment groups.
    int group = team.league_rank() % groups;
    int chunk = (team.league_rank() / groups) % chunks;
    int cell = team.league_rank() / (groups * chunks);
    int i = cell % nx, j = (cell / nx) % ny;
    int k = (cell / (nx * ny)) % nz, m = cell / (nx * ny * nz);
    int begin = chunk * chunk_cells;
    int end = (begin + chunk_cells < samples) ? begin + chunk_cells : samples;
    auto reduce_group = [=](auto group_start) {
      int first;
      if constexpr (powers == 1) { first = decltype(group_start)::value; }
      else { first = group_start; }
      Real r0, r1, r2, r3, r4, r5, r6, r7;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, begin, end),
      [=](const int offset, Real &s0, Real &s1, Real &s2, Real &s3,
          Real &s4, Real &s5, Real &s6, Real &s7) {
        int ii = is + i * fx + offset % fx;
        int jj = js + j * fy + (offset / fx) % fy;
        int kk = ks + k * fz + offset / (fx * fy);
        MHDCons1D state{};
        state.d = u(m,IDN,kk,jj,ii);
        state.mx = u(m,IM1,kk,jj,ii);
        state.my = u(m,IM2,kk,jj,ii);
        state.mz = u(m,IM3,kk,jj,ii);
        if constexpr (is_mhd) {
          state.e = u(m,IEN,kk,jj,ii);
          state.bx = bcc(m,IBX,kk,jj,ii);
          state.by = bcc(m,IBY,kk,jj,ii);
          state.bz = bcc(m,IBZ,kk,jj,ii);
        }
        auto moment = [=](int n) {
          if constexpr (is_mhd) {
            Real value = MHDSGSMoment(n / powers, state);
            Real result = value;
            for (int p=0; p<n % powers; ++p) { result *= value; }
            return result;
          } else {
            return HydroSGS3DMoment(n, state);
          }
        };
        s0 += moment(first);     s1 += moment(first + 1);
        s2 += moment(first + 2); s3 += moment(first + 3);
        s4 += moment(first + 4); s5 += moment(first + 5);
        s6 += moment(first + 6); s7 += moment(first + 7);
      }, r0, r1, r2, r3, r4, r5, r6, r7);
      Kokkos::single(Kokkos::PerTeam(team), [=]() {
        auto save = [=](int n, Real value) {
          if (n < nvars) {
            if (chunks == 1) { coarse(n,m,k,j,i) = value/samples; }
            else { partial(n,cell,chunk) = value; }
          }
        };
        save(first, r0);     save(first + 1, r1);
        save(first + 2, r2); save(first + 3, r3);
        save(first + 4, r4); save(first + 5, r5);
        save(first + 6, r6); save(first + 7, r7);
      });
    };
    // Resolve moment selectors once per team, so the fine-cell loop has fixed algebra.
    if constexpr (powers == 1) {
      switch (group) {
        case 0: reduce_group(std::integral_constant<int,0>{}); break;
        case 1: reduce_group(std::integral_constant<int,8>{}); break;
        default:
          if constexpr (is_mhd) {
            switch (group) {
              case 2: reduce_group(std::integral_constant<int,16>{}); break;
              case 3: reduce_group(std::integral_constant<int,24>{}); break;
              case 4: reduce_group(std::integral_constant<int,32>{}); break;
              case 5: reduce_group(std::integral_constant<int,40>{}); break;
              case 6: reduce_group(std::integral_constant<int,48>{}); break;
              case 7: reduce_group(std::integral_constant<int,56>{}); break;
            }
          }
      }
    } else {
      reduce_group(group * group_size);
    }
  });
  if (chunks > 1) {
    Kokkos::parallel_for("sgs_3d_combine",
      Kokkos::TeamPolicy<DevExeSpace>(nvars * cells, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamMember_t &team) {
      int cell = team.league_rank() % cells, n = team.league_rank() / cells;
      int i = cell % nx, j = (cell / nx) % ny;
      int k = (cell / (nx * ny)) % nz, m = cell / (nx * ny * nz);
      Real sum;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, chunks),
      [=](const int chunk, Real &subtotal) { subtotal += partial(n,cell,chunk); }, sum);
      Kokkos::single(Kokkos::PerTeam(team), [=]() {
        coarse(n,m,k,j,i) = sum/samples;
      });
    });
  }
  Kokkos::deep_copy(output, coarse);
}

}  // namespace

//----------------------------------------------------------------------------------------
// Constructor: also calls BaseTypeOutput base class constructor

CoarsenedBinaryOutput::CoarsenedBinaryOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int factor = out_params.coarsen_factor;
  if (factor < 1) {
    FatalCoarsenedBinaryError("coarsen_factor must be positive.");
  }
  if (out_params.include_gzs || out_params.slice1 || out_params.slice2 ||
      out_params.slice3 || out_params.gid >= 0) {
    FatalCoarsenedBinaryError(
        "Coarsened binary output requires the full domain without ghost zones, "
        "slices, or gid selection.");
  }
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  if ((indcs.nx1 > 1 && indcs.nx1 % factor != 0) ||
      (indcs.nx2 > 1 && indcs.nx2 % factor != 0) ||
      (indcs.nx3 > 1 && indcs.nx3 % factor != 0)) {
    FatalCoarsenedBinaryError(
        "Every active MeshBlock dimension must be divisible by coarsen_factor.");
  }
  if ((out_params.variable.compare("hydro_sgs_2d") == 0 ||
       out_params.variable.compare("hydro_sgs_3d") == 0) &&
      out_params.compute_moments) {
    FatalCoarsenedBinaryError(
        "hydro_sgs output does not support compute_moments=true.");
  }
  // create directories for outputs
  // useful for mpiio-based outputs because on some supercomputers you may need to
  // set different stripe counts depending on whether mpiio is used in order to
  // achieve the best performance and not to crash the filesystem
  std::string dir_name;
  dir_name.assign("cbin_");
  dir_name.append(out_params.file_id);
  dir_name.append("_");
  dir_name.append(std::to_string(out_params.coarsen_factor));
  mkdir(dir_name.c_str(),0775);
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    dir_name.append("/");
    dir_name.append(rank_dir);
    mkdir(dir_name.c_str(), 0775);
  }
}

//----------------------------------------------------------------------------------------
// BaseTypeOutput::LoadOutputData()
// create std::vector of HostArray3Ds containing data specified in <output> block for
// this output type

void CoarsenedBinaryOutput::LoadOutputData(Mesh *pm) {
  // out_data_ vector (indexed over # of output MBs) stores 4D array of variables
  // so start iteration over number of MeshBlocks
  // TODO(@user): get this working for multiple physics, which may be either defined/undef

  // With AMR, number and location of output MBs can change between output times.
  // So start with clean vector of output MeshBlock info, and re-compute
  outmbs.clear();

  // Use active-cell bounds; unsupported selections are rejected by the constructor.
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  auto &size  = pm->pmb_pack->pmb->mb_size;
  for (int m=0; m<(pm->pmb_pack->nmb_thispack); ++m) {
    int ois = indcs.is, oie = indcs.ie;
    int ojs = indcs.js, oje = indcs.je;
    int oks = indcs.ks, oke = indcs.ke;

    // set coordinate geometry information for MB
    Real x1min = size.h_view(m).x1min;
    Real x1max = size.h_view(m).x1max;
    Real x2min = size.h_view(m).x2min;
    Real x2max = size.h_view(m).x2max;
    Real x3min = size.h_view(m).x3min;
    Real x3max = size.h_view(m).x3max;

    int id = pm->pmb_pack->pmb->mb_gid.h_view(m);
    outmbs.emplace_back(id,ois,oie,ojs,oje,oks,oke,x1min,x1max,x2min,x2max,x3min,x3max);
  }

  // get number of output vars and MBs, then realloc outarray (HostArray)
  int nout_vars_with_moments;
  if (out_params.compute_moments) {
    nout_vars_with_moments = outvars.size() * 4;
  } else {
    nout_vars_with_moments = outvars.size();
  }
  int nout_vars = outvars.size();
  int nout_mbs = outmbs.size();
  // note that while ois,oie,etc. can be different on each MB, the number of cells output
  // on each MeshBlock, i.e. (ois-ois+1), etc. is the same.
  if (nout_mbs > 0) {
    int full_nout1 = outmbs[0].oie - outmbs[0].ois + 1;
    int full_nout2 = outmbs[0].oje - outmbs[0].ojs + 1;
    int full_nout3 = outmbs[0].oke - outmbs[0].oks + 1;
    int nout1 = full_nout1/ActiveCoarsenFactor(full_nout1, out_params.coarsen_factor);
    int nout2 = full_nout2/ActiveCoarsenFactor(full_nout2, out_params.coarsen_factor);
    int nout3 = full_nout3/ActiveCoarsenFactor(full_nout3, out_params.coarsen_factor);
    // NB: outarray stores all output data on Host
    // Degenerate dimensions remain one cell wide.
    Kokkos::realloc(outarray, nout_vars_with_moments, nout_mbs, nout3, nout2, nout1);
  }

  if (out_params.variable.compare("hydro_sgs_2d") == 0) {
    Kokkos::Profiling::ScopedRegion region("SGS2D/load");
    int nx = outarray.extent_int(4), ny = outarray.extent_int(3);
    int fx = ActiveCoarsenFactor(indcs.nx1, out_params.coarsen_factor);
    int fy = ActiveCoarsenFactor(indcs.nx2, out_params.coarsen_factor);
    int samples = fx * fy;
    // Bound each team's work so even a whole-block filter uses many GPU teams.
    constexpr int chunk_cells = 1024;
    int chunks = (samples + chunk_cells - 1) / chunk_cells;
    int cells = nout_mbs * ny * nx;
    int is = indcs.is, js = indcs.js, ks = indcs.ks;
    auto u = pm->pmb_pack->phydro->u0;
    DvceArray5D<Real> coarse(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                            "sgs_coarse"), 6, nout_mbs, 1, ny, nx);
    auto partial = coarse;
    if (chunks > 1) {
      partial = DvceArray5D<Real>(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                 "sgs_partial"), 6, nout_mbs, chunks, ny, nx);
    }

    Kokkos::parallel_for("sgs_2d_moments",
      Kokkos::TeamPolicy<DevExeSpace>(cells * chunks, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamMember_t &team) {
      int chunk = team.league_rank() % chunks;
      int cell = team.league_rank() / chunks;
      int i = cell % nx, j = (cell / nx) % ny, m = cell / (nx * ny);
      int begin = chunk * chunk_cells;
      int end = (begin + chunk_cells < samples) ? begin + chunk_cells : samples;
      Real r, x, y, xx, xy, yy;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, begin, end),
      [=](const int offset, Real &sr, Real &sx, Real &sy,
          Real &sxx, Real &sxy, Real &syy) {
        int ii = is + i * fx + offset % fx;
        int jj = js + j * fy + offset / fx;
        Real rho = u(m,IDN,ks,jj,ii);
        Real mx = u(m,IM1,ks,jj,ii), my = u(m,IM2,ks,jj,ii);
        sr += rho; sx += mx; sy += my;
        sxx += mx*mx/rho; sxy += mx*my/rho; syy += my*my/rho;
      }, r, x, y, xx, xy, yy);
      Kokkos::single(Kokkos::PerTeam(team), [=]() {
        Real norm = (chunks == 1) ? 1.0/samples : 1.0;
        partial(0,m,chunk,j,i) = norm*r;
        partial(1,m,chunk,j,i) = norm*x;
        partial(2,m,chunk,j,i) = norm*y;
        partial(3,m,chunk,j,i) = norm*xx;
        partial(4,m,chunk,j,i) = norm*xy;
        partial(5,m,chunk,j,i) = norm*yy;
      });
    });

    if (chunks > 1) {
      Kokkos::parallel_for("sgs_2d_combine",
        Kokkos::TeamPolicy<DevExeSpace>(cells, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamMember_t &team) {
        int cell = team.league_rank();
        int i = cell % nx, j = (cell / nx) % ny, m = cell / (nx * ny);
        Real r, x, y, xx, xy, yy;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, chunks),
        [=](const int chunk, Real &sr, Real &sx, Real &sy,
            Real &sxx, Real &sxy, Real &syy) {
          sr += partial(0,m,chunk,j,i);
          sx += partial(1,m,chunk,j,i);
          sy += partial(2,m,chunk,j,i);
          sxx += partial(3,m,chunk,j,i);
          sxy += partial(4,m,chunk,j,i);
          syy += partial(5,m,chunk,j,i);
        }, r, x, y, xx, xy, yy);
        Kokkos::single(Kokkos::PerTeam(team), [=]() {
          coarse(0,m,0,j,i) = r/samples;
          coarse(1,m,0,j,i) = x/samples;
          coarse(2,m,0,j,i) = y/samples;
          coarse(3,m,0,j,i) = xx/samples;
          coarse(4,m,0,j,i) = xy/samples;
          coarse(5,m,0,j,i) = yy/samples;
        });
      });
    }
    // One transfer per filter width, independent of variable and MeshBlock counts.
    Kokkos::deep_copy(outarray, coarse);
    for (int m=0; m<nout_mbs; ++m) {
      for (int j=0; j<ny; ++j) {
        for (int i=0; i<nx; ++i) {
          Real rho = outarray(0,m,0,j,i);
          if (!(rho > 0.0)) {
            FatalCoarsenedBinaryError(
                "hydro_sgs_2d encountered non-positive filtered density.");
          }
          Real mx = outarray(1,m,0,j,i), my = outarray(2,m,0,j,i);
          outarray(1,m,0,j,i) = mx/rho;
          outarray(2,m,0,j,i) = my/rho;
          outarray(3,m,0,j,i) -= mx*mx/rho;
          outarray(4,m,0,j,i) -= mx*my/rho;
          outarray(5,m,0,j,i) -= my*my/rho;
        }
      }
    }
    return;
  }

  if (out_params.variable.compare("mhd_sgs") == 0) {
    Kokkos::Profiling::ScopedRegion region("MHD_SGS/load");
    if (out_params.compute_moments) {
      CoarsenSGS<true,4>(pm, out_params.coarsen_factor, outarray);
    } else {
      CoarsenSGS<true>(pm, out_params.coarsen_factor, outarray);
    }
    return;
  }
  if (out_params.variable.compare("hydro_sgs_3d") == 0) {
    Kokkos::Profiling::ScopedRegion region("SGS3D/load");
    CoarsenSGS<false>(pm, out_params.coarsen_factor, outarray);
    for (int m=0; m<nout_mbs; ++m) {
      for (int k=0; k<outarray.extent_int(2); ++k) {
        for (int j=0; j<outarray.extent_int(3); ++j) {
          for (int i=0; i<outarray.extent_int(4); ++i) {
            Real rho = outarray(0,m,k,j,i);
            if (!(rho > 0.0)) {
              FatalCoarsenedBinaryError(
                  "hydro_sgs_3d encountered non-positive filtered density.");
            }
            Real mx = outarray(1,m,k,j,i);
            Real my = outarray(2,m,k,j,i);
            Real mz = outarray(3,m,k,j,i);
            outarray(1,m,k,j,i) = mx/rho;
            outarray(2,m,k,j,i) = my/rho;
            outarray(3,m,k,j,i) = mz/rho;
            outarray(4,m,k,j,i) -= mx*mx/rho;
            outarray(5,m,k,j,i) -= mx*my/rho;
            outarray(6,m,k,j,i) -= mx*mz/rho;
            outarray(7,m,k,j,i) -= my*my/rho;
            outarray(8,m,k,j,i) -= my*mz/rho;
            outarray(9,m,k,j,i) -= mz*mz/rho;
          }
        }
      }
    }
    return;
  }

  // Calculate derived variables, if required
  if (out_params.contains_derived) {
    ComputeDerivedVariable(out_params.variable, pm);
  }

  // Now copy data to host (outarray) over all variables and MeshBlocks
  for (int n=0; n<nout_vars; ++n) {
    for (int m=0; m<nout_mbs; ++m) {
      int mbi = pm->FindMeshBlockIndex(outmbs[m].mb_gid);
      std::pair<int,int> irange = std::make_pair(outmbs[m].ois, outmbs[m].oie+1);
      std::pair<int,int> jrange = std::make_pair(outmbs[m].ojs, outmbs[m].oje+1);
      std::pair<int,int> krange = std::make_pair(outmbs[m].oks, outmbs[m].oke+1);
      std::pair<int,int> moment_range;
      if (out_params.compute_moments) {
        moment_range = std::make_pair(n*4, n*4+4);
      } else {
        moment_range = std::make_pair(n, n+1);
      }
      int nout1 = (outmbs[0].oie - outmbs[0].ois + 1);
      int nout2 = (outmbs[0].oje - outmbs[0].ojs + 1);
      int nout3 = (outmbs[0].oke - outmbs[0].oks + 1);
      int coarsen1 = ActiveCoarsenFactor(nout1, out_params.coarsen_factor);
      int coarsen2 = ActiveCoarsenFactor(nout2, out_params.coarsen_factor);
      int coarsen3 = ActiveCoarsenFactor(nout3, out_params.coarsen_factor);
      int coarsened_nout1 = nout1/coarsen1;
      int coarsened_nout2 = nout2/coarsen2;
      int coarsened_nout3 = nout3/coarsen3;

      auto d_slice = Kokkos::subview(*(outvars[n].data_ptr), mbi, outvars[n].data_index,
                                     krange,jrange,irange);
      int number_of_moments = 1;
      if (out_params.compute_moments) {
        number_of_moments = 4;
      }
      DvceArray4D<Real> d_output_var_coarsened("d_output_var_coarsened",
        number_of_moments, coarsened_nout3, coarsened_nout2, coarsened_nout1);

      int coarsen_cells = coarsen1 * coarsen2 * coarsen3;

      if (nout1 % coarsen1 != 0 || nout2 % coarsen2 != 0 || nout3 % coarsen3 != 0) {
        FatalCoarsenedBinaryError(
            "Active output dimensions must be divisible by coarsen_factor.");
      }

      // One team sums each coarse cell, avoiding factor^D contended atomic adds.
      int coarse_cells = coarsened_nout3 * coarsened_nout2 * coarsened_nout1;
      Kokkos::parallel_for("coarsen_variable",
        Kokkos::TeamPolicy<DevExeSpace>(number_of_moments * coarse_cells, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamMember_t &team) {
        int idx = team.league_rank();
        int moment_idx = idx / coarse_cells;
        int k_c = (idx / (coarsened_nout2 * coarsened_nout1)) % coarsened_nout3;
        int j_c = (idx / coarsened_nout1) % coarsened_nout2;
        int i_c = idx % coarsened_nout1;
        Real sum;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, coarsen_cells),
        [=](const int offset, Real &subtotal) {
          int k = k_c * coarsen3 + offset / (coarsen2 * coarsen1);
          int j = j_c * coarsen2 + (offset / coarsen1) % coarsen2;
          int i = i_c * coarsen1 + offset % coarsen1;
          Real value = d_slice(k,j,i);
          Real moment = value;
          for (int p=0; p<moment_idx; ++p) { moment *= value; }
          subtotal += moment;
        }, sum);
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          d_output_var_coarsened(moment_idx,k_c,j_c,i_c) = sum/coarsen_cells;
        });
      });


      // Now, create a host mirror for the coarsened data.
      DvceArray4D<Real>::HostMirror h_output_var = Kokkos::create_mirror(
        d_output_var_coarsened
      );

      // Copy the coarsened data to the host mirror.
      Kokkos::deep_copy(h_output_var, d_output_var_coarsened);
      Kokkos::fence(); // Ensure complete copy before using h_output_var on the host

      // copy host mirror to 5D host View containing all output variables
      // if (out_params.compute_moments) {
      auto h_slice = Kokkos::subview(outarray,
        moment_range,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL
      );
      Kokkos::deep_copy(h_slice,h_output_var);
    }
  }

  // Each filter width owns an output object; do not retain fine-grid scratch per width.
  if (out_params.contains_derived) {
    derived_var = DvceArray5D<Real>();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CoarsenedBinaryOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in Coarsenedbinary format
//   All MeshBlocks are written to the same file.

void CoarsenedBinaryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  Kokkos::Profiling::ScopedRegion region("cbin/write");
  // create filename: "cbin_"+"file_id"+"_"+"coarsening_factor"+"/file_basename"
  // + "." + "file_id" + "." + XXXXX + ".cbin"
  // where XXXXX = 5-digit file_number
  bool single_file_per_rank = out_params.single_file_per_rank;

  std::string fname;
  char number[6];
  std::snprintf(number, sizeof(number), "%05d", out_params.file_number);

  fname.assign("cbin_");
  fname.append(out_params.file_id);
  fname.append("_");
  fname.append(std::to_string(out_params.coarsen_factor));
  fname.append("/");
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname.append(rank_dir);
  }
  fname.append(out_params.file_basename);
  fname.append(".");
  fname.append(out_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".cbin");

  IOWrapper cbinfile;
  std::size_t header_offset=0;
  cbinfile.Open(fname.c_str(), IOWrapper::FileMode::write, single_file_per_rank);

  int number_of_moments = 1;
  if (out_params.compute_moments) {
    number_of_moments = 4;
  }

  // Basic parts of the format:
  // 1. Size of the header
  // 2. Current time
  // 3. List of variables in the file
  // 4. Header (input file information)
  {std::stringstream msg;
  msg << "Athena binary output version=1.1" << std::endl
      // preheader size includes "size of preheader" line up to "number of variables"
      << "  size of preheader=7" << std::endl
      << "  time=" << pm->time << std::endl
      << "  cycle=" << pm->ncycle << std::endl
      << "  number of moments=" << number_of_moments << std::endl
      << "  coarsening factor=" << out_params.coarsen_factor << std::endl
      << "  size of location=" << sizeof(Real) << std::endl
      << "  size of variable=" << sizeof(float) << std::endl
      << "  number of variables=" << outvars.size()*number_of_moments << std::endl
      << "  variables:  ";
  if (out_params.compute_moments) {
    // need to write the label for each of the 4 moments
    for (int n=0; n<outvars.size(); n++) {
      msg << outvars[n].label.c_str() << "_1st  ";
      msg << outvars[n].label.c_str() << "_2nd  ";
      msg << outvars[n].label.c_str() << "_3rd  ";
      msg << outvars[n].label.c_str() << "_4th  ";
    }
  } else {
    for (int n=0; n<outvars.size(); n++) {
      msg << outvars[n].label.c_str() << "  ";
    }
  }
  msg << std::endl;
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    cbinfile.Write_any_type(msg.str().c_str(),msg.str().size(), "byte",
                            single_file_per_rank);
  }
  header_offset += msg.str().size();}
  {std::stringstream msg;
  // prepare the input parameters
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf=ost.str();
  msg << "  header offset=" << sbuf.size()*sizeof(char)  << std::endl;
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    cbinfile.Write_any_type(msg.str().c_str(),msg.str().size(), "byte",
                            single_file_per_rank);
    cbinfile.Write_any_type(sbuf.c_str(),sbuf.size(), "byte", single_file_per_rank);
  }
  header_offset += sbuf.size()*sizeof(char);
  header_offset += msg.str().size();}

  //  5. Data.  An arbitrary number of scalars and vectors can be written (every node
  //  in the OutputData doubly linked lists), all in binary floats format

  int nout_vars = outvars.size();
  if (out_params.compute_moments) {
    nout_vars *= 4;
  }
  int nout_mbs = outmbs.size();
  int full_nout1 = outmbs[0].oie - outmbs[0].ois + 1;
  int full_nout2 = outmbs[0].oje - outmbs[0].ojs + 1;
  int full_nout3 = outmbs[0].oke - outmbs[0].oks + 1;
  int nout1 = full_nout1/ActiveCoarsenFactor(full_nout1, out_params.coarsen_factor);
  int nout2 = full_nout2/ActiveCoarsenFactor(full_nout2, out_params.coarsen_factor);
  int nout3 = full_nout3/ActiveCoarsenFactor(full_nout3, out_params.coarsen_factor);
  int cells = nout1*nout2*nout3;


  // ois, oie, ojs, oje, oks, oke + il1, il2, il3, level +
  // x1min, x1max, x2min, x2max, x3min, x3max + data
  std::size_t data_size = 10*sizeof(int32_t) + 6*sizeof(Real)
                        + (cells*nout_vars)*sizeof(float);

  int ns_mbs = pm->gids_eachrank[global_variable::my_rank];
  int nb_mbs = pm->nmb_eachrank[global_variable::my_rank];

  // allocate 1D vector of floats used to convert and output data
  char *data = new char[nb_mbs*data_size];
  float *single_data = new float[cells];

  // Loop over MeshBlocks
  for (int m=0; m<nout_mbs; ++m) {
    char *pdata=&(data[m*data_size]);
    LogicalLocation loc = pm->lloc_eachmb[outmbs[m].mb_gid];
    // of the starting indexes maybe I need to subtract of nghost,
    // divide by coarsen factor, and then add nghost back in
    int ois = outmbs[m].ois;
    int oie = outmbs[m].ois+nout1-1;
    int ojs = outmbs[m].ojs;
    int oje = outmbs[m].ojs+nout2-1;
    int oks = outmbs[m].oks;
    int oke = outmbs[m].oks+nout3-1;

    // output indexing for MB
    int32_t nx = (int32_t)(ois);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(oie);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(ojs);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(oje);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(oks);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(oke);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);


    // TODO(@DBF): not sure how to shift these properly for the reduced grid
    // logical location lx1, lx2, lx3
    nx = (int32_t)(loc.lx1);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(loc.lx2);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(loc.lx3);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);

    // TODO(@DBF): This probably won't work for AMR
    // physical refinement level
    nx = (int32_t)(loc.level-pm->root_level);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);

    // coordinate location
    Real xv = outmbs[m].x1min;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);
    xv = outmbs[m].x1max;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);
    xv = outmbs[m].x2min;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);
    xv = outmbs[m].x2max;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);
    xv = outmbs[m].x3min;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);
    xv = outmbs[m].x3max;
    memcpy(pdata,&(xv),sizeof(xv));
    pdata+=sizeof(xv);

    // output variables
    float tmp_data;
    for (int n=0; n<nout_vars; n++) {
      int cnt=0;
      for (int k=oks; k<=oke; k++) {
        for (int j=ojs; j<=oje; j++) {
          for (int i=ois; i<=oie; i++) {
            tmp_data = static_cast<float>(outarray(n,m,k-oks,j-ojs,i-ois));
            single_data[cnt] = tmp_data;
            cnt++;
          }
        }
      }
      memcpy(pdata,single_data,cells*sizeof(float));
      pdata+=cells*sizeof(float);
    }
  }

  // now write Coarsenedbinary data
  // check if elements larger than 2^31
  if (data_size*nb_mbs<=2147483648) {
    // now write Coarsenedbinary data in parallel
    std::size_t myoffset=header_offset;
    if (!single_file_per_rank) {
      myoffset += data_size*ns_mbs;
    }
    if (cbinfile.Write_any_type_at_all(data,(data_size*nb_mbs),myoffset,"byte",
                                     single_file_per_rank) != data_size*nb_mbs) {
      FatalCoarsenedBinaryError("Coarsened binary data were not written completely.");
    }
  } else {
    // write data over each MeshBlock sequentially and in parallel
    // calculate max/min number of MeshBlocks across all ranks
    noutmbs_max = pm->nmb_eachrank[0];
    noutmbs_min = pm->nmb_eachrank[0];
    for (int i=0; i<(global_variable::nranks); ++i) {
      noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
      noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
    }
    for (int m=0;  m<noutmbs_max; ++m) {
      char *pdata=&(data[m*data_size]);
      std::size_t myoffset=header_offset + data_size*m;
      if (!single_file_per_rank) {
        myoffset += data_size*ns_mbs;
      }
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        if (cbinfile.Write_any_type_at_all(pdata,(data_size),myoffset,"byte",
                                            single_file_per_rank) != data_size) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "binary data not written correctly to binary file, "
              << "binary file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        if (cbinfile.Write_any_type_at(pdata,(data_size),myoffset,"byte",
                                        single_file_per_rank) != data_size) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
               << std::endl << "binary data not written correctly to binary file, "
               << "binary file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // close the output file and clean up ptrs to data
  cbinfile.Close(single_file_per_rank);
  delete [] data;
  delete [] single_data;

  // increment counters
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}
