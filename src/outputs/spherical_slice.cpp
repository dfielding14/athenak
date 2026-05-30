//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_slice.cpp
//! \brief writes an origin-centered spherical slice in binary analysis format

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"
#include "parameter_input.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \class SphericalSlice
//! \brief constructs an equal-area angular grid and trilinearly samples local cells

class SphericalSlice {
 public:
  SphericalSlice(MeshBlockPack *pack, Real r, int nt, int np)
      : pmy_pack(pack), radius(r), ntheta(nt), nphi(np), nangles(nt*np),
        interp_indcs("sphslice_indices", nt*np, 4),
        interp_wghts("sphslice_weights", nt*np, 3),
        interp_vals("sphslice_values", nt*np) {
    Rebuild();
  }

  void Rebuild() {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    auto &size = pmy_pack->pmb->mb_size;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    int is = indcs.is;
    int js = indcs.js;
    int ks = indcs.ks;
    int nmb = pmy_pack->nmb_thispack;
    auto indices = interp_indcs.h_view;
    auto weights = interp_wghts.h_view;

    for (int it = 0; it < ntheta; ++it) {
      Real costheta = -1.0 + 2.0*(static_cast<Real>(it) + 0.5)/ntheta;
      Real sintheta = std::sqrt(std::max(static_cast<Real>(0.0),
                                         1.0 - costheta*costheta));
      for (int ip = 0; ip < nphi; ++ip) {
        Real phi = 2.0*M_PI*(static_cast<Real>(ip) + 0.5)/nphi;
        Real x = radius*sintheta*std::cos(phi);
        Real y = radius*sintheta*std::sin(phi);
        Real z = radius*costheta;
        int a = it*nphi + ip;
        int owner = -1;
        for (int m = 0; m < nmb; ++m) {
          auto block = size.h_view(m);
          if (x >= block.x1min && x < block.x1max &&
              y >= block.x2min && y < block.x2max &&
              z >= block.x3min && z < block.x3max) {
            owner = m;
            break;
          }
        }
        indices(a, 0) = owner;
        if (owner < 0) {
          indices(a, 1) = indices(a, 2) = indices(a, 3) = -1;
          weights(a, 0) = weights(a, 1) = weights(a, 2) = 0.0;
          continue;
        }
        auto block = size.h_view(owner);
        Real fi = (x - block.x1min)/block.dx1 - 0.5;
        Real fj = (y - block.x2min)/block.dx2 - 0.5;
        Real fk = (z - block.x3min)/block.dx3 - 0.5;
        int i0 = static_cast<int>(std::floor(fi));
        int j0 = static_cast<int>(std::floor(fj));
        int k0 = static_cast<int>(std::floor(fk));
        indices(a, 1) = i0 + is;
        indices(a, 2) = j0 + js;
        indices(a, 3) = k0 + ks;
        weights(a, 0) = fi - i0;
        weights(a, 1) = fj - j0;
        weights(a, 2) = fk - k0;
      }
    }
    interp_indcs.template modify<HostMemSpace>();
    interp_wghts.template modify<HostMemSpace>();
    interp_indcs.template sync<DevExeSpace>();
    interp_wghts.template sync<DevExeSpace>();
    owned_angles.clear();
    for (int a = 0; a < nangles; ++a) {
      if (indices(a, 0) >= 0) {
        owned_angles.push_back(static_cast<std::int32_t>(a));
      }
    }
  }

  void Interpolate(int variable, DvceArray5D<Real> &data) {
    auto indices = interp_indcs.d_view;
    auto weights = interp_wghts.d_view;
    auto values = interp_vals.d_view;
    par_for("sphslice_interpolate", DevExeSpace(), 0, nangles-1,
    KOKKOS_LAMBDA(int a) {
      int m = indices(a, 0);
      if (m < 0) {
        values(a) = 0.0;
        return;
      }
      int i = indices(a, 1);
      int j = indices(a, 2);
      int k = indices(a, 3);
      Real wx = weights(a, 0);
      Real wy = weights(a, 1);
      Real wz = weights(a, 2);
      Real c00 = data(m, variable, k,   j,   i)*(1.0-wx)
                 + data(m, variable, k,   j,   i+1)*wx;
      Real c10 = data(m, variable, k,   j+1, i)*(1.0-wx)
                 + data(m, variable, k,   j+1, i+1)*wx;
      Real c01 = data(m, variable, k+1, j,   i)*(1.0-wx)
                 + data(m, variable, k+1, j,   i+1)*wx;
      Real c11 = data(m, variable, k+1, j+1, i)*(1.0-wx)
                 + data(m, variable, k+1, j+1, i+1)*wx;
      Real c0 = c00*(1.0-wy) + c10*wy;
      Real c1 = c01*(1.0-wy) + c11*wy;
      values(a) = c0*(1.0-wz) + c1*wz;
    });
    interp_vals.template modify<DevExeSpace>();
    interp_vals.template sync<HostMemSpace>();
  }

  MeshBlockPack *pmy_pack;
  Real radius;
  int ntheta;
  int nphi;
  int nangles;
  DualArray2D<int> interp_indcs;
  DualArray2D<Real> interp_wghts;
  DualArray1D<Real> interp_vals;
  std::vector<std::int32_t> owned_angles;
};

//----------------------------------------------------------------------------------------
// Constructor and destructor

SphericalSliceOutput::SphericalSliceOutput(ParameterInput *pin, Mesh *pm,
                                           OutputParameters op)
    : BaseTypeOutput(pin, pm, op), psph(nullptr) {
  if (pm->mesh_indcs.nx2 <= 1 || pm->mesh_indcs.nx3 <= 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sphslice output requires a 3D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real radius = pin->GetReal(op.block_name, "slice_r");
  int ntheta = pin->GetOrAddInteger(op.block_name, "ntheta", 64);
  int nphi = pin->GetOrAddInteger(op.block_name, "nphi", 128);
  Real max_radius = std::min({pm->mesh_size.x1max, -pm->mesh_size.x1min,
                              pm->mesh_size.x2max, -pm->mesh_size.x2min,
                              pm->mesh_size.x3max, -pm->mesh_size.x3min});
  if (!(radius > 0.0 && radius < max_radius)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sphslice slice_r=" << radius << " in block '"
              << op.block_name << "' must lie strictly inside the origin-centered domain"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ntheta < 2 || nphi < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sphslice requires ntheta>=2 and nphi>=2" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (out_params.contains_derived) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "sphslice variable '" << out_params.variable
              << "' in block '" << out_params.block_name
              << "' requires derived-field interpolation, which is not supported "
              << "until ghost-zone-safe sampling is implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  mkdir("bin", 0775);
  if (IsSharded(op.shard_mode)) {
    std::string shard_path = "bin/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    mkdir(shard_path.c_str(), 0775);
  }
  psph = new SphericalSlice(pm->pmb_pack, radius, ntheta, nphi);
}

SphericalSliceOutput::~SphericalSliceOutput() {
  delete psph;
}

//----------------------------------------------------------------------------------------
// Data collection

void SphericalSliceOutput::LoadOutputData(Mesh *pm) {
  if (pm->adaptive) {
    psph->Rebuild();
  }
  int nvars = outvars.size();
  int npoints = psph->nangles;
  shard_owned_angles.clear();
  shard_values.clear();

  if (out_params.contains_derived) {
    out_params.i_derived = 0;
    ComputeDerivedVariable(out_params.variable, pm);
  }

  if (out_params.shard_mode == FileShardMode::shared) {
    Kokkos::realloc(outarray, nvars, 1, 1, psph->ntheta, psph->nphi);
    for (int n = 0; n < nvars; ++n) {
      psph->Interpolate(outvars[n].data_index, *(outvars[n].data_ptr));
      for (int a = 0; a < npoints; ++a) {
        outarray(n, 0, 0, a/psph->nphi, a%psph->nphi) =
            psph->interp_vals.h_view(a);
      }
    }
#if MPI_PARALLEL_ENABLED
    int count = nvars*npoints;
    if (global_variable::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, outarray.data(), count, MPI_ATHENA_REAL,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(outarray.data(), outarray.data(), count, MPI_ATHENA_REAL,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    return;
  }

  shard_owned_angles = psph->owned_angles;
  int local_points = static_cast<int>(shard_owned_angles.size());
  shard_values.resize(static_cast<std::size_t>(nvars)*local_points);
  for (int n = 0; n < nvars; ++n) {
    psph->Interpolate(outvars[n].data_index, *(outvars[n].data_ptr));
    for (int q = 0; q < local_points; ++q) {
      shard_values[static_cast<std::size_t>(n)*local_points + q] =
          static_cast<float>(psph->interp_vals.h_view(shard_owned_angles[q]));
    }
  }

#if MPI_PARALLEL_ENABLED
  if (IsNodeSharded(out_params.shard_mode)) {
    std::vector<int> counts;
    if (global_variable::node_rank == 0) {
      counts.resize(global_variable::node_size);
    }
    MPI_Gather(&local_points, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
               global_variable::node_comm);
    std::vector<int> offsets;
    int node_points = 0;
    if (global_variable::node_rank == 0) {
      offsets.resize(global_variable::node_size, 0);
      for (int r = 0; r < global_variable::node_size; ++r) {
        offsets[r] = node_points;
        node_points += counts[r];
      }
    }
    std::vector<std::int32_t> node_angles;
    std::vector<float> node_values;
    if (global_variable::node_rank == 0) {
      node_angles.resize(node_points);
      node_values.resize(static_cast<std::size_t>(nvars)*node_points);
    }
    MPI_Gatherv(shard_owned_angles.data(), local_points, MPI_INT32_T, node_angles.data(),
                counts.data(), offsets.data(), MPI_INT32_T, 0, global_variable::node_comm);
    for (int n = 0; n < nvars; ++n) {
      const float *send_values = local_points > 0
          ? &(shard_values[static_cast<std::size_t>(n)*local_points]) : nullptr;
      float *recv_values = (global_variable::node_rank == 0 && node_points > 0)
          ? &(node_values[static_cast<std::size_t>(n)*node_points]) : nullptr;
      MPI_Gatherv(send_values, local_points, MPI_FLOAT, recv_values, counts.data(),
                  offsets.data(), MPI_FLOAT, 0, global_variable::node_comm);
    }
    if (global_variable::node_rank == 0) {
      shard_owned_angles.swap(node_angles);
      shard_values.swap(node_values);
    }
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \brief Writes slice metadata followed by dense values or sparse rank records.

void SphericalSliceOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  bool sharded = IsSharded(out_params.shard_mode);
  bool i_write = IsRankSharded(out_params.shard_mode) ||
      (IsNodeSharded(out_params.shard_mode) && global_variable::node_rank == 0) ||
      (out_params.shard_mode == FileShardMode::shared && global_variable::my_rank == 0);
  if (i_write) {
    char number[7];
    char radius_token[32];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    std::snprintf(radius_token, sizeof(radius_token), "r_%g",
                  static_cast<double>(psph->radius));
    std::string path = "bin/";
    if (sharded) {
      path += ShardDirectoryName(out_params.shard_mode, global_variable::my_rank,
                                 global_variable::node_id) + "/";
    }
    std::string filename = path + out_params.file_basename + "." + out_params.file_id
        + "." + radius_token + number + ".sph.bin";
    std::FILE *output = std::fopen(filename.c_str(), "wb");
    if (output == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Cannot open sphslice output '" << filename << "'"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int nvars = outvars.size();
    int npoints = sharded ? static_cast<int>(shard_owned_angles.size()) : psph->nangles;
    std::stringstream parameters;
    pin->ParameterDump(parameters);
    std::string parameter_dump = parameters.str();
    std::stringstream header;
    header << "Athena spherical slice version=1.0\n"
           << "  layout=" << (sharded ? "sparse_angles" : "dense") << "\n"
           << "  distribution=" << ShardDistributionName(out_params.shard_mode) << "\n"
           << "  rank=" << global_variable::my_rank << "\n"
           << "  time=" << pm->time << "\n"
           << "  cycle=" << pm->ncycle << "\n"
           << "  radius=" << psph->radius << "\n"
           << "  ntheta=" << psph->ntheta << "\n"
           << "  nphi=" << psph->nphi << "\n"
           << "  size of variable=" << sizeof(float) << "\n"
           << "  number of variables=" << nvars << "\n"
           << "  npoints=" << npoints << "\n";
    if (IsNodeSharded(out_params.shard_mode)) {
      header << "  node=" << global_variable::node_id << "\n";
    }
    header << "  variables: ";
    for (int n = 0; n < nvars; ++n) {
      header << outvars[n].label << " ";
    }
    header << "\n  header offset=" << parameter_dump.size() << "\n";
    std::string metadata = header.str();
    std::fwrite(metadata.data(), sizeof(char), metadata.size(), output);
    std::fwrite(parameter_dump.data(), sizeof(char), parameter_dump.size(), output);
    if (sharded) {
      std::fwrite(shard_owned_angles.data(), sizeof(std::int32_t),
                  shard_owned_angles.size(), output);
      std::fwrite(shard_values.data(), sizeof(float), shard_values.size(), output);
    } else {
      std::vector<float> values(static_cast<std::size_t>(nvars)*psph->nangles);
      std::size_t q = 0;
      for (int n = 0; n < nvars; ++n) {
        for (int it = 0; it < psph->ntheta; ++it) {
          for (int ip = 0; ip < psph->nphi; ++ip) {
            values[q++] = static_cast<float>(outarray(n, 0, 0, it, ip));
          }
        }
      }
      std::fwrite(values.data(), sizeof(float), values.size(), output);
    }
    std::fclose(output);
  }

  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
