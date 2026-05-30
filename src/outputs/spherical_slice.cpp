//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_slice.cpp
//! \brief writes an origin-centered spherical slice in binary analysis format

#include <sys/stat.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
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

namespace {

[[noreturn]] void FatalSphericalSliceError(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

#if MPI_PARALLEL_ENABLED
void CheckMpi(int mpi_error, const char *context) {
  if (mpi_error != MPI_SUCCESS) {
    FatalSphericalSliceError(std::string(context) + " failed.");
  }
}
#endif

std::size_t CheckedAdd(std::size_t left, std::size_t right, const char *context) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    FatalSphericalSliceError(std::string(context) + " size overflow.");
  }
  return left + right;
}

std::size_t CheckedProduct(std::size_t left, std::size_t right, const char *context) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max()/left) {
    FatalSphericalSliceError(std::string(context) + " size overflow.");
  }
  return left*right;
}

int CountAsInt(std::size_t count, const char *context) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    FatalSphericalSliceError(std::string(context) + " exceeds int range.");
  }
  return static_cast<int>(count);
}

int CheckedAngularCount(int ntheta, int nphi) {
  if (ntheta <= 0 || nphi <= 0) {
    FatalSphericalSliceError("sphslice angular dimensions must be positive.");
  }
  return CountAsInt(CheckedProduct(static_cast<std::size_t>(ntheta),
                                   static_cast<std::size_t>(nphi),
                                   "sphslice angular grid"),
                    "sphslice angular grid");
}

void ValidateOwnedAngles(const std::vector<std::int32_t> &angles, int nangles,
                         const char *context) {
  std::vector<std::int32_t> sorted(angles);
  std::sort(sorted.begin(), sorted.end());
  for (std::size_t q = 0; q < sorted.size(); ++q) {
    if (sorted[q] < 0 || sorted[q] >= nangles) {
      FatalSphericalSliceError(std::string(context) + " has out-of-range angle index " +
                               std::to_string(sorted[q]) + ".");
    }
    if (q > 0 && sorted[q] == sorted[q - 1]) {
      FatalSphericalSliceError(std::string(context) + " has duplicate angle index " +
                               std::to_string(sorted[q]) + ".");
    }
  }
}

void SortAndValidateShardRecords(std::vector<std::int32_t> &angles,
                                 std::vector<float> &values, int nvars,
                                 int nangles, const char *context) {
  std::size_t npoints = angles.size();
  std::size_t expected_values = CheckedProduct(
      static_cast<std::size_t>(nvars), npoints, "sphslice sparse values");
  if (values.size() != expected_values) {
    FatalSphericalSliceError(std::string(context) + " value count does not match its "
                             "angular index count.");
  }
  ValidateOwnedAngles(angles, nangles, context);
  std::vector<std::size_t> order(npoints);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&angles](std::size_t left, std::size_t right) {
    return angles[left] < angles[right];
  });
  std::vector<std::int32_t> sorted_angles(npoints);
  std::vector<float> sorted_values(expected_values);
  for (std::size_t q = 0; q < npoints; ++q) {
    sorted_angles[q] = angles[order[q]];
    for (int n = 0; n < nvars; ++n) {
      sorted_values[static_cast<std::size_t>(n)*npoints + q] =
          values[static_cast<std::size_t>(n)*npoints + order[q]];
    }
  }
  angles.swap(sorted_angles);
  values.swap(sorted_values);
}

void ValidateGlobalOwnership(const std::vector<std::int32_t> &angles, int nangles) {
  ValidateOwnedAngles(angles, nangles, "sphslice local ownership");
  std::vector<int> local_owners(nangles, 0);
  std::vector<int> global_owners(nangles, 0);
  for (std::int32_t angle : angles) {
    local_owners[angle] = 1;
  }
#if MPI_PARALLEL_ENABLED
  if (MPI_Allreduce(local_owners.data(), global_owners.data(), nangles, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    FatalSphericalSliceError("Could not reduce sphslice angular ownership.");
  }
#else
  global_owners.swap(local_owners);
#endif
  for (int angle = 0; angle < nangles; ++angle) {
    if (global_owners[angle] != 1) {
      FatalSphericalSliceError("sphslice angle index " + std::to_string(angle) +
                               " has " + std::to_string(global_owners[angle]) +
                               " global owners; expected exactly one.");
    }
  }
}

void CheckedFileWrite(std::FILE *output, const void *data, std::size_t element_size,
                      std::size_t count, const std::string &filename,
                      const char *context) {
  if (count == 0) {
    return;
  }
  if (std::fwrite(data, element_size, count, output) != count) {
    std::remove(filename.c_str());
    FatalSphericalSliceError(std::string(context) + " was not written completely to '" +
                             filename + "'.");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \class SphericalSlice
//! \brief constructs an equal-area angular grid and trilinearly samples local cells

class SphericalSlice {
 public:
  SphericalSlice(MeshBlockPack *pack, Real r, int nt, int np)
      : pmy_pack(pack), radius(r), ntheta(nt), nphi(np),
        nangles(CheckedAngularCount(nt, np)),
        interp_indcs("sphslice_indices", nangles, 4),
        interp_wghts("sphslice_weights", nangles, 3),
        interp_vals("sphslice_values", nangles) {
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
  ValidateGlobalOwnership(psph->owned_angles, psph->nangles);

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
    int count = CountAsInt(CheckedProduct(static_cast<std::size_t>(nvars),
                                          static_cast<std::size_t>(npoints),
                                          "dense sphslice values"),
                           "dense sphslice values");
    if (global_variable::my_rank == 0) {
      CheckMpi(MPI_Reduce(MPI_IN_PLACE, outarray.data(), count, MPI_ATHENA_REAL,
                          MPI_SUM, 0, MPI_COMM_WORLD),
               "MPI_Reduce for dense sphslice values");
    } else {
      CheckMpi(MPI_Reduce(outarray.data(), outarray.data(), count, MPI_ATHENA_REAL,
                          MPI_SUM, 0, MPI_COMM_WORLD),
               "MPI_Reduce for dense sphslice values");
    }
#endif
    return;
  }

  shard_owned_angles = psph->owned_angles;
  int local_points = CountAsInt(shard_owned_angles.size(), "local sphslice point count");
  shard_values.resize(CheckedProduct(static_cast<std::size_t>(nvars),
                                     static_cast<std::size_t>(local_points),
                                     "local sphslice values"));
  for (int n = 0; n < nvars; ++n) {
    psph->Interpolate(outvars[n].data_index, *(outvars[n].data_ptr));
    for (int q = 0; q < local_points; ++q) {
      shard_values[static_cast<std::size_t>(n)*local_points + q] =
          static_cast<float>(psph->interp_vals.h_view(shard_owned_angles[q]));
    }
  }
  SortAndValidateShardRecords(shard_owned_angles, shard_values, nvars, psph->nangles,
                              "local sphslice shard");

#if MPI_PARALLEL_ENABLED
  if (IsNodeSharded(out_params.shard_mode)) {
    std::vector<int> counts;
    if (global_variable::node_rank == 0) {
      counts.resize(global_variable::node_size);
    }
    CheckMpi(MPI_Gather(&local_points, 1, MPI_INT,
                        global_variable::node_rank == 0 ? counts.data() : nullptr,
                        1, MPI_INT, 0, global_variable::node_comm),
             "MPI_Gather for node sphslice point counts");
    std::vector<int> offsets;
    std::size_t node_points = 0;
    if (global_variable::node_rank == 0) {
      offsets.resize(global_variable::node_size, 0);
      for (int r = 0; r < global_variable::node_size; ++r) {
        offsets[r] = CountAsInt(node_points, "node sphslice point offset");
        node_points = CheckedAdd(node_points, static_cast<std::size_t>(counts[r]),
                                 "node sphslice point count");
      }
    }
    std::vector<std::int32_t> node_angles;
    std::vector<float> node_values;
    if (global_variable::node_rank == 0) {
      node_angles.resize(node_points);
      node_values.resize(CheckedProduct(static_cast<std::size_t>(nvars), node_points,
                                        "node sphslice values"));
    }
    CheckMpi(MPI_Gatherv(shard_owned_angles.data(), local_points, MPI_INT32_T,
                         global_variable::node_rank == 0 ? node_angles.data() : nullptr,
                         global_variable::node_rank == 0 ? counts.data() : nullptr,
                         global_variable::node_rank == 0 ? offsets.data() : nullptr,
                         MPI_INT32_T, 0, global_variable::node_comm),
             "MPI_Gatherv for node sphslice angles");
    for (int n = 0; n < nvars; ++n) {
      const float *send_values = local_points > 0
          ? &(shard_values[static_cast<std::size_t>(n)*local_points]) : nullptr;
      float *recv_values = (global_variable::node_rank == 0 && node_points > 0)
          ? &(node_values[static_cast<std::size_t>(n)*node_points]) : nullptr;
      CheckMpi(MPI_Gatherv(send_values, local_points, MPI_FLOAT, recv_values,
                           global_variable::node_rank == 0 ? counts.data() : nullptr,
                           global_variable::node_rank == 0 ? offsets.data() : nullptr,
                           MPI_FLOAT, 0, global_variable::node_comm),
               "MPI_Gatherv for node sphslice values");
    }
    if (global_variable::node_rank == 0) {
      shard_owned_angles.swap(node_angles);
      shard_values.swap(node_values);
      SortAndValidateShardRecords(shard_owned_angles, shard_values, nvars, psph->nangles,
                                  "node sphslice shard");
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
    int nvars = CountAsInt(outvars.size(), "sphslice variable count");
    if (sharded) {
      SortAndValidateShardRecords(shard_owned_angles, shard_values, nvars, psph->nangles,
                                  "published sphslice shard");
    }
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
    std::string temporary_filename = filename + ".tmp";
    std::FILE *output = std::fopen(temporary_filename.c_str(), "wb");
    if (output == nullptr) {
      FatalSphericalSliceError(
          "Cannot open sphslice output '" + temporary_filename + "'.");
    }
    int npoints = sharded ? CountAsInt(shard_owned_angles.size(), "sphslice point count")
                          : psph->nangles;
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
      header << "  node=" << global_variable::node_id << "\n"
             << "  number of nodes=" << global_variable::nnodes << "\n";
    } else if (IsRankSharded(out_params.shard_mode)) {
      header << "  number of ranks=" << global_variable::nranks << "\n";
    }
    header << "  variables: ";
    for (int n = 0; n < nvars; ++n) {
      header << outvars[n].label << " ";
    }
    header << "\n  header offset=" << parameter_dump.size() << "\n";
    std::string metadata = header.str();
    CheckedFileWrite(output, metadata.data(), sizeof(char), metadata.size(),
                     temporary_filename, "sphslice metadata");
    CheckedFileWrite(output, parameter_dump.data(), sizeof(char), parameter_dump.size(),
                     temporary_filename, "sphslice input header");
    if (sharded) {
      CheckedFileWrite(output, shard_owned_angles.data(), sizeof(std::int32_t),
                       shard_owned_angles.size(), temporary_filename,
                       "sphslice angular indices");
      CheckedFileWrite(output, shard_values.data(), sizeof(float), shard_values.size(),
                       temporary_filename, "sphslice sparse values");
    } else {
      std::vector<float> values(CheckedProduct(static_cast<std::size_t>(nvars),
                                               static_cast<std::size_t>(psph->nangles),
                                               "dense sphslice values"));
      std::size_t q = 0;
      for (int n = 0; n < nvars; ++n) {
        for (int it = 0; it < psph->ntheta; ++it) {
          for (int ip = 0; ip < psph->nphi; ++ip) {
            values[q++] = static_cast<float>(outarray(n, 0, 0, it, ip));
          }
        }
      }
      CheckedFileWrite(output, values.data(), sizeof(float), values.size(),
                       temporary_filename, "sphslice dense values");
    }
    if (std::fclose(output) != 0) {
      std::remove(temporary_filename.c_str());
      FatalSphericalSliceError("Could not close sphslice output '" + temporary_filename +
                               "'.");
    }
    if (std::rename(temporary_filename.c_str(), filename.c_str()) != 0) {
      int rename_errno = errno;
      std::remove(temporary_filename.c_str());
      FatalSphericalSliceError(
          "Could not atomically publish sphslice output '" + filename +
          "': " + std::strerror(rename_errno));
    }
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
