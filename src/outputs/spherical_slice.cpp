//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_slice.cpp
//! \brief writes an origin-centered spherical slice in binary analysis format

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"
#include "parameter_input.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

[[noreturn]] void FatalSphericalSliceError(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" + message);
}

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

std::size_t SphericalSlicePersistentBytes(int nangles) {
  std::size_t angles = static_cast<std::size_t>(nangles);
  std::size_t fixed_per_angle = 8*sizeof(int) + 8*sizeof(Real) +
                                sizeof(std::int32_t);
  return output_file_utils::CheckedSizeProduct(
      angles, fixed_per_angle,
      "sphslice allocation", FatalSphericalSliceError);
}

std::size_t SphericalSliceSparseBytes(std::size_t npoints, std::size_t nvars) {
  return CheckedAdd(
      CheckedProduct(npoints, sizeof(std::int32_t), "sphslice sparse angles"),
      CheckedProduct(CheckedProduct(nvars, npoints, "sphslice sparse values"),
                     sizeof(float), "sphslice sparse values"),
      "sphslice sparse staging");
}

std::size_t SphericalSliceDenseBytes(std::size_t npoints, std::size_t nvars,
                                     std::size_t value_size,
                                     const char *context) {
  return CheckedProduct(CheckedProduct(nvars, npoints, context), value_size, context);
}

std::size_t SphericalSliceSortBytes(std::size_t npoints, std::size_t nvars) {
  return CheckedAdd(
      CheckedAdd(
          CheckedProduct(CheckedProduct(2, npoints, "sphslice sort angles"),
                         sizeof(std::int32_t), "sphslice sort angles"),
          CheckedProduct(npoints, sizeof(std::size_t), "sphslice sort order"),
          "sphslice sort staging"),
      CheckedProduct(CheckedProduct(nvars, npoints, "sphslice sort values"),
                     sizeof(float), "sphslice sort values"),
      "sphslice sort staging");
}

void RequireSphericalSliceBudget(std::size_t persistent_bytes,
                                 std::size_t staging_bytes,
                                 std::size_t max_writer_allocation_bytes,
                                 const char *context) {
  output_file_utils::RequireAllocationBudget(
      CheckedAdd(persistent_bytes, staging_bytes, context),
      max_writer_allocation_bytes, context, FatalSphericalSliceError);
}

class CountingStreamBuffer : public std::streambuf {
 public:
  std::size_t size() const { return size_; }

 protected:
  std::streamsize xsputn(const char *, std::streamsize count) override {
    if (count < 0) {
      FatalSphericalSliceError("sphslice serialization received a negative count.");
    }
    size_ = CheckedAdd(size_, static_cast<std::size_t>(count),
                       "sphslice serialization staging");
    return count;
  }

  int_type overflow(int_type character) override {
    if (!traits_type::eq_int_type(character, traits_type::eof())) {
      size_ = CheckedAdd(size_, 1, "sphslice serialization staging");
    }
    return traits_type::not_eof(character);
  }

 private:
  std::size_t size_ = 0;
};

class FileStreamBuffer : public std::streambuf {
 public:
  explicit FileStreamBuffer(std::FILE *output) : output_(output) {}

 protected:
  std::streamsize xsputn(const char *data, std::streamsize count) override {
    if (count < 0) return 0;
    return static_cast<std::streamsize>(
        std::fwrite(data, sizeof(char), static_cast<std::size_t>(count), output_));
  }

  int_type overflow(int_type character) override {
    if (traits_type::eq_int_type(character, traits_type::eof())) {
      return traits_type::not_eof(character);
    }
    return std::fputc(traits_type::to_char_type(character), output_) == EOF
        ? traits_type::eof() : character;
  }

  int sync() override {
    return std::fflush(output_) == 0 ? 0 : -1;
  }

 private:
  std::FILE *output_;
};

std::size_t ParameterDumpSerializedSize(const ParameterInput &pin) {
  constexpr char marker[] =
      "#------------------------- PAR_DUMP -------------------------\n";
  std::size_t bytes = sizeof(marker) - 1;
  for (const auto &block : pin.block) {
    bytes = CheckedAdd(bytes, CheckedAdd(block.block_name.size(), 3,
                                         "sphslice parameter block"),
                       "sphslice parameter dump");
    for (const auto &line : block.line) {
      std::size_t line_bytes = CheckedAdd(block.max_len_parname, 1,
                                          "sphslice parameter name");
      line_bytes = CheckedAdd(line_bytes, 2, "sphslice parameter assignment");
      line_bytes = CheckedAdd(line_bytes,
                              CheckedAdd(block.max_len_parvalue, 1,
                                         "sphslice parameter value"),
                              "sphslice parameter line");
      line_bytes = CheckedAdd(line_bytes, line.param_comment.size(),
                              "sphslice parameter comment");
      bytes = CheckedAdd(bytes, CheckedAdd(line_bytes, 1, "sphslice parameter line"),
                         "sphslice parameter dump");
    }
  }
  bytes = CheckedAdd(bytes, sizeof(marker) - 1, "sphslice parameter dump");
  return CheckedAdd(bytes, sizeof("<par_end>\n") - 1, "sphslice parameter dump");
}

std::size_t StringStorageBytes(std::size_t characters, const char *context) {
  return CheckedAdd(characters, 1, context);
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
  for (float value : values) {
    if (!std::isfinite(value)) {
      FatalSphericalSliceError(std::string(context) + " contains a non-finite value.");
    }
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
  mpi_utils::CheckMpi(
      MPI_Allreduce(local_owners.data(), global_owners.data(), nangles, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD),
      "MPI_Allreduce for sphslice angular ownership");
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
    output_file_utils::DiscardOwnedPath(filename);
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
    owned_angles.reserve(nangles);
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
    : BaseTypeOutput(pin, pm, op), psph(nullptr),
      persistent_writer_allocation_bytes(0) {
  if (pm->mesh_indcs.nx2 <= 1 || pm->mesh_indcs.nx3 <= 1) {
    FatalSphericalSliceError("sphslice output requires a 3D mesh.");
  }
  Real radius = pin->GetReal(op.block_name, "slice_r");
  int ntheta = pin->GetOrAddInteger(op.block_name, "ntheta", 64);
  int nphi = pin->GetOrAddInteger(op.block_name, "nphi", 128);
  int configured_limit = pin->GetOrAddInteger(
      op.block_name, "max_writer_allocation_bytes",
      static_cast<int>(output_file_utils::kDefaultMaxWriterAllocationBytes));
  if (configured_limit <= 0) {
    FatalSphericalSliceError("sphslice max_writer_allocation_bytes must be positive.");
  }
  max_writer_allocation_bytes = static_cast<std::size_t>(configured_limit);
  Real max_radius = std::min({pm->mesh_size.x1max, -pm->mesh_size.x1min,
                              pm->mesh_size.x2max, -pm->mesh_size.x2min,
                              pm->mesh_size.x3max, -pm->mesh_size.x3min});
  if (!(radius > 0.0 && radius < max_radius)) {
    FatalSphericalSliceError("sphslice slice_r=" + std::to_string(radius) +
                             " in block '" + op.block_name +
                             "' must lie strictly inside the origin-centered domain.");
  }
  if (ntheta < 2 || nphi < 2) {
    FatalSphericalSliceError("sphslice requires ntheta>=2 and nphi>=2.");
  }
  if (out_params.contains_derived) {
    FatalSphericalSliceError(
        "sphslice variable '" + out_params.variable + "' in block '" +
        out_params.block_name +
        "' requires derived-field interpolation, which is not supported until "
        "ghost-zone-safe sampling is implemented.");
  }
  int nangles = CheckedAngularCount(ntheta, nphi);
  persistent_writer_allocation_bytes = SphericalSlicePersistentBytes(nangles);
  output_file_utils::RequireAllocationBudget(
      persistent_writer_allocation_bytes,
      max_writer_allocation_bytes, "sphslice allocation", FatalSphericalSliceError);
  output_file_utils::EnsureDirectory("bin", 0775, "sphslice output",
                                     FatalSphericalSliceError);
  if (IsSharded(op.shard_mode)) {
    std::string shard_path = "bin/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    output_file_utils::EnsureDirectory(shard_path, 0775, "sphslice output",
                                       FatalSphericalSliceError);
  }
  psph = std::make_unique<SphericalSlice>(pm->pmb_pack, radius, ntheta, nphi);
}

SphericalSliceOutput::~SphericalSliceOutput() = default;

//----------------------------------------------------------------------------------------
// Data collection

void SphericalSliceOutput::LoadOutputData(Mesh *pm) {
  std::size_t retained_output_bytes = CheckedAdd(
      CheckedProduct(outarray.size(), sizeof(Real), "retained dense sphslice values"),
      CheckedAdd(
          CheckedProduct(shard_owned_angles.capacity(), sizeof(std::int32_t),
                         "retained sparse sphslice angles"),
          CheckedProduct(shard_values.capacity(), sizeof(float),
                         "retained sparse sphslice values"),
          "retained sparse sphslice values"),
      "retained sphslice values");
  int npoints = psph->nangles;
  std::size_t angular_points = static_cast<std::size_t>(npoints);
  std::size_t ownership_staging_bytes = CheckedAdd(
      CheckedProduct(angular_points, sizeof(std::int32_t),
                     "sphslice ownership validation"),
      CheckedProduct(CheckedProduct(2, angular_points, "sphslice ownership validation"),
                     sizeof(int), "sphslice ownership validation"),
      "sphslice ownership validation");
  RequireSphericalSliceBudget(
      CheckedAdd(persistent_writer_allocation_bytes, retained_output_bytes,
                 "sphslice ownership validation"),
      ownership_staging_bytes, max_writer_allocation_bytes,
      "sphslice ownership validation");
  if (pm->adaptive) {
    psph->Rebuild();
  }
  int nvars = CountAsInt(outvars.size(), "sphslice variable count");
  shard_owned_angles.clear();
  shard_values.clear();
  ValidateGlobalOwnership(psph->owned_angles, psph->nangles);

  if (out_params.contains_derived) {
    out_params.i_derived = 0;
    ComputeDerivedVariable(out_params.variable, pm);
  }

  if (out_params.shard_mode == FileShardMode::shared) {
    RequireSphericalSliceBudget(
        persistent_writer_allocation_bytes,
        SphericalSliceDenseBytes(angular_points, static_cast<std::size_t>(nvars),
                                 sizeof(Real), "dense sphslice values"),
        max_writer_allocation_bytes, "dense sphslice values");
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
      mpi_utils::CheckMpi(MPI_Reduce(MPI_IN_PLACE, outarray.data(), count,
                                     MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD),
                          "MPI_Reduce for dense sphslice values");
    } else {
      mpi_utils::CheckMpi(MPI_Reduce(outarray.data(), outarray.data(), count,
                                     MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD),
                          "MPI_Reduce for dense sphslice values");
    }
#endif
    if (global_variable::my_rank == 0) {
      for (int n = 0; n < nvars; ++n) {
        for (int a = 0; a < npoints; ++a) {
          if (!std::isfinite(outarray(n, 0, 0, a/psph->nphi, a%psph->nphi))) {
            FatalSphericalSliceError("dense sphslice values contain a non-finite value.");
          }
        }
      }
    }
    return;
  }

  std::size_t local_points_size = psph->owned_angles.size();
  std::size_t local_sparse_bytes =
      SphericalSliceSparseBytes(local_points_size, static_cast<std::size_t>(nvars));
  RequireSphericalSliceBudget(
      persistent_writer_allocation_bytes,
      CheckedAdd(local_sparse_bytes,
                 SphericalSliceSortBytes(local_points_size,
                                         static_cast<std::size_t>(nvars)),
                 "local sphslice sort staging"),
      max_writer_allocation_bytes, "local sphslice sort staging");
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
      std::size_t node_metadata_bytes = CheckedProduct(
          CheckedProduct(2, static_cast<std::size_t>(global_variable::node_size),
                         "node sphslice metadata"),
          sizeof(int), "node sphslice metadata");
      RequireSphericalSliceBudget(
          persistent_writer_allocation_bytes,
          CheckedAdd(local_sparse_bytes, node_metadata_bytes,
                     "node sphslice metadata"),
          max_writer_allocation_bytes, "node sphslice metadata");
      counts.resize(global_variable::node_size);
    }
    mpi_utils::CheckMpi(
        MPI_Gather(&local_points, 1, MPI_INT,
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
      std::size_t node_metadata_bytes = CheckedProduct(
          CheckedProduct(2, static_cast<std::size_t>(global_variable::node_size),
                         "node sphslice metadata"),
          sizeof(int), "node sphslice metadata");
      std::size_t node_sparse_bytes =
          SphericalSliceSparseBytes(node_points, static_cast<std::size_t>(nvars));
      RequireSphericalSliceBudget(
          persistent_writer_allocation_bytes,
          CheckedAdd(
              CheckedAdd(CheckedAdd(local_sparse_bytes, node_metadata_bytes,
                                    "node sphslice sort staging"),
                         node_sparse_bytes, "node sphslice sort staging"),
              SphericalSliceSortBytes(node_points, static_cast<std::size_t>(nvars)),
              "node sphslice sort staging"),
          max_writer_allocation_bytes, "node sphslice sort staging");
      node_angles.resize(node_points);
      node_values.resize(CheckedProduct(static_cast<std::size_t>(nvars), node_points,
                                        "node sphslice values"));
    }
    mpi_utils::CheckMpi(
        MPI_Gatherv(shard_owned_angles.data(), local_points, MPI_INT32_T,
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
      mpi_utils::CheckMpi(
          MPI_Gatherv(send_values, local_points, MPI_FLOAT, recv_values,
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
      std::size_t sparse_bytes = SphericalSliceSparseBytes(
          shard_owned_angles.size(), static_cast<std::size_t>(nvars));
      RequireSphericalSliceBudget(
          persistent_writer_allocation_bytes,
          CheckedAdd(sparse_bytes,
                     SphericalSliceSortBytes(shard_owned_angles.size(),
                                             static_cast<std::size_t>(nvars)),
                     "published sphslice sort staging"),
          max_writer_allocation_bytes, "published sphslice sort staging");
      SortAndValidateShardRecords(shard_owned_angles, shard_values, nvars, psph->nangles,
                                  "published sphslice shard");
    }
    int npoints = sharded ? CountAsInt(shard_owned_angles.size(), "sphslice point count")
                          : psph->nangles;
    std::size_t parameter_dump_size = ParameterDumpSerializedSize(*pin);
    auto write_metadata = [&](std::ostream &header) {
      header << "Athena spherical slice version=1.0\n"
             << "  layout=" << (sharded ? "sparse_angles" : "dense") << "\n"
             << "  distribution=" << ShardDistributionName(out_params.shard_mode) << "\n"
             << "  rank=" << global_variable::my_rank << "\n"
             << "  time=" << pm->time << "\n"
             << "  cycle=" << pm->ncycle << "\n"
             << "  radius=" << std::scientific
             << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1)
             << psph->radius << std::defaultfloat << "\n"
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
      header << "\n  header offset=" << parameter_dump_size << "\n";
    };
    CountingStreamBuffer metadata_counter_buffer;
    std::ostream metadata_counter(&metadata_counter_buffer);
    write_metadata(metadata_counter);
    if (!metadata_counter) {
      FatalSphericalSliceError("Could not size sphslice metadata.");
    }
    std::size_t serialization_bytes =
        CheckedAdd(metadata_counter_buffer.size(), parameter_dump_size,
                   "sphslice serialization staging");
    constexpr std::size_t kMaxRadiusTokenCharacters = 64;
    constexpr std::size_t kMaxSequenceTokenCharacters = 10;
    constexpr std::size_t kShardDirectoryCharacters = sizeof("node_00000000/") - 1;
    std::size_t shard_directory_characters = sharded ? kShardDirectoryCharacters : 0;
    std::size_t filename_characters = sizeof("bin/") - 1;
    filename_characters = CheckedAdd(filename_characters, shard_directory_characters,
                                     "sphslice filename");
    filename_characters = CheckedAdd(filename_characters, out_params.file_basename.size(),
                                     "sphslice filename");
    filename_characters = CheckedAdd(filename_characters, out_params.file_id.size(),
                                     "sphslice filename");
    filename_characters = CheckedAdd(filename_characters, kMaxRadiusTokenCharacters,
                                     "sphslice filename");
    filename_characters = CheckedAdd(filename_characters, kMaxSequenceTokenCharacters,
                                     "sphslice filename");
    filename_characters = CheckedAdd(filename_characters, sizeof("....sph.bin") - 1,
                                     "sphslice filename");
    std::size_t path_staging_bytes = CheckedAdd(
        StringStorageBytes(filename_characters, "sphslice filename"),
        StringStorageBytes(CheckedAdd(filename_characters, sizeof(".tmp") - 1,
                                      "sphslice temporary filename"),
                           "sphslice temporary filename"),
        "sphslice path staging");
    path_staging_bytes = CheckedAdd(
        path_staging_bytes,
        CheckedAdd(StringStorageBytes(kMaxSequenceTokenCharacters,
                                      "sphslice sequence token"),
                   StringStorageBytes(kMaxRadiusTokenCharacters,
                                      "sphslice radius token"),
                   "sphslice token staging"),
        "sphslice path staging");
    path_staging_bytes = CheckedAdd(
        path_staging_bytes,
        StringStorageBytes(shard_directory_characters, "sphslice shard directory"),
        "sphslice path staging");
    serialization_bytes = CheckedAdd(serialization_bytes, path_staging_bytes,
                                     "sphslice serialization staging");
    if (sharded) {
      serialization_bytes = CheckedAdd(
          serialization_bytes,
          SphericalSliceSparseBytes(shard_owned_angles.size(),
                                    static_cast<std::size_t>(nvars)),
          "sphslice serialization staging");
    } else {
      std::size_t angular_points = static_cast<std::size_t>(psph->nangles);
      serialization_bytes = CheckedAdd(
          serialization_bytes,
          CheckedAdd(
              SphericalSliceDenseBytes(angular_points, static_cast<std::size_t>(nvars),
                                       sizeof(Real), "dense sphslice values"),
              SphericalSliceDenseBytes(angular_points, static_cast<std::size_t>(nvars),
                                       sizeof(float), "dense sphslice serialization"),
              "sphslice serialization staging"),
          "sphslice serialization staging");
    }
    RequireSphericalSliceBudget(persistent_writer_allocation_bytes,
                                serialization_bytes,
                                max_writer_allocation_bytes,
                                "sphslice serialization staging");
    std::vector<float> dense_values;
    if (!sharded) {
      dense_values.resize(CheckedProduct(static_cast<std::size_t>(nvars),
                                         static_cast<std::size_t>(psph->nangles),
                                         "dense sphslice values"));
      std::size_t q = 0;
      for (int n = 0; n < nvars; ++n) {
        for (int it = 0; it < psph->ntheta; ++it) {
          for (int ip = 0; ip < psph->nphi; ++ip) {
            float value = static_cast<float>(outarray(n, 0, 0, it, ip));
            if (!std::isfinite(value)) {
              FatalSphericalSliceError(
                  "dense sphslice values contain a non-finite serialized value.");
            }
            dense_values[q++] = value;
          }
        }
      }
    }
    std::string number = output_file_utils::FormatSequence(
        out_params.file_number, "sphslice output", FatalSphericalSliceError);
    std::string radius_token =
        output_file_utils::FormatSphericalSliceRadius(psph->radius);
    std::string shard_directory;
    if (sharded) {
      shard_directory = ShardDirectoryName(
          out_params.shard_mode, global_variable::my_rank, global_variable::node_id)
          + "/";
    }
    std::string filename = "bin/" + shard_directory + out_params.file_basename + "."
        + out_params.file_id + "." + radius_token + "." + number + ".sph.bin";
    std::string temporary_filename = output_file_utils::TemporaryPath(filename);
    std::FILE *output = std::fopen(temporary_filename.c_str(), "wb");
    if (output == nullptr) {
      FatalSphericalSliceError(
          "Cannot open sphslice output '" + temporary_filename + "'.");
    }
    FileStreamBuffer header_buffer(output);
    std::ostream serialized_header(&header_buffer);
    write_metadata(serialized_header);
    pin->ParameterDump(serialized_header);
    serialized_header.flush();
    if (!serialized_header) {
      std::fclose(output);
      output_file_utils::DiscardOwnedPath(temporary_filename);
      FatalSphericalSliceError("Could not write sphslice serialized header to '" +
                               temporary_filename + "'.");
    }
    if (sharded) {
      CheckedFileWrite(output, shard_owned_angles.data(), sizeof(std::int32_t),
                       shard_owned_angles.size(), temporary_filename,
                       "sphslice angular indices");
      CheckedFileWrite(output, shard_values.data(), sizeof(float), shard_values.size(),
                       temporary_filename, "sphslice sparse values");
    } else {
      CheckedFileWrite(output, dense_values.data(), sizeof(float), dense_values.size(),
                       temporary_filename, "sphslice dense values");
    }
    if (std::fclose(output) != 0) {
      output_file_utils::DiscardOwnedPath(temporary_filename);
      FatalSphericalSliceError("Could not close sphslice output '" + temporary_filename +
                               "'.");
    }
    output_file_utils::PublishTemporaryFile(
        temporary_filename, filename, "sphslice output", FatalSphericalSliceError);
  }

  // These buffers are write-only staging. Release them so the next output cycle
  // starts without carrying a previous dense or sparse allocation into admission.
  outarray = HostArray5D<Real>();
  std::vector<std::int32_t>().swap(shard_owned_angles);
  std::vector<float>().swap(shard_values);
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "sphslice output", FatalSphericalSliceError);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
