//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file projection.cpp
//! \brief writes Cartesian projections as additive native-AMR patches or guarded maps

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "file_sharding.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "outputs.hpp"
#include "parameter_input.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

Real AxisMinimum(const RegionSize &size, int axis) {
  if (axis == 1) return size.x1min;
  if (axis == 2) return size.x2min;
  return size.x3min;
}

Real AxisMaximum(const RegionSize &size, int axis) {
  if (axis == 1) return size.x1max;
  if (axis == 2) return size.x2max;
  return size.x3max;
}

int AxisCellCount(const RegionIndcs &indcs, int axis) {
  if (axis == 1) return indcs.nx1;
  if (axis == 2) return indcs.nx2;
  return indcs.nx3;
}

int AxisStartIndex(const RegionIndcs &indcs, int axis) {
  if (axis == 1) return indcs.is;
  if (axis == 2) return indcs.js;
  return indcs.ks;
}

Real BlockMinimum(const OutputMeshBlockInfo &mb, int axis) {
  if (axis == 1) return mb.x1min;
  if (axis == 2) return mb.x2min;
  return mb.x3min;
}

Real BlockMaximum(const OutputMeshBlockInfo &mb, int axis) {
  if (axis == 1) return mb.x1max;
  if (axis == 2) return mb.x2max;
  return mb.x3max;
}

Real CellLowerEdge(const OutputMeshBlockInfo &mb, int axis, int cell, int ncells) {
  return BlockMinimum(mb, axis) +
      static_cast<Real>(cell) * (BlockMaximum(mb, axis) - BlockMinimum(mb, axis)) /
      static_cast<Real>(ncells);
}

Real Overlap(Real amin, Real amax, Real bmin, Real bmax) {
  return std::max(static_cast<Real>(0.0),
                  std::min(amax, bmax) - std::max(amin, bmin));
}

bool ProjectionStatsEnabled() {
  const char *env = std::getenv("ATHENAK_OUTPUT_IO_STATS");
  return (env != nullptr && env[0] != '\0' && env[0] != '0');
}

template <typename T>
void AppendBytes(std::vector<char> &buffer, const T &value) {
  std::size_t start = buffer.size();
  buffer.resize(start + sizeof(T));
  std::memcpy(buffer.data() + start, &value, sizeof(T));
}

void AppendMomentValues(std::vector<char> &buffer, const std::vector<Real> &values) {
  for (Real value : values) {
    double stored = static_cast<double>(value);
    AppendBytes(buffer, stored);
  }
}

void CheckedFwrite(const void *ptr, std::size_t size, std::size_t count,
                   std::FILE *stream, const std::string &path) {
  if (count > 0 && std::fwrite(ptr, size, count, stream) != count) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Failed writing projection output '" << path << "'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// Constructor: reads projection geometry, layout, statistic, and weighting options.

ProjectionOutput::ProjectionOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
    BaseTypeOutput(pin, pm, op),
    projection_axes_{false, false, false},
    n_projection_axes_(0),
    n_image_axes_(0),
    projection_level_(-1),
    uniform_max_megabytes_(0),
    projection_min_{0.0, 0.0, 0.0},
    projection_max_{0.0, 0.0, 0.0},
    layout_(Layout::native_amr),
    layout_name_("native_amr"),
    weighting_(Weighting::integral),
    weighting_name_("integral"),
    statistic_(Statistic::value),
    statistic_name_("value"),
    projection_axes_name_(""),
    density_ptr_(nullptr),
    projected_data_("projection", 1, 1, 1),
    second_moment_("projection_second_moment", 1, 1, 1),
    normalization_("projection_normalization", 1, 1),
    local_patch_count_(0),
    shard_patch_count_(0),
    local_staged_cells_(0),
    shard_staged_cells_(0),
    local_payload_bytes_(0),
    communicated_bytes_(0),
    peak_output_bytes_(0),
    output_elapsed_seconds_(0.0) {
  mkdir("proj", 0775);
  if (op.file_shard_mode != FileShardMode::shared) {
    std::string shard_dir = std::string("proj/") + ShardDirectoryName(op.file_shard_mode);
    mkdir(shard_dir.c_str(), 0775);
  }

  if (out_params.include_gzs || out_params.gid >= 0 || out_params.slice1 ||
      out_params.slice2 || out_params.slice3) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' does not accept ghost_zones, gid, or slice_x* parameters."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  layout_name_ = pin->GetOrAddString(op.block_name, "projection_layout", "native_amr");
  if (layout_name_ == "native_amr") {
    layout_ = Layout::native_amr;
    if (pin->DoesParameterExist(op.block_name, "projection_level")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Native-AMR projection block '" << out_params.block_name
                << "' stores leaf patches and does not accept projection_level. "
                << "Choose composition resolution in the reader." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else if (layout_name_ == "uniform") {
    layout_ = Layout::uniform;
    if (out_params.file_shard_mode != FileShardMode::shared) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Guarded uniform projection block '"
                << out_params.block_name << "' is a shared quicklook path and cannot "
                << "set single_file_per_rank or single_file_per_node." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!pin->DoesParameterExist(op.block_name, "projection_level")) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Uniform projection block '" << out_params.block_name
                << "' requires an explicit projection_level." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection block '" << out_params.block_name
              << "' has invalid projection_layout='" << layout_name_
              << "'. Valid choices are native_amr and uniform." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (pin->DoesParameterExist(op.block_name, "projection_axes") &&
      pin->DoesParameterExist(op.block_name, "projection_axis")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' sets both projection_axis and projection_axes. Use only one."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string axes = pin->DoesParameterExist(op.block_name, "projection_axes")
      ? pin->GetString(op.block_name, "projection_axes")
      : pin->GetOrAddString(op.block_name, "projection_axis", "x3");
  std::stringstream axes_stream(axes);
  std::string axis_token;
  while (std::getline(axes_stream, axis_token, ',')) {
    axis_token.erase(std::remove_if(axis_token.begin(), axis_token.end(),
        [](unsigned char c) { return std::isspace(c); }), axis_token.end());
    int axis = 0;
    if (axis_token == "x1") axis = 1;
    if (axis_token == "x2") axis = 2;
    if (axis_token == "x3") axis = 3;
    if (axis == 0 || projection_axes_[axis - 1]) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Projection output block '" << out_params.block_name
                << "' has invalid projection axes specification '" << axes
                << "'. Use one or two unique axes chosen from x1,x2,x3." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    projection_axes_[axis - 1] = true;
    if (!projection_axes_name_.empty()) projection_axes_name_ += ",";
    projection_axes_name_ += axis_token;
    n_projection_axes_++;
  }
  if (n_projection_axes_ < 1 || n_projection_axes_ > 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' must reduce over one or two Cartesian axes." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bool has_bounded_reduction = false;
  for (int axis = 1; axis <= 3; ++axis) {
    if (!projection_axes_[axis - 1]) image_axes_[n_image_axes_++] = axis;
  }

  int finest_level = pm->max_level - pm->root_level;
  if (layout_ == Layout::uniform) {
    projection_level_ = pin->GetInteger(op.block_name, "projection_level");
    if (projection_level_ < 0 || projection_level_ > finest_level) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Projection block '" << out_params.block_name
                << "' requests projection_level=" << projection_level_
                << " outside physical AMR range [0, " << finest_level << "]."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  int refine_factor = (layout_ == Layout::uniform) ? (1 << projection_level_) : 1;
  for (int d = 0; d < n_image_axes_; ++d) {
    image_min_[d] = AxisMinimum(pm->mesh_size, image_axes_[d]);
    image_max_[d] = AxisMaximum(pm->mesh_size, image_axes_[d]);
    image_nx_[d] = AxisCellCount(pm->mesh_indcs, image_axes_[d]) * refine_factor;
    image_dx_[d] = (image_max_[d] - image_min_[d]) / static_cast<Real>(image_nx_[d]);
  }
  if (n_image_axes_ == 1) {
    image_axes_[1] = 0;
    image_min_[1] = 0.0;
    image_max_[1] = 1.0;
    image_nx_[1] = 1;
    image_dx_[1] = 1.0;
  }

  if (n_projection_axes_ > 1 &&
      (pin->DoesParameterExist(op.block_name, "projection_min") ||
       pin->DoesParameterExist(op.block_name, "projection_max"))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' reduces over multiple axes; use projection_xN_min/max bounds."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for (int axis = 1; axis <= 3; ++axis) {
    if (!projection_axes_[axis - 1]) continue;
    Real full_min = AxisMinimum(pm->mesh_size, axis);
    Real full_max = AxisMaximum(pm->mesh_size, axis);
    std::string min_name = "projection_x" + std::to_string(axis) + "_min";
    std::string max_name = "projection_x" + std::to_string(axis) + "_max";
    if (n_projection_axes_ == 1 &&
        !pin->DoesParameterExist(op.block_name, min_name) &&
        !pin->DoesParameterExist(op.block_name, max_name) &&
        (pin->DoesParameterExist(op.block_name, "projection_min") ||
         pin->DoesParameterExist(op.block_name, "projection_max"))) {
      projection_min_[axis - 1] = pin->GetOrAddReal(op.block_name, "projection_min",
                                                     full_min);
      projection_max_[axis - 1] = pin->GetOrAddReal(op.block_name, "projection_max",
                                                     full_max);
    } else {
      projection_min_[axis - 1] = pin->GetOrAddReal(op.block_name, min_name, full_min);
      projection_max_[axis - 1] = pin->GetOrAddReal(op.block_name, max_name, full_max);
    }
    if (projection_min_[axis - 1] < full_min ||
        projection_max_[axis - 1] > full_max ||
        projection_min_[axis - 1] >= projection_max_[axis - 1]) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Projection bounds [" << projection_min_[axis - 1]
                << ", " << projection_max_[axis - 1] << "] in block '"
                << out_params.block_name << "' are invalid on x" << axis << "."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (projection_min_[axis - 1] > full_min ||
        projection_max_[axis - 1] < full_max) {
      has_bounded_reduction = true;
    }
  }
  if (layout_ == Layout::uniform && !has_bounded_reduction) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Uniform projection block '" << out_params.block_name
              << "' is a bounded quicklook/reference mode; at least one projected "
              << "axis must use a strict sub-domain interval." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  weighting_name_ = pin->GetOrAddString(op.block_name, "weighting", "integral");
  if (weighting_name_ == "integral") {
    weighting_ = Weighting::integral;
  } else if (weighting_name_ == "volume") {
    weighting_ = Weighting::volume;
  } else if (weighting_name_ == "mass") {
    weighting_ = Weighting::mass;
    if (pm->pmb_pack->phydro != nullptr) {
      density_ptr_ = &(pm->pmb_pack->phydro->u0);
    } else if (pm->pmb_pack->pmhd != nullptr) {
      density_ptr_ = &(pm->pmb_pack->pmhd->u0);
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Mass-weighted projection block '"
                << out_params.block_name
                << "' requires a Hydro or MHD density field." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' has invalid weighting='" << weighting_name_
              << "'. Valid choices are integral, volume, and mass." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  statistic_name_ = pin->GetOrAddString(op.block_name, "statistic", "value");
  if (statistic_name_ == "value") {
    statistic_ = Statistic::value;
  } else if (statistic_name_ == "stddev" && weighting_ != Weighting::integral) {
    statistic_ = Statistic::stddev;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Projection output block '" << out_params.block_name
              << "' must use statistic=value, or statistic=stddev with "
              << "weighting=volume or mass." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (layout_ == Layout::uniform) {
    uniform_max_megabytes_ = pin->GetOrAddInteger(op.block_name,
                                                   "projection_max_megabytes", 256);
    if (uniform_max_megabytes_ <= 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Uniform projection block '" << out_params.block_name
                << "' requires projection_max_megabytes > 0." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::uint64_t arrays = static_cast<std::uint64_t>(outvars.size());
    if (weighting_ != Weighting::integral) arrays += 1;
    if (statistic_ == Statistic::stddev) arrays += outvars.size();
    std::uint64_t bytes = arrays * static_cast<std::uint64_t>(image_nx_[0]) *
        static_cast<std::uint64_t>(image_nx_[1]) * sizeof(Real);
    std::uint64_t limit = static_cast<std::uint64_t>(uniform_max_megabytes_) *
        1024ULL * 1024ULL;
    if (bytes > limit) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Uniform projection block '" << out_params.block_name
                << "' needs " << bytes << " working bytes per rank, exceeding "
                << "projection_max_megabytes=" << uniform_max_megabytes_
                << ". Use projection_layout=native_amr or reduce projection_level."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Form additive patches after staging only cells in intersecting MeshBlocks.

void ProjectionOutput::LoadOutputData(Mesh *pm) {
  Kokkos::Timer output_timer;
  patches_.clear();
  shard_payload_.clear();
  local_patch_count_ = 0;
  shard_patch_count_ = 0;
  local_staged_cells_ = 0;
  shard_staged_cells_ = 0;
  local_payload_bytes_ = 0;
  communicated_bytes_ = 0;
  peak_output_bytes_ = 0;

  if (out_params.contains_derived) {
    out_params.i_derived = 0;
    ComputeDerivedVariable(out_params.variable, pm);
  }

  const int nvars = static_cast<int>(outvars.size());
  const auto &indcs = pm->mb_indcs;
  const int ncell[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int begin_active[3] = {indcs.is, indcs.js, indcs.ks};
  auto &size = pm->pmb_pack->pmb->mb_size;
  std::size_t accumulated_patch_bytes = 0;
  std::size_t dense_bytes = 0;
  std::size_t uniform_limit = 0;
  if (layout_ == Layout::uniform) {
    std::size_t arrays = static_cast<std::size_t>(nvars);
    if (statistic_ == Statistic::stddev) arrays += nvars;
    if (weighting_ != Weighting::integral) arrays += 1;
    dense_bytes = arrays * static_cast<std::size_t>(image_nx_[0]) *
        image_nx_[1] * sizeof(Real);
    uniform_limit = static_cast<std::size_t>(uniform_max_megabytes_) *
        1024ULL * 1024ULL;
  }

  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    int gid = pm->pmb_pack->pmb->mb_gid.h_view(m);
    OutputMeshBlockInfo mb(gid, indcs.is, indcs.ie, indcs.js, indcs.je,
                           indcs.ks, indcs.ke, size.h_view(m).x1min,
                           size.h_view(m).x1max, size.h_view(m).x2min,
                           size.h_view(m).x2max, size.h_view(m).x3min,
                           size.h_view(m).x3max);
    bool intersects = true;
    int first[3] = {indcs.is, indcs.js, indcs.ks};
    int last[3] = {indcs.ie, indcs.je, indcs.ke};
    for (int axis = 1; axis <= 3; ++axis) {
      if (!projection_axes_[axis - 1]) continue;
      if (Overlap(BlockMinimum(mb, axis), BlockMaximum(mb, axis),
                  projection_min_[axis - 1], projection_max_[axis - 1]) <= 0.0) {
        intersects = false;
        break;
      }
      first[axis - 1] = begin_active[axis - 1] + ncell[axis - 1];
      last[axis - 1] = begin_active[axis - 1] - 1;
      for (int c = 0; c < ncell[axis - 1]; ++c) {
        Real lower = CellLowerEdge(mb, axis, c, ncell[axis - 1]);
        Real upper = CellLowerEdge(mb, axis, c + 1, ncell[axis - 1]);
        if (Overlap(lower, upper, projection_min_[axis - 1],
                    projection_max_[axis - 1]) > 0.0) {
          first[axis - 1] = std::min(first[axis - 1], begin_active[axis - 1] + c);
          last[axis - 1] = std::max(last[axis - 1], begin_active[axis - 1] + c);
        }
      }
    }
    if (!intersects) continue;

    int ni = last[0] - first[0] + 1;
    int nj = last[1] - first[1] + 1;
    int nk = last[2] - first[2] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) continue;
    local_staged_cells_ += static_cast<std::size_t>(ni) * nj * nk;
    std::size_t npixels = static_cast<std::size_t>(
        ncell[image_axes_[0] - 1]) *
        ((n_image_axes_ == 2) ? ncell[image_axes_[1] - 1] : 1);
    std::size_t patch_arrays = 1 + static_cast<std::size_t>(nvars);
    if (statistic_ == Statistic::stddev) patch_arrays += nvars;
    std::size_t patch_storage_bytes = sizeof(Patch) +
        patch_arrays * npixels * sizeof(Real);
    std::size_t staging_arrays = static_cast<std::size_t>(nvars) +
        ((weighting_ == Weighting::mass) ? 1 : 0);
    std::size_t staging_bytes = static_cast<std::size_t>(ni) * nj * nk *
        staging_arrays * sizeof(Real);
    if (layout_ == Layout::uniform) {
      std::size_t required = std::max(
          accumulated_patch_bytes + staging_bytes + patch_storage_bytes,
          accumulated_patch_bytes + patch_storage_bytes + dense_bytes);
      if (required > uniform_limit) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Uniform projection block '"
                  << out_params.block_name << "' needs " << required
                  << " bytes including selected-cell staging, exceeding "
                  << "projection_max_megabytes=" << uniform_max_megabytes_
                  << ". Use projection_layout=native_amr." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    HostArray4D<Real> values("projection_values", nvars, nk, nj, ni);
    for (int n = 0; n < nvars; ++n) {
      auto d_slice = Kokkos::subview(*(outvars[n].data_ptr), m, outvars[n].data_index,
          std::make_pair(first[2], last[2] + 1),
          std::make_pair(first[1], last[1] + 1),
          std::make_pair(first[0], last[0] + 1));
      auto h_slice = Kokkos::subview(values, n, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      Kokkos::deep_copy(h_slice, d_slice);
    }
    HostArray3D<Real> density;
    if (weighting_ == Weighting::mass) {
      density = HostArray3D<Real>("projection_density", nk, nj, ni);
      auto d_slice = Kokkos::subview(*density_ptr_, m, static_cast<int>(IDN),
          std::make_pair(first[2], last[2] + 1),
          std::make_pair(first[1], last[1] + 1),
          std::make_pair(first[0], last[0] + 1));
      Kokkos::deep_copy(density, d_slice);
    }

    Patch patch;
    patch.gid = gid;
    patch.level = pm->lloc_eachmb[gid].level - pm->root_level;
    patch.nx = ncell[image_axes_[0] - 1];
    patch.ny = (n_image_axes_ == 2) ? ncell[image_axes_[1] - 1] : 1;
    patch.xmin = BlockMinimum(mb, image_axes_[0]);
    patch.xmax = BlockMaximum(mb, image_axes_[0]);
    patch.ymin = (n_image_axes_ == 2) ? BlockMinimum(mb, image_axes_[1]) : 0.0;
    patch.ymax = (n_image_axes_ == 2) ? BlockMaximum(mb, image_axes_[1]) : 1.0;
    patch.weight.assign(npixels, 0.0);
    patch.first_moment.assign(static_cast<std::size_t>(nvars) * npixels, 0.0);
    if (statistic_ == Statistic::stddev) {
      patch.second_moment.assign(static_cast<std::size_t>(nvars) * npixels, 0.0);
    }

    for (int kk = 0; kk < nk; ++kk) {
      int actual_k = first[2] + kk;
      for (int jj = 0; jj < nj; ++jj) {
        int actual_j = first[1] + jj;
        for (int ii = 0; ii < ni; ++ii) {
          int actual_i = first[0] + ii;
          int actual[3] = {actual_i, actual_j, actual_k};
          Real reduced_measure = 1.0;
          for (int axis = 1; axis <= 3; ++axis) {
            if (!projection_axes_[axis - 1]) continue;
            int c = actual[axis - 1] - begin_active[axis - 1];
            reduced_measure *= Overlap(
                CellLowerEdge(mb, axis, c, ncell[axis - 1]),
                CellLowerEdge(mb, axis, c + 1, ncell[axis - 1]),
                projection_min_[axis - 1], projection_max_[axis - 1]);
          }
          if (reduced_measure <= 0.0) continue;
          int retained0 = actual[image_axes_[0] - 1] - begin_active[image_axes_[0] - 1];
          int retained1 = (n_image_axes_ == 2)
              ? actual[image_axes_[1] - 1] - begin_active[image_axes_[1] - 1] : 0;
          std::size_t p = static_cast<std::size_t>(retained1) * patch.nx + retained0;
          Real factor = reduced_measure;
          if (weighting_ == Weighting::mass) factor *= density(kk, jj, ii);
          patch.weight[p] += factor;
          for (int n = 0; n < nvars; ++n) {
            Real value = values(n, kk, jj, ii);
            std::size_t q = static_cast<std::size_t>(n) * npixels + p;
            patch.first_moment[q] += factor * value;
            if (statistic_ == Statistic::stddev) {
              patch.second_moment[q] += factor * value * value;
            }
          }
        }
      }
    }
    if (layout_ == Layout::native_amr) {
      std::int32_t stored_gid = patch.gid;
      std::int32_t stored_level = patch.level;
      std::int32_t stored_nx = patch.nx;
      std::int32_t stored_ny = patch.ny;
      AppendBytes(shard_payload_, stored_gid);
      AppendBytes(shard_payload_, stored_level);
      AppendBytes(shard_payload_, stored_nx);
      AppendBytes(shard_payload_, stored_ny);
      AppendBytes(shard_payload_, static_cast<double>(patch.xmin));
      AppendBytes(shard_payload_, static_cast<double>(patch.xmax));
      AppendBytes(shard_payload_, static_cast<double>(patch.ymin));
      AppendBytes(shard_payload_, static_cast<double>(patch.ymax));
      AppendMomentValues(shard_payload_, patch.weight);
      AppendMomentValues(shard_payload_, patch.first_moment);
      if (statistic_ == Statistic::stddev) {
        AppendMomentValues(shard_payload_, patch.second_moment);
      }
      local_patch_count_++;
      shard_patch_count_ = local_patch_count_;
      local_payload_bytes_ = shard_payload_.size();
      peak_output_bytes_ = std::max(
          peak_output_bytes_,
          shard_payload_.capacity() + patch_storage_bytes + staging_bytes);
    } else {
      patches_.push_back(std::move(patch));
      local_patch_count_++;
      accumulated_patch_bytes += patch_storage_bytes;
      peak_output_bytes_ = std::max(
          peak_output_bytes_,
          std::max(accumulated_patch_bytes + staging_bytes,
                   accumulated_patch_bytes + dense_bytes));
    }
  }

  if (layout_ == Layout::uniform) {
    int nx = image_nx_[0];
    int ny = image_nx_[1];
    Kokkos::realloc(projected_data_, nvars, ny, nx);
    Kokkos::deep_copy(projected_data_, 0.0);
    if (statistic_ == Statistic::stddev) {
      Kokkos::realloc(second_moment_, nvars, ny, nx);
      Kokkos::deep_copy(second_moment_, 0.0);
    }
    if (weighting_ != Weighting::integral) {
      Kokkos::realloc(normalization_, ny, nx);
      Kokkos::deep_copy(normalization_, 0.0);
    }

    for (const auto &patch : patches_) {
      Real pdx = (patch.xmax - patch.xmin) / static_cast<Real>(patch.nx);
      Real pdy = (patch.ymax - patch.ymin) / static_cast<Real>(patch.ny);
      std::size_t npixels = static_cast<std::size_t>(patch.nx) * patch.ny;
      for (int py = 0; py < patch.ny; ++py) {
        Real plo1 = patch.ymin + static_cast<Real>(py) * pdy;
        Real phi1 = plo1 + pdy;
        int y0 = (n_image_axes_ == 2) ? std::max(0, static_cast<int>(
            std::floor((plo1 - image_min_[1]) / image_dx_[1] + 1.0e-10))) : 0;
        int y1 = (n_image_axes_ == 2) ? std::min(ny - 1, static_cast<int>(
            std::ceil((phi1 - image_min_[1]) / image_dx_[1] - 1.0e-10)) - 1) : 0;
        for (int px = 0; px < patch.nx; ++px) {
          Real plo0 = patch.xmin + static_cast<Real>(px) * pdx;
          Real phi0 = plo0 + pdx;
          int x0 = std::max(0, static_cast<int>(
              std::floor((plo0 - image_min_[0]) / image_dx_[0] + 1.0e-10)));
          int x1 = std::min(nx - 1, static_cast<int>(
              std::ceil((phi0 - image_min_[0]) / image_dx_[0] - 1.0e-10)) - 1);
          std::size_t p = static_cast<std::size_t>(py) * patch.nx + px;
          for (int iy = y0; iy <= y1; ++iy) {
            Real fraction = (n_image_axes_ == 2)
                ? Overlap(plo1, phi1, image_min_[1] + iy * image_dx_[1],
                          image_min_[1] + (iy + 1) * image_dx_[1]) / image_dx_[1]
                : 1.0;
            for (int ix = x0; ix <= x1; ++ix) {
              fraction *= Overlap(plo0, phi0, image_min_[0] + ix * image_dx_[0],
                                  image_min_[0] + (ix + 1) * image_dx_[0]) /
                                  image_dx_[0];
              for (int n = 0; n < nvars; ++n) {
                std::size_t q = static_cast<std::size_t>(n) * npixels + p;
                projected_data_(n, iy, ix) += patch.first_moment[q] * fraction;
                if (statistic_ == Statistic::stddev) {
                  second_moment_(n, iy, ix) += patch.second_moment[q] * fraction;
                }
              }
              if (weighting_ != Weighting::integral) {
                normalization_(iy, ix) += patch.weight[p] * fraction;
              }
              fraction = (n_image_axes_ == 2)
                  ? Overlap(plo1, phi1, image_min_[1] + iy * image_dx_[1],
                            image_min_[1] + (iy + 1) * image_dx_[1]) / image_dx_[1]
                  : 1.0;
            }
          }
        }
      }
    }

#if MPI_PARALLEL_ENABLED
    int count = nvars * nx * ny;
    if (global_variable::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, projected_data_.data(), count, MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(projected_data_.data(), projected_data_.data(), count, MPI_ATHENA_REAL,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (statistic_ == Statistic::stddev) {
      if (global_variable::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, second_moment_.data(), count, MPI_ATHENA_REAL, MPI_SUM,
                   0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(second_moment_.data(), second_moment_.data(), count, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      }
    }
    if (weighting_ != Weighting::integral) {
      if (global_variable::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, normalization_.data(), nx * ny, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(normalization_.data(), normalization_.data(), nx * ny,
                   MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
      }
    }
    std::size_t communicated_arrays = static_cast<std::size_t>(nvars);
    if (statistic_ == Statistic::stddev) communicated_arrays += nvars;
    if (weighting_ != Weighting::integral) communicated_arrays += 1;
    communicated_bytes_ = communicated_arrays * nx * ny * sizeof(Real) *
        static_cast<std::size_t>(std::max(global_variable::nranks - 1, 0));
#endif
    if (global_variable::my_rank == 0 && weighting_ != Weighting::integral) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
          Real norm = normalization_(iy, ix);
          for (int n = 0; n < nvars; ++n) {
            if (norm <= 0.0) {
              projected_data_(n, iy, ix) = std::numeric_limits<Real>::quiet_NaN();
              continue;
            }
            Real mean = projected_data_(n, iy, ix) / norm;
            if (statistic_ == Statistic::stddev) {
              Real variance = second_moment_(n, iy, ix) / norm - mean * mean;
              projected_data_(n, iy, ix) =
                  std::sqrt(std::max(variance, static_cast<Real>(0.0)));
            } else {
              projected_data_(n, iy, ix) = mean;
            }
          }
        }
      }
    }
    std::size_t arrays = static_cast<std::size_t>(nvars);
    if (statistic_ == Statistic::stddev) arrays += nvars;
    if (weighting_ != Weighting::integral) arrays += 1;
    peak_output_bytes_ = std::max(peak_output_bytes_, accumulated_patch_bytes +
                                  arrays * static_cast<std::size_t>(nx) * ny *
                                  sizeof(Real));
    patches_.clear();
    patches_.shrink_to_fit();
  } else {
    shard_staged_cells_ = local_staged_cells_;

#if MPI_PARALLEL_ENABLED
    if (out_params.file_shard_mode != FileShardMode::per_rank) {
      MPI_Comm comm = (out_params.file_shard_mode == FileShardMode::per_node)
          ? global_variable::node_comm : MPI_COMM_WORLD;
      int local_index = (out_params.file_shard_mode == FileShardMode::per_node)
          ? global_variable::rank_in_node : global_variable::my_rank;
      int participants = (out_params.file_shard_mode == FileShardMode::per_node)
          ? global_variable::ranks_per_node : global_variable::nranks;
      if (local_payload_bytes_ >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Projection shard payload exceeds MPI_Gatherv int count; "
                  << "use single_file_per_rank=true." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      int local_bytes = static_cast<int>(local_payload_bytes_);
      int local_patches = shard_patch_count_;
      std::uint64_t local_staged = static_cast<std::uint64_t>(local_staged_cells_);
      std::uint64_t gathered_staged = local_staged;
      std::vector<int> counts;
      std::vector<int> patch_counts;
      if (local_index == 0) {
        counts.resize(participants, 0);
        patch_counts.resize(participants, 0);
      }
      MPI_Gather(&local_bytes, 1, MPI_INT, counts.empty() ? nullptr : counts.data(),
                 1, MPI_INT, 0, comm);
      MPI_Gather(&local_patches, 1, MPI_INT,
                 patch_counts.empty() ? nullptr : patch_counts.data(), 1, MPI_INT, 0,
                 comm);
      MPI_Reduce(&local_staged, &gathered_staged, 1, MPI_UINT64_T, MPI_SUM, 0, comm);
      std::vector<int> displs;
      std::vector<char> gathered;
      if (local_index == 0) {
        displs.resize(participants, 0);
        std::size_t total = 0;
        shard_patch_count_ = 0;
        for (int r = 0; r < participants; ++r) {
          if (total > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "Projection gathered shard exceeds MPI count; use "
                      << "single_file_per_rank=true." << std::endl;
            std::exit(EXIT_FAILURE);
          }
          displs[r] = static_cast<int>(total);
          total += counts[r];
          shard_patch_count_ += patch_counts[r];
        }
        gathered.resize(total);
        shard_staged_cells_ = static_cast<std::size_t>(gathered_staged);
        communicated_bytes_ = total - local_payload_bytes_;
        peak_output_bytes_ = std::max(
            peak_output_bytes_, shard_payload_.capacity() + gathered.capacity());
      }
      MPI_Gatherv(shard_payload_.data(), local_bytes, MPI_CHAR,
                  gathered.empty() ? nullptr : gathered.data(),
                  counts.empty() ? nullptr : counts.data(),
                  displs.empty() ? nullptr : displs.data(), MPI_CHAR, 0, comm);
      if (local_index == 0) {
        shard_payload_.swap(gathered);
      } else {
        shard_payload_.clear();
        shard_patch_count_ = 0;
      }
    }
#endif
  }

  output_elapsed_seconds_ = output_timer.seconds();
}

//----------------------------------------------------------------------------------------
//! \brief Write binary native-AMR patches or the guarded uniform ASCII quicklook.

void ProjectionOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  Kokkos::Timer write_timer;
  bool i_write = (layout_ == Layout::uniform)
      ? (global_variable::my_rank == 0) : IsShardWriter(out_params.file_shard_mode);
  std::string path_prefix = "proj/";
  if (layout_ == Layout::native_amr &&
      out_params.file_shard_mode != FileShardMode::shared) {
    path_prefix += ShardDirectoryName(out_params.file_shard_mode);
  }
  char number[7];
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  std::string fname = path_prefix + out_params.file_basename + "." +
      out_params.file_id + number + ".proj";

  if (i_write && layout_ == Layout::native_amr) {
    std::FILE *pfile = std::fopen(fname.c_str(), "wb");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::stringstream params;
    pin->ParameterDump(params);
    std::string pbuf = params.str();
    std::stringstream hdr;
    hdr << "AthenaK projection output version=3.0" << std::endl
        << "  layout=native_amr" << std::endl
        << "  distribution=" << ShardDistributionName(out_params.file_shard_mode)
        << std::endl
        << "  writer=" << ShardWriterId(out_params.file_shard_mode) << std::endl
        << "  time=" << pm->time << std::endl
        << "  cycle=" << pm->ncycle << std::endl
        << "  projection_axes=" << projection_axes_name_ << std::endl
        << "  weighting=" << weighting_name_ << std::endl
        << "  statistic=" << statistic_name_ << std::endl
        << "  image_dimension=" << n_image_axes_ << std::endl
        << "  image_axis0=x" << image_axes_[0] << std::endl
        << "  image_axis1=" << (n_image_axes_ == 2
            ? ("x" + std::to_string(image_axes_[1])) : "none") << std::endl
        << "  domain_xmin=" << image_min_[0] << std::endl
        << "  domain_xmax=" << image_max_[0] << std::endl
        << "  domain_ymin=" << image_min_[1] << std::endl
        << "  domain_ymax=" << image_max_[1] << std::endl
        << "  root_nx=" << AxisCellCount(pm->mesh_indcs, image_axes_[0]) << std::endl
        << "  root_ny=" << (n_image_axes_ == 2
            ? AxisCellCount(pm->mesh_indcs, image_axes_[1]) : 1) << std::endl
        << "  number_of_variables=" << outvars.size() << std::endl
        << "  number_of_moments=" << (statistic_ == Statistic::stddev ? 3 : 2)
        << std::endl
        << "  npatches=" << shard_patch_count_ << std::endl
        << "  size_of_moment=" << sizeof(double) << std::endl;
    for (int axis = 1; axis <= 3; ++axis) {
      if (projection_axes_[axis - 1]) {
        hdr << "  projection_x" << axis << "_min=" << projection_min_[axis - 1]
            << std::endl
            << "  projection_x" << axis << "_max=" << projection_max_[axis - 1]
            << std::endl;
      }
    }
    hdr << "  variables:";
    for (const auto &var : outvars) hdr << " " << var.label;
    hdr << std::endl << "  header_offset=" << pbuf.size() << std::endl;
    std::string hbuf = hdr.str();
    CheckedFwrite(hbuf.data(), 1, hbuf.size(), pfile, fname);
    CheckedFwrite(pbuf.data(), 1, pbuf.size(), pfile, fname);
    CheckedFwrite(shard_payload_.data(), 1, shard_payload_.size(), pfile, fname);
    if (std::fclose(pfile) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed closing projection output '" << fname << "'"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (i_write && layout_ == Layout::uniform) {
    std::FILE *pfile = std::fopen(fname.c_str(), "w");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fprintf(pfile, "# AthenaK projection output version=3.0\n");
    std::fprintf(pfile, "# layout=uniform distribution=shared time=%.17e cycle=%d\n",
                 static_cast<double>(pm->time), pm->ncycle);
    std::fprintf(pfile, "# projection_axes=%s weighting=%s statistic=%s "
                 "projection_level=%d\n", projection_axes_name_.c_str(),
                 weighting_name_.c_str(), statistic_name_.c_str(), projection_level_);
    for (int axis = 1; axis <= 3; ++axis) {
      if (projection_axes_[axis - 1]) {
        std::fprintf(pfile, "# projection_x%d_min=%.17e projection_x%d_max=%.17e\n",
                     axis, static_cast<double>(projection_min_[axis - 1]), axis,
                     static_cast<double>(projection_max_[axis - 1]));
      }
    }
    std::fprintf(pfile, "# image_dimension=%d image_axis0=x%d image_axis1=%s nx=%d "
                 "ny=%d xmin=%.17e xmax=%.17e ymin=%.17e ymax=%.17e\n",
                 n_image_axes_, image_axes_[0],
                 n_image_axes_ == 2 ? ("x" + std::to_string(image_axes_[1])).c_str()
                                    : "none",
                 image_nx_[0], image_nx_[1], static_cast<double>(image_min_[0]),
                 static_cast<double>(image_max_[0]), static_cast<double>(image_min_[1]),
                 static_cast<double>(image_max_[1]));
    std::fprintf(pfile, "# columns: i j x%dv %s", image_axes_[0],
                 n_image_axes_ == 2 ? ("x" + std::to_string(image_axes_[1]) + "v").c_str()
                                    : "dummy");
    for (const auto &var : outvars) std::fprintf(pfile, " %s", var.label.c_str());
    std::fprintf(pfile, "\n");
    for (int iy = 0; iy < image_nx_[1]; ++iy) {
      Real y = image_min_[1] + (static_cast<Real>(iy) + 0.5) * image_dx_[1];
      for (int ix = 0; ix < image_nx_[0]; ++ix) {
        Real x = image_min_[0] + (static_cast<Real>(ix) + 0.5) * image_dx_[0];
        std::fprintf(pfile, "%d %d %.17e %.17e", ix, iy, static_cast<double>(x),
                     static_cast<double>(y));
        for (int n = 0; n < static_cast<int>(outvars.size()); ++n) {
          std::fprintf(pfile, " %.17e", static_cast<double>(projected_data_(n, iy, ix)));
        }
        std::fprintf(pfile, "\n");
      }
    }
    if (std::fclose(pfile) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Failed closing projection output '" << fname << "'"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  output_elapsed_seconds_ += write_timer.seconds();
  if (i_write && ProjectionStatsEnabled()) {
    std::cout << "[output-io] type=proj"
              << " id=" << out_params.file_id
              << " layout=" << layout_name_
              << " mode=" << ShardDistributionName(out_params.file_shard_mode)
              << " writer=" << ShardWriterId(out_params.file_shard_mode)
              << " local_patches=" << local_patch_count_
              << " shard_patches=" << shard_patch_count_
              << " staged_cells=" << local_staged_cells_
              << " shard_staged_cells=" << shard_staged_cells_
              << " payload_bytes=" << local_payload_bytes_
              << " shard_payload_bytes=" << shard_payload_.size()
              << " communicated_bytes=" << communicated_bytes_
              << " peak_output_bytes=" << peak_output_bytes_
              << " elapsed=" << output_elapsed_seconds_ << std::endl;
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
