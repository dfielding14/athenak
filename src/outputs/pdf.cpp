//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file pdf.cpp
//! \brief writes versioned N-dimensional PDF output data

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "diagnostic_semantics.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"
#include "parameter_input.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

[[noreturn]] void FatalPDFError(const std::string &message);

std::size_t PDFCountAsSize(int value, const char *context) {
  if (value < 0) {
    FatalPDFError(std::string(context) + " is negative.");
  }
  return static_cast<std::size_t>(value);
}

[[noreturn]] void FatalPDFError(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" +
                        message);
}

std::size_t PDFPersistentAllocationBytes(const OutputParameters &op) {
  std::size_t edge_elements = 0;
  std::size_t max_edge_elements = 0;
  std::size_t total_bins = 1;
  for (int d = 0; d < op.pdf_ndim; ++d) {
    std::size_t bins = PDFCountAsSize(op.pdf_nbin[d], "PDF bin count");
    edge_elements = output_file_utils::CheckedSizeAdd(
        edge_elements, bins + 1, "PDF bin edges", FatalPDFError);
    max_edge_elements = std::max(max_edge_elements, bins + 1);
    total_bins = output_file_utils::CheckedSizeProduct(
        total_bins, bins + 2, "PDF histogram", FatalPDFError);
  }
  std::size_t edges = output_file_utils::CheckedSizeProduct(
      edge_elements, sizeof(Real), "PDF device bin edges", FatalPDFError);
  std::size_t edge_mirror = output_file_utils::CheckedSizeProduct(
      max_edge_elements, sizeof(Real), "PDF host bin-edge mirror", FatalPDFError);
  std::size_t histogram = output_file_utils::CheckedSizeProduct(
      total_bins, sizeof(Real), "PDF histogram result storage",
      FatalPDFError);
  return output_file_utils::CheckedSizeAdd(
      output_file_utils::CheckedSizeAdd(edges, edge_mirror,
                                        "PDF persistent allocation", FatalPDFError),
      histogram, "PDF persistent allocation", FatalPDFError);
}

std::size_t PDFResultMirrorBytes(int total_bins) {
  return output_file_utils::CheckedSizeProduct(
      PDFCountAsSize(total_bins, "PDF total bin count"), sizeof(Real),
      "PDF host result mirror", FatalPDFError);
}

void RequirePDFBudget(std::size_t bytes, std::size_t limit, const char *context) {
  output_file_utils::RequireAllocationBudget(bytes, limit, context, FatalPDFError);
}

void ValidatePDFResult(const DvceArray1D<Real> &result, int total_bins,
                       const char *context) {
  auto host_result = Kokkos::create_mirror_view(result);
  Kokkos::deep_copy(host_result, result);
  Kokkos::fence();
  for (int n = 0; n < total_bins; ++n) {
    if (!output_diagnostics::IsFinite(host_result(n))) {
      FatalPDFError(std::string(context) + " contains a nonfinite histogram bin.");
    }
  }
}

#if MPI_PARALLEL_ENABLED
void ReducePDFResultViaHost(DvceArray1D<Real> &result, int total_bins, MPI_Comm comm,
                            int root, const char *context) {
  auto host_result = Kokkos::create_mirror_view(result);
  Kokkos::deep_copy(host_result, result);
  Kokkos::fence();
  int comm_rank = 0;
  mpi_utils::CheckMpi(MPI_Comm_rank(comm, &comm_rank),
                      "MPI_Comm_rank for host-staged PDF reduction");
  if (comm_rank == root) {
    mpi_utils::CheckMpi(MPI_Reduce(MPI_IN_PLACE, host_result.data(), total_bins,
                                   MPI_ATHENA_REAL, MPI_SUM, root, comm), context);
    Kokkos::deep_copy(result, host_result);
    Kokkos::fence();
  } else {
    mpi_utils::CheckMpi(MPI_Reduce(host_result.data(), host_result.data(), total_bins,
                                   MPI_ATHENA_REAL, MPI_SUM, root, comm), context);
  }
}
#endif

void DiscardTemporaryPDFFile(std::FILE *output, const std::string &filename) {
  static_cast<void>(std::fclose(output));
  output_file_utils::DiscardOwnedPath(filename);
}

void CheckedPDFPrint(std::FILE *output, const std::string &filename, const char *text) {
  if (std::fputs(text, output) == EOF) {
    DiscardTemporaryPDFFile(output, filename);
    FatalPDFError("Could not write PDF header '" + filename + "'.");
  }
}

template <typename Arg, typename... Args>
void CheckedPDFPrint(std::FILE *output, const std::string &filename,
                     const char *format, Arg arg, Args... args) {
  if (std::fprintf(output, format, arg, args...) < 0) {
    DiscardTemporaryPDFFile(output, filename);
    FatalPDFError("Could not write PDF header '" + filename + "'.");
  }
}

void CheckedPDFWrite(std::FILE *output, const void *data, std::size_t element_size,
                     std::size_t count, const std::string &filename,
                     const char *context) {
  if (count == 0) {
    return;
  }
  if (std::fwrite(data, element_size, count, output) != count) {
    DiscardTemporaryPDFFile(output, filename);
    FatalPDFError(std::string(context) + " was not written completely to '" +
                  filename + "'.");
  }
}

void CheckedLegacyPDFPrint(std::FILE *output, const std::string &filename,
                           const char *text) {
  if (std::fputs(text, output) == EOF) {
    FatalPDFError("Could not append legacy PDF output '" + filename + "'.");
  }
}

template <typename Arg, typename... Args>
void CheckedLegacyPDFPrint(std::FILE *output, const std::string &filename,
                           const char *format, Arg arg, Args... args) {
  if (std::fprintf(output, format, arg, args...) < 0) {
    FatalPDFError("Could not append legacy PDF output '" + filename + "'.");
  }
}

void CheckedLegacyPDFClose(std::FILE *output, const std::string &filename) {
  if (std::fclose(output) != 0) {
    FatalPDFError("Could not close legacy PDF output '" + filename + "'.");
  }
}

void PublishTemporaryPDFFile(std::FILE *output, const std::string &temporary_filename,
                             const std::string &filename, const char *context) {
  if (std::fclose(output) != 0) {
    output_file_utils::DiscardOwnedPath(temporary_filename);
    FatalPDFError("Could not close " + std::string(context) + " '" +
                  temporary_filename + "'.");
  }
  output_file_utils::PublishTemporaryFile(temporary_filename, filename, context,
                                          FatalPDFError);
}

std::string PDFDirectory(const OutputParameters &op) {
  std::string directory = "pdf_" + op.file_id;
  for (int d = 1; d < op.pdf_ndim; ++d) {
    directory += "_" + op.pdf_variables[d];
  }
  return directory;
}

void AdvanceOutputCounters(OutputParameters &op, Mesh *pm, ParameterInput *pin) {
  op.file_number = output_file_utils::AdvanceFileNumber(
      op.file_number, "PDF output", FatalPDFError);
  if (op.last_time < 0.0) {
    op.last_time = pm->time;
  } else {
    op.last_time += op.dt;
  }
  pin->SetInteger(op.block_name, "file_number", op.file_number);
  pin->SetReal(op.block_name, "last_time", op.last_time);
}

}  // namespace

//----------------------------------------------------------------------------------------
// Constructor

PDFOutput::PDFOutput(ParameterInput *pin, Mesh *pm, OutputParameters op)
    : BaseTypeOutput(pin, pm, op) {
  int configured_limit = pin->GetOrAddInteger(
      op.block_name, "max_writer_allocation_bytes",
      static_cast<int>(output_file_utils::kDefaultMaxWriterAllocationBytes));
  if (configured_limit <= 0) {
    FatalPDFError("PDF max_writer_allocation_bytes must be positive.");
  }
  max_writer_allocation_bytes = static_cast<std::size_t>(configured_limit);
  persistent_writer_allocation_bytes = PDFPersistentAllocationBytes(op);
  RequirePDFBudget(persistent_writer_allocation_bytes, max_writer_allocation_bytes,
                   "PDF persistent allocation");
  if (op.include_gzs) {
    FatalPDFError("PDF output block '" + op.block_name +
                  "' cannot set ghost_zones=true; PDFs sample active zones only.");
  }
  if (op.pdf_weight == "mass" && pm->pmb_pack->pionn != nullptr) {
    FatalPDFError("Mass-weighted PDF output block '" + op.block_name +
                  "' is ambiguous for <ion-neutral> two-fluid runs.");
  }
  std::string directory = PDFDirectory(op);
  output_file_utils::EnsureDirectory(directory, 0775, "PDF output", FatalPDFError);
  if (IsSharded(op.shard_mode)) {
    directory += "/" + ShardDirectoryName(op.shard_mode, global_variable::my_rank,
                                           global_variable::node_id);
    output_file_utils::EnsureDirectory(directory, 0775, "PDF output", FatalPDFError);
  }

  pdf_data.Initialize(op.pdf_ndim, op.pdf_nbin, op.pdf_bin_min, op.pdf_bin_max,
                      op.pdf_scale, op.pdf_linthresh);
  pdf_data.PopulateBinEdges();

  int expected_vars = op.pdf_ndim + (op.pdf_weight == "variable" ? 1 : 0);
  if (outvars.size() != static_cast<std::size_t>(expected_vars)) {
    FatalPDFError("PDF output block '" + op.block_name +
                  "' requires one scalar output field per axis" +
                  (op.pdf_weight == "variable" ? " and for its variable weight" : "") +
                  ".");
  }
}

//----------------------------------------------------------------------------------------
//! \brief Computes an N-dimensional histogram over active zones.

void PDFOutput::LoadOutputData(Mesh *pm) {
  int weight_mode = 0;  // 0=volume, 1=mass, 2=cell variable times volume
  if (out_params.pdf_weight == "mass") weight_mode = 1;
  if (out_params.pdf_weight == "variable") weight_mode = 2;

  DvceArray5D<Real> density_data;
  if (weight_mode == 1) {
    if (pm->pmb_pack->phydro != nullptr) {
      density_data = pm->pmb_pack->phydro->u0;
    } else if (pm->pmb_pack->pmhd != nullptr) {
      density_data = pm->pmb_pack->pmhd->u0;
    } else {
      FatalPDFError("Mass-weighted PDF requires Hydro or MHD density.");
    }
  }

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  int nmb = pm->pmb_pack->nmb_thispack;
  int nx1 = indcs.nx1 + 2*indcs.ng;
  int nx2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
  int nx3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;

  std::size_t field_elements = output_file_utils::CheckedSizeProduct(
      outvars.size(), PDFCountAsSize(nmb, "PDF MeshBlock count"),
      "PDF copied fields", FatalPDFError);
  field_elements = output_file_utils::CheckedSizeProduct(
      field_elements, PDFCountAsSize(nx1, "PDF x1 extent"), "PDF copied fields",
      FatalPDFError);
  field_elements = output_file_utils::CheckedSizeProduct(
      field_elements, PDFCountAsSize(nx2, "PDF x2 extent"), "PDF copied fields",
      FatalPDFError);
  field_elements = output_file_utils::CheckedSizeProduct(
      field_elements, PDFCountAsSize(nx3, "PDF x3 extent"), "PDF copied fields",
      FatalPDFError);
  std::size_t field_bytes = output_file_utils::CheckedSizeProduct(
      field_elements, sizeof(Real), "PDF copied fields", FatalPDFError);
  std::size_t cell_elements = output_file_utils::CheckedSizeProduct(
      PDFCountAsSize(nmb, "PDF MeshBlock count"), PDFCountAsSize(nx1, "PDF x1 extent"),
      "PDF derived fields", FatalPDFError);
  cell_elements = output_file_utils::CheckedSizeProduct(
      cell_elements, PDFCountAsSize(nx2, "PDF x2 extent"), "PDF derived fields",
      FatalPDFError);
  cell_elements = output_file_utils::CheckedSizeProduct(
      cell_elements, PDFCountAsSize(nx3, "PDF x3 extent"), "PDF derived fields",
      FatalPDFError);
  std::size_t derived_bytes = output_file_utils::CheckedSizeProduct(
      output_file_utils::CheckedSizeProduct(
          cell_elements, PDFCountAsSize(out_params.n_derived, "PDF derived-field count"),
          "PDF derived fields", FatalPDFError),
      sizeof(Real), "PDF derived fields", FatalPDFError);
  std::size_t result_mirror_bytes = PDFResultMirrorBytes(pdf_data.total_bins);
  constexpr std::size_t metadata_bytes =
      2*PDFData::MAX_DIM*(3*sizeof(int) + 5*sizeof(Real)) + 2*sizeof(int);
  std::size_t load_bytes = output_file_utils::CheckedSizeAdd(
      persistent_writer_allocation_bytes, field_bytes, "PDF load allocation",
      FatalPDFError);
  load_bytes = output_file_utils::CheckedSizeAdd(
      load_bytes, derived_bytes, "PDF load allocation", FatalPDFError);
  load_bytes = output_file_utils::CheckedSizeAdd(
      load_bytes, result_mirror_bytes, "PDF load allocation", FatalPDFError);
  load_bytes = output_file_utils::CheckedSizeAdd(
      load_bytes, metadata_bytes, "PDF load allocation", FatalPDFError);
  RequirePDFBudget(load_bytes, max_writer_allocation_bytes, "PDF load allocation");

  if (out_params.contains_derived) {
    out_params.i_derived = 0;
    for (int d = 0; d < out_params.pdf_ndim; ++d) {
      ComputeDerivedVariable(out_params.pdf_variables[d], pm);
    }
    if (out_params.pdf_weight == "variable") {
      ComputeDerivedVariable(out_params.pdf_weight_variable, pm);
    }
  }

  DvceArray5D<Real> fields("pdf_fields", outvars.size(), nmb, nx3, nx2, nx1);
  for (std::size_t n = 0; n < outvars.size(); ++n) {
    auto source = Kokkos::subview(*(outvars[n].data_ptr), Kokkos::make_pair(0, nmb),
        outvars[n].data_index, Kokkos::make_pair(0, nx3), Kokkos::make_pair(0, nx2),
        Kokkos::make_pair(0, nx1));
    auto target = Kokkos::subview(fields, n, Kokkos::ALL(), Kokkos::ALL(),
                                  Kokkos::ALL(), Kokkos::ALL());
    Kokkos::deep_copy(target, source);
  }
  Kokkos::fence();

  Kokkos::View<int[PDFData::MAX_DIM]> d_nbin("pdf_nbin");
  Kokkos::View<int[PDFData::MAX_DIM]> d_stride("pdf_stride");
  Kokkos::View<int[PDFData::MAX_DIM]> d_scale("pdf_scale");
  Kokkos::View<Real[PDFData::MAX_DIM]> d_step("pdf_step");
  Kokkos::View<Real[PDFData::MAX_DIM]> d_min("pdf_min");
  Kokkos::View<Real[PDFData::MAX_DIM]> d_max("pdf_max");
  Kokkos::View<Real[PDFData::MAX_DIM]> d_transform_min("pdf_transform_min");
  Kokkos::View<Real[PDFData::MAX_DIM]> d_linthresh("pdf_linthresh");
  auto h_nbin = Kokkos::create_mirror_view(d_nbin);
  auto h_stride = Kokkos::create_mirror_view(d_stride);
  auto h_scale = Kokkos::create_mirror_view(d_scale);
  auto h_step = Kokkos::create_mirror_view(d_step);
  auto h_min = Kokkos::create_mirror_view(d_min);
  auto h_max = Kokkos::create_mirror_view(d_max);
  auto h_transform_min = Kokkos::create_mirror_view(d_transform_min);
  auto h_linthresh = Kokkos::create_mirror_view(d_linthresh);
  for (int d = 0; d < PDFData::MAX_DIM; ++d) {
    h_nbin(d) = pdf_data.nbin[d];
    h_stride(d) = pdf_data.stride[d];
    h_scale(d) = pdf_data.scale[d];
    h_step(d) = pdf_data.step_size[d];
    h_min(d) = pdf_data.bin_min[d];
    h_max(d) = pdf_data.bin_max[d];
    h_transform_min(d) = pdf_data.transformed_min[d];
    h_linthresh(d) = pdf_data.linthresh[d];
  }
  Kokkos::deep_copy(d_nbin, h_nbin);
  Kokkos::deep_copy(d_stride, h_stride);
  Kokkos::deep_copy(d_scale, h_scale);
  Kokkos::deep_copy(d_step, h_step);
  Kokkos::deep_copy(d_min, h_min);
  Kokkos::deep_copy(d_max, h_max);
  Kokkos::deep_copy(d_transform_min, h_transform_min);
  Kokkos::deep_copy(d_linthresh, h_linthresh);
  Kokkos::fence();

  auto result = pdf_data.result_;
  Kokkos::deep_copy(result, 0.0);
  Kokkos::fence();
  int ndim = pdf_data.ndim;
  int weight_index = (weight_mode == 2) ? ndim : -1;
  DvceArray1D<int> invalid("pdf_invalid_sample", 1);
  Kokkos::deep_copy(invalid, 0);

  par_for("pdf_nd", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    int flat_index = 0;
    for (int d = 0; d < ndim; ++d) {
      Real value = fields(d, m, k, j, i);
      int bin;
      if (!output_diagnostics::IsFinite(value)) {
        Kokkos::atomic_exchange(&invalid(0), 1);
        bin = 0;
      } else if (value < d_min(d)) {
        bin = 0;
      } else if (!(value < d_max(d))) {
        bin = d_nbin(d) + 1;
      } else {
        Real transformed = PDFTransformValue(value, d_scale(d), d_linthresh(d));
        Real position = (transformed - d_transform_min(d))/d_step(d);
        if (!output_diagnostics::IsFinite(transformed) ||
            !output_diagnostics::IsFinite(position)) {
          Kokkos::atomic_exchange(&invalid(0), 1);
          bin = 0;
        } else {
          bin = static_cast<int>(position) + 1;
          if (bin < 1) bin = 1;
          if (bin > d_nbin(d)) bin = d_nbin(d);
        }
      }
      flat_index += bin*d_stride(d);
    }

    Real weight = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    if (weight_mode == 1) {
      weight *= density_data(m, IDN, k, j, i);
    } else if (weight_mode == 2) {
      weight *= fields(weight_index, m, k, j, i);
    }
    if (!output_diagnostics::IsFinite(weight)) {
      Kokkos::atomic_exchange(&invalid(0), 1);
    }
    Kokkos::atomic_add(&result(flat_index), weight);
  });

  Kokkos::fence();
  auto host_invalid = Kokkos::create_mirror_view(invalid);
  Kokkos::deep_copy(host_invalid, invalid);
  Kokkos::fence();
  if (host_invalid(0) != 0) {
    FatalPDFError("PDF output encountered a nonfinite axis, transform, or weight.");
  }
  ValidatePDFResult(result, pdf_data.total_bins, "Local PDF result");

#if MPI_PARALLEL_ENABLED
  if (out_params.shard_mode == FileShardMode::shared) {
    ReducePDFResultViaHost(result, pdf_data.total_bins, MPI_COMM_WORLD, 0,
                           "MPI_Reduce for shared PDF output");
  } else if (IsNodeSharded(out_params.shard_mode)) {
    ReducePDFResultViaHost(result, pdf_data.total_bins, global_variable::node_comm, 0,
                           "MPI_Reduce for node PDF output");
  }
#endif
  if (out_params.shard_mode == FileShardMode::shared) {
    if (global_variable::my_rank == 0) {
      ValidatePDFResult(result, pdf_data.total_bins, "Reduced shared PDF result");
    }
  } else if (IsNodeSharded(out_params.shard_mode)) {
    if (global_variable::node_rank == 0) {
      ValidatePDFResult(result, pdf_data.total_bins, "Reduced node PDF result");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Writes an explicitly versioned PDF header and binary payload.
//!
//! V2 binary payload starts with: 8-byte magic "AKPDFV2", uint32 version,
//! uint32 layout (0=dense, 1=sparse COO), uint32 ndim, uint32 writer rank,
//! uint64 record count, double time, int64 cycle. Dense records are float64 values in
//! flattened row-major order. Sparse records are repeated (uint64 flat index, float64).

void PDFOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (out_params.pdf_legacy_layout) {
    if (global_variable::my_rank == 0) {
      RequirePDFBudget(output_file_utils::CheckedSizeAdd(
          persistent_writer_allocation_bytes, PDFResultMirrorBytes(pdf_data.total_bins),
          "legacy PDF write allocation", FatalPDFError),
          max_writer_allocation_bytes, "legacy PDF write allocation");
      std::string path = PDFDirectory(out_params) + "/";
      if (!pdf_data.bins_written) {
        std::string header_name = path + out_params.file_basename + ".bins.pdf";
        std::FILE *header = std::fopen(header_name.c_str(), "a");
        if (header == nullptr) {
          FatalPDFError("Cannot open legacy PDF header '" + header_name + "'.");
        }
        CheckedLegacyPDFPrint(header, header_name, "# pdf bins \n");
        CheckedLegacyPDFPrint(header, header_name, "# [1]= %.20s \n",
                              outvars[0].label.c_str());
        if (pdf_data.ndim == 2) {
          CheckedLegacyPDFPrint(header, header_name, "# [2]= %.20s \n",
                                outvars[1].label.c_str());
        }
        for (int d = 0; d < pdf_data.ndim; ++d) {
          auto edges = Kokkos::create_mirror_view(pdf_data.bin_edges[d]);
          Kokkos::deep_copy(edges, pdf_data.bin_edges[d]);
          Kokkos::fence();
          for (int n = 0; n <= pdf_data.nbin[d]; ++n) {
            CheckedLegacyPDFPrint(header, header_name, out_params.data_format.c_str(),
                                  edges(n));
          }
          CheckedLegacyPDFPrint(header, header_name, "\n");
        }
        CheckedLegacyPDFClose(header, header_name);
        pdf_data.bins_written = true;
      }

      std::string sequence = output_file_utils::FormatSequence(
          out_params.file_number, "legacy PDF output", FatalPDFError);
      std::string data_name = path + out_params.file_basename + "." + sequence + ".pdf";
      std::FILE *output = std::fopen(data_name.c_str(), "a");
      if (output == nullptr) {
        FatalPDFError("Cannot open legacy PDF data file '" + data_name + "'.");
      }
      auto values = Kokkos::create_mirror_view(pdf_data.result_);
      Kokkos::deep_copy(values, pdf_data.result_);
      Kokkos::fence();
      CheckedLegacyPDFPrint(output, data_name, "# time= ");
      CheckedLegacyPDFPrint(output, data_name, out_params.data_format.c_str(), pm->time);
      CheckedLegacyPDFPrint(output, data_name, "\n");
      int rows = (pdf_data.ndim == 2) ? pdf_data.nbin_with_overflow[1] : 1;
      for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < pdf_data.nbin_with_overflow[0]; ++x) {
          int flat_index = x*pdf_data.stride[0] + y;
          CheckedLegacyPDFPrint(output, data_name, out_params.data_format.c_str(),
                                values(flat_index));
        }
        CheckedLegacyPDFPrint(output, data_name, "\n");
      }
      CheckedLegacyPDFPrint(output, data_name, "\n");
      CheckedLegacyPDFClose(output, data_name);
    }
    AdvanceOutputCounters(out_params, pm, pin);
    return;
  }

  bool sharded = IsSharded(out_params.shard_mode);
  bool i_write = IsRankSharded(out_params.shard_mode) ||
      (IsNodeSharded(out_params.shard_mode) && global_variable::node_rank == 0) ||
      (out_params.shard_mode == FileShardMode::shared && global_variable::my_rank == 0);
  if (i_write) {
    std::size_t payload_staging = output_file_utils::CheckedSizeProduct(
        PDFCountAsSize(pdf_data.total_bins, "PDF total bin count"),
        sizeof(Real) + sizeof(std::uint64_t) + sizeof(double),
        "PDF payload staging", FatalPDFError);
    RequirePDFBudget(output_file_utils::CheckedSizeAdd(
        persistent_writer_allocation_bytes, payload_staging, "PDF write allocation",
        FatalPDFError), max_writer_allocation_bytes, "PDF write allocation");
    std::string path = PDFDirectory(out_params) + "/";
    if (sharded) {
      path += ShardDirectoryName(out_params.shard_mode, global_variable::my_rank,
                                 global_variable::node_id) + "/";
    }

    if (!pdf_data.bins_written) {
      std::string header_name = path + out_params.file_basename + ".header.pdf";
      std::string temporary_header_name = header_name + ".tmp";
      std::FILE *header = std::fopen(temporary_header_name.c_str(), "w");
      if (header == nullptr) {
        FatalPDFError("Cannot open PDF header '" + temporary_header_name + "'.");
      }
      CheckedPDFPrint(header, temporary_header_name, "# AthenaK PDF format version=2\n");
      CheckedPDFPrint(header, temporary_header_name, "binary_magic = AKPDFV2\n");
      CheckedPDFPrint(header, temporary_header_name, "layout = %s\n",
                      sharded ? "sparse_coo" : "dense");
      CheckedPDFPrint(header, temporary_header_name, "distribution = %s\n",
                      ShardDistributionName(out_params.shard_mode));
      if (IsNodeSharded(out_params.shard_mode)) {
        CheckedPDFPrint(header, temporary_header_name, "node = %d\n",
                        global_variable::node_id);
        CheckedPDFPrint(header, temporary_header_name, "number_of_nodes = %d\n",
                        global_variable::nnodes);
      } else if (IsRankSharded(out_params.shard_mode)) {
        CheckedPDFPrint(header, temporary_header_name, "rank = %d\n",
                        global_variable::my_rank);
        CheckedPDFPrint(header, temporary_header_name, "number_of_ranks = %d\n",
                        global_variable::nranks);
      }
      CheckedPDFPrint(header, temporary_header_name, "ndim = %d\n", pdf_data.ndim);
      CheckedPDFPrint(header, temporary_header_name, "total_bins = %d\n",
                      pdf_data.total_bins);
      CheckedPDFPrint(header, temporary_header_name, "weight = %s\n",
                      out_params.pdf_weight.c_str());
      if (out_params.pdf_weight == "variable") {
        CheckedPDFPrint(header, temporary_header_name, "weight_variable = %s\n",
                        out_params.pdf_weight_variable.c_str());
      }
      CheckedPDFPrint(header, temporary_header_name,
                      "symlog_transform = sign(x)*(abs(x)/linthresh if "
                      "abs(x)<=linthresh else 1+log10(abs(x)/linthresh))\n");
      for (int d = 0; d < pdf_data.ndim; ++d) {
        CheckedPDFPrint(header, temporary_header_name, "variable_%d = %s\n", d + 1,
                        out_params.pdf_variables[d].c_str());
        CheckedPDFPrint(header, temporary_header_name, "nbin%d = %d\n", d + 1,
                        pdf_data.nbin[d]);
        CheckedPDFPrint(header, temporary_header_name, "bin%d_min = %.17e\n", d + 1,
                        pdf_data.bin_min[d]);
        CheckedPDFPrint(header, temporary_header_name, "bin%d_max = %.17e\n", d + 1,
                        pdf_data.bin_max[d]);
        CheckedPDFPrint(header, temporary_header_name, "scale%d = %s\n", d + 1,
                        PDFScaleName(pdf_data.scale[d]));
        if (pdf_data.scale[d] == PDF_SCALE_SYMLOG) {
          CheckedPDFPrint(header, temporary_header_name, "linthresh%d = %.17e\n",
                          d + 1, pdf_data.linthresh[d]);
        }
        CheckedPDFPrint(header, temporary_header_name, "stride%d = %d\n", d + 1,
                        pdf_data.stride[d]);
        auto edges = Kokkos::create_mirror_view(pdf_data.bin_edges[d]);
        Kokkos::deep_copy(edges, pdf_data.bin_edges[d]);
        Kokkos::fence();
        CheckedPDFPrint(header, temporary_header_name, "bin_edges_%d =", d + 1);
        for (int n = 0; n <= pdf_data.nbin[d]; ++n) {
          CheckedPDFPrint(header, temporary_header_name, " %.17e", edges(n));
        }
        CheckedPDFPrint(header, temporary_header_name, "\n");
      }
      PublishTemporaryPDFFile(header, temporary_header_name, header_name, "PDF header");
      pdf_data.bins_written = true;
    }

    std::string sequence = output_file_utils::FormatSequence(
        out_params.file_number, "PDF output", FatalPDFError);
    std::string data_name = path + out_params.file_basename + "." + sequence + ".pdf";
    std::string temporary_data_name = data_name + ".tmp";
    std::FILE *output = std::fopen(temporary_data_name.c_str(), "wb");
    if (output == nullptr) {
      FatalPDFError("Cannot open PDF data file '" + temporary_data_name + "'.");
    }
    auto values = Kokkos::create_mirror_view(pdf_data.result_);
    Kokkos::deep_copy(values, pdf_data.result_);
    Kokkos::fence();
    std::vector<std::uint64_t> sparse_indices;
    std::vector<double> sparse_values;
    if (sharded) {
      std::size_t nonzero_count = 0;
      for (int n = 0; n < pdf_data.total_bins; ++n) {
        if (values(n) != 0.0) {
          ++nonzero_count;
        }
      }
      sparse_indices.resize(nonzero_count);
      sparse_values.resize(nonzero_count);
      std::size_t sparse_index = 0;
      for (int n = 0; n < pdf_data.total_bins; ++n) {
        if (values(n) != 0.0) {
          sparse_indices[sparse_index] = static_cast<std::uint64_t>(n);
          sparse_values[sparse_index] = static_cast<double>(values(n));
          ++sparse_index;
        }
      }
    }
    const char magic[8] = {'A', 'K', 'P', 'D', 'F', 'V', '2', '\0'};
    std::uint32_t version = 2;
    std::uint32_t layout = sharded ? 1 : 0;
    std::uint32_t ndim = static_cast<std::uint32_t>(pdf_data.ndim);
    std::uint32_t rank = static_cast<std::uint32_t>(global_variable::my_rank);
    std::uint64_t count = sharded ? sparse_indices.size() :
        static_cast<std::uint64_t>(pdf_data.total_bins);
    double time = static_cast<double>(pm->time);
    std::int64_t cycle = static_cast<std::int64_t>(pm->ncycle);
    CheckedPDFWrite(output, magic, sizeof(char), 8, temporary_data_name,
                    "PDF binary magic");
    CheckedPDFWrite(output, &version, sizeof(version), 1, temporary_data_name,
                    "PDF format version");
    CheckedPDFWrite(output, &layout, sizeof(layout), 1, temporary_data_name,
                    "PDF layout");
    CheckedPDFWrite(output, &ndim, sizeof(ndim), 1, temporary_data_name,
                    "PDF dimensionality");
    CheckedPDFWrite(output, &rank, sizeof(rank), 1, temporary_data_name,
                    "PDF writer rank");
    CheckedPDFWrite(output, &count, sizeof(count), 1, temporary_data_name,
                    "PDF record count");
    CheckedPDFWrite(output, &time, sizeof(time), 1, temporary_data_name,
                    "PDF simulation time");
    CheckedPDFWrite(output, &cycle, sizeof(cycle), 1, temporary_data_name,
                    "PDF simulation cycle");
    if (sharded) {
      for (std::size_t n = 0; n < sparse_indices.size(); ++n) {
        CheckedPDFWrite(output, &(sparse_indices[n]), sizeof(std::uint64_t), 1,
                        temporary_data_name, "PDF sparse index");
        CheckedPDFWrite(output, &(sparse_values[n]), sizeof(double), 1,
                        temporary_data_name, "PDF sparse value");
      }
    } else {
      for (int n = 0; n < pdf_data.total_bins; ++n) {
        double value = static_cast<double>(values(n));
        CheckedPDFWrite(output, &value, sizeof(value), 1, temporary_data_name,
                        "PDF dense value");
      }
    }
    PublishTemporaryPDFFile(output, temporary_data_name, data_name, "PDF data file");
  }

  AdvanceOutputCounters(out_params, pm, pin);
}
