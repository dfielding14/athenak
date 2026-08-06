//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coarsened_binary.cpp
//! \brief writes output data in binary format, which simply consists of each MeshBlock
//! written contiguously in order of "gid" in binary format.

#include <algorithm>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coarsened_binary_layout.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"

namespace {

std::string &ActiveCoarsenedBinaryTemporary() {
  static std::string path;
  return path;
}

void CleanupActiveCoarsenedBinaryTemporary() {
  output_file_utils::DiscardOwnedPath(ActiveCoarsenedBinaryTemporary());
}

[[noreturn]] void FatalCoarsenedBinaryError(const std::string &message) {
  CleanupActiveCoarsenedBinaryTemporary();
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

std::size_t CheckedAdd(std::size_t left, std::size_t right, const char *context) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    FatalCoarsenedBinaryError(std::string(context) + " size overflow.");
  }
  return left + right;
}

std::size_t CheckedProduct(std::size_t left, std::size_t right, const char *context) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max()/left) {
    FatalCoarsenedBinaryError(std::string(context) + " size overflow.");
  }
  return left*right;
}

std::size_t CountAsSize(int count, const char *context) {
  if (count < 0) {
    FatalCoarsenedBinaryError(std::string(context) + " is negative.");
  }
  return static_cast<std::size_t>(count);
}

std::size_t Count64AsSize(std::uint64_t count, const char *context) {
  if (count > std::numeric_limits<std::size_t>::max()) {
    FatalCoarsenedBinaryError(std::string(context) + " exceeds size_t range.");
  }
  return static_cast<std::size_t>(count);
}

int CountAsInt(std::size_t count, const char *context) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    FatalCoarsenedBinaryError(std::string(context) + " exceeds int range.");
  }
  return static_cast<int>(count);
}

std::size_t RankPrefixSum(const std::vector<int> &counts, int rank) {
  std::size_t prefix = 0;
  for (int r = 0; r < rank; ++r) {
    prefix = CheckedAdd(
        prefix, CountAsSize(counts[r], "coarsened-binary rank MeshBlock count"),
        "coarsened-binary rank MeshBlock prefix");
  }
  return prefix;
}

void CheckedWrite(IOWrapper &file, const void *data, std::size_t count,
                  const char *context, bool independent_file) {
  if (file.Write_any_type(data, count, "byte", independent_file) != count) {
    FatalCoarsenedBinaryError(std::string(context) + " was not written completely.");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// Constructor: also calls BaseTypeOutput base class constructor

CoarsenedBinaryOutput::CoarsenedBinaryOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  if (out_params.slice1 || out_params.slice2 || out_params.slice3) {
    FatalCoarsenedBinaryError(
        "Sliced coarsened-binary output is not supported.");
  }
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  if (pm->multilevel) {
    FatalCoarsenedBinaryError(
        "Coarsened-binary output supports uniform meshes only; static refinement "
        "and AMR are not supported.");
  }
  if (indcs.nx2 <= 1 || indcs.nx3 <= 1) {
    FatalCoarsenedBinaryError(
        "Coarsened-binary output supports three-dimensional meshes only.");
  }
  if (out_params.include_gzs) {
    FatalCoarsenedBinaryError(
        "Coarsened-binary output does not support ghost_zones=true.");
  }
  coarsened_binary_layout::KernelLayout::Build(
      indcs.nx1, indcs.nx2, indcs.nx3, out_params.coarsen_factor,
      out_params.compute_moments ? 4 : 1, FatalCoarsenedBinaryError);
  // create directories for outputs
  // useful for mpiio-based outputs because on some supercomputers you may need to
  // set different stripe counts depending on whether mpiio is used in order to
  // achieve the best performance and not to crash the filesystem
  std::string dir_name;
  dir_name.assign("cbin_");
  dir_name.append(out_params.file_id);
  dir_name.append("_");
  dir_name.append(std::to_string(out_params.coarsen_factor));
  output_file_utils::EnsureDirectory(dir_name, 0775, "coarsened-binary output",
                                     FatalCoarsenedBinaryError);
  if (IsSharded(op.shard_mode)) {
    dir_name.append("/");
    dir_name.append(ShardDirectoryName(op.shard_mode, global_variable::my_rank,
                                       global_variable::node_id));
    output_file_utils::EnsureDirectory(dir_name, 0775, "coarsened-binary output",
                                       FatalCoarsenedBinaryError);
  }
}

//----------------------------------------------------------------------------------------
// BaseTypeOutput::LoadOutputData()
// create std::vector of HostArray3Ds containing data specified in <output> block for
// this output type

void CoarsenedBinaryOutput::LoadOutputData(Mesh *pm) {
  // out_data_ vector (indexed over # of output MBs) stores 4D array of variables
  // so start iteration over number of MeshBlocks
  // Recompute the emitted MeshBlock inventory at every output time.
  outmbs.clear();

  // loop over all MeshBlocks
  // set size & starting indices of output arrays, adjusted accordingly if gz included
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  auto &size  = pm->pmb_pack->pmb->mb_size;
  for (int m=0; m<(pm->pmb_pack->nmb_thispack); ++m) {
    int id = pm->pmb_pack->pmb->mb_gid.h_view(m);
    if (out_params.gid >= 0 && id != out_params.gid) { continue; }

    int ois = indcs.is;
    int oie = indcs.ie;
    int ojs = indcs.js;
    int oje = indcs.je;
    int oks = indcs.ks;
    int oke = indcs.ke;

    // set coordinate geometry information for MB
    Real x1min = size.h_view(m).x1min;
    Real x1max = size.h_view(m).x1max;
    Real x2min = size.h_view(m).x2min;
    Real x2max = size.h_view(m).x2max;
    Real x3min = size.h_view(m).x3min;
    Real x3max = size.h_view(m).x3max;

    outmbs.emplace_back(id,ois,oie,ojs,oje,oks,oke,x1min,x1max,x2min,x2max,x3min,x3max);
  }

  std::fill(noutmbs.begin(), noutmbs.end(), 0);
  noutmbs[global_variable::my_rank] =
      CountAsInt(outmbs.size(), "coarsened-binary MeshBlock count");
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(
      MPI_Allreduce(MPI_IN_PLACE, noutmbs.data(), global_variable::nranks,
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD),
      "MPI_Allreduce for coarsened-binary MeshBlock counts");
#endif
  noutmbs_min = *std::min_element(noutmbs.begin(), noutmbs.end());
  noutmbs_max = *std::max_element(noutmbs.begin(), noutmbs.end());

  // get number of output vars and MBs, then realloc outarray (HostArray)
  int nout_vars_with_moments;
  if (out_params.compute_moments) {
    nout_vars_with_moments = CountAsInt(
        CheckedProduct(outvars.size(), 4, "coarsened-binary variable count"),
        "coarsened-binary variable count");
  } else {
    nout_vars_with_moments =
        CountAsInt(outvars.size(), "coarsened-binary variable count");
  }
  int nout_vars = CountAsInt(outvars.size(), "coarsened-binary variable count");
  int nout_mbs = CountAsInt(outmbs.size(), "coarsened-binary MeshBlock count");
  // note that while ois,oie,etc. can be different on each MB, the number of cells output
  // on each MeshBlock, i.e. (ois-ois+1), etc. is the same.
  if (nout_mbs > 0) {
    auto layout = coarsened_binary_layout::KernelLayout::Build(
        outmbs[0].oie - outmbs[0].ois + 1, outmbs[0].oje - outmbs[0].ojs + 1,
        outmbs[0].oke - outmbs[0].oks + 1, out_params.coarsen_factor,
        out_params.compute_moments ? 4 : 1, FatalCoarsenedBinaryError);
    std::size_t allocation_elements =
        coarsened_binary_layout::CheckedAllocationElements(
            nout_vars_with_moments, nout_mbs, layout, FatalCoarsenedBinaryError);
    coarsened_binary_layout::CheckedAllocationBytes(
        allocation_elements, sizeof(Real), FatalCoarsenedBinaryError);
    Kokkos::realloc(outarray, nout_vars_with_moments, nout_mbs,
                    layout.coarsened_nout3, layout.coarsened_nout2,
                    layout.coarsened_nout1);
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
      std::size_t input_cells = CheckedProduct(
          CheckedProduct(CountAsSize(nout1, "coarsened-binary x1 extent"),
                         CountAsSize(nout2, "coarsened-binary x2 extent"),
                         "coarsened-binary input allocation"),
          CountAsSize(nout3, "coarsened-binary x3 extent"),
          "coarsened-binary input allocation");
      CheckedProduct(input_cells, sizeof(Real), "coarsened-binary input allocation");
      // copy output variable to new device View
      DvceArray3D<Real> d_output_var("d_out_var",nout3,nout2,nout1);
      auto d_slice = Kokkos::subview(*(outvars[n].data_ptr), mbi, outvars[n].data_index,
                                     krange,jrange,irange);
      Kokkos::deep_copy(d_output_var,d_slice);


      int number_of_moments = 1;
      if (out_params.compute_moments) {
        number_of_moments = 4;
      }
      auto layout = coarsened_binary_layout::KernelLayout::Build(
          nout1, nout2, nout3, out_params.coarsen_factor, number_of_moments,
          FatalCoarsenedBinaryError);
      std::size_t coarsened_elements = CheckedProduct(
          CountAsSize(number_of_moments, "coarsened-binary moment count"),
          layout.coarsened_cells, "coarsened-binary temporary allocation");
      coarsened_binary_layout::CheckedAllocationBytes(
          coarsened_elements, sizeof(Real), FatalCoarsenedBinaryError);
      DvceArray4D<Real> d_output_var_coarsened("d_output_var_coarsened",
        number_of_moments, layout.coarsened_nout3, layout.coarsened_nout2,
        layout.coarsened_nout1);

      // Coarsen the d_slice and store the result in d_output_var
      // CoarsenVariable(d_output_var, d_output_var_coarsened, out_params.coarsen_factor);
      bool compute_moments = out_params.compute_moments;
      Kokkos::parallel_for("coarsen_variable",
       Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<std::int64_t>>(
           0, layout.coarsen_iterations),
      KOKKOS_LAMBDA(const std::int64_t idx) {
        // Calculate the 3D indices for the coarsened data
        std::int64_t total_coarsened_elements =
            static_cast<std::int64_t>(layout.coarsened_cells);
        std::int64_t k_c = (idx / layout.coarsened_plane)
            % layout.coarsened_nout3;
        std::int64_t j_c = (idx / layout.coarsened_nout1) % layout.coarsened_nout2;
        std::int64_t i_c = idx % layout.coarsened_nout1;

        // Calculate the offset within the coarsen_factor_cubed cube
        std::int64_t offset = idx / total_coarsened_elements;
        std::int64_t kk = offset / layout.coarsen_factor_squared;
        std::int64_t jj =
            (offset / layout.coarsen_factor) % layout.coarsen_factor;
        std::int64_t ii = offset % layout.coarsen_factor;

        // Calculate the corresponding indices in the full data
        std::int64_t k = k_c * layout.coarsen_factor + kk;
        std::int64_t j = j_c * layout.coarsen_factor + jj;
        std::int64_t i = i_c * layout.coarsen_factor + ii;

        // Perform the coarsening operation
        if(k < nout3 && j < nout2 && i < nout1) {
          Kokkos::atomic_add(&d_output_var_coarsened(0, k_c, j_c, i_c),
            d_output_var(k, j, i));
          if (compute_moments) {
            Kokkos::atomic_add(&d_output_var_coarsened(1, k_c, j_c, i_c),
              d_output_var(k, j, i)*d_output_var(k, j, i));
            Kokkos::atomic_add(&d_output_var_coarsened(2, k_c, j_c, i_c),
              d_output_var(k, j, i)*d_output_var(k, j, i)*d_output_var(k, j, i));
            Kokkos::atomic_add(&d_output_var_coarsened(3, k_c, j_c, i_c),
               d_output_var(k, j, i)*d_output_var(k, j, i)
              *d_output_var(k, j, i)*d_output_var(k, j, i));
          }
        }
      });
      // Normalize the coarsened data
      Kokkos::parallel_for("normalize_coarsened_variable",
        Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<std::int64_t>>(
            0, layout.normalize_iterations),
      KOKKOS_LAMBDA(const std::int64_t idx) {
        std::int64_t total_coarsened_elements =
            static_cast<std::int64_t>(layout.coarsened_cells);
        std::int64_t moment_idx = idx / total_coarsened_elements;
        std::int64_t k = (idx / layout.coarsened_plane)
            % layout.coarsened_nout3;
        std::int64_t j = (idx / layout.coarsened_nout1) % layout.coarsened_nout2;
        std::int64_t i = idx % layout.coarsened_nout1;

        d_output_var_coarsened(moment_idx, k, j, i) /=
            layout.coarsen_factor_cubed_kernel;
      });


      // Now, create a host mirror for the coarsened data.
      DvceArray4D<Real>::HostMirror h_output_var = Kokkos::create_mirror(
        d_output_var_coarsened
      );

      // Copy the coarsened data to the host mirror.
      Kokkos::deep_copy(h_output_var, d_output_var_coarsened);

      // copy host mirror to 5D host View containing all output variables
      // if (out_params.compute_moments) {
      auto h_slice = Kokkos::subview(outarray,
        moment_range,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL
      );
      Kokkos::deep_copy(h_slice,h_output_var);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CoarsenedBinaryOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in Coarsenedbinary format
//   All MeshBlocks are written to the same file.

void CoarsenedBinaryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // create filename: "cbin_"+"file_id"+"_"+"coarsening_factor"+"/file_basename"
  // + "." + "file_id" + "." + XXXXX + ".cbin"
  // where XXXXX = file_number with a minimum width of 5 digits
  FileShardMode shard_mode = out_params.shard_mode;
  bool independent_file = UsesIndependentFileIO(shard_mode);
  bool shard_writer = IsRankSharded(shard_mode) ||
      (IsNodeSharded(shard_mode) && global_variable::node_rank == 0) ||
      (shard_mode == FileShardMode::shared && global_variable::my_rank == 0);
  std::string published_fname;
  std::string number = output_file_utils::FormatSequence(
      out_params.file_number, "coarsened-binary output", FatalCoarsenedBinaryError);

  published_fname.assign("cbin_");
  published_fname.append(out_params.file_id);
  published_fname.append("_");
  published_fname.append(std::to_string(out_params.coarsen_factor));
  published_fname.append("/");
  if (IsSharded(shard_mode)) {
    published_fname.append(ShardDirectoryName(shard_mode, global_variable::my_rank,
                                              global_variable::node_id));
    published_fname.append("/");
  }
  published_fname.append(out_params.file_basename);
  published_fname.append(".");
  published_fname.append(out_params.file_id);
  published_fname.append(".");
  published_fname.append(number);
  published_fname.append(".cbin");
  std::string fname = output_file_utils::TemporaryPath(published_fname);
  ActiveCoarsenedBinaryTemporary() = shard_writer ? fname : "";
  mpi_utils::SetFatalCleanupHook(CleanupActiveCoarsenedBinaryTemporary);

  IOWrapper cbinfile;
  std::size_t header_offset=0;
#if MPI_PARALLEL_ENABLED
  if (IsNodeSharded(shard_mode)) {
    cbinfile.SetCommunicator(global_variable::node_comm);
  }
#endif
  cbinfile.Open(fname.c_str(), IOWrapper::FileMode::write, independent_file);

  int number_of_moments = 1;
  if (out_params.compute_moments) {
    number_of_moments = 4;
  }
  int nout_vars = CountAsInt(
      CheckedProduct(outvars.size(), CountAsSize(number_of_moments,
                                                 "coarsened-binary moment count"),
                     "coarsened-binary variable count"),
      "coarsened-binary variable count");

  // Basic parts of the format:
  // 1. Size of the header
  // 2. Current time
  // 3. List of variables in the file
  // 4. Header (input file information)
  int nout_mbs = CountAsInt(outmbs.size(), "coarsened-binary MeshBlock count");
  std::uint64_t shard_nout_mbs = IsNodeSharded(shard_mode) ?
      global_variable::NodeSum64(nout_mbs) : nout_mbs;
  {std::stringstream msg;
  msg << "Athena binary output version=1.1" << std::endl
      // preheader size includes "size of preheader" line up to "number of variables"
      << "  size of preheader=" << (IsNodeSharded(shard_mode) ? 11 : 7) << std::endl
      << "  time=" << pm->time << std::endl
      << "  cycle=" << pm->ncycle << std::endl
      << "  number of moments=" << number_of_moments << std::endl
      << "  coarsening factor=" << out_params.coarsen_factor << std::endl
      << "  size of location=" << sizeof(Real) << std::endl
      << "  size of variable=" << sizeof(float) << std::endl;
  if (IsNodeSharded(shard_mode)) {
    msg << "  distribution=node" << std::endl
        << "  node=" << global_variable::node_id << std::endl
        << "  number of nodes=" << global_variable::nnodes << std::endl
        << "  number of meshblocks=" << shard_nout_mbs << std::endl;
  }
  msg << "  number of variables=" << nout_vars << std::endl
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
  std::string metadata = msg.str();
  if (shard_writer) {
    CheckedWrite(cbinfile, metadata.data(), metadata.size(),
                 "coarsened-binary metadata", independent_file);
  }
  header_offset = CheckedAdd(header_offset, metadata.size(), "coarsened-binary header");}
  {std::stringstream msg;
  // prepare the input parameters
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf=ost.str();
  msg << "  header offset=" << sbuf.size()*sizeof(char)  << std::endl;
  std::string offset_metadata = msg.str();
  if (shard_writer) {
    CheckedWrite(cbinfile, offset_metadata.data(), offset_metadata.size(),
                 "coarsened-binary header-offset metadata", independent_file);
    CheckedWrite(cbinfile, sbuf.data(), sbuf.size(), "coarsened-binary input header",
                 independent_file);
  }
  header_offset = CheckedAdd(header_offset, sbuf.size(), "coarsened-binary header");
  header_offset = CheckedAdd(header_offset, offset_metadata.size(),
                             "coarsened-binary header");}

  //  5. Data.  An arbitrary number of scalars and vectors can be written (every node
  //  in the OutputData doubly linked lists), all in binary floats format

  int nout1 = 0;
  int nout2 = 0;
  int nout3 = 0;
  if (nout_mbs > 0) {
    nout1 = ((outmbs[0].oie - outmbs[0].ois + 1)/out_params.coarsen_factor);
    nout2 = ((outmbs[0].oje - outmbs[0].ojs + 1)/out_params.coarsen_factor);
    nout3 = ((outmbs[0].oke - outmbs[0].oks + 1)/out_params.coarsen_factor);
  }
  std::size_t cells = CheckedProduct(
      CheckedProduct(CountAsSize(nout1, "coarsened-binary x1 extent"),
                     CountAsSize(nout2, "coarsened-binary x2 extent"),
                     "coarsened-binary cell count"),
      CountAsSize(nout3, "coarsened-binary x3 extent"), "coarsened-binary cell count");
#if MPI_PARALLEL_ENABLED
  if (!IsRankSharded(shard_mode)) {
    std::uint64_t local_cells = cells;
    std::uint64_t shard_cells = 0;
    MPI_Comm comm =
        IsNodeSharded(shard_mode) ? global_variable::node_comm : MPI_COMM_WORLD;
    mpi_utils::CheckMpi(
        MPI_Allreduce(&local_cells, &shard_cells, 1, MPI_UINT64_T, MPI_MAX, comm),
        "MPI_Allreduce for coarsened-binary MeshBlock cell counts");
    if (shard_cells > std::numeric_limits<std::size_t>::max()) {
      FatalCoarsenedBinaryError(
          "coarsened-binary MeshBlock cell count exceeds size_t range.");
    }
    cells = static_cast<std::size_t>(shard_cells);
  }
#endif


  // ois, oie, ojs, oje, oks, oke + il1, il2, il3, level +
  // x1min, x1max, x2min, x2max, x3min, x3max + data
  std::size_t value_bytes = CheckedProduct(
      CheckedProduct(cells, CountAsSize(nout_vars, "coarsened-binary variable count"),
                     "coarsened-binary MeshBlock values"),
      sizeof(float), "coarsened-binary MeshBlock values");
  std::size_t data_size = CheckedAdd(10*sizeof(int32_t) + 6*sizeof(Real),
                                     value_bytes, "coarsened-binary MeshBlock record");

  std::size_t node_offset = IsNodeSharded(shard_mode) ? Count64AsSize(
      global_variable::NodePrefixSum64(nout_mbs),
      "coarsened-binary node MeshBlock prefix") : 0;
  std::size_t payload_bytes = CheckedProduct(
      CountAsSize(nout_mbs, "coarsened-binary MeshBlock count"), data_size,
      "coarsened-binary payload");

  // allocate 1D vector of floats used to convert and output data
  std::vector<char> data(payload_bytes);
  std::vector<float> single_data(cells);

  // Loop over MeshBlocks
  for (int m=0; m<nout_mbs; ++m) {
    char *pdata = data.data() + CheckedProduct(
        CountAsSize(m, "coarsened-binary local MeshBlock index"), data_size,
        "coarsened-binary local MeshBlock offset");
    LogicalLocation loc = pm->lloc_eachmb[outmbs[m].mb_gid];
    // The constructor restricts cbin to active-zone, full-volume, uniform-grid
    // output, so the reduced-grid record starts at the active-zone indices.
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


    // Preserve the uniform-grid logical location. AMR is rejected during
    // construction because reduced-grid refinement semantics are not defined.
    nx = (int32_t)(loc.lx1);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(loc.lx2);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);
    nx = (int32_t)(loc.lx3);
    memcpy(pdata,&(nx),sizeof(nx));
    pdata+=sizeof(nx);

    // Uniform-grid output retains the existing relative-level field, which is
    // always zero under the supported contract.
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
      std::size_t cnt=0;
      for (int k=oks; k<=oke; k++) {
        for (int j=ojs; j<=oje; j++) {
          for (int i=ois; i<=oie; i++) {
            tmp_data = static_cast<float>(outarray(n,m,k-oks,j-ojs,i-ois));
            single_data[cnt] = tmp_data;
            cnt++;
          }
        }
      }
      std::size_t variable_bytes = CheckedProduct(
          cells, sizeof(float), "coarsened-binary variable payload");
      memcpy(pdata, single_data.data(), variable_bytes);
      pdata += variable_bytes;
    }
  }

  // now write Coarsenedbinary data
  std::size_t block_offset = 0;
  if (shard_mode == FileShardMode::shared) {
    block_offset = RankPrefixSum(noutmbs, global_variable::my_rank);
  } else if (IsNodeSharded(shard_mode)) {
    block_offset = node_offset;
  }
  std::size_t myoffset = CheckedAdd(
      header_offset,
      CheckedProduct(data_size, block_offset, "coarsened-binary payload offset"),
      "coarsened-binary payload offset");
  char dummy = '\0';
  const char *payload = data.empty() ? &dummy : data.data();
  if (cbinfile.Write_any_type_at_all(payload, payload_bytes, myoffset, "byte",
                                     independent_file) != payload_bytes) {
    if (shard_writer) output_file_utils::DiscardOwnedPath(fname);
    FatalCoarsenedBinaryError("coarsened-binary payload was not written completely.");
  }

  // close the output file and clean up ptrs to data
  if (cbinfile.Close(independent_file) != 0) {
    if (shard_writer) output_file_utils::DiscardOwnedPath(fname);
    FatalCoarsenedBinaryError(
        "Could not close coarsened-binary output file '" + fname + "'.");
  }
  if (shard_writer) {
    output_file_utils::PublishTemporaryFile(
        fname, published_fname, "coarsened-binary output", FatalCoarsenedBinaryError);
  }
  ActiveCoarsenedBinaryTemporary().clear();
  mpi_utils::SetFatalCleanupHook(nullptr);
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                      "MPI_Barrier after coarsened-binary output publication");
#endif

  // increment counters
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "coarsened-binary output", FatalCoarsenedBinaryError);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}
