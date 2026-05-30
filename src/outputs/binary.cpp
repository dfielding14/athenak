//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file binary.cpp
//! \brief writes output data in binary format, which simply consists of each MeshBlock
//! written contiguously in order of "gid" in binary format.

#include <algorithm>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "outputs.hpp"

namespace {

std::string &ActiveBinaryTemporary() {
  static std::string path;
  return path;
}

void CleanupActiveBinaryTemporary() {
  output_file_utils::DiscardOwnedPath(ActiveBinaryTemporary());
}

[[noreturn]] void FatalBinaryError(const std::string &message) {
  CleanupActiveBinaryTemporary();
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

std::size_t CheckedAdd(std::size_t left, std::size_t right, const char *context) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    FatalBinaryError(std::string(context) + " size overflow.");
  }
  return left + right;
}

std::size_t CheckedProduct(std::size_t left, std::size_t right, const char *context) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max()/left) {
    FatalBinaryError(std::string(context) + " size overflow.");
  }
  return left*right;
}

std::size_t CountAsSize(int count, const char *context) {
  if (count < 0) {
    FatalBinaryError(std::string(context) + " is negative.");
  }
  return static_cast<std::size_t>(count);
}

std::size_t Count64AsSize(std::uint64_t count, const char *context) {
  if (count > std::numeric_limits<std::size_t>::max()) {
    FatalBinaryError(std::string(context) + " exceeds size_t range.");
  }
  return static_cast<std::size_t>(count);
}

int CountAsInt(std::size_t count, const char *context) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    FatalBinaryError(std::string(context) + " exceeds int range.");
  }
  return static_cast<int>(count);
}

std::size_t RankPrefixSum(const std::vector<int> &counts, int rank) {
  std::size_t prefix = 0;
  for (int r = 0; r < rank; ++r) {
    prefix = CheckedAdd(prefix, CountAsSize(counts[r], "binary rank MeshBlock count"),
                        "binary rank MeshBlock prefix");
  }
  return prefix;
}

void CheckedWrite(IOWrapper &file, const void *data, std::size_t count,
                  const char *context, bool independent_file) {
  if (file.Write_any_type(data, count, "byte", independent_file) != count) {
    FatalBinaryError(std::string(context) + " was not written completely.");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// Constructor: also calls BaseTypeOutput base class constructor

MeshBinaryOutput::MeshBinaryOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs
  // useful for mpiio-based outputs because on some supercomputers you may need to
  // set different stripe counts depending on whether mpiio is used in order to
  // achieve the best performance and not to crash the filesystem
  output_file_utils::EnsureDirectory("bin", 0775, "binary output", FatalBinaryError);
  if (IsSharded(op.shard_mode)) {
    std::string shard_dir = "bin/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    output_file_utils::EnsureDirectory(shard_dir, 0775, "binary output",
                                       FatalBinaryError);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBinaryOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in binary format
//   All MeshBlocks are written to the same file.

void MeshBinaryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  FileShardMode shard_mode = out_params.shard_mode;
  bool independent_file = UsesIndependentFileIO(shard_mode);
  bool shard_writer = IsRankSharded(shard_mode) ||
      (IsNodeSharded(shard_mode) && global_variable::node_rank == 0) ||
      (shard_mode == FileShardMode::shared && global_variable::my_rank == 0);

  // create filename: "bin/file_basename" + "." + "file_id" + "." + XXXXX + ".bin"
  // where XXXXX = file_number with a minimum width of 5 digits

  std::string sequence = output_file_utils::FormatSequence(
      out_params.file_number, "binary output", FatalBinaryError);
  std::string published_fname;
  if (IsSharded(shard_mode)) {
    published_fname = std::string("bin/") + ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id) + "/"
          + out_params.file_basename
          + "." + out_params.file_id + "." + sequence + ".bin";
  } else {
    published_fname = std::string("bin/") + out_params.file_basename
          + "." + out_params.file_id + "." + sequence + ".bin";
  }
  std::string fname = output_file_utils::TemporaryPath(published_fname);
  ActiveBinaryTemporary() = shard_writer ? fname : "";
  mpi_utils::SetFatalCleanupHook(CleanupActiveBinaryTemporary);

  IOWrapper binfile;
  std::size_t header_offset=0;
#if MPI_PARALLEL_ENABLED
  if (IsNodeSharded(shard_mode)) {
    binfile.SetCommunicator(global_variable::node_comm);
  }
#endif
  binfile.Open(fname.c_str(), IOWrapper::FileMode::write, independent_file);

  int nout_mbs = CountAsInt(outmbs.size(), "binary MeshBlock count");
  std::uint64_t shard_nout_mbs = nout_mbs;
  if (IsNodeSharded(shard_mode)) {
    shard_nout_mbs = global_variable::NodeSum64(nout_mbs);
  }

  // Basic parts of the format:
  // 1. Size of the header
  // 2. Current time
  // 3. List of variables in the file
  // 4. Header (input file information)
  {
    std::stringstream msg;
    const int time_precision = std::numeric_limits<Real>::max_digits10 - 1;
    msg << "Athena binary output version=1.1" << std::endl
        // preheader size includes "size of preheader" line up to "number of variables"
        << "  size of preheader=" << (IsNodeSharded(shard_mode) ? 9 : 5) << std::endl
        << std::scientific << std::setprecision(time_precision)
        << "  time=" << pm->time << std::endl
        << "  cycle=" << pm->ncycle << std::endl
        << "  size of location=" << sizeof(Real) << std::endl
        << "  size of variable=" << sizeof(float) << std::endl;
    if (IsNodeSharded(shard_mode)) {
      msg << "  distribution=node" << std::endl
          << "  node=" << global_variable::node_id << std::endl
          << "  number of nodes=" << global_variable::nnodes << std::endl
          << "  number of meshblocks=" << shard_nout_mbs << std::endl;
    }
    msg << "  number of variables=" << outvars.size() << std::endl
        << "  variables:  ";
    for (int n=0; n<outvars.size(); n++) {
      msg << outvars[n].label.c_str() << "  ";
    }
    msg << std::endl;
    std::string metadata = msg.str();
    if (shard_writer) {
      CheckedWrite(binfile, metadata.data(), metadata.size(), "binary metadata",
                   independent_file);
    }
    header_offset = CheckedAdd(header_offset, metadata.size(), "binary header");
  }
  {
    std::stringstream msg;
    // prepare the input parameters
    std::stringstream ost;
    pin->ParameterDump(ost);
    std::string sbuf=ost.str();
    msg << "  header offset=" << sbuf.size()*sizeof(char)  << std::endl;
    std::string offset_metadata = msg.str();
    if (shard_writer) {
      CheckedWrite(binfile, offset_metadata.data(), offset_metadata.size(),
                   "binary header-offset metadata", independent_file);
      CheckedWrite(binfile, sbuf.data(), sbuf.size(), "binary input header",
                   independent_file);
    }
    header_offset = CheckedAdd(header_offset, sbuf.size(), "binary header");
    header_offset = CheckedAdd(header_offset, offset_metadata.size(), "binary header");
  }

  //  5. Data.  An arbitrary number of scalars and vectors can be written (every element
  //  of the outvars vector), all in binary floats format

  int nout_vars = outvars.size();
  std::size_t cells = 0;
  if (nout_mbs > 0) {
    std::size_t nout1 = CountAsSize(outmbs[0].oie - outmbs[0].ois + 1,
                                    "binary x1 extent");
    std::size_t nout2 = CountAsSize(outmbs[0].oje - outmbs[0].ojs + 1,
                                    "binary x2 extent");
    std::size_t nout3 = CountAsSize(outmbs[0].oke - outmbs[0].oks + 1,
                                    "binary x3 extent");
    cells = CheckedProduct(CheckedProduct(nout1, nout2, "binary cell count"),
                           nout3, "binary cell count");
  }
#if MPI_PARALLEL_ENABLED
  if (!IsRankSharded(shard_mode)) {
    std::uint64_t local_cells = cells;
    std::uint64_t shard_cells = 0;
    MPI_Comm comm =
        IsNodeSharded(shard_mode) ? global_variable::node_comm : MPI_COMM_WORLD;
    mpi_utils::CheckMpi(
        MPI_Allreduce(&local_cells, &shard_cells, 1, MPI_UINT64_T, MPI_MAX, comm),
        "MPI_Allreduce for binary MeshBlock cell counts");
    if (shard_cells > std::numeric_limits<std::size_t>::max()) {
      FatalBinaryError("binary MeshBlock cell count exceeds size_t range.");
    }
    cells = static_cast<std::size_t>(shard_cells);
  }
#endif

  // ois, oie, ojs, oje, oks, oke + il1, il2, il3, level +
  // x1min, x1max, x2min, x2max, x3min, x3max + data
  std::size_t value_bytes = CheckedProduct(
      CheckedProduct(cells, CountAsSize(nout_vars, "binary variable count"),
                     "binary MeshBlock values"),
      sizeof(float), "binary MeshBlock values");
  std::size_t data_size = CheckedAdd(10*sizeof(int32_t) + 6*sizeof(Real),
                                     value_bytes, "binary MeshBlock record");

  std::size_t node_offset = IsNodeSharded(shard_mode) ? Count64AsSize(
      global_variable::NodePrefixSum64(nout_mbs), "binary node MeshBlock prefix") : 0;
  std::size_t payload_bytes = CheckedProduct(
      CountAsSize(nout_mbs, "binary MeshBlock count"), data_size, "binary payload");

  // allocate 1D vector of floats used to convert and output data
  std::vector<char> data(payload_bytes);
  std::vector<float> single_data(cells);

  // Loop over MeshBlocks
  for (int m=0; m<nout_mbs; ++m) {
    char *pdata = data.data() + CheckedProduct(
        CountAsSize(m, "binary local MeshBlock index"), data_size,
        "binary local MeshBlock offset");
    LogicalLocation loc = pm->lloc_eachmb[outmbs[m].mb_gid];
    int &ois = outmbs[m].ois;
    int &oie = outmbs[m].oie;
    int &ojs = outmbs[m].ojs;
    int &oje = outmbs[m].oje;
    int &oks = outmbs[m].oks;
    int &oke = outmbs[m].oke;

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
      std::size_t variable_bytes = CheckedProduct(cells, sizeof(float),
                                                  "binary variable payload");
      memcpy(pdata, single_data.data(), variable_bytes);
      pdata += variable_bytes;
    }
  }

  // now write binary data
  std::size_t block_offset = 0;
  if (shard_mode == FileShardMode::shared) {
    block_offset = RankPrefixSum(noutmbs, global_variable::my_rank);
  } else if (IsNodeSharded(shard_mode)) {
    block_offset = node_offset;
  }
  std::size_t myoffset = CheckedAdd(
      header_offset, CheckedProduct(data_size, block_offset, "binary payload offset"),
      "binary payload offset");
  char dummy = '\0';
  const char *payload = data.empty() ? &dummy : data.data();
  if (binfile.Write_any_type_at_all(payload, payload_bytes, myoffset, "byte",
                                    independent_file) != payload_bytes) {
    if (shard_writer) output_file_utils::DiscardOwnedPath(fname);
    FatalBinaryError("binary payload was not written completely.");
  }

  // close the output file and clean up ptrs to data
  if (binfile.Close(independent_file) != 0) {
    if (shard_writer) output_file_utils::DiscardOwnedPath(fname);
    FatalBinaryError("Could not close binary output file '" + fname + "'.");
  }
  if (shard_writer) {
    output_file_utils::PublishTemporaryFile(fname, published_fname, "binary output",
                                            FatalBinaryError);
  }
  ActiveBinaryTemporary().clear();
  mpi_utils::SetFatalCleanupHook(nullptr);
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                      "MPI_Barrier after binary output publication");
#endif

  // increment counters
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "binary output", FatalBinaryError);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}
