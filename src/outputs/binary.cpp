//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file binary.cpp
//! \brief writes output data in binary format, which simply consists of each MeshBlock
//! written contiguously in order of "gid" in binary format.

#include <sys/stat.h>  // mkdir

#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> // min

#include "athena.hpp"
#include "file_sharding.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

bool OutputIOStatsEnabled() {
  const char *env = std::getenv("ATHENAK_OUTPUT_IO_STATS");
  return (env != nullptr && env[0] != '\0' && env[0] != '0');
}

void PrintBinaryOutputStats(const char *type, FileShardMode mode, bool slice,
                            bool skipped_empty_shard, double open_time,
                            double header_time, double payload_time, double close_time) {
  if (!OutputIOStatsEnabled()) {
    return;
  }

  double local_times[4] = {open_time, header_time, payload_time, close_time};
  double max_times[4] = {open_time, header_time, payload_time, close_time};
  int local_skipped = skipped_empty_shard ? 1 : 0;
  int skipped_ranks = local_skipped;
#if MPI_PARALLEL_ENABLED
  MPI_Reduce(local_times, max_times, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_skipped, &skipped_ranks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (global_variable::my_rank != 0) {
    return;
  }
#endif

  std::cout << "[output-io] type=" << type
            << " mode=" << ShardDistributionName(mode)
            << " slice=" << (slice ? 1 : 0)
            << " skipped_empty_shard_ranks=" << skipped_ranks
            << " open=" << max_times[0]
            << " header=" << max_times[1]
            << " payload=" << max_times[2]
            << " close=" << max_times[3]
            << std::endl;
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
  mkdir("bin",0775);
  if (op.file_shard_mode != FileShardMode::shared) {
    std::string shard_dir = std::string("bin/") + ShardDirectoryName(op.file_shard_mode);
    mkdir(shard_dir.c_str(), 0775);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBinaryOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in binary format
//   All MeshBlocks are written to the same file.

void MeshBinaryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // check if slicing
  bool bin_slice = (out_params.slice1 || out_params.slice2 || out_params.slice3 ||
                    out_params.gid >= 0);
  auto AdvanceOutputCounters = [&]() {
    out_params.file_number++;
    if (out_params.last_time < 0.0) {
      out_params.last_time = pm->time;
    } else {
      out_params.last_time += out_params.dt;
    }
    pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  };

  // create filename: "bin/file_basename" + "." + "file_id" + "." + XXXXX + ".bin"
  // where XXXXX = 5-digit file_number

  FileShardMode shard_mode = out_params.file_shard_mode;
  std::vector<int> slice_shard_counts;
  bool skip_empty_slice_shard = false;
  if (bin_slice) {
    slice_shard_counts = GatherShardCounts(static_cast<int>(outmbs.size()), shard_mode);
    int slice_shard_total = std::accumulate(slice_shard_counts.begin(),
                                            slice_shard_counts.end(), 0);
    skip_empty_slice_shard = (shard_mode == FileShardMode::per_node &&
                              slice_shard_total == 0);
    if (skip_empty_slice_shard) {
      PrintBinaryOutputStats("bin", shard_mode, bin_slice, true, 0.0, 0.0, 0.0, 0.0);
      AdvanceOutputCounters();
      return;
    }
  }

  double open_time = 0.0;
  double header_time = 0.0;
  double payload_time = 0.0;
  double close_time = 0.0;

  std::string fname;
  char number[7];
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  if (shard_mode != FileShardMode::shared) {
    fname = std::string("bin/") + ShardDirectoryName(shard_mode) +
            out_params.file_basename + "." + out_params.file_id + number +
            ".bin";
  } else {
    // Existing behavior: single restart file
    fname = std::string("bin/") + out_params.file_basename + "." +
            out_params.file_id + number + ".bin";
  }

  IOWrapper binfile;
  std::size_t header_offset = 0;
#if MPI_PARALLEL_ENABLED
  if (shard_mode == FileShardMode::per_node) {
    binfile.SetCommunicator(global_variable::node_comm);
  }
#endif
  bool use_serial_io = UsesSerialIO(shard_mode);
  Kokkos::Timer phase_timer;
  binfile.Open(fname.c_str(), IOWrapper::FileMode::write, use_serial_io);
  open_time = phase_timer.seconds();

  // Basic parts of the format:
  // 1. Size of the header
  // 2. Current time
  // 3. List of variables in the file
  // 4. Header (input file information)
  phase_timer.reset();
  {
    std::stringstream msg;
    const int time_precision = std::numeric_limits<Real>::max_digits10 - 1;
    msg << "Athena binary output version=1.1" << std::endl
        // preheader size includes "size of preheader" line up to "number of variables"
        << "  size of preheader=5" << std::endl
        << std::scientific << std::setprecision(time_precision)
        << "  time=" << pm->time << std::endl
        << "  cycle=" << pm->ncycle << std::endl
        << "  size of location=" << sizeof(Real) << std::endl
        << "  size of variable=" << sizeof(float) << std::endl
        << "  number of variables=" << outvars.size() << std::endl
        << "  variables:  ";
    for (int n=0; n<outvars.size(); n++) {
      msg << outvars[n].label.c_str() << "  ";
    }
    msg << std::endl;
    if (IsShardWriter(shard_mode)) {
      binfile.Write_any_type(msg.str().c_str(),msg.str().size(),"byte",
                             use_serial_io);
    }
    header_offset += msg.str().size();
  }
  {
    std::stringstream msg;
    // prepare the input parameters
    std::stringstream ost;
    pin->ParameterDump(ost);
    std::string sbuf=ost.str();
    msg << "  header offset=" << sbuf.size()*sizeof(char)  << std::endl;
    if (IsShardWriter(shard_mode)) {
      binfile.Write_any_type(msg.str().c_str(),msg.str().size(),"byte",
                             use_serial_io);
      binfile.Write_any_type(sbuf.c_str(),sbuf.size(),"byte", use_serial_io);
    }
    header_offset += sbuf.size()*sizeof(char);
    header_offset += msg.str().size();
  }
  header_time = phase_timer.seconds();

  //  5. Data.  An arbitrary number of scalars and vectors can be written (every element
  //  of the outvars vector), all in binary floats format

  int nout_vars = outvars.size();
  int nout_mbs = outmbs.size();
  int cells = 0;
  if (nout_mbs > 0) {
    int nout1 = outmbs[0].oie - outmbs[0].ois + 1;
    int nout2 = outmbs[0].oje - outmbs[0].ojs + 1;
    int nout3 = outmbs[0].oke - outmbs[0].oks + 1;
    cells = nout1*nout2*nout3;
  }

  // ois, oie, ojs, oje, oks, oke + il1, il2, il3, level +
  // x1min, x1max, x2min, x2max, x3min, x3max + data
  std::size_t data_size = 10*sizeof(int32_t) + 6*sizeof(Real)
                        + (cells*nout_vars)*sizeof(float);

  int nb_mbs = pm->nmb_eachrank[global_variable::my_rank];

  int write_mbs = bin_slice ? nout_mbs : nb_mbs;
  // allocate 1D vector of floats used to convert and output data
  char *data = new char[write_mbs*data_size];
  float *single_data = new float[cells];

  // Loop over MeshBlocks
  for (int m=0; m<nout_mbs; ++m) {
    char *pdata=&(data[m*data_size]);
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

  // now write binary data
  phase_timer.reset();
  if (bin_slice) {
    int shard_prefix = PrefixCountBeforeMe(slice_shard_counts, shard_mode);
    std::size_t myoffset = header_offset + data_size*shard_prefix;
    binfile.Write_any_type_at_all(data,(data_size*nout_mbs),myoffset,"byte",
                                  use_serial_io);
  } else {
    std::vector<int> shard_counts = GatherShardCounts(nb_mbs, shard_mode);
    int shard_prefix = PrefixCountBeforeMe(shard_counts, shard_mode);
    int shard_min = *std::min_element(shard_counts.begin(), shard_counts.end());
    int shard_max = *std::max_element(shard_counts.begin(), shard_counts.end());
    // check if elements larger than 2^31
    if (data_size*nb_mbs<=2147483648) {
      // now write binary data in parallel
      std::size_t myoffset = header_offset + data_size*shard_prefix;
      binfile.Write_any_type_at_all(data,(data_size*nb_mbs),myoffset,"byte",
                                    use_serial_io);
    } else {
      // write data over each MeshBlock sequentially and in parallel
      // calculate max/min number of MeshBlocks across all ranks
      for (int m=0;  m<shard_max; ++m) {
        std::size_t myoffset = header_offset + data_size*(shard_prefix + m);
        // every rank has a MB to write, so write collectively
        if (m < shard_min) {
          char *pdata = &(data[m*data_size]);
          if (binfile.Write_any_type_at_all(pdata,(data_size),myoffset,"byte",
                                            use_serial_io) != data_size) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "binary data not written correctly to binary file, "
                << "binary file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
        // some ranks are finished writing, so use non-collective write
        } else if (m < pm->nmb_thisrank) {
          char *pdata = &(data[m*data_size]);
          if (binfile.Write_any_type_at(pdata,(data_size),myoffset,"byte",
                                        use_serial_io) != data_size) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                 << std::endl << "binary data not written correctly to binary file, "
                 << "binary file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }
  payload_time = phase_timer.seconds();

  // close the output file and clean up ptrs to data
  phase_timer.reset();
  binfile.Close(use_serial_io);
  close_time = phase_timer.seconds();
  delete [] data;
  delete [] single_data;

  // increment counters
  AdvanceOutputCounters();

  PrintBinaryOutputStats("bin", shard_mode, bin_slice, false, open_time, header_time,
                         payload_time, close_time);

  return;
}
