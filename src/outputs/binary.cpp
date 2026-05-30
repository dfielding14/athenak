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
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// Constructor: also calls BaseTypeOutput base class constructor

MeshBinaryOutput::MeshBinaryOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs
  // useful for mpiio-based outputs because on some supercomputers you may need to
  // set different stripe counts depending on whether mpiio is used in order to
  // achieve the best performance and not to crash the filesystem
  mkdir("bin",0775);
  if (IsSharded(op.shard_mode)) {
    std::string shard_dir = "bin/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    mkdir(shard_dir.c_str(), 0775);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBinaryOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in binary format
//   All MeshBlocks are written to the same file.

void MeshBinaryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // check if slicing
  bool bin_slice = (out_params.slice1 || out_params.slice2 || out_params.slice3);
  FileShardMode shard_mode = out_params.shard_mode;
  bool independent_file = UsesIndependentFileIO(shard_mode);
  bool shard_writer = IsRankSharded(shard_mode) ||
      (IsNodeSharded(shard_mode) && global_variable::node_rank == 0) ||
      (shard_mode == FileShardMode::shared && global_variable::my_rank == 0);

  // create filename: "bin/file_basename" + "." + "file_id" + "." + XXXXX + ".bin"
  // where XXXXX = 5-digit file_number

  std::string fname;
  if (IsSharded(shard_mode)) {
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("bin/") + ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id) + "/"
          + out_params.file_basename
          + "." + out_params.file_id + number + ".bin";
  } else {
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("bin/") + out_params.file_basename
          + "." + out_params.file_id + number + ".bin";
  }

  IOWrapper binfile;
  std::size_t header_offset=0;
#if MPI_PARALLEL_ENABLED
  if (IsNodeSharded(shard_mode)) {
    binfile.SetCommunicator(global_variable::node_comm);
  }
#endif
  binfile.Open(fname.c_str(), IOWrapper::FileMode::write, independent_file);

  int nout_mbs = outmbs.size();
  int shard_nout_mbs = nout_mbs;
  if (IsNodeSharded(shard_mode)) {
    shard_nout_mbs = global_variable::NodeSum(nout_mbs);
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
        << "  size of preheader=" << (IsNodeSharded(shard_mode) ? 8 : 5) << std::endl
        << std::scientific << std::setprecision(time_precision)
        << "  time=" << pm->time << std::endl
        << "  cycle=" << pm->ncycle << std::endl
        << "  size of location=" << sizeof(Real) << std::endl
        << "  size of variable=" << sizeof(float) << std::endl;
    if (IsNodeSharded(shard_mode)) {
      msg << "  distribution=node" << std::endl
          << "  node=" << global_variable::node_id << std::endl
          << "  number of meshblocks=" << shard_nout_mbs << std::endl;
    }
    msg << "  number of variables=" << outvars.size() << std::endl
        << "  variables:  ";
    for (int n=0; n<outvars.size(); n++) {
      msg << outvars[n].label.c_str() << "  ";
    }
    msg << std::endl;
    if (shard_writer) {
      binfile.Write_any_type(msg.str().c_str(),msg.str().size(),"byte",
                             independent_file);
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
    if (shard_writer) {
      binfile.Write_any_type(msg.str().c_str(),msg.str().size(),"byte",
                             independent_file);
      binfile.Write_any_type(sbuf.c_str(),sbuf.size(),"byte", independent_file);
    }
    header_offset += sbuf.size()*sizeof(char);
    header_offset += msg.str().size();
  }

  //  5. Data.  An arbitrary number of scalars and vectors can be written (every element
  //  of the outvars vector), all in binary floats format

  int nout_vars = outvars.size();
  int cells = 0;
  if (nout_mbs > 0) {
    int nout1 = outmbs[0].oie - outmbs[0].ois + 1;
    int nout2 = outmbs[0].oje - outmbs[0].ojs + 1;
    int nout3 = outmbs[0].oke - outmbs[0].oks + 1;
    cells = nout1*nout2*nout3;
  }
#if MPI_PARALLEL_ENABLED
  if (bin_slice && !IsRankSharded(shard_mode)) {
    int shard_cells = 0;
    MPI_Comm comm = IsNodeSharded(shard_mode) ? global_variable::node_comm : MPI_COMM_WORLD;
    MPI_Allreduce(&cells, &shard_cells, 1, MPI_INT, MPI_MAX, comm);
    cells = shard_cells;
  }
#endif

  // ois, oie, ojs, oje, oks, oke + il1, il2, il3, level +
  // x1min, x1max, x2min, x2max, x3min, x3max + data
  std::size_t data_size = 10*sizeof(int32_t) + 6*sizeof(Real)
                        + (cells*nout_vars)*sizeof(float);

  int ns_mbs = pm->gids_eachrank[global_variable::my_rank];
  int nb_mbs = pm->nmb_eachrank[global_variable::my_rank];
  int node_offset = IsNodeSharded(shard_mode) ? global_variable::NodePrefixSum(
      bin_slice ? nout_mbs : nb_mbs) : 0;

  // allocate 1D vector of floats used to convert and output data
  char *data = new char[nb_mbs*data_size];
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
  if (bin_slice) {
    std::vector<int> rank_offset(global_variable::nranks, 0);
    std::partial_sum(noutmbs.begin(),std::prev(noutmbs.end()),
                     std::next(rank_offset.begin()));
    std::size_t myoffset = header_offset+data_size*rank_offset[global_variable::my_rank];

    if (IsRankSharded(shard_mode)) {
      myoffset = header_offset;  // Reset offset for individual files
    } else if (IsNodeSharded(shard_mode)) {
      myoffset = header_offset + data_size*node_offset;
    }

    if (noutmbs_min > 0) {
      binfile.Write_any_type_at_all(data,(data_size*nout_mbs),myoffset,"byte",
                                    independent_file);
    } else {
      if (nout_mbs > 0) {
        binfile.Write_any_type_at(data,(data_size*nout_mbs),myoffset,"byte",
                                    independent_file);
      }
    }
  } else {
    // check if elements larger than 2^31
    if (data_size*nb_mbs<=2147483648) {
      // now write binary data in parallel
      std::size_t myoffset = header_offset;
      if (shard_mode == FileShardMode::shared) {
        myoffset += data_size*ns_mbs;
      } else if (IsNodeSharded(shard_mode)) {
        myoffset += data_size*node_offset;
      }
      binfile.Write_any_type_at_all(data,(data_size*nb_mbs),myoffset,"byte",
                                    independent_file);
    } else {
      // write data over each MeshBlock sequentially and in parallel
      // calculate max/min number of MeshBlocks across all ranks
      if (IsNodeSharded(shard_mode)) {
        noutmbs_max = global_variable::NodeMax(nb_mbs);
        noutmbs_min = global_variable::NodeMin(nb_mbs);
      } else {
        noutmbs_max = pm->nmb_eachrank[0];
        noutmbs_min = pm->nmb_eachrank[0];
        for (int i=0; i<(global_variable::nranks); ++i) {
          noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
          noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
        }
      }
      for (int m=0;  m<noutmbs_max; ++m) {
        char *pdata=&(data[m*data_size]);
        std::size_t myoffset = header_offset + data_size*m;
        if (shard_mode == FileShardMode::shared) {
          myoffset += data_size*ns_mbs;
        } else if (IsNodeSharded(shard_mode)) {
          myoffset += data_size*node_offset;
        }
        // every rank has a MB to write, so write collectively
        if (m < noutmbs_min) {
            if (binfile.Write_any_type_at_all(pdata,(data_size),myoffset,"byte",
                                              independent_file) != data_size) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "binary data not written correctly to binary file, "
                << "binary file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
        // some ranks are finished writing, so use non-collective write
        } else if (m < pm->nmb_thisrank) {
          if (binfile.Write_any_type_at(pdata,(data_size),myoffset,"byte",
                                          independent_file) != data_size) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                 << std::endl << "binary data not written correctly to binary file, "
                 << "binary file is broken." << std::endl;
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  // close the output file and clean up ptrs to data
  binfile.Close(independent_file);
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
