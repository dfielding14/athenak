//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <chrono>  // NOLINT(build/c++11)
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility> // make_pair
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "srcterms/turb_driver.hpp"
#include "outputs.hpp"
#include "restart_layout.hpp"
#include "restart_manifest.hpp"

namespace {

constexpr std::size_t kMaxGeneratedPayloadPathBytes = 1024;
constexpr IOWrapperSizeT kMaxNodeRestartManifestBytes = 64ULL*1024ULL*1024ULL;

bool IsSafeRestartLeaf(const std::string &leaf) {
  return !leaf.empty() && leaf != "." && leaf != ".." &&
      leaf.find('/') == std::string::npos &&
      leaf.find('\\') == std::string::npos;
}

[[noreturn]] void FailNodeRestartWrite(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(EXIT_FAILURE);
}

std::string NodeRestartRelativePayloadPath(const std::string &payload_name, int node_id) {
  std::string relative_path = ShardDirectoryName(FileShardMode::node, 0, node_id)
      + "/" + payload_name;
  if (relative_path.size() > kMaxGeneratedPayloadPathBytes) {
    FailNodeRestartWrite("node restart generated payload path exceeds the "
                         "1024-byte limit.");
  }
  return relative_path;
}

IOWrapperSizeT CheckedRestartAdd(IOWrapperSizeT left, IOWrapperSizeT right,
                                 const std::string &context) {
  return restart_layout::CheckedAdd(left, right, FailNodeRestartWrite, context);
}

IOWrapperSizeT CheckedRestartMultiply(IOWrapperSizeT left, IOWrapperSizeT right,
                                      const std::string &context) {
  return restart_layout::CheckedMultiply(left, right, FailNodeRestartWrite, context);
}

IOWrapperSizeT CheckedRestartCount(int value, const std::string &context) {
  return restart_layout::CheckedNonNegative(value, FailNodeRestartWrite, context);
}

int CheckedRestartInt(IOWrapperSizeT value, const std::string &context) {
  if (value > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    FailNodeRestartWrite(context + " exceeds INT_MAX.");
  }
  return static_cast<int>(value);
}

int CheckedRestartIntAdd(int left, int right, const std::string &context) {
  return CheckedRestartInt(
      CheckedRestartAdd(CheckedRestartCount(left, context),
                        CheckedRestartCount(right, context), context),
      context);
}

void CheckedRestartWrite(IOWrapper &file, const void *data, IOWrapperSizeT bytes,
                         const std::string &context, bool independent_file) {
  if (file.Write_any_type(data, bytes, "byte", independent_file) != bytes) {
    FailNodeRestartWrite(context + " was not written completely.");
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  mkdir("rst",0775);
  if (IsSharded(op.shard_mode)) {
    std::string shard_dir = "rst/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    mkdir(shard_dir.c_str(), 0775);
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx1, indcs.ng, true, FailNodeRestartWrite, "restart x1 extent"),
      "restart x1 extent");
  int nout2 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx2, indcs.ng, indcs.nx2 > 1, FailNodeRestartWrite, "restart x2 extent"),
      "restart x2 extent");
  int nout3 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx3, indcs.ng, indcs.nx3 > 1, FailNodeRestartWrite, "restart x3 extent"),
      "restart x3 extent");
  int nmb = pm->pmb_pack->nmb_thispack;

  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nadm=0, nz4c=0;
  if (phydro != nullptr) {
    nhydro = CheckedRestartIntAdd(phydro->nhydro, phydro->nscalars,
                                  "restart Hydro component count");
  }
  if (pmhd != nullptr) {
    nmhd = CheckedRestartIntAdd(pmhd->nmhd, pmhd->nscalars,
                                "restart MHD component count");
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  // if the spacetime is evolved, we do not need to checkpoint/recover the ADM variables
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  restart_layout::PayloadLayout layout = restart_layout::PayloadLayout::Build(
      {CheckedRestartCount(nout1, "restart x1 extent"),
       CheckedRestartCount(nout2, "restart x2 extent"),
       CheckedRestartCount(nout3, "restart x3 extent"),
       CheckedRestartCount(nhydro, "restart Hydro component count"),
       CheckedRestartCount(nmhd, "restart MHD component count"),
       CheckedRestartCount(nrad, "restart radiation component count"),
       pturb == nullptr ? 0 : CheckedRestartCount(nforce, "restart forcing count"),
       CheckedRestartCount(nz4c, "restart Z4c component count"),
       CheckedRestartCount(nadm, "restart ADM component count"),
       sizeof(Real)}, FailNodeRestartWrite);
  IOWrapperSizeT local_blocks = CheckedRestartCount(nmb, "restart local MeshBlock count");
  const auto preflight_local = [&](IOWrapperSizeT block_bytes,
                                   const std::string &context) {
    restart_layout::CheckedMemorySize(
        CheckedRestartMultiply(local_blocks, block_bytes, context),
        FailNodeRestartWrite, context);
  };
  int nout1f = nout1;
  int nout2f = nout2;
  int nout3f = nout3;
  if (pmhd != nullptr) {
    nout1f = CheckedRestartInt(CheckedRestartAdd(nout1, 1, "restart MHD x1 extent"),
                               "restart MHD x1 extent");
    nout2f = CheckedRestartInt(CheckedRestartAdd(nout2, 1, "restart MHD x2 extent"),
                               "restart MHD x2 extent");
    nout3f = CheckedRestartInt(CheckedRestartAdd(nout3, 1, "restart MHD x3 extent"),
                               "restart MHD x3 extent");
  }

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  if (phydro != nullptr) {
    preflight_local(layout.hydro_bytes, "restart local Hydro allocation bytes");
    Kokkos::realloc(outarray_hyd, nmb, nhydro, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    preflight_local(layout.mhd_bytes, "restart local MHD allocation bytes");
    preflight_local(layout.mhd_x1f_bytes, "restart local MHD x1-face allocation bytes");
    preflight_local(layout.mhd_x2f_bytes, "restart local MHD x2-face allocation bytes");
    preflight_local(layout.mhd_x3f_bytes, "restart local MHD x3-face allocation bytes");
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1f);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2f, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x3f, nmb, nout3f, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (prad != nullptr) {
    preflight_local(layout.radiation_bytes, "restart local radiation allocation bytes");
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pturb != nullptr) {
    preflight_local(layout.forcing_bytes, "restart local forcing allocation bytes");
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    preflight_local(layout.z4c_bytes, "restart local Z4c allocation bytes");
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    preflight_local(layout.adm_bytes, "restart local ADM allocation bytes");
    Kokkos::realloc(outarray_adm, nmb, nadm, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_adm, Kokkos::subview(padm->u_adm, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }

  // calculate max/min number of MeshBlocks across all ranks
  noutmbs_max = pm->nmb_eachrank[0];
  noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes everything to a single restart file

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx1, indcs.ng, true, FailNodeRestartWrite, "restart x1 extent"),
      "restart x1 extent");
  int nout2 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx2, indcs.ng, indcs.nx2 > 1, FailNodeRestartWrite, "restart x2 extent"),
      "restart x2 extent");
  int nout3 = CheckedRestartInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx3, indcs.ng, indcs.nx3 > 1, FailNodeRestartWrite, "restart x3 extent"),
      "restart x3 extent");
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nz4c=0, nadm=0;
  IOWrapperSizeT nco=0;
  if (phydro != nullptr) {
    nhydro = CheckedRestartIntAdd(phydro->nhydro, phydro->nscalars,
                                  "restart Hydro component count");
  }
  if (pmhd != nullptr) {
    nmhd = CheckedRestartIntAdd(pmhd->nmhd, pmhd->nscalars,
                                "restart MHD component count");
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
    nco = restart_layout::CheckedSizeT(pz4c->ptracker.size(), FailNodeRestartWrite,
                                       "restart compact-object tracker count");
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  FileShardMode shard_mode = out_params.shard_mode;
  bool independent_file = UsesIndependentFileIO(shard_mode);
  bool node_sharded = IsNodeSharded(shard_mode);
  bool shard_writer = IsRankSharded(shard_mode) ||
      (node_sharded && global_variable::node_rank == 0) ||
      (shard_mode == FileShardMode::shared && global_variable::my_rank == 0);
  std::string fname;
  std::string manifest_name;
  std::string payload_name;
  std::uint64_t generation = 0;
  char number[7];
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  if (IsRankSharded(shard_mode)) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    fname = std::string("rst/") + ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id) + "/"
      + out_params.file_basename
      + number + ".rst";
  } else if (node_sharded) {
    if (!IsSafeRestartLeaf(out_params.file_basename)) {
      FailNodeRestartWrite("Node restart output requires <job>/basename to be a "
                           "single safe path component.");
    }
    if (global_variable::my_rank == 0) {
      generation = static_cast<std::uint64_t>(
          std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&generation, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
#endif
    payload_name = out_params.file_basename + number + ".g"
        + std::to_string(generation) + ".payload.rst";
    for (int id = 0; id < global_variable::nnodes; ++id) {
      NodeRestartRelativePayloadPath(payload_name, id);
    }
    std::string shard_dir = ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id);
    fname = std::string("rst/") + shard_dir + "/" + payload_name + ".tmp";
    manifest_name = std::string("rst/") + out_params.file_basename + number + ".rst";
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
#if MPI_PARALLEL_ENABLED
  if (node_sharded) {
    resfile.SetCommunicator(global_variable::node_comm);
  }
#endif
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write, independent_file);
  if (shard_writer) {
    // output the input parameters (input file)
    CheckedRestartWrite(resfile, sbuf.c_str(), sbuf.size(), "restart parameter dump",
                        independent_file);
    if (node_sharded) {
      CheckedRestartWrite(resfile, kNodeRestartPayloadMarker,
                          kNodeRestartPayloadMarkerSize,
                          "node restart payload marker", independent_file);
    }

    // output Mesh information
    CheckedRestartWrite(resfile, &(pm->nmb_total), sizeof(int),
                        "restart total MeshBlock count", independent_file);
    CheckedRestartWrite(resfile, &(pm->root_level), sizeof(int),
                        "restart root level", independent_file);
    CheckedRestartWrite(resfile, &(pm->mesh_size), sizeof(RegionSize),
                        "restart mesh size", independent_file);
    CheckedRestartWrite(resfile, &(pm->mesh_indcs), sizeof(RegionIndcs),
                        "restart mesh indices", independent_file);
    CheckedRestartWrite(resfile, &(pm->mb_indcs), sizeof(RegionIndcs),
                        "restart MeshBlock indices", independent_file);
    CheckedRestartWrite(resfile, &(pm->time), sizeof(Real),
                        "restart time", independent_file);
    CheckedRestartWrite(resfile, &(pm->dt), sizeof(Real),
                        "restart timestep", independent_file);
    CheckedRestartWrite(resfile, &(pm->ncycle), sizeof(int),
                        "restart cycle", independent_file);
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (shard_writer) {
    CheckedRestartWrite(resfile, &(pm->lloc_eachmb[0]),
                        CheckedRestartMultiply(
                            CheckedRestartCount(pm->nmb_total,
                                                "restart total MeshBlock count"),
                            sizeof(LogicalLocation),
                                               "restart logical-location bytes"),
                        "restart logical locations", independent_file);
    CheckedRestartWrite(resfile, &(pm->cost_eachmb[0]),
                        CheckedRestartMultiply(
                            CheckedRestartCount(pm->nmb_total,
                                                "restart total MeshBlock count"),
                            sizeof(float),
                                               "restart MeshBlock-cost bytes"),
                        "restart MeshBlock costs", independent_file);
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (shard_writer) {
    // store z4c information
    if (pz4c != nullptr) {
      CheckedRestartWrite(resfile, &(pz4c->last_output_time), sizeof(Real),
                          "restart z4c output time", independent_file);
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        CheckedRestartWrite(resfile, pt->GetPos(), 3*sizeof(Real),
                            "restart puncture position", independent_file);
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      CheckedRestartWrite(resfile, &(pturb->rstate), sizeof(RNG_State),
                          "restart turbulence RNG state", independent_file);
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  restart_layout::PayloadLayout layout = restart_layout::PayloadLayout::Build(
      {CheckedRestartCount(nout1, "restart x1 extent"),
       CheckedRestartCount(nout2, "restart x2 extent"),
       CheckedRestartCount(nout3, "restart x3 extent"),
       CheckedRestartCount(nhydro, "restart Hydro component count"),
       CheckedRestartCount(nmhd, "restart MHD component count"),
       CheckedRestartCount(nrad, "restart radiation component count"),
       pturb == nullptr ? 0 : CheckedRestartCount(nforce, "restart forcing count"),
       CheckedRestartCount(nz4c, "restart Z4c component count"),
       CheckedRestartCount(nadm, "restart ADM component count"),
       sizeof(Real)}, FailNodeRestartWrite);
  IOWrapperSizeT data_size = layout.data_bytes;
  if (shard_writer) {
    CheckedRestartWrite(resfile, &(data_size), sizeof(IOWrapperSizeT),
                        "restart per-MeshBlock byte count", independent_file);
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT fixed_mesh_bytes = 3*sizeof(int) + 2*sizeof(Real)
      + sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  IOWrapperSizeT state_bytes = CheckedRestartMultiply(
      CheckedRestartMultiply(3, nco, "restart tracker state bytes"), sizeof(Real),
      "restart tracker state bytes");
  if (pz4c != nullptr) {
    state_bytes = CheckedRestartAdd(state_bytes, sizeof(Real),
                                    "restart Z4c state bytes");
  }
  if (pturb != nullptr) {
    state_bytes = CheckedRestartAdd(state_bytes, sizeof(RNG_State),
                                    "restart turbulence state bytes");
  }
  restart_layout::HeaderLayout header_layout = restart_layout::HeaderLayout::Build(
      {restart_layout::CheckedSizeT(sbuf.size(), FailNodeRestartWrite,
                                    "restart parameter dump bytes"),
       node_sharded ? kNodeRestartPayloadMarkerSize : 0,
       fixed_mesh_bytes,
       restart_layout::CheckedPositive(pm->nmb_total, FailNodeRestartWrite,
                                       "restart total MeshBlock count"),
       sizeof(LogicalLocation) + sizeof(float),
       state_bytes,
       sizeof(IOWrapperSizeT)}, FailNodeRestartWrite);

  // write cell-centered variables in parallel
  IOWrapperSizeT header_size = header_layout.header_bytes;
  IOWrapperSizeT offset_myrank = header_size;
  int payload_block_offset = 0;
  if (shard_mode == FileShardMode::shared) {
    offset_myrank = restart_layout::CheckedOffset(
        header_size, data_size,
        CheckedRestartCount(pm->gids_eachrank[global_variable::my_rank],
                            "restart shared rank block offset"),
        FailNodeRestartWrite, "restart shared rank byte offset");
  } else if (node_sharded) {
    payload_block_offset = global_variable::NodePrefixSum(pm->nmb_thisrank);
    offset_myrank = CheckedRestartAdd(
        offset_myrank,
        CheckedRestartMultiply(data_size,
                               CheckedRestartCount(payload_block_offset,
                                                   "node restart payload block offset"),
                               "node restart rank payload offset"),
        "node restart rank payload offset");
    noutmbs_min = global_variable::NodeMin(pm->nmb_thisrank);
    noutmbs_max = global_variable::NodeMax(pm->nmb_thisrank);
  }

  IOWrapperSizeT myoffset = offset_myrank;

  // write cell-centered variables, one MeshBlock at a time (but parallelized over all
  // ranks). MeshBlocks are written seperately to reduce number of data elements per write
  // call, to avoid exceeding 2^31 limit for very large grids per MPI rank.
  if (phydro != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart Hydro subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart Hydro subview count");
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.hydro_bytes,
                                      "restart field offset"); // hydro u0
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart MHD subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart MHD subview count");
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.mhd_bytes,
                                      "restart field offset"); // mhd u0
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        IOWrapperSizeT fldcnt = restart_layout::CheckedSizeT(
            x1fptr.size(), FailNodeRestartWrite, "restart MHD x1-face subview count");
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x1f_bytes,
                                     "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestartWrite, "restart MHD x2-face subview count");
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x2f_bytes,
                                     "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestartWrite, "restart MHD x3-face subview count");
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x3f_bytes,
                                     "restart MHD face offset");

        myoffset = CheckedRestartAdd(
            myoffset,
            layout.mhd_face_stride_remainder_bytes,
            "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        IOWrapperSizeT fldcnt = restart_layout::CheckedSizeT(
            x1fptr.size(), FailNodeRestartWrite, "restart MHD x1-face subview count");
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x1f_bytes,
                                     "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestartWrite, "restart MHD x2-face subview count");
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x2f_bytes,
                                     "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestartWrite, "restart MHD x3-face subview count");
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      independent_file) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x3f_bytes,
                                     "restart MHD face offset");

        myoffset = CheckedRestartAdd(
            myoffset,
            layout.mhd_face_stride_remainder_bytes,
            "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.mhd_x1f_bytes,
                                      "restart field offset"); // mhd b0.x1f
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.mhd_x2f_bytes,
                                      "restart field offset"); // mhd b0.x2f
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.mhd_x3f_bytes,
                                      "restart field offset"); // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart radiation subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart radiation subview count");
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.radiation_bytes,
                                      "restart field offset"); // radiation i0
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart forcing subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart forcing subview count");
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.forcing_bytes,
                                      "restart field offset"); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart Z4c subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart Z4c subview count");
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.z4c_bytes,
                                      "restart field offset"); // z4c u0
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart ADM subview count");
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestartWrite, "restart ADM subview count");
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      independent_file) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset = CheckedRestartAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    offset_myrank = CheckedRestartAdd(offset_myrank, layout.adm_bytes,
                                      "restart field offset"); // adm u_adm
    myoffset = offset_myrank;
  }

  // close file, clean up
  if (resfile.Close(independent_file) != 0) {
    FailNodeRestartWrite("restart payload could not be closed cleanly.");
  }

  if (node_sharded) {
    int node_blocks = global_variable::NodeSum(pm->nmb_thisrank);
    IOWrapperSizeT expected_size = CheckedRestartAdd(
        header_size,
        CheckedRestartMultiply(data_size,
                               CheckedRestartCount(node_blocks,
                                                   "node restart payload block count"),
                               "node restart payload size"),
        "node restart payload size");
    std::string completed_payload = fname.substr(0, fname.size() - 4);
    if (global_variable::node_rank == 0) {
      std::ifstream payload_check(fname, std::ios::binary | std::ios::ate);
      IOWrapperSizeT observed_size = payload_check.good()
          ? static_cast<IOWrapperSizeT>(payload_check.tellg()) : 0;
      payload_check.close();
      if (observed_size != expected_size ||
          std::rename(fname.c_str(), completed_payload.c_str()) != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Node restart payload '" << fname
                  << "' was not completed atomically; expected " << expected_size
                  << " bytes and found " << observed_size << "." << std::endl;
#if MPI_PARALLEL_ENABLED
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        std::exit(EXIT_FAILURE);
      }
    }

#if MPI_PARALLEL_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    std::vector<int> manifest_nodes;
    std::vector<int> manifest_offsets;
    if (global_variable::my_rank == 0) {
      manifest_nodes.resize(global_variable::nranks);
      manifest_offsets.resize(global_variable::nranks);
    }
#if MPI_PARALLEL_ENABLED
    MPI_Gather(&(global_variable::node_id), 1, MPI_INT, manifest_nodes.data(), 1,
               MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&payload_block_offset, 1, MPI_INT, manifest_offsets.data(), 1,
               MPI_INT, 0, MPI_COMM_WORLD);
#else
    manifest_nodes[0] = global_variable::node_id;
    manifest_offsets[0] = payload_block_offset;
#endif

    if (global_variable::my_rank == 0) {
      std::vector<int> blocks_per_node(global_variable::nnodes, 0);
      std::vector<int> next_payload_block(global_variable::nnodes, 0);
      int next_gid = 0;
      for (int r = 0; r < global_variable::nranks; ++r) {
        int id = manifest_nodes[r];
        int blocks = pm->nmb_eachrank[r];
        if (id < 0 || id >= global_variable::nnodes || blocks < 0 ||
            pm->gids_eachrank[r] != next_gid ||
            manifest_offsets[r] != next_payload_block[id]) {
          FailNodeRestartWrite("Node restart segment map is inconsistent and cannot "
                               "be published.");
        }
        blocks_per_node[id] = CheckedRestartIntAdd(
            blocks_per_node[id], blocks, "node restart blocks-per-node count");
        next_payload_block[id] = CheckedRestartIntAdd(
            next_payload_block[id], blocks, "node restart payload block offset");
        next_gid = CheckedRestartIntAdd(next_gid, blocks,
                                        "node restart global block count");
      }
      if (next_gid != pm->nmb_total) {
        FailNodeRestartWrite("Node restart segment map does not cover all mesh blocks.");
      }
      std::string temporary_manifest = manifest_name + ".tmp.g"
          + std::to_string(generation);
      std::ofstream manifest(temporary_manifest, std::ios::trunc);
      restart_layout::ManifestBudget manifest_budget{0, kMaxNodeRestartManifestBytes};
      const auto write_manifest_record = [&](const std::string &record) {
        manifest_budget.Add(CheckedRestartAdd(
            restart_layout::CheckedSizeT(record.size(), FailNodeRestartWrite,
                                         "node restart manifest record bytes"),
            1, "node restart manifest record bytes"), FailNodeRestartWrite);
        manifest << record << "\n";
      };
      write_manifest_record("AthenaK node restart manifest version=1");
      write_manifest_record("complete=1");
      write_manifest_record("payload_count=" + std::to_string(global_variable::nnodes));
      write_manifest_record("nmb_total=" + std::to_string(pm->nmb_total));
      write_manifest_record("header_size=" + std::to_string(header_size));
      write_manifest_record("data_size=" + std::to_string(data_size));
      for (int id = 0; id < global_variable::nnodes; ++id) {
        std::string relative_path = NodeRestartRelativePayloadPath(payload_name, id);
        IOWrapperSizeT payload_size = CheckedRestartAdd(
            header_size,
            CheckedRestartMultiply(data_size,
                                   CheckedRestartCount(blocks_per_node[id],
                                                       "node restart payload block "
                                                       "count"),
                                   "node restart manifest payload size"),
            "node restart manifest payload size");
        write_manifest_record("payload " + std::to_string(id) + " " +
                              std::to_string(blocks_per_node[id]) + " " +
                              std::to_string(payload_size) + " " + relative_path);
      }
      for (int r = 0; r < global_variable::nranks; ++r) {
        write_manifest_record("segment " + std::to_string(manifest_nodes[r]) + " " +
                              std::to_string(pm->gids_eachrank[r]) + " " +
                              std::to_string(pm->nmb_eachrank[r]) + " " +
                              std::to_string(manifest_offsets[r]));
      }
      write_manifest_record("end");
      std::streampos observed_manifest_bytes = manifest.tellp();
      manifest.close();
      if (!manifest.good() || observed_manifest_bytes < 0 ||
          static_cast<IOWrapperSizeT>(observed_manifest_bytes) != manifest_budget.bytes ||
          std::rename(temporary_manifest.c_str(), manifest_name.c_str()) != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Node restart manifest '" << manifest_name
                  << "' could not be published atomically." << std::endl;
#if MPI_PARALLEL_ENABLED
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        std::exit(EXIT_FAILURE);
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  return;
}
