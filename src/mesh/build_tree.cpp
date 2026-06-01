//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file build_tree.cpp
//! \brief Functions to build MeshBlockTreee, both for new runs and restarts

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cinttypes>
#include <limits> // numeric_limits<>
#include <memory> // make_unique<>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/restart_utils.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

bool SameRegionSize(const RegionSize &lhs, const RegionSize &rhs) {
  return lhs.x1min == rhs.x1min && lhs.x2min == rhs.x2min && lhs.x3min == rhs.x3min &&
         lhs.x1max == rhs.x1max && lhs.x2max == rhs.x2max && lhs.x3max == rhs.x3max &&
         lhs.dx1 == rhs.dx1 && lhs.dx2 == rhs.dx2 && lhs.dx3 == rhs.dx3;
}

bool SameRegionIndcs(const RegionIndcs &lhs, const RegionIndcs &rhs) {
  return lhs.ng == rhs.ng && lhs.nx1 == rhs.nx1 && lhs.nx2 == rhs.nx2 &&
         lhs.nx3 == rhs.nx3 && lhs.is == rhs.is && lhs.ie == rhs.ie &&
         lhs.js == rhs.js && lhs.je == rhs.je && lhs.ks == rhs.ks &&
         lhs.ke == rhs.ke && lhs.cnx1 == rhs.cnx1 && lhs.cnx2 == rhs.cnx2 &&
         lhs.cnx3 == rhs.cnx3 && lhs.cis == rhs.cis && lhs.cie == rhs.cie &&
         lhs.cjs == rhs.cjs && lhs.cje == rhs.cje && lhs.cks == rhs.cks &&
         lhs.cke == rhs.cke;
}

bool SameMeshRegionIndcs(const RegionIndcs &lhs, const RegionIndcs &rhs) {
  return lhs.ng == rhs.ng && lhs.nx1 == rhs.nx1 && lhs.nx2 == rhs.nx2 &&
         lhs.nx3 == rhs.nx3 && lhs.is == rhs.is && lhs.ie == rhs.ie &&
         lhs.js == rhs.js && lhs.je == rhs.je && lhs.ks == rhs.ks &&
         lhs.ke == rhs.ke;
}

bool SameLogicalLocation(const LogicalLocation &lhs, const LogicalLocation &rhs) {
  return lhs.lx1 == rhs.lx1 && lhs.lx2 == rhs.lx2 && lhs.lx3 == rhs.lx3 &&
         lhs.level == rhs.level;
}

bool CheckedMpiByteCount(const IOWrapperSizeT count, const IOWrapperSizeT width,
                         int &bytes) {
  constexpr IOWrapperSizeT max = static_cast<IOWrapperSizeT>(
      std::numeric_limits<int>::max());
  if (count != 0 && width > max/count) return false;
  bytes = static_cast<int>(count*width);
  return true;
}

int CheckedMaxRefinementLevel(ParameterInput *pin, const int root_level) {
  const int num_levels = pin->GetOrAddInteger("mesh_refinement", "num_levels", 1);
  if (root_level < 0 || root_level > 30 || num_levels < 1 ||
      num_levels > 31 - root_level) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Number of refinement levels must be between 1 and "
              << 31 - root_level << " for the configured root grid." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  return root_level + num_levels - 1;
}

int CheckedRootLevel(const int nmbmax) {
  if (nmbmax <= 0 || nmbmax > (1 << 30)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Root-grid MeshBlock count exceeds supported logical bounds."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  int level = 0;
  while ((1 << level) < nmbmax) ++level;
  return level;
}

int CheckedStaticRefinementLevel(const int phy_ref_lev, const int root_level,
                                 const int max_level) {
  if (phy_ref_lev < 1 || phy_ref_lev > 30 - root_level) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Refinement level exceeds supported logical bounds." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const int log_ref_lev = root_level + phy_ref_lev;
  if (log_ref_lev > max_level) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Refinement level exceeds maximum allowed ("
              << max_level << ")" << std::endl << "Reduce/specify 'num_levels' in "
              << "<mesh_refinement> input block if using AMR" << std::endl;
    restart_utils::AbortOnFatalError();
  }
  return log_ref_lev;
}

std::int32_t CheckedRefinedLogicalExtent(const int nmb_root, const int phy_ref_lev) {
  const std::int64_t factor = std::int64_t{1} << phy_ref_lev;
  if (nmb_root <= 0 ||
      static_cast<std::int64_t>(nmb_root) >
          std::numeric_limits<std::int32_t>::max()/factor) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Refined logical MeshBlock extent exceeds supported bounds."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  return static_cast<std::int32_t>(static_cast<std::int64_t>(nmb_root)*factor);
}

std::int32_t FirstRefinedLogicalIndex(const std::int32_t logical_extent,
                                     const Real mesh_min, const Real mesh_max,
                                     const Real region_min) {
  std::int32_t lower = 0;
  std::int32_t upper = logical_extent;
  while (lower < upper) {
    const std::int32_t middle = lower + (upper - lower)/2;
    if (LeftEdgeX(middle + 1, logical_extent, mesh_min, mesh_max) > region_min) {
      upper = middle;
    } else {
      lower = middle + 1;
    }
  }
  return lower;
}

std::int32_t LastRefinedLogicalIndex(const std::int32_t logical_extent,
                                    const std::int32_t first,
                                    const Real mesh_min, const Real mesh_max,
                                    const Real region_max) {
  std::int32_t lower = first;
  std::int32_t upper = logical_extent;
  while (lower < upper) {
    const std::int32_t middle = lower + (upper - lower)/2;
    if (LeftEdgeX(middle + 1, logical_extent, mesh_min, mesh_max) >= region_max) {
      upper = middle;
    } else {
      lower = middle + 1;
    }
  }
  return lower;
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromScratch():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for new runs (starting from scratch), using parameters
//! read from input file.  Also does initial load balance based on simple cost estimate.

void Mesh::BuildTreeFromScratch(ParameterInput *pin) {
  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;

  // find maximum number of MeshBlocks at root level in any dir
  int nmbmax = (nmb_rootx1 > nmb_rootx2) ? nmb_rootx1 : nmb_rootx2;
  nmbmax = (nmbmax > nmb_rootx3) ? nmbmax : nmb_rootx3;

  // Find smallest N such that 2^N >= max number of MeshBlocks in any dimension.
  root_level = CheckedRootLevel(nmbmax);
  int current_level = root_level;

  // Construct tree and create root grid
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    max_level = CheckedMaxRefinementLevel(pin, root_level);
  } else {
    max_level = 30;
  }

  // For meshes with refinement, construct new nodes for <refinement> blocks in input file

  if (multilevel) {
    // error check that number of cells in MeshBlock divisible by two
    if (mb_indcs.nx1 % 2 != 0 ||
       (mb_indcs.nx2 % 2 != 0 && multi_d) ||
       (mb_indcs.nx3 % 2 != 0 && three_d)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of cells in MeshBlock must be divisible by 2 "
                << "with SMR or AMR." << std::endl;
      restart_utils::AbortOnFatalError();
    }

    // cycle through ParameterInput list and find "refinement" blocks (SMR), extract data
    // Expand MeshBlockTree to include "refinement" regions specified in input file:
    for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
      if (it->block_name.compare(0, 10, "refinement") == 0) {
        RegionSize ref_size;
        ref_size.x1min = pin->GetReal(it->block_name, "x1min");
        ref_size.x1max = pin->GetReal(it->block_name, "x1max");
        if (multi_d) {
          ref_size.x2min = pin->GetReal(it->block_name, "x2min");
          ref_size.x2max = pin->GetReal(it->block_name, "x2max");
        } else {
          ref_size.x2min = mesh_size.x2min;
          ref_size.x2max = mesh_size.x2max;
        }
        if (three_d) {
          ref_size.x3min = pin->GetReal(it->block_name, "x3min");
          ref_size.x3max = pin->GetReal(it->block_name, "x3max");
        } else {
          ref_size.x3min = mesh_size.x3min;
          ref_size.x3max = mesh_size.x3max;
        }
        int phy_ref_lev = pin->GetInteger(it->block_name, "level");
        int log_ref_lev =
            CheckedStaticRefinementLevel(phy_ref_lev, root_level, max_level);
        if (log_ref_lev > current_level) current_level = log_ref_lev;

        // error check parameters in "refinement" blocks
        if (   ref_size.x1min > ref_size.x1max
            || ref_size.x2min > ref_size.x2max
            || ref_size.x3min > ref_size.x3max)  {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid refinement region (xmax < xmin in one direction)."
              << std::endl;
          restart_utils::AbortOnFatalError();
        }
        if (   ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
            || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
            || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Refinement region must be fully contained within root mesh"
              << std::endl;
          restart_utils::AbortOnFatalError();
        }

        // note: if following is too slow, it could be replaced with bi-section search.
        // Suppose entire root domain is tiled with MeshBlocks at the desired refinement
        // level. Find range of x1-integer indices of such MeshBlocks that cover the
        // refinement region
        std::int32_t lx1min = 0, lx1max = 0;
        std::int32_t lx2min = 0, lx2max = 0;
        std::int32_t lx3min = 0, lx3max = 0;
        std::int32_t lxmax = CheckedRefinedLogicalExtent(nmb_rootx1, phy_ref_lev);
        lx1min = FirstRefinedLogicalIndex(lxmax, mesh_size.x1min, mesh_size.x1max,
                                          ref_size.x1min);
        lx1max = LastRefinedLogicalIndex(lxmax, lx1min, mesh_size.x1min,
                                        mesh_size.x1max, ref_size.x1max);
        if (lx1min % 2 == 1) lx1min--;
        if (lx1max % 2 == 0) lx1max++;

        // Find range of x2-indices of such MeshBlocks that cover the refinement region
        if (multi_d) { // 2D or 3D
          lxmax = CheckedRefinedLogicalExtent(nmb_rootx2, phy_ref_lev);
          lx2min = FirstRefinedLogicalIndex(lxmax, mesh_size.x2min, mesh_size.x2max,
                                            ref_size.x2min);
          lx2max = LastRefinedLogicalIndex(lxmax, lx2min, mesh_size.x2min,
                                          mesh_size.x2max, ref_size.x2max);
          if (lx2min % 2 == 1) lx2min--;
          if (lx2max % 2 == 0) lx2max++;
        }

        // Find range of x3-indices of such MeshBlocks that cover the refinement region
        if (three_d) { // 3D
          lxmax = CheckedRefinedLogicalExtent(nmb_rootx3, phy_ref_lev);
          lx3min = FirstRefinedLogicalIndex(lxmax, mesh_size.x3min, mesh_size.x3max,
                                            ref_size.x3min);
          lx3max = LastRefinedLogicalIndex(lxmax, lx3min, mesh_size.x3min,
                                          mesh_size.x3max, ref_size.x3max);
          if (lx3min % 2 == 1) lx3min--;
          if (lx3max % 2 == 0) lx3max++;
        }

        // Now add nodes to the MeshBlockTree corresponding to these MeshBlocks
        if (one_d) {  // 1D
          for (std::int32_t i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nlloc;
            nlloc.level = log_ref_lev;
            nlloc.lx1 = i;
            nlloc.lx2 = 0;
            nlloc.lx3 = 0;
            int nnew = 0;
            ptree->AddNode(nlloc, nnew);
          }
        }
        if (two_d) {  // 2D
          for (std::int32_t j=lx2min; j<lx2max; j+=2) {
            for (std::int32_t i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nlloc;
              nlloc.level = log_ref_lev;
              nlloc.lx1 = i;
              nlloc.lx2 = j;
              nlloc.lx3 = 0;
              int nnew = 0;
              ptree->AddNode(nlloc, nnew);
            }
          }
        }
        if (three_d) {  // 3D
          for (std::int32_t k=lx3min; k<lx3max; k+=2) {
            for (std::int32_t j=lx2min; j<lx2max; j+=2) {
              for (std::int32_t i=lx1min; i<lx1max; i+=2) {
                LogicalLocation nlloc;
                nlloc.level = log_ref_lev;
                nlloc.lx1 = i;
                nlloc.lx2 = j;
                nlloc.lx3 = k;
                int nnew = 0;
                ptree->AddNode(nlloc, nnew);
              }
            }
          }
        }
      }
    }
  } // if (multilevel)

  if (!adaptive) max_level = current_level;

  // initial mesh hierarchy construction is completed here
  ptree->CountMeshBlocks(nmb_total);
  if (nmb_total <= 0 || nmb_total > kMaxSupportedMeshBlocks) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "MeshBlock count exceeds supported topology bounds." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  ptree->ResetMeshBlockCount(nmb_total);

  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // following returns LogicalLocation list sorted by Z-ordering, and total # of MBs
  ptree->CreateZOrderedLLList(lloc_eachmb, nullptr, nmb_total);

#if MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (nmb_total < global_variable::nranks) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
    restart_utils::AbortOnFatalError();
  }
#endif

  // initialize cost array with the simplest estimate; all the blocks are equal
  // TODO(@user): implement variable cost per MeshBlock as needed
  for (int i=0; i<nmb_total; i++) {cost_eachmb[i] = 1.0;}
  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  nmb_packs_thisrank = 1;
  pmb_pack->AddMeshBlocks(pin);
  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);

  // Fix maximum number of MeshBlocks per rank with AMR
  nmb_maxperrank = nmb_thisrank;
  if (adaptive) {
    if (pin->DoesParameterExist("mesh_refinement", "max_nmb_per_rank")) {
      nmb_maxperrank = pin->GetReal("mesh_refinement", "max_nmb_per_rank");
      if (nmb_maxperrank < nmb_thisrank) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "On rank=" << global_variable::my_rank << " Root grid requires "
          << "more MeshBlocks (nmb_thisrank=" << nmb_thisrank << ") than specified by "
          << "<mesh_refinement>/max_nmb_per_rank=" << nmb_maxperrank << std::endl;
        restart_utils::AbortOnFatalError();
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  if (nmb_maxperrank > (1 << (NUM_BITS_LID))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "Maximum number of MeshBlocks per rank cannot exceed 2^(NUM_BITS_LID) due to MPI"
      << " tag limits" << std::endl;
    restart_utils::AbortOnFatalError();
  }
#endif

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
  }

  // set initial time/cycle parameters, output diagnostics
  time = pin->GetOrAddReal("time", "start_time", 0.0);
  dt   = std::numeric_limits<float>::max();
  cfl_no = pin->GetReal("time", "cfl_number");
  ncycle = 0;
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromRestart():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for restart runs, using parameters and data read from
//! restart file.

void Mesh::BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
                                                     bool single_file_per_rank) {
  const RegionSize configured_mesh_size = mesh_size;
  const RegionIndcs configured_mesh_indcs = mesh_indcs;
  const RegionIndcs configured_mb_indcs = mb_indcs;

  // At this point, the restartfile is already open and the ParameterInput (input file)
  // data has already been read in main(). Thus the file pointer is set to after <par_end>
  IOWrapperSizeT headeroffset = resfile.GetPosition(single_file_per_rank);

  // following must be identical to calculation of headeroffset (excluding size of
  // ParameterInput data) in restart.cpp
  IOWrapperSizeT headersize = 4*sizeof(int) + 2*sizeof(Real)
    + sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  char *headerdata = new char[headersize];

  // the master process reads the header data if single_file_per_rank is false
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    IOWrapperSizeT read_size = resfile.Read_bytes(headerdata, 1, headersize,
                                                  single_file_per_rank);
    if (read_size != headersize) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Header size read from restart file is incorrect, "
                << "expected " << headersize << ", got " << read_size << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }

#if MPI_PARALLEL_ENABLED
  // then broadcast the header data
  if (!single_file_per_rank) {
    int mpi_err = MPI_Bcast(headerdata, headersize, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (mpi_err != MPI_SUCCESS) {
      char error_string[1024];
      int length_of_error_string;
      MPI_Error_string(mpi_err, error_string, &length_of_error_string);
      std::cout << "MPI_Bcast failed with error: " << error_string << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#endif

  // Now copy mesh data read from restart file into Mesh variables. Order of variables
  // set by Write()'s in restart.cpp
  // Note this overwrites size and indices initialized in Mesh constructor.
  IOWrapperSizeT hdos = 0;
  std::memcpy(&nmb_total, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos += sizeof(RegionSize);
  std::memcpy(&mesh_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&mb_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&(restart_meta.original_nranks), &(headerdata[hdos]), sizeof(int));
  delete [] headerdata;

  if (nmb_total <= 0 || nmb_total > kMaxSupportedMeshBlocks ||
      restart_meta.original_nranks <= 0 ||
      restart_meta.original_nranks > nmb_total) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "MeshBlock or rank count stored in restart file is invalid or exceeds "
              << "the supported restart limit."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!SameRegionSize(mesh_size, configured_mesh_size) ||
      !SameMeshRegionIndcs(mesh_indcs, configured_mesh_indcs) ||
      !SameRegionIndcs(mb_indcs, configured_mb_indcs)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart mesh geometry does not match configured mesh geometry."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (root_level < 0 || root_level > 30 ||
      mesh_indcs.nx1 % mb_indcs.nx1 != 0 ||
      mesh_indcs.nx2 % mb_indcs.nx2 != 0 ||
      mesh_indcs.nx3 % mb_indcs.nx3 != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart mesh geometry or root level is invalid."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const auto checked_add_product = [](IOWrapperSizeT &total, const IOWrapperSizeT count,
                                      const IOWrapperSizeT width) {
    constexpr IOWrapperSizeT max = std::numeric_limits<IOWrapperSizeT>::max();
    if ((count != 0 && width > max/count) || total > max - count*width) return false;
    total += count*width;
    return true;
  };
  IOWrapperSizeT minimum_size = headeroffset;
  const IOWrapperSizeT meshblocks = static_cast<IOWrapperSizeT>(nmb_total);
  const IOWrapperSizeT original_nranks =
      static_cast<IOWrapperSizeT>(restart_meta.original_nranks);
  bool valid_minimum_size = checked_add_product(minimum_size, 1, headersize);
  valid_minimum_size =
      valid_minimum_size &&
      checked_add_product(minimum_size, meshblocks,
                          sizeof(LogicalLocation) + sizeof(float));
  valid_minimum_size =
      valid_minimum_size && checked_add_product(minimum_size, meshblocks, sizeof(int));
  valid_minimum_size =
      valid_minimum_size &&
      checked_add_product(minimum_size, original_nranks, 2*sizeof(int));
  if (!valid_minimum_size || minimum_size > resfile.GetSize(single_file_per_rank)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart artifact is too small for serialized MeshBlock layout."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;
  if (nmb_rootx1 <= 0 || nmb_rootx2 <= 0 || nmb_rootx3 <= 0 ||
      nmb_rootx1 > nmb_total || nmb_rootx2 > nmb_total/nmb_rootx1 ||
      nmb_rootx3 > nmb_total/(nmb_rootx1*nmb_rootx2)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart root MeshBlock grid exceeds serialized MeshBlock count."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const int nmbmax = std::max(nmb_rootx1, std::max(nmb_rootx2, nmb_rootx3));
  int expected_root_level = 0;
  while ((1 << expected_root_level) < nmbmax) ++expected_root_level;
  if (root_level != expected_root_level) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart root level is inconsistent with configured mesh geometry."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  int current_level = root_level;

  int idlist_mpi_bytes = 0;
  int meshblock_int_mpi_bytes = 0;
  int rank_int_mpi_bytes = 0;
  if (!CheckedMpiByteCount(meshblocks, sizeof(LogicalLocation) + sizeof(float),
                           idlist_mpi_bytes) ||
      !CheckedMpiByteCount(meshblocks, sizeof(int), meshblock_int_mpi_bytes) ||
      !CheckedMpiByteCount(original_nranks, sizeof(int), rank_int_mpi_bytes)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart MeshBlock layout exceeds supported MPI broadcast counts."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    max_level = CheckedMaxRefinementLevel(pin, root_level);
  } else {
    max_level = 31;
  }

  // allocate memory for lists read from restart
  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // allocate idlist buffer and read list of logical locations and cost
  IOWrapperSizeT listsize = sizeof(LogicalLocation) + sizeof(float);
  char *idlist = new char[listsize*nmb_total];
  // only the master process reads the ID list
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(idlist,listsize,nmb_total,single_file_per_rank) !=
        static_cast<unsigned int>(nmb_total)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Incorrect number of MeshBlocks in restart file; "
                << "restart file is broken." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the ID list
  if (!single_file_per_rank) {
    MPI_Bcast(idlist, idlist_mpi_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

  // everyone sets the logical location and cost lists based on bradcasted data
  int os = 0;
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(lloc_eachmb[i]), &(idlist[os]), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
  }
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(cost_eachmb[i]), &(idlist[os]), sizeof(float));
    os += sizeof(float);
    if (lloc_eachmb[i].level > current_level) current_level = lloc_eachmb[i].level;
  }
  delete [] idlist;
  double total_cost = 0.0;
  for (int i=0; i<nmb_total; ++i) {
    const LogicalLocation &loc = lloc_eachmb[i];
    const bool valid_level = loc.level >= root_level && loc.level <= 30;
    const int level_offset = valid_level ? loc.level - root_level : 0;
    std::int64_t nx1 = 0;
    std::int64_t nx2 = 0;
    std::int64_t nx3 = 0;
    if (valid_level) {
      const std::int64_t scale = static_cast<std::int64_t>(1) << level_offset;
      nx1 = static_cast<std::int64_t>(nmb_rootx1)*scale;
      nx2 = static_cast<std::int64_t>(nmb_rootx2)*scale;
      nx3 = static_cast<std::int64_t>(nmb_rootx3)*scale;
    }
    const std::int64_t max_logical_extent = std::numeric_limits<std::int32_t>::max();
    if (!valid_level || nx1 > max_logical_extent || nx2 > max_logical_extent ||
        nx3 > max_logical_extent || loc.lx1 < 0 || loc.lx1 >= nx1 ||
        loc.lx2 < 0 || loc.lx2 >= nx2 || loc.lx3 < 0 || loc.lx3 >= nx3 ||
        (!multi_d && loc.lx2 != 0) || (!three_d && loc.lx3 != 0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart MeshBlock logical location is outside supported bounds."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (!std::isfinite(cost_eachmb[i]) || cost_eachmb[i] < 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart MeshBlock load-balance cost is invalid." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    total_cost += static_cast<double>(cost_eachmb[i]);
  }
  if (!std::isfinite(total_cost) || !(total_cost > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart MeshBlock load-balance total cost is invalid." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (!adaptive) max_level = current_level;
  if (adaptive && current_level > max_level) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart MeshBlock logical level exceeds configured AMR levels."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  restart_meta.rank_eachmb.assign(nmb_total, 0);
  if (restart_meta.original_nranks > 0) {
    restart_meta.gids_eachrank.assign(restart_meta.original_nranks, 0);
    restart_meta.nmb_eachrank.assign(restart_meta.original_nranks, 0);
  } else {
    restart_meta.gids_eachrank.clear();
    restart_meta.nmb_eachrank.clear();
  }

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(restart_meta.rank_eachmb.data(), sizeof(int), nmb_total,
                           single_file_per_rank)
        != static_cast<unsigned int>(nmb_total)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "MeshBlock rank list read from restart file is incorrect,"
                << " restart file is broken." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Bcast(restart_meta.rank_eachmb.data(), meshblock_int_mpi_bytes, MPI_CHAR, 0,
              MPI_COMM_WORLD);
  }
#endif

  if (restart_meta.original_nranks > 0) {
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(restart_meta.gids_eachrank.data(), sizeof(int),
                             restart_meta.original_nranks, single_file_per_rank)
          != static_cast<unsigned int>(restart_meta.original_nranks)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Rank gid table read from restart file is incorrect,"
                  << " restart file is broken." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(restart_meta.gids_eachrank.data(),
                rank_int_mpi_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif

    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(restart_meta.nmb_eachrank.data(), sizeof(int),
                             restart_meta.original_nranks, single_file_per_rank)
          != static_cast<unsigned int>(restart_meta.original_nranks)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Rank MeshBlock count read from restart file is incorrect,"
                  << " restart file is broken." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(restart_meta.nmb_eachrank.data(),
                rank_int_mpi_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
  }

  int partition_end = 0;
  for (int rank = 0; rank < restart_meta.original_nranks; ++rank) {
    const int start = restart_meta.gids_eachrank[rank];
    const int count = restart_meta.nmb_eachrank[rank];
    if (start != partition_end || count <= 0 || count > nmb_total - partition_end) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart rank layout is not a contiguous MeshBlock partition."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    partition_end += count;
  }
  if (partition_end != nmb_total) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart rank layout does not cover every MeshBlock." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  for (int gid = 0; gid < nmb_total; ++gid) {
    const int rank = restart_meta.rank_eachmb[gid];
    if (rank < 0 || rank >= restart_meta.original_nranks ||
        gid < restart_meta.gids_eachrank[rank] ||
        gid >= restart_meta.gids_eachrank[rank] + restart_meta.nmb_eachrank[rank]) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart MeshBlock rank assignment is inconsistent with rank layout."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }

  restart_meta.ncyc_since_ref.clear();
  restart_meta.checkpoint_nonce = 0;
  std::uint64_t mesh_metadata_magic = 0;
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    const IOWrapperSizeT extension_offset = resfile.GetPosition(single_file_per_rank);
    if (resfile.Read_bytes_at(&mesh_metadata_magic, 1, sizeof(mesh_metadata_magic),
                              extension_offset, single_file_per_rank)
        != sizeof(mesh_metadata_magic)) {
      mesh_metadata_magic = 0;
    }
    if (resfile.Seek(extension_offset, single_file_per_rank) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Unable to restore restart stream position after mesh-metadata peek."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Bcast(&mesh_metadata_magic, sizeof(mesh_metadata_magic), MPI_CHAR, 0,
              MPI_COMM_WORLD);
  }
#endif
  if (mesh_metadata_magic == restart_utils::kMeshMetadataMagic) {
    int mesh_metadata_version = 0;
    int has_refinement_cooldown = 0;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(&mesh_metadata_magic, 1, sizeof(mesh_metadata_magic),
                             single_file_per_rank) != sizeof(mesh_metadata_magic)
          || resfile.Read_bytes(&mesh_metadata_version, sizeof(int), 1,
                                single_file_per_rank) != 1
          || resfile.Read_bytes(&has_refinement_cooldown, sizeof(int), 1,
                                single_file_per_rank) != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Mesh restart metadata read is incomplete, restart file is broken."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&mesh_metadata_version, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&has_refinement_cooldown, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
#endif
    if ((mesh_metadata_version !=
         restart_utils::kMeshMetadataVersionWithoutCheckpointNonce &&
         mesh_metadata_version != restart_utils::kMeshMetadataVersion)
        || (has_refinement_cooldown != 0 && has_refinement_cooldown != 1)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Unsupported mesh restart metadata version="
                << mesh_metadata_version << " or refinement-cooldown flag="
                << has_refinement_cooldown << "." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (mesh_metadata_version >= restart_utils::kMeshMetadataVersion) {
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_bytes(&restart_meta.checkpoint_nonce,
                               sizeof(restart_meta.checkpoint_nonce), 1,
                               single_file_per_rank) != 1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Restart checkpoint nonce read is incomplete." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        MPI_Bcast(&restart_meta.checkpoint_nonce, sizeof(restart_meta.checkpoint_nonce),
                  MPI_CHAR, 0, MPI_COMM_WORLD);
      }
#endif
      if (restart_meta.checkpoint_nonce == 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Restart checkpoint nonce is invalid." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
    if (has_refinement_cooldown == 0) {
      // Keep the empty vector as the backward-compatible constructor signal.
    } else {
      restart_meta.ncyc_since_ref.assign(nmb_total, 0);
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_bytes(restart_meta.ncyc_since_ref.data(), sizeof(int),
                               nmb_total, single_file_per_rank)
            != static_cast<unsigned int>(nmb_total)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "MeshBlock refinement cooldown list read from restart file is "
                    << "incorrect, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        MPI_Bcast(restart_meta.ncyc_since_ref.data(), meshblock_int_mpi_bytes,
                  MPI_CHAR, 0, MPI_COMM_WORLD);
      }
#endif
      for (const int cooldown : restart_meta.ncyc_since_ref) {
        if (cooldown < 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "MeshBlock refinement cooldown list contains a negative value."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
    }
  }
  restart_meta.common_prefix_bytes =
      static_cast<std::uint64_t>(resfile.GetPosition(single_file_per_rank));

  // rebuild the MeshBlockTree
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();
  ptree->ResetRestartNodeBudget(nmb_total);
  for (int i=0; i<nmb_total; i++) {ptree->AddNodeWithoutRefinement(lloc_eachmb[i]);}

  // Check that the serialized leaves rebuild a complete physical tree in canonical order.
  {
    int nnb = 0;
    ptree->CountMeshBlocks(nnb);
    if (nnb != nmb_total || !ptree->IsRestartTreeComplete()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Tree reconstruction failed. Total number of blocks in "
        << "reconstructed tree=" << nnb << ", number in file=" << nmb_total << std::endl;
      restart_utils::AbortOnFatalError();
    }
    ptree->ResetMeshBlockCount(nnb);
    auto canonical_lloc_eachmb = std::make_unique<LogicalLocation[]>(nmb_total);
    ptree->CreateZOrderedLLList(canonical_lloc_eachmb.get(), nullptr, nnb);
    for (int gid=0; gid<nmb_total; ++gid) {
      if (!SameLogicalLocation(lloc_eachmb[gid], canonical_lloc_eachmb[gid])) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Restart MeshBlock logical locations are not in canonical gid order."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (nmb_total < global_variable::nranks) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
      << __LINE__ << std::endl
      << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
      << global_variable::nranks << ")" << std::endl;
    restart_utils::AbortOnFatalError();
  }
#endif

  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  pmb_pack->AddMeshBlocks(pin);
  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);

  // Fix maximum number of MeshBlocks per rank with AMR
  nmb_maxperrank = nmb_thisrank;
  if (adaptive) {
    if (pin->DoesParameterExist("mesh_refinement", "max_nmb_per_rank")) {
      nmb_maxperrank = pin->GetReal("mesh_refinement", "max_nmb_per_rank");
      if (nmb_maxperrank < nmb_thisrank) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "On rank=" << global_variable::my_rank << " Root grid requires "
          << "more MeshBlocks (nmb_thisrank=" << nmb_thisrank << ") than specified by "
          << "<mesh_refinement>/max_nmb_per_rank=" << nmb_maxperrank << std::endl;
        restart_utils::AbortOnFatalError();
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
  }

  // set remaining parameters, output diagnostics
  cfl_no = pin->GetReal("time", "cfl_number");
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::ValidateRestartShardCommonPrefix() const
//! \brief reject per-rank restart sets whose duplicated mesh metadata is inconsistent

void Mesh::ValidateRestartShardCommonPrefix() const {
  if (!restart_meta.single_file_per_rank || restart_meta.original_nranks <= 1) return;
  if (restart_meta.common_prefix_bytes == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart shard common-prefix size is invalid." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  constexpr std::size_t chunk_bytes = 64*1024;
  std::array<char, chunk_bytes> reference = {};
  std::array<char, chunk_bytes> candidate = {};
  const std::string base_prefix = restart_meta.base_dir.empty() ? "" :
      (restart_meta.base_dir == "/" ? "/" : restart_meta.base_dir + "/");
  const auto shard_path = [&](const int rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", rank);
    return base_prefix + rank_dir + "/" + restart_meta.file_name;
  };

  for (int rank=1; rank<restart_meta.original_nranks; ++rank) {
    std::ifstream rank_zero(shard_path(0), std::ios::binary);
    std::ifstream other(shard_path(rank), std::ios::binary);
    std::uint64_t remaining = restart_meta.common_prefix_bytes;
    while (remaining > 0) {
      const auto count = static_cast<std::streamsize>(
          std::min<std::uint64_t>(remaining, chunk_bytes));
      rank_zero.read(reference.data(), count);
      other.read(candidate.data(), count);
      if (rank_zero.gcount() != count || other.gcount() != count ||
          std::memcmp(reference.data(), candidate.data(),
                      static_cast<std::size_t>(count)) != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Restart shard common mesh metadata is inconsistent across files."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
      remaining -= static_cast<std::uint64_t>(count);
    }
  }
}
