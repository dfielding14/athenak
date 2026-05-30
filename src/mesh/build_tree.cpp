//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file build_tree.cpp
//! \brief Functions to build MeshBlock, both for new runs and restarts

#include <algorithm>
#include <cctype>
#include <cinttypes>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits> // numeric_limits<>
#include <memory> // make_unique<>
#include <set>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "restart_layout.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr int kMaxSafeLogicalLevel = 30;

[[noreturn]] void FailRestartLayout(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(EXIT_FAILURE);
}

bool ValidRestartAxis(int cells, int lower, int upper, int ghost_zones,
                      bool active) {
  if (cells <= 0) return false;
  if (!active) return cells == 1 && lower == 0 && upper == 0;
  return lower == ghost_zones &&
      static_cast<std::int64_t>(upper) ==
          static_cast<std::int64_t>(lower) + cells - 1;
}

int CheckedRootLevelForMeshBlocks(int nmb_root_max) {
  int root_level = 0;
  std::int64_t root_capacity = 1;
  while (root_capacity < nmb_root_max) {
    if (root_level == kMaxSafeLogicalLevel) {
      FailRestartLayout("MeshBlock grid exceeds the signed-safe logical-level range.");
    }
    ++root_level;
    root_capacity *= 2;
  }
  return root_level;
}

void ValidateRestartMeshHeader(const RegionSize &mesh_size,
                               const RegionIndcs &mesh_indcs,
                               const RegionIndcs &mb_indcs,
                               int root_level, bool multilevel,
                               bool multi_d, bool three_d) {
  bool finite_mesh_size = std::isfinite(mesh_size.x1min) &&
      std::isfinite(mesh_size.x1max) && std::isfinite(mesh_size.x2min) &&
      std::isfinite(mesh_size.x2max) && std::isfinite(mesh_size.x3min) &&
      std::isfinite(mesh_size.x3max) && std::isfinite(mesh_size.dx1) &&
      std::isfinite(mesh_size.dx2) && std::isfinite(mesh_size.dx3);
  if (!finite_mesh_size || mesh_size.x1max <= mesh_size.x1min ||
      mesh_size.x2max <= mesh_size.x2min || mesh_size.x3max <= mesh_size.x3min ||
      mesh_size.dx1 <= 0.0 || mesh_size.dx2 <= 0.0 || mesh_size.dx3 <= 0.0) {
    FailRestartLayout("restart mesh coordinate bounds are inconsistent.");
  }
  if (root_level < 0 || root_level > kMaxSafeLogicalLevel) {
    FailRestartLayout("restart root level must be between 0 and 30.");
  }
  if ((mesh_indcs.nx2 > 1) != multi_d || (mesh_indcs.nx3 > 1) != three_d ||
      (mb_indcs.nx2 > 1) != multi_d || (mb_indcs.nx3 > 1) != three_d) {
    FailRestartLayout("restart mesh dimensions disagree with input parameters.");
  }
  if (mesh_indcs.ng < 2 || mb_indcs.ng != mesh_indcs.ng ||
      (multilevel && mesh_indcs.ng % 2 != 0) ||
      !ValidRestartAxis(mesh_indcs.nx1, mesh_indcs.is, mesh_indcs.ie,
                        mesh_indcs.ng, true) ||
      !ValidRestartAxis(mesh_indcs.nx2, mesh_indcs.js, mesh_indcs.je,
                        mesh_indcs.ng, multi_d) ||
      !ValidRestartAxis(mesh_indcs.nx3, mesh_indcs.ks, mesh_indcs.ke,
                        mesh_indcs.ng, three_d) ||
      mesh_indcs.nx1 < 4 ||
      (multi_d && mesh_indcs.nx2 < 4) ||
      (three_d && mesh_indcs.nx3 < 4)) {
    FailRestartLayout("restart mesh indices are inconsistent.");
  }
  Real expected_dx1 = (mesh_size.x1max - mesh_size.x1min)/
      static_cast<Real>(mesh_indcs.nx1);
  Real expected_dx2 = (mesh_size.x2max - mesh_size.x2min)/
      static_cast<Real>(mesh_indcs.nx2);
  Real expected_dx3 = (mesh_size.x3max - mesh_size.x3min)/
      static_cast<Real>(mesh_indcs.nx3);
  if (mesh_size.dx1 != expected_dx1 || mesh_size.dx2 != expected_dx2 ||
      mesh_size.dx3 != expected_dx3) {
    FailRestartLayout("restart mesh spacing is inconsistent with its bounds.");
  }
  if (!ValidRestartAxis(mb_indcs.nx1, mb_indcs.is, mb_indcs.ie,
                        mb_indcs.ng, true) ||
      !ValidRestartAxis(mb_indcs.nx2, mb_indcs.js, mb_indcs.je,
                        mb_indcs.ng, multi_d) ||
      !ValidRestartAxis(mb_indcs.nx3, mb_indcs.ks, mb_indcs.ke,
                        mb_indcs.ng, three_d) ||
      mb_indcs.nx1 < 4 ||
      (multi_d && mb_indcs.nx2 < 4) ||
      (three_d && mb_indcs.nx3 < 4) ||
      (multilevel && (mb_indcs.nx1 % 2 != 0 ||
                      (multi_d && mb_indcs.nx2 % 2 != 0) ||
                      (three_d && mb_indcs.nx3 % 2 != 0))) ||
      mesh_indcs.nx1 % mb_indcs.nx1 != 0 ||
      mesh_indcs.nx2 % mb_indcs.nx2 != 0 ||
      mesh_indcs.nx3 % mb_indcs.nx3 != 0 ||
      mb_indcs.cnx1 != mb_indcs.nx1/2 ||
      mb_indcs.cnx2 != std::max(1, mb_indcs.nx2/2) ||
      mb_indcs.cnx3 != std::max(1, mb_indcs.nx3/2) ||
      mb_indcs.cis != mb_indcs.ng ||
      mb_indcs.cie != mb_indcs.cis + mb_indcs.cnx1 - 1 ||
      (multi_d && (mb_indcs.cjs != mb_indcs.ng ||
                   mb_indcs.cje != mb_indcs.cjs + mb_indcs.cnx2 - 1)) ||
      (!multi_d && (mb_indcs.cjs != 0 || mb_indcs.cje != 0)) ||
      (three_d && (mb_indcs.cks != mb_indcs.ng ||
                   mb_indcs.cke != mb_indcs.cks + mb_indcs.cnx3 - 1)) ||
      (!three_d && (mb_indcs.cks != 0 || mb_indcs.cke != 0))) {
    FailRestartLayout("restart MeshBlock indices are inconsistent.");
  }
  int nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  int nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  int nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;
  int nmb_root_max = std::max(nmb_rootx1, std::max(nmb_rootx2, nmb_rootx3));
  int expected_root_level = CheckedRootLevelForMeshBlocks(nmb_root_max);
  if (root_level != expected_root_level) {
    FailRestartLayout("restart root level is inconsistent with its MeshBlock grid.");
  }
}

bool ValidRestartLogicalAxis(std::int32_t location, int root_blocks, int level_delta,
                             bool active) {
  if (!active) return location == 0;
  if (level_delta < 0 || level_delta > kMaxSafeLogicalLevel ||
      root_blocks > (std::numeric_limits<std::int32_t>::max() >> level_delta)) {
    return false;
  }
  std::int64_t axis_blocks = static_cast<std::int64_t>(root_blocks) << level_delta;
  return location >= 0 && static_cast<std::int64_t>(location) < axis_blocks;
}

void ValidateRestartMeshBlockMetadata(const LogicalLocation &location, float cost,
                                      int root_level, int max_level,
                                      int nmb_rootx1, int nmb_rootx2,
                                      int nmb_rootx3, bool multilevel,
                                      bool multi_d, bool three_d) {
  if (location.level < root_level || location.level > max_level ||
      location.level > kMaxSafeLogicalLevel ||
      (!multilevel && location.level != root_level)) {
    FailRestartLayout("restart MeshBlock logical level is inconsistent.");
  }
  int level_delta = location.level - root_level;
  if (!ValidRestartLogicalAxis(location.lx1, nmb_rootx1, level_delta, true) ||
      !ValidRestartLogicalAxis(location.lx2, nmb_rootx2, level_delta, multi_d) ||
      !ValidRestartLogicalAxis(location.lx3, nmb_rootx3, level_delta, three_d)) {
    FailRestartLayout("restart MeshBlock logical location is inconsistent.");
  }
  if (!std::isfinite(cost) || cost <= 0.0) {
    FailRestartLayout("restart MeshBlock cost is inconsistent.");
  }
}

struct LogicalLocationLess {
  bool operator()(const LogicalLocation &left, const LogicalLocation &right) const {
    if (left.level != right.level) return left.level < right.level;
    if (left.lx1 != right.lx1) return left.lx1 < right.lx1;
    if (left.lx2 != right.lx2) return left.lx2 < right.lx2;
    return left.lx3 < right.lx3;
  }
};

LogicalLocation RestartParent(const LogicalLocation &location) {
  return {location.lx1 >> 1, location.lx2 >> 1, location.lx3 >> 1,
          location.level - 1};
}

LogicalLocation RestartChild(const LogicalLocation &location, int ox1, int ox2,
                             int ox3) {
  return {static_cast<std::int32_t>(location.lx1*2 + ox1),
          static_cast<std::int32_t>(location.lx2*2 + ox2),
          static_cast<std::int32_t>(location.lx3*2 + ox3),
          location.level + 1};
}

bool SameRestartLocation(const LogicalLocation &left, const LogicalLocation &right) {
  return left.lx1 == right.lx1 && left.lx2 == right.lx2 &&
      left.lx3 == right.lx3 && left.level == right.level;
}

bool RestartBranchIsComplete(
    const LogicalLocation &location,
    const std::set<LogicalLocation, LogicalLocationLess> &leaves,
    const std::set<LogicalLocation, LogicalLocationLess> &branches,
    bool multi_d, bool three_d) {
  if (leaves.count(location) != 0) return true;
  if (branches.count(location) == 0) return false;
  int nchild2 = multi_d ? 2 : 1;
  int nchild3 = three_d ? 2 : 1;
  for (int ox3 = 0; ox3 < nchild3; ++ox3) {
    for (int ox2 = 0; ox2 < nchild2; ++ox2) {
      for (int ox1 = 0; ox1 < 2; ++ox1) {
        if (!RestartBranchIsComplete(RestartChild(location, ox1, ox2, ox3),
                                     leaves, branches, multi_d, three_d)) {
          return false;
        }
      }
    }
  }
  return true;
}

void ValidateRestartMeshBlockInventory(
    const LogicalLocation *locations, const float *costs, int nmb_total,
    int root_level, int max_level, int nmb_rootx1, int nmb_rootx2,
    int nmb_rootx3, bool multilevel, bool multi_d, bool three_d) {
  std::set<LogicalLocation, LogicalLocationLess> leaves;
  std::set<LogicalLocation, LogicalLocationLess> branches;
  float total_cost = 0.0;
  for (int i = 0; i < nmb_total; ++i) {
    ValidateRestartMeshBlockMetadata(locations[i], costs[i], root_level, max_level,
                                     nmb_rootx1, nmb_rootx2, nmb_rootx3,
                                     multilevel, multi_d, three_d);
    if (!leaves.insert(locations[i]).second) {
      FailRestartLayout("restart MeshBlock inventory contains a duplicate location.");
    }
    if (branches.count(locations[i]) != 0) {
      FailRestartLayout("restart MeshBlock inventory contains overlapping levels.");
    }
    total_cost += costs[i];
    if (!std::isfinite(total_cost)) {
      FailRestartLayout("restart MeshBlock aggregate cost is inconsistent.");
    }
  }
  if (!std::isfinite(total_cost/static_cast<float>(global_variable::nranks)) ||
      total_cost/static_cast<float>(global_variable::nranks) <= 0.0) {
    FailRestartLayout("restart MeshBlock aggregate cost is inconsistent.");
  }
  for (int i = 0; i < nmb_total; ++i) {
    LogicalLocation ancestor = locations[i];
    while (ancestor.level > root_level) {
      ancestor = RestartParent(ancestor);
      if (leaves.count(ancestor) != 0) {
        FailRestartLayout("restart MeshBlock inventory contains overlapping levels.");
      }
      branches.insert(ancestor);
    }
  }
  for (int lx3 = 0; lx3 < nmb_rootx3; ++lx3) {
    for (int lx2 = 0; lx2 < nmb_rootx2; ++lx2) {
      for (int lx1 = 0; lx1 < nmb_rootx1; ++lx1) {
        LogicalLocation root_location{lx1, lx2, lx3, root_level};
        if (!RestartBranchIsComplete(root_location, leaves, branches, multi_d,
                                     three_d)) {
          FailRestartLayout("restart MeshBlock inventory has incomplete refinement.");
        }
      }
    }
  }
}

void ValidateRestartNeighborBalance(MeshBlockTree *tree,
                                    const LogicalLocation *locations,
                                    int nmb_total, bool multi_d, bool three_d) {
  int lower2 = multi_d ? -1 : 0;
  int upper2 = multi_d ? 1 : 0;
  int lower3 = three_d ? -1 : 0;
  int upper3 = three_d ? 1 : 0;
  for (int i = 0; i < nmb_total; ++i) {
    for (int ox3 = lower3; ox3 <= upper3; ++ox3) {
      for (int ox2 = lower2; ox2 <= upper2; ++ox2) {
        for (int ox1 = -1; ox1 <= 1; ++ox1) {
          if (ox1 != 0 || ox2 != 0 || ox3 != 0) {
            tree->FindNeighbor(locations[i], ox1, ox2, ox3);
          }
        }
      }
    }
  }
}

void ValidateRestartRankInventory(const int *gids_eachrank, const int *nmb_eachrank,
                                  int nranks, int nmb_total) {
  int next_gid = 0;
  for (int rank = 0; rank < nranks; ++rank) {
    if (gids_eachrank[rank] != next_gid || nmb_eachrank[rank] <= 0 ||
        nmb_eachrank[rank] > nmb_total - next_gid) {
      FailRestartLayout("restart load balancing produced an inconsistent rank map.");
    }
    next_gid += nmb_eachrank[rank];
  }
  if (next_gid != nmb_total) {
    FailRestartLayout("restart load balancing did not assign every MeshBlock.");
  }
}

int RestartMaxLevel(ParameterInput *pin, int root_level, bool adaptive) {
  if (!adaptive) return kMaxSafeLogicalLevel;
  std::string num_levels_text;
  if (pin->DoesParameterExist("mesh_refinement", "num_levels")) {
    num_levels_text = pin->GetString("mesh_refinement", "num_levels");
  } else {
    num_levels_text = std::to_string(
        pin->GetOrAddInteger("mesh_refinement", "num_levels", 1));
  }
  std::size_t parsed = 0;
  std::int64_t num_levels;
  try {
    num_levels = std::stoll(num_levels_text, &parsed, 10);
  } catch (const std::exception &) {
    FailRestartLayout("mesh_refinement/num_levels is not a valid integer.");
  }
  while (parsed < num_levels_text.size() &&
         std::isspace(static_cast<unsigned char>(num_levels_text[parsed]))) {
    ++parsed;
  }
  if (parsed != num_levels_text.size()) {
    FailRestartLayout("mesh_refinement/num_levels is not a valid integer.");
  }
  std::int64_t max_num_levels = kMaxSafeLogicalLevel - root_level + 1;
  if (num_levels <= 0 || num_levels > max_num_levels) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Number of refinement levels must be between 1 and "
              << max_num_levels << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::int64_t max_level = num_levels + root_level - 1;
  return static_cast<int>(max_level);
}

}  // namespace

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

  // Find smallest N such that 2^N >= max number of MeshBlocks in any dimension (nmbmax)
  // Then N is logical level of root grid.
  root_level = CheckedRootLevelForMeshBlocks(nmbmax);
  int current_level = root_level;

  // Construct tree and create root grid
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();

  // Error check properties of input paraemters for SMR/AMR meshes.
  max_level = RestartMaxLevel(pin, root_level, adaptive);

  // Read <refined_region> blocks and construct tree accordingly
  // These regions can be used with both SMR (in which case they will remain fixed) and
  // AMR (in which case they may be defined, unless the location refinement criteria used)
  if (multilevel) {
    // error check that number of cells in MeshBlock divisible by two
    if (mb_indcs.nx1 % 2 != 0 ||
       (mb_indcs.nx2 % 2 != 0 && multi_d) ||
       (mb_indcs.nx3 % 2 != 0 && three_d)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of cells in MeshBlock must be divisible by 2 "
                << "with SMR or AMR." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // cycle through ParameterInput list and find "refined_region" blocks, extract data
    // and expand MeshBlockTree
    for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
      if (it->block_name.compare(0, 14, "refined_region") == 0) {
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
        int log_ref_lev = phy_ref_lev + root_level;
        if (log_ref_lev > current_level) current_level = log_ref_lev;

        // error check parameters in "refinement" blocks
        if (phy_ref_lev < 1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl <<"<refined_region> level must be larger than 0 (root level=0)"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (log_ref_lev > max_level) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> level exceeds maximum allowed ("
              << max_level << ")" << std::endl << "Reduce/specify 'num_levels' in "
              << "<mesh_refinement> input block if using AMR" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min > ref_size.x1max
            || ref_size.x2min > ref_size.x2max
            || ref_size.x3min > ref_size.x3max)  {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid <refined_region> (xmax < xmin in one direction)."
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
            || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
            || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> must be fully contained within root mesh"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }

        // note: if following is too slow, it could be replaced with bi-section search.
        // Suppose entire root domain is tiled with MeshBlocks at the desired refinement
        // level. Find range of x1-integer indices of such MeshBlocks that cover the
        // refinement region
        std::int32_t lx1min = 0, lx1max = 0;
        std::int32_t lx2min = 0, lx2max = 0;
        std::int32_t lx3min = 0, lx3max = 0;
        std::int32_t lxmax = nmb_rootx1*(1<<phy_ref_lev);
        for (lx1min=0; lx1min<lxmax; lx1min++) {
          if (LeftEdgeX(lx1min+1,lxmax,mesh_size.x1min,mesh_size.x1max) > ref_size.x1min)
            break;
        }
        for (lx1max=lx1min; lx1max<lxmax; lx1max++) {
          if (LeftEdgeX(lx1max+1,lxmax,mesh_size.x1min,mesh_size.x1max) >= ref_size.x1max)
            break;
        }
        if (lx1min % 2 == 1) lx1min--;
        if (lx1max % 2 == 0) lx1max++;

        // Find range of x2-indices of such MeshBlocks that cover the refinement region
        if (multi_d) { // 2D or 3D
          lxmax = nmb_rootx2*(1<<phy_ref_lev);
          for (lx2min=0; lx2min<lxmax; lx2min++) {
            if (LeftEdgeX(lx2min+1, lxmax, mesh_size.x2min, mesh_size.x2max) >
                ref_size.x2min)
            break;
          }
          for (lx2max=lx2min; lx2max<lxmax; lx2max++) {
            if (LeftEdgeX(lx2max+1, lxmax, mesh_size.x2min, mesh_size.x2max) >=
                ref_size.x2max)
            break;
          }
          if (lx2min % 2 == 1) lx2min--;
          if (lx2max % 2 == 0) lx2max++;
        }

        // Find range of x3-indices of such MeshBlocks that cover the refinement region
        if (three_d) { // 3D
          lxmax = nmb_rootx3*(1<<phy_ref_lev);
          for (lx3min=0; lx3min<lxmax; lx3min++) {
            if (LeftEdgeX(lx3min+1, lxmax, mesh_size.x3min, mesh_size.x3max) >
                ref_size.x3min)
            break;
          }
          for (lx3max=lx3min; lx3max<lxmax; lx3max++) {
            if (LeftEdgeX(lx3max+1, lxmax, mesh_size.x3min, mesh_size.x3max) >=
                ref_size.x3max)
            break;
          }
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
            int nnew;
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
              int nnew;
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
                int nnew;
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
    std::exit(EXIT_FAILURE);
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
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  if (nmb_maxperrank > (1 << (NUM_BITS_LID))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "Maximum number of MeshBlocks per rank cannot exceed 2^(NUM_BITS_LID) due to MPI"
      << " tag limits" << std::endl;
    std::exit(EXIT_FAILURE);
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
  // At this point, the restartfile is already open and the ParameterInput (input file)
  // data has already been read in main(). Thus the file pointer is set to after <par_end>
  IOWrapperSizeT headeroffset = resfile.GetPosition(single_file_per_rank);

  // following must be identical to calculation of headeroffset (excluding size of
  // ParameterInput data) in restart.cpp
  IOWrapperSizeT headersize = 3*sizeof(int) + 2*sizeof(Real)
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
      exit(EXIT_FAILURE);
    }
  }

#if MPI_PARALLEL_ENABLED
  // then broadcast the header data
  if (!single_file_per_rank) {
    io_wrapper::BroadcastBytes(headerdata, headersize, 0, MPI_COMM_WORLD);
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
  delete [] headerdata;
  ValidateRestartMeshHeader(mesh_size, mesh_indcs, mb_indcs, root_level, multilevel,
                            multi_d, three_d);
  restart_layout::Size nmb_total_size =
      restart_layout::CheckedPositive(nmb_total, FailRestartLayout,
                                      "restart total MeshBlock count");

  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;
  int current_level = root_level;

  // Error check properties of input paraemters for SMR/AMR meshes.
  max_level = RestartMaxLevel(pin, root_level, adaptive);

  // allocate memory for lists read from restart
  restart_layout::CheckedMemorySize(
      restart_layout::CheckedMultiply(nmb_total_size, sizeof(float), FailRestartLayout,
                                      "restart MeshBlock-cost allocation bytes"),
      FailRestartLayout, "restart MeshBlock-cost allocation bytes");
  restart_layout::CheckedMemorySize(
      restart_layout::CheckedMultiply(nmb_total_size, sizeof(int), FailRestartLayout,
                                      "restart MeshBlock-rank allocation bytes"),
      FailRestartLayout, "restart MeshBlock-rank allocation bytes");
  restart_layout::CheckedMemorySize(
      restart_layout::CheckedMultiply(nmb_total_size, sizeof(LogicalLocation),
                                      FailRestartLayout,
                                      "restart logical-location allocation bytes"),
      FailRestartLayout, "restart logical-location allocation bytes");
  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // allocate idlist buffer and read list of logical locations and cost
  IOWrapperSizeT listsize = restart_layout::CheckedAdd(
      sizeof(LogicalLocation), sizeof(float), FailRestartLayout,
      "restart MeshBlock metadata record bytes");
  IOWrapperSizeT metadata_bytes = restart_layout::CheckedMultiply(
      listsize, nmb_total_size, FailRestartLayout,
      "restart MeshBlock metadata bytes");
  char *idlist = new char[restart_layout::CheckedMemorySize(
      metadata_bytes, FailRestartLayout, "restart MeshBlock metadata bytes")];
  // only the master process reads the ID list
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(idlist, 1, metadata_bytes, single_file_per_rank) !=
        metadata_bytes) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Incorrect MeshBlock metadata size in restart file; "
                << "restart file is broken." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the ID list
  if (!single_file_per_rank) {
    io_wrapper::BroadcastBytes(idlist, metadata_bytes, 0, MPI_COMM_WORLD);
  }
#endif

  // everyone sets the logical location and cost lists based on broadcast data
  IOWrapperSizeT os = 0;
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(lloc_eachmb[i]), &(idlist[restart_layout::CheckedMemorySize(
        os, FailRestartLayout, "restart MeshBlock metadata offset")]),
        sizeof(LogicalLocation));
    os = restart_layout::CheckedAdd(os, sizeof(LogicalLocation), FailRestartLayout,
                                    "restart MeshBlock metadata offset");
  }
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(cost_eachmb[i]), &(idlist[restart_layout::CheckedMemorySize(
        os, FailRestartLayout, "restart MeshBlock metadata offset")]), sizeof(float));
    os = restart_layout::CheckedAdd(os, sizeof(float), FailRestartLayout,
                                    "restart MeshBlock metadata offset");
    if (lloc_eachmb[i].level > current_level) current_level = lloc_eachmb[i].level;
  }
  delete [] idlist;
  ValidateRestartMeshBlockInventory(lloc_eachmb, cost_eachmb, nmb_total, root_level,
                                    max_level, nmb_rootx1, nmb_rootx2, nmb_rootx3,
                                    multilevel, multi_d, three_d);
  if (!adaptive) max_level = current_level;

  // rebuild the MeshBlockTree
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();
  for (int i=0; i<nmb_total; i++) {ptree->AddNodeWithoutRefinement(lloc_eachmb[i]);}
  ValidateRestartNeighborBalance(ptree.get(), lloc_eachmb, nmb_total, multi_d, three_d);

  // check the tree structure by making sure total # of MBs counted in tree same as the
  // number read from the restart file. Persisted records must remain in canonical order.
  {
    int nnb;
    LogicalLocation *zordered_lloc = new LogicalLocation[nmb_total];
    ptree->CreateZOrderedLLList(zordered_lloc, nullptr, nnb);
    if (nnb != nmb_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Tree reconstruction failed. Total number of blocks in "
        << "reconstructed tree=" << nnb << ", number in file=" << nmb_total << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nmb_total; ++i) {
      if (!SameRestartLocation(lloc_eachmb[i], zordered_lloc[i])) {
        FailRestartLayout("restart MeshBlock inventory is not in canonical order.");
      }
    }
    delete [] zordered_lloc;
  }

#ifdef MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (!single_file_per_rank) {
    if (nmb_total < global_variable::nranks) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
        << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#endif

  std::fill(gids_eachrank, gids_eachrank + global_variable::nranks, -1);
  std::fill(nmb_eachrank, nmb_eachrank + global_variable::nranks, 0);
  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);
  ValidateRestartRankInventory(gids_eachrank, nmb_eachrank, global_variable::nranks,
                               nmb_total);

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
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "With AMR maximum number of MeshBlocks per rank must be "
        << "specified in input file using <mesh_refinement>/max_nmb_per_rank"
        << std::endl;
      std::exit(EXIT_FAILURE);
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
