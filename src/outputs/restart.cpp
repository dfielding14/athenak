//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir
#include <sys/time.h>  // gettimeofday
#include <unistd.h>    // getpid

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <initializer_list>
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
#include "particles/particles.hpp"
#include "outputs/restart_utils.hpp"
//#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr IOWrapperSizeT kMaxSupportedRestartParticlePayloadBytesPerRank =
    static_cast<IOWrapperSizeT>(8) << 30;

void AbortRestartLayoutOverflow() {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl
            << "Restart output layout exceeds supported offset or allocation limits."
            << std::endl;
  restart_utils::AbortOnFatalError();
}

IOWrapperSizeT CheckedRestartProduct(IOWrapperSizeT count, IOWrapperSizeT width) {
  constexpr IOWrapperSizeT max = std::numeric_limits<IOWrapperSizeT>::max();
  if (count != 0 && width > max/count) AbortRestartLayoutOverflow();
  return count*width;
}

IOWrapperSizeT CheckedRestartProduct(std::initializer_list<IOWrapperSizeT> factors) {
  IOWrapperSizeT product = 1;
  for (const IOWrapperSizeT factor : factors) {
    product = CheckedRestartProduct(product, factor);
  }
  return product;
}

void AdvanceRestartOffset(IOWrapperSizeT &offset, IOWrapperSizeT count,
                          IOWrapperSizeT width) {
  constexpr IOWrapperSizeT max = std::numeric_limits<IOWrapperSizeT>::max();
  const IOWrapperSizeT increment = CheckedRestartProduct(count, width);
  if (offset > max - increment) AbortRestartLayoutOverflow();
  offset += increment;
}

void AdvanceRestartRealArrayOffset(IOWrapperSizeT &offset,
                                   std::initializer_list<int> extents) {
  IOWrapperSizeT count = 1;
  for (const int extent : extents) {
    if (extent < 0) AbortRestartLayoutOverflow();
    count = CheckedRestartProduct(count, static_cast<IOWrapperSizeT>(extent));
  }
  AdvanceRestartOffset(offset, count, sizeof(Real));
}

int CheckedRestartOutputExtent(int active_cells, int ghost_cells) {
  if (active_cells <= 0 || ghost_cells < 0) AbortRestartLayoutOverflow();
  IOWrapperSizeT extent = static_cast<IOWrapperSizeT>(active_cells);
  if (active_cells > 1) {
    AdvanceRestartOffset(extent, static_cast<IOWrapperSizeT>(ghost_cells), 2);
  }
  if (extent > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    AbortRestartLayoutOverflow();
  }
  return static_cast<int>(extent);
}

int CheckedRestartExtentPlusOne(int extent) {
  if (extent < 0 || extent == std::numeric_limits<int>::max()) {
    AbortRestartLayoutOverflow();
  }
  return extent + 1;
}

int CheckedRestartIntProduct(std::initializer_list<int> factors) {
  IOWrapperSizeT product = 1;
  for (const int factor : factors) {
    if (factor < 0) AbortRestartLayoutOverflow();
    product = CheckedRestartProduct(product, static_cast<IOWrapperSizeT>(factor));
  }
  if (product > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    AbortRestartLayoutOverflow();
  }
  return static_cast<int>(product);
}

int CheckedRestartIntCount(std::size_t count) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    AbortRestartLayoutOverflow();
  }
  return static_cast<int>(count);
}

int CheckedRestartIntAdd(int lhs, int rhs) {
  if ((rhs > 0 && lhs > std::numeric_limits<int>::max() - rhs) ||
      (rhs < 0 && lhs < std::numeric_limits<int>::min() - rhs)) {
    AbortRestartLayoutOverflow();
  }
  return lhs + rhs;
}

std::size_t CheckedRestartSizeT(IOWrapperSizeT count) {
  if (count > static_cast<IOWrapperSizeT>(std::numeric_limits<std::size_t>::max())) {
    AbortRestartLayoutOverflow();
  }
  return static_cast<std::size_t>(count);
}

std::size_t CheckedRestartAllocationSize(IOWrapperSizeT count, std::size_t width) {
  constexpr std::size_t max = std::numeric_limits<std::size_t>::max();
  if (width == 0 || count > static_cast<IOWrapperSizeT>(max/width)) {
    AbortRestartLayoutOverflow();
  }
  return static_cast<std::size_t>(count);
}

void AdvanceRestartPastFaceArrayRemainder(
    IOWrapperSizeT &offset, IOWrapperSizeT data_size,
    std::initializer_list<IOWrapperSizeT> face_counts) {
  IOWrapperSizeT face_bytes = 0;
  for (const IOWrapperSizeT count : face_counts) {
    AdvanceRestartOffset(face_bytes, count, sizeof(Real));
  }
  if (face_bytes > data_size) AbortRestartLayoutOverflow();
  AdvanceRestartOffset(offset, data_size - face_bytes, 1);
}

void ValidateRestartAllocation(std::initializer_list<int> extents) {
  std::size_t count = 1;
  constexpr std::size_t max = std::numeric_limits<std::size_t>::max();
  for (const int extent : extents) {
    if (extent < 0 || (extent != 0 &&
                       count > max/static_cast<std::size_t>(extent))) {
      AbortRestartLayoutOverflow();
    }
    count *= static_cast<std::size_t>(extent);
  }
  if (count > max/sizeof(Real)) AbortRestartLayoutOverflow();
}

std::string CheckedRestartSequenceSuffix(const int file_number) {
  if (file_number < 0 || file_number == std::numeric_limits<int>::max()) {
    AbortRestartLayoutOverflow();
  }
  std::ostringstream suffix;
  suffix << '.' << std::setfill('0') << std::setw(5) << file_number;
  return suffix.str();
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  if (mkdir("rst",0775) != 0 && errno != EEXIST) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unable to create restart output directory 'rst': "
              << std::strerror(errno) << std::endl;
    restart_utils::AbortOnFatalError();
  }
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rst/rank_%08d/", global_variable::my_rank);
    if (mkdir(rank_dir, 0775) != 0 && errno != EEXIST) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Unable to create restart rank directory '" << rank_dir
                << "': " << std::strerror(errno) << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedRestartOutputExtent(indcs.nx1, indcs.ng);
  int nout2 = CheckedRestartOutputExtent(indcs.nx2, indcs.ng);
  int nout3 = CheckedRestartOutputExtent(indcs.nx3, indcs.ng);
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
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
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

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  if (phydro != nullptr) {
    ValidateRestartAllocation({nmb, nhydro, nout3, nout2, nout1});
    Kokkos::realloc(outarray_hyd, nmb, nhydro, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    ValidateRestartAllocation({nmb, nmhd, nout3, nout2, nout1});
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    const int nout1p1 = CheckedRestartExtentPlusOne(nout1);
    const int nout2p1 = CheckedRestartExtentPlusOne(nout2);
    const int nout3p1 = CheckedRestartExtentPlusOne(nout3);
    ValidateRestartAllocation({nmb, nout3, nout2, nout1p1});
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1p1);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    ValidateRestartAllocation({nmb, nout3, nout2p1, nout1});
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2p1, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    ValidateRestartAllocation({nmb, nout3p1, nout2, nout1});
    Kokkos::realloc(outfield.x3f, nmb, nout3p1, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (prad != nullptr) {
    ValidateRestartAllocation({nmb, nrad, nout3, nout2, nout1});
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pturb != nullptr) {
    ValidateRestartAllocation({nmb, nforce, nout3, nout2, nout1});
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force_tmp1,
                      std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    ValidateRestartAllocation({nmb, nz4c, nout3, nout2, nout1});
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    ValidateRestartAllocation({nmb, nadm, nout3, nout2, nout1});
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
  if (pm->nmb_total <= 0 || pm->nmb_total > kMaxSupportedMeshBlocks) {
    AbortRestartLayoutOverflow();
  }
  if (pm->adaptive) {
    if (pm->pmr == nullptr) {
      AbortRestartLayoutOverflow();
    }
    if (pm->pmr->ncyc_since_ref.extent(0) <
        static_cast<std::size_t>(pm->nmb_total)) {
      AbortRestartLayoutOverflow();
    }
    for (int m=0; m<pm->nmb_total; ++m) {
      if (pm->pmr->ncyc_since_ref(m) < 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Adaptive MeshBlock cooldown state is negative." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
  }

  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedRestartOutputExtent(indcs.nx1, indcs.ng);
  int nout2 = CheckedRestartOutputExtent(indcs.nx2, indcs.ng);
  int nout3 = CheckedRestartOutputExtent(indcs.nx3, indcs.ng);
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nz4c=0, nadm=0, nco=0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
    nco = CheckedRestartIntCount(pz4c->ptracker.size());
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  bool single_file_per_rank = out_params.single_file_per_rank;
  std::string fname;
  const std::string number = CheckedRestartSequenceSuffix(out_params.file_number);
  if (single_file_per_rank) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = file_number padded to at least 5 digits
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("rst/") + std::string(rank_dir) + out_params.file_basename
      + number + ".rst";

    // Debugging output to check directory and filename
    // std::cout << "Rank " << global_variable::my_rank << " generated filename: "
    //           << fname << std::endl;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = file_number padded to at least 5 digits
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  const std::string partial_fname = fname + ".partial";
  const std::string manifest_fname =
      std::string("rst/") + out_params.file_basename + number + ".rst.manifest";
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  std::uint64_t checkpoint_nonce = 0;
  if (global_variable::my_rank == 0) {
    timeval wall_time = {};
    if (gettimeofday(&wall_time, nullptr) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Unable to create restart checkpoint nonce." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    checkpoint_nonce = static_cast<std::uint64_t>(wall_time.tv_sec)*1000000;
    checkpoint_nonce += static_cast<std::uint64_t>(wall_time.tv_usec);
    checkpoint_nonce ^= static_cast<std::uint64_t>(getpid()) << 16;
    checkpoint_nonce ^= static_cast<std::uint64_t>(pm->ncycle) << 32;
    checkpoint_nonce ^= static_cast<std::uint64_t>(out_params.file_number);
    if (checkpoint_nonce == 0) checkpoint_nonce = 1;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&checkpoint_nonce, sizeof(checkpoint_nonce), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
  resfile.Open(partial_fname.c_str(), IOWrapper::FileMode::write,
               single_file_per_rank);
  auto write_header_bytes = [&](const void *buffer, IOWrapperSizeT count) {
    if (resfile.Write_any_type(buffer, count, "byte", single_file_per_rank) != count) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to write restart header to partial artifact '"
                << partial_fname << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  };
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // output the input parameters (input file)
    write_header_bytes(sbuf.c_str(), sbuf.size());

    // output Mesh information
    write_header_bytes(&(pm->nmb_total), sizeof(int));
    write_header_bytes(&(pm->root_level), sizeof(int));
    write_header_bytes(&(pm->mesh_size), sizeof(RegionSize));
    write_header_bytes(&(pm->mesh_indcs), sizeof(RegionIndcs));
    write_header_bytes(&(pm->mb_indcs), sizeof(RegionIndcs));
    write_header_bytes(&(pm->time), sizeof(Real));
    write_header_bytes(&(pm->dt), sizeof(Real));
    write_header_bytes(&(pm->ncycle), sizeof(int));
    write_header_bytes(&(global_variable::nranks), sizeof(int));
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    write_header_bytes(&(pm->lloc_eachmb[0]),
                       CheckedRestartProduct(pm->nmb_total, sizeof(LogicalLocation)));
    write_header_bytes(&(pm->cost_eachmb[0]),
                       CheckedRestartProduct(pm->nmb_total, sizeof(float)));
    write_header_bytes(&(pm->rank_eachmb[0]),
                       CheckedRestartProduct(pm->nmb_total, sizeof(int)));
    if (global_variable::nranks > 0) {
      write_header_bytes(&(pm->gids_eachrank[0]),
                         CheckedRestartProduct(global_variable::nranks, sizeof(int)));
      write_header_bytes(&(pm->nmb_eachrank[0]),
                         CheckedRestartProduct(global_variable::nranks, sizeof(int)));
    }
    const std::uint64_t mesh_metadata_magic = restart_utils::kMeshMetadataMagic;
    const int mesh_metadata_version = restart_utils::kMeshMetadataVersion;
    const int has_refinement_cooldown = pm->adaptive ? 1 : 0;
    write_header_bytes(&mesh_metadata_magic, sizeof(mesh_metadata_magic));
    write_header_bytes(&mesh_metadata_version, sizeof(mesh_metadata_version));
    write_header_bytes(&has_refinement_cooldown, sizeof(has_refinement_cooldown));
    write_header_bytes(&checkpoint_nonce, sizeof(checkpoint_nonce));
    if (pm->adaptive) {
      write_header_bytes(pm->pmr->ncyc_since_ref.data(),
                         CheckedRestartProduct(pm->nmb_total, sizeof(int)));
    }
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // store z4c information
    if (pz4c != nullptr) {
      write_header_bytes(&(pz4c->last_output_time), sizeof(Real));
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        write_header_bytes(pt.GetPos(), 3*sizeof(Real));
      }
    }
    // Turbulence-driver OU accumulator metadata and internal RNG.  The raw
    // accumulator itself is stored in outarray_force below.
    if (pturb != nullptr) {
      TurbulenceRestartState state{1, pturb->n_turb_updates_yet};
      write_header_bytes(&state, sizeof(TurbulenceRestartState));
      write_header_bytes(&(pturb->rstate), sizeof(RNG_State));
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  IOWrapperSizeT data_size = 0;
  if (phydro != nullptr) {
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nhydro});
  }
  if (pmhd != nullptr) {
    const int nout1p1 = CheckedRestartExtentPlusOne(nout1);
    const int nout2p1 = CheckedRestartExtentPlusOne(nout2);
    const int nout3p1 = CheckedRestartExtentPlusOne(nout3);
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nmhd});
    AdvanceRestartRealArrayOffset(data_size, {nout1p1, nout2, nout3});
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2p1, nout3});
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3p1});
  }
  if (prad != nullptr) {
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nrad});
  }
  if (pturb != nullptr) {
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nforce});
  }
  if (pz4c != nullptr) {
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nz4c});
  } else if (padm != nullptr) {
    AdvanceRestartRealArrayOffset(data_size, {nout1, nout2, nout3, nadm});
  }
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    write_header_bytes(&(data_size), sizeof(IOWrapperSizeT));
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT step1size = 0;
  AdvanceRestartOffset(step1size, sbuf.size(), sizeof(char));
  AdvanceRestartOffset(step1size, 4, sizeof(int));
  AdvanceRestartOffset(step1size, 2, sizeof(Real));
  AdvanceRestartOffset(step1size, 1, sizeof(RegionSize));
  AdvanceRestartOffset(step1size, 2, sizeof(RegionIndcs));
  IOWrapperSizeT step2size = 0;
  AdvanceRestartOffset(step2size, pm->nmb_total, sizeof(LogicalLocation));
  AdvanceRestartOffset(step2size, pm->nmb_total, sizeof(float));
  AdvanceRestartOffset(step2size, pm->nmb_total, sizeof(int));
  AdvanceRestartOffset(step2size, global_variable::nranks, 2*sizeof(int));
  AdvanceRestartOffset(step2size, 1, sizeof(std::uint64_t));
  AdvanceRestartOffset(step2size, 2, sizeof(int));
  AdvanceRestartOffset(step2size, 1, sizeof(std::uint64_t));
  if (pm->adaptive) AdvanceRestartOffset(step2size, pm->nmb_total, sizeof(int));

  IOWrapperSizeT step3size = 0;
  AdvanceRestartOffset(step3size, nco, 3*sizeof(Real));
  if (pz4c != nullptr) AdvanceRestartOffset(step3size, 1, sizeof(Real));
  if (pturb != nullptr) {
    AdvanceRestartOffset(step3size, 1, sizeof(TurbulenceRestartState));
    AdvanceRestartOffset(step3size, 1, sizeof(RNG_State));
  }

  // write cell-centered variables in parallel
  IOWrapperSizeT offset_myrank = step1size;
  AdvanceRestartOffset(offset_myrank, step2size, 1);
  AdvanceRestartOffset(offset_myrank, step3size, 1);
  AdvanceRestartOffset(offset_myrank, 1, sizeof(IOWrapperSizeT));

  if (!single_file_per_rank) {
    AdvanceRestartOffset(offset_myrank,
                         pm->gids_eachrank[global_variable::my_rank], data_size);
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
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nhydro});
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nmhd});
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        auto fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        AdvanceRestartPastFaceArrayRemainder(
            myoffset, data_size, {x1fptr.size(), x2fptr.size(), x3fptr.size()});

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        auto fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, fldcnt, sizeof(Real));

        AdvanceRestartPastFaceArrayRemainder(
            myoffset, data_size, {x1fptr.size(), x2fptr.size(), x3fptr.size()});
      }
    }
    AdvanceRestartRealArrayOffset(
        offset_myrank, {CheckedRestartExtentPlusOne(nout1), nout2, nout3});
    AdvanceRestartRealArrayOffset(
        offset_myrank, {nout1, CheckedRestartExtentPlusOne(nout2), nout3});
    AdvanceRestartRealArrayOffset(
        offset_myrank, {nout1, nout2, CheckedRestartExtentPlusOne(nout3)});
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nrad});
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nforce});
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nz4c});
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceRestartOffset(myoffset, data_size, 1);
      }
    }
    AdvanceRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nadm});
    myoffset = offset_myrank;
  }

  //--- STEP 5. Optional particle restart state.
  particles::Particles *ppart = pm->pmb_pack->ppart;
  if (ppart != nullptr) {
    constexpr std::uint64_t kPicMagic = 0x5049435253543031ULL;
    constexpr int kPicVersion = particles::Particles::PIC_RESTART_SCHEMA_VERSION;

    const int nmb_local = pm->nmb_thisrank;
    const int gids_local = pm->gids_eachrank[global_variable::my_rank];
    const int npart = ppart->nprtcl_thispack;
    const int nrdata = ppart->nrdata;
    const int nidata = ppart->nidata;
    if (nmb_local < 0 || npart < 0 || nrdata <= 0 || nidata <= 0) {
      AbortRestartLayoutOverflow();
    }
    const IOWrapperSizeT packed_pr_count = CheckedRestartProduct(npart, nrdata);
    const IOWrapperSizeT packed_pi_count = CheckedRestartProduct(npart, nidata);
    const IOWrapperSizeT packed_pr_bytes =
        CheckedRestartProduct(packed_pr_count, sizeof(Real));
    const IOWrapperSizeT packed_pi_bytes =
        CheckedRestartProduct(packed_pi_count, sizeof(int));
    if (packed_pi_bytes > kMaxSupportedRestartParticlePayloadBytesPerRank ||
        packed_pr_bytes >
            kMaxSupportedRestartParticlePayloadBytesPerRank - packed_pi_bytes) {
      AbortRestartLayoutOverflow();
    }
    const int moment_cnt =
        CheckedRestartIntProduct({particles::Particles::NMOM, nout3, nout2, nout1});
    const bool has_moments = ppart->deposit_moments && (ppart->moments.size() > 0);
    const bool requires_moments = ppart->couple_moments_to_mhd && ppart->deposit_moments;
    const bool has_edge = (ppart->couple_j_to_efield_representation ==
                           CoupledCurrentRepresentation::edge_staggered) &&
                          (ppart->j_edge_x1e.size() > 0);
    const int nout1p1 = CheckedRestartExtentPlusOne(nout1);
    const int nout2p1 = CheckedRestartExtentPlusOne(nout2);
    const int nout3p1 = CheckedRestartExtentPlusOne(nout3);
    const int edge1_cnt = CheckedRestartIntProduct({nout3p1, nout2p1, nout1});
    const int edge2_cnt = CheckedRestartIntProduct({nout3p1, nout2, nout1p1});
    const int edge3_cnt = CheckedRestartIntProduct({nout3, nout2p1, nout1p1});
    const int state_kind = ppart->UsesRelativisticCRState() ? 1 : 0;
    const int physical_mode = static_cast<int>(ppart->pic_physical_mode);
    const Real cr_light_speed = ppart->pic_cr_light_speed;
    std::array<int, particles::Particles::NPIC_RESTART_MODEL_INTS> model_ints;
    std::array<Real, particles::Particles::NPIC_RESTART_MODEL_REALS> model_reals;
    ppart->FillRestartModelMetadata(model_ints, model_reals);

    if (requires_moments && !has_moments) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart requested coupled moments, but moments are not "
                << "allocated."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }

    auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
    auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_idata);

    std::vector<int> local_mb_counts(nmb_local, 0);
    for (int p=0; p<npart; ++p) {
      const int m = h_pi(PGID, p) - gids_local;
      if ((m < 0) || (m >= nmb_local)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Particle gid is not local at restart write time (gid="
                  << h_pi(PGID, p) << ", local gid range=[" << gids_local << ","
                  << (gids_local + nmb_local - 1) << "])." << std::endl;
        restart_utils::AbortOnFatalError();
      }
      if (ppart->particle_type == ParticleType::cosmic_ray) {
        const int sp = h_pi(PSP, p);
        if (sp < 0 || sp >= ppart->nspecies) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Particle species is out of range at restart write time "
                    << "(species=" << sp << ", nspecies=" << ppart->nspecies << ")."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
      if (local_mb_counts[m] == std::numeric_limits<int>::max()) {
        AbortRestartLayoutOverflow();
      }
      local_mb_counts[m] += 1;
    }

    std::vector<IOWrapperSizeT> local_mb_offsets(
        CheckedRestartExtentPlusOne(nmb_local), 0);
    for (int m=0; m<nmb_local; ++m) {
      local_mb_offsets[m + 1] = local_mb_offsets[m];
      AdvanceRestartOffset(local_mb_offsets[m + 1], local_mb_counts[m], 1);
    }

    std::vector<Real> packed_pr(
        CheckedRestartAllocationSize(packed_pr_count, sizeof(Real)), 0.0);
    std::vector<int> packed_pi(
        CheckedRestartAllocationSize(packed_pi_count, sizeof(int)), 0);
    std::vector<int> local_mb_cursor(nmb_local, 0);
    for (int p=0; p<npart; ++p) {
      const int m = h_pi(PGID, p) - gids_local;
      const int pinmb = local_mb_cursor[m]++;
      IOWrapperSizeT lp = local_mb_offsets[m];
      AdvanceRestartOffset(lp, pinmb, 1);
      IOWrapperSizeT real_index = CheckedRestartProduct(lp, nrdata);
      IOWrapperSizeT int_index = CheckedRestartProduct(lp, nidata);
      for (int n=0; n<nrdata; ++n) {
        packed_pr[CheckedRestartSizeT(real_index)] = h_pr(n, p);
        AdvanceRestartOffset(real_index, 1, 1);
      }
      for (int n=0; n<nidata; ++n) {
        packed_pi[CheckedRestartSizeT(int_index)] = h_pi(n, p);
        AdvanceRestartOffset(int_index, 1, 1);
      }
    }

    std::vector<int> mb_counts_section;
    if (single_file_per_rank) {
      mb_counts_section = local_mb_counts;
    } else {
      mb_counts_section.assign(pm->nmb_total, 0);
      for (int m=0; m<nmb_local; ++m) {
        mb_counts_section[CheckedRestartIntAdd(gids_local, m)] = local_mb_counts[m];
      }
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, mb_counts_section.data(), pm->nmb_total,
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    }

    const int nmb_section = CheckedRestartIntCount(mb_counts_section.size());
    std::vector<IOWrapperSizeT> mb_offsets_section(
        CheckedRestartExtentPlusOne(nmb_section), 0);
    for (int m=0; m<nmb_section; ++m) {
      if (mb_counts_section[m] < 0) AbortRestartLayoutOverflow();
      mb_offsets_section[m + 1] = mb_offsets_section[m];
      AdvanceRestartOffset(mb_offsets_section[m + 1], mb_counts_section[m], 1);
    }
    const IOWrapperSizeT npart_section = mb_offsets_section[nmb_section];

    IOWrapperSizeT step5offset = step1size;
    AdvanceRestartOffset(step5offset, step2size, 1);
    AdvanceRestartOffset(step5offset, step3size, 1);
    AdvanceRestartOffset(step5offset, 1, sizeof(IOWrapperSizeT));
    if (single_file_per_rank) {
      AdvanceRestartOffset(step5offset, pm->nmb_thisrank, data_size);
    } else {
      AdvanceRestartOffset(step5offset, pm->nmb_total, data_size);
    }

    IOWrapperSizeT section_offset = step5offset;
    const IOWrapperSizeT magic_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(std::uint64_t));
    const IOWrapperSizeT version_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nmb_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nrdata_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nidata_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nout1_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nout2_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT nout3_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT has_mom_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT has_edge_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT mom_cnt_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT edge1_cnt_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT edge2_cnt_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT edge3_cnt_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT state_kind_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT physical_mode_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(int));
    const IOWrapperSizeT cr_light_speed_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(Real));
    const IOWrapperSizeT model_ints_offset = section_offset;
    AdvanceRestartOffset(section_offset, model_ints.size(), sizeof(int));
    const IOWrapperSizeT model_reals_offset = section_offset;
    AdvanceRestartOffset(section_offset, model_reals.size(), sizeof(Real));
    const IOWrapperSizeT npart_offset = section_offset;
    AdvanceRestartOffset(section_offset, 1, sizeof(IOWrapperSizeT));
    const IOWrapperSizeT mb_count_offset = section_offset;
    AdvanceRestartOffset(section_offset, nmb_section, sizeof(int));
    const IOWrapperSizeT pr_real_offset = section_offset;
    AdvanceRestartOffset(section_offset, npart_section,
                         CheckedRestartProduct(nrdata, sizeof(Real)));
    const IOWrapperSizeT pr_int_offset = section_offset;
    AdvanceRestartOffset(section_offset, npart_section,
                         CheckedRestartProduct(nidata, sizeof(int)));
    const IOWrapperSizeT moments_offset = section_offset;
    if (has_moments) {
      AdvanceRestartOffset(section_offset, nmb_section,
                           CheckedRestartProduct(moment_cnt, sizeof(Real)));
    }
    const IOWrapperSizeT edge1_offset = section_offset;
    if (has_edge) {
      AdvanceRestartOffset(section_offset, nmb_section,
                           CheckedRestartProduct(edge1_cnt, sizeof(Real)));
    }
    const IOWrapperSizeT edge2_offset = section_offset;
    if (has_edge) {
      AdvanceRestartOffset(section_offset, nmb_section,
                           CheckedRestartProduct(edge2_cnt, sizeof(Real)));
    }
    const IOWrapperSizeT edge3_offset = section_offset;
    if (has_edge) {
      AdvanceRestartOffset(section_offset, nmb_section,
                           CheckedRestartProduct(edge3_cnt, sizeof(Real)));
    }

    if (global_variable::my_rank == 0 || single_file_per_rank) {
      int i_has_mom = has_moments ? 1 : 0;
      int i_has_edge = has_edge ? 1 : 0;
      if (resfile.Write_any_type_at(&kPicMagic, sizeof(std::uint64_t), magic_offset,
                                    "byte", single_file_per_rank)
            != sizeof(std::uint64_t) ||
          resfile.Write_any_type_at(&kPicVersion, 1, version_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nmb_section, 1, nmb_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nrdata, 1, nrdata_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nidata, 1, nidata_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout1, 1, nout1_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout2, 1, nout2_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout3, 1, nout3_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&i_has_mom, 1, has_mom_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&i_has_edge, 1, has_edge_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&moment_cnt, 1, mom_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge1_cnt, 1, edge1_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge2_cnt, 1, edge2_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge3_cnt, 1, edge3_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&state_kind, 1, state_kind_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&physical_mode, 1, physical_mode_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&cr_light_speed, 1, cr_light_speed_offset, "Real",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(model_ints.data(), model_ints.size(),
                                    model_ints_offset, "int",
                                    single_file_per_rank) != model_ints.size() ||
          resfile.Write_any_type_at(model_reals.data(), model_reals.size(),
                                    model_reals_offset, "Real",
                                    single_file_per_rank) != model_reals.size() ||
          resfile.Write_any_type_at(&npart_section, sizeof(IOWrapperSizeT),
                                    npart_offset, "byte", single_file_per_rank)
            != sizeof(IOWrapperSizeT)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart metadata." << std::endl;
        restart_utils::AbortOnFatalError();
      }
      if (nmb_section > 0) {
        if (resfile.Write_any_type_at(mb_counts_section.data(), nmb_section,
                                      mb_count_offset, "int",
                                      single_file_per_rank) != nmb_section) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle MeshBlock counts."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
    }

#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    for (int m=0; m<nmb_local; ++m) {
      const int cnt = local_mb_counts[m];
      if (cnt <= 0) continue;
      const int section_m =
          (single_file_per_rank ? m : CheckedRestartIntAdd(gids_local, m));
      const IOWrapperSizeT gstart = mb_offsets_section[section_m];
      const IOWrapperSizeT lstart = local_mb_offsets[m];
      IOWrapperSizeT pr_off = pr_real_offset;
      IOWrapperSizeT pi_off = pr_int_offset;
      AdvanceRestartOffset(pr_off, gstart,
                           CheckedRestartProduct(nrdata, sizeof(Real)));
      AdvanceRestartOffset(pi_off, gstart,
                           CheckedRestartProduct(nidata, sizeof(int)));
      const IOWrapperSizeT pr_start = CheckedRestartProduct(lstart, nrdata);
      const IOWrapperSizeT pi_start = CheckedRestartProduct(lstart, nidata);
      const IOWrapperSizeT pr_count = CheckedRestartProduct(cnt, nrdata);
      const IOWrapperSizeT pi_count = CheckedRestartProduct(cnt, nidata);
      if (resfile.Write_any_type_at(&(packed_pr[CheckedRestartSizeT(pr_start)]),
                                    pr_count, pr_off, "Real",
                                    single_file_per_rank) != pr_count) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart real data." << std::endl;
        restart_utils::AbortOnFatalError();
      }
      if (resfile.Write_any_type_at(&(packed_pi[CheckedRestartSizeT(pi_start)]),
                                    pi_count, pi_off, "int",
                                    single_file_per_rank) != pi_count) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart integer data." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }

    if (has_moments) {
      auto h_mom = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->moments);
      for (int m=0; m<nmb_local; ++m) {
        const int section_m =
            (single_file_per_rank ? m : CheckedRestartIntAdd(gids_local, m));
        auto mom_mb = Kokkos::subview(h_mom, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                      Kokkos::ALL);
        if (mom_mb.size() != static_cast<std::size_t>(moment_cnt)) {
          AbortRestartLayoutOverflow();
        }
        IOWrapperSizeT moff = moments_offset;
        AdvanceRestartOffset(moff, section_m,
                             CheckedRestartProduct(moment_cnt, sizeof(Real)));
        if (resfile.Write_any_type_at(mom_mb.data(), moment_cnt, moff, "Real",
                                      single_file_per_rank) != moment_cnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle moment restart data." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
    }

    if (has_edge) {
      auto h_x1e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x1e);
      auto h_x2e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x2e);
      auto h_x3e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x3e);

      for (int m=0; m<nmb_local; ++m) {
        const int section_m =
            (single_file_per_rank ? m : CheckedRestartIntAdd(gids_local, m));
        auto x1_mb = Kokkos::subview(h_x1e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto x2_mb = Kokkos::subview(h_x2e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto x3_mb = Kokkos::subview(h_x3e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        if (x1_mb.size() != static_cast<std::size_t>(edge1_cnt) ||
            x2_mb.size() != static_cast<std::size_t>(edge2_cnt) ||
            x3_mb.size() != static_cast<std::size_t>(edge3_cnt)) {
          AbortRestartLayoutOverflow();
        }
        IOWrapperSizeT x1off = edge1_offset;
        IOWrapperSizeT x2off = edge2_offset;
        IOWrapperSizeT x3off = edge3_offset;
        AdvanceRestartOffset(x1off, section_m,
                             CheckedRestartProduct(edge1_cnt, sizeof(Real)));
        AdvanceRestartOffset(x2off, section_m,
                             CheckedRestartProduct(edge2_cnt, sizeof(Real)));
        AdvanceRestartOffset(x3off, section_m,
                             CheckedRestartProduct(edge3_cnt, sizeof(Real)));
        if (resfile.Write_any_type_at(x1_mb.data(), edge1_cnt, x1off, "Real",
                                      single_file_per_rank) != edge1_cnt ||
            resfile.Write_any_type_at(x2_mb.data(), edge2_cnt, x2off, "Real",
                                      single_file_per_rank) != edge2_cnt ||
            resfile.Write_any_type_at(x3_mb.data(), edge3_cnt, x3off, "Real",
                                      single_file_per_rank) != edge3_cnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle edge-current restart data." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
    }
  }

  // Sync, close, and publish only complete artifacts. Previous checkpoints use
  // distinct sequence names and remain available until this promotion succeeds.
  const int sync_status = resfile.Sync(single_file_per_rank);
  const int close_status = resfile.Close(single_file_per_rank);
  if (sync_status != 0 || close_status != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Failed to sync or close restart partial artifact '" << partial_fname
              << "'." << std::endl;
    restart_utils::AbortOnFatalError();
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  restart_utils::FileDigest local_digest;
  if (single_file_per_rank || global_variable::my_rank == 0) {
    restart_utils::PublishRestartArtifact(partial_fname, fname);
    local_digest = restart_utils::ComputeFileDigest(fname);
    restart_utils::WriteCompletionMarker(fname, local_digest);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  std::vector<std::pair<std::string, restart_utils::FileDigest>> members;
  if (single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
    std::array<std::uint64_t, 2> local_values = {
        local_digest.size, local_digest.fnv1a64};
    std::vector<std::uint64_t> gathered;
    if (global_variable::my_rank == 0) {
      gathered.resize(CheckedRestartSizeT(CheckedRestartProduct(
          static_cast<IOWrapperSizeT>(global_variable::nranks), 2)));
    }
    MPI_Gather(local_values.data(), 2, MPI_UINT64_T, gathered.data(), 2, MPI_UINT64_T,
               0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      for (int rank = 0; rank < global_variable::nranks; ++rank) {
        char rank_dir[20];
        std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", rank);
        restart_utils::FileDigest digest;
        const std::size_t digest_offset = CheckedRestartSizeT(
            CheckedRestartProduct(static_cast<IOWrapperSizeT>(rank), 2));
        digest.size = gathered[digest_offset];
        digest.fnv1a64 = gathered[digest_offset + 1];
        members.emplace_back(std::string("rst/") + rank_dir +
                             out_params.file_basename + number + ".rst", digest);
      }
    }
#else
    members.emplace_back(fname, local_digest);
#endif
  } else if (global_variable::my_rank == 0) {
    members.emplace_back(fname, local_digest);
  }
  if (global_variable::my_rank == 0) {
    restart_utils::WriteRestartManifest(manifest_fname, members);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return;
}
