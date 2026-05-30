//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>  // NOLINT(build/c++11)
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <cstring>
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
#include "mpi_utils.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "srcterms/turb_driver.hpp"
#include "outputs.hpp"
#include "output_file_utils.hpp"
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

struct RestartAttemptCleanup {
  std::string temporary_payload;
  std::string published_payload;
  std::string temporary_manifest;
  std::string published_manifest;
  std::string reservation;
  bool owns_temporary_payload = false;
  bool owns_published_payload = false;
  bool owns_temporary_manifest = false;
  bool owns_published_manifest = false;
  bool owns_reservation = false;
};

RestartAttemptCleanup *active_restart_cleanup = nullptr;

void CleanupActiveRestartAttempt() {
  if (active_restart_cleanup == nullptr) return;
  if (active_restart_cleanup->owns_temporary_payload) {
    output_file_utils::DiscardOwnedPath(active_restart_cleanup->temporary_payload);
  }
  if (active_restart_cleanup->owns_published_payload) {
    output_file_utils::DiscardOwnedPath(active_restart_cleanup->published_payload);
  }
  if (active_restart_cleanup->owns_temporary_manifest) {
    output_file_utils::DiscardOwnedPath(active_restart_cleanup->temporary_manifest);
  }
  if (active_restart_cleanup->owns_published_manifest) {
    output_file_utils::DiscardOwnedPath(active_restart_cleanup->published_manifest);
  }
  if (active_restart_cleanup->owns_reservation) {
    output_file_utils::DiscardOwnedPath(active_restart_cleanup->reservation);
  }
}

[[noreturn]] void FailNodeRestartWrite(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" + message);
}

[[noreturn]] void FailNodeRestartWriteCoordinated(const std::string &message) {
  CleanupActiveRestartAttempt();
#if MPI_PARALLEL_ENABLED
  int barrier_error = MPI_Barrier(MPI_COMM_WORLD);
  if (barrier_error != MPI_SUCCESS) {
    std::cerr << "MPI_Barrier for coordinated node-restart cleanup failed with MPI "
              << "error: " << mpi_utils::MpiErrorString(barrier_error) << std::endl;
  }
#endif
  active_restart_cleanup = nullptr;
  mpi_utils::SetFatalCleanupHook(nullptr);
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" + message);
}

void RecordNodeRestartPayloadFailure(bool coordinate_failure, int &local_failure,
                                     std::string &local_error,
                                     const std::string &message) {
  if (!coordinate_failure) {
    FailNodeRestartWrite(message);
  }
  local_failure = 1;
  if (local_error.empty()) {
    local_error = message;
  }
}

std::uint64_t NodeRestartGenerationSeed() {
  const char *configured_seed = std::getenv("ATHENAK_TEST_NODE_RESTART_GENERATION");
  if (configured_seed == nullptr) {
    return static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }
  char *end = nullptr;
  errno = 0;
  auto value = std::strtoull(configured_seed, &end, 10);
  if (errno != 0 || end == configured_seed || *end != '\0') {
    FailNodeRestartWrite("ATHENAK_TEST_NODE_RESTART_GENERATION must be an unsigned "
                         "integer.");
  }
  return static_cast<std::uint64_t>(value);
}

bool InjectNodeRestartFailure(const std::string &stage) {
  const char *configured_stage = std::getenv("ATHENAK_TEST_NODE_RESTART_FAIL_STAGE");
  return configured_stage != nullptr && stage == configured_stage;
}

bool AnyWorldFailure(int local_failure, const std::string &context) {
#if MPI_PARALLEL_ENABLED
  int any_failure = 0;
  mpi_utils::CheckMpi(MPI_Allreduce(&local_failure, &any_failure, 1, MPI_INT, MPI_MAX,
                                    MPI_COMM_WORLD), context.c_str());
  return any_failure != 0;
#else
  return local_failure != 0;
#endif
}

bool BroadcastRootFailure(int root_failure, const std::string &context) {
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(MPI_Bcast(&root_failure, 1, MPI_INT, 0, MPI_COMM_WORLD),
                      context.c_str());
#endif
  return root_failure != 0;
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

bool ReserveRestartGeneration(const std::string &reservation_name) {
  int descriptor = open(reservation_name.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0600);
  if (descriptor >= 0) {
    if (close(descriptor) != 0) {
      output_file_utils::DiscardOwnedPath(reservation_name);
      FailNodeRestartWrite("Could not close node restart generation reservation '" +
                           reservation_name + "'.");
    }
    return true;
  }
  if (errno == EEXIST) {
    return false;
  }
  FailNodeRestartWrite("Could not reserve node restart generation '" + reservation_name +
                       "': " + std::strerror(errno) + ".");
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  output_file_utils::EnsureDirectory("rst", 0775, "restart output",
                                     FailNodeRestartWrite);
  if (IsSharded(op.shard_mode)) {
    std::string shard_dir = "rst/" + ShardDirectoryName(
        op.shard_mode, global_variable::my_rank, global_variable::node_id);
    output_file_utils::EnsureDirectory(shard_dir, 0775, "restart output",
                                       FailNodeRestartWrite);
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
  std::string published_fname;
  std::string manifest_name;
  std::string payload_name;
  std::string reservation_name;
  std::uint64_t generation = 0;
  RestartAttemptCleanup cleanup;
  active_restart_cleanup = &cleanup;
  mpi_utils::SetFatalCleanupHook(CleanupActiveRestartAttempt);
  std::string number = "." + output_file_utils::FormatSequence(
      out_params.file_number, "restart output", FailNodeRestartWrite);
  if (IsRankSharded(shard_mode)) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    published_fname = std::string("rst/") + ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id) + "/"
      + out_params.file_basename
      + number + ".rst";
    fname = output_file_utils::TemporaryPath(published_fname);
    cleanup.temporary_payload = fname;
    cleanup.published_payload = published_fname;
    cleanup.owns_temporary_payload = shard_writer;
  } else if (node_sharded) {
    if (!IsSafeRestartLeaf(out_params.file_basename)) {
      FailNodeRestartWrite("Node restart output requires <job>/basename to be a "
                           "single safe path component.");
    }
    if (global_variable::my_rank == 0) {
      generation = NodeRestartGenerationSeed();
      while (true) {
        payload_name = out_params.file_basename + number + ".g"
            + std::to_string(generation) + ".payload.rst";
        for (int id = 0; id < global_variable::nnodes; ++id) {
          NodeRestartRelativePayloadPath(payload_name, id);
        }
        reservation_name = std::string("rst/.") + out_params.file_basename + number
            + ".g" + std::to_string(generation) + ".reserve";
        if (ReserveRestartGeneration(reservation_name)) {
          cleanup.reservation = reservation_name;
          cleanup.owns_reservation = true;
          break;
        }
        if (generation == std::numeric_limits<std::uint64_t>::max()) {
          FailNodeRestartWrite("Could not reserve a node restart generation: the "
                               "generation counter is exhausted.");
        }
        ++generation;
      }
    }
#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(MPI_Bcast(&generation, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD),
                        "MPI_Bcast for node restart generation reservation");
#endif
    reservation_name = std::string("rst/.") + out_params.file_basename + number
        + ".g" + std::to_string(generation) + ".reserve";
    payload_name = out_params.file_basename + number + ".g"
        + std::to_string(generation) + ".payload.rst";
    for (int id = 0; id < global_variable::nnodes; ++id) {
      NodeRestartRelativePayloadPath(payload_name, id);
    }
    std::string shard_dir = ShardDirectoryName(
        shard_mode, global_variable::my_rank, global_variable::node_id);
    published_fname = std::string("rst/") + shard_dir + "/" + payload_name;
    fname = output_file_utils::TemporaryPath(published_fname);
    manifest_name = std::string("rst/") + out_params.file_basename + number + ".rst";
    cleanup.temporary_payload = fname;
    cleanup.published_payload = published_fname;
    cleanup.reservation = reservation_name;
    cleanup.owns_temporary_payload = global_variable::node_rank == 0;
    cleanup.owns_reservation = global_variable::my_rank == 0;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    published_fname = std::string("rst/") + out_params.file_basename + number + ".rst";
    fname = output_file_utils::TemporaryPath(published_fname);
    cleanup.temporary_payload = fname;
    cleanup.published_payload = published_fname;
    cleanup.owns_temporary_payload = shard_writer;
  }
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number = output_file_utils::AdvanceFileNumber(
      out_params.file_number, "restart output", FailNodeRestartWrite);
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
  int local_payload_io_failure = 0;
  std::string local_payload_io_error;
  auto checked_payload_write = [&](const void *data, IOWrapperSizeT bytes,
                                   const std::string &context) {
    if (resfile.Write_any_type(data, bytes, "byte", independent_file) != bytes) {
      RecordNodeRestartPayloadFailure(node_sharded, local_payload_io_failure,
                                      local_payload_io_error,
                                      context + " was not written completely.");
    }
  };
  if (shard_writer) {
    // output the input parameters (input file)
    checked_payload_write(sbuf.c_str(), sbuf.size(), "restart parameter dump");
    if (node_sharded) {
      checked_payload_write(kNodeRestartPayloadMarker, kNodeRestartPayloadMarkerSize,
                            "node restart payload marker");
    }

    // output Mesh information
    checked_payload_write(&(pm->nmb_total), sizeof(int), "restart total MeshBlock count");
    checked_payload_write(&(pm->root_level), sizeof(int), "restart root level");
    checked_payload_write(&(pm->mesh_size), sizeof(RegionSize), "restart mesh size");
    checked_payload_write(&(pm->mesh_indcs), sizeof(RegionIndcs), "restart mesh indices");
    checked_payload_write(&(pm->mb_indcs), sizeof(RegionIndcs),
                          "restart MeshBlock indices");
    checked_payload_write(&(pm->time), sizeof(Real), "restart time");
    checked_payload_write(&(pm->dt), sizeof(Real), "restart timestep");
    checked_payload_write(&(pm->ncycle), sizeof(int), "restart cycle");
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (shard_writer) {
    checked_payload_write(&(pm->lloc_eachmb[0]),
                        CheckedRestartMultiply(
                            CheckedRestartCount(pm->nmb_total,
                                                "restart total MeshBlock count"),
                            sizeof(LogicalLocation),
                                               "restart logical-location bytes"),
                        "restart logical locations");
    checked_payload_write(&(pm->cost_eachmb[0]),
                        CheckedRestartMultiply(
                            CheckedRestartCount(pm->nmb_total,
                                                "restart total MeshBlock count"),
                            sizeof(float),
                                               "restart MeshBlock-cost bytes"),
                        "restart MeshBlock costs");
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (shard_writer) {
    // store z4c information
    if (pz4c != nullptr) {
      checked_payload_write(&(pz4c->last_output_time), sizeof(Real),
                            "restart z4c output time");
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        checked_payload_write(pt->GetPos(), 3*sizeof(Real), "restart puncture position");
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      checked_payload_write(&(pturb->rstate), sizeof(RNG_State),
                            "restart turbulence RNG state");
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
    checked_payload_write(&(data_size), sizeof(IOWrapperSizeT),
                          "restart per-MeshBlock byte count");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered hydro data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered hydro data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered mhd data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered mhd data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x1f data not written correctly to rst file, restart file is broken.");
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x1f_bytes,
                                     "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestartWrite, "restart MHD x2-face subview count");
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          independent_file) != fldcnt) {
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x2f data not written correctly to rst file, restart file is broken.");
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x2f_bytes,
                                     "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestartWrite, "restart MHD x3-face subview count");
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          independent_file) != fldcnt) {
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x3f data not written correctly to rst file, restart file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x1f data not written correctly to rst file, restart file is broken.");
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x1f_bytes,
                                     "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestartWrite, "restart MHD x2-face subview count");
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      independent_file) != fldcnt) {
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x2f data not written correctly to rst file, restart file is broken.");
        }
        myoffset = CheckedRestartAdd(myoffset, layout.mhd_x2f_bytes,
                                     "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestartWrite, "restart MHD x3-face subview count");
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      independent_file) != fldcnt) {
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "b0.x3f data not written correctly to rst file, restart file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered rad data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered rad data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered turb data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered turb data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered z4c data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered z4c data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered adm data not written correctly to rst file, restart "
              "file is broken.");
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
          RecordNodeRestartPayloadFailure(
              node_sharded, local_payload_io_failure, local_payload_io_error,
              "cell-centered adm data not written correctly to rst file, restart "
              "file is broken.");
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
    RecordNodeRestartPayloadFailure(
        node_sharded, local_payload_io_failure, local_payload_io_error,
        "restart payload could not be closed cleanly.");
  }
  if (node_sharded && InjectNodeRestartFailure("after_payload_write") &&
      global_variable::my_rank == 0) {
    RecordNodeRestartPayloadFailure(
        true, local_payload_io_failure, local_payload_io_error,
        "Injected node restart failure after payload write.");
  }
  if (node_sharded &&
      AnyWorldFailure(local_payload_io_failure,
                      "MPI_Allreduce for node restart payload IO completion")) {
    FailNodeRestartWriteCoordinated(local_payload_io_error.empty()
        ? "A node restart payload could not be written completely."
        : local_payload_io_error);
  }

  if (!node_sharded) {
    if (shard_writer) {
      output_file_utils::PublishTemporaryFile(
          fname, published_fname, "restart output", FailNodeRestartWrite);
      cleanup.owns_temporary_payload = false;
    }
#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                        "MPI_Barrier after restart output publication");
#endif
    active_restart_cleanup = nullptr;
    mpi_utils::SetFatalCleanupHook(nullptr);
    return;
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
    int local_payload_failure = 0;
    std::string local_payload_error;
    if (global_variable::node_rank == 0) {
      std::ifstream payload_check(fname, std::ios::binary | std::ios::ate);
      IOWrapperSizeT observed_size = payload_check.good()
          ? static_cast<IOWrapperSizeT>(payload_check.tellg()) : 0;
      payload_check.close();
      if (observed_size != expected_size) {
        local_payload_failure = 1;
        local_payload_error = "Node restart payload '" + fname +
            "' could not be published: expected " + std::to_string(expected_size) +
            " bytes and found " + std::to_string(observed_size) + ".";
      } else {
        std::string publish_error;
        if (!output_file_utils::TryPublishTemporaryFile(
                fname, published_fname, "node restart payload", &publish_error)) {
          local_payload_failure = 1;
          local_payload_error = publish_error;
        } else {
          cleanup.owns_temporary_payload = false;
          cleanup.owns_published_payload = true;
        }
      }
    }
    if (AnyWorldFailure(local_payload_failure,
                        "MPI_Allreduce for node restart payload publication")) {
      FailNodeRestartWriteCoordinated(local_payload_error.empty()
          ? "A node restart payload could not be published."
          : local_payload_error);
    }

#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                        "MPI_Barrier after node restart payload publication");
#endif
    if (AnyWorldFailure(InjectNodeRestartFailure("after_payload_publication"),
                        "MPI_Allreduce for injected node restart failure")) {
      FailNodeRestartWriteCoordinated(
          "Injected node restart failure after payload publication.");
    }
    std::vector<int> manifest_nodes;
    std::vector<int> manifest_offsets;
    if (global_variable::my_rank == 0) {
      manifest_nodes.resize(global_variable::nranks);
      manifest_offsets.resize(global_variable::nranks);
    }
#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(
        MPI_Gather(&(global_variable::node_id), 1, MPI_INT,
                   global_variable::my_rank == 0 ? manifest_nodes.data() : nullptr,
                   1, MPI_INT, 0, MPI_COMM_WORLD),
        "MPI_Gather for node restart manifest node IDs");
    mpi_utils::CheckMpi(
        MPI_Gather(&payload_block_offset, 1, MPI_INT,
                   global_variable::my_rank == 0 ? manifest_offsets.data() : nullptr,
                   1, MPI_INT, 0, MPI_COMM_WORLD),
        "MPI_Gather for node restart manifest payload offsets");
#else
    manifest_nodes[0] = global_variable::node_id;
    manifest_offsets[0] = payload_block_offset;
#endif

    std::vector<int> blocks_per_node;
    std::string manifest_error;
    if (global_variable::my_rank == 0) {
      blocks_per_node.assign(global_variable::nnodes, 0);
      std::vector<int> next_payload_block(global_variable::nnodes, 0);
      int next_gid = 0;
      for (int r = 0; r < global_variable::nranks; ++r) {
        int id = manifest_nodes[r];
        int blocks = pm->nmb_eachrank[r];
        if (id < 0 || id >= global_variable::nnodes || blocks < 0 ||
            pm->gids_eachrank[r] != next_gid ||
            manifest_offsets[r] != next_payload_block[id]) {
          manifest_error = "Node restart segment map is inconsistent and cannot "
                           "be published.";
          break;
        }
        if (blocks > std::numeric_limits<int>::max() - blocks_per_node[id] ||
            blocks > std::numeric_limits<int>::max() - next_payload_block[id] ||
            blocks > std::numeric_limits<int>::max() - next_gid) {
          manifest_error = "Node restart segment map exceeds INT_MAX and cannot "
                           "be published.";
          break;
        }
        blocks_per_node[id] += blocks;
        next_payload_block[id] += blocks;
        next_gid += blocks;
      }
      if (manifest_error.empty() && next_gid != pm->nmb_total) {
        manifest_error = "Node restart segment map does not cover all mesh blocks.";
      }
    }
    if (BroadcastRootFailure(
            global_variable::my_rank == 0 && !manifest_error.empty(),
            "MPI_Bcast for node restart segment-map validation")) {
      FailNodeRestartWriteCoordinated(manifest_error.empty()
          ? "Node restart segment-map validation failed."
          : manifest_error);
    }

    if (global_variable::my_rank == 0) {
      std::string temporary_manifest =
          manifest_name + ".tmp.g" + std::to_string(generation);
      cleanup.temporary_manifest = temporary_manifest;
      cleanup.owns_temporary_manifest = true;
      std::ofstream manifest(temporary_manifest, std::ios::trunc);
      restart_layout::ManifestBudget manifest_budget{0, kMaxNodeRestartManifestBytes};
      const auto write_manifest_record = [&](const std::string &record) {
        if (!manifest_error.empty()) return;
        IOWrapperSizeT record_bytes = static_cast<IOWrapperSizeT>(record.size());
        if (record.size() > std::numeric_limits<IOWrapperSizeT>::max() ||
            record_bytes >= kMaxNodeRestartManifestBytes - manifest_budget.bytes) {
          manifest_error = "Node restart manifest exceeds the " +
              std::to_string(kMaxNodeRestartManifestBytes) + "-byte limit.";
          return;
        }
        manifest_budget.bytes += record_bytes + 1;
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
      if (manifest_error.empty() &&
          (!manifest.good() || observed_manifest_bytes < 0 ||
           static_cast<IOWrapperSizeT>(observed_manifest_bytes) !=
               manifest_budget.bytes)) {
        manifest_error = "Node restart manifest '" + manifest_name +
                         "' was not written completely.";
      }
      if (manifest_error.empty()) {
        std::string publish_error;
        if (!output_file_utils::TryPublishTemporaryFile(
                temporary_manifest, manifest_name, "node restart manifest",
                &publish_error)) {
          manifest_error = publish_error;
        } else {
          cleanup.owns_temporary_manifest = false;
          cleanup.published_manifest = manifest_name;
          cleanup.owns_published_manifest = true;
        }
      }
    }
    if (BroadcastRootFailure(
            global_variable::my_rank == 0 && !manifest_error.empty(),
            "MPI_Bcast for node restart manifest publication")) {
      FailNodeRestartWriteCoordinated(manifest_error.empty()
          ? "Node restart manifest publication failed."
          : manifest_error);
    }
    int reservation_remove_failure = 0;
    std::string reservation_remove_error;
    if (global_variable::my_rank == 0) {
      if (std::remove(reservation_name.c_str()) != 0) {
        reservation_remove_failure = 1;
        reservation_remove_error =
            "Could not remove node restart generation reservation '" +
            reservation_name + "'.";
      } else {
        cleanup.owns_reservation = false;
      }
    }
    if (BroadcastRootFailure(
            reservation_remove_failure,
            "MPI_Bcast for node restart reservation removal")) {
      FailNodeRestartWriteCoordinated(reservation_remove_error.empty()
          ? "Node restart reservation removal failed."
          : reservation_remove_error);
    }
#if MPI_PARALLEL_ENABLED
    mpi_utils::CheckMpi(MPI_Barrier(MPI_COMM_WORLD),
                        "MPI_Barrier after node restart manifest publication");
#endif
    cleanup.owns_published_payload = false;
    cleanup.owns_published_manifest = false;
  }

  active_restart_cleanup = nullptr;
  mpi_utils::SetFatalCleanupHook(nullptr);
  return;
}
