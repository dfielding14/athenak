//========================================================================================
// Athena astrophysical MHD code (Kokkos version)
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//////////////////////////////////// Athena Main Program /////////////////////////////////
//! \file main.cpp
//! \brief Athena main program
//!
//! Based on the Athena (Cambridge version) and Athena++ MHD codes. Athena originally was
//! written in 2002-2005 by Jim Stone, Tom Gardiner, and Peter Teuben, with many important
//! contributions by many other developers after that, i.e. 2005-2014.
//!
//! Athena++ was started in Jan 2014, with the core design developed 4-7/2014 during an
//! extended visit to the KITP at UCSB by J. Stone. GR was implemented by Chris White and
//! AMR by Kengo Tomida 2014-2016, with contributions from many others (esp. K. Felker)
//! continuing after that.
//!
//! Athena (Kokkos version) is an outgrowth of the Athena-Parthenon collaboration, and is
//! a completely new implementation based on the Kokkos performance-portability library
//! (which is an external dependency required for this version). It was started 6/2020
//! during the pandemic. As part of the keep-it-simple design, only a fraction of the
//! features of the C++ version are implemented.
//========================================================================================

// C/C++ headers
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdio> // sscanf
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>

// Athena headers
#include "athena.hpp"
#include "file_sharding.hpp"
#include "globals.hpp"
#include "utils/utils.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "outputs/outputs.hpp"
#include "driver/driver.hpp"
#include "particles/particles.hpp"

// MPI/OpenMP headers
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#if OPENMP_PARALLEL_ENABLED
#include <omp.h>
#endif

#if defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace {

namespace fs = std::filesystem;

struct RestartPathInfo {
  FileShardMode shard_mode = FileShardMode::shared;
  std::string normalized_open_path;
  std::string base_dir;
  std::string file_name;
};

constexpr int kConversionFields = 8;
constexpr long double kConversionTolerance = 5.0e-13L;
const char *kConversionFieldNames[kConversionFields] = {
  "density", "momentum_x1", "momentum_x2", "momentum_x3",
  "total_energy", "metallicity", "dust_small", "dust_large"
};

struct ConversionDiagnostics {
  long double source[kConversionFields] = {};
  long double source_abs[kConversionFields] = {};
  long double target[kConversionFields] = {};
  long double source_kinetic = 0.0;
  long double source_internal = 0.0;
  long double target_kinetic = 0.0;
  long double target_internal = 0.0;
  long double max_scaled_integral_error = 0.0;
  long long invalid_target_cells = 0;
};

bool HasShardPrefix(const std::string &name, const char *prefix) {
  return name.rfind(prefix, 0) == 0;
}

bool IsExactShardDirectoryName(const std::string &name, const char *prefix) {
  constexpr int kShardIdDigits = 8;
  std::size_t prefix_len = std::char_traits<char>::length(prefix);
  if (!HasShardPrefix(name, prefix) || name.size() != prefix_len + kShardIdDigits) {
    return false;
  }
  for (std::size_t i = prefix_len; i < name.size(); ++i) {
    if (!std::isdigit(static_cast<unsigned char>(name[i]))) {
      return false;
    }
  }
  return true;
}

std::string PathToString(const fs::path &path) {
  return path.empty() ? std::string() : path.generic_string();
}

RestartPathInfo NormalizeRestartPath(const std::string &restart_arg) {
  if (restart_arg.empty()) {
    throw std::runtime_error("Restart path is empty.");
  }

  fs::path normalized = fs::path(restart_arg).lexically_normal();
  if (normalized.filename().empty()) {
    throw std::runtime_error("Restart path must include a file name.");
  }

  FileShardMode embedded_mode = FileShardMode::shared;
  std::string embedded_name;
  for (const auto &part : normalized.parent_path()) {
    std::string name = part.string();
    if (name.empty() || name == "." || name == ".." || name == "/") {
      continue;
    }

    bool looks_like_node = HasShardPrefix(name, "node_");
    bool looks_like_rank = HasShardPrefix(name, "rank_");
    if (!looks_like_node && !looks_like_rank) {
      continue;
    }

    bool exact_node = IsExactShardDirectoryName(name, "node_");
    bool exact_rank = IsExactShardDirectoryName(name, "rank_");
    if (!exact_node && !exact_rank) {
      throw std::runtime_error("Malformed restart shard directory '" + name
                               + "'. Expected node_######## or rank_########.");
    }

    FileShardMode part_mode = exact_node ? FileShardMode::per_node : FileShardMode::per_rank;
    if (embedded_mode != FileShardMode::shared) {
      throw std::runtime_error("Restart path may contain only one shard directory.");
    }
    embedded_mode = part_mode;
    embedded_name = name;
  }

  RestartPathInfo info;
  info.shard_mode = embedded_mode;
  info.file_name = normalized.filename().string();

  std::error_code ec;
  if (embedded_mode != FileShardMode::shared) {
    std::string parent_name = normalized.parent_path().filename().string();
    if (parent_name != embedded_name) {
      throw std::runtime_error("Restart shard directory must be the direct parent of the "
                               "restart file.");
    }

    fs::path base_dir = normalized.parent_path().parent_path();
    info.base_dir = PathToString(base_dir);

    if (!fs::exists(normalized, ec) || ec) {
      throw std::runtime_error("Restart shard file '" + PathToString(normalized)
                               + "' does not exist.");
    }

    if (embedded_mode == FileShardMode::per_node) {
      fs::path manifest_path = base_dir / info.file_name;
      info.normalized_open_path = PathToString(manifest_path);
      ec.clear();
      if (!fs::exists(manifest_path, ec) || ec) {
        throw std::runtime_error("Per-node restart manifest '" + info.normalized_open_path
                                 + "' does not exist.");
      }
    } else {
      info.normalized_open_path = PathToString(normalized);
    }
    return info;
  }

  fs::path base_dir = normalized.parent_path();
  info.base_dir = PathToString(base_dir);
  info.normalized_open_path = PathToString(normalized);

  fs::path search_dir = base_dir.empty() ? fs::path(".") : base_dir;
  bool found_payload = false;
  ec.clear();
  if (fs::exists(search_dir, ec) && !ec && fs::is_directory(search_dir, ec)) {
    for (fs::directory_iterator it(search_dir, ec); !ec && it != fs::directory_iterator();
         it.increment(ec)) {
      if (ec) break;
      if (!it->is_directory(ec) || ec) {
        ec.clear();
        continue;
      }

      std::string name = it->path().filename().string();
      if (!IsExactShardDirectoryName(name, "node_")) {
        continue;
      }

      fs::path payload_path = it->path() / info.file_name;
      std::error_code exists_ec;
      if (fs::exists(payload_path, exists_ec) && !exists_ec) {
        found_payload = true;
        break;
      }
    }
  }

  if (found_payload) {
    info.shard_mode = FileShardMode::per_node;
    ec.clear();
    if (!fs::exists(normalized, ec) || ec) {
      throw std::runtime_error("Per-node restart manifest '" + info.normalized_open_path
                               + "' does not exist.");
    }
    return info;
  }

  ec.clear();
  if (!fs::exists(normalized, ec) || ec) {
    throw std::runtime_error("Restart file '" + info.normalized_open_path
                             + "' does not exist.");
  }

  return info;
}

ConversionDiagnostics MeasureRestartConversion(Mesh *mesh) {
  ConversionDiagnostics diagnostics;
  auto *pack = mesh->pmb_pack;
  auto &indcs = mesh->mb_indcs;
  auto target = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), Kokkos::subview(pack->phydro->u0,
                                     std::make_pair(0, pack->nmb_thispack),
                                     Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL));

  for (int n=0; n<kConversionFields; ++n) {
    diagnostics.source[n] = mesh->restart_meta.source_integral[n];
    diagnostics.source_abs[n] = mesh->restart_meta.source_abs_integral[n];
  }
  diagnostics.source_kinetic = mesh->restart_meta.source_kinetic;
  diagnostics.source_internal = mesh->restart_meta.source_internal;

  for (int m=0; m<pack->nmb_thispack; ++m) {
    int gid = pack->gids + m;
    long double cell_volume = static_cast<long double>(mesh->mesh_size.dx1) *
                              static_cast<long double>(mesh->mesh_size.dx2) *
                              static_cast<long double>(mesh->mesh_size.dx3);
    for (int level=mesh->root_level; level<mesh->lloc_eachmb[gid].level; ++level) {
      cell_volume *= 0.125L;
    }
    for (int k=indcs.ks; k<=indcs.ke; ++k) {
      for (int j=indcs.js; j<=indcs.je; ++j) {
        for (int i=indcs.is; i<=indcs.ie; ++i) {
          long double value[kConversionFields];
          bool finite = true;
          for (int n=0; n<kConversionFields; ++n) {
            value[n] = static_cast<long double>(target(m, n, k, j, i));
            diagnostics.target[n] += value[n]*cell_volume;
            finite = finite && std::isfinite(value[n]);
          }
          long double kinetic = 0.5L*(value[1]*value[1] + value[2]*value[2] +
                                      value[3]*value[3])/value[0];
          long double internal = value[4] - kinetic;
          diagnostics.target_kinetic += kinetic*cell_volume;
          diagnostics.target_internal += internal*cell_volume;
          if (!finite || value[0] <= 0.0 || internal <= 0.0 ||
              !std::isfinite(kinetic) || !std::isfinite(internal)) {
            ++diagnostics.invalid_target_cells;
          }
        }
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, diagnostics.source, kConversionFields, MPI_LONG_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, diagnostics.source_abs, kConversionFields, MPI_LONG_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, diagnostics.target, kConversionFields, MPI_LONG_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  long double energies[4] = {diagnostics.source_kinetic,
                             diagnostics.source_internal,
                             diagnostics.target_kinetic,
                             diagnostics.target_internal};
  MPI_Allreduce(MPI_IN_PLACE, energies, 4, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  diagnostics.source_kinetic = energies[0];
  diagnostics.source_internal = energies[1];
  diagnostics.target_kinetic = energies[2];
  diagnostics.target_internal = energies[3];
  MPI_Allreduce(MPI_IN_PLACE, &diagnostics.invalid_target_cells, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int n=0; n<kConversionFields; ++n) {
    long double scale = std::max(1.0L, diagnostics.source_abs[n]);
    diagnostics.max_scaled_integral_error = std::max(
        diagnostics.max_scaled_integral_error,
        std::fabs(diagnostics.target[n] - diagnostics.source[n])/scale);
  }
  return diagnostics;
}

void WriteRestartConversion(Driver *driver, Outputs *outputs, Mesh *mesh,
                            ParameterInput *pin, const fs::path &output_dir) {
  int floor_count = mesh->ecounter.neos_dfloor + mesh->ecounter.neos_efloor +
                    mesh->ecounter.neos_tfloor;
  // Populate target ghosts and primitives through the ordinary restart path, without
  // entering Driver::Initialize(), which would replace the stored timestep.
  driver->InitBoundaryValuesAndPrimitives(mesh);
  int new_floor_count = mesh->ecounter.neos_dfloor + mesh->ecounter.neos_efloor +
                        mesh->ecounter.neos_tfloor;
  if (new_floor_count != floor_count) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Primitive recovery applied a floor to the converted "
              << "state; no restart was written." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  ConversionDiagnostics diagnostics = MeasureRestartConversion(mesh);
  if (diagnostics.invalid_target_cells != 0 ||
      diagnostics.max_scaled_integral_error > kConversionTolerance) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Restart conversion failed its finite-state or conserved-"
              << "integral acceptance gate; no restart was written." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  BaseTypeOutput *fluid_output = nullptr;
  BaseTypeOutput *particle_output = nullptr;
  for (auto *output : outputs->pout_list) {
    if (output->out_params.file_type == "rst") {
      fluid_output = output;
    } else if (output->out_params.file_type == "rst_prtcl") {
      particle_output = output;
    }
  }
  if (fluid_output == nullptr || particle_output == nullptr ||
      fluid_output->out_params.file_shard_mode != FileShardMode::per_node) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Restart conversion requires one per-node fluid restart "
              << "output and one particle restart output." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string basename = pin->GetString("job", "basename");
  if (mesh->restart_meta.file_name.rfind(basename + ".", 0) == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Converted output basename must differ from the archived "
              << "source basename." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  fluid_output->out_params.file_number = 0;
  particle_output->out_params.file_number = 0;
  pin->SetInteger(fluid_output->out_params.block_name, "file_number", 1);
  pin->SetInteger(particle_output->out_params.block_name, "file_number", 1);
  pin->SetBoolean("restart_conversion", "enabled", false);

  fs::path fluid_path = output_dir / "rst" / (basename + ".00000.rst");
  fs::path particle_path = output_dir / "rst_prtcl" /
                           (basename + ".00000.rst_prtcl");
  pin->SetString("particles", "star_particle_file", PathToString(particle_path));

  // The sidecar is complete before its absolute path is embedded in the fluid restart.
  particle_output->LoadOutputData(mesh);
  Kokkos::fence();
  particle_output->WriteOutputFile(mesh, pin, true);
  Kokkos::fence();
  fluid_output->LoadOutputData(mesh);
  Kokkos::fence();
  fluid_output->WriteOutputFile(mesh, pin, true);
  Kokkos::fence();
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (global_variable::my_rank == 0) {
    std::ofstream manifest("conversion_manifest.tmp", std::ios::trunc);
    manifest << std::setprecision(21);
    manifest << "source_fluid=" << mesh->restart_meta.source_restart_path << '\n'
             << "source_particles=" << mesh->restart_meta.source_particle_path << '\n'
             << "converted_fluid=" << PathToString(fluid_path) << '\n'
             << "converted_particles=" << PathToString(particle_path) << '\n'
             << "output_shards=" << global_variable::nnodes << '\n'
             << "source_blocks=" << mesh->restart_meta.source_nmb_total << '\n'
             << "target_blocks=" << mesh->nmb_total << '\n'
             << "source_max_physical_level="
             << mesh->restart_meta.source_max_level - mesh->root_level << '\n'
             << "target_max_physical_level="
             << mesh->restart_meta.target_max_level - mesh->root_level << '\n'
             << "time=" << mesh->time << '\n'
             << "cycle=" << mesh->ncycle << '\n'
             << "dt=" << mesh->dt << '\n'
             << "timestep_policy=preserved from source restart\n"
             << "physical_source_invocations=0\n"
             << "particle_updates=0\n"
             << "primitive_recovery_floors=0\n"
             << "invalid_target_cells=" << diagnostics.invalid_target_cells << '\n'
             << "conservation_tolerance=" << kConversionTolerance << '\n'
             << "max_scaled_integral_error="
             << diagnostics.max_scaled_integral_error << '\n'
             << "turbulence_rng=preserved\n"
             << "turbulence_force=zeroed for inactive legacy driver\n"
             << "hashes=pending\n";
    for (int n=0; n<kConversionFields; ++n) {
      manifest << "source_integral." << kConversionFieldNames[n] << '='
               << diagnostics.source[n] << '\n'
               << "source_abs_integral." << kConversionFieldNames[n] << '='
               << diagnostics.source_abs[n] << '\n'
               << "target_integral." << kConversionFieldNames[n] << '='
               << diagnostics.target[n] << '\n';
    }
    manifest << "source_resolved_kinetic=" << diagnostics.source_kinetic << '\n'
             << "target_resolved_kinetic=" << diagnostics.target_kinetic << '\n'
             << "resolved_kinetic_change="
             << diagnostics.target_kinetic - diagnostics.source_kinetic << '\n'
             << "source_internal_energy=" << diagnostics.source_internal << '\n'
             << "target_internal_energy=" << diagnostics.target_internal << '\n'
             << "internal_energy_change="
             << diagnostics.target_internal - diagnostics.source_internal << '\n'
             << "energy_partition_closure="
             << (diagnostics.target_kinetic - diagnostics.source_kinetic) +
                (diagnostics.target_internal - diagnostics.source_internal) << '\n';
    manifest.close();
    if (!manifest || std::rename("conversion_manifest.tmp", "conversion_manifest.txt") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Unable to publish conversion manifest." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // Hashing every shard is a conversion-workflow concern.  Publish only a written
    // marker here; the focused finalizer promotes it to CONVERSION_COMPLETE.
    std::ofstream marker("CONVERSION_WRITTEN.tmp", std::ios::trunc);
    marker << "Restart pair and conversion manifest written successfully.\n";
    marker.close();
    if (!marker || std::rename("CONVERSION_WRITTEN.tmp", "CONVERSION_WRITTEN") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Unable to publish conversion written marker."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn int main(int argc, char *argv[])
//! \brief Athena main program

int main(int argc, char *argv[]) {
  std::string input_file, restart_file, run_dir;
  std::string restart_base_dir, restart_file_name;
  std::string source_restart_path, source_particle_path, conversion_output_dir;
  bool iarg_flag = false;  // set to true if -i <file> argument is on cmdline
  bool marg_flag = false;  // set to true if -m        argument is on cmdline
  bool narg_flag = false;  // set to true if -n        argument is on cmdline
  bool  res_flag = false;  // set to true if -r <file> argument is on cmdline
  Real wtlim = 0;

  //--- Step 1. --------------------------------------------------------------------------
  // Initialize environment (must initialize MPI first, then Kokkos)

#if MPI_PARALLEL_ENABLED
#if defined(KOKKOS_ENABLE_HIP)
  // JMF: This is a bizarre workaround to avoid segmentation faults on Frontier.
  // See OLCFDEV-1655: Occasional seg-fault during MPI_Init inside the Frontier
  // documentation.
  (void) hipInit(0);
#endif
#if OPENMP_PARALLEL_ENABLED
  int mpiprv;
  if (MPI_SUCCESS != MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpiprv)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI Initialization failed." << std::endl;
    return(0);
  }
  if (mpiprv != MPI_THREAD_MULTIPLE) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI_THREAD_MULTIPLE must be supported for hybrid parallelzation. "
              << MPI_THREAD_MULTIPLE << " : " << mpiprv
              << std::endl;
    MPI_Finalize();
    return(0);
  }
#else  // no OpenMP
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI Initialization failed." << std::endl;
    return(0);
  }
#endif  // OPENMP_PARALLEL_ENABLED
  // Get process id (rank) in MPI_COMM_WORLD
  if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &(global_variable::my_rank))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI_Comm_rank failed." << std::endl;
    MPI_Finalize();
    return(0);
  }

  // Get total number of MPI processes (ranks)
  if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &global_variable::nranks)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI_Comm_size failed." << std::endl;
    MPI_Finalize();
    return(0);
  }

  if (MPI_SUCCESS != MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,
                                         global_variable::my_rank, MPI_INFO_NULL,
                                         &global_variable::node_comm)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed." << std::endl;
    MPI_Finalize();
    return(0);
  }
  MPI_Comm_rank(global_variable::node_comm, &global_variable::rank_in_node);
  MPI_Comm_size(global_variable::node_comm, &global_variable::ranks_per_node);

  int node_leader_world_rank = global_variable::my_rank;
  MPI_Bcast(&node_leader_world_rank, 1, MPI_INT, 0, global_variable::node_comm);
  std::vector<int> node_leader_for_rank(global_variable::nranks, 0);
  MPI_Allgather(&node_leader_world_rank, 1, MPI_INT, node_leader_for_rank.data(), 1,
                MPI_INT, MPI_COMM_WORLD);

  std::unordered_map<int, int> leader_to_node;
  global_variable::rank_to_node.assign(global_variable::nranks, 0);
  for (int r = 0; r < global_variable::nranks; ++r) {
    auto [it, inserted] = leader_to_node.emplace(node_leader_for_rank[r],
                                                 static_cast<int>(leader_to_node.size()));
    global_variable::rank_to_node[r] = it->second;
  }
  global_variable::node_id = global_variable::rank_to_node[global_variable::my_rank];
  global_variable::nnodes = leader_to_node.size();
#else  // no MPI
  global_variable::my_rank = 0;
  global_variable::nranks  = 1;
  global_variable::node_id = 0;
  global_variable::rank_in_node = 0;
  global_variable::ranks_per_node = 1;
  global_variable::nnodes = 1;
  global_variable::rank_to_node.assign(1, 0);
#endif  // MPI_PARALLEL_ENABLED

  Kokkos::initialize(argc, argv);

  auto FinalizeMpi = []() {
#if MPI_PARALLEL_ENABLED
    if (global_variable::node_comm != MPI_COMM_NULL) {
      MPI_Comm_free(&global_variable::node_comm);
      global_variable::node_comm = MPI_COMM_NULL;
    }
    MPI_Finalize();
#endif
  };

  //--- Step 2. --------------------------------------------------------------------------
  // Check for command line options and respond.

  for (int i=1; i<argc; i++) {
    // If argv[i] is a 2 character string of the form "-?" then:
    if (*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0') {
      // check that command line options that require arguments actually have them:
      char opt_letter = *(argv[i]+1);
      switch(opt_letter) {
        case 'c':
        case 'h':
        case 'm':
        case 'n':
          break;
        default:
          if ((i+1 >= argc) // no argument after option
              || (*argv[i+1] == '-') ) { // option is followed by another option
            if (global_variable::my_rank == 0) {
              std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                        << std::endl << "-" << opt_letter
                        << " must be followed by a valid argument" << std::endl;
              Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
              FinalizeMpi();
#endif
              return(0);
            }
          }
      }

      // set arguments, flags, or execute tasks specified by options
      switch(*(argv[i]+1)) {
        case 'i':                      // -i <input_file>
          input_file.assign(argv[++i]);
          iarg_flag = true;
          break;
        case 'r':                      // -r <restart_file>
          restart_file.assign(argv[++i]);
          res_flag = true;
          break;
        case 'd':                      // -d <run_directory>
          run_dir.assign(argv[++i]);
          break;
        case 'n':
          narg_flag = true;
          break;
        case 'm':
          marg_flag = true;
          break;
        case 't':                      // -t <hh:mm:ss>
          int wth, wtm, wts;
          std::sscanf(argv[++i], "%d:%d:%d", &wth, &wtm, &wts);
          wtlim = static_cast<Real>(wth*3600 + wtm*60 + wts);
          break;
        case 'c':
          if (global_variable::my_rank == 0) ShowConfig();
          Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
          FinalizeMpi();
#endif
          return(0);
          break;
        case 'h':
        default:
          if (global_variable::my_rank == 0) {
            std::cout << "Athena v" << ATHENA_VERSION_MAJOR << "."
                                    << ATHENA_VERSION_MINOR << std::endl;
            std::cout << "Usage: " << argv[0] << " [options] [block/par=value ...]\n";
            std::cout << "Options:" << std::endl;
            std::cout << "  -i <file>       specify input file [athinput]\n";
            std::cout << "  -r <file>       restart with this file\n";
            std::cout << "  -d <directory>  specify run dir [current dir]\n";
            std::cout << "  -n              parse input file and quit\n";
            std::cout << "  -c              show configuration and quit\n";
            std::cout << "  -m              output mesh structure and quit\n";
            std::cout << "  -t hh:mm:ss     wall time limit for final output\n";
            std::cout << "  -h              this help\n";
            ShowConfig();
          }
          Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
          MPI_Finalize();
#endif
          return(0);
          break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }

  // print error if input or restart file not given
  if (restart_file.empty() && input_file.empty()) {
    // no input file is given
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Either an input or restart file must be specified." << std::endl
              << "See " << argv[0] << " -h for options and usage." << std::endl;
    Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
    MPI_Finalize();
#endif
    return(0);
  }

  // Start the wall clock timer. This is done here rather than in the Driver to ensure
  // that the time taken in ProblemGenerator is also captured.
  Kokkos::Timer timer;

  //--- Step 3. --------------------------------------------------------------------------
  // Construct ParameterInput object and load data either from restart or input file.
  // With MPI, the input is read by every rank in parallel using MPI-IO.

  ParameterInput* pinput = new ParameterInput;
  IOWrapper infile, restartfile;
  // read parameters from restart file
  FileShardMode restart_shard_mode = FileShardMode::shared;
  if (res_flag) {
    RestartPathInfo restart_path;
    try {
      restart_path = NormalizeRestartPath(restart_file);
    } catch (const std::runtime_error &err) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << err.what() << std::endl;
      Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
      FinalizeMpi();
#endif
      return(0);
    }

    restart_shard_mode = restart_path.shard_mode;
    restart_base_dir = restart_path.base_dir;
    restart_file_name = restart_path.file_name;
    restart_file = restart_path.normalized_open_path;

    // read parameters from restart file
    restartfile.Open(restart_file.c_str(), IOWrapper::FileMode::read,
                     UsesSerialIO(restart_shard_mode));
    pinput->LoadFromFile(restartfile, restart_shard_mode);
    IOWrapperSizeT headeroffset = restartfile.GetPosition(UsesSerialIO(restart_shard_mode));
  }

  // read parameters from input file.  If both -r and -i are specified, this will
  // override parameters from the restart file
  if (iarg_flag) {
    infile.Open(input_file.c_str(), IOWrapper::FileMode::read);
    pinput->LoadFromFile(infile);
    infile.Close();
  }
  pinput->ModifyFromCmdline(argc, argv);

  // Dump input parameters and quit if code was run with -n option.
  if (narg_flag) {
    if (global_variable::my_rank == 0) pinput->ParameterDump(std::cout);
    if (res_flag) restartfile.Close(UsesSerialIO(restart_shard_mode));
    delete pinput;
    Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
    FinalizeMpi();
#endif
    return(0);
  }

  bool restart_conversion = pinput->DoesBlockExist("restart_conversion") &&
      pinput->DoesParameterExist("restart_conversion", "enabled") &&
      pinput->GetBoolean("restart_conversion", "enabled");
  if (restart_conversion) {
    if (!res_flag || run_dir.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Restart conversion requires both -r and a separate "
                << "output directory supplied with -d." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    source_restart_path = PathToString(fs::absolute(restart_file).lexically_normal());
    source_particle_path = PathToString(fs::absolute(pinput->GetString(
        "restart_conversion", "particle_restart")).lexically_normal());
    // Particle initialization happens before ChangeRunDir(), but store an absolute path
    // so the exact source sidecar is also unambiguous in conversion provenance.
    pinput->SetInteger("particles", "star_prtcl_rst_flag", 1);
    pinput->SetString("particles", "star_particle_file", source_particle_path);
    conversion_output_dir = PathToString(fs::absolute(run_dir).lexically_normal());
    if (fs::exists(conversion_output_dir)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Restart conversion output directory already exists: "
                << conversion_output_dir << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  //--- Step 4. --------------------------------------------------------------------------
  // Construct Mesh.  Then build MeshBlockTree and add MeshBlockPack containing MeshBlocks
  // on this rank.  Latter cannot be performed in Mesh constructor since it requires
  // pointer to Mesh.

  Mesh* pmesh = new Mesh(pinput);
  if (res_flag) {
    pmesh->SetRestartFileInfo(restart_base_dir, restart_file_name, restart_shard_mode);
    pmesh->restart_meta.source_restart_path = source_restart_path;
    pmesh->restart_meta.source_particle_path = source_particle_path;
  }
  if (!res_flag) {
    pmesh->BuildTreeFromScratch(pinput);
  } else {
    pmesh->BuildTreeFromRestart(pinput, restartfile, restart_shard_mode);
  }

  //  If code was run with -m option, write mesh structure to file and quit.
  if (marg_flag) {
    if (global_variable::my_rank == 0) {pmesh->WriteMeshStructure();}
    if (res_flag) {restartfile.Close(UsesSerialIO(restart_shard_mode));}
    delete pmesh;
    delete pinput;
    Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
    FinalizeMpi();
#endif
    return(0);
  }

  //--- Step 5. --------------------------------------------------------------------------
  // Add coordinates and physics modules to MeshBlockPack, and set initial conditions.
  // Note these steps must occur after Mesh (including MeshBlocks and MeshBlockPack)
  // is fully constructed.

  pmesh->AddCoordinatesAndPhysics(pinput);
  if (!res_flag) {
    // set ICs using ProblemGenerator constructor for new runs
    pmesh->pgen = std::make_unique<ProblemGenerator>(pinput, pmesh);
    if (pmesh->pmb_pack->ppart != nullptr) {
      pmesh->pmb_pack->ppart->CreateParticleTags(pinput);
    }
  } else {
    // read ICs from restart file using ProblemGenerator constructor for restarts
    pmesh->pgen = std::make_unique<ProblemGenerator>(pinput,
                                                     pmesh,
                                                     restartfile,
                                                     restart_shard_mode);
    restartfile.Close(UsesSerialIO(restart_shard_mode));
  }
  //--- Step 6. --------------------------------------------------------------------------
  // Construct Driver and Outputs. Actual outputs (including initial conditions) are made
  // in Driver.Initialize(). Add wall clock timer to Driver if necessary.

  ChangeRunDir(run_dir);
  Driver* pdriver = new Driver(pinput, pmesh, wtlim, &timer);
  Outputs* pout = new Outputs(pinput, pmesh);



  //--- Step 7. --------------------------------------------------------------------------
  // Execute Driver.
  //    1. Initial conditions set in Driver::Initialize()
  //    2. TaskList(s) executed in Driver::Execute()
  //    3. Any final analysis or diagnostics run in Driver::Finalize()

  if (restart_conversion) {
    WriteRestartConversion(pdriver, pout, pmesh, pinput, conversion_output_dir);
    if (pmesh->pgen->pgen_final_func != nullptr) {
      (pmesh->pgen->pgen_final_func)(pinput, pmesh);
    }
  } else {
    pdriver->Initialize(pmesh, pinput, pout, res_flag);
    pdriver->Execute(pmesh, pinput, pout);
    pdriver->Finalize(pmesh, pinput, pout);
  }

  //--- Step 8. -------------------------------------------------------------------------
  // clean up, and terminate
  // Note anything containing a Kokkos::view must be deleted before Kokkos::finalize()

  delete pout;
  delete pdriver;
  delete pmesh;
  delete pinput;
  Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
  FinalizeMpi();
#endif
  return(0);
}
