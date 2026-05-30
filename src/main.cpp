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
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <string>
#include <memory>
#include <cstdio> // sscanf
#include <fstream>  // Include this for std::ifstream
#include <limits>
#include <sstream>
#include <vector>

// Athena headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "outputs/outputs.hpp"
#include "driver/driver.hpp"
#include "utils/utils.hpp"

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

constexpr const char *kNodeRestartMagic = "AthenaK node restart manifest version=1";

struct NodeRestartPayload {
  int node;
  int blocks;
  std::uint64_t bytes;
  std::string path;
};

struct NodeRestartSegment {
  int node;
  int gid_start;
  int count;
  int payload_block_start;
};

[[noreturn]] void FailNodeRestart(const std::string &message) {
  std::cerr << "### FATAL ERROR while reading node restart manifest: "
            << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(EXIT_FAILURE);
}

bool IsNodeRestartManifest(const std::string &path) {
  std::ifstream input(path);
  std::string first_line;
  return input.good() && std::getline(input, first_line) &&
      first_line == kNodeRestartMagic;
}

std::string ParentDirectory(const std::string &path) {
  std::size_t slash = path.rfind('/');
  return (slash == std::string::npos) ? std::string(".") : path.substr(0, slash);
}

bool CopyFileRange(std::ifstream &input, std::ofstream &output, std::uint64_t input_offset,
                   std::uint64_t output_offset, std::uint64_t bytes) {
  constexpr std::size_t kCopyBytes = 1024*1024;
  std::vector<char> buffer(kCopyBytes);
  input.clear();
  input.seekg(static_cast<std::streamoff>(input_offset), std::ios::beg);
  output.seekp(static_cast<std::streamoff>(output_offset), std::ios::beg);
  while (bytes > 0) {
    std::size_t amount = static_cast<std::size_t>(
        std::min<std::uint64_t>(bytes, buffer.size()));
    input.read(buffer.data(), static_cast<std::streamsize>(amount));
    if (input.gcount() != static_cast<std::streamsize>(amount)) return false;
    output.write(buffer.data(), static_cast<std::streamsize>(amount));
    if (!output.good()) return false;
    bytes -= amount;
  }
  return true;
}

bool ParseUnsignedField(const std::string &line, const std::string &prefix,
                        std::uint64_t &value) {
  if (line.rfind(prefix, 0) != 0) return false;
  std::string token = line.substr(prefix.size());
  if (token.empty() ||
      !std::all_of(token.begin(), token.end(),
                   [](char ch) { return ch >= '0' && ch <= '9'; })) {
    return false;
  }
  std::istringstream input(token);
  input >> value;
  std::string trailing;
  return input && !(input >> trailing);
}

bool ParseSignedField(const std::string &line, const std::string &prefix, int &value) {
  if (line.rfind(prefix, 0) != 0) return false;
  std::istringstream input(line.substr(prefix.size()));
  input >> value;
  std::string trailing;
  return input && !(input >> trailing);
}

std::string ManifestPayloadPrefix(const std::string &manifest_path) {
  std::size_t slash = manifest_path.rfind('/');
  std::string leaf = (slash == std::string::npos) ? manifest_path :
      manifest_path.substr(slash + 1);
  constexpr const char *suffix = ".rst";
  if (leaf.size() <= 4 || leaf.compare(leaf.size() - 4, 4, suffix) != 0 ||
      leaf.find('/') != std::string::npos || leaf.find('\\') != std::string::npos) {
    FailNodeRestart("manifest filename does not follow the restart leaf contract.");
  }
  return leaf.substr(0, leaf.size() - 4) + ".g";
}

std::string ValidatePayloadPath(const NodeRestartPayload &payload,
                                const std::string &payload_prefix) {
  if (payload.path.empty() || payload.path[0] == '/' ||
      payload.path.find('\\') != std::string::npos) {
    FailNodeRestart("payload path is not a relative node-shard path.");
  }
  std::size_t slash = payload.path.find('/');
  if (slash == std::string::npos || slash == 0 ||
      slash + 1 >= payload.path.size() ||
      payload.path.find('/', slash + 1) != std::string::npos) {
    FailNodeRestart("payload path must contain exactly one node-directory component.");
  }
  std::string directory = payload.path.substr(0, slash);
  std::string leaf = payload.path.substr(slash + 1);
  char expected_directory[32];
  std::snprintf(expected_directory, sizeof(expected_directory), "node_%08d", payload.node);
  if (directory != expected_directory || directory == "." || directory == ".." ||
      leaf == "." || leaf == "..") {
    FailNodeRestart("payload path does not match its declared node directory.");
  }
  constexpr const char *payload_suffix = ".payload.rst";
  if (leaf.rfind(payload_prefix, 0) != 0 ||
      leaf.size() <= payload_prefix.size() + 12 ||
      leaf.compare(leaf.size() - 12, 12, payload_suffix) != 0) {
    FailNodeRestart("payload leaf does not match the generated restart contract.");
  }
  std::string generation = leaf.substr(payload_prefix.size(),
      leaf.size() - payload_prefix.size() - 12);
  if (generation.empty() ||
      !std::all_of(generation.begin(), generation.end(),
                   [](char ch) { return ch >= '0' && ch <= '9'; })) {
    FailNodeRestart("payload leaf has an invalid generation token.");
  }
  return leaf;
}

std::uint64_t ExpectedPayloadBytes(std::uint64_t header_size, std::uint64_t data_size,
                                   int blocks) {
  if (blocks < 0 || (blocks > 0 && data_size >
      (std::numeric_limits<std::uint64_t>::max() - header_size) /
      static_cast<std::uint64_t>(blocks))) {
    FailNodeRestart("payload byte count overflows its restart inventory.");
  }
  return header_size + data_size*static_cast<std::uint64_t>(blocks);
}

std::string StageNodeRestart(const std::string &manifest_path) {
  std::ifstream manifest(manifest_path);
  std::string line;
  if (!std::getline(manifest, line) || line != kNodeRestartMagic) {
    FailNodeRestart("invalid manifest signature in '" + manifest_path + "'.");
  }
  bool complete = false;
  int payload_count = -1;
  int nmb_total = -1;
  std::uint64_t header_size = 0;
  std::uint64_t data_size = 0;
  std::vector<NodeRestartPayload> payloads;
  std::vector<NodeRestartSegment> segments;
  bool saw_complete = false;
  bool saw_payload_count = false;
  bool saw_nmb_total = false;
  bool saw_header_size = false;
  bool saw_data_size = false;
  bool saw_end = false;
  while (std::getline(manifest, line)) {
    if (line == "end") {
      saw_end = true;
      break;
    }
    if (line.rfind("complete=", 0) == 0) {
      if (saw_complete || line != "complete=1") {
        FailNodeRestart("invalid or duplicate completion record.");
      }
      complete = true;
      saw_complete = true;
    } else if (line.rfind("payload_count=", 0) == 0) {
      if (saw_payload_count || !ParseSignedField(line, "payload_count=", payload_count) ||
          payload_count <= 0) {
        FailNodeRestart("invalid or duplicate payload count.");
      }
      saw_payload_count = true;
    } else if (line.rfind("nmb_total=", 0) == 0) {
      if (saw_nmb_total || !ParseSignedField(line, "nmb_total=", nmb_total) ||
          nmb_total < 0) {
        FailNodeRestart("invalid or duplicate total mesh block count.");
      }
      saw_nmb_total = true;
    } else if (line.rfind("header_size=", 0) == 0) {
      if (saw_header_size || !ParseUnsignedField(line, "header_size=", header_size) ||
          header_size == 0) {
        FailNodeRestart("invalid or duplicate header byte count.");
      }
      saw_header_size = true;
    } else if (line.rfind("data_size=", 0) == 0) {
      if (saw_data_size || !ParseUnsignedField(line, "data_size=", data_size)) {
        FailNodeRestart("invalid or duplicate per-block byte count.");
      }
      saw_data_size = true;
    } else if (line.rfind("payload ", 0) == 0) {
      std::istringstream values(line);
      std::string tag;
      NodeRestartPayload payload;
      values >> tag >> payload.node >> payload.blocks >> payload.bytes >> payload.path;
      std::string trailing;
      if (!values || payload.node < 0 || payload.blocks < 0 || (values >> trailing)) {
        FailNodeRestart("malformed payload entry in '" + manifest_path + "'.");
      }
      if (payload.node != static_cast<int>(payloads.size())) {
        FailNodeRestart("payload inventory must be ordered by contiguous node id.");
      }
      payloads.push_back(payload);
    } else if (line.rfind("segment ", 0) == 0) {
      std::istringstream values(line);
      std::string tag;
      NodeRestartSegment segment;
      values >> tag >> segment.node >> segment.gid_start >> segment.count
             >> segment.payload_block_start;
      std::string trailing;
      if (!values || segment.node < 0 || segment.gid_start < 0 || segment.count < 0 ||
          segment.payload_block_start < 0 || (values >> trailing)) {
        FailNodeRestart("malformed segment entry in '" + manifest_path + "'.");
      }
      segments.push_back(segment);
    } else {
      FailNodeRestart("unrecognized inventory record in '" + manifest_path + "'.");
    }
  }
  while (std::getline(manifest, line)) {
    if (!line.empty()) {
      FailNodeRestart("records found after the manifest terminator.");
    }
  }
  if (!complete || !saw_end || !saw_payload_count || !saw_nmb_total ||
      !saw_header_size || !saw_data_size ||
      payload_count != static_cast<int>(payloads.size()) || payloads.empty()) {
    FailNodeRestart("incomplete inventory in '" + manifest_path + "'.");
  }

  std::string expected_payload_leaf;
  std::string payload_prefix = ManifestPayloadPrefix(manifest_path);
  std::vector<int> node_to_payload(static_cast<std::size_t>(payload_count), -1);
  for (std::size_t i = 0; i < payloads.size(); ++i) {
    if (payloads[i].node >= payload_count ||
        node_to_payload[payloads[i].node] != -1) {
      FailNodeRestart("node payload inventory is duplicated or non-contiguous.");
    }
    node_to_payload[payloads[i].node] = static_cast<int>(i);
    std::string leaf = ValidatePayloadPath(payloads[i], payload_prefix);
    if (expected_payload_leaf.empty()) {
      expected_payload_leaf = leaf;
    } else if (leaf != expected_payload_leaf) {
      FailNodeRestart("node payload inventory refers to mixed generations.");
    }
  }
  for (int node = 0; node < payload_count; ++node) {
    if (node_to_payload[node] < 0) {
      FailNodeRestart("node payload inventory omits node " + std::to_string(node) + ".");
    }
  }

  std::vector<int> covered(static_cast<std::size_t>(nmb_total), 0);
  std::vector<int> mapped_blocks(payloads.size(), 0);
  std::vector<int> next_payload_block(payloads.size(), 0);
  int next_gid = 0;
  ExpectedPayloadBytes(header_size, data_size, nmb_total);
  for (const auto &segment : segments) {
    int index = (segment.node < payload_count) ? node_to_payload[segment.node] : -1;
    if (index < 0 || segment.gid_start != next_gid ||
        segment.payload_block_start != next_payload_block[index] ||
        segment.count > nmb_total - segment.gid_start ||
        segment.count > payloads[index].blocks - segment.payload_block_start) {
      FailNodeRestart("segment ordering or node-local range is inconsistent.");
    }
    mapped_blocks[index] += segment.count;
    next_payload_block[index] = segment.payload_block_start + segment.count;
    next_gid = segment.gid_start + segment.count;
    for (int gid = segment.gid_start; gid < segment.gid_start + segment.count; ++gid) {
      if (++covered[gid] != 1) {
        FailNodeRestart("payload segments overlap at mesh block " + std::to_string(gid) + ".");
      }
    }
  }
  if (next_gid != nmb_total) {
    FailNodeRestart("payload segments do not cover the declared block range.");
  }
  for (int gid = 0; gid < nmb_total; ++gid) {
    if (covered[gid] != 1) {
      FailNodeRestart("payload segments do not cover mesh block " + std::to_string(gid) + ".");
    }
  }
  std::string directory = ParentDirectory(manifest_path);
  for (std::size_t i = 0; i < payloads.size(); ++i) {
    if (mapped_blocks[i] != payloads[i].blocks ||
        next_payload_block[i] != payloads[i].blocks ||
        payloads[i].bytes != ExpectedPayloadBytes(header_size, data_size,
                                                   payloads[i].blocks)) {
      FailNodeRestart("payload block count or byte count is inconsistent.");
    }
    std::ifstream payload(directory + "/" + payloads[i].path,
                          std::ios::binary | std::ios::ate);
    std::uint64_t bytes = payload.good()
        ? static_cast<std::uint64_t>(payload.tellg()) : 0;
    if (bytes != payloads[i].bytes) {
      FailNodeRestart("payload '" + payloads[i].path + "' is absent or incomplete.");
    }
  }

  std::string assembled_path = manifest_path + ".assembled";
  std::string temporary_path = assembled_path + ".tmp";
  std::ofstream assembled(temporary_path, std::ios::binary | std::ios::trunc);
  std::ifstream first_payload(directory + "/" + payloads[0].path, std::ios::binary);
  if (!assembled.good() || !first_payload.good() ||
      !CopyFileRange(first_payload, assembled, 0, 0, header_size)) {
    FailNodeRestart("could not stage shared restart header.");
  }
  for (const auto &segment : segments) {
    if (segment.count == 0) continue;
    int index = node_to_payload[segment.node];
    std::ifstream payload(directory + "/" + payloads[index].path, std::ios::binary);
    std::uint64_t bytes = data_size*static_cast<std::uint64_t>(segment.count);
    std::uint64_t source = header_size + data_size*segment.payload_block_start;
    std::uint64_t destination = header_size + data_size*segment.gid_start;
    if (!payload.good() ||
        !CopyFileRange(payload, assembled, source, destination, bytes)) {
      FailNodeRestart("could not stage payload data for node "
                      + std::to_string(segment.node) + ".");
    }
  }
  assembled.close();
  if (!assembled.good() || std::rename(temporary_path.c_str(), assembled_path.c_str()) != 0) {
    FailNodeRestart("could not publish transient assembled restart file.");
  }
  return assembled_path;
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn int main(int argc, char *argv[])
//! \brief Athena main program

int main(int argc, char *argv[]) {
  std::string input_file, restart_file, run_dir;
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
              << "MPI_THREAD_MULTIPLE must be supported for hybrid parallelization. "
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
#else  // no MPI
  global_variable::my_rank = 0;
  global_variable::nranks  = 1;
#endif  // MPI_PARALLEL_ENABLED

  Kokkos::initialize(argc, argv);

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
              MPI_Finalize();
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
          MPI_Finalize();
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
  bool staged_node_restart = false;
  std::string staged_restart_file;
  auto cleanup_staged_restart = [&]() {
    if (!staged_node_restart) return;
#if MPI_PARALLEL_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::remove(staged_restart_file.c_str());
    }
#if MPI_PARALLEL_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  };
  // read parameters from restart file
  bool single_file_per_rank = false; // DBF: flag for single_file_per_rank for rst files
  if (res_flag) {
    if (IsNodeRestartManifest(restart_file)) {
      staged_node_restart = true;
      global_variable::InitializeNodeCommunicator();
      staged_restart_file = restart_file + ".assembled";
      if (global_variable::my_rank == 0) {
        staged_restart_file = StageNodeRestart(restart_file);
      }
#if MPI_PARALLEL_ENABLED
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      restart_file = staged_restart_file;
    }
    // Check if the path contains "rank_" directory
    size_t rank_pos = restart_file.find("/rank_");
    single_file_per_rank = (rank_pos != std::string::npos);

    // If single_file_per_rank is true, modify the path for the current rank
    if (single_file_per_rank) {
        // Extract the base directory and file name
        size_t last_slash = restart_file.rfind('/');
        std::string base_dir = restart_file.substr(0, rank_pos);
        std::string file_name = restart_file.substr(last_slash + 1);

        // Construct the path for the current rank
        char rank_dir[20];
        std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", global_variable::my_rank);
        restart_file = base_dir + "/" + rank_dir + "/" + file_name;
    }

    // Now use restart_file for opening the file
    std::ifstream file_check(restart_file);
    if (!file_check.good()) {
        std::cerr << "Error: Unable to open restart file: " << restart_file << std::endl;
        // Handle the error (e.g., exit the program or use a default configuration)
    }

    // read parameters from restart file
    restartfile.Open(restart_file.c_str(),IOWrapper::FileMode::read,single_file_per_rank);
    pinput->LoadFromFile(restartfile, single_file_per_rank);
    IOWrapperSizeT headeroffset = restartfile.GetPosition(single_file_per_rank);
  }

  // read parameters from input file.  If both -r and -i are specified, this will
  // override parameters from the restart file
  if (iarg_flag) {
    infile.Open(input_file.c_str(), IOWrapper::FileMode::read);
    pinput->LoadFromFile(infile);
    infile.Close();
    pinput->CheckBlockNames();
  }
  pinput->ModifyFromCmdline(argc, argv);

  // Dump input parameters and quit if code was run with -n option.
  if (narg_flag) {
    if (global_variable::my_rank == 0) pinput->ParameterDump(std::cout);
    if (res_flag) restartfile.Close(single_file_per_rank);
    cleanup_staged_restart();
    delete pinput;
    Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
    global_variable::FinalizeNodeCommunicator();
    MPI_Finalize();
#endif
    return(0);
  }

  //--- Step 4. --------------------------------------------------------------------------
  // Construct Mesh.  Then build MeshBlockTree and add MeshBlockPack containing MeshBlocks
  // on this rank.  Latter cannot be performed in Mesh constructor since it requires
  // pointer to Mesh.

  Mesh* pmesh = new Mesh(pinput);
  if (!res_flag) {
    pmesh->BuildTreeFromScratch(pinput);
  } else {
    pmesh->BuildTreeFromRestart(pinput, restartfile, single_file_per_rank);
  }

  //  If code was run with -m option, write mesh structure to file and quit.
  if (marg_flag) {
    if (global_variable::my_rank == 0) {pmesh->WriteMeshStructure();}
    if (res_flag) {restartfile.Close(single_file_per_rank);}
    cleanup_staged_restart();
    delete pmesh;
    delete pinput;
    Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
    global_variable::FinalizeNodeCommunicator();
    MPI_Finalize();
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
  } else {
    // read ICs from restart file using ProblemGenerator constructor for restarts
    pmesh->pgen = std::make_unique<ProblemGenerator>(pinput,
                                                     pmesh,
                                                     restartfile,
                                                     single_file_per_rank);
    restartfile.Close(single_file_per_rank);
    cleanup_staged_restart();
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

  pdriver->Initialize(pmesh, pinput, pout, res_flag);
  pdriver->Execute(pmesh, pinput, pout);
  pdriver->Finalize(pmesh, pinput, pout);

  //--- Step 8. -------------------------------------------------------------------------
  // clean up, and terminate
  // Note anything containing a Kokkos::view must be deleted before Kokkos::finalize()

  delete pout;
  delete pdriver;
  delete pmesh;
  delete pinput;
  Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
  global_variable::FinalizeNodeCommunicator();
  MPI_Finalize();
#endif
  return(0);
}
