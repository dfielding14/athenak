//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file track_prtcl.cpp
//! \brief writes data for tracked particles in unformatted binary

#include <vector>

#include <algorithm>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "mpi_utils.hpp"
#include "mesh/mesh.hpp"
#include "output_file_utils.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

namespace {

[[noreturn]] void FatalTrackedParticleError(const std::string &message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" +
                        message);
}

IOWrapperSizeT CheckedTrackedParticleProduct(IOWrapperSizeT left,
                                             IOWrapperSizeT right,
                                             const std::string &context) {
  if (left != 0 && right > std::numeric_limits<IOWrapperSizeT>::max()/left) {
    FatalTrackedParticleError(context + " offset overflow.");
  }
  return left*right;
}

IOWrapperSizeT CheckedTrackedParticleAdd(IOWrapperSizeT left,
                                         IOWrapperSizeT right,
                                         const std::string &context) {
  if (right > std::numeric_limits<IOWrapperSizeT>::max() - left) {
    FatalTrackedParticleError(context + " offset overflow.");
  }
  return left + right;
}

}  // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  if (pm->pmb_pack->ppart == nullptr) {
    FatalTrackedParticleError(
        "Tracked-particle output requires a <particles> block.");
  }
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  if (ntrack < 0 || ntrack > pm->nprtcl_total) {
    FatalTrackedParticleError(
        "Tracked-particle output requires nparticles between 0 and the total "
        "particle count (" + std::to_string(pm->nprtcl_total) + "), but received " +
        std::to_string(ntrack) + ".");
  }

  // create new directory for this output. Comments in binary.cpp constructor explain why
  output_file_utils::EnsureDirectory("trk", 0775, "tracked particle output",
                                     FatalTrackedParticleError);
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  ntrack_thisrank = ntrack;
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  if (pm->pmb_pack->ppart == nullptr) {
    FatalTrackedParticleError(
        "Tracked-particle output requires a <particles> block.");
  }

  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank);
  DvceArray1D<int> counter("tracked_particle_counter", 1);
  Kokkos::deep_copy(counter, 0);
  int npart = pm->nprtcl_thisrank;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int tag = pi(PTAG,p);
    if (tag >= 0 && tag < ntrack) {
      int index = Kokkos::atomic_fetch_add(&counter(0),1);
      if (index < ntrack) {
        tracked_prtcl.d_view(index).tag = tag;
        tracked_prtcl.d_view(index).x   = pr(IPX,p);
        tracked_prtcl.d_view(index).y   = pr(IPY,p);
        tracked_prtcl.d_view(index).z   = pr(IPZ,p);
        tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
        tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
        tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
      }
    }
  });
  Kokkos::fence();
  auto host_counter = Kokkos::create_mirror_view(counter);
  Kokkos::deep_copy(host_counter, counter);
  npout = host_counter(0);
  if (npout > ntrack) {
    FatalTrackedParticleError(
        "Tracked-particle tag selection exceeded the requested particle count.");
  }

  // share number of tracked particles to be output across all ranks
  npout_eachrank[global_variable::my_rank] = npout;
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(
      MPI_Allgather(&npout, 1, MPI_INT, npout_eachrank.data(), 1, MPI_INT,
                    MPI_COMM_WORLD),
      "MPI_Allgather for tracked-particle output counts");
#endif

  std::uint64_t global_count = 0;
  for (int count : npout_eachrank) {
    global_count += static_cast<std::uint64_t>(count);
  }
  if (global_count != static_cast<std::uint64_t>(ntrack)) {
    FatalTrackedParticleError(
        "Tracked-particle tags below nparticles must form the dense range "
        "[0, nparticles); selected " + std::to_string(global_count) +
        " records for nparticles=" + std::to_string(ntrack) + ".");
  }

  // sync tracked particle device array with host
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

  // copy host view into host outpart array
  Kokkos::realloc(outpart, npout);
  auto tracked_host = Kokkos::subview(
      tracked_prtcl.h_view, std::make_pair(0, npout));
  Kokkos::deep_copy(outpart, tracked_host);

  // The fixed-width file layout requires exactly one owner for every tag slot.
  std::vector<int> tag_owners(ntrack, 0);
  for (int p=0; p<npout; ++p) {
    int tag = outpart(p).tag;
    if (tag < 0 || tag >= ntrack) {
      FatalTrackedParticleError(
          "Tracked-particle tag is outside the requested dense range.");
    }
    ++tag_owners[tag];
  }
#if MPI_PARALLEL_ENABLED
  if (ntrack > 0) {
    std::vector<int> global_tag_owners(ntrack, 0);
    mpi_utils::CheckMpi(
        MPI_Allreduce(tag_owners.data(), global_tag_owners.data(), ntrack,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD),
        "MPI_Allreduce for tracked-particle tag ownership");
    tag_owners.swap(global_tag_owners);
  }
#endif
  for (int tag=0; tag<ntrack; ++tag) {
    if (tag_owners[tag] != 1) {
      FatalTrackedParticleError(
          "Tracked-particle tags must form the dense range [0, nparticles); "
          "tag " + std::to_string(tag) + " has " +
          std::to_string(tag_owners[tag]) + " owners.");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // create filename: "trk/file_basename".trk
  std::string fname;
  fname.assign("trk/");
  fname.append(out_params.file_basename);
  fname.append(".trk");

  // Root process opens/creates file and appends string
  if (global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle
        << "  ntracked_prtcls=" << ntrack << std::endl;
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
      FatalTrackedParticleError("Output file '" + fname + "' could not be opened.");
    }
    if (std::fprintf(pfile,"%s \n",msg.str().c_str()) < 0) {
      FatalTrackedParticleError("Tracked-particle header could not be written to '" +
                                fname + "'.");
    }
    if (std::fclose(pfile) != 0) {
      FatalTrackedParticleError("Tracked-particle header could not be closed in '" +
                                fname + "'.");
    }
  }
#if MPI_PARALLEL_ENABLED
  mpi_utils::CheckMpi(
      MPI_Barrier(MPI_COMM_WORLD),
      "MPI_Barrier after tracked-particle header write");
#endif

  // Now all ranks open file and append data
  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
  IOWrapperSizeT header_offset = partfile.GetPosition();

  // Payload floats intentionally retain the native-endian legacy format.
  constexpr IOWrapperSizeT values_per_particle = 6;
  std::size_t staging_values = output_file_utils::CheckedSizeProduct(
      static_cast<std::size_t>(values_per_particle), static_cast<std::size_t>(npout),
      "Tracked-particle staging", FatalTrackedParticleError);
  std::vector<float> data(staging_values);
  // Loop over particles, load positions into data[]
  for (int p=0; p<npout; ++p) {
    std::size_t begin = static_cast<std::size_t>(values_per_particle)*p;
    data[begin    ] = static_cast<float>(outpart(p).x);
    data[begin + 1] = static_cast<float>(outpart(p).y);
    data[begin + 2] = static_cast<float>(outpart(p).z);
    data[begin + 3] = static_cast<float>(outpart(p).vx);
    data[begin + 4] = static_cast<float>(outpart(p).vy);
    data[begin + 5] = static_cast<float>(outpart(p).vz);
  }
  int npout_min = *std::min_element(npout_eachrank.begin(), npout_eachrank.end());
  IOWrapperSizeT record_bytes = CheckedTrackedParticleProduct(
      values_per_particle, sizeof(float), "Tracked-particle record");

  // Write tracked particle data collectively over minimum shared number of prtcls
  for (int p=0; p<npout_min; ++p) {
    IOWrapperSizeT myoffset = CheckedTrackedParticleAdd(
        header_offset,
        CheckedTrackedParticleProduct(
            record_bytes, static_cast<IOWrapperSizeT>(outpart(p).tag),
            "Tracked-particle record"),
        "Tracked-particle record");
    std::size_t begin = static_cast<std::size_t>(values_per_particle)*p;
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at_all(data.data() + begin, values_per_particle,
                                       myoffset, "float") != values_per_particle) {
      FatalTrackedParticleError(
          "Particle data not written correctly to tracked-particle file.");
    }
  }
  // Write particle positions individually for remaining particles on each rank
  for (int p=npout_min; p<npout; ++p) {
    IOWrapperSizeT myoffset = CheckedTrackedParticleAdd(
        header_offset,
        CheckedTrackedParticleProduct(
            record_bytes, static_cast<IOWrapperSizeT>(outpart(p).tag),
            "Tracked-particle record"),
        "Tracked-particle record");
    std::size_t begin = static_cast<std::size_t>(values_per_particle)*p;
    // Write particle positions individually for local remainder records.
    if (partfile.Write_any_type_at(data.data() + begin, values_per_particle,
                                   myoffset, "float") != values_per_particle) {
      FatalTrackedParticleError(
          "Particle data not written correctly to tracked-particle file.");
    }
  }

  // close the output file and clean up
  if (partfile.Close() != 0) {
    FatalTrackedParticleError("Tracked-particle output file could not be closed.");
  }

  // increment counters
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}
