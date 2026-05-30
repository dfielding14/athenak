//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file track_prtcl.cpp
//! \brief writes data for tracked particles in unformatted binary

#include <sys/stat.h>  // mkdir
#include <vector>

#include <algorithm>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  if (pm->pmb_pack->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Tracked-particle output requires an active <particles> block"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("trk",0775);
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  if (ntrack < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked-particle output nparticles must be nonnegative"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack;
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank);
  int npart = pm->nprtcl_thisrank;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  int ntrack_ = ntrack;
  int ntrack_thisrank_ = ntrack_thisrank;

  // Create device-side counters
  Kokkos::View<int> d_counter("counter");
  Kokkos::View<int> d_overflow("overflow");
  Kokkos::deep_copy(d_counter, 0);
  Kokkos::deep_copy(d_overflow, 0);

  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    const int tag = pi(PTAG,p);
    if ((tag >= 0) && (tag < ntrack_)) {
      int index = Kokkos::atomic_fetch_add(&d_counter(),1);
      if (index < ntrack_thisrank_) {
        tracked_prtcl.d_view(index).tag = tag;
        tracked_prtcl.d_view(index).x   = pr(IPX,p);
        tracked_prtcl.d_view(index).y   = pr(IPY,p);
        tracked_prtcl.d_view(index).z   = pr(IPZ,p);
        tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
        tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
        tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
      } else {
        Kokkos::atomic_fetch_add(&d_overflow(),1);
      }
    }
  });

  Kokkos::fence();
  Kokkos::deep_copy(npout, d_counter);
  int overflow = 0;
  Kokkos::deep_copy(overflow, d_overflow);
  if ((overflow > 0) || (npout > ntrack_thisrank)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "Tracked-particle output found " << npout
      << " particles with tags in [0," << ntrack << "), but this rank only allocated "
      << ntrack_thisrank << " output slots" << std::endl;
    exit(EXIT_FAILURE);
  }

  // share number of tracked particles to be output across all ranks
  npout_eachrank[global_variable::my_rank] = npout;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&npout, 1, MPI_INT, npout_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#endif

  // sync tracked particle device array with host
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

  // copy host view into host outpart array
  Kokkos::realloc(outpart, npout);
  if (npout > 0) {
    auto tracked_slice = Kokkos::subview(tracked_prtcl.h_view, std::make_pair(0, npout));
    Kokkos::deep_copy(outpart, tracked_slice);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  int big_end = IsBigEndian(); // =1 on big endian machine

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
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,"%s \n",msg.str().c_str());
    std::fclose(pfile);
  }

  // Now all ranks open file and append data
  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
  std::size_t header_offset = partfile.GetPosition();

  // allocate 1D vector of floats used to convert and output particle data
  float *data = new float[6*npout];
  // Loop over particles, load positions into data[]
  for (int p=0; p<npout; ++p) {
    data[ 6*p   ] = static_cast<float>(outpart(p).x);
    data[(6*p)+1] = static_cast<float>(outpart(p).y);
    data[(6*p)+2] = static_cast<float>(outpart(p).z);
    data[(6*p)+3] = static_cast<float>(outpart(p).vx);
    data[(6*p)+4] = static_cast<float>(outpart(p).vy);
    data[(6*p)+5] = static_cast<float>(outpart(p).vz);
  }
  if (!big_end) {
    for (int i=0; i<(6*npout); ++i) { Swap4Bytes(&data[i]); }
  }

  // calculate local data offset
  std::vector<int> rank_offset(global_variable::nranks, 0);
  int npout_min = npout_eachrank[0];
  for (int n=1; n<global_variable::nranks; ++n) {
    rank_offset[n] = rank_offset[n-1] + npout_eachrank[n-1];
    npout_min = std::min(npout_min, npout_eachrank[n]);
  }

  std::vector<int> tag_counts(ntrack, 0);
  for (int p=0; p<npout; ++p) {
    const int tag = outpart(p).tag;
    if (tag < 0 || tag >= ntrack) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Tracked-particle output encountered invalid tag "
                << tag << " for ntracked_prtcls=" << ntrack << std::endl;
      std::exit(EXIT_FAILURE);
    }
    tag_counts[tag] += 1;
  }
#if MPI_PARALLEL_ENABLED
  if (ntrack > 0) {
    std::vector<int> global_tag_counts(ntrack, 0);
    MPI_Allreduce(tag_counts.data(), global_tag_counts.data(), ntrack, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    tag_counts.swap(global_tag_counts);
  }
#endif
  for (int tag=0; tag<ntrack; ++tag) {
    if (tag_counts[tag] != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Tracked-particle output expected exactly one particle with tag "
                << tag << ", but found " << tag_counts[tag] << "." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Write tracked particle data collectively over minimum shared number of prtcls
  for (int p=0; p<npout_min; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset = header_offset + 6 * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at_all(&(data[6*p]),6,myoffset,"float") != 6) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Write particle positions individually for remaining particles on each rank
  for (int p=npout_min; p<npout; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset = header_offset + 6 * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at(&(data[6*p]),6,myoffset,"float") != 6) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // close the output file and clean up
  partfile.Close();
  delete[] data;

  // increment counters
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}
