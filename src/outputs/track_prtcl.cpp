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
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("trk",0775);
  particles_per_species = pin->GetOrAddInteger(op.block_name,"particles_per_species",0);
  if (particles_per_species < 0) {
    std::cout << "### FATAL ERROR: trk particles_per_species must be nonnegative"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (particles_per_species > 0) {
    auto pp = pm->pmb_pack->ppart;
    if (pp == nullptr || !pm->pmb_pack->frozen_mhd ||
        pp->pusher != ParticlesPusher::boris ||
        pin->GetOrAddString("particles","assign_tag","index_order") != "index_order" ||
        particles_per_species > pm->nprtcl_total/pp->nspecies ||
        pm->nprtcl_total > (1 << 24)) {
      std::cout << "### FATAL ERROR: rich trk requires frozen Boris particles, index_order, "
                << "particles_per_species <= available particles per bin, and "
                << "at most 2^24 particles for exact float32 IDs" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ntrack = particles_per_species*pp->nspecies;
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "trk/rank_%08d", global_variable::my_rank);
    mkdir(rank_dir,0775);
    std::string fname = std::string(rank_dir) + "/" + op.file_basename + ".trk";
    track_file = std::fopen(fname.c_str(), "ab");
    track_buffer.resize(1024*1024);
    if (track_file == nullptr ||
        std::setvbuf(track_file, track_buffer.data(), _IOFBF, track_buffer.size()) != 0) {
      std::cout << "### FATAL ERROR: cannot open buffered track output " << fname
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    return;
  }
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  if (pm->pmb_pack->ppart == nullptr || ntrack < 1 || ntrack > pm->nprtcl_total) {
    std::cout << "### FATAL ERROR: trk requires 1 <= nparticles <= particle count"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack;
}

TrackedParticleOutput::~TrackedParticleOutput() {
  if (track_file != nullptr && std::fclose(track_file) != 0) {
    std::cout << "### FATAL ERROR: failed to flush rich particle tracks" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  if (particles_per_species > 0) {
    pm->pmb_pack->ppart->LoadTrackedData(particles_per_species, track_records, npout);
    return;
  }
  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank);
  int npart = pm->pmb_pack->ppart->nprtcl_thispack;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  int ntrack_local = ntrack;
  DvceArray1D<int> counter("tracked_particle_counter", 1);
  Kokkos::deep_copy(counter, 0);
  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    if (pi(PTAG,p) >= 0 && pi(PTAG,p) < ntrack_local) {
      int index = Kokkos::atomic_fetch_add(&counter(0),1);
      if (index >= ntrack_local) {return;}
      tracked_prtcl.d_view(index).tag = pi(PTAG,p);
      tracked_prtcl.d_view(index).x   = pr(IPX,p);
      tracked_prtcl.d_view(index).y   = pr(IPY,p);
      tracked_prtcl.d_view(index).z   = pr(IPZ,p);
      tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
      tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
      tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
    }
  });
  auto host_counter = Kokkos::create_mirror_view_and_copy(HostMemSpace(), counter);
  npout = host_counter(0);
  if (npout > ntrack) {
    std::cout << "### FATAL ERROR: duplicate tracked particle tags" << std::endl;
    std::exit(EXIT_FAILURE);
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
  Kokkos::deep_copy(outpart, Kokkos::subview(tracked_prtcl.h_view,
                                           std::make_pair(0, npout)));
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (particles_per_species > 0) {
    if (last_output_cycle == pm->ncycle) {return;}
    std::vector<float> records(static_cast<std::size_t>(npout)*17);
    for (int p=0; p<npout; ++p) {
      for (int field=0; field<17; ++field) {
        records[17*p + field] = static_cast<float>(track_records.h_view(p,field));
      }
    }
    int header = std::fprintf(track_file,
        "\n# AthenaK tracked particle data at time= %.17g  nranks= %d  cycle=%d"
        "  ntracked_prtcls=%d  ntrack_per_species=%d  record_count=%d\n"
        "# trk_format=rich_v1  nfields=17  layout=per_rank  rank=%d\n"
        "# fields=tag,track_tag,species,mass,x,y,z,vx,vy,vz,bx,by,bz,dx,dy,dz,dpar\n\n",
        static_cast<double>(pm->time), global_variable::nranks, pm->ncycle,
        ntrack, particles_per_species, npout, global_variable::my_rank);
    if (header < 0 || std::fwrite(records.data(), sizeof(float), records.size(),
                                track_file) != records.size()) {
      std::cout << "### FATAL ERROR: failed to write rich particle tracks" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    last_output_cycle = pm->ncycle;
    out_params.last_time = (out_params.last_time < 0.0) ? pm->time :
                          out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

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
//  std::size_t header_offset = 0;

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
  // calculate local data offset
  int npout_min = npout_eachrank[0];
  for (int n=1; n<global_variable::nranks; ++n) {
    npout_min = std::min(npout_min, npout_eachrank[n]);
  }

  // Write tracked particle data collectively over minimum shared number of prtcls
  for (int p=0; p<npout_min; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset = header_offset +
                           6*sizeof(float)*static_cast<std::size_t>(outpart(p).tag);
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
    std::size_t myoffset = header_offset +
                           6*sizeof(float)*static_cast<std::size_t>(outpart(p).tag);
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
