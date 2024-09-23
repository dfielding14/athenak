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

ParticleDxHistOutput::ParticleDxHistOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;	  
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("dxh",0775);
  // allocate arrays
  nbin = pin->GetOrAddInteger(op.block_name,"nbin", 100);
  vmin = pin->GetOrAddReal(op.block_name, "vmin", -1.0);
  vmax = pin->GetOrAddReal(op.block_name, "vmax", 1.0);
  Kokkos::realloc(host_histogram, nspec*3, nbin);  
  dxhist_single_file_per_rank = pin->GetOrAddInteger(op.block_name,"dxhist_single_file_per_rank", 0);
  if (dxhist_single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "dxh/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir, 0775);
  }  
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void ParticleDxHistOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;

  // Load data for tracked particles on this rank into new device array
  int npart = pm->nprtcl_thisrank;
  npart = pm->pmb_pack->ppart->nprtcl_thispack;  
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  auto nspecies_ = pp->nspecies;
  auto vmin_ = vmin;
  auto vmax_ = vmax;
  auto nbin_ = nbin;
  auto nprtcl_total_ = pm->nprtcl_total;
  Kokkos::View<int**> local_histogram("local_hist", nspecies_*3, nbin_);
  Kokkos::deep_copy(local_histogram, 0);  
  par_for("part_hist",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
      int spec = pi(PSP,p);
      int ip = (pr(IPDX,p) - vmin_)/(vmax_-vmin_)*nbin_;  	
      ip = max(0,ip);
      ip = min(nbin_-1, ip);
      Kokkos::atomic_add(&local_histogram(spec*3, ip),1);
      ip = (pr(IPDY,p) - vmin_)/(vmax_-vmin_)*nbin_;
      ip = max(0,ip);
      ip = min(nbin_-1, ip);      
      Kokkos::atomic_add(&local_histogram(spec*3+1, ip),1);
      ip = (pr(IPDZ,p) - vmin_)/(vmax_-vmin_)*nbin_;
      ip = max(0,ip);
      ip = min(nbin_-1, ip);      
      Kokkos::atomic_add(&local_histogram(spec*3+2, ip),1);      
    });
  // Wait for all threads to finish
  Kokkos::fence();
  auto local_histogram_m = Kokkos::create_mirror_view(local_histogram);  
  Kokkos::deep_copy(local_histogram_m, local_histogram);
  Kokkos::deep_copy(host_histogram, local_histogram_m);
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void ParticleDxHistOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  int big_end = IsBigEndian(); // =1 on big endian machine

  // create filename: "trk/file_basename".trk
  std::string fname;
  fname.assign("dxh/");
  fname.append(out_params.file_basename);
  fname.append(".dxh");

  if(dxhist_single_file_per_rank){
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("dxh/") + std::string(rank_dir) + out_params.file_basename  + ".dxh";
      std::stringstream msg;
      msg << std::endl << "# AthenaK particle spatial displacement histogram at time= " << pm->time
        << "  rank= " << global_variable::my_rank
        << "  cycle=" << pm->ncycle <<std::endl;
      FILE *pfile;
      if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
        exit(EXIT_FAILURE);
      }
      std::fprintf(pfile,"%s \n",msg.str().c_str());
      std::fclose(pfile);  
  }

  else{
    // Root process opens/creates file and appends string
    if (global_variable::my_rank == 0) {
      std::stringstream msg;
      msg << std::endl << "# AthenaK particle spatial displacement histogram at time=" << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle <<std::endl;
      FILE *pfile;
      if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
        exit(EXIT_FAILURE);
      }
      std::fprintf(pfile,"%s \n",msg.str().c_str());
      std::fclose(pfile);
    }
  }
  // Now all ranks open file and append data
  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
  std::size_t header_offset = partfile.GetPosition();
//  std::size_t header_offset = 0;
  int nspecies_ = pm->pmb_pack->ppart->nspecies;
  // allocate 1D vector of floats used to convert and output particle data
  int nout = 3*nspecies_*nbin;
  int *data = new int[nout];
  // Loop over particles, load positions into data[]
  for (int sp=0; sp<nspecies_; ++sp) { 
    for (int bin=0; bin<nbin; ++bin){
      data[sp*3*nbin + bin] = host_histogram(sp*3, bin);
      data[sp*3*nbin + nbin +  bin] = host_histogram(sp*3+1, bin);
      data[sp*3*nbin + 2*nbin + bin] = host_histogram(sp*3+2, bin);
    }
  }
  // calculate local data offset
  std::vector<int> rank_offset(global_variable::nranks, 0);
  for (int n=1; n<global_variable::nranks; ++n) {
    rank_offset[n] = rank_offset[n-1] + nout;
  }
  
  std::size_t datasize = sizeof(int);
  std::size_t myoffset = header_offset;
  if (!dxhist_single_file_per_rank) myoffset +=  rank_offset[global_variable::my_rank] * datasize;    
  // Write particle positions collectively for minimum number of particles across ranks
  if (partfile.Write_any_type_at_all(&(data[0]),nout,myoffset,"int") != nout) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "particle data not written correctly to tracked particle file"
        << std::endl;
    exit(EXIT_FAILURE);
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
