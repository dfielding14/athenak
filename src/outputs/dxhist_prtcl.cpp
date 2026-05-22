//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dxhist_prtcl.cpp
//! \brief writes displacement histograms for cosmic-ray particles

#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

ParticleDxHistOutput::ParticleDxHistOutput(ParameterInput *pin, Mesh *pm,
                                           OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("dxh",0775);
  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  vmin = pin->GetOrAddReal(op.block_name,"vmin",-1.0);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",1.0);
  Kokkos::realloc(host_histogram, nspec*3, nbin);
  dxhist_single_file_per_rank = pin->GetOrAddInteger(op.block_name,
                                                     "dxhist_single_file_per_rank",0);
  reduce_histogram = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (dxhist_single_file_per_rank) {reduce_histogram = false;}
  if (dxhist_single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "dxh/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleDxHistOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  Real vmin_ = vmin;
  Real vmax_ = vmax;
  int nbin_ = nbin;

  DvceArray2D<int> local_histogram("local_dxh_hist", nspecies*3, nbin);
  Kokkos::deep_copy(local_histogram, 0);
  par_for("part_dxh_hist",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    Real dx[3] = {pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p)};
    for (int d=0; d<3; ++d) {
      int ip = static_cast<int>((dx[d] - vmin_)/(vmax_ - vmin_)*nbin_);
      ip = (ip < 0) ? 0 : ip;
      ip = (ip > nbin_ - 1) ? nbin_ - 1 : ip;
      Kokkos::atomic_add(&local_histogram(3*spec + d,ip),1);
    }
  });
  Kokkos::deep_copy(host_histogram, local_histogram);
#if MPI_PARALLEL_ENABLED
  if (reduce_histogram) {
    HostArray2D<int> global_histogram("global_dxh_hist", nspecies*3, nbin);
    MPI_Reduce(host_histogram.data(), global_histogram.data(), nspecies*3*nbin,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_histogram, global_histogram);
    }
  }
#endif
}

void ParticleDxHistOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname("dxh/");
  if (dxhist_single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + ".dxh";

  if (reduce_histogram && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  bool single = static_cast<bool>(dxhist_single_file_per_rank) || reduce_histogram;
  if (dxhist_single_file_per_rank || reduce_histogram ||
      global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << "\n# AthenaK particle displacement histogram at time= " << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle << std::endl;
    FILE *pfile = std::fopen(fname.c_str(),"a");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,"%s\n",msg.str().c_str());
    std::fclose(pfile);
  }

  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append, single);
  std::size_t header_offset = partfile.GetPosition(single);

  int nspecies = pm->pmb_pack->ppart->nspecies;
  int nout = 3*nspecies*nbin;
  int *data = new int[std::max(1,nout)];
  for (int sp=0; sp<nspecies; ++sp) {
    for (int bin=0; bin<nbin; ++bin) {
      data[sp*3*nbin + bin] = host_histogram(3*sp,bin);
      data[sp*3*nbin + nbin + bin] = host_histogram(3*sp + 1,bin);
      data[sp*3*nbin + 2*nbin + bin] = host_histogram(3*sp + 2,bin);
    }
  }

  std::size_t myoffset = header_offset;
  if (reduce_histogram) {
    if (partfile.Write_any_type(data,nout,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle displacement histogram output failed"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (!single) {
      myoffset += static_cast<std::size_t>(global_variable::my_rank)*nout*sizeof(int);
    }
    if (partfile.Write_any_type_at_all(data,nout,myoffset,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle displacement histogram output failed"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  partfile.Close(single);
  delete[] data;

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
