//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file df_prtcl.cpp
//! \brief writes pitch-angle histograms for cosmic-ray particles

#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

ParticleDFOutput::ParticleDFOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("df",0775);
  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  vmin = pin->GetOrAddReal(op.block_name,"vmin",-1.0);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",1.0);
  Kokkos::realloc(host_histogram, nspec, nbin);
  df_single_file_per_rank = pin->GetOrAddInteger(op.block_name,
                                                 "df_single_file_per_rank",0);
  reduce_histogram = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (df_single_file_per_rank) {reduce_histogram = false;}
  if (df_single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "df/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleDFOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  Real vmin_ = vmin;
  Real vmax_ = vmax;
  int nbin_ = nbin;

  DvceArray2D<int> local_histogram("local_df_hist", nspecies, nbin);
  Kokkos::deep_copy(local_histogram, 0);
  par_for("part_df_hist",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    Real bmag = sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real vmag = sqrt(SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p)));
    Real mu = 0.0;
    if (bmag > 0.0 && vmag > 0.0) {
      mu = (pr(IPVX,p)*pr(IPBX,p) + pr(IPVY,p)*pr(IPBY,p) +
            pr(IPVZ,p)*pr(IPBZ,p))/(vmag*bmag);
    }
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    int ip = static_cast<int>((mu - vmin_)/(vmax_ - vmin_)*nbin_);
    ip = (ip < 0) ? 0 : ip;
    ip = (ip > nbin_ - 1) ? nbin_ - 1 : ip;
    Kokkos::atomic_add(&local_histogram(spec,ip),1);
  });
  Kokkos::deep_copy(host_histogram, local_histogram);
#if MPI_PARALLEL_ENABLED
  if (reduce_histogram) {
    HostArray2D<int> global_histogram("global_df_hist", nspecies, nbin);
    MPI_Reduce(host_histogram.data(), global_histogram.data(), nspecies*nbin,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_histogram, global_histogram);
    }
  }
#endif
  pp->LogPerformance("df_output", 0, 0, 0, 0, npart);
}

void ParticleDFOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname("df/");
  if (df_single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + ".df";

  if (reduce_histogram && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  bool single = static_cast<bool>(df_single_file_per_rank) || reduce_histogram;
  if (df_single_file_per_rank || reduce_histogram ||
      global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << "\n# AthenaK particle distribution function at time= " << pm->time
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
  int nout = nspecies*nbin;
  int *data = new int[std::max(1,nout)];
  for (int sp=0; sp<nspecies; ++sp) {
    for (int bin=0; bin<nbin; ++bin) {
      data[sp*nbin + bin] = host_histogram(sp,bin);
    }
  }

  std::size_t myoffset = header_offset;
  if (reduce_histogram) {
    if (partfile.Write_any_type(data,nout,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle DF output failed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (!single) {
      myoffset += static_cast<std::size_t>(global_variable::my_rank)*nout*sizeof(int);
    }
    if (partfile.Write_any_type_at_all(data,nout,myoffset,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle DF output failed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  partfile.Close(single);
  delete[] data;

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}

namespace {

constexpr int kParticleMomentCount = 12;

enum ParticleMomentIndex {
  PMOM_COUNT=0,
  PMOM_MU=1,
  PMOM_MU2=2,
  PMOM_DX=3,
  PMOM_DY=4,
  PMOM_DZ=5,
  PMOM_DX2=6,
  PMOM_DY2=7,
  PMOM_DZ2=8,
  PMOM_DPAR=9,
  PMOM_DPAR2=10,
  PMOM_SPEED2=11
};

} // namespace

ParticleMomentsOutput::ParticleMomentsOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op),
  nmoment(kParticleMomentCount) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("pmom",0775);
  Kokkos::realloc(host_moments, nspec, nmoment);
  single_file_per_rank = pin->GetOrAddInteger(op.block_name,
                                              "pmom_single_file_per_rank",0);
  reduce_moments = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (single_file_per_rank) {reduce_moments = false;}
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "pmom/rank_%08d/",
                  global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleMomentsOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;

  DvceArray2D<Real> local_moments("local_particle_moments", nspecies, nmoment);
  Kokkos::deep_copy(local_moments, 0.0);
  par_for("part_moments",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    Real bmag = sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real vmag = sqrt(speed2);
    Real mu = 0.0;
    if (bmag > 0.0 && vmag > 0.0) {
      mu = (pr(IPVX,p)*pr(IPBX,p) + pr(IPVY,p)*pr(IPBY,p) +
            pr(IPVZ,p)*pr(IPBZ,p))/(vmag*bmag);
    }

    Kokkos::atomic_add(&local_moments(spec,PMOM_COUNT),1.0);
    Kokkos::atomic_add(&local_moments(spec,PMOM_MU),mu);
    Kokkos::atomic_add(&local_moments(spec,PMOM_MU2),mu*mu);
    Kokkos::atomic_add(&local_moments(spec,PMOM_DX),pr(IPDX,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DY),pr(IPDY,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DZ),pr(IPDZ,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DX2),SQR(pr(IPDX,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DY2),SQR(pr(IPDY,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DZ2),SQR(pr(IPDZ,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DPAR),pr(IPDB,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DPAR2),SQR(pr(IPDB,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_SPEED2),speed2);
  });
  Kokkos::deep_copy(host_moments, local_moments);
#if MPI_PARALLEL_ENABLED
  if (reduce_moments) {
    HostArray2D<Real> global_moments("global_particle_moments", nspecies, nmoment);
    MPI_Reduce(host_moments.data(), global_moments.data(), nspecies*nmoment,
               MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_moments, global_moments);
    }
  }
#endif
  pp->LogPerformance("pmom_output", 0, 0, 0, 0, npart);
}

void ParticleMomentsOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname("pmom/");
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + ".pmom";

  if (reduce_moments && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  FILE *pfile = std::fopen(fname.c_str(),"a");
  if (pfile == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Output file '" << fname << "' could not be opened"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::fprintf(pfile,
      "\n# AthenaK particle moments at time= %.17e  nranks= %d  cycle= %d\n",
      pm->time, global_variable::nranks, pm->ncycle);
  std::fprintf(pfile,
      "# columns: species count mean_mu mean_mu2 anisotropy mean_dx mean_dy mean_dz "
      "mean_dx2 mean_dy2 mean_dz2 mean_dparallel mean_dparallel2 mean_speed2\n");

  int nspecies = pm->pmb_pack->ppart->nspecies;
  for (int sp=0; sp<nspecies; ++sp) {
    Real count = host_moments(sp,PMOM_COUNT);
    Real inv_count = (count > 0.0) ? 1.0/count : 0.0;
    Real mean_mu = host_moments(sp,PMOM_MU)*inv_count;
    Real mean_mu2 = host_moments(sp,PMOM_MU2)*inv_count;
    std::fprintf(pfile,
        "%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e "
        "%.17e %.17e %.17e\n",
        sp, count, mean_mu, mean_mu2, 3.0*mean_mu2 - 1.0,
        host_moments(sp,PMOM_DX)*inv_count,
        host_moments(sp,PMOM_DY)*inv_count,
        host_moments(sp,PMOM_DZ)*inv_count,
        host_moments(sp,PMOM_DX2)*inv_count,
        host_moments(sp,PMOM_DY2)*inv_count,
        host_moments(sp,PMOM_DZ2)*inv_count,
        host_moments(sp,PMOM_DPAR)*inv_count,
        host_moments(sp,PMOM_DPAR2)*inv_count,
        host_moments(sp,PMOM_SPEED2)*inv_count);
  }
  std::fclose(pfile);

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
