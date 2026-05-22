//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file df_prtcl.cpp
//! \brief writes pitch-angle histograms for cosmic-ray particles

#include <sys/stat.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdint>
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

constexpr int kParticleMomentCount = 24;

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
  PMOM_SPEED2=11,
  PMOM_VX=12,
  PMOM_VY=13,
  PMOM_VZ=14,
  PMOM_VX2=15,
  PMOM_VY2=16,
  PMOM_VZ2=17,
  PMOM_VXVY=18,
  PMOM_VXVZ=19,
  PMOM_VYVZ=20,
  PMOM_DXDY=21,
  PMOM_DXDZ=22,
  PMOM_DYDZ=23
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
    Kokkos::atomic_add(&local_moments(spec,PMOM_VX),pr(IPVX,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VY),pr(IPVY,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VZ),pr(IPVZ,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VX2),SQR(pr(IPVX,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VY2),SQR(pr(IPVY,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VZ2),SQR(pr(IPVZ,p)));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VXVY),pr(IPVX,p)*pr(IPVY,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VXVZ),pr(IPVX,p)*pr(IPVZ,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_VYVZ),pr(IPVY,p)*pr(IPVZ,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DXDY),pr(IPDX,p)*pr(IPDY,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DXDZ),pr(IPDX,p)*pr(IPDZ,p));
    Kokkos::atomic_add(&local_moments(spec,PMOM_DYDZ),pr(IPDY,p)*pr(IPDZ,p));
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
  std::fprintf(pfile,
      "# transport_columns: mean_vx mean_vy mean_vz mean_vx2 mean_vy2 mean_vz2 "
      "mean_vxvy mean_vxvz mean_vyvz mean_dxdy mean_dxdz mean_dydz\n");

  int nspecies = pm->pmb_pack->ppart->nspecies;
  for (int sp=0; sp<nspecies; ++sp) {
    Real count = host_moments(sp,PMOM_COUNT);
    Real inv_count = (count > 0.0) ? 1.0/count : 0.0;
    Real mean_mu = host_moments(sp,PMOM_MU)*inv_count;
    Real mean_mu2 = host_moments(sp,PMOM_MU2)*inv_count;
    std::fprintf(pfile,
        "%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e "
        "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e "
        "%.17e %.17e %.17e %.17e\n",
        sp, count, mean_mu, mean_mu2, 3.0*mean_mu2 - 1.0,
        host_moments(sp,PMOM_DX)*inv_count,
        host_moments(sp,PMOM_DY)*inv_count,
        host_moments(sp,PMOM_DZ)*inv_count,
        host_moments(sp,PMOM_DX2)*inv_count,
        host_moments(sp,PMOM_DY2)*inv_count,
        host_moments(sp,PMOM_DZ2)*inv_count,
        host_moments(sp,PMOM_DPAR)*inv_count,
        host_moments(sp,PMOM_DPAR2)*inv_count,
        host_moments(sp,PMOM_SPEED2)*inv_count,
        host_moments(sp,PMOM_VX)*inv_count,
        host_moments(sp,PMOM_VY)*inv_count,
        host_moments(sp,PMOM_VZ)*inv_count,
        host_moments(sp,PMOM_VX2)*inv_count,
        host_moments(sp,PMOM_VY2)*inv_count,
        host_moments(sp,PMOM_VZ2)*inv_count,
        host_moments(sp,PMOM_VXVY)*inv_count,
        host_moments(sp,PMOM_VXVZ)*inv_count,
        host_moments(sp,PMOM_VYVZ)*inv_count,
        host_moments(sp,PMOM_DXDY)*inv_count,
        host_moments(sp,PMOM_DXDZ)*inv_count,
        host_moments(sp,PMOM_DYDZ)*inv_count);
  }
  std::fclose(pfile);

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}

namespace {

enum ParticleSpectrumQuantity {
  PSPEC_P=0,
  PSPEC_E=1,
  PSPEC_LOGE=2,
  PSPEC_MU=3,
  PSPEC_MAGNETIC_MOMENT=4
};

} // namespace

ParticleSpectrumOutput::ParticleSpectrumOutput(ParameterInput *pin, Mesh *pm,
                                               OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("pspec",0775);

  quantity_name = pin->GetOrAddString(op.block_name,"quantity","p");
  if (quantity_name == "p" || quantity_name == "dNdp") {
    spectrum_quantity = PSPEC_P;
    quantity_name = "p";
  } else if (quantity_name == "E" || quantity_name == "energy" ||
             quantity_name == "dNdE") {
    spectrum_quantity = PSPEC_E;
    quantity_name = "E";
  } else if (quantity_name == "logE" || quantity_name == "dNdlogE") {
    spectrum_quantity = PSPEC_LOGE;
    quantity_name = "logE";
  } else if (quantity_name == "mu") {
    spectrum_quantity = PSPEC_MU;
    quantity_name = "mu";
  } else if (quantity_name == "magnetic_moment") {
    spectrum_quantity = PSPEC_MAGNETIC_MOMENT;
    quantity_name = "magnetic_moment";
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown pspec quantity '" << quantity_name
              << "' in output block '" << op.block_name << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  Real default_vmin = 0.0;
  Real default_vmax = 2.0;
  if (spectrum_quantity == PSPEC_MU) {
    default_vmin = -1.0;
    default_vmax = 1.0;
  } else if (spectrum_quantity == PSPEC_LOGE) {
    default_vmin = -12.0;
    default_vmax = 1.0;
  } else if (spectrum_quantity == PSPEC_E) {
    default_vmax = 2.0;
  } else if (spectrum_quantity == PSPEC_MAGNETIC_MOMENT) {
    default_vmax = 4.0;
  }
  vmin = pin->GetOrAddReal(op.block_name,"vmin",default_vmin);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",default_vmax);
  if (vmax <= vmin) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "pspec output block '" << op.block_name
              << "' requires vmax > vmin" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Kokkos::realloc(host_histogram, nspec, nbin);
  single_file_per_rank = pin->GetOrAddInteger(op.block_name,
                                              "pspec_single_file_per_rank",0);
  reduce_histogram = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (single_file_per_rank) {reduce_histogram = false;}
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "pspec/rank_%08d/",
                  global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleSpectrumOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  Real vmin_ = vmin;
  Real vmax_ = vmax;
  int nbin_ = nbin;
  int quantity_ = spectrum_quantity;

  DvceArray2D<int> local_histogram("local_particle_spectrum", nspecies, nbin);
  Kokkos::deep_copy(local_histogram, 0);
  par_for("part_spectrum_hist",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;

    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real speed = sqrt(speed2);
    Real energy = 0.5*speed2;
    Real bmag = sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real dotvb = pr(IPVX,p)*pr(IPBX,p) + pr(IPVY,p)*pr(IPBY,p) +
                 pr(IPVZ,p)*pr(IPBZ,p);
    Real value = speed;
    if (quantity_ == PSPEC_E) {
      value = energy;
    } else if (quantity_ == PSPEC_LOGE) {
      Real safe_energy = (energy > 1.0e-30) ? energy : 1.0e-30;
      value = log10(safe_energy);
    } else if (quantity_ == PSPEC_MU) {
      value = (bmag > 0.0 && speed > 0.0) ? dotvb/(speed*bmag) : 0.0;
    } else if (quantity_ == PSPEC_MAGNETIC_MOMENT) {
      Real vpar = (bmag > 0.0) ? dotvb/bmag : 0.0;
      Real vperp2 = speed2 - vpar*vpar;
      vperp2 = (vperp2 > 0.0) ? vperp2 : 0.0;
      value = (bmag > 0.0) ? vperp2/bmag : 0.0;
    }

    int ip = static_cast<int>((value - vmin_)/(vmax_ - vmin_)*nbin_);
    ip = (ip < 0) ? 0 : ip;
    ip = (ip > nbin_ - 1) ? nbin_ - 1 : ip;
    Kokkos::atomic_add(&local_histogram(spec,ip),1);
  });
  Kokkos::deep_copy(host_histogram, local_histogram);
#if MPI_PARALLEL_ENABLED
  if (reduce_histogram) {
    HostArray2D<int> global_histogram("global_particle_spectrum",
                                      nspecies, nbin);
    MPI_Reduce(host_histogram.data(), global_histogram.data(), nspecies*nbin,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_histogram, global_histogram);
    }
  }
#endif
  pp->LogPerformance("pspec_output", 0, 0, 0, 0, npart);
}

void ParticleSpectrumOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname("pspec/");
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + ".pspec";

  if (reduce_histogram && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  bool single = static_cast<bool>(single_file_per_rank) || reduce_histogram;
  if (single_file_per_rank || reduce_histogram || global_variable::my_rank == 0) {
    FILE *pfile = std::fopen(fname.c_str(),"a");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,
        "\n# AthenaK particle spectrum at time= %.17e  nranks= %d  cycle= %d\n",
        pm->time, global_variable::nranks, pm->ncycle);
    std::fprintf(pfile,
        "# metadata: quantity=%s nbin=%d vmin=%.17e vmax=%.17e reduce=%d\n",
        quantity_name.c_str(), nbin, vmin, vmax, reduce_histogram ? 1 : 0);
    std::fprintf(pfile,"# layout: species-major int32 histogram\n");
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
                << std::endl << "Particle spectrum output failed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (!single) {
      myoffset += static_cast<std::size_t>(global_variable::my_rank)*nout*sizeof(int);
    }
    if (partfile.Write_any_type_at_all(data,nout,myoffset,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle spectrum output failed" << std::endl;
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

enum ParticleJointSpectrumQuantity {
  PSPEC2_MU_P=0,
  PSPEC2_MU_E=1,
  PSPEC2_VPAR_VPERP=2
};

} // namespace

ParticleJointSpectrumOutput::ParticleJointSpectrumOutput(ParameterInput *pin, Mesh *pm,
                                                         OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("pspec2",0775);

  quantity_name = pin->GetOrAddString(op.block_name,"quantity","mu_p");
  if (quantity_name == "mu_p" || quantity_name == "mu,p" ||
      quantity_name == "f_mu_p") {
    spectrum_quantity = PSPEC2_MU_P;
    quantity_name = "mu_p";
  } else if (quantity_name == "mu_E" || quantity_name == "mu,e" ||
             quantity_name == "f_mu_E") {
    spectrum_quantity = PSPEC2_MU_E;
    quantity_name = "mu_E";
  } else if (quantity_name == "vpar_vperp" ||
             quantity_name == "v_parallel_v_perp") {
    spectrum_quantity = PSPEC2_VPAR_VPERP;
    quantity_name = "vpar_vperp";
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown pspec2 quantity '" << quantity_name
              << "' in output block '" << op.block_name << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  nbin1 = pin->GetOrAddInteger(op.block_name,"nbin1",32);
  nbin2 = pin->GetOrAddInteger(op.block_name,"nbin2",32);
  if (nbin1 <= 0 || nbin2 <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "pspec2 output block '" << op.block_name
              << "' requires positive nbin1 and nbin2" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real default_vmin1 = -1.0;
  Real default_vmax1 = 1.0;
  Real default_vmin2 = 0.0;
  Real default_vmax2 = 2.0;
  if (spectrum_quantity == PSPEC2_VPAR_VPERP) {
    default_vmin1 = -2.0;
    default_vmax1 = 2.0;
  }
  vmin1 = pin->GetOrAddReal(op.block_name,"vmin1",default_vmin1);
  vmax1 = pin->GetOrAddReal(op.block_name,"vmax1",default_vmax1);
  vmin2 = pin->GetOrAddReal(op.block_name,"vmin2",default_vmin2);
  vmax2 = pin->GetOrAddReal(op.block_name,"vmax2",default_vmax2);
  if (vmax1 <= vmin1 || vmax2 <= vmin2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "pspec2 output block '" << op.block_name
              << "' requires vmax1 > vmin1 and vmax2 > vmin2" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Kokkos::realloc(host_histogram, nspec, nbin1, nbin2);
  single_file_per_rank = pin->GetOrAddInteger(op.block_name,
                                              "pspec2_single_file_per_rank",0);
  reduce_histogram = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (single_file_per_rank) {reduce_histogram = false;}
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "pspec2/rank_%08d/",
                  global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleJointSpectrumOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  Real vmin1_ = vmin1;
  Real vmax1_ = vmax1;
  Real vmin2_ = vmin2;
  Real vmax2_ = vmax2;
  int nbin1_ = nbin1;
  int nbin2_ = nbin2;
  int quantity_ = spectrum_quantity;

  DvceArray3D<int> local_histogram("local_particle_joint_spectrum",
                                   nspecies, nbin1, nbin2);
  Kokkos::deep_copy(local_histogram, 0);
  par_for("part_joint_spectrum_hist",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;

    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real speed = sqrt(speed2);
    Real energy = 0.5*speed2;
    Real bmag = sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real dotvb = pr(IPVX,p)*pr(IPBX,p) + pr(IPVY,p)*pr(IPBY,p) +
                 pr(IPVZ,p)*pr(IPBZ,p);
    Real mu = (bmag > 0.0 && speed > 0.0) ? dotvb/(speed*bmag) : 0.0;
    Real vpar = (bmag > 0.0) ? dotvb/bmag : 0.0;
    Real vperp2 = speed2 - vpar*vpar;
    vperp2 = (vperp2 > 0.0) ? vperp2 : 0.0;

    Real value1 = mu;
    Real value2 = speed;
    if (quantity_ == PSPEC2_MU_E) {
      value2 = energy;
    } else if (quantity_ == PSPEC2_VPAR_VPERP) {
      value1 = vpar;
      value2 = sqrt(vperp2);
    }

    int i1 = static_cast<int>((value1 - vmin1_)/(vmax1_ - vmin1_)*nbin1_);
    int i2 = static_cast<int>((value2 - vmin2_)/(vmax2_ - vmin2_)*nbin2_);
    i1 = (i1 < 0) ? 0 : i1;
    i1 = (i1 > nbin1_ - 1) ? nbin1_ - 1 : i1;
    i2 = (i2 < 0) ? 0 : i2;
    i2 = (i2 > nbin2_ - 1) ? nbin2_ - 1 : i2;
    Kokkos::atomic_add(&local_histogram(spec,i1,i2),1);
  });
  Kokkos::deep_copy(host_histogram, local_histogram);
#if MPI_PARALLEL_ENABLED
  if (reduce_histogram) {
    HostArray3D<int> global_histogram("global_particle_joint_spectrum",
                                      nspecies, nbin1, nbin2);
    MPI_Reduce(host_histogram.data(), global_histogram.data(),
               nspecies*nbin1*nbin2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_histogram, global_histogram);
    }
  }
#endif
  pp->LogPerformance("pspec2_output", 0, 0, 0, 0, npart);
}

void ParticleJointSpectrumOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname("pspec2/");
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + ".pspec2";

  if (reduce_histogram && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  bool single = static_cast<bool>(single_file_per_rank) || reduce_histogram;
  if (single_file_per_rank || reduce_histogram || global_variable::my_rank == 0) {
    FILE *pfile = std::fopen(fname.c_str(),"a");
    if (pfile == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,
        "\n# AthenaK particle joint spectrum at time= %.17e  nranks= %d  cycle= %d\n",
        pm->time, global_variable::nranks, pm->ncycle);
    std::fprintf(pfile,
        "# metadata: quantity=%s nbin1=%d nbin2=%d vmin1=%.17e vmax1=%.17e "
        "vmin2=%.17e vmax2=%.17e reduce=%d\n",
        quantity_name.c_str(), nbin1, nbin2, vmin1, vmax1, vmin2, vmax2,
        reduce_histogram ? 1 : 0);
    std::fprintf(pfile,"# layout: species-major int32 histogram, bin1 then bin2\n");
    std::fclose(pfile);
  }

  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append, single);
  std::size_t header_offset = partfile.GetPosition(single);

  int nspecies = pm->pmb_pack->ppart->nspecies;
  int nout = nspecies*nbin1*nbin2;
  int *data = new int[std::max(1,nout)];
  for (int sp=0; sp<nspecies; ++sp) {
    for (int bin1=0; bin1<nbin1; ++bin1) {
      for (int bin2=0; bin2<nbin2; ++bin2) {
        int n = (sp*nbin1 + bin1)*nbin2 + bin2;
        data[n] = host_histogram(sp,bin1,bin2);
      }
    }
  }

  std::size_t myoffset = header_offset;
  if (reduce_histogram) {
    if (partfile.Write_any_type(data,nout,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle joint spectrum output failed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (!single) {
      myoffset += static_cast<std::size_t>(global_variable::my_rank)*nout*sizeof(int);
    }
    if (partfile.Write_any_type_at_all(data,nout,myoffset,"int",single) !=
        static_cast<std::size_t>(nout)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle joint spectrum output failed" << std::endl;
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

enum ParticleSampleFieldKind {
  PSAMP_REAL_FIELD=0,
  PSAMP_DERIVED_FIELD=1
};

enum ParticleSampleDerivedIndex {
  PSAMP_SPEED=0,
  PSAMP_ENERGY=1,
  PSAMP_BMAG=2,
  PSAMP_MU=3,
  PSAMP_VPAR=4,
  PSAMP_VPERP=5,
  PSAMP_MAGNETIC_MOMENT=6
};

std::string NormalizeParticleSampleField(std::string name) {
  name.erase(std::remove_if(name.begin(), name.end(),
      [](unsigned char c) { return std::isspace(c); }), name.end());
  std::transform(name.begin(), name.end(), name.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return name;
}

uint64_t ParticleSampleHash(int species, int tag) {
  uint64_t x = (static_cast<uint64_t>(static_cast<uint32_t>(species)) << 32) ^
               static_cast<uint32_t>(tag);
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

bool ParticleSampleSelected(int species, int tag, int sample_species,
                            int sample_stride, int sample_offset) {
  if (sample_species >= 0 && species != sample_species) {
    return false;
  }
  return static_cast<int>(ParticleSampleHash(species, tag) %
                         static_cast<uint64_t>(sample_stride)) == sample_offset;
}

ParticleSampleField MakeParticleSampleField(const std::string &name,
                                            const std::string &block_name) {
  std::string field = NormalizeParticleSampleField(name);
  if (field == "x" || field == "x1") {
    return ParticleSampleField("x", PSAMP_REAL_FIELD, IPX);
  } else if (field == "y" || field == "x2") {
    return ParticleSampleField("y", PSAMP_REAL_FIELD, IPY);
  } else if (field == "z" || field == "x3") {
    return ParticleSampleField("z", PSAMP_REAL_FIELD, IPZ);
  } else if (field == "vx" || field == "v1") {
    return ParticleSampleField("vx", PSAMP_REAL_FIELD, IPVX);
  } else if (field == "vy" || field == "v2") {
    return ParticleSampleField("vy", PSAMP_REAL_FIELD, IPVY);
  } else if (field == "vz" || field == "v3") {
    return ParticleSampleField("vz", PSAMP_REAL_FIELD, IPVZ);
  } else if (field == "m" || field == "mass") {
    return ParticleSampleField("mass", PSAMP_REAL_FIELD, IPM);
  } else if (field == "bx" || field == "b1") {
    return ParticleSampleField("bx", PSAMP_REAL_FIELD, IPBX);
  } else if (field == "by" || field == "b2") {
    return ParticleSampleField("by", PSAMP_REAL_FIELD, IPBY);
  } else if (field == "bz" || field == "b3") {
    return ParticleSampleField("bz", PSAMP_REAL_FIELD, IPBZ);
  } else if (field == "dx") {
    return ParticleSampleField("dx", PSAMP_REAL_FIELD, IPDX);
  } else if (field == "dy") {
    return ParticleSampleField("dy", PSAMP_REAL_FIELD, IPDY);
  } else if (field == "dz") {
    return ParticleSampleField("dz", PSAMP_REAL_FIELD, IPDZ);
  } else if (field == "dpar" || field == "dparallel") {
    return ParticleSampleField("dpar", PSAMP_REAL_FIELD, IPDB);
  } else if (field == "speed" || field == "p") {
    return ParticleSampleField("speed", PSAMP_DERIVED_FIELD, PSAMP_SPEED);
  } else if (field == "energy" || field == "e") {
    return ParticleSampleField("energy", PSAMP_DERIVED_FIELD, PSAMP_ENERGY);
  } else if (field == "bmag") {
    return ParticleSampleField("bmag", PSAMP_DERIVED_FIELD, PSAMP_BMAG);
  } else if (field == "mu") {
    return ParticleSampleField("mu", PSAMP_DERIVED_FIELD, PSAMP_MU);
  } else if (field == "vpar" || field == "vparallel") {
    return ParticleSampleField("vpar", PSAMP_DERIVED_FIELD, PSAMP_VPAR);
  } else if (field == "vperp") {
    return ParticleSampleField("vperp", PSAMP_DERIVED_FIELD, PSAMP_VPERP);
  } else if (field == "magnetic_moment") {
    return ParticleSampleField("magnetic_moment", PSAMP_DERIVED_FIELD,
                               PSAMP_MAGNETIC_MOMENT);
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Unknown psamp field '" << name
            << "' in output block '" << block_name << "'" << std::endl;
  std::exit(EXIT_FAILURE);
}

Real ParticleSampleDerivedValue(const HostArray2D<Real> &pr, int p, int index) {
  Real vx = pr(IPVX,p);
  Real vy = pr(IPVY,p);
  Real vz = pr(IPVZ,p);
  Real bx = pr(IPBX,p);
  Real by = pr(IPBY,p);
  Real bz = pr(IPBZ,p);
  Real speed2 = SQR(vx) + SQR(vy) + SQR(vz);
  Real speed = std::sqrt(speed2);
  Real bmag = std::sqrt(SQR(bx) + SQR(by) + SQR(bz));
  Real dotvb = vx*bx + vy*by + vz*bz;
  Real vpar = (bmag > 0.0) ? dotvb/bmag : 0.0;
  Real vperp2 = speed2 - vpar*vpar;
  vperp2 = (vperp2 > 0.0) ? vperp2 : 0.0;

  if (index == PSAMP_SPEED) {
    return speed;
  } else if (index == PSAMP_ENERGY) {
    return 0.5*speed2;
  } else if (index == PSAMP_BMAG) {
    return bmag;
  } else if (index == PSAMP_MU) {
    return (bmag > 0.0 && speed > 0.0) ? dotvb/(speed*bmag) : 0.0;
  } else if (index == PSAMP_VPAR) {
    return vpar;
  } else if (index == PSAMP_VPERP) {
    return std::sqrt(vperp2);
  } else if (index == PSAMP_MAGNETIC_MOMENT) {
    return (bmag > 0.0) ? vperp2/bmag : 0.0;
  }
  return 0.0;
}

} // namespace

ParticleSampleOutput::ParticleSampleOutput(ParameterInput *pin, Mesh *pm,
                                           OutputParameters op) :
  BaseTypeOutput(pin, pm, op),
  npout_thisrank(0) {
  mkdir("psamp",0775);
  char rank_dir[32];
  std::snprintf(rank_dir, sizeof(rank_dir), "psamp/rank_%08d/",
                global_variable::my_rank);
  mkdir(rank_dir,0775);

  sample_species = pin->GetOrAddInteger(op.block_name,"species",-1);
  sample_stride = pin->GetOrAddInteger(op.block_name,"sample_stride",1);
  sample_offset = pin->GetOrAddInteger(op.block_name,"sample_offset",0);
  if (sample_stride <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "psamp output block '" << op.block_name
              << "' requires sample_stride > 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sample_offset < 0 || sample_offset >= sample_stride) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "psamp output block '" << op.block_name
              << "' requires 0 <= sample_offset < sample_stride" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string field_text = pin->GetOrAddString(op.block_name, "fields",
      "x,y,z,vx,vy,vz,bx,by,bz,dx,dy,dz,dpar,mass");
  std::stringstream field_stream(field_text);
  std::string field_name;
  while (std::getline(field_stream, field_name, ',')) {
    std::string normalized = NormalizeParticleSampleField(field_name);
    if (normalized.empty()) {
      continue;
    }
    if (normalized == "pgid" || normalized == "gid" ||
        normalized == "tag" || normalized == "species") {
      continue;
    }
    sample_fields.push_back(MakeParticleSampleField(field_name, op.block_name));
  }
  if (sample_fields.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "psamp output block '" << op.block_name
              << "' must include at least one real or derived field" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  out_params.last_time = pm->time;
}

void ParticleSampleOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pp->nprtcl_thispack;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  if (npout_thisrank > 0) {
    Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
    Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
  }
  pp->LogPerformance("psamp_output", 0, 0, 0, 0, npout_thisrank);
}

void ParticleSampleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  char rank_dir[32];
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
  std::string fname("psamp/");
  fname += std::string(rank_dir);
  fname += out_params.file_basename + ".psamp";

  FILE *pfile = std::fopen(fname.c_str(),"a");
  if (pfile == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Output file '" << fname << "' could not be opened"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int sample_count = 0;
  for (int p=0; p<npout_thisrank; ++p) {
    int spec = outpart_idata(PSP,p);
    int tag = outpart_idata(PTAG,p);
    if (ParticleSampleSelected(spec, tag, sample_species, sample_stride,
                               sample_offset)) {
      ++sample_count;
    }
  }

  std::fprintf(pfile,
      "\n# AthenaK particle sample at time= %.17e  nranks= %d  cycle= %d\n",
      pm->time, global_variable::nranks, pm->ncycle);
  std::fprintf(pfile,
      "# metadata: selected_species=%d sample_stride=%d sample_offset=%d "
      "sample_count=%d field_count=%zu\n",
      sample_species, sample_stride, sample_offset, sample_count,
      sample_fields.size());
  std::fprintf(pfile,
      "# selection: splitmix64(species,tag) %% sample_stride == sample_offset\n");
  std::fprintf(pfile, "# columns: rank pgid tag species");
  for (const auto &field : sample_fields) {
    std::fprintf(pfile, " %s", field.label.c_str());
  }
  std::fprintf(pfile, "\n");

  for (int p=0; p<npout_thisrank; ++p) {
    int pgid = outpart_idata(PGID,p);
    int tag = outpart_idata(PTAG,p);
    int spec = outpart_idata(PSP,p);
    if (!ParticleSampleSelected(spec, tag, sample_species, sample_stride,
                                sample_offset)) {
      continue;
    }
    std::fprintf(pfile, "%d %d %d %d", global_variable::my_rank, pgid, tag, spec);
    for (const auto &field : sample_fields) {
      Real value = (field.kind == PSAMP_REAL_FIELD)
                 ? outpart_rdata(field.index,p)
                 : ParticleSampleDerivedValue(outpart_rdata, p, field.index);
      std::fprintf(pfile, " %.17e", value);
    }
    std::fprintf(pfile, "\n");
  }
  std::fclose(pfile);

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
