//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dxhist_prtcl.cpp
//! \brief writes displacement histograms for cosmic-ray particles

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "particles/relativistic_state.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

bool IsRelativisticParticleOutput(const particles::Particles *pp) {
  return pp->pusher == ParticlesPusher::relativistic_hc;
}

void WriteRelativisticParticleMetadata(FILE *pfile, const particles::Particles *pp) {
  if (!IsRelativisticParticleOutput(pp)) {return;}
  std::fprintf(pfile,
      "# metadata: schema=akcr_particle_output_v1 mode=relativistic_hc "
      "units=code_model c_model=%.17e alpha_s=%.17e\n",
      pp->c_model, pp->alpha_s);
}

void RequireHistogramGeometry(const std::string &block_name, const std::string &kind,
                              const int nbin, const Real vmin, const Real vmax) {
  if (nbin <= 0 || !std::isfinite(vmin) || !std::isfinite(vmax) ||
      !std::isfinite(vmax - vmin) || vmax <= vmin) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << kind << " output block '" << block_name
              << "' requires positive nbin and finite non-overflowing vmax > vmin"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

KOKKOS_INLINE_FUNCTION
int HistogramBin(const Real value, const Real vmin, const Real vmax, const int nbin) {
  if (value <= vmin) {return 0;}
  if (value >= vmax) {return nbin - 1;}
  return static_cast<int>((value - vmin)/(vmax - vmin)*nbin);
}

[[noreturn]] void FatalParticleDiagnostic(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

ParticleDxHistOutput::ParticleDxHistOutput(ParameterInput *pin, Mesh *pm,
                                           OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("dxh",0775);
  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  vmin = pin->GetOrAddReal(op.block_name,"vmin",-1.0);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",1.0);
  RequireHistogramGeometry(op.block_name, "dxh", nbin, vmin, vmax);
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
  DvceArray1D<int> derived_failure("local_dxh_hist_failure", 1);
  Kokkos::deep_copy(local_histogram, 0);
  Kokkos::deep_copy(derived_failure, 0);
  par_for("part_dxh_hist",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    Real dx[3] = {pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p)};
    for (int d=0; d<3; ++d) {
      if (!particles::relativistic::IsFinite(dx[d])) {
        Kokkos::atomic_max(&derived_failure(0), 1);
        return;
      }
      int ip = HistogramBin(dx[d], vmin_, vmax_, nbin_);
      Kokkos::atomic_add(&local_histogram(3*spec + d,ip),1);
    }
  });
  auto derived_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), derived_failure);
  if (derived_failure_h(0) != 0) {
    FatalParticleDiagnostic("particle dxh derived quantity evaluation failed");
  }
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
  pp->LogPerformance("dxh_output", 0, 0, 0, 0, npart);
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
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (IsRelativisticParticleOutput(pp)) {
      WriteRelativisticParticleMetadata(pfile, pp);
      std::fprintf(pfile,
          "# metadata: quantity=dx_dy_dz quantity_units=code_length "
          "histogram_units=particle_count nspecies=%d nbin=%d "
          "vmin=%.17e vmax=%.17e reduce=%d\n",
          pp->nspecies, nbin, vmin, vmax, reduce_histogram ? 1 : 0);
      std::fprintf(pfile,
          "# layout: species-major int32 histograms ordered dx, dy, dz\n");
    }
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

ParticleScalarDxHistOutput::ParticleScalarDxHistOutput(ParameterInput *pin, Mesh *pm,
                                                       OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  use_parallel_displacement = (op.file_type.compare("dparh") == 0);
  output_dir = use_parallel_displacement ? "dparh" : "drh";
  header_marker = use_parallel_displacement ?
      "# AthenaK particle parallel displacement histogram" :
      "# AthenaK particle scalar displacement histogram";
  mkdir(output_dir.c_str(),0775);
  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  vmin = pin->GetOrAddReal(op.block_name,"vmin",
                           use_parallel_displacement ? -1.0 : 0.0);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",1.0);
  RequireHistogramGeometry(op.block_name, output_dir, nbin, vmin, vmax);
  Kokkos::realloc(host_histogram, nspec, nbin);
  single_file_per_rank = pin->GetOrAddInteger(
      op.block_name, (output_dir + "_single_file_per_rank").c_str(), 0);
  reduce_histogram = pin->GetOrAddBoolean(op.block_name,"reduce",true);
  if (single_file_per_rank) {reduce_histogram = false;}
  if (single_file_per_rank) {
    char rank_dir[64];
    std::snprintf(rank_dir, sizeof(rank_dir), "%s/rank_%08d/",
                  output_dir.c_str(), global_variable::my_rank);
    mkdir(rank_dir,0775);
  }
  out_params.last_time = pm->time;
}

void ParticleScalarDxHistOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  int npart = pp->nprtcl_thispack;
  int nspecies = pp->nspecies;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  Real vmin_ = vmin;
  Real vmax_ = vmax;
  int nbin_ = nbin;
  bool use_parallel = use_parallel_displacement;

  DvceArray2D<int> local_histogram("local_scalar_dxh_hist", nspecies, nbin);
  DvceArray1D<int> derived_failure("local_scalar_dxh_hist_failure", 1);
  Kokkos::deep_copy(local_histogram, 0);
  Kokkos::deep_copy(derived_failure, 0);
  par_for("part_scalar_dxh_hist",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    Real val = pr(IPDB,p);
    if (!use_parallel) {
      val = particles::relativistic::ScaledNorm3(
          {pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p)});
    }
    if (!particles::relativistic::IsFinite(val)) {
      Kokkos::atomic_max(&derived_failure(0), 1);
      return;
    }
    int ip = HistogramBin(val, vmin_, vmax_, nbin_);
    Kokkos::atomic_add(&local_histogram(spec,ip),1);
  });
  auto derived_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), derived_failure);
  if (derived_failure_h(0) != 0) {
    FatalParticleDiagnostic("particle scalar displacement evaluation failed");
  }
  Kokkos::deep_copy(host_histogram, local_histogram);
#if MPI_PARALLEL_ENABLED
  if (reduce_histogram) {
    HostArray2D<int> global_histogram("global_scalar_dxh_hist", nspecies, nbin);
    MPI_Reduce(host_histogram.data(), global_histogram.data(), nspecies*nbin,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      Kokkos::deep_copy(host_histogram, global_histogram);
    }
  }
#endif
  pp->LogPerformance(output_dir + "_output", 0, 0, 0, 0, npart);
}

void ParticleScalarDxHistOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname(output_dir + "/");
  if (single_file_per_rank) {
    char rank_dir[32];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname += std::string(rank_dir);
  }
  fname += out_params.file_basename + "." + output_dir;

  if (reduce_histogram && global_variable::my_rank != 0) {
    bool reset_time = (out_params.last_time < 0.0 ||
                       out_params.last_time == pm->time);
    out_params.last_time = reset_time ? pm->time : out_params.last_time + out_params.dt;
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  bool single = static_cast<bool>(single_file_per_rank) || reduce_histogram;
  if (single_file_per_rank || reduce_histogram || global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << "\n" << header_marker << " at time= " << pm->time
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
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (IsRelativisticParticleOutput(pp)) {
      WriteRelativisticParticleMetadata(pfile, pp);
      if (use_parallel_displacement) {
        std::fprintf(pfile,
            "# metadata: quantity=dparallel quantity_units=code_length "
            "histogram_units=particle_count "
            "definition=accumulated_midpoint_sum_dx_dot_Bhat "
            "nspecies=%d nbin=%d vmin=%.17e vmax=%.17e reduce=%d\n",
            pp->nspecies, nbin, vmin, vmax, reduce_histogram ? 1 : 0);
      } else {
        std::fprintf(pfile,
            "# metadata: quantity=displacement_norm quantity_units=code_length "
            "histogram_units=particle_count nspecies=%d nbin=%d "
            "vmin=%.17e vmax=%.17e reduce=%d\n",
            pp->nspecies, nbin, vmin, vmax, reduce_histogram ? 1 : 0);
      }
      std::fprintf(pfile,"# layout: species-major int32 histogram\n");
    }
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
                << std::endl << "Particle scalar displacement output failed"
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
                << std::endl << "Particle scalar displacement output failed"
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
