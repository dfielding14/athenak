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
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

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

constexpr Real kParticleSpectrumLogFloor = 1.0e-30;

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

[[noreturn]] void FatalParticleDiagnostic(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void RequireHistogramGeometry(const std::string &block_name, const std::string &kind,
                              const int nbin, const Real vmin, const Real vmax) {
  if (nbin <= 0 || !std::isfinite(vmin) || !std::isfinite(vmax) ||
      !std::isfinite(vmax - vmin) || vmax <= vmin) {
    FatalParticleDiagnostic(kind + " output block '" + block_name +
        "' requires positive nbin and finite non-overflowing vmax > vmin");
  }
}

KOKKOS_INLINE_FUNCTION
int HistogramBin(const Real value, const Real vmin, const Real vmax, const int nbin) {
  if (value <= vmin) {return 0;}
  if (value >= vmax) {return nbin - 1;}
  return static_cast<int>((value - vmin)/(vmax - vmin)*nbin);
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
particles::relativistic::Vector3 ParticleW(const RealView &pr, const int p) {
  return {pr(IPWX,p), pr(IPWY,p), pr(IPWZ,p)};
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticlePhysicalSpeed(const RealView &pr, const int p) {
  return particles::relativistic::ScaledNorm3(
      {pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p)});
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticleBmag(const RealView &pr, const int p) {
  return particles::relativistic::ScaledNorm3(
      {pr(IPBX,p), pr(IPBY,p), pr(IPBZ,p)});
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticleParallelVelocity(const RealView &pr, const int p, const Real bmag) {
  if (bmag <= 0.0) {return 0.0;}
  return pr(IPVX,p)*(pr(IPBX,p)/bmag) + pr(IPVY,p)*(pr(IPBY,p)/bmag) +
         pr(IPVZ,p)*(pr(IPBZ,p)/bmag);
}

KOKKOS_INLINE_FUNCTION
Real PerpendicularMagnitude(const particles::relativistic::Vector3 &value,
                            const Real bx, const Real by, const Real bz,
                            const Real bmag) {
  if (bmag <= 0.0) {return 0.0;}
  Real bhx = bx/bmag;
  Real bhy = by/bmag;
  Real bhz = bz/bmag;
  particles::relativistic::Vector3 cross{
    value.y*bhz - value.z*bhy,
    value.z*bhx - value.x*bhz,
    value.x*bhy - value.y*bhx
  };
  return particles::relativistic::ScaledNorm3(cross);
}

KOKKOS_INLINE_FUNCTION
Real PositiveSquareRatio(const Real numerator, const Real denominator) {
  if (!Kokkos::isfinite(numerator) || numerator < 0.0 ||
      !Kokkos::isfinite(denominator) || denominator <= 0.0) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  if (numerator == 0.0) {return 0.0;}
  Real numerator_exponent = Kokkos::logb(numerator);
  Real denominator_exponent = Kokkos::logb(denominator);
  Real numerator_scaled = numerator/Kokkos::exp2(numerator_exponent);
  Real denominator_scaled = denominator/Kokkos::exp2(denominator_exponent);
  Real scaled = numerator_scaled*numerator_scaled/denominator_scaled;
  Real scaled_exponent = Kokkos::logb(scaled);
  scaled /= Kokkos::exp2(scaled_exponent);
  return scaled*Kokkos::exp2(
      2.0*numerator_exponent - denominator_exponent + scaled_exponent);
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticlePitchAngleCosine(const RealView &pr, const int p, const Real speed,
                              const Real bmag) {
  if (speed <= 0.0 || bmag <= 0.0) {return 0.0;}
  return ParticleParallelVelocity(pr, p, bmag)/speed;
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticlePerpendicularVelocity(const RealView &pr, const int p, const Real bmag) {
  return PerpendicularMagnitude(
      {pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p)},
      pr(IPBX,p), pr(IPBY,p), pr(IPBZ,p), bmag);
}

template<typename RealView>
KOKKOS_INLINE_FUNCTION
Real ParticleVelocityMagneticMomentProxy(const RealView &pr, const int p,
                                         const Real bmag) {
  if (bmag <= 0.0) {return 0.0;}
  Real vperp = ParticlePerpendicularVelocity(pr, p, bmag);
  return PositiveSquareRatio(vperp, bmag);
}

} // namespace

ParticleDFOutput::ParticleDFOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  mkdir("df",0775);
  nbin = pin->GetOrAddInteger(op.block_name,"nbin",100);
  vmin = pin->GetOrAddReal(op.block_name,"vmin",-1.0);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",1.0);
  RequireHistogramGeometry(op.block_name, "df", nbin, vmin, vmax);
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
  bool relativistic_mode = IsRelativisticParticleOutput(pp);

  DvceArray2D<int> local_histogram("local_df_hist", nspecies, nbin);
  DvceArray1D<int> relativistic_failure("local_df_hist_failure", 1);
  Kokkos::deep_copy(local_histogram, 0);
  Kokkos::deep_copy(relativistic_failure, 0);
  par_for("part_df_hist",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    Real bmag = ParticleBmag(pr, p);
    Real vmag = ParticlePhysicalSpeed(pr, p);
    if (!particles::relativistic::IsFinite(bmag) ||
        !particles::relativistic::IsFinite(vmag) ||
        (relativistic_mode && (bmag <= 0.0 || vmag <= 0.0))) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }
    Real mu = ParticlePitchAngleCosine(pr, p, vmag, bmag);
    if (!particles::relativistic::IsFinite(mu)) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    int ip = HistogramBin(mu, vmin_, vmax_, nbin_);
    Kokkos::atomic_add(&local_histogram(spec,ip),1);
  });
  auto relativistic_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), relativistic_failure);
  if (relativistic_failure_h(0) != 0) {
    FatalParticleDiagnostic("particle df derived quantity evaluation failed");
  }
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
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (IsRelativisticParticleOutput(pp)) {
      WriteRelativisticParticleMetadata(pfile, pp);
      std::fprintf(pfile,
          "# metadata: quantity=mu quantity_units=dimensionless "
          "histogram_units=particle_count nspecies=%d nbin=%d "
          "vmin=%.17e vmax=%.17e reduce=%d\n",
          pp->nspecies, nbin, vmin, vmax, reduce_histogram ? 1 : 0);
      std::fprintf(pfile,
          "# layout: species-major int32 histogram; mu=dot(v,B)/(|v||B|)\n");
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
  bool relativistic_mode = IsRelativisticParticleOutput(pp);

  DvceArray2D<Real> local_moments("local_particle_moments", nspecies, nmoment);
  DvceArray1D<int> relativistic_failure("local_particle_moments_failure", 1);
  Kokkos::deep_copy(local_moments, 0.0);
  Kokkos::deep_copy(relativistic_failure, 0);
  par_for("part_moments",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;
    Real bmag = ParticleBmag(pr, p);
    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real vmag = ParticlePhysicalSpeed(pr, p);
    Real mu = ParticlePitchAngleCosine(pr, p, vmag, bmag);
    Real values[] = {
      mu, mu*mu,
      pr(IPDX,p), pr(IPDY,p), pr(IPDZ,p),
      SQR(pr(IPDX,p)), SQR(pr(IPDY,p)), SQR(pr(IPDZ,p)),
      pr(IPDB,p), SQR(pr(IPDB,p)), speed2,
      pr(IPVX,p), pr(IPVY,p), pr(IPVZ,p),
      SQR(pr(IPVX,p)), SQR(pr(IPVY,p)), SQR(pr(IPVZ,p)),
      pr(IPVX,p)*pr(IPVY,p), pr(IPVX,p)*pr(IPVZ,p),
      pr(IPVY,p)*pr(IPVZ,p), pr(IPDX,p)*pr(IPDY,p),
      pr(IPDX,p)*pr(IPDZ,p), pr(IPDY,p)*pr(IPDZ,p)
    };
    if (relativistic_mode) {
      if (!particles::relativistic::IsFinite(bmag) ||
          !particles::relativistic::IsFinite(vmag) ||
          bmag <= 0.0 || vmag <= 0.0) {
        Kokkos::atomic_max(&relativistic_failure(0), 1);
        return;
      }
      for (int n=0; n<kParticleMomentCount - 1; ++n) {
        if (!particles::relativistic::IsFinite(values[n])) {
          Kokkos::atomic_max(&relativistic_failure(0), 1);
          return;
        }
      }
    }

    Kokkos::atomic_add(&local_moments(spec,PMOM_COUNT),1.0);
    for (int n=1; n<kParticleMomentCount; ++n) {
      Kokkos::atomic_add(&local_moments(spec,n),values[n-1]);
    }
  });
  auto relativistic_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), relativistic_failure);
  if (relativistic_failure_h(0) != 0) {
    FatalParticleDiagnostic("relativistic_hc pmom derived quantity evaluation failed");
  }
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
  particles::Particles *pp = pm->pmb_pack->ppart;
  if (IsRelativisticParticleOutput(pp)) {
    WriteRelativisticParticleMetadata(pfile, pp);
    std::fprintf(pfile,
        "# metadata: quantity_basis=physical_velocity_shadow "
        "mu_units=dimensionless displacement_units=code_length "
        "displacement_second_moment_units=code_length_squared "
        "velocity_units=code_velocity "
        "velocity_second_moment_units=code_velocity_squared "
        "dparallel_definition=accumulated_midpoint_sum_dx_dot_Bhat\n");
  }
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
  PSPEC_SPEED=0,
  PSPEC_LEGACY_ENERGY=1,
  PSPEC_LEGACY_LOGE=2,
  PSPEC_MU=3,
  PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY=4,
  PSPEC_WMAG=5,
  PSPEC_KINETIC_ENERGY_MODEL=6,
  PSPEC_LOG10_KINETIC_ENERGY_MODEL=7
};

const char *ParticleSpectrumQuantityUnits(const int quantity) {
  if (quantity == PSPEC_MU) {return "dimensionless";}
  if (quantity == PSPEC_LEGACY_ENERGY ||
      quantity == PSPEC_KINETIC_ENERGY_MODEL) {
    return "code_velocity_squared";
  }
  if (quantity == PSPEC_LEGACY_LOGE ||
      quantity == PSPEC_LOG10_KINETIC_ENERGY_MODEL) {
    return "log10_code_velocity_squared";
  }
  if (quantity == PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY) {
    return "code_velocity_squared_per_code_B";
  }
  return "code_velocity";
}

} // namespace

ParticleSpectrumOutput::ParticleSpectrumOutput(ParameterInput *pin, Mesh *pm,
                                               OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  bool relativistic_mode =
      IsRelativisticParticleOutput(pm->pmb_pack->ppart);
  mkdir("pspec",0775);

  if (relativistic_mode &&
      !pin->DoesParameterExist(op.block_name,"quantity")) {
    FatalParticleDiagnostic("relativistic_hc pspec output block '" +
        op.block_name + "' requires explicit quantity");
  }
  quantity_name = relativistic_mode ?
      pin->GetString(op.block_name,"quantity") :
      pin->GetOrAddString(op.block_name,"quantity","p");
  if (relativistic_mode) {
    if (quantity_name == "speed") {
      spectrum_quantity = PSPEC_SPEED;
    } else if (quantity_name == "wmag") {
      spectrum_quantity = PSPEC_WMAG;
    } else if (quantity_name == "kinetic_energy_model") {
      spectrum_quantity = PSPEC_KINETIC_ENERGY_MODEL;
    } else if (quantity_name == "log10_kinetic_energy_model") {
      spectrum_quantity = PSPEC_LOG10_KINETIC_ENERGY_MODEL;
    } else if (quantity_name == "mu") {
      spectrum_quantity = PSPEC_MU;
    } else if (quantity_name == "velocity_magnetic_moment_proxy") {
      spectrum_quantity = PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY;
    } else if (quantity_name == "p" || quantity_name == "dNdp" ||
               quantity_name == "E" || quantity_name == "energy" ||
               quantity_name == "dNdE" || quantity_name == "logE" ||
               quantity_name == "dNdlogE" ||
               quantity_name == "magnetic_moment") {
      FatalParticleDiagnostic("relativistic_hc pspec output block '" +
          op.block_name + "' rejects ambiguous quantity alias '" +
          quantity_name + "'");
    } else {
      FatalParticleDiagnostic("Unknown relativistic_hc pspec quantity '" +
          quantity_name + "' in output block '" + op.block_name + "'");
    }
  } else if (quantity_name == "p" || quantity_name == "dNdp") {
    spectrum_quantity = PSPEC_SPEED;
    quantity_name = "p";
  } else if (quantity_name == "speed") {
    spectrum_quantity = PSPEC_SPEED;
  } else if (quantity_name == "E" || quantity_name == "energy" ||
             quantity_name == "dNdE") {
    spectrum_quantity = PSPEC_LEGACY_ENERGY;
    quantity_name = "E";
  } else if (quantity_name == "logE" || quantity_name == "dNdlogE") {
    spectrum_quantity = PSPEC_LEGACY_LOGE;
    quantity_name = "logE";
  } else if (quantity_name == "mu") {
    spectrum_quantity = PSPEC_MU;
    quantity_name = "mu";
  } else if (quantity_name == "magnetic_moment") {
    spectrum_quantity = PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY;
    quantity_name = "magnetic_moment";
  } else if (quantity_name == "velocity_magnetic_moment_proxy") {
    spectrum_quantity = PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY;
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
  } else if (spectrum_quantity == PSPEC_LEGACY_LOGE ||
             spectrum_quantity == PSPEC_LOG10_KINETIC_ENERGY_MODEL) {
    default_vmin = -12.0;
    default_vmax = 1.0;
  } else if (spectrum_quantity == PSPEC_LEGACY_ENERGY ||
             spectrum_quantity == PSPEC_KINETIC_ENERGY_MODEL) {
    default_vmax = 2.0;
  } else if (spectrum_quantity == PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY) {
    default_vmax = 4.0;
  }
  vmin = pin->GetOrAddReal(op.block_name,"vmin",default_vmin);
  vmax = pin->GetOrAddReal(op.block_name,"vmax",default_vmax);
  RequireHistogramGeometry(op.block_name, "pspec", nbin, vmin, vmax);

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
  bool relativistic_mode = IsRelativisticParticleOutput(pp);
  Real c_model = pp->c_model;
  Real log_floor = kParticleSpectrumLogFloor;

  DvceArray2D<int> local_histogram("local_particle_spectrum", nspecies, nbin);
  DvceArray1D<int> relativistic_failure("local_particle_spectrum_failure", 1);
  Kokkos::deep_copy(local_histogram, 0);
  Kokkos::deep_copy(relativistic_failure, 0);
  par_for("part_spectrum_hist",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;

    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real speed = relativistic_mode ? ParticlePhysicalSpeed(pr, p) : sqrt(speed2);
    Real energy = 0.5*speed2;
    Real bmag = relativistic_mode ? ParticleBmag(pr, p) :
                sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real vpar = ParticleParallelVelocity(pr, p, bmag);
    Real vperp = ParticlePerpendicularVelocity(pr, p, bmag);
    if (quantity_ == PSPEC_MU &&
        (!particles::relativistic::IsFinite(speed) ||
         !particles::relativistic::IsFinite(bmag) ||
         !particles::relativistic::IsFinite(vpar) ||
         (relativistic_mode && (speed <= 0.0 || bmag <= 0.0)))) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }
    if (quantity_ == PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY &&
        (!particles::relativistic::IsFinite(bmag) ||
         !particles::relativistic::IsFinite(vperp) ||
         (relativistic_mode && bmag <= 0.0))) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }
    Real value = speed;
    if (quantity_ == PSPEC_LEGACY_ENERGY) {
      value = energy;
    } else if (quantity_ == PSPEC_LEGACY_LOGE) {
      Real safe_energy = (energy > log_floor) ? energy : log_floor;
      value = log10(safe_energy);
    } else if (quantity_ == PSPEC_MU) {
      value = ParticlePitchAngleCosine(pr, p, speed, bmag);
    } else if (quantity_ == PSPEC_VELOCITY_MAGNETIC_MOMENT_PROXY) {
      value = ParticleVelocityMagneticMomentProxy(pr, p, bmag);
    } else if (quantity_ == PSPEC_WMAG) {
      value = particles::relativistic::ScaledNorm3(ParticleW(pr, p));
      if (!particles::relativistic::IsFinite(value)) {
        Kokkos::atomic_max(&relativistic_failure(0), 1);
        return;
      }
    } else if (quantity_ == PSPEC_KINETIC_ENERGY_MODEL ||
               quantity_ == PSPEC_LOG10_KINETIC_ENERGY_MODEL) {
      auto status = particles::relativistic::KineticEnergyFromW(
          ParticleW(pr, p), c_model, value);
      if (status != particles::relativistic::StateStatus::success) {
        Kokkos::atomic_max(&relativistic_failure(0), 1 + static_cast<int>(status));
        return;
      }
      if (quantity_ == PSPEC_LOG10_KINETIC_ENERGY_MODEL) {
        Real safe_energy = (value > log_floor) ? value : log_floor;
        value = log10(safe_energy);
      }
    }
    if (!particles::relativistic::IsFinite(value)) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }

    int ip = HistogramBin(value, vmin_, vmax_, nbin_);
    Kokkos::atomic_add(&local_histogram(spec,ip),1);
  });
  auto relativistic_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), relativistic_failure);
  if (relativistic_failure_h(0) != 0) {
    FatalParticleDiagnostic("relativistic_hc pspec derived quantity evaluation failed");
  }
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
        "# metadata: quantity=%s nspecies=%d nbin=%d vmin=%.17e vmax=%.17e "
        "reduce=%d\n",
        quantity_name.c_str(), pm->pmb_pack->ppart->nspecies, nbin, vmin, vmax,
        reduce_histogram ? 1 : 0);
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (IsRelativisticParticleOutput(pp)) {
      WriteRelativisticParticleMetadata(pfile, pp);
      std::fprintf(pfile,
          "# metadata: quantity_units=%s histogram_units=particle_count\n",
          ParticleSpectrumQuantityUnits(spectrum_quantity));
      if (spectrum_quantity == PSPEC_LOG10_KINETIC_ENERGY_MODEL) {
        std::fprintf(pfile,
            "# metadata: log10_floor=%.17e "
            "log10_floor_units=code_velocity_squared\n",
            kParticleSpectrumLogFloor);
      }
    }
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
  PSPEC2_MU_SPEED=0,
  PSPEC2_MU_LEGACY_ENERGY=1,
  PSPEC2_VPAR_VPERP=2,
  PSPEC2_MU_WMAG=3,
  PSPEC2_MU_KINETIC_ENERGY_MODEL=4
};

const char *ParticleJointSpectrumAxis1Units(const int quantity) {
  return (quantity == PSPEC2_VPAR_VPERP) ? "code_velocity" : "dimensionless";
}

const char *ParticleJointSpectrumAxis2Units(const int quantity) {
  return (quantity == PSPEC2_MU_KINETIC_ENERGY_MODEL) ?
      "code_velocity_squared" : "code_velocity";
}

} // namespace

ParticleJointSpectrumOutput::ParticleJointSpectrumOutput(ParameterInput *pin, Mesh *pm,
                                                         OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  int nspec = pm->pmb_pack->ppart->nspecies;
  bool relativistic_mode =
      IsRelativisticParticleOutput(pm->pmb_pack->ppart);
  mkdir("pspec2",0775);

  if (relativistic_mode &&
      !pin->DoesParameterExist(op.block_name,"quantity")) {
    FatalParticleDiagnostic("relativistic_hc pspec2 output block '" +
        op.block_name + "' requires explicit quantity");
  }
  quantity_name = relativistic_mode ?
      pin->GetString(op.block_name,"quantity") :
      pin->GetOrAddString(op.block_name,"quantity","mu_p");
  if (relativistic_mode) {
    if (quantity_name == "mu_speed") {
      spectrum_quantity = PSPEC2_MU_SPEED;
    } else if (quantity_name == "mu_wmag") {
      spectrum_quantity = PSPEC2_MU_WMAG;
    } else if (quantity_name == "mu_kinetic_energy_model") {
      spectrum_quantity = PSPEC2_MU_KINETIC_ENERGY_MODEL;
    } else if (quantity_name == "vpar_vperp") {
      spectrum_quantity = PSPEC2_VPAR_VPERP;
    } else if (quantity_name == "mu_p" || quantity_name == "mu,p" ||
               quantity_name == "f_mu_p" || quantity_name == "mu_E" ||
               quantity_name == "mu,e" || quantity_name == "f_mu_E" ||
               quantity_name == "v_parallel_v_perp") {
      FatalParticleDiagnostic("relativistic_hc pspec2 output block '" +
          op.block_name + "' rejects ambiguous quantity alias '" +
          quantity_name + "'");
    } else {
      FatalParticleDiagnostic("Unknown relativistic_hc pspec2 quantity '" +
          quantity_name + "' in output block '" + op.block_name + "'");
    }
  } else if (quantity_name == "mu_p" || quantity_name == "mu,p" ||
             quantity_name == "f_mu_p") {
    spectrum_quantity = PSPEC2_MU_SPEED;
    quantity_name = "mu_p";
  } else if (quantity_name == "mu_speed") {
    spectrum_quantity = PSPEC2_MU_SPEED;
  } else if (quantity_name == "mu_E" || quantity_name == "mu,e" ||
             quantity_name == "f_mu_E") {
    spectrum_quantity = PSPEC2_MU_LEGACY_ENERGY;
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
  if (!std::isfinite(vmin1) || !std::isfinite(vmax1) ||
      !std::isfinite(vmin2) || !std::isfinite(vmax2) ||
      !std::isfinite(vmax1 - vmin1) || !std::isfinite(vmax2 - vmin2) ||
      vmax1 <= vmin1 || vmax2 <= vmin2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "pspec2 output block '" << op.block_name
              << "' requires finite non-overflowing vmax1 > vmin1 and "
              << "vmax2 > vmin2" << std::endl;
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
  bool relativistic_mode = IsRelativisticParticleOutput(pp);
  Real c_model = pp->c_model;

  DvceArray3D<int> local_histogram("local_particle_joint_spectrum",
                                   nspecies, nbin1, nbin2);
  DvceArray1D<int> relativistic_failure("local_particle_joint_spectrum_failure", 1);
  Kokkos::deep_copy(local_histogram, 0);
  Kokkos::deep_copy(relativistic_failure, 0);
  par_for("part_joint_spectrum_hist",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    int spec = pi(PSP,p);
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > nspecies - 1) ? nspecies - 1 : spec;

    Real speed2 = SQR(pr(IPVX,p)) + SQR(pr(IPVY,p)) + SQR(pr(IPVZ,p));
    Real speed = relativistic_mode ? ParticlePhysicalSpeed(pr, p) : sqrt(speed2);
    Real energy = 0.5*speed2;
    Real bmag = relativistic_mode ? ParticleBmag(pr, p) :
                sqrt(SQR(pr(IPBX,p)) + SQR(pr(IPBY,p)) + SQR(pr(IPBZ,p)));
    Real vpar = ParticleParallelVelocity(pr, p, bmag);
    Real vperp = ParticlePerpendicularVelocity(pr, p, bmag);
    if (!particles::relativistic::IsFinite(speed) ||
        !particles::relativistic::IsFinite(bmag) ||
        !particles::relativistic::IsFinite(vpar) ||
        (relativistic_mode && (speed <= 0.0 || bmag <= 0.0))) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }
    Real mu = ParticlePitchAngleCosine(pr, p, speed, bmag);

    Real value1 = mu;
    Real value2 = speed;
    if (quantity_ == PSPEC2_MU_LEGACY_ENERGY) {
      if (!particles::relativistic::IsFinite(speed2)) {
        Kokkos::atomic_max(&relativistic_failure(0), 1);
        return;
      }
      value2 = energy;
    } else if (quantity_ == PSPEC2_MU_WMAG) {
      value2 = particles::relativistic::ScaledNorm3(ParticleW(pr, p));
      if (!particles::relativistic::IsFinite(value2)) {
        Kokkos::atomic_max(&relativistic_failure(0), 1);
        return;
      }
    } else if (quantity_ == PSPEC2_MU_KINETIC_ENERGY_MODEL) {
      auto status = particles::relativistic::KineticEnergyFromW(
          ParticleW(pr, p), c_model, value2);
      if (status != particles::relativistic::StateStatus::success) {
        Kokkos::atomic_max(&relativistic_failure(0), 1 + static_cast<int>(status));
        return;
      }
    } else if (quantity_ == PSPEC2_VPAR_VPERP) {
      value1 = vpar;
      value2 = vperp;
    }
    if (!particles::relativistic::IsFinite(value1) ||
        !particles::relativistic::IsFinite(value2)) {
      Kokkos::atomic_max(&relativistic_failure(0), 1);
      return;
    }

    int i1 = HistogramBin(value1, vmin1_, vmax1_, nbin1_);
    int i2 = HistogramBin(value2, vmin2_, vmax2_, nbin2_);
    Kokkos::atomic_add(&local_histogram(spec,i1,i2),1);
  });
  auto relativistic_failure_h =
      Kokkos::create_mirror_view_and_copy(HostMemSpace(), relativistic_failure);
  if (relativistic_failure_h(0) != 0) {
    FatalParticleDiagnostic("relativistic_hc pspec2 derived quantity evaluation failed");
  }
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
        "# metadata: quantity=%s nspecies=%d nbin1=%d nbin2=%d "
        "vmin1=%.17e vmax1=%.17e vmin2=%.17e vmax2=%.17e reduce=%d\n",
        quantity_name.c_str(), pm->pmb_pack->ppart->nspecies, nbin1, nbin2,
        vmin1, vmax1, vmin2, vmax2, reduce_histogram ? 1 : 0);
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (IsRelativisticParticleOutput(pp)) {
      WriteRelativisticParticleMetadata(pfile, pp);
      std::fprintf(pfile,
          "# metadata: axis1_units=%s axis2_units=%s histogram_units=particle_count\n",
          ParticleJointSpectrumAxis1Units(spectrum_quantity),
          ParticleJointSpectrumAxis2Units(spectrum_quantity));
    }
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
  PSAMP_LEGACY_ENERGY=1,
  PSAMP_BMAG=2,
  PSAMP_MU=3,
  PSAMP_VPAR=4,
  PSAMP_VPERP=5,
  PSAMP_VELOCITY_MAGNETIC_MOMENT_PROXY=6,
  PSAMP_WMAG=7,
  PSAMP_GAMMA=8,
  PSAMP_KINETIC_ENERGY_MODEL=9,
  PSAMP_R_LARMOR_OVER_DX_MIN=10
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
                                            const std::string &block_name,
                                            const bool relativistic_mode) {
  std::string field = NormalizeParticleSampleField(name);
  if (relativistic_mode &&
      (field == "p" || field == "energy" || field == "e" ||
       field == "mass" || field == "m" || field == "magnetic_moment")) {
    FatalParticleDiagnostic("relativistic_hc psamp output block '" + block_name +
        "' rejects ambiguous field alias '" + name + "'");
  }
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
    return ParticleSampleField("energy", PSAMP_DERIVED_FIELD, PSAMP_LEGACY_ENERGY);
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
                               PSAMP_VELOCITY_MAGNETIC_MOMENT_PROXY);
  } else if (field == "velocity_magnetic_moment_proxy") {
    return ParticleSampleField("velocity_magnetic_moment_proxy", PSAMP_DERIVED_FIELD,
                               PSAMP_VELOCITY_MAGNETIC_MOMENT_PROXY);
  } else if (relativistic_mode && (field == "wx" || field == "w1")) {
    return ParticleSampleField("wx", PSAMP_REAL_FIELD, IPWX);
  } else if (relativistic_mode && (field == "wy" || field == "w2")) {
    return ParticleSampleField("wy", PSAMP_REAL_FIELD, IPWY);
  } else if (relativistic_mode && (field == "wz" || field == "w3")) {
    return ParticleSampleField("wz", PSAMP_REAL_FIELD, IPWZ);
  } else if (relativistic_mode && (field == "cex" || field == "ce1")) {
    return ParticleSampleField("cex", PSAMP_REAL_FIELD, IPCEX);
  } else if (relativistic_mode && (field == "cey" || field == "ce2")) {
    return ParticleSampleField("cey", PSAMP_REAL_FIELD, IPCEY);
  } else if (relativistic_mode && (field == "cez" || field == "ce3")) {
    return ParticleSampleField("cez", PSAMP_REAL_FIELD, IPCEZ);
  } else if (relativistic_mode && field == "work") {
    return ParticleSampleField("work", PSAMP_REAL_FIELD, IPWORK);
  } else if (relativistic_mode && (field == "alpha_s" || field == "alpha")) {
    return ParticleSampleField("alpha_s", PSAMP_REAL_FIELD, IPALPHA);
  } else if (relativistic_mode && field == "wmag") {
    return ParticleSampleField("wmag", PSAMP_DERIVED_FIELD, PSAMP_WMAG);
  } else if (relativistic_mode && field == "gamma") {
    return ParticleSampleField("gamma", PSAMP_DERIVED_FIELD, PSAMP_GAMMA);
  } else if (relativistic_mode && field == "kinetic_energy_model") {
    return ParticleSampleField("kinetic_energy_model", PSAMP_DERIVED_FIELD,
                               PSAMP_KINETIC_ENERGY_MODEL);
  } else if (relativistic_mode && field == "r_larmor_over_dx_min") {
    return ParticleSampleField("r_larmor_over_dx_min", PSAMP_DERIVED_FIELD,
                               PSAMP_R_LARMOR_OVER_DX_MIN);
  }

  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Unknown psamp field '" << name
            << "' in output block '" << block_name << "'" << std::endl;
  std::exit(EXIT_FAILURE);
}

Real PositiveRatio3(const Real numerator, const Real denominator1,
                    const Real denominator2, const Real denominator3) {
  if (!std::isfinite(numerator) || numerator < 0.0 ||
      !std::isfinite(denominator1) || denominator1 <= 0.0 ||
      !std::isfinite(denominator2) || denominator2 <= 0.0 ||
      !std::isfinite(denominator3) || denominator3 <= 0.0) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  if (numerator == 0.0) {return 0.0;}
  int numerator_exponent = 0;
  int denominator1_exponent = 0;
  int denominator2_exponent = 0;
  int denominator3_exponent = 0;
  Real scaled = std::frexp(numerator, &numerator_exponent);
  scaled /= std::frexp(denominator1, &denominator1_exponent);
  scaled /= std::frexp(denominator2, &denominator2_exponent);
  scaled /= std::frexp(denominator3, &denominator3_exponent);
  int exponent = numerator_exponent - denominator1_exponent -
                 denominator2_exponent - denominator3_exponent;
  return std::ldexp(scaled, exponent);
}

Real ParticleSampleDerivedValue(const HostArray2D<Real> &pr, int p, int index,
                                Real c_model, Real dx_min) {
  Real vx = pr(IPVX,p);
  Real vy = pr(IPVY,p);
  Real vz = pr(IPVZ,p);
  Real bx = pr(IPBX,p);
  Real by = pr(IPBY,p);
  Real bz = pr(IPBZ,p);
  Real speed2 = SQR(vx) + SQR(vy) + SQR(vz);
  Real speed = particles::relativistic::ScaledNorm3({vx, vy, vz});
  Real bmag = particles::relativistic::ScaledNorm3({bx, by, bz});
  Real vpar = (bmag > 0.0) ?
      vx*(bx/bmag) + vy*(by/bmag) + vz*(bz/bmag) : 0.0;
  Real vperp = PerpendicularMagnitude({vx, vy, vz}, bx, by, bz, bmag);
  Real invalid = std::numeric_limits<Real>::quiet_NaN();

  if (index == PSAMP_SPEED) {
    return std::isfinite(speed) ? speed : invalid;
  } else if (index == PSAMP_LEGACY_ENERGY) {
    return std::isfinite(speed2) ? 0.5*speed2 : invalid;
  } else if (index == PSAMP_BMAG) {
    return std::isfinite(bmag) ? bmag : invalid;
  } else if (index == PSAMP_MU) {
    if (!std::isfinite(speed) || !std::isfinite(bmag) || !std::isfinite(vpar) ||
        speed <= 0.0 || bmag <= 0.0) {
      return invalid;
    }
    return (bmag > 0.0 && speed > 0.0) ? vpar/speed : 0.0;
  } else if (index == PSAMP_VPAR) {
    if (!std::isfinite(speed) || !std::isfinite(bmag) || !std::isfinite(vpar) ||
        bmag <= 0.0) {
      return invalid;
    }
    return vpar;
  } else if (index == PSAMP_VPERP) {
    if (!std::isfinite(bmag) || !std::isfinite(vperp) || bmag <= 0.0) {
      return invalid;
    }
    return vperp;
  } else if (index == PSAMP_VELOCITY_MAGNETIC_MOMENT_PROXY) {
    if (!std::isfinite(bmag) || !std::isfinite(vperp) || bmag <= 0.0) {
      return invalid;
    }
    return PositiveSquareRatio(vperp, bmag);
  } else if (index == PSAMP_WMAG) {
    return particles::relativistic::ScaledNorm3(ParticleW(pr, p));
  } else if (index == PSAMP_GAMMA) {
    Real gamma = 0.0;
    auto status = particles::relativistic::GammaFromW(ParticleW(pr, p), c_model, gamma);
    return (status == particles::relativistic::StateStatus::success) ?
        gamma : std::numeric_limits<Real>::quiet_NaN();
  } else if (index == PSAMP_KINETIC_ENERGY_MODEL) {
    Real energy = 0.0;
    auto status =
        particles::relativistic::KineticEnergyFromW(ParticleW(pr, p), c_model, energy);
    return (status == particles::relativistic::StateStatus::success) ?
        energy : std::numeric_limits<Real>::quiet_NaN();
  } else if (index == PSAMP_R_LARMOR_OVER_DX_MIN) {
    Real alpha = std::abs(pr(IPALPHA,p));
    Real wperp = PerpendicularMagnitude(ParticleW(pr, p), bx, by, bz, bmag);
    if (!std::isfinite(bmag) || bmag <= 0.0 ||
        !std::isfinite(alpha) || alpha <= 0.0 ||
        !std::isfinite(dx_min) || dx_min <= 0.0 ||
        !std::isfinite(wperp)) {
      return invalid;
    }
    return PositiveRatio3(wperp, bmag, alpha, dx_min);
  }
  return 0.0;
}

Real ParticleSampleDxMin(Mesh *pm, int pgid) {
  int m = pm->FindMeshBlockIndex(pgid);
  if (m < 0) {
    FatalParticleDiagnostic("psamp could not map local particle PGID to a MeshBlock");
  }
  const auto &size = pm->pmb_pack->pmb->mb_size.h_view(m);
  Real dx_min = std::min(size.dx1, std::min(size.dx2, size.dx3));
  if (!std::isfinite(dx_min) || dx_min <= 0.0) {
    FatalParticleDiagnostic("psamp MeshBlock cell widths must be finite and positive");
  }
  return dx_min;
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
  bool relativistic_mode =
      IsRelativisticParticleOutput(pm->pmb_pack->ppart);
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

  if (relativistic_mode &&
      !pin->DoesParameterExist(op.block_name,"fields")) {
    FatalParticleDiagnostic("relativistic_hc psamp output block '" +
        op.block_name + "' requires explicit fields");
  }
  std::string field_text = relativistic_mode ?
      pin->GetString(op.block_name,"fields") :
      pin->GetOrAddString(op.block_name, "fields",
          "x,y,z,vx,vy,vz,bx,by,bz,dx,dy,dz,dpar,mass");
  std::stringstream field_stream(field_text);
  std::string field_name;
  std::set<std::string> canonical_fields;
  while (std::getline(field_stream, field_name, ',')) {
    std::string normalized = NormalizeParticleSampleField(field_name);
    if (normalized.empty()) {
      continue;
    }
    if (normalized == "pgid" || normalized == "gid" ||
        normalized == "tag" || normalized == "species") {
      continue;
    }
    ParticleSampleField sample_field =
        MakeParticleSampleField(field_name, op.block_name, relativistic_mode);
    if (relativistic_mode &&
        !canonical_fields.insert(sample_field.label).second) {
      FatalParticleDiagnostic("relativistic_hc psamp output block '" +
          op.block_name + "' contains duplicate canonical field '" +
          sample_field.label + "'");
    }
    sample_fields.push_back(sample_field);
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
  particles::Particles *pp = pm->pmb_pack->ppart;
  if (IsRelativisticParticleOutput(pp)) {
    WriteRelativisticParticleMetadata(pfile, pp);
    std::fprintf(pfile,
        "# metadata: position_units=code_length velocity_units=code_velocity "
        "w_units=code_velocity cE_units=code_velocity_times_code_B "
        "gamma_units=dimensionless "
        "kinetic_energy_model_units=code_velocity_squared "
        "work_units=code_velocity_squared "
        "alpha_s_units=inverse_code_time_per_code_B "
        "r_larmor_over_dx_min_units=dimensionless\n");
  }
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
    std::vector<Real> sample_values;
    sample_values.reserve(sample_fields.size());
    for (const auto &field : sample_fields) {
      Real dx_min = (field.kind == PSAMP_DERIVED_FIELD &&
                     field.index == PSAMP_R_LARMOR_OVER_DX_MIN) ?
                    ParticleSampleDxMin(pm, pgid) : 0.0;
      Real value = (field.kind == PSAMP_REAL_FIELD)
                 ? outpart_rdata(field.index,p)
                 : ParticleSampleDerivedValue(outpart_rdata, p, field.index,
                                              pp->c_model, dx_min);
      if (IsRelativisticParticleOutput(pp) && !std::isfinite(value)) {
        FatalParticleDiagnostic("relativistic_hc psamp derived field evaluation failed");
      }
      sample_values.push_back(value);
    }
    std::fprintf(pfile, "%d %d %d %d", global_variable::my_rank, pgid, tag, spec);
    for (const Real value : sample_values) {
      std::fprintf(pfile, " %.17e", value);
    }
    std::fprintf(pfile, "\n");
  }
  std::fclose(pfile);

  out_params.last_time = (out_params.last_time < 0.0 || out_params.last_time == pm->time)
                       ? pm->time : out_params.last_time + out_params.dt;
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
