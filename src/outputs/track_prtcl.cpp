//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file track_prtcl.cpp
//! \brief writes data for tracked particles in unformatted binary

#include <algorithm>
#include <chrono>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>  // mkdir
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op),
  track_cache_probe_initialized(false),
  last_output_cycle(-1),
  track_cycles_buffered(0) {
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("trk",0775);
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  track_per_species = pin->GetOrAddBoolean(op.block_name,"track_per_species",true);
  track_cache_probe = pin->GetOrAddBoolean(op.block_name,"cache_probe",false);
  track_buffer_size = pin->GetOrAddInteger(op.block_name,"buffer_size",0);
  track_ncycle_buffer = pin->GetOrAddInteger(op.block_name,"ncycle",1);
  track_ncycle_buffer = std::max(track_ncycle_buffer, 1);
  int nspecies = pm->pmb_pack->ppart->nspecies;
  ntrack_total = track_per_species ? ntrack*nspecies : ntrack;
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack_total;
  if (track_buffer_size > 0 && global_variable::my_rank == 0) {
    track_buffer.reserve(static_cast<std::size_t>(track_buffer_size));
  }
  if (track_cache_probe) {
    Kokkos::realloc(track_cache_indices, ntrack_total);
    Kokkos::deep_copy(track_cache_indices, -1);
  }
  out_params.last_time = pm->time;
}

TrackedParticleOutput::~TrackedParticleOutput() {
  if (global_variable::my_rank == 0) {
    std::string fname("trk/");
    fname.append(out_params.file_basename);
    fname.append(".trk");
    FlushTrackBuffer(fname);
  }
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  double probe_ms = 0.0;
  int probe_prev_local = 0;
  int probe_hits = 0;
  int probe_stale_oob = 0;
  int probe_stale_mismatch = 0;
  bool have_probe_stats = false;

  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank);
  int npart = pm->nprtcl_thisrank;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  int ntrack_local = ntrack;
  int ntrack_total_local = ntrack_total;
  bool track_per_species_local = track_per_species;
  int nspecies = pm->pmb_pack->ppart->nspecies;
  int species_offset = (nspecies > 0) ? pm->nprtcl_total/nspecies : 0;
  bool cache_probe_local = track_cache_probe;
  auto cache_indices = track_cache_indices;
  if (cache_probe_local) {
    if (track_cache_probe_initialized) {
      auto probe_start = std::chrono::steady_clock::now();
      Kokkos::parallel_reduce("trk_cache_probe",
          Kokkos::RangePolicy<>(DevExeSpace(), 0, ntrack_total),
          KOKKOS_LAMBDA(const int t, int &prev, int &hits,
                        int &stale_oob, int &stale_mismatch) {
        int p = cache_indices(t);
        if (p < 0) {return;}
        prev += 1;
        if (p >= npart) {
          stale_oob += 1;
          return;
        }
        int tag = pi(PTAG,p);
        int spec = pi(PSP,p);
        int track_tag = tag;
        if (track_per_species_local) {
          track_tag = tag - spec*species_offset;
        }
        int output_tag = track_per_species_local ? spec*ntrack_local + track_tag : tag;
        if (track_tag >= 0 && track_tag < ntrack_local && output_tag == t) {
          hits += 1;
        } else {
          stale_mismatch += 1;
        }
      }, Kokkos::Sum<int>(probe_prev_local),
         Kokkos::Sum<int>(probe_hits),
         Kokkos::Sum<int>(probe_stale_oob),
         Kokkos::Sum<int>(probe_stale_mismatch));
      Kokkos::fence();
      auto probe_end = std::chrono::steady_clock::now();
      probe_ms = std::chrono::duration<double,std::milli>(probe_end - probe_start).count();
      have_probe_stats = true;
    }
    Kokkos::deep_copy(cache_indices, -1);
  }
  auto scan_start = std::chrono::steady_clock::now();
  DvceArray1D<int> counter("tracked_particle_counter",1);
  Kokkos::deep_copy(counter, 0);
  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int tag = pi(PTAG,p);
    int spec = pi(PSP,p);
    int track_tag = tag;
    if (track_per_species_local) {
      track_tag = tag - spec*species_offset;
    }
    if (track_tag >= 0 && track_tag < ntrack_local) {
      int output_tag = track_per_species_local ? spec*ntrack_local + track_tag : tag;
      if (output_tag < 0 || output_tag >= ntrack_total_local) {
        return;
      }
      if (cache_probe_local) {
        cache_indices(output_tag) = p;
      }
      int index = Kokkos::atomic_fetch_add(&counter(0),1);
      tracked_prtcl.d_view(index).tag = output_tag;
      tracked_prtcl.d_view(index).x   = pr(IPX,p);
      tracked_prtcl.d_view(index).y   = pr(IPY,p);
      tracked_prtcl.d_view(index).z   = pr(IPZ,p);
      tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
      tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
      tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
    }
  });
  auto counter_host = Kokkos::create_mirror_view_and_copy(HostMemSpace(), counter);
  npout = counter_host(0);
  if (npout > ntrack_thisrank) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle output found " << npout
              << " local records, but the allocation only permits "
              << ntrack_thisrank << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // share number of tracked particles to be output across all ranks
  npout_eachrank[global_variable::my_rank] = npout;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&npout, 1, MPI_INT, npout_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#endif
  tracked_prtcl.resize(npout);
  // sync tracked particle device array with host
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

  // copy host view into host outpart array
  Kokkos::realloc(outpart, npout);
  Kokkos::deep_copy(outpart, tracked_prtcl.h_view);
  if (cache_probe_local) {
    track_cache_probe_initialized = true;
    auto scan_end = std::chrono::steady_clock::now();
    double scan_ms = std::chrono::duration<double,std::milli>(scan_end - scan_start).count();
    LogTrackCacheProbe(pm, have_probe_stats, probe_prev_local, probe_hits,
                       probe_stale_oob, probe_stale_mismatch, probe_ms, scan_ms);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (last_output_cycle == pm->ncycle) {
    if (global_variable::my_rank == 0) {
      std::string fname("trk/");
      fname.append(out_params.file_basename);
      fname.append(".trk");
      FlushTrackBuffer(fname);
    }
    return;
  }

  // create filename: "trk/file_basename".trk
  std::string fname;
  fname.assign("trk/");
  fname.append(out_params.file_basename);
  fname.append(".trk");

  int nout_total = 0;
  for (int count : npout_eachrank) {
    nout_total += count;
  }
  if (nout_total != ntrack_total) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle output found " << nout_total
              << " global records, but expected " << ntrack_total << std::endl;
    std::exit(EXIT_FAILURE);
  }

#if MPI_PARALLEL_ENABLED
  std::vector<int> recv_counts;
  std::vector<int> displs;
  std::vector<TrackedParticleData> gathered;
  if (global_variable::my_rank == 0) {
    recv_counts.resize(global_variable::nranks);
    displs.resize(global_variable::nranks);
    int offset = 0;
    for (int n=0; n<global_variable::nranks; ++n) {
      recv_counts[n] = npout_eachrank[n]*static_cast<int>(sizeof(TrackedParticleData));
      displs[n] = offset;
      offset += recv_counts[n];
    }
    gathered.resize(nout_total);
  }
  int send_count = npout*static_cast<int>(sizeof(TrackedParticleData));
  void *send_buffer = (npout > 0) ? static_cast<void*>(outpart.data()) : nullptr;
  void *recv_buffer = gathered.empty() ? nullptr : static_cast<void*>(gathered.data());
  MPI_Gatherv(send_buffer, send_count, MPI_BYTE,
              recv_buffer, recv_counts.data(), displs.data(), MPI_BYTE,
              0, MPI_COMM_WORLD);
#else
  std::vector<TrackedParticleData> gathered(nout_total);
  for (int p=0; p<npout; ++p) {
    gathered[p] = outpart(p);
  }
#endif

  if (global_variable::my_rank == 0) {
    std::vector<float> data(6*ntrack_total,
                            std::numeric_limits<float>::quiet_NaN());
    std::vector<unsigned char> seen(ntrack_total, 0);
    for (const auto &record : gathered) {
      if (record.tag < 0 || record.tag >= ntrack_total) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Tracked particle tag " << record.tag
                  << " is outside [0," << ntrack_total << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (seen[record.tag] != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Duplicate tracked particle tag "
                  << record.tag << std::endl;
        std::exit(EXIT_FAILURE);
      }
      seen[record.tag] = 1;
      int base = 6*record.tag;
      data[base    ] = static_cast<float>(record.x);
      data[base + 1] = static_cast<float>(record.y);
      data[base + 2] = static_cast<float>(record.z);
      data[base + 3] = static_cast<float>(record.vx);
      data[base + 4] = static_cast<float>(record.vy);
      data[base + 5] = static_cast<float>(record.vz);
    }

    std::stringstream msg;
    msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle
        << "  ntracked_prtcls=" << ntrack_total
        << "  ntrack_per_species=" << ntrack
        << "  track_per_species=" << (track_per_species ? 1 : 0) << std::endl;
    std::string header = msg.str();
    header.append(" \n");
    track_buffer.insert(track_buffer.end(), header.begin(), header.end());
    const char *payload = reinterpret_cast<const char*>(data.data());
    track_buffer.insert(track_buffer.end(), payload,
                        payload + data.size()*sizeof(float));
    track_cycles_buffered += 1;
  }

  // increment counters
  float time_32 = static_cast<float>(pm->time);
  float next_32 = static_cast<float>(out_params.last_time + out_params.dt);
  bool final_forced_output = (out_params.dt > 0.0 && out_params.last_time >= 0.0 &&
                              time_32 < next_32);
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  last_output_cycle = pm->ncycle;
  if (global_variable::my_rank == 0 &&
      (final_forced_output ||
       track_buffer_size <= 0 ||
       static_cast<int>(track_buffer.size()) >= track_buffer_size ||
       track_cycles_buffered >= track_ncycle_buffer)) {
    FlushTrackBuffer(fname);
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}

void TrackedParticleOutput::FlushTrackBuffer(const std::string &fname) {
  if (track_buffer.empty()) {return;}
  FILE *pfile;
  if ((pfile = std::fopen(fname.c_str(),"ab")) == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (std::fwrite(track_buffer.data(), 1, track_buffer.size(), pfile) !=
      track_buffer.size()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle buffer not written correctly"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::fclose(pfile);
  track_buffer.clear();
  track_cycles_buffered = 0;
}

void TrackedParticleOutput::LogTrackCacheProbe(Mesh *pm, bool have_probe_stats,
                                               int prev_local, int hits,
                                               int stale_oob, int stale_mismatch,
                                               double probe_ms, double scan_ms) {
  long long local_stats[5] = {
    static_cast<long long>(prev_local),
    static_cast<long long>(hits),
    static_cast<long long>(stale_oob),
    static_cast<long long>(stale_mismatch),
    static_cast<long long>(npout)};
  long long global_stats[5] = {0, 0, 0, 0, 0};
  double local_times[2] = {probe_ms, scan_ms};
  double global_times[2] = {0.0, 0.0};
#if MPI_PARALLEL_ENABLED
  MPI_Reduce(local_stats, global_stats, 5, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_times, global_times, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else
  for (int n=0; n<5; ++n) {global_stats[n] = local_stats[n];}
  for (int n=0; n<2; ++n) {global_times[n] = local_times[n];}
#endif
  if (global_variable::my_rank != 0) {return;}

  long long prev_global = global_stats[0];
  long long hits_global = global_stats[1];
  long long stale_global = global_stats[2] + global_stats[3];
  long long current_global = global_stats[4];
  double hit_frac = (prev_global > 0) ?
      static_cast<double>(hits_global)/static_cast<double>(prev_global) : 0.0;
  double stale_frac = (prev_global > 0) ?
      static_cast<double>(stale_global)/static_cast<double>(prev_global) : 0.0;
  double arrival_frac = (current_global > 0) ?
      static_cast<double>(current_global - hits_global)/static_cast<double>(current_global) : 0.0;

  std::cout << "trk_cache_probe: output=" << (have_probe_stats ? "check" : "init")
            << " time=" << std::setprecision(12) << pm->time
            << " cycle=" << pm->ncycle
            << " tracked=" << ntrack_total
            << " prev=" << prev_global
            << " current=" << current_global
            << " hits=" << hits_global
            << " stale=" << stale_global
            << " stale_oob=" << global_stats[2]
            << " stale_mismatch=" << global_stats[3]
            << " hit_frac=" << std::setprecision(6) << hit_frac
            << " stale_frac=" << stale_frac
            << " arrival_or_moved_frac=" << arrival_frac
            << " probe_ms_max=" << global_times[0]
            << " scan_ms_max=" << global_times[1]
            << std::endl;
}
