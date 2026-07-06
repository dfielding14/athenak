//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file track_prtcl.cpp
//! \brief writes rich tracked-particle records in unformatted binary

#include <algorithm>
#include <chrono>
#include <cstdio>      // fwrite(), fclose(), fopen(), snprintf()
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>  // mkdir
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "outputs.hpp"
#include "particles/particles.hpp"

namespace {
constexpr int kTrackRecordFields = 18;
constexpr char kCompactTrackFileMagic[8] = {'A','K','T','R','K','2','F','\0'};
constexpr char kCompactTrackFrameMagic[8] = {'A','K','T','R','K','2','R','\0'};
constexpr std::uint16_t kCompactTrackVersion = 2;
constexpr std::uint16_t kCompactTrackPrologueBytes = 68;
constexpr std::uint16_t kCompactTrackFrameBytes = 40;
constexpr const char *kTrackFieldList =
    "tag,time,x,y,z,vx,vy,vz,bx,by,bz,k1,k2,k3,db1,db2,db3,jmag";

template <typename T>
void AppendScalar(std::vector<char> &buffer, const T &value) {
  const char *data = reinterpret_cast<const char*>(&value);
  buffer.insert(buffer.end(), data, data + sizeof(T));
}

void AppendBytes(std::vector<char> &buffer, const void *data, std::size_t bytes) {
  const char *ptr = reinterpret_cast<const char*>(data);
  buffer.insert(buffer.end(), ptr, ptr + bytes);
}

std::uint32_t CompactLayoutCode(FileShardMode mode) {
  switch (mode) {
    case FileShardMode::shared:
      return 0;
    case FileShardMode::per_node:
      return 1;
    case FileShardMode::per_rank:
      return 2;
    default:
      return 0;
  }
}

KOKKOS_INLINE_FUNCTION
Real SafeBmag(Real bx, Real by, Real bz) {
  return Kokkos::sqrt(bx*bx + by*by + bz*bz);
}

KOKKOS_INLINE_FUNCTION
void SafeBhat(Real bx, Real by, Real bz, Real &b1, Real &b2, Real &b3) {
  Real bmag = SafeBmag(bx, by, bz);
  if (bmag > 0.0) {
    Real inv_bmag = 1.0/bmag;
    b1 = bx*inv_bmag;
    b2 = by*inv_bmag;
    b3 = bz*inv_bmag;
  } else {
    b1 = 0.0;
    b2 = 0.0;
    b3 = 0.0;
  }
}

template <typename T>
bool WriteBytes(FILE *pfile, const T *data, std::size_t count,
                const char *error_context) {
  if (count == 0) {return true;}
  if (std::fwrite(data, sizeof(T), count, pfile) != count) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << error_context << " not written correctly"
              << std::endl;
    return false;
  }
  return true;
}
} // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op),
  track_cache_probe_initialized(false),
  last_output_cycle(-1),
  track_cycles_buffered(0) {
  if (pm->pmb_pack == nullptr || pm->pmb_pack->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle output requires particles"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pm->pmb_pack->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Rich tracked particle output requires MHD fields"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  mkdir("trk",0775);
  if (out_params.file_shard_mode != FileShardMode::shared && IsShardWriter(
          out_params.file_shard_mode)) {
    std::string dirname("trk/");
    dirname.append(ShardDirectoryName(out_params.file_shard_mode));
    mkdir(dirname.c_str(),0775);
  }

  npout_eachrank.resize(global_variable::nranks, 0);
  ntrack = pin->GetInteger(op.block_name,"nparticles");
  track_per_species = pin->GetOrAddBoolean(op.block_name,"track_per_species",true);
  track_cache_probe = pin->GetOrAddBoolean(op.block_name,"cache_probe",false);
  track_validate_global_tags = pin->GetOrAddBoolean(op.block_name,
                                                    "validate_global_tags",false);
  std::string header_format = pin->GetOrAddString(op.block_name, "trk_header_format",
                                                  "legacy");
  if (header_format.compare("legacy") == 0 || header_format.compare("rich_v1") == 0) {
    track_header_format = TrackHeaderFormat::legacy;
  } else if (header_format.compare("compact") == 0 ||
             header_format.compare("rich_v2") == 0) {
    track_header_format = TrackHeaderFormat::compact;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Unknown trk_header_format = '" << header_format
              << "' in output block '" << op.block_name << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  track_buffer_size = pin->GetOrAddInteger(op.block_name,"buffer_size",0);
  track_ncycle_buffer = pin->GetOrAddInteger(op.block_name,"ncycle",1);
  track_ncycle_buffer = std::max(track_ncycle_buffer, 1);
  int nspecies = pm->pmb_pack->ppart->nspecies;
  ntrack_total = track_per_species ? ntrack*nspecies : ntrack;
  ntrack_thisrank = ntrack_total;
  if (track_buffer_size > 0 && out_params.file_shard_mode != FileShardMode::shared &&
      IsShardWriter(out_params.file_shard_mode)) {
    track_buffer.reserve(static_cast<std::size_t>(track_buffer_size));
  }
  if (track_cache_probe) {
    Kokkos::realloc(track_cache_indices, ntrack_total);
    Kokkos::deep_copy(track_cache_indices, -1);
  }
  out_params.last_time = pm->time;
}

TrackedParticleOutput::~TrackedParticleOutput() {
  if (out_params.file_shard_mode != FileShardMode::shared &&
      IsShardWriter(out_params.file_shard_mode)) {
    FlushTrackBuffer(TrackFilename());
  }
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies rich data for tracked particles on this rank to host outpart array.

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  double probe_ms = 0.0;
  int probe_prev_local = 0;
  int probe_hits = 0;
  int probe_stale_oob = 0;
  int probe_stale_mismatch = 0;
  bool have_probe_stats = false;

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

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int nmb = pm->pmb_pack->nmb_thispack;
  const int gids = pm->pmb_pack->gids;
  auto &mbsize = pm->pmb_pack->pmb->mb_size;
  auto &bcc = pm->pmb_pack->pmhd->bcc0;

  par_for("part_trackout",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
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
      if (index >= ntrack_total_local) {return;}

      tracked_prtcl.d_view(index).tag = output_tag;
      tracked_prtcl.d_view(index).x   = pr(IPX,p);
      tracked_prtcl.d_view(index).y   = pr(IPY,p);
      tracked_prtcl.d_view(index).z   = pr(IPZ,p);
      tracked_prtcl.d_view(index).vx  = pr(IPVX,p);
      tracked_prtcl.d_view(index).vy  = pr(IPVY,p);
      tracked_prtcl.d_view(index).vz  = pr(IPVZ,p);
      tracked_prtcl.d_view(index).Bx  = pr(IPBX,p);
      tracked_prtcl.d_view(index).By  = pr(IPBY,p);
      tracked_prtcl.d_view(index).Bz  = pr(IPBZ,p);

      Real k1 = 0.0;
      Real k2 = 0.0;
      Real k3 = 0.0;
      Real db1 = 0.0;
      Real db2 = 0.0;
      Real db3 = 0.0;
      Real jmag = 0.0;

      int m = pi(PGID,p) - gids;
      if (m >= 0 && m < nmb) {
        int i = static_cast<int>((pr(IPX,p) - mbsize.d_view(m).x1min)/
                                 mbsize.d_view(m).dx1) + is;
        int j = static_cast<int>((pr(IPY,p) - mbsize.d_view(m).x2min)/
                                 mbsize.d_view(m).dx2) + js;
        int k = static_cast<int>((pr(IPZ,p) - mbsize.d_view(m).x3min)/
                                 mbsize.d_view(m).dx3) + ks;
        i = (i < is) ? is : ((i > ie) ? ie : i);
        j = (j < js) ? js : ((j > je) ? je : j);
        k = (k < ks) ? ks : ((k > ke) ? ke : k);

        Real b1c, b2c, b3c;
        SafeBhat(bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i), bcc(m,IBZ,k,j,i),
                 b1c, b2c, b3c);
        const bool have_x1_stencil = (i > 0 && i + 1 < ncells1 &&
                                      mbsize.d_view(m).dx1 > 0.0);
        const bool have_x2_stencil = (nx2 > 1 && j > 0 && j + 1 < ncells2 &&
                                      mbsize.d_view(m).dx2 > 0.0);
        const bool have_x3_stencil = (nx3 > 1 && k > 0 && k + 1 < ncells3 &&
                                      mbsize.d_view(m).dx3 > 0.0);

        Real b1_ip1 = b1c, b2_ip1 = b2c, b3_ip1 = b3c;
        Real b1_im1 = b1c, b2_im1 = b2c, b3_im1 = b3c;
        Real b1_jp1 = b1c, b2_jp1 = b2c, b3_jp1 = b3c;
        Real b1_jm1 = b1c, b2_jm1 = b2c, b3_jm1 = b3c;
        Real b1_kp1 = b1c, b2_kp1 = b2c, b3_kp1 = b3c;
        Real b1_km1 = b1c, b2_km1 = b2c, b3_km1 = b3c;
        if (have_x1_stencil) {
          SafeBhat(bcc(m,IBX,k,j,i+1), bcc(m,IBY,k,j,i+1),
                   bcc(m,IBZ,k,j,i+1), b1_ip1, b2_ip1, b3_ip1);
          SafeBhat(bcc(m,IBX,k,j,i-1), bcc(m,IBY,k,j,i-1),
                   bcc(m,IBZ,k,j,i-1), b1_im1, b2_im1, b3_im1);
        }
        if (have_x2_stencil) {
          SafeBhat(bcc(m,IBX,k,j+1,i), bcc(m,IBY,k,j+1,i),
                   bcc(m,IBZ,k,j+1,i), b1_jp1, b2_jp1, b3_jp1);
          SafeBhat(bcc(m,IBX,k,j-1,i), bcc(m,IBY,k,j-1,i),
                   bcc(m,IBZ,k,j-1,i), b1_jm1, b2_jm1, b3_jm1);
        }
        if (have_x3_stencil) {
          SafeBhat(bcc(m,IBX,k+1,j,i), bcc(m,IBY,k+1,j,i),
                   bcc(m,IBZ,k+1,j,i), b1_kp1, b2_kp1, b3_kp1);
          SafeBhat(bcc(m,IBX,k-1,j,i), bcc(m,IBY,k-1,j,i),
                   bcc(m,IBZ,k-1,j,i), b1_km1, b2_km1, b3_km1);
        }

        Real dbhat1_dx1 = 0.0, dbhat2_dx1 = 0.0, dbhat3_dx1 = 0.0;
        Real dbhat1_dx2 = 0.0, dbhat2_dx2 = 0.0, dbhat3_dx2 = 0.0;
        Real dbhat1_dx3 = 0.0, dbhat2_dx3 = 0.0, dbhat3_dx3 = 0.0;
        if (have_x1_stencil) {
          dbhat1_dx1 = (b1_ip1 - b1_im1)/(2.0*mbsize.d_view(m).dx1);
          dbhat2_dx1 = (b2_ip1 - b2_im1)/(2.0*mbsize.d_view(m).dx1);
          dbhat3_dx1 = (b3_ip1 - b3_im1)/(2.0*mbsize.d_view(m).dx1);
        }
        if (have_x2_stencil) {
          dbhat1_dx2 = (b1_jp1 - b1_jm1)/(2.0*mbsize.d_view(m).dx2);
          dbhat2_dx2 = (b2_jp1 - b2_jm1)/(2.0*mbsize.d_view(m).dx2);
          dbhat3_dx2 = (b3_jp1 - b3_jm1)/(2.0*mbsize.d_view(m).dx2);
        }
        if (have_x3_stencil) {
          dbhat1_dx3 = (b1_kp1 - b1_km1)/(2.0*mbsize.d_view(m).dx3);
          dbhat2_dx3 = (b2_kp1 - b2_km1)/(2.0*mbsize.d_view(m).dx3);
          dbhat3_dx3 = (b3_kp1 - b3_km1)/(2.0*mbsize.d_view(m).dx3);
        }

        k1 = b1c*dbhat1_dx1 + b2c*dbhat1_dx2 + b3c*dbhat1_dx3;
        k2 = b1c*dbhat2_dx1 + b2c*dbhat2_dx2 + b3c*dbhat2_dx3;
        k3 = b1c*dbhat3_dx1 + b2c*dbhat3_dx2 + b3c*dbhat3_dx3;

        Real bmag_c = SafeBmag(bcc(m,IBX,k,j,i), bcc(m,IBY,k,j,i),
                               bcc(m,IBZ,k,j,i));
        Real bmag_ip1 = bmag_c, bmag_im1 = bmag_c;
        Real bmag_jp1 = bmag_c, bmag_jm1 = bmag_c;
        Real bmag_kp1 = bmag_c, bmag_km1 = bmag_c;
        if (have_x1_stencil) {
          bmag_ip1 = SafeBmag(bcc(m,IBX,k,j,i+1), bcc(m,IBY,k,j,i+1),
                              bcc(m,IBZ,k,j,i+1));
          bmag_im1 = SafeBmag(bcc(m,IBX,k,j,i-1), bcc(m,IBY,k,j,i-1),
                              bcc(m,IBZ,k,j,i-1));
          db1 = (bmag_ip1 - bmag_im1)/(2.0*mbsize.d_view(m).dx1);
        }
        if (have_x2_stencil) {
          bmag_jp1 = SafeBmag(bcc(m,IBX,k,j+1,i), bcc(m,IBY,k,j+1,i),
                              bcc(m,IBZ,k,j+1,i));
          bmag_jm1 = SafeBmag(bcc(m,IBX,k,j-1,i), bcc(m,IBY,k,j-1,i),
                              bcc(m,IBZ,k,j-1,i));
          db2 = (bmag_jp1 - bmag_jm1)/(2.0*mbsize.d_view(m).dx2);
        }
        if (have_x3_stencil) {
          bmag_kp1 = SafeBmag(bcc(m,IBX,k+1,j,i), bcc(m,IBY,k+1,j,i),
                              bcc(m,IBZ,k+1,j,i));
          bmag_km1 = SafeBmag(bcc(m,IBX,k-1,j,i), bcc(m,IBY,k-1,j,i),
                              bcc(m,IBZ,k-1,j,i));
          db3 = (bmag_kp1 - bmag_km1)/(2.0*mbsize.d_view(m).dx3);
        }

        Real dBx_dx2 = 0.0, dBx_dx3 = 0.0;
        Real dBy_dx1 = 0.0, dBy_dx3 = 0.0;
        Real dBz_dx1 = 0.0, dBz_dx2 = 0.0;
        if (have_x1_stencil) {
          dBy_dx1 = (bcc(m,IBY,k,j,i+1) - bcc(m,IBY,k,j,i-1))/
                    (2.0*mbsize.d_view(m).dx1);
          dBz_dx1 = (bcc(m,IBZ,k,j,i+1) - bcc(m,IBZ,k,j,i-1))/
                    (2.0*mbsize.d_view(m).dx1);
        }
        if (have_x2_stencil) {
          dBx_dx2 = (bcc(m,IBX,k,j+1,i) - bcc(m,IBX,k,j-1,i))/
                    (2.0*mbsize.d_view(m).dx2);
          dBz_dx2 = (bcc(m,IBZ,k,j+1,i) - bcc(m,IBZ,k,j-1,i))/
                    (2.0*mbsize.d_view(m).dx2);
        }
        if (have_x3_stencil) {
          dBx_dx3 = (bcc(m,IBX,k+1,j,i) - bcc(m,IBX,k-1,j,i))/
                    (2.0*mbsize.d_view(m).dx3);
          dBy_dx3 = (bcc(m,IBY,k+1,j,i) - bcc(m,IBY,k-1,j,i))/
                    (2.0*mbsize.d_view(m).dx3);
        }

        Real j1 = dBz_dx2 - dBy_dx3;
        Real j2 = dBx_dx3 - dBz_dx1;
        Real j3 = dBy_dx1 - dBx_dx2;
        jmag = Kokkos::sqrt(j1*j1 + j2*j2 + j3*j3);
      }

      tracked_prtcl.d_view(index).K1 = k1;
      tracked_prtcl.d_view(index).K2 = k2;
      tracked_prtcl.d_view(index).K3 = k3;
      tracked_prtcl.d_view(index).dB1 = db1;
      tracked_prtcl.d_view(index).dB2 = db2;
      tracked_prtcl.d_view(index).dB3 = db3;
      tracked_prtcl.d_view(index).jmag = jmag;
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

  std::fill(npout_eachrank.begin(), npout_eachrank.end(), 0);
  npout_eachrank[global_variable::my_rank] = npout;
  tracked_prtcl.resize(npout);
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

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
//! \brief Writes rich tracked-particle records.  Production mode writes one file per rank.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (last_output_cycle == pm->ncycle) {
    if (out_params.file_shard_mode != FileShardMode::shared &&
        IsShardWriter(out_params.file_shard_mode)) {
      FlushTrackBuffer(TrackFilename());
    }
    return;
  }

  ValidateTrackedRecords();
  std::vector<float> records = PackLocalTrackRecords(pm);

  if (out_params.file_shard_mode == FileShardMode::per_rank) {
    AppendTrackBuffer(pm, records);
  } else if (out_params.file_shard_mode == FileShardMode::per_node) {
    WriteNodeTrackFrame(pm, records);
  } else {
    WriteSharedTrackFrame(pm, records);
  }

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
  if (out_params.file_shard_mode != FileShardMode::shared &&
      IsShardWriter(out_params.file_shard_mode) &&
      (final_forced_output ||
       track_buffer_size <= 0 ||
       static_cast<int>(track_buffer.size()) >= track_buffer_size ||
       track_cycles_buffered >= track_ncycle_buffer)) {
    FlushTrackBuffer(TrackFilename());
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}

std::string TrackedParticleOutput::TrackFilename() const {
  std::string fname("trk/");
  fname.append(ShardDirectoryName(out_params.file_shard_mode));
  fname.append(out_params.file_basename);
  fname.append(".trk");
  return fname;
}

std::string TrackedParticleOutput::TrackHeader(Mesh *pm, int record_count) const {
  std::stringstream msg;
  msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
      << "  nranks= " << global_variable::nranks
      << "  nnodes= " << global_variable::nnodes
      << "  cycle=" << pm->ncycle
      << "  ntracked_prtcls=" << ntrack_total
      << "  ntrack_per_species=" << ntrack
      << "  track_per_species=" << (track_per_species ? 1 : 0)
      << "  record_count=" << record_count
      << std::endl;
  msg << "# trk_format=rich_v1"
      << "  nfields=" << kTrackRecordFields
      << "  layout=" << ShardDistributionName(out_params.file_shard_mode)
      << "  rank=" << global_variable::my_rank
      << "  node=" << global_variable::node_id
      << "  ranks_per_node=" << global_variable::ranks_per_node
      << std::endl;
  msg << "# fields=" << kTrackFieldList << std::endl;
  msg << " " << std::endl;
  return msg.str();
}

std::vector<char> TrackedParticleOutput::CompactTrackPrologue() const {
  std::vector<char> buffer;
  const std::string fields(kTrackFieldList);
  buffer.reserve(kCompactTrackPrologueBytes + fields.size());
  AppendBytes(buffer, kCompactTrackFileMagic, sizeof(kCompactTrackFileMagic));
  AppendScalar(buffer, kCompactTrackVersion);
  AppendScalar(buffer, kCompactTrackPrologueBytes);
  AppendScalar(buffer, static_cast<std::uint32_t>(kTrackRecordFields));
  AppendScalar(buffer, static_cast<std::int64_t>(ntrack_total));
  AppendScalar(buffer, static_cast<std::int64_t>(ntrack));
  AppendScalar(buffer, static_cast<std::uint32_t>(track_per_species ? 1 : 0));
  AppendScalar(buffer, CompactLayoutCode(out_params.file_shard_mode));
  AppendScalar(buffer, static_cast<std::int32_t>(global_variable::my_rank));
  AppendScalar(buffer, static_cast<std::int32_t>(global_variable::node_id));
  AppendScalar(buffer, static_cast<std::int32_t>(global_variable::nranks));
  AppendScalar(buffer, static_cast<std::int32_t>(global_variable::nnodes));
  AppendScalar(buffer, static_cast<std::int32_t>(global_variable::ranks_per_node));
  AppendScalar(buffer, static_cast<std::uint32_t>(fields.size()));
  AppendScalar(buffer, static_cast<std::uint32_t>(0));
  AppendBytes(buffer, fields.data(), fields.size());
  return buffer;
}

std::vector<char> TrackedParticleOutput::CompactTrackFrame(
    Mesh *pm, const std::vector<float> &records) const {
  std::vector<char> buffer;
  const std::uint32_t record_count =
      static_cast<std::uint32_t>(records.size()/kTrackRecordFields);
  const std::uint32_t payload_bytes =
      static_cast<std::uint32_t>(records.size()*sizeof(float));
  buffer.reserve(kCompactTrackFrameBytes + payload_bytes);
  AppendBytes(buffer, kCompactTrackFrameMagic, sizeof(kCompactTrackFrameMagic));
  AppendScalar(buffer, kCompactTrackVersion);
  AppendScalar(buffer, kCompactTrackFrameBytes);
  AppendScalar(buffer, record_count);
  AppendScalar(buffer, static_cast<std::int64_t>(pm->ncycle));
  AppendScalar(buffer, static_cast<double>(pm->time));
  AppendScalar(buffer, payload_bytes);
  AppendScalar(buffer, static_cast<std::uint32_t>(0));
  AppendBytes(buffer, records.data(), records.size()*sizeof(float));
  return buffer;
}

std::vector<float> TrackedParticleOutput::PackLocalTrackRecords(Mesh *pm) const {
  std::vector<float> records(static_cast<std::size_t>(npout)*kTrackRecordFields);
  for (int p=0; p<npout; ++p) {
    const int base = p*kTrackRecordFields;
    records[base     ] = static_cast<float>(outpart(p).tag);
    records[base +  1] = static_cast<float>(pm->time);
    records[base +  2] = static_cast<float>(outpart(p).x);
    records[base +  3] = static_cast<float>(outpart(p).y);
    records[base +  4] = static_cast<float>(outpart(p).z);
    records[base +  5] = static_cast<float>(outpart(p).vx);
    records[base +  6] = static_cast<float>(outpart(p).vy);
    records[base +  7] = static_cast<float>(outpart(p).vz);
    records[base +  8] = static_cast<float>(outpart(p).Bx);
    records[base +  9] = static_cast<float>(outpart(p).By);
    records[base + 10] = static_cast<float>(outpart(p).Bz);
    records[base + 11] = static_cast<float>(outpart(p).K1);
    records[base + 12] = static_cast<float>(outpart(p).K2);
    records[base + 13] = static_cast<float>(outpart(p).K3);
    records[base + 14] = static_cast<float>(outpart(p).dB1);
    records[base + 15] = static_cast<float>(outpart(p).dB2);
    records[base + 16] = static_cast<float>(outpart(p).dB3);
    records[base + 17] = static_cast<float>(outpart(p).jmag);
  }
  return records;
}

void TrackedParticleOutput::AppendTrackBuffer(Mesh *pm,
                                              const std::vector<float> &records) {
  if (!IsShardWriter(out_params.file_shard_mode)) {return;}
  if (track_header_format == TrackHeaderFormat::compact) {
    std::vector<char> frame = CompactTrackFrame(pm, records);
    track_buffer.insert(track_buffer.end(), frame.begin(), frame.end());
  } else {
    int record_count = static_cast<int>(records.size()/kTrackRecordFields);
    std::string header = TrackHeader(pm, record_count);
    track_buffer.insert(track_buffer.end(), header.begin(), header.end());
    const char *payload = reinterpret_cast<const char*>(records.data());
    track_buffer.insert(track_buffer.end(), payload,
                        payload + records.size()*sizeof(float));
  }
  track_cycles_buffered += 1;
}

bool TrackedParticleOutput::TrackFileHasBytes(const std::string &fname) const {
  struct stat st;
  if (stat(fname.c_str(), &st) != 0) {return false;}
  return st.st_size > 0;
}

bool TrackedParticleOutput::TrackFileStartsWithCompactMagic(const std::string &fname) const {
  char magic[sizeof(kCompactTrackFileMagic)] = {};
  FILE *pfile = std::fopen(fname.c_str(), "rb");
  if (pfile == nullptr) {return false;}
  std::size_t nread = std::fread(magic, 1, sizeof(magic), pfile);
  std::fclose(pfile);
  return nread == sizeof(magic) &&
         std::memcmp(magic, kCompactTrackFileMagic, sizeof(magic)) == 0;
}

void TrackedParticleOutput::ValidateTrackFileAppendFormat(
    const std::string &fname) const {
  if (!TrackFileHasBytes(fname)) {return;}
  bool existing_compact = TrackFileStartsWithCompactMagic(fname);
  bool requested_compact = (track_header_format == TrackHeaderFormat::compact);
  if (existing_compact == requested_compact) {return;}
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Tracked particle output file '" << fname
            << "' already exists with a different trk_header_format. "
            << "Use a new basename or run directory before changing "
            << "trk_header_format." << std::endl;
  std::exit(EXIT_FAILURE);
}

void TrackedParticleOutput::FlushTrackBuffer(const std::string &fname) {
  if (track_buffer.empty()) {return;}
  ValidateTrackFileAppendFormat(fname);
  bool write_compact_prologue =
      (track_header_format == TrackHeaderFormat::compact && !TrackFileHasBytes(fname));
  FILE *pfile;
  if ((pfile = std::fopen(fname.c_str(),"ab")) == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (write_compact_prologue) {
    std::vector<char> prologue = CompactTrackPrologue();
    if (!WriteBytes(pfile, prologue.data(), prologue.size(),
                    "Compact tracked particle prologue")) {
      std::fclose(pfile);
      std::exit(EXIT_FAILURE);
    }
  }
  if (!WriteBytes(pfile, track_buffer.data(), track_buffer.size(),
                  "Tracked particle buffer")) {
    std::fclose(pfile);
    std::exit(EXIT_FAILURE);
  }
  std::fclose(pfile);
  track_buffer.clear();
  track_cycles_buffered = 0;
}

void TrackedParticleOutput::ValidateTrackedRecords() {
  int local_errors = 0;
  std::vector<int> local_seen;
  if (track_validate_global_tags) {
    local_seen.assign(ntrack_total, 0);
  }
  std::vector<int> local_tags;
  local_tags.reserve(npout);
  for (int p=0; p<npout; ++p) {
    int tag = outpart(p).tag;
    if (tag < 0 || tag >= ntrack_total) {
      local_errors += 1;
      continue;
    }
    local_tags.push_back(tag);
    if (track_validate_global_tags) {
      local_seen[tag] += 1;
    }
  }
  std::sort(local_tags.begin(), local_tags.end());
  for (std::size_t n=1; n<local_tags.size(); ++n) {
    if (local_tags[n] == local_tags[n-1]) {
      local_errors += 1;
    }
  }

  int nout_total = npout;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&npout, &nout_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  int fatal_error = 0;
  if (nout_total != ntrack_total) {
    if (global_variable::my_rank == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Tracked particle output found " << nout_total
                << " global records, but expected " << ntrack_total << std::endl;
    }
    fatal_error = 1;
  }

  int global_local_errors = local_errors;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&local_errors, &global_local_errors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  if (global_local_errors != 0) {
    if (global_variable::my_rank == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Tracked particle output found "
                << global_local_errors
                << " local duplicate or out-of-range tags" << std::endl;
    }
    fatal_error = 1;
  }

  if (track_validate_global_tags) {
    std::vector<int> global_seen(ntrack_total, 0);
#if MPI_PARALLEL_ENABLED
    MPI_Reduce(local_seen.data(), global_seen.data(), ntrack_total, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
#else
    global_seen = local_seen;
#endif
    if (global_variable::my_rank == 0) {
      int first_bad_tag = -1;
      int first_bad_count = 0;
      for (int t=0; t<ntrack_total; ++t) {
        if (global_seen[t] != 1) {
          first_bad_tag = t;
          first_bad_count = global_seen[t];
          break;
        }
      }
      if (first_bad_tag >= 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Tracked particle tag " << first_bad_tag
                  << " appears " << first_bad_count
                  << " times globally; expected exactly once" << std::endl;
        fatal_error = 1;
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&fatal_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  }

  if (fatal_error != 0) {
    std::exit(EXIT_FAILURE);
  }
}

void TrackedParticleOutput::WriteSharedTrackFrame(Mesh *pm,
                                                  const std::vector<float> &records) {
  const std::string fname = TrackFilename();
  std::vector<int> counts = GatherShardCounts(npout, FileShardMode::shared);
  int prefix_records = PrefixCountBeforeMe(counts, FileShardMode::shared);
  int record_count = 0;
  for (int count : counts) {
    record_count += count;
  }

#if MPI_PARALLEL_ENABLED
  if (global_variable::my_rank == 0) {
    ValidateTrackFileAppendFormat(fname);
    bool write_compact_prologue =
        (track_header_format == TrackHeaderFormat::compact && !TrackFileHasBytes(fname));
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"ab")) == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname
                << "' could not be opened" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (write_compact_prologue) {
      std::vector<char> prologue = CompactTrackPrologue();
      if (!WriteBytes(pfile, prologue.data(), prologue.size(),
                      "Compact tracked particle prologue")) {
        std::fclose(pfile);
        std::exit(EXIT_FAILURE);
      }
    }
    if (track_header_format == TrackHeaderFormat::compact) {
      std::vector<float> empty_records;
      std::vector<char> frame = CompactTrackFrame(pm, empty_records);
      std::uint32_t global_records = static_cast<std::uint32_t>(record_count);
      std::uint32_t payload_bytes = static_cast<std::uint32_t>(
          record_count*kTrackRecordFields*sizeof(float));
      std::memcpy(frame.data() + 12, &global_records, sizeof(global_records));
      std::memcpy(frame.data() + 32, &payload_bytes, sizeof(payload_bytes));
      if (!WriteBytes(pfile, frame.data(), kCompactTrackFrameBytes,
                      "Compact tracked particle frame header")) {
        std::fclose(pfile);
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::string header = TrackHeader(pm, record_count);
      if (!WriteBytes(pfile, header.data(), header.size(), "Tracked particle header")) {
        std::fclose(pfile);
        std::exit(EXIT_FAILURE);
      }
    }
    std::fclose(pfile);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File fh;
  int errcode = MPI_File_open(MPI_COMM_WORLD, fname.c_str(), MPI_MODE_WRONLY,
                              MPI_INFO_NULL, &fh);
  if (errcode != MPI_SUCCESS) {
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    Kokkos::printf("%.*s\n", resultlen, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Offset payload_offset = 0;
  errcode = MPI_File_get_size(fh, &payload_offset);
  if (errcode != MPI_SUCCESS) {
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    Kokkos::printf("%.*s\n", resultlen, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Offset myoffset = payload_offset +
      static_cast<MPI_Offset>(prefix_records)*kTrackRecordFields*sizeof(float);
  MPI_Status status;
  errcode = MPI_File_write_at_all(fh, myoffset,
                                  const_cast<float*>(records.data()),
                                  static_cast<int>(records.size()), MPI_FLOAT, &status);
  if (errcode != MPI_SUCCESS) {
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    Kokkos::printf("%.*s\n", resultlen, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_File_close(&fh);
#else
  ValidateTrackFileAppendFormat(fname);
  bool write_compact_prologue =
      (track_header_format == TrackHeaderFormat::compact && !TrackFileHasBytes(fname));
  FILE *pfile;
  if ((pfile = std::fopen(fname.c_str(),"ab")) == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Output file '" << fname
              << "' could not be opened" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (write_compact_prologue) {
    std::vector<char> prologue = CompactTrackPrologue();
    if (!WriteBytes(pfile, prologue.data(), prologue.size(),
                    "Compact tracked particle prologue")) {
      std::fclose(pfile);
      std::exit(EXIT_FAILURE);
    }
  }
  if (track_header_format == TrackHeaderFormat::compact) {
    std::vector<char> frame = CompactTrackFrame(pm, records);
    if (!WriteBytes(pfile, frame.data(), frame.size(), "Tracked particle data")) {
      std::fclose(pfile);
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::string header = TrackHeader(pm, record_count);
    if (!WriteBytes(pfile, header.data(), header.size(), "Tracked particle header") ||
        !WriteBytes(pfile, records.data(), records.size(), "Tracked particle data")) {
      std::fclose(pfile);
      std::exit(EXIT_FAILURE);
    }
  }
  std::fclose(pfile);
#endif
}

void TrackedParticleOutput::WriteNodeTrackFrame(Mesh *pm,
                                                const std::vector<float> &records) {
#if MPI_PARALLEL_ENABLED
  std::vector<int> counts = GatherShardCounts(npout, FileShardMode::per_node);
  std::vector<int> recv_counts(counts.size(), 0);
  std::vector<int> displs(counts.size(), 0);
  int record_count = 0;
  int float_count = 0;
  for (std::size_t n=0; n<counts.size(); ++n) {
    recv_counts[n] = counts[n]*kTrackRecordFields;
    displs[n] = float_count;
    float_count += recv_counts[n];
    record_count += counts[n];
  }

  std::vector<float> node_records;
  if (IsShardWriter(FileShardMode::per_node)) {
    node_records.resize(float_count);
  }
  MPI_Gatherv(const_cast<float*>(records.data()), static_cast<int>(records.size()),
              MPI_FLOAT,
              node_records.empty() ? nullptr : node_records.data(),
              recv_counts.data(), displs.data(), MPI_FLOAT, 0,
              global_variable::node_comm);
  if (IsShardWriter(FileShardMode::per_node)) {
    AppendTrackBuffer(pm, node_records);
    if (record_count != static_cast<int>(node_records.size()/kTrackRecordFields)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Node tracked-particle gather produced an "
                << "inconsistent record count" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#else
  AppendTrackBuffer(pm, records);
#endif
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
