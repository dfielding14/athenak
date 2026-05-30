//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rst_prtcl.cpp
//! \brief writes per-rank particle restart dumps

#include <sys/stat.h>

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particle_restart.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

uint64_t UpdateFNV1a(uint64_t hash, const void *data, std::size_t nbyte) {
  const unsigned char *bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i=0; i<nbyte; ++i) {
    hash ^= static_cast<uint64_t>(bytes[i]);
    hash *= 1099511628211ULL;
  }
  return hash;
}

void WriteBytesAndRename(const std::string &fname,
                         const std::vector<unsigned char> &bytes) {
  std::string tmp_fname = fname + ".tmp";
  FILE *file = std::fopen(tmp_fname.c_str(), "wb");
  bool published = file != nullptr;
  if (published) {
    published = std::fwrite(bytes.data(), 1, bytes.size(), file) == bytes.size();
    published = std::fclose(file) == 0 && published;
    if (published) {published = std::rename(tmp_fname.c_str(), fname.c_str()) == 0;}
  }
  if (!published) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Relativistic typed-v2 restart shard '" << fname
              << "' could not be published" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void PublishRelativisticManifest(
    const std::string &fname, const std::string &mesh_restart,
    const particles::restart::V2Header &local_header,
    const particles::restart::MeshWitness &mesh_witness) {
  std::vector<std::uint64_t> local_counts(global_variable::nranks, 0);
  std::vector<std::uint64_t> byte_counts(global_variable::nranks, 0);
  std::vector<std::uint64_t> payload_checksums(global_variable::nranks, 0);
  std::vector<std::uint64_t> header_checksums(global_variable::nranks, 0);
#if MPI_PARALLEL_ENABLED
  MPI_Gather(&local_header.local_count, 1, MPI_UINT64_T, local_counts.data(), 1,
             MPI_UINT64_T, 0, MPI_COMM_WORLD);
  std::uint64_t byte_count =
      particles::restart::kV2HeaderBytes + local_header.payload_bytes;
  MPI_Gather(&byte_count, 1, MPI_UINT64_T, byte_counts.data(), 1, MPI_UINT64_T,
             0, MPI_COMM_WORLD);
  MPI_Gather(&local_header.payload_checksum, 1, MPI_UINT64_T,
             payload_checksums.data(), 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Gather(&local_header.header_checksum, 1, MPI_UINT64_T,
             header_checksums.data(), 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#else
  local_counts[0] = local_header.local_count;
  byte_counts[0] = particles::restart::kV2HeaderBytes + local_header.payload_bytes;
  payload_checksums[0] = local_header.payload_checksum;
  header_checksums[0] = local_header.header_checksum;
#endif
  if (global_variable::my_rank == 0) {
    std::string manifest_fname =
        particles::restart::ManifestPathFromRankZeroShard(fname);
    std::string tmp_fname = manifest_fname + ".tmp";
    std::ofstream manifest(tmp_fname);
    manifest << std::setprecision(17);
    manifest << "magic AKPRST-MANIFEST\n";
    manifest << "version 2.0\n";
    manifest << "topology_policy reject_rank_count_change\n";
    manifest << "paired_mesh_checkpoint required\n";
    manifest << "checkpoint_id " << local_header.checkpoint_id << '\n';
    manifest << "saved_nranks " << local_header.saved_nranks << '\n';
    manifest << "global_count " << local_header.global_count << '\n';
    manifest << "mesh_cycle " << local_header.mesh_cycle << '\n';
    manifest << "mesh_time " << local_header.mesh_time << '\n';
    manifest << "mesh_dt " << local_header.mesh_dt << '\n';
    manifest << "particle_dtnew " << local_header.particle_dtnew << '\n';
    manifest << "config_fingerprint " << local_header.config_fingerprint << '\n';
    manifest << "mesh_byte_count " << mesh_witness.mesh_byte_count << '\n';
    manifest << "mesh_checksum " << mesh_witness.mesh_checksum << '\n';
    manifest << "mesh_topology_hash " << mesh_witness.mesh_topology_hash << '\n';
    manifest << "mesh_restart " << mesh_restart << '\n';
    for (int rank=0; rank<global_variable::nranks; ++rank) {
      char rank_dir[32];
      std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", rank);
      std::size_t slash = fname.rfind('/');
      std::string basename = fname.substr(slash + 1);
      manifest << "shard " << rank << ' ' << rank_dir << basename << ' '
               << byte_counts[rank] << ' ' << local_counts[rank] << ' '
               << payload_checksums[rank] << ' ' << header_checksums[rank] << '\n';
    }
    manifest.close();
    if (!manifest || std::rename(tmp_fname.c_str(), manifest_fname.c_str()) != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Relativistic typed-v2 restart manifest '"
                << manifest_fname << "' could not be published" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

} // namespace

ParticleRestartOutput::ParticleRestartOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  mkdir("prst",0775);
  char rank_dir[32];
  std::snprintf(rank_dir, sizeof(rank_dir), "prst/rank_%08d/", global_variable::my_rank);
  mkdir(rank_dir,0775);

  int prtcl_rst_flag = pin->GetOrAddInteger("problem","prtcl_rst_flag",0);
  if (prtcl_rst_flag) {
    std::string prst_fname = pin->GetString("problem","prtcl_res_file");
    particles::Particles *pp = pm->pmb_pack->ppart;
    if (pp != nullptr && pp->pusher == ParticlesPusher::relativistic_hc) {
      try {
        out_params.file_number = particles::restart::FileNumber(prst_fname) + 1;
        out_params.last_time = pm->time;
      } catch (const std::exception &error) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Relativistic typed-v2 restart output initialization "
                  << "rejected: " << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::size_t last_period = prst_fname.rfind('.');
      if (last_period != std::string::npos && last_period >= 5) {
        std::string outnumber_str = prst_fname.substr(last_period-5,5);
        out_params.file_number = std::stoi(outnumber_str) + 1;
        out_params.last_time = pm->time;
      }
    }
  }
}

void ParticleRestartOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  pp->CheckConsistency("particle restart output");
  if (pp->pusher == ParticlesPusher::relativistic_hc) {
    pp->ValidateRelativisticRestartState("typed-v2 particle restart output");
  }
  npout_thisrank = pp->nprtcl_thispack;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
}

void ParticleRestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  if (pp->pusher == ParticlesPusher::relativistic_hc &&
      (out_params.file_number < 0 || out_params.file_number > 99999)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Relativistic typed-v2 particle restart file_number must "
              << "remain in the five-digit range [0, 99999]" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  char rank_dir[32];
  char number[8];
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  std::string fname = "prst/" + std::string(rank_dir) + out_params.file_basename +
                      std::string(number) + ".prst";

  if (pp->pusher == ParticlesPusher::relativistic_hc) {
    try {
      std::string number_string(number);
      std::string mesh_restart =
          "rst/" + out_params.file_basename + number_string + ".rst";
      std::ifstream mesh_file(mesh_restart, std::ios::binary);
      particles::restart::Require(
          mesh_file.good(),
          "typed-v2 particle restart requires paired mesh checkpoint '" +
          mesh_restart + "' to be published first");
      std::vector<std::int32_t> idata(
          particles::restart::kV2IntegerFields*npout_thisrank);
      std::vector<double> rdata(particles::restart::kV2RealFields*npout_thisrank);
      for (int p=0; p<npout_thisrank; ++p) {
        idata[particles::restart::kV2IntegerFields*p] = outpart_idata(PGID,p);
        idata[particles::restart::kV2IntegerFields*p + 1] = outpart_idata(PTAG,p);
        idata[particles::restart::kV2IntegerFields*p + 2] = outpart_idata(PSP,p);
        for (std::uint32_t n=0; n<particles::restart::kV2RealFields; ++n) {
          rdata[particles::restart::kV2RealFields*p + n] = outpart_rdata(n,p);
        }
      }
      particles::restart::V2Header header;
      header.saved_nranks = global_variable::nranks;
      header.saved_rank = global_variable::my_rank;
      header.local_count = npout_thisrank;
      header.global_count = npout_total;
      header.mesh_cycle = pm->ncycle;
      header.mesh_time = pm->time;
      header.mesh_dt = pm->dt;
      header.particle_dtnew = pp->dtnew;
      std::uint32_t field_source =
          (pp->relativistic_field_source == RelativisticFieldSource::prescribed_test) ?
          1U : 2U;
      std::uint32_t temporal_sampling =
          (pp->relativistic_temporal_sampling == RelativisticTemporalSampling::none) ?
          0U : 1U;
      header.config_fingerprint =
          particles::restart::ComputeRelativisticConfigFingerprint(
              pp->c_model, pp->alpha_s, field_source, temporal_sampling,
              pp->subcycle, pp->subcycle_strict, pp->subcycle_max_steps,
              pp->subcycle_cell_fraction, pp->subcycle_meshblock_fraction,
              pp->subcycle_gyro_fraction, pp->subcycle_electric_kick_max);
      header.checkpoint_id = particles::restart::ComputeCheckpointID(
          pm->ncycle, pm->time, pm->dt, out_params.file_number,
          global_variable::nranks);
      particles::restart::MeshWitness mesh_witness =
          particles::restart::ReadMeshWitness(mesh_restart);
      particles::restart::Require(
          mesh_witness.checkpoint_id == header.checkpoint_id &&
          mesh_witness.saved_nranks == header.saved_nranks &&
          mesh_witness.mesh_cycle == header.mesh_cycle &&
          mesh_witness.mesh_time == header.mesh_time &&
          mesh_witness.mesh_dt == header.mesh_dt,
          "typed-v2 particle checkpoint does not match paired mesh witness");
      std::vector<unsigned char> bytes =
          particles::restart::EncodeShard(header, idata, rdata);
      WriteBytesAndRename(fname, bytes);
      PublishRelativisticManifest(fname, mesh_restart, header, mesh_witness);
    } catch (const std::exception &error) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Relativistic typed-v2 restart output rejected: "
                << error.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    out_params.file_number++;
    out_params.last_time = (out_params.last_time < 0.0) ? pm->time :
                           out_params.last_time + out_params.dt;
    pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
    pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
    return;
  }

  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::write, true);

  Real header[3] = {pm->time, pm->dt, static_cast<Real>(npout_thisrank)};
  if (partfile.Write_any_type(header,3,"Real",true) != 3) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle restart header output failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real *real_data = new Real[std::max(1,17*npout_thisrank)];
  for (int p=0; p<npout_thisrank; ++p) {
    real_data[17*p] = static_cast<Real>(outpart_idata(PGID,p));
    real_data[17*p + 1] = static_cast<Real>(outpart_idata(PTAG,p));
    real_data[17*p + 2] = static_cast<Real>(outpart_idata(PSP,p));
    real_data[17*p + 3] = outpart_rdata(IPX,p);
    real_data[17*p + 4] = outpart_rdata(IPY,p);
    real_data[17*p + 5] = outpart_rdata(IPZ,p);
    real_data[17*p + 6] = outpart_rdata(IPVX,p);
    real_data[17*p + 7] = outpart_rdata(IPVY,p);
    real_data[17*p + 8] = outpart_rdata(IPVZ,p);
    real_data[17*p + 9] = outpart_rdata(IPM,p);
    real_data[17*p + 10] = outpart_rdata(IPBX,p);
    real_data[17*p + 11] = outpart_rdata(IPBY,p);
    real_data[17*p + 12] = outpart_rdata(IPBZ,p);
    real_data[17*p + 13] = outpart_rdata(IPDX,p);
    real_data[17*p + 14] = outpart_rdata(IPDY,p);
    real_data[17*p + 15] = outpart_rdata(IPDZ,p);
    real_data[17*p + 16] = outpart_rdata(IPDB,p);
  }

  if (partfile.Write_any_type(real_data,17*npout_thisrank,"Real",true) !=
      static_cast<std::size_t>(17*npout_thisrank)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle restart data output failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  partfile.Close(true);

  uint64_t checksum = 1469598103934665603ULL;
  checksum = UpdateFNV1a(checksum, header, 3*sizeof(Real));
  checksum = UpdateFNV1a(checksum, real_data,
                         static_cast<std::size_t>(17*npout_thisrank)*sizeof(Real));
  std::string meta_fname = fname + ".pmeta";
  FILE *mfile = std::fopen(meta_fname.c_str(),"w");
  if (mfile == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle restart metadata file '" << meta_fname
              << "' could not be opened" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::fprintf(mfile, "magic AKPRST\n");
  std::fprintf(mfile, "version 1\n");
  std::fprintf(mfile, "binary_payload legacy_prst\n");
  std::fprintf(mfile, "real_size %zu\n", sizeof(Real));
  std::fprintf(mfile, "integer_size %zu\n", sizeof(int));
  std::fprintf(mfile, "species_count %d\n", pp->nspecies);
  std::fprintf(mfile, "real_field_count %d\n", pp->nrdata);
  std::fprintf(mfile, "integer_field_count %d\n", pp->nidata);
  std::fprintf(mfile, "record_real_count 17\n");
  std::fprintf(mfile, "particle_count %d\n", npout_thisrank);
  std::fprintf(mfile, "byte_count %zu\n",
               static_cast<std::size_t>(3 + 17*npout_thisrank)*sizeof(Real));
  std::fprintf(mfile, "checksum_fnv1a64 %" PRIu64 "\n", checksum);
  std::fclose(mfile);

  delete[] real_data;

  out_params.file_number++;
  out_params.last_time = (out_params.last_time < 0.0) ? pm->time :
                         out_params.last_time + out_params.dt;
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
