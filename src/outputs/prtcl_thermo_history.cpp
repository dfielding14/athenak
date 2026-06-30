//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file prtcl_thermo_history.cpp
//! \brief append-only binary thermodynamic history output for lagrangian_mc particles

#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "particles/tracer_fields.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

struct FileHeader {
  char magic[32];
  int version;
  int real_size;
  int nfields;
  int names_size;
};

struct BlockHeader {
  char magic[16];
  int version;
  int nrecords;
  int int_per_record;
  int real_per_record;
  int cycle;
  Real time;
};

std::string JoinFieldNames(const std::vector<std::string> &names) {
  std::string joined;
  for (std::size_t n=0; n<names.size(); ++n) {
    if (n > 0) joined += "\n";
    joined += names[n];
  }
  return joined;
}

void FatalHistoryOutput(const std::string &msg) {
  std::cout << "### FATAL ERROR in prtcl_thermo_history.cpp" << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void EnsureFileHeader(FILE *file, const std::vector<std::string> &field_names) {
  const std::string names_blob = JoinFieldNames(field_names);
  std::fseek(file, 0, SEEK_END);
  std::int64_t file_size = std::ftell(file);
  if (file_size == 0) {
    FileHeader header;
    std::memset(&header, 0, sizeof(header));
    std::strncpy(header.magic, "ATHK_PRTCL_THERMO_HISTORY", sizeof(header.magic)-1);
    header.version = 2;
    header.real_size = sizeof(Real);
    header.nfields = static_cast<int>(field_names.size());
    header.names_size = static_cast<int>(names_blob.size());
    if (std::fwrite(&header, sizeof(header), 1, file) != 1 ||
        (!names_blob.empty() &&
         std::fwrite(names_blob.data(), sizeof(char), names_blob.size(), file) !=
         names_blob.size())) {
      FatalHistoryOutput("failed to write particle thermo file header");
    }
    return;
  }

  std::fseek(file, 0, SEEK_SET);
  FileHeader header;
  if (std::fread(&header, sizeof(header), 1, file) != 1) {
    FatalHistoryOutput("existing particle thermo file has a truncated header");
  }
  if (std::strncmp(header.magic, "ATHK_PRTCL_THERMO_HISTORY", 25) != 0) {
    FatalHistoryOutput("existing particle thermo file has an unrecognized header");
  }
  if (header.version != 2) {
    FatalHistoryOutput("existing particle thermo file uses an older schema; "
                       "write to a clean output directory or remove the old file");
  }
  if (header.real_size != static_cast<int>(sizeof(Real)) ||
      header.nfields != static_cast<int>(field_names.size()) ||
      header.names_size != static_cast<int>(names_blob.size())) {
    FatalHistoryOutput("existing particle thermo file schema does not match this output");
  }
  std::string existing_names(header.names_size, '\0');
  if (header.names_size > 0 &&
      std::fread(existing_names.data(), sizeof(char), header.names_size, file) !=
      static_cast<std::size_t>(header.names_size)) {
    FatalHistoryOutput("existing particle thermo file has a truncated schema");
  }
  if (existing_names != names_blob) {
    FatalHistoryOutput("existing particle thermo file variable list does not match this "
                       "output");
  }
  std::fseek(file, 0, SEEK_END);
}

void WriteBlockHeader(FILE *file, const BlockHeader &header) {
  int ints[5] = {header.version, header.nrecords, header.int_per_record,
                 header.real_per_record, header.cycle};
  if (std::fwrite(header.magic, sizeof(char), sizeof(header.magic), file) !=
      sizeof(header.magic) ||
      std::fwrite(ints, sizeof(int), 5, file) != 5 ||
      std::fwrite(&header.time, sizeof(Real), 1, file) != 1) {
    FatalHistoryOutput("failed to write particle thermo history block header");
  }
}

} // namespace

//----------------------------------------------------------------------------------------
// ctor

ParticleThermoHistoryOutput::ParticleThermoHistoryOutput(ParameterInput *pin, Mesh *pm,
                                                         OutputParameters op) :
    BaseTypeOutput(pin, pm, op),
    npout_thisrank(0),
    npout_total(0),
    tracer_gamma(5.0/3.0) {
  mkdir("prtcl_thermo_history", 0775);
  if (pm->pmb_pack->ppart == nullptr ||
      !pm->pmb_pack->ppart->IsLagrangianMC()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "file_type=prtcl_thermo_history requires "
              << "particle_type=lagrangian_mc" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool has_mhd = (pm->pmb_pack->pmhd != nullptr);
  tracer_gamma = has_mhd ? pin->GetOrAddReal("mhd", "gamma", 5.0/3.0) :
                           pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);
  int nscalars = pm->pmb_pack->ppart->GetLagrangianMCScalarCount();
  std::string field_list = pin->DoesParameterExist(op.block_name, "variables") ?
                           pin->GetString(op.block_name, "variables") :
                           pin->GetOrAddString("particles", "track_variables", "default");
  tracer_fields = particles::ParseTracerFieldList(field_list, has_mhd, nscalars,
                                                  op.block_name + "/variables");
  tracer_field_names = particles::TracerFieldNames(tracer_fields);
}

//----------------------------------------------------------------------------------------
// LoadOutputData

void ParticleThermoHistoryOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pm->nprtcl_thisrank;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  Kokkos::realloc(outfield_data, tracer_fields.size(), npout_thisrank);
  if (npout_thisrank <= 0) return;

  Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);

  auto &indcs = pm->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int ncells1 = indcs.nx1 + 2*indcs.ng;
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb_alloc = std::max(pm->pmb_pack->nmb_thispack, pm->nmb_maxperrank);
  int nfluid = 0, nscalars = 0;
  bool has_mhd = (pm->pmb_pack->pmhd != nullptr);

  if (pm->pmb_pack->phydro != nullptr) {
    nfluid = pm->pmb_pack->phydro->nhydro;
    nscalars = pm->pmb_pack->phydro->nscalars;
  } else {
    nfluid = pm->pmb_pack->pmhd->nmhd;
    nscalars = pm->pmb_pack->pmhd->nscalars;
  }

  HostArray5D<Real> h_w0("prtcl_history_w0", nmb_alloc, nfluid+nscalars,
                         ncells3, ncells2, ncells1);
  HostArray5D<Real> h_bcc("prtcl_history_bcc", nmb_alloc, NMAG,
                          ncells3, ncells2, ncells1);
  if (pm->pmb_pack->phydro != nullptr) {
    Kokkos::deep_copy(h_w0, pm->pmb_pack->phydro->w0);
  } else {
    Kokkos::deep_copy(h_w0, pm->pmb_pack->pmhd->w0);
    Kokkos::deep_copy(h_bcc, pm->pmb_pack->pmhd->bcc0);
  }

  pm->pmb_pack->pmb->mb_size.template sync<HostMemSpace>();
  auto h_size = pm->pmb_pack->pmb->mb_size.h_view;
  int nmb = pm->pmb_pack->nmb_thispack;
  int gids = pm->pmb_pack->gids;

  for (int p=0; p<npout_thisrank; ++p) {
    int m = outpart_idata(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      FatalHistoryOutput("tracer particle references a MeshBlock not owned by this rank");
    }
    int i = static_cast<int>((outpart_rdata(LMCX,p) - h_size(m).x1min)/
                             h_size(m).dx1) + is;
    int j = pm->multi_d ?
      static_cast<int>((outpart_rdata(LMCY,p) - h_size(m).x2min)/h_size(m).dx2) + js : js;
    int k = pm->three_d ?
      static_cast<int>((outpart_rdata(LMCZ,p) - h_size(m).x3min)/h_size(m).dx3) + ks : ks;
    i = std::max(is, std::min(is + indcs.nx1 - 1, i));
    j = std::max(js, std::min(js + indcs.nx2 - 1, j));
    k = std::max(ks, std::min(ks + indcs.nx3 - 1, k));
    for (int n=0; n<static_cast<int>(tracer_fields.size()); ++n) {
      outfield_data(n,p) = particles::EvaluateTracerFieldHost(tracer_fields[n], h_w0,
                                                              h_bcc, has_mhd,
                                                              tracer_gamma,
                                                              nfluid, m, k, j, i);
    }
  }
}

//----------------------------------------------------------------------------------------
// WriteOutputFile

void ParticleThermoHistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  constexpr int int_per_record = 4;   // cycle, tag, seed_id, gid
  int real_per_record = 4 + static_cast<int>(tracer_fields.size()); // time, x, fields

  std::vector<int> local_i(int_per_record*npout_thisrank);
  std::vector<Real> local_r(real_per_record*npout_thisrank);
  for (int p=0; p<npout_thisrank; ++p) {
    local_i[int_per_record*p] = pm->ncycle;
    local_i[int_per_record*p + 1] = outpart_idata(PTAG,p);
    local_i[int_per_record*p + 2] = outpart_idata(PSEEDID,p);
    local_i[int_per_record*p + 3] = outpart_idata(PGID,p);

    local_r[real_per_record*p] = pm->time;
    local_r[real_per_record*p + 1] = outpart_rdata(LMCX,p);
    local_r[real_per_record*p + 2] = outpart_rdata(LMCY,p);
    local_r[real_per_record*p + 3] = outpart_rdata(LMCZ,p);
    for (int n=0; n<static_cast<int>(tracer_fields.size()); ++n) {
      local_r[real_per_record*p + 4 + n] = outfield_data(n,p);
    }
  }

  std::vector<int> all_i;
  std::vector<Real> all_r;
#if MPI_PARALLEL_ENABLED
  std::vector<int> counts(global_variable::nranks), idispl(global_variable::nranks),
                   rcounts(global_variable::nranks), rdispl(global_variable::nranks);
  int local_count = npout_thisrank;
  MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) {
    int total_i = 0, total_r = 0;
    for (int n=0; n<global_variable::nranks; ++n) {
      idispl[n] = total_i;
      rdispl[n] = total_r;
      total_i += counts[n]*int_per_record;
      total_r += counts[n]*real_per_record;
      rcounts[n] = counts[n]*real_per_record;
      counts[n] *= int_per_record;
    }
    all_i.resize(total_i);
    all_r.resize(total_r);
  }
  int local_i_count = static_cast<int>(local_i.size());
  int local_r_count = static_cast<int>(local_r.size());
  MPI_Gatherv(local_i.data(), local_i_count, MPI_INT,
              all_i.data(), counts.data(), idispl.data(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_r.data(), local_r_count, MPI_ATHENA_REAL,
              all_r.data(), rcounts.data(), rdispl.data(), MPI_ATHENA_REAL, 0,
              MPI_COMM_WORLD);
#else
  all_i = std::move(local_i);
  all_r = std::move(local_r);
#endif

  if (global_variable::my_rank == 0) {
    std::string fname = "prtcl_thermo_history/" + out_params.file_basename + "." +
                        out_params.file_id + ".thp";
    FILE *file = std::fopen(fname.c_str(), "ab+");
    if (file == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    EnsureFileHeader(file, tracer_field_names);

    BlockHeader header;
    std::memset(&header, 0, sizeof(header));
    std::strncpy(header.magic, "ATHKTHPBLK", sizeof(header.magic)-1);
    header.version = 2;
    header.nrecords = npout_total;
    header.int_per_record = int_per_record;
    header.real_per_record = real_per_record;
    header.cycle = pm->ncycle;
    header.time = pm->time;
    WriteBlockHeader(file, header);
    bool write_failed =
        (!all_i.empty() && std::fwrite(all_i.data(), sizeof(int), all_i.size(), file) !=
         all_i.size()) ||
        (!all_r.empty() && std::fwrite(all_r.data(), sizeof(Real), all_r.size(), file) !=
         all_r.size());
    if (write_failed) {
      FatalHistoryOutput("failed to write particle thermo history block");
    }
    std::fclose(file);
  }

  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
