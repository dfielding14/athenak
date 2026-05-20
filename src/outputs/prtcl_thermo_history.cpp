//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file prtcl_thermo_history.cpp
//! \brief append-only binary thermodynamic history output for lagrangian_mc particles

#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

struct FileHeader {
  char magic[32];
  int version;
  int real_size;
  int nscalars;
  int reserved;
};

struct BlockHeader {
  char magic[16];
  int version;
  int nrecords;
  int nscalars;
  int int_per_record;
  int real_per_record;
  int cycle;
  Real time;
};

void WriteFileHeaderIfNeeded(FILE *file, int nscalars) {
  std::fseek(file, 0, SEEK_END);
  if (std::ftell(file) != 0) return;

  FileHeader header;
  std::memset(&header, 0, sizeof(header));
  std::strncpy(header.magic, "ATHK_PRTCL_THERMO_HISTORY", sizeof(header.magic)-1);
  header.version = 1;
  header.real_size = sizeof(Real);
  header.nscalars = nscalars;
  if (std::fwrite(&header, sizeof(header), 1, file) != 1) {
    std::cout << "### FATAL ERROR: failed to write particle thermo file header"
              << std::endl;
    std::exit(EXIT_FAILURE);
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
    nscalars(0) {
  mkdir("prtcl_thermo_history", 0775);
  if (pm->pmb_pack->ppart == nullptr ||
      !pm->pmb_pack->ppart->IsLagrangianMC()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "file_type=prtcl_thermo_history requires "
              << "particle_type=lagrangian_mc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  nscalars = pm->pmb_pack->ppart->GetLagrangianMCScalarCount();
}

//----------------------------------------------------------------------------------------
// LoadOutputData

void ParticleThermoHistoryOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  pp->SampleThermodynamics();
  npout_thisrank = pm->nprtcl_thisrank;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  if (npout_thisrank > 0) {
    Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
    Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
  }
}

//----------------------------------------------------------------------------------------
// WriteOutputFile

void ParticleThermoHistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  constexpr int int_per_record = 4;   // cycle, tag, seed_id, gid
  int real_per_record = 12 + nscalars; // time, x, thermo, velocity, scalars

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
    local_r[real_per_record*p + 4] = outpart_rdata(LMC_RHO,p);
    local_r[real_per_record*p + 5] = outpart_rdata(LMC_PRESSURE,p);
    local_r[real_per_record*p + 6] = outpart_rdata(LMC_TEMPERATURE,p);
    local_r[real_per_record*p + 7] = outpart_rdata(LMC_ENTROPY,p);
    local_r[real_per_record*p + 8] = outpart_rdata(LMC_EINT,p);
    local_r[real_per_record*p + 9] = outpart_rdata(LMC_VX,p);
    local_r[real_per_record*p + 10] = outpart_rdata(LMC_VY,p);
    local_r[real_per_record*p + 11] = outpart_rdata(LMC_VZ,p);
    for (int n=0; n<nscalars; ++n) {
      local_r[real_per_record*p + 12 + n] = outpart_rdata(LMC_SCALAR0+n,p);
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
  MPI_Gatherv(local_i.data(), local_i.size(), MPI_INT,
              all_i.data(), counts.data(), idispl.data(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_r.data(), local_r.size(), MPI_ATHENA_REAL,
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
    WriteFileHeaderIfNeeded(file, nscalars);

    BlockHeader header;
    std::memset(&header, 0, sizeof(header));
    std::strncpy(header.magic, "ATHKTHPBLK", sizeof(header.magic)-1);
    header.version = 1;
    header.nrecords = npout_total;
    header.nscalars = nscalars;
    header.int_per_record = int_per_record;
    header.real_per_record = real_per_record;
    header.cycle = pm->ncycle;
    header.time = pm->time;
    std::fwrite(&header, sizeof(header), 1, file);
    if (!all_i.empty()) std::fwrite(all_i.data(), sizeof(int), all_i.size(), file);
    if (!all_r.empty()) std::fwrite(all_r.data(), sizeof(Real), all_r.size(), file);
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
