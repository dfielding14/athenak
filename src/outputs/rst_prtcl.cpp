//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rst_prtcl.cpp
//! \brief writes sidecar restart files for particles

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

void FatalParticleRestartWrite(const char *file, int line, const std::string &msg) {
  std::cout << "### FATAL ERROR in " << file << " at line " << line << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn ParticleRestartOutput::ParticleRestartOutput()
//! \brief Construct sidecar particle restart output.

ParticleRestartOutput::ParticleRestartOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
    BaseTypeOutput(pin, pm, op) {
  if (pm->pmb_pack->ppart == nullptr) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "file_type=rst_prtcl requires a <particles> block.");
  }
  if (pm->pmb_pack->ppart->particle_type != ParticleType::star) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "file_type=rst_prtcl currently supports particle_type=star only.");
  }
  if (mkdir("rst_prtcl", 0775) != 0 && errno != EEXIST) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "Unable to create rst_prtcl output directory.");
  }
}

//----------------------------------------------------------------------------------------
//! \fn ParticleRestartOutput::LoadOutputData()
//! \brief Copy particle arrays to host before writing a sidecar restart.

void ParticleRestartOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  nprtcl_thisrank = pp->nprtcl_thispack;
  nprtcl_total = pm->nprtcl_total;
  nrdata = pp->nrdata;
  nidata = pp->nidata;
  local_max_tag = -1;
  local_star_mass_formed_total = pp->star_mass_formed_total;
  local_star_mass_accreted_total = pp->star_mass_accreted_total;

  Kokkos::realloc(outpart_rdata, nrdata, nprtcl_thisrank);
  Kokkos::realloc(outpart_idata, nidata, nprtcl_thisrank);
  if (nprtcl_thisrank > 0) {
    Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
    Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
    for (int p=0; p<nprtcl_thisrank; ++p) {
      local_max_tag = std::max(local_max_tag, outpart_idata(PTAG,p));
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn ParticleRestartOutput::WriteOutputFile()
//! \brief Write a single global particle sidecar restart file.

void ParticleRestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::string fname;
  char number[7];
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  fname = std::string("rst_prtcl/") + out_params.file_basename + number +
      ".rst_prtcl";

  std::vector<int> nprtcl_eachrank(global_variable::nranks, 0);
  nprtcl_eachrank[global_variable::my_rank] = nprtcl_thisrank;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&nprtcl_thisrank, 1, MPI_INT, nprtcl_eachrank.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
#endif

  int total_particles = 0;
  int my_particle_offset = 0;
  for (int n=0; n<global_variable::nranks; ++n) {
    if (n < global_variable::my_rank) {my_particle_offset += nprtcl_eachrank[n];}
    total_particles += nprtcl_eachrank[n];
  }
  if (total_particles != nprtcl_total) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "Particle count mismatch while writing rst_prtcl sidecar.");
  }

  int nrdata_min = nrdata;
  int nrdata_max = nrdata;
  int nidata_min = nidata;
  int nidata_max = nidata;
  int max_tag = local_max_tag;
  Real formed_total = local_star_mass_formed_total;
  Real accreted_total = local_star_mass_accreted_total;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &nrdata_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nrdata_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nidata_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nidata_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max_tag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &formed_total, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &accreted_total, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  if (nrdata_min != nrdata_max || nidata_min != nidata_max) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "Particle restart cannot write ranks with different particle array sizes.");
  }

  particles::ParticleRestartHeader header;
  std::memset(&header, 0, sizeof(header));
  std::strncpy(header.magic, particles::kParticleRestartMagic,
               particles::kParticleRestartMagicSize - 1);
  header.version = particles::kParticleRestartVersion;
  header.header_size = sizeof(header);
  header.endian_marker = particles::kParticleRestartEndianMarker;
  header.real_size = sizeof(Real);
  header.int_size = sizeof(int);
  header.ncycle = pm->ncycle;
  header.total_particles = total_particles;
  header.nrdata = nrdata;
  header.nidata = nidata;
  header.max_tag = max_tag;
  header.time = pm->time;
  header.star_mass_formed_total = formed_total;
  header.star_mass_accreted_total = accreted_total;
  if (out_params.file_basename.size() >=
      static_cast<std::size_t>(particles::kParticleRestartBasenameSize)) {
    FatalParticleRestartWrite(__FILE__, __LINE__,
        "Particle restart basename is too long for the sidecar header.");
  }
  header.basename_length = static_cast<int>(out_params.file_basename.size());
  std::strncpy(header.basename, out_params.file_basename.c_str(),
               particles::kParticleRestartBasenameSize - 1);

  IOWrapper prtcl_file;
  prtcl_file.Open(fname.c_str(), IOWrapper::FileMode::write);
  if (global_variable::my_rank == 0) {
    if (prtcl_file.Write_any_type(&header, sizeof(header), "byte") != sizeof(header)) {
      FatalParticleRestartWrite(__FILE__, __LINE__,
          "Particle restart header was not written correctly.");
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (nprtcl_thisrank > 0) {
    std::vector<Real> rdata_transposed(nprtcl_thisrank*nrdata);
    std::vector<int> idata_transposed(nprtcl_thisrank*nidata);
    for (int p=0; p<nprtcl_thisrank; ++p) {
      for (int n=0; n<nrdata; ++n) {
        rdata_transposed[p*nrdata + n] = outpart_rdata(n,p);
      }
      for (int n=0; n<nidata; ++n) {
        idata_transposed[p*nidata + n] = outpart_idata(n,p);
      }
    }

    IOWrapperSizeT rdata_offset = sizeof(header) +
        static_cast<IOWrapperSizeT>(my_particle_offset)*nrdata*sizeof(Real);
    IOWrapperSizeT idata_offset = sizeof(header) +
        static_cast<IOWrapperSizeT>(total_particles)*nrdata*sizeof(Real) +
        static_cast<IOWrapperSizeT>(my_particle_offset)*nidata*sizeof(int);

    if (prtcl_file.Write_any_type_at(rdata_transposed.data(),
                                     nprtcl_thisrank*nrdata, rdata_offset,
                                     "Real") !=
        static_cast<std::size_t>(nprtcl_thisrank*nrdata)) {
      FatalParticleRestartWrite(__FILE__, __LINE__,
          "Particle real data were not written correctly to rst_prtcl sidecar.");
    }
    if (prtcl_file.Write_any_type_at(idata_transposed.data(),
                                     nprtcl_thisrank*nidata, idata_offset,
                                     "int") !=
        static_cast<std::size_t>(nprtcl_thisrank*nidata)) {
      FatalParticleRestartWrite(__FILE__, __LINE__,
          "Particle integer data were not written correctly to rst_prtcl sidecar.");
    }
  }
  prtcl_file.Close();

  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
