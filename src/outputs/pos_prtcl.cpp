//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pos_prtcl.cpp
//! \brief writes compact particle species and position dumps

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
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

namespace {

bool IsRelativisticParticleOutput(const particles::Particles *pp) {
  return pp->pusher == ParticlesPusher::relativistic_hc;
}

constexpr bool RelativisticPpdSpeciesExactlyRepresentable(const int species) {
  return species >= 0 &&
         static_cast<double>(static_cast<float>(species)) ==
             static_cast<double>(species);
}

static_assert(RelativisticPpdSpeciesExactlyRepresentable(16777216));
static_assert(!RelativisticPpdSpeciesExactlyRepresentable(16777217));
static_assert(RelativisticPpdSpeciesExactlyRepresentable(2147483520));

} // namespace

ParticlePositionsOutput::ParticlePositionsOutput(ParameterInput *pin, Mesh *pm,
                                                 OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  mkdir("ppd",0775);
}

void ParticlePositionsOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pp->nprtcl_thispack;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
}

void ParticlePositionsOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  char number[6];
  std::snprintf(number, sizeof(number), "%05d", out_params.file_number);
  std::string fname = "ppd/" + out_params.file_basename + "." +
                      std::string(number) + ".ppd";

  particles::Particles *pp = pm->pmb_pack->ppart;
  bool relativistic_mode = IsRelativisticParticleOutput(pp);
  float *data = new float[std::max(1,4*npout_thisrank)];
  int relativistic_failure = 0;
  for (int p=0; p<npout_thisrank; ++p) {
    int species = outpart_idata(PSP,p);
    data[4*p] = static_cast<float>(species);
    data[4*p + 1] = static_cast<float>(outpart_rdata(IPX,p));
    data[4*p + 2] = static_cast<float>(pm->multi_d ? outpart_rdata(IPY,p) :
                                       pm->mesh_size.x2min);
    data[4*p + 3] = static_cast<float>(pm->three_d ? outpart_rdata(IPZ,p) :
                                       pm->mesh_size.x3min);
    if (relativistic_mode &&
        (!RelativisticPpdSpeciesExactlyRepresentable(species) ||
         !std::isfinite(data[4*p]) || !std::isfinite(data[4*p + 1]) ||
         !std::isfinite(data[4*p + 2]) || !std::isfinite(data[4*p + 3]))) {
      relativistic_failure = 1;
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &relativistic_failure, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#endif
  if (relativistic_failure != 0) {
    delete[] data;
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "relativistic_hc ppd values must be finite and "
              << "identities exactly representable in float32" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::write);
  std::stringstream msg;
  msg << "# AthenaK particle position data at time= " << pm->time
      << "  nranks= " << global_variable::nranks
      << "  particles= " << pm->nprtcl_total
      << "  cycle=" << pm->ncycle << std::endl;
  if (relativistic_mode) {
    msg << "# metadata: schema=akcr_particle_output_v1 mode=relativistic_hc "
        << "units=code_model c_model=" << std::scientific << std::setprecision(17)
        << pp->c_model
        << " alpha_s=" << pp->alpha_s << std::endl
        << "# metadata: columns=species_x_y_z position_units=code_length "
        << "payload=float32" << std::endl;
  }
  if (global_variable::my_rank == 0) {
    partfile.Write_any_type(msg.str().c_str(),msg.str().size(),"byte");
  }
  std::size_t header_offset = msg.str().size();

  std::vector<int> rank_offset(global_variable::nranks,0);
  for (int n=1; n<global_variable::nranks; ++n) {
    rank_offset[n] = rank_offset[n-1] + pm->nprtcl_eachrank[n-1];
  }

  std::size_t myoffset = header_offset +
    static_cast<std::size_t>(4*rank_offset[global_variable::my_rank])*sizeof(float);
  if (partfile.Write_any_type_at_all(data,4*npout_thisrank,myoffset,"float") !=
      static_cast<std::size_t>(4*npout_thisrank)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle position output failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  partfile.Close();
  delete[] data;

  out_params.file_number++;
  out_params.last_time = (out_params.last_time < 0.0) ? pm->time :
                         out_params.last_time + out_params.dt;
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
