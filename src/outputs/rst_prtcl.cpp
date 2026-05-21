//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rst_prtcl.cpp
//! \brief writes per-rank particle restart dumps

#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

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
    std::size_t last_period = prst_fname.rfind('.');
    if (last_period != std::string::npos && last_period >= 5) {
      std::string outnumber_str = prst_fname.substr(last_period-5,5);
      out_params.file_number = std::stoi(outnumber_str) + 1;
      out_params.last_time = pm->time;
    }
  }
}

void ParticleRestartOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pp->nprtcl_thispack;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);
}

void ParticleRestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  char rank_dir[32];
  char number[8];
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
  std::string fname = "prst/" + std::string(rank_dir) + out_params.file_basename +
                      std::string(number) + ".prst";

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
  delete[] real_data;

  out_params.file_number++;
  out_params.last_time = (out_params.last_time < 0.0) ? pm->time :
                         out_params.last_time + out_params.dt;
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
}
