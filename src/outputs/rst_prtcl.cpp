#include <sys/stat.h>  // mkdir
#include <vector>

#include <algorithm>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor
// Checks compatibility options for VTK outputs

ParticleRestartOutput::ParticleRestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("prst",0775);
  // allocate arrays
  char rank_dir[20];
  std::snprintf(rank_dir, sizeof(rank_dir), "prst/rank_%08d/", global_variable::my_rank);
  mkdir(rank_dir, 0775);

  int prtcl_rst_flag = pin->GetOrAddInteger("problem","prtcl_rst_flag",0);
  if (prtcl_rst_flag) {
    std::string prst_fname = pin->GetString("problem","prtcl_res_file");
    size_t last_period = prst_fname.rfind('.');
    std::string outnumber_str = prst_fname.substr(last_period-5,5);
    int outnumber = std::stoi(outnumber_str);
    out_params.file_number = outnumber + 1; 
    out_params.last_time = pm->time;    
  }
}

//----------------------------------------------------------------------------------------
// ParticlePositionsOutput::LoadOutputData()
// Copies real and integer particle data to host for outputs

void ParticleRestartOutput::LoadOutputData(Mesh *pm) {
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pm->nprtcl_thisrank;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);

  // Create mirror view on device of host view of output particle real/int data
  auto d_outpart_rdata = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace(),
                                                    outpart_rdata);
  auto d_outpart_idata = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace(),
                                                    outpart_idata);
  // Copy particle positions into device mirrors
  Kokkos::deep_copy(d_outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(d_outpart_idata, pp->prtcl_idata);
  // Copy particle positions from device mirror to host output array
  Kokkos::deep_copy(outpart_rdata, d_outpart_rdata);
  Kokkos::deep_copy(outpart_idata, d_outpart_idata);
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlePositionsOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all particles and writes all particle positions to a binary file.

void ParticleRestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  int big_end = IsBigEndian(); // =1 on big endian machine

  // create filename: "ppd/file_basename"."file_id"."XXXXX".part.vtk
  // where XXXXX = 5-digit file_number
  std::string fname;
  char rank_dir[20];
  char number[7];
  std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);  
  std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
  fname = std::string("prst/") + std::string(rank_dir) + out_params.file_basename  + number + ".prst";
  FILE *pfile;
  if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    exit(EXIT_FAILURE);
  }
  std::fclose(pfile);


  IOWrapper partfile;
  std::size_t header_offset=0;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::write);

  // Write Part 5: Write particle data
  // allocate 1D vector of floats used to convert and output particle data
  //#int *int_data = new int[3*npout_thisrank];
  Real *real_data = new Real[17*npout_thisrank];
 


  // Loop over particles, load positions into data[]
  for (int p=0; p<npout_thisrank; ++p) {
    real_data[17*p] = static_cast<Real>(outpart_idata(PGID,p));
    real_data[17*p+1] = static_cast<Real>(outpart_idata(PTAG,p)); 
    real_data[17*p+2] = static_cast<Real>(outpart_idata(PSP,p));
    real_data[17*p+3] = outpart_rdata(IPX,p);
    real_data[17*p+4] = outpart_rdata(IPY,p);
    real_data[17*p+5] = outpart_rdata(IPZ,p);
    real_data[17*p+6] = outpart_rdata(IPVX,p);
    real_data[17*p+7] = outpart_rdata(IPVY,p);
    real_data[17*p+8] = outpart_rdata(IPVZ,p);
    real_data[17*p+9] = outpart_rdata(IPM,p);
    real_data[17*p+10] = outpart_rdata(IPBX,p);
    real_data[17*p+11] = outpart_rdata(IPBY,p);
    real_data[17*p+12] = outpart_rdata(IPBZ,p);
    real_data[17*p+13] = outpart_rdata(IPDX,p);
    real_data[17*p+14] = outpart_rdata(IPDY,p);
    real_data[17*p+15] = outpart_rdata(IPDZ,p);
    real_data[17*p+16] = outpart_rdata(IPDB,p);
  }

  Real npout_thisrank_r = static_cast<Real>(npout_thisrank);
  // output timestep/rank information
  partfile.Write_any_type(&(pm->time), (sizeof(Real)), "byte");
  partfile.Write_any_type(&(pm->dt), (sizeof(Real)), "byte");
  partfile.Write_any_type(&(npout_thisrank_r), (sizeof(Real)), "byte");

  // Write particle positions
  {
    std::size_t datasize = sizeof(int);
    std::size_t myoffset=header_offset +3*(sizeof(Real)); //+ (sizeof(int)) ; //+ 3*npout_thisrank*datasize;
    // collective writes for minimum number of particles across ranks
    //if (partfile.Write_any_type_at_all(&(int_data[0]),3*npout_thisrank,myoffset,"int")
    //      != 3*npout_thisrank) {
    //  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    //      << std::endl << "particle data not written correctly to particle restart file, "
    //      << "file is broken." << std::endl;
    //  exit(EXIT_FAILURE);
    //}
    //myoffset += 3*npout_thisrank*datasize;
    // individual writes for remaining particles on each rank
    if (partfile.Write_any_type_at_all(&(real_data[0]),17*npout_thisrank,myoffset,"Real")
          != 17*npout_thisrank) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to particle restart file, "
          << "file is broken." << std::endl;
      exit(EXIT_FAILURE);
    }    
  }

  // close the output file and clean up
  partfile.Close();
  delete[] real_data;
  //delete[] int_data;
  // increment counters
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  return;
}
