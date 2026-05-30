//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility> // make_pair
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "srcterms/turb_driver.hpp"
#include "particles/particles.hpp"
//#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  mkdir("rst",0775);
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rst/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir, 0775);
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb = pm->pmb_pack->nmb_thispack;

  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nadm=0, nz4c=0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  // if the spacetime is evolved, we do not need to checkpoint/recover the ADM variables
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  if (phydro != nullptr) {
    Kokkos::realloc(outarray_hyd, nmb, nhydro, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1+1);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2+1, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x3f, nmb, nout3+1, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (prad != nullptr) {
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pturb != nullptr) {
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    Kokkos::realloc(outarray_adm, nmb, nadm, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_adm, Kokkos::subview(padm->u_adm, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }

  // calculate max/min number of MeshBlocks across all ranks
  noutmbs_max = pm->nmb_eachrank[0];
  noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes everything to a single restart file

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nz4c=0, nadm=0, nco=0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
    nco = pz4c->ptracker.size();
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  bool single_file_per_rank = out_params.single_file_per_rank;
  std::string fname;
  if (single_file_per_rank) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    char rank_dir[20];
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("rst/") + std::string(rank_dir) + out_params.file_basename
      + number + ".rst";

    // Debugging output to check directory and filename
    // std::cout << "Rank " << global_variable::my_rank << " generated filename: "
    //           << fname << std::endl;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  // increment counters now so values for *next* dump are stored in restart file
  out_params.file_number++;
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetInteger(out_params.block_name, "file_number", out_params.file_number);
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write, single_file_per_rank);
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // output the input parameters (input file)
    resfile.Write_any_type(sbuf.c_str(), sbuf.size(), "byte", single_file_per_rank);

    // output Mesh information
    resfile.Write_any_type(&(pm->nmb_total), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->root_level), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_size), (sizeof(RegionSize)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mb_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->time), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->dt), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->ncycle), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(global_variable::nranks), (sizeof(int)), "byte",
                            single_file_per_rank);
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(pm->lloc_eachmb[0]),(pm->nmb_total)*sizeof(LogicalLocation),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->cost_eachmb[0]), (pm->nmb_total)*sizeof(float),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->rank_eachmb[0]), (pm->nmb_total)*sizeof(int),
                           "byte", single_file_per_rank);
    if (global_variable::nranks > 0) {
      resfile.Write_any_type(&(pm->gids_eachrank[0]),
                             (global_variable::nranks)*sizeof(int),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&(pm->nmb_eachrank[0]),
                             (global_variable::nranks)*sizeof(int),
                             "byte", single_file_per_rank);
    }
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // store z4c information
    if (pz4c != nullptr) {
      resfile.Write_any_type(&(pz4c->last_output_time), sizeof(Real), "byte",
                             single_file_per_rank);
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        resfile.Write_any_type(pt.GetPos(), 3*sizeof(Real), "byte",
                               single_file_per_rank);
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      resfile.Write_any_type(&(pturb->rstate), sizeof(RNG_State), "byte",
                             single_file_per_rank);
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  IOWrapperSizeT data_size = 0;
  if (phydro != nullptr) {
    data_size += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
  }
  if (pmhd != nullptr) {
    data_size += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    data_size += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    data_size += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    data_size += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
  }
  if (prad != nullptr) {
    data_size += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
  }
  if (pturb != nullptr) {
    data_size += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
  }
  if (pz4c != nullptr) {
    data_size += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
  } else if (padm != nullptr) {
    data_size += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
  }
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(data_size), sizeof(IOWrapperSizeT), "byte",
                            single_file_per_rank);
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 4*sizeof(int) + 2*sizeof(Real) +
                             sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  IOWrapperSizeT step2size = (pm->nmb_total)*(sizeof(LogicalLocation) + sizeof(float)
                           + sizeof(int))
                           + (global_variable::nranks)*2*sizeof(int);

  IOWrapperSizeT step3size = 3*nco*sizeof(Real);
  if (pz4c != nullptr) step3size += sizeof(Real);
  if (pturb != nullptr) step3size += sizeof(RNG_State);

  // write cell-centered variables in parallel
  IOWrapperSizeT offset_myrank = (step1size + step2size + step3size
                                  + sizeof(IOWrapperSizeT));

  if (!single_file_per_rank) {
    offset_myrank += data_size*(pm->gids_eachrank[global_variable::my_rank]);
  }

  IOWrapperSizeT myoffset = offset_myrank;

  // write cell-centered variables, one MeshBlock at a time (but parallelized over all
  // ranks). MeshBlocks are written seperately to reduce number of data elements per write
  // call, to avoid exceeding 2^31 limit for very large grids per MPI rank.
  if (phydro != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      }
    }
    offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nz4c*sizeof(Real); // z4c u0
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nadm*sizeof(Real); // adm u_adm
    myoffset = offset_myrank;
  }

  //--- STEP 5. Optional particle restart state.
  particles::Particles *ppart = pm->pmb_pack->ppart;
  if (ppart != nullptr) {
    constexpr std::uint64_t kPicMagic = 0x5049435253543031ULL;
    constexpr int kPicVersion = 1;

    const int nmb_local = pm->nmb_thisrank;
    const int gids_local = pm->gids_eachrank[global_variable::my_rank];
    const int npart = ppart->nprtcl_thispack;
    const int nrdata = ppart->nrdata;
    const int nidata = ppart->nidata;
    const int moment_cnt = particles::Particles::NMOM*nout3*nout2*nout1;
    const bool has_moments = ppart->deposit_moments && (ppart->moments.size() > 0);
    const bool requires_moments = ppart->couple_moments_to_mhd && ppart->deposit_moments;
    const bool has_edge = (ppart->couple_j_to_efield_representation ==
                           CoupledCurrentRepresentation::edge_staggered) &&
                          (ppart->j_edge_x1e.size() > 0);
    const int edge1_cnt = (nout3 + 1)*(nout2 + 1)*nout1;
    const int edge2_cnt = (nout3 + 1)*nout2*(nout1 + 1);
    const int edge3_cnt = nout3*(nout2 + 1)*(nout1 + 1);

    if (requires_moments && !has_moments) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart requested coupled moments, but moments are not "
                << "allocated."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
    auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_idata);

    std::vector<int> local_mb_counts(nmb_local, 0);
    for (int p=0; p<npart; ++p) {
      const int m = h_pi(PGID, p) - gids_local;
      if ((m < 0) || (m >= nmb_local)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Particle gid is not local at restart write time (gid="
                  << h_pi(PGID, p) << ", local gid range=[" << gids_local << ","
                  << (gids_local + nmb_local - 1) << "])." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ppart->particle_type == ParticleType::cosmic_ray) {
        const int sp = h_pi(PSP, p);
        if (sp < 0 || sp >= ppart->nspecies) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Particle species is out of range at restart write time "
                    << "(species=" << sp << ", nspecies=" << ppart->nspecies << ")."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      local_mb_counts[m] += 1;
    }

    std::vector<IOWrapperSizeT> local_mb_offsets(nmb_local + 1, 0);
    for (int m=0; m<nmb_local; ++m) {
      local_mb_offsets[m + 1] = local_mb_offsets[m] +
                                static_cast<IOWrapperSizeT>(local_mb_counts[m]);
    }

    std::vector<Real> packed_pr(npart*nrdata, 0.0);
    std::vector<int> packed_pi(npart*nidata, 0);
    std::vector<int> local_mb_cursor(nmb_local, 0);
    for (int p=0; p<npart; ++p) {
      const int m = h_pi(PGID, p) - gids_local;
      const int pinmb = local_mb_cursor[m]++;
      const IOWrapperSizeT lp = local_mb_offsets[m] + static_cast<IOWrapperSizeT>(pinmb);
      for (int n=0; n<nrdata; ++n) {
        packed_pr[lp*nrdata + n] = h_pr(n, p);
      }
      for (int n=0; n<nidata; ++n) {
        packed_pi[lp*nidata + n] = h_pi(n, p);
      }
    }

    std::vector<int> mb_counts_section;
    if (single_file_per_rank) {
      mb_counts_section = local_mb_counts;
    } else {
      mb_counts_section.assign(pm->nmb_total, 0);
      for (int m=0; m<nmb_local; ++m) {
        mb_counts_section[gids_local + m] = local_mb_counts[m];
      }
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, mb_counts_section.data(), pm->nmb_total,
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    }

    const int nmb_section = static_cast<int>(mb_counts_section.size());
    std::vector<IOWrapperSizeT> mb_offsets_section(nmb_section + 1, 0);
    for (int m=0; m<nmb_section; ++m) {
      mb_offsets_section[m + 1] = mb_offsets_section[m] +
                                  static_cast<IOWrapperSizeT>(mb_counts_section[m]);
    }
    const IOWrapperSizeT npart_section = mb_offsets_section[nmb_section];

    IOWrapperSizeT step5offset = step1size + step2size + step3size +
                                 sizeof(IOWrapperSizeT);
    if (single_file_per_rank) {
      step5offset += data_size*pm->nmb_thisrank;
    } else {
      step5offset += data_size*pm->nmb_total;
    }

    IOWrapperSizeT section_offset = step5offset;
    const IOWrapperSizeT magic_offset = section_offset;
    section_offset += sizeof(std::uint64_t);
    const IOWrapperSizeT version_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nmb_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nrdata_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nidata_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nout1_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nout2_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT nout3_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT has_mom_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT has_edge_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT mom_cnt_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT edge1_cnt_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT edge2_cnt_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT edge3_cnt_offset = section_offset;
    section_offset += sizeof(int);
    const IOWrapperSizeT npart_offset = section_offset;
    section_offset += sizeof(IOWrapperSizeT);
    const IOWrapperSizeT mb_count_offset = section_offset;
    section_offset += static_cast<IOWrapperSizeT>(nmb_section)*sizeof(int);
    const IOWrapperSizeT pr_real_offset = section_offset;
    section_offset += npart_section*nrdata*sizeof(Real);
    const IOWrapperSizeT pr_int_offset = section_offset;
    section_offset += npart_section*nidata*sizeof(int);
    const IOWrapperSizeT moments_offset = section_offset;
    if (has_moments) {
      section_offset += static_cast<IOWrapperSizeT>(nmb_section)*moment_cnt*sizeof(Real);
    }
    const IOWrapperSizeT edge1_offset = section_offset;
    if (has_edge) {
      section_offset += static_cast<IOWrapperSizeT>(nmb_section)*edge1_cnt*sizeof(Real);
    }
    const IOWrapperSizeT edge2_offset = section_offset;
    if (has_edge) {
      section_offset += static_cast<IOWrapperSizeT>(nmb_section)*edge2_cnt*sizeof(Real);
    }
    const IOWrapperSizeT edge3_offset = section_offset;
    if (has_edge) {
      section_offset += static_cast<IOWrapperSizeT>(nmb_section)*edge3_cnt*sizeof(Real);
    }

    if (global_variable::my_rank == 0 || single_file_per_rank) {
      int i_has_mom = has_moments ? 1 : 0;
      int i_has_edge = has_edge ? 1 : 0;
      if (resfile.Write_any_type_at(&kPicMagic, sizeof(std::uint64_t), magic_offset,
                                    "byte", single_file_per_rank)
            != sizeof(std::uint64_t) ||
          resfile.Write_any_type_at(&kPicVersion, 1, version_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nmb_section, 1, nmb_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nrdata, 1, nrdata_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nidata, 1, nidata_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout1, 1, nout1_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout2, 1, nout2_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&nout3, 1, nout3_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&i_has_mom, 1, has_mom_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&i_has_edge, 1, has_edge_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&moment_cnt, 1, mom_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge1_cnt, 1, edge1_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge2_cnt, 1, edge2_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&edge3_cnt, 1, edge3_cnt_offset, "int",
                                    single_file_per_rank) != 1 ||
          resfile.Write_any_type_at(&npart_section, sizeof(IOWrapperSizeT),
                                    npart_offset, "byte", single_file_per_rank)
            != sizeof(IOWrapperSizeT)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart metadata." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (nmb_section > 0) {
        if (resfile.Write_any_type_at(mb_counts_section.data(), nmb_section,
                                      mb_count_offset, "int",
                                      single_file_per_rank) != nmb_section) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle MeshBlock counts."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }

#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    for (int m=0; m<nmb_local; ++m) {
      const int cnt = local_mb_counts[m];
      if (cnt <= 0) continue;
      const int section_m = (single_file_per_rank ? m : gids_local + m);
      const IOWrapperSizeT gstart = mb_offsets_section[section_m];
      const IOWrapperSizeT lstart = local_mb_offsets[m];
      const IOWrapperSizeT pr_off = pr_real_offset + gstart*nrdata*sizeof(Real);
      const IOWrapperSizeT pi_off = pr_int_offset + gstart*nidata*sizeof(int);
      if (resfile.Write_any_type_at(&(packed_pr[lstart*nrdata]), cnt*nrdata, pr_off,
                                    "Real", single_file_per_rank) != cnt*nrdata) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart real data." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (resfile.Write_any_type_at(&(packed_pi[lstart*nidata]), cnt*nidata, pi_off,
                                    "int", single_file_per_rank) != cnt*nidata) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to write particle restart integer data." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    if (has_moments) {
      auto h_mom = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->moments);
      for (int m=0; m<nmb_local; ++m) {
        const int section_m = (single_file_per_rank ? m : gids_local + m);
        auto mom_mb = Kokkos::subview(h_mom, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                      Kokkos::ALL);
        const IOWrapperSizeT moff = moments_offset +
                                    static_cast<IOWrapperSizeT>(section_m)*moment_cnt*
                                    sizeof(Real);
        if (resfile.Write_any_type_at(mom_mb.data(), moment_cnt, moff, "Real",
                                      single_file_per_rank) != moment_cnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle moment restart data." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }

    if (has_edge) {
      auto h_x1e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x1e);
      auto h_x2e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x2e);
      auto h_x3e = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->j_edge_x3e);

      for (int m=0; m<nmb_local; ++m) {
        const int section_m = (single_file_per_rank ? m : gids_local + m);
        auto x1_mb = Kokkos::subview(h_x1e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto x2_mb = Kokkos::subview(h_x2e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto x3_mb = Kokkos::subview(h_x3e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        const IOWrapperSizeT x1off = edge1_offset +
                                     static_cast<IOWrapperSizeT>(section_m)*edge1_cnt*
                                     sizeof(Real);
        const IOWrapperSizeT x2off = edge2_offset +
                                     static_cast<IOWrapperSizeT>(section_m)*edge2_cnt*
                                     sizeof(Real);
        const IOWrapperSizeT x3off = edge3_offset +
                                     static_cast<IOWrapperSizeT>(section_m)*edge3_cnt*
                                     sizeof(Real);
        if (resfile.Write_any_type_at(x1_mb.data(), edge1_cnt, x1off, "Real",
                                      single_file_per_rank) != edge1_cnt ||
            resfile.Write_any_type_at(x2_mb.data(), edge2_cnt, x2off, "Real",
                                      single_file_per_rank) != edge2_cnt ||
            resfile.Write_any_type_at(x3_mb.data(), edge3_cnt, x3off, "Real",
                                      single_file_per_rank) != edge3_cnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to write particle edge-current restart data." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  // close file, clean up
  resfile.Close(single_file_per_rank);

  return;
}
