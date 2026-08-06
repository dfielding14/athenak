//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.cpp
//! \brief Implementation of constructors and functions in class ProblemGenerator.
//! Default constructor calls problem generator function, while  constructor for restarts
//! reads data from restart file, as well as re-initializing problem-specific data.

#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "athena.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "particles/particles.hpp"
#include "restart_layout.hpp"
#include "restart_manifest.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen.hpp"

namespace {

IOWrapperSizeT CheckedRestartReadAdd(IOWrapperSizeT left, IOWrapperSizeT right,
                                     const std::string &context) {
  return restart_layout::CheckedAdd(left, right, FailNodeRestart, context);
}

IOWrapperSizeT CheckedRestartReadCount(int value, const std::string &context) {
  return restart_layout::CheckedNonNegative(value, FailNodeRestart, context);
}

int CheckedRestartReadInt(IOWrapperSizeT value, const std::string &context) {
  if (value > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    FailNodeRestart(context + " exceeds INT_MAX.");
  }
  return static_cast<int>(value);
}

int CheckedRestartReadIntAdd(int left, int right, const std::string &context) {
  return CheckedRestartReadInt(
      CheckedRestartReadAdd(CheckedRestartReadCount(left, context),
                            CheckedRestartReadCount(right, context), context),
      context);
}

}  // namespace


//----------------------------------------------------------------------------------------
// default constructor, calls pgen function.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm) :
    user_bcs(false),
    user_srcs(false),
    user_dt(false),
    user_hist(false),
    pmy_mesh_(pm) {
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }

  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_dt = pin->GetOrAddBoolean("problem","user_dt",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);

  // second argument false since this IS NOT a restart
  CallProblemGenerator(pin, false);

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User BCs specified in <mesh> block, but not enrolled "
                << "by SetProblemData()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr && user_stage_srcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User SRCs specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined timestep outputs were enrolled if needed
  if (user_dt) {
    if (user_time_step_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User timestep specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User history output specified in <problem> block, but "
                << "not enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
// constructor for restarts
// When called, data needed to rebuild mesh has been read from restart file by
// Mesh::BuildTreeFromRestart() function. This constructor reads from the restart file and
// initializes all the dependent variables (u0,b0,etc) stored in each Physics class. It
// also calls ProblemGenerator::SetProblemData() function to set any user-defined BCs,
// and any data necessary for restart runs to continue correctly.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm, IOWrapper resfile,
                                   bool single_file_per_rank,
                                   const NodeRestartManifest *node_restart_manifest) :
    user_bcs(false),
    user_srcs(false),
    user_dt(false),
    user_hist(false),
    pmy_mesh_(pm) {
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }
  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_dt = pin->GetOrAddBoolean("problem","user_dt",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);
  if (pm->pmb_pack->ppart != nullptr &&
      pm->pmb_pack->ppart->IsLagrangianMC() &&
      (node_restart_manifest != nullptr ||
       (global_variable::nranks > 1 && !single_file_per_rank))) {
    FailNodeRestart(
        "lagrangian_mc particle restarts do not yet support node-sharded files; "
        "MPI runs require single_file_per_rank=true.");
  }

  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedRestartReadInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx1, indcs.ng, true, FailNodeRestart, "restart x1 extent"),
      "restart x1 extent");
  int nout2 = CheckedRestartReadInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx2, indcs.ng, indcs.nx2 > 1, FailNodeRestart, "restart x2 extent"),
      "restart x2 extent");
  int nout3 = CheckedRestartReadInt(restart_layout::CheckedExtentWithGhosts(
      indcs.nx3, indcs.ng, indcs.nx3 > 1, FailNodeRestart, "restart x3 extent"),
      "restart x3 extent");
  int nmb = pm->pmb_pack->nmb_thispack;
  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad=pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nrad = 0, nhydro = 0, nmhd = 0, nforce = 3, nadm = 0, nz4c = 0;
  if (phydro != nullptr) {
    nhydro = CheckedRestartReadIntAdd(phydro->nhydro, phydro->nscalars,
                                      "restart Hydro component count");
  }
  if (pmhd != nullptr) {
    nmhd = CheckedRestartReadIntAdd(pmhd->nmhd, pmhd->nscalars,
                                    "restart MHD component count");
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }

  // root process reads z4c last_output_time and tracker data
  if (pz4c != nullptr) {
    Real last_output_time;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_Reals(&last_output_time, 1,single_file_per_rank) != 1) {
        FailNodeRestart("z4c::last_output_time data size read from restart file is "
                        "incorrect, restart file is broken.");
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      io_wrapper::BroadcastBytes(&last_output_time, sizeof(Real), 0,
                                 MPI_COMM_WORLD);
    }
#endif
    pz4c->last_output_time = last_output_time;

    for (auto &pt : pz4c->ptracker) {
      Real pos[3];
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_Reals(&pos[0], 3, single_file_per_rank) != 3) {
          FailNodeRestart("compact object tracker data size read from restart file is "
                          "incorrect, restart file is broken.");
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        io_wrapper::BroadcastBytes(&pos[0], 3*sizeof(Real), 0, MPI_COMM_WORLD);
      }
#endif
      pt->SetPos(&pos[0]);
    }
  }

  if (pturb != nullptr) {
    // root process reads size the random seed
    char *rng_data = new char[sizeof(RNG_State)];
    // the master process reads the variables data
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(rng_data, 1, sizeof(RNG_State), single_file_per_rank)
          != sizeof(RNG_State)) {
        FailNodeRestart("RNG data size read from restart file is incorrect, restart "
                        "file is broken.");
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      // then broadcast the RNG information
      io_wrapper::BroadcastBytes(rng_data, sizeof(RNG_State), 0,
                                 MPI_COMM_WORLD);
    }
#endif
    std::memcpy(&(pturb->rstate), &(rng_data[0]), sizeof(RNG_State));
  }

  // root process reads size of CC and FC data arrays from restart file
  IOWrapperSizeT variablesize = sizeof(IOWrapperSizeT);
  char *variabledata = new char[variablesize];
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
        != variablesize) {
      FailNodeRestart("Variable data size read from restart file is incorrect, restart "
                      "file is broken.");
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the datasize information
  if (!single_file_per_rank) {
    io_wrapper::BroadcastBytes(variabledata, variablesize, 0, MPI_COMM_WORLD);
  }
#endif
  IOWrapperSizeT data_size;
  std::memcpy(&data_size, &(variabledata[0]), sizeof(IOWrapperSizeT));

  // calculate total number of CC variables
  IOWrapperSizeT headeroffset;
  // master process gets file offset
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    headeroffset = resfile.GetPosition(single_file_per_rank);
  }
#if MPI_PARALLEL_ENABLED
  // then broadcasts it
  if (!single_file_per_rank) {
    io_wrapper::BroadcastBytes(&headeroffset, sizeof(IOWrapperSizeT), 0,
                               MPI_COMM_WORLD);
  }
#endif

  restart_layout::PayloadLayout layout = restart_layout::PayloadLayout::Build(
      {CheckedRestartReadCount(nout1, "restart x1 extent"),
       CheckedRestartReadCount(nout2, "restart x2 extent"),
       CheckedRestartReadCount(nout3, "restart x3 extent"),
       CheckedRestartReadCount(nhydro, "restart Hydro component count"),
       CheckedRestartReadCount(nmhd, "restart MHD component count"),
       CheckedRestartReadCount(nrad, "restart radiation component count"),
       pturb == nullptr ? 0 : CheckedRestartReadCount(nforce, "restart forcing count"),
       CheckedRestartReadCount(nz4c, "restart Z4c component count"),
       CheckedRestartReadCount(nadm, "restart ADM component count"),
       sizeof(Real)}, FailNodeRestart);
  IOWrapperSizeT data_size_ = layout.data_bytes;
  IOWrapperSizeT local_blocks =
      CheckedRestartReadCount(nmb, "restart local MeshBlock count");
  const auto preflight_local = [&](IOWrapperSizeT block_bytes,
                                   const std::string &context) {
    restart_layout::CheckedMemorySize(
        restart_layout::CheckedMultiply(local_blocks, block_bytes, FailNodeRestart,
                                        context),
        FailNodeRestart, context);
  };
  preflight_local(layout.hydro_bytes, "restart local Hydro allocation bytes");
  preflight_local(layout.mhd_bytes, "restart local MHD allocation bytes");
  preflight_local(layout.mhd_x1f_bytes, "restart local MHD x1-face allocation bytes");
  preflight_local(layout.mhd_x2f_bytes, "restart local MHD x2-face allocation bytes");
  preflight_local(layout.mhd_x3f_bytes, "restart local MHD x3-face allocation bytes");
  preflight_local(layout.radiation_bytes, "restart local radiation allocation bytes");
  preflight_local(layout.forcing_bytes, "restart local forcing allocation bytes");
  preflight_local(layout.z4c_bytes, "restart local Z4c allocation bytes");
  preflight_local(layout.adm_bytes, "restart local ADM allocation bytes");
  int nout1f = nout1;
  int nout2f = nout2;
  int nout3f = nout3;
  if (pmhd != nullptr) {
    nout1f = CheckedRestartReadInt(
        CheckedRestartReadAdd(nout1, 1, "restart MHD x1 extent"),
        "restart MHD x1 extent");
    nout2f = CheckedRestartReadInt(
        CheckedRestartReadAdd(nout2, 1, "restart MHD x2 extent"),
        "restart MHD x2 extent");
    nout3f = CheckedRestartReadInt(
        CheckedRestartReadAdd(nout3, 1, "restart MHD x3 extent"),
        "restart MHD x3 extent");
  }

  if (data_size_ != data_size) {
    FailNodeRestart("CC data size read from restart file not equal to size of Hydro, "
                    "MHD, Rad, and/or Z4c arrays, restart file is broken.");
  }

  // read CC data into host array
  int mygids = pm->gids_eachrank[global_variable::my_rank];
  IOWrapperSizeT offset_myrank = headeroffset;
  if (!single_file_per_rank) {
    offset_myrank = restart_layout::CheckedOffset(
        headeroffset, data_size_,
        CheckedRestartReadCount(pm->gids_eachrank[global_variable::my_rank],
                                "restart shared rank block offset"),
        FailNodeRestart, "restart shared rank byte offset");
  }
  IOWrapperSizeT node_restart_virtual_base = offset_myrank;
  IOWrapperSizeT myoffset = offset_myrank;

  std::vector<char> node_restart_blocks;
  if (node_restart_manifest != nullptr) {
    if (single_file_per_rank) {
      FailNodeRestart("node restart manifests cannot use rank-sharded restart mode.");
    }
    if (pm->nmb_total != node_restart_manifest->NumMeshBlocks() ||
        headeroffset != node_restart_manifest->HeaderSize()) {
      FailNodeRestart("canonical payload header does not match the manifest inventory.");
    }
    node_restart_manifest->LoadLocalBlocks(mygids, nmb, data_size, &node_restart_blocks);
  }
  auto read_restart_reals_at =
      [&](void *buffer, IOWrapperSizeT count, IOWrapperSizeT offset, bool collective) {
    if (node_restart_manifest == nullptr) {
      return collective
          ? resfile.Read_Reals_at_all(buffer, count, offset, single_file_per_rank)
          : resfile.Read_Reals_at(buffer, count, offset, single_file_per_rank);
    }
    IOWrapperSizeT bytes = restart_layout::CheckedMultiply(
        count, sizeof(Real), FailNodeRestart, "local restart field byte count");
    if (offset < node_restart_virtual_base ||
        offset - node_restart_virtual_base > node_restart_blocks.size() ||
        bytes > node_restart_blocks.size() - (offset - node_restart_virtual_base)) {
      FailNodeRestart("local field range is outside the routed restart buffer.");
    }
    if (bytes > 0) {
      std::memcpy(buffer,
                  &(node_restart_blocks[offset - node_restart_virtual_base]), bytes);
    }
    return static_cast<std::size_t>(count);
  };

  HostArray5D<Real> ccin("rst-cc-in", 1, 1, 1, 1, 1);
  HostFaceFld4D<Real> fcin("rst-fc-in", 1, 1, 1, 1);

  // calculate max/min number of MeshBlocks across all ranks
  int noutmbs_max = pm->nmb_eachrank[0];
  int noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }

  if (phydro != nullptr) {
    Kokkos::realloc(ccin, nmb, nhydro, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart Hydro subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true)
            != mbcnt) {
          FailNodeRestart("CC hydro data not read correctly from rst file, restart "
                          "file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart Hydro subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false)
            != mbcnt) {
          FailNodeRestart("CC hydro data not read correctly from rst file, restart "
                          "file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.hydro_bytes,
                                          "restart field offset"); // hydro u0
    myoffset = offset_myrank;
  }

  if (pmhd != nullptr) {
    Kokkos::realloc(ccin, nmb, nmhd, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart MHD subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true)
            != mbcnt) {
          FailNodeRestart("CC mhd data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart MHD subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false)
            != mbcnt) {
          FailNodeRestart("CC mhd data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.mhd_bytes,
                                          "restart field offset"); // mhd u0
    myoffset = offset_myrank;

    Kokkos::realloc(fcin.x1f, nmb, nout3, nout2, nout1f);
    Kokkos::realloc(fcin.x2f, nmb, nout3, nout2f, nout1);
    Kokkos::realloc(fcin.x3f, nmb, nout3f, nout2, nout1);
    // read FC data into host array, again one MeshBlock at a time
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT fldcnt = restart_layout::CheckedSizeT(
            x1fptr.size(), FailNodeRestart, "restart MHD x1-face subview count");

        if (read_restart_reals_at(x1fptr.data(), fldcnt, myoffset, true) != fldcnt) {
          FailNodeRestart("Input b0.x1f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x1f_bytes,
                                         "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestart, "restart MHD x2-face subview count");

        if (read_restart_reals_at(x2fptr.data(), fldcnt, myoffset, true) != fldcnt) {
          FailNodeRestart("Input b0.x2f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x2f_bytes,
                                         "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestart, "restart MHD x3-face subview count");

        if (read_restart_reals_at(x3fptr.data(), fldcnt, myoffset, true) != fldcnt) {
          FailNodeRestart("Input b0.x3f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x3f_bytes,
                                         "restart MHD face offset");

        myoffset = CheckedRestartReadAdd(
            myoffset,
            layout.mhd_face_stride_remainder_bytes,
            "restart MeshBlock offset");
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT fldcnt = restart_layout::CheckedSizeT(
            x1fptr.size(), FailNodeRestart, "restart MHD x1-face subview count");

        if (read_restart_reals_at(x1fptr.data(), fldcnt, myoffset, false) != fldcnt) {
          FailNodeRestart("Input b0.x1f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x1f_bytes,
                                         "restart MHD face offset");

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x2fptr.size(), FailNodeRestart, "restart MHD x2-face subview count");

        if (read_restart_reals_at(x2fptr.data(), fldcnt, myoffset, false) != fldcnt) {
          FailNodeRestart("Input b0.x2f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x2f_bytes,
                                         "restart MHD face offset");

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = restart_layout::CheckedSizeT(
            x3fptr.size(), FailNodeRestart, "restart MHD x3-face subview count");

        if (read_restart_reals_at(x3fptr.data(), fldcnt, myoffset, false) != fldcnt) {
          FailNodeRestart("Input b0.x3f field not read correctly from rst file, "
                          "restart file is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, layout.mhd_x3f_bytes,
                                         "restart MHD face offset");

        myoffset = CheckedRestartReadAdd(
            myoffset,
            layout.mhd_face_stride_remainder_bytes,
            "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x1f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x2f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x3f);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.mhd_x1f_bytes,
                                          "restart field offset"); // mhd b0.x1f
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.mhd_x2f_bytes,
                                          "restart field offset"); // mhd b0.x2f
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.mhd_x3f_bytes,
                                          "restart field offset"); // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    Kokkos::realloc(ccin, nmb, nrad, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart radiation subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true) != mbcnt) {
          FailNodeRestart("CC rad data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart radiation subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false) != mbcnt) {
          FailNodeRestart("CC rad data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(prad->i0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.radiation_bytes,
                                          "restart field offset"); // radiation i0
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    Kokkos::realloc(ccin, nmb, nforce, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart forcing subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true) != mbcnt) {
          FailNodeRestart("CC turb data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart forcing subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false) != mbcnt) {
          FailNodeRestart("CC turb data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pturb->force, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.forcing_bytes,
                                          "restart field offset"); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    Kokkos::realloc(ccin, nmb, nz4c, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart Z4c subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true) != mbcnt) {
          FailNodeRestart("CC z4c data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart Z4c subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false) != mbcnt) {
          FailNodeRestart("CC z4c data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pz4c->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.z4c_bytes,
                                          "restart field offset"); // z4c u0
    myoffset = offset_myrank;

    // We also need to reinitialize the ADM data.
    pz4c->Z4cToADM(pmy_mesh_->pmb_pack);
  } else if (padm != nullptr) {
    Kokkos::realloc(ccin, nmb, nadm, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart ADM subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, true) != mbcnt) {
          FailNodeRestart("CC adm data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT mbcnt = restart_layout::CheckedSizeT(
            mbptr.size(), FailNodeRestart, "restart ADM subview count");
        if (read_restart_reals_at(mbptr.data(), mbcnt, myoffset, false) != mbcnt) {
          FailNodeRestart("CC adm data not read correctly from rst file, restart file "
                          "is broken.");
        }
        myoffset = CheckedRestartReadAdd(myoffset, data_size, "restart MeshBlock offset");
      }
    }
    Kokkos::deep_copy(Kokkos::subview(padm->u_adm, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    offset_myrank = CheckedRestartReadAdd(offset_myrank, layout.adm_bytes,
                                          "restart field offset"); // adm u_adm
    myoffset = offset_myrank;
  }

  if (pm->pmb_pack->ppart != nullptr && pm->pmb_pack->ppart->IsLagrangianMC()) {
    pm->pmb_pack->ppart->ReadRestartData(resfile, single_file_per_rank);
  }

  // call problem generator again to re-initialize data, fn ptrs, as needed
  // second argument true since this IS a restart
  CallProblemGenerator(pin, true);

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      FailNodeRestart("User BCs specified in <mesh> block, but not enrolled "
                      "during restart by SetProblemData().");
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr && user_stage_srcs_func == nullptr) {
      FailNodeRestart("User SRCs specified in <problem> block, but not "
                      "enrolled by UserProblem().");
    }
  }
  // Check that user defined timestep outputs were enrolled if needed
  if (user_dt) {
    if (user_time_step_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User timestep specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      FailNodeRestart("User history output specified in <problem> block, "
                      "but not enrolled by UserProblem().");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::OutputErrors()
//! \brief Generic function for computing the L1 and L-infty difference between solutions
//! stored in the u0 and u1 registers, and outputting them to an error file.  This is
//! used for linear wave convergence tests, for example.
//! Function requires appropriate solutions already stored in u0 and u1.

void ProblemGenerator::OutputErrors(ParameterInput *pin, Mesh *pm) {
  Real l1_err[16];
  Real linfty_err=0.0;
  int nvars=0,nprev=0;

  // capture class variables for kernel
  auto &indcs = pm->mb_indcs;
  int &nx1 = indcs.nx1;
  int &nx2 = indcs.nx2;
  int &nx3 = indcs.nx3;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

  // compute errors for Hydro  -----------------------------------------------------------
  if (pmbp->phydro != nullptr) {
    nvars = pmbp->phydro->nhydro;

    auto &is_ideal_ = pmbp->phydro->peos->eos_data.is_ideal;
    auto &u0_ = pmbp->phydro->u0;
    auto &u1_ = pmbp->phydro->u1;

    const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    array_sum::GlobalSum sum_this_mb;
    Kokkos::parallel_reduce("L1-err",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum, Real &max_err) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

      // conserved variables:
      array_sum::GlobalSum evars;
      evars.the_array[IDN] = vol*fabs(u0_(m,IDN,k,j,i) - u1_(m,IDN,k,j,i));
      max_err = fmax(max_err, evars.the_array[IDN]);
      evars.the_array[IM1] = vol*fabs(u0_(m,IM1,k,j,i) - u1_(m,IM1,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM1]);
      evars.the_array[IM2] = vol*fabs(u0_(m,IM2,k,j,i) - u1_(m,IM2,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM2]);
      evars.the_array[IM3] = vol*fabs(u0_(m,IM3,k,j,i) - u1_(m,IM3,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM3]);
      if (is_ideal_) {
        evars.the_array[IEN] = vol*fabs(u0_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i));
        max_err = fmax(max_err, evars.the_array[IEN]);
      }

      // fill rest of the_array with zeros, if narray < NREDUCTION_VARIABLES
      for (int n=nvars; n<NREDUCTION_VARIABLES; ++n) {
        evars.the_array[n] = 0.0;
      }

      // sum into parallel reduce
      mb_sum += evars;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb), Kokkos::Max<Real>(linfty_err));

    // store data into l1_err array
    for (int n=0; n<nvars; ++n) {
      l1_err[n] = sum_this_mb.the_array[n];
    }
    nprev += nvars;
  }

  // compute errors for MHD  -------------------------------------------------------------
  if (pmbp->pmhd != nullptr) {
    nvars = pmbp->pmhd->nmhd + 3;  // include 3-compts of cell-centered B in errors
    auto &is_ideal_ = pmbp->pmhd->peos->eos_data.is_ideal;

    int bindx;
    if (is_ideal_) {
      bindx = 5;
    } else {
      bindx = 4;
    }

    auto &u0_ = pmbp->pmhd->u0;
    auto &u1_ = pmbp->pmhd->u1;
    auto &b0_ = pmbp->pmhd->b0;
    auto &b1_ = pmbp->pmhd->b1;

    const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    array_sum::GlobalSum sum_this_mb;
    Kokkos::parallel_reduce("L1-err-Sums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum, Real &max_err) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

      // conserved variables:
      array_sum::GlobalSum evars;
      evars.the_array[IDN] = vol*fabs(u0_(m,IDN,k,j,i) - u1_(m,IDN,k,j,i));
      max_err = fmax(max_err, evars.the_array[IDN]);
      evars.the_array[IM1] = vol*fabs(u0_(m,IM1,k,j,i) - u1_(m,IM1,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM1]);
      evars.the_array[IM2] = vol*fabs(u0_(m,IM2,k,j,i) - u1_(m,IM2,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM2]);
      evars.the_array[IM3] = vol*fabs(u0_(m,IM3,k,j,i) - u1_(m,IM3,k,j,i));
      max_err = fmax(max_err, evars.the_array[IM3]);
      if (is_ideal_) {
        evars.the_array[IEN] = vol*fabs(u0_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i));
        max_err = fmax(max_err, evars.the_array[IEN]);
      }

      // cell-centered B
      Real bcc0 = 0.5*(b0_.x1f(m,k,j,i) + b0_.x1f(m,k,j,i+1));
      Real bcc1 = 0.5*(b1_.x1f(m,k,j,i) + b1_.x1f(m,k,j,i+1));
      evars.the_array[bindx] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+1]);

      bcc0 = 0.5*(b0_.x2f(m,k,j,i) + b0_.x2f(m,k,j+1,i));
      bcc1 = 0.5*(b1_.x2f(m,k,j,i) + b1_.x2f(m,k,j+1,i));
      evars.the_array[bindx+1] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+2]);

      bcc0 = 0.5*(b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j,i));
      bcc1 = 0.5*(b1_.x3f(m,k,j,i) + b1_.x3f(m,k+1,j,i));
      evars.the_array[bindx+2] = vol*fabs(bcc0 - bcc1);
      max_err = fmax(max_err, evars.the_array[IEN+3]);

      // fill rest of the_array with zeros, if narray < NREDUCTION_VARIABLES
      for (int n=nvars; n<NREDUCTION_VARIABLES; ++n) {
        evars.the_array[n] = 0.0;
      }

      // sum into parallel reduce
      mb_sum += evars;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb), Kokkos::Max<Real>(linfty_err));

    // store data into l1_err array
    for (int n=0; n<nvars; ++n) {
      l1_err[n+nprev] = sum_this_mb.the_array[n];
    }
    nprev += nvars;
  }

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &l1_err, nprev, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &linfty_err, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

  // normalize errors by number of cells
  Real vol=  (pmbp->pmesh->mesh_size.x1max - pmbp->pmesh->mesh_size.x1min)
            *(pmbp->pmesh->mesh_size.x2max - pmbp->pmesh->mesh_size.x2min)
            *(pmbp->pmesh->mesh_size.x3max - pmbp->pmesh->mesh_size.x3min);
  for (int i=0; i<nprev; ++i) l1_err[i] = l1_err[i]/vol;
  linfty_err /= vol;

  // compute rms error
  Real rms_err = 0.0;
  for (int i=0; i<nprev; ++i) {
    rms_err += SQR(l1_err[i]);
  }
  rms_err = std::sqrt(rms_err);

  // root process opens output file and writes out errors
  if (global_variable::my_rank == 0) {
    std::string fname;
    fname.assign(pin->GetString("job","basename"));
    fname.append("-errs.dat");
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" <<std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::fprintf(pfile, "# Nx1  Nx2  Nx3   Ncycle   RMS-L1       L-infty       ");
      if (pmbp->phydro != nullptr) {
        std::fprintf(pfile,"d_L1          M1_L1         M2_L1         M3_L1         ");
        if (pmbp->phydro->peos->eos_data.is_ideal) {
          std::fprintf(pfile,"E_L1          ");
        }
      }
      if (pmbp->pmhd != nullptr) {
        std::fprintf(pfile,"d_L1          M1_L1         M2_L1         M3_L1         ");
        if (pmbp->pmhd->peos->eos_data.is_ideal) {
          std::fprintf(pfile,"E_L1          ");
        }
        std::fprintf(pfile,"B1_L1         B2_L1         B3_L1");
      }
      std::fprintf(pfile, "\n");
    }

    // write errors
    std::fprintf(pfile, "%04d", pmbp->pmesh->mesh_indcs.nx1);
    std::fprintf(pfile, "  %04d", pmbp->pmesh->mesh_indcs.nx2);
    std::fprintf(pfile, "  %04d", pmbp->pmesh->mesh_indcs.nx3);
    std::fprintf(pfile, "  %05d  %e %e", pmbp->pmesh->ncycle, rms_err, linfty_err);
    for (int i=0; i<nprev; ++i) {
      std::fprintf(pfile, "  %e", l1_err[i]);
    }
    std::fprintf(pfile, "\n");
    std::fclose(pfile);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::CallProblemGenerator()
//! \brief selects one of the default problem generators compiled automatically with
//! the source code depending on input string in <problem> block ELSE selects a
//! user-defined problem generator function compiled with the code.

void ProblemGenerator::CallProblemGenerator(ParameterInput *pin, bool is_restart) {
#if USER_PROBLEM_ENABLED
  // call user-defined problem generator (if USER_PROBLEM_ENABLED macro defined at build)
  UserProblem(pin, is_restart);
#else
  // else read name of built-in pgen from <problem> block in input file, and call
  std::string pgen_fun_name = pin->GetOrAddString("problem", "pgen_name", "none");

  if (pgen_fun_name.compare("advection") == 0) {
    Advection(pin, is_restart);
  } else if (pgen_fun_name.compare("cpaw") == 0) {
    AlfvenWave(pin, is_restart);
  } else if (pgen_fun_name.compare("gr_bondi") == 0) {
    BondiAccretion(pin, is_restart);
  } else if (pgen_fun_name.compare("cshock") == 0) {
    CShock(pin, is_restart);
  } else if (pgen_fun_name.compare("divb_amr") == 0) {
    DivBAMR(pin, is_restart);
  } else if (pgen_fun_name.compare("linear_wave") == 0) {
    LinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("implode") == 0) {
    LWImplode(pin, is_restart);
  } else if (pgen_fun_name.compare("gr_monopole") == 0) {
    Monopole(pin, is_restart);
  } else if (pgen_fun_name.compare("mri3d") == 0) {
    MRI3d(pin, is_restart);
  } else if (pgen_fun_name.compare("orszag_tang") == 0) {
    OrszagTang(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_linear_wave") == 0) {
    RadiationLinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("rad_beam") == 0) {
    RadiationBeam(pin, is_restart);
  } else if (pgen_fun_name.compare("shock_tube") == 0) {
    ShockTube(pin, is_restart);
  } else if (pgen_fun_name.compare("shwave") == 0) {
    Shwave(pin, is_restart);
  } else if (pgen_fun_name.compare("z4c_boosted_puncture") == 0) {
    Z4cBoostedPuncture(pin, is_restart);
  } else if (pgen_fun_name.compare("z4c_linear_wave") == 0) {
    Z4cLinearWave(pin, is_restart);
  } else if (pgen_fun_name.compare("spherical_collapse") == 0) {
    SphericalCollapse(pin, is_restart);
  } else if (pgen_fun_name.compare("diffusion") == 0) {
    Diffusion(pin, is_restart);
  // else, name not set on command line or input file, print warning and quit
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Problem generator name could not be found in <problem> block in input file"
        << std::endl
        << "and it was not set by -D PROBLEM option on cmake command line during build"
        << std::endl
        << "Rerun cmake with -D PROBLEM=file to specify custom problem generator file"
        << std::endl;;
    std::exit(EXIT_FAILURE);
  }
#endif
  return;
}
