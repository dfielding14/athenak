//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file track_prtcl.cpp
//! \brief writes data for tracked particles in unformatted binary

#include <sys/stat.h>  // mkdir
#include <vector>

#include <algorithm>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "outputs.hpp"

namespace {
constexpr int NTRACK_VALUES = 14;
const char *TRACK_VARIABLES =
    "x y z vx vy vz rho press temp eint scalar0 gid level active";
} // namespace

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

TrackedParticleOutput::TrackedParticleOutput(ParameterInput *pin, Mesh *pm,
                                             OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create new directory for this output. Comments in binary.cpp constructor explain why
  mkdir("trk",0775);
  // allocate arrays
  npout_eachrank.resize(global_variable::nranks);
  if (pin->DoesParameterExist(op.block_name,"nparticles")) {
    ntrack = pin->GetInteger(op.block_name,"nparticles");
  } else if (pin->DoesParameterExist(op.block_name,"ntrack")) {
    ntrack = pin->GetInteger(op.block_name,"ntrack");
    if (global_variable::my_rank == 0) {
      std::cout << "### WARNING in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << op.block_name << ">/ntrack is deprecated; use "
                << "nparticles instead." << std::endl;
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle output requires 'nparticles'"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack;
}

//----------------------------------------------------------------------------------------
// TrackedParticleOutput::LoadOutputData()
// Copies data for tracked particles on this rank to host outpart array

void TrackedParticleOutput::LoadOutputData(Mesh *pm) {
  pm->CountParticles();

  // Load data for tracked particles on this rank into new device array
  DualArray1D<TrackedParticleData> tracked_prtcl("d_trked",ntrack_thisrank);
  int npart = pm->nprtcl_thisrank;
  auto &pr = pm->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pm->pmb_pack->ppart->prtcl_idata;
  int ntrack_ = ntrack;
  bool is_lagrangian_mc =
      (pm->pmb_pack->ppart->particle_type == ParticleType::lagrangian_mc);
  int gids = pm->pmb_pack->gids;
  int nmb = pm->pmb_pack->nmb_thispack;
  auto &mblev = pm->pmb_pack->pmb->mb_lev;
  auto &mbsize = pm->pmb_pack->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  bool multi_d = pm->multi_d;
  bool three_d = pm->three_d;
  DvceArray5D<Real> w0;
  int nfluid = 0;
  int nscalars = 0;
  EOS_Data eos;
  if (pm->pmb_pack->phydro != nullptr) {
    w0 = pm->pmb_pack->phydro->w0;
    nfluid = pm->pmb_pack->phydro->nhydro;
    nscalars = pm->pmb_pack->phydro->nscalars;
    eos = pm->pmb_pack->phydro->peos->eos_data;
  } else {
    w0 = pm->pmb_pack->pmhd->w0;
    nfluid = pm->pmb_pack->pmhd->nmhd;
    nscalars = pm->pmb_pack->pmhd->nscalars;
    eos = pm->pmb_pack->pmhd->peos->eos_data;
  }
  const Real missing = -std::numeric_limits<Real>::max();

  // Create device-side counter
  Kokkos::View<int> d_counter("counter");
  Kokkos::deep_copy(d_counter, 0);

  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    if (pi(PTAG,p) < ntrack_) {
      int index = Kokkos::atomic_fetch_add(&d_counter(),1);
      int m = pi(PGID,p) - gids;
      bool active = (m >= 0 && m < nmb);
      if (is_lagrangian_mc && pi(PLASTMOVE,p) < 0) {
        active = false;
      }
      Real x = is_lagrangian_mc ? pr(LMCX,p) : pr(IPX,p);
      Real y = is_lagrangian_mc ? pr(LMCY,p) : pr(IPY,p);
      Real z = is_lagrangian_mc ? pr(LMCZ,p) : pr(IPZ,p);
      tracked_prtcl.d_view(index).tag = pi(PTAG,p);
      tracked_prtcl.d_view(index).gid = pi(PGID,p);
      tracked_prtcl.d_view(index).level = active ? mblev.d_view(m) : -1;
      tracked_prtcl.d_view(index).active = active ? 1 : 0;
      tracked_prtcl.d_view(index).x = x;
      tracked_prtcl.d_view(index).y = y;
      tracked_prtcl.d_view(index).z = z;
      tracked_prtcl.d_view(index).vx = missing;
      tracked_prtcl.d_view(index).vy = missing;
      tracked_prtcl.d_view(index).vz = missing;
      tracked_prtcl.d_view(index).rho = missing;
      tracked_prtcl.d_view(index).press = missing;
      tracked_prtcl.d_view(index).temp = missing;
      tracked_prtcl.d_view(index).eint = missing;
      tracked_prtcl.d_view(index).scalar0 = missing;

      if (active) {
        int i = static_cast<int>((x - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1) + is;
        int j = js;
        int k = ks;
        if (multi_d) {
          j = static_cast<int>((y - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2) + js;
        }
        if (three_d) {
          k = static_cast<int>((z - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3) + ks;
        }
        i = (i < is) ? is : ((i > ie) ? ie : i);
        j = (j < js) ? js : ((j > je) ? je : j);
        k = (k < ks) ? ks : ((k > ke) ? ke : k);

        tracked_prtcl.d_view(index).vx = is_lagrangian_mc ? w0(m,IVX,k,j,i) : pr(IPVX,p);
        tracked_prtcl.d_view(index).vy = is_lagrangian_mc ? w0(m,IVY,k,j,i) : pr(IPVY,p);
        tracked_prtcl.d_view(index).vz = is_lagrangian_mc ? w0(m,IVZ,k,j,i) : pr(IPVZ,p);
        tracked_prtcl.d_view(index).rho = w0(m,IDN,k,j,i);
        if (eos.is_ideal) {
          tracked_prtcl.d_view(index).eint = w0(m,IEN,k,j,i);
          tracked_prtcl.d_view(index).press = eos.IdealGasPressure(w0(m,IEN,k,j,i));
          tracked_prtcl.d_view(index).temp =
              tracked_prtcl.d_view(index).press / tracked_prtcl.d_view(index).rho;
        } else {
          tracked_prtcl.d_view(index).temp = eos.iso_cs*eos.iso_cs;
          tracked_prtcl.d_view(index).press =
              tracked_prtcl.d_view(index).rho * tracked_prtcl.d_view(index).temp;
          tracked_prtcl.d_view(index).eint = missing;
        }
        if (nscalars > 0) {
          tracked_prtcl.d_view(index).scalar0 = w0(m,nfluid,k,j,i);
        }
      }
    }
  });
  
  Kokkos::fence();
  Kokkos::deep_copy(npout, d_counter);
  
  // share number of tracked particles to be output across all ranks
  npout_eachrank[global_variable::my_rank] = npout;
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(&npout, 1, MPI_INT, npout_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#endif

  tracked_prtcl.resize(npout);
  // sync tracked particle device array with host
  tracked_prtcl.template modify<DevExeSpace>();
  tracked_prtcl.template sync<HostMemSpace>();

  // copy host view into host outpart array
  Kokkos::realloc(outpart, npout);
  Kokkos::deep_copy(outpart, tracked_prtcl.h_view);

}

//----------------------------------------------------------------------------------------
//! \fn void TrackedParticleOutput:::WriteOutputFile(Mesh *pm)
//! \brief Cycles over all tracked particles on this rank and writes ouput data
//! With MPI, all particles are written to the same file.

void TrackedParticleOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  int big_end = IsBigEndian(); // =1 on big endian machine

  // create filename: "trk/file_basename".trk
  std::string fname;
  fname.assign("trk/");
  fname.append(out_params.file_basename);
  fname.append(".trk");

  // Root process opens/creates file and appends string
  if (global_variable::my_rank == 0) {
    std::stringstream msg;
    msg << std::endl << "# AthenaK tracked particle data at time= " << pm->time
        << "  nranks= " << global_variable::nranks
        << "  cycle=" << pm->ncycle
        << "  ntracked_prtcls=" << ntrack
        << "  nvalues=" << NTRACK_VALUES
        << "  variables=" << TRACK_VARIABLES << std::endl;
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      exit(EXIT_FAILURE);
    }
    std::fprintf(pfile,"%s",msg.str().c_str());
    std::fclose(pfile);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  struct stat file_stat;
  if (stat(fname.c_str(), &file_stat) != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "Output file '" << fname << "' could not be stat'ed" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Now all ranks open file and append data
  IOWrapper partfile;
  partfile.Open(fname.c_str(), IOWrapper::FileMode::append);
  std::size_t header_offset = static_cast<std::size_t>(file_stat.st_size);

  // allocate 1D vector of floats used to convert and output particle data
  float *data = new float[NTRACK_VALUES*npout];
  // Loop over particles, load positions into data[]
  for (int p=0; p<npout; ++p) {
    data[(NTRACK_VALUES*p)    ] = static_cast<float>(outpart(p).x);
    data[(NTRACK_VALUES*p) + 1] = static_cast<float>(outpart(p).y);
    data[(NTRACK_VALUES*p) + 2] = static_cast<float>(outpart(p).z);
    data[(NTRACK_VALUES*p) + 3] = static_cast<float>(outpart(p).vx);
    data[(NTRACK_VALUES*p) + 4] = static_cast<float>(outpart(p).vy);
    data[(NTRACK_VALUES*p) + 5] = static_cast<float>(outpart(p).vz);
    data[(NTRACK_VALUES*p) + 6] = static_cast<float>(outpart(p).rho);
    data[(NTRACK_VALUES*p) + 7] = static_cast<float>(outpart(p).press);
    data[(NTRACK_VALUES*p) + 8] = static_cast<float>(outpart(p).temp);
    data[(NTRACK_VALUES*p) + 9] = static_cast<float>(outpart(p).eint);
    data[(NTRACK_VALUES*p) +10] = static_cast<float>(outpart(p).scalar0);
    data[(NTRACK_VALUES*p) +11] = static_cast<float>(outpart(p).gid);
    data[(NTRACK_VALUES*p) +12] = static_cast<float>(outpart(p).level);
    data[(NTRACK_VALUES*p) +13] = static_cast<float>(outpart(p).active);
  }

  // calculate local data offset
  std::vector<int> rank_offset(global_variable::nranks, 0);
  int npout_min = npout_eachrank[0];
  for (int n=1; n<global_variable::nranks; ++n) {
    rank_offset[n] = rank_offset[n-1] + npout_eachrank[n-1];
    npout_min = std::min(npout_min, npout_eachrank[n]);
  }

  if (global_variable::my_rank == 0) {
    std::vector<float> inactive(NTRACK_VALUES*ntrack, 0.0);
    for (int p=0; p<ntrack; ++p) {
      inactive[(NTRACK_VALUES*p) + 13] = 0.0;
    }
    if (partfile.Write_any_type_at(inactive.data(), inactive.size(), header_offset, "float")
        != static_cast<int>(inactive.size())) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "inactive tracked particle records were not initialized"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Write tracked particle data collectively over minimum shared number of prtcls
  for (int p=0; p<npout_min; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset =
        header_offset + NTRACK_VALUES * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at_all(&(data[NTRACK_VALUES*p]), NTRACK_VALUES,
                                       myoffset, "float") != NTRACK_VALUES) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Write particle positions individually for remaining particles on each rank
  for (int p=npout_min; p<npout; ++p) {
    // offset computed assuming tags run 0...(ntrack-1) sequentially
    std::size_t myoffset =
        header_offset + NTRACK_VALUES * outpart(p).tag * sizeof(float);
    // Write particle positions collectively for minimum number of particles across ranks
    if (partfile.Write_any_type_at(&(data[NTRACK_VALUES*p]), NTRACK_VALUES,
                                   myoffset, "float") != NTRACK_VALUES) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "particle data not written correctly to tracked particle file"
          << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // close the output file and clean up
  partfile.Close();
  delete[] data;

  // increment counters
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}
