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
#include <cmath>
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
#include "units/units.hpp"
#include "outputs.hpp"

namespace {
constexpr int NTRACK_VALUES = 32;
constexpr int TRACK_TAG_INDEX = 0;
constexpr int TRACK_ACTIVE_INDEX = 14;
const char *TRACK_VARIABLES =
    "tag x y z vx vy vz rho press temp eint scalar0 gid level active "
    "edot_cool edot_heat edot_net dedt_rad_mass dTdt_rad tcool "
    "entropy ln_entropy divv dTdt_ad T_mix_scalar T_minus_T_mix_scalar "
    "T_label_mix T_minus_T_label_mix gradT_mag grad_scalar_mag strain_mag";
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
  const bool is_lagrangian_mc =
      (pm->pmb_pack->ppart->particle_type == ParticleType::lagrangian_mc);
  track_by_slot =
      is_lagrangian_mc && pm->pmb_pack->ppart->lmc_scheduled_tracking;
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
  if (ntrack <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Tracked particle output requires nparticles > 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (track_by_slot &&
      ntrack != pm->pmb_pack->ppart->lmc_track_nparticles_total) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "scheduled_injection tracking requires <"
              << op.block_name << ">/nparticles to equal "
              << "<particles>/track_nparticles_total" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // TODO(@user) improve guess below?
  ntrack_thisrank = ntrack;

  auto has_problem = [&](const char *name) {
    return pin->DoesParameterExist("problem", name);
  };
  const bool have_trml_cooling_inputs =
      has_problem("rho_0") && has_problem("pgas_0") &&
      (has_problem("density_contrast") || has_problem("contrast"));
  if (have_trml_cooling_inputs) {
    trml_rho_0 = pin->GetReal("problem", "rho_0");
    trml_pgas_0 = pin->GetReal("problem", "pgas_0");
    const Real contrast = has_problem("density_contrast") ?
        pin->GetReal("problem", "density_contrast") :
        pin->GetReal("problem", "contrast");
    trml_T_cold = trml_pgas_0/trml_rho_0/contrast;
    trml_T_hot = trml_pgas_0/trml_rho_0;
    trml_t_cool_start = has_problem("t_cool_start") ?
        pin->GetReal("problem", "t_cool_start") : 0.0;

    const Real velocity = has_problem("velocity") ? pin->GetReal("problem", "velocity") : 0.0;
    const Real velocity_abs = std::fabs(velocity);
    const bool is_mhd = (pm->pmb_pack->pmhd != nullptr);
    EOS_Data &eos = (is_mhd) ?
        pm->pmb_pack->pmhd->peos->eos_data : pm->pmb_pack->phydro->peos->eos_data;
    trml_cooling::ReferenceState ref;
    ref.rho_0 = trml_rho_0;
    ref.pgas_0 = trml_pgas_0;
    ref.contrast = contrast;
    ref.gamma = eos.gamma;
    ref.t_shear = (velocity_abs > 0.0) ? 1.0/velocity_abs :
                  std::numeric_limits<Real>::infinity();
    ref.allow_zero_shear = true;
    trml_cooling::Setup setup =
        trml_cooling::ReadInputs(pin, pm->pmb_pack->punit, ref, false, false);
    if (setup.valid) {
      trml_cooling_params = setup.params;
      trml_t_cool_0 = setup.t_cool_0;
      trml_cooling_diagnostics = true;
    }
  }
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
  bool track_by_slot_ = track_by_slot;
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
  const bool trml_cooling_diagnostics_ = trml_cooling_diagnostics;
  const trml_cooling::Params trml_cooling_params_ = trml_cooling_params;
  const Real trml_T_hot_ = trml_T_hot;
  const Real trml_T_cold_ = trml_T_cold;
  const Real trml_t_cool_start_ = trml_t_cool_start;
  const Real time_ = pm->time;

  // Create device-side counter
  Kokkos::View<int> d_counter("counter");
  Kokkos::deep_copy(d_counter, 0);

  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    const int tag = pi(PTAG,p);
    int track_index = -1;
    if (track_by_slot_) {
      track_index = pi(PTRACK,p);
    } else if (tag >= 0 && tag < ntrack_) {
      track_index = tag;
    }
    if (track_index >= 0 && track_index < ntrack_) {
      int index = Kokkos::atomic_fetch_add(&d_counter(),1);
      int m = pi(PGID,p) - gids;
      bool active = (m >= 0 && m < nmb);
      if (is_lagrangian_mc && pi(PLASTMOVE,p) < 0) {
        active = false;
      }
      Real x = is_lagrangian_mc ? pr(LMCX,p) : pr(IPX,p);
      Real y = is_lagrangian_mc ? pr(LMCY,p) : pr(IPY,p);
      Real z = is_lagrangian_mc ? pr(LMCZ,p) : pr(IPZ,p);
      tracked_prtcl.d_view(index).track_index = track_index;
      tracked_prtcl.d_view(index).tag = tag;
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
      tracked_prtcl.d_view(index).edot_cool = missing;
      tracked_prtcl.d_view(index).edot_heat = missing;
      tracked_prtcl.d_view(index).edot_net = missing;
      tracked_prtcl.d_view(index).dedt_rad_mass = missing;
      tracked_prtcl.d_view(index).dTdt_rad = missing;
      tracked_prtcl.d_view(index).tcool = missing;
      tracked_prtcl.d_view(index).entropy = missing;
      tracked_prtcl.d_view(index).ln_entropy = missing;
      tracked_prtcl.d_view(index).divv = missing;
      tracked_prtcl.d_view(index).dTdt_ad = missing;
      tracked_prtcl.d_view(index).T_mix_scalar = missing;
      tracked_prtcl.d_view(index).T_minus_T_mix_scalar = missing;
      tracked_prtcl.d_view(index).T_label_mix = missing;
      tracked_prtcl.d_view(index).T_minus_T_label_mix = missing;
      tracked_prtcl.d_view(index).gradT_mag = missing;
      tracked_prtcl.d_view(index).grad_scalar_mag = missing;
      tracked_prtcl.d_view(index).strain_mag = missing;

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

        const Real vx = is_lagrangian_mc ? w0(m,IVX,k,j,i) : pr(IPVX,p);
        const Real vy = is_lagrangian_mc ? w0(m,IVY,k,j,i) : pr(IPVY,p);
        const Real vz = is_lagrangian_mc ? w0(m,IVZ,k,j,i) : pr(IPVZ,p);
        const Real rho = w0(m,IDN,k,j,i);
        Real press = missing;
        Real temp = missing;
        Real eint = missing;
        tracked_prtcl.d_view(index).vx = vx;
        tracked_prtcl.d_view(index).vy = vy;
        tracked_prtcl.d_view(index).vz = vz;
        tracked_prtcl.d_view(index).rho = rho;
        if (eos.is_ideal) {
          eint = w0(m,IEN,k,j,i);
          press = eos.IdealGasPressure(eint);
          temp = press/rho;
        } else {
          temp = eos.iso_cs*eos.iso_cs;
          press = rho*temp;
          eint = missing;
        }
        tracked_prtcl.d_view(index).press = press;
        tracked_prtcl.d_view(index).temp = temp;
        tracked_prtcl.d_view(index).eint = eint;

        const int im = (i > is) ? i - 1 : i;
        const int ip = (i < ie) ? i + 1 : i;
        const int jm = (multi_d && j > js) ? j - 1 : j;
        const int jp = (multi_d && j < je) ? j + 1 : j;
        const int km = (three_d && k > ks) ? k - 1 : k;
        const int kp = (three_d && k < ke) ? k + 1 : k;
        const Real inv_dx = (ip > im) ? 1.0/((ip - im)*mbsize.d_view(m).dx1) : 0.0;
        const Real inv_dy = (jp > jm) ? 1.0/((jp - jm)*mbsize.d_view(m).dx2) : 0.0;
        const Real inv_dz = (kp > km) ? 1.0/((kp - km)*mbsize.d_view(m).dx3) : 0.0;

        const Real dvx_dx = (w0(m,IVX,k,j,ip) - w0(m,IVX,k,j,im))*inv_dx;
        const Real dvy_dx = (w0(m,IVY,k,j,ip) - w0(m,IVY,k,j,im))*inv_dx;
        const Real dvz_dx = (w0(m,IVZ,k,j,ip) - w0(m,IVZ,k,j,im))*inv_dx;
        const Real dvx_dy = (w0(m,IVX,k,jp,i) - w0(m,IVX,k,jm,i))*inv_dy;
        const Real dvy_dy = (w0(m,IVY,k,jp,i) - w0(m,IVY,k,jm,i))*inv_dy;
        const Real dvz_dy = (w0(m,IVZ,k,jp,i) - w0(m,IVZ,k,jm,i))*inv_dy;
        const Real dvx_dz = (w0(m,IVX,kp,j,i) - w0(m,IVX,km,j,i))*inv_dz;
        const Real dvy_dz = (w0(m,IVY,kp,j,i) - w0(m,IVY,km,j,i))*inv_dz;
        const Real dvz_dz = (w0(m,IVZ,kp,j,i) - w0(m,IVZ,km,j,i))*inv_dz;
        const Real divv = dvx_dx + dvy_dy + dvz_dz;
        const Real sxx = dvx_dx;
        const Real syy = dvy_dy;
        const Real szz = dvz_dz;
        const Real sxy = 0.5*(dvx_dy + dvy_dx);
        const Real sxz = 0.5*(dvx_dz + dvz_dx);
        const Real syz = 0.5*(dvy_dz + dvz_dy);
        tracked_prtcl.d_view(index).divv = divv;
        tracked_prtcl.d_view(index).strain_mag =
            sqrt(2.0*(SQR(sxx) + SQR(syy) + SQR(szz) +
                      2.0*(SQR(sxy) + SQR(sxz) + SQR(syz))));

        if (eos.is_ideal && rho > 0.0 && press > 0.0 && temp > 0.0) {
          const Real gamma = eos.gamma;
          const Real gm1 = gamma - 1.0;
          tracked_prtcl.d_view(index).entropy = press/pow(rho, gamma);
          tracked_prtcl.d_view(index).ln_entropy = log(press) - gamma*log(rho);
          tracked_prtcl.d_view(index).dTdt_ad = -gm1*temp*divv;

          Real temp_im = eos.IdealGasPressure(w0(m,IEN,k,j,im))/w0(m,IDN,k,j,im);
          Real temp_ip = eos.IdealGasPressure(w0(m,IEN,k,j,ip))/w0(m,IDN,k,j,ip);
          Real temp_jm = eos.IdealGasPressure(w0(m,IEN,k,jm,i))/w0(m,IDN,k,jm,i);
          Real temp_jp = eos.IdealGasPressure(w0(m,IEN,k,jp,i))/w0(m,IDN,k,jp,i);
          Real temp_km = eos.IdealGasPressure(w0(m,IEN,km,j,i))/w0(m,IDN,km,j,i);
          Real temp_kp = eos.IdealGasPressure(w0(m,IEN,kp,j,i))/w0(m,IDN,kp,j,i);
          const Real dT_dx = (temp_ip - temp_im)*inv_dx;
          const Real dT_dy = (temp_jp - temp_jm)*inv_dy;
          const Real dT_dz = (temp_kp - temp_km)*inv_dz;
          tracked_prtcl.d_view(index).gradT_mag =
              sqrt(SQR(dT_dx) + SQR(dT_dy) + SQR(dT_dz));

          if (trml_cooling_diagnostics_) {
            trml_cooling::Rates rates;
            if (time_ > trml_t_cool_start_) {
              rates = trml_cooling::Evaluate(trml_cooling_params_, rho, temp, eint);
            }
            tracked_prtcl.d_view(index).edot_cool = rates.edot_cool;
            tracked_prtcl.d_view(index).edot_heat = rates.edot_heat;
            tracked_prtcl.d_view(index).edot_net = rates.edot_net;
            tracked_prtcl.d_view(index).dedt_rad_mass = -rates.edot_net/rho;
            tracked_prtcl.d_view(index).dTdt_rad = -gm1*rates.edot_net/rho;
            if (rates.edot_net > 0.0) {
              tracked_prtcl.d_view(index).tcool = eint/rates.edot_net;
            }
          }
        }

        if (nscalars > 0) {
          const Real scalar0 = w0(m,nfluid,k,j,i);
          tracked_prtcl.d_view(index).scalar0 = scalar0;
          if (trml_cooling_diagnostics_) {
            const Real T_mix_scalar =
                (1.0 - scalar0)*trml_T_hot_ + scalar0*trml_T_cold_;
            tracked_prtcl.d_view(index).T_mix_scalar = T_mix_scalar;
            tracked_prtcl.d_view(index).T_minus_T_mix_scalar = temp - T_mix_scalar;
          }

          const Real dsc_dx = (w0(m,nfluid,k,j,ip) - w0(m,nfluid,k,j,im))*inv_dx;
          const Real dsc_dy = (w0(m,nfluid,k,jp,i) - w0(m,nfluid,k,jm,i))*inv_dy;
          const Real dsc_dz = (w0(m,nfluid,kp,j,i) - w0(m,nfluid,km,j,i))*inv_dz;
          tracked_prtcl.d_view(index).grad_scalar_mag =
              sqrt(SQR(dsc_dx) + SQR(dsc_dy) + SQR(dsc_dz));
        }
        if (trml_cooling_diagnostics_ && nscalars > 1) {
          const Real T_label_mix = w0(m,nfluid+1,k,j,i);
          tracked_prtcl.d_view(index).T_label_mix = T_label_mix;
          tracked_prtcl.d_view(index).T_minus_T_label_mix = temp - T_label_mix;
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
    data[(NTRACK_VALUES*p)    ] = static_cast<float>(outpart(p).tag);
    data[(NTRACK_VALUES*p) + 1] = static_cast<float>(outpart(p).x);
    data[(NTRACK_VALUES*p) + 2] = static_cast<float>(outpart(p).y);
    data[(NTRACK_VALUES*p) + 3] = static_cast<float>(outpart(p).z);
    data[(NTRACK_VALUES*p) + 4] = static_cast<float>(outpart(p).vx);
    data[(NTRACK_VALUES*p) + 5] = static_cast<float>(outpart(p).vy);
    data[(NTRACK_VALUES*p) + 6] = static_cast<float>(outpart(p).vz);
    data[(NTRACK_VALUES*p) + 7] = static_cast<float>(outpart(p).rho);
    data[(NTRACK_VALUES*p) + 8] = static_cast<float>(outpart(p).press);
    data[(NTRACK_VALUES*p) + 9] = static_cast<float>(outpart(p).temp);
    data[(NTRACK_VALUES*p) +10] = static_cast<float>(outpart(p).eint);
    data[(NTRACK_VALUES*p) +11] = static_cast<float>(outpart(p).scalar0);
    data[(NTRACK_VALUES*p) +12] = static_cast<float>(outpart(p).gid);
    data[(NTRACK_VALUES*p) +13] = static_cast<float>(outpart(p).level);
    data[(NTRACK_VALUES*p) +14] = static_cast<float>(outpart(p).active);
    data[(NTRACK_VALUES*p) +15] = static_cast<float>(outpart(p).edot_cool);
    data[(NTRACK_VALUES*p) +16] = static_cast<float>(outpart(p).edot_heat);
    data[(NTRACK_VALUES*p) +17] = static_cast<float>(outpart(p).edot_net);
    data[(NTRACK_VALUES*p) +18] = static_cast<float>(outpart(p).dedt_rad_mass);
    data[(NTRACK_VALUES*p) +19] = static_cast<float>(outpart(p).dTdt_rad);
    data[(NTRACK_VALUES*p) +20] = static_cast<float>(outpart(p).tcool);
    data[(NTRACK_VALUES*p) +21] = static_cast<float>(outpart(p).entropy);
    data[(NTRACK_VALUES*p) +22] = static_cast<float>(outpart(p).ln_entropy);
    data[(NTRACK_VALUES*p) +23] = static_cast<float>(outpart(p).divv);
    data[(NTRACK_VALUES*p) +24] = static_cast<float>(outpart(p).dTdt_ad);
    data[(NTRACK_VALUES*p) +25] = static_cast<float>(outpart(p).T_mix_scalar);
    data[(NTRACK_VALUES*p) +26] = static_cast<float>(outpart(p).T_minus_T_mix_scalar);
    data[(NTRACK_VALUES*p) +27] = static_cast<float>(outpart(p).T_label_mix);
    data[(NTRACK_VALUES*p) +28] = static_cast<float>(outpart(p).T_minus_T_label_mix);
    data[(NTRACK_VALUES*p) +29] = static_cast<float>(outpart(p).gradT_mag);
    data[(NTRACK_VALUES*p) +30] = static_cast<float>(outpart(p).grad_scalar_mag);
    data[(NTRACK_VALUES*p) +31] = static_cast<float>(outpart(p).strain_mag);
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
      inactive[(NTRACK_VALUES*p) + TRACK_TAG_INDEX] =
          track_by_slot ? -1.0f : static_cast<float>(p);
      inactive[(NTRACK_VALUES*p) + TRACK_ACTIVE_INDEX] = 0.0;
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
    std::size_t myoffset =
        header_offset + NTRACK_VALUES * outpart(p).track_index * sizeof(float);
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
    std::size_t myoffset =
        header_offset + NTRACK_VALUES * outpart(p).track_index * sizeof(float);
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
