//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file history.cpp
//  \brief writes history output data, volume-averaged quantities that are output
//         frequently in time to trace their evolution.

#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/viscosity.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/turb_driver.hpp"
#include "z4c/z4c.hpp"
#include "coordinates/adm.hpp"
#include "outputs.hpp"

namespace {
constexpr int nturb_history = 9;

// Positive viscous stress at the lower face normal to dir. Match
// IsotropicViscousFlux: face-averaged rho*nu and centered transverse derivatives.
KOKKOS_INLINE_FUNCTION
Kokkos::Array<Real, 3> ViscousFaceStress(const DvceArray5D<Real> &w, int m,
    int k, int j, int i, int dir, int ndim, const Real *idx, Real nu) {
  int r[3] = {i, j, k};
  int l[3] = {i, j, k};
  --l[dir];
  Kokkos::Array<Real, 3> stress;
  for (int n=0; n<3; ++n) {
    stress[n] = (w(m,IVX+n,r[2],r[1],r[0]) -
                 w(m,IVX+n,l[2],l[1],l[0]))*idx[dir];
  }
  Real div = stress[dir];
  stress[dir] *= 2.0;
  for (int d=0; d<ndim; ++d) {
    if (d == dir) continue;
    ++r[d]; ++l[d];
    Real shear = w(m,IVX+dir,r[2],r[1],r[0]) + w(m,IVX+dir,l[2],l[1],l[0]);
    Real normal = w(m,IVX+d,r[2],r[1],r[0]) + w(m,IVX+d,l[2],l[1],l[0]);
    r[d] -= 2; l[d] -= 2;
    shear -= w(m,IVX+dir,r[2],r[1],r[0]) + w(m,IVX+dir,l[2],l[1],l[0]);
    normal -= w(m,IVX+d,r[2],r[1],r[0]) + w(m,IVX+d,l[2],l[1],l[0]);
    ++r[d]; ++l[d];
    stress[d] += 0.25*idx[d]*shear;
    div += 0.25*idx[d]*normal;
  }
  stress[dir] -= (2.0/3.0)*div;
  Real mu = 0.5*nu*(w(m,IDN,r[2],r[1],r[0]) + w(m,IDN,l[2],l[1],l[0]));
  for (int n=0; n<3; ++n) stress[n] *= mu;
  return stress;
}
}  // namespace

//----------------------------------------------------------------------------------------
// Constructor: also calls BaseTypeOutput base class constructor

HistoryOutput::HistoryOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op),
  turbulence_history_(pin->GetOrAddBoolean(op.block_name, "turbulence", false)) {
  if (turbulence_history_) {
    auto *pc = pm->pmb_pack->pcoord;
    bool periodic = true;
    int ndim = pm->three_d ? 3 : (pm->two_d ? 2 : 1);
    for (int f=0; f<2*ndim; ++f) {
      periodic = periodic && (pm->mesh_bcs[f] == BoundaryFlag::periodic);
    }
    if (!periodic || pm->multilevel || pm->mb_indcs.ng < 2 ||
        pc->is_special_relativistic || pc->is_general_relativistic ||
        pc->is_dynamical_relativistic || op.user_hist_only ||
        (pm->pmb_pack->phydro == nullptr && pm->pmb_pack->pmhd == nullptr)) {
      std::cout << "### FATAL ERROR: turbulence history requires a periodic, uniform "
                << "Newtonian fluid grid with at least two ghost cells and "
                << "user_hist_only=false" << std::endl;
      exit(EXIT_FAILURE);
    }
    auto *ph = pm->pmb_pack->phydro;
    auto *pb = pm->pmb_pack->pmhd;
    if ((ph && (7 + ph->peos->eos_data.is_ideal + ph->nscalars + nturb_history >
                NHISTORY_VARIABLES)) ||
        (pb && (10 + pb->peos->eos_data.is_ideal + pb->nscalars + nturb_history >
                NHISTORY_VARIABLES))) {
      std::cout << "### FATAL ERROR: turbulence history plus fluid/scalar columns "
                << "exceeds NHISTORY_VARIABLES=" << NHISTORY_VARIABLES << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // cycle through physics modules and add HistoryData struct for each
  hist_data.clear();

  if (pm->pgen->user_hist && op.user_hist_only) {
    hist_data.emplace_back(PhysicsModule::UserDefined);
  } else {
    if (pm->pmb_pack->phydro != nullptr) {
      hist_data.emplace_back(PhysicsModule::HydroDynamics);
    }
    if (pm->pmb_pack->pmhd != nullptr) {
      hist_data.emplace_back(PhysicsModule::MagnetoHydroDynamics);
    }
    if (pm->pgen->user_hist) {
      hist_data.emplace_back(PhysicsModule::UserDefined);
    }
  }

  if (pm->pmb_pack->pz4c != nullptr) {
    hist_data.emplace_back(PhysicsModule::SpaceTimeDynamics);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void HistoryOutput::LoadOutputData()
//  \brief Wrapper function that cycles through hist_data vector and calls
//  appropriate LoadXXXData() function for that physics

void HistoryOutput::LoadOutputData(Mesh *pm) {
  for (auto &data : hist_data) {
    if (data.physics == PhysicsModule::HydroDynamics) {
      LoadHydroHistoryData(&data, pm);
    } else if (data.physics == PhysicsModule::MagnetoHydroDynamics) {
      LoadMHDHistoryData(&data, pm);
    } else if (data.physics == PhysicsModule::SpaceTimeDynamics) {
      LoadZ4cHistoryData(&data, pm);
    } else if (data.physics == PhysicsModule::UserDefined) {
      (pm->pgen->user_hist_func)(&data, pm);
    }
    if (turbulence_history_ && (data.physics == PhysicsModule::HydroDynamics ||
                               data.physics == PhysicsModule::MagnetoHydroDynamics)) {
      LoadTurbulenceHistoryData(&data, pm);
    }
  }
}

//----------------------------------------------------------------------------------------
// Instantaneous spatial-operator rates, evaluated only at history output. Periodic
// summation by parts gives dZ/dt|a = integral a.curl(omega); this avoids acceleration
// halos (forcing ghosts are not populated) and additional mesh-sized work arrays.
void HistoryOutput::LoadTurbulenceHistoryData(HistoryData *pdata, Mesh *pm) {
  const char *labels[nturb_history] = {"enstrophy", "drag-KE", "drag-enst", "visc-KE",
                                    "visc-enst", "force-KE", "force-enst", "p-dilat",
                                    "enst-comp"};
  auto *pack = pm->pmb_pack;
  bool hydro = pdata->physics == PhysicsModule::HydroDynamics;
  auto w = hydro ? pack->phydro->w0 : pack->pmhd->w0;
  auto eos = hydro ? pack->phydro->peos->eos_data : pack->pmhd->peos->eos_data;
  auto *visc = hydro ? pack->phydro->pvisc : pack->pmhd->pvisc;
  auto *src = hydro ? pack->phydro->psrc : pack->pmhd->psrc;
  Real nu = visc ? visc->nu_iso : 0.0;
  Real alpha = (src && src->linear_drag) ? src->drag_rate : 0.0;
  bool forced = pack->pturb != nullptr;
  DvceArray5D<Real> force;
  if (forced) force = pack->pturb->force;
  auto size = pack->pmb->mb_size;
  auto &ind = pm->mb_indcs;
  int is = ind.is, js = ind.js, ks = ind.ks;
  int nx1 = ind.nx1, nx2 = ind.nx2, nx3 = ind.nx3;
  int ndim = pm->three_d ? 3 : (pm->two_d ? 2 : 1);
  int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  array_sum::GlobalSum sums;
  Kokkos::parallel_reduce("TurbulenceHistory",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, pack->nmb_thispack*nkji),
      KOKKOS_LAMBDA(const int cell, array_sum::GlobalSum &sum) {
    int m = cell/nkji;
    int k = (cell-m*nkji)/nji + ks;
    int j = (cell-m*nkji)%nji/nx1 + js;
    int i = cell%nx1 + is;
    Real idx[3] = {1.0/size.d_view(m).dx1, 1.0/size.d_view(m).dx2,
                   1.0/size.d_view(m).dx3};
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real rho = w(m,IDN,k,j,i);
    Real grad[3][3] = {};  // grad[velocity component][derivative direction]
    Real visc_rhs[3] = {};
    for (int d=0; d<ndim; ++d) {
      int p[3] = {i,j,k}, q[3] = {i,j,k};
      ++p[d]; --q[d];
      for (int n=0; n<3; ++n) {
        grad[n][d] = 0.5*idx[d]*(w(m,IVX+n,p[2],p[1],p[0]) -
                                 w(m,IVX+n,q[2],q[1],q[0]));
      }
      if (nu != 0.0) {
        auto lo = ViscousFaceStress(w,m,k,j,i,d,ndim,idx,nu);
        auto hi = ViscousFaceStress(w,m,p[2],p[1],p[0],d,ndim,idx,nu);
        for (int n=0; n<3; ++n) visc_rhs[n] += (hi[n]-lo[n])*idx[d];
      }
    }
    Real omega[3] = {grad[2][1]-grad[1][2], grad[0][2]-grad[2][0],
                     grad[1][0]-grad[0][1]};
    Real omega2 = SQR(omega[0]) + SQR(omega[1]) + SQR(omega[2]);
    Real div = grad[0][0] + grad[1][1] + grad[2][2];
    Real pressure = eos.is_ideal ? (eos.gamma-1.0)*w(m,IEN,k,j,i) :
                                  SQR(eos.iso_cs)*rho;
    sum.the_array[0] += 0.5*vol*omega2;
    sum.the_array[2] += alpha*vol*omega2;
    sum.the_array[7] += vol*pressure*div;
    sum.the_array[8] -= 0.5*vol*omega2*div;
    for (int n=0; n<3; ++n) {
      Real v = w(m,IVX+n,k,j,i);
      Real a = forced ? force(m,n,k,j,i) : 0.0;
      sum.the_array[1] += alpha*vol*rho*v*v;
      sum.the_array[3] -= vol*v*visc_rhs[n];
      sum.the_array[5] += vol*rho*v*a;
      if (nu == 0.0 && !forced) continue;
      // curl(curl(u)) with the SAME centered derivative as omega above.
      Real curl_omega = 0.0;
      for (int d=0; d<ndim; ++d) {
        if (d == n) continue;
        int p[3] = {i,j,k}, q[3] = {i,j,k};
        p[d] += 2; q[d] -= 2;
        curl_omega -= 0.25*SQR(idx[d])*(w(m,IVX+n,p[2],p[1],p[0]) -
                                       2.0*v + w(m,IVX+n,q[2],q[1],q[0]));
        if (n >= ndim) continue;
        --p[d]; ++q[d]; ++p[n]; ++q[n];
        Real cross = w(m,IVX+d,p[2],p[1],p[0]) - w(m,IVX+d,q[2],q[1],q[0]);
        p[n] -= 2; q[n] -= 2;
        cross -= w(m,IVX+d,p[2],p[1],p[0]) - w(m,IVX+d,q[2],q[1],q[0]);
        curl_omega += 0.25*idx[d]*idx[n]*cross;
      }
      sum.the_array[4] += vol*(visc_rhs[n]/rho)*curl_omega;
      sum.the_array[6] += vol*a*curl_omega;
    }
  }, Kokkos::Sum<array_sum::GlobalSum>(sums));
  for (int n=0; n<nturb_history; ++n) {
    pdata->label[pdata->nhist+n] = labels[n];
    pdata->hdata[pdata->nhist+n] = sums.the_array[n];
  }
  pdata->nhist += nturb_history;
}

//----------------------------------------------------------------------------------------
//! \fn void HistoryOutput::LoadHydroHistoryData()
//  \brief Compute and store history data over all MeshBlocks on this rank
//  Data is stored in a Real array defined in derived class.

void HistoryOutput::LoadHydroHistoryData(HistoryData *pdata, Mesh *pm) {
  auto &eos_data = pm->pmb_pack->phydro->peos->eos_data;
  int &nhydro_ = pm->pmb_pack->phydro->nhydro;
  int &nscalars_ = pm->pmb_pack->phydro->nscalars;

  // set number of and names of history variables for hydro
  if (eos_data.is_ideal) {
    pdata->nhist = 8;
  } else {
    pdata->nhist = 7;
  }
  if (nscalars_>0) {
    pdata->nhist += nscalars_;
  }
  pdata->label[IDN] = "mass";
  pdata->label[IM1] = "1-mom";
  pdata->label[IM2] = "2-mom";
  pdata->label[IM3] = "3-mom";
  if (eos_data.is_ideal) {
    pdata->label[IEN] = "tot-E";
  }
  pdata->label[nhydro_  ] = "1-KE";
  pdata->label[nhydro_+1] = "2-KE";
  pdata->label[nhydro_+2] = "3-KE";
  for (int s=0; s<nscalars_; ++s) {
    std::ostringstream labelSS;
    labelSS << "scal-" << s;
    pdata->label[nhydro_+3+s] = labelSS.str();
  }

  // capture class variables for kernel
  auto &u0_ = pm->pmb_pack->phydro->u0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    // Hydro conserved variables:
    array_sum::GlobalSum hvars;
    hvars.the_array[IDN] = vol*u0_(m,IDN,k,j,i);
    hvars.the_array[IM1] = vol*u0_(m,IM1,k,j,i);
    hvars.the_array[IM2] = vol*u0_(m,IM2,k,j,i);
    hvars.the_array[IM3] = vol*u0_(m,IM3,k,j,i);
    if (eos_data.is_ideal) {
      hvars.the_array[IEN] = vol*u0_(m,IEN,k,j,i);
    }

    // Hydro KE
    hvars.the_array[nhydro_  ] = vol*0.5*SQR(u0_(m,IM1,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nhydro_+1] = vol*0.5*SQR(u0_(m,IM2,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nhydro_+2] = vol*0.5*SQR(u0_(m,IM3,k,j,i))/u0_(m,IDN,k,j,i);

    // Scalar masses
    for (int s=0; s<nscalars_; ++s) {
      hvars.the_array[nhydro_+3+s] = vol*u0_(m,nhydro_+s,k,j,i);
    }

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HistoryOutput::LoadZ4cHistoryData()
//  \brief Compute and store history data over all MeshBlocks on this rank
//  Data is stored in a Real array defined in derived class.

void HistoryOutput::LoadZ4cHistoryData(HistoryData *pdata, Mesh *pm) {
  // set number of and names of history variables for z4c
  pdata->nhist = 9;
  pdata->label[0] = "C-norm2";
  pdata->label[1] = "H-norm2";
  pdata->label[2] = "M-norm2";
  pdata->label[3] = "Z-norm2";
  pdata->label[4] = "Mx-norm2";
  pdata->label[5] = "My-norm2";
  pdata->label[6] = "Mz-norm2";
  pdata->label[7] = "Theta-norm2";
  pdata->label[8] = "Volume";

  // capture class variabels for kernel
  auto &u0_ = pm->pmb_pack->pz4c->u0;
  auto &u_con_ = pm->pmb_pack->pz4c->u_con;
  const int &I_Z4c_Theta_ =  pm->pmb_pack->pz4c->I_Z4C_THETA;
  auto &z4c = pm->pmb_pack->pz4c->z4c;
  auto &adm = pm->pmb_pack->padm->adm;

  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;
  auto &opt = pm->pmb_pack->pz4c->opt;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real detg = adm::SpatialDet(adm.g_dd(m,0,0,k,j,i), adm.g_dd(m,0,1,k,j,i),
                                adm.g_dd(m,0,2,k,j,i), adm.g_dd(m,1,1,k,j,i),
                                adm.g_dd(m,1,2,k,j,i), adm.g_dd(m,2,2,k,j,i));

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3
               * std::sqrt(std::abs(detg));

    // Excise the punctures based on chi
    array_sum::GlobalSum hvars;
    if (z4c.chi(m,k,j,i)>=opt.excise_chi) {
      hvars.the_array[0] = vol*u_con_(m,0,k,j,i); // ||C||^2 (comes already squared)
      hvars.the_array[1] = vol*SQR(u_con_(m,1,k,j,i)); //||H||^2
      hvars.the_array[2] = vol*u_con_(m,2,k,j,i); // ||M||^2 (comes already squared)
      hvars.the_array[3] = vol*u_con_(m,3,k,j,i); // ||Z||^2 (comes already squared)
      hvars.the_array[4] = vol*SQR(u_con_(m,4,k,j,i));      // ||Mx||^2
      hvars.the_array[5] = vol*SQR(u_con_(m,5,k,j,i));      // ||My||^2
      hvars.the_array[6] = vol*SQR(u_con_(m,6,k,j,i));      // ||Mz||^2
      hvars.the_array[7] = vol*SQR(u0_(m,I_Z4c_Theta_,k,j,i)); // ||Theta||^2
      hvars.the_array[8] = vol;
    } else {
      hvars.the_array[0] = 0;
      hvars.the_array[1] = 0;
      hvars.the_array[2] = 0;
      hvars.the_array[3] = 0;
      hvars.the_array[4] = 0;
      hvars.the_array[5] = 0;
      hvars.the_array[6] = 0;
      hvars.the_array[7] = 0;
      hvars.the_array[8] = 0;
    }

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HistoryOutput::LoadMHDHistoryData()
//  \brief Compute and store history data over all MeshBlocks on this rank
//  Data is stored in a Real array defined in derived class.

void HistoryOutput::LoadMHDHistoryData(HistoryData *pdata, Mesh *pm) {
  auto &eos_data = pm->pmb_pack->pmhd->peos->eos_data;
  int &nmhd_ = pm->pmb_pack->pmhd->nmhd;
  int &nscalars_ = pm->pmb_pack->pmhd->nscalars;

  // set number of and names of history variables for mhd
  if (eos_data.is_ideal) {
    pdata->nhist = 11;
  } else {
    pdata->nhist = 10;
  }
  if (nscalars_>0) {
    pdata->nhist += nscalars_;
  }
  pdata->label[IDN] = "mass";
  pdata->label[IM1] = "1-mom";
  pdata->label[IM2] = "2-mom";
  pdata->label[IM3] = "3-mom";
  if (eos_data.is_ideal) {
    pdata->label[IEN] = "tot-E";
  }
  pdata->label[nmhd_  ] = "1-KE";
  pdata->label[nmhd_+1] = "2-KE";
  pdata->label[nmhd_+2] = "3-KE";
  pdata->label[nmhd_+3] = "1-ME";
  pdata->label[nmhd_+4] = "2-ME";
  pdata->label[nmhd_+5] = "3-ME";

  for (int s=0; s<nscalars_; ++s) {
    std::ostringstream labelSS;
    labelSS << "scal-" << s;
    pdata->label[nmhd_+6+s] = labelSS.str();
  }

  // capture class variabels for kernel
  auto &u0_ = pm->pmb_pack->pmhd->u0;
  auto &bx1f = pm->pmb_pack->pmhd->b0.x1f;
  auto &bx2f = pm->pmb_pack->pmhd->b0.x2f;
  auto &bx3f = pm->pmb_pack->pmhd->b0.x3f;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    // MHD conserved variables:
    array_sum::GlobalSum hvars;
    hvars.the_array[IDN] = vol*u0_(m,IDN,k,j,i);
    hvars.the_array[IM1] = vol*u0_(m,IM1,k,j,i);
    hvars.the_array[IM2] = vol*u0_(m,IM2,k,j,i);
    hvars.the_array[IM3] = vol*u0_(m,IM3,k,j,i);
    if (eos_data.is_ideal) {
      hvars.the_array[IEN] = vol*u0_(m,IEN,k,j,i);
    }

    // MHD KE
    hvars.the_array[nmhd_  ] = vol*0.5*SQR(u0_(m,IM1,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nmhd_+1] = vol*0.5*SQR(u0_(m,IM2,k,j,i))/u0_(m,IDN,k,j,i);
    hvars.the_array[nmhd_+2] = vol*0.5*SQR(u0_(m,IM3,k,j,i))/u0_(m,IDN,k,j,i);

    // MHD ME
    hvars.the_array[nmhd_+3] = vol*0.25*(SQR(bx1f(m,k,j,i+1)) + SQR(bx1f(m,k,j,i)));
    hvars.the_array[nmhd_+4] = vol*0.25*(SQR(bx2f(m,k,j+1,i)) + SQR(bx2f(m,k,j,i)));
    hvars.the_array[nmhd_+5] = vol*0.25*(SQR(bx3f(m,k+1,j,i)) + SQR(bx3f(m,k,j,i)));

    // Scalar masses
    for (int s=0; s<nscalars_; ++s) {
      hvars.the_array[nmhd_+6+s] = vol*u0_(m,nmhd_+s,k,j,i);
    }

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HistoryOutput::WriteOutputFile()
//  \brief Cycles through hist_data vector and writes history file for each component

void HistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  for (auto &data : hist_data) {
    // first, perform in-place sum over all MPI ranks
#if MPI_PARALLEL_ENABLED
    if (global_variable::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, &(data.hdata[0]), data.nhist, MPI_ATHENA_REAL,
         MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(&(data.hdata[0]), &(data.hdata[0]), data.nhist,
         MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    // only the master rank writes the file
    if (global_variable::my_rank == 0) {
      // create filename: "file_basename" + ".physics" + ".hst"
      // There is no file number or id in history output filenames.
      std::string fname;
      fname.assign(out_params.file_basename);
      switch (data.physics) {
        case PhysicsModule::HydroDynamics:
          fname.append(".hydro");
          break;
        case PhysicsModule::MagnetoHydroDynamics:
          fname.append(".mhd");
          break;
        case PhysicsModule::SpaceTimeDynamics:
          fname.append(".z4c");
        case PhysicsModule::UserDefined:
          fname.append(".user");
          break;
        default:
          break;
      }
      fname.append(".hst");

      // open file for output
      FILE *pfile;
      if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
        exit(EXIT_FAILURE);
      }

      // Write header, if it has not been written already
      if (!(data.header_written)) {
        int iout = 1;
        std::fprintf(pfile,"# Athena++ history data\n");
        std::fprintf(pfile,"#  [%d]=time      ", iout++);
        std::fprintf(pfile,"[%d]=dt       ", iout++);
        for (int n=0; n<data.nhist; ++n) {
          std::fprintf(pfile,"[%d]=%.10s    ", iout++, data.label[n].c_str());
        }
        std::fprintf(pfile,"\n");                              // terminate line
        data.header_written = true;
      }

      // write history variables
      std::fprintf(pfile, out_params.data_format.c_str(), pm->time);
      std::fprintf(pfile, out_params.data_format.c_str(), pm->dt);
      for (int n=0; n<data.nhist; ++n)
        std::fprintf(pfile, out_params.data_format.c_str(), data.hdata[n]);
      std::fprintf(pfile,"\n"); // terminate line
      std::fclose(pfile);
    }
  } // End loop over hist_data vector

  // increment counters, clean up
  if (out_params.last_time < 0.0) {
    out_params.last_time = pm->time;
  } else {
    out_params.last_time += out_params.dt;
  }
  pin->SetReal(out_params.block_name, "last_time", out_params.last_time);
  return;
}
