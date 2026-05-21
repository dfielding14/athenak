//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file streaming_instability.cpp
//! \brief dust-gas drag relaxation and streaming-instability test problem

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

struct StreamingProblemData {
  Real rho0 = 1.0;
  Real epsilon = 1.0;
  Real amp = 1.0e-4;
  Real kx = 0.0;
  Real kz = 0.0;
  Real gas_vx = 0.0;
  Real gas_vy = 0.0;
  Real gas_vz = 0.0;
  Real part_vx = 0.0;
  Real part_vy = 0.0;
  Real part_vz = 0.0;
  bool use_eigenmode = false;
  Real rhog_re = 1.0;
  Real rhog_im = 0.0;
  Real vgx_re = 1.0;
  Real vgx_im = 0.0;
  Real vgy_re = 0.0;
  Real vgy_im = 0.0;
  Real vgz_re = 1.0;
  Real vgz_im = 0.0;
  Real rhod_re = 1.0;
  Real rhod_im = 0.0;
  Real vdx_re = 1.0;
  Real vdx_im = 0.0;
  Real vdy_re = 0.0;
  Real vdy_im = 0.0;
  Real vdz_re = 1.0;
  Real vdz_im = 0.0;
  bool orbital_terms = false;
  Real omega0 = 1.0;
  Real qshear = 1.5;
  Real pressure_accel = 0.0;
  Real user_drag_stopping_time = 0.1;
};

StreamingProblemData si;

} // namespace

void StreamingInstabilitySource(Mesh *pm, const Real bdt);
void StreamingInstabilityHistory(HistoryData *pdata, Mesh *pm);
void StreamingInstabilityUserDrag(MeshBlockPack *pmbp, const Real bdt);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::StreamingInstability

void ProblemGenerator::StreamingInstability(ParameterInput *pin, const bool restart) {
  user_hist_func = StreamingInstabilityHistory;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  bool missing_fluid = (pmbp->phydro == nullptr) && (pmbp->pmhd == nullptr);
  if (missing_fluid || (pmbp->ppart == nullptr)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "streaming_instability requires <hydro> or <mhd>, "
              << "plus <particles>" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bool use_mhd = (pmbp->pmhd != nullptr);

  std::string mode = pin->GetOrAddString("problem","mode","streaming");
  si.rho0 = pin->GetOrAddReal("problem","rho0",1.0);
  si.epsilon = pin->GetOrAddReal("problem","epsilon",1.0);
  si.amp = pin->GetOrAddReal("problem","amp",1.0e-4);
  if (pin->DoesBlockExist("drag_particles")) {
    std::string drag_model = pin->GetOrAddString("drag_particles","model","none");
    if (drag_model.compare("user") == 0) {
      Real default_tstop = pmbp->ppart->drag_stopping_time;
      si.user_drag_stopping_time = pin->GetOrAddReal("problem","user_drag_stopping_time",
                                                     default_tstop);
      if (si.user_drag_stopping_time <= 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<problem>/user_drag_stopping_time must be > 0"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      user_particle_drag_func = StreamingInstabilityUserDrag;
    }
  }

  Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real lz = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  int nwx = pin->GetOrAddInteger("problem","nwx",1);
  int nwz = pin->GetOrAddInteger("problem","nwz",1);
  si.kx = (lx > 0.0) ? 2.0*M_PI*static_cast<Real>(nwx)/lx : 0.0;
  si.kz = (lz > 0.0) ? 2.0*M_PI*static_cast<Real>(nwz)/lz : 0.0;

  if (mode.compare("relaxation") == 0) {
    si.orbital_terms = false;
    si.use_eigenmode = false;
    si.gas_vx = pin->GetOrAddReal("problem","gas_vx",0.0);
    si.gas_vy = pin->GetOrAddReal("problem","gas_vy",0.0);
    si.gas_vz = pin->GetOrAddReal("problem","gas_vz",0.0);
    si.part_vx = pin->GetOrAddReal("problem","particle_vx",0.1);
    si.part_vy = pin->GetOrAddReal("problem","particle_vy",0.0);
    si.part_vz = pin->GetOrAddReal("problem","particle_vz",0.0);
    si.amp = 0.0;
  } else {
    if (use_mhd) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "streaming_instability mode=streaming currently "
                << "requires <hydro>; use mode=relaxation for MHD drag smoke tests"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    user_srcs = true;
    Real tau = pin->GetOrAddReal("problem","tau_s",0.1);
    Real eta_vk = pin->GetOrAddReal("problem","eta_vk",0.05);
    si.omega0 = pin->GetOrAddReal("problem","omega0",1.0);
    si.qshear = pin->GetOrAddReal("problem","qshear",1.5);
    si.pressure_accel = pin->GetOrAddReal("problem","pressure_accel",
                                          2.0*eta_vk*si.omega0);
    si.orbital_terms = true;
    user_srcs_func = StreamingInstabilitySource;
    Real denom = (1.0 + si.epsilon)*(1.0 + si.epsilon) + tau*tau;
    si.gas_vx = 2.0*si.epsilon*tau*eta_vk/denom;
    si.gas_vy = -(1.0 + si.epsilon + tau*tau)*eta_vk/denom;
    si.gas_vz = 0.0;
    si.part_vx = -2.0*tau*eta_vk/denom;
    si.part_vy = -(1.0 + si.epsilon)*eta_vk/denom;
    si.part_vz = 0.0;

    std::string perturbation = pin->GetOrAddString("problem","perturbation","eigenmode");
    si.use_eigenmode = (perturbation.compare("eigenmode") == 0);
    if (si.use_eigenmode) {
      si.rhog_re = pin->GetOrAddReal("problem","eigen_rhog_re",
                                     -6.3941133490089087e-04);
      si.rhog_im = pin->GetOrAddReal("problem","eigen_rhog_im",
                                      4.8248529845820872e-04);
      si.vgx_re = pin->GetOrAddReal("problem","eigen_vgx_re",
                                    -7.0111836003118488e-02);
      si.vgx_im = pin->GetOrAddReal("problem","eigen_vgx_im",
                                    -2.2111103367351170e-02);
      si.vgy_re = pin->GetOrAddReal("problem","eigen_vgy_re",
                                     1.2150806113551665e-01);
      si.vgy_im = pin->GetOrAddReal("problem","eigen_vgy_im",
                                    -3.9156174090960547e-02);
      si.vgz_re = pin->GetOrAddReal("problem","eigen_vgz_re",
                                     3.5079348742086472e-02);
      si.vgz_im = pin->GetOrAddReal("problem","eigen_vgz_im",
                                     1.1032183616959859e-02);
      si.rhod_re = pin->GetOrAddReal("problem","eigen_rhod_re",1.0);
      si.rhod_im = pin->GetOrAddReal("problem","eigen_rhod_im",0.0);
      si.vdx_re = pin->GetOrAddReal("problem","eigen_vdx_re",
                                     5.7833082952576656e-03);
      si.vdx_im = pin->GetOrAddReal("problem","eigen_vdx_im",
                                    -2.3675841652271375e-02);
      si.vdy_re = pin->GetOrAddReal("problem","eigen_vdy_re",
                                     1.1725574099850000e-01);
      si.vdy_im = pin->GetOrAddReal("problem","eigen_vdy_im",
                                    -6.6509847792263046e-03);
      si.vdz_re = pin->GetOrAddReal("problem","eigen_vdz_re",
                                     2.9536786564495800e-02);
      si.vdz_im = pin->GetOrAddReal("problem","eigen_vdz_im",
                                     1.7505935642441284e-02);
    } else if (perturbation.compare("cosine") != 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "streaming_instability perturbation must be "
                << "eigenmode or cosine" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (restart) {return;}

  auto u0 = use_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  EOS_Data eos = use_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  Real rho0 = si.rho0;
  Real amp = si.amp;
  Real kx = si.kx;
  Real kz = si.kz;
  Real gvx = si.gas_vx;
  Real gvy = si.gas_vy;
  Real gvz = si.gas_vz;
  bool use_eigenmode = si.use_eigenmode;
  Real rhog_re = si.rhog_re;
  Real rhog_im = si.rhog_im;
  Real vgx_re = si.vgx_re;
  Real vgx_im = si.vgx_im;
  Real vgy_re = si.vgy_re;
  Real vgy_im = si.vgy_im;
  Real vgz_re = si.vgz_re;
  Real vgz_im = si.vgz_im;
  bool ideal = eos.is_ideal;
  Real gm1 = eos.gamma - 1.0;
  Real pres = pin->GetOrAddReal("problem","pressure",1.0);

  par_for("streaming_gas", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                         size.d_view(m).x1max);
    Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                         size.d_view(m).x3max);
    Real phase = kx*x + kz*z;
    Real cph = cos(phase);
    Real sph = sin(phase);
    Real drho = use_eigenmode ? (rhog_re*cph - rhog_im*sph) : cph;
    Real dvx = use_eigenmode ? (vgx_re*cph - vgx_im*sph) : cph;
    Real dvy = use_eigenmode ? (vgy_re*cph - vgy_im*sph) : 0.0;
    Real dvz = use_eigenmode ? (vgz_re*cph - vgz_im*sph) : cph;
    Real rho = rho0*(1.0 + amp*drho);
    Real vx = gvx + amp*dvx;
    Real vy = gvy + amp*dvy;
    Real vz = gvz + amp*dvz;

    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = rho*vx;
    u0(m,IM2,k,j,i) = rho*vy;
    u0(m,IM3,k,j,i) = rho*vz;
    if (ideal) {
      u0(m,IEN,k,j,i) = pres/gm1 + 0.5*rho*(vx*vx + vy*vy + vz*vz);
    }
  });

  if (use_mhd) {
    auto &b0 = pmbp->pmhd->b0;
    par_for("streaming_mhd_bzero", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = 0.0;
      if (i == ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j == je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k == ke) {b0.x3f(m,k+1,j,i) = 0.0;}
    });
  }

  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &pi = pmbp->ppart->prtcl_idata;
  int npart = pmbp->ppart->nprtcl_thispack;
  int nc1 = indcs.nx1;
  int nc2 = indcs.nx2;
  int nc3 = indcs.nx3;
  int ncells = nc1*nc2*nc3;
  int nmb = pmbp->nmb_thispack;
  int gids = pmbp->gids;
  Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
  if (ppc <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "streaming_instability requires particles/ppc > 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_mesh_->nprtcl_total <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "streaming_instability requires at least one particle"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real dust_density = si.epsilon*si.rho0;
  Real pvx = si.part_vx;
  Real pvy = si.part_vy;
  Real pvz = si.part_vz;
  Real rhod_re = si.rhod_re;
  Real rhod_im = si.rhod_im;
  Real vdx_re = si.vdx_re;
  Real vdx_im = si.vdx_im;
  Real vdy_re = si.vdy_re;
  Real vdy_im = si.vdy_im;
  Real vdz_re = si.vdz_re;
  Real vdz_im = si.vdz_im;

  par_for("streaming_particles", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    int local = (ncells > 0) ? p % (nmb*ncells) : 0;
    int m = local/ncells;
    int rem = local - m*ncells;
    int kk = rem/(nc1*nc2);
    rem -= kk*nc1*nc2;
    int jj = rem/nc1;
    int ii = rem - jj*nc1;

    Real x = CellCenterX(ii, nc1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real y = CellCenterX(jj, nc2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real z = CellCenterX(kk, nc3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phase = kx*x + kz*z;
    Real cph = cos(phase);
    Real sph = sin(phase);
    Real drho = use_eigenmode ? (rhod_re*cph - rhod_im*sph) : cph;
    Real dvx = use_eigenmode ? (vdx_re*cph - vdx_im*sph) : cph;
    Real dvy = use_eigenmode ? (vdy_re*cph - vdy_im*sph) : 0.0;
    Real dvz = use_eigenmode ? (vdz_re*cph - vdz_im*sph) : cph;
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real pmass0 = dust_density*vol/ppc;

    pi(PGID,p) = gids + m;
    pr(IPX,p) = x;
    pr(IPY,p) = y;
    pr(IPZ,p) = z;
    pr(IPVX,p) = pvx + amp*dvx;
    pr(IPVY,p) = pvy + amp*dvy;
    pr(IPVZ,p) = pvz + amp*dvz;
    pr(IPM,p) = pmass0*(1.0 + amp*drho);
  });

  Real vmax = std::max(std::abs(si.part_vx), std::abs(si.part_vy));
  vmax = std::max(vmax, std::abs(si.part_vz)) + si.amp;
  if (vmax > 0.0) {
    Real dx = std::min(size.h_view(0).dx1, size.h_view(0).dx3);
    pmbp->ppart->dtnew = std::min(pmbp->ppart->dtnew, 0.3*dx/vmax);
  }
}

//----------------------------------------------------------------------------------------
//! \fn StreamingInstabilitySource

void StreamingInstabilitySource(Mesh *pm, const Real bdt) {
  if (!si.orbital_terms) {return;}
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &w0 = pmbp->phydro->w0;
  auto &u0 = pmbp->phydro->u0;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  Real omega0 = si.omega0;
  Real qshear = si.qshear;
  Real pressure_accel = si.pressure_accel;

  par_for("streaming_orbital_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real rho = w0(m,IDN,k,j,i);
    Real vx = w0(m,IVX,k,j,i);
    Real vy = w0(m,IVY,k,j,i);
    u0(m,IM1,k,j,i) += bdt*rho*(2.0*omega0*vy + pressure_accel);
    u0(m,IM2,k,j,i) += bdt*rho*(-(2.0 - qshear)*omega0*vx);
  });
}

//----------------------------------------------------------------------------------------
//! \fn StreamingInstabilityUserDrag

void StreamingInstabilityUserDrag(MeshBlockPack *pmbp, const Real bdt) {
  auto *ppart = pmbp->ppart;
  if (ppart->nprtcl_thispack <= 0) {return;}

  auto &indcs = pmbp->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  bool multi_d = pmbp->pmesh->multi_d;
  bool three_d = pmbp->pmesh->three_d;
  bool use_mhd = (pmbp->pmhd != nullptr);
  auto w0 = use_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto u0 = use_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  EOS_Data eos_data = use_mhd ? pmbp->pmhd->peos->eos_data
                              : pmbp->phydro->peos->eos_data;
  auto &mbsize = pmbp->pmb->mb_size;
  auto pr = ppart->prtcl_rdata;
  auto pi = ppart->prtcl_idata;
  int npart = ppart->nprtcl_thispack;
  int gids = pmbp->gids;
  int nmb = pmbp->nmb_thispack;
  Real stopping_time = si.user_drag_stopping_time;
  Real default_mass = ppart->drag_particle_mass;
  bool backreaction = ppart->drag_backreaction;
  bool include_energy = ppart->drag_include_energy && eos_data.is_ideal;
  DragParticlesCoupling interpolation = ppart->drag_interpolation;
  DragParticlesCoupling deposition = ppart->drag_deposition;

  par_for("streaming_user_particle_drag", DevExeSpace(), 0, npart-1,
  KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    if ((m < 0) || (m >= nmb)) {return;}
    Real pmass = pr(IPM,p);
    if (pmass <= 0.0) {
      pmass = default_mass;
      pr(IPM,p) = pmass;
    }

    DragParticleStencil interp_stencil;
    BuildDragParticleStencil(interpolation, pr(IPX,p), pr(IPY,p), pr(IPZ,p),
                             mbsize.d_view(m).x1min, mbsize.d_view(m).x2min,
                             mbsize.d_view(m).x3min, mbsize.d_view(m).dx1,
                             mbsize.d_view(m).dx2, mbsize.d_view(m).dx3, is, ie, js,
                             je, ks, ke, multi_d, three_d, interp_stencil);
    Real vg1, vg2, vg3;
    InterpolateDragVelocity(w0, m, interp_stencil, multi_d, three_d, vg1, vg2, vg3);

    Real vx_old = pr(IPVX,p);
    Real vy_old = multi_d ? pr(IPVY,p) : 0.0;
    Real vz_old = three_d ? pr(IPVZ,p) : 0.0;
    Real drag_factor = 1.0 - exp(-bdt/stopping_time);
    Real dv1_drag = drag_factor*(vg1 - vx_old);
    Real dv2_drag = multi_d ? drag_factor*(vg2 - vy_old) : 0.0;
    Real dv3_drag = three_d ? drag_factor*(vg3 - vz_old) : 0.0;
    Real vx_new = vx_old + dv1_drag;
    Real vy_new = vy_old + dv2_drag;
    Real vz_new = vz_old + dv3_drag;

    pr(IPVX,p) = vx_new;
    if (multi_d) {pr(IPVY,p) = vy_new;}
    if (three_d) {pr(IPVZ,p) = vz_new;}

    if (backreaction) {
      Real vol = mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*mbsize.d_view(m).dx3;
      Real denergy = 0.0;
      if (include_energy) {
        Real ke_old = 0.5*pmass*(vx_old*vx_old + vy_old*vy_old + vz_old*vz_old);
        Real ke_new = 0.5*pmass*(vx_new*vx_new + vy_new*vy_new + vz_new*vz_new);
        denergy = -(ke_new - ke_old);
      }
      DragParticleStencil deposit_stencil;
      BuildDragParticleStencil(deposition, pr(IPX,p), pr(IPY,p), pr(IPZ,p),
                               mbsize.d_view(m).x1min, mbsize.d_view(m).x2min,
                               mbsize.d_view(m).x3min, mbsize.d_view(m).dx1,
                               mbsize.d_view(m).dx2, mbsize.d_view(m).dx3, is, ie, js,
                               je, ks, ke, multi_d, three_d, deposit_stencil);
      DepositDragBackreaction(u0, m, deposit_stencil, 1.0/vol, -pmass*dv1_drag,
                              -pmass*dv2_drag, -pmass*dv3_drag, denergy, multi_d,
                              three_d, include_energy);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn StreamingInstabilityHistory

void StreamingInstabilityHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 14;
  pdata->label[0] = "gmass";
  pdata->label[1] = "pmass";
  pdata->label[2] = "gmom1";
  pdata->label[3] = "pmom1";
  pdata->label[4] = "gmom2";
  pdata->label[5] = "pmom2";
  pdata->label[6] = "gmom3";
  pdata->label[7] = "pmom3";
  pdata->label[8] = "gmodec";
  pdata->label[9] = "gmodes";
  pdata->label[10] = "pmodec";
  pdata->label[11] = "pmodes";
  pdata->label[12] = "gtotE";
  pdata->label[13] = "pkinE";

  MeshBlockPack *pmbp = pm->pmb_pack;
  bool use_mhd = (pmbp->pmhd != nullptr);
  auto u0 = use_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  EOS_Data eos = use_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  auto &pr = pmbp->ppart->prtcl_rdata;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  int nkji = nx3*nx2*nx1;
  int nji = nx2*nx1;
  Real rho0 = si.rho0;
  Real kx = si.kx;
  Real kz = si.kz;
  bool ideal = eos.is_ideal;
  array_sum::GlobalSum gas_sum;

  Kokkos::parallel_reduce("streaming_gas_history",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phase = kx*x + kz*z;
    sum.the_array[0] += vol*u0(m,IDN,k,j,i);
    sum.the_array[2] += vol*u0(m,IM1,k,j,i);
    sum.the_array[4] += vol*u0(m,IM2,k,j,i);
    sum.the_array[6] += vol*u0(m,IM3,k,j,i);
    sum.the_array[8] += vol*(u0(m,IDN,k,j,i) - rho0)*cos(phase);
    sum.the_array[9] += vol*(u0(m,IDN,k,j,i) - rho0)*sin(phase);
    if (ideal) {sum.the_array[12] += vol*u0(m,IEN,k,j,i);}
  }, Kokkos::Sum<array_sum::GlobalSum>(gas_sum));

  int npart = pmbp->ppart->nprtcl_thispack;
  array_sum::GlobalSum part_sum;
  Kokkos::parallel_reduce("streaming_particle_history",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, npart),
  KOKKOS_LAMBDA(const int &p, array_sum::GlobalSum &sum) {
    Real phase = kx*pr(IPX,p) + kz*pr(IPZ,p);
    Real mass = pr(IPM,p);
    sum.the_array[1] += mass;
    sum.the_array[3] += mass*pr(IPVX,p);
    sum.the_array[5] += mass*pr(IPVY,p);
    sum.the_array[7] += mass*pr(IPVZ,p);
    sum.the_array[10] += mass*cos(phase);
    sum.the_array[11] += mass*sin(phase);
    sum.the_array[13] += 0.5*mass*(pr(IPVX,p)*pr(IPVX,p) + pr(IPVY,p)*pr(IPVY,p) +
                                    pr(IPVZ,p)*pr(IPVZ,p));
  }, Kokkos::Sum<array_sum::GlobalSum>(part_sum));

  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = gas_sum.the_array[n] + part_sum.the_array[n];
  }
}
