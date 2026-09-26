//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistive_tests.cpp
//! \brief Magnetic diffusion and single-sheet reconnection initial conditions.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "diffusion/current_limited_resistivity.hpp"
#include "diffusion/resistivity.hpp"
#include "eos/eos.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

namespace {

enum class ResistiveTest {forcefree, gaussian, sheet, harris};

current_limited::Parameters history_params;
bool history_harris = false;
Real history_xc, history_yc, history_b0, history_density;

[[noreturn]] void ResistiveTestFatal(const std::string &message) {
  std::cout << "### FATAL ERROR in resistive_tests: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real Sinc(const Real x) {
  return x == 0.0 ? 1.0 : sin(x)/x;
}

// Stable log(cosh(x)) also far outside a thin current sheet.
KOKKOS_INLINE_FUNCTION
Real LogCosh(const Real x) {
  return fabs(x) + log(1.0 + exp(-2.0*fabs(x))) - log(2.0);
}

KOKKOS_INLINE_FUNCTION
Real HarrisAz(const Real x, const Real y, const Real b0, const Real width,
              const Real flux, const Real perturb_width, const bool gem,
              const Real lx, const Real ly) {
  Real perturb = gem ? cos(2.0*M_PI*x/lx)*cos(M_PI*y/ly)
                     : exp(-(x*x + y*y)/(perturb_width*perturb_width));
  return b0*width*LogCosh(y/width) + flux*perturb;
}

void ResistiveHistory(HistoryData *pdata, Mesh *pm);

} // namespace

void ProblemGenerator::ResistiveTests(ParameterInput *pin, const bool restart) {
  auto *pack = pmy_mesh_->pmb_pack;
  if (pack->pmhd == nullptr || pack->pcoord->is_special_relativistic ||
      pack->pcoord->is_general_relativistic || pack->pcoord->is_dynamical_relativistic) {
    ResistiveTestFatal("requires Newtonian MHD");
  }
  if (!pack->pmhd->peos->eos_data.is_ideal) {
    ResistiveTestFatal("requires ideal-gas MHD to track magnetic heating");
  }
  const std::string name = pin->GetOrAddString("problem", "test", "forcefree");
  ResistiveTest test;
  if (name == "forcefree") test = ResistiveTest::forcefree;
  else if (name == "gaussian") test = ResistiveTest::gaussian;
  else if (name == "sheet") test = ResistiveTest::sheet;
  else if (name == "harris") test = ResistiveTest::harris;
  else ResistiveTestFatal("test must be forcefree, gaussian, sheet, or harris");

  const Real density = pin->GetOrAddReal("problem", "density", 1.0);
  const Real pressure = pin->GetOrAddReal("problem", "pressure", 1.0);
  const Real bamp = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real guide = pin->GetOrAddReal("problem", "guide_field", 0.0);
  const Real width = pin->GetOrAddReal("problem", "width", 0.1);
  const int n1 = pin->GetOrAddInteger("problem", "wave_n1", 1);
  const int n2 = pin->GetOrAddInteger("problem", "wave_n2", 0);
  const int n3 = pin->GetOrAddInteger("problem", "wave_n3", 0);
  const Real sheet_density = pin->GetOrAddReal("problem", "sheet_density", 1.0);
  const Real flux = pin->GetOrAddReal("problem", "perturbation_flux", 0.0);
  const Real perturb_width = pin->GetOrAddReal("problem", "perturbation_width", 0.5);
  const std::string perturb = pin->GetOrAddString("problem", "perturbation", "gaussian");
  const bool gem = perturb == "gem";
  if (!(density > 0.0) || !(pressure > 0.0) || !(bamp > 0.0) || !(width > 0.0) ||
      !(sheet_density >= 0.0) || !(perturb_width > 0.0) || !std::isfinite(density) ||
      !std::isfinite(pressure) || !std::isfinite(bamp) || !std::isfinite(width) ||
      !std::isfinite(sheet_density) || !std::isfinite(perturb_width) ||
      !std::isfinite(guide) || !std::isfinite(flux)) {
    ResistiveTestFatal("density, pressure, b0, widths must be positive; "
                       "sheet_density nonnegative and all field parameters finite");
  }
  if (perturb != "gaussian" && !gem) ResistiveTestFatal("unknown perturbation profile");
  if (test == ResistiveTest::forcefree &&
      ((n1 == 0 && n2 == 0 && n3 == 0) || (!pmy_mesh_->multi_d && n2 != 0) ||
       (!pmy_mesh_->three_d && n3 != 0) || guide != 0.0)) {
    ResistiveTestFatal("forcefree needs a nonzero resolved wavevector and zero guide_field");
  }
  if (test == ResistiveTest::harris &&
      (!pmy_mesh_->two_d || pmy_mesh_->mesh_indcs.nx1 % 2 != 0 ||
       pmy_mesh_->mesh_indcs.nx2 % 2 != 0)) {
    ResistiveTestFatal("harris requires 2D and even root-grid sizes so its central X "
                       "lies on a native electric-field edge");
  }
  auto &domain = pmy_mesh_->mesh_size;
  const Real xc = 0.5*(domain.x1min + domain.x1max);
  const Real yc = 0.5*(domain.x2min + domain.x2max);
  const Real zc = 0.5*(domain.x3min + domain.x3max);
  const Real lx = domain.x1max - domain.x1min;
  const Real ly = domain.x2max - domain.x2min;
  const Real kx = 2.0*M_PI*n1/lx;
  const Real ky = 2.0*M_PI*n2/ly;
  const Real kz = 2.0*M_PI*n3/(domain.x3max - domain.x3min);
  const Real kxy = sqrt(kx*kx + ky*ky);
  const Real knorm = sqrt(kx*kx + ky*ky + kz*kz);
  // Orthonormal polarization vectors, with e2 = khat cross e1.
  const Real e1x = kxy > 0.0 ? -ky/kxy : 1.0;
  const Real e1y = kxy > 0.0 ? kx/kxy : 0.0;
  const Real e2x = knorm > 0.0 ? -kz*e1y/knorm : 0.0;
  const Real e2y = knorm > 0.0 ? kz*e1x/knorm : 0.0;
  const Real e2z = knorm > 0.0 ? (kx*e1y-ky*e1x)/knorm : 0.0;

  // Re-enroll on restart. Histories use the same ion interpretation and fixed b_rec
  // as the closure. In constant/no-resistivity controls, d_i and b_rec affect only q.
  auto *resist = pack->pmhd->presist;
  if (resist != nullptr && resist->current_limited) {
    history_params = resist->current_limited_params;
  } else {
    history_params.eta0 = resist == nullptr ? 0.0 : resist->eta_ohm;
    history_params.eta_max = history_params.eta0;
    history_params.eta_star = history_params.eta0;
    history_params.q_star = 0.0;
    history_params.d_i = pin->GetOrAddReal("mhd", "d_i", 1.0);
    history_params.b_rec = pin->GetOrAddReal("mhd", "b_rec", bamp);
  }
  if (!(history_params.d_i > 0.0) || !(history_params.b_rec > 0.0)) {
    ResistiveTestFatal("diagnostic d_i and b_rec must be positive");
  }
  history_harris = test == ResistiveTest::harris;
  history_xc = xc;
  history_yc = yc;
  history_b0 = bamp;
  history_density = density;
  user_hist_func = ResistiveHistory;
  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto size = pack->pmb->mb_size;
  auto b = pack->pmhd->b0;
  const int nmb = pack->nmb_thispack;

  // Analytic face averages of the force-free wave are a discrete curl of its
  // line-averaged vector potential. This uses K_i=2 sin(k_i dx_i/2)/dx_i and
  // preserves discrete divergence, including coarse/fine flux matching.
  par_for("resistive_pgen_b1", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const auto s = size.d_view(m);
    const Real x = LeftEdgeX(i-is,nx1,s.x1min,s.x1max) - xc;
    const Real yf = LeftEdgeX(j-js,nx2,s.x2min,s.x2max) - yc;
    const Real y = CellCenterX(j-js,nx2,s.x2min,s.x2max) - yc;
    Real bx = 0.0;
    if (test == ResistiveTest::forcefree) {
      const Real z = CellCenterX(k-ks,nx3,s.x3min,s.x3max) - zc;
      const Real phase = kx*x + ky*y + kz*z;
      bx = bamp*Sinc(0.5*ky*s.dx2)*Sinc(0.5*kz*s.dx3)*
           (e1x*sin(phase) + e2x*cos(phase));
    } else if (test == ResistiveTest::harris) {
      bx = (HarrisAz(x,yf+s.dx2,bamp,width,flux,perturb_width,gem,lx,ly) -
            HarrisAz(x,yf,bamp,width,flux,perturb_width,gem,lx,ly))/s.dx2;
    }
    b.x1f(m,k,j,i) = bx;
  });
  par_for("resistive_pgen_b2", DevExeSpace(), 0,nmb-1,ks,ke,js,je+1,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const auto s = size.d_view(m);
    const Real xf = LeftEdgeX(i-is,nx1,s.x1min,s.x1max) - xc;
    const Real x = CellCenterX(i-is,nx1,s.x1min,s.x1max) - xc;
    const Real y = LeftEdgeX(j-js,nx2,s.x2min,s.x2max) - yc;
    Real by;
    if (test == ResistiveTest::forcefree) {
      const Real z = CellCenterX(k-ks,nx3,s.x3min,s.x3max) - zc;
      const Real phase = kx*x + ky*y + kz*z;
      by = bamp*Sinc(0.5*kx*s.dx1)*Sinc(0.5*kz*s.dx3)*
           (e1y*sin(phase) + e2y*cos(phase));
    } else if (test == ResistiveTest::gaussian) {
      by = bamp*exp(-0.5*SQR(x/width));
    } else if (test == ResistiveTest::sheet) {
      by = bamp*tanh(x/width);
    } else {
      by = -(HarrisAz(xf+s.dx1,y,bamp,width,flux,perturb_width,gem,lx,ly) -
             HarrisAz(xf,y,bamp,width,flux,perturb_width,gem,lx,ly))/s.dx1;
    }
    b.x2f(m,k,j,i) = by;
  });
  par_for("resistive_pgen_b3", DevExeSpace(), 0,nmb-1,ks,ke+1,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const auto s = size.d_view(m);
    const Real x = CellCenterX(i-is,nx1,s.x1min,s.x1max) - xc;
    const Real y = CellCenterX(j-js,nx2,s.x2min,s.x2max) - yc;
    const Real z = LeftEdgeX(k-ks,nx3,s.x3min,s.x3max) - zc;
    b.x3f(m,k,j,i) = test == ResistiveTest::forcefree ?
        bamp*Sinc(0.5*kx*s.dx1)*Sinc(0.5*ky*s.dx2)*e2z*cos(kx*x + ky*y + kz*z)
        : guide;
  });

  auto u = pack->pmhd->u0;
  const Real gm1 = pack->pmhd->peos->eos_data.gamma - 1.0;
  const int nscalars = pack->pmhd->nscalars;
  const int nmhd = pack->pmhd->nmhd;
  par_for("resistive_pgen_u", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const auto s = size.d_view(m);
    const Real y = CellCenterX(j-js,nx2,s.x2min,s.x2max) - yc;
    const Real bx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
    const Real by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
    const Real bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
    Real rho = density, pgas = pressure;
    if (test == ResistiveTest::sheet) {
      pgas += 0.5*(bamp*bamp - by*by);
    } else if (test == ResistiveTest::harris) {
      const Real t = tanh(y/width);
      rho += sheet_density*(1.0 - t*t);
      const Real yf = LeftEdgeX(j-js,nx2,s.x2min,s.x2max) - yc;
      const Real bx_unperturbed = bamp*width*
          (LogCosh((yf+s.dx2)/width) - LogCosh(yf/width))/s.dx2;
      pgas += 0.5*(bamp*bamp - bx_unperturbed*bx_unperturbed);
    }
    u(m,IDN,k,j,i) = rho;
    u(m,IM1,k,j,i) = 0.0;
    u(m,IM2,k,j,i) = 0.0;
    u(m,IM3,k,j,i) = 0.0;
    u(m,IEN,k,j,i) = pgas/gm1 + 0.5*(bx*bx + by*by + bz*bz);
    for (int n=0; n<nscalars; ++n) u(m,nmhd+n,k,j,i) = 0.0;
  });
}

namespace {

// q/eta maxima use actual closure edges. Fractions and eta J^2 are cell-centered
// reconstructed proxies, not the exact discrete magnetic-energy loss. The central
// X diagnostics assume the symmetric Harris X remains at the domain midpoint.
// Reference diagnostics use the positive-x boundary at the same y, not an O point.
void ResistiveHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = history_harris ? 12 : 6;
  const char *labels[] = {"q_max_edge", "eta_max_ed", "frac_qstar", "frac_q1",
                         "etaJ2_cc", "heat_frac", "x_q", "x_etaJ", "x_etaJz", "x_Ez",
                         "ref_etaJz", "ref_Ez"};
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->label[n] = labels[n];
    pdata->hdata[n] = 0.0;
  }
  auto *pack = pm->pmb_pack;
  auto b = pack->pmhd->b0;
  auto w = pack->pmhd->w0;
  auto size = pack->pmb->mb_size;
  const auto indcs = pm->mb_indcs;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nji = nx2*nx1, nkji = nx3*nji;
  const bool multi_d = pm->multi_d, three_d = pm->three_d;
  const bool harris = history_harris;
  auto *resist = pack->pmhd->presist;
  const auto params = resist != nullptr && resist->current_limited ?
      resist->current_limited_params : history_params;
  const Real dfloor = pack->pmhd->peos->eos_data.dfloor;
  const Real xc = history_xc, yc = history_yc, xref = pm->mesh_size.x1max;
  const Real enorm = history_b0*history_b0/sqrt(history_density);
  array_sum::GlobalSum sums;
  Real maxq = 0.0, maxeta = 0.0;
  Kokkos::parallel_reduce("resistive_history",
      Kokkos::RangePolicy<>(DevExeSpace(),0,pack->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum, Real &qmax, Real &etamax) {
    const int m = idx/nkji;
    const int k = (idx % nkji)/nji + indcs.ks;
    const int j = (idx % nji)/nx1 + indcs.js;
    const int i = idx % nx1 + indcs.is;
    const auto s = size.d_view(m);
    const Real vol = s.dx1*s.dx2*s.dx3;
    const auto cell = current_limited::CellState(b,w,s,params,dfloor,
                                               multi_d,three_d,m,k,j,i);
    const Real diss = cell.eta*(SQR(cell.j1) + SQR(cell.j2) + SQR(cell.j3))*vol;
    array_sum::GlobalSum local;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) local.the_array[n] = 0.0;
    local.the_array[0] = vol;
    local.the_array[1] = cell.q > params.q_star ? vol : 0.0;
    local.the_array[2] = cell.q > 1.0 ? vol : 0.0;
    local.the_array[3] = diss;
    local.the_array[4] = cell.q > params.q_star ? diss : 0.0;
    for (int c=0; c<3; ++c) {
      if (c == 0 && !multi_d) continue;
      for (int dk=0; dk<=(three_d && c!=2); ++dk) {
        for (int dj=0; dj<=(multi_d && c!=1); ++dj) {
          for (int di=0; di<=(c!=0); ++di) {
            const auto edge = current_limited::EdgeState(b,w,s,params,dfloor,
                multi_d,three_d,c,m,k+dk,j+dj,i+di);
            qmax = fmax(qmax,edge.q);
            etamax = fmax(etamax,edge.eta);
          }
        }
      }
    }
    if (harris) {
      const Real x = LeftEdgeX(i-indcs.is,nx1,s.x1min,s.x1max);
      const Real y = LeftEdgeX(j-indcs.js,nx2,s.x2min,s.x2max);
      for (int point=0; point<2; ++point) {
        const int ei = i + point;
        const bool at_point = point == 0 ? fabs(x-xc) < 1.0e-10*s.dx1 :
            i == indcs.ie && fabs(x+s.dx1-xref) < 1.0e-10*s.dx1;
        if (!at_point || fabs(y-yc) >= 1.0e-10*s.dx2) continue;
        const auto edge = current_limited::EdgeState(b,w,s,params,dfloor,
            multi_d,three_d,2,m,k,j,ei);
        Real vx = 0.0, vy = 0.0;
        for (int dj=-1; dj<=0; ++dj) {
          for (int di=-1; di<=0; ++di) {
            vx += 0.25*w(m,IVX,k,j+dj,ei+di);
            vy += 0.25*w(m,IVY,k,j+dj,ei+di);
          }
        }
        const Real bx = 0.5*(b.x1f(m,k,j-1,ei) + b.x1f(m,k,j,ei));
        const Real by = 0.5*(b.x2f(m,k,j,ei-1) + b.x2f(m,k,j,ei));
        const int offset = point == 0 ? 7 : 10;
        local.the_array[offset] = edge.eta*edge.j3/enorm;
        local.the_array[offset+1] = (vy*bx-vx*by+edge.eta*edge.j3)/enorm;
        local.the_array[offset+2] = 1.0;
        if (point == 0) {
          local.the_array[5] = edge.q;
          local.the_array[6] = edge.eta*sqrt(SQR(edge.j1)+SQR(edge.j2)+SQR(edge.j3))/enorm;
        }
      }
    }
    sum += local;
  }, Kokkos::Sum<array_sum::GlobalSum>(sums), Kokkos::Max<Real>(maxq),
     Kokkos::Max<Real>(maxeta));

  Real maxima[2] = {maxq,maxeta};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE,sums.the_array,13,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,maxima,2,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_WORLD);
#endif
  if (harris && (sums.the_array[9] == 0.0 || sums.the_array[12] == 0.0)) {
    ResistiveTestFatal("central X or boundary reference edge was not found for history output");
  }
  // HistoryOutput performs MPI_SUM afterward; publish global ratios/maxima on rank 0 only.
  if (global_variable::my_rank != 0) return;
  pdata->hdata[0] = maxima[0];
  pdata->hdata[1] = maxima[1];
  pdata->hdata[2] = sums.the_array[1]/sums.the_array[0];
  pdata->hdata[3] = sums.the_array[2]/sums.the_array[0];
  pdata->hdata[4] = sums.the_array[3];
  pdata->hdata[5] = sums.the_array[3] > 0.0 ? sums.the_array[4]/sums.the_array[3] : 0.0;
  if (harris && sums.the_array[9] > 0.0) {
    for (int n=6; n<10; ++n) pdata->hdata[n] = sums.the_array[n-1]/sums.the_array[9];
    for (int n=10; n<12; ++n) pdata->hdata[n] = sums.the_array[n]/sums.the_array[12];
  }
}

} // namespace
