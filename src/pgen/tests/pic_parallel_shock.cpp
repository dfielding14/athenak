//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pic_parallel_shock.cpp
//! \brief Parallel-shock MHD-PIC benchmark problem generator.
//!
//! Minimal AthenaK adaptation of the Section 5.4 parallel-shock setup in:
//!   Sun & Bai (2022), "The MHD-PIC Module in Athena++"
//! using:
//! - reflecting wall at inner x1
//! - inflow at outer x1
//! - background field B0 parallel to x
//! - eta-based shock-surface CR injection with conservative gas subtraction
//! - AMR refinement based on density/pressure curvature metrics.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "pgen/pgen.hpp"

namespace {

struct ShockCell {
  int m, k, j, i;
  int gid;
  Real x1c, x2c, x3c;
  Real dx1, dx2, dx3;
  Real area;
  Real vol;
};

struct InjectedParticle {
  int gid;
  int m, k, j, i;
  Real vol;
  Real x1, x2, x3;
  Real vx, vy, vz;
};

struct GasDelta {
  int m, k, j, i;
  Real dm;
  Real dmx;
  Real dmy;
  Real dmz;
  Real de;
};

// Runtime controls populated by ProblemGenerator::PICParallelShock().
Real ps_rho0 = 1.0;
Real ps_p0 = 1.0;
Real ps_u0 = 10.0;
Real ps_b0 = 1.0;
Real ps_eta = 1.0e-3;
Real ps_vinj_over_u0 = std::sqrt(10.0);
Real ps_shock_speed = 0.0;
Real ps_xshock0 = 0.0;
Real ps_inject_half_width_cells = 0.5;
Real ps_inject_t_start = 0.0;
Real ps_inject_t_stop = 1.0e99;
Real ps_seed_noise_amp = 0.0;
int ps_seed_noise_seed = 1234;
std::array<Real, 4> ps_seed_noise_phase_by = {{0.0, 0.0, 0.0, 0.0}};
std::array<Real, 4> ps_seed_noise_phase_bz = {{0.0, 0.0, 0.0, 0.0}};
Real ps_refine_curv = 1.0;
Real ps_derefine_curv = 0.1;
Real ps_rho_floor_frac = 1.0e-6;
Real ps_p_floor_frac = 1.0e-8;
bool ps_enable_injection = true;
bool ps_enable_subtraction = true;
bool ps_enable_curvature_amr = true;
bool ps_enable_frame_tracking = false;
bool ps_use_2d3v = false;
enum class PSFrameMode { velocity, recenter };
PSFrameMode ps_frame_mode = PSFrameMode::velocity;
Real ps_frame_t_start = 0.0;
Real ps_frame_t_ramp = 0.0;
Real ps_frame_vfrac = 1.0;
Real ps_frame_dv_max = 1.0e99;
bool ps_frame_apply_to_particles = true;
bool ps_frame_apply_to_inflow = true;
bool ps_frame_require_uniform = true;
int ps_frame_diag_dcycle = 200;
int ps_feedback_diag_dcycle = 0;
Real ps_recenter_x_target = 2.0;
Real ps_recenter_x_trigger = 3.0;
Real ps_recenter_dx1 = 0.0;
int ps_recenter_shift_cells = 2;
Real ps_recenter_vshock_model = -1.0;
int ps_inject_species = 0;
int ps_inject_seed = 1234;
Real ps_particle_mass = 1.0;
Real ps_particle_charge = 1.0;
Real ps_particle_q_over_m = 1.0;
Real ps_particle_macro_mass = 1.0;
Real ps_mass_reservoir_local = 0.0;
bool ps_tag_seeded = false;
std::int64_t ps_next_tag = 0;

inline bool FrameModeVelocity() {
  return ps_frame_mode == PSFrameMode::velocity;
}

inline bool FrameModeRecenter() {
  return ps_frame_mode == PSFrameMode::recenter;
}

inline Real FrameRampFactor(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  if (!ps_enable_frame_tracking) return 0.0;
  if (t <= ps_frame_t_start) return 0.0;
  if (ps_frame_t_ramp <= 0.0) return 1.0;

  const Real s = (t - ps_frame_t_start)/ps_frame_t_ramp;
  if (s >= 1.0) return 1.0;
  if (s <= 0.0) return 0.0;
  return 0.5*(1.0 - std::cos(M_PI*s));
}

inline Real FrameVelocityOffset(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  const Real vfollow = ps_frame_vfrac*ps_shock_speed;
  return -vfollow*FrameRampFactor(t);
}

inline Real FrameOffsetDisplacement(const Real t) {
  if (!FrameModeVelocity()) return 0.0;
  if (!ps_enable_frame_tracking) return 0.0;
  if (t <= ps_frame_t_start) return 0.0;

  const Real vfollow = ps_frame_vfrac*ps_shock_speed;
  if (ps_frame_t_ramp <= 0.0) {
    return -vfollow*(t - ps_frame_t_start);
  }

  const Real tr = ps_frame_t_ramp;
  const Real dt = t - ps_frame_t_start;
  if (dt >= tr) {
    const Real iramp = 0.5*tr;
    const Real itail = dt - tr;
    return -vfollow*(iramp + itail);
  }

  const Real s = dt/tr;
  const Real iramp = tr*(0.5*s - std::sin(M_PI*s)/(2.0*M_PI));
  return -vfollow*iramp;
}

inline Real ShockSurfaceUnshiftedX1(const Real t) {
  if (FrameModeRecenter()) {
    Real vmodel = ps_shock_speed;
    if (ps_recenter_vshock_model > 0.0) vmodel = ps_recenter_vshock_model;
    return ps_xshock0 + vmodel*t;
  }
  return ps_xshock0 + ps_shock_speed*t + FrameOffsetDisplacement(t);
}

inline int RecenterEventCount(const Real t) {
  if (!ps_enable_frame_tracking) return 0;
  if (!FrameModeRecenter()) return 0;
  if (ps_recenter_shift_cells < 1 || ps_recenter_dx1 <= 0.0) return 0;

  const Real xu = ShockSurfaceUnshiftedX1(t);
  const Real dx_shift = ps_recenter_shift_cells*ps_recenter_dx1;
  const Real trigger = ps_recenter_x_trigger;
  const Real target = ps_recenter_x_target;
  if (xu <= trigger) return 0;
  const Real overshoot = xu - target;
  if (overshoot <= 0.0) return 0;
  return static_cast<int>(std::ceil(overshoot/dx_shift));
}

inline Real RecenterDisplacement(const Real t) {
  if (!ps_enable_frame_tracking) return 0.0;
  if (!FrameModeRecenter()) return 0.0;
  const int nevents = RecenterEventCount(t);
  return static_cast<Real>(nevents*ps_recenter_shift_cells)*ps_recenter_dx1;
}

inline Real ShockSurfaceModelX1(const Real t) {
  return ShockSurfaceUnshiftedX1(t) - RecenterDisplacement(t);
}

inline std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline Real UniformFromUint64(std::uint64_t x) {
  constexpr Real inv = 1.0/static_cast<Real>(std::numeric_limits<std::uint64_t>::max());
  return static_cast<Real>(x)*inv;
}

inline void ConfigureSeedNoisePhases() {
  std::uint64_t state = static_cast<std::uint64_t>(ps_seed_noise_seed);
  for (int n = 0; n < 4; ++n) {
    state = SplitMix64(state + 0x9e3779b97f4a7c15ULL);
    ps_seed_noise_phase_by[n] = 2.0*M_PI*UniformFromUint64(state);
    state = SplitMix64(state + 0xbf58476d1ce4e5b9ULL);
    ps_seed_noise_phase_bz[n] = 2.0*M_PI*UniformFromUint64(state);
  }
}

inline Real EstimateShockSpeed(const Real gamma, const Real rho0, const Real p0,
                               const Real u0) {
  // Finite-Mach hydrodynamic estimate for a piston-driven shock reflected off
  // the inner wall. For B parallel to shock normal (Bx-only), magnetic terms
  // do not change the 1D compression ratio in this minimal setup.
  const Real cs2 = gamma*p0/rho0;
  if (cs2 <= 0.0) return 0.0;
  const Real ms2 = SQR(u0)/cs2;
  if (ms2 <= 1.0) return 0.0;
  const Real r = ((gamma + 1.0)*ms2)/((gamma - 1.0)*ms2 + 2.0);
  if (r <= 1.0) return 0.0;
  return u0/(r - 1.0);
}

void ApplyFrameShiftToFluid(Mesh *pm, const Real dvx) {
  if (dvx == 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto *pmhd = pmbp->pmhd;
  auto &u0 = pmhd->u0;
  auto &w0 = pmhd->w0;
  auto &b0 = pmhd->b0;
  auto &bcc0 = pmhd->bcc0;

  pmhd->peos->ConsToPrim(u0, b0, w0, bcc0, false, is, ie, js, je, ks, ke);
  par_for("ps_frame_shift_w0", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w0(m, IVX, k, j, i) += dvx;
  });
  pmhd->peos->PrimToCons(w0, bcc0, u0, is, ie, js, je, ks, ke);
}

void ApplyFrameShiftToParticles(Mesh *pm, const Real dvx) {
  if (dvx == 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;

  auto *ppart = pmbp->ppart;
  if (ppart->nprtcl_thispack <= 0) return;
  auto &prtcl_rdata = ppart->prtcl_rdata;
  par_for("ps_frame_shift_particles", DevExeSpace(), 0, ppart->nprtcl_thispack - 1,
  KOKKOS_LAMBDA(const int p) {
    prtcl_rdata(IPVX, p) += dvx;
  });
}

void UpdateOuterInflowState(Mesh *pm, const Real frame_vx) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;
  auto &u_in = pmbp->pmhd->pbval_u->u_in;
  auto &b_in = pmbp->pmhd->pbval_b->b_in;
  const Real gamma = pmbp->pmhd->peos->eos_data.gamma;
  const Real gm1 = gamma - 1.0;
  const Real vin = -ps_u0 + frame_vx;
  const Real ein = ps_p0/gm1 + 0.5*ps_rho0*SQR(vin) + 0.5*SQR(ps_b0);

  u_in.h_view(IDN, BoundaryFace::outer_x1) = ps_rho0;
  u_in.h_view(IM1, BoundaryFace::outer_x1) = ps_rho0*vin;
  u_in.h_view(IM2, BoundaryFace::outer_x1) = 0.0;
  u_in.h_view(IM3, BoundaryFace::outer_x1) = 0.0;
  u_in.h_view(IEN, BoundaryFace::outer_x1) = ein;
  b_in.h_view(IBX, BoundaryFace::outer_x1) = ps_b0;
  b_in.h_view(IBY, BoundaryFace::outer_x1) = 0.0;
  b_in.h_view(IBZ, BoundaryFace::outer_x1) = 0.0;
  u_in.template modify<HostMemSpace>();
  u_in.template sync<DevExeSpace>();
  b_in.template modify<HostMemSpace>();
  b_in.template sync<DevExeSpace>();
}

void ApplyRecenteringShiftToParticles(Mesh *pm, const Real xshift) {
  if (xshift <= 0.0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;

  auto *ppart = pmbp->ppart;
  const int np_old = ppart->nprtcl_thispack;
  if (np_old <= 0) return;

  const int nr = ppart->nrdata;
  const int ni = ppart->nidata;
  auto h_pr_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_rdata);
  auto h_pi_old = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       ppart->prtcl_idata);
  HostArray2D<Real> h_pr_new("ps_pr_recent", nr, np_old);
  HostArray2D<int> h_pi_new("ps_pi_recent", ni, np_old);

  const Real xmin = pm->mesh_size.x1min;
  const Real xmax = pm->mesh_size.x1max;
  int np_new = 0;
  for (int p = 0; p < np_old; ++p) {
    const Real xnew = h_pr_old(IPX, p) - xshift;
    if (!(xnew > xmin && xnew < xmax)) continue;

    for (int q = 0; q < nr; ++q) h_pr_new(q, np_new) = h_pr_old(q, p);
    for (int q = 0; q < ni; ++q) h_pi_new(q, np_new) = h_pi_old(q, p);
    h_pr_new(IPX, np_new) = xnew;
    ++np_new;
  }

  Kokkos::resize(ppart->prtcl_rdata, nr, np_new);
  Kokkos::resize(ppart->prtcl_idata, ni, np_new);
  if (np_new > 0) {
    auto r_sub = Kokkos::subview(ppart->prtcl_rdata, Kokkos::ALL,
                                 std::make_pair(0, np_new));
    auto i_sub = Kokkos::subview(ppart->prtcl_idata, Kokkos::ALL,
                                 std::make_pair(0, np_new));
    auto r_src = Kokkos::subview(h_pr_new, Kokkos::ALL,
                                 std::make_pair(0, np_new));
    auto i_src = Kokkos::subview(h_pi_new, Kokkos::ALL,
                                 std::make_pair(0, np_new));
    Kokkos::deep_copy(r_sub, r_src);
    Kokkos::deep_copy(i_sub, i_src);
  }
  ppart->nprtcl_thispack = np_new;
}

void ApplyRecenteringShift(Mesh *pm, const int nshift) {
  if (nshift <= 0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ng = indcs.ng;
  if (nshift > ng) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock recenter shift must be <= nghost."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto *pmhd = pmbp->pmhd;
  auto h_u0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->u0);
  auto h_x1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x1f);
  auto h_x2f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x2f);
  auto h_x3f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmhd->b0.x3f);
  const int nvar = static_cast<int>(pmhd->u0.extent_int(1));
  const Real vin = -ps_u0 + FrameVelocityOffset(pm->time + pm->dt);
  const Real gm1 = pmhd->peos->eos_data.gamma - 1.0;
  const Real ein = ps_p0/gm1 + 0.5*ps_rho0*SQR(vin) + 0.5*SQR(ps_b0);

  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          if (src <= ie + ng) {
            for (int n = 0; n < nvar; ++n) {
              h_u0(m, n, k, j, i) = h_u0(m, n, k, j, src);
            }
          } else {
            for (int n = 0; n < nvar; ++n) {
              h_u0(m, n, k, j, i) = h_u0(m, n, k, j, ie);
            }
            h_u0(m, IDN, k, j, i) = ps_rho0;
            h_u0(m, IM1, k, j, i) = ps_rho0*vin;
            h_u0(m, IM2, k, j, i) = 0.0;
            h_u0(m, IM3, k, j, i) = 0.0;
            h_u0(m, IEN, k, j, i) = ein;
          }
        }
      }
    }
  }

  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie + 1; ++i) {
          const int src = i + nshift;
          h_x1f(m, k, j, i) = (src <= ie + 1 + ng) ?
              h_x1f(m, k, j, src) : ps_b0;
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je + 1; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          h_x2f(m, k, j, i) = (src <= ie + ng) ? h_x2f(m, k, j, src) : 0.0;
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          const int src = i + nshift;
          h_x3f(m, k, j, i) = (src <= ie + ng) ? h_x3f(m, k, j, src) : 0.0;
        }
      }
    }
  }

  Kokkos::deep_copy(pmhd->u0, h_u0);
  Kokkos::deep_copy(pmhd->b0.x1f, h_x1f);
  Kokkos::deep_copy(pmhd->b0.x2f, h_x2f);
  Kokkos::deep_copy(pmhd->b0.x3f, h_x3f);

  auto &u0 = pmhd->u0;
  auto &b0 = pmhd->b0;
  auto &w0 = pmhd->w0;
  auto &bcc0 = pmhd->bcc0;
  pmhd->peos->ConsToPrim(u0, b0, w0, bcc0, false, is, ie, js, je, ks, ke);

  if (ps_frame_apply_to_particles) {
    const Real xshift = static_cast<Real>(nshift)*ps_recenter_dx1;
    ApplyRecenteringShiftToParticles(pm, xshift);
    pm->CountParticles();
  }
}

inline Real Uniform01(std::mt19937_64 &rng) {
  static std::uniform_real_distribution<Real> uni(0.0, 1.0);
  return uni(rng);
}

inline Real ClampInsideDomain(const Real x, const Real xmin, const Real xmax) {
  const Real span = xmax - xmin;
  const Real eps = std::max(static_cast<Real>(1.0e-12)*span,
                            static_cast<Real>(1.0e-14));
  return std::min(std::max(x, xmin + eps), xmax - eps);
}

void SeedNextTag(particles::Particles *ppart) {
  if (ps_tag_seeded) return;

  int local_max = -1;
  int npart = ppart->nprtcl_thispack;
  if (npart > 0) {
    auto h_pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_idata);
    for (int p = 0; p < npart; ++p) {
      local_max = std::max(local_max, h_pi(PTAG, p));
    }
  }

#if MPI_PARALLEL_ENABLED
  int global_max = -1;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
  int global_max = local_max;
#endif

  const std::int64_t start = static_cast<std::int64_t>(global_max) + 1;
  const std::int64_t nranks = static_cast<std::int64_t>(global_variable::nranks);
  std::int64_t rem = start % nranks;
  if (rem < 0) rem += nranks;
  std::int64_t shift = static_cast<std::int64_t>(global_variable::my_rank) - rem;
  if (shift < 0) shift += nranks;
  ps_next_tag = start + shift;
  ps_tag_seeded = true;
}

void ParallelShockSource(Mesh *pm, const Real bdt) {
  if (!ps_enable_injection) return;
  if (bdt <= 0.0) return;
  if (pm->time < ps_inject_t_start || pm->time > ps_inject_t_stop) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->pmhd == nullptr || pmbp->ppart == nullptr) return;

  auto *pmhd = pmbp->pmhd;
  auto *ppart = pmbp->ppart;
  SeedNextTag(ppart);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const bool three_d = pm->three_d;

  auto &mb_size = pmbp->pmb->mb_size;
  auto &mb_gid = pmbp->pmb->mb_gid;
  mb_size.template sync<HostMemSpace>();
  mb_gid.template sync<HostMemSpace>();

  const Real xshock = ShockSurfaceModelX1(pm->time);
  std::vector<ShockCell> cells;
  cells.reserve(static_cast<std::size_t>(pmbp->nmb_thispack)*indcs.nx2*indcs.nx3);

  Real area_total = 0.0;
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    const Real x1min = mb_size.h_view(m).x1min;
    const Real x2min = mb_size.h_view(m).x2min;
    const Real x3min = mb_size.h_view(m).x3min;
    const Real dx1 = mb_size.h_view(m).dx1;
    const Real dx2 = mb_size.h_view(m).dx2;
    const Real dx3 = mb_size.h_view(m).dx3;
    const Real half_width = ps_inject_half_width_cells*dx1;

    for (int k = ks; k <= ke; ++k) {
      Real x3c = x3min + (static_cast<Real>(k - ks) + 0.5)*dx3;
      if (!three_d) x3c = 0.0;
      for (int j = js; j <= je; ++j) {
        const Real x2c = x2min + (static_cast<Real>(j - js) + 0.5)*dx2;
        for (int i = is; i <= ie; ++i) {
          const Real x1c = x1min + (static_cast<Real>(i - is) + 0.5)*dx1;
          if (std::abs(x1c - xshock) > half_width) continue;

          ShockCell cell;
          cell.m = m;
          cell.k = k;
          cell.j = j;
          cell.i = i;
          cell.gid = mb_gid.h_view(m);
          cell.x1c = x1c;
          cell.x2c = x2c;
          cell.x3c = x3c;
          cell.dx1 = dx1;
          cell.dx2 = dx2;
          cell.dx3 = dx3;
          cell.area = dx2*dx3;
          cell.vol = dx1*dx2*dx3;
          area_total += cell.area;
          cells.push_back(cell);
        }
      }
    }
  }

  int ninj = 0;
  if (!cells.empty() && area_total > 0.0) {
    const Real swept_mass = ps_eta*ps_rho0*ps_shock_speed*bdt*area_total;
    const Real mass_budget = ps_mass_reservoir_local + swept_mass;
    if (mass_budget > 0.0) {
      ninj = static_cast<int>(std::floor(mass_budget/ps_particle_macro_mass));
      ps_mass_reservoir_local = mass_budget -
          static_cast<Real>(ninj)*ps_particle_macro_mass;
    }
  }

  int ninj_global = ninj;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&ninj, &ninj_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (ninj_global <= 0) return;

  if (ninj <= 0) {
    pm->CountParticles();
    return;
  }

  std::vector<Real> area_prefix(cells.size(), 0.0);
  Real running = 0.0;
  for (std::size_t n = 0; n < cells.size(); ++n) {
    running += cells[n].area;
    area_prefix[n] = running;
  }

  std::map<std::tuple<int, int, int, int>, GasDelta> gas_deltas;
  std::vector<InjectedParticle> injected;
  injected.reserve(static_cast<std::size_t>(ninj));

  const std::uint64_t seed =
      static_cast<std::uint64_t>(ps_inject_seed) ^
      (static_cast<std::uint64_t>(global_variable::my_rank + 1) * 1315423911ull) ^
      (static_cast<std::uint64_t>(pm->ncycle + 1) * 2654435761ull);
  std::mt19937_64 rng(seed);

  const Real frame_vx = FrameVelocityOffset(pm->time);
  const Real vinj = ps_vinj_over_u0*ps_u0;
  const auto &mesh_size = pm->mesh_size;
  const Real x1min = mesh_size.x1min;
  const Real x1max = mesh_size.x1max;
  const Real x2min = mesh_size.x2min;
  const Real x2max = mesh_size.x2max;
  const Real x3min = mesh_size.x3min;
  const Real x3max = mesh_size.x3max;
  for (int n = 0; n < ninj; ++n) {
    const Real draw = Uniform01(rng)*running;
    auto it = std::lower_bound(area_prefix.begin(), area_prefix.end(), draw);
    std::size_t idx = static_cast<std::size_t>(std::distance(area_prefix.begin(), it));
    if (idx >= cells.size()) idx = cells.size() - 1;
    const ShockCell &cell = cells[idx];

    Real dirx = 0.0;
    Real diry = 0.0;
    Real dirz = 0.0;
    if (three_d || ps_use_2d3v) {
      // Half-space isotropy toward +x keeps injected particles away from the
      // reflecting wall in this minimal benchmark implementation.
      const Real mu = Uniform01(rng);
      const Real phi = 2.0*M_PI*Uniform01(rng);
      const Real st = std::sqrt(std::max(static_cast<Real>(0.0), 1.0 - mu*mu));
      dirx = mu;
      diry = st*std::cos(phi);
      dirz = st*std::sin(phi);
    } else {
      const Real phi = M_PI*(Uniform01(rng) - 0.5);
      dirx = std::cos(phi);
      diry = std::sin(phi);
      dirz = 0.0;
    }

    InjectedParticle part;
    part.gid = cell.gid;
    part.m = cell.m;
    part.k = cell.k;
    part.j = cell.j;
    part.i = cell.i;
    part.vol = cell.vol;
    part.x1 = cell.x1c + (Uniform01(rng) - 0.5)*cell.dx1;
    part.x2 = cell.x2c + (Uniform01(rng) - 0.5)*cell.dx2;
    part.x3 = three_d ? (cell.x3c + (Uniform01(rng) - 0.5)*cell.dx3) : 0.0;
    part.x1 = ClampInsideDomain(part.x1, x1min, x1max);
    part.x2 = ClampInsideDomain(part.x2, x2min, x2max);
    if (three_d) part.x3 = ClampInsideDomain(part.x3, x3min, x3max);
    part.vx = ps_shock_speed + frame_vx + vinj*dirx;
    part.vy = vinj*diry;
    part.vz = vinj*dirz;
    injected.push_back(part);

    if (!ps_enable_subtraction) continue;
    const auto key = std::make_tuple(cell.m, cell.k, cell.j, cell.i);
    auto dit = gas_deltas.find(key);
    if (dit == gas_deltas.end()) {
      GasDelta delta;
      delta.m = cell.m;
      delta.k = cell.k;
      delta.j = cell.j;
      delta.i = cell.i;
      delta.dm = 0.0;
      delta.dmx = 0.0;
      delta.dmy = 0.0;
      delta.dmz = 0.0;
      delta.de = 0.0;
      dit = gas_deltas.emplace(key, delta).first;
    }
    const Real inv_vol = 1.0/part.vol;
    const Real mass_rho = ps_particle_macro_mass*inv_vol;
    dit->second.dm += mass_rho;
    dit->second.dmx += mass_rho*part.vx;
    dit->second.dmy += mass_rho*part.vy;
    dit->second.dmz += mass_rho*part.vz;
    dit->second.de += 0.5*mass_rho*(SQR(part.vx) + SQR(part.vy) + SQR(part.vz));
  }

  const std::int64_t nranks_ll = static_cast<std::int64_t>(global_variable::nranks);
  const std::int64_t nreq = nranks_ll*static_cast<std::int64_t>(injected.size() + 1);
  if (ps_next_tag >
      static_cast<std::int64_t>(std::numeric_limits<int>::max()) - nreq) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle tag range exhausted during pic_parallel_shock injection."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const int old_npart = ppart->nprtcl_thispack;
  const int ninject = static_cast<int>(injected.size());
  const int new_npart = old_npart + ninject;
  if (ninject <= 0) return;

  HostArray2D<int> h_pi_new("ps_pi_new", ppart->nidata, ninject);
  HostArray2D<Real> h_pr_new("ps_pr_new", ppart->nrdata, ninject);
  for (int n = 0; n < ninject; ++n) {
    for (int q = 0; q < ppart->nidata; ++q) h_pi_new(q, n) = 0;
    for (int q = 0; q < ppart->nrdata; ++q) h_pr_new(q, n) = 0.0;
  }

  for (int n = 0; n < ninject; ++n) {
    h_pi_new(PGID, n) = injected[n].gid;
    h_pi_new(PSP, n) = ps_inject_species;
    h_pi_new(PTAG, n) = static_cast<int>(ps_next_tag);
    ps_next_tag += nranks_ll;

    h_pr_new(IPX, n) = injected[n].x1;
    h_pr_new(IPY, n) = injected[n].x2;
    h_pr_new(IPZ, n) = injected[n].x3;
    h_pr_new(IPVX, n) = injected[n].vx;
    h_pr_new(IPVY, n) = injected[n].vy;
    h_pr_new(IPVZ, n) = injected[n].vz;
    h_pr_new(IPM, n) = ps_particle_q_over_m;
  }

  Kokkos::resize(ppart->prtcl_idata, ppart->nidata, new_npart);
  Kokkos::resize(ppart->prtcl_rdata, ppart->nrdata, new_npart);
  auto idata_new = Kokkos::subview(ppart->prtcl_idata, Kokkos::ALL,
                                   std::make_pair(old_npart, new_npart));
  auto rdata_new = Kokkos::subview(ppart->prtcl_rdata, Kokkos::ALL,
                                   std::make_pair(old_npart, new_npart));
  Kokkos::deep_copy(idata_new, h_pi_new);
  Kokkos::deep_copy(rdata_new, h_pr_new);
  ppart->nprtcl_thispack = new_npart;
  pm->CountParticles();

  if (!ps_enable_subtraction || gas_deltas.empty()) return;

  const int nsub = static_cast<int>(gas_deltas.size());
  HostArray1D<int> h_m("ps_m", nsub);
  HostArray1D<int> h_k("ps_k", nsub);
  HostArray1D<int> h_j("ps_j", nsub);
  HostArray1D<int> h_i("ps_i", nsub);
  HostArray1D<Real> h_dm("ps_dm", nsub);
  HostArray1D<Real> h_dmx("ps_dmx", nsub);
  HostArray1D<Real> h_dmy("ps_dmy", nsub);
  HostArray1D<Real> h_dmz("ps_dmz", nsub);
  HostArray1D<Real> h_de("ps_de", nsub);
  int idx_sub = 0;
  for (const auto &kv : gas_deltas) {
    const GasDelta &d = kv.second;
    h_m(idx_sub) = d.m;
    h_k(idx_sub) = d.k;
    h_j(idx_sub) = d.j;
    h_i(idx_sub) = d.i;
    h_dm(idx_sub) = d.dm;
    h_dmx(idx_sub) = d.dmx;
    h_dmy(idx_sub) = d.dmy;
    h_dmz(idx_sub) = d.dmz;
    h_de(idx_sub) = d.de;
    ++idx_sub;
  }

  auto d_m = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_m);
  auto d_k = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_k);
  auto d_j = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_j);
  auto d_i = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_i);
  auto d_dm = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dm);
  auto d_dmx = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmx);
  auto d_dmy = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmy);
  auto d_dmz = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_dmz);
  auto d_de = Kokkos::create_mirror_view_and_copy(DevExeSpace(), h_de);

  auto &u0 = pmhd->u0;
  auto &b0 = pmhd->b0;
  const Real rho_floor = ps_rho_floor_frac*ps_rho0;
  const Real p_floor = ps_p_floor_frac*ps_p0;
  const Real gm1 = pmhd->peos->eos_data.gamma - 1.0;
  par_for("ps_gas_subtract", DevExeSpace(), 0, nsub - 1,
  KOKKOS_LAMBDA(const int n) {
    const int m = d_m(n);
    const int k = d_k(n);
    const int j = d_j(n);
    const int i = d_i(n);
    const Real dm = d_dm(n);
    if (dm <= 0.0) return;

    const Real rho_old = u0(m, IDN, k, j, i);
    Real dm_eff = dm;
    if (rho_old - dm_eff < rho_floor) {
      dm_eff = fmax(rho_old - rho_floor, static_cast<Real>(0.0));
    }
    if (dm_eff <= 0.0) return;

    const Real frac = dm_eff/dm;
    u0(m, IDN, k, j, i) = rho_old - dm_eff;
    u0(m, IM1, k, j, i) -= frac*d_dmx(n);
    u0(m, IM2, k, j, i) -= frac*d_dmy(n);
    u0(m, IM3, k, j, i) -= frac*d_dmz(n);
    u0(m, IEN, k, j, i) -= frac*d_de(n);

    Real bx = 0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    Real by = 0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    Real bz = 0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
    Real emag = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    Real kin = 0.5*(SQR(u0(m, IM1, k, j, i)) +
                    SQR(u0(m, IM2, k, j, i)) +
                    SQR(u0(m, IM3, k, j, i)))/fmax(u0(m, IDN, k, j, i), rho_floor);
    Real efloor = p_floor/gm1 + kin + emag;
    if (u0(m, IEN, k, j, i) < efloor) {
      u0(m, IEN, k, j, i) = efloor;
    }
  });
}

void ParallelShockRefinement(MeshBlockPack *pmbp) {
  if (!ps_enable_curvature_amr) return;
  if (pmbp == nullptr || pmbp->pmhd == nullptr) return;
  if (!pmbp->pmesh->multilevel) return;
  if (!pmbp->pmesh->multi_d) return;

  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  auto &w0 = pmbp->pmhd->w0;
  auto &refine_flag = pmbp->pmesh->pmr->refine_flag;
  const int mbs = pmbp->pmesh->gids_eachrank[global_variable::my_rank];
  const Real refine_thresh = ps_refine_curv;
  const Real deref_thresh = ps_derefine_curv;
  const Real tiny = 1.0e-30;

  par_for("ps_curvature_amr", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
  KOKKOS_LAMBDA(const int m) {
    Real max_grho = 0.0;
    Real max_gprs = 0.0;
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real rho = fmax(w0(m, IDN, k, j, i), tiny);
          Real prs = fmax(w0(m, IEN, k, j, i), tiny);

          Real gx_rho = fabs((w0(m, IDN, k, j, i + 1) +
                              w0(m, IDN, k, j, i - 1))/rho - 2.0);
          Real gy_rho = fabs((w0(m, IDN, k, j + 1, i) +
                              w0(m, IDN, k, j - 1, i))/rho - 2.0);
          Real gx_prs = fabs((w0(m, IEN, k, j, i + 1) +
                              w0(m, IEN, k, j, i - 1))/prs - 2.0);
          Real gy_prs = fabs((w0(m, IEN, k, j + 1, i) +
                              w0(m, IEN, k, j - 1, i))/prs - 2.0);

          Real grho = gx_rho + gy_rho;
          Real gprs = gx_prs + gy_prs;
          max_grho = fmax(max_grho, grho);
          max_gprs = fmax(max_gprs, gprs);
        }
      }
    }

    if (fmax(max_grho, max_gprs) > refine_thresh) {
      refine_flag.d_view(m + mbs) = 1;
    } else if (max_grho < deref_thresh && max_gprs < deref_thresh) {
      refine_flag.d_view(m + mbs) = -1;
    }
  });
}

void MaybePrintFeedbackDiagnostics(Mesh *pm) {
  if (pm == nullptr) return;
  if (ps_feedback_diag_dcycle < 1) return;
  if ((pm->ncycle % ps_feedback_diag_dcycle) != 0) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp == nullptr || pmbp->ppart == nullptr) return;
  auto *ppart = pmbp->ppart;
  if (!ppart->deposit_moments) return;

  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  constexpr int nfield = 8;
  const int idx[nfield] = {
      particles::Particles::IMOM_JX,
      particles::Particles::IMOM_JY,
      particles::Particles::IMOM_JZ,
      particles::Particles::IMOM_DPXDT,
      particles::Particles::IMOM_DPYDT,
      particles::Particles::IMOM_DPZDT,
      particles::Particles::IMOM_DEDT,
      particles::Particles::IMOM_EBDOT,
  };

  Real sum_abs_local[nfield] = {0.0};
  Real sum_sq_local[nfield] = {0.0};
  Real max_abs_local[nfield] = {0.0};
  std::int64_t ncell_local = 0;
  constexpr int nprt_field = 8;
  const int pr_idx[nprt_field] = {
      IPDPX, IPDPY, IPDPZ, IPDE, IPEBDOT, IPBX, IPBY, IPBZ};
  Real pr_sum_abs_local[nprt_field] = {0.0};
  Real pr_sum_sq_local[nprt_field] = {0.0};
  Real pr_max_abs_local[nprt_field] = {0.0};

  auto h_mom = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->moments);
  for (int m = 0; m < pmbp->nmb_thispack; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          for (int n = 0; n < nfield; ++n) {
            const Real v = h_mom(m, idx[n], k, j, i);
            const Real av = std::abs(v);
            sum_abs_local[n] += av;
            sum_sq_local[n] += v*v;
            if (av > max_abs_local[n]) max_abs_local[n] = av;
          }
          ncell_local += 1;
        }
      }
    }
  }
  auto h_pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppart->prtcl_rdata);
  for (int p = 0; p < ppart->nprtcl_thispack; ++p) {
    for (int n = 0; n < nprt_field; ++n) {
      const Real v = h_pr(pr_idx[n], p);
      const Real av = std::abs(v);
      pr_sum_abs_local[n] += av;
      pr_sum_sq_local[n] += v*v;
      if (av > pr_max_abs_local[n]) pr_max_abs_local[n] = av;
    }
  }

  std::int64_t npart_local = static_cast<std::int64_t>(ppart->nprtcl_thispack);
  std::int64_t npart_global = npart_local;
  std::int64_t ncell_global = ncell_local;
  Real sum_abs_global[nfield] = {0.0};
  Real sum_sq_global[nfield] = {0.0};
  Real max_abs_global[nfield] = {0.0};
  Real pr_sum_abs_global[nprt_field] = {0.0};
  Real pr_sum_sq_global[nprt_field] = {0.0};
  Real pr_max_abs_global[nprt_field] = {0.0};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(sum_abs_local, sum_abs_global, nfield, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(sum_sq_local, sum_sq_global, nfield, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(max_abs_local, max_abs_global, nfield, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(pr_sum_abs_local, pr_sum_abs_global, nprt_field, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pr_sum_sq_local, pr_sum_sq_global, nprt_field, MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pr_max_abs_local, pr_max_abs_global, nprt_field, MPI_ATHENA_REAL,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&ncell_local, &ncell_global, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&npart_local, &npart_global, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
#else
  for (int n = 0; n < nfield; ++n) {
    sum_abs_global[n] = sum_abs_local[n];
    sum_sq_global[n] = sum_sq_local[n];
    max_abs_global[n] = max_abs_local[n];
  }
  for (int n = 0; n < nprt_field; ++n) {
    pr_sum_abs_global[n] = pr_sum_abs_local[n];
    pr_sum_sq_global[n] = pr_sum_sq_local[n];
    pr_max_abs_global[n] = pr_max_abs_local[n];
  }
#endif

  if (global_variable::my_rank != 0) return;
  const Real inv_ncell = (ncell_global > 0) ?
      1.0/static_cast<Real>(ncell_global) : 0.0;
  Real rms[nfield];
  for (int n = 0; n < nfield; ++n) {
    rms[n] = std::sqrt(std::max(sum_sq_global[n]*inv_ncell, 0.0));
  }
  const Real inv_npart = (npart_global > 0) ?
      1.0/static_cast<Real>(npart_global) : 0.0;
  Real pr_rms[nprt_field];
  for (int n = 0; n < nprt_field; ++n) {
    pr_rms[n] = std::sqrt(std::max(pr_sum_sq_global[n]*inv_npart, 0.0));
  }

  static bool printed_header = false;
  if (!printed_header) {
    const bool is_boris =
        (ppart->pusher == ParticlesPusher::boris_lin) ||
        (ppart->pusher == ParticlesPusher::boris_tsc);
    const bool use_delta_feedback =
        (is_boris &&
         (ppart->pic_feedback_mode == PICFeedbackMode::coupled) &&
         ppart->couple_moments_to_mhd &&
         (ppart->couple_moments_momentum_to_mhd || ppart->couple_moments_energy_to_mhd));
    std::cout << "pic_parallel_shock feedback_diag_cfg: "
              << "pusher=" << static_cast<int>(ppart->pusher)
              << " pic_feedback_mode=" << static_cast<int>(ppart->pic_feedback_mode)
              << " couple_moments_to_mhd=" << (ppart->couple_moments_to_mhd ? 1 : 0)
              << " mom_feedback=" << (ppart->couple_moments_momentum_to_mhd ? 1 : 0)
              << " eng_feedback=" << (ppart->couple_moments_energy_to_mhd ? 1 : 0)
              << " pic_enable_2d3v=" << (ppart->pic_enable_2d3v ? 1 : 0)
              << " use_delta_feedback=" << (use_delta_feedback ? 1 : 0)
              << std::endl;
    printed_header = true;
  }

  const Real step_mom_l1 = pm->dt*ppart->couple_moments_momentum_coeff*
      (sum_abs_global[3] + sum_abs_global[4] + sum_abs_global[5]);
  const Real step_eng_l1 = pm->dt*ppart->couple_moments_energy_coeff*
      sum_abs_global[6];
  auto old_flags = std::cout.flags();
  auto old_prec = std::cout.precision();
  std::cout << std::scientific << std::setprecision(12);
  std::cout << "pic_parallel_shock feedback_diag: cycle=" << pm->ncycle
            << " time=" << pm->time
            << " npart=" << npart_global
            << " ncell=" << ncell_global
            << " j_rms=(" << rms[0] << "," << rms[1] << "," << rms[2] << ")"
            << " dpdt_rms=(" << rms[3] << "," << rms[4] << "," << rms[5] << ")"
            << " dedt_rms=" << rms[6]
            << " ebdot_rms=" << rms[7]
            << " dpdt_l1=(" << sum_abs_global[3] << ","
            << sum_abs_global[4] << "," << sum_abs_global[5] << ")"
            << " dedt_l1=" << sum_abs_global[6]
            << " step_mom_l1=" << step_mom_l1
            << " step_eng_l1=" << step_eng_l1
            << " max_abs_dpdt=(" << max_abs_global[3] << ","
            << max_abs_global[4] << "," << max_abs_global[5] << ")"
            << " pr_dpdt_rms=(" << pr_rms[0] << "," << pr_rms[1] << ","
            << pr_rms[2] << ")"
            << " pr_dedt_rms=" << pr_rms[3]
            << " pr_ebdot_rms=" << pr_rms[4]
            << " pr_b_rms=(" << pr_rms[5] << "," << pr_rms[6] << ","
            << pr_rms[7] << ")"
            << " pr_dpdt_l1=(" << pr_sum_abs_global[0] << ","
            << pr_sum_abs_global[1] << "," << pr_sum_abs_global[2] << ")"
            << " pr_max_abs_dpdt=(" << pr_max_abs_global[0] << ","
            << pr_max_abs_global[1] << "," << pr_max_abs_global[2] << ")"
            << std::endl;
  std::cout.flags(old_flags);
  std::cout.precision(old_prec);
}

void ParallelShockWorkInLoop(Mesh *pm) {
  if (pm == nullptr || pm->dt <= 0.0) return;
  if (ps_frame_diag_dcycle < 1) ps_frame_diag_dcycle = 1;

  if (!ps_enable_frame_tracking) {
    MaybePrintFeedbackDiagnostics(pm);
    return;
  }

  if (FrameModeVelocity()) {
    const Real v_now = FrameVelocityOffset(pm->time);
    const Real v_next_target = FrameVelocityOffset(pm->time + pm->dt);
    Real dv = v_next_target - v_now;
    const bool clipped = (ps_frame_dv_max > 0.0 && std::abs(dv) > ps_frame_dv_max);
    if (ps_frame_dv_max > 0.0 && std::abs(dv) > ps_frame_dv_max) {
      dv = std::copysign(ps_frame_dv_max, dv);
    }
    if (std::abs(dv) <= 1.0e-15) return;

    ApplyFrameShiftToFluid(pm, dv);
    if (ps_frame_apply_to_particles) {
      ApplyFrameShiftToParticles(pm, dv);
    }
    if (ps_frame_apply_to_inflow) {
      UpdateOuterInflowState(pm, v_next_target);
    }

    if (global_variable::my_rank == 0 &&
        ((pm->ncycle % ps_frame_diag_dcycle) == 0)) {
      const Real x_model = ShockSurfaceModelX1(pm->time + pm->dt);
      std::cout << "pic_parallel_shock frame(velocity): cycle=" << pm->ncycle
                << " time=" << pm->time
                << " v_now=" << v_now
                << " v_target=" << v_next_target
                << " dv=" << dv
                << " clipped=" << (clipped ? 1 : 0)
                << " xshock_model=" << x_model << std::endl;
    }
    MaybePrintFeedbackDiagnostics(pm);
    return;
  }

  if (FrameModeRecenter()) {
    const int ev_now = RecenterEventCount(pm->time);
    const int ev_next = RecenterEventCount(pm->time + pm->dt);
    const int dev = ev_next - ev_now;
    int nshift = 0;
    if (dev > 0) {
      for (int n = 0; n < dev; ++n) {
        ApplyRecenteringShift(pm, ps_recenter_shift_cells);
      }
      nshift = dev*ps_recenter_shift_cells;
    }

    if (global_variable::my_rank == 0 &&
        ((pm->ncycle % ps_frame_diag_dcycle) == 0 || dev > 0)) {
      const Real x_model = ShockSurfaceModelX1(pm->time + pm->dt);
      const Real xshift = static_cast<Real>(nshift)*ps_recenter_dx1;
      std::cout << "pic_parallel_shock frame(recenter): cycle=" << pm->ncycle
                << " time=" << pm->time
                << " ev_now=" << ev_now
                << " ev_next=" << ev_next
                << " dev=" << dev
                << " nshift=" << nshift
                << " xshift=" << xshift
                << " xshock_model=" << x_model << std::endl;
    }
  }
  MaybePrintFeedbackDiagnostics(pm);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::PICParallelShock()
//! \brief Parallel-shock MHD-PIC benchmark setup.

void ProblemGenerator::PICParallelShock(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires an active <mhd> block." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmbp->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires an active <particles> block."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmbp->phydro != nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock is MHD-only; do not enable <hydro>."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!pmbp->pmhd->peos->eos_data.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ideal MHD EOS." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::reflect) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects <mesh>/ix1_bc=reflect." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const BoundaryFlag ox1_bc = pmy_mesh_->mesh_bcs[BoundaryFace::outer_x1];
  if (!(ox1_bc == BoundaryFlag::inflow || ox1_bc == BoundaryFlag::outflow)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects <mesh>/ox1_bc=inflow or outflow."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
      pmy_mesh_->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock expects periodic y boundaries." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Parse runtime controls.
  ps_rho0 = pin->GetOrAddReal("problem", "ps_rho0", 1.0);
  ps_p0 = pin->GetOrAddReal("problem", "ps_p0", 1.0);
  ps_u0 = pin->GetOrAddReal("problem", "ps_u0", 30.0);
  ps_b0 = pin->GetOrAddReal("problem", "ps_b0", 1.0);
  ps_eta = pin->GetOrAddReal("problem", "ps_eta", 1.0e-3);
  ps_vinj_over_u0 = pin->GetOrAddReal("problem", "ps_vinj_over_u0", std::sqrt(10.0));
  ps_inject_half_width_cells = pin->GetOrAddReal(
      "problem", "ps_inject_half_width_cells", 0.5);
  ps_inject_t_start = pin->GetOrAddReal("problem", "ps_inject_t_start", 0.0);
  ps_inject_t_stop = pin->GetOrAddReal("problem", "ps_inject_t_stop", 1.0e99);
  ps_seed_noise_amp = pin->GetOrAddReal("problem", "ps_seed_noise_amp", 0.0);
  ps_seed_noise_seed = pin->GetOrAddInteger("problem", "ps_seed_noise_seed", 1234);
  ps_refine_curv = pin->GetOrAddReal("problem", "ps_refine_curv", 1.0);
  ps_derefine_curv = pin->GetOrAddReal("problem", "ps_derefine_curv", 0.1);
  ps_rho_floor_frac = pin->GetOrAddReal("problem", "ps_rho_floor_frac", 1.0e-6);
  ps_p_floor_frac = pin->GetOrAddReal("problem", "ps_p_floor_frac", 1.0e-8);
  ps_enable_injection = pin->GetOrAddBoolean("problem", "ps_enable_injection", true);
  ps_enable_subtraction = pin->GetOrAddBoolean(
      "problem", "ps_enable_gas_subtraction", true);
  ps_enable_curvature_amr = pin->GetOrAddBoolean(
      "problem", "ps_enable_curvature_amr", true);
  ps_enable_frame_tracking = pin->GetOrAddBoolean(
      "problem", "ps_enable_frame_tracking", false);
  std::string frame_mode = pin->GetOrAddString("problem", "ps_frame_mode",
                                                "velocity");
  ps_frame_t_start = pin->GetOrAddReal("problem", "ps_frame_t_start", 0.0);
  ps_frame_t_ramp = pin->GetOrAddReal("problem", "ps_frame_t_ramp", 0.0);
  ps_frame_vfrac = pin->GetOrAddReal("problem", "ps_frame_vfrac", 1.0);
  ps_frame_dv_max = pin->GetOrAddReal("problem", "ps_frame_dv_max", 1.0e99);
  ps_frame_apply_to_particles = pin->GetOrAddBoolean(
      "problem", "ps_frame_apply_to_particles", true);
  ps_frame_apply_to_inflow = pin->GetOrAddBoolean(
      "problem", "ps_frame_apply_to_inflow", true);
  ps_frame_require_uniform = pin->GetOrAddBoolean(
      "problem", "ps_frame_require_uniform", true);
  ps_frame_diag_dcycle = pin->GetOrAddInteger(
      "problem", "ps_frame_diag_dcycle", 200);
  ps_feedback_diag_dcycle = pin->GetOrAddInteger(
      "problem", "ps_feedback_diag_dcycle", 0);
  ps_recenter_x_target = pin->GetOrAddReal("problem", "ps_recenter_x_target", 2.0);
  ps_recenter_x_trigger = pin->GetOrAddReal("problem", "ps_recenter_x_trigger", 3.0);
  ps_recenter_shift_cells = pin->GetOrAddInteger("problem",
                                                  "ps_recenter_shift_cells", 2);
  ps_recenter_vshock_model = pin->GetOrAddReal("problem",
                                                "ps_recenter_vshock_model", -1.0);
  ps_inject_species = pin->GetOrAddInteger("problem", "ps_inject_species", 0);
  ps_inject_seed = pin->GetOrAddInteger("problem", "ps_inject_seed", 1234);

  if (ps_rho0 <= 0.0 || ps_p0 <= 0.0 || ps_u0 <= 0.0 || ps_b0 <= 0.0 || ps_eta < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock parameters ps_rho0/ps_p0/ps_u0/ps_b0 must be > 0 "
              << "and ps_eta must be >= 0." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ps_seed_noise_amp < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_seed_noise_amp >= 0."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ps_frame_t_ramp < 0.0 || ps_frame_vfrac < 0.0 || ps_frame_dv_max < 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock frame controls require ps_frame_t_ramp >= 0, "
              << "ps_frame_vfrac >= 0, and ps_frame_dv_max >= 0." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ps_frame_diag_dcycle < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_frame_diag_dcycle >= 1."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ps_feedback_diag_dcycle < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires ps_feedback_diag_dcycle >= 0."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (frame_mode == "velocity") {
    ps_frame_mode = PSFrameMode::velocity;
  } else if (frame_mode == "recenter") {
    ps_frame_mode = PSFrameMode::recenter;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "ps_frame_mode must be 'velocity' or 'recenter'." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const int nspecies = pmbp->ppart->nspecies;
  if (ps_inject_species < 0 || ps_inject_species >= nspecies) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "ps_inject_species is out of range for configured species."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string species_block = "species" + std::to_string(ps_inject_species);
  ps_particle_mass = pin->GetOrAddReal(species_block, "mass", 1.0);
  ps_particle_charge = pin->GetOrAddReal(species_block, "charge", 1.0);
  if (ps_particle_mass <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Injected species mass must be > 0." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  ps_particle_q_over_m = ps_particle_charge/ps_particle_mass;
  Real qscale = pin->GetOrAddReal("particles", "deposit_qscale", 1.0);
  ps_particle_macro_mass = qscale*ps_particle_mass;
  if (ps_particle_macro_mass <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Computed injected macro-mass must be > 0." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Analytic finite-Mach shock-speed estimate in the reflecting-wall frame.
  Real gamma = pmbp->pmhd->peos->eos_data.gamma;
  ps_shock_speed = EstimateShockSpeed(gamma, ps_rho0, ps_p0, ps_u0);
  if (ps_shock_speed <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock requires supersonic inflow with positive "
              << "finite-Mach shock-speed estimate." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  ps_xshock0 = pmy_mesh_->mesh_size.x1min;
  ps_recenter_dx1 = pmy_mesh_->mesh_size.dx1;
  ps_use_2d3v = (pmy_mesh_->two_d && pmbp->ppart->pic_enable_2d3v);
  ps_mass_reservoir_local = 0.0;
  ps_tag_seeded = false;
  ps_next_tag = 0;
  ConfigureSeedNoisePhases();

  if (ps_enable_frame_tracking && ps_frame_require_uniform &&
      (pmy_mesh_->adaptive || pmy_mesh_->multilevel)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "pic_parallel_shock frame tracking currently supports uniform grids "
              << "only (no AMR/SMR)." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ps_enable_frame_tracking && FrameModeRecenter()) {
    if (!(ps_recenter_x_target > pmy_mesh_->mesh_size.x1min &&
          ps_recenter_x_target < pmy_mesh_->mesh_size.x1max)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_target must lie strictly inside the x1 domain."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(ps_recenter_x_trigger > pmy_mesh_->mesh_size.x1min)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_trigger must be greater than x1min."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(ps_recenter_x_trigger < pmy_mesh_->mesh_size.x1max)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_x_trigger must be less than x1max."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ps_recenter_dx1 <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "pic_parallel_shock recentering requires positive mesh dx1."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ps_recenter_shift_cells < 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_shift_cells must be >= 1 in recenter mode."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ps_recenter_shift_cells > pmy_mesh_->mb_indcs.ng) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_shift_cells must be <= nghost in recenter mode."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const Real vmodel = ps_recenter_vshock_model;
    const bool use_default = std::abs(vmodel + 1.0) < 1.0e-12;
    const bool use_override = vmodel > 0.0;
    if (!(use_default || use_override)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "ps_recenter_vshock_model must be > 0 or left at default -1."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Enroll benchmark callbacks.
  user_srcs = true;
  user_srcs_func = ParallelShockSource;
  user_ref_func = ParallelShockRefinement;
  user_work_in_loop = true;
  user_work_in_loop_func = ParallelShockWorkInLoop;

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mb_size;
  Real gm1 = gamma - 1.0;
  const bool add_seed_noise = (ps_seed_noise_amp > 0.0);
  const Real seed_noise_amp = ps_seed_noise_amp*ps_b0;
  const Real x1min_global = pmy_mesh_->mesh_size.x1min;
  const Real lx_global = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const Real seed_noise_norm = 1.0/1.875;
  const Real ph_by_1 = ps_seed_noise_phase_by[0];
  const Real ph_by_2 = ps_seed_noise_phase_by[1];
  const Real ph_by_3 = ps_seed_noise_phase_by[2];
  const Real ph_by_4 = ps_seed_noise_phase_by[3];
  const Real ph_bz_1 = ps_seed_noise_phase_bz[0];
  const Real ph_bz_2 = ps_seed_noise_phase_bz[1];
  const Real ph_bz_3 = ps_seed_noise_phase_bz[2];
  const Real ph_bz_4 = ps_seed_noise_phase_bz[3];

  // Set uniform upstream state (flow toward the reflecting wall).
  par_for("pgen_pic_parallel_shock", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    u0(m, IDN, k, j, i) = ps_rho0;
    u0(m, IM1, k, j, i) = -ps_rho0*ps_u0;
    u0(m, IM2, k, j, i) = 0.0;
    u0(m, IM3, k, j, i) = 0.0;

    Real dby = 0.0;
    Real dbz = 0.0;
    if (add_seed_noise) {
      const Real x1c = size.d_view(m).x1min +
          (static_cast<Real>(i - is) + 0.5)*size.d_view(m).dx1;
      const Real xfrac = (x1c - x1min_global)/lx_global;
      const Real xphase = 2.0*M_PI*xfrac;

      const Real nby = sin(1.0*xphase + ph_by_1) +
                       0.5*sin(2.0*xphase + ph_by_2) +
                       0.25*sin(4.0*xphase + ph_by_3) +
                       0.125*sin(8.0*xphase + ph_by_4);
      const Real nbz = sin(1.0*xphase + ph_bz_1) +
                       0.5*sin(2.0*xphase + ph_bz_2) +
                       0.25*sin(4.0*xphase + ph_bz_3) +
                       0.125*sin(8.0*xphase + ph_bz_4);
      dby = seed_noise_amp*seed_noise_norm*nby;
      dbz = seed_noise_amp*seed_noise_norm*nbz;
    }

    b0.x1f(m, k, j, i) = ps_b0;
    b0.x2f(m, k, j, i) = dby;
    b0.x3f(m, k, j, i) = dbz;
    if (i == ie) b0.x1f(m, k, j, i + 1) = ps_b0;
    if (j == je) b0.x2f(m, k, j + 1, i) = dby;
    if (k == ke) b0.x3f(m, k + 1, j, i) = dbz;
  });

  par_for("pgen_pic_parallel_shock_e", DevExeSpace(), 0, pmbp->nmb_thispack - 1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real bx = 0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i + 1));
    Real by = 0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j + 1, i));
    Real bz = 0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k + 1, j, i));
    Real ekin = 0.5*(SQR(u0(m, IM1, k, j, i)) +
                     SQR(u0(m, IM2, k, j, i)) +
                     SQR(u0(m, IM3, k, j, i)))/u0(m, IDN, k, j, i);
    Real emag = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    u0(m, IEN, k, j, i) = ps_p0/gm1 + ekin + emag;
  });

  // Outer x1 inflow reservoir set to upstream state.
  UpdateOuterInflowState(pmy_mesh_, FrameVelocityOffset(pmy_mesh_->time));

  // Keep size data synced for host-side injection cell selection.
  size.template modify<HostMemSpace>();
  size.template sync<DevExeSpace>();
}
