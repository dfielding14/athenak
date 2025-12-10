//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file scalar_mixing.cpp
//  \brief Problem generator for scalar mixing in frozen turbulent velocity field
//
//  Single passive scalar with mean gradient forcing: ds/dt = G * v_x
//  This simulates mixing of a background scalar gradient by turbulence.
//
//  When scalar_only=true in <hydro>, velocity is frozen and only the scalar evolves.

#include <iostream>
#include <cmath>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen.hpp"
#include "globals.hpp"
#include "utils/random.hpp"

#include <Kokkos_Random.hpp>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

// Function declarations
void MeanGradientForcing(Mesh* pm, const Real bdt);

// Store mean gradient parameter (read in ProblemGenerator, used in source term)
namespace {
  Real mean_gradient_G = 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn InitTurbulentVelocity
//  \brief Initialize a turbulent velocity field with specified power-law spectrum.

void InitTurbulentVelocity(MeshBlockPack *pmbp, ParameterInput *pin, Real den) {
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmbp->nmb_thispack;

  // Read turbulence parameters
  Real v_rms = pin->GetOrAddReal("problem", "turb_v_rms", 1.0);
  int nlow = pin->GetOrAddInteger("problem", "turb_nlow", 1);
  int nhigh = pin->GetOrAddInteger("problem", "turb_nhigh", 4);
  Real expo = pin->GetOrAddReal("problem", "turb_expo", 5.0/3.0);
  Real sol_frac = pin->GetOrAddReal("problem", "turb_sol_frac", 1.0);
  int rseed = pin->GetOrAddInteger("problem", "turb_rseed", 12345);

  // Domain size for computing wavenumbers
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  Real dkx = 2.0*M_PI/lx;
  Real dky = (nx2 > 1) ? 2.0*M_PI/ly : 0.0;
  Real dkz = (nx3 > 1) ? 2.0*M_PI/lz : 0.0;

  // Determine dimensionality
  bool multi_d = (nx2 > 1);
  bool three_d = (nx3 > 1);

  // Count modes within [nlow, nhigh]
  int nlow_sqr = nlow*nlow;
  int nhigh_sqr = nhigh*nhigh;
  int mode_count = 0;

  for (int nkx = 0; nkx <= nhigh; nkx++) {
    for (int nky = (multi_d ? -nhigh : 0); nky <= (multi_d ? nhigh : 0); nky++) {
      for (int nkz = (three_d ? -nhigh : 0); nkz <= (three_d ? nhigh : 0); nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) {
          mode_count++;
        }
      }
    }
  }

  if (mode_count == 0) {
    std::cout << "### WARNING: No turbulent modes in range [" << nlow << ", " << nhigh
              << "]. Using zero velocity." << std::endl;
    return;
  }

  if (global_variable::my_rank == 0) {
    std::cout << "Initializing turbulent velocity field:" << std::endl
              << "  v_rms=" << v_rms << ", nlow=" << nlow << ", nhigh=" << nhigh
              << ", expo=" << expo << ", sol_frac=" << sol_frac << std::endl
              << "  mode_count=" << mode_count << std::endl;
  }

  // Allocate arrays for Fourier coefficients
  std::vector<Real> kx_arr(mode_count), ky_arr(mode_count), kz_arr(mode_count);
  std::vector<Real> aka0(mode_count), aka1(mode_count), aka2(mode_count);
  std::vector<Real> akb0(mode_count), akb1(mode_count), akb2(mode_count);

  // Initialize RNG
  RNG_State rstate;
  rstate.idum = -rseed;

  // Generate mode amplitudes and phases
  int nmode = 0;
  for (int nkx = 0; nkx <= nhigh; nkx++) {
    for (int nky = (multi_d ? -nhigh : 0); nky <= (multi_d ? nhigh : 0); nky++) {
      for (int nkz = (three_d ? -nhigh : 0); nkz <= (three_d ? nhigh : 0); nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) {
          Real kx = dkx*nkx;
          Real ky = dky*nky;
          Real kz = dkz*nkz;
          Real kiso = sqrt(kx*kx + ky*ky + kz*kz);

          kx_arr[nmode] = kx;
          ky_arr[nmode] = ky;
          kz_arr[nmode] = kz;

          // Amplitude scaling for E(k) ~ k^{-expo}
          Real norm = (kiso > 1e-16) ? 1.0/pow(kiso, (expo+2.0)/2.0) : 0.0;

          // Generate random Fourier coefficients
          Real k[3] = {kx, ky, kz};
          Real a[3], b[3];
          Real ka = 0.0, kb = 0.0;

          for (int dir = 0; dir < 3; dir++) {
            a[dir] = norm * RanGaussianSt(&rstate);
            b[dir] = norm * RanGaussianSt(&rstate);
            ka += k[dir] * b[dir];
            kb += k[dir] * a[dir];
          }

          // Apply solenoidal projection
          Real kiso_sqr = kiso*kiso;
          for (int dir = 0; dir < 3; dir++) {
            Real diva = k[dir]*ka/kiso_sqr;
            Real divb = k[dir]*kb/kiso_sqr;
            Real curla = a[dir] - divb;
            Real curlb = b[dir] - diva;
            a[dir] = sol_frac*curla + (1.0-sol_frac)*divb;
            b[dir] = sol_frac*curlb + (1.0-sol_frac)*diva;
          }

          aka0[nmode] = a[0]; aka1[nmode] = a[1]; aka2[nmode] = a[2];
          akb0[nmode] = b[0]; akb1[nmode] = b[1]; akb2[nmode] = b[2];
          nmode++;
        }
      }
    }
  }

  // Copy coefficients to device
  DvceArray1D<Real> d_kx("kx", mode_count), d_ky("ky", mode_count), d_kz("kz", mode_count);
  DvceArray1D<Real> d_aka0("aka0", mode_count), d_aka1("aka1", mode_count), d_aka2("aka2", mode_count);
  DvceArray1D<Real> d_akb0("akb0", mode_count), d_akb1("akb1", mode_count), d_akb2("akb2", mode_count);

  auto h_kx = Kokkos::create_mirror_view(d_kx);
  auto h_ky = Kokkos::create_mirror_view(d_ky);
  auto h_kz = Kokkos::create_mirror_view(d_kz);
  auto h_aka0 = Kokkos::create_mirror_view(d_aka0);
  auto h_aka1 = Kokkos::create_mirror_view(d_aka1);
  auto h_aka2 = Kokkos::create_mirror_view(d_aka2);
  auto h_akb0 = Kokkos::create_mirror_view(d_akb0);
  auto h_akb1 = Kokkos::create_mirror_view(d_akb1);
  auto h_akb2 = Kokkos::create_mirror_view(d_akb2);

  for (int n = 0; n < mode_count; n++) {
    h_kx(n) = kx_arr[n]; h_ky(n) = ky_arr[n]; h_kz(n) = kz_arr[n];
    h_aka0(n) = aka0[n]; h_aka1(n) = aka1[n]; h_aka2(n) = aka2[n];
    h_akb0(n) = akb0[n]; h_akb1(n) = akb1[n]; h_akb2(n) = akb2[n];
  }

  Kokkos::deep_copy(d_kx, h_kx);
  Kokkos::deep_copy(d_ky, h_ky);
  Kokkos::deep_copy(d_kz, h_kz);
  Kokkos::deep_copy(d_aka0, h_aka0);
  Kokkos::deep_copy(d_aka1, h_aka1);
  Kokkos::deep_copy(d_aka2, h_aka2);
  Kokkos::deep_copy(d_akb0, h_akb0);
  Kokkos::deep_copy(d_akb1, h_akb1);
  Kokkos::deep_copy(d_akb2, h_akb2);

  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  int mode_count_ = mode_count;

  // Compute velocity field by summing Fourier modes
  par_for("turb_vel_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real vx = 0.0, vy = 0.0, vz = 0.0;
    for (int n = 0; n < mode_count_; n++) {
      Real kdotx = d_kx(n)*x1v + d_ky(n)*x2v + d_kz(n)*x3v;
      Real cosk = cos(kdotx);
      Real sink = sin(kdotx);

      vx += d_aka0(n)*cosk - d_akb0(n)*sink;
      vy += d_aka1(n)*cosk - d_akb1(n)*sink;
      vz += d_aka2(n)*cosk - d_akb2(n)*sink;
    }

    u0(m,IM1,k,j,i) = den * vx;
    u0(m,IM2,k,j,i) = den * vy;
    u0(m,IM3,k,j,i) = den * vz;
  });

  // Compute mean momentum and v_rms for normalization
  const int nmkji = nmb*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  Real sum_mass = 0.0, sum_mom1 = 0.0, sum_mom2 = 0.0, sum_mom3 = 0.0, sum_vsqr = 0.0;
  Kokkos::parallel_reduce("turb_vel_stats", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &s_mass, Real &s_mom1, Real &s_mom2, Real &s_mom3, Real &s_vsqr) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    Real rho = u0(m,IDN,k,j,i);
    Real mom1 = u0(m,IM1,k,j,i);
    Real mom2 = u0(m,IM2,k,j,i);
    Real mom3 = u0(m,IM3,k,j,i);

    s_mass += rho * vol;
    s_mom1 += mom1 * vol;
    s_mom2 += mom2 * vol;
    s_mom3 += mom3 * vol;
    s_vsqr += (mom1*mom1 + mom2*mom2 + mom3*mom3)/rho * vol;
  }, Kokkos::Sum<Real>(sum_mass), Kokkos::Sum<Real>(sum_mom1),
     Kokkos::Sum<Real>(sum_mom2), Kokkos::Sum<Real>(sum_mom3),
     Kokkos::Sum<Real>(sum_vsqr));

#if MPI_PARALLEL_ENABLED
  Real local_sums[5] = {sum_mass, sum_mom1, sum_mom2, sum_mom3, sum_vsqr};
  Real global_sums[5];
  MPI_Allreduce(local_sums, global_sums, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum_mass = global_sums[0];
  sum_mom1 = global_sums[1];
  sum_mom2 = global_sums[2];
  sum_mom3 = global_sums[3];
  sum_vsqr = global_sums[4];
#endif

  Real vmean1 = sum_mom1 / sum_mass;
  Real vmean2 = sum_mom2 / sum_mass;
  Real vmean3 = sum_mom3 / sum_mass;

  Real totvol = sum_mass / den;
  Real current_vsqr = sum_vsqr / totvol - (vmean1*vmean1 + vmean2*vmean2 + vmean3*vmean3);
  Real current_vrms = sqrt(current_vsqr);

  Real scale = (current_vrms > 1e-16) ? v_rms / current_vrms : 0.0;

  if (global_variable::my_rank == 0) {
    std::cout << "  vmean = (" << vmean1 << ", " << vmean2 << ", " << vmean3 << ")" << std::endl
              << "  current_vrms = " << current_vrms << ", scale = " << scale << std::endl;
  }

  // Apply mean subtraction and v_rms scaling
  par_for("turb_vel_normalize", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real rho = u0(m,IDN,k,j,i);
    u0(m,IM1,k,j,i) = rho * scale * (u0(m,IM1,k,j,i)/rho - vmean1);
    u0(m,IM2,k,j,i) = rho * scale * (u0(m,IM2,k,j,i)/rho - vmean2);
    u0(m,IM3,k,j,i) = rho * scale * (u0(m,IM3,k,j,i)/rho - vmean3);
  });
}

//----------------------------------------------------------------------------------------
//! \brief Problem Generator for scalar mixing
//
//  Initializes density, turbulent velocity, and a single passive scalar.
//  Scalar is initialized to scalar_init (default 0.5) to allow symmetric fluctuations.
//  Mean gradient forcing ds/dt = G*v_x is applied via UserSourceTerm.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  int nscalars = pmbp->phydro->nscalars;
  int nhydro = pmbp->phydro->nhydro;

  // Read mean gradient parameter
  mean_gradient_G = pin->GetOrAddReal("problem", "mean_gradient", 1.0);

  // Enroll user source term
  user_srcs_func = MeanGradientForcing;

  if (restart) return;

  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;

  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  // Read initial conditions
  Real den = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real temp = pin->GetOrAddReal("problem", "temp0", 1.0);
  Real scalar_init = pin->GetOrAddReal("problem", "scalar_init", 0.5);

  int nmb = pmbp->nmb_thispack;

  if (global_variable::my_rank == 0) {
    std::cout << "Scalar mixing problem:" << std::endl
              << "  mean_gradient G = " << mean_gradient_G << std::endl
              << "  scalar_init = " << scalar_init << std::endl
              << "  nscalars = " << nscalars << std::endl;
  }

  // Initialize density and energy
  par_for("pgen_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IDN,k,j,i) = den;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = den * temp / gm1;
  });

  // Initialize turbulent velocity field
  InitTurbulentVelocity(pmbp, pin, den);

  // Update energy with kinetic energy
  par_for("pgen_energy", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real ke = 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
                   SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
    u0(m,IEN,k,j,i) += ke;
  });

  // Initialize passive scalar to offset value (allows symmetric fluctuations)
  par_for("pgen_scalar", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real rho = u0(m, IDN, k, j, i);
    for (int ns = 0; ns < nscalars; ns++) {
      u0(m, nhydro+ns, k, j, i) = rho * scalar_init;
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeanGradientForcing
//  \brief Apply mean gradient forcing to passive scalar: ds/dt = G * v_x
//
//  This simulates the effect of a background mean scalar gradient in the x-direction.
//  Turbulence mixes this gradient, creating scalar fluctuations correlated with velocity.

void MeanGradientForcing(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  int nhydro = pmbp->phydro->nhydro;

  Real G = mean_gradient_G;

  par_for("mean_grad_forcing", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real density = w0(m, IDN, k, j, i);
    Real vx = w0(m, IVX, k, j, i);

    // Mean gradient forcing: ds/dt = G * v_x
    // In conservative form: d(rho*s)/dt = rho * G * v_x
    u0(m, nhydro, k, j, i) += density * G * vx * bdt;
  });

  return;
}
