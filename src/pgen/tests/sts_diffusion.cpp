//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sts_diffusion.cpp
//! \brief Validation and visualization initial conditions for parabolic STS operators.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

namespace {

bool create_reference_state = false;

void ScalarModeErrors(ParameterInput *pin, Mesh *pm);

void FatalCase(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << message << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \brief Set up scalar Fourier modes, a moving scalar blob, a conducting front, or a
//! temperature-dependent viscous shear layer.

void ProblemGenerator::STSDiffusion(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    FatalCase("sts_diffusion currently requires a <hydro> block.");
  }
  auto *phydro = pmbp->phydro;
  EOS_Data &eos = phydro->peos->eos_data;
  if (!eos.is_ideal) {
    FatalCase("sts_diffusion requires <hydro>/eos = ideal.");
  }

  const std::string sts_case = pin->GetOrAddString("problem", "sts_case", "scalar_modes");
  const bool scalar_modes = (sts_case == "scalar_modes");
  const bool scalar_blob = (sts_case == "scalar_blob");
  const bool thermal_front = (sts_case == "thermal_front");
  const bool viscous_shear = (sts_case == "viscous_shear");
  if (!(scalar_modes || scalar_blob || thermal_front || viscous_shear)) {
    FatalCase("Unknown <problem>/sts_case = '" + sts_case + "'.");
  }

  if ((scalar_modes && phydro->nscalars < 2) ||
      (scalar_blob && phydro->nscalars < 1)) {
    FatalCase("scalar_modes requires nscalars >= 2 and "
              "scalar_blob requires nscalars >= 1.");
  }
  if (scalar_modes) {
    pgen_final_func = ScalarModeErrors;
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nmb = pmbp->nmb_thispack;
  const int nhydro = phydro->nhydro;
  const int nscalars = phydro->nscalars;
  auto &size = pmbp->pmb->mb_size;
  auto &u = create_reference_state ? phydro->u1 : phydro->u0;

  const Real x1min = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  const Real x2min = pmy_mesh_->mesh_size.x2min;
  const Real x2max = pmy_mesh_->mesh_size.x2max;
  const Real lx = x1max - x1min;
  const Real ly = x2max - x2min;
  const Real gm1 = eos.gamma - 1.0;
  const Real p0 = pin->GetOrAddReal("problem", "pressure", 1.0);
  const Real amp = pin->GetOrAddReal("problem", "amp", 0.2);
  const Real mode = pin->GetOrAddReal("problem", "mode", 1.0);
  const Real kx = 2.0*M_PI*mode/lx;
  const Real final_time = create_reference_state ? pmy_mesh_->time : 0.0;

  const Real kappa_default = scalar_modes ?
      pin->GetReal("hydro", "scalar_diffusivity") : 0.0;
  const Real kappa0 = scalar_modes ?
      pin->GetOrAddReal("hydro", "scalar_diffusivity_0", kappa_default) : 0.0;
  const Real kappa1 = scalar_modes ?
      pin->GetOrAddReal("hydro", "scalar_diffusivity_1", kappa_default) : 0.0;
  const Real scalar_decay0 = std::exp(-kappa0*kx*kx*final_time);
  const Real scalar_decay1 = std::exp(-kappa1*kx*kx*final_time);

  const Real sigma = pin->GetOrAddReal("problem", "sigma", 0.08*lx);
  const Real xcenter = pin->GetOrAddReal("problem", "xcenter", x1min + 0.25*lx);
  const Real ycenter = pin->GetOrAddReal("problem", "ycenter", x2min + 0.5*ly);
  const Real vx = pin->GetOrAddReal("problem", "vx", scalar_blob ? 0.35 : 0.0);
  const Real vy = pin->GetOrAddReal("problem", "vy", scalar_blob ? 0.15 : 0.0);
  const Real tcold = pin->GetOrAddReal("problem", "tcold", 0.25);
  const Real thot = pin->GetOrAddReal("problem", "thot", 2.0);
  const Real front_width = pin->GetOrAddReal("problem", "front_width", 0.04*lx);

  par_for("sts_diffusion_init", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real rho = 1.0;
    Real temp = 1.0;
    Real v1 = vx;
    Real v2 = vy;
    Real s0 = 0.0;
    Real s1 = 0.0;

    if (scalar_modes) {
      s0 = 0.5 + amp*scalar_decay0*std::sin(kx*(x - x1min));
      s1 = 0.5 + amp*scalar_decay1*std::cos(kx*(x - x1min));
    } else if (scalar_blob) {
      Real dx = x - xcenter;
      Real dy = y - ycenter;
      if (dx > 0.5*lx) dx -= lx;
      if (dx < -0.5*lx) dx += lx;
      if (dy > 0.5*ly) dy -= ly;
      if (dy < -0.5*ly) dy += ly;
      s0 = std::exp(-0.5*(dx*dx + dy*dy)/(sigma*sigma));
    } else if (thermal_front) {
      temp = tcold + 0.5*(thot - tcold)*
          (1.0 - std::tanh((x - xcenter)/front_width));
      rho = p0/temp;
      v1 = 0.0;
      v2 = 0.0;
    } else if (viscous_shear) {
      temp = tcold + 0.5*(thot - tcold)*
          (1.0 + std::sin(kx*(x - x1min)));
      rho = p0/temp;
      v1 = 0.0;
      v2 = amp*std::sin(kx*(x - x1min));
    }

    u(m,IDN,k,j,i) = rho;
    u(m,IM1,k,j,i) = rho*v1;
    u(m,IM2,k,j,i) = rho*v2;
    u(m,IM3,k,j,i) = 0.0;
    u(m,IEN,k,j,i) = p0/gm1 + 0.5*rho*(v1*v1 + v2*v2);
    for (int n = 0; n < nscalars; ++n) {
      u(m,nhydro+n,k,j,i) = 0.0;
    }
    if (scalar_modes) {
      u(m,nhydro,k,j,i) = rho*s0;
      u(m,nhydro+1,k,j,i) = rho*s1;
    } else if (scalar_blob) {
      u(m,nhydro,k,j,i) = rho*s0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \brief Compare the two scalar modes with their exact diffusion decay solution.

namespace {

void ScalarModeErrors(ParameterInput *pin, Mesh *pm) {
  create_reference_state = true;
  pm->pgen->STSDiffusion(pin, false);
  create_reference_state = false;

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto *phydro = pmbp->phydro;
  const int ns0 = phydro->nhydro;
  const int ns1 = ns0 + 1;
  auto u0 = phydro->u0;
  auto u1 = phydro->u1;
  auto size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx1*nx2*nx3;
  const int nkji = nx1*nx2*nx3;
  const int nji = nx1*nx2;
  Real l1_s0 = 0.0, l1_s1 = 0.0, linf_s0 = 0.0, linf_s1 = 0.0;

  Kokkos::parallel_reduce("sts_scalar_l1_0",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    sum += vol*std::fabs(u0(m,ns0,k,j,i) - u1(m,ns0,k,j,i));
  }, Kokkos::Sum<Real>(l1_s0));
  Kokkos::parallel_reduce("sts_scalar_l1_1",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    sum += vol*std::fabs(u0(m,ns1,k,j,i) - u1(m,ns1,k,j,i));
  }, Kokkos::Sum<Real>(l1_s1));
  Kokkos::parallel_reduce("sts_scalar_linf_0",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &max_error) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    max_error = fmax(max_error, std::fabs(u0(m,ns0,k,j,i) - u1(m,ns0,k,j,i)));
  }, Kokkos::Max<Real>(linf_s0));
  Kokkos::parallel_reduce("sts_scalar_linf_1",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &max_error) {
    const int m = idx/nkji;
    const int rem_m = idx - m*nkji;
    const int k0 = rem_m/nji;
    const int rem_k = rem_m - k0*nji;
    const int j0 = rem_k/nx1;
    const int i = rem_k - j0*nx1 + is;
    const int j = j0 + js;
    const int k = k0 + ks;
    max_error = fmax(max_error, std::fabs(u0(m,ns1,k,j,i) - u1(m,ns1,k,j,i)));
  }, Kokkos::Max<Real>(linf_s1));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &l1_s0, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &l1_s1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &linf_s0, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &linf_s1, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

  const Real volume = (pm->mesh_size.x1max - pm->mesh_size.x1min)*
      (pm->mesh_size.x2max - pm->mesh_size.x2min)*
      (pm->mesh_size.x3max - pm->mesh_size.x3min);
  l1_s0 /= volume;
  l1_s1 /= volume;

  if (global_variable::my_rank == 0) {
    const std::string fname = pin->GetString("job", "basename") + "-scalar-errors.dat";
    FILE *pfile = std::fopen(fname.c_str(), "w");
    if (pfile == nullptr) {
      FatalCase("Could not write STS scalar error output.");
    }
    std::fprintf(pfile, "# Nx1 Ncycle L1_s0 L1_s1 Linf_s0 Linf_s1\n");
    std::fprintf(pfile, "%d %d %.16e %.16e %.16e %.16e\n",
                 pm->mesh_indcs.nx1, pm->ncycle, l1_s0, l1_s1, linf_s0, linf_s1);
    std::fclose(pfile);
  }
}

} // namespace
