//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_init.cpp
//  \brief Problem generator for a turbulent box
#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "utils/random.hpp"

//----------------------------------------------------------------------------------------
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  // Turbulence parameters
  int nlow   = pin->GetOrAddInteger("problem", "cgm_turb_nlow", 1);
  int nhigh  = pin->GetOrAddInteger("problem", "cgm_turb_nhigh", 8);
  Real expo  = pin->GetOrAddReal("problem", "cgm_turb_expo", 5.0/3.0);
  Real v_rms = pin->GetOrAddReal("problem", "cgm_turb_rms", 0.1);
  Real cgm_turb_xscale = pin->GetOrAddReal("problem", "cgm_turb_xscale", 0.01);
  Real cgm_turb_yscale = pin->GetOrAddReal("problem", "cgm_turb_yscale", 0.01);
  Real cgm_turb_zscale = pin->GetOrAddReal("problem", "cgm_turb_zscale", 0.01);

  // Initialize random state
  RNG_State rstate;
  rstate.idum = -1;

  // Domain size
  Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real ly = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  Real lz = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  Real dkx = 2.0*M_PI/lx;
  Real dky = 2.0*M_PI/ly;
  Real dkz = 2.0*M_PI/lz;

  // Count modes
  int nmodes = 0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
	  nmodes++;
	}
  }}}

  // Allocate arrays
  DualArray2D<Real> k_modes, aka, akb;
  Kokkos::realloc(k_modes, 3, nmodes);
  Kokkos::realloc(aka, 3, nmodes);
  Kokkos::realloc(akb, 3, nmodes);

  // Generate modes
  int nmode = 0;
  Real total_energy = 0.0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
          Real kx = dkx*nkx, ky = dky*nky, kz = dkz*nkz;
          Real kiso = sqrt(kx*kx + ky*ky + kz*kz);

	  k_modes.h_view(0, nmode) = kx;
          k_modes.h_view(1, nmode) = ky;
	  k_modes.h_view(2, nmode) = kz;

          Real norm = 1.0/pow(kiso, (expo+2.0)/2.0);

          Real aval[3], bval[3];
          for (int dir = 0; dir < 3; dir++) {
            aval[dir] = norm * RanGaussianSt(&rstate);
            bval[dir] = norm * RanGaussianSt(&rstate);
          }
 
	  Real k_dirs[3] = {kx, ky, kz};
          Real ka = kx*aval[0] + ky*aval[1] + kz*aval[2];
          Real kb = kx*bval[0] + ky*bval[1] + kz*bval[2];

          for (int dir = 0; dir < 3; dir++) {
            aval[dir] -= k_dirs[dir]*ka/(kiso*kiso);
            bval[dir] -= k_dirs[dir]*kb/(kiso*kiso);

	    aka.h_view(dir,nmode) = aval[dir];
            akb.h_view(dir,nmode) = bval[dir];
		
	    total_energy += 0.5*(aval[dir]*aval[dir] + bval[dir]*bval[dir]);
          }
          nmode++;	  
        }
      }
    }
  }

  Real v_norm = v_rms/sqrt(total_energy);

  k_modes.template modify<HostMemSpace>();
  k_modes.template sync<DevExeSpace>();
  aka.template modify<HostMemSpace>();
  aka.template sync<DevExeSpace>();
  akb.template modify<HostMemSpace>();
  akb.template sync<DevExeSpace>();

  // Optional seeds of cold gas, both at pressure p0: log-normal density perturbations,
  // ln(rho/rho0) a random field of rms drho_rms on box modes [drho_nlow, drho_nhigh] with
  // the velocity field's spectral slope; and a cloud of density cloud_chi*rho0 and radius
  // cloud_radius at the origin.
  Real drho_rms = pin->GetOrAddReal("problem", "drho_rms", 0.0);
  int drho_nlow = pin->GetOrAddInteger("problem", "drho_nlow", 4);
  int drho_nhigh = pin->GetOrAddInteger("problem", "drho_nhigh", 8);
  Real cloud_radius = pin->GetOrAddReal("problem", "cloud_radius", 0.0);
  Real cloud_chi = pin->GetOrAddReal("problem", "cloud_chi", 1.0);

  // Count density modes (none without perturbations)
  int ndmodes = 0;
  if (drho_rms > 0.0) {
    for (int nkx = -drho_nhigh; nkx <= drho_nhigh; nkx++) {
      for (int nky = -drho_nhigh; nky <= drho_nhigh; nky++) {
        for (int nkz = -drho_nhigh; nkz <= drho_nhigh; nkz++) {
          if (nkx == 0 && nky == 0 && nkz == 0) continue;
          int nsqr = nkx*nkx + nky*nky + nkz*nkz;
          if (nsqr >= drho_nlow*drho_nlow && nsqr <= drho_nhigh*drho_nhigh) ndmodes++;
    }}}
  }

  // Density modes come from their own random sequence, so the velocities don't change
  DualArray2D<Real> dk_modes, dka;  // dka(0,n), dka(1,n): cos and sin amplitudes
  Kokkos::realloc(dk_modes, 3, ndmodes);
  Kokkos::realloc(dka, 2, ndmodes);
  Real dvar = 0.0;
  if (drho_rms > 0.0) {
    RNG_State dstate;
    dstate.idum = -2;
    int nd = 0;
    for (int nkx = -drho_nhigh; nkx <= drho_nhigh; nkx++) {
      for (int nky = -drho_nhigh; nky <= drho_nhigh; nky++) {
        for (int nkz = -drho_nhigh; nkz <= drho_nhigh; nkz++) {
          if (nkx == 0 && nky == 0 && nkz == 0) continue;
          int nsqr = nkx*nkx + nky*nky + nkz*nkz;
          if (nsqr >= drho_nlow*drho_nlow && nsqr <= drho_nhigh*drho_nhigh) {
            Real kx = dkx*nkx, ky = dky*nky, kz = dkz*nkz;
            Real norm = 1.0/pow(sqrt(kx*kx + ky*ky + kz*kz), (expo+2.0)/2.0);
            dk_modes.h_view(0,nd) = kx;
            dk_modes.h_view(1,nd) = ky;
            dk_modes.h_view(2,nd) = kz;
            dka.h_view(0,nd) = norm*RanGaussianSt(&dstate);
            dka.h_view(1,nd) = norm*RanGaussianSt(&dstate);
            dvar += 0.5*(SQR(dka.h_view(0,nd)) + SQR(dka.h_view(1,nd)));
            nd++;
          }
    }}}
  }
  Real d_norm = (ndmodes > 0) ? drho_rms/sqrt(dvar) : 0.0;

  dk_modes.template modify<HostMemSpace>();
  dk_modes.template sync<DevExeSpace>();
  dka.template modify<HostMemSpace>();
  dka.template sync<DevExeSpace>();

  // Initialize density, pressure, velocities
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  // Passive scalar 0 is the metal mass density, as GOTHAM's u0(IZS) = Z*Zsol*rho
  Real zmet = 0.02*pin->GetOrAddReal("problem", "metallicity", 1.0/3.0);
  int nhyd = pmbp->phydro->nhydro;
  int nscal = pmbp->phydro->nscalars;

  int &is = indcs.is, &ie = indcs.ie;
  int &js = indcs.js, &je = indcs.je;
  int &ks = indcs.ks, &ke = indcs.ke;

  par_for("pgen_turb_init", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // Get cell coordinates
    Real &x1min = size.d_view(m).x1min, &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min, &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min, &x3max = size.d_view(m).x3max;

    Real x = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real y = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real z = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Compute turbulent velocities by summing Fourier modes
    Real vx = 0.0, vy = 0.0, vz = 0.0;
    for (int n = 0; n < nmodes; n++) {
      Real phase = k_modes.d_view(0,n)*x + k_modes.d_view(1,n)*y + k_modes.d_view(2,n)*z;
      Real cos_phase = cos(phase);
      Real sin_phase = sin(phase);

      vx += aka.d_view(0,n)*cos_phase - akb.d_view(0,n)*sin_phase;
      vy += aka.d_view(1,n)*cos_phase - akb.d_view(1,n)*sin_phase;
      vz += aka.d_view(2,n)*cos_phase - akb.d_view(2,n)*sin_phase;
    }

    // Attenuate in the center by 1 - Gaussian
    Real att = 1.0 - exp(-0.5 * ( SQR(x)/SQR(cgm_turb_xscale) 
 		                + SQR(y)/SQR(cgm_turb_yscale)
  		                + SQR(z)/SQR(cgm_turb_zscale)));

    // Normalize to desired RMS velocity
    vx *= v_norm*att; vy *= v_norm*att; vz *= v_norm*att;

    // Density: log-normal perturbations (mean rho0), then the cloud
    Real lnd = 0.0;
    for (int n = 0; n < ndmodes; n++) {
      Real phase = dk_modes.d_view(0,n)*x + dk_modes.d_view(1,n)*y
                 + dk_modes.d_view(2,n)*z;
      lnd += dka.d_view(0,n)*cos(phase) - dka.d_view(1,n)*sin(phase);
    }
    Real rho = rho0*exp(d_norm*lnd - 0.5*SQR(drho_rms));
    if (SQR(x) + SQR(y) + SQR(z) < SQR(cloud_radius)) rho = cloud_chi*rho0;

    // Set conserved variables; the pressure is p0 everywhere
    u0(m,IDN,k,j,i) = rho;
    if (nscal > 0) u0(m,nhyd,k,j,i) = rho*zmet;
    u0(m,IM1,k,j,i) = rho * vx;
    u0(m,IM2,k,j,i) = rho * vy;
    u0(m,IM3,k,j,i) = rho * vz;

    Real ekin = 0.5*rho*(vx*vx + vy*vy + vz*vz);
    u0(m,IEN,k,j,i) = p0/gm1 + ekin;
  });
}
