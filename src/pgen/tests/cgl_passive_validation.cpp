// Independent passive/isothermal validation initial data. Built-in regression pgen.
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

void ProblemGenerator::CGLPassiveValidation(ParameterInput *pin, const bool restart) {
  if (restart) return;
  auto *pack = pmy_mesh_->pmb_pack;
  if (pack->pmhd == nullptr || pack->phydro != nullptr) {
    std::cout << "passive_validation requires MHD only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const auto eos = pack->pmhd->peos->eos_data;
  if (eos.is_ideal && !eos.is_cgl) {
    std::cout << "passive_validation requires isothermal or CGL EOS" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const auto indcs = pmy_mesh_->mb_indcs;
  const int is=indcs.is, ie=indcs.ie, js=indcs.js, je=indcs.je, ks=indcs.ks, ke=indcs.ke;
  const int nmb = pack->nmb_thispack;
  const auto size = pack->pmb->mb_size;
  const auto w = pack->pmhd->w0;
  const auto bcc = pack->pmhd->bcc0;
  const auto b = pack->pmhd->b0;
  const int mode = pin->GetOrAddInteger("problem", "case", 0);
  const Real amplitude = pin->GetOrAddReal("problem", "amplitude", 0.1);
  const Real thermal_amplitude = pin->GetOrAddReal("problem", "thermal_amplitude", 0.2);
  const Real p0 = pin->GetOrAddReal("problem", "p0", 2.0);
  const Real advect = pin->GetOrAddReal("problem", "advect_velocity", 0.4);
  const Real bx = pin->GetOrAddReal("problem", "bx", 1.0);
  const Real by = pin->GetOrAddReal("problem", "by", mode == 2 ? 0.0 : 1.0);
  const Real bz = pin->GetOrAddReal("problem", "bz", 0.0);
  const Real x0 = pmy_mesh_->mesh_size.x1min;
  const Real kx = 2.0*M_PI/(pmy_mesh_->mesh_size.x1max-x0);
  par_for("passive_validation_prim", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real c = cos(kx*(x-x0)), s = sin(kx*(x-x0));
    Real rho=1.0, vx=0.0, vy=0.0, vz=0.0, pp=p0, pt=p0;
    if (mode == 0) {
      rho = 1.0 + amplitude*c;
      vx = amplitude*s;
      vy = 0.7*amplitude*c;
      vz = 0.3*amplitude*sin(2.0*kx*(x-x0));
      pp = p0*(1.0 + thermal_amplitude*c);
      pt = p0*(1.0 + 0.5*thermal_amplitude*s);
    } else if (mode == 1) {
      vy = amplitude*s;
      pp = p0 - TWO_3RDS*thermal_amplitude*c;
      pt = p0 + ONE_3RD*thermal_amplitude*c;
    } else if (mode == 2) {
      rho = 1.0 + amplitude*c;
      vx = eos.iso_cs*amplitude*c;
      pp = p0*(1.0 + 3.0*thermal_amplitude*c);
      pt = p0*(1.0 + thermal_amplitude*c);
    } else if (mode == 3) {
      vx = advect;
      pp = p0*(1.0 + thermal_amplitude*c);
      pt = p0*(1.0 + 0.5*thermal_amplitude*s);
    }
    w(m,IDN,k,j,i)=rho; w(m,IVX,k,j,i)=vx;
    w(m,IVY,k,j,i)=vy; w(m,IVZ,k,j,i)=vz;
    if (eos.is_cgl) { w(m,IPR,k,j,i)=pp; w(m,IPP,k,j,i)=pt; }
    bcc(m,IBX,k,j,i)=bx; bcc(m,IBY,k,j,i)=by; bcc(m,IBZ,k,j,i)=bz;
  });
  par_for("passive_validation_bx", DevExeSpace(), 0,nmb-1, ks,ke, js,je, is,ie+1,
  KOKKOS_LAMBDA(int m,int k,int j,int i) { b.x1f(m,k,j,i)=bx; });
  par_for("passive_validation_by", DevExeSpace(), 0,nmb-1, ks,ke, js,je+1, is,ie,
  KOKKOS_LAMBDA(int m,int k,int j,int i) { b.x2f(m,k,j,i)=by; });
  par_for("passive_validation_bz", DevExeSpace(), 0,nmb-1, ks,ke+1, js,je, is,ie,
  KOKKOS_LAMBDA(int m,int k,int j,int i) { b.x3f(m,k,j,i)=bz; });
  pack->pmhd->peos->PrimToCons(w,bcc,pack->pmhd->u0,is,ie,js,je,ks,ke);
}
