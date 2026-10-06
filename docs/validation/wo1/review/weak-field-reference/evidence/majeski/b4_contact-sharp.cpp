#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  const Real velocity = pin->GetReal("problem", "weak_field_velocity");
  auto *pack = pmy_mesh_->pmb_pack;
  auto *pmhd = pack->pmhd;
  const auto indcs = pmy_mesh_->mb_indcs;
  const int is=indcs.is, ie=indcs.ie, js=indcs.js, je=indcs.je;
  const int ks=indcs.ks, ke=indcs.ke, nmb=pack->nmb_thispack;
  auto w=pmhd->w0; auto bcc=pmhd->bcc0; auto b=pmhd->b0;
  par_for("b4_init", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m,int k,int j,int i) {
    const bool weak=(i <= (is+ie)/2) == (velocity > 0.0);
    w(m,IDN,k,j,i)=1.0;
    w(m,IVX,k,j,i)=velocity;
    w(m,IVY,k,j,i)=w(m,IVZ,k,j,i)=0.0;
    w(m,IPR,k,j,i)=w(m,IPP,k,j,i)=weak ? 1.5 : 1.0;
    bcc(m,IBX,k,j,i)=bcc(m,IBZ,k,j,i)=0.0;
    bcc(m,IBY,k,j,i)=weak ? 1.0e-12 : 1.0;
  });
  par_for("b4_b1", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie+1,
  KOKKOS_LAMBDA(int m,int k,int j,int i) { b.x1f(m,k,j,i)=0.0; });
  par_for("b4_b2", DevExeSpace(), 0,nmb-1,ks,ke,js,je+1,is,ie,
  KOKKOS_LAMBDA(int m,int k,int j,int i) {
    const bool weak=(i <= (is+ie)/2) == (velocity > 0.0);
    b.x2f(m,k,j,i)=weak ? 1.0e-12 : 1.0;
  });
  par_for("b4_b3", DevExeSpace(), 0,nmb-1,ks,ke+1,js,je,is,ie,
  KOKKOS_LAMBDA(int m,int k,int j,int i) { b.x3f(m,k,j,i)=0.0; });
  pmhd->peos->PrimToCons(w,bcc,pmhd->u0,is,ie,js,je,ks,ke);
}
