from pathlib import Path
root=Path('/tmp/cgl-b4-transition-20261001/source')
(root/'src/mhd/b4_trace.hpp').write_text(r'''#pragma once
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
inline int b4_cycle = -1;
inline int b4_stage = -1;
inline bool B4TraceOn() { return std::getenv("B4_TRACE") != nullptr; }
inline void B4Snapshot(MeshBlockPack *pack, mhd::MHD *p, const char *phase, int stage) {
  if (!B4TraceOn() || pack->pmesh->ncycle >= 5) return;
  Kokkos::fence();
  auto u=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->u0);
  auto u1=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->u1);
  auto w=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->w0);
  auto bx=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->b0.x1f);
  auto by=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->b0.x2f);
  auto bz=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->b0.x3f);
  auto by1=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->b1.x2f);
  auto bc=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->bcc0);
  auto f=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->uflx.x1f);
  auto e=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),p->e3x1);
  static std::ofstream out("b4-state.csv");
  static bool header=false;
  if (!header) { out << "cycle,stage,phase,cell,dt,rho,mx,E,A,rho1,A1,By1,Bx,By,Bz,Bccy,ppar,pperp,FrhoL,FrhoR,FAL,FAR,FEL,FER,FByL,FByR\n"; header=true; }
  const auto &d=pack->pmesh->mb_indcs;
  out << std::setprecision(17);
  for (int i=d.is;i<=d.ie;++i) {
    int m=0,k=d.ks,j=d.js;
    out << pack->pmesh->ncycle << ',' << stage << ',' << phase << ',' << i-d.is << ',' << pack->pmesh->dt;
    for (int n : {IDN,IM1,IEN,IAN}) out << ',' << u(m,n,k,j,i);
    out << ',' << u1(m,IDN,k,j,i) << ',' << u1(m,IAN,k,j,i) << ',' << by1(m,k,j,i);
    out << ',' << .5*(bx(m,k,j,i)+bx(m,k,j,i+1)) << ',' << .5*(by(m,k,j,i)+by(m,k,j+1,i)) << ',' << .5*(bz(m,k,j,i)+bz(m,k+1,j,i));
    out << ',' << bc(m,IBY,k,j,i) << ',' << w(m,IPR,k,j,i) << ',' << w(m,IPP,k,j,i);
    for (int n : {IDN,IAN,IEN}) out << ',' << f(m,n,k,j,i) << ',' << f(m,n,k,j,i+1);
    out << ',' << -e(m,k,j,i) << ',' << -e(m,k,j,i+1) << '\n';
  }
  out.flush();
}
''')
p=root/'src/mhd/mhd_tasks.cpp';s=p.read_text().replace('#include "mhd/mhd.hpp"','#include "mhd/mhd.hpp"\n#include "mhd/b4_trace.hpp"')
s=s.replace('TaskStatus MHD::Fluxes(Driver *pdrive, int stage) {','TaskStatus MHD::Fluxes(Driver *pdrive, int stage) {\n  b4_cycle=pmy_pack->pmesh->ncycle; b4_stage=stage;')
s=s.replace('  // call FOFC if necessary','  B4Snapshot(pmy_pack,this,"flux_pre_fofc",stage);\n  // call FOFC if necessary')
a=s.index('TaskStatus MHD::Fluxes(');b=s.index('  return TaskStatus::complete;',a);s=s[:b]+'  B4Snapshot(pmy_pack,this,"flux_post_fofc",stage);\n'+s[b:]
s=s.replace('  peos->ConsToPrim(u0, b0, w0, bcc0, false, 0, n1m1, 0, n2m1, 0, n3m1);','  B4Snapshot(pmy_pack,this,"c2p_pre",stage);\n  peos->ConsToPrim(u0, b0, w0, bcc0, false, 0, n1m1, 0, n2m1, 0, n3m1);\n  B4Snapshot(pmy_pack,this,"c2p_post",stage);')
s=s.replace('  peos->Collisions(w0, bcc0, u0, pmy_pack->pmesh->dt, mode,','  B4Snapshot(pmy_pack,this,"coll_pre",stage);\n  peos->Collisions(w0, bcc0, u0, pmy_pack->pmesh->dt, mode,')
s=s.replace('  // Rates or walls can change fine-cell anisotropy after the last restriction.','  B4Snapshot(pmy_pack,this,"coll_post",stage);\n  // Rates or walls can change fine-cell anisotropy after the last restriction.')
p.write_text(s)
for filename,name in [('mhd_update.cpp','RKUpdate'),('mhd_ct.cpp','CT')]:
 p=root/'src/mhd'/filename;s=p.read_text().replace('#include "mhd.hpp"','#include "mhd.hpp"\n#include "mhd/b4_trace.hpp"')
 a=s.index(f'TaskStatus MHD::{name}(');b=s.index('  return TaskStatus::complete;',a)
 s=s[:b]+f'  B4Snapshot(pmy_pack,this,"{name.lower()}_post",stage);\n'+s[b:]
 p.write_text(s)
p=root/'src/eos/cgl_mhd.cpp';s=p.read_text().replace('#include "mhd/mhd.hpp"','#include "mhd/mhd.hpp"\n#include "mhd/b4_trace.hpp"')
s=s.replace('    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);',r'''    const Real b4_old_a=u.mu;
    const Real b4_old_e=u.e;
    SingleC2P_CGLMHD(u, eos, w, dfloor_used, efloor_used, tfloor_used, bfloor_used);
    if (B4TraceOn() && b4_cycle < 5 && i>=3 && i<131) {
      printf("B4EOS %d %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %d %d %d %d\n",
             b4_cycle,b4_stage,only_testfloors,i-3,u.d,u.by,b4_old_a,u.mu,b4_old_e,w.e,w.pp,
             dfloor_used,efloor_used,tfloor_used,bfloor_used);
    }''')
p.write_text(s)
p=root/'src/mhd/rsolvers/hlle_cgl.hpp';s=p.read_text().replace('#include "mhd/mhd.hpp"','#include "mhd/mhd.hpp"\n#include "mhd/b4_trace.hpp"')
s=s.replace('    Real fmutmp = fdtmp*anis_upwind;',r'''    Real fmutmp = fdtmp*anis_upwind;
    if (B4TraceOn() && b4_cycle < 5 && i>=61 && i<=78) {
      printf("B4FACE %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
             b4_cycle,b4_stage,i-3,wl_idn,wr_idn,wl_ipr,wl_ipp,wr_ipr,wr_ipp,bmagl,bmagr,fdtmp,fmutmp);
    }''')
p.write_text(s)
print('Instrumented scratch source only:',root)
