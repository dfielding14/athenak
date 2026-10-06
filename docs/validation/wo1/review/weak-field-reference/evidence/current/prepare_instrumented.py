from pathlib import Path
import subprocess,difflib,json
root=Path(__file__).parent;src=root/'source'
repo=Path('/Users/dbf75/.codex/worktrees/bf22/athenak-DF')
archive=repo/'docs/validation/wo1/review/weak-field-transition'
subprocess.run(['patch','-p1','-i',str(archive/'instrumentation.patch')],cwd=src,check=True)
def edit(name,old,new,count=1):
 p=src/name;s=p.read_text();assert s.count(old)==count,(name,s.count(old),old[:45]);p.write_text(s.replace(old,new))
edit('src/mhd/b4_trace.hpp','#include <cstdlib>','#include <cstdlib>\n#include <cstring>')
edit('src/mhd/b4_trace.hpp','  if (!B4TraceOn() || pack->pmesh->ncycle >= 5) return;', '  if (!B4TraceOn()) return;\n  if (pack->pmesh->ncycle >= 5 && std::strcmp(phase,"coll_post") != 0) return;')
edit('src/mhd/b4_trace.hpp','cycle,stage,phase,cell,dt,rho','cycle,stage,phase,cell,time,dt,rho')
edit('src/mhd/b4_trace.hpp',"<< i-d.is << ',' << pack->pmesh->dt;", "<< i-d.is << ',' << pack->pmesh->time << ',' << pack->pmesh->dt;")
# Expand diagnostic index range to each test's active cells. The face trace window
# remains around the128-cell discontinuity; full cell/flux snapshots cover everycell.
edit('src/eos/cgl_mhd.cpp','b4_cycle < 5 && i>=3 && i<131','b4_cycle < 5 && i>=3 && i<ni-3')
controls='''#pragma once
#include <cstdlib>
// SCRATCH-ONLY CPU comparison controls. Unset selects unchanged WO1 physics.
inline bool B4OriginalFace() { return std::getenv("B4_ORIGINAL_FACE") != nullptr; }
inline bool B4NoWall() { return std::getenv("B4_NO_WALL") != nullptr; }
'''
(src/'eos').mkdir(exist_ok=True) if False else None
(src/'src/eos/b4_controls.hpp').write_text(controls)
edit('src/mhd/rsolvers/hlle_cgl.hpp','#include "athena.hpp"','#include "athena.hpp"\n#include "eos/b4_controls.hpp"')
edit('src/mhd/rsolvers/llf_mhd_singlestate.hpp','#include "coordinates/cartesian_ks.hpp"','#include "coordinates/cartesian_ks.hpp"\n#include "eos/b4_controls.hpp"')
edit('src/eos/ideal_c2p_mhd.hpp','#include "eos/cgl_physics.hpp"','#include "eos/cgl_physics.hpp"\n#include "eos/b4_controls.hpp"')
edit('src/eos/ideal_c2p_mhd.hpp','                                const bool backup) {\n  const Real bsqr', '                                const bool backup) {\n  if (B4NoWall() && !backup) return anisotropy;\n  const Real bsqr')
edit('src/eos/ideal_c2p_mhd.hpp','void SingleCollWalls_CGLMHD(MHDPrim1D &w, const EOS_Data &eos, const bool backup) {','void SingleCollWalls_CGLMHD(MHDPrim1D &w, const EOS_Data &eos, const bool backup) {\n  if (B4NoWall() && !backup) return;')
p=src/'src/mhd/rsolvers/hlle_cgl.hpp';s=p.read_text();start=s.index('    //--- Step 2. Isotropize');end=s.index('    //--- Step 3.',start);original=s[start:end]
s=s[:start]+'''    const bool original_face = B4OriginalFace();
    Real original_al = 0.0, original_ar = 0.0;
    Real fhl = 1.0, fhr = 1.0;
    if (original_face) {
      fhl = 1.0 + (wl_ipp - wl_ipr)/(2.0*pbl);
      fhr = 1.0 + (wr_ipp - wr_ipr)/(2.0*pbr);
      original_al = wl_idn*log(wl_ipp/wl_ipr*SQR(wl_idn)/(bmagl*SQR(bmagl)));
      original_ar = wr_idn*log(wr_ipp/wr_ipr*SQR(wr_idn)/(bmagr*SQR(bmagr)));
      if (bmagl < bfloor || bmagr < bfloor) {
        fhl = fhr = 1.0;
        wl_ipr = TWO_3RDS*wl_ipp + ONE_3RD*wl_ipr;
        wl_ipp = wl_ipr;
        wr_ipr = TWO_3RDS*wr_ipp + ONE_3RD*wr_ipr;
        wr_ipp = wr_ipr;
        original_al = wl_idn*log(SQR(wl_idn)/(bfloor*SQR(bfloor)));
        original_ar = wr_idn*log(SQR(wr_idn)/(bfloor*SQR(bfloor)));
      }
    } else {
'''+original.replace('    Real fhl = 1.0, fhr = 1.0;\n','')+'''    }

'''+s[end:]
old='''    Real anis_upwind;
    if (fdtmp >= 0.0) {'''
new='''    Real anis_upwind;
    if (original_face) {
      original_al /= wl_idn;
      original_ar /= wr_idn;
      anis_upwind = fdtmp >= 0.0 ? original_al : original_ar;
    } else if (fdtmp >= 0.0) {'''
assert s.count(old)==1;s=s.replace(old,new);p.write_text(s)
edit('src/mhd/rsolvers/llf_mhd_singlestate.hpp','  if (sqrt(SQR(bxi) + SQR(upwind.by) + SQR(upwind.bz)) <= eos.bfloor) {','  if (!B4OriginalFace() &&\n      sqrt(SQR(bxi) + SQR(upwind.by) + SQR(upwind.bz)) <= eos.bfloor) {')
# Optional resolved profile for a later fixed-width refinement study.
p=src/'src/pgen/tests/cgl_fofc.cpp';s=p.read_text();key='  const bool weak_field = mode == "weak_field_transport";';assert s.count(key)==1;s=s.replace(key,key+'\n  const Real profile_width = pin->GetOrAddReal("problem", "profile_width", 0.0);')
key='  auto b0 = pmhd->b0;\n\n  par_for("cgl_fofc_e2e_init"';new='''  auto b0 = pmhd->b0;
  const Real x_min = pmy_mesh_->mesh_size.x1min;
  const Real x_max = pmy_mesh_->mesh_size.x1max;
  const Real x_step = (x_max-x_min)/indcs.nx1;

  par_for("cgl_fofc_e2e_init"''';assert s.count(key)==1;s=s.replace(key,new)
key='''      bcc0(m,IBY,k,j,i) = weak_side ? 1.0e-12 : 1.0;''';new=key+'''
      if (profile_width > 0.0) {
        const Real xcoord = x_min+(i-is+0.5)*x_step;
        const Real sign = weak_velocity >= 0.0 ? 1.0 : -1.0;
        const Real field = 1.0e-12+(1.0-1.0e-12)*0.5*
            (1.0+tanh(sign*(xcoord-0.5*(x_min+x_max))/profile_width));
        bcc0(m,IBY,k,j,i) = field;
        w0(m,IPR,k,j,i) = w0(m,IPP,k,j,i) = 1.5-0.5*SQR(field);
      }''';assert s.count(key)==1;s=s.replace(key,new)
key='''      b0.x2f(m,k,j,i) = weak_side ? 1.0e-12 : 1.0;''';new=key+'''
      if (profile_width > 0.0) {
        const Real xcoord = x_min+(i-is+0.5)*x_step;
        const Real sign = weak_velocity >= 0.0 ? 1.0 : -1.0;
        b0.x2f(m,k,j,i) = 1.0e-12+(1.0-1.0e-12)*0.5*
            (1.0+tanh(sign*(xcoord-0.5*(x_min+x_max))/profile_width));
      }''';assert s.count(key)==1;s=s.replace(key,new);p.write_text(s)
print('Prepared scratchinstrumentation and4variant controls.')
