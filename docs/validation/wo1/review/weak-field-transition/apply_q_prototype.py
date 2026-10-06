from pathlib import Path
root=Path(__file__).parent/'source'
def edit(name, old, new, count=1):
 p=root/name; s=p.read_text(); assert s.count(old)==count,(name,s.count(old),old[:50]); p.write_text(s.replace(old,new))
edit('src/eos/ideal_c2p_mhd.hpp', '  const Real log_isotropic_anisotropy = 2.0*log(rho) - 3.0*log(bmag);\n  return rho*(log(p_perp) - log(p_parallel) + log_isotropic_anisotropy);', '  // SCRATCH Q prototype: IAN stores rho*log(p_perp/p_parallel).\n  return rho*(log(p_perp) - log(p_parallel));')
edit('src/eos/ideal_c2p_mhd.hpp', '''  const Real log_isotropic_anisotropy = 2.0*log(rho) - 3.0*log(bmag);
  Real log_p_ratio = anisotropy/rho - log_isotropic_anisotropy;
  // Preserve exact round trips for isotropic states without a tolerance or limiter.
  if (Kokkos::isfinite(anisotropy) &&
      anisotropy == rho*log_isotropic_anisotropy) {
    log_p_ratio = 0.0;
  }''', '  const Real log_p_ratio = anisotropy/rho;')
edit('src/eos/ideal_c2p_mhd.hpp', '    const Real log_ratio = log_exp - (2.0*log(original_density) - 3.0*log(bmag));', '    const Real log_ratio = log_exp;')
p=root/'src/mhd/rsolvers/hlle_cgl.hpp'; s=p.read_text(); start=s.index('    // Evaluate only the selected side,'); end=s.index('    Real fmutmp =',start); s=s[:start]+'''    // SCRATCH Q prototype: physical log pressure ratio has no B reference.
    const Real anis_upwind = (fdtmp >= 0.0) ? log(wl_ipp/wl_ipr)
                                          : log(wr_ipp/wr_ipr);
'''+s[end:]; p.write_text(s)
edit('src/mhd/rsolvers/llf_mhd_singlestate.hpp', '    u.mu = w.d*log(pperp/ppar*SQR(w.d)/(bmag*bsq));', '    u.mu = w.d*log(pperp/ppar);')
edit('src/mhd/rsolvers/llf_mhd_singlestate.hpp', '    u.mu = w.d*log(SQR(w.d)/(eos.bfloor*SQR(eos.bfloor)));', '    u.mu = 0.0;')
p=root/'src/mhd/rsolvers/llf_mhd_singlestate.hpp'; s=p.read_text(); start=s.index('  const MHDPrim1D &upwind ='); end=s.index('  flux.mu =',start); s=s[:start]+'''  Real specific_anis = (flux.d >= 0.0) ? ul.mu/ul.d : ur.mu/ur.d;
'''+s[end:]; p.write_text(s)
helper='''// SCRATCH prototype: 1D-only physical Q source using centered velocity gradients.
KOKKOS_INLINE_FUNCTION
Real LogRatioSource1D(const DvceArray5D<Real> &w,
                     const DvceArray5D<Real> &b,
                     const int m, const int k, const int j, const int i,
                     const Real dx, const Real bfloor) {
  const Real bx = b(m,IBX,k,j,i), by = b(m,IBY,k,j,i), bz = b(m,IBZ,k,j,i);
  const Real bsqr = SQR(bx) + SQR(by) + SQR(bz);
  if (bsqr <= SQR(bfloor)) return 0.0;
  const Real dvx = (w(m,IVX,k,j,i+1) - w(m,IVX,k,j,i-1))/(2.0*dx);
  const Real dvy = (w(m,IVY,k,j,i+1) - w(m,IVY,k,j,i-1))/(2.0*dx);
  const Real dvz = (w(m,IVZ,k,j,i+1) - w(m,IVZ,k,j,i-1))/(2.0*dx);
  return w(m,IDN,k,j,i)*(3.0*bx*(bx*dvx+by*dvy+bz*dvz)/bsqr-dvx);
}

'''
edit('src/eos/cgl_physics.hpp', '} // namespace cgl', helper+'} // namespace cgl')
for name in ['src/mhd/mhd_update.cpp','src/mhd/mhd_fofc.cpp']:
 edit(name, '#include "eos/eos.hpp"', '#include "eos/eos.hpp"\n#include "eos/cgl_physics.hpp"')
edit('src/mhd/mhd_update.cpp', '  auto &mbsize = pmy_pack->pmb->mb_size;\n\n  if (diagnose_nonfinite_rk_update)', '''  auto &mbsize = pmy_pack->pmb->mb_size;
  const auto q_w = w0;
  const auto q_b = bcc0;
  const auto q_eos = peos->eos_data;
  if (q_eos.is_cgl && (multi_d || has_cgl_lf_split)) {
    std::cerr << "Q prototype supports only 1D pure CGL" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (diagnose_nonfinite_rk_update)''')
edit('src/mhd/mhd_update.cpp', '      u0_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - beta_dt*divf(i);', '''      u0_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - beta_dt*divf(i);
      if (q_eos.is_cgl && n == IAN) {
        u0_(m,n,k,j,i) += beta_dt*cgl::LogRatioSource1D(
            q_w, q_b, m, k, j, i, mbsize.d_view(m).dx1, q_eos.bfloor);
      }''')
edit('src/mhd/mhd_fofc.cpp', '    auto &b1_ = b1;\n', '    auto &b1_ = b1;\n    const auto q_w = w0;\n    const auto q_eos = peos->eos_data;\n')
edit('src/mhd/mhd_fofc.cpp', '        utest_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - divf;', '''        utest_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - divf;
        if (q_eos.is_cgl && n == IAN) {
          utest_(m,n,k,j,i) += beta_dt*cgl::LogRatioSource1D(
              q_w, bcc0_, m, k, j, i, size.d_view(m).dx1, q_eos.bfloor);
        }''')
print('Patched scratch source only')
